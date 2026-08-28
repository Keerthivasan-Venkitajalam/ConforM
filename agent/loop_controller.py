"""Closed-loop controller: observe -> decide -> execute -> update -> repeat.

Safety features (master prompt #21):
  - MAX_ITERATIONS (default 5)
  - compute budget in seconds
  - duplicate action detection via deterministic input hashes
  - per-tool timeouts (enforced inside the tool wrappers via subprocess timeout)
  - artifact validation after every step
  - failures are recorded into scientific memory, never silently recovered
"""
from __future__ import annotations

import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from agent.discovery_score import DiscoveryScoreInputs, DiscoveryScoreWeights, compute_discovery_score
from agent.policies import Decision, Policy, PolicyConfig
from agent.state import Action, ExperimentState, action_hash
from db.repository import Repository
from pipelines import engines
from tools import reinvent_tool
from tools.bioemu_tool import cuda_available

# Actions that constitute a real closed-loop experiment cycle (they consume
# significant compute and produce new ligand evidence). Setup stages do not
# count against MAX_ITERATIONS -- see ExperimentState.iteration.
EXPERIMENT_ACTIONS = {Action.SCREEN_LIGANDS, Action.OPTIMIZE_LIGAND, Action.VALIDATE_POSE}

# Hard cap on total actions, independent of iteration accounting, so a policy
# bug can never produce an unbounded loop.
MAX_TOTAL_STEPS = 20


class ClosedLoopAgent:
    def __init__(self, config_path: Path, ligand_csv: Path, out_root: Path = Path("artifacts"),
                 mode: str = "conform-agent", repo: Repository | None = None,
                 max_iterations: int | None = None, budget_seconds: float = 3600.0):
        self.cfg = yaml.safe_load(Path(config_path).read_text())
        self.config_path = Path(config_path)
        self.ligand_csv = Path(ligand_csv)
        self.mode = mode
        self.experiment_id = f"{self.cfg['target']['name'].lower()}_{mode}_{int(time.time())}"
        self.dir = Path(out_root) / self.experiment_id
        self.dirs = {name: self.dir / name for name in
                     ("input", "structures", "pockets", "ligands", "docking", "metrics", "report")}
        for d in self.dirs.values():
            d.mkdir(parents=True, exist_ok=True)

        self.repo = repo or Repository()
        self.repo.create_experiment(self.experiment_id, self.cfg["target"]["name"], mode, self.cfg)

        self.policy = Policy(PolicyConfig(
            max_iterations=max_iterations or self.cfg.get("agent", {}).get("max_iterations", 5)))
        self.state = ExperimentState(target=self.cfg["target"]["name"], budget_seconds=budget_seconds)
        self.log: list[str] = []
        self.weights = DiscoveryScoreWeights(**self.cfg["discovery_score"]["weights"])
        self.ensemble_provider = None
        self.optimizer_mode = None

    # ------------------------------------------------------------------
    def event(self, msg: str):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        self.log.append(line)
        print(line, flush=True)

    def run(self) -> dict:
        t_start = time.time()
        self.event(f"Target loaded: {self.state.target} ({self.cfg['target']['mutation']}) "
                   f"| mode={self.mode} | experiment={self.experiment_id}")

        while True:
            if self.state.step_count >= MAX_TOTAL_STEPS:
                self.event(f"Hard step cap ({MAX_TOTAL_STEPS}) reached; terminating loop.")
                break
            decision = self.policy.decide(self.state)
            self.event(f"DECISION (iter {self.state.iteration}): {decision.action.value} "
                       f"-- {decision.rationale}")

            if decision.action is Action.STOP:
                self.repo.log_step(self.experiment_id, self.state.iteration, Action.STOP.value,
                                    action_hash(Action.STOP, {}),
                                    interpretation=decision.rationale, next_action=None)
                break

            h = action_hash(decision.action, decision.params)
            if self.repo.has_completed(self.experiment_id, h):
                self.event(f"  duplicate action detected ({h}); refusing to re-run. Stopping.")
                break

            t0 = time.time()
            try:
                metrics, interpretation = self.execute(decision)
                failure = None
            except Exception as exc:  # noqa: BLE001
                metrics, interpretation, failure = {}, None, f"{type(exc).__name__}: {exc}"
                self.event(f"  ACTION FAILED: {failure}")
            runtime = time.time() - t0
            self.state.consumed_seconds += runtime

            self.repo.log_step(
                self.experiment_id, self.state.iteration, decision.action.value, h,
                tool=metrics.get("tool"), params=decision.params, metrics=metrics,
                interpretation=interpretation, failure=failure, runtime_seconds=runtime)
            self.state.record(decision.action, decision.params, metrics, failure)

            if failure is not None:
                self.event("  recording failure in scientific memory and stopping this branch.")
                break

            self.state.step_count += 1
            if decision.action in EXPERIMENT_ACTIONS:
                self.state.iteration += 1

        manifest = self.build_manifest(time.time() - t_start)
        (self.dirs["metrics"] / "experiment_manifest.json").write_text(json.dumps(manifest, indent=2))
        (self.dir / "agent_log.txt").write_text("\n".join(self.log))
        self.repo.finish_experiment(self.experiment_id, manifest)
        self.event(f"Experiment complete. Artifacts: {self.dir}")
        return manifest

    # ------------------------------------------------------------------
    def execute(self, decision: Decision) -> tuple[dict, str]:
        action = decision.action

        if action is Action.GENERATE_ENSEMBLE:
            ens = engines.generate_ensemble(self.cfg, self.dirs["structures"])
            if not ens.structures:
                raise RuntimeError("Ensemble provider returned zero structures")
            for p in ens.structures:
                if not Path(p).exists() or Path(p).stat().st_size == 0:
                    raise RuntimeError(f"Ensemble artifact missing/empty: {p}")
            self.state.ensemble = [str(p) for p in ens.structures]
            self.state.ensemble_provider = ens.provider
            self.state.ensemble_is_equilibrium = ens.is_equilibrium_sample
            self.state.baseline_structure = str(ens.structures[0])
            self.ensemble_provider = ens
            self.event(f"  provider={ens.provider} n_states={len(ens.structures)} "
                       f"equilibrium_sample={ens.is_equilibrium_sample}")
            return ({"tool": "bioemu_tool", "n_states": len(ens.structures),
                     "provider": ens.provider, "is_equilibrium_sample": ens.is_equilibrium_sample},
                    f"Obtained {len(ens.structures)} conformational states via {ens.provider}. "
                    f"{'' if ens.is_equilibrium_sample else 'These are NOT equilibrium samples; state populations are uniform placeholders.'}")

        if action is Action.ANALYZE_ENSEMBLE:
            a = engines.analyze_ensemble([Path(p) for p in self.state.ensemble])
            self.state.conformational_states = {
                "pdb_ids": a.pdb_ids,
                "n_common_residues": len(a.common_resids),
                "max_rmsf": float(a.rmsf_per_residue.max()),
                "mean_rmsf": float(a.rmsf_per_residue.mean()),
                "max_pairwise_rmsd": float(a.rmsd_matrix.max()),
                "pca_explained_variance_ratio": [float(x) for x in a.pca_explained_variance_ratio],
                "pca_coords": a.pca_coords.tolist(),
                "rmsd_matrix": a.rmsd_matrix.tolist(),
                "rmsf_per_residue": a.rmsf_per_residue.tolist(),
                "common_resids": a.common_resids,
                "dimensionality_reduction": "PCA (TICA invalid: ensemble has no temporal ordering)",
            }
            (self.dirs["metrics"] / "ensemble_analysis.json").write_text(
                json.dumps(self.state.conformational_states, indent=2))
            self.event(f"  max_RMSF={a.rmsf_per_residue.max():.2f} A, "
                       f"max_pairwise_RMSD={a.rmsd_matrix.max():.2f} A, "
                       f"PC1_var={a.pca_explained_variance_ratio[0]:.3f}")
            return ({"tool": "structural_analysis", "max_rmsf": float(a.rmsf_per_residue.max()),
                     "max_pairwise_rmsd": float(a.rmsd_matrix.max())},
                    f"Ensemble spans {a.rmsd_matrix.max():.2f} A maximum pairwise Ca RMSD with peak "
                    f"per-residue RMSF {a.rmsf_per_residue.max():.2f} A, indicating a mobile region "
                    "worth probing for transient cavities.")

        if action is Action.FIND_POCKETS:
            gt = self.cfg["target"]["ground_truth_pocket_residues"]
            candidates = engines.find_pockets([Path(p) for p in self.state.ensemble],
                                               self.dirs["pockets"], gt)
            if not candidates:
                raise RuntimeError("fpocket returned no pocket candidates")
            self.state.pocket_candidates = candidates
            self.repo.save_pockets(self.experiment_id, candidates)
            self.event(f"  {len(candidates)} raw pocket candidates across "
                       f"{len(self.state.ensemble)} states")
            return ({"tool": "mdpocket_tool(fpocket)", "n_candidates": len(candidates)},
                    f"Detected {len(candidates)} candidate cavities across the ensemble using "
                    "fpocket Voronoi tessellation.")

        if action is Action.SELECT_POCKET:
            baseline_id = self.cfg["target"]["baseline_pdb_id"]
            families = engines.cluster_pocket_families(
                self.state.pocket_candidates, baseline_id, len(self.state.ensemble))
            ranked = engines.rank_pocket_families(families)
            self.state.evidence["pocket_families"] = [{
                "representative": f["representative"].key,
                "persistence": f["persistence"],
                "baseline_volume": f["baseline_volume"],
                "max_volume": f["max_volume"],
                "cryptic_volume_gain": f["cryptic_volume_gain"],
                "novelty": f["novelty"],
                "max_druggability": f["max_druggability"],
                "rank_score": f["rank_score"],
                "ground_truth_overlap": f["representative"].ground_truth_overlap,
            } for f in ranked]
            (self.dirs["metrics"] / "pocket_families.json").write_text(
                json.dumps(self.state.evidence["pocket_families"], indent=2))

            top_family = ranked[0]
            top = top_family["representative"]
            self.state.selected_pocket = top
            self.state.evidence["selected_family"] = self.state.evidence["pocket_families"][0]
            self.event(f"  {len(ranked)} cross-state pocket families; selected {top.key}")
            self.event(f"  volume={top.volume:.1f}A^3 druggability={top.druggability:.2f} "
                       f"novelty={top_family['novelty']:.2f} "
                       f"persistence={top_family['persistence']:.2f} "
                       f"baseline_volume={top_family['baseline_volume']:.1f}A^3")
            self.event(f"  post-hoc ground-truth overlap={top.ground_truth_overlap:.2f} "
                       "(NOT used in ranking)")
            return ({"tool": "rank_pocket_families", "selected": top.key,
                     "n_families": len(ranked), "volume": top.volume,
                     "druggability": top.druggability, "novelty": top_family["novelty"],
                     "persistence": top_family["persistence"],
                     "cryptic_volume_gain": top_family["cryptic_volume_gain"],
                     "ground_truth_overlap": top.ground_truth_overlap},
                    f"Grouped {len(self.state.pocket_candidates)} per-state cavities into "
                    f"{len(ranked)} cross-state families and ranked them blind on "
                    f"druggability + novelty-vs-baseline + volume. Committed to {top.key}: "
                    f"{top.volume:.0f} A^3, druggability {top.druggability:.2f}, novelty "
                    f"{top_family['novelty']:.2f} (baseline volume "
                    f"{top_family['baseline_volume']:.0f} A^3), present in "
                    f"{top_family['n_states_present']}/{len(self.state.ensemble)} states. "
                    f"Post-hoc, it covers {top.ground_truth_overlap*100:.0f}% of documented "
                    "Switch-II ground-truth residues.")

        if action is Action.SCREEN_LIGANDS:
            pairs = engines.load_ligand_library(self.ligand_csv)
            validated = engines.validate_ligands(pairs, self.dirs["ligands"])
            n_valid = sum(1 for v in validated if v.valid and v.embedded_3d)
            self.event(f"  {n_valid}/{len(validated)} ligands passed RDKit validation + 3D embedding")
            results = engines.dock_ligands(
                validated, Path(self.state.ensemble[self._state_index(self.state.selected_pocket)]),
                self.state.selected_pocket, self.dirs["pockets"], self.dirs["docking"],
                exhaustiveness=self.cfg["docking"]["exhaustiveness"],
                n_poses=self.cfg["docking"]["num_modes"], on_event=self.event)
            if not results:
                raise RuntimeError("No ligand docked successfully")
            self.state.docking_results.extend(results)
            self._score_and_persist(results, origin="library")
            best = min(r["best_affinity_kcal"] for r in results)
            return ({"tool": "docking_tool(vina)", "n_docked": len(results), "best_affinity": best},
                    f"Screened {len(results)} validated ligands into {self.state.selected_pocket.key}; "
                    f"best Vina affinity {best:.2f} kcal/mol. Computational score only -- not "
                    "evidence of binding.")

        if action is Action.OPTIMIZE_LIGAND:
            best_record = min(self.state.docking_results, key=lambda r: r["best_affinity_kcal"])
            optimizer, mode = reinvent_tool.get_optimizer(prefer_reinvent=True)
            self.optimizer_mode = mode
            self.event(f"  optimizer mode: {mode} (seed={best_record['ligand_name']})")
            analogs = optimizer.generate(best_record["ligand_name"], best_record["smiles"])
            if not analogs:
                raise RuntimeError("Optimizer produced no valid analogs")
            self.event(f"  generated {len(analogs)} validated analogs")

            validated = engines.validate_ligands([(a.name, a.smiles) for a in analogs],
                                                  self.dirs["ligands"])
            results = engines.dock_ligands(
                validated, Path(self.state.ensemble[self._state_index(self.state.selected_pocket)]),
                self.state.selected_pocket, self.dirs["pockets"], self.dirs["docking"],
                exhaustiveness=self.cfg["docking"]["exhaustiveness"],
                n_poses=self.cfg["docking"]["num_modes"], on_event=self.event)
            if not results:
                raise RuntimeError("No analog docked successfully")
            for r in results:
                r["parent"] = best_record["ligand_name"]
            self.state.docking_results.extend(results)
            self.state.optimization_results.extend(results)
            self._score_and_persist(results, origin="optimized")

            best_analog = min(results, key=lambda r: r["best_affinity_kcal"])
            delta = best_analog["best_affinity_kcal"] - best_record["best_affinity_kcal"]
            improved = delta < 0
            self.event(f"  best analog {best_analog['ligand_name']}: "
                       f"{best_analog['best_affinity_kcal']:.2f} kcal/mol "
                       f"(delta {delta:+.2f} vs seed)")
            return ({"tool": f"reinvent_tool({mode})", "n_analogs": len(results),
                     "seed": best_record["ligand_name"], "delta_kcal": delta, "improved": improved},
                    f"Generated {len(analogs)} analogs of {best_record['ligand_name']} via {mode} "
                    f"and re-docked them. Best analog changed affinity by {delta:+.2f} kcal/mol "
                    f"({'improvement' if improved else 'no improvement'}).")

        raise NotImplementedError(f"Action {action} has no executor")

    # ------------------------------------------------------------------
    def _state_index(self, pocket) -> int:
        for i, p in enumerate(self.state.ensemble):
            if Path(p).stem == pocket.state_pdb_id:
                return i
        raise RuntimeError(f"Receptor structure for pocket {pocket.key} not found in ensemble")

    def _score_and_persist(self, results: list[dict], origin: str):
        """Discovery Score is computed ONLY here, by deterministic code."""
        ref = self.cfg["discovery_score"]["affinity_reference"]
        best_aff, worst_aff = ref["best_kcal"], ref["worst_kcal"]
        max_vol = max((c.volume for c in self.state.pocket_candidates), default=1.0)
        pocket = self.state.selected_pocket

        family = self.state.evidence.get("selected_family", {})
        # True cryptic novelty (volume absent from the apo baseline), not druggability.
        pocket_novelty = family.get("novelty", 0.0)
        # State novelty: a cavity present in few states is the transient one.
        state_frequency = family.get("persistence", 1.0 / max(len(self.state.ensemble), 1))

        for r in results:
            ds = compute_discovery_score(DiscoveryScoreInputs(
                pocket_volume=pocket.volume, max_observed_volume=max_vol,
                pocket_druggability=pocket_novelty,
                state_frequency=state_frequency,
                binding_affinity_kcal=r["best_affinity_kcal"],
                best_possible_affinity_kcal=best_aff, worst_possible_affinity_kcal=worst_aff,
                ligand_qed=r.get("qed") or 0.0,
                lipinski_violations=r.get("lipinski_violations") or 0,
            ), self.weights)
            r["discovery_score"] = ds["discovery_score"]
            r["discovery_components"] = ds["components"]
            r["origin"] = origin
            self.state.best_discovery_score = max(self.state.best_discovery_score, ds["discovery_score"])

        self.repo.save_docking(self.experiment_id, results)

    # ------------------------------------------------------------------
    def build_manifest(self, runtime: float) -> dict:
        ranked = sorted((r for r in self.state.docking_results if "discovery_score" in r),
                        key=lambda r: r["discovery_score"], reverse=True)
        return {
            "experiment_id": self.experiment_id,
            "mode": self.mode,
            "target": self.state.target,
            "reproducibility": self.provenance(),
            "ensemble": {
                "provider": self.state.ensemble_provider,
                "n_states": len(self.state.ensemble),
                "is_equilibrium_sample": self.state.ensemble_is_equilibrium,
                "structures": self.state.ensemble,
            },
            "ensemble_analysis": {k: v for k, v in self.state.conformational_states.items()
                                   if k not in ("pca_coords", "rmsd_matrix", "rmsf_per_residue", "common_resids")},
            "n_pocket_candidates": len(self.state.pocket_candidates),
            "selected_pocket": (self.state.selected_pocket.__dict__
                                 if self.state.selected_pocket else None),
            "optimizer_mode": self.optimizer_mode,
            "n_ligands_docked": len(self.state.docking_results),
            "best_discovery_score": (None if self.state.best_discovery_score == float("-inf")
                                      else self.state.best_discovery_score),
            "ranked_results": ranked,
            "closed_loop_iterations_executed": self.state.iteration,
            "total_agent_steps": self.state.step_count,
            "agent_decisions": self.state.history,
            "runtime_seconds": round(runtime, 1),
            "compute_budget_seconds": self.state.budget_seconds,
            "tools_executed": ["RDKit", "fpocket 4.0", "AutoDock Vina 1.2.7", "MDAnalysis", "OpenBabel 3.1.0"],
            "tools_unavailable_fallback_used": self._fallbacks_used(),
        }

    def _fallbacks_used(self) -> dict:
        """Only report a fallback for a tool this run actually invoked and that
        actually fell back -- do not claim a tool "was used instead" of one this
        closed loop never attempted to call."""
        fallbacks = {}
        if self.state.ensemble_provider and self.state.ensemble_provider != "bioemu":
            reason = "no CUDA GPU detected" if not cuda_available() else "BioEmu unavailable"
            fallbacks["BioEmu"] = f"{reason}; {self.state.ensemble_provider} used instead"
        if self.optimizer_mode and self.optimizer_mode != "reinvent4":
            fallbacks["REINVENT4"] = f"not installed; optimizer fallback used: {self.optimizer_mode}"
        return fallbacks

    def provenance(self) -> dict:
        try:
            commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                     text=True, timeout=10).stdout.strip()
        except Exception:  # noqa: BLE001
            commit = "unknown"
        import numpy
        import rdkit
        return {
            "git_commit": commit,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": numpy.__version__,
            "rdkit": rdkit.__version__,
            "cuda": "available" if cuda_available() else "unavailable (no nvidia-smi on this host)",
            "random_seeds": {"vina": 42, "rdkit_embed": 42},
            "config_file": str(self.config_path),
            "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }


if __name__ == "__main__":
    agent = ClosedLoopAgent(Path("configs/kras_g12d.yaml"), Path("data/ligands_kras.csv"))
    m = agent.run()
    print(json.dumps({k: v for k, v in m.items()
                      if k not in ("ranked_results", "agent_decisions")}, indent=2))
