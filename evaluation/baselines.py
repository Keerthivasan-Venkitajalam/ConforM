"""Baseline modes for the Cryptic Recovery Challenge.

Modes (master prompt #26):
  static                 single baseline (apo) structure + docking; no ensemble
  random                 ensemble, but a RANDOMLY chosen pocket (no ranking logic)
  no-pocket-guidance     ensemble, but dock into EVERY state's top cavity uniformly
  no-ligand-optimization ConforM-Agent with the optimization step disabled
  conform-agent          full system

Baseline 3 from the research plan (1 microsecond classical MD) is NOT
implemented: it requires ~10,000 GPU-hours by the plan's own estimate and is
infeasible here. Its absence is reported rather than approximated, so no
fabricated MD numbers enter the comparison.

Every mode below runs the SAME real tools (fpocket + Vina + RDKit) so the
comparison is apples-to-apples; only the agent's decision logic differs.
"""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from agent.discovery_score import DiscoveryScoreInputs, DiscoveryScoreWeights, compute_discovery_score
from agent.loop_controller import ClosedLoopAgent
from agent.policies import PolicyConfig
from db.repository import Repository
from pipelines import engines

MODES = ["static", "random", "no-pocket-guidance", "no-ligand-optimization", "conform-agent"]


def _score_all(results, pocket, families_entry, cfg, n_states, max_vol):
    ref = cfg["discovery_score"]["affinity_reference"]
    weights = DiscoveryScoreWeights(**cfg["discovery_score"]["weights"])
    for r in results:
        ds = compute_discovery_score(DiscoveryScoreInputs(
            pocket_volume=pocket.volume, max_observed_volume=max_vol,
            pocket_druggability=families_entry.get("novelty", 0.0),
            state_frequency=families_entry.get("persistence", 1.0 / max(n_states, 1)),
            binding_affinity_kcal=r["best_affinity_kcal"],
            best_possible_affinity_kcal=ref["best_kcal"],
            worst_possible_affinity_kcal=ref["worst_kcal"],
            ligand_qed=r.get("qed") or 0.0,
            lipinski_violations=r.get("lipinski_violations") or 0,
        ), weights)
        r["discovery_score"] = ds["discovery_score"]
        r["discovery_components"] = ds["components"]
    return results


def run_simple_mode(mode: str, config_path: Path, ligand_csv: Path, out_root: Path,
                    seed: int = 42) -> dict:
    """Runs the non-agentic baselines with the same real tooling."""
    t0 = time.time()
    cfg = yaml.safe_load(Path(config_path).read_text())
    exp_id = f"{cfg['target']['name'].lower()}_{mode}_{int(t0)}"
    exp_dir = out_root / exp_id
    dirs = {n: exp_dir / n for n in ("structures", "pockets", "ligands", "docking", "metrics")}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    log = []

    def event(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        log.append(line)
        print(line, flush=True)

    baseline_id = cfg["target"]["baseline_pdb_id"]
    gt = cfg["target"]["ground_truth_pocket_residues"]

    if mode == "static":
        # Only the apo baseline structure: no conformational sampling at all.
        # Must also force bioemu off, not just override ensemble_pdb_ids --
        # ensemble_pdb_ids only feeds the ExperimentalEnsembleProvider fallback,
        # and get_ensemble() tries BioEmu first whenever ensemble.bioemu.enabled
        # is true, which it inherits from cfg. Without this, "static" silently
        # samples the full real BioEmu ensemble like every other mode, defeating
        # the entire point of a no-sampling control (see RESEARCH_CORRECTIONS.md).
        sub_cfg = {
            **cfg,
            "target": {**cfg["target"], "ensemble_pdb_ids": [baseline_id]},
            "ensemble": {**cfg["ensemble"],
                         "bioemu": {**cfg["ensemble"].get("bioemu", {}), "enabled": False}},
        }
        ens = engines.generate_ensemble(sub_cfg, dirs["structures"])
        event(f"STATIC baseline: single structure {baseline_id}, no ensemble generated")
    else:
        ens = engines.generate_ensemble(cfg, dirs["structures"])
        event(f"{mode}: ensemble of {len(ens.structures)} states")

    candidates = engines.find_pockets(list(ens.structures), dirs["pockets"], gt)
    event(f"{len(candidates)} raw pocket candidates")
    families = engines.cluster_pocket_families(candidates, baseline_id, len(ens.structures))

    rng = random.Random(seed)
    if mode == "random":
        chosen_family = rng.choice(families)
        event(f"RANDOM baseline: pocket chosen at random (seed={seed}), no ranking applied")
    elif mode == "no-pocket-guidance":
        # No geometric guidance: take the first pocket of the first state.
        chosen_family = next(f for f in families
                             if f["representative"].state_pdb_id == Path(ens.structures[0]).stem)
        event("NO-POCKET-GUIDANCE: first cavity of first state, no volumetric ranking")
    else:
        chosen_family = engines.rank_pocket_families(families)[0]
        event(f"Ranked {len(families)} families on druggability+novelty+volume")

    pocket = chosen_family["representative"]
    event(f"selected {pocket.key} vol={pocket.volume:.1f} drug={pocket.druggability:.2f} "
          f"novelty={chosen_family['novelty']:.2f} "
          f"(post-hoc ground-truth overlap={pocket.ground_truth_overlap:.2f})")

    pairs = engines.load_ligand_library(ligand_csv)
    validated = engines.validate_ligands(pairs, dirs["ligands"])
    receptor = next(p for p in ens.structures if Path(p).stem == pocket.state_pdb_id)
    results = engines.dock_ligands(validated, Path(receptor), pocket, dirs["pockets"],
                                    dirs["docking"], exhaustiveness=cfg["docking"]["exhaustiveness"],
                                    n_poses=cfg["docking"]["num_modes"], on_event=event)
    max_vol = max((c.volume for c in candidates), default=1.0)
    results = _score_all(results, pocket, chosen_family, cfg, len(ens.structures), max_vol)
    ranked = sorted(results, key=lambda r: r["discovery_score"], reverse=True)

    manifest = {
        "experiment_id": exp_id, "mode": mode, "target": cfg["target"]["name"],
        "ensemble": {"provider": ens.provider, "n_states": len(ens.structures),
                      "is_equilibrium_sample": ens.is_equilibrium_sample,
                      "structures": [str(p) for p in ens.structures]},
        "n_pocket_candidates": len(candidates),
        "selected_pocket": pocket.__dict__,
        "selected_pocket_novelty": chosen_family["novelty"],
        "selected_pocket_persistence": chosen_family["persistence"],
        "n_ligands_docked": len(results),
        "best_discovery_score": max((r["discovery_score"] for r in ranked), default=None),
        "ranked_results": ranked,
        "closed_loop_iterations_executed": 1,
        "optimizer_mode": None,
        "runtime_seconds": round(time.time() - t0, 1),
    }
    (dirs["metrics"] / "experiment_manifest.json").write_text(json.dumps(manifest, indent=2))
    (exp_dir / "agent_log.txt").write_text("\n".join(log))
    return manifest


def run_mode(mode: str, config_path: Path, ligand_csv: Path, out_root: Path = Path("artifacts"),
             repo: Repository | None = None) -> dict:
    if mode in ("static", "random", "no-pocket-guidance"):
        return run_simple_mode(mode, config_path, ligand_csv, out_root)

    agent = ClosedLoopAgent(config_path, ligand_csv, out_root, mode=mode, repo=repo)
    if mode == "no-ligand-optimization":
        # Ablate the optimization step by stopping once library evidence exists.
        agent.policy.config = PolicyConfig(max_iterations=1, stop_score=-999.0)
    return agent.run()
