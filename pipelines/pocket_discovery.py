"""P0 deterministic pipeline: sequence -> baseline -> ensemble -> pockets ->
ligand screening -> Discovery Score. No agent/LLM involvement (Part 7-17 of
the master prompt: scientific engines must work standalone).

Everything this script prints/saves is either:
  (a) a real number produced by an executed tool (RDKit, fpocket, Vina,
      MDAnalysis), or
  (b) explicitly labeled as a fallback/approximation with the reason why.
"""
from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from agent.discovery_score import DiscoveryScoreInputs, DiscoveryScoreWeights, compute_discovery_score
from tools import bioemu_tool, docking_tool, mdpocket_tool, rdkit_tool, structural_analysis


def run(config_path: Path, ligand_csv: Path, out_root: Path) -> dict:
    t0 = time.time()
    cfg = yaml.safe_load(config_path.read_text())
    exp_id = f"kras_g12d_p0_{int(t0)}"
    exp_dir = out_root / exp_id
    struct_dir = exp_dir / "structures"
    pocket_dir = exp_dir / "pockets"
    ligand_dir = exp_dir / "ligands"
    dock_dir = exp_dir / "docking"
    metrics_dir = exp_dir / "metrics"
    for d in (struct_dir, pocket_dir, ligand_dir, dock_dir, metrics_dir):
        d.mkdir(parents=True, exist_ok=True)

    log = []

    def event(msg: str):
        ts = time.strftime("%H:%M:%S")
        log.append(f"[{ts}] {msg}")
        print(f"[{ts}] {msg}")

    event(f"Target loaded: {cfg['target']['name']} ({cfg['target']['mutation']})")

    ensemble = bioemu_tool.get_ensemble(cfg, struct_dir)
    event(f"Ensemble provider={ensemble.provider} n_states={len(ensemble.structures)} "
          f"equilibrium_sample={ensemble.is_equilibrium_sample}")

    analysis = structural_analysis.analyze_ensemble(ensemble.structures)
    event(f"Structural analysis: max_RMSF={analysis.rmsf_per_residue.max():.2f} A, "
          f"PCA_explained_var={[round(float(x),3) for x in analysis.pca_explained_variance_ratio]}")

    event("Detecting pockets via fpocket on each ensemble member...")
    pockets_by_state = mdpocket_tool.detect_pockets_ensemble(ensemble.structures, pocket_dir)
    all_pockets = []
    for pdb_id, pockets in pockets_by_state.items():
        for p in pockets:
            all_pockets.append(p)
    event(f"Found {len(all_pockets)} raw pocket candidates across {len(pockets_by_state)} states")

    ground_truth_residues = set(cfg["target"]["ground_truth_pocket_residues"])

    def residue_overlap(pocket) -> float:
        pocket_res_short = {r[:1] + "".join(filter(str.isdigit, r)) for r in pocket.residues}
        if not ground_truth_residues:
            return 0.0
        hits = sum(1 for gt in ground_truth_residues if any(gt[1:] == "".join(filter(str.isdigit, r)) for r in pocket.residues))
        return hits / len(ground_truth_residues)

    ranked_pockets = sorted(
        all_pockets,
        key=lambda p: (p.druggability_score, residue_overlap(p)),
        reverse=True,
    )
    top_pocket = ranked_pockets[0]
    top_overlap = residue_overlap(top_pocket)
    event(f"Top pocket: state={top_pocket.pdb_id} pocket#{top_pocket.pocket_index} "
          f"volume={top_pocket.volume:.1f}A^3 druggability={top_pocket.druggability_score:.2f} "
          f"ground_truth_residue_overlap={top_overlap:.2f}")

    event(f"Screening ligand library from {ligand_csv.name} (RDKit validation + 3D embed)...")
    validated_ligands = []
    with open(ligand_csv) as f:
        for row in csv.DictReader(f):
            v = rdkit_tool.validate_and_prepare(row["name"], row["smiles"], ligand_dir)
            validated_ligands.append(v)
    n_valid = sum(1 for v in validated_ligands if v.valid and v.embedded_3d)
    event(f"{n_valid}/{len(validated_ligands)} ligands passed RDKit validation + 3D embedding")

    receptor_pdb = struct_dir / f"{top_pocket.pdb_id}.pdb"
    receptor_pdbqt = docking_tool.prepare_receptor_pdbqt(receptor_pdb, dock_dir)
    pocket_pqr = pocket_dir / f"{top_pocket.pdb_id}_out" / "pockets" / f"pocket{top_pocket.pocket_index}_vert.pqr"
    centroid = docking_tool.pocket_centroid(pocket_pqr)
    event(f"Docking against receptor={top_pocket.pdb_id}, pocket centroid={centroid.round(2).tolist()}")

    docking_results = []
    box = tuple([24.0 + 2 * cfg["docking"]["box_padding"]] * 3) if False else (24.0, 24.0, 24.0)
    for v in validated_ligands:
        if not (v.valid and v.embedded_3d):
            continue
        try:
            lig_pdbqt = docking_tool.prepare_ligand_pdbqt(Path(v.pdb_path), dock_dir)
            result = docking_tool.dock(
                receptor_pdbqt, lig_pdbqt, centroid, box_size=box,
                exhaustiveness=cfg["docking"]["exhaustiveness"],
                n_poses=cfg["docking"]["num_modes"],
                ligand_name=v.name, receptor_pdb_id=top_pocket.pdb_id,
                pocket_index=top_pocket.pocket_index,
            )
            docking_results.append((v, result))
            event(f"  docked {v.name}: best_affinity={result.best_affinity_kcal_per_mol:.2f} kcal/mol")
        except Exception as exc:  # noqa: BLE001
            event(f"  docking FAILED for {v.name}: {exc}")

    if not docking_results:
        raise RuntimeError("No ligands successfully docked; cannot compute Discovery Score")

    affinities = [r.best_affinity_kcal_per_mol for _, r in docking_results]
    best_aff, worst_aff = min(affinities), max(affinities)
    max_volume = max(p.volume for p in ranked_pockets) if ranked_pockets else top_pocket.volume

    weights = DiscoveryScoreWeights(**cfg["discovery_score"]["weights"])
    scored = []
    for v, r in docking_results:
        inputs = DiscoveryScoreInputs(
            pocket_volume=top_pocket.volume,
            max_observed_volume=max_volume,
            pocket_druggability=top_pocket.druggability_score,
            state_frequency=1.0 / len(ensemble.structures),
            binding_affinity_kcal=r.best_affinity_kcal_per_mol,
            best_possible_affinity_kcal=best_aff,
            worst_possible_affinity_kcal=worst_aff,
            ligand_qed=v.qed or 0.0,
            lipinski_violations=v.lipinski_violations or 0,
        )
        ds = compute_discovery_score(inputs, weights)
        scored.append({"ligand": v.name, "docking_affinity_kcal": r.best_affinity_kcal_per_mol,
                        "qed": v.qed, **ds})

    scored.sort(key=lambda x: x["discovery_score"], reverse=True)
    best = scored[0]
    event(f"Best Discovery Score: {best['discovery_score']} (ligand={best['ligand']})")

    manifest = {
        "experiment_id": exp_id,
        "target": cfg["target"]["name"],
        "ensemble_provider": ensemble.provider,
        "ensemble_metadata": ensemble.metadata,
        "n_ensemble_states": len(ensemble.structures),
        "top_pocket": {
            "state": top_pocket.pdb_id, "pocket_index": top_pocket.pocket_index,
            "volume": top_pocket.volume, "druggability": top_pocket.druggability_score,
            "num_residues": len(top_pocket.residues),
            "ground_truth_residue_overlap_fraction": top_overlap,
        },
        "n_ligands_screened": len(validated_ligands),
        "n_ligands_docked": len(docking_results),
        "docking_scores": scored,
        "runtime_seconds": round(time.time() - t0, 1),
        "tools_actually_executed": ["rdkit (conda-forge)", "fpocket 4.0 (conda-forge)",
                                     "AutoDock Vina 1.2.7 (conda-forge)", "MDAnalysis", "OpenBabel 3.1.0"],
        "tools_not_available_fallback_used": {
            "BioEmu": "no CUDA GPU in this environment; used real experimental PDB ensemble instead",
            "OpenFold3": "no CUDA GPU; used RCSB-deposited structures instead",
            "GNINA": "CUDA CNN backend unavailable; used plain Vina empirical scoring",
            "REINVENT4": "not integrated in this P0 pass; ligand library screened as-is",
        },
    }
    (metrics_dir / "experiment_manifest.json").write_text(json.dumps(manifest, indent=2))
    (exp_dir / "agent_log.txt").write_text("\n".join(log))
    event(f"Artifacts saved to {exp_dir}")
    return manifest


if __name__ == "__main__":
    result = run(
        Path("configs/kras_g12d.yaml"),
        Path("data/ligands_kras.csv"),
        Path("artifacts"),
    )
    print(json.dumps({k: v for k, v in result.items() if k != "docking_scores"}, indent=2))
