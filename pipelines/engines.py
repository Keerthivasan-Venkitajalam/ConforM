"""Scientific engines: independently callable, deterministic, no LLM involved.

Each function performs one real scientific operation and returns structured
results. The agent (agent/loop_controller.py) orchestrates these; it never
computes scientific quantities itself.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from agent.state import LigandRecord, PocketCandidate
from tools import bioemu_tool, docking_tool, mdpocket_tool, rdkit_tool, structural_analysis


def residue_number(res: str) -> str:
    return "".join(filter(str.isdigit, res))


def ground_truth_overlap(residues: list[str], ground_truth: list[str]) -> float:
    """Fraction of documented ground-truth pocket residues present in this pocket."""
    if not ground_truth:
        return 0.0
    pocket_numbers = {residue_number(r) for r in residues}
    hits = sum(1 for gt in ground_truth if residue_number(gt) in pocket_numbers)
    return hits / len(ground_truth)


def generate_ensemble(cfg: dict, out_dir: Path):
    return bioemu_tool.get_ensemble(cfg, out_dir)


def analyze_ensemble(structures: list[Path]):
    return structural_analysis.analyze_ensemble(structures)


def find_pockets(structures: list[Path], work_dir: Path, ground_truth: list[str],
                 min_volume: float = 0.0) -> list[PocketCandidate]:
    by_state = mdpocket_tool.detect_pockets_ensemble(structures, work_dir)
    candidates = []
    for _pdb_id, pockets in by_state.items():
        for p in pockets:
            if p.volume < min_volume:
                continue
            candidates.append(PocketCandidate(
                state_pdb_id=p.pdb_id, pocket_index=p.pocket_index, volume=p.volume,
                druggability=p.druggability_score, residues=p.residues,
                ground_truth_overlap=ground_truth_overlap(p.residues, ground_truth),
            ))
    return candidates


def _jaccard(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def cluster_pocket_families(candidates: list[PocketCandidate], baseline_pdb_id: str,
                            n_states: int, similarity_threshold: float = 0.5) -> list[dict]:
    """Group per-state pockets into cross-state 'families' by residue overlap.

    This is the ensemble-level analysis that true `mdpocket` would do via
    volumetric grid persistence; since we run fpocket per structure (see
    docs/RESEARCH_CORRECTIONS.md), we establish cross-state correspondence
    by Jaccard similarity of the lining-residue sets instead.

    Per family we compute the quantities the Discovery Score actually needs:
      persistence         fraction of ensemble states where the cavity is present
      baseline_volume     its volume in the apo/baseline structure (0.0 if absent)
      cryptic_volume_gain max_volume - baseline_volume  (Part 8 of the plan:
                          "volume that does not exist in the apo baseline")
      novelty             cryptic_volume_gain / max_volume, in [0, 1];
                          1.0 == cavity entirely absent from the baseline
    """
    families: list[list[PocketCandidate]] = []
    for c in sorted(candidates, key=lambda x: x.volume, reverse=True):
        res_c = {residue_number(r) for r in c.residues}
        placed = False
        for fam in families:
            res_f = {residue_number(r) for r in fam[0].residues}
            if _jaccard(res_c, res_f) >= similarity_threshold:
                fam.append(c)
                placed = True
                break
        if not placed:
            families.append([c])

    summaries = []
    for fam in families:
        states_present = {m.state_pdb_id for m in fam}
        baseline_members = [m for m in fam if m.state_pdb_id == baseline_pdb_id]
        baseline_volume = max((m.volume for m in baseline_members), default=0.0)
        best = max(fam, key=lambda m: (m.druggability, m.volume))
        max_volume = max(m.volume for m in fam)
        gain = max(max_volume - baseline_volume, 0.0)
        summaries.append({
            "representative": best,
            "members": fam,
            "n_states_present": len(states_present),
            "persistence": len(states_present) / max(n_states, 1),
            "baseline_volume": baseline_volume,
            "max_volume": max_volume,
            "cryptic_volume_gain": gain,
            "novelty": (gain / max_volume) if max_volume > 0 else 0.0,
            "max_druggability": max(m.druggability for m in fam),
        })
    return summaries


def rank_pocket_families(families: list[dict], weights: dict | None = None) -> list[dict]:
    """Deterministic ranking targeting CRYPTIC cavities.

    Ground-truth residue overlap is deliberately NOT an input -- that would
    leak the answer. Ranking uses only blind descriptors: druggability,
    novelty vs. the apo baseline, and normalized volume. Transient pockets
    (low persistence, high novelty) are rewarded, which is the entire point
    of the system: a cavity that is always open in the baseline is not a
    cryptic-pocket discovery.
    """
    w = weights or {"druggability": 0.40, "novelty": 0.40, "volume": 0.20}
    max_vol = max((f["max_volume"] for f in families), default=1.0) or 1.0

    def score(f: dict) -> float:
        return (w["druggability"] * f["max_druggability"]
                + w["novelty"] * f["novelty"]
                + w["volume"] * (f["max_volume"] / max_vol))

    for f in families:
        f["rank_score"] = round(score(f), 4)
    return sorted(families, key=lambda f: f["rank_score"], reverse=True)


def load_ligand_library(csv_path: Path) -> list[tuple[str, str]]:
    with open(csv_path) as f:
        return [(row["name"], row["smiles"]) for row in csv.DictReader(f)]


def validate_ligands(pairs: list[tuple[str, str]], out_dir: Path):
    return [rdkit_tool.validate_and_prepare(name, smiles, out_dir) for name, smiles in pairs]


def dock_ligands(validated, receptor_pdb: Path, pocket: PocketCandidate, pocket_dir: Path,
                 work_dir: Path, exhaustiveness: int = 8, n_poses: int = 5,
                 box_size=(24.0, 24.0, 24.0), on_event=None) -> list[dict]:
    receptor_pdbqt = docking_tool.prepare_receptor_pdbqt(receptor_pdb, work_dir)
    pqr = pocket_dir / f"{pocket.state_pdb_id}_out" / "pockets" / f"pocket{pocket.pocket_index}_vert.pqr"
    centroid = docking_tool.pocket_centroid(pqr)

    results = []
    for v in validated:
        if not (v.valid and v.embedded_3d):
            if on_event:
                on_event(f"  skipped {v.name}: failed RDKit validation ({v.error})")
            continue
        try:
            lig_pdbqt = docking_tool.prepare_ligand_pdbqt(Path(v.pdb_path), work_dir)
            r = docking_tool.dock(
                receptor_pdbqt, lig_pdbqt, centroid, box_size=box_size,
                exhaustiveness=exhaustiveness, n_poses=n_poses, ligand_name=v.name,
                receptor_pdb_id=pocket.state_pdb_id, pocket_index=pocket.pocket_index,
            )
            results.append({
                "ligand_name": v.name, "smiles": v.canonical_smiles,
                "receptor_pdb_id": pocket.state_pdb_id, "pocket_index": pocket.pocket_index,
                "engine": "vina", "best_affinity_kcal": r.best_affinity_kcal_per_mol,
                "poses": r.poses_affinity_kcal_per_mol, "qed": v.qed,
                "lipinski_violations": v.lipinski_violations,
            })
            if on_event:
                on_event(f"  docked {v.name}: {r.best_affinity_kcal_per_mol:.2f} kcal/mol")
        except Exception as exc:  # noqa: BLE001
            if on_event:
                on_event(f"  docking FAILED for {v.name}: {exc}")
    return results


def to_ligand_records(docking_results: list[dict], origin: str = "library") -> list[LigandRecord]:
    return [LigandRecord(
        name=d["ligand_name"], smiles=d.get("smiles", ""), qed=d.get("qed"),
        lipinski_violations=d.get("lipinski_violations"),
        best_affinity_kcal=d.get("best_affinity_kcal"), origin=origin,
        parent=d.get("parent"),
    ) for d in docking_results]
