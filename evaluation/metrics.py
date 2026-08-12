"""Objective evaluation metrics. All deterministic, all computed from real
pipeline outputs -- never from an LLM's summary of them."""
from __future__ import annotations


def residue_number(res: str) -> str:
    return "".join(filter(str.isdigit, res))


def pocket_residue_recovery(predicted_residues: list[str], ground_truth_residues: list[str]) -> dict:
    """Recall/precision/Jaccard of a predicted pocket's lining residues against
    the documented ground-truth cryptic-site residues.

    NOTE: this is a residue-level proxy for the Discretized Volume Overlap
    (DVO) named in the research plan. True DVO requires voxelizing both the
    predicted cavity grid and the experimental ligand-occupied volume; we do
    not compute it here (see docs/LIMITATIONS.md), and this metric is reported
    under its own name so it is not mistaken for DVO.
    """
    pred = {residue_number(r) for r in predicted_residues}
    truth = {residue_number(r) for r in ground_truth_residues}
    if not truth:
        return {"recall": 0.0, "precision": 0.0, "jaccard": 0.0, "n_hits": 0}
    hits = pred & truth
    return {
        "recall": len(hits) / len(truth),
        "precision": len(hits) / len(pred) if pred else 0.0,
        "jaccard": len(hits) / len(pred | truth) if (pred | truth) else 0.0,
        "n_hits": len(hits),
        "hit_residues": sorted(hits, key=int),
    }


def enrichment(docking_results: list[dict], actives: set[str], top_n: int = 3) -> dict:
    """Fraction of designated 'active' ligands recovered in the top-N by affinity.

    With a 10-ligand toy library this is a smoke-test statistic, not a
    statistically meaningful enrichment factor. Reported with n so the reader
    can judge.
    """
    ranked = sorted((r for r in docking_results if r.get("best_affinity_kcal") is not None),
                    key=lambda r: r["best_affinity_kcal"])
    top = {r["ligand_name"] for r in ranked[:top_n]}
    recovered = top & actives
    return {
        "top_n": top_n,
        "n_ligands": len(ranked),
        "n_actives": len(actives),
        "actives_in_top_n": len(recovered),
        "recovered": sorted(recovered),
        "enrichment_ratio": (
            (len(recovered) / top_n) / (len(actives) / len(ranked))
            if ranked and actives and top_n else None
        ),
        "caveat": "toy library; not a statistically powered enrichment factor",
    }


def best_affinity(docking_results: list[dict]) -> float | None:
    vals = [r["best_affinity_kcal"] for r in docking_results
            if r.get("best_affinity_kcal") is not None]
    return min(vals) if vals else None


def summarize_run(manifest: dict, ground_truth_residues: list[str],
                  actives: set[str] | None = None) -> dict:
    pocket = manifest.get("selected_pocket") or {}
    results = manifest.get("ranked_results", [])
    recovery = pocket_residue_recovery(pocket.get("residues", []), ground_truth_residues)
    out = {
        "mode": manifest.get("mode"),
        "experiment_id": manifest.get("experiment_id"),
        "ensemble_provider": manifest.get("ensemble", {}).get("provider"),
        "n_states": manifest.get("ensemble", {}).get("n_states"),
        "n_pocket_candidates": manifest.get("n_pocket_candidates"),
        "selected_pocket": pocket.get("state_pdb_id"),
        "selected_pocket_volume": pocket.get("volume"),
        "selected_pocket_druggability": pocket.get("druggability"),
        "cryptic_residue_recall": recovery["recall"],
        "cryptic_residue_precision": recovery["precision"],
        "cryptic_residue_jaccard": recovery["jaccard"],
        "n_ligands_docked": manifest.get("n_ligands_docked"),
        "best_affinity_kcal": best_affinity(results),
        "best_discovery_score": manifest.get("best_discovery_score"),
        "closed_loop_iterations": manifest.get("closed_loop_iterations_executed"),
        "runtime_seconds": manifest.get("runtime_seconds"),
        "gpu_hours": 0.0,
        "cpu_seconds": manifest.get("runtime_seconds"),
    }
    if actives:
        out["enrichment"] = enrichment(results, actives)
    return out
