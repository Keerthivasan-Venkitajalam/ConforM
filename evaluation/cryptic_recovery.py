"""Cryptic pocket recovery evaluation against documented ground truth.

CryptoBench integration status: the CryptoBench dataset (skrhakv/CryptoBench,
MIT) is NOT downloaded or evaluated in this build. Doing it properly means
running the full pipeline over ~1,107 apo-holo pairs, which is far beyond the
compute available here, and a partial run would produce a number that looks
like a benchmark result without being one. Instead this module evaluates
recovery against the explicitly documented KRAS G12D Switch-II ground truth
(PDB 7RPZ / MRTX1133 site) declared in configs/kras_g12d.yaml.

`--subset` / `--target` style scoping over CryptoBench is left as documented
future work in docs/BENCHMARKS.md.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from evaluation.metrics import pocket_residue_recovery


def evaluate_manifest(manifest: dict, ground_truth_residues: list[str]) -> dict:
    pocket = manifest.get("selected_pocket") or {}
    recovery = pocket_residue_recovery(pocket.get("residues", []), ground_truth_residues)
    return {
        "mode": manifest.get("mode"),
        "selected_pocket": f"{pocket.get('state_pdb_id')}:pocket{pocket.get('pocket_index')}",
        "ground_truth_residues": ground_truth_residues,
        "recovered_residues": recovery["hit_residues"],
        "recall": recovery["recall"],
        "precision": recovery["precision"],
        "jaccard": recovery["jaccard"],
        "recovered_cryptic_site": recovery["recall"] >= 0.6,
        "metric_note": "residue-level recovery, not Discretized Volume Overlap (see metrics.py)",
    }


def main(manifest_path: Path, config_path: Path = Path("configs/kras_g12d.yaml")):
    cfg = yaml.safe_load(config_path.read_text())
    manifest = json.loads(Path(manifest_path).read_text())
    result = evaluate_manifest(manifest, cfg["target"]["ground_truth_pocket_residues"])
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    main(Path(sys.argv[1]))
