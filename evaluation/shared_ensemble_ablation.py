"""Re-runs the ablation table with every non-static mode sharing ONE real
BioEmu ensemble, to eliminate independent-sampling variance as a confound.

Background (see docs/RESEARCH_CORRECTIONS.md #8): the original ablation table
had each mode independently call BioEmu, so random/no-pocket-guidance/
no-ligand-optimization/conform-agent each evaluated a DIFFERENT real
generative draw. That confounds any cross-mode comparison of which pocket
gets selected -- a recall difference between two modes could be genuine
algorithmic behavior, or could just be two different draws of the diffusion
model. This script removes that confound by reusing the same, already
GPU-verified 99-state ensemble (kras_g12d_conform-agent_1787945235) for every
mode. `static` is intentionally excluded from the shared ensemble -- it is
supposed to be a true single-structure control, and forcing it onto a
99-state ensemble would defeat that purpose (that was bug #3 in
RESEARCH_CORRECTIONS.md #8, already fixed separately).

No new BioEmu sampling is performed by this script -- it reuses the
already-generated, already GPU-verified structures on disk. Only the
deterministic pocket-detection/ranking/docking steps are re-run per mode.

Usage:
    python evaluation/shared_ensemble_ablation.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from evaluation.baselines import run_mode, run_simple_mode
from evaluation.cryptic_recovery import evaluate_manifest
from evaluation.metrics import summarize_run

MASTER_ENSEMBLE_DIR = Path(
    "artifacts/kras_g12d_conform-agent_1787945235/structures/bioemu_run/frames")
SHARED_MODES = ["random", "no-pocket-guidance", "no-ligand-optimization", "conform-agent"]
ACTIVES = {"Hypothetical_piperazine_pyrimidine_scaffold"}


def print_table(rows: list[dict]):
    header = (f"{'mode':<24}{'states':>7}{'pockets':>9}{'recall':>8}"
              f"{'best_kcal':>11}{'discovery':>11}{'iters':>7}{'sec':>8}")
    print("\n" + "=" * len(header))
    print("SHARED-ENSEMBLE ABLATION (identical structural input across all non-static modes)")
    print("=" * len(header))
    print(header)
    print("-" * len(header))
    for r in rows:
        if r.get("failed"):
            print(f"{r['mode']:<24}{'FAILED':>7}  {r.get('error', '')[:50]}")
            continue
        print(f"{r['mode']:<24}{r.get('n_states') or 0:>7}"
              f"{r.get('n_pocket_candidates') or 0:>9}"
              f"{r.get('cryptic_residue_recall') or 0:>8.2f}"
              f"{r.get('best_affinity_kcal') or 0:>11.2f}"
              f"{r.get('best_discovery_score') or 0:>11.3f}"
              f"{r.get('closed_loop_iterations') or 0:>7}"
              f"{r.get('runtime_seconds') or 0:>8.0f}")
    print("=" * len(header))


def main():
    config_path = Path("configs/kras_g12d.yaml")
    ligand_csv = Path("data/ligands_kras.csv")
    out_root = Path("artifacts")
    cfg = yaml.safe_load(config_path.read_text())
    gt = cfg["target"]["ground_truth_pocket_residues"]

    frames = sorted(MASTER_ENSEMBLE_DIR.glob("frame_*.pdb"))
    if not frames:
        raise SystemExit(f"No frames found under {MASTER_ENSEMBLE_DIR} -- "
                          "run the real BioEmu experiment first (see RESEARCH_CORRECTIONS.md #8).")
    print(f"Master shared ensemble: {len(frames)} real BioEmu states from "
          f"{MASTER_ENSEMBLE_DIR}\n")

    rows = []

    # static stays a true single-structure control -- NOT part of the shared ensemble.
    print("=" * 70, "\nMODE: static (true single-structure control, unchanged)\n", "=" * 70,
          flush=True)
    manifest = run_mode("static", config_path, ligand_csv, out_root)
    summary = summarize_run(manifest, gt, actives=ACTIVES)
    summary["cryptic_recovery"] = evaluate_manifest(manifest, gt)
    summary["failed"] = False
    rows.append(summary)

    for mode in SHARED_MODES:
        print(f"\n{'=' * 70}\nMODE: {mode} (shared master ensemble)\n{'=' * 70}", flush=True)
        t0 = time.time()
        try:
            manifest = run_mode(mode, config_path, ligand_csv, out_root,
                                 shared_structures=frames)
            summary = summarize_run(manifest, gt, actives=ACTIVES)
            summary["cryptic_recovery"] = evaluate_manifest(manifest, gt)
            summary["failed"] = False
        except Exception as exc:  # noqa: BLE001
            print(f"MODE {mode} FAILED: {type(exc).__name__}: {exc}", flush=True)
            summary = {"mode": mode, "failed": True,
                       "error": f"{type(exc).__name__}: {exc}",
                       "runtime_seconds": round(time.time() - t0, 1)}
        rows.append(summary)

    report = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "target": cfg["target"]["name"],
        "ground_truth_residues": gt,
        "shared_ensemble_source": str(MASTER_ENSEMBLE_DIR),
        "shared_ensemble_n_states": len(frames),
        "modes": rows,
        "note": ("static intentionally excluded from the shared ensemble -- it is a "
                 "true single-structure control by design, not a bug."),
    }
    path = out_root / "ablation_report_shared_ensemble.json"
    path.write_text(json.dumps(report, indent=2))
    print_table(rows)
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main()
