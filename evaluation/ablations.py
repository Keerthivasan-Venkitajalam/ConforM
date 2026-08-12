"""Runs every mode and produces a real comparison table.

Usage:
    python evaluation/ablations.py                    # all modes
    python evaluation/ablations.py --mode static      # one mode
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from evaluation.baselines import MODES, run_mode
from evaluation.cryptic_recovery import evaluate_manifest
from evaluation.metrics import summarize_run

ACTIVES = {"Hypothetical_piperazine_pyrimidine_scaffold"}


def run_all(config_path: Path, ligand_csv: Path, out_root: Path, modes: list[str]) -> dict:
    cfg = yaml.safe_load(config_path.read_text())
    gt = cfg["target"]["ground_truth_pocket_residues"]
    rows = []
    for mode in modes:
        print(f"\n{'='*70}\nMODE: {mode}\n{'='*70}", flush=True)
        t0 = time.time()
        try:
            manifest = run_mode(mode, config_path, ligand_csv, out_root)
            summary = summarize_run(manifest, gt, actives=ACTIVES)
            summary["cryptic_recovery"] = evaluate_manifest(manifest, gt)
            summary["failed"] = False
        except Exception as exc:  # noqa: BLE001
            print(f"MODE {mode} FAILED: {type(exc).__name__}: {exc}", flush=True)
            summary = {"mode": mode, "failed": True,
                       "error": f"{type(exc).__name__}: {exc}",
                       "runtime_seconds": round(time.time() - t0, 1)}
        rows.append(summary)

    out_root.mkdir(parents=True, exist_ok=True)
    report = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "target": cfg["target"]["name"],
        "ground_truth_residues": gt,
        "modes": rows,
        "omitted_baseline": {
            "name": "Baseline 3 (1 microsecond classical MD)",
            "reason": "requires ~10,000 GPU-hours by the research plan's own estimate; "
                      "not run, and deliberately not approximated",
        },
    }
    path = out_root / "ablation_report.json"
    path.write_text(json.dumps(report, indent=2))
    print_table(rows)
    print(f"\nWrote {path}")
    return report


def print_table(rows: list[dict]):
    header = (f"{'mode':<24}{'states':>7}{'pockets':>9}{'recall':>8}"
              f"{'best_kcal':>11}{'discovery':>11}{'iters':>7}{'sec':>8}")
    print("\n" + "=" * len(header))
    print("ABLATION / BASELINE COMPARISON (all numbers from real executed runs)")
    print("=" * len(header))
    print(header)
    print("-" * len(header))
    for r in rows:
        if r.get("failed"):
            print(f"{r['mode']:<24}{'FAILED':>7}  {r.get('error','')[:50]}")
            continue
        print(f"{r['mode']:<24}{r.get('n_states') or 0:>7}"
              f"{r.get('n_pocket_candidates') or 0:>9}"
              f"{r.get('cryptic_residue_recall') or 0:>8.2f}"
              f"{r.get('best_affinity_kcal') or 0:>11.2f}"
              f"{r.get('best_discovery_score') or 0:>11.3f}"
              f"{r.get('closed_loop_iterations') or 0:>7}"
              f"{r.get('runtime_seconds') or 0:>8.0f}")
    print("=" * len(header))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=MODES, help="run a single mode")
    parser.add_argument("--config", default="configs/kras_g12d.yaml")
    parser.add_argument("--ligands", default="data/ligands_kras.csv")
    parser.add_argument("--out", default="artifacts")
    args = parser.parse_args()
    run_all(Path(args.config), Path(args.ligands), Path(args.out),
            [args.mode] if args.mode else MODES)
