"""conform-agent CLI.

    python scripts/run_experiment.py setup
    python scripts/run_experiment.py validate
    python scripts/run_experiment.py run --target kras-g12d
    python scripts/run_experiment.py run --target kras-g12d --closed-loop
    python scripts/run_experiment.py benchmark --target kras-g12d
    python scripts/run_experiment.py ablate --mode static
    python scripts/run_experiment.py report [--experiment ID]
    python scripts/run_experiment.py dashboard
    python scripts/run_experiment.py history
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

TARGETS = {"kras-g12d": Path("configs/kras_g12d.yaml")}
DEFAULT_LIGANDS = Path("data/ligands_kras.csv")

REQUIRED_BINARIES = ["fpocket", "obabel"]
REQUIRED_MODULES = ["rdkit", "vina", "MDAnalysis", "numpy", "yaml"]


def cmd_setup(_args):
    print("ConforM-Agent setup\n" + "=" * 40)
    print("Create the scientific environment with:\n")
    print("    conda env create -f environment.yml")
    print("    conda activate conform\n")
    print("Then verify with:  python scripts/run_experiment.py validate")
    return 0


def cmd_validate(_args):
    print("Validating environment\n" + "=" * 40)
    ok = True
    for binary in REQUIRED_BINARIES:
        found = shutil.which(binary)
        print(f"  {'PASS' if found else 'FAIL'}  binary {binary:10s} {found or 'NOT FOUND'}")
        ok &= bool(found)
    for module in REQUIRED_MODULES:
        try:
            __import__(module)
            print(f"  PASS  module {module}")
        except ImportError as exc:
            print(f"  FAIL  module {module}: {exc}")
            ok = False
    try:
        subprocess.run(["nvidia-smi"], capture_output=True, timeout=10, check=True)
        print("  INFO  CUDA GPU detected -- BioEmu/GNINA integration could be enabled")
    except Exception:
        print("  INFO  no CUDA GPU: BioEmu/OpenFold3/GNINA/DiffDock fallbacks will be used")
    print("\nSTATUS:", "READY" if ok else "INCOMPLETE")
    return 0 if ok else 1


def cmd_run(args):
    from agent.loop_controller import ClosedLoopAgent
    from evaluation.baselines import run_simple_mode

    config = TARGETS[args.target]
    ligands = Path(args.ligands)
    if args.closed_loop:
        agent = ClosedLoopAgent(config, ligands, Path(args.out), mode="conform-agent",
                                 max_iterations=args.max_iterations)
        manifest = agent.run()
    else:
        manifest = run_simple_mode("conform-agent", config, ligands, Path(args.out))
    print(f"\nExperiment: {manifest['experiment_id']}")
    print(f"Best Discovery Score: {manifest.get('best_discovery_score')}")
    return 0


def cmd_benchmark(args):
    from evaluation.ablations import run_all
    from evaluation.baselines import MODES

    run_all(TARGETS[args.target], Path(args.ligands), Path(args.out), MODES)
    return 0


def cmd_ablate(args):
    from evaluation.ablations import run_all
    from evaluation.baselines import MODES

    run_all(TARGETS[args.target], Path(args.ligands), Path(args.out),
            [args.mode] if args.mode else MODES)
    return 0


def cmd_report(args):
    from scripts.generate_report import build_report, latest_experiment

    root = Path(args.out)
    d = (Path(args.experiment) if args.experiment and Path(args.experiment).exists()
         else (root / args.experiment) if args.experiment else latest_experiment(root))
    print(f"Report written to: {build_report(d)}")
    return 0


def cmd_dashboard(_args):
    return subprocess.call([sys.executable, "-m", "streamlit", "run", "visualization/dashboard.py"])


def cmd_history(args):
    from db.repository import Repository

    repo = Repository()
    for exp in repo.list_experiments():
        print(f"{exp['created_at'][:19]}  {exp['mode']:<22} {exp['experiment_id']:<45} {exp['status']}")
        if args.verbose:
            for s in repo.steps(exp["experiment_id"]):
                status = s["failure"] or "ok"
                print(f"    iter {s['iteration']}  {s['action']:<20} {status}")
    return 0


def main():
    ap = argparse.ArgumentParser(prog="conform-agent", description="ConforM-Agent CLI")
    sub = ap.add_subparsers(dest="command", required=True)

    def common(p):
        p.add_argument("--target", default="kras-g12d", choices=list(TARGETS))
        p.add_argument("--ligands", default=str(DEFAULT_LIGANDS))
        p.add_argument("--out", default="artifacts")

    sub.add_parser("setup").set_defaults(func=cmd_setup)
    sub.add_parser("validate").set_defaults(func=cmd_validate)

    p = sub.add_parser("run"); common(p)
    p.add_argument("--closed-loop", action="store_true", help="run the autonomous agent loop")
    p.add_argument("--samples", type=int, help="(reserved; BioEmu sample count when a GPU is available)")
    p.add_argument("--max-iterations", type=int, default=None)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("benchmark"); common(p); p.set_defaults(func=cmd_benchmark)

    p = sub.add_parser("ablate"); common(p)
    p.add_argument("--mode", default=None)
    p.set_defaults(func=cmd_ablate)

    p = sub.add_parser("report")
    p.add_argument("--experiment", default=None)
    p.add_argument("--out", default="artifacts")
    p.set_defaults(func=cmd_report)

    sub.add_parser("dashboard").set_defaults(func=cmd_dashboard)

    p = sub.add_parser("history")
    p.add_argument("--verbose", action="store_true")
    p.set_defaults(func=cmd_history)

    args = ap.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
