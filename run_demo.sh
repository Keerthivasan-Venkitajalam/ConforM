#!/usr/bin/env bash
# One-command ConforM-Agent demo: verifies deps, runs a real lightweight
# closed-loop experiment on KRAS G12D, generates the report, launches the
# dashboard.
#
# This demo runs REAL tools (fpocket, AutoDock Vina, RDKit) on REAL
# experimental structures. It uses documented fallbacks for GPU-only models
# (BioEmu, OpenFold3, GNINA, REINVENT4) and labels them as such everywhere.
# It does not ship precomputed results presented as live inference.
set -euo pipefail
cd "$(dirname "$0")"

PY="${CONFORM_PYTHON:-python}"
export PYTHONPATH="${PYTHONPATH:-.}"

echo "========================================"
echo "ConforM-Agent demo"
echo "========================================"
echo

echo "[1/5] Verifying dependencies..."
if ! $PY scripts/run_experiment.py validate; then
    echo
    echo "Environment incomplete. Create it with:"
    echo "    conda env create -f environment.yml && conda activate conform"
    exit 1
fi

echo
echo "[2/5] Running closed-loop experiment on KRAS G12D..."
echo "      (downloads real PDB structures on first run; needs network)"
$PY scripts/run_experiment.py run --target kras-g12d --closed-loop

echo
echo "[3/5] Evaluating cryptic-pocket recovery..."
$PY - <<'PYEOF'
import json
from pathlib import Path
import yaml
from evaluation.cryptic_recovery import evaluate_manifest
from scripts.generate_report import latest_experiment
cfg = yaml.safe_load(Path("configs/kras_g12d.yaml").read_text())
m = json.loads((latest_experiment() / "metrics" / "experiment_manifest.json").read_text())
print(json.dumps(evaluate_manifest(m, cfg["target"]["ground_truth_pocket_residues"]), indent=2))
PYEOF

echo
echo "[4/5] Generating scientific report..."
$PY scripts/run_experiment.py report

echo
echo "[5/5] Launching Streamlit dashboard (Ctrl-C to stop)..."
exec $PY -m streamlit run visualization/dashboard.py
