#!/usr/bin/env bash
# One-shot GPU rental session for ConforM-Agent.
#
# Run this ONCE on a rented CUDA machine (Vast.ai / RunPod RTX 4090
# recommended -- see docs/DEPENDENCIES.md) to do EVERY compute-requiring
# step in a single sitting:
#   1. Install real BioEmu and generate a genuine equilibrium ensemble
#      from the apo KRAS G12D sequence (no crystal structures as input).
#   2. Run the full closed-loop pipeline against that real ensemble.
#   3. Install GNINA (prebuilt binary, no build) and CNN-rescore the top
#      Vina hits from that run.
#   4. Package everything into one tarball to pull back to the laptop.
#
# Usage:
#   git clone <repo> && cd ConforM
#   bash scripts/gpu_session.sh
#
# Expected total GPU time for KRAS G12D (169 residues, 1000 samples):
# a few minutes for BioEmu (official benchmark: ~4 min/1000 samples at
# 100 residues on an A100; 4090 is comparable) + docking/GNINA overhead.
# Budget ~1-2 hours wall-clock including environment setup; well inside a
# 200 CNY (~$30) rental budget even on a per-hour billed instance.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "================================================================"
echo "ConforM-Agent GPU session"
echo "================================================================"

echo "[0/6] Verifying CUDA GPU is present..."
if ! nvidia-smi >/dev/null 2>&1; then
    echo "ERROR: nvidia-smi not found. This script must run on a CUDA GPU instance."
    exit 1
fi
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

echo "[1/6] Setting up Python environment..."
if command -v conda >/dev/null 2>&1; then
    conda env create -f environment.yml -n conform 2>/dev/null || true
    eval "$(conda shell.bash hook)"
    conda activate conform
else
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -q -r <(python3 - <<'PY'
import yaml
env = yaml.safe_load(open("environment.yml"))
for dep in env.get("dependencies", []):
    if isinstance(dep, str) and dep not in ("python",) and "=" not in dep:
        print(dep)
PY
)
fi
export PYTHONPATH=.

echo "[2/6] Installing BioEmu (CUDA build)..."
pip install -q "bioemu[cuda]"
python -c "from tools.bioemu_tool import cuda_available; assert cuda_available(), 'CUDA not detected by bioemu_tool'"
echo "  CUDA confirmed available to the pipeline."

echo "[3/6] Enabling BioEmu in configs/kras_g12d.yaml..."
python3 - <<'PY'
import yaml
from pathlib import Path
p = Path("configs/kras_g12d.yaml")
cfg = yaml.safe_load(p.read_text())
cfg["ensemble"]["provider"] = "bioemu"
cfg["ensemble"]["bioemu"]["enabled"] = True
cfg["ensemble"]["bioemu"].setdefault("num_samples", 1000)
p.write_text(yaml.dump(cfg, sort_keys=False))
print("  bioemu.enabled=True, num_samples=", cfg["ensemble"]["bioemu"]["num_samples"])
PY

echo "[4/6] Running the closed-loop pipeline against a REAL BioEmu ensemble..."
echo "  (this is the single most important run in the whole project: it is"
echo "   the first genuine test of generative sampling from the apo sequence,"
echo "   with no ligand-bound structure anywhere in the input.)"
python scripts/run_experiment.py run --target kras-g12d --closed-loop 2>&1 | tee /tmp/gpu_session_run.log

EXPERIMENT_DIR=$(python3 - <<'PY'
from pathlib import Path
from scripts.generate_report import latest_experiment
print(latest_experiment())
PY
)
echo "  Experiment: $EXPERIMENT_DIR"

echo "[5/6] Installing GNINA (prebuilt binary) and CNN-rescoring top hits..."
if [ ! -f ./gnina ]; then
    wget -q https://github.com/gnina/gnina/releases/download/v1.3.2/gnina -O ./gnina
    chmod +x ./gnina
fi
export PATH="$PWD:$PATH"
python scripts/gnina_rescore.py --experiment "$EXPERIMENT_DIR" --top-n 5 \
    2>&1 | tee /tmp/gpu_session_gnina.log || \
    echo "  GNINA rescoring failed or unavailable -- Vina-only results are still valid and complete."

echo "[6/6] Generating report and packaging results..."
python scripts/generate_report.py --experiment "$EXPERIMENT_DIR"
python evaluation/ablations.py 2>&1 | tee /tmp/gpu_session_ablations.log

TARBALL="conform_gpu_session_$(date +%Y%m%d_%H%M%S).tar.gz"
tar czf "$TARBALL" \
    "$EXPERIMENT_DIR" \
    artifacts/ablation_report.json \
    /tmp/gpu_session_run.log \
    /tmp/gpu_session_gnina.log \
    /tmp/gpu_session_ablations.log \
    configs/kras_g12d.yaml

echo
echo "================================================================"
echo "DONE. Pull this file back to your laptop and stop the instance:"
echo "  $TARBALL"
echo "================================================================"
