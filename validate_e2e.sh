#!/usr/bin/env bash
# ConforM-Agent end-to-end validation.
# Prints PASS only when the corresponding check actually succeeded.
set -uo pipefail
cd "$(dirname "$0")"

PY="${CONFORM_PYTHON:-python}"
declare -a NAMES=() STATES=()
FAILED=0

check() {  # check <name> <command...>
    local name="$1"; shift
    if "$@" >/tmp/conform_e2e_step.log 2>&1; then
        NAMES+=("$name"); STATES+=("PASS")
    else
        NAMES+=("$name"); STATES+=("FAIL")
        FAILED=1
        echo "--- $name FAILED ---"
        tail -15 /tmp/conform_e2e_step.log
    fi
}

echo "Running ConforM-Agent E2E validation (this executes real tools; takes a few minutes)..."
echo

check "Environment"       $PY scripts/run_experiment.py validate
check "Imports"           $PY -c "import agent.loop_controller, pipelines.engines, evaluation.ablations, visualization.plots, scripts.generate_report"
check "Unit tests"        $PY -m pytest tests/ -q
check "Structure"         $PY -c "
from pathlib import Path
from tools.structure_tool import StructureProvider
r = StructureProvider().get_baseline('4DST', Path('data/structures'))
assert r.path.exists() and r.path.stat().st_size > 1000
"
check "Ensemble"          $PY -c "
import yaml
from pathlib import Path
from tools.bioemu_tool import get_ensemble
cfg = yaml.safe_load(Path('configs/kras_g12d.yaml').read_text())
e = get_ensemble(cfg, Path('data/structures'))
assert len(e.structures) >= 2
"
check "Structural analysis" $PY -c "
import yaml
from pathlib import Path
from tools.bioemu_tool import get_ensemble
from tools.structural_analysis import analyze_ensemble
cfg = yaml.safe_load(Path('configs/kras_g12d.yaml').read_text())
e = get_ensemble(cfg, Path('data/structures'))
a = analyze_ensemble(e.structures)
assert a.rmsd_matrix.max() > 0 and len(a.common_resids) > 50
"
check "Pocket detection"  $PY -c "
from pathlib import Path
from tools.mdpocket_tool import run_fpocket
p = run_fpocket(Path('data/structures/7RPZ.pdb'), Path('artifacts/_e2e_pockets'))
assert len(p) > 0 and p[0].volume > 0
"
check "RDKit"             $PY -c "
from tools.rdkit_tool import validate_and_prepare
r = validate_and_prepare('benzene','c1ccccc1')
assert r.valid and r.qed > 0
assert not validate_and_prepare('junk','!!!not_smiles(((').valid
"
check "Ligand optimizer"  $PY -c "
from tools.reinvent_tool import get_optimizer
o, mode = get_optimizer(prefer_reinvent=True)
a = o.generate('seed','COc1ccc(N2CCNCC2)cc1Nc1ncccn1')
assert len(a) > 0, 'no analogs generated'
print('optimizer mode:', mode)
"
check "Discovery Score"   $PY -c "
from agent.discovery_score import DiscoveryScoreInputs, compute_discovery_score
i = DiscoveryScoreInputs(800,1000,0.9,0.25,-9.6,-12.0,-4.0,0.6,0)
a, b = compute_discovery_score(i), compute_discovery_score(i)
assert a == b, 'Discovery Score is not deterministic'
"
check "Memory / DB"       $PY -c "
from db.repository import Repository
r = Repository()
r.create_experiment('e2e_probe','KRAS','test',{})
r.log_step('e2e_probe',0,'FIND_POCKETS','hash_probe')
assert r.has_completed('e2e_probe','hash_probe'), 'duplicate detection broken'
assert not r.has_completed('e2e_probe','other_hash')
"
check "Docking"           $PY -c "
from pathlib import Path
import numpy as np
from tools import docking_tool
from tools.rdkit_tool import validate_and_prepare
w = Path('artifacts/_e2e_dock'); w.mkdir(parents=True, exist_ok=True)
rec = docking_tool.prepare_receptor_pdbqt(Path('data/structures/7RPZ.pdb'), w)
v = validate_and_prepare('benzene','c1ccccc1', w)
lig = docking_tool.prepare_ligand_pdbqt(Path(v.pdb_path), w)
c = docking_tool.pocket_centroid(Path('artifacts/_e2e_pockets/7RPZ_out/pockets/pocket1_vert.pqr'))
r = docking_tool.dock(rec, lig, c, exhaustiveness=2, n_poses=2)
assert r.best_affinity_kcal_per_mol < 0, 'implausible docking score'
"
check "Closed loop"       $PY -c "
from pathlib import Path
from agent.loop_controller import ClosedLoopAgent
a = ClosedLoopAgent(Path('configs/kras_g12d.yaml'), Path('data/ligands_e2e.csv'),
                    Path('artifacts'), mode='e2e-validation')
m = a.run()
assert m['closed_loop_iterations_executed'] >= 2, 'closed loop did not complete 2 iterations'
assert m['selected_pocket'] is not None
assert m['best_discovery_score'] is not None
"
check "Evaluation"        $PY -c "
import json
from pathlib import Path
from evaluation.cryptic_recovery import evaluate_manifest
from evaluation.metrics import summarize_run
from scripts.generate_report import latest_experiment
m = json.loads((latest_experiment() / 'metrics' / 'experiment_manifest.json').read_text())
gt = ['H95','Y96','Q99','V9','D69']
s = summarize_run(m, gt); e = evaluate_manifest(m, gt)
assert 'cryptic_residue_recall' in s and 'recall' in e
"
check "Report"            $PY -c "
from scripts.generate_report import build_report, latest_experiment
p = build_report(latest_experiment())
assert p.exists() and p.stat().st_size > 5000
"
check "Dashboard"         $PY -c "
import ast, pathlib
ast.parse(pathlib.Path('visualization/dashboard.py').read_text())
import streamlit, visualization.molecular_viewer, visualization.plots
"

echo
echo "========================================"
echo "ConforM-Agent E2E VALIDATION"
echo "========================================"
echo
for i in "${!NAMES[@]}"; do
    printf "%-22s %s\n" "${NAMES[$i]}" "${STATES[$i]}"
done
echo
if [ "$FAILED" -eq 0 ]; then
    echo "STATUS: READY"
else
    echo "STATUS: FAILED"
fi
echo "========================================"
exit $FAILED
