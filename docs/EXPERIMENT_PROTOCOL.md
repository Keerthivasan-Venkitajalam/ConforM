# Experiment Protocol

## Reproducing the KRAS G12D experiment

```bash
conda env create -f environment.yml
conda activate conform
export PYTHONPATH=.

python scripts/run_experiment.py validate                    # environment check
python scripts/run_experiment.py run --target kras-g12d --closed-loop
python scripts/run_experiment.py report
python scripts/run_experiment.py dashboard
```

Full validation (executes every real tool, ~5 min):
```bash
./validate_e2e.sh
```

Benchmark across all baselines/ablations (~10 min):
```bash
python evaluation/ablations.py
python evaluation/ablations.py --mode static     # single mode
```

## Protocol parameters (configs/kras_g12d.yaml)

| Stage | Setting |
|---|---|
| Ensemble | 4 experimental KRAS G12D structures: 4DST (baseline/apo reference), 5US4, 7RPZ, 5XCO |
| Structural analysis | Kabsch superposition on common Cα atoms; PCA (not TICA) |
| Pocket detection | fpocket 4.0, default parameters, run per structure |
| Cross-state grouping | residue-set Jaccard ≥ 0.50 |
| Pocket ranking | 0.40·druggability + 0.40·novelty + 0.20·normalized volume |
| Ligand prep | RDKit sanitize → ETKDG embed (seed 42) → MMFF optimize (500 iters) |
| Docking | AutoDock Vina 1.2.7, box 24×24×24 Å at pocket alpha-sphere centroid, exhaustiveness 8, 5 poses, seed 42, 4 CPU |
| Affinity normalization | absolute anchors −12.0 → 1.0, −4.0 → 0.0 kcal/mol |
| Optimization | RDKit R-group enumeration (7 groups on aromatic C–H), max 12 analogs, MW ≤ 600, ≤1 Lipinski violation |
| Agent | MAX_ITERATIONS 5 (experiment cycles), stop_score 0.85, hard cap 20 total steps |

## Determinism
Fixed seeds: Vina 42, RDKit embedding 42, random baseline 42. Re-running
produces identical pocket selection and near-identical affinities (Vina's
Monte Carlo search is seeded but small numerical variation across runs of
fpocket volume calculation has been observed, e.g. 813–836 Å³ for the same
pocket, because fpocket alpha-sphere placement is sensitive to input atom
ordering).

## Artifacts produced
```
artifacts/<experiment_id>/
    structures/   real PDB files
    pockets/      fpocket output (_info.txt, pockets/*.pqr, *_atm.pdb)
    ligands/      RDKit 3D-embedded PDB per ligand
    docking/      PDBQT receptor/ligand inputs
    metrics/      experiment_manifest.json, ensemble_analysis.json, pocket_families.json
    report/       experiment_report.html
    agent_log.txt
artifacts/conform_memory.db     SQLite scientific memory (all experiments)
artifacts/ablation_report.json  benchmark comparison
```

`experiment_manifest.json` carries the reproducibility block: git commit,
Python/NumPy/RDKit versions, platform, seeds, config path, UTC timestamp.

## Reproducing the generalization tests (ABL kinase, PRMT5)

```bash
python scripts/compute_ground_truth.py data/structures/1IEP.pdb STI --chain A   # ABL DFG-out contacts
python scripts/compute_ground_truth.py data/structures/6UXX.pdb QL1             # PRMT5 EE-loop contacts

python scripts/run_experiment.py run --target abl-kinase --closed-loop
python scripts/run_experiment.py run --target prmt5 --closed-loop
```

`tests/test_frozen_config.py` enforces that `configs/abl_kinase.yaml` and
`configs/prmt5.yaml` share byte-identical `pocket_detection` / `docking` /
`discovery_score` / `agent` sections with `configs/kras_g12d.yaml`, and that
all three targets use the same ligand library — the generalization claim
fails CI if anyone tunes a per-target weight. See
[GENERALIZATION.md](GENERALIZATION.md) for results and interpretation.
