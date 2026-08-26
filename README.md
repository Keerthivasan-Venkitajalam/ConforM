# ConforM-Agent

**A closed-loop scientific agent that explores protein conformational state
space to find transient (cryptic) binding pockets and evaluate ligand
hypotheses against them.**

Demonstration target: **KRAS G12D**, validated against the documented
Switch-II cryptic pocket.

---

## ⚠️ Scientific honesty notice

This is a **computational hypothesis-generation platform**. It does not
demonstrate that any molecule binds KRAS G12D.

- A favorable docking score means *"this ligand produced a favorable
  computational docking score under the evaluated protocol"* — nothing more.
- Recovering the Switch-II pocket is **known ground-truth recovery**, which
  validates the method. It is **not** a novel biological discovery: that
  pocket is well characterized and already drugged (MRTX1133).
- Several models from the original design (BioEmu, OpenFold3, GNINA,
  DiffDock-Pocket, REINVENT 4) **could not be executed on this CPU-only
  machine**. Documented fallbacks were used and are labeled as such in every
  manifest, report, and dashboard view. See [docs/LIMITATIONS.md](docs/LIMITATIONS.md).

---

## What it does

```
KRAS G12D
   ↓  structure acquisition        (OpenFold3 → ESMFold → RCSB)      [RCSB used]
   ↓  conformational ensemble      (BioEmu → experimental ensemble)  [experimental used]
   ↓  structural analysis          RMSD / RMSF / PCA
   ↓  pocket detection             (mdpocket → fpocket → P2Rank)     [fpocket used]
   ↓  cross-state pocket families  persistence + novelty vs. apo baseline
   ↓  blind pocket ranking         druggability + novelty + volume
   ↓  ligand validation            RDKit sanitize / 3D embed / Lipinski / QED
   ↓  docking                      (GNINA → Vina)                    [Vina used]
   ↓  Discovery Score              deterministic Python, never LLM-generated
   ↓  agent decision               screen? optimize? stop?
   ↺  next experiment
```

The agent decides what to run next from an explicit state machine
([agent/policies.py](agent/policies.py)) — it is not a hard-coded linear
pipeline. Scientific engines are independently callable without any agent or
LLM.

## Key result (real executed run, apo-only ensemble)

Blind ranking — with ground-truth overlap explicitly **excluded** from the
objective — selects a cavity that is:

- **absent from the apo baseline structure** (novelty 1.00, baseline volume 0 Å³)
- **present in only 1 of 4 sampled states** (persistence 0.25 — genuinely transient)
- **partially covering the documented Switch-II ground-truth residues** (H95, Y96 recovered; Q99, V9, D69 not — recall 0.40)

The `static` baseline, given only the apo structure, instead selects the
always-open nucleotide site: ground-truth recall **0.00**. This is the
central experimental claim, and it is reproducible with `./validate_e2e.sh`.

**This is a partial, honest result, not a clean win.** An earlier ensemble
that included the MRTX1133-bound structure (7RPZ) reported 100% recovery —
but that structure's own bound inhibitor is what holds the pocket open, so
that result was circular and has been retracted (see
[docs/RESEARCH_CORRECTIONS.md](docs/RESEARCH_CORRECTIONS.md) #7). With 7RPZ
removed and the ensemble restricted to structures with no synthetic ligand
bound anywhere, the ranking-guided modes still clearly outperform every
ablated baseline (0.40 recall vs. 0.00–0.20), but do not fully recover the
site from static apo crystal heterogeneity alone. Closing that gap is
exactly the job real BioEmu generative sampling is intended to do — see
[docs/LIMITATIONS.md](docs/LIMITATIONS.md).

> Design note: the first version of this ranking used druggability + volume
> only and also selected the nucleotide site. Crypticity had to be *in the
> objective* for the system to work. See [docs/SCIENTIFIC_METHOD.md](docs/SCIENTIFIC_METHOD.md).

## Generalization: does it only work on KRAS?

The frozen KRAS pocket-ranking config (identical weights, thresholds, and
ligand library, no per-target tuning — enforced by
[tests/test_frozen_config.py](tests/test_frozen_config.py)) was pointed at
two proteins never used during development: **ABL1 kinase** (DFG-out/
imatinib pocket) and **PRMT5:MEP50** (EE-loop allosteric pocket). On ABL
kinase it produced a **documented negative result** — the ensemble was too
thin (2 apo structures) for the novelty signal to be meaningful, and the
policy correctly refused to dock into the resulting low-quality pocket
rather than fabricate a positive outcome. Full writeup, including exactly
why it failed and what that reveals about the method's requirements, in
[docs/GENERALIZATION.md](docs/GENERALIZATION.md).

## Quick start

```bash
conda env create -f environment.yml
conda activate conform
export PYTHONPATH=.

python scripts/run_experiment.py validate
./run_demo.sh
```

Individual commands:

```bash
python scripts/run_experiment.py run --target kras-g12d --closed-loop
python scripts/run_experiment.py run --target abl-kinase --closed-loop   # generalization test
python scripts/run_experiment.py run --target prmt5 --closed-loop        # generalization test
python scripts/run_experiment.py benchmark --target kras-g12d
python scripts/run_experiment.py ablate --mode static
python scripts/run_experiment.py report
python scripts/run_experiment.py dashboard
python scripts/run_experiment.py history --verbose
```

## Validation

```bash
./validate_e2e.sh      # 16 checks, executes every real tool
python -m pytest tests/ -q
```

`validate_e2e.sh` prints PASS only when the check actually succeeded — it
runs real fpocket, real Vina docking, a real 2-iteration closed loop, and
real report generation.

## Requirements

**No GPU required.** Everything currently implemented runs on CPU. A full
KRAS closed-loop run takes ~4 minutes on a laptop. BioEmu / OpenFold3 /
GNINA / DiffDock-Pocket / REINVENT 4 would each need CUDA; the provider
interfaces for them exist and are documented in
[docs/DEPENDENCIES.md](docs/DEPENDENCIES.md).

## Project structure

```
agent/          state.py · policies.py · loop_controller.py · discovery_score.py
tools/          structure · bioemu · structural_analysis · mdpocket · rdkit · docking · reinvent
pipelines/      engines.py (deterministic engine layer) · pocket_discovery.py (P0 runner)
db/             schema.sql · repository.py   (SQLite default, Postgres/pgvector declared)
evaluation/     metrics · baselines · ablations · cryptic_recovery
visualization/  dashboard.py (Streamlit) · molecular_viewer.py (py3Dmol) · plots.py (SVG)
scripts/        run_experiment.py (CLI) · generate_report.py
configs/        kras_g12d.yaml
tests/          34 tests
docs/           ARCHITECTURE · SCIENTIFIC_METHOD · EXPERIMENT_PROTOCOL · BENCHMARKS
                LIMITATIONS · RESEARCH_CORRECTIONS · DEPENDENCIES · IMPLEMENTATION_STATUS
```

## Documentation

| Document | Contents |
|---|---|
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | Layering, closed loop, why the agent can't fake results |
| [SCIENTIFIC_METHOD.md](docs/SCIENTIFIC_METHOD.md) | Hypothesis, Discovery Score, blind-evaluation discipline |
| [EXPERIMENT_PROTOCOL.md](docs/EXPERIMENT_PROTOCOL.md) | Exact parameters, reproduction steps, artifacts |
| [BENCHMARKS.md](docs/BENCHMARKS.md) | Baselines, ablations, what each metric does and does not mean |
| [LIMITATIONS.md](docs/LIMITATIONS.md) | What did not run and the scientific consequences |
| [RESEARCH_CORRECTIONS.md](docs/RESEARCH_CORRECTIONS.md) | Where the original plan was wrong or infeasible |
| [DEPENDENCIES.md](docs/DEPENDENCIES.md) | Versions, licenses, GPU-worker isolation strategy |
| [IMPLEMENTATION_STATUS.md](docs/IMPLEMENTATION_STATUS.md) | Component-by-component status table |

## Licenses

All dependencies are open source. Note that Open Babel and MDAnalysis are
GPL-licensed and are used here as a library/subprocess — see
[docs/DEPENDENCIES.md](docs/DEPENDENCIES.md) for the full audit.

## Citation

If this architecture is useful, cite the underlying tools — RDKit, fpocket
(Le Guilloux et al. 2009), AutoDock Vina (Eberhardt et al. 2021), MDAnalysis,
Open Babel — and the KRAS G12D / MRTX1133 structural literature (PDB 7RPZ).
