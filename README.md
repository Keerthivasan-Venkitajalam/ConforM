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
- The core pipeline is CPU-only and reproducible without a GPU. Two of the
  originally GPU-only models, **BioEmu and GNINA, have since been executed
  for real** on a CUDA machine (2026-08-29) — see the GPU-verified result
  below and [docs/RESEARCH_CORRECTIONS.md](docs/RESEARCH_CORRECTIONS.md) #8,
  which also documents two bugs that first produced a false GPU result and
  were caught before being reported anywhere. **OpenFold3, DiffDock-Pocket,
  and REINVENT 4 still have not been run**; documented fallbacks are used
  for those and labeled as such in every manifest, report, and dashboard
  view. See [docs/LIMITATIONS.md](docs/LIMITATIONS.md).

---

## What it does

```
KRAS G12D
   ↓  structure acquisition        (OpenFold3 → ESMFold → RCSB)      [RCSB used]
   ↓  conformational ensemble      (BioEmu → experimental ensemble)  [BioEmu verified on GPU; experimental fallback on CPU]
   ↓  structural analysis          RMSD / RMSF / PCA
   ↓  pocket detection             (mdpocket → fpocket → P2Rank)     [fpocket used]
   ↓  cross-state pocket families  persistence + novelty vs. apo baseline
   ↓  blind pocket ranking         druggability + novelty + volume
   ↓  ligand validation            RDKit sanitize / 3D embed / Lipinski / QED
   ↓  docking                      (Vina, GNINA CNN rescore on top)  [Vina always; GNINA verified on GPU]
   ↓  Discovery Score              deterministic Python, never LLM-generated
   ↓  agent decision               screen? optimize? stop?
   ↺  next experiment
```

The agent decides what to run next from an explicit state machine
([agent/policies.py](agent/policies.py)) — it is not a hard-coded linear
pipeline. Scientific engines are independently callable without any agent or
LLM.

## Key result (real GPU-verified run, real BioEmu ensemble, 2026-08-29)

On a real CUDA GPU, `tools/bioemu_tool.py` ran genuine BioEmu diffusion
inference from the apo KRAS G12D sequence alone — no experimental structure
as input. Blind ranking — with ground-truth overlap explicitly **excluded**
from the objective — selects a cavity that is:

- **absent from the apo baseline structure** (novelty 1.00)
- **present in only a handful of a 94–100-state real equilibrium sample**
  (persistence 0.04–0.16 — genuinely transient, not a coin flip the way a
  4-structure ensemble risks being)
- **covering 60% of the documented Switch-II ground-truth residues**
  (replicated independently across two separate BioEmu runs)

Full corrected 5-mode ablation table (static / random / no-pocket-guidance /
no-ligand-optimization / conform-agent, all real BioEmu runs) in
[docs/RESEARCH_CORRECTIONS.md](docs/RESEARCH_CORRECTIONS.md) #8. GNINA CNN
rescoring was also run for real on this verified result (top hit CNNscore
0.683).

**Update (2026-08-30): the #8 cross-mode sampling-variance caveat is now
fixed.** Every non-`static` mode was re-run against the *same* shared,
already-verified 99-state ensemble (no new GPU sampling needed) — see
[docs/RESEARCH_CORRECTIONS.md](docs/RESEARCH_CORRECTIONS.md) #9. With the
confound removed, `no-ligand-optimization` and `conform-agent` converge on
the identical selected pocket (0.60 recall both), and a 10,000-permutation
test on the real, unmodified ranking of all 149 detected pocket families
found the algorithm places the best-overlap family available at rank 2 of
149 — a rank that good occurs in only 1.44% of random relabelings
(**p = 0.0144**, not pre-specified).

**This supersedes an earlier false GPU result** (`n_states=1` for every
mode, a `conform-agent` recall of 0.80) that was caught and retracted before
being used anywhere — a trajectory-extraction bug silently discarded ~99 of
100 real BioEmu samples down to 1 reference frame. Full root-cause writeup
in RESEARCH_CORRECTIONS.md #8.

## Earlier result (CPU-only, no GPU required, apo-only crystal ensemble)

Kept as a separately valid, no-GPU-required reproduction path — this is what
`./validate_e2e.sh` exercises. Blind ranking on 4 real KRAS G12D crystal
structures selects a cavity **partially covering the documented Switch-II
ground-truth residues** (H95, Y96 recovered; Q99, V9, D69 not — recall
0.40), present in only 1 of 4 sampled states (persistence 0.25). The
`static` baseline, given only the apo structure, instead selects the
always-open nucleotide site: ground-truth recall **0.00**.

**This is a partial, honest result, not a clean win.** An earlier ensemble
that included the MRTX1133-bound structure (7RPZ) reported 100% recovery —
but that structure's own bound inhibitor is what holds the pocket open, so
that result was circular and has been retracted (see
[docs/RESEARCH_CORRECTIONS.md](docs/RESEARCH_CORRECTIONS.md) #7). With 7RPZ
removed and the ensemble restricted to structures with no synthetic ligand
bound anywhere, the ranking-guided modes still clearly outperform every
ablated baseline (0.40 recall vs. 0.00–0.20), but do not fully recover the
site from static apo crystal heterogeneity alone. Closing that gap is
exactly what the GPU-verified real BioEmu result above does — see
[docs/LIMITATIONS.md](docs/LIMITATIONS.md).

> Design note: the first version of this ranking used druggability + volume
> only and also selected the nucleotide site. Crypticity had to be *in the
> objective* for the system to work. See [docs/SCIENTIFIC_METHOD.md](docs/SCIENTIFIC_METHOD.md).

## Generalization: does it only work on KRAS?

The frozen KRAS pocket-ranking config (identical weights, thresholds, and
ligand library, no per-target tuning — enforced by
[tests/test_frozen_config.py](tests/test_frozen_config.py)) was pointed at
two proteins never used during development: **ABL1 kinase** (DFG-out/
imatinib pocket) and **PRMT5:MEP50** (EE-loop allosteric pocket).

| Target | Recall | Outcome |
|---|---|---|
| KRAS G12D | 0.40 | Partial recovery; beats all ablated baselines |
| ABL1 kinase | 0.00 | **Documented negative result** — 2-structure ensemble too thin for the novelty signal to be meaningful; agent correctly refused to dock rather than fabricate a result |
| PRMT5:MEP50 | 0.25 | Real partial recovery on a structurally unrelated fold (druggability 0.998, best Discovery Score 0.735) |

Consistent partial signal, an honest failure mode instead of a fabricated
success, and zero per-target tuning across three unrelated protein folds —
full writeup, root-cause analysis, and per-target manifests in
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

**No GPU required for the core pipeline.** Everything in `./validate_e2e.sh`
runs on CPU; a full KRAS closed-loop run takes ~4 minutes on a laptop with
the CPU-only experimental-ensemble fallback (`bioemu.enabled: false`, or no
CUDA GPU detected).

**A GPU unlocks the real BioEmu and GNINA path** (`bioemu.enabled: true` in
`configs/kras_g12d.yaml`, plus `scripts/gnina_rescore.py`) — verified for
real on an RTX 4060 via WSL2, see the GPU-verified result above and
[scripts/gpu_session.sh](scripts/gpu_session.sh) for a scripted end-to-end
GPU session. A single closed-loop run with real BioEmu (100 samples) plus
GNINA rescoring takes ~13 minutes. OpenFold3, DiffDock-Pocket, and
REINVENT 4 would also each need CUDA and have not been run; the provider
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
tests/          49 tests
docs/           ARCHITECTURE · SCIENTIFIC_METHOD · EXPERIMENT_PROTOCOL · BENCHMARKS
                LIMITATIONS · RESEARCH_CORRECTIONS · DEPENDENCIES · IMPLEMENTATION_STATUS
                GENERALIZATION
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
| [GENERALIZATION.md](docs/GENERALIZATION.md) | Zero-shot evaluation on ABL1 kinase and PRMT5:MEP50; real-BioEmu KRAS addendum |

## Licenses

All dependencies are open source. Note that Open Babel and MDAnalysis are
GPL-licensed and are used here as a library/subprocess — see
[docs/DEPENDENCIES.md](docs/DEPENDENCIES.md) for the full audit.

## Citation

If this architecture is useful, cite the underlying tools — RDKit, fpocket
(Le Guilloux et al. 2009), AutoDock Vina (Eberhardt et al. 2021), MDAnalysis,
Open Babel — and the KRAS G12D / MRTX1133 structural literature (PDB 7RPZ).
