# Architecture

## Layering principle
Three strictly separated layers. Scientific engines work with no LLM and no
agent; the agent only decides *what to run* and *how to interpret it*; the
presentation layer only reads artifacts.

```
                     ┌─────────────────────────────┐
                     │   Agent orchestration       │
                     │   policies.py (state machine)│
                     │   loop_controller.py         │
                     └──────────────┬──────────────┘
                                    │ Action + params
                     ┌──────────────▼──────────────┐
                     │   pipelines/engines.py      │   ← deterministic, testable
                     └──────────────┬──────────────┘
          ┌─────────────────┬───────┴───────┬──────────────────┐
          ▼                 ▼               ▼                  ▼
   structure_tool     bioemu_tool    mdpocket_tool      docking_tool
   (RCSB / OF3 /      (BioEmu /      (mdpocket /        (GNINA /
    ESMFold)           experimental)  fpocket / P2Rank)  Vina)
          │                 │               │                  │
          └─────────────────┴───────┬───────┴──────────────────┘
                                    ▼
                     ┌─────────────────────────────┐
                     │ agent/discovery_score.py    │  ← ONLY place a score
                     │ (pure deterministic Python) │     is ever computed
                     └──────────────┬──────────────┘
                                    ▼
                     ┌─────────────────────────────┐
                     │ db/repository.py            │  SQLite (default)
                     │ scientific memory           │  PostgreSQL/pgvector (declared)
                     └──────────────┬──────────────┘
                                    ▼
              artifacts/<experiment_id>/{structures,pockets,ligands,
                                         docking,metrics,report}
                                    ▼
              visualization/dashboard.py · scripts/generate_report.py
```

## The closed loop
```
observe state → policy.decide() → execute engine → validate artifacts
     ▲                                                      │
     └──────── update state + persist to memory ◄───────────┘
```

Terminates on: confidence threshold met, budget exhausted, MAX_ITERATIONS,
duplicate action detected, or no useful action remaining.

## Iteration accounting
`ExperimentState.iteration` counts **closed-loop experiment cycles**
(SCREEN_LIGANDS / OPTIMIZE_LIGAND / VALIDATE_POSE). Setup stages
(GENERATE_ENSEMBLE, ANALYZE_ENSEMBLE, FIND_POCKETS, SELECT_POCKET) are
one-time prerequisites tracked by `step_count` and do not consume the
MAX_ITERATIONS budget. A separate hard cap (`MAX_TOTAL_STEPS = 20`) bounds
total actions so a policy bug cannot loop forever.

## Why the agent cannot fake results
- Discovery Score is computed only in `agent/discovery_score.py`, from
  numbers produced by executed tools.
- Every action's output artifacts are existence/size-validated before the
  state is updated.
- Failures are written to scientific memory and terminate the branch; there
  is no "recover by assuming success" path.
- Ground-truth residue overlap is computed but explicitly **excluded** from
  pocket ranking, so the system cannot score well by peeking at the answer.

## Fallback hierarchy
| Stage | Preferred | Fallback chain actually used here |
|---|---|---|
| Structure | OpenFold3 | ESMFold → **RCSB deposited structure** |
| Ensemble | BioEmu | **experimental multi-structure ensemble** |
| Pockets | mdpocket | **fpocket per structure** → P2Rank |
| Docking | GNINA | **AutoDock Vina** |
| Ligand opt. | REINVENT4 | **RDKit R-group enumeration** → library only |

Bold = what ran in this build. The mode actually used is recorded in every
manifest, report, and dashboard view.
