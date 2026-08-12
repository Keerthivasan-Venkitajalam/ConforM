# Implementation Status

Environment: macOS, **no CUDA GPU** (`nvidia-smi` absent). Internet available.
All numbers below come from actually executed runs; see
`artifacts/<experiment_id>/metrics/experiment_manifest.json` and
`artifacts/ablation_report.json`.

Verify everything with `./validate_e2e.sh` (16 checks, all currently PASS).

| Component | Status | Implementation | Missing work | Test status | Fallback used |
|---|---|---|---|---|---|
| Structure acquisition | **Done** (fallback tier) | `tools/structure_tool.py` — provider chain, fetches real RCSB structures | OpenFold3/ESMFold inference | E2E PASS | Yes: RCSB instead of OpenFold3/ESMFold |
| Conformational ensemble | **Done** (fallback tier) | `tools/bioemu_tool.py` — 4 real KRAS G12D structures (4DST, 5US4, 7RPZ, 5XCO) | Real BioEmu diffusion sampling | E2E PASS | Yes: experimental ensemble; NOT equilibrium samples |
| Structural analysis | **Done** (real) | `tools/structural_analysis.py` — Kabsch, RMSD/RMSF, PCA | TICA (deliberately excluded, invalid here) | E2E PASS | No |
| Pocket detection | **Done** (real, fallback tier) | `tools/mdpocket_tool.py` — real fpocket 4.0 per structure | True mdpocket trajectory mode | E2E PASS | Yes: per-structure fpocket |
| Cross-state pocket families | **Done** (real) | `pipelines/engines.py` — Jaccard clustering, persistence, novelty vs. apo baseline | Volumetric grid persistence | 4 unit tests PASS | No |
| Pocket ranking | **Done** (real) | Blind: druggability + novelty + volume; ground truth excluded | Weights not yet split to own YAML | Unit test PASS | No |
| Ligand validation | **Done** (real) | `tools/rdkit_tool.py` — sanitize, Lipinski, QED, ETKDG+MMFF | ChEMBL-scale library | 4 unit tests PASS | No |
| Docking | **Done** (real, fallback tier) | `tools/docking_tool.py` — real AutoDock Vina 1.2.7 | GNINA CNN rescoring, DiffDock-Pocket | E2E PASS | Yes: Vina empirical scoring |
| Ligand optimization | **Done** (fallback tier) | `tools/reinvent_tool.py` — RDKit R-group enumeration | REINVENT 4 RL loop | E2E PASS | Yes; **and it did not improve results — see BENCHMARKS.md** |
| Discovery Score | **Done** (real) | `agent/discovery_score.py` — deterministic, absolute affinity normalization | — | 6 unit tests PASS | No |
| Agent state machine | **Done** (real) | `agent/state.py`, `agent/policies.py` | LLM proposal layer on top of policy | 10 unit tests PASS | No |
| Closed loop | **Done** (real) | `agent/loop_controller.py` — 2 real iterations executed | Multi-pocket exploration branches | E2E PASS (asserts ≥2 iterations) | No |
| Scientific memory | **Done** (real) | `db/repository.py` + `schema.sql`, SQLite | PostgreSQL/pgvector code path | 5 unit tests PASS | SQLite instead of Postgres |
| Duplicate detection | **Done** (real) | SHA-256 action hashes, enforced in loop | — | Unit + E2E PASS | No |
| Evaluation / baselines / ablations | **Done** (real) | `evaluation/` — 5 modes all executed | CryptoBench full dataset | E2E PASS | Baseline 3 (1 µs MD) omitted, not approximated |
| Report generator | **Done** (real) | `scripts/generate_report.py` → HTML with SVG plots + py3Dmol | — | E2E PASS | No |
| Dashboard | **Done** (real) | `visualization/dashboard.py` — Streamlit, consumes real artifacts | — | E2E PASS (import + parse) | No |
| CLI | **Done** | `scripts/run_experiment.py` — setup/validate/run/benchmark/ablate/report/dashboard/history | — | Manually verified | No |
| Docker | **Partial** | `docker-compose.yml` + `docker/Dockerfile` written | **Never built or run** on this host | Not tested | — |
| MCP layer | **Not started** | — | Pydantic tool schemas, MCP server | — | Local pipeline works without it |
| CryptoBench | **Not started** | — | Dataset download + evaluation harness | — | KRAS-specific ground truth used |
| DiffDock-Pocket / OpenMM | **Not started** | — | P5 items | — | — |

## Headline result (real)
Blind ranking selects **7RPZ:pocket1** — novelty 1.00 (absent from apo
baseline), persistence 0.25 (present in 1 of 4 states), **100% overlap with
documented Switch-II ground-truth residues**. The `static` baseline scores
0.00 recall. Full table in [BENCHMARKS.md](BENCHMARKS.md).

## Honest caveats
- Ligand optimization produced a 0.07 kcal/mol change — inside noise. It did
  **not** reproduce the research plan's expected ablation effect.
- The ensemble is 4 crystal structures, not a Boltzmann sample; state
  populations are uniform placeholders.
- The library is 10 molecules; enrichment is a smoke test, not a benchmark.
- Docking is Vina empirical scoring, unvalidated against CNN or experimental data.
- `docker-compose.yml` and the Dockerfile are written but were never built.
