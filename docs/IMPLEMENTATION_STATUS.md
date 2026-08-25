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

## Headline result (real, apo-only ensemble, 2026-08-25)
Blind ranking selects **5XCO:pocket2** — novelty 1.00 (absent from apo
baseline), persistence 0.25 (present in 1 of 4 states), **40% recall of
documented Switch-II ground-truth residues** (H95, Y96 recovered; Q99, V9,
D69 not). The `static` baseline scores 0.00 recall. Full table in
[BENCHMARKS.md](BENCHMARKS.md).

**This supersedes an earlier 100%-recall result** obtained when PDB 7RPZ
(MRTX1133-bound) was included in the discovery ensemble — that inhibitor
itself forces the pocket open, making that result circular. It has been
retracted; see [RESEARCH_CORRECTIONS.md](RESEARCH_CORRECTIONS.md) #7 for the
full HETATM audit and fix. The corrected ensemble contains only structures
with no synthetic ligand bound anywhere.

## Honest caveats
- **The corrected ensemble only partially recovers the cryptic site** (40%,
  below this project's own 60% "recovered" threshold). Static apo crystal
  heterogeneity alone is not sufficient to fully open Switch-II — real
  generative sampling (BioEmu) is the intended way to close this gap and has
  not yet been run (no GPU on this machine; `tools/bioemu_tool.py` now has a
  real, CUDA-gated implementation ready to run on rented compute).
- Ligand optimization in the corrected run improved affinity by 0.76
  kcal/mol (a different pocket/seed ligand than the earlier retracted run,
  so not directly comparable) — still a single run on 12 analogs, not a
  validated claim.
- The ensemble is 4 crystal structures, not a Boltzmann sample; state
  populations are uniform placeholders.
- The library is 10 molecules; enrichment is a smoke test, not a benchmark.
- Docking is Vina empirical scoring, unvalidated against CNN or experimental data.
- `docker-compose.yml` and the Dockerfile are written but were never built.
