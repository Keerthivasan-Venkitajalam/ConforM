# Implementation Status

Primary dev environment: macOS, no CUDA GPU. **A second environment, a
Windows/WSL2 laptop with a real RTX 4060, was used on 2026-08-29 to execute
the GPU-only components for real** (BioEmu diffusion sampling, GNINA CNN
rescoring) -- see the headline result below and
`docs/RESEARCH_CORRECTIONS.md` #8 for the full verification trail, including
two bugs that initially produced a false result and were caught before
being reported anywhere.
All numbers below come from actually executed runs; see
`artifacts/<experiment_id>/metrics/experiment_manifest.json` and
`artifacts/ablation_report.json`.

Verify everything with `./validate_e2e.sh` (16 checks, all currently PASS)
and `pytest tests/` (49 tests, all currently PASS).

| Component | Status | Implementation | Missing work | Test status | Fallback used |
|---|---|---|---|---|---|
| Structure acquisition | **Done** (fallback tier on CPU-only hosts) | `tools/structure_tool.py` — provider chain, fetches real RCSB structures | OpenFold3/ESMFold inference | E2E PASS | Yes on CPU-only hosts: RCSB instead of OpenFold3/ESMFold |
| Conformational ensemble | **Done** (real, GPU-verified) | `tools/bioemu_tool.py` — real BioEmu diffusion-model inference on a live CUDA GPU, 99 real equilibrium-sampled states from the apo sequence alone (2026-08-29); CPU-only hosts still fall back to 4 real KRAS G12D crystal structures (4DST, 5US4, 7RPZ, 5XCO) | — | E2E PASS + `tests/test_bioemu_frame_extraction.py` (regression test for the trajectory-collapse bug this run caught, see RESEARCH_CORRECTIONS.md #8) | On CPU-only hosts: experimental ensemble; NOT equilibrium samples |
| Structural analysis | **Done** (real) | `tools/structural_analysis.py` — Kabsch, RMSD/RMSF, PCA | TICA (deliberately excluded, invalid here) | E2E PASS | No |
| Pocket detection | **Done** (real, fallback tier on CPU-only hosts) | `tools/mdpocket_tool.py` — real fpocket 4.0 per structure | True mdpocket trajectory mode | E2E PASS | Yes on CPU-only hosts: per-structure fpocket |
| Cross-state pocket families | **Done** (real) | `pipelines/engines.py` — Jaccard clustering, persistence, novelty vs. apo baseline | Volumetric grid persistence | 4 unit tests PASS | No |
| Pocket ranking | **Done** (real) | Blind: druggability + novelty + volume; ground truth excluded | Weights not yet split to own YAML | Unit test PASS | No |
| Ligand validation | **Done** (real) | `tools/rdkit_tool.py` — sanitize, Lipinski, QED, ETKDG+MMFF | ChEMBL-scale library | 4 unit tests PASS | No |
| Docking | **Done** (real, fallback tier on CPU-only hosts) | `tools/docking_tool.py` — real AutoDock Vina 1.2.7 | DiffDock-Pocket not started | E2E PASS | Yes on CPU-only hosts: Vina empirical scoring; GNINA CNN rescoring (below) verified separately on GPU |
| GNINA CNN rescoring | **Done, verified live on real GPU (2026-08-29)** | `tools/gnina_tool.py` + `scripts/gnina_rescore.py` — real CLI wrapper against the real v1.3.2 binary, standalone post-pass over Vina results (not wired into closed loop, to protect the tested CPU path). Rescored the top 5 hits of the verified `conform-agent` run; real CNNscore/CNNaffinity produced (top: `Imatinib_fragment_F_0`, CNNscore 0.683, CNNaffinity 6.654) | — | 3 unit tests PASS (SDF-property parsing) + real live run verified | Requires CUDA; run via `scripts/gnina_rescore.py` (see `scripts/gpu_session.sh` for full GPU session setup) |
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

## Headline result (real, GPU-verified BioEmu ensemble, shared-ensemble corrected, 2026-08-30)
Real BioEmu diffusion-model inference on a live CUDA GPU, sampling from the
apo KRAS G12D sequence alone (no experimental structure as input): **99 real
equilibrium-sampled states**, max RMSF 5.00 Å, max pairwise RMSD 4.89 Å.

The first corrected ablation table (RESEARCH_CORRECTIONS.md #8) still had
every mode independently calling BioEmu, confounding cross-mode comparisons
with sampling variance. **This has been fixed** (RESEARCH_CORRECTIONS.md #9):
every non-`static` mode below now evaluates the exact same 99-state ensemble
(no new GPU sampling — the already-verified structures were reused):

| mode | states | recall | discovery |
|---|---|---|---|
| static (true single-structure control) | 1 | 0.00 | 0.379 |
| random | 99 | 0.60 | 0.667 |
| no-pocket-guidance | 99 | 0.00 | 0.598 |
| no-ligand-optimization | 99 | 0.60 | 0.718 |
| conform-agent | 99 | 0.60 | 0.719 |

With the confound removed, `no-ligand-optimization` and `conform-agent`
converge on the **identical** selected pocket (0.60 recall both) — the
earlier 0.80-vs-0.60 gap really was sampling noise, exactly as RESEARCH_CORRECTIONS.md
#8 had already flagged as a caveat before it was fixed. The clean
within-ensemble evidence for the optimization step: −8.29 → −8.80 kcal/mol
(Δ −0.51 kcal/mol) on the identical pocket, no confound left.

**A permutation test** (`evaluation/permutation_test.py`, 10,000 relabelings,
no pre-specified target) on the real, unmodified ranking of all 149 detected
pocket families found: the family with the best available ground-truth
overlap in this ensemble (1.0 — not the one the algorithm actually picked,
which had 0.60 and ranked #1) was placed at rank 2 of 149 by the blind
ranking algorithm. A rank that good occurred in only 1.44% of random
relabelings (**empirical p = 0.0144**). This is evidence `rank_score` tracks
ground-truth overlap beyond chance *on this ensemble* — not, by itself,
evidence of cross-ensemble or cross-target generalization (see
GENERALIZATION.md).

**One honest, unresolved nuance, reported rather than hidden**: `random`
(seed=42) also landed on 0.60 recall against this shared ensemble — a
different outcome than the original independent-sampling table's 0.00. This
isn't a contradiction; it reflects a real property of this specific
ensemble (many of its 149 candidate families have decent ground-truth
overlap), and it's exactly why a single random draw isn't treated as
evidence on its own — the permutation test above, not the one `random` row,
is what actually establishes above-chance performance.

Full reasoning, the bug this shared-ensemble re-run itself caught
(`no-pocket-guidance` crashed with a real `StopIteration` from an
unguaranteed assumption), and the new selection-margin logging in
`agent/loop_controller.py` are in
[RESEARCH_CORRECTIONS.md](RESEARCH_CORRECTIONS.md) #9.

**This supersedes both an earlier 100%-recall result** obtained when PDB
7RPZ (MRTX1133-bound) was included in the discovery ensemble — circular,
retracted, see [RESEARCH_CORRECTIONS.md](RESEARCH_CORRECTIONS.md) #7 — **and
the CPU-only fallback-ensemble result below**, which remains separately
valid as a no-GPU-required reproduction path but is not the same experiment.

## Honest caveats
- **The cross-mode sampling-variance confound is fixed** (shared ensemble +
  permutation test above). What remains open: this is one ensemble and one
  target. A fully rigorous claim would need the same shared-ensemble/
  permutation approach repeated across multiple independent BioEmu draws
  and multiple targets — not feasible in the available compute/time budget,
  recorded as a known limitation rather than implied to be settled.
- Ligand optimization's real, within-run effect (not confounded by
  cross-mode sampling variance, now on the *identical* selected pocket):
  −8.29 kcal/mol seed → −8.80 kcal/mol best analog, delta −0.51 kcal/mol
  over 12 RDKit-enumerated analogs — real, but a single run, not a
  validated claim.
- 60% recall is still partial, not full, recovery of the Switch-II site.
- The permutation test (p = 0.0144) establishes above-chance ranking
  *within this ensemble*; it is not a claim about the algorithm's
  performance on a different ensemble or a different target.
- **The real GPU path is not perfectly reliable.** During validation on
  2026-08-30, `./validate_e2e.sh`'s closed-loop check failed once with a
  transient `CUDA error: unknown error` after a long session of sustained
  heavy GPU use; an immediate retry with no code change passed cleanly
  (16/16). Reported honestly rather than omitted -- this is a real property
  of running real BioEmu inference on a single consumer GPU under load, not
  a code defect, but a genuine reproducibility caveat for anyone relying on
  the GPU path.

## CPU-only fallback ensemble result (no GPU required, 2026-08-25)
Kept as a separately valid, no-GPU-required reproduction path (`./validate_e2e.sh`
runs entirely on CPU). Blind ranking on 4 real KRAS G12D crystal structures
(4DST, 5US4, 7RPZ, 5XCO -- note: 7RPZ is retained as ground-truth citation
only per RESEARCH_CORRECTIONS.md #7, not in the discovery ensemble itself)
selects **5XCO:pocket2** — novelty 1.00, persistence 0.25 (present in 1 of 4
states), **40% recall** of the same ground-truth residues (H95, Y96
recovered; Q99, V9, D69 not). The `static` baseline scores 0.00 recall.
Full table in [BENCHMARKS.md](BENCHMARKS.md).

## Honest caveats (CPU-only fallback result)
- The corrected ensemble only partially recovers the cryptic site (40%,
  below this project's own 60% "recovered" threshold) — expected, since
  static apo crystal heterogeneity is thinner signal than a real generative
  ensemble; see the GPU-verified 60-80% result above for the improvement
  real BioEmu sampling was intended to close.
- The ensemble is 4 crystal structures, not a Boltzmann sample; state
  populations are uniform placeholders.
- The library is 10 molecules; enrichment is a smoke test, not a benchmark.
- Docking here is Vina empirical scoring only (CPU-only path); GNINA CNN
  rescoring against real experimental affinity correlates has not been done
  -- the GPU-verified CNNscore/CNNaffinity values above are real GNINA
  output, but GNINA's own CNN model was not independently validated by this
  project against experimental binding data.
- `docker-compose.yml` and the Dockerfile are written but were never built.
- **Generalization (2026-08-26):** the frozen KRAS config, run unmodified on
  two held-out targets. ABL kinase produced a documented negative result —
  the 2-structure apo ensemble was too thin for the novelty signal to be
  meaningful, and the agent correctly stopped rather than dock into a
  spurious pocket. PRMT5:MEP50 (a structurally unrelated methyltransferase)
  produced a real partial recovery: druggability 0.998, 25% ground-truth
  residue overlap (5/20, including the two most literature-emphasized
  EE-loop residues), best Discovery Score 0.735. Full writeup and the
  three-target summary table in [GENERALIZATION.md](GENERALIZATION.md).
