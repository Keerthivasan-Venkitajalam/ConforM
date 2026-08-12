# Benchmarks

Run with `python evaluation/ablations.py`. Results are written to
`artifacts/ablation_report.json`. **The table below is regenerated from real
executed runs — see that file for the authoritative current numbers.**

## Results (real executed run, 2026-08-12, CPU-only)

```
mode                     states  pockets  recall  best_kcal  discovery  iters     sec
-------------------------------------------------------------------------------------
static                        1        9    0.00     -10.65      0.439      1      64
random                        4       49    0.00      -8.55      0.659      1      62
no-pocket-guidance            4       49    0.00     -10.65      0.551      1      64
no-ligand-optimization        4       49    1.00      -9.60      0.777      1      68
conform-agent                 4       49    1.00      -9.67      0.776      2     237
```

### What this supports
**Conformational sampling + crypticity-aware ranking is what finds the
cryptic site.** Only the two modes with both (`no-ligand-optimization` and
`conform-agent`) recover the Switch-II ground truth (recall 1.00). All three
degraded modes score 0.00:
- `static` has no ensemble, so the cryptic state is never sampled.
- `no-pocket-guidance` has the ensemble but no ranking, and defaults to the
  always-open nucleotide site.
- `random` has the ensemble but no ranking logic, and lands on an
  undruggable cavity (druggability 0.00).

**Raw docking score is a misleading metric, exactly as warned below.**
`static` posts the *best* affinity in the whole table (−10.65 kcal/mol) while
completely failing the actual task. It achieves this by docking into the
large, always-open nucleotide pocket. Anyone comparing these systems on
affinity alone would rank the worst method first.

### Negative result: ligand optimization did not help
`conform-agent` (0.776) scores **marginally lower** than
`no-ligand-optimization` (0.777), and the optimization round improved
affinity by only **0.07 kcal/mol** (−9.60 → −9.67) — far inside docking noise
— at ~3.5× the runtime (237 s vs 68 s).

This does **not** support the research plan's Ablation 3 expectation that
removing iterative refinement would "demonstrate a significantly lower final
binding affinity." On this target, with this 10-molecule library and RDKit
R-group enumeration standing in for REINVENT 4, the optimization loop bought
essentially nothing. Plausible reasons: the fallback optimizer only decorates
aromatic C–H positions rather than performing reward-driven scaffold
optimization; the seed scaffold may already be near the ceiling of what Vina
scores in this pocket; and a single round on 12 analogs is a very small
search. This should be re-tested with real REINVENT 4 on a GPU before any
claim about closed-loop optimization is made.

## Modes compared

| Mode | What it ablates |
|---|---|
| `static` | No conformational sampling: apo baseline structure only (research plan Baseline 1) |
| `random` | Ensemble generated, but pocket chosen at random — no ranking logic (Baseline 2) |
| `no-pocket-guidance` | Ensemble generated, but no volumetric ranking: first cavity of first state |
| `no-ligand-optimization` | Full agent with the optimization iteration disabled |
| `conform-agent` | Full closed-loop system |

## Baseline deliberately NOT run
**Baseline 3 (1 µs classical MD).** The research plan's own estimate is
>10,000 GPU-hours to sample Switch-II opening by unbiased MD. It is not run
and is **not approximated** — an invented MD number would be worse than a
missing one. Its absence is recorded in `ablation_report.json` under
`omitted_baseline`.

## Metrics reported
- `cryptic_residue_recall` — fraction of documented Switch-II ground-truth
  residues present in the selected pocket's lining residues.
  **This is a residue-level proxy, not the Discretized Volume Overlap (DVO)**
  named in the research plan; true DVO requires voxelizing the predicted
  cavity grid against the experimental ligand-occupied volume, which is not
  implemented here.
- `best_affinity_kcal` — best AutoDock Vina empirical score achieved.
- `best_discovery_score` — deterministic Discovery Score (absolute affinity
  normalization, so values are comparable across modes).
- `closed_loop_iterations`, `runtime_seconds`, `gpu_hours` (0.0 — CPU-only).
- `enrichment` — smoke-test only; the library is 10 molecules and this is not
  a statistically powered enrichment factor.

## Interpretation guidance
The meaningful comparison is **`cryptic_residue_recall`**, not docking score.
A mode that docks into the large always-open nucleotide site can post a
respectable affinity while completely failing the actual task (finding the
cryptic site). Comparing the modes on affinity alone would be misleading.

Do not claim superiority over any baseline beyond what
`artifacts/ablation_report.json` actually shows, and note that with a single
target and a 10-molecule library these results demonstrate that the
*mechanism works*, not that it generalizes.

## CryptoBench
Not integrated. Evaluating the ~1,107 apo-holo pairs properly is far beyond
the compute available here, and a partial run would produce a number that
looks like a benchmark result without being one. `evaluation/cryptic_recovery.py`
evaluates against the explicitly documented KRAS G12D Switch-II ground truth
instead. Full CryptoBench integration with `--subset` / `--target` scoping
remains future work.
