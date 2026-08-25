# Benchmarks

Run with `python evaluation/ablations.py`. Results are written to
`artifacts/ablation_report.json`. **The table below is regenerated from real
executed runs — see that file for the authoritative current numbers.**

## Results (real executed run, 2026-08-25, CPU-only, apo-only ensemble)

**⚠ Superseded numbers:** an earlier run (2026-08-12) reported 1.00 recall
for `no-ligand-optimization`/`conform-agent`. That ensemble included PDB
7RPZ, which is co-crystallized with MRTX1133 — the inhibitor that itself
holds the Switch-II pocket open. That result was circular (see
`docs/RESEARCH_CORRECTIONS.md` #7) and is retracted. The numbers below are
from the corrected ensemble containing **only structures with no synthetic
ligand bound anywhere** (5US4, 5XCO, 7F0W, 7EYX — verified apo-like by
HETATM audit).

```
mode                     states  pockets  recall  best_kcal  discovery  iters     sec
-------------------------------------------------------------------------------------
static                        1       19    0.00      -9.85      0.384      1      58
random                        4       47    0.20      -6.94      0.605      1      66
no-pocket-guidance            4       47    0.00      -9.85      0.482      1      63
no-ligand-optimization        4       47    0.40      -8.31      0.763      1      72
conform-agent                 4       47    0.40      -9.07      0.785      2     335
```

### What this supports
The ranking-guided modes (`no-ligand-optimization`, `conform-agent`) still
outperform every ablated mode on cryptic-residue recall (0.40 vs 0.00–0.20),
and land on the same pocket (5XCO:pocket2, 891 Å³, novelty 1.00 — fully
absent from the apo baseline). `random` gets partial credit (0.20) by luck,
landing on 7F0W's Switch-I-open pocket. `static` and `no-pocket-guidance`
both land on 0.00 by defaulting to a pocket at the apo baseline itself
(no novelty).

**The pocket is only partially recovered: H95 and Y96, not Q99, V9, or D69.**
Recall 0.40 is below the `recovered_cryptic_site` threshold (0.6) — by this
project's own criterion, **the corrected, honest ensemble does not fully
recover the Switch-II site from apo crystal heterogeneity alone.** This is
the expected and scientifically believable outcome: a handful of static apo
crystal forms happen to show partial pre-organization of the pocket lip
(consistent with a conformational-selection mechanism), but full opening to
the ligand-bound geometry requires either (a) a ligand physically forcing it
open — which is what the retracted run was secretly measuring — or (b)
broader conformational sampling than 4 discrete crystal forms can provide.
This is exactly the gap real BioEmu sampling is intended to close; see
`docs/LIMITATIONS.md` for status.

**Raw docking score is still a misleading metric.** `static` posts the best
affinity in the table (−9.85 kcal/mol) while scoring 0.00 recall — it docks
into the large, always-open nucleotide pocket. Ranking methods by affinity
alone would pick the worst-performing mode first, in both the old and the
corrected results.

### Ligand optimization: small real improvement, still inconclusive
`conform-agent` improved best affinity by **0.76 kcal/mol** over the library
baseline (−8.31 → −9.07), a larger and more genuine gain than the 0.07
kcal/mol noise-level change seen in the earlier (contaminated) run — but on
a different pocket and different seed ligand (Imatinib_fragment rather than
the piperazine scaffold), so the two results are not directly comparable.
`conform-agent`'s Discovery Score (0.785) now exceeds
`no-ligand-optimization` (0.763), reversing the earlier (also
noise-level) direction. With a single run, single target, and 12 analogs,
neither result should be treated as a validated claim that optimization
helps — it should be re-tested with real REINVENT 4 on a GPU, and ideally
across more than one seed ligand, before any claim is made either way.

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
