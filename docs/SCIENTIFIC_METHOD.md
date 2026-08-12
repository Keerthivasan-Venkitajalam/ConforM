# Scientific Method

## Hypothesis
A system that samples multiple protein conformational states and ranks
cavities by *novelty relative to the apo baseline* will identify cryptic
binding sites that a single static structure does not expose.

## Why crypticity must be in the objective
The first working version of this pipeline ranked pockets by druggability
and volume alone. On KRAS G12D it selected the **nucleotide (GDP) site** —
the largest, most druggable cavity in the protein, and one that is always
open. Ground-truth overlap with the Switch-II cryptic site was 0.00.

That is the correct behavior for that objective, and the wrong objective for
this problem. A cavity that is open in the apo structure is not a cryptic
pocket discovery. The ranking function was changed to include:

```
novelty = (max_volume_across_states - volume_in_apo_baseline) / max_volume
```

with cross-state pocket correspondence established by Jaccard similarity of
lining-residue sets (`pipelines/engines.py:cluster_pocket_families`). With
novelty weighted at 0.40, blind ranking selects the Switch-II pocket.

This correction is recorded because it is the substantive scientific content
of the system: the result depends entirely on the objective function, not on
the tooling.

## Discovery Score
```
DiscoveryScore = w_pocket  · pocket_novelty          (volume absent from apo baseline)
               + w_volume  · normalized_volume
               + w_binding · normalized_binding      (absolute Vina reference scale)
               + w_state   · state_novelty           (1 − persistence across ensemble)
               + w_ligand  · ligand_quality          (QED, Lipinski-penalized)
               − w_invalid · structural_penalty
               − w_cost    · computational_cost
```
Weights in `configs/kras_g12d.yaml`. Implementation and range documented in
`agent/discovery_score.py`.

**Binding normalization is absolute, not run-relative.** Anchors are
−12.0 kcal/mol → 1.0 and −4.0 kcal/mol → 0.0. Run-relative min/max
normalization was rejected because it awards the best ligand of *every* run a
score near 1.0 regardless of actual affinity, which makes Discovery Scores
incomparable across baselines and invites reward hacking.

## Blind evaluation discipline
`ground_truth_overlap` is computed for every pocket and carried through the
pipeline, but it is **never** an input to ranking or scoring. It is used only
post-hoc, in `evaluation/`, to measure recovery. Any change that lets
ground truth influence selection invalidates the benchmark.

## Known ground-truth recovery vs. novel hypothesis
The competition claim must distinguish these:

- **Ground-truth recovery (validation).** The system blindly re-identifies
  the documented Switch-II pocket in KRAS G12D. This validates the method.
  It is not a discovery — the pocket has been known and drugged (MRTX1133).
- **Novel hypothesis generation (future potential).** Ligand rankings
  produced against that pocket are computational hypotheses requiring
  experimental validation. No binding, potency, or efficacy claim follows
  from a docking score.

## Statistical honesty
- The fallback ensemble is 4 crystal structures, not an equilibrium sample.
  State "populations" are uniform placeholders; no free-energy or
  thermodynamic claim can be made from them.
- The ligand library is 10 molecules. Enrichment computed on it is a smoke
  test, not a statistically powered enrichment factor, and is reported as
  such in `evaluation/metrics.py`.
- Docking scores are AutoDock Vina empirical-function values, unvalidated
  against GNINA CNN scores or experimental affinity data in this build.
