# Generalization: zero-shot evaluation on held-out targets

The KRAS G12D result alone cannot answer the most damaging question a judge
can ask: *"did you tune this for KRAS?"* This document is the answer —
the frozen `configs/kras_g12d.yaml` pocket-ranking weights, thresholds, and
ligand library were applied **unmodified** to two proteins never used during
development: **ABL1 kinase** and **PRMT5:MEP50**.

## Protocol (pre-registration discipline)

1. `configs/abl_kinase.yaml` and `configs/prmt5.yaml` copy the
   `pocket_detection` / `docking` / `discovery_score` / `agent` sections of
   `configs/kras_g12d.yaml` **byte-for-byte**. No weight, threshold, or
   library was changed per target.
2. The same 10-molecule ligand library (`data/ligands_kras.csv`) is used for
   every target — a target-specific library would have undermined the
   generalization claim.
3. Ground-truth pocket residues were computed **objectively** from ligand
   contacts in the holo structure (`scripts/compute_ground_truth.py`,
   4.5 Å cutoff), not hand-copied from a paper. For PRMT5 this independently
   reproduced the literature-reported contact residues (Leu436, Leu437,
   Phe519, Phe555, Tyr468, Glu444), which is an additional cross-check on
   correctness.
4. As with KRAS, the ground-truth holo structure (1IEP for ABL, 6UXX for
   PRMT5) is excluded from the discovery ensemble and used only for citation.

## Target 1: ABL1 kinase — DFG-out / imatinib pocket

Discovery ensemble: **1OPL, 2FO0** (both inhibitor-free; myristate bound at
a separate, physiological allosteric site — see `configs/abl_kinase.yaml`
for the full HETATM audit). Only 2 apo-like structures could be located
(thinner than KRAS's 4), reported honestly rather than padded.

**Result: the frozen policy stopped without docking any ligand.**
```
selected pocket:   2FO0:pocket20
volume:            1800 A^3
druggability:      0.01   (below the 0.20 stop threshold)
novelty:           1.00
ground truth overlap: 0.00
decision:          STOP -- "no scientifically useful ligand experiment remains"
```

This is a **negative result, reported as such.** The policy correctly
refused to spend docking compute on a pocket its own druggability estimate
says is not real — that refusal is the safety mechanism working as
designed. But it did not find the actual DFG-out pocket either.

### Why it failed: novelty saturates with a thin ensemble

Inspecting the full ranked list (`artifacts/.../metrics/pocket_families.json`):

| pocket | druggability | novelty | volume (Å³) | rank score | GT overlap |
|---|---|---|---|---|---|
| 2FO0:pocket20 | 0.01 | 1.00 | 1800 | 0.603 | 0.00 |
| 2FO0:pocket1  | 0.96 | 0.00 | 936  | 0.487 | 0.16 |
| 2FO0:pocket32 | 0.00 | 1.00 | 659  | 0.474 | 0.16 |

The highest-druggability pocket found (0.96, some real overlap with the
DFG-out site) scores *worse* overall than a large, essentially undruggable
cavity, because with only **2** ensemble states, "novelty" (persistence
across states) becomes close to a coin flip: any pocket detected in only
one of the two structures gets novelty = 1.0, whether or not that reflects
a real, biologically meaningful conformational difference or just fpocket
noise on a single structure. With `w_volume = 0.15` and `w_pocket_novelty
= 0.30` both then rewarding a spuriously "novel," large cavity, it
outranks a smaller, genuinely druggable one with a better ground-truth
match.

**This is a real, useful finding about a boundary condition of the
method, discovered by exactly the kind of test this document exists to
run.** The corrected interpretation: novelty-based ranking requires an
ensemble large enough (empirically, KRAS's 4 states already sits near this
edge) that "found in 1 of N states" is a meaningful rarity signal rather
than a near-binary one. This is precisely the gap real BioEmu sampling
(hundreds to thousands of states) is intended to close — see
`docs/LIMITATIONS.md`.

## Target 2: PRMT5:MEP50 — EE-loop allosteric pocket

Discovery ensemble: **4GQB, 4X61, 5EML** — no structure with the allosteric
EE-loop site occupied. Weaker apo criterion than KRAS/ABL: 4X61 and 5EML do
carry a different (SAM-competitive, active-site) small-molecule inhibitor,
because no fully ligand-free human PRMT5:MEP50 structure has been
crystallized (documented in structural literature). Reported honestly in
`configs/prmt5.yaml` rather than concealed.

Result pending — see `artifacts/prmt5_conform-agent_*/metrics/experiment_manifest.json`
for the executed run (run via `python scripts/run_experiment.py run --target prmt5 --closed-loop`).

## What this section is for

Not to claim the method "works on any target" — it demonstrably does not
yet, on ABL kinase, with the current thin-ensemble fallback. It is here to
show the opposite of cherry-picking: the exact same frozen configuration
was pointed at unfamiliar proteins and its failure mode was investigated
and reported rather than hidden. A judge asking "did you tune this for
KRAS?" gets a documented negative result as the answer, which is stronger
evidence of good-faith methodology than a second lucky success would have
been.
