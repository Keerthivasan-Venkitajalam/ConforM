"""Permutation test: does the pocket-ranking algorithm beat chance at picking
the true cryptic pocket, on a single fixed structural ensemble?

This directly answers the causal question the shared-ensemble ablation
(evaluation/shared_ensemble_ablation.py) sets up but doesn't itself test:
given the exact same computed rank_score for every candidate pocket family,
is the family with the best ground-truth overlap ranked where it is because
the score is genuinely informative, or could a ranking with no real signal
have landed it there by chance?

Method (standard label-permutation / enrichment-significance test, the same
family of test used to assess virtual-screening enrichment factors):
  1. Take the REAL rank_score ordering of pocket families exactly as computed
     by the real, deterministic ranking algorithm -- never modified here.
  2. Find the real observed rank R_obs of the family with the best
     ground-truth residue overlap (the best geometric match available in
     this ensemble, wherever the algorithm placed it).
  3. Randomly reshuffle WHICH family carries the ground-truth-overlap
     values, N times. Each shuffle asks: "if ground truth had been handed
     out at random instead of reflecting real pocket geometry, how good a
     rank would the top-labeled family land under this same real score
     ordering, purely by chance?"
  4. Empirical p-value = fraction of shuffles at least as good as R_obs.

No p-value is pre-specified as a target; whatever comes out is reported.

Usage:
    python evaluation/permutation_test.py <experiment_dir> [--n-permutations 10000]
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def run_permutation_test(pocket_families: list[dict], n_permutations: int = 10000,
                          seed: int = 42) -> dict:
    n = len(pocket_families)
    if n < 2:
        raise ValueError(f"Need at least 2 pocket families for a meaningful test, got {n}")

    # Real ranks 1..n by the real, unmodified rank_score, descending.
    ordered = sorted(pocket_families, key=lambda f: f["rank_score"], reverse=True)
    real_rank_score = [f["rank_score"] for f in ordered]  # fixed throughout
    real_overlaps = [f["ground_truth_overlap"] for f in ordered]

    best_overlap = max(real_overlaps)
    # Real observed rank: best (lowest-numbered) rank among families tied for
    # the best ground-truth overlap actually available in this ensemble.
    r_obs = min(i + 1 for i, ov in enumerate(real_overlaps) if ov == best_overlap)

    rng = random.Random(seed)
    indices = list(range(n))
    null_ranks = []
    at_least_as_good = 0
    for _ in range(n_permutations):
        rng.shuffle(indices)
        shuffled_overlaps = [real_overlaps[i] for i in indices]
        perm_best = max(shuffled_overlaps)
        perm_rank = min(i + 1 for i, ov in enumerate(shuffled_overlaps) if ov == perm_best)
        null_ranks.append(perm_rank)
        if perm_rank <= r_obs:
            at_least_as_good += 1

    p_value = at_least_as_good / n_permutations

    return {
        "n_pocket_families": n,
        "n_permutations": n_permutations,
        "seed": seed,
        "real_observed_rank_of_best_overlap_family": r_obs,
        "best_ground_truth_overlap_available": best_overlap,
        "null_distribution_mean_rank": sum(null_ranks) / len(null_ranks),
        "null_distribution_median_rank": sorted(null_ranks)[len(null_ranks) // 2],
        "empirical_p_value": p_value,
        "interpretation": (
            f"Under {n_permutations} random relabelings of which pocket family carries the "
            f"ground-truth overlap, a rank as good as the real ranking's rank {r_obs} occurred "
            f"in {p_value * 100:.2f}% of permutations. This is a test of whether rank_score is "
            f"associated with ground-truth overlap beyond chance on THIS fixed ensemble; it does "
            f"not by itself establish cross-ensemble generalization (see the shared-ensemble "
            f"ablation table for that comparison)."
        ),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment_dir", help="experiment directory containing "
                                            "metrics/pocket_families.json")
    ap.add_argument("--n-permutations", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    families_path = Path(args.experiment_dir) / "metrics" / "pocket_families.json"
    families = json.loads(families_path.read_text())
    result = run_permutation_test(families, n_permutations=args.n_permutations, seed=args.seed)
    result["experiment_dir"] = str(args.experiment_dir)
    result["source"] = str(families_path)

    out_path = Path(args.experiment_dir) / "metrics" / "permutation_test.json"
    out_path.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    print(f"\nWrote {out_path}")
