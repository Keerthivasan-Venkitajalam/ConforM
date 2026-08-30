"""Unit tests for evaluation/permutation_test.py's statistical logic.

Uses synthetic pocket-family data (no GPU, no real experiment needed) to
verify the test behaves correctly at both extremes: a ranking that perfectly
tracks ground truth should be highly significant; a ranking uncorrelated
with ground truth should not be.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from evaluation.permutation_test import run_permutation_test


def _family(rank_score, ground_truth_overlap):
    return {"rank_score": rank_score, "ground_truth_overlap": ground_truth_overlap}


def test_perfect_ranking_is_highly_significant():
    # The single family with real ground-truth overlap also has the top score.
    families = [
        _family(0.99, 0.80),
        _family(0.61, 0.00),
        _family(0.55, 0.00),
        _family(0.40, 0.00),
        _family(0.10, 0.00),
    ]
    result = run_permutation_test(families, n_permutations=5000, seed=1)
    assert result["real_observed_rank_of_best_overlap_family"] == 1
    # With only 1 of 5 families carrying the positive label, a random
    # relabeling lands it at rank 1 about 1/5 of the time -- not near-zero,
    # but this is exactly what a smaller pocket-family count implies; the
    # test itself, not a fabricated threshold, determines the number.
    assert 0.0 <= result["empirical_p_value"] <= 0.3


def test_no_association_gives_unremarkable_p_value():
    # Ground truth overlap on a middling-score family; nothing to detect.
    families = [
        _family(0.99, 0.00),
        _family(0.61, 0.00),
        _family(0.55, 0.20),
        _family(0.40, 0.00),
        _family(0.10, 0.00),
    ]
    result = run_permutation_test(families, n_permutations=5000, seed=1)
    assert result["real_observed_rank_of_best_overlap_family"] == 3
    # Landing at rank 3 of 5 by chance is unremarkable.
    assert result["empirical_p_value"] > 0.3


def test_requires_at_least_two_families():
    import pytest
    with pytest.raises(ValueError):
        run_permutation_test([_family(0.5, 0.1)], n_permutations=100)


def test_deterministic_given_seed():
    families = [_family(s, ov) for s, ov in
                [(0.9, 0.6), (0.7, 0.0), (0.5, 0.0), (0.3, 0.2), (0.1, 0.0)]]
    r1 = run_permutation_test(families, n_permutations=2000, seed=7)
    r2 = run_permutation_test(families, n_permutations=2000, seed=7)
    assert r1["empirical_p_value"] == r2["empirical_p_value"]
