"""Ground-truth leakage firewall.

This is the automated check referenced in docs/RESEARCH_CORRECTIONS.md #7.
It exists because the single easiest way for this project to produce a fake
result is to let `ground_truth_overlap` (how much a candidate pocket matches
the literature Switch-II residue list) influence which pocket gets selected
or how ligands get scored. If that ever happens, every "the agent discovered
X" claim collapses into "the agent was told the answer."

These tests do not merely inspect source text (a rename would defeat that).
They construct pocket families that are IDENTICAL on every blind descriptor
(druggability, novelty, volume) but DIFFERENT on ground_truth_overlap, and
assert the ranking function is completely insensitive to that difference.
If someone later wires ground_truth_overlap into the score, these tests fail.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import inspect

from agent.discovery_score import DiscoveryScoreInputs, compute_discovery_score
from agent.state import PocketCandidate
from pipelines.engines import rank_pocket_families


def _family(rank_key: str, ground_truth_overlap: float, druggability=0.8,
           novelty=0.9, volume=800.0) -> dict:
    rep = PocketCandidate(
        state_pdb_id=rank_key, pocket_index=1, volume=volume, druggability=druggability,
        residues=["ALA1", "GLY2"], ground_truth_overlap=ground_truth_overlap,
    )
    return {
        "representative": rep, "members": [rep], "n_states_present": 1,
        "persistence": 0.25, "baseline_volume": 0.0, "max_volume": volume,
        "cryptic_volume_gain": volume, "novelty": novelty, "max_druggability": druggability,
    }


def test_ranking_identical_for_zero_and_perfect_ground_truth_overlap():
    """Two families, blind-identical, differing ONLY in ground_truth_overlap
    (0.0 vs 1.0) -- ranking and scores must be byte-identical."""
    zero_gt = [_family("stateA", ground_truth_overlap=0.0)]
    full_gt = [_family("stateA", ground_truth_overlap=1.0)]
    ranked_zero = rank_pocket_families(zero_gt)
    ranked_full = rank_pocket_families(full_gt)
    assert ranked_zero[0]["rank_score"] == ranked_full[0]["rank_score"]


def test_low_ground_truth_pocket_can_outrank_high_ground_truth_pocket():
    """A pocket with WORSE ground-truth overlap but BETTER blind descriptors
    must win. If ground truth were leaking in, the perfect-overlap pocket
    with weaker geometry would win instead."""
    families = [
        _family("weak_geometry_perfect_match", ground_truth_overlap=1.0,
               druggability=0.3, novelty=0.2, volume=200.0),
        _family("strong_geometry_no_match", ground_truth_overlap=0.0,
               druggability=0.9, novelty=0.95, volume=900.0),
    ]
    ranked = rank_pocket_families(families)
    assert ranked[0]["representative"].state_pdb_id == "strong_geometry_no_match"


def test_rank_pocket_families_signature_excludes_ground_truth():
    """Static guard: the ranking function must not even accept a ground-truth
    parameter. Catches a future refactor that threads it in via new kwargs."""
    sig = inspect.signature(rank_pocket_families)
    for name in sig.parameters:
        assert "ground_truth" not in name.lower(), (
            f"rank_pocket_families gained a ground-truth parameter: {name}")


def test_discovery_score_inputs_have_no_ground_truth_field():
    """Static guard on the Discovery Score dataclass itself."""
    for field_name in DiscoveryScoreInputs.__dataclass_fields__:
        assert "ground_truth" not in field_name.lower(), (
            f"DiscoveryScoreInputs gained a ground-truth field: {field_name}")


def test_discovery_score_unaffected_by_ground_truth_regardless_of_pocket_object():
    """compute_discovery_score takes no pocket/family object at all -- only
    scalar geometric/binding/ligand quantities -- so ground truth structurally
    cannot reach it. This test locks that call signature."""
    sig = inspect.signature(compute_discovery_score)
    param_names = list(sig.parameters)
    assert param_names[0] == "inputs"
    for name in DiscoveryScoreInputs.__dataclass_fields__:
        assert "ground_truth" not in name.lower()


def test_pocket_candidate_carries_ground_truth_but_find_pockets_does_not_rank_on_it():
    """PocketCandidate legitimately STORES ground_truth_overlap for post-hoc
    evaluation (evaluation/metrics.py) -- storing it is fine and necessary.
    This test documents that boundary: stored != used-for-ranking, which the
    two tests above already enforce functionally."""
    p = PocketCandidate(state_pdb_id="x", pocket_index=1, volume=1.0,
                        druggability=1.0, residues=[], ground_truth_overlap=0.73)
    assert p.ground_truth_overlap == 0.73  # stored for evaluation/, not scoring
