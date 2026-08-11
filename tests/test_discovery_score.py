import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from agent.discovery_score import DiscoveryScoreInputs, DiscoveryScoreWeights, compute_discovery_score, normalized_binding_score


def test_normalized_binding_score_best_is_one():
    assert normalized_binding_score(-10.0, best=-10.0, worst=-2.0) == 1.0


def test_normalized_binding_score_worst_is_zero():
    assert normalized_binding_score(-2.0, best=-10.0, worst=-2.0) == 0.0


def test_normalized_binding_score_handles_equal_best_worst():
    assert normalized_binding_score(-5.0, best=-5.0, worst=-5.0) == 0.5


def test_discovery_score_is_deterministic():
    inputs = DiscoveryScoreInputs(
        pocket_volume=800.0, max_observed_volume=1000.0, pocket_druggability=0.9,
        state_frequency=0.25, binding_affinity_kcal=-9.6,
        best_possible_affinity_kcal=-9.6, worst_possible_affinity_kcal=-2.8,
        ligand_qed=0.6, lipinski_violations=0,
    )
    r1 = compute_discovery_score(inputs)
    r2 = compute_discovery_score(inputs)
    assert r1 == r2
    assert 0.0 <= r1["discovery_score"] <= 1.0


def test_discovery_score_penalizes_structural_invalidity():
    base = DiscoveryScoreInputs(
        pocket_volume=800.0, max_observed_volume=1000.0, pocket_druggability=0.9,
        state_frequency=0.25, binding_affinity_kcal=-9.6,
        best_possible_affinity_kcal=-9.6, worst_possible_affinity_kcal=-2.8,
        ligand_qed=0.6, lipinski_violations=0,
    )
    penalized = DiscoveryScoreInputs(**{**base.__dict__, "structural_clash_penalty": 1.0})
    s_base = compute_discovery_score(base)["discovery_score"]
    s_penalized = compute_discovery_score(penalized)["discovery_score"]
    assert s_penalized < s_base


def test_weights_are_configurable():
    w = DiscoveryScoreWeights(w_binding=0.9, w_pocket_novelty=0.0, w_volume=0.0,
                               w_state_novelty=0.0, w_ligand_quality=0.0,
                               w_invalid_penalty=0.0, w_cost_penalty=0.0)
    inputs = DiscoveryScoreInputs(
        pocket_volume=0.0, max_observed_volume=1.0, pocket_druggability=0.0,
        state_frequency=1.0, binding_affinity_kcal=-9.6,
        best_possible_affinity_kcal=-9.6, worst_possible_affinity_kcal=-2.8,
        ligand_qed=0.0, lipinski_violations=4,
    )
    result = compute_discovery_score(inputs, w)
    assert abs(result["discovery_score"] - 0.9) < 1e-6
