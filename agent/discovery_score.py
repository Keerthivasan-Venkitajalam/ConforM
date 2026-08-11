"""Deterministic Discovery Score (research plan, Part 8).

DiscoveryScore =
      w_pocket   * pocket_novelty
    + w_volume   * normalized_volume
    + w_binding  * normalized_binding_score
    + w_state    * state_novelty
    + w_ligand   * ligand_quality
    - w_invalid  * structural_penalty
    - w_cost     * computational_cost

All inputs are normalized to [0, 1] before weighting. This function is pure,
deterministic Python -- it is the ONLY place a Discovery Score may be
computed. No LLM/agent code may invent or override this value (research plan
Part 17 / master prompt rule #17).
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class DiscoveryScoreInputs:
    pocket_volume: float              # Angstrom^3, from fpocket
    max_observed_volume: float        # Angstrom^3, normalization reference across candidate pockets
    pocket_druggability: float        # fpocket druggability score, already in [0, 1]
    state_frequency: float            # fractional population of the conformational state (uniform-fallback if not equilibrium)
    binding_affinity_kcal: float      # Vina best affinity (more negative = better)
    best_possible_affinity_kcal: float  # normalization reference, e.g. best score seen this run
    worst_possible_affinity_kcal: float  # normalization reference, e.g. weakest score seen this run
    ligand_qed: float                 # RDKit QED in [0, 1]
    lipinski_violations: int          # 0-4
    structural_clash_penalty: float = 0.0   # in [0, 1]; 0 = no fallback-provider structural validation performed
    computational_cost_normalized: float = 0.0  # in [0, 1]


@dataclass
class DiscoveryScoreWeights:
    w_pocket_novelty: float = 0.30
    w_volume: float = 0.15
    w_binding: float = 0.30
    w_state_novelty: float = 0.15
    w_ligand_quality: float = 0.10
    w_invalid_penalty: float = 0.20
    w_cost_penalty: float = 0.05


def _clip01(x: float) -> float:
    return max(0.0, min(1.0, x))


def normalized_binding_score(affinity: float, best: float, worst: float) -> float:
    """Map Vina affinity (more negative = better) onto [0, 1], 1 = best seen."""
    if best == worst:
        return 0.5
    return _clip01((worst - affinity) / (worst - best))


def compute_discovery_score(inputs: DiscoveryScoreInputs,
                             weights: DiscoveryScoreWeights = DiscoveryScoreWeights()) -> dict:
    normalized_volume = _clip01(inputs.pocket_volume / inputs.max_observed_volume) \
        if inputs.max_observed_volume > 0 else 0.0
    pocket_novelty = _clip01(inputs.pocket_druggability)
    state_novelty = _clip01(1.0 - inputs.state_frequency)
    binding_norm = normalized_binding_score(
        inputs.binding_affinity_kcal, inputs.best_possible_affinity_kcal, inputs.worst_possible_affinity_kcal)
    lipinski_penalty = inputs.lipinski_violations / 4.0
    ligand_quality = _clip01(inputs.ligand_qed * (1.0 - 0.5 * lipinski_penalty))
    structural_penalty = _clip01(inputs.structural_clash_penalty)
    cost_penalty = _clip01(inputs.computational_cost_normalized)

    score = (
        weights.w_pocket_novelty * pocket_novelty
        + weights.w_volume * normalized_volume
        + weights.w_binding * binding_norm
        + weights.w_state_novelty * state_novelty
        + weights.w_ligand_quality * ligand_quality
        - weights.w_invalid_penalty * structural_penalty
        - weights.w_cost_penalty * cost_penalty
    )

    return {
        "discovery_score": round(score, 4),
        "components": {
            "pocket_novelty": round(pocket_novelty, 4),
            "normalized_volume": round(normalized_volume, 4),
            "binding_norm": round(binding_norm, 4),
            "state_novelty": round(state_novelty, 4),
            "ligand_quality": round(ligand_quality, 4),
            "structural_penalty": round(structural_penalty, 4),
            "cost_penalty": round(cost_penalty, 4),
        },
        "range": "theoretical range approx [-0.25, 1.0] given default weights; not a probability",
    }
