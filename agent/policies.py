"""Deterministic decision policy for the closed-loop agent.

This is the agent's "what should I do next" reasoning, implemented as an
explicit, auditable state machine rather than an LLM call. An LLM may later
be layered on top to *propose* actions, but this policy remains the arbiter
so that no scientific quantity or experiment selection is hallucinated
(master prompt rules #17, #18).

Termination conditions (rule #18):
  - confidence threshold reached (Discovery Score >= stop_score)
  - budget exhausted
  - MAX_ITERATIONS reached
  - no scientifically useful next action exists
"""
from __future__ import annotations

from dataclasses import dataclass

from agent.state import Action, ExperimentState, action_hash


@dataclass
class PolicyConfig:
    max_iterations: int = 5
    stop_score: float = 0.85
    min_pocket_druggability: float = 0.20
    optimize_if_best_affinity_worse_than: float = -9.0


@dataclass
class Decision:
    action: Action
    params: dict
    rationale: str


class Policy:
    def __init__(self, config: PolicyConfig | None = None):
        self.config = config or PolicyConfig()

    def decide(self, state: ExperimentState) -> Decision:
        cfg = self.config

        if state.iteration >= cfg.max_iterations:
            return Decision(Action.STOP, {}, f"Maximum iterations ({cfg.max_iterations}) reached.")

        if state.budget_remaining() <= 0:
            return Decision(Action.STOP, {},
                            f"Compute budget exhausted ({state.consumed_seconds:.0f}s used).")

        if not state.ensemble:
            return Decision(Action.GENERATE_ENSEMBLE, {"target": state.target},
                            "No conformational ensemble exists yet; a state-space sample is the "
                            "prerequisite for transient pocket detection.")

        if not state.conformational_states:
            return Decision(Action.ANALYZE_ENSEMBLE, {"n_states": len(state.ensemble)},
                            "Ensemble exists but has not been characterized; need RMSD/RMSF/PCA "
                            "to describe the sampled state space.")

        if not state.pocket_candidates:
            return Decision(Action.FIND_POCKETS, {"n_states": len(state.ensemble)},
                            "State space characterized but no cavity detection has been run.")

        if state.selected_pocket is None:
            return Decision(Action.SELECT_POCKET, {"n_candidates": len(state.pocket_candidates)},
                            "Pocket candidates exist; rank them and commit to the most druggable "
                            "cavity before spending docking compute.")

        if state.selected_pocket.druggability < cfg.min_pocket_druggability:
            return Decision(Action.STOP, {},
                            f"Best available pocket has druggability "
                            f"{state.selected_pocket.druggability:.2f} < "
                            f"{cfg.min_pocket_druggability}; no scientifically useful ligand "
                            "experiment remains for this ensemble.")

        if not state.docking_results:
            return Decision(Action.SCREEN_LIGANDS, {"pocket": state.selected_pocket.key},
                            "A druggable cavity is selected but no ligand evidence exists yet; "
                            "screen the library to establish a baseline affinity distribution.")

        if state.best_discovery_score >= cfg.stop_score:
            return Decision(Action.STOP, {},
                            f"Confidence threshold met: Discovery Score "
                            f"{state.best_discovery_score:.3f} >= {cfg.stop_score}.")

        best_affinity = min((d["best_affinity_kcal"] for d in state.docking_results
                             if d.get("best_affinity_kcal") is not None), default=0.0)
        already_optimized = action_hash(Action.OPTIMIZE_LIGAND, {"pocket": state.selected_pocket.key}) \
            in state.completed_actions()

        if not already_optimized:
            return Decision(
                Action.OPTIMIZE_LIGAND, {"pocket": state.selected_pocket.key},
                f"Best library affinity is {best_affinity:.2f} kcal/mol and Discovery Score "
                f"{state.best_discovery_score:.3f} is below threshold; generate and re-dock "
                "analogs of the best scaffold to test whether affinity can be improved.")

        return Decision(Action.STOP, {},
                        "Library screening and one optimization round are complete; no further "
                        "distinct, non-duplicate experiment is available under the current policy.")
