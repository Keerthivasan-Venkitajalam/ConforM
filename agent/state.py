"""Explicit experiment state and action vocabulary for the closed-loop agent.

The state object is the single source of truth about what the agent knows.
It is serializable so it can be persisted to (and rehydrated from) the
scientific-memory database in `db/repository.py`.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any


class Action(str, Enum):
    GENERATE_ENSEMBLE = "GENERATE_ENSEMBLE"
    ANALYZE_ENSEMBLE = "ANALYZE_ENSEMBLE"
    FIND_POCKETS = "FIND_POCKETS"
    SELECT_POCKET = "SELECT_POCKET"
    SCREEN_LIGANDS = "SCREEN_LIGANDS"
    OPTIMIZE_LIGAND = "OPTIMIZE_LIGAND"
    VALIDATE_POSE = "VALIDATE_POSE"
    STOP = "STOP"


@dataclass
class PocketCandidate:
    state_pdb_id: str
    pocket_index: int
    volume: float
    druggability: float
    residues: list[str]
    ground_truth_overlap: float = 0.0

    @property
    def key(self) -> str:
        return f"{self.state_pdb_id}:pocket{self.pocket_index}"


@dataclass
class LigandRecord:
    name: str
    smiles: str
    qed: float | None = None
    lipinski_violations: int | None = None
    best_affinity_kcal: float | None = None
    discovery_score: float | None = None
    origin: str = "library"        # "library" | "optimized"
    parent: str | None = None


@dataclass
class ExperimentState:
    target: str
    baseline_structure: str | None = None
    ensemble: list[str] = field(default_factory=list)          # paths
    ensemble_provider: str | None = None
    ensemble_is_equilibrium: bool = False
    conformational_states: dict[str, Any] = field(default_factory=dict)
    pocket_candidates: list[PocketCandidate] = field(default_factory=list)
    selected_pocket: PocketCandidate | None = None
    ligand_candidates: list[LigandRecord] = field(default_factory=list)
    docking_results: list[dict] = field(default_factory=list)
    optimization_results: list[dict] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)
    # `iteration` counts closed-loop EXPERIMENT cycles (ligand screening /
    # optimization / pose validation) -- the things that consume real compute
    # and generate new evidence. Setup stages (ensemble generation, analysis,
    # pocket detection/selection) are one-time pipeline prerequisites and are
    # tracked by `step_count` instead, so they do not burn the iteration
    # budget that MAX_ITERATIONS is meant to bound.
    iteration: int = 0
    step_count: int = 0
    budget_seconds: float = 3600.0
    consumed_seconds: float = 0.0
    history: list[dict] = field(default_factory=list)
    best_discovery_score: float = float("-inf")

    def budget_remaining(self) -> float:
        return self.budget_seconds - self.consumed_seconds

    def record(self, action: Action, params: dict, outcome: dict, failure: str | None = None):
        self.history.append({
            "iteration": self.iteration,
            "action": action.value,
            "params": params,
            "input_hash": action_hash(action, params),
            "outcome": outcome,
            "failure": failure,
        })

    def completed_actions(self) -> set[str]:
        return {h["input_hash"] for h in self.history if h["failure"] is None}

    def to_dict(self) -> dict:
        d = asdict(self)
        d["selected_pocket"] = asdict(self.selected_pocket) if self.selected_pocket else None
        d["best_discovery_score"] = (
            None if self.best_discovery_score == float("-inf") else self.best_discovery_score
        )
        return d


def action_hash(action: Action, params: dict) -> str:
    """Deterministic key for duplicate-experiment detection."""
    payload = json.dumps({"action": action.value, "params": params}, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]
