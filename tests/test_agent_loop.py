"""Agent policy / state-machine tests using mock scientific data (no real tools)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from agent.policies import Policy, PolicyConfig
from agent.state import Action, ExperimentState, PocketCandidate, action_hash


def pocket(druggability=0.9):
    return PocketCandidate(state_pdb_id="7RPZ", pocket_index=1, volume=800.0,
                            druggability=druggability, residues=["HIS95", "TYR96"])


def test_first_action_is_generate_ensemble():
    s = ExperimentState(target="KRAS")
    assert Policy().decide(s).action is Action.GENERATE_ENSEMBLE


def test_progresses_through_setup_stages():
    p = Policy()
    s = ExperimentState(target="KRAS")
    s.ensemble = ["a.pdb", "b.pdb"]
    assert p.decide(s).action is Action.ANALYZE_ENSEMBLE
    s.conformational_states = {"max_rmsf": 5.0}
    assert p.decide(s).action is Action.FIND_POCKETS
    s.pocket_candidates = [pocket()]
    assert p.decide(s).action is Action.SELECT_POCKET
    s.selected_pocket = pocket()
    assert p.decide(s).action is Action.SCREEN_LIGANDS


def test_stops_when_max_iterations_reached():
    s = ExperimentState(target="KRAS")
    s.iteration = 5
    d = Policy(PolicyConfig(max_iterations=5)).decide(s)
    assert d.action is Action.STOP
    assert "Maximum iterations" in d.rationale


def test_stops_when_budget_exhausted():
    s = ExperimentState(target="KRAS", budget_seconds=10.0)
    s.consumed_seconds = 11.0
    d = Policy().decide(s)
    assert d.action is Action.STOP
    assert "budget" in d.rationale.lower()


def test_stops_when_no_druggable_pocket():
    s = ExperimentState(target="KRAS")
    s.ensemble = ["a.pdb"]
    s.conformational_states = {"x": 1}
    s.pocket_candidates = [pocket(druggability=0.01)]
    s.selected_pocket = pocket(druggability=0.01)
    d = Policy().decide(s)
    assert d.action is Action.STOP
    assert "no scientifically useful" in d.rationale.lower()


def test_optimizes_when_score_below_threshold():
    s = ExperimentState(target="KRAS")
    s.ensemble = ["a.pdb"]
    s.conformational_states = {"x": 1}
    s.pocket_candidates = [pocket()]
    s.selected_pocket = pocket()
    s.docking_results = [{"best_affinity_kcal": -7.0}]
    s.best_discovery_score = 0.5
    assert Policy(PolicyConfig(stop_score=0.85)).decide(s).action is Action.OPTIMIZE_LIGAND


def test_stops_when_confidence_threshold_met():
    s = ExperimentState(target="KRAS")
    s.ensemble = ["a.pdb"]
    s.conformational_states = {"x": 1}
    s.pocket_candidates = [pocket()]
    s.selected_pocket = pocket()
    s.docking_results = [{"best_affinity_kcal": -11.0}]
    s.best_discovery_score = 0.95
    d = Policy(PolicyConfig(stop_score=0.85)).decide(s)
    assert d.action is Action.STOP
    assert "threshold" in d.rationale.lower()


def test_does_not_repeat_optimization():
    s = ExperimentState(target="KRAS")
    s.ensemble = ["a.pdb"]
    s.conformational_states = {"x": 1}
    s.pocket_candidates = [pocket()]
    s.selected_pocket = pocket()
    s.docking_results = [{"best_affinity_kcal": -7.0}]
    s.best_discovery_score = 0.5
    s.record(Action.OPTIMIZE_LIGAND, {"pocket": s.selected_pocket.key}, {})
    assert Policy().decide(s).action is Action.STOP


def test_action_hash_is_deterministic_and_order_insensitive():
    a = action_hash(Action.FIND_POCKETS, {"x": 1, "y": 2})
    b = action_hash(Action.FIND_POCKETS, {"y": 2, "x": 1})
    c = action_hash(Action.FIND_POCKETS, {"x": 1, "y": 3})
    assert a == b
    assert a != c


def test_completed_actions_excludes_failures():
    s = ExperimentState(target="KRAS")
    s.record(Action.FIND_POCKETS, {"a": 1}, {}, failure="boom")
    s.record(Action.SELECT_POCKET, {"b": 2}, {})
    assert action_hash(Action.FIND_POCKETS, {"a": 1}) not in s.completed_actions()
    assert action_hash(Action.SELECT_POCKET, {"b": 2}) in s.completed_actions()
