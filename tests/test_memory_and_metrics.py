"""Scientific memory (duplicate detection, persistence) and evaluation metrics."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest

from agent.state import PocketCandidate
from db.repository import Repository
from evaluation.metrics import enrichment, pocket_residue_recovery, summarize_run
from pipelines.engines import cluster_pocket_families, ground_truth_overlap, rank_pocket_families


@pytest.fixture
def repo(tmp_path):
    r = Repository(tmp_path / "mem.db")
    yield r
    r.close()


def test_duplicate_detection(repo):
    repo.create_experiment("e1", "KRAS", "test", {})
    assert not repo.has_completed("e1", "h1")
    repo.log_step("e1", 0, "FIND_POCKETS", "h1")
    assert repo.has_completed("e1", "h1")


def test_failed_step_is_not_treated_as_completed(repo):
    repo.create_experiment("e1", "KRAS", "test", {})
    repo.log_step("e1", 0, "FIND_POCKETS", "h1", failure="tool crashed")
    assert not repo.has_completed("e1", "h1")


def test_steps_are_retrievable(repo):
    repo.create_experiment("e1", "KRAS", "test", {})
    repo.log_step("e1", 0, "GENERATE_ENSEMBLE", "h0", interpretation="made 4 states")
    repo.log_step("e1", 1, "FIND_POCKETS", "h1")
    steps = repo.steps("e1")
    assert len(steps) == 2
    assert steps[0]["scientific_interpretation"] == "made 4 states"


def test_docking_results_persist(repo):
    repo.create_experiment("e1", "KRAS", "test", {})
    repo.save_docking("e1", [{"ligand_name": "benzene", "best_affinity_kcal": -5.0,
                               "engine": "vina", "discovery_score": 0.4}])
    rows = repo.docking_results("e1")
    assert len(rows) == 1 and rows[0]["ligand_name"] == "benzene"


def test_pocket_residue_recovery_perfect():
    r = pocket_residue_recovery(["HIS95", "TYR96", "GLN99"], ["H95", "Y96", "Q99"])
    assert r["recall"] == 1.0


def test_pocket_residue_recovery_none():
    r = pocket_residue_recovery(["ALA11", "GLY13"], ["H95", "Y96"])
    assert r["recall"] == 0.0


def test_ground_truth_overlap_matches_by_residue_number():
    assert ground_truth_overlap(["HIS95", "TYR96"], ["H95", "Y96"]) == 1.0
    assert ground_truth_overlap(["ALA11"], ["H95", "Y96"]) == 0.0


def test_enrichment_finds_active_in_top_n():
    results = [{"ligand_name": "active", "best_affinity_kcal": -10.0},
               {"ligand_name": "decoy1", "best_affinity_kcal": -5.0},
               {"ligand_name": "decoy2", "best_affinity_kcal": -4.0}]
    e = enrichment(results, {"active"}, top_n=1)
    assert e["actives_in_top_n"] == 1


def _pc(pdb, idx, vol, drug, residues):
    return PocketCandidate(state_pdb_id=pdb, pocket_index=idx, volume=vol,
                            druggability=drug, residues=residues)


def test_cryptic_family_absent_from_baseline_has_novelty_one():
    # A cavity that appears only in a non-baseline state is fully cryptic.
    candidates = [
        _pc("BASE", 1, 900.0, 0.5, ["ALA11", "GLY13", "LYS16"]),
        _pc("OTHER", 1, 800.0, 0.9, ["HIS95", "TYR96", "GLN99"]),
    ]
    fams = cluster_pocket_families(candidates, "BASE", n_states=2)
    cryptic = next(f for f in fams if f["representative"].state_pdb_id == "OTHER")
    assert cryptic["baseline_volume"] == 0.0
    assert cryptic["novelty"] == 1.0
    assert cryptic["persistence"] == 0.5


def test_pocket_present_in_baseline_has_low_novelty():
    candidates = [
        _pc("BASE", 1, 900.0, 0.8, ["ALA11", "GLY13", "LYS16"]),
        _pc("OTHER", 1, 950.0, 0.8, ["ALA11", "GLY13", "LYS16"]),
    ]
    fams = cluster_pocket_families(candidates, "BASE", n_states=2)
    assert len(fams) == 1
    assert fams[0]["novelty"] < 0.1
    assert fams[0]["persistence"] == 1.0


def test_ranking_prefers_cryptic_over_always_open_pocket():
    """The whole point of the system: an always-open cavity must not outrank
    a comparably druggable cavity that is absent from the apo baseline."""
    candidates = [
        _pc("BASE", 1, 1000.0, 0.85, ["ALA11", "GLY13", "LYS16"]),
        _pc("OTHER", 1, 1000.0, 0.85, ["ALA11", "GLY13", "LYS16"]),   # always open
        _pc("OTHER", 2, 850.0, 0.88, ["HIS95", "TYR96", "GLN99"]),    # cryptic
    ]
    fams = cluster_pocket_families(candidates, "BASE", n_states=2)
    ranked = rank_pocket_families(fams)
    assert ranked[0]["representative"].residues == ["HIS95", "TYR96", "GLN99"]


def test_summarize_run_handles_empty_manifest():
    s = summarize_run({}, ["H95"])
    assert s["cryptic_residue_recall"] == 0.0
