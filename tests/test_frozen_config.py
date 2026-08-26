"""Enforces the generalization-test pre-registration claim in
docs/GENERALIZATION.md: every held-out target's pocket_detection/docking/
discovery_score/agent sections must be byte-identical to KRAS's. If someone
edits one target's weights without the others, this test catches it -- the
whole point of the generalization result is that nothing was tuned per
target.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import yaml

FROZEN_SECTIONS = ["pocket_detection", "docking", "discovery_score", "agent"]
CONFIGS_DIR = Path(__file__).resolve().parent.parent / "configs"


def _load(name: str) -> dict:
    return yaml.safe_load((CONFIGS_DIR / name).read_text())


def test_abl_kinase_config_matches_kras_on_frozen_sections():
    kras = _load("kras_g12d.yaml")
    abl = _load("abl_kinase.yaml")
    for section in FROZEN_SECTIONS:
        assert abl[section] == kras[section], f"abl_kinase.yaml differs from kras_g12d.yaml in '{section}'"


def test_prmt5_config_matches_kras_on_frozen_sections():
    kras = _load("kras_g12d.yaml")
    prmt5 = _load("prmt5.yaml")
    for section in FROZEN_SECTIONS:
        assert prmt5[section] == kras[section], f"prmt5.yaml differs from kras_g12d.yaml in '{section}'"


def test_all_targets_registered_in_cli():
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import run_experiment
    assert set(run_experiment.TARGETS) == {"kras-g12d", "abl-kinase", "prmt5"}
    for name, path in run_experiment.TARGETS.items():
        assert (Path(__file__).resolve().parent.parent / path).exists(), f"{name} config missing: {path}"


def test_same_ligand_library_used_for_every_target():
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import run_experiment
    # DEFAULT_LIGANDS is target-agnostic by construction (module-level constant
    # used as the argparse default regardless of --target); this test locks
    # that design decision against a future per-target special-case.
    assert run_experiment.DEFAULT_LIGANDS.name == "ligands_kras.csv"
