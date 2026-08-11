import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tools.rdkit_tool import validate_and_prepare


def test_valid_smiles_passes():
    r = validate_and_prepare("ethanol", "CCO")
    assert r.valid
    assert r.canonical_smiles is not None
    assert r.mol_weight is not None and r.mol_weight > 0


def test_invalid_smiles_is_rejected():
    r = validate_and_prepare("garbage", "not_a_smiles!!!(((")
    assert not r.valid
    assert r.error is not None


def test_lipinski_violations_for_large_molecule():
    # A molecule engineered to be well outside Lipinski bounds.
    huge_smiles = "C" * 60
    r = validate_and_prepare("huge_alkane", huge_smiles)
    assert r.valid
    assert r.lipinski_violations >= 1


def test_3d_embedding_writes_pdb(tmp_path):
    r = validate_and_prepare("benzene", "c1ccccc1", tmp_path)
    assert r.embedded_3d
    assert r.pdb_path is not None
    assert Path(r.pdb_path).exists()
