import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tools.mdpocket_tool import _parse_info_txt, _parse_pocket_residues

FAKE_INFO_TXT = """Pocket 1 :
\tScore : \t0.702
\tDruggability Score : \t0.834
\tNumber of Alpha Spheres : \t145
\tVolume : \t1088.510

Pocket 2 :
\tScore : \t0.349
\tDruggability Score : \t0.308
\tNumber of Alpha Spheres : \t85
\tVolume : \t770.683
"""

FAKE_ATM_PDB = """HEADER
ATOM      1  CA  HIS A  95      10.000  20.000  30.000  1.00  0.00           C
ATOM      2  CA  TYR A  96      11.000  21.000  31.000  1.00  0.00           C
ATOM      3  CA  GLN A  99      12.000  22.000  32.000  1.00  0.00           C
"""


def test_parse_info_txt_extracts_both_pockets(tmp_path):
    p = tmp_path / "info.txt"
    p.write_text(FAKE_INFO_TXT)
    pockets = _parse_info_txt(p)
    assert len(pockets) == 2
    assert pockets[0]["Volume"] == 1088.510
    assert pockets[0]["Druggability Score"] == 0.834
    assert pockets[1]["Score"] == 0.349


def test_parse_pocket_residues(tmp_path):
    p = tmp_path / "pocket1_atm.pdb"
    p.write_text(FAKE_ATM_PDB)
    residues = _parse_pocket_residues(p)
    assert residues == ["GLN99", "HIS95", "TYR96"]
