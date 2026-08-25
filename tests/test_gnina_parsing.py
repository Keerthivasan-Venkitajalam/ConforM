"""Unit tests for GNINA SDF-output parsing (tools/gnina_tool.py).

The GNINA binary itself needs a CUDA GPU and cannot run in this environment
(see docs/LIMITATIONS.md), so this only tests the parsing logic against a
hand-written SDF fixture using GNINA's documented property tags (CNNscore,
CNNaffinity, CNNvariance, minimizedAffinity). The full rescore() path is
exercised for real on rented GPU compute via scripts/gpu_session.sh.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tools.gnina_tool import _parse_sdf_poses, gnina_available

FIXTURE_SDF = """benzene_pose1
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.3968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2098    0.6984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2098   -0.6984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.3968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2098   -0.6984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2098    0.6984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
>  <CNNscore>
0.812

>  <CNNaffinity>
6.34

>  <CNNvariance>
0.021

>  <minimizedAffinity>
-7.85

$$$$
benzene_pose2
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.3968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2098    0.6984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2098   -0.6984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.3968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2098   -0.6984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2098    0.6984    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
>  <CNNscore>
0.553

>  <CNNaffinity>
5.90

>  <CNNvariance>
0.034

>  <minimizedAffinity>
-6.20

$$$$
"""


def test_parse_sdf_poses_extracts_all_properties(tmp_path):
    sdf = tmp_path / "test.sdf"
    sdf.write_text(FIXTURE_SDF)
    poses = _parse_sdf_poses(sdf)
    assert len(poses) == 2
    assert poses[0].cnn_score == 0.812
    assert poses[0].cnn_affinity == 6.34
    assert poses[0].cnn_variance == 0.021
    assert poses[0].vina_affinity_kcal == -7.85


def test_parse_sdf_poses_second_pose_has_lower_score(tmp_path):
    sdf = tmp_path / "test.sdf"
    sdf.write_text(FIXTURE_SDF)
    poses = _parse_sdf_poses(sdf)
    assert poses[1].cnn_score < poses[0].cnn_score


def test_gnina_available_returns_bool():
    # No GPU/binary on this dev machine -- must return False cleanly, not raise.
    assert gnina_available() is False
