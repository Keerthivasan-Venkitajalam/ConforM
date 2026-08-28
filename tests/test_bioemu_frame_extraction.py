"""Regression test for the BioEmu trajectory-collapse bug.

A real GPU run on 2026-08-28 showed every ablation mode reporting n_states=1
despite requesting num_samples=100. Root cause: bioemu.sample.main() writes
its ensemble as samples.xtc + topology.pdb (a trajectory), and
BioEmuProvider.generate() was only globbing for *.pdb / frame_*.pdb /
samples_*.pdb -- none of which match that format, so it silently collapsed
the whole ensemble down to the single reference topology frame. This does
not require bioemu or a GPU to test: it only needs MDAnalysis (a base CPU
dependency) to write a synthetic multi-frame trajectory and confirm
BioEmuProvider._extract_frames() recovers every frame, not just one.
"""
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tools.bioemu_tool import BioEmuProvider

N_ATOMS = 12
N_FRAMES = 7


def _write_synthetic_trajectory(out_dir: Path) -> tuple[Path, Path]:
    universe = mda.Universe.empty(N_ATOMS, trajectory=True)
    universe.add_TopologyAttr("name", [f"C{i}" for i in range(N_ATOMS)])
    universe.add_TopologyAttr("resname", ["MOL"])
    universe.add_TopologyAttr("resid", [1])
    topology_path = out_dir / "topology.pdb"
    universe.atoms.positions = np.zeros((N_ATOMS, 3))
    universe.atoms.write(str(topology_path))

    traj_path = out_dir / "samples.xtc"
    with mda.Writer(str(traj_path), n_atoms=N_ATOMS) as writer:
        rng = np.random.default_rng(seed=42)
        for _ in range(N_FRAMES):
            universe.atoms.positions = rng.random((N_ATOMS, 3)).astype(np.float32)
            writer.write(universe.atoms)
    return traj_path, topology_path


def test_extract_frames_recovers_full_ensemble_not_just_topology(tmp_path):
    traj_path, topology_path = _write_synthetic_trajectory(tmp_path)

    frames = BioEmuProvider._extract_frames(traj_path, topology_path, tmp_path)

    assert len(frames) == N_FRAMES, (
        f"expected {N_FRAMES} states extracted from the trajectory, got "
        f"{len(frames)} -- this is the exact collapse-to-1 bug if it regresses"
    )
    assert all(p.exists() for p in frames)


def test_extract_frames_produces_distinct_coordinates(tmp_path):
    traj_path, topology_path = _write_synthetic_trajectory(tmp_path)
    frames = BioEmuProvider._extract_frames(traj_path, topology_path, tmp_path)

    coords = [mda.Universe(str(p)).atoms.positions.copy() for p in frames]
    assert not np.allclose(coords[0], coords[1]), (
        "extracted frames must have distinct coordinates -- identical frames "
        "would mean the trajectory wasn't really being read frame-by-frame"
    )
