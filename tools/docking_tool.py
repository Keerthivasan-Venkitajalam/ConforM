"""Docking wrapper: real AutoDock Vina (Python bindings, conda-forge build).

Tiered docking per the research plan (fast Vina screen -> GNINA CNN rescore
-> DiffDock-Pocket refinement). GNINA needs a CUDA build of its forked Caffe
backend and DiffDock-Pocket needs a GPU diffusion model; neither is
available in this CPU-only environment (see docs/LIMITATIONS.md). Only the
Vina tier is real/executed here. Scores below are unmodified Vina empirical
scoring-function output -- never fabricated or rescaled to look like a CNN
score.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from vina import Vina


@dataclass
class DockingResult:
    ligand_name: str
    receptor_pdb_id: str
    pocket_index: int
    engine: str
    center: list[float]
    box_size: list[float]
    poses_affinity_kcal_per_mol: list[float]
    best_affinity_kcal_per_mol: float


def prepare_receptor_pdbqt(receptor_pdb: Path, work_dir: Path) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    protein_only = work_dir / f"{receptor_pdb.stem}_protein.pdb"
    with open(receptor_pdb) as fin, open(protein_only, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "TER", "END")):
                fout.write(line)
    pdbqt_path = work_dir / f"{receptor_pdb.stem}.pdbqt"
    if not pdbqt_path.exists():
        subprocess.run(["obabel", str(protein_only), "-O", str(pdbqt_path), "-xr"],
                        check=True, capture_output=True, timeout=120)
    return pdbqt_path


def prepare_ligand_pdbqt(ligand_pdb: Path, work_dir: Path) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_path = work_dir / f"{ligand_pdb.stem}.pdbqt"
    if not pdbqt_path.exists():
        subprocess.run(["obabel", str(ligand_pdb), "-O", str(pdbqt_path)],
                        check=True, capture_output=True, timeout=60)
    return pdbqt_path


def pocket_centroid(pocket_vert_pqr: Path) -> np.ndarray:
    coords = []
    for line in pocket_vert_pqr.read_text().splitlines():
        if line.startswith("ATOM"):
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not coords:
        raise RuntimeError(f"No alpha-sphere coordinates found in {pocket_vert_pqr}")
    return np.array(coords).mean(axis=0)


def dock(receptor_pdbqt: Path, ligand_pdbqt: Path, center: np.ndarray,
         box_size=(24.0, 24.0, 24.0), exhaustiveness=8, n_poses=5,
         ligand_name="ligand", receptor_pdb_id="receptor", pocket_index=0) -> DockingResult:
    if shutil.which("obabel") is None:
        raise RuntimeError("Open Babel not found; required for PDBQT preparation.")
    v = Vina(sf_name="vina", cpu=4, seed=42, verbosity=0)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=list(map(float, center)), box_size=list(box_size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    energies = v.energies()
    affinities = [float(row[0]) for row in energies]
    return DockingResult(
        ligand_name=ligand_name, receptor_pdb_id=receptor_pdb_id, pocket_index=pocket_index,
        engine="vina", center=list(map(float, center)), box_size=list(box_size),
        poses_affinity_kcal_per_mol=affinities,
        best_affinity_kcal_per_mol=min(affinities) if affinities else float("nan"),
    )
