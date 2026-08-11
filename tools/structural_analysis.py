"""Structural ensemble analysis: RMSD, RMSF, PCA over Calpha coordinates.

Uses PCA rather than TICA. TICA assumes a meaningfully ordered time series
(a trajectory) so it can estimate a lag-time correlation; our fallback
ensemble is a pool of independent experimental crystal structures with no
temporal ordering, so TICA would be scientifically invalid here (see
docs/RESEARCH_CORRECTIONS.md, Part 4 of the research plan). PCA on the
Calpha coordinate covariance is the appropriate choice for an unordered
structural ensemble.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms


@dataclass
class EnsembleAnalysis:
    pdb_ids: list[str]
    common_resids: list[int]
    rmsd_matrix: np.ndarray          # (N, N) pairwise CA RMSD in Angstrom, vs first structure as ref for alignment
    rmsf_per_residue: np.ndarray     # (n_residues,) CA RMSF in Angstrom
    pca_coords: np.ndarray           # (N, n_components)
    pca_explained_variance_ratio: np.ndarray


def _load_chain_a_ca(pdb_path: Path) -> mda.Universe:
    u = mda.Universe(str(pdb_path))
    protein = u.select_atoms("protein and segid A or (protein and chainID A)")
    if len(protein) == 0:
        protein = u.select_atoms("protein")
    return u, protein


def analyze_ensemble(pdb_paths: list[Path]) -> EnsembleAnalysis:
    universes = []
    resid_sets = []
    for p in pdb_paths:
        u, protein = _load_chain_a_ca(p)
        ca = protein.select_atoms("name CA")
        universes.append((u, ca))
        resid_sets.append(set(ca.resids.tolist()))

    common_resids = sorted(set.intersection(*resid_sets))
    if len(common_resids) < 10:
        raise RuntimeError(f"Too few common CA residues across ensemble: {len(common_resids)}")

    coords = []
    for u, ca in universes:
        sub = ca.select_atoms(f"resid {' '.join(str(r) for r in common_resids)}")
        by_resid = {}
        for a in sub:
            # Keep first atom per resid (handles altlocs/duplicate CA records).
            by_resid.setdefault(a.resid, a)
        xyz = np.array([by_resid[r].position for r in common_resids])
        coords.append(xyz)
    coords = np.stack(coords)  # (N, n_res, 3)

    # Superpose every structure onto structure 0 using Kabsch algorithm.
    ref = coords[0] - coords[0].mean(axis=0)
    aligned = [ref]
    for i in range(1, coords.shape[0]):
        mobile = coords[i] - coords[i].mean(axis=0)
        R, _rmsd = align_kabsch(mobile, ref)
        aligned.append(mobile @ R.T)
    aligned = np.stack(aligned)  # (N, n_res, 3)

    n = aligned.shape[0]
    rmsd_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            diff = aligned[i] - aligned[j]
            rmsd_matrix[i, j] = np.sqrt((diff ** 2).sum(axis=1).mean())

    mean_struct = aligned.mean(axis=0)
    rmsf = np.sqrt(((aligned - mean_struct) ** 2).sum(axis=2).mean(axis=0))

    flat = aligned.reshape(n, -1)
    flat_centered = flat - flat.mean(axis=0)
    n_components = min(3, n - 1) if n > 1 else 1
    if n_components >= 1 and flat_centered.shape[0] > 1:
        U, S, Vt = np.linalg.svd(flat_centered, full_matrices=False)
        pca_coords = U[:, :n_components] * S[:n_components]
        total_var = (S ** 2).sum()
        explained = (S[:n_components] ** 2) / total_var if total_var > 0 else np.zeros(n_components)
    else:
        pca_coords = np.zeros((n, 1))
        explained = np.zeros(1)

    return EnsembleAnalysis(
        pdb_ids=[p.stem for p in pdb_paths],
        common_resids=common_resids,
        rmsd_matrix=rmsd_matrix,
        rmsf_per_residue=rmsf,
        pca_coords=pca_coords,
        pca_explained_variance_ratio=explained,
    )


def align_kabsch(mobile: np.ndarray, ref: np.ndarray) -> tuple[np.ndarray, float]:
    """Standard Kabsch superposition; returns rotation matrix and RMSD."""
    H = mobile.T @ ref
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    rotated = mobile @ R.T
    rmsd = np.sqrt(((rotated - ref) ** 2).sum(axis=1).mean())
    return R, rmsd


if __name__ == "__main__":
    import sys

    paths = [Path(p) for p in sys.argv[1:]]
    result = analyze_ensemble(paths)
    print("pdb_ids:", result.pdb_ids)
    print("n_common_residues:", len(result.common_resids))
    print("rmsd_matrix (A):\n", np.round(result.rmsd_matrix, 2))
    print("max RMSF (A):", round(float(result.rmsf_per_residue.max()), 2))
    print("PCA explained variance ratio:", np.round(result.pca_explained_variance_ratio, 3))
