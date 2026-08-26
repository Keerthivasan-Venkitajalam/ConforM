"""Computes ground-truth pocket residues objectively from a holo structure's
ligand contacts, rather than hand-typing residue numbers from literature
(which risks transcription error and cannot be independently re-derived).

Usage:
    python scripts/compute_ground_truth.py <pdb_path> <ligand_resname> [--cutoff 4.5] [--chain A]

Used to generate the ground_truth_pocket_residues lists in configs/*.yaml
for KRAS G12D (7RPZ/6IC), PRMT5 (6UXX/QL1), and ABL kinase (1IEP/STI). Every
number reported by this script is directly checkable against the cited PDB
file -- there is no hidden literature lookup step to get wrong.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E",
    "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def hetatm_coords(pdb_path: Path, resname: str) -> np.ndarray:
    coords = [
        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        for line in pdb_path.read_text().splitlines()
        if line.startswith("HETATM") and line[17:20].strip() == resname
    ]
    if not coords:
        raise ValueError(f"No HETATM records for '{resname}' found in {pdb_path}")
    return np.array(coords)


def contact_residues(pdb_path: Path, ligand_coords: np.ndarray, cutoff: float = 4.5,
                     chain: str | None = None) -> list[str]:
    residues = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if chain and line[21] != chain:
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if np.min(np.linalg.norm(ligand_coords - xyz, axis=1)) <= cutoff:
            resname = line[17:20].strip()
            resid = line[22:26].strip()
            one = THREE_TO_ONE.get(resname, "X")
            residues.add(f"{one}{resid}")
    return sorted(residues, key=lambda r: int(r[1:]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb_path", type=Path)
    ap.add_argument("ligand_resname")
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--chain", default=None)
    args = ap.parse_args()

    lig = hetatm_coords(args.pdb_path, args.ligand_resname)
    residues = contact_residues(args.pdb_path, lig, args.cutoff, args.chain)
    print(f"{len(lig)} ligand atoms, {len(residues)} contact residues within {args.cutoff} A:")
    print(residues)
