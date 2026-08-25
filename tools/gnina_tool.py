"""GNINA CNN-scoring docking wrapper (CUDA required).

GNINA (gnina/gnina, GPL-2.0/Apache-2.0) forks AutoDock Vina and replaces/
augments its empirical scoring function with an ensemble of 3D CNNs. It
needs a CUDA-capable GPU; this project's local dev machine has none, so
this module has never been executed end-to-end here (see
docs/LIMITATIONS.md). It is written against GNINA's documented CLI
(https://github.com/gnina/gnina, v1.3.x) and is intended to run unmodified
the first time on rented GPU compute via `scripts/gpu_session.sh`.

Install (no build required as of v1.3.2 -- a prebuilt binary is published):
    wget https://github.com/gnina/gnina/releases/download/v1.3.2/gnina -O gnina
    chmod +x gnina
    # requires CUDA >= 12.0 on the host

CLI reference used here (verified against the official README, NOT guessed):
    gnina -r receptor.pdbqt -l ligand.pdbqt \
        --center_x X --center_y Y --center_z Z \
        --size_x SX --size_y SY --size_z SZ \
        --cnn_scoring rescore --exhaustiveness N --num_modes N \
        -o out.sdf --log log.txt
GNINA writes CNNscore/CNNaffinity/CNNvariance and minimizedAffinity (Vina-
style empirical score) as SDF per-pose properties on the output file.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

GNINA_BIN_ENV_DEFAULT = "gnina"


@dataclass
class GninaPose:
    cnn_score: float           # CNN pose-quality probability, [0, 1], higher = better
    cnn_affinity: float        # CNN predicted binding affinity (pK units, higher = better)
    cnn_variance: float | None
    vina_affinity_kcal: float | None  # minimizedAffinity, Vina-style empirical kcal/mol


@dataclass
class GninaResult:
    ligand_name: str
    receptor_pdb_id: str
    pocket_index: int
    engine: str = "gnina"
    cnn_scoring_mode: str = "rescore"
    poses: list[GninaPose] = field(default_factory=list)
    best_cnn_score: float = float("nan")
    best_cnn_affinity: float = float("nan")
    raw_output_path: str | None = None


def gnina_available(binary: str = GNINA_BIN_ENV_DEFAULT) -> bool:
    return shutil.which(binary) is not None


def _parse_sdf_poses(sdf_path: Path) -> list[GninaPose]:
    """Parse GNINA's per-pose SDF output properties.

    Uses RDKit's SDMolSupplier (already a project dependency) rather than
    hand-rolled SDF parsing. Property names (CNNscore, CNNaffinity,
    CNNvariance, minimizedAffinity) are GNINA's documented output tags.
    """
    from rdkit import Chem

    poses = []
    supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=False)
    for mol in supplier:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        poses.append(GninaPose(
            cnn_score=float(props.get("CNNscore", "nan")),
            cnn_affinity=float(props.get("CNNaffinity", "nan")),
            cnn_variance=float(props["CNNvariance"]) if "CNNvariance" in props else None,
            vina_affinity_kcal=float(props["minimizedAffinity"])
                if "minimizedAffinity" in props else None,
        ))
    return poses


def rescore(receptor_pdbqt: Path, ligand_pdbqt: Path, center: np.ndarray,
           box_size=(24.0, 24.0, 24.0), exhaustiveness: int = 8, num_modes: int = 5,
           cnn_scoring: str = "rescore", binary: str = GNINA_BIN_ENV_DEFAULT,
           work_dir: Path = Path("artifacts/gnina"), ligand_name: str = "ligand",
           receptor_pdb_id: str = "receptor", pocket_index: int = 0,
           timeout: int = 600) -> GninaResult:
    """Runs real GNINA docking/CNN-rescoring. Raises if the binary is absent
    or the subprocess fails -- never fabricates scores on failure."""
    if not gnina_available(binary):
        raise RuntimeError(
            f"GNINA binary '{binary}' not found on PATH. Install the prebuilt "
            "release (see module docstring) or pass an explicit binary path. "
            "This tool requires a CUDA-capable GPU."
        )
    work_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = work_dir / f"{ligand_name}_gnina.sdf"
    log_path = work_dir / f"{ligand_name}_gnina.log"

    cmd = [
        binary,
        "-r", str(receptor_pdbqt),
        "-l", str(ligand_pdbqt),
        "--center_x", str(float(center[0])),
        "--center_y", str(float(center[1])),
        "--center_z", str(float(center[2])),
        "--size_x", str(float(box_size[0])),
        "--size_y", str(float(box_size[1])),
        "--size_z", str(float(box_size[2])),
        "--cnn_scoring", cnn_scoring,
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", "42",
        "-o", str(out_sdf),
        "--log", str(log_path),
    ]
    subprocess.run(cmd, check=True, capture_output=True, timeout=timeout)

    if not out_sdf.exists():
        raise RuntimeError(f"GNINA reported success but produced no output at {out_sdf}")

    poses = _parse_sdf_poses(out_sdf)
    if not poses:
        raise RuntimeError(f"GNINA output {out_sdf} contained no parseable poses")

    best = max(poses, key=lambda p: p.cnn_score)
    return GninaResult(
        ligand_name=ligand_name, receptor_pdb_id=receptor_pdb_id, pocket_index=pocket_index,
        cnn_scoring_mode=cnn_scoring, poses=poses,
        best_cnn_score=best.cnn_score, best_cnn_affinity=best.cnn_affinity,
        raw_output_path=str(out_sdf),
    )


if __name__ == "__main__":
    import sys

    print(f"gnina binary available: {gnina_available()}")
    if len(sys.argv) >= 3 and gnina_available():
        # Smoke test: gnina --autobox_ligand mode against a known pocket-adjacent
        # ligand, using pocket_centroid from docking_tool for the box center.
        from tools.docking_tool import pocket_centroid
        receptor, ligand, pqr = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
        center = pocket_centroid(pqr)
        result = rescore(receptor, ligand, center, ligand_name="smoke_test")
        print(f"best_cnn_score={result.best_cnn_score:.3f} "
              f"best_cnn_affinity={result.best_cnn_affinity:.3f} "
              f"n_poses={len(result.poses)}")
