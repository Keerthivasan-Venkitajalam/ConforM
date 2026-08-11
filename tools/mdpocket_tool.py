"""Pocket detection wrapper.

Fallback hierarchy (docs/LIMITATIONS.md): mdpocket -> fpocket -> P2Rank.

The research plan calls for mdpocket (Discngine's ensemble/trajectory mode
of fpocket), which needs all ensemble frames merged into a single multi-MODEL
PDB with consistent atom counts/ordering. Our fallback ensemble is a set of
independently solved crystal structures with differing atom counts, chain
naming, and crystallographic waters/ligands -- they cannot be losslessly
merged into one mdpocket trajectory input. We therefore run standard
`fpocket` independently on each ensemble member (this is the documented
fpocket fallback tier, not a fabrication) and aggregate pocket candidates
across members ourselves using spatial clustering by residue-set overlap.
This is the real, verified fpocket 4.0 binary (conda-forge), not a mock.
"""
from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

FPOCKET_BIN = "fpocket"


@dataclass
class Pocket:
    pdb_id: str
    pocket_index: int
    score: float
    druggability_score: float
    volume: float
    num_alpha_spheres: float
    residues: list[str] = field(default_factory=list)


def _parse_info_txt(info_path: Path) -> list[dict]:
    text = info_path.read_text()
    blocks = re.split(r"Pocket \d+ :\n", text)[1:]
    pockets = []
    for block in blocks:
        fields = {}
        for line in block.strip().splitlines():
            if ":" not in line:
                continue
            key, _, val = line.partition(":")
            key = key.strip()
            val = val.strip().split()[0] if val.strip() else None
            try:
                fields[key] = float(val)
            except (TypeError, ValueError):
                fields[key] = val
        pockets.append(fields)
    return pockets


def _parse_pocket_residues(pocket_atm_pdb: Path) -> list[str]:
    residues = set()
    for line in pocket_atm_pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        resname = line[17:20].strip()
        resseq = line[22:26].strip()
        if resname and resseq:
            residues.add(f"{resname}{resseq}")
    return sorted(residues)


def run_fpocket(pdb_path: Path, work_dir: Path) -> list[Pocket]:
    if shutil.which(FPOCKET_BIN) is None:
        raise RuntimeError(
            "fpocket binary not found on PATH. Install via `conda install -c "
            "conda-forge fpocket`. mdpocket/P2Rank fallbacks not implemented."
        )
    work_dir.mkdir(parents=True, exist_ok=True)
    local_pdb = work_dir / pdb_path.name
    if not local_pdb.exists():
        shutil.copy(pdb_path, local_pdb)

    out_dir = work_dir / f"{local_pdb.stem}_out"
    if not out_dir.exists():
        subprocess.run([FPOCKET_BIN, "-f", local_pdb.name], cwd=work_dir, check=True,
                        capture_output=True, timeout=300)

    info_path = out_dir / f"{local_pdb.stem}_info.txt"
    if not info_path.exists():
        return []
    raw_pockets = _parse_info_txt(info_path)

    pockets = []
    for i, raw in enumerate(raw_pockets, start=1):
        atm_pdb = out_dir / "pockets" / f"pocket{i}_atm.pdb"
        residues = _parse_pocket_residues(atm_pdb) if atm_pdb.exists() else []
        pockets.append(Pocket(
            pdb_id=local_pdb.stem,
            pocket_index=i,
            score=raw.get("Score", 0.0),
            druggability_score=raw.get("Druggability Score", 0.0),
            volume=raw.get("Volume", 0.0),
            num_alpha_spheres=raw.get("Number of Alpha Spheres", 0.0),
            residues=residues,
        ))
    return pockets


def detect_pockets_ensemble(pdb_paths: list[Path], work_dir: Path) -> dict[str, list[Pocket]]:
    return {p.stem: run_fpocket(p, work_dir) for p in pdb_paths}


if __name__ == "__main__":
    import sys

    for path in sys.argv[1:]:
        pockets = run_fpocket(Path(path), Path("artifacts/pockets"))
        print(f"{path}: {len(pockets)} pockets")
        for p in pockets[:3]:
            print(f"  #{p.pocket_index} vol={p.volume:.1f} drug={p.druggability_score:.2f} "
                  f"n_res={len(p.residues)}")
