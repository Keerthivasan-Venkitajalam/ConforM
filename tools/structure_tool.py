"""Structure acquisition with an explicit provider fallback hierarchy.

Fallback order (per docs/IMPLEMENTATION_STATUS.md):
    OpenFold3 -> ESMFold -> ExistingStructureProvider (RCSB PDB fetch)

OpenFold3 and ESMFold require GPU inference and multi-GB model weights that
are not available in this environment (see docs/LIMITATIONS.md). Only
ExistingStructureProvider is real/executed here. It downloads actual
experimentally-deposited coordinates from RCSB -- it never fabricates atoms.
"""
from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


@dataclass
class StructureResult:
    pdb_id: str
    path: Path
    provider: str


class ExistingStructureProvider:
    """Fetches a real, experimentally deposited structure from RCSB PDB."""

    name = "rcsb_existing_structure"

    def fetch(self, pdb_id: str, out_dir: Path) -> StructureResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        dest = out_dir / f"{pdb_id.upper()}.pdb"
        if not dest.exists():
            url = RCSB_URL.format(pdb_id=pdb_id.upper())
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = resp.read()
            if not data.startswith(b"HEADER") and b"ATOM" not in data[:2000]:
                raise RuntimeError(f"RCSB fetch for {pdb_id} did not return a valid PDB file")
            dest.write_bytes(data)
        return StructureResult(pdb_id=pdb_id.upper(), path=dest, provider=self.name)


class OpenFold3Provider:
    """Not runnable in this environment: needs CUDA GPU + multi-GB weights."""

    name = "openfold3"

    def fetch(self, *_args, **_kwargs):
        raise NotImplementedError(
            "OpenFold3 requires a CUDA-capable GPU (RTX 4090+) and model weights "
            "not available in this environment. See docs/LIMITATIONS.md."
        )


class ESMFoldProvider:
    """Not installed in this environment (no GPU); documented as unavailable."""

    name = "esmfold"

    def fetch(self, *_args, **_kwargs):
        raise NotImplementedError(
            "ESMFold inference was not installed/executed in this environment. "
            "See docs/LIMITATIONS.md."
        )


class StructureProvider:
    """Tries providers in fallback order, records which one actually ran."""

    def __init__(self):
        self._chain = [OpenFold3Provider(), ESMFoldProvider(), ExistingStructureProvider()]

    def get_baseline(self, pdb_id: str, out_dir: Path) -> StructureResult:
        errors = []
        for provider in self._chain:
            try:
                return provider.fetch(pdb_id, out_dir)
            except NotImplementedError as exc:
                errors.append(f"{provider.name}: {exc}")
                continue
        raise RuntimeError("All structure providers failed:\n" + "\n".join(errors))


if __name__ == "__main__":
    import sys

    sp = StructureProvider()
    result = sp.get_baseline(sys.argv[1] if len(sys.argv) > 1 else "4DST", Path("data/structures"))
    print(f"provider={result.provider} path={result.path}")
