"""Conformational ensemble generation.

Fallback hierarchy (docs/IMPLEMENTATION_STATUS.md, docs/LIMITATIONS.md):
    BioEmu -> experimental ensemble -> generated conformer set

BioEmu (microsoft/bioemu) requires a CUDA GPU for diffusion-model inference.
BioEmuProvider auto-detects CUDA via torch; when present it runs REAL
`bioemu.sample.main()` inference (pip install bioemu[cuda]) and returns
actual generated conformers -- this is genuine generative sampling from the
apo sequence alone, not a stand-in. When no GPU is available it raises and
the pipeline falls back to ExperimentalEnsembleProvider, which uses only
real, experimentally deposited KRAS G12D structures (multiple independent
apo-like crystal forms) downloaded from RCSB via structure_tool.py.

IMPORTANT SCIENTIFIC CAVEAT for the fallback path: unlike true BioEmu
samples, these crystal-structure snapshots are NOT i.i.d. draws from a
Boltzmann equilibrium distribution, so per-state population/frequency has no
rigorous thermodynamic meaning. The pipeline still computes a "frequency"
field for interface compatibility with the Discovery Score, but it is
defined as 1/N_states (uniform) and is labeled
`population_is_uniform_fallback=True` everywhere it is used so downstream
consumers do not mistake it for equilibrium statistics. True BioEmu output
sets `is_equilibrium_sample=True` and real per-state populations are
computable from the returned samples.
"""
from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from tools.structure_tool import StructureProvider


@dataclass
class EnsembleResult:
    provider: str
    structures: list[Path]
    is_equilibrium_sample: bool
    metadata: dict = field(default_factory=dict)


def cuda_available() -> bool:
    """Detect a usable CUDA GPU via nvidia-smi.

    Deliberately does NOT fall back to `import torch; torch.cuda.is_available()`.
    On some local dev machines (observed on this project's macOS/conda setup)
    importing torch alongside rdkit triggers a double-initialized-OpenMP
    runtime abort (SIGABRT / exit 134) at the C level, which Python's
    exception handling cannot catch -- it would take down the whole process,
    not just this check. nvidia-smi is the standard, dependency-free signal
    of real GPU hardware and is reliably present on any CUDA-enabled rental
    instance (RunPod/Vast.ai/Lambda all ship it), so it is sufficient on its
    own without needing torch imported this early.
    """
    try:
        subprocess.run(["nvidia-smi"], capture_output=True, timeout=10, check=True)
        return True
    except Exception:  # noqa: BLE001
        return False


class BioEmuProvider:
    """Real microsoft/bioemu diffusion-model inference (requires CUDA).

    Install: pip install bioemu[cuda]   (Python 3.10+; ~3.5GB AlphaFold2
    weights auto-download to ~/.cache/colabfold/ on first use)
    """

    name = "bioemu"

    def generate(self, sequence: str, out_dir: Path, num_samples: int = 1000,
                 model_name: str = "bioemu-v1.1", filter_samples: bool = True,
                 seed: int | None = None) -> EnsembleResult:
        if not cuda_available():
            raise NotImplementedError(
                "BioEmu requires a CUDA-capable GPU for diffusion inference; none "
                "detected in this environment (torch.cuda.is_available() is False "
                "and nvidia-smi is unavailable). See docs/LIMITATIONS.md."
            )
        try:
            from bioemu.sample import main as bioemu_sample
        except ImportError as exc:
            raise NotImplementedError(
                "bioemu package not installed. Run: pip install bioemu[cuda]"
            ) from exc

        out_dir = Path(out_dir) / "bioemu_run"
        out_dir.mkdir(parents=True, exist_ok=True)

        if seed is not None:
            import random
            import numpy as np
            random.seed(seed)
            np.random.seed(seed)
            try:
                import torch
                torch.manual_seed(seed)
            except ImportError:
                pass

        bioemu_sample(
            sequence=sequence,
            num_samples=num_samples,
            output_dir=str(out_dir),
            filter_samples=filter_samples,
            model_name=model_name,
        )

        # bioemu.sample writes samples.xtc + topology.pdb (multi-frame trajectory)
        # and/or per-frame PDBs depending on version; collect whatever it produced.
        pdb_files = sorted(out_dir.glob("*.pdb"))
        frame_files = sorted(out_dir.glob("frame_*.pdb")) or sorted(out_dir.glob("samples_*.pdb"))
        structures = frame_files or pdb_files
        if not structures:
            raise RuntimeError(
                f"BioEmu ran but produced no PDB output in {out_dir}; check "
                f"samples.xtc/topology.pdb for trajectory-format output requiring "
                f"MDAnalysis extraction to per-frame PDBs."
            )

        return EnsembleResult(
            provider=self.name,
            structures=structures,
            is_equilibrium_sample=True,
            metadata={
                "num_states": len(structures),
                "population_is_uniform_fallback": False,
                "model_name": model_name,
                "requested_samples": num_samples,
                "filter_samples": filter_samples,
                "seed": seed,
                "note": (
                    "Real BioEmu diffusion-model samples: an approximated equilibrium "
                    "Boltzmann distribution generated from the apo sequence alone, "
                    "with no experimental structure (including ligand-bound ones) "
                    "used as input."
                ),
            },
        )


class ExperimentalEnsembleProvider:
    """Real fallback: downloads and pools multiple deposited PDB conformers."""

    name = "fallback_experimental_ensemble"

    def __init__(self):
        self._structure_provider = StructureProvider()

    def generate(self, pdb_ids: list[str], out_dir: Path) -> EnsembleResult:
        paths = []
        for pdb_id in pdb_ids:
            result = self._structure_provider.get_baseline(pdb_id, out_dir)
            paths.append(result.path)
        return EnsembleResult(
            provider=self.name,
            structures=paths,
            is_equilibrium_sample=False,
            metadata={
                "num_states": len(paths),
                "population_is_uniform_fallback": True,
                "note": (
                    "NOT BioEmu output. These are independent experimental crystal "
                    "structures, not equilibrium samples from a generative model."
                ),
            },
        )


def get_ensemble(config: dict, out_dir: Path) -> EnsembleResult:
    bioemu_cfg = config.get("ensemble", {}).get("bioemu", {})
    if bioemu_cfg.get("enabled"):
        try:
            sequence = config["target"]["sequence"].replace("\n", "").strip()
            return BioEmuProvider().generate(
                sequence=sequence,
                out_dir=out_dir,
                num_samples=bioemu_cfg.get("num_samples", 1000),
                model_name=bioemu_cfg.get("model_name", "bioemu-v1.1"),
                filter_samples=bioemu_cfg.get("filter_samples", True),
                seed=bioemu_cfg.get("seed"),
            )
        except NotImplementedError as exc:
            print(f"[bioemu_tool] BioEmu unavailable, falling back: {exc}")
    pdb_ids = config["target"]["ensemble_pdb_ids"]
    return ExperimentalEnsembleProvider().generate(pdb_ids, out_dir)


if __name__ == "__main__":
    import yaml

    cfg = yaml.safe_load(Path("configs/kras_g12d.yaml").read_text())
    res = get_ensemble(cfg, Path("data/structures"))
    print(f"provider={res.provider} n_states={len(res.structures)} equilibrium={res.is_equilibrium_sample}")
