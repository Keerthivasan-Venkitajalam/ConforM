"""Conformational ensemble generation.

Fallback hierarchy (docs/IMPLEMENTATION_STATUS.md, docs/LIMITATIONS.md):
    BioEmu -> experimental ensemble -> generated conformer set

BioEmu requires a CUDA GPU for diffusion-model inference. This machine has
none (`nvidia-smi` unavailable), so BioEmuProvider.generate() always raises
and the pipeline falls back to ExperimentalEnsembleProvider, which uses only
real, experimentally deposited KRAS G12D structures (multiple independent
crystal forms/ligand states) downloaded from RCSB via structure_tool.py.

IMPORTANT SCIENTIFIC CAVEAT: unlike BioEmu samples, these snapshots are NOT
i.i.d. draws from a Boltzmann equilibrium distribution, so per-state
population/frequency has no rigorous thermodynamic meaning here. The
pipeline still computes a "frequency" field for interface compatibility with
the Discovery Score, but it is defined as 1/N_states (uniform) and is
labeled `population_is_uniform_fallback=True` everywhere it is used so
downstream consumers do not mistake it for equilibrium statistics.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from tools.structure_tool import StructureProvider


@dataclass
class EnsembleResult:
    provider: str
    structures: list[Path]
    is_equilibrium_sample: bool
    metadata: dict = field(default_factory=dict)


class BioEmuProvider:
    name = "bioemu"

    def generate(self, *_args, **_kwargs) -> EnsembleResult:
        raise NotImplementedError(
            "BioEmu requires a CUDA-capable GPU for diffusion inference; none "
            "detected in this environment. See docs/LIMITATIONS.md."
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
    if config.get("ensemble", {}).get("bioemu", {}).get("enabled"):
        try:
            return BioEmuProvider().generate()
        except NotImplementedError as exc:
            print(f"[bioemu_tool] BioEmu unavailable, falling back: {exc}")
    pdb_ids = config["target"]["ensemble_pdb_ids"]
    return ExperimentalEnsembleProvider().generate(pdb_ids, out_dir)


if __name__ == "__main__":
    import yaml

    cfg = yaml.safe_load(Path("configs/kras_g12d.yaml").read_text())
    res = get_ensemble(cfg, Path("data/structures"))
    print(f"provider={res.provider} n_states={len(res.structures)} equilibrium={res.is_equilibrium_sample}")
