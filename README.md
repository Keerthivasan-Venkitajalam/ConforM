# ConforM-Agent (P0 slice)

A closed-loop scientific pipeline for discovering transient/cryptic druggable
pockets on KRAS G12D and screening ligand hypotheses against them.

**Status: P0 deterministic pipeline only.** The autonomous agent, ligand
optimization loop, database, MCP layer, and dashboard described in the full
research plan are **not yet built** — see `docs/IMPLEMENTATION_STATUS.md` for
the exact breakdown of what runs today vs. what remains.

## Scientific honesty note
This is a computational hypothesis-generation tool. It never claims a ligand
"binds" KRAS — only that it produced a favorable score under the evaluated
computational protocol. See `docs/LIMITATIONS.md` before interpreting any
number below.

## What actually runs today
`sequence → real RCSB structures (fallback for OpenFold3/BioEmu, no GPU
here) → RMSD/RMSF/PCA (MDAnalysis) → fpocket pocket detection → RDKit ligand
validation → AutoDock Vina docking → deterministic Discovery Score`

This recovers the documented Switch-II ground-truth pocket residues
(H95/Y96/Q99) from PDB 7RPZ with 100% overlap, and docks a small ligand set
into it with real Vina affinities. See `docs/RESEARCH_CORRECTIONS.md` for
where this deviates from the original research plan and why.

## Quick start

```bash
conda env create -f environment.yml
conda activate conform
python pipelines/pocket_discovery.py
```

Artifacts (structures, pockets, ligands, docking poses, experiment manifest,
agent log) are written to `artifacts/kras_g12d_p0_<timestamp>/`.

## Tests

```bash
conda activate conform
python -m pytest tests/ -q
```

## Repository layout

```
agent/          discovery_score.py (deterministic scoring — no LLM/agent loop yet)
tools/          structure_tool, bioemu_tool, structural_analysis, mdpocket_tool, rdkit_tool, docking_tool
pipelines/      pocket_discovery.py (P0 end-to-end runner)
configs/        kras_g12d.yaml (target, fallback ensemble PDB IDs, weights)
data/           ligands_kras.csv (small test ligand library)
docs/           IMPLEMENTATION_STATUS.md, LIMITATIONS.md, RESEARCH_CORRECTIONS.md
tests/          unit tests for discovery_score, RDKit validation, pocket parsing
artifacts/      per-experiment outputs (gitignored)
```

## GPU requirements
None of the currently implemented code requires a GPU. BioEmu, OpenFold3,
GNINA, and DiffDock-Pocket (from the full research plan) all require CUDA
and are not yet integrated — see `docs/LIMITATIONS.md`.
