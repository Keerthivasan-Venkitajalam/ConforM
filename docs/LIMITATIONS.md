# Limitations

## Hardware
This build was developed and executed on a Mac with **no CUDA-capable GPU**
(`nvidia-smi` is not present). The research plan's MVP assumes a single
RTX 4090 (24GB VRAM). Every GPU-dependent model below could not be run here:

| Tool | Why it can't run here | What ran instead |
|---|---|---|
| BioEmu | Diffusion-model inference requires CUDA | 4 real experimental KRAS G12D crystal structures used as an ensemble (`tools/bioemu_tool.py`) |
| OpenFold3 | AlphaFold3-scale model, requires high-VRAM GPU | Structures fetched directly from RCSB (`tools/structure_tool.py`) |
| ESMFold | Runs on GPU or very slow CPU; not installed | Not exercised — RCSB fetch used directly |
| GNINA | 3D CNN scoring backend needs a CUDA-built fork of Caffe | Plain AutoDock Vina empirical scoring only |
| DiffDock-Pocket | Diffusion pose model, GPU required | Not run |
| REINVENT4 | RL training loop is GPU-accelerated in practice, and was not integrated in this pass at all | Not run; ligand library screened as a fixed set |

## Consequences for scientific interpretation
- The "ensemble" used here is **not** a Boltzmann-weighted equilibrium sample. It is
  4 independently solved crystallographic states of KRAS G12D. State "frequency"
  in the Discovery Score is therefore a uniform placeholder (`1/N`), not a real
  population estimate. Any claim about "rare state discovery" from this ensemble
  is not statistically meaningful until BioEmu (or MD) is actually run.
- Pocket detection used standalone `fpocket` per structure rather than true
  `mdpocket` ensemble/trajectory mode, because the crystal structures are not
  atom-identical frames of one trajectory (different resolved residues, ligands,
  crystallographic waters). Cross-state pocket "persistence" was not computed;
  only per-structure pocket volume/druggability was.
- Docking scores are raw AutoDock Vina empirical-function affinities
  (kcal/mol). They have not been validated against GNINA's CNN scores or any
  experimental Ki/Kd data in this pass, so treat them as an unvalidated,
  fast-screening-tier signal only, per the tiered docking design in the
  research plan.
- No claim is made that any docked ligand is a real binder. Per the master
  prompt's scientific-honesty rule, all outputs should be read as "produced a
  favorable computational docking score under the evaluated protocol," not as
  evidence of biological activity.

## Not yet built
Agent/closed-loop state machine, scientific memory database, MCP tool
schemas, REINVENT4 integration, GNINA/DiffDock-Pocket tiers, Streamlit
dashboard, CryptoBench-based benchmarking/ablations, Docker/Postgres
orchestration. See `docs/IMPLEMENTATION_STATUS.md` for the itemized table
and `docs/RESEARCH_CORRECTIONS.md` for corrections to the original plan's
specific technical claims.
