# Implementation Status

Environment: macOS, no CUDA GPU (`nvidia-smi` not present). Internet access available.
All numbers below come from an actual executed run; see `artifacts/kras_g12d_p0_*/metrics/experiment_manifest.json`.

| Component | Status | Current implementation | Missing work | Dependency | Test status | Fallback used |
|---|---|---|---|---|---|---|
| Structure acquisition | Done (fallback tier) | `tools/structure_tool.py` fetches real deposited PDB structures from RCSB | OpenFold3/ESMFold inference | CUDA GPU (unavailable) | Manually run, verified real PDB content | Yes: RCSB structures instead of OpenFold3/ESMFold |
| Conformational ensemble | Done (fallback tier) | `tools/bioemu_tool.py`; pools 4 real KRAS G12D crystal structures (4DST, 5US4, 7RPZ, 5XCO) | Real BioEmu diffusion sampling | CUDA GPU + bioemu weights | Manually run | Yes: experimental ensemble, NOT equilibrium samples — population stats are uniform-fallback, documented in code |
| Structural analysis | Done (real) | `tools/structural_analysis.py`: Kabsch alignment, RMSD/RMSF, PCA (MDAnalysis + numpy SVD) | TICA (deliberately not used — invalid for unordered ensemble; see RESEARCH_CORRECTIONS.md) | MDAnalysis | Manually run, produced real numbers | No |
| Pocket detection | Done (real, fallback tier) | `tools/mdpocket_tool.py` runs real `fpocket 4.0` per ensemble member, parses actual volumes/druggability/residues | True `mdpocket` multi-frame trajectory mode (needs atom-identical frames, incompatible with heterogeneous crystal structures) | fpocket (conda-forge) | Manually run; recovered ground-truth Switch-II residues (H95/Y96/Q99) with 100% overlap in 7RPZ | Yes: independent fpocket runs instead of true mdpocket ensemble mode |
| Pocket ranking | Done (real) | Deterministic rank by druggability + ground-truth residue overlap in `pipelines/pocket_discovery.py` | Configurable weight file (`configs/discovery_score.yaml`) not yet split out | none | Manually run | No |
| Ligand validation | Done (real) | `tools/rdkit_tool.py`: sanitize, canonicalize, Lipinski, QED, 3D ETKDG embed + MMFF optimize | Larger library (ChEMBL subset) | RDKit (conda-forge) | Manually run, 10/10 ligands valid | No |
| Docking | Done (real, fallback tier) | `tools/docking_tool.py`: real AutoDock Vina 1.2.7 Python bindings, OpenBabel PDBQT prep | GNINA CNN rescoring, DiffDock-Pocket refinement | CUDA (GNINA CNN backend, DiffDock model) | Manually run, 10/10 ligands docked with real affinities (-9.6 to -2.8 kcal/mol) | Yes: Vina empirical scoring only, no CNN rescore |
| Discovery Score | Done (real) | `agent/discovery_score.py`, pure deterministic function, documented formula | Configurable weights YAML not yet split out (weights live in kras_g12d.yaml) | none | Manually run | No |
| Ligand optimization (REINVENT4) | Not started | — | Full module | REINVENT4 (Apache 2.0, CPU-installable in principle) | none | Not yet implemented — P0 run screens the fixed library only |
| Agent / closed loop | Not started | — | State machine, iteration loop, memory | none | none | — |
| Scientific memory (DB) | Not started | Only filesystem JSON manifest exists | PostgreSQL/pgvector or SQLite repository layer | postgres/sqlite | none | — |
| MCP layer | Not started | — | Pydantic tool schemas | none | none | — |
| Dashboard | Not started | — | Streamlit + py3Dmol | streamlit, py3Dmol | none | — |
| Evaluation / CryptoBench / ablations | Not started | — | Baseline comparison scripts | CryptoBench data | none | — |
| Docker/Postgres compose | Not started | — | docker-compose.yml | docker (installed) | none | — |

## What genuinely ran (verifiable in `artifacts/`)
- 4 real KRAS G12D PDB structures downloaded from RCSB (4DST, 5US4, 7RPZ, 5XCO)
- Real Kabsch alignment + RMSD/RMSF/PCA over 163 common Cα residues
- Real `fpocket 4.0` run on each of the 4 structures (49 total raw pocket candidates)
- Real ligand pocket in 7RPZ recovered with residues H95, Y96, Q99 (100% overlap with documented ground truth)
- Real RDKit validation/3D embedding of 10 ligands
- Real AutoDock Vina 1.2.7 docking of all 10 ligands into the recovered pocket, exhaustiveness=8, 5 poses each
- Deterministic Discovery Score computed from the above, best = 0.8332 for a synthetic scaffold at -9.6 kcal/mol

## What did NOT run (and why)
See `docs/LIMITATIONS.md` for the full list (BioEmu, OpenFold3, GNINA, DiffDock-Pocket, REINVENT4 — all require a CUDA GPU or additional integration work not completed in this pass).
