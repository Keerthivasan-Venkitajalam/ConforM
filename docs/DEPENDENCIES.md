# Dependencies and Licenses

## Installed and actually executed in this build
| Package | Version | License | Install path that worked |
|---|---|---|---|
| RDKit | conda-forge | BSD-3-Clause | `conda install -c conda-forge rdkit` |
| fpocket | 4.0 | MIT / BSD-3-Clause | `conda install -c conda-forge fpocket` |
| AutoDock Vina | 1.2.7 | Apache-2.0 | `conda install -c conda-forge vina` |
| Open Babel | 3.1.0 | GPL-2.0 | `conda install -c conda-forge openbabel` |
| MDAnalysis | conda-forge | GPL-2.0+ | `conda install -c conda-forge mdanalysis` |
| NumPy / SciPy / pandas / scikit-learn | conda-forge | BSD-3-Clause | conda |
| Streamlit | 1.61.1 | Apache-2.0 | pip/conda |
| py3Dmol | latest | BSD-3-Clause | pip/conda |
| pytest | latest | MIT | pip/conda |

**Install note:** `pip install vina` fails on macOS without a system Boost
installation (`ValueError: Boost library location was not found!`). The
conda-forge binary build is the working path and is what `environment.yml`
uses.

**GPL note:** Open Babel (GPL-2.0) and MDAnalysis (GPL-2.0+) are invoked as
libraries/subprocesses. If this project is distributed as a combined work,
GPL obligations apply to the distribution. GNINA is also GPL-2.0. This is
compatible with an open-source competition submission but should be stated
in the submission's license section.

## Declared in the plan but NOT installed or executed
| Package | License | Why not | Consequence |
|---|---|---|---|
| BioEmu | MIT / Apache-2.0 | CUDA GPU required; none on host | Experimental-structure ensemble fallback |
| OpenFold3 | Apache-2.0 | High-VRAM GPU + weights | RCSB structures used |
| ESMFold | MIT | GPU (or impractically slow CPU) | Not exercised |
| GNINA | GPL-2.0 / Apache-2.0 | CUDA-built Caffe fork required | Vina empirical scoring only |
| DiffDock-Pocket | MIT | GPU diffusion model | Not run |
| REINVENT 4 | Apache-2.0 | Prior checkpoints + GPU RL loop; not integrated | RDKit R-group enumeration fallback |
| CryptoBench | MIT | Full benchmark needs ~1,107 apo-holo pairs | KRAS-specific ground truth used instead |
| PostgreSQL / pgvector | PostgreSQL License | Backend declared, code path not exercised | SQLite default |

## Dependency isolation strategy
The GPU tools have mutually conflicting CUDA/PyTorch pins and must **not**
share one environment. The intended layout when a GPU host is available:

```
core env (this environment.yml)   ← RDKit, fpocket, Vina, MDAnalysis, app
  ├── bioemu worker  (own env, own CUDA/PyTorch pin)
  ├── gnina          (binary or container)
  └── reinvent4      (own env)
```

Each is reached through the provider interfaces already in `tools/`
(`StructureProvider`, `bioemu_tool.get_ensemble`, `docking_tool.dock`,
`reinvent_tool.get_optimizer`), so enabling a GPU tool means implementing one
provider class — no pipeline changes.
