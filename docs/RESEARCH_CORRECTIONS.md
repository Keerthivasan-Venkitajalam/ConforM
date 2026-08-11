# Corrections to the Research Plan

1. **PDB 6ARK is not a G12D structure.** The original ensemble list included
   6ARK ("KRAS G12D-SOS1 complex") but on inspection its actual title is
   "Compound 10 covalently bound to K-RAS G12C". It was dropped from
   `configs/kras_g12d.yaml` to keep the ensemble mutant-consistent (4DST,
   5US4, 7RPZ, 5XCO retained, all confirmed G12D by their RCSB `TITLE`/`COMPND`
   records).

2. **mdpocket vs. fpocket.** The plan specifies `mdpocket`, Discngine's
   multi-frame/trajectory extension of fpocket, which requires all ensemble
   frames to be atom-identical (same trajectory). Our fallback ensemble
   consists of independently solved crystal structures with differing
   resolved residues, crystallographic ligands, and waters — they cannot be
   merged into one mdpocket trajectory file without lossy renumbering. We run
   standalone `fpocket` (same underlying Voronoi-tessellation algorithm,
   same binary family, BSD-3-Clause) independently per structure instead.
   True mdpocket ensemble-persistence tracking should be revisited once a
   real BioEmu (or MD) trajectory with consistent topology is available.

3. **TICA replaced with PCA.** Part 4 of the plan recommends TICA for
   dimensionality reduction. TICA estimates a lag-time correlation and
   therefore requires a time-ordered trajectory. Our ensemble (whether the
   BioEmu-unavailable fallback, or even true BioEmu output, which draws
   i.i.d. equilibrium samples rather than a continuous trajectory) has no
   meaningful temporal ordering, so TICA is not statistically valid here.
   PCA on Cα coordinates is used instead, consistent with the plan's own
   guidance elsewhere ("If TICA is not scientifically valid ... use PCA").

4. **GNINA CNN scoring not executed.** The plan treats GNINA as the P0
   primary docking tier. GNINA's CNN scoring backend requires a CUDA build;
   this environment has no GPU. Plain AutoDock Vina (empirical scoring
   function) was substituted as the fast-screening tier. This is explicitly
   the plan's own fallback path ("GNINA -> smina/Vina" in the failure-mode
   table), not a silent substitution.

5. **vina Python package install path.** `pip install vina` fails on macOS
   without a pre-installed Boost library (`Boost library location was not
   found!`). The working installation path used here was the conda-forge
   binary build (`conda install -c conda-forge vina`), which ships prebuilt
   and required no Boost setup.

6. **No PostgreSQL/pgvector, MCP, or REINVENT4 in this pass.** These are
   listed in the plan as P1/P3 priorities to be added after the P0
   deterministic pipeline works (master prompt rule #42). Per that explicit
   priority ordering, this pass stopped at a working, real P0 pipeline and
   did not start P1+ components; see `docs/IMPLEMENTATION_STATUS.md`.
