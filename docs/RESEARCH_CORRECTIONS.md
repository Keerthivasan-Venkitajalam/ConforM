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

7. **Discovery-ensemble circularity: 7RPZ removed from the discovery set
   (2026-08-25).** The original ensemble (4DST, 5US4, 7RPZ, 5XCO) included
   7RPZ, the KRAS G12D structure co-crystallized with MRTX1133 -- the very
   inhibitor that holds the Switch-II pocket open. Ranking candidate pockets
   across an ensemble that already contains the ligand-forced-open state is
   circular: it demonstrates that the pipeline can *recognize* a pocket a
   drug has already carved out, not that it can *discover* one from
   unliganded conformations. This would not withstand scrutiny from a
   structural biologist reviewer.

   Fix: every ensemble member was re-audited against its RCSB HETATM
   records. The corrected discovery ensemble contains only structures whose
   sole heteroatoms are the physiological GDP/Mg2+ cofactor plus ordinary
   crystallographic solvent (glycerol, ethylene glycol, water, backbone
   capping groups) -- no synthetic small-molecule ligand at any site:

   | PDB  | HETATM records          | Note |
   |------|--------------------------|------|
   | 5US4 | GDP, GOL, HOH, MG        | apo-like |
   | 5XCO | GDP, ACE, EDO, HOH, NH2  | apo-like (ACE/NH2 are peptide caps) |
   | 7F0W | GDP, HOH, MG             | apo-like; Switch-I open conformation |
   | 7EYX | GDP, HOH                 | apo-like; Mg-free |

   7RPZ (MRTX1133) and 4DST (ligand 9LI at a distinct, non-Switch-II
   pocket) were dropped from the discovery ensemble. 7RPZ is retained only
   as the literature citation for `ground_truth_pocket_residues` -- the
   residue list itself is hardcoded from the McCarthy et al. Switch-II SAR
   literature, and no code path opens `7RPZ.pdb` during discovery or
   ranking. See `configs/kras_g12d.yaml` for the annotated ensemble
   definition and `tests/test_no_leakage.py` for the automated check that
   ground-truth residues cannot reach a ranking function.

   As a separate, independent verification that pocket detection is not an
   artifact of the bound ligand's own volume: stripping all HETATM records
   from 7RPZ (including MRTX1133 itself) and re-running fpocket on the bare
   protein coordinates still recovers the cavity at 820 Å³ with H95/Y96/Q99
   lining residues -- confirming the geometry is defined by the protein
   backbone/side chains, not by ligand-occupied volume.
