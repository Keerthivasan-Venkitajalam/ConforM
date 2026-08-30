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

8. **First real GPU/BioEmu run initially reported a false `n_states=1`
   result for every ablation mode -- retracted and re-verified
   (2026-08-29).** A GPU session on this project's own RTX 4060 (WSL2,
   `nvidia-smi` confirmed live) produced an ablation table where every
   mode, including `conform-agent`, showed `n_states=1` and `conform-agent`
   scored ground-truth recall 0.80 -- suspicious for a `num_samples=100`
   BioEmu request. That number was never published anywhere and was caught
   before use. Root cause, fully traced:

   - **`tools/bioemu_tool.py`, `BioEmuProvider.generate()`.** Real BioEmu
     inference genuinely ran on real CUDA (34 batch `.npz` files plus a
     `samples.xtc` trajectory were on disk, ~100 real GPU-generated
     conformers). But the code only collected structures via
     `out_dir.glob("*.pdb")` / `frame_*.pdb` / `samples_*.pdb` -- none of
     which match BioEmu's `samples.xtc` + `topology.pdb` trajectory output
     format. Only `topology.pdb` (the single reference frame) matched, so
     the entire real ensemble was silently discarded and every downstream
     step (pocket detection, ranking, docking) ran against 1 structure
     while believing it had sampled 100. Fixed by extracting per-frame
     PDBs from the trajectory via MDAnalysis (`_extract_frames`), covered
     by `tests/test_bioemu_frame_extraction.py`, which fails if this
     regresses.
   - **`agent/loop_controller.py`, `build_manifest()`/`provenance()`.**
     `"cuda": "unavailable (no nvidia-smi on this host)"` and the entire
     `tools_unavailable_fallback_used` dict were hardcoded string literals
     written into every manifest unconditionally, left over from before
     GPU support existed. They directly contradicted the real `.npz`/`.xtc`
     evidence sitting in the same experiment directory. Fixed: `cuda` now
     reflects a live `cuda_available()` check, and
     `tools_unavailable_fallback_used` only reports a tool as "fallback
     used" if this specific run's own state shows it actually fell back
     (`BioEmu` only appears if `ensemble_provider != "bioemu"`; `REINVENT4`
     only if `optimizer_mode != "reinvent4"`). `GNINA`/`OpenFold3`/
     `DiffDock-Pocket` were removed from this per-run dict entirely --
     `loop_controller.py`'s closed loop never invokes any of the three
     (confirmed against its own `tools_executed` list), so claiming a
     fallback for tools never attempted was misleading regardless of GPU
     status.
   - **`evaluation/baselines.py`, `static` mode.** A third, independent bug
     surfaced once the first two were fixed and the ablation table was
     re-run: the `static` baseline (meant to be "single apo structure, no
     conformational sampling at all" -- the naive control) only overrode
     `target.ensemble_pdb_ids`, which feeds the *fallback* ensemble
     provider. It never set `ensemble.bioemu.enabled = False`, and
     `get_ensemble()` tries BioEmu first whenever that flag is true
     (inherited unchanged from the main config). So with real BioEmu
     enabled, `static` silently sampled the same full ~100-state ensemble
     as every other mode and then ran through the exact same
     `rank_pocket_families()` ranking as `conform-agent` -- defeating the
     entire purpose of a no-sampling control. This is why the first
     corrected table still showed `static` tied with the real pipeline at
     0.80 recall. Fixed by explicitly forcing
     `ensemble.bioemu.enabled = False` in `static`'s sub-config.

   **Verified real result** (`kras_g12d_conform-agent_1787945235`,
   git commit at run time `753753f`, re-run after both fixes): real BioEmu
   inference, **99 real states**, `is_equilibrium_sample=True`, max RMSF
   5.00 Å, max pairwise RMSD 4.89 Å, PC1 explained variance 0.199 (all
   0.0 in the false run). Ground-truth pocket residue recall **0.60**,
   best Discovery Score 0.7194. GNINA CNN rescoring, run for the first time
   against a real, verified pocket/ligand set, produced real CNNscore/
   CNNaffinity values (top hit `Imatinib_fragment_F_0`: CNNscore 0.683,
   CNNaffinity 6.654) -- see `docs/IMPLEMENTATION_STATUS.md`.

   **Corrected 5-mode ablation table** (all real BioEmu runs, one flagship
   `static` re-run after its fix, four other modes' original real runs
   reused unchanged since that bug did not affect them):

   | mode | states | pockets | recall | best kcal | discovery | iters | sec |
   |---|---|---|---|---|---|---|---|
   | static | 1 | 19 | 0.00 | -9.81 | 0.383 | 1 | 48 |
   | random | 99 | 1192 | 0.00 | -6.85 | 0.595 | 1 | 526 |
   | no-pocket-guidance | 98 | 1210 | 0.00 | -7.90 | 0.628 | 1 | 519 |
   | no-ligand-optimization | 100 | 1225 | 0.80 | -9.41 | 0.781 | 1 | 531 |
   | conform-agent | 94 | 1185 | 0.60 | -9.07 | 0.756 | 2 | 767 |

   **Honest caveat this run surfaced, reported rather than hidden:**
   `no-ligand-optimization` and `conform-agent` run the *identical* pocket
   selection algorithm (`rank_pocket_families`), yet scored 0.80 vs. 0.60.
   Because each mode draws its **own independent** BioEmu sample rather
   than sharing one ensemble, that 0.20 gap is more likely sampling
   variance between two separate real diffusion draws than a real effect
   of the optimization step -- pocket selection happens in iteration 0,
   before the optimization branch is even reached, so optimization cannot
   causally change which pocket gets picked. Weak supporting evidence: an
   earlier standalone `conform-agent` run
   (`kras_g12d_conform-agent_1787945235`, the one GNINA-rescored above)
   independently landed on the same 0.60 recall -- two of two
   `conform-agent` draws at 0.60 vs. one `no-ligand-optimization` draw at
   0.80 -- but n=2 is far too small to separate a real effect from noise.
   The clean, *within-run* evidence for the optimization step is instead
   `conform-agent`'s own before/after: best library affinity -8.75 kcal/mol
   -> best analog -9.07 kcal/mol (delta -0.32 kcal/mol, 12 analogs
   generated and re-docked). A statistically meaningful cross-mode
   ablation would need multiple BioEmu seeds per mode; that was not
   feasible in the remaining compute/time budget and is recorded here as a
   known limitation rather than glossed over -- see
   `docs/GENERALIZATION.md` for how this interacts with the earlier
   CPU-fallback three-target generalization result, which used a
   different, non-stochastic ensemble source and is not affected by this
   issue.

9. **The #8 caveat resolved with a shared-ensemble re-run and a permutation
   test (2026-08-30).** #8 flagged, but did not fix, the independent-sampling
   confound in the ablation table. This entry fixes it, using the already
   GPU-verified 99-state ensemble from #8's `conform-agent` run
   (`kras_g12d_conform-agent_1787945235`) as a shared master ensemble --
   **no new BioEmu sampling was performed**; only the deterministic
   pocket-detection/ranking/docking steps were re-run per mode
   (`evaluation/shared_ensemble_ablation.py`).

   | mode | states | pockets | recall | best kcal | discovery | iters | sec |
   |---|---|---|---|---|---|---|---|
   | static (true single-structure control) | 1 | 19 | 0.00 | -9.81 | 0.379 | 1 | 49 |
   | random | 99 | 1228 | 0.60 | -7.24 | 0.667 | 1 | 101 |
   | no-pocket-guidance | 99 | 1228 | 0.00 | -6.98 | 0.598 | 1 | 82 |
   | no-ligand-optimization | 99 | 1228 | 0.60 | -8.29 | 0.718 | 1 | 86 |
   | conform-agent | 99 | 1228 | 0.60 | -8.80 | 0.719 | 2 | 310 |

   With every non-static mode now evaluating the **identical** 99-state
   ensemble, `no-ligand-optimization` and `conform-agent` converge on the
   exact same selected pocket family (0.60 recall both) -- confirming #8's
   caveat was correct: the earlier 0.80-vs-0.60 gap was sampling variance
   between two independent BioEmu draws, not a real effect of the
   optimization step. The clean within-ensemble evidence for the
   optimization step is unchanged from #8: -8.29 -> -8.80 kcal/mol
   (Δ -0.51 kcal/mol) on the identical pocket, with no confound left to
   explain it away.

   **A second real bug surfaced running this**: `no-pocket-guidance`
   crashed with `StopIteration`. Its selection logic
   (`evaluation/baselines.py`) assumed the ensemble's first structure
   always has at least one detected pocket family -- untrue in general, and
   false for `frame_0000` of this specific real ensemble. Fixed to walk the
   ensemble in order and take the first structure that actually has a
   detection, still applying no volumetric/novelty guidance.

   **An honest, unresolved-by-construction nuance**: `random` (seed=42)
   also landed on 0.60 recall this time -- a different outcome from the
   original independent-sampling table, where `random` scored 0.00. This
   is not a contradiction; it is a real property of *this* ensemble: enough
   of its 149 detected pocket families apparently have decent ground-truth
   overlap that one random draw can land well. A single random seed cannot
   distinguish "the ranking algorithm is better than chance" from "this
   ensemble is easy" -- that is exactly what the permutation test below is
   for, and a single `random` draw is not treated as evidence on its own.

   **Permutation test** (`evaluation/permutation_test.py`, unit-tested in
   `tests/test_permutation_test.py`): holding the real, unmodified
   `rank_score` ordering of all 149 candidate pocket families fixed,
   10,000 random relabelings of which family carries the ground-truth
   overlap were generated. The real ranking placed the best-overlap family
   available in the ensemble (overlap 1.0 -- notably not the family the
   algorithm actually selected, which had overlap 0.60 and ranked #1) at
   rank 2 of 149. A rank that good occurred in **1.44% of random
   relabelings (empirical p = 0.0144, seed=42)**. No p-value was
   pre-specified; this is what the test returned. This is evidence that
   `rank_score` is associated with ground-truth overlap beyond chance on
   this specific ensemble -- it is not, by itself, evidence that the method
   generalizes across ensembles or targets (see `docs/GENERALIZATION.md`
   for that separate question).

   **Selection-margin logging added** (`agent/loop_controller.py`): every
   `SELECT_POCKET` decision now also logs the score gap between the chosen
   pocket family and its runner-up (e.g. "0.998 vs. runner-up 0.978, margin
   +0.020" for this run) -- a small, low-confidence-visible signal about
   how decisive a given selection was, logged from data the ranking already
   computes. This is not a claim of probabilistic reasoning; it is exactly
   what it says: a margin between two already-computed scores.
