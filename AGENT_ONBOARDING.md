# TITAN Codebase Onboarding for AI Agents and New Contributors

This file is the fast-path map for understanding TITAN's code structure,
runtime workflow, and performance constraints before making changes.

Read this first together with:
- `AGENT_RULES.md` for hard editing/build rules.
- `TITAN_RULES.md` for architectural constraints inside the engine.

================================================================================
1. What TITAN is
================================================================================

TITAN is a high-performance molecular trajectory analysis engine centered on:
- interaction energy / force analysis / subgroup decomposition (`-sg`)
- GBSA / ALPB / screening-style reweighting helpers
- hydrogen bonds, HBond/salt/water kinetics, and geometry observables
- RMSD / RMSF / DSSP / Rg2 / LJ condensation / structure factor `S(q)`
- DCCM / network / communities / PCA / allostery / corr-kNN / LSP-MD
- TICA / MSM / Validation / DHAM / MBAR / UMAP / TPT / ZTS
- LLPS / GIST / IFP / AWH RC-binned IFP

The project is performance-sensitive. Most heavy loops are expected to be
streaming, cache-aware, SIMD-friendly, and conservative about allocations.

================================================================================
2. First files to read
================================================================================

If you only have a few minutes, read these in order:

1. `arguments.h`
   Central runtime configuration and grouped defaults (`MsmArgs`, `DhamArgs`,
   `MbarArgs`, `ValidationArgs`, `CorrKnnArgs`, `SqArgs`,
   `HBondKineticsArgs`, `LlpsArgs`, `AwhArgs`, PBC knobs, subgroup `-sg`
   knobs, etc.).

2. `parse_arguments.cpp`
   CLI parsing, module implication rules, option normalization, and validation
   glue for new modules such as `-corr-knn`, `-validation`, `-mbar`, `-sq`,
   `-hbk`, and `-sg`.

3. `impl_engine.cpp`
   Main orchestration file. This is where setup, streaming analyzers (`HBK`,
   `SQ`, `Rg2`, `IFP`, `AWH`, subgroup decomposition), and postproc launch all
   meet.

4. `chunk_streamer.h` and `chunk_streamer.cpp`
   Input streaming layer for PDB and XTC/TRR trajectories.

5. `engine_postproc.h` and `engine_postproc.cpp`
   Post-loop solve/export orchestration for TICA, MSM, Validation, DHAM, MBAR,
   TPT, UMAP, LLPS, corr-kNN, and LSP-MD.

6. `dccm_network.cpp`
   DCCM cache, allostery, community export, aligned-trace ownership, and the
   LSP-MD implementation used downstream.

7. `titan_kinetics.*`, `titan_msm.*`, `titan_dham.*`, `titan_validation.*`
   The core kinetics / reweighting / uncertainty pipeline.

8. `titan_corr_knn.*`, `titan_sq.*`, `titan_hbond_kinetics.*`,
   `titan_llps.*`, `titan_llps_cluster.*`
   Exact continuous MI/GC/TE, RDF/S(q), streaming HBond/salt/water
   kinetics, and the upgraded LLPS stack (macro-contact core, DBSCAN/fuzzy
   morphology, lineage, dense-phase mobility, and RNA/topology layers).

9. `titan_ifp.*`, `titan_awh_ifp.*`, `titan_interaction_primitives.*`,
   `titan_selection_utils.*`, `titan_sg_streaming.*`,
   `titan_energy_decomposition.*`
   Interaction fingerprinting, AWH binning, shared interaction helpers, and
   the subgroup cache-vs-stream decomposition path.

10. `manual.txt`, `help.cpp`, `info.h`, `compile.bat`, `Makefile`,
    `solv_pbc_pbctopo.cpp`
    User docs, output filename source of truth, real build flow, and the main
    topology-aware PBC hotspot.

================================================================================
3. High-level runtime workflow
================================================================================

The overall execution model is:

1. Parse CLI / arg file into `arguments`.
2. Resolve defaults, validate flags, and apply implied modules.
3. Initialize topology, groups, force-field data, and reference frame(s).
4. Open a `ChunkStreamer` for either:
   - text PDB/GRO-like chunks, or
   - XDR streaming (`.xtc` / `.trr`) with a structure reference from `-s`.
5. Stream frames chunk-by-chunk through the main loop in `impl_engine.cpp`.
6. Accumulate online or cached data for active modules.
7. After the trajectory loop finishes, run post-processing in
   `engine_postproc.cpp`.
8. Export results using `info.h` filename getters and log through `TitanUI`.

Important pattern:
- TITAN tries to do cheap, streaming work inside the trajectory loop.
- Heavier solves and exports happen after the loop.

================================================================================
4. Core architecture map
================================================================================

CLI and configuration:
- `arguments.h`: defaults and grouped config structs (`MsmArgs`, `DhamArgs`,
  `MbarArgs`, `ValidationArgs`, `CorrKnnArgs`, `SqArgs`,
  `HBondKineticsArgs`, `LlpsArgs`, `GistArgs`, `UmapArgs`, `AwhArgs`, subgroup
  `-sg` settings, etc.).
- `parse_arguments.cpp`: parsing and module dependency activation.
- `arguments.cpp`: validation and post-processing of options.

Input, topology, and setup:
- `initialization.cpp`: startup path and initial frame handling.
- `pre_proc_pdb.cpp`: frame-0 preprocessing and topology-aware setup.
- `take_topo.cpp`, `take_para.cpp`, `take_charge.cpp`, `take_coordinates.cpp`:
  topology, parameters, charges, and coordinate loading.

Streaming:
- `chunk_streamer.cpp`: chunked PDB/XDR reads, frame splitting into receptor /
  ligand / water, and jthread-based coordinate partitioning.

Main engine:
- `impl_engine.cpp`: central runtime coordinator.
  This file decides module initialization, memory strategy, chunk size, and
  per-frame work routing.

Post-processing:
- `engine_postproc.cpp`: heavy solves and exports after accumulation.

Analysis modules:
- `energy.cpp`, `energy_intrin.cpp`, `force.cpp`, `force_intrin.cpp`
- `dccm_network.cpp`
- `titan_rms.cpp`
- `titan_rg.cpp`
- `titan_dssp.cpp`
- `titan_lj.cpp`
- `titan_sq.cpp`, `titan_sq_core.cpp`
- `titan_hbond_kinetics.cpp`
- `titan_pca.cpp`
- `titan_kinetics.cpp` for TICA
- `titan_msm.cpp` for MSM
- `titan_dham.cpp` for DHAM
- `titan_validation.cpp`, `titan_validation_core.cpp`
- `titan_tpt.cpp` for TPT
- `titan_umap.cpp`, `titan_umap_knn.cpp`
- `titan_corr_knn.cpp`, `titan_corr_knn_core.cpp`, `titan_corr_knn_shared.cpp`
- `titan_llps.cpp`, `titan_llps_contacts.cpp`, `titan_llps_async.cpp`,
  `titan_llps_dynamics.cpp`, `titan_llps_export.cpp`
- `titan_llps_cluster.cpp`, `titan_llps_cluster_dbscan.cpp`,
  `titan_llps_cluster_dynamics.cpp`,
  `titan_llps_cluster_rna_compaction.cpp`,
  `titan_llps_cluster_rna_information.cpp`,
  `titan_llps_cluster_rna_topology.cpp`,
  `titan_llps_cluster_export.cpp`, `titan_llps_stats.*`
- `titan_gist.cpp`
- `titan_ifp.cpp`
- `titan_awh_ifp.cpp`
- `titan_sg_streaming.cpp`, `titan_energy_decomposition.cpp`
- `titan_interaction_primitives.cpp`, `titan_selection_utils.cpp`
- `cas_analysis.cpp`
- `water_analysis.cpp`

LLPS split of responsibility:
- `titan_llps.cpp`: macro-node facade, config, frame orchestration, C(t)
  finalize, and per-frame accumulation flow.
- `titan_llps_contacts.cpp`: macro contacts, protein-RNA contact search,
  pair packing, PBC distance helpers, GCC Rg, and LLPS kernel dispatch.
- `titan_llps_async.cpp`: bounded async frame queue and worker lifecycle.
- `titan_llps_dynamics.cpp`: contact lifetimes, rewiring, interface turnover,
  and detailed contact traces.
- `titan_llps_export.cpp`: macro summary, density/performance/memory reports,
  sticker-hotspot export, output audit, and log.txt LLPS diagnostics.
- `titan_llps_cluster*.cpp`: bead-level DBSCAN/fuzzy morphology, lineage,
  droplet-shape timeseries, genealogy/FRAP proxies, dense-phase MSD/liquidity,
  RNA compaction/information/topology, and DBSCAN/RNA export families.
- `titan_llps_stats.*`: shared LLPS percentile, fitting, entropy, binning,
  CSV-formatting, and small packing helpers.

Shared math / reusable numerics:
- `titan_gyration.cpp`: reusable `Rg`, `Rg^2`, principal gyration axes, shape
  descriptors, and AVX2 SoA point-stat accumulation
- `titan_validation_core.cpp`: CK/VAMP/bootstrap math shared by Validation and
  TPT uncertainty flows
- `titan_corr_knn_shared.cpp`: pair filtering / config-sharing logic reused by
  corr-kNN and allostery kNN mode
- `titan_sq_core.cpp`: RDF/S(q) kernel helpers and isotropic Fourier path
- `pbctopo_hull.cpp`: reusable convex-hull volume and hull-rescue support for
  pbctopo-style reconstruction paths

Output and UI:
- `info.h`: output filename getters
- `titan_ui.h` and related implementation: user-facing logging/progress

PBC and topology-aware reconstruction:
- `solv_pbc_pbctopo.cpp`: VIBE-pbctopo (`-pbctopo`, Virtual Intact Block
  Extraction), the topology-aware PBC path used for intact-block recovery on
  difficult chain-aware systems. This module keeps its detailed timing/debug
  output in `pbctopo_debug.log`, separate from the main `log.txt`.

================================================================================
5. Module dependency graph
================================================================================

Important implied flags from parsing:
- `-pca` implies `-network`
- `-allostery` implies `-network`
- `-corr-knn` implies `-network`
- `-tica` implies `-network`
- `-umap` implies `-tica` and `-network`
- `-msm` implies `-tica` and `-network`
- `-dham` implies `-msm`, `-tica`, and `-network`
- `-validation` implies `-msm`, `-tica`, and `-network`
- `-mbar` is standalone REMD/HREMD reweighting; `domain=tica` implies
  `-msm/-tica/-network`, `domain=pca` implies `-pca/-network`, and
  `domain=both` enables both paths
- `-tpt` implies `-msm`, `-tica`, and currently defaults to `-network`
- `-dssp` enables a secondary-structure feature module; do not assume it is
  only a side effect of protein workflows
- `-lspmd` is standalone and does not need `-network`
- `-rms`, `-lj`, `-rg2`, and `-sq` auto-enable a safe PBC path when the user
  did not choose an explicit PBC mode
- `-hbk`, `-ifp`, `-awh-ifp`, and `-sg` stay in the main streaming loop and do
  not depend on the DCCM/TICA stack
- `-lj` used with `-tica` / `-msm` / `-dham` / `-umap` / `-tpt` may disable an
  auto-enabled `-network` path and use permutation-invariant LJ phase-2
  features instead

Interpretation:
- `network` owns the aligned Ca/backbone trajectory used downstream by
  allostery, corr-kNN, TICA, and optional equilibrium overlays.
- `corr-knn` reuses the aligned DCCM trace cache and can also back
  `-allostery-estimator=knn`.
- `tica` consumes aligned network-derived features.
- `lj` can also act as a trace-free phase-2 feature source for kinetics-style
  workflows.
- `dssp` is its own feature pipeline and should be reviewed as such.
- `msm` consumes TICA projections.
- `dham` consumes MSM states plus temperature/energy metadata.
- `validation` reuses solved TICA/MSM/TPT context; it should not trigger a
  second raw-trajectory collection path.
- `mbar` is a standalone reweighting solver whose frame weights only become
  default downstream weights when `-eq-weight-source mbar` resolves.
- `tpt` consumes the built MSM TPM plus source/sink basin definitions.
- `sq`, `hbk`, `ifp`, `awh-ifp`, and `sg` are streaming analyzers with their
  own outputs and caches.
- `llps` is a hybrid path: per-frame macro contact accumulation starts in
  `TitanLLPS`, while richer morphology/dynamics/RNA exports complete through
  `TitanLLPSCluster` finalize/export work after streaming.
  Detailed chain/RNA traces are streamed when possible; DBSCAN shape output is
  written as `llps_droplet_shape_timeseries.csv`.
- `pca` can consume DHAM frame weights if DHAM converges.

================================================================================
6. TICA -> MSM -> DHAM workflow
================================================================================

This is one of the most important flows in the repo.

Stage A: DCCM / alignment cache
- `dccm_network.cpp` provides the aligned residue-level trajectory.
- TICA uses pairwise inter-Ca distances as features, not raw xyz, to remove
  rigid-body artifacts.

Stage B: TICA accumulation
- `impl_engine.cpp` initializes `TitanKinetics`.
- The chosen feature dimension is `N_ca * (N_ca - 1) / 2`.
- Lag buffers are reset at trajectory boundaries so lagged pairs are never
  constructed across replicas/chunks.

Stage C: TICA backend selection
- Implemented in `impl_engine.cpp` + `titan_kinetics.cpp`.
- Available backends:
  - `exact_dense`: dense covariance solve for smaller feature spaces
  - `exact_dual`: exact reduced solve with in-memory feature cache
  - `randomized`: file-backed reduced solve for very large systems
- `auto` chooses based on estimated memory.
- If per-frame projections are needed by MSM/DHAM/UMAP, exact dense replay is
  not enough, so TITAN promotes to a cache-preserving backend when necessary.

Stage D: TICA solve and projection replay
- `engine_postproc.cpp::run_tica_msm_postproc()` calls `solve_tica()`.
- TITAN writes:
  - `tica_eigenvalues.dat`
  - `tica_eigenvectors.dat`
  - `tica_projections.csv`
- Projection replay comes from:
  - MSM's pending feature cache when available, or
  - DCCM cached trajectory as a fallback rebuild path.

Stage E: MSM microstates
- `TitanMSM` clusters TICA projections with K-means++.
- It builds a discrete trajectory of microstate IDs.
- Transition counting is boundary-aware by default.

Stage F: MSM model build
- `build_markov_model()` constructs:
  - count matrix
  - TPM
  - stationary distribution
- `export_results()` can additionally perform:
  - implied-timescale scan
  - crisp PCCA macrostate lumping

Stage G: DHAM
- `TitanDHAM` consumes:
  - MSM state trajectory
  - replica IDs / timestamps
  - temperatures
  - GROMACS energy `.xvg` files
- If DHAM converges, TITAN emits reweighted populations and passes frame
  weights forward into PCA/FES and TICA-FES.
- If DHAM fails, the rest of the workflow still runs in uniform-weight mode.

Stage H: MBAR
- Implemented in `engine_postproc.cpp`.
- MBAR v1 is standalone temperature-ladder reweighting and writes `mbar_*`
  outputs.
- `domain=tica|pca|both` decides which solved metadata path(s) it reuses.
- MBAR only changes downstream default equilibrium weighting when
  `-eq-weight-source mbar` is selected.

Stage I: Validation / TPT overlays
- `TitanValidation::RuntimeCache` is allocated in `impl_engine.cpp` and filled
  during postproc.
- Validation can preselect a better MSM model (`-validation-auto-pick`) before
  the final build, then run CK/VAMP/bootstrap after the model is solved.
- TPT bootstrap uncertainty can be driven either from legacy TPT flags or from
  the shared validation cache when the requested bootstrap settings match.

Files to inspect for this pipeline:
- `impl_engine.cpp`
- `engine_postproc.cpp`
- `titan_kinetics.h` and `titan_kinetics.cpp`
- `titan_msm.h` and `titan_msm.cpp`
- `titan_dham.h` and `titan_dham.cpp`
- `titan_validation.h`, `titan_validation.cpp`, and `titan_validation_core.cpp`
- `titan_tpt.cpp`
- `titan_pca.cpp`

================================================================================
7. Streaming and memory model
================================================================================

TITAN is not a naive load-everything-into-RAM program.

Key rules:
- PDB inputs may run in-memory if small enough.
- XTC/TRR inputs are streamed conservatively.
- Chunk size is chosen adaptively from estimated frame size and available RAM.
- Some modules keep only the latest state.
- Other modules intentionally cache compact derived data rather than full atom
  trajectories.

Examples:
- `ChunkStreamer` handles PDB and XDR chunk reads.
- DCCM may store the trajectory required by downstream network/PCA/TICA work.
- TICA uses either dense matrices, in-memory reduced caches, or on-disk spool
  files depending on backend.
- Dense pair caches are guarded in `impl_engine.cpp`; do not assume O(NxM)
  interaction caches are always safe or always enabled.

If you add a new module:
- Prefer streaming accumulation over full-frame storage.
- Reuse existing caches if they already contain the needed reduced data.
- Avoid a second full-trajectory copy unless there is a strong reason.

================================================================================
8. Build and toolchain constraints
================================================================================

Current build intent:
- Windows: `compile.bat` -> MSYS2 UCRT64 `mingw32-make.exe` -> `clang++.exe`
- Linux: either `g++` (GCC 13-era target) or `TOOLCHAIN=clang` with LLVM 22 +
  libc++ as configured in `Makefile`

Important practical rule:
- Treat Linux support as "C++23 syntax with GCC 13-era reality".
- Do not assume full latest-standard library coverage just because `-std=c++23`
  is enabled.
- `PROFILE=native` uses `-march=native`; `PROFILE=portable` currently uses
  generic AVX2/FMA tuning (`-mavx2 -mfma -mtune=generic`) in `Makefile`.
- The local native profile is good for workstation tuning, but it is not a
  portability contract.

Allowed / encouraged:
- `std::string::contains`
- `std::expected`
- `std::jthread`
- modern ranges when support is reliable in the target compiler setup

Avoid or treat as unsupported:
- `<mdspan>`
- C++20 modules
- `<print>`

Build workflow:
- Use `compile.bat [native|portable] [all|release|clean|print-config]` on
  Windows for project builds.
- If compilation behavior needs changing, edit `Makefile`, not ad-hoc compile
  commands.
- The batch wrapper intentionally calls `mingw32-make ... -j1` today to cap
  peak RAM. Do not assume higher `-j` is safe just because the repo is using
  incremental object builds.
- `compile.bat all` cleans the selected profile first, then snapshots the built
  binary into `test/bin/<profile>/`.
- `release` target adds LTO flags.
- Keep AVX2/FMA assumptions paired with scalar fallback paths; do not assume
  every supported build machine is identical to the local workstation.
- For Linux `TOOLCHAIN=clang`, prefer the `Makefile`-managed LLVM 22/libc++
  path rather than hand-written compiler/linker commands.
- The codebase now carries small compatibility shims because libc++ / GCC 13
  feature coverage still differs in practice (`titan_thread_compat.h`,
  `titan_charconv_compat.h`, etc.).
- On Windows, `titan.exe` may occasionally be file-locked by an active run.
  In that case, linking an alternate target name is acceptable for validation,
  but keep `Makefile` as the source of truth.
- `test/` is the local regression workspace. Keep new generated outputs in
  dedicated scenario subdirectories there instead of spraying files into the
  repo root.
- `test/bin/native` and `test/bin/portable` are the intended snapshot locations
  for profile-specific binaries used by local regression scripts.

When changing headers or build rules:
- Recheck `Makefile`.
- Be aware that stale object reuse can still bite in real workflows even with
  dependency files enabled.

================================================================================
9. Known hotspots and current architectural debt
================================================================================

These are current-behavior notes, not blanket invitations to rewrite
everything. Touch them carefully and profile before expanding scope.

- `chunk_streamer.cpp`
  Current behavior: the split-frame path copies `init_pro.m_mol[0]` and
  `init_lig.m_mol[0]` into each frame, then overwrites coordinates.
  Recommended direction: avoid additional deep copies on streaming paths unless
  a topology-preservation fix truly requires them.
- `impl_engine.cpp` + `cache_sigma_epsilon.cpp`
  Current behavior: dense pair caches are O(NxM) in memory and only enabled
  behind explicit guards.
  Recommended direction: preserve the guard logic and prefer sparse or
  streaming accumulation when adding new pairwise analyses.
- `dccm_network.cpp`
  Current behavior: several heavy paths still use `std::vector<std::vector<double>>`
  and dense pair lists.
  Recommended direction: flatten hot buffers before introducing even more pair-
  dense work.
- `titan_pca.cpp`
  Current behavior: PCA still builds a dense `3N x 3N` covariance matrix.
  Recommended direction: treat memory growth as a first-class constraint for
  large proteins and weighted workflows.
- `titan_tpt.cpp`
  Current behavior: committor solves and flux matrices are dense O(K^2)/O(K^3)
  in the number of MSM states.
  Recommended direction: do not casually increase state counts or add extra
  dense post-processing without checking solver cost and conditioning.
- `solv_pbc_pbctopo.cpp`
  Current behavior: VIBE-pbctopo is a layered heuristic solver with projection
  cuts, component-graph refinement, and anchor-based topology metrics. The
  current hotspot is usually the coarse/exact anchor scoring path inside local
  refine, not the CLI glue around `-pbctopo`. Recent AVX2 work already covers
  the gyration accumulation path, the cut-candidate edge scan, continuity
  checks, and the anchor MST relax/select loops; remaining wins are now more
  likely in component-search state management than in the already-vectorized
  geometry kernels.
  Recommended direction: profile with `pbctopo_debug.log`, preserve output
  equivalence, and treat anchor sampling / prune thresholds as benchmark-driven
  changes rather than cosmetic rewrites.

================================================================================
10. HPC coding rules
================================================================================

This repo should be treated as HPC C++, not generic app C++.

Data layout:
- Prefer Structure-of-Arrays and flat contiguous buffers.
- Avoid `std::vector<std::vector<T>>` in hot paths.
- Reuse `AlignedDoubleVector`, `AlignedIntVector`, and row-major Eigen matrices
  where that matches the existing design.

Alignment:
- 32-byte alignment is the baseline SIMD expectation.
- See `aligned_allocator.h`.
- On Windows the code uses `_aligned_malloc`.
- On non-Windows paths use `std::aligned_alloc` with padded size rounded to a
  multiple of the alignment.

SIMD:
- AVX2 is a first-class target in heavy numeric code.
- Preserve `#ifdef __AVX__` / `#if defined(__AVX2__)` dispatch and always keep a
  scalar fallback path.
- Prefer contiguous loads/stores and branch-light kernels.
- `simd_compat.h` provides FMA-compatible wrappers.

Threads:
- Prefer `std::jthread` for scoped worker lifetimes.
- TITAN already uses jthreads in streamer, DCCM, MSM, PCA, LLPS, DHAM, GIST,
  and other modules.
- The preferred style is:
  - static chunk partitioning
  - no detached threads
  - RAII join on scope exit
  - minimal allocation inside worker loops

Logging:
- In engine code, prefer `TitanUI` and log files over ad-hoc console prints.

================================================================================
11. HPC review checklist
================================================================================

When reviewing or extending a compute-heavy module, quickly verify:

- SoA / flat row-major storage is preserved in hot loops.
- 32-byte alignment assumptions still hold for AVX2 paths.
- No avoidable `std::vector<std::vector<T>>` or deep topology/frame copies were
  introduced into hot paths.
- `std::jthread` usage is amortized over meaningful work; do not respawn tiny
  worker batches inside already-small loops.
- AVX2 kernels still have scalar remainder handling and a non-AVX fallback.
- Numerical guards exist where needed for `NaN` / `Inf`, zero row sums, zero
  distances, and ill-conditioned dense solves.

================================================================================
12. Safe modification checklist
================================================================================

If you change CLI or defaults:
- update `arguments.h`
- update `parse_arguments.cpp`
- update validation in `arguments.cpp` if needed
- update `manual.txt`
- update the embedded help text in `help.cpp`

If you change output filenames:
- update `info.h`
- search for docs/examples that mention the old names

If you change TICA/MSM/DHAM behavior:
- inspect `impl_engine.cpp`
- inspect `engine_postproc.cpp`
- inspect `titan_kinetics.cpp`
- inspect `titan_msm.cpp`
- inspect `titan_dham.cpp`
- inspect `titan_validation.cpp` / `titan_validation_core.cpp` when the change
  affects model selection, CK, bootstrap, or TPT uncertainty
- inspect `titan_pca.cpp` for downstream weight coupling

If you change corr-kNN / allostery-kNN / MBAR behavior:
- inspect `engine_postproc.cpp`
- inspect `titan_corr_knn.cpp`, `titan_corr_knn_core.cpp`, and
  `titan_corr_knn_shared.cpp`
- inspect `dccm_network.cpp` because corr-kNN reuses the aligned trace cache
- inspect `titan_validation.cpp` if the change affects equilibrium overlays or
  shared bootstrap/model-selection context
- inspect `manual.txt`, `help.cpp`, and `info.h` when pair-score or summary
  outputs change

If you change trajectory streaming:
- inspect `chunk_streamer.cpp`
- inspect XDR and PDB path handling in `impl_engine.cpp`
- verify replica boundaries still propagate correctly
- recheck whether the change increases frame-copy pressure in the split path

If you change RMS / DSSP / LJ / TPT behavior:
- inspect `titan_rms.cpp`
- inspect `titan_rg.cpp` if the change affects trace-atom geometry observables
- inspect `titan_dssp.cpp`
- inspect `titan_lj.cpp`
- inspect `titan_sq.cpp` / `titan_sq_core.cpp` if the change touches trace
  atom or condensed-phase channel semantics
- inspect `titan_tpt.cpp`
- inspect `arguments.cpp` and `parse_arguments.cpp` for implied-module and PBC
  semantics

If you change IFP / HBK / AWH / subgroup decomposition behavior:
- inspect `titan_ifp.cpp`
- inspect `titan_hbond_kinetics.cpp`
- inspect `titan_awh_ifp.cpp`
- inspect `titan_interaction_primitives.cpp` and `titan_selection_utils.cpp`
- inspect `titan_sg_streaming.cpp` and `titan_energy_decomposition.cpp`
- inspect `info.h`, `parse_arguments.cpp`, `manual.txt`, and `help.cpp` if
  filenames, flags, or representative exports move

If you change LLPS behavior:
- inspect the relevant `titan_llps*.cpp` or `titan_llps_cluster*.cpp` split
  file, not only the facade files
- inspect `arguments.h` and `parse_arguments.cpp` because LLPS option coupling
  now spans density policy, kernel policy, RNA analysis, and energy correlation
- inspect `info.h`, `manual.txt`, and `help.cpp` if any LLPS output family,
  filename, or flag wording changes
- keep LLPS run metadata in `log.txt` through `TitanUI`; do not add a manifest
  sidecar unless the user explicitly asks for one
- benchmark both `-llps-kernel=scalar` and `-llps-kernel=avx2` sensitive paths
  when touching contact kernels, density accumulation, or gyration reuse
- recheck no-RNA, `-nocal`, and protein-only fallback paths when touching
  protein-RNA, compaction, topology, or energy-correlation logic

If you change gyration / shape-metric behavior:
- inspect `titan_gyration.cpp`
- inspect `titan_rg.cpp`
- inspect `titan_llps.cpp` for GCC/Rg reuse
- inspect `solv_pbc_pbctopo.cpp` because pbctopo scoring also consumes these
  metrics
- inspect `info.h`, `help.cpp`, and `manual.txt` if output names or units move

If you change VIBE-pbctopo behavior:
- inspect `solv_pbc_pbctopo.cpp`
- inspect `parse_arguments.cpp`, `arguments.cpp`, `manual.txt`, and `help.cpp`
  for `-pbctopo` semantics and user-facing wording
- keep `pbctopo_debug.log` as the detailed debug/timing sink; do not spill hot
  path diagnostics back into `log.txt`
- benchmark on a representative trajectory chunk and compare both timing and
  output hashes before keeping heuristic changes

================================================================================
13. Common pitfalls
================================================================================

- Do not break module implication logic in `parse_arguments.cpp`.
- Do not hardcode filenames when `info.h` already owns them.
- Do not let `manual.txt` and the embedded raw-string manual in `help.cpp`
  drift apart.
- Do not add raw `std::cout` in hot engine paths if `TitanUI` is the intended
  logging interface.
- Do not accidentally disable boundary-aware lag handling in TICA/MSM.
- Do not add memory-hungry caches without considering XDR and multi-replica
  workflows.
- Do not assume `compile.bat` is conservative by accident; the current `-j1`
  wrapper is there to control memory pressure.
- Do not bypass dense pair-cache guards just because a workstation build has
  enough RAM today.
- Do not move detailed VIBE-pbctopo diagnostics back into `log.txt`; use
  `pbctopo_debug.log` for that path.
- Do not assume `exact_dense` TICA is safe for workflows that need projection
  replay.
- Do not assume MBAR mutates DHAM/TICA/PCA/TPT default behavior unless
  `-eq-weight-source mbar` is actually active.
- Do not recompute large alignment/feature caches inside Validation when the
  `TitanValidation::RuntimeCache` is already the intended reuse path.
- Do not treat `-sg` as "just an export layer"; it now has real cache-vs-stream
  memory policy and kernel-selection logic.
- Do not treat LJ phase-2 kinetics, DSSP features, RMS traces, and network-
  derived TICA features as interchangeable; the runtime chooses among them
  based on flags.
- Do not replace jthread-based scoped parallelism with unmanaged threads.
- Do not introduce shiny C++23 library features that GCC 13 likely does not
  support.

================================================================================
14. Current 2026-04 Status
================================================================================

This section is intentionally current-state and should be kept in sync with
recent engine work.

Current user-facing stack additions:
- Exact corr-kNN is live in `titan_corr_knn.*` and can be reused by
  `-allostery-estimator=knn` without recollecting the aligned DCCM trace cache.
- Validation is now a first-class postproc layer (`titan_validation.*`) with
  VAMP-2 sweeps, CK checks, optional bootstrap, and optional auto-pick reuse of
  the winning MSM model.
- MBAR v1 is standalone reweighting in `engine_postproc.cpp`; it writes
  `mbar_*` outputs and only changes downstream default equilibrium overlays when
  `-eq-weight-source mbar` is selected.
- Streaming analyzers now include `-hbk`, `-sq`, `-rg2`, `-ifp`, `-awh`
  RC-binned IFP, and subgroup decomposition (`-sg`). Most are initialized in
  `impl_engine.cpp` and export via `info.h`.
- `-sg` now routes through `SgStreamingAnalyzer` plus
  `EnergyDecompositionProvider` and chooses cache/stream internals from memory
  heuristics. Treat changes there as benchmark-sensitive.
- LLPS is no longer just a coarse phase detector. The current stack adds
  `-llps-density`, `-llps-kernel`, DBSCAN/fuzzy droplet morphology, lineage and
  fusion/fission tracking, dense-phase MSD/liquidity, protein-RNA contact
  turnover, sticker hotspot ranking, droplet-shape/genealogy summaries,
  FRAP-like recovery proxies, detailed performance reporting, and `log.txt`
  run-context/output audit sections instead of a separate manifest file.
- RNA-aware LLPS analysis now has three higher-level layers behind
  `titan_llps_cluster.*`: compaction/binding coupling, information-style state
  metrics, and topology outputs (bridge/clamp/loop closure/local competition),
  with optional detailed trace exports.
- `test/` doubles as the local regression matrix plus binary snapshot area.
  Scenario directories are intentionally numerous; keep new outputs isolated in
  dedicated subdirectories.

Reusable gyration / `-rg2`:
- `titan_gyration.*` is now the shared source of truth for:
  - `Rg^2`
  - `Rg`
  - principal gyration-axis eigenvalues
  - `kappa2`, asphericity, and acylindricity
- `titan_rg.*` is a standalone analyzer, not folded into RMS code.
- `-rg2` exports:
  - `rg2*.dat`
  - `rg2_summary*.dat`
  - `rg2_hist*.dat`
- For multi-input runs it also exports aggregate files:
  - `rg2_all.dat`
  - `rg2_summary_all.dat`
  - `rg2_hist_all.dat`
- Aggregate timeseries columns are:
  - `traj_index frame time_ps rg2_A2 rg_A`
- Aggregate summary contains one row per trajectory plus an `ALL` row.

Current AVX2 coverage:
- `titan_gyration.cpp`
  - AVX2 SoA accumulation is implemented for hot `PointStats` reduction.
  - The 3x3 eigensolve and shape-descriptor formulas remain scalar by design.
- `titan_llps.cpp`
  - GCC `Rg` now reuses the same SoA gyration path instead of scalar
    point-by-point accumulation.
- `solv_pbc_pbctopo.cpp`
  - AVX2 currently covers:
    - cut-candidate edge scan
    - continuity metric scan
    - anchor MST relax step
    - anchor MST next-node selection
  - Component-search hot loops were also tightened with cheaper sorting and
    less state copying, but they are not "fully SIMD-ized".

Current pbctopo behavior:
- VIBE-pbctopo now relies on shared gyration metrics from `titan_gyration.*`.
- Split-contacting-owner rejection and rescue hierarchy are already in place.
- QuickHull was moved into reusable `pbctopo_hull.*`.
- Current regression expectation on the recent local validation windows:
  - dimer window: no final split, no rescue, no pbci fallback
  - tetramer window: one split candidate may be rejected internally, but final
    output still has no split and no pbci fallback

Representative local-only validation fixtures:
- Difficult dimer case:
  - `D:\Data_MD\Npro\dimerunphos`
- Extended tetramer case:
  - `D:\Data_MD\Npro\tetramer_conf0_unphos`
- These are workstation-local paths, not portable repo fixtures. Treat them as
  regression references when they exist, not as guaranteed CI assets.

================================================================================
15. Quick mental model
================================================================================

If you need one sentence per layer:

- `arguments.*`: what the user asked for
- `impl_engine.cpp`: what TITAN will actually run
- `chunk_streamer.*`: how frames arrive
- module files: what gets accumulated per frame
- `engine_postproc.*`: what gets solved after streaming
- `info.h`: where outputs are written
- `manual.txt` and `help.cpp`: what users think the program does

Keep those layers in sync.
