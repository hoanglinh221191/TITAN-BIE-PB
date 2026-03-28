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
- interaction energy / force analysis
- GBSA / ALPB
- hydrogen bonds and geometry observables
- RMSD / RMSF / DSSP / LJ order parameters
- DCCM / network / PCA / allostery
- TICA / MSM / DHAM / UMAP / TPT
- LLPS / GIST / IFP / LSP-MD

The project is performance-sensitive. Most heavy loops are expected to be
streaming, cache-aware, SIMD-friendly, and conservative about allocations.

================================================================================
2. First files to read
================================================================================

If you only have a few minutes, read these in order:

1. `arguments.h`
   Central runtime configuration and default values.

2. `parse_arguments.cpp`
   CLI parsing, module implication rules, and option semantics.

3. `impl_engine.cpp`
   Main orchestration file. This is where setup, streaming, module dispatch,
   and data collection all meet.

4. `chunk_streamer.h` and `chunk_streamer.cpp`
   Input streaming layer for PDB and XTC/TRR trajectories.

5. `engine_postproc.h` and `engine_postproc.cpp`
   Post-loop analysis orchestration for TICA, MSM, DHAM, PCA, UMAP, LLPS, etc.

6. `info.h`
   Output filename source of truth.

7. `manual.txt` and `help.cpp`
   User-facing documentation. If behavior changes, these often need updates.

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
  `LlpsArgs`, `GistArgs`, `UmapArgs`, etc.).
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
- `titan_dssp.cpp`
- `titan_lj.cpp`
- `titan_pca.cpp`
- `titan_kinetics.cpp` for TICA
- `titan_msm.cpp` for MSM
- `titan_dham.cpp` for DHAM
- `titan_tpt.cpp` for TPT
- `titan_umap.cpp`, `titan_umap_knn.cpp`
- `titan_llps.cpp`, `titan_llps_cluster.cpp`
- `titan_gist.cpp`
- `titan_ifp.cpp`
- `cas_analysis.cpp`
- `water_analysis.cpp`

Output and UI:
- `info.h`: output filename getters
- `titan_ui.h` and related implementation: user-facing logging/progress

================================================================================
5. Module dependency graph
================================================================================

Important implied flags from parsing:
- `-pca` implies `-network`
- `-allostery` implies `-network`
- `-tica` implies `-network`
- `-umap` implies `-tica` and `-network`
- `-msm` implies `-tica` and `-network`
- `-dham` implies `-msm`, `-tica`, and `-network`
- `-tpt` implies `-msm`, `-tica`, and currently defaults to `-network`
- `-dssp` enables a secondary-structure feature module; do not assume it is
  only a side effect of protein workflows
- `-lspmd` is standalone and does not need `-network`
- `-rms` and `-lj` auto-enable `-pbcn` when the user did not choose an
  explicit PBC mode
- `-lj` used with `-tica` / `-msm` / `-dham` / `-umap` / `-tpt` may disable an
  auto-enabled `-network` path and use permutation-invariant LJ phase-2
  features instead

Interpretation:
- `network` owns the aligned Ca/backbone trajectory used downstream.
- `tica` consumes aligned network-derived features.
- `lj` can also act as a trace-free phase-2 feature source for kinetics-style
  workflows.
- `dssp` is its own feature pipeline and should be reviewed as such.
- `msm` consumes TICA projections.
- `dham` consumes MSM states plus temperature/energy metadata.
- `tpt` consumes the built MSM TPM plus source/sink basin definitions.
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

Files to inspect for this pipeline:
- `impl_engine.cpp`
- `engine_postproc.cpp`
- `titan_kinetics.h` and `titan_kinetics.cpp`
- `titan_msm.h` and `titan_msm.cpp`
- `titan_dham.h` and `titan_dham.cpp`
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
- Windows: `clang++` via `compile.bat` and `mingw32-make.exe`
- Linux: `g++` with `-std=c++23`

Important practical rule:
- Treat Linux support as "C++23 syntax with GCC 13-era reality".
- Do not assume full latest-standard library coverage just because `-std=c++23`
  is enabled.
- Current local build policy uses `-march=native` in `Makefile`.
  This is good for workstation tuning, but it is not a portability contract.

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
- Use `compile.bat` on Windows for project builds.
- If compilation behavior needs changing, edit `Makefile`, not ad-hoc compile
  commands.
- `release` target adds LTO flags.
- Keep AVX2/FMA assumptions paired with scalar fallback paths; do not assume
  every supported build machine is identical to the local workstation.

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
- inspect `titan_pca.cpp` for downstream weight coupling

If you change trajectory streaming:
- inspect `chunk_streamer.cpp`
- inspect XDR and PDB path handling in `impl_engine.cpp`
- verify replica boundaries still propagate correctly
- recheck whether the change increases frame-copy pressure in the split path

If you change RMS / DSSP / LJ / TPT behavior:
- inspect `titan_rms.cpp`
- inspect `titan_dssp.cpp`
- inspect `titan_lj.cpp`
- inspect `titan_tpt.cpp`
- inspect `arguments.cpp` and `parse_arguments.cpp` for implied-module and PBC
  semantics

================================================================================
13. Common pitfalls
================================================================================

- Do not break module implication logic in `parse_arguments.cpp`.
- Do not hardcode filenames when `info.h` already owns them.
- Do not add raw `std::cout` in hot engine paths if `TitanUI` is the intended
  logging interface.
- Do not accidentally disable boundary-aware lag handling in TICA/MSM.
- Do not add memory-hungry caches without considering XDR and multi-replica
  workflows.
- Do not bypass dense pair-cache guards just because a workstation build has
  enough RAM today.
- Do not assume `exact_dense` TICA is safe for workflows that need projection
  replay.
- Do not treat LJ phase-2 kinetics, DSSP features, RMS traces, and network-
  derived TICA features as interchangeable; the runtime chooses among them
  based on flags.
- Do not replace jthread-based scoped parallelism with unmanaged threads.
- Do not introduce shiny C++23 library features that GCC 13 likely does not
  support.

================================================================================
14. Quick mental model
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
