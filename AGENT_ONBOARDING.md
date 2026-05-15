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
- `titan_pbe.h` and `titan_pbe/` for frame-based PBE/PBSA electrostatics
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
  `titan_llps_cluster_export.cpp`, `titan_llps_paper.cpp`,
  `titan_llps_stats.*`
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
- `titan_llps_paper.cpp`: opt-in `-llps-paper` Martini protein/RNA paper
  bundle, using DBSCAN dense labels for phase partitioning, radial/slab
  profiles, composition time series, convergence summaries, and log audit.
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
  written as `llps_droplet_shape_timeseries.csv`. The `-llps-paper` path is
  opt-in and should keep default `-llps` output unchanged; it reuses DBSCAN
  dense labels instead of rescanning contacts globally.
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
- PBE is intentionally a frame analyzer, not a standalone one-shot tool. Keep
  `pbe_energy.csv` long-table output and memory preflight behavior intact when
  extending nonlinear, Gaussian, or focusing variants.
- PBE public API is now a small PIMPL shim in `titan_pbe.h`; implementation
  lives under `titan_pbe/`. Keep `impl_engine.cpp` including `titan_pbe.h` and
  put new PBE internals in the subdirectory. `Makefile` already builds
  `titan_pbe/*.cpp`; `compile.bat` should not need separate source-list edits.
  Phase 3 splitting has begun: map/cube I/O is in `audit_io.cpp`, PBSA
  aggregation/writers are in `pbsa.cpp`, focus/probe helpers are in
  `focus_probe.cpp`, and Gaussian/RPB audit helpers plus source diagnostics are
  in `gaussian_rpb.cpp`. Grid/radius/surface setup is now in
  `grid_surface.cpp`, and charge assignment/operator setup plus oracle replay
  helpers are now in `charge_operator.cpp`. Linear boundary/SOR/focus routing
  is now in `linear_solver.cpp`, and linear V-cycle hierarchy is now in
  `multigrid.cpp`. Continue moving one small subsystem at a time, then run
  native build plus focused PBE/RPB/PBSA smokes before the next slice.
- PBE Phase 2 exposes VDW/SAS surface selection with `-pbe-surface`; default
  behavior remains `vdw`, and `-pbe-probe` only changes the SAS radius inflate.
- PBE Surface Phase adds `-pbe-surface=ses`: it marks probe-excluded voxels
  and then flood-fills solvent-accessible grid points from the outer boundary.
  The final solute dielectric is the complement of accessible solvent, so
  buried probe-inaccessible cavities are handled without changing the solver.
- PBE Accuracy Parity Phase 3 adds `-pbe-surface=delphi-ses`,
  `-pbe-ion-access=stern`, `-pbe-ion-radius`, and `-pbe-map=ion`. The
  DelPhi-like SES mode keeps the SES flood-fill center classification but
  refines dielectric midpoint faces with 2x2 samples. Stern ion accessibility
  uses a separate ion probe so salt screening is applied only in
  ion-accessible solvent; the default remains `-pbe-surface=vdw` with
  `-pbe-ion-access=dielectric` for legacy regression stability.
- PBE nonlinear parity adds `-pbe-model=linear|nonlinear`. Keep the default
  `linear`; nonlinear currently seeds from the selected linear off/focus/vcycle
  solve and then runs damped Picard/SOR or the newer damped Newton path on the
  final fine grid. Salt-free nonlinear requests intentionally fall back to the
  linear effective solve with `nl_status=linearized_no_salt`; vacuum endpoints
  in `-pbe-dgsolv` and `-pbsa` remain linear/no-salt. Public nonlinear knobs
  include `-pbe-nl-solver=picard|newton`, `-pbe-nl-inner=auto|sor|vcycle`,
  `-pbe-nl-multigrid=off|fas`, `-pbe-nl-max-iter`,
  `-pbe-nl-inner-iter`, `-pbe-nl-stagnation`, `-pbe-nl-block`,
  `-pbe-nl-dphi`, `-pbe-nl-relax`, and `-pbe-nl-clip`.
- PBE Accuracy Parity Phase 4 adds discretization-lock diagnostics and knobs:
  `-pbe-boundary=coulombic|dipolar|zero`,
  `-pbe-charge-assignment=trilinear|delphi`, `-pbe-surface=delphi-vdw`,
  `-pbe-nl-continuation`, `-pbe-nl-line-search`, and `-pbe-nl-residual`.
  Keep the defaults legacy-compatible (`coulombic`, `trilinear`, `vdw`,
  `linear`). `pbe_energy.csv` appends grid origin/center/spacing, charge
  conservation/moment, boundary RMS/max, linear residual, nonlinear residual,
  continuation steps, and backtracks for DelPhi parity debugging.
- Hidden DelPhi PQR grids with `-pbe-test-grid-size` use pyDelPhi
  `COORDLIMITS` auto-centering, i.e. the radius-expanded boundary-box center,
  to match DelPhi `acenter=NA`; public adaptive grids still use
  radius-expanded bounds. `delphi-vdw` uses DelPhi-compatible discrete media
  indices at positive x/y/z half-grid faces: each midpoint is classified by VDW
  atom radius, then the PQR parity path applies the deterministic VDW-to-MS
  midpoint-media remap controlled by hidden `-pbe-test-vdwms-remap=auto|off`.
  For PQR parity, VDW media assignment uses the explicit radius read from
  PQR/SIZ and does not apply small-atom floors; Bondi floors remain a fallback
  for non-PQR `delphi-bondi` paths without explicit parity radii.
- `pbe_energy.csv` appends DelPhi-parity corrected RF diagnostics:
  solvent/vacuum phase RF energies, atom/grid corrected RF energies, the active
  reaction-energy source, and analytic Coulomb energy.
- PBE Accuracy Parity Roadmap adds opt-in Gaussian/RPB/focus-probe plumbing:
  `-pbe-biomodel=pbe|rpbe`, `-pbe-dielectric=two-dielectric|gaussian`,
  `-pbe-surface=gaussian|gaussian-cutoff`, Gaussian density knobs, and
  `-pbe-radius-convention=auto|native|delphi-bondi`. Defaults remain legacy.
  `-pbe-focus=delphi` is a parity/test path that writes `pbe_probe_sites.csv`
  when probe residue/site fixture options are supplied.
- PBE Accuracy Parity Phase 5 moves Newton from a diagonal-only update to a
  full-stencil Jacobian correction solved by SOR smoothing. The CSV appends
  `newton_inner_iterations`, residual-before/after, step norm, line-search
  status, and Gaussian density/epsilon audit columns. Barnase `nwt` is still a
  strict-parity diagnostic target; do not mark it OK unless `nl_status=ok` and
  the residual/DelPhi tolerance gates pass.
- PBE Accuracy Parity Phase 5b separates outer Newton steps
  (`-pbe-nl-max-iter`) from correction sweeps (`-pbe-nl-inner-iter`) and adds
  `-pbe-nl-stagnation`. The correction solver now reports residual RMS/max for
  `J delta = residual`, Jacobian diagonal min/max, accepted lambda, rejected
  line-search steps, active nonlinear points, and `newton_stop_reason`.
  Expected not-converged statuses include explicit suffixes such as
  `nonlinear_not_converged_stagnated` and
  `nonlinear_not_converged_line_search_floor`.
- PBE Accuracy Parity Phase 5c makes `-pbe-nl-inner=vcycle` a real Newton
  correction solver for `J delta = residual`. The inner V-cycle uses zero
  correction boundary values, the same epsilon off-diagonal stencil, and a
  Jacobian diagonal of `eps_sum + screen*cosh(alpha*phi)`. CSV/summary/log
  append `newton_inner_levels`, `newton_inner_cycles`,
  `newton_inner_residual_before/after`, and `newton_inner_effective`. Keep
  Gaussian/RPB/focus probes report-only until barnase `nwt` passes strict.
- PBE Accuracy Parity Phase 5d adds hidden PQR-only root-cause diagnostics:
  `-pbe-parity-dump=off|summary|points|all`, `-pbe-test-nl-source`,
  `-pbe-test-newton-derivative`, and `-pbe-test-energy-convention`. These are
  valid only with `-pbe-test-pqr`; do not expose them as production workflow.
  Phase 5e makes `-pbe-test-nl-source=delphi-dimensionless` a full solve-units
  convention, not a source-only toggle: nonlinear solves use
  `psi=alpha*phi`, an alpha-scaled RHS, `sinh(psi)`, and convert back to
  physical `phi` before interpolation, map output, and reaction-energy
  selection. Without external DelPhi stage TSV evidence, keep this as a hidden
  parity/debug path rather than changing public nonlinear defaults.
  Phase 5f keeps that decision and makes stage-oracle TSV comparison the
  acceptance boundary: `pbe_delphi_reference_smoke.ps1 -StageTsv <file>` checks
  grid/charge, epsilon/screen, boundary, linear residual, nonlinear source,
  Newton correction/derivative, and energy fields in order, reporting the first
  failing stage before any final energy verdict.
  `pbe_parity_stage.csv` records grid/charge/epsilon/screen/source/residual,
  Newton correction, derivative, and reaction-energy convention fields.
  `pbe_parity_points.csv` is only for `points|all` and dumps top residual
  points plus deterministic sample points.
- PBE Accuracy Parity Phase 5g keeps nonlinear/Newton work paused until the
  linear barnase mismatch is explained. When `-pbe-test-pqr` points at barnase
  and `-pbe-parity-dump` is enabled, TITAN now writes hidden forensic CSVs:
  `pbe_barnase_forensic_inputs.csv`, `pbe_barnase_surface_audit.csv`,
  `pbe_barnase_operator_audit.csv`, and `pbe_atom_phi_audit.csv`. The new
  `test/pbe_barnase_discretization_matrix.ps1` varies charge assignment,
  surface, radius convention, and boundary mode for salt 0/0.15/2.0 as a
  report-only root-cause matrix; do not treat mismatch as a failing gate
  without a real DelPhi stage TSV oracle.
- PBE Accuracy Parity Phase 5h adds the hidden PQR-only
  `-pbe-test-energy-convention=delphi-iscm` path for DelPhi two-dielectric
  induced-surface-charge reaction-field energy. It records ISCM boundary
  counts, induced-charge norms, water/vacuum/corrected RF energies in
  `pbe_energy.csv`, `pbe_parity_stage.csv`, and barnase operator audit output;
  `-pbe-parity-dump=all` also writes `pbe_iscm_surface_charge_samples.csv`.
  The selected corrected RF follows DelPhi's VDW water-phase ISCM convention;
  vacuum grid RF is retained as an audit field rather than subtracted in that
  hidden mode.
  Keep this as a DelPhi oracle/parity convention, not the PBSA/native default.
  The companion scripts are `test/pbe_iscm_energy_smoke.ps1` and
  `test/pbe_barnase_field_oracle_smoke.ps1`.
- PBE Accuracy Parity Phase 5i keeps nonlinear/Gaussian/FAS/focus paused and
  adds hidden barnase field-oracle tooling. `-pbe-parity-map-dump=off|samples|full`
  is valid only with `-pbe-test-pqr`; `samples` writes
  `pbe_field_samples.csv`, `pbe_boundary_samples.csv`, and
  `pbe_charge_grid_audit.csv`. The test-only replay knobs
  `-pbe-test-charge-grid=native|oracle`, `-pbe-test-surface-map=native|oracle`,
  `-pbe-test-midpoint-epsilon=native|oracle`, and
  `-pbe-test-boundary-map=native|oracle` read `pbe_oracle_*.csv` files from the
  run directory and isolate charge-grid, surface, midpoint-epsilon, and
  boundary-map drift one stage at a time. Do not expose these as production
  PBE/PBSA behavior. The companion scripts are
  `test/pbe_barnase_map_oracle_smoke.ps1` and
  `test/pbe_barnase_stage_replay_smoke.ps1`.
- PBE Accuracy Parity Phase 5j adds `test/pbe_delphi_export_maps.ps1`, which
  calls the local `.codex_tmp/pyDelPhi` checkout as an external oracle for
  barnase linear salt 0/0.15/2.0 and writes optional oracle CSVs under
  `test/tmp`. Phase 5k hardens that extractor: pyDelPhi VDW `surfmap` is
  normalized as `1=solvent, 0=solute`, the replayable
  `pbe_oracle_surface_map.csv` is a sparse list of `inside=1` solute rows, and
  sampled boundary data is written as `pbe_oracle_boundary_samples.csv` rather
  than being advertised as a complete replay map. If pyDelPhi cannot expose the
  needed maps or its dependencies are unavailable, the extractor writes a
  report-only CSV instead of fabricating golden data.
  Phase 5l switches `-pbe-surface=delphi-vdw` from the old fractional/harmonic
  face approximation to DelPhi-style discrete positive-face midpoint epsilon.
  `test/pbe_delphi_export_maps.ps1` writes
  `pbe_oracle_midpoint_epsilon_samples.csv` by default and a full replayable
  `pbe_oracle_midpoint_epsilon.csv` for one salt, default 0.15 M, with
  `-FullMidpointEpsilon` and optional `-FullMidpointSalt`; stage replay can
  then run `-pbe-test-midpoint-epsilon=oracle`.
  Phase 5m locks radius/media assignment: barnase forensic output now includes
  `pbe_barnase_radius_media_audit.csv` plus per-axis midpoint face counts in
  `pbe_barnase_surface_audit.csv`. The DelPhi PDB+CRG+SIZ test adapter strips
  inline `!`/`#` comments before taking SIZ radii, matching pyDelPhi rows such
  as `HG SER 0.6000 ! ... 0.0001`; without this, generated PQR files assign
  tiny placeholder H radii and poison midpoint-epsilon comparisons. Phase 5n
  adds a clean-room VDW-to-MS remap layer and VDWMS audit fields
  (`vdwms_enabled`, boundary counts, exposed-point count, and midpoint RMS).
  Phase 5o adds `test/pbe_midpoint_oracle_roundtrip_smoke.ps1` and the
  `midpoint_eps_sum_*` audit fields, separating oracle layout mistakes from
  native remap mismatch. The full pyDelPhi midpoint replay now roundtrips, but
  native `delphi-vdw` still has a `midpoint_epsilon` mismatch on barnase, so
  keep barnase strict parity disabled until the VDW-to-MS remap itself matches
  the oracle face counts.
  Phase 5q adds hidden `-pbe-test-midpoint-diff=off|samples|full`, which writes
  `pbe_midpoint_face_diff.csv` from `pbe_oracle_midpoint_epsilon.csv` with only
  the differing positive faces plus raw/final media owner, signed encoded
  media/entity state, queue iteration, remap action, nearest atom distance,
  and exposed-point distance. Phase 5r moves the `delphi-vdw` remap internals
  to a signed discrete epsilon state, but barnase still needs the exact
  pyDelPhi exposed/contact remap geometry before strict midpoint parity can
  pass. Use
  `test/pbe_midpoint_face_diff_smoke.ps1` report-only first; `-Strict` is for
  the moment native face diffs are expected to be zero.
  Phase 5s adds VDWMS state capture on both sides: the pyDelPhi exporter accepts
  `-FullVdwmsState`, and native full midpoint diff writes
  `pbe_native_midpoint_state.csv` plus `pbe_vdwms_action_histogram.csv`.
  `test/pbe_barnase_stage_replay_smoke.ps1` uses those artifacts to report the
  first failing VDWMS stage as `raw_vdw`, `exposed_sas`, `boundary_queue`,
  `signed_remap`, or the later field/energy stage.
  `test/pbe_delphi_reference_smoke.ps1 -StageTsv <file>` now also resolves
  `test_center=auto` from `grid_center_x/y/z` in the oracle stage row; without
  that center lock, barnase native midpoint comparisons can report a false
  `raw_vdw` mismatch even when the remap algorithm is correct.
  `test/pbe_barnase_field_oracle_smoke.ps1 -OracleRoot <dir>` now compares
  normalized surface/epsilon, boundary, water/vacuum/reaction phi, ISCM
  samples, and writes `pbe_barnase_root_cause_report.csv`.
  `test/pbe_barnase_timing_compare.ps1` is a report-only TITAN vs pyDelPhi
  timing harness for barnase linear salt 0/0.15/2.0. Keep these oracle CSVs
  hidden PQR parity inputs only. Arg files use `-flag=value`; direct CLI
  invocations use `-flag value`.
- PBE Accuracy Parity Phase 5t focuses on DelPhi boundary condition and field
  units after midpoint epsilon is locked. Hidden PQR-only flags
  `-pbe-test-boundary-convention=native-real|delphi-grid|oracle` and
  `-pbe-test-field-units=native|delphi-dimensionless` isolate boundary
  sampling from solver units. `delphi-grid` uses DelPhi atom-grid distance,
  Debye-Huckel scaling, and exdi/indi divided by `epkt`; `oracle` replays
  `pbe_oracle_boundary_map.csv`. The pyDelPhi exporter accepts
  `-FullBoundaryOracle`, and the new reports are
  `test/pbe_boundary_formula_matrix.ps1`,
  `test/pbe_boundary_replay_smoke.ps1`, `pbe_boundary_diff.csv`, and
  `pbe_operator_stage.csv`. Do not change PBSA/native defaults from these
  diagnostics unless boundary and water-phi oracle stages both justify it.
- PBE Accuracy Parity Phase 5u moves the oracle comparator to water-phi
  operator setup. The pyDelPhi exporter accepts `-FullOperatorOracle` and writes
  `pbe_oracle_operator_stage.csv`, `pbe_oracle_sor_trace.csv`, and
  `pbe_oracle_solver_summary.csv`. Hidden PQR flags
  `-pbe-test-operator-convention=native|delphi`,
  `-pbe-test-rhs-map=native|oracle`, `-pbe-test-screen-map=native|oracle`, and
  `-pbe-test-denom-map=native|oracle` isolate RHS, screening, and denominator
  mismatches before SOR changes. Use `test/pbe_water_phi_operator_matrix.ps1`
  and `test/pbe_sor_trace_smoke.ps1` for this stage; keep arg-file syntax as
  `-flag=value` and direct CLI syntax as `-flag value`.
- PBE Accuracy Parity Phase 5v isolates charge-grid assignment before SOR
  update semantics. The pyDelPhi exporter accepts `-FullChargeGrid` and writes
  `pbe_oracle_charge_grid_full.csv`; TITAN hidden PQR diagnostics add
  `-pbe-test-charge-grid-diff=off|samples|full`, `pbe_charge_grid_diff.csv`,
  and `pbe_charge_grid_diff_summary.csv`. `pbe_sor_trace.csv` appends
  `rhs_delta_vs_oracle` and `charge_grid_source`. Use
  `test/pbe_charge_grid_full_diff_smoke.ps1`,
  `test/pbe_charge_grid_replay_smoke.ps1`,
  `test/pbe_sor_update_lock_smoke.ps1`, `test/pbe_water_phi_lock_smoke.ps1`,
  and `test/pbe_iscm_transition_smoke.ps1` for this stage. Do not change ISCM
  energy convention until charge-grid replay and water-phi samples pass.
- PBE Accuracy Parity Phase 5w locks the water field oracle precision before
  touching ISCM or final energy. `test/pbe_delphi_export_maps.ps1` accepts
  `-FullFieldOracle` and writes
  `pbe_oracle_field_samples_fullprecision.csv` with
  `field_source=solver_capture|ctx|cube_fallback` and
  `field_precision_status`. `test/pbe_field_precision_smoke.ps1 -Strict`
  must pass before using `test/pbe_water_phi_lock_smoke.ps1 -Strict`.
  Current barnase evidence moves the first failing stage past `water_phi_raw`
  and `ISCM` samples to `final_energy`, so the next parity work should isolate
  final reaction-energy selection/convention rather than changing the solver.
- PBE Accuracy Parity Phase 5x locks the final DelPhi reaction-energy
  convention for hidden PQR parity. The pyDelPhi extractor accepts
  `-FullIscmEnergyOracle` and writes `pbe_oracle_energy_components.csv`,
  `pbe_oracle_iscm_full.csv`, and `pbe_oracle_iscm_energy_terms.csv`.
  TITAN `delphi-iscm`/`delphi-corrected` now build DelPhi-style scaled
  surface-charge positions natively; `pbe_oracle_iscm_full.csv` is only a
  compare/audit oracle (`iscm_position_rms_vs_oracle`, sample diffs), not an
  input to final ISCM accumulation. Use
  `test/pbe_final_energy_convention_smoke.ps1 -Strict` with a full oracle root
  before claiming barnase final-energy parity. Public PBE/PBSA defaults remain
  on the native atom/grid convention.
- PBE nonlinear/Gaussian/RPB/focus lock resumes after accepting the small
  native ISCM position residual as diagnostic-only. `-pbe-model=nonlinear
  -pbe-nl-solver=newton` now keeps DelPhi PQR parity rows on SOR when
  `-pbe-nl-inner=auto`, and a rejected Newton correction falls back to a native
  NWT SOR block using the DelPhi linearization
  `denom += screen*cosh(phi)` and `rhs += screen*(phi*cosh(phi)-sinh(phi))`
  instead of the older Picard fallback. Hidden PQR parity dump `all` also writes
  `pbe_nwt_trace.csv`, and `pbe_energy.csv`/`pbe_parity_stage.csv` append
  nonlinear electrostatic-stress, osmotic-pressure, and total-correction
  columns for pyDelPhi NWT energy-component comparison. Gaussian density for
  `-pbe-dielectric=gaussian`, `-pbe-surface=gaussian|gaussian-cutoff`,
  `-pbe-charge-assignment=gaussian`, and `-pbe-biomodel=rpbe` uses the
  multiplicative DelPhi-style union at grid centers and positive midpoint
  faces. Phase 6b extends the pyDelPhi exporter with
  `pbe_oracle_nwt_trace.csv` and `pbe_oracle_nwt_solver_summary.csv`; use
  `test/pbe_nonlinear_nwt_trace_compare.ps1` to fail in fixed order
  `nwt_units -> nwt_rhs -> nwt_screen -> nwt_denom -> nwt_update` before
  tuning residuals. Phase 6c makes the DelPhi NWT parity path use pyDelPhi
  block-stop semantics directly: `-pbe-nl-block` is the RMS/max-dphi block
  size, `pbe_nwt_solver_progress.csv` records per-block convergence, and
  `test/pbe_nwt_solver_progress_smoke.ps1` compares that curve against the
  oracle before nonlinear energy strict checks. Do not reintroduce TITAN
  line-search/stagnation early stops in `delphi_nwt_parity_active()`. New
  Phase-6 gates include
  `test/pbe_nonlinear_nwt_trace_smoke.ps1`,
  `test/pbe_nwt_solver_progress_smoke.ps1`,
  `test/pbe_nonlinear_energy_components_smoke.ps1`,
  `test/pbe_gaussian_density_oracle_smoke.ps1`,
  `test/pbe_rpbe_source_smoke.ps1`, and
  `test/pbe_focus_boundary_oracle_smoke.ps1`. Keep RPB and focus/probe strict
  gates tied to the pyDelPhi oracle smoke scripts until their stage diagnostics
  pass. Phase 6d adds `-FullGaussianOracle`, `-FullRpbOracle`,
  `-FullFocusProbeOracle`, and `-Fixture` to
  `test/pbe_delphi_export_maps.ps1`, plus `test/pbe_timing_compare.ps1` for
  report-only TITAN vs pyDelPhi timing rows. Gaussian density/charge oracle is
  strict-capable on the 5tif Gaussian-density fixture. RPB now computes a
  native analytic Coulomb/grad-Coulomb/source audit instead of using the vacuum
  field, and `test/pbe_rpbe_source_smoke.ps1 -Strict` passes after locking
  `test_center=auto` to pyDelPhi's radius-expanded `COORDLIMITS` grid center.
  With parity dump enabled, `pbe_rpbe_source_top.csv` records the strongest
  native source cells for native-vs-oracle diagnostics.
  Focus/probe oracle export now builds a parent phimap and uses pyDelPhi
  `focusing` boundary syntax; probe energy may still be missing in the oracle
  and should not be claimed strict without `pbe_oracle_focus_probe.csv`
  containing finite probe energy.
  Phase 6e adds default-off hidden timing diagnostics:
  `-pbe-timing=off|summary|detail` writes `pbe_timing.csv` with setup,
  surface, epsilon, boundary, solver, energy, map I/O, RPB source, focus setup,
  probe, total, atom count, thread count, workflow scope, iteration, block, and
  status buckets. Detail mode also records focus/probe sub-buckets for parent
  solve, child setup, boundary interpolation, child surface, child epsilon,
  child solver, probe energy, and audit output. Public/default PBE and PBSA
  outputs are unchanged. `test/pbe_timing_compare.ps1` now consumes that file,
  reports `speed_ratio_titan_vs_pydelphi`, records `titan_threads`,
  `pydelphi_threads`, and records the `bottleneck_bucket` per fixture. The
  local golden CSV has harness-only
  `surface_cutoff`, `pbe_timing`, and `parity_gate` columns; new Gaussian,
  RPB, and focus/probe matrix rows start as `report-only` until their
  source/operator/final-energy oracle data is finite. Use
  `test/pbe_parity_difference_matrix.ps1 -Fixtures <a,b>` or `-MaxRows` for
  targeted expansion checks. When passing multiple fixtures to the exporter,
  prefer current-shell array syntax such as
  `& .\test\pbe_delphi_export_maps.ps1 ... -Fixture barnase_nwt,barnase_nwt_salt2`;
  nested `powershell -File` can pass the comma list as one string. A full
  all-oracle refresh can exceed 15 minutes locally, so regenerate targeted
  oracle groups when only one strict file is missing.
  Phase 6f closes the timing harness linear row by using
  `test_energy_convention=delphi-iscm` and `temp_K=297.15` for
  `barnase_linear_salt015`; the row is now `pass` instead of
  `mismatch_report_only`. `test/pbe_parity_difference_matrix.ps1 -Fixtures`
  now honors explicitly requested strict rows, not only report-only matrix
  rows. `test/pbe_timing_compare.ps1 -RepeatCount N` runs repeat subfolders and
  writes `pbe_timing_compare_repeats.csv` plus a median
  `pbe_timing_compare.csv`. The first optimization pass keeps formulas intact:
  red/black SOR/NWT sweeps step through same-color cells directly, Gaussian
  operator setup reuses the center density computed during surface
  classification, and RPB source audit precomputes per-atom Gaussian factors.
  The refreshed timing snapshot is: linear `19.50s` and pass, NWT `276.03s`,
  Gaussian-cutoff `30.85s`, RPB `216.25s`, focus/probe `45.27s`.
  Phase 6g records fair thread/scope metadata and closes the locked RPB source
  bottleneck. The pyDelPhi exporter writes `pydelphi_threads=2`; RPB source
  audit now uses the caller executor/thread count, pre-skips Gaussian-gradient
  contributions outside the locked density cutoff, and reuses solvent RPB
  audit/source/Coulomb maps for energy/top diagnostics. The refreshed timing
  snapshot is: linear `19.20s` vs pyDelPhi `21.87s`, NWT `295.45s` vs
  `142.22s`, Gaussian-cutoff `42.71s` vs `48.95s`, RPB `42.96s` vs `82.96s`,
  and focus/probe full workflow `27.14s` vs pyDelPhi child oracle `4.45s`.
  The focus/probe row remains scope-mismatched: TITAN times 1697 atoms and the
  full parent/child/probe workflow, while pyDelPhi times a 599-atom child oracle;
  TITAN detail reports parent solve `14.35s` and child solver `0.65s`.
  Phase 6h/6i adds hidden NWT/timing controls:
  `-pbe-nwt-progress=auto|off|strict`,
  `-pbe-nwt-kernel=auto|legacy|delphi-half`, and
  `-pbe-nwt-homo-kernel=auto|packed-list|row-span|avx2|salt-tiled`. Strict gates keep
  progress output by default; timing rows can use `pbe_nwt_progress=off` to
  avoid charging per-block source/residual/progress diagnostics to TITAN. The
  `delphi-half` NWT kernel now has a persistent workspace across blocks: split
  RHS/epsilon/screen caches, half-neighbor interior-cell packs, separate
  homogeneous/heterogeneous update loops, row-span/fixed-offset homogeneous
  streaming, and a final flush back to `phi_`. `auto` selects this path for
  pyDelPhi/NWT fair timing and now promotes `salt-tiled` after strict/kernel
  smokes pass. That path tiles homogeneous salt spans, hoists color-specific
  offsets/pointers, and uses a double-precision polynomial `sinh`/`cosh` fast
  branch for small clipped potentials; `packed-list`, `row-span`, and `avx2`
  remain forced bisect modes. Strict progress rows keep the full-grid trace
  path. The latest full median snapshot has barnase NWT at `106.64s` total /
  `89.09s` solver versus pyDelPhi `124.09s` / `94.54s`, with
  `nwt_half_update=57.45s`, `nwt_homo_update=41.93s`, and
  `nwt_homo_ns_per_update=4.24`.
  `pbe_timing.csv` records NWT detail buckets (`nwt_half_update`,
  `nwt_coeff_pack`, `nwt_homo_update`, `nwt_hetero_update`, `nwt_reduce`,
  residual/source/progress I/O, `nwt_homo_stencil`, `nwt_homo_salt_math`,
  `nwt_homo_write`, `nwt_thread_wait`, `nwt_block_overhead`), homogeneous counters (`nwt_homo_cells`,
  `nwt_homo_linear_cells`, `nwt_homo_salt_cells`, `nwt_homo_updates`,
  `nwt_dispatch_calls`, `nwt_homo_ns_per_update`), and Gaussian epsilon detail buckets
  (`epsilon_center_density`, `epsilon_midpoint_density`,
  `epsilon_ion_exclusion`, `epsilon_finalize`). The current 5tif
  gaussian-cutoff full row is `23.62s` total / `8.10s` epsilon versus
  pyDelPhi `57.50s`; Gaussian epsilon is no longer the row bottleneck after
  atom-SoA plus parallel center/midpoint/ion/finalize passes.
  `test/pbe_timing_compare.ps1` also emits a report-only
  `titan_child_scope` row for focus/probe, now using the true child atom count.
  The current full workflow is `32.85s` with parent solve `16.31s`, while the
  child-scope estimate is `2.31s` for 599 child atoms versus pyDelPhi child
  oracle `8.66s`; use this only for scope-aware diagnostics, not as a
  replacement for production full-workflow timing. The timing harness default
  is now `RepeatCount=2` median with a comparable-row speed gate of 1.05.
  Oracle exporter/timing scripts append black-box provenance to
  `test/pbe_oracle_manifest.csv`; production TITAN must never read that file.
  `-pbe-linear-kernel=half-packed` exists as a forced hidden bisect knob, but
  `auto` intentionally remains legacy because the packed linear focus path was
  slower on this workload.
- PBE snapshot stability diagnostics are default-off and separate from parity
  oracle inputs. `-pbe-stability-audit=off|summary|detail` writes
  `pbe_stability_audit.csv`; detail mode also writes
  `pbe_stability_atom_samples.csv` with nearest atom-to-integer/half-grid-plane
  samples. `-pbe-grid-anchor=auto|per-frame|first-frame|fixed-center` is a
  hidden trajectory-stability bisect knob, not a public output change.
  `-pbe-stability-mode=off|diagnose|anchor|hysteresis|ensemble` is also hidden:
  `diagnose` forces summary audit only, `anchor` uses first-frame anchoring,
  `hysteresis` applies opt-in solvent-medium hysteresis near atom-surface
  thresholds, and `ensemble` writes deterministic jitter mean/std diagnostics
  without replacing `pbe_energy.csv`.
  `test/pbe_snapshot_stability_scan.ps1` generates a synthetic grid-boundary
  smoke when no snapshots are supplied and can scan user PQR snapshots with
  `-PqrPaths`. Use its `suspect_stage`/`scan_label` columns to distinguish
  atom-on-grid-plane candidates from grid recentering, surface/midpoint,
  ion-access, or ISCM-boundary jumps.
  `test/pbe_grid_boundary_pydelphi_compare.ps1` adds a report-only pyDelPhi
  synthetic compare when the local pyDelPhi environment is present; keep
  pyDelPhi out of production TITAN paths.
- The strict DelPhi scripts are diagnostic gates, not a reason to loosen
  tolerances. `test/pbe_delphi_strict_linear_smoke.ps1` and
  `test/pbe_barnase_nonlinear_smoke.ps1` can run non-strict to emit reports,
  and `-Strict` turns those reports into failing gates once the linear
  discretization is locked.
- PBE Accuracy Parity Phase starts with nested-grid focusing via
  `-pbe-focus=auto`, `-pbe-multigrid=focus`, `-pbe-focus-levels`, and
  `-pbe-focus-margin`. This is a coarse-to-fine focusing path that interpolates
  the previous level potential into the next level boundary/initial guess.
- PBE Accuracy Parity Phase 2 adds `-pbe-multigrid=vcycle` plus
  `-pbe-mg-levels`, `-pbe-mg-pre`, `-pbe-mg-post`, and
  `-pbe-mg-coarse-iter`. The V-cycle path builds a geometric cuboidal
  hierarchy, runs residual restriction/correction prolongation, then finishes
  with the same fine-grid SOR convergence polish used by the baseline solver.
  Keep `-pbe-multigrid=auto` mapped to nested focusing unless the public CLI
  contract is deliberately changed.
- PBE Phase 3 regression lives in `test/pbe_regression.ps1` with PQR fixtures
  under `test/pbe_inputs_local/`. The hidden `-pbe-test-pqr` adapter is for
  this harness only; do not document it as a public PBE input path.
- PBE Accuracy Parity Phase 3 adds `test/pbe_stern_smoke.ps1`,
  `test/pbe_delphi_surface_smoke.ps1`, and
  `test/pbe_delphi_parity_matrix.ps1`. The matrix writes internal parity
  results and treats external DelPhi/pyDelPhi reference energies as optional
  `delphi_E_rxn_kJmol` columns when available.
- The real DelPhi/pyDelPhi golden harness is
  `test/pbe_delphi_generate_golden.ps1` plus
  `test/pbe_delphi_reference_smoke.ps1`. It reads the external local checkout
  under `.codex_tmp/pyDelPhi` or a DelPhi TSV and writes TITAN golden CSV only
  when called with `-UpdateReference`; do not vendor or copy pyDelPhi source
  into TITAN. `PBE_CLEAN_ROOM_POLICY.md` is the forward-looking source of truth:
  production PBE is self-contained and may only be validated against black-box
  oracle outputs, timing rows, hashes, and stage CSVs. Hidden parity-only knobs
  `-pbe-test-grid-size` and
  `-pbe-test-center` are valid only with `-pbe-test-pqr` for reproducing
  DelPhi scale/grid/center settings. Barnase root-cause reports now include
  `test/pbe_barnase_linear_stage_smoke.ps1`,
  `test/pbe_barnase_energy_convention_smoke.ps1`,
  `test/pbe_barnase_forensic_smoke.ps1`,
  `test/pbe_barnase_atom_phi_smoke.ps1`,
  `test/pbe_barnase_discretization_matrix.ps1`,
  `test/pbe_barnase_parity_stage_report.ps1`,
  `test/pbe_delphi_export_maps.ps1`,
  `test/pbe_nonlinear_source_convention_smoke.ps1`, and
  `test/pbe_nonlinear_units_matrix_smoke.ps1`, and
  `test/pbe_newton_derivative_smoke.ps1`.
- PBE backend routing is explicit. Bare `-pbe`, `-pbe=grid`, and aliases
  `-pbe=delphi|pydelphi|fd` use the established finite-difference grid solver.
  `-pbe=tabi` is a separate clean-room C++23 TABI v1 backend under
  `titan_pbe/`; it is linear PB only and rejects nonlinear/RPB, Gaussian,
  focus/multigrid, PBSA/decomposition, and grid-map modes. Keep TABI validation
  black-box-only against papers/oracle outputs; do not vendor or call external
  TABI-PB/NanoShaper at runtime. TABI now defaults to the JCP 2013 coupled
  boundary-integral formulation with a clean-room triangle surface
  (`-pbe-tabi-geometry=triangle`). Self-panel analytic singular integration is
  available through the hidden `-pbe-tabi-near-integration` validation knob;
  one-ring Duffy/CPV is now a guarded validation path that must not worsen the
  Born matrix, and `auto` remains centroid until broader analytic gates pass.
  The old scalar `titan-v1` formulation and the
  shared `TitanSurface` point-panel path remain hidden debug paths.
  `jcp2013-centroid` is the direct paper-style path. `jcp2013-treecode` runs
  the coupled JCP 2013 treecode operator and should remain gated by
  direct-vs-treecode accuracy tests before becoming the default auto path.
  `tabi2-nodepatch` is still an experimental/report mode. The legacy scalar
  treecode remains available when explicitly pinning
  `-pbe-tabi-formulation=titan-v1`.
  TABI timing now separates surface-neighbor, mesh, solid-angle, RHS, GMRES,
  matvec, preconditioner, BIE formulation, and tree moment/traversal buckets in
  `pbe_timing.csv`; `-pbe-timing=summary|detail` also writes
  `pbe_memory_profile.csv` for grid-vs-TABI RAM estimates/RSS samples, including
  both clamped and signed RSS deltas, and
  `detail` writes `pbe_tabi_profile.csv` with one aggregate row per solve.
  `test/pbe_born_ion_matrix.ps1` is the primary TABI analytic gate; it records
  Born status, max error by surface density, convergence slope, and fallback
  use across charge/radius/salt/density cases. `test/pbe_grid_born_convergence.ps1`
  is the separate grid diagnostic for Born convergence, vacuum-reference
  reporting, grid center anchoring, time, and RAM.
  `test/pbe_tabi_vs_grid_compare.ps1 -ProfileMemory` compares energy/time/RAM
  for grid versus TABI; use `-FixtureSet barnase` or
  `test/pbe_tabi_vs_grid_barnase_compare.ps1` for the report-only barnase
  rows built from the local pyDelPhi fixture cache.
  Every TABI solve writes
  `pbe_tabi_formulation_audit.csv` with formulation, energy projection,
  unknown count, BIE system, diagonal constants, kernel convention, source and
  solution norms, dielectric contrast, kappa, boundary energy, Born reference,
  surface_generator, discretization, near_integration/self/one-ring ASI
  counters, near_kernel_policy, and preconditioner columns. Hidden
  accuracy-bisect knobs are
  `-pbe-tabi-formulation=titan-v1|jcp2013|jcp2013-centroid|jcp2013-treecode|tabi2-nodepatch|calibrated-debug`,
  `-pbe-tabi-energy-projection=boundary|born-check|charge-potential-debug`,
  `-pbe-tabi-near-integration=auto|centroid|analytic_singular|duffy_adaptive`,
  `-pbe-tabi-precond=auto|none|diagonal|block`, and
  `-pbe-tabi-surface-points=<24..4096>`.
  The TABI energy source must be explicit (`tabi_boundary_integral`,
  `tabi_born_check`, `tabi_charge_potential_debug`, `tabi_calibrated_debug`, or
  `tabi_born_fallback`); do not reintroduce silent Born fallback behavior.
- PBE `-pbe-dgsolv` is opt-in because it runs six solves per frame for
  two-molecule systems: complex, mol1, and mol2 each get solvent and vacuum
  solves. Keep `deltaG_solv` and component energies in both kJ/mol and
  kcal/mol when changing this path.
- `-pbsa` is the public MM/PBSA wrapper on top of PBE. It implies `-pbe` and
  `-pbe-dgsolv`, uses the native force-field elec/vdw provider plus
  Shrake-Rupley-style SASA and optional interaction entropy. PBSA IE uses the
  shared Duan-Zhang interaction entropy core also used by the `-entropy`
  module, then adds PBSA-specific ESS/block/filter diagnostics. It writes
  `mm_pbsa_energy.csv` plus `mm_pbsa_summary.dat` in both kJ/mol and kcal/mol.
  `-pbsa-decomp=residue` writes additive component residue decomposition from
  the existing six endpoint solves; it is not residue deletion/mutation PBSA.
  PBSA also writes `mm_pbsa_by_trajectory.csv`,
  `mm_pbsa_ie_diagnostics.csv`, and `mm_pbsa_performance.dat` unless disabled
  with `-pbsa-perf=off`.
- GBSA, PBSA, and TABI share `titan_surface.*` for Shrake-Rupley exposed
  surface workspaces. `titan_sasa.*` remains the compatibility wrapper for
  per-atom SASA/SA nonpolar callers, so the same probe, 96-point sphere, and
  gamma conversion still drive GBSA/PBSA workflows. TABI consumes the exposed
  panel list `{position, normal, area, atom_index}` only when the hidden
  point-panel debug path is selected; the default TABI path owns triangle mesh
  reconstruction so it can compute local solid-angle/diagonal diagnostics.
  `TitanSurface::Options`
  has an internal `acceleration` selector (`auto|legacy|reordered_padded`);
  `auto` uses the Morton-reordered padded-neighbor backend for larger target
  sets while keeping legacy as a parity oracle. Keep residue decomposition in
  the caller: GBSA injects per-atom nonpolar deltas into its existing residue
  aggregator, while PBSA groups the shared per-atom SASA deltas by residue
  alongside PB reaction-field terms.
- PBE still requires charges. Zero-radius polar hydrogens from supported force
  fields get a small Bondi-style surface radius fallback so AMBER-style
  zero-LJ hydrogens do not block PBSA.
- TITAN arg files use `-flag=value`; direct CLI invocations use
  `-flag value`. Keep examples in that convention, especially for decimal PBE
  values.

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
- for `-llps-paper`, keep the single-run bundle opt-in, emit `LLPS Paper Output
  Audit` in `log.txt`, and keep multi-condition comparisons in
  `plot_llps_paper.py --ensemble-manifest` rather than running multiple
  trajectories in the engine
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
  FRAP-like recovery proxies, opt-in `-llps-paper` manuscript summaries,
  detailed performance reporting, and `log.txt` run-context/output audit
  sections instead of a separate manifest file.
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

** C:\Users\ADMIN\OneDrive\Documents\antigravity\tools use rg in this folder
