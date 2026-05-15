# TITAN PBE TABI Roadmap

This roadmap tracks the remaining work needed to move the clean-room TABI
backend from the current validated core toward paper-level TABI 2013 and
TABI-PB 2.0 coverage.

Last updated: 2026-05-14

## Current Status

TABI is now useful as an internal linear-PB backend, but it is not yet a full
implementation of either TABI 2013 or TABI-PB 2.0.

Implemented and gated:

- CLI routing: `-pbe=tabi` is separate from the grid backend.
- Core formulation: JCP 2013-style coupled BIE unknowns `[phi, dphi/dn]`.
- Direct JCP path: `jcp2013` / `jcp2013-centroid`.
- Analytic diagonal constants and boundary-integral energy projection.
- Triangle surface path using TITAN-generated sphere/patch geometry.
- Self-panel analytic singular integration path.
- One-ring CPV/Duffy path now reports per-block acceptance/rejection for
  `K1/K2/K3/K4`, so guarded centroid fallbacks are auditable by block.
- One-ring-specific counters now split self decisions from one-ring decisions:
  `near_one_ring_k1_accepted/rejected` through
  `near_one_ring_k4_accepted/rejected`.
- `-pbe-timing=detail` now emits `pbe_tabi_near_block_audit.csv`, aggregated
  by `block=K1..K4`, `pair_type=self|edge|vertex|smooth_near`, and
  `orientation=same|flipped|mixed`, with local-reference error summaries and
  `candidate_worse_than_centroid` counts.
- The same audit now records explicit rejection reasons:
  `reference_nonconverged`, `orientation_sign_mismatch`, and
  `strict_ref_reject`; one-ring rejection should no longer be diagnosed only
  as a silent `centroid_guarded` fallback.
- `test/pbe_tabi_near_block_rank.ps1` ranks guarded buckets by
  `candidate_worse_than_centroid`, K2 orientation diagnostics, reference
  error, and rejection reasons.
- K2 source-normal CPV now has a separate guarded path: Laplace source-normal
  from oriented solid angle plus Yukawa smooth subtraction from Duffy. K2 is
  accepted only when the local reference converges, the orientation pair check
  passes, the candidate improves centroid by margin, and
  `candidate_ref_rel <= 1e-8`.
- Surface-quality diagnostics now exist for TABI triangle surfaces:
  `pbe_tabi_surface_quality.csv` is emitted under `-pbe-timing=detail`, and
  the hidden mesh selector `-pbe-tabi-sphere-mesh=auto|latlong|icosphere`
  keeps `auto=latlong`.
- Current priority remains K2 vertex/edge same-orientation because it is still
  the top ranked guarded bucket on mini-protein. K1/K3/K4 remain
  conservatively strict-guarded until Born/Kirkwood gates improve or stay flat.
- Focused near-kernel block smoke exists for self, edge-neighbor,
  vertex-neighbor, smooth near, and orientation-flipped canonical triangle
  pairs; it verifies finite CPV decisions while allowing guarded rejection when
  global analytic gates would otherwise worsen.
- JCP treecode experimental path: `jcp2013-treecode`.
- TABI timing/profile output and grid-vs-TABI memory profile.
- Primary Born/Kirkwood analytic gate:
  - surface density 96: max relative error <= 5%.
  - surface density 384: max relative error <= 1.5%.
  - surface density 768: max relative error <= 1%.
- Grid Born diagnostic is separated from TABI correctness.
- `test/pbe_kirkwood_oracle_selfcheck.ps1` cross-checks centered screened
  Kirkwood/Born closed form against the spherical-Bessel coefficient for
  radii `1.5, 2, 3, 5 A`, charges `+/-1,+/-2`, and salts
  `0, 0.05, 0.15, 0.5 M`.
- Off-center Kirkwood sphere rows now enforce centered salted cases after that
  selfcheck passes. Salted off-center and salted multi-charge rows remain
  report-only until their series tail/convergence checks are hardened.
- Off-center and multi-charge Kirkwood scripts now report adaptive oracle
  truncation diagnostics: `l_used`, `oracle_tail_abs`, `oracle_tail_rel`, and
  `oracle_status`. Salted multi-charge rows are enforced when tail and symmetry
  checks pass; salted off-center rows are enforced case-by-case when both oracle
  and solver tolerance pass.

Known limitations:

- TABI surface is still TITAN-generated triangle/sphere patch geometry, not a
  full MSMS/NanoShaper-equivalent SES pipeline.
- One-ring near-field CPV/Duffy is audit-visible and guarded; individual
  finite-part blocks may stay centroid when the CPV candidate fails local
  consistency guards or worsens analytic gates.
- Current one-ring policy is deliberately conservative: only block-local
  candidates that do not worsen Born/Kirkwood gates may be accepted. K1/K3/K4
  are intentionally strict-guarded in this slice after Born-gate checks showed
  tiny but real `duffy_adaptive` regressions when they were accepted.
- K2 is also still conservative: canonical local block tests accept K2, but
  protein-like buckets remain mostly guarded unless the stronger
  local-reference margin passes. This kept Born/Kirkwood analytic gates from
  regressing.
- The targeted salted off-center case `radius=5 A`, `offset=0.35`,
  `salt=0.05 M` is now classified as `surface_mesh_limited` for the legacy
  lat/long sphere mesh. Lat/long stayed above threshold even at `1536`
  requested points, while `icosphere` passed at `384` requested points
  (`rel_error_vs_analytic ~= 1.22%`) and improved further at `1536`
  (`~= 0.60%`). Surface-quality audit explains the split:
  lat/long `triangle_aspect_max` grows from about `6.86` to `12.82`, while
  `icosphere` stays near `1.17`.
- `auto` near integration still resolves to `centroid`.
- `jcp2013-treecode` exists but should not become the default until it passes
  direct-vs-treecode energy gates with the final near-field integration.
- `tabi2-nodepatch` is experimental/report-only.
- Block preconditioning is not yet the real TABI-PB 2.0 block preconditioner.
- 1A63 and protein paper comparisons remain report-only unless stable local
  PQR/mesh fixtures are available.

## Clean-Room Boundary

Production TITAN must remain self-contained C++23 code.

Allowed:

- Use the local TABI 2013 and TABI-PB 2.0 PDFs as algorithm references.
- Use external APBS, DelPhi, pyDelPhi, NanoShaper, or TABI-PB outputs as
  report-only oracle data in tests.
- Cache fixed input fixtures when their provenance is documented.

Not allowed:

- Vendor TABI-PB/NanoShaper source into TITAN.
- Call TABI-PB/NanoShaper/APBS/DelPhi/pyDelPhi at production runtime.
- Use Born fallback as accepted TABI energy.

## Stage 1 - Finish Near-Field ASI

Goal: remove the biggest remaining mathematical weakness in TABI: crude
near-field and one-ring treatment.

Tasks:

- Implement closed-form triangle integral for Laplace single-layer `G0`.
- Implement Yukawa `Gk` with singularity subtraction:
  `Gk = G0 + smooth(Gk - G0)`.
- Implement source-normal, target-normal, and mixed-normal one-ring kernels
  using Duffy adaptive quadrature with principal-value guards.
- Use the same near-field code path for:
  - direct JCP operator.
  - treecode leaf near-field.
  - energy projection near terms, where applicable.
- Keep JCP jump/diagonal constants on the LHS only; do not double-count them
  inside self integrals.
- Add per-kernel audit counters:
  - `near_kernel_policy`
  - `near_self_pairs`
  - `near_one_ring_pairs`
  - `near_duffy_pairs`
  - `near_centroid_pairs`
  - `near_fallbacks`
  - `duffy_subdivisions_max`
  - `duffy_subdivisions_mean`
  - `near_k1_accepted/rejected`, `near_k2_accepted/rejected`
  - `near_k3_accepted/rejected`, `near_k4_accepted/rejected`
  - `near_one_ring_k1_accepted/rejected`
  - `near_one_ring_k2_accepted/rejected`
  - `near_one_ring_k3_accepted/rejected`
  - `near_one_ring_k4_accepted/rejected`
  - `near_k1_policy`, `near_k2_policy`, `near_k3_policy`, `near_k4_policy`
  - detail audit only: `k2_laplace_sign_ok`,
    `k2_yukawa_smooth_ref_rel`, `k2_orientation_pair_ref_rel`

Acceptance:

- `test/pbe_born_ion_matrix.ps1 -SkipGrid` passes.
- `duffy_adaptive` is never worse than `analytic_singular` on the Born matrix.
- Off-center single-charge sphere test passes after it is added.
- Multi-charge Kirkwood sphere test passes after it is added.
- `test/pbe_tabi_near_kernel_blocks_smoke.ps1` passes and confirms that
  canonical one-ring pairs create per-block K1-K4 decisions.
- Direct and treecode leaf near-field report the same near integration policy.

Default policy:

- Keep `-pbe-tabi-near-integration=auto` mapped to `centroid` until the full
  ASI/Duffy path passes all analytic gates.
- After it passes, consider switching `auto` to:
  - self panels: `analytic_singular`.
  - one-ring panels: `duffy_adaptive`.
  - far non-neighbor panels: centroid or Gaussian triangle quadrature.

## Stage 2 - Analytic Validation Beyond Born Ion

Goal: make analytic correctness stronger than a centered Born ion.

Tasks:

- Add off-center charge Kirkwood sphere fixture.
- Add multi-charge Kirkwood sphere fixture.
- Add charge sign and charge magnitude symmetry checks.
- Add salt-dependence checks for multiple radii and ionic strengths.
- Keep salted off-center/multi-charge rows report-only until the centered
  screened-sphere oracle is promoted from sanity check to enforced gate.
- Add convergence slope reports for surface density and near-integration mode.
- Add `density_status` for off-center rows so remaining salted rows can be
  separated into `pass_threshold`, `density_refinement_needed`,
  `above_threshold`, or `oracle_not_converged`.

Acceptance:

- Centered Born ion remains the primary fast gate.
- Off-center and multi-charge fixtures pass within documented tolerances.
- `reaction_energy_source` must be `tabi_boundary_integral`.
- Any fallback or debug projection fails acceptance tests.

Test targets:

- `test/pbe_born_ion_matrix.ps1`
- `test/pbe_tabi_near_kernel_blocks_smoke.ps1`
- `test/pbe_tabi_near_block_rank.ps1`
- `test/pbe_kirkwood_offcenter_matrix.ps1`
- `test/pbe_kirkwood_multicharge_matrix.ps1`

## Stage 3 - Paper 2013 Surface And Protein Fixture Closure

Goal: make comparisons with the 2013 TABI paper meaningful rather than only
algorithmic.

Tasks:

- Add a stable clean-room 1A63 PQR fixture under `test/pbe_inputs_local` when
  available.
- Add optional report-only mesh fixture support if a paper-like mesh is
  available and provenance is clear.
- Improve TITAN SES mesh generation:
  - solvent-excluded surface mode.
  - deterministic vertex welding.
  - triangle quality checks.
  - nonmanifold edge audit.
  - patch valence diagnostics.
- Report surface provenance:
  - `surface_generator`
  - `discretization`
  - `mesh_vertices`
  - `mesh_triangles`
  - `boundary_nodes`
  - `nonmanifold_edges`
  - `solid_angle_fallbacks`

Acceptance:

- `test/pbe_tabi_jcp2013_paper_compare.ps1` passes.
- `test/pbe_tabi_1a63_paper_smoke.ps1` passes when fixture exists, otherwise
  skips with a clear reason.
- Exact paper values are report-only unless the same PQR, dielectric settings,
  salt, and surface mesh are used.

## Stage 4 - Production JCP Treecode

Goal: make `jcp2013-treecode` a trustworthy acceleration path.

Tasks:

- Ensure all 4 coupled BIE operator blocks use the final ASI near-field path.
- Reuse interaction lists across GMRES matvecs when geometry is unchanged.
- Keep analytic Taylor derivatives as the default.
- Keep finite-difference derivatives only as debug/reference.
- Add error gates:
  - fixed-density matvec relative error.
  - full energy relative error.
  - GMRES residual behavior.
- Tune tree parameters:
  - opening angle.
  - leaf size.
  - Taylor order 3 default.
  - Taylor order 4 debug/check path.
- Preserve deterministic results for `nt=1` and parallel runs.

Acceptance:

- Direct-vs-treecode matvec relative error <= `1e-4`.
- Direct-vs-treecode energy relative error <= `5e-4`.
- Order 4 is not worse than order 3.
- Treecode is faster than direct on medium protein fixtures.
- `-pbe-tabi-sum=auto` may choose treecode only after these gates pass.

Test targets:

- `test/pbe_tabi_treecode_accuracy.ps1`
- `test/pbe_tabi_jcp2013_treecode_full_smoke.ps1`
- `test/pbe_tabi_vs_grid_compare.ps1 -Mode perf`

## Stage 5 - Solver Stability And Preconditioning

Goal: reduce GMRES iterations and improve robustness on proteins.

Tasks:

- Keep restarted GMRES with Givens rotations as the default solver.
- Add stronger residual and stagnation diagnostics.
- Implement true block-Jacobi preconditioner:
  - leaf block or near-patch block grouping.
  - small dense factorization.
  - fallback to diagonal on singular blocks.
- Add preconditioner timing and memory estimates.
- Add acceptance gates on:
  - iteration count reduction.
  - no energy drift versus no-preconditioner direct.
  - no NaN/Inf in density, potential, or energy.

Acceptance:

- `-pbe-tabi-precond=auto` picks diagonal for small cases.
- `-pbe-tabi-precond=auto` picks block only when it reduces iterations on
  medium fixtures.
- Block preconditioner never changes accepted energy beyond tolerance.

## Stage 6 - TABI-PB 2.0 Node-Patch Mode

Goal: turn `tabi2-nodepatch` from report-only into a real experimental solver.

Tasks:

- Implement node-patch discretization:
  - node area from one-ring patch.
  - node normal from patch geometry.
  - local patch reconstruction for low-quality triangles.
- Implement TABI-PB 2.0-style block preconditioning.
- Add BLDTT/FMM-style evaluator as clean-room CPU implementation.
- Compare node-patch energy to centroid direct on the same geometry.
- Track node-patch counts and patch quality diagnostics.

Acceptance:

- `test/pbe_tabi2_nodepatch_smoke.ps1` no longer expects pending block status.
- Node-patch direct agrees with centroid direct within documented tolerance on
  Born/Kirkwood fixtures.
- Block preconditioner reduces iterations on mini-protein and 1A63 when
  fixture exists.
- BLDTT/FMM-style evaluator agrees with node-patch direct within tolerance.

## Stage 7 - Grid, Protein, And Experimental Comparisons

Goal: provide honest comparison reports without making grid parity the TABI
correctness target.

Tasks:

- Keep `TABI vs analytic` as the primary correctness gate.
- Keep `grid vs analytic` as a convergence diagnostic.
- Keep `TABI vs grid` as a secondary report for:
  - energy difference.
  - time.
  - RAM.
  - panels/grid points.
  - GMRES iterations.
  - surface/discretization provenance.
- Expand fixture sets:
  - sphere.
  - two-atoms.
  - mini-protein.
  - charged peptide.
  - compact mini-protein.
  - barnase.
  - 1A63, if fixture exists.
- Keep pKa experimental benchmark report-only until titration-state workflow is
  implemented.

Acceptance:

- `test/pbe_tabi_vs_grid_compare.ps1 -ProfileMemory` writes usable CSVs.
- Barnase compare runs with `-FixtureSet barnase` when local fixture cache is
  present.
- Reports explain remaining differences as surface/discretization/model
  differences instead of forcing TABI to match grid bitwise.

## Stage 8 - Performance And Memory Optimization

Goal: make TABI competitive for protein use without hiding accuracy failures.

Stage 2H/8A implementation status:

- Hidden treecode tuning knobs are available for performance sweeps:
  `-pbe-tabi-tree-leaf-size=<8..512>` and
  `-pbe-tabi-opening-theta=<0.05..2.0>`. Defaults remain leaf size `32`
  and opening theta `0.45`.
- `test/pbe_tabi_perf_matrix.ps1` writes timing/RAM/interaction-ratio CSVs
  for grid, TABI direct, and TABI treecode. It records derived
  `singular_near_pairs`, `leaf_near_pairs`, per-matvec near/far counts,
  seconds per million near interactions, ratios vs grid, and
  `performance_classification`.
- Current small-fixture smoke still classifies mini-protein treecode as
  `operator_limited`: leaf-near interactions dominate, so speed work should
  target interaction-list pruning and leaf-near kernels before more GMRES
  tuning.
- A mini-protein tuning smoke over leaf size `32/64/128` and opening theta
  `0.6/0.8/1.0` did not remove the bottleneck: best observed treecode row was
  still about `1.24s` at `24` surface points, with about `2.09M` near
  interactions per matvec and almost no far-node work. This confirms the next
  optimization should change the leaf-near algorithm, not just constants.
- Public defaults are unchanged; these knobs remain hidden/debug and are used
  by perf harnesses only.

Stage 2I/8B implementation status:

- JCP treecode now builds deterministic interaction lists once per solve and
  reuses them across all GMRES matvecs:
  `singular_near`, `leaf_near`, and `far_nodes`.
- `singular_near` keeps self/edge/vertex/one-ring panels on the existing
  ASI/Duffy/centroid path. Non-singular `leaf_near` uses a dedicated centroid
  JCP hot loop, so treecode leaf evaluation no longer pays near-patch
  classification cost per matvec.
- Timing/profile CSVs now report:
  `tabi_singular_near_eval_seconds`, `tabi_leaf_near_eval_seconds`,
  `tabi_singular_near_interactions`, `tabi_leaf_near_interactions`, and
  `tabi_far_node_interactions`.
- `test/pbe_tabi_treecode_interaction_list_smoke.ps1` verifies list validity,
  deterministic `nt=1` vs parallel output, interaction-build reuse, and
  singular/leaf/far count reporting.
- Mini-protein `24` panels/atom tuning improved best treecode runtime from
  about `1.24s` to about `0.47s` (`leaf=128`, `theta=0.8`) while keeping the
  classification `operator_limited`. The remaining bottleneck is real
  leaf-near work: about `2.09M` leaf-near interactions per matvec.
- A representative barnase treecode run (`24` panels/atom, `leaf=128`,
  `theta=0.8`) still timed out at `15` minutes. Grid barnase completed in
  about `0.83s` in that harness, so the next protein-scale slice must reduce
  leaf-near interaction count and memory traffic, not only scalar math cost.
- Analytic gates still pass: Born ion, off-center Kirkwood, multi-charge
  Kirkwood, and direct-vs-treecode accuracy. Grid/PBSA/SASA and clean-room
  regressions also pass.

Stage 2J/8C implementation status:

- Added hidden policy `-pbe-tabi-leaf-near-policy=auto|exact|smooth-node`.
  `auto` still resolves to `exact`, so production/default behavior is
  unchanged.
- JCP treecode interaction lists now split leaf work into
  `singular_near`, `leaf_panel_pairs`, `smooth_leaf_nodes`, and `far_nodes`.
  Smooth leaf nodes are accepted only for non-singular leaves that pass a
  conservative opening test; self/edge/vertex/one-ring pairs continue to use
  the shared ASI/Duffy/centroid near path.
- Timing/profile CSVs now report
  `tabi_leaf_panel_interactions`,
  `tabi_smooth_leaf_node_interactions`,
  `tabi_leaf_pruned_interactions`, and `tabi_leaf_near_policy`. Memory
  estimates include the split interaction-list entries.
- Added `test/pbe_tabi_leaf_near_policy_smoke.ps1` to verify `exact` vs
  `smooth-node` matvec agreement, deterministic `nt=1` vs parallel, and count
  sanity. `test/pbe_tabi_perf_matrix.ps1` can sweep leaf-near policy during
  tree tuning.
- Mini-protein `24` panels/atom tuning now shows best `exact` treecode around
  `0.37s` (`leaf=64`, `theta=0.6`), improving over the previous `~0.47s`.
  However, `smooth-node` is currently slower despite pruning some leaf panels:
  the Taylor node evaluator is more expensive than the panel loop for these
  small leaf buckets. The correct classification remains `operator_limited`.
- A representative barnase Stage 2J run
  (`24` panels/atom, `leaf=128`, `theta=0.8`, `smooth-node`) still timed out
  at `15` minutes. The same harness completed grid barnase in about `0.76s`
  with peak working set about `50MB`. This reinforces that the next slice must
  optimize leaf-panel and smooth-node evaluation before spending more time on
  barnase sweeps.
- The latest focused mini-protein run showed `smooth-node` pruning only about
  `78k` to `466k` panel interactions per matvec depending on leaf/theta, while
  leaving `~1.2M` to `2.0M` leaf-panel interactions per matvec. This is useful
  diagnostic data, but not enough to change defaults.
- Analytic gates still pass after the split: Born ion, off-center Kirkwood,
  multi-charge Kirkwood, direct-vs-treecode accuracy, profile/memory smokes,
  grid/PBSA/SASA regressions, and clean-room audit.

Stage 2K/8D implementation status:

- Added a `TreePanelCache` SoA geometry cache for TABI treecode and a
  per-matvec unknown cache for leaf-panel evaluation. The cache is internal to
  TABI treecode and does not change production CLI/defaults.
- Added a leaf-panel kernel dispatcher with `scalar` and forced `avx2-gather`
  paths. `auto` currently resolves to `scalar` because the AVX2 gather
  implementation is accurate but slower on mini-protein; forced AVX2 remains
  available for smoke testing and future tile work.
- Timing/profile CSVs now report:
  `tabi_leaf_panel_eval_seconds`,
  `tabi_smooth_leaf_node_eval_seconds`,
  `tabi_leaf_panel_kernel`, and
  `tabi_leaf_panel_pairs_per_second`.
- `test/pbe_tabi_leaf_near_policy_smoke.ps1` now compares scalar vs AVX2
  leaf-panel output directly and gates relative matvec difference at
  `<=1e-12`.
- Mini-protein `24` panels/atom tuning improved best exact treecode from the
  previous `~0.371s` to `~0.324s` (`leaf=64`, `theta=0.6`, scalar cache path).
  It is still not below the aspirational `<0.25s` target.
- The perf matrix now distinguishes `leaf_kernel_memory_limited` from generic
  `operator_limited`. The best mini-protein row reported about `88.7M`
  leaf-panel pairs/s but still around `2.09M` leaf-panel pairs/matvec.
- `smooth-node` remains report-only: it prunes some panels but is consistently
  slower, and rows are classified as `smooth_node_eval_limited`.
- Analytic gates and regressions still pass: Born ion, off-center Kirkwood,
  multi-charge Kirkwood, direct-vs-treecode accuracy, profile/memory smokes,
  grid/PBSA/SASA regressions, and clean-room audit.

Tasks:

- Profile hot paths:
  - surface mesh build.
  - adjacency/solid-angle/near-patch build.
  - tree build.
  - moment build.
  - traversal.
  - near eval.
  - far eval.
  - GMRES orthogonalization.
  - preconditioner apply.
- Move hot arrays to aligned SoA/AoSoA layouts where useful.
- AVX2 optimize:
  - centroid direct source loop.
  - near Gaussian/Duffy smooth quadrature batches.
  - far Taylor node evaluation.
- Avoid nested `TitanExecutor` use.
- Keep memory estimates and observed working set in the same CSV.

Acceptance:

- No regression in Born/Kirkwood gates.
- No regression in direct-vs-treecode accuracy.
- Barnase or medium protein treecode path is significantly faster than direct.
- TABI memory remains below comparable grid runs for large enough fixtures.

## Stage 9 - Documentation And Release Readiness

Goal: make TABI behavior understandable for users and future agents.

Tasks:

- Keep `manual.txt`, `help.cpp`, and `AGENT_ONBOARDING.md` in sync.
- Document hidden/debug knobs separately from public CLI.
- Keep public default simple:
  - `-pbe=tabi`
  - linear PB only.
  - no nonlinear/RPB/PBSA/focus support until explicitly implemented.
- Maintain clear rejection messages for unsupported combinations.
- Keep clean-room audit passing.

Release criteria:

- Native build passes.
- Core TABI analytic gates pass.
- Grid/PBSA/SASA regressions pass.
- Clean-room audit passes.
- Timing/memory compare CSVs are generated and documented.

## Stage 2L/8E - Contiguous Leaf Tiles

Implemented:

- Added hidden `-pbe-tabi-leaf-panel-kernel=auto|scalar|avx2-gather|tile`.
  Public/default behavior remains unchanged: `auto` resolves to `scalar`.
- Added exact contiguous tile interaction lists for non-singular complete leaf
  nodes. Leaves with self/edge/vertex/one-ring exclusions keep the residual
  flat panel list so singular near-field still goes through the same ASI/CPV
  path as direct mode.
- Added ordered Morton SoA cache for panel geometry and per-matvec ordered
  unknown channels. The tile evaluator reads contiguous leaf ranges instead of
  gathering source panel indices.
- Added timing/profile fields:
  - `tabi_leaf_tile_interactions`
  - `tabi_leaf_tile_panel_interactions`
  - `tabi_leaf_residual_panel_interactions`
  - `tabi_leaf_tile_eval_seconds`
  - `tabi_leaf_residual_eval_seconds`
  - `tabi_leaf_tile_status`
- Added `test/pbe_tabi_leaf_tile_kernel_smoke.ps1`. Current smoke result:
  scalar flat vs tile exact relative matvec error `<=1e-12`, deterministic
  `nt=1` vs `nt=4`, and active tile counts.
- Updated `test/pbe_tabi_perf_matrix.ps1` so `-TuneTreecode` sweeps
  `leaf_panel_kernel=scalar|tile` and can classify `tile_fast`,
  `tile_math_limited`, or `tile_memory_limited`.

Still open:

- The tile kernel is exact and contiguous but still scalar inside the tile.
  It removes gather/random index overhead first; the next slice should add
  aligned AVX2 loads over tile ranges if perf data says tile is memory-bound.
- Mini-protein/barnase perf matrix must be rerun after native build to decide
  whether tile hits the `<0.25s` target or remains `tile_memory_limited`.
- `smooth-node` remains report-only until exact tile is clearly faster.

Measured after native build:

- Focused mini-protein `24` panels/atom tuning completed.
- Best exact tile row: `leaf=64`, `theta=0.6`, `0.273s`, active tile list,
  about `101.7M` leaf-panel pairs/s, classified `tile_memory_limited`.
- Best exact scalar comparison in the same sweep: `leaf=64`, `theta=0.6`,
  `0.287s`, about `102.4M` pairs/s. Tile reduced random flat residual work
  (`~23.0M` tile panels plus `~2.1M` residual panels total across matvecs) but
  did not yet reach the `<0.25s` target.
- `leaf=128`, `theta=0.6`, tile was classified `tile_fast` versus scalar
  (`0.274s` vs `0.290s`), but overall best time is still just above target.

## Recommended Next Slice

Stage 2M/8F has now landed:

- Ordered tile geometry and ordered per-matvec unknown arrays use TITAN's
  `AlignedDoubleVector` from `aligned_allocator.h`.
- Added debug kernel choices `tile-scalar` and `tile-avx2`; the existing
  `tile` mode runtime-dispatches to `tile-avx2` when AVX2 is available and
  otherwise falls back to `tile-scalar`.
- `tile-avx2` uses aligned `_mm256_load_pd` over ordered tile ranges, with
  scalar prologue/tail for non-4-lane boundaries and lane-wise `std::exp` to
  preserve accuracy.
- `test/pbe_tabi_leaf_tile_kernel_smoke.ps1` now compares `tile-scalar` vs
  `tile-avx2` at relative matvec error `<=1e-12`.
- `test/pbe_tabi_perf_matrix.ps1 -TuneTreecode` now sweeps
  `leaf_panel_kernel=scalar|tile-scalar|tile-avx2` and reports
  `tile_avx2_fast`, `tile_avx2_math_limited`, or
  `tile_avx2_memory_limited`.

Measured Stage 2M/8F mini-protein result:

- Best exact `tile-avx2`: `leaf=128`, `theta=0.8`, `0.200s`,
  `~153.9M` leaf-panel pairs/s, classified `tile_avx2_fast`.
- Best exact `tile-scalar`: `leaf=128`, `theta=0.8`, `0.351s`,
  `~80.5M` pairs/s.
- Best exact scalar flat comparison in the same sweep: `0.377s`,
  `~78.8M` pairs/s.
- The `<0.25s` mini-protein target is met.

Stage 2N/8G has now landed:

- Added `test/pbe_tabi_barnase_tile_avx2_perf.ps1`, a focused barnase wrapper
  for `leaf=128`, `theta=0.8`, `leaf_near_policy=exact`, and
  `leaf_panel_kernel=tile-avx2`.
- Added `-CaseTimeoutSeconds` to `test/pbe_tabi_perf_matrix.ps1`. A timeout
  now records a CSV row instead of hanging the whole gate. If no partial
  `pbe_timing.csv` exists, future timeout rows are classified as
  `timeout_no_profile`; if partial timing exists, the row is classified from the
  emitted counters.
- Added barnase-specific classification for completed tile-AVX2 runs:
  `barnase_tile_avx2_candidate`, `barnase_residual_panel_limited`,
  `barnase_tile_eval_limited`, or `barnase_tile_avx2_slow`.

Measured Stage 2N/8G barnase result:

- Grid barnase completed in `1.235s`, peak RAM `50,417,664` bytes, estimated
  bytes `58,280,744`, energy `-14893.83415 kJ/mol`.
- TABI treecode barnase with exact tile-AVX2 timed out at `900.09s`. It did not
  reach normal solver shutdown, so no TABI `pbe_timing.csv` or memory profile
  was emitted for that row.
- This means the current barnase bottleneck is not yet measurable from the final
  profile CSV; the next performance slice must add progress/partial timing
  inside the TABI solve before another long barnase run.

Measured medium fixture result with the same tile-AVX2 policy:

- `charged-peptide`: grid `0.00205s`; TABI treecode `0.804s`, peak RAM ratio
  `1.57x`, `1201` panels, `13` GMRES iterations, `86.6M` leaf-panel pairs/s,
  `15.1M` tile-panel interactions and `3.56M` residual-panel interactions.
  Classification: `leaf_kernel_memory_limited`.
- `compact-mini-protein`: grid `0.00232s`; TABI treecode `6.728s`, peak RAM
  ratio `2.15x`, `2241` panels, `38` GMRES iterations, `98.1M` leaf-panel
  pairs/s, `166.4M` tile-panel interactions and `19.3M` residual-panel
  interactions. Classification: `leaf_kernel_memory_limited`.
- Tile-AVX2 still helps mini-protein, but medium fixtures show the exact
  leaf-panel path is memory/work-count limited before barnase scale.

Stage 2O/8H has now landed:

- Added hidden `-pbe-tabi-leaf-near-policy=large-protein`.
- `auto` remains unchanged. `large-protein` is opt-in and keeps
  self/edge/vertex/one-ring pairs on the exact singular near path.
- Added masked tile interactions: leaves with a small singular exclusion set
  now store a contiguous leaf tile plus an exclusion list instead of expanding
  all non-singular panels into the flat residual list.
- Added conservative smooth-node v2 inside `large-protein`: only non-singular
  leaves with `node.radius / distance < min(0.25, 0.35 * opening_theta)` are
  promoted to Taylor node evaluation.
- Added `test/pbe_tabi_large_protein_policy_smoke.ps1`. Current smoke result:
  masked tile reduced residual interactions from `57,660` to `5,032` and
  matched exact residual-flat output at relative matvec error `<=1e-12`.
- Added `test/pbe_tabi_large_protein_perf.ps1`, a focused wrapper for
  `surface_points=24`, `leaf=128`, `theta=0.8`,
  `leaf_panel_kernel=tile-avx2`, `leaf_near_policy=large-protein`, and
  `precond=block`.

Measured Stage 2O/8H medium result:

- `charged-peptide`: TABI large-protein `0.597s` vs previous exact tile-AVX2
  `~0.804s`; residual panels fell from `3.56M` to `3,528`; tile panels were
  `18.67M`; GMRES stayed `13` iterations. Classification:
  `large_protein_residual_pruned`.
- `compact-mini-protein`: TABI large-protein `7.920s` vs previous exact
  tile-AVX2 `~6.728s`; residual panels fell from `19.3M` to `0`, but tile
  panels rose to `185.39M`, smooth nodes were `2,680`, and GMRES stayed `38`
  iterations. Classification: `large_protein_residual_pruned`, but wall time
  regressed.
- Block preconditioner remains report-only in this slice:
  `-pbe-tabi-precond=block` still audits as the existing pending/diagonal path,
  so no GMRES iteration reduction was observed.
- Full Born matrix was attempted with `-SkipGrid` but did not finish within the
  300 second turn budget; no correctness tolerance was changed.

The next high-value slice is Stage 2P/8I:

1. Keep `large-protein` opt-in; do not promote it to `auto`.
2. Add an error/perf gate comparing `exact` vs `large-protein` on the same
   medium fixture rows before running barnase again.
3. Add a policy to disable smooth-node v2 when tile-panel work grows more than
   residual work shrinks; compact-mini-protein currently shows that failure
   mode.
4. Implement a real micro-block preconditioner or explicitly remove `block`
   from perf wrappers until it reduces GMRES iterations.
5. If barnase is rerun before partial progress profiling exists, use
   `large-protein` only after charged/compact both improve wall time.

Stage 2P/8I has now landed:

- Replaced the `large-protein` masked tile runtime exclusion scan with
  `leaf_tile_span` lists. Each target leaf now stores up to four contiguous
  ordered-panel spans with singular panels removed ahead of time. Residual flat
  panels are kept only when spans are too fragmented or too small.
- Disabled smooth-node v2 inside `large-protein`; `smooth-node` remains a
  separate hidden/debug policy.
- Added a real callable GMRES preconditioner path. `-pbe-tabi-precond=block`
  now builds micro-block Jacobi factors in Morton order with a cap of 32 panels
  per block. Invalid blocks fall back to diagonal, and the audit reports
  `preconditioner_effective=block` only when at least one LU block is valid.
- `test/pbe_tabi_large_protein_policy_smoke.ps1` now checks exact vs span-tile
  `large-protein` at relative matvec error `<=1e-12`, deterministic parallel
  behavior, and span/residual sanity. Current smoke result: `848` panels,
  residual interactions reduced from `57,660` to `7,217`, with `25,125` span
  panels.
- `test/pbe_tabi_large_protein_perf.ps1` now runs exact and `large-protein`
  rows side-by-side for medium fixtures, but for barnase it runs only grid plus
  optimized `large-protein` to avoid two 900 second TABI attempts.
- `test/pbe_tabi_perf_matrix.ps1` now records preconditioner requested/effective
  and classifies `large-protein` against the matching exact row as
  `large_protein_candidate`, `span_tile_regressed`,
  `span_tile_work_regressed`, `block_precond_regressed`, or
  `span_tile_not_enough`.

Measured Stage 2P/8I medium result with block preconditioner:

- `charged-peptide`: exact tile-AVX2 `0.713s`, `13` GMRES iterations,
  `18.67M` leaf-panel interactions, `3.56M` residual panels, peak RAM
  `15.2MB`. `large-protein` span-tile `0.722s`, `13` iterations, same total
  leaf-panel work, residual panels reduced to `0.954M`, peak RAM `12.8MB`.
  Classification: `large_protein_candidate` in the latest run, but repeated
  runs fluctuate around exact because the wall-time delta is small.
- `compact-mini-protein`: exact tile-AVX2 `7.739s`, `35` GMRES iterations,
  `171.80M` leaf-panel interactions, `17.84M` residual panels, peak RAM
  `22.3MB`. `large-protein` span-tile `7.933s`, `35` iterations, same total
  leaf-panel work, residual panels reduced to `5.01M`, peak RAM `16.8MB`.
  Classification: `span_tile_regressed`.
- Preconditioner comparison: on `compact-mini-protein`, diagonal preconditioning
  needed `38` GMRES iterations and `6.954s` exact tile time; block
  preconditioning reduced this to `35` iterations and `5.528s`. On
  `charged-peptide`, block did not reduce iterations (`13` either way).

Stage 2P decision:

- Barnase long TABI was not rerun in this slice because the medium gate did not
  pass for both fixtures. This follows the Stage 2P rule to avoid another
  900 second barnase attempt until the optimized policy is no slower than exact
  on `charged-peptide` and `compact-mini-protein`.
- A fresh barnase grid-only baseline was run for comparison: grid completed in
  `0.995s`, peak RAM `50.4MB`, estimated bytes `58.3MB`, `145` iterations,
  and energy `-14893.83415 kJ/mol`.
- The current optimization is useful for RAM and residual-list shrinkage, but
  not yet a wall-time win because total leaf-panel work is unchanged. The next
  slice should reduce evaluated leaf work, not only move residual flat panels
  into spans.

Recommended Stage 2Q/8J:

1. Keep exact tile-AVX2 as the performance baseline; keep `large-protein`
   hidden/report-only.
2. Add an adaptive span gate: use span tile only when measured residual-flat
   throughput is worse than tile throughput for that target/node bucket, or when
   RAM pressure is the explicit goal.
3. Revisit smooth-node with a stricter local error estimator that reduces total
   leaf-panel work instead of merely shifting residual panels into tile panels.
4. Extend block preconditioner grouping to real tree leaves/one-ring patches and
   accept it only where it reduces GMRES iterations by at least 20%.
5. Run barnase grid plus optimized TABI only after medium fixtures show no wall
   time regression; until then, keep the existing Stage 2N grid vs timeout
   barnase result as the current large-protein benchmark.

Stage 2Q/8J has now landed:

- Added hidden `-pbe-tabi-leaf-near-policy=adaptive-large`; `auto`,
  `exact`, `large-protein`, and `smooth-node` defaults remain unchanged.
- `adaptive-large` keeps singular self/edge/vertex/one-ring pairs exact, uses a
  stricter span-tile gate for singular leaf exclusions, and can promote
  non-singular leaves to smooth-node Taylor only when a conservative geometric
  p4 proxy passes.
- Added adaptive counters to timing/profile CSVs:
  `adaptive_smooth_leaf_accepted/rejected`,
  `adaptive_span_tile_accepted/rejected`,
  `adaptive_exact_leaf_fallbacks`, reject-reason counters, and
  `leaf_work_reduction_ratio`.
- Added `test/pbe_tabi_adaptive_leaf_policy_smoke.ps1`. Current smoke result:
  `1058` panels, `span_reject=34690`, `smooth_accept=492`,
  `smooth_reject=16090`, and `leaf_work_reduction_ratio=0.0240`.
- Block preconditioner grouping now follows true tree leaves and chunks leaf
  ranges at 32 panels/block; `preconditioner_effective=block` requires valid LU
  blocks to cover at least 70% of panels.
- `test/pbe_tabi_large_protein_perf.ps1` now compares `exact`,
  `large-protein`, and `adaptive-large` on medium fixtures, while barnase mode
  runs only optimized TABI to avoid duplicate long attempts.

Measured Stage 2Q/8J medium result:

- `charged-peptide`: exact tile-AVX2 `0.531s`, `13` GMRES iterations,
  `18.67M` leaf-panel interactions, peak RAM `15.3MB`. `adaptive-large`
  `0.580s`, `13` iterations, `18.67M` leaf-panel interactions, no pruned leaf
  work, peak RAM `13.2MB`. Classification: `adaptive_error_limited`.
- `compact-mini-protein`: exact tile-AVX2 `7.524s`, `35` GMRES iterations,
  `171.80M` leaf-panel interactions, peak RAM `22.4MB`. `adaptive-large`
  `7.119s`, `35` iterations, `171.64M` leaf-panel interactions,
  `165,760` pruned panels, `leaf_work_reduction_ratio=0.000965`, peak RAM
  `18.7MB`. Classification: `adaptive_error_limited`.
- `large-protein` remains the best hidden medium policy for
  `compact-mini-protein` in this run (`6.755s`, `17.0MB`), but it still does
  not reduce total leaf-panel work.

Stage 2Q decision:

- Barnase TABI was not rerun because `adaptive-large` did not meet the medium
  gate. The current barnase comparison remains: grid completes in about `1s`
  and the previous TABI exact tile-AVX2 run timed out at `900s`.
- The adaptive estimator is now instrumented but too conservative to affect
  protein-scale work. The next performance slice should target a stronger
  mathematically bounded approximation or a different solver-level reduction.

Recommended Stage 2R/8K:

1. Keep `adaptive-large` hidden/report-only; do not promote it.
2. Add a target-cluster or leaf-pair batched Taylor evaluator so one accepted
   smooth leaf removes many panel interactions at once with lower overhead.
3. Add a configurable report-only sweep for the adaptive p4 proxy threshold
   (`strict`, `balanced`, `aggressive`) to find the accuracy/perf knee on
   medium fixtures before touching barnase.
4. Improve block preconditioner strength: include cross-chunk one-ring couplings
   or switch to a small Schwarz-style additive preconditioner; accept only if
   GMRES iterations fall by at least 20%.
5. Consider a production-safe early barnase progress snapshot only after medium
   leaf work reduction exceeds 25%.

Stage 2R/8K has now landed:

- Added hidden `-pbe-tabi-adaptive-level=strict|balanced|aggressive`. The
  public/default path is unchanged; `auto` still stays on the safe/exact path.
- Added target-cluster counters to timing/profile CSVs:
  `target_cluster_interactions`, `target_cluster_targets`,
  `target_cluster_panels_pruned`, `target_cluster_rejected`, and
  `target_cluster_eval_seconds`.
- `adaptive-large` now uses the level to relax or tighten the existing p4
  proxy for non-singular smooth leaf pruning. This is still report-only.
- Added `test/pbe_tabi_target_cluster_policy_smoke.ps1`. Current smoke result:
  `192` synthetic clustered panels, `384` target-cluster interactions,
  `24,576` pruned panels, `0` rejects, and relative matvec error
  `1.75e-15`.
- `test/pbe_tabi_large_protein_perf.ps1` now sweeps `strict`, `balanced`, and
  `aggressive` on medium fixtures. Barnase mode uses a single aggressive
  adaptive probe to avoid multiplying timeout risk.
- `test/pbe_tabi_perf_matrix.ps1` now classifies adaptive rows as
  `cluster_leaf_pruned`, `adaptive_level_too_aggressive`,
  `target_cluster_error_limited`, or the existing exact/tile classes.

Measured Stage 2R/8K medium result:

- `charged-peptide`: exact tile-AVX2 strict `0.561s`, `13` GMRES iterations,
  `18.67M` leaf-panel interactions, peak RAM `14.6MB`. `adaptive-large`
  strict stayed exact (`0` pruned panels) and slowed to `0.855s`.
  `adaptive-large` balanced pruned `0.532M` panels
  (`leaf_work_reduction_ratio=0.0285`) but slowed to `1.238s`.
  `adaptive-large` aggressive pruned `1.976M` panels
  (`leaf_work_reduction_ratio=0.1058`) but slowed to `2.056s`.
  Classification: strict `adaptive_error_limited`, balanced/aggressive
  `adaptive_level_too_aggressive`.
- `compact-mini-protein`: exact tile-AVX2 strict `8.523s`, `35` iterations,
  `171.80M` leaf-panel interactions, peak RAM `21.3MB`. `large-protein`
  strict is still the best hidden medium row in this run at `7.347s` and
  `16.4MB`, but it does not reduce total leaf work. `adaptive-large` strict
  pruned only `165,760` panels and took `8.268s`. Balanced pruned `8.335M`
  panels but took `13.919s`; aggressive pruned `28.83M` panels but took
  `29.278s`. Classifications: strict `target_cluster_error_limited`,
  balanced/aggressive `adaptive_level_too_aggressive`.

Measured Stage 2R/8K barnase probe:

- Barnase grid completed in `1.769s`, peak RAM `48.4MB`, `145` grid
  iterations, energy `-14893.83415 kJ/mol`.
- Barnase TABI `adaptive-large/aggressive`, `leaf=128`, `theta=0.8`,
  `tile-avx2`, `precond=block`, `surface_points=24` timed out after
  `180.37s` before producing a timing row. Classification:
  `barnase_probe_timeout`.

Stage 2R decision:

- Do not promote `adaptive-large` or `adaptive-level` tuning. The smoother
  target-cluster path reduces mathematical work on medium fixtures, but its
  evaluator overhead dominates wall time.
- Barnase full `900s` should still not be run until medium fixtures beat exact
  wall time and reduce total leaf-panel work meaningfully.
- The current bottleneck is no longer just residual list construction; it is
  the cost of per-target Taylor/smooth-node evaluation and unchanged GMRES
  iteration count.

Recommended Stage 2S/8L:

1. Keep exact tile-AVX2 as the baseline and keep `large-protein` as the only
   hidden policy worth comparing on compact fixtures.
2. Replace per-target smooth-node evaluation with true target-leaf/source-leaf
   pair batching: one leaf-pair accept should evaluate all target panels in a
   contiguous batch using shared node moments and vectorized target loops.
3. Add a cheap pre-estimator before entering the smooth evaluator so rejected
   leaves do not pay setup cost.
4. Strengthen solver reduction instead of tree-only work reduction: try
   restarted GMRES parameter sweep and a Schwarz/near-patch preconditioner that
   demonstrably reduces iterations on `charged-peptide` and
   `compact-mini-protein`.
5. Add a lightweight progress snapshot only for long TABI runs, so barnase
   timeouts can report panels, GMRES iteration, and interaction counts before
   being killed.

Stage 2S/8L has now landed:

- `test/pbe_tabi_perf_matrix.ps1` now supports `-Threads` instead of hard-coded
  `-nt=1`, and emits `threads`, `speedup_vs_nt1`,
  `parallel_efficiency`, and `serial_bottleneck_class`.
- Added `test/pbe_tabi_thread_scaling.ps1` for focused scaling runs.
- Added hidden parser support for
  `-pbe-tabi-cluster-policy=off|leaf-pair-strict|leaf-pair-balanced`; it is a
  scaffold only in this slice and remains `off` in all default/perf rows.
- Tree moments are now built with a level-parallel bottom-up pass. The JCP
  treecode path also computes the four required moment channels
  (`dphi/dn`, `phi*nx`, `phi*ny`, `phi*nz`) in one packed traversal instead of
  four separate serial traversals.
- GMRES vector operations now use `TitanExecutor` for large vectors:
  dot/norm, basis normalization, orthogonalization updates, residual updates,
  and solution updates. The threshold is intentionally high to avoid slowing
  small/medium vectors.
- Block preconditioner setup now factors independent micro-blocks in parallel,
  and block preconditioner apply runs independent blocks in parallel when
  active.
- Added `test/pbe_tabi_parallel_moments_smoke.ps1`. Current result:
  charged-peptide exact treecode `nt=1` and `nt=4` produce identical energy
  (`delta=0`).

Measured Stage 2S/8L thread scaling result:

- `charged-peptide`, exact tile-AVX2: `nt=1` `0.907s`, `nt=4` `0.443s`
  (`2.05x`), `nt=8` `0.282s` (`3.22x`). Energy unchanged:
  `-801.2530074 kJ/mol`.
- `charged-peptide`, `large-protein`: `nt=8` `0.303s`, lower RAM than exact
  but slower than exact at this thread count.
- `compact-mini-protein`, exact tile-AVX2: `nt=1` `9.972s`, `nt=4`
  `3.600s` (`2.77x`), `nt=8` `1.931s` (`5.16x`). Energy unchanged:
  `-1294.750935 kJ/mol`.
- `compact-mini-protein`, `large-protein`: `nt=4` `3.505s`, `nt=8`
  `2.240s`; useful for RAM (`~17.7MB` vs exact `~23.2MB`) but no longer the
  fastest wall-time row.

Measured Stage 2S barnase result:

- Barnase grid with `nt=8`: `0.939s`, peak RAM `50.5MB`, energy
  `-14893.83415 kJ/mol`.
- Barnase TABI exact tile-AVX2 with `nt=8`, `leaf=128`, `theta=0.8`,
  `precond=block`, and `surface_points=24`: timeout after `900.13s` before
  producing a timing row. Classification remains `timeout_no_profile`.

Stage 2S decision:

- Parallelization is now materially better for medium protein fixtures. The
  `compact-mini-protein` target is met (`<4s` at `nt=4`, `<2s` at `nt=8`), and
  `charged-peptide` meets the `<0.35s` target at `nt=8`.
- Barnase still needs an algorithmic reduction, not just more threads. The
  missing piece is true leaf-pair/source-leaf batching or a stronger solver
  preconditioner that reduces GMRES work.
- Do not promote `large-protein`, `adaptive-large`, or cluster-policy paths.
  Exact tile-AVX2 is currently the fastest hidden TABI performance baseline.

Recommended Stage 2T/8M:

1. Implement real leaf-pair batching rather than the current cluster-policy
   scaffold: source leaf accepted once, target leaf evaluated as a contiguous
   batch, shared Taylor setup amortized across all target panels.
2. Add a lightweight long-run progress snapshot before another barnase attempt:
   write panels, current GMRES iteration, matvec count, and interaction counts
   periodically so timeout rows are diagnosable.
3. Add a GMRES/preconditioner sweep for restart size and Schwarz/near-patch
   block coverage; accept only if iterations drop on both medium fixtures.
4. Keep exact tile-AVX2 `nt=8` as the performance baseline for all future
   protein-large comparisons.

Stage 2T/8M has now landed:

- Added hidden TABI progress controls:
  `-pbe-tabi-progress=off|auto|detail` and
  `-pbe-tabi-progress-interval=<1..300>`. `auto` is enabled only when
  `-pbe-timing=detail` reaches the TABI backend.
- TABI now writes and flushes `pbe_tabi_longrun_progress.csv` during long
  solves. Rows are emitted for surface, tree build, preconditioner, GMRES
  start, timed GMRES iteration snapshots, restart boundaries, and solve
  completion.
- Added `test/pbe_tabi_longrun_progress_smoke.ps1` and
  `test/pbe_tabi_timeout_progress_smoke.ps1`. Both pass; timeout runs now
  preserve a readable progress CSV.
- Added `test/pbe_tabi_barnase_timing_profile.ps1`, which writes
  `pbe_tabi_barnase_timing_profile.csv` by combining the grid/TABI perf matrix
  with the last TABI progress row when TABI times out.

Measured Stage 2T Barnase profile:

- Barnase grid with `nt=8`: `0.928s`, peak RAM `50.5MB`, energy
  `-14893.83415 kJ/mol`.
- Barnase TABI exact tile-AVX2 with `nt=8`, `leaf=128`, `theta=0.8`,
  `precond=block`, `surface_points=24`: timeout after `900.11s`.
- The last progress row was at `893.59s`, phase `gmres_iter`, GMRES iteration
  `49`, matvec calls `50`, residual `2.0923e-03`, panels `39011`, unknowns
  `78022`, tree nodes `350`, RSS `390.6MB`.
- Cumulative work by the last row: `42.18B` leaf-panel interactions,
  `41.70B` tile-panel interactions, `473.8M` residual-panel interactions,
  `32.0M` far-node interactions, and `17.5M` singular-near interactions.
  Classification: `matvec_limited`.

Stage 2T decision:

- Barnase no longer needs timing infrastructure work before optimization; the
  bottleneck is now measured. Surface (`0.086s`), tree build (`~3.0s`), and
  preconditioner setup (`0.024s`) are not the blockers.
- The next performance slice must reduce per-matvec leaf-panel work by orders
  of magnitude or reduce GMRES iterations substantially. More outer-thread
  scaling alone will not solve Barnase.

Recommended Stage 2U/8N:

1. Implement true leaf-pair batching for JCP treecode: accept source-leaf /
   target-leaf interactions once, evaluate a contiguous target leaf batch, and
   amortize Taylor setup across target panels.
2. Add a near/far acceptance sweep that targets Barnase leaf-panel work, with
   medium fixtures as accuracy gates and Barnase as report/profile.
3. Build a stronger Schwarz/near-patch preconditioner and accept only if GMRES
   iterations drop on both medium fixtures and Barnase progress rows show a
   lower residual per matvec.
4. Keep `pbe_tabi_longrun_progress.csv` enabled for all Barnase profile runs
   until TABI finishes under the timeout.

Stage 2U/8N has now landed:

- Hidden TABI surface density controls now allow
  `-pbe-tabi-surface-points=8..4096`; normal defaults remain unchanged.
- Timing/profile rows now report
  `tabi_surface_point_count_requested`,
  `tabi_panel_density_per_A2`, and `tabi_panel_budget_status`.
- Added `test/pbe_tabi_low_density_surface_smoke.ps1`, covering
  `surface_points=8,12,16,24`, deterministic `nt=1` vs `nt=4`, finite panel
  geometry, and long-run progress output.
- Added `test/pbe_tabi_barnase_surface_budget_probe.ps1`, which runs the grid
  reference plus TABI treecode surface/tree sweeps and summarizes progress rows
  even when TABI times out or finishes with non-converged status.
- `test/pbe_tabi_perf_matrix.ps1` now has report-only
  `-AllowFailedStatus` for performance probes; regular tests still fail on
  non-`ok` energy status by default.

Measured Stage 2U Barnase focused sweep:

- Grid reference with `nt=8`: `0.936s`, peak RAM `50.5MB`, energy
  `-14893.83415 kJ/mol`.
- Focused `90s` sweep used `surface_points=8,12,16`, `leaf=32,64`,
  `theta=0.8`, `tile-avx2`, `precond=block`, and `nt=8`.
- `surface_points=8` and `12` both produced `15129` panels and `30258`
  unknowns; this is effectively at the requested protein panel budget
  (`panel_density_per_A2=1.5067`, `panel_budget_status=target_density` on a
  completed row).
- Best progress row in the short sweep was `surface_points=8`, `leaf=64`:
  `24` matvec calls in `90.16s`, residual `2.4068e-03`, RSS `151.9MB`,
  and `144.1M` leaf-panel interactions per matvec. This is a `5.85x`
  reduction from the Stage 2T baseline of `843.5M` per matvec.
- `surface_points=16` produced `25672` panels and remained too dense for the
  immediate Barnase target; short-sweep work was `215M-327M` leaf-panel
  interactions per matvec depending on leaf size.

Measured Stage 2U best-candidate rerun:

- Candidate: `surface_points=8`, `leaf=64`, `theta=0.8`,
  `leaf_near_policy=exact`, `leaf_panel_kernel=tile-avx2`,
  `precond=block`, `nt=8`.
- TABI finished before the `900s` timeout in `363.85s`, with finite energy
  `-12558.31419 kJ/mol`, `15129` panels, `30258` unknowns, and final RSS
  `123.9MB` from progress (`269.7MB` peak working set in final memory CSV).
- The run finished with `status=failed` because GMRES stopped at
  residual `7.0733e-03` after `96` iterations / `99` matvec calls. The best
  pre-restart residual reached `1.3061e-03`, then the final restart row reset
  the reported residual upward.
- Cumulative work was `14.27B` leaf-panel interactions, or `144.1M` per
  matvec; tile work still dominates (`14.08B` tile-panel interactions).
  Classification: `barnase_tile_eval_limited` plus
  `leaf_work_reduced_5x`.

Stage 2U decision:

- The panel-budget comment was directionally correct: lowering the TABI
  surface density is the first large Barnase win. It reduced panels by
  `2.58x` (`39011 -> 15129`), unknowns by `2.58x`, RSS by roughly `3x` in
  progress rows, and per-matvec leaf-panel work by `5.85x`.
- This alone is not enough for production Barnase: the optimized coarse run
  finishes under timeout but is still hundreds of seconds and fails the current
  GMRES convergence gate.
- `leaf=64` is better than `leaf=32` for this coarse Barnase setup despite
  more leaf work per matvec, because it advanced more GMRES iterations within
  the same wall time.
- `precond=block` is active but not strong enough. It does not prevent
  `96` GMRES iterations or the final non-converged status on Barnase.

Recommended Stage 2V/8O:

1. Add a Barnase convergence-focused sweep over GMRES tolerance/restart and
   `precond=none|diagonal|block` at the coarse `surface_points=8` budget.
   Record residual-before-restart and final residual separately so restart
   bookkeeping cannot hide the best achieved residual.
2. Make the panel-budget mode explicit for protein perf tests:
   choose `surface_points=8` as the hidden Barnase coarse candidate, and keep
   `surface_points=24+` for analytic/sphere correctness gates.
3. Optimize the tile-dominated path next: tile exponent/math throughput,
   larger target chunks, and true leaf-pair/source-leaf batching. The current
   dominant count is still `14.08B` tile-panel interactions on the completed
   coarse Barnase run.
4. Keep grid comparison in the same CSV. Grid remains much faster for this
   fixture (`<1s`), while TABI's current advantage is lower asymptotic memory
   after surface density is controlled, not wall time.

Stage 2V/8O has now landed:

- Added hidden `-pbe-tabi-gmres-restart=<8..128>`. The production default
  remains `32`.
- TABI timing/profile/progress rows now distinguish requested/effective restart
  and best/final residual:
  `tabi_gmres_restart_requested`, `tabi_gmres_restart_effective`,
  `tabi_best_residual`, `tabi_best_residual_iter`, and
  `tabi_final_residual`.
- `pbe_tabi_longrun_progress.csv` now flushes the same best/final residual
  fields, so timeout rows no longer lose the best residual reached before a
  restart boundary.
- Direct TABI with `precond=block` now builds the JCP tree workspace needed for
  block grouping without using treecode matvec. This lets the truth matrix test
  direct+block as an actual preconditioner branch instead of silently behaving
  like diagonal-only.
- Added `test/pbe_tabi_gmres_restart_smoke.ps1`, covering restart `32/64/100`,
  progress/timing columns, and deterministic energy with `nt=1` vs `nt=4`.
- Added `test/pbe_tabi_barnase_truth_matrix.ps1`. It emits grid, direct, and
  treecode rows with timeout-safe progress summaries and classifies the run as
  `treecode_noise_limited`, `nearfield_or_bie_limited`, `restart_limited`,
  `block_precond_regressed`, `safe_treecode_candidate`, or `inconclusive`.

Measured Stage 2V quick Barnase truth sanity:

- This was a harness sanity run, not the full `900s` truth matrix:
  `-Quick -CaseTimeoutSeconds 60 -ProfileMemory`.
- Grid reference completed in `0.859s`, peak RAM `50.5MB`, energy
  `-14879.32072 kJ/mol`.
- Direct TABI, `surface_points=8`, `precond=diagonal`, `restart=32`,
  `nt=8`: timeout at `60.04s`, `15129` panels, `30258` unknowns, `14`
  matvec calls, best/final residual `4.3236e-03`, progress RSS `22.4MB`.
- Treecode TABI, `theta=0.5`, `leaf=32`, `precond=diagonal`, `restart=64`,
  `nt=8`: timeout at `60.05s`, `15129` panels, `30258` unknowns, `10`
  matvec calls, best/final residual `5.3049e-03`, cumulative
  `1.782B` leaf-panel interactions, about `178.2M` per matvec, progress RSS
  `272.1MB`.
- Quick classification: `nearfield_or_bie_limited`. This is not yet a full
  convergence verdict, but it is enough to show that a direct coarse Barnase
  truth row does not rapidly converge below `1e-4` in the first minute.

Verification after Stage 2V:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`.
- `test/pbe_tabi_gmres_restart_smoke.ps1`: PASS.
- `test/pbe_tabi_barnase_truth_matrix.ps1 -Quick -CaseTimeoutSeconds 60
  -ProfileMemory`: completed and wrote
  `test/tmp/pbe_tabi_barnase_truth_matrix/pbe_tabi_barnase_truth_matrix.csv`.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS
  (`panels=1422`, `rel3=9.38761e-07`, `rel4=1.35948e-07`).
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_tabi_profile_smoke.ps1`: PASS
  (`panels=5386`, `matvec_calls=30`).
- `test/pbe_clean_room_audit.ps1`: PASS.

Recommended Stage 2W/8P:

1. Run the full Barnase truth matrix only when a long run is acceptable:
   direct diagonal/block with restart `32/64/100`, then treecode
   `theta=0.5`, `leaf=32/64`, restart `64/100`. The quick run proves the
   harness works but does not replace the full branch decision.
2. If full direct still stagnates above `1e-3`, pause treecode optimization and
   strengthen the BIE near-field/preconditioner path: one-ring ASI/CPV,
   near-patch Schwarz or leaf-local block coupling, and preconditioner
   diagnostics versus diagonal.
3. If full direct converges but treecode does not, tune treecode accuracy before
   speed: smaller MAC, higher Taylor order, and direct-vs-treecode residual
   drift on the same Barnase coarse surface.
4. Keep `surface_points=8` as the hidden Barnase coarse performance candidate,
   but do not use it for analytic sphere gates.
5. Keep grid in every Barnase CSV for time/RAM/energy context; grid remains the
   diagnostic reference, not the TABI correctness oracle.

Stage 2W/8P has now landed:

- Added hidden geometric near-buffer control
  `-pbe-tabi-quad-buffer=auto|off|gauss7|adaptive`.
  Public/default behavior stays safe: `auto` resolves to `off`. Validation
  probes must pin `gauss7` or `adaptive` explicitly.
- TABI now classifies non-topological close triangle pairs between one-ring
  CPV and centroid. The classifier uses a cheap centroid/radius prefilter and
  then triangle-triangle distance; pairs with
  `triangle_gap <= 0.75 * max_edge_length` enter the buffer.
- Direct JCP and treecode leaf near-field share the same buffer-aware
  `accumulate_jcp_near_ordered(...)` path. Tree interaction build now keeps
  geometric-buffer pairs out of exact leaf tiles/residual centroid loops, so
  close cleft pairs cannot bypass quadrature in treecode mode.
- Buffer kernels:
  - `gauss7`: fixed 7-point triangle quadrature for close non-one-ring pairs.
  - `adaptive`: subdivided 7-point quadrature for the closest bucket, fixed
    7-point otherwise.
  - Self and one-ring ASI/CPV stay unchanged; block preconditioner uses the
    same near operator path.
- Timing/profile now includes:
  `tabi_quad_buffer`, `tabi_near_geometric_buffer_pairs`,
  `tabi_near_geometric_buffer_gauss_pairs`,
  `tabi_near_geometric_buffer_adaptive_pairs`,
  `tabi_near_geometric_buffer_centroid_pairs`,
  `tabi_near_geometric_buffer_min_gap`,
  `tabi_near_geometric_buffer_mean_gap_ratio`, and
  `tabi_near_geometric_buffer_max_edge_mean`.
- Added detailed audit file for timing detail runs:
  `pbe_tabi_geometric_near_audit.csv`, bucketed by gap ratio, center ratio,
  block label, and kernel. This file is also listed in `pbe_summary.dat` and
  the log audit.
- Added `test/pbe_tabi_geometric_near_buffer_smoke.ps1`. It verifies a
  non-one-ring pair of parallel triangles separated by `0.3 A` is detected,
  that `gauss7` and `adaptive` counters are recorded, and that far
  non-one-ring pairs still use centroid.
- Extended `test/pbe_tabi_treecode_accuracy.ps1` so a geometric-buffer pair is
  counted through the treecode near path instead of being swallowed by the tile
  kernel.

Measured Stage 2W Barnase quick probe:

- Command:
  `test/pbe_tabi_barnase_truth_matrix.ps1 -Quick -CaseTimeoutSeconds 60
  -ProfileMemory -NearIntegration duffy_adaptive -QuadBuffer adaptive`.
- Grid reference completed in `1.606s`, peak RAM `50.6MB`, energy
  `-14879.32072 kJ/mol`.
- Direct TABI, `surface_points=8`, `precond=diagonal`, `restart=32`,
  `near=duffy_adaptive`, `quad_buffer=adaptive`, `nt=8`:
  completed in `43.46s`, `15129` panels, `30258` unknowns, `6` matvec calls,
  best residual `8.3519e-05` at iteration `5`, final residual
  `2.2824e-03`, progress RSS `12.9MB`, peak working set `22.4MB`, energy
  `-13678.44097 kJ/mol`.
- The direct profile recorded `2,156,892` geometric-buffer pairs over the solve:
  `1,644,240` `gauss7` pairs and `512,652` subdivided adaptive pairs. This is
  strong evidence that Barnase has many close non-one-ring cleft pairs that the
  previous centroid path missed.
- Treecode TABI, `theta=0.5`, `leaf=32`, diagonal preconditioner, restart `64`,
  with the same near/buffer settings:
  completed in `54.62s`, `6` matvec calls, best residual `8.3515e-05`, final
  residual `2.2823e-03`, energy `-13678.73778 kJ/mol`, but used about
  `177.9M` leaf-panel interactions per matvec and peak working set
  `551.4MB`.
- Quick classification remains `treecode_noise_limited` because the harness
  currently judges final residual after restart, but the new best-residual
  fields show the direct and treecode rows both crossed the `1e-4` target
  before the final restart bookkeeping row. This should be handled explicitly
  in the next solver-reporting slice.

Stage 2W decision:

- The protein-cleft comment was correct. Born/Kirkwood spheres did not expose
  this because convex surfaces have fewer close non-one-ring triangle pairs.
  Barnase does expose it, and the geometric buffer materially changes the
  convergence picture.
- Near-field is no longer the obvious primary failure for the coarse Barnase
  truth row: direct reached best residual below `1e-4` in only `6` matvecs with
  the buffer enabled.
- A Born ion gate immediately caught that auto-enabling the buffer for every
  `near=duffy_adaptive` run can worsen at least one enforced sphere row
  (`q=1`, `r=1.5 A`, `salt=0.05 M`, `points=96`, latlong). The implementation
  therefore keeps `quad_buffer=auto -> off`; Barnase/protein probes must use
  explicit `-pbe-tabi-quad-buffer=adaptive` until Stage 2X analytic matrix
  checks pass.
- The next bottlenecks are now:
  1. GMRES/status semantics: preserve and accept the best residual when it
     crosses tolerance before a restart row worsens the final reported residual.
  2. Treecode memory/work: the buffer-aware treecode row is finite and close in
     energy to direct, but still has huge leaf-panel work and high peak memory.
  3. Decide whether `quad_buffer=adaptive` can become a hidden Barnase/protein
     perf candidate after Born/Kirkwood analytic gates confirm it does not
     worsen enforced rows.

Verification after Stage 2W:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`
  completed after an initial long full-rebuild timeout and an incremental
  rerun.
- `test/pbe_tabi_geometric_near_buffer_smoke.ps1`: PASS
  (`gap_ratio=0.212132`, `center_ratio=0.212132`).
- `test/pbe_tabi_near_kernel_blocks_smoke.ps1`: PASS
  (`k1_accept=0`, `k2_accept=3`, `k3_accept=0`, `k4_accept=0`).
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS
  (`panels=1422`, `rel3=9.38761e-07`, `rel4=1.35948e-07`).
- `test/pbe_tabi_profile_smoke.ps1`: PASS
  after keeping `quad_buffer=auto -> off`
  (`panels=5386`, `matvec_calls=30`, `E_kJ=-1113.71`).
- `test/pbe_tabi_barnase_truth_matrix.ps1 -Quick -CaseTimeoutSeconds 60
  -ProfileMemory -NearIntegration duffy_adaptive -QuadBuffer adaptive`:
  completed and wrote
  `test/tmp/pbe_tabi_barnase_truth_matrix/pbe_tabi_barnase_truth_matrix.csv`.
- Focused Born row after the `auto -> off` correction:
  `test/pbe_born_ion_matrix.ps1 -SkipGrid -RadiiA 1.5 -Charges 1.0
  -SaltM 0.05 -SurfacePointCounts 96`: PASS. A full default Born matrix was
  started but hit the local `600s` tool timeout before export, so Stage 2X
  should still run the full analytic matrix explicitly.
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.

Recommended Stage 2X/8Q:

1. Fix GMRES success/status semantics for TABI: if `best_residual <= tolerance`
   before restart bookkeeping, mark the solve converged, retain the best
   solution vector, and report both best and final residual without failing the
   frame.
2. Run enforced Born/Kirkwood matrices with explicit
   `near=duffy_adaptive` and `quad_buffer=adaptive`. If adaptive worsens sphere
   rows, add geometry-aware gating so convex/sphere cases keep centroid/ASI
   while protein clefts still get buffer quadrature.
3. Re-run Barnase quick/full truth with direct diagonal/block and restart
   `64/100`. If direct remains best-residual converged, switch the main
   performance work back to treecode leaf work/memory.
4. For treecode, reduce buffer-aware leaf work before more AVX2 math: build
   source-leaf/target-leaf pair batching that still routes close cleft pairs
   through the geometric buffer instead of expanding everything into panel
   interactions.

## Stage 2X/8Q Update - True Residual And On-The-Fly Treecode

Implemented in this slice:

- Added hidden `-pbe-tabi-interaction-mode=auto|cached|onthefly`.
  `auto` still resolves to cached behavior; `onthefly` is explicit debug/perf
  mode only.
- Added true-residual GMRES semantics for TABI:
  - Krylov iterations may still use the preconditioned system internally.
  - Solver status and `result.residual` now use the true residual of the
    original BIE system, `||b - A*x|| / ||b||`.
  - TABI stores the best solution vector by true residual and can return it
    with `convergence_reason=converged_true_residual_best` if a later restart
    row is worse.
  - Added separate preconditioned residual diagnostics:
    `tabi_best_preconditioned_residual`,
    `tabi_final_preconditioned_residual`,
    `tabi_true_residual_abs`, and `tabi_true_residual_rel`.
- Extended long-run progress output with the same true/preconditioned residual
  columns so timeout rows can distinguish physical convergence from
  preconditioned residual progress.
- Added on-the-fly JCP treecode traversal:
  - Reuses the tree nodes, Morton order, panel cache, and moments.
  - Does not build/store the full cached per-target interaction lists.
  - Routes self/one-ring/geometric-buffer pairs through the existing near path,
    so protein cleft quadrature is not bypassed.
  - Emits aggregate counters with `tabi_interaction_mode=onthefly`.
- Kept `quad_buffer=auto -> off`. The new diagnostic
  `tabi_quad_buffer_auto_candidate` is `off` for auto/off and `explicit` when
  the user pins `gauss7|adaptive`; no auto-enable policy was activated.

Measured Stage 2X smoke results:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS. Confirms timing,
  profile, and progress CSVs contain true/preconditioned residual columns and
  that `status=ok` is gated by true residual.
- `test/pbe_tabi_onthefly_treecode_smoke.ps1`: PASS. Cached and on-the-fly
  JCP treecode agree on the small sphere fixture, and on-the-fly is
  deterministic for `nt=1` vs `nt=4`.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS
  (`panels=1422`, `rel3=9.38761e-07`, `rel4=1.35948e-07`, `far=199`,
  `jcp_duffy=19416`).
- Focused Born row with `quad_buffer=auto -> off`: PASS for
  `r=1.5 A`, `q=1`, `salt=0.05 M`, `points=96`.
- `test/pbe_tabi_profile_smoke.ps1`: PASS after narrowing the smoke to the
  small sphere fixture so it validates profile schema without depending on a
  protein row that now correctly fails true-residual convergence
  (`panels=120`, `matvec_calls=5`, `E_kJ=-237.854`).
- `test/pbe_memory_profile_smoke.ps1`: PASS after narrowing to the sphere
  memory schema check. Protein memory/perf rows now belong in the truth/perf
  harness because true-residual semantics expose real TABI failures instead of
  accepting preconditioned-residual-only success.
- Regression gates unchanged:
  `test/pbe_regression.ps1`, `test/pbsa_smoke.ps1`,
  `test/sasa_surface_smoke.ps1`, and `test/pbe_clean_room_audit.ps1`: PASS.

Barnase note after Stage 2X:

- A quick 60s Barnase truth run with explicit
  `near=duffy_adaptive`, `quad_buffer=adaptive` now writes the new true and
  preconditioned residual columns, but the TABI rows timed out before a
  restart/true-residual recompute produced a useful best row. This is expected
  with true-residual logging because in-window GMRES progress cannot know the
  true residual until the current window is applied.
- Next Barnase diagnostic should use a longer direct truth window or a smaller
  surface before classifying. The important change is that we no longer mark
  a protein solve `ok` just because the preconditioned residual looked good.

Next recommended Stage 2Y:

1. Add an optional low-cost true-residual check inside long GMRES windows for
   Barnase/protein diagnostics, e.g. every `N` iterations or every progress
   interval, without making it the default for tiny cases.
2. Re-run Barnase direct truth with `surface_points=8`,
   `near=duffy_adaptive`, `quad_buffer=adaptive`, restart `64/100`, and enough
   timeout to reach at least one full restart after the first update.
3. Compare cached vs on-the-fly treecode on a medium protein fixture for RAM
   delta. If memory drops materially and energy stays within tolerance, use
   on-the-fly for Barnase profiling before further leaf-pair batching.
4. If direct true residual still converges but treecode stalls, tune
   `theta=0.5/0.6/0.7`, `leaf=64/128/256`, and Taylor order with the new true
   residual columns instead of relying on preconditioned residual status.

## Stage 2Y/8R Update - In-Window True Residual Probe

Implemented in this slice:

- Added hidden diagnostic probe controls:
  - `-pbe-tabi-true-residual-probe=restart|interval|off`
  - `-pbe-tabi-true-residual-probe-interval=<1..128>`
- Defaults remain unchanged: `restart` mode with interval `8`. The extra
  in-window matvecs run only when the user explicitly pins `interval`.
- Extended TABI GMRES so `interval` mode builds a temporary Krylov-window
  solution, applies the original BIE matvec, computes the true residual
  `||b - A*x|| / ||b||`, updates the best true solution, and can return
  `convergence_reason=converged_true_residual_probe` if the probe crosses
  tolerance.
- Extended CSV/progress output:
  `tabi_true_residual_probe_mode`,
  `tabi_true_residual_probe_interval`,
  `tabi_true_residual_probe_calls`, and
  `tabi_true_residual_probe_seconds` in `pbe_timing.csv`; matching non-`tabi_`
  columns in `pbe_tabi_profile.csv` and `pbe_tabi_longrun_progress.csv`.
- Fixed long-run progress semantics so probe counters/seconds are synchronized
  while GMRES is still running. Timeout rows now preserve probe work instead of
  showing `0` until solve completion.
- Updated `test/pbe_tabi_barnase_truth_matrix.ps1`:
  - quick rows now use direct diagonal/block restart `64`, plus cached and
    on-the-fly treecode rows;
  - Barnase/protein validation pins `near=duffy_adaptive`,
    `quad_buffer=adaptive`, `true_residual_probe=interval`, interval `4`;
  - added classifications including `true_probe_converged`,
    `probe_too_expensive`, `phantom_preconditioned`, and
    `true_residual_failed`.
- Updated `test/pbe_tabi_perf_matrix.ps1` with `-InteractionModes` so medium
  and Barnase probes can sweep `cached|onthefly`; default remains `cached`.

Measured Stage 2Y results:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS
  after an incremental rebuild of `titan_pbe/tabi_solver.cpp`.
- `test/pbe_tabi_true_residual_probe_smoke.ps1`: PASS. Confirms interval mode
  emits `gmres_true_probe` rows and nonzero probe call counts.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS. Confirms the default
  probe policy remains `restart` and interval `8`.
- `test/pbe_tabi_profile_smoke.ps1`: PASS
  (`panels=120`, `matvec_calls=5`, `E_kJ=-237.854`).
- `test/pbe_tabi_onthefly_treecode_smoke.ps1`: PASS.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS
  (`panels=1422`, `rel3=9.38761e-07`, `rel4=1.35948e-07`, `far=199`,
  `jcp_duffy=19416`).
- Focused Born row with `quad_buffer=auto -> off`: PASS for
  `r=1.5 A`, `q=1`, `salt=0.05 M`, `points=96`.
- `test/pbe_tabi_perf_matrix.ps1` tiny interaction-mode sweep:
  sphere treecode rows completed for both requested `cached` and `onthefly`.
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- Regression gates:
  `test/pbe_regression.ps1`, `test/pbsa_smoke.ps1`,
  `test/sasa_surface_smoke.ps1`, and `test/pbe_clean_room_audit.ps1`: PASS.

Barnase quick probe after Stage 2Y:

- Command:
  `test/pbe_tabi_barnase_truth_matrix.ps1 -Quick -CaseTimeoutSeconds 60
  -NearIntegration duffy_adaptive -QuadBuffer adaptive -ProfileMemory`
- Output:
  `test/tmp/pbe_tabi_barnase_truth_matrix/pbe_tabi_barnase_truth_matrix.csv`
- Grid reference completed:
  `1.12 s`, `E=-14879.32072 kJ/mol`, peak working set about `50.6 MB`.
- TABI direct diagonal, `surface_points=8`, `15129 panels`,
  `30258 unknowns`, timeout at `60.06 s`:
  `best_true_residual=1.423942550642e-01`,
  `best_preconditioned_residual=2.760553305897e-03`,
  `true_residual_probe_calls=1`,
  `true_residual_probe_seconds=6.877925`.
- TABI direct block, same surface, timeout at `60.41 s`:
  `best_true_residual=1.461759937450e-01`,
  `best_preconditioned_residual=2.918109817025e-03`,
  `true_residual_probe_calls=1`.
  The observed RSS row reached about `1.60 GB`, so current block setup is not a
  Barnase candidate.
- Treecode `theta=0.5`, `leaf=32`, diagonal, cached, timeout at `60.12 s`:
  `best_true_residual=1.423939605251e-01`,
  `leaf_panel_interactions=1.422935016e9`,
  RSS about `281.6 MB`,
  `true_residual_probe_calls=1`.
- Treecode `theta=0.5`, `leaf=32`, diagonal, on-the-fly, timeout at `60.06 s`:
  `best_true_residual=1.423939605251e-01`,
  `leaf_panel_interactions=1.067201262e9`,
  RSS about `33.0 MB`,
  `true_residual_probe_calls=1`.
- Decision classification for the 60s quick probe:
  `true_residual_failed`. This is an early-window diagnostic, not a final
  900s convergence result. It does show that true residual remains far above
  tolerance while the preconditioned residual is already around `2.8e-3`.

Recommended Stage 2Z:

1. Reduce the cost of true-residual probe for Barnase by probing against a
   sampled or blockwise residual first, then confirming full true residual only
   when the sampled residual drops. Full probe currently costs about `7 s` per
   call at `15129 panels`.
2. Disable or heavily gate block preconditioner for Barnase until the setup/apply
   memory blow-up is fixed; diagonal is the safer diagnostic baseline.
3. Keep on-the-fly tree traversal as the preferred Barnase memory strategy: the
   quick probe cut RSS from about `281.6 MB` cached to about `33.0 MB` while
   preserving the same early true residual.
4. The next performance bottleneck is still leaf work, not interaction-list RAM:
   on-the-fly reduces memory but still performs about `1.07B` leaf-panel
   interactions in 60s. Stage 2Z should implement target/source leaf-pair
   pruning or a stronger admissibility test for non-singular leaf pairs.
5. Re-run a longer diagonal-direct truth row after the lower-cost probe exists.
   If direct true residual trends downward below `1e-3`, return to treecode/MAC
   tuning; if it stagnates, continue near-field/preconditioner investigation.

## Stage 2Z/8S Update - Preconditioner Truth And Tree MAC Audit

Implemented in this slice:

- Added hidden TABI controls:
  - `-pbe-tabi-precond-block-near=singular|buffer-aware`
  - `-pbe-tabi-node-bounds=moment|tight-aabb`
- Kept defaults unchanged:
  - `precond=auto` still does not auto-enable block;
  - `precond-block-near` defaults to `singular`;
  - `quad_buffer=auto` remains `off`;
  - `interaction-mode=auto` remains the safe cached path.
- Refactored near-pair evaluation so block preconditioner can use the same
  self/one-ring/geometric-buffer kernel policy as direct/treecode near-field.
  `buffer-aware` block mode now includes geometric-buffer pairs through the
  `gauss7|subdivided` quadrature path instead of centroid.
- Added a block-preconditioner memory guard:
  - for large panel counts, block setup is skipped before allocating cached
    interaction lists or dense block data;
  - direct+block uses on-the-fly tree workspace when the tree is needed only
    for leaf grouping, preventing cached interaction-list RSS blow-up.
- Added tree MAC diagnostics:
  `tabi_node_bounds_mode`, `tabi_mac_accept_rate`,
  `tabi_node_radius_max`, `tabi_node_radius_mean`,
  `tabi_leaf_panel_work_per_matvec`, and
  `tabi_far_node_work_per_matvec` in timing/profile CSVs.
- Added tight-AABB tree node bounds mode for report-only MAC experiments. It
  computes node center/radius from a tight AABB over triangle vertices where
  available and recomputes moments around that center.
- Added/extended harnesses:
  - `test/pbe_tabi_block_precond_buffer_smoke.ps1`
  - `test/pbe_tabi_barnase_precond_truth.ps1`
  - `test/pbe_tabi_barnase_tree_mac_probe.ps1`
  - `test/pbe_tabi_barnase_truth_matrix.ps1` now supports `none`,
    `block-buffer-aware`, preconditioner block-near labels, tree node-bounds
    sweeps, and tree-probe-only mode.

Measured Stage 2Z smoke/probe results:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- `test/pbe_tabi_block_precond_buffer_smoke.ps1`: PASS. Confirms a
  non-topological close triangle pair is classified as geometric-buffer and
  the buffer-aware block option reaches `pbe_timing.csv`.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_onthefly_treecode_smoke.ps1`: PASS.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS
  (`panels=1422`, `rel3=9.38761e-07`, `rel4=1.35948e-07`, `far=199`,
  `jcp_duffy=19416`).
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.
- Regression gates: `test/pbe_regression.ps1`, `test/pbsa_smoke.ps1`, and
  `test/sasa_surface_smoke.ps1`: PASS.
- Barnase preconditioner truth probe:
  `test/pbe_tabi_barnase_precond_truth.ps1 -CaseTimeoutSeconds 20 -ProfileMemory`
  completed the grid row and wrote timeout rows for
  `none`, `diagonal`, `block`, and `block-buffer-aware`.
  - Grid reference: about `1.63 s`.
  - TABI direct rows used `15129 panels`, `30258 unknowns`.
  - With the new prebuild/on-the-fly guard, block rows no longer inflate RSS
    to ~`1.6 GB`; the short probe rows stayed around `32-41 MB`.
  - The 20s window is too short to reach the interval true-residual probe for
    Barnase; decision remains `true_residual_failed` for this short diagnostic.
- Barnase tree MAC probe syntax/progress check:
  `test/pbe_tabi_barnase_tree_mac_probe.ps1 -CaseTimeoutSeconds 1`: PASS.
  This validates the sweep harness and CSV plumbing for
  `leaf=64/128/256`, `theta=0.5/0.7/1.0/1.2`,
  `node_bounds=moment|tight-aabb`, and `interaction-mode=onthefly`.

Current interpretation:

- The comment was right about the block-preconditioner blind spot: singular-only
  block matrices missed geometric-buffer cleft pairs. The code now has an
  explicit buffer-aware path, but Barnase should still use diagonal/none as the
  truth baseline until a longer run proves block improves true residual.
- The previous Barnase RSS blow-up was caused largely by direct+block building
  a cached tree interaction workspace even when only leaf grouping was needed.
  That path now builds on-the-fly tree metadata instead.
- The 20s Barnase preconditioner probe is not a convergence result. It is only
  a safety/progress check showing memory is now controlled and the new matrix
  cases run.

Next recommended Stage 3A:

1. Run a longer Barnase diagonal-direct truth case with
   `surface_points=8`, `near=duffy_adaptive`, `quad_buffer=adaptive`,
   `restart=64`, and `true_residual_probe=interval` long enough to get at
   least one probe row. Use `none` as a secondary baseline.
2. If diagonal/none true residual drops below `1e-3`, run the tree MAC sweep
   for `onthefly`, `leaf=64/128/256`, `theta=0.5/0.7/1.0/1.2`, and
   `node_bounds=moment|tight-aabb`.
3. If direct true residual still stalls above `1e-3`, return to BIE
   near-field/CPV formulation and protein cleft quadrature before tuning
   treecode further.
4. Only revisit block or buffer-aware block after it improves true residual or
   iteration count versus diagonal without exceeding the memory guard.

## Stage 3A/8T - Sampled True-Residual Probe

Implemented:

- Added hidden TABI GMRES probe controls:
  - `-pbe-tabi-true-residual-probe-scope=full|sampled|adaptive`
  - `-pbe-tabi-true-residual-sample-count=<128..4096>`
- Default remains unchanged: probe mode `restart`, probe scope `full`, sample
  count `1024`.
- Direct JCP GMRES can now compute deterministic sampled true-residual probes
  over stratified unknown rows without running a full matvec. Sampled probes
  use the same direct near operator policy as the solve.
- Sampled probes are diagnostic only. They update sampled probe counters and
  progress rows but never set `status=ok`; convergence still requires full
  true residual at restart/end or full adaptive confirmation.
- Added CSV/progress fields:
  `tabi_true_residual_probe_scope`, `tabi_true_residual_sample_count`,
  `tabi_sampled_true_residual_rel`, `tabi_sampled_true_residual_rows`,
  `tabi_sampled_probe_calls`, `tabi_sampled_probe_seconds`,
  `tabi_full_true_residual_confirm_calls`, and `tabi_full_probe_seconds`.
- Added `gmres_sampled_true_probe` rows to `pbe_tabi_longrun_progress.csv`.
- Added/extended harnesses:
  - `test/pbe_tabi_sampled_true_residual_smoke.ps1`
  - `test/pbe_tabi_barnase_direct_truth_long.ps1`
  - `test/pbe_tabi_barnase_truth_matrix.ps1` now accepts sampled/adaptive
    probe scope and sample count.

Measured Stage 3A results:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
  - Sampled scope produced `6` sampled probes over `128` rows.
  - In-window full confirm calls stayed `0`.
  - Solver status was still gated by full true residual semantics
    (`converged_true_residual_after_preconditioned_check`), not by sampled
    residual.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS. Confirms default probe
  scope remains `full`.
- `test/pbe_tabi_true_residual_probe_smoke.ps1`: PASS. Confirms legacy
  interval full-probe behavior still works.
- `test/pbe_tabi_profile_smoke.ps1`: PASS.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS
  (`panels=1422`, `rel3=9.38761e-07`, `rel4=1.35948e-07`, `far=199`,
  `jcp_duffy=19416`).
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.
- Regression gates: `test/pbe_regression.ps1`, `test/pbsa_smoke.ps1`, and
  `test/sasa_surface_smoke.ps1`: PASS.
- Focused Born row:
  `test/pbe_born_ion_matrix.ps1 -SkipGrid -RadiiA 2.0 -Charges 1.0 -SaltM 0.15 -SurfacePointCounts 96 -NearIntegrationModes analytic_singular`:
  PASS.
- Full Born matrix with `-SkipGrid` exceeded a `300 s` local timeout in this
  run; no tolerance was changed.
- Barnase sampled direct-truth smoke:
  `test/pbe_tabi_barnase_direct_truth_long.ps1 -CaseTimeoutSeconds 20 -ProfileMemory`
  completed grid plus direct timeout rows.
  - Grid reference: about `0.99 s`.
  - TABI direct rows used `15129 panels`, `30258 unknowns`.
  - The 20s window reached only `2` GMRES iterations, so interval `4` did not
    yet emit sampled rows. This validates harness plumbing, not convergence.

Current interpretation:

- Stage 3A landed the diagnostic mechanism needed for Barnase, but the short
  20s Barnase run is still too short to classify direct convergence.
- The next decisive run should use the new wrapper with a longer timeout so at
  least two `gmres_sampled_true_probe` rows are written before stopping.
- Tree MAC tuning should remain gated behind a direct diagonal/none truth row
  that shows full true residual trending below `1e-3`.

Next recommended Stage 3B:

1. Run `test/pbe_tabi_barnase_direct_truth_long.ps1 -CaseTimeoutSeconds 300 -ProfileMemory`
   or a single-case diagonal/none variant if wall time must be capped.
2. If sampled residual trends down and triggers full confirm below `1e-3`,
   run `test/pbe_tabi_barnase_tree_mac_probe.ps1 -CaseTimeoutSeconds 300 -ProfileMemory`.
3. If sampled/full residual stalls above `1e-3`, return to BIE
   near-field/CPV conditioning for protein clefts before further treecode
   optimization.

## Stage 3B/8U - Top-Down Octree + Diagonal Protein Path

Implemented:

- Added hidden tree-builder control:
  - `-pbe-tabi-tree-builder=grouped-morton|topdown-aabb`
  - Default remains `grouped-morton`.
- Added `topdown-aabb` tree construction:
  - recursively partitions panel order by real AABB octants;
  - preserves contiguous panel spans per leaf;
  - stores child indices through `TreeWorkspace::child_indices`, so traversal
    and moment build no longer assume child nodes are contiguous in `nodes`.
- Updated scalar treecode, JCP treecode, cached interaction-list build,
  on-the-fly traversal, and packed moment translation to use the shared
  child-index helper.
- Added diagnostics to timing/profile data:
  `tabi_tree_builder`, `tabi_tree_depth_max`,
  `tabi_mac_reject_radius_limited`, plus the existing MAC accept rate and
  node-radius metrics.
- Updated Barnase MAC probe plumbing:
  `test/pbe_tabi_barnase_tree_mac_probe.ps1` now sweeps
  `tree_builder=grouped-morton|topdown-aabb` with on-the-fly traversal,
  `leaf=64/128/256`, `theta=0.5/0.7/1.0/1.2`, and diagonal preconditioning.
- Added a conservative block-preconditioner pivot-ratio guard. Ill-conditioned
  micro-block LU factors now fall back to diagonal with
  `preconditioner_effective=diagonal_block_ill_conditioned` instead of being
  applied to GMRES.
- Added `test/pbe_tabi_tree_builder_smoke.ps1`.

Measured Stage 3B results:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- `test/pbe_tabi_tree_builder_smoke.ps1`: PASS.
  - Synthetic spatial-jump fixture: `512 panels`.
  - `grouped_nodes=73`, `topdown_nodes=73`.
  - `grouped_radius_mean=0.665403`, `topdown_radius_mean=0.602967`.
  - `topdown_depth=2`.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS.
  - Default builder: `rel3=9.38761e-07`, `rel4=1.35948e-07`.
  - `topdown-aabb`: `topdown_rel3=4.61951e-05`, within the existing gate.
- `test/pbe_tabi_block_precond_buffer_smoke.ps1`: PASS after the pivot-ratio
  guard.
- CSV/profile sanity:
  - `test/pbe_tabi_profile_smoke.ps1`: PASS.
  - `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
  - `test/pbe_memory_profile_smoke.ps1`: PASS.

Current interpretation:

- The current `grouped-morton` default remains available and unchanged.
- `topdown-aabb` is now a safe hidden candidate: it preserves treecode
  accuracy on the existing smoke fixture and lowers mean node radius on the
  synthetic spatial-jump case.
- Barnase/protein optimized probes should pin `precond=diagonal` first. Block
  preconditioning remains report-only until it proves lower full true residual
  than diagonal without conditioning or memory regressions.

Next recommended Stage 3C:

1. Run the Barnase direct diagonal/none truth gate long enough to capture
   sampled true-residual rows:
   `test/pbe_tabi_barnase_direct_truth_long.ps1 -CaseTimeoutSeconds 300 -ProfileMemory`.
2. If direct residual trends below `1e-3`, run:
   `test/pbe_tabi_barnase_tree_mac_probe.ps1 -CaseTimeoutSeconds 300 -ProfileMemory`
   and compare `grouped-morton` vs `topdown-aabb` using
   `tabi_mac_accept_rate`, `tabi_mac_reject_radius_limited`,
   `tabi_leaf_panel_work_per_matvec`, RSS, and residual trend.
3. If direct still stalls above `1e-3`, continue with BIE near-field/CPV
   conditioning for protein clefts before spending more time on treecode.

## Stage 3C/8V - GMRES True-Residual Semantics + Barnase Focused Tree Probe

Implemented:

- Fixed TABI GMRES stop semantics:
  - preconditioned residual below tolerance now triggers a diagnostic true-residual
    check but no longer terminates the solve unless the original-system true
    residual also passes;
  - in-window GMRES continues after a phantom preconditioned event, allowing the
    sampled/full true residual trend to keep improving instead of stopping early.
- Made Barnase truth wrappers narrower and safer for iterative diagnostics:
  - `test/pbe_tabi_barnase_truth_matrix.ps1` can limit direct preconditioners and
    restart values;
  - `test/pbe_tabi_barnase_direct_truth_long.ps1` forwards those direct controls;
  - `test/pbe_tabi_barnase_tree_mac_probe.ps1` now accepts explicit leaf/theta,
    bounds, builder, interaction-mode, and sampled true-residual probe controls.

Measured Stage 3C results:

- Native rebuild after the GMRES fix:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Focused smoke tests:
  - `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS.
  - `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
  - Focused Barnase tree MAC wrapper plumbing with `CaseTimeoutSeconds=1`: PASS.
- Barnase direct truth, diagonal preconditioner, restart `64`, `surface_points=8`,
  `near=duffy_adaptive`, `quad_buffer=adaptive`, `nt=8`, timeout `300s`:
  - before the GMRES semantics fix, the same case stopped at `24` iterations with
    preconditioned residual `9.94e-5` but true residual `2.496e-3`;
  - after the fix, it kept running to the timeout and reached best full true
    residual `1.850e-4` at iteration `32`;
  - latest sampled true residual was `1.449e-4`;
  - panels `15129`, unknowns `30258`, RSS around `31 MB`;
  - harness classification: `diagonal_truth_candidate`.
- Focused Barnase treecode probe, `interaction-mode=onthefly`,
  `node_bounds=tight-aabb`, `leaf=128`, timeout `120s`:
  - `grouped-morton theta=0.7`: `4.599B` leaf interactions, sampled residual
    `2.32e-3`;
  - `topdown-aabb theta=0.7`: `0.734B` leaf interactions, sampled residual
    `5.84e-2`;
  - `grouped-morton theta=1.0`: `2.951B` leaf interactions, sampled residual
    `2.33e-3`;
  - `topdown-aabb theta=1.0`: `0.431B` leaf interactions, sampled residual
    `4.47e-2`;
  - `topdown-aabb theta=0.5`: `1.327B` leaf interactions, sampled residual
    `5.83e-2`.

Current interpretation:

- Direct Barnase is no longer a core BIE failure at this coarse density: with
  diagonal preconditioning and the geometric near buffer, true residual trends
  into the `1e-4` range.
- The previous `phantom_preconditioned` diagnosis was partly caused by solver
  stop semantics. Preconditioned residual remains useful as a progress signal,
  but full true residual must remain the acceptance source.
- `topdown-aabb` fixes a real work problem: it reduces Barnase leaf-panel work by
  roughly `4x-7x` in the focused probes. However, the current topdown treecode
  rows show much worse sampled residual than grouped-morton, even at `theta=0.5`.
  That points to a treecode consistency problem for the topdown path or its
  moment/near routing on protein-scale geometry, not just MAC looseness.

Next recommended Stage 3D:

1. Build a medium-size direct-vs-treecode diagnostic that compares
   `grouped-morton` and `topdown-aabb` with the same `duffy_adaptive` +
   `quad_buffer=adaptive` near policy, and reports sampled/full true residual,
   matvec relative error, and energy drift.
2. Audit topdown moment translation and near routing:
   - verify every child span/panel is visited exactly once;
   - verify singular, one-ring, and geometric-buffer pairs are not being routed
     into far-node Taylor for topdown;
   - compare far-node accepted pairs between grouped and topdown on the same
     target sample.
3. Keep Barnase optimized probes on `precond=diagonal`.
4. Do not promote `topdown-aabb` yet despite the work reduction; first fix the
   sampled-residual regression relative to grouped-morton.

## Stage 3D/8W - Topdown Far-Field Leakage Audit + Physical-Gap MAC

Implemented:

- Added hidden topdown MAC safety controls:
  - `-pbe-tabi-mac-safety=off|audit|physical-gap`, default `off`;
  - `-pbe-tabi-mac-gap-factor=<0.5..4.0>`, default `1.5`;
  - `-pbe-tabi-moment-audit=off|sampled|full-small`, default `off`.
- Added geometric far-field leakage accounting in tree traversal:
  - `audit` records far-node candidates whose AABB gap is too small relative to
    the target/source triangle scale;
  - `physical-gap` rejects those candidates and continues traversal, so
    geometric-buffer/near paths can handle cleft-like pairs.
- Added tree/node diagnostics to timing/profile outputs:
  - `tabi_mac_safety`, `tabi_mac_gap_factor`;
  - `tabi_far_geometric_unsafe_count`,
    `tabi_far_geometric_unsafe_panels`,
    `tabi_far_geometric_min_gap`,
    `tabi_far_geometric_gap_ratio_min/mean`,
    `tabi_far_geometric_audit_samples`,
    `tabi_mac_physical_gap_rejects`;
  - `tabi_moment_audit`, `tabi_moment_audit_nodes`,
    `tabi_moment_translation_rel_error_max/mean`.
- Added detail CSVs:
  - `pbe_tabi_far_geometric_audit.csv`;
  - `pbe_tabi_moment_audit.csv`.
- Updated Barnase tree MAC harnesses to forward `MacSafety`, `MacGapFactor`, and
  `MomentAudit`.
- Added `test/pbe_tabi_topdown_mac_audit.ps1`, a synthetic cleft smoke where
  topdown far-field can steal near work and `physical-gap` must block it.

Measured Stage 3D results:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- `test/pbe_tabi_topdown_mac_audit.ps1`: PASS.
  - Synthetic cleft fixture: `unsafe=392`, `rejects=392`;
  - far-node interactions dropped from `392` in audit mode to `0` with
    `physical-gap`.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS.
  - Default builder: `rel3=9.38761e-07`, `rel4=1.35948e-07`;
  - `topdown-aabb`: `topdown_rel3=4.61951e-05`.
- `test/pbe_tabi_profile_smoke.ps1`: PASS.
  - `panels=120`, `matvec_calls=9`, energy finite.
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- Short Barnase topdown MAC wrapper check:
  - command used `leaf=128`, `theta=0.7`, `node_bounds=tight-aabb`,
    `tree_builder=topdown-aabb`, `mac_safety=audit,physical-gap`,
    timeout `20s`;
  - wrapper completed and wrote `pbe_tabi_barnase_tree_mac_probe.csv` with
    timeout rows, confirming option plumbing.

Current interpretation:

- The synthetic cleft test confirms the suspected failure mode exists in the
  code path: a topdown far-node can be accepted even when the physical gap is
  inside the geometric-near safety zone, and `physical-gap` correctly prevents
  that acceptance.
- This does not yet prove Barnase topdown residual is fixed. The short Barnase
  run timed out before final timing rows, so the next useful run must be a
  longer Barnase/medium probe that compares sampled residual and leaf-work with
  `mac-safety=audit` versus `physical-gap`.
- Moment audit is now wired as a sampled diagnostic. If physical-gap improves
  unsafe far counts but residual remains bad, the next branch is
  `moment_translation_limited`.

Next recommended Stage 3E:

1. Run a medium fixture matrix with `topdown-aabb`, `tight-aabb`,
   `mac-safety=off|audit|physical-gap`, and compare direct-vs-treecode matvec
   error, energy drift, unsafe far counts, and sampled true residual.
2. Run Barnase `120-300s` probes with `physical-gap` after the medium fixture
   passes, keeping `precond=diagonal`, `interaction-mode=onthefly`,
   `near=duffy_adaptive`, and `quad_buffer=adaptive`.
3. If `physical-gap` fixes residual but restores too much leaf work, optimize
   leaf-pair batching under the new MAC safety.
4. If residual remains bad while unsafe far counts drop, turn on
   `-pbe-tabi-moment-audit=sampled` and fix moment translation/direct moment
   construction before more tree tuning.
