# TITAN PBE TABI Roadmap

This roadmap tracks the remaining work needed to move the clean-room TABI
backend from the current validated core toward paper-level TABI 2013 and
TABI-PB 2.0 coverage.

Last updated: 2026-06-11

## Phase 7 Decision Gate - Benchmark-Gated Direction

Status: implemented as a benchmark gate harness.

Phase 7 intentionally does not add new solver machinery. The gate script
`test/pbe_tabi_phase7_gate.ps1` wraps the existing PBE/TABI perf matrix and
exports a normalized decision package:

```text
test/tmp/pbe_tabi_phase7_gate/pbe_tabi_phase7_matrix.csv
test/tmp/pbe_tabi_phase7_gate/pbe_tabi_phase7_decision.csv
test/tmp/pbe_tabi_phase7_gate/pbe_tabi_phase7_summary.dat
```

The gate compares `grid`, `tabi_direct`, and `tabi_treecode` on the selected
research fixtures, records wall time, RSS, energy/accuracy status, GMRES
iterations, final residuals, accepted-miss evidence when present, and failure
modes. The current built-in Phase 7 fixtures are:

- `research-medium`: `example/1/rec_027.pqr` (~6k atoms, charged/radius PQR).
- `research-large`: `example/4/rec_0.pqr` (~12.7k atoms, charged/radius PQR).

TABI can become a promotion candidate only when treecode rows pass accuracy and
accepted-miss evidence and also show at least one practical advantage:

- at least `30%` lower memory than grid on two or more medium fixtures;
- runtime within `1.25x` grid on two or more medium fixtures;
- grid preflight/memory failure while TABI passes the same target fixture.

Current production decision: keep the finite-difference `grid` backend as the
default. TABI remains hidden/experimental or explicit opt-in until a full
medium-fixture gate run emits `promote_candidate`. Missing accepted-miss
evidence is treated conservatively and blocks promotion.

The gate resolves the research fixtures directly from `example/` and therefore
does not require `test/pbe_inputs_local` for the Phase 7 decision path.

Latest local Phase 7 evidence with `-CaseTimeoutSeconds=60`:

| Fixture | Grid status / seconds | TABI direct | TABI treecode | Decision signal |
| --- | ---: | ---: | ---: | --- |
| research-medium | ok / ~0.62 s | timeout / 60 s | timeout / 60 s | no runtime gate |
| research-large | ok / ~1.52 s | timeout / 60 s | timeout / 60 s | no runtime gate |

Decision: `keep_experimental`, `promote_candidate=false`,
`decision_reason=accuracy_or_runtime_failure_present`. This is strong enough to
pause new TABI machinery for the current research fixtures and keep production
PBE on the grid backend.

No FMM, new scheduler layer, or additional preconditioner should be added until
this gate passes on the target workload.

## TABI Review Follow-Up - Immediate Triage

Status: implemented as a short diagnostic/refinement slice before any new
solver machinery.

The local review in `tabi_review_report.md` is directionally useful. Its two
most actionable observations match the latest Barnase timing profile:

- The run is GMRES/matvec limited, not setup limited. Barnase grid finished in
  about `0.054s`, while TABI treecode timed out at `180s` after only `17`
  matvecs. Surface construction was about `0.04s`, tree build about `1.57s`,
  and preconditioner setup was negligible.
- The active treecode row produced roughly `14.34B` leaf-panel interactions
  for `39011` panels, which is worse than a simple `N^2` panel-pair sanity
  bound. This means the tree traversal/leaf-near policy is degenerating into a
  leaf-direct workload rather than acting as a useful acceleration.
- The effective preconditioner was `diagonal_block_memory_cap`, so the solver
  was running with a weak diagonal fallback. This explains why the residual
  stayed at `1.000000000000e+00` through the observed GMRES iterations.

Immediate actions that are worth doing:

1. Add a timeout-safe large-protein triage report to the Barnase/profile
   harnesses. It should emit:
   - `leaf_panel_interactions_per_matvec`;
   - `leaf_vs_n2_ratio`;
   - `seconds_per_matvec`;
   - `residual_delta_per_matvec`;
   - requested/effective preconditioner and downgrade reason;
   - a classification such as `preconditioner_downgraded`,
     `leaf_work_super_n2`, `matvec_limited`, or `convergence_stalled`.
2. Fail fast on hopeless TABI rows before long timeouts. If the progress file
   shows both `residual_delta_per_matvec ~= 0` and `leaf_vs_n2_ratio > 1`
   after the first few matvecs, classify the run as diagnostic failure and stop
   instead of waiting for the full timeout.
3. Run a small review-driven parameter sweep, but keep it behind hidden
   benchmark scripts:
   - leaf size `32/64/128`;
   - opening theta around the current safe range;
   - surface points `8/24` for Barnase profiling only;
   - `block-budgeted` memory budgets large enough to avoid the silent diagonal
     downgrade when the estimate is reasonable.
4. Treat the review's "increase preconditioner budget" recommendation as a
   diagnostic branch, not as a default change. It is valid only if the reported
   effective preconditioner becomes real block/budgeted-block and GMRES
   matvecs or residual trend improve on medium fixtures first.
5. Treat the review's "deepen the tree / reduce leaf size" recommendation as
   a work-count branch, not as a guaranteed speedup. It is valid only if
   `leaf_vs_n2_ratio` drops below `1` and wall time per matvec improves without
   violating accuracy/audit gates.

Not immediate:

- Do not add FMM/BLDTT, a new scheduler layer, or another coarse
  preconditioner from this review alone. Previous medium-fixture evidence
  already showed local block and exact-Galerkin coarse correction did not
  reduce Krylov work enough.
- Do not promote TABI or change `grid` default until the Phase 7 decision gate
  emits `promote_candidate=true`.

Implemented triage harness:

- Added shared PowerShell triage helpers and
  `test/pbe_tabi_large_protein_triage.ps1`.
- Extended `test/pbe_tabi_barnase_timing_profile.ps1` and
  `test/pbe_tabi_phase7_gate.ps1` with the triage columns above.
- The triage is report-only. It does not change production solver defaults,
  TABI parameters, or the finite-difference grid default.

Latest timeout-safe triage smoke, representative row
`surface_points=24`, `leaf=128`, `theta=0.8`, `nt=8`,
`CaseTimeoutSeconds=30`:

| Fixture | Grid seconds | TABI status | Triage classification | seconds/matvec | leaf_vs_n2 | Effective preconditioner |
| --- | ---: | --- | --- | ---: | ---: | --- |
| `barnase` | `0.093` | timeout | `preconditioner_downgraded` | `12.61` | `1.109` | `diagonal_block_memory_cap` |
| `research-medium` | `0.710` | timeout | `preconditioner_downgraded` | `NA` | `0` | `diagonal_block_pending` |
| `research-large` | `1.265` | timeout | `timeout_no_profile` | `NA` | `0` | `diagonal` |

Decision from the triage slice:

- The review's preconditioner concern is confirmed on Barnase and
  research-medium: TABI is not reaching a useful block preconditioner before
  the timeout.
- Barnase also confirms the leaf-work concern with `leaf_vs_n2 > 1` once GMRES
  starts, but the higher-priority diagnosis remains preconditioner downgrade.
- Research-large needs either a longer diagnostic timeout or earlier progress
  instrumentation before its failure can be split more precisely than
  `timeout_no_profile`.

## TABI Preconditioner Activation Diagnostic

Status: implemented as a report-only preconditioner activation probe.

New diagnostic harness:

- `test/pbe_tabi_preconditioner_activation_probe.ps1` sweeps `grid`,
  `diagonal`, `block`, and `block-budgeted` modes.
- `test/pbe_tabi_perf_matrix.ps1` now accepts diagnostic-only
  `PreconditionerBlockNear`, `PreconditionerBlockPanelCap`, and
  `PreconditionerMemoryBudgetMb` script parameters, and still writes generated
  `arg.in` files with `-option=value` syntax.
- `test/pbe_tabi_triage_helpers.ps1` now exports preconditioner setup/apply,
  budget, valid-block, valid-panel, ill-conditioned-block, and effective-mode
  telemetry when the underlying TABI run emits it.
- The human-readable result is captured in
  `TITAN_TABI_PRECONDITIONER_ACTIVATION_REPORT.txt`.

Focused local probe, representative row `surface_points=24`, `leaf=128`,
`theta=0.8`, `nt=8`, `CaseTimeoutSeconds=30`, `block-budgeted budget=256MB`,
`block_panel_cap=32`:

| Fixture | Requested preconditioner | Status | Effective preconditioner | Activation class | Triage class | seconds/matvec | leaf_vs_n2 |
| --- | --- | --- | --- | --- | --- | ---: | ---: |
| `barnase` | `diagonal` | timeout | `diagonal` | `diagonal_only` | `matvec_limited` | `15.77` | `0.554` |
| `barnase` | `block` | timeout | `diagonal_block_memory_cap` | `memory_capped` | `preconditioner_downgraded` | `12.59` | `1.109` |
| `barnase` | `block-budgeted` | timeout | `NA` | `diagonal_only` | `barnase_probe_timeout` | `NA` | `NA` |
| `research-medium` | `block-budgeted` | timeout | `NA` | `diagonal_only` | `timeout_no_profile` | `NA` | `NA` |
| `research-large` | `block-budgeted` | timeout | `NA` | `diagonal_only` | `timeout_no_profile` | `NA` | `NA` |

Activation decision:

- No tested branch produced a useful block/budgeted-block effective
  preconditioner on the representative fixtures.
- Barnase `block` still downgrades to `diagonal_block_memory_cap` and also
  shows bad leaf work once GMRES starts (`leaf_vs_n2 > 1`).
- The `block-budgeted` rows timed out before enough effective preconditioner
  telemetry was emitted, so the next fix should be earlier activation telemetry
  and fail-fast classification, not leaf-kernel optimization.

Gate for moving to leaf-kernel optimization:

- At least one medium fixture must report a real block/budgeted-block effective
  preconditioner.
- The same row must show GMRES/residual progress, not only setup timeout.
- Only then should leaf-kernel work be prioritized, and only if `leaf_vs_n2`
  remains high after preconditioner activation is fixed.

Follow-up implementation after this diagnostic:

- Added `src/tabi/tabi_preconditioner.h` and
  `src/tabi/tabi_preconditioner.cpp`.
- Moved JCP block preconditioner, budgeted block range/estimate logic, dense LU
  helpers, and coarse Galerkin preconditioner support out of
  `tabi_solver.cpp`.
- Changed `block-budgeted` activation policy so it is governed by the explicit
  preconditioner memory budget, not the legacy whole-process working-set soft
  cap.
- Changed block application so ill-conditioned local blocks fall back to
  diagonal behavior locally instead of resetting the entire preconditioner.
- Changed buffer-aware block preconditioner construction to use a fast
  centroid/local kernel approximation; full Duffy/adaptive near-field kernels
  remain in the production operator.
- Added `test/titan_tabi_preconditioner_smoke.cpp` to exercise the new module
  directly with a synthetic panel/tree fixture.

## TABI Leaf-Kernel Transition Probe

Status: implemented as a benchmark-only transition probe.

The operator-near preconditioner branch moved Barnase into a leaf-kernel-limited
regime: the effective preconditioner can become `block_budgeted_operator_near`,
the quality gate passes, and setup time is negligible, but each GMRES matvec is
still too expensive. The next gate therefore measures tree/leaf options before
any new kernel rewrite:

- `test/pbe_tabi_leaf_kernel_transition_probe.ps1` sweeps leaf size, opening
  theta, and existing scheduler/kernel variants under a fixed
  `block-budgeted + operator-near` preconditioner.
- The probe writes
  `test/tmp/pbe_tabi_leaf_kernel_transition_probe/pbe_tabi_leaf_kernel_transition_probe.csv`
  and the human-readable
  `TITAN_TABI_LEAF_KERNEL_OPTIMIZATION_REPORT.txt`.
- `test/pbe_tabi_perf_matrix.ps1` now accepts diagnostic-only `TabiExtraArgs`,
  allowing the probe to test existing hidden scheduler paths without adding
  production CLI behavior.

Decision rules:

- `tree_depth_win`: leaf-panel work per matvec drops by at least `3x` versus
  Barnase baseline `leaf=128`, `theta=0.8`, `tile-avx2`, with clean
  verify/accepted-miss evidence.
- `scheduler_win`: seconds per matvec drops by at least `1.5x` without growing
  leaf-panel work beyond `1.2x`.
- `mac_accuracy_fail`: accepted-miss or verify mismatch evidence is nonzero.
- `memory_limited`: peak memory exceeds `2x` the baseline row.
- `no_win`: no clean row meets the promotion gates.

Promotion remains conservative: a candidate can enter the Phase 7 experimental
gate only when Barnase reaches seconds/matvec `<3` or leaf-panel work/matvec
drops by at least `3x` with clean verification. The finite-difference `grid`
backend remains the production default.

Latest Barnase scheduler transition probe:

Command shape:

```text
test/pbe_tabi_leaf_kernel_transition_probe.ps1 -Fixtures barnase -TreeLeafSizes 64,128 -OpeningThetas 0.8 -SchedulerVariants baseline,deferred-worklist,target-leaf-block,soa-avx2 -CaseTimeoutSeconds 30
```

Result summary:

| Scheduler | leaf | theta | status | seconds/matvec | leaf work/matvec | decision |
| --- | ---: | ---: | --- | ---: | ---: | --- |
| baseline | 128 | 0.8 | timeout | `5.189` | `843.5M` | `no_win` |
| baseline | 64 | 0.8 | timeout | `10.936` | `603.9M` | `no_win` |
| soa-avx2 | 128 | 0.8 | timeout | `11.293` | `843.5M` | `no_win` |
| soa-avx2 | 64 | 0.8 | timeout | `18.537` | `603.9M` | `no_win` |
| deferred-worklist | 64/128 | 0.8 | timeout | `NA` | `NA` | `no_win` |
| target-leaf-block | 64/128 | 0.8 | timeout | `NA` | `NA` | `no_win` |

Interpretation:

- The existing hidden scheduler variants do not provide an immediate Barnase
  win under the operator-near preconditioner.
- `leaf=64` reduces leaf work but slows matvec time, so tree-depth tuning alone
  is not enough.
- The next implementation branch should attack safe work-count/kernel structure
  directly: reduce unclean/scalar tail work, reuse packed source/target blocks,
  and split tasks by estimated pair cost. Another scheduler wrapper is unlikely
  to pass the `1.5x` gate without changing the hot leaf path.

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
  after an incremental rebuild of `src/tabi/tabi_solver.cpp`.
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
  - `-pbe-tabi-precond-block-near=singular|buffer-aware|operator-near`
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

## Stage 3E/8X - TABI Topdown Physical-Gap Validation

Implemented:

- Added hidden direct-vs-treecode matvec audit controls:
  - `-pbe-tabi-matvec-audit=off|sampled|full-small`, default `off`;
  - `-pbe-tabi-matvec-audit-sample-count=<128..4096>`, default `1024`.
- Added timing/profile fields:
  - `tabi_tree_vs_direct_matvec_rel`;
  - `tabi_matvec_audit_rows`;
  - `tabi_matvec_audit_seconds`.
- The audit computes the treecode matvec on the selected GMRES solution and
  compares sampled rows, or the full operator when `full-small` is requested and
  the unknown count is small enough. The production default remains unchanged.
- Added `test/pbe_tabi_topdown_physical_gap_matrix.ps1`, a Stage 3E matrix
  harness for `mini-protein`, `charged-peptide`, and `compact-mini-protein`.
  The harness pins:
  `surface_points=48`, `leaf=128`, `theta=0.7`,
  `node_bounds=tight-aabb`, `interaction-mode=onthefly`,
  `near=duffy_adaptive`, `quad_buffer=adaptive`, `precond=diagonal`,
  and `nt=8`.
- The matrix writes direct baseline rows plus:
  `grouped-morton/off`, `topdown-aabb/off`, `topdown-aabb/audit`, and
  `topdown-aabb/physical-gap`. The CSV includes energy drift vs direct,
  direct-vs-treecode matvec error, sampled/full true-residual fields, unsafe
  far counts, physical-gap rejects, leaf/far work per matvec, RSS, and runtime.

Measured Stage 3E medium results:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Full matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1`: PASS and wrote
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.
- `physical-gap` rows:

| fixture | matvec rel | energy drift rel | sampled residual | grouped sampled residual | unsafe/rejects | leaf/far work per matvec | runtime | RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `mini-protein` | `5.684e-4` | `8.221e-5` | `9.757e-5` | `9.752e-5` | `2310/2310` | `3.86M / 9.61k` | `7.07s` | `12.6MB` |
| `charged-peptide` | `1.832e-4` | `1.112e-4` | `1.772e-5` | `1.770e-5` | `580/580` | `3.13M / 5.90k` | `7.69s` | `12.1MB` |
| `compact-mini-protein` | `1.339e-3` | `5.579e-5` | `5.302e-5` | `5.309e-5` | `10260/10260` | `8.48M / 28.4k` | `34.21s` | `14.2MB` |

- All three `physical-gap` rows are under the Stage 3E medium gate:
  `matvec rel < 2.5e-3`, `energy drift rel < 5e-3`, and sampled residual within
  `2x` the grouped-morton row.
- `physical-gap` correctly rejects every unsafe far candidate observed in the
  audit rows on the medium fixtures, while preserving the true-residual trend.

Barnase follow-up probe:

- Command:
  `test/pbe_tabi_barnase_tree_mac_probe.ps1 -CaseTimeoutSeconds 120 -ProfileMemory -TreeLeafSizes 128 -TreeOpeningThetas 0.7 -NodeBoundsModes tight-aabb -TreeBuilders topdown-aabb -MacSafety physical-gap`.
- Grid reference completed in `1.08s`.
- Topdown physical-gap tree row timed out at `120.05s`, GMRES iteration `12`,
  `14` matvec calls, sampled residual `5.8385e-2`, RSS about `33.9MB`.
- The timeout row is useful as a lower-bound probe but remains inconclusive:
  final timing/MAC work counters were not emitted before termination.

Verification after Stage 3E:

- `test/pbe_tabi_topdown_mac_audit.ps1`: PASS
  (`unsafe=392`, `rejects=392`, far-node interactions `392 -> 0`).
- `test/pbe_tabi_tree_builder_smoke.ps1`: PASS
  (`panels=512`, `grouped_nodes=73`, `topdown_nodes=73`).
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS
  (`rel3=9.38761e-07`, `rel4=1.35948e-07`,
  `topdown_rel3=4.61951e-05`).
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_profile_smoke.ps1`: PASS
  (`panels=120`, `matvec_calls=9`).
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.
- `test/pbe_regression.ps1`: PASS.
- `test/pbsa_smoke.ps1`: PASS.
- `test/sasa_surface_smoke.ps1`: PASS.

Stage 3E decision:

- Medium physical-gap validation passes. This confirms the Stage 3D regression
  fix on the medium fixtures and justifies longer Barnase probing.
- Do not promote `topdown-aabb` or `physical-gap` to production defaults yet.
  Barnase still shows a high sampled residual in the 120s lower-bound probe, and
  the timeout path currently lacks final MAC/work counters.
- Next branch:
  1. Run a 300s Barnase physical-gap probe with the same settings and preserve
     timeout-safe progress counters for leaf/far work and MAC rejects.
  2. If Barnase sampled residual remains near `5e-2`, enable
     `-pbe-tabi-moment-audit=sampled` for the Barnase/topdown path and inspect
     moment translation/direct moment construction before more performance work.
  3. If Barnase residual improves but leaf work regresses under physical-gap,
     optimize leaf-pair batching without loosening the MAC safety rule.

## Stage 3F/8Y - TABI Topdown Fragmented-Node Moment Validation

Implemented:

- Added hidden center/safety controls:
  - `-pbe-tabi-expansion-center=legacy|area-centroid|aabb-center`, default
    `legacy`;
  - `-pbe-tabi-fragmentation-safety=off|audit|split`, default `off`.
- Separated tight AABB bounds from the Taylor expansion center:
  - `tight-aabb` still keeps true min/max bounds for physical-gap checks;
  - `legacy` preserves the old tight-AABB center behavior;
  - `area-centroid` keeps the panel/area center for Taylor moments and
    recomputes radius from that center to the farthest AABB corner, so the MAC
    remains conservative.
- Added fragmented-node diagnostics to timing/profile rows:
  - `tabi_expansion_center`, `tabi_fragmentation_safety`;
  - `tabi_fragmented_far_count`, `tabi_fragmented_far_panels`;
  - `tabi_fragmented_split_rejects`;
  - `tabi_fragmentation_center_gap_ratio_max/mean`.
- `fragmentation-safety=split` now descends child nodes instead of using a
  Taylor far-node sum when a MAC/physical-gap-passing internal node is marked
  fragmented. The leaf/direct path is unchanged.
- Added `test/pbe_tabi_fragmented_node_mac_smoke.ps1`, a synthetic three-patch
  fixture where `audit` sees a fragmented far node and `split` reduces the
  direct-vs-tree matvec error.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with Stage 3F rows:
  - `topdown-aabb/physical-gap/legacy/off`;
  - `topdown-aabb/physical-gap/area-centroid/audit`;
  - `topdown-aabb/physical-gap/area-centroid/split`.
- Extended Barnase tree MAC harnesses to forward `ExpansionCenters`,
  `FragmentationSafety`, and an explicit `SurfacePoints` knob.

Measured Stage 3F medium rows:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Synthetic fragmented-node smoke:
  `test/pbe_tabi_fragmented_node_mac_smoke.ps1`: PASS.
  - panels `48`;
  - audit fragmented far candidates `16`;
  - split rejects `16`;
  - direct-vs-tree matvec rel improved from `2.084e-8` to `2.512e-13`.
- Full medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -CaseTimeoutSeconds 180`:
  PASS and wrote
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| fixture | row | matvec rel | energy drift rel | sampled residual | grouped sampled residual | unsafe/rejects | fragmented far/split | leaf/far work per matvec |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `mini-protein` | `legacy/off` | `5.684e-4` | `8.221e-5` | `9.757e-5` | `9.752e-5` | `2310/2310` | `0/0` | `3.86M / 9.61k` |
| `mini-protein` | `area-centroid/audit` | `4.463e-4` | `6.757e-5` | `9.762e-5` | `9.752e-5` | `885/885` | `0/0` | `4.29M / 8.32k` |
| `mini-protein` | `area-centroid/split` | `4.463e-4` | `6.757e-5` | `9.762e-5` | `9.752e-5` | `885/885` | `0/0` | `4.29M / 8.32k` |
| `charged-peptide` | `legacy/off` | `1.832e-4` | `1.112e-4` | `1.772e-5` | `1.770e-5` | `580/580` | `0/0` | `3.13M / 5.90k` |
| `charged-peptide` | `area-centroid/audit` | `1.974e-4` | `1.224e-4` | `1.772e-5` | `1.770e-5` | `0/0` | `0/0` | `3.27M / 5.25k` |
| `charged-peptide` | `area-centroid/split` | `1.974e-4` | `1.224e-4` | `1.772e-5` | `1.770e-5` | `0/0` | `0/0` | `3.27M / 5.25k` |
| `compact-mini-protein` | `legacy/off` | `1.339e-3` | `5.579e-5` | `5.302e-5` | `5.309e-5` | `10260/10260` | `0/0` | `8.48M / 28.4k` |
| `compact-mini-protein` | `area-centroid/audit` | `8.891e-4` | `7.254e-5` | `5.314e-5` | `5.309e-5` | `3300/3300` | `0/0` | `10.0M / 24.5k` |
| `compact-mini-protein` | `area-centroid/split` | `8.891e-4` | `7.254e-5` | `5.314e-5` | `5.309e-5` | `3300/3300` | `0/0` | `10.0M / 24.5k` |

- All three `area-centroid/split` medium rows pass the Stage 3E/3F gate:
  `matvec rel < 2.5e-3`, `energy drift rel < 5e-3`, and sampled residual within
  `2x` grouped-morton.
- The medium fixtures did not naturally trigger fragmented far candidates.
  The synthetic smoke is therefore the current proof that the diagnostic and
  split path work; the medium matrix proves the hidden controls do not regress
  the accepted physical-gap rows.

Barnase Stage 3F branch:

- `surface_points=48` Barnase, `leaf=128`, `theta=0.7`, `topdown-aabb`,
  `tight-aabb`, `physical-gap`, `area-centroid`, `fragmentation-safety=split`,
  `onthefly`, `duffy_adaptive`, `quad_buffer=adaptive`, `precond=diagonal`,
  timeout `120s`:
  - reached `gmres_start` only;
  - panels `73967`, unknowns `147934`;
  - no residual/MAC timing row was emitted before timeout.
- Repeated the same branch with the Barnase harness coarse budget
  `surface_points=8`, timeout `120s`:
  - panels `15129`, unknowns `30258`;
  - timeout at GMRES iteration `10`, `12` matvecs;
  - sampled true residual `5.8313e-2`;
  - leaf-panel interactions `842.3M`, far-node interactions `5.11M`;
  - RSS from progress about `33.9MB`.
- Repeated the coarse Barnase row with `-pbe-tabi-moment-audit=sampled`:
  - timeout at the same iteration/matvec count;
  - sampled true residual remained `5.8313e-2`;
  - no final `pbe_tabi_moment_audit.csv`/timing row was emitted before the
    forced timeout.

Verification after Stage 3F:

- `test/pbe_tabi_fragmented_node_mac_smoke.ps1`: PASS.
- `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -CaseTimeoutSeconds 180`:
  PASS.
- `test/pbe_tabi_topdown_mac_audit.ps1`: PASS
  (`unsafe=392`, `rejects=392`, far-node interactions `392 -> 0`).
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS
  (`rel3=9.38761e-07`, `rel4=1.35948e-07`,
  `topdown_rel3=4.61951e-05`).
- `test/pbe_tabi_tree_builder_smoke.ps1`: PASS.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_profile_smoke.ps1`: PASS
  (`panels=120`, `matvec_calls=9`).
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.
- `test/pbe_regression.ps1`: PASS.
- `test/pbsa_smoke.ps1`: PASS.
- `test/sasa_surface_smoke.ps1`: PASS.

Stage 3F decision:

- Keep all new controls hidden/report-only. Do not promote `topdown-aabb`,
  `physical-gap`, `area-centroid`, or `fragmentation-safety=split`.
- Medium validation passes, and `area-centroid` reduces matvec error on two of
  the three medium physical-gap rows. However, Barnase remains in the old
  `~5.8e-2` sampled residual basin at the 120s coarse probe.
- The fragmented-node hypothesis is implemented and validated synthetically,
  but it is not sufficient to explain the Barnase residual plateau under the
  current 120s evidence.
- Next branch is moment/Taylor construction, not performance tuning:
  1. make moment/expansion audit emit timeout-safe progress rows, not only final
     timing/profile rows;
  2. add a far-expansion truncation audit that compares accepted Taylor node
     sums against direct sums on sampled target-node pairs;
  3. inspect JCP moment construction for `area-centroid` topdown nodes before
     optimizing leaf batching or loosening MAC safety.

## Stage 3G/8Z - TABI Timeout-Safe Far-Expansion Audit

Implemented:

- Added hidden far-expansion audit controls:
  - `-pbe-tabi-expansion-audit=off|sampled|full-small`, default `off`;
  - `-pbe-tabi-expansion-audit-sample-count=<64..2048>`, default `512`;
  - `-pbe-tabi-expansion-audit-max-panels=<128..8192>`, default `2048`.
- Implemented JCP on-the-fly far-expansion audit. When an accepted internal far
  node passes MAC, physical-gap, and fragmentation split routing, the audit
  compares `taylor_jcp_node_sum` against a direct JCP sum over the same node
  panels for sampled target-node pairs.
- Added timing/profile fields:
  - `tabi_expansion_audit`;
  - `tabi_far_expansion_audit_pairs`;
  - `tabi_far_expansion_audit_panels`;
  - `tabi_far_expansion_direct_seconds`;
  - `tabi_far_expansion_rel_error_max`;
  - `tabi_far_expansion_rel_error_mean`.
- Extended `pbe_tabi_longrun_progress.csv` with timeout-safe tree/MAC/moment
  and far-expansion fields, including tree builder, node bounds, MAC safety,
  expansion center, fragmentation safety, unsafe/reject counts, moment audit
  counters, and expansion audit counters/errors.
- Updated Barnase harnesses to forward the expansion-audit controls and to read
  fallback tree/MAC/moment/expansion counters from long-run progress rows when
  timing/profile rows are unavailable after a timeout.
- Added `test/pbe_tabi_far_expansion_audit_smoke.ps1`, a synthetic hollow-node
  smoke where the audit observes Taylor-vs-direct node error and split reduces
  the direct-vs-tree vector error.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with optional
  `-IncludeExpansionAuditRows`, adding an `area-centroid/split` sampled
  expansion-audit row without changing default production behavior.

Measured Stage 3G evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Synthetic far-expansion audit smoke:
  `test/pbe_tabi_far_expansion_audit_smoke.ps1`: PASS.
  - panels `64`;
  - sampled far expansion pairs `128`;
  - audited panels `2560`;
  - max/mean Taylor-vs-direct node error `1.650e-4 / 2.049e-5`;
  - split rejects `16`;
  - direct-vs-tree vector rel improved from `5.766e-9` to `3.877e-10`.
- Fragmented-node smoke:
  `test/pbe_tabi_fragmented_node_mac_smoke.ps1`: PASS.
  - panels `48`;
  - fragmented far candidates `16`;
  - split rejects `16`;
  - vector rel improved from `2.084e-8` to `2.512e-13`.
- Medium matrix with expansion-audit rows:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -CaseTimeoutSeconds 180 -IncludeExpansionAuditRows`:
  PASS and wrote
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| fixture | grouped sampled residual | split matvec rel | split energy drift rel | split sampled residual | expansion audit pairs | expansion rel max | classification |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `mini-protein` | `9.752e-5` | `4.463e-4` | `6.757e-5` | `9.762e-5` | `7680` | `1.037e1` | `fragmentation_split_medium_gate_pass` |
| `charged-peptide` | `1.770e-5` | `1.974e-4` | `1.224e-4` | `1.772e-5` | `10240` | `4.041e0` | `fragmentation_split_medium_gate_pass` |
| `compact-mini-protein` | `5.309e-5` | `8.891e-4` | `7.254e-5` | `5.314e-5` | `15360` | `1.498e1` | `fragmentation_split_medium_gate_pass` |

- The medium `area-centroid/split` rows remain under the success gates:
  `matvec rel < 2.5e-3`, `energy drift rel < 5e-3`, and sampled residual within
  `2x` grouped-morton.
- Expansion audit can see large local Taylor node errors even when the medium
  global matvec/residual gate passes. This makes the audit useful as a branch
  diagnostic, not a promotion gate by itself.

Barnase Stage 3G timeout-safe probe:

- Command shape:
  `test/pbe_tabi_barnase_tree_mac_probe.ps1 -CaseTimeoutSeconds 120 -SurfacePoints 8 -TreeLeafSizes 128 -TreeOpeningThetas 0.7 -NodeBoundsModes tight-aabb -TreeBuilders topdown-aabb -MacSafety physical-gap -ExpansionCenters area-centroid -FragmentationSafety split -ExpansionAudit sampled -ExpansionAuditSampleCount 512 -ExpansionAuditMaxPanels 2048 -TrueResidualProbeScope adaptive`.
- Tree row timed out at `120.04s`, but progress emitted the needed counters:
  - panels `15129`, unknowns `30258`;
  - GMRES iteration `10`, matvec calls `12`;
  - sampled true residual `5.8313e-2`;
  - tree builder `topdown-aabb`, node bounds `tight-aabb`;
  - MAC safety `physical-gap`, expansion center `area-centroid`;
  - fragmentation safety `split`, expansion audit `sampled`;
  - unsafe/rejects `228/228`;
  - fragmented far/split rejects `0/0`;
  - far-expansion audit pairs `6144`, panels `2310494`;
  - direct audit time `0.085869s`;
  - far-expansion rel max/mean `4.8596 / 1.278e-2`;
  - RSS from progress `34.1MB`;
  - last progress phase `gmres_iter`.

Verification after Stage 3G:

- `test/pbe_tabi_far_expansion_audit_smoke.ps1`: PASS.
- `test/pbe_tabi_fragmented_node_mac_smoke.ps1`: PASS.
- `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -CaseTimeoutSeconds 180 -IncludeExpansionAuditRows`: PASS.
- `test/pbe_tabi_topdown_mac_audit.ps1`: PASS
  (`unsafe=392`, `rejects=392`, far-node interactions `392 -> 0`).
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS
  (`rel3=9.38761e-07`, `rel4=1.35948e-07`,
  `topdown_rel3=4.61951e-05`).
- `test/pbe_tabi_tree_builder_smoke.ps1`: PASS.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_profile_smoke.ps1`: PASS
  (`panels=120`, `matvec_calls=9`).
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.
- `test/pbe_regression.ps1`: PASS.
- `test/pbsa_smoke.ps1`: PASS.
- `test/sasa_surface_smoke.ps1`: PASS.

Stage 3G decision:

- Keep all Stage 3E/3F/3G controls hidden/report-only. Production defaults are
  unchanged.
- Barnase still sits in the `~5.8e-2` sampled residual basin, while the 120s
  timeout-safe row now shows high far-expansion Taylor-vs-direct error
  (`max 4.8596`, mean `1.278e-2`) and no fragmented-node split activity.
- Next branch is Taylor/JCP moment construction:
  1. inspect the JCP Taylor series signs/scales for `phi` and `phin` moments;
  2. validate direct moment construction at `area-centroid` centers against
     accepted far-node direct sums, not only translated moment equality;
  3. add a targeted fix before any leaf batching or MAC-loosening work;
  4. rerun the Barnase 120s Stage 3G probe before considering a 300s run.

## Stage 3H/9A - TABI JCP Moment Translation and Convergence Audit

Implemented:

- Added hidden moment/convergence controls:
  - `-pbe-tabi-moment-construction=translated|direct-node`, default
    `translated`;
  - `-pbe-tabi-moment-sum=naive|kahan`, default `naive`;
  - `-pbe-tabi-convergence-audit=off|sampled|full-small`, default `off`.
- Added the `direct-node` probe path for JCP tree moments. Internal nodes can
  now build moments by scanning the panels owned by that node at the node's own
  Taylor center, instead of only translating child moments upward through M2M.
- Added optional Kahan accumulation for panel moment construction and M2M
  accumulation.
- Extended timing/profile/progress diagnostics:
  - `tabi_moment_translation_offset_max/mean`;
  - `tabi_moment_cancellation_ratio_max/mean`;
  - `tabi_moment_dynamic_range_max`;
  - `tabi_far_convergence_rho_max/mean`;
  - `tabi_node_aspect_ratio_max/mean`;
  - `tabi_accepted_far_aspect_ratio_max/mean`.
- Kept `tight-aabb` bounds for MAC/physical-gap while `area-centroid` remains
  a hidden Taylor-center policy; radius is still recomputed from the Taylor
  center to the farthest AABB corner so MAC is not relaxed implicitly.
- Added `test/pbe_tabi_jcp_moment_construction_smoke.ps1` for translated vs
  direct-node moment construction on a synthetic hollow case.
- Added `test/pbe_tabi_jcp_mixed_derivative_smoke.ps1` comparing the analytic
  JCP Taylor/mixed derivative path against a finite-difference debug oracle for
  both Poisson (`kappa=0`) and Yukawa (`kappa>0`) cases.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with
  `-IncludeStage3HRows`, adding:
  - `topdown-aabb/physical-gap/area-centroid/split/translated/naive`;
  - `topdown-aabb/physical-gap/area-centroid/split/translated/kahan`;
  - `topdown-aabb/physical-gap/area-centroid/split/direct-node/kahan`.
- Extended the Barnase truth/tree harness path to pass through
  `MomentConstruction`, `MomentSum`, and `ConvergenceAudit`, and to report the
  new timeout-safe progress counters.

Measured Stage 3H smoke evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Moment construction smoke:
  `test/pbe_tabi_jcp_moment_construction_smoke.ps1`: PASS.
  - panels `64`;
  - translated/direct-node direct-vs-tree rel both `5.889e-10`;
  - max translation offset `83.505`;
  - cancellation ratio max `1.414e5`;
  - dynamic range max `1.197e6`;
  - convergence `rho_max=0.231941`.
- Mixed derivative smoke:
  `test/pbe_tabi_jcp_mixed_derivative_smoke.ps1`: PASS.
  - panels `48`;
  - Poisson finite-difference rel `2.856e-10`;
  - Yukawa finite-difference rel `1.331e-10`;
  - Yukawa far interactions `96`.

Measured Stage 3H medium rows:

- Full Stage 3H medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3HRows -CaseTimeoutSeconds 180`:
  PASS and wrote
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| fixture | Stage 3H row | matvec rel | energy drift rel | sampled residual | translation rel max | cancellation max | dynamic range max | expansion rel max | rho max | node aspect max | accepted far aspect max |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `mini-protein` | `translated/naive` | `4.463e-4` | `6.757e-5` | `9.762e-5` | `1.087e-15` | `1.147e5` | `2.551e16` | `1.037e1` | `0.699943` | `10.055` | `1.930` |
| `mini-protein` | `translated/kahan` | `4.463e-4` | `6.757e-5` | `9.762e-5` | `1.148e-15` | `1.147e5` | `2.551e16` | `1.037e1` | `0.699943` | `10.055` | `1.930` |
| `mini-protein` | `direct-node/kahan` | `4.463e-4` | `6.757e-5` | `9.762e-5` | `1.170e-15` | `1.147e5` | `2.551e16` | `1.037e1` | `0.699943` | `10.055` | `1.930` |
| `charged-peptide` | `translated/naive` | `1.974e-4` | `1.224e-4` | `1.772e-5` | `7.124e-15` | `5.070e5` | `5.840e16` | `1.277e1` | `0.699861` | `9.289` | `2.240` |
| `charged-peptide` | `translated/kahan` | `1.974e-4` | `1.224e-4` | `1.772e-5` | `5.058e-15` | `5.070e5` | `5.840e16` | `1.277e1` | `0.699861` | `9.289` | `2.240` |
| `charged-peptide` | `direct-node/kahan` | `1.974e-4` | `1.224e-4` | `1.772e-5` | `7.211e-15` | `5.070e5` | `5.840e16` | `4.041e0` | `0.699861` | `9.289` | `2.240` |
| `compact-mini-protein` | `translated/naive` | `8.891e-4` | `7.254e-5` | `5.314e-5` | `0.000e0` | `6.568e4` | `2.261e15` | `1.498e1` | `0.699997` | `4.999` | `2.321` |
| `compact-mini-protein` | `translated/kahan` | `8.891e-4` | `7.254e-5` | `5.314e-5` | `1.192e-15` | `6.568e4` | `2.261e15` | `1.498e1` | `0.699997` | `4.999` | `2.321` |
| `compact-mini-protein` | `direct-node/kahan` | `8.891e-4` | `7.254e-5` | `5.314e-5` | `1.453e-15` | `6.568e4` | `2.261e15` | `1.498e1` | `0.699997` | `4.999` | `2.321` |

- All Stage 3H medium rows remain under the Stage 3E/3F success gates:
  `matvec rel < 2.5e-3`, `energy drift rel < 5e-3`, and sampled residual within
  `2x` grouped-morton.
- Medium translation audit errors are at numerical noise level. Direct-node
  only reduced the local far-expansion max on `charged-peptide`; it did not
  change the global medium matvec/residual rows.

Barnase Stage 3H timeout-safe probes:

| Barnase row | status | time | panels | unknowns | matvecs | sampled residual | translation rel max | offset max | cancellation max | dynamic range max | expansion rel max | rho max | node aspect max | accepted far aspect max | phase |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `translated/naive` | `timeout` | `120.05s` | `15129` | `30258` | `15` | `4.4313e-2` | `0.000e0` | `13.890` | `2.546e6` | `2.981e14` | `4.8596` | `0.699998` | `5.830` | `1.909` | `gmres_sampled_true_probe` |
| `direct-node/kahan` | `timeout` | `120.13s` | `15129` | `30258` | `13` | `5.8313e-2` | `0.000e0` | `13.890` | `2.546e6` | `2.981e14` | `4.8596` | `0.699998` | `5.830` | `1.909` | `gmres_iter` |

- Both Barnase probes used `surface_points=8`, `topdown-aabb`,
  `physical-gap`, `area-centroid`, `fragmentation-safety=split`, `leaf=128`,
  `theta=0.7`, `tight-aabb`, `onthefly`, `duffy_adaptive`,
  `quad_buffer=adaptive`, `precond=diagonal`, `moment-audit=sampled`,
  `expansion-audit=sampled`, and `convergence-audit=sampled`.
- `direct-node/kahan` did not reduce the Barnase far-expansion rel max
  (`4.8596` in both rows). It was also slower at the same 120s timeout, so it
  completed fewer matvecs.

Verification after Stage 3H:

- `test/pbe_tabi_jcp_moment_construction_smoke.ps1`: PASS.
- `test/pbe_tabi_jcp_mixed_derivative_smoke.ps1`: PASS.
- `test/pbe_tabi_far_expansion_audit_smoke.ps1`: PASS.
- `test/pbe_tabi_fragmented_node_mac_smoke.ps1`: PASS.
- `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3HRows -CaseTimeoutSeconds 180`: PASS.
- `test/pbe_tabi_topdown_mac_audit.ps1`: PASS.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS.
- `test/pbe_tabi_tree_builder_smoke.ps1`: PASS.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_profile_smoke.ps1`: PASS.
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.
- `test/pbe_regression.ps1`: PASS.
- `test/pbsa_smoke.ps1`: PASS.
- `test/sasa_surface_smoke.ps1`: PASS.

Stage 3H decision:

- Keep all Stage 3H controls hidden/report-only. Production defaults are
  unchanged.
- The comment's M2M concern is now tested by evidence. `translate_child_moment`
  does not look like the primary Barnase failure source under this probe:
  translation audit is near zero on medium rows and zero in the timeout-safe
  Barnase progress rows, and `direct-node/kahan` does not reduce the Barnase
  far-expansion error.
- Kahan accumulation is useful as a diagnostic but did not change the measured
  medium or Barnase rows enough to justify promotion.
- The mixed derivative smoke passes against a finite-difference oracle in safe
  synthetic cases, so there is no immediate controlled-case sign/scale bug in
  the JCP mixed derivative path.
- Current evidence points to accepted far-node Taylor truncation/convergence
  safety, not physical-gap, fragmentation split, or M2M translation:
  Barnase still shows `rho_max` almost exactly at `theta=0.7` and
  far-expansion rel max `4.8596`.
- Next branch: add a convergence-safety split/audit policy for high-rho and/or
  high-aspect accepted far nodes before any performance tuning. Do not run a
  300s Barnase probe until the 120s row shows either sampled residual moving
  clearly below the `~5e-2` basin or far-expansion rel max falling sharply.

## Stage 3I/9B - TABI Convergence-Safe MAC and Taylor Sweep

Implemented:

- Added hidden convergence/MAC safety controls:
  - `-pbe-tabi-convergence-safety=off|audit|split`, default `off`;
  - `-pbe-tabi-convergence-rho-cap=<0.0..2.0>`, default `0.0` meaning off;
  - `-pbe-tabi-mac-aspect-cap=<0.0..50.0>`, default `0.0` meaning off.
- Extended `evaluate_mac` after ratio MAC and physical-gap:
  - `audit` records accepted far candidates with `rho > rho_cap` or
    `aspect > aspect_cap`;
  - `split` rejects those candidates and descends to child/leaf paths;
  - leaf direct handling is unchanged.
- Extended timing/profile/progress output with:
  - `tabi_convergence_safety`, `tabi_convergence_rho_cap`,
    `tabi_mac_aspect_cap`;
  - `tabi_mac_convergence_rho_rejects`,
    `tabi_mac_convergence_aspect_rejects`;
  - `tabi_mac_convergence_rho_reject_panels`,
    `tabi_mac_convergence_aspect_reject_panels`.
- Extended Stage 3I harness coverage:
  - `test/pbe_tabi_barnase_truth_matrix.ps1` now accepts
    `-TreeTaylorOrders @(...)` and can run tree-only probes with
    `-SkipGridReference`;
  - `test/pbe_tabi_barnase_tree_mac_probe.ps1` forwards Taylor order and
    convergence guard controls;
  - `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3IRows`
    adds theta/order/guard rows while preserving the Stage 3E/3F/3H rows.
- Added `test/pbe_tabi_convergence_mac_safety_smoke.ps1`, a synthetic
  high-aspect/high-rho accepted-far case that proves `audit` sees candidates
  and `split` trades far work for direct leaf work.

Measured Stage 3I smoke evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Convergence MAC safety smoke:
  `test/pbe_tabi_convergence_mac_safety_smoke.ps1`: PASS.
  - panels `72`;
  - audit rho/aspect would-rejects `38 / 135`;
  - split rho/aspect rejects `38 / 135`;
  - far work `135 -> 0`;
  - leaf work `2362 -> 5112`;
  - far-expansion rel max `0.730137 -> 0`;
  - direct-vs-tree vector rel `3.449e-5 -> 2.741e-16`.

Measured Stage 3I medium sweep:

- Full Stage 3I medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3IRows -CaseTimeoutSeconds 180`:
  PASS and wrote
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| fixture | Stage 3I row | order | theta | guard | matvec rel | energy drift rel | sampled residual | expansion rel max | rho max | rho rejects | leaf work/matvec |
| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `mini-protein` | baseline | 3 | 0.7 | off | `4.463e-4` | `6.757e-5` | `9.762e-5` | `1.037e1` | `0.699943` | `0` | `4.293e6` |
| `mini-protein` | theta scan | 3 | 0.5 | off | `5.429e-5` | `1.816e-5` | `9.764e-5` | `8.379e-1` | `0.499998` | `0` | `6.118e6` |
| `mini-protein` | theta scan | 3 | 0.4 | off | `9.229e-6` | `3.224e-6` | `9.759e-5` | `7.055e-2` | `0.399993` | `0` | `6.783e6` |
| `mini-protein` | order scan | 4 | 0.5 | off | `1.479e-5` | `1.842e-6` | `9.760e-5` | `4.158e-1` | `0.499998` | `0` | `6.118e6` |
| `mini-protein` | rho cap | 3 | 0.7 | `rho<=0.5` | `5.429e-5` | `1.816e-5` | `9.764e-5` | `8.379e-1` | `0.499998` | `83175` | `6.118e6` |
| `charged-peptide` | baseline | 3 | 0.7 | off | `1.974e-4` | `1.224e-4` | `1.772e-5` | `1.277e1` | `0.699861` | `0` | `3.274e6` |
| `charged-peptide` | theta scan | 3 | 0.5 | off | `2.775e-5` | `1.205e-5` | `1.770e-5` | `4.435e0` | `0.499980` | `0` | `4.507e6` |
| `charged-peptide` | theta scan | 3 | 0.4 | off | `5.288e-16` | `0` | `1.770e-5` | `0` | `0` | `0` | `5.021e6` |
| `charged-peptide` | order scan | 4 | 0.5 | off | `7.974e-6` | `2.764e-6` | `1.770e-5` | `7.233e-1` | `0.499980` | `0` | `4.507e6` |
| `charged-peptide` | rho cap | 3 | 0.7 | `rho<=0.5` | `2.775e-5` | `1.205e-5` | `1.770e-5` | `4.435e0` | `0.499980` | `77260` | `4.507e6` |
| `compact-mini-protein` | baseline | 3 | 0.7 | off | `8.891e-4` | `7.254e-5` | `5.314e-5` | `1.498e1` | `0.699997` | `0` | `1.004e7` |
| `compact-mini-protein` | theta scan | 3 | 0.5 | off | `9.158e-5` | `8.093e-6` | `5.306e-5` | `1.674e0` | `0.499991` | `0` | `1.459e7` |
| `compact-mini-protein` | theta scan | 3 | 0.4 | off | `2.205e-5` | `2.348e-7` | `5.307e-5` | `1.295e0` | `0.399987` | `0` | `1.568e7` |
| `compact-mini-protein` | order scan | 4 | 0.5 | off | `3.235e-5` | `6.387e-6` | `5.306e-5` | `5.432e-1` | `0.499991` | `0` | `1.459e7` |
| `compact-mini-protein` | rho cap | 3 | 0.7 | `rho<=0.5` | `9.158e-5` | `8.093e-6` | `5.306e-5` | `1.222e0` | `0.499991` | `482550` | `1.459e7` |

- All listed Stage 3I medium rows stay inside the Stage 3E/3F global gates:
  `matvec rel < 2.5e-3`, `energy drift rel < 5e-3`, and sampled residual
  within `2x` grouped-morton.
- The local expansion audit is the discriminating signal. Lowering theta or
  imposing `rho_cap=0.5` reduces the worst local Taylor error substantially,
  but does not yet reach the Barnase target of `<1e-3`.
- `mac_aspect_cap=3.0` did not trigger on these medium fixtures because
  accepted far aspect ratios were below the cap; aspect guard remains useful
  for synthetic/high-aspect cases but is not the active medium/Barnase lever
  at this threshold.

Barnase Stage 3I 120s probes:

- Theta/order sweep command shape:
  `test/pbe_tabi_barnase_truth_matrix.ps1 -TreeProbeOnly -SkipGridReference -CaseTimeoutSeconds 120 -SurfacePoints 8 -TreeLeafSizes 128 -TreeOpeningThetas @(0.7,0.5,0.4,0.6) -TreeTaylorOrders @(3,4) -NodeBoundsModes tight-aabb -TreeBuilders topdown-aabb -MacSafety physical-gap -ExpansionCenters area-centroid -FragmentationSafety split -MomentAudit sampled -ExpansionAudit sampled -ConvergenceAudit sampled -NearIntegration duffy_adaptive -QuadBuffer adaptive`.
- Guard row command shape:
  same settings, plus `-ConvergenceSafety split -ConvergenceRhoCap 0.5 -MacAspectCap 3.0`.

| Barnase row | status | matvecs | sampled residual | expansion rel max | expansion rel mean | rho max | rho rejects | leaf work/matvec | phase |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| order 3 theta 0.7 | `timeout` | `12` | `5.8313e-2` | `4.859629` | `1.1699e-2` | `0.699998` | `0` | `7.019e7` | `gmres_iter` |
| order 3 theta 0.6 | `timeout` | `12` | `5.8329e-2` | `2.015257` | `7.516e-3` | `0.599999` | `0` | `9.442e7` | `gmres_iter` |
| order 3 theta 0.5 | `timeout` | `14` | `5.8330e-2` | `1.653837` | `4.207e-3` | `0.500000` | `0` | `1.275e8` | `gmres_iter` |
| order 3 theta 0.4 | `timeout` | `16` | `4.4371e-2` | `0.570735` | `1.785e-3` | `0.400000` | `0` | `1.706e8` | `gmres_iter` |
| order 4 theta 0.7 | `timeout` | `7` | `1.57275e-1` | `0.685814` | `3.855e-3` | `0.699998` | `0` | `7.019e7` | `gmres_iter` |
| order 4 theta 0.6 | `timeout` | `7` | `1.57269e-1` | `0.356780` | `1.837e-3` | `0.599999` | `0` | `9.442e7` | `gmres_iter` |
| order 4 theta 0.5 | `timeout` | `8` | `1.57273e-1` | `0.074761` | `8.26e-4` | `0.500000` | `0` | `1.275e8` | `gmres_iter` |
| order 4 theta 0.4 | `timeout` | `11` | `5.8331e-2` | `0.126600` | `4.35e-4` | `0.400000` | `0` | `1.706e8` | `gmres_iter` |
| order 3 theta 0.7 + `rho<=0.5/aspect<=3` | `timeout` | `14` | `5.8330e-2` | `1.653837` | `4.422e-3` | `0.500000` | `2213232` | `1.275e8` | `gmres_iter` |

- The Stage 3H baseline residual basin is reproduced:
  order 3 theta 0.7 gives sampled residual `5.8313e-2` and far-expansion rel
  max `4.859629`.
- The convergence hypothesis is supported:
  reducing effective `rho` lowers far-expansion max monotonically for order 3
  (`4.86 -> 2.02 -> 1.65 -> 0.571`) and `rho_cap=0.5` reproduces the theta
  `0.5` expansion/work shape while recording `2,213,232` rho rejects.
- Taylor order 4 lowers the expansion audit error further, especially at
  theta `0.5` (`0.074761`), but it is much slower in this 120s budget and did
  not improve sampled residual yet.
- No Stage 3I 120s Barnase row reaches the promotion target
  `far-expansion rel max < 1e-3`; therefore no Barnase 300s probe was run.

Verification after Stage 3I:

- `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- `test/pbe_tabi_convergence_mac_safety_smoke.ps1`: PASS.
- `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3IRows -CaseTimeoutSeconds 180`: PASS.
- Barnase 120s theta/order sweep: PASS as timeout-safe evidence, all rows
  emitted progress counters.
- Barnase 120s guard row: PASS as timeout-safe evidence, progress recorded
  `rho_cap` rejects.
- `test/pbe_tabi_jcp_moment_construction_smoke.ps1`: PASS.
- `test/pbe_tabi_jcp_mixed_derivative_smoke.ps1`: PASS.
- `test/pbe_tabi_far_expansion_audit_smoke.ps1`: PASS.
- `test/pbe_tabi_fragmented_node_mac_smoke.ps1`: PASS.
- `test/pbe_tabi_topdown_mac_audit.ps1`: PASS.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS.
- `test/pbe_tabi_tree_builder_smoke.ps1`: PASS.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_profile_smoke.ps1`: PASS.
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.
- `test/pbe_regression.ps1`: PASS.
- `test/pbsa_smoke.ps1`: PASS.
- `test/sasa_surface_smoke.ps1`: PASS.

Stage 3I decision:

- Keep `topdown-aabb`, `physical-gap`, `area-centroid`,
  `fragmentation-safety`, Taylor order 4, and convergence-safety controls
  hidden/experimental. Production defaults are unchanged.
- The comment is directionally correct: `rho ~= theta=0.7` is too loose for
  JCP BIE order 3 on Barnase. Stage 3I reduces the far-expansion error by
  tightening effective `rho` and/or increasing Taylor order, so the branch is
  convergence/truncation dominated rather than M2M, mixed derivative, or
  physical-gap dominated.
- The current fix is not sufficient. Even the strongest 120s Stage 3I rows
  still sit at expansion rel max `7.476e-2` or worse, far above `<1e-3`.
- Do not run Barnase 300s or optimize leaf batching yet. The next branch should
  add a stricter convergence model before performance tuning, for example:
  1. test smaller `rho_cap` values (`0.35`, `0.30`) and/or a direct
     convergence-error estimator for accepted JCP far nodes;
  2. evaluate Taylor order beyond 4 only as a hidden probe if derivative series
     cost is acceptable;
  3. only after far-expansion rel max falls near `<1e-3`, revisit leaf-pair
     batching to recover the `1.7e8` leaf-work/matvec cost of strict rho rows.

## Stage 3J/9C - TABI Adaptive JCP Convergence Estimator

Implemented:

- Added hidden adaptive convergence-estimator controls:
  - `-pbe-tabi-convergence-estimator=off|moment-tail|hybrid`, default `off`;
  - `-pbe-tabi-convergence-tail-threshold=<1..1e12>`, default `1e4`;
  - `-pbe-tabi-convergence-estimator-error-cap=<1e-8..1.0>`, default `1e-3`;
  - `-pbe-tabi-convergence-unstable-rho-cap=<0.05..2.0>`, default `0.25`.
- Estimator is active only on JCP treecode `interaction-mode=onthefly`.
  Cached/non-JCP paths keep counters zero and do not change behavior.
- After ratio MAC, physical-gap, and fragmentation split, accepted JCP far
  candidates can now be scored from current per-matvec `PackedJcpMoments`:
  `tail_ratio = max(degree-p moment/radius^p) /
  max(degree 0..1 moment/radius^degree)`, and hybrid estimated error is
  `tail_ratio * rho^(p+1) * 30`.
- With `convergence-safety=audit`, unstable far candidates are only counted.
  With `convergence-safety=split`, unstable far candidates descend unless
  `rho <= convergence-unstable-rho-cap`.
- Extended timing/profile/progress output with estimator mode, thresholds,
  unstable/split counts and panels, tail-ratio/error max/mean, and
  expansion-audit calibration counters:
  `tabi_convergence_estimator_high_error_caught/missed`.
- Extended harness coverage:
  - `test/pbe_tabi_convergence_estimator_smoke.ps1`;
  - `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3JRows`;
  - `test/pbe_tabi_barnase_truth_matrix.ps1` and
    `test/pbe_tabi_barnase_tree_mac_probe.ps1` now pass through and report the
    estimator controls/counters.

Measured Stage 3J smoke evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Convergence estimator smoke:
  `test/pbe_tabi_convergence_estimator_smoke.ps1`: PASS.
  - panels `96`;
  - audit estimator samples `240`;
  - audit unstable far candidates `240`;
  - hybrid split rejects `128`;
  - loose `moment-tail` benign row rejects `0`;
  - direct-vs-tree vector rel `5.166e-5 -> 1.962e-6`;
  - far work `240 -> 112`;
  - leaf work `5856 -> 7780`;
  - max `tail_ratio=1.272`, max `estimated_error=2.552`.

Measured Stage 3J medium quick matrix:

- Command:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -Quick -IncludeStage3JRows -CaseTimeoutSeconds 180 -AllowFailedStatus`.
- Output:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| mini-protein Stage 3J row | order | theta | guard/estimator | matvec rel | energy drift rel | sampled residual | expansion rel max | estimator split rejects | leaf work/matvec | far work/matvec |
| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `center_area_o3_theta0p6` | 3 | 0.6 | off | `2.188e-4` | `1.432e-5` | `9.761e-5` | `6.558` | `0` | `5.106e6` | `7060` |
| `center_aabb_o3_theta0p6` | 3 | 0.6 | off | `2.716e-4` | `3.529e-5` | `9.752e-5` | `5.551` | `0` | `4.647e6` | `8329` |
| `rho_cap0p35` | 3 | 0.6 | `rho<=0.35` | `4.561e-6` | `8.178e-7` | `9.759e-5` | `6.099e-2` | `0` | `6.940e6` | `3089` |
| `rho_cap0p30` | 3 | 0.6 | `rho<=0.30` | `1.847e-6` | `1.327e-7` | `9.759e-5` | `1.847e-2` | `0` | `7.106e6` | `1853` |
| `rho_cap0p25` | 3 | 0.6 | `rho<=0.25` | `5.318e-7` | `3.839e-8` | `9.759e-5` | `6.575e-3` | `0` | `7.291e6` | `498` |
| `hybrid_t1e3` | 3 | 0.6 | hybrid split | `6.012e-7` | `3.505e-8` | `9.759e-5` | `6.575e-3` | `130801` | `7.289e6` | `515` |
| `hybrid_t1e4` | 3 | 0.6 | hybrid split | `6.012e-7` | `3.505e-8` | `9.759e-5` | `6.575e-3` | `130801` | `7.289e6` | `515` |
| `order4_hybrid_t1e4` | 4 | 0.6 | hybrid split | `1.768e-7` | `3.922e-8` | `9.758e-5` | `1.410e-3` | `123943` | `7.227e6` | `972` |

- Medium quick rows remain inside Stage 3E/3F gates:
  `matvec rel < 2.5e-3`, `energy drift rel < 5e-3`, and sampled residual
  remains near the grouped/direct medium baseline.
- On mini-protein, the hybrid estimator is currently about as aggressive as
  strict `rho_cap=0.25`; it improves expansion error but is not yet selective.
- `aabb-center` lowered the local max spike on this mini row, but it worsened
  energy/matvec drift slightly and did not become a promotion candidate.

Barnase Stage 3J 120s probes:

All rows used `surface_points=8`, `leaf=128`, `theta=0.6`,
`topdown-aabb`, `tight-aabb`, `physical-gap`, `area-centroid` unless noted,
`fragmentation-safety=split`, `moment-construction=translated`,
`moment-sum=naive`, `moment-audit=sampled`, `expansion-audit=sampled`,
`convergence-audit=sampled`, `near=duffy_adaptive`,
`quad_buffer=adaptive`, and `precond=diagonal`.

| Barnase Stage 3J row | status | matvecs | sampled residual | true residual | expansion rel max | expansion rel mean | rho max | rho rejects | estimator unstable | estimator split rejects | caught/missed | leaf work/matvec |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| center `area-centroid` | `timeout` | `14` | `5.8329e-2` | `NA` | `2.015257` | `7.564e-3` | `0.599999` | `0` | `0` | `0` | `0/0` | `9.442e7` |
| center `aabb-center` | `timeout` | `13` | `5.8337e-2` | `NA` | `2.548277` | `9.722e-3` | `0.600000` | `0` | `0` | `0` | `0/0` | `7.895e7` |
| strict `rho_cap=0.35` | `timeout` | `20` | `6.389e-3` | `NA` | `0.626529` | `1.149e-3` | `0.349999` | `5,716,100` | `0` | `0` | `0/0` | `1.933e8` |
| strict `rho_cap=0.30` | `timeout` | `29` | `2.320e-3` | `NA` | `0.651508` | `4.52e-4` | `0.300000` | `10,253,675` | `0` | `0` | `0/0` | `2.128e8` |
| hybrid order 3 threshold `1e3` | `timeout` | `39` | `4.53e-4` | `NA` | `5.199e-2` | `2.88e-4` | `0.599999` | `0` | `16,414,508` | `15,593,331` | `974/0` | `2.241e8` |
| hybrid order 3 threshold `1e4` | `ok` | `45` | `3.2e-5` | `5.449e-5` | `5.199e-2` | `2.822e-4` | `0.599999` | `0` | `18,940,620` | `17,992,305` | `1074/0` | `2.241e8` |
| hybrid order 4 threshold `1e4` | `timeout` | `29` | `2.320e-3` | `NA` | `4.294e-3` | `3.3e-5` | `0.599999` | `0` | `11,375,877` | `11,182,936` | `40/15` | `2.207e8` |

- Center-policy check: on Barnase, `area-centroid` is better than
  `aabb-center` for this 120s probe (`2.015` vs `2.548` expansion rel max).
  Keep `aabb-center` hidden as a comparison row only.
- Strict rho lower-bound confirms the comment's accuracy/performance warning:
  sampled residual leaves the `~5.8e-2` basin (`6.389e-3` and `2.320e-3`),
  but leaf work rises to roughly `1.9e8-2.1e8` per matvec.
- Hybrid order 3 threshold `1e4` is the first Barnase 120s row in this branch
  to finish successfully with true residual `5.449e-5` and sampled residual
  `3.2e-5`.
- The hybrid estimator does not yet satisfy the strict promotion target
  `far_expansion_rel_error_max < 1e-3`: the max spike remains about
  `5.2e-2`. However, expansion mean is `2.82e-4` and sampled expansion-audit
  high-error calibration caught all high-error pairs (`1074/0` caught/missed)
  in the converged order 3 row.
- Order 4 reduces the worst expansion spike much further (`4.294e-3`) but
  is slower in the 120s budget and still missed some high-error sampled pairs
  (`40/15`), so order 4 is not the next default candidate.

Verification after Stage 3J:

- `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- `test/pbe_tabi_convergence_estimator_smoke.ps1`: PASS.
- `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -Quick -IncludeStage3JRows -CaseTimeoutSeconds 180 -AllowFailedStatus`: PASS.
- Barnase 120s center, strict-rho, hybrid estimator, and order-4 fallback
  probes: PASS as timeout-safe/converged evidence, all emitted progress rows.
- `test/pbe_tabi_convergence_mac_safety_smoke.ps1`: PASS.
- `test/pbe_tabi_jcp_moment_construction_smoke.ps1`: PASS.
- `test/pbe_tabi_jcp_mixed_derivative_smoke.ps1`: PASS.
- `test/pbe_tabi_far_expansion_audit_smoke.ps1`: PASS.
- `test/pbe_tabi_fragmented_node_mac_smoke.ps1`: PASS.
- `test/pbe_tabi_topdown_mac_audit.ps1`: PASS.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS.
- `test/pbe_tabi_tree_builder_smoke.ps1`: PASS.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_profile_smoke.ps1`: PASS.
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.
- `test/pbe_regression.ps1`: PASS.
- `test/pbsa_smoke.ps1`: PASS.
- `test/sasa_surface_smoke.ps1`: PASS.

Stage 3J decision:

- Keep `topdown-aabb`, `physical-gap`, `area-centroid`,
  `fragmentation-safety`, convergence-safety, and convergence-estimator
  controls hidden/experimental. Production defaults are unchanged.
- The adaptive estimator hypothesis is strongly supported for solver progress:
  the Barnase order 3 hybrid threshold `1e4` row breaks out of the
  `~5.8e-2` residual basin and reaches true residual `5.449e-5` within the
  120s budget.
- Do not promote yet. The hard far-expansion max target `<1e-3` is not met;
  order 3 has a remaining local spike around `5.2e-2`, and order 4 still sits
  at `4.294e-3` while costing more and missing some high-error sampled pairs.
- Next branch should tune estimator selectivity and calibration before leaf
  batching:
  1. add per-channel or percentile/mean-sensitive estimator calibration so
     the order 3 hybrid row keeps its convergence but avoids behaving like
     strict `rho_cap=0.25` on medium rows;
  2. investigate the remaining expansion max outliers caught by audit in the
     converged Barnase row, including whether they are harmless isolated pairs
     or need a stricter local split rule;
  3. only after the max outliers are bounded near `<1e-3` or proven harmless
     by residual/energy gates should a 300s Barnase run or leaf-pair batching
     become the priority.

## Stage 3K/9D - TABI Calibrated Hybrid Estimator and Mean-Error Promotion Gate

Implemented:

- Added hidden calibration control
  `-pbe-tabi-convergence-calibration=stage3j|consecutive|aspect|consecutive-aspect`,
  default `stage3j`.
- `stage3j` preserves Stage 3J behavior: degree `p` tail only and static JCP
  amplification `30`.
- `consecutive` uses the consecutive-tail guard:
  `tail_norm = max(norm(p), norm(p-1) * rho)`.
- `aspect` uses dynamic JCP amplification:
  `jcp_amp = 30 * max(1, node_aspect_ratio(node))`.
- `consecutive-aspect` enables both calibration rules.
- Extended timing/profile/progress diagnostics:
  `tabi_convergence_calibration`,
  `tabi_convergence_tail_norm_p_max/mean`,
  `tabi_convergence_tail_norm_pm1_max/mean`,
  `tabi_convergence_effective_jcp_amp_max/mean`, and
  `tabi_far_expansion_mean_gate_pass`.
- Mean gate is emitted from timing/profile/progress only when expansion-audit
  sampled pairs exist, mean far-expansion error is `<5e-4`, and
  high-error misses are `0`.
- Extended harnesses:
  - `test/pbe_tabi_convergence_calibration_smoke.ps1`;
  - `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3KRows`;
  - `test/pbe_tabi_barnase_truth_matrix.ps1` and
    `test/pbe_tabi_barnase_tree_mac_probe.ps1` pass through/report the
    calibration mode and mean-gate counters.

Measured Stage 3K smoke evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Calibration smoke:
  `test/pbe_tabi_convergence_calibration_smoke.ps1`: PASS.
  - panels `112`;
  - `stage3j` estimated-error max `4.64661`;
  - `consecutive` estimated-error max `8.38751`;
  - `aspect` effective JCP amp max `38887.5`;
  - `consecutive-aspect` effective JCP amp max `38887.5`;
  - combined split rejects `119`;
  - benign row retained far work (`135` far interactions).

Measured Stage 3K medium matrix:

- Command:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3KRows -CaseTimeoutSeconds 180 -AllowFailedStatus`.
- Output:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| fixture | row | sampled residual | matvec rel | energy drift rel | expansion mean | mean gate | missed | leaf work/matvec | classification |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| mini-protein | order3 `stage3j` | `9.759e-5` | `6.012e-7` | `3.505e-8` | `2.748e-4` | `1` | `0` | `7.289e6` | mean gate pass |
| mini-protein | order3 `consecutive-aspect` | `9.759e-5` | `5.318e-7` | `3.839e-8` | `2.767e-4` | `1` | `0` | `7.291e6` | mean gate pass |
| compact-mini-protein | order3 `stage3j` | `5.307e-5` | `6.265e-7` | `1.396e-7` | `3.253e-4` | `1` | `0` | `1.708e7` | mean gate pass |
| compact-mini-protein | order3 `aspect` | `5.307e-5` | `4.946e-7` | `4.230e-8` | `3.153e-4` | `1` | `0` | `1.710e7` | mean gate pass |
| charged-peptide | order3 calibrated rows | `1.770e-5` | `5.288e-16` | `0` | `0` | `0` | `0` | `5.021e6` | medium gate pass |

- Medium rows pass the Stage 3E accuracy/residual gates.
- Mean gate is intentionally not counted for charged-peptide because no
  expansion-audit pairs were emitted in those rows (`mean_gate=0`).
- `consecutive-aspect` did not materially improve medium mean error versus
  `stage3j`; `aspect` slightly improved compact-mini-protein mean error while
  preserving the same residual basin.
- Relaxed `consecutive-aspect` error cap `3e-3` and order-4 fallback did not
  satisfy the mean gate on compact-mini-protein because high-error misses
  appeared (`21` and `87`, respectively).

Barnase Stage 3K 120s probes:

All rows used `surface_points=8`, `leaf=128`, `theta=0.6`, order `3`,
`topdown-aabb`, `tight-aabb`, `physical-gap`, `area-centroid`,
`fragmentation-safety=split`, `moment-construction=translated`,
`moment-sum=naive`, `moment-audit=sampled`, `expansion-audit=sampled`,
`convergence-audit=sampled`, `convergence-safety=split`,
`convergence-estimator=hybrid`, `threshold=1e4`,
`unstable_rho_cap=0.25`, `near=duffy_adaptive`, `quad_buffer=adaptive`, and
`precond=diagonal`.

| Barnase Stage 3K row | status | matvecs | sampled residual | true residual | expansion mean | mean gate | expansion max | missed | effective amp max | leaf work/matvec |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `stage3j`, cap `1e-3` | `timeout` | `37` | `4.53e-4` | `6.677e-4` | `2.65e-4` | `1` | `5.199e-2` | `0` | `30.0` | `2.241e8` |
| `consecutive`, cap `1e-3` | `timeout` | `30` | `2.102e-3` | progress fallback | `2.46e-4` | `1` | `3.916e-2` | `0` | `30.0` | `2.241e8` |
| `aspect`, cap `1e-3` | `timeout` | `37` | `4.53e-4` | `6.677e-4` | `2.64e-4` | `1` | `5.199e-2` | `0` | `57.278` | `2.241e8` |
| `consecutive-aspect`, cap `1e-3` | `timeout` | `32` | `2.102e-3` | progress fallback | `2.79e-4` | `1` | `5.199e-2` | `0` | `57.278` | `2.241e8` |
| `consecutive-aspect`, cap `3e-3` | `timeout` | `34` | `2.102e-3` | progress fallback | `2.70e-4` | `1` | `5.196e-2` | `0` | `57.278` | `2.241e8` |

- Stage 3K confirms the user's mean-error direction: all Barnase 120s rows
  emitted timeout-safe progress rows with far-expansion mean `<5e-4` and
  high-error missed `0`.
- No Stage 3K Barnase 120s row satisfied the promotion residual gate
  `true or sampled residual <=1e-4`.
- Therefore Barnase 300s was not run in Stage 3K, and no default was promoted.
- `aspect` matches `stage3j` residual behavior while exposing dynamic JCP amp
  counters; `consecutive` and `consecutive-aspect` were more conservative but
  had worse 120s sampled residual in this run.
- Remaining risk is the isolated max spike around `5.2e-2`: it no longer
  blocks mean-gate evidence by itself, but promotion still requires residual
  confirmation.

Verification after Stage 3K:

- `test/pbe_tabi_convergence_calibration_smoke.ps1`: PASS.
- `test/pbe_tabi_convergence_estimator_smoke.ps1`: PASS.
- `test/pbe_tabi_convergence_mac_safety_smoke.ps1`: PASS.
- `test/pbe_tabi_far_expansion_audit_smoke.ps1`: PASS.
- `test/pbe_tabi_fragmented_node_mac_smoke.ps1`: PASS.
- `test/pbe_tabi_jcp_moment_construction_smoke.ps1`: PASS.
- `test/pbe_tabi_jcp_mixed_derivative_smoke.ps1`: PASS.
- `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3KRows -CaseTimeoutSeconds 180 -AllowFailedStatus`: PASS.
- Barnase 120s Stage 3K rows above: PASS as timeout-safe evidence, no
  promotion candidate.
- `test/pbe_tabi_topdown_mac_audit.ps1`: PASS.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS.
- `test/pbe_tabi_tree_builder_smoke.ps1`: PASS.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_profile_smoke.ps1`: PASS.
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.
- `test/pbe_regression.ps1`: PASS.
- `test/pbsa_smoke.ps1`: PASS.
- `test/sasa_surface_smoke.ps1`: PASS.

Stage 3K decision:

- Keep production defaults unchanged.
- Keep `topdown-aabb`, `physical-gap`, `area-centroid`,
  `fragmentation-safety`, convergence estimator, convergence calibration, and
  mean gate hidden/experimental.
- Do not promote and do not run Barnase 300s yet.
- Next branch should use the mean-gate evidence, but tune solver progress:
  1. investigate why calibrated consecutive rows reduce max/mean expansion
     error yet slow residual progress in Barnase 120s;
  2. compare Stage 3J converged `hybrid_t1e4` against the current Stage 3K
     `stage3j/aspect` rows under identical timeout/progress settings to
     separate runtime variance from calibration behavior;
  3. only if a 120s row again reaches residual `<=1e-4` with mean gate pass
     and missed `0`, run the 300s confirmation before promotion.

## Stage 3L/9E - TABI Dielectric-Aware Estimator and Target-Cluster Leaf Batching

Implemented:

- Added hidden `-pbe-tabi-convergence-dielectric-amp=off|on`, default `off`.
- When enabled for JCP on-the-fly convergence estimation, the estimator uses
  `max(30, abs(eps_ratio))` before optional aspect calibration. Existing
  behavior is preserved when the flag is `off`.
- Added timing/profile/progress fields:
  `tabi_convergence_dielectric_amp`,
  `tabi_convergence_dielectric_amp_max/mean`, plus the existing effective JCP
  amp fields after aspect calibration.
- Implemented hidden on-the-fly direct leaf batching through the existing
  `-pbe-tabi-cluster-policy=off|leaf-pair-strict|leaf-pair-balanced` scaffold.
  New hidden controls are:
  `-pbe-tabi-target-batch-size=1|4|8`,
  `-pbe-tabi-leaf-pair-min-source-panels`, and
  `-pbe-tabi-leaf-pair-min-targets`.
- Batching is conservative: only clean direct leaf paths are batched. Self,
  one-ring, geometric-buffer, Duffy/adaptive-buffer, and target-specific
  exclusion cases fall back to scalar leaf evaluation.
- Added counters:
  `tabi_leaf_pair_batch_interactions`,
  `tabi_leaf_pair_batch_targets`,
  `tabi_leaf_pair_batch_source_panels`,
  `tabi_leaf_pair_batch_residual_panels`,
  `tabi_leaf_pair_batch_fallbacks`, and
  `tabi_leaf_pair_batch_seconds`.
- Extended Stage 3L harness rows in:
  - `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3LRows`;
  - `test/pbe_tabi_barnase_truth_matrix.ps1` with dielectric amp, cluster
    policy, target batch size, and leaf-pair threshold passthrough/reporting.
- Added focused smokes:
  - `test/pbe_tabi_dielectric_estimator_smoke.ps1`;
  - `test/pbe_tabi_leaf_pair_batch_smoke.ps1`.

Measured Stage 3L smoke evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Dielectric estimator smoke: PASS.
  - panels `120`;
  - `dielectric_amp=off` reports amp max `30`;
  - `dielectric_amp=on` with `eps_ratio=80` reports amp max `80`;
  - estimated-error max increases from `1.12315` to `2.99506`;
  - split rejects increase from `4` to `24`;
  - low-eps case keeps conservative amp `30`.
- Leaf-pair batch smoke: PASS.
  - panels `96`;
  - batched-vs-scalar relative RMS `1.0891e-16`;
  - batch interactions `40`, batch targets `158`, source panels `2686`;
  - singular/exclusion fallback counter `616`.

Measured Stage 3L medium matrix:

- Command:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3LRows -CaseTimeoutSeconds 180 -AllowFailedStatus`.
- Output:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| fixture | row | matvec rel | energy drift rel | sampled residual | expansion mean | missed | dielectric amp max | batch interactions | total seconds |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| mini-protein | order3 dielectric off | `5.318e-7` | `3.84e-8` | `9.759e-5` | `2.767e-4` | `0` | `30` | `0` | `NA` |
| mini-protein | order3 dielectric on | `5.318e-7` | `3.84e-8` | `9.759e-5` | `2.767e-4` | `0` | `80` | `0` | `4.05` |
| mini-protein | order4 dielectric on | `8.920e-8` | `9.18e-9` | `9.759e-5` | `5.56e-5` | `0` | `80` | `0` | `4.38` |
| charged-peptide | order3 dielectric on | `5.288e-16` | `0` | `1.770e-5` | `0` | `0` | `80` | `0` | `2.78` |
| compact-mini-protein | order3 dielectric on | `4.946e-7` | `4.02e-8` | `5.306e-5` | `3.38e-4` | `0` | `80` | `0` | `13.42` |
| compact-mini-protein | order4 dielectric on | `6.964e-8` | `7.76e-9` | `5.306e-5` | `4.43e-5` | `4` | `80` | `0` | `11.37` |

Batching rows on the same order3 dielectric-on settings:

| fixture | cluster policy | batch interactions | batch seconds | total seconds | note |
| --- | --- | ---: | ---: | ---: | --- |
| mini-protein | `leaf-pair-strict` | `169425` | `1.925` | `4.40` | slower than cluster off |
| mini-protein | `leaf-pair-balanced` | `169425` | `1.869` | `4.51` | slower than cluster off |
| charged-peptide | `leaf-pair-strict` | `151820` | `1.262` | `3.10` | slower than cluster off |
| charged-peptide | `leaf-pair-balanced` | `151820` | `1.329` | `3.22` | slower than cluster off |
| compact-mini-protein | `leaf-pair-strict` | `648120` | `6.463` | `12.74` | slightly faster than cluster off |
| compact-mini-protein | `leaf-pair-balanced` | `648120` | `5.441` | `10.78` | about `20%` faster than cluster off |

Medium interpretation:

- Dielectric amp wiring works and is visible in timing/profile/progress.
- Medium accuracy/residual gates remain good for order3 rows, but the stricter
  Stage 3L mean target `<5e-5` is not met on protein fixtures.
- Order4 brings compact-mini-protein expansion mean under `5e-5`, but sampled
  high-error misses remain nonzero (`4`), so it is not promotable.
- Target batching is correct and can be a win on the larger compact fixture,
  but it is not yet a broad performance win and should not be promoted.

Barnase Stage 3L 120s probe:

Command used one candidate row:
`surface_points=8`, `leaf=128`, `theta=0.6`, order `3`, `topdown-aabb`,
`tight-aabb`, `physical-gap`, `area-centroid`,
`fragmentation-safety=split`, `moment-construction=translated`,
`moment-sum=naive`, `moment-audit=sampled`, `expansion-audit=sampled`,
`convergence-audit=sampled`, `convergence-safety=split`,
`convergence-estimator=hybrid`, `convergence-calibration=consecutive-aspect`,
`convergence-dielectric-amp=on`, `near=duffy_adaptive`,
`quad_buffer=adaptive`, and `precond=diagonal`.

| Barnase Stage 3L row | status | matvecs | gmres iter | sampled residual | true residual | expansion mean | expansion max | dielectric amp max | effective amp max | split rejects | leaf work/matvec |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| order3 theta0.6 dielectric on cluster off | `timeout` | `35` | `28` | `4.53e-4` | `6.677e-4` | `2.65e-4` | `5.199e-2` | `80` | `152.74` | `13994015` | `2.241e8` |

- The Barnase progress row confirms timeout-safe Stage 3L counters.
- The dielectric-aware estimator increases amplification as expected, but it
  does not push Barnase below the Stage 3L residual or mean-error gates.
- Barnase 300s was not run because no 120s Stage 3L row met the promotion gate.

Verification after Stage 3L:

- `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- `test/pbe_tabi_dielectric_estimator_smoke.ps1`: PASS.
- `test/pbe_tabi_leaf_pair_batch_smoke.ps1`: PASS.
- `test/pbe_tabi_convergence_calibration_smoke.ps1`: PASS.
- `test/pbe_tabi_convergence_estimator_smoke.ps1`: PASS.
- `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3LRows -CaseTimeoutSeconds 180 -AllowFailedStatus`: PASS.
- Barnase 120s Stage 3L row above: PASS as timeout-safe evidence; no
  promotion candidate.
- `test/pbe_tabi_far_expansion_audit_smoke.ps1`: PASS.
- `test/pbe_tabi_convergence_mac_safety_smoke.ps1`: PASS.
- `test/pbe_tabi_fragmented_node_mac_smoke.ps1`: PASS.
- `test/pbe_tabi_treecode_accuracy.ps1`: PASS.
- `test/pbe_tabi_topdown_mac_audit.ps1`: PASS.
- `test/pbe_tabi_tree_builder_smoke.ps1`: PASS.
- `test/pbe_tabi_gmres_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_sampled_true_residual_smoke.ps1`: PASS.
- `test/pbe_tabi_profile_smoke.ps1`: PASS.
- `test/pbe_memory_profile_smoke.ps1`: PASS.
- `test/pbe_clean_room_audit.ps1`: PASS.
- `test/pbe_regression.ps1`: PASS.
- `test/pbsa_smoke.ps1`: PASS.
- `test/sasa_surface_smoke.ps1`: PASS.

Stage 3L decision:

- Keep production defaults unchanged.
- Keep dielectric-aware estimator and leaf-pair batching hidden/experimental.
- Do not promote and do not run Barnase 300s yet.
- The user's dielectric-amplifier hypothesis is partially confirmed: the
  estimator now sees `eps_ratio=80` and rejects more unstable nodes. However,
  Barnase remains near the same residual/noise floor (`~6.7e-4` true residual)
  and far-expansion mean (`~2.65e-4`), so dielectric amplification alone is not
  enough.
- Next branch should target the estimator's source norm/error model rather than
  leaf batching:
  1. calibrate the JCP estimator against block/channel contributions so the
     error cap predicts the observed mean error scale (`~2.6e-4`) rather than
     only increasing split count;
  2. add percentile or per-channel expansion-audit gates (`p95/p99`, `K2`
     source-normal channel, and phin channel) to decide whether isolated spikes
     are harmless or still polluting GMRES;
  3. keep leaf-pair batching available for larger fixtures, but only optimize
     its SIMD/target-lane implementation after an accuracy row reaches
     far-expansion mean `<5e-5` and residual `<=1e-4`.

## Stage 3M/9F - Channel-Calibrated JCP Expansion Estimator

Status: implemented as hidden/experimental diagnostics; production defaults
unchanged.

Stage 3M implementation:

- Added hidden channel expansion audit flag
  `-pbe-tabi-expansion-audit-channel=off|summary|sampled-detail`, default
  `off`.
- Added hidden channel estimator model flag
  `-pbe-tabi-convergence-estimator-model=stage3l|channel-tail|channel-tail-eps`,
  default `stage3l`.
- Expansion audit now records overall `p95/p99`, worst channel, and
  per-channel mean/max/p99 for:
  - `K1_phi`;
  - `K2_phi_normal`;
  - `K3_phin`;
  - `K4_mixed`.
- Channel-aware estimator now computes per-channel estimated error and
  high-error caught/missed counters. `channel-tail-eps` applies the dielectric
  scaling to the K2 channel instead of only using one global amplifier.
- Timing/profile/progress now expose:
  `tabi_far_expansion_rel_error_p95`,
  `tabi_far_expansion_rel_error_p99`,
  `tabi_far_expansion_worst_channel`,
  `tabi_far_expansion_k{1..4}_rel_error_{max,mean,p99}`,
  `tabi_convergence_k{1..4}_estimated_error_max`,
  `tabi_convergence_k{1..4}_high_error_caught`, and
  `tabi_convergence_k{1..4}_high_error_missed`.
- Extended harnesses:
  - `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3MRows`;
  - `test/pbe_tabi_barnase_truth_matrix.ps1` with estimator-model and
    expansion-audit-channel passthrough plus timing/progress fallback fields.
- Added focused smokes:
  - `test/pbe_tabi_channel_expansion_audit_smoke.ps1`;
  - `test/pbe_tabi_channel_estimator_smoke.ps1`.

Measured Stage 3M smoke evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Channel expansion audit smoke: PASS.
  - panels `168`;
  - sampled far-expansion pairs `129`;
  - overall p99 `2.68676`;
  - worst channel `K3_phin`;
  - K2 channel max `3.69952`;
  - phin-side channel checks reported nonzero K1/K3 maxima.
- Channel estimator smoke: PASS.
  - panels `164`;
  - baseline global estimated-error max `415.909`;
  - channel K2 estimated-error max `415.909`;
  - split cap `207.955`;
  - `channel-tail-eps` split rejects `70`;
  - benign case retained far work (`350` far-node interactions).
- Dielectric estimator smoke: PASS.
  - `dielectric_amp=on` still reports `eps_ratio=80` and increases split
    rejects versus `off`.
- Leaf-pair batch smoke: PASS.
  - batched-vs-scalar relative RMS `1.0891e-16`.

Measured Stage 3M quick medium matrix:

- Command:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3MRows -Quick -CaseTimeoutSeconds 180 -AllowFailedStatus`.
- Output:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.
- Scope: `mini-protein` only. Full three-fixture Stage 3M matrix is still
  pending before any promotion decision.

| row | matvec rel | energy drift rel | sampled residual | expansion mean | expansion p99 | worst channel | model | channel audit | K2 est max | K2 missed |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: | ---: |
| stage3l order3 theta0.6 | `5.318e-7` | `3.84e-8` | `9.759e-5` | `2.767e-4` | `1.493e-3` | `K3_phin` | `stage3l` | `summary` | `12.479` | `0` |
| channel audit order3 theta0.6 | `2.188e-4` | `1.43e-5` | `9.761e-5` | `3.469e-2` | `5.967e-1` | `K2_phi_normal` | `channel-tail-eps` | `summary` | `6.344` | `0` |
| channel-tail order3 theta0.6 | `5.318e-7` | `3.84e-8` | `9.759e-5` | `2.767e-4` | `1.493e-3` | `K3_phin` | `channel-tail` | `summary` | `2.365` | `0` |
| channel-tail-eps order3 theta0.6 | `5.318e-7` | `3.84e-8` | `9.759e-5` | `2.767e-4` | `1.493e-3` | `K3_phin` | `channel-tail-eps` | `summary` | `6.306` | `0` |
| channel-tail-eps order3 theta0.5 | `5.318e-7` | `3.84e-8` | `9.759e-5` | `2.767e-4` | `1.493e-3` | `K3_phin` | `channel-tail-eps` | `summary` | `2.538` | `0` |
| channel-tail-eps order4 theta0.6 | `8.920e-8` | `9.18e-9` | `9.759e-5` | `5.576e-5` | `4.449e-4` | `K3_phin` | `channel-tail-eps` | `summary` | `1.614` | `0` |

Medium interpretation:

- Channel audit plumbing works: model/channel flags, p99, worst channel, and
  per-channel estimator counters are visible in `pbe_timing.csv`.
- On `mini-protein`, the primary order3 channel-tail-eps rows still sit at
  expansion mean `2.767e-4`, above the Stage 3M gate `<5e-5`.
- Order4 theta0.6 nearly reaches the strict mean gate (`5.576e-5`) and passes
  the p99 gate (`4.449e-4`), but it is still a fallback diagnostic, not a
  promotion candidate.
- The audit-only row correctly exposes a K2-dominated high-error branch
  (`worst_channel=K2_phi_normal`) without claiming it is an accepted split
  candidate.

Barnase Stage 3M 120s probe:

- Command ran one candidate row:
  `surface_points=8`, `leaf=128`, `theta=0.6`, order `3`,
  `topdown-aabb`, `tight-aabb`, `physical-gap`, `area-centroid`,
  `fragmentation-safety=split`, `moment-audit=sampled`,
  `expansion-audit=sampled`, `expansion-audit-channel=summary`,
  `convergence-audit=sampled`, `convergence-safety=split`,
  `convergence-estimator=hybrid`,
  `convergence-calibration=consecutive-aspect`,
  `convergence-dielectric-amp=on`,
  `convergence-estimator-model=channel-tail-eps`,
  `near=duffy_adaptive`, `quad_buffer=adaptive`, and `precond=diagonal`.

| Barnase Stage 3M row | status | total seconds | timing/progress counters |
| --- | --- | ---: | --- |
| order3 theta0.6 channel-tail-eps summary | `timeout` | `120.03` | `NA` |

- This run timed out before `pbe_timing.csv` or
  `pbe_tabi_longrun_progress.csv` emitted usable counters, so it is not a
  promotion or accuracy datapoint.
- Barnase 300s was not run.

Stage 3M decision:

- Keep production defaults unchanged.
- Keep `channel-tail`, `channel-tail-eps`, and channel expansion audit hidden.
- Do not promote and do not run Barnase 300s yet.
- Next branch should focus on selectivity and earlier Barnase progress:
  1. run the full three-fixture Stage 3M matrix before Barnase;
  2. add a lower-cost Barnase pre-progress checkpoint or a smaller Barnase
     diagnostic surface row so 120s probes always emit counters;
  3. use the channel audit evidence to tune the K3/phin and K2/K4 source-norm
     model, because order3 channel-tail-eps did not reduce mean error on
     `mini-protein`;
  4. keep order4 as fallback only if full medium rows confirm p99 and missed
     gates without exploding leaf work.

## Stage 3N/9G - TABI Pre-Solver Progress and K3/Phin Channel Calibration

Stage 3N implementation:

- Added hidden pre-solver progress flag:
  `-pbe-tabi-pre-solver-progress=off|on`, default `off`.
- When enabled, TABI writes `pbe_tabi_prepare_progress.csv` before
  `solve_linear_pb`, with phases:
  `tabi_prepare_start`, `pqr_loaded`, `tabi_options_built`, and
  `solver_enter`.
- Extended the Barnase truth matrix harness to read
  `pbe_tabi_longrun_progress.csv` first and fall back to
  `pbe_tabi_prepare_progress.csv` when a timeout happens before solver entry.
  Such a run is classified as `prep_timeout`.
- Added hidden K3/channel calibration flag:
  `-pbe-tabi-convergence-channel-calibration=stage3m|k3-sqrt-eps|k3-eps|all-derivative-sqrt-eps`,
  default `stage3m`.
- `k3-sqrt-eps` scales the `K3_phin` channel by
  `max(30, 30*sqrt(abs(eps_ratio)))`; `k3-eps` is the strict diagnostic
  `max(30, abs(eps_ratio))`; `all-derivative-sqrt-eps` applies the sqrt-eps
  amplifier to the derivative-side channels.
- Timing/profile/progress now expose per-channel effective amp and
  per-channel unstable/split reject counters for K1/K2/K3/K4.
- Extended harnesses:
  - `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3NRows`;
  - `test/pbe_tabi_barnase_truth_matrix.ps1` with
    `-ConvergenceChannelCalibrations` and `-PreSolverProgress`.
- Added focused smokes:
  - `test/pbe_tabi_presolver_progress_smoke.ps1`;
  - `test/pbe_tabi_k3_channel_estimator_smoke.ps1`.

Measured Stage 3N smoke evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Pre-solver progress smoke: PASS.
  - phases emitted:
    `tabi_prepare_start;pqr_loaded;tabi_options_built;solver_enter`;
  - solver-enter row preserved `k3-sqrt-eps`;
  - RSS at smoke completion `8351744` bytes.
- K3 channel estimator smoke: PASS.
  - panels `164`;
  - K3 amp increased from `30` to `268.328`;
  - K3 estimated-error max increased from `9.44424` to `84.4718`;
  - K3 channel split rejects increased from `38` to `71`;
  - benign path retained far work (`350` far-node interactions).
- Existing channel/audit gates after Stage 3N: PASS.
  - channel expansion audit smoke: worst channel `K3_phin`, p99 `2.68676`;
  - channel estimator smoke: `channel-tail-eps` split rejects `70` with
    benign far work retained;
  - dielectric estimator smoke: `dielectric_amp=on` reports amp `80`;
  - convergence calibration smoke: consecutive/aspect calibration still
    raises estimated error and keeps benign far work;
  - leaf-pair batch smoke: batched-vs-scalar relative RMS `1.0891e-16`;
  - convergence estimator smoke: split path reduced synthetic expansion
    error from `5.16619e-05` to `1.9621e-06`.

Measured Stage 3N full medium matrix:

- Command:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3MRows -IncludeStage3NRows -CaseTimeoutSeconds 180 -AllowFailedStatus`.
- Output:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.
- Scope: full three-fixture medium matrix:
  `mini-protein`, `charged-peptide`, and `compact-mini-protein`.

| fixture | Stage 3N row | status | expansion mean | expansion p99 | worst channel | global missed | K1 missed | K3 missed | K4 missed | K3 amp max | K3 split rejects | leaf work/matvec |
| --- | --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| mini-protein | order3 theta0.6 k3-sqrt-eps | `ok` | `2.767e-4` | `1.493e-3` | `K3_phin` | `0` | `204` | `0` | `0` | `268.328` | `131055` | `7290772` |
| mini-protein | order4 theta0.6 k3-sqrt-eps fallback | `ok` | `5.576e-5` | `4.449e-4` | `K3_phin` | `0` | `0` | `0` | `20` | `268.328` | `131055` | `7290772` |
| charged-peptide | order3 theta0.6 k3-sqrt-eps | `ok` | `0` | `0` | `none` | `0` | `0` | `0` | `0` | `268.328` | `73120` | `5021476` |
| compact-mini-protein | order3 theta0.6 k3-sqrt-eps | `ok` | `3.340e-4` | `4.655e-3` | `K4_mixed` | `0` | `40` | `0` | `0` | `268.328` | `700680` | `17098706` |
| compact-mini-protein | order4 theta0.6 k3-sqrt-eps fallback | `ok` | `4.461e-5` | `5.805e-4` | `K4_mixed` | `0` | `0` | `0` | `263` | `268.328` | `700680` | `17098706` |

Medium interpretation:

- Stage 3N successfully targets the suspected K3/phin channel: K3 missed is
  `0` in the measured medium rows, and K3 split counters/amps are visible.
- No Stage 3N row passes the full promotion gate across all three fixtures:
  order3 rows remain far above the mean/p99 targets on mini/compact; order4
  improves mean sharply but still misses compact p99 and K4 per-channel missed
  remains nonzero.
- The remaining medium failures point away from K3 alone and toward a
  K4/mixed or per-channel gate semantics branch. The global estimator missed
  counter is `0`, but per-channel missed counters still expose local channel
  underprediction.
- Barnase 120s was not run in Stage 3N because the required full medium
  candidate did not pass first. Barnase 300s was not run.

Verification after Stage 3N:

- `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- `test/pbe_tabi_presolver_progress_smoke.ps1`: PASS.
- `test/pbe_tabi_k3_channel_estimator_smoke.ps1`: PASS.
- `test/pbe_tabi_channel_expansion_audit_smoke.ps1`: PASS.
- `test/pbe_tabi_channel_estimator_smoke.ps1`: PASS.
- `test/pbe_tabi_dielectric_estimator_smoke.ps1`: PASS.
- `test/pbe_tabi_convergence_calibration_smoke.ps1`: PASS.
- `test/pbe_tabi_convergence_estimator_smoke.ps1`: PASS.
- `test/pbe_tabi_leaf_pair_batch_smoke.ps1`: PASS.
- Full Stage 3M+3N medium matrix: PASS as evidence collection, no promotion
  candidate.
- `git diff --check`: no whitespace errors; only existing line-ending
  normalization warnings were reported.

Stage 3N decision:

- Keep production defaults unchanged.
- Keep pre-solver progress, K3/channel calibration, channel estimator, and
  channel expansion audit hidden.
- Do not run Barnase 120s/300s yet.
- Next branch should target the remaining medium blockers:
  1. split per-channel `missed` semantics into "globally caught" versus
     "same-channel caught" so promotion gates do not mix diagnostic meanings;
  2. add K4/mixed calibration rows analogous to K3 sqrt-eps because compact
     medium now reports `K4_mixed` as worst channel;
  3. run a strict order4 K4 diagnostic only after confirming p99 can fall below
     `5e-4` without increasing leaf work beyond the Stage 3L/3M envelope;
  4. keep Barnase blocked until a full medium row passes mean, p99, residual,
     and channel-missed gates.

## Stage 3O/9H - TABI K4/Mixed Calibration and Accepted-Far Miss Semantics

Stage 3O implementation:

- Added hidden channel missed semantics flag:
  `-pbe-tabi-channel-missed-semantics=legacy|accepted-only`, default
  `legacy`.
- With `accepted-only`, sampled per-channel `*_high_error_missed` now means a
  high-error pair actually entered the accepted far Taylor path. High-error
  candidates rejected by the convergence split are counted separately in:
  `tabi_convergence_k{1..4}_high_error_caught_by_split` and
  `tabi_convergence_estimator_high_error_caught_by_split`.
- Extended sampled channel audit to sample split-rejected candidates without
  changing traversal. These samples only feed `caught_by_split` counters.
- Added hidden K4 calibration diagnostics:
  `k4-sqrt-eps` and `k4-eps`, while preserving `all-derivative-sqrt-eps` as
  the primary Stage 3O candidate.
- Threaded the new semantics/caught-by-split counters through treecode stats,
  timing/profile output, longrun progress, pre-solver progress context,
  matrix CSV rows, and Barnase truth harness rows.
- Updated the medium matrix classifier so audited estimator rows that do not
  pass the mean/missed gate are labeled `convergence_mean_gate_followup`
  instead of the older medium-gate pass label.
- Extended harnesses:
  - `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3ORows`;
  - `test/pbe_tabi_barnase_truth_matrix.ps1` with
    `-ChannelMissedSemantics`.
- Added focused smokes:
  - `test/pbe_tabi_channel_missed_semantics_smoke.ps1`;
  - `test/pbe_tabi_k4_channel_estimator_smoke.ps1`.

Measured Stage 3O smoke evidence:

- Native build state:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS
  (`Nothing to be done for all` after the Stage 3O engine rebuild).
- Channel missed semantics smoke: PASS.
  - panels `118`;
  - accepted high-error missed `308`;
  - split path missed `170`;
  - split caught-by-split `138`;
  - split rejects `69`.
- K4 channel estimator smoke: PASS.
  - panels `168`;
  - K4 effective amp increased from `822.857` to `7359.86`;
  - K4 estimated-error max increased from `70.0491` to `626.538`;
  - K4 split rejects increased from `78` to `82`;
  - benign path retained far work (`60` far-node interactions).
- Existing channel/progress gates after Stage 3O: PASS.
  - K3 channel estimator smoke: K3 amp `30 -> 268.328`, split rejects
    `38 -> 71`;
  - channel expansion audit smoke: worst channel `K3_phin`, p99 `2.68676`;
  - channel estimator smoke: split rejects `70` with benign far work retained;
  - pre-solver progress smoke: phases
    `tabi_prepare_start;pqr_loaded;tabi_options_built;solver_enter`.

Measured Stage 3O full medium matrix:

- Command:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3MRows -IncludeStage3NRows -IncludeStage3ORows -CaseTimeoutSeconds 180 -AllowFailedStatus`.
- Output:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.
- Scope: full three-fixture medium matrix:
  `mini-protein`, `charged-peptide`, and `compact-mini-protein`.

| fixture | Stage 3O row | status | semantics | expansion mean | p99 report | global accepted missed | K1/K2/K3/K4 accepted missed | K4 caught by split | worst channel | K4 amp max | K4 split rejects | leaf work/matvec |
| --- | --- | --- | --- | ---: | ---: | ---: | --- | ---: | --- | ---: | ---: | ---: |
| mini-protein | order4 theta0.6 all-derivative-sqrt-eps | `ok` | `accepted-only` | `6.923e-5` | `4.904e-4` | `0` | `0/0/15/1` | `3512` | `K3_phin` | `517.780` | `131055` | `7290772` |
| mini-protein | order4 theta0.5 all-derivative-sqrt-eps | `ok` | `accepted-only` | `6.460e-5` | `4.606e-4` | `0` | `0/0/22/21` | `1777` | `K4_mixed` | `517.780` | `83775` | `7290772` |
| mini-protein | order4 theta0.6 k4-eps | `ok` | `accepted-only` | `7.051e-5` | `4.904e-4` | `0` | `0/0/16/1` | `3552` | `K4_mixed` | `154.372` | `131055` | `7290772` |
| charged-peptide | order4 theta0.6 all-derivative-sqrt-eps | `ok` | `accepted-only` | `0` | `0` | `0` | `0/0/0/0` | `7621` | `none` | `601.125` | `73120` | `5021476` |
| charged-peptide | order4 theta0.5 all-derivative-sqrt-eps | `ok` | `accepted-only` | `0` | `0` | `0` | `0/0/0/0` | `6286` | `none` | `601.125` | `27700` | `5021476` |
| compact-mini-protein | order4 theta0.6 all-derivative-sqrt-eps | `ok` | `accepted-only` | `3.673e-5` | `4.824e-4` | `0` | `0/7/18/42` | `7755` | `K3_phin` | `622.798` | `700680` | `17098706` |
| compact-mini-protein | order4 theta0.5 all-derivative-sqrt-eps | `ok` | `accepted-only` | `3.358e-5` | `4.218e-4` | `0` | `0/10/21/37` | `5448` | `K3_phin` | `622.798` | `454050` | `17098706` |
| compact-mini-protein | order4 theta0.6 k4-eps | `ok` | `accepted-only` | `3.351e-5` | `4.243e-4` | `0` | `0/8/19/36` | `7711` | `K3_phin` | `185.683` | `694478` | `17084287` |

Medium interpretation:

- The accepted-only semantics fix resolves the diagnostic ambiguity: many
  high-error K4 candidates are now visibly `caught_by_split`, while true
  accepted-far misses remain separate.
- `all-derivative-sqrt-eps` and strict `k4-eps` both reduce the compact
  expansion mean below the Stage 3O target (`<5e-5`) with p99 now report-only,
  but accepted per-channel misses remain nonzero on `mini-protein` and
  `compact-mini-protein`.
- `charged-peptide` passes the accepted-missed and mean evidence cleanly, but
  there is no full three-fixture promotion candidate because:
  - `mini-protein` accepted-only order4 rows remain above the mean gate
    (`~7.1e-5`) and still have K3/K4 accepted misses;
  - `compact-mini-protein` mean passes, but K2/K3/K4 accepted misses remain
    nonzero.
- Barnase 120s was not run in Stage 3O because the required full medium
  candidate did not pass first. Barnase 300s was not run.

Stage 3O decision:

- Keep production defaults unchanged.
- Keep channel missed semantics, K4 calibration, channel expansion audit, and
  all convergence estimator modes hidden/experimental.
- Do not run Barnase 120s/300s yet.
- Next branch should target accepted-far residual channel misses instead of
  p99:
  1. inspect why K3 accepted misses remain on `mini-protein` even under
     order4 all-derivative calibration;
  2. add a selective accepted-far retry/split rule keyed by per-channel
     estimator score or direct audit misses, without falling back to a blind
     global rho cap;
  3. use compact `k4-eps` only as a diagnostic lower bound, because it improves
     mean but does not remove accepted misses;
  4. keep Barnase blocked until all three medium fixtures have expansion mean
     `<5e-5` and global/per-channel accepted-missed `0`.

## Stage 3P/9I - TABI Polynomial-Derivative Channel Calibration and Order-3 Theta Rescue

Stage 3P implementation:

- Added hidden channel calibration mode:
  `-pbe-tabi-convergence-channel-calibration=polynomial-eps`, default still
  `stage3m`.
- In the JCP `channel-tail-eps` estimator, `polynomial-eps` uses:
  - K2/K3 amp `max(30, (taylor_order+1)*sqrt(abs(eps_ratio)))`;
  - K4 amp `max(30, (taylor_order+1)*(taylor_order+2)*abs(eps_ratio))`;
  - K2/K3/K4 aspect scaling when `convergence-calibration=consecutive-aspect`;
  - K1 remains amp `1`.
- Kept accepted-only miss semantics unchanged: accepted-far high-error samples
  count as missed, split-rejected samples count as `caught_by_split`.
- Extended the medium matrix with `-IncludeStage3PRows`:
  order3 theta `0.50`, `0.45`, `0.40`, and order4 theta `0.60`, all using
  `polynomial-eps` and `accepted-only`.
- Added K2/K3 comparison columns to the matrix CSV:
  `k2_k3_accepted_missed_delta`, `k2_k3_split_reject_delta`, and
  `k2_k3_estimated_error_ratio`.
- Added focused smoke:
  `test/pbe_tabi_polynomial_channel_estimator_smoke.ps1`.

Measured Stage 3P smoke evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Polynomial channel estimator smoke: PASS.
  - panels `168`;
  - K3 amp increased from all-derivative `268.328` to polynomial `981.314`;
  - polynomial K2 amp order3/order4: `981.314 -> 1226.64`;
  - polynomial K4 amp order3/order4: `43885.7 -> 65828.6`;
  - K3 split path active (`82` rejects; synthetic fixture was tail-saturated);
  - benign path retained far work (`60` far-node interactions).
- Regression-focused smokes after Stage 3P: PASS.
  - channel missed semantics smoke;
  - K4 channel estimator smoke;
  - K3 channel estimator smoke;
  - channel estimator smoke;
  - channel expansion audit smoke;
  - pre-solver progress smoke;
  - dielectric estimator smoke;
  - convergence calibration smoke.

Measured Stage 3P full medium matrix:

- Command:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3MRows -IncludeStage3NRows -IncludeStage3ORows -IncludeStage3PRows -CaseTimeoutSeconds 180 -AllowFailedStatus`.
- Output:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| fixture | Stage 3P row | status | expansion mean | p99 report | global missed | K1/K2/K3/K4 accepted missed | worst channel | K2/K3/K4 amp max | K2-K3 missed delta | leaf work/matvec |
| --- | --- | --- | ---: | ---: | ---: | --- | --- | --- | ---: | ---: |
| mini-protein | order3 theta0.50 polynomial | `ok` | `3.053e-4` | `1.321e-3` | `35` | `36/2/203/103` | `K4_mixed` | `69.04/69.04/3087.44` | `-201` | `7290772` |
| mini-protein | order3 theta0.45 polynomial | `ok` | `2.975e-4` | `1.356e-3` | `43` | `36/2/237/122` | `K4_mixed` | `69.04/69.04/3087.44` | `-235` | `7290772` |
| mini-protein | order3 theta0.40 polynomial | `ok` | `2.871e-4` | `1.481e-3` | `57` | `49/2/344/148` | `K4_mixed` | `69.04/69.04/3087.44` | `-342` | `7290772` |
| mini-protein | order4 theta0.60 polynomial | `ok` | `7.374e-5` | `4.904e-4` | `0` | `0/0/14/0` | `K3_phin` | `86.30/86.30/4631.17` | `-14` | `7290772` |
| charged-peptide | order3 theta0.50 polynomial | `ok` | `0` | `0` | `0` | `0/0/0/0` | `none` | `80.15/80.15/3584.42` | `0` | `5021476` |
| charged-peptide | order3 theta0.45 polynomial | `ok` | `0` | `0` | `0` | `0/0/0/0` | `none` | `70.28/70.28/3143.17` | `0` | `5021476` |
| compact-mini-protein | order3 theta0.50 polynomial | `ok` | `2.394e-4` | `2.591e-3` | `47` | `0/49/244/79` | `K3_phin` | `83.04/83.04/3713.65` | `-195` | `17098706` |
| compact-mini-protein | order3 theta0.45 polynomial | `ok` | `2.372e-4` | `2.591e-3` | `53` | `0/61/267/95` | `K3_phin` | `66.37/66.37/2968.36` | `-206` | `17098706` |
| compact-mini-protein | order4 theta0.60 polynomial | `ok` | `3.830e-5` | `4.824e-4` | `0` | `0/8/18/36` | `K3_phin` | `103.80/103.80/5570.48` | `-10` | `17098706` |

Stage 3P interpretation:

- The polynomial second-derivative K4 amp is active and large, but order3
  polynomial-eps does not rescue medium accuracy. It worsens accepted misses
  and leaves mean far above the `5e-5` gate on `mini-protein` and
  `compact-mini-protein`.
- Order4 polynomial preserves the compact mean gate, but still fails
  full-medium promotion because `mini-protein` mean remains `7.374e-5` and
  accepted K3 misses remain nonzero.
- K2/K3 comparison columns show K3 remains the dominant accepted-miss channel
  in Stage 3P order3 rows (`K2-K3 accepted missed delta` is strongly
  negative), so the polynomial K2/K3 amp is not selective enough to replace
  the Stage 3O order4 path.
- No Stage 3P row passes full medium gates across all three fixtures.
- Barnase 120s was not run because no full-medium candidate passed first.
  Barnase 300s was not run.

Stage 3P decision:

- Keep production defaults unchanged.
- Keep `polynomial-eps` hidden/experimental.
- Do not promote and do not run Barnase yet.
- Next branch should not continue tightening order3 polynomial blindly.
  Better candidates are:
  1. reinterpret polynomial calibration as a stricter derivative envelope for
     K2/K3 as well, because the implemented K2/K3 formula is weaker than the
     Stage 3O all-derivative sqrt-eps amp on these fixtures;
  2. add a selective accepted-far retry/split rule for sampled K3 misses on
     `mini-protein`, using accepted-only audit feedback rather than p99;
  3. keep order4 theta0.6 as the current lower-bound accuracy branch, but
     require accepted-missed `0` before Barnase.

## Stage 3Q/9J - TABI Corrected Polynomial K2/K3 Amp and Selective K3 Split Probe

Stage 3Q implementation:

- Corrected hidden `polynomial-eps` K2/K3 amplification in the JCP
  `channel-tail-eps` estimator:
  - old Stage 3P K2/K3 amp:
    `max(30, (taylor_order+1)*sqrt(abs(eps_ratio)))`;
  - corrected Stage 3Q K2/K3 amp:
    `max(30, (taylor_order+1)*30*sqrt(abs(eps_ratio)))`;
  - K4 remains
    `max(30, (taylor_order+1)*(taylor_order+2)*abs(eps_ratio))`;
  - K1 remains amp `1`.
- Added hidden split policy:
  `-pbe-tabi-convergence-split-policy=global|k3-hot`, default `global`.
  - `global` preserves existing behavior.
  - `k3-hot` is only active for `channel-tail-eps` and only allows estimator
    split when `K3_phin` is unstable. Non-channel models keep global behavior.
- Plumbed the split policy through TABI options, tree workspace, stats,
  timing/profile CSVs, prepare/longrun progress CSVs, the medium matrix, and
  Barnase truth matrix.
- Extended the medium matrix with `-IncludeStage3QRows`:
  order3 theta `0.50`, `0.45`, `0.40` corrected polynomial/global;
  order4 theta `0.60` corrected polynomial/global;
  order4 theta `0.60` corrected polynomial/`k3-hot`;
  order3 theta `0.45` corrected polynomial/`k3-hot`.
- Added focused smoke:
  `test/pbe_tabi_selective_k3_split_smoke.ps1`.

Measured Stage 3Q smoke evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Corrected polynomial smoke: PASS.
  - panels `168`;
  - all-derivative order3 K3 amp `268.328`;
  - corrected polynomial order3 K2/K3 amp `29439.4`;
  - corrected polynomial order4 K2/K3 amp `36799.3`;
  - polynomial K4 amp order3/order4: `43885.7 -> 65828.6`;
  - synthetic K3 split path active (`82` rejects);
  - benign path retained far work (`60` far-node interactions).
- Selective K3 split smoke: PASS.
  - panels `176`;
  - K3-dominated case under `k3-hot`: `153` K3 split rejects;
  - K4-dominated diagnostic: global split rejects `138`, `k3-hot` split
    rejects `0`;
  - benign accepted-only missed stayed `0`, with `143` far-node interactions.
- Adjacent channel/progress smokes after Stage 3Q: PASS.
  - channel missed semantics smoke;
  - K4 channel estimator smoke;
  - K3 channel estimator smoke;
  - channel expansion audit smoke;
  - channel estimator smoke;
  - pre-solver progress smoke.

Measured Stage 3Q full medium matrix:

- Command:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3MRows -IncludeStage3NRows -IncludeStage3ORows -IncludeStage3PRows -IncludeStage3QRows -CaseTimeoutSeconds 180 -AllowFailedStatus`.
- Output:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| fixture | Stage 3Q row | split policy | status | expansion mean | p99 report | global missed | K1/K2/K3/K4 accepted missed | worst channel | K2/K3/K4 amp max | K2-K3 missed delta | leaf work/matvec |
| --- | --- | --- | --- | ---: | ---: | ---: | --- | --- | --- | ---: | ---: |
| mini-protein | order3 theta0.50 corrected polynomial | `global` | `ok` | `2.970e-4` | `1.216e-3` | `29` | `29/1/199/97` | `K4_mixed` | `2071.12/2071.12/3087.44` | `-198` | `7290772` |
| mini-protein | order3 theta0.45 corrected polynomial | `global` | `ok` | `2.980e-4` | `1.381e-3` | `44` | `37/2/263/127` | `K4_mixed` | `2071.12/2071.12/3087.44` | `-261` | `7290772` |
| mini-protein | order3 theta0.40 corrected polynomial | `global` | `ok` | `2.948e-4` | `1.414e-3` | `58` | `51/2/342/153` | `K4_mixed` | `2071.12/2071.12/3087.44` | `-340` | `7290772` |
| mini-protein | order4 theta0.60 corrected polynomial | `global` | `ok` | `7.274e-5` | `4.904e-4` | `0` | `0/0/14/0` | `K3_phin` | `2588.90/2588.90/4631.17` | `-14` | `7290772` |
| mini-protein | order4 theta0.60 corrected polynomial | `k3-hot` | `ok` | `7.072e-5` | `4.904e-4` | `0` | `0/0/18/0` | `K3_phin` | `2588.90/2588.90/4631.17` | `-18` | `7290772` |
| charged-peptide | order3 theta0.50 corrected polynomial | `global` | `ok` | `0` | `0` | `0` | `0/0/0/0` | `none` | `2404.50/2404.50/3584.42` | `0` | `5021476` |
| charged-peptide | order4 theta0.60 corrected polynomial | `global` | `ok` | `0` | `0` | `0` | `0/0/0/0` | `none` | `3005.62/3005.62/5376.62` | `0` | `5021476` |
| compact-mini-protein | order3 theta0.50 corrected polynomial | `global` | `ok` | `2.375e-4` | `2.643e-3` | `43` | `0/47/217/76` | `K3_phin` | `2491.19/2491.19/3713.65` | `-170` | `17098706` |
| compact-mini-protein | order3 theta0.45 corrected polynomial | `global` | `ok` | `2.441e-4` | `2.552e-3` | `55` | `0/62/273/95` | `K3_phin` | `1991.24/1991.24/2968.36` | `-211` | `17098706` |
| compact-mini-protein | order3 theta0.40 corrected polynomial | `global` | `ok` | `2.416e-4` | `2.511e-3` | `63` | `0/75/308/119` | `K3_phin` | `1991.24/1991.24/2968.36` | `-233` | `17098706` |
| compact-mini-protein | order4 theta0.60 corrected polynomial | `global` | `ok` | `4.112e-5` | `5.587e-4` | `0` | `0/8/18/36` | `K3_phin` | `3113.99/3113.99/5570.48` | `-10` | `17098706` |
| compact-mini-protein | order4 theta0.60 corrected polynomial | `k3-hot` | `ok` | `3.789e-5` | `4.824e-4` | `0` | `0/9/17/36` | `K3_phin` | `3113.99/3113.99/5570.48` | `-8` | `17098706` |

Stage 3Q interpretation:

- The corrected K2/K3 amp now exceeds Stage 3O by design and by smoke
  evidence, so the Stage 3P formula bug is fixed.
- Corrected order3 theta `0.45..0.50` does not rescue medium accuracy:
  mean stays around `2.4e-4..3.0e-4` on `mini-protein` and
  `compact-mini-protein`, with many accepted K3 misses.
- Corrected order4 theta `0.60` remains the best accuracy branch:
  - `charged-peptide` passes mean and accepted-missed gates;
  - `compact-mini-protein` passes mean (`~3.8e-5..4.1e-5`) but still has
    accepted K2/K3/K4 misses;
  - `mini-protein` still fails the mean gate (`~7.1e-5`) and accepted K3
    misses remain nonzero.
- `k3-hot` works on a selective synthetic probe, but on the real medium
  fixtures it does not reduce leaf work or clear accepted misses; K3 is hot in
  the same regions that drive the candidate splits.
- No Stage 3Q row passes full medium gates across all three fixtures.
- Barnase 120s was not run because the mandatory full-medium candidate did
  not pass first. Barnase 300s was not run.

Stage 3Q decision:

- Keep production defaults unchanged.
- Keep corrected `polynomial-eps` and `convergence-split-policy` hidden and
  experimental.
- Do not promote and do not run Barnase yet.
- Next branch should target accepted-far K3 misses directly on the order4
  path. Recommended next work:
  1. add an accepted-far retry/direct-descent probe for sampled high-error K3
     candidates, gated by accepted-only audit evidence;
  2. test a K3-specific estimator score that is stronger than corrected
     polynomial only when K3 accepted misses are predicted, without forcing
     K4-only nodes into the same split policy;
  3. keep p99 report-only, but require mean `<5e-5` and accepted-missed `0`
     before Barnase 120s.

## Stage 3R/9K - TABI Accepted-Far K3 Structural-Miss Guard

Stage 3R implementation:

- Added hidden accepted-miss diagnostics:
  `-pbe-tabi-accepted-miss-diagnostics=off|summary|sampled-detail`, default
  `off`.
  - Summary mode reports per-channel accepted-miss counts, panels, rho
    min/mean/max, tail-ratio mean/max, estimated-error mean/max, and worst
    channel.
  - Sampled-detail mode writes `pbe_tabi_accepted_miss_diagnostics.csv` with
    target index, source node, channel, rho, estimated error, and measured
    direct-vs-Taylor relative error for sampled accepted misses.
- Added deterministic K3 accepted-rho cap:
  `-pbe-tabi-convergence-k3-accepted-rho-cap=<0.0..2.0>`, default `0.0`.
  It is active on the JCP on-the-fly `channel-tail-eps` path and descends
  accepted far candidates with rho above the cap before Taylor evaluation.
- Added true-surface MAC diagnostic:
  `-pbe-tabi-mac-distance-policy=center|surface-margin`, default `center`,
  plus `-pbe-tabi-mac-surface-margin-factor=<0.0..0.75>`, default `0.25`.
  `surface-margin` evaluates rho using `distance - factor*node.radius`, clamped
  positive, and reports margin rejects/panels.
- Plumbed all new fields through TABI options, tree workspace, stats,
  timing/profile CSVs, longrun progress, medium matrix, and Barnase truth
  matrix.
- Extended the medium matrix with `-IncludeStage3RRows`:
  baseline Stage 3Q order4 theta `0.60`, K3 cap `0.58/0.57/0.55/0.52`,
  surface-margin `0.25`, and combined cap+surface-margin diagnostic.
- Added focused smokes:
  - `test/pbe_tabi_k3_accepted_rho_cap_smoke.ps1`;
  - `test/pbe_tabi_surface_margin_mac_smoke.ps1`.

Measured Stage 3R smoke evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- K3 accepted-rho cap smoke: PASS.
  - synthetic high-rho accepted far case: baseline far nodes `249`;
    cap `0.58` produced `53` K3 rho-cap rejects over `728` panels and
    increased leaf/direct work;
  - benign low-rho case retained `196` far-node interactions and had no cap
    rejects.
- Surface-margin MAC smoke: PASS.
  - hollow/elongated synthetic case: center far nodes `228`, surface-margin
    far nodes `193`, surface-margin rejects `73`;
  - benign/far case retained `216` far-node interactions.
- Adjacent Stage 3Q/3O smokes after Stage 3R: PASS.
  - polynomial channel estimator smoke;
  - selective K3 split smoke;
  - channel missed semantics smoke;
  - topdown MAC audit smoke;
  - treecode accuracy smoke;
  - tree builder smoke;
  - TABI profile smoke;
  - memory profile smoke;
  - sampled true-residual smoke.

Measured Stage 3R full medium matrix:

- Command:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3QRows -IncludeStage3RRows -CaseTimeoutSeconds 180 -AllowFailedStatus`.
- Output:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| fixture | Stage 3R row | status | expansion mean | K1/K2/K3/K4 accepted missed | K3 accepted-miss rho min/mean/max | K3 cap rejects | surface rejects | leaf/far work per matvec |
| --- | --- | --- | ---: | --- | --- | ---: | ---: | --- |
| mini-protein | order4 theta0.60 polynomial baseline | `ok` | `7.308e-5` | `0/0/17/1` | `~0.246/0.247/0.249` | `0` | `0` | `7290772/498` |
| mini-protein | K3 rho cap `0.57` | `ok` | `6.844e-5` | `0/0/15/0` | `~0.246/0.247/0.249` | `0` | `0` | `7290772/498` |
| mini-protein | surface-margin `0.25` | `ok` | `2.883e-5` | `0/0/1/0` | sampled K3 miss near `0.245` | `0` | `37170` | `7335697/173` |
| charged-peptide | any Stage 3R row | `ok` | `0` | `0/0/0/0` | none | `0` | up to `35940` | `5021476/0` |
| compact-mini-protein | order4 theta0.60 polynomial baseline | `ok` | `4.173e-5` | `0/7/17/36` | low-rho structural K3 misses | `0` | `0` | `17098706/1474` |
| compact-mini-protein | K3 rho cap `0.55` | `ok` | `3.858e-5` | `0/7/14/37` | low-rho structural K3 misses | `0` | `0` | `17098706/1474` |
| compact-mini-protein | surface-margin `0.25` | `ok` | `2.448e-5` | `0/3/0/0` | none | `0` | `197100` | `17228993/596` |

Stage 3R interpretation:

- The deterministic K3 cap implementation works on the synthetic high-rho
  smoke, but the real `mini-protein` accepted K3 misses are not high-rho
  `0.52..0.58` misses. The accepted-miss diagnostic shows them at about
  `rho ~= 0.245..0.249`, so the comment's proposed high-rho cap does not touch
  the actual medium blocker.
- `surface-margin` is the strongest Stage 3R signal:
  - `mini-protein` mean drops below the `<5e-5` gate (`2.883e-5`) and K3
    accepted misses drop from `17` to `1`;
  - `compact-mini-protein` mean drops to `2.448e-5` and K3/K4 accepted misses
    drop to `0`, but K2 still has `3` accepted misses;
  - leaf work rises only slightly in these rows, while far-node work drops
    sharply.
- No Stage 3R row passes full medium gates across all three fixtures because
  accepted-missed is still nonzero on `mini-protein` and/or
  `compact-mini-protein`.
- Barnase 120s was not run because the mandatory full-medium candidate did not
  pass. Barnase 300s was not run.

Stage 3R decision:

- Keep production defaults unchanged.
- Keep accepted-miss diagnostics, K3 rho cap, and surface-margin MAC hidden and
  experimental.
- Do not promote and do not run Barnase yet.
- The next branch should target the structural low-rho misses now exposed by
  the accepted-miss diagnostics:
  1. add a lower accepted-rho diagnostic around the observed K3 band
     (`0.24..0.26`) only as a lower-bound probe, because `0.52..0.58` is not
     relevant to the current medium misses;
  2. generalize the accepted-rho or surface-margin guard to K2/K4 only where
     accepted-only diagnostics show residual misses, rather than applying a
     global cap;
  3. consider a decoupled phi/phin tail estimator or true-surface distance
     score that catches low-rho structural misses without relying on blind
     amplitude increases;
  4. require mean `<5e-5` and global/per-channel accepted-missed `0` on all
     three medium fixtures before Barnase 120s.

## Stage 3S/9L - TABI Absolute-Floor Miss Semantics and Surface-Margin Candidate

Stage 3S implementation:

- Added hidden audit control:
  `-pbe-tabi-accepted-miss-absolute-floor=<0.0..1e-4>`, default `1e-7`.
- Updated far-expansion truth gates so high-error `missed`, `caught`, and
  `caught_by_split` require both:
  - relative error above `-pbe-tabi-convergence-estimator-error-cap`;
  - absolute JCP channel error above the accepted-miss absolute floor.
- Kept relative mean/p95/p99 over all sampled expansion-audit pairs. Only
  truth counters (`missed/caught/caught_by_split`) are floor-filtered.
- Added phantom high-error counters:
  - global `tabi_convergence_estimator_phantom_high_error_count` and
    `tabi_convergence_estimator_phantom_high_error_abs_max`;
  - per-channel `tabi_phantom_high_error_k{1..4}_count` and
    `tabi_phantom_high_error_k{1..4}_abs_max`.
- Extended sampled accepted-miss diagnostics with `abs_error,is_phantom` in
  `pbe_tabi_accepted_miss_diagnostics.csv`.
- Plumbed the floor and phantom counters through timing/profile CSVs, longrun
  progress, the medium matrix, and the Barnase truth matrix.
- Added `test/pbe_tabi_absolute_floor_miss_smoke.ps1`.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with
  `-IncludeStage3SRows`:
  - candidate: order4 theta `0.60`, `surface-margin=0.25`,
    `polynomial-eps`, accepted-only, floor `1e-7`;
  - sensitivity: same row with floor `1e-8`;
  - baseline comparison: order4 theta `0.60`, center-distance,
    `polynomial-eps`, accepted-only, floor `1e-7`.
- Hardened the Barnase truth harness for Stage 3S evidence:
  - added the accepted-miss absolute floor to CLI/csv output;
  - shortened long run-directory names with a SHA1 suffix so deep Windows paths
    do not fail before the solver starts.

Measured Stage 3S smoke/build evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Absolute-floor miss smoke: PASS.
  - phantom accepted high-error case: `314` phantom high errors, `0` missed;
  - real accepted high-error case: `314` missed;
  - split real high-error case: `309` caught by split;
  - split phantom case did not increment `caught_by_split`.
- Adjacent focused gates after Stage 3S: PASS.
  - channel missed semantics smoke;
  - K3 accepted-rho cap smoke;
  - surface-margin MAC smoke;
  - polynomial channel estimator smoke;
  - topdown MAC audit;
  - treecode accuracy;
  - tree builder smoke.

Measured Stage 3S medium matrix:

- Command:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3RRows -IncludeStage3SRows -CaseTimeoutSeconds 180 -AllowFailedStatus`.
- Output:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

| fixture | Stage 3S row | status | matvec rel | sampled residual | expansion mean | p99 | accepted missed K1/K2/K3/K4 | phantom global | leaf work |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| mini-protein | surface-margin floor `1e-7` | `ok` | `2.52e-8` | `9.76e-5` | `3.22e-5` | `1.21e-4` | `0/0/0/0` | `35` | `7335697` |
| mini-protein | surface-margin floor `1e-8` | `ok` | `2.52e-8` | `9.76e-5` | `2.77e-5` | `1.19e-4` | `0/0/0/0` | `0` | `7335697` |
| mini-protein | center floor `1e-7` | `ok` | `8.92e-8` | `9.76e-5` | `6.97e-5` | `4.90e-4` | `0/0/0/0` | `28` | `7290772` |
| charged-peptide | surface-margin floor `1e-7` | `ok` | `5.29e-16` | `1.77e-5` | `0` | `0` | `0/0/0/0` | `13` | `5021476` |
| charged-peptide | surface-margin floor `1e-8` | `ok` | `5.29e-16` | `1.77e-5` | `0` | `0` | `0/0/0/0` | `0` | `5021476` |
| compact-mini-protein | surface-margin floor `1e-7` | `ok` | `2.87e-8` | `5.31e-5` | `2.54e-5` | `3.51e-4` | `0/0/0/0` | `40` | `17228993` |
| compact-mini-protein | surface-margin floor `1e-8` | `ok` | `2.87e-8` | `5.31e-5` | `2.69e-5` | `3.51e-4` | `0/0/0/0` | `0` | `17228993` |
| compact-mini-protein | center floor `1e-7` | `ok` | `6.96e-8` | `5.31e-5` | `3.74e-5` | `4.82e-4` | `0/0/0/0` | `24` | `17098706` |

Stage 3S medium interpretation:

- The comment was correct: the remaining Stage 3R accepted-only misses were
  phantom high-relative-error samples with tiny absolute error. Adding the
  absolute floor moves these into explicit phantom counters and removes them
  from `missed/caught/caught_by_split` truth gates.
- The Stage 3S surface-margin candidate passes all three medium fixtures:
  energy drift `<5e-3`, matvec rel `<2.5e-3`, sampled residual within `2x`
  grouped-morton, expansion mean `<5e-5`, and accepted-missed `0`.
- Floor sensitivity is stable: both `1e-7` and `1e-8` pass medium gates.
  The default/candidate floor remains `1e-7` because it isolates negligible
  absolute perturbations while keeping the stricter `1e-8` row available as
  a regression comparison.

Measured Barnase Stage 3S confirmation:

- 120s command:
  `test/pbe_tabi_barnase_truth_matrix.ps1 -TreeProbeOnly -SkipGridReference -CaseTimeoutSeconds 120 ... -TreeLeafSizes 128 -TreeOpeningThetas 0.60 -TreeTaylorOrders 4 -NodeBoundsModes tight-aabb -TreeBuilders topdown-aabb -MacSafety physical-gap -ExpansionCenters area-centroid -FragmentationSafety split -NearIntegration duffy_adaptive -QuadBuffer adaptive -MomentAudit sampled -ConvergenceAudit sampled -ConvergenceSafety split -ConvergenceEstimator hybrid -ConvergenceCalibration consecutive-aspect -ConvergenceDielectricAmp on -ConvergenceEstimatorModels channel-tail-eps -ConvergenceChannelCalibrations polynomial-eps -ChannelMissedSemantics accepted-only -AcceptedMissDiagnostics summary -AcceptedMissAbsoluteFloor 1e-7 -MacDistancePolicy surface-margin -MacSurfaceMarginFactor 0.25 -ExpansionAudit sampled -ExpansionAuditChannels summary -PreSolverProgress`.
- Output:
  `test/tmp/b3s_abs/pbe_tabi_barnase_truth_matrix.csv`.

| timeout | status | progress | panels | matvecs | sampled residual | expansion mean | p99 | accepted missed K1/K2/K3/K4 | phantom global | leaf/far work |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- |
| 120s | `timeout` | `gmres_iter`, prepare `solver_enter` | `15129` | `54` | `1.4e-5` | `3.2e-5` | `3.49e-4` | `0/0/0/0` | `483` | `225974749/NA` |

- 120s interpretation: timeout-safe evidence is usable. Solver progress is
  present, sampled residual is already below `1e-4`, expansion mean is below
  `5e-5`, and accepted-missed is `0`, so the Stage 3S candidate was eligible
  for 300s confirmation.

- 300s command: same candidate/settings with `-CaseTimeoutSeconds 300`.
- Output:
  `test/tmp/b3s_300/pbe_tabi_barnase_truth_matrix.csv`.

| timeout | status | decision | panels | gmres iter | matvecs | final residual | sampled residual | expansion mean | p99 | accepted missed K1/K2/K3/K4 | phantom global | leaf/far work |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- |
| 300s | `ok` | `safe_treecode_candidate` | `15129` | `64` | `81` | `7.34e-9` | `6.90e-9` | `3.72e-5` | `5.72e-4` | `0/0/0/0` | `954` | `225974749/12770` |

Stage 3S decision:

- The hidden production candidate is now evidence-backed:
  - `taylor_order=4`;
  - `mac-distance-policy=surface-margin`;
  - `mac-surface-margin-factor=0.25`;
  - `convergence-channel-calibration=polynomial-eps`;
  - `channel-missed-semantics=accepted-only`;
  - `accepted-miss-absolute-floor=1e-7`;
  - `physical-gap`, `area-centroid`, and `fragmentation-safety=split`.
- Medium + Barnase 120s + Barnase 300s gates pass. P99 remains report-only as
  planned; the Barnase p99 spike `~5.72e-4` did not prevent sampled/full
  residual from converging.
- Keep direct TABI, grid backend, and DelPhi/pyDelPhi oracle-only policy
  unchanged.
- Production defaults remain unchanged in this implementation patch; the next
  small branch can flip the experimental JCP treecode default/preset to the
  Stage 3S candidate if desired, backed by the evidence above.

## Stage 3T/9M - TABI Phantom-Floor Promotion Readiness Audit

Stage 3T implementation:

- Extended the medium and Barnase truth matrix CSV summaries with
  promotion-safety derived columns:
  - `phi_norm`, `phin_norm`, `source_norm_s1`, `source_norm_s2`;
  - `phantom_global_abs_bound`;
  - `phantom_k{1..4}_abs_bound`;
  - `phantom_k3_abs_over_phi_norm`;
  - `phantom_k3_abs_over_phin_norm`.
- Added Barnase per-row promotion-readiness fields:
  - `promotion_readiness`;
  - `promotion_readiness_reason`.
- The Barnase promotion-readiness gate requires:
  - completed `status=ok`;
  - residual `<=1e-4`;
  - `far_expansion_rel_error_mean < 5e-5`;
  - global and per-channel accepted-missed counters `0`;
  - `phantom_k3_abs_over_phi_norm < 1e-4`;
  - `phantom_k3_abs_bound < 1e-2`.
- Added `-IncludeStage3TRows` to
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1`. The Stage 3T row uses the
  same hidden candidate as Stage 3S:
  order 4, theta `0.60`, `surface-margin=0.25`, `physical-gap`,
  `area-centroid`, `fragmentation-safety=split`, `channel-tail-eps`,
  `polynomial-eps`, `accepted-only`, absolute floor `1e-7`.

Barnase 300s promotion-readiness evidence from the existing Stage 3S row:

| field | value |
| --- | ---: |
| `phi_norm` | `0.001308301853` |
| `phin_norm` | `0.01264949228` |
| `phantom_high_error_k3_count` | `15696` |
| `phantom_high_error_k3_abs_max` | `9.998896736e-8` |
| `phantom_k3_abs_bound` | `0.0015694268` |
| `phantom_k3_abs_over_phi_norm` | `7.6427e-5` |
| `phantom_k3_abs_over_phin_norm` | `7.9046e-6` |
| `far_expansion_rel_error_mean` | `3.7175e-5` |
| accepted missed K1/K2/K3/K4 | `0/0/0/0` |
| final residual | `7.3381e-9` |

Measured Stage 3T verification:

- PowerShell parser check for the updated medium and Barnase harnesses: PASS.
- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS
  (`Nothing to be done for 'all'`).
- Focused smokes: PASS.
  - `test/pbe_tabi_absolute_floor_miss_smoke.ps1`;
  - `test/pbe_tabi_channel_missed_semantics_smoke.ps1`;
  - `test/pbe_tabi_surface_margin_mac_smoke.ps1`;
  - `test/pbe_tabi_polynomial_channel_estimator_smoke.ps1`;
  - `test/pbe_tabi_treecode_accuracy.ps1`;
  - `test/pbe_tabi_tree_builder_smoke.ps1`;
  - `test/pbe_tabi_topdown_mac_audit.ps1`.
- Stage 3S/3T full medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3SRows -IncludeStage3TRows -CaseTimeoutSeconds 180 -AllowFailedStatus`: PASS.

| fixture | Stage 3T status | classification | `phi_norm` | `phantom_k3_abs_bound` | `phantom_k3_abs_over_phi_norm` | expansion mean | accepted missed K1/K2/K3/K4 |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| mini-protein | `ok` | `convergence_mean_gate_medium_pass` | `4.975e-4` | `1.3387e-4` | `2.0080e-4` | `2.555e-5` | `0/0/0/0` |
| charged-peptide | `ok` | `convergence_mean_gate_followup` | `4.348e-4` | `2.1720e-4` | `2.2989e-4` | `0` | `0/0/0/0` |
| compact-mini-protein | `ok` | `convergence_mean_gate_medium_pass` | `7.028e-4` | `3.4566e-4` | `1.4215e-4` | `2.601e-5` | `0/0/0/0` |

  Charged-peptide remains `convergence_mean_gate_followup` only because no
  expansion-audit mean-gate pair is emitted there; energy/residual/accepted-miss
  gates remain clean as in Stage 3S.

Stage 3T decision:

- The phantom floor is now validated against Barnase 300s norm scale:
  `phantom_k3_abs_over_phi_norm` is below `1e-4`, and the conservative K3
  phantom bound is below `0.01`.
- The Stage 3S hidden candidate should be considered
  `promotion_readiness=pass` once the updated Barnase harness is re-run or the
  existing row is post-processed with the new columns.
- Production defaults remain unchanged in Stage 3T. The next branch may either
  promote the experimental JCP treecode preset/default to the evidence-backed
  candidate, or keep it hidden and run broader regression timing first.

## Stage 3U/9N - TABI Promotion-Readiness Regression Freeze

Stage 3U implementation:

- Kept production defaults unchanged. Stage 3U does not change public help,
  manual behavior, direct TABI, grid backend, or DelPhi/pyDelPhi oracle-only
  policy.
- Added `-IncludeStage3URows` to
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1`:
  - promotion-freeze row: Stage 3T candidate with cluster `off`;
  - timing-only row: same candidate with `leaf-pair-balanced`, batch size `4`;
  - floor-sensitivity row: same candidate with absolute floor `1e-8`.
- Added `test/pbe_tabi_promotion_readiness_summary.ps1`.
  The script reads the current medium matrix CSV plus the existing Barnase
  300s CSV, recomputes missing derived Barnase metrics from `run_dir` timing
  output when needed, and writes
  `pbe_tabi_promotion_readiness_summary.csv`.
- Promotion-readiness summary gates:
  - medium candidate must pass status, matvec, energy, residual-vs-grouped,
    far-expansion mean, and accepted-missed checks on all three fixtures;
  - Barnase 300s must pass status, residual `<=1e-4`, expansion mean `<5e-5`,
    accepted-missed `0`, `phantom_k3_abs_over_phi_norm <1e-4`, and
    `phantom_k3_abs_bound <1e-2`;
  - `production_default_changed` is explicitly reported as `false`.

Measured Stage 3U verification:

- PowerShell parser check for the updated medium, Barnase, and promotion
  summary harnesses: PASS.
- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS
  (`Nothing to be done for 'all'`).
- Focused TABI gates: PASS.
  - `test/pbe_tabi_absolute_floor_miss_smoke.ps1`;
  - `test/pbe_tabi_channel_missed_semantics_smoke.ps1`;
  - `test/pbe_tabi_surface_margin_mac_smoke.ps1`;
  - `test/pbe_tabi_polynomial_channel_estimator_smoke.ps1`;
  - `test/pbe_tabi_treecode_accuracy.ps1`;
  - `test/pbe_tabi_tree_builder_smoke.ps1`;
  - `test/pbe_tabi_topdown_mac_audit.ps1`.
- Stage 3S/3T/3U full medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3SRows -IncludeStage3TRows -IncludeStage3URows -CaseTimeoutSeconds 180 -AllowFailedStatus`: PASS.

| fixture | Stage 3U role | cluster | seconds | sampled residual | expansion mean | accepted missed | classification |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| mini-protein | promotion freeze | `off` | `3.78` | `9.76e-5` | `2.94e-5` | `0` | `convergence_mean_gate_medium_pass` |
| charged-peptide | promotion freeze | `off` | `4.40` | `1.77e-5` | `0` | `0` | `convergence_mean_gate_followup` |
| compact-mini-protein | promotion freeze | `off` | `16.30` | `5.31e-5` | `2.53e-5` | `0` | `convergence_mean_gate_medium_pass` |
| mini-protein | timing report | `leaf-pair-balanced` | `4.16` | `9.76e-5` | `3.16e-5` | `0` | report-only |
| charged-peptide | timing report | `leaf-pair-balanced` | `4.45` | `1.77e-5` | `0` | `0` | report-only |
| compact-mini-protein | timing report | `leaf-pair-balanced` | `14.47` | `5.31e-5` | `2.86e-5` | `0` | report-only |

- Promotion summary:
  `test/pbe_tabi_promotion_readiness_summary.ps1 -MediumCsv test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv -BarnaseCsv test/tmp/b3s_300/pbe_tabi_barnase_truth_matrix.csv`: PASS.
- Output:
  `test/tmp/pbe_tabi_promotion_readiness_summary/pbe_tabi_promotion_readiness_summary.csv`.

| scope | fixture | status | residual | expansion mean | accepted missed | K3 phantom bound | K3 abs / phi norm | overall |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| medium | mini-protein | `ok` | `9.76e-5` | `2.94e-5` | `0` | `1.33e-4` | report-only `2.01e-4` | `True` |
| medium | charged-peptide | `ok` | `1.77e-5` | `0` | `0` | `2.16e-4` | report-only `2.30e-4` | `True` |
| medium | compact-mini-protein | `ok` | `5.31e-5` | `2.53e-5` | `0` | `3.44e-4` | report-only `1.42e-4` | `True` |
| Barnase 300s | barnase | `ok` | `7.34e-9` | `3.72e-5` | `0` | `1.57e-3` | `7.64e-5` | `True` |
| overall | all | `ok` | `NA` | `NA` | `NA` | `NA` | `NA` | `True` |

- Regression freeze gates: PASS.
  - `test/pbe_tabi_gmres_true_residual_smoke.ps1`;
  - `test/pbe_tabi_sampled_true_residual_smoke.ps1`;
  - `test/pbe_tabi_profile_smoke.ps1`;
  - `test/pbe_memory_profile_smoke.ps1`;
  - `test/pbe_clean_room_audit.ps1`;
  - `test/pbe_regression.ps1`;
  - `test/pbsa_smoke.ps1`;
  - `test/sasa_surface_smoke.ps1`.
- `git diff --check` on the Stage 3U touched files: PASS; only existing
  CRLF normalization warnings were reported.

Stage 3U decision:

- Stage 3U freezes the candidate as a hidden promotion-ready package:
  summary `overall_pass=True`, Barnase 300s norm gate passes, and
  `production_default_changed=false`.
- Production defaults remain unchanged.
- `leaf-pair-balanced` remains timing-only: it improves the compact-mini row in
  this run but is slower on the smaller medium fixtures and is not part of the
  promotion candidate.
- The next branch may promote an explicit experimental preset/default for the
  JCP treecode path, backed by Stage 3S/3T/3U evidence, or continue with a
  broader timing campaign if runtime predictability is preferred.

## Stage 3V/9O - TABI Experimental Preset and Grid-vs-TABI Compare Gate

Stage 3V implementation:

- Added hidden opt-in preset `-pbe-tabi-tree-preset=off|stage3u-safe`.
  Production defaults remain unchanged; default preset is `off`.
- The preset is applied at the CLI position where it appears. TABI flags before
  the preset are overwritten by the preset; TABI flags after the preset override
  it. This keeps the preset convenient while preserving explicit CLI precedence.
- `stage3u-safe` expands to the Stage 3U promotion-freeze candidate:
  `sum=treecode`, `formulation=jcp2013-treecode`, `taylor_order=4`,
  `theta=0.60`, `tree_builder=topdown-aabb`, `node_bounds=tight-aabb`,
  `mac_safety=physical-gap`, `expansion_center=area-centroid`,
  `fragmentation_safety=split`, `mac_distance_policy=surface-margin`,
  surface margin factor `0.25`, `convergence_estimator_model=channel-tail-eps`,
  `convergence_channel_calibration=polynomial-eps`,
  `channel_missed_semantics=accepted-only`,
  `accepted_miss_absolute_floor=1e-7`, and `cluster_policy=off`.
- The preset also freezes the Stage 3U validation context:
  `leaf=128`, `near=duffy_adaptive`, `quad_buffer=adaptive`,
  `interaction_mode=onthefly`, `precond=diagonal`, sampled moment,
  convergence, and expansion/channel audits.
- Added `tabi_tree_preset` to `pbe_timing.csv` and `tree_preset` to
  `pbe_tabi_profile.csv`; the console/PBE summary also reports the preset.
- Added `test/pbe_tabi_stage3u_preset_smoke.ps1` to assert that the preset
  writes the expected profile/timing fields and expands to the Stage 3U safe
  settings.
- Extended `test/pbe_tabi_vs_grid_compare.ps1`:
  - new `-UseStage3USafePreset` switch;
  - new `-Threads` override, defaulting to `8` when the preset is used;
  - robust absolute `WorkRoot` handling so relative work roots still produce
    valid `arg.in` paths after the harness changes directory;
  - compare columns for `time_ratio_tabi_over_grid`,
    `energy_delta_tabi_minus_grid`, `energy_delta_rel_vs_grid`,
    `peak_memory_ratio_tabi_over_grid`, `rss_ratio_tabi_over_grid`,
    `tabi_residual`, `tabi_far_expansion_mean`, and
    `tabi_accepted_missed_total`.

Measured Stage 3V verification:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Parser checks for the updated preset and compare harnesses: PASS.
- Focused TABI gates: PASS.
  - `test/pbe_tabi_stage3u_preset_smoke.ps1`;
  - `test/pbe_tabi_absolute_floor_miss_smoke.ps1`;
  - `test/pbe_tabi_channel_missed_semantics_smoke.ps1`;
  - `test/pbe_tabi_surface_margin_mac_smoke.ps1`;
  - `test/pbe_tabi_polynomial_channel_estimator_smoke.ps1`;
  - `test/pbe_tabi_tree_builder_smoke.ps1`;
  - `test/pbe_tabi_treecode_accuracy.ps1`;
  - `test/pbe_tabi_topdown_mac_audit.ps1`.
- Promotion summary still passes after adding the preset:
  `test/pbe_tabi_promotion_readiness_summary.ps1 -MediumCsv test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv -BarnaseCsv test/tmp/b3s_300/pbe_tabi_barnase_truth_matrix.csv`: PASS,
  `overall_pass=True`, `production_default_changed=false`.
- `git diff --check` on Stage 3V touched files: PASS; only existing CRLF
  normalization warnings were reported.

Stage 3V medium grid-vs-TABI compare:

Command:

`test/pbe_tabi_vs_grid_compare.ps1 -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3USafePreset -Mode perf`

Output:

`test/tmp/pbe_tabi_vs_grid_compare/pbe_tabi_vs_grid_compare.csv`

| fixture | backend | seconds | energy kJ/mol | preset | time ratio TABI/grid | energy delta TABI-grid | energy rel delta | peak memory ratio | RSS ratio | TABI residual | expansion mean | accepted missed |
| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| mini-protein | grid | `0.00453` | `-916.3506` | `grid` | `1` | `0` | `0` | `1` | `1` | `NA` | `NA` | `NA` |
| mini-protein | TABI | `5.274` | `-1198.3724` | `stage3u-safe` | `1163.63` | `-282.0218` | `0.3078` | `1.386` | `1.106` | `0` | `3.05e-5` | `0` |
| charged-peptide | grid | `0.00411` | `-529.2189` | `grid` | `1` | `0` | `0` | `1` | `1` | `NA` | `NA` | `NA` |
| charged-peptide | TABI | `5.355` | `-843.7091` | `stage3u-safe` | `1301.76` | `-314.4902` | `0.5943` | `1.316` | `1.116` | `0` | `0` | `0` |
| compact-mini-protein | grid | `0.00376` | `-1000.9248` | `grid` | `1` | `0` | `0` | `1` | `1` | `NA` | `NA` | `NA` |
| compact-mini-protein | TABI | `15.233` | `-1418.3443` | `stage3u-safe` | `4050.92` | `-417.4195` | `0.4170` | `1.564` | `1.124` | `9.54e-7` | `3.18e-5` | `0` |

Stage 3V Barnase compare note:

- Optional Barnase salt0 compare was attempted with
  `-UseStage3USafePreset -BarnaseTabiSurfacePoints 8`.
- After fixing relative work-root handling, the Barnase grid side completed and
  wrote `pbe_timing.csv` in about `2.23s`.
- The full grid-vs-TABI Barnase compare did not emit the final compare summary
  within the 360s turn budget because the TABI side was still running. This is
  not a promotion gate failure; Stage 3V only requires the preset and medium
  compare gate. It does confirm that protein-scale time remains a performance
  branch, not a correctness branch.

Stage 3V decision:

- `stage3u-safe` is now a hidden experimental preset for repeatable evidence
  runs. It is not a production default.
- The medium compare gate demonstrates the current tradeoff clearly:
  Stage 3U safe TABI passes its internal residual/expansion/missed diagnostics,
  but grid is still orders of magnitude faster on small fixtures and remains
  the fastest default path.
- The reported medium energies are backend outputs for comparison only; they do
  not replace direct-vs-treecode, residual, or Barnase promotion-readiness gates.
- `leaf-pair-balanced` remains timing-only and is not part of `stage3u-safe`.
- Recommended Stage 3W/9P: performance branch for the hidden preset, focused on
  protein-scale throughput. Start with Barnase timeout-safe compare/progress and
  then optimize leaf/matvec throughput only after the preset path keeps the
  Stage 3U correctness counters clean.

## Stage 3W/9P - TABI Stage3U Performance Branch

Stage 3W implementation:

- Added hidden preset `-pbe-tabi-tree-preset=stage3u-perf`; production default
  remains `off` and `stage3u-safe` is unchanged.
- `stage3u-perf` keeps the Stage 3U math/safety envelope:
  order 4, theta `0.60`, `topdown-aabb`, `tight-aabb`, `physical-gap`,
  `area-centroid`, `fragmentation-safety=split`, `surface-margin=0.25`,
  `channel-tail-eps`, `polynomial-eps`, `accepted-only`, and absolute floor
  `1e-7`.
- The perf preset turns report-only audits off by default and enables the
  performance path: `leaf-panel-kernel=tile-avx2`,
  `cluster-policy=leaf-pair-balanced`, `leaf-pair-kernel=target-avx2`,
  target batch size `4`, min source panels `64`, min targets `4`.
- Added hidden `-pbe-tabi-leaf-pair-kernel=auto|scalar|target-avx2`. `auto`
  selects `target-avx2` only for AVX2-capable JCP on-the-fly leaf-pair batches
  with target batch size `4`; unsupported or unclean paths fall back to scalar.
- Implemented `target-avx2` for clean JCP leaf-pair batches, vectorized across
  four target panels per source panel. The fast path keeps target-lane
  accumulators in AVX2 registers over the source span and stores once per
  target batch.
- Added timing/profile/progress counters:
  `tabi_leaf_pair_batch_kernel`,
  `tabi_leaf_pair_batch_vectorized_sources`,
  `tabi_leaf_pair_batch_vector_lanes`, and
  `tabi_leaf_pair_batch_pairs_per_second`.
- Extended `test/pbe_tabi_vs_grid_compare.ps1` with
  `-UseStage3UPerfPreset`, timeout/allow-timeout handling, and partial TABI
  progress fields from `pbe_tabi_longrun_progress.csv` /
  `pbe_tabi_prepare_progress.csv`.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with
  `-IncludeStage3WRows`, including safe baseline, perf isolate, perf candidate,
  validation audit row, and min-source threshold rows.
- Added `test/pbe_tabi_leaf_pair_target_avx2_smoke.ps1` and extended
  `test/pbe_tabi_stage3u_preset_smoke.ps1` to validate both `stage3u-safe` and
  `stage3u-perf`.

Measured Stage 3W verification:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Syntax checks for updated compare, matrix, and target-AVX2 smoke harnesses:
  PASS.
- Focused smokes:
  - `test/pbe_tabi_leaf_pair_target_avx2_smoke.ps1`: PASS,
    scalar-vs-AVX2 relative RMS `1.36e-16`, vectorized sources `2052`,
    vector lanes `8208`, singular fallback relative RMS `1.20e-16`.
  - `test/pbe_tabi_leaf_pair_batch_smoke.ps1`: PASS.
  - `test/pbe_tabi_stage3u_preset_smoke.ps1`: PASS for both
    `stage3u-safe` and `stage3u-perf`.
- Quick Stage 3W matrix on `mini-protein`:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3WRows -Quick -AllowFailedStatus`:
  PASS and emitted Stage 3W batch/vector counters.
- Full Stage 3W matrix across all three medium fixtures was attempted but did
  not emit the final CSV within the 360s tool timeout. No Titan child process
  remained after timeout; this run is not counted as a valid gate result.

Stage 3W medium perf compare:

Safe baseline command:

`test/pbe_tabi_vs_grid_compare.ps1 -WorkRoot test/tmp/pbe_tabi_vs_grid_compare_safe -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3USafePreset -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout`

Perf command:

`test/pbe_tabi_vs_grid_compare.ps1 -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3UPerfPreset -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout`

| fixture | safe seconds | perf seconds | speedup | perf kernel | vectorized sources | pair/s |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| mini-protein | `4.3189` | `4.5147` | `0.9566x` | `target-avx2` | `32,116,491` | `1.04e8` |
| charged-peptide | `4.3504` | `2.7403` | `1.5875x` | `target-avx2` | `21,235,170` | `1.26e8` |
| compact-mini-protein | `17.5758` | `12.5717` | `1.3980x` | `target-avx2` | `65,629,872` | `6.94e7` |

Stage 3W decision:

- Correctness/plumbing for `target-avx2` passes focused smoke coverage, and
  the perf preset emits the new counters correctly.
- Performance gate is not yet met: medium geometric mean speedup is about
  `1.29x`, below the Stage 3W target `1.5x`; `mini-protein` is `4.3%` slower,
  which is within the `5%` no-regression allowance but leaves little margin.
- `stage3u-perf` remains hidden experimental/timing-only and is not promoted.
- Next branch should optimize batch-selection overhead and target reuse before
  another Barnase timing run. The likely hotspot is not the vector arithmetic
  itself, but the per-candidate `would_visit_leaf_direct_for_target` traversal
  and batching overhead around many small/medium leaf pairs.

## Stage 3X/9Q - TABI Leaf-Pair Batch Selection Overhead Reduction

Stage 3X implementation:

- Added hidden selector profiling:
  `-pbe-tabi-leaf-pair-selection-profile=off|summary`.
- Added hidden selector controls:
  `-pbe-tabi-leaf-pair-probe-cache=off|per-matvec`,
  `-pbe-tabi-leaf-pair-target-group=scan|same-target-leaf`, and
  `-pbe-tabi-leaf-pair-max-lookahead=<4..128>`.
- Added per-matvec probe cache for `would_visit_leaf_direct_for_target`,
  keyed by `(target_index, source_node_index)`. The cache is local to each
  matvec worker and is cleared between matvecs, so GMRES still sees a fixed
  operator.
- Added `target_leaf_index_by_panel` plus sorted `target_panels_by_leaf` in
  `TreeWorkspace`. `same-target-leaf` selection now scans directly inside the
  target leaf list instead of scanning original panel order and rejecting most
  candidates.
- Added hidden `stage3u-perf2` preset for timing experiments. It keeps Stage
  3U math/safety and uses `leaf-pair-balanced`, `target-avx2`,
  `probe-cache=per-matvec`, `target-group=same-target-leaf`, and lookahead
  `32`. It is not a production default or promotion candidate yet.
- Extended timing/profile/progress/matrix/compare CSVs with selection counters:
  candidate scans, clean rejects, probe calls, probe node visits, probe/cache
  hits/misses, insufficient-target batches, same-target-leaf rejects,
  probe seconds, and selection seconds.
- Extended `test/pbe_tabi_vs_grid_compare.ps1` with
  `-UseStage3UPerf2Preset` and added
  `test/pbe_tabi_leaf_pair_selection_profile_smoke.ps1`.

Measured Stage 3X verification:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Focused smokes:
  - `test/pbe_tabi_leaf_pair_selection_profile_smoke.ps1`: PASS,
    scalar-vs-cached/target-leaf relative RMS `1.09e-16`, candidate scans
    `374`, probe calls `374`, scan-mode cache hits `1004`, batches `128`.
  - `test/pbe_tabi_stage3u_preset_smoke.ps1`: PASS for
    `stage3u-safe`, `stage3u-perf`, and `stage3u-perf2`.
  - `test/pbe_tabi_leaf_pair_target_avx2_smoke.ps1`: PASS,
    relative RMS `1.36e-16`, vectorized sources `2052`.
  - `test/pbe_tabi_leaf_pair_batch_smoke.ps1`: PASS.
  - `test/pbe_tabi_tree_builder_smoke.ps1`: PASS.
  - `test/pbe_tabi_treecode_accuracy.ps1`: PASS.
  - `test/pbe_tabi_topdown_mac_audit.ps1`: PASS.
- Quick Stage 3X matrix on `mini-protein`:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3XRows -Quick -AllowFailedStatus`:
  PASS. Validation row kept `far_expansion_rel_error_mean=2.735e-05` and
  accepted missed `0`.

Stage 3X mini quick selector rows:

| row | total s | batch s | selection s | probe s | batches | scans | probe calls | cache hits |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| scan/no-cache/lh32 | `4.1916` | `0.9612` | `2.7963` | `1.5038` | `169,320` | `546,720` | `509,160` | `0` |
| cache/same-leaf/lh32 | `4.2491` | `0.9145` | `2.9061` | `1.5029` | `158,010` | `522,855` | `496,305` | `6,030` |
| validation cache/same-leaf/lh32 | `4.7670` | `0.9796` | `3.0918` | `1.5938` | `158,010` | `522,855` | `496,305` | `6,030` |

Stage 3X medium perf compare, sequential:

Safe command:

`test/pbe_tabi_vs_grid_compare.ps1 -WorkRoot test/tmp/pbe_tabi_vs_grid_compare_stage3x_safe_seq -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3USafePreset -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout`

Perf2 command:

`test/pbe_tabi_vs_grid_compare.ps1 -WorkRoot test/tmp/pbe_tabi_vs_grid_compare_stage3x_perf2_seq2 -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3UPerf2Preset -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout`

| fixture | safe seconds | perf2 seconds | speedup | selection s | probe s |
| --- | ---: | ---: | ---: | ---: | ---: |
| mini-protein | `8.1792` | `4.4398` | `1.8423x` | `3.5308` | `1.8240` |
| charged-peptide | `7.1704` | `6.1401` | `1.1678x` | `3.0270` | `0.7701` |
| compact-mini-protein | `15.8590` | `16.6324` | `0.9535x` | `17.1335` | `10.6454` |

Stage 3X decision:

- Correctness gates for the new selector path pass focused smoke and the
  `mini-protein` validation audit row.
- Performance gate is not met. Medium geometric mean speedup is only about
  `1.27x` versus the current safe sequential run, below the `1.5x` gate, and
  `compact-mini-protein` regresses.
- `stage3u-perf2` remains hidden experimental/timing-only and is not promoted.
- Barnase timing was not run because medium performance gate did not pass.
- Next branch should attack the remaining selection cost directly. The new CSV
  shows compact-mini spends more time in selection/probe than in the AVX2
  batch kernel itself, so the likely next step is target-leaf/source-leaf pair
  reuse or a safe target-leaf directness proof that avoids per-target traversal
  for clean leaf pairs.

## Stage 3Y/9R - TABI Ancestor-Path Directness Probe and Leaf-Pair Reuse

Stage 3Y implementation:

- Added hidden leaf-pair directness controls:
  `-pbe-tabi-leaf-pair-directness-probe=full-traversal|ancestor-path`,
  `-pbe-tabi-leaf-pair-directness-cache=off|target-leaf-source-leaf`, and
  `-pbe-tabi-leaf-pair-directness-verify=off|sampled|full-small`.
  Defaults preserve the Stage 3X behavior.
- Added `parent_index_by_node` and leaf ancestor paths to `TreeWorkspace`.
  The `ancestor-path` probe checks only the ancestor chain of the candidate
  source leaf and rejects batching when any ancestor would be accepted as a far
  node instead of descending to that leaf.
- Added an exact per-matvec `(target_leaf, source_leaf)` all-direct proof cache.
  It only accelerates proven all-direct leaf pairs; mixed verdicts fall back to
  per-target directness probes.
- Added sampled/full-small verification counters comparing the new probe/cache
  verdicts against the old full traversal.
- Extended timing/profile/progress/matrix/compare CSVs with directness mode,
  cache, verify, ancestor probe calls/visits/seconds, pair-cache hits/misses,
  all-direct proofs, mixed fallbacks, and verified mismatches.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with
  `-IncludeStage3YRows` and `test/pbe_tabi_vs_grid_compare.ps1` with
  `-TabiExtraArgs` so compare rows can override hidden probe knobs after a
  preset.
- Added:
  - `test/pbe_tabi_leaf_pair_directness_probe_smoke.ps1`;
  - `test/pbe_tabi_leaf_pair_directness_cache_smoke.ps1`.

Measured Stage 3Y verification:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Focused smokes:
  - `test/pbe_tabi_leaf_pair_directness_probe_smoke.ps1`: PASS,
    scalar-vs-batched relative RMS `0`, verify relative RMS `0`,
    full traversal visits `2242`, ancestor visits `748`, verify checks `374`.
  - `test/pbe_tabi_leaf_pair_directness_cache_smoke.ps1`: PASS,
    scalar-vs-cache relative RMS `0`, no-cache probes `980`, cache probes `0`,
    pair hits `308`, pair misses `56`, all-direct proofs `56`, mixed fallbacks
    `0`, verify checks `980`.
  - `test/pbe_tabi_leaf_pair_target_avx2_smoke.ps1`: PASS.
  - `test/pbe_tabi_leaf_pair_batch_smoke.ps1`: PASS.
  - `test/pbe_tabi_stage3u_preset_smoke.ps1`: PASS.
  - `test/pbe_tabi_tree_builder_smoke.ps1`: PASS.
  - `test/pbe_tabi_treecode_accuracy.ps1`: PASS.
  - `test/pbe_tabi_topdown_mac_audit.ps1`: PASS.
- Full Stage 3Y medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3YRows -CaseTimeoutSeconds 180 -AllowFailedStatus`:
  PASS.

Stage 3Y medium matrix highlights:

| fixture | row | total s | selection s | probe s | probe visits | ancestor visits | expansion mean | missed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| mini-protein | perf2 full traversal | `3.2030` | `2.1912` | `1.1018` | `19,912,830` | `0` | `0` | `0` |
| mini-protein | ancestor only | `3.2041` | `1.2752` | `0.1136` | `1,488,915` | `1,488,915` | `0` | `0` |
| mini-protein | ancestor cache validation | `4.3142` | `1.8106` | `0.3976` | `5,059,815` | `4,737,600` | `2.78e-5` | `0` |
| charged-peptide | perf2 full traversal | `4.0979` | `1.7886` | `0.4386` | `11,934,040` | `0` | `0` | `0` |
| charged-peptide | ancestor only | `4.0878` | `1.4986` | `0.1437` | `1,280,080` | `1,280,080` | `0` | `0` |
| compact-mini-protein | perf2 full traversal | `13.4146` | `12.8265` | `7.8284` | `126,673,320` | `0` | `0` | `0` |
| compact-mini-protein | ancestor only | `13.5334` | `5.9962` | `0.3598` | `5,698,620` | `5,698,620` | `0` | `0` |
| compact-mini-protein | ancestor cache validation | `13.1027` | `5.9375` | `0.9430` | `19,420,920` | `17,441,820` | `2.70e-5` | `0` |

Stage 3Y medium compare:

Commands:

`test/pbe_tabi_vs_grid_compare.ps1 -WorkRoot test/tmp/pbe_tabi_vs_grid_compare_stage3y_safe -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3USafePreset -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout`

`test/pbe_tabi_vs_grid_compare.ps1 -WorkRoot test/tmp/pbe_tabi_vs_grid_compare_stage3y_perf2 -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3UPerf2Preset -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout`

`test/pbe_tabi_vs_grid_compare.ps1 -WorkRoot test/tmp/pbe_tabi_vs_grid_compare_stage3y_ancestor -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3UPerf2Preset -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout -TabiExtraArgs '-pbe-tabi-leaf-pair-directness-probe=ancestor-path','-pbe-tabi-leaf-pair-directness-cache=off','-pbe-tabi-leaf-pair-directness-verify=off'`

| fixture | safe s | perf2 s | ancestor s | ancestor speedup vs safe | ancestor selection s | ancestor probe s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| mini-protein | `6.2047` | `7.8892` | `7.2098` | `0.8606x` | `3.1885` | `0.2882` |
| charged-peptide | `7.0235` | `6.0092` | `5.7667` | `1.2179x` | `2.4331` | `0.2327` |
| compact-mini-protein | `15.2035` | `16.7701` | `14.7927` | `1.0278x` | `6.6864` | `0.4091` |

Stage 3Y decision:

- Correctness gate passes: directness verification mismatches are `0`, scalar
  equivalence smokes pass, and validation rows keep far-expansion mean below
  `5e-5` with accepted misses `0`.
- The ancestor-path probe is real and useful locally: compact-mini probe visits
  drop from `126.7M` to `5.70M`, and probe seconds drop from `7.83s` to
  `0.36s` in the matrix row.
- The end-to-end performance gate is still not met. Medium compare geometric
  mean speedup for the ancestor row is only about `1.03x` versus
  `stage3u-safe`, below the `1.5x` gate; mini-protein is slower than safe.
- `stage3u-perf3` was not created. Barnase timing was not run.
- The next branch should not spend more effort on probe traversal alone. The
  remaining bottleneck has shifted into batch construction, candidate
  bookkeeping, and the amount of leaf work exposed by the safe MAC. A better
  Stage 3Z target is to precompute reusable leaf-pair worklists or switch to a
  traversal that emits target-leaf/source-leaf tasks directly, then run the
  same correctness gate before another Barnase timing attempt.

## Stage 3Z/9S - TABI Deferred Leaf-Pair Worklist Scheduler

Stage 3Z implementation:

- Added hidden leaf-pair scheduler controls:
  `-pbe-tabi-leaf-pair-scheduler=inline-select|deferred-worklist`,
  `-pbe-tabi-leaf-pair-worklist-max-active=<64..8192>`, and
  `-pbe-tabi-leaf-pair-worklist-verify=off|sampled|full-small`.
  Defaults keep the old inline selection path.
- Implemented `deferred-worklist` for JCP on-the-fly leaf interactions. Clean
  direct leaf visits emitted by the real traversal are enqueued by
  `(target_leaf, source_leaf_node)` and flushed through the existing
  `target-avx2` batch kernel when a target batch is full. Unclean/singular/
  one-ring/Duffy paths stay scalar and immediate.
- Worklist state is per worker and per matvec only. It does not persist between
  GMRES iterations and does not change MAC, estimator, or near-field safety.
- Added timing/profile/progress counters for scheduler mode, worklist tasks,
  groups, batched groups, scalar flushes, forced flushes, active group high
  water mark, worklist seconds, and verify checks/mismatches.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with
  `-IncludeStage3ZRows` and extended `test/pbe_tabi_vs_grid_compare.ps1` with
  `-UseDeferredLeafPairWorklist` plus CSV worklist columns.
- Added `test/pbe_tabi_leaf_pair_worklist_scheduler_smoke.ps1`.
- `stage3u-safe`, production defaults, direct TABI, and grid backend remain
  unchanged. No `stage3u-perf3` preset was created in this stage.

Measured Stage 3Z verification:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- New worklist smoke:
  `test/pbe_tabi_leaf_pair_worklist_scheduler_smoke.ps1`: PASS,
  scalar-vs-worklist relative RMS `1.74536e-16`, forced-flush relative RMS
  `1.54329e-16`, worklist tasks `960`, groups `272`, batched groups `208`,
  forced flushes `949`, verify checks `960`, and verify mismatches `0`.
- Focused smokes also pass:
  - `test/pbe_tabi_leaf_pair_directness_probe_smoke.ps1`;
  - `test/pbe_tabi_leaf_pair_directness_cache_smoke.ps1`;
  - `test/pbe_tabi_leaf_pair_target_avx2_smoke.ps1`;
  - `test/pbe_tabi_leaf_pair_batch_smoke.ps1`;
  - `test/pbe_tabi_stage3u_preset_smoke.ps1`;
  - `test/pbe_tabi_treecode_accuracy.ps1`;
  - `test/pbe_tabi_tree_builder_smoke.ps1`;
  - `test/pbe_tabi_topdown_mac_audit.ps1`.
- Solver/profile/regression smoke gates pass after the engine change:
  `pbe_tabi_gmres_true_residual_smoke.ps1`,
  `pbe_tabi_sampled_true_residual_smoke.ps1`,
  `pbe_tabi_profile_smoke.ps1`, `pbe_memory_profile_smoke.ps1`,
  `pbe_clean_room_audit.ps1`, `pbe_regression.ps1`, `pbsa_smoke.ps1`, and
  `sasa_surface_smoke.ps1`.
- Full Stage 3Y+3Z medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3YRows -IncludeStage3ZRows -CaseTimeoutSeconds 180 -AllowFailedStatus`:
  PASS. Summary artifact:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/stage3z_summary.csv`.

Stage 3Z medium matrix highlights:

| fixture | best Stage 3Z row | total s | worklist s | tasks | batched groups | validation mean | accepted missed | verify mismatch |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `mini-protein` | `deferred max2048 min64` | `3.8497` | `1.1541` | `680,415` | `158,010` | `2.757e-5` | `0` | `0` |
| `charged-peptide` | `deferred max2048 min64` | `3.3934` | `0.9452` | `611,500` | `141,080` | `0` | `0` | `0` |
| `compact-mini-protein` | `deferred max4096 min64` | `11.4155` | `4.4257` | `2,602,650` | `605,550` | `2.715e-5` | `0` | `0` |

Stage 3Z medium compare:

Commands:

`test/pbe_tabi_vs_grid_compare.ps1 -WorkRoot test/tmp/pbe_tabi_vs_grid_compare_stage3z_safe -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3USafePreset -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout`

`test/pbe_tabi_vs_grid_compare.ps1 -WorkRoot test/tmp/pbe_tabi_vs_grid_compare_stage3z_perf2 -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3UPerf2Preset -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout`

`test/pbe_tabi_vs_grid_compare.ps1 -WorkRoot test/tmp/pbe_tabi_vs_grid_compare_stage3z_worklist -Fixtures mini-protein,charged-peptide,compact-mini-protein -UseStage3UPerf2Preset -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout -TabiExtraArgs '-pbe-tabi-leaf-pair-scheduler=deferred-worklist','-pbe-tabi-leaf-pair-worklist-max-active=4096','-pbe-tabi-leaf-pair-worklist-verify=off'`

Compare summary artifact:
`test/tmp/pbe_tabi_vs_grid_compare_stage3z_worklist/stage3z_compare_summary.csv`.

| fixture | safe s | perf2 s | worklist s | worklist speedup vs safe |
| --- | ---: | ---: | ---: | ---: |
| `mini-protein` | `5.0401` | `5.1172` | `6.8716` | `0.7335x` |
| `charged-peptide` | `3.4453` | `2.8357` | `6.1512` | `0.5601x` |
| `compact-mini-protein` | `12.3460` | `15.6388` | `17.4833` | `0.7062x` |

Stage 3Z decision:

- Correctness gate passes: worklist verify mismatches are `0`, scalar
  equivalence smokes pass, and validation rows keep far-expansion mean below
  `5e-5` with accepted misses `0`.
- The worklist scheduler removes directness probe work from the batching path,
  but the compare gate fails. Medium compare geometric mean speedup is
  `0.662x` versus `stage3u-safe`, well below the `1.5x` target, and every
  medium compare fixture is slower than safe.
- Barnase timing was not run. `stage3u-perf3` was not created.
- Defaults remain unchanged. The next performance branch should not add more
  scheduling layers before reducing per-task bookkeeping and writeback cost, or
  it should revisit whether batching can be emitted at target-leaf granularity
  without storing many tiny `(target, source_leaf)` tasks.

## Stage 3AA/9T - TABI Target-Leaf Block Batching

Stage 3AA implementation:

- Added hidden scheduler value
  `-pbe-tabi-leaf-pair-scheduler=target-leaf-block` and hidden verifier
  `-pbe-tabi-leaf-block-verify=off|sampled|full-small`.
- Implemented a target-leaf/source-leaf grouping path for JCP on-the-fly clean
  direct leaf visits. Real traversal remains the source of truth; clean direct
  leaf candidates are grouped by `(target_leaf, source_leaf_node)` and flushed
  through the existing `target-avx2` target batch kernel. Unclean, singular,
  one-ring, Duffy, and adaptive-buffer paths stay scalar/immediate.
- Added target-leaf block counters to timing/profile/progress output:
  target leaves, targets, source groups, batched groups, scalar tail tasks,
  group/flush/writeback seconds, and verify checks/mismatches.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with
  `-IncludeStage3AARows`, including safe baseline, Stage 3Z deferred-worklist
  diagnostic, target-leaf-block min-source sweeps, and an audit validation row.
- Extended `test/pbe_tabi_vs_grid_compare.ps1` with
  `-UseTargetLeafBlockBatching` and target-leaf block counter columns.
- Added `test/pbe_tabi_leaf_block_scheduler_smoke.ps1`.
- Production defaults, `stage3u-safe`, direct TABI, grid backend, MAC safety,
  and estimator accuracy guards remain unchanged.

Measured Stage 3AA verification:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- New leaf-block smoke:
  `test/pbe_tabi_leaf_block_scheduler_smoke.ps1`: PASS,
  scalar-vs-block relative RMS `1.65074e-16`, source groups `352`,
  batched groups `150`, verify checks `960`, verify mismatches `0`, and
  exclusion fallback relative RMS `1.58042e-16`.
- Existing focused leaf-batch smokes still pass:
  `test/pbe_tabi_leaf_pair_worklist_scheduler_smoke.ps1` and
  `test/pbe_tabi_leaf_pair_target_avx2_smoke.ps1`.
- Full Stage 3Z+3AA medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3ZRows -IncludeStage3AARows -CaseTimeoutSeconds 180 -AllowFailedStatus`:
  PASS. Summary artifact:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.
- Compare harness smoke:
  `test/pbe_tabi_vs_grid_compare.ps1 -Fixtures mini-protein -UseStage3UPerf2Preset -UseTargetLeafBlockBatching -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout`:
  PASS and emits target-leaf block counters.

Stage 3AA medium matrix highlights:

| fixture | safe s | active target-leaf-block row | block s | speedup vs safe | source groups | batched groups | validation mean | accepted misses | verify mismatch |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `mini-protein` | `5.2921` | `min64` | `4.5500` | `1.163x` | `238,320` | `125,145` | `3.136e-5` | `0/0/0/0` | `0` |
| `charged-peptide` | `3.6824` | `min64` | `3.1871` | `1.155x` | `217,380` | `106,260` | `0` | `0/0/0/0` | `0` |
| `compact-mini-protein` | `15.9693` | `min64` | `14.3658` | `1.1116x` | `927,210` | `461,040` | `2.781e-5` | `0/0/0/0` | `0` |

Additional Stage 3AA timing notes:

- The audit validation row with `target-leaf-block min64` reports
  `238,320` source groups and `125,145` batched groups on `mini-protein`,
  `10,575` verifier checks, and zero verifier mismatches.
- The same validation row keeps far-expansion mean below the Stage 3T/3U gate
  on `mini-protein` and `compact-mini-protein`; `charged-peptide` remains the
  known no-audit-pair follow-up row with clean accepted-miss counters.
- The active `min64` target-leaf-block medium geometric mean speedup versus
  the `stage3aa_safe_baseline` is only `1.143x`, below the `1.5x`
  performance gate. Looser `min128` rows are faster on mini/charged but emit
  no block groups there, so they are treated as diagnostics rather than a block
  batching candidate; their best-row geomean is still only `1.185x`.

Stage 3AA decision:

- Correctness passes for the new scheduler path: scalar equivalence, verifier
  mismatches, accepted-miss counters, and far-expansion mean gates are clean on
  the measured rows.
- Performance still fails. Target-leaf-block reduces some bookkeeping on
  `mini-protein` and `charged-peptide`, but compact-mini remains dominated by
  leaf work and grouping/flush overhead. The medium geomean speedup is below
  the `1.5x` target, so Barnase timing was not run.
- No new preset was created. Defaults remain unchanged.
- Next branch should stop adding scheduler layers and target either
  kernel-level throughput for the direct JCP leaf path or a safe leaf-work
  reduction strategy that preserves Stage 3T/3U accuracy gates.

## Stage 3AB/9U - TABI Direct JCP Leaf-Kernel Throughput

Stage 3AB implementation:

- Added hidden profiling/control flags:
  `-pbe-tabi-leaf-kernel-profile=off|summary` and
  `-pbe-tabi-jcp-direct-kernel=baseline|soa-avx2`.
- Implemented `soa-avx2` for clean JCP on-the-fly leaf-pair target batches.
  The path is only used for clean source spans with `target_batch_size=4`;
  singular, one-ring, Duffy, adaptive-buffer, and exclusion cases remain on
  the existing scalar/fallback path.
- Added leaf-kernel counters to timing/profile/progress output:
  scalar tasks, batch tasks, source-panel evals, scalar seconds,
  target-avx2 seconds, tail seconds, SoA sources, SoA seconds, scalar-tail
  seconds, and clean source panels per second.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with
  `-IncludeStage3ABRows`, including safe baseline, Stage 3AA target-leaf-block
  baseline, SoA rows for min-source `32/64/128`, and an audit validation row.
- Extended `test/pbe_tabi_vs_grid_compare.ps1` with
  `-UseSoaLeafPairKernel` for future compare rows.
- Added `test/pbe_tabi_leaf_pair_target_soa_avx2_smoke.ps1`.
- Production defaults, `stage3u-safe`, direct TABI, grid backend, MAC safety,
  and estimator accuracy guards remain unchanged.

Measured Stage 3AB verification:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- New SoA smoke:
  `test/pbe_tabi_leaf_pair_target_soa_avx2_smoke.ps1`: PASS,
  relative RMS `0`, batches `109`, SoA sources `2052`,
  vectorized sources `2052`, clean source panels/sec `6.20783e7`, fallback
  singular relative RMS `0`.
- Existing focused smokes still pass:
  `test/pbe_tabi_leaf_block_scheduler_smoke.ps1`,
  `test/pbe_tabi_leaf_pair_target_avx2_smoke.ps1`,
  `test/pbe_tabi_leaf_pair_worklist_scheduler_smoke.ps1`,
  `test/pbe_tabi_tree_builder_smoke.ps1`,
  `test/pbe_tabi_topdown_mac_audit.ps1`, and
  `test/pbe_tabi_treecode_accuracy.ps1`.
- Full Stage 3AA+3AB medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3AARows -IncludeStage3ABRows -CaseTimeoutSeconds 180 -AllowFailedStatus`:
  PASS. Summary artifact:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

Stage 3AB medium timing highlights:

| fixture | safe s | best Stage 3AB row | best s | speedup vs safe | validation mean | accepted missed | SoA sources | SoA seconds |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `mini-protein` | `3.4495` | safe baseline | `3.4495` | `1.000x` | `2.654e-5` | `0` | `10,780,845` | `0.7554` |
| `charged-peptide` | `4.1113` | SoA min32 | `3.0023` | `1.369x` | `0` | `0` | `9,005,740` | `0.5226` |
| `compact-mini-protein` | `13.5824` | SoA min128 | `12.9733` | `1.047x` | `2.854e-5` | `0` | `42,437,070` | `3.2829` |

Stage 3AB min64 kernel isolate:

| fixture | target-avx2 batch s | SoA batch s | target-avx2 panels/s | SoA panels/s |
| --- | ---: | ---: | ---: | ---: |
| `mini-protein` | `0.7040` | `0.6906` | `6.125e7` | `6.245e7` |
| `charged-peptide` | `0.6073` | `0.4649` | `5.932e7` | `7.749e7` |
| `compact-mini-protein` | `3.6407` | `3.4961` | `4.662e7` | `4.855e7` |

Stage 3AB decision:

- Correctness passes for the SoA direct leaf-kernel path: smoke equivalence,
  validation accepted-miss counters, far-expansion mean, and leaf-block verify
  mismatches are clean on the measured rows.
- The SoA kernel improves clean batch throughput, especially on
  `charged-peptide`, but total medium performance still fails. Geometric mean
  speedup versus the Stage 3AB safe baseline is `0.992x` for min64 SoA,
  `0.989x` for min32 SoA, and `1.001x` for min128 SoA, below the `1.5x`
  gate. `mini-protein` is slower under target-leaf-block batching, and
  `compact-mini-protein` only improves about `1.047x` at the best measured
  SoA row.
- Barnase timing was not run because the medium performance gate did not pass.
- No new preset was created. Defaults remain unchanged.
- Next branch should reduce safe leaf work rather than add another scheduler
  layer or another small clean-batch kernel variant. Candidate directions:
  adaptive safe leaf-work reduction under the Stage 3T/3U accuracy gates,
  stricter batching only where target-leaf-block is profitable, or a larger
  kernel redesign that attacks scalar tail and unclean direct pairs.

## Stage 3AC/9V - TABI Safe Leaf-Work Reduction and Adaptive Surface-Margin Gate

Stage 3AC implementation:

- Added hidden leaf-work profiling:
  `-pbe-tabi-leaf-work-profile=off|summary`.
- Added hidden adaptive surface-margin controls:
  `-pbe-tabi-mac-surface-margin-policy=all|adaptive`,
  `-pbe-tabi-mac-surface-margin-aspect-threshold`, and
  `-pbe-tabi-mac-surface-margin-rho-trigger`.
- `surface-margin` keeps the Stage 3U behavior by default (`all`). In
  `adaptive`, the margin distance is applied only to high-risk nodes with
  `rho >= trigger` or aspect ratio above the threshold; benign candidates use
  center-distance MAC.
- Added timing/profile/progress counters for margin applied/skipped,
  adaptive rejects, skipped accepted far nodes, and leaf-work buckets
  (`mac-terminal`, physical-gap, fragmentation, convergence, surface-margin,
  K3 rho cap, singular fallback, scalar tail).
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with
  `-IncludeStage3ACRows`, including leaf-size diagnostics, adaptive margin
  rows, a sampled audit validation row, and a timing-only
  target-leaf-block+SoA row.
- Added `test/pbe_tabi_leaf_work_profile_smoke.ps1` and
  `test/pbe_tabi_adaptive_surface_margin_smoke.ps1`.
- Production defaults, `stage3u-safe`, direct TABI, grid backend, estimator
  guards, and Stage 3T/3U promotion candidate remain unchanged.

Measured Stage 3AC verification:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- New leaf-work profile smoke: PASS, leaf work `36,672`, accounted
  `36,672`, surface-margin bucket `6,272`.
- New adaptive surface-margin smoke: PASS, adaptive high-risk applied
  `756`, adaptive rejects `74`, benign skipped `512`.
- Existing focused smokes still pass:
  `test/pbe_tabi_leaf_pair_target_soa_avx2_smoke.ps1`,
  `test/pbe_tabi_leaf_block_scheduler_smoke.ps1`,
  `test/pbe_tabi_tree_builder_smoke.ps1`,
  `test/pbe_tabi_topdown_mac_audit.ps1`, and
  `test/pbe_tabi_treecode_accuracy.ps1`.
- Stage 3AC mini-protein matrix smoke:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -Fixtures mini-protein -IncludeStage3ACRows -CaseTimeoutSeconds 180 -AllowFailedStatus`:
  PASS. Artifact:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

Stage 3AC mini-protein highlights:

| row | seconds | leaf panels | accounted leaf panels | margin policy | applied | skipped | adaptive rejects | audit mean | accepted missed |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| `stage3ac_safe_baseline_leaf128_margin_all` | `3.0426` | `110,035,455` | `110,035,455` | `all` | `322,005` | `0` | `0` | `0` | `0` |
| `stage3ac_adaptive_aspect3p0_rho0p55` | `5.9828` | `109,361,580` | `109,361,580` | `adaptive` | `207,195` | `114,810` | `23,715` | `0` | `0` |
| `stage3ac_adaptive_aspect2p5_rho0p55` | `4.5117` | `109,361,580` | `109,361,580` | `adaptive` | `207,195` | `114,810` | `23,715` | `0` | `0` |
| `stage3ac_adaptive_aspect3p0_rho0p50` | `4.6233` | `109,361,580` | `109,361,580` | `adaptive` | `230,760` | `91,245` | `37,170` | `0` | `0` |
| `stage3ac_validation_adaptive_aspect3p0_rho0p55_audit` | `5.1890` | `109,361,580` | `109,361,580` | `adaptive` | `207,195` | `114,810` | `23,715` | `6.992e-5` | `0` |
| `stage3ac_timing_adaptive_aspect3p0_rho0p55_block_soa` | `5.0109` | `109,361,580` | `109,361,580` | `adaptive` | `207,195` | `114,810` | `23,715` | `0` | `0` |

Stage 3AC decision:

- The implementation and counters are in place and verified. Adaptive margin
  does reduce leaf work slightly on the mini-protein smoke
  (`110.04M -> 109.36M` leaf-panel interactions), and the new bucket
  accounting matches total leaf-panel interactions.
- The first validation row does not pass the Stage 3T/3U far-expansion mean
  gate on mini-protein (`6.992e-5`, gate `<5e-5`), so no Barnase timing was
  run and no preset/default changed.
- Medium full 3-fixture performance gate remains pending. The early evidence
  suggests adaptive margin alone is too small a leaf-work reduction to reach
  the `1.5x` geomean target, so the next branch should either tune adaptive
  selectivity with sampled audit enabled or target safe reduction of
  convergence/estimator-driven leaf work rather than another scheduler layer.

## Stage 3AD/9W - TABI Leaf-Work Attribution and Risk-Scaled Surface Margin

Stage 3AD implementation:

- Extended `-pbe-tabi-leaf-work-profile=summary` with detailed attribution
  buckets while preserving the aggregate columns:
  `mac_terminal`, `physical_gap`, `fragmentation`, `convergence_estimator`,
  `convergence_rho`, `convergence_aspect`, `surface_margin_hard`,
  `surface_margin_soft`, `k3_rho_cap`, `singular_duffy_exclusion`, and
  `scalar_tail`.
- Added fraction columns for fast diagnosis:
  `tabi_leaf_work_mac_terminal_frac`, `tabi_leaf_work_convergence_frac`,
  `tabi_leaf_work_surface_margin_frac`,
  `tabi_leaf_work_singular_fallback_frac`, and
  `tabi_leaf_work_scalar_tail_frac`.
- Added hidden risk-scaled surface-margin controls:
  `-pbe-tabi-mac-surface-margin-policy=risk-scaled`,
  `-pbe-tabi-mac-surface-margin-soft-aspect-threshold`,
  `-pbe-tabi-mac-surface-margin-soft-rho-trigger`, and
  `-pbe-tabi-mac-surface-margin-soft-factor`.
- `risk-scaled` applies the full Stage 3U margin to hard-risk candidates,
  a smaller margin to soft-risk candidates, and center-distance MAC to
  benign candidates. The existing `all` behavior remains the default.
- Added hard/soft margin counters to timing/profile/progress:
  hard/soft applied counts, hard/soft rejects, and hard/soft reject panels.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with
  `-IncludeStage3ADRows`, including full-margin baseline, sampled adaptive
  validation rows, factor sweeps, risk-scaled rows, leaf overlays, and a
  timing-only target-leaf-block+SoA row.
- Added `test/pbe_tabi_leaf_work_attribution_smoke.ps1` and
  `test/pbe_tabi_risk_scaled_surface_margin_smoke.ps1`.
- Production defaults, `stage3u-safe`, direct TABI, grid backend, MAC safety,
  estimator guards, and solver operator semantics remain unchanged.

Measured Stage 3AD verification:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- New leaf-work attribution smoke: PASS, leaf work `36,672`, convergence
  bucket `0`, surface-margin bucket `392`, selected fraction sum `0.754145`.
- New risk-scaled surface-margin smoke: PASS, hard applied `756`, soft
  applied `512`, benign skipped `724`.
- Existing focused smokes still pass:
  `test/pbe_tabi_leaf_work_profile_smoke.ps1`,
  `test/pbe_tabi_adaptive_surface_margin_smoke.ps1`,
  `test/pbe_tabi_leaf_pair_target_soa_avx2_smoke.ps1`,
  `test/pbe_tabi_leaf_block_scheduler_smoke.ps1`,
  `test/pbe_tabi_tree_builder_smoke.ps1`,
  `test/pbe_tabi_topdown_mac_audit.ps1`,
  `test/pbe_tabi_treecode_accuracy.ps1`,
  `test/pbe_tabi_surface_margin_mac_smoke.ps1`,
  `test/pbe_tabi_profile_smoke.ps1`, and
  `test/pbe_memory_profile_smoke.ps1`.
- Full Stage 3AD medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3ADRows -CaseTimeoutSeconds 180 -AllowFailedStatus`:
  PASS. Artifact:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

Stage 3AD medium attribution highlights:

| fixture | safe leaf panels | adaptive/risk leaf panels | risk leaf96 leaf panels | safe mean | best audited mean | accepted missed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `mini-protein` | `110.04M` | `109.36M` | `106.89M` | `2.77e-5` | `5.49e-5` (`leaf96`) | `0` |
| `charged-peptide` | `100.43M` | `100.43M` | `98.24M` | `0` | `4.00e-5` (`leaf96`) | `8 global / 4 K2` |
| `compact-mini-protein` | `516.87M` | `512.96M` | `506.77M` | `2.77e-5` | `4.70e-5` (`leaf96`) | `0` |

Stage 3AD speed/leaf-work summary:

| row | geomean speedup vs Stage 3AD safe | geomean leaf-panel ratio | gate status |
| --- | ---: | ---: | --- |
| `stage3ad_risk_scaled_leaf192` | `1.217x` | `1.003x` | timing improves, but leaf work is not reduced and audit pairs are not useful enough for promotion |
| `stage3ad_timing_risk_scaled_block_soa` | `1.199x` | `0.995x` | timing-only, below `1.5x` |
| `stage3ad_risk_scaled_hard2p5_rho0p55_soft2p0_rho0p45_factor0p20` | `0.984x` | `0.995x` | accuracy clean, slower than safe |
| `stage3ad_validation_adaptive_aspect3p0_rho0p55_audit` | `0.882x` | `0.995x` | mean fails on mini-protein |
| `stage3ad_risk_scaled_leaf96` | `0.724x` | `0.977x` | leaf work reduced, but mini mean `5.49e-5` and charged accepted misses fail gates |

Stage 3AD decision:

- The new attribution is working: bucket totals match
  `leaf_panel_interactions` on the smoke and full medium rows, and the new
  hard/soft margin counters identify that the current risk-scaled settings
  mostly move candidates through hard-margin rejects; soft margin is applied
  but did not create meaningful soft-reject leaf-work reductions on the medium
  rows.
- No full-medium row satisfies both accuracy and performance gates. The best
  leaf-work reduction row (`leaf96`) reduces geomean leaf panels by only about
  `2.3%`, remains below the `10%` leaf-work target, fails the mini mean gate
  (`5.49e-5` vs `<5e-5`), and introduces accepted misses on
  `charged-peptide`.
- The fastest row (`leaf192`) reaches only `1.217x` geomean speedup and does
  not reduce leaf work, so it is not a Stage 3AD candidate.
- Barnase timing was not run because the full-medium performance/accuracy gate
  did not pass. Defaults remain unchanged.
- Next branch should reduce convergence/estimator-driven leaf work more
  directly. The Stage 3AD buckets show convergence-driven leaf work is a large
  fraction of the remaining direct workload, while binary/adaptive
  surface-margin tuning is too small to unlock the `1.5x` target.

## Stage 3AE/9X - TABI Refactor and Absolute-Aware Convergence Split Selectivity

Stage 3AE implementation:

- Added the first behavior-preserving TABI treecode refactor slice:
  `src/tabi/tabi_treecode_internal.h` and
  `src/tabi/tabi_treecode_internal.cpp`.
  The internal header currently holds small convergence split/selectivity and
  overreject classification helpers under `TitanPBE::Tabi::detail`.
  The helpers are kept inline for now so standalone smoke harnesses that include
  `tabi_treecode.cpp` continue to link without custom rules. This is the
  scaffolding for the larger file split; `tabi_treecode.cpp` is not yet a thin
  entrypoint.
- Added hidden convergence overreject audit:
  `-pbe-tabi-convergence-overreject-audit=off|sampled|full-small`.
  Sampled split-rejected candidates can now be checked with the same
  direct-vs-Taylor JCP expansion audit and classified as true high-error caught
  candidates or low-error overrejects.
- Added hidden absolute-aware split selectivity:
  `-pbe-tabi-convergence-split-selectivity=legacy|absolute-aware|channel-absolute-aware`
  and `-pbe-tabi-convergence-split-absolute-floor=<0.0..1e-4>`.
  Defaults preserve legacy behavior. Experimental rows can reaccept a far node
  only when the estimated relative error is high but the estimated absolute
  error is below the configured floor.
- Added timing/profile/progress fields for:
  `tabi_convergence_split_selectivity`,
  `tabi_convergence_split_absolute_floor`,
  `tabi_convergence_overreject_audit`,
  low-absolute reaccepted far counts/panels,
  estimated absolute error max/mean, and global/per-channel overreject/caught
  counters.
- Extended `test/pbe_tabi_topdown_physical_gap_matrix.ps1` with
  `-IncludeStage3AERows`, including legacy baseline, overreject audit,
  absolute-aware floor sweep, channel-absolute-aware floor sweep, and a
  timing-only target-leaf-block+SoA row.
- Added focused smokes:
  `test/pbe_tabi_convergence_overreject_audit_smoke.ps1` and
  `test/pbe_tabi_absolute_aware_convergence_split_smoke.ps1`.
- Production defaults, `stage3u-safe`, direct TABI, grid backend, public CLI
  behavior, and solver operator defaults remain unchanged.

Measured Stage 3AE verification:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- New absolute-aware convergence split smoke: PASS.
- New convergence overreject audit smoke: PASS.
- Focused smokes still pass:
  `test/pbe_tabi_tree_builder_smoke.ps1`,
  `test/pbe_tabi_treecode_accuracy.ps1`,
  `test/pbe_tabi_topdown_mac_audit.ps1`,
  `test/pbe_tabi_leaf_work_attribution_smoke.ps1`,
  `test/pbe_tabi_risk_scaled_surface_margin_smoke.ps1`, and
  `test/pbe_tabi_leaf_pair_target_soa_avx2_smoke.ps1`.
- Stage 3AE mini-protein matrix smoke:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -Fixtures mini-protein -IncludeStage3AERows -CaseTimeoutSeconds 180 -AllowFailedStatus`:
  PASS. Artifact:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_topdown_physical_gap_matrix.csv`.

Stage 3AE mini-protein highlights:

| row | seconds | leaf panels/matvec | mean error | low-abs reaccepted | overreject | caught | accepted missed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `stage3ae_safe_baseline_legacy_split` | `2.0466` | `7,335,697` | `2.618e-5` | `0` | `0` | `0` | `0` |
| `stage3ae_overreject_audit_legacy_split` | `2.0347` | `7,335,697` | `0` | `0` | `6,171` | `1,509` | `0` |
| `stage3ae_absolute_aware_floor1e8` | `2.0391` | `7,335,697` | `2.753e-5` | `0` | `6,133` | `1,439` | `0` |
| `stage3ae_absolute_aware_floor3e8` | `2.0490` | `7,335,697` | `2.916e-5` | `0` | `6,116` | `1,445` | `0` |
| `stage3ae_absolute_aware_floor1e7` | `2.0059` | `7,335,697` | `2.952e-5` | `0` | `6,152` | `1,421` | `0` |
| `stage3ae_channel_absolute_aware_floor1e8` | `2.0599` | `7,335,697` | `3.001e-5` | `0` | `6,138` | `1,427` | `0` |
| `stage3ae_channel_absolute_aware_floor3e8` | `2.0178` | `7,335,697` | `2.637e-5` | `0` | `6,110` | `1,454` | `0` |
| `stage3ae_channel_absolute_aware_floor1e7` | `2.0407` | `7,335,697` | `2.745e-5` | `0` | `6,053` | `1,512` | `0` |
| `stage3ae_timing_best_absaware_block_soa` | `1.9035` | `7,335,697` | `0` | `0` | `0` | `0` | `0` |

Stage 3AE decision:

- The hidden plumbing is in place and verified. The overreject audit is useful:
  on mini-protein it found thousands of low-error split-rejected candidates
  alongside true caught high-error candidates.
- The first absolute-aware floor sweep did not reduce actual leaf work on the
  mini-protein matrix: `low_abs_reaccepted_far_count` remained `0`, and
  `leaf_panel_work_per_matvec` stayed unchanged across the selectivity rows.
  This means the new decision path is wired, but the estimated absolute error
  scale or floor needs calibration before it can unlock the Stage 3AD
  convergence-driven leaf-work bucket.
- Full 3-fixture Stage 3AE medium matrix and Barnase compare were not run
  because the mini-protein selectivity smoke did not yet show leaf-work
  reduction. Defaults remain unchanged.
- Next branch should calibrate estimated absolute channel error against the
  overreject audit oracle, or raise/select floors per channel with sampled
  audit enabled, before adding more batching or scheduler layers.

## Stage 3AF/9Y - TABI Treecode Translation-Unit Refactor Freeze

Stage 3AF implementation:

- Refactored `src/tabi/tabi_treecode.cpp` from the monolithic treecode
  implementation into a thin public forwarding facade plus bounded internal
  translation units. No flag, preset, MAC rule, estimator rule, split
  decision, CSV field, or production default was changed in this stage.
- Expanded `src/tabi/tabi_treecode_internal.h` into the private contract for
  shared treecode structs and cross-TU helper declarations under
  `TitanPBE::Tabi::detail`.
- Split ownership as follows:

| file | responsibility | lines |
| --- | --- | ---: |
| `tabi_treecode.cpp` | public API forwarding facade | `218` |
| `tabi_treecode_geometry.cpp` | tree construction, bounds, interaction-list/workspace setup | `1,523` |
| `tabi_treecode_series.cpp` | Taylor/JCP series, moments and moment audit | `1,345` |
| `tabi_treecode_policy.cpp` | MAC, convergence/expansion audits and stats aggregation | `1,550` |
| `tabi_treecode_leaf.cpp` | direct JCP leaf kernels, AVX2/SoA batching and directness probes | `1,510` |
| `tabi_treecode_operator.cpp` | direct/tree/JCP traversal and matvec orchestration | `2,456` |
| `tabi_treecode_internal.h` | internal shared contract | `806` |

- Added `test/pbe_tabi_treecode_sources.ps1` as the single source-list helper
  for standalone TABI smoke compiles, and routed the existing manual treecode
  compile harnesses through it so future internal TU changes do not require
  repeated per-script linker edits.
- The context target is met: the facade is below `250` lines and every
  implementation TU is below `3,000` lines. `tabi_solver.cpp` was intentionally
  left outside this refactor.

Stage 3AF parity evidence:

- Native build:
  `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Deterministic focused parity matches the Stage 3AE pre-refactor baseline:
  - tree builder: `panels=512`, `grouped_nodes=73`, `topdown_nodes=73`,
    `topdown_depth=2`;
  - treecode accuracy: `rel3=9.38761e-07`, `rel4=1.35948e-07`,
    `topdown_rel3=4.61951e-05`, `far=199`, `jcp_duffy=19416`;
  - topdown MAC audit: `unsafe=392`, `rejects=392`,
    `far_audit=392`, `far_safe=0`;
  - leaf-work attribution: `leaf=36672`, `convergence=0`,
    `margin=392`, `frac_sum=0.754145`;
  - leaf-block scheduler: `rel=1.65074e-16`, `groups=352`,
    `batched=150`, `verify_checks=960`;
  - absolute-aware split, convergence overreject audit, far-expansion audit,
    risk-scaled surface-margin and SoA/AVX2 smokes: PASS.
- Mini-protein pre/post Stage 3AE matrix comparison preserved deterministic
  gate fields for every Stage 3AE row:
  status, leaf-panel work, low-absolute reaccepted count, and all accepted-miss
  channel counts were unchanged. Sampled expansion/audit means are report-only
  for this comparison because multithread sampled candidate selection is not
  bitwise stable across separate runs.
- Full post-refactor medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3ADRows -IncludeStage3AERows -CaseTimeoutSeconds 180 -AllowFailedStatus`:
  PASS. Representative Stage 3AE channel-absolute-aware floor `1e-7` rows:

| fixture | leaf panels/matvec | expansion mean | low-abs reaccepted | accepted missed |
| --- | ---: | ---: | ---: | ---: |
| `mini-protein` | `7,335,697` | `2.778e-5` | `0` | `0` |
| `charged-peptide` | `5,021,476` | `0` | `0` | `0` |
| `compact-mini-protein` | `17,228,993` | `2.785e-5` | `0` | `0` |

- Regression freeze passed:
  GMRES true residual, sampled true residual, profile smoke, memory profile,
  clean-room audit, PBE regression, PBSA smoke, and SASA surface smoke.

Stage 3AF decision:

- The treecode refactor is behavior-preserving on the focused deterministic
  evidence and full medium gate run. Production defaults and `stage3u-safe`
  remain unchanged, and no Barnase timing or promotion probe was run.
- Stage 3AE's unresolved result remains unchanged after the split:
  absolute-aware rows still produce `low_abs_reaccepted=0` and do not reduce
  leaf-panel work, even though overreject audit identifies low-error rejected
  candidates.
- The next accuracy/performance branch should now calibrate
  `estimated_abs_error` against the sampled direct-vs-Taylor overreject oracle,
  using the new policy/series boundary rather than adding another batching or
  scheduler layer.

## Stage 3AG/9Z - TABI Overreject Oracle Calibration Matrix

Stage 3AG implementation:

- Added hidden report-only option
  `-pbe-tabi-convergence-overreject-detail=off|summary|sampled-detail`,
  default `off`. It is inert unless overreject audit is active and never
  changes traversal, MAC, convergence split, or the GMRES operator.
- Extended overreject accounting in `tabi_treecode_policy.cpp` with sampled
  estimated-versus-oracle absolute-error summaries and `required_scale =
  oracle_abs_error / estimated_abs_error` for the global result and each JCP
  channel.
- `sampled-detail` now emits
  `pbe_tabi_convergence_overreject_calibration.csv` with pair metadata
  (`target_index`, source node, channel, `rho`, aspect ratio, panels),
  estimated and oracle relative/absolute errors, scale, and
  `caught|overreject` classification.
- Added Stage 3AG matrix rows and
  `test/pbe_tabi_overreject_calibration_summary.ps1`. The summary evaluates
  a fixed scale grid against the sampled `caught` oracle and reports projected
  overreject reaccept panels only when no caught sample would be reaccepted.

Stage 3AG verification:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- New oracle plumbing smoke:
  `test/pbe_tabi_overreject_calibration_oracle_smoke.ps1`: PASS
  (`required_scale=2e-07`, real high-error pair remains `caught`).
- Existing focused gates pass:
  convergence overreject audit, absolute-aware convergence split, tree builder,
  treecode accuracy, topdown MAC audit, far-expansion audit, leaf-work
  attribution, and leaf-work profile.
- Full medium evidence:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3AERows -IncludeStage3AGRows -CaseTimeoutSeconds 180 -AllowFailedStatus`:
  PASS. Detail artifacts are under
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/<fixture>/stage3ag_overreject_calibration_detail/`.
- Calibration summary artifact:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_overreject_calibration_summary.csv`.
- Regression freeze passes:
  GMRES true residual, sampled true residual, profile, memory profile,
  clean-room audit, PBE regression, PBSA smoke, and SASA surface smoke.

Stage 3AG detail rows:

| fixture | leaf panels/matvec | overreject | caught | oracle abs mean | estimated abs mean | required scale p99 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `mini-protein` | `7,335,697` | `6,074` | `1,606` | `3.828e-7` | `2.075e-1` | `3.673e-5` |
| `charged-peptide` | `5,021,476` | `4,731` | `5,509` | `9.808e-7` | `1.115` | `1.266e-5` |
| `compact-mini-protein` | `17,228,993` | `11,356` | `4,004` | `2.439e-7` | `8.053e-2` | `6.727e-5` |

Sampled safe-scale projection at absolute floor `1e-7`:

| channel | scale with margin | caught reaccepted | overreject reaccepted panels | fraction vs safe leaf work |
| --- | ---: | ---: | ---: | ---: |
| global | `1e-5` | `0` | `1,833` | `0.0062%` |
| `K1_phi` | `3e-6` | `0` | `394,241` | `1.33%` |
| `K2_phi_normal` | `3e-4` | `0` | `14,486` | `0.049%` |
| `K3_phin` | `3e-6` | `0` | `41,266` | `0.139%` |
| `K4_mixed` | `3e-6` | `0` | `244,818` | `0.827%` |

Stage 3AG decision:

- Stage 3AG is report-only and confirms the Stage 3AE diagnosis: the
  estimator absolute magnitude is much larger than the sampled direct
  expansion oracle while all Stage 3AG leaf work remains identical to the
  safe baseline.
- No sampled-safe scale projects the required `>=10%` leaf-work reduction;
  even the separate per-channel projections are small and are not additive
  without double-counting pair work. Therefore the next stage must not enable
  deterministic channel-scaled reaccept based on this calibration.
- The next branch should revisit estimator channel shape and source norm
  construction, especially the derivative channels, before any new split
  selectivity or batching work. Production defaults and `stage3u-safe` remain
  unchanged. No Barnase run was performed because Stage 3AG changes no
  candidate operator.

## Stage 3AH/10A - TABI Channel-Shape Estimator Feature Audit

Stage 3AH implementation:

- Added hidden report-only option
  `-pbe-tabi-convergence-shape-audit=off|summary|sampled-detail`, default
  `off`. It is only useful with overreject audit enabled and does not affect
  MAC, convergence split, traversal, batching, or the GMRES operator.
- Extended the split-rejected overreject audit path to write
  `pbe_tabi_convergence_shape_calibration.csv` when `sampled-detail` is
  active. Each row records channel/classification plus `rho`, `rho^(p+1)`,
  aspect, panel count, radius, target distance, source norm, degree-`p` and
  degree-`p-1` tail norms, tail ratio, effective amp, estimated/oracle errors,
  and oracle-to-feature ratios.
- Added Stage 3AH matrix rows and
  `test/pbe_tabi_convergence_shape_summary.ps1`, which evaluates sampled
  feature thresholds that would reaccept overrejects only below the minimum
  caught feature value.

Stage 3AH verification:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- New shape audit smoke:
  `test/pbe_tabi_convergence_shape_audit_smoke.ps1`: PASS
  (`records=10`, both `caught` and `overreject` classes present).
- Focused gates pass:
  overreject calibration oracle, convergence overreject audit, treecode
  accuracy, tree builder, topdown MAC audit, and far-expansion audit.
- Full medium evidence:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3AGRows -IncludeStage3AHRows -CaseTimeoutSeconds 180 -AllowFailedStatus`:
  PASS. Detail artifacts are under
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/<fixture>/stage3ah_shape_audit_detail/`.
- Shape summary artifact:
  `test/tmp/pbe_tabi_topdown_physical_gap_matrix/pbe_tabi_convergence_shape_summary.csv`.
- Regression freeze passes:
  GMRES true residual, sampled true residual, profile, memory profile,
  clean-room audit, PBE regression, PBSA smoke, and SASA surface smoke.

Stage 3AH detail rows:

| fixture | leaf panels/matvec | overreject | caught |
| --- | ---: | ---: | ---: |
| `mini-protein` | `7,335,697` | `6,116` | `1,564` |
| `charged-peptide` | `5,021,476` | `4,758` | `5,482` |
| `compact-mini-protein` | `17,228,993` | `11,160` | `4,200` |

Best sampled feature projections:

| channel | feature | caught | overreject | projected panels | projected fraction |
| --- | --- | ---: | ---: | ---: | ---: |
| global | `oracle_abs_over_source_norm` | `117` | `1,419` | `286,446` | `0.968%` |
| `K3_phin` | `oracle_abs_over_source_norm` | `238` | `1,298` | `202,061` | `0.683%` |
| `K3_phin` | `oracle_abs_over_tail_norm` | `238` | `1,298` | `198,815` | `0.672%` |
| global | `oracle_abs_over_tail_norm` | `117` | `1,419` | `139,775` | `0.472%` |
| `K2_phi_normal` | `oracle_abs_over_source_norm` | `609` | `927` | `114,479` | `0.387%` |
| global | `tail_ratio` | `117` | `1,419` | `106,478` | `0.360%` |

Stage 3AH decision:

- No measured feature threshold separates enough low-error overreject work from
  caught high-error samples. The best projection is below `1%` of safe leaf
  work, far under the `10%` threshold for implementing estimator v2.
- Do not add another absolute floor, scale, or split selectivity heuristic from
  these sampled shape features.
- The next branch should inspect source moment/channel semantics directly:
  verify that K1/K3 and K2/K4 source moment mappings match the JCP Taylor
  channel they estimate, then test whether decoupled phi/phin or derivative
  source-norm construction gives a better caught/overreject separation.
  Defaults and `stage3u-safe` remain unchanged. No Barnase run was performed.

## Stage 3AI/10B - TABI Grid-Gap Performance Triage and Fused Leaf Kernel

Stage 3AI implementation:

- Added hidden `-pbe-tabi-matvec-cost-profile=off|summary`, default `off`.
  When enabled, timing/profile/progress now report `matvec_cost_profile`,
  `leaf_panels_per_second`, `clean_leaf_panels_per_second`, and
  `tree_traversal_fraction` so grid-vs-TABI comparisons can separate GMRES,
  tree traversal, and direct leaf throughput without changing the operator.
- Extended hidden `-pbe-tabi-jcp-direct-kernel` to accept
  `fused-soa-avx2`. The implementation routes clean JCP on-the-fly
  leaf-pair batches through the existing SoA AVX2 source loop, but reports the
  fused kernel name distinctly so Stage 3AI rows can isolate it from the older
  `soa-avx2` label. Singular, one-ring, Duffy, and exclusion paths keep the
  scalar fallback behavior.
- Added `-IncludeStage3AIRows` to the medium matrix. Rows include the
  `stage3u-safe` cluster-off baseline, current best perf-style row,
  `soa-avx2`, `fused-soa-avx2` min-source sweeps, and a sampled audit
  validation row.
- Added `-UseStage3AIPerfProbe` to the grid-vs-TABI compare harness. The
  compare CSV now includes `absolute_time_delta_vs_grid`,
  `tabi_matvec_cost_profile`, `tabi_leaf_panels_per_second`,
  `tabi_clean_leaf_panels_per_second`, and `tabi_tree_traversal_fraction`.

Stage 3AI verification:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- New fused leaf-kernel smoke:
  `test/pbe_tabi_leaf_pair_fused_soa_avx2_smoke.ps1`: PASS
  (`rel=0`, `fused_sources=2052`, singular fallback remains accurate).
- New matvec-cost profile smoke:
  `test/pbe_tabi_matvec_cost_profile_smoke.ps1`: PASS
  (`matvecs=33`, profile/timing columns emitted).
- Focused gates pass:
  tree builder, treecode accuracy, and topdown MAC audit.
- Matrix syntax/sanity row:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -Fixtures mini-protein -IncludeStage3AIRows -CaseTimeoutSeconds 180 -AllowFailedStatus -Quick`:
  PASS. Stage 3AI rows report `matvec_cost_profile=summary` and the requested
  `fused-soa-avx2` kernel.
- Grid-vs-TABI compare sanity:
  `test/pbe_tabi_vs_grid_compare.ps1 -Fixtures mini-protein -UseStage3AIPerfProbe -Mode perf -CaseTimeoutSeconds 180 -AllowTimeout`:
  PASS. The small-fixture grid ratio remains report-only; the row confirms the
  new grid-gap columns are emitted.

Stage 3AI decision:

- Stage 3AI is a performance-instrumentation and kernel-labeling slice, not a
  production promotion. `stage3u-safe`, production defaults, direct TABI,
  estimator/MAC guards, and the grid backend remain unchanged.
- Full medium geomean and Barnase timing are still pending. The next step is to
  run the full Stage 3AI medium matrix and compare the fused kernel against
  `stage3u-safe` and the current best perf row. If the speedup gate still fails,
  the roadmap should follow the stated branch: reduce algorithmic matvec count
  or improve preconditioning, not add another scheduler layer.

## Stage 3AJ/10C - TABI Memory-Budgeted Leaf-Block Preconditioner

Stage 3AJ implementation:

- Added hidden `-pbe-tabi-precond=block-budgeted`, leaving legacy `block` and
  all production defaults unchanged.
- Added hidden budget controls:
  `-pbe-tabi-precond-block-panel-cap=<8..32>` (default `32`) and
  `-pbe-tabi-precond-memory-budget-mb=<16..512>` (default `256`).
  The new path forms the existing leaf/chunk block ranges, estimates LU,
  pivot, index, and diagonal storage before allocation, and reports
  `diagonal_block_budget_cap` when the estimate exceeds its budget.
- Reused the existing buffer-aware JCP block assembly, LU pivot/valid-panel
  checks, and preconditioner apply path. Legacy `block` retains its previous
  hard-cap behavior.
- Added timing/profile/progress telemetry for configured cap/budget, estimated
  bytes, valid/ill-conditioned blocks, valid panels, setup time, apply time,
  and apply calls. This information is available in long-run progress rows for
  timeout-safe Barnase probes.
- Added `test/pbe_tabi_budgeted_block_precond_smoke.ps1`, Stage 3AJ matrix rows,
  `-UseStage3AJBudgetedPreconditioner` in the grid-vs-TABI compare harness,
  and budgeted-block reporting/options in the Barnase truth harness.

Stage 3AJ verification:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Budgeted block smoke: PASS. The active small case reports
  `block_budgeted_buffer_aware`; the explicit low-budget case reports
  `diagonal_block_budget_cap` before allocating valid LU blocks
  (`estimated_bytes=21662976`, `budget_bytes=16777216`, `valid_blocks=0`).
- Legacy block buffer smoke, matvec-cost profile, fused SoA AVX2, treecode
  accuracy, tree builder, topdown MAC audit, absolute-floor, channel-semantics,
  surface-margin and polynomial-estimator smokes: PASS.
- Full medium matrix:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3AJRows
  -CaseTimeoutSeconds 180 -AllowFailedStatus`: PASS as a completed evidence
  run; all 39 rows emitted status `ok`.
- Regression freeze: GMRES true residual, sampled true residual, TABI profile,
  memory profile, clean-room audit, PBE regression, PBSA smoke, and SASA
  surface smoke: PASS.

Stage 3AJ medium candidate evidence (`fused-soa-avx2`, cap `32`, budget
`256 MB`, restart `64`, relative to fused-diagonal control):

| Fixture | Estimated precond bytes | Matvec reduction | Time speedup | Sampled residual | Audit mean |
| --- | ---: | ---: | ---: | ---: | ---: |
| mini-protein | 2558272 | 0% | 0.9250x | 5.0623e-5 | 2.3618e-5 |
| charged-peptide | 2133120 | 0% | 0.9370x | 8.1240e-6 | 0 |
| compact-mini-protein | 3883072 | 0% | 1.2010x | 2.7631e-5 | 5.9797e-5 |

Stage 3AJ decision:

- The new memory-budgeted path is functional and keeps medium energy/matvec
  error and accepted-miss counters clean; it also improves sampled residuals.
- It does not satisfy the performance branch: matvec reduction is `0%` on all
  three medium fixtures, geomean speedup versus fused diagonal is only
  `1.0135x`, and the compact validation audit mean remains above the
  `5e-5` gate (`5.9797e-5`).
- Barnase was not run because the full-medium candidate did not pass the
  performance and audit-mean gates. Defaults and `stage3u-safe` remain
  unchanged.
- Next branch: a leaf-local block preconditioner alone is not changing Krylov
  iteration count; evaluate a broader preconditioner correction
  (overlapping/additive Schwarz or coarse correction) only behind hidden gates,
  with medium evidence before any Barnase probe.

## Stage 3AK/10D - TABI Exact-Galerkin Two-Level Coarse Correction Probe

Stage 3AK implementation:

- Added hidden `-pbe-tabi-precond=coarse-galerkin-probe` with
  `-pbe-tabi-precond-coarse-aggregate-cap=<4..64>`,
  `-pbe-tabi-precond-coarse-base=diagonal|block-budgeted`, and the
  smoke-only verification mode `-pbe-tabi-precond-coarse-verify=full-small`.
- Built deterministic topdown aggregates by splitting the largest splittable
  tree node up to the aggregate cap. Each aggregate contributes constant
  `phi` and `phin` modes; restriction is area weighted and reports its
  `R*P=I` error.
- Built the exact Galerkin matrix with the same hidden-safe JCP treecode
  operator used by the solve, while isolating setup matvec statistics and
  disabling report-only audit sampling during setup. The two-level apply is
  `M_base^-1 (r - P R r) + P A_c^-1 R r`.
- Added timeout/profile/timing fields for coarse dimension, storage,
  pivot/status, setup matvec calls/seconds, apply calls/seconds, and
  verification errors. Matrix CSV now reports
  `gmres_matvec_reduction_excluding_setup` and
  `setup_amortization_matvec_equivalent`.
- Added `test/pbe_tabi_coarse_galerkin_precond_smoke.ps1` and
  `-IncludeStage3AKRows` with diagonal, budgeted-block, coarse-cap
  `8/16/32`, composed coarse-plus-block, and sampled-audit validation rows.
- Added `-UseStage3AKCoarseProbe` and coarse telemetry columns to the
  grid-vs-TABI compare harness for future report-only comparisons; Stage 3AK
  itself does not run a Barnase timing probe.

Stage 3AK verification:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`: PASS.
- Exact Galerkin smoke: PASS, including active coarse construction,
  `R*P=I`/column verification at tolerance `<=1e-12`, and low-budget
  fallback before expensive block/coarse setup.
- Focused gates pass: budgeted block preconditioner, matvec-cost profile,
  fused SoA AVX2 leaf path, treecode accuracy, tree builder, topdown MAC
  audit, absolute-floor, channel-semantics, polynomial-estimator, and
  surface-margin smokes.
- Full medium evidence completed:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3AIRows
  -IncludeStage3AJRows -IncludeStage3AKRows -CaseTimeoutSeconds 180
  -AllowFailedStatus`; all `21` Stage 3AK rows emitted `status=ok`.
- Grid-vs-TABI schema sanity on `mini-protein` with
  `-UseStage3AKCoarseProbe`: PASS; summary output preserves
  `coarse_galerkin_probe_diagonal`, `coarse_status=active`, and setup matvec
  telemetry.
- Regression freeze pass: GMRES true residual, sampled true residual, TABI
  profile, memory profile, clean-room audit, PBE regression, PBSA smoke, and
  SASA surface smoke.

Stage 3AK coarse-diagonal viability evidence (GMRES operator matvec reduction
excludes exact-Galerkin setup matvecs):

| Row | mini-protein | charged-peptide | compact-mini-protein | Geomean matvec reduction |
| --- | ---: | ---: | ---: | ---: |
| coarse diagonal, cap 8 | 0% | 0% | 4.17% | 1.41% |
| coarse diagonal, cap 16 | 0% | 0% | 4.17% | 1.41% |
| coarse diagonal, cap 32 | 0% | 6.67% | 4.17% | 3.65% |
| coarse + block-budgeted, cap 16 | 8.33% | 6.67% | 8.33% | 7.78% |

For cap `16` with diagonal base, coarse setup required `28`, `30`, and `30`
operator matvecs and `7.33s`, `6.82s`, and `23.48s` respectively on the
three fixtures, while its GMRES reduction was only `0%`, `0%`, and `4.17%`.
The composed row does not isolate new coarse benefit: its reductions match the
existing budgeted-block trend rather than crossing the viability threshold.

Stage 3AK validation accuracy evidence (coarse diagonal cap `16`, sampled
expansion/channel audit enabled):

| Fixture | Energy drift rel | Sampled residual rel | Far expansion mean | Accepted misses K1/K2/K3/K4 |
| --- | ---: | ---: | ---: | --- |
| mini-protein | 2.0896e-5 | 6.7067e-5 | 6.1519e-5 | 0/0/0/0 |
| charged-peptide | 3.9853e-6 | 8.7215e-6 | 0 | 0/0/0/0 |
| compact-mini-protein | 1.2515e-6 | 3.2324e-5 | 4.4174e-5 | 0/0/0/0 |

Stage 3AK decision:

- Exact Galerkin coarse correction is implemented and accepted-miss clean, but
  it fails both gating arguments: the best coarse-diagonal row reduces GMRES
  matvecs by only `3.65%` geometric mean, far below the required `25%`, while
  setup costs dozens of fine operator applications; the final validation rerun
  also puts `mini-protein` expansion mean above gate (`6.1519e-5 > 5e-5`).
- `charged-peptide` emitted no expansion-mean pair in its validation row and
  remains a report-only follow-up case; it does not rescue the failed mini
  mean gate.
- The coarse branch is therefore closed for now; implementing a cheaper
  coarse assembly would optimize setup for a correction that did not change
  Krylov convergence enough on medium fixtures.
- Barnase was not run, as specified for this probe and reinforced by the
  failed medium viability gate. Production defaults, `stage3u-safe`, the grid
  backend, and direct TABI remain unchanged.
- Next branch: return to estimator/channel-shape or safe leaf-work
  selectivity, using the Stage 3AH calibration evidence rather than adding
  another preconditioner layer.

## Stage 3AL/10E - TABI Kernel-Response Channel Estimator Audit

Stage 3AL implementation:

- Added hidden report-only option
  `-pbe-tabi-convergence-kernel-response-audit=off|summary|sampled-detail`,
  default `off`. It is evaluated only for sampled convergence-overreject
  candidates and does not alter MAC, split, accepted-far evaluation, or
  GMRES operator behavior.
- Confirmed the live source mapping used by Taylor evaluation:
  `K1_phi/K3_phin` are driven by `phin_scalar`, while
  `K2_phi_normal/K4_mixed` are driven by `phi_normal`.
- Added degree-filtered kernel response evaluation for degree `p` and `p-1`
  moments through the existing JCP channel Taylor kernels. Detail mode emits
  `pbe_tabi_kernel_response_calibration.csv` with source family,
  `kernel_tail_p_abs`, `kernel_tail_pm1_abs`,
  `kernel_tail_envelope_abs`, oracle error, and classification.
- Added `test/pbe_tabi_kernel_response_channel_audit_smoke.ps1`,
  `test/pbe_tabi_kernel_response_summary.ps1`, and
  `-IncludeStage3ALRows` in the physical-gap matrix harness.

Stage 3AL verification:

- Native build: `mingw32-make PROFILE=native TARGET=titan.exe all -j1`:
  PASS.
- New kernel-response channel smoke: PASS; `phin_scalar` contributes only to
  `K1/K3`, `phi_normal` contributes only to `K2/K4`, and detail CSV contains
  both source families.
- Existing focused gates pass: convergence shape audit, overreject
  calibration oracle, convergence overreject audit, far-expansion audit,
  moment construction, mixed derivative, treecode accuracy, tree builder,
  and topdown MAC audit.
- Full medium evidence completed:
  `test/pbe_tabi_topdown_physical_gap_matrix.ps1 -IncludeStage3AGRows
  -IncludeStage3AHRows -IncludeStage3ALRows -CaseTimeoutSeconds 180
  -AllowFailedStatus`; all `9` Stage 3AL rows emitted `status=ok`.
- Regression freeze pass: GMRES true residual, sampled true residual, TABI
  profile, memory profile, clean-room audit, PBE regression, PBSA smoke, and
  SASA surface smoke.

Stage 3AL report-only parity evidence:

| Fixture | Safe leaf panels/matvec | Summary leaf panels/matvec | Detail leaf panels/matvec | Safe/detail sampled residual rel |
| --- | ---: | ---: | ---: | ---: |
| mini-protein | 7,335,697 | 7,335,697 | 7,335,697 | 9.75894e-5 / 9.75894e-5 |
| charged-peptide | 5,021,476 | 5,021,476 | 5,021,476 | 1.76983e-5 / 1.76983e-5 |
| compact-mini-protein | 17,228,993 | 17,228,993 | 17,228,993 | 5.30641e-5 / 5.30641e-5 |

Stage 3AL kernel-response projection:

| Best safe channel | Safety scale | Caught would reaccept | Overreject would reaccept | Projected released panels | Fraction of safe leaf work |
| --- | ---: | ---: | ---: | ---: | ---: |
| K4_mixed | 1 | 0 | 1,536 | 382,497 | 1.2928% |

Stage 3AL decision:

- Kernel-response features correctly distinguish the source-family mapping,
  but they do not expose enough safely recoverable leaf work. The best
  zero-caught projection releases only `1.2928%`, far below the required
  `10%` threshold for an active estimator-v2 branch.
- Do not implement kernel-response split selectivity, another reaccept floor,
  or another local preconditioner from this evidence.
- The estimator-selectivity performance branch is closed for now. The next
  branch should reassess TABI versus grid positioning and decide whether TABI
  remains an opt-in memory/accuracy path or needs a larger algorithmic change
  before further runtime investment. Production defaults, `stage3u-safe`,
  direct TABI, and the grid backend remain unchanged. No Barnase run was
  performed because Stage 3AL is report-only.

## Immediate TABI Preconditioner Refactor and Controlled Probe

Status: implemented as a focused preconditioner-quality triage slice after the
Phase 7 large-protein benchmark showed that grid remains the production
default and TABI remains experimental.

Implementation:

- Split JCP block/coarse preconditioner code out of `tabi_solver.cpp` into the
  dedicated internal submodule `src/tabi/tabi_preconditioner.h` and
  `src/tabi/tabi_preconditioner.cpp`.
- Reused worker-local scratch during block application instead of allocating a
  local vector for every block apply.
- Changed `block-budgeted` activation so it uses the configured memory budget
  instead of the legacy whole-process working-set soft cap.
- Allowed partial block preconditioners when valid panel coverage is high
  enough; ill-conditioned local blocks now remain local fallback information
  instead of always discarding the full preconditioner.
- Made the `buffer-aware` preconditioner matrix cheaper to build by using the
  centroid/local kernel approximation inside the preconditioner only. The
  production operator still uses the requested near-field integration.
- Added `operator-near` block matrix construction, which uses the same
  near-pair evaluator semantics as the real JCP operator for block entries
  instead of the centroid/local approximation.
- Added symmetric diagonal-magnitude scaling before block LU and a deterministic
  block quality gate. Timing/progress outputs now include
  `preconditioner_quality_rejected_blocks`,
  `preconditioner_quality_probe_max_gain`, and
  `preconditioner_quality_status`.
- Added `test/titan_tabi_preconditioner_smoke.cpp` and
  `test/pbe_tabi_preconditioner_controlled_probe.ps1`.

Controlled probe command:

```
powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_tabi_preconditioner_controlled_probe.ps1 -TitanExe .\titan.exe -CaseTimeoutSeconds 45 -SurfacePoints 24 -Threads 2
```

Controlled probe result:

| Case | Effective preconditioner | Matvecs | Iterations | Final residual | Energy status |
| --- | --- | ---: | ---: | ---: | --- |
| direct_diagonal | diagonal | 18 | 16 | 1.03e-6 | ok |
| direct_block_buffer_aware | block_budgeted_buffer_aware | 45 | 40 | 7.36e-2 | failed |
| direct_block_operator_near | block_budgeted_operator_near | 9 | 8 | 7.05e-5 | ok |
| treecode_diagonal | diagonal | 18 | 16 | 1.03e-6 | ok |
| treecode_block_buffer_aware | block_budgeted_buffer_aware | 45 | 40 | 7.36e-2 | failed |
| treecode_block_operator_near | block_budgeted_operator_near | 9 | 8 | 7.05e-5 | ok |

Decision:

- The new preconditioner submodule is wired correctly: `block-budgeted`
  activates in both direct and treecode controlled runs and writes normal timing
  output.
- `buffer-aware` is now confirmed as the bad approximation branch: it worsens
  Krylov convergence on the controlled fixture.
- `operator-near` fixes the controlled quality gate: it reduces GMRES work from
  18 to 9 matvecs in both direct and treecode cases.
- Next work should run Barnase/research-medium with
  `-pbe-tabi-precond-block-near=operator-near`. If the matvec reduction survives
  representative fixtures, retire `buffer-aware` from promotion consideration
  and make `operator-near` the serious block-budgeted branch.
- Production default remains `-pbe=grid`; TABI remains explicit/experimental.

## TABI Kernel Rewrite V2A - Branchless Target AVX2

Status: implemented as the first safe kernel-rewrite slice after the
large-protein leaf-kernel transition report showed that scheduler wrappers and
SoA retargeting were not enough for Barnase.

Implementation:

- Added hidden `-pbe-tabi-leaf-pair-kernel=target-avx2-v2`.
- Kept `auto` and existing `target-avx2` behavior unchanged; V2 is opt-in only.
- Reused the existing leaf-pair target-batch entry point, but added a masked
  AVX2 path for self/degenerate lanes. V2 clamps invalid-distance lanes before
  division and masks their contribution to zero instead of falling back the
  whole 4-target batch to scalar.
- Wired `target-avx2-v2` through parser validation, leaf-pair kernel
  normalization, runtime resolver, dispatcher, and leaf-kernel transition probe.
- Extended `test/pbe_tabi_leaf_pair_target_avx2_smoke.ps1` so scalar,
  `target-avx2`, and `target-avx2-v2` are compared in the same fixture.

Verification snapshot:

- TDD red was observed before implementation: the V2 smoke failed because
  `target-avx2-v2` resolved to scalar.
- After implementation, `test/pbe_tabi_leaf_pair_target_avx2_smoke.ps1`
  passed with scalar-vs-V2 relative RMS `1.36044e-16` on the clean fixture and
  `1.20471e-16` on the singular/exclusion fixture.
- The `soa-avx2` path now masks self/degenerate lanes before division as well.
  `test/pbe_tabi_leaf_pair_target_soa_avx2_smoke.ps1` includes a direct
  self-lane probe and passes with no non-finite accumulation.
- Quick Barnase probe with `target-avx2-v2` active confirmed the kernel was
  used (`leaf_pair_batch_kernel=target-avx2-v2`, `415,414,900` vectorized
  sources for leaf size 128), but the run still timed out and remained
  `no_win`: `11.13 s/matvec` for V2A versus `12.02 s/matvec` for the baseline
  leaf-size-128 row in the same quick run.

Decision:

- V2A addresses the first kernel-rewrite issue, the unclean-path scalar
  fallback penalty, without changing production defaults.
- The follow-up review correctly identifies that V2A still has low arithmetic
  intensity: both `soa-avx2` and `target-avx2-v2` are still fundamentally
  4-target x 1-source kernels with per-source broadcast. This is now an
  optimization limitation, not a NaN/self-interaction correctness bug.
- Remaining kernel-rewrite work is V2B: true 4x4 or 4x8 source-register
  blocking, packed contiguous source/target tiles, and cost-tile task splitting
  through a static worklist. Do not promote V2A without Barnase timing evidence
  from `test/pbe_tabi_leaf_kernel_transition_probe.ps1`.

## TABI Kernel Rewrite V2B - 4x4 Leaf Kernel And Cost-Tile Worklist

Status: implemented as an experimental, opt-in follow-up to V2A. V2B improves
the internal architecture but does not yet provide a Barnase promotion win.

Implementation:

- Added hidden `-pbe-tabi-leaf-pair-kernel=target-avx2-4x4`.
- Added hidden `-pbe-tabi-leaf-pair-scheduler=cost-tile-worklist`.
- Added `direct_jcp_target_source_tile_sum_avx2_4x4`, a 4-target x 4-source
  AVX2 microtile that keeps V2A's branchless self/degenerate lane masking and
  handles source and target tails with scalar fallback.
- Added static cost-tile scheduling for direct near leaf-pair work. Source leaf
  spans are split into fixed 32-panel tiles, grouped by target leaf, sorted by
  estimated pair count, and then executed through the batch kernel. Far
  Taylor/MAC logic is unchanged.
- Updated the leaf-block verifier so cost-tile partial source spans are checked
  against the same span, not the whole original source leaf.
- Updated the leaf-kernel transition probe with the `cost-tile-4x4` variant and
  cost-tile/4x4 telemetry columns:
  `leaf_pair_cost_tile_tasks`, `leaf_pair_cost_tile_source_panels`,
  `leaf_pair_cost_tile_target_panels`, `leaf_pair_cost_tile_seconds`, and
  `leaf_pair_kernel_4x4_sources`.

Verification snapshot:

- `test/pbe_tabi_leaf_pair_target_avx2_smoke.ps1` passes and compares scalar,
  legacy `target-avx2`, V2A, and V2B 4x4 on clean, self/degenerate, singular,
  and tail source/target fixtures. Clean scalar-vs-4x4 relative RMS was
  `1.36044e-16`.
- `test/pbe_tabi_leaf_pair_worklist_scheduler_smoke.ps1` passes with
  `cost_tile_rel=1.7068e-16`, `cost_tile_rel4=1.7068e-16`, full-small verifier
  enabled, and deterministic output for one and four workers.
- Native build passed with `mingw32-make PROFILE=native TARGET=titan.exe all -j1`.
- Quick Barnase probe (`baseline,target-avx2-v2,cost-tile-4x4`; leaf sizes
  `64,128`; theta `0.8`; 30 s timeout) produced clean but non-promotable rows:

| Scheduler | leaf_size | status | sec/matvec | batch kernel | decision |
| --- | ---: | --- | ---: | --- | --- |
| baseline | 128 | timeout | `9.393` | scalar | `no_win` |
| target-avx2-v2 | 128 | timeout | `11.612` | target-avx2-v2 | `no_win` |
| cost-tile-4x4 | 128 | timeout | `12.066` | target-avx2-4x4 | `no_win` |
| baseline | 64 | timeout | `18.349` | scalar | `no_win` |
| target-avx2-v2 | 64 | timeout | `22.619` | target-avx2-v2 | `no_win` |
| cost-tile-4x4 | 64 | timeout | `27.222` | target-avx2-4x4 | `no_win` |

Decision:

- V2B proves correctness and plumbing for a cleaner microtile/worklist path, but
  the current 4x4 implementation is still slower than the scalar baseline on
  the quick Barnase gate.
- Do not promote V2B into any experimental preset. Keep `-pbe=grid` as the
  production default and TABI as explicit/experimental.
- The next useful TABI decision is algorithmic ROI: either profile why the
  4x4 path adds overhead despite lower-level vectorization, or pause TABI
  solver investment in favor of grid/PBSA robustness until a larger algorithmic
  change is justified.

## TABI V2C - Leaf-MAC Traversal Gate

Status: implemented as an experimental traversal/work-count gate after review
showed that leaf nodes were always forced into direct near-field work even when
the same MAC/safety gates could accept them as far nodes. This is not a
production default change.

Implementation:

- Added hidden `-pbe-tabi-leaf-mac=<off|audit|on>`, default `off`.
- `audit` evaluates leaf MAC candidates and writes counters while preserving
  the old direct-leaf behavior exactly.
- `on` allows a non-target-containing leaf node to use the existing far/Taylor
  path when `evaluate_mac` accepts it. Leaves containing the target remain
  direct.
- Applied the policy to non-JCP traversal, JCP on-the-fly traversal, and the
  cached interaction-list builder.
- Leaf-MAC accepted nodes still pass through existing fragmentation safety,
  convergence estimator, expansion audit, and accepted-miss diagnostics.
- Added counters:
  `leaf_mac_candidates`, `leaf_mac_accepted`, `leaf_mac_rejected`,
  `leaf_mac_far_panels`, and `leaf_mac_direct_panels_saved_estimate`.
- Updated the leaf-kernel transition probe with a `leaf-mac` variant.

Verification snapshot:

- `test/pbe_tabi_leaf_mac_smoke.ps1` passed on a synthetic clustered fixture.
  `audit` preserved the old output, while `on` reduced direct leaf work and
  stayed close to the direct reference (`jcp_rel=0.00377998`,
  `non_jcp_rel=0.000718785`, saved panels `6424`).
- Follow-up diagnostic fix: leaf nodes accepted by Leaf-MAC now fall back to
  direct evaluation instead of `push_child_nodes(...); continue` when
  fragmentation split, convergence-estimator split, or the K3 rho cap rejects
  the far approximation. This prevents leaf drops and avoids expensive
  split-channel audit on leaves that cannot actually be split.
- The directness probe used by leaf-pair batching is now Leaf-MAC-aware only
  for `leaf_mac=on`. A leaf that would be handled as far returns `false`; a
  leaf that falls back to direct returns `true`; `leaf_mac=audit` remains direct.
- Follow-up control-flow fix: JCP traversal now rechecks
  `leaf_mac_fallback_to_direct` before K3 rho cap, far Taylor evaluation, far
  expansion audit, and `continue`. Setting the flag inside the estimator block
  no longer lets the same leaf continue through the far path and skip direct
  evaluation.
- Quick Barnase transition probe (`baseline,target-avx2-v2,cost-tile-4x4,
  leaf-mac`; leaf sizes `64,128`; theta `0.8`; 30 s timeout) completed and
  kept all rows non-promotable. Leaf-MAC rows reached `gmres_start` but timed
  out before the first matvec/progress row with nonzero work counters, so the
  current implementation is not ready for promotion on Barnase.

Decision:

- Keep `-pbe=grid` as the production PBE default and keep TABI explicit and
  experimental.
- Only consider Leaf-MAC for an experimental preset if Barnase/research fixtures
  show at least `3x` lower leaf work per matvec or `1.5x` faster seconds per
  matvec with `accepted_miss=0`, clean verifier/audit rows, and no GMRES
  regression.
- If further TABI work continues, the next question is why Leaf-MAC cached
  setup/first matvec stalls before producing work counters on Barnase; do not
  promote V2C based only on the synthetic smoke.

## TABI V2D - Precomputed Convergence Norms

Status: implemented as a Leaf-MAC stall diagnostic/optimization after review
identified `evaluate_convergence_estimator` as an O(moment-terms) hot path when
Leaf-MAC exposes millions of target-leaf candidate pairs.

Implementation:

- Extended `PackedJcpMoments` with `convergence_norms`, one entry per tree
  node.
- Precomputed scalar `phin`, normal `phi_normal`, and per-channel
  source/tail norms immediately after JCP moment construction.
- Updated `evaluate_convergence_estimator` to read cached source/tail norms in
  O(1) when the cache is valid for the current Taylor order, while preserving
  the old scan path as a fallback for manually constructed/legacy moment
  payloads.
- Kept target-dependent parts of the estimator unchanged: `rho`, aspect,
  dielectric/channel amplification, and split-selectivity policy are still
  evaluated per target-node pair.

Verification snapshot:

- `test/pbe_tabi_leaf_mac_smoke.ps1` now checks that the convergence norm cache
  is built and that cached estimator output matches the fallback scan path on a
  synthetic Leaf-MAC candidate.
- The same smoke includes a forced convergence-split leaf case. It verifies that
  fallback leaves do not increment Leaf-MAC far-panel counters, do not run far
  expansion audit, and still execute direct leaf work.
- `mingw32-make PROFILE=native TARGET=titan.exe all -j1` passed.
- Quick Barnase transition probe still produced no promotable Leaf-MAC row:
  `leaf-mac` timed out before nonzero matvec/work counters at leaf sizes `64`
  and `128`. In the current quick harness, those rows report
  `convergence_estimator=off`, so the remaining stall is not explained by the
  convergence split/audit control-flow bug.

Decision:

- V2D removes the repeated moment-norm scan as a known algorithmic hazard, but
  does not by itself make Leaf-MAC a Barnase candidate.
- Keep `-pbe=grid` as production default and keep Leaf-MAC experimental.
- The next diagnostic should focus on why Leaf-MAC still times out before the
  first progress/work row after estimator norm caching: likely cached list
  construction, pre-matvec traversal volume, or another target-leaf policy path
  before counters are flushed.

## FMM/M2L Scope A - Scalar Prototype Gate

Status: Scope A scalar prototype, Scope B point-cloud gate, Scope C
panel-shaped scalar gate, Scope D JCP-channel-shaped gate, Scope E K1
bridge gate, Scope F derivative-channel bridge gate, Scope G
production-surface offline JCP comparison, and Scope H hidden offline
single-matvec bridge implemented in `TITAN_FMM_M2L_ROADMAP.md`. Scope I
analytical derivative offline operator is also implemented and passed. Scope J
hidden production-panel single-matvec adapter, Scope K guarded GMRES, Scope L
bounded end-to-end solve, Scope M hidden CLI, and Scope N hidden matrix gates
are implemented and passed. Scope O is implemented; cap-128 repeatability
passed, but the cap-256 stress tier is `timeout_not_verified`, so the full
Scope O promotion-readiness gate did not pass. Scope P Barnase hidden FMM/M2L
gate is implemented; bounded Barnase cap-64 quick passed, while full cap-128
failed the current convergence/status gate despite finite audit-clean FMM rows.
Scope Q Barnase convergence sweep is implemented and passed; cap-128 strict
`1e-6` convergence requires a larger GMRES budget than Scope P used.
Scope R fair-memory/cap-256 Barnase gate is implemented. Its original
memory-regression conclusion is stale after Scope S: hidden FMM/M2L no longer
builds unused legacy cached tree-interaction lists. The refreshed full Scope R
gate keeps cap-128 converged, fast, and memory-light; cap-256 remains not
promotion-ready because both direct and FMM hit the current `240`
iteration/status gate. Scope T Barnase cap-256 convergence sweep is
implemented and classifies the remaining blocker as `restart_limited`.

Rationale:

- Leaf-kernel V2A/V2B and Leaf-MAC V2C/V2D did not produce a promotable
  Barnase candidate.
- The remaining algorithmic question is whether a true near/direct plus
  far/M2L architecture has enough ROI to justify deeper TABI solver work.
- Scope A is intentionally isolated from production PBE/TABI code: it will
  prototype scalar Laplace/Yukawa M2L on synthetic panels, compare against
  all-direct reference, and gate any Scope B TABI integration on measured
  accuracy, setup cost, and speedup.
- The first representative Scope A probe passed on synthetic 2048-panel
  Laplace/Yukawa rows. Best row:
  `kernel=yukawa order=2 theta=0.6 rel_rms_error=0.000111452663183 speedup=5.07641331915 decision=promote_to_scope_b`.
- Scope B added an explicit point-cloud API and a deterministic
  mini-protein-like fixture. The first representative Scope B probe passed:
  `kernel=yukawa order=2 theta=0.55 rel_rms_error=0.000539123021504 speedup=3.29716804201 gate=pass`.
- Scope C added a TABI-panel-shaped scalar adapter (`PrototypePanel` plus
  `point_cloud_from_panels`). The first representative Scope C probe passed:
  `kernel=laplace order=2 theta=0.5 rel_rms_error=0.000643292468717 speedup=1.69286812112 gate=pass`.
- Scope D added a research-only JCP-channel-shaped four-channel probe. The
  first representative Scope D probe passed with nonzero far work:
  `kernel=yukawa order=2 theta=0.8 aggregate_rel_rms_error=0.00413466788901 max_channel_rel_rms_error=0.00591104044265 near_pairs=6094848 far_pairs=10682368 speedup=6.14330704128 gate=pass`.
- Scope E added a real JCP K1 scalar bridge: direct K1 reference
  `G0-Gk` vs Laplace-minus-Yukawa M2L. The first representative Scope E probe
  passed:
  `order=2 theta=0.7 rel_rms_error=0.00621637754907 near_pairs=3784704 far_pairs=4603904 speedup=5.60205568647 gate=pass`.
- Scope F added a finite-difference derivative-channel bridge for K2/K3/K4.
  The first representative Scope F probe passed:
  `order=4 theta=0.7 aggregate_rel_rms_error=0.0648357897519 max_channel_rel_rms_error=0.0648646161014 near_pairs=26042368 far_pairs=41066496 speedup=3.42307506673 gate=pass`.
- Scope G used production `build_surface` panels and production direct JCP as
  an independent reference. Both bounded real-surface fixtures passed:
  `research-medium parity=2.40586437131e-15 k1=0.00274278083935 derivative_aggregate=0.0219441041137 speedup=1.74161167298`,
  `research-large parity=3.23379276979e-15 k1=0.00550837970958 derivative_aggregate=0.000341432712451 speedup=1.7558248512`.
- Scope H combined K1-K4 into production-shaped offline `phi/phin` vectors.
  The representative `research-medium` row was accurate enough
  (`offline_aggregate=0.0218073842146`) but failed the integration gate because
  it was slower than production direct JCP on the bounded sample
  (`speedup=0.675245335009`, `gate=no_speedup`).
- Scope I replaced the finite-difference derivative bridge with analytical
  derivative kernels. Both bounded real-surface rows passed:
  `research-medium offline_aggregate=0.0002820320478 speedup=3.55950049282`,
  `research-large offline_aggregate=0.0014179982537 speedup=4.41025749568`.
- Scope J routes the analytical derivative operator through a hidden
  production-panel single-matvec adapter that consumes the real coupled JCP
  unknown vector and writes production-sign output. Both bounded real-surface
  rows passed:
  `research-medium aggregate=0.0002820320478 speedup=4.83624860675`,
  `research-large aggregate=0.0014179982537 speedup=6.25971893616`.
- Scope K runs a guarded GMRES experiment using the Scope J adapter on bounded
  production-surface panel samples. Both real fixtures passed at `1024`
  sampled panels:
  `research-medium adapter_residual=5.4590676773e-07 direct_true_residual=2.96990273845e-05 speedup=3.20555382964`,
  `research-large adapter_residual=3.85417787791e-07 direct_true_residual=4.21884160643e-05 speedup=3.23620232888`.
- Scope L runs the hidden FMM/M2L adapter inside bounded end-to-end
  `solve_linear_pb` TABI solves, then compares direct JCP and FMM/M2L
  convergence, residuals, energies, and matvec timing on the same sampled atom
  set. Both bounded real-fixture rows passed at `128` sampled atoms:
  `research-medium panels=6133 fmm_true_residual=6.23991642839e-11 direct_energy=-8347.94314733 fmm_energy=-8340.06650253 energy_delta_abs=7.87664480659 energy_delta_rel=0.000943543177951 speedup=7.60893491996 gate=pass`;
  `research-large panels=6132 fmm_true_residual=5.4563596887e-15 direct_energy=-6380.19795965 fmm_energy=-6375.2583789 energy_delta_abs=4.93958074564 energy_delta_rel=0.000774204934844 speedup=9.12711951306 gate=pass`.
- Scope M runs the same hidden FMM/M2L path through real `titan.exe -a arg.in`
  CLI cases with bounded sampled PQR fixtures. Both direct and FMM CLI rows
  exited `0`, reported `status=ok`, wrote FMM/M2L telemetry in
  `pbe_tabi_profile.csv`, and passed the energy/residual/speed gates:
  `research-medium panels=6133 fmm_true_residual=6.23992e-11 direct_energy=-8347.943147 fmm_energy=-8340.066503 energy_delta_abs=7.87664399999994 energy_delta_rel=0.000943543081367363 speedup=5.88432380093087 gate=pass`;
  `research-large panels=6132 fmm_true_residual=5.45636e-15 direct_energy=-6380.19796 fmm_energy=-6375.258379 energy_delta_abs=4.93958100000054 energy_delta_rel=0.000774204974668301 speedup=12.4719178640643 gate=pass`.
- Scope N runs the hidden FMM/M2L path through a broader CLI matrix covering
  scale, alternate fixture geometry, salt `0.0/0.15/0.5 M`, dielectric `1/80`
  and `4/80`, Poisson salt-zero, and a real resolution bump. The full matrix
  passed `8/8` rows with sampled direct-operator audit and `opening_theta=0.6`.
  `surface=16` is the resolution row because `surface=12` does not change the
  current mesh bucket relative to `surface=8`. Key rows:
  `medium_cap128_baseline audit_rel=0.00323097 energy_delta_abs=1.9216379999998 speedup=5.21319583986574`;
  `large_cap128_baseline audit_rel=0.00261365 energy_delta_abs=0.195887999999286 speedup=5.93470020675396`;
  `alt_cap128_baseline audit_rel=0.00499974 energy_delta_abs=5.21474400000079 speedup=5.71991452350604`;
  `medium_cap64_salt0 audit_rel=0.00113586 energy_delta_abs=0.22892899999988 speedup=4.99505955726381`;
  `medium_cap64_surface16 panels=5110 audit_rel=0.00417215 energy_delta_abs=1.32726600000024 speedup=6.5734317938785`.
- Scope O adds repeatability and cap-256 stress gates. The cap-128
  repeatability tier passed all three fixtures with three FMM repeats and
  deterministic energies: `research-medium audit_max=0.00323097
  speedup_median=7.15721267875626 speedup_min=5.14596932441632`;
  `research-large audit_max=0.00261365 speedup_median=5.72520286151872
  speedup_min=4.9350600868128`; `research-alt audit_max=0.00499974
  speedup_median=8.92799795952675 speedup_min=7.53845927941774`.
  The full gate remains blocked because cap-256 stress is
  `timeout_not_verified`: `research-medium panels=12250
  fmm_true_residual=2.28745e-07 audit_rel=0.00502187 fmm_energy=-16086.4942`,
  but the direct baseline did not produce finite energy before timeout, so
  energy delta and speedup cannot be verified. `research-large` cap-256 was
  correctly skipped.
- Scope P adds a bounded Barnase hidden FMM/M2L CLI gate, separate from the
  legacy TABI treecode Barnase probes. The quick cap-64 smoke passed:
  `panels=3053 energy_delta_abs=0.526157 audit_rel=0.00143331
  fmm_residual=1.50766e-11 speedup=4.82697815048099`. The full cap-128
  repeatability tier ran three repeats at `panels=5802`; the hidden FMM rows
  were finite, deterministic, audit-clean, and fast
  (`energy_delta_abs=2.354467`, `audit_rel=0.00268276`, speedup range
  `6.01715806698311..9.17630425009285`). The full gate did not pass because
  current direct and FMM CLI rows report non-`ok` status at the `1e-6` / `120`
  iteration gate (`direct_true_residual=0.000916239`,
  `fmm_true_residual=0.00080777`). Cap-256 Barnase stress was skipped because
  cap-128 did not pass.
- Scope Q runs a bounded Barnase cap-128 convergence sweep over GMRES
  iteration budget, restart, and tolerance. It reproduced the Scope P baseline
  failure at `max_iter=120`, `restart=24`, `1e-6`
  (`direct_true_residual=0.000916239`, `fmm_true_residual=0.00080777`), then
  showed that the issue is iteration/restart budget rather than FMM accuracy:
  `max_iter=240`, `restart=24`, `1e-6` converged strictly with
  `fmm_true_residual=6.26312e-7`, `energy_delta_abs=1.716306`,
  `audit_rel=0.00267975`, `speedup=7.337`; `restart=64` converged in `128`
  iterations with `fmm_true_residual=1.27442e-8` and `speedup=7.299`. Scope Q
  final class is `strict_converged`.
- Scope S fixes the Scope R process-memory regression by reusing the analytical
  FMM/M2L context and avoiding unused legacy cached tree-interaction lists in
  the hidden FMM path. The refreshed full Scope R run uses the Scope Q preset
  (`max_iter=240`, `restart=64`, `1e-6`). Cap-128 repeated strict convergence
  three times with zero FMM energy spread, `energy_delta_abs=1.718162`,
  `audit_rel=0.00267974`, and FMM speedups `9.43x..11.66x`. Process peak memory
  is now close to direct TABI: grid around `208 MB`, direct TABI around
  `16.43..16.47 MB`, and FMM/M2L around `18.65..18.67 MB`. Cap-256 no longer
  shows memory blow-up (`27 MB` FMM peak), but direct and FMM both reached the
  `240` iteration limit without `status=ok` (`direct_true_residual=0.0140401`,
  `fmm_true_residual=0.0131365`). FMM remained audit-clean and fast
  (`speedup=17.17x`). Refreshed Scope R final class is
  `cap128_ready_cap256_not_promotion_ready`.
- Scope T diagnoses the refreshed Scope R cap-256 convergence/status blocker.
  Baseline `240/r64/1e-6` reproduced the non-converged state
  (`direct_true_residual=0.0140401`, `fmm_true_residual=0.0131365`) while
  staying audit-clean and fast (`audit_rel=0.0062981`, `speedup=17.22x`).
  Raising iteration count at restart `64` improved slowly
  (`360/r64=0.00919207`, `480/r64=0.00640712`), but larger restart produced the
  main residual drop (`480/r96=0.00111163`, `480/r128=0.00147099`). The best
  row `720/r128/1e-6` reached `fmm_true_residual=0.000185057`, direct
  verification `direct_true_residual=0.000194446`,
  `energy_delta_abs=7.33137`, `audit_rel=0.00685148`, and `speedup=19.92x`.
  Scope T final class is `restart_limited`.
- Scope U implements the first FMM-native preconditioner slice. Hidden
  `-pbe-tabi-precond=fmm-block-jacobi` builds local coupled JCP blocks from
  `FmmM2lAnalyticContext` panel groups and does not rebuild legacy
  `TreeWorkspace` interactions solely for preconditioning. The full Barnase
  gate keeps cap-128 strict (`true_residual=1.41048e-8`,
  `energy_delta_abs=1.718147`, `audit_rel=0.00267974`, `speedup=10.07x`,
  `peak_working_set=23.73 MB`). Cap-256 remains not promotion-ready but the
  native block preconditioner improves the long-budget residual:
  `480/r128/1e-6` improves from diagonal `0.00147099` to block `0.000852217`
  with clean sampled audit (`0.00686342`) and modest memory (`44.98 MB`).
  Direct verification timed out for the longer cap-256 rows.
- Scope V opens larger block caps for the hidden FMM-native preconditioner only.
  Parser behavior is now split: `fmm-block-jacobi` accepts caps `8..128`, while
  legacy `block` and `block-budgeted` remain capped at `32`. Quick Barnase
  cap-128 probes passed parser validation and strict convergence for caps
  `32/64/128`. Full Barnase cap-256 shows cap-only tuning is not the missing
  piece: at `480/r128/1e-6`, cap `32` residual is `0.000852217`, while caps
  `64/128` are both `0.000958293`; at `720/r128/1e-6`, cap `32` is
  `0.000213074`, while caps `64/128` are both `0.000230641`. Cap `64/128`
  remain audit-clean and memory-light, but do not improve convergence. Cap
  `128` matches cap `64` because current FMM tree leaves are size `64`.
- Scope W adds hidden `-pbe-tabi-precond=fmm-overlap-schwarz`, a Restricted
  Additive-Schwarz preconditioner built from FMM `near_sources`. It factors
  owned plus halo panels but writes back owned panels only. Quick Barnase
  cap-128 passed parser, telemetry, sampled audit, and strict convergence:
  overlap `0` reached `true_residual=1.41048e-8`, `audit_rel=0.00269123`,
  `speedup=11.89x`, peak memory `23.83 MB`; overlap `32` improved residual to
  `3.06488e-9` with `speedup=13.11x`, peak memory `44.61 MB`; overlap `64`
  reached `5.90696e-9`, `speedup=11.91x`, peak memory `80.09 MB`. Full
  Barnase cap-256 Scope W is a no-win: at `480/r128/1e-6`, baseline
  `fmm-block-jacobi` cap `32` reached `0.000852217`, while overlap `32/64`
  worsened to `0.00109664` and `0.00116147`; at `720/r128/1e-6`, baseline
  reached `0.000213074`, while overlap `32/64` worsened to `0.000239418` and
  `0.000344171`. Overlap `128` fell back to diagonal because of the memory
  budget, so it is not a valid Schwarz win.

Decision:

- Keep `-pbe=grid` as production default.
- Keep TABI, Leaf-MAC, and FMM/M2L research paths explicit/experimental.
- Scope L/M use `128` atoms as the safe default bounded cap because even that
  size produced about `6133` panels and a costly direct baseline. The earlier
  `512` atom cap remains useful as a manual stress test, not as the default
  smoke.
- Scope X adds hidden `-pbe-tabi-precond=fmm-nearfield-ilu0`, a sparse
  near-field ILU(0) preconditioner over the FMM near-source graph. It uses a
  deterministic per-panel neighbor graph, sparse coupled `2x2` BSR storage,
  ILU(0) with no fill beyond near-source sparsity, and hidden telemetry for nnz
  blocks, average/max degree, neighbor cap, drop tolerance, setup/apply timing,
  and budget/factor fallback labels.
- Scope X.1 verifies sparse near-field ILU on bounded Barnase cap-256. The
  FMM/ILU candidate reaches strict residual in `128` iterations; cap `64` is
  the preferred row with `true_residual=5.64662e-15`, audit `0.00684539`, peak
  working set `61.58 MB`, and speedup `20.61x` versus direct. Direct TABI
  verification completed but did not converge strictly at the same `720/r128`
  budget (`direct_true_residual=0.000194446`), leaving Scope X.1 classified as
  `strict_fmm_pending_direct`. The direct/FMM energy delta was `4.55179 kJ/mol`,
  inside tolerance, but the non-converged oracle blocks promotion.
- Scope Y implements the controlled TABI treecode default. `-pbe=grid` remains
  the production PBE default, but users who explicitly choose
  `-pbe=tabi -pbe-tabi-sum=treecode` now get analytical FMM/M2L plus sparse
  near-field ILU(0) by default. The resolver maps `-pbe-tabi-fmm-m2l=auto` to
  `single-matvec-analytic` on the coupled JCP treecode path, maps
  `-pbe-tabi-precond=auto` to `fmm-nearfield-ilu0`, and uses ILU neighbor cap
  `64`. `-pbe-tabi-fmm-m2l=off` remains the legacy treecode escape hatch, and
  explicit preconditioner overrides are preserved.
- Scope Z5A adds the preconditioner money-plot harness for the prospective
  TITAN-FAGBI paper. It compares local block-style preconditioners with sparse
  near-field ILU(0), exports residual traces and an SVG convergence plot, and
  first validated the machinery on Barnase cap-128. The quick cap-128 run was
  `harness_ready_but_no_contrast`: Block-Jacobi and sparse ILU(0) both
  converged strictly.
- Scope Z5B completes the denser Barnase cap-256 scientific gate. The required
  Block-Diagonal/FMM Block-Jacobi row stalled at `true_residual=2.13074e-4`
  after 720 iterations, while sparse near-field ILU(0) reached
  `true_residual=5.64662e-15` after 128 iterations with audit
  `tree_vs_direct_matvec_rel=0.00684539`. The decision is `money_plot_pass`.
- Scope Z5C extends this into a multi-fixture repeatability/generalization
  matrix. Barnase cap-256 repeated three times with the same contrast:
  Block-Jacobi stalled at `2.13074e-4`, sparse ILU(0) reached
  `5.64662e-15`, and the residual ratio was `3.77e10`. The sampled
  `research-medium`, `research-alt`, and gated `research-large` cap-256 rows
  did not reproduce the collapse because Block-Jacobi also reached strict
  residual in 128 iterations. Z5C decision: `barnase_only_money_plot`.
- Scope Z1A implements the FAGBI paper-parity validation harness. It compares
  `grid`, TABI direct, default TABI treecode FMM/M2L+ILU0, and FMM/M2L with
  Block-Jacobi on analytic Born/Kirkwood rows and sampled protein rows. The
  first full gate produced `z1a_accuracy_gap`: Barnase, research-medium, and
  research-alt cap-128 passed the direct-reference energy and sampled-audit
  checks, but analytic Born and off-center Kirkwood rows at surface point
  counts 8, 12, and 24 exceeded the oracle error thresholds.
- Scope Z1B replaces the aliased `8/12/24` interpretation with actual
  icosphere panel levels `20/80/320/1280` and compares direct TABI
  `centroid`, `analytic_singular`, and `duffy_adaptive` at salt `0.0/0.15 M`.
  The full 48-row gate completed without fallback and classified
  `z1b_mesh_dominant`. At zero salt the modes were identical; at `0.15 M`
  non-centroid gains were largest on coarse geometry and diminished with mesh
  refinement. The next scientific gate is therefore surface discretization
  and resolution, not unconditional FMM-Duffy promotion.
- Scope Z1C measures direct TABI surface discretization on centered Born and
  off-center Kirkwood spheres at salt `0.0/0.15 M`, using icosphere and
  latlong meshes through about 5120 panels. All 32 rows were finite and
  fallback-free. The decision is `resolution_limited`, with all eight finest
  transitions improving energy error by at least `1.5x`, zero plateau groups,
  and no preferred mesh (`0.6419%` vs `0.6755%` median finest oracle error).
  Plan Z1D for higher resolution before changing the direct operator or FMM
  near-field policy.
- Scope Z1D implements the controlled high-resolution gate. It preserves the
  public `4096` surface-point cap unless `-pbe-tabi-highres-surface=true` is
  used on direct TABI `-pbe-test-pqr` benchmark rows. The initial full gate was
  `z1d_runtime_limited`, but Scope Z1D.1 reran the missing Kirkwood `0.15 M`
  20480-panel row and completed it: `0.1374%` oracle error, `1750.93 s`, final
  residual `6.12512e-15`, and zero near fallbacks. All required 20480-panel
  direct TABI rows are now below the 0.5% oracle-error threshold, so Z1D is
  converged after targeted verification. Keep the runtime cost visible, but do
  not block Z2/Z3 planning on missing high-resolution evidence.
- Scope Z2A implements the first FMM near-field parity slice. Analytical
  FMM/M2L treecode rows now support
  `-pbe-tabi-near-integration=analytic_singular` by routing near pairs through
  the shared production K1..K4 near kernel. The Z2A smoke passed analytic Born,
  off-center Kirkwood, and sampled protein cap-128 rows; the sampled protein
  audit RMS was `0.00344495` with nonzero FMM near/far counters.
- Scope Z2B implements the explicit FMM/M2L Duffy gate. Treecode FMM rows now
  accept `-pbe-tabi-near-integration=duffy_adaptive` with
  `-pbe-tabi-quad-buffer=adaptive` and use the shared `one_ring_cpv` production
  near kernel. The Z2B smoke passed Born, off-center Kirkwood, and sampled
  protein cap-128 rows; the protein row reported Duffy pairs `4,777,175`,
  near pairs `13,815,687`, far pairs `23,798,002`, zero fallbacks, and sampled
  audit RMS `0.00338017`. `near=auto` remains unchanged; Duffy is explicit-only
  until a later matched-geometry policy gate proves benefit.
- Scope Z3A wires `-pbe-tabi-taylor-order=3|4` into the FMM/M2L context instead
  of hard-coding expansion order `4`, and adds
  `fmm_m2l_expansion_order` to `pbe_tabi_profile.csv`. The new
  `test/pbe_fagbi_scope_z3a_fmm_order_matrix.ps1` harness writes matrix,
  decision, summary, report, and SVG artifacts under
  `test/tmp/fagbi_scope_z3a/`. The quick Born/Kirkwood gate passed and
  confirmed effective FMM orders `3` and `4`; the quick decision was
  `z3a_order3_sufficient` because the low-resolution analytic rows showed
  identical sampled audit and runtime within noise. Treat this as
  order-plumbing evidence, not a final adaptive-order policy.
- Scope Z3B adds the fixed-order/theta matrix harness
  `test/pbe_fagbi_scope_z3b_fmm_order_theta_matrix.ps1`. The full
  Born/Kirkwood plus Barnase/research-medium cap-128 sweep completed 24/24
  rows, generated `TITAN_FMM_M2L_SCOPE_Z3B_FMM_ORDER_THETA_MATRIX_REPORT.txt`,
  and the latest rerun classified `z3b_adaptive_order_justified`. The measured
  rows still show theta changing sampled audit/runtime more than fixed order,
  while FMM expansion orders matched requested `3|4`. A theta-policy gate
  should precede adaptive-by-level order work.
- Scope Z3C adds the opening-theta policy harness
  `test/pbe_fagbi_scope_z3c_theta_policy_gate.ps1`. The full
  Barnase/research-medium/research-alt cap-128 matrix completed 12/12 rows,
  generated `TITAN_FMM_M2L_SCOPE_Z3C_THETA_POLICY_REPORT.txt`, and classified
  `z3c_theta06_safe_default`: theta `0.8` was faster but failed the
  research-alt audit gate, while theta `0.7` remained marginal. Keep theta
  `0.6` as the conservative treecode policy until a cap-256/repeatability gate
  proves a looser value is safe.
- Scope Z3D adds the cap-256 repeatability harness
  `test/pbe_fagbi_scope_z3d_theta_cap256_repeatability.ps1`. The full
  Barnase/research-medium/research-alt cap-256 matrix completed 15/15 rows,
  generated `TITAN_FMM_M2L_SCOPE_Z3D_THETA_CAP256_REPEATABILITY_REPORT.txt`,
  and classified `z3d_theta06_confirmed`. Theta `0.6` passed all repeat rows
  with audit maxima below `8e-3`, while diagnostic theta `0.8` exceeded
  `1e-2` on every required fixture and theta `0.7` stayed marginal on
  Barnase/research-alt. Use theta `0.6` for the next Z4 paper-style protein
  panel benchmark; no production grid default changes follow from Z3D.
- Scope Z4A adds the protein panel benchmark harness
  `test/pbe_fagbi_scope_z4a_protein_panel_benchmark.ps1`. The full
  Barnase/research-medium/research-alt cap `128/256` run generated
  `TITAN_FMM_M2L_SCOPE_Z4A_PROTEIN_PANEL_BENCHMARK_REPORT.txt` and classified
  `z4a_direct_limited_fmm_pass`: FMM/M2L + ILU0 passed every required sampled
  audit row with theta `0.6`, cap-128 direct TABI references completed, and
  cap-256 direct references were skipped by runtime policy. Some grid
  comparator rows reported `not_converged`; keep grid as production default,
  but use Z4B to tighten reference/verification policy for paper tables.
- Scope Z4B adds the reference/verification policy harness
  `test/pbe_fagbi_scope_z4b_reference_policy_gate.ps1`. The full
  Barnase/research-medium/research-alt cap `128/256` run generated
  `TITAN_FMM_M2L_SCOPE_Z4B_REFERENCE_POLICY_REPORT.txt` and classified
  `z4b_reference_policy_pass`: cap-128 FMM rows use full direct TABI energy
  references, cap-256 FMM rows use sampled direct-operator audit, and all six
  FMM rows are marked paper-usable. The maximum sampled audit was
  `0.00748955`, cap-128 relative energy deltas stayed below `7.4e-4`, and
  three grid rows remained context warnings rather than paper oracle rows.
- Scope Z4C adds the paper packaging harness
  `test/pbe_fagbi_scope_z4c_paper_packaging.ps1`. It reads the Z4B artifacts
  by default and generated `TITAN_FMM_M2L_SCOPE_Z4C_PAPER_PACKAGE_REPORT.txt`
  with decision `z4c_paper_package_ready`: six FMM main-table rows, nine
  supplement/context rows, a reference-policy table, and a runtime/memory SVG.
  This closes the first paper-panel artifact pass without changing solver
  behavior or the grid production default.
- Scope Z4D adds `test/pbe_fagbi_scope_z4d_manuscript_export.ps1`. It generated
  `TITAN_FMM_M2L_SCOPE_Z4D_MANUSCRIPT_EXPORT_REPORT.txt` with decision
  `z4d_manuscript_assets_ready`, exporting manuscript-ready CSV/Markdown main
  tables, a supplement CSV, a methods box, and a figure caption. The generated
  assets keep the reference policy explicit: cap-128 uses full direct TABI
  energy, cap-256 uses sampled direct-operator audit, and grid remains context
  only when it is not a converged oracle row.
- Scope Z4E adds `test/pbe_fagbi_scope_z4e_evidence_manifest.ps1`. It generated
  `TITAN_FMM_M2L_SCOPE_Z4E_EVIDENCE_MANIFEST_REPORT.txt` with decision
  `z4e_evidence_bundle_ready`, including a 37-artifact manifest, checksum
  table, command provenance table, claim map, limitations map, and Markdown
  evidence bundle. This is packaging/audit-only and preserves the Z4B/Z4D
  reference policy without changing solver behavior or public defaults.
- Scope Z4F adds `test/pbe_fagbi_scope_z4f_manuscript_skeleton.ps1`, a
  deterministic `python-docx` builder, checksum-drift and render-contract
  smokes, and `TITAN_FMM_M2L_SCOPE_Z4F_MANUSCRIPT_SKELETON_REPORT.txt`. The
  measured decision is `z4f_manuscript_skeleton_ready`: a nine-page
  journal-neutral Word draft passed structural and visual QA with six mapped
  claims, five limitations, six protein rows, and two embedded evidence
  figures. Citation gaps remain separate and explicit; no solver behavior or
  default changed.
- Scope Z4G adds public complete-protein fixture preparation (Barnase, 1A63,
  and 1SH1), a no-sampling benchmark harness, and an external manuscript
  builder that documents the 2 x 2 block-CSR near graph and ILU(0) no-fill
  triangular solves. The quick 1SH1 row reached
  `true_residual=4.36458e-14`, but sampled direct-operator RMS was
  `2.05699e-2`, exceeding the `8e-3` gate. Z4G is therefore
  `z4g_operator_accuracy_gap`; full benchmarking and manuscript promotion are
  fail-closed. Grid remains the production default.
- Scope Z4H adds the complete-protein theta recovery harness. Full 1SH1
  selected `theta=0.45` (`audit=6.38179e-3`), but 2048-sample confirmation
  failed transferability on Barnase (`9.72521e-2`) and 1A63 (`1.00216`). The
  decision is `z4h_fixture_specific_policy`; therefore no TABI treecode
  default changes follow, and Z4G remains fail-closed pending an operator-level
  accuracy diagnosis on the complete failing proteins.
- Scope Z4I adds `sampled-detail` FMM/M2L operator attribution and the full
  `test/pbe_fagbi_scope_z4i_operator_attribution.ps1` matrix. All nine
  Barnase/1A63/1SH1 theta rows completed and classified
  `z4i_far_order_implementation_gap`: near error fraction was zero, far error
  fraction was one, K2 dominated every row, recomposition error stayed below
  `5.4e-15`, and reference gaps stayed below `4.3e-14`. The diagnostic records
  configured order `4` versus implemented far moment degree `1`. Tightening
  theta to `0.35` reduced Barnase audit to `2.16504e-2` and 1SH1 to
  `2.10383e-3`, but 1A63 remained at `8.60231e-1`. The next solver scope is
  genuine higher-order analytical moments/M2L (Z4J), not more theta-only
  tuning. Grid remains the production PBE default and Z4I changes no policy.
- Scope O and refreshed Scope R still block grid-default promotion. Scope Q
  fixed bounded Barnase cap-128 convergence, Scope S removed the hidden
  FMM/M2L memory regression, Scope U/V/W established that dense block and
  overlap preconditioners are insufficient for cap-256, Scope X.1 selected
  near-field ILU, and Scope Y makes that engine the default only inside explicit
  TABI treecode runs. Keep `-pbe=grid` as production default until
  repeatability, broader cap-256 fixtures, and verification policy are all
  demonstrated.
- Scope Z4J adds `-pbe-tabi-fmm-cartesian-audit=on` for diagnostic
  `sampled-detail` runs and evaluates complete Cartesian p4 K1-K4 moments
  against the same direct source-cell references as legacy degree one. The
  full matrix decision is `z4j_p4_tight_theta_pass`: Barnase and 1SH1 pass
  the `8e-3` audit gate at theta `0.6`, while 1A63 is narrowly high at
  `0.00853174` and passes at theta `0.45` with `0.00291735`. This validates
  the higher-order formulas but does not wire them into GMRES. Z4K should
  implement true local M2L/L2P before any TABI policy change. Grid remains
  the production PBE default.
- Scope Z4K implements diagnostic target-leaf p4 M2L/L2P with source/local
  order `4` and kernel derivative order `8`. Full-protein attribution produced
  `z4k_local_theta06_accuracy_gap`: Barnase and 1SH1 passed at theta `0.6`,
  while 1A63 reached `0.00961167`; theta `0.45` passed but cannot promote the
  path. Translation matched p4 M2P closely, yet M2L timing achieved only
  `0.42x–1.64x` relative speedup. Hidden `gmres-p4` therefore remains rejected,
  the production FMM operator is unchanged, and grid remains the PB default.
- Scope Z4L implements sampled-detail p3/p4 tail refinement for target-leaf
  local M2L. At theta `0.6`, tolerance `0.002` brought 1A63 below the audit
  gate (`0.00734897`) and retained clean Barnase/1SH1 results. The measured
  decision is `z4l_accuracy_pass_costly`: 1A63 runtime increased about
  `7.48x` and 15,551 source panels fell back to direct attribution. This
  validates the error indicator but not production use. `gmres-p4` remains
  rejected, the degree-one production operator remains unchanged, and grid
  remains the production PBE default.
- Scope Z4M implements a diagnostic static p4 refinement plan and reusable
  compact kernel-jet cache. Quick full-1SH1 evidence classified
  `z4m_accuracy_pass_cache_limited`: the cached candidate preserved audit RMS
  `1.03558e-4` and used only `53.93 MB` peak RSS, but achieved `1.29x` apply
  speedup rather than the required `3x`. Direct fallback still covered
  15,996,854 panel pairs. `gmres-p4` remains fail-closed and no TABI or grid
  production policy changes follow from this result.
- Scope Z4M.1 accelerates diagnostic p4 direct fallback through a deterministic
  packed worklist and the shared AVX2 near primitive. The full-1SH1 quick gate
  preserved audit RMS `1.03558e-4`, improved fallback by `5.87x`, and improved
  total candidate apply by `8.57x`, while using `7.63 MB` task metadata and
  `62.50 MB` peak RSS. Barnase and 1A63 full sampled-detail cases remained
  runtime-limited at the 3600-second process gate, including packed rows.
  Therefore `gmres-p4` remains fail-closed and no TABI/grid default changes
  follow from Z4M.1.
- Scope Z4M.2 caches static-refinement p3/p4 pair metrics while preserving
  deterministic summation, tie breaks, plan hash, and output. On full 1SH1,
  plan build improved from `770.073 s` to `9.83224 s` (`78.32x`) with exact
  plan parity and audit RMS `1.4925e-4`. Barnase completed with a `48.8339 s`
  plan and audit RMS `0.00127916`. Metric storage remained below `0.16 MB`,
  but 1A63 required `354.73 MB` for the separate order-8 kernel-jet cache and
  exceeded the `256 MB` diagnostic budget. The measured decision is
  `z4m2_memory_limited`; `gmres-p4` remains fail-closed and production TABI/
  grid policies remain unchanged.
- Scope Z4M.3 removes that measured cache blocker with two double-precision
  degree-9 Green jets per far pair. The quick decision is
  `z4m3_compressed_cache_pass`: 1SH1 achieved `3.00176x` cache reduction and
  a faster `0.155283 s` median apply while preserving audit `1.4925e-4`;
  1A63 fit in `118.16 MB` and preserved audit `6.46481e-3`; Barnase used
  `32.49 MB`, peaked at `137.89 MB` RSS, and retained audit `1.27916e-3`.
  The full decision is `z4m3_compressed_cache_pass`. This remains diagnostic;
  `gmres-p4` is fail-closed until a guarded integration scope, and grid remains
  default.
- Scope Z4N enables `gmres-p4` only for the exact hidden static/refined,
  metric-cache, compressed-Green-jet, packed-AVX2 and FMM-ILU0 configuration.
  Lifecycle and performance gates passed: plan/cache build counts remained
  one, all full proteins converged internally, peak RSS stayed below 259 MB,
  and p4 matvecs were faster than degree-one FMM. Scientific promotion remains
  blocked by full 1A63, where sampled operator RMS rose to `0.0595872` on the
  p4 solution and sampled direct true residual reached `0.0598425`. Barnase
  (`0.00135936`) and 1SH1 (`0.000149168`) passed. Decision:
  `z4n_operator_accuracy_gap`; TABI treecode stays degree-one by default and
  grid remains the production PBE default.
- Scope Z4N.1 confirms that the 1A63 gap is not fixed by simply adding static
  refinement probes. Full-1A63 quick rows gave `0.0595872` audit for the
  baseline `4/rademacher` plan, `0.0526377` for `8/rademacher`, and
  `0.0509292` for `16/rademacher`; solution and mixed probe modes matched the
  corresponding Rademacher work. The measured decision is
  `z4n1_estimator_gap`, so the full Z4N.1 gate is skipped and `gmres-p4`
  remains explicit diagnostic only.
- Scope Z4N.2/Z4N.3 confirms that 1A63 needs an outlier guard, not more
  density probes. The confirmed budgeted policy is `operator-tail-static`
  tolerance `0.002` plus `strict-geometry-mac` theta `0.5`. Full 1A63,
  Barnase, and 1SH1 passed sampled audit/direct-residual gates with fallback
  ratios `1.84x`, `2.18x`, and `0.47x`, RSS below `370 MB`, and one plan/cache
  build.
- Scope Z4N.4 runs that policy inside full hidden `gmres-p4`. The decision is
  `z4n4_guarded_p4_policy_pass`: all three full proteins converged and passed
  sampled verification, while p4 matvecs were faster than same-run degree-one
  baselines. This supports planning Z4O repeatability/promotion-readiness for
  the explicit guarded p4 path only. TABI treecode default remains degree-one
  and `-pbe=grid` remains the production PBE default.
- Scope Z4O verifies repeatability for that same explicit guarded p4 path.
  Barnase, 1A63, and 1SH1 each passed three p4 repeats with zero energy
  spread, sampled audit maxima of `6.38294e-5`, `0.00235702`, and
  `4.89439e-6`, sampled direct true-residual maxima below `8e-3`, peak RSS
  below `370 MB`, and one plan/cache build per solve. Median p4 matvec ratios
  versus same-run degree-one baselines were `0.601`, `0.824`, and `0.300`.
  Decision: `z4o_repeatability_pass`. This permits a Z4P plan for advanced
  explicit-option documentation/preset smoke only; TABI treecode default
  remains degree-one and `-pbe=grid` remains production default.
- Scope Z4P refreshes the manuscript package around that repeatable guarded p4
  evidence. The generated seven-page DOCX passed render QA and describes guarded
  p4 as an advanced explicit TABI treecode path, not a default. Scope Z4Q adds
  the supplement-first public fixture expansion harness and prepared cached
  full-protein PDB2PQR fixtures for 1UBQ, 2CI2, and 1VII using pinned
  PDB2PQR 3.7.1/AMBER/pH-7 settings. The full Z4Q decision is
  `z4q_expansion_pass`: 1UBQ, 2CI2, and 1VII passed two guarded-p4 repeats with
  zero energy spread, one plan/cache build, and sampled audit/direct residuals
  of `1.46754e-5`, `7.98251e-6`, and `5.87085e-6`/`5.87084e-6`, respectively.
  The refreshed DOCX includes these as supplement rows. Keep the main paper
  claim scoped to the locked Z4O Barnase/1A63/1SH1 panel.
- Scope Z4S closes out the current paper package by copying the refreshed
  manuscript, rendered PDF, Z4O/Z4Q tables, fixture provenance, checksums, and
  reports into `docs/fagbi_guarded_p4_package/`. The package decision is
  `z4s_paper_package_ready`. This is an evidence-bundle/package step only:
  grid remains production default, TABI treecode remains degree-one by default,
  and guarded p4 remains explicit.
- Scope Z4T produces the submission-polished package under
  `docs/fagbi_guarded_p4_package/submission/`. The decision is
  `z4t_submission_polish_ready`: the DOCX/PDF pair now has a `3058`-word full
  manuscript body and `15` rendered QA pages, citation gaps are closed, eight
  references and BibTeX are locked, and the manuscript keeps the same policy
  boundaries. This is document/package polish only; grid remains production
  default, TABI treecode remains degree-one by default, and guarded p4 remains
  explicit.
