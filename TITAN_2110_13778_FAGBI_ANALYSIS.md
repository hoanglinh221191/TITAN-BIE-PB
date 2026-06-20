# TITAN Analysis Of arXiv 2110.13778

Source reviewed: `2110.13778v1.pdf`, "A Cartesian FMM-accelerated Galerkin
boundary integral Poisson-Boltzmann solver" by Chen, Tausch, and Geng.

## Paper In One Page

The paper describes a Galerkin boundary integral Poisson-Boltzmann solver
accelerated by Cartesian FMM. The solver targets the same core problem as
TITAN's TABI/PBE work: solve a dense, nonsymmetric boundary-integral PB system
with GMRES, where a direct matrix-vector product is too expensive.

Core ingredients:

- Well-conditioned derivative-form boundary integral PB formulation with
  coupled surface unknowns: potential and normal derivative.
- Galerkin discretization on a triangulated molecular surface.
- Four coupled kernels `K1..K4` from Coulomb, screened Coulomb, and their
  normal derivatives.
- Duffy/cubature treatment for singular or near-singular panel integrals.
- Cartesian FMM for the boundary matrix-vector product.
- Treecode for final solvation-energy projection to atomic charge locations.
- Block diagonal preconditioner built from local direct interactions in tree
  cubes/clusters.
- Adaptive FMM order by level: coarser levels use higher Taylor order so
  FMM truncation error tracks discretization error.

The important empirical messages:

- For the Kirkwood sphere, Galerkin + better singular treatment improves
  convergence versus older collocation/removal approaches.
- Direct matrix-vector products scale like `O(N^2)`, while FMM shows `O(N)`
  behavior in CPU and memory.
- On protein 1A63, at `N=132,133` triangles, direct matrix-vector solve took
  about 36,530 seconds while FMM took about 37 seconds, with similar solvation
  energy.
- Block diagonal preconditioning substantially reduces GMRES iterations for
  hard protein cases without changing energy meaningfully.

## Mapping To TITAN Today

TITAN is already aligned with the paper in several important ways.

Implemented or partially implemented in TITAN:

- TABI/PBE backend with boundary-integral linear PB machinery.
- Coupled JCP 2013 treecode path with `K1..K4` semantics.
- Analytical FMM/M2L matvec path for production panels and production unknowns.
- `FmmM2lAnalyticContext` cache so geometry/work are not rebuilt every GMRES
  iteration.
- Sparse near-field ILU(0) over the FMM near-source graph.
- Controlled Scope Y behavior: `-pbe=grid` remains the production default, but
  explicit `-pbe=tabi -pbe-tabi-sum=treecode` defaults to analytical FMM/M2L +
  `fmm-nearfield-ilu0`.
- Extensive hidden probe matrix for Barnase, bounded cap-128/cap-256, memory,
  repeatability, and direct/FMM parity.

This means TITAN has moved from the old TABI treecode bottleneck toward the
paper's main algorithmic direction: an FMM-like boundary-integral matvec with
local near-field preconditioning.

## Major Gaps Versus The Paper

### 1. TITAN FMM Is Still A Matvec Adapter, Not Full FAGBI

The paper's solver is a full Cartesian FMM boundary-element algorithm with
near-field direct blocks, moment calculation, upward pass, interaction phase,
downward pass, and final evaluation. TITAN's current FMM/M2L path is a strong
analytical single-matvec adapter integrated into GMRES, but it should be audited
as a full FMM pipeline:

- persistent moment storage per node;
- explicit M2M/M2L/L2L/L2P accounting;
- adaptive order by tree level;
- error estimate tied to discretization and surface resolution.

### 2. Near-Field Integration Mismatch

The paper emphasizes Duffy/cubature singular treatment as a key reason Galerkin
accuracy improves. TITAN currently has advanced near-integration paths, but the
FMM/M2L single-matvec treecode route has historically constrained the FMM path
to centroid-compatible near integration in some code paths.

This is the biggest scientific gap before claiming paper-level parity.

### 3. Galerkin Discretization Parity

The paper is Galerkin with piecewise constant basis and panel integrals.
TITAN's TABI path has JCP/treecode, direct analytical kernels, and surface
panels, but the exact Galerkin-vs-collocation equivalence should be made
explicit in one contract:

- what each matrix entry represents;
- whether entries are panel-integrated or centroid-evaluated;
- which quadrature/singularity policy is used per `K1..K4`;
- how solvation energy projection is computed.

### 4. Preconditioner Direction

The paper's block diagonal preconditioner is local-cluster direct interaction.
TITAN tested that family and found it insufficient for bounded Barnase cap-256.
TITAN's sparse near-field ILU(0) is stronger and appears more appropriate for
the current implementation.

Interpretation: TITAN should keep ILU0 as the practical default for TABI
treecode, but keep block-Jacobi in the benchmark suite because it is the closest
paper analogue.

### 5. Validation Set

The paper validates on:

- Kirkwood sphere;
- protein 1A63;
- a 27-protein set with iteration and CPU comparisons.

TITAN currently has a strong Barnase-centered gate and research-medium/large
sampled gates, but it needs a paper-facing validation matrix:

- Kirkwood analytic sphere over surface resolutions;
- 1A63 or paper-like PDB/PQR fixture;
- a small multi-protein set with `grid`, direct TABI, FMM/M2L+ILU0, and
  optional block-Jacobi comparison.

## Recommended TITAN Roadmap From This Paper

Scopes Z1-Z4 close parity gaps with the FAGBI paper. Scope Z5 is different:
it is the TITAN novelty axis. Its purpose is to show that TITAN is not merely
reimplementing FAGBI, but advancing the preconditioning strategy for dense
protein boundary-integral PB systems.

### Scope Z1 - Paper-Parity Validation Matrix

Goal: make TITAN's TABI treecode default defensible against the paper's
benchmarks.

Tasks:

- Add a Kirkwood sphere FMM/TABI matrix with surface resolution sweep.
- Add 1A63 if fixture can be generated from PDB2PQR-compatible inputs.
- Report solvation energy, GMRES iterations, true residual, sampled direct
  matvec parity, CPU time, memory, and preconditioner effective mode.
- Compare `grid`, TABI direct, FMM/M2L+ILU0, and FMM/M2L+block-Jacobi.

Acceptance:

- FMM/M2L energy tracks direct TABI within the existing audit tolerance.
- Iteration count and wall time are stable across repeated runs.
- Grid remains default unless the broader user-facing gate later says
  otherwise.

Implementation status:

- Z1A harness implemented:
  `test/pbe_fagbi_scope_z1a_parity_matrix.ps1`.
- Z1A full result: `z1a_accuracy_gap`. Barnase, research-medium, and
  research-alt cap-128 protein rows passed TABI direct energy and sampled
  matvec-audit checks, but the analytic centered Born and off-center Kirkwood
  sphere rows at surface point counts 8, 12, and 24 exceeded the oracle error
  thresholds. This is useful paper-parity evidence: the current treecode
  default is operational on protein sampled rows, but analytic accuracy parity
  still needs Z2/Z3 work before claiming FAGBI equivalence.
- Z1A report:
  `TITAN_FAGBI_SCOPE_Z1A_PARITY_VALIDATION_REPORT.txt`.
- Z1B near-field attribution harness implemented:
  `test/pbe_fagbi_scope_z1b_near_field_attribution.ps1`.
- Z1B full result: `z1b_mesh_dominant`. The direct TABI matrix compared
  `centroid`, `analytic_singular`, and `duffy_adaptive` on centered Born and
  off-center Kirkwood spheres at `20/80/320/1280` actual icosphere panels and
  salt `0.0/0.15 M`. All 48 rows completed without near fallback. At zero salt
  the three methods produced identical energies. At `0.15 M`, the largest
  non-centroid improvement was `1.90x` on the coarsest off-center sphere and
  shrank to about `1.03-1.12x` at 1280 panels. Mesh refinement reduced the
  oracle error consistently for both fixtures, reaching about `0.83-1.50%`.
  Therefore the Z1A oracle gap is primarily a surface-resolution/discretization
  issue; centroid near quadrature is a secondary contribution rather than the
  demonstrated sole cause.
- Z1B report:
  `TITAN_FMM_M2L_SCOPE_Z1B_NEAR_FIELD_ATTRIBUTION_REPORT.txt`.
- Z1C surface-discretization harness implemented:
  `test/pbe_fagbi_scope_z1c_surface_discretization.ps1`.
- Z1C full result: `resolution_limited`, with `preferred_mesh=none`. All 32
  direct TABI rows completed without fallback across centered Born and
  off-center Kirkwood spheres, salt `0.0/0.15 M`, `icosphere/latlong`, and
  approximately `80/320/1280/5120` panels. Every finest-mesh transition
  improved the energy error by at least `1.5x`, while no transition met the
  operator-plateau pattern. At the finest level, median oracle error was
  `0.6419%` for icosphere and `0.6755%` for latlong; planar area error was
  `0.1195%` and `0.1586%`, respectively. The two meshes are therefore a quality
  tie at this resolution, and the current evidence supports a Z1D refinement
  study beyond about 5120 panels rather than an immediate collocation rewrite.
- Z1C report:
  `TITAN_FMM_M2L_SCOPE_Z1C_SURFACE_DISCRETIZATION_REPORT.txt`.
- Z1D high-resolution harness implemented:
  `test/pbe_fagbi_scope_z1d_high_resolution_surface.ps1`.
- Z1D full result is now `z1d_resolution_converged_after_targeted_verify`. The
  full matrix first completed seven of eight rows, then Scope Z1D.1 reran the
  missing off-center Kirkwood `0.15 M` / 20480-panel row with the same direct
  TABI `analytic_singular` settings and a longer timeout. The targeted row
  completed in `1750.93 s`, with `0.1374%` oracle error, final residual
  `6.12512e-15`, and zero near fallbacks. All four 20480-panel rows are now
  below the 0.5% oracle threshold (`0.3705%`, `0.2682%`, `0.3874%`, and
  `0.1374%`). Treat high-resolution direct runtime as expensive but no longer
  a missing-evidence blocker for Z1.
- Z1D report:
  `TITAN_FMM_M2L_SCOPE_Z1D_HIGH_RESOLUTION_REPORT.txt`.
- Z1D.1 targeted report:
  `TITAN_FMM_M2L_SCOPE_Z1D_TARGETED_VERIFY_REPORT.txt`.

### Scope Z2 - Duffy/Cubature Near-Field For FMM Path

Goal: integrate higher-order near quadrature into FMM only after the Z1B
attribution result and a mesh-resolution gate justify its cost.

Status:

- Scope Z2A implemented the first parity step:
  `-pbe-tabi-fmm-m2l=single-matvec-analytic` now accepts
  `-pbe-tabi-near-integration=analytic_singular` and routes FMM/M2L near
  pairs through the shared production K1..K4 near kernel.
- Scope Z2B implemented the explicit Duffy/adaptive gate:
  `-pbe-tabi-near-integration=duffy_adaptive` with
  `-pbe-tabi-quad-buffer=adaptive` now runs through the FMM/M2L near-source
  path and reports the shared `one_ring_cpv` near-kernel policy.
- Centroid remains the conservative `auto`/default near policy for treecode
  rows; Z2B does not promote Duffy into `auto`.
- Z2A smoke result:
  `TITAN_FMM_M2L_SCOPE_Z2A_FMM_NEAR_PARITY_REPORT.txt`. Born and off-center
  Kirkwood analytic rows completed with FMM/M2L `analytic_singular`; a
  sampled protein cap-128 row completed with near pairs `13,815,687`, far
  pairs `23,798,002`, and sampled audit RMS `0.00344495`.
- Z2B smoke result:
  `TITAN_FMM_M2L_SCOPE_Z2B_DUFFY_GATE_REPORT.txt`. Born and off-center
  Kirkwood direct/FMM Duffy rows completed with `one_ring_cpv`, zero
  fallbacks, and finite outputs. A sampled protein cap-128 FMM row completed
  with Duffy pairs `4,777,175`, near pairs `13,815,687`, far pairs
  `23,798,002`, and sampled audit RMS `0.00338017`.

Tasks:

- Z1C completed the controlled geometry study through about 5120 panels and
  classified the remaining gap as resolution-limited. Z1D/Z1D.1 then verified
  that all required 20480-panel direct TABI rows cross below 0.5% oracle error.
  The remaining caveat is runtime cost, not observed surface-convergence
  accuracy.
- Implemented for `analytic_singular`: the FMM/M2L centroid-only restriction
  is narrowed, and self-singular near pairs use the same production kernel
  used by direct JCP.
- Implemented for explicit `duffy_adaptive`: FMM/M2L uses the same production
  Duffy/CPV near kernel as direct JCP. Promotion into `near=auto` remains a
  future policy decision and should require matched-geometry benefit, not just
  safety.
- Add per-kernel `K1..K4` audit rows for FMM near-field acceptance/fallback.

Acceptance:

- FMM/M2L with Duffy near-field remains finite and audit-clean on small
  fixtures.
- Kirkwood energy convergence improves versus centroid at matched geometry, or
  the integration is not promoted.
- Performance regression is bounded and explained.

### Scope Z3 - Adaptive FMM Order By Level

Goal: adopt the paper's accuracy/performance tradeoff instead of one global
Taylor order.

Status:

- Scope Z3A implemented fixed-order FMM plumbing and a validation harness:
  `-pbe-tabi-taylor-order=3|4` now controls the FMM/M2L context instead of
  leaving the analytical context hard-coded at order `4`.
- `pbe_tabi_profile.csv` now reports `fmm_m2l_expansion_order`, allowing
  order-specific matrix rows to distinguish legacy treecode order from the
  effective FMM/M2L expansion order.
- Z3A quick result:
  `TITAN_FMM_M2L_SCOPE_Z3A_FMM_ORDER_MATRIX_REPORT.txt`. Born and off-center
  Kirkwood rows at theta `0.6` completed for orders `3` and `4`; reported FMM
  expansion orders matched the request, and the quick decision was
  `z3a_order3_sufficient`. This is order-plumbing evidence only; it does not
  yet replace the planned adaptive-by-level matrix.
- Scope Z3B adds the fixed-order/theta matrix harness
  `test/pbe_fagbi_scope_z3b_fmm_order_theta_matrix.ps1`. The full
  Born/Kirkwood + Barnase/research-medium cap-128 matrix completed all 24 rows
  with effective FMM orders matching requested `3|4`; the latest rerun
  classified `z3b_adaptive_order_justified`. The measured rows still show theta
  controls the main audit/runtime tradeoff more than fixed order in this
  matrix. Adaptive-by-level should therefore wait for a theta-policy decision
  rather than being implemented only because order `4` exists.
- Scope Z3C adds `test/pbe_fagbi_scope_z3c_theta_policy_gate.ps1`. The full
  Barnase/research-medium/research-alt cap-128 matrix completed all 12 rows
  with decision `z3c_theta06_safe_default`: theta `0.8` was about `2.23x`
  faster than theta `0.6` by median matvec time, but research-alt exceeded the
  sampled audit gate (`0.0137388`), and theta `0.7` was still marginal on
  research-alt (`0.00901184`). Keep theta `0.6` as the conservative policy
  until a later cap-256/repeatability confirmation proves otherwise.
- Scope Z3D adds `test/pbe_fagbi_scope_z3d_theta_cap256_repeatability.ps1`.
  The full Barnase/research-medium/research-alt cap-256 gate completed all 15
  rows with decision `z3d_theta06_confirmed`: theta `0.6` repeated cleanly
  three times per fixture with audit maxima below `8e-3`, theta `0.8` exceeded
  `1e-2` on every required fixture, and theta `0.7` remained marginal on
  Barnase/research-alt. This confirms theta `0.6` as the fixed policy for the
  next Z4 paper-style protein panel benchmark.
- Scope Z4A adds `test/pbe_fagbi_scope_z4a_protein_panel_benchmark.ps1`. The
  full Barnase/research-medium/research-alt cap `128/256` protein panel
  completed with decision `z4a_direct_limited_fmm_pass`: FMM/ILU0 used theta
  `0.6`, order `4`, and `analytic_singular` near integration, passed all
  sampled audit gates (`max=0.00748955`), and matched completed cap-128 direct
  references within the energy gate. Direct cap-256 remains runtime-limited,
  so the next paper-evidence step should refine reference policy rather than
  adding adaptive order immediately.
- Scope Z4B adds `test/pbe_fagbi_scope_z4b_reference_policy_gate.ps1` and
  `TITAN_FMM_M2L_SCOPE_Z4B_REFERENCE_POLICY_REPORT.txt`. The full
  Barnase/research-medium/research-alt cap `128/256` gate completed with
  decision `z4b_reference_policy_pass`: cap-128 FMM rows use full direct TABI
  energy references, cap-256 FMM rows use sampled direct-operator audit, and
  all six FMM rows are paper-reference usable. The maximum FMM audit remained
  `0.00748955`; cap-128 relative energy deltas were `0.000216675`,
  `0.000230417`, and `0.000728069` for Barnase, research-medium, and
  research-alt. Three grid rows were context warnings, so grid should stay a
  runtime/memory comparator rather than a paper oracle for this panel.
- Scope Z4C adds `test/pbe_fagbi_scope_z4c_paper_packaging.ps1` and
  `TITAN_FMM_M2L_SCOPE_Z4C_PAPER_PACKAGE_REPORT.txt`. It packages the Z4B
  reference-policy artifacts into a six-row FMM main table, nine-row
  supplement/context table, reference-policy table, and runtime/memory SVG.
  The package classified `z4c_paper_package_ready`: all six FMM rows are
  paper-usable, cap-128 rows retain full direct TABI energy references, and
  cap-256 rows retain sampled direct-operator audit references. This is a
  paper-artifact step only; it does not change TABI/FMM defaults.
- Scope Z4D adds `test/pbe_fagbi_scope_z4d_manuscript_export.ps1` and
  `TITAN_FMM_M2L_SCOPE_Z4D_MANUSCRIPT_EXPORT_REPORT.txt`. It converts the
  Z4C package into manuscript-ready CSV/Markdown main tables, a supplement CSV,
  a methods box, and a figure caption. The export classified
  `z4d_manuscript_assets_ready` with six FMM main rows, nine supplement rows,
  and maximum sampled audit `0.00748955`; cap-128 rows cite full direct TABI
  energy while cap-256 rows cite sampled direct-operator audit.
- Scope Z4E adds `test/pbe_fagbi_scope_z4e_evidence_manifest.ps1` and
  `TITAN_FMM_M2L_SCOPE_Z4E_EVIDENCE_MANIFEST_REPORT.txt`. The evidence
  manifest classified `z4e_evidence_bundle_ready`, hashing `37` required
  artifacts, recording command provenance, mapping six paper-facing claims to
  their Z1-Z5 evidence, and recording five limitations. Z4E is packaging/audit
  only: it introduces no new numerical evidence and keeps grid as production
  default.
- Scope Z4F adds `test/pbe_fagbi_scope_z4f_manuscript_skeleton.ps1` plus a
  bundled-Python DOCX builder. It classified `z4f_manuscript_skeleton_ready`
  after SHA256 verification, structural checks, and visual inspection of all
  nine rendered pages. The journal-neutral Word draft contains a 204-word
  abstract, six protein rows, two evidence figures, a supplementary context
  table, a six-claim trace, and all five Z4E limitations. Bibliographic gaps
  remain explicit in a separate CSV; Z4F does not invent citations or add new
  numerical evidence.
- Scope Z4G replaces sampled internal labels with cached public complete-
  protein fixtures for Barnase, 1A63, and 1SH1. Fixture preparation is pinned
  to PDB2PQR 3.7.1, AMBER, pH 7.0, model 1/chain A; non-positive AMBER
  hydrogen radii are deterministically floored to 0.6 A and recorded in the
  manifest without removing atoms. The quick complete-1SH1 gate converged to
  `true_residual=4.36458e-14`, but its sampled direct-operator RMS was
  `2.05699e-2`, above the required `8e-3`. Z4G is therefore currently
  `z4g_operator_accuracy_gap`; the full matrix and external manuscript are
  intentionally blocked rather than packaging the old atom-sampled results as
  complete-protein evidence.

Tasks:

- Implement level-dependent order `p_l = p_L + L - l` behind an advanced flag
  only after fixed-order protein/order/theta evidence justifies it.
- Record per-level order, term count, and timing.
- Gate on Kirkwood and protein fixtures.

Acceptance:

- Accuracy improves at similar timing, or timing improves at similar accuracy.
- No hidden default change until matrix evidence is good.

### Scope Z4 - Paper-Style Protein Panel Benchmark

Goal: produce a direct analogue of the paper's 1A63 and 27-protein tables.

Tasks:

- Build reusable PQR/surface fixture workflow.
- Add deterministic cap/resolution tiers so CI can run small rows and local
  developer machines can run full rows.
- Export `titan_pbe_fagbi_benchmark.csv`.

Acceptance:

- Benchmark table includes energy, iterations, timing, memory, and error/audit
  columns.
- Rows are reproducible and use `-option=value` arg files.

### Scope Z5 - Preconditioner Collapse And Sparse ILU Breakthrough

Goal: make TITAN's main scientific difference visible: paper-style local
block-diagonal preconditioning can stall on dense protein meshes, while TITAN's
sparse near-field ILU(0) over the FMM near-source graph restores Krylov
convergence.

This should become the "money plot" for a TITAN-FAGBI paper.

Tasks:

- Build an ultra-dense benchmark case, starting with bounded Barnase cap-256.
- Run the same FMM/M2L operator and GMRES budget with:
  - diagonal preconditioner;
  - FMM-native Block-Jacobi as the closest paper analogue;
  - restricted Additive-Schwarz as an intermediate local-overlap control;
  - TITAN sparse near-field ILU(0), default neighbor cap `64`.
- Export per-iteration convergence traces:
  `iteration,preconditioner,true_residual,gmres_residual,matvec_ms,apply_ms`.
- Export final ranking table with iterations to `1e-6`, final true residual,
  wall time, peak memory, preconditioner setup/apply time, and direct/FMM energy
  delta when direct verification finishes.
- Generate a publication-ready log-residual plot:
  `diagonal` versus `block-diagonal` versus `sparse ILU(0)`.

Acceptance:

- Block-Jacobi remains non-strict on the ultra-dense case, documenting the
  local-block limitation.
- Sparse ILU(0) reaches strict residual on the same FMM/M2L operator without
  violating sampled matvec audit or energy gates.
- Peak memory stays within the existing Scope X.1 envelope and is reported
  beside convergence, so the gain is not hiding a memory blow-up.
- The plot is reproducible from committed CSV outputs and can be regenerated
  without manual spreadsheet work.

Implementation status:

- Z5A harness implemented:
  `test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1`.
- Z5A quick Barnase cap-128 result: `harness_ready_but_no_contrast`;
  both Block-Jacobi and sparse ILU(0) reached strict residual, so the paper
  contrast must be measured in Z5B on Barnase cap-256.
- Z5B full Barnase cap-256 result: `money_plot_pass`. The required
  Block-Diagonal/FMM Block-Jacobi row stalled at
  `true_residual=2.13074e-4` after 720 iterations, while TITAN sparse
  near-field ILU(0) reached `true_residual=5.64662e-15` after 128 iterations.
  Both rows were audit-clean (`tree_vs_direct_matvec_rel` about `6.85e-3`).
- Z5C multi-fixture result: `barnase_only_money_plot`. Barnase cap-256
  repeated the same money plot in three repeats, but sampled
  `research-medium`, `research-alt`, and gated `research-large` cap-256 did not
  show Block-Diagonal collapse because Block-Jacobi also reached strict
  residual in 128 iterations. This narrows the publication claim to dense or
  high-coupling Barnase-style meshes unless future fixtures reproduce the same
  contrast.
- Z5 report:
  `TITAN_FMM_M2L_SCOPE_Z5_PRECONDITIONER_MONEY_PLOT_REPORT.txt`.
- Z5C report:
  `TITAN_FMM_M2L_SCOPE_Z5C_MULTI_FIXTURE_MONEY_PLOT_REPORT.txt`.
- Generated paper-figure artifacts:
  `test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_residual_trace.csv` and
  `test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_money_plot.svg`.
  Z5C generalization artifacts live under `test/tmp/fmm_m2l_scope_z5c/`.

Publication claim:

> Purely local block-diagonal preconditioning is insufficient for dense
> protein-scale boundary-integral PB systems because strong near-field coupling
> crosses tree-cell boundaries. TITAN's sparse near-field ILU(0), built from the
> FMM near-source graph, captures that coupling while preserving the accelerated
> FMM/M2L matvec.

## Practical Conclusion For TITAN

The paper strongly supports the direction TITAN has already moved toward:
treecode alone is not enough; a real FMM-style boundary matvec plus local
near-field preconditioning is the right architecture.

However, the paper also raises the bar for scientific parity. TITAN should not
claim "FAGBI-equivalent" yet. The current best description is:

> TITAN TABI treecode now defaults to an analytical FMM/M2L matvec with sparse
> near-field ILU(0), inspired by the FMM-accelerated boundary integral PB
> direction. Full Galerkin/Duffy/adaptive-order paper parity remains a staged
> validation and implementation roadmap.

The strongest independent paper angle is still Scope Z5: the sparse near-field
ILU(0) preconditioner is not a catch-up feature. It is TITAN's differentiating
algorithmic contribution over paper-style block-diagonal preconditioning on the
Barnase cap-256 dense/high-coupling mesh. Z5C shows this should be presented as
a targeted breakthrough first, not yet as a universal protein-fixture result.

## Scope Z4H Full-Protein Accuracy Recovery

Z4H tested whether tightening the FMM opening criterion could recover the
complete-protein operator-audit gap exposed by Z4G. On full 1SH1, decreasing
theta from `0.6` to `0.45` reduced sampled direct-operator RMS from
`2.05699e-2` to `6.38179e-3`, at roughly a threefold wall-time cost in the
full sweep. The highest passing 1SH1 value was therefore `theta=0.45`.

That value did not transfer to the other complete proteins. With 2048 audit
samples, candidate RMS was `9.72521e-2` on Barnase and `1.00216` on 1A63,
while 1SH1 remained within gate at `6.52454e-3`. All GMRES solves reported
small residuals, so this is an operator-parity failure rather than a linear
solver convergence failure. The measured decision is
`z4h_fixture_specific_policy`.

Consequences:

- do not replace the current TABI treecode policy with `theta=0.45`;
- keep Z4G external manuscript promotion fail-closed;
- do not present complete-protein energies for Barnase or 1A63 as validated;
- the next numerical scope must diagnose geometry/interaction-list or
  expansion-accuracy behavior on the failing full proteins, rather than
  extending the manuscript pipeline.

## Scope Z4I Full-Protein Operator Attribution

Z4I implemented `sampled-detail` attribution for the analytical FMM/M2L
operator and completed the full Barnase, 1A63, and 1SH1 matrix at opening
theta `0.6`, `0.45`, and `0.35`. The diagnostic decomposes the sampled
operator error into near/far regions, coupled K1-K4 channels, source-tree
depth, and MAC-ratio bins for the converged solution and deterministic
`phi_only`/`phin_only` probes.

The measured decision is `z4i_far_order_implementation_gap`:

- all nine rows completed with finite output;
- the shared analytical-singular near path reproduced the direct near
  contribution exactly in the sampled rows (`near_error_fraction=0`);
- accepted far cells accounted for the complete measured channel error
  (`far_error_fraction=1`), and K2 was the dominant channel in every row;
- attribution recomposition stayed below `5.4e-15`, while the bridge to the
  existing sampled direct-operator audit stayed below `4.3e-14`;
- the configured expansion order was `4`, but the production far-cell
  operator reported implemented moment degree `1` (monopole/dipole data only).

Tightening theta reduced Barnase audit RMS from `0.288126` to `0.0918591`
and `0.0216504`, but did not reach the `8e-3` gate. On 1A63, the same sweep
gave `0.967063`, `1.00371`, and `0.860231`; the tightest value improved by
only about 11 percent and remained unusable. 1SH1 was less demanding and
reached `0.00638179` at theta `0.45` and `0.00210383` at theta `0.35`.

Therefore, Z4I rules out the near kernel and the direct-reference verifier as
the primary causes of the full-protein gap. The next numerical scope should
be Z4J: implement genuine higher-order analytical source moments and M2L/local
evaluation, beginning with the K2 translation, then repeat the same
attribution gate. Further theta-only tuning is not a sufficient remedy.

## Scope Z4J Cartesian p4 M2P Accuracy Gate

Z4J implemented complete Cartesian source moments through degree four for
K1-K4 and evaluated them only as a `sampled-detail` attribution candidate.
The production GMRES operator remains the legacy degree-one far evaluator.
The full Barnase, 1A63, and 1SH1 matrix at theta `0.6/0.45` classified
`z4j_p4_tight_theta_pass`.

At theta `0.6`, Cartesian p4 reduced sampled operator RMS from `0.288126` to
`0.00366315` on Barnase and from `0.0205699` to `0.000339393` on 1SH1. On
1A63 it reduced RMS from `0.967063` to `0.00853174`, a large correction but
slightly above the strict `8e-3` gate. Tightening theta to `0.45` reduced the
1A63 p4 RMS to `0.00291735`; Barnase and 1SH1 reached `0.00030948` and
`0.000032918`, respectively.

K2 remained the largest residual channel, but its squared error fell by
orders of magnitude on every required row; no channel migration or reference
bridge failure was observed. The p4 evaluator is intentionally expensive
because it performs source-cell M2P independently for sampled targets. The
next scope is Z4K: translate these validated degree-four moments into a true
target-local M2L/L2P path and test whether the accuracy gain can be retained
at production matvec cost. No production policy changes follow from Z4J.

## Scope Z4K Target-Leaf Local M2L/L2P Gate

Z4K implemented p4 source-to-local translation at target-leaf centers using
dedicated order-8 kernel jets (165 terms), followed by p4 L2P. The candidate
remained attribution-only and the full Barnase/1A63/1SH1 theta `0.6/0.45`
matrix classified `z4k_local_theta06_accuracy_gap`; guarded GMRES was not
enabled.

Local-vs-p4-M2P relative RMS stayed between about `1.0e-6` and `1.15e-4`,
confirming the translation and K1-K4 L2P formulas. Barnase and 1SH1 passed the
theta-0.6 audit gate at `0.00420459` and `0.000399779`, respectively, but
1A63 reached `0.00961167`, above the fixed `8e-3` threshold. Theta `0.45`
reduced the 1A63 local audit to `0.00337859` but is diagnostic-only.

The local path also missed the performance gate: measured p4-M2P/local
speedup ranged from `0.42x` to `1.64x`, versus the required `4x`. L2P cost was
negligible; order-8 M2L jet construction dominated. Further work should
optimize or precompute translation operators before reconsidering GMRES, and
must separately address the theta-0.6 1A63 accuracy margin.

## Scope Z4L Adaptive Far-Cell Accuracy Refinement

Z4L added a diagnostic-only p3/p4 local-tail estimator that selectively
splits risky accepted far cells while leaving theta fixed at `0.6`. The full
256-target matrix classified `z4l_accuracy_pass_costly`. On 1A63, tolerance
`0.002` reduced sampled audit RMS from `0.00961167` to `0.00734897`, crossing
the fixed `8e-3` gate. Barnase and 1SH1 remained within gate at `0.00333376`
and `0.00034668`.

The estimator therefore correlates with the direct operator error, but the
winning 1A63 row refined 770 cells, moved 15,551 source panels to direct
attribution, and increased wall time from `261.35 s` to `1954.01 s` (about
`7.48x`). This is an accuracy result, not a production-performance result.
GMRES-p4 remains disabled. Any Z4M work should first reduce direct fallback
and order-8 translation cost without relaxing the theta-0.6 accuracy gate.

## Scope Z4M Static Refinement And Cached p4 Translation Gate

Z4M replaced solution-dependent refinement with a density-independent plan
built from four normalized Rademacher operator probes and added compact
order-8 kernel-jet caching. The quick full-1SH1 comparison classified
`z4m_accuracy_pass_cache_limited`: cached and uncached candidates had the
same sampled audit RMS (`1.03558e-4`), zero recomposition/reference gap, and
the same plan hash. The cache occupied `17.25 MB`, increased peak RSS from
`44.46 MB` to `53.93 MB`, and reduced median full-vector apply time from
`0.760 s` to `0.587 s` (`1.29x`), below the required `3x` gate.

The static plan refined 8,846 cells and left 15,996,854 direct-fallback panel
pairs. The remaining cost is therefore fallback work rather than repeated
kernel-jet construction. Full three-protein verification remains pending,
`gmres-p4` stays fail-closed, and Z4N is not justified by Z4M quick evidence.

## Scope Z4M.1 Packed AVX2 Direct-Fallback Acceleration

Z4M.1 retained the static refinement plan and p4 far operator while packing
direct fallback work into deterministic four-target by 32-source tasks. On
the full 1SH1 quick gate, packed AVX2 preserved the sampled audit RMS
(`1.03558e-4`) and zero recomposition/reference gap. The fallback phase fell
from `0.814163 s` to `0.138782 s` (`5.87x`), while median candidate apply fell
from the same-run uncached `1.35455 s` to `0.157966 s` (`8.57x`). Worklist
metadata used `7.63 MB` and peak RSS remained `62.50 MB`.

The full gate is classified `z4m1_runtime_limited`: Barnase and 1A63 scalar,
uncached, and packed CLI rows all reached the 3600-second process timeout
before sampled-detail attribution completed. This does not invalidate the
1SH1 kernel result, but it blocks Z4N guarded GMRES planning. The next work
must separate candidate-apply timing from the expensive full-solve/direct-
reference envelope or reduce that diagnostic cost without relaxing accuracy.
