# TITAN FMM/M2L Preconditioner Design

## Current Decision

Scope Y changes the operating decision: analytical FMM/M2L plus sparse
near-field ILU(0) is now the default engine for users who explicitly choose
`-pbe=tabi -pbe-tabi-sum=treecode`.

This is not a grid-default promotion. `-pbe=grid` remains the production PBE
default, and TABI direct/non-treecode modes keep their existing behavior.
Legacy TABI treecode remains available as an explicit regression escape hatch
with `-pbe-tabi-fmm-m2l=off`.

Scope U added an FMM-native Block-Jacobi preconditioner. Scope V showed that
larger leaf-local block caps do not solve cap-256 convergence. Scope W showed
restricted Additive-Schwarz is useful on cap-128 but a no-win on cap-256. Scope
X/X.1 showed sparse near-field ILU over the FMM near-source graph is the best
bounded Barnase tradeoff, with cap `64` as the default treecode ILU neighbor
cap.

## Scope S Clarification

Scope S did not prove that block preconditioning is impossible or that the old
preconditioner should be restored wholesale. It showed that hidden FMM/M2L must
not depend on the legacy cached `TreeWorkspace` interaction lists when those
lists are not needed by the FMM matvec.

That memory fix avoided unused legacy tree/worklist allocation and brought the
Barnase cap-128 hidden FMM/M2L peak working set close to direct TABI. The tradeoff
is architectural: the current block/coarse preconditioners are still tied to
legacy `TreeWorkspace` geometry, while the FMM/M2L matvec now owns a separate
`FmmM2lAnalyticContext`.

The correct fix is to decouple preconditioning from the legacy tree and build a
preconditioner directly from the FMM/M2L context.

## Recommended Preconditioner Path

### 1. FMM-Native Block-Jacobi

Status: implemented in Scope U; cap sweep completed in Scope V.

- Build local dense blocks from FMM target leaf/chunk panel groups.
- Use the same coupled JCP local operator convention as the production matrix:
  LHS diagonal plus local direct K1/K2/K3/K4 interactions.
- Factor each block once during setup and apply through the existing GMRES
  `custom_preconditioner` hook.
- Reuse existing block LU, scaling, quality gate, telemetry, and memory-budget
  conventions where possible.
- Do not build the legacy cached `TreeWorkspace` solely for this preconditioner.

Observed value: it improves some cap-256 residuals without regressing memory,
but leaf-local cap tuning alone is not enough for strict cap-256 convergence.

### 2. Overlap / Additive-Schwarz Block-Jacobi

Status: implemented in Scope W; cap-256 gate completed as no-win.

- Extend each FMM leaf block with selected near-source neighbor panels.
- Keep the solved variables restricted to the target leaf while using overlap to
  represent cross-leaf strong coupling.
- Reuse the same dense block factorization and memory-budget telemetry.
- Gate on Barnase cap-256 residual improvement versus Scope V cap `32`.

### 3. Near-Field ILU(0)

Status: implemented in Scope X as the sparse near-field preconditioner used by
the controlled TABI treecode default after Scope Y;
Barnase cap-128/cap-256 gates are measured by the Scope X harness.

- Build a sparse/block-sparse `2x2` BSR matrix from the FMM near-source lists.
- Factor with ILU(0), using no fill beyond the FMM near-source graph.
- Expect better convergence but higher implementation risk and memory cost.

### 4. FMM Multigrid / Coarse Galerkin

This is future research work, not the next engineering slice.

- Use the FMM hierarchy as a coarse space.
- Potentially strong asymptotic behavior, but much larger mathematical and code
  surface area.

## Scope U Contract

Scope U implements the first slice only:

- Keep `-pbe=grid` as production default.
- Keep `-pbe-tabi-fmm-m2l=single-matvec-analytic` hidden/experimental.
- Add hidden `-pbe-tabi-precond=fmm-block-jacobi`.
- Build FMM-native block groups from `FmmM2lAnalyticContext`.
- Gate on Barnase cap-128 quick and cap-256 convergence diagnostics.

FMM/M2L should only be considered for default TABI treecode routing after a
later scope proves cap-256 convergence, accuracy, memory, and repeatability.

## Scope U Implementation Status

Status: implemented and verified as a hidden experimental path.

Implemented:

- Hidden parser value `-pbe-tabi-precond=fmm-block-jacobi`.
- FMM-native block plan from `FmmM2lAnalyticContext` leaf/chunk panel groups.
- Dense local coupled JCP block construction using the existing block LU,
  scaling, quality gate, and telemetry path.
- Solver wiring that builds the FMM context before GMRES and does not build
  legacy `TreeWorkspace` interaction lists solely for the preconditioner.

Full Barnase Scope U result:

- Cap-128: strict converged and audit-clean with
  `true_residual=1.41048e-8`, `energy_delta_abs=1.718147`,
  `audit_rel=0.00267974`, `speedup=10.07x`, and peak working set `23.73 MB`.
- Cap-256: not promotion-ready, but the native block preconditioner improves
  the long-budget residual. At `480/r128/1e-6`, diagonal FMM reached
  `0.00147099`, while `fmm-block-jacobi` reached `0.000852217` with
  `audit_rel=0.00686342` and peak working set `44.98 MB`.

Decision:

- Keep `-pbe=grid` as production default.
- Keep hidden FMM/M2L and `fmm-block-jacobi` experimental.
- Scope V supersedes the original next-step recommendation by completing the
  cap-256 long-budget/block-cap comparison. Do not promote the FMM/M2L treecode
  path until cap-256 strict convergence and verification are solved.

## Scope V Implementation Status

Status: implemented and verified as a hidden experimental cap sweep.

Implemented:

- `fmm-block-jacobi` accepts `-pbe-tabi-precond-block-panel-cap=8..128`.
- Legacy `block` and `block-budgeted` remain capped at `8..32`.
- Scope V CLI gate validates parser behavior, generated `arg.in` syntax, and
  Barnase cap-128/cap-256 telemetry.

Full Barnase Scope V result:

- Quick cap-128: caps `32/64/128` all strict-converged; parser accepted FMM cap
  `128`, rejected legacy cap `128`, and rejected FMM cap `129`.
- Full cap-256: cap `64/128` did not improve over cap `32`.
  At `480/r128/1e-6`, cap `32` residual was `0.000852217`, while cap `64/128`
  were `0.000958293`. At `720/r128/1e-6`, cap `32` was `0.000213074`, while
  cap `64/128` were `0.000230641`.
- Cap `128` matches cap `64` because the current FMM tree leaf size is `64`.

Decision:

- Keep FMM/M2L and `fmm-block-jacobi` hidden.
- Do not pursue larger leaf-local caps as the main convergence strategy.
- Scope W showed cross-leaf overlap/additive-Schwarz is not sufficient on
  cap-256; Scope X now targets sparse near-field ILU over the FMM near-source
  graph.

## Scope W Implementation Status

Status: implemented as a hidden experimental Restricted Additive-Schwarz path;
Barnase cap-128 quick verified; Barnase cap-256 full gate completed as no-win.

Implemented:

- Hidden parser value `-pbe-tabi-precond=fmm-overlap-schwarz`.
- Hidden halo option `-pbe-tabi-precond-overlap-panel-cap=0..128`.
- FMM context overlap groups with owned panels from target leaf/chunks and halo
  panels selected from FMM `near_sources`.
- Restricted dense block apply: factor `owned + halo`, but write back owned
  panels only.
- Shared preconditioner builder now supports `owned/factor` panel groups while
  legacy block paths keep `owned=factor`.
- Scope W CLI harness and report outputs.

Quick Barnase Scope W result:

- Cap-128, `240/r64/1e-6`, owned cap `32`.
- Overlap `0`: `true_residual=1.41048e-8`, `audit_rel=0.00269123`,
  `speedup=11.89x`, peak working set `23.83 MB`.
- Overlap `32`: `true_residual=3.06488e-9`, `audit_rel=0.00269123`,
  `speedup=13.11x`, peak working set `44.61 MB`.
- Overlap `64`: `true_residual=5.90696e-9`, `audit_rel=0.00269123`,
  `speedup=11.91x`, peak working set `80.09 MB`.

Full Barnase Scope W cap-256 result:

- Baseline `fmm-block-jacobi` cap `32` at `480/r128/1e-6`:
  `true_residual=0.000852217`, `audit_rel=0.00686342`.
- Overlap `32/64` at `480/r128/1e-6` worsened residual to
  `0.00109664` and `0.00116147`.
- Overlap `128` exceeded the current memory budget and fell back to diagonal,
  giving `0.00147099`.
- Baseline `fmm-block-jacobi` cap `32` at `720/r128/1e-6`:
  `true_residual=0.000213074`, `audit_rel=0.00685375`.
- Overlap `32/64` at `720/r128/1e-6` worsened residual to
  `0.000239418` and `0.000344171`.
- Overlap `128` again fell back to diagonal, giving `0.000185057`; this is not
  a valid overlap-Schwarz win because `preconditioner_effective` was
  `diagonal_fmm_overlap_budget_cap`.

Decision:

- Scope W verifies the restricted overlap machinery and telemetry.
- Do not promote FMM/M2L or overlap Schwarz publicly.
- Move to sparse near-field ILU instead of further leaf-local or overlap cap
  tuning.

## Scope X Implementation Status

Status: implemented as a hidden experimental Sparse Near-Field ILU(0) path;
Barnase cap-256 full gate completed as `strict_fmm_pending_direct`.

Implemented:

- Hidden parser value `-pbe-tabi-precond=fmm-nearfield-ilu0`.
- Hidden ILU controls:
  `-pbe-tabi-precond-ilu-neighbor-cap=32..1024` and
  `-pbe-tabi-precond-ilu-drop-tol=0..1e-2`.
- Deterministic FMM near-source graph export from `FmmM2lAnalyticContext`.
- Sparse coupled `2x2` BSR ILU(0) assembly/factor/apply using local analytical
  JCP `K1/K2/K3/K4` blocks plus LHS diagonal.
- Solver wiring that builds ILU after the FMM context and before GMRES without
  rebuilding legacy `TreeWorkspace`.
- Hidden sparse telemetry: nnz blocks, average/max degree, neighbor cap, drop
  tolerance, estimated bytes, setup/apply timings, and explicit fallback labels.
- Scope X CLI harness and report outputs.

Measured Barnase Scope X.1 result:

- Cap-128 quick: all ILU rows were strict-verified against direct. Cap `64`
  was the best practical row with `true_residual=2.44098e-15`, audit
  `0.00269123`, energy delta `1.718295 kJ/mol`, speedup `8.91x`, and peak
  working set `34.08 MB`.
- Cap-256 full: baseline `fmm-block-jacobi` at `720/r128/1e-6` remained
  non-strict (`true_residual=0.000213074`), while `fmm-nearfield-ilu0`
  reached strict residual in `128` iterations for every screened cap. Cap `64`
  is the best default candidate because it reached `true_residual=5.64662e-15`,
  audit `0.00684539`, `644416` sparse blocks, peak working set `61.58 MB`,
  and speedup `20.61x` versus the direct verifier.
- Direct TABI verification at `720/r128/1e-6` completed without timeout but did
  not report `status=ok` (`direct_true_residual=0.000194446`). The FMM/direct
  energy delta for cap-64 ILU was `4.55179 kJ/mol`, inside the planned energy
  tolerance, but the non-converged direct oracle keeps the official class at
  `strict_fmm_pending_direct` rather than `strict_verified`.

Decision policy after Scope Y:

- Keep `-pbe=grid` as production default.
- Use analytical FMM/M2L plus `fmm-nearfield-ilu0` by default when users
  explicitly choose TABI treecode.
- Keep `-pbe-tabi-fmm-m2l=off` as the legacy treecode regression escape hatch.
- Treat local near-field preconditioning as solved for bounded Barnase
  cap-256. The next gates should validate repeatability and broader cap-256
  fixtures with ILU cap `64`, while keeping grid as the production default until
  an independent verification policy replaces or supplements the slow direct
  oracle.

## Scope Y Implementation Status

Status: implemented as the controlled TABI treecode default.

Implemented:

- `-pbe-tabi-fmm-m2l` now accepts `auto|off|single-matvec-analytic` and
  defaults to `auto`.
- For `-pbe=tabi -pbe-tabi-sum=treecode`, `auto` resolves to
  `single-matvec-analytic` on the coupled JCP treecode path.
- If FMM/M2L is active and `-pbe-tabi-precond` is `auto`, the resolver selects
  `fmm-nearfield-ilu0` and default ILU neighbor cap `64`.
- Explicit overrides are preserved: `-pbe-tabi-fmm-m2l=off` keeps the legacy
  treecode path, and `-pbe-tabi-precond=diagonal` keeps diagonal
  preconditioning even with the FMM/M2L matvec.

Smoke evidence:

- `test/pbe_fmm_m2l_scope_y_treecode_default_probe.ps1 -Quick` verifies that
  the omitted-FMM treecode arg file resolves to
  `fmm_m2l_mode=single-matvec-analytic`,
  `preconditioner_effective=fmm_nearfield_ilu0`, and ILU neighbor cap `64`.
- The same smoke verifies the legacy escape row reports `fmm_m2l_mode=off` and
  direct TABI does not auto-enable FMM/M2L.
