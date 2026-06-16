# TITAN FMM/M2L Scope A Roadmap

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a small, benchmark-gated scalar FMM/M2L prototype for TABI research before wiring anything into the production PBE solver.

**Architecture:** Scope A is an isolated research prototype: it builds a cell tree over synthetic boundary panels, keeps near interactions exact/direct, and evaluates far interactions with a cached M2L-style operator for scalar Laplace/Yukawa kernels. It does not change `-pbe=grid`, TABI defaults, JCP channels, PBSA outputs, or production runtime behavior.

**Tech Stack:** C++23, existing TITAN TABI geometry conventions, `TitanExecutor`, aligned vectors where useful, PowerShell benchmark harnesses, and CSV/report outputs under `test/tmp`.

---

## Why This Exists

Recent TABI work showed that micro-kernel and Leaf-MAC changes do not produce a promotable Barnase candidate. The current Leaf-MAC path tries to push millions of near leaf interactions into per-target analytical Taylor evaluation, which is not the usual FMM shape and remains too expensive or fragile for real workloads.

Scope A changes the question from "can the current treecode leaf path be patched?" to "can a true M2L architecture reduce far-field work with acceptable accuracy and setup cost?" The first prototype must be scalar and synthetic so the answer is not hidden behind JCP channel complexity, GMRES behavior, or protein-specific preprocessing.

## Non-Goals

- Do not change `-pbe=grid` as the production default.
- Do not promote TABI or any FMM path to a public preset.
- Do not wire the prototype into the JCP TABI matvec in Scope A.
- Do not add a third-party FMM runtime dependency.
- Do not approximate near-field panels with M2L. Near/direct remains exact.
- Do not optimize AVX2/FMA kernels before the algorithmic gate passes.

## Scientific And Numerical Contract

- Prototype kernels:
  - Laplace: `G(r) = 1 / r`.
  - Yukawa: `G(r) = exp(-kappa * r) / r`.
- Near-field direct reference is the source of truth.
- Far-field approximation must be compared against all-direct for synthetic fixtures.
- M2L acceptance is based on well-separated cell-cell pairs, not target-leaf opportunistic MAC alone.
- The prototype must report both absolute and relative errors:
  - potential RMS error;
  - max absolute error;
  - relative RMS error against direct;
  - accepted far pair count;
  - direct near pair count;
  - M2L pair count;
  - setup seconds;
  - apply seconds.

## Files

- Create: `src/tabi/tabi_fmm_m2l_prototype.h`
  - Internal prototype types: synthetic panel source/target arrays, tree cells, interaction lists, operator cache keys, and run configuration.
- Create: `src/tabi/tabi_fmm_m2l_prototype.cpp`
  - Scalar direct reference, tree build, near/far list build, M2L apply, validation metrics, and benchmark helper entry points.
- Create: `test/titan_fmm_m2l_prototype_smoke.cpp`
  - C++ smoke for direct reference, tree interaction partition, M2L accuracy, Yukawa mode, deterministic threading, and setup/apply counters.
- Create: `test/fmm_m2l_scope_a_smoke.ps1`
  - PowerShell wrapper that builds/runs the smoke and writes `test/tmp/fmm_m2l_scope_a`.
- Create: `test/pbe_fmm_m2l_scope_a_probe.ps1`
  - Synthetic benchmark harness with small/medium panel counts and generated CSV/report.
- Create: `TITAN_FMM_M2L_SCOPE_A_REPORT.txt`
  - Generated report from the probe, committed only when it contains a useful decision snapshot.
- Modify: `.gitignore`
  - Whitelist any new smoke scripts or C++ smoke files under `test/` if the broad `/test/` rule hides them.
- Modify: `TITAN_PBE_TABI_ROADMAP.md`
  - Point from the TABI Leaf-MAC/TABI closure section to this FMM/M2L Scope A roadmap.

## Phase A1 - Prototype Boundary And Tests

- [ ] Add `src/tabi/tabi_fmm_m2l_prototype.h` with a small API that cannot be called from production TABI code by accident.

  Required public-internal types:

  ```cpp
  namespace titan::tabi::fmm_m2l {

  enum class KernelKind {
    Laplace,
    Yukawa,
  };

  struct PrototypeConfig {
    KernelKind kernel = KernelKind::Laplace;
    double kappa = 0.0;
    int max_leaf_size = 64;
    int expansion_order = 4;
    double separation_theta = 0.6;
    int threads = 1;
  };

  struct PrototypeMetrics {
    std::size_t sources = 0;
    std::size_t targets = 0;
    std::size_t cells = 0;
    std::size_t near_pairs = 0;
    std::size_t far_pairs = 0;
    std::size_t m2l_operators = 0;
    double setup_seconds = 0.0;
    double apply_seconds = 0.0;
    double rms_abs_error = 0.0;
    double max_abs_error = 0.0;
    double rel_rms_error = 0.0;
  };

  PrototypeMetrics run_synthetic_probe(const PrototypeConfig& config,
                                       std::size_t panel_count,
                                       unsigned seed);

  }  // namespace titan::tabi::fmm_m2l
  ```

- [ ] Add `test/titan_fmm_m2l_prototype_smoke.cpp` with a direct-only failing test first.

  Required first checks:

  ```cpp
  // Build two deterministic point clouds.
  // Compute all-direct Laplace potential.
  // Run prototype with theta so strict that all interactions remain near/direct.
  // Assert rel_rms_error <= 1e-14 and far_pairs == 0.
  ```

- [ ] Run the smoke before implementation and confirm it fails because the prototype symbols do not exist.

  Command:

  ```powershell
  mingw32-make PROFILE=native TARGET=test\titan_fmm_m2l_prototype_smoke.exe test\titan_fmm_m2l_prototype_smoke.exe -j1
  ```

  Expected: compile failure for missing `tabi_fmm_m2l_prototype` symbols.

## Phase A2 - Direct Reference And Tree Partition

- [ ] Implement deterministic synthetic panel generation and all-direct scalar reference.

  Requirements:
  - frame-independent deterministic RNG from `seed`;
  - SoA arrays for `x/y/z/q`;
  - no random values inside the timed apply loop;
  - `r2 == 0` contributes zero for self-interaction.

- [ ] Implement a simple octree or balanced AABB tree for synthetic points.

  Requirements:
  - leaf size controlled by `PrototypeConfig::max_leaf_size`;
  - deterministic child/cell ordering;
  - cell center, radius, first/last point span;
  - no production TABI tree mutation.

- [ ] Implement near/far cell pair partition.

  Rule:

  ```text
  accept_far if (source_radius + target_radius) / center_distance <= separation_theta
  otherwise descend or mark near leaf pair
  ```

  Safety:
  - identical or overlapping cells are always near;
  - empty cells are skipped;
  - near/far pair counts are exported in metrics.

- [ ] Extend the smoke to assert that a two-cluster fixture produces both near and far pairs.

  Expected:
  - `near_pairs > 0`;
  - `far_pairs > 0`;
  - all targets receive finite potentials.

## Phase A3 - Scalar M2L Prototype

- [ ] Implement a simple interpolation/equivalent-point M2L approximation.

  Minimum acceptable Scope A design:
  - place fixed equivalent source points on a cube or sphere around each source cell;
  - solve small dense least-squares or interpolation weights once per source cell;
  - translate equivalent representation to target cell check points;
  - evaluate targets from target-local representation.

  Constraint:
  - if this becomes too large, use a deterministic low-order multipole expansion for Laplace first and keep the same API. The report must name the exact approximation used.

- [ ] Cache M2L operators by relative cell geometry bucket.

  Required key fields:
  - expansion order;
  - kernel kind;
  - rounded relative offset;
  - rounded source/target radius ratio bucket.

  Metrics:
  - `m2l_operators` is the number of unique cached operators;
  - setup time includes operator creation;
  - apply time excludes direct reference.

- [ ] Add Yukawa support with `kappa > 0`.

  Required smoke:
  - compare all-direct Yukawa reference and M2L output;
  - ensure all outputs finite;
  - require `rel_rms_error < 1e-2` for a loose first gate.

## Phase A4 - Parallel Apply And Determinism

- [ ] Parallelize independent target-cell or target-block apply work with `TitanExecutor`.

  Rules:
  - no nested executor calls;
  - worker-local accumulation when two tasks can touch the same target;
  - deterministic reduce order;
  - `threads < 1` resolves to `1`.

- [ ] Add deterministic threading smoke.

  Required assertion:

  ```text
  nt=1 and nt=4 produce identical metrics and potentials within 1e-13 relative RMS.
  ```

- [ ] Keep AVX2 out of Scope A unless scalar M2L clearly passes the algorithmic gate.

## Phase A5 - Probe Harness And Decision Gate

- [ ] Add `test/pbe_fmm_m2l_scope_a_probe.ps1`.

  Required options:

  ```powershell
  param(
    [string]$TitanExe = ".\titan.exe",
    [switch]$Quick,
    [int[]]$PanelCounts = @(2048, 8192, 32768),
    [int[]]$Orders = @(2, 4, 6),
    [double[]]$Theta = @(0.45, 0.6, 0.8),
    [int]$Threads = 8
  )
  ```

  Output:

  ```text
  test/tmp/fmm_m2l_scope_a/fmm_m2l_scope_a_matrix.csv
  test/tmp/fmm_m2l_scope_a/fmm_m2l_scope_a_summary.dat
  TITAN_FMM_M2L_SCOPE_A_REPORT.txt
  ```

- [ ] CSV columns must include:

  ```text
  panel_count,kernel,order,theta,status,setup_seconds,apply_seconds,
  direct_seconds,near_pairs,far_pairs,m2l_operators,rms_abs_error,
  max_abs_error,rel_rms_error,speedup_vs_direct,decision_class
  ```

- [ ] Decision classes:

  ```text
  accuracy_fail
  setup_limited
  no_speedup
  promising_scalar
  promote_to_scope_b
  ```

- [ ] Scope A promotion rule:

  `promote_to_scope_b=true` only when all are true on at least one medium synthetic fixture:

  - `rel_rms_error <= 1e-3` for Laplace;
  - `rel_rms_error <= 3e-3` for Yukawa;
  - `apply_seconds <= 0.5 * direct_seconds`;
  - `setup_seconds <= direct_seconds`;
  - no non-finite potentials;
  - deterministic `nt=1` vs `nt=4`.

## Phase A6 - Documentation And Handoff

- [ ] Update `TITAN_PBE_TABI_ROADMAP.md` with a section:

  ```markdown
  ## FMM/M2L Scope A - Scalar Prototype Gate

  Status: planned.
  Rationale: Leaf-MAC and leaf-kernel micro-optimizations did not produce a
  Barnase candidate; Scope A tests whether a true near/direct plus far/M2L
  architecture has algorithmic ROI before touching production TABI JCP.
  ```

- [ ] Update this roadmap with verification evidence after Scope A runs.

  Required evidence:
  - native build command and status;
  - C++ smoke status;
  - PowerShell probe status;
  - best row table;
  - final decision: stop, tune prototype, or move to Scope B.

## Scope B Entry Criteria

Do not start Scope B unless Scope A produces `promote_to_scope_b=true` or the user explicitly overrides the gate.

If Scope B starts, its goal is not a production solver. It should wire the scalar M2L architecture into a controlled TABI research path on one mini-protein fixture, then decide whether JCP channel support is scientifically and computationally justified.

## Scope B - Mini-Protein-Like Point-Cloud Gate

Goal: prove the Scope A scalar M2L engine can consume an externally supplied
panel/point cloud instead of only its internal two-cluster synthetic fixture.

Boundary:

- This is still research-only and does not call the production TABI/JCP matvec.
- The integration object is an explicit point cloud API:
  `run_point_cloud_probe(config, PrototypePointCloud)`.
- The first fixture is deterministic `mini-protein-like`: a folded/helical
  charged surface-like point cloud. It is not a scientific protein surface
  replacement, but it exercises elongated geometry and uneven charge patterns
  before JCP channels are considered.

Implementation checklist:

- [x] Add public-internal `PrototypePointCloud`.
- [x] Add `run_point_cloud_probe(...)`.
- [x] Add deterministic `make_mini_protein_like_point_cloud(...)`.
- [x] Extend the C++ smoke to validate mini-protein-like finite output,
  near/far partitioning, and low-order approximation quality.
- [x] Add a Scope B probe/report:
  `test/pbe_fmm_m2l_scope_b_probe.ps1` and
  `TITAN_FMM_M2L_SCOPE_B_REPORT.txt`.

Scope B gate:

- `rel_rms_error <= 5e-2` for mini-protein-like Laplace and Yukawa.
- `apply_seconds < direct_seconds`.
- `setup_seconds <= direct_seconds`.
- finite outputs and deterministic smoke pass.
- If the gate passes, the next step may design Scope C: controlled TABI surface
  panel extraction or one mini-protein JCP research matvec. Scope C still must
  be opt-in and hidden.

## Scope C - TABI-Panel-Shaped Scalar Gate

Goal: consume TABI `Panel`-shaped data through a scalar M2L adapter before any
JCP-channel or GMRES integration is attempted.

Boundary:

- Scope C is still research-only.
- It does not call `apply_jcp_treecode_operator`.
- It does not replace or alter production TABI treecode traversal.
- It treats a panel as a scalar source at its center with source strength
  derived from `area * scalar_density`. This is a geometry and work-count gate,
  not a scientific PB boundary operator.

Implementation checklist:

- [x] Add a lightweight public-internal `PrototypePanel` with center, normal,
  area, and scalar density fields.
- [x] Add `point_cloud_from_panels(...)`.
- [x] Add deterministic `make_tabi_panel_like_point_cloud(...)` or equivalent
  panel fixture with uneven areas and normals.
- [x] Extend the C++ smoke to validate the panel adapter.
- [x] Add `test/pbe_fmm_m2l_scope_c_probe.ps1` and
  `TITAN_FMM_M2L_SCOPE_C_REPORT.txt`.

Scope C gate:

- `rel_rms_error <= 5e-2` for Laplace and Yukawa on the panel-like fixture.
- `apply_seconds < direct_seconds`.
- `setup_seconds <= direct_seconds`.
- finite outputs.
- If Scope C passes, the next step is a separate Scope D design for one
  controlled TABI JCP research matvec. Scope D must define source-channel
  semantics before code changes.

## Scope D - JCP-Channel-Shaped Research Matvec Gate

Goal: test whether the scalar M2L engine can handle a multi-channel
JCP-shaped source layout on TABI-panel-like geometry before production JCP
treecode is touched.

Boundary:

- Scope D is still research-only.
- It does not call production `apply_jcp_treecode_operator`.
- It does not replace direct leaf kernels, GMRES, or preconditioners.
- It defines four channel-shaped scalar source arrays inspired by the current
  JCP direct path:
  - `k1_phi_from_phin`
  - `k2_phi_from_phi_normal`
  - `k3_phin_from_phin`
  - `k4_phin_from_phi`
- Each channel is evaluated with the same near/direct plus far/M2L machinery
  and compared against all-direct scalar reference. This is a channel-layout
  and work-count gate, not a full PB boundary-operator validation.

Implementation checklist:

- [x] Add `JcpChannelPrototypeMetrics`.
- [x] Add deterministic channel source construction from `PrototypePanel`.
- [x] Add `run_jcp_channel_probe(...)`.
- [x] Extend C++ smoke to verify finite channel outputs, deterministic metrics,
  and channel aggregate error.
- [x] Add `test/pbe_fmm_m2l_scope_d_probe.ps1` and
  `TITAN_FMM_M2L_SCOPE_D_REPORT.txt`.

Scope D gate:

- aggregate channel relative RMS `<= 5e-2`.
- every finite channel output must be finite.
- far-pair work must be nonzero; `no_far_work` is not a valid pass.
- total channel apply seconds `<` total direct seconds.
- setup seconds `<=` total direct seconds.
- If Scope D passes, the next step is Scope E: a controlled opt-in prototype
  that extracts one real TABI surface and unknown vector from a mini-protein
  solve and runs an offline JCP-channel comparison. Scope E must still avoid
  changing production defaults.

## Scope E - Offline JCP K1 Bridge Gate

Goal: compare an exact direct JCP channel against the scalar M2L engine where
the mathematics genuinely matches. Scope E starts with the K1 phi-from-phin
channel because it is the scalar kernel difference `G0 - Gk`.

Boundary:

- Scope E is still research-only.
- It does not call production GMRES or replace the TABI treecode operator.
- It does not pretend K2/K3/K4 are scalar M2L channels. Those derivative
  channels need a separate derivative-M2L design in Scope F.
- The K1 direct reference is:
  `sum_j area_j * (G0(r_ij) - Gk(r_ij)) * phin_j`.
- The M2L approximation is built from two scalar passes:
  Laplace with source `area * phin / (4*pi)` minus Yukawa with the same source.

Implementation checklist:

- [x] Add `JcpK1BridgeMetrics`.
- [x] Add deterministic `phin` construction from `PrototypePanel`.
- [x] Add `run_jcp_k1_bridge_probe(...)`.
- [x] Extend the C++ smoke to verify finite K1 bridge outputs, nonzero far work,
  and bounded relative error.
- [x] Add `test/pbe_fmm_m2l_scope_e_probe.ps1` and
  `TITAN_FMM_M2L_SCOPE_E_REPORT.txt`.

Scope E gate:

- K1 relative RMS `<= 5e-2`.
- far-pair work must be nonzero.
- M2L apply seconds `<` direct K1 seconds.
- setup seconds `<=` direct K1 seconds.
- If Scope E passes, Scope F may design derivative M2L for K2/K3/K4. If Scope E
  fails, do not move to derivative channels.

## Scope F - Derivative-Channel Bridge Gate

Goal: test K2/K3/K4 derivative-shaped channels without touching production
TABI. Scope F uses finite-difference dipole source/target probes around the
existing scalar M2L engine, then compares against direct analytic derivative
reference formulas from the JCP leaf path.

Boundary:

- Scope F is still research-only.
- It does not call production GMRES or replace the TABI treecode operator.
- It is not the final derivative M2L implementation. It is a bridge to decide
  whether derivative-channel work has enough numerical and performance ROI.
- K2 source-normal derivative is represented by plus/minus scalar sources
  shifted along the source normal.
- K3 target-normal derivative is represented by evaluating scalar potentials
  at target points shifted along the target normal.
- K4 mixed derivative uses both source-normal plus/minus sources and
  target-normal finite differences.

Implementation checklist:

- [x] Add `JcpDerivativeBridgeMetrics`.
- [x] Add direct analytic K2/K3/K4 reference evaluation.
- [x] Add finite-difference M2L approximation for K2/K3/K4.
- [x] Extend C++ smoke to verify finite output, nonzero far work, and bounded
  aggregate error.
- [x] Add `test/pbe_fmm_m2l_scope_f_probe.ps1` and
  `TITAN_FMM_M2L_SCOPE_F_REPORT.txt`.

Scope F gate:

- aggregate derivative-channel relative RMS `<= 1e-1`.
- max per-channel relative RMS `<= 2e-1`.
- far-pair work must be nonzero.
- M2L apply seconds `<` direct derivative seconds.
- setup seconds `<=` direct derivative seconds.
- If Scope F passes, Scope G may design a hidden offline real-surface JCP
  operator comparison. Production defaults still remain unchanged.

## Scope G - Production-Surface Offline JCP Gate

Goal: compare production direct JCP output with the K1-K4 FMM/M2L research
prototype on panels created by TITAN's production TABI surface builder from a
real PQR fixture.

Boundary:

- Scope G is offline and research-only.
- It uses `TitanPBE::Tabi::build_surface` and the production direct JCP fallback.
- It does not call GMRES, replace the production matvec, or add a public CLI.
- Direct comparisons use a deterministic bounded panel sample to keep the
  quadratic reference tractable.

Implementation checklist:

- [x] Add production direct JCP parity metrics to the prototype.
- [x] Add a C++ real-surface probe using production surface panels.
- [x] Add `test/pbe_fmm_m2l_scope_g_probe.ps1` and
  `TITAN_FMM_M2L_SCOPE_G_REPORT.txt`.
- [x] Gate `research-large` behind a passing `research-medium` result.

Scope G gate:

- production direct parity aggregate relative RMS `<= 1e-10`;
- K1 relative RMS `<= 5e-2`;
- derivative aggregate relative RMS `<= 1e-1`;
- derivative maximum relative RMS `<= 2e-1` for channels whose reference RMS
  is at least `0.1%` of the strongest derivative channel;
- far-pair work must be nonzero;
- combined M2L apply seconds `<` combined direct seconds.

## Scope H - Hidden Offline Single-Matvec Bridge

Goal: expose one research-only offline operator result that combines the Scope
G K1-K4 prototype channels into production-shaped `phi/phin` vectors, then
compare that single matvec against production direct JCP on bounded real-surface
samples.

Boundary:

- Scope H is still offline and hidden; it does not run GMRES.
- It does not add a public CLI or production backend.
- It reuses the Scope G production direct JCP reference and the same bounded
  `research-medium` / `research-large` fixtures.
- A pass means "continue to a controlled hidden TABI matvec experiment", not
  "promote FMM/M2L".

Implementation checklist:

- [x] Add a prototype API returning approximate `phi` and `phin` vectors from
  the current K1-K4 M2L bridge.
- [x] Add `test/titan_fmm_m2l_scope_h_probe.cpp` with a generated-surface smoke
  and real-fixture CSV mode.
- [x] Add `test/pbe_fmm_m2l_scope_h_probe.ps1` producing matrix, summary, and
  report files.
- [x] Gate `research-large` behind a passing `research-medium` result.

Scope H gate:

- production direct parity aggregate relative RMS `<= 1e-10`;
- offline single-matvec aggregate relative RMS `<= 1e-1`;
- offline single-matvec significant-channel maximum relative RMS `<= 2e-1`;
- far-pair work must be nonzero;
- prototype apply seconds `<` production direct seconds on the same sample.

## Scope I - Analytical Derivative Offline Operator

Goal: replace Scope H's finite-difference derivative bridge with a research-only
analytical derivative operator that computes K1/K2/K3/K4 in one tree traversal.

Boundary:

- Scope I is offline and hidden; it does not run GMRES.
- It does not add a public CLI or production backend.
- Near leaf pairs remain direct analytical K1/K2/K3/K4.
- Accepted far leaf pairs use cell moments for scalar `phin` and source-normal
  `phi` derivative terms.
- A pass means "plan a controlled hidden production matvec adapter"; a fail
  means do not wire FMM/M2L into the solver.

Scope I gate:

- production direct parity aggregate relative RMS `<= 1e-10`;
- analytical offline aggregate relative RMS `<= 1e-1`;
- analytical offline `max(phi_rel, phin_rel) <= 2e-1`;
- far-pair work must be nonzero;
- analytical apply seconds faster than production direct JCP, preferably
  `speedup > 1.2`.

## Current Status

Status: Scope A, Scope B point-cloud, Scope C panel-shaped, Scope D
JCP-channel-shaped, Scope E K1 bridge, Scope F derivative-channel bridge, Scope
G production-surface offline JCP, and Scope H hidden offline single-matvec
bridge are implemented. Scope H did not pass the performance gate. Scope I
analytical derivative offline operator is implemented and passed the offline
gate. Scope J hidden production-panel single-matvec adapter, Scope K guarded
GMRES probe, Scope L bounded end-to-end solve gate, Scope M hidden CLI gate,
and Scope N hidden physics/resolution/scale matrix are implemented and passed.
Scope O repeatability/cap-256 stress gate is implemented; its cap-128
repeatability tier passed, but the full gate is blocked by cap-256 timeout
classification and is not promotion-ready. Scope P Barnase hidden CLI gate is
implemented. The quick bounded Barnase cap-64 tier passed, but the full cap-128
repeatability tier did not pass because both direct and FMM rows remained in
`energy=failed;profile=failed` non-converged status at the current `1e-6`
iteration budget. The bounded FMM rows were finite, audit-clean, and fast, so
the result is a convergence-gate failure rather than a timeout or NaN failure.
Scope Q Barnase convergence sweep is implemented and passed: bounded Barnase
cap-128 reaches strict `1e-6` convergence when the GMRES iteration budget is
raised to `240`, with clean sampled audit and strong hidden FMM/M2L speedup.
Scope R fair-memory/cap-256 Barnase gate is implemented. Its original
memory-regression conclusion is stale after Scope S: hidden FMM/M2L no longer
builds the unused legacy cached `TreeWorkspace` interaction lists. The refreshed
full Scope R gate confirms cap-128 repeatability with process-level memory close
to direct TABI, while cap-256 remains not promotion-ready because direct and FMM
both hit the current `240` iteration/status gate. Scope T Barnase cap-256
convergence sweep is implemented; it classifies the cap-256 blocker as
`restart_limited`.

Implemented:

- Added `src/tabi/tabi_fmm_m2l_prototype.{h,cpp}` as an isolated research
  prototype. It is compiled with the TABI sources but is not called by the
  production PBE/TABI solver.
- Added deterministic synthetic panel generation, a balanced AABB tree,
  near/direct exact interaction lists, far low-order multipole/M2L-like
  evaluation, Laplace and Yukawa kernels, and deterministic `TitanExecutor`
  apply over target leaves.
- Added `test/titan_fmm_m2l_prototype_smoke.cpp` plus
  `test/fmm_m2l_scope_a_smoke.ps1`.
- Added `test/pbe_fmm_m2l_scope_a_probe.ps1` and generated
  `TITAN_FMM_M2L_SCOPE_A_REPORT.txt`.
- Added explicit point-cloud integration for Scope B:
  `PrototypePointCloud`, `run_point_cloud_probe(...)`, and deterministic
  `make_mini_protein_like_point_cloud(...)`.
- Added `test/pbe_fmm_m2l_scope_b_probe.ps1` and generated
  `TITAN_FMM_M2L_SCOPE_B_REPORT.txt`.
- Added TABI-panel-shaped scalar integration for Scope C:
  `PrototypePanel`, `point_cloud_from_panels(...)`, deterministic
  `make_tabi_panel_like_fixture(...)`, and `test/pbe_fmm_m2l_scope_c_probe.ps1`.
- Added JCP-channel-shaped research integration for Scope D:
  `JcpChannelPrototypeMetrics`, deterministic four-channel source construction,
  `run_jcp_channel_probe(...)`, and `test/pbe_fmm_m2l_scope_d_probe.ps1`.
- Added real JCP K1 scalar bridge for Scope E:
  `JcpK1BridgeMetrics`, direct K1 reference evaluation, Laplace-minus-Yukawa
  M2L approximation, and `test/pbe_fmm_m2l_scope_e_probe.ps1`.
- Added derivative-channel bridge for Scope F:
  `JcpDerivativeBridgeMetrics`, direct analytic K2/K3/K4 reference evaluation,
  finite-difference source/target normal scalar M2L probes, per-channel error
  diagnostics, and `test/pbe_fmm_m2l_scope_f_probe.ps1`.
- Added production-surface offline comparison for Scope G:
  `JcpDirectParityMetrics`, production `build_surface` integration, production
  direct JCP parity, bounded deterministic real-surface panel sampling, and
  `test/pbe_fmm_m2l_scope_g_probe.ps1`.
- Added Scope H hidden offline single-matvec bridge:
  `JcpOfflineOperatorResult`, `evaluate_jcp_offline_operator(...)`,
  production-shaped `phi/phin` output vectors, and
  `test/pbe_fmm_m2l_scope_h_probe.ps1`.
- Added Scope I analytical derivative offline operator:
  `JcpAnalyticOperatorResult`,
  `evaluate_jcp_offline_operator_analytic(...)`, direct near K1/K2/K3/K4,
  far cell moments for scalar `phin` and source-normal `phi`, and
  `test/pbe_fmm_m2l_scope_i_probe.ps1`.
- Added Scope J hidden production-panel single-matvec adapter:
  `-pbe-tabi-fmm-m2l=single-matvec-analytic` as a hidden option,
  `FmmM2lSingleMatvecStats`,
  `apply_jcp_fmm_m2l_single_matvec(...)`, production-sign output
  `lhs*unknown - integral_sum`, production-panel/unknown-vector conversion,
  hidden routing in `apply_jcp_treecode_operator`, and
  `test/pbe_fmm_m2l_scope_j_probe.ps1`.
- Added hidden CLI matrix/repeatability harnesses through Scope O, including
  `test/pbe_fmm_m2l_cli_helpers.ps1`,
  `test/pbe_fmm_m2l_scope_n_matrix.ps1`, and
  `test/pbe_fmm_m2l_scope_o_repeatability_probe.ps1`.
- Added Scope P bounded Barnase hidden CLI gate:
  `test/pbe_fmm_m2l_scope_p_barnase_gate.ps1`, Barnase fixture resolution via
  `test/input/barnase.pqr` or pyDelPhi fixture conversion, and
  `TITAN_FMM_M2L_SCOPE_P_BARNASE_REPORT.txt`.
- Added Scope Q bounded Barnase convergence sweep:
  `test/pbe_fmm_m2l_scope_q_barnase_convergence_sweep.ps1` and
  `TITAN_FMM_M2L_SCOPE_Q_BARNASE_CONVERGENCE_REPORT.txt`.
- Added Scope R bounded Barnase fair-memory/cap-256 gate:
  `test/pbe_fmm_m2l_scope_r_barnase_memory_cap256_gate.ps1`, process-level
  peak working-set capture in the shared CLI helper, and
  `TITAN_FMM_M2L_SCOPE_R_BARNASE_MEMORY_CAP256_REPORT.txt`.
- Added Scope S FMM/M2L memory fix:
  cached analytical FMM/M2L geometry context reuse, caller-owned
  `TitanExecutor` reuse, and hidden FMM/M2L solver routing that avoids building
  legacy cached tree interactions when no preconditioner requires them.
- Added Scope T bounded Barnase cap-256 convergence sweep:
  `test/pbe_fmm_m2l_scope_t_cap256_convergence_sweep.ps1` and
  `TITAN_FMM_M2L_SCOPE_T_CAP256_CONVERGENCE_REPORT.txt`.

Verification snapshot:

- `powershell -NoProfile -ExecutionPolicy Bypass -File test\fmm_m2l_scope_a_smoke.ps1`
  passed with `PASS titan_fmm_m2l_prototype_smoke`.
- Quick probe passed and wrote `test/tmp/fmm_m2l_scope_a/fmm_m2l_scope_a_matrix.csv`.
- Representative probe:

  ```powershell
  powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_a_probe.ps1 -PanelCounts 2048 -Orders 2 -Theta 0.6 -Threads 4
  ```

  Best row from `TITAN_FMM_M2L_SCOPE_A_REPORT.txt`:

  ```text
  panels=2048 kernel=yukawa order=2 theta=0.6 rel_rms_error=0.000111452663183 speedup=5.07641331915 decision=promote_to_scope_b
  ```
- Scope B representative probe:

  ```powershell
  powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_b_probe.ps1 -PanelCounts 2048 -Orders 2 -Theta 0.55 -Threads 4
  ```

  Best row from `TITAN_FMM_M2L_SCOPE_B_REPORT.txt`:

  ```text
  panels=2048 kernel=yukawa order=2 theta=0.55 rel_rms_error=0.000539123021504 speedup=3.29716804201 gate=pass
  ```
- Scope C representative probe:

  ```powershell
  powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_c_probe.ps1 -PanelCounts 2048 -Orders 2 -Theta 0.5 -Threads 4
  ```

  Best row from `TITAN_FMM_M2L_SCOPE_C_REPORT.txt`:

  ```text
  panels=2048 kernel=laplace order=2 theta=0.5 rel_rms_error=0.000643292468717 speedup=1.69286812112 gate=pass
  ```
- Scope D representative probe:

  ```powershell
  powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_d_probe.ps1 -PanelCounts 2048 -Orders 2 -Theta 0.8 -Threads 4
  ```

  Best row from `TITAN_FMM_M2L_SCOPE_D_REPORT.txt`:

  ```text
  panels=2048 kernel=yukawa order=2 theta=0.8 aggregate_rel_rms_error=0.00413466788901 max_channel_rel_rms_error=0.00591104044265 near_pairs=6094848 far_pairs=10682368 speedup=6.14330704128 gate=pass
  ```
- Scope E representative probe:

  ```powershell
  powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_e_probe.ps1 -PanelCounts 2048 -Orders 2 -Theta 0.7 -Threads 4
  ```

  Best row from `TITAN_FMM_M2L_SCOPE_E_REPORT.txt`:

  ```text
  panels=2048 order=2 theta=0.7 rel_rms_error=0.00621637754907 near_pairs=3784704 far_pairs=4603904 speedup=5.60205568647 gate=pass
  ```
- Scope F representative probe:

  ```powershell
  powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_f_probe.ps1
  ```

  Best row from `TITAN_FMM_M2L_SCOPE_F_REPORT.txt`:

  ```text
  panels=2048 order=4 theta=0.7 aggregate_rel_rms_error=0.0648357897519 max_channel_rel_rms_error=0.0648646161014 k2=0.0648646161014 k3=0.00587158026322 k4=0.0505188958911 near_pairs=26042368 far_pairs=41066496 speedup=3.42307506673 gate=pass
  ```
- Scope G production-surface probe:

  ```powershell
  powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_g_probe.ps1
  ```

  Representative rows from `TITAN_FMM_M2L_SCOPE_G_REPORT.txt`:

  ```text
  research-medium atoms=6063 surface_panels=55062 sampled_panels=2048 parity=2.40586437131e-15 k1=0.00274278083935 derivative_aggregate=0.0219441041137 speedup=1.74161167298 gate=pass
  research-large atoms=12752 surface_panels=113468 sampled_panels=2048 parity=3.23379276979e-15 k1=0.00550837970958 derivative_aggregate=0.000341432712451 speedup=1.7558248512 gate=pass
  ```
- Scope H hidden offline single-matvec bridge:

  ```powershell
  powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_h_probe.ps1 -Repeats 5
  ```

  Representative row from `TITAN_FMM_M2L_SCOPE_H_REPORT.txt`:

  ```text
  research-medium atoms=6063 surface_panels=55062 sampled_panels=2048 parity=2.40586437131e-15 offline_aggregate=0.0218073842146 offline_phi=0.0218061086095 offline_phin=0.0298849275179 near_pairs=52838400 far_pairs=22659072 speedup=0.675245335009 gate=no_speedup
  ```
- Scope I analytical derivative offline operator:

  ```powershell
  powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_i_probe.ps1 -Repeats 5
  ```

  Representative rows from `TITAN_FMM_M2L_SCOPE_I_REPORT.txt`:

  ```text
  research-medium atoms=6063 surface_panels=55062 sampled_panels=2048 parity=2.40586437131e-15 offline_aggregate=0.0002820320478 offline_phi=0.000273385086965 offline_phin=0.00601041008738 near_pairs=3272704 far_pairs=921600 speedup=3.55950049282 gate=pass
  research-large atoms=12752 surface_panels=113468 sampled_panels=2048 parity=3.23379276979e-15 offline_aggregate=0.0014179982537 offline_phi=0.00140552481985 offline_phin=0.00925772972043 near_pairs=2953216 far_pairs=1241088 speedup=4.41025749568 gate=pass
  ```
- Scope J hidden production-panel single-matvec adapter:

  ```powershell
  powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_j_probe.ps1 -Repeats 5
  ```

  Representative rows from `TITAN_FMM_M2L_SCOPE_J_REPORT.txt`:

  ```text
  research-medium atoms=6063 surface_panels=55062 sampled_panels=2048 aggregate=0.0002820320478 phi=0.000273385086965 phin=0.00601041008738 near_pairs=3272704 far_pairs=921600 setup_ms=16.7886 apply_ms=47.3712 speedup=4.83624860675 gate=pass
  research-large atoms=12752 surface_panels=113468 sampled_panels=2048 aggregate=0.0014179982537 phi=0.00140552481985 phin=0.00925772972043 near_pairs=2953216 far_pairs=1241088 setup_ms=12.0007 apply_ms=27.3461 speedup=6.25971893616 gate=pass
  ```

Decision:

- Scope A passes its first scalar synthetic gate on the 2048-panel
  representative probe.
- Scope B passes the explicit point-cloud mini-protein-like gate on the
  2048-panel representative probe.
- Scope C passes the TABI-panel-shaped scalar geometry/work-count gate on the
  2048-panel representative probe.
- Scope D passes the JCP-channel-shaped work-count gate on the 2048-panel
  representative probe, with nonzero far-pair work.
- Scope E passes the real JCP K1 scalar bridge gate on the 2048-panel
  representative probe, with nonzero far-pair work.
- Scope F passes the finite-difference derivative-channel bridge gate for
  K2/K3/K4 on the 2048-panel representative probe, with nonzero far-pair work.
- Scope G passes production direct parity plus K1-K4 M2L accuracy/performance
  gates on bounded 2048-panel samples from both real research fixtures.
- Scope H confirms that the combined offline `phi/phin` operator is accurate
  enough on `research-medium` (`~2.18%` aggregate relative RMS), but it is
  slower than production direct JCP on the bounded sample (`speedup=0.675`) and
  therefore fails the integration gate.
- Scope I fixes the Scope H pass-count bottleneck with analytical derivative
  kernels. Both bounded real-surface fixtures pass with aggregate relative RMS
  below `0.2%` and speedup above `3.5x`.
- Scope J wires the analytical derivative operator behind a hidden
  production-panel single-matvec adapter. It preserves production sign
  convention, consumes the real coupled JCP unknown vector, passes the
  bounded real-fixture accuracy gate (`<=1e-2` aggregate relative RMS), and
  exceeds the `1.2x` speed gate on `research-medium`.
- Scope K runs a guarded GMRES probe on bounded production-surface panel
  samples with `b = A_direct * deterministic_x_true`. Both real fixtures pass
  at `1024` sampled panels: `research-medium adapter_residual=5.4590676773e-07`,
  `direct_true_residual=2.96990273845e-05`, `speedup=3.20555382964`; and
  `research-large adapter_residual=3.85417787791e-07`,
  `direct_true_residual=4.21884160643e-05`, `speedup=3.23620232888`.
- Scope L runs the hidden FMM/M2L adapter inside bounded end-to-end
  `solve_linear_pb` TABI solves and compares convergence, residual, energy,
  and matvec timing against the direct JCP baseline. Both real fixtures passed
  at `128` deterministically sampled atoms with `surface_point_count=8`:
  `research-medium sampled_atoms=128 panels=6133 fmm_true_residual=6.23991642839e-11 energy_delta_abs=7.87664480659 energy_delta_rel=0.000943543177951 speedup=7.60893491996 gate=pass`;
  `research-large sampled_atoms=128 panels=6132 fmm_true_residual=5.4563596887e-15 energy_delta_abs=4.93958074564 energy_delta_rel=0.000774204934844 speedup=9.12711951306 gate=pass`.
- Scope M runs the hidden FMM/M2L path through real `titan.exe -a arg.in`
  CLI cases using bounded sampled PQR fixtures. Both direct and FMM CLI runs
  exited cleanly, wrote `pbe_energy.csv`, `pbe_timing.csv`, and
  `pbe_tabi_profile.csv`, and passed the energy/residual/speed gates:
  `research-medium sampled_atoms=128 panels=6133 fmm_true_residual=6.23992e-11 energy_delta_abs=7.87664399999994 energy_delta_rel=0.000943543081367363 speedup=5.88432380093087 gate=pass`;
  `research-large sampled_atoms=128 panels=6132 fmm_true_residual=5.45636e-15 energy_delta_abs=4.93958100000054 energy_delta_rel=0.000774204974668301 speedup=12.4719178640643 gate=pass`.
- Scope N runs a broader hidden CLI matrix across scale, fixture, salt,
  dielectric, Poisson salt-zero, and surface-resolution conditions using
  sampled direct-operator audit. The full matrix passed `8/8` required rows
  with `opening_theta=0.6`; `surface=16` was used for the resolution row
  because `surface=12` maps to the same mesh bucket as `surface=8` in the
  current surface builder. Representative rows from
  `TITAN_FMM_M2L_SCOPE_N_REPORT.txt`:
  `medium_cap128_baseline panels=6133 energy_delta_abs=1.9216379999998 audit_rel=0.00323097 speedup=5.21319583986574 gate=pass`;
  `large_cap128_baseline panels=6132 energy_delta_abs=0.195887999999286 audit_rel=0.00261365 speedup=5.93470020675396 gate=pass`;
  `alt_cap128_baseline panels=6120 energy_delta_abs=5.21474400000079 audit_rel=0.00499974 speedup=5.71991452350604 gate=pass`;
  `medium_cap64_salt0 panels=3065 energy_delta_abs=0.22892899999988 audit_rel=0.00113586 speedup=4.99505955726381 gate=pass`;
  `medium_cap64_surface16 panels=5110 energy_delta_abs=1.32726600000024 audit_rel=0.00417215 speedup=6.5734317938785 gate=pass`.
- Scope O implements repeatability plus cap-256 stress gating. The repeatability
  tier passed on all three cap-128 fixtures with deterministic energies across
  three FMM repeats and sampled audit below `1e-2`: `research-medium
  audit_max=0.00323097 speedup_median=7.15721267875626 speedup_min=5.14596932441632`;
  `research-large audit_max=0.00261365 speedup_median=5.72520286151872
  speedup_min=4.9350600868128`; `research-alt audit_max=0.00499974
  speedup_median=8.92799795952675 speedup_min=7.53845927941774`.
  The full Scope O gate did not pass because the cap-256 stress tier is
  `timeout_not_verified`: `research-medium atom_cap=256 panels=12250
  fmm_true_residual=2.28745e-07 audit_rel=0.00502187 fmm_energy=-16086.4942`,
  but the direct baseline did not produce finite energy before the timeout, so
  energy delta and speedup cannot be verified; `research-large` cap-256 was
  correctly skipped. This blocks promotion planning without classifying the
  hidden FMM solver as scientifically failed.
- Scope P fills the Barnase gap for the hidden FMM/M2L path. Legacy TABI
  treecode Barnase rows were non-promotable, but Scope P tests only bounded
  sampled Barnase through the hidden analytical FMM/M2L CLI path. The quick
  cap-64 smoke passed with `panels=3053`, `energy_delta_abs=0.526157`,
  `audit_rel=0.00143331`, `fmm_residual=1.50766e-11`, and
  `speedup=4.82697815048099`. The full cap-128 repeatability tier ran three
  repeats at `panels=5802`; all FMM repeats were finite and deterministic with
  `energy_delta_abs=2.354467`, `audit_rel=0.00268276`, and speedups
  `7.37443711383391`, `6.01715806698311`, and `9.17630425009285`. The full
  gate still failed as `convergence_fail;row_fail` because the direct and FMM
  CLI rows did not report `status=ok` under the current `1e-6` / `120`
  iteration gate (`direct_true_residual=0.000916239`,
  `fmm_true_residual=0.00080777`). Cap-256 Barnase stress was correctly
  skipped because cap-128 did not pass.
- Scope Q diagnoses the Scope P cap-128 convergence failure. The baseline
  `max_iter=120`, `restart=24`, `1e-6` row reproduced the non-converged status
  (`direct_true_residual=0.000916239`, `fmm_true_residual=0.00080777`) while
  remaining audit-clean and fast (`audit_rel=0.00268276`, `speedup=7.199`).
  Raising `max_iter` to `240` at the same restart produced strict convergence:
  `direct_true_residual=9.31161e-7`, `fmm_true_residual=6.26312e-7`,
  `energy_delta_abs=1.716306`, `audit_rel=0.00267975`, `speedup=7.337`.
  Restart `48` and `64` improved the same strict row further:
  `restart=48` converged in `144` iterations with
  `fmm_true_residual=3.60205e-8`, `speedup=5.802`; `restart=64` converged in
  `128` iterations with `fmm_true_residual=1.27442e-8`, `speedup=7.299`.
  The practical `1e-4` row also converged in `96` iterations. Scope Q final
  class is `strict_converged`.
- Scope S fixed the Scope R process-memory regression by reusing the analytical
  FMM/M2L context and avoiding unused legacy cached tree-interaction lists in
  the hidden FMM path. The refreshed full Scope R run uses the same Scope Q
  preset (`max_iter=240`, `restart=64`, `1e-6`). Cap-128 repeated three times
  with deterministic energy and strict convergence: FMM
  `true_residual=1.27442e-8`, `energy_delta_abs=1.718162`,
  `audit_rel=0.00267974`, speedup range `9.43x..11.66x`, and zero FMM energy
  spread. Process peak working set is now close to direct TABI: grid rows peak
  around `208 MB`, direct TABI around `16.43..16.47 MB`, and FMM/M2L around
  `18.65..18.67 MB` (`~1.13x` direct, `~0.09x` grid). Cap-256 no longer shows
  memory blow-up (`27 MB` FMM peak), but direct/FMM both reached the `240`
  iteration limit without `status=ok` (`direct_true_residual=0.0140401`,
  `fmm_true_residual=0.0131365`). FMM remained audit-clean
  (`audit_rel=0.0062981`) and fast (`speedup=17.17x`). Refreshed Scope R final
  class is `cap128_ready_cap256_not_promotion_ready`.
- Scope T diagnoses the refreshed Scope R cap-256 convergence/status blocker
  with an FMM-first GMRES budget sweep plus direct verification for the baseline
  and best candidate. Baseline cap-256 `240/r64/1e-6` reproduced the Scope R
  non-converged state (`direct_true_residual=0.0140401`,
  `fmm_true_residual=0.0131365`, `audit_rel=0.0062981`,
  `speedup=17.22x`). Increasing iterations at restart `64` improved FMM
  residual gradually (`360/r64=0.00919207`, `480/r64=0.00640712`), while larger
  restart produced the main gain (`480/r96=0.00111163`,
  `480/r128=0.00147099`). The best row `720/r128/1e-6` reached
  `fmm_true_residual=0.000185057` with direct verification
  `direct_true_residual=0.000194446`, `energy_delta_abs=7.33137`,
  `audit_rel=0.00685148`, `speedup=19.92x`, and no timeout. Scope T final
  class is `restart_limited`.
- Scope U adds a hidden FMM-native Block-Jacobi preconditioner
  (`-pbe-tabi-precond=fmm-block-jacobi`) built from
  `FmmM2lAnalyticContext` leaf/block geometry instead of legacy
  `TreeWorkspace` interaction lists. The full Barnase gate keeps cap-128
  strict and audit-clean: `true_residual=1.41048e-8`,
  `energy_delta_abs=1.718147`, `audit_rel=0.00267974`, `speedup=10.07x`,
  `preconditioner_valid_blocks=256`, and peak working set `23.73 MB`.
  Cap-256 remains hidden/research-only but shows useful convergence movement:
  `240/r64` improves residual from diagonal `0.0131365` to block
  `0.0121808`, and the longer `480/r128` row improves from diagonal
  `0.00147099` to block `0.000852217` with audit still clean
  (`audit_rel=0.00686342`) and no memory regression (`44.98 MB` peak). The
  direct verifier timed out for the longer cap-256 rows, so cap-256 is still
  not promotion-ready.
- Scope V opens the FMM-native block cap to `128` while keeping legacy
  `block`/`block-budgeted` capped at `32`, then runs a Barnase cap-256 block-cap
  sweep. Quick cap-128 parser/CLI probes passed: `fmm-block-jacobi` accepts cap
  `128`, legacy `block-budgeted` rejects cap `128`, and FMM rejects cap `129`.
  Full cap-256 shows the cap-only hypothesis is a no-win: at `480/r128/1e-6`,
  cap `32` residual is `0.000852217`, while cap `64/128` are both
  `0.000958293`; at `720/r128/1e-6`, cap `32` is `0.000213074`, while
  cap `64/128` are both `0.000230641`. Cap `64/128` remain audit-clean
  (`~0.006853`) and memory-light (`~48.6 MB` peak), but do not improve
  convergence. Cap `128` matches cap `64` because the current FMM tree leaf size
  is `64`, so the block cannot grow beyond the leaf-local span.
- Scope W adds a hidden FMM-native Restricted Additive-Schwarz preconditioner
  (`-pbe-tabi-precond=fmm-overlap-schwarz`) with deterministic halo panels from
  FMM `near_sources` and restricted writeback to owned panels. The quick
  Barnase cap-128 gate validates parser behavior, generated `arg.in` syntax,
  overlap telemetry, and strict convergence. With owned cap `32`, overlap `0`
  gives `true_residual=1.41048e-8`, `audit_rel=0.00269123`,
  `speedup=11.89x`, peak memory `23.83 MB`; overlap `32` improves residual to
  `3.06488e-9` with `speedup=13.11x` and peak memory `44.61 MB`; overlap `64`
  gives `5.90696e-9`, `speedup=11.91x`, and peak memory `80.09 MB`. This
  proves the cross-leaf restricted-overlap machinery works on cap-128. The full
  Barnase cap-256 gate is a no-win: at `480/r128/1e-6`, baseline
  `fmm-block-jacobi` cap `32` reached `0.000852217`, while overlap `32/64`
  worsened to `0.00109664` and `0.00116147`; at `720/r128/1e-6`, baseline
  reached `0.000213074`, while overlap `32/64` worsened to `0.000239418` and
  `0.000344171`. Overlap `128` exceeded the memory budget and fell back to
  diagonal (`diagonal_fmm_overlap_budget_cap`), so it is not a valid Schwarz
  win.
- Production default is unchanged. `-pbe=grid` remains the production default,
  and TABI/FMM work remains explicit research-only.
- Scope L/M use `128` atoms as the safe default bounded cap after a measured
  `128`-atom direct baseline already produced about `6133` panels and a costly
  direct solve. The original `512`-atom idea remains a manual stress override,
  not the default smoke size.
- Scope X adds the hidden sparse near-field ILU(0) diagnostic path:
  `-pbe-tabi-precond=fmm-nearfield-ilu0`. It exports a deterministic FMM
  near-source graph, assembles a sparse coupled `2x2` BSR near-field matrix,
  factors it with ILU(0), and applies it through the existing GMRES custom
  preconditioner hook. Hidden controls are
  `-pbe-tabi-precond-ilu-neighbor-cap=32..1024` and
  `-pbe-tabi-precond-ilu-drop-tol=0..1e-2`; default drop tolerance is `0`.
  Sparse telemetry records nnz blocks, average/max degree, cap/drop controls,
  setup/apply timing, and explicit budget/factor fallback labels.
- Scope X.1 completes the Barnase cap-256 full gate. The block baseline at
  `720/r128/1e-6` remained non-strict (`true_residual=0.000213074`), while
  sparse ILU(0) converged all screened caps to machine-level residual in `128`
  iterations. The cap-64 row is the best practical candidate:
  `true_residual=5.64662e-15`, `audit_rel=0.00684539`, `nnz_blocks=644416`,
  peak working set `61.58 MB`, and speedup `20.61x` versus the direct verifier.
  Direct TABI at `720/r128/1e-6` completed but remained non-converged
  (`direct_true_residual=0.000194446`), so the final Scope X.1 decision is
  `strict_fmm_pending_direct`. Energy delta for cap-64 ILU was
  `4.55179 kJ/mol`, inside the planned tolerance, but direct non-convergence
  blocks `strict_verified`.
- Scope Y implements the controlled TABI treecode default. `-pbe=grid` remains
  the production PBE default, but `-pbe=tabi -pbe-tabi-sum=treecode` now
  resolves `-pbe-tabi-fmm-m2l=auto` to `single-matvec-analytic` and resolves
  `-pbe-tabi-precond=auto` to `fmm-nearfield-ilu0` with ILU neighbor cap `64`.
  Explicit overrides remain intact: `-pbe-tabi-fmm-m2l=off` keeps legacy
  treecode behavior and `-pbe-tabi-precond=diagonal` keeps diagonal
  preconditioning. Scope Y quick smoke verified the default treecode row
  matches the explicit FMM+ILU0 reference, the legacy escape reports
  `fmm_m2l_mode=off`, and direct TABI does not auto-enable FMM/M2L.
- Scope Z5A adds the preconditioner money-plot harness:
  `test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1`. This harness
  exports matrix, residual trace, summary, report, and SVG plot artifacts for
  diagonal/block/overlap/ILU preconditioner comparisons. The quick Barnase
  cap-128 run produced `harness_ready_but_no_contrast`: Block-Jacobi reached
  `true_residual=3.89211e-14`, sparse ILU(0) reached
  `true_residual=2.85424e-15`, and both rows were audit-clean. Z5A is
  harness-first.
- Scope Z5B completes the full Barnase cap-256 money-plot gate. The required
  rows were run first: Block-Diagonal / FMM Block-Jacobi stalled at
  `true_residual=2.13074e-4` after 720 iterations, while sparse near-field
  ILU(0) reached `true_residual=5.64662e-15` after 128 iterations with
  `tree_vs_direct_matvec_rel=0.00684539`. Peak process memory stayed modest
  (`45.2 MB` for Block-Jacobi, `61.6 MB` for ILU0). The decision is
  `money_plot_pass`, and the SVG/CSV artifacts are under
  `test/tmp/fmm_m2l_scope_z5/`.
- Scope Z5C adds the multi-fixture repeatability/generalization matrix:
  `test/pbe_fmm_m2l_scope_z5c_multi_fixture_money_plot.ps1`. The full run
  repeated Barnase cap-256 three times and reproduced the money plot each time:
  Block-Jacobi stalled at `2.13074e-4`, sparse ILU(0) reached
  `5.64662e-15`, and the ILU-vs-block residual ratio was `3.77e10`.
  However, sampled `research-medium`, `research-alt`, and gated
  `research-large` cap-256 rows all had Block-Jacobi and sparse ILU(0) strict
  at 128 iterations; the decision is `barnase_only_money_plot`. The paper claim
  should therefore be framed as dense/high-coupling Barnase-style mesh evidence,
  not yet a broad protein-fixture generalization.
- Scope Z1A adds the FAGBI paper-parity validation harness:
  `test/pbe_fagbi_scope_z1a_parity_matrix.ps1`. The full gate writes matrix,
  decision, summary, report, and SVG artifacts under `test/tmp/fagbi_scope_z1a/`.
  Its first full result is `z1a_accuracy_gap`: Barnase, research-medium, and
  research-alt cap-128 FMM/ILU0 rows passed direct-energy and sampled-audit
  checks, but the analytic Born and off-center Kirkwood sphere rows at surface
  point counts 8, 12, and 24 exceeded the oracle thresholds. This points the
  next scientific work toward Z2/Z3 accuracy parity rather than another
  preconditioner slice.
- Scope Z1B resolves the first Z1A attribution question with a direct-only
  `20/80/320/1280` icosphere matrix across centered Born, off-center Kirkwood,
  salt `0.0/0.15 M`, and the three near modes. The measured decision is
  `z1b_mesh_dominant`: all 48 rows were finite and fallback-free, but no
  `0.15 M` group reached the planned `2x` non-centroid improvement. At zero
  salt all three near modes were energetically identical. Refinement reduced
  oracle error consistently, while the non-centroid gain decayed at finer
  meshes. Do not treat FMM Duffy integration as the automatic next fix; first
  extend the surface-resolution/discretization study, then integrate only the
  near policy with demonstrated marginal value.
- Scope Z1C adds corrected surface-quality telemetry and a direct-only
  geometry-convergence matrix. The full 32-row gate completed without fallback
  and classified `resolution_limited`: all eight `1280 -> approximately 5120`
  transitions improved energy error by at least `1.5x`, while none showed an
  operator-discretization plateau. Finest median oracle error was `0.6419%` for
  icosphere and `0.6755%` for latlong, so `preferred_mesh=none` and
  `mesh_quality_tie`. The next accuracy scope is Z1D refinement beyond the
  current approximately 5120-panel range, not an immediate FMM-Duffy or
  collocation rewrite.
- Scope Z1D adds the high-resolution direct TABI gate:
  `test/pbe_fagbi_scope_z1d_high_resolution_surface.ps1`, with the hidden
  `-pbe-tabi-highres-surface=true` guard. The full run initially reached
  `z1d_runtime_limited`: 7/8 rows completed, and the completed 20480-panel rows
  were below 0.5% oracle error (`0.3705%`, `0.2682%`, `0.3874%`). Scope Z1D.1
  added `test/pbe_fagbi_scope_z1d_targeted_missing_row.ps1` and verified the
  missing Kirkwood `0.15 M` row at 20480 panels: `0.1374%` oracle error,
  `1750.93 s`, final residual `6.12512e-15`, and zero near fallbacks. Treat
  Z1D as high-resolution converged after targeted verification; the remaining
  constraint is direct-runtime cost, not missing surface-convergence evidence.
- Next step is broader controlled-treecode validation, not grid promotion.
  Scope Q removed the Barnase cap-128 convergence blocker, Scope S removed the
  hidden FMM/M2L memory regression, Scope U/V/W eliminated dense local
  preconditioner variants, and Scope X.1 selected sparse near-field ILU as the
  treecode default preconditioner. Keep `-pbe=grid` as production default until
  repeatability, broader fixtures, and a robust verification oracle are all
  demonstrated.
