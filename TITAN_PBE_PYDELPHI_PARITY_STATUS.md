# TITAN PBE / pyDelPhi Parity Status

Date: 2026-05-07

This note records the current PBE parity state so later phases can continue from
evidence instead of rediscovering the same mismatch chain.

## Runtime Independence

TITAN PBE is intended to run independently. Production TITAN does not require
DelPhi or pyDelPhi at runtime, and should not read oracle CSV files for normal
PBE, PBSA, Gaussian, RPB, focus, or probe calculations.

TITAN PBE is a self-contained implementation validated against black-box
DelPhi/pyDelPhi oracle outputs. pyDelPhi/DelPhi is used only as a local test
oracle under `.codex_tmp` and `test/tmp` to export traces, maps, component
energies, and timing references. Those oracle files are for regression/debug
comparison only.

Clean-room boundary: `PBE_CLEAN_ROOM_POLICY.md` is now the source of truth.
The guarantee is forward-looking from 2026-05-07 because older parity work did
inspect a local pyDelPhi checkout. Future production PBE changes must not read,
copy, translate, import, link, or call DelPhi/pyDelPhi source/runtime; only
black-box oracle outputs, timing rows, hashes, stage CSVs, internal specs, and
public documentation may be used.

## Source Layout

- The public API is `src/titan_pbe.h`; it now exposes `FrameResult` and a PIMPL
  `TitanPBE::Analyzer` wrapper so engine call sites keep including the same
  header.
- The implementation lives under `src/tabi/`. Phase 3 module splitting has
  moved map/cube I/O to `audit_io.cpp`, PBSA aggregation/writers to `pbsa.cpp`,
  focus/probe helpers to `focus_probe.cpp`, and Gaussian/RPB audit helpers plus
  RPB source CSV/top-cell diagnostics to `gaussian_rpb.cpp`. Grid/radius/surface
  setup now lives in `grid_surface.cpp`, and charge assignment/operator setup
  plus oracle replay helpers now live in `charge_operator.cpp`. Linear
  boundary/SOR/focus routing now lives in `linear_solver.cpp`, and linear
  V-cycle hierarchy now lives in `multigrid.cpp`. The remaining nonlinear,
  energy, reaction, and frame orchestration bodies are still mostly in
  `analyzer_impl.cpp` and should be moved one or two subsystems at a time with
  native build/smoke gates after each slice.
- `titan_surface.*` is the shared Shrake-Rupley surface workspace for GBSA,
  PBSA, and TABI. `titan_sasa.*` is now a compatibility wrapper for legacy
  per-atom SASA callers, while TABI consumes exposed panels with position,
  normal, area, and atom index. It now has an internal
  `auto|legacy|reordered_padded` acceleration selector with CSR cell lists, Morton atom
  ordering, and padded neighbor slots; public GBSA/PBSA/TABI outputs remain
  ordered by original atom index. TABI treecode acceleration is internal C++23
  code using Morton-ordered linear octree arrays plus a direct fallback.
- `Makefile` builds `src/tabi/*.cpp`; `compile.bat` remains unchanged because
  it delegates to the Makefile.

## What Is Currently Locked Or Usable

- Linear PBE, two-dielectric, CPU double path is usable.
- DelPhi-style VDW parity path has locked the main stage chain used during the
  long barnase debug: midpoint epsilon, boundary condition, charge grid,
  RHS/screen/denominator, SOR update, and water field samples.
- ISCM / DelPhi-corrected reaction-energy path is good enough for current
  parity work; residual native ISCM position mismatch is accepted as diagnostic,
  not a blocker for the next stages.
- Gaussian and gaussian-cutoff paths have strict smoke coverage against the
  current pyDelPhi oracle fixtures.
- Gaussian charge, RPB final-energy/operator, and RPB analytic source-audit
  paths have strict smoke coverage against the current oracle fixtures. The
  RPB source lock depends on matching pyDelPhi `COORDLIMITS` test-grid
  auto-centering, i.e. the radius-expanded boundary-box center, not the raw
  coordinate-extrema center.
- Nonlinear NWT barnase final convergence and energy now pass the current
  strict acceptance rows for salt `0.15 M` and `2.0 M`:
  - `barnase_nwt`: `nl_status=ok`, `nl_residual=1.231546642e-07`,
    `delta=4.36591291189961 kJ/mol`, tolerance `254.751661045595 kJ/mol`.
  - `barnase_nwt_salt2`: `nl_status=ok`, `nl_residual=2.436165942e-08`,
    `delta=4.66465514742413 kJ/mol`, tolerance `254.955489057371 kJ/mol`.
  - The NWT progress comparator now passes the strict block curve gate for the
    current oracle: `barnase_nwt` and `barnase_nwt_salt2` both match block
    counts and final iterations, with max block deltas near roundoff.
- DelPhi focus/probe strict smoke currently passes:
  - TITAN `E_rxn_kJmol = -6.96234747329947`
  - pyDelPhi `E_rxn_kJmol = -6.96695999427574`
  - delta `0.00461252097626907 kJ/mol`
  - tolerance `0.348347999713787 kJ/mol`
- PBSA smoke continues to pass with the current PBE backend changes.

## What Should Not Yet Be Claimed As Fully Complete

Do not claim TITAN is a complete pyDelPhi replacement yet. The current state is
best described as: important pyDelPhi-compatible PBE paths are implemented and
tested, but full feature parity still needs more gates.

Remaining or partially locked areas:

- Nonlinear PB / DelPhi NWT final convergence and current block-by-block
  progress gates pass, but keep broader NWT parameter coverage before claiming
  all nonlinear pyDelPhi paths.
- Nonlinear FAS / true nonlinear multigrid remains report-only unless a later
  phase explicitly locks it.
- Gaussian/RPB is covered by current fixtures, but should be expanded across
  more pyDelPhi parameter combinations before claiming broad parity.
- Surface/radius conventions beyond the locked DelPhi VDW/Gaussian paths need
  wider fixture coverage.
- Focus/probe now passes the current smoke, but broader residue/site/group
  semantics still need more fixtures.
- Full DelPhi workflow parity, map formats, and uncommon pyDelPhi options are
  not guaranteed.

## Snapshot Stability Diagnostics

Snapshot-to-snapshot energy jumps are now diagnosable with default-off
stability tooling. Atom coordinates exactly on, or very near, integer or
half-grid planes are plausible jump triggers because charge assignment,
surface/midpoint epsilon, ion accessibility, and ISCM boundary sets are all
discrete grid decisions. The new audit labels the first suspicious stage instead
of smoothing or changing default energy behavior.

- `-pbe-stability-audit=off|summary|detail` writes
  `pbe_stability_audit.csv` with per-frame grid center/origin/spacing, energy
  delta, min atom-to-integer/half-grid-plane distance, near-plane atom counts at
  `1e-6`, `1e-4`, and `1e-3` grid units, stage flip candidates,
  surface/ion/ISCM count deltas, `suspect_stage`, anchor status/margin,
  hysteresis counters, and native-vs-stable/ensemble energy diagnostics.
- Detail mode also writes `pbe_stability_atom_samples.csv` with the nearest
  atom/residue/site samples, fractional grid coordinates, distance to the
  nearest integer/half-grid plane, and a stage hint.
- `-pbe-stability-mode=off|diagnose|anchor|hysteresis|ensemble` is hidden and
  default-off. `diagnose` only forces the summary audit. `anchor` promotes
  auto/per-frame grids to first-frame anchoring for trajectory stability.
  `hysteresis` additionally keeps previous solvent medium decisions for cells
  that remain inside a narrow atom-surface release band. `ensemble` runs a
  deterministic `+-1e-7 A` coordinate-jitter probe and reports mean/std in the
  audit without replacing `pbe_energy.csv`.
- `-pbe-grid-anchor=auto|per-frame|first-frame|fixed-center` is a hidden
  stability/bisect knob. `auto` and `per-frame` keep current behavior.
  `first-frame` reuses the first grid shape for matching scale/margin/focus
  grid keys inside a trajectory; `fixed-center` anchors generated grids to
  `-pbe-test-center` or `-pbe-focus-center`.
- `test/pbe_snapshot_stability_scan.ps1` writes
  `test/tmp/pbe_snapshot_stability_scan/pbe_snapshot_stability_summary.csv`.
  Without user snapshots it generates synthetic one-atom cases on the grid
  plane, offset by `1e-7 A` and `1e-5 A`, plus a half-grid control. The current
  synthetic gate passes and verifies near-plane counts without changing default
  energies.
- `test/pbe_grid_boundary_pydelphi_compare.ps1` writes
  `test/tmp/pbe_grid_boundary_pydelphi_compare/pbe_grid_boundary_pydelphi_compare.csv`.
  It keeps TITAN as the strict side of the synthetic scan and records pyDelPhi
  report-only rows when the local `pydelphi_env` is available. pyDelPhi still
  behaves as a deterministic grid discretizer here; no smoothing is assumed for
  parity paths.

## Current Performance Snapshot

Timing data below is from the latest full `test/pbe_timing_compare.ps1`
median snapshot. The timing harness now defaults to `RepeatCount=2`, writes raw
repeats plus median rows, and enforces the comparable-row speed gate on the
median. The script supports comma-separated `-Fixtures`; these are wall-clock
snapshots on the current machine, not portable benchmarks.

The hidden `-pbe-timing=off|summary|detail` flag is default-off. When enabled,
TITAN writes `pbe_timing.csv` with buckets for `setup`, `surface`, `epsilon`,
`boundary`, `solver`, `energy`, `map_io`, `rpb_source`, `focus_setup`, `probe`,
and `total`, plus atom count, configured thread count, workflow scope, iteration,
block, and status fields. Detail columns now also record focus/probe
sub-buckets plus NWT sub-buckets (`nwt_half_update`, `nwt_coeff_pack`,
`nwt_homo_update`, `nwt_hetero_update`, `nwt_reduce`,
`nwt_operator_residual`, `nwt_source_update`, `nwt_progress_io`,
`nwt_homo_stencil`, `nwt_homo_salt_math`, `nwt_homo_write`,
`nwt_thread_wait`, `nwt_block_overhead`), NWT homogeneous-kernel counters
(`nwt_homo_cells`, `nwt_homo_linear_cells`, `nwt_homo_salt_cells`,
`nwt_homo_updates`, `nwt_dispatch_calls`, `nwt_homo_ns_per_update`), and
Gaussian epsilon sub-buckets (`epsilon_center_density`, `epsilon_midpoint_density`,
`epsilon_ion_exclusion`, `epsilon_finalize`). Production output is unchanged
unless the hidden flag is explicitly enabled.

Oracle/timing scripts append black-box provenance to
`test/pbe_oracle_manifest.csv` with fixture, command class, input/output hashes,
Python executable, pyDelPhi root/commit when available, thread count, status,
and output path. This manifest is a test artifact only; production TITAN must
not read it.

Hidden NWT controls added in this pass:

- `-pbe-nwt-progress=auto|off|strict`: `auto` keeps the strict progress CSV
  whenever parity dump is enabled; timing rows can set `off` to avoid charging
  TITAN for per-block source/residual/progress diagnostics that pyDelPhi does
  not include in its solver timing.
- `-pbe-nwt-kernel=auto|legacy|delphi-half`: `auto` now uses the packed
  DelPhi half-grid workspace for the pyDelPhi/NWT fair-timing path when
  `pbe_nwt_progress=off`; strict progress rows still use the legacy full-grid
  trace path unless `delphi-half` is forced explicitly.
- `-pbe-nwt-homo-kernel=auto|packed-list|row-span|avx2|salt-tiled`: hidden
  bisect knob for the homogeneous half-grid NWT update. `auto` now selects the
  `salt-tiled` path after strict/kernel-mode gates pass. This path keeps the
  packed-list and row-span kernels available for bisection, tiles homogeneous
  salt spans, hoists color-specific offsets/pointers, and uses a double
  precision polynomial `sinh`/`cosh` fast branch for small clipped potentials.
  The code is compiled without fast-math and remains guarded by strict NWT
  progress and kernel-mode smokes.
- `-pbe-linear-kernel=auto|legacy|half-packed`: hidden bisect knob for linear
  SOR. The half-packed path passes focus strict gates but is slower on the
  current focus parent solve, so `auto` intentionally remains legacy.

| Workload | TITAN seconds | pyDelPhi seconds | Speed ratio | Bottleneck bucket | Current read |
| --- | ---: | ---: | ---: | --- | --- |
| Barnase linear | 18.84 | 32.05 | 1.702 | epsilon | `pass`; both rows report 2 threads, and TITAN remains faster in the median profile. |
| Barnase NWT | 106.64 | 124.09 | 1.164 | solver | Strict parity pass. Auto half-grid fair timing now uses the `salt-tiled` homogeneous update: `solver=89.09s`, `nwt_half_update=57.45s`, `nwt_coeff_pack=0.16s`, `nwt_homo_update=41.93s`, `nwt_hetero_update=2.27s`, `nwt_homo_ns_per_update=4.24`. The two raw repeats showed scheduler noise (`1.042x` and `1.251x`), so the locked fair profile is median `RepeatCount=2`. |
| 5tif gaussian-cutoff | 23.62 | 57.50 | 2.434 | epsilon | Strict parity pass. Gaussian remains faster overall; epsilon is `8.10s` and solver is `7.03s` in the median profile. |
| 5tif RPB | 63.13 | 136.86 | 2.168 | rpb_source | Strict parity pass; `rpb_source=35.39s`, still faster overall than pyDelPhi. |
| Barnase focus/probe | 32.85 | 8.66 | 0.264 | solver | Strict parity pass, but still not apples-to-apples: TITAN times the full 1697-atom parent/child/probe workflow, while pyDelPhi times a 599-atom child oracle probe. A synthetic TITAN child-scope report-only row is `2.31s` versus pyDelPhi `8.66s`, with true child atom count `599`, child solver `0.84s`, and parent solve `16.31s`. |

Short version: thread/scope accounting is explicit, the locked RPB source row is
no longer the slowest hotspot, NWT fair timing no longer pays strict diagnostic
passes, and Gaussian epsilon now reports the density/ion/finalize split. The
remaining performance work should follow the measured buckets: on comparable
rows TITAN is now faster than pyDelPhi at the same recorded thread count, while
full focus/probe remains a workflow-scope comparison and the report-only child
scope is faster than the pyDelPhi child oracle. Gaussian midpoint-density
construction and the locked RPB source row are no longer the rows to chase
first.

## Difference Matrix Expansion

`test/pbe_inputs_local/pbe_delphi_golden.csv` now includes hidden harness-only
columns `surface_cutoff`, `pbe_timing`, and `parity_gate`. Existing locked rows
remain `strict`; the expanded Gaussian/RPB/focus-probe rows start as
`report-only` until their source/operator/final-energy oracle data is finite and
stable enough for strict gates.

New report-only coverage includes Gaussian salt `0` and `2.0`, Gaussian-charge
salt `0` and `2.0`, explicit gaussian-cutoff surface cutoff, density cutoff
`0.50`, RPB salt `0` and `2.0`, RPB Gaussian-charge/surface-cutoff variants,
and focus/probe residue/site/no-map variants. Use
`test/pbe_parity_difference_matrix.ps1` with `-Fixtures` or `-MaxRows` for
targeted runs before promoting any row from report-only to strict.

`test/pbe_parity_difference_matrix.ps1 -Fixtures barnase_linear_salt015` now
includes explicitly requested strict rows instead of only report-only rows. The
current row reports `strict_pass` with `delta=7.04121652685171 kJ/mol` and
`tolerance=152.839878735806 kJ/mol`.

## Verification Evidence To Keep Handy

Commands that have recently passed in native mode:

```powershell
mingw32-make PROFILE=native
powershell -ExecutionPolicy Bypass -File test\pbe_delphi_export_maps.ps1 -PythonExe C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe -FullOperatorOracle -FullFieldOracle -NwtOnly -FullMidpointSalt 0.15
powershell -ExecutionPolicy Bypass -File test\pbe_delphi_export_maps.ps1 -PythonExe C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe -FullGaussianOracle -FullRpbOracle
powershell -ExecutionPolicy Bypass -File test\pbe_delphi_export_maps.ps1 -PythonExe C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe -FullFocusProbeOracle
powershell -ExecutionPolicy Bypass -File test\pbe_regression.ps1
powershell -ExecutionPolicy Bypass -File test\pbe_iscm_energy_smoke.ps1
powershell -ExecutionPolicy Bypass -File test\pbsa_smoke.ps1
powershell -ExecutionPolicy Bypass -File test\pbe_snapshot_stability_scan.ps1
powershell -ExecutionPolicy Bypass -File test\pbe_nwt_solver_progress_smoke.ps1 -Strict -StageTsv test\tmp\pbe_delphi_export_maps\barnase_stage_oracle.tsv
powershell -ExecutionPolicy Bypass -File test\pbe_nwt_kernel_modes_smoke.ps1
powershell -ExecutionPolicy Bypass -File test\pbe_gaussian_density_oracle_smoke.ps1 -OracleRoot test\tmp\pbe_delphi_export_maps -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_gaussian_parity_smoke.ps1 -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_rpbe_source_smoke.ps1 -OracleRoot test\tmp\pbe_delphi_export_maps -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_rpbe_smoke.ps1 -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_focus_boundary_oracle_smoke.ps1 -OracleRoot test\tmp\pbe_delphi_export_maps -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_focus_probe_smoke.ps1 -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_timing_compare.ps1 -PythonExe C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe
powershell -ExecutionPolicy Bypass -File test\pbe_timing_compare.ps1 -PythonExe C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe -Fixtures barnase_nwt
powershell -ExecutionPolicy Bypass -File test\pbe_timing_compare.ps1 -PythonExe C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe -Fixtures 5tif_gaussian_cutoff_salt015
powershell -ExecutionPolicy Bypass -File test\pbe_timing_compare.ps1 -PythonExe C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe -Fixtures barnase_focus_probe
powershell -ExecutionPolicy Bypass -File test\pbe_parity_difference_matrix.ps1 -Fixtures barnase_linear_salt015
$env:TITAN_NO_PAUSE='1'; .\compile.bat portable all
```

Phase 2 RPB source-audit lock:

```powershell
powershell -ExecutionPolicy Bypass -File test\pbe_rpbe_source_smoke.ps1 -OracleRoot test\tmp\pbe_delphi_export_maps -Strict
```

This now passes. The key fix was to make hidden PQR `test_center=auto` reproduce
pyDelPhi's `COORDLIMITS` grid center from the radius-expanded boundary box. The
current 5tif RPB source row reports TITAN `source_rms=1341.232978` and
`source_max=39239.94855`, inside the strict oracle gate. With parity dump
enabled, TITAN also writes `pbe_rpbe_source_top.csv` to compare the strongest
native source cells against oracle top-source diagnostics.

The 2026-05-06 timing pass also changed the production-side source audit to use
the caller executor/thread count instead of a private worker cap, pre-skip
Gaussian-gradient contributions outside the locked density cutoff, and reuse the
solvent RPB audit/source/Coulomb maps for energy/top diagnostics when available.
The current locked 5tif full timing row has `rpb_source=49.87s`, still down
substantially from the older 164s snapshot and still faster overall than
pyDelPhi, while strict source/final-energy gates still pass.

The NWT solver pass keeps strict block-progress semantics for oracle gates, but
allows fair timing rows to set `pbe_nwt_progress=off`. The half-grid workspace
now persists across blocks, caches split RHS/epsilon/screen coefficients, packs
interior cells into half-neighbor indices, separates homogeneous from
heterogeneous coefficient updates, and routes homogeneous cells through a
salt-tiled fixed-offset kernel. The `auto` homogeneous kernel keeps
`packed-list`, `row-span`, and `avx2` available for bisection, but now selects
`salt-tiled` for fair timing after strict gates pass. That path tiles
homogeneous salt spans, hoists color-specific pointers/offsets, and uses a
double-precision polynomial `sinh`/`cosh` fast branch for small clipped
potentials without enabling fast-math. On the current median full snapshot,
TITAN is `106.64s` total with `solver=89.09s`, down from the row-span/AVX2
`149.25s` / `134.82s` row, the older `179.02s` / `168.18s`, the older
`232.62s` / `221.11s` persistent-half row, and the older `273.50s` / `253.35s`
fair-timing row. pyDelPhi reports `124.09s` total with `solver=94.54s` at
`pydelphi_threads=2`, so the comparable NWT row is now faster than pyDelPhi in
the locked `RepeatCount=2` profile.

The linear half-packed SOR path was also implemented behind
`-pbe-linear-kernel=half-packed`, but the focus parent solve was slower
(`17.81s` median parent solve with half-packed versus `13.95s` in the final
legacy-auto full snapshot). It remains a forced bisect/perf lab path; default
`auto` stays legacy to avoid focus/probe regressions.

The Gaussian epsilon pass hoists invariant Gaussian density parameters, uses a
small atom SoA, parallelizes center/midpoint density, ion exclusion, and the
initial epsilon finalize pass through `TitanExecutor`, and keeps strict
Gaussian/RPB oracle gates green. On the current 5tif gaussian-cutoff full row,
TITAN is `23.62s` total with `epsilon=8.10s`, still faster than pyDelPhi's
`57.50s` full-row timing.

The focus/probe timing compare now emits a synthetic
`titan_child_scope` report-only row with the true child atom count. The current
full workflow is `32.85s` because parent solve is `16.31s` and child epsilon is
`1.44s`, but the child-scope estimate is `2.31s` versus pyDelPhi child oracle
`8.66s`. Keep using
the full row for production cost and the child row only for apples-to-apples
diagnostics.

Exporter note: `test\pbe_delphi_export_maps.ps1` now cleans only the oracle
files for the requested export group. A NWT-only export no longer deletes
Gaussian/RPB/focus oracle files in the same `OracleRoot`.

Command to refresh the timing snapshot:

```powershell
powershell -ExecutionPolicy Bypass -File test\pbe_timing_compare.ps1 -PythonExe C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe
```

Exporter caveat from the 2026-05-06 timing/difference pass: a full all-oracle
refresh with `-FullOperatorOracle -FullFieldOracle -FullGaussianOracle
-FullRpbOracle -FullFocusProbeOracle` exceeded a 900 second local timeout. The
strict files were restored by targeted Gaussian, RPB, focus/probe, and NWT
exports. Prefer targeted exporter groups when only one strict gate is missing.

## Recommended Next Work Order

1. Use `test/pbe_snapshot_stability_scan.ps1 -PqrPaths <snapshots>` on any
   real trajectory subset that shows sudden PBE energy jumps. Treat
   `grid_boundary_candidate` as evidence for atom/grid-plane sensitivity, but
   check `grid_center`, `surface_midpoint`, `ion_access`, and `iscm_boundary`
   labels before changing solver or energy conventions.
2. Keep the strict gates green while promoting the new report-only fixtures one
   stage at a time: grid/radius, surface/midpoint, charge, boundary, operator,
   field, ISCM/energy, then nonlinear components.
3. Keep NWT `salt-tiled` guarded by strict/kernel-mode smokes; the next
   comparable-row performance work should target focus/probe same-scope
   harnessing and parent-solve/setup reuse, not another NWT rewrite unless a
   fresh median profile shows regression.
4. Add a second RPB optimization pass only after newly promoted RPB matrix rows
   show source construction as a remaining bottleneck; the current locked 5tif
   RPB row is already faster than pyDelPhi at the same recorded thread count.
5. Continue the real module split only after the current timing/parity snapshot
   is stable. Next good slices are nonlinear/NWT helpers or ISCM/energy audit
   helpers, one or two subsystems at a time.
6. Refresh `test/tmp/pbe_timing_compare/pbe_timing_compare.csv` after each
   performance change with median `RepeatCount=2`; single-run ratios can swing
   enough to mislead the gate.
7. Run portable build only after native parity/timing stabilizes.
