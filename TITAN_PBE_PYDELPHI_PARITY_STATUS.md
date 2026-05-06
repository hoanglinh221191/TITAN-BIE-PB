# TITAN PBE / pyDelPhi Parity Status

Date: 2026-05-05

This note records the current PBE parity state so later phases can continue from
evidence instead of rediscovering the same mismatch chain.

## Runtime Independence

TITAN PBE is intended to run independently. Production TITAN does not require
DelPhi or pyDelPhi at runtime, and should not read oracle CSV files for normal
PBE, PBSA, Gaussian, RPB, focus, or probe calculations.

pyDelPhi/DelPhi is used only as a local test oracle under `.codex_tmp` and
`test/tmp` to export traces, maps, component energies, and timing references.
Those oracle files are for regression/debug comparison only.

## Source Layout

- The public API is `titan_pbe.h`; it now exposes `FrameResult` and a PIMPL
  `TitanPBE::Analyzer` wrapper so engine call sites keep including the same
  header.
- The implementation lives under `titan_pbe/`. Phase 3 module splitting has
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
- `Makefile` builds `titan_pbe/*.cpp`; `compile.bat` remains unchanged because
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

## Current Performance Snapshot

Timing data below comes from the refreshed
`test/tmp/pbe_timing_compare/pbe_timing_compare.csv`. It is a wall-clock
snapshot on the current machine, not a portable benchmark.

| Workload | TITAN seconds | pyDelPhi seconds | Current read |
| --- | ---: | ---: | --- |
| Barnase linear | 23.15 | 24.24 | Similar speed, TITAN about 1.05x faster; timing harness still marks this row `mismatch_report_only`. |
| Barnase NWT | 215.68 | 135.18 | TITAN about 1.59x slower; parity status pass with NWT parity settings. |
| 5tif gaussian-cutoff | 30.00 | 76.64 | TITAN about 2.55x faster; parity status pass. |
| 5tif RPB | 178.62 | 136.59 | TITAN about 1.31x slower; parity status pass. |
| Barnase focus/probe | 31.41 | 6.90 | TITAN about 4.55x slower; parity status pass. |

Short version: TITAN is already faster than pyDelPhi on gaussian-cutoff, roughly
similar on this linear barnase timing row, and slower on NWT, RPB, and
focus/probe. The next performance work should target NWT solver cost, RPB
source/solver construction, and focus/probe map I/O/setup after preserving the
current parity gates.

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
powershell -ExecutionPolicy Bypass -File test\pbe_nwt_solver_progress_smoke.ps1 -Strict -StageTsv test\tmp\pbe_delphi_export_maps\barnase_stage_oracle.tsv
powershell -ExecutionPolicy Bypass -File test\pbe_gaussian_density_oracle_smoke.ps1 -OracleRoot test\tmp\pbe_delphi_export_maps -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_gaussian_parity_smoke.ps1 -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_rpbe_source_smoke.ps1 -OracleRoot test\tmp\pbe_delphi_export_maps -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_rpbe_smoke.ps1 -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_focus_boundary_oracle_smoke.ps1 -OracleRoot test\tmp\pbe_delphi_export_maps -Strict
powershell -ExecutionPolicy Bypass -File test\pbe_focus_probe_smoke.ps1 -Strict
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

Exporter note: `test\pbe_delphi_export_maps.ps1` now cleans only the oracle
files for the requested export group. A NWT-only export no longer deletes
Gaussian/RPB/focus oracle files in the same `OracleRoot`.

Command to refresh the timing snapshot:

```powershell
powershell -ExecutionPolicy Bypass -File test\pbe_timing_compare.ps1 -PythonExe C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe
```

## Recommended Next Work Order

1. Continue the real module split: PBSA, focus/probe, Gaussian/RPB,
   grid/surface, charge/operator, linear solver, and linear multigrid slices are
   now separated. Next good slices are nonlinear/NWT helpers or ISCM/energy
   audit helpers. Move only one or two subsystems at a time, with a native
   build/smoke after each slice.
2. Expand Gaussian/RPB fixtures beyond the current oracle set.
3. Add more focus/probe fixtures covering residue/site/group corner cases.
4. Optimize NWT, RPB, and focus/probe after keeping parity strict; those are
   currently slower than pyDelPhi in the refreshed timing snapshot.
5. Run portable build only after native parity/timing stabilizes.
