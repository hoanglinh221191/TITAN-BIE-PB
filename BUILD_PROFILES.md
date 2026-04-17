# TITAN Build Profiles and Benchmark Suites

This document defines the Phase 0 build/baseline workflow.

## Build profiles

TITAN supports two build profiles through `Makefile` and `compile.bat`.

- `native`: enables `-march=native` for local benchmarking on the current host.
- `portable`: avoids `-march=native`, uses
  `-mavx2 -mfma -mtune=generic`, and stays within the C++23 subset/codegen
  assumptions we want to keep compatible with `g++13/libstdc++13`.
  This is not scalar-only portability yet; Phase 4 currently centralizes
  runtime-dispatched AVX2/FMA hotspots while keeping the portable baseline
  itself on AVX2/FMA.

Phase 4 slice 1 also centralizes build/runtime capability reporting:

- build profile is exported to the binary as `native` or `portable`
- object files and dependency files live under profile-specific directories:
  `build/native/` and `build/portable/`
- startup logs now report:
  - build profile
  - `g++13/libstdc++13` compatibility floor
  - compiled SIMD level
  - runtime SIMD level
  - whether the current CPU satisfies the binary's compiled SIMD requirement

This does not change numeric logic; it only makes profile/ISA assumptions
explicit and easier to validate on another machine.

Phase 4 slice 2 extends this with shared runtime-dispatch helpers for selected
hot kernels:

- reductions used by LJ and RMS
- row-major squared-distance/add-row helpers used by MSM and UMAP
- dot-product and 4-distance helpers used by DCCM/LSP-MD

These helpers keep scalar fallbacks in one place and let `native` /
`portable` builds share the same dispatch policy. They do not change the
scientific algorithm or output contract.

Phase 4 slice 3 extends the same dispatch layer into classic contact-screening
paths that previously carried their own AVX blocks:

- `titan_ifp.cpp` now uses the shared point-vs-4 cutoff helper for ligand atom
  screening
- `h_bond.cpp` now uses the shared point-vs-4 gathered cutoff helper for
  hydrogen-bond group prefiltering

Regression for this slice is tracked under:

- `test/phase4_slice3_ifp_baseline_native`
- `test/phase4_slice3_ifp_native`
- `test/phase4_slice3_ifp_portable`

Phase 4 slice 4 extends the same approach into the LLPS autocorrelation
popcount path:

- `titan_llps.cpp` now uses a shared XOR-overlap popcount helper from
  `titan_simd_dispatch.h` instead of carrying a local AVX2-only block
- the helper keeps the exact integer counting logic for overlap and reference
  populations, with scalar and AVX2 implementations selected through the same
  runtime-dispatch layer used by other Phase 4 slices
- this does not change the LLPS model or file contract; it only centralizes the
  implementation of the hot popcount kernel

Regression for this slice is tracked under:

- `test/phase4_slice4_llps_baseline_native`
- `test/phase4_slice4_llps_native`
- `test/phase4_slice4_llps_portable`

The current LLPS fixture does not emit `llps_autocorr.dat` because there are
not enough post-percolation frames for `C(t)`; regression for this slice checks
that the file remains consistently absent and that the remaining LLPS outputs
stay hash-identical.

Phase 4 slice 5 extends the shared dispatch layer into DSSP CA-distance
prefiltering:

- `titan_dssp.cpp` now uses the shared `point_vs_four_cutoff_mask(...)` helper
  instead of carrying its own AVX register math for the 9 A backbone CA
  candidate screen
- the downstream hydrogen-bond energy evaluation and residue-pair exclusion
  logic are unchanged; this slice only centralizes the hot 4-lane cutoff test
- the shared helper still uses AVX2 at runtime when available, so this keeps
  the fast path while avoiding another local ISA-specific block

Regression for this slice is tracked under:

- `test/phase4_slice5_dssp_baseline_native`
- `test/phase4_slice5_dssp_native`
- `test/phase4_slice5_dssp_portable`

The accepted hash gates for this fixture are:

- `dssp_mol1.dat`
- `dssp_mol1_content.dat`

Phase 4 slice 6 extends the same dispatch layer into the GBSA Shrake-Rupley
SASA burial filter:

- `gbsa_intrin.cpp` now routes the 96-point burial mask through the shared
  `point_vs_four_cutoff_mask(...)` helper instead of carrying a local AVX2
  register block in `compute_sasa_shrake_rupley()`
- the scientific logic is unchanged: the burial test still uses strict
  `< Rj^2`, and the downstream GBSA/SASA math remains exactly the same
- this slice intentionally leaves the larger GBSA Born-radii and energy kernels
  alone; it only centralizes the hot 4-lane cutoff path with runtime dispatch

Regression for this slice is tracked under:

- `test/phase4_slice6_gbsa_baseline_native`
- `test/phase4_slice6_gbsa_native`
- `test/phase4_slice6_gbsa_portable`

The accepted hash gates for this fixture are:

- `gbsa_summary.dat`
- `gbsa_timeseries_avg.dat`
- `traj_0_gbsa_snapshots.dat`
- `gbsa_decomposition_ligand.dat`
- `gbsa_decomposition_receptor.dat`

Phase 4 slice 7 extends the same dispatch layer into the GBSA effective
Born-radii helper:

- `gbsa_intrin.cpp` no longer carries a local AVX horizontal-sum + Born-radii
  register loop inside `cal_effective_born_radii_avx2()`
- the hot `rho^3 / 3 / r^4` accumulation now goes through shared
  `gbsa_inverse_r4_volume_sum(...)` helpers in `titan_simd_dispatch.h`, with
  scalar, AVX2, and AVX2+FMA implementations selected through the same runtime
  dispatch layer used by earlier Phase 4 slices
- the scientific logic is unchanged: the self/near-zero exclusion still uses
  `r2 > 1e-6`, the volume term still uses exact `/ 3.0`, and the downstream
  GBSA Born-radius formulas remain untouched

Regression for this slice is tracked under:

- `test/phase4_slice7_gbsa_baseline_native`
- `test/phase4_slice7_gbsa_native`
- `test/phase4_slice7_gbsa_portable`

The accepted hash gates for this fixture are:

- `gbsa_summary.dat`
- `gbsa_timeseries_avg.dat`
- `traj_0_gbsa_snapshots.dat`
- `gbsa_decomposition_ligand.dat`
- `gbsa_decomposition_receptor.dat`

Phase 4 slice 8 extends the same dispatch layer into the GBSA LJ energy helper
path used by the AMBER/OPLS-style GBSA kernels:

- `gbsa_intrin.cpp` no longer carries duplicated local AVX register math for
  the mixed-sigma/mixed-epsilon four-lane LJ path inside `cal_energy_gbsa_avx2()`
- the hot `4 * epsilon * (sr12 - sr6)` evaluation for mixed sigma/epsilon now
  goes through shared
  `gbsa_lj_mixed_sigma_epsilon_r6_fourlanes(...)` helpers in
  `titan_simd_dispatch.h`, with scalar and AVX2 implementations selected
  through the same runtime dispatch layer used by earlier Phase 4 slices
- the scientific logic is unchanged: OPLS still uses geometric sigma mixing,
  the AMBER-style path still uses arithmetic sigma mixing, epsilon still uses
  geometric mixing, the `r^6` floor remains `1e-12`, and the existing
  `r_nm >= 0.05` post-mask is applied exactly as before

Regression for this slice is tracked under:

- `test/phase4_slice8_gbsa_baseline_native`
- `test/phase4_slice8_gbsa_native`
- `test/phase4_slice8_gbsa_portable`

The accepted hash gates for this fixture are:

- `gbsa_summary.dat`
- `gbsa_timeseries_avg.dat`
- `traj_0_gbsa_snapshots.dat`
- `gbsa_decomposition_ligand.dat`
- `gbsa_decomposition_receptor.dat`

Phase 4 slice 9 extends the same dispatch layer into the AMBER energy helper
paths used by the filtered and no-filter intrinsics:

- `energy_intrin.cpp` and `energy_intrin_nofilter.cpp` no longer carry
  duplicated local AVX register math for the four-lane mixed-sigma /
  mixed-epsilon Lennard-Jones energy block
- the hot `epsilon * (c12 - c6)` evaluation now goes through shared
  `lj_mixed_sigma_epsilon_distance_fourlanes(...)` helpers in
  `titan_simd_dispatch.h`, with scalar and AVX2 implementations selected
  through the same runtime dispatch layer used by earlier Phase 4 slices
- the scientific logic is unchanged: AMBER sigma mixing remains arithmetic,
  epsilon mixing remains `sqrt(16 * eps_i * eps_j)`, the screening path still
  uses `safe_d` for the LJ helper while preserving the existing self-filter on
  the original distance, and the CHARMM branches are untouched

Regression for this slice is tracked under:

- `test/phase4_slice9_energy_baseline_native`
- `test/phase4_slice9_energy_native`
- `test/phase4_slice9_energy_portable`
- `test/phase4_slice9_energy_screening_baseline_native`
- `test/phase4_slice9_energy_screening_native`
- `test/phase4_slice9_energy_screening_portable`

The accepted hash gates for this fixture are all generated outputs except
`arg.in` and `log.txt`; the baseline/native/portable directories are
hash-identical for both the default and screened runs.

Phase 4 slice 10 extends the shared dispatch layer into the AMBER/CHARMM
force helper path:

- `force_intrin.cpp` no longer carries duplicated local AVX register math for
  the four-lane mixed-sigma / mixed-epsilon Lennard-Jones block in
  `force_energy_add_intrin(...)` for the AMBER/CHARMM branch
- the hot LJ path now goes through
  `lj_mixed_sigma_epsilon_distance_c6_fourlanes(...)` in
  `titan_simd_dispatch.h`, returning both `c6` and energy with scalar and AVX2
  implementations selected through the shared runtime dispatch layer
- the scientific logic is unchanged: AMBER sigma mixing remains arithmetic,
  epsilon mixing remains `sqrt(eps_i * eps_j)`, the existing `c6`-based
  self-interaction guard is preserved, and the OPLS branch is intentionally
  left untouched in this slice to keep the blast radius small

Regression for this slice is tracked under:

- `test/phase4_slice10_force_baseline_native`
- `test/phase4_slice10_force_native`
- `test/phase4_slice10_force_portable`

The accepted hash gates for this fixture are all generated outputs except
`arg.in`, `log.txt`, and `run_stdout.txt`; the baseline/native/portable
directories are hash-identical.

Independent compile-speed note: `titan_rf_model.cpp` is now forced to compile
with `-O0` by the Makefile in every profile, instead of relying on a source
pragma. This keeps the file on the portable `g++13/libstdc++13` floor and
avoids compiler-specific pragma handling.

Phase 4 slice 11 extends the shared dispatch layer into the legacy OPLS force
helper path:

- `force_intrin.cpp` no longer carries duplicated local AVX register math for
  the four-lane OPLS Lennard-Jones block in `force_energy_add_intrin(...)`
- the hot OPLS path now goes through
  `lj_sigma_product_distance_c6_fourlanes(...)` in `titan_simd_dispatch.h`,
  returning both `c6` and energy with scalar and AVX2 implementations selected
  through the shared runtime dispatch layer
- the scientific logic is unchanged relative to the pre-slice AVX path: the
  helper preserves the legacy OPLS lane math order based on sigma-product
  lanes, keeps the existing `c6`-based self-interaction guard, and leaves the
  AMBER/CHARMM branch from slice 10 untouched

Regression for this slice is tracked under:

- `test/phase4_slice11_force_opls_baseline_native`
- `test/phase4_slice11_force_opls_native`
- `test/phase4_slice11_force_opls_portable`

The accepted hash gates for this fixture are all generated outputs except
`arg.in`, `log.txt`, and `run_stdout.txt`; the baseline/native/portable
directories are hash-identical.

Phase 4 slice 12 closes the remaining low-risk AVX-local cleanup and records
the audit boundary for the rest:

- `dccm_network.cpp` no longer carries duplicated local AVX distance helpers or
  stale dot-product fallback code; LSP-MD now calls
  `titan_simd_dispatch::compute_4dist(...)` and the shared dot-product helper
  directly
- `MOL.h` no longer carries a local `_mm256` reduction loop for thread-local
  force accumulation; it now uses `add_row_inplace(...)` from the shared
  dispatch layer for the hot vectorizable prefix and keeps the scalar tail
- `find_pull_direction.cpp` no longer carries a local AVX translate loop for
  scratch coordinates; it now uses `translate_xyz(...)` from the shared
  dispatch layer
- after this slice, the remaining `_mm256` call sites in the tree are
  intentional:
  `energy_intrin.cpp`, `energy_intrin_nofilter.cpp`, `force_intrin.cpp`,
  `gbsa_intrin.cpp`, `titan_gist.cpp`, and `cal_sub_groups.cpp` are retained as
  full compute kernels or highly specialized microkernels;
  `avx_math.h`, `BoxMatrix.h`, `functions.h`, and `simd_compat.h` remain SIMD
  support headers; `find_pull_direction.cpp` still has `__AVX__` gating only to
  select pre-existing intrinsic backends, not to carry local register math

Regression for this slice is tracked under:

- `test/phase4_slice12_native_lspmd_example2`
- `test/phase4_slice12_portable_lspmd_example2`
- `test/phase4_slice12_force_opls_native`
- `test/phase4_slice12_force_opls_portable`
- `test/phase4_slice12_pull_native`
- `test/phase4_slice12_pull_native_rerun`
- `test/phase4_slice12_pull_portable`

The accepted hash gates for this slice are:

- LSP-MD: `lspmd_edges.csv`, `lspmd_degree_centrality.csv`,
  `lspmd_betweenness_centrality.csv`; native matches the existing phase-4
  baseline and portable matches native
- OPLS force: `cal_force.txt`, `cal_force.xvg`, `total_force.dat`,
  `total_energy.dat`, `electrostatic_force.dat`, `vanderwaals_force.dat`, and
  the main decomposition outputs; portable matches native and the remaining
  text-only diffs versus the older baseline are limited to path formatting in a
  few stats files
- pull direction: `DE_out_vector.dat`, `DE_vector.pdb`, `best_conf.pdb`, and
  the aligned snapshot outputs are deterministic across native reruns and
  portable matches native; `DE_log.dat` / `DE_report.dat` may still differ
  because they embed execution-profile/report text rather than scientific state

## Windows commands

Build the default native profile:

```bat
set TITAN_NO_PAUSE=1
compile.bat native
```

Build the portable profile:

```bat
set TITAN_NO_PAUSE=1
compile.bat portable
```

Each successful `all` build snapshots the executable into:

- `test/bin/native/titan.exe`
- `test/bin/portable/titan.exe`

You can override the output executable name for local testing without touching
the default artifact path by setting `TITAN_TARGET`, for example:

```bat
set TITAN_NO_PAUSE=1
set TITAN_TARGET=titan_stage.exe
compile.bat portable
```

## Benchmark suites

`test/benchmark_harness.ps1` supports:

- `-Suite small`
- `-Suite medium`
- `-Suite large`
- `-Suite all`

It also records `-BuildProfile native|portable|unknown` in the markdown summary.

### Small

Use for quick smoke/regression after hot-path edits.

- screened electrostatics valid/invalid guards
- mini XTC GBSA path
- single-thread GIST sanity

### Medium

Use for data-layout refactors and normal post-processing regression.

- chunk streamer full example
- PCA example
- LLPS example
- multi-thread GIST sanity

### Large

Use for memory-sensitive work and native/portable comparisons.

- XTC TICA + UMAP local case
- optional heavier all-atom case via explicit `-Cases allatom_8cbf_ba1`

## Example benchmark commands

Run the small native suite:

```powershell
powershell -ExecutionPolicy Bypass -File .\test\benchmark_harness.ps1 `
  -Suite small `
  -BuildProfile native `
  -TitanExePath .\test\bin\native\titan.exe
```

Run the small portable suite:

```powershell
powershell -ExecutionPolicy Bypass -File .\test\benchmark_harness.ps1 `
  -Suite small `
  -BuildProfile portable `
  -TitanExePath .\test\bin\portable\titan.exe
```

Run an opt-in heavy case:

```powershell
powershell -ExecutionPolicy Bypass -File .\test\benchmark_harness.ps1 `
  -Cases allatom_8cbf_ba1 `
  -BuildProfile native `
  -TitanExePath .\test\bin\native\titan.exe
```
