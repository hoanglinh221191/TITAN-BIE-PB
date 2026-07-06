# TITAN manuscript-release build instructions

These commands document the build and smoke-test environment used by the manuscript release.

## Windows native build

Command:
mingw32-make PROFILE=native TARGET=titan.exe all -j2

The manuscript benchmarks were executed with native double-precision CPU code and AVX2-capable kernels when available.

## Minimal guarded p4 smoke

Commands:
powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_y_treecode_default_probe.ps1 -Quick -TitanExe .\titan.exe
powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fagbi_scope_z4xd_p_order_tiled_cache_gate.ps1 -TitanExe .\titan.exe -Quick -AllowFailedDecision

## Release manifest regeneration

Command:
powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fagbi_release_readiness_manifest.ps1

## APBS context

Supplementary Figure S1 uses APBS input and log files archived under docs/fagbi_guarded_p4_package/release_payload/apbs_context/. APBS itself is external to TITAN; the archived .in files, stdout/stderr logs, and grid-size settings CSV document the exact context rows used in the manuscript supplement.

## Notes

- Production PBE default remains -pbe=grid.
- Guarded local p4 TABI/FMM is an explicit advanced path and is not enabled by default.
- Full-protein dense direct TABI is not required for the manuscript tables; sampled dense direct-operator audits and small dense-oracle rows are archived instead.
