# TITAN PBE Clean-Room Policy

Date: 2026-05-07

This policy defines the forward-looking clean-room boundary for TITAN PBE work.
It is intentionally precise about what is guaranteed and what is not.

## Claim

TITAN PBE is a self-contained C++ implementation validated against black-box
DelPhi/pyDelPhi oracle outputs and public TABI papers. Production TITAN must
not import, link, call, vendor, or read DelPhi/pyDelPhi/TABI-PB/NanoShaper
runtime artifacts.

This is a forward-looking engineering policy. Earlier parity work inspected a
local pyDelPhi checkout while diagnosing behavior, so do not claim strict
historical legal clean-room provenance for all existing PBE code. From this
phase onward, production PBE changes must be based on TITAN code, internal
specs, public documentation, and black-box output comparisons only.

## Allowed

- Production C++ may use compatibility labels such as `delphi-vdw`,
  `delphi-grid`, or `delphi-half` to describe behavior modes.
- Test harnesses under `test/` may run DelPhi/pyDelPhi as external black-box
  oracles and may write CSV/TSV diagnostics under `test/tmp`.
- Benchmark scripts may compare TITAN and pyDelPhi timing, energy, stage, and
  hash outputs.
- TABI-PB may be used only as an optional external black-box oracle or
  benchmark; production `-pbe=tabi` must remain self-contained C++.
- Oracle maintainers may refresh fixture outputs and provenance manifests.

## Forbidden In Production

- No source includes, imports, subprocess launches, dynamic loading, or linking
  against DelPhi, pyDelPhi, TABI-PB, or NanoShaper.
- No production reads from `.codex_tmp/pyDelPhi`, `pydelphi_env`, DelPhi
  executables, Python oracle scripts, pyDelPhi source files, TABI-PB source
  trees, or NanoShaper meshes.
- No production fallback that silently replays `pbe_oracle_*.csv` outside
  hidden test/parity diagnostics.
- No copying or translating DelPhi/pyDelPhi source code into TITAN.

## Oracle Provenance

Oracle scripts should append provenance rows to `test/pbe_oracle_manifest.csv`
with fixture, command class, input hash, output hash, pyDelPhi root/commit when
available, Python executable, timestamp, thread count, status, and output path.
The manifest is a test artifact only; production TITAN must never read it.

## Required Gate

Run the clean-room audit before claiming a PBE change preserves this boundary:

```powershell
powershell -ExecutionPolicy Bypass -File test\pbe_clean_room_audit.ps1
```

The audit fails on production runtime/source references to DelPhi/pyDelPhi and
writes inventories for compatibility labels and test-only oracle usage.
