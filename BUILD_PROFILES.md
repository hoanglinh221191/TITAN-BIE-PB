# TITAN Build Profiles and Benchmark Suites

This document defines the Phase 0 build/baseline workflow.

## Build profiles

TITAN supports two build profiles through `Makefile` and `compile.bat`.

- `native`: enables `-march=native` for local benchmarking on the current host.
- `portable`: avoids `-march=native`, uses
  `-mavx2 -mfma -mtune=generic`, and stays within the C++23 subset/codegen
  assumptions we want to keep compatible with `g++13/libstdc++13`.
  This is not scalar-only portability yet; full runtime dispatch remains a
  later roadmap item.

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
