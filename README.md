# TITAN biomolecular Poisson-Boltzmann solver

This repository is the compact, manuscript-facing source release for TITAN's
boundary-integral Poisson-Boltzmann workflow. It contains the C++23 source,
the command reference, and reproducible example runs used to demonstrate:

- watertight analytical solvent-excluded surfaces;
- node-patch boundary-integral discretization;
- degree-one analytical FMM/M2L matrix-vector products;
- sparse near-field block ILU(0) preconditioning;
- polar solvation and binding-energy output; and
- export of surface electrostatic potentials for visualization.

The broader private TITAN development repository contains additional molecular
simulation and analysis modules. This public tree is intentionally limited to
the source snapshot and examples needed to build and audit the manuscript code.

## Build

TITAN requires a C++23 compiler and OpenMP. Eigen is vendored under
`src/eigen`.

### Windows (MSYS2 UCRT64)

```powershell
mingw32-make PROFILE=native TARGET=titan.exe all -j4
```

The Makefile defaults to `C:/msys64/ucrt64/bin/clang++.exe`. Override
`MSYS2_UCRT64_BIN` if MSYS2 is installed elsewhere.

### Linux

```bash
make PROFILE=native TOOLCHAIN=gcc all -j4
```

GCC 13 or newer is recommended. Use `PROFILE=portable` for an AVX2/FMA build
that does not use `-march=native`.

## Examples

Start with [examples/pbe/README.md](examples/pbe/README.md). On Windows, after
building `titan.exe`, run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File examples\pbe\01_1vii_solvation\run.ps1
```

Each example executes in its own ignored `run/` directory and leaves the
curated inputs unchanged.

## Command reference

See [manual.txt](manual.txt) or run `titan.exe -h`.

## License

TITAN is released under the [MIT License](LICENSE).
