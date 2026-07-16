# TITAN biomolecular Poisson-Boltzmann solver

This repository is the manuscript-facing source release for TITAN
(Thermodynamic & Interaction Trajectory ANalyzer). It contains the C++23
source, command reference, and reproducible examples supporting the paper
"Stabilized Boundary-Integral Poisson-Boltzmann Solvers for Biomolecular
Solvation Electrostatics."

The release centers on two contributions:

1. a sparse near-field block ILU(0) preconditioner that stabilizes GMRES for
   the coupled boundary-integral Poisson-Boltzmann operator; and
2. an end-to-end TITAN workflow that constructs and audits a watertight
   solvent-excluded surface, solves the boundary-integral system, and exports
   polar solvation energies, PB binding-cycle terms, MM/PBSA decomposition,
   and surface electrostatic potentials.

The implemented workflow includes:

- watertight analytical solvent-excluded surfaces;
- node-patch boundary-integral discretization;
- degree-one analytical FMM/M2L matrix-vector products;
- sparse near-field block ILU(0) preconditioning;
- polar solvation and binding-energy output; and
- export of surface electrostatic potentials for visualization.

The default production PB solver remains the grid backend. The analytical
TABI/FMM path documented in the manuscript is selected explicitly. This public
tree contains the source snapshot and examples needed to build, run, and audit
the reported workflow.

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

Start with [example/pbe/README.md](example/pbe/README.md). On Windows, after
building `titan.exe`, run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File example\pbe\01_1vii_solvation\run.ps1
```

Additional examples are provided in:

- `example/mm_pbsa_tabi/`: trajectory-native complex/receptor/ligand TABI
  MM/PBSA with residue/atom decomposition and reference outputs.
- `example/compact_mm_pbsa/`: a compact two-model submission companion.

The MM/PBSA path reads PDB/XTC/TRR trajectories plus index groups, reconstructs
periodic molecules, forms the complex from the two selected molecular groups,
and assigns charges/radii from a supplied topology or the bundled force-field
library. Separately prepared endpoint PQR files are not required.

## Command reference

See [manual.txt](manual.txt) or run `titan.exe -h`.

## License

TITAN is released under the [MIT License](LICENSE).
