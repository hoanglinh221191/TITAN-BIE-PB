# TITAN boundary-integral Poisson-Boltzmann solver

This repository is the source release of TITAN (Thermodynamic & Interaction Trajectory ANalyzer) for the paper
"A Stabilized Boundary-Integral Poisson-Boltzmann Solver for Biomolecular Solvation Electrostatics". It contains
the C++23 source of the whole program, the user manual of the boundary-integral solver, runnable examples, and the
data and scripts behind the paper.

The boundary-integral Poisson-Boltzmann (BIE-PB) solver of TITAN

- builds the analytic solvent-excluded surface of the molecule, refined towards charges close to the surface;
- discretizes it with curved node-patch cells built from the analytic surface normals, and integrates the
  singular and near-singular kernels over them (polar Duffy-type quadrature of the self cells, converged to a set
  tolerance; geometrically adaptive subdivision of nearby cells);
- applies the far field with a hierarchical Cartesian multipole expansion;
- solves the coupled system with GMRES and a sparse block ILU(0) preconditioner built on the near-field graph;
- checks the accelerated operator against direct summation on the converged solution;
- reports polar solvation energies, polar binding terms and surface electrostatic potentials.

The BIE solver is selected with `-pbe=bie`. TITAN's default Poisson-Boltzmann solver remains the grid solver.

## Start here

- [MANUAL_BIE_PB.md](MANUAL_BIE_PB.md): building TITAN, preparing a protein from its PDB file, running the
  paper's protocol, reading and checking the results.
- [example/pbe](example/pbe/README.md): solvation energy, surface potential map, polar binding cycle, and the
  27 proteins of the paper from the deposited PDB files.
- [validation/rsos_v06](validation/rsos_v06/README.md): the results behind every table and figure of the
  paper, with the scripts that produced them.
- [manual.txt](manual.txt): the complete command-line reference of TITAN (`titan -h`).

## Build

TITAN requires a C++23 compiler and OpenMP; Eigen is included under `src/eigen`.

Windows (MSYS2 UCRT64):

```bash
export PATH=/c/msys64/ucrt64/bin:$PATH
mingw32-make all -j8
```

Linux (GCC 13 or newer):

```bash
make PROFILE=native TOOLCHAIN=gcc all -j8
```

`PROFILE=portable` builds for AVX2/FMA without `-march=native`.

## Quick start

```bash
cd example/pbe
python run_example.py 01_1vii_solvation --titan ../../titan.exe
```

prints the node patches, GMRES iterations, residuals and ΔG_solv = -891.625 kcal/mol for the villin headpiece
(PDB 1VII), as in the paper.

## Changes from titan-bie-pb-v0.5

- Curved point-normal cells on the analytic SES, with analytic vertex normals, replace the flat node patches.
- The self cells are integrated by polar Duffy-type quadrature with a convergence test and the nearby cells by
  geometrically adaptive subdivision; a cell that cannot be resolved stops the calculation.
- The analytic SES is refined towards charges close to the surface (factor 0.35, minimum edge 0.1 Å).
- The routine far field (p1) is a degree-one particle-cluster treecode; the degree-four expansion (p4) is the
  accuracy control.
- The sampled operator audit runs on the converged solution (512 targets, 1,024 rows).
- New: `MANUAL_BIE_PB.md`, the 27-protein example with the PDB-to-PQR preparation, the surface-map script, and
  the evidence of the revised paper (`validation/rsos_v06`), which replaces the v0.5 matrix archive.
- The BIE MM/PBSA examples of v0.5 are not part of this release; the MM/PBSA options remain documented in
  `manual.txt`.

## License and citation

TITAN is released under the [MIT License](LICENSE). If you use this release, please cite the paper and this
versioned archive ([CITATION.cff](CITATION.cff)).
