# Poisson-Boltzmann examples

These examples use public structures and the manuscript-facing TITAN protocol.
All paths are local to the cloned repository; no download occurs at run time.

## Protocol shared by the protein examples

- linear Poisson-Boltzmann equation at 298.15 K and 0.15 M monovalent salt;
- solute/solvent dielectric constants 1/80;
- probe radius 1.4 Angstrom;
- watertight analytical solvent-excluded surface, spacing 1.25 Angstrom;
- node-patch boundary-integral discretization;
- degree-one analytical FMM/M2L, opening parameter 0.45;
- analytical singular near integration;
- sparse FMM-near-graph block ILU(0), neighbor cap 64;
- GMRES true-residual tolerance 1e-4; and
- four CPU threads.

The option `-pbe-tabi-taylor-order=4` remains in the argument files because it
is part of the current parser contract. With `-pbe-tabi-fmm-local=off`, the
production analytical far operator is the degree-one path described above.

## Included runs

1. `01_1vii_solvation`: polar solvation energy for the villin headpiece
   subdomain (PDB 1VII).
2. `02_1vii_surface_potential`: the same solve with CSV and VTK surface fields
   for visualization.
3. `03_6lzg_binding_cycle`: complex/receptor/ligand polar electrostatic cycle
   for PDB 6LZG chains A and B. This is a software-validation example on a
   fixed experimental structure, not a biological binding claim.

Each `run.ps1` accepts `-TitanExe` and recreates a local `run/` directory. The
curated `reference_*.csv` files record expected outputs from the tagged source
revision; small timing differences across processors are normal.

