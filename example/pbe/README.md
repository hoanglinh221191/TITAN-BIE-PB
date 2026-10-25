# Poisson-Boltzmann examples

These examples run the boundary-integral Poisson-Boltzmann protocol of the paper. Their inputs are the PQR files
the paper used: deposited PDB structures protonated with PDB2PQR 3.7.1 (CHARMM charges and radii, PROPKA at pH 7),
as described in `MANUAL_BIE_PB.md`, section 3. Each `arg.in` is the argument file of the same calculation in the
paper, with relative paths. Nothing is downloaded at run time.

| Directory | What it shows |
|---|---|
| `01_1vii_solvation` | ΔG_solv of the villin headpiece subdomain (PDB 1VII) |
| `02_1vii_surface_potential` | The same solve, with the surface potential written as CSV and VTK (`-pbe-bie-surface-output=all`) |
| `03_6lzg_binding_cycle` | Complex, receptor and ligand solves and the polar binding term of the ACE2-RBD complex (PDB 6LZG, chains A and B), Table 6 of the paper |
| `04_fagbi27_protocol` | The 27 proteins of Figure 5: deposited PDB files, their preparation, the protocol, and the paper's results |
| `common` | The neutral four-atom driver pair (`driver.pdb`, `driver.top`) and the minimal force-field library every run needs |

Run an example from this directory:

```bash
python run_example.py 01_1vii_solvation --titan ../../titan.exe
```

or, on Windows, `.\01_1vii_solvation\run.ps1 -TitanExe ..\..\titan.exe`. The run writes into `<example>/run/`.
`reference_results.csv` in each example holds the results of this release; built from the same source, the
energies agree to the printed digits, and other compilers or processors may differ in the last digits.

The protocol, in brief: linear Poisson-Boltzmann equation at 298.15 K and 0.15 M monovalent salt; dielectric 1
inside and 80 outside; the analytic solvent-excluded surface of a 1.4 Å probe at spacing h = 1.25 Å, refined
towards charges near the surface; curved node-patch cells with convergence-tested self quadrature; a degree-one
Cartesian treecode with opening parameter 0.45; block ILU(0) on the near-field graph with 64 blocks per row;
GMRES to a true relative residual of 10^-4; four threads. `MANUAL_BIE_PB.md`, section 5, explains every option.
