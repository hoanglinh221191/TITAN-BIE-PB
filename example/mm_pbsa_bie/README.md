# TITAN BIE MM/PBSA Example

This self-contained example applies TITAN's boundary-integral Poisson-Boltzmann
backend to the complete protein-ligand trajectory from `example/1`. The input
contains 52 `MODEL` records from 89.8 to 100.0 ns. TITAN uses the first model
as the coordinate/topology seed and analyzes the remaining 51 frames from
90.0 to 100.0 ns. Before any energy evaluation, every analyzed frame is made
whole with TITAN's
FixBox-style periodic-boundary module (`-pbcn`). TITAN then evaluates the
complex, receptor, and ligand with identical electrostatic settings and writes
the conventional MM/PBSA estimate

`Delta G = E_MM + DeltaDelta G_solv,PB + Delta G_nonpolar,SASA`.

Interaction entropy is evaluated from the full trajectory and is reported both
separately and in the final binding-energy estimate. The example is intended to
accompany the manuscript. It uses the production
BIE protocol: a watertight analytical solvent-excluded surface, node-patch
discretization, hierarchical Cartesian multipole acceleration, sparse near-field
ILU(0), a 0.45 opening parameter, and a `1e-4` true-residual target.

The native MM/PBSA path does not require separately prepared complex,
receptor, ligand, or frame-by-frame PQR files. The two molecular components
are selected from `index.ndx` with `-g=<group1>-<group2>`. For every frame,
TITAN extracts those two groups as molecule 1 and molecule 2 and forms their
union as the complex endpoint. The bundled example uses a multi-model PDB,
but the same workflow accepts a reference structure and an XTC/TRR trajectory:

```text
-s=complex.pdb
-traj=complex.xtc
-n=index.ndx
-g=<group1>-<group2>
```

Charges, radii, and Lennard-Jones parameters can be read from a user-supplied
GROMACS/AMBER topology (`-common-topo`, or `-topo` and `-topo2` for separate
components). When no topology is supplied, `-fflib-auto=on` assigns the
parameters from the force-field data in `lib_ff`. The same assigned parameters
are retained for the complex and both isolated endpoint calculations.

## Run

From the repository root:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File example\mm_pbsa_bie\run_mm_pbsa_bie.ps1 `
  -TitanExe .\titan.exe `
  -Threads 4
```

Use `-Force` to replace an existing `results` directory. A complete run checks
that all 51 sampled frames associated with the 52-model input were evaluated,
the PBC stage appeared in the run log, every endpoint cycle succeeded, and the
decomposition closes to the frame-level MM/PBSA total.

## Decomposition semantics

`-pbsa-decomp=hybrid` reports protein contributions by residue and ligand
contributions by atom. `-edecomp-source=pbsa -edecomp-level=both` additionally
writes residue and atom maps/rankings for both endpoints. Polar contributions
are differences of per-atom reaction-field energies between the complex and
isolated endpoint. Molecular-mechanics pair energies are divided equally
between their two partners, and the nonpolar term is distributed from the
per-atom change in solvent-accessible surface area. The output includes closure
checks against the total no-entropy MM/PBSA energy.

These rows are an additive endpoint decomposition, not a unique pairwise
decomposition of the nonlinear binding free energy. Interaction entropy is a
trajectory-level correction and is therefore not assigned uniquely to residues
or atoms.

## Main outputs

- `mm_pbsa_energy.csv`: frame-level complex/receptor/ligand energy cycle.
- `mm_pbsa_summary.dat`: mean energies and run protocol.
- `mm_pbsa_unit_summary.csv`: hybrid residue/ligand-atom contributions.
- `edecomp_mmpbsa_residue_{map,ranking}.csv`: residue decomposition.
- `edecomp_mmpbsa_atom_{map,ranking}.csv`: atom decomposition.
- `mm_pbsa_performance.dat`: endpoint and decomposition timing.
- `pbe_tabi_profile.csv`: BIE solver, multipole, preconditioner, and audit
  telemetry. The `pbe_tabi_*` stem is retained for output-schema compatibility.

All paths in `arg_bie_mm_pbsa.in` are relative to the generated `results`
directory, so the example remains portable within the repository layout.
