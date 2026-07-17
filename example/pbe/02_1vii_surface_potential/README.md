# 1VII surface electrostatic-potential export

This example uses the same physical and numerical protocol as the solvation
example, with `-pbe-bie-surface-output=all`. TITAN writes:

- `pbe_tabi_surface_solution.csv` for analysis and plotting;
- `pbe_tabi_surface_solution.vtk` for PyMOL/ParaView-compatible workflows;
- `pbe_atom_reaction_potential.csv` for per-atom reaction potentials; and
- the standard energy, timing, memory, surface-quality, and solver profiles.

The exported potential is a TITAN boundary-solution field. Cross-software
surface comparisons should use a common geometry and a common color range.

