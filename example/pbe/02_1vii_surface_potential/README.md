# 1VII surface electrostatic potential

The solvation calculation of `01_1vii_solvation` with `-pbe-bie-surface-output=all`. Besides the usual energy,
timing and profile tables, TITAN writes

- `pbe_tabi_surface_solution.csv`: for every surface node its position, normal and cell area, the boundary
  unknowns (`phi`, `phin`) and the total, Coulomb and reaction potentials (the `total_potential_kT_e` column is
  the potential in kT/e used for the maps of the paper);
- `pbe_tabi_surface_solution.vtk`: the same field on the surface mesh, for ParaView or PyMOL.

To draw the map in PyMOL on the TITAN surface mesh, coloured red (negative) to blue (positive) on the ±5 kT/e
scale of Figure 6:

```bash
pymol view_surface_pymol.py -- run/pbe_tabi_surface_solution.vtk 1vii_full_charmm_ph7.pqr
pymol -cq view_surface_pymol.py -- run/pbe_tabi_surface_solution.vtk 1vii_full_charmm_ph7.pqr 5 1vii_map.png
```

In the boundary-integral method the surface potential is an unknown of the solve, so no interpolation is
involved. When maps from different programs are compared, use one colour scale for all of them.
