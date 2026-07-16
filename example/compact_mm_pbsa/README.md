# Compact TITAN MM/PBSA example

This two-model trajectory is a compact submission companion derived from the
full example distributed in the TITAN repository. It demonstrates the native
workflow without prebuilt PQR endpoint files.

`index.ndx` and `-g=1-13` select the protein and ligand. TITAN extracts the two
molecules for each analyzed frame and forms their union as the complex. Atomic
parameters are read from the supplied GROMACS topology and force-field data.
The first PDB model seeds the system and the second model is analyzed.

From this directory:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File run.ps1 -TitanExe C:\path\to\titan.exe
```

For XTC/TRR data, replace the multi-model input with `-s=complex.pdb` and
`-traj=trajectory.xtc` while retaining the same index groups and topology.
