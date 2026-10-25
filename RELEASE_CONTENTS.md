# Contents of titan-bie-pb-v0.6

- `src/`, `Makefile`, `VERSION`, `tools/`: the TITAN C++23 source and its build (Eigen included).
- `manual.txt`: the command-line reference of TITAN.
- `MANUAL_BIE_PB.md`: the user manual of the boundary-integral Poisson-Boltzmann solver.
- `example/pbe/`: solvation energy (1VII), surface potential map (1VII), polar binding cycle (6LZG), and the 27
  proteins of the paper from the deposited PDB files, each with its inputs and reference results.
- `validation/rsos_v06/`: the result tables behind the paper's tables and figures, the campaign and comparison
  scripts, the APBS and PAFMPB comparison inputs and notes, and `SHA256SUMS`.

No executable, build object or temporary run directory is included.
