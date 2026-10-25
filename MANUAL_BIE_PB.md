# TITAN BIE-PB: user manual

This manual covers the boundary-integral Poisson-Boltzmann (BIE-PB) solver of TITAN as used in the paper
"A Stabilized Boundary-Integral Poisson-Boltzmann Solver for Biomolecular Solvation Electrostatics". It explains
how to prepare a protein, run the paper's protocol, read the results, and check them. The complete command-line
reference of TITAN, including its molecular-dynamics analyses and the grid PB solver, is `manual.txt`.

Contents

1. What the solver computes
2. Building TITAN
3. From a PDB file to the PB input
4. Running a calculation
5. The protocol of the paper, option by option
6. Results
7. Accuracy checks and failures
8. Variations of the protocol
9. Examples
10. Reproducing the paper

## 1. What the solver computes

TITAN solves the linearized Poisson-Boltzmann equation for a solute of dielectric 1 in an ionic solvent of
dielectric 80. The unknowns are the electrostatic potential and its normal derivative on the solvent-excluded
surface (SES) of a 1.4 Å probe; they satisfy the derivative boundary-integral formulation of Juffer et al., as
in TABI. From the solution TITAN reports

- the polar solvation energy, ΔG_solv = ½ Σ_i q_i φ_reac,i, summed over the atomic charges;
- for a complex and its two components, the polar binding term
  ΔG_pol,bind = G_pol(complex) - G_pol(receptor) - G_pol(ligand);
- optionally, the potential on the surface nodes for visualization.

The surface is the analytic SES of the atoms, traced on a grid of spacing h and refined towards charges close to
the surface. Each mesh vertex carries one node patch, lifted to a curved cell built from the analytic surface
normals. Self and near interactions are integrated over the curved cells by quadrature; far interactions use a
hierarchical Cartesian multipole expansion of degree one (p1). GMRES solves the coupled system with a sparse
block ILU(0) preconditioner built on the near-field graph.

## 2. Building TITAN

TITAN needs a C++23 compiler with OpenMP. Eigen is included under `src/eigen`.

Windows (MSYS2 UCRT64, the build used for the paper):

```bash
export PATH=/c/msys64/ucrt64/bin:$PATH
mingw32-make all -j8
```

Linux:

```bash
make PROFILE=native TOOLCHAIN=gcc all -j8
```

GCC 13 or newer is recommended. `PROFILE=portable` builds for AVX2/FMA without `-march=native`. The build writes
`titan.exe` (Windows) or `titan` (Linux) in the repository root. `titan -v` prints the version.

## 3. From a PDB file to the PB input

TITAN takes the charges and radii of the solute from a PQR file. The paper prepared every protein in the same two
steps, and `example/pbe/04_fagbi27_protocol/prepare_pqr.py` repeats them:

1. Keep the first model and its `ATOM` records (protein chains; `HETATM` records such as waters and ligands are
   dropped).
2. Run PDB2PQR 3.7.1, which adds hydrogens, assigns protonation states with PROPKA at pH 7.0 and writes CHARMM
   charges and radii:

   ```bash
   pdb2pqr --ff=CHARMM --ffout=CHARMM --keep-chain --titration-state-method=propka --with-ph=7.0 \
       protein_model1.pdb protein.pqr
   ```

```bash
pip install pdb2pqr==3.7.1
python example/pbe/04_fagbi27_protocol/prepare_pqr.py 1abc.pdb 1abc.pqr
python example/pbe/04_fagbi27_protocol/prepare_pqr.py 1abc.pdb 1abc.pqr --policy clean
```

The default policy keeps the `ATOM` records as deposited. If PDB2PQR rejects the structure because of alternate
conformers, `--policy clean` removes the input hydrogens and keeps one conformer per residue (highest occupancy,
then altloc A); the paper prepared 1CBN this way. `prepare_pqr.py --verify` rebuilds the 27 PQR files of the
paper from the deposited PDB files and checks that they are byte-identical to the shipped ones.

For a complex, prepare the complex, and cut the receptor and ligand PQR files from the complex PQR by chain, so
that the three share coordinates, charges and radii (example 03).

## 4. Running a calculation

TITAN reads its options from an argument file:

```bash
titan -a arg.in
```

Inside the file each option is written `-option=value`, one per line; on the command line the same option is
written `-option value`. Relative paths are resolved from the working directory, and TITAN writes its outputs
there.

TITAN is a general analysis program whose input is a structure and a topology. For a PB calculation from a PQR
file it is given a neutral four-atom driver structure (`driver.pdb`, `driver.top`, with a minimal force-field
library in `lib_ff/`), and the PQR is the source of the PB problem (`-pbe-test-pqr`). The driver contributes
nothing to the electrostatics. `example/pbe/common` holds these files; every example copies them into its run
directory.

A minimal session, with the paper's protocol:

```bash
mkdir run && cd run
cp ../example/pbe/common/{driver.pdb,driver.top,index.ndx,cal.txt} .
cp -r ../example/pbe/common/lib_ff .
cp /path/to/protein.pqr .
sed 's/{PQR}/protein.pqr/' ../example/pbe/04_fagbi27_protocol/protocol_arg.in.template > arg.in
../titan -a arg.in
```

On Windows, `example/pbe/<example>/run.ps1 -TitanExe <path>` does the same for each example.

## 5. The protocol of the paper, option by option

`example/pbe/04_fagbi27_protocol/protocol_arg.in.template` is the argument file of every protein calculation in
the paper (Figure 5, Supplementary Table S7); `{PQR}` stands for the input file.

| Option | Value | Meaning |
|---|---|---|
| `-f`, `-common_topo`, `-ff`, `-g`, `-1pdb`, `-path-a`, `-path-n` | driver files | The neutral driver structure and its force-field library (section 4). |
| `-pbe-test-pqr` | `protein.pqr` | The PQR file: atom positions, charges (e) and radii (Å). |
| `-pbe` | `bie` | Boundary-integral PB solver (linear PB). |
| `-nt` | `4` | OpenMP threads. |
| `-pbe-int-diel`, `-pbe-ext-diel` | `1`, `80` | Solute and solvent dielectric constants. |
| `-pbe-salt`, `-pbe-temp` | `0.15`, `298.15` | Monovalent salt (M) and temperature (K); the Debye screening applies throughout the solvent outside the SES. |
| `-pbe-probe` | `1.4` | Probe radius (Å) of the SES. |
| `-pbe-bie-primary-surface` | `analytic-ses` | The analytic solvent-excluded surface (`grid-ses` is the alternative grid construction of Section 2.2). |
| `-pbe-bie-analytic-ses-engine` | `patch-ray-resolution` | Builds the analytic SES by casting rays on the grid of spacing h and remeshing on the analytic patches. |
| `-pbe-bie-ses-grid-spacing` | `1.25` | The nominal spacing h (Å). Spacings coarser than the probe radius are not recommended (Section 3.3). |
| `-pbe-bie-ses-charge-refine` | `0.35` | Charge-graded refinement: near a charge at depth d below the surface the target edge is 0.35 d. |
| `-pbe-bie-ses-charge-refine-min-edge` | `0.1` | Smallest target edge (Å) of the refinement. |
| `-pbe-surface` | `vdw` | Grid-solver setting, kept in the protocol file; the BIE surface is the one above. |
| `-pbe-bie-discretization` | `node-patch` | One unknown pair per mesh vertex, on the curved cell around it. |
| `-pbe-bie-near-integration` | `analytic_singular` | Selects the near-field quadrature of Section 2.2: polar Duffy-type quadrature of the self cell, converged to 10^-6 per kernel (Gauss order 8 to 128), and geometrically adaptive subdivision of nearby cells. The name is historical. |
| `-pbe-bie-quad-buffer` | `off` | No quadrature buffer zone beyond the near field. |
| `-pbe-bie-formulation` | `jcp2013-treecode` | The coupled derivative formulation (Juffer et al.; TABI). |
| `-pbe-bie-sum`, `-pbe-bie-fast-operator` | `treecode`, `p1` | Far field by the degree-one Cartesian expansion of each accepted source cell (particle-cluster treecode). |
| `-pbe-bie-tree-leaf-size` | `64` | At most 64 panels per tree leaf. |
| `-pbe-bie-opening-theta` | `0.45` | Opening criterion: a source region is expanded only when its radius divided by its distance from the target is below θ. |
| `-pbe-bie-precond` | `fmm-nearfield-ilu0` | Block ILU(0) on the sparse near-field graph. |
| `-pbe-bie-precond-ilu-neighbor-cap` | `64` | At most 64 blocks per row of the graph. |
| `-pbe-bie-precond-ilu-drop-tol` | `0` | Zero fill: no entry outside the graph is kept. |
| `-pbe-bie-precond-ilu-regularization`, `...-max` | `adaptive`, `0.01` | A diagonal shift is added only if the factor fails its quality gate; every reported protein ran without a shift. |
| `-pbe-bie-precond-memory-budget-mb` | `4096` | Memory budget of the preconditioner (MB). Under a smaller budget a large surface can fall back to a weaker preconditioner; keep it large enough. |
| `-pbe-bie-gmres-restart` | `128` | GMRES restart length. |
| `-pbe-rms`, `-pbe-dphi` | `1e-4` | GMRES stops when the true relative residual is below 10^-4. |
| `-pbe-max-iter` | `720` | GMRES iteration limit. |
| `-pbe-bie-matvec-audit`, `...-sample-count` | `targets`, `512` | After the solve, 512 target panels (1,024 rows) of the accelerated operator are compared with the direct sum on the solution (section 7). |
| `-pbe-timing` | `detail` | Writes the timing and profile tables. |

## 6. Results

TITAN writes its results into the working directory.

| File | Contents |
|---|---|
| `pbe_energy.csv` | `E_rxn_kcalmol` is ΔG_solv (kcal/mol); `iterations_solvent` is the GMRES iteration count; `solute_points` is the number of node patches; `status`. |
| `pbe_timing.csv` | `tabi_true_residual_rel` is the final true relative residual; `tabi_tree_vs_direct_matvec_rel` is the sampled operator discrepancy e_op,S; stage times. |
| `pbe_tabi_profile.csv` | Mesh, preconditioner and solver profile, and the node-patch certificate (`pn_star_normals`, `pn_corner_vertices`, `self_max_order`, `self_fallback`, `near_unresolved`, ...). |
| `pbe_summary.dat` | A readable summary of the run. |
| `pbe_binding_electrostatics.csv` | With `-pbe-binding-cycle=pqr-components`: `Gpol_complex_kJmol`, `Gpol_receptor_kJmol`, `Gpol_ligand_kJmol` and `Delta_Gpol_bind_kJmol` (divide kJ/mol by 4.184 for kcal/mol). |
| `pbe_tabi_surface_solution.csv`, `.vtk` | With `-pbe-bie-surface-output=all`: the surface mesh with, at every node, the boundary unknowns and the total, Coulomb and reaction potentials (`total_potential_kT_e`, kT/e). `example/pbe/02_1vii_surface_potential/view_surface_pymol.py` draws the map in PyMOL on the ±5 kT/e scale of Figure 6; ParaView reads the VTK file directly. |

## 7. Accuracy checks and failures

A converged GMRES solve is not by itself an accurate one. The protocol therefore checks the operator on the
converged solution: 512 target panels of the accelerated operator are recomputed by direct summation with the
same near-field quadrature, and their relative L2 discrepancy is e_op,S (`tabi_tree_vs_direct_matvec_rel`). The
paper accepts a protein solution only if e_op,S ≤ 8 × 10^-3. On 26 of the 27 proteins of the paper e_op,S is at
most 4.4 × 10^-3. 451C fails the gate at h = 1.25 Å because of an oscillating mode on very small cells in a deep
pocket; at h = 0.75 or 0.5 Å the mode is absent and the gate passes (Section 3.4). If a protein fails the gate,
repeat it at a finer spacing.

TITAN stops with an error rather than return a doubtful result when

- the SES fails its watertightness or orientation checks;
- a self cell does not converge under the polar rule up to order 128 and a second rule
  (`curved_self_quadrature_unconverged`);
- a nearby cell cannot be resolved geometrically;
- GMRES does not reach the tolerance within the iteration limit (the energy is then not reported).

The curved-cell certificate in `pbe_tabi_profile.csv` records how often the surface needed a fallback: vertices
whose normal was replaced by the normal of their triangle star, and corner vertices collocated at a smooth point
inside a patch. On the paper's 27 proteins these were at most 3.4% and 1.9% of the vertices.

## 8. Variations of the protocol

- Finer surface: `-pbe-bie-ses-grid-spacing=1.0`, `0.75`, `0.6` or `0.5`. Cost and memory grow roughly with the
  number of node patches; at h = 1.25 Å the protocol is within 0.16% (0.15 M) and 0.19% (10^-5 M) of the finest
  calculation of the paper.
- Degree-four control (p4), used as the reference in the paper: add

  ```text
  -pbe-bie-fmm-local=gmres-p4
  -pbe-bie-fmm-local-refine=operator-tail-static
  -pbe-bie-fmm-local-refine-tol=0.002
  -pbe-bie-fmm-local-cache=green-jet-compressed
  -pbe-bie-fmm-local-cache-budget-mb=4096
  -pbe-bie-fmm-local-fallback=packed-avx2
  -pbe-bie-fmm-local-plan-builder=metric-cache
  -pbe-bie-fmm-local-refine-probe-count=4
  -pbe-bie-fmm-local-refine-probe-mode=rademacher
  -pbe-bie-fmm-local-geometry-guard=strict-geometry-mac
  -pbe-bie-fmm-local-geometry-theta=0.5
  -pbe-bie-taylor-order=4
  ```

  On the analytic SES p4 lowered e_op,S by a factor of 360 to 660 at 1.5 to 2.1 times the p1 wall time.
- Low salt: `-pbe-salt=1e-5` (Debye length about 960 Å), as in the comparison with APBS.
- Other preconditioners, for comparison only: `-pbe-bie-precond=diagonal`; leaf-block Jacobi
  `-pbe-bie-precond=fmm-block-jacobi -pbe-bie-precond-block-panel-cap=64`; restricted additive Schwarz
  `-pbe-bie-precond=fmm-overlap-schwarz -pbe-bie-precond-block-panel-cap=32 -pbe-bie-precond-overlap-panel-cap=64`;
  ILU(0) on a distance graph `-pbe-bie-precond=distance-knn-ilu0`.
- Threads: `-nt=N`.

## 9. Examples

| Directory | What it shows | Time (4 threads) |
|---|---|---|
| `example/pbe/01_1vii_solvation` | ΔG_solv of 1VII with the protocol | seconds |
| `example/pbe/02_1vii_surface_potential` | The same, with the surface potential written for visualization | seconds |
| `example/pbe/03_6lzg_binding_cycle` | Polar binding term of the ACE2-RBD complex 6LZG (Table 6) | minutes |
| `example/pbe/04_fagbi27_protocol` | The 27 proteins of Figure 5, from the deposited PDB files to ΔG_solv | about 10 minutes |

Each directory holds its inputs, its `arg.in` and reference results from this release. In
`04_fagbi27_protocol`, `python run_all.py --titan <path to titan>` runs all 27 proteins and compares them with the
paper's values.

## 10. Reproducing the paper

`validation/rsos_v06` holds the results behind every table and figure of the paper as CSV files, the scripts
that ran and summarized the calculations, the APBS and PAFMPB comparison inputs, and a SHA-256 manifest. Its
README maps each table and figure to its files.
