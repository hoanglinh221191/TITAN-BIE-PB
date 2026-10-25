# Evidence for the revised paper

This directory holds the results behind the tables and figures of "A Stabilized Boundary-Integral
Poisson-Boltzmann Solver for Biomolecular Solvation Electrostatics" (revised version), the scripts that ran and
summarized the calculations, the inputs that are not already under `example/pbe`, and `SHA256SUMS`.

The TITAN results were computed with builds of the revised code (MSYS2 UCRT64 on Windows, four threads, one
workstation); the timed stages ran on the source of this release. Rebuilt from this release, TITAN reproduces the
27 protein energies, iteration counts and node-patch counts of Figure 5 to the printed digits
(`example/pbe/04_fagbi27_protocol/run_all.py`, starting from the deposited PDB files).

## Where each result is

| Paper | Files |
|---|---|
| Tables 3-5, Figures 2-3, Supplementary Tables S8-S9 (preconditioner matrix: five geometries, three surfaces, six preconditioners) | `data/titan/matrix.csv` |
| Figure 4 (surface generators and threads) | `data/titan/xgen.csv`, `data/titan/threads.csv` |
| Figure 5, Supplementary Table S7 (27 proteins) | `data/titan/proteins.csv`, `data/titan/control451c.csv`, `data/apbs/apbs_27_proteins_0.5A.csv`, `data/fagbi/fagbi27_published.csv` |
| Table 6 (polar binding terms) | `data/titan/binding.csv`, `data/apbs/apbs_binding_0.5A.csv` |
| Figure 6 (surface potential of 6LZG) | `data/figure6/6lzg_surface_comparison_summary.csv`, `scripts/figure6/` |
| Section 3.4, cost-accuracy curve | `data/titan/presets.csv`, `data/titan/presets_lowsalt.csv`, `data/titan/presets27.csv`, `data/compare/presets_summary.csv` |
| Supplementary Table S1 | `data/titan/s1.csv` |
| Supplementary Tables S2 and S6 (sphere-patch stress ladder) | `data/titan/s2.csv`, `data/titan/s6.csv` |
| Supplementary Table S3 (ILU(0) neighbour cap) | `data/titan/cap.csv` |
| Supplementary Table S4 (dense calculation) | `data/titan/dense.csv` |
| Supplementary Table S5 (spacing, 1VII) | `data/titan/s5.csv` |
| Supplementary Table S11 (p1 against p4) | `data/titan/p1p4.csv` |
| Supplementary Table S12 (spacing refinement) | `data/titan/spacing.csv` |
| Supplementary Table S13 (coarse correction) | `data/titan/coarse.csv` |
| Supplementary Table S14 (spacing convergence, five proteins, two salt concentrations, p1 and p4) | `data/titan/hsweep.csv`, `data/titan/hsweep_p4.csv`, `data/titan/hsweep_lowsalt.csv`, `data/titan/hsweep_lowsalt_p4.csv` |
| Supplementary Table S15 (cost against APBS) | `data/compare/equal_accuracy_levels.csv`, `data/compare/equal_accuracy_targets.csv`, `data/compare/balanced.csv`, `data/apbs/sweep_*.csv` |
| Supplementary Table S16 (Kirkwood spheres) | `data/titan/kirkwood_R2.csv`, `data/titan/kirkwood_R10.csv` and their `_p4` versions (energies), `data/kirkwood_surface/` (surface potential), `data/pafmpb/kirkwood_pafmpb_off.csv` (PAFMPB on the same nodes) |
| Supplementary Table S17 (PAFMPB on 1SH1 and 1A63) | `data/pafmpb/afmpb_bench.csv` |

Each row of `data/titan/*.csv` is one TITAN run: case, spacing, status, ΔG_solv (`energy_kcalmol`), node
patches, GMRES iterations, true residual, sampled operator discrepancy (`op_audit_rel`, e_op,S) and direct
residual (`direct_residual_rel`, r_dir,S), preconditioner, times, peak memory, and the curved-cell certificate.
The `run_dir` column names the directory on the machine that ran the campaign.

## APBS

APBS 3.4.1 (`mg-auto`, `lpbe`, `bcfl sdh`, `srfm smol` at 0.15 M with ion radius 2.0 Å, `srfm mol` with ion
radius 0 at 10^-5 M, `chgm spl2`, `sdens 10`, `srad 1.4`, `swin 0.3`, dielectric 1/80, 298.15 K) ran on the same
PQR files as TITAN; `scripts/campaign/tabi_rsos_apbs_sweep.py` writes its input files. For the comparison at
10^-5 M, the cavities enclosed by the protein, which TITAN assigns to the solute, were filled in the APBS input
with uncharged atoms of probe radius (`inputs/apbs_cavities_filled/`, made by
`scripts/checks/apbs_cavity_box_check.py`: probe-centre classification on a 0.15 Å grid with 26-connectivity).
`data/apbs/coarse_box_padding_check.csv` is the outer-boundary check (coarse-grid padding 16 and 48 Å per side).

## PAFMPB

PAFMPB (Zhang et al. 2015, AEGB v2.0) was compiled on Windows with GCC and OpenMP in place of the Intel compilers
and Cilk Plus it was released for. `scripts/pafmpb/PATCH_NOTES.txt` lists every change, including two setup
loops run serially to remove data races; `patch_win.py` and `compat/cilk/` make those changes. The PAFMPB source
itself is not redistributed here.

## Scripts

`scripts/campaign/tabi_rsos_campaign.py` defines every stage (its inputs, spacings and options) and runs it with
`tabi_rsos_runner.py`; `run_night_v06_timing.sh` is the sequence that produced the timed results, and
`tabi_rsos_compare.py` and `tabi_rsos_presets_curve.py` derive the comparisons with APBS. The scripts keep the
absolute paths of the machine that ran them: point them at this release and at the inputs here and in
`example/pbe` before running them elsewhere.

## Inputs

- The 27 proteins and 6LZG: `example/pbe/04_fagbi27_protocol/inputs` and `example/pbe/03_6lzg_binding_cycle`.
- `inputs/preconditioner_matrix/`: the five test geometries of the preconditioner matrix (Barnase, 1A63, 1SH1,
  1A7M, 1UBQ; AMBER PQR files and their source structures).
- `inputs/binding_1brs/`: Barnase-Barstar (1BRS), complex and components.
- `inputs/apbs_cavities_filled/`: the APBS inputs with enclosed cavities filled.

The Kirkwood spheres are generated by `tabi_rsos_campaign.py`.
