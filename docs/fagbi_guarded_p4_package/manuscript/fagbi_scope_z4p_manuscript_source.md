# TITAN: Analytical FMM/M2L Boundary-Integral Poisson-Boltzmann Solves with Sparse Near-Field ILU(0) Preconditioning

## Abstract

Boundary-integral Poisson-Boltzmann solvers avoid volumetric grids but require accurate acceleration and robust Krylov preconditioning. We present a guarded fourth-order analytical FMM/M2L path in TITAN that remains an advanced explicit TABI treecode option. The path uses static far-cell refinement, a strict geometric guard, compressed Green-function jets, packed AVX2 direct fallback, and sparse near-field ILU(0). Full prepared proteins were evaluated without atom sampling: Barnase (PyDelPhi benchmark), RNA-binding domain of E. coli Rho factor (1A63), Stichodactyla helianthus neurotoxin I (1SH1). Across three repeats per protein, the guarded path was finite, repeatable, and passed sampled direct-operator checks. Grid-based Poisson-Boltzmann remains TITAN's production default.

## Methods

The guarded p4 path is selected explicitly. It is not the automatic TABI treecode policy. The configuration uses Taylor order 4, opening theta 0.6, static tail tolerance 0.002, strict-geometry MAC theta 0.5, compressed degree-9 Green jets, packed AVX2 direct fallback, analytic-singular near integration, and sparse near-field ILU(0) with 64 neighbors.

For panel i, the coupled unknown is u_i = (phi_i, partial_n phi_i)^T. The preconditioner retains a 2 x 2 block on edges from the FMM near-source graph and applies ILU(0) without fill. The operator is verified with sampled direct-operator RMS and sampled direct true residuals, rather than requiring full direct TABI energies for large systems.

## Results

- Barnase (PyDelPhi benchmark): atoms=1700, panels=14942, repeats=3, audit_max=6.38294E-05, direct_residual_max=6.3829E-05, energy_spread=0, median_matvec_ratio=0.45016110118065.
- RNA-binding domain of E. coli Rho factor (1A63): atoms=2065, panels=17250, repeats=3, audit_max=0.00235702, direct_residual_max=0.00235711, energy_spread=0, median_matvec_ratio=0.674389927487744.
- Stichodactyla helianthus neurotoxin I (1SH1): atoms=696, panels=6673, repeats=3, audit_max=4.89439E-06, direct_residual_max=4.89439E-06, energy_spread=0, median_matvec_ratio=0.383251316898367.

## Limitations

- The guarded p4 path is an advanced explicit option, not a default.
- Grid remains the production Poisson-Boltzmann default.
- Large-system verification uses sampled direct-operator checks.
- Full direct TABI energies are not required for the largest rows.


## Supplementary public-protein expansion

The same guarded p4 gate was repeated on Ubiquitin (1UBQ), Villin headpiece subdomain (1VII), Chymotrypsin inhibitor 2 (2CI2). These rows are treated as supplement evidence because they were added after the primary three-protein panel was locked.
