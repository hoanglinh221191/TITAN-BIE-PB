# TITAN: Analytical FMM/M2L Boundary-Integral Poisson-Boltzmann Solves with Sparse Near-Field ILU(0) Preconditioning

## Abstract

Boundary-integral Poisson-Boltzmann solvers can avoid volumetric grids but require accurate acceleration and robust Krylov preconditioning. We present a guarded fourth-order analytical FMM/M2L path in TITAN, built around analytic singular near-field integration, deterministic far-cell refinement, compressed Green-function jets, packed AVX2 direct fallback, and sparse near-field ILU(0). The method is evaluated as an advanced explicit TABI treecode option, not as a default backend. The primary full-protein panel includes Barnase, the RNA-binding domain of E. coli Rho factor (1A63), and Stichodactyla helianthus neurotoxin I (1SH1). All guarded p4 repeats converged with finite energies, sampled direct-operator residuals below 8e-3, and one plan/cache build. Additional public proteins, Ubiquitin (1UBQ), Chymotrypsin inhibitor 2 (2CI2), and Villin headpiece subdomain (1VII), passed the same supplement gate. A dense Barnase convergence profile shows that sparse near-field ILU(0) resolves a high-coupling case where local block-diagonal preconditioning stalls. The production Poisson-Boltzmann default remains the grid backend, and the default TABI treecode path remains degree-one.

## Methods

TITAN solves a coupled boundary-integral Poisson-Boltzmann system over triangulated molecular surfaces. The accelerated path follows the FMM motivation of Greengard and Rokhlin [2] and the Cartesian boundary-integral PB context of Chen, Tausch, and Geng [1]. GMRES is used for the nonsymmetric coupled system [3]. The guarded p4 configuration is explicit: Taylor order 4, opening theta 0.6, static tail tolerance 0.002, strict-geometry MAC theta 0.5, analytic-singular near integration, compressed degree-9 Green jets, packed AVX2 direct fallback, and sparse near-field ILU(0).

For panel i, the coupled unknown is u_i = (phi_i, partial_n phi_i)^T. The sparse preconditioner keeps a 2 x 2 block on FMM near-source graph edges, including self edges and at most 64 neighbors. ILU(0) keeps the original graph sparsity and applies block forward/back substitution, following the no-fill incomplete-factorization principle [4]. Duffy-adaptive quadrature [5] remains an explicit near-field option but is not selected automatically.

Protein structures are prepared with PDB2PQR/APBS-family tooling [6], [7] and public structure identifiers from the RCSB Protein Data Bank [8]. Large-system verification uses sampled direct-operator RMS and sampled direct true residuals because full direct TABI reaction energies are not practical for the largest rows.

## Results

The primary guarded p4 panel passed all repeatability gates. Supplementary full-protein systems also passed, broadening the public benchmark context without changing the main claim. The Barnase convergence profile is interpreted narrowly as a dense/high-coupling preconditioner contrast.

## Limitations

Grid remains TITAN's production Poisson-Boltzmann default. Guarded p4 is an advanced explicit TABI treecode option and is not selected automatically. Full direct TABI reaction energies are not used as the large-system oracle. Duffy-adaptive near integration is explicit-only. The Barnase preconditioner contrast should not be generalized to all proteins.
