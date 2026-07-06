July 2, 2026

Editors
Journal of Computational Physics

Dear Editors,

We are pleased to submit our manuscript, "TITAN: Guarded Cartesian FMM/M2L Boundary-Integral Poisson-Boltzmann Solves with Sparse Near-Field ILU(0) Preconditioning," for consideration as an original research article in the Journal of Computational Physics.

The manuscript addresses a numerical bottleneck in boundary-integral Poisson-Boltzmann electrostatics. Boundary-integral formulations reduce the problem to the molecular surface, but the resulting coupled potential/normal-derivative operator is dense, derivative-rich, and difficult for Krylov solvers when near-field couplings are not represented in the preconditioner. Our study asks whether fast matrix-vector evaluation, graph-derived sparse near-field ILU(0), and guarded high-order Cartesian FMM/M2L can be combined into an auditable CPU workflow for protein-scale surface discretizations.

The main contributions are:

1. A sparse near-field ILU(0) preconditioner built directly from the FMM near-source graph, with 2 x 2 coupled local operator blocks and no formation of the full dense matrix.
2. A guarded Cartesian p4 FMM/M2L path with compressed Green-function jets, deterministic geometric fallback to direct evaluation for risky far interactions, and packed AVX2 direct-fallback kernels.
3. A verification protocol that separates internal GMRES convergence from agreement with sampled dense direct-operator rows.

Across the six-protein benchmark panel, weakly preconditioned degree-one FMM configurations did not meet the manuscript acceptance criteria, whereas the near-field ILU(0) rows reached the GMRES stopping criterion in 128 reported iterations. The guarded p4 + sparse ILU(0) path retained this convergence behavior while passing sampled direct-operator audits; the largest sampled audit RMS reported in the manuscript is 2.534e-3, and the largest peak resident memory is 369.6 MB. The supplementary material includes p-order, geometric-guard, ILU neighbor-cap, APBS context, small dense-direct calibration, and surface-density stress checks that bound the claims without presenting the method as a universal parameter optimum.

We believe the manuscript is well suited to the Journal of Computational Physics because it presents a numerical method for a physical modeling problem, emphasizes robustness and reproducibility, and reports both algorithmic evidence and failure-boundary diagnostics. The work is intended as a computational physics contribution: it does not claim superiority over GPU implementations or production grid solvers, but instead provides a guarded, auditable boundary-integral pathway for Poisson-Boltzmann solves on molecular surfaces.

This manuscript is original, has not been published previously, and is not under consideration elsewhere. All authors have approved the submission. The source code, benchmark manifests, scripts, figures, and reference outputs corresponding to the reported tables will be released under the MIT License in a versioned public repository and archived with a persistent DOI before final submission.

Thank you for considering our manuscript. We would be grateful for the opportunity to have it reviewed by the Journal of Computational Physics.

Sincerely,

Hoang Linh Nguyen, Ph.D.

Corresponding author

Institute of Fundamental and Applied Sciences, Duy Tan University, Ho Chi Minh City, Vietnam

School of Engineering and Technology, Duy Tan University, Da Nang, Vietnam

Email: nguyenhoanglinh9@duytan.edu.vn
