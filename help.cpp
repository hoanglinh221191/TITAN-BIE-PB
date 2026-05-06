#include "arguments.h"
#include <iostream>

int arguments::help_print() {
  if (help == false)
    return 1;
  else {
    cout << R"MANUAL(
================================================================================
                       TITAN ENGINE - USER MANUAL
                 author:nguyenhoanglinh9@duytan.edu.vn
================================================================================

SYNOPSIS
    titan.exe [OPTIONS]
    titan.exe -a arg.in

  DESCRIPTION
    Titan is a high-performance molecular analysis engine for computing
    electrostatic and van der Waals interactions, GBSA/ALPB solvation,
    hydrogen bonds, DCCM network analysis, PCA/FES, allostery, corr-kNN,
    TICA, MSM,
    ZTS minimum-energy paths, DHAM/MBAR reweighting, TPT pathway analysis,
    LJ condensation, LLPS, GIST water thermodynamics, IFP occupancy,
    LSP-MD, and related observables from PDB/XTC/TRR trajectories.

================================================================================
1. CORE INPUT / OUTPUT
================================================================================

SUPPORTED INPUT MODES

  1) Single structure / single-PDB trajectory
     Use:
       -f <file>
       -1pdb
     Typical case:
       all frames are stored in one PDB/GRO-like text file.

  2) XTC/TRR trajectory with topology reference
     Use:
       -s <structure.pdb>
       -traj <traj.xtc|traj.trr>
     Notes:
       -s sets the reference structure and enables XDR mode.
       Coordinates from XTC/TRR are converted from nm to Angstrom.

  3) Multi-trajectory batch mode
     Use:
       -i <input.txt>
     Notes:
       each line in input.txt is one trajectory file. Titan can process
       multiple replicas or multiple trajectory chunks in one run.

  4) Argument-file mode
     Use:
       -a arg.in
     Notes:
       options can be stored in a text file. For single-token entries,
       arg.in commonly uses key=value syntax, for example:
         -traj=md.xtc
         -ff=AMBER99SB
         -g=1&2

CORE DEPENDENCY CHAIN

  -pca        implies -network
  -allostery  implies -network
  -corr-knn   implies -network
  -tica       implies -network
  -zts        requires -pca or -tica; if both are active, set -zts-domain
  -umap       implies -tica and -network
  -msm        implies -tica and -network
  -dham       implies -msm and -tica; for protein-like traces it also uses
              -network, while LJ jobs switch to the LJ phase-2 feature backend
  -mbar       standalone REMD/HREMD temperature-ladder reweighting;
              domain=tica implies -msm/-tica/-network, domain=pca implies
              -pca/-network, domain=both enables both paths
  -tpt        implies -msm and -tica; phase 0 uses the built MSM TPM +
              stationary distribution directly
  -community=louvain implies -network
  -graph-source=dccm implies -network
  -graph-source=lspmd implies -lspmd
  -graph-source=allostery-gc / allostery-te imply -allostery
  -graph-source=tpt implies -tpt
  -lspmd      runs as its own network-style module and does not need -network

GENERAL OUTPUT BEHAVIOR

  Common files:
    log.txt               Main runtime log
    cal.txt / cal.xvg     Total interaction energy summaries/time-series
    out_*.pdb             Optional coordinate output from -opdb / -opdbsep

  Module-specific output files are created only when the corresponding module
  is enabled. See Chapter 2 for module input/output details.

IMPORTANT NOTES

  - When DHAM is not enabled, or DHAM does not converge, the rest of the
    pipeline still works normally. PCA, FES, basin extraction, TICA, MSM,
    LLPS, GIST, and other active modules continue in their standard mode.

  - Default equilibrium weighting uses `-eq-weight-source auto`:
    if DHAM converges, PCA/FES uses DHAM frame weights at the requested
    target temperature; otherwise TITAN falls back to uniform weights.

  - MBAR v1 is standalone-only. It writes `mbar_*` outputs, reweights to a
    target temperature from REMD/HREMD energy ladders, and does not overwrite
    the standard PCA/TICA/DHAM outputs unless you explicitly opt in with
    `-eq-weight-source mbar`.

  - PCA basin extraction writes an auto-generated shell script:
      extract_basins.sh
    Titan does not check whether GROMACS is installed. Please make sure
    `gmx` is available in your environment before running that script.

  - Representative basin structures are chosen as the frames closest to the
    FES sink in PCA space. They are not guaranteed to be the lowest-energy
    frames inside the basin.

================================================================================
2. MODULE INPUTS / OUTPUTS
================================================================================

2.1 Interaction Energy / Force Analysis

  Purpose:
    Compute electrostatic and van der Waals interaction energies, force maps,
    residue decomposition, COM-distance projections, and related observables.

  Main inputs:
    -f / -s / -traj / -i
    -g <group1&group2>
    -ff <forcefield>
    -topo / -topo2 / -common-topo when needed

  Typical outputs:
    cal.txt
    cal.xvg
    average_emap_pairs.tsv
    average_emap_matrix.tsv
    out_dis_com_energy.dat
    out_dis_com.dat

2.2 GBSA / ALPB

  Purpose:
    Add implicit-solvent free-energy terms on top of interaction analysis.

  Main inputs:
    -gbsa or -alpb
    -igb, -ext_diel, -int_diel, -saltcon, -epsilon_in

  Typical outputs:
    gbsa_snapshots.dat
    gbsa_timeseries_avg.dat
    gbsa_summary.dat
    gbsa_decomposition.dat
    gbsa_average_emap_pairs.tsv
    gbsa_total_energy_maptype.tsv
    out_entropy.dat   (auto-enabled by -gbsa for free-energy completeness)

2.2b PBE

  Purpose:
    Solve Poisson-Boltzmann electrostatics on each processed frame.
    The v1 implementation uses CPU double precision, a two-dielectric model,
    VDW/SAS/SES/DelPhi-like SES/DelPhi-like VDW surface classification,
    Stern-aware linear ion screening, selectable Coulombic/dipolar/zero
    boundaries, linear PB by default, optional nonlinear PB Picard/SOR
    polishing on the fine grid, and reaction-field energy reporting.
    `-pbsa` wraps this path for production two-molecule MM/PBSA: it enables
    PBE deltaG_solv endpoint solves, force-field elec/vdw interaction energy,
    SASA nonpolar energy, optional interaction entropy, and per-frame binding
    energy CSV output. The SASA/SA nonpolar core is shared with GBSA, while
    each workflow keeps its own residue grouping and output files.

  Main inputs:
    -pbe
    -pbe-scale, -pbe-margin, -pbe-gridbox
    -pbe-model, -pbe-biomodel, -pbe-dielectric
    -pbe-int-diel, -pbe-ext-diel, -pbe-salt, -pbe-temp
    -pbe-surface, -pbe-probe, -pbe-boundary, -pbe-charge-assignment,
    -pbe-radius-convention
    -pbe-ion-access, -pbe-ion-radius
    -pbe-dgsolv
    -pbsa
    -pbsa-components, -pbsa-decomp, -pbsa-residue-frames
    -pbsa-ie, -pbsa-ie-block, -pbsa-ie-min-samples
    -pbsa-ie-filter-sigma, -pbsa-perf
    -pbe-max-iter, -pbe-rms, -pbe-dphi, -pbe-omega
    -pbe-nl-max-iter, -pbe-nl-inner-iter, -pbe-nl-stagnation
    -pbe-nl-block, -pbe-nl-dphi, -pbe-nl-relax, -pbe-nl-clip
    -pbe-nl-continuation, -pbe-nl-line-search, -pbe-nl-residual
    -pbe-nl-solver, -pbe-nl-inner, -pbe-nl-multigrid
    -pbe-gaussian-exp, -pbe-gaussian-sigma, -pbe-density-cutoff,
    -pbe-surface-cutoff, -pbe-gap-diel
    -pbe-focus, -pbe-multigrid, -pbe-focus-levels, -pbe-focus-margin
    -pbe-mg-levels, -pbe-mg-pre, -pbe-mg-post, -pbe-mg-coarse-iter
    -pbe-map, -pbe-map-frames, -pbe-max-grid, -pbe-kernel
    hidden PQR parity only: -pbe-parity-dump, -pbe-parity-map-dump,
    -pbe-test-charge-grid, -pbe-test-charge-grid-diff,
    -pbe-test-surface-map,
    -pbe-test-midpoint-epsilon, -pbe-test-midpoint-diff,
    -pbe-test-vdwms-remap,
    -pbe-test-boundary-map, -pbe-test-boundary-convention,
    -pbe-test-field-units, -pbe-test-operator-convention,
    -pbe-test-rhs-map, -pbe-test-screen-map, -pbe-test-denom-map,
    -pbe-test-nl-source, -pbe-test-newton-derivative,
    -pbe-test-energy-convention
    In arg files use -flag=value; on the CLI use -flag value.
    v1 requires normal charge/radius initialization; do not combine with -nocal
    Zero-radius polar hydrogens from supported force fields use a small
    Bondi-style PBE surface fallback radius while charges remain mandatory.

  Typical outputs:
    pbe_energy.csv
    pbe_summary.dat
    mm_pbsa_energy.csv, mm_pbsa_summary.dat, mm_pbsa_by_trajectory.csv,
    mm_pbsa_ie_diagnostics.csv, and mm_pbsa_performance.dat when -pbsa is
    enabled
    optional: mm_pbsa_decomposition.dat and mm_pbsa_residue_*.csv when
    -pbsa-decomp is molecule/residue
    energies are reported in kJ/mol and kcal/mol; kT is retained for continuity
    optional: pbe_phi_*.cube, pbe_epsilon_*.cube, and pbe_ion_*.cube
    optional: pbe_probe_sites.csv when DelPhi focus/probe fixture knobs are used
    hidden PQR parity only: pbe_parity_stage.csv and pbe_parity_points.csv
    barnase PQR forensic mode additionally writes pbe_barnase_forensic_inputs.csv,
    pbe_barnase_surface_audit.csv, pbe_barnase_operator_audit.csv, and
    pbe_atom_phi_audit.csv; ISCM parity with -pbe-parity-dump=all writes
    pbe_iscm_surface_charge_samples.csv and pbe_iscm_energy_terms.csv;
    -pbe-parity-map-dump=samples writes
    pbe_field_samples.csv, pbe_boundary_samples.csv,
    pbe_operator_stage.csv, pbe_charge_grid_audit.csv,
    pbe_nwt_trace.csv, and pbe_nwt_solver_progress.csv for DelPhi NWT
    nonlinear parity rows. The hidden
    DelPhi map extractor script writes
    oracle-side pbe_oracle_field_samples.csv, pbe_oracle_boundary_samples.csv,
    pbe_oracle_boundary_map.csv, pbe_oracle_boundary_inputs.csv,
    pbe_oracle_nwt_trace.csv, pbe_oracle_nwt_solver_summary.csv,
    pbe_oracle_midpoint_epsilon_samples.csv,
    pbe_oracle_atom_radius_audit.csv, a replayable sparse
    pbe_oracle_surface_map.csv with inside semantics, and
    barnase_stage_oracle.tsv under test/tmp for stage comparison. With
    -FullMidpointEpsilon it writes a full replayable midpoint-epsilon oracle
    for one salt only, default 0.15 M, selectable with -FullMidpointSalt.
    With -FullIscmEnergyOracle it also writes
    pbe_oracle_energy_components.csv, pbe_oracle_iscm_full.csv, and
    pbe_oracle_iscm_energy_terms.csv for final DelPhi ISCM convention checks.
    Phase 6d extends the hidden exporter with -FullGaussianOracle,
    -FullRpbOracle, -FullFocusProbeOracle, and -Fixture. The Gaussian density
    oracle path is strict-capable for the 5tif Gaussian-density fixture; RPB
    source comparison is report-first until the remaining grad-epsilon/source
    convention is locked. test/pbe_timing_compare.ps1 writes
    pbe_timing_compare.csv for TITAN vs pyDelPhi wall-clock/profile rows.
    test/pbe_nonlinear_nwt_trace_compare.ps1 compares TITAN NWT rows against
    the pyDelPhi oracle in fixed stage order, and
    test/pbe_nwt_solver_progress_smoke.ps1 compares block RMS/max-dphi progress
    before Gaussian/RPB/focus strict.

2.3 Hydrogen Bonds

  Purpose:
    Detect hydrogen bonds between the selected molecular groups.

  Main inputs:
    -hb
    optional: -hb side-chain, -ahb, -dhb

  Typical outputs:
    stat_hydrogenbonds.dat
    stat_average_HB.dat
    out_list_res_HB_<group>.dat
    out_hydrogen_bonds.dat
    out_hydrogen_bonds_receptor.dat
    out_hydrogen_bonds_ligand.dat

2.3k HBond / Salt / Water Kinetics

  Purpose:
    Track primitive hydrogen bonds, salt bridges, and water bridges as
    time-resolved on/off events, then export residence / lifetime statistics
    plus continuous and intermittent correlation functions.

  Main inputs:
    -hbk
    -hbk-kind <all|hbond|salt|water>
    -hbk-domain <inter|intra|both>
    -hbk-maxlag <N>
    -hbk-hbond-cutoff <A>
    -hbk-angle <deg>
    -hbk-salt-cutoff <A>
    -hbk-water-cutoff <A>

  Typical outputs:
    hbond_kinetics_summary.dat
    hbond_kinetics_hbond_pairs.csv
    hbond_kinetics_salt_pairs.csv
    hbond_kinetics_water_pairs.csv
    hbond_kinetics_hbond_ct.dat
    hbond_kinetics_salt_ct.dat
    hbond_kinetics_water_ct.dat
    hbond_kinetics_hbond_lifetimes.dat
    hbond_kinetics_salt_lifetimes.dat
    hbond_kinetics_water_lifetimes.dat

  Notes:
    - `-hbk` is a new kinetics module. It does not replace legacy `-hb`,
      `-ifp`, or the classic `global_water_bridges_*` outputs.
    - V1 supports `inter`, `intra`, and `both`; the default is `both`.
    - `inter` means molecule 1 <-> molecule 2.
    - `intra` means molecule 1, plus molecule 2 only outside `-nocal`.
    - Primitive keys are not collapsed to residue-pairs:
      hbond = donor / hydrogen / acceptor
      salt = acidic atom / basic atom
      water = water index / solute atom A / solute atom B
    - Water-bridge kinetics are exported only when a water trajectory is
      available. If water is requested but unavailable, TITAN warns and writes
      header-only water files instead of failing the whole run.
    - Canonical syntax:
      arg.in uses `-hbk-kind=water`
      CLI uses `-hbk-kind water`

2.3a DSSP Secondary Structure

  Purpose:
    Compute DSSP-like secondary-structure states per residue and aggregate
    them across trajectories as mean +/- SEM, plus a consensus text string
    suitable for martinize2-style workflows.

  Main inputs:
    -dssp                enable DSSP output
    -no-dssp             disable DSSP output
    -dssp-hbond-cutoff <val>
    protein/protein-user molecule types

  Typical outputs:
    dssp_mol1.dat
    dssp_mol1_content.dat
    dssp_mol2.dat        only outside -nocal and only if molecule 2 is
                         protein/protein-user
    dssp_mol2_content.dat
                         only outside -nocal and only if molecule 2 is
                         protein/protein-user

  Notes:
    - In -nocal mode, DSSP is evaluated only for molecule 1.
    - TITAN auto-skips DSSP if topology lacks a DSSP-compatible protein
      backbone signature, which avoids meaningless output on LJ/monoatomic jobs.
    - For each residue, Titan writes the dominant state, state fractions,
      SEM across input trajectories, and one `martinize2_dssp` consensus
      string built from the highest-population state per residue.
    - The *_content.dat file reports, for each secondary-structure state,
      1. protein-average content
      2. per-residue occupancy +/- SEM
    - Trajectory aggregation is two-stage:
      1. frame fractions are averaged within each input trajectory
      2. Titan reports mean +/- SEM across trajectories

2.3b RMSD / RMSF

  Purpose:
    Compute a PBC-safe trace-atom RMSD time series and per-trace-atom RMSF
    profile after Kabsch alignment.

  Main inputs:
    -rms                enable RMSD + RMSF
    aliases: -rmsd, -rmsf

  Trace atoms:
    atomistic protein/protein-user: CA
    atomistic RNA/DNA: C4' (fallback C4*)
    MARTINI: BB (fallback BB2)

  Typical outputs:
    rmsd.dat            one file per trajectory; indexed if multiple inputs
                        share one output folder
    rmsf.dat            one file per trajectory; RMSD is combined over all
                        active trace atoms, RMSF is split by molecule

  Notes:
    - In -nocal mode, only molecule 1 is analyzed.
    - Outside -nocal, Titan builds one joint fit set from all supported
      trace atoms in molecule 1 + molecule 2.
    - RMSD is reported for the combined fit set.
    - RMSF is reported separately for molecule 1 and molecule 2 in the same
      RMSF file when both contribute trace atoms.
    - If no explicit PBC mode is requested, Titan auto-enables -pbci 0 for
      RMSD/RMSF and logs that decision.

2.3c Rg2

  Purpose:
    Compute a PBC-safe trace-atom radius of gyration time series for the same
    trace atom set used by RMSD/RMSF.

  Main inputs:
    -rg2                enable Rg2 analysis
    -rg2-hist-bins <N>  histogram bin count for rg2_hist.dat. Default: 100

  Outputs:
    rg2.dat             per-frame frame/time_ps/Rg^2/Rg time series
    rg2_summary.dat     mean/std/min/max for Rg^2 and Rg
    rg2_hist.dat        histogram over Rg^2
    rg2_all.dat         aggregated per-frame time series across multi-input runs
    rg2_summary_all.dat per-trajectory summaries plus a global ALL row
    rg2_hist_all.dat    histogram over Rg^2 using all frames from all inputs

  Units:
    - Rg^2 is reported in A^2.
    - Rg is reported in A.

  Notes:
    - In -nocal mode, only molecule 1 is analyzed.
    - Outside -nocal, Titan builds one joint trace set from molecule 1 +
      molecule 2, matching RMSD/RMSF trace selection.
    - If no explicit PBC mode is requested, Titan auto-enables -pbci 0 for
      Rg2 and logs that decision.

2.3d LJ Condensation / Cluster Density

  Purpose:
    Analyze monoatomic or generic Lennard-Jones-like systems through
    cluster detection on molecule 1, then report condensation observables
    such as largest-cluster Rg and radial density profile rho(r).

  Main inputs:
    -lj
    -lj-cut <val>          cluster cutoff in Angstrom; default 5.0
    -lj-min-cluster <N>    minimum size counted as a cluster; default 2
    -lj-rho-bin <val>      rho(r) bin width in Angstrom; default 1.0
    -lj-rho-max <val>      rho(r) max radius in Angstrom; default auto
    -lj-solid-q6 <val>     q6 threshold for solid-like atoms; default 0.40
    -lj-solid-min-cn <N>   minimum CN for dense/solid-like atoms; default 6
    -lj-crystal-tol <val>  q4_bar/q6_bar crystal-class tolerance; default 0.08

  Typical outputs:
    lj_timeseries.dat      one file per trajectory; indexed if multiple inputs
                           share one output folder
    lj_cluster_rho_r.dat   trajectory-averaged radial density profile around
                           the largest-cluster centroid
    lj_features.dat        phase-2 feature matrix used by TICA/MSM for LJ jobs
    lj_summary.dat         trajectory summary with mean +/- SEM observables

  Notes:
    - LJ analysis always reads molecule 1 only, both in -nocal and outside it.
    - All atoms of molecule 1 are used; no protein/RNA trace-atom filter is
      applied.
    - TITAN auto-skips DSSP in LJ/monoatomic jobs when no protein backbone
      signature is found in the topology.
    - Timeseries columns include bulk density, number of clusters,
      largest-cluster size/fraction, Rg, max cluster radius, dense-region
      q4/q6, Lechner-Dellago q4_bar/q6_bar, solid-like fraction, largest
      solid-cluster fraction, and FCC/HCP/BCC/ICO-like fractions.
    - rho(r) is reported in atoms / nm^3 as a function of radius in Angstrom.
    - If -lj is combined with -tica/-msm and no DCCM trace atoms exist, TITAN
      automatically switches TICA input to LJ phase-2 features:
      bulk density, largest-cluster fraction, normalized Rg/Rmax, CN
      statistics, dense-region q4/q6 + q4_bar/q6_bar, solid-like fraction,
      FCC/HCP fractions, and near-core shell enrichments.
    - If no explicit PBC mode is requested, Titan auto-enables -pbci 0 for LJ
      analysis and logs that decision.

  2.3e Structure Factor / Pair Correlation

  Purpose:
    Build streaming RDF / pair-correlation histograms on molecule 1 and derive
    isotropic structure factors S(q) at finalize time.

  Main inputs:
    -sq
    -sq-channel <all|trace|backbone|sticker>
    -sq-rdf-bin <A>
    -sq-rdf-max <A|auto>
    -sq-q-bin <A^-1>
    -sq-q-max <A^-1|auto>

  Typical outputs:
    sq_summary.dat
    sq_rdf_<channel>.dat
    sq_structure_factor_<channel>.dat
    sq_rdf_<channelA>__<channelB>.dat
    sq_structure_factor_<channelA>__<channelB>.dat

  Notes:
    - SQ analysis always reads molecule 1 only.
    - `all` uses all molecule-1 atoms/beads.
    - `trace` reuses the RMS/Rg trace-selection rules.
    - `backbone` uses CA, C4', C4*, BB, BB1, BB2, and BB3 atoms/beads.
    - `sticker` reuses the -llps-sticker bead list.
    - TITAN only exports RDF / S(q) from valid periodic condensed-phase frames.
      Invalid / non-periodic boxes are skipped with a warning.
    - S(q) is derived from RDF using an isotropic Fourier transform with a
      Lorch window to reduce finite-cutoff ripple.
    - If no explicit PBC mode is requested, Titan auto-enables -pbci 0 for SQ
      analysis and logs that decision.

  2.4 DCCM Network

  Purpose:
    Build an aligned CA/backbone trajectory and compute dynamic
    cross-correlation networks.

  Main inputs:
    -network
    a protein-like trajectory with stable residue mapping

  Typical outputs:
    network_edges.csv
    network_nodes.csv
    dccm_equilibrium_context.dat     (when -eq-weight-source mbar resolves)
    allostery_equilibrium_ranking.csv (when -allostery and -eq-weight-source mbar resolve)

  Integration note:
    `-eq-weight-source mbar` only adds equilibrium/context sidecars for DCCM
    and allostery. It does not alter correlation, MI, GC, or TE calculations.

2.4a Community Detection

  Purpose:
    Partition the DCCM graph into residue communities using Louvain modularity
    optimization.

  Main inputs:
    -community=louvain
    implies -network

  Typical outputs:
    network_communities.csv
    network_community_summary.csv
    network_supergraph.csv
    network_community_id_bfactor.pdb
    network_participation_coeff_bfactor.pdb

2.4b Unified Sparse Graph Analytics

  Purpose:
    Run one shared sparse-graph analytics layer on top of internal TITAN graph
    sources from DCCM, LSP-MD, thresholded allostery GC/TE, and TPT.

  Main inputs:
    -graph
    -graph-source <dccm|lspmd|allostery-gc|allostery-te|tpt>
    -graph-community <auto|off|louvain>
    -graph-filter-topk <N>
    -graph-filter-component <all|largest>

  Typical outputs:
    graph_summary.dat
    graph_<source>_edges.csv
    graph_<source>_nodes.csv
    graph_<source>_communities.csv
    graph_<source>_community_summary.csv
    graph_<source>_supergraph.csv
    graph_<source>_bottleneck_edges.csv
    graph_<source>_bottleneck_nodes.csv

  Notes:
    - V1 consumes only in-memory graph data already produced inside TITAN.
      It does not read standalone CSV edge lists.
    - `dccm`, `lspmd`, and `allostery-gc` can run Louvain communities.
    - Directed sources (`allostery-te`, `tpt`) skip community detection in V1
      and export bottleneck-style node/edge rankings instead.
    - `-community=louvain` is treated as a compatibility alias onto the new
      graph core for the DCCM source.
    - Canonical syntax:
      arg.in uses `-graph-source=dccm`
      CLI uses `-graph-source dccm`

2.5 PCA / FES / Basin Extraction

  Purpose:
    Project aligned CA motion onto principal components, build a 2D free
    energy surface from PC1 vs PC2, cluster FES basins, and export
    representative structures metadata.

  Main inputs:
    -pca
    optional: -fes-grid=<N>
    optional: -pca-solver=<full|partial>, -pca-topk=<N>, -pca-max-k-full=<N>
    optional: -pca-symmetry-aware
    optional DHAM weights if -dham converges

  Typical outputs:
    pca_fel.csv
    pca_fes_2d.dat
    cluster_representatives.dat
    extract_basins.sh
    pca_zts_path.csv        (when -zts -zts-domain pca is enabled)

  Notes:
    - Without DHAM, FES is built from uniform frame counts.
    - With converged DHAM, FES and basin populations are reweighted to
      -dham-target-temp automatically.
    - `-pca-symmetry-aware` canonicalizes equivalent chain blocks before
      global PCA alignment. This helps symmetric dimers/multimers whose
      subunits exchange roles across the trajectory.
    - `extract_basins.sh` uses stored trajectory metadata
      (replica id + time in ps) to help recover full-atom basin structures.
    - `-zts` on PCA exports a minimum-energy path between two FES basins.

2.6 Allostery

  Purpose:
    Estimate nonlinear coupling and directed information flow between residues
    using mutual-information-derived Generalized Correlation (GC) and
    Transfer Entropy (TE).

  Main inputs:
    -allostery
    optional: -allostery-tau=<N>, -allostery-bins=<N>
    optional: -allostery-knn=<N>, -allostery-knn-prior=<frame|contact>
    optional: -allostery-knn-contact-cut=<A>
    optional: -allostery-cut=<A>, -allostery-prune-window=<avg|first>
    optional: -allostery-no-te
    optional: -allostery-estimator=<hist|knn|kde>
    optional: -allostery-estimator-k=<N>
    optional: -allostery-kde-bandwidth=<val>

  Typical outputs:
    raw_mi_edges.csv
    network_gc_edges.csv
    network_te_edges.csv

  Notes:
    - TITAN computes pairwise mutual information (MI) internally and exports it
      to `raw_mi_edges.csv` for every candidate pair that survives the pair
      filter.
    - GC is derived from MI using:
        GC = sqrt(1 - exp(-2 * MI))
    - `network_gc_edges.csv` remains a thresholded coupling graph, while
      `raw_mi_edges.csv` stores the unthresholded MI values for the same
      candidate-pair set.
    - `-allostery-estimator=knn` routes through the exact continuous 3D
      corr-kNN backend. `-allostery-estimator-k` controls the estimator
      neighborhood size, while `-allostery-knn` still means pair-pruning top-k.
    - `-allostery-bins` and `-allostery-kde-bandwidth` are ignored in
      knn mode and TITAN logs that decision.

2.6a Mutual Information / Generalized Correlation

  Interpretation:
    - `raw_mi_edges.csv` reports the raw pairwise MI estimated from the
      discretized allostery projections.
    - `network_gc_edges.csv` reports the MI-derived GC transform, which maps
      nonlinear dependence onto the bounded interval [0, 1).
    - `network_te_edges.csv` remains the directed lagged-information export.

  2.6b Corr-kNN

  Purpose:
    Export raw continuous-estimator MI / GC / TE tables from aligned
    residue-trace fluctuations using exact kNN estimators.

  Main inputs:
    -corr-knn
    optional: -corr-knn-k=<N>, -corr-knn-tau=<N>
    optional: -corr-knn-pair-cut=<A>, -corr-knn-pair-topk=<N>
    optional: -corr-knn-pair-prior=<frame|contact>
    optional: -corr-knn-pair-contact-cut=<A>
    optional: -corr-knn-prune-window=<avg|first>
    optional: -corr-knn-no-te

  Typical outputs:
    corr_knn_summary.dat
    corr_knn_mi_pairs.csv
    corr_knn_gc_pairs.csv
    corr_knn_te_pairs.csv
    corr_knn_equilibrium_ranking.csv   (when -eq-weight-source mbar resolves)

  Notes:
    - `-corr-knn` implies `-network` and reuses the aligned DCCM trace cache.
    - MI uses exact KSG type-1 on 3D fluctuation vectors.
    - TE uses continuous conditional-MI (Frenzel-Pompe style) on
      (source[t-tau], target[t], target[t-tau]).
    - GC is derived from MI using the 3D vector form:
        GC = sqrt(1 - exp(-2 * MI / 3))
    - Outputs are raw tables only; thresholded graph policy stays with the
      legacy allostery/network path.

  2.7 TICA

  Purpose:
    Extract slow collective coordinates from time-lagged covariance and
    optionally build TICA-space FES / ZTS pathways.

  Main inputs:
    -tica
    optional: -tica-lag=<N>, -tica-components=<N>
    optional: -tica-backend=<auto|exact_dense|exact_dual|randomized>
    optional: -tica-rank=<N>

  Typical outputs:
    tica_eigenvalues.dat
    tica_eigenvectors.dat
    tica_projections.csv
    tica_fes_2d.dat         (when TICA FES postproc runs)
    tica_cluster_representatives.dat
    tica_zts_path.csv       (when -zts -zts-domain tica is enabled)

  Notes:
    - TICA-space ZTS uses the projection-cache capable backend path.
    - TITAN auto-promotes `-tica-backend exact_dense` when projection caching
      is required by MSM, DHAM, UMAP, or ZTS.

2.8 MSM

  Purpose:
    Build a Markov State Model on top of TICA projections via K-means
    clustering and lagged transition counting.

  Main inputs:
    -msm
    optional:
      -msm-clusters
      -msm-lag
      -msm-nts
      -msm-kmeans-iter
      -msm-seed
      -msm-macrostates
      -msm-no-boundary-split
      -msm-no-pcca
      -msm-no-its
      -msm-its-lags
      -msm-its-solver
      -msm-its-topk
      -msm-its-max-k-full

  Typical outputs:
    msm_tpm.dat
    msm_discrete_traj.dat
    msm_timescales.dat
    msm_centroids.dat
    msm_its_scan.dat
    msm_macrostates.dat
    msm_macrostate_traj.dat
    msm_count_matrix.dat
    msm_stationary.dat
    msm_equilibrium_summary.dat       (when -eq-weight-source mbar resolves)
    msm_its_tpm/

2.8a Validation

  Purpose:
    Run Chapman-Kolmogorov checks, VAMP-2 candidate ranking, and optional
    bootstrap uncertainty on the current TICA -> MSM -> TPT workflow.

  Main inputs:
    -validation
    optional:
      -validation-ck-multiples <m1,m2,...>
      -validation-vamp-components <c1,c2,...>
      -validation-vamp-lags <l1,l2,...>
      -validation-vamp-clusters <k1,k2,...>
      -validation-vamp-folds <N>
      -validation-auto-pick
      -validation-bootstrap <N>
      -validation-bootstrap-seed <N>

  Typical outputs:
    validation_summary.dat
    validation_ck_summary.dat
    validation_ck_states.dat
    validation_vamp2_ranking.csv
    validation_selected_model.dat
    validation_msm_bootstrap_summary.dat
    validation_msm_bootstrap_states.dat
    validation_tpt_bootstrap_summary.dat
    validation_tpt_bootstrap_states.dat

  Notes:
    - `-validation` implies `-msm` and reuses the solved TICA embedding plus
      cached MSM segments.
    - VAMP-2 candidate sweeps only vary component-prefix, MSM lag, and cluster
      count inside the already solved TICA space.
    - `-validation-auto-pick` is opt-in; when off, TITAN only reports ranking
      and leaves the current MSM/TPT pipeline unchanged.
    - `-validation-bootstrap` covers MSM uncertainty and, when `-tpt` is also
      active, TPT uncertainty. Legacy `-tpt-bootstrap` files remain available.

2.9 UMAP

  Purpose:
    Build a nonlinear 2D/3D embedding from cached TICA projections using
    Titan's internal UMAP solver.

  Main inputs:
    -umap
    optional:
      -umap-neighbors=<N>
      -umap-components=<2|3>
      -umap-epochs=<N>
      -umap-min-dist=<val>
      -umap-spread=<val>
      -umap-seed=<N>
      -umap-knn-mode=<exact|approx|auto>
      -umap-knn-benchmark
      -umap-verbose

  Required context:
    UMAP depends on network -> TICA and reuses the same cached TICA
    projection matrix used by MSM/DHAM workflows. Enabling -umap does not
    automatically enable MSM clustering.

  Typical outputs:
    titan_umap_embedding.csv
    umap_fes_2d.dat                  (when -eq-weight-source mbar resolves)
    umap_cluster_representatives.csv (when -eq-weight-source mbar resolves)
    umap_equilibrium_summary.dat     (when -eq-weight-source mbar resolves)

  Integration note:
    `-eq-weight-source mbar` adds equilibrium-only UMAP sidecars without
    changing the embedding fit itself. Default UMAP behavior is unchanged
    unless MBAR is requested explicitly.

2.10 DHAM Reweighting

  Purpose:
    Reweight MSM state populations using per-replica temperatures and energy
    time-series from GROMACS .xvg energy outputs.

  Main inputs:
    -dham
    -dham-temps <temps.txt>
    -dham-ener <energy_list.txt>
    optional:
      -dham-target-temp <K>
      -dham-ensemble <nvt|npt>
      -dham-pressure-bar <bar>
      -dham-max-iter=<N>
      -dham-tol=<val>
      -dham-n-replicas=<N>
      -dham-dt-energy <ps>
      -dham-window-match <auto|off|strict>

  Required context:
    DHAM depends on TICA -> MSM and requires a valid MSM lag.
    Protein-like jobs use network -> TICA -> MSM.
    LJ jobs use the LJ phase-2 feature backend -> TICA -> MSM.
    Energy list files should point to GROMACS .xvg energy trajectories.
    In NPT mode, TITAN combines the first numeric .xvg column with trajectory
    box volumes using H = E_col1 + P*V.
    The pressure term is the ensemble set-point, not instantaneous pressure(t).
    By default, TITAN trims each replica to the common trajectory/energy
    overlap before DHAM (`-dham-window-match auto`).

  Typical outputs:
    dham_stationary_distribution.dat
    dham_free_energy_profile.dat
    dham_unbiased_tpm.dat
    dham_window_audit.dat

  Integration note:
    `-eq-weight-source auto` prefers converged DHAM frame weights for the
    default PCA/TICA FES outputs.
    If DHAM does not converge, Titan prints a warning and PCA/FES falls back
    to uniform weights while other active modules continue normally.
    NVT mode uses energy-based reweighting.
    NPT mode uses an enthalpy-like observable assembled from the aligned .xvg
    column plus trajectory volume metadata.

2.10a MBAR Reweighting

  Purpose:
    Solve temperature-ladder MBAR for REMD/HREMD-like energy series, export
    per-frame weights, and generate standalone weighted FES outputs in TICA
    and/or PCA space at a chosen target temperature.

  Main inputs:
    -mbar
    -mbar-temps <temps.txt>
    -mbar-ener <energy_list.txt>
    optional:
      -mbar-target-temp <K>
      -mbar-domain <tica|pca|both>
      -mbar-ensemble <nvt|npt>
      -mbar-window-match <auto|off|strict>
      -mbar-temp-mode <fixed|per_frame>
      -mbar-temp-col <col>
      -mbar-demux-map <file>
      -mbar-pressure-bar <bar>
      -mbar-dt-energy <ps>
      -mbar-bootstrap <N>
      -mbar-scan-points <N>

  Required context:
    MBAR v1 is temperature-ladder only and expects one GROMACS .xvg
    observable trajectory per replica/window.
    `-mbar-temp-mode fixed` keeps the original ladder-per-file workflow and
    expects `-mbar-temps`.
    `-mbar-temp-mode per_frame` supports raw GROMACS T-REMD-style XVG files
    where temperature changes with time inside each series.
    In `per_frame` mode, `-mbar-temps` is optional: TITAN can infer the
    ladder from replica temperatures directly, or from `-mbar-demux-map`
    plus the XVG temperature column.
    `-mbar-demux-map` is the authoritative state-assignment route for raw
    replica outputs; if it is omitted, TITAN maps each frame to the nearest
    inferred ladder temperature using the XVG temperature column.
    `-mbar-ensemble nvt` treats the XVG observable as energy.
    `-mbar-ensemble npt` treats the XVG observable as enthalpy-like.
    `-mbar-domain tica` reuses the TICA/MSM frame metadata path.
    `-mbar-domain pca` reuses PCA trajectory metadata and the standard
    `pca_fel.csv` projection export.
    `-mbar-domain both` runs both standalone weighted exports.
    MBAR uses DHAM-style window alignment controls, but does not replace DHAM
    automatically.
    `-eq-weight-source mbar` can opt MBAR weights into the default PCA/TICA
    FES path plus equilibrium-only MSM/TPT summaries and UMAP/DCCM/LSP-MD
    context sidecars.
    It does not reweight MSM kinetics, DHAM, or TPT core TPM/pi.
    If a requested projection domain is unavailable, TITAN warns and skips it.

  Typical outputs:
    mbar_frame_weights.dat
    mbar_summary.dat
    mbar_window_audit.dat
    mbar_cv.dat                       (when -mbar-scan-points > 1)
    mbar_tica_fes_2d.dat              (when domain includes tica)
    mbar_tica_cluster_representatives.csv
    mbar_pca_fes_2d.dat               (when domain includes pca)
    mbar_pca_cluster_representatives.csv

  Integration note:
    MBAR outputs are intentionally standalone in v1.
    Standard `tica_fes_2d.dat`, `pca_fes_2d.dat`, and DHAM/TPT outputs keep
    their original behavior unless `-eq-weight-source mbar` is requested.

2.10b TPT

  Purpose:
    Run Transition Path Theory on the built MSM to estimate committors,
    reactive flux edges, and dominant pathways between user-defined source
    and sink microstate basins.

  Main inputs:
    -tpt
    -tpt-source <micro:list>
    -tpt-sink <micro:list>
    optional:
      -tpt-mode <ensemble|kinetic>
      -tpt-top-paths <N>
      -tpt-flux-cut <val>

  Required context:
    TPT phase 0 depends on a completed TICA -> MSM workflow.
    Source/sink basins currently support:
      `micro:0,1-3`
      `macro:0`
      `tica_rect:x1,y1,x2,y2`
      `lj:q6_bar>=0.172`
      `ljpreset:gas`
      `protein:major_macro`
    `lj:` basin specs are only meaningful for LJ jobs. Protein jobs
    normally use `micro:`, `macro:`, or `tica_rect:`.
    Named LJ presets currently include `gas`, `dense`, `liquid`, and `solid`.
    Protein presets currently include `major_macro`, `minor_macro`,
    `low_pc1`, `high_pc1`, `low_pc2`, `high_pc2`,
    `rmsd_low`, `rmsd_high`, `helical`, `strand`, `coil`, `turn`, and `bend`.
    `protein:rmsd_*` presets require `-rms`; DSSP-informed presets require
    `-dssp`.
    Phase 0 uses the MSM TPM + stationary distribution directly.
    In ensemble mode, TITAN can use the converged DHAM-reweighted TPM/pi
    when `-tpt-use-dham yes` (default). If DHAM does not converge, TITAN
    falls back to the built MSM model automatically.
    `-eq-weight-source mbar` only adds an equilibrium overlay summary from
    MBAR populations. It does not change committors, rates, MFPTs, dominant
    paths, or the TPT model source.
    In kinetic mode, TITAN always uses the built MSM model and writes
    `tpt_kinetics.dat` diagnostics; it still does not certify kinetic
    validity automatically. If `-tpt-demux-map <file>` is supplied, TITAN
    rebuilds a kinetic MSM from the demuxed microstate segments before TPT.
    Kinetic diagnostics now include CK / Markov sanity metrics from the shared
    validation core.
    `-tpt-bootstrap <N>` enables bootstrap uncertainty for both MSM-backed and
    DHAM-backed TPT models. When `-validation-bootstrap` matches the legacy
    TPT bootstrap config, TITAN computes once and writes both file sets.

  Typical outputs:
    tpt_summary.dat
    tpt_committor.dat
    tpt_flux_edges.dat
    tpt_paths.dat
    tpt_state_labels.dat
    tpt_kinetics.dat
    tpt_reactive_fes_2d.dat
    tpt_flux_map_2d.dat
    tpt_bootstrap_summary.dat
    tpt_bootstrap_states.dat
    tpt_all_pairs_kinetics.dat
    tpt_path_clusters.dat
    tpt_bottleneck_edges.dat
    tpt_bottleneck_states.dat
    tpt_equilibrium_overlay.dat       (when -eq-weight-source mbar resolves)
    `tpt_summary.dat` now appends a compact top-bottleneck table.
    `tpt_state_labels.dat` always contains generic state columns.
    For LJ jobs, TITAN appends LJ microstate-mean order parameters.
    If TICA/MSM centroids are at least 2D, TITAN also exports a smoothed
    reactive free-energy surface and a 2D centroid-to-centroid net-flux map.
    Bottleneck tables rank reactive edges and microstates by flux throughput
    near the transition region; dominant-path support is reported as an
    auxiliary signal.

2.11 LLPS

  Purpose:
    Detect dense and dilute phases, percolation, spatial clustering, and
    time-autocorrelation for condensate-like systems.
    The current LLPS stack also covers contact dynamics, dense-phase mobility,
    protein-RNA coupling, RNA compaction/information metrics, and RNA topology
    signatures such as bridges and clamps. Paper-facing sidecars include
    sticker hotspot ranking, droplet shape anisotropy/sphericity proxies,
    droplet genealogy summaries, and FRAP-like residence/recovery proxies.
    `-llps-paper` is an opt-in Martini protein/RNA paper bundle that reuses
    DBSCAN dense labels to export phase partitioning, composition profiles,
    convergence checks, and manuscript-friendly summaries without changing
    default `-llps` output volume.
    Run metadata and diagnostics are written into the central TITAN log
    (`log.txt`) under LLPS Run Context / Output Audit sections; TITAN does not
    create a separate LLPS manifest sidecar.

  Main inputs:
    -llps
    optional:
      -llps-rcut
      -llps-perc
      -llps-fuzzy
      -llps-dynamics
      -llps-dynamics-detailed
      -llps-maxlag
      -llps-dimer-size
      -llps-contacts
      -llps-sticker
      -llps-voxel
      -llps-dbscan-eps
      -llps-dbscan-minpts
      -llps-lineage-overlap
      -llps-prna-rcut
      -llps-msd-fit-start-ns
      -llps-msd-fit-end-ns
      -llps-rna-compaction
      -llps-rna-information
      -llps-rna-topology
      -llps-rna-topology-detailed
      -llps-rna-intra-rcut
      -llps-rna-bridge-geom-min
      -llps-rna-patch-gap-max
      -llps-rna-clamp-require-loop
      -llps-rna-info-rg-bins
      -llps-rna-info-longrange-bins
      -llps-rna-longrange-min
      -llps-energy-correlation
      -llps-energy-detailed
      -llps-paper
      -llps-paper-profile-bin
      -llps-paper-blocks
      -llps-paper-geometry

  Typical outputs:
    llps_timeseries.dat
    llps_autocorr.dat
    llps_degree_dist.dat
    llps_contact_map.dat
    llps_sticker_hotspots.csv
    llps_performance_report.dat
    llps_summary.dat
    llps_density_protein.dx
    llps_density_rna.dx
    polymer_scaling.dat
    llps_dbscan_droplets.pdb
    llps_dbscan_timeseries.csv
    llps_droplet_shape_timeseries.csv
    llps_membership_autocorr_protein.dat
    llps_membership_autocorr_rna.dat
    llps_residence_times.csv
    llps_entry_exit_rates.dat
    llps_cluster_lineage.csv
    llps_fusion_fission.dat
    llps_droplet_genealogy_summary.dat
    llps_frap_recovery_proxy.dat
    llps_sticker_contact_lifetimes.dat
    llps_protein_rna_contact_lifetimes.dat
    llps_neighbor_rewiring.dat
    llps_interface_turnover.dat
    llps_dense_msd_protein.dat
    llps_dense_msd_rna.dat
    llps_liquidity_summary.dat
    llps_rna_compaction_binding_summary.dat
    llps_rna_compaction_binding_binned.csv
    llps_rna_chain_compaction_summary.csv
    llps_rna_intrachain_contacts_summary.csv
    llps_rna_compaction_trace.csv (with -llps-dynamics-detailed=true)
    llps_rna_energy_compaction_summary.dat
    llps_rna_energy_compaction_binned.csv
    llps_rna_binding_energy_trace.csv
    llps_prna_pair_energy_trace.csv (with -llps-energy-detailed=true)
    llps_rna_information_summary.dat
    llps_rna_information_per_chain.csv
    llps_rna_compaction_state_histogram.csv
    llps_rna_partner_state_histogram.csv
    llps_rna_transition_information.csv
    llps_rna_binding_compaction_information.csv
    llps_rna_partner_state_trace.csv (with -llps-dynamics-detailed=true)
    llps_rna_topology_summary.dat
    llps_rna_topology_per_chain.csv
    llps_rna_topology_binding_binned.csv
    llps_rna_topology_state_histogram.csv
    llps_rna_bridge_lifetimes.dat
    llps_rna_clamp_lifetimes.dat
    llps_rna_loop_closure_efficiency.dat
    llps_rna_local_competition_summary.dat
    llps_rna_topology_trace.csv
    llps_rna_bridge_trace.csv
    llps_rna_clamp_trace.csv
    with -llps-paper:
      llps_phase_partition_summary.dat
      llps_radial_profiles.csv or llps_slab_profiles.csv
      llps_composition_timeseries.csv
      llps_contact_taxonomy.csv
      llps_contact_taxonomy_summary.dat
      llps_convergence_report.dat
      llps_paper_summary.md
      llps_paper_summary.json
    log.txt sections:
      LLPS Run Context
      LLPS Macro Summary
      LLPS Output Audit
      LLPS DBSCAN Output Audit
      LLPS Dynamics Output Audit
      LLPS Paper Output Audit

2.12 GIST

  Purpose:
    Map solvent occupancy, energy, and entropy onto a 3D grid around the
    solute to identify favorable and unfavorable hydration regions.

  Main inputs:
    -gist
    optional:
      -gist-resolution=<A>
      -gist-padding=<A>
      -gist-temperature=<K>
      -gist-water-model=<name>
      -gist-min-count=<N>
      -gist-cutoff=<A>

  Typical outputs:
    gist_*.dx
    (Titan exports a set of DX maps using the prefix `gist`)

2.12 IFP (Interaction Fingerprint)

  Purpose:
    Compute per-contact occupancy fingerprints for H-bond, salt bridge,
    hydrophobic, pi-pi, cation-pi, and halogen contacts.
    Supports protein-protein and protein-RNA/DNA systems while keeping
    ifp_occupancy.csv backward compatible.

  Main inputs:
    -ifp
    optional:
      -ifp-hbond
      -ifp-salt
      -ifp-hphob
      -ifp-pipi
      -ifp-cpi
      -ifp-hal
      -ifp-out <file>

  Typical outputs:
    ifp_occupancy.csv
    ifp_rna_annotations.csv   (only when nucleic-acid contacts are detected)
    nucleic_acid_interactions.csv
      coarse RNA/DNA recognition summary, separate from residue-pair IFP

2.12a AWH RC-Binned IFP

  Purpose:
    For GROMACS AWH trajectories, read walker-aligned RC time series from
    .xvg, bin snapshots along the biased reaction coordinate, then report
    IFP patterns and representative structures along that coordinate.
    This workflow is intended for structural mechanism mapping on a biased
    ensemble, not for absolute kinetics or equilibrium populations.

  Main inputs:
    -awh-ifp
    -awh-xvg <filelist>
    optional:
      -awh-xvg-col <N>
      -awh-bins <N>
      -awh-out-prefix <prefix>

  Typical outputs:
    awh_rc_bin_summary.dat
    awh_walker_diagnostics.dat
    awh_rc_overlap_diagnostics.dat
    awh_walker_stability_by_bin.dat
    awh_walker_pair_overlap.dat
    awh_ifp_by_bin.csv
    awh_ifp_heatmap.dat
    awh_bin_structural_summary.dat
    awh_bin_top_contacts.csv
    awh_rc_frame_map.dat
    awh_frame_assignments.dat
    awh_exploratory_context.dat
    awh_rep_bin_*.pdb

2.13 LSP-MD

  Purpose:
    Build a Local Spatial Pattern network from residue geometry without
    requiring the full DCCM workflow.

  Main inputs:
    -lspmd
    optional: -lspmd_cut=<A>

  Typical outputs:
    lspmd_edges.csv
    lspmd_degree_centrality.csv
    lspmd_betweenness_centrality.csv
    dccm_equilibrium_context.dat     (when -eq-weight-source mbar resolves)
    lspmd_equilibrium_ranking.csv    (when -eq-weight-source mbar resolves)

  Integration note:
    `-eq-weight-source mbar` only adds contextual ranking/overlay sidecars.
    TITAN does not rewrite core LSP-MD graph construction or centrality values.

2.14 Other Modules

  CAS:
    Alanine-scanning style hotspot analysis.
    Outputs include:
      cas_receptor.dat
      cas_ligand.dat
      cas_hotspots_b_factor.pdb

  TDA:
    Topological data analysis of structure/geometry.

  ML:
    Simple ML affinity prediction.
    Outputs include:
      titan_ml_features.csv
      titan_rf_prediction.txt

================================================================================
3. EXAMPLES
================================================================================

3.1 Basic energy decomposition from a PDB

  titan.exe -f complex.pdb -g 1&2 -ff AMBER99SB

3.2 GBSA / ALPB free-energy style analysis

  titan.exe -f complex.pdb -g 1&2 -gbsa -igb 5

  titan.exe -f complex.pdb -g 1&2 -alpb -igb 5

3.3 Single XTC trajectory with topology reference

  titan.exe -s reference.pdb -traj md.xtc -g 1&2 -ff AMBER99SB

3.4 Multi-trajectory batch mode from input list

  titan.exe -s reference.pdb -i input.txt -g 1&2 -ff AMBER99SB

  Example input.txt:
    rep0.xtc
    rep1.xtc
    rep2.xtc

3.5 DCCM network + PCA/FES

  titan.exe -s reference.pdb -traj md.xtc -network -community=louvain \
      -pca -fes-grid=80

3.6 Allostery analysis

  titan.exe -s reference.pdb -traj md.xtc -allostery \
      -allostery-tau=10 -allostery-bins=48

  titan.exe -s reference.pdb -traj md.xtc -allostery \
      -allostery-cut=12 -allostery-prune-window=avg -allostery-no-te

  titan.exe -s reference.pdb -traj md.xtc -allostery \
      -allostery-cut=12 -allostery-knn=4

  titan.exe -s reference.pdb -traj md.xtc -allostery \
      -allostery-knn=4 -allostery-knn-prior=contact \
      -allostery-knn-contact-cut=8 -allostery-no-te

  titan.exe -s reference.pdb -traj md.xtc -allostery \
      -allostery-estimator=knn -allostery-estimator-k=4 \
      -allostery-knn=6 -allostery-cut=12

  Inspect:
    raw_mi_edges.csv for unthresholded mutual information
    network_gc_edges.csv for the derived GC network
    network_te_edges.csv for directed transfer entropy

  3.6a Corr-kNN raw information export

  titan.exe -s reference.pdb -traj md.xtc -corr-knn \
      -corr-knn-k=4 -corr-knn-tau=10 -corr-knn-pair-cut=12

  titan.exe -s reference.pdb -traj md.xtc -corr-knn \
      -corr-knn-pair-topk=6 -corr-knn-pair-prior=contact \
      -corr-knn-pair-contact-cut=8 -corr-knn-no-te

  3.6b HBond kinetics

  CLI:
    titan.exe -s reference.pdb -traj md.xtc -hbk \
        -hbk-kind hbond -hbk-kind salt -hbk-domain both -hbk-maxlag 250

  arg.in:
    -hbk
    -hbk-kind=hbond
    -hbk-kind=salt
    -hbk-domain=both
    -hbk-maxlag=250
    -hbk-hbond-cutoff=3.5
    -hbk-angle=135
    -hbk-salt-cutoff=4.0

  3.7 TICA + MSM workflow

  titan.exe -s reference.pdb -i input.txt -tica -tica-lag=10 \
      -tica-components=5 -msm -msm-clusters 150 -msm-lag 10

3.8 DHAM + weighted FES from REMD-style data

  titan.exe -s reference.pdb -i xtc_list.txt -network -tica -msm -dham \
      -msm-clusters 150 -msm-lag 10 \
      -dham-temps temps.txt \
      -dham-ener energy_xvg_list.txt \
      -dham-target-temp 300 \
      -dham-max-iter=5000 \
      -dham-tol=1e-6

  Example temps.txt:
    300
    320
    340
    360

  Example energy_xvg_list.txt:
    rep0_energy.xvg
    rep1_energy.xvg
    rep2_energy.xvg
    rep3_energy.xvg

 Result:
    if DHAM converges, PCA/FES is reweighted automatically to the requested
    target temperature. Basin extraction still works even if DHAM is off.

3.8a Protein NPT DHAM workflow

  titan.exe -s reference.pdb -i xtc_list.txt -network -tica -msm -dham \
      -dham-temps temps.txt \
      -dham-ener potential_xvg_list.txt \
      -dham-target-temp 300 \
      -dham-ensemble npt \
      -dham-pressure-bar 1.0

  Notes:
    Use the ensemble set pressure with -dham-pressure-bar.
    TITAN does not need instantaneous pressure(t) extracted from .edr.
    In NPT mode, TITAN builds H = E_col1 + P*V from the .xvg column and
    per-frame trajectory box volumes.

3.8b LJ REMD / Argon DHAM workflow

  titan.exe -s conf.pdb -i xtc_list.txt -nocal -lj -tica -msm -dham \
      -common-topo topol.top \
      -dham-temps temps.txt \
      -dham-ener energy_xvg_list.txt \
      -dham-target-temp 94 \
      -dham-ensemble npt \
      -dham-pressure-bar 1.0

  Notes:
    Each trajectory file is treated as one replica.
    TITAN uses LJ phase-2 features instead of DCCM trace atoms.
    `-dham-ensemble nvt` keeps the REMD-style energy workflow.
    `-dham-ensemble npt` uses H = E_col1 + P*V with volume taken from the
    trajectory box at each frame.

3.9 MBAR standalone weighted TICA/PCA FES from REMD/HREMD data

  titan.exe -s reference.pdb -i xtc_list.txt -network -tica -msm -pca -mbar \
      -mbar-temps temps.txt \
      -mbar-ener energy_xvg_list.txt \
      -mbar-target-temp 300 \
      -mbar-domain both \
      -mbar-temp-mode per_frame \
      -mbar-window-match auto \
      -mbar-scan-points 24

  Notes:
    MBAR writes separate `mbar_*` outputs and keeps the default TICA/PCA FES
    files untouched.
    Raw GROMACS T-REMD NVT/NPT energy files with a fluctuating temperature
    column can be handled via `-mbar-temp-mode per_frame`.
    If you have a demux/state-assignment map, prefer `-mbar-demux-map`;
    otherwise TITAN falls back to nearest-ladder assignment from the XVG
    temperature column.
    `-mbar-domain tica` is the lightest-weight integration path.
    `-mbar-domain pca` reuses PCA projection exports and is useful when you
    want a target-temperature PCA free-energy map without changing the default
    PCA analysis.

3.9 PCA basin extraction after FES clustering

  titan.exe -s reference.pdb -traj md.xtc -pca

  After the run:
    1. inspect cluster_representatives.dat
    2. make sure GROMACS `gmx` is available in your shell
    3. run:
         bash extract_basins.sh

  PCA ZTS example:
    titan.exe -s reference.pdb -traj md.xtc -pca -zts -zts-domain pca \
        -zts-images 64 -zts-source 0 -zts-sink 1

  TICA ZTS example:
    titan.exe -s reference.pdb -traj md.xtc -tica -zts -zts-domain tica \
        -zts-images 64 -zts-source 0 -zts-sink 1

3.10 LLPS analysis (Martini)

  titan.exe -s system.pdb -traj md.xtc -ff MARTINI -llps

  titan.exe -s system.pdb -traj md.xtc -ff MARTINI -llps \
      -llps-dbscan-eps 12 -llps-dbscan-minpts 8 -llps-voxel 0.6 -llps-fuzzy

3.11 GIST water thermodynamics

  titan.exe -s complex.pdb -traj md.xtc -gist \
      -gist-resolution=0.5 \
      -gist-padding=8.0 \
      -gist-temperature=300 \
      -gist-water-model=TIP3P

3.12 IFP occupancy analysis

  titan.exe -s complex.pdb -traj md.xtc -ifp \
      -ifp-hbond 3.5 -ifp-salt 4.0 -ifp-hphob 4.0

  Protein-RNA/DNA runs keep ifp_occupancy.csv unchanged and may also
  write ifp_rna_annotations.csv with base/sugar/phosphate contact bias.

3.13 LSP-MD network analysis

  titan.exe -s reference.pdb -traj md.xtc -lspmd -lspmd_cut=12.0

3.14 Using arg.in

  titan.exe -a arg.in

  Example arg.in:
    -s=reference.pdb
    -traj=md.xtc
    -g=1&2
    -ff=AMBER99SB
    -network
    -pca
    -tica
    -tica-lag=10
    -msm
    -msm-lag=10

================================================================================
4. OPTION REFERENCE
================================================================================

4.1 Input / Output

  -f <file>             Single PDB/GRO input file (reference structure).
                        Default: complex.pdb

  -s <file>             Structure file for XTC/TRR mode.
                        Sets the reference structure and enables XDR mode.

  -traj <file>          XTC or TRR trajectory file.
                        Auto-detects format from extension.

  -i <file>             List file with multiple trajectory filenames.
                        Default: input.txt

  -1pdb                 Use single-PDB trajectory mode.

  -a <file>             Read arguments from file. Default: arg.in

  -opdb=<name.pdb>      Write output PDB.
  -opdbsep=<name.pdb>   Write per-trajectory/separated output PDB.

  -n <file>             Index file. Default: index.ndx

  -topo <file>          Topology for molecule 1. Supports .itp/.top/.parm7
  -topo2 <file>         Topology for molecule 2
  -common-topo <file>   Use one topology for both molecules
  -common_topo <file>   Alias of -common-topo

  -ff <name>            Force field. Default: AMBER99SB
                        Supported:
                        AMBER99SB
                        AMBER99SB-ILDN
                        AMBER99SB-BSC1
                        AMBER99SB-DISP
                        AMBER14SB
                        AMBER19SB
                        GROMOS43A1
                        GROMOS53A6
                        GROMOS54A7
                        CHARMM27
                        CHARMM36
                        OPLSAA
                        MARTINI
  -fflib-auto on|off    Auto-assign charge/VDW from lib_ff when no topology is
                        provided. Default: on
  -fflib-missing warn|strict
                        Missing lib_ff atom/template policy. Default: warn
  -fflib-hydrogens reference|none
                        Add missing hydrogens for PDB/reference inputs only.
                        Default: reference

  -pattern-input <pat> [b] [e] [step]
                        Pattern input mode; <pat> must contain '*'

4.2 Groups / Molecules / Time

  -g <G1>&<G2>          Interaction groups
  -sg <SG1>-<SG2>       CLI sub-groups within main groups
  -sg= <SG1>-<SG2>      arg.in sub-groups within main groups
  -sg-mode=auto|stream|cache
                        Sub-group energy mode. stream avoids dense pair cache.
  -sg-cache=auto|off|on Dense pair-cache preference for -sg
  -sg-memory-gb=auto|N  RAM budget for -sg auto mode
  -sg-timeseries=auto|full|summary-only
                        Control -sg time-series output volume
  -sg-kernel=auto|scalar|avx2
                        Streaming -sg compute kernel preference
  -sg-cutoff=auto|VALUE|none
                        Streaming -sg cutoff in nm; auto uses energy cutoff
  -mol-1-type <type>    protein, complex, lig, rna, dna
  -mol-2-type <type>    protein, protein-user, complex, lig, rna, dna
  -mol2 <N>             Starting residue/atom index for molecule 2
  -n_wat <N>            Number of waters to track; -1 disables

  -b <time>             Begin time (ps)
  -e <time>             End time (ps)
  -dt <N>               Frame step interval
  -nt <N>               Number of threads

4.3 Energy / Solvation / Entropy

  -gbsa                 Enable GBSA solvation analysis
  -alpb                 Enable ALPB and imply -gbsa
  -igb <N>              GB model. Default: 5
  -ext_diel <val>       External dielectric. Default: 78.5
  -int_diel <val>       Internal dielectric. Default: 1.0
  -saltcon <val>        Salt concentration (M). Default: 0.15
  -epsilon_in <val>     Internal epsilon for solvation models

  -pbe                  Enable frame-based PB electrostatics
  -pbe-model <mode>     linear or nonlinear. Default: linear
                        nonlinear uses the current linear off/focus/vcycle
                        solve as a seed, then runs Picard/SOR or damped Newton
                        on the final fine grid. Vacuum endpoints remain linear.
  -pbe-biomodel <mode>  pbe or rpbe. Default: pbe
  -pbe-dielectric <m>   two-dielectric or gaussian. Default: two-dielectric
  -pbe-scale <val>      Grid scale in points/A. Default: 2.0
  -pbe-margin <val>     Solvent margin around VDW surface (A). Default: 15.0
  -pbe-gridbox <mode>   cuboidal or cubic. Default: cuboidal
  -pbe-int-diel <val>   PBE internal dielectric. Default: 1.0
  -pbe-ext-diel <val>   PBE external dielectric. Default: 80.0
  -pbe-salt <val>       PBE salt concentration (M). Default: 0.15
  -pbe-temp <val>       PBE temperature (K). Default: 298.15
  -pbe-surface <mode>   vdw, sas, ses, delphi-ses, delphi-vdw, gaussian, or
                        gaussian-cutoff. Default: vdw
                        ses uses probe-excluded voxels plus boundary flood fill
                        to identify solvent-accessible grid points.
                        delphi-ses keeps the SES flood fill and refines
                        dielectric faces by 2x2 midpoint-face sampling.
                        delphi-vdw uses DelPhi-compatible discrete media
                        indices at positive x/y/z half-grid faces: each
                        midpoint is classified by VDW atom radius, then the
                        PQR DelPhi-parity path applies a deterministic
                        VDW-to-MS boundary remap before producing final
                        midpoint epsilon. VDW media assignment uses explicit
                        PQR/SIZ radii without small-atom floors; Bondi floors
                        remain a fallback for non-PQR delphi-bondi paths.
                        DelPhi-parity diagnostics append phase RF energies
                        (solvent/vacuum), atom/grid corrected RF energy, and
                        analytic Coulomb energy to pbe_energy.csv.
  -pbe-boundary <mode>  coulombic, dipolar, or zero. Default: coulombic
                        DelPhi parity harnesses can switch boundary mode
                        without changing legacy default behavior.
  -pbe-charge-assignment <mode>
                        trilinear, delphi, or gaussian. Default: trilinear
                        delphi enables strict-grid charge diagnostics for
                        hidden PQR/DelPhi regression fixtures.
  -pbe-radius-convention <mode>
                        auto, native, or delphi-bondi. Default: auto
  -pbe-probe <val>      SAS/SES probe radius (A). Default: 1.4
  -pbe-ion-access <m>   dielectric or stern. Default: dielectric
                        stern applies Debye screening only in ion-accessible
                        solvent after excluding atom radius + ion radius.
  -pbe-ion-radius <val> Ion probe radius for Stern accessibility (A). Default: 2.0
  -pbe-dgsolv [bool]    For two-molecule systems, compute deltaG_solv by
                        six solves: complex, mol1, mol2 each solvent/vacuum.
                        Arg file: -pbe-dgsolv=true; CLI: -pbe-dgsolv
  -pbsa [bool]          Enable MM/PBSA summary for two-molecule systems.
                        Implies -pbe and -pbe-dgsolv; writes
                        mm_pbsa_energy.csv and mm_pbsa_summary.dat with PB
                        endpoint polar, SASA nonpolar, elec/vdw,
                        interaction entropy, and final dG.
                        Arg file: -pbsa=true; CLI: -pbsa true
  -pbsa-components <m>  auto or complex,receptor,ligand. Default: auto
  -pbsa-decomp <mode>   off, molecule, or residue. Default: off
                        residue writes component decomposition without
                        re-solving PB per residue.
  -pbsa-ie <bool>       Enable interaction entropy correction. Default: true
                        Uses the same Duan-Zhang interaction entropy core as
                        TITAN's -entropy module, with PBSA diagnostics.
                        Arg file: -pbsa-ie=true; CLI: -pbsa-ie true
  -pbsa-residue-frames <mode>
                        off or all. Default: off. all writes
                        mm_pbsa_residue_frame.csv.
  -pbsa-ie-block <N>    auto or positive block size. Default: auto
  -pbsa-ie-min-samples <N>
                        Minimum samples before IE correction is applied.
                        Default: 10
  -pbsa-ie-filter-sigma <val>
                        Sigma filter for IE samples. Default: 2.0
  -pbsa-perf <mode>     off, summary, or detail. Default: summary
  -pbe-max-iter <N>     SOR iteration cap. Default: 10000
  -pbe-rms <val>        SOR RMS convergence tolerance. Default: 1e-6
  -pbe-dphi <val>       SOR max delta-phi tolerance. Default: 1e-5
  -pbe-omega <val>      auto or SOR relaxation value in (0,2). Default: auto
  -pbe-nl-max-iter <N>  Nonlinear Picard/SOR sweep cap or Newton work budget.
                        Default: 80
  -pbe-nl-inner-iter <N>
                        Newton correction SOR sweeps or V-cycle cap per outer
                        step.
                        Default: 40
  -pbe-nl-stagnation <N>
                        Newton residual stagnation window. Default: 8
  -pbe-nl-block <N>     Nonlinear diagnostic block size. Default: 5
                        In DelPhi NWT parity rows this is the true block-stop
                        size used for RMS/max-dphi checks.
  -pbe-nl-dphi <val>    Nonlinear max delta-phi tolerance. Default: 1e-4
  -pbe-nl-relax <val>   auto or Picard damping in [0.05,1.0]. Default: auto
  -pbe-nl-clip <val>    Clip alpha*phi before sinh. Default: 20.0
  -pbe-nl-continuation <mode>
                        off or auto. Default: auto
  -pbe-nl-line-search <mode>
                        off or auto. Default: auto
  -pbe-nl-residual <val>
                        Nonlinear residual tolerance. Default: 1e-5
  -pbe-nl-solver <mode> picard or newton. Default: picard
  -pbe-nl-inner <mode>  auto, sor, or vcycle. Default: auto
                        newton solves a stencil Jacobian correction with SOR
                        or geometric V-cycle and appends correction residual,
                        inner levels/cycles, Jacobian, lambda, and stop-reason
                        diagnostics. In DelPhi PQR parity rows, auto keeps the
                        strict Newton path on SOR first; if line search rejects
                        a correction, TITAN falls back to a DelPhi NWT SOR
                        block with denominator += screen*cosh(phi) and RHS
                        += screen*(phi*cosh(phi)-sinh(phi)). PQR parity dump
                        mode writes pbe_nwt_trace.csv and
                        pbe_nwt_solver_progress.csv; pbe_energy.csv appends
                        nonlinear stress, osmotic, and total correction columns.
  -pbe-nl-multigrid <m> off or fas. Default: off
  -pbe-gaussian-exp <v> Gaussian density exponent. Default: 2
  -pbe-gaussian-sigma <v>
                        Gaussian density sigma scale. Default: 1.0
  -pbe-density-cutoff <v>
                        Gaussian cutoff density. Default: 0.759
                        Gaussian dielectric now uses DelPhi-style
                        multiplicative density union at centers and positive
                        midpoint faces instead of additive clamp.
  -pbe-surface-cutoff <v>
                        auto or explicit Gaussian surface cutoff. Default: auto
  -pbe-gap-diel <v>     Gaussian-cutoff gap dielectric. Default: 22.0
  -pbe-max-grid <N>     Maximum grid points before preflight rejection
  -pbe-map <mode>       off, phi, epsilon, ion, or all. Default: off
  -pbe-map-frames <m>   none, first, last, or all. Default: none
  -pbe-kernel <mode>    auto, scalar, or avx2. Default: auto
                        v1 records the request but uses scalar-threaded SOR
                        Arg file form: -pbe-kernel=auto
                        CLI form:      -pbe-kernel auto
  -pbe-focus <mode>     off, auto, or delphi. Default: off
                        auto runs nested-grid focusing: coarse solve first,
                        then interpolated coarse potential seeds fine-grid
                        boundary/initial values.
                        delphi uses the parent phimap as the child boundary
                        path and exposes hidden probe/focus fixture
                        diagnostics.
  -pbe-multigrid <mode> off, focus, auto, or vcycle. Default: off
                        focus/auto keep the nested focusing path; vcycle runs
                        a geometric V-cycle correction before final fine-grid
                        convergence polish.
  -pbe-focus-levels <N> Number of nested grids, 1..4. Values >1 imply auto.
                        Default effective value: 1 when off, 2 when auto.
  -pbe-focus-margin <A> Extra coarse-grid outer margin in Angstrom. Default: 8.0
  -pbe-mg-levels <N>    auto or 1..4. Default: auto
  -pbe-mg-pre <N>       V-cycle pre-smoothing sweeps. Default: 2
  -pbe-mg-post <N>      V-cycle post-smoothing sweeps. Default: 2
  -pbe-mg-coarse-iter <N>
                        Coarsest-grid SOR iterations per cycle. Default: 200
  Hidden PQR parity diagnostics:
  -pbe-parity-dump <m>  off, summary, points, or all. Requires -pbe-test-pqr.
                        Writes pbe_parity_stage.csv and optionally
                        pbe_parity_points.csv for DelPhi root-cause reports.
                        Barnase PQR rows also emit hidden forensic input,
                        surface, operator, and per-atom phi audits so linear
                        field/source/energy drift can be isolated before
                        nonlinear Newton tuning.
  -pbe-parity-map-dump <m>
                        off, samples, or full. Requires -pbe-test-pqr.
                        samples writes pbe_field_samples.csv,
                        pbe_boundary_samples.csv, and
                        pbe_charge_grid_audit.csv for barnase field-oracle
                        comparison. full is intentionally large. Use
                        test/pbe_delphi_export_maps.ps1 to generate optional
                        DelPhi oracle CSVs outside production paths. The
                        extractor normalizes pyDelPhi VDW surfmap semantics
                        as 1=solvent and 0=solute and writes a sparse
                        pbe_oracle_surface_map.csv containing inside=1 rows
                        for replay. It also writes midpoint-epsilon samples by
                        default, or a full pbe_oracle_midpoint_epsilon.csv for
                        one salt (default 0.15 M, set with -FullMidpointSalt)
                        when called with -FullMidpointEpsilon. If local pyDelPhi
                        cannot export maps, the script remains report-only and
                        strict barnase parity stays disabled.
                        test/pbe_midpoint_oracle_roundtrip_smoke.ps1 verifies
                        the full midpoint oracle layout, axis order, positive
                        face convention, and midpoint six-face sum before
                        native remap changes are trusted.
                        test/pbe_barnase_field_oracle_smoke.ps1
                        writes pbe_barnase_root_cause_report.csv, and
                        test/pbe_barnase_timing_compare.ps1 compares TITAN and
                        pyDelPhi wall-clock timing on the same barnase rows.
  -pbe-test-charge-grid <m>
                        native or oracle. Requires -pbe-test-pqr. oracle reads
                        pbe_oracle_charge_grid_full.csv or
                        pbe_oracle_charge_grid.csv from the run directory.
  -pbe-test-charge-grid-diff <m>
                        off, samples, or full. Requires -pbe-test-pqr and a
                        pyDelPhi charge-grid oracle. Writes
                        pbe_charge_grid_diff.csv and
                        pbe_charge_grid_diff_summary.csv to isolate local
                        charge assignment before SOR update debugging.
  -pbe-test-surface-map <m>
                        native or oracle. Requires -pbe-test-pqr. oracle reads
                        pbe_oracle_surface_map.csv with i,j,k,inside columns.
  -pbe-test-midpoint-epsilon <m>
                        native or oracle. Requires -pbe-test-pqr. oracle reads
                        pbe_oracle_midpoint_epsilon.csv with
                        i,j,k,eps_x,eps_y,eps_z columns and overrides
                        positive-face epsilon before denominator/solver setup.
  -pbe-test-midpoint-diff <m>
                        off, samples, or full. Requires -pbe-test-pqr and
                        pbe_oracle_midpoint_epsilon.csv. Writes only native
                        vs oracle positive-face mismatches to
                        pbe_midpoint_face_diff.csv with remap action, raw face
                        state, signed encoded media/entity state, owner, and
                        nearest atom/exposed-point metrics. full also writes
                        pbe_native_midpoint_state.csv and
                        pbe_vdwms_action_histogram.csv for VDWMS stage
                        root-cause analysis.
  -pbe-test-vdwms-remap <m>
                        auto or off. Requires -pbe-test-pqr. auto enables the
                        DelPhi-parity VDW-to-MS midpoint-media remap; off keeps
                        the raw VDW half-grid media map for forensic runs.
                        pbe_parity_stage.csv and
                        pbe_barnase_surface_audit.csv append
                        midpoint_eps_sum_* columns to distinguish solver
                        denominator audits from raw midpoint-epsilon sums.
  -pbe-test-boundary-map <m>
                        native or oracle. Requires -pbe-test-pqr. oracle reads
                        a complete pbe_oracle_boundary_map.csv boundary phi
                        map and fails clearly if any boundary point is absent.
                        Oracle CSVs are hidden PQR parity inputs only; arg
                        files use -flag=value, while direct CLI usage uses
                        -flag value.
  -pbe-test-boundary-convention <m>
                        native-real, delphi-grid, or oracle. Requires
                        -pbe-test-pqr. delphi-grid samples boundary distances
                        from DelPhi atom-grid coordinates and uses DelPhi
                        exdi/indi scaled by epkt; oracle implies boundary-map
                        replay.
  -pbe-test-field-units <m>
                        native or delphi-dimensionless. Requires
                        -pbe-test-pqr. delphi-dimensionless stores the field
                        in pyDelPhi raw units, scales RHS by epkt, and writes
                        raw/native phi columns in field, boundary, and
                        operator-stage diagnostics.
  -pbe-test-operator-convention <m>
                        native or delphi. Requires -pbe-test-pqr. delphi keeps
                        native field units but uses DelPhi's solvent screening
                        convention in the linear operator audit path.
  -pbe-test-rhs-map <m>, -pbe-test-screen-map <m>, -pbe-test-denom-map <m>
                        native or oracle. Requires -pbe-test-pqr and a full
                        pbe_oracle_operator_stage.csv. These replay RHS,
                        screen, or denominator maps to isolate water_phi
                        mismatch before changing the solver.
  -pbe-test-nl-source <m>
                        titan or delphi-dimensionless. Requires -pbe-test-pqr.
                        delphi-dimensionless solves psi=alpha*phi with alpha-
                        scaled RHS/source, then converts maps/energy back to
                        physical phi for parity diagnostics.
  -pbe-test-newton-derivative [bool]
                        Finite-difference check for the Newton Jacobian.
                        Requires -pbe-test-pqr.
  -pbe-test-energy-convention <m>
                        auto, atom, grid, delphi-corrected, or delphi-iscm.
                        Requires -pbe-test-pqr. Production default remains auto.
                        delphi-iscm is a hidden DelPhi two-dielectric parity
                        convention for induced-surface-charge water RF energy;
                        its final corrected RF follows the DelPhi VDW
                        water-phase ISCM convention while keeping vacuum grid
                        RF as an audit field.
                        TITAN now builds DelPhi-style scaled surface-charge
                        positions natively for delphi-iscm/delphi-corrected.
                        When pbe_oracle_iscm_full.csv is present in a hidden
                        PQR parity run, it is used only for compare/audit
                        fields such as iscm_position_rms_vs_oracle, not for
                        final ISCM accumulation.
                        native/PBSA defaults still use the normal atom/grid
                        convention unless this test path is selected.
                        Use the DelPhi stage-oracle smoke scripts, not this
                        hidden knob alone, to decide whether an energy mismatch
                        is field/operator or reaction-energy convention drift.

  -emap [type]          Energy map type
  -emap-level <level>   per-res or all
  -fmap [type]          Force map type
  -fmap-level <level>   per-res or all
  -energy-res           Per-residue energy decomposition
  -force-res            Per-residue force decomposition
  -new-emap             New energy map format
  -com_energy <axis>    x, y, z, all, or default euclidean
  -interchain-energy    Pairwise inter-chain energy analysis
  -useForceOptionforEnergy
                        Reuse force options for energy output
  -denfo <val>          Distance-energy cutoff
  -dsc <val>            Standard distance cutoff. Default: 6.5 A

  -entropy              Interaction entropy
  -entropy-quasi <com> <lig> <rec>
                        Quasi-harmonic entropy

4.4 PBC

  -pbc                  Basic PBC correction
  -pbci <chain_ref>     Internal-reference PBC
  -pbci2 <chain_ref>    More aggressive internal-reference PBC
  -pbci3 <chain_ref>    Most aggressive internal-reference PBC
  -pbca <chain_ref>     Aggressive PBC mode
  -pbcn                 FixBox-style PBC mode
  -pbctopo              pbctopo (VIBE): Virtual Intact Block Extraction PBC mode
  -pbctopo-anchor-backend=<name>
                        auto | prim | spatial-all. Default: auto
  -pbcclust <axes>      Clustered PBC mode; axes like x, xy, xyz
  -pbcm [chain_ref]     Membrane-oriented PBC mode
  -refpbc <file>        Separate reference PDB for PBC correction
  -pbc-envi <chains...> Environment chains for PBC handling
  -no-center            Disable COM centering during PBC correction

4.5 LJ, RMS, Network, PCA, TICA, MSM, DHAM, LLPS, GIST, IFP, LSP-MD

  -network              DCCM network analysis
  -community=louvain    Louvain community detection on the DCCM graph

  -pca                  PCA and FES analysis
  -fes-grid=<N>         FES grid size. Default: 60
  -pca-solver=<mode>    PCA spectral solver: full | partial (default: full)
  -pca-topk=<N>         Top spectrum size for partial solver; 0 = auto from FES needs
  -pca-max-k-full=<N>   Full solver cutoff K; if K > cutoff, auto-switch full -> partial
  -pca-symmetry-aware   Canonicalize equivalent chain blocks before PCA alignment

  -allostery            MI-derived GC + Transfer Entropy analysis
  -allostery-tau=<N>    TE lag in frames. Default: 5
  -allostery-bins=<N>   Histogram bins. Default: auto (-1)
  -allostery-knn=<N>    Keep top-N ranked partners per residue. Default: 0
  -allostery-knn-prior=<frame|contact>
                        frame = nearest-neighbor votes; contact = Cbeta/N contact frequency
  -allostery-knn-contact-cut=<A>
                        Anchor-contact cutoff for contact prior. Default: 8 A
  -allostery-cut=<A>    Distance-prune candidate residue pairs. Default: off
  -allostery-prune-window=<avg|first>
                        avg scans all stored frames; first only frame 0 for pruning
  -allostery-cut-mode=<avg|first>
                        Deprecated alias for -allostery-prune-window
  -allostery-no-te      Skip Transfer Entropy stage and emit header-only TE CSV
  -allostery-estimator=<mode>
                        hist | knn | kde. Default: hist
  -allostery-estimator-k=<N>
                        Exact kNN estimator neighborhood size. Default: 4
  -allostery-kde-bandwidth=<val>
                        KDE bandwidth override for estimator=kde. Default: 0 (auto)

  -mem-budget-gb=<auto|N>
                        Shared post-processing RAM budget. Default: auto
                        (60% currently available RAM; lower in low-memory mode)
  -low-memory=<true|false>
                        Prefer RAM-saving fallbacks when estimates are large
  -network-mode=<auto|dense|sparse>
                        DCCM network memory policy. Default: auto
  -network-sparse-cut=<A>
                        Sparse candidate cutoff in Angstrom. Default: 30
  -network-sparse-k=<N>
                        Sparse top-N nearest candidates per node. Default: 24
  -network-cache-precision=<auto|fp32|fp64>
                        Cache precision policy metadata for large network runs
  -graph                Unified sparse graph analytics
  -graph-source <name>  dccm | lspmd | allostery-gc | allostery-te | tpt
                        Repeatable. Default: auto (attach to active producers)
  -graph-community <m>  auto | off | louvain. Default: auto
  -graph-filter-topk <N>
                        Keep only top-N weighted edges per node after source
                        construction. Default: 0
  -graph-filter-component <m>
                        all | largest. Default: all

  -corr-knn             Raw continuous kNN MI / GC / TE export
  -corr-knn-k=<N>       Exact kNN estimator neighborhood size. Default: 4
  -corr-knn-tau=<N>     TE lag in frames. Default: 5
  -corr-knn-pair-cut=<A>
                        Distance-prune candidate residue pairs. Default: off
  -corr-knn-pair-topk=<N>
                        Keep top-N ranked partners per residue. Default: 0
  -corr-knn-pair-prior=<frame|contact>
                        frame = nearest-neighbor votes; contact = Cbeta/N contact frequency
  -corr-knn-pair-contact-cut=<A>
                        Anchor-contact cutoff for contact prior. Default: 8 A
  -corr-knn-prune-window=<avg|first>
                         avg scans all stored frames; first only frame 0 for pruning
  -corr-knn-no-te       Skip Transfer Entropy stage and emit header-only TE CSV

  -hbk                  HBond / salt / water kinetics analysis
  -hbk-kind <mode>      all | hbond | salt | water. Repeatable. Default: all
  -hbk-domain <mode>    inter | intra | both. Default: both
  -hbk-maxlag <N>       Max lag in frames for C(t). Default: 500
  -hbk-hbond-cutoff <A>
                        HBond donor-acceptor cutoff. Default: effective -dhb or 3.5
  -hbk-angle <deg>      HBond angle cutoff. Default: effective -ahb or 135
  -hbk-salt-cutoff <A>  Salt-bridge cutoff. Default: 4.0
  -hbk-water-cutoff <A>
                        Water-bridge O...solute cutoff. Default: effective
                        -dis-water or 3.5

  -tica                 TICA analysis
  -tica-lag=<N>         TICA lag. Default: 10
  -tica-components=<N>  Components kept. Default: 5
  -tica-backend=<name>  auto | exact_dense | exact_dual | randomized
  -tica-rank=<N>        Reduced-rank target for large-system TICA. Default: 64

  -zts                  Zero-Temperature String minimum-energy path on PCA/TICA FES
  -zts-domain=<mode>    pca | tica. Default: pca unless both domains are active
  -zts-images=<N>       Number of string images. Default: 50
  -zts-iter=<N>         Maximum ZTS iterations. Default: 2000
  -zts-step=<val>       Gradient-descent step size. Default: 0.01
  -zts-tol=<val>        Convergence tolerance. Default: 1e-5
  -zts-source=<id>      Source basin sink_id. Default: auto lowest basin
  -zts-sink=<id>        Sink basin sink_id. Default: auto second-lowest basin
  -zts-sigma=<val>      Gaussian smoothing radius in grid units. Default: 1.0
  -zts-no-smooth        Disable pre-smoothing of the FES grid
  -zts-init=<mode>      Initialization mode. V1 supports: linear

  -msm                  MSM analysis
  -msm-clusters <N>     K-means cluster count. Default: 25
  -msm-lag <N>          MSM lag. Default: 10
  -msm-nts <N>          Number of implied timescales. Default: 3
  -msm-kmeans-iter <N>  K-means iterations. Default: 300
  -msm-seed <N>         Random seed. Default: 42
  -msm-macrostates <N>  Number of macrostates. Default: -1 (auto)
  -msm-no-boundary-split
                        Disable boundary splitting
  -msm-no-pcca          Disable PCCA
  -msm-no-its           Disable implied-timescale scan
  -msm-its-lags <list>  Comma-separated lag list
  -msm-its-solver <m>   ITS spectral solver: full | schur | partial (default: full)
  -msm-its-topk <N>     Top spectrum size for partial solver; 0 = auto from -msm-nts
  -msm-its-max-k-full <N>
                         Full solver cutoff K; if K > cutoff, auto-switch full -> schur
  -eq-weight-source <mode>
                        auto | uniform | dham | mbar. Default: auto
                        Controls which equilibrium frame weights feed the
                        default PCA/TICA FES path and equilibrium-only
                        MSM/TPT summaries. `auto` preserves legacy behavior.

  -dham                 DHAM reweighting
  -dham-temps <file>    Replica temperatures in K
  -dham-ener <file>     File containing list of GROMACS .xvg energy files
  -dham-target-temp <K> Target reweighting temperature. Default: 300
  -dham-ensemble <mode> NVT or NPT DHAM observable. Default: nvt
  -dham-pressure-bar <bar>
                        Target pressure for NPT enthalpy term. Default: 1.0
  -dham-max-iter=<N>    Maximum DHAM iterations. Default: 5000
  -dham-tol=<val>       DHAM convergence tolerance. Default: 1e-6
  -dham-n-replicas=<N>  Explicit replica count validation
  -dham-dt-energy <ps>  Energy sampling interval in ps; 0 = auto
  -dham-window-match <mode>
                        auto | off | strict. Default: auto
                        Raw DHAM .xvg input also honors global -b/-e/-dt

  -mbar                 MBAR reweighting for REMD/HREMD temperature ladders
-mbar-temps <file>    Optional ladder temperatures in K; required for fixed mode,
                      optional override for per-frame mode
  -mbar-ener <file>     File containing list of GROMACS .xvg energy files
  -mbar-target-temp <K> Target reweighting temperature. Default: 300
  -mbar-domain <mode>   tica | pca | both. Default: tica
  -mbar-ensemble <mode> nvt | npt. Default: nvt
  -mbar-window-match <mode>
                         auto | off | strict. Default: auto
  -mbar-temp-mode <mode>
                        fixed | per_frame. Default: fixed
  -mbar-temp-col <col>  XVG temperature-column index for raw per-frame mode.
                        Default: 2
  -mbar-demux-map <file>
                        Optional per-frame state-assignment map for raw
                        replica outputs
  -mbar-pressure-bar <bar>
                        Target/reporting pressure for NPT MBAR. Default: 1.0
  -mbar-dt-energy <ps>  Energy sampling interval in ps; 0 = auto
  -mbar-bootstrap <N>   Bootstrap replicates for Cv/<E> uncertainty. Default: 0
  -mbar-scan-points <N> Number of temperature scan points; >1 writes mbar_cv.dat
                        Raw MBAR .xvg input also honors global -b/-e/-dt

  -tpt                  Transition Path Theory on the built MSM
  -tpt-mode <mode>      ensemble | kinetic. Default: ensemble
  -tpt-source <list>    Basin spec, e.g. micro:0,1-3 | macro:0 |
                        tica_rect:-1,-1,1,1 | lj:q6_bar>=0.172 |
                        ljpreset:gas | protein:major_macro |
                        protein:rmsd_low | protein:helical
  -tpt-sink <list>      Basin spec, same syntax as -tpt-source
  -tpt-top-paths <N>    Dominant pathways to export. Default: 5
  -tpt-flux-cut <val>   Absolute flux cutoff. Default: 0.0
  -tpt-use-dham <yes/no>
                        Use converged DHAM TPM/pi when available.
                        Default: yes
  -tpt-demux-map <file>
                        Optional demux replica-map for kinetic-TPT. Each row
                        is one time step; each column is one demuxed segment.
  -tpt-bootstrap <N>    Segment bootstrap replicates for uncertainty. Default: 0
  -tpt-bootstrap-seed <N>
                        RNG seed for TPT bootstrap. Default: 42
                        `-eq-weight-source mbar` writes
                        tpt_equilibrium_overlay.dat without changing kinetics

  -llps                 LLPS analysis
  -llps-rcut <val>      Contact cutoff (nm). Default: 0.7
  -llps-perc <val>      Percolation threshold. Default: 0.5
  -llps-fuzzy           Fuzzy EM refinement after DBSCAN
  -llps-dynamics        Enable liquid-like dynamics outputs
  -llps-dynamics-detailed
                        Export detailed chain/pair transition traces
  -llps-maxlag <N>      Max lag for C(t). Default: 500
  -llps-dimer-size <N>  Chains per macro-node. Default: 1
  -llps-contacts        Export contact map data
  -llps-sticker <beads> Sticker bead names
  -llps-voxel <val>     Voxel size for acceleration. Default: 0.5
  -llps-lineage-overlap <val>
                        Jaccard threshold for lineage matching. Default: 0.30
  -llps-prna-rcut <val> Protein-RNA chain contact cutoff (nm). Default: llps-rcut
  -llps-density=<auto|full|projection|off>
                        Density accumulation policy. Default: auto
  -llps-kernel=<auto|scalar|avx2>
                        LLPS hot-kernel policy for benchmark/regression. Default: auto
  -llps-rna-compaction=<auto|true|false>
                        RNA compaction vs protein-binding metrics. Default: auto
  -llps-energy-correlation=<auto|true|false>
                        Contact-gated protein-RNA energy vs RNA compaction
  -llps-energy-detailed=<true|false>
                        Export protein-chain/RNA-chain energy trace. Default: false
  -llps-paper=<true|false>
                        Opt-in Martini protein/RNA paper bundle. Enables LLPS,
                        contacts, dynamics, DBSCAN paper summaries, composition
                        profiles, partition proxies, and convergence report.
                        Default: false
  -llps-paper-profile-bin <nm>
                        Profile bin width for radial/slab paper profiles.
                        Default: 0.2
  -llps-paper-blocks <N|auto>
                        Block count for paper SEM/convergence summaries.
                        Default: auto
  -llps-paper-geometry=<auto|droplet|slab>
                        Profile geometry policy. Default: auto
  -llps-rna-information=<auto|true|false>
                        RNA information-theoretic entropy and coupling metrics
  -llps-rna-topology=<auto|true|false>
                        RNA bridge/clamp/topology metrics in dense phase. Default: auto
  -llps-rna-topology-detailed=<true|false>
                        Export detailed topology/bridge/clamp traces. Default: true
  -llps-rna-intra-rcut <val>
                        RNA intra-chain contact cutoff (nm). Default: llps-rcut
  -llps-rna-longrange-min <N>
                        Minimum residue span counted as long-range RNA contact. Default: 10
  -llps-rna-bridge-geom-min <val|auto>
                        Minimum geometric bridge span (nm). Default: auto
  -llps-rna-patch-gap-max <N>
                        Max residue gap merged into an RNA patch. Default: 1
  -llps-rna-clamp-require-loop <true|false>
                        Require loop support when calling RNA clamps. Default: true
  -llps-rna-info-rg-bins <N>
                        Number of Rg bins for RNA information states. Default: 16
  -llps-rna-info-longrange-bins <N>
                        Number of long-range-contact bins for RNA information states. Default: 16
  -llps-msd-fit-start-ns <val|auto>
                        Dense-phase MSD fit window start. Default: auto
  -llps-msd-fit-end-ns <val|auto>
                        Dense-phase MSD fit window end. Default: auto
  -llps-dbscan-eps <val>
                        DBSCAN epsilon in Angstrom. Default: 15.0
  -llps-dbscan-minpts <N>
                        DBSCAN min neighbours. Default: 5

  -gist                 GIST water thermodynamics
  -gist-resolution=<val>
                        Grid resolution in Angstrom. Default: 0.5
  -gist-padding=<val>   Solute padding in Angstrom. Default: 8.0
  -gist-temperature=<val>
                        Temperature in K. Default: 300.0
  -gist-water-model=<name>
                        Water model. Default: TIP3P
  -gist-min-count=<N>   Min voxel count. Default: 10
  -gist-cutoff=<val>    Interaction cutoff in Angstrom. Default: 9.0

  -ifp                  Interaction fingerprint analysis for protein-protein
                        and protein-RNA/DNA systems
  -ifp-hbond <val>      H-bond cutoff. Default: 3.5
  -ifp-salt <val>       Salt-bridge cutoff. Default: 4.0
  -ifp-hphob <val>      Hydrophobic cutoff. Default: 4.0
  -ifp-pipi <val>       Pi-pi cutoff. Default: 5.5
  -ifp-cpi <val>        Cation-pi cutoff. Default: 6.0
  -ifp-hal <val>        Halogen cutoff. Default: 4.0
  -ifp-out <file>       Custom IFP output filename
                        Extra RNA/DNA annotations still go to
                        ifp_rna_annotations.csv when applicable
  -awh-ifp              RC-binned IFP analysis for GROMACS AWH trajectories
                        Auto-enables -ifp and keeps classic IFP output intact
                        Writes awh_frame_assignments.dat and
                        awh_exploratory_context.dat for exploratory-only
                        biased-ensemble overlays when downstream context exists
  -awh-xvg <file>       Text file listing one walker RC .xvg path per line
  -awh-xvg-col <N>      RC value column index in the .xvg file. Default: 1
  -awh-bins <N>         Number of equal-width RC bins. Default: 20
  -awh-out-prefix <p>   Prefix for AWH representative PDB filenames

  -dssp                 DSSP-like secondary-structure analysis
                        Default: ON
                        Protein/protein-user only
                        In -nocal mode, only molecule 1 is analyzed
  -no-dssp              Disable DSSP analysis
  -dssp-hbond-cutoff <val>
                        DSSP H-bond energy threshold. Default: -0.5

  -lj                   LJ condensation / cluster-density analysis
                        Uses all atoms of molecule 1
                        If no PBC mode is set, Titan auto-enables -pbci 0
  -lj-cut <val>         Cluster cutoff in Angstrom. Default: 5.0
  -lj-min-cluster <N>   Minimum size counted as a cluster. Default: 2
  -lj-rho-bin <val>     rho(r) bin width in Angstrom. Default: 1.0
  -lj-rho-max <val>     rho(r) max radius in Angstrom. Default: auto
  -lj-solid-q6 <val>    q6 threshold for solid-like atoms. Default: 0.40
  -lj-solid-min-cn <N>  Minimum CN for dense/solid-like atoms. Default: 6
  -lj-crystal-tol <val> Crystal classification tolerance. Default: 0.08
  -sq                   Structure-factor / RDF analysis on molecule 1
                        Channels: all, trace, backbone, sticker
                        If no PBC mode is set, Titan auto-enables -pbci 0
  -sq-channel <name>    Repeated selection channel. Default: all
  -sq-rdf-bin <val>     RDF bin width in Angstrom. Default: 1.0
  -sq-rdf-max <val>     RDF max radius in Angstrom. Default: auto
  -sq-q-bin <val>       S(q) bin width in A^-1. Default: 0.05
  -sq-q-max <val>       S(q) max q in A^-1. Default: auto

  -rms                  RMSD + RMSF analysis
                        Atomistic protein/protein-user uses CA
                        Atomistic RNA/DNA uses C4' (fallback C4*)
                        MARTINI uses BB (fallback BB2)
                        In -nocal mode, only molecule 1 is analyzed
                        If no PBC mode is set, Titan auto-enables -pbci 0
  -rmsd                 Alias for -rms
  -rmsf                 Alias for -rms
  -rg2                  Radius-of-gyration analysis on the RMS trace set
                        Outputs Rg^2 in A^2 and Rg in A
                        In -nocal mode, only molecule 1 is analyzed
                        If no PBC mode is set, Titan auto-enables -pbci 0
  -rg2-hist-bins <N>    Histogram bin count for rg2_hist.dat. Default: 100

  -lspmd                Local Spatial Pattern MD analysis
  -lspmd_cut=<val>      C-alpha cutoff in Angstrom. Default: 12.0

  -cas                  Computational alanine scanning
  -tda                  Topological data analysis
  -tda-max-points=<N>   Reservoir cap for TDA point-cloud CSV. Default: 200000
  -ml                   ML affinity prediction

4.6 Hydrogen Bonds / Water / Geometry

  -hb                   Hydrogen-bond analysis
  -hb side-chain        Side-chain-only hydrogen bonds
  -ahb <val>            H-bond angle cutoff. Default: 135.0
  -dhb <val>            H-bond distance cutoff. Default: 3.5

  -water-type <type>    Water mode: surface or sphere
  -dis-water <val>      Water distance cutoff. Default: 3.5

  -charged_atom_dis [N] Charged-atom distance distribution
  -charged_distances_min <val>
                        Min charged distance
  -charged_distances_max <val>
                        Max charged distance
  -angle_dipole_efield [N]
                        Dipole-electric-field angle distribution
  -cal-mindis           Minimum distance per frame
  -resdis2com           Residue-to-COM distance analysis
  -com-dis <axis>       x, y, z, or default euclidean

4.7 Miscellaneous

  -nocal                No-energy mode; useful for PBC/utility workflows
                        DSSP in this mode analyzes molecule 1 only
                        LJ in this mode analyzes molecule 1 only
                        RMSD/RMSF in this mode analyzes molecule 1 only
  -cscd                 Side-chain contact distances
  -pqr                  Write PQR-style output
  -cal-force            Calculate atomic forces
  -cal-volume           Molecular volume
  -max-force            Report maximum force
  -cache                Enable cache
  -cache -no            Disable cache
  -use-CI               Use confidence intervals
  -elec-screening       Electrostatic screening
  -dipole-moment        Dipole moment analysis
  -split                Split results per trajectory
  -all-subfolders       Process all subfolders and split outputs
  -cat-result           Concatenate results from multiple inputs
  -all-respair          Output all residue pairs
  -contact-map no       Disable contact map output
  -act                  Use all side-chain contacts
  -new_al               New alignment in single-PDB mode
  -new-al               Alias of -new_al
  -side_chain_contact <mode>
                        full-res or other
  -noprint-side_chain_contact
                        Suppress side-chain contact printing

4.8 Pull Direction Optimizers

  -pull-direct          Pulling-direction search on the first input only.
                        Uses the last snapshot of the first pdb/xtc/trr input.
                        If -e is set, uses the last snapshot with time <= -e.
                        TITAN exits after this module finishes and ignores
                        remaining inputs for pull-direct.
  -pull-opt <mode>      Optimizer mode. Current build supports: de, de+pattern,
                        de-shade, de-l-shade, de-jso, de-jso+pattern,
                        de-nl-shade-lbc, de-nl-shade-lbc+pattern, direct,
                        direct+pattern
  -pull-max-evals <N>   Global evaluation budget. Default: 0 (unlimited)
  -pull-seed <N>        Deterministic seed for pull-direction search. Default: 42
  -de-diff-weight <val> Differential weight for legacy de only. Default: 0.9
  -de-cross-prob <val>  Crossover probability for legacy de only. Default: 0.9
  -de-F <val>           Scaling factor for legacy de only. Default: 0.5
  -de-scan-radius <val> Scan radius. Default: 2.0
  -de-min-scan-radius <val>
                        Minimum scan radius. Default: 0.5
  -de-score-function <N>
                        Score function index
  -de-pop-num <N>       Population size. Default: 6480
  -de-ite-num <N>       Maximum iterations. Default: 100
  -de-algorithm <N>     Legacy de mutation variant
  -de-angle <val>       Restrain angle. Default: 90.0
  -de-auto-term         Auto-termination
  -de-adaptive-F        Adaptive F scaling for legacy de only
  -de-stagnation-limit <N>
                        Stagnation limit. Default: 20
  -de-memory-size <N>   SHADE-family memory size. Default: 10
  -de-p-best-rate <val> SHADE-family p-best rate. Defaults: 0.2 (de-shade,
                        de-l-shade, de-nl-shade-lbc), 0.11 (de-jso,
                        de-jso+pattern)
  -de-archive-scale <v> SHADE-family archive size scale. Default: 1.0
  -de-min-pop-num <N>   Minimum population for de-l-shade / de-jso /
                        de-jso+pattern / de-nl-shade-lbc /
                        de-nl-shade-lbc+pattern. Default: 4
  -de-pathway-group <N> Pathway group on reference side (subset of g1)
  -de-pathway-radius <A>
                        Geometric rescue radius for auto-ordering pathway nodes. Default: 4.0
                        In pathway mode, TITAN reuses DE_out_vector.dat as the
                        previous pull direction, enforces angle(new, old) <=
                        -de-angle, and writes pull_direction_status.txt with
                        STOP/CONTINUE for Titan+Gromacs chaining.
                        Pathway mode supports DE-family plus direct and
                        direct+pattern in the current build.
                        Recovery mode prioritizes reattaching mol2 to the
                        nearest reachable forward point on the pathway before
                        energy. Advance mode prioritizes energy first and then
                        forward progress.
                        Legacy outputs DE_out_vector.dat / DE_log.dat /
                        DE_report.dat / DE_vector.pdb / best_conf.pdb stay in
                        the original input frame. Only
                        pull_direction_aligned_snapshot.pdb and
                        pull_direction_system_info.txt are rotated so that the
                        best pull vector points to +Z.
                        de-nl-shade-lbc is the recommended DE-family branch in
                        the current build. direct+pattern remains the strongest
                        score-focused competitor on the benchmark case.
                        direct can be faster than de on small/medium cases
                        because it often reaches a good answer with far fewer
                        energy evaluations.

4.9 Help / Version

  -h, --help            Print help message
  -v, --version         Print build number

================================================================================
                               END OF MANUAL
================================================================================
)MANUAL";
    return 0;
  }
}
