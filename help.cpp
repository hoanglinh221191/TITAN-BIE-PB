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
    hydrogen bonds, DCCM network analysis, PCA/FES, allostery, TICA, MSM,
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

  - If DHAM converges, PCA/FES automatically uses DHAM frame weights at the
    requested target temperature. This enables weighted FES without changing
    the rest of the PCA workflow.

  - MBAR v1 is standalone-only. It writes `mbar_*` outputs, reweights to a
    target temperature from REMD/HREMD energy ladders, and does not overwrite
    the standard PCA/TICA/DHAM outputs.

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
    average_emap.dat
    average_emap_matrix.dat
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
    gbsa_average_emap.dat
    gbsa_total_energy_maptype.dat
    out_entropy.dat   (auto-enabled by -gbsa for free-energy completeness)

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

2.6a Mutual Information / Generalized Correlation

  Interpretation:
    - `raw_mi_edges.csv` reports the raw pairwise MI estimated from the
      discretized allostery projections.
    - `network_gc_edges.csv` reports the MI-derived GC transform, which maps
      nonlinear dependence onto the bounded interval [0, 1).
    - `network_te_edges.csv` remains the directed lagged-information export.

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
    msm_its_tpm/

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
    If DHAM converges, PCA/FES uses DHAM frame weights automatically.
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
    and does not feed weights back into MSM/TPT/PCA/TICA defaults automatically.
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
    their original behavior.

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
    In kinetic mode, TITAN always uses the built MSM model and writes
    `tpt_kinetics.dat` diagnostics; it still does not certify kinetic
    validity automatically. If `-tpt-demux-map <file>` is supplied, TITAN
    rebuilds a kinetic MSM from the demuxed microstate segments before TPT.
    Kinetic diagnostics now include a first-pass CK / Markov sanity check.
    `-tpt-bootstrap <N>` enables segment bootstrap uncertainty for MSM-backed
    TPT models. DHAM-reweighted bootstrap is not implemented yet.

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

  Main inputs:
    -llps
    optional:
      -llps-rcut
      -llps-perc
      -llps-fuzzy
      -llps-maxlag
      -llps-dimer-size
      -llps-contacts
      -llps-sticker
      -llps-voxel
      -llps-dbscan-eps
      -llps-dbscan-minpts

  Typical outputs:
    llps_timeseries.dat
    llps_autocorr.dat
    llps_degree_dist.dat
    llps_contact_map.dat
    llps_summary.dat
    llps_density_protein.dx
    llps_density_rna.dx
    polymer_scaling.dat
    llps_dbscan_droplets.pdb
    llps_dbscan_timeseries.csv

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

  Inspect:
    raw_mi_edges.csv for unthresholded mutual information
    network_gc_edges.csv for the derived GC network
    network_te_edges.csv for directed transfer entropy

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

  -pattern-input <pat> [b] [e] [step]
                        Pattern input mode; <pat> must contain '*'

4.2 Groups / Molecules / Time

  -g <G1>&<G2>          Interaction groups
  -sg <SG1>-<SG2>       Sub-groups within main groups
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
  -allostery-kde-bandwidth=<val>
                        KDE bandwidth override for estimator=kde. Default: 0 (auto)

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

  -llps                 LLPS analysis
  -llps-rcut <val>      Contact cutoff (nm). Default: 0.7
  -llps-perc <val>      Percolation threshold. Default: 0.5
  -llps-fuzzy           Fuzzy EM refinement after DBSCAN
  -llps-maxlag <N>      Max lag for C(t). Default: 500
  -llps-dimer-size <N>  Chains per macro-node. Default: 1
  -llps-contacts        Export contact map data
  -llps-sticker <beads> Sticker bead names
  -llps-voxel <val>     Voxel size for acceleration. Default: 0.5
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
