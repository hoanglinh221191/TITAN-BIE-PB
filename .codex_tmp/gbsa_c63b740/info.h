#ifndef INFO_H
#define INFO_H
#include "arguments.h"
#include "groups.h"
#include <string>

using std::string;
class info {
  string Build, PdbFileName, TotalEnergyFileName, TotalEnergyTimeFileName,
      SideChainName, LogFile, EmapLgFileName, EmapRecFileName, EntropyFile,
      EntropyQuasiFile, group_1, group_2, stat_hbonds, stat_average_HB,
      stat_maxforce_res, emap, dis_energy_2mol, res_list_hb_rec, res_list_hb_lg,
      side_chain_out, side_chain_time, out_energy_res_motif, dis_com_2mol,
      dipole_filename;
  double current_begin_time, current_end_time;
  double current_dt;

public:
  info(const arguments &Arguments)
      : Build{Arguments.get_build_number()}, PdbFileName{Arguments.pdb_name},
        TotalEnergyFileName{"cal.txt"}, TotalEnergyTimeFileName{"cal.xvg"},
        SideChainName{"out_sc.txt"}, LogFile{"log.txt"},
        EmapLgFileName{"ener_lg.txt"}, EmapRecFileName{"ener_rec.txt"},
        EntropyFile{"out_entropy.dat"},
        EntropyQuasiFile{"out_entropy_quasiHarmonic.dat"}, group_1{"mol1"},
        group_2{"mol2"}, stat_hbonds{"stat_hydrogenbonds.dat"},
        stat_average_HB{"stat_average_HB.dat"},
        stat_maxforce_res{"_out_stat_maxforce_res.dat"}, emap{"out_emap_"},
        dis_energy_2mol{"out_dis_com_energy.dat"}, res_list_hb_rec{},
        res_list_hb_lg{}, side_chain_out{"out_sc.txt"},
        side_chain_time{"out_sc_time.dat"},
        out_energy_res_motif{"_stat_energy_res.dat"},
        dis_com_2mol{"out_dis_com.dat"}, current_begin_time{1000000000000.},
        current_end_time{}, current_dt{} {}
  string getBuild() const { return Build; }
  int get_info_group(const groups_cal &groups) {
    if (groups.get_parse_index_flag()) {
      group_1 = groups.group_name(groups.group_1_index());
      group_2 = groups.group_name(groups.group_2_index());
    } else {
    }
    return 1;
  }

  string get_pca_fel_filename() const { return "pca_fel.csv"; }
  string get_fes_2d_filename() const { return "fes_2d.dat"; }
  string get_cluster_representatives_filename() const {
    return "cluster_representatives.dat";
  }

  void setPdbName(const arguments &Arguments) {
    PdbFileName = Arguments.pdb_name;
  }
  string get_group_1_name() const { return group_1; }
  string get_group_2_name() const { return group_2; }
  string get_stathb_name() const { return stat_hbonds; }
  string get_stat_ave_hb_name() const { return stat_average_HB; }
  string get_stat_maxforce_name() const { return stat_maxforce_res; }
  string get_emap() const { return emap; }
  string getPdbFileName() const { return PdbFileName; }
  string get_pdboutname(const arguments &Arguments) const {
    return Arguments.get_pdbOutput();
  }
  string getEntropyQuasi() const { return EntropyQuasiFile; }
  string get_res_list_hb_rec_filename() {
    return "out_list_res_HB_" + group_1 + ".dat";
  }
  string get_res_list_hb_lg_filename() {
    return "out_list_res_HB_" + group_2 + ".dat";
  }
  string get_res_list_sc_contact_rec_filename() {
    return "out_list_res_sidechain_contact_" + group_1 + ".dat";
  }
  string get_res_list_sc_contact_lg_filename() {
    return "out_list_res_sidechain_contact_" + group_2 + ".dat";
  }
  string getDisCom2Mol() const { return dis_com_2mol; }
  string getDisEner2Mol() const { return dis_energy_2mol; }
  string getDipoleFilename() const { return "_dipole_moment.dat"; }
  string getDistributionDistanceCharged() const {
    return "_distributionChargedAtoms.dat";
  }
  string getDistributionAngleDipoleEfield() const {
    return "_distributionDipoleEfield.dat";
  }
  string getOutputDistanceRes2COM() const { return "_distances_res_2COM.dat"; }
  string getOuputAverageContactMap() const { return "average_contact_map.dat"; }
  string getOuputEnerMaptypeAve() const { return "_energy_maptype.dat"; }
  string getOutputAverageEmap() const { return "average_emap.dat"; }
  string getOutputAverageEmapMatrix() const {
    return "average_emap_matrix.dat";
  }
  string getOutputEnerPdbtypeAve() const { return "_energy_pdbtype.pdb"; }
  string get_gbsa_snapshots_filename() const { return "gbsa_snapshots.dat"; }
  string get_gbsa_timeseries_avg_filename() const {
    return "gbsa_timeseries_avg.dat";
  }
  string get_gbsa_summary_filename() const { return "gbsa_summary.dat"; }
  string get_gbsa_decomposition_filename() const {
    return "gbsa_decomposition.dat";
  }

  string get_DE_output_filename() { return "DE_out_vector.dat"; }
  string get_DE_log_filename() { return "DE_log.dat"; }
  string get_DE_pdb_vector_filename() { return "DE_vector.pdb"; }

  string get_cas_receptor_filename() const { return "cas_receptor.dat"; }
  string get_cas_ligand_filename() const { return "cas_ligand.dat"; }
  string get_cas_hotspot_pdb_filename() const {
    return "cas_hotspots_b_factor.pdb";
  }

  string get_global_water_bridges_timeseries() const {
    return "global_water_bridges_timeseries.csv";
  }
  string get_global_water_bridges_summary() const {
    return "global_water_bridges_summary.txt";
  }
  string get_global_water_survival() const {
    return "global_water_survival.csv";
  }
  string get_water_pdb_filename(const arguments &Arguments,
                                const std::vector<string> &input_line,
                                int input_name) const {
    if (Arguments.flag_sep_opdb == true && !input_line.empty()) {
      size_t posfile = input_name - Arguments.b_input;
      return Arguments.get_pdbOutput().substr(
                 0, Arguments.get_pdbOutput().length() - 4) +
             "_waters_" + std::to_string(posfile) + ".pdb";
    }
    return "waters_out" +
           (input_line.empty()
                ? ""
                : "_" + input_line[input_name - Arguments.b_input]) +
           ".pdb";
  }

  string get_pca_filename(const arguments &Arguments,
                          const std::vector<string> &input_line,
                          int input_name) const {
    if (Arguments.flag_sep_opdb == true && !input_line.empty()) {
      size_t posfile = input_name - Arguments.b_input;
      return "pca_fel_" + std::to_string(posfile) + ".csv";
    }
    return "pca_fel" +
           (input_line.empty()
                ? ""
                : "_" + input_line[input_name - Arguments.b_input]) +
           ".csv";
  }

  string getSCout() const { return side_chain_out; }
  string getContmap() const { return "contact_map_traj_"; }
  string getSCtimeout() const { return side_chain_time; }
  void set_current_begin_time(double t) { current_begin_time = t; }
  double get_current_begin_time() const { return current_begin_time; }
  void set_current_end_time(double t) { current_end_time = t; }
  double get_current_end_time() const { return current_end_time; }
  void set_current_dt(double dt) { current_dt = dt; }
  double get_current_dt() const { return current_dt; }
  void set_out_energy_res_motif(string motif) { out_energy_res_motif = motif; }
  string get_out_energy_res_motif() { return out_energy_res_motif; }
  string get_ml_features_filename() const { return "titan_ml_features.csv"; }
  string get_ml_prediction_filename() const {
    return "titan_rf_prediction.txt";
  }
  string get_tica_eigenvalues_filename() const {
    return "tica_eigenvalues.dat";
  }
  string get_tica_eigenvectors_filename() const {
    return "tica_eigenvectors.dat";
  }
  string get_tica_projections_filename() const {
    return "tica_projections.csv";
  }
  // --- MSM Phase 2 output files ---
  string get_msm_tpm_filename() const { return "msm_tpm.dat"; }
  string get_msm_discrete_traj_filename() const {
    return "msm_discrete_traj.dat";
  }
  string get_msm_timescales_filename() const { return "msm_timescales.dat"; }
  string get_msm_centroids_filename() const { return "msm_centroids.dat"; }
  string get_msm_its_filename() const { return "msm_its_scan.dat"; }
  string get_msm_macrostates_filename() const { return "msm_macrostates.dat"; }
  // Python-friendly extended outputs
  string get_msm_macrostate_traj_filename() const {
    return "msm_macrostate_traj.dat";
  }
  string get_msm_count_matrix_filename() const {
    return "msm_count_matrix.dat";
  }
  string get_msm_stationary_filename() const { return "msm_stationary.dat"; }
  string get_msm_its_tpm_dir() const { return "msm_its_tpm"; }
  string get_out_mindis_subgroups() { return "mininum_distance_"; }
  string get_out_histogram_mindis() { return "histo_mindistance_"; }
  // LLPS output files
  string get_llps_timeseries_filename() const { return "llps_timeseries.dat"; }
  string get_llps_ct_filename() const { return "llps_autocorr.dat"; }
  string get_llps_degree_dist_filename() const {
    return "llps_degree_dist.dat";
  }
  string get_llps_contact_map_filename() const {
    return "llps_contact_map.dat";
  }
  string get_llps_summary_filename() const { return "llps_summary.dat"; }
  string get_llps_density_pro_filename() const {
    return "llps_density_protein.dx";
  }
  string get_llps_density_rna_filename() const { return "llps_density_rna.dx"; }
  string get_llps_polymer_scaling_filename() const {
    return "polymer_scaling.dat";
  }
  // GIST output prefix
  string get_gist_prefix() const { return "gist"; }
  // LLPS DBSCAN output files
  string get_llps_dbscan_droplets_filename() const { return "llps_dbscan_droplets.pdb"; }
  string get_llps_dbscan_timeseries_filename() const { return "llps_dbscan_timeseries.csv"; }
  // IFP output file
  string get_ifp_filename() const { return "ifp_occupancy.csv"; }
  // DHAM output files
  string get_dham_pi_filename() const { return "dham_stationary_distribution.dat"; }
  string get_dham_fes_filename() const { return "dham_free_energy_profile.dat"; }
  string get_dham_tpm_filename() const { return "dham_unbiased_tpm.dat"; }
};
#endif
