#ifndef ARGUMENTS_H
#define ARGUMENTS_H
#include <string>
#include <vector>
using namespace std;

// ── Domain sub-structs ──────────────────────────────────────
struct PbcArgs {
  int solv_pbc{0};
  bool solv_pbcn{false};
  bool solv_pbctopo{false};
  bool auto_enabled_pbcn{false};
  bool auto_enabled_pbci{false};
  int solv_pbc_int{0}, solv_pbc_int2{0}, solv_pbc_int3{0};
  bool pbc_aggressive{false};
  bool pbc_mem{false};
  bool no_center{false};
  int chain_ref{-1};
  bool solv_pbc_clust{false};
  vector<char> pbc_clust_axis;
  bool pbc_envi{false};
  vector<string> pbc_envi_chain;
  int pbc_mol2_aggressive{0};
  bool flag_pbc_ref{false};
  string ref_pdb_name{};
};

struct MsmArgs {
  bool flag_msm{false};
  bool flag_tica{false};
  // Default/cap for MSM microstates when the user does not pass
  // -msm-clusters explicitly. The engine may auto-reduce this value after
  // TICA finishes if the trajectory is too short for a stable MSM/DHAM fit.
  int msm_clusters{25};
  bool msm_clusters_user_set{false};
  int msm_lag{10};
  int msm_nts{3};
  int msm_kmeans_iter{300};
  int msm_seed{42};
  bool msm_auto_boundary{true};
  int msm_macrostates{-1};
  bool msm_run_pcca{true};
  bool msm_run_its{true};
  vector<int> msm_its_lags{};
  int tica_lag{10}, tica_components{5};
  std::string tica_backend{"auto"};
  int tica_rank{64};
};

struct GistArgs {
  bool flag_gist{false};
  double gist_resolution{0.5};
  double gist_padding{8.0};
  double gist_temperature{300.0};
  std::string gist_water_model{"TIP3P"};
  int gist_min_count{10};
  double gist_cutoff{9.0};
};

struct DeArgs {
  bool flag_find_pulling_direct{false};
  double de_diff_weight{0.9};
  double de_cross_prob{0.9};
  double de_F{0.5};
  double de_scan_radius{2.};
  double de_min_scan_radius{0.5};
  double de_restrain_angle{65.};
  int de_flag_score_function{1};
  size_t de_pop_num{6480};
  size_t de_algorithm{};
  size_t de_ite{100};
  bool de_adaptive_F{true};
  int de_stagnation_limit{20};
  bool flag_de_auto_term{false};
  double dis_along_path{4.};
  bool flag_pathway{false};
  size_t group_pathway{};
};

struct LlpsArgs {
  bool flag_llps{false};
  float llps_rcut{0.7f};
  float llps_perc_thr{0.5f};
  int llps_maxlag{500};
  int llps_dimer_size{1};
  bool llps_contacts{false};
  std::string llps_sticker_beads{"BB,SC1,SC2,BB1,BB2,BB3,SC3,SC4"};
  float llps_voxel_size{0.5f};
  // DBSCAN bead-level morphology parameters (unit: Angstrom)
  double llps_dbscan_eps{15.0};    // epsilon (Å); default 15 Å ≈ 1.5 nm
  int    llps_dbscan_minpts{5};    // min neighbours (self excluded)
  bool   llps_fuzzy{false};        // enable Fuzzy EM refinement after DBSCAN
};

struct IfpArgs {
  bool   flag_ifp             {false};
  double cutoff_hbond         {3.5};   // D···A distance (Å)
  double cutoff_salt_bridge   {4.0};   // charged atom pair (Å)
  double cutoff_hydrophobic   {4.0};   // C···C contact (Å)
  double cutoff_pi_pi         {5.5};   // ring COM–COM (Å)
  double cutoff_cation_pi     {6.0};   // cation–ring COM (Å)
  double cutoff_halogen       {4.0};   // X···acceptor (Å)
  std::string output_filename {};  // empty → uses info.get_ifp_filename()
};

struct DhamArgs {
  bool   flag_dham       {false};
  std::string temps_file {};       // file with replica temperatures (K)
  std::string ener_list_file {};   // file listing paths to .xvg energy files
  std::string ensemble   {"nvt"};  // nvt | npt
  std::string window_match {"auto"}; // auto | off | strict
  double target_temp     {300.0};  // target reweighting temperature (K)
  double target_pressure_bar {1.0}; // target pressure for NPT reweighting (bar)
  int    max_iter        {5000};
  double tolerance       {1e-6};
  int    n_replicas      {0};      // explicit validation (0 = auto from temps)
  double dt_energy_ps    {0.0};    // energy .xvg output freq (ps); 0 = auto from file
};

struct TptArgs {
  bool flag_tpt{false};
  std::string mode{"ensemble"};   // ensemble | kinetic
  std::string source_spec{};      // microstate list/range, e.g. micro:0,1-3
  std::string sink_spec{};        // microstate list/range, e.g. 8-12
  int top_paths{5};
  double flux_cut{0.0};
  bool use_dham{true};
  std::string demux_map_file{};   // optional demux replica-map for kinetic mode
  int bootstrap_samples{0};       // 0 = off
  unsigned int bootstrap_seed{42};
};

struct UmapArgs {
  bool flag_umap{false};
  int umap_neighbors{15};
  int umap_components{2};
  int umap_epochs{200};
  double umap_min_dist{0.1};
  double umap_spread{1.0};
  unsigned int umap_seed{42};
  string umap_knn_mode{"exact"};
  bool umap_knn_benchmark{false};
  bool umap_verbose{false};
};

struct DsspArgs {
  bool flag_dssp{false};
  double hbond_energy_cutoff{-0.5};
};

struct RmsArgs {
  bool flag_rms{false};
};

struct LjArgs {
  bool flag_lj{false};
  double lj_cluster_cutoff{5.0};
  int lj_min_cluster_size{2};
  double lj_rho_bin{1.0};
  double lj_rho_max{0.0};
  double lj_solid_q6{0.40};
  int lj_solid_min_cn{6};
  double lj_crystal_tol{0.08};
};

// ── Main arguments struct ───────────────────────────────────
struct arguments {
private:
  string Buildnumber;
  string set_Buildnumber();
  bool parse_input_file_is_open();
  std::vector<std::string> tokenize_arg_file(const std::string& filename);

public:
  // Domain sub-structs
  PbcArgs pbc;
  MsmArgs msm;
  GistArgs gist;
  DeArgs de;
  LlpsArgs llps;
  IfpArgs ifp;
  DhamArgs dham;
  TptArgs tpt;
  UmapArgs umap;
  DsspArgs dssp;
  RmsArgs rms;
  LjArgs lj;

  // LSP-MD (Kornev, JCTC 2025) — Local Spatial Pattern network analysis
  struct LSPMDArgs {
    bool   flag_lspmd    = false;  // -lspmd: enable module
    double lspmd_cutoff  = 12.0;  // -lspmd_cut: Cα–Cα cutoff (Å)
    bool   lspmd_use_avg = true;  // use average distance for cutoff check
  };
  LSPMDArgs lspmd;

  // Core flags
  bool flag_gbsa, flag_alpb;
  bool cal_dis_sc_contact, opdb, file_par, new_al_1pdb, entropy, entropy_quasi,
      flag_index, flag_topo_1, flag_topo_2, flag_1_pdb, flag_listfile,
      flag_use_common_topo, flag_nocal, flag_pattern_input, flag_range_input,
      help, flag_network, flag_network_explicit, flag_pca, flag_tda, flag_ml,
      flag_cas,
      flag_allostery;
  bool flag_all_folders, flag_separation_results, flag_sidechain_hbond,
      flag_cal_force, flag_use_CI, print_force_res, flag_begin_set,
      flag_end_set, flag_energy_same_force, flag_maxforce, flag_new_emap,
      print_energy_res, self_int, flag_sidechain, flag_print_sidechain,
      flag_all_respair;
  bool flag_sub_groups, flag_cal_contact,
      flag_contact_map_explicit,
      flag_cal_dis_res_2com, flag_set_type_1, flag_set_type_2,
      flag_cal_dipole_moment, flag_charged_atom_distance_dis,
      flag_angle_dipole_efield, flag_sep_opdb, heat_capacity_cp, flag_cal_vol,
      flag_elec_screening, flag_interchain_energy,
      flag_cal_mindis, mol_1_type_explicit, mol_2_type_explicit,
      flag_cat_result_multi_infile, opdb_name_explicit;
  int igb;
  int flag_pqr, flag_all_side_chain_ct, begin_mol_2,
      energy_map, force_map, argc, group_1, group_2,
      flag_h_bond, group_entropy_com, group_entropy_rec,
      group_entropy_lg, k_neighbor, dt, b_input, e_input, step_input,
      flag_cache, com_dis_energy_2mol, com_dis_2mol,
      n_wat;
  int fes_grid;
  int allostery_tau;
  int allostery_bins;
  int allostery_knn;
  double allostery_cutoff;
  bool allostery_use_avg_cutoff;
  bool allostery_compute_te;
  // XDR trajectory (XTC/TRR)
  string traj_file;
  bool flag_xdr;
  string xdr_type;
  vector<int> sub_group_1, sub_group_2;
  size_t number_of_files, nt, number_charged_distances_bin,
      number_bin_angle_dipole_efield, number_mindis_bin;
  double ext_diel, int_diel, saltcon;
  double d_standard, begin_time, end_time, angle_h, dis_h, limit_dis, eu_par,
      fu_par, dis_clust, d_crit_water, distance_energy_cutoff{-1.1},
      d_crit_minmax_mol{4.}, charged_distances_max, charged_distances_min,
      salt_concentration, epsilon_in;
  string topo_name, topo_name_2, pdb_name, out_pdb_name, index_name, force_f,
      command, path_amino_acid, path_nonbonded, name_out_cal_dis, arg_filename,
      list_filename, emap_type, fmap_type, input_pattern, emap_level,
      fmap_level, out_motif, eu, fu, ref_file, tu, input_file_type, mol_1_type,
      mol_2_type, water_type;

  arguments()
      : Buildnumber{set_Buildnumber()},
        // Core bools
        flag_gbsa{false}, flag_alpb{false},
        cal_dis_sc_contact{false}, opdb{false}, file_par{false},
        new_al_1pdb{false}, entropy{false}, entropy_quasi{false},
        flag_index{false}, flag_topo_1{false}, flag_topo_2{false},
        flag_1_pdb{false}, flag_listfile{false}, flag_use_common_topo{false},
        flag_nocal{false}, flag_pattern_input{false}, flag_range_input{false},
        help{false}, flag_network{false}, flag_network_explicit{false},
        flag_pca{false},
        flag_tda{false}, flag_ml{false}, flag_cas{false},
        flag_allostery{false},
        flag_all_folders{false}, flag_separation_results{false},
        flag_sidechain_hbond{false}, flag_cal_force{false}, flag_use_CI{false},
        print_force_res{false}, flag_begin_set{false}, flag_end_set{false},
        flag_energy_same_force{false}, flag_maxforce{false},
        flag_new_emap{false}, print_energy_res{false}, self_int{false},
        flag_sidechain{true}, flag_print_sidechain{true},
        flag_all_respair{false},
        flag_sub_groups{false}, flag_cal_contact{true},
        flag_contact_map_explicit{false},
        flag_cal_dis_res_2com{false}, flag_set_type_1{false},
        flag_set_type_2{false}, flag_cal_dipole_moment{false},
        flag_charged_atom_distance_dis{false}, flag_angle_dipole_efield{false},
        flag_sep_opdb{false}, heat_capacity_cp{false}, flag_cal_vol{false},
        flag_elec_screening{false}, flag_interchain_energy{false},
        flag_cal_mindis{false}, mol_1_type_explicit{false},
        mol_2_type_explicit{false},
        flag_cat_result_multi_infile{false}, opdb_name_explicit{false},
        // Core ints
        igb{5},
        flag_pqr{2}, flag_all_side_chain_ct{1}, begin_mol_2{},
        energy_map{2}, force_map{2}, argc{1}, group_1{}, group_2{},
        flag_h_bond{0},
        group_entropy_com{-1}, group_entropy_rec{-1}, group_entropy_lg{-1},
        k_neighbor{4}, dt{0}, b_input{}, e_input{100000000},
        step_input{1}, flag_cache{1}, com_dis_energy_2mol{4}, com_dis_2mol{0},
        n_wat{-1},
        fes_grid{60},
        allostery_tau{5}, allostery_bins{-1},
        allostery_knn{0},
        allostery_cutoff{0.0}, allostery_use_avg_cutoff{true},
        allostery_compute_te{true},
        // XDR
        traj_file{}, flag_xdr{false}, xdr_type{},
        // size_t
        number_of_files{1ull}, nt{}, number_charged_distances_bin{400},
        number_bin_angle_dipole_efield{200}, number_mindis_bin{100},
        // doubles
        ext_diel{78.5}, int_diel{1.0}, saltcon{0.15},
        d_standard{6.5}, begin_time{0.}, end_time{100000000000.}, angle_h{135.},
        dis_h{3.5}, limit_dis{10e10}, eu_par{1.0}, fu_par{1.66}, dis_clust{6.5},
        d_crit_water{3.5}, distance_energy_cutoff{-1.1}, d_crit_minmax_mol{4.},
        charged_distances_max{50.}, charged_distances_min{0.},
        salt_concentration{0.15}, epsilon_in{1.00000000},
        // strings
        topo_name{"lig.itp"}, topo_name_2{"no"}, pdb_name{"complex.pdb"},
        out_pdb_name{"out_complex.pdb"}, index_name{"index.ndx"},
        force_f{"AMBER99SB"}, command{"-"}, path_amino_acid{"lib_ff/"},
        path_nonbonded{"lib_ff/"}, name_out_cal_dis{"out_cal_dis_com_"},
        arg_filename{"arg.in"}, list_filename{"input.txt"},
        emap_type{"residue-atom"}, fmap_type{"residue-atom"}, input_pattern{},
        emap_level{"per-res"}, fmap_level{"per-res"}, out_motif{"*"},
        eu{"kJ/mol"}, fu{"pN"}, ref_file{}, tu{"ps"}, input_file_type{".pdb"},
        mol_1_type{"protein"}, mol_2_type{"lig"}, water_type{"surface"} {}

  int post_process();
  string get_build_number() const;
  int parse_arguments(int argc, char **argv);
  int parse_arguments_interactive(int argc);
  int parse_arg_file();
  bool check_topo(string NameofTopoFile, int Topo);
  int help_print();
  void set_number_of_files(size_t num) { number_of_files = num; }
  void setPdbname(string name) { pdb_name = name; }
  string get_pdbOutput() const;
  int recheck_and_change_pdb_input();
  string get_input_file_type() const { return input_file_type; }
  int add_input_file_type() {
    if (pdb_name.length() <= 4 ||
        pdb_name.substr(pdb_name.length() - 4) != get_input_file_type()) {
      pdb_name += get_input_file_type();
      return 0;
    } else
      return 1;
  }
};
#endif
