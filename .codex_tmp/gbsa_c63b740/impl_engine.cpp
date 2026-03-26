#include "MOL_para.h"
#include "cas_analysis.h"
#include "chunk_streamer.h"
#include "dccm_network.h"
#include "differential_evolution.h"
#include "functions.h"
#include "gbsa_intrin.h"
#include "groups.h"
#include "memory_utils.h"
#include "miscellaneous.h" // Added as per instruction's Code Edit
#include "ml_extractor.h"
#include "nucleic_profiler.h"
#include "tda.h"
#include "titan_ui.h"

#include "engine_postproc.h"
#include "titan_alignment.h"
#include "titan_gist.h"
#include "titan_kinetics.h"
#include "titan_llps.h"
#include "titan_llps_cluster.h"
#include "titan_msm.h"
#include "titan_ifp.h"
#include "titan_pca.h" // Ensured inclusion as per instruction
#include "water_analysis.h"


#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <fcntl.h>
#include <io.h>
#include <windows.h>

#endif

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#pragma warning(disable : 4996)
using namespace std;
using namespace chrono;
int engine(arguments &Arguments) {
#ifdef _WIN32
  _setmode(_fileno(stdout), _O_TEXT);
  _setmode(_fileno(stderr), _O_TEXT);
  setvbuf(stdout, NULL, _IONBF, 0);
  setvbuf(stderr, NULL, _IONBF, 0);
  SetConsoleOutputCP(CP_UTF8);
  SetConsoleCP(CP_UTF8);

  // KHÔNG bật ENABLE_VIRTUAL_TERMINAL_PROCESSING ở đây nữa
  // Windows Terminal tự bật VT cho mọi app hiện đại
  // Việc bật thủ công kích hoạt DA response từ ConPTY
#endif
  TitanUI::print_logo();
  Arguments.post_process();
  /************************************/
  folders paths;
  vector<string> list_paths;
  de_arg DE_arguments;
  if (Arguments.de.flag_find_pulling_direct == true)
    DE_arguments.parse_input(Arguments);
  paths.add_paths(Arguments, list_paths);
  paths.create_folders(Arguments);

  int sc_flag{Arguments.flag_all_side_chain_ct};
  auto total_run_time{chrono::steady_clock::now()};
  chrono::milliseconds time_energymaptype{};
  string force_f{Arguments.force_f};

  if (Arguments.energy_map > 0) {
    fstream lg_ener_inter("ener_lg.txt", ios::out);
    if (Arguments.mol_2_type.find("protein") != string::npos)
      lg_ener_inter << "#Res" << setw(11) << "Electro." << setw(12)
                    << "Van der W." << '\n';
    else
      lg_ener_inter << "#Atom/Res" << setw(11) << "Electro." << setw(12)
                    << "Van der W." << '\n';
    fstream rec_ener_inter("ener_rec.txt", ios::out);
    rec_ener_inter << "#Res" << setw(12) << "Electro." << setw(12)
                   << "Van der W." << '\n';
    lg_ener_inter.close();
    rec_ener_inter.close();
  }

  if (Arguments.flag_cal_force && Arguments.force_map > 0) {
    fstream lg_force_inter("force_lg.txt", ios::out);
    if (Arguments.mol_2_type.find("protein") != string::npos)
      lg_force_inter << "#Res" << setw(11) << "Electro." << setw(12)
                     << "Van der W." << '\n';
    else
      lg_force_inter << "#Atom/Res" << setw(11) << "Electro." << setw(12)
                     << "Van der W." << '\n';
    fstream rec_force_inter("force_rec.txt", ios::out);
    rec_force_inter << "#Res" << setw(12) << "Electro." << setw(12)
                    << "Van der W." << '\n';
    lg_force_inter.close();
    rec_force_inter.close();
  }

  /*************************************/
  double d_standard{Arguments.d_standard};
  int num_amine{};
  size_t num_threads{Arguments.nt}, tot_snap_allfiles{};
  string name_pro{}, name_lig{};
  para parameter;

  vector<double> E_v, E_e;
  vector<size_t> size_snap;
  vector<vector<double>> E_v_vector, E_e_vector, Volumes;
  vector<force> Force_e, Force_v;
  vector<vector<force>> Force_e_vector, Force_v_vector;
  vector<energy_res> Force_res_stat_mol1, Force_res_stat_mol2;
  vector<vector<energy_res>> Energy_res_stat_1, Energy_res_stat_2;
  vector<vector<energy_groups>> EF_subgroups;
  vector<vector<map_pos>> average_contact_matrix;
  vector<mindistance_group> MIN_DIS;
  vector<vector<size_t>> Hbonds_vector, sidechain_contact_num;
  amino_acid amine;
  vector<double> time_of_snap;
  vector<string> input_line, cscd_out_name;
  all_cscd cscd_out_dis;
  info_ener_res infoener_res;
  bool topo_2_is_open{Arguments.flag_topo_2},
      topo_is_open(Arguments.flag_topo_1), status{true};

  groups_cal Groups_from_index_file;

  info info_output{Arguments};
  vector<traj_hb> res_hb_rec, res_hb_lg;
  vector<traj_sc> res_sc_rec, res_sc_lg;

  fstream out_sc(info_output.getSCout(), ios::out);
  /************************************/
  fstream cal_out;
  cal_out.open("cal.txt", ios::out);
  /*************************************/
  fstream log;
  log.open("log.txt", ios::out);
  log.close();
  log.open("log.txt", ios::out | ios::app);
  /*************************************/
  /******************************************/

  if (Arguments.com_dis_energy_2mol != 4) {
    for (size_t i{}; i < paths.get_size(); ++i) {
      fstream output_tot_dis(
          paths.get_path(i) + "/" + info_output.getDisCom2Mol(), ios::out);
      fstream output_elec_dis(paths.get_path(i) + "/electrostatic_" +
                                  info_output.getDisCom2Mol(),
                              ios::out);
      fstream output_vdw_dis(paths.get_path(i) + "/vanderWaals_" +
                                 info_output.getDisCom2Mol(),
                             ios::out);

      output_elec_dis.close();
      output_tot_dis.close();
      output_vdw_dis.close();
    }
  }

  if (Arguments.flag_h_bond == 1) {
    fstream out_h_bonds("out_hydrogen_bonds.dat", ios::out);
    out_h_bonds.close();
    fstream out_h_bonds_rec("out_hydrogen_bonds_receptor.dat", ios::out);
    out_h_bonds_rec << "# atom of receptor | atom of ligand\n";
    out_h_bonds_rec.close();
    fstream out_h_bonds_lig("out_hydrogen_bonds_ligand.dat", ios::out);
    out_h_bonds_lig << "# atom of ligand | atom of receptor\n";
    out_h_bonds_lig.close();
  }

  TitanUI::log_section(log, "TITAN ENGINE INITIALIZATION");
  TitanUI::log_kv(log, "Software Build", info_output.getBuild());

  time_t timeRun{time(nullptr)};
  char time_str[100];
  strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S",
           localtime(&timeRun));
  TitanUI::log_kv(log, "Execution Time", std::string(time_str));
  TitanUI::log_kv(log, "Force Field", force_f);
  TitanUI::log_kv(log, "Number of Threads", std::to_string(num_threads));

  TitanUI::log_section(log, "SYSTEM CONFIGURATION");
  TitanUI::log_kv(log, "Input File (-f)", Arguments.pdb_name);
  if (Arguments.flag_topo_1)
    TitanUI::log_kv(log, "Topology File (-topo)", Arguments.topo_name);
  if (Arguments.flag_topo_2)
    TitanUI::log_kv(log, "Topology File 2 (-topo2)", Arguments.topo_name_2);

  TitanUI::log_kv(log, "Molecule 1 Type", Arguments.mol_1_type);
  TitanUI::log_kv(log, "Molecule 2 Type", Arguments.mol_2_type);

  std::string group_str = std::to_string(Arguments.group_1) + " - " +
                          std::to_string(Arguments.group_2);
  TitanUI::log_kv(log, "Interaction Groups (-g)", group_str);

  size_t available_ram = get_available_ram();
  if (available_ram > 0) {
    TitanUI::log_kv(log, "Available RAM",
                    std::to_string(available_ram / (1024 * 1024)) + " MB");
  }

  // Strip trailing zeros for doubles using a lambda or just stringstream
  auto format_double = [](double v) {
    std::ostringstream oss;
    oss << v;
    return oss.str();
  };

  TitanUI::log_kv(log, "Distance Cutoff (-dsc)",
                  format_double(Arguments.d_standard) + " A");

  if (Arguments.flag_gbsa) {
    TitanUI::log_section(log, "GBSA CONFIGURATION");
    TitanUI::log_kv(log, "IGB Model", std::to_string(Arguments.igb));
    TitanUI::log_kv(log, "Internal Dielectric",
                    format_double(Arguments.int_diel));
    TitanUI::log_kv(log, "External Dielectric",
                    format_double(Arguments.ext_diel));
    TitanUI::log_kv(log, "Salt Concentration",
                    format_double(Arguments.saltcon) + " M");
  }

  std::string active_modules = "";
  if (Arguments.flag_cas)
    active_modules += "CAS ";
  if (Arguments.flag_network)
    active_modules += "DCCM/Network ";
  if (Arguments.flag_pca)
    active_modules += "PCA_FEL ";
  if (Arguments.flag_tda)
    active_modules += "TDA ";
  if (Arguments.msm.flag_tica)
    active_modules += "TICA ";
  if (Arguments.msm.flag_msm)
    active_modules += "MSM ";
  if (Arguments.flag_h_bond)
    active_modules += "H-Bonds ";
  if (Arguments.n_wat >= 0)
    active_modules += "Water ";
  if (Arguments.flag_new_emap)
    active_modules += "EMAP ";
  if (Arguments.flag_energy_same_force) {
    active_modules += "Energy/Force(=) ";
  }
  if (Arguments.flag_allostery)
    active_modules += "Allostery ";
  if (Arguments.lspmd.flag_lspmd)
    active_modules += "LSP-MD ";

  if (!active_modules.empty()) {
    TitanUI::log_section(log, "ACTIVE MODULES");
    TitanUI::log_kv(log, "Enabled Features", active_modules);
    // PBC algorithm in use
    if (Arguments.pbc.solv_pbcn)
      TitanUI::log_kv(log, "PBC Algorithm", "FixBox (Baptista et al. 2022)");
    else if (Arguments.pbc.solv_pbc_int == 1)
      TitanUI::log_kv(log, "PBC Algorithm",
                      "COM-heuristic (chain ref=" +
                          std::to_string(Arguments.pbc.chain_ref) + ")");
    else if (Arguments.pbc.solv_pbc == 1)
      TitanUI::log_kv(log, "PBC Algorithm", "Prob (legacy)");
    else
      TitanUI::log_kv(log, "PBC Algorithm", "None");
  }

  take_para(parameter, force_f, Arguments);
  take_amino_acids(amine, num_amine, force_f, Arguments);
  take_ff_para(amine, parameter);

  int input_name{};
  gen_input_list(Arguments, input_line, paths, log, input_name);

  // ── Auto-detect XDR trajectory from -i file list when -s was used ──
  // If -s set flag_xdr but no -traj was given, detect from input_line
  if (Arguments.flag_xdr && Arguments.traj_file.empty() &&
      !input_line.empty()) {
    const std::string &first_file = input_line[0];
    if (first_file.length() > 4) {
      std::string ext = first_file.substr(first_file.length() - 4);
      if (ext == ".xtc" || ext == ".trr") {
        Arguments.traj_file = first_file;
        Arguments.xdr_type = (ext == ".xtc") ? "xtc" : "trr";
        TitanUI::log_kv(log, "Auto-detect XDR",
                        Arguments.xdr_type + " from " + first_file);
      }
    }
  }

  GBSA_Aggregator global_gbsa_agg(Arguments, info_output);
  CASAnalyzer cas_analyzer;
  NucleicProfiler na_profiler;
  TDAAnalyzer tda_analyzer;
  bool global_has_nucleic_acid = false;

  // Network Analysis (DCCM) Variables
  DCCMNetwork dccm_analyzer;
  std::chrono::milliseconds time_gbsa{0}, time_cas{0}, time_network{0},
      time_tda{0}, time_llps{0}, time_tica{0};
  std::unique_ptr<TitanLLPSCluster> dbscan_llps;
  std::chrono::milliseconds time_llps_dbscan{0};

  WaterAnalyzer water_analyzer;
  MLFeatureExtractor ml_extractor;
  std::unique_ptr<TitanKinetics> tica_analyzer;
  std::unique_ptr<TitanMSM> msm_analyzer;
  std::unique_ptr<TitanLLPS> llps_analyzer;
  std::unique_ptr<TitanGIST> gist_analyzer;

  double time_new_emap{};
  std::chrono::milliseconds time_cal{}, time_emap{}, time_fmap{},
      time_subgroups_cal{}, time_solve_pbc{}, time_sc{}, time_hb{},
      time_other{}, time_print_ener_res{}, time_print_force{}, time_maxforce{},
      time_new_emap_other{}, time_avg_emap{}, time_contact{}, time_pdb{},
      time_cscd{0}, time_water{0}, time_na{0};
  vector<thread> threads(num_threads);
  std::vector<std::thread> threads_maxF(Arguments.nt);

  Multi_MOL pro_in(force_f), lig(force_f), water_traj(force_f);
  bool topo_initialized = false;
  bool charge_initialized = false;

  // --- IFP (Interaction Fingerprint) Initialization ---
  std::unique_ptr<TitanIFP> ifp_analyzer;
  if (Arguments.ifp.flag_ifp) {
    ifp_analyzer = std::make_unique<TitanIFP>();
    TitanUI::print_info("[IFP] Interaction Fingerprint module enabled.");
  }
  std::chrono::milliseconds time_ifp{0};

  // === TICA workspace — hoisted outside outer loop (Fix 1: single reference
  // frame shared across all replicas for H-REMD cross-replica consistency) ===
  int tica_n_ca_global = 0, tica_feat_dim_global = 0;
  Eigen::VectorXd tica_features_global;
  Eigen::MatrixXd tica_mobile_global, tica_ref_centered_global;
  Eigen::Matrix3d tica_R_global = Eigen::Matrix3d::Identity();
  bool tica_ref_initialized_global = false;

  while (true) {

    if (!input_line.empty()) {
      if (Arguments.flag_all_folders == false) {
        name_lig = "lg_" + input_line[input_name - Arguments.b_input] +
                   Arguments.get_input_file_type();
        name_pro = "rec_" + input_line[input_name - Arguments.b_input] +
                   Arguments.get_input_file_type();
      } else {

        name_lig = fs::path(input_line[input_name - Arguments.b_input])
                       .parent_path()
                       .string() +
                   "/lg_" +
                   fs::path(input_line[input_name - Arguments.b_input])
                       .stem()
                       .string() +
                   Arguments.get_input_file_type();
        name_pro = fs::path(input_line[input_name - Arguments.b_input])
                       .parent_path()
                       .string() +
                   "/rec_" +
                   fs::path(input_line[input_name - Arguments.b_input])
                       .stem()
                       .string() +
                   Arguments.get_input_file_type();
      }
    } else if (Arguments.flag_1_pdb == false &&
               Arguments.flag_all_folders == false && Arguments.flag_listfile) {
      name_lig = Arguments.input_pattern;
      name_pro = Arguments.input_pattern;
    } else {
      name_lig =
          "lg_" + to_string(input_name) + Arguments.get_input_file_type();
      name_pro =
          "rec_" + to_string(input_name) + Arguments.get_input_file_type();
    }

    h_groups groups_1, groups_2;
    vector<energy_res> Energy_res_mol1, Energy_res_mol2, Force_res_mol1,
        Force_res_mol2;
    unordered_map<pair<string, size_t>, size_t, pair_hash> map_energy_res1,
        map_energy_res2;
    // vector<vector<com> > COM_vec_rec, COM_vec_lg;
    vector<energy_groups> Energy_groups_out;
    vector<fstream> out_cal_dis_files;
    traj_hb hb_rec_1traj, hb_lg_1traj;
    traj_sc sc_rec_1traj, sc_lg_1traj;

    // Water Bridge Tracking
    if (Arguments.n_wat >= 0 && input_name == Arguments.b_input) {
      std::ofstream wb_out(info_output.get_global_water_bridges_timeseries());
      wb_out << "Traj_ID,Frame,Time_ps,Water_ID,Protein_Atom,Ligand_Atom,"
                "Distance_A,Angle_deg\n";
    }
    water_analyzer.reset_trajectory_state();

    multi_snap_emap rec_ener_inter;
    vector<vector<map_pos>> cont_matrix;
    vector<vector<int>> pos_lookup_1, pos_lookup_2;
    auto start = std::chrono::steady_clock::now();
    auto start_pdb{std::chrono::steady_clock::now()};

    // --- Total File Size ---
    size_t total_traj_bytes = 1;
    std::error_code ec;
    std::string ref_file = Arguments.flag_1_pdb ? Arguments.pdb_name : name_pro;
    if (std::filesystem::exists(ref_file, ec)) {
      total_traj_bytes = std::filesystem::file_size(ref_file, ec);
      if (total_traj_bytes == 0)
        total_traj_bytes = 1;
    }

    // Initialize the topology (Frame 0) and setup groups
    if (!topo_initialized) {
      auto init_result = initialize_topology_streaming(
          Arguments, log, name_pro, name_lig, force_f, input_name, input_line,
          topo_2_is_open, topo_is_open, pro_in, lig, water_traj,
          Groups_from_index_file, parameter, amine, info_output);

      if (!init_result) {
        if (Arguments.flag_nocal) {
          TitanUI::print_info("Nocal mode is turned on!");
          return 0; // nocal: init parsed what it could, done
        }
        if (input_name > Arguments.b_input)
          --input_name;
        status = false;
        TitanUI::print_error("Initialization failed: " + init_result.error());
        break;
      }
      topo_initialized = true;

    } else {
      // Subsequent files: ONLY reload Frame 0 coordinates, skip topo parsing
      pro_in.m_mol.clear();
      lig.m_mol.clear();
      water_traj.m_mol.clear();

      // Reload coordinates-only from new file (flag_raw=1 = coords only, no
      // topo)
      if (Arguments.flag_1_pdb) {
        const int r = pre_proc_pdb(force_f, pro_in, lig, water_traj,
                                   Groups_from_index_file, -1, Arguments, amine,
                                   info_output);
        if (r != 0) {
          TitanUI::print_error("Failed to reload coordinates for: " +
                               Arguments.pdb_name);
          break;
        }
      } else {
        pro_in.take_coordinate(name_pro, force_f, Arguments, amine, info_output,
                               1, 1);
        lig.take_coordinate(name_lig, force_f, Arguments, amine, info_output, 1,
                            1);
      }

      if (pro_in.m_mol.empty() || lig.m_mol.empty()) {
        TitanUI::print_error("Coordinate reload failed: empty mol for " +
                             name_pro);
        break;
      }

      if (Arguments.flag_1_pdb && !input_line.empty()) {
        Arguments.pdb_name = input_line[input_name - Arguments.b_input];
        Arguments.add_input_file_type();
      }

      Arguments.setPdbname(Arguments.flag_1_pdb
                               ? Arguments.pdb_name
                               : name_lig.substr(name_lig.find("_")));
      info_output.setPdbName(Arguments);

      TitanUI::log_kv(
          log, "Topo Cache",
          "Reusing topology, reloaded coords from: " +
              (Arguments.flag_1_pdb ? Arguments.pdb_name : name_pro));
    }

    info_output.get_info_group(Groups_from_index_file);
    process_time_of_infiles(pro_in, Arguments, info_output, time_of_snap);
    set_chainnum(pro_in, lig, Arguments);

    if (Arguments.opdb == true && Arguments.flag_sep_opdb == false) {
      fstream out_pdb(Arguments.get_pdbOutput(), ios::out);
      out_pdb.close();
    }

    if (Arguments.com_dis_2mol != 0) {
      fstream output_dis_2mol_file(
          paths.get_path(input_name - Arguments.b_input) + "/" +
              info_output.getDisCom2Mol(),
          ios::out);
      output_dis_2mol_file.close();
    }

    // --- CHARGE INITIALIZATION (Frame 0) ---
    if (!charge_initialized) {
      if (Arguments.flag_use_common_topo) {
        // Khi dùng common topo cho protein-protein, cả 2 phải dùng cờ 1
        // (chuẩn topology)
        pro_in.take_charge(amine, parameter, pro_in.m_mol[0].size(), 0ull,
                           force_f, 1);
        lig.take_charge(amine, parameter, lig.m_mol[0].size(), 0ull, force_f,
                        1);
      } else if (topo_is_open && topo_2_is_open) {
        pro_in.take_charge(amine, parameter, pro_in.m_mol[0].size(), 0ull,
                           force_f, 1);
        lig.take_charge(amine, parameter, lig.m_mol[0].size(), 0ull, force_f,
                        1);
      } else if (topo_is_open && Arguments.flag_topo_1 && !topo_2_is_open) {
        pro_in.take_charge(amine, parameter, pro_in.m_mol[0].size(), 0ull,
                           force_f, 0);
        lig.take_charge(amine, parameter, lig.m_mol[0].size(), 0ull, force_f,
                        1);
      } else if (!topo_is_open && topo_2_is_open) {
        pro_in.take_charge(amine, parameter, pro_in.m_mol[0].size(), 0ull,
                           force_f, 1);
        lig.take_charge(amine, parameter, lig.m_mol[0].size(), 0ull, force_f,
                        0);
      } else {
        pro_in.take_charge(amine, parameter, pro_in.m_mol[0].size(), 0ull,
                           force_f, 0);
        lig.take_charge(amine, parameter, lig.m_mol[0].size(), 0ull, force_f,
                        0);
      }
      charge_initialized = true;
    }

    time_pdb = std::chrono::duration_cast<milliseconds>(
        std::chrono::steady_clock::now() - start_pdb);

    cache_sigma_epsilon CachedData{Arguments.distance_energy_cutoff};
    // Martini CG systems are too large for O(N×M) pair cache (~GB).
    // Energy functions handle empty cache via elec_ener.empty() checks.
    // Downstream emap/force_map are guarded to skip when cache is empty.
    if (force_f != "MARTINI" &&
        (Arguments.flag_cache == 1 || Arguments.force_map == 1 ||
         Arguments.force_map == 3 || Arguments.flag_new_emap == true ||
         Arguments.energy_map == 1 || Arguments.energy_map == 3 ||
         Arguments.energy_map == 4)) {
      CachedData.set_size(pro_in.m_mol[0].size() * lig.m_mol[0].size(),
                          Arguments, pro_in.m_mol[0].size(),
                          lig.m_mol[0].size());
    }

    if (lig.m_mol.empty())
      lig.m_mol.resize(1);

    // --- MARTINI TOPOLOGY ASSIGNMENT (charges, LJ params, chain IDs) ---
    // Must run BEFORE total_charge() so that parameters[].charge is set.
    if (force_f == "MARTINI") {
      pro_in.parse_martini_topology(Arguments, log);
      lig.parse_martini_topology(Arguments, log);
    }

    double chg_pro = pro_in.total_charge();
    double chg_lig = lig.total_charge();

    TitanUI::log_section(log, "SYSTEM TOPOLOGY & TRAJECTORY");
    TitanUI::log_kv(log, "Receptor Group",
                    "Group 1 [" + std::to_string(Arguments.group_1) + "] (" +
                        std::to_string(pro_in.info.size()) + " atoms)");
    TitanUI::log_kv(log, "Ligand Group",
                    "Group 2 [" + std::to_string(Arguments.group_2) + "] (" +
                        std::to_string(lig.info.size()) + " atoms)");

    auto format_double = [](double v) {
      std::ostringstream oss;
      oss << v;
      return oss.str();
    };
    TitanUI::log_kv(log, "Receptor Charge", format_double(chg_pro) + " e");
    TitanUI::log_kv(log, "Ligand Charge", format_double(chg_lig) + " e");

    if (Arguments.de.flag_find_pulling_direct) {
      find_pull_direct(pro_in, lig, CachedData, DE_arguments, Arguments,
                       info_output, Groups_from_index_file);
      return 0;
    }

    const AlignedDoubleVector &rec_vdw_radii = pro_in.parameters.radius;
    const AlignedDoubleVector &lg_vdw_radii = lig.parameters.radius;

    if (Arguments.flag_cas)
      cas_analyzer.init(pro_in, lig, Arguments.mol_2_type);

    if (Arguments.flag_network || Arguments.flag_pca) {
      dccm_analyzer.init(pro_in, lig, Arguments);
      TitanUI::log_kv(log, "Network/PCA Nodes",
                      std::to_string(dccm_analyzer.get_resnrs().size()));
    }
    // LSP-MD standalone: init DCCM data collector without DCCM/PCA
    if (Arguments.lspmd.flag_lspmd && !Arguments.flag_network &&
        !Arguments.flag_pca) {
      dccm_analyzer.init(pro_in, lig, Arguments);
      TitanUI::log_kv(log, "LSP-MD Nodes",
                      std::to_string(dccm_analyzer.get_resnrs().size()));
    }

    // --- STREAMING ARCHITECTURE INITIALIZATION ---
    ChunkStreamer streamer(Arguments, info_output, Groups_from_index_file,
                           force_f, num_threads);
    // Save file positions after Frame 0 was consumed by init
    std::streampos pos_pro = -1, pos_lig = -1;
    {
      std::ifstream probe_pro(
          Arguments.flag_1_pdb ? Arguments.pdb_name : name_pro, std::ios::in);
      if (probe_pro.is_open()) {
        std::string tmp;
        while (std::getline(probe_pro, tmp)) {
          if (tmp.starts_with("ENDMDL") || tmp.starts_with("END")) {
            pos_pro = probe_pro.tellg();
            break;
          }
        }
      }
    }
    if (!Arguments.flag_1_pdb) {
      std::ifstream probe_lig(name_lig, std::ios::in);
      if (probe_lig.is_open()) {
        std::string tmp;
        while (std::getline(probe_lig, tmp)) {
          if (tmp.starts_with("ENDMDL") || tmp.starts_with("END")) {
            pos_lig = probe_lig.tellg();
            break;
          }
        }
      }
    }
    streamer.set_frame0_positions(pos_pro, pos_lig);

    auto stream_ok = streamer.open_streams(name_pro, name_lig, pro_in, lig);
    if (!stream_ok) {
      TitanUI::print_error(stream_ok.error());
      return 1;
    }

    size_t available_ram_bytes = get_available_ram();
    size_t CHUNK_SIZE = 100; // fallback default

    if (total_traj_bytes > 0 && available_ram_bytes > 0 &&
        total_traj_bytes <= available_ram_bytes * 0.5 && !Arguments.flag_xdr) {
      // File fits in memory AND is PDB — load everything at once
      CHUNK_SIZE = std::numeric_limits<size_t>::max();
      TitanUI::log_kv(log, "Memory Strategy",
                      "In-Memory Execution (File < 50% RAM)");
    } else if (Arguments.flag_xdr) {
      // XDR: use small fixed chunk — topology copy is expensive per frame,
      // and atom counts may not be available yet. XDR sequential read has
      // no penalty for small chunks.
      CHUNK_SIZE = 10;
      TitanUI::log_kv(log, "Memory Strategy",
                      "XDR Streaming (Chunk = 10 frames)");
    } else {
      // Adaptive chunk: target ~256 MB per chunk, clamped [10, 5000] frames
      size_t atoms_per_frame = pro_in.info.size() + lig.info.size();
      size_t bytes_per_frame = std::max(atoms_per_frame * 80, size_t(1));
      constexpr size_t TARGET_CHUNK_BYTES = 256ULL * 1024 * 1024; // 256 MB
      CHUNK_SIZE =
          std::max(size_t(10), std::min(size_t(5000),
                                        TARGET_CHUNK_BYTES / bytes_per_frame));
      TitanUI::log_kv(log, "Memory Strategy",
                      "Out-of-Core Streaming (Adaptive Chunk = " +
                          std::to_string(CHUNK_SIZE) + " frames, " +
                          std::to_string(bytes_per_frame / 1024) +
                          " KB/frame)");
    }

    // --- TICA INITIALIZATION (after DCCM init, uses Cα count) ---
    if (Arguments.msm.flag_tica && !tica_analyzer) {
      int n_ca = static_cast<int>(dccm_analyzer.get_resnrs().size());
      // Science fix: pairwise inter-Cα distances (rotation+translation
      // invariant) D = N_ca*(N_ca-1)/2 — eliminates 6 spurious DOF from raw xyz
      int tica_feature_dim = n_ca * (n_ca - 1) / 2;
      double tica_dt_ns = info_output.get_current_dt() > 0
                              ? info_output.get_current_dt() / 1000.0
                              : 0.001; // fallback 1 ps
      tica_analyzer = std::make_unique<TitanKinetics>(
          tica_feature_dim, Arguments.msm.tica_lag,
          Arguments.msm.tica_components, tica_dt_ns);
      TitanUI::log_kv(log, "TICA Cα Atoms", std::to_string(n_ca));
      TitanUI::log_kv(log, "TICA Feature Dim (pairwise dist)",
                      std::to_string(tica_feature_dim));
      TitanUI::log_kv(log, "TICA Lag Time",
                      std::to_string(Arguments.msm.tica_lag) + " frames");
      TitanUI::log_kv(log, "TICA Components",
                      std::to_string(Arguments.msm.tica_components));
    }

    // --- LLPS INITIALIZATION (after Martini topo parsed, chainnum set) ---
    if (Arguments.llps.flag_llps && !llps_analyzer) {
      // Build chain_start / chain_len from pro_in.info (chainnum assigned by
      // impl_martini)
      const auto &mol_info_vec = pro_in.info;
      const int n_info = static_cast<int>(mol_info_vec.size());
      int n_chains_llps = 0;
      for (int i = 0; i < n_info; ++i)
        n_chains_llps =
            std::max(n_chains_llps,
                     static_cast<int>(mol_info_vec[i].get_chainnum()) + 1);
      std::vector<int> llps_chain_start(n_chains_llps, -1);
      std::vector<int> llps_chain_len(n_chains_llps, 0);
      for (int i = 0; i < n_info; ++i) {
        int c = static_cast<int>(mol_info_vec[i].get_chainnum());
        if (c < 0 || c >= n_chains_llps)
          continue;
        if (llps_chain_start[c] < 0)
          llps_chain_start[c] = i;
        llps_chain_len[c]++;
      }
      llps_analyzer = std::make_unique<TitanLLPS>(
          n_chains_llps, llps_chain_start, llps_chain_len,
          Arguments.llps.llps_rcut, Arguments.llps.llps_perc_thr,
          Arguments.llps.llps_dimer_size, Arguments.llps.llps_contacts,
          Arguments.llps.llps_voxel_size);
      llps_analyzer->set_sticker_beads(mol_info_vec,
                                       Arguments.llps.llps_sticker_beads);
      TitanUI::log_kv(log, "LLPS Chains", std::to_string(n_chains_llps));
      TitanUI::log_kv(log, "LLPS Macro-Nodes",
                      std::to_string(llps_analyzer->n_macronodes()));
      TitanUI::log_kv(log, "LLPS Dimer Size",
                      std::to_string(Arguments.llps.llps_dimer_size));
      TitanUI::log_kv(log, "LLPS Rcut (nm)",
                      std::to_string(Arguments.llps.llps_rcut));
      // Start dedicated LLPS worker thread (pipeline mode)
      llps_analyzer->start_async();
    }

    // --- LLPS DBSCAN (bead-level morphology) ---
    if (Arguments.llps.flag_llps && !dbscan_llps) {
      dbscan_llps = std::make_unique<TitanLLPSCluster>();
      const std::string pdb_out =
          info_output.get_llps_dbscan_droplets_filename();
      const std::string csv_out =
          info_output.get_llps_dbscan_timeseries_filename();
      // [BUG 1 fix] Use user-supplied eps (Å) and min_pts; default 15 Å / 5
      // beads. Pass csv_out so data is flushed per-frame (crash-safe, ISSUE 4
      // fix).
      TitanUI::log_kv(log, "LLPS-DBSCAN eps",
                      std::to_string(Arguments.llps.llps_dbscan_eps) + " A");
      TitanUI::log_kv(log, "LLPS-DBSCAN min_pts",
                      std::to_string(Arguments.llps.llps_dbscan_minpts));
      dbscan_llps->init(Arguments.llps.llps_dbscan_eps,
                        Arguments.llps.llps_dbscan_minpts, pdb_out, csv_out);
      if (Arguments.llps.llps_fuzzy) {
        dbscan_llps->enable_fuzzy(true);
        TitanUI::print_info("[LLPS-Fuzzy] Fuzzy EM refinement enabled.");
      }
    }

    // --- GIST Water Thermodynamics Initialization ---
    if (Arguments.gist.flag_gist) {
      // Auto-enable water tracking if not explicitly set
      if (Arguments.n_wat < 0) {
        Arguments.n_wat = 0;
        TitanUI::print_info("[GIST] Auto-enabling water tracking (-n_wat 0) "
                            "for GIST analysis.");
      }

      if (!gist_analyzer && !lig.m_mol.empty() &&
          lig.m_mol[0].get_pos_x().size() > 0) {
        gist_analyzer = std::make_unique<TitanGIST>();
        gist_analyzer->init(
            lig.m_mol[0], Arguments.gist.gist_padding,
            Arguments.gist.gist_resolution, Arguments.gist.gist_temperature,
            Arguments.gist.gist_cutoff, Arguments.nt > 0 ? Arguments.nt : 1,
            Arguments.gist.gist_water_model, pro_in.parameters, pro_in.info,
            lig.parameters, lig.info);
        // BUG2: pass DCCM's exact Cα indices for TICA centroid
        if (!dccm_analyzer.get_ca_indices().empty())
          gist_analyzer->set_tica_ca_indices(dccm_analyzer.get_ca_indices());
        TitanUI::log_kv(log, "GIST Water Model",
                        Arguments.gist.gist_water_model);
        TitanUI::log_kv(log, "GIST Resolution",
                        std::to_string(Arguments.gist.gist_resolution) + " A");
        TitanUI::log_kv(log, "GIST Voxels",
                        std::to_string(gist_analyzer->n_voxels()));
        TitanUI::log_kv(log, "GIST Cutoff",
                        std::to_string(Arguments.gist.gist_cutoff) + " A");
      }
    }

    // --- MSM INITIALIZATION (after TICA init, shares n_components) ---
    if (Arguments.msm.flag_msm && tica_analyzer && !msm_analyzer) {
      double msm_dt_ns = info_output.get_current_dt() > 0
                             ? info_output.get_current_dt() / 1000.0
                             : 0.001;
      msm_analyzer = std::make_unique<TitanMSM>(
          Arguments.msm.msm_clusters, tica_analyzer->n_components(), msm_dt_ns);
      TitanUI::log_kv(log, "MSM Clusters",
                      std::to_string(Arguments.msm.msm_clusters));
      TitanUI::log_kv(log, "MSM Lag",
                      std::to_string(Arguments.msm.msm_lag) + " frames");
      TitanUI::log_kv(log, "MSM PCCA",
                      Arguments.msm.msm_run_pcca ? "enabled" : "disabled");
    }

    // Fix 5: hoist TICA workspace OUTSIDE streaming loops
    // (pairwise distance vector + Kabsch mobile/reference matrices — zero heap
    // per frame)
    // Fix 1: these are aliases to global variables hoisted outside while(true)
    // so all replicas share a single Kabsch reference frame.
    if (tica_n_ca_global == 0 && Arguments.msm.flag_tica && tica_analyzer) {
      tica_n_ca_global = static_cast<int>(dccm_analyzer.get_resnrs().size());
      tica_feat_dim_global = tica_analyzer->feature_dim();
      if (tica_n_ca_global > 0) {
        tica_features_global.resize(tica_feat_dim_global);
        tica_mobile_global.resize(tica_n_ca_global, 3);
        tica_ref_centered_global.resize(tica_n_ca_global, 3);
      }
    }
    const int& tica_n_ca       = tica_n_ca_global;
    const int& tica_feat_dim   = tica_feat_dim_global;
    auto& tica_features        = tica_features_global;
    auto& tica_mobile          = tica_mobile_global;
    auto& tica_ref_centered    = tica_ref_centered_global;
    auto& tica_R               = tica_R_global;
    bool& tica_ref_initialized = tica_ref_initialized_global;


    size_t global_frame_counter = 0;
    bool is_streaming = true;
    size_t NUM = 0;

    // ========================================================================
    // Save Frame 0 topology+coords for post-loop functions that need m_mol
    // (e.g. print_energy_pdbtype which writes B-factor PDB).
    // The streaming loop clears pro_in/lig.m_mol each iteration.
    // ========================================================================
    const std::vector<MOL> saved_pro_frame0 = pro_in.m_mol;
    const std::vector<MOL> saved_lig_frame0 = lig.m_mol;

    // ========================================================================
    // MAIN STREAMING LOOP
    // ========================================================================
    while (is_streaming) {
      auto expected_frames =
          streamer.read_next_chunk(CHUNK_SIZE, pro_in, lig, water_traj);
      if (!expected_frames) {
        TitanUI::print_error(expected_frames.error());
        break;
      }

      size_t chunk_num_snaps = expected_frames.value();
      if (chunk_num_snaps == 0) {
        is_streaming = false;
        break;
      }

      NUM = chunk_num_snaps;

      // --- PBC & WATER SHIFT LOGIC (Applied Per-Chunk) ---
      struct shift_vector {
        double x, y, z;
      };
      std::vector<shift_vector> old_cogs(pro_in.m_mol.size(), {0.0, 0.0, 0.0});
      if (Arguments.n_wat >= 0) {
        for (size_t snap = 0; snap < pro_in.m_mol.size(); ++snap) {
          double old_cx = 0.0, old_cy = 0.0, old_cz = 0.0;
          size_t p_size = pro_in.m_mol[snap].pos.x.size();
          if (p_size > 0) {
            for (size_t i = 0; i < p_size; ++i) {
              old_cx += pro_in.m_mol[snap].pos.x[i];
              old_cy += pro_in.m_mol[snap].pos.y[i];
              old_cz += pro_in.m_mol[snap].pos.z[i];
            }
            old_cx /= p_size;
            old_cy /= p_size;
            old_cz /= p_size;
          }
          old_cogs[snap] = {old_cx, old_cy, old_cz};
        }
      }

      auto start_pbc = std::chrono::steady_clock::now();
      if (Arguments.pbc.solv_pbcn) {
        pro_in.solv_pbc_pbcn(lig, Arguments);
        align_complex(pro_in, lig, 0, 0);
      } else if (Arguments.pbc.solv_pbc_int == 1) {
        if (Arguments.pbc.chain_ref == -2) {
          for (size_t i{}; i < num_threads; ++i)
            threads[i] =
                thread(&Multi_MOL::solve_pbc_int, &pro_in, i, num_threads);
          for (auto &t : threads)
            if (t.joinable())
              t.join();
        }
        if (Arguments.pbc.chain_ref > -2) {
          pro_in.solve_pbc_new(lig, Arguments, amine);
          align_complex(pro_in, lig, 0, 0);
        }
      }
      if (Arguments.pbc.solv_pbc == 1) {
        for (size_t i{}; i < num_threads; ++i)
          threads[i] =
              thread(&Multi_MOL::pbc_prob, &pro_in, ref(lig), i, num_threads);
        for (auto &t : threads)
          if (t.joinable())
            t.join();
        align_complex(pro_in, lig, 0, 0);
      }
      time_solve_pbc += std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::steady_clock::now() - start_pbc);

      // ── NOCAL MODE: PBC only → write PDB/PQR → skip energy ──────────
      if (Arguments.flag_nocal) {
        // Use align_to_box (position atoms within box) instead of
        // align_complex (which aligns two molecules for interaction calc)
        if (Arguments.pbc.solv_pbcn || Arguments.pbc.solv_pbc_int)
          pro_in.align_to_box();

        // PBC clustering (non-int mode)
        if (Arguments.pbc.solv_pbc_clust && !Arguments.pbc.solv_pbc_int)
          pro_in.solve_pbc_cluster(Arguments, info_output);

        TitanUI::log_kv(
            log, "Nocal Chunk",
            "Writing frames " + std::to_string(global_frame_counter + 1) +
                " to " +
                std::to_string(global_frame_counter + chunk_num_snaps));

        for (size_t chunk_snap_idx = 0; chunk_snap_idx < chunk_num_snaps;
             ++chunk_snap_idx) {
          // COM distance between molecules
          if (Arguments.com_dis_2mol != 0) {
            fstream output_dis_2mol_file(
                paths.get_path(input_name - Arguments.b_input) + "/" +
                    info_output.getDisCom2Mol(),
                ios::out | ios::app);
            pro_in.impl_com_dis_2mol(lig, Arguments, info_output,
                                     output_dis_2mol_file, chunk_snap_idx);
          }

          // Dipole moment
          if (Arguments.flag_cal_dipole_moment)
            pro_in.cal_dipole_moment(
                info_output, Arguments,
                Arguments.pdb_name.substr(0, Arguments.pdb_name.length() - 4));

          // Charged atom distance distribution
          if (Arguments.flag_charged_atom_distance_dis)
            pro_in.cal_distribution_distance_charged_atoms(Arguments,
                                                           info_output);

          // Volume calculation
          if (Arguments.flag_cal_vol) {
            vector<double> volumes_temp;
            pro_in.cal_volumes(volumes_temp);
            Volumes.emplace_back(std::move(volumes_temp));
          }

          // H-bond (self-interaction mode for nocal)
          if (Arguments.flag_h_bond == 1) {
            h_groups groups_1_nc, groups_2_nc;
            MOL temp_mol(pro_in.m_mol[chunk_snap_idx]);
            find_h_bonds(pro_in.m_mol[chunk_snap_idx], pro_in.m_mol[0],
                         temp_mol, temp_mol, groups_1_nc, groups_2_nc,
                         Arguments, chunk_snap_idx, pro_in.info, pro_in.info,
                         Hbonds_vector);
          }

          // PQR output
          if (Arguments.flag_pqr == 1)
            write_pqr(pro_in, lig, input_line, input_name, chunk_snap_idx,
                      Arguments, paths);

          // PDB output (critical for nocal)
          write_pdb(pro_in, lig, Arguments, chunk_snap_idx, paths, input_line,
                    input_name);
        }

        tot_snap_allfiles += chunk_num_snaps;
        global_frame_counter += chunk_num_snaps;
        continue; // skip energy calculations, read next chunk
      }
      // ── END NOCAL ───────────────────────────────────────────────────────

      if (Arguments.entropy_quasi) {
        if (Arguments.pbc.solv_pbc == 0 && Arguments.pbc.chain_ref == -1)
          align_complex(pro_in, lig, 0, 0);
        cal_entropy(pro_in, lig, Groups_from_index_file, Arguments, paths,
                    static_cast<size_t>(input_name), info_output);
      }

      if (Arguments.n_wat >= 0 && !water_traj.m_mol.empty()) {
        size_t limit = std::min(pro_in.m_mol.size(), water_traj.m_mol.size());
        for (size_t snap = 0; snap < limit; ++snap) {
          size_t p_size = pro_in.m_mol[snap].pos.x.size();
          if (p_size > 0) {
            double new_cx = 0.0, new_cy = 0.0, new_cz = 0.0;
            for (size_t i = 0; i < p_size; ++i) {
              new_cx += pro_in.m_mol[snap].pos.x[i];
              new_cy += pro_in.m_mol[snap].pos.y[i];
              new_cz += pro_in.m_mol[snap].pos.z[i];
            }
            new_cx /= p_size;
            new_cy /= p_size;
            new_cz /= p_size;
            double shift_x = new_cx - old_cogs[snap].x;
            double shift_y = new_cy - old_cogs[snap].y;
            double shift_z = new_cz - old_cogs[snap].z;
            size_t w_size = water_traj.m_mol[snap].pos.x.size();
            for (size_t i = 0; i < w_size; ++i) {
              water_traj.m_mol[snap].pos.x[i] += shift_x;
              water_traj.m_mol[snap].pos.y[i] += shift_y;
              water_traj.m_mol[snap].pos.z[i] += shift_z;
            }
          }
        }
      }

      TitanUI::log_kv(
          log, "Processing Chunk",
          "Frames " + std::to_string(global_frame_counter + 1) + " to " +
              std::to_string(global_frame_counter + chunk_num_snaps));

      // --- INNER FRAME LOOP ---
      for (size_t chunk_snap_idx = 0; chunk_snap_idx < chunk_num_snaps;
           ++chunk_snap_idx) {
        size_t global_snap_idx = global_frame_counter + chunk_snap_idx;

        // --- LLPS per-frame analysis (async pipeline — overlaps energy calc)
        // ---
        if (Arguments.llps.flag_llps && llps_analyzer) {
          auto start_llps = std::chrono::steady_clock::now();
          const auto &mol_frame = pro_in.m_mol[chunk_snap_idx];
          // Box dimensions from CRYST1 (Å → nm)
          const float bx_nm = static_cast<float>(mol_frame.Box.get_a() / 10.0);
          const float by_nm = static_cast<float>(mol_frame.Box.get_b() / 10.0);
          const float bz_nm = static_cast<float>(mol_frame.Box.get_c() / 10.0);
          const float t_ns =
              static_cast<float>(mol_frame.time_of_snap / 1000.0);
          // Pass RNA (lig) data if mol_2_type is rna
          const MOL *rna_ptr = nullptr;
          const std::vector<mol_info> *rna_info_ptr = nullptr;
          if (Arguments.mol_2_type.find("rna") != std::string::npos ||
              Arguments.mol_2_type.find("RNA") != std::string::npos) {
            if (chunk_snap_idx < lig.m_mol.size()) {
              rna_ptr = &lig.m_mol[chunk_snap_idx];
              rna_info_ptr = &lig.info;
            }
          }
          // Submit to async queue — returns immediately, runs on worker thread
          llps_analyzer->process_frame_async(
              mol_frame, pro_in.info, bx_nm, by_nm, bz_nm, t_ns,
              static_cast<int>(global_snap_idx), rna_ptr, rna_info_ptr);
          time_llps += std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::steady_clock::now() - start_llps);
        }

        // --- LLPS DBSCAN per-frame ---
        if (dbscan_llps) {
          auto t0_dbscan = std::chrono::steady_clock::now();
          const MOL *rna_mol = nullptr;
          const std::vector<mol_info> *rna_info = nullptr;
          if (chunk_snap_idx < lig.m_mol.size()) {
            rna_mol = &lig.m_mol[chunk_snap_idx];
            rna_info = &lig.info;
          }
          static const MOL empty_mol;
          static const std::vector<mol_info> empty_info;
          // [BUG 2 note] pro.Box is used as the shared periodic box for all
          // nodes (protein + RNA/lig). This is correct when both molecules come
          // from the same trajectory frame. If lig were ever loaded from a
          // separate file with a different box, MIC for RNA nodes would be
          // wrong — document that here.
          dbscan_llps->process_frame(pro_in.m_mol[chunk_snap_idx], pro_in.info,
                                     rna_mol ? *rna_mol : empty_mol,
                                     rna_info ? *rna_info : empty_info,
                                     pro_in.m_mol[chunk_snap_idx].time_of_snap,
                                     static_cast<int>(global_snap_idx + 1));
          time_llps_dbscan +=
              std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - t0_dbscan);
        }

        std::string active_in_loop = "(";
        if (Arguments.flag_gbsa)
          active_in_loop += "GBSA ";
        if (Arguments.flag_cas)
          active_in_loop += "CAS ";
        if (global_has_nucleic_acid)
          active_in_loop += "RNA_Profile ";
        if (Arguments.flag_network)
          active_in_loop += "DCCM ";
        if (Arguments.flag_pca)
          active_in_loop += "PCA ";
        if (Arguments.flag_tda)
          active_in_loop += "TDA ";
        if (Arguments.msm.flag_tica)
          active_in_loop += "TICA ";
        if (Arguments.dham.flag_dham)
          active_in_loop += "DHAM ";
        if (Arguments.llps.flag_llps)
          active_in_loop += "LLPS ";
        if (Arguments.ifp.flag_ifp)
          active_in_loop += "IFP ";
        if (Arguments.lspmd.flag_lspmd)
          active_in_loop += "LSP-MD ";
        std::string extra_info =
            active_in_loop + ") | Total frame: " +
            std::to_string((tot_snap_allfiles + global_snap_idx + 1));

        size_t atoms_per_frame = pro_in.info.size() + lig.info.size();
        size_t estimated_bytes_per_frame = atoms_per_frame * 80;
        size_t estimated_current_bytes =
            global_snap_idx * estimated_bytes_per_frame;

        size_t display_current =
            std::min(estimated_current_bytes,
                     total_traj_bytes > 1 ? total_traj_bytes - 1 : 0);

        std::filesystem::path current_pdb =
            Arguments.flag_1_pdb ? Arguments.pdb_name : name_pro;
        std::string module_name = current_pdb.filename().string();
        if (module_name.length() > 20)
          module_name = module_name.substr(0, 17) + "...";

        if (Arguments.flag_xdr || Arguments.flag_begin_set ||
            Arguments.flag_end_set || Arguments.dt > 0) {
          // XDR or filtered: total frames unknown, use continuous counter
          TitanUI::show_continuous_progress(module_name, global_snap_idx + 1,
                                            extra_info);
        } else {
          TitanUI::show_progress(
              module_name, display_current, total_traj_bytes,
              "Frame: " + std::to_string(global_snap_idx + 1) + " " +
                  extra_info);
        }

        vector<d_min> dis_min(num_threads);
        vector<cache_energy> Energy(num_threads);
        vector<cache_force> Force(num_threads);

        double E_V_{}, E_E_{};
        force f_e_{}, f_v_{};

        // Print formatted timestamps in a neat grid (10 columns)
        if (log.is_open()) {
          log << std::fixed << std::setprecision(1) << std::setw(10)
              << std::right << pro_in.m_mol[chunk_snap_idx].time_of_snap;
          if ((chunk_snap_idx + 1) % 10 == 0 || chunk_snap_idx + 1 == NUM) {
            log << "\n";
          }
        }

        AlignedDoubleVector rec_free_radii;
        AlignedDoubleVector lg_free_radii;
        AlignedDoubleVector rec_complex_radii;
        AlignedDoubleVector lg_complex_radii;
        size_t lg_snap = (lig.m_mol.size() > 1) ? chunk_snap_idx : 0;
        const MOL &active_lg_mol = lig.m_mol[lg_snap];
        const vector<mol_info> &active_lg_info = lig.info;

        auto start_energy_eval = std::chrono::steady_clock::now();
        GBSA_energy_snap frame_gbsa_total;
        if (Arguments.flag_gbsa) {
#ifdef __AVX__
          rec_free_radii.resize(pro_in.m_mol[chunk_snap_idx].size());
          lg_free_radii.resize(active_lg_mol.size());
          rec_complex_radii.resize(pro_in.m_mol[chunk_snap_idx].size());
          lg_complex_radii.resize(active_lg_mol.size());

          std::thread t1(cal_effective_born_radii_avx2,
                         std::ref(pro_in.m_mol[chunk_snap_idx]),
                         std::ref(rec_vdw_radii), std::ref(rec_free_radii),
                         Arguments.igb, 0, pro_in.m_mol[chunk_snap_idx].size(),
                         nullptr, nullptr);
          std::thread t2;
          std::thread t3;
          std::thread t4;

          if (active_lg_mol.size() > 0) {
            t2 = std::thread(cal_effective_born_radii_avx2,
                             std::ref(active_lg_mol), std::ref(lg_vdw_radii),
                             std::ref(lg_free_radii), Arguments.igb, 0,
                             active_lg_mol.size(), nullptr, nullptr);
            t3 = std::thread(cal_effective_born_radii_avx2,
                             std::ref(pro_in.m_mol[chunk_snap_idx]),
                             std::ref(rec_vdw_radii),
                             std::ref(rec_complex_radii), Arguments.igb, 0,
                             pro_in.m_mol[chunk_snap_idx].size(),
                             &active_lg_mol, &lg_vdw_radii);
            t4 = std::thread(cal_effective_born_radii_avx2,
                             std::ref(active_lg_mol), std::ref(lg_vdw_radii),
                             std::ref(lg_complex_radii), Arguments.igb, 0,
                             active_lg_mol.size(),
                             &pro_in.m_mol[chunk_snap_idx], &rec_vdw_radii);
          }

          if (t1.joinable())
            t1.join();
          if (active_lg_mol.size() > 0) {
            if (t2.joinable())
              t2.join();
            if (t3.joinable())
              t3.join();
            if (t4.joinable())
              t4.join();
          }
#else
          std::cerr << "Error: AVX2 is required for GBSA!" << std::endl;
          exit(1);
#endif

          std::vector<GBSA_Aggregator> local_aggs(
              num_threads, GBSA_Aggregator(Arguments, info_output));
          std::vector<GBSA_energy_snap> gbsa_energies(num_threads);
          std::vector<info_ener_res> infoener_res_v(num_threads);

          size_t active_lg_size = active_lg_mol.size();
          size_t rec_size = pro_in.m_mol[chunk_snap_idx].size();
          for (size_t i{}; i < num_threads; ++i) {
#ifdef __AVX__
            threads[i] = thread(
                cal_energy_gbsa_avx2, std::ref(pro_in.m_mol[chunk_snap_idx]),
                std::ref(active_lg_mol), std::ref(pro_in.info),
                std::ref(active_lg_info), std::ref(CachedData),
                rec_size * i / num_threads, rec_size * (i + 1) / num_threads,
                active_lg_size * i / num_threads,
                active_lg_size * (i + 1) / num_threads,
                std::ref(gbsa_energies[i]), chunk_snap_idx, std::ref(Arguments),
                std::ref(pro_in.parameters), std::ref(lig.parameters),
                std::ref(pro_in.lcpo_p1), std::ref(pro_in.lcpo_p2),
                std::ref(lig.lcpo_p1), std::ref(lig.lcpo_p2),
                std::ref(rec_free_radii), std::ref(lg_free_radii),
                std::ref(rec_complex_radii), std::ref(lg_complex_radii),
                std::ref(local_aggs[i]), &infoener_res_v[i]);
#endif
          }
          for (auto &t : threads)
            if (t.joinable())
              if (t.joinable())
                t.join();

          if (Arguments.flag_energy_same_force == false &&
              (Arguments.energy_map == 1 || Arguments.energy_map == 3 ||
               Arguments.energy_map == 4)) {
            for (auto &ii : infoener_res_v) {
              for (auto &ele : ii.position_map) {
                infoener_res.add_energy(
                    ele.first, ii.emap_ee[ele.second],
                    ii.emap_ev[ele.second]);
              }
            }
            infoener_res.wrapper_impl_ener_maptype(
                pro_in, lig, info_output, Arguments, chunk_snap_idx + 1);
          }

          // Merge thread-local GBSA results back to the main pipeline
          frame_gbsa_total.time = pro_in.m_mol[chunk_snap_idx].time_of_snap;

          // double current_E_vdw = 0.0;
          // double current_E_elec = 0.0;

          GBSA_Aggregator frame_decomp(
              Arguments,
              info_output); // GÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â»ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¾ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢p
                            // cÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡c
                            // luÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â»ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¦ÃƒÂ¢Ã¢â€šÂ¬Ã…â€œng
                            // cÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â»ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â§a
                            // frame
                            // nÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â y

          for (size_t i{}; i < num_threads; ++i) {
            Energy[i].vdw_ener = gbsa_energies[i].E_vdw;
            Energy[i].elec_ener = gbsa_energies[i].E_elec;

            frame_gbsa_total.E_vdw += gbsa_energies[i].E_vdw;
            frame_gbsa_total.E_elec += gbsa_energies[i].E_elec;
            frame_gbsa_total.G_polar += gbsa_energies[i].G_polar;
            frame_gbsa_total.G_nonpolar += gbsa_energies[i].G_nonpolar;
            frame_gbsa_total.total_energy += gbsa_energies[i].total_energy;

            frame_decomp.merge(
                local_aggs
                    [i]); // GÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â»ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¾ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¢p
                          // Thread -> Frame
          }

          global_gbsa_agg.commit_frame_decompositions(
              frame_decomp, tot_snap_allfiles + global_snap_idx);
          global_gbsa_agg.add_snapshot(input_name - Arguments.b_input,
                                       frame_gbsa_total);

          time_gbsa += std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::steady_clock::now() - start_energy_eval);
        } else {

          if (Arguments.force_f.find("CHARMM") != string::npos ||
              Arguments.force_f.find("AMBER") != string::npos) {
            for (size_t i{}; i < num_threads; ++i) {
              if (Arguments.flag_elec_screening == false) {
                if (Arguments.flag_cal_force == false) {
#ifdef __AVX__
                  threads[i] =
                      thread(cal_energy_amber_intrin,
                             ref(pro_in.m_mol[chunk_snap_idx]),
                             ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                             ref(lig.info), ref(CachedData),
                             lig.m_mol[0].size() * i / num_threads,
                             lig.m_mol[0].size() * (i + 1) / num_threads,
                             ref(Energy[i]), chunk_snap_idx, ref(Arguments),
                             ref(pro_in.parameters), ref(lig.parameters));
#else
                  threads[i] = thread(
                      cal_energy_amber, ref(pro_in.m_mol[chunk_snap_idx]),
                      ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                      ref(lig.info), ref(CachedData),
                      lig.m_mol[0].size() * i / num_threads,
                      lig.m_mol[0].size() * (i + 1) / num_threads,
                      ref(Energy[i]), chunk_snap_idx, ref(Arguments),
                      ref(pro_in.parameters), ref(lig.parameters));
#endif
                } else
                  threads[i] = thread(
                      cal_energy_force_amber, ref(pro_in.m_mol[chunk_snap_idx]),
                      ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                      ref(lig.info), ref(CachedData),
                      lig.m_mol[0].size() * i / num_threads,
                      lig.m_mol[0].size() * (i + 1) / num_threads,
                      ref(Energy[i]), ref(Force[i]), chunk_snap_idx,
                      ref(Arguments), ref(pro_in.parameters),
                      ref(lig.parameters));
              } else {
                if (Arguments.flag_cal_force == false)

#ifdef __AVX__
                  threads[i] =
                      thread(cal_energy_amber_screening_intrin,
                             ref(pro_in.m_mol[chunk_snap_idx]),
                             ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                             ref(lig.info), ref(CachedData),
                             lig.m_mol[0].size() * i / num_threads,
                             lig.m_mol[0].size() * (i + 1) / num_threads,
                             ref(Energy[i]), chunk_snap_idx, ref(Arguments),
                             ref(pro_in.parameters), ref(lig.parameters));
#else
                  threads[i] =
                      thread(cal_energy_amber_screening,
                             ref(pro_in.m_mol[chunk_snap_idx]),
                             ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                             ref(lig.info), ref(CachedData),
                             lig.m_mol[0].size() * i / num_threads,
                             lig.m_mol[0].size() * (i + 1) / num_threads,
                             ref(Energy[i]), chunk_snap_idx, ref(Arguments),
                             ref(pro_in.parameters), ref(lig.parameters));
#endif
                else
                  threads[i] = thread(
                      cal_energy_force_amber, ref(pro_in.m_mol[chunk_snap_idx]),
                      ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                      ref(lig.info), ref(CachedData),
                      lig.m_mol[0].size() * i / num_threads,
                      lig.m_mol[0].size() * (i + 1) / num_threads,
                      ref(Energy[i]), ref(Force[i]), chunk_snap_idx,
                      ref(Arguments), ref(pro_in.parameters),
                      ref(lig.parameters));
              }
            }
            for (auto &t : threads)
              if (t.joinable())
                t.join();

          } else if (Arguments.force_f.find("OPLSAA") != string::npos) {
            for (size_t i{}; i < num_threads; ++i) {
              if (Arguments.flag_cal_force == false)
#ifdef __AVX__
                threads[i] = thread(
                    cal_energy_opls_intrin, ref(pro_in.m_mol[chunk_snap_idx]),
                    ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                    ref(lig.info), ref(CachedData),
                    lig.m_mol[0].size() * i / num_threads,
                    lig.m_mol[0].size() * (i + 1) / num_threads, ref(Energy[i]),
                    chunk_snap_idx, ref(Arguments), ref(pro_in.parameters),
                    ref(lig.parameters));
#else
                threads[i] =
                    thread(cal_energy_opls, ref(pro_in.m_mol[chunk_snap_idx]),
                           ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                           ref(lig.info), ref(CachedData),
                           lig.m_mol[0].size() * i / num_threads,
                           lig.m_mol[0].size() * (i + 1) / num_threads,
                           ref(Energy[i]), chunk_snap_idx, ref(Arguments),
                           ref(pro_in.parameters), ref(lig.parameters));
#endif
              else
                threads[i] = thread(
                    cal_energy_force_opls, ref(pro_in.m_mol[chunk_snap_idx]),
                    ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                    ref(lig.info), ref(CachedData),
                    lig.m_mol[0].size() * i / num_threads,
                    lig.m_mol[0].size() * (i + 1) / num_threads, ref(Energy[i]),
                    ref(Force[i]), chunk_snap_idx, ref(Arguments),
                    ref(pro_in.parameters), ref(lig.parameters));
            }
            for (auto &t : threads)
              if (t.joinable())
                if (t.joinable())
                  t.join();
          } else if (Arguments.force_f.find("MARTINI") != string::npos) {
            for (size_t i{}; i < num_threads; ++i) {
              if (Arguments.flag_cal_force == false)
#ifdef __AVX__
                threads[i] =
                    thread(cal_energy_martini_intrin,
                           ref(pro_in.m_mol[chunk_snap_idx]),
                           ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                           ref(lig.info), ref(CachedData),
                           lig.m_mol[0].size() * i / num_threads,
                           lig.m_mol[0].size() * (i + 1) / num_threads,
                           ref(Energy[i]), chunk_snap_idx, ref(Arguments),
                           ref(pro_in.parameters), ref(lig.parameters));
#else
                threads[i] = thread(
                    cal_energy_martini, ref(pro_in.m_mol[chunk_snap_idx]),
                    ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                    ref(lig.info), ref(CachedData),
                    lig.m_mol[0].size() * i / num_threads,
                    lig.m_mol[0].size() * (i + 1) / num_threads, ref(Energy[i]),
                    chunk_snap_idx, ref(Arguments), ref(pro_in.parameters),
                    ref(lig.parameters));
#endif
              else
                threads[i] = thread(
                    cal_energy_martini, ref(pro_in.m_mol[chunk_snap_idx]),
                    ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                    ref(lig.info), ref(CachedData),
                    lig.m_mol[0].size() * i / num_threads,
                    lig.m_mol[0].size() * (i + 1) / num_threads, ref(Energy[i]),
                    chunk_snap_idx, ref(Arguments), ref(pro_in.parameters),
                    ref(lig.parameters));
              // NOTE: cal_energy_force_martini not yet implemented.
              // Fall back to scalar for flag_cal_force=true path.
            }
            for (auto &t : threads)
              if (t.joinable())
                t.join();
          } else {
            for (size_t i{}; i < num_threads; ++i) {
              if (Arguments.flag_cal_force == false)
#ifdef __AVX__
                threads[i] = thread(
                    cal_energy_gromos_intrin, ref(pro_in.m_mol[chunk_snap_idx]),
                    ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                    ref(lig.info), ref(CachedData),
                    lig.m_mol[0].size() * i / num_threads,
                    lig.m_mol[0].size() * (i + 1) / num_threads, ref(Energy[i]),
                    chunk_snap_idx, ref(Arguments), ref(pro_in.parameters),
                    ref(lig.parameters));
#else
                threads[i] =
                    thread(cal_energy_gromos, ref(pro_in.m_mol[chunk_snap_idx]),
                           ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                           ref(lig.info), ref(CachedData),
                           lig.m_mol[0].size() * i / num_threads,
                           lig.m_mol[0].size() * (i + 1) / num_threads,
                           ref(Energy[i]), chunk_snap_idx, ref(Arguments),
                           ref(pro_in.parameters), ref(lig.parameters));
#endif
              else
                threads[i] = thread(
                    cal_energy_force_gromos, ref(pro_in.m_mol[chunk_snap_idx]),
                    ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                    ref(lig.info), ref(CachedData),
                    lig.m_mol[0].size() * i / num_threads,
                    lig.m_mol[0].size() * (i + 1) / num_threads, ref(Energy[i]),
                    ref(Force[i]), chunk_snap_idx, ref(Arguments),
                    ref(pro_in.parameters), ref(lig.parameters));
            }
            for (auto &t : threads)
              if (t.joinable())
                if (t.joinable())
                  t.join();
          }
          time_cal += std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::steady_clock::now() - start_energy_eval);
          // cout << "done cal snap " << chunk_snap_idx << " ";
        } // Close else block for standard energy eval

        if (Arguments.self_int == true) {
          remove_double_counting(CachedData, Energy);
        }

        for (size_t i{}; i < Energy.size(); ++i) {
          E_E_ += Energy[i].elec_ener;
          E_V_ += Energy[i].vdw_ener;
          if (Arguments.flag_cal_force == true) {
            f_e_ += Force[i].forces[0];
            f_v_ += Force[i].forces[1];
          }
        }

        process_adding_energies_forces(
            Arguments, input_name, E_v, E_e, E_E_, E_V_, f_e_, f_v_, E_v_vector,
            E_e_vector, Force_e, Force_v, Force_e_vector, Force_v_vector,
            size_snap, global_snap_idx);

        // Keep time_of_snap in sync with E_e for streaming/XTC mode
        while (time_of_snap.size() < E_e.size()) {
          double t = pro_in.m_mol[chunk_snap_idx].time_of_snap;
          if (Arguments.tu == "ns")
            t /= 1000.0;
          time_of_snap.push_back(t);
        }

        auto start_sc{chrono::steady_clock::now()};

        size_t num_threads_rest(
            pro_in.Index.size() - 1 > 0 ? pro_in.Index.size() - 1 : 1);
        size_t ite_thread{};
        while (num_threads_rest > 0) {
          if (num_threads_rest <= num_threads) {
            for (size_t i{}; i < num_threads_rest; ++i) {
              threads[i] =
                  thread(&MOL::take_com, &pro_in.m_mol[chunk_snap_idx],
                         ref(pro_in.parameters), i + ite_thread,
                         pro_in.Index[i + ite_thread + 1].get_res_pos_e(0),
                         ref(pro_in.Index), ref(pro_in.info));
            }
            for_each(threads.begin(), threads.begin() + num_threads_rest,
                     mem_fn(&thread::join));
            num_threads_rest = 0ull;
          } else {
            for (size_t i{}; i < num_threads; ++i) {
              threads[i] =
                  thread(&MOL::take_com, &pro_in.m_mol[chunk_snap_idx],
                         ref(pro_in.parameters), i + ite_thread,
                         pro_in.Index[i + ite_thread + 1].get_res_pos_e(0),
                         ref(pro_in.Index), ref(pro_in.info));
            }
            for (auto &t : threads)
              if (t.joinable())
                if (t.joinable())
                  t.join();
            num_threads_rest -= num_threads;
            ite_thread += num_threads;
          }
        }

        num_threads_rest = lig.Index.size() - 1 > 0 ? lig.Index.size() - 1 : 1;
        ite_thread = 0;
        while (num_threads_rest > 0) {
          if (num_threads_rest <= num_threads) {
            for (size_t i{}; i < num_threads_rest; ++i) {
              threads[i] =
                  thread(&MOL::take_com, &lig.m_mol[chunk_snap_idx],
                         ref(lig.parameters), i + ite_thread,
                         lig.Index[i + ite_thread + 1].get_res_pos_e(0),
                         ref(lig.Index), ref(lig.info));
            }
            for_each(threads.begin(), threads.begin() + num_threads_rest,
                     mem_fn(&thread::join));
            num_threads_rest = 0ull;
          } else {
            for (size_t i{}; i < num_threads; ++i) {
              threads[i] =
                  thread(&MOL::take_com, &lig.m_mol[chunk_snap_idx],
                         ref(lig.parameters), i + ite_thread,
                         lig.Index[i + ite_thread + 1].get_res_pos_e(0),
                         ref(lig.Index), ref(lig.info));
            }
            for (auto &t : threads)
              if (t.joinable())
                if (t.joinable())
                  t.join();
            num_threads_rest -= num_threads;
            ite_thread += num_threads;
          }
        }

        impl_add_coms_2_mols(lig, pro_in, Arguments, chunk_snap_idx);
        /*****************************************************************************************************************/

        for (size_t i{}; i < num_threads; ++i) {
          threads[i] = thread(
              take_min, ref(pro_in.m_mol[chunk_snap_idx].Com),
              ref(lig.m_mol[chunk_snap_idx].Com), ref(dis_min[i]),
              i * pro_in.m_mol[chunk_snap_idx].Com.size() / num_threads,
              (i + 1) * pro_in.m_mol[chunk_snap_idx].Com.size() / num_threads,
              0, lig.m_mol[chunk_snap_idx].Com.size(), ref(Arguments), 0,
              d_standard);
        }
        for (auto &t : threads)
          if (t.joinable())
            if (t.joinable())
              t.join();

        time_sc += chrono::duration_cast<chrono::milliseconds>(
            chrono::steady_clock::now() - start_sc);
        // cout << "\nsnap " << chunk_snap_idx << " mol1 consists of " <<
        // pro_in.m_mol[chunk_snap_idx].Com.size() << " mol2 consists of " <<
        // lig.m_mol[chunk_snap_idx].Com.size() << " residues\n"; cout <<
        // "\nsnap:
        // "<<chunk_snap_idx<< ": numof dismin " << dis_min.size() << " "; for
        // (auto& d : dis_min) 	cout << d.base.size() << ' '; cout <<
        // "mem size: " << Arguments.nt * (sizeof(std::vector<dmin_base>) +
        // (sizeof(dmin_base) * dis_min[0].base.size())); cout << '\n';
        print_sidechain_contact(num_threads, dis_min, Arguments, out_sc,
                                sidechain_contact_num, sc_rec_1traj,
                                sc_lg_1traj, input_name - Arguments.b_input);

        // --- IFP per-frame processing ---
        if (Arguments.ifp.flag_ifp && ifp_analyzer) {
          auto start_ifp = std::chrono::steady_clock::now();
          size_t lg_snap_ifp = (lig.m_mol.size() > 1) ? chunk_snap_idx : 0;
          ifp_analyzer->process_frame(
              pro_in.m_mol[chunk_snap_idx],
              lig.m_mol[lg_snap_ifp],
              pro_in.info, lig.info,
              pro_in.Index, lig.Index,
              pro_in.info.size(),
              lig.info.size(),
              dis_min,
              Arguments.ifp
          );
          time_ifp += std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::steady_clock::now() - start_ifp);
        }

        if (Arguments.cal_dis_sc_contact == true) {
          auto start_cscd = std::chrono::steady_clock::now();
          impl_cscd(pro_in, lig, Arguments, log, threads, chunk_snap_idx,
                    sc_flag, input_name, cscd_out_name, cscd_out_dis, paths,
                    out_cal_dis_files);
          time_cscd += std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::steady_clock::now() - start_cscd);
        }

        if (global_has_nucleic_acid) {
          auto start_na = std::chrono::steady_clock::now();
          na_profiler.process_frame(
              pro_in.m_mol[chunk_snap_idx], lig.m_mol[chunk_snap_idx],
              pro_in.info, lig.info, pro_in.Index, lig.Index, chunk_snap_idx,
              pro_in.m_mol[chunk_snap_idx].time_of_snap);
          time_na += std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::steady_clock::now() - start_na);
        }

        if (Arguments.flag_pqr == 1) {
          write_pqr(pro_in, lig, input_line, input_name, chunk_snap_idx,
                    Arguments, paths);
        }
        write_pdb(pro_in, lig, Arguments, chunk_snap_idx, paths, input_line,
                  input_name);

        if (Arguments.flag_cas) {
          auto start_cas = std::chrono::steady_clock::now();

          cas_analyzer.process_frame(
              pro_in, lig, chunk_snap_idx, num_threads, dis_min, CachedData,
              Arguments, info_output, rec_free_radii, lg_free_radii,
              rec_complex_radii, lg_complex_radii, frame_gbsa_total);

          time_cas += std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::steady_clock::now() - start_cas);
        }

        // Network Analysis (DCCM) Data Collection
        auto start_network = std::chrono::steady_clock::now();
        if (Arguments.flag_network || Arguments.flag_pca ||
            Arguments.lspmd.flag_lspmd) {
          dccm_analyzer.process_frame(pro_in, lig, chunk_snap_idx);
        }
        time_network += std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_network);

        // --- TICA per-frame accumulation (pairwise Cα distances + Kabsch
        // align) ---
        auto start_tica = std::chrono::steady_clock::now();
        if (Arguments.msm.flag_tica && tica_analyzer && tica_n_ca > 0) {
          size_t nf = dccm_analyzer.get_total_frames();
          if (nf > 0) {
            size_t last = nf - 1;
            const auto &tx = dccm_analyzer.get_traj_x()[last];
            const auto &ty = dccm_analyzer.get_traj_y()[last];
            const auto &tz = dccm_analyzer.get_traj_z()[last];

            // Fill tica_mobile (pre-allocated N_ca×3)
            for (int a = 0; a < tica_n_ca; ++a) {
              tica_mobile(a, 0) = tx[static_cast<size_t>(a)];
              tica_mobile(a, 1) = ty[static_cast<size_t>(a)];
              tica_mobile(a, 2) = tz[static_cast<size_t>(a)];
            }

            if (!tica_ref_initialized) {
              // Frame 0: center and store as Kabsch reference
              tica_ref_centered = tica_mobile;
              tica_ref_centered.rowwise() -=
                  tica_ref_centered.colwise().mean().eval();
              tica_mobile = tica_ref_centered;
              tica_ref_initialized = true;
            } else {
              // Frame > 0: zero-alloc Kabsch alignment
              TitanAlignment::align_frame_inplace(tica_ref_centered,
                                                  tica_mobile, tica_R);
            }

            // Science fix: pairwise inter-Cα distances (rotation+translation
            // invariant) D = N_ca*(N_ca-1)/2 — no spurious global
            // rotation/translation modes
            int feat_idx = 0;
            for (int a = 0; a < tica_n_ca; ++a) {
              for (int b = a + 1; b < tica_n_ca; ++b) {
                double dx = tica_mobile(a, 0) - tica_mobile(b, 0);
                double dy = tica_mobile(a, 1) - tica_mobile(b, 1);
                double dz = tica_mobile(a, 2) - tica_mobile(b, 2);
                tica_features[feat_idx++] =
                    std::sqrt(dx * dx + dy * dy + dz * dz);
              }
            }

            tica_analyzer->accumulate_frame(tica_features);

            // --- MSM Phase 2A: collect TICA projection in RAM ---
            if (Arguments.msm.flag_msm && msm_analyzer &&
                tica_analyzer->is_solved()) {
              Eigen::VectorXd proj =
                  tica_analyzer->project_frame(tica_features);
              int current_replica = input_name - Arguments.b_input;
              msm_analyzer->add_projected_frame(proj, current_replica);
              msm_analyzer->record_frame_time(
                  pro_in.m_mol[chunk_snap_idx].time_of_snap, current_replica);
            }
          }
        }
        time_tica += std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_tica);

        auto start_tda = std::chrono::steady_clock::now();
        if (Arguments.flag_tda) {
          if (!Arguments.flag_sidechain) {
            TitanUI::print_info(
                "TDA is running without sidechain filtering (extracting ALL "
                "heavy atoms). Enable flag -sc for active-only.");
            std::vector<int>
                target_idx_1; // Empty means "Extract ALL Heavy Atoms"
            std::vector<int> target_idx_2;
            std::string tda_out = "tda_global_pointcloud.csv";
            tda_analyzer.export_TDA_pointcloud(
                pro_in, lig, input_name, chunk_snap_idx, target_idx_1,
                target_idx_2, Arguments, tda_out);
          } else {
            std::unordered_set<int> active_1;
            std::unordered_set<int> active_2;

            bool mol1_lig = (Arguments.mol_1_type == "ligand" ||
                             Arguments.mol_1_type == "lig");
            bool mol2_lig = (Arguments.mol_2_type == "ligand" ||
                             Arguments.mol_2_type == "lig");

            // Directly grab the active residue/atom numbers forming contacts
            // from d_min array
            for (size_t t{}; t < num_threads; ++t) {
              for (const auto &contact : dis_min[t].base) {
                if (mol1_lig) {
                  active_1.insert(contact.atom_num_rec -
                                  1); // 0-indexed atom num
                } else {
                  active_1.insert(contact.resnr_rec);
                }

                if (mol2_lig) {
                  active_2.insert(contact.atom_num - 1); // 0-indexed atom num
                } else {
                  active_2.insert(contact.resnr_ab);
                }
              }
            }

            std::vector<int> target_idx_1(active_1.begin(), active_1.end());
            std::vector<int> target_idx_2(active_2.begin(), active_2.end());

            std::string tda_out = "tda_global_pointcloud.csv";
            tda_analyzer.export_TDA_pointcloud(
                pro_in, lig, input_name, chunk_snap_idx, target_idx_1,
                target_idx_2, Arguments, tda_out);
          }
        }
        time_tda += std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_tda);

        auto start_emap = std::chrono::steady_clock::now();
        if (Arguments.flag_energy_same_force == false &&
            (Arguments.energy_map == 1 || Arguments.energy_map == 3 ||
             Arguments.energy_map == 4) &&
            !CachedData.elec_ener.empty()) {
          wrap_cal_energy_map(threads, pro_in, lig, Arguments, chunk_snap_idx,
                              dis_min, CachedData, Energy_res_mol1,
                              Energy_res_mol2, info_output, infoener_res);
        }
        time_emap += std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_emap);
        // cout << "check emap..";
        // cout << "Ener resmol1 size: " << Energy_res_mol1.size() << " size 1
        // ele: " << sizeof(Energy_res_mol1[0]) << '\n'; cout << "Ener resmol2
        // size: " << Energy_res_mol2.size() << " size 1 ele: " <<
        // sizeof(Energy_res_mol2[0]) << '\n';
        auto start_fmap{std::chrono::steady_clock::now()};
        if (Arguments.flag_cal_force == true &&
            (Arguments.force_map == 1 || Arguments.force_map == 3) &&
            !CachedData.elec_ener.empty()) {
          wrap_cal_force_map(threads, pro_in, lig, Arguments, chunk_snap_idx,
                             dis_min, CachedData, Force_res_mol1,
                             Force_res_mol2, info_output);
        }
        time_fmap += std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_fmap);

        auto start_hb{std::chrono::steady_clock::now()};

        if (Arguments.flag_h_bond == 1) {
          find_h_bonds(pro_in.m_mol[chunk_snap_idx], pro_in.m_mol[0],
                       lig.m_mol[chunk_snap_idx], lig.m_mol[0], groups_1,
                       groups_2, Arguments, chunk_snap_idx, pro_in.info,
                       lig.info, Hbonds_vector);

          if (Arguments.nt >= 2) {
            thread t_hb_1(add_hb_res, ref(pro_in.m_mol[chunk_snap_idx]),
                          ref(lig.m_mol[chunk_snap_idx]), ref(pro_in.info),
                          ref(lig.info), Arguments, ref(hb_rec_1traj));
            thread t_hb_2(add_hb_res, ref(lig.m_mol[chunk_snap_idx]),
                          ref(pro_in.m_mol[chunk_snap_idx]), ref(lig.info),
                          ref(pro_in.info), Arguments, ref(hb_lg_1traj));
            if (t_hb_1.joinable())
              if (t_hb_1.joinable())
                t_hb_1.join();
            if (t_hb_2.joinable())
              if (t_hb_2.joinable())
                t_hb_2.join();
          } else {
            add_hb_res(pro_in.m_mol[chunk_snap_idx], lig.m_mol[chunk_snap_idx],
                       pro_in.info, lig.info, Arguments, hb_rec_1traj);
            add_hb_res(lig.m_mol[chunk_snap_idx], pro_in.m_mol[chunk_snap_idx],
                       lig.info, pro_in.info, Arguments, hb_lg_1traj);
          }
        } // End of flag_h_bond block

        time_hb += std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_hb);

        if (Arguments.n_wat >= 0 && !water_traj.m_mol.empty()) {
          if (chunk_snap_idx >= water_traj.m_mol.size()) {
            if (chunk_snap_idx == 1) {
              std::cout
                  << "\n[WARNING] Water group has only "
                  << water_traj.m_mol.size() << " frame(s), but trajectory has "
                  << NUM << " frames.\n"
                  << "Water analysis will only process available frames!\n";
            }
          } else {
            // Optimization: Distill heavy atom indices just once
            static std::vector<size_t> rec_heavy_atoms;
            static std::vector<size_t> lig_heavy_atoms;

            if (chunk_snap_idx == 0) {
              rec_heavy_atoms.clear();
              lig_heavy_atoms.clear();

              for (size_t i = 0; i < pro_in.info.size(); ++i) {
                if (pro_in.info[i].get_atom()[0] != 'H')
                  rec_heavy_atoms.emplace_back(i);
              }
              for (size_t i = 0; i < lig.info.size(); ++i) {
                if (lig.info[i].get_atom()[0] != 'H')
                  lig_heavy_atoms.emplace_back(i);
              }
            }

            const MOL &w_frame = water_traj.m_mol[chunk_snap_idx];
            const MOL &p_frame = pro_in.m_mol[chunk_snap_idx];
            size_t l_snap = (lig.m_mol.size() > 1) ? chunk_snap_idx : 0;
            const MOL &l_frame = lig.m_mol[l_snap];

            /*if (chunk_snap_idx == 0) {
              std::cout << "[DEBUG WB] Snapshot 0 - rec heavy atoms: "
                        << rec_heavy_atoms.size()
                        << ", lig heavy atoms: " << lig_heavy_atoms.size()
                        << ", water atoms: " << w_frame.pos.x.size() << "\n";
            }*/

            // GÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â»ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Âi
            // class WaterAnalyzer
            // ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¾ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¹Ãƒâ€¦Ã¢â‚¬Å“ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â£
            // ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¾ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¹Ãƒâ€¦Ã¢â‚¬Å“ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€šÃ‚Â ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â°ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â»ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â£c
            // tÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ch
            // file
            // cÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â»ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±c
            // kÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â»ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â³
            // gÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â»ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Ân
            // gÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â ng!
            auto start_water = std::chrono::steady_clock::now();
            water_analyzer.process_frame(
                p_frame, l_frame, w_frame, rec_heavy_atoms, lig_heavy_atoms,
                pro_in.info, lig.info, input_name, chunk_snap_idx,
                info_output.get_global_water_bridges_timeseries(),
                info_output.get_water_pdb_filename(Arguments, input_line,
                                                   input_name));
            time_water += std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start_water);

            // --- GIST per-frame ---
            if (Arguments.gist.flag_gist && gist_analyzer) {
              const bool gist_use_tica =
                  Arguments.msm.flag_tica && tica_ref_initialized;
              gist_analyzer->process_frame(
                  p_frame, l_frame, w_frame, pro_in.info,
                  gist_use_tica ? tica_ref_centered : Eigen::MatrixXd{},
                  gist_use_tica ? tica_R : Eigen::Matrix3d::Identity(),
                  (global_snap_idx == 0), gist_use_tica);
            }
          }
        }
        // cout << "check hbond..";
        clock_t start_new_emap{clock()};
        if (Arguments.flag_new_emap == true && !CachedData.elec_ener.empty()) {
          imp_cal_energy_rec_res(rec_ener_inter, pro_in, lig, chunk_snap_idx, 0,
                                 0, dis_min, CachedData, Arguments);
        }
        time_new_emap += clock() - start_new_emap;
        // cout << "check new emap..";
        auto start_subgroups = std::chrono::steady_clock::now();
        if (Arguments.flag_sub_groups == true && Arguments.flag_gbsa == false &&
            !CachedData.elec_ener.empty()) {
          Energy_force_subgroups(Groups_from_index_file, pro_in, lig,
                                 CachedData, EF_subgroups, Energy_groups_out,
                                 chunk_snap_idx, NUM, Arguments, threads);
          if (Arguments.flag_cal_mindis == true)
            mindis_subgroups(Groups_from_index_file, pro_in, lig,
                             chunk_snap_idx, NUM, Arguments, threads, MIN_DIS);
        }
        time_subgroups_cal +=
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start_subgroups);
        // cout << "check sub groups..";
        auto start_dis2commol = std::chrono::steady_clock::now();
        if (Arguments.com_dis_energy_2mol != 4)
          dis_com_2mol(pro_in, lig, chunk_snap_idx, NUM, Arguments, input_name,
                       paths, E_E_, E_V_, info_output);
        if (Arguments.com_dis_2mol != 0) {
          fstream output_dis_2mol_file(
              paths.get_path(input_name - Arguments.b_input) + "/" +
                  info_output.getDisCom2Mol(),
              ios::out | ios::app);
          pro_in.impl_com_dis_2mol(lig, Arguments, info_output,
                                   output_dis_2mol_file, chunk_snap_idx);
        }
        // cout << "begin check contact ";
        if (Arguments.flag_cal_contact == true) {
          calc_contact_map(pro_in, lig, num_threads, dis_min, Arguments,
                           chunk_snap_idx, cont_matrix, pos_lookup_1,
                           pos_lookup_2);
        }
        // cout << "checked contact..";
        time_other += std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_dis2commol);

        // --- ML Feature Extraction ---
        if (Arguments.flag_ml) {
          ml_extractor.process_frame(
              pro_in.m_mol[chunk_snap_idx], lig.m_mol[chunk_snap_idx],
              pro_in.info, lig.info, lig.parameters, dis_min, E_E_ + E_V_, 6.0);
        }

        // In chunk streaming, we do not clear individual MOLs here, we clear
        // the entire vector before reading the next chunk.

      } // End of inner chunk loop

      // We update the sizes for the next iteration
      // INITIATE THREAD-LOCAL FORCES FOR OMP LOOP
      for (size_t f = 0; f < pro_in.m_mol.size(); ++f) {
        if (pro_in.compute_forces)
          pro_in.m_mol[f].initialize_forces(num_threads);
      }
      for (size_t f = 0; f < lig.m_mol.size(); ++f) {
        if (lig.compute_forces)
          lig.m_mol[f].initialize_forces(num_threads);
      }

      global_frame_counter += chunk_num_snaps;

    } // End of outer chunk while loop

    tot_snap_allfiles += global_frame_counter;

    // --- IFP post-trajectory export ---
    // (aggregate across all files at the very end, handled after all-files loop)

    // Reset TICA lag buffer at trajectory boundary
    if (Arguments.msm.flag_tica && tica_analyzer) {
      tica_analyzer->reset_lag_buffer();
    }
    // Insert MSM trajectory break at each file boundary
    if (Arguments.msm.flag_msm && msm_analyzer &&
        Arguments.msm.msm_auto_boundary) {
      msm_analyzer->add_trajectory_break();
      TitanUI::log_kv(log, "MSM trajectory break at frame",
                      std::to_string(msm_analyzer->n_frames()));
    }

    if (Arguments.n_wat >= 0 && !water_traj.m_mol.empty()) {
      water_analyzer.finalize_trajectory();
    }

    if (!out_cal_dis_files.empty()) {
      for (auto &file : out_cal_dis_files)
        if (file.is_open())
          file.close();
    }

    auto start_print_ener = std::chrono::steady_clock::now();
    if (Arguments.print_energy_res == true) {
      if (Arguments.flag_energy_same_force == true) {
        print_energy_res(Force_res_mol1, info_output.get_group_1_name(), paths,
                         input_name - Arguments.b_input, Arguments);
        print_energy_res(Force_res_mol2, info_output.get_group_2_name(), paths,
                         input_name - Arguments.b_input, Arguments);
      } else {
        print_energy_res(Energy_res_mol1, info_output.get_group_1_name(), paths,
                         input_name - Arguments.b_input, Arguments);
        print_energy_res(Energy_res_mol2, info_output.get_group_2_name(), paths,
                         input_name - Arguments.b_input, Arguments);
      }

      Energy_res_stat_1.emplace_back(Energy_res_mol1);
      Energy_res_stat_2.emplace_back(Energy_res_mol2);
    }
    //<< std::endl;
    time_print_ener_res +=
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_print_ener);

    auto start_print_force = std::chrono::steady_clock::now();
    if (Arguments.flag_cal_force == true && Arguments.print_force_res == true) {
      print_force_res(Force_res_mol1, info_output.get_group_1_name(), paths,
                      input_name - Arguments.b_input, Arguments);
      print_force_res(Force_res_mol2, info_output.get_group_2_name(), paths,
                      input_name - Arguments.b_input, Arguments);
    }
    time_print_force += std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_print_force);

    auto start_maxforce = std::chrono::steady_clock::now();
    if (Arguments.flag_maxforce == true) {
      find_maxforce_res(Force_res_mol1, info_output.get_group_1_name(), paths,
                        input_name - Arguments.b_input, Arguments,
                        Force_res_stat_mol1);
      find_maxforce_res(Force_res_mol2, info_output.get_group_2_name(), paths,
                        input_name - Arguments.b_input, Arguments,
                        Force_res_stat_mol2);
    }
    time_maxforce += std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_maxforce);

    auto start_new_emap_other = std::chrono::steady_clock::now();
    if (Arguments.flag_new_emap == true) {
      rec_ener_inter.ana_emap(Arguments, info_output, paths, input_name);
    }
    time_new_emap_other +=
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_new_emap_other);

    auto start_avg_emap = std::chrono::steady_clock::now();
    if (Arguments.energy_map == 1 && tot_snap_allfiles > 0) {
      infoener_res.print_average_emap(pro_in, lig, tot_snap_allfiles,
                                      info_output, Arguments);
    }
    time_avg_emap += std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_avg_emap);

    if (Arguments.flag_h_bond == 1) {
      hb_lg_1traj.snap_num = global_frame_counter;
      hb_rec_1traj.snap_num = global_frame_counter;
      res_hb_lg.emplace_back(hb_lg_1traj);
      res_hb_rec.emplace_back(hb_rec_1traj);
    }

    auto start_contact = std::chrono::steady_clock::now();
    if (Arguments.flag_cal_contact == true) {
      sc_lg_1traj.snap_num = global_frame_counter;
      sc_rec_1traj.snap_num = global_frame_counter;

      sc_rec_1traj.clear_map();
      sc_lg_1traj.clear_map();

      res_sc_lg.emplace_back(sc_lg_1traj);
      res_sc_rec.emplace_back(sc_rec_1traj);

      string cont_map_filename{info_output.getContmap() +
                               to_string(input_name) + ".dat"};
      fstream cont_map_file(cont_map_filename, ios::out);

      print_cont_map(cont_matrix, cont_map_file, global_frame_counter);
    }
    time_contact += std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_contact);

    /**************************************************************************/
    TitanUI::log_section(log, "PERFORMANCE & TIMING REPORT");
    double total_traj_time = std::chrono::duration_cast<milliseconds>(
                                 std::chrono::steady_clock::now() - start)
                                 .count() /
                             1000.0;

    TitanUI::log_time(log, "Parsing PDB", time_pdb.count() / 1000.0);
    TitanUI::log_time(log, "Solving PBC", time_solve_pbc.count() / 1000.0);
    TitanUI::log_time(log, "GBSA Calculation", time_gbsa.count() / 1000.0);
    TitanUI::log_time(log, "CAS Calculation", time_cas.count() / 1000.0);
    TitanUI::log_time(log, "Network Data Collection",
                      time_network.count() / 1000.0);
    TitanUI::log_time(log, "Standard Elec+vdW Calc", time_cal.count() / 1000.0);
    TitanUI::log_time(log, "Sidechain Contact", time_sc.count() / 1000.0);
    TitanUI::log_time(log, "Sidechain CSCD", time_cscd.count() / 1000.0);
    TitanUI::log_time(log, "Hydrogen Bonds", time_hb.count() / 1000.0);
    TitanUI::log_time(log, "Energy Map (EMAP)", time_emap.count() / 1000.0);
    TitanUI::log_time(log, "Force Map (FMAP)", time_fmap.count() / 1000.0);
    TitanUI::log_time(log, "Contact Map & Dis2Mol",
                      time_other.count() / 1000.0);
    TitanUI::log_time(log, "Water Bridges Analysis",
                      time_water.count() / 1000.0);
    TitanUI::log_time(log, "Subgroups Calc",
                      time_subgroups_cal.count() / 1000.0);
    TitanUI::log_time(log, "Nucleic Acid Profiler", time_na.count() / 1000.0);
    TitanUI::log_time(log, "LLPS Analysis", time_llps.count() / 1000.0);
    TitanUI::log_time(log, "TICA/MSM Accumulation", time_tica.count() / 1000.0);
    TitanUI::log_time(log, "IFP Analysis", time_ifp.count() / 1000.0);

    log << "-----------------------------------------------------------------"
           "--"
           "-------------\n";
    TitanUI::log_time(log, "TOTAL TRAJECTORY TIME", total_traj_time);
    std::cout << std::endl;
    if (Arguments.flag_cal_contact == true) {
      if (average_contact_matrix.empty()) {
        average_contact_matrix = vector<vector<map_pos>>(cont_matrix.size());
        for (size_t a_i{}; a_i < average_contact_matrix.size(); ++a_i) {
          average_contact_matrix[a_i] =
              vector<map_pos>(cont_matrix[a_i].size(), map_pos());
          for (size_t a_j{}; a_j < average_contact_matrix[a_i].size(); ++a_j)
            average_contact_matrix[a_i][a_j] = cont_matrix[a_i][a_j];
        }
      } else {
        for (size_t a_i{}; a_i < average_contact_matrix.size(); ++a_i) {
          for (size_t a_j{}; a_j < average_contact_matrix[a_i].size(); ++a_j) {
            if (a_i < cont_matrix.size() && a_j < cont_matrix[a_i].size())
              average_contact_matrix[a_i][a_j].cont_num +=
                  cont_matrix[a_i][a_j].cont_num;
          }
        }
      }
    }

    //<< std::endl;
    //<< std::endl;
    if (check_end_cal(Arguments, input_line, input_name) == true) {
      //<< Arguments.energy_map << std::endl;
      if (Arguments.energy_map == 1 || Arguments.energy_map == 4) {
        auto start_enermaptype = chrono::steady_clock::now();
        TitanUI::log_kv(log, "Status", "Building energy maptype...");
        //<< std::endl;
        // Restore Frame 0 data so print_energy_pdbtype can write coords
        if (pro_in.m_mol.empty())
          pro_in.m_mol = saved_pro_frame0;
        if (lig.m_mol.empty())
          lig.m_mol = saved_lig_frame0;

        infoener_res.wrapper_impl_ener_maptype(pro_in, lig, info_output,
                                               Arguments, tot_snap_allfiles);
        TitanUI::log_kv(log, "Status", "Energy maptype done.");
        time_energymaptype += chrono::duration_cast<milliseconds>(
            chrono::steady_clock::now() - start_enermaptype);
      }

      break;
    }

    std::cout << std::endl;
  }

  TitanUI::print_info(
      "\033[3mFinalizing module calculations and aggregating data...\033[0m");

  // --- IFP: Aggregate occupancy and export CSV ---
  if (Arguments.ifp.flag_ifp && ifp_analyzer) {
    ifp_analyzer->aggregate_occupancy();
    std::string ifp_outfile = Arguments.ifp.output_filename.empty()
                                  ? info_output.get_ifp_filename()
                                  : Arguments.ifp.output_filename;
    ifp_analyzer->export_csv(ifp_outfile);
  }

  out_sc.close();
  if (status == false && input_name == Arguments.b_input) {
    TitanUI::print_error("Error occurred during processing!");
    return 1;
  }

  TitanUI::print_module_header("Final Results Analysis & Summary");

  try {
    out(E_e, E_v, time_of_snap, cal_out, Arguments, input_name, E_e_vector,
        E_v_vector, size_snap);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] out(): " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] out() crashed!\n";
  }

  try {
    paths.out_to_folders(E_e, E_v, time_of_snap, E_e_vector, E_v_vector,
                         Arguments);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] paths.out_to_folders(): " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] paths.out_to_folders() crashed!\n";
  }

  if (Arguments.flag_cal_force == true) {
    try {
      out_force(Force_e, Force_v, time_of_snap, Arguments, input_name,
                Force_e_vector, Force_v_vector, size_snap);
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] out_force(): " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] out_force() crashed!\n";
    }

    try {
      paths.out_force_to_folders(Force_e, Force_v, time_of_snap, Force_e_vector,
                                 Force_v_vector, Arguments);
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] paths.out_force_to_folders(): " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] paths.out_force_to_folders() crashed!\n";
    }
  }

  if (Arguments.flag_maxforce == true) {
    try {
      stat_maxforce_res(Force_res_stat_mol1, info_output.get_group_1_name(),
                        info_output, Arguments);
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] stat_maxforce_res() mol1: " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] stat_maxforce_res() mol1 crashed!\n";
    }

    try {
      stat_maxforce_res(Force_res_stat_mol2, info_output.get_group_2_name(),
                        info_output, Arguments);
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] stat_maxforce_res() mol2: " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] stat_maxforce_res() mol2 crashed!\n";
    }
  }

  try {
    cal_sidechain_contact_time(Arguments, time_of_snap, sidechain_contact_num,
                               info_output, paths);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] cal_sidechain_contact_time(): " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] cal_sidechain_contact_time() crashed!\n";
  }

  try {
    impl_entropy_duan_zhang(paths, E_e_vector, E_v_vector, time_of_snap,
                            Arguments, E_e, E_v);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] impl_entropy_duan_zhang(): " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] impl_entropy_duan_zhang() crashed!\n";
  }

  try {
    print_cscd(Arguments, cscd_out_dis, cscd_out_name);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] print_cscd(): " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] print_cscd() crashed!\n";
  }

  try {
    paths.out_volume(time_of_snap, Volumes, Arguments);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] paths.out_volume(): " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] paths.out_volume() crashed!\n";
  }

  if (Arguments.flag_sub_groups == true) {
    try {
      post_process_Energy_force_subgroups(time_of_snap, EF_subgroups, paths,
                                          Arguments, info_output);
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] post_process_Energy_force_subgroups(): " << e.what()
                << "\n";
    } catch (...) {
      std::cerr << "[ERROR] post_process_Energy_force_subgroups() crashed!\n";
    }

    if (Arguments.flag_cal_mindis == true) {
      try {
        post_proc_mindis(Groups_from_index_file, Arguments, info_output,
                         MIN_DIS);
      } catch (const std::exception &e) {
        std::cerr << "[ERROR] post_proc_mindis(): " << e.what() << "\n";
      } catch (...) {
        std::cerr << "[ERROR] post_proc_mindis() crashed!\n";
      }
    }
  }

  try {
    stat_hbonds(Hbonds_vector, info_output, Arguments, time_of_snap,
                input_line);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] stat_hbonds(): " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] stat_hbonds() crashed!\n";
  }

  if (Arguments.flag_h_bond == 1) {
    if (num_threads > 1) {
      thread s_hbres_1(stat_hb_res, ref(res_hb_rec),
                       info_output.get_res_list_hb_rec_filename());
      thread s_hbres_2(stat_hb_res, ref(res_hb_lg),
                       info_output.get_res_list_hb_lg_filename());
      if (s_hbres_1.joinable())
        if (s_hbres_1.joinable())
          s_hbres_1.join();
      if (s_hbres_2.joinable())
        if (s_hbres_2.joinable())
          s_hbres_2.join();
    } else {
      try {
        stat_hb_res(res_hb_rec, info_output.get_res_list_hb_rec_filename());
      } catch (const std::exception &e) {
        std::cerr << "[ERROR] stat_hb_res() rec: " << e.what() << "\n";
      } catch (...) {
        std::cerr << "[ERROR] stat_hb_res() rec crashed!\n";
      }

      try {
        stat_hb_res(res_hb_lg, info_output.get_res_list_hb_lg_filename());
      } catch (const std::exception &e) {
        std::cerr << "[ERROR] stat_hb_res() lg: " << e.what() << "\n";
      } catch (...) {
        std::cerr << "[ERROR] stat_hb_res() lg crashed!\n";
      }
    }
  }

  if (Arguments.flag_cal_contact == true) {
    try {
      stat_sidechain_contact(
          res_sc_rec, info_output.get_res_list_sc_contact_rec_filename());
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] stat_sidechain_contact() rec: " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] stat_sidechain_contact() rec crashed!\n";
    }

    try {
      stat_sidechain_contact(res_sc_lg,
                             info_output.get_res_list_sc_contact_lg_filename());
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] stat_sidechain_contact() lg: " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] stat_sidechain_contact() lg crashed!\n";
    }
  }

  if (Arguments.print_energy_res == true) {
    try {
      stat_energy_res(Energy_res_stat_1, 1, Arguments, info_output,
                      info_output.get_group_1_name());
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] stat_energy_res() 1: " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] stat_energy_res() 1 crashed!\n";
    }

    try {
      stat_energy_res(Energy_res_stat_2, 1, Arguments, info_output,
                      info_output.get_group_2_name());
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] stat_energy_res() 2: " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] stat_energy_res() 2 crashed!\n";
    }
  }

  if (Arguments.flag_cal_contact == true) {
    try {
      fstream ave_cont_file(info_output.getOuputAverageContactMap(), ios::out);
      print_cont_map(average_contact_matrix, ave_cont_file, tot_snap_allfiles);
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] print_cont_map(): " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] print_cont_map() crashed!\n";
    }
  } else {
    log << "Note: Average contact map output is skipped because the option "
           "is not enabled in input.\n";
  }

  // --- Post-loop analysis (extracted to engine_postproc.cpp) ---------------
  auto start_postproc = std::chrono::steady_clock::now();

  auto start_pp_gbsa = std::chrono::steady_clock::now();
  run_gbsa_postproc(Arguments, tot_snap_allfiles, ml_extractor, global_gbsa_agg,
                    cas_analyzer, na_profiler, global_has_nucleic_acid,
                    info_output, log);
  auto time_pp_gbsa = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_pp_gbsa);

  auto start_pp_tica = std::chrono::steady_clock::now();
  run_tica_msm_postproc(Arguments, tica_analyzer, msm_analyzer, num_threads,
                        info_output, log);
  auto time_pp_tica = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_pp_tica);

  run_dham_postproc(Arguments, msm_analyzer, num_threads, info_output, log);

  // Drain LLPS worker queue before finalize
  auto start_pp_llps = std::chrono::steady_clock::now();
  if (Arguments.llps.flag_llps && llps_analyzer)
    llps_analyzer->flush_async();
  run_llps_postproc(Arguments, llps_analyzer, num_threads, info_output, log);
  auto time_pp_llps = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_pp_llps);

  auto start_pp_dccm = std::chrono::steady_clock::now();
  run_dccm_pca_postproc(Arguments, tot_snap_allfiles, dccm_analyzer,
                        num_threads, time_network, info_output, log);
  auto time_pp_dccm = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_pp_dccm);

  auto start_pp_lspmd = std::chrono::steady_clock::now();
  run_lspmd_postproc(Arguments, dccm_analyzer, num_threads, info_output, log);
  auto time_pp_lspmd = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_pp_lspmd);

  auto start_pp_water = std::chrono::steady_clock::now();
  if (Arguments.n_wat >= 0) {
    TitanUI::print_module_header(
        "Multi-Trajectory Global Water Bridges Analysis");
    water_analyzer.export_global_summary(
        info_output.get_global_water_bridges_summary(),
        info_output.get_global_water_survival());
  }
  auto time_pp_water = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_pp_water);

  auto start_pp_gist = std::chrono::steady_clock::now();
  if (Arguments.gist.flag_gist && gist_analyzer) {
    TitanUI::print_module_header("GIST Water Thermodynamics");
    gist_analyzer->export_dx_maps(info_output.get_gist_prefix(),
                                  static_cast<int>(tot_snap_allfiles),
                                  Arguments.gist.gist_min_count);
  }
  auto time_pp_gist = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_pp_gist);

  // --- LLPS DBSCAN export ---
  if (dbscan_llps) {
    TitanUI::print_module_header("LLPS DBSCAN Cluster Analysis");
    dbscan_llps->export_timeseries(
        info_output.get_llps_dbscan_timeseries_filename());
    TitanUI::log_time(log, "LLPS DBSCAN", time_llps_dbscan.count() / 1000.0);
  }

  auto total_postproc = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_postproc);

  // --- Post-processing timing report ---
  TitanUI::log_section(log, "POST-PROCESSING TIMING");
  TitanUI::log_time(log, "GBSA/CAS/ML Postproc", time_pp_gbsa.count() / 1000.0);
  TitanUI::log_time(log, "TICA/MSM Solve+Export",
                    time_pp_tica.count() / 1000.0);
  TitanUI::log_time(log, "LLPS Flush+Postproc", time_pp_llps.count() / 1000.0);
  TitanUI::log_time(log, "DCCM/PCA Postproc", time_pp_dccm.count() / 1000.0);
  TitanUI::log_time(log, "LSP-MD Postproc", time_pp_lspmd.count() / 1000.0);
  TitanUI::log_time(log, "Water Bridges Export",
                    time_pp_water.count() / 1000.0);
  TitanUI::log_time(log, "GIST Thermodynamics", time_pp_gist.count() / 1000.0);
  log << "--------------------------------------------------------------------"
         "----------\n";
  TitanUI::log_time(log, "TOTAL POST-PROCESSING",
                    total_postproc.count() / 1000.0);

  TitanUI::log_section(log, "TITAN ENGINE TERMINATED SUCCESSFULLY");
  TitanUI::print_info(
      "\n\033[1;32mTITAN ENGINE TERMINATED SUCCESSFULLY\033[0m");
  double totalruntime = chrono::duration_cast<milliseconds>(
                            chrono::steady_clock::now() - total_run_time)
                            .count() /
                        1000.0;
  TitanUI::log_time(log, "GLOBAL TOTAL RUN TIME", totalruntime);
  log << "\nContact: nguyenhoanglinh9@duytan.edu.vn\n";
  log.close();
  return 0;
}
