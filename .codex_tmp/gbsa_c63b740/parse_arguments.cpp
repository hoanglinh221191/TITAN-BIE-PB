#include "MOL.h"
#include "arguments.h"
#include "functions.h"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
int arguments::parse_arguments(int argc, char **argv) {
  argc = argc;
  for (int i = 1; i < argc; ++i) {
    if (string(argv[i]) == "-gbsa") {
      flag_gbsa = true;
    }
    if (string(argv[i]) == "-alpb") {
      flag_alpb = true;
      flag_gbsa = true;
    }
    if (string(argv[i]) == "-cas") {
      flag_cas = true;
    }
    if (string(argv[i]) == "-network") {
      flag_network = true;
    }
    if (string(argv[i]) == "-pca") {
      flag_pca = true;
      flag_network = true; // PCA depends on CA data from network module
    }
    if (string(argv[i]) == "-ml") {
      flag_ml = true;
    }
    if (string(argv[i]).find("--fes-grid=") == 0) {
      fes_grid = stoi(string(argv[i]).substr(11));
    }
    if (string(argv[i]) == "-allostery") {
      flag_allostery = true;
      flag_network = true; // Allostery depends on DCCM trajectory data
    }
    if (string(argv[i]) == "-lspmd") {
      lspmd.flag_lspmd = true;
      // NOTE: does NOT auto-enable flag_network — LSP-MD runs standalone
      // without triggering DCCM matrix calculation.
    }
    if (string(argv[i]).find("-lspmd_cut=") == 0) {
      lspmd.lspmd_cutoff = stod(string(argv[i]).substr(11));
    }
    if (string(argv[i]).find("--allostery-tau=") == 0) {
      allostery_tau = stoi(string(argv[i]).substr(16));
    }
    if (string(argv[i]).find("--allostery-bins=") == 0) {
      allostery_bins = stoi(string(argv[i]).substr(17));
    }
    if (string(argv[i]) == "-tda") {
      flag_tda = true;
    }
    if (string(argv[i]) == "-tica") {
      msm.flag_tica = true;
      flag_network = true; // TICA depends on CA data from network module
    }
    if (string(argv[i]).find("--tica-lag=") == 0) {
      msm.tica_lag = stoi(string(argv[i]).substr(11));
    }
    if (string(argv[i]).find("--tica-components=") == 0) {
      msm.tica_components = stoi(string(argv[i]).substr(18));
    }
    // --- MSM Phase 2 flags ---
    if (string(argv[i]) == "-msm") {
      msm.flag_msm = true;
      msm.flag_tica = true;    // MSM depends on TICA projection
      flag_network = true; // TICA depends on CA data from network module
    }
    if (string(argv[i]) == "-msm-clusters" && i + 1 < argc) {
      msm.msm_clusters = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]) == "-msm-lag" && i + 1 < argc) {
      msm.msm_lag = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]) == "-msm-nts" && i + 1 < argc) {
      msm.msm_nts = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]) == "-msm-kmeans-iter" && i + 1 < argc) {
      msm.msm_kmeans_iter = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]) == "-msm-seed" && i + 1 < argc) {
      msm.msm_seed = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]) == "-msm-macrostates" && i + 1 < argc) {
      msm.msm_macrostates = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]) == "-msm-no-boundary-split") {
      msm.msm_auto_boundary = false;
    }
    if (string(argv[i]) == "-msm-no-pcca") {
      msm.msm_run_pcca = false;
    }
    if (string(argv[i]) == "-msm-no-its") {
      msm.msm_run_its = false;
    }
    if (string(argv[i]) == "-msm-its-lags" && i + 1 < argc) {
      msm.msm_its_lags.clear();
      istringstream iss{string(argv[i + 1])};
      for (string tok; getline(iss, tok, ',');)
        if (!tok.empty())
          msm.msm_its_lags.push_back(stoi(tok));
      ++i;
    }
    // --- LLPS flags ---
    if (string(argv[i]) == "-llps") {
      llps.flag_llps = true;
    }
    if (string(argv[i]) == "-llps-rcut" && i + 1 < argc) {
      llps.llps_rcut = stof(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-llps-perc" && i + 1 < argc) {
      llps.llps_perc_thr = stof(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-llps-maxlag" && i + 1 < argc) {
      llps.llps_maxlag = stoi(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-llps-dimer-size" && i + 1 < argc) {
      llps.llps_dimer_size = stoi(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-llps-contacts") {
      llps.llps_contacts = true;
    }
    if (string(argv[i]) == "-llps-sticker" && i + 1 < argc) {
      llps.llps_sticker_beads = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-llps-voxel" && i + 1 < argc) {
      llps.llps_voxel_size = stof(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]) == "-llps-dbscan-eps" && i + 1 < argc) {
      llps.llps_dbscan_eps = stod(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-llps-dbscan-minpts" && i + 1 < argc) {
      llps.llps_dbscan_minpts = stoi(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-llps-fuzzy") {
      llps.llps_fuzzy = true;
    }
    // --- GIST Water Thermodynamics flags ---
    if (string(argv[i]) == "-gist") {
      gist.flag_gist = true;
    }
    if (string(argv[i]).find("--gist-resolution=") == 0) {
      gist.gist_resolution = stod(string(argv[i]).substr(18));
    }
    if (string(argv[i]).find("--gist-padding=") == 0) {
      gist.gist_padding = stod(string(argv[i]).substr(15));
    }
    if (string(argv[i]).find("--gist-temperature=") == 0) {
      gist.gist_temperature = stod(string(argv[i]).substr(19));
    }
    if (string(argv[i]).find("--gist-water-model=") == 0) {
      gist.gist_water_model = string(argv[i]).substr(19);
    }
    if (string(argv[i]).find("--gist-min-count=") == 0) {
      gist.gist_min_count = stoi(string(argv[i]).substr(17));
    }
    if (string(argv[i]).find("--gist-cutoff=") == 0) {
      gist.gist_cutoff = stod(string(argv[i]).substr(14));
    }
    // --- IFP (Interaction Fingerprint) flags ---
    if (string(argv[i]) == "-ifp") {
      ifp.flag_ifp = true;
    }
    if (string(argv[i]) == "-ifp-hbond" && i + 1 < argc) {
      ifp.cutoff_hbond = stod(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-ifp-salt" && i + 1 < argc) {
      ifp.cutoff_salt_bridge = stod(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-ifp-hphob" && i + 1 < argc) {
      ifp.cutoff_hydrophobic = stod(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-ifp-pipi" && i + 1 < argc) {
      ifp.cutoff_pi_pi = stod(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-ifp-cpi" && i + 1 < argc) {
      ifp.cutoff_cation_pi = stod(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-ifp-hal" && i + 1 < argc) {
      ifp.cutoff_halogen = stod(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-ifp-out" && i + 1 < argc) {
      ifp.output_filename = string(argv[i + 1]);
      ++i;
    }
    // --- XDR trajectory (XTC/TRR) ---
    if (string(argv[i]) == "-traj" && i + 1 < argc) {
      traj_file = string(argv[i + 1]);
      flag_xdr = true;
      flag_1_pdb = true; // XDR always uses single-topology reference
      if (traj_file.length() > 4) {
        string ext = traj_file.substr(traj_file.length() - 4);
        if (ext == ".xtc")
          xdr_type = "xtc";
        else if (ext == ".trr")
          xdr_type = "trr";
      }
      ++i;
    }
    if (string(argv[i]) == "-igb" && i + 1 < argc) {
      igb = stoi(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-ext_diel" && i + 1 < argc) {
      ext_diel = stod(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-int_diel" && i + 1 < argc) {
      int_diel = stod(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-saltcon" && i + 1 < argc) {
      saltcon = stod(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-n_wat" && i + 1 < argc) {
      n_wat = stoi(argv[i + 1]);
      ++i;
    } else if (string(argv[i]).find("-n_wat=") == 0) {
      n_wat = stoi(string(argv[i]).substr(7));
    }
    if (string(argv[i]) == "-cscd") {
      cal_dis_sc_contact = true;
    }
    if (string(argv[i]) == "-opdb" || string(argv[i]) == "-opdbsep") {
      opdb = true;
      if (string(argv[i]) == "-opdbsep")
        flag_separation_results = true;
      if (i + 1 < argc) {
        if (string(argv[i + 1]).find(".pdb"))
          out_pdb_name = string(argv[i + 1]);
      }
    }
    if (string(argv[i]) == "-pqr") {
      flag_pqr = 1;
    }
    if (string(argv[i]) == "-interchain-energy") {
      flag_interchain_energy = true;
    }
    if (string(argv[i]) == "-dipole-moment") {
      flag_cal_dipole_moment = true;
    }
    if (string(argv[i]) == "-act") {
      flag_all_side_chain_ct = 0;
    }
    if (string(argv[i]) == "-useForceOptionforEnergy") {
      flag_energy_same_force = true;
    }
    if (string(argv[i]) == "-h" || string(argv[i]) == "--help") {
      help = true;
      if (help_print() == 0)
        return 2;
    }
    if (string(argv[i]) == "-v" || string(argv[i]) == "--version") {
      cout << "Build: " << get_build_number() << '\n';
      return 2;
    }
    if (string(argv[i]) == "-pbc") {
      pbc.solv_pbc = 1;
    }
    if (string(argv[i]) == "-cat-result") {
      flag_cat_result_multi_infile = true;
    }
    if (string(argv[i]) == "-pbci") {
      pbc.solv_pbc_int = 1;
      pbc.solv_pbc = 1;
      if (i + 1 < argc && isdigit(argv[i + 1][0])) {
        pbc.solv_pbc = 0;
        pbc.chain_ref = atoi(argv[i + 1]);
      } else {
        cout << "WARNING!! The option -pbci lacks reference chain!\n";
        pbc.solv_pbc_int = 0;
        pbc.solv_pbc = 0;
      }
    }
    if (string(argv[i]) == "-pbcn") {
      pbc.solv_pbcn = true;
    }
    if (string(argv[i]) == "-pbci2") {
      pbc.solv_pbc_int = 1;
      pbc.solv_pbc_int2 = 1;
      pbc.solv_pbc = 1;
      if (i + 1 < argc && isdigit(argv[i + 1][0])) {
        pbc.solv_pbc = 0;
        pbc.chain_ref = atoi(argv[i + 1]);
      } else {
        cout << "WARNING!! The option -pbci2 lacks reference chain!\n";
        pbc.solv_pbc_int = 0;
        pbc.solv_pbc_int2 = 0;
        pbc.solv_pbc = 0;
      }
    }
    if (string(argv[i]) == "-pbci3") {
      pbc.solv_pbc_int = 1;
      pbc.solv_pbc_int2 = 1;
      pbc.solv_pbc = 1;
      pbc.solv_pbc_int3 = 1;
      if (i + 1 < argc && isdigit(argv[i + 1][0])) {
        pbc.solv_pbc = 0;
        pbc.chain_ref = atoi(argv[i + 1]);
      } else {
        cout << "WARNING!! The option -pbci3 lacks reference chain!\n";
        pbc.solv_pbc_int = 0;
        pbc.solv_pbc_int2 = 0;
        pbc.solv_pbc_int3 = 0;
        pbc.solv_pbc = 0;
      }
    }
    if (string(argv[i]) == "-pbca") {
      pbc.solv_pbc_int = 1;
      pbc.solv_pbc = 1;
      pbc.pbc_aggressive = true;
      if (i + 1 < argc && isdigit(argv[i + 1][0])) {
        pbc.solv_pbc = 0;
        pbc.chain_ref = atoi(argv[i + 1]);
        if (pbc.chain_ref < 0)
          pbc.chain_ref = 0;
      } else {
        pbc.solv_pbc_int = 0;
        pbc.solv_pbc = 0;
      }
    }
    if (string(argv[i]) == "-pbcclust") {
      pbc.solv_pbc_clust = true;
      if (i + 1 < argc && string(argv[i + 1]).find("-") == string::npos) {
        string temp{string(argv[i + 1])};
        for (size_t i{}; i < temp.size(); ++i)
          pbc.pbc_clust_axis.emplace_back(temp[i]);
      } else
        pbc.solv_pbc_clust = false;
    }
    if (string(argv[i]) == "-1pdb") {
      flag_1_pdb = true;
    }
    if (string(argv[i]) == "-hb") {
      flag_h_bond = 1;
      if (argc > i + 1 && string(argv[i + 1]) == "side-chain")
        flag_sidechain_hbond = true;
    }
    if (string(argv[i]) == "-entropy") {
      entropy = true;
    }
    if (string(argv[i]) == "-nocal") {
      flag_nocal = true;
      flag_1_pdb = true;
    }
    if (string(argv[i]) == "-all-subfolders") {
      flag_all_folders = true;
      flag_separation_results = true;
    }
    if (string(argv[i]) == "-split") {
      flag_separation_results = true;
    }
    if (string(argv[i]) == "-cal-force") {
      flag_cal_force = true;
    }
    if (string(argv[i]) == "-cal-volume") {
      flag_cal_vol = true;
    }
    if (string(argv[i]) == "-cache") {
      flag_cache = 1;
      if (i + 1 < argc &&
          (string(argv[i + 1]) == "-no" || string(argv[i + 1]) == "no"))
        flag_cache = 2;
    }
    if (string(argv[i]) == "-max-force") {
      flag_maxforce = true;
    }
    if (string(argv[i]) == "-new-emap") {
      flag_new_emap = true;
    }
    if (string(argv[i]) == "-energy-res") {
      print_energy_res = true;
      energy_map = 1;
    }
    if (string(argv[i]) == "-force-res") {
      print_force_res = true;
      force_map = 1;
    }
    if (string(argv[i]) == "-use-CI") {
      flag_use_CI = true;
    }
    if (string(argv[i]) == "-elec-screening") {
      flag_elec_screening = true;
    }
    if (string(argv[i]) == "-charged_atom_dis") {
      flag_charged_atom_distance_dis = true;
      if (i + 1 < argc && argv[i + 1][0] != '-' && isdigit(argv[i + 1][0]) != 0)
        number_charged_distances_bin = atoi(argv[i + 1]);
    }
    if (string(argv[i]) == "-charged_distances_min") {
      flag_charged_atom_distance_dis = true;
      if (i + 1 < argc && argv[i + 1][0] != '-' && isdigit(argv[i + 1][0]) != 0)
        charged_distances_min = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-charged_distances_max") {
      flag_charged_atom_distance_dis = true;
      if (i + 1 < argc && argv[i + 1][0] != '-' && isdigit(argv[i + 1][0]) != 0)
        charged_distances_max = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-resdis2com") {
      flag_cal_dis_res_2com = true;
    }
    if (string(argv[i]) == "-angle_dipole_efield") {
      flag_angle_dipole_efield = true;
      if (i + 1 < argc && argv[i + 1][0] != '-' && isdigit(argv[i + 1][0]) != 0)
        number_bin_angle_dipole_efield = atoi(argv[i + 1]);
    }
    if (string(argv[i]) == "-side_chain_contact" && i + 1 < argc) {
      if (string(argv[i + 1]) == "full-res")
        flag_sidechain = false;
      else
        flag_sidechain = true;
    }
    if (string(argv[i]) == "-noprint-side_chain_contact" && i + 1 < argc) {
      flag_print_sidechain = false;
    }
    if (string(argv[i]) == "-mol-1-type" && i + 1 < argc) {
      std::string t = string(argv[i + 1]);
      if (t == "protein" || t == "complex" || t == "lig" || t == "rna" ||
          t == "dna")
        mol_1_type = t;
    }
    if (string(argv[i]) == "-mol-2-type" && i + 1 < argc) {
      std::string t = string(argv[i + 1]);
      if (t == "protein" || t == "protein-user" || t == "complex" ||
          t == "lig" || t == "rna" || t == "dna")
        mol_2_type = t;
    }
    if (string(argv[i]) == "-com_energy" && i + 1 < argc) {
      if (string(argv[i + 1]) == "x")
        com_dis_energy_2mol = 1;
      else if (string(argv[i + 1]) == "y")
        com_dis_energy_2mol = 2;
      else if (string(argv[i + 1]) == "z")
        com_dis_energy_2mol = 3;
      else if (string(argv[i + 1]) == "all")
        com_dis_energy_2mol = 0;
      else
        com_dis_energy_2mol = 4;
    }
    if (string(argv[i]) == "-b" && i + 1 < argc) {
      begin_time = atof(argv[i + 1]);
      flag_begin_set = true;
    }
    if (string(argv[i]) == "-e" && i + 1 < argc) {
      end_time = atof(argv[i + 1]);
      flag_end_set = true;
    }
    if (string(argv[i]) == "-dt" && i + 1 < argc) {
      dt = atoi(argv[i + 1]);
    }
    if (string(argv[i]) == "-mol2" && i + 1 < argc) {
      begin_mol_2 = atoi(argv[i + 1]);
    }
    if (string(argv[i]) == "-topo" && i + 1 < argc) {
      topo_name = argv[i + 1];
      if (check_topo(topo_name, 1) == true)
        flag_topo_1 = true;
    }
    if (string(argv[i]) == "-topo2" && i + 1 < argc) {
      topo_name_2 = argv[i + 1];
      if (check_topo(topo_name_2, 2) == true)
        flag_topo_2 = true;
    }
    if (string(argv[i]) == "-f" && i + 1 < argc) {
      pdb_name = argv[i + 1];
      flag_1_pdb = true;
    }
    // -s = structure file (topology reference for XDR trajectory)
    // Sets pdb_name as reference and enables XDR mode.
    // XDR type is auto-detected from -traj or -i file list.
    if (string(argv[i]) == "-s" && i + 1 < argc) {
      pdb_name = argv[i + 1];
      flag_1_pdb = true;
      flag_xdr = true;
    }
    if (string(argv[i]) == "-n" && i + 1 < argc) {
      index_name = argv[i + 1];
      flag_index = true;
    }
    if (string(argv[i]) == "-eu" && i + 1 < argc) {
      eu = argv[i + 1];
    }
    if (string(argv[i]) == "-fu" && i + 1 < argc) {
      fu = argv[i + 1];
    }

    if (string(argv[i]) == "-ff" && i + 1 < argc) {
      force_f = argv[i + 1];
      transform(force_f.begin(), force_f.end(), force_f.begin(), ::toupper);
    }
    if (string(argv[i]) == "-epsilon_in" && i + 1 < argc) {
      epsilon_in = atof(argv[i + 1]);
      ++argc;
    }
    if (string(argv[i]) == "-g" && i + 1 < argc) {
      auto found = string(argv[i + 1]).find("&");
      auto found_2 = string(argv[i + 1]).find("-");
      if (found != std::string::npos) {
        group_1 = stoi(string(argv[i + 1]).substr(0, found));
        group_2 = stoi(string(argv[i + 1]).substr(found + 1));
      } else if (found_2 != std::string::npos) {
        group_1 = stoi(string(argv[i + 1]).substr(0, found_2));
        group_2 = stoi(string(argv[i + 1]).substr(found_2 + 1));
      } else {
        cout << "Your argument -g is invalid!\n";
      }
      if (group_1 == group_2)
        self_int = true;
    }
    if (string(argv[i]) == "-sg" && i + 1 < argc) {
      // Format: -sg=<mol1_groups>-<mol2_groups> [extra_mol2_groups...]
      // or:     -sg=<mol1_groups>&<mol2_groups> [extra_mol2_groups...]
      // The separator between mol1 and mol2 groups is '-' or '&'.
      // Subsequent non-flag tokens are additional mol2 groups.
      string spec = string(argv[i + 1]);
      auto sep = spec.find("&");
      if (sep == string::npos)
        sep = spec.find("-");
      if (sep == string::npos) {
        cout << "Your argument -sg is invalid!\n"
             << "Use - or & to separate groups of two molecules!\n";
      } else {
        string mol1_spec = spec.substr(0, sep);
        string mol2_first = spec.substr(sep + 1);

        regex reg{"(&?)([[:blank:]]*)(\\d*)([[:blank:]]*)(&*)"};

        // Parse mol1 groups
        sregex_iterator pos(mol1_spec.cbegin(), mol1_spec.cend(), reg), end;
        for (; pos != end; ++pos) {
          if (!pos->str(3).empty() && isdigit(pos->str(3)[0]))
            sub_group_1.emplace_back(stoi(pos->str(3)));
        }

        // Parse mol2 groups from the first token (after separator)
        sregex_iterator pos_2(mol2_first.cbegin(), mol2_first.cend(), reg);
        for (; pos_2 != end; ++pos_2) {
          if (!pos_2->str(3).empty() && isdigit(pos_2->str(3)[0]))
            sub_group_2.emplace_back(stoi(pos_2->str(3)));
        }

        // Consume additional mol2 groups from subsequent non-flag tokens
        while (i + 2 < argc && argv[i + 2][0] != '-') {
          ++i;
          try {
            sub_group_2.emplace_back(stoi(string(argv[i + 1])));
          } catch (...) {
            break;
          }
        }


        flag_sub_groups = true;
      }
    }
    if (string(argv[i]) == "-emap") {
      energy_map = 1;
      if (i + 1 < argc && (string(argv[i + 1]) == "residue-atom" ||
                           string(argv[i + 1]) == "residue-residue" ||
                           string(argv[i + 1]) == "atom-atom" ||
                           string(argv[i + 1]) == "atom-residue"))
        emap_type = string(argv[i + 1]);
    }
    if (string(argv[i]) == "-emap-level" && i + 1 < argc) {
      if (string(argv[1 + 1]) == "all" || string(argv[i + 1]) == "per-res")
        emap_level = string(argv[i + 1]);
    }
    if (string(argv[i]) == "-fmap") {
      force_map = 1;
      if (i + 1 < argc && (string(argv[i + 1]) == "residue-atom" ||
                           string(argv[i + 1]) == "residue-residue" ||
                           string(argv[i + 1]) == "atom-atom" ||
                           string(argv[i + 1]) == "atom-residue"))
        fmap_type = string(argv[i + 1]);
    }
    if (string(argv[i]) == "-fmap-level" && i + 1 < argc) {
      if (string(argv[1 + 1]) == "all" || string(argv[i + 1]) == "per-res")
        fmap_level = string(argv[i + 1]);
    }
    if (string(argv[i]) == "-nt" && i + 1 < argc) {
      if (stoi(string(argv[i + 1])) < 0)
        nt = 0;
      else
        nt = stoi(string(argv[i + 1]));
    }
    if (string(argv[i]) == "-dsc" && i + 1 < argc) {
      d_standard = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-denfo" && i + 1 < argc) {
      distance_energy_cutoff = atof(argv[i + 1]);
      distance_energy_cutoff *= 0.1;
    }
    if (string(argv[i]) == "-path-a" && i + 1 < argc) {
      path_amino_acid = argv[i + 1];
      if (path_amino_acid[path_amino_acid.length() - 1] != '/')
        path_amino_acid += '/';
    }
    if (string(argv[i]) == "-path-n" && i + 1 < argc) {
      path_nonbonded = argv[i + 1];
      if (path_nonbonded[path_nonbonded.length() - 1] != '/')
        path_nonbonded += '/';
    }
    if (string(argv[i]) == "-pattern-input" && i + 1 < argc) {
      flag_pattern_input = true;
      flag_1_pdb = true;
      if (string(argv[i + 1]).find("*") != string::npos) {
        input_pattern = string(argv[i + 1]);
        if (i + 2 < argc && argv[i + 2][0] != '-')
          b_input = atoi(argv[i + 2]);
        if (i + 3 < argc && argv[i + 3][0] != '-')
          e_input = atoi(argv[i + 3]);
        if (i + 4 < argc && argv[i + 4][0] != '-')
          step_input = atoi(argv[i + 4]);
        if (e_input != 10000000)
          flag_range_input = true;
      } else
        cout << "Can not parse -pattern-input option, it must contain * "
                "character!\n";
    }
    if (string(argv[i]) == "-ahb" && i + 1 < argc) {
      angle_h = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-dhb" && i + 1 < argc) {
      dis_h = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-entropy-quasi" && i + 3 < argc) {
      entropy_quasi = true;
      if (isdigit(argv[i + 1][0]) && isdigit(argv[i + 2][0]) &&
          isdigit(argv[i + 3][0])) {
        group_entropy_com = atoi(argv[i + 1]);
        group_entropy_lg = atoi(argv[i + 2]);
        group_entropy_rec = atoi(argv[i + 3]);
      } else {
        cout << "-entropy-quasi option is invalid\n";
        entropy_quasi = false;
      }
    }
    if (string(argv[i]) == "-a" && string(argv[i]) != "-alpb" && i + 1 < argc) {
      arg_filename = string(argv[i + 1]);
      parse_arg_file();
      return 0;
    }
    if (string(argv[i]) == "-i" && i + 1 < argc) {
      list_filename = string(argv[i + 1]);
      fstream test_listFile(list_filename, ios::in);
      if (test_listFile.is_open())
        flag_listfile = true;
      else {
        flag_listfile = false;
        cout << "Cannot open " << list_filename << " file!\n";
        exit(EXIT_FAILURE);
      }
      test_listFile.close();
    }
    if (string(argv[i]) == "-water-type" && i + 1 < argc) {
      if (string(argv[i + 1]) == "surface" && string(argv[i + 1]) == "sphere")
        water_type = string(argv[i + 1]);
      else {
        cout << "WARNING!!! -water-type option is incorrect!\n";
      }
    }
    if (string(argv[i]) == "-dis-water" && i + 1 < argc) {
      d_crit_water = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-all-respair" && i + 1 < argc) {
      flag_all_respair = true;
    }
    if (string(argv[i]) == "-pull-direct") {
      de.flag_find_pulling_direct = true;
    }
    if (string(argv[i]) == "-de-diff-weight" && i + 1 < argc) {
      de.de_diff_weight = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-de-cross-prob" && i + 1 < argc) {
      de.de_cross_prob = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-de-F" && i + 1 < argc) {
      de.de_F = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-de-scan-radius" && i + 1 < argc) {
      de.de_scan_radius = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-de-min-scan-radius" && i + 1 < argc) {
      de.de_min_scan_radius = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-de-score-function" && i + 1 < argc) {
      de.de_flag_score_function = atoi(argv[i + 1]);
    }
    if (string(argv[i]) == "-de-pop-num" && i + 1 < argc) {
      de.de_pop_num = stoi(string(argv[i + 1]));
    }
    if (string(argv[i]) == "-de-ite-num" && i + 1 < argc) {
      de.de_ite = stoi(string(argv[i + 1]));
    }
    if (string(argv[i]) == "-de-algorithm" && i + 1 < argc) {
      de.de_algorithm = stoi(string(argv[i + 1]));
    }
    if (string(argv[i]) == "-de-angle" && i + 1 < argc) {
      de.de_restrain_angle = atof(argv[i + 1]);
    }
    if (string(argv[i]) == "-de-auto-term") {
      de.flag_de_auto_term = true;
    }
    if (string(argv[i]) == "-de-adaptive-F") {
      de.de_adaptive_F = true;
    }
    if (string(argv[i]) == "-de-stagnation-limit" && i + 1 < argc) {
      de.de_stagnation_limit = stoi(string(argv[i + 1]));
    }
    // --- Flags previously only in arg-file parser, now unified ---
    if (string(argv[i]) == "-new_al" || string(argv[i]) == "-new-al") {
      new_al_1pdb = true;
    }
    if (string(argv[i]) == "-refpbc" && i + 1 < argc) {
      pbc.solv_pbc_int = 1;
      pbc.chain_ref = 0;
      pbc.ref_pdb_name = string(argv[i + 1]);
      if (pbc.ref_pdb_name.length() <= 4 ||
          pbc.ref_pdb_name.substr(pbc.ref_pdb_name.length() - 4) != ".pdb")
        pbc.ref_pdb_name += ".pdb";
      pbc.flag_pbc_ref = true;
      ++i;
    }
    if (string(argv[i]) == "-pbcm") {
      pbc.solv_pbc_int = 1;
      pbc.solv_pbc = 1;
      pbc.pbc_mem = true;
      if (i + 1 < argc && isdigit(argv[i + 1][0])) {
        pbc.solv_pbc = 0;
        pbc.chain_ref = atoi(argv[i + 1]);
        ++i;
      }
    }
    if (string(argv[i]) == "-pbc-envi") {
      while (i + 1 < argc && argv[i + 1][0] != '-') {
        string t = string(argv[++i]);
        if (t == "@")
          t = " ";
        pbc.pbc_envi_chain.emplace_back(t);
      }
      if (!pbc.pbc_envi_chain.empty())
        pbc.pbc_envi = true;
    }
    if (string(argv[i]) == "-com-dis" && i + 1 < argc) {
      string axis = string(argv[i + 1]);
      if (axis == "x")
        com_dis_2mol = 1;
      else if (axis == "y")
        com_dis_2mol = 2;
      else if (axis == "z")
        com_dis_2mol = 3;
      else
        com_dis_2mol = 4;
      ++i;
    }
    if (string(argv[i]) == "-cal-mindis") {
      flag_cal_mindis = true;
    }
    if ((string(argv[i]) == "-common_topo" ||
         string(argv[i]) == "-common-topo") &&
        i + 1 < argc) {
      flag_use_common_topo = true;
      if (flag_topo_1 == false && flag_topo_2 == false) {
        topo_name = string(argv[i + 1]);
        topo_name_2 = string(argv[i + 1]);
        flag_topo_1 = true;
        flag_topo_2 = true;
        mol_1_type = "protein-user";
        mol_2_type = "protein-user";
      }
      ++i;
    }
    if (string(argv[i]) == "-contact-map" && i + 1 < argc) {
      string val = string(argv[i + 1]);
      if (val == "no")
        flag_cal_contact = false;
      ++i;
    }
    // --- DHAM (Dynamic Histogram Analysis Method) flags ---
    if (string(argv[i]) == "-dham") {
      dham.flag_dham = true;
      msm.flag_msm = true;      // DHAM depends on MSM (K-means + dtrajs)
      msm.flag_tica = true;     // MSM depends on TICA projection
      flag_network = true;      // TICA depends on CA data from network module
    }
    if (string(argv[i]) == "-dham-temps" && i + 1 < argc) {
      dham.temps_file = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-dham-ener" && i + 1 < argc) {
      dham.ener_list_file = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-dham-target-temp" && i + 1 < argc) {
      dham.target_temp = stod(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--dham-max-iter=") == 0) {
      dham.max_iter = stoi(string(argv[i]).substr(16));
    }
    if (string(argv[i]).find("--dham-tol=") == 0) {
      dham.tolerance = stod(string(argv[i]).substr(11));
    }
    if (string(argv[i]).find("--dham-n-replicas=") == 0) {
      dham.n_replicas = stoi(string(argv[i]).substr(18));
    }
    if (string(argv[i]) == "-dham-dt-energy" && i + 1 < argc) {
      dham.dt_energy_ps = stod(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]) == "-no-center") {
      pbc.no_center = true;
    }
    command += argv[i];
    command += " ";
  }

  if (flag_topo_1 && flag_topo_2 && topo_name == topo_name_2) {
    flag_use_common_topo = true;
  }

  return 0;
}
// ─── Tokenizer: convert arg-file lines into CLI-style tokens ───
std::vector<std::string>
arguments::tokenize_arg_file(const std::string &filename) {
  std::vector<std::string> tokens;
  std::fstream in_file(filename, std::ios::in);
  if (!in_file.is_open()) {
    std::cout << "Parameter file " << filename << " is not found!\n";
    return tokens;
  }
  for (std::string line; std::getline(in_file, line);) {
    // Strip comments
    auto hash_pos = line.find('#');
    if (hash_pos != std::string::npos)
      line = line.substr(0, hash_pos);
    // Strip CR
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    // Trim whitespace
    size_t start = line.find_first_not_of(" \t");
    if (start == std::string::npos)
      continue;
    line = line.substr(start);
    size_t end = line.find_last_not_of(" \t");
    if (end != std::string::npos)
      line = line.substr(0, end + 1);
    if (line.empty())
      continue;

    auto eq_pos = line.find('=');
    if (eq_pos != std::string::npos && line[0] == '-') {
      std::string flag = line.substr(0, eq_pos);
      std::string value = line.substr(eq_pos + 1);
      // Trim value
      size_t vs = value.find_first_not_of(" \t");
      if (vs != std::string::npos)
        value = value.substr(vs);
      size_t ve = value.find_last_not_of(" \t");
      if (ve != std::string::npos)
        value = value.substr(0, ve + 1);

      // Split flag and value into separate tokens (single-dash only in arg.in)
      tokens.push_back(flag);
      if (!value.empty()) {
        std::istringstream iss(value);
        std::string word;
        while (iss >> word)
          tokens.push_back(word);
      }
    } else {
      // No '=', push as-is (boolean flags like -gbsa, -entropy)
      tokens.push_back(line);
    }
  }
  in_file.close();
  return tokens;
}

// ─── Unified arg-file parser: tokenize → fake argv → parse_arguments() ───
int arguments::parse_arg_file() {
  auto tokens = tokenize_arg_file(arg_filename);
  if (tokens.empty())
    return 1;

  file_par = true;

  // Build fake argv[]
  std::vector<std::string> argv_storage;
  argv_storage.push_back("titan"); // argv[0]
  for (auto &t : tokens)
    argv_storage.push_back(t);

  std::vector<char *> argv_ptrs;
  for (auto &s : argv_storage)
    argv_ptrs.push_back(const_cast<char *>(s.c_str()));

  int result = parse_arguments(static_cast<int>(argv_ptrs.size()),
                               argv_ptrs.data());

  // Resolve relative paths for topo/index files relative to arg file dir
  std::string arg_dir;
  if (arg_filename.find('\\') != std::string::npos ||
      arg_filename.find('/') != std::string::npos) {
    arg_dir = arg_filename.substr(0, arg_filename.find_last_of("\\/"));
  }
  if (!arg_dir.empty()) {
    auto resolve = [&](std::string &path) {
      if (!path.empty() && path != "no" &&
          path.find('\\') == std::string::npos &&
          path.find('/') == std::string::npos) {
        path = arg_dir + "/" + path;
      }
    };
    if (flag_topo_1)
      resolve(topo_name);
    if (flag_topo_2)
      resolve(topo_name_2);
    if (flag_index)
      resolve(index_name);
  }

  return result;
}

int arguments::parse_arguments_interactive(int argc_in) {
  int flag;
  string force_f;
  do {
    cout << "Please choose the force field that you want to use:" << '\n';
    cout << "1. AMBER99SB" << '\n';
    cout << "2. AMBER99SB-ILDN" << '\n';
    cout << "3. AMBER99SB-BSC1" << '\n';
    cout << "4. AMBER99SB-DISP" << '\n';
    cout << "5. AMBER14SB" << '\n';
    cout << "6. AMBER19SB" << '\n';
    cout << "7. GROMOS43A1" << '\n';
    cout << "8. GROMOS53A6" << '\n';
    cout << "9. GROMOS54A7" << '\n';
    cout << "10. CHARMM27" << '\n';
    cout << "11. CHARMM36" << '\n';
    cout << "12. OPLSAA" << '\n';
    getline(cin, force_f);
    istringstream st(force_f);
    st >> flag;
  } while (flag > 12 || flag < 1);
  if (flag == 1) {
    force_f = "AMBER99SB";
  } else if (flag == 2) {
    force_f = "AMBER99SB-ILDN";
  } else if (flag == 3) {
    force_f = "AMBER99SB-BSC1";
  } else if (flag == 4) {
    force_f = "AMBER99SB-DISP";
  } else if (flag == 5) {
    force_f = "AMBER14SB";
  } else if (flag == 6) {
    force_f = "AMBER19SB";
  } else if (flag == 7) {
    force_f = "GROMOS43A1";
  } else if (flag == 8) {
    force_f = "GROMOS53A6";
  } else if (flag == 9) {
    force_f = "GROMOS54A7";
  } else if (flag == 10) {
    force_f = "CHARMM27";
  } else if (flag == 11) {
    force_f = "CHARMM36";
  } else
    force_f = "OPLSAA";

  do {
    cout
        << "Do you want to create energy map? (yes = 1, no = 2, or custom = 3)";
    if (!(cin >> energy_map)) {
      cin.clear();
      cin.ignore(10000, '\n');
    }
  } while (energy_map != 1 && energy_map != 2 && energy_map != 3);
  if (energy_map == 3) {
    cout << "Please specify type of energy map calculation: (residue-atom | "
            "atom-atom | residue-residue | atom-residue)";
    string temp{};
    getline(cin, temp);
    string::size_type pos{temp.find("-")};
    if (pos != string::npos &&
        (temp.substr(0, pos) == "residue" || temp.substr(0, pos) == "atom") &&
        (temp.substr(pos + 1) == "residue" || temp.substr(pos + 1) == "atom")) {
      emap_type = temp;
    } else {
      cout << "Cannot parse! will use default value (residue-atom)\n";
    }
  }

  do {
    cout << "Do you want to take only minimun sc-ct distance? (Y/N=1/2)";
    if (!(cin >> flag_all_side_chain_ct)) {
      cin.clear();
      cin.ignore(10000, '\n');
    }
  } while (flag_all_side_chain_ct != 1 && flag_all_side_chain_ct != 2);
  do {
    cout << "Do you want to creat .PQR file? (Y/N=1/2)";
    cin >> flag_pqr;
  } while (flag_pqr != 1 && flag_pqr != 2);

  cout << "Do you want to use topo file (no if you don't want to use topo): ";
  cin >> topo_name;
  if (check_topo(topo_name, 1) == true)
    flag_topo_1 = true;

  cout << "Number of threads: ";
  cin >> nt;
  return 0;
}
bool arguments::check_topo(string NameofTopoFile, int Topo) {
  if (NameofTopoFile == "no")
    return false;
  else {
    fstream input_topo_protein(NameofTopoFile, ios::in);
    bool topo_is_open{input_topo_protein.is_open()};
    input_topo_protein.close();
    if (topo_is_open == true)
      return true;
    else {
      while (!topo_is_open) {
        cout << "Cannot open " << NameofTopoFile
             << ". Please type your lig topo file (no if you don't want to use "
                "topo): ";
        cin >> NameofTopoFile;
        cout << NameofTopoFile << "\n";
        topo_name = NameofTopoFile;
        if (NameofTopoFile == "no")
          return false;
        else {
          input_topo_protein.open(NameofTopoFile, ios::in);
          if (input_topo_protein.is_open()) {
            input_topo_protein.close();
            return true;
          } else {
            input_topo_protein.close();
            continue;
          }
        }
      }
      return false;
    }
  }
}
