#include "MOL.h"
#include "arguments.h"
#include "functions.h"
#include "titan_ui.h"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
int arguments::parse_arguments(int argc, char **argv) {
  argc = argc;
  auto normalize_mol_type = [](std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (value == "ligand")
      value = "lig";
    return value;
  };

  auto assign_mol_type = [&](const std::string &value, bool mol1) {
    const std::string t = normalize_mol_type(value);
    if (t != "protein" && t != "protein-user" && t != "complex" &&
        t != "lig" && t != "rna" && t != "dna")
      return false;

    if (mol1) {
      if (t == "protein-user")
        return false;
      mol_1_type = t;
      mol_1_type_explicit = true;
    } else {
      mol_2_type = t;
      mol_2_type_explicit = true;
    }
    return true;
  };
  auto parse_allostery_cut_mode = [&](std::string raw_mode) {
    std::transform(raw_mode.begin(), raw_mode.end(), raw_mode.begin(),
                   [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (raw_mode == "avg" || raw_mode == "all" || raw_mode == "mean") {
      allostery_use_avg_cutoff = true;
      return;
    }
    if (raw_mode == "first" || raw_mode == "frame0" || raw_mode == "instant") {
      allostery_use_avg_cutoff = false;
      return;
    }
    TitanUI::print_error(
        "Allostery: -allostery-cut-mode must be 'avg' or 'first'.");
    exit(EXIT_FAILURE);
  };

  for (int i = 1; i < argc; ++i) {
    const std::string arg = string(argv[i]);

    if (arg.rfind("-mol-1-type=", 0) == 0) {
      assign_mol_type(arg.substr(12), true);
      continue;
    }
    if (arg.rfind("-mol-2-type=", 0) == 0) {
      assign_mol_type(arg.substr(12), false);
      continue;
    }

    if (arg.rfind("-community=", 0) == 0) {
      community.flag_community = true;
      community.method = arg.substr(11);
      flag_network = true;
      continue;
    }
    if (arg == "-community") {
      TitanUI::print_error(
          "Community: use -community=louvain (single token, no spaces).");
      exit(EXIT_FAILURE);
    }
    if (arg.rfind("--community", 0) == 0) {
      TitanUI::print_error(
          "Community: use -community=louvain (single dash, no spaces).");
      exit(EXIT_FAILURE);
    }

    if (arg == "-gbsa") {
      flag_gbsa = true;
    }
    if (arg == "-alpb") {
      flag_alpb = true;
      flag_gbsa = true;
    }
    if (arg == "-cas") {
      flag_cas = true;
    }
    if (arg == "-dssp") {
      dssp.flag_dssp = true;
    }
    if (arg == "-no-dssp") {
      dssp.flag_dssp = false;
    }
    if (arg == "-rms" || arg == "-rmsd" || arg == "-rmsf") {
      rms.flag_rms = true;
    }
    if (arg == "-lj") {
      lj.flag_lj = true;
    }
    if (arg.rfind("-lj-cut=", 0) == 0 ||
        arg.rfind("--lj-cut=", 0) == 0) {
      lj.lj_cluster_cutoff = stod(arg.substr(arg.find('=') + 1));
    }
    if ((arg == "-lj-cut" || arg == "--lj-cut") && i + 1 < argc) {
      lj.lj_cluster_cutoff = stod(string(argv[i + 1]));
      ++i;
    }
    if (arg.rfind("-lj-min-cluster=", 0) == 0 ||
        arg.rfind("--lj-min-cluster=", 0) == 0) {
      lj.lj_min_cluster_size = stoi(arg.substr(arg.find('=') + 1));
    }
    if ((arg == "-lj-min-cluster" || arg == "--lj-min-cluster") &&
        i + 1 < argc) {
      lj.lj_min_cluster_size = stoi(string(argv[i + 1]));
      ++i;
    }
    if (arg.rfind("-lj-rho-bin=", 0) == 0 ||
        arg.rfind("--lj-rho-bin=", 0) == 0) {
      lj.lj_rho_bin = stod(arg.substr(arg.find('=') + 1));
    }
    if ((arg == "-lj-rho-bin" || arg == "--lj-rho-bin") &&
        i + 1 < argc) {
      lj.lj_rho_bin = stod(string(argv[i + 1]));
      ++i;
    }
    if (arg.rfind("-lj-rho-max=", 0) == 0 ||
        arg.rfind("--lj-rho-max=", 0) == 0) {
      lj.lj_rho_max = stod(arg.substr(arg.find('=') + 1));
    }
    if ((arg == "-lj-rho-max" || arg == "--lj-rho-max") &&
        i + 1 < argc) {
      lj.lj_rho_max = stod(string(argv[i + 1]));
      ++i;
    }
    if (arg.rfind("-lj-solid-q6=", 0) == 0 ||
        arg.rfind("--lj-solid-q6=", 0) == 0) {
      lj.lj_solid_q6 = stod(arg.substr(arg.find('=') + 1));
    }
    if ((arg == "-lj-solid-q6" || arg == "--lj-solid-q6") && i + 1 < argc) {
      lj.lj_solid_q6 = stod(string(argv[i + 1]));
      ++i;
    }
    if (arg.rfind("-lj-solid-min-cn=", 0) == 0 ||
        arg.rfind("--lj-solid-min-cn=", 0) == 0) {
      lj.lj_solid_min_cn = stoi(arg.substr(arg.find('=') + 1));
    }
    if ((arg == "-lj-solid-min-cn" || arg == "--lj-solid-min-cn") &&
        i + 1 < argc) {
      lj.lj_solid_min_cn = stoi(string(argv[i + 1]));
      ++i;
    }
    if (arg.rfind("-lj-crystal-tol=", 0) == 0 ||
        arg.rfind("--lj-crystal-tol=", 0) == 0) {
      lj.lj_crystal_tol = stod(arg.substr(arg.find('=') + 1));
    }
    if ((arg == "-lj-crystal-tol" || arg == "--lj-crystal-tol") &&
        i + 1 < argc) {
      lj.lj_crystal_tol = stod(string(argv[i + 1]));
      ++i;
    }
    if (arg.rfind("-dssp-hbond-cutoff=", 0) == 0 ||
        arg.rfind("--dssp-hbond-cutoff=", 0) == 0) {
      dssp.hbond_energy_cutoff = stod(arg.substr(arg.find('=') + 1));
    }
    if ((arg == "-dssp-hbond-cutoff" || arg == "--dssp-hbond-cutoff") &&
        i + 1 < argc) {
      dssp.hbond_energy_cutoff = stod(string(argv[i + 1]));
      ++i;
    }
    if (arg == "-network") {
      flag_network = true;
      flag_network_explicit = true;
    }
    if (string(argv[i]) == "-pca") {
      flag_pca = true;
      flag_network = true; // PCA depends on CA data from network module
    }
    if (string(argv[i]) == "-ml") {
      flag_ml = true;
    }
    if (string(argv[i]).find("--fes-grid=") == 0 ||
        string(argv[i]).find("-fes-grid=") == 0) {
      fes_grid = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
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
    if (string(argv[i]).find("--allostery-tau=") == 0 ||
        string(argv[i]).find("-allostery-tau=") == 0) {
      allostery_tau =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]).find("--allostery-bins=") == 0 ||
        string(argv[i]).find("-allostery-bins=") == 0) {
      allostery_bins =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]).find("--allostery-knn=") == 0 ||
        string(argv[i]).find("-allostery-knn=") == 0) {
      allostery_knn =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--allostery-knn" ||
         string(argv[i]) == "-allostery-knn") &&
        i + 1 < argc) {
      allostery_knn = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--allostery-cut=") == 0 ||
        string(argv[i]).find("-allostery-cut=") == 0) {
      allostery_cutoff =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--allostery-cut" ||
         string(argv[i]) == "-allostery-cut") &&
        i + 1 < argc) {
      allostery_cutoff = stod(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--allostery-cut-mode=") == 0 ||
        string(argv[i]).find("-allostery-cut-mode=") == 0) {
      parse_allostery_cut_mode(
          string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--allostery-cut-mode" ||
         string(argv[i]) == "-allostery-cut-mode") &&
        i + 1 < argc) {
      parse_allostery_cut_mode(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]) == "--allostery-no-te" ||
        string(argv[i]) == "-allostery-no-te") {
      allostery_compute_te = false;
    }
    if (string(argv[i]) == "-tda") {
      flag_tda = true;
    }
    if (string(argv[i]) == "-tica") {
      msm.flag_tica = true;
      flag_network = true; // TICA depends on CA data from network module
    }
    if (string(argv[i]).find("--tica-lag=") == 0 ||
        string(argv[i]).find("-tica-lag=") == 0) {
      msm.tica_lag =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--tica-lag" || string(argv[i]) == "-tica-lag") &&
        i + 1 < argc) {
      msm.tica_lag = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--tica-components=") == 0 ||
        string(argv[i]).find("-tica-components=") == 0) {
      msm.tica_components =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--tica-components" ||
         string(argv[i]) == "-tica-components") &&
        i + 1 < argc) {
      msm.tica_components = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--tica-backend=") == 0 ||
        string(argv[i]).find("-tica-backend=") == 0) {
      msm.tica_backend =
          string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if ((string(argv[i]) == "--tica-backend" ||
         string(argv[i]) == "-tica-backend") &&
        i + 1 < argc) {
      msm.tica_backend = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]).find("--tica-rank=") == 0 ||
        string(argv[i]).find("-tica-rank=") == 0) {
      msm.tica_rank =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--tica-rank" ||
         string(argv[i]) == "-tica-rank") &&
        i + 1 < argc) {
      msm.tica_rank = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]) == "-umap") {
      umap.flag_umap = true;
      msm.flag_tica = true;
      flag_network = true;
    }
    if (string(argv[i]).find("--umap-neighbors=") == 0 ||
        string(argv[i]).find("-umap-neighbors=") == 0) {
      umap.umap_neighbors =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--umap-neighbors" ||
         string(argv[i]) == "-umap-neighbors") &&
        i + 1 < argc) {
      umap.umap_neighbors = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--umap-components=") == 0 ||
        string(argv[i]).find("-umap-components=") == 0) {
      umap.umap_components =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--umap-components" ||
         string(argv[i]) == "-umap-components") &&
        i + 1 < argc) {
      umap.umap_components = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--umap-epochs=") == 0 ||
        string(argv[i]).find("-umap-epochs=") == 0) {
      umap.umap_epochs =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--umap-epochs" ||
         string(argv[i]) == "-umap-epochs") &&
        i + 1 < argc) {
      umap.umap_epochs = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--umap-min-dist=") == 0 ||
        string(argv[i]).find("-umap-min-dist=") == 0) {
      umap.umap_min_dist =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--umap-min-dist" ||
         string(argv[i]) == "-umap-min-dist") &&
        i + 1 < argc) {
      umap.umap_min_dist = stod(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--umap-spread=") == 0 ||
        string(argv[i]).find("-umap-spread=") == 0) {
      umap.umap_spread =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--umap-spread" ||
         string(argv[i]) == "-umap-spread") &&
        i + 1 < argc) {
      umap.umap_spread = stod(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--umap-seed=") == 0 ||
        string(argv[i]).find("-umap-seed=") == 0) {
      umap.umap_seed = static_cast<unsigned int>(
          stoul(string(argv[i]).substr(string(argv[i]).find('=') + 1)));
    }
    if ((string(argv[i]) == "--umap-seed" ||
         string(argv[i]) == "-umap-seed") &&
        i + 1 < argc) {
      umap.umap_seed = static_cast<unsigned int>(stoul(string(argv[i + 1])));
      ++i;
    }
    if (string(argv[i]).find("--umap-knn-mode=") == 0 ||
        string(argv[i]).find("-umap-knn-mode=") == 0) {
      umap.umap_knn_mode =
          string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if ((string(argv[i]) == "--umap-knn-mode" ||
         string(argv[i]) == "-umap-knn-mode") &&
        i + 1 < argc) {
      umap.umap_knn_mode = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]) == "-umap-knn-benchmark" ||
        string(argv[i]) == "--umap-knn-benchmark") {
      umap.umap_knn_benchmark = true;
    }
    if (string(argv[i]) == "-umap-verbose" ||
        string(argv[i]) == "--umap-verbose") {
      umap.umap_verbose = true;
    }
    // --- MSM Phase 2 flags ---
    if (string(argv[i]) == "-msm") {
      msm.flag_msm = true;
      msm.flag_tica = true;    // MSM depends on TICA projection
      flag_network = true; // TICA depends on CA data from network module
    }
    if (string(argv[i]) == "-msm-clusters" && i + 1 < argc) {
      msm.msm_clusters = stoi(string(argv[i + 1]));
      msm.msm_clusters_user_set = true;
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
    if (string(argv[i]).find("--gist-resolution=") == 0 ||
        string(argv[i]).find("-gist-resolution=") == 0) {
      gist.gist_resolution =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]).find("--gist-padding=") == 0 ||
        string(argv[i]).find("-gist-padding=") == 0) {
      gist.gist_padding =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]).find("--gist-temperature=") == 0 ||
        string(argv[i]).find("-gist-temperature=") == 0) {
      gist.gist_temperature =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]).find("--gist-water-model=") == 0 ||
        string(argv[i]).find("-gist-water-model=") == 0) {
      gist.gist_water_model =
          string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if (string(argv[i]).find("--gist-min-count=") == 0 ||
        string(argv[i]).find("-gist-min-count=") == 0) {
      gist.gist_min_count =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]).find("--gist-cutoff=") == 0 ||
        string(argv[i]).find("-gist-cutoff=") == 0) {
      gist.gist_cutoff =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
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
      salt_concentration = saltcon;
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
    if (string(argv[i]) == "-opdb" || string(argv[i]) == "-opdbsep" ||
        string(argv[i]).rfind("-opdb=", 0) == 0 ||
        string(argv[i]).rfind("-opdbsep=", 0) == 0) {
      opdb = true;
      const string arg(argv[i]);
      if (arg == "-opdbsep" || arg.rfind("-opdbsep=", 0) == 0)
        flag_separation_results = true;

      if (arg.rfind("-opdb=", 0) == 0 || arg.rfind("-opdbsep=", 0) == 0) {
        const size_t eq_pos = arg.find('=');
        if (eq_pos != string::npos && eq_pos + 1 < arg.size()) {
          out_pdb_name = arg.substr(eq_pos + 1);
          opdb_name_explicit = true;
        }
      } else if (i + 1 < argc) {
        if (string(argv[i + 1]).find(".pdb")) {
          out_pdb_name = string(argv[i + 1]);
          opdb_name_explicit = true;
        }
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
    if (string(argv[i]) == "-pbctopo") {
      pbc.solv_pbctopo = true;
    }
    if (string(argv[i]).find("-pbctopo-anchor-backend=") == 0) {
      pbc.pbctopo_anchor_backend =
          string(argv[i]).substr(string(argv[i]).find('=') + 1);
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
    if (arg == "-mol-1-type" && i + 1 < argc) {
      assign_mol_type(string(argv[i + 1]), true);
    }
    if (arg == "-mol-2-type" && i + 1 < argc) {
      assign_mol_type(string(argv[i + 1]), false);
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
      ++i;
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
      string spec = string(argv[i + 1]);
      int consumed_tokens = 1;
      while (i + consumed_tokens + 1 < argc) {
        const string next_token = string(argv[i + consumed_tokens + 1]);
        if (!next_token.empty() && next_token[0] == '-' && next_token != "-" &&
            next_token != "&") {
          break;
        }
        spec += " " + next_token;
        ++consumed_tokens;
      }

      auto sep = spec.find('&');
      if (sep == string::npos)
        sep = spec.find('-');

      auto parse_group_list = [](const string &text, vector<int> &out_groups) {
        regex reg{"(\\d+)"};
        sregex_iterator pos(text.cbegin(), text.cend(), reg), end;
        for (; pos != end; ++pos)
          out_groups.emplace_back(stoi(pos->str(1)));
      };

      vector<int> mol1_groups, mol2_groups;
      if (sep == string::npos) {
        cout << "Your argument -sg is invalid!\n"
             << "Use - or & to separate groups of two molecules!\n";
      } else {
        parse_group_list(spec.substr(0, sep), mol1_groups);
        parse_group_list(spec.substr(sep + 1), mol2_groups);
        if (mol1_groups.empty() || mol2_groups.empty()) {
          cout << "Your argument -sg is invalid!\n"
               << "Mol 1 and Mol 2 subgroup lists must both be non-empty.\n";
        } else {
          sub_group_1 = std::move(mol1_groups);
          sub_group_2 = std::move(mol2_groups);
          flag_sub_groups = true;
        }
      }
      i += consumed_tokens;
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
    if (string(argv[i]).find("-contact-map=") == 0) {
      string val = string(argv[i]).substr(string("-contact-map=").size());
      flag_contact_map_explicit = true;
      if (val == "no")
        flag_cal_contact = false;
      else if (val == "yes")
        flag_cal_contact = true;
    }
    if (string(argv[i]) == "-contact-map" && i + 1 < argc) {
      string val = string(argv[i + 1]);
      flag_contact_map_explicit = true;
      if (val == "no")
        flag_cal_contact = false;
      else if (val == "yes")
        flag_cal_contact = true;
      ++i;
    }
    // --- DHAM (Dynamic Histogram Analysis Method) flags ---
    if (string(argv[i]) == "-dham") {
      dham.flag_dham = true;
      msm.flag_msm = true;      // DHAM depends on MSM (K-means + dtrajs)
      msm.flag_tica = true;     // MSM depends on TICA projection
      flag_network = true;      // TICA depends on CA data from network module
    }
    if (string(argv[i]).find("-dham-temps=") == 0 ||
        string(argv[i]).find("--dham-temps=") == 0) {
      dham.temps_file = string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if (string(argv[i]) == "-dham-temps" && i + 1 < argc) {
      dham.temps_file = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]).find("-dham-ener=") == 0 ||
        string(argv[i]).find("--dham-ener=") == 0) {
      dham.ener_list_file =
          string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if (string(argv[i]) == "-dham-ener" && i + 1 < argc) {
      dham.ener_list_file = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]).find("-dham-target-temp=") == 0 ||
        string(argv[i]).find("--dham-target-temp=") == 0) {
      dham.target_temp =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]) == "-dham-target-temp" && i + 1 < argc) {
      dham.target_temp = stod(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("-dham-ensemble=") == 0 ||
        string(argv[i]).find("--dham-ensemble=") == 0) {
      dham.ensemble = string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if (string(argv[i]) == "-dham-ensemble" && i + 1 < argc) {
      dham.ensemble = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]).find("-dham-window-match=") == 0 ||
        string(argv[i]).find("--dham-window-match=") == 0) {
      dham.window_match =
          string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if (string(argv[i]) == "-dham-window-match" && i + 1 < argc) {
      dham.window_match = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]).find("-dham-pressure-bar=") == 0 ||
        string(argv[i]).find("--dham-pressure-bar=") == 0) {
      dham.target_pressure_bar =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]) == "-dham-pressure-bar" && i + 1 < argc) {
      dham.target_pressure_bar = stod(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--dham-max-iter=") == 0 ||
        string(argv[i]).find("-dham-max-iter=") == 0) {
      dham.max_iter =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]).find("--dham-tol=") == 0 ||
        string(argv[i]).find("-dham-tol=") == 0) {
      dham.tolerance =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]).find("--dham-n-replicas=") == 0 ||
        string(argv[i]).find("-dham-n-replicas=") == 0) {
      dham.n_replicas =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]).find("-dham-dt-energy=") == 0 ||
        string(argv[i]).find("--dham-dt-energy=") == 0) {
      dham.dt_energy_ps =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if (string(argv[i]) == "-dham-dt-energy" && i + 1 < argc) {
      dham.dt_energy_ps = stod(string(argv[i + 1]));
      ++i;
    }
    // --- TPT (Transition Path Theory) flags ---
    if (string(argv[i]) == "-tpt") {
      tpt.flag_tpt = true;
      msm.flag_msm = true;  // TPT depends on MSM microstates / TPM
      msm.flag_tica = true; // MSM depends on TICA projection
      flag_network = true;  // protein workflows still use network -> TICA -> MSM
    }
    if (string(argv[i]).find("--tpt-mode=") == 0 ||
        string(argv[i]).find("-tpt-mode=") == 0) {
      tpt.mode = string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if ((string(argv[i]) == "--tpt-mode" || string(argv[i]) == "-tpt-mode") &&
        i + 1 < argc) {
      tpt.mode = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]).find("--tpt-source=") == 0 ||
        string(argv[i]).find("-tpt-source=") == 0) {
      tpt.source_spec = string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if ((string(argv[i]) == "--tpt-source" ||
         string(argv[i]) == "-tpt-source") &&
        i + 1 < argc) {
      tpt.source_spec = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]).find("--tpt-sink=") == 0 ||
        string(argv[i]).find("-tpt-sink=") == 0) {
      tpt.sink_spec = string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if ((string(argv[i]) == "--tpt-sink" || string(argv[i]) == "-tpt-sink") &&
        i + 1 < argc) {
      tpt.sink_spec = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]).find("--tpt-top-paths=") == 0 ||
        string(argv[i]).find("-tpt-top-paths=") == 0) {
      tpt.top_paths =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--tpt-top-paths" ||
         string(argv[i]) == "-tpt-top-paths") &&
        i + 1 < argc) {
      tpt.top_paths = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--tpt-flux-cut=") == 0 ||
        string(argv[i]).find("-tpt-flux-cut=") == 0) {
      tpt.flux_cut =
          stod(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--tpt-flux-cut" ||
         string(argv[i]) == "-tpt-flux-cut") &&
        i + 1 < argc) {
      tpt.flux_cut = stod(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--tpt-use-dham=") == 0 ||
        string(argv[i]).find("-tpt-use-dham=") == 0) {
      std::string value =
          string(argv[i]).substr(string(argv[i]).find('=') + 1);
      std::transform(value.begin(), value.end(), value.begin(),
                     [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                     });
      tpt.use_dham =
          !(value == "no" || value == "false" || value == "0" ||
            value == "off");
    }
    if (string(argv[i]).find("--tpt-demux-map=") == 0 ||
        string(argv[i]).find("-tpt-demux-map=") == 0) {
      tpt.demux_map_file =
          string(argv[i]).substr(string(argv[i]).find('=') + 1);
    }
    if ((string(argv[i]) == "--tpt-demux-map" ||
         string(argv[i]) == "-tpt-demux-map") &&
        i + 1 < argc) {
      tpt.demux_map_file = string(argv[i + 1]);
      ++i;
    }
    if (string(argv[i]).find("--tpt-bootstrap=") == 0 ||
        string(argv[i]).find("-tpt-bootstrap=") == 0) {
      tpt.bootstrap_samples =
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
    }
    if ((string(argv[i]) == "--tpt-bootstrap" ||
         string(argv[i]) == "-tpt-bootstrap") &&
        i + 1 < argc) {
      tpt.bootstrap_samples = stoi(string(argv[i + 1]));
      ++i;
    }
    if (string(argv[i]).find("--tpt-bootstrap-seed=") == 0 ||
        string(argv[i]).find("-tpt-bootstrap-seed=") == 0) {
      tpt.bootstrap_seed = static_cast<unsigned int>(
          stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1)));
    }
    if ((string(argv[i]) == "--tpt-bootstrap-seed" ||
         string(argv[i]) == "-tpt-bootstrap-seed") &&
        i + 1 < argc) {
      tpt.bootstrap_seed = static_cast<unsigned int>(stoi(string(argv[i + 1])));
      ++i;
    }
    if ((string(argv[i]) == "--tpt-use-dham" ||
         string(argv[i]) == "-tpt-use-dham") &&
        i + 1 < argc) {
      std::string value = string(argv[i + 1]);
      std::transform(value.begin(), value.end(), value.begin(),
                     [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                     });
      tpt.use_dham =
          !(value == "no" || value == "false" || value == "0" ||
            value == "off");
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

  if (pbc.solv_pbctopo) {
    const bool conflict_with_reference =
        pbc.solv_pbc == 1 || pbc.solv_pbc_int == 1 || pbc.solv_pbcn ||
        pbc.pbc_aggressive || pbc.pbc_mem || pbc.flag_pbc_ref ||
        pbc.pbc_envi || pbc.solv_pbc_clust;
    if (conflict_with_reference) {
      cout << "ERROR!! The option -pbctopo can not be combined with other "
              "reference/legacy PBC modes.\n";
      return 2;
    }
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

      if (flag == "-pbctopo-anchor-backend" || flag == "-community") {
        tokens.push_back(flag + "=" + value);
      } else {
        // Split flag and value into separate tokens (single-dash only in arg.in)
        tokens.push_back(flag);
        if (!value.empty()) {
          std::istringstream iss(value);
          std::string word;
          while (iss >> word)
            tokens.push_back(word);
        }
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

  namespace fs = std::filesystem;

  // Resolve relative file/directory arguments against the location of arg.in
  // before parse_arguments(). Some flags (e.g. -topo) validate paths
  // immediately while parsing, so post-parse fixups are too late.
  std::string arg_dir;
  if (arg_filename.find('\\') != std::string::npos ||
      arg_filename.find('/') != std::string::npos) {
    arg_dir = arg_filename.substr(0, arg_filename.find_last_of("\\/"));
  }
  if (!arg_dir.empty()) {
    const fs::path arg_dir_path(arg_dir);
    const auto is_path_flag = [](const std::string &flag) {
      return flag == "-topo" || flag == "-topo2" || flag == "-common_topo" ||
             flag == "-common-topo" || flag == "-f" || flag == "-s" ||
             flag == "-n" || flag == "-i" || flag == "-traj" ||
             flag == "-dham-temps" || flag == "-dham-ener" ||
             flag == "-path-a" || flag == "-path-n" || flag == "-pbc-ref";
    };

    for (size_t i = 0; i + 1 < tokens.size(); ++i) {
      if (!is_path_flag(tokens[i]))
        continue;
      if (tokens[i + 1].empty() || tokens[i + 1] == "no")
        continue;

      fs::path candidate(tokens[i + 1]);
      if (candidate.is_relative()) {
        tokens[i + 1] = (arg_dir_path / candidate).lexically_normal().string();
      }
      ++i;
    }
  }

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

  // Resolve relative paths for file arguments relative to arg file dir
  if (!arg_dir.empty()) {
    const fs::path arg_dir_path(arg_dir);
    auto resolve = [&](std::string &path) {
      if (!path.empty() && path != "no") {
        fs::path candidate(path);
        if (candidate.is_relative()) {
          path = (arg_dir_path / candidate).lexically_normal().string();
        }
      }
    };
    if (flag_topo_1)
      resolve(topo_name);
    if (flag_topo_2)
      resolve(topo_name_2);
    if (flag_1_pdb)
      resolve(pdb_name);
    if (flag_index)
      resolve(index_name);
    if (flag_listfile)
      resolve(list_filename);
    if (flag_xdr)
      resolve(traj_file);
    if (pbc.flag_pbc_ref)
      resolve(pbc.ref_pdb_name);
    resolve(path_amino_acid);
    resolve(path_nonbonded);
    if (dham.flag_dham) {
      resolve(dham.temps_file);
      resolve(dham.ener_list_file);
    }
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
    auto assign_topo_name = [&](const std::string &name) {
      if (Topo == 2)
        topo_name_2 = name;
      else
        topo_name = name;
    };

    fstream input_topo_protein(NameofTopoFile, ios::in);
    bool topo_is_open{input_topo_protein.is_open()};
    input_topo_protein.close();
    if (topo_is_open == true)
      return true;
    else {
      int attempts = 0;
      constexpr int MAX_ATTEMPTS = 5;
      while (!topo_is_open) {
        cout << "Cannot open " << NameofTopoFile
             << ". Please type your lig topo file (no if you don't want to use "
                 "topo): ";
        cin >> NameofTopoFile;
        cout << NameofTopoFile << "\n";
        assign_topo_name(NameofTopoFile);
        if (NameofTopoFile == "no")
          return false;
        else {
          input_topo_protein.open(NameofTopoFile, ios::in);
          if (input_topo_protein.is_open()) {
            input_topo_protein.close();
            return true;
          } else {
            input_topo_protein.close();
            ++attempts;
            if (attempts >= MAX_ATTEMPTS) {
              TitanUI::print_error("Too many invalid attempts (" +
                                   std::to_string(MAX_ATTEMPTS) +
                                   "). Exiting.");
              return false;
            }
            cout << "File not found. Attempts remaining: "
                 << (MAX_ATTEMPTS - attempts) << '\n';
            continue;
          }
        }
      }
      return false;
    }
  }
}
