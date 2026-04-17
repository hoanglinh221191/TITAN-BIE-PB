#include "arguments.h"
#include "titan_tpt.h"
#include "titan_ui.h"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <thread>

namespace fs = std::filesystem;

bool arguments::parse_input_file_is_open() {
  bool check_file_open{false};
  fstream test_file(pdb_name, ios::in);
  if (test_file.is_open() || flag_all_folders == true)
    check_file_open = true;
  test_file.close();
  return check_file_open;
}

int arguments::recheck_and_change_pdb_input() {
  add_input_file_type();

  if (flag_1_pdb == true && parse_input_file_is_open() == false &&
      flag_pattern_input == false && flag_listfile == false) {
    bool check_file_open{false};
    int attempts = 0;
    constexpr int MAX_ATTEMPTS = 5;
    do {
      pdb_name = string();
      cout << "Please type your input " << get_input_file_type() << " file:\n";
      cin >> pdb_name;
      add_input_file_type();
      cout << "Your input: " << pdb_name << '\n';
      if (pdb_name.length() < 1)
        return 1;
      flag_1_pdb = true;
      check_file_open = parse_input_file_is_open();
      if (!check_file_open) {
        ++attempts;
        if (attempts >= MAX_ATTEMPTS) {
          TitanUI::print_error("Too many invalid attempts (" +
                               std::to_string(MAX_ATTEMPTS) + "). Exiting.");
          return 1;
        }
        cout << "File not found. Attempts remaining: "
             << (MAX_ATTEMPTS - attempts) << '\n';
      }
    } while (check_file_open == false);
  } else
    return 0;
  return 0;
}

string arguments::set_Buildnumber() { return "20260406"; }
string arguments::get_build_number() const { return "20260406"; }

string arguments::get_pdbOutput() const {
  if (opdb_name_explicit && !out_pdb_name.empty())
    return out_pdb_name;
  fs::path p_temp = pdb_name;
  return "out_" + p_temp.stem().string() + ".pdb";
}

int arguments::post_process() {
  recheck_and_change_pdb_input();

  if (pdb_name.length() > 4 &&
      pdb_name.substr(pdb_name.length() - 4) != get_input_file_type()) {
    cout << "Your input file " << pdb_name << " is invalid, input file must be "
         << get_input_file_type() << " file!\n";
    exit(EXIT_FAILURE);
  }
  if (force_f != "AMBER99SB" && force_f != "AMBER99SB-ILDN" &&
      force_f != "AMBER99SB-BSC1" && force_f != "AMBER99SB-DISP" &&
      force_f != "AMBER14SB" && force_f != "AMBER19SB" &&
      force_f != "CHARMM27" && force_f != "CHARMM36" &&
      force_f != "GROMOS43A1" && force_f != "GROMOS53A6" &&
      force_f != "GROMOS54A7" && force_f != "OPLSAA" && force_f != "MARTINI") {
    cout << "Force field " << force_f << " is not supported!\n"
         << "This tool supports AMBER99SB, AMBER99SB-ILDN, AMBER99SB-BSC1, "
            "AMBER99SB-DISP, AMBER14SB, AMBER19SB, GROMOS43A1, GROMOS53A6, "
            "GROMOS54A7, CHARMM27, CHARMM36, OPLSAA and MARTINI\n"
            "Please check your option.\n";
    exit(EXIT_FAILURE);
  }
  // Martini CG uses relative permittivity epsilon_r = 15
  if (force_f == "MARTINI" && epsilon_in == 1.0) {
    epsilon_in = 15.0;
    TitanUI::print_info(
        "[MARTINI] Auto-set epsilon_in = 15.0 (Martini relative "
        "permittivity). Override with -epsilon_in if needed.");
  }
  if (flag_elec_screening) {
    if (salt_concentration <= 0.0) {
      TitanUI::print_error(
          "Electrostatic screening: salt concentration must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (epsilon_in <= 0.0) {
      TitanUI::print_error(
          "Electrostatic screening: -epsilon_in must be > 0.");
      exit(EXIT_FAILURE);
    }
  }
  if (begin_time < 0. || end_time < 0.) {
    cout << "Your begin time and end time must not be negative!\n";
    exit(EXIT_FAILURE);
  }
  if (dt < 0) {
    cout << "Time interval must not be negative!\n";
    exit(EXIT_FAILURE);
  }
  if (mbar.flag_mbar) {
    std::transform(mbar.domain.begin(), mbar.domain.end(), mbar.domain.begin(),
                   [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    std::transform(mbar.ensemble.begin(), mbar.ensemble.end(),
                   mbar.ensemble.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    std::transform(mbar.window_match.begin(), mbar.window_match.end(),
                   mbar.window_match.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    std::transform(mbar.temp_mode.begin(), mbar.temp_mode.end(),
                   mbar.temp_mode.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    flag_network = true;
    if (mbar.domain == "tica" || mbar.domain == "both") {
      msm.flag_tica = true;
      msm.flag_msm = true;
    }
    if (mbar.domain == "pca" || mbar.domain == "both") {
      flag_pca = true;
    }
  }
  if ((msm.flag_msm || dham.flag_dham ||
       (mbar.flag_mbar && (mbar.domain == "tica" || mbar.domain == "both"))) &&
      msm.msm_lag <= 0) {
    TitanUI::print_error("MSM/DHAM/MBAR: -msm-lag must be > 0.");
    exit(EXIT_FAILURE);
  }
  if (dham.flag_dham) {
    std::transform(dham.ensemble.begin(), dham.ensemble.end(),
                   dham.ensemble.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    std::transform(dham.window_match.begin(), dham.window_match.end(),
                   dham.window_match.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (dham.ensemble != "nvt" && dham.ensemble != "npt") {
      TitanUI::print_error("DHAM: -dham-ensemble must be either nvt or npt.");
      exit(EXIT_FAILURE);
    }
    if (dham.window_match != "auto" && dham.window_match != "off" &&
        dham.window_match != "strict") {
      TitanUI::print_error(
          "DHAM: -dham-window-match must be auto, off, or strict.");
      exit(EXIT_FAILURE);
    }
    if (dham.target_temp <= 0.0) {
      TitanUI::print_error("DHAM: -dham-target-temp must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (dham.target_pressure_bar < 0.0) {
      TitanUI::print_error("DHAM: -dham-pressure-bar must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (dham.max_iter <= 0) {
      TitanUI::print_error("DHAM: --dham-max-iter must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (dham.tolerance <= 0.0) {
      TitanUI::print_error("DHAM: --dham-tol must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (dham.dt_energy_ps < 0.0) {
      TitanUI::print_error("DHAM: -dham-dt-energy must be >= 0.");
      exit(EXIT_FAILURE);
    }
  }
  if (mbar.flag_mbar) {
    if (mbar.domain != "tica" && mbar.domain != "pca" &&
        mbar.domain != "both") {
      TitanUI::print_error("MBAR: -mbar-domain must be tica, pca, or both.");
      exit(EXIT_FAILURE);
    }
    if (mbar.ensemble != "nvt" && mbar.ensemble != "npt") {
      TitanUI::print_error("MBAR: -mbar-ensemble must be either nvt or npt.");
      exit(EXIT_FAILURE);
    }
    if (mbar.window_match != "auto" && mbar.window_match != "off" &&
        mbar.window_match != "strict") {
      TitanUI::print_error(
          "MBAR: -mbar-window-match must be auto, off, or strict.");
      exit(EXIT_FAILURE);
    }
    if (mbar.temp_mode != "fixed" && mbar.temp_mode != "per_frame") {
      TitanUI::print_error(
          "MBAR: -mbar-temp-mode must be fixed or per_frame.");
      exit(EXIT_FAILURE);
    }
    if (mbar.temps_file.empty() && mbar.temp_mode != "per_frame") {
      TitanUI::print_error(
          "MBAR: omitting -mbar-temps is only supported with "
          "-mbar-temp-mode per_frame.");
      exit(EXIT_FAILURE);
    }
    if (mbar.target_temp <= 0.0) {
      TitanUI::print_error("MBAR: -mbar-target-temp must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (mbar.target_pressure_bar < 0.0) {
      TitanUI::print_error("MBAR: -mbar-pressure-bar must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (mbar.dt_energy_ps < 0.0) {
      TitanUI::print_error("MBAR: -mbar-dt-energy must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (mbar.temp_column < 1) {
      TitanUI::print_error("MBAR: -mbar-temp-col must be >= 1.");
      exit(EXIT_FAILURE);
    }
    if (mbar.temp_mode == "per_frame" && mbar.demux_map_file.empty() &&
        mbar.temp_column < 2) {
      TitanUI::print_error(
          "MBAR: raw per-frame mode expects a 3-column XVG layout "
          "(time, observable, temperature), so -mbar-temp-col must be >= 2.");
      exit(EXIT_FAILURE);
    }
    if (mbar.bootstrap_samples < 0) {
      TitanUI::print_error("MBAR: -mbar-bootstrap must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (mbar.scan_points < 0) {
      TitanUI::print_error("MBAR: -mbar-scan-points must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (!mbar.demux_map_file.empty() && mbar.temp_mode != "per_frame") {
      TitanUI::print_info(
          "MBAR: -mbar-demux-map is only used with -mbar-temp-mode per_frame. "
          "Fixed-temperature mode will ignore it.");
    }
  }
  if (tpt.flag_tpt) {
    std::transform(tpt.mode.begin(), tpt.mode.end(), tpt.mode.begin(),
                   [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (tpt.mode != "ensemble" && tpt.mode != "kinetic") {
      TitanUI::print_error(
          "TPT: -tpt-mode must be either ensemble or kinetic.");
      exit(EXIT_FAILURE);
    }
    if (tpt.source_spec.empty()) {
      TitanUI::print_error(
          "TPT: -tpt-source is required when -tpt is enabled.");
      exit(EXIT_FAILURE);
    }
    if (tpt.sink_spec.empty()) {
      TitanUI::print_error(
          "TPT: -tpt-sink is required when -tpt is enabled.");
      exit(EXIT_FAILURE);
    }
    if (tpt.top_paths <= 0) {
      TitanUI::print_error("TPT: -tpt-top-paths must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (tpt.flux_cut < 0.0) {
      TitanUI::print_error("TPT: -tpt-flux-cut must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (tpt.bootstrap_samples < 0) {
      TitanUI::print_error("TPT: -tpt-bootstrap must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (tpt.use_dham && !dham.flag_dham) {
      TitanUI::print_info(
          "TPT: -tpt-use-dham requested without -dham. TPT will fall back to MSM.");
    }
    if (!tpt.demux_map_file.empty() && tpt.mode != "kinetic") {
      TitanUI::print_info(
          "TPT: -tpt-demux-map is only used in kinetic mode. Ensemble TPT will ignore it.");
    }

    auto basin_prefix = [](const std::string &spec) {
      const auto colon_pos = spec.find(':');
      if (colon_pos == std::string::npos)
        return std::string("micro");
      std::string prefix = spec.substr(0, colon_pos);
      std::transform(prefix.begin(), prefix.end(), prefix.begin(),
                     [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                     });
      return prefix;
    };
    const std::string source_prefix = basin_prefix(tpt.source_spec);
    const std::string sink_prefix = basin_prefix(tpt.sink_spec);
    const auto prefix_supported = [](const std::string &prefix) {
      return prefix == "micro" || prefix == "state" || prefix == "states" ||
             prefix == "macro" || prefix == "tica_rect" || prefix == "lj" ||
             prefix == "ljpreset" || prefix == "lj_preset" ||
             prefix == "protein" || prefix == "prot" ||
             prefix == "protein_preset";
    };
    if (!prefix_supported(source_prefix)) {
      TitanUI::print_error("TPT: unsupported -tpt-source prefix '" +
                           source_prefix + "'.");
      exit(EXIT_FAILURE);
    }
    if (!prefix_supported(sink_prefix)) {
      TitanUI::print_error("TPT: unsupported -tpt-sink prefix '" +
                           sink_prefix + "'.");
      exit(EXIT_FAILURE);
    }
    if (source_prefix == "micro" || source_prefix == "state" ||
        source_prefix == "states") {
      std::vector<int> source_states;
      std::string parse_error;
      if (!TitanTPT::parse_state_spec(tpt.source_spec, -1, source_states,
                                      &parse_error)) {
        TitanUI::print_error("TPT: invalid -tpt-source: " + parse_error);
        exit(EXIT_FAILURE);
      }
    }
    if (sink_prefix == "micro" || sink_prefix == "state" ||
        sink_prefix == "states") {
      std::vector<int> sink_states;
      std::string parse_error;
      if (!TitanTPT::parse_state_spec(tpt.sink_spec, -1, sink_states,
                                      &parse_error)) {
        TitanUI::print_error("TPT: invalid -tpt-sink: " + parse_error);
        exit(EXIT_FAILURE);
      }
    }
    if ((source_prefix == "micro" || source_prefix == "state" ||
         source_prefix == "states") &&
        (sink_prefix == "micro" || sink_prefix == "state" ||
         sink_prefix == "states")) {
      std::vector<int> source_states;
      std::vector<int> sink_states;
      std::string parse_error;
      TitanTPT::parse_state_spec(tpt.source_spec, -1, source_states,
                                 &parse_error);
      TitanTPT::parse_state_spec(tpt.sink_spec, -1, sink_states,
                                 &parse_error);
      std::vector<unsigned char> basin_mask;
      int max_state = -1;
      for (int state : source_states)
        max_state = std::max(max_state, state);
      for (int state : sink_states)
        max_state = std::max(max_state, state);
      basin_mask.assign(static_cast<size_t>(max_state + 1), 0);
      for (int state : source_states)
        basin_mask[static_cast<size_t>(state)] = 1;
      for (int state : sink_states) {
        if (basin_mask[static_cast<size_t>(state)] != 0) {
          TitanUI::print_error(
              "TPT: source and sink basins must be disjoint.");
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if (awh.flag_awh_ifp) {
    if (!ifp.flag_ifp) {
      TitanUI::print_error(
          "AWH: -awh-ifp requires IFP to be enabled.");
      exit(EXIT_FAILURE);
    }
    if (awh.xvg_list_file.empty()) {
      TitanUI::print_error(
          "AWH: -awh-xvg <filelist> is required when -awh-ifp is enabled.");
      exit(EXIT_FAILURE);
    }
    if (awh.xvg_value_column < 1) {
      TitanUI::print_error("AWH: -awh-xvg-col must be >= 1.");
      exit(EXIT_FAILURE);
    }
    if (awh.n_bins <= 0) {
      TitanUI::print_error("AWH: -awh-bins must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (flag_nocal) {
      TitanUI::print_error(
          "AWH: -awh-ifp currently requires two-group analysis and does not support -nocal.");
      exit(EXIT_FAILURE);
    }
    if (awh.output_prefix.empty()) {
      awh.output_prefix = "awh";
    }
  }

  if (community.flag_community) {
    std::transform(community.method.begin(), community.method.end(),
                   community.method.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (community.method != "louvain") {
      TitanUI::print_error("Community: unsupported method '" +
                           community.method +
                           "'. V1 supports only -community=louvain.");
      exit(EXIT_FAILURE);
    }
  }
  if (msm.flag_msm || dham.flag_dham) {
    if (msm.msm_nts <= 0) {
      TitanUI::print_error("MSM: -msm-nts must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (msm.msm_its_topk < 0) {
      TitanUI::print_error("MSM: -msm-its-topk must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (msm.msm_its_max_k_full <= 0) {
      TitanUI::print_error("MSM: -msm-its-max-k-full must be > 0.");
      exit(EXIT_FAILURE);
    }
  }
  if (flag_pca) {
    if (msm.pca_topk < 0) {
      TitanUI::print_error("PCA: -pca-topk must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (msm.pca_max_k_full <= 0) {
      TitanUI::print_error("PCA: -pca-max-k-full must be > 0.");
      exit(EXIT_FAILURE);
    }
  }
  if (msm.flag_tica) {
    if (msm.tica_lag <= 0) {
      TitanUI::print_error("TICA: -tica-lag must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (msm.tica_components <= 0) {
      TitanUI::print_error("TICA: -tica-components must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (msm.tica_rank <= 0) {
      TitanUI::print_error("TICA: -tica-rank must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (msm.tica_rank < msm.tica_components) {
      TitanUI::print_error("TICA: -tica-rank must be >= -tica-components.");
      exit(EXIT_FAILURE);
    }
    std::string backend = msm.tica_backend;
    std::transform(backend.begin(), backend.end(), backend.begin(),
                   [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (backend != "auto" && backend != "exact_dense" && backend != "dense" &&
        backend != "exact_dual" && backend != "dual" &&
        backend != "randomized" && backend != "rand") {
      TitanUI::print_error(
          "TICA: -tica-backend must be one of auto, exact_dense, exact_dual, randomized.");
      exit(EXIT_FAILURE);
    }
    if (backend == "dense") {
      backend = "exact_dense";
    } else if (backend == "dual") {
      backend = "exact_dual";
    } else if (backend == "rand") {
      backend = "randomized";
    }
    const bool need_tica_projection_cache =
        msm.flag_msm || dham.flag_dham || umap.flag_umap ||
        (zts.flag_zts && zts.zts_domain == "tica");
    if (need_tica_projection_cache && backend == "exact_dense") {
      TitanUI::print_warning(
          "TICA: exact_dense does not retain per-frame projection history "
          "needed by MSM/DHAM/UMAP/ZTS. Switching backend to exact_dual.");
      backend = "exact_dual";
    }
    msm.tica_backend = backend;
  }
  if (zts.flag_zts) {
    std::transform(zts.zts_domain.begin(), zts.zts_domain.end(),
                   zts.zts_domain.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    std::transform(zts.zts_init_mode.begin(), zts.zts_init_mode.end(),
                   zts.zts_init_mode.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (flag_pca && msm.flag_tica && !zts.zts_domain_explicit) {
      TitanUI::print_error(
          "ZTS: when both -pca and -tica are enabled, -zts-domain is required.");
      exit(EXIT_FAILURE);
    }
    if (zts.zts_domain != "pca" && zts.zts_domain != "tica") {
      TitanUI::print_error("ZTS: -zts-domain must be either pca or tica.");
      exit(EXIT_FAILURE);
    }
    if (zts.zts_domain == "pca" && !flag_pca) {
      TitanUI::print_error("ZTS: -zts-domain pca requires -pca.");
      exit(EXIT_FAILURE);
    }
    if (zts.zts_domain == "tica" && !msm.flag_tica) {
      TitanUI::print_error("ZTS: -zts-domain tica requires -tica.");
      exit(EXIT_FAILURE);
    }
    if (zts.zts_images < 3) {
      TitanUI::print_error("ZTS: -zts-images must be >= 3.");
      exit(EXIT_FAILURE);
    }
    if (zts.zts_step <= 0.0) {
      TitanUI::print_error("ZTS: -zts-step must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (zts.zts_tol <= 0.0) {
      TitanUI::print_error("ZTS: -zts-tol must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (zts.zts_max_iter <= 0) {
      TitanUI::print_error("ZTS: -zts-iter must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (zts.zts_fes_sigma <= 0.0) {
      TitanUI::print_error("ZTS: -zts-sigma must be > 0.");
      exit(EXIT_FAILURE);
    }
    const bool source_explicit = zts.zts_source >= 0;
    const bool sink_explicit = zts.zts_sink >= 0;
    if (source_explicit != sink_explicit) {
      TitanUI::print_error(
          "ZTS: -zts-source and -zts-sink must be set together, or both left at auto.");
      exit(EXIT_FAILURE);
    }
    if (source_explicit && zts.zts_source == zts.zts_sink) {
      TitanUI::print_error(
          "ZTS: -zts-source and -zts-sink must refer to different basins.");
      exit(EXIT_FAILURE);
    }
    if (zts.zts_init_mode != "linear") {
      TitanUI::print_error(
          "ZTS: -zts-init currently supports only 'linear' in this build.");
      exit(EXIT_FAILURE);
    }
  }
  if (umap.flag_umap) {
    std::transform(umap.umap_knn_mode.begin(), umap.umap_knn_mode.end(),
                   umap.umap_knn_mode.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (umap.umap_neighbors < 2) {
      TitanUI::print_error("UMAP: -umap-neighbors must be >= 2.");
      exit(EXIT_FAILURE);
    }
    if (umap.umap_components != 2 && umap.umap_components != 3) {
      TitanUI::print_error("UMAP: -umap-components must be 2 or 3.");
      exit(EXIT_FAILURE);
    }
    if (umap.umap_epochs <= 0) {
      TitanUI::print_error("UMAP: -umap-epochs must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (umap.umap_min_dist <= 0.0) {
      TitanUI::print_error("UMAP: -umap-min-dist must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (umap.umap_spread <= 0.0) {
      TitanUI::print_error("UMAP: -umap-spread must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (umap.umap_knn_mode != "exact" && umap.umap_knn_mode != "approx" &&
        umap.umap_knn_mode != "auto") {
      TitanUI::print_error(
          "UMAP: -umap-knn-mode must be one of exact, approx, auto.");
      exit(EXIT_FAILURE);
    }
  }
  if (group_1 < 0 || group_2 < 0) {
    cout << "Groups for calculating energy are invalid! ";
    cout << "Your index group(s) must not be negative!\n";
    exit(EXIT_FAILURE);
  }
  if (de.flag_find_pulling_direct) {
    std::transform(de.pull_opt.begin(), de.pull_opt.end(), de.pull_opt.begin(),
                   [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (de.pull_opt.empty())
      de.pull_opt = "de";
    if (de.pull_opt == "de-pattern" || de.pull_opt == "de_pattern")
      de.pull_opt = "de+pattern";
    if (de.pull_opt == "de-jso-pattern" || de.pull_opt == "de_jso_pattern")
      de.pull_opt = "de-jso+pattern";
    if (de.pull_opt == "de-nl-shade-lbc-pattern" ||
        de.pull_opt == "de_nl_shade_lbc_pattern")
      de.pull_opt = "de-nl-shade-lbc+pattern";
    if (de.pull_opt == "direct-pattern" || de.pull_opt == "direct_pattern")
      de.pull_opt = "direct+pattern";
    if (de.pull_opt != "de" && de.pull_opt != "de+pattern" &&
        de.pull_opt != "de-shade" && de.pull_opt != "de-l-shade" &&
        de.pull_opt != "de-jso" && de.pull_opt != "de-jso+pattern" &&
        de.pull_opt != "de-nl-shade-lbc" &&
        de.pull_opt != "de-nl-shade-lbc+pattern" && de.pull_opt != "direct" &&
        de.pull_opt != "direct+pattern") {
      TitanUI::print_error(
          "Pull-direction optimizer: supported values are de, de+pattern, de-shade, de-l-shade, de-jso, de-jso+pattern, de-nl-shade-lbc, de-nl-shade-lbc+pattern, direct, and direct+pattern.");
      exit(EXIT_FAILURE);
    }
    const bool use_legacy_de =
        (de.pull_opt == "de" || de.pull_opt == "de+pattern");
    const bool use_shade_family =
        (de.pull_opt == "de-shade" || de.pull_opt == "de-l-shade" ||
         de.pull_opt == "de-jso" || de.pull_opt == "de-jso+pattern" ||
         de.pull_opt == "de-nl-shade-lbc" ||
         de.pull_opt == "de-nl-shade-lbc+pattern");
    const bool use_de_family = use_legacy_de || use_shade_family;
    const bool use_direct =
        (de.pull_opt == "direct" || de.pull_opt == "direct+pattern");
    if ((de.pull_opt == "de-jso" || de.pull_opt == "de-jso+pattern") &&
        !de.user_set_de_p_best_rate)
      de.de_p_best_rate = 0.11;
    if (use_de_family && de.de_pop_num < 4) {
      TitanUI::print_error(
          "Pull-direction DE: -de-pop-num must be >= 4.");
      exit(EXIT_FAILURE);
    }
    if (use_de_family && de.de_ite == 0) {
      TitanUI::print_error(
          "Pull-direction DE: -de-ite-num must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (de.de_cross_prob < 0.0 || de.de_cross_prob > 1.0) {
      TitanUI::print_error(
          "Pull-direction DE: -de-cross-prob must be in [0, 1].");
      exit(EXIT_FAILURE);
    }
    if (de.de_scan_radius <= 0.0 || de.de_min_scan_radius <= 0.0) {
      TitanUI::print_error(
          "Pull-direction DE: scan radii must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (de.de_scan_radius < de.de_min_scan_radius) {
      TitanUI::print_error(
          "Pull-direction DE: -de-scan-radius must be >= -de-min-scan-radius.");
      exit(EXIT_FAILURE);
    }
    if (use_de_family && de.de_stagnation_limit < 1) {
      TitanUI::print_error(
          "Pull-direction DE: -de-stagnation-limit must be >= 1.");
      exit(EXIT_FAILURE);
    }
    if (use_shade_family && de.de_memory_size == 0) {
      TitanUI::print_error(
          "Pull-direction SHADE-family: -de-memory-size must be >= 1.");
      exit(EXIT_FAILURE);
    }
    if (use_shade_family &&
        (de.de_p_best_rate <= 0.0 || de.de_p_best_rate > 1.0)) {
      TitanUI::print_error(
          "Pull-direction SHADE-family: -de-p-best-rate must be in (0, 1].");
      exit(EXIT_FAILURE);
    }
    if (use_shade_family && de.de_archive_scale < 0.0) {
      TitanUI::print_error(
          "Pull-direction SHADE-family: -de-archive-scale must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if ((de.pull_opt == "de-l-shade" || de.pull_opt == "de-jso" ||
         de.pull_opt == "de-jso+pattern" ||
         de.pull_opt == "de-nl-shade-lbc" ||
         de.pull_opt == "de-nl-shade-lbc+pattern") &&
        de.de_min_pop_num < 4) {
      TitanUI::print_error(
          "Pull-direction SHADE-family: -de-min-pop-num must be >= 4.");
      exit(EXIT_FAILURE);
    }
    if ((de.pull_opt == "de-l-shade" || de.pull_opt == "de-jso" ||
         de.pull_opt == "de-jso+pattern" ||
         de.pull_opt == "de-nl-shade-lbc" ||
         de.pull_opt == "de-nl-shade-lbc+pattern") &&
        de.de_min_pop_num > de.de_pop_num) {
      TitanUI::print_error(
          "Pull-direction SHADE-family: -de-min-pop-num must be <= -de-pop-num.");
      exit(EXIT_FAILURE);
    }
    if (use_de_family && de.pull_max_evals != 0 &&
        de.pull_max_evals < de.de_pop_num) {
      TitanUI::print_error(
          "Pull-direction DE: -pull-max-evals must be 0 or >= -de-pop-num.");
      exit(EXIT_FAILURE);
    }
    if (use_direct && de.pull_max_evals == 0 && de.de_ite == 0) {
      TitanUI::print_error(
          "Pull-direction DIRECT: provide -pull-max-evals or keep -de-ite-num > 0 for fallback budget.");
      exit(EXIT_FAILURE);
    }
    if (de.flag_pathway && de.dis_along_path <= 0.0) {
      TitanUI::print_error(
          "Pull-direction DE pathway: -de-pathway-radius must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (use_shade_family &&
        (de.user_set_de_diff_weight || de.user_set_de_cross_prob ||
         de.user_set_de_F || de.user_set_de_algorithm ||
         de.user_set_de_adaptive_F)) {
      TitanUI::print_error(
          "Pull-direction SHADE-family: legacy flags -de-diff-weight, -de-cross-prob, -de-F, -de-algorithm, and -de-adaptive-F are not supported.");
      exit(EXIT_FAILURE);
    }
    if ((use_legacy_de || use_direct) &&
        (de.user_set_de_memory_size || de.user_set_de_p_best_rate ||
         de.user_set_de_archive_scale || de.user_set_de_min_pop_num)) {
      TitanUI::print_error(
          "Pull-direction optimizer: SHADE-family flags -de-memory-size, -de-p-best-rate, -de-archive-scale, and -de-min-pop-num are only valid with de-shade, de-l-shade, de-jso, de-jso+pattern, de-nl-shade-lbc, or de-nl-shade-lbc+pattern.");
      exit(EXIT_FAILURE);
    }
  }
  if (entropy_quasi == true &&
      (group_entropy_com < 0 || group_entropy_rec < 0 ||
       group_entropy_lg < 0)) {
    cout << "Entropy-quasi option is invalid! ";
    cout << "Your entropy-quasi group(s) must not be negative!\n";
    exit(EXIT_FAILURE);
  }
  if (b_input < 0 || e_input < 0) {
    cout << "Begin/end input must be positive!\n";
    exit(EXIT_FAILURE);
  }
  if (b_input > e_input) {
    cout << "Begin input must be smaller than end!\n";
    exit(EXIT_FAILURE);
  }
  if (step_input <= 0) {
    cout << "Step input must be positive!\n";
    exit(EXIT_FAILURE);
  }
  if (nt == 0 || nt > std::thread::hardware_concurrency()) {
    nt = std::thread::hardware_concurrency();
    cout << "The number of threads is invalid. Will use " << nt
         << " threads!\n";
  }

  if (flag_gbsa) {
    if (entropy == false) {
      TitanUI::print_info(
          "Note: -gbsa option is on, will automatically activate -entropy "
          "option for Free Energy completeness!");
      entropy = true;
    }
    if (igb != 1 && igb != 2 && igb != 5) {
      cout << "GBSA Error: Invalid igb model! Only igb=1 (HCT), igb=2 (OBC1), "
              "and igb=5 (OBC2) are supported.\n";
      exit(EXIT_FAILURE);
    }
    if (ext_diel <= 0.0 || int_diel <= 0.0) {
      cout << "GBSA Error: Dielectric constants (ext_diel and int_diel) must "
              "be strictly positive!\n";
      exit(EXIT_FAILURE);
    }
    if (saltcon < 0.0) {
      cout << "GBSA Error: Salt concentration (saltcon) cannot be negative!\n";
      exit(EXIT_FAILURE);
    }
  }

  if (flag_listfile == true) {
    fstream list(list_filename, ios::in);
    if (!list.is_open()) {
      cout << "Cannot open " << list_filename << "!\n";
      flag_listfile = false;
      return 1;
    }
    list.close();
  }

  if (flag_energy_same_force == true) {
    emap_level = fmap_level;
    emap_type = fmap_type;
    print_energy_res = print_force_res;
    cout << "Note: Options for energy and force are same!\n";
  }

  if (lj.flag_lj) {
    if (lj.lj_cluster_cutoff <= 0.0) {
      TitanUI::print_error("LJ: -lj-cut must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (lj.lj_min_cluster_size <= 0) {
      TitanUI::print_error("LJ: -lj-min-cluster must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (lj.lj_rho_bin <= 0.0) {
      TitanUI::print_error("LJ: -lj-rho-bin must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (lj.lj_rho_max < 0.0) {
      TitanUI::print_error("LJ: -lj-rho-max must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (lj.lj_solid_q6 < 0.0 || lj.lj_solid_q6 > 1.0) {
      TitanUI::print_error("LJ: -lj-solid-q6 must be in [0, 1].");
      exit(EXIT_FAILURE);
    }
    if (lj.lj_solid_min_cn <= 0) {
      TitanUI::print_error("LJ: -lj-solid-min-cn must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (lj.lj_crystal_tol <= 0.0) {
      TitanUI::print_error("LJ: -lj-crystal-tol must be > 0.");
      exit(EXIT_FAILURE);
    }
  }

  if (flag_allostery) {
    if (allostery_tau <= 0) {
      TitanUI::print_error("Allostery: -allostery-tau must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (allostery_bins == 0) {
      TitanUI::print_error(
          "Allostery: -allostery-bins must be > 0 or -1 for auto.");
      exit(EXIT_FAILURE);
    }
    if (allostery_bins < -1) {
      TitanUI::print_error(
          "Allostery: -allostery-bins must be >= 1 or -1 for auto.");
      exit(EXIT_FAILURE);
    }
    if (allostery_cutoff < 0.0) {
      TitanUI::print_error("Allostery: -allostery-cut must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (allostery_knn < 0) {
      TitanUI::print_error("Allostery: -allostery-knn must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (allostery_knn_use_contact_prior && allostery_knn_contact_cutoff <= 0.0) {
      TitanUI::print_error(
          "Allostery: -allostery-knn-contact-cut must be > 0.");
      exit(EXIT_FAILURE);
    }
    std::string estimator = allostery_estimator;
    std::transform(estimator.begin(), estimator.end(), estimator.begin(),
                   [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (estimator == "histogram")
      estimator = "hist";
    if (estimator != "hist" && estimator != "knn" && estimator != "kde") {
      TitanUI::print_error(
          "Allostery: -allostery-estimator must be one of hist, knn, kde.");
      exit(EXIT_FAILURE);
    }
    allostery_estimator = estimator;
    if (allostery_kde_bandwidth < 0.0) {
      TitanUI::print_error("Allostery: -allostery-kde-bandwidth must be >= 0.");
      exit(EXIT_FAILURE);
    }
  }

  const bool lj_trace_free_kinetics =
      lj.flag_lj && (msm.flag_tica || msm.flag_msm || dham.flag_dham ||
                     tpt.flag_tpt ||
                     umap.flag_umap) &&
      !flag_network_explicit && !flag_pca && !flag_allostery &&
      !community.flag_community && !lspmd.flag_lspmd;
  if (lj_trace_free_kinetics && flag_network) {
    flag_network = false;
    TitanUI::print_info(
        "[LJ] Using permutation-invariant LJ phase-2 features for "
        "TICA/MSM/UMAP/DHAM; skipping auto-enabled DCCM network setup.");
  }

  const bool has_explicit_pbc_mode =
      pbc.solv_pbcn || pbc.solv_pbctopo || pbc.solv_pbc == 1 ||
      pbc.solv_pbc_int == 1 ||
      pbc.solv_pbc_clust || pbc.flag_pbc_ref || pbc.pbc_mem ||
      pbc.pbc_envi;
  std::transform(pbc.pbctopo_anchor_backend.begin(),
                 pbc.pbctopo_anchor_backend.end(),
                 pbc.pbctopo_anchor_backend.begin(), [](unsigned char c) {
                   return static_cast<char>(std::tolower(c));
                 });
  if (pbc.pbctopo_anchor_backend != "auto" &&
      pbc.pbctopo_anchor_backend != "prim" &&
      pbc.pbctopo_anchor_backend != "spatial-all") {
    TitanUI::print_error(
        "PBCTOPO: -pbctopo-anchor-backend must be one of auto, prim, spatial-all.");
    exit(EXIT_FAILURE);
  }
  if (rg.flag_rg2 && rg.hist_bins <= 0) {
    TitanUI::print_error("RG2: -rg2-hist-bins must be a positive integer.");
    exit(EXIT_FAILURE);
  }
  if ((flag_nocal || rms.flag_rms || rg.flag_rg2 || lj.flag_lj) &&
      !has_explicit_pbc_mode) {
    pbc.solv_pbc_int = 1;
    pbc.solv_pbc = 0;
    pbc.chain_ref = 0;
    pbc.auto_enabled_pbci = true;
    vector<string> modules;
    if (flag_nocal)
      modules.push_back("Nocal");
    if (rms.flag_rms)
      modules.push_back("RMSD/RMSF");
    if (rg.flag_rg2)
      modules.push_back("Rg2");
    if (lj.flag_lj)
      modules.push_back("LJ");
    std::string reason;
    for (size_t i = 0; i < modules.size(); ++i) {
      if (i > 0)
        reason += (i + 1 == modules.size()) ? " and " : ", ";
      reason += modules[i];
    }
    if (modules.size() == 1 && modules.front() == "LJ")
      reason += " analysis";
    reason += " requested without an explicit PBC mode. Auto-enabling -pbci 0.";
    TitanUI::print_warning(reason);
  }

  if (pbc.solv_pbc_int == 1 && pbc.solv_pbc_clust == true) {
    cout << "WARNING!! Both intensive pbc solver and clustered pbc solver are "
            "turned on, will turn off clustered pbc solver!\n";
    pbc.solv_pbc_clust = false;
  }
  /*if (flag_cal_force == true && print_energy_res == false)
  {
          cout << "WARNING! -cal-force option is on, will activate -energy-res
  option\n"; print_energy_res = true;
  }*/

  if (flag_maxforce == true && flag_cal_force == false) {
    cout << "Note: -max-force option is on, will activate -cal-force option!\n";
    flag_cal_force = true;
  }

  if (flag_topo_1 == false)
    mol_2_type = "protein";

  if (mol_2_type.find("protein") != string::npos &&
      mol_1_type.find("protein") != string::npos)
    emap_type = "residue-residue";

  if (group_1 == group_2)
    self_int = true;

  transform(eu.begin(), eu.end(), eu.begin(), ::toupper);
  if (eu != "KJ/MOL" && eu != "KCAL/MOL") {
    cout << "Energy unit is invalid, only kcal/mol and kJ/mol are supported. "
            "Will use the default unit (kJ/mol)\n";
    eu = "kJ/mol";
  } else if (eu == "KJ/MOL")
    eu = "kJ/mol";
  else {
    eu = "kcal/mol";
    eu_par = 4.184;
  }

  transform(fu.begin(), fu.end(), fu.begin(), ::toupper);
  if (fu != "PN" && fu != "KJ/MOL/NM") {
    cout << "Force unit is invalid, only pN and kJ/mol/nm are supported. Will "
            "use the default unit (pN)\n";
    fu = "pN";
  } else if (fu == "PN") {
    fu = "pN";
  } else {
    fu = "kJ/mol/nm";
    fu_par = 1.;
  }

  transform(tu.begin(), tu.end(), tu.begin(), ::toupper);
  if (tu != "PS" && tu != "NS") {
    cout << "Time unit is invalid, only ps and ns are supported. Will use the "
            "default unit (ps)\n";
  } else if (tu == "PS")
    tu = "ps";
  else
    tu = "ns";

  if (tu == "ns") {
    begin_time *= 1000.;
    end_time *= 1000.;
    dt *= 1000;
  }

  if (flag_index == false) {
    fstream index_file("index.ndx", ios::in);
    if (index_file.is_open())
      flag_index = true;
  }

  return 0;
}
