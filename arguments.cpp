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
  std::transform(eq_weight_source.begin(), eq_weight_source.end(),
                 eq_weight_source.begin(), [](unsigned char c) {
                   return static_cast<char>(std::tolower(c));
                 });
  if (eq_weight_source != "auto" && eq_weight_source != "uniform" &&
      eq_weight_source != "dham" && eq_weight_source != "mbar") {
    TitanUI::print_error(
        "-eq-weight-source must be auto, uniform, dham, or mbar.");
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
  if (graph.flag_graph || community.flag_community) {
    graph.flag_graph = true;
    if (graph.filter_topk < 0) {
      TitanUI::print_error("Graph: -graph-filter-topk must be >= 0.");
      exit(EXIT_FAILURE);
    }

    std::transform(graph.community_method.begin(), graph.community_method.end(),
                   graph.community_method.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (graph.community_method.empty())
      graph.community_method = "auto";
    if (graph.community_method != "auto" && graph.community_method != "off" &&
        graph.community_method != "louvain") {
      TitanUI::print_error(
          "Graph: -graph-community must be one of auto, off, louvain.");
      exit(EXIT_FAILURE);
    }
    if (community.flag_community && graph.community_method == "auto")
      graph.community_method = community.method;

    std::transform(graph.filter_component.begin(), graph.filter_component.end(),
                   graph.filter_component.begin(), [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (graph.filter_component.empty())
      graph.filter_component = "all";
    if (graph.filter_component != "all" &&
        graph.filter_component != "largest") {
      TitanUI::print_error(
          "Graph: -graph-filter-component must be one of all or largest.");
      exit(EXIT_FAILURE);
    }

    std::vector<std::string> normalized_sources;
    normalized_sources.reserve(graph.sources.size());
    bool use_auto_sources = graph.sources.empty();
    for (const std::string &raw_source : graph.sources) {
      std::string source = raw_source;
      std::transform(source.begin(), source.end(), source.begin(),
                     [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                     });
      if (source.empty())
        continue;
      if (source == "auto") {
        use_auto_sources = true;
        continue;
      }
      if (source == "network")
        source = "dccm";
      else if (source == "allostery_gc" || source == "allostery-gc-graph")
        source = "allostery-gc";
      else if (source == "allostery_te" || source == "allostery-te-graph")
        source = "allostery-te";
      if (source != "dccm" && source != "lspmd" && source != "allostery-gc" &&
          source != "allostery-te" && source != "tpt") {
        TitanUI::print_error(
            "Graph: -graph-source must be one of dccm, lspmd, allostery-gc, allostery-te, tpt.");
        exit(EXIT_FAILURE);
      }
      if (std::find(normalized_sources.begin(), normalized_sources.end(),
                    source) == normalized_sources.end()) {
        normalized_sources.push_back(source);
      }
    }

    if (!use_auto_sources && community.flag_community &&
        std::find(normalized_sources.begin(), normalized_sources.end(),
                  "dccm") == normalized_sources.end()) {
      normalized_sources.insert(normalized_sources.begin(), "dccm");
    }

    if (use_auto_sources || normalized_sources.empty()) {
      graph.sources = {"auto"};
    } else {
      graph.sources = std::move(normalized_sources);
    }

    const auto graph_has_source = [&](const char *name) {
      return std::find(graph.sources.begin(), graph.sources.end(), name) !=
             graph.sources.end();
    };
    const bool graph_auto_sources =
        graph.sources.size() == 1 && graph.sources.front() == "auto";

    if (!graph_auto_sources) {
      if (graph_has_source("dccm")) {
        flag_network = true;
        flag_network_explicit = true;
      }
      if (graph_has_source("allostery-gc") || graph_has_source("allostery-te")) {
        flag_allostery = true;
        flag_network = true;
      }
      if (graph_has_source("lspmd")) {
        lspmd.flag_lspmd = true;
      }
      if (graph_has_source("tpt")) {
        tpt.flag_tpt = true;
      }
    }

    const bool dccm_source_active =
        flag_network_explicit || flag_pca || flag_allostery ||
        corr_knn.flag_corr_knn || community.flag_community;
    const bool allostery_source_active = flag_allostery;
    const bool lspmd_source_active = lspmd.flag_lspmd;
    const bool tpt_source_active = tpt.flag_tpt;
    if (graph_auto_sources && !dccm_source_active && !allostery_source_active &&
        !lspmd_source_active && !tpt_source_active) {
      TitanUI::print_error(
          "Graph: -graph with auto sources requires at least one active source producer (dccm/network, allostery, lspmd, or tpt).");
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

  if (pbe.flag_pbe) {
    if (pbe.pbsa) {
      pbe.dgsolv = true;
      pbe.flag_pbe = true;
    }
    if (flag_nocal) {
      TitanUI::print_error(
          "PBE: -pbe requires charge/radius initialization; remove -nocal for "
          "the v1 frame-based solver.");
      exit(EXIT_FAILURE);
    }
    if (pbe.scale <= 0.0) {
      TitanUI::print_error("PBE: -pbe-scale must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.margin < 0.0) {
      TitanUI::print_error("PBE: -pbe-margin must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.int_diel <= 0.0 || pbe.ext_diel <= 0.0) {
      TitanUI::print_error(
          "PBE: -pbe-int-diel and -pbe-ext-diel must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.salt < 0.0) {
      TitanUI::print_error("PBE: -pbe-salt must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.temp <= 0.0) {
      TitanUI::print_error("PBE: -pbe-temp must be > 0 K.");
      exit(EXIT_FAILURE);
    }
    if (pbe.model != "linear" && pbe.model != "nonlinear") {
      TitanUI::print_error("PBE: -pbe-model must be linear or nonlinear.");
      exit(EXIT_FAILURE);
    }
    if (pbe.biomodel != "pbe" && pbe.biomodel != "rpbe") {
      TitanUI::print_error("PBE: -pbe-biomodel must be pbe or rpbe.");
      exit(EXIT_FAILURE);
    }
    if (pbe.dielectric_model != "two-dielectric" &&
        pbe.dielectric_model != "gaussian") {
      TitanUI::print_error(
          "PBE: -pbe-dielectric must be two-dielectric or gaussian.");
      exit(EXIT_FAILURE);
    }
    if (pbe.boundary != "coulombic" && pbe.boundary != "dipolar" &&
        pbe.boundary != "zero") {
      TitanUI::print_error(
          "PBE: -pbe-boundary must be coulombic, dipolar, or zero.");
      exit(EXIT_FAILURE);
    }
    if (pbe.charge_assignment != "trilinear" &&
        pbe.charge_assignment != "delphi" &&
        pbe.charge_assignment != "gaussian") {
      TitanUI::print_error(
          "PBE: -pbe-charge-assignment must be trilinear, delphi, or gaussian.");
      exit(EXIT_FAILURE);
    }
    if (pbe.surface != "vdw" && pbe.surface != "sas" &&
        pbe.surface != "ses" && pbe.surface != "delphi-ses" &&
        pbe.surface != "delphi-vdw" && pbe.surface != "gaussian" &&
        pbe.surface != "gaussian-cutoff") {
      TitanUI::print_error(
          "PBE: -pbe-surface must be vdw, sas, ses, delphi-ses, delphi-vdw, gaussian, or gaussian-cutoff.");
      exit(EXIT_FAILURE);
    }
    if (pbe.radius_convention != "auto" &&
        pbe.radius_convention != "native" &&
        pbe.radius_convention != "delphi-bondi") {
      TitanUI::print_error(
          "PBE: -pbe-radius-convention must be auto, native, or delphi-bondi.");
      exit(EXIT_FAILURE);
    }
    if (pbe.probe_radius < 0.0) {
      TitanUI::print_error("PBE: -pbe-probe must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.ion_access != "dielectric" && pbe.ion_access != "stern") {
      TitanUI::print_error(
          "PBE: -pbe-ion-access must be dielectric or stern.");
      exit(EXIT_FAILURE);
    }
    if (pbe.ion_radius < 0.0) {
      TitanUI::print_error("PBE: -pbe-ion-radius must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_max_iter <= 0) {
      TitanUI::print_error("PBE: -pbe-nl-max-iter must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_inner_iter <= 0) {
      TitanUI::print_error("PBE: -pbe-nl-inner-iter must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_stagnation <= 0) {
      TitanUI::print_error("PBE: -pbe-nl-stagnation must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_block <= 0) {
      TitanUI::print_error("PBE: -pbe-nl-block must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_dphi <= 0.0) {
      TitanUI::print_error("PBE: -pbe-nl-dphi must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_relax != "auto") {
      if (pbe.nl_relax_value < 0.05 || pbe.nl_relax_value > 1.0) {
        TitanUI::print_error(
            "PBE: -pbe-nl-relax must be auto or a value in [0.05, 1.0].");
        exit(EXIT_FAILURE);
      }
    }
    if (pbe.nl_clip <= 0.0) {
      TitanUI::print_error("PBE: -pbe-nl-clip must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_continuation != "off" && pbe.nl_continuation != "auto") {
      TitanUI::print_error(
          "PBE: -pbe-nl-continuation must be off or auto.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_line_search != "off" && pbe.nl_line_search != "auto") {
      TitanUI::print_error(
          "PBE: -pbe-nl-line-search must be off or auto.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_residual <= 0.0) {
      TitanUI::print_error("PBE: -pbe-nl-residual must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_solver != "picard" && pbe.nl_solver != "newton") {
      TitanUI::print_error("PBE: -pbe-nl-solver must be picard or newton.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_inner != "auto" && pbe.nl_inner != "sor" &&
        pbe.nl_inner != "vcycle") {
      TitanUI::print_error("PBE: -pbe-nl-inner must be auto, sor, or vcycle.");
      exit(EXIT_FAILURE);
    }
    if (pbe.nl_multigrid != "off" && pbe.nl_multigrid != "fas") {
      TitanUI::print_error("PBE: -pbe-nl-multigrid must be off or fas.");
      exit(EXIT_FAILURE);
    }
    if (pbe.gaussian_exp <= 0.0 || pbe.gaussian_sigma <= 0.0) {
      TitanUI::print_error(
          "PBE: -pbe-gaussian-exp and -pbe-gaussian-sigma must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.density_cutoff < 0.0 || pbe.density_cutoff > 1.0) {
      TitanUI::print_error("PBE: -pbe-density-cutoff must be in [0, 1].");
      exit(EXIT_FAILURE);
    }
    if (pbe.surface_cutoff != "auto" &&
        (pbe.surface_cutoff_value < 0.0 || pbe.surface_cutoff_value > 1.0)) {
      TitanUI::print_error(
          "PBE: -pbe-surface-cutoff must be auto or a value in [0, 1].");
      exit(EXIT_FAILURE);
    }
    if (pbe.gap_diel <= 0.0) {
      TitanUI::print_error("PBE: -pbe-gap-diel must be > 0.");
      exit(EXIT_FAILURE);
    }
    if ((pbe.test_grid_size > 0 || pbe.test_center_user) &&
        pbe.test_pqr.empty()) {
      TitanUI::print_error(
          "PBE: -pbe-test-grid-size and -pbe-test-center require "
          "-pbe-test-pqr.");
      exit(EXIT_FAILURE);
    }
    if (pbe.test_grid_size != 0 && pbe.test_grid_size < 3) {
      TitanUI::print_error("PBE: -pbe-test-grid-size must be >= 3.");
      exit(EXIT_FAILURE);
    }
    if (pbe.parity_dump != "off" && pbe.parity_dump != "summary" &&
        pbe.parity_dump != "points" && pbe.parity_dump != "all") {
      TitanUI::print_error(
          "PBE: -pbe-parity-dump must be off, summary, points, or all.");
      exit(EXIT_FAILURE);
    }
    if (pbe.parity_map_dump != "off" && pbe.parity_map_dump != "samples" &&
        pbe.parity_map_dump != "full") {
      TitanUI::print_error(
          "PBE: -pbe-parity-map-dump must be off, samples, or full.");
      exit(EXIT_FAILURE);
    }
    if (pbe.test_charge_grid != "native" &&
        pbe.test_charge_grid != "oracle") {
      TitanUI::print_error(
          "PBE: -pbe-test-charge-grid must be native or oracle.");
      exit(EXIT_FAILURE);
    }
    if (pbe.test_surface_map != "native" &&
        pbe.test_surface_map != "oracle") {
      TitanUI::print_error(
          "PBE: -pbe-test-surface-map must be native or oracle.");
      exit(EXIT_FAILURE);
    }
    if (pbe.test_midpoint_epsilon != "native" &&
        pbe.test_midpoint_epsilon != "oracle") {
      TitanUI::print_error(
          "PBE: -pbe-test-midpoint-epsilon must be native or oracle.");
      exit(EXIT_FAILURE);
    }
    if (pbe.test_boundary_map != "native" &&
        pbe.test_boundary_map != "oracle") {
      TitanUI::print_error(
          "PBE: -pbe-test-boundary-map must be native or oracle.");
      exit(EXIT_FAILURE);
    }
    if (pbe.test_nl_source != "titan" &&
        pbe.test_nl_source != "delphi-dimensionless") {
      TitanUI::print_error(
          "PBE: -pbe-test-nl-source must be titan or delphi-dimensionless.");
      exit(EXIT_FAILURE);
    }
    if (pbe.test_energy_convention != "auto" &&
        pbe.test_energy_convention != "atom" &&
        pbe.test_energy_convention != "grid" &&
        pbe.test_energy_convention != "delphi-corrected" &&
        pbe.test_energy_convention != "delphi-iscm") {
      TitanUI::print_error(
          "PBE: -pbe-test-energy-convention must be auto, atom, grid, delphi-corrected, or delphi-iscm.");
      exit(EXIT_FAILURE);
    }
    if ((pbe.parity_dump != "off" || pbe.parity_map_dump != "off" ||
         pbe.test_charge_grid != "native" ||
         pbe.test_surface_map != "native" ||
         pbe.test_midpoint_epsilon != "native" ||
         pbe.test_boundary_map != "native" ||
         pbe.test_nl_source != "titan" || pbe.test_newton_derivative ||
         pbe.test_energy_convention != "auto") &&
        pbe.test_pqr.empty()) {
      TitanUI::print_error(
          "PBE: hidden parity diagnostics require -pbe-test-pqr.");
      exit(EXIT_FAILURE);
    }
    if (pbe.dgsolv && !pbe.test_pqr.empty()) {
      TitanUI::print_error(
          "PBE: -pbe-dgsolv is not compatible with the hidden PQR fixture "
          "adapter; run it on a native two-molecule TITAN input.");
      exit(EXIT_FAILURE);
    }
    if (pbe.pbsa && !pbe.test_pqr.empty()) {
      TitanUI::print_error(
          "PBSA: -pbsa is not compatible with the hidden PQR fixture adapter; "
          "run it on a native two-molecule TITAN input.");
      exit(EXIT_FAILURE);
    }
    if (pbe.pbsa_components != "auto" &&
        pbe.pbsa_components != "complex,receptor,ligand") {
      TitanUI::print_error(
          "PBSA: -pbsa-components must be auto or complex,receptor,ligand.");
      exit(EXIT_FAILURE);
    }
    if (pbe.pbsa_decomp != "off" && pbe.pbsa_decomp != "molecule" &&
        pbe.pbsa_decomp != "residue") {
      TitanUI::print_error(
          "PBSA: -pbsa-decomp must be off, molecule, or residue.");
      exit(EXIT_FAILURE);
    }
    if (pbe.pbsa_residue_frames != "off" &&
        pbe.pbsa_residue_frames != "all") {
      TitanUI::print_error(
          "PBSA: -pbsa-residue-frames must be off or all.");
      exit(EXIT_FAILURE);
    }
    if (pbe.pbsa_ie_block != "auto") {
      try {
        const int block_size = std::stoi(pbe.pbsa_ie_block);
        if (block_size <= 0) {
          TitanUI::print_error(
              "PBSA: -pbsa-ie-block must be auto or a positive integer.");
          exit(EXIT_FAILURE);
        }
      } catch (...) {
        TitanUI::print_error(
            "PBSA: -pbsa-ie-block must be auto or a positive integer.");
        exit(EXIT_FAILURE);
      }
    }
    if (pbe.pbsa_ie_min_samples <= 0) {
      TitanUI::print_error("PBSA: -pbsa-ie-min-samples must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.pbsa_ie_filter_sigma < 0.0) {
      TitanUI::print_error("PBSA: -pbsa-ie-filter-sigma must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.pbsa_perf != "off" && pbe.pbsa_perf != "summary" &&
        pbe.pbsa_perf != "detail") {
      TitanUI::print_error("PBSA: -pbsa-perf must be off, summary, or detail.");
      exit(EXIT_FAILURE);
    }
    if (pbe.max_iter <= 0) {
      TitanUI::print_error("PBE: -pbe-max-iter must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.rms_tol <= 0.0 || pbe.dphi_tol <= 0.0) {
      TitanUI::print_error("PBE: convergence tolerances must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.max_grid_points < 27) {
      TitanUI::print_error("PBE: -pbe-max-grid must allow at least 27 points.");
      exit(EXIT_FAILURE);
    }
    if (pbe.omega != "auto" &&
        (pbe.omega_value <= 0.0 || pbe.omega_value >= 2.0)) {
      TitanUI::print_error("PBE: -pbe-omega must be auto or in (0, 2).");
      exit(EXIT_FAILURE);
    }
    if (pbe.focus != "off" && pbe.focus != "auto" &&
        pbe.focus != "delphi") {
      TitanUI::print_error("PBE: -pbe-focus must be off, auto, or delphi.");
      exit(EXIT_FAILURE);
    }
    if (pbe.multigrid != "off" && pbe.multigrid != "focus" &&
        pbe.multigrid != "auto" && pbe.multigrid != "vcycle") {
      TitanUI::print_error(
          "PBE: -pbe-multigrid must be off, focus, auto, or vcycle.");
      exit(EXIT_FAILURE);
    }
    if (pbe.focus_levels < 1 || pbe.focus_levels > 4) {
      TitanUI::print_error("PBE: -pbe-focus-levels must be in [1, 4].");
      exit(EXIT_FAILURE);
    }
    if (pbe.focus_margin < 0.0) {
      TitanUI::print_error("PBE: -pbe-focus-margin must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (pbe.mg_levels != "auto") {
      try {
        const int parsed_mg_levels = std::stoi(pbe.mg_levels);
        if (parsed_mg_levels < 1 || parsed_mg_levels > 4) {
          TitanUI::print_error("PBE: -pbe-mg-levels must be auto or in [1, 4].");
          exit(EXIT_FAILURE);
        }
      } catch (...) {
        TitanUI::print_error("PBE: -pbe-mg-levels must be auto or an integer.");
        exit(EXIT_FAILURE);
      }
    }
    if (pbe.mg_pre < 0 || pbe.mg_post < 0 || pbe.mg_coarse_iter <= 0) {
      TitanUI::print_error(
          "PBE: -pbe-mg-pre/post must be >= 0 and -pbe-mg-coarse-iter "
          "must be > 0.");
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
  if (hbk.flag_hbk) {
    if (hbk.hbk_maxlag <= 0) {
      TitanUI::print_error("HBond kinetics: -hbk-maxlag must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (!hbk.hbk_hbond_cutoff_explicit)
      hbk.hbk_hbond_cutoff = dis_h;
    if (!hbk.hbk_angle_explicit)
      hbk.hbk_angle_deg = angle_h;
    if (!hbk.hbk_water_cutoff_explicit)
      hbk.hbk_water_cutoff = d_crit_water;
    if (hbk.hbk_hbond_cutoff <= 0.0) {
      TitanUI::print_error("HBond kinetics: -hbk-hbond-cutoff must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (hbk.hbk_angle_deg <= 0.0 || hbk.hbk_angle_deg >= 180.0) {
      TitanUI::print_error("HBond kinetics: -hbk-angle must be in (0, 180).");
      exit(EXIT_FAILURE);
    }
    if (hbk.hbk_salt_cutoff <= 0.0) {
      TitanUI::print_error("HBond kinetics: -hbk-salt-cutoff must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (hbk.hbk_water_cutoff <= 0.0) {
      TitanUI::print_error("HBond kinetics: -hbk-water-cutoff must be > 0.");
      exit(EXIT_FAILURE);
    }

    std::string domain = hbk.hbk_domain;
    std::transform(domain.begin(), domain.end(), domain.begin(),
                   [](unsigned char c) {
                     return static_cast<char>(std::tolower(c));
                   });
    if (domain != "inter" && domain != "intra" && domain != "both") {
      TitanUI::print_error(
          "HBond kinetics: -hbk-domain must be one of inter, intra, both.");
      exit(EXIT_FAILURE);
    }
    hbk.hbk_domain = domain;

    std::vector<std::string> normalized_kinds;
    normalized_kinds.reserve(hbk.hbk_kinds.size());
    bool use_all = hbk.hbk_kinds.empty();
    for (const std::string &raw_kind : hbk.hbk_kinds) {
      std::string kind = raw_kind;
      std::transform(kind.begin(), kind.end(), kind.begin(),
                     [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                     });
      if (kind.empty())
        continue;
      if (kind == "all") {
        use_all = true;
        break;
      }
      if (kind != "hbond" && kind != "salt" && kind != "water") {
        TitanUI::print_error(
            "HBond kinetics: -hbk-kind must be one of all, hbond, salt, water.");
        exit(EXIT_FAILURE);
      }
      if (std::find(normalized_kinds.begin(), normalized_kinds.end(), kind) ==
          normalized_kinds.end()) {
        normalized_kinds.push_back(kind);
      }
    }
    if (use_all)
      normalized_kinds = {"hbond", "salt", "water"};
    hbk.hbk_kinds = std::move(normalized_kinds);
  }
  if (sq.flag_sq) {
    if (sq.sq_rdf_bin <= 0.0) {
      TitanUI::print_error("SQ: -sq-rdf-bin must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (!sq.sq_rdf_max_auto && sq.sq_rdf_max <= 0.0) {
      TitanUI::print_error("SQ: -sq-rdf-max must be > 0 or auto.");
      exit(EXIT_FAILURE);
    }
    if (!sq.sq_rdf_max_auto && sq.sq_rdf_max < sq.sq_rdf_bin) {
      TitanUI::print_error("SQ: -sq-rdf-max must be >= -sq-rdf-bin.");
      exit(EXIT_FAILURE);
    }
    if (sq.sq_q_bin <= 0.0) {
      TitanUI::print_error("SQ: -sq-q-bin must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (!sq.sq_q_max_auto && sq.sq_q_max <= 0.0) {
      TitanUI::print_error("SQ: -sq-q-max must be > 0 or auto.");
      exit(EXIT_FAILURE);
    }
    if (!sq.sq_q_max_auto && sq.sq_q_max < sq.sq_q_bin) {
      TitanUI::print_error("SQ: -sq-q-max must be >= -sq-q-bin.");
      exit(EXIT_FAILURE);
    }

    std::vector<std::string> normalized_channels;
    normalized_channels.reserve(sq.sq_channels.size());
    for (const std::string &raw_channel : sq.sq_channels) {
      std::string channel = raw_channel;
      std::transform(channel.begin(), channel.end(), channel.begin(),
                     [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                     });
      if (channel.empty())
        continue;
      if (channel != "all" && channel != "trace" && channel != "backbone" &&
          channel != "sticker") {
        TitanUI::print_error(
            "SQ: -sq-channel must be one of all, trace, backbone, sticker.");
        exit(EXIT_FAILURE);
      }
      if (std::find(normalized_channels.begin(), normalized_channels.end(),
                    channel) == normalized_channels.end()) {
        normalized_channels.push_back(channel);
      }
    }
    if (normalized_channels.empty())
      normalized_channels.push_back("all");
    sq.sq_channels = std::move(normalized_channels);
  }

  if (corr_knn.flag_corr_knn) {
    flag_network = true;
    if (corr_knn.estimator_k <= 0) {
      TitanUI::print_error("Corr-kNN: -corr-knn-k must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (corr_knn.tau_lag <= 0) {
      TitanUI::print_error("Corr-kNN: -corr-knn-tau must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (corr_knn.pair_cutoff_ang < 0.0) {
      TitanUI::print_error("Corr-kNN: -corr-knn-pair-cut must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (corr_knn.pair_topk < 0) {
      TitanUI::print_error("Corr-kNN: -corr-knn-pair-topk must be >= 0.");
      exit(EXIT_FAILURE);
    }
    if (corr_knn.pair_use_contact_prior &&
        corr_knn.pair_contact_cutoff_ang <= 0.0) {
      TitanUI::print_error(
          "Corr-kNN: -corr-knn-pair-contact-cut must be > 0.");
        exit(EXIT_FAILURE);
    }
  }

  if (validation.flag_validation) {
    flag_network = true;
    msm.flag_tica = true;
    msm.flag_msm = true;

    auto normalize_positive_int_list =
        [](std::vector<int> &values, const char *flag_name) {
          std::vector<int> normalized;
          normalized.reserve(values.size());
          for (int value : values) {
            if (value <= 0) {
              TitanUI::print_error(std::string("Validation: ") + flag_name +
                                   " entries must be > 0.");
              exit(EXIT_FAILURE);
            }
            if (std::find(normalized.begin(), normalized.end(), value) ==
                normalized.end()) {
              normalized.push_back(value);
            }
          }
          values = std::move(normalized);
        };

    normalize_positive_int_list(validation.ck_multiples,
                                "-validation-ck-multiples");
    if (!validation.vamp_components.empty()) {
      normalize_positive_int_list(validation.vamp_components,
                                  "-validation-vamp-components");
    }
    if (!validation.vamp_lags.empty()) {
      normalize_positive_int_list(validation.vamp_lags,
                                  "-validation-vamp-lags");
    }
    if (!validation.vamp_clusters.empty()) {
      normalize_positive_int_list(validation.vamp_clusters,
                                  "-validation-vamp-clusters");
    }
    if (validation.vamp_folds <= 0) {
      TitanUI::print_error("Validation: -validation-vamp-folds must be > 0.");
      exit(EXIT_FAILURE);
    }
    if (validation.bootstrap_samples < 0) {
      TitanUI::print_error(
          "Validation: -validation-bootstrap must be >= 0.");
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
    if (allostery_estimator == "knn" && allostery_estimator_k <= 0) {
      TitanUI::print_error("Allostery: -allostery-estimator-k must be > 0.");
      exit(EXIT_FAILURE);
    }
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
      !corr_knn.flag_corr_knn &&
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
  if ((flag_nocal || rms.flag_rms || rg.flag_rg2 || lj.flag_lj || sq.flag_sq) &&
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
    if (sq.flag_sq)
      modules.push_back("SQ");
    std::string reason;
    for (size_t i = 0; i < modules.size(); ++i) {
      if (i > 0)
        reason += (i + 1 == modules.size()) ? " and " : ", ";
      reason += modules[i];
    }
    if (modules.size() == 1 &&
        (modules.front() == "LJ" || modules.front() == "SQ"))
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
