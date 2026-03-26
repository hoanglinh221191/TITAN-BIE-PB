// engine_postproc.cpp — Post-loop analysis orchestration
// Extracted from impl_engine.cpp for faster incremental builds.
// Each function manages its own timing via std::chrono and logs via fstream.
#include "engine_postproc.h"
#include "titan_ui.h"
#include <Eigen/Dense>
#include <mutex>
#include <thread>

// Clang 21 COMDAT workaround: force emission of Matrix destructor
namespace {
void _eigen_force_dtor() {
  Eigen::MatrixXd m;
  (void)m;
  Eigen::VectorXd v;
  (void)v;
}
} // namespace

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

// ============================================================================
// 1. GBSA + ML + CAS + Nucleic Acid post-processing
// ============================================================================
void run_gbsa_postproc(const arguments &args, size_t tot_snap_allfiles,
                       MLFeatureExtractor &ml_extractor,
                       GBSA_Aggregator &global_gbsa_agg,
                       CASAnalyzer &cas_analyzer, NucleicProfiler &na_profiler,
                       bool global_has_nucleic_acid, info &info_output,
                       std::fstream &log) {
  if (args.flag_gbsa) {
    if (args.flag_ml) {
      if (tot_snap_allfiles > 0) {
        TitanUI::print_module_header("Machine Learning Feature Extraction");
        std::string ml_out_file = info_output.get_ml_features_filename();
        ml_extractor.export_features(ml_out_file);
        log << "Exported ML features to " << ml_out_file << "\n";
        TitanUI::print_info("Successfully exported 5 ML features to " +
                            ml_out_file);

        double predicted_pKD = ml_extractor.predict_affinity();
        double predicted_KD = std::pow(10.0, -predicted_pKD);

        std::ostringstream ss_out;
        ss_out << std::fixed << std::setprecision(4) << predicted_pKD
               << " (KD = " << std::scientific << std::setprecision(2)
               << predicted_KD << " M)";
        std::string pKD_str = ss_out.str();

        std::string pred_file = info_output.get_ml_prediction_filename();
        std::ofstream pred_out(pred_file);
        pred_out << "[Titan ML] Predicted pKD = " << pKD_str << "\n";
        pred_out.close();

        log << "[Titan ML] Predicted pKD = " << pKD_str << " (saved in "
            << pred_file << ")\n";
        TitanUI::print_info("[Titan ML] Predicted pKD = " + pKD_str);
      }
    }
    try {
      if (tot_snap_allfiles > 0) {
        log << "Printing GBSA statistics & decompositions...\n";
        TitanUI::print_module_header("GBSA Statistics & Decompositions");
        std::string out_dir = ".";
        global_gbsa_agg.print_all_stats(out_dir);
        TitanUI::print_info(
            "GBSA Decomposition files successfully written to " + out_dir);
      }
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] Exception in GBSA: " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] Unknown exception in GBSA!\n";
    }

    if (args.flag_cas) {
      try {
        if (tot_snap_allfiles > 0) {
          cas_analyzer.export_results(info_output.get_cas_receptor_filename(),
                                      info_output.get_cas_ligand_filename(),
                                      args.eu_par);
          cas_analyzer.export_hotspot_pdb(
              info_output.get_cas_hotspot_pdb_filename(), args.eu_par);
        }
      } catch (const std::exception &e) {
        std::cerr << "[ERROR] Exception in CAS: " << e.what() << "\n";
      } catch (...) {
        std::cerr << "[ERROR] Unknown exception in CAS!\n";
      }
    }
  }

  if (global_has_nucleic_acid) {
    try {
      if (tot_snap_allfiles > 0) {
        na_profiler.export_results("nucleic_acid_interactions.csv");
      }
    } catch (const std::exception &e) {
      std::cerr << "[ERROR] Exception in NA Profiler: " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ERROR] Unknown exception in NA Profiler!\n";
    }
  }
}

// ============================================================================
// 2. TICA + MSM post-processing
// ============================================================================
void run_tica_msm_postproc(const arguments &args,
                           std::unique_ptr<TitanKinetics> &tica_analyzer,
                           std::unique_ptr<TitanMSM> &msm_analyzer,
                           size_t num_threads, info &info_output,
                           std::fstream &log) {
  if (!(args.msm.flag_tica && tica_analyzer))
    return;

  try {
    auto start_tica = std::chrono::steady_clock::now();
    TitanUI::print_module_header("TICA Analysis");
    TitanUI::log_kv(log, "TICA Total Frames",
                    std::to_string(tica_analyzer->num_frames()));
    TitanUI::log_kv(log, "TICA Lagged Pairs",
                    std::to_string(tica_analyzer->num_pairs()));

    tica_analyzer->solve_tica();
    tica_analyzer->write_results(info_output.get_tica_eigenvalues_filename(),
                                 info_output.get_tica_eigenvectors_filename(),
                                 log);

    // Write TICA projections — proj_mat_ rows are ALREADY in TICA subspace.
    if (args.msm.flag_msm && msm_analyzer && msm_analyzer->n_frames() > 0) {
      tica_analyzer->write_projections(
          msm_analyzer->proj_mat_ref(), msm_analyzer->n_frames(),
          info_output.get_tica_projections_filename(), log);
    } else {
      TitanUI::log_kv(log, "TICA Projections",
                      "skipped (MSM not enabled or no frames)");
    }

    // --- MSM Phase 2 ---
    if (args.msm.flag_msm && msm_analyzer) {
      TitanUI::print_module_header("MSM Phase 2 — Collecting TICA Projections");
      if (msm_analyzer->n_frames() == 0) {
        TitanUI::print_info(
            "  [MSM] No projected frames collected during first pass "
            "(TICA was not pre-solved). "
            "MSM requires a pre-solved TICA model loaded from file. "
            "Skipping.");
      } else {
        TitanUI::log_kv(log, "MSM Projected Frames",
                        std::to_string(msm_analyzer->n_frames()));

        msm_analyzer->train_kmeans(args.msm.msm_kmeans_iter,
                                   static_cast<unsigned>(args.msm.msm_seed),
                                   static_cast<int>(num_threads));

        msm_analyzer->build_markov_model(args.msm.msm_lag);

        msm_analyzer->export_results(info_output, log, args.msm.msm_nts,
                                     args.msm.msm_run_its, args.msm.msm_run_pcca,
                                     args.msm.msm_macrostates, args.msm.msm_its_lags);
      }
    }

    auto time_tica = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_tica);
    TitanUI::log_time(log, "TICA + MSM Analysis", time_tica.count() / 1000.0);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] Exception in TICA/MSM: " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] Unknown exception in TICA/MSM!\n";
  }
}

// ============================================================================
// 3. LLPS post-processing
// ============================================================================
void run_llps_postproc(const arguments &args,
                       std::unique_ptr<TitanLLPS> &llps_analyzer,
                       size_t num_threads, info &info_output,
                       std::fstream &log) {
  if (!(args.llps.flag_llps && llps_analyzer &&
        llps_analyzer->n_frames_processed() > 0))
    return;

  try {
    auto start_llps = std::chrono::steady_clock::now();
    TitanUI::print_module_header("LLPS Analysis");
    llps_analyzer->finalize(static_cast<int>(num_threads), args.llps.llps_maxlag);
    llps_analyzer->export_results(info_output, log);
    auto time_llps = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_llps);
    TitanUI::log_time(log, "LLPS Analysis", time_llps.count() / 1000.0);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] Exception in LLPS: " << e.what() << "\n";
  }
}

// ============================================================================
// 4. DCCM + PCA post-processing
// ============================================================================
void run_dccm_pca_postproc(const arguments &args, size_t tot_snap_allfiles,
                           DCCMNetwork &dccm_analyzer, size_t num_threads,
                           std::chrono::milliseconds &time_network,
                           info &info_output, std::fstream &log) {
  if (!(args.flag_network || args.flag_pca))
    return;

  try {
    if (tot_snap_allfiles > 0) {
      auto start_dccm = std::chrono::steady_clock::now();
      TitanUI::print_module_header("DCCM Network Analysis");

      dccm_analyzer.calculate_and_export(num_threads);

      time_network += std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::steady_clock::now() - start_dccm);
      TitanUI::log_time(log, "Network Output", time_network.count() / 1000.0);

      // Allosteric Information Theory (GC + Transfer Entropy)
      if (args.flag_allostery) {
        auto start_allost = std::chrono::steady_clock::now();
        TitanUI::print_module_header("Allosteric Information Theory (GC + TE)");
        dccm_analyzer.calculate_allostery_information(
            num_threads, args.allostery_tau, args.allostery_bins);
        auto time_allost =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start_allost);
        TitanUI::log_time(log, "Allostery GC+TE",
                          time_allost.count() / 1000.0);
      }

      if (args.flag_pca) {
        auto ca_traj_aos = dccm_analyzer.get_ca_traj_aos();
        TitanPCA::execute_pca_fes_analysis(
            args, ca_traj_aos, dccm_analyzer.get_resnames(),
            dccm_analyzer.get_resnrs(), dccm_analyzer.get_chains(),
            dccm_analyzer.get_is_lig(), info_output, log);
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] Exception in DCCM/PCA: " << e.what() << "\n";
  }
}

// ============================================================================
// 5. LSP-MD post-processing (Kornev, JCTC 2025)
// ============================================================================
void run_lspmd_postproc(const arguments &args, DCCMNetwork &dccm_analyzer,
                        size_t num_threads, info &info_output,
                        std::fstream &log) {
  if (!args.lspmd.flag_lspmd)
    return;
  if (dccm_analyzer.get_total_frames() == 0) {
    TitanUI::print_warning("LSP-MD: No frames collected, skipping.");
    return;
  }

  try {
    auto start_lspmd = std::chrono::steady_clock::now();
    TitanUI::print_module_header(
        "LSP-MD: Local Spatial Pattern Network Analysis");
    TitanUI::log_kv(log, "LSP-MD Frames",
                    std::to_string(dccm_analyzer.get_total_frames()));
    TitanUI::log_kv(
        log, "LSP-MD Cutoff",
        std::to_string(args.lspmd.lspmd_cutoff) + " A");

    DCCMNetwork::LSPMDConfig cfg;
    cfg.cutoff_ang = args.lspmd.lspmd_cutoff;
    cfg.use_avg_cutoff = args.lspmd.lspmd_use_avg;

    DCCMNetwork::LSPMDResult result;
    dccm_analyzer.calculate_lspmd(num_threads, cfg, result);
    dccm_analyzer.export_lspmd(result);

    TitanUI::log_kv(log, "LSP-MD Edges",
                    std::to_string(result.edges.size()));
    auto time_lspmd = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_lspmd);
    TitanUI::log_time(log, "LSP-MD Analysis", time_lspmd.count() / 1000.0);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] Exception in LSP-MD: " << e.what() << "\n";
  }
}

// ============================================================================
// 6. DHAM post-processing
// ============================================================================
void run_dham_postproc(const arguments &args,
                       std::unique_ptr<TitanMSM> &msm_analyzer,
                       size_t num_threads, info &info_output,
                       std::fstream &log) {
  if (!args.dham.flag_dham || !msm_analyzer)
    return;

  auto start_dham = std::chrono::steady_clock::now();
  TitanUI::print_module_header("DHAM Post-Processing");

  // 1. Read temperatures
  if (args.dham.temps_file.empty()) {
    TitanUI::print_error("DHAM: -dham-temps not set. Skipping DHAM.");
    return;
  }
  auto temps = TitanDHAM::read_temps(args.dham.temps_file);
  if (temps.empty()) {
    TitanUI::print_error("DHAM: No temperatures read. Skipping.");
    return;
  }

  // 2. Read energy file list
  if (args.dham.ener_list_file.empty()) {
    TitanUI::print_error("DHAM: -dham-ener not set. Skipping DHAM.");
    return;
  }
  auto ener_paths = TitanDHAM::read_file_list(args.dham.ener_list_file);
  if (ener_paths.empty()) {
    TitanUI::print_error("DHAM: No energy files read. Skipping.");
    return;
  }

  const int n_replicas = static_cast<int>(temps.size());
  if (static_cast<int>(ener_paths.size()) != n_replicas) {
    TitanUI::print_error(
        "DHAM: #temps (" + std::to_string(n_replicas) +
        ") != #energy_files (" + std::to_string(ener_paths.size()) + ").");
    return;
  }

  // Validate n_replicas if explicitly set
  if (args.dham.n_replicas > 0 && args.dham.n_replicas != n_replicas) {
    TitanUI::print_warning(
        "DHAM: --dham-n-replicas=" + std::to_string(args.dham.n_replicas) +
        " != temps file count (" + std::to_string(n_replicas) +
        "). Using temps file count.");
  }

  TitanUI::print_info("DHAM: " + std::to_string(n_replicas) + " replicas, " +
                       "target T = " + std::to_string(args.dham.target_temp) + " K");
  // Fix C: Document the 1-file-per-replica assumption
  TitanUI::print_warning(
      "DHAM: Each input trajectory file is treated as one replica. "
      "Use 'gmx trjcat' to merge multi-part trajectories per replica "
      "before running DHAM.");

  // 3. Extract per-replica data from MSM
  if (!msm_analyzer->is_clustered()) {
    TitanUI::print_error("DHAM: MSM K-means has not run. "
        "Ensure -msm flag is set and enough frames are collected.");
    return;
  }

  auto per_replica_dtrajs     = msm_analyzer->get_per_replica_dtrajs();
  auto per_replica_timestamps = msm_analyzer->get_per_replica_timestamps();

  if (per_replica_dtrajs.empty()) {
    TitanUI::print_error("DHAM: No per-replica discrete trajectories. "
                          "Was replica_id passed to add_projected_frame?");
    return;
  }

  // 4. Parse + align energy .xvg files (Opt 3: threaded I/O + align)
  double effective_dt_energy = args.dham.dt_energy_ps;

  // Pre-detect dt_energy from first file (sequential, fast single read)
  if (effective_dt_energy <= 0.0) {
    auto probe = TitanDHAM::parse_xvg(ener_paths[0]);
    if (probe.size() >= 2) {
      effective_dt_energy = probe[1].first - probe[0].first;
      if (effective_dt_energy <= 0.0) effective_dt_energy = 2.0;
    } else {
      effective_dt_energy = 2.0;
    }
    TitanUI::print_info("DHAM: Auto-detected dt_energy = " +
        std::to_string(effective_dt_energy) + " ps from " + ener_paths[0]);
  }

  // Threaded parse + align: each thread handles a chunk of replicas
  std::unordered_map<int, std::vector<double>> per_replica_energies;
  std::mutex ener_mutex;
  int n_io_threads = std::min(static_cast<int>(num_threads), n_replicas);
  {
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(n_io_threads));
    int chunk = (n_replicas + n_io_threads - 1) / n_io_threads;
    for (int t = 0; t < n_io_threads; ++t) {
      int ks = t * chunk;
      int ke = std::min(ks + chunk, n_replicas);
      if (ks >= ke) continue;
      workers.emplace_back([&, ks, ke, effective_dt_energy]() {
        for (int k = ks; k < ke; ++k) {
          auto xvg_data = TitanDHAM::parse_xvg(
              ener_paths[static_cast<size_t>(k)]);
          if (xvg_data.empty()) continue;
          auto it_ts = per_replica_timestamps.find(k);
          if (it_ts == per_replica_timestamps.end()) continue;
          auto aligned = TitanDHAM::align_energy_to_dtraj(
              it_ts->second, xvg_data, effective_dt_energy);
          std::lock_guard<std::mutex> lock(ener_mutex);
          per_replica_energies[k] = std::move(aligned);
        }
      });
    }
    for (auto &w : workers) w.join();
  }

  // 5. Determine number of microstates from MSM
  const int n_states = msm_analyzer->n_clusters();

  // 6. Compute per-state mean energies (Opt 4: multi-threaded)
  auto state_energies = TitanDHAM::compute_state_energies(
      per_replica_dtrajs, per_replica_energies, n_states, n_replicas,
      static_cast<int>(num_threads));

  // 7. Create TitanDHAM and load data
  TitanDHAM dham(n_states, n_replicas, temps, args.dham.target_temp);

  // Load state energies (skip NaN = unvisited state/replica pairs)
  for (int k = 0; k < n_replicas; ++k) {
    for (int i = 0; i < n_states; ++i) {
      double e = state_energies(i, k);
      if (!std::isnan(e))
        dham.set_state_energy(k, i, e);
      // else: u_ik stays 0 — c_out will be 0 so no effect on solver
    }
  }

  // Build count matrices from per-replica dtrajs (use same lag as MSM)
  // Fix A: Skip transitions that cross trajectory breaks (time discontinuity)
  const int dham_lag = args.msm.msm_lag;
  TitanUI::print_info("DHAM: Using MSM lag = " + std::to_string(dham_lag) +
                       " frames for count matrix construction.");
  size_t total_counts = 0, skipped_breaks = 0;
  for (int k = 0; k < n_replicas; ++k) {
    auto it = per_replica_dtrajs.find(k);
    if (it == per_replica_dtrajs.end()) continue;
    const auto &dtraj = it->second;
    auto it_ts = per_replica_timestamps.find(k);
    const bool have_times = (it_ts != per_replica_timestamps.end() &&
                              it_ts->second.size() == dtraj.size());
    // Estimate expected gap = lag * dt_frame
    double dt_frame = 1.0, expected_gap = dham_lag;
    if (have_times && it_ts->second.size() > 1) {
      dt_frame = it_ts->second[1] - it_ts->second[0];
      expected_gap = dham_lag * dt_frame;
    }
    for (size_t t = 0; t + static_cast<size_t>(dham_lag) < dtraj.size(); ++t) {
      if (have_times) {
        double actual_gap = it_ts->second[t + static_cast<size_t>(dham_lag)]
                          - it_ts->second[t];
        if (std::abs(actual_gap - expected_gap) > 0.1 * std::abs(expected_gap)) {
          ++skipped_breaks;
          continue; // trajectory break — skip phantom transition
        }
      }
      dham.accumulate_counts(k, dtraj[t], dtraj[t + static_cast<size_t>(dham_lag)]);
      ++total_counts;
    }
  }
  if (skipped_breaks > 0) {
    TitanUI::print_warning("DHAM: Skipped " + std::to_string(skipped_breaks) +
        " transitions across trajectory breaks (kept " +
        std::to_string(total_counts) + ").");
  }

  // 8. Solve DHAM
  bool converged = dham.solve(args.dham.max_iter, args.dham.tolerance);
  TitanUI::log_kv(log, "DHAM Converged", converged ? "yes" : "no");

  // 9. Export results
  dham.write_results(info_output.get_dham_pi_filename(),
                     info_output.get_dham_fes_filename(),
                     info_output.get_dham_tpm_filename(), log,
                     args.eu_par);

  auto time_dham = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_dham);
  TitanUI::log_time(log, "DHAM Analysis", time_dham.count() / 1000.0);
}
