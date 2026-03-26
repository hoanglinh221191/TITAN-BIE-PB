// engine_postproc.cpp — Post-loop analysis orchestration
// Extracted from impl_engine.cpp for faster incremental builds.
// Each function manages its own timing via std::chrono and logs via fstream.
#include "engine_postproc.h"
#include "titan_alignment.h"
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

namespace {

bool is_trajectory_boundary(const std::vector<int> &replicas,
                            const std::vector<double> &times_ps, size_t prev,
                            size_t curr) {
  if (replicas.size() > curr && replicas[prev] != replicas[curr])
    return true;
  if (times_ps.size() > curr && times_ps[curr] <= times_ps[prev])
    return true;
  return false;
}

int select_effective_msm_clusters(const arguments &args,
                                  const TitanKinetics &tica_analyzer) {
  const int requested_k = std::max(2, args.msm.msm_clusters);
  if (args.msm.msm_clusters_user_set) {
    return requested_k;
  }

  // Auto mode is intentionally conservative for DHAM-facing workflows.
  // We only reduce the default cluster count when the TICA dataset is too
  // short to support a stable microstate model.
  const int64_t n_frames = tica_analyzer.num_frames();
  const int64_t n_pairs = tica_analyzer.num_pairs();
  const int64_t by_frames = std::max<int64_t>(2, n_frames / 5);
  const int64_t by_pairs =
      (n_pairs > 0) ? std::max<int64_t>(2, n_pairs / 3) : by_frames;
  const int64_t auto_cap =
      std::min<int64_t>(requested_k, std::min(by_frames, by_pairs));
  return std::clamp(static_cast<int>(auto_cap), 2, requested_k);
}

std::string describe_msm_cluster_selection(const arguments &args,
                                           int effective_k,
                                           const TitanKinetics &tica_analyzer) {
  const std::string mode =
      args.msm.msm_clusters_user_set ? "user" : "auto";
  return mode + " (requested=" + std::to_string(args.msm.msm_clusters) +
         ", effective=" + std::to_string(effective_k) +
         ", frames=" + std::to_string(tica_analyzer.num_frames()) +
         ", lagged_pairs=" + std::to_string(tica_analyzer.num_pairs()) + ")";
}

void rebuild_msm_projections_from_dccm(const arguments &args,
                                       const DCCMNetwork &dccm_analyzer,
                                       const TitanKinetics &tica_analyzer,
                                       TitanMSM &msm_analyzer,
                                       std::fstream &log) {
  const size_t n_frames = dccm_analyzer.get_total_frames();
  const int64_t expected_frames = tica_analyzer.num_frames();
  const int n_ca = static_cast<int>(dccm_analyzer.get_resnrs().size());
  const int feat_dim = tica_analyzer.feature_dim();

  if (n_frames == 0 || n_ca <= 0 || feat_dim <= 0) {
    TitanUI::print_warning(
        "MSM: Cannot rebuild TICA projections from cached DCCM trajectory.");
    return;
  }
  if (expected_frames > 0 && static_cast<int64_t>(n_frames) != expected_frames) {
    TitanUI::print_warning(
        "MSM: DCCM cache holds " + std::to_string(n_frames) +
        " frames, but TICA accumulated " + std::to_string(expected_frames) +
        ". Refusing partial projection rebuild.");
    return;
  }

  const auto &traj_x = dccm_analyzer.get_traj_x();
  const auto &traj_y = dccm_analyzer.get_traj_y();
  const auto &traj_z = dccm_analyzer.get_traj_z();
  const auto &traj_times = dccm_analyzer.get_traj_times();
  const auto &traj_replicas = dccm_analyzer.get_traj_replicas();
  const auto &global_indices = dccm_analyzer.get_global_frame_indices();

  if (traj_x.size() != n_frames || traj_y.size() != n_frames ||
      traj_z.size() != n_frames) {
    TitanUI::print_warning(
        "MSM: DCCM trajectory cache is inconsistent. Skipping projection rebuild.");
    return;
  }

  msm_analyzer.clear_projected();

  Eigen::VectorXd features(feat_dim);
  Eigen::MatrixXd mobile(n_ca, 3);
  Eigen::MatrixXd ref_centered(n_ca, 3);
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  bool ref_initialized = false;
  size_t n_breaks = 0;

  for (size_t f = 0; f < n_frames; ++f) {
    if (traj_x[f].size() < static_cast<size_t>(n_ca) ||
        traj_y[f].size() < static_cast<size_t>(n_ca) ||
        traj_z[f].size() < static_cast<size_t>(n_ca)) {
      TitanUI::print_warning(
          "MSM: Encountered truncated cached CA coordinates. Aborting projection rebuild.");
      msm_analyzer.clear_projected();
      return;
    }

    for (int a = 0; a < n_ca; ++a) {
      const size_t ai = static_cast<size_t>(a);
      mobile(a, 0) = traj_x[f][ai];
      mobile(a, 1) = traj_y[f][ai];
      mobile(a, 2) = traj_z[f][ai];
    }

    if (!ref_initialized) {
      ref_centered = mobile;
      ref_centered.rowwise() -= ref_centered.colwise().mean().eval();
      mobile = ref_centered;
      ref_initialized = true;
    } else {
      TitanAlignment::align_frame_inplace(ref_centered, mobile, R);
    }

    int feat_idx = 0;
    for (int a = 0; a < n_ca; ++a) {
      for (int b = a + 1; b < n_ca; ++b) {
        const double dx = mobile(a, 0) - mobile(b, 0);
        const double dy = mobile(a, 1) - mobile(b, 1);
        const double dz = mobile(a, 2) - mobile(b, 2);
        features[feat_idx++] = std::sqrt(dx * dx + dy * dy + dz * dz);
      }
    }

    const int replica_id =
        (traj_replicas.size() > f) ? traj_replicas[f] : -1;
    const int64_t global_idx =
        (global_indices.size() > f) ? global_indices[f]
                                    : static_cast<int64_t>(f);

    msm_analyzer.add_projected_frame(
        tica_analyzer.project_frame(features), replica_id, global_idx);
    if (traj_times.size() > f)
      msm_analyzer.record_frame_time(traj_times[f], replica_id);

    if (args.msm.msm_auto_boundary && f + 1 < n_frames &&
        is_trajectory_boundary(traj_replicas, traj_times, f, f + 1)) {
      msm_analyzer.add_trajectory_break();
      ++n_breaks;
    }
  }

  TitanUI::log_kv(log, "MSM Projection Rebuild", "from cached DCCM trajectory");
  TitanUI::log_kv(log, "MSM Projected Frames",
                  std::to_string(msm_analyzer.n_frames()));
  if (args.msm.msm_auto_boundary) {
    TitanUI::log_kv(log, "MSM Rebuilt Trajectory Breaks",
                    std::to_string(n_breaks));
  }
}

} // namespace

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
                           DCCMNetwork &dccm_analyzer,
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
    TitanUI::log_kv(log, "TICA Backend",
                    tica_backend_to_string(tica_analyzer->backend()));
    TitanUI::log_kv(log, "TICA Rank",
                    std::to_string(tica_analyzer->randomized_rank()));

    tica_analyzer->solve_tica();
    TitanUI::log_kv(log, "TICA Retained Components",
                    std::to_string(tica_analyzer->retained_components()));
    tica_analyzer->write_results(info_output.get_tica_eigenvalues_filename(),
                                 info_output.get_tica_eigenvectors_filename(),
                                 log);

    if (msm_analyzer) {
      TitanUI::print_module_header("TICA Projection Cache");
      if (msm_analyzer->n_frames() == 0 &&
          msm_analyzer->pending_feature_frames() > 0) {
        const int64_t pending_frames = msm_analyzer->pending_feature_frames();
        msm_analyzer->project_pending_features(*tica_analyzer);
        TitanUI::log_kv(log, "TICA Projection Rebuild",
                        "from TitanKinetics feature cache");
        TitanUI::log_kv(log, "TICA Pending Feature Frames",
                        std::to_string(pending_frames));
        TitanUI::log_kv(log, "TICA Projected Frames",
                        std::to_string(msm_analyzer->n_frames()));
      }

      if (msm_analyzer->n_frames() == 0) {
        rebuild_msm_projections_from_dccm(args, dccm_analyzer, *tica_analyzer,
                                          *msm_analyzer, log);
      }

      if (msm_analyzer->n_frames() > 0) {
        tica_analyzer->write_projections(
            msm_analyzer->proj_mat_ref(), msm_analyzer->n_frames(),
            info_output.get_tica_projections_filename(), log);
      } else {
        TitanUI::log_kv(log, "TICA Projections",
                        "skipped (no projected frames available)");
      }
    } else {
      TitanUI::log_kv(log, "TICA Projections",
                      "skipped (projection cache not enabled)");
    }

    // --- MSM Phase 2 ---
    if (args.msm.flag_msm && msm_analyzer) {
      const int effective_k =
          select_effective_msm_clusters(args, *tica_analyzer);
      msm_analyzer->reconfigure_clusters(effective_k);
      TitanUI::log_kv(log, "MSM Cluster Selection",
                      describe_msm_cluster_selection(args, effective_k,
                                                     *tica_analyzer));
      TitanUI::print_module_header("MSM Phase 2 — Collecting TICA Projections");
      if (msm_analyzer->n_frames() == 0 &&
          msm_analyzer->pending_feature_frames() > 0) {
        const int64_t pending_frames = msm_analyzer->pending_feature_frames();
        msm_analyzer->project_pending_features(*tica_analyzer);
        TitanUI::log_kv(log, "MSM Projection Rebuild",
                        "from TitanKinetics feature cache");
        TitanUI::log_kv(log, "MSM Pending Feature Frames",
                        std::to_string(pending_frames));
        TitanUI::log_kv(log, "MSM Projected Frames",
                        std::to_string(msm_analyzer->n_frames()));
      }

      if (msm_analyzer->n_frames() == 0) {
        rebuild_msm_projections_from_dccm(args, dccm_analyzer, *tica_analyzer,
                                          *msm_analyzer, log);
      }

      if (msm_analyzer->n_frames() > 0) {
        msm_analyzer->train_kmeans(args.msm.msm_kmeans_iter,
                                   static_cast<unsigned>(args.msm.msm_seed),
                                   static_cast<int>(num_threads));

        if (msm_analyzer->is_clustered())
          msm_analyzer->build_markov_model(args.msm.msm_lag);

        if (msm_analyzer->is_built()) {
          msm_analyzer->export_results(
              info_output, log, args.msm.msm_nts, args.msm.msm_run_its,
              args.msm.msm_run_pcca, args.msm.msm_macrostates,
              args.msm.msm_its_lags);
        }
      } else {
        TitanUI::log_kv(log, "TICA Projections",
                        "skipped (no projected frames available)");
      }
    } else {
      TitanUI::log_kv(log, "MSM Analysis", "skipped (module not enabled)");
    }

    auto time_tica = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_tica);
    TitanUI::log_time(log, args.msm.flag_msm ? "TICA + MSM Analysis"
                                             : "TICA Analysis",
                      time_tica.count() / 1000.0);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] Exception in TICA/MSM: " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] Unknown exception in TICA/MSM!\n";
  }
}

// ============================================================================
// 3. TICA 2D FES post-processing
// ============================================================================
void run_tica_fes_postproc(const arguments &args,
                           std::unique_ptr<TitanMSM> &msm_analyzer,
                           info &info_output, std::fstream &log,
                           const std::vector<double> &frame_weights) {
  if (!(args.msm.flag_tica && args.msm.flag_msm && msm_analyzer))
    return;
  if (msm_analyzer->n_frames() <= 0)
    return;

  try {
    auto start_tica_fes = std::chrono::steady_clock::now();
    TitanUI::print_module_header("TICA 2D Free Energy Surface");
    TitanPCA::execute_tica_fes_analysis(
        args, msm_analyzer->proj_mat_ref(), msm_analyzer->n_frames(),
        msm_analyzer->get_frame_global_indices(), info_output, log,
        frame_weights);
    auto time_tica_fes = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_tica_fes);
    TitanUI::log_time(log, "TICA 2D FES", time_tica_fes.count() / 1000.0);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] Exception in TICA 2D FES: " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] Unknown exception in TICA 2D FES!\n";
  }
}

void run_umap_postproc(const arguments &args,
                       std::unique_ptr<TitanMSM> &msm_analyzer,
                       std::unique_ptr<TitanUMAP> &umap_analyzer,
                       size_t num_threads,
                       const info &info_output, std::fstream &log) {
  if (!(args.umap.flag_umap && msm_analyzer && msm_analyzer->n_frames() > 0))
    return;

  try {
    auto start_umap = std::chrono::steady_clock::now();
    TitanUI::print_module_header("UMAP Analysis");

    TitanUMAP::Config cfg;
    cfg.n_neighbors = args.umap.umap_neighbors;
    cfg.n_components = args.umap.umap_components;
    cfg.n_epochs = args.umap.umap_epochs;
    cfg.min_dist = args.umap.umap_min_dist;
    cfg.spread = args.umap.umap_spread;
    cfg.seed = args.umap.umap_seed;
    cfg.n_threads = static_cast<int>(std::max<size_t>(1, num_threads));
    if (args.umap.umap_knn_mode == "approx") {
      cfg.knn_mode = TitanUMAP::KnnMode::Approx;
    } else if (args.umap.umap_knn_mode == "auto") {
      cfg.knn_mode = TitanUMAP::KnnMode::Auto;
    } else {
      cfg.knn_mode = TitanUMAP::KnnMode::Exact;
    }
    cfg.benchmark_knn = args.umap.umap_knn_benchmark;
    cfg.verbose = args.umap.umap_verbose;

    TitanUI::log_kv(log, "UMAP Input Frames",
                    std::to_string(msm_analyzer->n_frames()));
    TitanUI::log_kv(log, "UMAP Input Dimensions",
                    std::to_string(msm_analyzer->proj_mat_ref().cols()));
    TitanUI::log_kv(log, "UMAP Components",
                    std::to_string(cfg.n_components));
    TitanUI::log_kv(log, "UMAP Neighbors",
                    std::to_string(cfg.n_neighbors));
    TitanUI::log_kv(log, "UMAP Epochs", std::to_string(cfg.n_epochs));
    TitanUI::log_kv(log, "UMAP Min Dist", std::to_string(cfg.min_dist));
    TitanUI::log_kv(log, "UMAP Spread", std::to_string(cfg.spread));
    TitanUI::log_kv(log, "UMAP Seed", std::to_string(cfg.seed));
    TitanUI::log_kv(log, "UMAP Threads", std::to_string(cfg.n_threads));
    TitanUI::log_kv(log, "UMAP kNN Mode", args.umap.umap_knn_mode);
    TitanUI::log_kv(log, "UMAP kNN Benchmark",
                    cfg.benchmark_knn ? "on" : "off");

    umap_analyzer = std::make_unique<TitanUMAP>(cfg);
    umap_analyzer->fit_transform(msm_analyzer->proj_mat_ref(),
                                 msm_analyzer->n_frames());
    if (umap_analyzer->is_fitted()) {
      umap_analyzer->write_embedding(
          info_output.get_umap_embedding_filename(), log);
    }

    auto time_umap = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_umap);
    TitanUI::log_time(log, "UMAP Analysis", time_umap.count() / 1000.0);
  } catch (const std::exception &e) {
    std::cerr << "[ERROR] Exception in UMAP: " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ERROR] Unknown exception in UMAP!\n";
  }
}

// ============================================================================
// 4. LLPS post-processing
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
// 5. DCCM + PCA post-processing
// ============================================================================
void run_dccm_pca_postproc(const arguments &args, size_t tot_snap_allfiles,
                           DCCMNetwork &dccm_analyzer, size_t num_threads,
                           std::chrono::milliseconds &time_network,
                           info &info_output, std::fstream &log,
                           const std::vector<double> &frame_weights,
                           const std::vector<std::string> &replica_traj_files) {
  if (!(args.flag_network_explicit || args.flag_pca || args.flag_allostery))
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
        dccm_analyzer.calculate_allostery_information_optimized(
            num_threads, args.allostery_tau, args.allostery_bins);
        auto time_allost =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start_allost);
        TitanUI::log_time(log, "Allostery GC+TE",
                          time_allost.count() / 1000.0);
      }

      if (args.flag_pca) {
        auto ca_traj_aos = dccm_analyzer.get_ca_traj_aos();
        std::vector<double> pca_frame_weights;
        if (!frame_weights.empty()) {
          const auto &dccm_global_indices = dccm_analyzer.get_global_frame_indices();
          if (dccm_global_indices.size() != ca_traj_aos.size()) {
            TitanUI::print_warning(
                "PCA: DCCM global frame index count (" +
                std::to_string(dccm_global_indices.size()) +
                ") != PCA frame count (" + std::to_string(ca_traj_aos.size()) +
                "). Falling back to uniform weights.");
          } else {
            pca_frame_weights.assign(ca_traj_aos.size(), 0.0);
            size_t matched_frames = 0;
            size_t weighted_frames = 0;
            double mapped_weight_sum = 0.0;
            for (size_t f = 0; f < dccm_global_indices.size(); ++f) {
              int64_t global_idx = dccm_global_indices[f];
              if (global_idx < 0 ||
                  static_cast<size_t>(global_idx) >= frame_weights.size()) {
                continue;
              }
              const double w = frame_weights[static_cast<size_t>(global_idx)];
              pca_frame_weights[f] = w;
              ++matched_frames;
              if (w > 0.0)
                ++weighted_frames;
              mapped_weight_sum += w;
            }
            TitanUI::print_info(
                "PCA: Mapped DHAM weights onto " +
                std::to_string(weighted_frames) + "/" +
                std::to_string(ca_traj_aos.size()) +
                " frames (matched indices: " + std::to_string(matched_frames) +
                ", weight sum = " + std::to_string(mapped_weight_sum) + ").");
          }
        }
        TitanPCA::execute_pca_fes_analysis(
            args, ca_traj_aos, dccm_analyzer.get_resnames(),
            dccm_analyzer.get_resnrs(), dccm_analyzer.get_chains(),
            dccm_analyzer.get_is_lig(), info_output, log, pca_frame_weights,
            dccm_analyzer.get_traj_times(), dccm_analyzer.get_traj_replicas(),
            replica_traj_files);
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
                       std::fstream &log,
                       std::vector<double> &out_frame_weights) {
  out_frame_weights.clear();
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
  if (!msm_analyzer->is_built()) {
    TitanUI::print_error(
        "DHAM: MSM transition model has not been built successfully.");
    return;
  }

  const int dham_lag = msm_analyzer->msm_lag_time();
  if (dham_lag <= 0) {
    TitanUI::print_error("DHAM: invalid MSM lag time. Aborting DHAM.");
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
  // Let each replica auto-detect its own cadence unless the user pinned one.
  // Threaded parse + align: each thread handles a chunk of replicas.
  std::unordered_map<int, std::vector<double>> per_replica_energies;
  std::mutex ener_mutex;
  int n_io_threads = std::min(static_cast<int>(num_threads), n_replicas);
  {
    std::vector<std::jthread> workers;
    workers.reserve(static_cast<size_t>(n_io_threads));
    int chunk = (n_replicas + n_io_threads - 1) / n_io_threads;
    for (int t = 0; t < n_io_threads; ++t) {
      int ks = t * chunk;
      int ke = std::min(ks + chunk, n_replicas);
      if (ks >= ke) continue;
      workers.emplace_back([&, ks, ke]() {
        for (int k = ks; k < ke; ++k) {
          auto xvg_data = TitanDHAM::parse_xvg(
              ener_paths[static_cast<size_t>(k)]);
          if (xvg_data.empty()) continue;
          auto it_ts = per_replica_timestamps.find(k);
          if (it_ts == per_replica_timestamps.end()) continue;
          const double dt_energy_ps =
              (args.dham.dt_energy_ps > 0.0) ? args.dham.dt_energy_ps : 0.0;
          auto aligned = TitanDHAM::align_energy_to_dtraj(
              it_ts->second, xvg_data, dt_energy_ps);
          std::lock_guard<std::mutex> lock(ener_mutex);
          per_replica_energies[k] = std::move(aligned);
        }
      });
    }
  }
  if (args.dham.dt_energy_ps > 0.0) {
    TitanUI::print_info("DHAM: Using user-specified dt_energy = " +
                        std::to_string(args.dham.dt_energy_ps) + " ps.");
  } else {
    TitanUI::print_info(
        "DHAM: Auto-detecting dt_energy independently for each energy file.");
  }

  // 5. Determine number of microstates from MSM
  const int n_states = msm_analyzer->n_clusters();

  // 6. Compute per-state mean energies (Opt 4: multi-threaded)
  auto state_energies = TitanDHAM::compute_state_energies(
      per_replica_dtrajs, per_replica_energies, n_states, n_replicas,
      static_cast<int>(num_threads));

  // 7. Create TitanDHAM and load data
  TitanDHAM dham(n_states, n_replicas, temps, args.dham.target_temp);
  std::vector<std::vector<unsigned char>> valid_energy(
      static_cast<size_t>(n_replicas),
      std::vector<unsigned char>(static_cast<size_t>(n_states), 0));
  size_t valid_state_energy_pairs = 0;

  // Load state energies (skip NaN = unvisited state/replica pairs)
  for (int k = 0; k < n_replicas; ++k) {
    for (int i = 0; i < n_states; ++i) {
      double e = state_energies(i, k);
      if (std::isfinite(e)) {
        valid_energy[static_cast<size_t>(k)][static_cast<size_t>(i)] = 1;
        ++valid_state_energy_pairs;
        dham.set_state_energy(k, i, e);
      }
    }
  }
  const size_t total_state_energy_pairs =
      static_cast<size_t>(n_states) * static_cast<size_t>(n_replicas);
  const size_t missing_state_energy_pairs =
      total_state_energy_pairs - valid_state_energy_pairs;
  if (missing_state_energy_pairs > 0) {
    TitanUI::print_warning(
        "DHAM: Missing aligned energies for " +
        std::to_string(missing_state_energy_pairs) + "/" +
        std::to_string(total_state_energy_pairs) +
        " state/replica pairs. Transitions touching those states will be skipped.");
  }

  // Build count matrices from per-replica dtrajs (use same lag as MSM)
  // Skip transitions that cross trajectory breaks or states with missing energy.
  TitanUI::print_info("DHAM: Using MSM lag = " + std::to_string(dham_lag) +
                        " frames for count matrix construction.");
  size_t total_counts = 0;
  size_t skipped_breaks = 0;
  size_t skipped_missing_energy = 0;
  size_t skipped_invalid_states = 0;
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
      const int state_i = dtraj[t];
      const int state_j = dtraj[t + static_cast<size_t>(dham_lag)];
      if (state_i < 0 || state_i >= n_states || state_j < 0 || state_j >= n_states) {
        ++skipped_invalid_states;
        continue;
      }
      if (!valid_energy[static_cast<size_t>(k)][static_cast<size_t>(state_i)] ||
          !valid_energy[static_cast<size_t>(k)][static_cast<size_t>(state_j)]) {
        ++skipped_missing_energy;
        continue;
      }
      dham.accumulate_counts(k, state_i, state_j);
      ++total_counts;
    }
  }
  if (skipped_breaks > 0) {
    TitanUI::print_warning("DHAM: Skipped " + std::to_string(skipped_breaks) +
        " transitions across trajectory breaks (kept " +
        std::to_string(total_counts) + ").");
  }
  if (skipped_missing_energy > 0) {
    TitanUI::print_warning("DHAM: Skipped " +
        std::to_string(skipped_missing_energy) +
        " transitions involving states without aligned energies.");
  }
  if (skipped_invalid_states > 0) {
    TitanUI::print_warning("DHAM: Skipped " +
        std::to_string(skipped_invalid_states) +
        " transitions with invalid microstate indices.");
  }

  // 8. Solve DHAM
  bool converged = dham.solve(args.dham.max_iter, args.dham.tolerance);
  TitanUI::log_kv(log, "DHAM Converged", converged ? "yes" : "no");
  if (converged) {
    auto per_replica_weights = dham.get_frame_weights(per_replica_dtrajs);
    const auto &frame_replica_ids = msm_analyzer->get_frame_replica_ids();
    const auto &frame_global_indices = msm_analyzer->get_frame_global_indices();

    if (frame_replica_ids.size() != frame_global_indices.size()) {
      TitanUI::print_warning(
          "DHAM: projected frame metadata mismatch: replica_ids=" +
          std::to_string(frame_replica_ids.size()) + ", global_indices=" +
          std::to_string(frame_global_indices.size()) +
          ". Some frame weights may be dropped.");
    }

    int64_t max_global_index = -1;
    for (int64_t global_idx : frame_global_indices) {
      if (global_idx > max_global_index)
        max_global_index = global_idx;
    }
    if (max_global_index >= 0) {
      out_frame_weights.assign(static_cast<size_t>(max_global_index) + 1, 0.0);
    }

    std::unordered_map<int, size_t> per_replica_offsets;
    per_replica_offsets.reserve(per_replica_weights.size());
    std::vector<unsigned char> assigned(out_frame_weights.size(), 0);
    const size_t tracked_projected_frames =
        std::min(frame_replica_ids.size(), frame_global_indices.size());
    size_t mapped_frames = 0;
    size_t dropped_missing_replica_weights = 0;
    size_t dropped_invalid_global_index = 0;
    size_t duplicate_global_indices = 0;

    for (size_t frame_idx = 0; frame_idx < tracked_projected_frames; ++frame_idx) {
      const int replica_id = frame_replica_ids[frame_idx];
      auto it_weights = per_replica_weights.find(replica_id);
      size_t &replica_offset = per_replica_offsets[replica_id];
      if (it_weights == per_replica_weights.end() ||
          replica_offset >= it_weights->second.size()) {
        ++dropped_missing_replica_weights;
        continue;
      }

      const int64_t global_idx = frame_global_indices[frame_idx];
      if (global_idx < 0 ||
          static_cast<size_t>(global_idx) >= out_frame_weights.size()) {
        ++dropped_invalid_global_index;
        ++replica_offset;
        continue;
      }

      const double w = it_weights->second[replica_offset++];
      const size_t global_slot = static_cast<size_t>(global_idx);
      if (assigned[global_slot]) {
        ++duplicate_global_indices;
      }
      out_frame_weights[global_slot] = w;
      assigned[global_slot] = 1;
      ++mapped_frames;
    }

    double weight_sum = 0.0;
    for (double w : out_frame_weights)
      weight_sum += w;
    TitanUI::print_info(
        "DHAM: Computed " + std::to_string(mapped_frames) +
        " global frame weights for reweighted FES.");
    TitanUI::print_info("DHAM: Weight sum = " + std::to_string(weight_sum) +
                        " (expected ~1.0)");
    if (dropped_missing_replica_weights > 0) {
      TitanUI::print_warning(
          "DHAM: Dropped " +
          std::to_string(dropped_missing_replica_weights) +
          " projected frames while rebuilding global weights.");
    }
    if (dropped_invalid_global_index > 0) {
      TitanUI::print_warning(
          "DHAM: Dropped " + std::to_string(dropped_invalid_global_index) +
          " projected frames with invalid global indices.");
    }
    if (duplicate_global_indices > 0) {
      TitanUI::print_warning(
          "DHAM: Encountered " + std::to_string(duplicate_global_indices) +
          " duplicate global frame indices while assigning weights.");
    }

    // 9. Export results
    dham.write_results(info_output.get_dham_pi_filename(),
                       info_output.get_dham_fes_filename(),
                       info_output.get_dham_tpm_filename(), log,
                       args.eu_par);
  } else {
    TitanUI::print_warning(
        "DHAM: solve() did not converge - results were not written and "
        "PCA/TICA FES will use uniform weights. Increase --dham-max-iter "
        "or check the aligned energy data.");
  }

  auto time_dham = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_dham);
  TitanUI::log_time(log, "DHAM Analysis", time_dham.count() / 1000.0);
}
