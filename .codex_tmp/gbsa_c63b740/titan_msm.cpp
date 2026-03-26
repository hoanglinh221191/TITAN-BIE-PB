// ============================================================================
// TitanMSM — Markov State Model (Phase 2) Implementation
// ============================================================================
#include "titan_msm.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// ============================================================================
// Constructor
// ============================================================================
TitanMSM::TitanMSM(int n_clusters, int n_tica_components, double dt_ns)
    : n_k_(n_clusters), n_comp_(n_tica_components), dt_ns_(dt_ns),
      msm_lag_time_(0),
      centroids_(Eigen::MatrixXd::Zero(n_clusters, n_tica_components)),
      count_matrix_(Eigen::MatrixXd::Zero(n_clusters, n_clusters)),
      tpm_(Eigen::MatrixXd::Zero(n_clusters, n_clusters)),
      stationary_dist_(
          Eigen::VectorXd::Constant(n_clusters, 1.0 / n_clusters)) {}

// ============================================================================
// Trajectory accumulation
// ============================================================================
void TitanMSM::add_projected_frame(const Eigen::VectorXd &proj_x,
                                   int replica_id) {
  // Fix 3: grow proj_mat_ with doubling strategy (amortized O(1), contiguous)
  if (n_proj_frames_ >= proj_mat_.rows()) {
    int64_t new_rows = std::max(int64_t(1024), proj_mat_.rows() * 2);
    proj_mat_.conservativeResize(new_rows, n_comp_);
    // Reserve replica tracking vectors to match
    frame_replica_ids_.reserve(static_cast<size_t>(new_rows));
    frame_timestamps_.reserve(static_cast<size_t>(new_rows));
  }
  proj_mat_.row(n_proj_frames_++) = proj_x.transpose();
  frame_replica_ids_.push_back(replica_id);
  frame_timestamps_.push_back(0.0);  // overwritten by record_frame_time
}

void TitanMSM::add_trajectory_break() {
  traj_break_indices_.push_back(static_cast<int>(n_proj_frames_));
}

void TitanMSM::clear_projected() {
  n_proj_frames_ = 0;
  // proj_mat_ capacity is kept to avoid re-allocation on reuse
  traj_break_indices_.clear();
  discrete_traj_.clear();
  clustered_ = false;
  msm_built_ = false;
  eigen_solved_ = false;
  frame_replica_ids_.clear();
  frame_timestamps_.clear();
}

void TitanMSM::record_frame_time(double time_ps, int /*replica_id*/) {
  if (!frame_timestamps_.empty())
    frame_timestamps_.back() = time_ps;
}

std::unordered_map<int, std::vector<int>>
TitanMSM::get_per_replica_dtrajs() const {
  std::unordered_map<int, std::vector<int>> result;
  for (size_t i = 0; i < discrete_traj_.size(); ++i) {
    int rid = (i < frame_replica_ids_.size())
              ? frame_replica_ids_[i] : -1;
    result[rid].push_back(discrete_traj_[i]);
  }
  return result;
}

std::unordered_map<int, std::vector<double>>
TitanMSM::get_per_replica_timestamps() const {
  std::unordered_map<int, std::vector<double>> result;
  for (size_t i = 0; i < frame_timestamps_.size(); ++i) {
    int rid = (i < frame_replica_ids_.size())
              ? frame_replica_ids_[i] : -1;
    result[rid].push_back(frame_timestamps_[i]);
  }
  return result;
}

// ============================================================================
// train_kmeans — K-Means++ init + Lloyd iteration  (Fix 3: RowMajor proj_mat_)
// n_threads: number of std::threads for the assignment step
// ============================================================================
void TitanMSM::train_kmeans(int max_iter, unsigned int seed, int n_threads) {
  n_threads_ = std::max(1, n_threads);
  // --- Preconditions -------------------------------------------------------
  if (n_proj_frames_ == 0) {
    TitanUI::print_error(
        "TitanMSM::train_kmeans: no projected frames collected. "
        "Call add_projected_frame() first.");
    return;
  }
  if (n_proj_frames_ < static_cast<int64_t>(n_k_)) {
    TitanUI::print_error(
        "TitanMSM::train_kmeans: fewer frames than clusters (" +
        std::to_string(n_proj_frames_) + " < " + std::to_string(n_k_) +
        "). Reduce msm_clusters.");
    return;
  }

  TitanUI::print_module_header("K-Means Clustering (MSM Phase 2A)");

  const int64_t N = n_proj_frames_;
  std::mt19937 rng(seed);

  // -------------------------------------------------------------------------
  // K-Means++ Initialization  (single-threaded, touches whole proj_mat_)
  // -------------------------------------------------------------------------
  {
    std::uniform_int_distribution<int64_t> uni(0, N - 1);
    centroids_.row(0) = proj_mat_.row(uni(rng));
  }

  // -------------------------------------------------------------------------
  // K-Means++ Initialization — parallel weight computation per centroid step
  // -------------------------------------------------------------------------
  std::vector<double> weights(N, 0.0);

  auto compute_weights_parallel = [&](int up_to_c) {
    const int64_t chunk = (N + n_threads_ - 1) / n_threads_;
    std::vector<std::jthread> wk;
    wk.reserve(n_threads_);

    for (int t = 0; t < n_threads_; ++t) {
      int64_t start = t * chunk;
      int64_t end = std::min(start + chunk, N);
      if (start >= end)
        continue;

      wk.emplace_back([&, start, end, up_to_c]() {
        Eigen::VectorXd d_local(n_comp_);
        for (int64_t i = start; i < end; ++i) {
          double min_d2 = std::numeric_limits<double>::max();
          for (int j = 0; j < up_to_c; ++j) {
            d_local.noalias() =
                proj_mat_.row(i).transpose() - centroids_.row(j).transpose();
            double d2 = d_local.squaredNorm();
            if (d2 < min_d2)
              min_d2 = d2;
          }
          weights[static_cast<size_t>(i)] = min_d2;
        }
      });
    }
    // jthread RAII join
  };

  for (int c = 1; c < n_k_; ++c) {
    compute_weights_parallel(c);
    std::discrete_distribution<int64_t> weighted(weights.begin(),
                                                 weights.end());
    centroids_.row(c) = proj_mat_.row(weighted(rng));
  }

  TitanUI::print_info("  K-Means++ initialization complete (" +
                      std::to_string(n_k_) + " centroids).");

  // -------------------------------------------------------------------------
  // K-Means Lloyd Iteration  (multi-threaded assignment via std::thread)
  // -------------------------------------------------------------------------
  discrete_traj_.assign(static_cast<size_t>(N), 0);
  Eigen::MatrixXd new_centroids(n_k_, n_comp_);
  std::vector<int> counts(n_k_);

  // Per-thread partial accumulators (to avoid mutex contention in update step)
  std::vector<Eigen::MatrixXd> thread_sums(
      n_threads_, Eigen::MatrixXd::Zero(n_k_, n_comp_));
  std::vector<std::vector<int>> thread_counts(n_threads_,
                                              std::vector<int>(n_k_, 0));

  double convergence = std::numeric_limits<double>::max();
  int iter = 0;

  for (iter = 0; iter < max_iter && convergence > 1e-6; ++iter) {

    // ----- Assignment step (parallel over frames) -------------------------
    auto assign_worker = [&](int tid, int64_t row_start, int64_t row_end) {
      Eigen::VectorXd d(n_comp_); // thread-local, no heap per frame
      for (int64_t i = row_start; i < row_end; ++i) {
        double best_d2 = std::numeric_limits<double>::max();
        int best_c = 0;
        for (int c = 0; c < n_k_; ++c) {
          d.noalias() =
              proj_mat_.row(i).transpose() - centroids_.row(c).transpose();
          double d2 = d.squaredNorm();
          if (d2 < best_d2) {
            best_d2 = d2;
            best_c = c;
          }
        }
        discrete_traj_[static_cast<size_t>(i)] = best_c;
      }
    };

    {
      std::vector<std::thread> workers;
      workers.reserve(n_threads_);
      int64_t chunk = (N + n_threads_ - 1) / n_threads_;
      for (int t = 0; t < n_threads_; ++t) {
        int64_t start = t * chunk;
        int64_t end = std::min(start + chunk, N);
        if (start < end)
          workers.emplace_back(assign_worker, t, start, end);
      }
      for (auto &w : workers)
        w.join();
    }

    // ----- Update step (parallel partial sums, serial merge) --------------
    for (int t = 0; t < n_threads_; ++t) {
      thread_sums[t].setZero();
      std::fill(thread_counts[t].begin(), thread_counts[t].end(), 0);
    }

    auto update_worker = [&](int tid, int64_t row_start, int64_t row_end) {
      for (int64_t i = row_start; i < row_end; ++i) {
        int c = discrete_traj_[static_cast<size_t>(i)];
        thread_sums[tid].row(c) += proj_mat_.row(i);
        thread_counts[tid][c]++;
      }
    };
    {
      std::vector<std::thread> workers;
      workers.reserve(n_threads_);
      int64_t chunk = (N + n_threads_ - 1) / n_threads_;
      for (int t = 0; t < n_threads_; ++t) {
        int64_t start = t * chunk;
        int64_t end = std::min(start + chunk, N);
        if (start < end)
          workers.emplace_back(update_worker, t, start, end);
      }
      for (auto &w : workers)
        w.join();
    }

    // Merge partial sums
    new_centroids.setZero();
    std::fill(counts.begin(), counts.end(), 0);
    for (int t = 0; t < n_threads_; ++t) {
      new_centroids += thread_sums[t];
      for (int c = 0; c < n_k_; ++c)
        counts[c] += thread_counts[t][c];
    }

    for (int c = 0; c < n_k_; ++c) {
      if (counts[c] > 0) {
        new_centroids.row(c) /= static_cast<double>(counts[c]);
      } else {
        // Empty cluster: reinit to farthest point (single-thread scan)
        TitanUI::print_info("  [WARNING] K-Means: cluster " +
                            std::to_string(c) +
                            " is empty. Reinitializing to farthest point.");
        double max_d2 = -1.0;
        int far_idx = 0;
        Eigen::VectorXd fd(n_comp_);
        for (int64_t i = 0; i < N; ++i) {
          int near_c = discrete_traj_[static_cast<size_t>(i)];
          fd.noalias() =
              proj_mat_.row(i).transpose() - centroids_.row(near_c).transpose();
          double d2 = fd.squaredNorm();
          if (d2 > max_d2) {
            max_d2 = d2;
            far_idx = static_cast<int>(i);
          }
        }
        new_centroids.row(c) = proj_mat_.row(far_idx);
        counts[c] = 1;
      }
    }

    // Convergence check
    convergence = 0.0;
    for (int c = 0; c < n_k_; ++c) {
      double move = (new_centroids.row(c) - centroids_.row(c)).norm();
      if (move > convergence)
        convergence = move;
    }
    centroids_ = new_centroids;
  }

  // Final assignment
  {
    auto assign_worker = [&](int tid, int64_t row_start, int64_t row_end) {
      Eigen::VectorXd d(n_comp_);
      for (int64_t i = row_start; i < row_end; ++i) {
        double best_d2 = std::numeric_limits<double>::max();
        int best_c = 0;
        for (int c = 0; c < n_k_; ++c) {
          d.noalias() =
              proj_mat_.row(i).transpose() - centroids_.row(c).transpose();
          double d2 = d.squaredNorm();
          if (d2 < best_d2) {
            best_d2 = d2;
            best_c = c;
          }
        }
        discrete_traj_[static_cast<size_t>(i)] = best_c;
      }
    };
    std::vector<std::thread> workers;
    int64_t chunk = (N + n_threads_ - 1) / n_threads_;
    for (int t = 0; t < n_threads_; ++t) {
      int64_t start = t * chunk;
      int64_t end = std::min(start + chunk, N);
      if (start < end)
        workers.emplace_back(assign_worker, t, start, end);
    }
    for (auto &w : workers)
      w.join();
  }

  // Inertia
  double inertia = 0.0;
  {
    Eigen::VectorXd d(n_comp_);
    for (int64_t i = 0; i < N; ++i) {
      d.noalias() =
          proj_mat_.row(i).transpose() -
          centroids_.row(discrete_traj_[static_cast<size_t>(i)]).transpose();
      inertia += d.squaredNorm();
    }
  }

  clustered_ = true;
  TitanUI::print_info("  Converged in " + std::to_string(iter) +
                      " iterations (max centroid shift = " +
                      std::to_string(convergence) + ").");
  TitanUI::print_info("  K-Means inertia: " + std::to_string(inertia));
}

// ============================================================================
// build_markov_model — boundary-aware count matrix → TPM → stationary dist
// ============================================================================
void TitanMSM::build_markov_model(int msm_lag_time) {
  // --- Preconditions -------------------------------------------------------
  if (!clustered_) {
    TitanUI::print_error(
        "TitanMSM::build_markov_model: train_kmeans() must run first.");
    return;
  }
  if (msm_lag_time <= 0) {
    TitanUI::print_error(
        "TitanMSM::build_markov_model: msm_lag_time must be > 0.");
    return;
  }
  const int N_frames = static_cast<int>(n_proj_frames_);
  if (N_frames < msm_lag_time * 10) {
    TitanUI::print_info(
        "[WARNING] TitanMSM: trajectory length (" + std::to_string(N_frames) +
        " frames) < 10 * lag_time (" + std::to_string(msm_lag_time * 10) +
        "). MSM statistics may be unreliable.");
  }

  msm_lag_time_ = msm_lag_time;

  // -------------------------------------------------------------------------
  // Step 1 — Count matrix (boundary-aware)
  // Build segments from traj_break_indices_. The last segment ends at N_frames.
  // -------------------------------------------------------------------------
  count_matrix_ = Eigen::MatrixXd::Zero(n_k_, n_k_);

  // Build segment boundary list: [0, b0, b1, ..., N_frames]
  std::vector<int> boundaries;
  boundaries.push_back(0);
  for (int b : traj_break_indices_) {
    if (b > 0 && b < N_frames)
      boundaries.push_back(b);
  }
  // If no breaks we declared, or all frames form one segment,
  // the last boundary is always N_frames.
  if (boundaries.back() != N_frames)
    boundaries.push_back(N_frames);

  for (std::size_t s = 0; s + 1 < boundaries.size(); ++s) {
    const int seg_start = boundaries[s];
    const int seg_end = boundaries[s + 1]; // exclusive
    for (int i = seg_start; i < seg_end - msm_lag_time; ++i) {
      count_matrix_(discrete_traj_[i], discrete_traj_[i + msm_lag_time]) += 1.0;
    }
  }

  // -------------------------------------------------------------------------
  // Step 2 — Symmetrize COUNT matrix (enforce detailed balance)
  // -------------------------------------------------------------------------
  Eigen::MatrixXd C_sym = 0.5 * (count_matrix_ + count_matrix_.transpose());

  // -------------------------------------------------------------------------
  // Step 3 — Row-normalize → row-stochastic TPM
  // -------------------------------------------------------------------------
  tpm_ = Eigen::MatrixXd::Zero(n_k_, n_k_);
  for (int i = 0; i < n_k_; ++i) {
    double row_sum = C_sym.row(i).sum();
    if (row_sum < 1e-10) {
      TitanUI::print_info(
          "[WARNING] TitanMSM: microstate " + std::to_string(i) +
          " is unvisited (zero row in count matrix). Setting uniform row. "
          "Consider pruning unvisited states.");
      tpm_.row(i) = Eigen::VectorXd::Constant(n_k_, 1.0 / n_k_).transpose();
    } else {
      tpm_.row(i) = C_sym.row(i) / row_sum;
    }
  }

  // -------------------------------------------------------------------------
  // Step 4 — Stationary distribution via power iteration
  // pi (1×k) = pi * T  until convergence (max 1000 steps)
  // Initialize pi as uniform.
  // -------------------------------------------------------------------------
  Eigen::RowVectorXd pi = Eigen::RowVectorXd::Constant(n_k_, 1.0 / n_k_);
  Eigen::RowVectorXd pi_new(n_k_);

  for (int it = 0; it < 1000; ++it) {
    pi_new.noalias() = pi * tpm_;
    double max_change = (pi_new - pi).cwiseAbs().maxCoeff();
    pi = pi_new;
    if (max_change < 1e-10)
      break;
  }

  // Normalize
  double pi_sum = pi.sum();
  if (pi_sum > 1e-30)
    pi /= pi_sum;

  stationary_dist_ = pi.transpose();
  msm_built_ = true;

  // Fix 4: cache eigendecomposition NOW (once, used by ITS + PCCA)
  {
    Eigen::EigenSolver<Eigen::MatrixXd> es(tpm_, /*computeEigenvectors=*/true);
    if (es.info() == Eigen::Success) {
      cached_eigenvalues_ = es.eigenvalues();
      cached_eigenvectors_ = es.eigenvectors();
      eigen_solved_ = true;
    } else {
      TitanUI::print_info(
          "[WARNING] TitanMSM: EigenSolver failed to converge — "
          "ITS and PCCA will be unavailable.");
      eigen_solved_ = false;
    }
  }
}

// ============================================================================
// compute_implied_timescales — eigendecomposition of TPM
// ============================================================================
std::vector<double>
TitanMSM::compute_implied_timescales(int n_timescales) const {
  if (!msm_built_) {
    TitanUI::print_error("TitanMSM::compute_implied_timescales: "
                         "build_markov_model() must run first.");
    return {};
  }
  if (!eigen_solved_)
    return {};

  // Fix 4: use cached eigenvalues from build_markov_model
  const auto &evals = cached_eigenvalues_;

  // Sort indices by descending |eigenvalue|
  const int k = n_k_;
  std::vector<int> idx(k);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int a, int b) {
    return std::abs(evals[a]) > std::abs(evals[b]);
  });

  std::vector<double> timescales;
  timescales.reserve(n_timescales);

  for (int rank = 1;
       rank < k && static_cast<int>(timescales.size()) < n_timescales; ++rank) {
    int ei = idx[rank];
    double abs_lambda = std::abs(evals[ei]);
    double imag_part = std::abs(evals[ei].imag());

    if (abs_lambda >= 1.0) {
      continue;
    }
    if (imag_part > 1e-6) {
      continue;
    }
    if (abs_lambda < 1e-15) {
      continue;
    }

    double t_i =
        -static_cast<double>(msm_lag_time_) * dt_ns_ / std::log(abs_lambda);
    timescales.push_back(t_i);
  }

  return timescales;
}

// ============================================================================
// generate_its_lags — geometric sequence for ITS scan
// ============================================================================
std::vector<int> TitanMSM::generate_its_lags() const {
  const int max_usable = static_cast<int>(n_proj_frames_) / 10;
  if (max_usable < 5) {
    TitanUI::print_info(
        "[WARNING] TitanMSM: trajectory too short for reliable ITS scan "
        "(max_usable_lag=" +
        std::to_string(max_usable) +
        " < 5). "
        "Skipping ITS scan.");
    return {};
  }

  // Geometric sequence: 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, ...
  static const int geo_steps[] = {1,    2,    5,     10,    20,
                                  50,   100,  200,   500,   1000,
                                  2000, 5000, 10000, 20000, 50000};
  std::vector<int> lags;
  for (int lag : geo_steps) {
    if (lag > max_usable)
      break;
    lags.push_back(lag);
  }
  return lags;
}

// ============================================================================
// run_pcca_crisp — argmax of top n_macro eigenvectors per microstate
// ============================================================================
std::vector<int> TitanMSM::run_pcca_crisp(int n_macro) const {
  // Fix 4: use cached eigenvectors from build_markov_model
  if (!eigen_solved_) {
    TitanUI::print_info("[WARNING] PCCA: no cached eigenvalues. Skipping.");
    return std::vector<int>(n_k_, 0);
  }

  // Sort by descending |eigenvalue|
  const int k = n_k_;
  std::vector<int> idx(k);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int a, int b) {
    return std::abs(cached_eigenvalues_[a]) > std::abs(cached_eigenvalues_[b]);
  });

  // Build V: top n_macro right eigenvectors (real parts)
  Eigen::MatrixXd V(k, n_macro);
  for (int j = 0; j < n_macro; ++j) {
    V.col(j) = cached_eigenvectors_.col(idx[j]).real();
  }

  // Crisp PCCA: macrostate[i] = argmax_j V(i, j)
  std::vector<int> macrostates(k, 0);
  for (int i = 0; i < k; ++i) {
    Eigen::Index best_j;
    V.row(i).maxCoeff(&best_j);
    macrostates[i] = static_cast<int>(best_j);
  }

  return macrostates;
}

// ============================================================================
// auto_detect_macrostates
// ============================================================================
int TitanMSM::auto_detect_macrostates(
    const std::vector<double> &timescales) const {
  double threshold = 2.0 * static_cast<double>(msm_lag_time_) * dt_ns_;
  int n_slow = 0;
  for (double t : timescales) {
    if (t > threshold)
      ++n_slow;
  }
  // n_macro = n_slow_processes + 1 (includes ground state)
  int n_macro = n_slow + 1;
  // Clamp: at least 2, at most k/2
  n_macro = std::clamp(n_macro, 2, std::max(2, n_k_ / 2));
  return n_macro;
}

// ============================================================================
// export_results — write all output files + run ITS scan + PCCA
// ============================================================================
void TitanMSM::export_results(const info &info_output, std::fstream &log,
                              int n_timescales, bool msm_run_its,
                              bool msm_run_pcca, int msm_macrostates,
                              const std::vector<int> &its_lags) const {
  if (!msm_built_) {
    TitanUI::print_error(
        "TitanMSM::export_results: build_markov_model() must run first.");
    return;
  }

  TitanUI::print_module_header("MSM Export Results (Phase 2D)");

  // =========================================================================
  // 1. Compute reference implied timescales (at the stored lag)
  // =========================================================================
  auto ref_timescales = compute_implied_timescales(n_timescales);

  // =========================================================================
  // 2. Write TPM
  // =========================================================================
  {
    std::ofstream f(info_output.get_msm_tpm_filename());
    if (!f.is_open()) {
      TitanUI::print_error("Cannot open " + info_output.get_msm_tpm_filename());
    } else {
      f << "# k=" << n_k_ << "  lag=" << msm_lag_time_ << "  dt=" << dt_ns_
        << " ns\n";
      f << std::fixed;
      for (int i = 0; i < n_k_; ++i) {
        for (int j = 0; j < n_k_; ++j) {
          if (j > 0)
            f << ' ';
          f << tpm_(i, j);
        }
        f << '\n';
      }
      TitanUI::print_info("  Written: " + info_output.get_msm_tpm_filename());
    }
  }

  // =========================================================================
  // 2b. Write raw count matrix (before symmetrization, for PyEMMA/msmtools MLE)
  // =========================================================================
  {
    std::ofstream f(info_output.get_msm_count_matrix_filename());
    if (!f.is_open()) {
      TitanUI::print_error("Cannot open " +
                           info_output.get_msm_count_matrix_filename());
    } else {
      f << "# Raw count matrix k=" << n_k_ << " lag=" << msm_lag_time_ << "\n";
      f << std::fixed << std::setprecision(1);
      for (int i = 0; i < n_k_; ++i) {
        for (int j = 0; j < n_k_; ++j) {
          if (j > 0)
            f << ' ';
          f << count_matrix_(i, j);
        }
        f << '\n';
      }
      TitanUI::print_info("  Written: " +
                          info_output.get_msm_count_matrix_filename());
    }
  }

  // =========================================================================
  // 2c. Write stationary distribution π (for deeptime/msmtools weighted PCCA+)
  // =========================================================================
  {
    std::ofstream f(info_output.get_msm_stationary_filename());
    if (!f.is_open()) {
      TitanUI::print_error("Cannot open " +
                           info_output.get_msm_stationary_filename());
    } else {
      f << "# microstate  stationary_prob\n";
      f << std::scientific << std::setprecision(8);
      for (int i = 0; i < n_k_; ++i)
        f << i << "  " << stationary_dist_[i] << "\n";
      TitanUI::print_info("  Written: " +
                          info_output.get_msm_stationary_filename());
    }
  }

  // =========================================================================
  // 3. Write discrete trajectory
  // =========================================================================
  {
    std::ofstream f(info_output.get_msm_discrete_traj_filename());
    if (!f.is_open()) {
      TitanUI::print_error("Cannot open " +
                           info_output.get_msm_discrete_traj_filename());
    } else {
      for (int s : discrete_traj_)
        f << s << '\n';
      TitanUI::print_info("  Written: " +
                          info_output.get_msm_discrete_traj_filename());
    }
  }

  // =========================================================================
  // 4. Write cluster centroids
  // =========================================================================
  {
    std::ofstream f(info_output.get_msm_centroids_filename());
    if (!f.is_open()) {
      TitanUI::print_error("Cannot open " +
                           info_output.get_msm_centroids_filename());
    } else {
      f << "# k=" << n_k_ << "  n_tica_components=" << n_comp_ << '\n';
      f << std::fixed;
      for (int i = 0; i < n_k_; ++i) {
        for (int j = 0; j < n_comp_; ++j) {
          if (j > 0)
            f << ' ';
          f << centroids_(i, j);
        }
        f << '\n';
      }
      TitanUI::print_info("  Written: " +
                          info_output.get_msm_centroids_filename());
    }
  }

  // =========================================================================
  // 5. Write implied timescales (reference lag)
  // =========================================================================
  {
    std::ofstream f(info_output.get_msm_timescales_filename());
    if (!f.is_open()) {
      TitanUI::print_error("Cannot open " +
                           info_output.get_msm_timescales_filename());
    } else {
      // Fix 4: use cached eigenvalues from build_markov_model for the timescale
      // table
      if (!eigen_solved_) {
        TitanUI::print_info(
            "[WARNING] Eigenvalues not cached; skipping table.");
      } else {
        const auto &evals = cached_eigenvalues_;
        std::vector<int> idx(n_k_);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
          return std::abs(evals[a]) > std::abs(evals[b]);
        });

        int written = 0;
        for (int rank = 1; rank < n_k_ && written < n_timescales; ++rank) {
          int ei = idx[rank];
          double abs_lambda = std::abs(evals[ei]);
          double imag_part = std::abs(evals[ei].imag());
          if (abs_lambda >= 1.0 || imag_part > 1e-6 || abs_lambda < 1e-15)
            continue;
          double t_i = -static_cast<double>(msm_lag_time_) * dt_ns_ /
                       std::log(abs_lambda);
          f << std::fixed << written << "  " << t_i << "  " << abs_lambda
            << '\n';
          ++written;
        }
      }
      TitanUI::print_info("  Written: " +
                          info_output.get_msm_timescales_filename());
    }
  }

  // =========================================================================
  // 6. ITS scan (multi-lag)
  // =========================================================================
  if (msm_run_its) {
    std::vector<int> lags_to_use =
        its_lags.empty() ? generate_its_lags() : its_lags;

    if (!lags_to_use.empty()) {
      std::ofstream f(info_output.get_msm_its_filename());
      if (!f.is_open()) {
        TitanUI::print_error("Cannot open " +
                             info_output.get_msm_its_filename());
      } else {
        f << "# lag_frames  lag_ns";
        for (int t = 0; t < n_timescales; ++t)
          f << "  t" << (t + 1) << "_ns";
        f << '\n';

        // Per-lag result struct (computed by workers, written in order)
        struct LagResult {
          int lag = 0;
          bool ok = false;
          Eigen::MatrixXd T_tmp;
          std::vector<std::string> its_values;
        };

        const int n_lags = static_cast<int>(lags_to_use.size());
        std::vector<LagResult> results(n_lags);

        // Pre-build segment boundaries (read-only, shared)
        const int N_frames = static_cast<int>(n_proj_frames_);
        std::vector<int> boundaries;
        boundaries.push_back(0);
        for (int b : traj_break_indices_) {
          if (b > 0 && b < N_frames)
            boundaries.push_back(b);
        }
        if (boundaries.back() != N_frames)
          boundaries.push_back(N_frames);

        // Cap worker count at n_lags
        const int n_workers = std::min(n_lags, n_threads_);

        // Launch jthread workers — each processes a shard of lags
        {
          std::vector<std::jthread> workers;
          workers.reserve(n_workers);

          for (int w = 0; w < n_workers; ++w) {
            int shard_begin = n_lags * w / n_workers;
            int shard_end = n_lags * (w + 1) / n_workers;
            if (shard_begin >= shard_end)
              continue;

            workers.emplace_back([&, shard_begin, shard_end]() {
              Eigen::MatrixXd C_tmp(n_k_, n_k_);
              Eigen::MatrixXd T_tmp_local(n_k_, n_k_);

              for (int li = shard_begin; li < shard_end; ++li) {
                int lag = lags_to_use[static_cast<size_t>(li)];
                LagResult &res = results[static_cast<size_t>(li)];
                res.lag = lag;

                // Build count matrix for this lag
                C_tmp.setZero();
                for (std::size_t s = 0; s + 1 < boundaries.size(); ++s) {
                  const int seg_start = boundaries[s];
                  const int seg_end = boundaries[s + 1];
                  for (int i = seg_start; i < seg_end - lag; ++i) {
                    C_tmp(discrete_traj_[i], discrete_traj_[i + lag]) += 1.0;
                  }
                }

                // Symmetrize + row-normalize
                Eigen::MatrixXd C_sym = 0.5 * (C_tmp + C_tmp.transpose());
                T_tmp_local.setZero();
                for (int i = 0; i < n_k_; ++i) {
                  double rs = C_sym.row(i).sum();
                  if (rs < 1e-10)
                    T_tmp_local.row(i) =
                        Eigen::VectorXd::Constant(n_k_, 1.0 / n_k_).transpose();
                  else
                    T_tmp_local.row(i) = C_sym.row(i) / rs;
                }

                // EigenSolver (thread-local instance)
                Eigen::EigenSolver<Eigen::MatrixXd> es_tmp(T_tmp_local);
                if (es_tmp.info() != Eigen::Success) {
                  res.its_values.assign(n_timescales, "  nan");
                  res.T_tmp = T_tmp_local;
                  res.ok = true;
                  continue;
                }
                auto evals_tmp = es_tmp.eigenvalues();

                std::vector<int> idx_tmp(n_k_);
                std::iota(idx_tmp.begin(), idx_tmp.end(), 0);
                std::sort(idx_tmp.begin(), idx_tmp.end(), [&](int a, int b) {
                  return std::abs(evals_tmp[a]) > std::abs(evals_tmp[b]);
                });

                int written = 0;
                for (int rank = 1; rank < n_k_ && written < n_timescales;
                     ++rank) {
                  int ei = idx_tmp[rank];
                  double abs_l = std::abs(evals_tmp[ei]);
                  double im_l = std::abs(evals_tmp[ei].imag());
                  if (abs_l >= 1.0 || im_l > 1e-6 || abs_l < 1e-15) {
                    res.its_values.push_back("  nan");
                  } else {
                    double t_i =
                        -static_cast<double>(lag) * dt_ns_ / std::log(abs_l);
                    std::ostringstream oss;
                    oss << std::fixed << "  " << t_i;
                    res.its_values.push_back(oss.str());
                  }
                  ++written;
                }
                while (written < n_timescales) {
                  res.its_values.push_back("  nan");
                  ++written;
                }

                res.T_tmp = T_tmp_local;
                res.ok = true;
              }
            }); // end jthread lambda
          }
          // jthread RAII: all join here
        }

        // Write results IN ORDER (single-threaded, after all workers joined)
        std::filesystem::create_directories(info_output.get_msm_its_tpm_dir());

        for (int li = 0; li < n_lags; ++li) {
          const LagResult &res = results[static_cast<size_t>(li)];
          if (!res.ok)
            continue;

          int lag = res.lag;

          // Write per-lag TPM file
          {
            std::string tpm_path = info_output.get_msm_its_tpm_dir() +
                                   "/tpm_lag" + std::to_string(lag) + ".dat";
            std::ofstream ft(tpm_path);
            if (ft.is_open()) {
              ft << "# TPM lag=" << lag << "  dt=" << dt_ns_
                 << " ns  k=" << n_k_ << "\n";
              ft << std::fixed;
              for (int i = 0; i < n_k_; ++i) {
                for (int j = 0; j < n_k_; ++j) {
                  if (j > 0)
                    ft << ' ';
                  ft << res.T_tmp(i, j);
                }
                ft << '\n';
              }
            }
          }

          // Write ITS scan row
          double lag_ns = static_cast<double>(lag) * dt_ns_;
          f << std::fixed << lag << "  " << lag_ns;
          for (const auto &s : res.its_values)
            f << s;
          f << '\n';
        }

        TitanUI::print_info("  Written: " + info_output.get_msm_its_filename());
        TitanUI::log_kv(log, "ITS Scan Lags",
                        std::to_string(lags_to_use.size()));
      }
    }
  } else {
    TitanUI::log_kv(log, "ITS Scan", "disabled");
  }

  // =========================================================================
  // 7. PCCA crisp macrostate lumping
  // =========================================================================
  int n_macro_used = 0;
  if (msm_run_pcca) {
    // Determine n_macro
    if (msm_macrostates == -1) {
      n_macro_used = auto_detect_macrostates(ref_timescales);
      TitanUI::log_kv(log, "Auto-detected macrostates",
                      std::to_string(n_macro_used));
    } else {
      n_macro_used = std::clamp(msm_macrostates, 2, n_k_);
    }

    auto macrostates = run_pcca_crisp(n_macro_used);

    // Write macrostate assignments
    std::ofstream f(info_output.get_msm_macrostates_filename());
    if (!f.is_open()) {
      TitanUI::print_error("Cannot open " +
                           info_output.get_msm_macrostates_filename());
    } else {
      f << "# microstate_id  macrostate_id\n";
      for (int i = 0; i < n_k_; ++i) {
        f << i << "  " << macrostates[i] << '\n';
      }
      TitanUI::print_info("  Written: " +
                          info_output.get_msm_macrostates_filename());
    }

    // Write macrostate trajectory (frame → macrostate, for Python scatter
    // plots)
    {
      std::ofstream fm(info_output.get_msm_macrostate_traj_filename());
      if (!fm.is_open()) {
        TitanUI::print_error("Cannot open " +
                             info_output.get_msm_macrostate_traj_filename());
      } else {
        fm << "# frame  microstate  macrostate\n";
        for (size_t i = 0; i < discrete_traj_.size(); ++i) {
          int micro = discrete_traj_[i];
          int macro = (micro < static_cast<int>(macrostates.size()))
                          ? macrostates[micro]
                          : -1;
          fm << i << "  " << micro << "  " << macro << "\n";
        }
        TitanUI::print_info("  Written: " +
                            info_output.get_msm_macrostate_traj_filename());
      }
    }

    TitanUI::log_kv(log, "PCCA Macrostates", std::to_string(n_macro_used));
  } else {
    TitanUI::log_kv(log, "PCCA Macrostates", "disabled");
  }

  // =========================================================================
  // 8. Log summary
  // =========================================================================
  TitanUI::log_kv(log, "MSM States", std::to_string(n_k_));
  TitanUI::log_kv(log, "MSM Lag", std::to_string(msm_lag_time_));

  if (ref_timescales.size() >= 1)
    TitanUI::log_kv(log, "Slowest Timescale (ns)",
                    std::to_string(ref_timescales[0]));
  if (ref_timescales.size() >= 2)
    TitanUI::log_kv(log, "Second Timescale (ns)",
                    std::to_string(ref_timescales[1]));
}
