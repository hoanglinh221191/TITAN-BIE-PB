#include "titan_kinetics.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

// ============================================================================
// Constructor — pre-allocate ALL matrices
// ============================================================================
TitanKinetics::TitanKinetics(int feature_dim, int lag_time, int n_components,
                             double dt_ns)
    : dim_(feature_dim), lag_(lag_time),
      n_comp_(std::min(n_components, feature_dim)), dt_ns_(dt_ns),
      sum_x_(Eigen::VectorXd::Zero(feature_dim)),
      sum_xx_(Eigen::MatrixXd::Zero(feature_dim, feature_dim)),
      sum_x_lagged_(Eigen::VectorXd::Zero(feature_dim)),
      sum_x_paired_(Eigen::VectorXd::Zero(feature_dim)),
      sum_xtau_(Eigen::MatrixXd::Zero(feature_dim, feature_dim)), n_frames_(0),
      n_pairs_(0), solved_(false), mean_(Eigen::VectorXd::Zero(feature_dim)),
      eigenvalues_(Eigen::VectorXd::Zero(n_comp_)),
      eigenvectors_(Eigen::MatrixXd::Zero(feature_dim, n_comp_)),
      timescales_(Eigen::VectorXd::Zero(n_comp_)),
      // Fix 2: pre-allocate circular lag buffer (lag_ × dim_), no heap in hot
      // loop
      circ_buf_(Eigen::MatrixXd::Zero(lag_time, feature_dim)), circ_head_(0),
      circ_count_(0) {}

// ============================================================================
// accumulate_frame — HOT LOOP (zero heap allocation)
// ============================================================================
void TitanKinetics::accumulate_frame(const Eigen::VectorXd &x) {
  // Fix 2 — write into circular buffer slot (zero heap allocation)
  int write_slot = (circ_head_ + circ_count_) % lag_;
  circ_buf_.row(write_slot) = x.transpose(); // Eigen memcpy, no alloc

  // Fix 1 — update only upper triangle of sum_xx_ (~2× fewer FMAs)
  sum_x_ += x;
  sum_xx_.selfadjointView<Eigen::Upper>().rankUpdate(x);
  ++n_frames_;

  if (circ_count_ == lag_) {
    // Buffer full: oldest entry is at circ_head_
    Eigen::VectorXd x_lagged = circ_buf_.row(circ_head_).transpose();

    sum_x_lagged_ += x_lagged;
    sum_x_paired_ += x;
    sum_xtau_.noalias() += x * x_lagged.transpose();
    ++n_pairs_;

    // Advance ring head (evict oldest)
    circ_head_ = (circ_head_ + 1) % lag_;
  } else {
    ++circ_count_; // still warming up
  }
}

// reset_lag_buffer — trajectory boundary: reset ring state (no dealloc)
void TitanKinetics::reset_lag_buffer() {
  circ_head_ = 0;
  circ_count_ = 0;
  // circ_buf_ data intentionally NOT zeroed — overwritten on next write
}

// ============================================================================
// solve_tica — center, build covariance matrices, solve eigenproblem
// ============================================================================
void TitanKinetics::solve_tica() {
  if (n_frames_ < 2) {
    TitanUI::print_error("TICA: insufficient frames (" +
                         std::to_string(n_frames_) + "). Need at least 2.");
    return;
  }
  if (n_pairs_ < 1) {
    TitanUI::print_error(
        "TICA: no time-lagged pairs formed. Increase trajectory "
        "length or decrease --tica-lag (current: " +
        std::to_string(lag_) + ").");
    return;
  }

  TitanUI::print_info("TICA: " + std::to_string(n_frames_) + " frames, " +
                      std::to_string(n_pairs_) + " lagged pairs, D=" +
                      std::to_string(dim_) + ", tau=" + std::to_string(lag_));

  // ---- Step 1: Compute means -----------------------------------------------
  mean_ = sum_x_ / static_cast<double>(n_frames_);
  Eigen::VectorXd mean_lagged = sum_x_lagged_ / static_cast<double>(n_pairs_);
  Eigen::VectorXd mean_paired = sum_x_paired_ / static_cast<double>(n_pairs_);

  // ---- Step 2: Build C0 (instantaneous covariance) -------------------------
  // Fix 1: sum_xx_ only has upper triangle filled — mirror it first
  Eigen::MatrixXd C0 =
      Eigen::MatrixXd(sum_xx_.selfadjointView<Eigen::Upper>()) /
          static_cast<double>(n_frames_) -
      mean_ * mean_.transpose();

  // ---- Step 3: Build Ctau (time-lagged covariance) -------------------------
  Eigen::MatrixXd Ctau = sum_xtau_ / static_cast<double>(n_pairs_) -
                         mean_paired * mean_lagged.transpose();

  // ---- Step 4: Symmetrize --------------------------------------------------
  C0 = 0.5 * (C0 + C0.transpose());
  Ctau = 0.5 * (Ctau + Ctau.transpose());

  // ---- Step 5: Regularize C0 -----------------------------------------------
  double epsilon = 1e-6 * C0.trace() / dim_;
  C0.diagonal().array() += epsilon;

  TitanUI::print_info("TICA: C0 trace = " + std::to_string(C0.trace()) +
                      ", regularization eps = " + std::to_string(epsilon));

  // ---- Step 6: Solve generalized eigenvalue problem Ctau * v = lambda * C0 *
  // v
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(Ctau, C0);
  if (solver.info() != Eigen::Success) {
    TitanUI::print_error(
        "TICA: GeneralizedSelfAdjointEigenSolver failed to converge.");
    return;
  }

  // ---- Step 7: Sort by DESCENDING eigenvalue (largest = slowest) -----------
  Eigen::VectorXd all_evals = solver.eigenvalues();
  Eigen::MatrixXd all_evecs = solver.eigenvectors();

  // Create index array and sort descending
  std::vector<int> idx(all_evals.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&all_evals](int a, int b) { return all_evals[a] > all_evals[b]; });

  // ---- Step 8: Retain top n_components -------------------------------------
  int actual_comp = std::min(n_comp_, static_cast<int>(all_evals.size()));
  eigenvalues_.resize(actual_comp);
  eigenvectors_.resize(dim_, actual_comp);
  timescales_.resize(actual_comp);

  for (int i = 0; i < actual_comp; ++i) {
    eigenvalues_[i] = all_evals[idx[i]];
    eigenvectors_.col(i) = all_evecs.col(idx[i]);

    // ---- Step 9: Implied timescales ----------------------------------------
    double abs_lambda = std::abs(eigenvalues_[i]);
    if (abs_lambda >= 1.0 || abs_lambda < 1e-15) {
      timescales_[i] = std::numeric_limits<double>::infinity();
    } else {
      timescales_[i] =
          -static_cast<double>(lag_) * dt_ns_ / std::log(abs_lambda);
    }
  }

  solved_ = true;
  TitanUI::print_info("TICA: solved successfully. Top " +
                      std::to_string(actual_comp) + " components retained.");
}

// ============================================================================
// project_frame — project a single frame onto the TICA subspace
// ============================================================================
Eigen::VectorXd TitanKinetics::project_frame(const Eigen::VectorXd &x) const {
  if (!solved_) {
    return Eigen::VectorXd::Zero(n_comp_);
  }
  return (x - mean_).transpose() * eigenvectors_;
}

// ============================================================================
// write_results — output eigenvalues, timescales, and eigenvectors
// ============================================================================
void TitanKinetics::write_results(const std::string &eval_file,
                                  const std::string &evec_file,
                                  std::fstream &log) const {
  if (!solved_) {
    TitanUI::print_warning("TICA: solve_tica() was not successful, skipping "
                           "output.");
    return;
  }

  // --- Eigenvalues + implied timescales ---
  {
    std::ofstream ofs(eval_file);
    if (!ofs.is_open()) {
      TitanUI::print_error("TICA: cannot open " + eval_file + " for writing.");
      return;
    }
    ofs << "# TICA Eigenvalues and Implied Timescales\n";
    ofs << "# lag_time = " << lag_ << " frames, dt = " << dt_ns_ << " ns\n";
    ofs << "# n_frames = " << n_frames_ << ", n_pairs = " << n_pairs_ << "\n";
    ofs << "# Component  Eigenvalue  Implied_Timescale(ns)\n";
    for (int i = 0; i < eigenvalues_.size(); ++i) {
      ofs << std::setw(4) << (i + 1) << "  " << std::fixed
          << std::setprecision(8) << eigenvalues_[i] << "  "
          << std::setprecision(4) << timescales_[i] << "\n";
    }
    ofs.close();
    TitanUI::print_info("TICA: eigenvalues written to " + eval_file);
    TitanUI::log_kv(log, "TICA Eigenvalues File", eval_file);
  }

  // --- Eigenvectors ---
  {
    std::ofstream ofs(evec_file);
    if (!ofs.is_open()) {
      TitanUI::print_error("TICA: cannot open " + evec_file + " for writing.");
      return;
    }
    ofs << "# TICA Eigenvectors (D=" << dim_ << " rows x "
        << eigenvectors_.cols() << " components)\n";
    ofs << "# Each column is one eigenvector, sorted by descending "
           "eigenvalue\n";
    for (int r = 0; r < eigenvectors_.rows(); ++r) {
      for (int c = 0; c < eigenvectors_.cols(); ++c) {
        if (c > 0)
          ofs << " ";
        ofs << std::fixed << std::setprecision(8) << eigenvectors_(r, c);
      }
      ofs << "\n";
    }
    ofs.close();
    TitanUI::print_info("TICA: eigenvectors written to " + evec_file);
    TitanUI::log_kv(log, "TICA Eigenvectors File", evec_file);
  }
}

// ============================================================================
// write_projections — write PRE-PROJECTED matrix rows as CSV for Python
// Input: proj_mat rows are ALREADY in TICA subspace (PC1..PCn per frame)
// DO NOT call project_frame() here — data is already projected.
// ============================================================================
void TitanKinetics::write_projections(const Eigen::MatrixXd &proj_mat,
                                      int64_t n_actual_frames,
                                      const std::string &proj_file,
                                      std::fstream &log) const {
  if (!solved_) {
    TitanUI::print_warning(
        "TICA: solve_tica() not successful — skipping projections.");
    return;
  }

  std::ofstream ofs(proj_file);
  if (!ofs.is_open()) {
    TitanUI::print_error("TICA: cannot open " + proj_file + " for writing.");
    return;
  }

  // Header: # frame  PC1  PC2  ...
  ofs << "# frame";
  const int n_cols = static_cast<int>(proj_mat.cols()); // actual stored cols
  for (int c = 0; c < n_cols; ++c)
    ofs << "  PC" << (c + 1);
  ofs << "\n";

  const int64_t N = n_actual_frames; // real frames, not allocated capacity
  ofs << std::fixed << std::setprecision(6);
  for (int64_t i = 0; i < N; ++i) {
    ofs << i;
    for (int c = 0; c < n_cols; ++c)
      ofs << "  " << proj_mat(i, c);
    ofs << "\n";
  }
  ofs.close();

  TitanUI::print_info("TICA: projections written to " + proj_file + " (" +
                      std::to_string(N) + " frames, " + std::to_string(n_cols) +
                      " components)");
  TitanUI::log_kv(log, "TICA Projections File", proj_file);
}
