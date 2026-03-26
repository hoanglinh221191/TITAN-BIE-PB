#include "titan_kinetics.h"

#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/SVD>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace {

using RowVectorXd = Eigen::Matrix<double, 1, Eigen::Dynamic>;

std::string lowercase_copy(const std::string &text) {
  std::string out = text;
  std::transform(out.begin(), out.end(), out.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return out;
}

double compute_tica_timescale(double eigenvalue, int lag, double dt_ns) {
  const double abs_lambda = std::abs(eigenvalue);
  if (abs_lambda >= 1.0 || abs_lambda < 1e-15) {
    return std::numeric_limits<double>::infinity();
  }
  return -static_cast<double>(lag) * dt_ns / std::log(abs_lambda);
}

int singular_rank(const Eigen::VectorXd &singular_values, int dim_hint) {
  if (singular_values.size() == 0) {
    return 0;
  }
  const double max_sv = singular_values[0];
  if (max_sv <= 0.0) {
    return 0;
  }
  const double tol = std::max(1, dim_hint) *
                     std::numeric_limits<double>::epsilon() * max_sv;
  int rank = 0;
  for (int i = 0; i < singular_values.size(); ++i) {
    if (singular_values[i] > tol) {
      ++rank;
    }
  }
  return rank;
}

} // namespace

std::string tica_backend_to_string(TicaBackend backend) {
  switch (backend) {
  case TicaBackend::Auto:
    return "auto";
  case TicaBackend::ExactDense:
    return "exact_dense";
  case TicaBackend::ExactDual:
    return "exact_dual";
  case TicaBackend::Randomized:
    return "randomized";
  }
  return "unknown";
}

bool tica_backend_from_string(const std::string &text, TicaBackend &backend) {
  const std::string value = lowercase_copy(text);
  if (value == "auto") {
    backend = TicaBackend::Auto;
    return true;
  }
  if (value == "exact_dense" || value == "dense") {
    backend = TicaBackend::ExactDense;
    return true;
  }
  if (value == "exact_dual" || value == "dual") {
    backend = TicaBackend::ExactDual;
    return true;
  }
  if (value == "randomized" || value == "rand") {
    backend = TicaBackend::Randomized;
    return true;
  }
  return false;
}

TitanKinetics::TitanKinetics(int feature_dim, int lag_time, int n_components,
                             double dt_ns, TicaBackend backend,
                             int randomized_rank, size_t feature_cache_cap_bytes,
                             bool allow_backend_promotion)
    : dim_(feature_dim), lag_(lag_time),
      n_comp_(std::max(1, std::min(n_components, feature_dim))), dt_ns_(dt_ns),
      requested_backend_(backend), backend_(backend),
      randomized_rank_(std::max(1, randomized_rank)),
      feature_cache_cap_bytes_(feature_cache_cap_bytes),
      allow_backend_promotion_(allow_backend_promotion),
      sum_x_(Eigen::VectorXd::Zero(feature_dim)),
      sum_x_lagged_(Eigen::VectorXd::Zero(feature_dim)),
      sum_x_paired_(Eigen::VectorXd::Zero(feature_dim)),
      n_frames_(0), n_pairs_(0),
      circ_buf_(FrameMatrixRM::Zero(std::max(1, lag_time), feature_dim)),
      solved_(false), mean_(Eigen::VectorXd::Zero(feature_dim)),
      eigenvalues_(Eigen::VectorXd::Zero(n_comp_)),
      eigenvectors_(Eigen::MatrixXd::Zero(feature_dim, n_comp_)),
      timescales_(Eigen::VectorXd::Constant(
          n_comp_, std::numeric_limits<double>::infinity())) {
  if (backend_ == TicaBackend::ExactDense) {
    sum_xx_ = Eigen::MatrixXd::Zero(feature_dim, feature_dim);
    sum_xtau_ = Eigen::MatrixXd::Zero(feature_dim, feature_dim);
  }
  if (backend_ == TicaBackend::Randomized) {
    open_spool_for_write_();
  }
}

TitanKinetics::~TitanKinetics() {
  close_spool_writer_();
  remove_spool_file_();
}

void TitanKinetics::open_spool_for_write_() {
  if (spool_writer_.is_open()) {
    return;
  }

  std::error_code ec;
  std::filesystem::path base_dir = std::filesystem::temp_directory_path(ec);
  if (ec || base_dir.empty()) {
    base_dir = std::filesystem::current_path(ec);
  }
  if (ec || base_dir.empty()) {
    throw std::runtime_error("TICA randomized: cannot determine spool directory.");
  }

  const auto stamp =
      std::chrono::steady_clock::now().time_since_epoch().count();
  spool_path_ =
      base_dir / ("titan_tica_spool_" + std::to_string(stamp) + ".bin");
  spool_writer_.open(spool_path_, std::ios::binary | std::ios::trunc);
  if (!spool_writer_.is_open()) {
    throw std::runtime_error("TICA randomized: cannot open spool file " +
                             spool_path_.string());
  }
}

void TitanKinetics::close_spool_writer_() {
  if (spool_writer_.is_open()) {
    spool_writer_.flush();
    spool_writer_.close();
  }
}

void TitanKinetics::remove_spool_file_() {
  if (spool_path_.empty()) {
    return;
  }
  std::error_code ec;
  std::filesystem::remove(spool_path_, ec);
  spool_path_.clear();
}

void TitanKinetics::append_cached_frame_(const Eigen::VectorXd &x) {
  if (backend_ == TicaBackend::Randomized) {
    if (!spool_writer_.is_open()) {
      open_spool_for_write_();
    }
    spool_writer_.write(reinterpret_cast<const char *>(x.data()),
                        static_cast<std::streamsize>(dim_ * sizeof(double)));
    if (!spool_writer_) {
      throw std::runtime_error("TICA randomized: failed writing spool file.");
    }
    return;
  }

  if (backend_ != TicaBackend::ExactDual) {
    return;
  }

  if (n_frames_ >= frame_cache_.rows()) {
    const int64_t new_rows =
        std::max<int64_t>(int64_t(1024), frame_cache_.rows() * 2);
    if (allow_backend_promotion_ && feature_cache_cap_bytes_ > 0) {
      const long double estimated_bytes =
          static_cast<long double>(new_rows) * static_cast<long double>(dim_) *
          static_cast<long double>(sizeof(double));
      if (estimated_bytes >
          static_cast<long double>(feature_cache_cap_bytes_)) {
        backend_ = TicaBackend::Randomized;
        open_spool_for_write_();
        for (int64_t i = 0; i < n_frames_; ++i) {
          spool_writer_.write(
              reinterpret_cast<const char *>(frame_cache_.row(i).data()),
              static_cast<std::streamsize>(dim_ * sizeof(double)));
        }
        frame_cache_.resize(0, dim_);
        TitanUI::print_warning(
            "TICA: feature cache exceeded auto memory budget. "
            "Promoting backend from exact_dual to randomized.");
        spool_writer_.write(reinterpret_cast<const char *>(x.data()),
                            static_cast<std::streamsize>(dim_ * sizeof(double)));
        if (!spool_writer_) {
          throw std::runtime_error(
              "TICA randomized: failed writing promoted spool file.");
        }
        return;
      }
    }
    frame_cache_.conservativeResize(new_rows, dim_);
  }

  frame_cache_.row(n_frames_) = x.transpose();
}

std::vector<int64_t> TitanKinetics::build_segment_boundaries_() const {
  std::vector<int64_t> boundaries;
  boundaries.reserve(traj_break_indices_.size() + 2);
  boundaries.push_back(0);
  for (int64_t b : traj_break_indices_) {
    if (b > 0 && b < n_frames_) {
      boundaries.push_back(b);
    }
  }
  if (boundaries.back() != n_frames_) {
    boundaries.push_back(n_frames_);
  }
  return boundaries;
}

void TitanKinetics::accumulate_frame(const Eigen::VectorXd &x) {
  if (x.size() != dim_) {
    throw std::runtime_error("TICA: feature dimension mismatch in accumulate_frame().");
  }

  append_cached_frame_(x);

  const int write_slot = (circ_head_ + circ_count_) % std::max(1, lag_);
  circ_buf_.row(write_slot) = x.transpose();

  sum_x_ += x;
  if (backend_ == TicaBackend::ExactDense) {
    sum_xx_.selfadjointView<Eigen::Upper>().rankUpdate(x);
  }
  ++n_frames_;

  if (circ_count_ == lag_) {
    Eigen::Map<const Eigen::VectorXd> x_lagged(circ_buf_.row(circ_head_).data(),
                                               dim_);
    sum_x_lagged_ += x_lagged;
    sum_x_paired_ += x;
    if (backend_ == TicaBackend::ExactDense) {
      sum_xtau_.noalias() += x * x_lagged.transpose();
    }
    ++n_pairs_;
    circ_head_ = (circ_head_ + 1) % lag_;
  } else {
    ++circ_count_;
  }
}

void TitanKinetics::reset_lag_buffer() {
  if (n_frames_ > 0 &&
      (traj_break_indices_.empty() || traj_break_indices_.back() != n_frames_)) {
    traj_break_indices_.push_back(n_frames_);
  }
  circ_head_ = 0;
  circ_count_ = 0;
}

void TitanKinetics::solve_exact_dense_() {
  const Eigen::VectorXd mean_lagged =
      sum_x_lagged_ / static_cast<double>(n_pairs_);
  const Eigen::VectorXd mean_paired =
      sum_x_paired_ / static_cast<double>(n_pairs_);

  Eigen::MatrixXd C0 =
      Eigen::MatrixXd(sum_xx_.selfadjointView<Eigen::Upper>()) /
          static_cast<double>(n_frames_) -
      mean_ * mean_.transpose();

  Eigen::MatrixXd Ctau =
      sum_xtau_ / static_cast<double>(n_pairs_) -
      mean_paired * mean_lagged.transpose();

  C0 = 0.5 * (C0 + C0.transpose());
  Ctau = 0.5 * (Ctau + Ctau.transpose());

  const double epsilon = 1e-6 * C0.trace() / std::max(1, dim_);
  C0.diagonal().array() += epsilon;

  TitanUI::print_info("TICA: C0 trace = " + std::to_string(C0.trace()) +
                      ", regularization eps = " + std::to_string(epsilon));

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(Ctau, C0);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error(
        "TICA: GeneralizedSelfAdjointEigenSolver failed to converge.");
  }

  const Eigen::VectorXd all_evals = solver.eigenvalues();
  const Eigen::MatrixXd all_evecs = solver.eigenvectors();

  std::vector<int> idx(all_evals.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&all_evals](int a, int b) { return all_evals[a] > all_evals[b]; });

  eigenvalues_.setZero();
  timescales_.setConstant(std::numeric_limits<double>::infinity());
  eigenvectors_.setZero();
  retained_components_ = std::min<int>(n_comp_, static_cast<int>(all_evals.size()));

  for (int i = 0; i < retained_components_; ++i) {
    eigenvalues_[i] = all_evals[idx[i]];
    eigenvectors_.col(i) = all_evecs.col(idx[i]);
    timescales_[i] =
        compute_tica_timescale(eigenvalues_[i], lag_, dt_ns_);
  }

  solved_ = true;
}

Eigen::MatrixXd TitanKinetics::build_exact_dual_basis_() {
  if (n_frames_ == 0) {
    return {};
  }

  Eigen::MatrixXd centered = frame_cache_.topRows(n_frames_);
  centered.rowwise() -= mean_.transpose();

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(centered,
                                        Eigen::ComputeThinU | Eigen::ComputeThinV);
  const int rank = singular_rank(svd.singularValues(), dim_);
  if (rank <= 0) {
    throw std::runtime_error("TICA exact_dual: centered feature matrix is rank-deficient.");
  }
  TitanUI::print_info("TICA exact_dual: retained PCA basis rank = " +
                      std::to_string(rank));
  return svd.matrixV().leftCols(rank);
}

Eigen::MatrixXd TitanKinetics::build_randomized_basis_() {
  close_spool_writer_();

  const int target_rank =
      std::min({dim_, static_cast<int>(n_frames_),
                std::max(n_comp_, randomized_rank_)});
  if (target_rank <= 0) {
    throw std::runtime_error("TICA randomized: invalid target rank.");
  }

  constexpr int oversample = 8;
  const int sketch_rank =
      std::min({dim_, static_cast<int>(n_frames_), target_rank + oversample});

  std::mt19937 rng(12345u);
  std::normal_distribution<double> normal(0.0, 1.0);
  Eigen::MatrixXd omega =
      Eigen::MatrixXd::NullaryExpr(dim_, sketch_rank,
                                   [&]() { return normal(rng); });

  RowVectorXd mean_row = mean_.transpose();
  RowVectorXd centered(dim_);
  Eigen::MatrixXd Y(n_frames_, sketch_rank);
  std::vector<double> row_buffer(static_cast<size_t>(dim_));

  {
    std::ifstream in(spool_path_, std::ios::binary);
    if (!in.is_open()) {
      throw std::runtime_error("TICA randomized: cannot open spool file for sketch pass.");
    }
    for (int64_t i = 0; i < n_frames_; ++i) {
      in.read(reinterpret_cast<char *>(row_buffer.data()),
              static_cast<std::streamsize>(dim_ * sizeof(double)));
      if (!in) {
        throw std::runtime_error("TICA randomized: short read in sketch pass.");
      }
      Eigen::Map<const RowVectorXd> row(row_buffer.data(), dim_);
      centered = row;
      centered -= mean_row;
      Y.row(i).noalias() = centered * omega;
    }
  }

  Eigen::HouseholderQR<Eigen::MatrixXd> qr(Y);
  Eigen::MatrixXd Q =
      qr.householderQ() * Eigen::MatrixXd::Identity(n_frames_, sketch_rank);

  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(sketch_rank, dim_);
  {
    std::ifstream in(spool_path_, std::ios::binary);
    if (!in.is_open()) {
      throw std::runtime_error(
          "TICA randomized: cannot open spool file for compression pass.");
    }
    for (int64_t i = 0; i < n_frames_; ++i) {
      in.read(reinterpret_cast<char *>(row_buffer.data()),
              static_cast<std::streamsize>(dim_ * sizeof(double)));
      if (!in) {
        throw std::runtime_error(
            "TICA randomized: short read in compression pass.");
      }
      Eigen::Map<const RowVectorXd> row(row_buffer.data(), dim_);
      centered = row;
      centered -= mean_row;
      B.noalias() += Q.row(i).transpose() * centered;
    }
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(B,
                                        Eigen::ComputeThinU | Eigen::ComputeThinV);
  const int available_rank = singular_rank(svd.singularValues(), dim_);
  const int basis_rank = std::min(target_rank, available_rank);
  if (basis_rank <= 0) {
    throw std::runtime_error("TICA randomized: randomized PCA returned zero rank.");
  }
  TitanUI::print_info("TICA randomized: sketch rank = " +
                      std::to_string(sketch_rank) +
                      ", retained PCA basis rank = " +
                      std::to_string(basis_rank));
  return svd.matrixV().leftCols(basis_rank);
}

Eigen::MatrixXd
TitanKinetics::project_all_frames_to_basis_(const Eigen::MatrixXd &basis) const {
  if (basis.cols() == 0) {
    return {};
  }

  Eigen::MatrixXd projected(n_frames_, basis.cols());
  RowVectorXd mean_row = mean_.transpose();
  RowVectorXd centered(dim_);

  if (backend_ == TicaBackend::ExactDual) {
    for (int64_t i = 0; i < n_frames_; ++i) {
      centered = frame_cache_.row(i);
      centered -= mean_row;
      projected.row(i).noalias() = centered * basis;
    }
    return projected;
  }

  std::ifstream in(spool_path_, std::ios::binary);
  if (!in.is_open()) {
    throw std::runtime_error("TICA randomized: cannot open spool file for projection pass.");
  }

  std::vector<double> row_buffer(static_cast<size_t>(dim_));
  for (int64_t i = 0; i < n_frames_; ++i) {
    in.read(reinterpret_cast<char *>(row_buffer.data()),
            static_cast<std::streamsize>(dim_ * sizeof(double)));
    if (!in) {
      throw std::runtime_error("TICA randomized: short read in projection pass.");
    }
    Eigen::Map<const RowVectorXd> row(row_buffer.data(), dim_);
    centered = row;
    centered -= mean_row;
    projected.row(i).noalias() = centered * basis;
  }

  return projected;
}

void TitanKinetics::solve_reduced_from_basis_(const Eigen::MatrixXd &basis,
                                              const std::string &label) {
  if (basis.cols() == 0) {
    throw std::runtime_error("TICA " + label + ": reduction basis is empty.");
  }

  const int reduced_dim = static_cast<int>(basis.cols());
  Eigen::MatrixXd C0 = Eigen::MatrixXd::Zero(reduced_dim, reduced_dim);
  const Eigen::VectorXd mean_lagged =
      sum_x_lagged_ / static_cast<double>(n_pairs_);
  const Eigen::VectorXd mean_paired =
      sum_x_paired_ / static_cast<double>(n_pairs_);
  const RowVectorXd delta_current =
      (mean_ - mean_paired).transpose() * basis;
  const RowVectorXd delta_lagged =
      (mean_ - mean_lagged).transpose() * basis;

  Eigen::MatrixXd Ctau = Eigen::MatrixXd::Zero(reduced_dim, reduced_dim);
  const int ring_size = std::max(1, lag_);
  FrameMatrixRM lag_buffer = FrameMatrixRM::Zero(ring_size, reduced_dim);
  RowVectorXd z_row(reduced_dim);
  RowVectorXd z_current(reduced_dim);
  RowVectorXd z_lagged(reduced_dim);
  const auto boundaries = build_segment_boundaries_();
  size_t next_boundary = 1;
  int ring_head = 0;
  int ring_count = 0;
  RowVectorXd mean_row = mean_.transpose();
  RowVectorXd centered(dim_);

  const auto process_projected_row = [&](int64_t frame_idx) {
    while (next_boundary < boundaries.size() &&
           boundaries[next_boundary] == frame_idx) {
      ring_head = 0;
      ring_count = 0;
      ++next_boundary;
    }

    C0.noalias() += z_row.transpose() * z_row;

    if (ring_count == lag_) {
      z_current = z_row;
      z_current += delta_current;
      z_lagged = lag_buffer.row(ring_head);
      z_lagged += delta_lagged;
      Ctau.noalias() += z_current.transpose() * z_lagged;

      lag_buffer.row(ring_head) = z_row;
      ring_head = (ring_head + 1) % ring_size;
    } else {
      const int write_slot = (ring_head + ring_count) % ring_size;
      lag_buffer.row(write_slot) = z_row;
      ++ring_count;
    }
  };

  if (backend_ == TicaBackend::ExactDual) {
    for (int64_t i = 0; i < n_frames_; ++i) {
      centered = frame_cache_.row(i);
      centered -= mean_row;
      z_row.noalias() = centered * basis;
      process_projected_row(i);
    }
  } else {
    std::ifstream in(spool_path_, std::ios::binary);
    if (!in.is_open()) {
      throw std::runtime_error("TICA randomized: cannot open spool file for reduced solve.");
    }

    std::vector<double> row_buffer(static_cast<size_t>(dim_));
    for (int64_t i = 0; i < n_frames_; ++i) {
      in.read(reinterpret_cast<char *>(row_buffer.data()),
              static_cast<std::streamsize>(dim_ * sizeof(double)));
      if (!in) {
        throw std::runtime_error("TICA randomized: short read in reduced solve.");
      }
      Eigen::Map<const RowVectorXd> row(row_buffer.data(), dim_);
      centered = row;
      centered -= mean_row;
      z_row.noalias() = centered * basis;
      process_projected_row(i);
    }
  }

  C0 /= static_cast<double>(n_frames_);
  const double epsilon = 1e-6 * C0.trace() / std::max(1, dim_);
  C0.diagonal().array() += epsilon;
  Ctau /= static_cast<double>(n_pairs_);
  Ctau = 0.5 * (Ctau + Ctau.transpose());

  TitanUI::print_info("TICA (" + label + "): reduced C0 trace = " +
                      std::to_string(C0.trace()) +
                      ", regularization eps = " + std::to_string(epsilon));

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(Ctau, C0);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("TICA " + label +
                             ": reduced generalized eigensolver failed.");
  }

  const Eigen::VectorXd all_evals = solver.eigenvalues();
  const Eigen::MatrixXd all_evecs = solver.eigenvectors();
  std::vector<int> idx(all_evals.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&all_evals](int a, int b) { return all_evals[a] > all_evals[b]; });

  eigenvalues_.setZero();
  timescales_.setConstant(std::numeric_limits<double>::infinity());
  eigenvectors_.setZero();
  retained_components_ = std::min<int>(n_comp_, static_cast<int>(all_evals.size()));

  for (int i = 0; i < retained_components_; ++i) {
    eigenvalues_[i] = all_evals[idx[i]];
    eigenvectors_.col(i).noalias() = basis * all_evecs.col(idx[i]);
    timescales_[i] =
        compute_tica_timescale(eigenvalues_[i], lag_, dt_ns_);
  }

  solved_ = true;
}

void TitanKinetics::solve_exact_dual_() {
  Eigen::MatrixXd basis = build_exact_dual_basis_();
  solve_reduced_from_basis_(basis, "exact_dual");
}

void TitanKinetics::solve_randomized_() {
  Eigen::MatrixXd basis = build_randomized_basis_();
  solve_reduced_from_basis_(basis, "randomized");
}

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

  mean_ = sum_x_ / static_cast<double>(n_frames_);

  TitanUI::print_info("TICA: " + std::to_string(n_frames_) + " frames, " +
                      std::to_string(n_pairs_) + " lagged pairs, D=" +
                      std::to_string(dim_) + ", tau=" + std::to_string(lag_) +
                      ", backend=" + tica_backend_to_string(backend_));

  try {
    switch (backend_) {
    case TicaBackend::ExactDense:
      solve_exact_dense_();
      break;
    case TicaBackend::ExactDual:
      solve_exact_dual_();
      break;
    case TicaBackend::Randomized:
      solve_randomized_();
      break;
    case TicaBackend::Auto:
      throw std::runtime_error(
          "TICA: backend 'auto' must be resolved before TitanKinetics construction.");
    }
  } catch (const std::exception &e) {
    TitanUI::print_error(e.what());
    solved_ = false;
    retained_components_ = 0;
    return;
  }

  if (retained_components_ < n_comp_) {
    TitanUI::print_warning("TICA: only " + std::to_string(retained_components_) +
                           " components retained out of requested " +
                           std::to_string(n_comp_) + ".");
  }

  TitanUI::print_info("TICA: solved successfully. Top " +
                      std::to_string(retained_components_) +
                      " components retained.");
}

Eigen::VectorXd TitanKinetics::project_frame(const Eigen::VectorXd &x) const {
  if (!solved_) {
    return Eigen::VectorXd::Zero(n_comp_);
  }
  return (x - mean_).transpose() * eigenvectors_;
}

void TitanKinetics::write_results(const std::string &eval_file,
                                  const std::string &evec_file,
                                  std::fstream &log) const {
  if (!solved_) {
    TitanUI::print_warning("TICA: solve_tica() was not successful, skipping "
                           "output.");
    return;
  }

  {
    std::ofstream ofs(eval_file);
    if (!ofs.is_open()) {
      TitanUI::print_error("TICA: cannot open " + eval_file + " for writing.");
      return;
    }
    ofs << "# TICA Eigenvalues and Implied Timescales\n";
    ofs << "# backend = " << tica_backend_to_string(backend_) << "\n";
    ofs << "# lag_time = " << lag_ << " frames, dt = " << dt_ns_ << " ns\n";
    ofs << "# n_frames = " << n_frames_ << ", n_pairs = " << n_pairs_ << "\n";
    ofs << "# Component  Eigenvalue  Implied_Timescale(ns)\n";
    for (int i = 0; i < retained_components_; ++i) {
      ofs << std::setw(4) << (i + 1) << "  " << std::fixed
          << std::setprecision(8) << eigenvalues_[i] << "  "
          << std::setprecision(4) << timescales_[i] << "\n";
    }
    ofs.close();
    TitanUI::print_info("TICA: eigenvalues written to " + eval_file);
    TitanUI::log_kv(log, "TICA Eigenvalues File", eval_file);
  }

  {
    std::ofstream ofs(evec_file);
    if (!ofs.is_open()) {
      TitanUI::print_error("TICA: cannot open " + evec_file + " for writing.");
      return;
    }
    ofs << "# TICA Eigenvectors (D=" << dim_ << " rows x "
        << retained_components() << " retained components)\n";
    ofs << "# backend = " << tica_backend_to_string(backend_) << "\n";
    ofs << "# Each column is one back-projected eigenvector\n";
    for (int r = 0; r < eigenvectors_.rows(); ++r) {
      for (int c = 0; c < retained_components_; ++c) {
        if (c > 0) {
          ofs << " ";
        }
        ofs << std::fixed << std::setprecision(8) << eigenvectors_(r, c);
      }
      ofs << "\n";
    }
    ofs.close();
    TitanUI::print_info("TICA: eigenvectors written to " + evec_file);
    TitanUI::log_kv(log, "TICA Eigenvectors File", evec_file);
  }
}

void TitanKinetics::write_projections(const MatrixXdAnyConstRef &proj_mat,
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

  ofs << "# frame";
  const int n_cols = static_cast<int>(proj_mat.cols());
  for (int c = 0; c < n_cols; ++c) {
    ofs << "  PC" << (c + 1);
  }
  ofs << "\n";

  ofs << std::fixed << std::setprecision(6);
  for (int64_t i = 0; i < n_actual_frames; ++i) {
    ofs << i;
    for (int c = 0; c < n_cols; ++c) {
      ofs << "  " << proj_mat(i, c);
    }
    ofs << "\n";
  }
  ofs.close();

  TitanUI::print_info("TICA: projections written to " + proj_file + " (" +
                      std::to_string(n_actual_frames) + " frames, " +
                      std::to_string(n_cols) + " components)");
  TitanUI::log_kv(log, "TICA Projections File", proj_file);
}

TitanKinetics::FrameMatrixRM TitanKinetics::project_all_frames() const {
  if (!solved_ || n_frames_ <= 0) {
    return {};
  }
  if (backend_ == TicaBackend::ExactDense) {
    throw std::runtime_error(
        "TICA exact_dense does not retain per-frame history for projection replay.");
  }

  const int basis_cols = std::max(0, retained_components_);
  FrameMatrixRM projected =
      FrameMatrixRM::Zero(n_frames_, std::max(1, n_comp_));
  if (basis_cols == 0) {
    return projected;
  }

  const Eigen::MatrixXd dense_projected =
      project_all_frames_to_basis_(eigenvectors_.leftCols(basis_cols));
  if (dense_projected.rows() != projected.rows()) {
    throw std::runtime_error(
        "TICA: project_all_frames() row-count mismatch.");
  }
  projected.leftCols(basis_cols) = dense_projected;
  return projected;
}
