#ifndef TITAN_KINETICS_H
#define TITAN_KINETICS_H

#include "titan_ui.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// ============================================================================
// TitanKinetics — Time-lagged Independent Component Analysis (TICA)
//
// Uses the "accumulate-then-center" online algorithm.
// All D×D matrices are pre-allocated in the constructor.
// ZERO heap allocation inside accumulate_frame().
// ============================================================================
class TitanKinetics {
public:
  // ---- Construction --------------------------------------------------------
  // feature_dim : dimensionality D of the feature vector (e.g. 3*N_atoms)
  // lag_time    : tau in frames for the time-lagged correlation
  // n_components: number of slow TICA components to retain
  // dt_ns       : physical timestep between frames in nanoseconds
  TitanKinetics(int feature_dim, int lag_time, int n_components, double dt_ns);

  // ---- Hot-loop methods (called once per frame) ----------------------------
  // Push a single feature vector and update all accumulators.
  // MUST NOT allocate heap memory.
  void accumulate_frame(const Eigen::VectorXd &x);

  // Flush the lag buffer at trajectory boundaries to prevent cross-trajectory
  // lag pairs. The global accumulators are NOT reset.
  void reset_lag_buffer();

  // ---- Post-loop methods ---------------------------------------------------
  // Build covariance matrices, solve generalised eigenvalue problem,
  // sort eigenvectors by descending eigenvalue, compute implied timescales.
  void solve_tica();

  // Project a single frame onto the TICA subspace.
  // Returns an n_components-dimensional vector.
  Eigen::VectorXd project_frame(const Eigen::VectorXd &x) const;

  // Write eigenvalues + implied timescales and eigenvectors to disk.
  void write_results(const std::string &eval_file, const std::string &evec_file,
                     std::fstream &log) const;

  // Write TICA projections from a PRE-PROJECTED matrix (rows = frames, cols =
  // PCs). Use with msm_analyzer->proj_mat_ref() which already holds TICA-space
  // coordinates. n_actual_frames: real row count (proj_mat may be
  // over-allocated with doubling). Each row [0..n_actual_frames) is written
  // as-is — no re-projection applied.
  void write_projections(const Eigen::MatrixXd &proj_mat,
                         int64_t n_actual_frames, const std::string &proj_file,
                         std::fstream &log) const;

  // ---- Accessors -----------------------------------------------------------
  int feature_dim() const { return dim_; }
  int lag_time() const { return lag_; }
  int n_components() const { return n_comp_; }
  int64_t num_frames() const { return n_frames_; }
  int64_t num_pairs() const { return n_pairs_; }
  bool is_solved() const { return solved_; }

  const Eigen::VectorXd &eigenvalues() const { return eigenvalues_; }
  const Eigen::MatrixXd &eigenvectors() const { return eigenvectors_; }
  const Eigen::VectorXd &implied_timescales() const { return timescales_; }
  const Eigen::VectorXd &mean() const { return mean_; }

private:
  // Dimensions
  int dim_;
  int lag_;
  int n_comp_;
  double dt_ns_;

  // Online accumulators (all pre-allocated)
  Eigen::VectorXd sum_x_;        // (D×1) sum of features at time t
  Eigen::MatrixXd sum_xx_;       // (D×D) sum of x_t * x_t^T
  Eigen::VectorXd sum_x_lagged_; // (D×1) sum of features at time t-tau
  Eigen::VectorXd sum_x_paired_; // (D×1) sum of x_t for paired frames only
  Eigen::MatrixXd sum_xtau_;     // (D×D) sum of x_t * x_{t-tau}^T
  int64_t n_frames_;
  int64_t n_pairs_;

  // Lag buffer — circular, pre-allocated (ZERO heap alloc per frame)
  Eigen::MatrixXd circ_buf_; // (lag_ × dim_) pre-allocated ring buffer
  int circ_head_ = 0;        // index of oldest entry
  int circ_count_ = 0;       // current fill level (< lag_ during warm-up)

  // Results (populated by solve_tica)
  bool solved_;
  Eigen::VectorXd mean_;
  Eigen::VectorXd eigenvalues_;
  Eigen::MatrixXd eigenvectors_; // (D × n_components)
  Eigen::VectorXd timescales_;
};

#endif // TITAN_KINETICS_H
