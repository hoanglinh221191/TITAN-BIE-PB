// ============================================================================
// TitanMSM — Markov State Model after TICA (Phase 2)
//
// Workflow:
//   1. Call add_projected_frame() for every frame from the TICA projection.
//   2. Call add_trajectory_break() after every trajectory file boundary.
//   3. train_kmeans()        : K-Means++ → discrete microstate trajectory.
//   4. build_markov_model()  : boundary-aware count matrix → TPM + stationary
//   dist.
//   5. compute_implied_timescales() : eigendecomposition for ITS.
//   6. export_results()      : write all output files + PCCA + ITS scan.
//
// Eigen only for linear algebra. No std::cout / printf anywhere.
// HOT LOOPS: zero heap allocation per K-Means iteration.
// ============================================================================
#ifndef TITAN_MSM_H
#define TITAN_MSM_H

#include "info.h"
#include "titan_ui.h"
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

class TitanMSM {
public:
  // -----------------------------------------------------------------------
  // Construction
  // n_clusters        : number of K-Means microstates (k)
  // n_tica_components : dimensionality of the projected space
  // dt_ns             : physical timestep between frames (nanoseconds)
  // -----------------------------------------------------------------------
  TitanMSM(int n_clusters, int n_tica_components, double dt_ns);

  // -----------------------------------------------------------------------
  // Trajectory accumulation
  // -----------------------------------------------------------------------
  // Push one TICA-projected frame into RAM.
  void add_projected_frame(const Eigen::VectorXd &proj_x,
                           int replica_id = -1);

  // Record the simulation time (ps) of the most recently added frame.
  void record_frame_time(double time_ps, int replica_id = -1);

  // Per-replica accessors (populated when replica_id != -1).
  const std::vector<int>&    get_frame_replica_ids() const { return frame_replica_ids_; }
  const std::vector<double>& get_frame_timestamps()  const { return frame_timestamps_; }

  // Per-replica discrete trajectory extraction (call after train_kmeans).
  std::unordered_map<int, std::vector<int>>    get_per_replica_dtrajs() const;
  std::unordered_map<int, std::vector<double>> get_per_replica_timestamps() const;

  // Mark a trajectory file boundary. The next add_projected_frame() call
  // starts a new segment. Must be called AFTER the last frame of each file.
  void add_trajectory_break();

  // Reset all trajectory data (projected frames + breaks) without touching
  // trained centroids. Useful when re-running multiple analyses.
  void clear_projected();

  // -----------------------------------------------------------------------
  // Phase 2A — K-Means clustering
  // max_iter : maximum EM iterations
  // seed     : RNG seed for K-Means++ initialisation
  // -----------------------------------------------------------------------
  void train_kmeans(int max_iter = 300, unsigned int seed = 42,
                    int n_threads = 1);

  // -----------------------------------------------------------------------
  // Phase 2B — Build Markov State Model
  // msm_lag_time : lag time in frames (must be > 0)
  // -----------------------------------------------------------------------
  void build_markov_model(int msm_lag_time);

  // -----------------------------------------------------------------------
  // Phase 2C — Implied timescales (eigendecomposition of TPM)
  // Returns a vector of timescales in nanoseconds (length <= n_timescales).
  // -----------------------------------------------------------------------
  std::vector<double> compute_implied_timescales(int n_timescales) const;

  // -----------------------------------------------------------------------
  // Phase 2D — Export all results
  // Runs ITS scan (if msm_run_its=true) and PCCA (if msm_run_pcca=true).
  // n_timescales      : number of timescales to report
  // msm_run_its       : perform multi-lag ITS scan
  // msm_run_pcca      : perform crisp PCCA macrostate lumping
  // msm_macrostates   : -1 = auto-detect from timescales
  // its_lags          : lag list for ITS scan (empty = auto-generate)
  // -----------------------------------------------------------------------
  void export_results(const info &info_output, std::fstream &log,
                      int n_timescales = 3, bool msm_run_its = true,
                      bool msm_run_pcca = true, int msm_macrostates = -1,
                      const std::vector<int> &its_lags = {}) const;

  // -----------------------------------------------------------------------
  // Accessors
  // -----------------------------------------------------------------------
  int n_clusters() const { return n_k_; }
  int n_components() const { return n_comp_; }
  int msm_lag_time() const { return msm_lag_time_; }
  bool is_clustered() const { return clustered_; }
  bool is_built() const { return msm_built_; }
  // Fix 3: n_frames reads n_proj_frames_ (proj_mat_ actual row count)
  int64_t n_frames() const { return n_proj_frames_; }
  int n_threads() const { return n_threads_; }
  // Accessor to the contiguous RowMajor projection matrix for
  // write_projections()
  const Eigen::MatrixXd &proj_mat_ref() const { return proj_mat_; }

  const Eigen::MatrixXd &centroids() const { return centroids_; }
  const Eigen::MatrixXd &tpm() const { return tpm_; }
  const Eigen::VectorXd &stationary_dist() const { return stationary_dist_; }
  const std::vector<int> &discrete_traj() const { return discrete_traj_; }

private:
  // -----------------------------------------------------------------------
  // Configuration
  // -----------------------------------------------------------------------
  int n_k_;          // number of clusters
  int n_comp_;       // TICA projected dimensions
  double dt_ns_;     // physical timestep (ns)
  int msm_lag_time_; // set by build_markov_model()

  // -----------------------------------------------------------------------
  // Trajectory data (Fix 3: RowMajor MatrixXd instead of vector<VectorXd>)
  // proj_mat_ rows are contiguous n_comp_ doubles → cache-friendly + SIMD
  // -----------------------------------------------------------------------
  // Fix 3: RowMajor layout: each row is one frame, rows are contiguous
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      proj_mat_;
  int64_t n_proj_frames_ = 0;           // number of frames actually stored
  std::vector<int> traj_break_indices_; // segment boundary frame indices
  std::vector<int> discrete_traj_;      // [N_frames] cluster assignments
  int n_threads_ = 1;                   // C++ thread count for K-Means assign

  // DHAM replica tracking
  std::vector<int>    frame_replica_ids_;   // replica_id per proj_mat_ row
  std::vector<double> frame_timestamps_;    // time_ps per proj_mat_ row

  // -----------------------------------------------------------------------
  // K-Means results
  // -----------------------------------------------------------------------
  Eigen::MatrixXd centroids_; // [k x n_comp_]
  bool clustered_ = false;

  // -----------------------------------------------------------------------
  // MSM results
  // -----------------------------------------------------------------------
  Eigen::MatrixXd count_matrix_;    // [k x k] double (for symmetrization)
  Eigen::MatrixXd tpm_;             // [k x k] row-stochastic TPM
  Eigen::VectorXd stationary_dist_; // [k] left eigenvector (pi)
  bool msm_built_ = false;

  // -----------------------------------------------------------------------
  // Private helpers
  // -----------------------------------------------------------------------

  // Generate a geometric ITS lag sequence: 1,2,5,10,20,50,...≤ max_usable.
  std::vector<int> generate_its_lags() const;

  // Run crisp PCCA: returns macrostate index per microstate (size = k).
  // Fix 4: uses cached_eigenvectors_ from build_markov_model().
  std::vector<int> run_pcca_crisp(int n_macro) const;

  // Compute n_macro from timescales (count t_i > 2 * lag * dt, +1, clamp).
  int auto_detect_macrostates(const std::vector<double> &timescales) const;

  // Fix 4 — cached eigendecomposition of tpm_ (computed once in
  // build_markov_model)
  mutable Eigen::VectorXcd cached_eigenvalues_;
  mutable Eigen::MatrixXcd cached_eigenvectors_;
  mutable bool eigen_solved_ = false;
};

// =======================================================================
// === ENGINE INTEGRATION CONCEPT (conceptual block for impl_engine.cpp) =
// =======================================================================
//
//  // --- Construct MSM analyzer (after TICA is solved) ---
//  TitanMSM msm_analyzer(Arguments.msm.msm_clusters,
//                         tica_analyzer->n_components(),
//                         dt_ns);
//
//  // --- Second Pass: Project + Collect ---
//  msm_analyzer.clear_projected();
//  int last_file_idx = -1;
//  while (chunk_streamer has chunks) {
//      if (chunk_streamer.current_file_idx() != last_file_idx) {
//          if (last_file_idx != -1 && Arguments.msm.msm_auto_boundary)
//              msm_analyzer.add_trajectory_break();
//          last_file_idx = chunk_streamer.current_file_idx();
//      }
//      for each frame in chunk:
//          auto proj = tica_analyzer->project_frame(feature_vector(frame));
//          msm_analyzer.add_projected_frame(proj);
//  }
//  // Final boundary for the last file
//  if (Arguments.msm.msm_auto_boundary)
//      msm_analyzer.add_trajectory_break();
//
//  // --- Build MSM ---
//  msm_analyzer.train_kmeans(Arguments.msm.msm_kmeans_iter, Arguments.msm.msm_seed);
//  msm_analyzer.build_markov_model(Arguments.msm.msm_lag);
//
//  // --- Export (ITS scan + PCCA inside) ---
//  msm_analyzer.export_results(info_output, log_file,
//                               Arguments.msm.msm_nts,
//                               Arguments.msm.msm_run_its,
//                               Arguments.msm.msm_run_pcca,
//                               Arguments.msm.msm_macrostates,
//                               Arguments.msm.msm_its_lags);

#endif // TITAN_MSM_H
