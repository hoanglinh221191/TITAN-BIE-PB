#include "titan_pca.h"
#include "titan_alignment.h"
#include "titan_ui.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>

namespace TitanPCA {

void execute_pca_fes_analysis(
    const arguments &args,
    const std::vector<std::vector<std::array<double, 3>>> &ca_traj,
    const std::vector<std::string> &ca_resname,
    const std::vector<int> &ca_resnrs, const std::vector<std::string> &ca_chain,
    const std::vector<int> &ca_is_lig, info &info_out, std::fstream &log) {

  if (ca_traj.empty() || ca_traj[0].empty()) {
    TitanUI::print_error("PCA ERROR: No valid trace coordinates found!");
    return;
  }

  size_t num_frames = ca_traj.size();
  size_t num_ca = ca_traj[0].size();
  size_t dimensions = num_ca * 3;

  auto t_start = std::chrono::steady_clock::now();
  TitanUI::log_kv(log, "PCA Snapshot Trace", std::to_string(num_frames));
  TitanUI::log_kv(log, "PCA Degrees of Freedom", std::to_string(dimensions));

  // Create mutable aligned trajectory copy
  std::vector<std::vector<std::array<double, 3>>> aligned_traj = ca_traj;

  // 1. Setup Reference Frame (Frame 0)
  Eigen::MatrixXd ref_mat(num_ca, 3);
  Eigen::RowVector3d ref_centroid(0, 0, 0);
  for (size_t i = 0; i < num_ca; ++i) {
    ref_centroid += Eigen::RowVector3d(
        aligned_traj[0][i][0], aligned_traj[0][i][1], aligned_traj[0][i][2]);
  }
  ref_centroid /= num_ca;

  for (size_t i = 0; i < num_ca; ++i) {
    ref_mat(i, 0) = aligned_traj[0][i][0] - ref_centroid(0);
    ref_mat(i, 1) = aligned_traj[0][i][1] - ref_centroid(1);
    ref_mat(i, 2) = aligned_traj[0][i][2] - ref_centroid(2);
    // Overwrite traj[0] with centered coords
    aligned_traj[0][i][0] = ref_mat(i, 0);
    aligned_traj[0][i][1] = ref_mat(i, 1);
    aligned_traj[0][i][2] = ref_mat(i, 2);
  }

  // 2. Kabsch Alignment for all frames using shared TitanAlignment module
  //    align_frame() centers mobile and rotates in-place; ref_mat stays
  //    unchanged.
  for (size_t f = 1; f < num_frames; ++f) {
    Eigen::MatrixXd frame_mat(num_ca, 3);
    for (size_t i = 0; i < num_ca; ++i) {
      frame_mat(i, 0) = aligned_traj[f][i][0];
      frame_mat(i, 1) = aligned_traj[f][i][1];
      frame_mat(i, 2) = aligned_traj[f][i][2];
    }
    TitanAlignment::align_frame(ref_mat,
                                frame_mat); // centers + rotates in-place
    for (size_t i = 0; i < num_ca; ++i) {
      aligned_traj[f][i][0] = frame_mat(i, 0);
      aligned_traj[f][i][1] = frame_mat(i, 1);
      aligned_traj[f][i][2] = frame_mat(i, 2);
    }
  }

  auto t_rmsd = std::chrono::steady_clock::now();
  TitanUI::log_time(
      log, "RMSD Global Alignment",
      std::chrono::duration_cast<std::chrono::milliseconds>(t_rmsd - t_start)
              .count() /
          1000.0);

  // 3. Compute Mean from aligned_traj
  std::vector<double> mean_coords(dimensions, 0.0);
  for (const auto &coords : aligned_traj) {
    for (size_t i = 0; i < num_ca; ++i) {
      mean_coords[i * 3 + 0] += coords[i][0];
      mean_coords[i * 3 + 1] += coords[i][1];
      mean_coords[i * 3 + 2] += coords[i][2];
    }
  }
  for (size_t i = 0; i < dimensions; ++i) {
    mean_coords[i] /= static_cast<double>(num_frames);
  }

  // 4. Compute Covariance Matrix (Flattened 3N x 3N natively)
  Eigen::MatrixXd cov_matrix = Eigen::MatrixXd::Zero(dimensions, dimensions);

  for (const auto &coords : aligned_traj) {
    std::vector<double> dev(dimensions);
    for (size_t i = 0; i < num_ca; ++i) {
      dev[i * 3 + 0] = coords[i][0] - mean_coords[i * 3 + 0];
      dev[i * 3 + 1] = coords[i][1] - mean_coords[i * 3 + 1];
      dev[i * 3 + 2] = coords[i][2] - mean_coords[i * 3 + 2];
    }

    for (size_t i = 0; i < dimensions; ++i) {
      for (size_t j = i; j < dimensions; ++j) {
        double val = dev[i] * dev[j] / static_cast<double>(num_frames - 1);
        cov_matrix(i, j) += val;
      }
    }
  }

  // Mirror the upper triangle to the lower triangle
  cov_matrix.triangularView<Eigen::StrictlyLower>() = cov_matrix.transpose();

  auto t_cov = std::chrono::steady_clock::now();
  TitanUI::log_time(
      log, "Covariance Building",
      std::chrono::duration_cast<std::chrono::milliseconds>(t_cov - t_rmsd)
              .count() /
          1000.0);

  // 5. Diagonalization natively executing
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(cov_matrix);
  if (es.info() != Eigen::Success) {
    TitanUI::print_error("PCA ERROR: Eigenvalue Decomposition failed!");
    return;
  }

  Eigen::VectorXd eigenvalues = es.eigenvalues();
  Eigen::MatrixXd eigenvectors = es.eigenvectors();

  // Eigen solver returns them in ascending order, we need descending (largest
  // variance first)
  std::vector<double> ev_sorted;
  std::vector<Eigen::VectorXd> evec_sorted;

  for (int i = dimensions - 1; i >= 0; --i) {
    ev_sorted.emplace_back(eigenvalues(i));
    evec_sorted.emplace_back(eigenvectors.col(i));
  }

  auto t_eigen = std::chrono::steady_clock::now();
  TitanUI::log_time(
      log, "Eigendecomposition",
      std::chrono::duration_cast<std::chrono::milliseconds>(t_eigen - t_cov)
              .count() /
          1000.0);

  // 6. Project Trajectory onto PC1 and PC2
  std::vector<std::pair<double, double>> projections(num_frames);
  const Eigen::VectorXd &pc1_vec = evec_sorted[0];
  const Eigen::VectorXd &pc2_vec = evec_sorted[1];

  double min_pc1 = 1e9, max_pc1 = -1e9;
  double min_pc2 = 1e9, max_pc2 = -1e9;

  for (size_t f = 0; f < num_frames; ++f) {
    double p1 = 0.0, p2 = 0.0;
    for (size_t i = 0; i < num_ca; ++i) {
      for (int j = 0; j < 3; ++j) {
        int dim_idx = i * 3 + j;
        double diff = aligned_traj[f][i][j] - mean_coords[dim_idx];
        p1 += diff * pc1_vec(dim_idx);
        p2 += diff * pc2_vec(dim_idx);
      }
    }
    projections[f] = {p1, p2};

    if (p1 < min_pc1)
      min_pc1 = p1;
    if (p1 > max_pc1)
      max_pc1 = p1;
    if (p2 < min_pc2)
      min_pc2 = p2;
    if (p2 > max_pc2)
      max_pc2 = p2;
  }

  // Export PCA projections smoothly natively avoiding OpenMP overlapping
  // instances
  std::ofstream pca_out(info_out.get_pca_fel_filename());
  pca_out << "Frame,PC1,PC2\n";
  for (size_t f = 0; f < num_frames; ++f) {
    pca_out << f << "," << std::fixed << std::setprecision(4)
            << projections[f].first << "," << projections[f].second << "\n";
  }
  pca_out.close();

  // 7. 2D FES Calculation explicitly extracting Histogram Density Peak mapping
  // cleanly
  int grid_size = args.fes_grid;
  if (grid_size <= 0)
    grid_size = 60; // Failsafe

  double pc1_range = max_pc1 - min_pc1;
  double pc2_range = max_pc2 - min_pc2;

  if (pc1_range < 1e-6)
    pc1_range = 1e-6;
  if (pc2_range < 1e-6)
    pc2_range = 1e-6;

  double step_pc1 = pc1_range / grid_size;
  double step_pc2 = pc2_range / grid_size;

  std::vector<std::vector<FES_Grid>> fes_map(grid_size,
                                             std::vector<FES_Grid>(grid_size));

  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      fes_map[i][j] = {min_pc1 + step_pc1 * (i + 0.5),
                       min_pc2 + step_pc2 * (j + 0.5), 0, 20.0, -1};
    }
  }

  // Populate Histogram
  int max_count = 0;
  for (size_t f = 0; f < num_frames; ++f) {
    int i = std::floor((projections[f].first - min_pc1) / step_pc1);
    int j = std::floor((projections[f].second - min_pc2) / step_pc2);

    if (i >= grid_size)
      i = grid_size - 1;
    if (i < 0)
      i = 0;
    if (j >= grid_size)
      j = grid_size - 1;
    if (j < 0)
      j = 0;

    fes_map[i][j].count++;
    if (fes_map[i][j].count > max_count) {
      max_count = fes_map[i][j].count;
    }
  }

  // Boltzmann Inversion explicitly scaling energies preventing zeros natively
  // kB*T at ~300K: 0.596 kcal/mol = 2.494 kJ/mol
  // eu_par = 4.184 when eu=kcal/mol (divisor to convert kJ->kcal)
  // eu_par = 1.0   when eu=kJ/mol  (no conversion, keep kJ)
  // Use kbT in kJ/mol (2.494) then divide by eu_par for final unit.
  double kbT = 2.494; // kJ/mol at ~300K
  double scale_factor = 1.0 / args.eu_par;

  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      if (fes_map[i][j].count > 0) {
        double prob = static_cast<double>(fes_map[i][j].count) / max_count;
        fes_map[i][j].free_energy = -kbT * scale_factor * std::log(prob);
      }
    }
  }

  // Export 2D FES
  std::ofstream fes_out(info_out.get_fes_2d_filename());
  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      fes_out << std::fixed << std::setprecision(4) << fes_map[i][j].pc1 << " "
              << fes_map[i][j].pc2 << " " << fes_map[i][j].free_energy << "\n";
    }
    fes_out << "\n"; // Gnuplot block separation natively
  }
  fes_out.close();

  // 6. Density-Aware Basin Clustering
  std::vector<ClusterNode> clusters;
  int current_cluster_id = 0;

  // Find Sinks (local minima) explicitly searching neighbors cleanly
  for (int i = 1; i < grid_size - 1; ++i) {
    for (int j = 1; j < grid_size - 1; ++j) {
      if (fes_map[i][j].count == 0)
        continue;

      double e = fes_map[i][j].free_energy;
      bool is_minimum = true;

      for (int di = -1; di <= 1 && is_minimum; ++di) {
        for (int dj = -1; dj <= 1 && is_minimum; ++dj) {
          if (di == 0 && dj == 0)
            continue;
          if (fes_map[i + di][j + dj].free_energy < e) {
            is_minimum = false;
          }
        }
      }

      if (is_minimum) {
        fes_map[i][j].sink_id = current_cluster_id++;
        ClusterNode c;
        c.id = fes_map[i][j].sink_id;
        c.min_pc1 = fes_map[i][j].pc1;
        c.min_pc2 = fes_map[i][j].pc2;
        c.min_energy = e;
        c.representative_frame_idx = -1;
        c.min_euclidean_distance = 1e9;
        c.population = 0;
        c.percentage = 0.0;
        clusters.emplace_back(c);
      }
    }
  }

  // Assign frames to nearest sink strictly matching nearest nodes linearly
  // isolating noisy bins
  for (size_t f = 0; f < num_frames; ++f) {
    double p1 = projections[f].first;
    double p2 = projections[f].second;

    int best_cluster = -1;
    double min_dist = 1e9;

    for (size_t c = 0; c < clusters.size(); ++c) {
      double dist = std::sqrt(std::pow(p1 - clusters[c].min_pc1, 2) +
                              std::pow(p2 - clusters[c].min_pc2, 2));
      if (dist < min_dist) {
        min_dist = dist;
        best_cluster = c;
      }
    }

    if (best_cluster != -1) {
      clusters[best_cluster].population++;

      // Track representative tracking shortest bounds solving extraction goals
      if (min_dist < clusters[best_cluster].min_euclidean_distance) {
        clusters[best_cluster].min_euclidean_distance = min_dist;
        clusters[best_cluster].representative_frame_idx = f;
      }
    }
  }

  // Reject noise clusters (< 2.0% frames) explicitly enforcing physics bounds
  // resolving anomalies
  double threshold = 0.02 * num_frames;
  std::vector<ClusterNode> valid_clusters;

  for (auto &c : clusters) {
    c.percentage = 100.0 * static_cast<double>(c.population) / num_frames;
    if (c.population >= threshold) {
      valid_clusters.emplace_back(c);
    }
  }

  // Export Cluster Representatives explicitly logging representations cleanly
  std::ofstream rep_out(info_out.get_cluster_representatives_filename());
  rep_out << "ClusterID,Population%,DeltaG_Sink,Representative_Frame,PC1,PC2\n";

  for (const auto &c : valid_clusters) {
    rep_out << c.id << "," << std::fixed << std::setprecision(2) << c.percentage
            << "%," << std::setprecision(4) << c.min_energy << ","
            << c.representative_frame_idx << "," << c.min_pc1 << ","
            << c.min_pc2 << "\n";
  }
  rep_out.close();

  TitanUI::log_kv(log, "Macrostate Clusters Extracted",
                  std::to_string(valid_clusters.size()));

  auto t_end = std::chrono::steady_clock::now();
  TitanUI::log_time(
      log, "PCA + FES Clustering",
      std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_eigen)
          .count());
}

} // namespace TitanPCA
