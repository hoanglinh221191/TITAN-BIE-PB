#include "titan_pca.h"
#include "titan_alignment.h"
#include "titan_executor.h"
#include "titan_ui.h"
#include "titan_zts.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <span>
#include <set>
#include <sstream>
#include <thread>

#ifndef _WIN32
#include <sys/stat.h>
#endif

namespace {

constexpr double kBoltzmannKjPerMolK = 0.008314463;
constexpr double kMinimumAxisRange = 1e-6;
constexpr double kDefaultFreeEnergy = 20.0;
using TitanPCA::FES2DResult;

struct FrameWeightResolution {
  bool use_weights{false};
  AlignedDoubleVector analysis_weights{};
  double analysis_weight_sum{0.0};
};

struct FlatCoordinateMatrix {
  static constexpr size_t kSimdLaneWidth =
      std::max<size_t>(1, kTitanSimdAlignment / sizeof(double));

  FlatCoordinateMatrix() = default;
  FlatCoordinateMatrix(size_t rows, size_t cols, double value = 0.0) {
    assign(rows, cols, value);
  }

  void assign(size_t rows, size_t cols, double value = 0.0) {
    rows_ = rows;
    cols_ = cols;
    stride_ = cols == 0 ? 0
                        : ((cols + kSimdLaneWidth - 1) / kSimdLaneWidth) *
                              kSimdLaneWidth;
    data_.assign(rows * stride_, value);
  }

  double &operator()(size_t row, size_t col) {
    return data_[row * stride_ + col];
  }

  const double &operator()(size_t row, size_t col) const {
    return data_[row * stride_ + col];
  }

  double *row_data(size_t row) { return data_.data() + row * stride_; }
  const double *row_data(size_t row) const {
    return data_.data() + row * stride_;
  }

private:
  size_t rows_{0};
  size_t cols_{0};
  size_t stride_{0};
  AlignedDoubleVector data_{};
};

size_t choose_worker_count(size_t requested_threads, size_t work_items,
                           size_t min_items_per_worker) {
  if (work_items == 0)
    return 1;
  const size_t capped_threads = std::max<size_t>(1, requested_threads);
  const size_t max_useful =
      std::max<size_t>(1, work_items / std::max<size_t>(1, min_items_per_worker));
  return std::max<size_t>(1, std::min(capped_threads, max_useful));
}

FrameWeightResolution resolve_frame_weights(const std::string &label,
                                            size_t num_frames,
                                            std::span<const double> weights) {
  FrameWeightResolution resolved;
  resolved.analysis_weights.assign(num_frames, 1.0);
  resolved.analysis_weight_sum = static_cast<double>(num_frames);

  if (weights.empty())
    return resolved;

  if (weights.size() != num_frames) {
    TitanUI::print_warning(
        label + ": frame_weights.size() (" + std::to_string(weights.size()) +
        ") != num_frames (" + std::to_string(num_frames) +
        "). Falling back to uniform weights. Ensure DHAM and this analysis "
        "use the same trajectory and frame stride.");
    return resolved;
  }

  resolved.analysis_weights.assign(weights.begin(), weights.end());
  resolved.analysis_weight_sum =
      std::accumulate(resolved.analysis_weights.begin(),
                      resolved.analysis_weights.end(), 0.0);
  if (resolved.analysis_weight_sum <= 0.0) {
    TitanUI::print_warning(label +
                           ": frame weights sum to zero. Falling back to "
                           "uniform weights.");
    resolved.analysis_weights.assign(num_frames, 1.0);
    resolved.analysis_weight_sum = static_cast<double>(num_frames);
    return resolved;
  }

  resolved.use_weights = true;
  return resolved;
}

AlignedDoubleVector map_global_weights_to_local_frames(
    const std::string &label, int64_t n_frames,
    const std::vector<int64_t> &frame_global_indices,
    const std::vector<double> &global_frame_weights) {
  if (global_frame_weights.empty())
    return {};
  if (n_frames <= 0)
    return {};

  const size_t local_frames = static_cast<size_t>(n_frames);
  if (frame_global_indices.size() != local_frames) {
    TitanUI::print_warning(
        label + ": frame_global_indices.size() (" +
        std::to_string(frame_global_indices.size()) + ") != n_frames (" +
        std::to_string(local_frames) +
        "). Falling back to uniform weights.");
    return {};
  }

  AlignedDoubleVector mapped(local_frames, 0.0);
  size_t matched_frames = 0;
  size_t weighted_frames = 0;
  size_t invalid_indices = 0;
  double weight_sum = 0.0;
  for (size_t f = 0; f < local_frames; ++f) {
    const int64_t global_idx = frame_global_indices[f];
    if (global_idx < 0 ||
        static_cast<size_t>(global_idx) >= global_frame_weights.size()) {
      ++invalid_indices;
      continue;
    }
    const double w = global_frame_weights[static_cast<size_t>(global_idx)];
    mapped[f] = w;
    ++matched_frames;
    if (w > 0.0)
      ++weighted_frames;
    weight_sum += w;
  }

  TitanUI::print_info(label + ": Mapped DHAM weights onto " +
                      std::to_string(weighted_frames) + "/" +
                      std::to_string(local_frames) +
                      " projected frames (matched indices: " +
                      std::to_string(matched_frames) + ", weight sum = " +
                      std::to_string(weight_sum) + ").");
  if (invalid_indices > 0) {
    TitanUI::print_warning(label + ": Dropped " +
                           std::to_string(invalid_indices) +
                           " projected frames with invalid global indices "
                           "during DHAM weight mapping.");
  }

  return mapped;
}

double resolve_target_temperature(const arguments &args, bool use_weights) {
  if (use_weights && args.dham.flag_dham && args.dham.target_temp > 0.0)
    return args.dham.target_temp;
  return 300.0;
}

FES2DResult compute_2d_fes_grid(
    const std::vector<std::pair<double, double>> &projections,
    const AlignedDoubleVector &analysis_weights, int requested_grid_size,
    double target_temp, double scale_factor, TitanExecutor &executor) {
  FES2DResult result;
  if (projections.empty())
    return result;

  result.grid_size = (requested_grid_size > 0) ? requested_grid_size : 60;
  result.total_weight =
      std::accumulate(analysis_weights.begin(), analysis_weights.end(), 0.0);

  double min_x = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::lowest();
  for (const auto &[x, y] : projections) {
    min_x = std::min(min_x, x);
    max_x = std::max(max_x, x);
    min_y = std::min(min_y, y);
    max_y = std::max(max_y, y);
  }

  const double range_x = std::max(max_x - min_x, kMinimumAxisRange);
  const double range_y = std::max(max_y - min_y, kMinimumAxisRange);
  result.min_x = min_x;
  result.min_y = min_y;
  result.step_x = range_x / static_cast<double>(result.grid_size);
  result.step_y = range_y / static_cast<double>(result.grid_size);

  const size_t grid_cells = static_cast<size_t>(result.grid_size) *
                            static_cast<size_t>(result.grid_size);
  result.bins.resize(grid_cells);
  for (int ix = 0; ix < result.grid_size; ++ix) {
    for (int iy = 0; iy < result.grid_size; ++iy) {
      auto &cell = result.bin(ix, iy);
      cell = {result.min_x + result.step_x * (static_cast<double>(ix) + 0.5),
              result.min_y + result.step_y * (static_cast<double>(iy) + 0.5), 0,
              0.0, kDefaultFreeEnergy, -1};
    }
  }

  const size_t worker_count =
      choose_worker_count(executor.capacity(), projections.size(), 2048);
  std::vector<AlignedIntVector> local_counts(worker_count,
                                             AlignedIntVector(grid_cells, 0));
  std::vector<AlignedDoubleVector> local_weights(
      worker_count, AlignedDoubleVector(grid_cells, 0.0));
  const double inv_step_x = 1.0 / result.step_x;
  const double inv_step_y = 1.0 / result.step_y;

  executor.parallel_do(worker_count, [&](size_t worker) {
    const size_t begin = projections.size() * worker / worker_count;
    const size_t end = projections.size() * (worker + 1) / worker_count;
    auto &counts = local_counts[worker];
    auto &weights = local_weights[worker];
    for (size_t frame = begin; frame < end; ++frame) {
      const auto &[x, y] = projections[frame];
      int ix = static_cast<int>(std::floor((x - result.min_x) * inv_step_x));
      int iy = static_cast<int>(std::floor((y - result.min_y) * inv_step_y));
      ix = std::clamp(ix, 0, result.grid_size - 1);
      iy = std::clamp(iy, 0, result.grid_size - 1);
      const size_t bin_index =
          static_cast<size_t>(ix) * static_cast<size_t>(result.grid_size) +
          static_cast<size_t>(iy);
      ++counts[bin_index];
      weights[bin_index] += analysis_weights[frame];
    }
  });

  for (size_t cell = 0; cell < grid_cells; ++cell) {
    int count = 0;
    double weight = 0.0;
    for (size_t worker = 0; worker < worker_count; ++worker) {
      count += local_counts[worker][cell];
      weight += local_weights[worker][cell];
    }
    result.bins[cell].count = count;
    result.bins[cell].count_w = weight;
    result.max_bin_weight = std::max(result.max_bin_weight, weight);
  }

  if (result.total_weight <= 0.0 || result.max_bin_weight <= 0.0)
    return result;

  const double kbT = kBoltzmannKjPerMolK * target_temp;
  const size_t fe_workers =
      choose_worker_count(executor.capacity(), grid_cells, 1024);
  executor.parallel_do(fe_workers, [&](size_t worker) {
    const size_t begin = grid_cells * worker / fe_workers;
    const size_t end = grid_cells * (worker + 1) / fe_workers;
    for (size_t cell = begin; cell < end; ++cell) {
      auto &bin = result.bins[cell];
      if (bin.count_w > 0.0) {
        const double prob = bin.count_w / result.max_bin_weight;
        bin.free_energy = -kbT * scale_factor * std::log(prob);
      }
    }
  });

  return result;
}

bool write_2d_fes_grid(const std::string &label, const std::string &path,
                       const FES2DResult &fes_result, std::fstream &log) {
  if (fes_result.total_weight <= 0.0 || fes_result.max_bin_weight <= 0.0) {
    TitanUI::print_warning(label +
                           ": histogram has zero total weight. Skipping FES "
                           "construction.");
    return false;
  }

  std::ofstream fes_out(path);
  if (!fes_out.is_open()) {
    TitanUI::print_error(label + ": cannot open " + path +
                         " for writing 2D FES.");
    return false;
  }

  for (int ix = 0; ix < fes_result.grid_size; ++ix) {
    for (int iy = 0; iy < fes_result.grid_size; ++iy) {
      const auto &bin = fes_result.bin(ix, iy);
      fes_out << std::fixed << std::setprecision(4) << bin.pc1 << " " << bin.pc2
              << " " << bin.free_energy << "\n";
    }
    fes_out << "\n";
  }

  TitanUI::print_info(label + ": 2D FES written to " + path);
  TitanUI::log_kv(log, label + " 2D FES File", path);
  return true;
}

std::vector<TitanPCA::ClusterNode> detect_fes_clusters(
    FES2DResult &fes_result,
    const std::vector<std::pair<double, double>> &projections,
    const AlignedDoubleVector &analysis_weights, bool use_weights,
    TitanExecutor &executor) {
  if (fes_result.grid_size <= 2 || projections.empty())
    return {};

  const size_t grid_cells = static_cast<size_t>(fes_result.grid_size) *
                            static_cast<size_t>(fes_result.grid_size);
  AlignedIntVector minima_flags(grid_cells, 0);
  const size_t worker_count =
      choose_worker_count(executor.capacity(), grid_cells, 256);
  executor.parallel_do(worker_count, [&](size_t worker) {
    const size_t begin = grid_cells * worker / worker_count;
    const size_t end = grid_cells * (worker + 1) / worker_count;
    for (size_t linear = begin; linear < end; ++linear) {
      const int ix = static_cast<int>(linear /
                                      static_cast<size_t>(fes_result.grid_size));
      const int iy = static_cast<int>(linear %
                                      static_cast<size_t>(fes_result.grid_size));
      if (ix <= 0 || iy <= 0 || ix >= fes_result.grid_size - 1 ||
          iy >= fes_result.grid_size - 1) {
        continue;
      }
      if (fes_result.bin(ix, iy).count_w <= 0.0)
        continue;

      const double e = fes_result.bin(ix, iy).free_energy;
      bool is_minimum = true;
      for (int dx = -1; dx <= 1 && is_minimum; ++dx) {
        for (int dy = -1; dy <= 1 && is_minimum; ++dy) {
          if (dx == 0 && dy == 0)
            continue;
          if (fes_result.bin(ix + dx, iy + dy).free_energy < e)
            is_minimum = false;
        }
      }
      minima_flags[linear] = is_minimum ? 1 : 0;
    }
  });

  std::vector<TitanPCA::ClusterNode> clusters;
  int current_cluster_id = 0;
  for (size_t linear = 0; linear < grid_cells; ++linear) {
    if (minima_flags[linear] == 0)
      continue;
    const int ix = static_cast<int>(linear /
                                    static_cast<size_t>(fes_result.grid_size));
    const int iy = static_cast<int>(linear %
                                    static_cast<size_t>(fes_result.grid_size));
    auto &bin = fes_result.bin(ix, iy);
    bin.sink_id = current_cluster_id++;
    TitanPCA::ClusterNode cluster{};
    cluster.id = bin.sink_id;
    cluster.min_pc1 = bin.pc1;
    cluster.min_pc2 = bin.pc2;
    cluster.min_energy = bin.free_energy;
    cluster.representative_frame_idx = -1;
    cluster.min_euclidean_distance = 1e9;
    cluster.population = 0;
    cluster.population_w = 0.0;
    cluster.percentage = 0.0;
    cluster.percentage_w = 0.0;
    clusters.emplace_back(cluster);
  }

  for (size_t frame = 0; frame < projections.size(); ++frame) {
    const double p1 = projections[frame].first;
    const double p2 = projections[frame].second;
    int best_cluster = -1;
    double min_dist = 1e9;
    for (size_t c = 0; c < clusters.size(); ++c) {
      const double dx = p1 - clusters[c].min_pc1;
      const double dy = p2 - clusters[c].min_pc2;
      const double dist = std::sqrt(dx * dx + dy * dy);
      if (dist < min_dist) {
        min_dist = dist;
        best_cluster = static_cast<int>(c);
      }
    }
    if (best_cluster < 0)
      continue;
    auto &cluster = clusters[static_cast<size_t>(best_cluster)];
    cluster.population++;
    cluster.population_w += analysis_weights[frame];
    if (min_dist < cluster.min_euclidean_distance) {
      cluster.min_euclidean_distance = min_dist;
      cluster.representative_frame_idx = static_cast<int>(frame);
    }
  }

  const double threshold =
      0.02 * (use_weights ? fes_result.total_weight
                          : static_cast<double>(projections.size()));
  std::vector<TitanPCA::ClusterNode> valid_clusters;
  for (auto &cluster : clusters) {
    cluster.percentage =
        100.0 * static_cast<double>(cluster.population) / projections.size();
    cluster.percentage_w =
        (fes_result.total_weight > 0.0)
            ? (100.0 * cluster.population_w / fes_result.total_weight)
            : 0.0;
    if ((use_weights && cluster.population_w >= threshold) ||
        (!use_weights &&
         static_cast<double>(cluster.population) >= threshold)) {
      valid_clusters.emplace_back(cluster);
    }
  }
  return valid_clusters;
}

std::string bash_quote(const std::string &value) {
  std::string quoted = "'";
  for (char ch : value) {
    if (ch == '\'')
      quoted += "'\"'\"'";
    else
      quoted.push_back(ch);
  }
  quoted.push_back('\'');
  return quoted;
}

std::string format_time_ps(double time_ps, int precision) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(precision) << time_ps;
  return oss.str();
}

std::string make_time_label(double time_ps) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(1) << time_ps;
  std::string label = oss.str();
  std::replace(label.begin(), label.end(), '.', 'p');
  std::replace(label.begin(), label.end(), '-', 'm');
  return label;
}

bool write_cluster_representatives_file(
    const std::string &label, const std::string &space_name,
    const std::string &filename,
    const std::vector<TitanPCA::ClusterNode> &clusters, bool use_weights,
    const std::vector<double> &traj_times, const std::vector<int> &traj_replicas,
    std::fstream &log) {
  std::ofstream rep_out(filename);
  if (!rep_out.is_open()) {
    TitanUI::print_error(label + ": cannot open " + filename +
                         " for writing basin metadata.");
    return false;
  }

  rep_out << "# Representative frame = closest frame to basin minimum in "
          << space_name << " space\n";
  rep_out << "# Not necessarily the lowest-energy frame in this basin\n";
  rep_out << "ClusterID,Population%,FrameCount,WeightSum,DeltaG_Sink,"
             "Representative_Frame,ReplicaID,Time_ps,PC1,PC2\n";
  for (const auto &cluster : clusters) {
    const double reported_percentage =
        use_weights ? cluster.percentage_w : cluster.percentage;
    int replica_id = -1;
    double time_ps = 0.0;
    const int rep_idx = cluster.representative_frame_idx;
    if (rep_idx >= 0 &&
        static_cast<size_t>(rep_idx) < traj_replicas.size()) {
      replica_id = traj_replicas[static_cast<size_t>(rep_idx)];
    }
    if (rep_idx >= 0 && static_cast<size_t>(rep_idx) < traj_times.size()) {
      time_ps = traj_times[static_cast<size_t>(rep_idx)];
    }

    rep_out << cluster.id << "," << std::fixed << std::setprecision(2)
            << reported_percentage << "%," << cluster.population << ","
            << std::setprecision(6) << cluster.population_w << ","
            << std::setprecision(4) << cluster.min_energy << ","
            << cluster.representative_frame_idx << "," << replica_id << ","
            << std::setprecision(1) << time_ps << "," << cluster.min_pc1 << ","
            << cluster.min_pc2 << "\n";
  }
  TitanUI::log_kv(log, label + " Basin File", filename);
  return true;
}

void write_pca_extraction_script(const arguments &args,
                                 const std::vector<TitanPCA::ClusterNode> &clusters,
                                 const std::vector<double> &traj_times,
                                 const std::vector<int> &traj_replicas,
                                 const std::vector<std::string> &replica_traj_files) {
  std::ofstream bash_out("extract_basins.sh");
  if (!bash_out.is_open()) {
    TitanUI::print_warning("PCA: failed to open extract_basins.sh for writing.");
    return;
  }

  bash_out << "#!/bin/bash\n";
  bash_out << "# Titan Engine: Auto-generated FES Basin Extraction Script\n";
  bash_out << "# Representative frame = closest frame to basin minimum in PCA "
              "space\n";
  bash_out << "# Not necessarily the lowest-energy frame in this basin\n";
  bash_out << "# Review TOPOLOGY, trajectory paths, and trjconv group selection "
              "before running.\n";
  bash_out << "# Default group selection below uses index group 0; edit the "
              "heredoc if needed.\n\n";

  const std::string topology_hint =
      !args.pbc.ref_pdb_name.empty()
          ? args.pbc.ref_pdb_name
          : (!args.pdb_name.empty() ? args.pdb_name : "UPDATE_ME_STRUCTURE.pdb");
  bash_out << "TOPOLOGY=" << bash_quote(topology_hint) << "\n";

  std::set<int> used_replicas;
  for (const auto &cluster : clusters) {
    const int rep_idx = cluster.representative_frame_idx;
    if (rep_idx >= 0 && static_cast<size_t>(rep_idx) < traj_replicas.size()) {
      used_replicas.insert(traj_replicas[static_cast<size_t>(rep_idx)]);
    }
  }

  if (!used_replicas.empty()) {
    bash_out << "# Trajectory files inferred from Titan input:\n";
    for (int replica_id : used_replicas) {
      std::string traj_value = "replica_" + std::to_string(replica_id) + ".xtc";
      if (replica_id >= 0 &&
          static_cast<size_t>(replica_id) < replica_traj_files.size() &&
          !replica_traj_files[static_cast<size_t>(replica_id)].empty()) {
        traj_value = replica_traj_files[static_cast<size_t>(replica_id)];
      }
      bash_out << "TRAJ_REP_" << replica_id << "=" << bash_quote(traj_value)
               << "\n";
    }
  } else {
    bash_out << "# No replica metadata was available; update trajectory paths "
                "manually if you add extraction commands below.\n";
  }
  bash_out << "\n";

  for (const auto &cluster : clusters) {
    const int rep_idx = cluster.representative_frame_idx;
    int replica_id = -1;
    double time_ps = 0.0;
    if (rep_idx >= 0 && static_cast<size_t>(rep_idx) < traj_replicas.size()) {
      replica_id = traj_replicas[static_cast<size_t>(rep_idx)];
    }
    if (rep_idx >= 0 && static_cast<size_t>(rep_idx) < traj_times.size()) {
      time_ps = traj_times[static_cast<size_t>(rep_idx)];
    }
    if (replica_id < 0) {
      bash_out << "# Basin " << cluster.id
               << ": skipped automatic extraction because frame metadata is "
                  "unavailable.\n\n";
      continue;
    }

    const double tol_ps = 1.0;
    const std::string time_begin_str = format_time_ps(time_ps - tol_ps, 3);
    const std::string time_end_str = format_time_ps(time_ps + tol_ps, 3);
    const std::string out_name = "FES_Basin_" + std::to_string(cluster.id) +
                                 "_Rep" + std::to_string(replica_id) + "_" +
                                 make_time_label(time_ps) + "ps.pdb";
    bash_out << "echo \"-> Extracting Basin " << cluster.id << " (Replica "
             << replica_id << " at " << format_time_ps(time_ps, 1)
             << " ps)...\"\n";
    bash_out << "gmx trjconv -f $TRAJ_REP_" << replica_id
             << " -s $TOPOLOGY -b " << time_begin_str << " -e "
             << time_end_str << " -o " << bash_quote(out_name)
             << " <<< \"0\" > /dev/null 2>&1\n\n";
  }
  bash_out.close();

#ifndef _WIN32
  chmod("extract_basins.sh",
        S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
#endif

  TitanUI::print_info("PCA: Generated 'extract_basins.sh'.");
  TitanUI::print_info(
      "PCA: Review TOPOLOGY, trajectory paths, and trjconv group selection "
      "before running the script.");
  TitanUI::print_info(
      "PCA: Run the script in an environment where GROMACS 'gmx' is already "
      "available in PATH.");
}

} // namespace

namespace TitanPCA {

void execute_pca_fes_analysis(
    const arguments &args, const CoordinateRows &traj_x,
    const CoordinateRows &traj_y, const CoordinateRows &traj_z,
    const std::vector<std::string> &ca_resname,
    const std::vector<int> &ca_resnrs, const std::vector<std::string> &ca_chain,
    const std::vector<int> &ca_is_lig, TitanExecutor &executor, info &info_out,
    std::fstream &log,
    const std::vector<double> &frame_weights,
    const std::vector<double> &traj_times,
    const std::vector<int> &traj_replicas,
    const std::vector<std::string> &replica_traj_files,
    FesPostprocMetrics *metrics) {

  (void)ca_resname;
  (void)ca_resnrs;
  (void)ca_chain;
  (void)ca_is_lig;
  if (metrics)
    *metrics = {};

  if (traj_x.empty() || traj_x[0].empty()) {
    TitanUI::print_error("PCA ERROR: No valid trace coordinates found!");
    return;
  }

  const size_t num_frames = traj_x.size();
  const size_t num_ca = traj_x[0].size();
  const size_t dimensions = num_ca * 3;
  if (traj_y.size() != num_frames || traj_z.size() != num_frames) {
    TitanUI::print_error(
        "PCA ERROR: DCCM trajectory component frame counts do not match.");
    return;
  }
  for (size_t f = 0; f < num_frames; ++f) {
    if (traj_x[f].size() != num_ca || traj_y[f].size() != num_ca ||
        traj_z[f].size() != num_ca) {
      TitanUI::print_error(
          "PCA ERROR: DCCM trajectory component atom counts do not match.");
      return;
    }
  }
  const bool have_frame_metadata =
      (traj_times.size() == num_frames && traj_replicas.size() == num_frames);

  if ((!traj_times.empty() || !traj_replicas.empty()) && !have_frame_metadata) {
    TitanUI::print_warning(
        "PCA: frame metadata size mismatch. Representative basin extraction "
        "metadata will be incomplete.");
  }

  const auto weight_state =
      resolve_frame_weights("PCA", num_frames, frame_weights);
  const bool use_weights = weight_state.use_weights;
  const auto &analysis_weights = weight_state.analysis_weights;
  const double analysis_weight_sum = weight_state.analysis_weight_sum;

  auto t_start = std::chrono::steady_clock::now();
  TitanUI::log_kv(log, "PCA Snapshot Trace", std::to_string(num_frames));
  TitanUI::log_kv(log, "PCA Degrees of Freedom", std::to_string(dimensions));
  TitanUI::log_kv(log, "PCA Weight Mode",
                  use_weights ? "DHAM-weighted" : "uniform");

  FlatCoordinateMatrix aligned_traj(num_frames, dimensions, 0.0);
  for (size_t f = 0; f < num_frames; ++f) {
    const auto &x_row = traj_x[f];
    const auto &y_row = traj_y[f];
    const auto &z_row = traj_z[f];
    double *aligned_row = aligned_traj.row_data(f);
    for (size_t i = 0; i < num_ca; ++i) {
      const size_t base = i * 3;
      aligned_row[base + 0] = x_row[i];
      aligned_row[base + 1] = y_row[i];
      aligned_row[base + 2] = z_row[i];
    }
  }

  // 1. Setup Reference Frame (Frame 0)
  Eigen::MatrixXd ref_mat(num_ca, 3);
  Eigen::RowVector3d ref_centroid(0, 0, 0);
  const double *ref_row = aligned_traj.row_data(0);
  for (size_t i = 0; i < num_ca; ++i) {
    const size_t base = i * 3;
    ref_centroid +=
        Eigen::RowVector3d(ref_row[base + 0], ref_row[base + 1], ref_row[base + 2]);
  }
  ref_centroid /= num_ca;

  for (size_t i = 0; i < num_ca; ++i) {
    const size_t base = i * 3;
    ref_mat(i, 0) = ref_row[base + 0] - ref_centroid(0);
    ref_mat(i, 1) = ref_row[base + 1] - ref_centroid(1);
    ref_mat(i, 2) = ref_row[base + 2] - ref_centroid(2);
    aligned_traj(0, base + 0) = ref_mat(i, 0);
    aligned_traj(0, base + 1) = ref_mat(i, 1);
    aligned_traj(0, base + 2) = ref_mat(i, 2);
  }

  // 2. Kabsch Alignment for all frames using shared TitanAlignment module
  //    align_frame() centers mobile and rotates in-place; ref_mat stays
  //    unchanged.
  Eigen::MatrixXd frame_mat(num_ca, 3);
  for (size_t f = 1; f < num_frames; ++f) {
    const double *frame_row = aligned_traj.row_data(f);
    for (size_t i = 0; i < num_ca; ++i) {
      const size_t base = i * 3;
      frame_mat(i, 0) = frame_row[base + 0];
      frame_mat(i, 1) = frame_row[base + 1];
      frame_mat(i, 2) = frame_row[base + 2];
    }
    TitanAlignment::align_frame(ref_mat,
                                frame_mat); // centers + rotates in-place
    for (size_t i = 0; i < num_ca; ++i) {
      const size_t base = i * 3;
      aligned_traj(f, base + 0) = frame_mat(i, 0);
      aligned_traj(f, base + 1) = frame_mat(i, 1);
      aligned_traj(f, base + 2) = frame_mat(i, 2);
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
  for (size_t f = 0; f < num_frames; ++f) {
    const double *coords = aligned_traj.row_data(f);
    const double w = analysis_weights[f];
    for (size_t dim = 0; dim < dimensions; ++dim)
      mean_coords[dim] += w * coords[dim];
  }
  for (size_t i = 0; i < dimensions; ++i) {
    mean_coords[i] /= analysis_weight_sum;
  }

  // 4. Compute Covariance Matrix (Flattened 3N x 3N natively)
  // DHAM frame weights behave like a probability mass over frames, so the
  // weighted covariance uses sum(w) as its normalization. The legacy
  // unweighted path keeps the original (N - 1) scaling.
  Eigen::MatrixXd cov_matrix = Eigen::MatrixXd::Zero(dimensions, dimensions);
  const double cov_norm =
      use_weights ? analysis_weight_sum
                  : static_cast<double>(std::max<size_t>(num_frames - 1, 1));
  std::vector<double> dev(dimensions);

  for (size_t f = 0; f < num_frames; ++f) {
    const double *coords = aligned_traj.row_data(f);
    const double w = analysis_weights[f];
    for (size_t dim = 0; dim < dimensions; ++dim)
      dev[dim] = coords[dim] - mean_coords[dim];

    for (size_t i = 0; i < dimensions; ++i) {
      for (size_t j = i; j < dimensions; ++j) {
        double val = (w * dev[i] * dev[j]) / cov_norm;
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

  for (size_t f = 0; f < num_frames; ++f) {
    double p1 = 0.0, p2 = 0.0;
    const double *coords = aligned_traj.row_data(f);
    for (size_t dim = 0; dim < dimensions; ++dim) {
      const double diff = coords[dim] - mean_coords[dim];
      p1 += diff * pc1_vec(static_cast<Eigen::Index>(dim));
      p2 += diff * pc2_vec(static_cast<Eigen::Index>(dim));
    }
    projections[f] = {p1, p2};
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

  const double target_temp = resolve_target_temperature(args, use_weights);
  auto fes_result = compute_2d_fes_grid(
      projections, analysis_weights, args.fes_grid, target_temp,
      1.0 / args.eu_par, executor);
  if (!write_2d_fes_grid("PCA", info_out.get_pca_fes_2d_filename(), fes_result,
                         log)) {
    return;
  }

  const auto valid_clusters =
      detect_fes_clusters(fes_result, projections, analysis_weights, use_weights,
                          executor);
  const std::vector<double> empty_times;
  const std::vector<int> empty_replicas;
  write_cluster_representatives_file(
      "PCA", "PCA", info_out.get_cluster_representatives_filename(),
      valid_clusters, use_weights, have_frame_metadata ? traj_times : empty_times,
      have_frame_metadata ? traj_replicas : empty_replicas, log);
  write_pca_extraction_script(args, valid_clusters, traj_times, traj_replicas,
                              replica_traj_files);
  TitanUI::log_kv(log, "Macrostate Clusters Extracted",
                  std::to_string(valid_clusters.size()));

  if (args.zts.flag_zts && args.zts.zts_domain == "pca") {
    auto start_zts = std::chrono::steady_clock::now();
    TitanUI::print_module_header("Zero-Temperature String Method (PCA)");
    const auto zts_result =
        TitanZTS::run_zts(args.zts, fes_result, valid_clusters, executor, log);
    const auto zts_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_zts);
    if (!zts_result.path.empty()) {
      TitanZTS::write_zts_path(info_out.get_pca_zts_path_filename(), zts_result);
      TitanUI::log_kv(log, "PCA ZTS Path File",
                      info_out.get_pca_zts_path_filename());
      TitanUI::log_kv(log, "PCA ZTS Barrier Fwd",
                      std::to_string(zts_result.barrier_fwd) + " " + args.eu);
      TitanUI::log_kv(log, "PCA ZTS Barrier Rev",
                      std::to_string(zts_result.barrier_rev) + " " + args.eu);
      TitanUI::log_kv(log, "PCA ZTS Converged",
                      zts_result.converged ? "yes" : "no (max iter reached)");
      if (metrics) {
        metrics->zts_ran = true;
        metrics->zts_time = zts_time;
      }
    }
    TitanUI::log_time(log, "PCA ZTS", zts_time.count() / 1000.0);
  }

  auto t_end = std::chrono::steady_clock::now();
  TitanUI::log_time(
      log, "PCA + FES Clustering",
      std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_eigen)
              .count() /
          1000.0);
}

void execute_tica_fes_analysis(
    const arguments &args, const MatrixXdAnyConstRef &tica_proj_mat,
    int64_t n_actual_frames, const std::vector<int64_t> &frame_global_indices,
    TitanExecutor &executor, info &info_out, std::fstream &log,
    const std::vector<double> &frame_weights, FesPostprocMetrics *metrics) {
  if (metrics)
    *metrics = {};
  if (n_actual_frames <= 0) {
    TitanUI::print_warning("TICA: no projected frames available for 2D FES.");
    return;
  }
  if (tica_proj_mat.cols() < 2) {
    TitanUI::print_warning(
        "TICA: fewer than 2 projection components are available. "
        "Skipping 2D FES.");
    return;
  }

  const int64_t clamped_frames =
      std::min<int64_t>(n_actual_frames, tica_proj_mat.rows());
  if (clamped_frames <= 0) {
    TitanUI::print_warning("TICA: projection matrix is empty. Skipping 2D FES.");
    return;
  }

  const size_t num_frames = static_cast<size_t>(clamped_frames);
  std::vector<std::pair<double, double>> projections(num_frames);
  for (size_t frame = 0; frame < num_frames; ++frame) {
    projections[frame] = {tica_proj_mat(static_cast<Eigen::Index>(frame), 0),
                          tica_proj_mat(static_cast<Eigen::Index>(frame), 1)};
  }

  const auto local_frame_weights = map_global_weights_to_local_frames(
      "TICA", clamped_frames, frame_global_indices, frame_weights);
  const auto weight_state =
      resolve_frame_weights("TICA", num_frames, local_frame_weights);
  const bool use_weights = weight_state.use_weights;

  TitanUI::log_kv(log, "TICA FES Frames", std::to_string(num_frames));
  TitanUI::log_kv(log, "TICA Weight Mode",
                  use_weights ? "DHAM-weighted" : "uniform");

  const double target_temp = resolve_target_temperature(args, use_weights);
  auto fes_result = compute_2d_fes_grid(
      projections, weight_state.analysis_weights, args.fes_grid, target_temp,
      1.0 / args.eu_par, executor);
  if (!write_2d_fes_grid("TICA", info_out.get_tica_fes_2d_filename(), fes_result,
                         log)) {
    return;
  }

  const auto valid_clusters =
      detect_fes_clusters(fes_result, projections, weight_state.analysis_weights,
                          use_weights, executor);
  const std::vector<double> empty_times;
  const std::vector<int> empty_replicas;
  write_cluster_representatives_file(
      "TICA", "TICA", info_out.get_tica_cluster_representatives_filename(),
      valid_clusters, use_weights, empty_times, empty_replicas, log);
  TitanUI::log_kv(log, "TICA Macrostate Clusters",
                  std::to_string(valid_clusters.size()));

  if (args.zts.flag_zts && args.zts.zts_domain == "tica") {
    auto start_zts = std::chrono::steady_clock::now();
    TitanUI::print_module_header("Zero-Temperature String Method (TICA)");
    const auto zts_result =
        TitanZTS::run_zts(args.zts, fes_result, valid_clusters, executor, log);
    const auto zts_time = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_zts);
    if (!zts_result.path.empty()) {
      TitanZTS::write_zts_path(info_out.get_tica_zts_path_filename(), zts_result);
      TitanUI::log_kv(log, "TICA ZTS Path File",
                      info_out.get_tica_zts_path_filename());
      TitanUI::log_kv(log, "TICA ZTS Barrier Fwd",
                      std::to_string(zts_result.barrier_fwd) + " " + args.eu);
      TitanUI::log_kv(log, "TICA ZTS Barrier Rev",
                      std::to_string(zts_result.barrier_rev) + " " + args.eu);
      TitanUI::log_kv(log, "TICA ZTS Converged",
                      zts_result.converged ? "yes" : "no (max iter reached)");
      if (metrics) {
        metrics->zts_ran = true;
        metrics->zts_time = zts_time;
      }
    }
    TitanUI::log_time(log, "TICA ZTS", zts_time.count() / 1000.0);
  }
}

} // namespace TitanPCA
