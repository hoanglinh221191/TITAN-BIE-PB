#include "dccm_network.h"
#include "simd_compat.h"
#include "titan_alignment.h"
#include "titan_neighbor_grid.h"
#include "titan_sparse_graph.h"
#include "titan_ui.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <thread>


#ifdef __AVX__
#include <immintrin.h>
#endif

namespace {

size_t coord_count(const MOL &mol) noexcept {
  return std::min({mol.pos.x.size(), mol.pos.y.size(), mol.pos.z.size()});
}

inline bool is_aligned_32(const void *ptr) noexcept {
  return (reinterpret_cast<std::uintptr_t>(ptr) & 31u) == 0u;
}

template <typename T> struct FlatAlignedMatrixT {
  static constexpr size_t kSimdLaneWidth =
      std::max<size_t>(1, kTitanSimdAlignment / sizeof(T));

  FlatAlignedMatrixT() = default;
  FlatAlignedMatrixT(size_t rows, size_t cols, T value = T{})
      : rows_(rows), cols_(cols),
        stride_(cols == 0 ? 0
                          : ((cols + kSimdLaneWidth - 1) / kSimdLaneWidth) *
                                kSimdLaneWidth),
        data_(rows * stride_, value) {}

  void assign(size_t rows, size_t cols, T value = T{}) {
    rows_ = rows;
    cols_ = cols;
    stride_ = cols == 0 ? 0
                        : ((cols + kSimdLaneWidth - 1) / kSimdLaneWidth) *
                              kSimdLaneWidth;
    data_.assign(rows * stride_, value);
  }

  T &operator()(size_t row, size_t col) { return data_[row * stride_ + col]; }

  const T &operator()(size_t row, size_t col) const {
    return data_[row * stride_ + col];
  }

  T *row_data(size_t row) { return data_.data() + row * stride_; }
  const T *row_data(size_t row) const { return data_.data() + row * stride_; }

  size_t rows() const noexcept { return rows_; }
  size_t cols() const noexcept { return cols_; }
  size_t stride() const noexcept { return stride_; }

  void fill(T value) { std::fill(data_.begin(), data_.end(), value); }

private:
  size_t rows_ = 0;
  size_t cols_ = 0;
  size_t stride_ = 0;
  AlignedVector<T> data_;
};

using FlatAlignedMatrix = FlatAlignedMatrixT<double>;

inline int resolve_allostery_bin_count(size_t total_frames, int num_bins) {
  if (num_bins > 0)
    return num_bins;
  return std::max(8, static_cast<int>(std::cbrt(total_frames / 5.0)));
}

inline std::vector<std::string>
build_allostery_node_labels(const std::vector<std::string> &ca_resname,
                            const std::vector<int> &ca_resnrs,
                            const std::vector<std::string> &ca_chain) {
  const size_t num_ca = ca_resname.size();
  std::vector<std::string> node_labels(num_ca);
  for (size_t i = 0; i < num_ca; ++i)
    node_labels[i] =
        ca_resname[i] + std::to_string(ca_resnrs[i]) + "_" + ca_chain[i];
  return node_labels;
}

inline size_t upper_triangle_pair_index(size_t i, size_t j,
                                        size_t node_count) noexcept {
  return i * (2 * node_count - i - 1) / 2 + (j - i - 1);
}

std::vector<std::pair<size_t, size_t>>
build_dense_upper_triangle_pairs(size_t node_count) {
  std::vector<std::pair<size_t, size_t>> pairs;
  if (node_count < 2)
    return pairs;

  pairs.reserve(node_count * (node_count - 1) / 2);
  for (size_t i = 0; i < node_count; ++i) {
    for (size_t j = i + 1; j < node_count; ++j)
      pairs.emplace_back(i, j);
  }
  return pairs;
}

std::vector<std::pair<size_t, size_t>> build_lspmd_candidate_pairs(
    const DCCMNetwork::TrajectoryRows &traj_x,
    const DCCMNetwork::TrajectoryRows &traj_y,
    const DCCMNetwork::TrajectoryRows &traj_z, size_t num_ca,
    size_t total_frames, double cutoff, bool use_avg_cutoff) {
  std::vector<std::pair<size_t, size_t>> pairs;
  if (num_ca < 2 || total_frames == 0 || cutoff <= 0.0)
    return pairs;

  const size_t dense_pair_count = num_ca * (num_ca - 1) / 2;
  const size_t frame_limit = use_avg_cutoff ? total_frames : 1;
  std::vector<unsigned char> seen(dense_pair_count, 0);
  titan_neighbor::HalfShellGrid neighbor_grid;
  pairs.reserve(std::min(dense_pair_count, num_ca * static_cast<size_t>(32)));

  for (size_t f = 0; f < frame_limit; ++f) {
    if (traj_x[f].size() < num_ca || traj_y[f].size() < num_ca ||
        traj_z[f].size() < num_ca) {
      return build_dense_upper_triangle_pairs(num_ca);
    }

    const auto &x_row = traj_x[f];
    const auto &y_row = traj_y[f];
    const auto &z_row = traj_z[f];
    const auto [xmin, xmax] =
        std::minmax_element(x_row.begin(), x_row.begin() + num_ca);
    const auto [ymin, ymax] =
        std::minmax_element(y_row.begin(), y_row.begin() + num_ca);
    const auto [zmin, zmax] =
        std::minmax_element(z_row.begin(), z_row.begin() + num_ca);

    titan_neighbor::OrthorhombicGridSpec spec;
    spec.cutoff = cutoff;
    spec.shift_x = -(*xmin) + cutoff;
    spec.shift_y = -(*ymin) + cutoff;
    spec.shift_z = -(*zmin) + cutoff;
    spec.len_x = (*xmax - *xmin) + 2.0 * cutoff;
    spec.len_y = (*ymax - *ymin) + 2.0 * cutoff;
    spec.len_z = (*zmax - *zmin) + 2.0 * cutoff;
    spec.periodic = false;

    neighbor_grid.for_each_pair(
        num_ca, spec,
        [&](std::size_t atom_index, double &x, double &y, double &z) {
          x = x_row[atom_index];
          y = y_row[atom_index];
          z = z_row[atom_index];
        },
        [&](std::size_t atom_index) { return atom_index < num_ca; },
        [&](std::size_t lhs, std::size_t rhs) {
          const double dx = x_row[lhs] - x_row[rhs];
          const double dy = y_row[lhs] - y_row[rhs];
          const double dz = z_row[lhs] - z_row[rhs];
          if (dx * dx + dy * dy + dz * dz > cutoff * cutoff)
            return;

          const size_t i = std::min(lhs, rhs);
          const size_t j = std::max(lhs, rhs);
          const size_t pair_index = upper_triangle_pair_index(i, j, num_ca);
          if (seen[pair_index] != 0)
            return;
          seen[pair_index] = 1;
          pairs.emplace_back(i, j);
        });
  }

  return pairs;
}

std::vector<std::pair<size_t, size_t>> build_allostery_candidate_pairs(
    const DCCMNetwork::TrajectoryRows &traj_x,
    const DCCMNetwork::TrajectoryRows &traj_y,
    const DCCMNetwork::TrajectoryRows &traj_z, size_t num_ca,
    size_t total_frames, const DCCMNetwork::AllosteryConfig &cfg) {
  const size_t dense_pair_count =
      num_ca < 2 ? 0 : num_ca * (num_ca - 1) / 2;
  if (cfg.cutoff_ang <= 0.0 && cfg.knn <= 0)
    return build_dense_upper_triangle_pairs(num_ca);

  std::vector<std::pair<size_t, size_t>> pairs;
  std::vector<unsigned char> seen(dense_pair_count, 0);
  auto mark_pair = [&](size_t i, size_t j) {
    const size_t a = std::min(i, j);
    const size_t b = std::max(i, j);
    const size_t pair_index = upper_triangle_pair_index(a, b, num_ca);
    if (seen[pair_index] != 0)
      return false;
    seen[pair_index] = 1;
    pairs.emplace_back(a, b);
    return true;
  };

  if (cfg.cutoff_ang > 0.0) {
    auto cutoff_pairs = build_lspmd_candidate_pairs(
        traj_x, traj_y, traj_z, num_ca, total_frames, cfg.cutoff_ang,
        cfg.use_avg_cutoff);
    for (const auto &[i, j] : cutoff_pairs)
      mark_pair(i, j);
  }

  if (cfg.knn > 0 && num_ca > 1 && total_frames > 0) {
    const size_t k = std::min<size_t>(static_cast<size_t>(cfg.knn), num_ca - 1);
    const size_t frame_limit = cfg.use_avg_cutoff ? total_frames : 1;
    bool can_rank = true;
    for (size_t f = 0; f < frame_limit; ++f) {
      if (traj_x[f].size() < num_ca || traj_y[f].size() < num_ca ||
          traj_z[f].size() < num_ca) {
        can_rank = false;
        break;
      }
    }
    if (!can_rank) {
      if (pairs.empty())
        return build_dense_upper_triangle_pairs(num_ca);
      std::sort(pairs.begin(), pairs.end());
      return pairs;
    }

    struct RankedNeighbor {
      size_t partner = 0;
      std::uint32_t votes = 0;
      float min_dist2 = std::numeric_limits<float>::infinity();
    };

    std::vector<std::uint32_t> vote_counts(dense_pair_count, 0);
    std::vector<float> min_dist2(dense_pair_count,
                                 std::numeric_limits<float>::infinity());
    for (size_t f = 0; f < frame_limit; ++f) {
      const auto &x_row = traj_x[f];
      const auto &y_row = traj_y[f];
      const auto &z_row = traj_z[f];

      std::vector<std::pair<float, size_t>> distances;
      distances.reserve(num_ca - 1);
      auto rank_cmp = [](const auto &lhs, const auto &rhs) {
        if (lhs.first != rhs.first)
          return lhs.first < rhs.first;
        return lhs.second < rhs.second;
      };

      for (size_t i = 0; i < num_ca; ++i) {
        distances.clear();
        const double xi = x_row[i];
        const double yi = y_row[i];
        const double zi = z_row[i];
        for (size_t j = 0; j < num_ca; ++j) {
          if (i == j)
            continue;
          const double dx = xi - x_row[j];
          const double dy = yi - y_row[j];
          const double dz = zi - z_row[j];
          distances.emplace_back(static_cast<float>(dx * dx + dy * dy + dz * dz),
                                 j);
        }

        if (distances.size() > k) {
          std::nth_element(distances.begin(),
                           distances.begin() +
                               static_cast<std::ptrdiff_t>(k),
                           distances.end(), rank_cmp);
          std::sort(distances.begin(),
                    distances.begin() + static_cast<std::ptrdiff_t>(k),
                    rank_cmp);
        } else {
          std::sort(distances.begin(), distances.end(), rank_cmp);
        }

        for (size_t n = 0; n < k; ++n) {
          const size_t partner = distances[n].second;
          const size_t a = std::min(i, partner);
          const size_t b = std::max(i, partner);
          const size_t pair_index = upper_triangle_pair_index(a, b, num_ca);
          if (vote_counts[pair_index] !=
              std::numeric_limits<std::uint32_t>::max()) {
            ++vote_counts[pair_index];
          }
          min_dist2[pair_index] =
              std::min(min_dist2[pair_index], distances[n].first);
        }
      }
    }

    std::vector<std::vector<RankedNeighbor>> ranked_neighbors(num_ca);
    for (auto &neighbors : ranked_neighbors)
      neighbors.reserve(k * 2);
    for (size_t i = 0; i < num_ca; ++i) {
      for (size_t j = i + 1; j < num_ca; ++j) {
        const size_t pair_index = upper_triangle_pair_index(i, j, num_ca);
        if (vote_counts[pair_index] == 0)
          continue;
        const RankedNeighbor ij{j, vote_counts[pair_index],
                                min_dist2[pair_index]};
        const RankedNeighbor ji{i, vote_counts[pair_index],
                                min_dist2[pair_index]};
        ranked_neighbors[i].push_back(ij);
        ranked_neighbors[j].push_back(ji);
      }
    }

    auto better_neighbor = [](const RankedNeighbor &lhs,
                              const RankedNeighbor &rhs) {
      if (lhs.votes != rhs.votes)
        return lhs.votes > rhs.votes;
      if (lhs.min_dist2 != rhs.min_dist2)
        return lhs.min_dist2 < rhs.min_dist2;
      return lhs.partner < rhs.partner;
    };

    for (size_t i = 0; i < num_ca; ++i) {
      auto &neighbors = ranked_neighbors[i];
      if (neighbors.empty())
        continue;
      if (neighbors.size() > k) {
        std::nth_element(neighbors.begin(),
                         neighbors.begin() + static_cast<std::ptrdiff_t>(k),
                         neighbors.end(), better_neighbor);
        std::sort(neighbors.begin(),
                  neighbors.begin() + static_cast<std::ptrdiff_t>(k),
                  better_neighbor);
      } else {
        std::sort(neighbors.begin(), neighbors.end(), better_neighbor);
      }
      const size_t keep = std::min(k, neighbors.size());
      for (size_t n = 0; n < keep; ++n)
        mark_pair(std::min(i, neighbors[n].partner),
                  std::max(i, neighbors[n].partner));
    }
  }

  std::sort(pairs.begin(), pairs.end());
  return pairs;
}

template <typename BinType>
size_t calculate_allostery_information_flat_impl(
    const DCCMNetwork::TrajectoryRows &traj_x,
    const DCCMNetwork::TrajectoryRows &traj_y,
    const DCCMNetwork::TrajectoryRows &traj_z,
    const std::vector<std::string> &node_labels, size_t total_frames,
    size_t num_threads, const DCCMNetwork::AllosteryConfig &cfg, int B,
    size_t *dense_pair_count_out = nullptr) {
  if (total_frames == 0 || node_labels.empty())
    return 0;
  const size_t num_ca = node_labels.size();
  const size_t tau = static_cast<size_t>(std::max(1, cfg.tau_lag));
  if (tau >= total_frames) {
    std::cerr << "[Allostery] tau_lag (" << tau << ") >= total_frames ("
              << total_frames << "). Skipping.\n";
    return 0;
  }

  const size_t B_sz = static_cast<size_t>(B);
  if (B_sz > static_cast<size_t>(std::numeric_limits<BinType>::max())) {
    TitanUI::print_warning(
        "Allostery bin count exceeds selected bin storage. Skipping.");
    return 0;
  }

  const size_t dense_pair_count = num_ca * (num_ca - 1) / 2;
  if (dense_pair_count_out != nullptr)
    *dense_pair_count_out = dense_pair_count;
  const auto pairs = build_allostery_candidate_pairs(traj_x, traj_y, traj_z,
                                                     num_ca, total_frames, cfg);
  const size_t n_pairs = pairs.size();
  std::ofstream gc_out("network_gc_edges.csv");
  gc_out << "Source,Target,GC\n";
  std::ofstream te_out("network_te_edges.csv");
  te_out << "Source,Target,Net_TE\n";
  if (n_pairs == 0)
    return 0;

  std::vector<double> mean_x(num_ca, 0.0), mean_y(num_ca, 0.0),
      mean_z(num_ca, 0.0);
  for (size_t f = 0; f < total_frames; ++f) {
    const auto &x_row = traj_x[f];
    const auto &y_row = traj_y[f];
    const auto &z_row = traj_z[f];
    for (size_t a = 0; a < num_ca; ++a) {
      mean_x[a] += x_row[a];
      mean_y[a] += y_row[a];
      mean_z[a] += z_row[a];
    }
  }
  const double inv_f = 1.0 / static_cast<double>(total_frames);
  for (size_t a = 0; a < num_ca; ++a) {
    mean_x[a] *= inv_f;
    mean_y[a] *= inv_f;
    mean_z[a] *= inv_f;
  }

  FlatAlignedMatrixT<BinType> proj_bins(num_ca, total_frames,
                                        static_cast<BinType>(0));
  std::vector<double> H1(num_ca, 0.0);
  std::vector<double> H_past;
  std::vector<double> H_auto;
  if (cfg.compute_te) {
    H_past.assign(num_ca, 0.0);
    H_auto.assign(num_ca, 0.0);
  }
  const double norm1_full =
      0.5 * static_cast<double>(B) + static_cast<double>(total_frames);
  const double norm1_tau =
      0.5 * static_cast<double>(B) +
      static_cast<double>(total_frames - tau);
  const double norm2_full =
      0.5 * static_cast<double>(B_sz * B_sz) +
      static_cast<double>(total_frames);
  const double norm2_tau =
      0.5 * static_cast<double>(B_sz * B_sz) +
      static_cast<double>(total_frames - tau);
  const double norm3_tau =
      0.5 * static_cast<double>(B_sz * B_sz * B_sz) +
      static_cast<double>(total_frames - tau);

  auto entropy_from_hist = [](const std::vector<double> &hist,
                              double norm) -> double {
    double H = 0.0;
    for (double count : hist) {
      const double p = count / norm;
      if (p > 0.0)
        H -= p * std::log(p);
    }
    return H;
  };

  {
    const size_t atom_worker_count = std::max<size_t>(1, num_threads);
    std::vector<std::jthread> atom_workers;
    atom_workers.reserve(atom_worker_count);
    const size_t atoms_per_worker =
        (num_ca + atom_worker_count - 1) / atom_worker_count;
    for (size_t t = 0; t < atom_worker_count; ++t) {
      const size_t a_begin = t * atoms_per_worker;
      const size_t a_end = std::min(a_begin + atoms_per_worker, num_ca);
      if (a_begin >= a_end)
        continue;
      atom_workers.emplace_back([&, a_begin, a_end]() {
        AlignedDoubleVector atom_proj(total_frames, 0.0);
        std::vector<double> hist1(B_sz, 0.0);
        std::vector<double> hist2(cfg.compute_te ? B_sz * B_sz : 0, 0.0);
        for (size_t a = a_begin; a < a_end; ++a) {
          double c_xx = 0.0, c_xy = 0.0, c_xz = 0.0;
          double c_yy = 0.0, c_yz = 0.0, c_zz = 0.0;
          const double avg_x = mean_x[a];
          const double avg_y = mean_y[a];
          const double avg_z = mean_z[a];
          for (size_t f = 0; f < total_frames; ++f) {
            const double dx = traj_x[f][a] - avg_x;
            const double dy = traj_y[f][a] - avg_y;
            const double dz = traj_z[f][a] - avg_z;
            c_xx += dx * dx;
            c_xy += dx * dy;
            c_xz += dx * dz;
            c_yy += dy * dy;
            c_yz += dy * dz;
            c_zz += dz * dz;
          }

          Eigen::Matrix3d cov3;
          cov3 << c_xx, c_xy, c_xz, c_xy, c_yy, c_yz, c_xz, c_yz, c_zz;
          cov3 /= static_cast<double>(total_frames);

          Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov3);
          const Eigen::Vector3d e1 = es.eigenvectors().col(2);
          const double ex = e1(0);
          const double ey = e1(1);
          const double ez = e1(2);

          double p_min = std::numeric_limits<double>::max();
          double p_max = std::numeric_limits<double>::lowest();
          for (size_t f = 0; f < total_frames; ++f) {
            const double proj_val = (traj_x[f][a] - avg_x) * ex +
                                    (traj_y[f][a] - avg_y) * ey +
                                    (traj_z[f][a] - avg_z) * ez;
            atom_proj[f] = proj_val;
            p_min = std::min(p_min, proj_val);
            p_max = std::max(p_max, proj_val);
          }

          double p_range = p_max - p_min;
          if (p_range < 1e-15)
            p_range = 1e-15;

          BinType *bins = proj_bins.row_data(a);
          std::fill(hist1.begin(), hist1.end(), 0.5);
          for (size_t f = 0; f < total_frames; ++f) {
            const double scaled = (atom_proj[f] - p_min) / p_range;
            int bin = static_cast<int>(scaled * static_cast<double>(B));
            if (bin < 0)
              bin = 0;
            else if (bin >= B)
              bin = B - 1;
            bins[f] = static_cast<BinType>(bin);
            hist1[static_cast<size_t>(bins[f])] += 1.0;
          }
          H1[a] = entropy_from_hist(hist1, norm1_full);

          if (cfg.compute_te) {
            std::fill(hist1.begin(), hist1.end(), 0.5);
            for (size_t f = tau; f < total_frames; ++f)
              hist1[static_cast<size_t>(bins[f - tau])] += 1.0;
            H_past[a] = entropy_from_hist(hist1, norm1_tau);

            std::fill(hist2.begin(), hist2.end(), 0.5);
            for (size_t f = tau; f < total_frames; ++f) {
              const size_t b_now = static_cast<size_t>(bins[f]);
              const size_t b_past = static_cast<size_t>(bins[f - tau]);
              hist2[b_now * B_sz + b_past] += 1.0;
            }
            H_auto[a] = entropy_from_hist(hist2, norm2_tau);
          }
        }
      });
    }
  }

  std::vector<double> gc_vals(n_pairs, 0.0);
  std::vector<double> net_te_vals(cfg.compute_te ? n_pairs : 0, 0.0);

  auto worker = [&](size_t p_start, size_t p_end) {
    std::vector<double> hist2d(B_sz * B_sz, 0.0);
    std::vector<double> hist2d_b(cfg.compute_te ? B_sz * B_sz : 0, 0.0);
    std::vector<double> hist3d(cfg.compute_te ? B_sz * B_sz * B_sz : 0, 0.0);

    for (size_t p = p_start; p < p_end; ++p) {
      const size_t i = pairs[p].first;
      const size_t j = pairs[p].second;
      const BinType *bins_i = proj_bins.row_data(i);
      const BinType *bins_j = proj_bins.row_data(j);

      std::fill(hist2d.begin(), hist2d.end(), 0.5);
      if (cfg.compute_te)
        std::fill(hist2d_b.begin(), hist2d_b.end(), 0.5);
      for (size_t f = 0; f < total_frames; ++f) {
        const size_t bi = static_cast<size_t>(bins_i[f]);
        const size_t bj = static_cast<size_t>(bins_j[f]);
        hist2d[bi * B_sz + bj] += 1.0;
        if (cfg.compute_te && f >= tau) {
          const size_t bi_past = static_cast<size_t>(bins_i[f - tau]);
          const size_t bj_past = static_cast<size_t>(bins_j[f - tau]);
          hist2d_b[bi_past * B_sz + bj_past] += 1.0;
        }
      }

      const double H_ij = entropy_from_hist(hist2d, norm2_full);
      double MI = H1[i] + H1[j] - H_ij;
      if (MI < 0.0)
        MI = 0.0;
      gc_vals[p] = std::sqrt(1.0 - std::exp(-2.0 * MI));

      if (!cfg.compute_te)
        continue;

      const double H_past_ij = entropy_from_hist(hist2d_b, norm2_tau);
      auto calc_h3 = [&](const BinType *src_bins,
                         const BinType *dst_bins) -> double {
        std::fill(hist3d.begin(), hist3d.end(), 0.5);
        for (size_t f = tau; f < total_frames; ++f) {
          const size_t b1 = static_cast<size_t>(dst_bins[f]);
          const size_t b2 = static_cast<size_t>(dst_bins[f - tau]);
          const size_t b3 = static_cast<size_t>(src_bins[f - tau]);
          hist3d[(b1 * B_sz + b2) * B_sz + b3] += 1.0;
        }
        return entropy_from_hist(hist3d, norm3_tau);
      };

      const double H3_j_from_i = calc_h3(bins_i, bins_j);
      const double H3_i_from_j = calc_h3(bins_j, bins_i);
      const double te_ij = H_auto[j] + H_past_ij - H_past[j] - H3_j_from_i;
      const double te_ji = H_auto[i] + H_past_ij - H_past[i] - H3_i_from_j;
      net_te_vals[p] = te_ij - te_ji;
    }
  };

  const size_t worker_count = std::max<size_t>(1, num_threads);
  {
    std::vector<std::jthread> workers;
    workers.reserve(worker_count);
    const size_t chunk = (n_pairs + worker_count - 1) / worker_count;
    for (size_t t = 0; t < worker_count; ++t) {
      const size_t start = t * chunk;
      const size_t end = std::min(start + chunk, n_pairs);
      if (start < end)
        workers.emplace_back(worker, start, end);
    }
  }

  for (size_t p = 0; p < n_pairs; ++p) {
    if (gc_vals[p] > 0.35) {
      const size_t i = pairs[p].first;
      const size_t j = pairs[p].second;
      gc_out << node_labels[i] << "," << node_labels[j] << "," << gc_vals[p]
             << "\n";
    }
  }
  gc_out.close();

  if (cfg.compute_te) {
    for (size_t p = 0; p < n_pairs; ++p) {
      if (std::abs(net_te_vals[p]) > 0.05) {
        size_t drv = pairs[p].first;
        size_t tgt = pairs[p].second;
        if (net_te_vals[p] <= 0.0)
          std::swap(drv, tgt);
        te_out << node_labels[drv] << "," << node_labels[tgt] << ","
               << std::abs(net_te_vals[p]) << "\n";
      }
    }
  }
  te_out.close();
  return n_pairs;
}

} // namespace

// ── Helper: find Cβ (or N for Gly) atom index within a residue ────────────────
int DCCMNetwork::find_cb_for_residue(const Multi_MOL &mol, int resnr,
                                     const std::string &chain, bool is_gly) {
  const std::string target = is_gly ? "N" : "CB";
  for (size_t i = 0; i < mol.info.size(); ++i) {
    std::string at = mol.info[i].get_atom();
    // trim whitespace
    auto p1 = at.find_first_not_of(" \t");
    auto p2 = at.find_last_not_of(" \t");
    if (p1 != std::string::npos)
      at = at.substr(p1, p2 - p1 + 1);
    if (at == target &&
        mol.info[i].get_resnr() == resnr &&
        mol.info[i].get_chain() == chain)
      return static_cast<int>(i);
  }
  return -1; // not found
}

// Hàm dot_product tối ưu bằng tập lệnh SIMD AVX2
inline double fast_dot_product(const double *a, const double *b, size_t n) {
  double total = 0.0;
  size_t i = 0;
#ifdef __AVX__
  __m256d sum_vec = _mm256_setzero_pd();
  for (; i + 3 < n; i += 4) {
    // Dùng loadu (unaligned) an toàn cho std::vector
    __m256d va = _mm256_loadu_pd(&a[i]);
    __m256d vb = _mm256_loadu_pd(&b[i]);
    sum_vec = simd_compat::fmadd_pd(va, vb, sum_vec);
  }
  double buffer[4];
  _mm256_storeu_pd(buffer, sum_vec);
  total = buffer[0] + buffer[1] + buffer[2] + buffer[3];
#endif
  // Phần dư còn lại nếu n không chia hết cho 4
  for (; i < n; ++i) {
    total += a[i] * b[i];
  }
  return total;
}

inline double fast_dot_product_aligned(const double *a, const double *b,
                                       size_t n) {
  double total = 0.0;
  size_t i = 0;
#ifdef __AVX__
  __m256d sum_vec = _mm256_setzero_pd();
  const bool aligned = is_aligned_32(a) && is_aligned_32(b);
  for (; i + 3 < n; i += 4) {
    const __m256d va =
        aligned ? _mm256_load_pd(&a[i]) : _mm256_loadu_pd(&a[i]);
    const __m256d vb =
        aligned ? _mm256_load_pd(&b[i]) : _mm256_loadu_pd(&b[i]);
    sum_vec = simd_compat::fmadd_pd(va, vb, sum_vec);
  }
  alignas(32) double buffer[4];
  _mm256_store_pd(buffer, sum_vec);
  total = buffer[0] + buffer[1] + buffer[2] + buffer[3];
#endif
  for (; i < n; ++i) {
    total += a[i] * b[i];
  }
  return total;
}

void DCCMNetwork::init(const Multi_MOL &pro_in, const Multi_MOL &lig,
                       const arguments &args) {
  // CRITICAL FIX: Clear all state from previous trajectories (Batch Mode
  // Safety)
  ca_indices.clear();
  ca_resnrs.clear();
  ca_is_lig.clear();
  ca_chain.clear();
  ca_resname.clear();
  cb_indices.clear();
  traj_x.clear();
  traj_y.clear();
  traj_z.clear();
  traj_cbx.clear();
  traj_cby.clear();
  traj_cbz.clear();
  frame_global_indices_.clear();
  traj_times_.clear();
  traj_replicas_.clear();
  latest_traj_x_.clear();
  latest_traj_y_.clear();
  latest_traj_z_.clear();
  total_frames = 0;
  kabsch_aligned_ = false;
  store_full_traj_ =
      args.flag_network_explicit || args.flag_pca || args.flag_allostery ||
      args.lspmd.flag_lspmd;
  store_cb_traj_ = args.flag_allostery || args.lspmd.flag_lspmd;

  const size_t pro_coord_count =
      pro_in.m_mol.empty() ? 0 : coord_count(pro_in.m_mol.front());
  const size_t lig_coord_count =
      lig.m_mol.empty() ? 0 : coord_count(lig.m_mol.front());
  size_t skipped_rec_nodes = 0;
  size_t skipped_lig_nodes = 0;
  size_t cb_fallback_nodes = 0;
  const bool include_lig_trace_nodes = !args.flag_nocal;

  bool is_lig_rec = (args.mol_1_type == "ligand" || args.mol_1_type == "lig");
  bool is_lig_lg = (args.mol_2_type == "ligand" || args.mol_2_type == "lig");

  // Collect Receptor Nodes
  for (size_t i = 0; i < pro_in.info.size(); ++i) {
    if (i >= pro_coord_count) {
      ++skipped_rec_nodes;
      continue;
    }
    std::string atom_name = pro_in.info[i].get_atom();
    if (is_lig_rec || atom_name == "CA" || atom_name == "CA " ||
        atom_name == "BB" || atom_name == "BB2" || atom_name == "C4'" ||
        atom_name == "C4*") {
      ca_indices.emplace_back((int)i);
      ca_resnrs.emplace_back(is_lig_rec ? (i + 1) : pro_in.info[i].get_resnr());
      ca_is_lig.emplace_back(0); // 0 = Receptor
      ca_chain.emplace_back(pro_in.info[i].get_chain());
      ca_resname.emplace_back(is_lig_rec ? pro_in.info[i].get_atom()
                                         : pro_in.info[i].get_res_name());
    }
  }
  // Collect Ligand Nodes
  if (include_lig_trace_nodes) {
    for (size_t i = 0; i < lig.info.size(); ++i) {
      if (i >= lig_coord_count) {
        ++skipped_lig_nodes;
        continue;
      }
      std::string atom_name = lig.info[i].get_atom();
      if (is_lig_lg || atom_name == "CA" || atom_name == "CA " ||
          atom_name == "BB" || atom_name == "BB2" || atom_name == "C4'" ||
          atom_name == "C4*") {
        ca_indices.emplace_back((int)i);
        ca_resnrs.emplace_back(is_lig_lg ? (i + 1) : lig.info[i].get_resnr());
        ca_is_lig.emplace_back(1); // 1 = Ligand
        ca_chain.emplace_back(lig.info[i].get_chain());
        ca_resname.emplace_back(is_lig_lg ? lig.info[i].get_atom()
                                          : lig.info[i].get_res_name());
      }
    }
  }

  // Only allostery needs the heavier Cβ/N side trace cache.
  if (store_cb_traj_) {
    cb_indices.resize(ca_indices.size());
    for (size_t i = 0; i < ca_indices.size(); ++i) {
      std::string rn = ca_resname[i];
      auto p1 = rn.find_first_not_of(" \t");
      auto p2 = rn.find_last_not_of(" \t");
      if (p1 != std::string::npos)
        rn = rn.substr(p1, p2 - p1 + 1);
      bool is_gly = (rn == "GLY");

      const Multi_MOL &mol_ref = (ca_is_lig[i] == 0) ? pro_in : lig;
      const size_t coord_limit =
          (ca_is_lig[i] == 0) ? pro_coord_count : lig_coord_count;
      const int cb_idx =
          find_cb_for_residue(mol_ref, ca_resnrs[i], ca_chain[i], is_gly);
      if (cb_idx >= 0 && static_cast<size_t>(cb_idx) < coord_limit) {
        cb_indices[i] = cb_idx;
      } else {
        cb_indices[i] = -1;
        if (cb_idx >= 0) {
          ++cb_fallback_nodes;
        }
      }
    }
  } else {
    cb_indices.clear();
  }

  total_frames = 0;
  frame_global_indices_.clear();
  if (skipped_rec_nodes > 0 || skipped_lig_nodes > 0 || cb_fallback_nodes > 0) {
    std::ostringstream oss;
    oss << "DCCM/TICA bounds filter active:";
    if (skipped_rec_nodes > 0) {
      oss << " skipped " << skipped_rec_nodes
          << " receptor topology atoms without coordinates";
    }
    if (skipped_lig_nodes > 0) {
      oss << (skipped_rec_nodes > 0 ? "," : "") << " skipped "
          << skipped_lig_nodes
          << " ligand topology atoms without coordinates";
    }
    if (cb_fallback_nodes > 0) {
      oss << ((skipped_rec_nodes > 0 || skipped_lig_nodes > 0) ? "," : "")
          << " fallback-to-CA for " << cb_fallback_nodes
          << " missing Cbeta/N nodes";
    }
    TitanUI::print_warning(oss.str() + ".");
  }
  if (ca_indices.empty()) {
    TitanUI::print_warning(
        "DCCM/TICA: no valid trace nodes remain after coordinate-bounds "
        "filtering.");
  }
}

void DCCMNetwork::process_frame(const Multi_MOL &pro_in, const Multi_MOL &lig,
                                size_t chunk_snap_idx,
                                size_t global_frame_idx,
                                double time_ps,
                                int replica_id) {
  static bool warned_bad_frame = false;
  if (ca_indices.empty())
    return;

  if (chunk_snap_idx >= pro_in.m_mol.size()) {
    if (!warned_bad_frame) {
      TitanUI::print_warning(
          "DCCM/TICA: skipping frame because receptor frame index exceeds "
          "loaded chunk.");
      warned_bad_frame = true;
    }
    return;
  }

  AlignedDoubleVector cx(ca_indices.size(), 0.0),
      cy(ca_indices.size(), 0.0), cz(ca_indices.size(), 0.0);
  AlignedDoubleVector cbx, cby, cbz;
  if (store_cb_traj_) {
    cbx.assign(ca_indices.size(), 0.0);
    cby.assign(ca_indices.size(), 0.0);
    cbz.assign(ca_indices.size(), 0.0);
  }
  size_t l_snap = (lig.m_mol.size() > 1) ? chunk_snap_idx : 0;

  if (l_snap >= lig.m_mol.size()) {
    if (!warned_bad_frame) {
      TitanUI::print_warning(
          "DCCM/TICA: skipping frame because ligand frame index exceeds "
          "loaded chunk.");
      warned_bad_frame = true;
    }
    return;
  }

  const size_t pro_coord_count = coord_count(pro_in.m_mol[chunk_snap_idx]);
  const size_t lig_coord_count = coord_count(lig.m_mol[l_snap]);

  for (size_t i = 0; i < ca_indices.size(); ++i) {
    const size_t coord_limit =
        (ca_is_lig[i] == 0) ? pro_coord_count : lig_coord_count;
    const int idx = ca_indices[i];
    if (idx < 0 || static_cast<size_t>(idx) >= coord_limit) {
      if (!warned_bad_frame) {
        TitanUI::print_warning(
            "DCCM/TICA: skipping at least one frame because a trace atom index "
            "is outside the available coordinate arrays. This usually means "
            "topology and trajectory atom counts differ.");
        warned_bad_frame = true;
      }
      return;
    }
  }

  for (size_t i = 0; i < ca_indices.size(); ++i) {
    int idx = ca_indices[i];
    if (ca_is_lig[i] == 0) {
      cx[i] = pro_in.m_mol[chunk_snap_idx].pos.x[idx];
      cy[i] = pro_in.m_mol[chunk_snap_idx].pos.y[idx];
      cz[i] = pro_in.m_mol[chunk_snap_idx].pos.z[idx];
    } else {
      cx[i] = lig.m_mol[l_snap].pos.x[idx];
      cy[i] = lig.m_mol[l_snap].pos.y[idx];
      cz[i] = lig.m_mol[l_snap].pos.z[idx];
    }

    if (store_cb_traj_) {
      // Cβ/N coord: use cb_indices if valid, else fallback to CA
      int cb_idx = cb_indices[i];
      const size_t coord_limit =
          (ca_is_lig[i] == 0) ? pro_coord_count : lig_coord_count;
      if (cb_idx >= 0 && static_cast<size_t>(cb_idx) < coord_limit) {
        if (ca_is_lig[i] == 0) {
          cbx[i] = pro_in.m_mol[chunk_snap_idx].pos.x[cb_idx];
          cby[i] = pro_in.m_mol[chunk_snap_idx].pos.y[cb_idx];
          cbz[i] = pro_in.m_mol[chunk_snap_idx].pos.z[cb_idx];
        } else {
          cbx[i] = lig.m_mol[l_snap].pos.x[cb_idx];
          cby[i] = lig.m_mol[l_snap].pos.y[cb_idx];
          cbz[i] = lig.m_mol[l_snap].pos.z[cb_idx];
        }
      } else {
        cbx[i] = cx[i];
        cby[i] = cy[i];
        cbz[i] = cz[i];
      }
    }
  }
  if (store_full_traj_) {
    traj_x.emplace_back(std::move(cx));
    traj_y.emplace_back(std::move(cy));
    traj_z.emplace_back(std::move(cz));
    if (store_cb_traj_) {
      traj_cbx.emplace_back(std::move(cbx));
      traj_cby.emplace_back(std::move(cby));
      traj_cbz.emplace_back(std::move(cbz));
    }
    frame_global_indices_.push_back(static_cast<int64_t>(global_frame_idx));
    traj_times_.push_back(time_ps);
    traj_replicas_.push_back(replica_id);
  } else {
    latest_traj_x_ = std::move(cx);
    latest_traj_y_ = std::move(cy);
    latest_traj_z_ = std::move(cz);
  }
  total_frames++;
}

void DCCMNetwork::calculate_and_export(size_t num_threads) {
  if (total_frames == 0 || ca_indices.empty())
    return;
  size_t num_ca = ca_indices.size();

  // ── STEP 0: Kabsch align (guarded — skip if already aligned) ───────────────
  if (!kabsch_aligned_) {
    Eigen::MatrixXd ref(num_ca, 3);
    for (size_t a = 0; a < num_ca; ++a) {
      ref(a, 0) = traj_x[0][a];
      ref(a, 1) = traj_y[0][a];
      ref(a, 2) = traj_z[0][a];
    }
    ref.rowwise() -= ref.colwise().mean().eval();

    Eigen::MatrixXd mobile(num_ca, 3);
    Eigen::Matrix3d R;
    for (size_t f = 0; f < total_frames; ++f) {
      for (size_t a = 0; a < num_ca; ++a) {
        mobile(a, 0) = traj_x[f][a];
        mobile(a, 1) = traj_y[f][a];
        mobile(a, 2) = traj_z[f][a];
      }
      Eigen::RowVector3d ca_centroid = mobile.colwise().mean();
      TitanAlignment::align_frame_inplace(ref, mobile, R);
      for (size_t a = 0; a < num_ca; ++a) {
        traj_x[f][a] = mobile(a, 0);
        traj_y[f][a] = mobile(a, 1);
        traj_z[f][a] = mobile(a, 2);
      }
      if (!traj_cbx.empty()) {
        for (size_t a = 0; a < num_ca; ++a) {
          Eigen::RowVector3d cb(traj_cbx[f][a], traj_cby[f][a],
                                traj_cbz[f][a]);
          cb -= ca_centroid;
          cb = (R * cb.transpose()).transpose();
          traj_cbx[f][a] = cb(0);
          traj_cby[f][a] = cb(1);
          traj_cbz[f][a] = cb(2);
        }
      }
    }
    kabsch_aligned_ = true;
  }

  // 1. Pivot Data [atom][frame] for cache locality using aligned rows.
  std::vector<double> mean_x(num_ca, 0.0), mean_y(num_ca, 0.0),
      mean_z(num_ca, 0.0);
  for (size_t f = 0; f < total_frames; ++f) {
    for (size_t a = 0; a < num_ca; ++a) {
      mean_x[a] += traj_x[f][a];
      mean_y[a] += traj_y[f][a];
      mean_z[a] += traj_z[f][a];
    }
  }
  const double inv_total_frames = 1.0 / static_cast<double>(total_frames);
  for (size_t a = 0; a < num_ca; ++a) {
    mean_x[a] *= inv_total_frames;
    mean_y[a] *= inv_total_frames;
    mean_z[a] *= inv_total_frames;
  }

  FlatAlignedMatrix dev_x(num_ca, total_frames, 0.0);
  FlatAlignedMatrix dev_y(num_ca, total_frames, 0.0);
  FlatAlignedMatrix dev_z(num_ca, total_frames, 0.0);

  const size_t prep_worker_count = std::max<size_t>(1, num_threads);
  {
    std::vector<std::jthread> prep_workers;
    prep_workers.reserve(prep_worker_count);
    const size_t atoms_per_worker =
        (num_ca + prep_worker_count - 1) / prep_worker_count;
    for (size_t t = 0; t < prep_worker_count; ++t) {
      const size_t a_begin = t * atoms_per_worker;
      const size_t a_end = std::min(a_begin + atoms_per_worker, num_ca);
      if (a_begin >= a_end)
        continue;
      prep_workers.emplace_back([&, a_begin, a_end]() {
        for (size_t a = a_begin; a < a_end; ++a) {
          const double avg_x = mean_x[a];
          const double avg_y = mean_y[a];
          const double avg_z = mean_z[a];
          double *dx = dev_x.row_data(a);
          double *dy = dev_y.row_data(a);
          double *dz = dev_z.row_data(a);
          for (size_t f = 0; f < total_frames; ++f) {
            dx[f] = traj_x[f][a] - avg_x;
            dy[f] = traj_y[f][a] - avg_y;
            dz[f] = traj_z[f][a] - avg_z;
          }
        }
      });
    }
  }

  // 2. Tính Ma trận Hiệp phương sai (Multithreaded + AVX2)
  //    FIX: pair-based distribution for balanced workload,
  //         no mirror-write in workers to avoid data race.
  FlatAlignedMatrix cov(num_ca, num_ca, 0.0);

  // Build upper-triangle pair list (including diagonal)
  std::vector<std::pair<size_t, size_t>> pairs;
  pairs.reserve(num_ca * (num_ca + 1) / 2);
  for (size_t i = 0; i < num_ca; ++i)
    for (size_t j = i; j < num_ca; ++j)
      pairs.emplace_back(i, j);

  // Worker: each thread processes a contiguous slice of pairs
  auto calc_cov_worker = [&](size_t p_start, size_t p_end) {
    for (size_t p = p_start; p < p_end; ++p) {
      auto [i, j] = pairs[p];
      const double *dx_i = dev_x.row_data(i);
      const double *dx_j = dev_x.row_data(j);
      const double *dy_i = dev_y.row_data(i);
      const double *dy_j = dev_y.row_data(j);
      const double *dz_i = dev_z.row_data(i);
      const double *dz_j = dev_z.row_data(j);
      double dot_x = fast_dot_product_aligned(dx_i, dx_j, total_frames);
      double dot_y = fast_dot_product_aligned(dy_i, dy_j, total_frames);
      double dot_z = fast_dot_product_aligned(dz_i, dz_j, total_frames);

      cov(i, j) = (dot_x + dot_y + dot_z) / static_cast<double>(total_frames);
    }
  };

  // Distribute pairs equally among threads
  const size_t cov_worker_count = std::max<size_t>(1, num_threads);
  {
    std::vector<std::jthread> workers;
    workers.reserve(cov_worker_count);
    size_t total_pairs = pairs.size();
    size_t chunk = (total_pairs + cov_worker_count - 1) / cov_worker_count;
    for (size_t t = 0; t < cov_worker_count; ++t) {
      size_t start = t * chunk;
      size_t end = std::min(start + chunk, total_pairs);
      if (start < end) {
        workers.emplace_back(calc_cov_worker, start, end);
      }
    }
  }

  // Single-threaded symmetry fill (safe, no race)
  for (size_t i = 0; i < num_ca; ++i)
    for (size_t j = i + 1; j < num_ca; ++j)
      cov(j, i) = cov(i, j);

  // 3. Tính Correlation và Xuất File
  std::ofstream edges_out("network_edges.csv");
  edges_out << "Source,Target,Correlation\n";
  for (size_t i = 0; i < num_ca; ++i) {
    for (size_t j = i + 1; j < num_ca; ++j) {
      double denom = std::sqrt(cov(i, i) * cov(j, j));
      if (denom > 1e-12) {
        double corr = cov(i, j) / denom;
        if (std::abs(corr) > 0.4) {
          edges_out << ca_resname[i] << ca_resnrs[i] << "_" << ca_chain[i]
                    << "," << ca_resname[j] << ca_resnrs[j] << "_"
                    << ca_chain[j] << "," << corr << "\n";
        }
      }
    }
  }
  edges_out.close();

  std::ofstream nodes_out("network_nodes.csv");
  nodes_out << "Id,Label,Type,AtomName\n";
  for (size_t i = 0; i < num_ca; ++i) {
    nodes_out << ca_resname[i] << ca_resnrs[i] << "_" << ca_chain[i] << ","
              << ca_resnrs[i] << "," << (ca_is_lig[i] ? "Ligand" : "Receptor")
              << "," << ca_resname[i] << "\n";
  }
  nodes_out.close();
}

// =============================================================================
// Information-Theory Allosteric Analysis: GC + Transfer Entropy
// =============================================================================
void DCCMNetwork::calculate_allostery_information(size_t num_threads,
                                                  int tau_lag, int num_bins) {
  AllosteryConfig cfg;
  cfg.tau_lag = tau_lag;
  cfg.num_bins = num_bins;
  calculate_allostery_information(num_threads, cfg);
}

void DCCMNetwork::calculate_allostery_information(
    size_t num_threads, const AllosteryConfig &cfg) {
  last_allostery_candidate_pairs_ = 0;
  last_allostery_dense_pairs_ = 0;
  const int resolved_bin_count =
      resolve_allostery_bin_count(total_frames, cfg.num_bins);
  const auto dispatch_node_labels =
      build_allostery_node_labels(ca_resname, ca_resnrs, ca_chain);
  last_allostery_candidate_pairs_ =
      calculate_allostery_information_flat_impl<std::uint32_t>(
          traj_x, traj_y, traj_z, dispatch_node_labels, total_frames,
          num_threads, cfg, resolved_bin_count, &last_allostery_dense_pairs_);
}

void DCCMNetwork::calculate_allostery_information_optimized(
    size_t num_threads, int tau_lag, int num_bins) {
  AllosteryConfig cfg;
  cfg.tau_lag = tau_lag;
  cfg.num_bins = num_bins;
  calculate_allostery_information_optimized(num_threads, cfg);
}

void DCCMNetwork::calculate_allostery_information_optimized(
    size_t num_threads, const AllosteryConfig &cfg) {
  last_allostery_candidate_pairs_ = 0;
  last_allostery_dense_pairs_ = 0;
  const int resolved_bin_count =
      resolve_allostery_bin_count(total_frames, cfg.num_bins);
  const auto dispatch_node_labels =
      build_allostery_node_labels(ca_resname, ca_resnrs, ca_chain);
  const size_t resolved_bin_count_sz =
      static_cast<size_t>(resolved_bin_count);
  if (resolved_bin_count_sz <=
      static_cast<size_t>(std::numeric_limits<std::uint8_t>::max())) {
    last_allostery_candidate_pairs_ =
        calculate_allostery_information_flat_impl<std::uint8_t>(
            traj_x, traj_y, traj_z, dispatch_node_labels, total_frames,
            num_threads, cfg, resolved_bin_count,
            &last_allostery_dense_pairs_);
  } else if (resolved_bin_count_sz <=
             static_cast<size_t>(std::numeric_limits<std::uint16_t>::max())) {
    last_allostery_candidate_pairs_ =
        calculate_allostery_information_flat_impl<std::uint16_t>(
            traj_x, traj_y, traj_z, dispatch_node_labels, total_frames,
            num_threads, cfg, resolved_bin_count,
            &last_allostery_dense_pairs_);
  } else {
    last_allostery_candidate_pairs_ =
        calculate_allostery_information_flat_impl<std::uint32_t>(
            traj_x, traj_y, traj_z, dispatch_node_labels, total_frames,
            num_threads, cfg, resolved_bin_count,
            &last_allostery_dense_pairs_);
  }
}
// =============================================================================
// LSP-MD: Local Spatial Pattern (Kornev, JCTC 2025, 21, 8699)
// =============================================================================

// ── Brandes weighted betweenness centrality (Dijkstra variant) ────────────────
namespace {
void brandes_weighted_bc(
    const titan_sparse::WeightedCSRGraph<double> &graph,
    std::vector<double> &bc_out,
    size_t num_threads = 1) {
  const size_t n = static_cast<size_t>(graph.node_count());
  const size_t worker_count = std::max<size_t>(1, num_threads);
  bc_out.assign(n, 0.0);

  // Per-thread partial sums (avoids false sharing / contention)
  FlatAlignedMatrix partial(worker_count, n, 0.0);

  auto worker = [&](size_t t, size_t s_start, size_t s_end) {
    double *bc_local = partial.row_data(t);
    // Allocate once per thread, reuse across source iterations
    std::vector<double> dist(n), sigma(n), delta(n);
    std::vector<std::vector<size_t>> pred(n);
    for (size_t node = 0; node < n; ++node)
      pred[node].reserve(
          static_cast<size_t>(std::max(1, graph.degree(static_cast<int>(node)))));
    std::vector<size_t> order;
    order.reserve(n);

    for (size_t s = s_start; s < s_end; ++s) {
      std::fill(dist.begin(), dist.end(), 1e30);
      std::fill(sigma.begin(), sigma.end(), 0.0);
      std::fill(delta.begin(), delta.end(), 0.0);
      for (auto &p : pred) p.clear(); // clear keeps capacity
      order.clear();

      using PQEntry = std::pair<double, size_t>;
      std::priority_queue<PQEntry, std::vector<PQEntry>,
                          std::greater<PQEntry>> pq;

      dist[s] = 0.0;
      sigma[s] = 1.0;
      pq.push({0.0, s});

      while (!pq.empty()) {
        auto [d_u, u] = pq.top();
        pq.pop();
        if (d_u > dist[u])
          continue;
        order.push_back(u);
        for (int edge_idx = graph.row_begin(static_cast<int>(u));
             edge_idx < graph.row_end(static_cast<int>(u)); ++edge_idx) {
          const size_t v =
              static_cast<size_t>(graph.col_index(edge_idx));
          const double w = graph.value(edge_idx);
          double alt = dist[u] + w;
          if (alt < dist[v] - 1e-12) {
            dist[v] = alt;
            sigma[v] = sigma[u];
            pred[v].clear();
            pred[v].push_back(u);
            pq.push({alt, v});
          } else if (std::abs(alt - dist[v]) < 1e-12) {
            sigma[v] += sigma[u];
            pred[v].push_back(u);
          }
        }
      }

      // Back-propagation of dependencies
      for (auto it = order.rbegin(); it != order.rend(); ++it) {
        size_t w = *it;
        for (size_t p : pred[w]) {
          double coeff = (sigma[p] / sigma[w]) * (1.0 + delta[w]);
          delta[p] += coeff;
        }
        if (w != s) {
          bc_local[w] += delta[w];
        }
      }
    }
  };

  // Distribute source nodes across threads
  {
    std::vector<std::jthread> threads;
    size_t chunk = (n + worker_count - 1) / worker_count;
    for (size_t t = 0; t < worker_count; ++t) {
      size_t s = t * chunk;
      size_t e = std::min(s + chunk, n);
      if (s < e)
        threads.emplace_back(worker, t, s, e);
    }
  }

  // Reduce partial sums
  for (size_t t = 0; t < worker_count; ++t) {
    const double *partial_row = partial.row_data(t);
    for (size_t i = 0; i < n; ++i)
      bc_out[i] += partial_row[i];
  }

  // Normalize for undirected graph
  if (n > 2) {
    double norm = static_cast<double>((n - 1) * (n - 2)) / 2.0;
    for (auto &v : bc_out)
      v /= norm;
  }
}
} // anonymous namespace

// ── AVX2: compute 4 distances in parallel ────────────────────────────────────
namespace {
#ifdef __AVX__
inline void compute_4dist_avx(
    double cax_i, double cay_i, double caz_i,
    double cbx_i, double cby_i, double cbz_i,
    double cax_j, double cay_j, double caz_j,
    double cbx_j, double cby_j, double cbz_j,
    double &d1, double &d2, double &d3, double &d4) {
  // _mm256_set_pd order: arg3→lane3, arg2→lane2, arg1→lane1, arg0→lane0
  // buf[0]=d1(CA-CA), buf[1]=d2(CA-CB), buf[2]=d3(CB-CA), buf[3]=d4(CB-CB)
  __m256d dx = _mm256_set_pd(cbx_i - cbx_j, cbx_i - cax_j,
                              cax_i - cbx_j, cax_i - cax_j);
  __m256d dy = _mm256_set_pd(cby_i - cby_j, cby_i - cay_j,
                              cay_i - cby_j, cay_i - cay_j);
  __m256d dz = _mm256_set_pd(cbz_i - cbz_j, cbz_i - caz_j,
                              caz_i - cbz_j, caz_i - caz_j);
  __m256d r2 = simd_compat::fmadd_pd(
      dz, dz, simd_compat::fmadd_pd(dy, dy, _mm256_mul_pd(dx, dx)));
  __m256d r = _mm256_sqrt_pd(r2);
  double buf[4];
  _mm256_storeu_pd(buf, r); // unaligned store — safe on stack
  d1 = buf[0]; d2 = buf[1]; d3 = buf[2]; d4 = buf[3];
}
#endif

[[maybe_unused]] inline void compute_4dist_scalar(
    double cax_i, double cay_i, double caz_i,
    double cbx_i, double cby_i, double cbz_i,
    double cax_j, double cay_j, double caz_j,
    double cbx_j, double cby_j, double cbz_j,
    double &d1, double &d2, double &d3, double &d4) {
  double dx, dy, dz;
  dx = cax_i - cax_j; dy = cay_i - cay_j; dz = caz_i - caz_j;
  d1 = std::sqrt(dx*dx + dy*dy + dz*dz);
  dx = cax_i - cbx_j; dy = cay_i - cby_j; dz = caz_i - cbz_j;
  d2 = std::sqrt(dx*dx + dy*dy + dz*dz);
  dx = cbx_i - cax_j; dy = cby_i - cay_j; dz = cbz_i - caz_j;
  d3 = std::sqrt(dx*dx + dy*dy + dz*dz);
  dx = cbx_i - cbx_j; dy = cby_i - cby_j; dz = cbz_i - cbz_j;
  d4 = std::sqrt(dx*dx + dy*dy + dz*dz);
}
} // anonymous namespace

void DCCMNetwork::calculate_lspmd(size_t num_threads,
                                   const LSPMDConfig &cfg,
                                   LSPMDResult &result) {
  if (total_frames == 0 || ca_indices.empty())
    return;

  if (total_frames < 10) {
    std::cerr << "[LSP-MD] WARNING: Too few frames (" << total_frames
              << " < 10), results may be unreliable.\n";
  }

  size_t num_ca = ca_indices.size();
  double cutoff = cfg.cutoff_ang;
  const auto has_full_cb_cache = [&]() {
    if (traj_cbx.size() != total_frames || traj_cby.size() != total_frames ||
        traj_cbz.size() != total_frames) {
      return false;
    }
    for (size_t f = 0; f < total_frames; ++f) {
      if (traj_cbx[f].size() < num_ca || traj_cby[f].size() < num_ca ||
          traj_cbz[f].size() < num_ca) {
        return false;
      }
    }
    return true;
  };
  const bool use_cb_cache = has_full_cb_cache();
  if (!use_cb_cache) {
    TitanUI::print_warning(
        "LSP-MD: Cbeta/N trajectory cache is unavailable or incomplete. "
        "Falling back to Calpha coordinates for side-chain anchors.");
  }

  // ── STEP 0: Kabsch alignment (guarded — skip if already aligned) ───────────
  if (!kabsch_aligned_) {
    Eigen::MatrixXd ref(num_ca, 3);
    for (size_t a = 0; a < num_ca; ++a) {
      ref(a, 0) = traj_x[0][a];
      ref(a, 1) = traj_y[0][a];
      ref(a, 2) = traj_z[0][a];
    }
    ref.rowwise() -= ref.colwise().mean().eval();

    Eigen::MatrixXd mobile(num_ca, 3);
    Eigen::Matrix3d R;
    for (size_t f = 0; f < total_frames; ++f) {
      for (size_t a = 0; a < num_ca; ++a) {
        mobile(a, 0) = traj_x[f][a];
        mobile(a, 1) = traj_y[f][a];
        mobile(a, 2) = traj_z[f][a];
      }
      Eigen::RowVector3d ca_centroid = mobile.colwise().mean();
      TitanAlignment::align_frame_inplace(ref, mobile, R);
      for (size_t a = 0; a < num_ca; ++a) {
        traj_x[f][a] = mobile(a, 0);
        traj_y[f][a] = mobile(a, 1);
        traj_z[f][a] = mobile(a, 2);
      }
      if (!traj_cbx.empty()) {
        for (size_t a = 0; a < num_ca; ++a) {
          Eigen::RowVector3d cb(traj_cbx[f][a], traj_cby[f][a],
                                traj_cbz[f][a]);
          cb -= ca_centroid;
          cb = (R * cb.transpose()).transpose();
          traj_cbx[f][a] = cb(0);
          traj_cby[f][a] = cb(1);
          traj_cbz[f][a] = cb(2);
        }
      }
    }
    kabsch_aligned_ = true;
  }

  // ── Pre-transpose to atom-major layout [atom][frame] for cache locality ────
  size_t F = total_frames;
  FlatAlignedMatrix ax(num_ca, F, 0.0);
  FlatAlignedMatrix ay(num_ca, F, 0.0);
  FlatAlignedMatrix az(num_ca, F, 0.0);
  FlatAlignedMatrix bx(num_ca, F, 0.0);
  FlatAlignedMatrix by(num_ca, F, 0.0);
  FlatAlignedMatrix bz(num_ca, F, 0.0);

  for (size_t f = 0; f < F; ++f)
    for (size_t a = 0; a < num_ca; ++a) {
      const double x = traj_x[f][a];
      const double y = traj_y[f][a];
      const double z = traj_z[f][a];
      ax(a, f) = x;
      ay(a, f) = y;
      az(a, f) = z;
      if (use_cb_cache) {
        bx(a, f) = traj_cbx[f][a];
        by(a, f) = traj_cby[f][a];
        bz(a, f) = traj_cbz[f][a];
      } else {
        bx(a, f) = ax(a, f);
        by(a, f) = ay(a, f);
        bz(a, f) = az(a, f);
      }
    }

  // ── STEP 1: Build upper-triangle pair list and compute W_ij ────────────────
  std::vector<std::pair<size_t, size_t>> pairs = build_lspmd_candidate_pairs(
      traj_x, traj_y, traj_z, num_ca, total_frames, cutoff,
      cfg.use_avg_cutoff);
  if (pairs.empty()) {
    result.edges.clear();
    result.dc.assign(num_ca, 0.0);
    result.bc.assign(num_ca, 0.0);
    return;
  }

  const size_t n_pairs = pairs.size();
  std::vector<double> W_vals(n_pairs,
      std::numeric_limits<double>::quiet_NaN());

  auto worker = [&](size_t p_start, size_t p_end) {
    AlignedDoubleVector d1v(F, 0.0), d2v(F, 0.0), d3v(F, 0.0), d4v(F, 0.0);
    for (size_t p = p_start; p < p_end; ++p) {
      size_t i = pairs[p].first;
      size_t j = pairs[p].second;

      // Pointers into transposed arrays — sequential access in inner loop
      const double *axi = ax.row_data(i), *ayi = ay.row_data(i),
                   *azi = az.row_data(i);
      const double *bxi = bx.row_data(i), *byi = by.row_data(i),
                   *bzi = bz.row_data(i);
      const double *axj = ax.row_data(j), *ayj = ay.row_data(j),
                   *azj = az.row_data(j);
      const double *bxj = bx.row_data(j), *byj = by.row_data(j),
                   *bzj = bz.row_data(j);

      // ── Single pass: compute distances + cache for deviation ──
      double sum_d1 = 0.0, sum_d2 = 0.0, sum_d3 = 0.0, sum_d4 = 0.0;
      for (size_t f = 0; f < F; ++f) {
#ifdef __AVX__
        compute_4dist_avx(axi[f], ayi[f], azi[f], bxi[f], byi[f], bzi[f],
                          axj[f], ayj[f], azj[f], bxj[f], byj[f], bzj[f],
                          d1v[f], d2v[f], d3v[f], d4v[f]);
#else
        compute_4dist_scalar(axi[f], ayi[f], azi[f], bxi[f], byi[f], bzi[f],
                             axj[f], ayj[f], azj[f], bxj[f], byj[f], bzj[f],
                             d1v[f], d2v[f], d3v[f], d4v[f]);
#endif
        sum_d1 += d1v[f]; sum_d2 += d2v[f];
        sum_d3 += d3v[f]; sum_d4 += d4v[f];
      }

      double inv_f = 1.0 / static_cast<double>(F);
      double avg_d1 = sum_d1 * inv_f;
      double avg_d2 = sum_d2 * inv_f;
      double avg_d3 = sum_d3 * inv_f;
      double avg_d4 = sum_d4 * inv_f;

      // ── Cutoff check ──
      if (cfg.use_avg_cutoff && avg_d1 > cutoff)
        continue;

      // ── Reuse cached distances for deviation ──
      double sum_Z2 = 0.0;
      for (size_t f = 0; f < F; ++f) {
        double dd1 = d1v[f] - avg_d1, dd2 = d2v[f] - avg_d2;
        double dd3 = d3v[f] - avg_d3, dd4 = d4v[f] - avg_d4;
        sum_Z2 += dd1*dd1 + dd2*dd2 + dd3*dd3 + dd4*dd4;
      }

      double Z_avg = std::sqrt(sum_Z2 * inv_f);
      W_vals[p] = std::exp(-Z_avg);
    }
  };

  // ── Launch multithreaded pair computation ──
  {
    std::vector<std::jthread> workers;
    const size_t worker_count = std::max<size_t>(1, num_threads);
    size_t chunk = (n_pairs + worker_count - 1) / worker_count;
    for (size_t t = 0; t < worker_count; ++t) {
      size_t start = t * chunk;
      size_t end = std::min(start + chunk, n_pairs);
      if (start < end)
        workers.emplace_back(worker, start, end);
    }
  }

  // ── Collect valid edges (NaN = skipped) ──
  result.edges.clear();
  for (size_t p = 0; p < n_pairs; ++p) {
    if (!std::isnan(W_vals[p])) {
      result.edges.emplace_back(static_cast<int>(pairs[p].first),
                                static_cast<int>(pairs[p].second),
                                W_vals[p]);
    }
  }

  // ── STEP 2: Degree Centrality ──────────────────────────────────────────────
  std::sort(result.edges.begin(), result.edges.end(),
            [](const auto &lhs, const auto &rhs) {
              if (std::get<0>(lhs) != std::get<0>(rhs))
                return std::get<0>(lhs) < std::get<0>(rhs);
              if (std::get<1>(lhs) != std::get<1>(rhs))
                return std::get<1>(lhs) < std::get<1>(rhs);
              return std::get<2>(lhs) > std::get<2>(rhs);
            });

  result.dc.assign(num_ca, 0.0);
  for (auto &[i, j, w] : result.edges) {
    result.dc[i] += w;
    result.dc[j] += w;
  }

  // ── STEP 3: Betweenness Centrality (Brandes + Dijkstra) ────────────────────
  // Build adjacency list with D_ij = -ln(W_ij) as distance
  std::vector<std::tuple<int, int, double>> weighted_edges;
  weighted_edges.reserve(result.edges.size());
  for (auto &[i, j, w] : result.edges) {
    const double d_weight = (w > 1e-15) ? -std::log(w) : 35.0;
    weighted_edges.emplace_back(i, j, d_weight);
  }
  const auto weighted_graph =
      titan_sparse::WeightedCSRGraph<double>::from_undirected_edges(
          static_cast<int>(num_ca), weighted_edges);

  brandes_weighted_bc(weighted_graph, result.bc, num_threads);
}

void DCCMNetwork::export_lspmd(const LSPMDResult &result,
                                double weight_threshold) {
  const auto &resname = ca_resname;
  const auto &resnrs = ca_resnrs;
  const auto &chains = ca_chain;
  size_t N = resname.size();

  // --- 1. Edges ---
  {
    std::ofstream ofs("lspmd_edges.csv");
    ofs << "Source,Target,Z_avg,Weight,D_ij\n";
    ofs << std::fixed << std::setprecision(6);
    for (const auto &[i, j, W] : result.edges) {
      if (W < weight_threshold)
        continue;
      double Z = -std::log(W > 1e-15 ? W : 1e-15);
      double D = Z; // D_ij = -ln(W_ij) = <Z>
      ofs << resname[i] << resnrs[i] << "_" << chains[i] << ","
          << resname[j] << resnrs[j] << "_" << chains[j] << "," << Z << ","
          << W << "," << D << "\n";
    }
  }

  // --- 2. Degree Centrality ---
  {
    double dc_max = *std::max_element(result.dc.begin(), result.dc.end());
    if (dc_max < 1e-15)
      dc_max = 1.0;
    std::ofstream ofs("lspmd_degree_centrality.csv");
    ofs << "Residue,ResNr,Chain,IsLig,DC,DC_norm\n";
    ofs << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < N; ++i) {
      ofs << resname[i] << "," << resnrs[i] << "," << chains[i] << ","
          << ca_is_lig[i] << "," << result.dc[i] << ","
          << result.dc[i] / dc_max << "\n";
    }
  }

  // --- 3. Betweenness Centrality ---
  {
    double bc_max = *std::max_element(result.bc.begin(), result.bc.end());
    if (bc_max < 1e-15)
      bc_max = 1.0;
    std::ofstream ofs("lspmd_betweenness_centrality.csv");
    ofs << "Residue,ResNr,Chain,IsLig,BC,BC_norm\n";
    ofs << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < N; ++i) {
      ofs << resname[i] << "," << resnrs[i] << "," << chains[i] << ","
          << ca_is_lig[i] << "," << result.bc[i] << ","
          << result.bc[i] / bc_max << "\n";
    }
  }

  TitanUI::print_info("[LSP-MD] Exported: lspmd_edges.csv, "
                      "lspmd_degree_centrality.csv, "
                      "lspmd_betweenness_centrality.csv");
}
