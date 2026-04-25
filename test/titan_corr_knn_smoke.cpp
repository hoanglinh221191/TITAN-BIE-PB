#include "titan_corr_knn_core.h"
#include "titan_corr_knn_shared.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using titan_corr_knn::PairPruningConfig;
using titan_corr_knn::TrajectoryRows;

[[noreturn]] void fail(const char *message) {
  std::cerr << "[titan_corr_knn_smoke] " << message << "\n";
  std::exit(1);
}

TrajectoryRows make_rows(std::size_t frame_count, std::size_t node_count) {
  TrajectoryRows rows(frame_count);
  for (std::size_t frame = 0; frame < frame_count; ++frame) {
    rows[frame].assign(node_count, 0.0);
  }
  return rows;
}

bool has_pair(const std::vector<std::pair<std::size_t, std::size_t>> &pairs,
              std::size_t lhs, std::size_t rhs) {
  for (const auto &pair : pairs) {
    if (pair.first == lhs && pair.second == rhs)
      return true;
  }
  return false;
}

double next_uniform(std::uint64_t &state) {
  state = state * 6364136223846793005ULL + 1442695040888963407ULL;
  const std::uint64_t mantissa = (state >> 11) & ((1ULL << 53) - 1ULL);
  const double unit =
      static_cast<double>(mantissa) / static_cast<double>(1ULL << 53);
  return 2.0 * unit - 1.0;
}

void fill_independent_and_correlated_series(TrajectoryRows &traj_x,
                                            TrajectoryRows &traj_y,
                                            TrajectoryRows &traj_z) {
  std::uint64_t rng_a = 0x123456789abcdef0ULL;
  std::uint64_t rng_b = 0x0fedcba987654321ULL;
  std::uint64_t rng_c = 0x13579bdf2468ace0ULL;
  for (std::size_t frame = 0; frame < traj_x.size(); ++frame) {
    const double ax = next_uniform(rng_a);
    const double ay = next_uniform(rng_a);
    const double az = next_uniform(rng_a);

    const double bx = next_uniform(rng_b);
    const double by = next_uniform(rng_b);
    const double bz = next_uniform(rng_b);

    const double nx = 0.08 * next_uniform(rng_c);
    const double ny = 0.08 * next_uniform(rng_c);
    const double nz = 0.08 * next_uniform(rng_c);

    traj_x[frame][0] = ax;
    traj_y[frame][0] = ay;
    traj_z[frame][0] = az;

    traj_x[frame][1] = bx;
    traj_y[frame][1] = by;
    traj_z[frame][1] = bz;

    traj_x[frame][2] = ax + nx;
    traj_y[frame][2] = ay + ny;
    traj_z[frame][2] = az + nz;

    traj_x[frame][3] = 0.0;
    traj_y[frame][3] = 0.0;
    traj_z[frame][3] = 0.0;
  }
}

void fill_te_driver_series(TrajectoryRows &traj_x, TrajectoryRows &traj_y,
                           TrajectoryRows &traj_z) {
  std::uint64_t rng_src = 0x4242424242424242ULL;
  std::uint64_t rng_tgt = 0xabcdef0123456789ULL;
  double prev_src_x = 0.0;
  double prev_src_y = 0.0;
  double prev_src_z = 0.0;
  double prev_tgt_x = 0.0;
  double prev_tgt_y = 0.0;
  double prev_tgt_z = 0.0;

  for (std::size_t frame = 0; frame < traj_x.size(); ++frame) {
    const double src_x = 0.55 * prev_src_x + 0.45 * next_uniform(rng_src);
    const double src_y = 0.55 * prev_src_y + 0.45 * next_uniform(rng_src);
    const double src_z = 0.55 * prev_src_z + 0.45 * next_uniform(rng_src);

    const double tgt_x =
        0.80 * prev_src_x + 0.15 * prev_tgt_x + 0.05 * next_uniform(rng_tgt);
    const double tgt_y =
        0.80 * prev_src_y + 0.15 * prev_tgt_y + 0.05 * next_uniform(rng_tgt);
    const double tgt_z =
        0.80 * prev_src_z + 0.15 * prev_tgt_z + 0.05 * next_uniform(rng_tgt);

    traj_x[frame][0] = src_x;
    traj_y[frame][0] = src_y;
    traj_z[frame][0] = src_z;
    traj_x[frame][1] = tgt_x;
    traj_y[frame][1] = tgt_y;
    traj_z[frame][1] = tgt_z;

    prev_src_x = src_x;
    prev_src_y = src_y;
    prev_src_z = src_z;
    prev_tgt_x = tgt_x;
    prev_tgt_y = tgt_y;
    prev_tgt_z = tgt_z;
  }
}

} // namespace

int main() {
  {
    auto traj_x = make_rows(2, 3);
    auto traj_y = make_rows(2, 3);
    auto traj_z = make_rows(2, 3);
    auto anchor_x = make_rows(2, 3);
    auto anchor_y = make_rows(2, 3);
    auto anchor_z = make_rows(2, 3);

    traj_x[0][0] = 0.0;
    traj_x[0][1] = 0.8;
    traj_x[0][2] = 4.0;
    traj_x[1][0] = 4.0;
    traj_x[1][1] = 0.9;
    traj_x[1][2] = 1.7;
    anchor_x = traj_x;

    PairPruningConfig cfg;
    cfg.pair_cutoff_ang = 1.1;
    cfg.use_avg_cutoff = true;
    auto pairs = titan_corr_knn::build_candidate_pairs(
        traj_x, traj_y, traj_z, anchor_x, anchor_y, anchor_z, 3, 2, cfg,
        nullptr);
    if (pairs.size() != 2 || !has_pair(pairs, 0, 1) || !has_pair(pairs, 1, 2)) {
      fail("avg cutoff pruning should keep union across frames");
    }

    cfg.use_avg_cutoff = false;
    pairs = titan_corr_knn::build_candidate_pairs(
        traj_x, traj_y, traj_z, anchor_x, anchor_y, anchor_z, 3, 2, cfg,
        nullptr);
    if (pairs.size() != 1 || !has_pair(pairs, 0, 1)) {
      fail("first-frame cutoff pruning should only keep frame-0 pair");
    }
  }

  {
    auto traj_x = make_rows(2, 3);
    auto traj_y = make_rows(2, 3);
    auto traj_z = make_rows(2, 3);
    auto anchor_x = make_rows(2, 3);
    auto anchor_y = make_rows(2, 3);
    auto anchor_z = make_rows(2, 3);

    traj_x[0][0] = 0.0;
    traj_x[0][1] = 0.7;
    traj_x[0][2] = 1.4;
    traj_x[1][0] = 0.0;
    traj_x[1][1] = 0.6;
    traj_x[1][2] = 1.3;
    anchor_x = traj_x;

    PairPruningConfig cfg;
    cfg.pair_topk = 1;
    cfg.use_avg_cutoff = true;
    auto frame_pairs = titan_corr_knn::build_candidate_pairs(
        traj_x, traj_y, traj_z, anchor_x, anchor_y, anchor_z, 3, 2, cfg,
        nullptr);
    if (frame_pairs.size() != 2 || !has_pair(frame_pairs, 0, 1) ||
        !has_pair(frame_pairs, 1, 2)) {
      fail("frame-voted top-k pruning mismatch");
    }

    anchor_x[0][0] = 0.0;
    anchor_x[0][1] = 4.0;
    anchor_x[0][2] = 0.2;
    anchor_x[1][0] = 0.0;
    anchor_x[1][1] = 4.0;
    anchor_x[1][2] = 0.2;
    cfg.pair_use_contact_prior = true;
    cfg.pair_contact_cutoff_ang = 0.5;
    auto contact_pairs = titan_corr_knn::build_candidate_pairs(
        traj_x, traj_y, traj_z, anchor_x, anchor_y, anchor_z, 3, 2, cfg,
        nullptr);
    if (contact_pairs.size() != 1 || !has_pair(contact_pairs, 0, 2)) {
      fail("contact-prior top-k pruning mismatch");
    }
  }

  TitanExecutor executor(2);

  {
    auto traj_x = make_rows(96, 4);
    auto traj_y = make_rows(96, 4);
    auto traj_z = make_rows(96, 4);
    fill_independent_and_correlated_series(traj_x, traj_y, traj_z);

    titan_corr_knn::ExactEstimatorConfig cfg;
    cfg.estimator_k = 2;
    cfg.tau_lag = 1;
    cfg.compute_te = true;
    const auto result = titan_corr_knn::compute_exact_scores(
        traj_x, traj_y, traj_z, {"A", "B", "C", "D"},
        {{0, 1}, {0, 2}, {0, 3}}, 6, cfg, executor);
    if (!result.mi_available)
      fail("exact MI should be available for synthetic smoke data");
    if (result.mi_values.size() != 3 || result.gc_values.size() != 3)
      fail("unexpected MI/GC vector sizes");
    if (result.mi_values[1] <= result.mi_values[0] + 0.02)
      fail("correlated pair should have larger MI than independent pair");
    if (!(result.gc_values[1] > 0.0) || !std::isfinite(result.gc_values[1]))
      fail("correlated GC should be finite and positive");
    if (result.mi_values[2] != 0.0 || result.gc_values[2] != 0.0)
      fail("degenerate node should emit zero MI/GC");
    if (result.degenerate_node_count == 0 || result.degenerate_pair_count == 0)
      fail("degenerate inputs should be reported in the summary result");
  }

  {
    auto traj_x = make_rows(128, 2);
    auto traj_y = make_rows(128, 2);
    auto traj_z = make_rows(128, 2);
    fill_te_driver_series(traj_x, traj_y, traj_z);

    titan_corr_knn::ExactEstimatorConfig cfg;
    cfg.estimator_k = 2;
    cfg.tau_lag = 1;
    cfg.compute_te = true;
    const auto result = titan_corr_knn::compute_exact_scores(
        traj_x, traj_y, traj_z, {"Src", "Tgt"}, {{0, 1}}, 1, cfg, executor);
    if (!result.te_available || result.te_forward_values.size() != 1 ||
        result.te_reverse_values.size() != 1) {
      fail("exact TE should be available for the driver-process smoke data");
    }
    if (result.te_forward_values[0] <= result.te_reverse_values[0] + 0.01) {
      fail("synthetic driver process should satisfy TE(X->Y) > TE(Y->X)");
    }
  }

  std::cout << "titan_corr_knn_smoke PASS\n";
  return 0;
}
