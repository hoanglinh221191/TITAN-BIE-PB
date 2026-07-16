#include "component_graph_math.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>

namespace titan_pbctopo {

void finalize_normalized_pair_shift_score(PbctopoPairShiftScore &score,
                                          std::size_t lhs_atom_count,
                                          std::size_t rhs_atom_count,
                                          double contact_cutoff) {
  (void)rhs_atom_count;
  const double lhs_count = static_cast<double>(std::max<std::size_t>(1, lhs_atom_count));
  score.contact_fraction = static_cast<double>(score.contact_atoms) / lhs_count;
  score.soft_penalty_mean =
      score.contact_pairs > 0
          ? score.soft_penalty / static_cast<double>(score.contact_pairs)
          : 0.0;
  const double contact_scale2 = std::max(1.0e-12, contact_cutoff * contact_cutoff);
  const double normalized_bbox =
      std::isfinite(score.bbox_d2) ? score.bbox_d2 / contact_scale2 : 1.0;
  score.quality_score = score.contact_fraction * 1000.0 -
                        score.soft_penalty_mean * 50.0 -
                        normalized_bbox * 5.0 - score.shift_norm2;
}

namespace {

struct Dsu {
  std::vector<std::size_t> parent;
  std::vector<unsigned char> rank;

  explicit Dsu(std::size_t n) : parent(n), rank(n, 0) {
    std::iota(parent.begin(), parent.end(), 0);
  }

  std::size_t find(std::size_t node) {
    if (parent[node] != node)
      parent[node] = find(parent[node]);
    return parent[node];
  }

  bool unite(std::size_t lhs, std::size_t rhs) {
    lhs = find(lhs);
    rhs = find(rhs);
    if (lhs == rhs)
      return false;
    if (rank[lhs] < rank[rhs])
      std::swap(lhs, rhs);
    parent[rhs] = lhs;
    if (rank[lhs] == rank[rhs])
      ++rank[lhs];
    return true;
  }
};

Int3 negate(Int3 value) { return {-value.x, -value.y, -value.z}; }

double residual_norm2(const Int3 &value) {
  return static_cast<double>(value.x * value.x + value.y * value.y +
                             value.z * value.z);
}

} // namespace

PbctopoCycleSyncMetrics compute_fundamental_cycle_residuals(
    std::size_t component_count,
    const std::vector<PbctopoPairShiftCandidateSet> &pair_candidates) {
  PbctopoCycleSyncMetrics metrics;
  if (component_count < 2 || pair_candidates.empty())
    return metrics;

  std::vector<std::size_t> order(pair_candidates.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs,
                                                   std::size_t rhs) {
    const double lhs_weight = pair_candidates[lhs].candidates.empty()
                                  ? 0.0
                                  : pair_candidates[lhs].candidates.front().weight;
    const double rhs_weight = pair_candidates[rhs].candidates.empty()
                                  ? 0.0
                                  : pair_candidates[rhs].candidates.front().weight;
    return lhs_weight > rhs_weight;
  });

  struct TreeEdge {
    std::size_t node = 0;
    Int3 delta;
  };
  std::vector<std::vector<TreeEdge>> tree(component_count);
  std::vector<char> is_tree_edge(pair_candidates.size(), 0);
  Dsu dsu(component_count);
  for (std::size_t edge_idx : order) {
    const auto &set = pair_candidates[edge_idx];
    if (set.a >= component_count || set.b >= component_count ||
        set.candidates.empty())
      continue;
    if (!dsu.unite(set.a, set.b))
      continue;
    const Int3 delta = set.candidates.front().delta;
    tree[set.a].push_back({set.b, delta});
    tree[set.b].push_back({set.a, negate(delta)});
    is_tree_edge[edge_idx] = 1;
  }

  std::vector<Int3> potential(component_count);
  std::vector<char> visited(component_count, 0);
  for (std::size_t root = 0; root < component_count; ++root) {
    if (visited[root])
      continue;
    std::queue<std::size_t> queue;
    queue.push(root);
    visited[root] = 1;
    while (!queue.empty()) {
      const std::size_t node = queue.front();
      queue.pop();
      for (const auto &edge : tree[node]) {
        if (visited[edge.node])
          continue;
        potential[edge.node] = potential[node] + edge.delta;
        visited[edge.node] = 1;
        queue.push(edge.node);
      }
    }
  }

  double residual_sum2 = 0.0;
  for (std::size_t edge_idx = 0; edge_idx < pair_candidates.size(); ++edge_idx) {
    if (is_tree_edge[edge_idx])
      continue;
    const auto &set = pair_candidates[edge_idx];
    if (set.a >= component_count || set.b >= component_count ||
        set.candidates.empty() || dsu.find(set.a) != dsu.find(set.b))
      continue;
    const Int3 tree_delta = potential[set.b] - potential[set.a];
    double best_residual2 = std::numeric_limits<double>::infinity();
    for (const auto &candidate : set.candidates)
      best_residual2 =
          std::min(best_residual2, residual_norm2(tree_delta - candidate.delta));
    if (!std::isfinite(best_residual2))
      continue;
    ++metrics.fundamental_cycle_count;
    residual_sum2 += best_residual2;
    metrics.fundamental_cycle_residual_max =
        std::max(metrics.fundamental_cycle_residual_max,
                 std::sqrt(best_residual2));
  }
  if (metrics.fundamental_cycle_count > 0) {
    metrics.fundamental_cycle_residual_rms =
        std::sqrt(residual_sum2 /
                  static_cast<double>(metrics.fundamental_cycle_count));
  }
  return metrics;
}

} // namespace titan_pbctopo
