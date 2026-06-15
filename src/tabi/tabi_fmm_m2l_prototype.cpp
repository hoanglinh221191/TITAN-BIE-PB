#include "tabi/tabi_fmm_m2l_prototype.h"
#include "tabi/tabi_types.h"

#include "aligned_allocator.h"
#include "titan_executor.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace titan::tabi::fmm_m2l {
namespace {

using Clock = std::chrono::steady_clock;

struct PointSoA {
  AlignedDoubleVector x;
  AlignedDoubleVector y;
  AlignedDoubleVector z;
  AlignedDoubleVector q;
  std::vector<std::size_t> panel_id;
};

struct Cell {
  std::size_t begin = 0;
  std::size_t end = 0;
  int left = -1;
  int right = -1;
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;
  double radius = 0.0;
  double charge = 0.0;
  double qx = 0.0;
  double qy = 0.0;
  double qz = 0.0;
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;

  [[nodiscard]] bool is_leaf() const noexcept { return left < 0 && right < 0; }
  [[nodiscard]] std::size_t count() const noexcept { return end - begin; }
};

struct TargetLeafWork {
  int target_node = -1;
  std::vector<int> near_sources;
  std::vector<int> far_sources;
};

struct PrototypeEvaluation {
  PrototypeMetrics metrics;
  AlignedDoubleVector reference;
  AlignedDoubleVector approx;
};

struct Tree {
  std::vector<std::size_t> order;
  std::vector<std::size_t> panel_id;
  std::vector<Cell> cells;
  std::vector<int> target_leaf_nodes;
  std::vector<TargetLeafWork> work;
};

double elapsed_seconds(Clock::time_point begin, Clock::time_point end) {
  return std::chrono::duration<double>(end - begin).count();
}

double kernel_value(KernelKind kernel, double kappa, double r2) {
  constexpr double kTiny = 1.0e-28;
  if (!(r2 > kTiny)) {
    return 0.0;
  }
  const double r = std::sqrt(r2);
  if (kernel == KernelKind::Yukawa) {
    return std::exp(-std::max(0.0, kappa) * r) / r;
  }
  return 1.0 / r;
}

PointSoA make_synthetic_points(std::size_t n, unsigned seed) {
  PointSoA points;
  points.x.resize(n);
  points.y.resize(n);
  points.z.resize(n);
  points.q.resize(n);
  points.panel_id.resize(n);

  std::uint64_t state = (static_cast<std::uint64_t>(seed) << 32) ^ 0x9e3779b97f4a7c15ull;
  auto next01 = [&]() {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    const std::uint64_t bits = (state >> 11) | 0x3ff0000000000000ull;
    double value = 0.0;
    static_assert(sizeof(value) == sizeof(bits));
    std::memcpy(&value, &bits, sizeof(value));
    return value - 1.0;
  };

  for (std::size_t i = 0; i < n; ++i) {
    const double cluster = (i < n / 2) ? -1.0 : 1.0;
    const double shell = static_cast<double>((i % 17) + 1) / 18.0;
    const double theta = 6.2831853071795864769 * next01();
    const double phi = std::acos(std::clamp(2.0 * next01() - 1.0, -1.0, 1.0));
    const double jitter = 0.08 * (next01() - 0.5);
    const double r = 0.35 + 0.35 * shell + jitter;
    points.x[i] = cluster * 2.0 + r * std::sin(phi) * std::cos(theta);
    points.y[i] = 0.35 * cluster + r * std::sin(phi) * std::sin(theta);
    points.z[i] = r * std::cos(phi);
    points.q[i] = ((i % 2) == 0 ? 1.0 : -0.7) + 0.05 * (next01() - 0.5);
    points.panel_id[i] = i;
  }
  return points;
}

PointSoA from_public_point_cloud(const PrototypePointCloud &input) {
  const std::size_t n = input.x.size();
  if (input.y.size() != n || input.z.size() != n || input.q.size() != n) {
    throw std::invalid_argument(
        "FMM/M2L prototype point cloud has inconsistent coordinate sizes");
  }
  if (n == 0) {
    throw std::invalid_argument("FMM/M2L prototype requires at least one point");
  }
  PointSoA points;
  points.x.assign(input.x.begin(), input.x.end());
  points.y.assign(input.y.begin(), input.y.end());
  points.z.assign(input.z.begin(), input.z.end());
  points.q.assign(input.q.begin(), input.q.end());
  if (input.panel_id.empty()) {
    points.panel_id.resize(n);
    std::iota(points.panel_id.begin(), points.panel_id.end(), std::size_t{0});
  } else {
    if (input.panel_id.size() != n) {
      throw std::invalid_argument(
          "FMM/M2L prototype point cloud has inconsistent panel ids");
    }
    points.panel_id = input.panel_id;
  }
  return points;
}

void compute_cell_geometry(Tree &tree, const PointSoA &points, int node_index) {
  Cell &cell = tree.cells[static_cast<std::size_t>(node_index)];
  double min_x = std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double min_z = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  double max_z = -std::numeric_limits<double>::infinity();
  double sx = 0.0;
  double sy = 0.0;
  double sz = 0.0;
  double q_sum = 0.0;
  double qx = 0.0;
  double qy = 0.0;
  double qz = 0.0;

  for (std::size_t pos = cell.begin; pos < cell.end; ++pos) {
    const std::size_t idx = tree.order[pos];
    const double x = points.x[idx];
    const double y = points.y[idx];
    const double z = points.z[idx];
    min_x = std::min(min_x, x);
    min_y = std::min(min_y, y);
    min_z = std::min(min_z, z);
    max_x = std::max(max_x, x);
    max_y = std::max(max_y, y);
    max_z = std::max(max_z, z);
    sx += x;
    sy += y;
    sz += z;
    q_sum += points.q[idx];
    qx += points.q[idx] * x;
    qy += points.q[idx] * y;
    qz += points.q[idx] * z;
  }

  const double inv_count = cell.count() > 0 ? 1.0 / static_cast<double>(cell.count()) : 0.0;
  cell.cx = sx * inv_count;
  cell.cy = sy * inv_count;
  cell.cz = sz * inv_count;
  cell.radius = 0.0;
  for (std::size_t pos = cell.begin; pos < cell.end; ++pos) {
    const std::size_t idx = tree.order[pos];
    const double dx = points.x[idx] - cell.cx;
    const double dy = points.y[idx] - cell.cy;
    const double dz = points.z[idx] - cell.cz;
    cell.radius = std::max(cell.radius, std::sqrt(dx * dx + dy * dy + dz * dz));
  }
  cell.charge = q_sum;
  cell.qx = 0.5 * (min_x + max_x);
  cell.qy = 0.5 * (min_y + max_y);
  cell.qz = 0.5 * (min_z + max_z);
  for (std::size_t pos = cell.begin; pos < cell.end; ++pos) {
    const std::size_t idx = tree.order[pos];
    cell.px += points.q[idx] * (points.x[idx] - cell.qx);
    cell.py += points.q[idx] * (points.y[idx] - cell.qy);
    cell.pz += points.q[idx] * (points.z[idx] - cell.qz);
  }
  (void)qx;
  (void)qy;
  (void)qz;
}

int build_node(Tree &tree, const PointSoA &points, std::size_t begin,
               std::size_t end, std::size_t max_leaf) {
  const int node_index = static_cast<int>(tree.cells.size());
  Cell cell;
  cell.begin = begin;
  cell.end = end;
  tree.cells.push_back(cell);
  compute_cell_geometry(tree, points, node_index);

  if (end - begin <= max_leaf) {
    tree.target_leaf_nodes.push_back(node_index);
    return node_index;
  }

  std::array<double, 3> min_v{std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::infinity()};
  std::array<double, 3> max_v{-std::numeric_limits<double>::infinity(),
                              -std::numeric_limits<double>::infinity(),
                              -std::numeric_limits<double>::infinity()};
  for (std::size_t pos = begin; pos < end; ++pos) {
    const std::size_t idx = tree.order[pos];
    min_v[0] = std::min(min_v[0], points.x[idx]);
    min_v[1] = std::min(min_v[1], points.y[idx]);
    min_v[2] = std::min(min_v[2], points.z[idx]);
    max_v[0] = std::max(max_v[0], points.x[idx]);
    max_v[1] = std::max(max_v[1], points.y[idx]);
    max_v[2] = std::max(max_v[2], points.z[idx]);
  }
  const std::array<double, 3> span{max_v[0] - min_v[0], max_v[1] - min_v[1],
                                  max_v[2] - min_v[2]};
  const int axis = (span[1] > span[0] && span[1] >= span[2]) ? 1
                   : (span[2] > span[0] && span[2] > span[1]) ? 2
                                                               : 0;
  const std::size_t mid = begin + (end - begin) / 2;
  auto coord = [&](std::size_t idx) {
    if (axis == 0) {
      return points.x[idx];
    }
    if (axis == 1) {
      return points.y[idx];
    }
    return points.z[idx];
  };
  std::stable_sort(tree.order.begin() + static_cast<std::ptrdiff_t>(begin),
                   tree.order.begin() + static_cast<std::ptrdiff_t>(end),
                   [&](std::size_t lhs, std::size_t rhs) {
                     const double cl = coord(lhs);
                     const double cr = coord(rhs);
                     if (cl != cr) {
                       return cl < cr;
                     }
                     return lhs < rhs;
                   });

  const int left = build_node(tree, points, begin, mid, max_leaf);
  const int right = build_node(tree, points, mid, end, max_leaf);
  tree.cells[static_cast<std::size_t>(node_index)].left = left;
  tree.cells[static_cast<std::size_t>(node_index)].right = right;
  return node_index;
}

Tree build_tree(const PointSoA &points, std::size_t max_leaf) {
  Tree tree;
  tree.order.resize(points.x.size());
  std::iota(tree.order.begin(), tree.order.end(), std::size_t{0});
  tree.panel_id = points.panel_id;
  tree.cells.reserve(points.x.size() / std::max<std::size_t>(1, max_leaf) * 2 + 8);
  build_node(tree, points, 0, points.x.size(), std::max<std::size_t>(1, max_leaf));
  tree.work.resize(tree.cells.size());
  for (int node : tree.target_leaf_nodes) {
    tree.work[static_cast<std::size_t>(node)].target_node = node;
  }
  return tree;
}

bool accept_far(const Tree &target_tree, const Tree &source_tree, int target_node,
                int source_node, double theta, bool exclude_same_index) {
  if (!(theta > 0.0)) {
    return false;
  }
  const Cell &t = target_tree.cells[static_cast<std::size_t>(target_node)];
  const Cell &s = source_tree.cells[static_cast<std::size_t>(source_node)];
  if (exclude_same_index) {
    for (std::size_t tpos = t.begin; tpos < t.end; ++tpos) {
      const std::size_t ti = target_tree.order[tpos];
      for (std::size_t spos = s.begin; spos < s.end; ++spos) {
        const std::size_t si = source_tree.order[spos];
        if (target_tree.panel_id[ti] == source_tree.panel_id[si]) {
          return false;
        }
      }
    }
  }
  const double dx = t.cx - s.cx;
  const double dy = t.cy - s.cy;
  const double dz = t.cz - s.cz;
  const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (!(dist > 1.0e-12)) {
    return false;
  }
  return (t.radius + s.radius) / dist <= theta;
}

void build_work(Tree &target_tree, const Tree &source_tree, int target_node,
                int source_node, double theta, bool exclude_same_index) {
  const Cell &target = target_tree.cells[static_cast<std::size_t>(target_node)];
  if (!target.is_leaf()) {
    build_work(target_tree, source_tree, target.left, source_node, theta,
               exclude_same_index);
    build_work(target_tree, source_tree, target.right, source_node, theta,
               exclude_same_index);
    return;
  }

  TargetLeafWork &work = target_tree.work[static_cast<std::size_t>(target_node)];
  const Cell &source = source_tree.cells[static_cast<std::size_t>(source_node)];
  if (accept_far(target_tree, source_tree, target_node, source_node, theta,
                 exclude_same_index)) {
    work.far_sources.push_back(source_node);
    return;
  }
  if (source.is_leaf()) {
    work.near_sources.push_back(source_node);
    return;
  }
  build_work(target_tree, source_tree, target_node, source.left, theta,
             exclude_same_index);
  build_work(target_tree, source_tree, target_node, source.right, theta,
             exclude_same_index);
}

void sort_work(Tree &tree) {
  for (int node : tree.target_leaf_nodes) {
    auto &work = tree.work[static_cast<std::size_t>(node)];
    std::sort(work.near_sources.begin(), work.near_sources.end());
    work.near_sources.erase(std::unique(work.near_sources.begin(), work.near_sources.end()),
                            work.near_sources.end());
    std::sort(work.far_sources.begin(), work.far_sources.end());
    work.far_sources.erase(std::unique(work.far_sources.begin(), work.far_sources.end()),
                           work.far_sources.end());
  }
}

std::string operator_key(const Tree &target_tree, const Tree &source_tree,
                         int target_node, int source_node,
                         const PrototypeConfig &config) {
  const Cell &t = target_tree.cells[static_cast<std::size_t>(target_node)];
  const Cell &s = source_tree.cells[static_cast<std::size_t>(source_node)];
  const double scale = std::max(1.0e-12, t.radius + s.radius);
  const auto bucket = [&](double value) {
    return static_cast<long long>(std::llround(value / scale * 16.0));
  };
  const long long ratio =
      static_cast<long long>(std::llround((s.radius / std::max(1.0e-12, t.radius)) * 8.0));
  return std::to_string(static_cast<int>(config.kernel)) + ":" +
         std::to_string(config.expansion_order) + ":" + std::to_string(bucket(s.cx - t.cx)) +
         ":" + std::to_string(bucket(s.cy - t.cy)) + ":" +
         std::to_string(bucket(s.cz - t.cz)) + ":" + std::to_string(ratio);
}

void compute_counts_and_operators(const Tree &target_tree, const Tree &source_tree,
                                  const PrototypeConfig &config,
                                  PrototypeMetrics &metrics) {
  std::set<std::string> keys;
  for (int tnode : target_tree.target_leaf_nodes) {
    const TargetLeafWork &work = target_tree.work[static_cast<std::size_t>(tnode)];
    const Cell &target = target_tree.cells[static_cast<std::size_t>(tnode)];
    for (int source_node : work.near_sources) {
      metrics.near_pairs +=
          target.count() * source_tree.cells[static_cast<std::size_t>(source_node)].count();
    }
    for (int source_node : work.far_sources) {
      metrics.far_pairs +=
          target.count() * source_tree.cells[static_cast<std::size_t>(source_node)].count();
      keys.insert(operator_key(target_tree, source_tree, tnode, source_node, config));
    }
  }
  metrics.m2l_operators = keys.size();
}

void apply_direct_reference(const PointSoA &sources, const PointSoA &targets,
                            const PrototypeConfig &config,
                            AlignedDoubleVector &out,
                            bool exclude_same_index = false) {
  out.assign(targets.x.size(), 0.0);
  for (std::size_t i = 0; i < targets.x.size(); ++i) {
    double sum = 0.0;
    for (std::size_t j = 0; j < sources.x.size(); ++j) {
      if (exclude_same_index && targets.panel_id[i] == sources.panel_id[j]) {
        continue;
      }
      const double dx = targets.x[i] - sources.x[j];
      const double dy = targets.y[i] - sources.y[j];
      const double dz = targets.z[i] - sources.z[j];
      const double r2 = dx * dx + dy * dy + dz * dz;
      sum += sources.q[j] * kernel_value(config.kernel, config.kappa, r2);
    }
    out[i] = sum;
  }
}

void apply_prototype(const PointSoA &sources, const Tree &source_tree,
                     const PointSoA &targets, const Tree &target_tree,
                     const PrototypeConfig &config, AlignedDoubleVector &out,
                     bool exclude_same_index = false) {
  out.assign(targets.x.size(), 0.0);
  TitanExecutor executor(static_cast<std::size_t>(std::max(1, config.threads)));
  const std::size_t workers = std::max<std::size_t>(1, executor.capacity());
  executor.parallel_for(target_tree.target_leaf_nodes.size(), workers, 1,
                        [&](std::size_t, std::size_t begin, std::size_t end) {
                          for (std::size_t leaf_pos = begin; leaf_pos < end; ++leaf_pos) {
                            const int target_node = target_tree.target_leaf_nodes[leaf_pos];
                            const Cell &target =
                                target_tree.cells[static_cast<std::size_t>(target_node)];
                            const TargetLeafWork &work =
                                target_tree.work[static_cast<std::size_t>(target_node)];
                            for (std::size_t tpos = target.begin; tpos < target.end; ++tpos) {
                              const std::size_t ti = target_tree.order[tpos];
                              double sum = 0.0;
                              for (int source_node : work.near_sources) {
                                const Cell &source =
                                    source_tree.cells[static_cast<std::size_t>(source_node)];
                                for (std::size_t spos = source.begin; spos < source.end; ++spos) {
                                  const std::size_t si = source_tree.order[spos];
                                  if (exclude_same_index &&
                                      targets.panel_id[ti] == sources.panel_id[si]) {
                                    continue;
                                  }
                                  const double dx = targets.x[ti] - sources.x[si];
                                  const double dy = targets.y[ti] - sources.y[si];
                                  const double dz = targets.z[ti] - sources.z[si];
                                  const double r2 = dx * dx + dy * dy + dz * dz;
                                  sum += sources.q[si] *
                                         kernel_value(config.kernel, config.kappa, r2);
                                }
                              }
                              for (int source_node : work.far_sources) {
                                const Cell &source =
                                    source_tree.cells[static_cast<std::size_t>(source_node)];
                                const double dx = targets.x[ti] - source.qx;
                                const double dy = targets.y[ti] - source.qy;
                                const double dz = targets.z[ti] - source.qz;
                                const double r2 = dx * dx + dy * dy + dz * dz;
                                const double g =
                                    kernel_value(config.kernel, config.kappa, r2);
                                double far = source.charge * g;
                                if (config.expansion_order >= 2 && r2 > 1.0e-28) {
                                  const double r = std::sqrt(r2);
                                  const double dot =
                                      source.px * dx + source.py * dy + source.pz * dz;
                                  const double dipole_scale =
                                      (config.kernel == KernelKind::Yukawa)
                                          ? std::exp(-std::max(0.0, config.kappa) * r) *
                                                (std::max(0.0, config.kappa) / (r * r) +
                                                 1.0 / (r * r * r))
                                          : 1.0 / (r * r * r);
                                  far += dot * dipole_scale;
                                }
                                sum += far;
                              }
                              out[ti] = sum;
                            }
                          }
                        });
}

void compute_error(const AlignedDoubleVector &reference,
                   const AlignedDoubleVector &approx, PrototypeMetrics &metrics) {
  double err2 = 0.0;
  double ref2 = 0.0;
  double max_abs = 0.0;
  bool finite = true;
  for (std::size_t i = 0; i < reference.size(); ++i) {
    const double err = approx[i] - reference[i];
    err2 += err * err;
    ref2 += reference[i] * reference[i];
    max_abs = std::max(max_abs, std::abs(err));
    finite = finite && std::isfinite(reference[i]) && std::isfinite(approx[i]);
  }
  const double inv = reference.empty() ? 0.0 : 1.0 / static_cast<double>(reference.size());
  metrics.rms_abs_error = std::sqrt(err2 * inv);
  metrics.max_abs_error = max_abs;
  metrics.rel_rms_error = std::sqrt(err2 / std::max(ref2, 1.0e-30));
  metrics.finite_outputs = finite;
}

bool finite_vector(const AlignedDoubleVector &values) {
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); });
}

struct PanelSourceMoments {
  double q_phin = 0.0;
  double phin_px = 0.0;
  double phin_py = 0.0;
  double phin_pz = 0.0;
  double phi_nx = 0.0;
  double phi_ny = 0.0;
  double phi_nz = 0.0;
};

double jcp_phin_for_panel(const PrototypePanel &panel, std::size_t index) {
  const double t = static_cast<double>(index + 1);
  return 0.65 * std::sin(0.017 * t) + 0.35 * panel.nz;
}

PointSoA make_panel_geometry_points(const std::vector<PrototypePanel> &panels) {
  PointSoA points;
  points.x.reserve(panels.size());
  points.y.reserve(panels.size());
  points.z.reserve(panels.size());
  points.q.assign(panels.size(), 1.0);
  points.panel_id.reserve(panels.size());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    points.x.push_back(panels[i].x);
    points.y.push_back(panels[i].y);
    points.z.push_back(panels[i].z);
    points.panel_id.push_back(i);
  }
  return points;
}

std::vector<PanelSourceMoments>
compute_panel_source_moments(const Tree &tree,
                             const std::vector<PrototypePanel> &panels,
                             const std::vector<double> &phi_values,
                             const std::vector<double> &phin_values) {
  std::vector<PanelSourceMoments> moments(tree.cells.size());
  for (std::size_t node_index = 0; node_index < tree.cells.size(); ++node_index) {
    const Cell &cell = tree.cells[node_index];
    PanelSourceMoments acc;
    for (std::size_t pos = cell.begin; pos < cell.end; ++pos) {
      const std::size_t source_index = tree.order[pos];
      const PrototypePanel &source = panels[source_index];
      const double phin = phin_values[source_index];
      const double q_phin = source.area * phin;
      acc.q_phin += q_phin;
      acc.phin_px += q_phin * (source.x - cell.qx);
      acc.phin_py += q_phin * (source.y - cell.qy);
      acc.phin_pz += q_phin * (source.z - cell.qz);
      const double q_phi = source.area * phi_values[source_index];
      acc.phi_nx += q_phi * source.nx;
      acc.phi_ny += q_phi * source.ny;
      acc.phi_nz += q_phi * source.nz;
    }
    moments[node_index] = acc;
  }
  return moments;
}

void add_jcp_direct_pair(const PrototypePanel &target,
                         const PrototypePanel &source, double source_phi,
                         double source_phin,
                         double eps_ratio, double kappa, double &phi_out,
                         double &phin_out) {
  constexpr double inv_four_pi = 0.079577471545947667884;
  const double safe_eps_ratio = std::max(1.0e-12, eps_ratio);
  const double inv_eps_ratio = 1.0 / safe_eps_ratio;
  const double kk = std::max(0.0, kappa);
  const double ux = target.x - source.x;
  const double uy = target.y - source.y;
  const double uz = target.z - source.z;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (!(r2 > 1.0e-28)) {
    return;
  }
  const double r = std::sqrt(r2);
  const double inv_r = 1.0 / r;
  const double inv_r2 = inv_r * inv_r;
  const double inv_r3 = inv_r2 * inv_r;
  const double expkr = std::exp(-kk * r);
  const double one_plus_kr = 1.0 + kk * r;
  const double source_projection =
      source.nx * ux + source.ny * uy + source.nz * uz;
  const double target_projection =
      target.nx * ux + target.ny * uy + target.nz * uz;
  const double phin = source_phin;
  const double g0 = inv_four_pi * inv_r;
  const double gk = expkr * g0;
  const double source_d0 = inv_four_pi * source_projection * inv_r3;
  const double source_dk = expkr * one_plus_kr * source_d0;
  const double target_d0 = -inv_four_pi * target_projection * inv_r3;
  const double target_dk = expkr * one_plus_kr * target_d0;
  const double nxu = target_projection * inv_r;
  const double nyu = source_projection * inv_r;
  const double nxny =
      target.nx * source.nx + target.ny * source.ny + target.nz * source.nz;
  const double nxu_nyu = nxu * nyu;
  const double tangent_metric = nxny - nxu_nyu;
  const double mixed0 = inv_four_pi * inv_r3 * (3.0 * nxu_nyu - nxny);
  const double fprime_k_over_r = -expkr * one_plus_kr * inv_four_pi * inv_r3;
  const double fsecond_k =
      expkr * inv_four_pi *
      (kk * kk * inv_r + 2.0 * kk * inv_r2 + 2.0 * inv_r3);
  const double mixedk = fsecond_k * nxu_nyu + fprime_k_over_r * tangent_metric;
  const double phi = source_phi;
  phi_out += source.area * ((g0 - gk) * phin +
                            (safe_eps_ratio * source_dk - source_d0) * phi);
  phin_out += source.area *
              ((target_d0 - inv_eps_ratio * target_dk) * phin +
               (mixedk - mixed0) * phi);
}

void add_jcp_direct_pair(const PrototypePanel &target,
                         const PrototypePanel &source, double source_phin,
                         double eps_ratio, double kappa, double &phi_out,
                         double &phin_out) {
  add_jcp_direct_pair(target, source, source.scalar_density, source_phin,
                      eps_ratio, kappa, phi_out, phin_out);
}

void add_jcp_far_cell(const PrototypePanel &target, const Cell &source_cell,
                      const PanelSourceMoments &moments, double eps_ratio,
                      double kappa, double &phi_out, double &phin_out) {
  constexpr double inv_four_pi = 0.079577471545947667884;
  const double safe_eps_ratio = std::max(1.0e-12, eps_ratio);
  const double inv_eps_ratio = 1.0 / safe_eps_ratio;
  const double kk = std::max(0.0, kappa);
  const double ux = target.x - source_cell.qx;
  const double uy = target.y - source_cell.qy;
  const double uz = target.z - source_cell.qz;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (!(r2 > 1.0e-28)) {
    return;
  }
  const double r = std::sqrt(r2);
  const double inv_r = 1.0 / r;
  const double inv_r2 = inv_r * inv_r;
  const double inv_r3 = inv_r2 * inv_r;
  const double expkr = std::exp(-kk * r);
  const double one_plus_kr = 1.0 + kk * r;
  const double g0 = inv_four_pi * inv_r;
  const double gk = expkr * g0;
  const double dipole_dot =
      moments.phin_px * ux + moments.phin_py * uy + moments.phin_pz * uz;
  const double lap_dipole_scale = inv_four_pi * inv_r3;
  const double yuk_dipole_scale = expkr * inv_four_pi *
                                  (kk * inv_r2 + inv_r3);
  const double k1 = moments.q_phin * (g0 - gk) +
                    dipole_dot * (lap_dipole_scale - yuk_dipole_scale);

  const double source_normal_projection =
      moments.phi_nx * ux + moments.phi_ny * uy + moments.phi_nz * uz;
  const double source_d0 = inv_four_pi * source_normal_projection * inv_r3;
  const double source_dk = expkr * one_plus_kr * source_d0;
  const double k2 = safe_eps_ratio * source_dk - source_d0;

  const double target_projection =
      target.nx * ux + target.ny * uy + target.nz * uz;
  const double target_d0 =
      -inv_four_pi * target_projection * inv_r3 * moments.q_phin;
  const double target_dk = expkr * one_plus_kr * target_d0;
  const double k3 = target_d0 - inv_eps_ratio * target_dk;

  const double nxu = target_projection * inv_r;
  const double nyu = source_normal_projection * inv_r;
  const double nxny =
      target.nx * moments.phi_nx + target.ny * moments.phi_ny +
      target.nz * moments.phi_nz;
  const double nxu_nyu = nxu * nyu;
  const double tangent_metric = nxny - nxu_nyu;
  const double mixed0 = inv_four_pi * inv_r3 * (3.0 * nxu_nyu - nxny);
  const double fprime_k_over_r = -expkr * one_plus_kr * inv_four_pi * inv_r3;
  const double fsecond_k =
      expkr * inv_four_pi *
      (kk * kk * inv_r + 2.0 * kk * inv_r2 + 2.0 * inv_r3);
  const double mixedk = fsecond_k * nxu_nyu + fprime_k_over_r * tangent_metric;
  const double k4 = mixedk - mixed0;

  phi_out += k1 + k2;
  phin_out += k3 + k4;
}

} // namespace

namespace {

PrototypeConfig normalize_jcp_analytic_config(PrototypeConfig config) {
  config.max_leaf_size = std::max(1, config.max_leaf_size);
  config.threads = std::max(1, config.threads);
  return config;
}

bool same_double(double a, double b) noexcept {
  return std::abs(a - b) <=
         1.0e-12 * std::max({1.0, std::abs(a), std::abs(b)});
}

} // namespace

struct FmmM2lAnalyticContext::Impl {
  PrototypeConfig config;
  std::vector<PrototypePanel> panels;
  PointSoA points;
  Tree source_tree;
  Tree target_tree;
  std::size_t near_pairs{0};
  std::size_t far_pairs{0};
  std::size_t work_items{0};
  double eps_ratio{1.0};
  double kappa{0.0};
  double setup_seconds{0.0};
};

FmmM2lAnalyticContext::FmmM2lAnalyticContext(
    const PrototypeConfig &input_config, std::vector<PrototypePanel> panels,
    double eps_ratio, double kappa)
    : impl_(std::make_unique<Impl>()) {
  if (panels.empty()) {
    throw std::invalid_argument(
        "JCP analytic FMM/M2L context requires at least one panel");
  }
  if (!std::isfinite(kappa) || kappa < 0.0) {
    throw std::invalid_argument(
        "JCP analytic FMM/M2L context requires finite kappa >= 0");
  }
  const auto setup_begin = Clock::now();
  impl_->config = normalize_jcp_analytic_config(input_config);
  impl_->panels = std::move(panels);
  impl_->eps_ratio = eps_ratio;
  impl_->kappa = kappa;
  impl_->points = make_panel_geometry_points(impl_->panels);
  impl_->source_tree = build_tree(
      impl_->points, static_cast<std::size_t>(impl_->config.max_leaf_size));
  impl_->target_tree = build_tree(
      impl_->points, static_cast<std::size_t>(impl_->config.max_leaf_size));
  if (!impl_->source_tree.cells.empty() && !impl_->target_tree.cells.empty()) {
    build_work(impl_->target_tree, impl_->source_tree, 0, 0,
               impl_->config.separation_theta, true);
  }
  sort_work(impl_->target_tree);
  PrototypeMetrics counts;
  compute_counts_and_operators(impl_->target_tree, impl_->source_tree,
                               impl_->config, counts);
  impl_->near_pairs = counts.near_pairs;
  impl_->far_pairs = counts.far_pairs;
  for (int target_node : impl_->target_tree.target_leaf_nodes) {
    const TargetLeafWork &work =
        impl_->target_tree.work[static_cast<std::size_t>(target_node)];
    impl_->work_items += work.near_sources.size() + work.far_sources.size();
  }
  impl_->setup_seconds = elapsed_seconds(setup_begin, Clock::now());
}

FmmM2lAnalyticContext::~FmmM2lAnalyticContext() = default;

FmmM2lAnalyticContext::FmmM2lAnalyticContext(
    FmmM2lAnalyticContext &&) noexcept = default;

FmmM2lAnalyticContext &
FmmM2lAnalyticContext::operator=(FmmM2lAnalyticContext &&) noexcept = default;

bool FmmM2lAnalyticContext::matches(const PrototypeConfig &input_config,
                                    std::size_t panel_count, double eps_ratio,
                                    double kappa) const noexcept {
  if (impl_ == nullptr) {
    return false;
  }
  const PrototypeConfig config = normalize_jcp_analytic_config(input_config);
  return panel_count == impl_->panels.size() &&
         config.kernel == impl_->config.kernel &&
         config.max_leaf_size == impl_->config.max_leaf_size &&
         config.expansion_order == impl_->config.expansion_order &&
         same_double(config.separation_theta, impl_->config.separation_theta) &&
         same_double(eps_ratio, impl_->eps_ratio) &&
         same_double(kappa, impl_->kappa);
}

std::size_t FmmM2lAnalyticContext::panel_count() const noexcept {
  return impl_ == nullptr ? 0 : impl_->panels.size();
}

std::size_t FmmM2lAnalyticContext::context_nodes() const noexcept {
  return impl_ == nullptr
             ? 0
             : impl_->source_tree.cells.size() + impl_->target_tree.cells.size();
}

std::size_t FmmM2lAnalyticContext::context_work_items() const noexcept {
  return impl_ == nullptr ? 0 : impl_->work_items;
}

double FmmM2lAnalyticContext::setup_seconds() const noexcept {
  return impl_ == nullptr ? 0.0 : impl_->setup_seconds;
}

std::vector<std::vector<std::size_t>>
FmmM2lAnalyticContext::block_panel_ids(std::size_t max_block_size) const {
  if (impl_ == nullptr) {
    throw std::logic_error("JCP analytic FMM/M2L context is not initialized");
  }
  const std::size_t cap = std::max<std::size_t>(1, max_block_size);
  std::vector<std::vector<std::size_t>> blocks;
  for (const int node_index : impl_->target_tree.target_leaf_nodes) {
    if (node_index < 0 ||
        static_cast<std::size_t>(node_index) >= impl_->target_tree.cells.size()) {
      continue;
    }
    const Cell &cell =
        impl_->target_tree.cells[static_cast<std::size_t>(node_index)];
    const std::size_t end = std::min(cell.end, impl_->target_tree.order.size());
    for (std::size_t begin = cell.begin; begin < end; begin += cap) {
      const std::size_t block_end = std::min(end, begin + cap);
      std::vector<std::size_t> block;
      block.reserve(block_end - begin);
      for (std::size_t pos = begin; pos < block_end; ++pos) {
        const std::size_t panel = impl_->target_tree.order[pos];
        if (panel < impl_->panels.size()) {
          block.push_back(panel);
        }
      }
      if (!block.empty()) {
        blocks.push_back(std::move(block));
      }
    }
  }
  return blocks;
}

std::vector<FmmM2lOverlapPanelGroup>
FmmM2lAnalyticContext::overlap_panel_groups(
    std::size_t owned_block_size, std::size_t overlap_panel_cap) const {
  if (impl_ == nullptr) {
    throw std::logic_error("JCP analytic FMM/M2L context is not initialized");
  }
  const std::size_t owned_cap = std::max<std::size_t>(1, owned_block_size);
  const std::size_t halo_cap = overlap_panel_cap;
  std::vector<FmmM2lOverlapPanelGroup> groups;
  const std::size_t n = impl_->panels.size();

  for (const int node_index : impl_->target_tree.target_leaf_nodes) {
    if (node_index < 0 ||
        static_cast<std::size_t>(node_index) >= impl_->target_tree.cells.size()) {
      continue;
    }
    const Cell &cell =
        impl_->target_tree.cells[static_cast<std::size_t>(node_index)];
    const std::size_t end = std::min(cell.end, impl_->target_tree.order.size());
    const TargetLeafWork &work =
        impl_->target_tree.work[static_cast<std::size_t>(node_index)];

    for (std::size_t begin = cell.begin; begin < end; begin += owned_cap) {
      const std::size_t block_end = std::min(end, begin + owned_cap);
      FmmM2lOverlapPanelGroup group;
      group.owned_panels.reserve(block_end - begin);
      group.factor_panels.reserve(block_end - begin + halo_cap);
      std::vector<unsigned char> owned_mask(n, 0);
      std::vector<unsigned char> factor_mask(n, 0);

      double cx = 0.0;
      double cy = 0.0;
      double cz = 0.0;
      for (std::size_t pos = begin; pos < block_end; ++pos) {
        const std::size_t panel = impl_->target_tree.order[pos];
        if (panel >= n || owned_mask[panel])
          continue;
        owned_mask[panel] = 1;
        factor_mask[panel] = 1;
        group.owned_panels.push_back(panel);
        group.factor_panels.push_back(panel);
        cx += impl_->panels[panel].x;
        cy += impl_->panels[panel].y;
        cz += impl_->panels[panel].z;
      }
      if (group.owned_panels.empty())
        continue;
      const double inv_owned =
          1.0 / static_cast<double>(group.owned_panels.size());
      cx *= inv_owned;
      cy *= inv_owned;
      cz *= inv_owned;

      if (halo_cap > 0) {
        std::vector<std::pair<double, std::size_t>> candidates;
        for (const int source_node : work.near_sources) {
          if (source_node < 0 ||
              static_cast<std::size_t>(source_node) >=
                  impl_->source_tree.cells.size()) {
            continue;
          }
          const Cell &source_cell =
              impl_->source_tree.cells[static_cast<std::size_t>(source_node)];
          const std::size_t source_end =
              std::min(source_cell.end, impl_->source_tree.order.size());
          for (std::size_t pos = source_cell.begin; pos < source_end; ++pos) {
            const std::size_t panel = impl_->source_tree.order[pos];
            if (panel >= n || owned_mask[panel] || factor_mask[panel])
              continue;
            factor_mask[panel] = 1;
            const double dx = impl_->panels[panel].x - cx;
            const double dy = impl_->panels[panel].y - cy;
            const double dz = impl_->panels[panel].z - cz;
            candidates.push_back({dx * dx + dy * dy + dz * dz, panel});
          }
        }
        std::sort(candidates.begin(), candidates.end(),
                  [](const auto &lhs, const auto &rhs) {
                    if (lhs.first != rhs.first)
                      return lhs.first < rhs.first;
                    return lhs.second < rhs.second;
                  });
        const std::size_t selected = std::min(halo_cap, candidates.size());
        group.halo_count = selected;
        for (std::size_t i = 0; i < selected; ++i)
          group.factor_panels.push_back(candidates[i].second);
      }
      groups.push_back(std::move(group));
    }
  }
  return groups;
}

std::vector<std::vector<std::size_t>>
FmmM2lAnalyticContext::near_panel_graph(std::size_t neighbor_cap) const {
  if (impl_ == nullptr) {
    throw std::logic_error("JCP analytic FMM/M2L context is not initialized");
  }
  const std::size_t n = impl_->panels.size();
  const std::size_t cap = std::max<std::size_t>(1, neighbor_cap);
  std::vector<std::vector<std::size_t>> graph(n);

  for (const int node_index : impl_->target_tree.target_leaf_nodes) {
    if (node_index < 0 ||
        static_cast<std::size_t>(node_index) >= impl_->target_tree.cells.size()) {
      continue;
    }
    const Cell &target_cell =
        impl_->target_tree.cells[static_cast<std::size_t>(node_index)];
    const std::size_t target_end =
        std::min(target_cell.end, impl_->target_tree.order.size());
    const TargetLeafWork &work =
        impl_->target_tree.work[static_cast<std::size_t>(node_index)];

    for (std::size_t tpos = target_cell.begin; tpos < target_end; ++tpos) {
      const std::size_t target_panel = impl_->target_tree.order[tpos];
      if (target_panel >= n)
        continue;
      std::vector<unsigned char> seen(n, 0);
      std::vector<std::pair<double, std::size_t>> candidates;
      auto add_candidate = [&](std::size_t source_panel) {
        if (source_panel >= n || seen[source_panel])
          return;
        seen[source_panel] = 1;
        const double dx = impl_->panels[source_panel].x -
                          impl_->panels[target_panel].x;
        const double dy = impl_->panels[source_panel].y -
                          impl_->panels[target_panel].y;
        const double dz = impl_->panels[source_panel].z -
                          impl_->panels[target_panel].z;
        candidates.push_back({dx * dx + dy * dy + dz * dz, source_panel});
      };
      add_candidate(target_panel);
      for (const int source_node : work.near_sources) {
        if (source_node < 0 ||
            static_cast<std::size_t>(source_node) >=
                impl_->source_tree.cells.size()) {
          continue;
        }
        const Cell &source_cell =
            impl_->source_tree.cells[static_cast<std::size_t>(source_node)];
        const std::size_t source_end =
            std::min(source_cell.end, impl_->source_tree.order.size());
        for (std::size_t spos = source_cell.begin; spos < source_end; ++spos)
          add_candidate(impl_->source_tree.order[spos]);
      }
      std::sort(candidates.begin(), candidates.end(),
                [](const auto &lhs, const auto &rhs) {
                  if (lhs.first != rhs.first)
                    return lhs.first < rhs.first;
                  return lhs.second < rhs.second;
                });
      std::vector<std::size_t> row;
      const std::size_t selected = std::min(cap, candidates.size());
      row.reserve(selected);
      for (std::size_t i = 0; i < selected; ++i)
        row.push_back(candidates[i].second);
      std::sort(row.begin(), row.end());
      row.erase(std::unique(row.begin(), row.end()), row.end());
      graph[target_panel] = std::move(row);
    }
  }
  for (std::size_t panel = 0; panel < n; ++panel) {
    if (graph[panel].empty())
      graph[panel].push_back(panel);
  }
  return graph;
}

JcpAnalyticOperatorResult
FmmM2lAnalyticContext::apply(const std::vector<double> &phi_values,
                             const std::vector<double> &phin_values,
                             TitanExecutor &executor) const {
  if (impl_ == nullptr) {
    throw std::logic_error("JCP analytic FMM/M2L context is not initialized");
  }
  if (phi_values.size() != impl_->panels.size() ||
      phin_values.size() != impl_->panels.size()) {
    throw std::invalid_argument(
        "JCP analytic FMM/M2L context density vectors must match panels");
  }

  JcpAnalyticOperatorResult result;
  result.panels = impl_->panels.size();
  result.near_pairs = impl_->near_pairs;
  result.far_pairs = impl_->far_pairs;
  result.setup_seconds = impl_->setup_seconds;
  result.context_nodes = context_nodes();
  result.context_work_items = context_work_items();
  result.phi.assign(impl_->panels.size(), 0.0);
  result.phin.assign(impl_->panels.size(), 0.0);

  const std::vector<PanelSourceMoments> source_moments =
      compute_panel_source_moments(impl_->source_tree, impl_->panels,
                                   phi_values, phin_values);

  const auto apply_begin = Clock::now();
  const std::size_t workers = std::max<std::size_t>(1, executor.capacity());
  executor.parallel_for(
      impl_->target_tree.target_leaf_nodes.size(), workers, 1,
      [&](std::size_t, std::size_t begin, std::size_t end) {
        for (std::size_t leaf_pos = begin; leaf_pos < end; ++leaf_pos) {
          const int target_node = impl_->target_tree.target_leaf_nodes[leaf_pos];
          const Cell &target_cell =
              impl_->target_tree.cells[static_cast<std::size_t>(target_node)];
          const TargetLeafWork &work =
              impl_->target_tree.work[static_cast<std::size_t>(target_node)];
          for (std::size_t tpos = target_cell.begin; tpos < target_cell.end;
               ++tpos) {
            const std::size_t ti = impl_->target_tree.order[tpos];
            double phi_sum = 0.0;
            double phin_sum = 0.0;
            const PrototypePanel &target = impl_->panels[ti];
            for (int source_node : work.near_sources) {
              const Cell &source_cell =
                  impl_->source_tree.cells[static_cast<std::size_t>(source_node)];
              for (std::size_t spos = source_cell.begin; spos < source_cell.end;
                   ++spos) {
                const std::size_t si = impl_->source_tree.order[spos];
                if (ti == si) {
                  continue;
                }
                add_jcp_direct_pair(target, impl_->panels[si], phi_values[si],
                                    phin_values[si], impl_->eps_ratio,
                                    impl_->kappa, phi_sum, phin_sum);
              }
            }
            for (int source_node : work.far_sources) {
              const std::size_t source_index =
                  static_cast<std::size_t>(source_node);
              add_jcp_far_cell(target, impl_->source_tree.cells[source_index],
                               source_moments[source_index], impl_->eps_ratio,
                               impl_->kappa, phi_sum, phin_sum);
            }
            result.phi[ti] = phi_sum;
            result.phin[ti] = phin_sum;
          }
        }
      });
  result.apply_seconds = elapsed_seconds(apply_begin, Clock::now());

  result.finite_outputs =
      std::all_of(result.phi.begin(), result.phi.end(),
                  [](double value) { return std::isfinite(value); }) &&
      std::all_of(result.phin.begin(), result.phin.end(),
                  [](double value) { return std::isfinite(value); });
  return result;
}

const char *kernel_kind_name(KernelKind kernel) noexcept {
  switch (kernel) {
  case KernelKind::Laplace:
    return "laplace";
  case KernelKind::Yukawa:
    return "yukawa";
  }
  return "unknown";
}

PrototypeEvaluation evaluate_source_target_cloud(
    const PrototypeConfig &input_config, const PrototypePointCloud &source_points,
    const PrototypePointCloud &target_points, bool exclude_same_index = false,
    bool compute_direct_reference = true) {
  const PointSoA sources = from_public_point_cloud(source_points);
  const PointSoA targets = from_public_point_cloud(target_points);
  PrototypeConfig config = input_config;
  config.max_leaf_size = std::max(1, config.max_leaf_size);
  config.threads = std::max(1, config.threads);
  if (config.kernel == KernelKind::Yukawa && !(config.kappa > 0.0)) {
    config.kappa = 0.125;
  }

  PrototypeEvaluation evaluation;
  PrototypeMetrics &metrics = evaluation.metrics;
  metrics.sources = sources.x.size();
  metrics.targets = targets.x.size();

  const auto setup_begin = Clock::now();
  Tree source_tree = build_tree(sources, static_cast<std::size_t>(config.max_leaf_size));
  Tree target_tree = build_tree(targets, static_cast<std::size_t>(config.max_leaf_size));
  if (!source_tree.cells.empty() && !target_tree.cells.empty()) {
    build_work(target_tree, source_tree, 0, 0, config.separation_theta,
               exclude_same_index);
  }
  sort_work(target_tree);
  metrics.cells = source_tree.cells.size() + target_tree.cells.size();
  compute_counts_and_operators(target_tree, source_tree, config, metrics);
  const auto setup_end = Clock::now();
  metrics.setup_seconds = elapsed_seconds(setup_begin, setup_end);

  if (compute_direct_reference) {
    const auto direct_begin = Clock::now();
    apply_direct_reference(sources, targets, config, evaluation.reference,
                           exclude_same_index);
    const auto direct_end = Clock::now();
    metrics.direct_seconds = elapsed_seconds(direct_begin, direct_end);
  }

  const auto apply_begin = Clock::now();
  apply_prototype(sources, source_tree, targets, target_tree, config,
                  evaluation.approx, exclude_same_index);
  const auto apply_end = Clock::now();
  metrics.apply_seconds = elapsed_seconds(apply_begin, apply_end);

  if (compute_direct_reference) {
    compute_error(evaluation.reference, evaluation.approx, metrics);
  } else {
    metrics.finite_outputs = finite_vector(evaluation.approx);
  }
  return evaluation;
}

PrototypeEvaluation evaluate_point_cloud(const PrototypeConfig &input_config,
                                         const PrototypePointCloud &input_points) {
  return evaluate_source_target_cloud(input_config, input_points, input_points);
}

PrototypeMetrics run_point_cloud_probe(const PrototypeConfig &input_config,
                                       const PrototypePointCloud &input_points) {
  return evaluate_point_cloud(input_config, input_points).metrics;
}

PrototypeMetrics run_synthetic_probe(const PrototypeConfig &config,
                                     std::size_t panel_count, unsigned seed) {
  if (panel_count == 0) {
    throw std::invalid_argument("FMM/M2L prototype requires at least one panel");
  }
  const PointSoA generated = make_synthetic_points(panel_count, seed);
  PrototypePointCloud public_points;
  public_points.x.assign(generated.x.begin(), generated.x.end());
  public_points.y.assign(generated.y.begin(), generated.y.end());
  public_points.z.assign(generated.z.begin(), generated.z.end());
  public_points.q.assign(generated.q.begin(), generated.q.end());
  return run_point_cloud_probe(config, public_points);
}

PrototypePointCloud make_mini_protein_like_point_cloud(std::size_t panel_count,
                                                       unsigned seed) {
  if (panel_count == 0) {
    throw std::invalid_argument("mini-protein-like fixture requires at least one point");
  }
  PrototypePointCloud points;
  points.x.resize(panel_count);
  points.y.resize(panel_count);
  points.z.resize(panel_count);
  points.q.resize(panel_count);

  std::uint64_t state = (static_cast<std::uint64_t>(seed) << 32) ^ 0xd1b54a32d192ed03ull;
  auto next01 = [&]() {
    state = state * 2862933555777941757ull + 3037000493ull;
    const std::uint64_t bits = (state >> 11) | 0x3ff0000000000000ull;
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value - 1.0;
  };

  for (std::size_t i = 0; i < panel_count; ++i) {
    const double t = 12.5663706143591729539 *
                     static_cast<double>(i) / static_cast<double>(panel_count);
    const double strand = static_cast<double>(i % 5) - 2.0;
    const double shell = 0.55 + 0.25 * next01();
    const double wobble = 0.15 * (next01() - 0.5);
    points.x[i] = 2.8 * std::cos(0.5 * t) + shell * std::cos(t) + 0.18 * strand;
    points.y[i] = 1.4 * std::sin(0.5 * t) + shell * std::sin(t) + wobble;
    points.z[i] = 0.22 * static_cast<double>(i % 31) + 0.45 * std::sin(1.5 * t);
    points.q[i] = ((i % 7) == 0 ? 1.0 : ((i % 11) == 0 ? -1.0 : 0.12 * (next01() - 0.5)));
  }
  const double z_mean =
      std::accumulate(points.z.begin(), points.z.end(), 0.0) /
      static_cast<double>(panel_count);
  for (double &z : points.z) {
    z -= z_mean;
  }
  return points;
}

std::vector<PrototypePanel>
make_tabi_panel_like_fixture(std::size_t panel_count, unsigned seed) {
  if (panel_count == 0) {
    throw std::invalid_argument("TABI-panel-like fixture requires at least one panel");
  }
  std::vector<PrototypePanel> panels(panel_count);
  std::uint64_t state = (static_cast<std::uint64_t>(seed) << 32) ^ 0x94d049bb133111ebull;
  auto next01 = [&]() {
    state = state * 3202034522624059733ull + 1ull;
    const std::uint64_t bits = (state >> 11) | 0x3ff0000000000000ull;
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value - 1.0;
  };

  const double golden = 2.3999632297286533222;
  for (std::size_t i = 0; i < panel_count; ++i) {
    const double u =
        -1.0 + 2.0 * (static_cast<double>(i) + 0.5) / static_cast<double>(panel_count);
    const double r = std::sqrt(std::max(0.0, 1.0 - u * u));
    const double theta = golden * static_cast<double>(i);
    const double ripple = 1.0 + 0.08 * std::sin(5.0 * theta) + 0.03 * (next01() - 0.5);
    PrototypePanel panel;
    panel.x = 2.2 * ripple * r * std::cos(theta);
    panel.y = 1.4 * ripple * r * std::sin(theta);
    panel.z = 3.0 * ripple * u;
    double nx = panel.x / (2.2 * 2.2);
    double ny = panel.y / (1.4 * 1.4);
    double nz = panel.z / (3.0 * 3.0);
    const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (norm > 0.0) {
      nx /= norm;
      ny /= norm;
      nz /= norm;
    }
    panel.nx = nx;
    panel.ny = ny;
    panel.nz = nz;
    panel.area = 4.0 * 3.14159265358979323846 * 2.2 * 1.4 /
                 static_cast<double>(panel_count) *
                 (0.8 + 0.4 * std::abs(nz));
    panel.scalar_density =
        ((i % 9) == 0 ? 1.0 : ((i % 13) == 0 ? -0.85 : 0.08 * std::sin(theta)));
    panels[i] = panel;
  }
  return panels;
}

PrototypePointCloud point_cloud_from_panels(const std::vector<PrototypePanel> &panels) {
  if (panels.empty()) {
    throw std::invalid_argument("FMM/M2L panel adapter requires at least one panel");
  }
  PrototypePointCloud points;
  points.x.reserve(panels.size());
  points.y.reserve(panels.size());
  points.z.reserve(panels.size());
  points.q.reserve(panels.size());
  for (const PrototypePanel &panel : panels) {
    if (!(panel.area > 0.0)) {
      throw std::invalid_argument("FMM/M2L panel adapter requires positive panel area");
    }
    points.x.push_back(panel.x);
    points.y.push_back(panel.y);
    points.z.push_back(panel.z);
    points.q.push_back(panel.area * panel.scalar_density);
  }
  return points;
}

std::array<PrototypePointCloud, 4>
make_jcp_channel_point_clouds(const std::vector<PrototypePanel> &panels) {
  if (panels.empty()) {
    throw std::invalid_argument("JCP-channel prototype requires at least one panel");
  }
  std::array<PrototypePointCloud, 4> channels;
  for (auto &cloud : channels) {
    cloud.x.reserve(panels.size());
    cloud.y.reserve(panels.size());
    cloud.z.reserve(panels.size());
    cloud.q.reserve(panels.size());
  }

  for (std::size_t i = 0; i < panels.size(); ++i) {
    const PrototypePanel &panel = panels[i];
    if (!(panel.area > 0.0)) {
      throw std::invalid_argument("JCP-channel prototype requires positive panel area");
    }
    const double t = static_cast<double>(i + 1);
    const double phi = panel.scalar_density;
    const double phin = 0.65 * std::sin(0.017 * t) + 0.35 * panel.nz;
    const double normal_mix = 0.45 * panel.nx - 0.25 * panel.ny + 0.15 * panel.nz;
    const std::array<double, 4> q{
        panel.area * phin,
        panel.area * phi * normal_mix,
        panel.area * phin * (1.0 + 0.2 * panel.nz),
        panel.area * phi * (0.5 + 0.5 * std::abs(panel.nz)),
    };
    for (std::size_t c = 0; c < channels.size(); ++c) {
      channels[c].x.push_back(panel.x);
      channels[c].y.push_back(panel.y);
      channels[c].z.push_back(panel.z);
      channels[c].q.push_back(q[c]);
    }
  }
  return channels;
}

PrototypePointCloud make_jcp_k1_source_cloud(const std::vector<PrototypePanel> &panels) {
  constexpr double inv_four_pi = 0.079577471545947667884;
  PrototypePointCloud cloud;
  cloud.x.reserve(panels.size());
  cloud.y.reserve(panels.size());
  cloud.z.reserve(panels.size());
  cloud.q.reserve(panels.size());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const PrototypePanel &panel = panels[i];
    cloud.x.push_back(panel.x);
    cloud.y.push_back(panel.y);
    cloud.z.push_back(panel.z);
    cloud.q.push_back(panel.area * jcp_phin_for_panel(panel, i) * inv_four_pi);
  }
  return cloud;
}

PrototypePointCloud make_panel_target_cloud(const std::vector<PrototypePanel> &panels,
                                            double normal_shift) {
  PrototypePointCloud cloud;
  cloud.x.reserve(panels.size());
  cloud.y.reserve(panels.size());
  cloud.z.reserve(panels.size());
  cloud.q.assign(panels.size(), 0.0);
  cloud.panel_id.reserve(panels.size());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const PrototypePanel &panel = panels[i];
    cloud.x.push_back(panel.x + normal_shift * panel.nx);
    cloud.y.push_back(panel.y + normal_shift * panel.ny);
    cloud.z.push_back(panel.z + normal_shift * panel.nz);
    cloud.panel_id.push_back(i);
  }
  return cloud;
}

PrototypePointCloud make_source_normal_dipole_cloud(
    const std::vector<PrototypePanel> &panels, bool use_phi, double h) {
  constexpr double inv_four_pi = 0.079577471545947667884;
  PrototypePointCloud cloud;
  cloud.x.reserve(panels.size() * 2);
  cloud.y.reserve(panels.size() * 2);
  cloud.z.reserve(panels.size() * 2);
  cloud.q.reserve(panels.size() * 2);
  cloud.panel_id.reserve(panels.size() * 2);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const PrototypePanel &panel = panels[i];
    const double value = use_phi ? panel.scalar_density : jcp_phin_for_panel(panel, i);
    const double q = panel.area * value * inv_four_pi / (2.0 * h);
    cloud.x.push_back(panel.x + h * panel.nx);
    cloud.y.push_back(panel.y + h * panel.ny);
    cloud.z.push_back(panel.z + h * panel.nz);
    cloud.q.push_back(q);
    cloud.panel_id.push_back(i);
    cloud.x.push_back(panel.x - h * panel.nx);
    cloud.y.push_back(panel.y - h * panel.ny);
    cloud.z.push_back(panel.z - h * panel.nz);
    cloud.q.push_back(-q);
    cloud.panel_id.push_back(i);
  }
  return cloud;
}

AlignedDoubleVector scalar_diff(const AlignedDoubleVector &lhs,
                                const AlignedDoubleVector &rhs, double scale) {
  AlignedDoubleVector out(lhs.size(), 0.0);
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    out[i] = (lhs[i] - rhs[i]) * scale;
  }
  return out;
}

struct DerivativeReference {
  AlignedDoubleVector k2;
  AlignedDoubleVector k3;
  AlignedDoubleVector k4;
};

DerivativeReference compute_direct_derivative_reference(
    const std::vector<PrototypePanel> &panels, double eps_ratio, double kappa) {
  constexpr double inv_four_pi = 0.079577471545947667884;
  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  const double kk = std::max(0.0, kappa);
  DerivativeReference ref;
  ref.k2.assign(panels.size(), 0.0);
  ref.k3.assign(panels.size(), 0.0);
  ref.k4.assign(panels.size(), 0.0);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const PrototypePanel &target = panels[i];
    for (std::size_t j = 0; j < panels.size(); ++j) {
      if (i == j) {
        continue;
      }
      const PrototypePanel &source = panels[j];
      const double ux = target.x - source.x;
      const double uy = target.y - source.y;
      const double uz = target.z - source.z;
      const double r2 = ux * ux + uy * uy + uz * uz;
      if (!(r2 > 1.0e-28)) {
        continue;
      }
      const double r = std::sqrt(r2);
      const double inv_r = 1.0 / r;
      const double inv_r2 = inv_r * inv_r;
      const double inv_r3 = inv_r2 * inv_r;
      const double source_projection = source.nx * ux + source.ny * uy + source.nz * uz;
      const double target_projection = target.nx * ux + target.ny * uy + target.nz * uz;
      const double source_d0 = inv_four_pi * source_projection * inv_r3;
      const double expkr = std::exp(-kk * r);
      const double one_plus_kr = 1.0 + kk * r;
      const double source_dk = expkr * one_plus_kr * source_d0;
      const double target_d0 = -inv_four_pi * target_projection * inv_r3;
      const double target_dk = expkr * one_plus_kr * target_d0;
      const double nxu = target_projection * inv_r;
      const double nyu = source_projection * inv_r;
      const double nxny =
          target.nx * source.nx + target.ny * source.ny + target.nz * source.nz;
      const double nxu_nyu = nxu * nyu;
      const double tangent_metric = nxny - nxu_nyu;
      const double mixed0 = inv_four_pi * inv_r3 * (3.0 * nxu_nyu - nxny);
      const double fprime_k_over_r = -expkr * one_plus_kr * inv_four_pi * inv_r3;
      const double fsecond_k =
          expkr * inv_four_pi * (kk * kk * inv_r + 2.0 * kk * inv_r2 + 2.0 * inv_r3);
      const double mixedk = fsecond_k * nxu_nyu + fprime_k_over_r * tangent_metric;
      const double phi = source.scalar_density;
      const double phin = jcp_phin_for_panel(source, j);
      ref.k2[i] += source.area * (eps_ratio * source_dk - source_d0) * phi;
      ref.k3[i] += source.area * (target_d0 - inv_eps_ratio * target_dk) * phin;
      ref.k4[i] += source.area * (mixedk - mixed0) * phi;
    }
  }
  return ref;
}

AlignedDoubleVector compute_direct_jcp_k1_reference(
    const std::vector<PrototypePanel> &panels, double kappa) {
  constexpr double inv_four_pi = 0.079577471545947667884;
  AlignedDoubleVector out(panels.size(), 0.0);
  const double kk = std::max(0.0, kappa);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    double sum = 0.0;
    const PrototypePanel &target = panels[i];
    for (std::size_t j = 0; j < panels.size(); ++j) {
      if (i == j) {
        continue;
      }
      const PrototypePanel &source = panels[j];
      const double dx = target.x - source.x;
      const double dy = target.y - source.y;
      const double dz = target.z - source.z;
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (!(r2 > 1.0e-28)) {
        continue;
      }
      const double r = std::sqrt(r2);
      const double phin = jcp_phin_for_panel(source, j);
      const double g0 = inv_four_pi / r;
      const double gk = std::exp(-kk * r) * g0;
      sum += source.area * (g0 - gk) * phin;
    }
    out[i] = sum;
  }
  return out;
}

JcpDirectParityMetrics compare_jcp_direct_reference(
    const std::vector<PrototypePanel> &panels,
    const std::vector<double> &production_phi,
    const std::vector<double> &production_phin, double eps_ratio,
    double kappa) {
  if (panels.empty()) {
    throw std::invalid_argument("JCP direct parity requires at least one panel");
  }
  if (production_phi.size() != panels.size() ||
      production_phin.size() != panels.size()) {
    throw std::invalid_argument(
        "JCP direct parity production vectors do not match panel count");
  }

  const AlignedDoubleVector k1 =
      compute_direct_jcp_k1_reference(panels, kappa);
  const DerivativeReference derivative =
      compute_direct_derivative_reference(panels, eps_ratio, kappa);
  AlignedDoubleVector analytic_phi(panels.size(), 0.0);
  AlignedDoubleVector analytic_phin(panels.size(), 0.0);
  AlignedDoubleVector production_phi_aligned(production_phi.begin(),
                                             production_phi.end());
  AlignedDoubleVector production_phin_aligned(production_phin.begin(),
                                              production_phin.end());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    analytic_phi[i] = k1[i] + derivative.k2[i];
    analytic_phin[i] = derivative.k3[i] + derivative.k4[i];
  }

  PrototypeMetrics phi_error;
  PrototypeMetrics phin_error;
  compute_error(production_phi_aligned, analytic_phi, phi_error);
  compute_error(production_phin_aligned, analytic_phin, phin_error);

  double err2 = 0.0;
  double ref2 = 0.0;
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const double phi_diff = analytic_phi[i] - production_phi_aligned[i];
    const double phin_diff = analytic_phin[i] - production_phin_aligned[i];
    err2 += phi_diff * phi_diff + phin_diff * phin_diff;
    ref2 += production_phi_aligned[i] * production_phi_aligned[i] +
            production_phin_aligned[i] * production_phin_aligned[i];
  }

  JcpDirectParityMetrics metrics;
  metrics.panels = panels.size();
  metrics.phi_rel_rms_error = phi_error.rel_rms_error;
  metrics.phin_rel_rms_error = phin_error.rel_rms_error;
  metrics.aggregate_rel_rms_error =
      std::sqrt(err2 / std::max(ref2, 1.0e-30));
  metrics.max_abs_error =
      std::max(phi_error.max_abs_error, phin_error.max_abs_error);
  metrics.finite_outputs = phi_error.finite_outputs && phin_error.finite_outputs &&
                           std::isfinite(metrics.aggregate_rel_rms_error);
  return metrics;
}

JcpK1BridgeMetrics
run_jcp_k1_bridge_probe(const PrototypeConfig &input_config,
                        const std::vector<PrototypePanel> &panels) {
  if (panels.empty()) {
    throw std::invalid_argument("JCP K1 bridge requires at least one panel");
  }
  PrototypeConfig laplace = input_config;
  laplace.kernel = KernelKind::Laplace;
  PrototypeConfig yukawa = input_config;
  yukawa.kernel = KernelKind::Yukawa;
  if (!(yukawa.kappa > 0.0)) {
    yukawa.kappa = 0.2;
  }

  const PrototypePointCloud source_cloud = make_jcp_k1_source_cloud(panels);
  JcpK1BridgeMetrics metrics;
  metrics.panels = panels.size();

  const auto direct_begin = Clock::now();
  const AlignedDoubleVector direct =
      compute_direct_jcp_k1_reference(panels, yukawa.kappa);
  const auto direct_end = Clock::now();
  metrics.direct_seconds = elapsed_seconds(direct_begin, direct_end);

  PrototypeEvaluation lap_eval = evaluate_point_cloud(laplace, source_cloud);
  PrototypeEvaluation yuk_eval = evaluate_point_cloud(yukawa, source_cloud);
  metrics.setup_seconds = lap_eval.metrics.setup_seconds + yuk_eval.metrics.setup_seconds;
  metrics.apply_seconds = lap_eval.metrics.apply_seconds + yuk_eval.metrics.apply_seconds;
  metrics.near_pairs = lap_eval.metrics.near_pairs + yuk_eval.metrics.near_pairs;
  metrics.far_pairs = lap_eval.metrics.far_pairs + yuk_eval.metrics.far_pairs;
  metrics.finite_outputs =
      lap_eval.metrics.finite_outputs && yuk_eval.metrics.finite_outputs;

  AlignedDoubleVector approx(direct.size(), 0.0);
  for (std::size_t i = 0; i < direct.size(); ++i) {
    approx[i] = lap_eval.approx[i] - yuk_eval.approx[i];
  }
  PrototypeMetrics err_metrics;
  compute_error(direct, approx, err_metrics);
  metrics.rel_rms_error = err_metrics.rel_rms_error;
  metrics.rms_abs_error = err_metrics.rms_abs_error;
  metrics.max_abs_error = err_metrics.max_abs_error;
  metrics.finite_outputs = metrics.finite_outputs && err_metrics.finite_outputs;
  return metrics;
}

std::string classify_scope_e_decision(const JcpK1BridgeMetrics &metrics) {
  if (!metrics.finite_outputs || !std::isfinite(metrics.rel_rms_error)) {
    return "accuracy_fail";
  }
  if (metrics.rel_rms_error > 5.0e-2) {
    return "accuracy_fail";
  }
  if (metrics.far_pairs == 0) {
    return "no_far_work";
  }
  if (metrics.setup_seconds > metrics.direct_seconds) {
    return "setup_limited";
  }
  if (!(metrics.apply_seconds < metrics.direct_seconds)) {
    return "no_speedup";
  }
  return "pass";
}

JcpDerivativeBridgeMetrics
run_jcp_derivative_bridge_probe(const PrototypeConfig &input_config,
                                const std::vector<PrototypePanel> &panels) {
  if (panels.empty()) {
    throw std::invalid_argument("JCP derivative bridge requires at least one panel");
  }
  PrototypeConfig laplace = input_config;
  laplace.kernel = KernelKind::Laplace;
  PrototypeConfig yukawa = input_config;
  yukawa.kernel = KernelKind::Yukawa;
  if (!(yukawa.kappa > 0.0)) {
    yukawa.kappa = 0.2;
  }
  constexpr double eps_ratio = 80.0;
  constexpr double inv_eps_ratio = 1.0 / eps_ratio;
  constexpr double h = 1.0e-3;

  JcpDerivativeBridgeMetrics metrics;
  metrics.panels = panels.size();
  metrics.channels = 3;

  const auto direct_begin = Clock::now();
  const DerivativeReference direct =
      compute_direct_derivative_reference(panels, eps_ratio, yukawa.kappa);
  const auto direct_end = Clock::now();
  metrics.direct_seconds = elapsed_seconds(direct_begin, direct_end);

  const PrototypePointCloud target_base = make_panel_target_cloud(panels, 0.0);
  const PrototypePointCloud target_plus = make_panel_target_cloud(panels, h);
  const PrototypePointCloud target_minus = make_panel_target_cloud(panels, -h);
  const PrototypePointCloud phi_source_dipole =
      make_source_normal_dipole_cloud(panels, true, h);
  const PrototypePointCloud phin_source = make_jcp_k1_source_cloud(panels);

  auto add_eval_metrics = [&](const PrototypeEvaluation &eval) {
    metrics.setup_seconds += eval.metrics.setup_seconds;
    metrics.apply_seconds += eval.metrics.apply_seconds;
    metrics.near_pairs += eval.metrics.near_pairs;
    metrics.far_pairs += eval.metrics.far_pairs;
    metrics.finite_outputs = metrics.finite_outputs && eval.metrics.finite_outputs;
  };

  const PrototypeEvaluation k2_lap =
      evaluate_source_target_cloud(laplace, phi_source_dipole, target_base, true);
  const PrototypeEvaluation k2_yuk =
      evaluate_source_target_cloud(yukawa, phi_source_dipole, target_base, true);
  add_eval_metrics(k2_lap);
  add_eval_metrics(k2_yuk);

  const PrototypeEvaluation k3_lap_plus =
      evaluate_source_target_cloud(laplace, phin_source, target_plus, true);
  const PrototypeEvaluation k3_lap_minus =
      evaluate_source_target_cloud(laplace, phin_source, target_minus, true);
  const PrototypeEvaluation k3_yuk_plus =
      evaluate_source_target_cloud(yukawa, phin_source, target_plus, true);
  const PrototypeEvaluation k3_yuk_minus =
      evaluate_source_target_cloud(yukawa, phin_source, target_minus, true);
  add_eval_metrics(k3_lap_plus);
  add_eval_metrics(k3_lap_minus);
  add_eval_metrics(k3_yuk_plus);
  add_eval_metrics(k3_yuk_minus);

  const PrototypeEvaluation k4_lap_plus =
      evaluate_source_target_cloud(laplace, phi_source_dipole, target_plus, true);
  const PrototypeEvaluation k4_lap_minus =
      evaluate_source_target_cloud(laplace, phi_source_dipole, target_minus, true);
  const PrototypeEvaluation k4_yuk_plus =
      evaluate_source_target_cloud(yukawa, phi_source_dipole, target_plus, true);
  const PrototypeEvaluation k4_yuk_minus =
      evaluate_source_target_cloud(yukawa, phi_source_dipole, target_minus, true);
  add_eval_metrics(k4_lap_plus);
  add_eval_metrics(k4_lap_minus);
  add_eval_metrics(k4_yuk_plus);
  add_eval_metrics(k4_yuk_minus);

  AlignedDoubleVector approx_k2(panels.size(), 0.0);
  AlignedDoubleVector approx_k3(panels.size(), 0.0);
  AlignedDoubleVector approx_k4(panels.size(), 0.0);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    approx_k2[i] = eps_ratio * k2_yuk.approx[i] - k2_lap.approx[i];
    const double lap_target_deriv =
        (k3_lap_plus.approx[i] - k3_lap_minus.approx[i]) / (2.0 * h);
    const double yuk_target_deriv =
        (k3_yuk_plus.approx[i] - k3_yuk_minus.approx[i]) / (2.0 * h);
    approx_k3[i] = lap_target_deriv - inv_eps_ratio * yuk_target_deriv;
    const double lap_mixed =
        (k4_lap_plus.approx[i] - k4_lap_minus.approx[i]) / (2.0 * h);
    const double yuk_mixed =
        (k4_yuk_plus.approx[i] - k4_yuk_minus.approx[i]) / (2.0 * h);
    approx_k4[i] = lap_mixed - yuk_mixed;
  }

  PrototypeMetrics e2, e3, e4;
  compute_error(direct.k2, approx_k2, e2);
  compute_error(direct.k3, approx_k3, e3);
  compute_error(direct.k4, approx_k4, e4);
  metrics.k2_rel_rms_error = e2.rel_rms_error;
  metrics.k3_rel_rms_error = e3.rel_rms_error;
  metrics.k4_rel_rms_error = e4.rel_rms_error;
  auto vector_rms = [](const AlignedDoubleVector &values) {
    double sum2 = 0.0;
    for (double value : values) {
      sum2 += value * value;
    }
    return values.empty()
               ? 0.0
               : std::sqrt(sum2 / static_cast<double>(values.size()));
  };
  metrics.k2_reference_rms = vector_rms(direct.k2);
  metrics.k3_reference_rms = vector_rms(direct.k3);
  metrics.k4_reference_rms = vector_rms(direct.k4);
  metrics.max_channel_rel_rms_error =
      std::max({e2.rel_rms_error, e3.rel_rms_error, e4.rel_rms_error});
  metrics.max_abs_error = std::max({e2.max_abs_error, e3.max_abs_error, e4.max_abs_error});
  metrics.finite_outputs =
      metrics.finite_outputs && e2.finite_outputs && e3.finite_outputs && e4.finite_outputs;
  const double err2 = e2.rms_abs_error * e2.rms_abs_error * panels.size() +
                      e3.rms_abs_error * e3.rms_abs_error * panels.size() +
                      e4.rms_abs_error * e4.rms_abs_error * panels.size();
  const double ref2 =
      (e2.rms_abs_error * e2.rms_abs_error * panels.size()) /
          std::max(e2.rel_rms_error * e2.rel_rms_error, 1.0e-30) +
      (e3.rms_abs_error * e3.rms_abs_error * panels.size()) /
          std::max(e3.rel_rms_error * e3.rel_rms_error, 1.0e-30) +
      (e4.rms_abs_error * e4.rms_abs_error * panels.size()) /
          std::max(e4.rel_rms_error * e4.rel_rms_error, 1.0e-30);
  metrics.aggregate_rel_rms_error = std::sqrt(err2 / std::max(ref2, 1.0e-30));
  return metrics;
}

std::string
classify_scope_f_decision(const JcpDerivativeBridgeMetrics &metrics) {
  if (!metrics.finite_outputs || !std::isfinite(metrics.aggregate_rel_rms_error)) {
    return "accuracy_fail";
  }
  if (metrics.aggregate_rel_rms_error > 1.0e-1 ||
      metrics.max_channel_rel_rms_error > 2.0e-1) {
    return "accuracy_fail";
  }
  if (metrics.far_pairs == 0) {
    return "no_far_work";
  }
  if (metrics.setup_seconds > metrics.direct_seconds) {
    return "setup_limited";
  }
  if (!(metrics.apply_seconds < metrics.direct_seconds)) {
    return "no_speedup";
  }
  return "pass";
}

JcpOfflineOperatorResult evaluate_jcp_offline_operator(
    const PrototypeConfig &input_config, const std::vector<PrototypePanel> &panels,
    double eps_ratio, double kappa) {
  if (panels.empty()) {
    throw std::invalid_argument("JCP offline operator requires at least one panel");
  }
  PrototypeConfig laplace = input_config;
  laplace.kernel = KernelKind::Laplace;
  PrototypeConfig yukawa = input_config;
  yukawa.kernel = KernelKind::Yukawa;
  yukawa.kappa = kappa > 0.0 ? kappa : input_config.kappa;
  if (!(yukawa.kappa > 0.0)) {
    yukawa.kappa = 0.2;
  }
  constexpr double h = 1.0e-3;
  const double safe_eps_ratio = std::max(1.0e-12, eps_ratio);
  const double inv_eps_ratio = 1.0 / safe_eps_ratio;

  JcpOfflineOperatorResult result;
  result.panels = panels.size();
  result.phi.assign(panels.size(), 0.0);
  result.phin.assign(panels.size(), 0.0);

  const PrototypePointCloud target_base = make_panel_target_cloud(panels, 0.0);
  const PrototypePointCloud target_plus = make_panel_target_cloud(panels, h);
  const PrototypePointCloud target_minus = make_panel_target_cloud(panels, -h);
  const PrototypePointCloud phin_source = make_jcp_k1_source_cloud(panels);
  const PrototypePointCloud phi_source_dipole =
      make_source_normal_dipole_cloud(panels, true, h);

  auto add_eval_metrics = [&](const PrototypeEvaluation &eval) {
    result.setup_seconds += eval.metrics.setup_seconds;
    result.apply_seconds += eval.metrics.apply_seconds;
    result.near_pairs += eval.metrics.near_pairs;
    result.far_pairs += eval.metrics.far_pairs;
    result.finite_outputs =
        result.finite_outputs && eval.metrics.finite_outputs;
  };

  const PrototypeEvaluation k1_lap =
      evaluate_source_target_cloud(laplace, phin_source, phin_source, false, false);
  const PrototypeEvaluation k1_yuk =
      evaluate_source_target_cloud(yukawa, phin_source, phin_source, false, false);
  add_eval_metrics(k1_lap);
  add_eval_metrics(k1_yuk);

  const PrototypeEvaluation k2_lap =
      evaluate_source_target_cloud(laplace, phi_source_dipole, target_base, true,
                                   false);
  const PrototypeEvaluation k2_yuk =
      evaluate_source_target_cloud(yukawa, phi_source_dipole, target_base, true,
                                   false);
  add_eval_metrics(k2_lap);
  add_eval_metrics(k2_yuk);

  const PrototypeEvaluation k3_lap_plus =
      evaluate_source_target_cloud(laplace, phin_source, target_plus, true, false);
  const PrototypeEvaluation k3_lap_minus =
      evaluate_source_target_cloud(laplace, phin_source, target_minus, true, false);
  const PrototypeEvaluation k3_yuk_plus =
      evaluate_source_target_cloud(yukawa, phin_source, target_plus, true, false);
  const PrototypeEvaluation k3_yuk_minus =
      evaluate_source_target_cloud(yukawa, phin_source, target_minus, true, false);
  add_eval_metrics(k3_lap_plus);
  add_eval_metrics(k3_lap_minus);
  add_eval_metrics(k3_yuk_plus);
  add_eval_metrics(k3_yuk_minus);

  const PrototypeEvaluation k4_lap_plus =
      evaluate_source_target_cloud(laplace, phi_source_dipole, target_plus, true,
                                   false);
  const PrototypeEvaluation k4_lap_minus =
      evaluate_source_target_cloud(laplace, phi_source_dipole, target_minus, true,
                                   false);
  const PrototypeEvaluation k4_yuk_plus =
      evaluate_source_target_cloud(yukawa, phi_source_dipole, target_plus, true,
                                   false);
  const PrototypeEvaluation k4_yuk_minus =
      evaluate_source_target_cloud(yukawa, phi_source_dipole, target_minus, true,
                                   false);
  add_eval_metrics(k4_lap_plus);
  add_eval_metrics(k4_lap_minus);
  add_eval_metrics(k4_yuk_plus);
  add_eval_metrics(k4_yuk_minus);

  for (std::size_t i = 0; i < panels.size(); ++i) {
    const double k1 = k1_lap.approx[i] - k1_yuk.approx[i];
    const double k2 = safe_eps_ratio * k2_yuk.approx[i] - k2_lap.approx[i];
    const double lap_target_deriv =
        (k3_lap_plus.approx[i] - k3_lap_minus.approx[i]) / (2.0 * h);
    const double yuk_target_deriv =
        (k3_yuk_plus.approx[i] - k3_yuk_minus.approx[i]) / (2.0 * h);
    const double k3 = lap_target_deriv - inv_eps_ratio * yuk_target_deriv;
    const double lap_mixed =
        (k4_lap_plus.approx[i] - k4_lap_minus.approx[i]) / (2.0 * h);
    const double yuk_mixed =
        (k4_yuk_plus.approx[i] - k4_yuk_minus.approx[i]) / (2.0 * h);
    const double k4 = lap_mixed - yuk_mixed;
    result.phi[i] = k1 + k2;
    result.phin[i] = k3 + k4;
    result.finite_outputs =
        result.finite_outputs && std::isfinite(result.phi[i]) &&
        std::isfinite(result.phin[i]);
  }
  return result;
}

JcpAnalyticOperatorResult evaluate_jcp_offline_operator_analytic(
    const PrototypeConfig &input_config, const std::vector<PrototypePanel> &panels,
    double eps_ratio, double kappa) {
  std::vector<double> phi_values(panels.size(), 0.0);
  std::vector<double> phin_values(panels.size(), 0.0);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    phi_values[i] = panels[i].scalar_density;
    phin_values[i] = jcp_phin_for_panel(panels[i], i);
  }
  return evaluate_jcp_offline_operator_analytic(input_config, panels, phi_values,
                                                phin_values, eps_ratio, kappa);
}

JcpAnalyticOperatorResult evaluate_jcp_offline_operator_analytic(
    const PrototypeConfig &input_config, const std::vector<PrototypePanel> &panels,
    const std::vector<double> &phi_values,
    const std::vector<double> &phin_values, double eps_ratio, double kappa) {
  PrototypeConfig config = normalize_jcp_analytic_config(input_config);
  TitanExecutor executor(static_cast<std::size_t>(std::max(1, config.threads)));
  return evaluate_jcp_offline_operator_analytic(config, panels, phi_values,
                                                phin_values, eps_ratio, kappa,
                                                executor);
}

JcpAnalyticOperatorResult evaluate_jcp_offline_operator_analytic(
    const PrototypeConfig &input_config, const std::vector<PrototypePanel> &panels,
    const std::vector<double> &phi_values,
    const std::vector<double> &phin_values, double eps_ratio, double kappa,
    TitanExecutor &executor) {
  if (panels.empty()) {
    throw std::invalid_argument(
        "JCP analytic offline operator requires at least one panel");
  }
  if (phi_values.size() != panels.size() || phin_values.size() != panels.size()) {
    throw std::invalid_argument(
        "JCP analytic offline operator density vectors must match panels");
  }
  FmmM2lAnalyticContext context(input_config, panels, eps_ratio, kappa);
  return context.apply(phi_values, phin_values, executor);
}

namespace {

constexpr double kScopeKGmresLhsPhi = 10.0;
constexpr double kScopeKGmresLhsPhin = 10.0;

double l2_norm(const std::vector<double> &values) {
  double sum = 0.0;
  for (double value : values) {
    sum += value * value;
  }
  return std::sqrt(sum);
}

std::vector<double> subtract_vectors(const std::vector<double> &a,
                                     const std::vector<double> &b) {
  if (a.size() != b.size()) {
    throw std::invalid_argument("Scope K vector size mismatch");
  }
  std::vector<double> out(a.size(), 0.0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    out[i] = a[i] - b[i];
  }
  return out;
}

double relative_norm(const std::vector<double> &residual,
                     const std::vector<double> &rhs) {
  return l2_norm(residual) / std::max(l2_norm(rhs), 1.0e-30);
}

bool finite_vector(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); });
}

std::vector<double> make_scope_k_truth_vector(std::size_t panels) {
  std::vector<double> x(2 * panels, 0.0);
  for (std::size_t i = 0; i < panels; ++i) {
    const double t = static_cast<double>(i + 1);
    x[i] = 0.7 * std::cos(0.013 * t) + 0.3 * std::sin(0.029 * t);
    x[panels + i] = 0.65 * std::sin(0.017 * t) + 0.2 * std::cos(0.011 * t);
  }
  return x;
}

std::vector<double> direct_jcp_scope_k_matvec(
    const std::vector<PrototypePanel> &panels, const std::vector<double> &x,
    double eps_ratio, double kappa) {
  const std::size_t n = panels.size();
  if (x.size() != 2 * n) {
    throw std::invalid_argument("Scope K direct matvec requires a 2N vector");
  }
  std::vector<double> out(2 * n, 0.0);
  for (std::size_t ti = 0; ti < n; ++ti) {
    double phi_sum = 0.0;
    double phin_sum = 0.0;
    for (std::size_t si = 0; si < n; ++si) {
      if (ti == si) {
        continue;
      }
      PrototypePanel source = panels[si];
      source.scalar_density = x[si];
      add_jcp_direct_pair(panels[ti], source, x[n + si], eps_ratio, kappa,
                          phi_sum, phin_sum);
    }
    out[ti] = kScopeKGmresLhsPhi * x[ti] - phi_sum;
    out[n + ti] = kScopeKGmresLhsPhin * x[n + ti] - phin_sum;
  }
  return out;
}

struct AdapterMatvecOutcome {
  std::vector<double> out;
  std::size_t near_pairs = 0;
  std::size_t far_pairs = 0;
  double setup_ms = 0.0;
  double apply_ms = 0.0;
  bool finite_outputs = true;
};

AdapterMatvecOutcome adapter_jcp_scope_k_matvec(
    const PrototypeConfig &config, const std::vector<PrototypePanel> &panels,
    const std::vector<double> &x, double eps_ratio, double kappa) {
  const std::size_t n = panels.size();
  if (x.size() != 2 * n) {
    throw std::invalid_argument("Scope K adapter matvec requires a 2N vector");
  }
  std::vector<double> phi(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(n));
  std::vector<double> phin(x.begin() + static_cast<std::ptrdiff_t>(n), x.end());
  const JcpAnalyticOperatorResult result =
      evaluate_jcp_offline_operator_analytic(config, panels, phi, phin,
                                             eps_ratio, kappa);
  AdapterMatvecOutcome out;
  out.out.assign(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    out.out[i] = kScopeKGmresLhsPhi * x[i] - result.phi[i];
    out.out[n + i] = kScopeKGmresLhsPhin * x[n + i] - result.phin[i];
  }
  out.near_pairs = result.near_pairs;
  out.far_pairs = result.far_pairs;
  out.setup_ms = result.setup_seconds * 1000.0;
  out.apply_ms = result.apply_seconds * 1000.0;
  out.finite_outputs = result.finite_outputs && finite_vector(out.out);
  return out;
}

std::vector<double> solve_dense_system(std::vector<std::vector<double>> a,
                                       std::vector<double> b) {
  const std::size_t n = b.size();
  for (std::size_t col = 0; col < n; ++col) {
    std::size_t pivot = col;
    double best = std::abs(a[col][col]);
    for (std::size_t row = col + 1; row < n; ++row) {
      const double value = std::abs(a[row][col]);
      if (value > best) {
        best = value;
        pivot = row;
      }
    }
    if (!(best > 1.0e-30)) {
      a[col][col] += 1.0e-12;
      best = std::abs(a[col][col]);
    }
    if (pivot != col) {
      std::swap(a[pivot], a[col]);
      std::swap(b[pivot], b[col]);
    }
    const double diag = a[col][col];
    for (std::size_t row = col + 1; row < n; ++row) {
      const double factor = a[row][col] / diag;
      if (factor == 0.0) {
        continue;
      }
      for (std::size_t k = col; k < n; ++k) {
        a[row][k] -= factor * a[col][k];
      }
      b[row] -= factor * b[col];
    }
  }
  std::vector<double> x(n, 0.0);
  for (std::size_t ii = 0; ii < n; ++ii) {
    const std::size_t i = n - 1 - ii;
    double rhs = b[i];
    for (std::size_t j = i + 1; j < n; ++j) {
      rhs -= a[i][j] * x[j];
    }
    x[i] = rhs / a[i][i];
  }
  return x;
}

struct ScopeKGmresResult {
  std::vector<double> x;
  int iterations = 0;
  double initial_residual = 0.0;
  double final_residual = 0.0;
  double matvec_ms = 0.0;
  bool converged = false;
  bool finite_outputs = true;
};

ScopeKGmresResult run_scope_k_gmres(
    const std::function<std::vector<double>(const std::vector<double> &)> &matvec,
    const std::vector<double> &rhs, int restart, int max_iterations,
    double tolerance) {
  if (rhs.empty()) {
    throw std::invalid_argument("Scope K GMRES requires a non-empty RHS");
  }
  restart = std::max(1, restart);
  max_iterations = std::max(1, max_iterations);
  tolerance = tolerance > 0.0 ? tolerance : 1.0e-6;
  const std::size_t dim = rhs.size();
  ScopeKGmresResult result;
  result.x.assign(dim, 0.0);
  result.initial_residual = 1.0;
  result.final_residual = 1.0;
  const double rhs_norm = std::max(l2_norm(rhs), 1.0e-30);

  while (result.iterations < max_iterations && !result.converged) {
    const std::vector<double> ax = matvec(result.x);
    std::vector<double> residual = subtract_vectors(rhs, ax);
    const double beta = l2_norm(residual);
    const double relative = beta / rhs_norm;
    if (result.iterations == 0) {
      result.initial_residual = relative;
    }
    result.final_residual = relative;
    if (relative <= tolerance) {
      result.converged = true;
      break;
    }

    const int cycle = std::min(restart, max_iterations - result.iterations);
    std::vector<std::vector<double>> v(static_cast<std::size_t>(cycle) + 1,
                                       std::vector<double>(dim, 0.0));
    for (std::size_t i = 0; i < dim; ++i) {
      v[0][i] = residual[i] / beta;
    }
    std::vector<std::vector<double>> h(static_cast<std::size_t>(cycle) + 1,
                                       std::vector<double>(
                                           static_cast<std::size_t>(cycle),
                                           0.0));
    std::vector<double> best_x = result.x;
    double best_residual = relative;
    int best_j = -1;

    for (int j = 0; j < cycle && result.iterations < max_iterations; ++j) {
      const auto mv_begin = Clock::now();
      std::vector<double> w = matvec(v[static_cast<std::size_t>(j)]);
      const auto mv_end = Clock::now();
      result.matvec_ms += elapsed_seconds(mv_begin, mv_end) * 1000.0;
      result.finite_outputs = result.finite_outputs && finite_vector(w);

      for (int i = 0; i <= j; ++i) {
        double dot = 0.0;
        for (std::size_t k = 0; k < dim; ++k) {
          dot += w[k] * v[static_cast<std::size_t>(i)][k];
        }
        h[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = dot;
        for (std::size_t k = 0; k < dim; ++k) {
          w[k] -= dot * v[static_cast<std::size_t>(i)][k];
        }
      }
      const double h_next = l2_norm(w);
      h[static_cast<std::size_t>(j + 1)][static_cast<std::size_t>(j)] =
          h_next;
      if (h_next > 1.0e-30 && j + 1 < cycle + 1) {
        for (std::size_t k = 0; k < dim; ++k) {
          v[static_cast<std::size_t>(j + 1)][k] = w[k] / h_next;
        }
      }

      const std::size_t m = static_cast<std::size_t>(j + 1);
      std::vector<std::vector<double>> normal(m, std::vector<double>(m, 0.0));
      std::vector<double> normal_rhs(m, 0.0);
      for (std::size_t row = 0; row <= m; ++row) {
        const double g = row == 0 ? beta : 0.0;
        for (std::size_t col = 0; col < m; ++col) {
          normal_rhs[col] += h[row][col] * g;
          for (std::size_t col2 = 0; col2 < m; ++col2) {
            normal[col][col2] += h[row][col] * h[row][col2];
          }
        }
      }
      std::vector<double> y = solve_dense_system(std::move(normal),
                                                 std::move(normal_rhs));
      std::vector<double> candidate = result.x;
      for (std::size_t col = 0; col < m; ++col) {
        for (std::size_t k = 0; k < dim; ++k) {
          candidate[k] += y[col] * v[col][k];
        }
      }
      const auto residual_begin = Clock::now();
      const std::vector<double> candidate_ax = matvec(candidate);
      const auto residual_end = Clock::now();
      result.matvec_ms +=
          elapsed_seconds(residual_begin, residual_end) * 1000.0;
      result.finite_outputs =
          result.finite_outputs && finite_vector(candidate) &&
          finite_vector(candidate_ax);
      const double candidate_residual =
          relative_norm(subtract_vectors(rhs, candidate_ax), rhs);
      ++result.iterations;
      if (candidate_residual < best_residual) {
        best_residual = candidate_residual;
        best_x = std::move(candidate);
        best_j = j;
      }
      if (candidate_residual <= tolerance) {
        result.x = best_x;
        result.final_residual = candidate_residual;
        result.converged = true;
        break;
      }
      if (!(h_next > 1.0e-30)) {
        break;
      }
    }
    if (!result.converged) {
      if (best_j < 0) {
        break;
      }
      result.x = std::move(best_x);
      result.final_residual = best_residual;
    }
  }
  return result;
}

} // namespace

JcpGuardedGmresMetrics run_jcp_fmm_m2l_guarded_gmres_probe(
    const PrototypeConfig &input_config, const std::vector<PrototypePanel> &panels,
    double eps_ratio, double kappa, int restart, int max_iterations,
    double tolerance) {
  if (panels.empty()) {
    throw std::invalid_argument("Scope K GMRES probe requires panels");
  }
  PrototypeConfig config = input_config;
  config.kernel = KernelKind::Yukawa;
  config.kappa = kappa > 0.0 ? kappa : input_config.kappa;
  config.max_leaf_size = std::max(1, config.max_leaf_size);
  config.threads = std::max(1, config.threads);

  JcpGuardedGmresMetrics metrics;
  metrics.panels = panels.size();
  const std::vector<double> x_true = make_scope_k_truth_vector(panels.size());
  const auto direct_matvec = [&](const std::vector<double> &x) {
    return direct_jcp_scope_k_matvec(panels, x, eps_ratio, config.kappa);
  };
  AdapterMatvecOutcome representative_adapter;
  const auto adapter_matvec = [&](const std::vector<double> &x) {
    AdapterMatvecOutcome outcome =
        adapter_jcp_scope_k_matvec(config, panels, x, eps_ratio, config.kappa);
    if (representative_adapter.out.empty()) {
      representative_adapter = outcome;
    }
    metrics.finite_outputs = metrics.finite_outputs && outcome.finite_outputs;
    return outcome.out;
  };

  const std::vector<double> rhs = direct_matvec(x_true);
  const ScopeKGmresResult direct =
      run_scope_k_gmres(direct_matvec, rhs, restart, max_iterations, tolerance);
  const ScopeKGmresResult adapter =
      run_scope_k_gmres(adapter_matvec, rhs, restart, max_iterations, tolerance);

  const std::vector<double> adapter_direct_ax = direct_matvec(adapter.x);
  const std::vector<double> solution_delta = subtract_vectors(adapter.x, x_true);
  metrics.direct_iterations = direct.iterations;
  metrics.adapter_iterations = adapter.iterations;
  metrics.direct_initial_residual = direct.initial_residual;
  metrics.direct_final_residual = direct.final_residual;
  metrics.adapter_initial_residual = adapter.initial_residual;
  metrics.adapter_final_residual = adapter.final_residual;
  metrics.direct_true_residual =
      relative_norm(subtract_vectors(rhs, adapter_direct_ax), rhs);
  metrics.solution_rel_error =
      l2_norm(solution_delta) / std::max(l2_norm(x_true), 1.0e-30);
  metrics.direct_matvec_ms =
      direct.iterations > 0 ? direct.matvec_ms / static_cast<double>(direct.iterations)
                            : direct.matvec_ms;
  metrics.adapter_matvec_ms =
      adapter.iterations > 0
          ? adapter.matvec_ms / static_cast<double>(adapter.iterations)
          : adapter.matvec_ms;
  metrics.speedup_vs_direct_matvec =
      metrics.adapter_matvec_ms > 0.0
          ? metrics.direct_matvec_ms / metrics.adapter_matvec_ms
          : 0.0;
  metrics.direct_converged = direct.converged;
  metrics.adapter_converged = adapter.converged;
  metrics.adapter_residual_decreased =
      adapter.final_residual < 0.99 * std::max(adapter.initial_residual, 1.0e-30);
  metrics.finite_outputs =
      metrics.finite_outputs && direct.finite_outputs && adapter.finite_outputs &&
      finite_vector(rhs) && finite_vector(adapter_direct_ax);
  if (!representative_adapter.out.empty()) {
    metrics.near_pairs = representative_adapter.near_pairs;
    metrics.far_pairs = representative_adapter.far_pairs;
    metrics.adapter_setup_ms = representative_adapter.setup_ms;
    metrics.adapter_apply_ms = representative_adapter.apply_ms;
  }
  return metrics;
}

JcpChannelPrototypeMetrics
run_jcp_channel_probe(const PrototypeConfig &config,
                      const std::vector<PrototypePanel> &panels) {
  const auto channels = make_jcp_channel_point_clouds(panels);
  JcpChannelPrototypeMetrics out;
  out.panels = panels.size();
  out.channels = channels.size();

  double weighted_err2 = 0.0;
  double weighted_ref2 = 0.0;
  for (const PrototypePointCloud &cloud : channels) {
    const PrototypeMetrics m = run_point_cloud_probe(config, cloud);
    out.near_pairs += m.near_pairs;
    out.far_pairs += m.far_pairs;
    out.setup_seconds += m.setup_seconds;
    out.apply_seconds += m.apply_seconds;
    out.direct_seconds += m.direct_seconds;
    out.max_channel_rel_rms_error =
        std::max(out.max_channel_rel_rms_error, m.rel_rms_error);
    out.max_abs_error = std::max(out.max_abs_error, m.max_abs_error);
    out.finite_outputs = out.finite_outputs && m.finite_outputs;

    const double samples = static_cast<double>(std::max<std::size_t>(1, m.targets));
    const double err2 = m.rms_abs_error * m.rms_abs_error * samples;
    const double ref2 = err2 / std::max(m.rel_rms_error * m.rel_rms_error, 1.0e-30);
    weighted_err2 += err2;
    weighted_ref2 += ref2;
  }
  out.aggregate_rel_rms_error =
      std::sqrt(weighted_err2 / std::max(weighted_ref2, 1.0e-30));
  return out;
}

std::string classify_scope_d_decision(const JcpChannelPrototypeMetrics &metrics) {
  if (!metrics.finite_outputs || !std::isfinite(metrics.aggregate_rel_rms_error)) {
    return "accuracy_fail";
  }
  if (metrics.aggregate_rel_rms_error > 5.0e-2 ||
      metrics.max_channel_rel_rms_error > 1.0e-1) {
    return "accuracy_fail";
  }
  if (metrics.far_pairs == 0) {
    return "no_far_work";
  }
  if (metrics.setup_seconds > metrics.direct_seconds) {
    return "setup_limited";
  }
  if (!(metrics.apply_seconds < metrics.direct_seconds)) {
    return "no_speedup";
  }
  return "pass";
}

std::string classify_scope_a_decision(const PrototypeMetrics &metrics,
                                      const PrototypeConfig &config) {
  if (!metrics.finite_outputs || !std::isfinite(metrics.rel_rms_error)) {
    return "accuracy_fail";
  }
  const double rel_gate =
      config.kernel == KernelKind::Yukawa ? 3.0e-3 : 1.0e-3;
  if (metrics.rel_rms_error > rel_gate) {
    return "accuracy_fail";
  }
  if (metrics.setup_seconds > metrics.direct_seconds) {
    return "setup_limited";
  }
  if (!(metrics.apply_seconds < metrics.direct_seconds)) {
    return "no_speedup";
  }
  if (metrics.apply_seconds <= 0.5 * metrics.direct_seconds) {
    return "promote_to_scope_b";
  }
  return "promising_scalar";
}

} // namespace titan::tabi::fmm_m2l

namespace TitanPBE::Tabi {

std::shared_ptr<titan::tabi::fmm_m2l::FmmM2lAnalyticContext>
ensure_jcp_fmm_m2l_context(
    const std::vector<Panel> &panels, double eps_ratio, double kappa,
    std::size_t leaf_size, double opening_theta, TitanExecutor &executor,
    std::shared_ptr<titan::tabi::fmm_m2l::FmmM2lAnalyticContext>
        *context_cache) {
  titan::tabi::fmm_m2l::PrototypeConfig config;
  config.kernel = titan::tabi::fmm_m2l::KernelKind::Yukawa;
  config.kappa = kappa;
  config.max_leaf_size =
      static_cast<int>(std::clamp<std::size_t>(leaf_size == 0 ? 64 : leaf_size,
                                               8, 512));
  config.expansion_order = 4;
  config.separation_theta =
      std::clamp(std::isfinite(opening_theta) && opening_theta > 0.0
                     ? opening_theta
                     : 0.8,
                 0.05, 2.0);
  config.threads =
      static_cast<int>(std::max<std::size_t>(1, executor.capacity()));

  std::shared_ptr<titan::tabi::fmm_m2l::FmmM2lAnalyticContext> local_context;
  std::shared_ptr<titan::tabi::fmm_m2l::FmmM2lAnalyticContext> &context =
      context_cache != nullptr ? *context_cache : local_context;
  if (context == nullptr ||
      !context->matches(config, panels.size(), eps_ratio, kappa)) {
    std::vector<titan::tabi::fmm_m2l::PrototypePanel> prototype;
    prototype.reserve(panels.size());
    for (const Panel &panel : panels) {
      titan::tabi::fmm_m2l::PrototypePanel converted;
      converted.x = panel.x;
      converted.y = panel.y;
      converted.z = panel.z;
      converted.nx = panel.nx;
      converted.ny = panel.ny;
      converted.nz = panel.nz;
      converted.area = panel.area;
      prototype.push_back(converted);
    }
    context = std::make_shared<titan::tabi::fmm_m2l::FmmM2lAnalyticContext>(
        config, std::move(prototype), eps_ratio, kappa);
  }
  return context;
}

void apply_jcp_fmm_m2l_single_matvec(
    const std::vector<Panel> &panels, const std::vector<double> &unknowns,
    double eps_ratio, double kappa, double lhs_phi, double lhs_phin,
    std::size_t leaf_size, double opening_theta, TitanExecutor &executor,
    std::vector<double> &out, FmmM2lSingleMatvecStats *stats, bool coupled,
    std::shared_ptr<titan::tabi::fmm_m2l::FmmM2lAnalyticContext>
        *context_cache) {
  const std::size_t n = panels.size();
  if (n == 0) {
    out.clear();
    if (stats != nullptr) {
      *stats = FmmM2lSingleMatvecStats{};
    }
    return;
  }
  if (unknowns.size() < (coupled ? 2 * n : n)) {
    throw std::invalid_argument(
        "PBE TABI FMM/M2L single-matvec received an incomplete JCP unknown vector");
  }

  std::vector<double> phi(n, 0.0);
  std::vector<double> phin(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    phi[i] = unknowns[i];
    phin[i] = coupled ? unknowns[n + i] : 0.0;
  }

  const std::shared_ptr<titan::tabi::fmm_m2l::FmmM2lAnalyticContext> context =
      ensure_jcp_fmm_m2l_context(panels, eps_ratio, kappa, leaf_size,
                                 opening_theta, executor, context_cache);

  const titan::tabi::fmm_m2l::JcpAnalyticOperatorResult result =
      context->apply(phi, phin, executor);

  out.assign(coupled ? 2 * n : n, 0.0);
  bool finite = result.finite_outputs;
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = lhs_phi * unknowns[i] - result.phi[i];
    finite = finite && std::isfinite(out[i]);
    if (coupled) {
      out[n + i] = lhs_phin * unknowns[n + i] - result.phin[i];
      finite = finite && std::isfinite(out[n + i]);
    }
  }

  if (stats != nullptr) {
    stats->near_pairs = result.near_pairs;
    stats->far_pairs = result.far_pairs;
    stats->setup_ms = 1000.0 * result.setup_seconds;
    stats->apply_ms = 1000.0 * result.apply_seconds;
    stats->cached = context_cache != nullptr;
    stats->context_nodes = result.context_nodes;
    stats->context_work_items = result.context_work_items;
    stats->finite_outputs = finite;
  }
}

} // namespace TitanPBE::Tabi
