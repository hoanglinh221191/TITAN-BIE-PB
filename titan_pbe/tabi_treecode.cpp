#include "tabi_types.h"

#include "../titan_executor.h"
#include "../titan_simd_dispatch.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace TitanPBE::Tabi {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr std::size_t kLeafSize = 32;
constexpr std::size_t kBranchFactor = 8;
constexpr double kOpeningTheta = 0.45;
constexpr std::uint32_t kMortonMax = (1u << 21) - 1u;
constexpr int kTaylorOrderDefault = 3;
constexpr int kTaylorOrderMax = 4;
constexpr std::size_t kTaylorMomentMax = 35;

struct MomentIndex {
  int ax{0};
  int ay{0};
  int az{0};
  double inv_factorial{1.0};
};

struct NodeMoments {
  std::array<double, kTaylorMomentMax> values{};
};

std::size_t effective_workers(const TitanExecutor &executor,
                              std::size_t work_items) noexcept {
  if (work_items == 0)
    return 1;
  return std::max<std::size_t>(
      1, std::min<std::size_t>(executor.capacity(), work_items));
}

double elapsed_ms(std::chrono::steady_clock::time_point begin,
                  std::chrono::steady_clock::time_point end) noexcept {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

double factorial(int value) noexcept {
  double out = 1.0;
  for (int i = 2; i <= value; ++i)
    out *= static_cast<double>(i);
  return out;
}

const std::vector<MomentIndex> &moment_indices(int order) {
  static const std::vector<MomentIndex> indices = [] {
    std::vector<MomentIndex> out;
    out.reserve(kTaylorMomentMax);
    for (int degree = 0; degree <= kTaylorOrderMax; ++degree) {
      for (int ax = degree; ax >= 0; --ax) {
        for (int ay = degree - ax; ay >= 0; --ay) {
          const int az = degree - ax - ay;
          out.push_back({ax, ay, az,
                         1.0 / (factorial(ax) * factorial(ay) *
                                factorial(az))});
        }
      }
    }
    return out;
  }();
  const int clamped = std::clamp(order, 0, kTaylorOrderMax);
  static thread_local std::vector<MomentIndex> active;
  active.clear();
  for (const MomentIndex &index : indices) {
    if (index.ax + index.ay + index.az <= clamped)
      active.push_back(index);
  }
  return active;
}

double pow_int(double value, int power) noexcept {
  double out = 1.0;
  for (int i = 0; i < power; ++i)
    out *= value;
  return out;
}

double binomial(int n, int k) noexcept {
  if (k < 0 || k > n)
    return 0.0;
  if (k == 0 || k == n)
    return 1.0;
  double out = 1.0;
  for (int i = 1; i <= k; ++i)
    out = out * static_cast<double>(n - k + i) / static_cast<double>(i);
  return out;
}

int moment_lookup_index(int ax, int ay, int az) noexcept {
  constexpr int side = kTaylorOrderMax + 1;
  return (ax * side + ay) * side + az;
}

const std::array<int, (kTaylorOrderMax + 1) * (kTaylorOrderMax + 1) *
                         (kTaylorOrderMax + 1)> &
moment_lookup() {
  static const auto lookup = [] {
    std::array<int, (kTaylorOrderMax + 1) * (kTaylorOrderMax + 1) *
                        (kTaylorOrderMax + 1)>
        out{};
    out.fill(-1);
    const auto &indices = moment_indices(kTaylorOrderMax);
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const MomentIndex &index = indices[i];
      out[static_cast<std::size_t>(
          moment_lookup_index(index.ax, index.ay, index.az))] =
          static_cast<int>(i);
    }
    return out;
  }();
  return lookup;
}

double point_kernel(const Panel &target, double sx, double sy, double sz,
                    double kappa) noexcept {
  const double dx = target.x - sx;
  const double dy = target.y - sy;
  const double dz = target.z - sz;
  const double r2 = dx * dx + dy * dy + dz * dz;
  if (r2 <= 1.0e-16)
    return 0.0;
  const double r = std::sqrt(r2);
  const double normal_projection =
      (target.nx * dx + target.ny * dy + target.nz * dz) / r;
  const double screened = std::exp(-std::max(0.0, kappa) * r);
  return screened * normal_projection /
         (4.0 * kPi * std::max(1.0e-12, r2));
}

double panel_kernel(const Panel &target, const Panel &source,
                    double kappa) noexcept {
  return source.area * point_kernel(target, source.x, source.y, source.z,
                                    kappa);
}

double finite_difference_derivative(const Panel &target, double cx, double cy,
                                    double cz, int ax, int ay, int az,
                                    double h, double kappa) noexcept {
  if (ax == 0 && ay == 0 && az == 0)
    return point_kernel(target, cx, cy, cz, kappa);
  if (ax > 0) {
    return (finite_difference_derivative(target, cx + h, cy, cz, ax - 1, ay,
                                         az, h, kappa) -
            finite_difference_derivative(target, cx - h, cy, cz, ax - 1, ay,
                                         az, h, kappa)) /
           (2.0 * h);
  }
  if (ay > 0) {
    return (finite_difference_derivative(target, cx, cy + h, cz, ax, ay - 1,
                                         az, h, kappa) -
            finite_difference_derivative(target, cx, cy - h, cz, ax, ay - 1,
                                         az, h, kappa)) /
           (2.0 * h);
  }
  return (finite_difference_derivative(target, cx, cy, cz + h, ax, ay, az - 1,
                                       h, kappa) -
          finite_difference_derivative(target, cx, cy, cz - h, ax, ay, az - 1,
                                       h, kappa)) /
         (2.0 * h);
}

double node_kernel_taylor(const Panel &target, const TreeNode &node,
                          const NodeMoments &moments, int order,
                          double kappa) noexcept {
  const double dx = target.x - node.cx;
  const double dy = target.y - node.cy;
  const double dz = target.z - node.cz;
  const double distance =
      std::sqrt(std::max(1.0e-16, dx * dx + dy * dy + dz * dz));
  const double h = std::clamp(0.05 * distance, 1.0e-4,
                              std::max(1.0e-3, 0.20 * node.radius));
  const auto &indices = moment_indices(order);
  double value = 0.0;
  for (std::size_t m = 0; m < indices.size(); ++m) {
    const MomentIndex &index = indices[m];
    const double derivative = finite_difference_derivative(
        target, node.cx, node.cy, node.cz, index.ax, index.ay, index.az, h,
        kappa);
    value += moments.values[m] * derivative * index.inv_factorial;
  }
  return std::isfinite(value) ? value : 0.0;
}

double direct_sum_scalar(const Panel &target, std::size_t target_index,
                         const std::vector<Panel> &panels,
                         const std::vector<double> &density, double kappa,
                         std::size_t begin, std::size_t end) noexcept {
  double coupled = 0.0;
  for (std::size_t j = begin; j < end; ++j) {
    if (j == target_index)
      continue;
    coupled += panel_kernel(target, panels[j], kappa) * density[j];
  }
  return coupled;
}

#if TITAN_SIMD_X86
TITAN_TARGET_AVX2
double direct_sum_avx2(const Panel &target, std::size_t target_index,
                       const std::vector<Panel> &panels,
                       const std::vector<double> &density, double kappa,
                       std::size_t begin, std::size_t end) noexcept {
  double coupled = 0.0;
  const __m256d tx = _mm256_set1_pd(target.x);
  const __m256d ty = _mm256_set1_pd(target.y);
  const __m256d tz = _mm256_set1_pd(target.z);
  const __m256d tnx = _mm256_set1_pd(target.nx);
  const __m256d tny = _mm256_set1_pd(target.ny);
  const __m256d tnz = _mm256_set1_pd(target.nz);
  const __m256d four_pi = _mm256_set1_pd(4.0 * kPi);
  const __m256d tiny_r2 = _mm256_set1_pd(1.0e-16);
  const __m256d min_denom = _mm256_set1_pd(1.0e-12);
  const double kk = std::max(0.0, kappa);

  std::size_t j = begin;
  alignas(32) double base_lane[4];
  alignas(32) double r_lane[4];
  for (; j + 4 <= end; j += 4) {
    const __m256d sx =
        _mm256_set_pd(panels[j + 3].x, panels[j + 2].x, panels[j + 1].x,
                      panels[j].x);
    const __m256d sy =
        _mm256_set_pd(panels[j + 3].y, panels[j + 2].y, panels[j + 1].y,
                      panels[j].y);
    const __m256d sz =
        _mm256_set_pd(panels[j + 3].z, panels[j + 2].z, panels[j + 1].z,
                      panels[j].z);
    const __m256d area =
        _mm256_set_pd(panels[j + 3].area, panels[j + 2].area,
                      panels[j + 1].area, panels[j].area);
    const __m256d dx = _mm256_sub_pd(tx, sx);
    const __m256d dy = _mm256_sub_pd(ty, sy);
    const __m256d dz = _mm256_sub_pd(tz, sz);
    __m256d r2 = _mm256_add_pd(_mm256_mul_pd(dx, dx),
                               _mm256_add_pd(_mm256_mul_pd(dy, dy),
                                             _mm256_mul_pd(dz, dz)));
    const __m256d valid = _mm256_cmp_pd(r2, tiny_r2, _CMP_GT_OQ);
    const __m256d r = _mm256_sqrt_pd(r2);
    const __m256d dot =
        _mm256_add_pd(_mm256_mul_pd(tnx, dx),
                      _mm256_add_pd(_mm256_mul_pd(tny, dy),
                                    _mm256_mul_pd(tnz, dz)));
    const __m256d normal = _mm256_div_pd(dot, r);
    const __m256d denom = _mm256_mul_pd(four_pi, _mm256_max_pd(r2, min_denom));
    __m256d base = _mm256_div_pd(_mm256_mul_pd(area, normal), denom);
    base = _mm256_and_pd(base, valid);
    _mm256_store_pd(base_lane, base);
    _mm256_store_pd(r_lane, r);
    for (int lane = 0; lane < 4; ++lane) {
      const std::size_t index = j + static_cast<std::size_t>(lane);
      if (index == target_index)
        continue;
      coupled += base_lane[lane] * std::exp(-kk * r_lane[lane]) *
                 density[index];
    }
  }
  for (; j < end; ++j) {
    if (j == target_index)
      continue;
    coupled += panel_kernel(target, panels[j], kappa) * density[j];
  }
  return coupled;
}
#endif

double direct_sum(const Panel &target, std::size_t target_index,
                  const std::vector<Panel> &panels,
                  const std::vector<double> &density, double kappa,
                  std::size_t begin, std::size_t end) noexcept {
  using Kernel = double (*)(const Panel &, std::size_t,
                            const std::vector<Panel> &,
                            const std::vector<double> &, double, std::size_t,
                            std::size_t) noexcept;
  static const Kernel kernel = []() noexcept -> Kernel {
#if TITAN_SIMD_X86 && TITAN_RUNTIME_TARGET_ATTRS
    return titan_simd_dispatch::runtime_has_avx2() ? &direct_sum_avx2
                                                   : &direct_sum_scalar;
#else
    return &direct_sum_scalar;
#endif
  }();
  return kernel(target, target_index, panels, density, kappa, begin, end);
}

double direct_ordered_sum(const Panel &target, std::size_t target_index,
                          const std::vector<Panel> &panels,
                          const std::vector<double> &density, double kappa,
                          const std::vector<std::size_t> &order,
                          std::size_t begin, std::size_t count,
                          std::size_t &near_count) noexcept {
  double coupled = 0.0;
  const std::size_t end = begin + count;
  for (std::size_t pos = begin; pos < end; ++pos) {
    const std::size_t j = order[pos];
    if (j == target_index)
      continue;
    coupled += panel_kernel(target, panels[j], kappa) * density[j];
    ++near_count;
  }
  return coupled;
}

std::uint64_t expand_bits_21(std::uint32_t value) noexcept {
  std::uint64_t out = 0;
  for (std::uint32_t bit = 0; bit < 21; ++bit)
    out |= static_cast<std::uint64_t>((value >> bit) & 1u) << (3u * bit);
  return out;
}

std::uint64_t morton_code(std::uint32_t x, std::uint32_t y,
                          std::uint32_t z) noexcept {
  return (expand_bits_21(x) << 2u) | (expand_bits_21(y) << 1u) |
         expand_bits_21(z);
}

struct MortonEntry {
  std::uint64_t code{0};
  std::size_t index{0};
};

std::vector<MortonEntry> radix_sort_morton(std::vector<MortonEntry> entries) {
  constexpr std::size_t kBits = 11;
  constexpr std::size_t kRadix = 1u << kBits;
  constexpr std::uint64_t kMask = kRadix - 1u;
  std::vector<MortonEntry> scratch(entries.size());
  for (std::size_t shift = 0; shift < 63; shift += kBits) {
    std::array<std::size_t, kRadix> counts{};
    for (const MortonEntry &entry : entries)
      ++counts[static_cast<std::size_t>((entry.code >> shift) & kMask)];
    std::size_t offset = 0;
    for (std::size_t bucket = 0; bucket < kRadix; ++bucket) {
      const std::size_t count = counts[bucket];
      counts[bucket] = offset;
      offset += count;
    }
    for (const MortonEntry &entry : entries) {
      const std::size_t bucket =
          static_cast<std::size_t>((entry.code >> shift) & kMask);
      scratch[counts[bucket]++] = entry;
    }
    entries.swap(scratch);
  }
  return entries;
}

TreeNode make_leaf_node(const std::vector<Panel> &panels,
                        const std::vector<std::size_t> &order,
                        std::size_t first_panel,
                        std::size_t panel_count) noexcept {
  TreeNode node;
  node.first_panel = first_panel;
  node.panel_count = panel_count;
  double weight = 0.0;
  for (std::size_t pos = first_panel; pos < first_panel + panel_count; ++pos) {
    const Panel &panel = panels[order[pos]];
    const double area = std::max(0.0, panel.area);
    node.cx += area * panel.x;
    node.cy += area * panel.y;
    node.cz += area * panel.z;
    node.area_sum += area;
    weight += area;
  }
  if (weight > 0.0) {
    node.cx /= weight;
    node.cy /= weight;
    node.cz /= weight;
  } else if (panel_count > 0) {
    for (std::size_t pos = first_panel; pos < first_panel + panel_count; ++pos) {
      const Panel &panel = panels[order[pos]];
      node.cx += panel.x;
      node.cy += panel.y;
      node.cz += panel.z;
    }
    const double inv_count = 1.0 / static_cast<double>(panel_count);
    node.cx *= inv_count;
    node.cy *= inv_count;
    node.cz *= inv_count;
  }
  for (std::size_t pos = first_panel; pos < first_panel + panel_count; ++pos) {
    const Panel &panel = panels[order[pos]];
    const double dx = panel.x - node.cx;
    const double dy = panel.y - node.cy;
    const double dz = panel.z - node.cz;
    node.radius = std::max(node.radius, std::sqrt(dx * dx + dy * dy + dz * dz));
  }
  return node;
}

TreeNode make_parent_node(const std::vector<TreeNode> &nodes,
                          std::size_t first_child,
                          std::size_t child_count) noexcept {
  TreeNode node;
  node.first_child = first_child;
  node.child_count = child_count;
  node.first_panel = nodes[first_child].first_panel;
  double weight = 0.0;
  for (std::size_t child = 0; child < child_count; ++child) {
    const TreeNode &child_node = nodes[first_child + child];
    const double area = std::max(0.0, child_node.area_sum);
    node.cx += area * child_node.cx;
    node.cy += area * child_node.cy;
    node.cz += area * child_node.cz;
    node.area_sum += area;
    node.panel_count += child_node.panel_count;
    weight += area;
  }
  if (weight > 0.0) {
    node.cx /= weight;
    node.cy /= weight;
    node.cz /= weight;
  } else if (child_count > 0) {
    node.cx = 0.0;
    node.cy = 0.0;
    node.cz = 0.0;
    for (std::size_t child = 0; child < child_count; ++child) {
      const TreeNode &child_node = nodes[first_child + child];
      node.cx += child_node.cx;
      node.cy += child_node.cy;
      node.cz += child_node.cz;
    }
    const double inv_count = 1.0 / static_cast<double>(child_count);
    node.cx *= inv_count;
    node.cy *= inv_count;
    node.cz *= inv_count;
  }
  for (std::size_t child = 0; child < child_count; ++child) {
    const TreeNode &child_node = nodes[first_child + child];
    const double dx = child_node.cx - node.cx;
    const double dy = child_node.cy - node.cy;
    const double dz = child_node.cz - node.cz;
    node.radius =
        std::max(node.radius,
                 std::sqrt(dx * dx + dy * dy + dz * dz) + child_node.radius);
  }
  return node;
}

std::vector<NodeMoments>
compute_node_moments(const std::vector<Panel> &panels,
                     const std::vector<double> &density,
                     const TreeWorkspace &tree, int taylor_order) {
  std::vector<NodeMoments> moments(tree.nodes.size());
  const auto &indices = moment_indices(taylor_order);
  const auto &lookup = moment_lookup();
  for (std::size_t node_index = 0; node_index < tree.nodes.size();
       ++node_index) {
    const TreeNode &node = tree.nodes[node_index];
    if (node.child_count == 0) {
      for (std::size_t pos = node.first_panel;
           pos < node.first_panel + node.panel_count; ++pos) {
        const std::size_t panel_index = tree.order[pos];
        const Panel &panel = panels[panel_index];
        const double weight = panel.area * density[panel_index];
        const double dx = panel.x - node.cx;
        const double dy = panel.y - node.cy;
        const double dz = panel.z - node.cz;
        for (std::size_t m = 0; m < indices.size(); ++m) {
          const MomentIndex &index = indices[m];
          moments[node_index].values[m] +=
              weight * pow_int(dx, index.ax) * pow_int(dy, index.ay) *
              pow_int(dz, index.az);
        }
      }
    } else {
      for (std::size_t child = 0; child < node.child_count; ++child) {
        const std::size_t child_index = node.first_child + child;
        const TreeNode &child_node = tree.nodes[child_index];
        const NodeMoments &child_moments = moments[child_index];
        const double ox = child_node.cx - node.cx;
        const double oy = child_node.cy - node.cy;
        const double oz = child_node.cz - node.cz;
        for (std::size_t alpha_index = 0; alpha_index < indices.size();
             ++alpha_index) {
          const MomentIndex &alpha = indices[alpha_index];
          double translated = 0.0;
          for (int bx = 0; bx <= alpha.ax; ++bx) {
            for (int by = 0; by <= alpha.ay; ++by) {
              for (int bz = 0; bz <= alpha.az; ++bz) {
                const int beta_index = lookup[static_cast<std::size_t>(
                    moment_lookup_index(bx, by, bz))];
                if (beta_index < 0)
                  continue;
                const double coeff =
                    binomial(alpha.ax, bx) * binomial(alpha.ay, by) *
                    binomial(alpha.az, bz) *
                    pow_int(ox, alpha.ax - bx) *
                    pow_int(oy, alpha.ay - by) *
                    pow_int(oz, alpha.az - bz);
                translated += coeff *
                              child_moments.values[static_cast<std::size_t>(
                                  beta_index)];
              }
            }
          }
          moments[node_index].values[alpha_index] += translated;
        }
      }
    }
  }
  return moments;
}

} // namespace

TreeWorkspace build_tree_workspace(const std::vector<Panel> &panels,
                                   TitanExecutor &executor) {
  const auto build_begin = std::chrono::steady_clock::now();
  TreeWorkspace tree;
  tree.taylor_order = kTaylorOrderDefault;
  if (panels.empty())
    return tree;

  double min_x = std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double min_z = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  double max_z = -std::numeric_limits<double>::infinity();
  for (const Panel &panel : panels) {
    min_x = std::min(min_x, panel.x);
    min_y = std::min(min_y, panel.y);
    min_z = std::min(min_z, panel.z);
    max_x = std::max(max_x, panel.x);
    max_y = std::max(max_y, panel.y);
    max_z = std::max(max_z, panel.z);
  }
  const double extent =
      std::max({max_x - min_x, max_y - min_y, max_z - min_z, 1.0e-9});
  auto quantize = [&](double value, double origin) noexcept {
    const double scaled =
        ((value - origin) / extent) * static_cast<double>(kMortonMax);
    const double clamped =
        std::clamp(scaled, 0.0, static_cast<double>(kMortonMax));
    return static_cast<std::uint32_t>(clamped);
  };

  std::vector<MortonEntry> entries(panels.size());
  const std::size_t workers = effective_workers(executor, panels.size());
  executor.parallel_for(
      panels.size(), workers, 256,
      [&](std::size_t, std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
          const Panel &panel = panels[i];
          entries[i].code =
              morton_code(quantize(panel.x, min_x), quantize(panel.y, min_y),
                          quantize(panel.z, min_z));
          entries[i].index = i;
        }
      });
  entries = radix_sort_morton(std::move(entries));

  tree.order.reserve(entries.size());
  tree.morton.reserve(entries.size());
  for (const MortonEntry &entry : entries) {
    tree.order.push_back(entry.index);
    tree.morton.push_back(entry.code);
  }

  std::vector<std::size_t> level_indices;
  const std::size_t leaf_count =
      (panels.size() + kLeafSize - 1) / kLeafSize;
  tree.nodes.reserve(leaf_count * 2 + 8);
  level_indices.reserve(leaf_count);
  for (std::size_t leaf = 0; leaf < leaf_count; ++leaf) {
    const std::size_t first_panel = leaf * kLeafSize;
    const std::size_t count =
        std::min(kLeafSize, panels.size() - first_panel);
    level_indices.push_back(tree.nodes.size());
    tree.nodes.push_back(make_leaf_node(panels, tree.order, first_panel, count));
  }
  tree.leaf_count = tree.nodes.size();

  while (level_indices.size() > 1) {
    std::vector<std::size_t> next_level;
    next_level.reserve((level_indices.size() + kBranchFactor - 1) /
                       kBranchFactor);
    for (std::size_t begin = 0; begin < level_indices.size();
         begin += kBranchFactor) {
      const std::size_t first_child = level_indices[begin];
      const std::size_t child_count =
          std::min(kBranchFactor, level_indices.size() - begin);
      next_level.push_back(tree.nodes.size());
      tree.nodes.push_back(
          make_parent_node(tree.nodes, first_child, child_count));
    }
    level_indices.swap(next_level);
  }

  tree.root_index = level_indices.empty() ? 0 : level_indices.front();
  tree.build_ms = elapsed_ms(build_begin, std::chrono::steady_clock::now());
  return tree;
}

void apply_direct_operator(const std::vector<Panel> &panels,
                           const std::vector<double> &density,
                           double contrast, double kappa,
                           TitanExecutor &executor,
                           std::vector<double> &out,
                           TreecodeStats *stats) {
  const std::size_t n = panels.size();
  out.assign(n, 0.0);
  if (n == 0)
    return;
  const std::size_t workers = effective_workers(executor, n);
  std::vector<std::size_t> local_near(workers, 0);
  executor.parallel_for(
      n, workers, 16,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        std::size_t near_count = 0;
        for (std::size_t i = begin; i < end; ++i) {
          const double coupled =
              direct_sum(panels[i], i, panels, density, kappa, 0, n);
          near_count += n > 0 ? n - 1 : 0;
          out[i] = panels[i].diagonal * density[i] + contrast * coupled;
        }
        local_near[worker] += near_count;
      });
  if (stats != nullptr) {
    stats->matvec_mode = "direct";
    stats->taylor_order = 0;
    stats->near_interactions +=
        std::accumulate(local_near.begin(), local_near.end(), std::size_t{0});
  }
}

void apply_treecode_operator(const std::vector<Panel> &panels,
                             const std::vector<double> &density,
                             double contrast, double kappa,
                             const TreeWorkspace &tree,
                             TitanExecutor &executor,
                             std::vector<double> &out,
                             TreecodeStats *stats) {
  const std::size_t n = panels.size();
  if (n == 0) {
    out.clear();
    return;
  }
  if (tree.nodes.empty() || tree.order.size() != n) {
    apply_direct_operator(panels, density, contrast, kappa, executor, out,
                          stats);
    return;
  }

  const auto eval_begin = std::chrono::steady_clock::now();
  const int taylor_order = std::clamp(tree.taylor_order, 0, kTaylorOrderMax);
  const auto moment_begin = std::chrono::steady_clock::now();
  const std::vector<NodeMoments> moments =
      compute_node_moments(panels, density, tree, taylor_order);
  const double moment_ms =
      elapsed_ms(moment_begin, std::chrono::steady_clock::now());
  const auto interaction_begin = std::chrono::steady_clock::now();
  std::vector<std::size_t> rank(n, 0);
  for (std::size_t pos = 0; pos < tree.order.size(); ++pos)
    rank[tree.order[pos]] = pos;
  const double interaction_ms =
      elapsed_ms(interaction_begin, std::chrono::steady_clock::now());

  out.assign(n, 0.0);
  const auto traversal_begin = std::chrono::steady_clock::now();
  const std::size_t workers = effective_workers(executor, n);
  std::vector<std::size_t> local_near(workers, 0);
  std::vector<std::size_t> local_far(workers, 0);
  executor.parallel_for(
      n, workers, 16,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        std::vector<std::size_t> stack;
        stack.reserve(64);
        std::size_t near_count = 0;
        std::size_t far_count = 0;
        for (std::size_t i = begin; i < end; ++i) {
          const Panel &target = panels[i];
          const std::size_t target_rank = rank[i];
          double coupled = 0.0;
          stack.clear();
          stack.push_back(tree.root_index);
          while (!stack.empty()) {
            const std::size_t node_index = stack.back();
            stack.pop_back();
            const TreeNode &node = tree.nodes[node_index];
            const bool leaf = node.child_count == 0;
            const bool contains_target =
                target_rank >= node.first_panel &&
                target_rank < node.first_panel + node.panel_count;
            if (!leaf && !contains_target) {
              const double dx = target.x - node.cx;
              const double dy = target.y - node.cy;
              const double dz = target.z - node.cz;
              const double distance =
                  std::sqrt(std::max(1.0e-16, dx * dx + dy * dy + dz * dz));
              if (node.radius / distance < kOpeningTheta) {
                coupled += node_kernel_taylor(target, node,
                                              moments[node_index], taylor_order,
                                              kappa);
                ++far_count;
                continue;
              }
            }
            if (leaf) {
              coupled += direct_ordered_sum(
                  target, i, panels, density, kappa, tree.order,
                  node.first_panel, node.panel_count, near_count);
            } else {
              for (std::size_t child = 0; child < node.child_count; ++child)
                stack.push_back(node.first_child + child);
            }
          }
          out[i] = panels[i].diagonal * density[i] + contrast * coupled;
        }
        local_near[worker] += near_count;
        local_far[worker] += far_count;
      });
  const double traversal_ms =
      elapsed_ms(traversal_begin, std::chrono::steady_clock::now());

  if (stats != nullptr) {
    const std::size_t near_total =
        std::accumulate(local_near.begin(), local_near.end(), std::size_t{0});
    const std::size_t far_total =
        std::accumulate(local_far.begin(), local_far.end(), std::size_t{0});
    const double interaction_total =
        static_cast<double>(near_total + far_total);
    stats->matvec_mode = "treecode";
    stats->taylor_order = taylor_order;
    stats->tree_moment_ms += moment_ms;
    stats->interaction_build_ms += interaction_ms;
    stats->tree_traversal_ms += traversal_ms;
    stats->tree_eval_ms +=
        elapsed_ms(eval_begin, std::chrono::steady_clock::now());
    if (interaction_total > 0.0) {
      stats->near_eval_ms += traversal_ms *
                             static_cast<double>(near_total) /
                             interaction_total;
      stats->far_eval_ms += traversal_ms *
                            static_cast<double>(far_total) /
                            interaction_total;
    }
    stats->near_interactions += near_total;
    stats->far_interactions += far_total;
  }
}

} // namespace TitanPBE::Tabi
