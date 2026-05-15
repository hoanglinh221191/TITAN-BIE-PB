#include "tabi_types.h"
#include "tabi_near.h"

#include "../titan_executor.h"
#include "../titan_simd_dispatch.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace TitanPBE::Tabi {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr std::size_t kBranchFactor = 8;
constexpr std::size_t kDefaultLeafSize = 32;
constexpr double kDefaultOpeningTheta = 0.45;
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

constexpr int kSeriesSide = kTaylorOrderMax + 1;
constexpr std::size_t kSeriesSize =
    kSeriesSide * kSeriesSide * kSeriesSide;

struct TaylorSeries {
  std::array<double, kSeriesSize> coeff{};
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

double point_distance(double ax, double ay, double az, double bx, double by,
                      double bz) noexcept {
  const double dx = ax - bx;
  const double dy = ay - by;
  const double dz = az - bz;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double tree_panel_max_edge_length(const Panel &panel) noexcept {
  if (!panel.has_triangle_geometry)
    return std::sqrt(std::max(0.0, panel.area));
  const double e01 = point_distance(panel.v0x, panel.v0y, panel.v0z, panel.v1x,
                                    panel.v1y, panel.v1z);
  const double e12 = point_distance(panel.v1x, panel.v1y, panel.v1z, panel.v2x,
                                    panel.v2y, panel.v2z);
  const double e20 = point_distance(panel.v2x, panel.v2y, panel.v2z, panel.v0x,
                                    panel.v0y, panel.v0z);
  return std::max({e01, e12, e20, std::sqrt(std::max(0.0, panel.area))});
}

void panel_aabb(const Panel &panel, double &min_x, double &min_y, double &min_z,
                double &max_x, double &max_y, double &max_z) noexcept {
  if (!panel.has_triangle_geometry) {
    min_x = max_x = panel.x;
    min_y = max_y = panel.y;
    min_z = max_z = panel.z;
    return;
  }
  min_x = std::min({panel.v0x, panel.v1x, panel.v2x});
  min_y = std::min({panel.v0y, panel.v1y, panel.v2y});
  min_z = std::min({panel.v0z, panel.v1z, panel.v2z});
  max_x = std::max({panel.v0x, panel.v1x, panel.v2x});
  max_y = std::max({panel.v0y, panel.v1y, panel.v2y});
  max_z = std::max({panel.v0z, panel.v1z, panel.v2z});
}

double aabb_gap(double amin_x, double amin_y, double amin_z, double amax_x,
                double amax_y, double amax_z, double bmin_x, double bmin_y,
                double bmin_z, double bmax_x, double bmax_y,
                double bmax_z) noexcept {
  const double dx = std::max({bmin_x - amax_x, amin_x - bmax_x, 0.0});
  const double dy = std::max({bmin_y - amax_y, amin_y - bmax_y, 0.0});
  const double dz = std::max({bmin_z - amax_z, amin_z - bmax_z, 0.0});
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double panel_node_aabb_gap(const Panel &target,
                           const TreeNode &node) noexcept {
  double min_x = 0.0;
  double min_y = 0.0;
  double min_z = 0.0;
  double max_x = 0.0;
  double max_y = 0.0;
  double max_z = 0.0;
  panel_aabb(target, min_x, min_y, min_z, max_x, max_y, max_z);
  return aabb_gap(min_x, min_y, min_z, max_x, max_y, max_z, node.min_x,
                  node.min_y, node.min_z, node.max_x, node.max_y, node.max_z);
}

struct MacAuditTally {
  std::size_t unsafe_count{0};
  std::size_t unsafe_panels{0};
  double min_gap{0.0};
  double min_gap_ratio{0.0};
  double sum_gap_ratio{0.0};
  std::size_t samples{0};
  std::size_t physical_gap_rejects{0};
};

void record_far_geometric_audit(MacAuditTally &tally, const TreeNode &node,
                                double gap, double gap_ratio,
                                bool unsafe) noexcept {
  ++tally.samples;
  if (!unsafe)
    return;
  ++tally.unsafe_count;
  tally.unsafe_panels += node.panel_count;
  if (tally.unsafe_count == 1 || gap < tally.min_gap)
    tally.min_gap = gap;
  if (tally.unsafe_count == 1 || gap_ratio < tally.min_gap_ratio)
    tally.min_gap_ratio = gap_ratio;
  tally.sum_gap_ratio += gap_ratio;
}

void merge_far_geometric_audit(MacAuditTally &dst,
                               const MacAuditTally &src) noexcept {
  if (src.samples == 0)
    return;
  dst.samples += src.samples;
  dst.physical_gap_rejects += src.physical_gap_rejects;
  if (src.unsafe_count == 0)
    return;
  if (dst.unsafe_count == 0 || src.min_gap < dst.min_gap)
    dst.min_gap = src.min_gap;
  if (dst.unsafe_count == 0 || src.min_gap_ratio < dst.min_gap_ratio)
    dst.min_gap_ratio = src.min_gap_ratio;
  dst.unsafe_count += src.unsafe_count;
  dst.unsafe_panels += src.unsafe_panels;
  dst.sum_gap_ratio += src.sum_gap_ratio;
}

void copy_far_geometric_audit_to_stats(TreecodeStats &stats,
                                       const MacAuditTally &tally) noexcept {
  stats.far_geometric_audit_samples += tally.samples;
  stats.mac_physical_gap_rejects += tally.physical_gap_rejects;
  if (tally.unsafe_count == 0)
    return;
  if (stats.far_geometric_unsafe_count == 0 ||
      tally.min_gap < stats.far_geometric_min_gap)
    stats.far_geometric_min_gap = tally.min_gap;
  if (stats.far_geometric_unsafe_count == 0 ||
      tally.min_gap_ratio < stats.far_geometric_gap_ratio_min)
    stats.far_geometric_gap_ratio_min = tally.min_gap_ratio;
  const double old_count =
      static_cast<double>(stats.far_geometric_unsafe_count);
  const double new_count = old_count + static_cast<double>(tally.unsafe_count);
  const double old_sum = stats.far_geometric_gap_ratio_mean * old_count;
  stats.far_geometric_unsafe_count += tally.unsafe_count;
  stats.far_geometric_unsafe_panels += tally.unsafe_panels;
  stats.far_geometric_gap_ratio_mean =
      new_count > 0.0 ? (old_sum + tally.sum_gap_ratio) / new_count : 0.0;
}

struct MacDecision {
  bool ratio_accept{false};
  bool accepted{false};
  bool unsafe{false};
  double gap{0.0};
  double gap_ratio{0.0};
};

MacDecision evaluate_mac(const Panel &target, const TreeNode &node,
                         const TreeWorkspace &tree,
                         double distance) noexcept {
  MacDecision decision;
  decision.ratio_accept = node.radius / distance < tree.opening_theta;
  decision.accepted = decision.ratio_accept;
  if (!decision.ratio_accept || tree.mac_safety == "off")
    return decision;
  const double target_edge = tree_panel_max_edge_length(target);
  const double guard_edge =
      std::max({target_edge, node.max_edge_length, 1.0e-12});
  decision.gap = panel_node_aabb_gap(target, node);
  decision.gap_ratio = decision.gap / guard_edge;
  decision.unsafe = decision.gap <= tree.mac_gap_factor * guard_edge;
  if (tree.mac_safety == "physical-gap" && decision.unsafe)
    decision.accepted = false;
  return decision;
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

template <typename Fn> void for_each_multi_index(int order, Fn &&fn) {
  const int clamped = std::clamp(order, 0, kTaylorOrderMax);
  for (int degree = 0; degree <= clamped; ++degree) {
    for (int ax = degree; ax >= 0; --ax) {
      for (int ay = degree - ax; ay >= 0; --ay) {
        const int az = degree - ax - ay;
        fn(ax, ay, az);
      }
    }
  }
}

TaylorSeries series_constant(double value) noexcept {
  TaylorSeries out;
  out.coeff[static_cast<std::size_t>(moment_lookup_index(0, 0, 0))] = value;
  return out;
}

TaylorSeries series_variable(double value, int axis) noexcept {
  TaylorSeries out = series_constant(value);
  const int ax = axis == 0 ? 1 : 0;
  const int ay = axis == 1 ? 1 : 0;
  const int az = axis == 2 ? 1 : 0;
  out.coeff[static_cast<std::size_t>(moment_lookup_index(ax, ay, az))] = 1.0;
  return out;
}

TaylorSeries series_add(const TaylorSeries &a, const TaylorSeries &b,
                        int order) noexcept {
  TaylorSeries out;
  for_each_multi_index(order, [&](int ax, int ay, int az) {
    const std::size_t idx =
        static_cast<std::size_t>(moment_lookup_index(ax, ay, az));
    out.coeff[idx] = a.coeff[idx] + b.coeff[idx];
  });
  return out;
}

TaylorSeries series_sub(const TaylorSeries &a, const TaylorSeries &b,
                        int order) noexcept {
  TaylorSeries out;
  for_each_multi_index(order, [&](int ax, int ay, int az) {
    const std::size_t idx =
        static_cast<std::size_t>(moment_lookup_index(ax, ay, az));
    out.coeff[idx] = a.coeff[idx] - b.coeff[idx];
  });
  return out;
}

TaylorSeries series_scale(const TaylorSeries &a, double scale,
                          int order) noexcept {
  TaylorSeries out;
  for_each_multi_index(order, [&](int ax, int ay, int az) {
    const std::size_t idx =
        static_cast<std::size_t>(moment_lookup_index(ax, ay, az));
    out.coeff[idx] = scale * a.coeff[idx];
  });
  return out;
}

TaylorSeries series_mul(const TaylorSeries &a, const TaylorSeries &b,
                        int order) noexcept {
  TaylorSeries out;
  for_each_multi_index(order, [&](int ax, int ay, int az) {
    double value = 0.0;
    for (int bx = 0; bx <= ax; ++bx) {
      for (int by = 0; by <= ay; ++by) {
        for (int bz = 0; bz <= az; ++bz) {
          const std::size_t lhs =
              static_cast<std::size_t>(moment_lookup_index(bx, by, bz));
          const std::size_t rhs = static_cast<std::size_t>(
              moment_lookup_index(ax - bx, ay - by, az - bz));
          value += a.coeff[lhs] * b.coeff[rhs];
        }
      }
    }
    out.coeff[static_cast<std::size_t>(moment_lookup_index(ax, ay, az))] =
        value;
  });
  return out;
}

TaylorSeries series_reciprocal(const TaylorSeries &a, int order) noexcept {
  TaylorSeries out;
  const std::size_t zero =
      static_cast<std::size_t>(moment_lookup_index(0, 0, 0));
  const double a0 = a.coeff[zero];
  if (!std::isfinite(a0) || std::abs(a0) <= 1.0e-30)
    return out;
  out.coeff[zero] = 1.0 / a0;
  for_each_multi_index(order, [&](int ax, int ay, int az) {
    if (ax == 0 && ay == 0 && az == 0)
      return;
    double sum = 0.0;
    for (int bx = 0; bx <= ax; ++bx) {
      for (int by = 0; by <= ay; ++by) {
        for (int bz = 0; bz <= az; ++bz) {
          if (bx == 0 && by == 0 && bz == 0)
            continue;
          const std::size_t lhs =
              static_cast<std::size_t>(moment_lookup_index(bx, by, bz));
          const std::size_t rhs = static_cast<std::size_t>(
              moment_lookup_index(ax - bx, ay - by, az - bz));
          sum += a.coeff[lhs] * out.coeff[rhs];
        }
      }
    }
    out.coeff[static_cast<std::size_t>(moment_lookup_index(ax, ay, az))] =
        -sum / a0;
  });
  return out;
}

TaylorSeries series_sqrt(const TaylorSeries &a, int order) noexcept {
  TaylorSeries out;
  const std::size_t zero =
      static_cast<std::size_t>(moment_lookup_index(0, 0, 0));
  const double a0 = a.coeff[zero];
  if (!std::isfinite(a0) || a0 <= 0.0)
    return out;
  const double root = std::sqrt(a0);
  out.coeff[zero] = root;
  for_each_multi_index(order, [&](int ax, int ay, int az) {
    if (ax == 0 && ay == 0 && az == 0)
      return;
    double sum = 0.0;
    for (int bx = 0; bx <= ax; ++bx) {
      for (int by = 0; by <= ay; ++by) {
        for (int bz = 0; bz <= az; ++bz) {
          if ((bx == 0 && by == 0 && bz == 0) ||
              (bx == ax && by == ay && bz == az))
            continue;
          const std::size_t lhs =
              static_cast<std::size_t>(moment_lookup_index(bx, by, bz));
          const std::size_t rhs = static_cast<std::size_t>(
              moment_lookup_index(ax - bx, ay - by, az - bz));
          sum += out.coeff[lhs] * out.coeff[rhs];
        }
      }
    }
    const std::size_t idx =
        static_cast<std::size_t>(moment_lookup_index(ax, ay, az));
    out.coeff[idx] = (a.coeff[idx] - sum) / (2.0 * root);
  });
  return out;
}

TaylorSeries series_exp(const TaylorSeries &a, int order) noexcept {
  TaylorSeries out;
  const std::size_t zero =
      static_cast<std::size_t>(moment_lookup_index(0, 0, 0));
  out.coeff[zero] = std::exp(a.coeff[zero]);
  for_each_multi_index(order, [&](int ax, int ay, int az) {
    if (ax == 0 && ay == 0 && az == 0)
      return;
    int axis = 0;
    if (ax == 0)
      axis = ay > 0 ? 1 : 2;
    const int gx = ax - (axis == 0 ? 1 : 0);
    const int gy = ay - (axis == 1 ? 1 : 0);
    const int gz = az - (axis == 2 ? 1 : 0);
    double rhs = 0.0;
    for (int bx = 0; bx <= gx; ++bx) {
      for (int by = 0; by <= gy; ++by) {
        for (int bz = 0; bz <= gz; ++bz) {
          const int dx = gx - bx + (axis == 0 ? 1 : 0);
          const int dy = gy - by + (axis == 1 ? 1 : 0);
          const int dz = gz - bz + (axis == 2 ? 1 : 0);
          const int axis_degree =
              axis == 0 ? dx : (axis == 1 ? dy : dz);
          const std::size_t lhs =
              static_cast<std::size_t>(moment_lookup_index(bx, by, bz));
          const std::size_t deriv =
              static_cast<std::size_t>(moment_lookup_index(dx, dy, dz));
          rhs += out.coeff[lhs] * static_cast<double>(axis_degree) *
                 a.coeff[deriv];
        }
      }
    }
    const int alpha_axis = axis == 0 ? ax : (axis == 1 ? ay : az);
    out.coeff[static_cast<std::size_t>(moment_lookup_index(ax, ay, az))] =
        rhs / static_cast<double>(alpha_axis);
  });
  return out;
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

struct KernelSeriesBase {
  TaylorSeries ux;
  TaylorSeries uy;
  TaylorSeries uz;
  TaylorSeries r;
  TaylorSeries inv_r;
  TaylorSeries inv_r2;
  TaylorSeries inv_r3;
  TaylorSeries screened;
};

KernelSeriesBase make_kernel_series_base(double tx, double ty, double tz,
                                         double sx, double sy, double sz,
                                         double kappa, int order) noexcept {
  KernelSeriesBase base;
  const TaylorSeries source_x = series_variable(sx, 0);
  const TaylorSeries source_y = series_variable(sy, 1);
  const TaylorSeries source_z = series_variable(sz, 2);
  base.ux = series_sub(series_constant(tx), source_x, order);
  base.uy = series_sub(series_constant(ty), source_y, order);
  base.uz = series_sub(series_constant(tz), source_z, order);
  const TaylorSeries r2 = series_add(
      series_add(series_mul(base.ux, base.ux, order),
                 series_mul(base.uy, base.uy, order), order),
      series_mul(base.uz, base.uz, order), order);
  const double r2_zero =
      r2.coeff[static_cast<std::size_t>(moment_lookup_index(0, 0, 0))];
  if (!std::isfinite(r2_zero) || r2_zero <= 1.0e-28)
    return base;
  base.r = series_sqrt(r2, order);
  base.inv_r = series_reciprocal(base.r, order);
  base.inv_r2 = series_mul(base.inv_r, base.inv_r, order);
  base.inv_r3 = series_mul(base.inv_r2, base.inv_r, order);
  base.screened = series_exp(series_scale(base.r, -std::max(0.0, kappa),
                                          order),
                             order);
  return base;
}

TaylorSeries green_series(double tx, double ty, double tz, double sx,
                          double sy, double sz, double kappa,
                          int order) noexcept {
  const KernelSeriesBase base =
      make_kernel_series_base(tx, ty, tz, sx, sy, sz, kappa, order);
  return series_scale(series_mul(base.screened, base.inv_r, order),
                      1.0 / (4.0 * kPi), order);
}

TaylorSeries source_normal_component_green_series(double tx, double ty,
                                                  double tz, double sx,
                                                  double sy, double sz,
                                                  int component, double kappa,
                                                  int order) noexcept {
  const double kk = std::max(0.0, kappa);
  const KernelSeriesBase base =
      make_kernel_series_base(tx, ty, tz, sx, sy, sz, kk, order);
  const TaylorSeries one_plus_kr =
      series_add(series_constant(1.0), series_scale(base.r, kk, order), order);
  const TaylorSeries &uc =
      component == 0 ? base.ux : (component == 1 ? base.uy : base.uz);
  TaylorSeries value =
      series_mul(series_mul(base.screened, one_plus_kr, order),
                 series_mul(uc, base.inv_r3, order), order);
  return series_scale(value, 1.0 / (4.0 * kPi), order);
}

TaylorSeries target_normal_green_series(const Panel &target, double sx,
                                        double sy, double sz, double kappa,
                                        int order) noexcept {
  const double kk = std::max(0.0, kappa);
  const KernelSeriesBase base = make_kernel_series_base(
      target.x, target.y, target.z, sx, sy, sz, kk, order);
  const TaylorSeries normal_dot_u = series_add(
      series_add(series_scale(base.ux, target.nx, order),
                 series_scale(base.uy, target.ny, order), order),
      series_scale(base.uz, target.nz, order), order);
  const TaylorSeries one_plus_kr =
      series_add(series_constant(1.0), series_scale(base.r, kk, order), order);
  TaylorSeries value =
      series_mul(series_mul(base.screened, one_plus_kr, order),
                 series_mul(normal_dot_u, base.inv_r3, order), order);
  return series_scale(value, -1.0 / (4.0 * kPi), order);
}

TaylorSeries mixed_normal_component_green_series(const Panel &target,
                                                 double sx, double sy,
                                                 double sz, int component,
                                                 double kappa,
                                                 int order) noexcept {
  const double kk = std::max(0.0, kappa);
  const KernelSeriesBase base = make_kernel_series_base(
      target.x, target.y, target.z, sx, sy, sz, kk, order);
  const TaylorSeries normal_dot_u = series_add(
      series_add(series_scale(base.ux, target.nx, order),
                 series_scale(base.uy, target.ny, order), order),
      series_scale(base.uz, target.nz, order), order);
  const TaylorSeries nxu = series_mul(normal_dot_u, base.inv_r, order);
  const TaylorSeries &uc =
      component == 0 ? base.ux : (component == 1 ? base.uy : base.uz);
  const TaylorSeries nyu = series_mul(uc, base.inv_r, order);
  const double nxny =
      component == 0 ? target.nx : (component == 1 ? target.ny : target.nz);
  const TaylorSeries one_plus_kr =
      series_add(series_constant(1.0), series_scale(base.r, kk, order), order);
  const TaylorSeries fprime = series_scale(
      series_mul(series_mul(base.screened, one_plus_kr, order), base.inv_r2,
                 order),
      -1.0 / (4.0 * kPi), order);
  TaylorSeries fsecond_radial = series_add(
      series_add(series_scale(base.inv_r, kk * kk, order),
                 series_scale(base.inv_r2, 2.0 * kk, order), order),
      series_scale(base.inv_r3, 2.0, order), order);
  const TaylorSeries fsecond = series_scale(
      series_mul(base.screened, fsecond_radial, order), 1.0 / (4.0 * kPi),
      order);
  const TaylorSeries nxu_nyu = series_mul(nxu, nyu, order);
  const TaylorSeries metric =
      series_sub(series_constant(nxny), nxu_nyu, order);
  return series_add(series_mul(fsecond, nxu_nyu, order),
                    series_mul(series_mul(fprime, base.inv_r, order), metric,
                               order),
                    order);
}

TaylorSeries point_kernel_series(const Panel &target, double sx, double sy,
                                 double sz, double kappa,
                                 int order) noexcept {
  const KernelSeriesBase base = make_kernel_series_base(
      target.x, target.y, target.z, sx, sy, sz, kappa, order);
  const TaylorSeries normal_dot_u = series_add(
      series_add(series_scale(base.ux, target.nx, order),
                 series_scale(base.uy, target.ny, order), order),
      series_scale(base.uz, target.nz, order), order);
  return series_scale(series_mul(base.screened,
                                 series_mul(normal_dot_u, base.inv_r3, order),
                                 order),
                      1.0 / (4.0 * kPi), order);
}

double taylor_sum_series(const NodeMoments &moments, int order,
                         const TaylorSeries &series) noexcept {
  const auto &indices = moment_indices(order);
  double value = 0.0;
  for (std::size_t m = 0; m < indices.size(); ++m) {
    const MomentIndex &index = indices[m];
    const std::size_t coeff =
        static_cast<std::size_t>(moment_lookup_index(index.ax, index.ay,
                                                     index.az));
    value += moments.values[m] * series.coeff[coeff];
  }
  return std::isfinite(value) ? value : 0.0;
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
                          double kappa,
                          const std::string &derivative_mode) noexcept {
  if (derivative_mode != "finite-diff-debug") {
    return taylor_sum_series(
        moments, order,
        point_kernel_series(target, node.cx, node.cy, node.cz, kappa, order));
  }
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

double green_value(double tx, double ty, double tz, double sx, double sy,
                   double sz, double kappa) noexcept {
  const double dx = tx - sx;
  const double dy = ty - sy;
  const double dz = tz - sz;
  const double r2 = dx * dx + dy * dy + dz * dz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  return std::exp(-std::max(0.0, kappa) * r) / (4.0 * kPi * r);
}

double source_normal_component_green(double tx, double ty, double tz,
                                     double sx, double sy, double sz,
                                     int component,
                                     double kappa) noexcept {
  const double ux = tx - sx;
  const double uy = ty - sy;
  const double uz = tz - sz;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  const double n_dot_u = component == 0 ? ux : (component == 1 ? uy : uz);
  return std::exp(-std::max(0.0, kappa) * r) *
         (1.0 + std::max(0.0, kappa) * r) * n_dot_u /
         (4.0 * kPi * r2 * r);
}

double target_normal_green(const Panel &target, double sx, double sy,
                           double sz, double kappa) noexcept {
  const double ux = target.x - sx;
  const double uy = target.y - sy;
  const double uz = target.z - sz;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  const double n_dot_u = target.nx * ux + target.ny * uy + target.nz * uz;
  const double kk = std::max(0.0, kappa);
  return -std::exp(-kk * r) * (1.0 + kk * r) * n_dot_u /
         (4.0 * kPi * r2 * r);
}

double mixed_normal_component_green(const Panel &target, double sx, double sy,
                                    double sz, int component,
                                    double kappa) noexcept {
  const double ux = target.x - sx;
  const double uy = target.y - sy;
  const double uz = target.z - sz;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  const double inv_r = 1.0 / r;
  const double nxu =
      (target.nx * ux + target.ny * uy + target.nz * uz) * inv_r;
  const double nyu =
      (component == 0 ? ux : (component == 1 ? uy : uz)) * inv_r;
  const double nxny =
      component == 0 ? target.nx : (component == 1 ? target.ny : target.nz);
  const double kk = std::max(0.0, kappa);
  const double expkr = std::exp(-kk * r);
  const double inv_r2 = inv_r * inv_r;
  const double fprime =
      -expkr * (1.0 + kk * r) * inv_r2 / (4.0 * kPi);
  const double fsecond =
      expkr *
      (kk * kk * inv_r + 2.0 * kk * inv_r2 + 2.0 * inv_r2 * inv_r) /
      (4.0 * kPi);
  return fsecond * nxu * nyu +
         (fprime * inv_r) * (nxny - nxu * nyu);
}

template <typename Kernel>
double finite_difference_source_derivative(const Kernel &kernel, double sx,
                                           double sy, double sz, int ax,
                                           int ay, int az, double h) noexcept {
  if (ax == 0 && ay == 0 && az == 0)
    return kernel(sx, sy, sz);
  if (ax > 0) {
    return (finite_difference_source_derivative(kernel, sx + h, sy, sz,
                                                ax - 1, ay, az, h) -
            finite_difference_source_derivative(kernel, sx - h, sy, sz,
                                                ax - 1, ay, az, h)) /
           (2.0 * h);
  }
  if (ay > 0) {
    return (finite_difference_source_derivative(kernel, sx, sy + h, sz, ax,
                                                ay - 1, az, h) -
            finite_difference_source_derivative(kernel, sx, sy - h, sz, ax,
                                                ay - 1, az, h)) /
           (2.0 * h);
  }
  return (finite_difference_source_derivative(kernel, sx, sy, sz + h, ax, ay,
                                              az - 1, h) -
          finite_difference_source_derivative(kernel, sx, sy, sz - h, ax, ay,
                                              az - 1, h)) /
         (2.0 * h);
}

template <typename Kernel>
double taylor_sum_kernel(const TreeNode &node, const NodeMoments &moments,
                         int order, double h, const Kernel &kernel) noexcept {
  const auto &indices = moment_indices(order);
  double value = 0.0;
  for (std::size_t m = 0; m < indices.size(); ++m) {
    const MomentIndex &index = indices[m];
    const double derivative = finite_difference_source_derivative(
        kernel, node.cx, node.cy, node.cz, index.ax, index.ay, index.az, h);
    value += moments.values[m] * derivative * index.inv_factorial;
  }
  return std::isfinite(value) ? value : 0.0;
}

struct JcpAccum {
  double phi{0.0};
  double phin{0.0};
};

struct TreeUnknownCache {
  std::vector<double> phi;
  std::vector<double> phin;
  AlignedDoubleVector ordered_phi;
  AlignedDoubleVector ordered_phin;
};

TreePanelCache build_tree_panel_cache(const std::vector<Panel> &panels,
                                      const std::vector<std::size_t> &order) {
  TreePanelCache cache;
  const std::size_t n = panels.size();
  cache.x.resize(n);
  cache.y.resize(n);
  cache.z.resize(n);
  cache.nx.resize(n);
  cache.ny.resize(n);
  cache.nz.resize(n);
  cache.area.resize(n);
  cache.ordered_x.resize(n);
  cache.ordered_y.resize(n);
  cache.ordered_z.resize(n);
  cache.ordered_nx.resize(n);
  cache.ordered_ny.resize(n);
  cache.ordered_nz.resize(n);
  cache.ordered_area.resize(n);
  cache.ordered_original_index.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const Panel &panel = panels[i];
    cache.x[i] = panel.x;
    cache.y[i] = panel.y;
    cache.z[i] = panel.z;
    cache.nx[i] = panel.nx;
    cache.ny[i] = panel.ny;
    cache.nz[i] = panel.nz;
    cache.area[i] = panel.area;
  }
  for (std::size_t pos = 0; pos < n; ++pos) {
    const std::size_t original = pos < order.size() ? order[pos] : pos;
    const Panel &panel = panels[original];
    cache.ordered_x[pos] = panel.x;
    cache.ordered_y[pos] = panel.y;
    cache.ordered_z[pos] = panel.z;
    cache.ordered_nx[pos] = panel.nx;
    cache.ordered_ny[pos] = panel.ny;
    cache.ordered_nz[pos] = panel.nz;
    cache.ordered_area[pos] = panel.area;
    cache.ordered_original_index[pos] = original;
  }
  return cache;
}

TreeUnknownCache build_tree_unknown_cache(const std::vector<Panel> &panels,
                                          const std::vector<double> &unknowns,
                                          bool coupled,
                                          const std::vector<std::size_t> &order) {
  TreeUnknownCache cache;
  const std::size_t n = panels.size();
  cache.phi.assign(n, 0.0);
  cache.phin.assign(n, 0.0);
  cache.ordered_phi.assign(n, 0.0);
  cache.ordered_phin.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    cache.phi[i] = i < unknowns.size() ? unknowns[i] : 0.0;
    if (coupled && n + i < unknowns.size())
      cache.phin[i] = unknowns[n + i];
  }
  for (std::size_t pos = 0; pos < n; ++pos) {
    const std::size_t original = pos < order.size() ? order[pos] : pos;
    cache.ordered_phi[pos] = original < cache.phi.size() ? cache.phi[original] : 0.0;
    cache.ordered_phin[pos] =
        original < cache.phin.size() ? cache.phin[original] : 0.0;
  }
  return cache;
}

inline bool accumulate_jcp_leaf_source(
    JcpAccum &acc, double tx, double ty, double tz, double tnx, double tny,
    double tnz, double sx, double sy, double sz, double snx, double sny,
    double snz, double area, double phi_j, double phin_j, double eps_ratio,
    double inv_eps_ratio, double poisson_scale, double inv_four_pi,
    double kk, bool coupled) noexcept {
  const double ux = tx - sx;
  const double uy = ty - sy;
  const double uz = tz - sz;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return false;
  const double r = std::sqrt(r2);
  const double inv_r = 1.0 / r;
  const double inv_r2 = inv_r * inv_r;
  const double inv_r3 = inv_r2 * inv_r;
  const double source_projection = snx * ux + sny * uy + snz * uz;
  const double target_projection = tnx * ux + tny * uy + tnz * uz;
  const double source_d0 = inv_four_pi * source_projection * inv_r3;
  if (coupled) {
    const double expkr = std::exp(-kk * r);
    const double one_plus_kr = 1.0 + kk * r;
    const double source_dk = expkr * one_plus_kr * source_d0;
    const double g0 = inv_four_pi * inv_r;
    const double gk = expkr * g0;
    const double target_d0 = -inv_four_pi * target_projection * inv_r3;
    const double target_dk = expkr * one_plus_kr * target_d0;
    const double nxu = target_projection * inv_r;
    const double nyu = source_projection * inv_r;
    const double nxny = tnx * snx + tny * sny + tnz * snz;
    const double nxu_nyu = nxu * nyu;
    const double tangent_metric = nxny - nxu_nyu;
    const double mixed0 = inv_four_pi * inv_r3 * (3.0 * nxu_nyu - nxny);
    const double fprime_k_over_r =
        -expkr * one_plus_kr * inv_four_pi * inv_r3;
    const double fsecond_k =
        expkr * inv_four_pi *
        (kk * kk * inv_r + 2.0 * kk * inv_r2 + 2.0 * inv_r3);
    const double mixedk =
        fsecond_k * nxu_nyu + fprime_k_over_r * tangent_metric;
    acc.phi += area * ((g0 - gk) * phin_j +
                       (eps_ratio * source_dk - source_d0) * phi_j);
    acc.phin += area * ((target_d0 - inv_eps_ratio * target_dk) * phin_j +
                        (mixedk - mixed0) * phi_j);
  } else {
    acc.phi += area * poisson_scale * source_d0 * phi_j;
  }
  return true;
}

JcpAccum direct_jcp_ordered_sum_scalar(
    const Panel &target, std::size_t target_index,
    const TreePanelCache &panel_cache, const TreeUnknownCache &unknown_cache,
    double eps_ratio, double kappa, bool coupled,
    const std::vector<std::size_t> &order, std::size_t begin,
    std::size_t count, std::size_t &near_count) noexcept {
  JcpAccum acc;
  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  const double poisson_scale = eps_ratio - 1.0;
  const double inv_four_pi = 1.0 / (4.0 * kPi);
  const double kk = std::max(0.0, kappa);
  const std::size_t end = begin + count;
  for (std::size_t pos = begin; pos < end; ++pos) {
    const std::size_t j = order[pos];
    if (j == target_index)
      continue;
    if (accumulate_jcp_leaf_source(
        acc, target.x, target.y, target.z, target.nx, target.ny, target.nz,
        panel_cache.x[j], panel_cache.y[j], panel_cache.z[j],
        panel_cache.nx[j], panel_cache.ny[j], panel_cache.nz[j],
        panel_cache.area[j], unknown_cache.phi[j], unknown_cache.phin[j],
        eps_ratio, inv_eps_ratio, poisson_scale, inv_four_pi, kk, coupled))
      ++near_count;
  }
  return acc;
}

inline bool accumulate_jcp_ordered_tile_position(
    JcpAccum &acc, const Panel &target, std::size_t target_index,
    const TreePanelCache &panel_cache, const TreeUnknownCache &unknown_cache,
    double eps_ratio, double inv_eps_ratio, double poisson_scale,
    double inv_four_pi, double kk, bool coupled, std::size_t pos) noexcept {
  const std::size_t original =
      pos < panel_cache.ordered_original_index.size()
          ? panel_cache.ordered_original_index[pos]
          : pos;
  if (original == target_index)
    return false;
  return accumulate_jcp_leaf_source(
      acc, target.x, target.y, target.z, target.nx, target.ny, target.nz,
      panel_cache.ordered_x[pos], panel_cache.ordered_y[pos],
      panel_cache.ordered_z[pos], panel_cache.ordered_nx[pos],
      panel_cache.ordered_ny[pos], panel_cache.ordered_nz[pos],
      panel_cache.ordered_area[pos], unknown_cache.ordered_phi[pos],
      unknown_cache.ordered_phin[pos], eps_ratio, inv_eps_ratio,
      poisson_scale, inv_four_pi, kk, coupled);
}

bool tile_position_is_excluded(const TreePanelCache &panel_cache,
                               std::size_t target_index, std::size_t pos,
                               const std::vector<std::size_t> *exclusions,
                               std::size_t exclusion_begin,
                               std::size_t exclusion_count) noexcept {
  const std::size_t original =
      pos < panel_cache.ordered_original_index.size()
          ? panel_cache.ordered_original_index[pos]
          : pos;
  if (original == target_index)
    return true;
  if (exclusions == nullptr || exclusion_count == 0)
    return false;
  const std::size_t end = exclusion_begin + exclusion_count;
  for (std::size_t k = exclusion_begin; k < end; ++k) {
    if ((*exclusions)[k] == original)
      return true;
  }
  return false;
}

inline bool accumulate_jcp_ordered_tile_position_masked(
    JcpAccum &acc, const Panel &target, std::size_t target_index,
    const TreePanelCache &panel_cache, const TreeUnknownCache &unknown_cache,
    double eps_ratio, double inv_eps_ratio, double poisson_scale,
    double inv_four_pi, double kk, bool coupled, std::size_t pos,
    const std::vector<std::size_t> *exclusions, std::size_t exclusion_begin,
    std::size_t exclusion_count) noexcept {
  if (tile_position_is_excluded(panel_cache, target_index, pos, exclusions,
                                exclusion_begin, exclusion_count))
    return false;
  return accumulate_jcp_ordered_tile_position(
      acc, target, target_index, panel_cache, unknown_cache, eps_ratio,
      inv_eps_ratio, poisson_scale, inv_four_pi, kk, coupled, pos);
}

JcpAccum direct_jcp_ordered_tile_sum_scalar(
    const Panel &target, std::size_t target_index,
    const TreePanelCache &panel_cache, const TreeUnknownCache &unknown_cache,
    double eps_ratio, double kappa, bool coupled, const TreeNode &node,
    std::size_t &near_count,
    const std::vector<std::size_t> *exclusions = nullptr,
    std::size_t exclusion_begin = 0,
    std::size_t exclusion_count = 0) noexcept {
  JcpAccum acc;
  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  const double poisson_scale = eps_ratio - 1.0;
  const double inv_four_pi = 1.0 / (4.0 * kPi);
  const double kk = std::max(0.0, kappa);
  const std::size_t end = node.first_panel + node.panel_count;
  for (std::size_t pos = node.first_panel; pos < end; ++pos) {
    if (accumulate_jcp_ordered_tile_position_masked(
            acc, target, target_index, panel_cache, unknown_cache, eps_ratio,
            inv_eps_ratio, poisson_scale, inv_four_pi, kk, coupled, pos,
            exclusions, exclusion_begin, exclusion_count))
      ++near_count;
  }
  return acc;
}

#if TITAN_SIMD_X86
TITAN_TARGET_AVX2
JcpAccum direct_jcp_ordered_tile_sum_avx2(
    const Panel &target, std::size_t target_index,
    const TreePanelCache &panel_cache, const TreeUnknownCache &unknown_cache,
    double eps_ratio, double kappa, bool coupled, const TreeNode &node,
    std::size_t &near_count,
    const std::vector<std::size_t> *exclusions = nullptr,
    std::size_t exclusion_begin = 0,
    std::size_t exclusion_count = 0) noexcept {
  JcpAccum acc;
  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  const double poisson_scale = eps_ratio - 1.0;
  const double inv_four_pi = 1.0 / (4.0 * kPi);
  const double kk = std::max(0.0, kappa);
  const std::size_t end = node.first_panel + node.panel_count;
  std::size_t pos = node.first_panel;
  while (pos < end && (pos & std::size_t{3}) != 0) {
    if (accumulate_jcp_ordered_tile_position_masked(
            acc, target, target_index, panel_cache, unknown_cache, eps_ratio,
            inv_eps_ratio, poisson_scale, inv_four_pi, kk, coupled, pos,
            exclusions, exclusion_begin, exclusion_count))
      ++near_count;
    ++pos;
  }

  const __m256d tx = _mm256_set1_pd(target.x);
  const __m256d ty = _mm256_set1_pd(target.y);
  const __m256d tz = _mm256_set1_pd(target.z);
  const __m256d tnx = _mm256_set1_pd(target.nx);
  const __m256d tny = _mm256_set1_pd(target.ny);
  const __m256d tnz = _mm256_set1_pd(target.nz);
  const __m256d eps_v = _mm256_set1_pd(eps_ratio);
  const __m256d inv_eps_v = _mm256_set1_pd(inv_eps_ratio);
  const __m256d poisson_scale_v = _mm256_set1_pd(poisson_scale);
  const __m256d inv_four_pi_v = _mm256_set1_pd(inv_four_pi);
  const __m256d kk_v = _mm256_set1_pd(kk);
  const __m256d one_v = _mm256_set1_pd(1.0);
  const __m256d two_v = _mm256_set1_pd(2.0);
  const __m256d three_v = _mm256_set1_pd(3.0);
  const __m256d neg_one_v = _mm256_set1_pd(-1.0);
  alignas(32) double r_lane[4];
  alignas(32) double exp_lane[4];
  alignas(32) double phi_lane[4];
  alignas(32) double phin_lane[4];

  for (; pos + 4 <= end; pos += 4) {
    bool masked_lane = false;
    for (std::size_t lane = 0; lane < 4; ++lane) {
      if (tile_position_is_excluded(panel_cache, target_index, pos + lane,
                                    exclusions, exclusion_begin,
                                    exclusion_count)) {
        masked_lane = true;
        break;
      }
    }
    if (masked_lane) {
      for (std::size_t lane = 0; lane < 4; ++lane) {
        if (accumulate_jcp_ordered_tile_position_masked(
                acc, target, target_index, panel_cache, unknown_cache,
                eps_ratio, inv_eps_ratio, poisson_scale, inv_four_pi, kk,
                coupled, pos + lane, exclusions, exclusion_begin,
                exclusion_count))
          ++near_count;
      }
      continue;
    }
    const __m256d sx = _mm256_load_pd(panel_cache.ordered_x.data() + pos);
    const __m256d sy = _mm256_load_pd(panel_cache.ordered_y.data() + pos);
    const __m256d sz = _mm256_load_pd(panel_cache.ordered_z.data() + pos);
    const __m256d snx = _mm256_load_pd(panel_cache.ordered_nx.data() + pos);
    const __m256d sny = _mm256_load_pd(panel_cache.ordered_ny.data() + pos);
    const __m256d snz = _mm256_load_pd(panel_cache.ordered_nz.data() + pos);
    const __m256d area =
        _mm256_load_pd(panel_cache.ordered_area.data() + pos);
    const __m256d phi_j =
        _mm256_load_pd(unknown_cache.ordered_phi.data() + pos);
    const __m256d phin_j =
        _mm256_load_pd(unknown_cache.ordered_phin.data() + pos);
    const __m256d ux = _mm256_sub_pd(tx, sx);
    const __m256d uy = _mm256_sub_pd(ty, sy);
    const __m256d uz = _mm256_sub_pd(tz, sz);
    const __m256d r2 = _mm256_add_pd(
        _mm256_add_pd(_mm256_mul_pd(ux, ux), _mm256_mul_pd(uy, uy)),
        _mm256_mul_pd(uz, uz));
    const __m256d r = _mm256_sqrt_pd(r2);
    const __m256d inv_r = _mm256_div_pd(one_v, r);
    const __m256d inv_r2 = _mm256_mul_pd(inv_r, inv_r);
    const __m256d inv_r3 = _mm256_mul_pd(inv_r2, inv_r);
    const __m256d source_projection = _mm256_add_pd(
        _mm256_add_pd(_mm256_mul_pd(snx, ux), _mm256_mul_pd(sny, uy)),
        _mm256_mul_pd(snz, uz));
    const __m256d target_projection = _mm256_add_pd(
        _mm256_add_pd(_mm256_mul_pd(tnx, ux), _mm256_mul_pd(tny, uy)),
        _mm256_mul_pd(tnz, uz));
    const __m256d source_d0 = _mm256_mul_pd(
        inv_four_pi_v, _mm256_mul_pd(source_projection, inv_r3));
    if (coupled) {
      _mm256_store_pd(r_lane, r);
      for (int lane = 0; lane < 4; ++lane)
        exp_lane[lane] = std::exp(-kk * r_lane[lane]);
      const __m256d expkr = _mm256_load_pd(exp_lane);
      const __m256d one_plus_kr =
          _mm256_add_pd(one_v, _mm256_mul_pd(kk_v, r));
      const __m256d source_dk =
          _mm256_mul_pd(_mm256_mul_pd(expkr, one_plus_kr), source_d0);
      const __m256d g0 = _mm256_mul_pd(inv_four_pi_v, inv_r);
      const __m256d gk = _mm256_mul_pd(expkr, g0);
      const __m256d target_d0 = _mm256_mul_pd(
          neg_one_v,
          _mm256_mul_pd(inv_four_pi_v,
                        _mm256_mul_pd(target_projection, inv_r3)));
      const __m256d target_dk =
          _mm256_mul_pd(_mm256_mul_pd(expkr, one_plus_kr), target_d0);
      const __m256d nxu = _mm256_mul_pd(target_projection, inv_r);
      const __m256d nyu = _mm256_mul_pd(source_projection, inv_r);
      const __m256d nxny = _mm256_add_pd(
          _mm256_add_pd(_mm256_mul_pd(tnx, snx), _mm256_mul_pd(tny, sny)),
          _mm256_mul_pd(tnz, snz));
      const __m256d nxu_nyu = _mm256_mul_pd(nxu, nyu);
      const __m256d tangent_metric = _mm256_sub_pd(nxny, nxu_nyu);
      const __m256d mixed0 = _mm256_mul_pd(
          inv_four_pi_v,
          _mm256_mul_pd(inv_r3,
                        _mm256_sub_pd(_mm256_mul_pd(three_v, nxu_nyu),
                                      nxny)));
      const __m256d fprime_k_over_r = _mm256_mul_pd(
          neg_one_v,
          _mm256_mul_pd(
              _mm256_mul_pd(expkr, one_plus_kr),
              _mm256_mul_pd(inv_four_pi_v, inv_r3)));
      const __m256d fsecond_k = _mm256_mul_pd(
          expkr,
          _mm256_mul_pd(
              inv_four_pi_v,
              _mm256_add_pd(
                  _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(kk_v, kk_v),
                                              inv_r),
                                _mm256_mul_pd(_mm256_mul_pd(two_v, kk_v),
                                              inv_r2)),
                  _mm256_mul_pd(two_v, inv_r3))));
      const __m256d mixedk = _mm256_add_pd(
          _mm256_mul_pd(fsecond_k, nxu_nyu),
          _mm256_mul_pd(fprime_k_over_r, tangent_metric));
      const __m256d phi_term = _mm256_mul_pd(
          area,
          _mm256_add_pd(
              _mm256_mul_pd(_mm256_sub_pd(g0, gk), phin_j),
              _mm256_mul_pd(
                  _mm256_sub_pd(_mm256_mul_pd(eps_v, source_dk), source_d0),
                  phi_j)));
      const __m256d phin_term = _mm256_mul_pd(
          area,
          _mm256_add_pd(
              _mm256_mul_pd(
                  _mm256_sub_pd(target_d0,
                                _mm256_mul_pd(inv_eps_v, target_dk)),
                  phin_j),
              _mm256_mul_pd(_mm256_sub_pd(mixedk, mixed0), phi_j)));
      _mm256_store_pd(phi_lane, phi_term);
      _mm256_store_pd(phin_lane, phin_term);
      acc.phi += phi_lane[0] + phi_lane[1] + phi_lane[2] + phi_lane[3];
      acc.phin +=
          phin_lane[0] + phin_lane[1] + phin_lane[2] + phin_lane[3];
    } else {
      const __m256d phi_term = _mm256_mul_pd(
          area, _mm256_mul_pd(poisson_scale_v,
                              _mm256_mul_pd(source_d0, phi_j)));
      _mm256_store_pd(phi_lane, phi_term);
      acc.phi += phi_lane[0] + phi_lane[1] + phi_lane[2] + phi_lane[3];
    }
    near_count += 4;
  }

  for (; pos < end; ++pos) {
    if (accumulate_jcp_ordered_tile_position_masked(
            acc, target, target_index, panel_cache, unknown_cache, eps_ratio,
            inv_eps_ratio, poisson_scale, inv_four_pi, kk, coupled, pos,
            exclusions, exclusion_begin, exclusion_count))
      ++near_count;
  }
  return acc;
}
#endif

JcpAccum direct_jcp_ordered_tile_sum(
    const Panel &target, std::size_t target_index,
    const TreePanelCache &panel_cache, const TreeUnknownCache &unknown_cache,
    double eps_ratio, double kappa, bool coupled, const TreeNode &node,
    const std::string &kernel, std::size_t &near_count,
    const std::vector<std::size_t> *exclusions = nullptr,
    std::size_t exclusion_begin = 0,
    std::size_t exclusion_count = 0) noexcept {
#if TITAN_SIMD_X86 && TITAN_RUNTIME_TARGET_ATTRS
  if (kernel == "tile-avx2") {
    return direct_jcp_ordered_tile_sum_avx2(
        target, target_index, panel_cache, unknown_cache, eps_ratio, kappa,
        coupled, node, near_count, exclusions, exclusion_begin,
        exclusion_count);
  }
#endif
  return direct_jcp_ordered_tile_sum_scalar(
      target, target_index, panel_cache, unknown_cache, eps_ratio, kappa,
      coupled, node, near_count, exclusions, exclusion_begin,
      exclusion_count);
}

#if TITAN_SIMD_X86
TITAN_TARGET_AVX2
JcpAccum direct_jcp_ordered_sum_avx2(
    const Panel &target, std::size_t target_index,
    const TreePanelCache &panel_cache, const TreeUnknownCache &unknown_cache,
    double eps_ratio, double kappa, bool coupled,
    const std::vector<std::size_t> &order, std::size_t begin,
    std::size_t count, std::size_t &near_count) noexcept {
  JcpAccum acc;
  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  const double poisson_scale = eps_ratio - 1.0;
  const double inv_four_pi = 1.0 / (4.0 * kPi);
  const double kk = std::max(0.0, kappa);
  const std::size_t end = begin + count;
  const __m256d tx = _mm256_set1_pd(target.x);
  const __m256d ty = _mm256_set1_pd(target.y);
  const __m256d tz = _mm256_set1_pd(target.z);
  alignas(32) double sx_lane[4];
  alignas(32) double sy_lane[4];
  alignas(32) double sz_lane[4];
  alignas(32) double snx_lane[4];
  alignas(32) double sny_lane[4];
  alignas(32) double snz_lane[4];
  alignas(32) double area_lane[4];
  alignas(32) double phi_lane[4];
  alignas(32) double phin_lane[4];
  std::size_t pos = begin;
  for (; pos + 4 <= end; pos += 4) {
    const std::int64_t j0 =
        static_cast<std::int64_t>(order[pos]);
    const std::int64_t j1 =
        static_cast<std::int64_t>(order[pos + 1]);
    const std::int64_t j2 =
        static_cast<std::int64_t>(order[pos + 2]);
    const std::int64_t j3 =
        static_cast<std::int64_t>(order[pos + 3]);
    const __m256i idx = _mm256_set_epi64x(j3, j2, j1, j0);
    _mm256_store_pd(sx_lane,
                    _mm256_i64gather_pd(panel_cache.x.data(), idx, 8));
    _mm256_store_pd(sy_lane,
                    _mm256_i64gather_pd(panel_cache.y.data(), idx, 8));
    _mm256_store_pd(sz_lane,
                    _mm256_i64gather_pd(panel_cache.z.data(), idx, 8));
    _mm256_store_pd(snx_lane,
                    _mm256_i64gather_pd(panel_cache.nx.data(), idx, 8));
    _mm256_store_pd(sny_lane,
                    _mm256_i64gather_pd(panel_cache.ny.data(), idx, 8));
    _mm256_store_pd(snz_lane,
                    _mm256_i64gather_pd(panel_cache.nz.data(), idx, 8));
    _mm256_store_pd(area_lane,
                    _mm256_i64gather_pd(panel_cache.area.data(), idx, 8));
    _mm256_store_pd(phi_lane,
                    _mm256_i64gather_pd(unknown_cache.phi.data(), idx, 8));
    _mm256_store_pd(phin_lane,
                    _mm256_i64gather_pd(unknown_cache.phin.data(), idx, 8));
    const __m256d sx = _mm256_load_pd(sx_lane);
    const __m256d sy = _mm256_load_pd(sy_lane);
    const __m256d sz = _mm256_load_pd(sz_lane);
    const __m256d ux = _mm256_sub_pd(tx, sx);
    const __m256d uy = _mm256_sub_pd(ty, sy);
    const __m256d uz = _mm256_sub_pd(tz, sz);
    alignas(32) double ux_lane[4];
    alignas(32) double uy_lane[4];
    alignas(32) double uz_lane[4];
    _mm256_store_pd(ux_lane, ux);
    _mm256_store_pd(uy_lane, uy);
    _mm256_store_pd(uz_lane, uz);
    const std::size_t source_index[4] = {
        static_cast<std::size_t>(j0), static_cast<std::size_t>(j1),
        static_cast<std::size_t>(j2), static_cast<std::size_t>(j3)};
    for (int lane = 0; lane < 4; ++lane) {
      if (source_index[lane] == target_index)
        continue;
      if (accumulate_jcp_leaf_source(
          acc, target.x, target.y, target.z, target.nx, target.ny, target.nz,
          sx_lane[lane], sy_lane[lane], sz_lane[lane], snx_lane[lane],
          sny_lane[lane], snz_lane[lane], area_lane[lane], phi_lane[lane],
          phin_lane[lane], eps_ratio, inv_eps_ratio, poisson_scale,
          inv_four_pi, kk, coupled))
        ++near_count;
    }
  }
  const JcpAccum tail = direct_jcp_ordered_sum_scalar(
      target, target_index, panel_cache, unknown_cache, eps_ratio, kappa,
      coupled, order, pos, end - pos, near_count);
  acc.phi += tail.phi;
  acc.phin += tail.phin;
  return acc;
}
#endif

std::string resolve_leaf_panel_kernel(const std::string &requested) {
  if (requested == "scalar" || requested == "auto")
    return "scalar";
  if (requested == "tile-scalar")
    return "tile-scalar";
#if TITAN_SIMD_X86 && TITAN_RUNTIME_TARGET_ATTRS
  if ((requested == "avx2-gather" || requested == "avx2") &&
      titan_simd_dispatch::runtime_has_avx2())
    return "avx2-gather";
  if ((requested == "tile" || requested == "tile-avx2") &&
      titan_simd_dispatch::runtime_has_avx2())
    return "tile-avx2";
#endif
  if (requested == "tile" || requested == "tile-avx2")
    return "tile-scalar";
  return "scalar";
}

JcpAccum direct_jcp_ordered_sum(const Panel &target, std::size_t target_index,
                                const TreePanelCache &panel_cache,
                                const TreeUnknownCache &unknown_cache,
                                double eps_ratio, double kappa, bool coupled,
                                const std::vector<std::size_t> &order,
                                std::size_t begin, std::size_t count,
                                const std::string &kernel,
                                std::size_t &near_count) noexcept {
#if TITAN_SIMD_X86 && TITAN_RUNTIME_TARGET_ATTRS
  if (kernel == "avx2-gather" || kernel == "avx2") {
    return direct_jcp_ordered_sum_avx2(
        target, target_index, panel_cache, unknown_cache, eps_ratio, kappa,
        coupled, order, begin, count, near_count);
  }
#endif
  return direct_jcp_ordered_sum_scalar(
      target, target_index, panel_cache, unknown_cache, eps_ratio, kappa,
      coupled, order, begin, count, near_count);
}

JcpAccum taylor_jcp_node_sum(const Panel &target, const TreeNode &node,
                             const NodeMoments &phin_scalar,
                             const std::array<NodeMoments, 3> &phi_normal,
                             double eps_ratio, double kappa, bool coupled,
                             int order,
                             const std::string &derivative_mode) noexcept {
  JcpAccum acc;
  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  const double poisson_scale = eps_ratio - 1.0;
  if (derivative_mode != "finite-diff-debug") {
    if (coupled) {
      acc.phi += taylor_sum_series(
          phin_scalar, order,
          series_sub(green_series(target.x, target.y, target.z, node.cx,
                                  node.cy, node.cz, 0.0, order),
                     green_series(target.x, target.y, target.z, node.cx,
                                  node.cy, node.cz, kappa, order),
                     order));
      acc.phin += taylor_sum_series(
          phin_scalar, order,
          series_sub(target_normal_green_series(target, node.cx, node.cy,
                                                node.cz, 0.0, order),
                     series_scale(target_normal_green_series(
                                      target, node.cx, node.cy, node.cz, kappa,
                                      order),
                                  inv_eps_ratio, order),
                     order));
    }
    for (int c = 0; c < 3; ++c) {
      const TaylorSeries d0 = source_normal_component_green_series(
          target.x, target.y, target.z, node.cx, node.cy, node.cz, c, 0.0,
          order);
      if (!coupled) {
        acc.phi += taylor_sum_series(
            phi_normal[static_cast<std::size_t>(c)], order,
            series_scale(d0, poisson_scale, order));
      } else {
        const TaylorSeries dk = source_normal_component_green_series(
            target.x, target.y, target.z, node.cx, node.cy, node.cz, c, kappa,
            order);
        acc.phi += taylor_sum_series(
            phi_normal[static_cast<std::size_t>(c)], order,
            series_sub(series_scale(dk, eps_ratio, order), d0, order));
        acc.phin += taylor_sum_series(
            phi_normal[static_cast<std::size_t>(c)], order,
            series_sub(mixed_normal_component_green_series(
                           target, node.cx, node.cy, node.cz, c, kappa,
                           order),
                       mixed_normal_component_green_series(
                           target, node.cx, node.cy, node.cz, c, 0.0, order),
                       order));
      }
    }
    return acc;
  }

  const double dx = target.x - node.cx;
  const double dy = target.y - node.cy;
  const double dz = target.z - node.cz;
  const double distance =
      std::sqrt(std::max(1.0e-16, dx * dx + dy * dy + dz * dz));
  const double h = std::clamp(0.05 * distance, 1.0e-4,
                              std::max(1.0e-3, 0.20 * node.radius));
  if (coupled) {
    acc.phi += taylor_sum_kernel(node, phin_scalar, order, h,
                                 [&](double sx, double sy, double sz) {
                                   return green_value(target.x, target.y,
                                                      target.z, sx, sy, sz,
                                                      0.0) -
                                          green_value(target.x, target.y,
                                                      target.z, sx, sy, sz,
                                                      kappa);
                                 });
    acc.phin += taylor_sum_kernel(node, phin_scalar, order, h,
                                  [&](double sx, double sy, double sz) {
                                    return target_normal_green(target, sx, sy,
                                                               sz, 0.0) -
                                           inv_eps_ratio *
                                               target_normal_green(
                                                   target, sx, sy, sz, kappa);
                                  });
  }
  for (int c = 0; c < 3; ++c) {
    acc.phi += taylor_sum_kernel(
        node, phi_normal[static_cast<std::size_t>(c)], order, h,
        [&](double sx, double sy, double sz) {
          const double d0 = source_normal_component_green(
              target.x, target.y, target.z, sx, sy, sz, c, 0.0);
          if (!coupled)
            return poisson_scale * d0;
          const double dk = source_normal_component_green(
              target.x, target.y, target.z, sx, sy, sz, c, kappa);
          return eps_ratio * dk - d0;
        });
    if (coupled) {
      acc.phin += taylor_sum_kernel(
          node, phi_normal[static_cast<std::size_t>(c)], order, h,
          [&](double sx, double sy, double sz) {
            return mixed_normal_component_green(target, sx, sy, sz, c, kappa) -
                   mixed_normal_component_green(target, sx, sy, sz, c, 0.0);
          });
    }
  }
  return acc;
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

void include_bound_point(TreeNode &node, double x, double y, double z,
                         bool &initialized) noexcept {
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
    return;
  if (!initialized) {
    node.min_x = node.max_x = x;
    node.min_y = node.max_y = y;
    node.min_z = node.max_z = z;
    initialized = true;
    return;
  }
  node.min_x = std::min(node.min_x, x);
  node.min_y = std::min(node.min_y, y);
  node.min_z = std::min(node.min_z, z);
  node.max_x = std::max(node.max_x, x);
  node.max_y = std::max(node.max_y, y);
  node.max_z = std::max(node.max_z, z);
}

void include_panel_bounds(TreeNode &node, const Panel &panel,
                          bool &initialized) noexcept {
  if (panel.has_triangle_geometry) {
    include_bound_point(node, panel.v0x, panel.v0y, panel.v0z, initialized);
    include_bound_point(node, panel.v1x, panel.v1y, panel.v1z, initialized);
    include_bound_point(node, panel.v2x, panel.v2y, panel.v2z, initialized);
  }
  include_bound_point(node, panel.x, panel.y, panel.z, initialized);
}

void finalize_tight_bounds(TreeNode &node, bool initialized) noexcept {
  if (!initialized) {
    node.min_x = node.max_x = node.cx;
    node.min_y = node.max_y = node.cy;
    node.min_z = node.max_z = node.cz;
  }
  node.cx = 0.5 * (node.min_x + node.max_x);
  node.cy = 0.5 * (node.min_y + node.max_y);
  node.cz = 0.5 * (node.min_z + node.max_z);
  const double hx = 0.5 * (node.max_x - node.min_x);
  const double hy = 0.5 * (node.max_y - node.min_y);
  const double hz = 0.5 * (node.max_z - node.min_z);
  node.radius = std::sqrt(std::max(0.0, hx * hx + hy * hy + hz * hz));
}

void set_moment_bounds_from_radius(TreeNode &node) noexcept {
  node.min_x = node.cx - node.radius;
  node.min_y = node.cy - node.radius;
  node.min_z = node.cz - node.radius;
  node.max_x = node.cx + node.radius;
  node.max_y = node.cy + node.radius;
  node.max_z = node.cz + node.radius;
}

TreeNode make_leaf_node(const std::vector<Panel> &panels,
                        const std::vector<std::size_t> &order,
                        std::size_t first_panel, std::size_t panel_count,
                        bool tight_bounds) noexcept {
  TreeNode node;
  node.first_panel = first_panel;
  node.panel_count = panel_count;
  bool bounds_initialized = false;
  double weight = 0.0;
  for (std::size_t pos = first_panel; pos < first_panel + panel_count; ++pos) {
    const Panel &panel = panels[order[pos]];
    include_panel_bounds(node, panel, bounds_initialized);
    node.max_edge_length =
        std::max(node.max_edge_length, tree_panel_max_edge_length(panel));
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
  if (tight_bounds)
    finalize_tight_bounds(node, bounds_initialized);
  else
    set_moment_bounds_from_radius(node);
  return node;
}

TreeNode make_parent_node(const std::vector<TreeNode> &nodes,
                          const std::vector<std::size_t> &children,
                          bool tight_bounds) noexcept {
  TreeNode node;
  node.child_count = children.size();
  if (!children.empty())
    node.first_panel = nodes[children.front()].first_panel;
  bool bounds_initialized = false;
  double weight = 0.0;
  for (std::size_t child_index : children) {
    const TreeNode &child_node = nodes[child_index];
    node.max_edge_length =
        std::max(node.max_edge_length, child_node.max_edge_length);
    node.depth = std::max(node.depth, child_node.depth + 1);
    include_bound_point(node, child_node.min_x, child_node.min_y,
                        child_node.min_z, bounds_initialized);
    include_bound_point(node, child_node.max_x, child_node.max_y,
                        child_node.max_z, bounds_initialized);
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
  } else if (!children.empty()) {
    node.cx = 0.0;
    node.cy = 0.0;
    node.cz = 0.0;
    for (std::size_t child_index : children) {
      const TreeNode &child_node = nodes[child_index];
      node.cx += child_node.cx;
      node.cy += child_node.cy;
      node.cz += child_node.cz;
    }
    const double inv_count = 1.0 / static_cast<double>(children.size());
    node.cx *= inv_count;
    node.cy *= inv_count;
    node.cz *= inv_count;
  }
  for (std::size_t child_index : children) {
    const TreeNode &child_node = nodes[child_index];
    const double dx = child_node.cx - node.cx;
    const double dy = child_node.cy - node.cy;
    const double dz = child_node.cz - node.cz;
    node.radius =
        std::max(node.radius,
                 std::sqrt(dx * dx + dy * dy + dz * dz) + child_node.radius);
  }
  if (tight_bounds)
    finalize_tight_bounds(node, bounds_initialized);
  else
    set_moment_bounds_from_radius(node);
  return node;
}

std::size_t child_node_index(const TreeWorkspace &tree, const TreeNode &node,
                             std::size_t child) noexcept {
  const std::size_t pos = node.first_child + child;
  if (pos < tree.child_indices.size())
    return tree.child_indices[pos];
  return static_cast<std::size_t>(-1);
}

void attach_children(TreeWorkspace &tree, TreeNode &node,
                     const std::vector<std::size_t> &children) {
  node.first_child = tree.child_indices.size();
  node.child_count = children.size();
  tree.child_indices.insert(tree.child_indices.end(), children.begin(),
                            children.end());
}

struct TopDownTempNode {
  std::size_t first_panel{0};
  std::size_t panel_count{0};
  std::array<std::size_t, kBranchFactor> children{};
  std::size_t child_count{0};
  std::size_t depth{0};
};

void panel_centroid_bounds(const std::vector<Panel> &panels,
                           const std::vector<std::size_t> &order,
                           std::size_t first, std::size_t count,
                           double &min_x, double &min_y, double &min_z,
                           double &max_x, double &max_y,
                           double &max_z) noexcept {
  min_x = min_y = min_z = std::numeric_limits<double>::infinity();
  max_x = max_y = max_z = -std::numeric_limits<double>::infinity();
  for (std::size_t pos = first; pos < first + count; ++pos) {
    const Panel &panel = panels[order[pos]];
    min_x = std::min(min_x, panel.x);
    min_y = std::min(min_y, panel.y);
    min_z = std::min(min_z, panel.z);
    max_x = std::max(max_x, panel.x);
    max_y = std::max(max_y, panel.y);
    max_z = std::max(max_z, panel.z);
  }
  if (!std::isfinite(min_x)) {
    min_x = min_y = min_z = 0.0;
    max_x = max_y = max_z = 0.0;
  }
}

std::size_t build_topdown_temp_node(
    const std::vector<Panel> &panels, std::vector<std::size_t> &order,
    std::size_t first, std::size_t count, std::size_t leaf_size,
    std::size_t depth, std::vector<TopDownTempNode> &temp_nodes) {
  const std::size_t node_index = temp_nodes.size();
  TopDownTempNode node;
  node.first_panel = first;
  node.panel_count = count;
  node.depth = depth;
  temp_nodes.push_back(node);
  if (count <= leaf_size || count <= 1 || depth >= 32)
    return node_index;

  double min_x = 0.0;
  double min_y = 0.0;
  double min_z = 0.0;
  double max_x = 0.0;
  double max_y = 0.0;
  double max_z = 0.0;
  panel_centroid_bounds(panels, order, first, count, min_x, min_y, min_z,
                        max_x, max_y, max_z);
  const double mid_x = 0.5 * (min_x + max_x);
  const double mid_y = 0.5 * (min_y + max_y);
  const double mid_z = 0.5 * (min_z + max_z);
  std::array<std::vector<std::size_t>, kBranchFactor> buckets;
  for (std::size_t pos = first; pos < first + count; ++pos) {
    const Panel &panel = panels[order[pos]];
    std::size_t bucket = 0;
    if (panel.x >= mid_x)
      bucket |= 4u;
    if (panel.y >= mid_y)
      bucket |= 2u;
    if (panel.z >= mid_z)
      bucket |= 1u;
    buckets[bucket].push_back(order[pos]);
  }
  std::size_t nonempty = 0;
  std::size_t largest = 0;
  for (const auto &bucket : buckets) {
    if (!bucket.empty()) {
      ++nonempty;
      largest = std::max(largest, bucket.size());
    }
  }
  if (nonempty <= 1 || largest == count)
    return node_index;

  std::size_t write = first;
  std::array<std::pair<std::size_t, std::size_t>, kBranchFactor> ranges{};
  for (std::size_t bucket = 0; bucket < kBranchFactor; ++bucket) {
    ranges[bucket] = {write, buckets[bucket].size()};
    for (std::size_t panel_index : buckets[bucket])
      order[write++] = panel_index;
  }

  for (std::size_t bucket = 0; bucket < kBranchFactor; ++bucket) {
    const auto [child_first, child_count] = ranges[bucket];
    if (child_count == 0)
      continue;
    const std::size_t child = build_topdown_temp_node(
        panels, order, child_first, child_count, leaf_size, depth + 1,
        temp_nodes);
    TopDownTempNode &stored = temp_nodes[node_index];
    stored.children[stored.child_count++] = child;
  }
  return node_index;
}

void flatten_topdown_tree(const std::vector<Panel> &panels,
                          TreeWorkspace &tree,
                          const std::vector<TopDownTempNode> &temp_nodes,
                          std::size_t temp_root, bool tight_bounds) {
  std::vector<std::size_t> temp_to_final(temp_nodes.size(),
                                         static_cast<std::size_t>(-1));
  std::size_t max_depth = 0;
  for (const TopDownTempNode &node : temp_nodes)
    max_depth = std::max(max_depth, node.depth);

  tree.nodes.clear();
  tree.child_indices.clear();
  tree.moment_levels.clear();
  tree.nodes.reserve(temp_nodes.size());
  for (std::size_t temp = 0; temp < temp_nodes.size(); ++temp) {
    const TopDownTempNode &node = temp_nodes[temp];
    if (node.child_count != 0)
      continue;
    temp_to_final[temp] = tree.nodes.size();
    tree.nodes.push_back(make_leaf_node(panels, tree.order, node.first_panel,
                                        node.panel_count, tight_bounds));
  }
  tree.leaf_count = tree.nodes.size();
  tree.moment_levels.push_back({0, tree.leaf_count});

  for (std::size_t reverse_depth = max_depth + 1; reverse_depth-- > 0;) {
    const std::size_t level_begin = tree.nodes.size();
    for (std::size_t temp = 0; temp < temp_nodes.size(); ++temp) {
      const TopDownTempNode &node = temp_nodes[temp];
      if (node.child_count == 0 || node.depth != reverse_depth)
        continue;
      std::vector<std::size_t> children;
      children.reserve(node.child_count);
      bool mapped = true;
      for (std::size_t child = 0; child < node.child_count; ++child) {
        const std::size_t temp_child = node.children[child];
        if (temp_child >= temp_to_final.size() ||
            temp_to_final[temp_child] == static_cast<std::size_t>(-1)) {
          mapped = false;
          break;
        }
        children.push_back(temp_to_final[temp_child]);
      }
      if (!mapped || children.empty())
        continue;
      TreeNode parent = make_parent_node(tree.nodes, children, tight_bounds);
      parent.first_panel = node.first_panel;
      parent.panel_count = node.panel_count;
      attach_children(tree, parent, children);
      temp_to_final[temp] = tree.nodes.size();
      tree.nodes.push_back(parent);
    }
    if (tree.nodes.size() > level_begin)
      tree.moment_levels.push_back({level_begin, tree.nodes.size()});
    if (reverse_depth == 0)
      break;
  }
  tree.root_index =
      temp_root < temp_to_final.size() ? temp_to_final[temp_root] : 0;
  tree.depth_max = max_depth;
}

bool valid_vertex_index(std::size_t index) noexcept {
  return index != static_cast<std::size_t>(-1);
}

bool panel_has_vertex(const Panel &panel, std::size_t vertex) noexcept {
  if (!valid_vertex_index(vertex))
    return false;
  return panel.v0_index == vertex || panel.v1_index == vertex ||
         panel.v2_index == vertex;
}

int shared_vertex_count_treecode(const Panel &a, const Panel &b) noexcept {
  if (!a.has_triangle_geometry || !b.has_triangle_geometry)
    return 0;
  int count = 0;
  if (panel_has_vertex(b, a.v0_index))
    ++count;
  if (a.v1_index != a.v0_index && panel_has_vertex(b, a.v1_index))
    ++count;
  if (a.v2_index != a.v0_index && a.v2_index != a.v1_index &&
      panel_has_vertex(b, a.v2_index))
    ++count;
  return count;
}

bool is_singular_near_pair(const Panel &target, std::size_t target_index,
                           const Panel &source,
                           std::size_t source_index) noexcept {
  if (source_index == target_index)
    return true;
  return shared_vertex_count_treecode(target, source) > 0;
}

std::string normalize_leaf_near_policy(const std::string &policy) {
  if (policy == "adaptive-large")
    return "adaptive-large";
  if (policy == "large-protein")
    return "large-protein";
  if (policy == "smooth-node")
    return "smooth-node";
  return "exact";
}

std::string normalize_leaf_panel_kernel(const std::string &kernel) {
  if (kernel == "tile" || kernel == "tile-scalar" || kernel == "tile-avx2")
    return kernel;
  if (kernel == "avx2-gather" || kernel == "avx2")
    return "avx2-gather";
  return "scalar";
}

bool is_tile_leaf_panel_kernel(const std::string &kernel) noexcept {
  return kernel == "tile" || kernel == "tile-scalar" ||
         kernel == "tile-avx2";
}

bool is_effective_tile_leaf_panel_kernel(const std::string &kernel) noexcept {
  return kernel == "tile-scalar" || kernel == "tile-avx2";
}

struct AdaptivePruneParams {
  double smooth_theta_cap{0.24};
  double smooth_theta_factor{0.30};
  double p4_proxy_limit{0.0035};
};

AdaptivePruneParams adaptive_prune_params(const std::string &level) noexcept {
  if (level == "aggressive")
    return {0.45, 0.60, 0.0400};
  if (level == "balanced")
    return {0.34, 0.45, 0.0120};
  return {};
}

std::vector<std::pair<std::size_t, std::size_t>>
make_leaf_spans(const TreeNode &node, const std::vector<std::size_t> &order,
                const std::vector<std::size_t> &excluded_original) {
  std::vector<std::pair<std::size_t, std::size_t>> spans;
  const std::size_t end = node.first_panel + node.panel_count;
  std::size_t span_begin = end;
  std::size_t span_count = 0;
  auto is_excluded = [&](std::size_t original) noexcept {
    for (std::size_t value : excluded_original) {
      if (value == original)
        return true;
    }
    return false;
  };
  for (std::size_t pos = node.first_panel; pos < end; ++pos) {
    const std::size_t original = pos < order.size() ? order[pos] : pos;
    if (is_excluded(original)) {
      if (span_count > 0)
        spans.push_back({span_begin, span_count});
      span_begin = end;
      span_count = 0;
      continue;
    }
    if (span_count == 0)
      span_begin = pos;
    ++span_count;
  }
  if (span_count > 0)
    spans.push_back({span_begin, span_count});
  return spans;
}

std::string tile_fallback_status(const std::string &kernel,
                                 std::size_t tile_nodes) {
  if (tile_nodes > 0)
    return "active";
  if (is_effective_tile_leaf_panel_kernel(kernel))
    return "fallback";
  return "off";
}

TreeInteractionLists build_tree_interactions(const std::vector<Panel> &panels,
                                             const TreeWorkspace &tree,
                                             TitanExecutor &executor) {
  const auto begin_time = std::chrono::steady_clock::now();
  TreeInteractionLists lists;
  const std::size_t n = panels.size();
  lists.singular_offsets.assign(n + 1, 0);
  lists.leaf_offsets.assign(n + 1, 0);
  lists.leaf_tile_offsets.assign(n + 1, 0);
  lists.smooth_leaf_offsets.assign(n + 1, 0);
  lists.far_offsets.assign(n + 1, 0);
  if (n == 0 || tree.nodes.empty() || tree.order.size() != n) {
    lists.build_ms =
        elapsed_ms(begin_time, std::chrono::steady_clock::now());
    return lists;
  }

  std::vector<std::size_t> rank(n, 0);
  for (std::size_t pos = 0; pos < tree.order.size(); ++pos)
    rank[tree.order[pos]] = pos;

  std::vector<std::vector<std::size_t>> singular_by_target(n);
  std::vector<std::vector<std::size_t>> leaf_by_target(n);
  std::vector<std::vector<std::size_t>> leaf_tile_by_target(n);
  std::vector<std::vector<std::size_t>> leaf_tile_span_counts_by_target(n);
  std::vector<std::vector<std::size_t>> leaf_tile_span_firsts_by_target(n);
  std::vector<std::vector<std::size_t>> leaf_tile_span_lengths_by_target(n);
  std::vector<std::vector<std::size_t>> smooth_leaf_by_target(n);
  std::vector<std::vector<std::size_t>> far_by_target(n);
  std::vector<std::size_t> leaf_pruned_by_target(n, 0);
  std::vector<std::size_t> adaptive_smooth_accepted_by_target(n, 0);
  std::vector<std::size_t> adaptive_smooth_rejected_by_target(n, 0);
  std::vector<std::size_t> adaptive_span_accepted_by_target(n, 0);
  std::vector<std::size_t> adaptive_span_rejected_by_target(n, 0);
  std::vector<std::size_t> adaptive_exact_fallbacks_by_target(n, 0);
  std::vector<std::size_t> adaptive_reject_one_ring_by_target(n, 0);
  std::vector<std::size_t> adaptive_reject_error_by_target(n, 0);
  std::vector<std::size_t> adaptive_reject_p4_by_target(n, 0);
  std::vector<std::size_t> adaptive_reject_too_close_by_target(n, 0);
  std::vector<std::size_t> target_cluster_interactions_by_target(n, 0);
  std::vector<std::size_t> target_cluster_targets_by_target(n, 0);
  std::vector<std::size_t> target_cluster_panels_by_target(n, 0);
  std::vector<std::size_t> target_cluster_rejected_by_target(n, 0);
  std::vector<std::size_t> mac_reject_by_target(n, 0);
  std::vector<MacAuditTally> mac_audit_by_target(n);
  const AdaptivePruneParams adaptive_params =
      adaptive_prune_params(tree.adaptive_level);
  const std::size_t workers = effective_workers(executor, n);
  executor.parallel_for(
      n, workers, 8,
      [&](std::size_t, std::size_t target_begin, std::size_t target_end) {
        std::vector<std::size_t> stack;
        stack.reserve(64);
        for (std::size_t i = target_begin; i < target_end; ++i) {
          const Panel &target = panels[i];
          const std::size_t target_rank = rank[i];
          auto &singular = singular_by_target[i];
          auto &leaf = leaf_by_target[i];
          auto &leaf_tile = leaf_tile_by_target[i];
          auto &leaf_tile_span_counts = leaf_tile_span_counts_by_target[i];
          auto &leaf_tile_span_firsts = leaf_tile_span_firsts_by_target[i];
          auto &leaf_tile_span_lengths = leaf_tile_span_lengths_by_target[i];
          auto &smooth_leaf = smooth_leaf_by_target[i];
          auto &far = far_by_target[i];
          stack.clear();
          stack.push_back(tree.root_index);
          while (!stack.empty()) {
            const std::size_t node_index = stack.back();
            stack.pop_back();
            const TreeNode &node = tree.nodes[node_index];
            const bool is_leaf = node.child_count == 0;
            const bool contains_target =
                target_rank >= node.first_panel &&
                target_rank < node.first_panel + node.panel_count;
            if (!is_leaf && !contains_target) {
              const double dx = target.x - node.cx;
              const double dy = target.y - node.cy;
              const double dz = target.z - node.cz;
              const double distance =
                  std::sqrt(std::max(1.0e-16, dx * dx + dy * dy + dz * dz));
              const MacDecision mac = evaluate_mac(target, node, tree, distance);
              if (mac.ratio_accept && tree.mac_safety != "off") {
                record_far_geometric_audit(mac_audit_by_target[i], node,
                                           mac.gap, mac.gap_ratio,
                                           mac.unsafe);
                if (!mac.accepted)
                  ++mac_audit_by_target[i].physical_gap_rejects;
              }
              if (mac.accepted) {
                far.push_back(node_index);
                continue;
              }
              ++mac_reject_by_target[i];
            }
            if (is_leaf) {
              bool has_singular_panel = false;
              std::vector<std::size_t> singular_candidates;
              std::vector<std::size_t> residual_candidates;
              singular_candidates.reserve(4);
              residual_candidates.reserve(node.panel_count);
              for (std::size_t pos = node.first_panel;
                   pos < node.first_panel + node.panel_count; ++pos) {
                const std::size_t source_index = tree.order[pos];
                const Panel &source = panels[source_index];
                if (is_singular_near_pair(target, i, source, source_index) ||
                    classify_geometric_near_buffer_pair(
                        target, source, tree.quad_buffer)
                        .active) {
                  has_singular_panel = true;
                  singular_candidates.push_back(source_index);
                } else {
                  residual_candidates.push_back(source_index);
                }
              }
              bool accepts_smooth_leaf = false;
              const bool adaptive_policy =
                  tree.leaf_near_policy == "adaptive-large";
              if (tree.leaf_near_policy == "smooth-node" &&
                  !contains_target && node.panel_count > 0 &&
                  !has_singular_panel) {
                  const double dx = target.x - node.cx;
                  const double dy = target.y - node.cy;
                  const double dz = target.z - node.cz;
                  const double distance = std::sqrt(
                      std::max(1.0e-16, dx * dx + dy * dy + dz * dz));
                  const double smooth_opening_theta =
                      std::min(0.35, 0.5 * tree.opening_theta);
                  accepts_smooth_leaf =
                      node.radius / distance < smooth_opening_theta;
              }
              if (adaptive_policy && !contains_target && node.panel_count > 0) {
                if (has_singular_panel) {
                  ++adaptive_smooth_rejected_by_target[i];
                  ++adaptive_reject_one_ring_by_target[i];
                } else {
                  const double dx = target.x - node.cx;
                  const double dy = target.y - node.cy;
                  const double dz = target.z - node.cz;
                  const double distance = std::sqrt(
                      std::max(1.0e-16, dx * dx + dy * dy + dz * dz));
                  const double ratio = node.radius / distance;
                  const double p4_proxy = ratio * ratio * ratio * ratio;
                  const double kAdaptiveSmoothTheta =
                      std::min(adaptive_params.smooth_theta_cap,
                               adaptive_params.smooth_theta_factor *
                                   tree.opening_theta);
                  const double kAdaptiveP4ProxyLimit =
                      adaptive_params.p4_proxy_limit;
                  if (ratio < kAdaptiveSmoothTheta &&
                      p4_proxy < kAdaptiveP4ProxyLimit) {
                    accepts_smooth_leaf = true;
                    ++adaptive_smooth_accepted_by_target[i];
                    ++target_cluster_interactions_by_target[i];
                    ++target_cluster_targets_by_target[i];
                    target_cluster_panels_by_target[i] += node.panel_count;
                  } else {
                    ++adaptive_smooth_rejected_by_target[i];
                    ++target_cluster_rejected_by_target[i];
                    if (ratio >= 0.5) {
                      ++adaptive_reject_too_close_by_target[i];
                    } else if (p4_proxy >= kAdaptiveP4ProxyLimit) {
                      ++adaptive_reject_p4_by_target[i];
                    } else {
                      ++adaptive_reject_error_by_target[i];
                    }
                  }
                }
              }
              if (accepts_smooth_leaf) {
                smooth_leaf.push_back(node_index);
                leaf_pruned_by_target[i] += node.panel_count;
                continue;
              }
              if ((tree.leaf_near_policy == "large-protein" ||
                   adaptive_policy) &&
                  is_tile_leaf_panel_kernel(tree.leaf_panel_kernel) &&
                  !residual_candidates.empty() && node.panel_count >= 16 &&
                  singular_candidates.size() * 4 <= node.panel_count) {
                const auto spans =
                    make_leaf_spans(node, tree.order, singular_candidates);
                std::size_t total_span_panels = 0;
                for (const auto &span : spans)
                  total_span_panels += span.second;
                const bool span_ok =
                    adaptive_policy
                        ? (spans.size() <= 4 && total_span_panels >= 64 &&
                           total_span_panels >= 32 * spans.size())
                        : (spans.size() <= 4 && total_span_panels >= 16);
                if (span_ok) {
                  if (adaptive_policy)
                    ++adaptive_span_accepted_by_target[i];
                  singular.insert(singular.end(), singular_candidates.begin(),
                                  singular_candidates.end());
                  leaf_tile.push_back(node_index);
                  leaf_tile_span_counts.push_back(spans.size());
                  for (const auto &span : spans) {
                    leaf_tile_span_firsts.push_back(span.first);
                    leaf_tile_span_lengths.push_back(span.second);
                  }
                  continue;
                }
                if (adaptive_policy)
                  ++adaptive_span_rejected_by_target[i];
              }
              if (adaptive_policy)
                ++adaptive_exact_fallbacks_by_target[i];
              if (is_tile_leaf_panel_kernel(tree.leaf_panel_kernel) &&
                  !contains_target && !has_singular_panel &&
                  node.panel_count > 0) {
                leaf_tile.push_back(node_index);
                leaf_tile_span_counts.push_back(1);
                leaf_tile_span_firsts.push_back(node.first_panel);
                leaf_tile_span_lengths.push_back(node.panel_count);
                continue;
              }
              singular.insert(singular.end(), singular_candidates.begin(),
                              singular_candidates.end());
              leaf.insert(leaf.end(), residual_candidates.begin(),
                          residual_candidates.end());
            } else {
              for (std::size_t child = 0; child < node.child_count; ++child) {
                const std::size_t child_index =
                    child_node_index(tree, node, child);
                if (child_index < tree.nodes.size())
                  stack.push_back(child_index);
              }
            }
          }
        }
      });

  for (std::size_t i = 0; i < n; ++i) {
    lists.singular_offsets[i + 1] =
        lists.singular_offsets[i] + singular_by_target[i].size();
    lists.leaf_offsets[i + 1] =
        lists.leaf_offsets[i] + leaf_by_target[i].size();
    lists.leaf_tile_offsets[i + 1] =
        lists.leaf_tile_offsets[i] + leaf_tile_by_target[i].size();
    lists.smooth_leaf_offsets[i + 1] =
        lists.smooth_leaf_offsets[i] + smooth_leaf_by_target[i].size();
    lists.far_offsets[i + 1] = lists.far_offsets[i] + far_by_target[i].size();
    lists.leaf_pruned_count += leaf_pruned_by_target[i];
    lists.adaptive_smooth_leaf_accepted +=
        adaptive_smooth_accepted_by_target[i];
    lists.adaptive_smooth_leaf_rejected +=
        adaptive_smooth_rejected_by_target[i];
    lists.adaptive_span_tile_accepted += adaptive_span_accepted_by_target[i];
    lists.adaptive_span_tile_rejected += adaptive_span_rejected_by_target[i];
    lists.adaptive_exact_leaf_fallbacks += adaptive_exact_fallbacks_by_target[i];
    lists.adaptive_reject_one_ring_overlap +=
        adaptive_reject_one_ring_by_target[i];
    lists.adaptive_reject_error_bound_failed +=
        adaptive_reject_error_by_target[i];
    lists.adaptive_reject_p4_check_failed += adaptive_reject_p4_by_target[i];
    lists.adaptive_reject_too_close += adaptive_reject_too_close_by_target[i];
    lists.target_cluster_interactions +=
        target_cluster_interactions_by_target[i];
    lists.target_cluster_targets += target_cluster_targets_by_target[i];
    lists.target_cluster_panels_pruned += target_cluster_panels_by_target[i];
    lists.target_cluster_rejected += target_cluster_rejected_by_target[i];
    lists.mac_reject_radius_limited += mac_reject_by_target[i];
    MacAuditTally tally;
    tally.unsafe_count = lists.far_geometric_unsafe_count;
    tally.unsafe_panels = lists.far_geometric_unsafe_panels;
    tally.min_gap = lists.far_geometric_min_gap;
    tally.min_gap_ratio = lists.far_geometric_gap_ratio_min;
    tally.sum_gap_ratio = lists.far_geometric_gap_ratio_sum;
    tally.samples = lists.far_geometric_audit_samples;
    tally.physical_gap_rejects = lists.mac_physical_gap_rejects;
    merge_far_geometric_audit(tally, mac_audit_by_target[i]);
    lists.far_geometric_unsafe_count = tally.unsafe_count;
    lists.far_geometric_unsafe_panels = tally.unsafe_panels;
    lists.far_geometric_min_gap = tally.min_gap;
    lists.far_geometric_gap_ratio_min = tally.min_gap_ratio;
    lists.far_geometric_gap_ratio_sum = tally.sum_gap_ratio;
    lists.far_geometric_audit_samples = tally.samples;
    lists.mac_physical_gap_rejects = tally.physical_gap_rejects;
  }
  lists.singular_indices.reserve(lists.singular_offsets.back());
  lists.leaf_indices.reserve(lists.leaf_offsets.back());
  lists.leaf_tile_node_indices.reserve(lists.leaf_tile_offsets.back());
  lists.smooth_leaf_node_indices.reserve(lists.smooth_leaf_offsets.back());
  lists.far_node_indices.reserve(lists.far_offsets.back());
  lists.leaf_tile_exclusion_offsets.clear();
  lists.leaf_tile_exclusion_offsets.reserve(lists.leaf_tile_offsets.back() + 1);
  lists.leaf_tile_exclusion_offsets.push_back(0);
  lists.leaf_tile_exclusion_indices.clear();
  lists.leaf_tile_span_offsets.clear();
  lists.leaf_tile_span_offsets.reserve(lists.leaf_tile_offsets.back() + 1);
  lists.leaf_tile_span_offsets.push_back(0);
  lists.leaf_tile_span_firsts.clear();
  lists.leaf_tile_span_counts.clear();
  for (std::size_t i = 0; i < n; ++i) {
    lists.singular_indices.insert(lists.singular_indices.end(),
                                  singular_by_target[i].begin(),
                                  singular_by_target[i].end());
    lists.leaf_indices.insert(lists.leaf_indices.end(),
                              leaf_by_target[i].begin(),
                              leaf_by_target[i].end());
    lists.leaf_tile_node_indices.insert(
        lists.leaf_tile_node_indices.end(), leaf_tile_by_target[i].begin(),
        leaf_tile_by_target[i].end());
    const auto &span_counts = leaf_tile_span_counts_by_target[i];
    const auto &span_firsts = leaf_tile_span_firsts_by_target[i];
    const auto &span_lengths = leaf_tile_span_lengths_by_target[i];
    std::size_t span_pos = 0;
    for (std::size_t entry = 0; entry < leaf_tile_by_target[i].size();
         ++entry) {
      lists.leaf_tile_exclusion_offsets.push_back(
          lists.leaf_tile_exclusion_indices.size());
      const std::size_t count =
          entry < span_counts.size() ? span_counts[entry] : 0;
      if (count > 0 && span_pos + count <= span_firsts.size() &&
          span_pos + count <= span_lengths.size()) {
        lists.leaf_tile_span_firsts.insert(
            lists.leaf_tile_span_firsts.end(),
            span_firsts.begin() + static_cast<std::ptrdiff_t>(span_pos),
            span_firsts.begin() + static_cast<std::ptrdiff_t>(span_pos + count));
        lists.leaf_tile_span_counts.insert(
            lists.leaf_tile_span_counts.end(),
            span_lengths.begin() + static_cast<std::ptrdiff_t>(span_pos),
            span_lengths.begin() +
                static_cast<std::ptrdiff_t>(span_pos + count));
      }
      span_pos += count;
      lists.leaf_tile_span_offsets.push_back(
          lists.leaf_tile_span_firsts.size());
    }
    lists.smooth_leaf_node_indices.insert(
        lists.smooth_leaf_node_indices.end(),
        smooth_leaf_by_target[i].begin(), smooth_leaf_by_target[i].end());
    lists.far_node_indices.insert(lists.far_node_indices.end(),
                                  far_by_target[i].begin(),
                                  far_by_target[i].end());
  }
  lists.singular_pair_count = lists.singular_indices.size();
  lists.leaf_pair_count = lists.leaf_indices.size();
  lists.leaf_tile_node_count = lists.leaf_tile_node_indices.size();
  lists.leaf_tile_panel_count = 0;
  for (std::size_t entry = 0; entry < lists.leaf_tile_node_indices.size();
       ++entry) {
    const std::size_t node_index = lists.leaf_tile_node_indices[entry];
    if (node_index < tree.nodes.size()) {
      const std::size_t span_begin =
          entry + 1 < lists.leaf_tile_span_offsets.size()
              ? lists.leaf_tile_span_offsets[entry]
              : 0;
      const std::size_t span_end =
          entry + 1 < lists.leaf_tile_span_offsets.size()
              ? lists.leaf_tile_span_offsets[entry + 1]
              : span_begin;
      for (std::size_t span = span_begin; span < span_end; ++span) {
        if (span < lists.leaf_tile_span_counts.size())
          lists.leaf_tile_panel_count += lists.leaf_tile_span_counts[span];
      }
    }
  }
  lists.leaf_pair_count += lists.leaf_tile_panel_count;
  lists.smooth_leaf_node_count = lists.smooth_leaf_node_indices.size();
  lists.far_node_count = lists.far_node_indices.size();
  lists.build_ms = elapsed_ms(begin_time, std::chrono::steady_clock::now());
  return lists;
}

void accumulate_leaf_moment(NodeMoments &moment, const Panel &panel,
                            const TreeNode &node, double weight,
                            const std::vector<MomentIndex> &indices) {
  const double dx = panel.x - node.cx;
  const double dy = panel.y - node.cy;
  const double dz = panel.z - node.cz;
  for (std::size_t m = 0; m < indices.size(); ++m) {
    const MomentIndex &index = indices[m];
    moment.values[m] += weight * pow_int(dx, index.ax) *
                        pow_int(dy, index.ay) * pow_int(dz, index.az);
  }
}

void translate_child_moment(NodeMoments &parent, const NodeMoments &child,
                            const TreeNode &parent_node,
                            const TreeNode &child_node,
                            const std::vector<MomentIndex> &indices,
                            const std::array<int, kSeriesSize> &lookup) {
  const double ox = child_node.cx - parent_node.cx;
  const double oy = child_node.cy - parent_node.cy;
  const double oz = child_node.cz - parent_node.cz;
  for (std::size_t alpha_index = 0; alpha_index < indices.size();
       ++alpha_index) {
    const MomentIndex &alpha = indices[alpha_index];
    double translated = 0.0;
    for (int bx = 0; bx <= alpha.ax; ++bx) {
      for (int by = 0; by <= alpha.ay; ++by) {
        for (int bz = 0; bz <= alpha.az; ++bz) {
          const int beta_index =
              lookup[static_cast<std::size_t>(moment_lookup_index(bx, by, bz))];
          if (beta_index < 0)
            continue;
          const double coeff =
              binomial(alpha.ax, bx) * binomial(alpha.ay, by) *
              binomial(alpha.az, bz) * pow_int(ox, alpha.ax - bx) *
              pow_int(oy, alpha.ay - by) * pow_int(oz, alpha.az - bz);
          translated +=
              coeff * child.values[static_cast<std::size_t>(beta_index)];
        }
      }
    }
    parent.values[alpha_index] += translated;
  }
}

std::vector<NodeMoments>
compute_node_moments(const std::vector<Panel> &panels,
                     const std::vector<double> &density,
                     const TreeWorkspace &tree, int taylor_order,
                     TitanExecutor &executor) {
  std::vector<NodeMoments> moments(tree.nodes.size());
  const auto &indices = moment_indices(taylor_order);
  const auto &lookup = moment_lookup();
  if (tree.moment_levels.empty()) {
    for (std::size_t node_index = 0; node_index < tree.nodes.size();
         ++node_index) {
      const TreeNode &node = tree.nodes[node_index];
      if (node.child_count == 0) {
        for (std::size_t pos = node.first_panel;
             pos < node.first_panel + node.panel_count; ++pos) {
          const std::size_t panel_index = tree.order[pos];
          const Panel &panel = panels[panel_index];
          accumulate_leaf_moment(moments[node_index], panel, node,
                                 panel.area * density[panel_index], indices);
        }
      } else {
        for (std::size_t child = 0; child < node.child_count; ++child) {
          const std::size_t child_index = child_node_index(tree, node, child);
          if (child_index >= tree.nodes.size())
            continue;
          translate_child_moment(moments[node_index], moments[child_index],
                                 node, tree.nodes[child_index], indices,
                                 lookup);
        }
      }
    }
    return moments;
  }

  const auto leaf_level = tree.moment_levels.front();
  const std::size_t leaf_count = leaf_level.second - leaf_level.first;
  const std::size_t leaf_workers = effective_workers(executor, leaf_count);
  executor.parallel_for(
      leaf_count, leaf_workers, 16,
      [&](std::size_t, std::size_t begin, std::size_t end) {
        for (std::size_t offset = begin; offset < end; ++offset) {
          const std::size_t node_index = leaf_level.first + offset;
          const TreeNode &node = tree.nodes[node_index];
          for (std::size_t pos = node.first_panel;
               pos < node.first_panel + node.panel_count; ++pos) {
            const std::size_t panel_index = tree.order[pos];
            const Panel &panel = panels[panel_index];
            accumulate_leaf_moment(moments[node_index], panel, node,
                                   panel.area * density[panel_index], indices);
          }
        }
      });

  for (std::size_t level = 1; level < tree.moment_levels.size(); ++level) {
    const auto range = tree.moment_levels[level];
    const std::size_t count = range.second - range.first;
    const std::size_t workers = effective_workers(executor, count);
    executor.parallel_for(
        count, workers, 16,
        [&](std::size_t, std::size_t begin, std::size_t end) {
          for (std::size_t offset = begin; offset < end; ++offset) {
            const std::size_t node_index = range.first + offset;
            const TreeNode &node = tree.nodes[node_index];
            for (std::size_t child = 0; child < node.child_count; ++child) {
              const std::size_t child_index =
                  child_node_index(tree, node, child);
              if (child_index >= tree.nodes.size())
                continue;
              translate_child_moment(moments[node_index],
                                     moments[child_index], node,
                                     tree.nodes[child_index], indices, lookup);
            }
          }
        });
  }
  return moments;
}

struct PackedJcpMoments {
  std::vector<NodeMoments> phin;
  std::array<std::vector<NodeMoments>, 3> phi_normal;
};

PackedJcpMoments compute_jcp_node_moments_packed(
    const std::vector<Panel> &panels, const std::vector<double> &unknowns,
    bool coupled, const TreeWorkspace &tree, int taylor_order,
    TitanExecutor &executor) {
  PackedJcpMoments packed;
  packed.phin.assign(tree.nodes.size(), NodeMoments{});
  for (auto &channel : packed.phi_normal)
    channel.assign(tree.nodes.size(), NodeMoments{});
  const auto &indices = moment_indices(taylor_order);
  const auto &lookup = moment_lookup();
  if (tree.moment_levels.empty())
    return packed;

  const auto leaf_level = tree.moment_levels.front();
  const std::size_t n = panels.size();
  const std::size_t leaf_count = leaf_level.second - leaf_level.first;
  const std::size_t leaf_workers = effective_workers(executor, leaf_count);
  executor.parallel_for(
      leaf_count, leaf_workers, 16,
      [&](std::size_t, std::size_t begin, std::size_t end) {
        for (std::size_t offset = begin; offset < end; ++offset) {
          const std::size_t node_index = leaf_level.first + offset;
          const TreeNode &node = tree.nodes[node_index];
          for (std::size_t pos = node.first_panel;
               pos < node.first_panel + node.panel_count; ++pos) {
            const std::size_t panel_index = tree.order[pos];
            const Panel &panel = panels[panel_index];
            const double phi = panel_index < unknowns.size()
                                   ? unknowns[panel_index]
                                   : 0.0;
            const double area_phi = panel.area * phi;
            accumulate_leaf_moment(packed.phi_normal[0][node_index], panel,
                                   node, area_phi * panel.nx, indices);
            accumulate_leaf_moment(packed.phi_normal[1][node_index], panel,
                                   node, area_phi * panel.ny, indices);
            accumulate_leaf_moment(packed.phi_normal[2][node_index], panel,
                                   node, area_phi * panel.nz, indices);
            if (coupled && n + panel_index < unknowns.size()) {
              accumulate_leaf_moment(packed.phin[node_index], panel, node,
                                     panel.area * unknowns[n + panel_index],
                                     indices);
            }
          }
        }
      });

  for (std::size_t level = 1; level < tree.moment_levels.size(); ++level) {
    const auto range = tree.moment_levels[level];
    const std::size_t count = range.second - range.first;
    const std::size_t workers = effective_workers(executor, count);
    executor.parallel_for(
        count, workers, 16,
        [&](std::size_t, std::size_t begin, std::size_t end) {
          for (std::size_t offset = begin; offset < end; ++offset) {
            const std::size_t node_index = range.first + offset;
            const TreeNode &node = tree.nodes[node_index];
            for (std::size_t child = 0; child < node.child_count; ++child) {
              const std::size_t child_index =
                  child_node_index(tree, node, child);
              if (child_index >= tree.nodes.size())
                continue;
              const TreeNode &child_node = tree.nodes[child_index];
              translate_child_moment(packed.phin[node_index],
                                     packed.phin[child_index], node,
                                     child_node, indices, lookup);
              for (int c = 0; c < 3; ++c) {
                translate_child_moment(
                    packed.phi_normal[static_cast<std::size_t>(c)]
                                      [node_index],
                    packed.phi_normal[static_cast<std::size_t>(c)]
                                      [child_index],
                    node, child_node, indices, lookup);
              }
            }
          }
        });
  }
  return packed;
}

struct MomentAuditSummary {
  std::size_t nodes{0};
  double max_rel_error{0.0};
  double mean_rel_error{0.0};
};

double moment_channel_rel_error(const NodeMoments &direct,
                                const NodeMoments &translated,
                                std::size_t count) noexcept {
  double diff = 0.0;
  double denom = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    diff = std::max(diff, std::abs(direct.values[i] - translated.values[i]));
    denom = std::max(denom, std::abs(direct.values[i]));
  }
  return diff / std::max(denom, 1.0e-30);
}

MomentAuditSummary audit_jcp_moment_translation(
    const std::vector<Panel> &panels, const std::vector<double> &unknowns,
    bool coupled, const TreeWorkspace &tree, int taylor_order,
    const PackedJcpMoments &packed) {
  MomentAuditSummary summary;
  if (tree.moment_audit == "off" || tree.nodes.empty())
    return summary;
  const auto &indices = moment_indices(taylor_order);
  if (indices.empty())
    return summary;
  const std::size_t moment_count = indices.size();
  const std::size_t n = panels.size();
  const bool full =
      tree.moment_audit == "full-small" && tree.nodes.size() <= 4096;
  const std::size_t sample_cap = full ? tree.nodes.size() : 64;
  const std::size_t stride =
      full ? 1 : std::max<std::size_t>(1, tree.nodes.size() / sample_cap);
  double rel_sum = 0.0;
  for (std::size_t node_index = 0; node_index < tree.nodes.size();
       node_index += stride) {
    if (!full && summary.nodes >= sample_cap)
      break;
    const TreeNode &node = tree.nodes[node_index];
    if (node.panel_count == 0)
      continue;
    NodeMoments direct_phin;
    std::array<NodeMoments, 3> direct_phi_normal{};
    const std::size_t end =
        std::min(tree.order.size(), node.first_panel + node.panel_count);
    for (std::size_t pos = node.first_panel; pos < end; ++pos) {
      const std::size_t panel_index = tree.order[pos];
      if (panel_index >= panels.size())
        continue;
      const Panel &panel = panels[panel_index];
      const double phi = panel_index < unknowns.size() ? unknowns[panel_index]
                                                       : 0.0;
      const double area_phi = panel.area * phi;
      accumulate_leaf_moment(direct_phi_normal[0], panel, node,
                             area_phi * panel.nx, indices);
      accumulate_leaf_moment(direct_phi_normal[1], panel, node,
                             area_phi * panel.ny, indices);
      accumulate_leaf_moment(direct_phi_normal[2], panel, node,
                             area_phi * panel.nz, indices);
      if (coupled && n + panel_index < unknowns.size()) {
        accumulate_leaf_moment(direct_phin, panel, node,
                               panel.area * unknowns[n + panel_index],
                               indices);
      }
    }
    double rel = 0.0;
    if (node_index < packed.phin.size())
      rel = std::max(rel, moment_channel_rel_error(
                              direct_phin, packed.phin[node_index],
                              moment_count));
    for (int c = 0; c < 3; ++c) {
      const std::size_t channel = static_cast<std::size_t>(c);
      if (node_index < packed.phi_normal[channel].size()) {
        rel = std::max(rel, moment_channel_rel_error(
                                direct_phi_normal[channel],
                                packed.phi_normal[channel][node_index],
                                moment_count));
      }
    }
    summary.max_rel_error = std::max(summary.max_rel_error, rel);
    rel_sum += rel;
    ++summary.nodes;
  }
  if (summary.nodes > 0)
    summary.mean_rel_error = rel_sum / static_cast<double>(summary.nodes);
  return summary;
}

} // namespace

TreeWorkspace build_tree_workspace(const std::vector<Panel> &panels,
                                   TitanExecutor &executor,
                                   std::size_t leaf_size,
                                   double opening_theta,
                                   const std::string &leaf_near_policy,
                                   const std::string &leaf_panel_kernel,
                                   const std::string &adaptive_level,
                                   const std::string &quad_buffer,
                                   const std::string &interaction_mode,
                                   const std::string &node_bounds_mode,
                                   const std::string &tree_builder,
                                   const std::string &mac_safety,
                                   double mac_gap_factor,
                                   const std::string &moment_audit) {
  const auto build_begin = std::chrono::steady_clock::now();
  TreeWorkspace tree;
  tree.taylor_order = kTaylorOrderDefault;
  const std::size_t effective_leaf_size =
      std::clamp<std::size_t>(leaf_size == 0 ? kDefaultLeafSize : leaf_size, 8,
                              512);
  tree.leaf_size = effective_leaf_size;
  tree.opening_theta =
      std::clamp(std::isfinite(opening_theta) && opening_theta > 0.0
                     ? opening_theta
                     : kDefaultOpeningTheta,
                 0.05, 2.0);
  tree.leaf_near_policy = normalize_leaf_near_policy(leaf_near_policy);
  tree.leaf_panel_kernel = normalize_leaf_panel_kernel(leaf_panel_kernel);
  tree.adaptive_level =
      (adaptive_level == "balanced" || adaptive_level == "aggressive")
          ? adaptive_level
          : std::string("strict");
  tree.interaction_mode =
      interaction_mode == "onthefly" ? std::string("onthefly")
                                     : std::string("cached");
  tree.node_bounds_mode =
      node_bounds_mode == "tight-aabb" ? std::string("tight-aabb")
                                       : std::string("moment");
  tree.tree_builder = tree_builder == "topdown-aabb"
                          ? std::string("topdown-aabb")
                          : std::string("grouped-morton");
  tree.mac_safety =
      (mac_safety == "audit" || mac_safety == "physical-gap")
          ? mac_safety
          : std::string("off");
  tree.mac_gap_factor =
      std::clamp(std::isfinite(mac_gap_factor) ? mac_gap_factor : 1.5, 0.5,
                 4.0);
  tree.moment_audit =
      (moment_audit == "sampled" || moment_audit == "full-small")
          ? moment_audit
          : std::string("off");
  tree.quad_buffer =
      (quad_buffer == "gauss7" || quad_buffer == "adaptive") ? quad_buffer
                                                             : "off";
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

  std::vector<std::uint64_t> panel_morton(panels.size(), 0);
  for (const MortonEntry &entry : entries)
    panel_morton[entry.index] = entry.code;
  tree.order.reserve(entries.size());
  tree.morton.reserve(entries.size());
  for (const MortonEntry &entry : entries) {
    tree.order.push_back(entry.index);
    tree.morton.push_back(entry.code);
  }

  const bool tight_bounds = tree.node_bounds_mode == "tight-aabb";
  if (tree.tree_builder == "topdown-aabb") {
    std::vector<TopDownTempNode> temp_nodes;
    temp_nodes.reserve((panels.size() + effective_leaf_size - 1) /
                           effective_leaf_size * 2 +
                       8);
    const std::size_t temp_root = build_topdown_temp_node(
        panels, tree.order, 0, panels.size(), effective_leaf_size, 0,
        temp_nodes);
    flatten_topdown_tree(panels, tree, temp_nodes, temp_root, tight_bounds);
  } else {
    std::vector<std::size_t> level_indices;
    const std::size_t leaf_count =
        (panels.size() + effective_leaf_size - 1) / effective_leaf_size;
    tree.nodes.reserve(leaf_count * 2 + 8);
    level_indices.reserve(leaf_count);
    for (std::size_t leaf = 0; leaf < leaf_count; ++leaf) {
      const std::size_t first_panel = leaf * effective_leaf_size;
      const std::size_t count =
          std::min(effective_leaf_size, panels.size() - first_panel);
      level_indices.push_back(tree.nodes.size());
      tree.nodes.push_back(
          make_leaf_node(panels, tree.order, first_panel, count, tight_bounds));
    }
    tree.leaf_count = tree.nodes.size();
    tree.moment_levels.clear();
    tree.moment_levels.push_back({0, tree.leaf_count});

    while (level_indices.size() > 1) {
      std::vector<std::size_t> next_level;
      next_level.reserve((level_indices.size() + kBranchFactor - 1) /
                         kBranchFactor);
      const std::size_t level_begin = tree.nodes.size();
      for (std::size_t begin = 0; begin < level_indices.size();
           begin += kBranchFactor) {
        const std::size_t child_count =
            std::min(kBranchFactor, level_indices.size() - begin);
        std::vector<std::size_t> children;
        children.reserve(child_count);
        for (std::size_t child = 0; child < child_count; ++child)
          children.push_back(level_indices[begin + child]);
        TreeNode parent = make_parent_node(tree.nodes, children, tight_bounds);
        attach_children(tree, parent, children);
        next_level.push_back(tree.nodes.size());
        tree.nodes.push_back(parent);
      }
      tree.moment_levels.push_back({level_begin, tree.nodes.size()});
      level_indices.swap(next_level);
    }
    tree.root_index = level_indices.empty() ? 0 : level_indices.front();
    tree.depth_max =
        tree.moment_levels.empty() ? 0 : tree.moment_levels.size() - 1;
  }
  tree.morton.clear();
  tree.morton.reserve(tree.order.size());
  for (std::size_t panel_index : tree.order)
    tree.morton.push_back(panel_index < panel_morton.size()
                              ? panel_morton[panel_index]
                              : std::uint64_t{0});

  if (!tree.nodes.empty()) {
    double radius_sum = 0.0;
    for (const TreeNode &node : tree.nodes) {
      tree.node_radius_max = std::max(tree.node_radius_max, node.radius);
      radius_sum += node.radius;
    }
    tree.node_radius_mean = radius_sum / static_cast<double>(tree.nodes.size());
  }
  tree.build_ms = elapsed_ms(build_begin, std::chrono::steady_clock::now());
  tree.panel_cache = build_tree_panel_cache(panels, tree.order);
  if (tree.interaction_mode != "onthefly")
    tree.interactions = build_tree_interactions(panels, tree, executor);
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
    stats->derivative_mode = "none";
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
      compute_node_moments(panels, density, tree, taylor_order, executor);
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
  std::vector<std::size_t> local_mac_reject(workers, 0);
  std::vector<MacAuditTally> local_mac_audit(workers);
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
              const MacDecision mac = evaluate_mac(target, node, tree, distance);
              if (mac.ratio_accept && tree.mac_safety != "off") {
                record_far_geometric_audit(local_mac_audit[worker], node,
                                           mac.gap, mac.gap_ratio,
                                           mac.unsafe);
                if (!mac.accepted)
                  ++local_mac_audit[worker].physical_gap_rejects;
              }
              if (mac.accepted) {
                coupled += node_kernel_taylor(target, node,
                                               moments[node_index], taylor_order,
                                               kappa, tree.derivative_mode);
                ++far_count;
                continue;
              }
              ++local_mac_reject[worker];
            }
            if (leaf) {
              coupled += direct_ordered_sum(
                  target, i, panels, density, kappa, tree.order,
                  node.first_panel, node.panel_count, near_count);
            } else {
              for (std::size_t child = 0; child < node.child_count; ++child) {
                const std::size_t child_index =
                    child_node_index(tree, node, child);
                if (child_index < tree.nodes.size())
                  stack.push_back(child_index);
              }
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
    MacAuditTally mac_audit_total;
    for (const MacAuditTally &worker_audit : local_mac_audit)
      merge_far_geometric_audit(mac_audit_total, worker_audit);
    stats->matvec_mode = "treecode";
    stats->taylor_order = taylor_order;
    stats->derivative_mode = tree.derivative_mode;
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
    stats->mac_reject_radius_limited += std::accumulate(
        local_mac_reject.begin(), local_mac_reject.end(), std::size_t{0});
    stats->tree_builder = tree.tree_builder;
    stats->mac_safety = tree.mac_safety;
    stats->mac_gap_factor = tree.mac_gap_factor;
    stats->moment_audit = tree.moment_audit;
    copy_far_geometric_audit_to_stats(*stats, mac_audit_total);
    stats->tree_depth_max = tree.depth_max;
  }
}

void apply_jcp_treecode_operator(const std::vector<Panel> &panels,
                                 const std::vector<double> &unknowns,
                                 double eps_ratio, double kappa,
                                 double lhs_phi, double lhs_phin,
                                 bool coupled, const TreeWorkspace &tree,
                                 const std::string &near_integration,
                                 const std::string &quad_buffer,
                                 TitanExecutor &executor,
                                 std::vector<double> &out,
                                 TreecodeStats *stats) {
  const std::size_t n = panels.size();
  if (n == 0) {
    out.clear();
    return;
  }
  if (tree.nodes.empty() || tree.order.size() != n ||
      unknowns.size() < (coupled ? 2 * n : n)) {
    out.assign(coupled ? 2 * n : n, 0.0);
    const std::vector<std::size_t> order = [&] {
      std::vector<std::size_t> identity(n);
      std::iota(identity.begin(), identity.end(), std::size_t{0});
      return identity;
    }();
    const std::size_t workers = effective_workers(executor, n);
    std::vector<std::size_t> local_near(workers, 0);
    std::vector<TreecodeStats> local_near_stats(workers);
    const auto near_begin = std::chrono::steady_clock::now();
    executor.parallel_for(
        n, workers, 8,
        [&](std::size_t worker, std::size_t begin, std::size_t end) {
          for (std::size_t i = begin; i < end; ++i) {
            TreecodeStats &worker_stats = local_near_stats[worker];
            const JcpNearAccum coupled_sum = accumulate_jcp_near_ordered(
                panels[i], i, panels, unknowns, eps_ratio, kappa, coupled,
                order, 0, n, near_integration, quad_buffer, &worker_stats);
            out[i] = lhs_phi * unknowns[i] - coupled_sum.phi;
            if (coupled)
              out[n + i] = lhs_phin * unknowns[n + i] - coupled_sum.phin;
          }
          local_near[worker] += local_near_stats[worker].near_centroid_pairs +
                                local_near_stats[worker].near_duffy_pairs +
                                local_near_stats[worker]
                                    .near_geometric_buffer_pairs;
        });
    if (stats != nullptr) {
      stats->matvec_mode = "jcp2013_direct_fallback";
      stats->taylor_order = 0;
      stats->derivative_mode = "none";
      stats->leaf_near_policy = "exact";
      stats->leaf_panel_kernel = "scalar";
      stats->near_integration = near_integration;
      stats->quad_buffer = quad_buffer;
      const double near_ms =
          elapsed_ms(near_begin, std::chrono::steady_clock::now());
      stats->near_integration_ms += near_ms;
      stats->near_eval_ms += near_ms;
      std::size_t near_total = 0;
      for (const TreecodeStats &worker_stats : local_near_stats) {
        const std::size_t previous_duffy = stats->near_duffy_pairs;
        stats->near_kernel_policy = worker_stats.near_kernel_policy;
        stats->near_self_pairs += worker_stats.near_self_pairs;
        stats->near_one_ring_pairs += worker_stats.near_one_ring_pairs;
        stats->near_duffy_pairs += worker_stats.near_duffy_pairs;
        stats->near_centroid_pairs += worker_stats.near_centroid_pairs;
        stats->near_fallbacks += worker_stats.near_fallbacks;
        merge_near_kernel_block_stats(*stats, worker_stats);
        stats->duffy_subdivisions_max =
            std::max(stats->duffy_subdivisions_max,
                     worker_stats.duffy_subdivisions_max);
        if (worker_stats.near_duffy_pairs > 0) {
          const double old_count = static_cast<double>(previous_duffy);
          const double new_count = static_cast<double>(stats->near_duffy_pairs);
          stats->duffy_subdivisions_mean =
              new_count <= 0.0
                  ? 0.0
                  : (stats->duffy_subdivisions_mean * old_count +
                     worker_stats.duffy_subdivisions_mean *
                         static_cast<double>(worker_stats.near_duffy_pairs)) /
                        new_count;
        }
        near_total += worker_stats.near_centroid_pairs +
                      worker_stats.near_duffy_pairs +
                      worker_stats.near_geometric_buffer_pairs;
      }
      stats->near_interactions += near_total;
    }
    return;
  }

  const auto eval_begin = std::chrono::steady_clock::now();
  const int taylor_order = std::clamp(tree.taylor_order, 0, kTaylorOrderMax);
  TreePanelCache fallback_panel_cache;
  const TreePanelCache *panel_cache = &tree.panel_cache;
  if (!panel_cache->valid_for(n)) {
    fallback_panel_cache = build_tree_panel_cache(panels, tree.order);
    panel_cache = &fallback_panel_cache;
  }
  const TreeUnknownCache unknown_cache =
      build_tree_unknown_cache(panels, unknowns, coupled, tree.order);
  const std::string leaf_panel_kernel =
      resolve_leaf_panel_kernel(tree.leaf_panel_kernel);
  const auto moment_begin = std::chrono::steady_clock::now();
  const PackedJcpMoments jcp_moments = compute_jcp_node_moments_packed(
      panels, unknowns, coupled, tree, taylor_order, executor);
  const MomentAuditSummary moment_audit_summary = audit_jcp_moment_translation(
      panels, unknowns, coupled, tree, taylor_order, jcp_moments);
  const double moment_ms =
      elapsed_ms(moment_begin, std::chrono::steady_clock::now());

  if (tree.interaction_mode == "onthefly") {
    const auto interaction_begin = std::chrono::steady_clock::now();
    std::vector<std::size_t> rank(n, 0);
    for (std::size_t pos = 0; pos < tree.order.size(); ++pos)
      rank[tree.order[pos]] = pos;
    const double interaction_ms =
        elapsed_ms(interaction_begin, std::chrono::steady_clock::now());

    out.assign(coupled ? 2 * n : n, 0.0);
    const std::size_t workers = effective_workers(executor, n);
    std::vector<std::size_t> local_leaf_residual(workers, 0);
    std::vector<std::size_t> local_far(workers, 0);
    std::vector<std::size_t> local_mac_reject(workers, 0);
    std::vector<MacAuditTally> local_mac_audit(workers);
    std::vector<TreecodeStats> local_near_stats(workers);
    const auto traversal_begin = std::chrono::steady_clock::now();
    executor.parallel_for(
        n, workers, 8,
        [&](std::size_t worker, std::size_t begin, std::size_t end) {
          TreecodeStats &worker_near_stats = local_near_stats[worker];
          std::vector<std::size_t> stack;
          std::vector<std::size_t> singular_order;
          std::vector<std::size_t> leaf_order;
          stack.reserve(64);
          singular_order.reserve(tree.leaf_size);
          leaf_order.reserve(tree.leaf_size);
          std::size_t leaf_count = 0;
          std::size_t far_count = 0;
          for (std::size_t i = begin; i < end; ++i) {
            const Panel &target = panels[i];
            const std::size_t target_rank = rank[i];
            double sum_phi_i = 0.0;
            double sum_phin_i = 0.0;
            stack.clear();
            stack.push_back(tree.root_index);
            while (!stack.empty()) {
              const std::size_t node_index = stack.back();
              stack.pop_back();
              if (node_index >= tree.nodes.size())
                continue;
              const TreeNode &node = tree.nodes[node_index];
              const bool is_leaf = node.child_count == 0;
              const bool contains_target =
                  target_rank >= node.first_panel &&
                  target_rank < node.first_panel + node.panel_count;
              if (!is_leaf && !contains_target) {
                const double dx = target.x - node.cx;
                const double dy = target.y - node.cy;
                const double dz = target.z - node.cz;
                const double distance = std::sqrt(
                    std::max(1.0e-16, dx * dx + dy * dy + dz * dz));
                const MacDecision mac =
                    evaluate_mac(target, node, tree, distance);
                if (mac.ratio_accept && tree.mac_safety != "off") {
                  record_far_geometric_audit(local_mac_audit[worker], node,
                                             mac.gap, mac.gap_ratio,
                                             mac.unsafe);
                  if (!mac.accepted)
                    ++local_mac_audit[worker].physical_gap_rejects;
                }
                if (mac.accepted) {
                  const std::array<NodeMoments, 3> normal_moments = {
                      jcp_moments.phi_normal[0][node_index],
                      jcp_moments.phi_normal[1][node_index],
                      jcp_moments.phi_normal[2][node_index]};
                  const JcpAccum far = taylor_jcp_node_sum(
                      target, node, jcp_moments.phin[node_index],
                      normal_moments, eps_ratio, kappa, coupled,
                      taylor_order, tree.derivative_mode);
                  sum_phi_i += far.phi;
                  if (coupled)
                    sum_phin_i += far.phin;
                  ++far_count;
                  continue;
                }
                ++local_mac_reject[worker];
              }
              if (is_leaf) {
                singular_order.clear();
                leaf_order.clear();
                for (std::size_t pos = node.first_panel;
                     pos < node.first_panel + node.panel_count; ++pos) {
                  const std::size_t source_index = tree.order[pos];
                  const Panel &source = panels[source_index];
                  if (is_singular_near_pair(target, i, source, source_index) ||
                      classify_geometric_near_buffer_pair(
                          target, source, tree.quad_buffer)
                          .active) {
                    singular_order.push_back(source_index);
                  } else {
                    leaf_order.push_back(source_index);
                  }
                }
                if (!singular_order.empty()) {
                  const JcpNearAccum near = accumulate_jcp_near_ordered(
                      target, i, panels, unknowns, eps_ratio, kappa, coupled,
                      singular_order, 0, singular_order.size(),
                      near_integration, quad_buffer, &worker_near_stats);
                  sum_phi_i += near.phi;
                  if (coupled)
                    sum_phin_i += near.phin;
                }
                if (!leaf_order.empty()) {
                  const JcpAccum near = direct_jcp_ordered_sum(
                      target, i, *panel_cache, unknown_cache, eps_ratio, kappa,
                      coupled, leaf_order, 0, leaf_order.size(),
                      leaf_panel_kernel, leaf_count);
                  sum_phi_i += near.phi;
                  if (coupled)
                    sum_phin_i += near.phin;
                }
              } else {
                for (std::size_t child = 0; child < node.child_count; ++child) {
                  const std::size_t child_index =
                      child_node_index(tree, node, child);
                  if (child_index < tree.nodes.size())
                    stack.push_back(child_index);
                }
              }
            }
            out[i] = lhs_phi * unknowns[i] - sum_phi_i;
            if (coupled)
              out[n + i] = lhs_phin * unknowns[n + i] - sum_phin_i;
          }
          local_leaf_residual[worker] += leaf_count;
          local_far[worker] += far_count;
        });
    const double traversal_ms =
        elapsed_ms(traversal_begin, std::chrono::steady_clock::now());

    if (stats != nullptr) {
      std::size_t singular_evaluated = 0;
      for (const TreecodeStats &worker_stats : local_near_stats)
        singular_evaluated += worker_stats.near_centroid_pairs +
                               worker_stats.near_duffy_pairs +
                               worker_stats.near_geometric_buffer_pairs;
      const std::size_t leaf_total =
          std::accumulate(local_leaf_residual.begin(),
                          local_leaf_residual.end(), std::size_t{0});
      const std::size_t far_total =
          std::accumulate(local_far.begin(), local_far.end(), std::size_t{0});
      MacAuditTally mac_audit_total;
      for (const MacAuditTally &worker_audit : local_mac_audit)
        merge_far_geometric_audit(mac_audit_total, worker_audit);
      stats->matvec_mode = "jcp2013_treecode";
      stats->interaction_mode = "onthefly";
      stats->taylor_order = taylor_order;
      stats->derivative_mode = tree.derivative_mode;
      stats->leaf_near_policy = tree.leaf_near_policy;
      stats->leaf_panel_kernel = leaf_panel_kernel;
      stats->leaf_tile_status = "off";
      stats->near_integration = near_integration;
      stats->quad_buffer = quad_buffer;
      stats->tree_moment_ms += moment_ms;
      if (stats->matvec_calls == 0)
        stats->interaction_build_ms += interaction_ms;
      stats->tree_traversal_ms += traversal_ms;
      stats->tree_eval_ms +=
          elapsed_ms(eval_begin, std::chrono::steady_clock::now());
      stats->leaf_panel_eval_ms += traversal_ms;
      stats->leaf_residual_eval_ms += traversal_ms;
      stats->leaf_near_eval_ms += traversal_ms;
      stats->near_eval_ms += traversal_ms;
      for (const TreecodeStats &worker_stats : local_near_stats) {
        const std::size_t previous_duffy = stats->near_duffy_pairs;
        stats->near_kernel_policy = worker_stats.near_kernel_policy;
        stats->near_self_pairs += worker_stats.near_self_pairs;
        stats->near_one_ring_pairs += worker_stats.near_one_ring_pairs;
        stats->near_duffy_pairs += worker_stats.near_duffy_pairs;
        stats->near_centroid_pairs += worker_stats.near_centroid_pairs;
        stats->near_fallbacks += worker_stats.near_fallbacks;
        merge_near_kernel_block_stats(*stats, worker_stats);
        stats->duffy_subdivisions_max =
            std::max(stats->duffy_subdivisions_max,
                     worker_stats.duffy_subdivisions_max);
        if (worker_stats.near_duffy_pairs > 0) {
          const double old_count = static_cast<double>(previous_duffy);
          const double new_count =
              static_cast<double>(stats->near_duffy_pairs);
          stats->duffy_subdivisions_mean =
              new_count <= 0.0
                  ? 0.0
                  : (stats->duffy_subdivisions_mean * old_count +
                     worker_stats.duffy_subdivisions_mean *
                         static_cast<double>(worker_stats.near_duffy_pairs)) /
                        new_count;
        }
      }
      stats->near_centroid_pairs += leaf_total;
      stats->singular_near_interactions += singular_evaluated;
      stats->leaf_near_interactions += leaf_total;
      stats->leaf_panel_interactions += leaf_total;
      stats->leaf_residual_panel_interactions += leaf_total;
      stats->far_node_interactions += far_total;
      stats->mac_reject_radius_limited += std::accumulate(
          local_mac_reject.begin(), local_mac_reject.end(), std::size_t{0});
      stats->near_interactions += singular_evaluated + leaf_total;
      stats->far_interactions += far_total;
      stats->leaf_panel_pairs_per_second =
          stats->leaf_panel_eval_ms > 0.0
              ? static_cast<double>(stats->leaf_panel_interactions) /
                    (stats->leaf_panel_eval_ms / 1000.0)
              : 0.0;
      stats->node_bounds_mode = tree.node_bounds_mode;
      stats->tree_builder = tree.tree_builder;
      stats->mac_safety = tree.mac_safety;
      stats->mac_gap_factor = tree.mac_gap_factor;
      stats->moment_audit = tree.moment_audit;
      stats->moment_audit_nodes += moment_audit_summary.nodes;
      stats->moment_translation_rel_error_max =
          std::max(stats->moment_translation_rel_error_max,
                   moment_audit_summary.max_rel_error);
      if (moment_audit_summary.nodes > 0) {
        const double old_count =
            static_cast<double>(stats->moment_audit_nodes -
                                moment_audit_summary.nodes);
        const double new_count = static_cast<double>(stats->moment_audit_nodes);
        stats->moment_translation_rel_error_mean =
            new_count > 0.0
                ? (stats->moment_translation_rel_error_mean * old_count +
                   moment_audit_summary.mean_rel_error *
                       static_cast<double>(moment_audit_summary.nodes)) /
                      new_count
                : 0.0;
      }
      copy_far_geometric_audit_to_stats(*stats, mac_audit_total);
      stats->tree_depth_max = tree.depth_max;
      stats->node_radius_max = tree.node_radius_max;
      stats->node_radius_mean = tree.node_radius_mean;
      const double calls_after = static_cast<double>(stats->matvec_calls + 1);
      stats->leaf_panel_work_per_matvec =
          static_cast<double>(stats->leaf_panel_interactions) / calls_after;
      stats->far_node_work_per_matvec =
          static_cast<double>(stats->far_node_interactions) / calls_after;
      const double mac_denom =
          static_cast<double>(stats->far_node_interactions +
                              stats->leaf_panel_interactions +
                              stats->leaf_tile_panel_interactions +
                              stats->leaf_residual_panel_interactions);
      stats->mac_accept_rate =
          mac_denom > 0.0
              ? static_cast<double>(stats->far_node_interactions) / mac_denom
              : 0.0;
      ++stats->matvec_calls;
    }
    return;
  }

  const TreeInteractionLists &interactions = tree.interactions;
  const double interaction_ms = interactions.valid_for(n) ? interactions.build_ms
                                                          : 0.0;
  out.assign(coupled ? 2 * n : n, 0.0);
  const std::size_t workers = effective_workers(executor, n);
  std::vector<double> sum_phi(n, 0.0);
  std::vector<double> sum_phin(coupled ? n : 0, 0.0);
  std::vector<std::size_t> local_leaf_residual(workers, 0);
  std::vector<std::size_t> local_leaf_tile_nodes(workers, 0);
  std::vector<std::size_t> local_leaf_tile_panels(workers, 0);
  std::vector<std::size_t> local_smooth_leaf(workers, 0);
  std::vector<std::size_t> local_far(workers, 0);
  std::vector<TreecodeStats> local_near_stats(workers);

  const auto singular_begin = std::chrono::steady_clock::now();
  executor.parallel_for(
      n, workers, 8,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        TreecodeStats &worker_near_stats = local_near_stats[worker];
        for (std::size_t i = begin; i < end; ++i) {
          const Panel &target = panels[i];
          const std::size_t singular_begin_pos =
              interactions.singular_offsets[i];
          const std::size_t singular_count =
              interactions.singular_offsets[i + 1] - singular_begin_pos;
          if (singular_count == 0)
            continue;
          const JcpNearAccum near = accumulate_jcp_near_ordered(
              target, i, panels, unknowns, eps_ratio, kappa, coupled,
              interactions.singular_indices, singular_begin_pos,
              singular_count, near_integration, quad_buffer,
              &worker_near_stats);
          sum_phi[i] += near.phi;
          if (coupled)
            sum_phin[i] += near.phin;
        }
      });
  const double singular_ms =
      elapsed_ms(singular_begin, std::chrono::steady_clock::now());

  const auto leaf_begin = std::chrono::steady_clock::now();
  executor.parallel_for(
      n, workers, 8,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        std::size_t residual_count = 0;
        for (std::size_t i = begin; i < end; ++i) {
          const std::size_t leaf_begin_pos = interactions.leaf_offsets[i];
          const std::size_t leaf_count_i =
              interactions.leaf_offsets[i + 1] - leaf_begin_pos;
          if (leaf_count_i == 0)
            continue;
          const JcpAccum near = direct_jcp_ordered_sum(
              panels[i], i, *panel_cache, unknown_cache, eps_ratio, kappa,
              coupled, interactions.leaf_indices, leaf_begin_pos,
              leaf_count_i, leaf_panel_kernel, residual_count);
          sum_phi[i] += near.phi;
          if (coupled)
            sum_phin[i] += near.phin;
        }
        local_leaf_residual[worker] += residual_count;
      });
  const double leaf_residual_ms =
      elapsed_ms(leaf_begin, std::chrono::steady_clock::now());

  const auto tile_begin = std::chrono::steady_clock::now();
  executor.parallel_for(
      n, workers, 8,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        std::size_t tile_node_count = 0;
        std::size_t tile_panel_count = 0;
        for (std::size_t i = begin; i < end; ++i) {
          const std::size_t tile_begin_pos =
              interactions.leaf_tile_offsets[i];
          const std::size_t tile_count_i =
              interactions.leaf_tile_offsets[i + 1] - tile_begin_pos;
          for (std::size_t pos = tile_begin_pos;
               pos < tile_begin_pos + tile_count_i; ++pos) {
            const std::size_t node_index =
                interactions.leaf_tile_node_indices[pos];
            if (node_index >= tree.nodes.size())
              continue;
            const std::size_t span_begin =
                pos + 1 < interactions.leaf_tile_span_offsets.size()
                    ? interactions.leaf_tile_span_offsets[pos]
                    : 0;
            const std::size_t span_end =
                pos + 1 < interactions.leaf_tile_span_offsets.size()
                    ? interactions.leaf_tile_span_offsets[pos + 1]
                    : span_begin;
            for (std::size_t span = span_begin; span < span_end; ++span) {
              if (span >= interactions.leaf_tile_span_firsts.size() ||
                  span >= interactions.leaf_tile_span_counts.size())
                continue;
              TreeNode span_node = tree.nodes[node_index];
              span_node.first_panel = interactions.leaf_tile_span_firsts[span];
              span_node.panel_count = interactions.leaf_tile_span_counts[span];
              std::size_t counted_panels = 0;
              const JcpAccum near = direct_jcp_ordered_tile_sum(
                  panels[i], i, *panel_cache, unknown_cache, eps_ratio, kappa,
                  coupled, span_node, leaf_panel_kernel, counted_panels);
              sum_phi[i] += near.phi;
              if (coupled)
                sum_phin[i] += near.phin;
              tile_panel_count += counted_panels;
            }
            ++tile_node_count;
          }
        }
        local_leaf_tile_nodes[worker] += tile_node_count;
        local_leaf_tile_panels[worker] += tile_panel_count;
      });
  const double leaf_tile_ms =
      elapsed_ms(tile_begin, std::chrono::steady_clock::now());
  const double leaf_ms = leaf_residual_ms + leaf_tile_ms;

  const auto smooth_leaf_begin = std::chrono::steady_clock::now();
  executor.parallel_for(
      n, workers, 8,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        std::size_t smooth_count = 0;
        for (std::size_t i = begin; i < end; ++i) {
          const Panel &target = panels[i];
          const std::size_t smooth_begin_pos =
              interactions.smooth_leaf_offsets[i];
          const std::size_t smooth_count_i =
              interactions.smooth_leaf_offsets[i + 1] - smooth_begin_pos;
          for (std::size_t pos = smooth_begin_pos;
               pos < smooth_begin_pos + smooth_count_i; ++pos) {
            const std::size_t node_index =
                interactions.smooth_leaf_node_indices[pos];
            const TreeNode &node = tree.nodes[node_index];
            const std::array<NodeMoments, 3> normal_moments = {
                jcp_moments.phi_normal[0][node_index],
                jcp_moments.phi_normal[1][node_index],
                jcp_moments.phi_normal[2][node_index]};
            const JcpAccum smooth = taylor_jcp_node_sum(
                target, node, jcp_moments.phin[node_index], normal_moments,
                eps_ratio, kappa, coupled, taylor_order,
                tree.derivative_mode);
            sum_phi[i] += smooth.phi;
            if (coupled)
              sum_phin[i] += smooth.phin;
            ++smooth_count;
          }
        }
        local_smooth_leaf[worker] += smooth_count;
      });
  const double smooth_leaf_ms =
      elapsed_ms(smooth_leaf_begin, std::chrono::steady_clock::now());

  const auto far_begin = std::chrono::steady_clock::now();
  executor.parallel_for(
      n, workers, 8,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        std::size_t far_count = 0;
        for (std::size_t i = begin; i < end; ++i) {
          const Panel &target = panels[i];
          const std::size_t far_begin_pos = interactions.far_offsets[i];
          const std::size_t far_count_i =
              interactions.far_offsets[i + 1] - far_begin_pos;
          for (std::size_t pos = far_begin_pos;
               pos < far_begin_pos + far_count_i; ++pos) {
            const std::size_t node_index = interactions.far_node_indices[pos];
            const TreeNode &node = tree.nodes[node_index];
            const std::array<NodeMoments, 3> normal_moments = {
                jcp_moments.phi_normal[0][node_index],
                jcp_moments.phi_normal[1][node_index],
                jcp_moments.phi_normal[2][node_index]};
            const JcpAccum far = taylor_jcp_node_sum(
                target, node, jcp_moments.phin[node_index], normal_moments,
                eps_ratio, kappa, coupled, taylor_order,
                tree.derivative_mode);
            sum_phi[i] += far.phi;
            if (coupled)
              sum_phin[i] += far.phin;
            ++far_count;
          }
          out[i] = lhs_phi * unknowns[i] - sum_phi[i];
          if (coupled)
            out[n + i] = lhs_phin * unknowns[n + i] - sum_phin[i];
        }
        local_far[worker] += far_count;
        });
  const double far_ms = elapsed_ms(far_begin, std::chrono::steady_clock::now());
  const double traversal_ms = singular_ms + leaf_ms + smooth_leaf_ms + far_ms;

  if (stats != nullptr) {
    std::size_t singular_evaluated = 0;
    for (const TreecodeStats &worker_stats : local_near_stats)
      singular_evaluated += worker_stats.near_centroid_pairs +
                             worker_stats.near_duffy_pairs +
                             worker_stats.near_geometric_buffer_pairs;
    const std::size_t leaf_residual_total =
        std::accumulate(local_leaf_residual.begin(), local_leaf_residual.end(),
                        std::size_t{0});
    const std::size_t leaf_tile_node_total =
        std::accumulate(local_leaf_tile_nodes.begin(),
                        local_leaf_tile_nodes.end(), std::size_t{0});
    const std::size_t leaf_tile_panel_total =
        std::accumulate(local_leaf_tile_panels.begin(),
                        local_leaf_tile_panels.end(), std::size_t{0});
    const std::size_t leaf_total =
        leaf_residual_total + leaf_tile_panel_total;
    const std::size_t smooth_leaf_total =
        std::accumulate(local_smooth_leaf.begin(), local_smooth_leaf.end(),
                        std::size_t{0});
    const std::size_t near_total =
        singular_evaluated + leaf_total + smooth_leaf_total;
    const std::size_t far_total =
        std::accumulate(local_far.begin(), local_far.end(), std::size_t{0});
    stats->matvec_mode = "jcp2013_treecode";
    stats->taylor_order = taylor_order;
    stats->derivative_mode = tree.derivative_mode;
    stats->leaf_near_policy = tree.leaf_near_policy;
    stats->leaf_panel_kernel = leaf_panel_kernel;
    stats->interaction_mode = tree.interaction_mode;
    stats->near_integration = near_integration;
    stats->quad_buffer = quad_buffer;
    stats->tree_moment_ms += moment_ms;
    if (stats->matvec_calls == 0)
      stats->interaction_build_ms += interaction_ms;
    stats->tree_traversal_ms += traversal_ms;
    stats->tree_eval_ms +=
        elapsed_ms(eval_begin, std::chrono::steady_clock::now());
    stats->singular_near_eval_ms += singular_ms;
    stats->leaf_panel_eval_ms += leaf_ms;
    stats->leaf_tile_eval_ms += leaf_tile_ms;
    stats->leaf_residual_eval_ms += leaf_residual_ms;
    stats->smooth_leaf_node_eval_ms += smooth_leaf_ms;
    stats->leaf_near_eval_ms += leaf_ms + smooth_leaf_ms;
    stats->near_eval_ms += singular_ms + leaf_ms + smooth_leaf_ms;
    stats->far_eval_ms += far_ms;
    stats->near_integration_ms += singular_ms;
    for (const TreecodeStats &worker_stats : local_near_stats) {
      const std::size_t previous_duffy = stats->near_duffy_pairs;
      stats->near_kernel_policy = worker_stats.near_kernel_policy;
      stats->near_self_pairs += worker_stats.near_self_pairs;
      stats->near_one_ring_pairs += worker_stats.near_one_ring_pairs;
      stats->near_duffy_pairs += worker_stats.near_duffy_pairs;
      stats->near_centroid_pairs += worker_stats.near_centroid_pairs;
      stats->near_fallbacks += worker_stats.near_fallbacks;
      merge_near_kernel_block_stats(*stats, worker_stats);
      stats->duffy_subdivisions_max =
          std::max(stats->duffy_subdivisions_max,
                   worker_stats.duffy_subdivisions_max);
      if (worker_stats.near_duffy_pairs > 0) {
        const double old_count = static_cast<double>(previous_duffy);
        const double new_count = static_cast<double>(stats->near_duffy_pairs);
        stats->duffy_subdivisions_mean =
            new_count <= 0.0
                ? 0.0
                : (stats->duffy_subdivisions_mean * old_count +
                   worker_stats.duffy_subdivisions_mean *
                       static_cast<double>(worker_stats.near_duffy_pairs)) /
                      new_count;
      }
    }
    stats->near_centroid_pairs += leaf_total;
    stats->singular_near_interactions +=
        interactions.valid_for(n) ? interactions.singular_pair_count : 0;
    stats->leaf_near_interactions += leaf_total;
    stats->leaf_panel_interactions += leaf_total;
    stats->leaf_tile_interactions += leaf_tile_node_total;
    stats->leaf_tile_panel_interactions += leaf_tile_panel_total;
    stats->leaf_residual_panel_interactions += leaf_residual_total;
    stats->leaf_tile_status =
        tile_fallback_status(leaf_panel_kernel, leaf_tile_node_total);
    stats->smooth_leaf_node_interactions += smooth_leaf_total;
    stats->leaf_pruned_interactions +=
        interactions.valid_for(n) ? interactions.leaf_pruned_count : 0;
    stats->adaptive_smooth_leaf_accepted +=
        interactions.valid_for(n) ? interactions.adaptive_smooth_leaf_accepted
                                  : 0;
    stats->adaptive_smooth_leaf_rejected +=
        interactions.valid_for(n) ? interactions.adaptive_smooth_leaf_rejected
                                  : 0;
    stats->adaptive_span_tile_accepted +=
        interactions.valid_for(n) ? interactions.adaptive_span_tile_accepted
                                  : 0;
    stats->adaptive_span_tile_rejected +=
        interactions.valid_for(n) ? interactions.adaptive_span_tile_rejected
                                  : 0;
    stats->adaptive_exact_leaf_fallbacks +=
        interactions.valid_for(n) ? interactions.adaptive_exact_leaf_fallbacks
                                  : 0;
    stats->adaptive_reject_one_ring_overlap +=
        interactions.valid_for(n)
            ? interactions.adaptive_reject_one_ring_overlap
            : 0;
    stats->adaptive_reject_error_bound_failed +=
        interactions.valid_for(n)
            ? interactions.adaptive_reject_error_bound_failed
            : 0;
    stats->adaptive_reject_p4_check_failed +=
        interactions.valid_for(n) ? interactions.adaptive_reject_p4_check_failed
                                  : 0;
    stats->adaptive_reject_too_close +=
        interactions.valid_for(n) ? interactions.adaptive_reject_too_close : 0;
    stats->target_cluster_interactions +=
        interactions.valid_for(n) ? interactions.target_cluster_interactions
                                  : 0;
    stats->target_cluster_targets +=
        interactions.valid_for(n) ? interactions.target_cluster_targets : 0;
    stats->target_cluster_panels_pruned +=
        interactions.valid_for(n) ? interactions.target_cluster_panels_pruned
                                  : 0;
    stats->target_cluster_rejected +=
        interactions.valid_for(n) ? interactions.target_cluster_rejected : 0;
    stats->target_cluster_eval_ms += smooth_leaf_ms;
    stats->adaptive_level = tree.adaptive_level;
    const double possible_leaf_work =
        static_cast<double>(stats->leaf_panel_interactions +
                            stats->leaf_pruned_interactions);
    stats->leaf_work_reduction_ratio =
        possible_leaf_work > 0.0
            ? static_cast<double>(stats->leaf_pruned_interactions) /
                  possible_leaf_work
            : 0.0;
    stats->leaf_near_policy = tree.leaf_near_policy;
    stats->far_node_interactions += far_total;
    stats->mac_reject_radius_limited +=
        interactions.valid_for(n) ? interactions.mac_reject_radius_limited : 0;
    stats->near_interactions += near_total;
    stats->far_interactions += far_total;
    stats->leaf_panel_pairs_per_second =
        stats->leaf_panel_eval_ms > 0.0
            ? static_cast<double>(stats->leaf_panel_interactions) /
                  (stats->leaf_panel_eval_ms / 1000.0)
            : 0.0;
    stats->node_bounds_mode = tree.node_bounds_mode;
    stats->tree_builder = tree.tree_builder;
    stats->mac_safety = tree.mac_safety;
    stats->mac_gap_factor = tree.mac_gap_factor;
    stats->moment_audit = tree.moment_audit;
    stats->moment_audit_nodes += moment_audit_summary.nodes;
    stats->moment_translation_rel_error_max =
        std::max(stats->moment_translation_rel_error_max,
                 moment_audit_summary.max_rel_error);
    if (moment_audit_summary.nodes > 0) {
      const double old_count =
          static_cast<double>(stats->moment_audit_nodes -
                              moment_audit_summary.nodes);
      const double new_count = static_cast<double>(stats->moment_audit_nodes);
      stats->moment_translation_rel_error_mean =
          new_count > 0.0
              ? (stats->moment_translation_rel_error_mean * old_count +
                 moment_audit_summary.mean_rel_error *
                     static_cast<double>(moment_audit_summary.nodes)) /
                    new_count
              : 0.0;
    }
    if (interactions.valid_for(n)) {
      MacAuditTally cached_audit;
      cached_audit.unsafe_count = interactions.far_geometric_unsafe_count;
      cached_audit.unsafe_panels = interactions.far_geometric_unsafe_panels;
      cached_audit.min_gap = interactions.far_geometric_min_gap;
      cached_audit.min_gap_ratio = interactions.far_geometric_gap_ratio_min;
      cached_audit.sum_gap_ratio = interactions.far_geometric_gap_ratio_sum;
      cached_audit.samples = interactions.far_geometric_audit_samples;
      cached_audit.physical_gap_rejects =
          interactions.mac_physical_gap_rejects;
      copy_far_geometric_audit_to_stats(*stats, cached_audit);
    }
    stats->tree_depth_max = tree.depth_max;
    stats->node_radius_max = tree.node_radius_max;
    stats->node_radius_mean = tree.node_radius_mean;
    const double calls_after = static_cast<double>(stats->matvec_calls + 1);
    stats->leaf_panel_work_per_matvec =
        static_cast<double>(stats->leaf_panel_interactions) / calls_after;
    stats->far_node_work_per_matvec =
        static_cast<double>(stats->far_node_interactions) / calls_after;
    const double mac_denom =
        static_cast<double>(stats->far_node_interactions +
                            stats->leaf_panel_interactions +
                            stats->leaf_tile_panel_interactions +
                            stats->leaf_residual_panel_interactions);
    stats->mac_accept_rate =
        mac_denom > 0.0
            ? static_cast<double>(stats->far_node_interactions) / mac_denom
            : 0.0;
    ++stats->matvec_calls;
  }
}

} // namespace TitanPBE::Tabi
