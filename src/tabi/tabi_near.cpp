#include "tabi_near.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace TitanPBE::Tabi {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct KernelValues {
  double k1{0.0};
  double k2{0.0};
  double k3{0.0};
  double k4{0.0};
};

enum class RejectReason {
  None,
  CandidateWorseThanCentroid,
  ReferenceNonconverged,
  OrientationSignMismatch,
  StrictReferenceReject,
};

struct KernelBlockDecision {
  bool accepted{false};
  const char *accepted_policy{"centroid"};
  const char *rejected_policy{"centroid_guarded"};
  RejectReason reject_reason{RejectReason::None};
};

struct AcceptedKernelValues {
  KernelValues kernel{};
  std::array<KernelBlockDecision, 4> block{};
  std::array<double, 4> candidate_ref_rel{};
  std::array<double, 4> centroid_ref_rel{};
  std::array<bool, 4> candidate_worse_than_centroid{};
  std::array<bool, 4> k2_laplace_sign_ok{};
  std::array<double, 4> k2_yukawa_smooth_ref_rel{};
  std::array<double, 4> k2_orientation_pair_ref_rel{};
};

struct QuadPoint {
  double l0{0.0};
  double l1{0.0};
  double l2{0.0};
  double weight{0.0};
};

struct Quad1D {
  double x{0.0};
  double weight{0.0};
};

struct Vec3 {
  double x{0.0};
  double y{0.0};
  double z{0.0};
};

struct RawNearIntegrals {
  double g0_minus_gk{0.0};
  double d0_src{0.0};
  double dk_minus_d0_src{0.0};
  double d0_tgt{0.0};
  double dk_tgt{0.0};
  double mixedk_minus_mixed0{0.0};
};

struct K2CandidateInfo {
  bool laplace_sign_ok{true};
  double laplace_integral{0.0};
  double yukawa_smooth_integral{0.0};
};

double green(double r, double kappa) noexcept {
  if (r <= 1.0e-14)
    return 0.0;
  const double kk = std::max(0.0, kappa);
  return std::exp(-kk * r) / (4.0 * kPi * r);
}

double source_normal_green(double tx, double ty, double tz, double sx,
                           double sy, double sz, double snx, double sny,
                           double snz, double kappa) noexcept {
  const double ux = tx - sx;
  const double uy = ty - sy;
  const double uz = tz - sz;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  const double ndotu = snx * ux + sny * uy + snz * uz;
  const double kk = std::max(0.0, kappa);
  return std::exp(-kk * r) * (1.0 + kk * r) * ndotu /
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
  const double ndotu = target.nx * ux + target.ny * uy + target.nz * uz;
  const double kk = std::max(0.0, kappa);
  return -std::exp(-kk * r) * (1.0 + kk * r) * ndotu /
         (4.0 * kPi * r2 * r);
}

double mixed_normal_green(const Panel &target, double sx, double sy, double sz,
                          double snx, double sny, double snz,
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
  const double nyu = (snx * ux + sny * uy + snz * uz) * inv_r;
  const double nxny = target.nx * snx + target.ny * sny + target.nz * snz;
  const double kk = std::max(0.0, kappa);
  const double expkr = std::exp(-kk * r);
  const double inv_r2 = inv_r * inv_r;
  const double fprime =
      -expkr * (1.0 + kk * r) * inv_r2 / (4.0 * kPi);
  const double fsecond =
      expkr *
      (kk * kk * inv_r + 2.0 * kk * inv_r2 + 2.0 * inv_r2 * inv_r) /
      (4.0 * kPi);
  return fsecond * nxu * nyu + (fprime * inv_r) * (nxny - nxu * nyu);
}

KernelValues point_kernel(const Panel &target, const Panel &source, double sx,
                          double sy, double sz, double eps_ratio,
                          double kappa) noexcept {
  const double dx = target.x - sx;
  const double dy = target.y - sy;
  const double dz = target.z - sz;
  const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (r <= 1.0e-14)
    return {};
  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  const double g0 = green(r, 0.0);
  const double gk = green(r, kappa);
  const double d0_src =
      source_normal_green(target.x, target.y, target.z, sx, sy, sz, source.nx,
                          source.ny, source.nz, 0.0);
  const double dk_src =
      source_normal_green(target.x, target.y, target.z, sx, sy, sz, source.nx,
                          source.ny, source.nz, kappa);
  const double d0_tgt = target_normal_green(target, sx, sy, sz, 0.0);
  const double dk_tgt = target_normal_green(target, sx, sy, sz, kappa);
  const double d2_0 =
      mixed_normal_green(target, sx, sy, sz, source.nx, source.ny, source.nz,
                         0.0);
  const double d2_k =
      mixed_normal_green(target, sx, sy, sz, source.nx, source.ny, source.nz,
                         kappa);
  return {g0 - gk, eps_ratio * dk_src - d0_src,
          d0_tgt - inv_eps_ratio * dk_tgt, d2_k - d2_0};
}

Vec3 make_vec(double x, double y, double z) noexcept { return {x, y, z}; }

Vec3 operator+(const Vec3 &a, const Vec3 &b) noexcept {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator-(const Vec3 &a, const Vec3 &b) noexcept {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 operator*(double s, const Vec3 &a) noexcept {
  return {s * a.x, s * a.y, s * a.z};
}

double dot(const Vec3 &a, const Vec3 &b) noexcept {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3 &a, const Vec3 &b) noexcept {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

double norm(const Vec3 &a) noexcept { return std::sqrt(dot(a, a)); }

double distance_sq(const Vec3 &a, const Vec3 &b) noexcept {
  const Vec3 d = a - b;
  return dot(d, d);
}

Vec3 panel_vertex(const Panel &panel, int index) noexcept {
  switch (index) {
  case 0:
    return make_vec(panel.v0x, panel.v0y, panel.v0z);
  case 1:
    return make_vec(panel.v1x, panel.v1y, panel.v1z);
  default:
    return make_vec(panel.v2x, panel.v2y, panel.v2z);
  }
}

Vec3 closest_point_on_triangle(const Vec3 &p, const Vec3 &a, const Vec3 &b,
                               const Vec3 &c) noexcept {
  const Vec3 ab = b - a;
  const Vec3 ac = c - a;
  const Vec3 ap = p - a;
  const double d1 = dot(ab, ap);
  const double d2 = dot(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0)
    return a;

  const Vec3 bp = p - b;
  const double d3 = dot(ab, bp);
  const double d4 = dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3)
    return b;

  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double v = d1 / (d1 - d3);
    return a + v * ab;
  }

  const Vec3 cp = p - c;
  const double d5 = dot(ab, cp);
  const double d6 = dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6)
    return c;

  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double w = d2 / (d2 - d6);
    return a + w * ac;
  }

  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return b + w * (c - b);
  }

  const double denom = 1.0 / (va + vb + vc);
  const double v = vb * denom;
  const double w = vc * denom;
  return a + v * ab + w * ac;
}

double point_triangle_distance_sq(const Vec3 &p, const Vec3 &a, const Vec3 &b,
                                  const Vec3 &c) noexcept {
  return distance_sq(p, closest_point_on_triangle(p, a, b, c));
}

double segment_segment_distance_sq(const Vec3 &p1, const Vec3 &q1,
                                   const Vec3 &p2, const Vec3 &q2) noexcept {
  constexpr double eps = 1.0e-14;
  const Vec3 d1 = q1 - p1;
  const Vec3 d2 = q2 - p2;
  const Vec3 r = p1 - p2;
  const double a = dot(d1, d1);
  const double e = dot(d2, d2);
  const double f = dot(d2, r);
  double s = 0.0;
  double t = 0.0;
  if (a <= eps && e <= eps)
    return distance_sq(p1, p2);
  if (a <= eps) {
    t = std::clamp(f / e, 0.0, 1.0);
  } else {
    const double c = dot(d1, r);
    if (e <= eps) {
      s = std::clamp(-c / a, 0.0, 1.0);
    } else {
      const double b = dot(d1, d2);
      const double denom = a * e - b * b;
      if (std::abs(denom) > eps)
        s = std::clamp((b * f - c * e) / denom, 0.0, 1.0);
      const double tnom = b * s + f;
      if (tnom < 0.0) {
        t = 0.0;
        s = std::clamp(-c / a, 0.0, 1.0);
      } else if (tnom > e) {
        t = 1.0;
        s = std::clamp((b - c) / a, 0.0, 1.0);
      } else {
        t = tnom / e;
      }
    }
  }
  const Vec3 c1 = p1 + s * d1;
  const Vec3 c2 = p2 + t * d2;
  return distance_sq(c1, c2);
}

double panel_max_edge_length(const Panel &panel) noexcept {
  if (!panel.has_triangle_geometry)
    return 0.0;
  const Vec3 a = panel_vertex(panel, 0);
  const Vec3 b = panel_vertex(panel, 1);
  const Vec3 c = panel_vertex(panel, 2);
  return std::sqrt(std::max({distance_sq(a, b), distance_sq(b, c),
                             distance_sq(c, a)}));
}

double panel_centroid_radius(const Panel &panel) noexcept {
  if (!panel.has_triangle_geometry)
    return 0.0;
  const Vec3 center = make_vec(panel.x, panel.y, panel.z);
  return std::sqrt(std::max({distance_sq(center, panel_vertex(panel, 0)),
                             distance_sq(center, panel_vertex(panel, 1)),
                             distance_sq(center, panel_vertex(panel, 2))}));
}

double triangle_triangle_distance_sq(const Panel &a_panel,
                                     const Panel &b_panel) noexcept {
  const std::array<Vec3, 3> a = {panel_vertex(a_panel, 0),
                                panel_vertex(a_panel, 1),
                                panel_vertex(a_panel, 2)};
  const std::array<Vec3, 3> b = {panel_vertex(b_panel, 0),
                                panel_vertex(b_panel, 1),
                                panel_vertex(b_panel, 2)};
  double best = std::numeric_limits<double>::infinity();
  for (const Vec3 &v : a)
    best = std::min(best, point_triangle_distance_sq(v, b[0], b[1], b[2]));
  for (const Vec3 &v : b)
    best = std::min(best, point_triangle_distance_sq(v, a[0], a[1], a[2]));
  for (int i = 0; i < 3; ++i) {
    const Vec3 a0 = a[static_cast<std::size_t>(i)];
    const Vec3 a1 = a[static_cast<std::size_t>((i + 1) % 3)];
    for (int j = 0; j < 3; ++j) {
      const Vec3 b0 = b[static_cast<std::size_t>(j)];
      const Vec3 b1 = b[static_cast<std::size_t>((j + 1) % 3)];
      best = std::min(best, segment_segment_distance_sq(a0, a1, b0, b1));
    }
  }
  return std::isfinite(best) ? best : 0.0;
}

double oriented_solid_angle(const Panel &target,
                            const Panel &source) noexcept {
  const Vec3 x = make_vec(target.x, target.y, target.z);
  const Vec3 a = panel_vertex(source, 0) - x;
  const Vec3 b = panel_vertex(source, 1) - x;
  const Vec3 c = panel_vertex(source, 2) - x;
  const double la = norm(a);
  const double lb = norm(b);
  const double lc = norm(c);
  if (la <= 1.0e-14 || lb <= 1.0e-14 || lc <= 1.0e-14)
    return 0.0;
  const double numerator = dot(a, cross(b, c));
  const double denominator =
      la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
  double omega = 2.0 * std::atan2(numerator, denominator);
  const Vec3 v0 = panel_vertex(source, 0);
  const Vec3 v1 = panel_vertex(source, 1);
  const Vec3 v2 = panel_vertex(source, 2);
  const Vec3 tri_normal = cross(v1 - v0, v2 - v0);
  const Vec3 source_normal = make_vec(source.nx, source.ny, source.nz);
  if (dot(tri_normal, source_normal) < 0.0)
    omega = -omega;
  return omega;
}

bool finite_kernel(const KernelValues &kernel) noexcept {
  return std::isfinite(kernel.k1) && std::isfinite(kernel.k2) &&
         std::isfinite(kernel.k3) && std::isfinite(kernel.k4);
}

bool relative_guard_ok(double candidate, double reference,
                       double rel_limit) noexcept {
  if (!std::isfinite(candidate) || !std::isfinite(reference))
    return false;
  const double scale = std::max({1.0e-14, std::abs(reference),
                                 0.25 * std::abs(candidate)});
  return std::abs(candidate - reference) <= rel_limit * scale;
}

bool same_sign_for_significant_values(double a, double b) noexcept {
  const double scale = std::max({1.0e-14, std::abs(a), std::abs(b)});
  if (std::abs(a) <= 1.0e-12 * scale || std::abs(b) <= 1.0e-12 * scale)
    return true;
  return a * b >= 0.0;
}

double relative_error_value(double value, double reference) noexcept {
  if (!std::isfinite(value) || !std::isfinite(reference))
    return std::numeric_limits<double>::infinity();
  const double scale =
      std::max({1.0e-14, std::abs(reference), 0.25 * std::abs(value)});
  return std::abs(value - reference) / scale;
}

std::size_t near_block_bucket_index(std::size_t block, std::size_t pair_type,
                                    std::size_t orientation) noexcept {
  return (block * kNearPairTypeCount + pair_type) * kNearOrientationCount +
         orientation;
}

std::size_t classify_pair_type(bool self_pair, int shared_vertices) noexcept {
  if (self_pair)
    return 0;
  if (shared_vertices >= 2)
    return 1;
  if (shared_vertices == 1)
    return 2;
  return 3;
}

std::size_t classify_orientation(const Panel &target,
                                 const Panel &source) noexcept {
  const double dotn = target.nx * source.nx + target.ny * source.ny +
                      target.nz * source.nz;
  if (dotn >= 0.25)
    return 0;
  if (dotn <= -0.25)
    return 1;
  return 2;
}

bool finite_kernel_values(const KernelValues &kernel) noexcept {
  return std::isfinite(kernel.k1) && std::isfinite(kernel.k2) &&
         std::isfinite(kernel.k3) && std::isfinite(kernel.k4);
}

KernelValues centroid_integral(const Panel &target, const Panel &source,
                               double eps_ratio, double kappa) noexcept {
  const double sx = source.has_source_quadrature_center ? source.source_x
                                                        : source.x;
  const double sy = source.has_source_quadrature_center ? source.source_y
                                                        : source.y;
  const double sz = source.has_source_quadrature_center ? source.source_z
                                                        : source.z;
  KernelValues kernel = point_kernel(target, source, sx, sy, sz, eps_ratio,
                                     kappa);
  kernel.k1 *= source.area;
  kernel.k2 *= source.area;
  kernel.k3 *= source.area;
  kernel.k4 *= source.area;
  return kernel;
}

[[maybe_unused]] bool share_vertex(const Panel &a, const Panel &b) noexcept {
  if (!a.has_triangle_geometry || !b.has_triangle_geometry)
    return false;
  const std::array<std::size_t, 3> av = {a.v0_index, a.v1_index, a.v2_index};
  const std::array<std::size_t, 3> bv = {b.v0_index, b.v1_index, b.v2_index};
  for (std::size_t ai : av) {
    for (std::size_t bi : bv) {
      if (ai != static_cast<std::size_t>(-1) && ai == bi)
        return true;
    }
  }
  return false;
}

int shared_vertex_count(const Panel &a, const Panel &b) noexcept {
  if (!a.has_triangle_geometry || !b.has_triangle_geometry)
    return 0;
  const std::array<std::size_t, 3> av = {a.v0_index, a.v1_index, a.v2_index};
  const std::array<std::size_t, 3> bv = {b.v0_index, b.v1_index, b.v2_index};
  int count = 0;
  for (std::size_t ai : av) {
    if (ai == static_cast<std::size_t>(-1))
      continue;
    for (std::size_t bi : bv) {
      if (ai == bi) {
        ++count;
        break;
      }
    }
  }
  return count;
}

double self_green_disk_integral(double area, double kappa) noexcept {
  const double radius = std::sqrt(std::max(0.0, area) / kPi);
  if (radius <= 0.0)
    return 0.0;
  const double kk = std::max(0.0, kappa);
  if (kk <= 1.0e-14)
    return 0.5 * radius;
  return (1.0 - std::exp(-kk * radius)) / (2.0 * kk);
}

[[maybe_unused]] KernelValues self_singular_integral(const Panel &source,
                                                     double kappa) noexcept {
  // The jump terms remain on the analytic JCP LHS.  For the collocation
  // source triangle itself, integrate only the single-layer singular part
  // with an equal-area disk surrogate and use flat-panel principal values for
  // the derivative kernels.
  const double g0 = self_green_disk_integral(source.area, 0.0);
  const double gk = self_green_disk_integral(source.area, kappa);
  return {g0 - gk, 0.0, 0.0, 0.0};
}

AcceptedKernelValues self_singular_decision(const Panel &source,
                                            double kappa) noexcept {
  AcceptedKernelValues out;
  out.kernel = self_singular_integral(source, kappa);
  out.block = {KernelBlockDecision{true, "self_asi", "centroid_guarded"},
               KernelBlockDecision{true, "self_flat_pv", "centroid_guarded"},
               KernelBlockDecision{true, "self_flat_pv", "centroid_guarded"},
               KernelBlockDecision{true, "self_flat_pv", "centroid_guarded"}};
  return out;
}

void accumulate_raw_near(RawNearIntegrals &raw, const Panel &target,
                         const Panel &source, const Vec3 &point,
                         double weight, double kappa) noexcept {
  const double dx = target.x - point.x;
  const double dy = target.y - point.y;
  const double dz = target.z - point.z;
  const double r2 = dx * dx + dy * dy + dz * dz;
  if (r2 <= 1.0e-28)
    return;
  const double r = std::sqrt(r2);
  const double g0 = green(r, 0.0);
  const double gk = green(r, kappa);
  const double d0_src =
      source_normal_green(target.x, target.y, target.z, point.x, point.y,
                          point.z, source.nx, source.ny, source.nz, 0.0);
  const double dk_src =
      source_normal_green(target.x, target.y, target.z, point.x, point.y,
                          point.z, source.nx, source.ny, source.nz, kappa);
  const double d0_tgt = target_normal_green(target, point.x, point.y, point.z,
                                            0.0);
  const double dk_tgt =
      target_normal_green(target, point.x, point.y, point.z, kappa);
  const double mixed0 =
      mixed_normal_green(target, point.x, point.y, point.z, source.nx,
                         source.ny, source.nz, 0.0);
  const double mixedk =
      mixed_normal_green(target, point.x, point.y, point.z, source.nx,
                         source.ny, source.nz, kappa);
  raw.g0_minus_gk += weight * (g0 - gk);
  raw.d0_src += weight * d0_src;
  raw.dk_minus_d0_src += weight * (dk_src - d0_src);
  raw.d0_tgt += weight * d0_tgt;
  raw.dk_tgt += weight * dk_tgt;
  raw.mixedk_minus_mixed0 += weight * (mixedk - mixed0);
}

KernelValues duffy_cpv_integral(const Panel &target, const Panel &source,
                                double eps_ratio, double kappa,
                                bool self_pair,
                                K2CandidateInfo *k2_info = nullptr) noexcept {
  static constexpr std::array<Quad1D, 8> q1 = {
      Quad1D{0.0198550717512319, 0.0506142681451881},
      Quad1D{0.1016667612931876, 0.1111905172266872},
      Quad1D{0.2372337950418355, 0.1568533229389436},
      Quad1D{0.4082826787521751, 0.1813418916891810},
      Quad1D{0.5917173212478249, 0.1813418916891810},
      Quad1D{0.7627662049581645, 0.1568533229389436},
      Quad1D{0.8983332387068124, 0.1111905172266872},
      Quad1D{0.9801449282487681, 0.0506142681451881}};

  const Vec3 v0 = panel_vertex(source, 0);
  const Vec3 v1 = panel_vertex(source, 1);
  const Vec3 v2 = panel_vertex(source, 2);
  const Vec3 target_point = make_vec(target.x, target.y, target.z);
  const Vec3 singular =
      self_pair ? target_point : closest_point_on_triangle(target_point, v0, v1, v2);
  const std::array<Vec3, 3> verts = {v0, v1, v2};
  RawNearIntegrals raw;
  for (int edge = 0; edge < 3; ++edge) {
    const Vec3 a = verts[edge];
    const Vec3 b = verts[(edge + 1) % 3];
    const Vec3 pa = a - singular;
    const Vec3 pb = b - singular;
    const double jac_base = norm(cross(pa, pb));
    if (jac_base <= 1.0e-18)
      continue;
    for (const Quad1D &u_quad : q1) {
      const double u = u_quad.x;
      for (const Quad1D &v_quad : q1) {
        const double v = v_quad.x;
        const Vec3 radial = (1.0 - v) * pa + v * pb;
        const Vec3 point = singular + u * radial;
        const double weight = jac_base * u * u_quad.weight * v_quad.weight;
        accumulate_raw_near(raw, target, source, point, weight, kappa);
      }
    }
  }

  double d0_src = raw.d0_src;
  bool laplace_sign_ok = true;
  if (self_pair) {
    d0_src = 0.0;
  } else {
    const double centroid_d0 =
        source_normal_green(target.x, target.y, target.z, source.x, source.y,
                            source.z, source.nx, source.ny, source.nz, 0.0) *
        source.area;
    const double analytic_d0 = -oriented_solid_angle(target, source) /
                               (4.0 * kPi);
    laplace_sign_ok =
        std::isfinite(analytic_d0) &&
        (same_sign_for_significant_values(analytic_d0, raw.d0_src) ||
         same_sign_for_significant_values(analytic_d0, centroid_d0));
    d0_src = std::isfinite(analytic_d0) ? analytic_d0 : centroid_d0;
    if (!laplace_sign_ok && relative_guard_ok(analytic_d0, centroid_d0, 0.20))
      d0_src = analytic_d0;
  }
  if (k2_info != nullptr) {
    k2_info->laplace_sign_ok = laplace_sign_ok;
    k2_info->laplace_integral = d0_src;
    k2_info->yukawa_smooth_integral = raw.dk_minus_d0_src;
  }

  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  const double dk_src = d0_src + raw.dk_minus_d0_src;
  return {raw.g0_minus_gk, eps_ratio * dk_src - d0_src,
          raw.d0_tgt - inv_eps_ratio * raw.dk_tgt,
          raw.mixedk_minus_mixed0};
}

Panel make_subtriangle_panel(const Panel &source, const Vec3 &a, const Vec3 &b,
                             const Vec3 &c) noexcept;
KernelValues triangle_gauss_subdivided_integral(const Panel &target,
                                                const Panel &source,
                                                double eps_ratio,
                                                double kappa) noexcept;

bool is_node_patch_panel(const Panel &panel) noexcept {
  return panel.node_patch_geometry != nullptr &&
         panel.node_patch_index != static_cast<std::size_t>(-1);
}

bool node_patch_adjacent(const Panel &target, const Panel &source) noexcept {
  if (!is_node_patch_panel(target) || !is_node_patch_panel(source) ||
      target.node_patch_geometry.get() != source.node_patch_geometry.get())
    return false;
  const NodePatchGeometry &geometry = *source.node_patch_geometry;
  const std::size_t source_node = source.node_patch_index;
  if (source_node + 1 >= geometry.adjacent_offsets.size())
    return false;
  const auto begin = geometry.adjacent_nodes.begin() +
                     geometry.adjacent_offsets[source_node];
  const auto end = geometry.adjacent_nodes.begin() +
                   geometry.adjacent_offsets[source_node + 1];
  return std::binary_search(begin, end, target.node_patch_index);
}

Vec3 node_patch_triangle_vertex(const NodePatchGeometry &geometry,
                                const NodePatchTriangle &triangle,
                                int local) noexcept {
  const auto &xyz = geometry.vertex_xyz[
      triangle.vertices[static_cast<std::size_t>(local)]];
  return Vec3{xyz[0], xyz[1], xyz[2]};
}

KernelValues node_patch_near_integral(const Panel &target,
                                      const Panel &source, double eps_ratio,
                                      double kappa, bool self_pair) {
  const NodePatchGeometry &geometry = *source.node_patch_geometry;
  const std::size_t node = source.node_patch_index;
  if (node + 1 >= geometry.incident_offsets.size())
    throw std::runtime_error(
        "PBE TABI: invalid node-patch incident-triangle topology.");
  KernelValues out;
  const std::size_t begin = geometry.incident_offsets[node];
  const std::size_t end = geometry.incident_offsets[node + 1];
  for (std::size_t pos = begin; pos < end; ++pos) {
    const NodePatchTriangle &triangle =
        geometry.triangles[geometry.incident_triangles[pos]];
    int local = -1;
    for (int k = 0; k < 3; ++k) {
      if (triangle.vertices[static_cast<std::size_t>(k)] == node) {
        local = k;
        break;
      }
    }
    if (local < 0)
      throw std::runtime_error(
          "PBE TABI: node-patch incidence does not reference its node.");
    const Vec3 vertex =
        node_patch_triangle_vertex(geometry, triangle, local);
    const Vec3 next =
        node_patch_triangle_vertex(geometry, triangle, (local + 1) % 3);
    const Vec3 previous =
        node_patch_triangle_vertex(geometry, triangle, (local + 2) % 3);
    const Vec3 midpoint_next = 0.5 * (vertex + next);
    const Vec3 midpoint_previous = 0.5 * (previous + vertex);
    const Vec3 centroid = (1.0 / 3.0) * (vertex + next + previous);
    Panel triangle_source = source;
    triangle_source.nx = triangle.nx;
    triangle_source.ny = triangle.ny;
    triangle_source.nz = triangle.nz;
    triangle_source.node_patch_geometry.reset();
    triangle_source.node_patch_index = static_cast<std::size_t>(-1);
    const std::array<Panel, 2> pieces = {
        make_subtriangle_panel(triangle_source, vertex, midpoint_next,
                               centroid),
        make_subtriangle_panel(triangle_source, vertex, centroid,
                               midpoint_previous)};
    for (const Panel &piece : pieces) {
      const KernelValues value =
          self_pair ? duffy_cpv_integral(target, piece, eps_ratio, kappa, true)
                    : triangle_gauss_subdivided_integral(target, piece,
                                                         eps_ratio, kappa);
      out.k1 += value.k1;
      out.k2 += value.k2;
      out.k3 += value.k3;
      out.k4 += value.k4;
    }
  }
  return out;
}

[[maybe_unused]] KernelValues
triangle_gauss_integral(const Panel &target, const Panel &source,
                        double eps_ratio, double kappa) noexcept {
  static constexpr std::array<QuadPoint, 7> points = {
      QuadPoint{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.2250000000000000},
      QuadPoint{0.4701420641051151, 0.4701420641051151, 0.0597158717897698,
                0.1323941527885062},
      QuadPoint{0.4701420641051151, 0.0597158717897698, 0.4701420641051151,
                0.1323941527885062},
      QuadPoint{0.0597158717897698, 0.4701420641051151, 0.4701420641051151,
                0.1323941527885062},
      QuadPoint{0.1012865073234563, 0.1012865073234563, 0.7974269853530873,
                0.1259391805448272},
      QuadPoint{0.1012865073234563, 0.7974269853530873, 0.1012865073234563,
                0.1259391805448272},
      QuadPoint{0.7974269853530873, 0.1012865073234563, 0.1012865073234563,
                0.1259391805448272}};
  KernelValues out;
  if (!source.has_triangle_geometry || source.area <= 0.0)
    return centroid_integral(target, source, eps_ratio, kappa);
  double g0_integral = 0.0;
  double gk_integral = 0.0;
  double d0_src_integral = 0.0;
  double dk_src_integral = 0.0;
  double d0_tgt_integral = 0.0;
  double dk_tgt_integral = 0.0;
  double mixed0_integral = 0.0;
  double mixedk_integral = 0.0;
  for (const QuadPoint &q : points) {
    const double sx = q.l0 * source.v0x + q.l1 * source.v1x +
                      q.l2 * source.v2x;
    const double sy = q.l0 * source.v0y + q.l1 * source.v1y +
                      q.l2 * source.v2y;
    const double sz = q.l0 * source.v0z + q.l1 * source.v1z +
                      q.l2 * source.v2z;
    const double w = source.area * q.weight;
    const double dx = target.x - sx;
    const double dy = target.y - sy;
    const double dz = target.z - sz;
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (r <= 1.0e-14)
      continue;
    g0_integral += w * green(r, 0.0);
    gk_integral += w * green(r, kappa);
    d0_src_integral +=
        w * source_normal_green(target.x, target.y, target.z, sx, sy, sz,
                                source.nx, source.ny, source.nz, 0.0);
    dk_src_integral +=
        w * source_normal_green(target.x, target.y, target.z, sx, sy, sz,
                                source.nx, source.ny, source.nz, kappa);
    d0_tgt_integral += w * target_normal_green(target, sx, sy, sz, 0.0);
    dk_tgt_integral += w * target_normal_green(target, sx, sy, sz, kappa);
    mixed0_integral +=
        w * mixed_normal_green(target, sx, sy, sz, source.nx, source.ny,
                               source.nz, 0.0);
    mixedk_integral +=
        w * mixed_normal_green(target, sx, sy, sz, source.nx, source.ny,
                               source.nz, kappa);
  }
  auto solid_angle = [](const Panel &target, const Panel &source) noexcept {
    const double ax = source.v0x - target.x;
    const double ay = source.v0y - target.y;
    const double az = source.v0z - target.z;
    const double bx = source.v1x - target.x;
    const double by = source.v1y - target.y;
    const double bz = source.v1z - target.z;
    const double cx = source.v2x - target.x;
    const double cy = source.v2y - target.y;
    const double cz = source.v2z - target.z;
    const double la = std::sqrt(ax * ax + ay * ay + az * az);
    const double lb = std::sqrt(bx * bx + by * by + bz * bz);
    const double lc = std::sqrt(cx * cx + cy * cy + cz * cz);
    if (la <= 1.0e-14 || lb <= 1.0e-14 || lc <= 1.0e-14)
      return 0.0;
    const double cross_x = by * cz - bz * cy;
    const double cross_y = bz * cx - bx * cz;
    const double cross_z = bx * cy - by * cx;
    const double numerator = ax * cross_x + ay * cross_y + az * cross_z;
    const double ab = ax * bx + ay * by + az * bz;
    const double bc = bx * cx + by * cy + bz * cz;
    const double ca = cx * ax + cy * ay + cz * az;
    const double denominator = la * lb * lc + ab * lc + bc * la + ca * lb;
    return 2.0 * std::atan2(numerator, denominator);
  };
  const double omega = solid_angle(target, source);
  if (std::isfinite(omega) && std::abs(omega) > 1.0e-14) {
    double analytic_d0_src = -omega / (4.0 * kPi);
    if (analytic_d0_src * d0_src_integral < 0.0)
      analytic_d0_src = -analytic_d0_src;
    d0_src_integral = analytic_d0_src;
  }
  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  out.k1 = g0_integral - gk_integral;
  out.k2 = eps_ratio * dk_src_integral - d0_src_integral;
  out.k3 = d0_tgt_integral - inv_eps_ratio * dk_tgt_integral;
  out.k4 = mixedk_integral - mixed0_integral;
  return out;
}

Panel make_subtriangle_panel(const Panel &source, const Vec3 &a, const Vec3 &b,
                             const Vec3 &c) noexcept {
  Panel sub = source;
  sub.v0x = a.x;
  sub.v0y = a.y;
  sub.v0z = a.z;
  sub.v1x = b.x;
  sub.v1y = b.y;
  sub.v1z = b.z;
  sub.v2x = c.x;
  sub.v2y = c.y;
  sub.v2z = c.z;
  const Vec3 normal_area = cross(b - a, c - a);
  sub.area = 0.5 * norm(normal_area);
  sub.x = (a.x + b.x + c.x) / 3.0;
  sub.y = (a.y + b.y + c.y) / 3.0;
  sub.z = (a.z + b.z + c.z) / 3.0;
  sub.has_triangle_geometry = source.has_triangle_geometry && sub.area > 0.0;
  sub.v0_index = static_cast<std::size_t>(-1);
  sub.v1_index = static_cast<std::size_t>(-1);
  sub.v2_index = static_cast<std::size_t>(-1);
  return sub;
}

KernelValues triangle_gauss_subdivided_integral(const Panel &target,
                                                const Panel &source,
                                                double eps_ratio,
                                                double kappa) noexcept {
  if (!source.has_triangle_geometry || source.area <= 0.0)
    return triangle_gauss_integral(target, source, eps_ratio, kappa);
  const Vec3 a = panel_vertex(source, 0);
  const Vec3 b = panel_vertex(source, 1);
  const Vec3 c = panel_vertex(source, 2);
  const Vec3 ab = 0.5 * (a + b);
  const Vec3 bc = 0.5 * (b + c);
  const Vec3 ca = 0.5 * (c + a);
  const std::array<Panel, 4> first_level = {
      make_subtriangle_panel(source, a, ab, ca),
      make_subtriangle_panel(source, ab, b, bc),
      make_subtriangle_panel(source, ca, bc, c),
      make_subtriangle_panel(source, ab, bc, ca)};
  KernelValues out;
  for (const Panel &sub : first_level) {
    if (!sub.has_triangle_geometry)
      continue;
    const Vec3 sa = panel_vertex(sub, 0);
    const Vec3 sb = panel_vertex(sub, 1);
    const Vec3 sc = panel_vertex(sub, 2);
    const Vec3 sab = 0.5 * (sa + sb);
    const Vec3 sbc = 0.5 * (sb + sc);
    const Vec3 sca = 0.5 * (sc + sa);
    const std::array<Panel, 4> second_level = {
        make_subtriangle_panel(sub, sa, sab, sca),
        make_subtriangle_panel(sub, sab, sb, sbc),
        make_subtriangle_panel(sub, sca, sbc, sc),
        make_subtriangle_panel(sub, sab, sbc, sca)};
    for (const Panel &leaf : second_level) {
      if (!leaf.has_triangle_geometry)
        continue;
      const KernelValues k =
          triangle_gauss_integral(target, leaf, eps_ratio, kappa);
      out.k1 += k.k1;
      out.k2 += k.k2;
      out.k3 += k.k3;
      out.k4 += k.k4;
    }
  }
  return out;
}

void add_stats(TreecodeStats *stats, bool self_pair, bool one_ring,
               bool duffy, bool centroid, bool fallback,
               std::size_t subdivisions) noexcept {
  if (stats == nullptr)
    return;
  if (self_pair)
    ++stats->near_self_pairs;
  if (one_ring)
    ++stats->near_one_ring_pairs;
  if (duffy)
    ++stats->near_duffy_pairs;
  if (centroid)
    ++stats->near_centroid_pairs;
  if (fallback)
    ++stats->near_fallbacks;
  if (subdivisions > 0) {
    const double count = static_cast<double>(stats->near_duffy_pairs);
    stats->duffy_subdivisions_mean =
        count <= 1.0
            ? static_cast<double>(subdivisions)
            : (stats->duffy_subdivisions_mean * (count - 1.0) +
               static_cast<double>(subdivisions)) /
                  count;
    stats->duffy_subdivisions_max =
        std::max(stats->duffy_subdivisions_max, subdivisions);
  }
}

std::size_t geometric_gap_bucket(double gap_ratio) noexcept {
  if (gap_ratio < 0.25)
    return 0;
  if (gap_ratio < 0.50)
    return 1;
  if (gap_ratio <= 0.75)
    return 2;
  return 3;
}

std::size_t geometric_center_bucket(double center_ratio) noexcept {
  if (center_ratio < 1.0)
    return 0;
  if (center_ratio < 2.0)
    return 1;
  return 2;
}

std::size_t geometric_kernel_bucket(const std::string &kernel) noexcept {
  if (kernel == "gauss7")
    return 1;
  if (kernel == "subdivided")
    return 2;
  return 0;
}

std::size_t geometric_audit_index(double gap_ratio, double center_ratio,
                                  const std::string &kernel) noexcept {
  const std::size_t gap = geometric_gap_bucket(gap_ratio);
  const std::size_t center = geometric_center_bucket(center_ratio);
  const std::size_t kernel_index = geometric_kernel_bucket(kernel);
  return (gap * kGeometricCenterBucketCount + center) *
             kGeometricKernelCount +
         kernel_index;
}

void add_geometric_buffer_stats(TreecodeStats *stats,
                                const GeometricNearDecision &decision,
                                const std::string &kernel) noexcept {
  if (stats == nullptr || !decision.active)
    return;
  const std::size_t old_count = stats->near_geometric_buffer_pairs;
  ++stats->near_geometric_buffer_pairs;
  if (kernel == "subdivided")
    ++stats->near_geometric_buffer_adaptive_pairs;
  else if (kernel == "gauss7")
    ++stats->near_geometric_buffer_gauss_pairs;
  else
    ++stats->near_geometric_buffer_centroid_pairs;
  if (old_count == 0 || stats->near_geometric_buffer_min_gap <= 0.0) {
    stats->near_geometric_buffer_min_gap = decision.gap;
  } else {
    stats->near_geometric_buffer_min_gap =
        std::min(stats->near_geometric_buffer_min_gap, decision.gap);
  }
  const double total = static_cast<double>(old_count + 1);
  stats->near_geometric_buffer_mean_gap_ratio =
      (stats->near_geometric_buffer_mean_gap_ratio *
           static_cast<double>(old_count) +
       decision.gap_ratio) /
      total;
  stats->near_geometric_buffer_max_edge_mean =
      (stats->near_geometric_buffer_max_edge_mean *
           static_cast<double>(old_count) +
       decision.max_edge) /
      total;
  const std::size_t index = geometric_audit_index(
      decision.gap_ratio, decision.center_ratio, kernel);
  if (index < stats->geometric_near_audit.size()) {
    auto &bucket = stats->geometric_near_audit[index];
    if (bucket.pairs == 0 || bucket.min_gap <= 0.0)
      bucket.min_gap = decision.gap;
    else
      bucket.min_gap = std::min(bucket.min_gap, decision.gap);
    bucket.max_gap_ratio =
        std::max(bucket.max_gap_ratio, decision.gap_ratio);
    bucket.sum_gap_ratio += decision.gap_ratio;
    bucket.sum_center_ratio += decision.center_ratio;
    ++bucket.pairs;
  }
}

void record_block_decisions(TreecodeStats *stats,
                            const AcceptedKernelValues &accepted,
                            bool one_ring, std::size_t pair_type,
                            std::size_t orientation) noexcept {
  if (stats == nullptr)
    return;
  auto record = [](std::size_t &accepted_count, std::size_t &rejected_count,
                   std::string &policy,
                   const KernelBlockDecision &decision) noexcept {
    if (decision.accepted) {
      ++accepted_count;
    } else {
      ++rejected_count;
    }
    if (accepted_count > 0 && rejected_count > 0)
      policy = "mixed_guarded";
    else
      policy = decision.accepted ? decision.accepted_policy
                                 : decision.rejected_policy;
  };
  record(stats->near_k1_accepted, stats->near_k1_rejected,
         stats->near_k1_policy, accepted.block[0]);
  record(stats->near_k2_accepted, stats->near_k2_rejected,
         stats->near_k2_policy, accepted.block[1]);
  record(stats->near_k3_accepted, stats->near_k3_rejected,
         stats->near_k3_policy, accepted.block[2]);
  record(stats->near_k4_accepted, stats->near_k4_rejected,
         stats->near_k4_policy, accepted.block[3]);
  for (std::size_t block = 0; block < accepted.block.size(); ++block) {
    const std::size_t bucket_index =
        near_block_bucket_index(block, pair_type, orientation);
    if (bucket_index >= stats->near_block_audit.size())
      continue;
    auto &bucket = stats->near_block_audit[bucket_index];
    if (accepted.block[block].accepted)
      ++bucket.accepted;
    else
      ++bucket.rejected;
    if (accepted.candidate_worse_than_centroid[block])
      ++bucket.candidate_worse_than_centroid;
    if (!accepted.block[block].accepted) {
      switch (accepted.block[block].reject_reason) {
      case RejectReason::ReferenceNonconverged:
        ++bucket.reference_nonconverged;
        break;
      case RejectReason::OrientationSignMismatch:
        ++bucket.orientation_sign_mismatch;
        break;
      case RejectReason::StrictReferenceReject:
        ++bucket.strict_ref_reject;
        break;
      case RejectReason::CandidateWorseThanCentroid:
        break;
      case RejectReason::None:
        ++bucket.strict_ref_reject;
        break;
      }
    }
    if (accepted.candidate_ref_rel[block] > bucket.max_candidate_ref_rel)
      bucket.max_candidate_ref_rel = accepted.candidate_ref_rel[block];
    bucket.sum_candidate_ref_rel += accepted.candidate_ref_rel[block];
    if (accepted.centroid_ref_rel[block] > bucket.max_centroid_ref_rel)
      bucket.max_centroid_ref_rel = accepted.centroid_ref_rel[block];
    if (block == 1) {
      if (accepted.k2_laplace_sign_ok[block])
        ++bucket.k2_laplace_sign_ok;
      if (accepted.k2_yukawa_smooth_ref_rel[block] >
          bucket.max_k2_yukawa_smooth_ref_rel)
        bucket.max_k2_yukawa_smooth_ref_rel =
            accepted.k2_yukawa_smooth_ref_rel[block];
      if (accepted.k2_orientation_pair_ref_rel[block] >
          bucket.max_k2_orientation_pair_ref_rel)
        bucket.max_k2_orientation_pair_ref_rel =
            accepted.k2_orientation_pair_ref_rel[block];
    }
    ++bucket.sample_count;
  }
  if (!one_ring)
    return;
  auto record_one_ring = [](std::size_t &accepted_count,
                            std::size_t &rejected_count,
                            const KernelBlockDecision &decision) noexcept {
    if (decision.accepted)
      ++accepted_count;
    else
      ++rejected_count;
  };
  record_one_ring(stats->near_one_ring_k1_accepted,
                  stats->near_one_ring_k1_rejected, accepted.block[0]);
  record_one_ring(stats->near_one_ring_k2_accepted,
                  stats->near_one_ring_k2_rejected, accepted.block[1]);
  record_one_ring(stats->near_one_ring_k3_accepted,
                  stats->near_one_ring_k3_rejected, accepted.block[2]);
  record_one_ring(stats->near_one_ring_k4_accepted,
                  stats->near_one_ring_k4_rejected, accepted.block[3]);
}

void apply_kernel(JcpNearAccum &acc, const KernelValues &kernel,
                  const std::vector<double> &unknowns, std::size_t n,
                  std::size_t j, bool coupled) noexcept {
  const double phi_j = unknowns[j];
  if (!coupled) {
    acc.phi += kernel.k2 * phi_j;
    return;
  }
  const double phin_j = unknowns[n + j];
  acc.phi += kernel.k1 * phin_j + kernel.k2 * phi_j;
  acc.phin += kernel.k3 * phin_j + kernel.k4 * phi_j;
}

bool has_all_triangle_geometry(const std::vector<Panel> &panels) noexcept {
  if (panels.empty())
    return false;
  for (const Panel &panel : panels) {
    if (!panel.has_triangle_geometry)
      return false;
  }
  return true;
}

std::string kernel_policy_for_mode(const std::string &near_integration) {
  if (near_integration == "analytic_singular")
    return "self_asi";
  if (near_integration == "duffy_adaptive")
    return "one_ring_cpv";
  return "centroid";
}

AcceptedKernelValues guarded_one_ring_cpv_integral_decision(
    const Panel &target, const Panel &source, double eps_ratio,
    double kappa) noexcept;

AcceptedKernelValues guarded_one_ring_cpv_integral_decision(
    const Panel &target, const Panel &source, double eps_ratio,
    double kappa) noexcept {
  const KernelValues centroid =
      centroid_integral(target, source, eps_ratio, kappa);
  K2CandidateInfo k2_info;
  const KernelValues duffy =
      duffy_cpv_integral(target, source, eps_ratio, kappa, false, &k2_info);
  const KernelValues reference_low =
      triangle_gauss_integral(target, source, eps_ratio, kappa);
  const KernelValues reference =
      triangle_gauss_subdivided_integral(target, source, eps_ratio, kappa);
  Panel flipped_source = source;
  flipped_source.nx = -flipped_source.nx;
  flipped_source.ny = -flipped_source.ny;
  flipped_source.nz = -flipped_source.nz;
  const KernelValues flipped_reference =
      triangle_gauss_subdivided_integral(target, flipped_source, eps_ratio,
                                         kappa);
  const double k2_orientation_pair_ref_rel =
      relative_error_value(reference.k2, -flipped_reference.k2);
  const double k2_laplace = k2_info.laplace_integral;
  const double k2_reference_smooth =
      (reference.k2 - (eps_ratio - 1.0) * k2_laplace) /
      std::max(1.0e-12, eps_ratio);
  const double k2_yukawa_smooth_ref_rel =
      relative_error_value(k2_info.yukawa_smooth_integral,
                           k2_reference_smooth);
  const std::array<double, 4> candidate_values = {duffy.k1, duffy.k2,
                                                  duffy.k3, duffy.k4};
  const std::array<double, 4> centroid_values = {centroid.k1, centroid.k2,
                                                 centroid.k3, centroid.k4};
  const std::array<double, 4> reference_values = {
      reference.k1, reference.k2, reference.k3, reference.k4};
  const std::array<double, 4> reference_low_values = {
      reference_low.k1, reference_low.k2, reference_low.k3, reference_low.k4};
  const std::array<double, 4> reference_convergence_limits = {0.10, 0.50,
                                                             0.50, 0.50};
  std::array<bool, 4> reference_converged{};
  std::array<bool, 4> orientation_mismatch{};
  std::array<bool, 4> accept{};
  for (std::size_t block = 0; block < accept.size(); ++block) {
    const double ref_ref_rel =
        relative_error_value(reference_low_values[block],
                             reference_values[block]);
    reference_converged[block] =
        finite_kernel_values(reference_low) && finite_kernel_values(reference) &&
        std::isfinite(ref_ref_rel) &&
        ref_ref_rel <= reference_convergence_limits[block];
    orientation_mismatch[block] =
        block == 1 &&
        (!k2_info.laplace_sign_ok ||
         !same_sign_for_significant_values(candidate_values[block],
                                           reference_values[block]) ||
         !relative_guard_ok(reference_values[block],
                            -flipped_reference.k2, 0.02));
    const bool candidate_worse =
        relative_error_value(candidate_values[block], reference_values[block]) >
        relative_error_value(centroid_values[block], reference_values[block]) +
            1.0e-12;
    const double candidate_ref_rel =
        relative_error_value(candidate_values[block], reference_values[block]);
    const double centroid_ref_rel =
        relative_error_value(centroid_values[block], reference_values[block]);
    const bool candidate_improves_centroid =
        candidate_ref_rel <= 0.95 * centroid_ref_rel;
    accept[block] =
        block == 1 && reference_converged[block] &&
        !orientation_mismatch[block] && !candidate_worse &&
        candidate_improves_centroid &&
        candidate_ref_rel <= 1.0e-8 &&
        relative_guard_ok(candidate_values[block], reference_values[block],
                          0.02);
  }
  AcceptedKernelValues out;
  out.kernel = {accept[0] ? duffy.k1 : centroid.k1,
                accept[1] ? duffy.k2 : centroid.k2,
                accept[2] ? duffy.k3 : centroid.k3,
                accept[3] ? duffy.k4 : centroid.k4};
  const auto reject_reason_for_block =
      [&](std::size_t block) noexcept -> RejectReason {
    if (accept[block])
      return RejectReason::None;
    if (!reference_converged[block])
      return RejectReason::ReferenceNonconverged;
    if (orientation_mismatch[block])
      return RejectReason::OrientationSignMismatch;
    if (relative_error_value(candidate_values[block], reference_values[block]) >
        relative_error_value(centroid_values[block], reference_values[block]) +
            1.0e-12) {
      return RejectReason::CandidateWorseThanCentroid;
    }
    return RejectReason::StrictReferenceReject;
  };
  const auto reject_policy_for_reason = [](RejectReason reason) noexcept {
    switch (reason) {
    case RejectReason::ReferenceNonconverged:
      return "reject_reference_nonconverged";
    case RejectReason::OrientationSignMismatch:
      return "reject_orientation_sign_mismatch";
    case RejectReason::CandidateWorseThanCentroid:
      return "reject_candidate_worse_than_centroid";
    case RejectReason::StrictReferenceReject:
      return "reject_strict_reference";
    case RejectReason::None:
      return "centroid_guarded";
    }
    return "centroid_guarded";
  };
  const std::array<RejectReason, 4> reasons = {
      reject_reason_for_block(0), reject_reason_for_block(1),
      reject_reason_for_block(2), reject_reason_for_block(3)};
  out.block = {
      KernelBlockDecision{accept[0], "duffy_g0_gk",
                          reject_policy_for_reason(reasons[0]), reasons[0]},
      KernelBlockDecision{accept[1], "duffy_source_normal",
                          reject_policy_for_reason(reasons[1]), reasons[1]},
      KernelBlockDecision{accept[2], "duffy_target_normal",
                          reject_policy_for_reason(reasons[2]), reasons[2]},
      KernelBlockDecision{accept[3], "duffy_mixed_normal",
                          reject_policy_for_reason(reasons[3]), reasons[3]}};
  for (std::size_t block = 0; block < out.block.size(); ++block) {
    out.candidate_ref_rel[block] =
        relative_error_value(candidate_values[block], reference_values[block]);
    out.centroid_ref_rel[block] =
        relative_error_value(centroid_values[block], reference_values[block]);
    out.candidate_worse_than_centroid[block] =
        out.candidate_ref_rel[block] > out.centroid_ref_rel[block] + 1.0e-12;
  }
  out.k2_laplace_sign_ok[1] = k2_info.laplace_sign_ok;
  out.k2_yukawa_smooth_ref_rel[1] = k2_yukawa_smooth_ref_rel;
  out.k2_orientation_pair_ref_rel[1] = k2_orientation_pair_ref_rel;
  return out;
}

} // namespace

std::string resolve_near_integration_mode(const std::string &requested,
                                          const std::vector<Panel> &panels) {
  const std::string mode = requested.empty() ? std::string("auto") : requested;
  const bool node_patch_ready =
      !panels.empty() &&
      std::all_of(panels.begin(), panels.end(), is_node_patch_panel);
  const bool triangle_ready =
      has_all_triangle_geometry(panels) || node_patch_ready;
  if (mode == "auto")
    return "centroid";
  if (mode == "centroid")
    return mode;
  if ((mode == "analytic_singular" || mode == "duffy_adaptive") &&
      !triangle_ready) {
    throw std::runtime_error(
        "PBE TABI: -pbe-tabi-near-integration requires triangle surface "
        "geometry; use -pbe-tabi-geometry=triangle or near-integration=auto.");
  }
  if (mode != "analytic_singular" && mode != "duffy_adaptive")
    throw std::runtime_error(
        "PBE TABI: invalid near integration mode; expected auto, centroid, "
        "analytic_singular, or duffy_adaptive.");
  return mode;
}

std::string resolve_quad_buffer_mode(const std::string &requested,
                                     const std::string &near_integration,
                                     const std::vector<Panel> &panels) {
  const std::string mode = requested.empty() ? std::string("auto") : requested;
  if (mode == "off")
    return "off";
  (void)near_integration;
  if (mode == "auto")
    return "off";
  if (mode != "gauss7" && mode != "adaptive")
    throw std::runtime_error(
        "PBE TABI: invalid quad-buffer mode; expected auto, off, gauss7, "
        "or adaptive.");
  const bool node_patch_ready =
      !panels.empty() &&
      std::all_of(panels.begin(), panels.end(), is_node_patch_panel);
  if (!has_all_triangle_geometry(panels) && !node_patch_ready) {
    throw std::runtime_error(
        "PBE TABI: -pbe-tabi-quad-buffer requires triangle surface geometry; "
        "use -pbe-tabi-geometry=triangle or quad-buffer=auto/off.");
  }
  return mode;
}

GeometricNearDecision classify_geometric_near_buffer_pair(
    const Panel &target, const Panel &source,
    const std::string &quad_buffer) noexcept {
  GeometricNearDecision decision;
  if (quad_buffer == "off" || quad_buffer.empty() ||
      !target.has_triangle_geometry || !source.has_triangle_geometry) {
    return decision;
  }
  const double target_edge = panel_max_edge_length(target);
  const double source_edge = panel_max_edge_length(source);
  const double max_edge = std::max(target_edge, source_edge);
  if (!std::isfinite(max_edge) || max_edge <= 1.0e-12)
    return decision;
  const double threshold = 0.75 * max_edge;
  const Vec3 target_center = make_vec(target.x, target.y, target.z);
  const Vec3 source_center = make_vec(source.x, source.y, source.z);
  const double center_distance =
      std::sqrt(std::max(0.0, distance_sq(target_center, source_center)));
  const double radius_sum =
      panel_centroid_radius(target) + panel_centroid_radius(source);
  if (center_distance - radius_sum > threshold)
    return decision;
  const double gap = std::sqrt(std::max(
      0.0, triangle_triangle_distance_sq(target, source)));
  const double gap_ratio = gap / max_edge;
  if (!std::isfinite(gap_ratio) || gap_ratio > 0.75)
    return decision;
  decision.active = true;
  decision.gap = gap;
  decision.gap_ratio = gap_ratio;
  decision.center_ratio = center_distance / max_edge;
  decision.max_edge = max_edge;
  decision.adaptive = quad_buffer == "adaptive" && gap_ratio <= 0.35;
  return decision;
}

JcpNearPairKernel evaluate_jcp_near_pair_kernel(
    const Panel &target, std::size_t target_index, const Panel &source,
    std::size_t source_index, double eps_ratio, double kappa, bool coupled,
    const std::string &near_integration, const std::string &quad_buffer,
    TreecodeStats *stats) {
  JcpNearPairKernel selected;
  selected.self_pair = source_index == target_index;
  const bool node_patch_pair =
      is_node_patch_panel(target) && is_node_patch_panel(source);
  const int shared_vertices =
      selected.self_pair ? 3 : shared_vertex_count(target, source);
  selected.one_ring =
      !selected.self_pair &&
      (node_patch_pair ? node_patch_adjacent(target, source)
                       : shared_vertices > 0);
  KernelValues kernel;
  bool used_duffy = false;
  bool used_centroid = false;
  bool fallback = false;
  std::size_t subdivisions = 0;
  GeometricNearDecision geometric_decision;
  std::string geometric_kernel;
  AcceptedKernelValues block_decision;
  bool has_block_decision = false;

  if (node_patch_pair && !coupled) {
    if (selected.self_pair)
      return selected;
    kernel = selected.one_ring
                 ? node_patch_near_integral(target, source, eps_ratio, 0.0,
                                            false)
                 : centroid_integral(target, source, eps_ratio, 0.0);
    used_duffy = selected.one_ring;
    used_centroid = !selected.one_ring;
    subdivisions = selected.one_ring ? 32 : 0;
  } else if (node_patch_pair && selected.self_pair) {
    kernel = self_singular_integral(source, kappa);
    used_duffy = true;
    subdivisions = 1;
  } else if (!coupled) {
    if (selected.self_pair)
      return selected;
    kernel = centroid_integral(target, source, eps_ratio, 0.0);
    used_centroid = true;
  } else if (near_integration == "centroid") {
    if (selected.self_pair)
      return selected;
    kernel = centroid_integral(target, source, eps_ratio, kappa);
    used_centroid = true;
  } else if (selected.self_pair) {
    if (!target.has_triangle_geometry || !source.has_triangle_geometry)
      throw std::runtime_error(
          "PBE TABI: analytic near integration requires triangle geometry "
          "for self panels.");
    block_decision = self_singular_decision(source, kappa);
    kernel = block_decision.kernel;
    has_block_decision = true;
    used_duffy = true;
    subdivisions = 3;
  } else if (selected.one_ring && near_integration == "duffy_adaptive") {
    if (!target.has_triangle_geometry || !source.has_triangle_geometry)
      throw std::runtime_error(
          "PBE TABI: Duffy one-ring near integration requires triangle "
          "geometry.");
    block_decision =
        guarded_one_ring_cpv_integral_decision(target, source, eps_ratio,
                                               kappa);
    kernel = block_decision.kernel;
    has_block_decision = true;
    used_duffy = true;
    subdivisions = 8;
  } else if (selected.one_ring) {
    kernel = centroid_integral(target, source, eps_ratio, kappa);
    used_centroid = true;
  } else if ((geometric_decision = classify_geometric_near_buffer_pair(
                  target, source, quad_buffer))
                 .active) {
    selected.geometric_buffer = true;
    if (geometric_decision.adaptive) {
      kernel =
          triangle_gauss_subdivided_integral(target, source, eps_ratio, kappa);
      geometric_kernel = "subdivided";
      subdivisions = 16;
    } else {
      kernel = triangle_gauss_integral(target, source, eps_ratio, kappa);
      geometric_kernel = "gauss7";
      subdivisions = 1;
    }
  } else {
    kernel = centroid_integral(target, source, eps_ratio, kappa);
    used_centroid = true;
  }

  if (!finite_kernel(kernel))
    throw std::runtime_error(
        "PBE TABI: non-finite analytic near-integration kernel.");
  selected.k1 = kernel.k1;
  selected.k2 = kernel.k2;
  selected.k3 = kernel.k3;
  selected.k4 = kernel.k4;
  selected.used_duffy = used_duffy;
  selected.used_centroid = used_centroid;
  selected.valid = true;
  selected.policy = selected.geometric_buffer
                        ? geometric_kernel
                        : (used_duffy ? std::string("duffy")
                                      : (used_centroid ? std::string("centroid")
                                                       : std::string("near")));
  add_stats(stats, selected.self_pair, selected.one_ring, used_duffy,
            used_centroid, fallback, subdivisions);
  if (geometric_decision.active)
    add_geometric_buffer_stats(stats, geometric_decision, geometric_kernel);
  if (has_block_decision) {
    const std::size_t pair_type =
        classify_pair_type(selected.self_pair, shared_vertices);
    const std::size_t orientation =
        selected.self_pair ? 0 : classify_orientation(target, source);
    record_block_decisions(stats, block_decision, selected.one_ring, pair_type,
                           orientation);
  }
  if (stats != nullptr) {
    stats->near_integration = near_integration;
    stats->quad_buffer = quad_buffer;
    stats->near_kernel_policy = kernel_policy_for_mode(near_integration);
  }
  return selected;
}

JcpNearAccum accumulate_jcp_near_ordered(
    const Panel &target, std::size_t target_index,
    const std::vector<Panel> &panels, const std::vector<double> &unknowns,
    double eps_ratio, double kappa, bool coupled,
    const std::vector<std::size_t> &order, std::size_t begin,
    std::size_t count, const std::string &near_integration,
    const std::string &quad_buffer,
    TreecodeStats *stats) {
  JcpNearAccum acc;
  const std::size_t n = panels.size();
  const std::size_t end = begin + count;
  for (std::size_t pos = begin; pos < end; ++pos) {
    const std::size_t j = order[pos];
    const Panel &source = panels[j];
    const JcpNearPairKernel selected = evaluate_jcp_near_pair_kernel(
        target, target_index, source, j, eps_ratio, kappa, coupled,
        near_integration, quad_buffer, stats);
    if (!selected.valid)
      continue;
    KernelValues kernel;
    kernel.k1 = selected.k1;
    kernel.k2 = selected.k2;
    kernel.k3 = selected.k3;
    kernel.k4 = selected.k4;
    apply_kernel(acc, kernel, unknowns, n, j, coupled);
  }
  if (stats != nullptr) {
    stats->near_integration = near_integration;
    stats->quad_buffer = quad_buffer;
    stats->near_kernel_policy = kernel_policy_for_mode(near_integration);
  }
  return acc;
}

} // namespace TitanPBE::Tabi
