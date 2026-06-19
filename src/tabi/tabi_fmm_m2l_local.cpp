#include "tabi_fmm_m2l_local.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace titan::tabi::fmm_m2l {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr int kJetSide = kCartesianLocalKernelOrder + 1;
constexpr std::size_t kJetSize =
    static_cast<std::size_t>(kJetSide * kJetSide * kJetSide);

struct MultiIndex {
  int x{0};
  int y{0};
  int z{0};
};

struct Jet8 {
  std::array<double, kJetSize> coeff{};
};

std::size_t cube_index(int x, int y, int z) noexcept {
  return static_cast<std::size_t>((x * kJetSide + y) * kJetSide + z);
}

const std::vector<MultiIndex> &indices_through(int order) {
  static const std::vector<MultiIndex> all = [] {
    std::vector<MultiIndex> result;
    result.reserve(165);
    for (int degree = 0; degree <= kCartesianLocalKernelOrder; ++degree) {
      for (int x = degree; x >= 0; --x) {
        for (int y = degree - x; y >= 0; --y)
          result.push_back({x, y, degree - x - y});
      }
    }
    return result;
  }();
  static thread_local std::vector<MultiIndex> active;
  active.clear();
  for (const MultiIndex &index : all) {
    if (index.x + index.y + index.z <= order)
      active.push_back(index);
  }
  return active;
}

double binomial(int n, int k) noexcept {
  if (k < 0 || k > n)
    return 0.0;
  double value = 1.0;
  for (int i = 1; i <= k; ++i)
    value = value * static_cast<double>(n - k + i) / static_cast<double>(i);
  return value;
}

double power(double value, int exponent) noexcept {
  double result = 1.0;
  for (int i = 0; i < exponent; ++i)
    result *= value;
  return result;
}

Jet8 constant(double value) noexcept {
  Jet8 result;
  result.coeff[cube_index(0, 0, 0)] = value;
  return result;
}

Jet8 variable(double value, int axis) noexcept {
  Jet8 result = constant(value);
  result.coeff[cube_index(axis == 0 ? 1 : 0, axis == 1 ? 1 : 0,
                           axis == 2 ? 1 : 0)] = 1.0;
  return result;
}

Jet8 add(const Jet8 &a, const Jet8 &b) noexcept {
  Jet8 result;
  for (const MultiIndex &index : indices_through(kCartesianLocalKernelOrder)) {
    const std::size_t i = cube_index(index.x, index.y, index.z);
    result.coeff[i] = a.coeff[i] + b.coeff[i];
  }
  return result;
}

Jet8 subtract(const Jet8 &a, const Jet8 &b) noexcept {
  Jet8 result;
  for (const MultiIndex &index : indices_through(kCartesianLocalKernelOrder)) {
    const std::size_t i = cube_index(index.x, index.y, index.z);
    result.coeff[i] = a.coeff[i] - b.coeff[i];
  }
  return result;
}

Jet8 scale(const Jet8 &a, double factor) noexcept {
  Jet8 result;
  for (const MultiIndex &index : indices_through(kCartesianLocalKernelOrder)) {
    const std::size_t i = cube_index(index.x, index.y, index.z);
    result.coeff[i] = factor * a.coeff[i];
  }
  return result;
}

Jet8 multiply(const Jet8 &a, const Jet8 &b) noexcept {
  Jet8 result;
  for (const MultiIndex &index : indices_through(kCartesianLocalKernelOrder)) {
    double value = 0.0;
    for (int x = 0; x <= index.x; ++x) {
      for (int y = 0; y <= index.y; ++y) {
        for (int z = 0; z <= index.z; ++z) {
          value += a.coeff[cube_index(x, y, z)] *
                   b.coeff[cube_index(index.x - x, index.y - y,
                                      index.z - z)];
        }
      }
    }
    result.coeff[cube_index(index.x, index.y, index.z)] = value;
  }
  return result;
}

Jet8 reciprocal(const Jet8 &a) noexcept {
  Jet8 result;
  const double a0 = a.coeff[cube_index(0, 0, 0)];
  if (!std::isfinite(a0) || std::abs(a0) <= 1.0e-30)
    return result;
  result.coeff[cube_index(0, 0, 0)] = 1.0 / a0;
  for (const MultiIndex &index : indices_through(kCartesianLocalKernelOrder)) {
    if (index.x == 0 && index.y == 0 && index.z == 0)
      continue;
    double sum = 0.0;
    for (int x = 0; x <= index.x; ++x) {
      for (int y = 0; y <= index.y; ++y) {
        for (int z = 0; z <= index.z; ++z) {
          if (x == 0 && y == 0 && z == 0)
            continue;
          sum += a.coeff[cube_index(x, y, z)] *
                 result.coeff[cube_index(index.x - x, index.y - y,
                                          index.z - z)];
        }
      }
    }
    result.coeff[cube_index(index.x, index.y, index.z)] = -sum / a0;
  }
  return result;
}

Jet8 square_root(const Jet8 &a) noexcept {
  Jet8 result;
  const double a0 = a.coeff[cube_index(0, 0, 0)];
  if (!std::isfinite(a0) || a0 <= 0.0)
    return result;
  const double root = std::sqrt(a0);
  result.coeff[cube_index(0, 0, 0)] = root;
  for (const MultiIndex &index : indices_through(kCartesianLocalKernelOrder)) {
    if (index.x == 0 && index.y == 0 && index.z == 0)
      continue;
    double sum = 0.0;
    for (int x = 0; x <= index.x; ++x) {
      for (int y = 0; y <= index.y; ++y) {
        for (int z = 0; z <= index.z; ++z) {
          if ((x == 0 && y == 0 && z == 0) ||
              (x == index.x && y == index.y && z == index.z))
            continue;
          sum += result.coeff[cube_index(x, y, z)] *
                 result.coeff[cube_index(index.x - x, index.y - y,
                                          index.z - z)];
        }
      }
    }
    const std::size_t i = cube_index(index.x, index.y, index.z);
    result.coeff[i] = (a.coeff[i] - sum) / (2.0 * root);
  }
  return result;
}

Jet8 exponential(const Jet8 &a) noexcept {
  Jet8 result;
  result.coeff[cube_index(0, 0, 0)] =
      std::exp(a.coeff[cube_index(0, 0, 0)]);
  for (const MultiIndex &index : indices_through(kCartesianLocalKernelOrder)) {
    if (index.x == 0 && index.y == 0 && index.z == 0)
      continue;
    int axis = index.x > 0 ? 0 : (index.y > 0 ? 1 : 2);
    const int gx = index.x - (axis == 0 ? 1 : 0);
    const int gy = index.y - (axis == 1 ? 1 : 0);
    const int gz = index.z - (axis == 2 ? 1 : 0);
    double rhs = 0.0;
    for (int x = 0; x <= gx; ++x) {
      for (int y = 0; y <= gy; ++y) {
        for (int z = 0; z <= gz; ++z) {
          const int dx = gx - x + (axis == 0 ? 1 : 0);
          const int dy = gy - y + (axis == 1 ? 1 : 0);
          const int dz = gz - z + (axis == 2 ? 1 : 0);
          const int axis_degree = axis == 0 ? dx : (axis == 1 ? dy : dz);
          rhs += result.coeff[cube_index(x, y, z)] *
                 static_cast<double>(axis_degree) *
                 a.coeff[cube_index(dx, dy, dz)];
        }
      }
    }
    const int degree = axis == 0 ? index.x : (axis == 1 ? index.y : index.z);
    result.coeff[cube_index(index.x, index.y, index.z)] =
        rhs / static_cast<double>(degree);
  }
  return result;
}

struct KernelBase {
  Jet8 ux;
  Jet8 uy;
  Jet8 uz;
  Jet8 r;
  Jet8 inv_r;
  Jet8 inv_r2;
  Jet8 inv_r3;
  Jet8 screened;
};

KernelBase kernel_base(double tx, double ty, double tz, double sx, double sy,
                       double sz, double kappa) noexcept {
  KernelBase base;
  const Jet8 source_x = variable(sx, 0);
  const Jet8 source_y = variable(sy, 1);
  const Jet8 source_z = variable(sz, 2);
  base.ux = subtract(constant(tx), source_x);
  base.uy = subtract(constant(ty), source_y);
  base.uz = subtract(constant(tz), source_z);
  const Jet8 r2 = add(add(multiply(base.ux, base.ux),
                          multiply(base.uy, base.uy)),
                      multiply(base.uz, base.uz));
  if (r2.coeff[cube_index(0, 0, 0)] <= 1.0e-28)
    return base;
  base.r = square_root(r2);
  base.inv_r = reciprocal(base.r);
  base.inv_r2 = multiply(base.inv_r, base.inv_r);
  base.inv_r3 = multiply(base.inv_r2, base.inv_r);
  base.screened = exponential(scale(base.r, -std::max(0.0, kappa)));
  return base;
}

Jet8 green_jet(double tx, double ty, double tz, double sx, double sy,
               double sz, double kappa) noexcept {
  const KernelBase base = kernel_base(tx, ty, tz, sx, sy, sz, kappa);
  return scale(multiply(base.screened, base.inv_r), 1.0 / (4.0 * kPi));
}

Jet8 source_normal_component_jet(double tx, double ty, double tz, double sx,
                                 double sy, double sz, int component,
                                 double kappa) noexcept {
  const double safe_kappa = std::max(0.0, kappa);
  const KernelBase base =
      kernel_base(tx, ty, tz, sx, sy, sz, safe_kappa);
  const Jet8 one_plus_kr =
      add(constant(1.0), scale(base.r, safe_kappa));
  const Jet8 &u = component == 0 ? base.ux : (component == 1 ? base.uy : base.uz);
  return scale(multiply(multiply(base.screened, one_plus_kr),
                        multiply(u, base.inv_r3)),
               1.0 / (4.0 * kPi));
}

void translate_moments(const std::array<double, kCartesianJcpMomentMax> &moments,
                       const Jet8 &kernel,
                       std::array<double, kCartesianLocalTermMax> &local) {
  const std::vector<MultiIndex> &source =
      indices_through(kCartesianLocalSourceOrder);
  const std::vector<MultiIndex> &target =
      indices_through(kCartesianLocalTargetOrder);
  for (std::size_t b = 0; b < target.size(); ++b) {
    const MultiIndex &beta = target[b];
    const double sign = ((beta.x + beta.y + beta.z) % 2) == 0 ? 1.0 : -1.0;
    double value = 0.0;
    for (std::size_t a = 0; a < source.size(); ++a) {
      const MultiIndex &alpha = source[a];
      const int gx = alpha.x + beta.x;
      const int gy = alpha.y + beta.y;
      const int gz = alpha.z + beta.z;
      const double factor =
          sign * binomial(gx, alpha.x) * binomial(gy, alpha.y) *
          binomial(gz, alpha.z);
      value += moments[a] * factor * kernel.coeff[cube_index(gx, gy, gz)];
    }
    local[b] += value;
  }
}

double evaluate_scalar(const std::array<double, kCartesianLocalTermMax> &local,
                       double dx, double dy, double dz) noexcept {
  const std::vector<MultiIndex> &indices =
      indices_through(kCartesianLocalTargetOrder);
  double value = 0.0;
  for (std::size_t i = 0; i < indices.size(); ++i) {
    const MultiIndex &index = indices[i];
    value += local[i] * power(dx, index.x) * power(dy, index.y) *
             power(dz, index.z);
  }
  return value;
}

std::array<double, 3>
evaluate_gradient(const std::array<double, kCartesianLocalTermMax> &local,
                  double dx, double dy, double dz) noexcept {
  const std::vector<MultiIndex> &indices =
      indices_through(kCartesianLocalTargetOrder);
  std::array<double, 3> gradient{};
  for (std::size_t i = 0; i < indices.size(); ++i) {
    const MultiIndex &index = indices[i];
    if (index.x > 0)
      gradient[0] += static_cast<double>(index.x) * local[i] *
                     power(dx, index.x - 1) * power(dy, index.y) *
                     power(dz, index.z);
    if (index.y > 0)
      gradient[1] += static_cast<double>(index.y) * local[i] *
                     power(dx, index.x) * power(dy, index.y - 1) *
                     power(dz, index.z);
    if (index.z > 0)
      gradient[2] += static_cast<double>(index.z) * local[i] *
                     power(dx, index.x) * power(dy, index.y) *
                     power(dz, index.z - 1);
  }
  return gradient;
}

} // namespace

std::size_t cartesian_local_kernel_term_count() noexcept { return 165; }

std::size_t cartesian_local_term_count() noexcept { return 35; }

CartesianJcpLocalP4 translate_cartesian_jcp_m2l_p4(
    const CartesianJcpMoments &moments, double source_center_x,
    double source_center_y, double source_center_z, double target_center_x,
    double target_center_y, double target_center_z, double eps_ratio,
    double kappa) {
  CartesianJcpLocalP4 local;
  if (moments.order != kCartesianLocalSourceOrder ||
      moments.term_count != kCartesianLocalTermMax) {
    local.finite_outputs = false;
    return local;
  }
  const double safe_eps = std::max(1.0e-12, eps_ratio);
  const Jet8 green0 = green_jet(target_center_x, target_center_y,
                                target_center_z, source_center_x,
                                source_center_y, source_center_z, 0.0);
  const Jet8 greenk = green_jet(target_center_x, target_center_y,
                                target_center_z, source_center_x,
                                source_center_y, source_center_z, kappa);
  translate_moments(moments.phin, subtract(green0, greenk),
                    local.fields[static_cast<std::size_t>(JcpLocalField::K1)]);
  translate_moments(
      moments.phin, subtract(green0, scale(greenk, 1.0 / safe_eps)),
      local.fields[static_cast<std::size_t>(JcpLocalField::K3Potential)]);
  for (int component = 0; component < 3; ++component) {
    const Jet8 source0 = source_normal_component_jet(
        target_center_x, target_center_y, target_center_z, source_center_x,
        source_center_y, source_center_z, component, 0.0);
    const Jet8 sourcek = source_normal_component_jet(
        target_center_x, target_center_y, target_center_z, source_center_x,
        source_center_y, source_center_z, component, kappa);
    translate_moments(
        moments.phi_normal[static_cast<std::size_t>(component)],
        subtract(scale(sourcek, safe_eps), source0),
        local.fields[static_cast<std::size_t>(JcpLocalField::K2)]);
    translate_moments(
        moments.phi_normal[static_cast<std::size_t>(component)],
        subtract(source0, sourcek),
        local.fields[static_cast<std::size_t>(JcpLocalField::K4Potential)]);
  }
  for (const auto &field : local.fields) {
    for (const double value : field)
      local.finite_outputs = local.finite_outputs && std::isfinite(value);
  }
  return local;
}

void accumulate_cartesian_jcp_local(
    CartesianJcpLocalP4 &destination,
    const CartesianJcpLocalP4 &source) noexcept {
  for (std::size_t field = 0; field < destination.fields.size(); ++field) {
    for (std::size_t term = 0; term < destination.fields[field].size();
         ++term)
      destination.fields[field][term] += source.fields[field][term];
  }
  destination.finite_outputs =
      destination.finite_outputs && source.finite_outputs;
}

double evaluate_cartesian_jcp_local_scalar(
    const CartesianJcpLocalP4 &local, JcpLocalField field, double dx,
    double dy, double dz) noexcept {
  return evaluate_scalar(local.fields[static_cast<std::size_t>(field)], dx,
                         dy, dz);
}

std::array<double, 4> evaluate_cartesian_jcp_l2p(
    const PrototypePanel &target, double target_center_x,
    double target_center_y, double target_center_z,
    const CartesianJcpLocalP4 &local) noexcept {
  const double dx = target.x - target_center_x;
  const double dy = target.y - target_center_y;
  const double dz = target.z - target_center_z;
  std::array<double, 4> channels{};
  channels[0] = evaluate_scalar(
      local.fields[static_cast<std::size_t>(JcpLocalField::K1)], dx, dy, dz);
  channels[1] = evaluate_scalar(
      local.fields[static_cast<std::size_t>(JcpLocalField::K2)], dx, dy, dz);
  const auto k3_gradient = evaluate_gradient(
      local.fields[static_cast<std::size_t>(JcpLocalField::K3Potential)], dx,
      dy, dz);
  const auto k4_gradient = evaluate_gradient(
      local.fields[static_cast<std::size_t>(JcpLocalField::K4Potential)], dx,
      dy, dz);
  channels[2] = target.nx * k3_gradient[0] + target.ny * k3_gradient[1] +
                target.nz * k3_gradient[2];
  channels[3] = target.nx * k4_gradient[0] + target.ny * k4_gradient[1] +
                target.nz * k4_gradient[2];
  return channels;
}

} // namespace titan::tabi::fmm_m2l
