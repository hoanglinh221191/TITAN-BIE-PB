#include "tabi/tabi_fmm_m2l_cartesian.h"
#include "tabi/tabi_fmm_m2l_local.h"

#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace fmm = titan::tabi::fmm_m2l;

double relative_rms(const std::vector<double> &lhs,
                    const std::vector<double> &rhs) {
  double diff_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    const double delta = lhs[i] - rhs[i];
    diff_sq += delta * delta;
    ref_sq += rhs[i] * rhs[i];
  }
  return std::sqrt(diff_sq / std::max(1.0e-300, ref_sq));
}

int main() {
  if (fmm::cartesian_local_kernel_term_count() != 165)
    throw std::runtime_error("order-8 kernel jet must contain 165 terms");
  if (fmm::cartesian_local_term_count() != 35)
    throw std::runtime_error("p4 local expansion must contain 35 terms");

  std::vector<fmm::PrototypePanel> sources(12);
  std::vector<double> phi(sources.size(), 0.0);
  std::vector<double> phin(sources.size(), 0.0);
  std::vector<std::size_t> source_ids(sources.size(), 0);
  for (std::size_t i = 0; i < sources.size(); ++i) {
    const double t = 0.47 * static_cast<double>(i + 1);
    sources[i].x = 0.18 * std::cos(t);
    sources[i].y = 0.15 * std::sin(t);
    sources[i].z = 0.11 * std::cos(1.7 * t);
    const double norm = std::sqrt(sources[i].x * sources[i].x +
                                  sources[i].y * sources[i].y +
                                  sources[i].z * sources[i].z);
    sources[i].nx = sources[i].x / norm;
    sources[i].ny = sources[i].y / norm;
    sources[i].nz = sources[i].z / norm;
    sources[i].area = 0.08 + 0.003 * static_cast<double>(i);
    phi[i] = 0.3 * std::sin(0.3 * static_cast<double>(i + 1)) + 0.2;
    phin[i] = 0.2 * std::cos(0.2 * static_cast<double>(i + 1)) - 0.1;
    source_ids[i] = i;
  }

  const auto moments = fmm::build_cartesian_jcp_moments(
      sources, phi, phin, std::span<const std::size_t>(source_ids), 0.0, 0.0,
      0.0, 4);
  constexpr double tx = 5.5;
  constexpr double ty = 0.35;
  constexpr double tz = -0.2;
  constexpr double eps_ratio = 80.0;
  constexpr double kappa = 0.125;
  const auto local = fmm::translate_cartesian_jcp_m2l_p4(
      moments, 0.0, 0.0, 0.0, tx, ty, tz, eps_ratio, kappa);

  fmm::PrototypePanel center;
  center.x = tx;
  center.y = ty;
  center.z = tz;
  center.nx = 0.3;
  center.ny = -0.4;
  center.nz = std::sqrt(0.75);
  const auto center_m2p = fmm::evaluate_cartesian_jcp_m2p(
      center, 0.0, 0.0, 0.0, moments, eps_ratio, kappa);
  const auto center_l2p = fmm::evaluate_cartesian_jcp_l2p(
      center, tx, ty, tz, local);
  for (std::size_t channel = 0; channel < 4; ++channel) {
    if (std::abs(center_m2p[channel] - center_l2p[channel]) >
        1.0e-11 * std::max(1.0, std::abs(center_m2p[channel]))) {
      std::cerr << "center channel " << channel << " m2p="
                << center_m2p[channel] << " l2p=" << center_l2p[channel]
                << '\n';
      throw std::runtime_error("M2L center value does not match p4 M2P");
    }
  }

  std::vector<double> local_values;
  std::vector<double> m2p_values;
  for (int i = 0; i < 9; ++i) {
    fmm::PrototypePanel target = center;
    target.x += 0.025 * static_cast<double>((i % 3) - 1);
    target.y += 0.02 * static_cast<double>(((i / 3) % 3) - 1);
    target.z += 0.015 * static_cast<double>((i % 2) ? 1 : -1);
    const auto m2p = fmm::evaluate_cartesian_jcp_m2p(
        target, 0.0, 0.0, 0.0, moments, eps_ratio, kappa);
    const auto l2p =
        fmm::evaluate_cartesian_jcp_l2p(target, tx, ty, tz, local);
    for (std::size_t channel = 0; channel < 4; ++channel) {
      m2p_values.push_back(m2p[channel]);
      local_values.push_back(l2p[channel]);
    }
  }
  if (relative_rms(local_values, m2p_values) > 1.0e-7)
    throw std::runtime_error("p4 local expansion does not reproduce p4 M2P");

  const double h = 1.0e-5;
  const auto field = fmm::JcpLocalField::K3Potential;
  const double plus = fmm::evaluate_cartesian_jcp_local_scalar(
      local, field, h, 0.0, 0.0);
  const double minus = fmm::evaluate_cartesian_jcp_local_scalar(
      local, field, -h, 0.0, 0.0);
  fmm::PrototypePanel x_normal = center;
  x_normal.nx = 1.0;
  x_normal.ny = 0.0;
  x_normal.nz = 0.0;
  const auto derivative_channels = fmm::evaluate_cartesian_jcp_l2p(
      x_normal, tx, ty, tz, local);
  const double finite_difference = (plus - minus) / (2.0 * h);
  if (std::abs(derivative_channels[2] - finite_difference) > 1.0e-8)
    throw std::runtime_error("K3 target-normal gradient is inconsistent");

  std::cout << "titan_fmm_m2l_scope_z4k_local_smoke ok\n";
  return 0;
}
