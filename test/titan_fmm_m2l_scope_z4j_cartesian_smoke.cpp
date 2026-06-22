#include "tabi/tabi_fmm_m2l_cartesian.h"
#include "tabi/tabi_near.h"
#include "tabi/tabi_types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using titan::tabi::fmm_m2l::CartesianJcpMoments;
using titan::tabi::fmm_m2l::PrototypePanel;
using titan::tabi::fmm_m2l::build_cartesian_jcp_moments;
using titan::tabi::fmm_m2l::cartesian_moment_term_count;
using titan::tabi::fmm_m2l::evaluate_cartesian_jcp_m2p;
using TitanPBE::Tabi::Panel;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_fmm_m2l_scope_z4j_cartesian_smoke: " << message
            << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

Panel make_target() {
  Panel target;
  target.x = 8.0;
  target.y = 1.1;
  target.z = -0.7;
  const double norm = std::sqrt(1.0 + 0.3 * 0.3 + 0.2 * 0.2);
  target.nx = 1.0 / norm;
  target.ny = 0.3 / norm;
  target.nz = -0.2 / norm;
  target.area = 0.4;
  return target;
}

std::vector<Panel> make_sources(std::size_t count) {
  std::vector<Panel> panels(count);
  constexpr double golden = 2.3999632297286533222;
  for (std::size_t i = 0; i < count; ++i) {
    const double u = -1.0 + 2.0 * (static_cast<double>(i) + 0.5) /
                                  static_cast<double>(count);
    const double ring = std::sqrt(std::max(0.0, 1.0 - u * u));
    const double angle = golden * static_cast<double>(i);
    Panel panel;
    panel.x = 0.55 * ring * std::cos(angle);
    panel.y = 0.42 * ring * std::sin(angle);
    panel.z = 0.48 * u;
    const double nx = panel.x / (0.55 * 0.55);
    const double ny = panel.y / (0.42 * 0.42);
    const double nz = panel.z / (0.48 * 0.48);
    const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
    panel.nx = nx / norm;
    panel.ny = ny / norm;
    panel.nz = nz / norm;
    panel.area = 0.09 + 0.002 * static_cast<double>(i % 5);
    panels[i] = panel;
  }
  return panels;
}

PrototypePanel to_prototype(const Panel &panel) {
  PrototypePanel out;
  out.x = panel.x;
  out.y = panel.y;
  out.z = panel.z;
  out.nx = panel.nx;
  out.ny = panel.ny;
  out.nz = panel.nz;
  out.area = panel.area;
  return out;
}

std::array<double, 4> direct_channels(
    const Panel &target, const std::vector<Panel> &sources,
    const std::vector<double> &phi, const std::vector<double> &phin,
    double eps_ratio, double kappa) {
  std::array<double, 4> out{};
  for (std::size_t i = 0; i < sources.size(); ++i) {
    const auto kernel = TitanPBE::Tabi::evaluate_jcp_near_pair_kernel(
        target, sources.size(), sources[i], i, eps_ratio, kappa, true,
        "centroid", "off", nullptr);
    expect(kernel.valid, "direct oracle kernel must be valid");
    out[0] += kernel.k1 * phin[i];
    out[1] += kernel.k2 * phi[i];
    out[2] += kernel.k3 * phin[i];
    out[3] += kernel.k4 * phi[i];
  }
  return out;
}

double squared_error(const std::array<double, 4> &lhs,
                     const std::array<double, 4> &rhs) {
  double sum = 0.0;
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    const double delta = lhs[i] - rhs[i];
    sum += delta * delta;
  }
  return sum;
}

void test_term_counts() {
  expect(cartesian_moment_term_count(3) == 20,
         "order three must expose 20 Cartesian moments");
  expect(cartesian_moment_term_count(4) == 35,
         "order four must expose 35 Cartesian moments");
}

void test_p4_improves_all_channel_far_accuracy() {
  const Panel target = make_target();
  const std::vector<Panel> sources = make_sources(48);
  std::vector<PrototypePanel> prototype;
  std::vector<std::size_t> indices(sources.size());
  std::vector<double> phi(sources.size());
  std::vector<double> phin(sources.size());
  prototype.reserve(sources.size());
  std::iota(indices.begin(), indices.end(), std::size_t{0});
  for (std::size_t i = 0; i < sources.size(); ++i) {
    prototype.push_back(to_prototype(sources[i]));
    const double t = static_cast<double>(i + 1);
    phi[i] = 0.7 * std::sin(0.31 * t) + 0.2 * std::cos(0.11 * t);
    phin[i] = 0.6 * std::cos(0.23 * t) - 0.15 * std::sin(0.19 * t);
  }
  constexpr double eps_ratio = 1.0 / 80.0;
  constexpr double kappa = 0.12;
  const std::array<double, 4> direct =
      direct_channels(target, sources, phi, phin, eps_ratio, kappa);
  const CartesianJcpMoments p1 = build_cartesian_jcp_moments(
      prototype, phi, phin, indices, 0.0, 0.0, 0.0, 1);
  const CartesianJcpMoments p4 = build_cartesian_jcp_moments(
      prototype, phi, phin, indices, 0.0, 0.0, 0.0, 4);
  const std::array<double, 4> approx1 = evaluate_cartesian_jcp_m2p(
      to_prototype(target), 0.0, 0.0, 0.0, p1, eps_ratio, kappa);
  const std::array<double, 4> approx4 = evaluate_cartesian_jcp_m2p(
      to_prototype(target), 0.0, 0.0, 0.0, p4, eps_ratio, kappa);
  const double p1_error = squared_error(approx1, direct);
  const double p4_error = squared_error(approx4, direct);
  expect(std::isfinite(p4_error), "p4 channel error must be finite");
  expect(p4_error * 10.0 < p1_error,
         "p4 must reduce aggregate far-channel squared error by at least 10x");
  expect(std::abs(approx4[1] - direct[1]) <
             0.25 * std::max(1.0e-14, std::abs(approx1[1] - direct[1])),
         "p4 must materially reduce the dominant K2 error");
}

void test_zeroed_higher_moments_reduce_to_degree_one() {
  const Panel target = make_target();
  const std::vector<Panel> sources = make_sources(16);
  std::vector<PrototypePanel> prototype;
  std::vector<std::size_t> indices(sources.size());
  std::vector<double> phi(sources.size(), 0.35);
  std::vector<double> phin(sources.size(), -0.2);
  prototype.reserve(sources.size());
  std::iota(indices.begin(), indices.end(), std::size_t{0});
  for (const Panel &panel : sources)
    prototype.push_back(to_prototype(panel));
  CartesianJcpMoments p1 = build_cartesian_jcp_moments(
      prototype, phi, phin, indices, 0.0, 0.0, 0.0, 1);
  CartesianJcpMoments p4 = p1;
  p4.order = 4;
  p4.term_count = cartesian_moment_term_count(4);
  const auto a = evaluate_cartesian_jcp_m2p(
      to_prototype(target), 0.0, 0.0, 0.0, p1, 1.0 / 80.0, 0.12);
  const auto b = evaluate_cartesian_jcp_m2p(
      to_prototype(target), 0.0, 0.0, 0.0, p4, 1.0 / 80.0, 0.12);
  expect(squared_error(a, b) <= 1.0e-28,
         "zero higher moments must reproduce the degree-one evaluator");
}

} // namespace

int main() {
  test_term_counts();
  test_p4_improves_all_channel_far_accuracy();
  test_zeroed_higher_moments_reduce_to_degree_one();
  std::cout << "titan_fmm_m2l_scope_z4j_cartesian_smoke ok\n";
  return 0;
}
