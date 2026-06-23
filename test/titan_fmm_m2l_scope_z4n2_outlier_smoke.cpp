#include "tabi/tabi_fmm_m2l_prototype.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {
namespace fmm = titan::tabi::fmm_m2l;
using TitanPBE::Tabi::Panel;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_fmm_m2l_scope_z4n2_outlier_smoke: " << message
            << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

std::vector<Panel> make_panels() {
  std::vector<Panel> panels(96);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const std::size_t cluster = i / 32;
    const double center = 2.2 * static_cast<double>(cluster);
    const double t = static_cast<double>(i % 32);
    panels[i].x = center + 0.44 * std::cos(0.41 * t);
    panels[i].y = 0.37 * std::sin(0.53 * t);
    panels[i].z = 0.31 * std::cos(0.29 * t);
    const double norm = std::sqrt(panels[i].x * panels[i].x +
                                  panels[i].y * panels[i].y +
                                  panels[i].z * panels[i].z);
    panels[i].nx = norm > 1.0e-12 ? panels[i].x / norm : 1.0;
    panels[i].ny = norm > 1.0e-12 ? panels[i].y / norm : 0.0;
    panels[i].nz = norm > 1.0e-12 ? panels[i].z / norm : 0.0;
    panels[i].area = 0.13 + 0.003 * static_cast<double>(i % 9);
    panels[i].diagonal = 1.0;
    const double ax = std::abs(panels[i].nx) < 0.8 ? 1.0 : 0.0;
    const double ay = std::abs(panels[i].nx) < 0.8 ? 0.0 : 1.0;
    const double az = 0.0;
    double ux = panels[i].ny * az - panels[i].nz * ay;
    double uy = panels[i].nz * ax - panels[i].nx * az;
    double uz = panels[i].nx * ay - panels[i].ny * ax;
    const double u_norm = std::sqrt(ux * ux + uy * uy + uz * uz);
    ux /= u_norm;
    uy /= u_norm;
    uz /= u_norm;
    const double vx = panels[i].ny * uz - panels[i].nz * uy;
    const double vy = panels[i].nz * ux - panels[i].nx * uz;
    const double vz = panels[i].nx * uy - panels[i].ny * ux;
    const double radius = std::sqrt(4.0 * panels[i].area /
                                    (3.0 * std::sqrt(3.0)));
    panels[i].v0x = panels[i].x + radius * ux;
    panels[i].v0y = panels[i].y + radius * uy;
    panels[i].v0z = panels[i].z + radius * uz;
    panels[i].v1x = panels[i].x +
                    radius * (-0.5 * ux + 0.8660254037844386 * vx);
    panels[i].v1y = panels[i].y +
                    radius * (-0.5 * uy + 0.8660254037844386 * vy);
    panels[i].v1z = panels[i].z +
                    radius * (-0.5 * uz + 0.8660254037844386 * vz);
    panels[i].v2x = panels[i].x +
                    radius * (-0.5 * ux - 0.8660254037844386 * vx);
    panels[i].v2y = panels[i].y +
                    radius * (-0.5 * uy - 0.8660254037844386 * vy);
    panels[i].v2z = panels[i].z +
                    radius * (-0.5 * uz - 0.8660254037844386 * vz);
    panels[i].has_triangle_geometry = true;
  }
  return panels;
}

std::vector<fmm::PrototypePanel>
to_prototype(const std::vector<Panel> &panels) {
  std::vector<fmm::PrototypePanel> converted;
  converted.reserve(panels.size());
  for (const Panel &panel : panels) {
    fmm::PrototypePanel out;
    out.x = panel.x;
    out.y = panel.y;
    out.z = panel.z;
    out.nx = panel.nx;
    out.ny = panel.ny;
    out.nz = panel.nz;
    out.area = panel.area;
    converted.push_back(out);
  }
  return converted;
}

void make_density(std::size_t count, std::vector<double> &phi,
                  std::vector<double> &phin) {
  phi.resize(count);
  phin.resize(count);
  for (std::size_t i = 0; i < count; ++i) {
    const double t = static_cast<double>(i + 1);
    phi[i] = 0.35 * std::sin(0.19 * t) + 0.11;
    phin[i] = 0.27 * std::cos(0.23 * t) - 0.07;
  }
}

fmm::PrototypeConfig make_config(int threads) {
  fmm::PrototypeConfig cfg;
  cfg.kernel = fmm::KernelKind::Yukawa;
  cfg.kappa = 0.1;
  cfg.max_leaf_size = 8;
  cfg.expansion_order = 4;
  cfg.separation_theta = 0.6;
  cfg.threads = threads;
  return cfg;
}

double rel_rms(const std::vector<double> &lhs, const std::vector<double> &rhs) {
  expect(lhs.size() == rhs.size(), "vector size mismatch");
  double diff_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    const double diff = lhs[i] - rhs[i];
    diff_sq += diff * diff;
    ref_sq += rhs[i] * rhs[i];
  }
  return std::sqrt(diff_sq / std::max(1.0e-300, ref_sq));
}

} // namespace

int main() {
  const std::vector<Panel> panels = make_panels();
  const auto prototype = to_prototype(panels);
  std::vector<double> phi, phin;
  make_density(panels.size(), phi, phin);
  TitanExecutor one(1);
  TitanExecutor four(4);

  fmm::FmmM2lAnalyticContext baseline_context(
      make_config(1), prototype, 1.0 / 80.0, 0.1);
  const auto baseline = baseline_context.apply_local_p4_candidate(
      phi, phin, one, true, "analytic_singular", "off", &panels,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "packed-avx2", "metric-cache");

  fmm::FmmM2lAnalyticContext linf_context(make_config(4), prototype,
                                          1.0 / 80.0, 0.1);
  const auto linf = linf_context.apply_local_p4_candidate(
      phi, phin, four, true, "analytic_singular", "off", &panels,
      "operator-tail-linf", 2.0e-3, "green-jet-compressed", 256,
      "packed-avx2", "metric-cache", 4, "rademacher", "off", 0.0);
  expect(linf.refine_mode == "operator-tail-linf",
         "Linf refinement mode must be reported");
  expect(linf.linf_outlier_panels > 0,
         "Linf refinement must record panel outliers");
  expect(linf.linf_max_relative_tail > 0.0,
         "Linf refinement must report max relative tail");
  expect(linf.direct_fallback_pairs >= baseline.direct_fallback_pairs,
         "Linf refinement must not reduce fallback work on the stress fixture");
  expect(linf.recomposition_rel_error <= 1.0e-10,
         "Linf refinement must recompose the candidate operator");

  fmm::FmmM2lAnalyticContext geometry_context(make_config(4), prototype,
                                              1.0 / 80.0, 0.1);
  const auto geometry = geometry_context.apply_local_p4_candidate(
      phi, phin, four, true, "analytic_singular", "off", &panels,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "packed-avx2", "metric-cache", 4, "rademacher",
      "strict-geometry-mac", 0.35);
  expect(geometry.geometry_guard_rejected_cells > 0,
         "strict geometry guard must reject at least one source cell");
  expect(geometry.geometry_guard_direct_fallback_pairs > 0,
         "strict geometry guard must report fallback pairs");
  expect(geometry.geometry_ratio_max > 0.0,
         "strict geometry guard must report a max ratio");

  fmm::FmmM2lAnalyticContext combined_context(make_config(4), prototype,
                                              1.0 / 80.0, 0.1);
  const auto combined = combined_context.apply_local_p4_candidate(
      phi, phin, four, true, "analytic_singular", "off", &panels,
      "operator-tail-linf", 2.0e-3, "green-jet-compressed", 256,
      "packed-avx2", "metric-cache", 4, "rademacher",
      "strict-geometry-mac", 0.35);
  expect(combined.direct_fallback_pairs >= linf.direct_fallback_pairs,
         "combined guard must include Linf fallback work");
  expect(combined.direct_fallback_pairs >= geometry.direct_fallback_pairs,
         "combined guard must include geometry fallback work");

  fmm::FmmM2lAnalyticContext combined_one_context(make_config(1), prototype,
                                                  1.0 / 80.0, 0.1);
  const auto combined_one = combined_one_context.apply_local_p4_candidate(
      phi, phin, one, true, "analytic_singular", "off", &panels,
      "operator-tail-linf", 2.0e-3, "green-jet-compressed", 256,
      "packed-avx2", "metric-cache", 4, "rademacher",
      "strict-geometry-mac", 0.35);
  expect(combined.plan_hash == combined_one.plan_hash,
         "combined outlier plan hash must be deterministic for nt=1/4");
  expect(rel_rms(combined.phi, combined_one.phi) <= 1.0e-12 &&
             rel_rms(combined.phin, combined_one.phin) <= 1.0e-12,
         "combined outlier output must be deterministic for nt=1/4");

  std::cout << "titan_fmm_m2l_scope_z4n2_outlier_smoke ok\n";
  return 0;
}
