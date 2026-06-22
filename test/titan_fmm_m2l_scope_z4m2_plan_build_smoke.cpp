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
  std::cerr << "FAIL titan_fmm_m2l_scope_z4m2_plan_build_smoke: "
            << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

std::vector<Panel> make_panels() {
  std::vector<Panel> panels(128);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const double cluster = 4.0 * static_cast<double>(i / 32);
    const double t = static_cast<double>(i % 32);
    panels[i].x = cluster + 0.38 * std::cos(0.31 * t);
    panels[i].y = 0.34 * std::sin(0.47 * t);
    panels[i].z = 0.29 * std::cos(0.23 * t);
    const double norm = std::sqrt(panels[i].x * panels[i].x +
                                  panels[i].y * panels[i].y +
                                  panels[i].z * panels[i].z);
    panels[i].nx = norm > 1.0e-12 ? panels[i].x / norm : 1.0;
    panels[i].ny = norm > 1.0e-12 ? panels[i].y / norm : 0.0;
    panels[i].nz = norm > 1.0e-12 ? panels[i].z / norm : 0.0;
    panels[i].area = 0.12 + 0.002 * static_cast<double>(i % 11);
    panels[i].diagonal = 1.0;
  }
  return panels;
}

std::vector<fmm::PrototypePanel>
to_prototype(const std::vector<Panel> &panels) {
  std::vector<fmm::PrototypePanel> result;
  result.reserve(panels.size());
  for (const Panel &panel : panels) {
    fmm::PrototypePanel item;
    item.x = panel.x;
    item.y = panel.y;
    item.z = panel.z;
    item.nx = panel.nx;
    item.ny = panel.ny;
    item.nz = panel.nz;
    item.area = panel.area;
    result.push_back(item);
  }
  return result;
}

fmm::PrototypeConfig make_config(int threads, double theta = 0.6) {
  fmm::PrototypeConfig config;
  config.kernel = fmm::KernelKind::Yukawa;
  config.kappa = 0.1;
  config.max_leaf_size = 8;
  config.expansion_order = 4;
  config.separation_theta = theta;
  config.threads = threads;
  return config;
}

void make_density(std::size_t count, double phase, std::vector<double> &phi,
                  std::vector<double> &phin) {
  phi.resize(count);
  phin.resize(count);
  for (std::size_t i = 0; i < count; ++i) {
    const double t = static_cast<double>(i + 1);
    phi[i] = std::sin(0.19 * t + phase);
    phin[i] = 0.65 * std::cos(0.27 * t - phase);
  }
}

double relative_rms(const std::vector<double> &lhs,
                    const std::vector<double> &rhs) {
  expect(lhs.size() == rhs.size(), "vector sizes differ");
  double difference_sq = 0.0;
  double reference_sq = 0.0;
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    const double difference = lhs[i] - rhs[i];
    difference_sq += difference * difference;
    reference_sq += rhs[i] * rhs[i];
  }
  return std::sqrt(difference_sq / std::max(1.0e-300, reference_sq));
}
} // namespace

int main() {
  const std::vector<Panel> panels = make_panels();
  const auto prototype = to_prototype(panels);
  std::vector<double> phi_a, phin_a, phi_b, phin_b;
  make_density(panels.size(), 0.0, phi_a, phin_a);
  make_density(panels.size(), 0.7, phi_b, phin_b);
  TitanExecutor one(1);
  TitanExecutor four(4);
  constexpr double stress_tolerance = 1.0e-4;

  fmm::FmmM2lAnalyticContext legacy_context(
      make_config(1), prototype, 1.0 / 80.0, 0.1);
  const auto legacy = legacy_context.apply_local_p4_candidate(
      phi_a, phin_a, one, true, "centroid", "off", &panels,
      "operator-tail-static", stress_tolerance, "kernel-jet", 256, "scalar",
      "legacy-rescan");

  fmm::FmmM2lAnalyticContext metric_context(
      make_config(4), prototype, 1.0 / 80.0, 0.1);
  const auto metric_a = metric_context.apply_local_p4_candidate(
      phi_a, phin_a, four, true, "centroid", "off", &panels,
      "operator-tail-static", stress_tolerance, "kernel-jet", 256, "scalar",
      "metric-cache");
  const auto metric_b = metric_context.apply_local_p4_candidate(
      phi_b, phin_b, four, true, "centroid", "off", &panels,
      "operator-tail-static", stress_tolerance, "kernel-jet", 256, "scalar",
      "metric-cache");

  expect(metric_a.plan_builder == "metric-cache",
         "metric-cache builder must be reported");
  expect(legacy.plan_hash == metric_a.plan_hash,
         "metric-cache must preserve the legacy plan hash");
  expect(legacy.refined_cells == metric_a.refined_cells &&
             legacy.split_children == metric_a.split_children &&
             legacy.direct_fallback_pairs == metric_a.direct_fallback_pairs &&
             legacy.far_pairs == metric_a.far_pairs,
         "metric-cache must preserve refinement counters");
  expect(relative_rms(metric_a.phi, legacy.phi) <= 1.0e-12 &&
             relative_rms(metric_a.phin, legacy.phin) <= 1.0e-12,
         "metric-cache must preserve candidate output");
  expect(metric_a.plan_pair_metric_evaluations > 0 &&
             metric_a.plan_pair_metric_cache_hits > 0,
         "metric-cache must evaluate and reuse pair metrics");
  expect(metric_a.plan_refinement_iterations > 0,
         "fixture must exercise iterative refinement");
  expect(metric_a.plan_metric_cache_peak_bytes > 0,
         "metric-cache memory telemetry must be nonzero");
  expect(metric_a.plan_parallel_target_leaves > 0,
         "target-leaf parallel telemetry must be nonzero");
  expect(metric_b.plan_hash == metric_a.plan_hash &&
             metric_b.plan_build_count == 1,
         "density changes must not rebuild the metric-cache plan");

  fmm::FmmM2lAnalyticContext metric_one_context(
      make_config(1), prototype, 1.0 / 80.0, 0.1);
  const auto metric_one = metric_one_context.apply_local_p4_candidate(
      phi_a, phin_a, one, true, "centroid", "off", &panels,
      "operator-tail-static", stress_tolerance, "kernel-jet", 256, "scalar",
      "metric-cache");
  expect(metric_one.plan_hash == metric_a.plan_hash &&
             relative_rms(metric_one.phi, metric_a.phi) <= 1.0e-12 &&
             relative_rms(metric_one.phin, metric_a.phin) <= 1.0e-12,
         "metric-cache plan and output must be deterministic for nt=1/4");

  fmm::FmmM2lAnalyticContext all_near_context(
      make_config(1, 0.0), prototype, 1.0 / 80.0, 0.1);
  const auto all_near = all_near_context.apply_local_p4_candidate(
      phi_a, phin_a, one, true, "centroid", "off", &panels,
      "operator-tail-static", 2.0e-3, "kernel-jet", 256, "scalar",
      "metric-cache");
  expect(all_near.far_pairs == 0 &&
             all_near.plan_pair_metric_evaluations == 0 &&
             all_near.plan_metric_cache_peak_bytes == 0,
         "all-near work must not create pair metrics");

  std::cout << "titan_fmm_m2l_scope_z4m2_plan_build_smoke ok\n";
  return 0;
}
