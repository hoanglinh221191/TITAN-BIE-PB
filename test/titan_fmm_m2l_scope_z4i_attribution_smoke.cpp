#include "tabi/tabi_fmm_m2l_prototype.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

using titan::tabi::fmm_m2l::FmmM2lAnalyticContext;
using titan::tabi::fmm_m2l::KernelKind;
using titan::tabi::fmm_m2l::PrototypeConfig;
using titan::tabi::fmm_m2l::PrototypePanel;
using TitanPBE::Tabi::Panel;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_fmm_m2l_scope_z4i_attribution_smoke: "
            << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

std::vector<Panel> make_panels(std::size_t count, bool separated) {
  std::vector<Panel> panels(count);
  for (std::size_t i = 0; i < count; ++i) {
    const double cluster = separated && i >= count / 2 ? 8.0 : 0.0;
    const double t = static_cast<double>(i % (count / 2));
    panels[i].x = cluster + 0.12 * t;
    panels[i].y = 0.21 * std::sin(0.47 * t);
    panels[i].z = 0.18 * std::cos(0.31 * t);
    const double norm = std::sqrt(panels[i].x * panels[i].x +
                                  panels[i].y * panels[i].y +
                                  panels[i].z * panels[i].z);
    panels[i].nx = norm > 1.0e-12 ? panels[i].x / norm : 1.0;
    panels[i].ny = norm > 1.0e-12 ? panels[i].y / norm : 0.0;
    panels[i].nz = norm > 1.0e-12 ? panels[i].z / norm : 0.0;
    panels[i].area = 0.35 + 0.01 * static_cast<double>(i % 5);
    panels[i].diagonal = 1.0;
  }
  return panels;
}

std::vector<PrototypePanel>
to_prototype(const std::vector<Panel> &panels) {
  std::vector<PrototypePanel> converted;
  converted.reserve(panels.size());
  for (const Panel &panel : panels) {
    PrototypePanel out;
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
    phi[i] = 0.4 * std::sin(0.13 * t) + 0.1;
    phin[i] = 0.3 * std::cos(0.17 * t) - 0.05;
  }
}

PrototypeConfig config(double theta) {
  PrototypeConfig result;
  result.kernel = KernelKind::Yukawa;
  result.kappa = 0.1;
  result.max_leaf_size = 4;
  result.expansion_order = 4;
  result.separation_theta = theta;
  result.threads = 4;
  return result;
}

void test_all_near_recomposes_exactly() {
  const std::vector<Panel> panels = make_panels(24, false);
  const PrototypeConfig cfg = config(0.0);
  FmmM2lAnalyticContext context(cfg, to_prototype(panels), 1.0 / 80.0,
                                cfg.kappa);
  std::vector<double> phi;
  std::vector<double> phin;
  make_density(panels.size(), phi, phin);
  TitanExecutor executor(4);
  const auto result = context.attribute(phi, phin, 8, "solution", executor,
                                        true, "centroid", "off",
                                        &panels);
  expect(result.finite_outputs, "all-near attribution must be finite");
  expect(result.sampled_targets.size() == 8,
         "all-near attribution sample count mismatch");
  expect(result.far_records == 0, "theta zero should not emit far records");
  expect(result.near_records == result.sampled_targets.size(),
         "all-near attribution should aggregate one near record per target");
  expect(result.recomposition_rel_error <= 1.0e-12,
         "all-near attribution must recompose the current operator");
  for (const auto &record : result.records) {
    if (record.region != "near")
      continue;
    for (std::size_t channel = 0; channel < 4; ++channel)
      expect(std::abs(record.approximate[channel] - record.direct[channel]) <=
                 1.0e-13,
             "shared near kernel must have zero attribution error");
  }
}

void test_far_attribution_reports_current_moment_degree() {
  const std::vector<Panel> panels = make_panels(48, true);
  const PrototypeConfig cfg = config(0.6);
  FmmM2lAnalyticContext context(cfg, to_prototype(panels), 1.0 / 80.0,
                                cfg.kappa);
  std::vector<double> phi;
  std::vector<double> phin;
  make_density(panels.size(), phi, phin);
  TitanExecutor executor(1);
  const auto result = context.attribute(phi, phin, 12, "phi_only", executor,
                                        true, "centroid", "off",
                                        &panels);
  expect(result.finite_outputs, "far attribution must be finite");
  expect(result.configured_expansion_order == 4,
         "configured expansion order telemetry mismatch");
  expect(result.implemented_far_moment_degree == 1,
         "current far operator must report monopole/dipole degree one");
  expect(result.far_records > 0,
         "separated clusters should produce far attribution records");
  expect(result.records.size() == result.near_records + result.far_records,
         "attribution record counters must recompose record storage");
  for (const auto &record : result.records) {
    expect(record.target_panel < panels.size(),
           "attribution target panel out of range");
    if (record.region == "far") {
      expect(record.source_depth >= 0, "far record must report tree depth");
      expect(record.source_panels > 0,
             "far record must report source panel count");
      expect(std::isfinite(record.mac_ratio),
             "far record must report finite MAC ratio");
      expect(record.mac_ratio <= cfg.separation_theta + 1.0e-12,
             "far record MAC ratio must match the accepted interaction metric");
      expect(std::isfinite(record.expansion_mac_ratio),
             "far record must report finite expansion-center MAC ratio");
      expect(record.source_expansion_radius + record.source_center_offset +
                     1.0e-12 >=
                 record.source_radius,
             "expansion radius plus center offset must enclose source panels");
    }
  }
}

void test_target_sampling_is_thread_deterministic() {
  const std::vector<Panel> panels = make_panels(40, true);
  const PrototypeConfig cfg = config(0.45);
  FmmM2lAnalyticContext context(cfg, to_prototype(panels), 1.0 / 80.0,
                                cfg.kappa);
  std::vector<double> phi;
  std::vector<double> phin;
  make_density(panels.size(), phi, phin);
  TitanExecutor one(1);
  TitanExecutor four(4);
  const auto serial = context.attribute(phi, phin, 10, "solution", one, true,
                                        "centroid", "off", &panels);
  const auto parallel = context.attribute(phi, phin, 10, "solution", four,
                                          true, "centroid", "off",
                                          &panels);
  expect(serial.sampled_targets == parallel.sampled_targets,
         "target sampling must be thread deterministic");
  expect(serial.records.size() == parallel.records.size(),
         "thread count must not change attribution record count");
  expect(std::abs(serial.recomposition_rel_error -
                  parallel.recomposition_rel_error) <= 1.0e-15,
         "thread count must not change attribution recomposition");
}

} // namespace

int main() {
  test_all_near_recomposes_exactly();
  test_far_attribution_reports_current_moment_degree();
  test_target_sampling_is_thread_deterministic();
  std::cout << "titan_fmm_m2l_scope_z4i_attribution_smoke ok\n";
  return 0;
}
