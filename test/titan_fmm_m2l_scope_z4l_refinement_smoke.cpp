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
  std::cerr << "FAIL titan_fmm_m2l_scope_z4l_refinement_smoke: " << message
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
    const std::size_t cluster_index = i / 32;
    const double center = 2.5 * static_cast<double>(cluster_index);
    const double t = static_cast<double>(i % 32);
    panels[i].x = center + 0.42 * std::cos(0.41 * t);
    panels[i].y = 0.36 * std::sin(0.57 * t);
    panels[i].z = 0.28 * std::cos(0.33 * t);
    const double norm = std::sqrt(panels[i].x * panels[i].x +
                                  panels[i].y * panels[i].y +
                                  panels[i].z * panels[i].z);
    panels[i].nx = norm > 1.0e-12 ? panels[i].x / norm : 1.0;
    panels[i].ny = norm > 1.0e-12 ? panels[i].y / norm : 0.0;
    panels[i].nz = norm > 1.0e-12 ? panels[i].z / norm : 0.0;
    panels[i].area = 0.16 + 0.004 * static_cast<double>(i % 7);
    panels[i].diagonal = 1.0;
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

fmm::PrototypeConfig config() {
  fmm::PrototypeConfig cfg;
  cfg.kernel = fmm::KernelKind::Yukawa;
  cfg.kappa = 0.1;
  cfg.max_leaf_size = 8;
  cfg.expansion_order = 4;
  cfg.separation_theta = 0.8;
  cfg.threads = 4;
  return cfg;
}

void expect_same_records(const fmm::JcpOperatorAttributionResult &lhs,
                         const fmm::JcpOperatorAttributionResult &rhs) {
  expect(lhs.records.size() == rhs.records.size(),
         "record count changed unexpectedly");
  for (std::size_t i = 0; i < lhs.records.size(); ++i) {
    expect(lhs.records[i].target_panel == rhs.records[i].target_panel,
           "target ordering changed");
    for (std::size_t channel = 0; channel < 4; ++channel) {
      expect(std::abs(lhs.records[i].approximate[channel] -
                      rhs.records[i].approximate[channel]) <= 1.0e-13,
             "thread count changed refined attribution values");
    }
  }
}

} // namespace

int main() {
  const std::vector<Panel> panels = make_panels();
  const fmm::PrototypeConfig cfg = config();
  fmm::FmmM2lAnalyticContext context(cfg, to_prototype(panels), 1.0 / 80.0,
                                     cfg.kappa);
  std::vector<double> phi;
  std::vector<double> phin;
  make_density(panels.size(), phi, phin);
  TitanExecutor one(1);
  TitanExecutor four(4);

  const auto baseline = context.attribute_local_p4(
      phi, phin, 24, "solution", one, true, "centroid", "off", &panels);
  const auto explicit_off = context.attribute_local_p4(
      phi, phin, 24, "solution", one, true, "centroid", "off", &panels,
      "off", 4.0e-3);
  expect_same_records(baseline, explicit_off);
  expect(explicit_off.local_refine_mode == "off",
         "explicit refinement-off telemetry mismatch");
  expect(explicit_off.local_refined_cells == 0,
         "refinement-off unexpectedly changed work");

  const auto loose = context.attribute_local_p4(
      phi, phin, 24, "solution", one, true, "centroid", "off", &panels,
      "p3-p4-tail", 8.0e-3);
  const auto tight = context.attribute_local_p4(
      phi, phin, 24, "solution", one, true, "centroid", "off", &panels,
      "p3-p4-tail", 1.0e-4);
  expect(tight.finite_outputs, "tight refinement must remain finite");
  expect(tight.local_refined_cells > 0,
         "tight tolerance must refine the risky separated fixture");
  expect(tight.local_direct_fallback_pairs > 0,
         "tight refinement must count direct fallback panel pairs");
  expect(tight.local_far_pairs_after != tight.local_far_pairs_before ||
             tight.local_direct_fallback_panels > 0,
         "tight refinement must change far work or use direct fallback");
  expect(tight.local_tail_ratio_after <= tight.local_tail_ratio_before,
         "refinement must not increase the estimated tail ratio");
  expect(tight.local_refined_cells >= loose.local_refined_cells,
         "tighter tolerance must not reduce refinement work");
  expect(tight.recomposition_rel_error <= 1.0e-10,
         "refined local/direct work must recompose the candidate operator");

  const auto tight_parallel = context.attribute_local_p4(
      phi, phin, 24, "solution", four, true, "centroid", "off", &panels,
      "p3-p4-tail", 1.0e-4);
  expect_same_records(tight, tight_parallel);
  expect(tight.local_refined_cells == tight_parallel.local_refined_cells &&
             tight.local_split_children ==
                 tight_parallel.local_split_children &&
             tight.local_direct_fallback_panels ==
                 tight_parallel.local_direct_fallback_panels,
         "refinement telemetry must be thread deterministic");

  std::cout << "titan_fmm_m2l_scope_z4l_refinement_smoke ok\n";
  return 0;
}
