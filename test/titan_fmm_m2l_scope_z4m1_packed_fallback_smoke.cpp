#include "tabi/tabi_fmm_m2l_fallback.h"
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
  std::cerr << "FAIL titan_fmm_m2l_scope_z4m1_packed_fallback_smoke: "
            << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

Panel make_panel(std::size_t index) {
  Panel panel;
  const double t = static_cast<double>(index);
  panel.x = 1.4 * t + 0.1 * std::sin(t);
  panel.y = 0.3 * std::cos(0.7 * t);
  panel.z = 0.2 * std::sin(0.4 * t);
  panel.nx = 0.25;
  panel.ny = 0.35;
  panel.nz = std::sqrt(1.0 - panel.nx * panel.nx - panel.ny * panel.ny);
  panel.area = 0.16 + 0.002 * t;
  panel.diagonal = 1.0;
  panel.atom_index = index;
  panel.has_triangle_geometry = true;
  panel.v0x = panel.x - 0.12;
  panel.v0y = panel.y - 0.08;
  panel.v0z = panel.z;
  panel.v1x = panel.x + 0.12;
  panel.v1y = panel.y - 0.08;
  panel.v1z = panel.z;
  panel.v2x = panel.x;
  panel.v2y = panel.y + 0.16;
  panel.v2z = panel.z;
  panel.v0_index = 3 * index;
  panel.v1_index = 3 * index + 1;
  panel.v2_index = 3 * index + 2;
  return panel;
}

std::vector<Panel> make_panels() {
  std::vector<Panel> panels;
  for (std::size_t i = 0; i < 13; ++i)
    panels.push_back(make_panel(i));
  panels[2].v0_index = panels[1].v1_index;
  panels[5].v0_index = panels[4].v2_index;
  return panels;
}

std::vector<fmm::PackedFallbackLeafInput> make_leaf_work() {
  std::vector<fmm::PackedFallbackLeafInput> leaves(2);
  leaves[0].target_panels = {0, 1, 2, 3, 4};
  leaves[0].source_spans = {{0, 1, false}, {1, 2, true}, {3, 3, true}};
  leaves[1].target_panels = {5, 6, 7, 8, 9, 10, 11};
  leaves[1].source_spans = {{4, 3, false}, {7, 2, true}, {9, 1, true}};
  return leaves;
}

void make_density(std::size_t n, double phase, std::vector<double> &phi,
                  std::vector<double> &phin) {
  phi.resize(n);
  phin.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double t = static_cast<double>(i + 1);
    phi[i] = 0.4 + std::sin(0.31 * t + phase);
    phin[i] = -0.2 + 0.7 * std::cos(0.27 * t - phase);
  }
}

double rel_rms(const std::vector<double> &a, const std::vector<double> &b) {
  expect(a.size() == b.size(), "output sizes differ");
  double diff_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const double d = a[i] - b[i];
    diff_sq += d * d;
    ref_sq += b[i] * b[i];
  }
  return std::sqrt(diff_sq / std::max(1.0e-300, ref_sq));
}
} // namespace

int main() {
  const auto panels = make_panels();
  std::vector<std::size_t> source_order(panels.size());
  for (std::size_t i = 0; i < source_order.size(); ++i)
    source_order[i] = i;
  const auto leaves = make_leaf_work();
  auto plan = fmm::build_packed_jcp_fallback_plan(
      panels, source_order, leaves, 4, 32);
  const auto metrics = plan.metrics();
  expect(metrics.build_count == 1, "packed plan must build once");
  expect(metrics.worklist_hash != 0, "packed plan hash must be nonzero");
  expect(metrics.target_batch_size == 4 && metrics.source_tile_size == 32,
         "packed tile geometry mismatch");
  expect(metrics.task_count > 0 && metrics.metadata_bytes > 0,
         "packed worklist must contain tasks");
  expect(metrics.base_near_pairs + metrics.refined_fallback_pairs ==
             metrics.total_pairs,
         "base/refined pair accounting mismatch");
  expect(metrics.self_pairs == 10,
         "fixture must identify exactly ten covered self pairs");

  std::vector<double> phi_a, phin_a, phi_b, phin_b;
  make_density(panels.size(), 0.0, phi_a, phin_a);
  make_density(panels.size(), 0.9, phi_b, phin_b);
  TitanExecutor one(1);
  TitanExecutor four(4);
  const auto scalar = fmm::apply_packed_jcp_fallback(
      plan, panels, phi_a, phin_a, 1.0 / 80.0, 0.1, true,
      "analytic_singular", "off", "scalar", one);
  const auto packed = fmm::apply_packed_jcp_fallback(
      plan, panels, phi_a, phin_a, 1.0 / 80.0, 0.1, true,
      "analytic_singular", "off", "packed-avx2", four);
  expect(scalar.finite_outputs && packed.finite_outputs,
         "fallback outputs must be finite");
  expect(rel_rms(packed.phi, scalar.phi) <= 1.0e-12 &&
             rel_rms(packed.phin, scalar.phin) <= 1.0e-12,
         "packed fallback must match scalar fallback");
  expect(packed.scalar_self_pairs == metrics.self_pairs,
         "self terms must be evaluated exactly once");
  expect(packed.vectorized_lanes > 0 && packed.vectorized_sources > 0,
         "packed path must exercise AVX2 lanes");

  const auto packed_one = fmm::apply_packed_jcp_fallback(
      plan, panels, phi_a, phin_a, 1.0 / 80.0, 0.1, true,
      "analytic_singular", "off", "packed-avx2", one);
  expect(rel_rms(packed_one.phi, packed.phi) <= 1.0e-12 &&
             rel_rms(packed_one.phin, packed.phin) <= 1.0e-12,
         "packed fallback must be deterministic for nt=1/4");
  const auto changed = fmm::apply_packed_jcp_fallback(
      plan, panels, phi_b, phin_b, 1.0 / 80.0, 0.1, true,
      "analytic_singular", "off", "packed-avx2", four);
  expect(rel_rms(changed.phi, packed.phi) > 1.0e-6,
         "density refresh must change packed output");
  expect(plan.metrics().build_count == 1 &&
             plan.metrics().worklist_hash == metrics.worklist_hash,
         "density refresh must not rebuild packed work");

  std::cout << "titan_fmm_m2l_scope_z4m1_packed_fallback_smoke ok\n";
  return 0;
}
