#include "tabi/tabi_fmm_m2l_local.h"
#include "tabi/tabi_fmm_m2l_prototype.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace fmm = titan::tabi::fmm_m2l;
using TitanPBE::Tabi::Panel;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_fmm_m2l_scope_z4m_static_cache_smoke: "
            << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

std::vector<Panel> make_panels() {
  std::vector<Panel> panels(96);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const double cluster = 3.0 * static_cast<double>(i / 32);
    const double t = static_cast<double>(i % 32);
    panels[i].x = cluster + 0.4 * std::cos(0.37 * t);
    panels[i].y = 0.35 * std::sin(0.51 * t);
    panels[i].z = 0.3 * std::cos(0.29 * t);
    const double r = std::sqrt(panels[i].x * panels[i].x +
                               panels[i].y * panels[i].y +
                               panels[i].z * panels[i].z);
    panels[i].nx = r > 1.0e-12 ? panels[i].x / r : 1.0;
    panels[i].ny = r > 1.0e-12 ? panels[i].y / r : 0.0;
    panels[i].nz = r > 1.0e-12 ? panels[i].z / r : 0.0;
    panels[i].area = 0.14 + 0.003 * static_cast<double>(i % 9);
    panels[i].diagonal = 1.0;
  }
  return panels;
}

std::vector<fmm::PrototypePanel>
to_prototype(const std::vector<Panel> &panels) {
  std::vector<fmm::PrototypePanel> out;
  out.reserve(panels.size());
  for (const Panel &panel : panels) {
    fmm::PrototypePanel item;
    item.x = panel.x;
    item.y = panel.y;
    item.z = panel.z;
    item.nx = panel.nx;
    item.ny = panel.ny;
    item.nz = panel.nz;
    item.area = panel.area;
    out.push_back(item);
  }
  return out;
}

fmm::PrototypeConfig config(int threads) {
  fmm::PrototypeConfig cfg;
  cfg.kernel = fmm::KernelKind::Yukawa;
  cfg.kappa = 0.1;
  cfg.max_leaf_size = 8;
  cfg.expansion_order = 4;
  cfg.separation_theta = 0.6;
  cfg.threads = threads;
  return cfg;
}

void make_density(std::size_t n, double phase, std::vector<double> &phi,
                  std::vector<double> &phin) {
  phi.resize(n);
  phin.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double t = static_cast<double>(i + 1);
    phi[i] = std::sin(0.17 * t + phase);
    phin[i] = 0.7 * std::cos(0.23 * t - phase);
  }
}

double relative_rms(const std::vector<double> &lhs,
                    const std::vector<double> &rhs) {
  expect(lhs.size() == rhs.size(), "vector sizes differ");
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
  const auto panels = make_panels();
  const auto prototype = to_prototype(panels);
  std::vector<double> phi_a, phin_a, phi_b, phin_b;
  make_density(panels.size(), 0.0, phi_a, phin_a);
  make_density(panels.size(), 0.8, phi_b, phin_b);
  TitanExecutor one(1);
  TitanExecutor four(4);

  fmm::FmmM2lAnalyticContext uncached_context(
      config(1), prototype, 1.0 / 80.0, 0.1);
  const auto uncached = uncached_context.apply_local_p4_candidate(
      phi_a, phin_a, one, true, "centroid", "off", &panels,
      "operator-tail-static", 2.0e-3, "off", 256);
  expect(uncached.static_probe_count == 4,
         "static refinement must use four operator probes");
  expect(uncached.plan_hash != 0, "static plan hash must be nonzero");
  expect(uncached.plan_build_count == 1, "static plan must build once");
  fmm::FmmM2lAnalyticContext loose_context(
      config(1), prototype, 1.0 / 80.0, 0.1);
  const auto loose = loose_context.apply_local_p4_candidate(
      phi_a, phin_a, one, true, "centroid", "off", &panels,
      "operator-tail-static", 8.0e-3, "off", 256);
  expect(uncached.refined_cells >= loose.refined_cells,
         "tighter static tolerance must not reduce refinement");

  fmm::FmmM2lAnalyticContext cached_context(
      config(4), prototype, 1.0 / 80.0, 0.1);
  const auto cached_a = cached_context.apply_local_p4_candidate(
      phi_a, phin_a, four, true, "centroid", "off", &panels,
      "operator-tail-static", 2.0e-3, "kernel-jet", 256);
  const auto cached_b = cached_context.apply_local_p4_candidate(
      phi_b, phin_b, four, true, "centroid", "off", &panels,
      "operator-tail-static", 2.0e-3, "kernel-jet", 256);
  expect(cached_a.plan_hash == cached_b.plan_hash,
         "plan hash must be density independent");
  expect(cached_b.plan_build_count == 1,
         "repeated applies must not rebuild the plan");
  expect(cached_b.cache_build_count == 1,
         "repeated applies must not rebuild the cache");
  expect(cached_b.cache_hits > cached_a.cache_hits,
         "second apply must reuse cached operators");
  expect(cached_a.cache_operator_count > 0 && cached_a.cache_bytes > 0,
         "separated fixture must cache far operators");
  expect(relative_rms(cached_a.phi, uncached.phi) <= 1.0e-12 &&
             relative_rms(cached_a.phin, uncached.phin) <= 1.0e-12,
         "cached and uncached candidates must agree");

  fmm::FmmM2lAnalyticContext one_context(
      config(1), prototype, 1.0 / 80.0, 0.1);
  const auto one_result = one_context.apply_local_p4_candidate(
      phi_a, phin_a, one, true, "centroid", "off", &panels,
      "operator-tail-static", 2.0e-3, "kernel-jet", 256);
  expect(one_result.plan_hash == cached_a.plan_hash,
         "plan hash must be thread independent");
  expect(relative_rms(one_result.phi, cached_a.phi) <= 1.0e-12 &&
             relative_rms(one_result.phin, cached_a.phin) <= 1.0e-12,
         "candidate output must be thread deterministic");

  bool budget_failed = false;
  try {
    fmm::FmmM2lAnalyticContext budget_context(
        config(1), prototype, 1.0 / 80.0, 0.1);
    (void)budget_context.apply_local_p4_candidate(
        phi_a, phin_a, one, true, "centroid", "off", &panels,
        "operator-tail-static", 2.0e-3, "kernel-jet", 0);
  } catch (const std::runtime_error &error) {
    budget_failed = std::string(error.what()).find("cache_budget_exceeded") !=
                    std::string::npos;
  }
  expect(budget_failed, "cache budget overflow must fail explicitly");

  auto all_near_config = config(1);
  all_near_config.separation_theta = 0.0;
  fmm::FmmM2lAnalyticContext all_near_context(
      all_near_config, prototype, 1.0 / 80.0, 0.1);
  const auto all_near = all_near_context.apply_local_p4_candidate(
      phi_a, phin_a, one, true, "centroid", "off", &panels,
      "operator-tail-static", 2.0e-3, "kernel-jet", 256);
  expect(all_near.far_pairs == 0 && all_near.cache_operator_count == 0 &&
             all_near.cache_bytes == 0,
         "all-near work must not allocate an M2L cache");

  std::cout << "titan_fmm_m2l_scope_z4m_static_cache_smoke ok\n";
  return 0;
}
