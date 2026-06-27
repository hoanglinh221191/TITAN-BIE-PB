#include "tabi/tabi_fmm_m2l_cartesian.h"
#include "tabi/tabi_fmm_m2l_local.h"
#include "tabi/tabi_fmm_m2l_prototype.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {
namespace fmm = titan::tabi::fmm_m2l;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_fmm_m2l_scope_z4xd_tiled_cache_smoke: "
            << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

double relative_rms(const std::vector<double> &lhs,
                    const std::vector<double> &rhs) {
  double diff_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    const double diff = lhs[i] - rhs[i];
    diff_sq += diff * diff;
    ref_sq += rhs[i] * rhs[i];
  }
  return std::sqrt(diff_sq / std::max(1.0e-300, ref_sq));
}

std::vector<TitanPBE::Tabi::Panel>
make_production_panels(const std::vector<fmm::PrototypePanel> &panels) {
  std::vector<TitanPBE::Tabi::Panel> production(panels.size());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    production[i].x = panels[i].x;
    production[i].y = panels[i].y;
    production[i].z = panels[i].z;
    production[i].nx = panels[i].nx;
    production[i].ny = panels[i].ny;
    production[i].nz = panels[i].nz;
    production[i].area = panels[i].area;
    production[i].diagonal = 1.0;
  }
  return production;
}

} // namespace

int main() {
  expect(fmm::cartesian_compressed_green_operator_bytes_order(4) == 3520,
         "p4 compressed operator byte count must remain compact");
  expect(fmm::cartesian_compressed_green_operator_bytes_order(5) == 5824,
         "p5 compressed operator byte count must be 5824");

  auto panels = fmm::make_tabi_panel_like_fixture(192, 41);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    panels[i].x += 5.0 * static_cast<double>(i / 48);
    panels[i].y += 1.75 * static_cast<double>((i / 24) % 2);
  }
  const auto production = make_production_panels(panels);

  std::vector<double> phi(panels.size());
  std::vector<double> phin(panels.size());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    phi[i] = std::sin(0.071 * static_cast<double>(i + 1));
    phin[i] = std::cos(0.097 * static_cast<double>(i + 1));
  }

  fmm::PrototypeConfig config;
  config.kernel = fmm::KernelKind::Yukawa;
  config.kappa = 0.17;
  config.max_leaf_size = 6;
  config.expansion_order = 5;
  config.separation_theta = 0.6;
  config.threads = 4;
  TitanExecutor executor(4);

  fmm::FmmM2lAnalyticContext uncached_context(config, panels, 1.0 / 80.0,
                                              0.17);
  const auto uncached = uncached_context.apply_local_p4_candidate(
      phi, phin, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "off", 256, "scalar",
      "metric-cache");

  fmm::FmmM2lAnalyticContext global_context(config, panels, 1.0 / 80.0,
                                            0.17);
  const auto global = global_context.apply_local_p4_candidate(
      phi, phin, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "scalar", "metric-cache");

  fmm::FmmM2lAnalyticContext tiled_context(config, panels, 1.0 / 80.0,
                                           0.17);
  const auto tiled = tiled_context.apply_local_p4_candidate(
      phi, phin, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-tiled", 1, "scalar",
      "metric-cache");

  expect(tiled.local_source_order == 5 && tiled.local_target_order == 5,
         "p5 tiled candidate must report effective order 5");
  expect(tiled.cache_layout == "green-jet-tiled",
         "tiled cache layout label is incorrect");
  expect(tiled.cache_operator_bytes == 5824,
         "p5 tiled cache operator byte count is incorrect");
  expect(tiled.cache_bytes == 4 * 5824,
         "p5 tiled cache bytes must be bounded by worker-local tiles");
  expect(global.cache_bytes > tiled.cache_bytes * 10,
         "global cache must materialize substantially more storage");
  expect(tiled.cache_tile_count == 4 && tiled.cache_tile_builds == tiled.far_pairs,
         "tiled cache telemetry must report worker tiles and per-pair builds");
  expect(tiled.plan_hash == global.plan_hash &&
             tiled.far_pairs == global.far_pairs &&
             tiled.direct_fallback_pairs == global.direct_fallback_pairs,
         "tiled cache must not change the static refinement plan");
  expect(relative_rms(tiled.phi, uncached.phi) <= 1.0e-12 &&
             relative_rms(tiled.phin, uncached.phin) <= 1.0e-12,
         "tiled and uncached p5 outputs must agree");
  expect(relative_rms(tiled.phi, global.phi) <= 1.0e-12 &&
             relative_rms(tiled.phin, global.phin) <= 1.0e-12,
         "tiled and global p5 outputs must agree");

  std::vector<double> phi_b(phi.size());
  std::vector<double> phin_b(phin.size());
  for (std::size_t i = 0; i < phi.size(); ++i) {
    phi_b[i] = std::sin(0.131 * static_cast<double>(i + 1));
    phin_b[i] = std::cos(0.149 * static_cast<double>(i + 1));
  }
  const auto tiled_b = tiled_context.apply_local_p4_candidate(
      phi_b, phin_b, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-tiled", 1, "scalar",
      "metric-cache");
  expect(tiled_b.plan_build_count == 1 && tiled_b.cache_build_count <= 1,
         "repeated tiled applies must not rebuild plan/global cache");
  expect(tiled_b.plan_hash == tiled.plan_hash &&
             tiled_b.cache_bytes == tiled.cache_bytes,
         "tiled cache metadata must be density independent");

  auto one_config = config;
  one_config.threads = 1;
  TitanExecutor one_executor(1);
  fmm::FmmM2lAnalyticContext one_context(one_config, panels, 1.0 / 80.0,
                                         0.17);
  const auto one = one_context.apply_local_p4_candidate(
      phi, phin, one_executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-tiled", 1, "scalar",
      "metric-cache");
  expect(one.plan_hash == tiled.plan_hash &&
             relative_rms(one.phi, tiled.phi) <= 1.0e-12 &&
             relative_rms(one.phin, tiled.phin) <= 1.0e-12,
         "tiled p5 output must be deterministic for nt=1/4");

  std::cout << "titan_fmm_m2l_scope_z4xd_tiled_cache_smoke ok\n";
  return 0;
}
