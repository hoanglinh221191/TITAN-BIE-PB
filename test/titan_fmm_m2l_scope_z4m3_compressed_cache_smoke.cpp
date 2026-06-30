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

namespace {
namespace fmm = titan::tabi::fmm_m2l;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_fmm_m2l_scope_z4m3_compressed_cache_smoke: "
            << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

double relative_rms(const fmm::CartesianJcpLocalP4 &lhs,
                    const fmm::CartesianJcpLocalP4 &rhs) {
  double diff_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t field = 0; field < lhs.fields.size(); ++field) {
    for (std::size_t term = 0; term < lhs.fields[field].size(); ++term) {
      const double diff = lhs.fields[field][term] - rhs.fields[field][term];
      diff_sq += diff * diff;
      ref_sq += rhs.fields[field][term] * rhs.fields[field][term];
    }
  }
  return std::sqrt(diff_sq / std::max(1.0e-300, ref_sq));
}
} // namespace

int main() {
  expect(fmm::cartesian_compressed_green_jet_term_count() == 364,
         "max degree-11 compact jet storage must contain 364 coefficients");
  expect(fmm::cartesian_compressed_green_jet_term_count_order(4) == 220,
         "p4 compressed Green jet must contain 220 active coefficients");
  expect(fmm::cartesian_compressed_green_operator_bytes_order(4) == 3520,
         "p4 compressed cache operator byte count must remain 3520");
  expect(sizeof(fmm::CartesianJcpCompressedGreenOperatorP4) == 5824,
         "max compressed operator wrapper must contain two 364-double jets");

  auto panels = fmm::make_tabi_panel_like_fixture(96, 29);
  for (std::size_t i = 0; i < panels.size(); ++i)
    panels[i].x += 4.0 * static_cast<double>(i / 32);
  std::vector<double> phi(panels.size());
  std::vector<double> phin(panels.size());
  std::vector<std::size_t> indices(panels.size());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    phi[i] = std::sin(0.13 * static_cast<double>(i + 1));
    phin[i] = std::cos(0.19 * static_cast<double>(i + 1));
    indices[i] = i;
  }
  const auto moments = fmm::build_cartesian_jcp_moments(
      panels, phi, phin, indices, 0.1, -0.2, 0.05, 4);
  const auto legacy = fmm::build_cartesian_jcp_m2l_operator_p4(
      0.1, -0.2, 0.05, 4.2, 1.4, -2.1, 1.0 / 80.0, 0.17);
  const auto compressed = fmm::build_cartesian_jcp_compressed_green_operator_p4(
      0.1, -0.2, 0.05, 4.2, 1.4, -2.1, 0.17);

  const auto legacy_local =
      fmm::apply_cartesian_jcp_m2l_operator_p4_order(moments, legacy, 4);
  const auto compressed_local =
      fmm::apply_cartesian_jcp_compressed_green_operator_p4_order(
          moments, compressed, 1.0 / 80.0, 4);
  expect(legacy.finite_outputs && compressed_local.finite_outputs,
         "operators and local fields must be finite");
  expect(relative_rms(compressed_local, legacy_local) <= 1.0e-12,
         "compressed Green jets must reproduce legacy K1-K4 locals");

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
  fmm::PrototypeConfig config;
  config.kernel = fmm::KernelKind::Yukawa;
  config.kappa = 0.17;
  config.max_leaf_size = 8;
  config.expansion_order = 4;
  config.separation_theta = 0.6;
  config.threads = 4;
  TitanExecutor executor(4);

  fmm::FmmM2lAnalyticContext legacy_context(
      config, panels, 1.0 / 80.0, 0.17);
  const auto legacy_candidate = legacy_context.apply_local_p4_candidate(
      phi, phin, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "kernel-jet", 256, "scalar",
      "metric-cache");
  fmm::FmmM2lAnalyticContext compressed_context(
      config, panels, 1.0 / 80.0, 0.17);
  const auto compressed_a = compressed_context.apply_local_p4_candidate(
      phi, phin, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "scalar", "metric-cache");
  std::vector<double> phi_b(phi.size());
  std::vector<double> phin_b(phin.size());
  for (std::size_t i = 0; i < phi.size(); ++i) {
    phi_b[i] = std::sin(0.31 * static_cast<double>(i + 1));
    phin_b[i] = std::cos(0.27 * static_cast<double>(i + 1));
  }
  const auto compressed_b = compressed_context.apply_local_p4_candidate(
      phi_b, phin_b, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "scalar", "metric-cache");

  auto vector_relative_rms = [](const std::vector<double> &lhs,
                                const std::vector<double> &rhs) {
    double diff_sq = 0.0;
    double ref_sq = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
      const double diff = lhs[i] - rhs[i];
      diff_sq += diff * diff;
      ref_sq += rhs[i] * rhs[i];
    }
    return std::sqrt(diff_sq / std::max(1.0e-300, ref_sq));
  };
  expect(compressed_a.plan_hash == legacy_candidate.plan_hash &&
             compressed_a.far_pairs == legacy_candidate.far_pairs &&
             compressed_a.direct_fallback_pairs ==
                 legacy_candidate.direct_fallback_pairs,
         "cache layout must not change the static refinement plan");
  expect(vector_relative_rms(compressed_a.phi, legacy_candidate.phi) <=
                 1.0e-12 &&
             vector_relative_rms(compressed_a.phin,
                                 legacy_candidate.phin) <= 1.0e-12,
         "compressed and legacy candidate vectors must agree");
  expect(compressed_a.cache_layout == "green-jet-compressed",
         "compressed cache layout label is incorrect");
  expect(compressed_a.cache_operator_bytes == 3520,
         "compressed cache operator byte count is incorrect");
  expect(compressed_a.cache_compression_ratio >= 2.8,
         "compressed cache ratio is below 2.8x");
  expect(compressed_a.cache_far_pair_offsets > 0,
         "compressed cache must expose deterministic far-pair offsets");
  expect(compressed_b.cache_build_count == 1 &&
             compressed_b.cache_hits > compressed_a.cache_hits,
         "repeated applies must reuse the compressed cache");

  auto one_config = config;
  one_config.threads = 1;
  TitanExecutor one_executor(1);
  fmm::FmmM2lAnalyticContext one_context(
      one_config, panels, 1.0 / 80.0, 0.17);
  const auto one_result = one_context.apply_local_p4_candidate(
      phi, phin, one_executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "scalar", "metric-cache");
  expect(one_result.plan_hash == compressed_a.plan_hash &&
             one_result.cache_far_pair_offsets ==
                 compressed_a.cache_far_pair_offsets &&
             vector_relative_rms(one_result.phi, compressed_a.phi) <=
                 1.0e-12 &&
             vector_relative_rms(one_result.phin, compressed_a.phin) <=
                 1.0e-12,
         "compressed cache layout and output must be deterministic for nt=1/4");

  bool budget_failed = false;
  try {
    fmm::FmmM2lAnalyticContext budget_context(
        config, panels, 1.0 / 80.0, 0.17);
    (void)budget_context.apply_local_p4_candidate(
        phi, phin, executor, true, "centroid", "off", &production,
        "operator-tail-static", 2.0e-3, "green-jet-compressed", 0,
        "scalar", "metric-cache");
  } catch (const std::runtime_error &error) {
    budget_failed = std::string(error.what()).find("cache_budget_exceeded") !=
                    std::string::npos;
  }
  expect(budget_failed, "compressed budget overflow must fail explicitly");

  config.separation_theta = 0.0;
  fmm::FmmM2lAnalyticContext all_near_context(
      config, panels, 1.0 / 80.0, 0.17);
  const auto all_near = all_near_context.apply_local_p4_candidate(
      phi, phin, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "scalar", "metric-cache");
  expect(all_near.far_pairs == 0 && all_near.cache_bytes == 0 &&
             all_near.cache_operator_count == 0,
         "all-near geometry must allocate no compressed cache");

  std::cout << "titan_fmm_m2l_scope_z4m3_compressed_cache_smoke ok\n";
  return 0;
}
