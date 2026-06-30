#include "tabi/tabi_fmm_m2l_cartesian.h"
#include "tabi/tabi_fmm_m2l_local.h"
#include "tabi/tabi_fmm_m2l_prototype.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {
namespace fmm = titan::tabi::fmm_m2l;

[[noreturn]] void fail(const std::string &message) {
  std::cerr << "FAIL titan_fmm_m2l_scope_z4xc_p5_smoke: " << message
            << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

double local_relative_rms(const fmm::CartesianJcpLocalP4 &lhs,
                          const fmm::CartesianJcpLocalP4 &rhs,
                          std::size_t term_count) {
  double diff_sq = 0.0;
  double ref_sq = 0.0;
  for (std::size_t field = 0; field < lhs.fields.size(); ++field) {
    for (std::size_t term = 0; term < term_count; ++term) {
      const double diff = lhs.fields[field][term] - rhs.fields[field][term];
      diff_sq += diff * diff;
      ref_sq += rhs.fields[field][term] * rhs.fields[field][term];
    }
  }
  return std::sqrt(diff_sq / std::max(1.0e-300, ref_sq));
}

double vector_relative_rms(const std::vector<double> &lhs,
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

} // namespace

int main() {
  expect(fmm::cartesian_moment_term_count(5) == 56,
         "p5 moment/local term count must be 56");
  expect(fmm::cartesian_local_term_count() == 56,
         "local storage must reserve 56 terms for p5");
  expect(fmm::cartesian_local_kernel_term_count() == 286,
         "degree-10 local kernel jet must contain 286 coefficients");
  expect(fmm::cartesian_compressed_green_jet_term_count() == 364,
         "degree-11 compressed Green jet must contain 364 coefficients");
  expect(fmm::cartesian_compressed_green_jet_term_count_order(4) == 220 &&
             fmm::cartesian_compressed_green_operator_bytes_order(4) == 3520,
         "p4 compressed cache storage must remain compact");
  expect(fmm::cartesian_compressed_green_jet_term_count_order(5) == 364 &&
             fmm::cartesian_compressed_green_operator_bytes_order(5) == 5824,
         "p5 compressed cache storage must use degree-11 Green jets");
  expect(sizeof(fmm::CartesianJcpCompressedGreenOperatorP4) == 5824,
         "compressed p5 operator must contain exactly two 364-double jets");

  auto panels = fmm::make_tabi_panel_like_fixture(112, 37);
  for (std::size_t i = 0; i < panels.size(); ++i)
    panels[i].x += 3.5 * static_cast<double>(i / 28);
  std::vector<double> phi(panels.size());
  std::vector<double> phin(panels.size());
  std::vector<std::size_t> indices(panels.size());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    phi[i] = std::sin(0.11 * static_cast<double>(i + 1));
    phin[i] = std::cos(0.17 * static_cast<double>(i + 1));
    indices[i] = i;
  }

  const auto moments = fmm::build_cartesian_jcp_moments(
      panels, phi, phin, indices, 0.1, -0.2, 0.05, 5);
  const auto kernel = fmm::build_cartesian_jcp_m2l_operator_p4(
      0.1, -0.2, 0.05, 5.2, 1.4, -2.1, 1.0 / 80.0, 0.17);
  const auto compressed = fmm::build_cartesian_jcp_compressed_green_operator_p4(
      0.1, -0.2, 0.05, 5.2, 1.4, -2.1, 0.17);
  const auto kernel_local =
      fmm::apply_cartesian_jcp_m2l_operator_p4_order(moments, kernel, 5);
  const auto compressed_local =
      fmm::apply_cartesian_jcp_compressed_green_operator_p4_order(
          moments, compressed, 1.0 / 80.0, 5);
  expect(kernel.finite_outputs && kernel_local.finite_outputs &&
             compressed_local.finite_outputs,
         "p5 kernel and compressed locals must be finite");
  expect(local_relative_rms(compressed_local, kernel_local, 56) <= 1.0e-12,
         "compressed p5 Green jets must reproduce p5 K1-K4 locals");

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
  config.expansion_order = 5;
  config.separation_theta = 0.6;
  config.threads = 4;
  TitanExecutor executor(4);

  fmm::FmmM2lAnalyticContext context(config, panels, 1.0 / 80.0, 0.17);
  const auto first = context.apply_local_p4_candidate(
      phi, phin, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "scalar", "metric-cache");
  std::vector<double> phi_b(phi.size());
  std::vector<double> phin_b(phin.size());
  for (std::size_t i = 0; i < phi.size(); ++i) {
    phi_b[i] = std::sin(0.31 * static_cast<double>(i + 1));
    phin_b[i] = std::cos(0.23 * static_cast<double>(i + 1));
  }
  const auto second = context.apply_local_p4_candidate(
      phi_b, phin_b, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "scalar", "metric-cache");
  expect(first.local_source_order == 5 && first.local_target_order == 5 &&
             first.local_kernel_derivative_order == 10,
         "p5 candidate must report source/target/kernel orders");
  expect(first.cache_operator_bytes == 5824 &&
             first.cache_layout == "green-jet-compressed",
         "p5 compressed cache byte accounting is incorrect");
  expect(second.plan_build_count == 1 && second.cache_build_count == 1 &&
             second.cache_hits > first.cache_hits,
         "p5 repeated applies must reuse static plan and compressed cache");
  expect(vector_relative_rms(first.phi, first.phi) <= 1.0e-15 &&
             vector_relative_rms(second.phin, second.phin) <= 1.0e-15,
         "p5 candidate vectors must be finite");

  auto one_config = config;
  one_config.threads = 1;
  TitanExecutor one_executor(1);
  fmm::FmmM2lAnalyticContext one_context(one_config, panels, 1.0 / 80.0,
                                         0.17);
  const auto one = one_context.apply_local_p4_candidate(
      phi, phin, one_executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "scalar", "metric-cache");
  expect(one.plan_hash == first.plan_hash &&
             vector_relative_rms(one.phi, first.phi) <= 1.0e-12 &&
             vector_relative_rms(one.phin, first.phin) <= 1.0e-12,
         "p5 static plan/output must be deterministic for nt=1/4");

  std::cout << "titan_fmm_m2l_scope_z4xc_p5_smoke ok\n";
  return 0;
}
