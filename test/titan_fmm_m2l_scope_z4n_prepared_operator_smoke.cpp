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
  std::cerr << "FAIL titan_fmm_m2l_scope_z4n_prepared_operator_smoke: "
            << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

double relative_rms(const std::vector<double> &lhs,
                    const std::vector<double> &rhs) {
  if (lhs.size() != rhs.size())
    return INFINITY;
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
production_panels(const std::vector<fmm::PrototypePanel> &panels) {
  std::vector<TitanPBE::Tabi::Panel> out(panels.size());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    out[i].x = panels[i].x;
    out[i].y = panels[i].y;
    out[i].z = panels[i].z;
    out[i].nx = panels[i].nx;
    out[i].ny = panels[i].ny;
    out[i].nz = panels[i].nz;
    out[i].area = panels[i].area;
    out[i].diagonal = 1.0;
  }
  return out;
}

fmm::JcpLocalP4OperatorOptions options() {
  fmm::JcpLocalP4OperatorOptions value;
  value.near_integration = "centroid";
  value.quad_buffer = "off";
  value.refine_tolerance = 2.0e-3;
  value.cache_mode = "green-jet-compressed";
  value.cache_budget_mb = 256;
  value.fallback_mode = "scalar";
  value.plan_builder = "metric-cache";
  return value;
}
} // namespace

int main() {
  auto panels = fmm::make_tabi_panel_like_fixture(96, 41);
  for (std::size_t i = 0; i < panels.size(); ++i)
    panels[i].x += 4.0 * static_cast<double>(i / 32);
  const auto production = production_panels(panels);

  fmm::PrototypeConfig config;
  config.kernel = fmm::KernelKind::Yukawa;
  config.kappa = 0.17;
  config.max_leaf_size = 8;
  config.expansion_order = 4;
  config.separation_theta = 0.6;
  config.threads = 4;
  TitanExecutor executor(4);

  fmm::FmmM2lAnalyticContext prepared_context(
      config, panels, 1.0 / 80.0, 0.17);
  bool unprepared_failed = false;
  try {
    std::vector<double> zero(panels.size(), 0.0);
    (void)prepared_context.apply_prepared_local_p4_operator(
        zero, zero, executor, true);
  } catch (const std::logic_error &error) {
    unprepared_failed =
        std::string(error.what()).find("not_prepared") != std::string::npos;
  }
  expect(unprepared_failed, "apply before prepare must fail clearly");

  const auto prepare = prepared_context.prepare_local_p4_operator(
      production, executor, options());
  expect(prepare.plan_build_count == 1 && prepare.cache_build_count == 1,
         "prepare must build plan and cache exactly once");
  expect(prepare.cache_layout == "green-jet-compressed" &&
             prepare.cache_bytes > 0,
         "prepare must create the compressed Green-jet cache");

  std::vector<double> phi(panels.size());
  std::vector<double> phin(panels.size());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    phi[i] = std::sin(0.13 * static_cast<double>(i + 1));
    phin[i] = std::cos(0.19 * static_cast<double>(i + 1));
  }

  const auto prepared = prepared_context.apply_prepared_local_p4_operator(
      phi, phin, executor, true);
  fmm::FmmM2lAnalyticContext one_shot_context(
      config, panels, 1.0 / 80.0, 0.17);
  const auto one_shot = one_shot_context.apply_local_p4_candidate(
      phi, phin, executor, true, "centroid", "off", &production,
      "operator-tail-static", 2.0e-3, "green-jet-compressed", 256,
      "scalar", "metric-cache");
  expect(relative_rms(prepared.phi, one_shot.phi) <= 1.0e-12 &&
             relative_rms(prepared.phin, one_shot.phin) <= 1.0e-12,
         "prepared and one-shot p4 outputs must agree");

  std::vector<double> unknowns(2 * phi.size(), 0.0);
  for (std::size_t i = 0; i < phi.size(); ++i) {
    unknowns[i] = phi[i];
    unknowns[phi.size() + i] = phin[i];
  }
  std::vector<double> system_out;
  fmm::assemble_jcp_local_p4_system_matvec(
      prepared, unknowns, 0.5, 0.5 / 80.0, true, system_out);
  std::vector<double> expected_system(system_out.size(), 0.0);
  for (std::size_t i = 0; i < phi.size(); ++i) {
    expected_system[i] = 0.5 * phi[i] - prepared.phi[i];
    expected_system[phi.size() + i] =
        (0.5 / 80.0) * phin[i] - prepared.phin[i];
  }
  expect(relative_rms(system_out, expected_system) <= 1.0e-15,
         "coupled p4 system assembly must preserve LHS and sign convention");

  auto phi_b = phi;
  auto phin_b = phin;
  auto phi_ab = phi;
  auto phin_ab = phin;
  for (std::size_t i = 0; i < panels.size(); ++i) {
    phi_b[i] = std::sin(0.31 * static_cast<double>(i + 1));
    phin_b[i] = std::cos(0.27 * static_cast<double>(i + 1));
    phi_ab[i] = 0.25 * phi[i] - 0.75 * phi_b[i];
    phin_ab[i] = 0.25 * phin[i] - 0.75 * phin_b[i];
  }
  const auto result_a = prepared_context.apply_prepared_local_p4_operator(
      phi, phin, executor, true);
  const auto result_b = prepared_context.apply_prepared_local_p4_operator(
      phi_b, phin_b, executor, true);
  const auto result_ab = prepared_context.apply_prepared_local_p4_operator(
      phi_ab, phin_ab, executor, true);
  std::vector<double> expected_phi(phi.size());
  std::vector<double> expected_phin(phin.size());
  for (std::size_t i = 0; i < phi.size(); ++i) {
    expected_phi[i] = 0.25 * result_a.phi[i] - 0.75 * result_b.phi[i];
    expected_phin[i] = 0.25 * result_a.phin[i] - 0.75 * result_b.phin[i];
  }
  expect(relative_rms(result_ab.phi, expected_phi) <= 1.0e-12 &&
             relative_rms(result_ab.phin, expected_phin) <= 1.0e-12,
         "prepared p4 operator must remain linear for GMRES");

  fmm::JcpLocalP4CandidateResult last = result_ab;
  for (int repeat = 0; repeat < 50; ++repeat) {
    const double scale = 1.0 + 0.01 * static_cast<double>(repeat);
    for (std::size_t i = 0; i < phi.size(); ++i) {
      phi_ab[i] = scale * phi[i];
      phin_ab[i] = scale * phin[i];
    }
    last = prepared_context.apply_prepared_local_p4_operator(
        phi_ab, phin_ab, executor, true);
  }
  expect(last.plan_build_count == 1 && last.cache_build_count == 1,
         "varying GMRES densities must not rebuild plan or cache");
  expect(last.cache_hits > prepared.cache_hits,
         "prepared cache hits must grow across repeated matvecs");

  auto expanded_options = options();
  expanded_options.refine_probe_count = 8;
  fmm::FmmM2lAnalyticContext expanded_context(
      config, panels, 1.0 / 80.0, 0.17);
  const auto expanded_prepare = expanded_context.prepare_local_p4_operator(
      production, executor, expanded_options);
  expect(expanded_prepare.static_probe_count == 8,
         "expanded rademacher plan must report requested probe count");
  expect(expanded_prepare.plan_hash != prepare.plan_hash,
         "changing the static probe count must rebuild a distinct plan");

  auto mixed_options = options();
  mixed_options.refine_probe_count = 8;
  mixed_options.refine_probe_mode = "mixed";
  fmm::FmmM2lAnalyticContext mixed_context(
      config, panels, 1.0 / 80.0, 0.17);
  const auto mixed_prepare = mixed_context.prepare_local_p4_operator(
      production, executor, mixed_options);
  expect(mixed_prepare.refine_probe_mode == "mixed" &&
             mixed_prepare.static_probe_count == 8 &&
             mixed_prepare.solution_probe_used,
         "mixed mode must report solution-probe diagnostics");
  const auto mixed_apply = mixed_context.apply_prepared_local_p4_operator(
      phi, phin, executor, true);
  expect(mixed_apply.plan_build_count == 1 && mixed_apply.cache_build_count == 1,
         "mixed prepared operator must reuse its plan and cache");

  std::cout << "titan_fmm_m2l_scope_z4n_prepared_operator_smoke ok\n";
  return 0;
}
