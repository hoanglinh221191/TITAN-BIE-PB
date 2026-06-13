#include "tabi/tabi_fmm_m2l_prototype.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using titan::tabi::fmm_m2l::JcpDirectParityMetrics;
using titan::tabi::fmm_m2l::JcpAnalyticOperatorResult;
using titan::tabi::fmm_m2l::PrototypeConfig;
using titan::tabi::fmm_m2l::PrototypePanel;
using titan::tabi::fmm_m2l::compare_jcp_direct_reference;
using titan::tabi::fmm_m2l::evaluate_jcp_offline_operator_analytic;

using Clock = std::chrono::steady_clock;

void require_true(bool value, const std::string &message) {
  if (!value) {
    throw std::runtime_error(message);
  }
}

double phi_for_panel(std::size_t index) {
  const double t = static_cast<double>(index + 1);
  return 0.7 * std::cos(0.013 * t) + 0.3 * std::sin(0.029 * t);
}

double phin_for_panel(const TitanPBE::Tabi::Panel &panel, std::size_t index) {
  const double t = static_cast<double>(index + 1);
  return 0.65 * std::sin(0.017 * t) + 0.35 * panel.nz;
}

std::vector<TitanPBE::Tabi::Panel> build_small_surface() {
  std::vector<TitanPBE::Tabi::Atom> atoms;
  for (int i = 0; i < 14; ++i) {
    const double t = static_cast<double>(i);
    atoms.push_back({-8.0 + 1.2 * t, 0.55 * std::sin(0.7 * t),
                     0.45 * std::cos(0.5 * t),
                     i % 2 == 0 ? 0.35 : -0.32,
                     1.45 + 0.05 * static_cast<double>(i % 5)});
  }
  TitanPBE::Tabi::Options options;
  options.surface_mode = "vdw";
  options.surface_geometry = "triangle";
  options.diagonal_mode = "solid-angle";
  options.surface_point_count = 12;
  TitanPBE::Tabi::Result result;
  std::vector<TitanPBE::Tabi::Panel> panels =
      TitanPBE::Tabi::build_surface(atoms, options, &result);
  require_true(panels.size() >= 64, "production surface is unexpectedly small");
  return panels;
}

std::vector<PrototypePanel>
to_prototype_panels(const std::vector<TitanPBE::Tabi::Panel> &panels) {
  std::vector<PrototypePanel> out;
  out.reserve(panels.size());
  for (std::size_t i = 0; i < panels.size(); ++i) {
    const auto &panel = panels[i];
    PrototypePanel converted;
    converted.x = panel.x;
    converted.y = panel.y;
    converted.z = panel.z;
    converted.nx = panel.nx;
    converted.ny = panel.ny;
    converted.nz = panel.nz;
    converted.area = panel.area;
    converted.scalar_density = phi_for_panel(i);
    out.push_back(converted);
  }
  return out;
}

std::vector<TitanPBE::Tabi::Atom> read_pqr(const std::string &path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("failed to open PQR fixture: " + path);
  }
  std::vector<TitanPBE::Tabi::Atom> atoms;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    if (!(line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0)) {
      continue;
    }
    std::istringstream stream(line);
    std::vector<std::string> fields;
    for (std::string field; stream >> field;) {
      fields.push_back(field);
    }
    if (fields.size() < 10) {
      throw std::runtime_error("malformed PQR row at line " +
                               std::to_string(line_number));
    }
    try {
      const std::size_t n = fields.size();
      TitanPBE::Tabi::Atom atom;
      atom.x = std::stod(fields[n - 5]);
      atom.y = std::stod(fields[n - 4]);
      atom.z = std::stod(fields[n - 3]);
      atom.q = std::stod(fields[n - 2]);
      atom.radius = std::stod(fields[n - 1]);
      if (!(atom.radius > 0.0)) {
        atom.radius = 0.6;
      }
      atoms.push_back(atom);
    } catch (const std::exception &) {
      throw std::runtime_error("invalid PQR numeric field at line " +
                               std::to_string(line_number));
    }
  }
  if (atoms.empty()) {
    throw std::runtime_error("PQR fixture contains no atoms: " + path);
  }
  return atoms;
}

std::vector<TitanPBE::Tabi::Panel>
sample_panels(const std::vector<TitanPBE::Tabi::Panel> &panels,
              std::size_t max_panels) {
  if (max_panels == 0 || panels.size() <= max_panels) {
    return panels;
  }
  std::vector<TitanPBE::Tabi::Panel> sampled;
  sampled.reserve(max_panels);
  for (std::size_t i = 0; i < max_panels; ++i) {
    const std::size_t index = i * panels.size() / max_panels;
    sampled.push_back(panels[index]);
  }
  return sampled;
}

struct ProductionReference {
  std::vector<double> phi;
  std::vector<double> phin;
  double seconds = 0.0;
};

ProductionReference production_direct_reference(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<PrototypePanel> &prototype, double eps_ratio,
    double kappa, int threads) {
  const std::size_t n = panels.size();
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    unknowns[i] = prototype[i].scalar_density;
    unknowns[n + i] = phin_for_panel(panels[i], i);
  }
  TitanExecutor executor(std::max(1, threads));
  TitanPBE::Tabi::TreeWorkspace direct_fallback;
  std::vector<double> output;
  TitanPBE::Tabi::TreecodeStats stats;
  const auto begin = Clock::now();
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, kappa, 0.0, 0.0, true, direct_fallback,
      "centroid", "off", executor, output, &stats);
  const auto end = Clock::now();
  require_true(output.size() == 2 * n, "production JCP output size");
  ProductionReference reference;
  reference.phi.resize(n);
  reference.phin.resize(n);
  reference.seconds = std::chrono::duration<double>(end - begin).count();
  for (std::size_t i = 0; i < n; ++i) {
    reference.phi[i] = -output[i];
    reference.phin[i] = -output[n + i];
  }
  return reference;
}

struct ErrorMetrics {
  double phi_rel = 0.0;
  double phin_rel = 0.0;
  double aggregate_rel = 0.0;
  double max_abs = 0.0;
  bool finite = true;
};

ErrorMetrics compare_vectors(const std::vector<double> &ref_phi,
                             const std::vector<double> &ref_phin,
                             const std::vector<double> &got_phi,
                             const std::vector<double> &got_phin) {
  require_true(ref_phi.size() == got_phi.size(), "phi vector size mismatch");
  require_true(ref_phin.size() == got_phin.size(), "phin vector size mismatch");
  double phi_err2 = 0.0;
  double phi_ref2 = 0.0;
  double phin_err2 = 0.0;
  double phin_ref2 = 0.0;
  double err2 = 0.0;
  double ref2 = 0.0;
  ErrorMetrics out;
  for (std::size_t i = 0; i < ref_phi.size(); ++i) {
    const double pe = got_phi[i] - ref_phi[i];
    const double ne = got_phin[i] - ref_phin[i];
    phi_err2 += pe * pe;
    phi_ref2 += ref_phi[i] * ref_phi[i];
    phin_err2 += ne * ne;
    phin_ref2 += ref_phin[i] * ref_phin[i];
    err2 += pe * pe + ne * ne;
    ref2 += ref_phi[i] * ref_phi[i] + ref_phin[i] * ref_phin[i];
    out.max_abs = std::max({out.max_abs, std::abs(pe), std::abs(ne)});
    out.finite = out.finite && std::isfinite(ref_phi[i]) &&
                 std::isfinite(ref_phin[i]) && std::isfinite(got_phi[i]) &&
                 std::isfinite(got_phin[i]);
  }
  out.phi_rel = std::sqrt(phi_err2 / std::max(phi_ref2, 1.0e-30));
  out.phin_rel = std::sqrt(phin_err2 / std::max(phin_ref2, 1.0e-30));
  out.aggregate_rel = std::sqrt(err2 / std::max(ref2, 1.0e-30));
  return out;
}

std::string classify_scope_i(const JcpDirectParityMetrics &parity,
                             const ErrorMetrics &offline,
                             const JcpAnalyticOperatorResult &result,
                             double production_direct_seconds) {
  if (!parity.finite_outputs || !offline.finite || !result.finite_outputs) {
    return "nonfinite";
  }
  if (parity.aggregate_rel_rms_error > 1.0e-10) {
    return "production_parity_fail";
  }
  if (offline.aggregate_rel > 1.0e-1 ||
      std::max(offline.phi_rel, offline.phin_rel) > 2.0e-1) {
    return "single_matvec_accuracy_fail";
  }
  if (result.far_pairs == 0) {
    return "no_far_work";
  }
  if (!(result.apply_seconds < production_direct_seconds)) {
    return "no_speedup";
  }
  return "pass";
}

double median(std::vector<double> values) {
  if (values.empty()) {
    return 0.0;
  }
  std::sort(values.begin(), values.end());
  const std::size_t middle = values.size() / 2;
  return values.size() % 2 == 0
             ? 0.5 * (values[middle - 1] + values[middle])
             : values[middle];
}

void run_generated_surface_smoke() {
  constexpr double eps_ratio = 80.0;
  constexpr double kappa = 0.2;
  std::vector<TitanPBE::Tabi::Panel> panels = build_small_surface();
  std::vector<PrototypePanel> prototype = to_prototype_panels(panels);
  const ProductionReference production =
      production_direct_reference(panels, prototype, eps_ratio, kappa, 4);
  const JcpDirectParityMetrics parity = compare_jcp_direct_reference(
      prototype, production.phi, production.phin, eps_ratio, kappa);
  PrototypeConfig config;
  config.expansion_order = 4;
  config.separation_theta = 0.7;
  config.kappa = kappa;
  config.threads = 4;
  const JcpAnalyticOperatorResult result =
      evaluate_jcp_offline_operator_analytic(config, prototype, eps_ratio, kappa);
  const ErrorMetrics offline =
      compare_vectors(production.phi, production.phin, result.phi, result.phin);
  require_true(parity.finite_outputs, "direct parity outputs must be finite");
  require_true(parity.aggregate_rel_rms_error <= 1.0e-10,
               "production direct parity must match analytic JCP reference");
  require_true(offline.finite, "offline result must be finite");
  require_true(result.phi.size() == panels.size() &&
                   result.phin.size() == panels.size(),
               "offline single-matvec result must preserve vector shape");
  require_true(result.near_pairs + result.far_pairs > 0,
               "offline single-matvec result must do work");
}

int run_fixture_probe(int argc, char **argv) {
  std::string fixture = "real-surface";
  std::string pqr_path;
  int surface_points = 8;
  std::size_t max_panels = 1024;
  int order = 2;
  double theta = 0.7;
  int threads = 8;
  int repeats = 5;
  bool header = false;
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    auto value = [&](const std::string &prefix) {
      return arg.rfind(prefix, 0) == 0 ? arg.substr(prefix.size()) : std::string{};
    };
    if (arg == "--csv-header") {
      header = true;
    } else if (const std::string v = value("--fixture="); !v.empty()) {
      fixture = v;
    } else if (const std::string v = value("--pqr="); !v.empty()) {
      pqr_path = v;
    } else if (const std::string v = value("--surface-points="); !v.empty()) {
      surface_points = std::stoi(v);
    } else if (const std::string v = value("--max-panels="); !v.empty()) {
      max_panels = static_cast<std::size_t>(std::stoull(v));
    } else if (const std::string v = value("--order="); !v.empty()) {
      order = std::stoi(v);
    } else if (const std::string v = value("--theta="); !v.empty()) {
      theta = std::stod(v);
    } else if (const std::string v = value("--threads="); !v.empty()) {
      threads = std::stoi(v);
    } else if (const std::string v = value("--repeats="); !v.empty()) {
      repeats = std::stoi(v);
    } else {
      throw std::invalid_argument("unknown argument: " + arg);
    }
  }
  require_true(!pqr_path.empty(), "--pqr is required");

  const std::vector<TitanPBE::Tabi::Atom> atoms = read_pqr(pqr_path);
  TitanPBE::Tabi::Options options;
  options.surface_mode = "vdw";
  options.surface_geometry = "triangle";
  options.diagonal_mode = "solid-angle";
  options.surface_point_count = std::max(4, surface_points);
  TitanPBE::Tabi::Result surface_result;
  TitanExecutor surface_executor(std::max(1, threads));
  const auto surface_begin = Clock::now();
  const std::vector<TitanPBE::Tabi::Panel> full_panels =
      TitanPBE::Tabi::build_surface(atoms, options, &surface_result,
                                    &surface_executor);
  const auto surface_end = Clock::now();
  const std::vector<TitanPBE::Tabi::Panel> panels =
      sample_panels(full_panels, max_panels);
  require_true(panels.size() >= 64, "sampled real surface is unexpectedly small");
  const std::vector<PrototypePanel> prototype = to_prototype_panels(panels);

  constexpr double eps_ratio = 80.0;
  constexpr double kappa = 0.2;
  const ProductionReference production =
      production_direct_reference(panels, prototype, eps_ratio, kappa, threads);
  const JcpDirectParityMetrics parity = compare_jcp_direct_reference(
      prototype, production.phi, production.phin, eps_ratio, kappa);
  PrototypeConfig config;
  config.max_leaf_size = 64;
  config.expansion_order = order;
  config.separation_theta = theta;
  config.kappa = kappa;
  config.threads = threads;
  require_true(repeats >= 1, "--repeats must be at least 1");

  JcpAnalyticOperatorResult result;
  std::vector<double> setup_times;
  std::vector<double> apply_times;
  for (int repeat = 0; repeat < repeats; ++repeat) {
    JcpAnalyticOperatorResult run =
        evaluate_jcp_offline_operator_analytic(config, prototype, eps_ratio, kappa);
    if (repeat == 0) {
      result = run;
    }
    setup_times.push_back(run.setup_seconds);
    apply_times.push_back(run.apply_seconds);
  }
  result.setup_seconds = median(setup_times);
  result.apply_seconds = median(apply_times);
  const ErrorMetrics offline =
      compare_vectors(production.phi, production.phin, result.phi, result.phin);
  const std::string gate =
      classify_scope_i(parity, offline, result, production.seconds);
  const double speedup = result.apply_seconds > 0.0
                             ? production.seconds / result.apply_seconds
                             : 0.0;
  const double surface_seconds =
      std::chrono::duration<double>(surface_end - surface_begin).count();

  if (header) {
    std::cout
        << "fixture,atoms,surface_panels,sampled_panels,surface_seconds,"
           "production_direct_seconds,order,theta,repeats,parity_aggregate_rel_rms,"
           "offline_phi_rel_rms,offline_phin_rel_rms,offline_aggregate_rel_rms,"
           "offline_max_abs,near_pairs,far_pairs,setup_seconds,apply_seconds,"
           "speedup_vs_production_direct,scope_i_gate\n";
  }
  std::cout << std::setprecision(12) << fixture << ',' << atoms.size() << ','
            << full_panels.size() << ',' << panels.size() << ','
            << surface_seconds << ',' << production.seconds << ',' << order
            << ',' << theta << ',' << repeats << ','
            << parity.aggregate_rel_rms_error << ',' << offline.phi_rel << ','
            << offline.phin_rel << ',' << offline.aggregate_rel << ','
            << offline.max_abs << ',' << result.near_pairs << ','
            << result.far_pairs << ',' << result.setup_seconds << ','
            << result.apply_seconds << ',' << speedup << ',' << gate << '\n';
  return gate == "pass" ? 0 : 2;
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc > 1 && std::string(argv[1]) == "--fixture-probe") {
      return run_fixture_probe(argc, argv);
    }
    run_generated_surface_smoke();
    std::cout << "PASS titan_fmm_m2l_scope_i_probe\n";
    return 0;
  } catch (const std::exception &ex) {
    std::cerr << "FAIL titan_fmm_m2l_scope_i_probe: " << ex.what() << '\n';
    return 1;
  }
}
