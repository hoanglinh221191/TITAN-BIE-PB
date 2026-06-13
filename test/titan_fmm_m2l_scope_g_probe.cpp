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
using titan::tabi::fmm_m2l::JcpDerivativeBridgeMetrics;
using titan::tabi::fmm_m2l::JcpK1BridgeMetrics;
using titan::tabi::fmm_m2l::PrototypeConfig;
using titan::tabi::fmm_m2l::PrototypePanel;
using titan::tabi::fmm_m2l::compare_jcp_direct_reference;
using titan::tabi::fmm_m2l::run_jcp_derivative_bridge_probe;
using titan::tabi::fmm_m2l::run_jcp_k1_bridge_probe;

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

std::string classify_scope_g(const JcpDirectParityMetrics &parity,
                             const JcpK1BridgeMetrics &k1,
                             const JcpDerivativeBridgeMetrics &derivative) {
  if (!parity.finite_outputs || !k1.finite_outputs ||
      !derivative.finite_outputs) {
    return "nonfinite";
  }
  if (parity.aggregate_rel_rms_error > 1.0e-10) {
    return "production_parity_fail";
  }
  if (k1.rel_rms_error > 5.0e-2) {
    return "k1_accuracy_fail";
  }
  const double max_reference = std::max(
      {derivative.k2_reference_rms, derivative.k3_reference_rms,
       derivative.k4_reference_rms});
  const double significance_floor = std::max(1.0e-12, 1.0e-3 * max_reference);
  double significant_max_relative = 0.0;
  if (derivative.k2_reference_rms >= significance_floor) {
    significant_max_relative =
        std::max(significant_max_relative, derivative.k2_rel_rms_error);
  }
  if (derivative.k3_reference_rms >= significance_floor) {
    significant_max_relative =
        std::max(significant_max_relative, derivative.k3_rel_rms_error);
  }
  if (derivative.k4_reference_rms >= significance_floor) {
    significant_max_relative =
        std::max(significant_max_relative, derivative.k4_rel_rms_error);
  }
  if (derivative.aggregate_rel_rms_error > 1.0e-1 ||
      significant_max_relative > 2.0e-1) {
    return "derivative_accuracy_fail";
  }
  if (k1.far_pairs + derivative.far_pairs == 0) {
    return "no_far_work";
  }
  if (!(k1.apply_seconds + derivative.apply_seconds <
        k1.direct_seconds + derivative.direct_seconds)) {
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
  require_true(parity.finite_outputs, "direct parity outputs must be finite");
  require_true(parity.aggregate_rel_rms_error <= 1.0e-10,
               "production direct parity must match analytic JCP reference");
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
  const ProductionReference production = production_direct_reference(
      panels, prototype, eps_ratio, kappa, threads);
  const JcpDirectParityMetrics parity = compare_jcp_direct_reference(
      prototype, production.phi, production.phin, eps_ratio, kappa);
  PrototypeConfig config;
  config.max_leaf_size = 64;
  config.expansion_order = order;
  config.separation_theta = theta;
  config.kappa = kappa;
  config.threads = threads;
  require_true(repeats >= 1, "--repeats must be at least 1");
  JcpK1BridgeMetrics k1;
  JcpDerivativeBridgeMetrics derivative;
  std::vector<double> k1_setup;
  std::vector<double> k1_apply;
  std::vector<double> k1_direct;
  std::vector<double> derivative_setup;
  std::vector<double> derivative_apply;
  std::vector<double> derivative_direct;
  for (int repeat = 0; repeat < repeats; ++repeat) {
    JcpK1BridgeMetrics k1_run = run_jcp_k1_bridge_probe(config, prototype);
    JcpDerivativeBridgeMetrics derivative_run =
        run_jcp_derivative_bridge_probe(config, prototype);
    if (repeat == 0) {
      k1 = k1_run;
      derivative = derivative_run;
    }
    k1_setup.push_back(k1_run.setup_seconds);
    k1_apply.push_back(k1_run.apply_seconds);
    k1_direct.push_back(k1_run.direct_seconds);
    derivative_setup.push_back(derivative_run.setup_seconds);
    derivative_apply.push_back(derivative_run.apply_seconds);
    derivative_direct.push_back(derivative_run.direct_seconds);
  }
  k1.setup_seconds = median(k1_setup);
  k1.apply_seconds = median(k1_apply);
  k1.direct_seconds = median(k1_direct);
  derivative.setup_seconds = median(derivative_setup);
  derivative.apply_seconds = median(derivative_apply);
  derivative.direct_seconds = median(derivative_direct);
  const double max_reference = std::max(
      {derivative.k2_reference_rms, derivative.k3_reference_rms,
       derivative.k4_reference_rms});
  const double significance_floor = std::max(1.0e-12, 1.0e-3 * max_reference);
  double significant_max_relative = 0.0;
  if (derivative.k2_reference_rms >= significance_floor) {
    significant_max_relative =
        std::max(significant_max_relative, derivative.k2_rel_rms_error);
  }
  if (derivative.k3_reference_rms >= significance_floor) {
    significant_max_relative =
        std::max(significant_max_relative, derivative.k3_rel_rms_error);
  }
  if (derivative.k4_reference_rms >= significance_floor) {
    significant_max_relative =
        std::max(significant_max_relative, derivative.k4_rel_rms_error);
  }
  const double combined_apply = k1.apply_seconds + derivative.apply_seconds;
  const double combined_direct = k1.direct_seconds + derivative.direct_seconds;
  const double speedup = combined_apply > 0.0 ? combined_direct / combined_apply : 0.0;
  const std::string gate = classify_scope_g(parity, k1, derivative);

  if (header) {
    std::cout << "fixture,atoms,surface_panels,sampled_panels,surface_seconds,"
                 "production_direct_seconds,order,theta,repeats,parity_phi_rel_rms,"
                 "parity_phin_rel_rms,parity_aggregate_rel_rms,k1_rel_rms,"
                 "derivative_aggregate_rel_rms,derivative_max_channel_rel_rms,"
                 "derivative_significant_max_rel_rms,k2_rel_rms,k3_rel_rms,"
                 "k4_rel_rms,k2_reference_rms,k3_reference_rms,k4_reference_rms,"
                 "near_pairs,far_pairs,setup_seconds,apply_seconds,"
                 "prototype_direct_seconds,speedup_vs_direct,scope_g_gate\n";
  }
  std::cout << fixture << ',' << atoms.size() << ',' << full_panels.size() << ','
            << panels.size() << ',' << std::setprecision(12)
            << std::chrono::duration<double>(surface_end - surface_begin).count()
            << ',' << production.seconds << ',' << order << ',' << theta << ','
            << repeats << ','
            << parity.phi_rel_rms_error << ',' << parity.phin_rel_rms_error << ','
            << parity.aggregate_rel_rms_error << ',' << k1.rel_rms_error << ','
            << derivative.aggregate_rel_rms_error << ','
            << derivative.max_channel_rel_rms_error << ','
            << significant_max_relative << ','
            << derivative.k2_rel_rms_error << ','
            << derivative.k3_rel_rms_error << ','
            << derivative.k4_rel_rms_error << ','
            << derivative.k2_reference_rms << ','
            << derivative.k3_reference_rms << ','
            << derivative.k4_reference_rms << ','
            << (k1.near_pairs + derivative.near_pairs) << ','
            << (k1.far_pairs + derivative.far_pairs) << ','
            << (k1.setup_seconds + derivative.setup_seconds) << ','
            << combined_apply << ',' << combined_direct << ',' << speedup << ','
            << gate << '\n';
  return gate == "pass" ? 0 : 2;
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc > 1 && std::string(argv[1]) == "--fixture-probe") {
      return run_fixture_probe(argc, argv);
    }
    run_generated_surface_smoke();
    std::cout << "PASS titan_fmm_m2l_scope_g_probe\n";
    return 0;
  } catch (const std::exception &ex) {
    std::cerr << "FAIL titan_fmm_m2l_scope_g_probe: " << ex.what() << '\n';
    return 1;
  }
}
