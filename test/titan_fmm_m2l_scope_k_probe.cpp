#include "tabi/tabi_fmm_m2l_prototype.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using titan::tabi::fmm_m2l::JcpGuardedGmresMetrics;
using titan::tabi::fmm_m2l::PrototypeConfig;
using titan::tabi::fmm_m2l::PrototypePanel;
using titan::tabi::fmm_m2l::run_jcp_fmm_m2l_guarded_gmres_probe;

void require_true(bool value, const std::string &message) {
  if (!value) {
    throw std::runtime_error(message);
  }
}

std::vector<PrototypePanel> make_generated_panels(std::size_t count) {
  std::vector<PrototypePanel> panels;
  panels.reserve(count);
  const double golden = 2.3999632297286533222;
  for (std::size_t i = 0; i < count; ++i) {
    const double u = -1.0 + 2.0 * (static_cast<double>(i) + 0.5) /
                                static_cast<double>(count);
    const double r = std::sqrt(std::max(0.0, 1.0 - u * u));
    const double theta = golden * static_cast<double>(i);
    PrototypePanel panel;
    panel.x = 2.2 * r * std::cos(theta);
    panel.y = 1.5 * r * std::sin(theta);
    panel.z = 2.8 * u;
    double nx = panel.x / (2.2 * 2.2);
    double ny = panel.y / (1.5 * 1.5);
    double nz = panel.z / (2.8 * 2.8);
    const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
    panel.nx = norm > 0.0 ? nx / norm : 0.0;
    panel.ny = norm > 0.0 ? ny / norm : 0.0;
    panel.nz = norm > 0.0 ? nz / norm : 1.0;
    panel.area = 4.0 * 3.14159265358979323846 * 2.2 * 1.5 /
                 static_cast<double>(count) *
                 (0.8 + 0.4 * std::abs(panel.nz));
    panels.push_back(panel);
  }
  return panels;
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
    const std::size_t n = fields.size();
    TitanPBE::Tabi::Atom atom;
    atom.x = std::stod(fields[n - 5]);
    atom.y = std::stod(fields[n - 4]);
    atom.z = std::stod(fields[n - 3]);
    atom.q = std::stod(fields[n - 2]);
    atom.radius = std::max(0.6, std::stod(fields[n - 1]));
    atoms.push_back(atom);
  }
  if (atoms.empty()) {
    throw std::runtime_error("PQR fixture contains no atoms: " + path);
  }
  return atoms;
}

std::vector<PrototypePanel>
sample_panels(const std::vector<TitanPBE::Tabi::Panel> &panels,
              std::size_t max_panels) {
  const std::size_t n =
      max_panels == 0 ? panels.size() : std::min(max_panels, panels.size());
  std::vector<PrototypePanel> sampled;
  sampled.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const TitanPBE::Tabi::Panel &src = panels[i * panels.size() / n];
    PrototypePanel dst;
    dst.x = src.x;
    dst.y = src.y;
    dst.z = src.z;
    dst.nx = src.nx;
    dst.ny = src.ny;
    dst.nz = src.nz;
    dst.area = src.area;
    sampled.push_back(dst);
  }
  return sampled;
}

std::string classify_scope_k(const JcpGuardedGmresMetrics &m,
                             double direct_residual_gate) {
  if (!m.finite_outputs) {
    return "nonfinite";
  }
  if (!m.adapter_residual_decreased) {
    return "no_residual_decrease";
  }
  if (m.adapter_final_residual > 1.0e-6) {
    return "adapter_residual_fail";
  }
  if (m.direct_true_residual > direct_residual_gate) {
    return "direct_true_residual_fail";
  }
  if (m.far_pairs == 0) {
    return "no_far_work";
  }
  if (m.solution_rel_error > 5.0e-2) {
    return "solution_error_fail";
  }
  if (m.speedup_vs_direct_matvec < 1.2) {
    return "no_speedup";
  }
  return "pass";
}

void run_generated_smoke() {
  PrototypeConfig config;
  config.kernel = titan::tabi::fmm_m2l::KernelKind::Yukawa;
  config.kappa = 0.2;
  config.max_leaf_size = 32;
  config.expansion_order = 4;
  config.separation_theta = 0.8;
  config.threads = 4;
  const std::vector<PrototypePanel> panels = make_generated_panels(96);
  const JcpGuardedGmresMetrics metrics =
      run_jcp_fmm_m2l_guarded_gmres_probe(config, panels, 80.0, 0.2, 16, 48,
                                          1.0e-6);
  require_true(metrics.finite_outputs, "generated GMRES outputs finite");
  require_true(metrics.direct_converged, "generated direct GMRES converged");
  require_true(metrics.adapter_residual_decreased,
               "generated adapter residual decreased");
}

int run_fixture_probe(int argc, char **argv) {
  std::string fixture = "research-medium";
  std::string pqr_path;
  int surface_points = 4;
  std::size_t max_panels = 1024;
  std::size_t leaf_size = 64;
  double theta = 0.8;
  int threads = 4;
  int restart = 24;
  int max_iterations = 72;
  double tolerance = 1.0e-6;
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    auto value = [&](const std::string &prefix) {
      return arg.rfind(prefix, 0) == 0 ? arg.substr(prefix.size())
                                       : std::string{};
    };
    if (const std::string v = value("--fixture="); !v.empty()) {
      fixture = v;
    } else if (const std::string v = value("--pqr="); !v.empty()) {
      pqr_path = v;
    } else if (const std::string v = value("--surface-points="); !v.empty()) {
      surface_points = std::stoi(v);
    } else if (const std::string v = value("--max-panels="); !v.empty()) {
      max_panels = static_cast<std::size_t>(std::stoull(v));
    } else if (const std::string v = value("--leaf-size="); !v.empty()) {
      leaf_size = static_cast<std::size_t>(std::stoull(v));
    } else if (const std::string v = value("--theta="); !v.empty()) {
      theta = std::stod(v);
    } else if (const std::string v = value("--threads="); !v.empty()) {
      threads = std::stoi(v);
    } else if (const std::string v = value("--restart="); !v.empty()) {
      restart = std::stoi(v);
    } else if (const std::string v = value("--max-iterations="); !v.empty()) {
      max_iterations = std::stoi(v);
    } else if (const std::string v = value("--tolerance="); !v.empty()) {
      tolerance = std::stod(v);
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
  TitanExecutor executor(std::max(1, threads));
  const auto surface_begin = std::chrono::steady_clock::now();
  const std::vector<TitanPBE::Tabi::Panel> production_panels =
      TitanPBE::Tabi::build_surface(atoms, options, &surface_result, &executor);
  const auto surface_end = std::chrono::steady_clock::now();
  const std::vector<PrototypePanel> panels =
      sample_panels(production_panels, max_panels);
  require_true(panels.size() >= 64, "sampled surface is too small");

  PrototypeConfig config;
  config.kernel = titan::tabi::fmm_m2l::KernelKind::Yukawa;
  config.kappa = 0.2;
  config.max_leaf_size = static_cast<int>(leaf_size);
  config.expansion_order = 4;
  config.separation_theta = theta;
  config.threads = threads;
  const JcpGuardedGmresMetrics metrics =
      run_jcp_fmm_m2l_guarded_gmres_probe(config, panels, 80.0, 0.2, restart,
                                          max_iterations, tolerance);
  const double surface_seconds =
      std::chrono::duration<double>(surface_end - surface_begin).count();
  const double direct_gate = fixture == "research-large" ? 2.0e-2 : 5.0e-3;
  const std::string gate = classify_scope_k(metrics, direct_gate);
  std::cout
      << "fixture,atoms,surface_panels,sampled_panels,surface_seconds,"
         "leaf_size,theta,restart,max_iterations,tolerance,direct_iterations,"
         "adapter_iterations,direct_final_residual,adapter_final_residual,"
         "direct_true_residual,solution_rel_error,direct_matvec_ms,"
         "adapter_matvec_ms,speedup_vs_direct_matvec,near_pairs,far_pairs,"
         "adapter_setup_ms,adapter_apply_ms,finite_outputs,scope_k_gate\n";
  std::cout << std::setprecision(12) << fixture << ',' << atoms.size() << ','
            << production_panels.size() << ',' << panels.size() << ','
            << surface_seconds << ',' << leaf_size << ',' << theta << ','
            << restart << ',' << max_iterations << ',' << tolerance << ','
            << metrics.direct_iterations << ',' << metrics.adapter_iterations
            << ',' << metrics.direct_final_residual << ','
            << metrics.adapter_final_residual << ','
            << metrics.direct_true_residual << ',' << metrics.solution_rel_error
            << ',' << metrics.direct_matvec_ms << ','
            << metrics.adapter_matvec_ms << ','
            << metrics.speedup_vs_direct_matvec << ',' << metrics.near_pairs
            << ',' << metrics.far_pairs << ',' << metrics.adapter_setup_ms
            << ',' << metrics.adapter_apply_ms << ','
            << (metrics.finite_outputs ? "true" : "false") << ',' << gate
            << '\n';
  return gate == "pass" ? 0 : 2;
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc > 1 && std::string(argv[1]) == "--fixture-probe") {
      return run_fixture_probe(argc, argv);
    }
    run_generated_smoke();
    std::cout << "PASS titan_fmm_m2l_scope_k_probe\n";
    return 0;
  } catch (const std::exception &ex) {
    std::cerr << "FAIL titan_fmm_m2l_scope_k_probe: " << ex.what() << '\n';
    return 1;
  }
}
