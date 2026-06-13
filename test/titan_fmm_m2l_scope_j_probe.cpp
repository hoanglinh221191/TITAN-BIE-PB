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

std::vector<TitanPBE::Tabi::Panel> make_panels(std::size_t count) {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(count);
  const double golden = 2.3999632297286533222;
  for (std::size_t i = 0; i < count; ++i) {
    const double u = -1.0 + 2.0 * (static_cast<double>(i) + 0.5) /
                                static_cast<double>(count);
    const double r = std::sqrt(std::max(0.0, 1.0 - u * u));
    const double theta = golden * static_cast<double>(i);
    TitanPBE::Tabi::Panel panel;
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
                 static_cast<double>(count) * (0.8 + 0.4 * std::abs(panel.nz));
    panels.push_back(panel);
  }
  return panels;
}

std::vector<double> make_unknowns(const std::vector<TitanPBE::Tabi::Panel> &panels) {
  std::vector<double> unknowns(2 * panels.size(), 0.0);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    unknowns[i] = phi_for_panel(i);
    unknowns[panels.size() + i] = phin_for_panel(panels[i], i);
  }
  return unknowns;
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
    atom.radius = std::stod(fields[n - 1]);
    if (!(atom.radius > 0.0)) {
      atom.radius = 0.6;
    }
    atoms.push_back(atom);
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
    sampled.push_back(panels[i * panels.size() / max_panels]);
  }
  return sampled;
}

struct ErrorMetrics {
  double rel_rms = 0.0;
  double phi_rel = 0.0;
  double phin_rel = 0.0;
  double max_abs = 0.0;
  bool finite = true;
};

ErrorMetrics compare_output(const std::vector<double> &reference,
                            const std::vector<double> &candidate) {
  require_true(reference.size() == candidate.size(), "output size mismatch");
  require_true(reference.size() % 2 == 0, "JCP output should be 2N");
  const std::size_t n = reference.size() / 2;
  ErrorMetrics out;
  double err2 = 0.0;
  double ref2 = 0.0;
  double phi_err2 = 0.0;
  double phi_ref2 = 0.0;
  double phin_err2 = 0.0;
  double phin_ref2 = 0.0;
  for (std::size_t i = 0; i < reference.size(); ++i) {
    const double e = candidate[i] - reference[i];
    err2 += e * e;
    ref2 += reference[i] * reference[i];
    out.max_abs = std::max(out.max_abs, std::abs(e));
    out.finite = out.finite && std::isfinite(reference[i]) &&
                 std::isfinite(candidate[i]);
    if (i < n) {
      phi_err2 += e * e;
      phi_ref2 += reference[i] * reference[i];
    } else {
      phin_err2 += e * e;
      phin_ref2 += reference[i] * reference[i];
    }
  }
  out.rel_rms = std::sqrt(err2 / std::max(ref2, 1.0e-30));
  out.phi_rel = std::sqrt(phi_err2 / std::max(phi_ref2, 1.0e-30));
  out.phin_rel = std::sqrt(phin_err2 / std::max(phin_ref2, 1.0e-30));
  return out;
}

double compare_scalar_output(const std::vector<double> &reference,
                             const std::vector<double> &candidate) {
  require_true(reference.size() == candidate.size(),
               "scalar output size mismatch");
  double err2 = 0.0;
  double ref2 = 0.0;
  for (std::size_t i = 0; i < reference.size(); ++i) {
    require_true(std::isfinite(reference[i]) && std::isfinite(candidate[i]),
                 "scalar output must be finite");
    const double error = candidate[i] - reference[i];
    err2 += error * error;
    ref2 += reference[i] * reference[i];
  }
  return std::sqrt(err2 / std::max(ref2, 1.0e-30));
}

double median(std::vector<double> values) {
  if (values.empty()) {
    return 0.0;
  }
  std::sort(values.begin(), values.end());
  const std::size_t middle = values.size() / 2;
  return values.size() % 2 == 0 ? 0.5 * (values[middle - 1] + values[middle])
                                : values[middle];
}

struct MatvecRun {
  std::vector<double> out;
  TitanPBE::Tabi::FmmM2lSingleMatvecStats fmm_stats;
  TitanPBE::Tabi::TreecodeStats tree_stats;
  double seconds = 0.0;
};

MatvecRun run_direct(const std::vector<TitanPBE::Tabi::Panel> &panels,
                     const std::vector<double> &unknowns, double eps_ratio,
                     double kappa, int threads) {
  TitanExecutor executor(std::max(1, threads));
  TitanPBE::Tabi::TreeWorkspace empty_tree;
  MatvecRun run;
  const auto begin = Clock::now();
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, kappa, 0.0, 0.0, true, empty_tree,
      "centroid", "off", executor, run.out, &run.tree_stats);
  const auto end = Clock::now();
  run.seconds = std::chrono::duration<double>(end - begin).count();
  return run;
}

MatvecRun run_adapter(const std::vector<TitanPBE::Tabi::Panel> &panels,
                      const std::vector<double> &unknowns, double eps_ratio,
                      double kappa, std::size_t leaf_size, double theta,
                      int threads) {
  TitanExecutor executor(std::max(1, threads));
  MatvecRun run;
  const auto begin = Clock::now();
  TitanPBE::Tabi::apply_jcp_fmm_m2l_single_matvec(
      panels, unknowns, eps_ratio, kappa, 0.0, 0.0, leaf_size, theta, executor,
      run.out, &run.fmm_stats);
  const auto end = Clock::now();
  run.seconds = std::chrono::duration<double>(end - begin).count();
  return run;
}

MatvecRun run_hidden_route(const std::vector<TitanPBE::Tabi::Panel> &panels,
                           const std::vector<double> &unknowns,
                           double eps_ratio, double kappa,
                           std::size_t leaf_size, double theta, int threads) {
  TitanExecutor executor(std::max(1, threads));
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(panels, executor, leaf_size, theta);
  tree.fmm_m2l = "single-matvec-analytic";
  MatvecRun run;
  const auto begin = Clock::now();
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, kappa, 0.0, 0.0, true, tree, "centroid",
      "off", executor, run.out, &run.tree_stats);
  const auto end = Clock::now();
  run.seconds = std::chrono::duration<double>(end - begin).count();
  return run;
}

std::string classify_scope_j(const ErrorMetrics &err, const MatvecRun &adapter,
                             double direct_seconds, double threshold) {
  if (!err.finite || !adapter.fmm_stats.finite_outputs) {
    return "nonfinite";
  }
  if (err.rel_rms > threshold) {
    return "accuracy_fail";
  }
  if (adapter.fmm_stats.far_pairs == 0) {
    return "no_far_work";
  }
  const double speedup =
      adapter.fmm_stats.apply_ms > 0.0
          ? direct_seconds / (adapter.fmm_stats.apply_ms * 1.0e-3)
          : 0.0;
  if (speedup < 1.2) {
    return "no_speedup";
  }
  return "pass";
}

void run_generated_surface_smoke() {
  constexpr double eps_ratio = 80.0;
  constexpr double kappa = 0.2;
  const std::vector<TitanPBE::Tabi::Panel> panels = make_panels(192);
  const std::vector<double> unknowns = make_unknowns(panels);
  const MatvecRun direct = run_direct(panels, unknowns, eps_ratio, kappa, 4);
  const MatvecRun adapter = run_adapter(panels, unknowns, eps_ratio, kappa, 32,
                                        0.8, 4);
  const MatvecRun hidden = run_hidden_route(panels, unknowns, eps_ratio, kappa,
                                            32, 0.8, 4);
  const ErrorMetrics adapter_err = compare_output(direct.out, adapter.out);
  const ErrorMetrics hidden_err = compare_output(direct.out, hidden.out);
  require_true(adapter.out.size() == unknowns.size(), "adapter output size");
  require_true(adapter_err.finite, "adapter output finite");
  require_true(hidden_err.finite, "hidden route output finite");
  require_true(adapter_err.rel_rms <= 1.0e-3,
               "generated adapter/direct aggregate RMS gate");
  require_true(hidden_err.rel_rms <= 1.0e-3,
               "generated hidden route/direct aggregate RMS gate");
  require_true(hidden.tree_stats.fmm_m2l_mode == "single-matvec-analytic",
               "hidden route should fill FMM/M2L stats");
}

void run_generated_poisson_smoke() {
  constexpr double eps_ratio = 80.0;
  constexpr double kappa = 0.0;
  constexpr double lhs_phi = 40.5;
  const std::vector<TitanPBE::Tabi::Panel> panels = make_panels(192);
  std::vector<double> unknowns(panels.size(), 0.0);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    unknowns[i] = phi_for_panel(i);
  }
  TitanExecutor executor(4);
  TitanPBE::Tabi::TreeWorkspace direct_tree;
  std::vector<double> direct;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, kappa, lhs_phi, 0.0, false, direct_tree,
      "centroid", "off", executor, direct, nullptr);

  TitanPBE::Tabi::TreeWorkspace fmm_tree =
      TitanPBE::Tabi::build_tree_workspace(panels, executor, 32, 0.8);
  fmm_tree.fmm_m2l = "single-matvec-analytic";
  std::vector<double> fmm;
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, kappa, lhs_phi, 0.0, false, fmm_tree,
      "centroid", "off", executor, fmm, &stats);
  require_true(fmm.size() == panels.size(),
               "Poisson FMM/M2L output should be N");
  require_true(compare_scalar_output(direct, fmm) <= 1.0e-3,
               "generated Poisson FMM/M2L direct parity gate");
  require_true(stats.fmm_m2l_mode == "single-matvec-analytic",
               "Poisson hidden route should fill FMM/M2L stats");
}

int run_fixture_probe(int argc, char **argv) {
  std::string fixture = "real-surface";
  std::string pqr_path;
  int surface_points = 4;
  std::size_t max_panels = 2048;
  std::size_t leaf_size = 64;
  double theta = 0.8;
  int threads = 4;
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
    } else if (const std::string v = value("--leaf-size="); !v.empty()) {
      leaf_size = static_cast<std::size_t>(std::stoull(v));
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
  require_true(repeats >= 1, "--repeats must be at least 1");

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

  constexpr double eps_ratio = 80.0;
  constexpr double kappa = 0.2;
  const std::vector<double> unknowns = make_unknowns(panels);
  const MatvecRun direct = run_direct(panels, unknowns, eps_ratio, kappa, threads);

  MatvecRun adapter;
  std::vector<double> adapter_total_times;
  std::vector<double> adapter_setup_ms;
  std::vector<double> adapter_apply_ms;
  for (int repeat = 0; repeat < repeats; ++repeat) {
    MatvecRun run =
        run_adapter(panels, unknowns, eps_ratio, kappa, leaf_size, theta, threads);
    if (repeat == 0) {
      adapter = run;
    }
    adapter_total_times.push_back(run.seconds);
    adapter_setup_ms.push_back(run.fmm_stats.setup_ms);
    adapter_apply_ms.push_back(run.fmm_stats.apply_ms);
  }
  adapter.seconds = median(adapter_total_times);
  adapter.fmm_stats.setup_ms = median(adapter_setup_ms);
  adapter.fmm_stats.apply_ms = median(adapter_apply_ms);
  const ErrorMetrics err = compare_output(direct.out, adapter.out);
  const double threshold = fixture == "generated" ? 1.0e-3 : 1.0e-2;
  const std::string gate = classify_scope_j(err, adapter, direct.seconds,
                                            threshold);
  const double speedup =
      adapter.fmm_stats.apply_ms > 0.0
          ? direct.seconds / (adapter.fmm_stats.apply_ms * 1.0e-3)
          : 0.0;
  const double surface_seconds =
      std::chrono::duration<double>(surface_end - surface_begin).count();

  if (header) {
    std::cout
        << "fixture,atoms,surface_panels,sampled_panels,surface_seconds,"
           "production_direct_seconds,leaf_size,theta,repeats,"
           "adapter_phi_rel_rms,adapter_phin_rel_rms,adapter_aggregate_rel_rms,"
           "adapter_max_abs,near_pairs,far_pairs,setup_ms,apply_ms,"
           "adapter_total_seconds,speedup_vs_production_direct,scope_j_gate\n";
  }
  std::cout << std::setprecision(12) << fixture << ',' << atoms.size() << ','
            << full_panels.size() << ',' << panels.size() << ','
            << surface_seconds << ',' << direct.seconds << ',' << leaf_size
            << ',' << theta << ',' << repeats << ',' << err.phi_rel << ','
            << err.phin_rel << ',' << err.rel_rms << ',' << err.max_abs << ','
            << adapter.fmm_stats.near_pairs << ',' << adapter.fmm_stats.far_pairs
            << ',' << adapter.fmm_stats.setup_ms << ','
            << adapter.fmm_stats.apply_ms << ',' << adapter.seconds << ','
            << speedup << ',' << gate << '\n';
  return gate == "pass" ? 0 : 2;
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc > 1 && std::string(argv[1]) == "--fixture-probe") {
      return run_fixture_probe(argc, argv);
    }
    run_generated_surface_smoke();
    run_generated_poisson_smoke();
    std::cout << "PASS titan_fmm_m2l_scope_j_probe\n";
    return 0;
  } catch (const std::exception &ex) {
    std::cerr << "FAIL titan_fmm_m2l_scope_j_probe: " << ex.what() << '\n';
    return 1;
  }
}
