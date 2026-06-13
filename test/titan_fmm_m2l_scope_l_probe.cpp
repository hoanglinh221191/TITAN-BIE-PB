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

std::vector<TitanPBE::Tabi::Atom> make_generated_atoms(std::size_t count) {
  std::vector<TitanPBE::Tabi::Atom> atoms;
  atoms.reserve(count);
  const double golden = 2.3999632297286533222;
  for (std::size_t i = 0; i < count; ++i) {
    const double u = -1.0 + 2.0 * (static_cast<double>(i) + 0.5) /
                                static_cast<double>(count);
    const double r = std::sqrt(std::max(0.0, 1.0 - u * u));
    const double theta = golden * static_cast<double>(i);
    TitanPBE::Tabi::Atom atom;
    atom.x = 5.0 * r * std::cos(theta);
    atom.y = 4.0 * r * std::sin(theta);
    atom.z = 6.0 * u;
    atom.q = (i % 2 == 0 ? 0.15 : -0.15);
    atom.radius = 1.5;
    atoms.push_back(atom);
  }
  return atoms;
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

std::vector<TitanPBE::Tabi::Atom>
sample_atoms(const std::vector<TitanPBE::Tabi::Atom> &atoms,
             std::size_t atom_cap) {
  const std::size_t n =
      atom_cap == 0 ? atoms.size() : std::min(atom_cap, atoms.size());
  std::vector<TitanPBE::Tabi::Atom> sampled;
  sampled.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    sampled.push_back(atoms[i * atoms.size() / n]);
  }
  return sampled;
}

TitanPBE::Tabi::Options make_options(bool fmm, int surface_points,
                                     int max_iter, int restart,
                                     double tolerance) {
  TitanPBE::Tabi::Options options;
  options.formulation = "jcp2013-treecode";
  options.sum_mode = fmm ? "treecode" : "direct";
  options.fmm_m2l = fmm ? "single-matvec-analytic" : "off";
  options.surface_mode = "vdw";
  options.surface_geometry = "triangle";
  options.diagonal_mode = "solid-angle";
  options.surface_point_count = surface_points;
  options.tree_leaf_size = 64;
  options.tree_opening_theta = 0.8;
  options.near_integration = "centroid";
  options.quad_buffer = "off";
  options.preconditioner = "diagonal";
  options.max_iter = max_iter;
  options.gmres_restart = restart;
  options.residual_tol = tolerance;
  options.true_residual_probe = "restart";
  options.true_residual_probe_scope = "sampled";
  options.true_residual_sample_count = 512;
  return options;
}

struct SolveRow {
  std::string fixture;
  std::size_t input_atoms = 0;
  std::size_t sampled_atoms = 0;
  double direct_seconds = 0.0;
  double fmm_seconds = 0.0;
  TitanPBE::Tabi::Result direct;
  TitanPBE::Tabi::Result fmm;
  double energy_delta_kj_mol = 0.0;
  double energy_delta_abs_kj_mol = 0.0;
  double energy_delta_rel = 0.0;
  double matvec_speedup = 0.0;
  std::string gate;
};

std::string classify_scope_l(const SolveRow &row) {
  const bool finite =
      std::isfinite(row.direct.energy_kj_mol) &&
      std::isfinite(row.fmm.energy_kj_mol) &&
      std::isfinite(row.direct.final_residual) &&
      std::isfinite(row.fmm.final_residual);
  if (!finite) {
    return "nonfinite";
  }
  if (row.direct.status != "ok" || row.fmm.status != "ok") {
    return "status_fail";
  }
  if (row.fmm.fmm_m2l_mode != "single-matvec-analytic") {
    return "telemetry_missing";
  }
  if (row.fmm.fmm_m2l_far_pairs == 0) {
    return "no_far_work";
  }
  const double residual_gate =
      std::max(1.0e-4, 10.0 * std::max(row.direct.true_residual_rel, 0.0));
  if (row.fmm.true_residual_rel > residual_gate) {
    return "true_residual_fail";
  }
  const double abs_energy = std::abs(row.direct.energy_kj_mol);
  const double energy_gate =
      row.fixture == "research-large" ? std::max(5.0, 0.08 * abs_energy)
                                      : std::max(2.0, 0.05 * abs_energy);
  if (row.energy_delta_abs_kj_mol > energy_gate) {
    return "energy_delta_fail";
  }
  if (row.matvec_speedup < 1.2) {
    return "no_speedup";
  }
  return "pass";
}

SolveRow run_solve_pair(const std::string &fixture,
                        const std::vector<TitanPBE::Tabi::Atom> &input_atoms,
                        std::size_t atom_cap, int surface_points,
                        int max_iter, int restart, double tolerance,
                        int threads) {
  SolveRow row;
  row.fixture = fixture;
  row.input_atoms = input_atoms.size();
  const std::vector<TitanPBE::Tabi::Atom> atoms =
      sample_atoms(input_atoms, atom_cap);
  row.sampled_atoms = atoms.size();
  TitanExecutor direct_executor(std::max(1, threads));
  TitanExecutor fmm_executor(std::max(1, threads));
  const auto direct_begin = Clock::now();
  row.direct = TitanPBE::Tabi::solve_linear_pb(
      atoms, make_options(false, surface_points, max_iter, restart, tolerance),
      direct_executor);
  const auto direct_end = Clock::now();
  const auto fmm_begin = Clock::now();
  row.fmm = TitanPBE::Tabi::solve_linear_pb(
      atoms, make_options(true, surface_points, max_iter, restart, tolerance),
      fmm_executor);
  const auto fmm_end = Clock::now();
  row.direct_seconds =
      std::chrono::duration<double>(direct_end - direct_begin).count();
  row.fmm_seconds =
      std::chrono::duration<double>(fmm_end - fmm_begin).count();
  row.energy_delta_kj_mol = row.fmm.energy_kj_mol - row.direct.energy_kj_mol;
  row.energy_delta_abs_kj_mol = std::abs(row.energy_delta_kj_mol);
  row.energy_delta_rel = row.energy_delta_abs_kj_mol /
                         std::max(1.0e-30, std::abs(row.direct.energy_kj_mol));
  row.matvec_speedup =
      row.fmm.matvec_ms_avg > 0.0 ? row.direct.matvec_ms_avg / row.fmm.matvec_ms_avg
                                  : 0.0;
  row.gate = classify_scope_l(row);
  return row;
}

void write_header() {
  std::cout
      << "fixture,input_atoms,sampled_atoms,direct_status,fmm_status,"
         "direct_panels,fmm_panels,direct_iterations,fmm_iterations,"
         "direct_final_residual,fmm_final_residual,direct_true_residual,"
         "fmm_true_residual,direct_energy_kj_mol,fmm_energy_kj_mol,"
         "energy_delta_kj_mol,energy_delta_abs_kj_mol,energy_delta_rel,"
         "direct_matvec_ms_avg,fmm_matvec_ms_avg,matvec_speedup,"
         "fmm_m2l_mode,fmm_near_pairs,fmm_far_pairs,fmm_setup_ms,"
         "fmm_apply_ms,direct_seconds,fmm_seconds,scope_l_gate\n";
}

void write_row(const SolveRow &row) {
  std::cout << std::setprecision(12) << row.fixture << ',' << row.input_atoms
            << ',' << row.sampled_atoms << ',' << row.direct.status << ','
            << row.fmm.status << ',' << row.direct.panel_count << ','
            << row.fmm.panel_count << ',' << row.direct.iterations << ','
            << row.fmm.iterations << ',' << row.direct.final_residual << ','
            << row.fmm.final_residual << ',' << row.direct.true_residual_rel
            << ',' << row.fmm.true_residual_rel << ','
            << row.direct.energy_kj_mol << ',' << row.fmm.energy_kj_mol
            << ',' << row.energy_delta_kj_mol << ','
            << row.energy_delta_abs_kj_mol << ',' << row.energy_delta_rel
            << ',' << row.direct.matvec_ms_avg << ','
            << row.fmm.matvec_ms_avg << ',' << row.matvec_speedup << ','
            << row.fmm.fmm_m2l_mode << ',' << row.fmm.fmm_m2l_near_pairs
            << ',' << row.fmm.fmm_m2l_far_pairs << ','
            << row.fmm.fmm_m2l_setup_ms << ',' << row.fmm.fmm_m2l_apply_ms
            << ',' << row.direct_seconds << ',' << row.fmm_seconds << ','
            << row.gate << '\n';
}

void run_generated_smoke() {
  const SolveRow row =
      run_solve_pair("generated", make_generated_atoms(12), 12, 8, 120, 16,
                     1.0e-6, 4);
  require_true(row.direct.status == "ok", "generated direct solve ok");
  require_true(row.fmm.status == "ok", "generated FMM solve ok");
  require_true(row.fmm.fmm_m2l_mode == "single-matvec-analytic",
               "generated FMM telemetry must be bridged");
}

int run_fixture_probe(int argc, char **argv) {
  std::string fixture = "research-medium";
  std::string pqr_path;
  std::size_t atom_cap = 128;
  int surface_points = 8;
  int max_iter = 120;
  int restart = 24;
  double tolerance = 1.0e-6;
  int threads = 4;
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
    } else if (const std::string v = value("--atom-cap="); !v.empty()) {
      atom_cap = static_cast<std::size_t>(std::stoull(v));
    } else if (const std::string v = value("--surface-points="); !v.empty()) {
      surface_points = std::stoi(v);
    } else if (const std::string v = value("--max-iter="); !v.empty()) {
      max_iter = std::stoi(v);
    } else if (const std::string v = value("--restart="); !v.empty()) {
      restart = std::stoi(v);
    } else if (const std::string v = value("--tolerance="); !v.empty()) {
      tolerance = std::stod(v);
    } else if (const std::string v = value("--threads="); !v.empty()) {
      threads = std::stoi(v);
    } else {
      throw std::invalid_argument("unknown argument: " + arg);
    }
  }
  require_true(!pqr_path.empty(), "--pqr is required");
  const SolveRow row = run_solve_pair(fixture, read_pqr(pqr_path), atom_cap,
                                      surface_points, max_iter, restart,
                                      tolerance, threads);
  write_header();
  write_row(row);
  return row.gate == "pass" ? 0 : 2;
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc > 1 && std::string(argv[1]) == "--fixture-probe") {
      return run_fixture_probe(argc, argv);
    }
    run_generated_smoke();
    std::cout << "PASS titan_fmm_m2l_scope_l_probe\n";
    return 0;
  } catch (const std::exception &ex) {
    std::cerr << "FAIL titan_fmm_m2l_scope_l_probe: " << ex.what() << '\n';
    return 1;
  }
}
