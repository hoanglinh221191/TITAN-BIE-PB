#include "titan_contact_geometry.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

mol_info make_atom(std::size_t resnr, const std::string &resname,
                   const std::string &atom, bool sidechain) {
  mol_info info(resnr, atom, atom, "A", resname);
  info.set_chainnum(0);
  info.set_respos(resnr);
  info.set_internal_resnr(resnr);
  info.set_side_chain(sidechain);
  return info;
}

Multi_MOL make_molecule(const std::vector<mol_info> &atoms,
                        const std::vector<double> &xs,
                        const std::vector<double> &masses) {
  Multi_MOL mol;
  mol.info = atoms;
  mol.parameters.resize(atoms.size());
  for (std::size_t i = 0; i < masses.size(); ++i)
    mol.parameters[i].set_mass(masses[i]);

  MOL frame;
  frame.pos.x.assign(xs.begin(), xs.end());
  frame.pos.y.assign(xs.size(), 0.0);
  frame.pos.z.assign(xs.size(), 0.0);
  mol.m_mol.push_back(std::move(frame));
  return mol;
}

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_contact_geometry_cache_smoke: " << message
              << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > 1.0e-9) {
    std::cerr << "FAIL titan_contact_geometry_cache_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

std::vector<titan_contact_geometry::ContactPoint>
reference_sidechain_com_points(const Multi_MOL &mol) {
  std::vector<titan_contact_geometry::ContactPoint> points;
  require(!mol.m_mol.empty(), "reference molecule has a frame");
  const MOL &frame = mol.m_mol.front();
  for (std::size_t atom_idx = 0; atom_idx < mol.info.size(); ++atom_idx) {
    const mol_info &atom_info = mol.info[atom_idx];
    if (atom_idx > 0 &&
        titan_contact::same_residue(points.back().residue,
                                    titan_contact::residue_key_from_atom(
                                        atom_info, 0, 1))) {
      continue;
    }
    const titan_contact::ResidueKey residue =
        titan_contact::residue_key_from_atom(atom_info, 0, 1);
    std::vector<std::size_t> all_atoms;
    std::vector<std::size_t> sidechain_atoms;
    for (std::size_t scan_idx = 0; scan_idx < mol.info.size(); ++scan_idx) {
      const titan_contact::ResidueKey scan_residue =
          titan_contact::residue_key_from_atom(mol.info[scan_idx], 0, 1);
      if (!(scan_residue == residue))
        continue;
      all_atoms.push_back(scan_idx);
      if (mol.info[scan_idx].is_side_chain())
        sidechain_atoms.push_back(scan_idx);
    }
    const auto &indices =
        sidechain_atoms.empty() ? all_atoms : sidechain_atoms;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double total_mass = 0.0;
    double unit_x = 0.0;
    double unit_y = 0.0;
    double unit_z = 0.0;
    for (const std::size_t idx : indices) {
      unit_x += frame.pos.x[idx];
      unit_y += frame.pos.y[idx];
      unit_z += frame.pos.z[idx];
      const double mass = mol.parameters.mass[idx];
      if (mass > 0.0) {
        x += frame.pos.x[idx] * mass;
        y += frame.pos.y[idx] * mass;
        z += frame.pos.z[idx] * mass;
        total_mass += mass;
      }
    }
    if (!(total_mass > 0.0)) {
      const double denom = static_cast<double>(indices.size());
      x = unit_x / denom;
      y = unit_y / denom;
      z = unit_z / denom;
      total_mass = 1.0;
    }
    points.push_back(titan_contact_geometry::ContactPoint{
        residue, x / total_mass, y / total_mass, z / total_mass});
  }
  return points;
}

void run_mass_weighted_sidechain_com() {
  Multi_MOL mol = make_molecule(
      {make_atom(1, "ALA", "N", false), make_atom(1, "ALA", "CB", true),
       make_atom(1, "ALA", "CG", true)},
      {0.0, 2.0, 8.0}, {9.0, 1.0, 3.0});

  titan_contact_geometry::ContactGeometryCache cache;
  require(cache.build(mol, 1, titan_contact::ContactMode::SidechainCOM),
          "side-chain COM cache build succeeds");
  std::vector<titan_contact_geometry::ContactPoint> points;
  require(cache.build_points(mol, 0, points), "side-chain COM points build");
  require(points.size() == 1, "one residue point emitted");
  require_close(points[0].x, 6.5, "COM uses positive side-chain masses");
}

void run_zero_mass_unit_weight_fallback() {
  Multi_MOL mol = make_molecule(
      {make_atom(1, "ALA", "CB", true), make_atom(1, "ALA", "CG", true)},
      {2.0, 8.0}, {0.0, 0.0});

  titan_contact_geometry::ContactGeometryCache cache;
  require(cache.build(mol, 1, titan_contact::ContactMode::SidechainCOM),
          "zero-mass cache build succeeds");
  std::vector<titan_contact_geometry::ContactPoint> points;
  require(cache.build_points(mol, 0, points), "zero-mass points build");
  require_close(points[0].x, 5.0, "zero-mass COM uses unit weights");
}

void run_glycine_all_atom_fallback() {
  Multi_MOL mol = make_molecule(
      {make_atom(2, "GLY", "N", false), make_atom(2, "GLY", "CA", false)},
      {10.0, 20.0}, {0.0, 0.0});

  titan_contact_geometry::ContactGeometryCache cache;
  require(cache.build(mol, 1, titan_contact::ContactMode::SidechainCOM),
          "glycine cache build succeeds");
  std::vector<titan_contact_geometry::ContactPoint> points;
  require(cache.build_points(mol, 0, points), "glycine points build");
  require_close(points[0].x, 15.0,
                "glycine/no-sidechain falls back to all atoms");
}

void run_atom_mode_selection() {
  Multi_MOL mol = make_molecule(
      {make_atom(1, "ALA", "N", false), make_atom(1, "ALA", "CA", false),
       make_atom(1, "ALA", "CB", true)},
      {1.0, 2.0, 3.0}, {1.0, 1.0, 1.0});

  titan_contact_geometry::ContactGeometryCache cache;
  require(cache.build(mol, 1, titan_contact::ContactMode::CA),
          "CA atom-mode cache build succeeds");
  std::vector<titan_contact_geometry::ContactPoint> points;
  require(cache.build_points(mol, 0, points), "CA atom-mode points build");
  require(points.size() == 1, "CA atom mode emits only selected atom");
  require_close(points[0].x, 2.0, "CA atom mode uses CA coordinate");
}

void run_invalid_metadata_fails() {
  Multi_MOL mol = make_molecule({make_atom(1, "ALA", "CB", true)}, {1.0},
                                {1.0});
  mol.parameters.resize(0);

  titan_contact_geometry::ContactGeometryCache cache;
  require(!cache.build(mol, 1, titan_contact::ContactMode::SidechainCOM),
          "side-chain COM cache rejects missing mass metadata");
}

void run_optional_benchmark() {
  if (std::getenv("TITAN_CONTACT_GEOMETRY_BENCH") == nullptr)
    return;

  const std::size_t residues = 1000;
  std::vector<mol_info> atoms;
  std::vector<double> xs;
  std::vector<double> masses;
  atoms.reserve(residues * 2);
  xs.reserve(residues * 2);
  masses.reserve(residues * 2);
  for (std::size_t i = 0; i < residues; ++i) {
    atoms.push_back(make_atom(i + 1, "ALA", "CB", true));
    atoms.push_back(make_atom(i + 1, "ALA", "CG", true));
    xs.push_back(static_cast<double>(i % 100));
    xs.push_back(static_cast<double>(i % 100) + 2.0);
    masses.push_back(1.0);
    masses.push_back(3.0);
  }
  Multi_MOL mol = make_molecule(atoms, xs, masses);

  const auto ref_start = std::chrono::steady_clock::now();
  auto reference = reference_sidechain_com_points(mol);
  const auto ref_end = std::chrono::steady_clock::now();

  titan_contact_geometry::ContactGeometryCache cache;
  require(cache.build(mol, 1, titan_contact::ContactMode::SidechainCOM),
          "benchmark cache build succeeds");
  std::vector<titan_contact_geometry::ContactPoint> cached;
  const auto cached_start = std::chrono::steady_clock::now();
  require(cache.build_points(mol, 0, cached), "benchmark cached build points");
  const auto cached_end = std::chrono::steady_clock::now();

  require(reference.size() == cached.size(), "benchmark point counts match");
  for (std::size_t i = 0; i < reference.size(); ++i)
    require_close(cached[i].x, reference[i].x,
                  "benchmark cached point matches reference");

  const auto ref_us = std::chrono::duration_cast<std::chrono::microseconds>(
                          ref_end - ref_start)
                          .count();
  const auto cached_us =
      std::chrono::duration_cast<std::chrono::microseconds>(cached_end -
                                                            cached_start)
          .count();
  std::cout << "BENCH titan_contact_geometry reference_us=" << ref_us
            << " cached_us=" << cached_us << " residues=" << residues << '\n';
}

} // namespace

int main() {
  run_mass_weighted_sidechain_com();
  run_zero_mass_unit_weight_fallback();
  run_glycine_all_atom_fallback();
  run_atom_mode_selection();
  run_invalid_metadata_fails();
  run_optional_benchmark();
  std::cout << "PASS titan_contact_geometry_cache_smoke\n";
  return 0;
}
