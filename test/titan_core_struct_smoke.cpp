#include "MOL.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_core_struct_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-9) {
    std::cerr << "FAIL titan_core_struct_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

mol_info make_atom(std::size_t resnr, const std::string &resname,
                   const std::string &atom, bool sidechain) {
  mol_info info(resnr, atom, atom, "A", resname);
  info.set_chainnum(0);
  info.set_respos(resnr);
  info.set_internal_resnr(resnr);
  info.set_side_chain(sidechain);
  return info;
}

index make_chain_index(std::size_t first_res_pos,
                       std::size_t second_res_pos) {
  index chain(std::string("A"), std::size_t{0}, first_res_pos,
              std::size_t{1}, 0.0, 0.0, 0.0, std::string("ALA"));
  chain.add_base(second_res_pos, std::size_t{2}, 0.0, 0.0, 0.0,
                 std::string("GLY"));
  return chain;
}

void run_sidechain_com_unit_weight_fallback() {
  MOL frame;
  frame.pos.x = {2.0, 8.0, 10.0, 20.0};
  frame.pos.y = {0.0, 0.0, 0.0, 0.0};
  frame.pos.z = {0.0, 0.0, 0.0, 0.0};

  std::vector<mol_info> info{
      make_atom(1, "ALA", "CB", true), make_atom(1, "ALA", "CG", true),
      make_atom(2, "GLY", "N", false), make_atom(2, "GLY", "CA", false)};

  Multi_MOL_parameters params;
  params.resize(info.size());
  for (std::size_t i = 0; i < info.size(); ++i)
    params[i].set_mass(0.0);

  std::vector<index> indices;
  indices.push_back(make_chain_index(1, 3));

  frame.take_com(params, 0, 5, indices, info);

  require_close(indices[0].get_res_x_e(0), 5.0,
                "zero-mass side-chain COM uses unit weights");
  require_close(indices[0].get_res_x_e(1), 15.0,
                "zero-mass all-atom fallback uses unit weights");
}

void run_sidechain_com_keeps_mass_weighting() {
  MOL frame;
  frame.pos.x = {2.0, 8.0};
  frame.pos.y = {0.0, 0.0};
  frame.pos.z = {0.0, 0.0};

  std::vector<mol_info> info{make_atom(1, "ALA", "CB", true),
                             make_atom(1, "ALA", "CG", true)};
  Multi_MOL_parameters params;
  params.resize(info.size());
  params[0].set_mass(1.0);
  params[1].set_mass(3.0);

  index chain(std::string("A"), std::size_t{0}, std::size_t{1},
              std::size_t{1}, 0.0, 0.0, 0.0, std::string("ALA"));
  std::vector<index> indices{chain};

  frame.take_com(params, 0, 3, indices, info);
  require_close(indices[0].get_res_x_e(0), 6.5,
                "positive masses keep mass-weighted COM");
}

void run_all_atom_com_unit_weight_fallback() {
  MOL frame;
  frame.pos.x = {1.0, 7.0};
  frame.pos.y = {2.0, 6.0};
  frame.pos.z = {0.0, 0.0};

  std::vector<mol_info> info{make_atom(1, "ALA", "N", false),
                             make_atom(1, "ALA", "CA", false)};
  Multi_MOL_parameters params;
  params.resize(info.size());
  params[0].set_mass(0.0);
  params[1].set_mass(0.0);

  index chain(std::string("A"), std::size_t{0}, std::size_t{1},
              std::size_t{1}, 0.0, 0.0, 0.0, std::string("ALA"));
  std::vector<index> indices{chain};
  multi_points out;

  frame.take_com_all_atoms(params, 0, 3, indices, info, out);
  require(out.points.size() == 1, "all-atom COM emits one point");
  require_close(out.points[0].x, 4.0, "all-atom zero mass x geometric center");
  require_close(out.points[0].y, 4.0, "all-atom zero mass y geometric center");
  require_close(out.mass[0], 2.0, "all-atom zero mass records effective mass");
}

void run_mol_invariants() {
  MOL frame;
  frame.pos.x = {1.0, 2.0};
  frame.pos.y = {1.0};
  frame.pos.z = {1.0, 2.0};
  require(!frame.coordinates_consistent(),
          "coordinate invariant detects x/y/z size mismatch");

  frame.pos.y.push_back(2.0);
  require(frame.coordinates_consistent(),
          "coordinate invariant accepts matched x/y/z sizes");

  const std::size_t old_capacity = frame.pos.x.capacity();
  frame.clear_MOL();
  require(frame.pos.x.capacity() >= old_capacity,
          "clear_MOL keeps coordinate capacity for reuse");
  frame.release_memory();
  require(frame.pos.x.capacity() == 0,
          "release_memory explicitly returns coordinate storage");
}

void run_force_accessors_guard_release_builds() {
  MOL frame;
  frame.pos.x = {1.0, 2.0};
  frame.pos.y = {1.0, 2.0};
  frame.pos.z = {1.0, 2.0};

  require(!frame.force_enabled(), "force buffers start disabled");
  try {
    (void)frame.get_force_x();
    require(false, "get_force_x throws before force initialization");
  } catch (const std::logic_error &err) {
    require(std::string(err.what()) == "MOL force buffers are not initialized",
            "get_force_x emits stable runtime guard message");
  }

  frame.initialize_forces(2);
  require(frame.force_enabled(), "initialize_forces enables force buffers");
  require(frame.get_force_x().size() == frame.pos.x.size(),
          "force x span matches atom count");
  require(frame.get_thread_force_x(1).size() == frame.pos.x.size(),
          "thread-local force span matches atom count");
  frame.clear_MOL();
  require(!frame.force_enabled(), "clear_MOL disables force buffers");
}

} // namespace

int main() {
  run_sidechain_com_unit_weight_fallback();
  run_sidechain_com_keeps_mass_weighting();
  run_all_atom_com_unit_weight_fallback();
  run_mol_invariants();
  run_force_accessors_guard_release_builds();
  std::cout << "PASS titan_core_struct_smoke\n";
  return 0;
}
