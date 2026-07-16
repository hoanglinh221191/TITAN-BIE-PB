#include "certification.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace titan_pbctopo {

PbctopoImageEquationResiduals evaluate_pbctopo_image_equations(
    std::span<const PbctopoAtom> wrapped_atoms,
    std::span<const Vec3> placed_points,
    std::span<const PbctopoLocalEdge> local_edges,
    const PbctopoLatticeMetric &lattice, double image_tolerance) {
  PbctopoImageEquationResiduals result;
  if (wrapped_atoms.empty() || placed_points.size() != wrapped_atoms.size() ||
      !std::isfinite(image_tolerance) || image_tolerance <= 0.0) {
    result.atom_image_residuals = 1;
    return result;
  }

  std::vector<Int3> atom_images(wrapped_atoms.size());
  std::vector<unsigned char> atom_image_valid(wrapped_atoms.size(), 0);
  for (std::size_t atom = 0; atom < wrapped_atoms.size(); ++atom) {
    const auto &wrapped = wrapped_atoms[atom];
    const auto &placed = placed_points[atom];
    if (!std::isfinite(placed.x) || !std::isfinite(placed.y) ||
        !std::isfinite(placed.z)) {
      ++result.atom_image_residuals;
      continue;
    }
    const auto scaled = lattice.cartesian_to_scaled(
        placed.x - wrapped.x, placed.y - wrapped.y, placed.z - wrapped.z);
    if (!std::isfinite(scaled[0]) || !std::isfinite(scaled[1]) ||
        !std::isfinite(scaled[2])) {
      ++result.atom_image_residuals;
      continue;
    }
    const Int3 image{static_cast<int>(std::llround(scaled[0])),
                     static_cast<int>(std::llround(scaled[1])),
                     static_cast<int>(std::llround(scaled[2]))};
    const double error = std::max(
        {std::fabs(scaled[0] - static_cast<double>(image.x)),
         std::fabs(scaled[1] - static_cast<double>(image.y)),
         std::fabs(scaled[2] - static_cast<double>(image.z))});
    result.max_atom_image_error = std::max(result.max_atom_image_error, error);
    if (error > image_tolerance) {
      ++result.atom_image_residuals;
      continue;
    }
    atom_images[atom] = image;
    atom_image_valid[atom] = 1;
  }

  for (const auto &edge : local_edges) {
    if (edge.kind != PbctopoLocalEdgeKind::TopologyBond)
      continue;
    ++result.topology_edges_checked;
    if (edge.a >= wrapped_atoms.size() || edge.b >= wrapped_atoms.size() ||
        edge.a == edge.b || atom_image_valid[edge.a] == 0 ||
        atom_image_valid[edge.b] == 0) {
      ++result.topology_edge_residuals;
      continue;
    }
    const auto &a = wrapped_atoms[edge.a];
    const auto &b = wrapped_atoms[edge.b];
    const auto nearest = lattice.nearest_image_scaled(
        b.sx - a.sx, b.sy - a.sy, b.sz - a.sz);
    if (!nearest.valid) {
      ++result.topology_edge_residuals;
      continue;
    }
    const Int3 expected{nearest.image_x, nearest.image_y, nearest.image_z};
    const Int3 actual = atom_images[edge.b] - atom_images[edge.a];
    if (!same_delta(actual, expected))
      ++result.topology_edge_residuals;
  }
  result.input_valid = true;
  return result;
}

PbctopoRequiredContactResiduals evaluate_pbctopo_required_contacts(
    std::span<const PbctopoAtom> atoms,
    std::span<const Vec3> placed_points,
    std::span<const PbctopoLocalEdge> local_edges) {
  PbctopoRequiredContactResiduals result;
  if (atoms.empty() || placed_points.size() != atoms.size())
    return result;

  for (const auto &edge : local_edges) {
    if (edge.kind != PbctopoLocalEdgeKind::TopologyBond)
      continue;
    if (edge.a >= atoms.size() || edge.b >= atoms.size() || edge.a == edge.b)
      continue;
    if (atoms[edge.a].owner == atoms[edge.b].owner)
      continue;

    ++result.contacts_checked;
    const Vec3 delta = placed_points[edge.a] - placed_points[edge.b];
    const double d2 = delta.x * delta.x + delta.y * delta.y +
                      delta.z * delta.z;
    if (!std::isfinite(d2) || !std::isfinite(edge.break_cutoff_d2) ||
        d2 > edge.break_cutoff_d2 + 1.0e-6) {
      ++result.contacts_lost;
    }
  }
  result.input_valid = true;
  return result;
}

PbctopoHardFeasibility evaluate_pbctopo_hard_feasibility(
    const PbctopoCutCandidate &candidate,
    const PbctopoImageEquationResiduals &image_equations,
    std::size_t cycle_mismatches,
    PbctopoRequiredContactResiduals required_contacts) {
  PbctopoHardFeasibility hard;
  hard.evaluated = candidate.valid;
  hard.assignment_available = candidate.assignment.valid &&
                              image_equations.input_valid;
  hard.topology_edges_checked = image_equations.topology_edges_checked;
  hard.hard_edge_residuals = image_equations.topology_edge_residuals;
  hard.atom_image_residuals = image_equations.atom_image_residuals;
  hard.max_atom_image_error = image_equations.max_atom_image_error;
  hard.cycle_mismatches = cycle_mismatches;
  hard.hard_clashes = candidate.ultra_clash ? 1 : 0;
  hard.required_contacts_checked = required_contacts.contacts_checked;
  hard.required_contacts_lost = required_contacts.contacts_lost;
  return hard;
}

bool certify_pbctopo_candidate(PbctopoCutCandidate &candidate,
                               const PbctopoImageEquationResiduals &image_equations,
                               std::size_t cycle_mismatches,
                               PbctopoRequiredContactResiduals required_contacts) {
  candidate.hard_feasibility = evaluate_pbctopo_hard_feasibility(
      candidate, image_equations, cycle_mismatches, required_contacts);
  return candidate.hard_feasibility.certified();
}

} // namespace titan_pbctopo
