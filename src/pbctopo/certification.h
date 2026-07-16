#ifndef TITAN_PBCTOPO_CERTIFICATION_H
#define TITAN_PBCTOPO_CERTIFICATION_H

#include "lattice_math.h"
#include "types.h"

#include <cstddef>
#include <span>

namespace titan_pbctopo {

struct PbctopoImageEquationResiduals {
  std::size_t topology_edges_checked = 0;
  std::size_t topology_edge_residuals = 0;
  std::size_t atom_image_residuals = 0;
  double max_atom_image_error = 0.0;
  bool input_valid = false;
};

struct PbctopoRequiredContactResiduals {
  std::size_t contacts_checked = 0;
  std::size_t contacts_lost = 0;
  bool input_valid = false;
};

PbctopoImageEquationResiduals evaluate_pbctopo_image_equations(
    std::span<const PbctopoAtom> wrapped_atoms,
    std::span<const Vec3> placed_points,
    std::span<const PbctopoLocalEdge> local_edges,
    const PbctopoLatticeMetric &lattice,
    double image_tolerance = 1.0e-6);

PbctopoRequiredContactResiduals evaluate_pbctopo_required_contacts(
    std::span<const PbctopoAtom> atoms,
    std::span<const Vec3> placed_points,
    std::span<const PbctopoLocalEdge> local_edges);

PbctopoHardFeasibility evaluate_pbctopo_hard_feasibility(
    const PbctopoCutCandidate &candidate,
    const PbctopoImageEquationResiduals &image_equations,
    std::size_t cycle_mismatches = 0,
    PbctopoRequiredContactResiduals required_contacts = {});

bool certify_pbctopo_candidate(PbctopoCutCandidate &candidate,
                               const PbctopoImageEquationResiduals &image_equations,
                               std::size_t cycle_mismatches = 0,
                               PbctopoRequiredContactResiduals required_contacts = {});

} // namespace titan_pbctopo

#endif
