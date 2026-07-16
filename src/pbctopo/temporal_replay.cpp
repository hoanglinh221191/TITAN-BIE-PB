#include "temporal_replay.h"

#include "lattice_math.h"

namespace titan_pbctopo {

bool replay_pbctopo_assignment_delta(
    std::span<const Vec3> applied_points,
    std::span<const std::uint32_t> atom_component_index,
    const PbctopoImageAssignment &applied_assignment,
    const PbctopoImageAssignment &selected_assignment, box Box,
    std::vector<Vec3> &replayed_points) {
  replayed_points.clear();
  if (!applied_assignment.valid || !selected_assignment.valid ||
      applied_points.empty() ||
      applied_points.size() != atom_component_index.size() ||
      applied_assignment.component_offsets.size() !=
          selected_assignment.component_offsets.size()) {
    return false;
  }

  const std::size_t component_count =
      applied_assignment.component_offsets.size();
  std::vector<Vec3> translations(component_count);
  const PbctopoLatticeMetric lattice(Box);
  for (std::size_t comp_idx = 0; comp_idx < component_count; ++comp_idx) {
    const Int3 delta = selected_assignment.component_offsets[comp_idx] -
                       applied_assignment.component_offsets[comp_idx];
    const auto cart = lattice.scaled_to_cartesian(
        static_cast<double>(delta.x), static_cast<double>(delta.y),
        static_cast<double>(delta.z));
    translations[comp_idx] = {cart[0], cart[1], cart[2]};
  }

  replayed_points.resize(applied_points.size());
  for (std::size_t atom_idx = 0; atom_idx < applied_points.size(); ++atom_idx) {
    const std::size_t comp_idx = atom_component_index[atom_idx];
    if (comp_idx >= translations.size()) {
      replayed_points.clear();
      return false;
    }
    replayed_points[atom_idx] = {
        applied_points[atom_idx].x + translations[comp_idx].x,
        applied_points[atom_idx].y + translations[comp_idx].y,
        applied_points[atom_idx].z + translations[comp_idx].z};
  }
  return true;
}

} // namespace titan_pbctopo
