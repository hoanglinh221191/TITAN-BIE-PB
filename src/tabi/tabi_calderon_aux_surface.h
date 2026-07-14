#ifndef TITAN_PBE_TABI_CALDERON_AUX_SURFACE_H
#define TITAN_PBE_TABI_CALDERON_AUX_SURFACE_H

#include "tabi_calderon_trace_space.h"
#include "tabi_types.h"

#include <span>
#include <string>

class TitanExecutor;

namespace TitanPBE::Tabi::calderon {

// A closed implicit surface mesh and its topology audit. The same deterministic
// extractor serves the auxiliary Calderon trace and the opt-in primary SES.
struct AuxiliarySurfaceResult {
  SurfaceMesh mesh;
  TraceSpaces trace_spaces;
  std::string status{"not_built"};
  double build_ms{0.0};
  std::size_t grid_vertices{0};
  std::size_t odd_ray_rows{0};
};

[[nodiscard]] AuxiliarySurfaceResult build_implicit_vdw_auxiliary_surface(
    std::span<const Atom> atoms, double grid_spacing,
    TitanExecutor *executor = nullptr);

[[nodiscard]] AuxiliarySurfaceResult build_implicit_ses_surface(
    std::span<const Atom> atoms, double grid_spacing, double probe_radius,
    TitanExecutor *executor = nullptr,
    std::size_t max_grid_vertices = 4000000,
    bool validate_topology = true);

[[nodiscard]] AuxiliarySurfaceResult build_analytic_ses_surface(
    std::span<const Atom> atoms, double grid_spacing, double probe_radius,
    TitanExecutor *executor = nullptr,
    std::size_t max_grid_vertices = 4000000,
    bool validate_topology = true);

} // namespace TitanPBE::Tabi::calderon

#endif // TITAN_PBE_TABI_CALDERON_AUX_SURFACE_H
