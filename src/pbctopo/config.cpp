#include "pbctopo/config.h"

#include "arguments.h"

#include <string>

namespace titan_pbctopo {

RuntimeCutoffs resolve_runtime_cutoffs(const arguments &Arguments) {
  RuntimeCutoffs cutoffs;
  if (Arguments.force_f == "MARTINI") {
    cutoffs.bond_cutoff = 5.5;
    cutoffs.clash_cutoff = 3.5;
    cutoffs.pair_contact_cutoff = 8.0;
    cutoffs.continuity_cutoff = 20.0;
  }
  if (Arguments.pbc.pbctopo_bond_cutoff_set)
    cutoffs.bond_cutoff = Arguments.pbc.pbctopo_bond_cutoff;
  if (Arguments.pbc.pbctopo_clash_cutoff_set)
    cutoffs.clash_cutoff = Arguments.pbc.pbctopo_clash_cutoff;
  if (Arguments.pbc.pbctopo_contact_cutoff_set)
    cutoffs.pair_contact_cutoff = Arguments.pbc.pbctopo_contact_cutoff;
  if (Arguments.pbc.pbctopo_continuity_cutoff_set)
    cutoffs.continuity_cutoff = Arguments.pbc.pbctopo_continuity_cutoff;
  return cutoffs;
}

SolveMode resolve_solve_mode(const arguments &Arguments) {
  if (Arguments.pbc.pbctopo_mode == "llps-local")
    return SolveMode::LlpsLocal;
  return SolveMode::Vibe;
}

const char *solve_mode_name(SolveMode mode) {
  switch (mode) {
  case SolveMode::LlpsLocal:
    return "llps-local";
  case SolveMode::Vibe:
  default:
    return "vibe";
  }
}

AnchorBackend resolve_anchor_backend(const arguments &Arguments) {
  const std::string &backend = Arguments.pbc.pbctopo_anchor_backend;
  if (backend == "prim")
    return AnchorBackend::Prim;
  if (backend == "spatial-all")
    return AnchorBackend::SpatialAll;
  return AnchorBackend::Auto;
}

const char *anchor_backend_name(AnchorBackend backend) {
  switch (backend) {
  case AnchorBackend::Prim:
    return "prim";
  case AnchorBackend::SpatialAll:
    return "spatial-all";
  case AnchorBackend::Auto:
  default:
    return "auto";
  }
}

} // namespace titan_pbctopo
