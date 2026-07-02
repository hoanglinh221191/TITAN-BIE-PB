#ifndef TITAN_PBCTOPO_COMPONENT_SEARCH_H
#define TITAN_PBCTOPO_COMPONENT_SEARCH_H

#include "config.h"
#include "types.h"
#include "../box_type.h"

#include <cstddef>
#include <vector>

namespace titan_pbctopo {

struct PbctopoComponentSearchCallbacks {
  using BuildSeedFn = PbctopoComponentCandidate (*)(
      const PbctopoCutCandidate &, const std::vector<std::vector<std::size_t>> &,
      const std::vector<Vec3> &);
  using ScoreAssignmentFn = PbctopoComponentCandidate (*)(
      const std::vector<PbctopoAtom> &, const std::vector<Vec3> &,
      const std::vector<std::vector<std::size_t>> &, double, std::size_t,
      std::size_t, const std::vector<PbctopoContinuityPair> &,
      const std::vector<std::size_t> &, double, double,
      const PbctopoComponentCandidate *, PbctopoTimingStats *, AnchorBackend,
      const PbctopoAnchorGraphContext *);
  using BetterCandidateFn = bool (*)(const PbctopoComponentCandidate &,
                                     const PbctopoComponentCandidate &);

  BuildSeedFn build_seed = nullptr;
  ScoreAssignmentFn score_assignment = nullptr;
  BetterCandidateFn better_candidate = nullptr;

  bool ready() const {
    return build_seed != nullptr && score_assignment != nullptr &&
           better_candidate != nullptr;
  }
};

PbctopoComponentCandidate search_component_shift_candidate(
    const std::vector<PbctopoAtom> &atoms, box Box,
    const std::vector<std::vector<std::size_t>> &components,
    const std::vector<PbctopoShift> &shifts,
    const std::vector<PbctopoContinuityPair> &continuity_pairs,
    const std::vector<std::size_t> &anchor_indices,
    const RuntimeCutoffs &cutoffs, std::vector<Vec3> &best_points_out,
    const PbctopoComponentSearchCallbacks &callbacks,
    const std::vector<Int3> *warm_component_offsets = nullptr,
    const PbctopoCutCandidate *incumbent_cut = nullptr,
    const std::vector<Vec3> *incumbent_points = nullptr,
    PbctopoTimingStats *timing = nullptr,
    AnchorBackend anchor_backend = AnchorBackend::Auto,
    const PbctopoAnchorGraphContext *anchor_context = nullptr);

PbctopoComponentCandidate search_component_shift_candidate_beam(
    const std::vector<PbctopoAtom> &atoms, box Box,
    const std::vector<std::vector<std::size_t>> &components,
    const std::vector<PbctopoShift> &shifts,
    const std::vector<PbctopoContinuityPair> &continuity_pairs,
    const std::vector<std::size_t> &anchor_indices,
    const RuntimeCutoffs &cutoffs, std::vector<Vec3> &best_points_out,
    std::size_t beam_width, const PbctopoComponentSearchCallbacks &callbacks,
    const std::vector<Int3> *warm_component_offsets = nullptr,
    const PbctopoCutCandidate *incumbent_cut = nullptr,
    const std::vector<Vec3> *incumbent_points = nullptr,
    PbctopoTimingStats *timing = nullptr,
    AnchorBackend anchor_backend = AnchorBackend::Auto,
    const PbctopoAnchorGraphContext *anchor_context = nullptr);

} // namespace titan_pbctopo

#endif
