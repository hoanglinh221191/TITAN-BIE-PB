#ifndef TITAN_PBCTOPO_TEMPORAL_PATH_H
#define TITAN_PBCTOPO_TEMPORAL_PATH_H

#include "types.h"

#include <array>
#include <cstddef>
#include <limits>
#include <vector>

namespace titan_pbctopo {

struct PbctopoTemporalCandidate {
  std::array<double, 3> center{1.5, 1.5, 1.5};
  double emission_cost = 0.0;
  std::vector<Int3> component_offsets;
  std::uint64_t assignment_signature = 0;
  bool certified = false;
};

struct PbctopoTemporalPathResult {
  bool valid = false;
  std::vector<std::size_t> selected;
  double total_cost = 0.0;
  std::size_t transition_count = 0;
};

double pbctopo_temporal_center_jump2(const std::array<double, 3> &lhs,
                                     const std::array<double, 3> &rhs);
double pbctopo_temporal_assignment_jump2(
    const PbctopoTemporalCandidate &lhs,
    const PbctopoTemporalCandidate &rhs);

PbctopoTemporalPathResult solve_pbctopo_temporal_viterbi(
    const std::vector<std::vector<PbctopoTemporalCandidate>> &frames,
    double transition_penalty);

} // namespace titan_pbctopo

#endif
