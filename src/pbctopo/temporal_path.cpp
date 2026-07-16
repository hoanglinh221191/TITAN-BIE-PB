#include "temporal_path.h"

#include <algorithm>
#include <cmath>

namespace titan_pbctopo {

double pbctopo_temporal_center_jump2(const std::array<double, 3> &lhs,
                                     const std::array<double, 3> &rhs) {
  double jump2 = 0.0;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    double delta = rhs[axis] - lhs[axis];
    delta -= std::round(delta);
    jump2 += delta * delta;
  }
  return jump2;
}

double pbctopo_temporal_assignment_jump2(
    const PbctopoTemporalCandidate &lhs,
    const PbctopoTemporalCandidate &rhs) {
  if (lhs.certified && rhs.certified && !lhs.component_offsets.empty() &&
      lhs.component_offsets.size() == rhs.component_offsets.size()) {
    double jump2 = 0.0;
    for (std::size_t idx = 0; idx < lhs.component_offsets.size(); ++idx) {
      const Int3 delta = rhs.component_offsets[idx] -
                         lhs.component_offsets[idx];
      jump2 += static_cast<double>(delta.x * delta.x + delta.y * delta.y +
                                   delta.z * delta.z);
    }
    return jump2 / static_cast<double>(lhs.component_offsets.size());
  }
  return pbctopo_temporal_center_jump2(lhs.center, rhs.center);
}

PbctopoTemporalPathResult solve_pbctopo_temporal_viterbi(
    const std::vector<std::vector<PbctopoTemporalCandidate>> &frames,
    double transition_penalty) {
  PbctopoTemporalPathResult result;
  const std::size_t invalid = std::numeric_limits<std::size_t>::max();
  result.selected.assign(frames.size(), invalid);
  transition_penalty = std::max(0.0, transition_penalty);

  std::size_t segment_begin = 0;
  while (segment_begin < frames.size()) {
    while (segment_begin < frames.size() && frames[segment_begin].empty())
      ++segment_begin;
    if (segment_begin == frames.size())
      break;
    std::size_t segment_end = segment_begin;
    while (segment_end + 1 < frames.size() &&
           !frames[segment_end + 1].empty())
      ++segment_end;

    std::vector<double> previous(frames[segment_begin].size(), 0.0);
    for (std::size_t state = 0; state < previous.size(); ++state)
      previous[state] = frames[segment_begin][state].emission_cost;
    std::vector<std::vector<std::size_t>> backtrace(segment_end - segment_begin + 1);

    for (std::size_t frame = segment_begin + 1; frame <= segment_end; ++frame) {
      std::vector<double> current(frames[frame].size(),
                                  std::numeric_limits<double>::infinity());
      backtrace[frame - segment_begin].assign(frames[frame].size(), 0);
      for (std::size_t state = 0; state < frames[frame].size(); ++state) {
        for (std::size_t prev = 0; prev < frames[frame - 1].size(); ++prev) {
          const double cost =
              previous[prev] + frames[frame][state].emission_cost +
              transition_penalty * pbctopo_temporal_assignment_jump2(
                                       frames[frame - 1][prev],
                                       frames[frame][state]);
          if (cost < current[state] ||
              (cost == current[state] &&
               prev < backtrace[frame - segment_begin][state])) {
            current[state] = cost;
            backtrace[frame - segment_begin][state] = prev;
          }
        }
      }
      previous = std::move(current);
    }

    const auto best_it = std::min_element(previous.begin(), previous.end());
    std::size_t state = static_cast<std::size_t>(best_it - previous.begin());
    result.total_cost += *best_it;
    for (std::size_t frame = segment_end;; --frame) {
      result.selected[frame] = state;
      if (frame == segment_begin)
        break;
      const std::size_t prior = backtrace[frame - segment_begin][state];
      if (pbctopo_temporal_assignment_jump2(
              frames[frame - 1][prior], frames[frame][state]) > 1.0e-18)
        ++result.transition_count;
      state = prior;
    }
    segment_begin = segment_end + 1;
  }
  result.valid = std::any_of(result.selected.begin(), result.selected.end(),
                             [&](std::size_t value) { return value != invalid; });
  return result;
}

} // namespace titan_pbctopo
