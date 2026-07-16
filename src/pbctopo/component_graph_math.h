#ifndef TITAN_PBCTOPO_COMPONENT_GRAPH_MATH_H
#define TITAN_PBCTOPO_COMPONENT_GRAPH_MATH_H

#include "types.h"

#include <cstddef>
#include <vector>

namespace titan_pbctopo {

void finalize_normalized_pair_shift_score(PbctopoPairShiftScore &score,
                                          std::size_t lhs_atom_count,
                                          std::size_t rhs_atom_count,
                                          double contact_cutoff);

PbctopoCycleSyncMetrics compute_fundamental_cycle_residuals(
    std::size_t component_count,
    const std::vector<PbctopoPairShiftCandidateSet> &pair_candidates);

} // namespace titan_pbctopo

#endif
