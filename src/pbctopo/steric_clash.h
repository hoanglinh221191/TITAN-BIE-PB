#ifndef TITAN_PBCTOPO_STERIC_CLASH_H
#define TITAN_PBCTOPO_STERIC_CLASH_H

#include "types.h"

#include <vector>

namespace titan_pbctopo {

struct PbctopoPairClashThreshold {
  double cutoff = 0.0;
  double cutoff2 = 0.0;
  bool used_vdw_radii = false;
};

PbctopoPairClashThreshold pbctopo_pair_clash_threshold(
    const PbctopoAtom &lhs, const PbctopoAtom &rhs, double fallback_cutoff);

void prepare_pbctopo_steric_radii(std::vector<PbctopoAtom> &atoms,
                                  double vdw_clash_scale);

double pbctopo_max_clash_cutoff(const std::vector<PbctopoAtom> &atoms,
                                double fallback_cutoff);

bool pbctopo_atoms_clash(double d2, const PbctopoAtom &lhs,
                         const PbctopoAtom &rhs, double fallback_cutoff,
                         PbctopoTimingStats *timing = nullptr);

} // namespace titan_pbctopo

#endif
