#ifndef TITAN_ENERGY_VALIDATION_H
#define TITAN_ENERGY_VALIDATION_H

#include "MOL.h"

#include <cstddef>
#include <optional>
#include <string_view>

class TitanExecutor;

namespace titan_energy_validation {

inline constexpr double kMinimumPairDistanceNm = 1.0e-10;
inline constexpr std::string_view kOverlappingAtomsError =
    "invalid_overlapping_atoms";

struct OverlappingAtomPair {
  std::size_t mol1_atom{};
  std::size_t mol2_atom{};
  double distance_nm{};
};

[[nodiscard]] std::optional<OverlappingAtomPair>
find_overlapping_atoms(const MOL &mol1, const MOL &mol2,
                       bool identical_selections, TitanExecutor &executor);

} // namespace titan_energy_validation

#endif
