#include "image_assignment.h"

#include <algorithm>
#include <cstdint>

namespace titan_pbctopo {
namespace {

constexpr std::uint64_t kFnvOffset = 1469598103934665603ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;

void hash_word(std::uint64_t &hash, std::uint32_t word) {
  for (unsigned shift = 0; shift < 32; shift += 8) {
    hash ^= static_cast<std::uint8_t>(word >> shift);
    hash *= kFnvPrime;
  }
}

} // namespace

PbctopoImageAssignment make_pbctopo_image_assignment(
    std::span<const Int3> component_offsets,
    std::span<const std::vector<std::size_t>> components) {
  PbctopoImageAssignment assignment;
  if (component_offsets.empty() ||
      component_offsets.size() != components.size()) {
    return assignment;
  }

  std::size_t root = 0;
  for (std::size_t idx = 1; idx < components.size(); ++idx) {
    if (components[idx].size() > components[root].size())
      root = idx;
  }

  assignment.global_gauge = component_offsets[root];
  assignment.component_offsets.reserve(component_offsets.size());
  std::uint64_t signature = kFnvOffset;
  hash_word(signature, static_cast<std::uint32_t>(component_offsets.size()));
  for (std::size_t idx = 0; idx < component_offsets.size(); ++idx) {
    const Int3 canonical = component_offsets[idx] - assignment.global_gauge;
    assignment.component_offsets.push_back(canonical);
    hash_word(signature, static_cast<std::uint32_t>(components[idx].size()));
    hash_word(signature, static_cast<std::uint32_t>(canonical.x));
    hash_word(signature, static_cast<std::uint32_t>(canonical.y));
    hash_word(signature, static_cast<std::uint32_t>(canonical.z));
  }
  assignment.signature = signature == 0 ? kFnvOffset : signature;
  assignment.valid = true;
  return assignment;
}

bool same_pbctopo_image_assignment(const PbctopoImageAssignment &lhs,
                                    const PbctopoImageAssignment &rhs) {
  if (!lhs.valid || !rhs.valid || lhs.signature != rhs.signature ||
      lhs.component_offsets.size() != rhs.component_offsets.size()) {
    return false;
  }
  for (std::size_t idx = 0; idx < lhs.component_offsets.size(); ++idx) {
    if (!same_delta(lhs.component_offsets[idx], rhs.component_offsets[idx]))
      return false;
  }
  return true;
}

} // namespace titan_pbctopo
