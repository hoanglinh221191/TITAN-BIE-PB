#ifndef TITAN_PBCTOPO_COMPONENT_SEARCH_PAIR_TABLE_H
#define TITAN_PBCTOPO_COMPONENT_SEARCH_PAIR_TABLE_H

#include "types.h"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace titan_pbctopo {

class PbctopoShiftMask {
public:
  PbctopoShiftMask() = default;
  explicit PbctopoShiftMask(std::size_t shift_count);

  static PbctopoShiftMask full(std::size_t shift_count);
  static PbctopoShiftMask only(std::size_t shift_count, std::size_t shift_idx);

  void reset(std::size_t shift_count);
  void set_all();
  void clear_all();
  void set_only(std::size_t shift_idx);
  void set(std::size_t shift_idx);
  void clear(std::size_t shift_idx);
  bool test(std::size_t shift_idx) const;
  bool empty() const;
  std::size_t count() const;
  bool intersect_with(const PbctopoShiftMask &other);
  std::size_t shift_count() const { return shift_count_; }
  std::size_t storage_bytes_estimate() const;

  template <typename Fn> void for_each_set_bit(Fn &&fn) const {
    if (uses_small_storage()) {
      std::uint32_t mask = small_mask_;
      while (mask != 0u) {
        const std::size_t bit = static_cast<std::size_t>(std::countr_zero(mask));
        mask &= (mask - 1u);
        if (bit < shift_count_)
          fn(bit);
      }
      return;
    }

    for (std::size_t word_idx = 0; word_idx < words_.size(); ++word_idx) {
      std::uint64_t word = words_[word_idx];
      while (word != 0ull) {
        const std::size_t bit =
            static_cast<std::size_t>(std::countr_zero(word));
        word &= (word - 1ull);
        const std::size_t shift_idx = word_idx * kBitsPerWord + bit;
        if (shift_idx < shift_count_)
          fn(shift_idx);
      }
    }
  }

private:
  static constexpr std::size_t kSmallShiftLimit = 31;
  static constexpr std::size_t kBitsPerWord = 64;

  bool uses_small_storage() const { return shift_count_ <= kSmallShiftLimit; }
  static std::size_t word_count_for(std::size_t shift_count);
  void trim_unused_bits();

  std::size_t shift_count_ = 0;
  std::uint32_t small_mask_ = 0;
  std::vector<std::uint64_t> words_;
};

struct PbctopoSparseShiftPair {
  std::size_t rhs_shift = 0;
  PbctopoComponentPairShiftEntry entry;
};

class PbctopoComponentPairShiftTable {
public:
  std::size_t component_count = 0;
  std::size_t shift_count = 0;

  void reset(std::size_t new_component_count, std::size_t new_shift_count);
  void store_directed_pair(std::size_t lhs_comp, std::size_t lhs_shift,
                           std::size_t rhs_comp, std::size_t rhs_shift,
                           const PbctopoComponentPairShiftEntry &entry);
  const PbctopoComponentPairShiftEntry &
  entry_or_default(std::size_t lhs_comp, std::size_t lhs_shift,
                   std::size_t rhs_comp, std::size_t rhs_shift) const;
  const PbctopoComponentPairShiftEntry &
  at(std::size_t lhs_comp, std::size_t lhs_shift, std::size_t rhs_comp,
     std::size_t rhs_shift) const {
    return entry_or_default(lhs_comp, lhs_shift, rhs_comp, rhs_shift);
  }

  const PbctopoShiftMask &compat_mask(std::size_t lhs_comp,
                                      std::size_t lhs_shift,
                                      std::size_t rhs_comp) const;
  void add_bbox_prefiltered(std::size_t count) { bbox_prefiltered_ += count; }

  std::size_t dense_entries_equivalent() const;
  std::size_t materialized_entries() const { return materialized_entries_; }
  std::size_t bbox_prefiltered() const { return bbox_prefiltered_; }
  std::size_t sparse_rows() const;
  std::size_t peak_bytes_estimate() const;

private:
  struct Row {
    std::vector<std::vector<PbctopoSparseShiftPair>> by_lhs_shift;
    std::vector<PbctopoShiftMask> compat_by_lhs_shift;
    std::size_t materialized_entries = 0;
  };

  std::size_t row_index(std::size_t lhs_comp, std::size_t rhs_comp) const;
  Row &row(std::size_t lhs_comp, std::size_t rhs_comp);
  const Row &row(std::size_t lhs_comp, std::size_t rhs_comp) const;

  std::vector<Row> rows_;
  PbctopoComponentPairShiftEntry default_entry_;
  PbctopoShiftMask empty_mask_;
  std::size_t materialized_entries_ = 0;
  std::size_t bbox_prefiltered_ = 0;
};

} // namespace titan_pbctopo

#endif
