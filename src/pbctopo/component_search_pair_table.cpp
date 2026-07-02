#include "component_search_pair_table.h"

#include <algorithm>
#include <limits>
#include <utility>

namespace titan_pbctopo {

PbctopoShiftMask::PbctopoShiftMask(std::size_t shift_count) {
  reset(shift_count);
}

PbctopoShiftMask PbctopoShiftMask::full(std::size_t shift_count) {
  PbctopoShiftMask mask(shift_count);
  mask.set_all();
  return mask;
}

PbctopoShiftMask PbctopoShiftMask::only(std::size_t shift_count,
                                        std::size_t shift_idx) {
  PbctopoShiftMask mask(shift_count);
  mask.set_only(shift_idx);
  return mask;
}

void PbctopoShiftMask::reset(std::size_t shift_count) {
  shift_count_ = shift_count;
  small_mask_ = 0;
  if (uses_small_storage()) {
    words_.clear();
  } else {
    words_.assign(word_count_for(shift_count_), 0ull);
  }
}

void PbctopoShiftMask::set_all() {
  if (shift_count_ == 0) {
    clear_all();
    return;
  }
  if (uses_small_storage()) {
    small_mask_ = (1u << shift_count_) - 1u;
    return;
  }
  std::fill(words_.begin(), words_.end(), ~0ull);
  trim_unused_bits();
}

void PbctopoShiftMask::clear_all() {
  small_mask_ = 0;
  std::fill(words_.begin(), words_.end(), 0ull);
}

void PbctopoShiftMask::set_only(std::size_t shift_idx) {
  clear_all();
  set(shift_idx);
}

void PbctopoShiftMask::set(std::size_t shift_idx) {
  if (shift_idx >= shift_count_)
    return;
  if (uses_small_storage()) {
    small_mask_ |= (1u << shift_idx);
    return;
  }
  words_[shift_idx / kBitsPerWord] |=
      (1ull << (shift_idx % kBitsPerWord));
}

void PbctopoShiftMask::clear(std::size_t shift_idx) {
  if (shift_idx >= shift_count_)
    return;
  if (uses_small_storage()) {
    small_mask_ &= ~(1u << shift_idx);
    return;
  }
  words_[shift_idx / kBitsPerWord] &=
      ~(1ull << (shift_idx % kBitsPerWord));
}

bool PbctopoShiftMask::test(std::size_t shift_idx) const {
  if (shift_idx >= shift_count_)
    return false;
  if (uses_small_storage())
    return (small_mask_ & (1u << shift_idx)) != 0u;
  return (words_[shift_idx / kBitsPerWord] &
          (1ull << (shift_idx % kBitsPerWord))) != 0ull;
}

bool PbctopoShiftMask::empty() const {
  if (uses_small_storage())
    return small_mask_ == 0u;
  return std::all_of(words_.begin(), words_.end(),
                     [](std::uint64_t word) { return word == 0ull; });
}

std::size_t PbctopoShiftMask::count() const {
  if (uses_small_storage())
    return static_cast<std::size_t>(std::popcount(small_mask_));
  std::size_t total = 0;
  for (const std::uint64_t word : words_)
    total += static_cast<std::size_t>(std::popcount(word));
  return total;
}

bool PbctopoShiftMask::intersect_with(const PbctopoShiftMask &other) {
  bool changed = false;
  if (uses_small_storage() && other.uses_small_storage()) {
    const std::uint32_t previous = small_mask_;
    small_mask_ &= other.small_mask_;
    return previous != small_mask_;
  }

  if (words_.size() != other.words_.size()) {
    PbctopoShiftMask filtered(shift_count_);
    for_each_set_bit([&](std::size_t shift_idx) {
      if (other.test(shift_idx))
        filtered.set(shift_idx);
    });
    changed = filtered.count() != count();
    *this = std::move(filtered);
    return changed;
  }

  for (std::size_t idx = 0; idx < words_.size(); ++idx) {
    const std::uint64_t previous = words_[idx];
    words_[idx] &= other.words_[idx];
    changed = changed || previous != words_[idx];
  }
  return changed;
}

std::size_t PbctopoShiftMask::storage_bytes_estimate() const {
  return sizeof(PbctopoShiftMask) + words_.capacity() * sizeof(std::uint64_t);
}

std::size_t PbctopoShiftMask::word_count_for(std::size_t shift_count) {
  return (shift_count + kBitsPerWord - 1) / kBitsPerWord;
}

void PbctopoShiftMask::trim_unused_bits() {
  if (uses_small_storage() || shift_count_ == 0 || words_.empty())
    return;
  const std::size_t used_bits = shift_count_ % kBitsPerWord;
  if (used_bits == 0)
    return;
  const std::uint64_t mask = (1ull << used_bits) - 1ull;
  words_.back() &= mask;
}

void PbctopoComponentPairShiftTable::reset(
    std::size_t new_component_count, std::size_t new_shift_count) {
  component_count = new_component_count;
  shift_count = new_shift_count;
  materialized_entries_ = 0;
  bbox_prefiltered_ = 0;
  empty_mask_.reset(shift_count);
  rows_.clear();
  rows_.resize(component_count * component_count);
  for (Row &entry_row : rows_) {
    entry_row.by_lhs_shift.assign(shift_count, {});
    entry_row.compat_by_lhs_shift.clear();
    entry_row.compat_by_lhs_shift.reserve(shift_count);
    for (std::size_t shift_idx = 0; shift_idx < shift_count; ++shift_idx)
      entry_row.compat_by_lhs_shift.push_back(PbctopoShiftMask::full(shift_count));
    entry_row.materialized_entries = 0;
  }
}

void PbctopoComponentPairShiftTable::store_directed_pair(
    std::size_t lhs_comp, std::size_t lhs_shift, std::size_t rhs_comp,
    std::size_t rhs_shift, const PbctopoComponentPairShiftEntry &entry) {
  if (lhs_comp >= component_count || rhs_comp >= component_count ||
      lhs_shift >= shift_count || rhs_shift >= shift_count) {
    return;
  }

  Row &entry_row = row(lhs_comp, rhs_comp);
  auto &bucket = entry_row.by_lhs_shift[lhs_shift];
  for (PbctopoSparseShiftPair &sparse_pair : bucket) {
    if (sparse_pair.rhs_shift == rhs_shift) {
      sparse_pair.entry = entry;
      if (entry.hard_clash)
        entry_row.compat_by_lhs_shift[lhs_shift].clear(rhs_shift);
      else
        entry_row.compat_by_lhs_shift[lhs_shift].set(rhs_shift);
      return;
    }
  }

  bucket.push_back({rhs_shift, entry});
  ++entry_row.materialized_entries;
  ++materialized_entries_;
  if (entry.hard_clash)
    entry_row.compat_by_lhs_shift[lhs_shift].clear(rhs_shift);
}

const PbctopoComponentPairShiftEntry &
PbctopoComponentPairShiftTable::entry_or_default(
    std::size_t lhs_comp, std::size_t lhs_shift, std::size_t rhs_comp,
    std::size_t rhs_shift) const {
  if (lhs_comp >= component_count || rhs_comp >= component_count ||
      lhs_shift >= shift_count || rhs_shift >= shift_count) {
    return default_entry_;
  }
  const Row &entry_row = row(lhs_comp, rhs_comp);
  const auto &bucket = entry_row.by_lhs_shift[lhs_shift];
  for (const PbctopoSparseShiftPair &sparse_pair : bucket) {
    if (sparse_pair.rhs_shift == rhs_shift)
      return sparse_pair.entry;
  }
  return default_entry_;
}

const PbctopoShiftMask &PbctopoComponentPairShiftTable::compat_mask(
    std::size_t lhs_comp, std::size_t lhs_shift, std::size_t rhs_comp) const {
  if (lhs_comp >= component_count || rhs_comp >= component_count ||
      lhs_shift >= shift_count) {
    return empty_mask_;
  }
  return row(lhs_comp, rhs_comp).compat_by_lhs_shift[lhs_shift];
}

std::size_t PbctopoComponentPairShiftTable::dense_entries_equivalent() const {
  return component_count * component_count * shift_count * shift_count;
}

std::size_t PbctopoComponentPairShiftTable::sparse_rows() const {
  std::size_t non_empty = 0;
  for (const Row &entry_row : rows_) {
    if (entry_row.materialized_entries > 0)
      ++non_empty;
  }
  return non_empty;
}

std::size_t PbctopoComponentPairShiftTable::peak_bytes_estimate() const {
  std::size_t bytes = sizeof(PbctopoComponentPairShiftTable);
  bytes += rows_.capacity() * sizeof(Row);
  for (const Row &entry_row : rows_) {
    bytes += entry_row.by_lhs_shift.capacity() *
             sizeof(std::vector<PbctopoSparseShiftPair>);
    bytes += entry_row.compat_by_lhs_shift.capacity() *
             sizeof(PbctopoShiftMask);
    for (const auto &bucket : entry_row.by_lhs_shift)
      bytes += bucket.capacity() * sizeof(PbctopoSparseShiftPair);
    for (const PbctopoShiftMask &mask : entry_row.compat_by_lhs_shift)
      bytes += mask.storage_bytes_estimate();
  }
  return bytes;
}

std::size_t PbctopoComponentPairShiftTable::row_index(
    std::size_t lhs_comp, std::size_t rhs_comp) const {
  return lhs_comp * component_count + rhs_comp;
}

PbctopoComponentPairShiftTable::Row &
PbctopoComponentPairShiftTable::row(std::size_t lhs_comp,
                                    std::size_t rhs_comp) {
  return rows_[row_index(lhs_comp, rhs_comp)];
}

const PbctopoComponentPairShiftTable::Row &
PbctopoComponentPairShiftTable::row(std::size_t lhs_comp,
                                    std::size_t rhs_comp) const {
  return rows_[row_index(lhs_comp, rhs_comp)];
}

} // namespace titan_pbctopo
