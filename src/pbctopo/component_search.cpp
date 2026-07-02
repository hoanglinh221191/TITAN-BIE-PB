#include "component_search.h"
#include "component_graph.h"
#include "component_refine.h"
#include "component_search_pair_table.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <utility>

namespace titan_pbctopo {
namespace {

constexpr double kDistanceEps = 1.0e-12;
constexpr int kPackedCellBias = 1 << 20;
constexpr std::uint64_t kPackedCellMask = (1ull << 21) - 1ull;
constexpr std::size_t kPbctopoBeamFinalRescore = 4;

std::uint64_t pack_cell_key(int x, int y, int z) noexcept {
  const std::uint64_t ux =
      static_cast<std::uint32_t>(x + kPackedCellBias) & kPackedCellMask;
  const std::uint64_t uy =
      static_cast<std::uint32_t>(y + kPackedCellBias) & kPackedCellMask;
  const std::uint64_t uz =
      static_cast<std::uint32_t>(z + kPackedCellBias) & kPackedCellMask;
  return (ux << 42) ^ (uy << 21) ^ uz;
}

struct PbctopoSearchContinuityPrefix {
  std::vector<std::vector<std::size_t>> pairs_closing_at_depth;
};

PbctopoSearchContinuityPrefix build_search_continuity_prefix(
    const std::vector<PbctopoContinuityPair> &continuity_pairs,
    const std::vector<std::size_t> &atom_component_index,
    const std::vector<std::size_t> &search_order, std::size_t component_count) {
  PbctopoSearchContinuityPrefix prefix;
  prefix.pairs_closing_at_depth.assign(search_order.size(), {});
  if (continuity_pairs.empty() || atom_component_index.empty() ||
      search_order.empty())
    return prefix;

  std::vector<std::size_t> component_depth(
      component_count, std::numeric_limits<std::size_t>::max());
  for (std::size_t depth = 0; depth < search_order.size(); ++depth) {
    if (search_order[depth] < component_depth.size())
      component_depth[search_order[depth]] = depth;
  }

  for (std::size_t pair_idx = 0; pair_idx < continuity_pairs.size();
       ++pair_idx) {
    const auto &pair = continuity_pairs[pair_idx];
    if (pair.a >= atom_component_index.size() ||
        pair.b >= atom_component_index.size())
      continue;
    const std::size_t comp_a = atom_component_index[pair.a];
    const std::size_t comp_b = atom_component_index[pair.b];
    if (comp_a >= component_depth.size() || comp_b >= component_depth.size())
      continue;
    const std::size_t depth_a = component_depth[comp_a];
    const std::size_t depth_b = component_depth[comp_b];
    if (depth_a == std::numeric_limits<std::size_t>::max() ||
        depth_b == std::numeric_limits<std::size_t>::max())
      continue;
    const std::size_t close_depth = std::max(depth_a, depth_b);
    if (close_depth < prefix.pairs_closing_at_depth.size())
      prefix.pairs_closing_at_depth[close_depth].push_back(pair_idx);
  }
  return prefix;
}

bool consume_search_continuity_prefix(
    const std::vector<Vec3> &placed_points,
    const std::vector<PbctopoContinuityPair> &continuity_pairs,
    const std::vector<std::size_t> &pair_indices, double hard_cutoff2,
    const PbctopoComponentCandidate *incumbent, double &path_sum,
    double &max_d2, PbctopoTimingStats *timing) {
  if (pair_indices.empty())
    return true;

  using clock = std::chrono::steady_clock;
  const auto started = clock::now();
  auto finish = [&](bool ok) {
    if (timing != nullptr) {
      timing->component_search_prefix_continuity_s +=
          std::chrono::duration<double>(clock::now() - started).count();
    }
    return ok;
  };

  const bool has_incumbent =
      incumbent != nullptr && incumbent->valid &&
      std::isfinite(incumbent->continuity_path) &&
      std::isfinite(incumbent->continuity_max_d2);

  for (std::size_t pair_idx : pair_indices) {
    if (pair_idx >= continuity_pairs.size())
      continue;
    const auto &pair = continuity_pairs[pair_idx];
    if (pair.a >= placed_points.size() || pair.b >= placed_points.size())
      continue;
    const Vec3 &point_a = placed_points[pair.a];
    const Vec3 &point_b = placed_points[pair.b];
    const double dx = point_a.x - point_b.x;
    const double dy = point_a.y - point_b.y;
    const double dz = point_a.z - point_b.z;
    const double d2 = dx * dx + dy * dy + dz * dz;
    if (timing != nullptr)
      ++timing->component_search_prefix_continuity_checks;
    if (d2 > max_d2)
      max_d2 = d2;
    if (d2 > hard_cutoff2) {
      if (timing != nullptr) {
        ++timing->component_search_prefix_continuity_prunes;
        ++timing->component_search_prefix_continuity_hard_rejects;
      }
      return finish(false);
    }
    path_sum += std::sqrt(std::max(d2, 0.0));
    if (has_incumbent &&
        max_d2 > incumbent->continuity_max_d2 + 1.0e-6 &&
        path_sum > incumbent->continuity_path + 1.0e-6) {
      if (timing != nullptr)
        ++timing->component_search_prefix_continuity_prunes;
      return finish(false);
    }
  }
  return finish(true);
}

PbctopoComponentPairShiftEntry compute_component_pair_shift_entry(
    const std::vector<std::size_t> &lhs_component,
    const PbctopoShiftedComponent &lhs,
    const std::vector<std::size_t> &rhs_component,
    const PbctopoShiftedComponent &rhs, double clash_cutoff) {
  PbctopoComponentPairShiftEntry entry;
  entry.bbox_d2 = shifted_component_bbox_distance2(lhs, rhs);
  if (!component_bbox_may_clash(lhs, rhs, clash_cutoff))
    return entry;

  const double cutoff2 = clash_cutoff * clash_cutoff;
  const std::vector<std::size_t> *scan_component = &lhs_component;
  const std::vector<std::size_t> *grid_component = &rhs_component;
  const PbctopoShiftedComponent *scan_shifted = &lhs;
  const PbctopoShiftedComponent *grid_shifted = &rhs;
  if (lhs_component.size() > rhs_component.size()) {
    scan_component = &rhs_component;
    grid_component = &lhs_component;
    scan_shifted = &rhs;
    grid_shifted = &lhs;
  }

  for (std::size_t local_idx = 0; local_idx < scan_component->size();
       ++local_idx) {
    const std::size_t atom_idx = (*scan_component)[local_idx];
    const auto &point = scan_shifted->points[local_idx];
    const auto &key = scan_shifted->clash_keys[local_idx];
    const int cx = key.x;
    const int cy = key.y;
    const int cz = key.z;
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          const auto it = grid_shifted->clash_grid.find(
              pack_cell_key(cx + dx, cy + dy, cz + dz));
          if (it == grid_shifted->clash_grid.end())
            continue;
          for (std::size_t other_local_idx : it->second) {
            const std::size_t other_idx = (*grid_component)[other_local_idx];
            const auto &other = grid_shifted->points[other_local_idx];
            const double ddx = point.x - other.x;
            const double ddy = point.y - other.y;
            const double ddz = point.z - other.z;
            const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (d2 < entry.min_pair_d2) {
              entry.min_pair_d2 = d2;
              entry.clash_atom_a = atom_idx;
              entry.clash_atom_b = other_idx;
            }
            if (d2 + kDistanceEps < cutoff2) {
              entry.hard_clash = true;
              return entry;
            }
          }
        }
      }
    }
  }
  return entry;
}

PbctopoComponentPairShiftTable build_component_pair_shift_table(
    const std::vector<std::vector<std::size_t>> &components,
    const std::vector<std::vector<PbctopoShiftedComponent>> &shifted_components,
    double clash_cutoff, PbctopoTimingStats *timing = nullptr) {
  using clock = std::chrono::steady_clock;
  const auto started = clock::now();
  PbctopoComponentPairShiftTable table;
  const std::size_t component_count = components.size();
  const std::size_t shift_count =
      shifted_components.empty() ? 0 : shifted_components.front().size();
  table.reset(component_count, shift_count);

  for (std::size_t lhs_comp = 0; lhs_comp < component_count; ++lhs_comp) {
    for (std::size_t rhs_comp = lhs_comp + 1; rhs_comp < component_count;
         ++rhs_comp) {
      for (std::size_t lhs_shift = 0; lhs_shift < shift_count; ++lhs_shift) {
        for (std::size_t rhs_shift = 0; rhs_shift < shift_count; ++rhs_shift) {
          const auto &lhs_shifted = shifted_components[lhs_comp][lhs_shift];
          const auto &rhs_shifted = shifted_components[rhs_comp][rhs_shift];
          if (!component_bbox_may_clash(lhs_shifted, rhs_shifted,
                                        clash_cutoff)) {
            table.add_bbox_prefiltered(2);
            continue;
          }
          const auto entry = compute_component_pair_shift_entry(
              components[lhs_comp], lhs_shifted, components[rhs_comp],
              rhs_shifted, clash_cutoff);
          if (entry.hard_clash || std::isfinite(entry.min_pair_d2)) {
            table.store_directed_pair(lhs_comp, lhs_shift, rhs_comp, rhs_shift,
                                      entry);
            table.store_directed_pair(rhs_comp, rhs_shift, lhs_comp, lhs_shift,
                                      entry);
          }
        }
      }
    }
  }

  if (timing != nullptr) {
    ++timing->component_pair_table_builds;
    timing->component_pair_table_s +=
        std::chrono::duration<double>(clock::now() - started).count();
    timing->component_pair_table_component_count =
        std::max(timing->component_pair_table_component_count, component_count);
    timing->component_pair_table_shift_count =
        std::max(timing->component_pair_table_shift_count, shift_count);
    timing->component_pair_table_dense_entries_equivalent =
        std::max(timing->component_pair_table_dense_entries_equivalent,
                 table.dense_entries_equivalent());
    timing->component_pair_table_materialized_entries =
        std::max(timing->component_pair_table_materialized_entries,
                 table.materialized_entries());
    timing->component_pair_table_bbox_prefiltered += table.bbox_prefiltered();
    timing->component_pair_table_sparse_rows =
        std::max(timing->component_pair_table_sparse_rows, table.sparse_rows());
    timing->component_pair_table_peak_bytes_estimate =
        std::max(timing->component_pair_table_peak_bytes_estimate,
                 table.peak_bytes_estimate());
  }
  return table;
}

std::vector<std::size_t> build_component_search_order(
    const std::vector<std::vector<std::size_t>> &components,
    const std::vector<std::size_t> &anchor_indices,
    const std::vector<std::size_t> &atom_component_index,
    const PbctopoComponentPairShiftTable &pair_table,
    std::size_t identity_shift_pos) {
  std::vector<std::size_t> order;
  if (components.empty())
    return order;
  order.reserve(components.size());
  order.push_back(0);
  if (components.size() == 1)
    return order;

  std::vector<std::size_t> anchor_counts(components.size(), 0);
  for (std::size_t anchor_idx : anchor_indices) {
    if (anchor_idx < atom_component_index.size()) {
      const std::size_t comp_idx = atom_component_index[anchor_idx];
      if (comp_idx < anchor_counts.size())
        ++anchor_counts[comp_idx];
    }
  }

  std::vector<std::size_t> remaining;
  remaining.reserve(components.size() - 1);
  for (std::size_t comp_idx = 1; comp_idx < components.size(); ++comp_idx)
    remaining.push_back(comp_idx);

  std::stable_sort(remaining.begin(), remaining.end(),
                   [&](std::size_t lhs, std::size_t rhs) {
                     if (components[lhs].size() != components[rhs].size())
                       return components[lhs].size() > components[rhs].size();
                     if (anchor_counts[lhs] != anchor_counts[rhs])
                       return anchor_counts[lhs] > anchor_counts[rhs];

                     const std::size_t lhs_compatible_shifts =
                         pair_table.compat_mask(0, identity_shift_pos, lhs)
                             .count();
                     const std::size_t rhs_compatible_shifts =
                         pair_table.compat_mask(0, identity_shift_pos, rhs)
                             .count();
                     double lhs_best_min_pair =
                         std::numeric_limits<double>::infinity();
                     double rhs_best_min_pair =
                         std::numeric_limits<double>::infinity();
                     for (std::size_t shift_idx = 0;
                          shift_idx < pair_table.shift_count; ++shift_idx) {
                       const auto &lhs_entry =
                           pair_table.at(0, identity_shift_pos, lhs, shift_idx);
                       const auto &rhs_entry =
                           pair_table.at(0, identity_shift_pos, rhs, shift_idx);
                       lhs_best_min_pair =
                           std::min(lhs_best_min_pair, lhs_entry.min_pair_d2);
                       rhs_best_min_pair =
                           std::min(rhs_best_min_pair, rhs_entry.min_pair_d2);
                     }
                     if (lhs_compatible_shifts != rhs_compatible_shifts)
                       return lhs_compatible_shifts < rhs_compatible_shifts;
                     if (std::fabs(lhs_best_min_pair - rhs_best_min_pair) >
                         1.0e-9)
                       return lhs_best_min_pair < rhs_best_min_pair;
                     return lhs < rhs;
                   });

  order.insert(order.end(), remaining.begin(), remaining.end());
  return order;
}

struct PbctopoComponentSearchShiftOption {
  std::size_t shift_idx = 0;
  double next_min_pair_d2 = std::numeric_limits<double>::infinity();
  double ordering_min_pair_d2 = std::numeric_limits<double>::infinity();
  double shift_norm2 = std::numeric_limits<double>::infinity();
  std::size_t future_remaining_shifts =
      std::numeric_limits<std::size_t>::max();
  std::size_t clash_atom_a = std::numeric_limits<std::size_t>::max();
  std::size_t clash_atom_b = std::numeric_limits<std::size_t>::max();
  bool warm_priority = false;
};

struct PbctopoBeamState {
  std::vector<std::size_t> selected_shift;
  std::vector<PbctopoShiftMask> allowed_masks;
  double shift_norm2 = 0.0;
  double min_pair_d2 = std::numeric_limits<double>::infinity();
  std::size_t clash_atom_a = std::numeric_limits<std::size_t>::max();
  std::size_t clash_atom_b = std::numeric_limits<std::size_t>::max();
  std::size_t future_remaining_shifts =
      std::numeric_limits<std::size_t>::max();
  std::size_t warm_hits = 0;
};

std::vector<Int3> normalize_search_warm_offsets(
    const std::vector<Int3> &warm_offsets,
    const std::vector<std::vector<std::size_t>> &components) {
  if (warm_offsets.empty())
    return {};
  std::vector<Int3> rebased = warm_offsets;
  rebase_component_offsets(rebased, components);
  std::vector<Int3> normalized;
  normalized.reserve(rebased.size());
  const Int3 root = rebased.front();
  for (const auto &offset : rebased) {
    normalized.push_back(
        {offset.x - root.x, offset.y - root.y, offset.z - root.z});
  }
  return normalized;
}

std::vector<PbctopoComponentSearchShiftOption> build_component_search_shift_order(
    const std::vector<PbctopoShift> &shifts, std::size_t actual_comp_idx,
    std::size_t depth, const std::vector<std::size_t> &search_order,
    const std::vector<std::size_t> &selected_shift,
    const std::vector<PbctopoShiftMask> &allowed_masks,
    const PbctopoComponentPairShiftTable &pair_table,
    double current_min_pair_d2, std::size_t clash_atom_a,
    std::size_t clash_atom_b,
    const std::vector<Int3> *normalized_warm_offsets = nullptr,
    PbctopoTimingStats *timing = nullptr) {
  std::vector<PbctopoComponentSearchShiftOption> ordered;
  if (actual_comp_idx >= allowed_masks.size())
    return ordered;
  const PbctopoShiftMask current_mask = allowed_masks[actual_comp_idx];
  if (current_mask.empty()) {
    if (timing != nullptr) {
      ++timing->component_pair_table_prunes;
      ++timing->component_search_zero_current_mask_prunes;
    }
    return ordered;
  }
  ordered.reserve(current_mask.count());

  Int3 warm_delta;
  const bool has_warm_delta =
      (depth > 0 && normalized_warm_offsets != nullptr &&
       actual_comp_idx < normalized_warm_offsets->size())
          ? (warm_delta = (*normalized_warm_offsets)[actual_comp_idx], true)
          : false;

  current_mask.for_each_set_bit([&](std::size_t shift_idx) {
    if (shift_idx >= shifts.size())
      return;
    PbctopoComponentSearchShiftOption option;
    option.shift_idx = shift_idx;
    option.next_min_pair_d2 = current_min_pair_d2;
    option.ordering_min_pair_d2 = std::numeric_limits<double>::infinity();
    option.future_remaining_shifts = 0;
    option.shift_norm2 = static_cast<double>(
        shifts[shift_idx].ix * shifts[shift_idx].ix +
        shifts[shift_idx].iy * shifts[shift_idx].iy +
        shifts[shift_idx].iz * shifts[shift_idx].iz);
    option.clash_atom_a = clash_atom_a;
    option.clash_atom_b = clash_atom_b;
    option.warm_priority =
        has_warm_delta && shifts[shift_idx].ix == warm_delta.x &&
        shifts[shift_idx].iy == warm_delta.y &&
        shifts[shift_idx].iz == warm_delta.z;

    for (std::size_t prev_depth = 0; prev_depth < depth; ++prev_depth) {
      const std::size_t prev_comp_idx = search_order[prev_depth];
      const std::size_t prev_shift_idx = selected_shift[prev_comp_idx];
      if (prev_shift_idx >= pair_table.shift_count)
        continue;
      const auto &entry =
          pair_table.at(actual_comp_idx, shift_idx, prev_comp_idx,
                        prev_shift_idx);
      option.ordering_min_pair_d2 =
          std::min(option.ordering_min_pair_d2, entry.min_pair_d2);
      if (entry.min_pair_d2 < option.next_min_pair_d2) {
        option.next_min_pair_d2 = entry.min_pair_d2;
        option.clash_atom_a = entry.clash_atom_a;
        option.clash_atom_b = entry.clash_atom_b;
      }
    }

    bool future_zero_mask = false;
    for (std::size_t future_depth = depth + 1;
         future_depth < search_order.size(); ++future_depth) {
      const std::size_t future_comp_idx = search_order[future_depth];
      if (future_comp_idx >= allowed_masks.size())
        continue;
      PbctopoShiftMask next_mask = allowed_masks[future_comp_idx];
      next_mask.intersect_with(
          pair_table.compat_mask(actual_comp_idx, shift_idx, future_comp_idx));
      if (next_mask.empty()) {
        future_zero_mask = true;
        break;
      }
      option.future_remaining_shifts += next_mask.count();
    }
    if (future_zero_mask) {
      if (timing != nullptr)
        ++timing->component_pair_table_prunes;
      if (timing != nullptr)
        ++timing->component_search_zero_future_mask_prunes;
      return;
    }
    ordered.push_back(option);
  });

  std::sort(ordered.begin(), ordered.end(),
            [](const PbctopoComponentSearchShiftOption &lhs,
               const PbctopoComponentSearchShiftOption &rhs) {
              if (lhs.warm_priority != rhs.warm_priority)
                return lhs.warm_priority;
              if (lhs.future_remaining_shifts != rhs.future_remaining_shifts)
                return lhs.future_remaining_shifts < rhs.future_remaining_shifts;
              if (std::fabs(lhs.ordering_min_pair_d2 -
                            rhs.ordering_min_pair_d2) > 1.0e-9)
                return lhs.ordering_min_pair_d2 > rhs.ordering_min_pair_d2;
              if (std::fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-9)
                return lhs.shift_norm2 < rhs.shift_norm2;
              return lhs.shift_idx < rhs.shift_idx;
            });
  if (!ordered.empty() && ordered.front().warm_priority && timing != nullptr)
    ++timing->component_search_shift_warm_hits;
  return ordered;
}

bool better_component_beam_state(const PbctopoBeamState &lhs,
                                 const PbctopoBeamState &rhs) {
  const bool lhs_has_pair = std::isfinite(lhs.min_pair_d2);
  const bool rhs_has_pair = std::isfinite(rhs.min_pair_d2);
  if (lhs_has_pair != rhs_has_pair)
    return lhs_has_pair;
  if (lhs_has_pair && std::fabs(lhs.min_pair_d2 - rhs.min_pair_d2) > 1.0e-9)
    return lhs.min_pair_d2 > rhs.min_pair_d2;
  if (lhs.warm_hits != rhs.warm_hits)
    return lhs.warm_hits > rhs.warm_hits;
  if (lhs.future_remaining_shifts != rhs.future_remaining_shifts)
    return lhs.future_remaining_shifts < rhs.future_remaining_shifts;
  if (std::fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-9)
    return lhs.shift_norm2 < rhs.shift_norm2;
  return lhs.selected_shift < rhs.selected_shift;
}

bool build_component_points_from_selected_shifts(
    const std::vector<PbctopoAtom> &atoms,
    const std::vector<std::vector<std::size_t>> &components,
    const std::vector<std::vector<PbctopoShiftedComponent>> &shifted_components,
    const std::vector<std::size_t> &selected_shift,
    std::vector<Vec3> &placed_points, double &shift_norm2) {
  if (components.size() != selected_shift.size())
    return false;
  placed_points.resize(atoms.size());
  shift_norm2 = 0.0;
  for (std::size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
    const std::size_t shift_idx = selected_shift[comp_idx];
    if (shift_idx >= shifted_components[comp_idx].size())
      return false;
    const auto &component = components[comp_idx];
    const auto &shifted = shifted_components[comp_idx][shift_idx];
    if (shifted.points.size() != component.size())
      return false;
    shift_norm2 += shifted.weighted_shift_norm2;
    for (std::size_t local_idx = 0; local_idx < component.size(); ++local_idx)
      placed_points[component[local_idx]] = shifted.points[local_idx];
  }
  return true;
}

std::size_t identity_shift_index(const std::vector<PbctopoShift> &shifts) {
  for (std::size_t shift_idx = 0; shift_idx < shifts.size(); ++shift_idx) {
    if (shifts[shift_idx].ix == 0 && shifts[shift_idx].iy == 0 &&
        shifts[shift_idx].iz == 0) {
      return shift_idx;
    }
  }
  return 0;
}

} // namespace

PbctopoComponentCandidate search_component_shift_candidate_beam(
    const std::vector<PbctopoAtom> &atoms, box Box,
    const std::vector<std::vector<std::size_t>> &components,
    const std::vector<PbctopoShift> &shifts,
    const std::vector<PbctopoContinuityPair> &continuity_pairs,
    const std::vector<std::size_t> &anchor_indices,
    const RuntimeCutoffs &cutoffs, std::vector<Vec3> &best_points_out,
    std::size_t beam_width, const PbctopoComponentSearchCallbacks &callbacks,
    const std::vector<Int3> *warm_component_offsets,
    const PbctopoCutCandidate *incumbent_cut,
    const std::vector<Vec3> *incumbent_points, PbctopoTimingStats *timing,
    AnchorBackend anchor_backend,
    const PbctopoAnchorGraphContext *anchor_context) {
  (void)Box;
  beam_width = std::max<std::size_t>(
      1, std::min(beam_width, static_cast<std::size_t>(kMaxBeamWidth)));
  PbctopoComponentCandidate best;
  if (components.empty() || !callbacks.ready())
    return best;

  std::vector<Vec3> placed_points(atoms.size());
  best_points_out.clear();
  best_points_out.reserve(atoms.size());

  const auto shifted_components =
      precompute_shifted_components(atoms, components, shifts,
                                    cutoffs.clash_cutoff,
                                    cutoffs.pair_contact_cutoff);

  const std::size_t identity_shift_pos = identity_shift_index(shifts);

  if (incumbent_cut != nullptr && incumbent_cut->valid &&
      !incumbent_cut->ultra_clash && incumbent_points != nullptr &&
      incumbent_points->size() == atoms.size()) {
    best = callbacks.build_seed(*incumbent_cut, components, *incumbent_points);
    if (timing != nullptr)
      ++timing->component_search_seed_graph;
  } else if (timing != nullptr) {
    ++timing->component_search_seed_none;
  }

  const PbctopoComponentPairShiftTable pair_table =
      build_component_pair_shift_table(components, shifted_components,
                                       cutoffs.clash_cutoff, timing);
  const std::vector<std::size_t> atom_component_index =
      build_atom_component_index(atoms.size(), components);
  const std::vector<std::size_t> search_order = build_component_search_order(
      components, anchor_indices, atom_component_index, pair_table,
      identity_shift_pos);
  const std::vector<Int3> normalized_warm_offsets =
      (warm_component_offsets != nullptr &&
       warm_component_offsets->size() == components.size())
          ? normalize_search_warm_offsets(*warm_component_offsets, components)
          : std::vector<Int3>{};
  const std::vector<Int3> *normalized_warm_offsets_ptr =
      normalized_warm_offsets.empty() ? nullptr : &normalized_warm_offsets;

  const PbctopoShiftMask full_shift_mask = PbctopoShiftMask::full(shifts.size());
  std::vector<PbctopoShiftMask> initial_allowed_masks(components.size(),
                                                      full_shift_mask);
  std::vector<std::size_t> initial_selected_shift(components.size(),
                                                  shifts.size());
  if (!search_order.empty() && identity_shift_pos < shifts.size()) {
    const std::size_t root_comp_idx = search_order.front();
    initial_allowed_masks[root_comp_idx] =
        PbctopoShiftMask::only(shifts.size(), identity_shift_pos);
    initial_selected_shift[root_comp_idx] = identity_shift_pos;
    for (std::size_t depth = 1; depth < search_order.size(); ++depth) {
      const std::size_t future_comp_idx = search_order[depth];
      initial_allowed_masks[future_comp_idx].intersect_with(
          pair_table.compat_mask(root_comp_idx, identity_shift_pos,
                                 future_comp_idx));
      if (initial_allowed_masks[future_comp_idx].empty() && timing != nullptr) {
        ++timing->component_pair_table_prunes;
        ++timing->component_search_zero_future_mask_prunes;
      }
    }
  }

  std::vector<PbctopoBeamState> beam_states;
  PbctopoBeamState root_state;
  root_state.selected_shift = initial_selected_shift;
  root_state.allowed_masks = initial_allowed_masks;
  root_state.future_remaining_shifts = 0;
  beam_states.push_back(std::move(root_state));
  if (timing != nullptr)
    ++timing->component_search_depth_visits[0];

  for (std::size_t depth = search_order.empty() ? 0 : 1;
       depth < search_order.size(); ++depth) {
    std::vector<PbctopoBeamState> next_states;
    next_states.reserve(beam_states.size() *
                        std::min<std::size_t>(beam_width, shifts.size()));
    const std::size_t actual_comp_idx = search_order[depth];
    for (const auto &state : beam_states) {
      if (timing != nullptr &&
          depth < timing->component_search_depth_visits.size()) {
        ++timing->component_search_depth_visits[depth];
      }
      const std::vector<PbctopoComponentSearchShiftOption> ordered_shifts =
          build_component_search_shift_order(
              shifts, actual_comp_idx, depth, search_order, state.selected_shift,
              state.allowed_masks, pair_table, state.min_pair_d2,
              state.clash_atom_a, state.clash_atom_b,
              normalized_warm_offsets_ptr, timing);
      for (const auto &option : ordered_shifts) {
        if (timing != nullptr) {
          ++timing->component_shift_trials;
          ++timing->component_large_mode_beam_expanded;
        }
        PbctopoBeamState next_state = state;
        next_state.selected_shift[actual_comp_idx] = option.shift_idx;
        if (actual_comp_idx < next_state.allowed_masks.size() &&
            option.shift_idx < shifts.size()) {
          next_state.allowed_masks[actual_comp_idx] =
              PbctopoShiftMask::only(shifts.size(), option.shift_idx);
        }
        for (std::size_t future_depth = depth + 1;
             future_depth < search_order.size(); ++future_depth) {
          const std::size_t future_comp_idx = search_order[future_depth];
          next_state.allowed_masks[future_comp_idx].intersect_with(
              pair_table.compat_mask(actual_comp_idx, option.shift_idx,
                                     future_comp_idx));
        }
        next_state.shift_norm2 +=
            shifted_components[actual_comp_idx][option.shift_idx]
                .weighted_shift_norm2;
        next_state.min_pair_d2 = option.next_min_pair_d2;
        next_state.clash_atom_a = option.clash_atom_a;
        next_state.clash_atom_b = option.clash_atom_b;
        next_state.future_remaining_shifts = option.future_remaining_shifts;
        next_state.warm_hits += option.warm_priority ? 1u : 0u;
        next_states.push_back(std::move(next_state));
      }
    }
    if (next_states.empty())
      return best;
    if (next_states.size() > beam_width) {
      std::partial_sort(
          next_states.begin(),
          next_states.begin() + static_cast<std::ptrdiff_t>(beam_width),
          next_states.end(), better_component_beam_state);
      next_states.resize(beam_width);
    } else {
      std::sort(next_states.begin(), next_states.end(),
                better_component_beam_state);
    }
    if (timing != nullptr)
      timing->component_large_mode_beam_kept += next_states.size();
    beam_states = std::move(next_states);
  }

  if (timing != nullptr)
    timing->component_large_mode_beam_complete_assignments += beam_states.size();
  const std::size_t exact_rescore_count =
      std::min(kPbctopoBeamFinalRescore, beam_states.size());
  if (beam_states.size() > exact_rescore_count) {
    std::partial_sort(
        beam_states.begin(),
        beam_states.begin() + static_cast<std::ptrdiff_t>(exact_rescore_count),
        beam_states.end(), better_component_beam_state);
  } else {
    std::sort(beam_states.begin(), beam_states.end(),
              better_component_beam_state);
  }
  if (timing != nullptr)
    timing->component_large_mode_beam_final_rescored += exact_rescore_count;
  for (std::size_t state_idx = 0; state_idx < exact_rescore_count;
       ++state_idx) {
    double candidate_shift_norm2 = 0.0;
    if (!build_component_points_from_selected_shifts(
            atoms, components, shifted_components,
            beam_states[state_idx].selected_shift, placed_points,
            candidate_shift_norm2)) {
      continue;
    }
    PbctopoComponentCandidate candidate = callbacks.score_assignment(
        atoms, placed_points, components, beam_states[state_idx].min_pair_d2,
        beam_states[state_idx].clash_atom_a,
        beam_states[state_idx].clash_atom_b, continuity_pairs, anchor_indices,
        cutoffs.continuity_cutoff, candidate_shift_norm2,
        best.valid ? &best : nullptr, timing, anchor_backend, anchor_context);
    if (callbacks.better_candidate(candidate, best)) {
      best = candidate;
      best_points_out = placed_points;
    }
  }

  return best;
}

PbctopoComponentCandidate search_component_shift_candidate(
    const std::vector<PbctopoAtom> &atoms, box Box,
    const std::vector<std::vector<std::size_t>> &components,
    const std::vector<PbctopoShift> &shifts,
    const std::vector<PbctopoContinuityPair> &continuity_pairs,
    const std::vector<std::size_t> &anchor_indices,
    const RuntimeCutoffs &cutoffs, std::vector<Vec3> &best_points_out,
    const PbctopoComponentSearchCallbacks &callbacks,
    const std::vector<Int3> *warm_component_offsets,
    const PbctopoCutCandidate *incumbent_cut,
    const std::vector<Vec3> *incumbent_points, PbctopoTimingStats *timing,
    AnchorBackend anchor_backend,
    const PbctopoAnchorGraphContext *anchor_context) {
  (void)Box;
  PbctopoComponentCandidate best;
  if (components.empty() || !callbacks.ready())
    return best;

  std::vector<Vec3> placed_points(atoms.size());
  best_points_out.clear();
  best_points_out.reserve(atoms.size());
  std::vector<std::size_t> selected_shift(components.size(), shifts.size());
  const auto shifted_components =
      precompute_shifted_components(atoms, components, shifts,
                                    cutoffs.clash_cutoff,
                                    cutoffs.pair_contact_cutoff);

  const std::size_t identity_shift_pos = identity_shift_index(shifts);

  for (const auto &component : components) {
    for (std::size_t atom_idx : component) {
      placed_points[atom_idx] = {atoms[atom_idx].x, atoms[atom_idx].y,
                                 atoms[atom_idx].z};
    }
  }

  if (incumbent_cut != nullptr && incumbent_cut->valid &&
      !incumbent_cut->ultra_clash && incumbent_points != nullptr &&
      incumbent_points->size() == atoms.size()) {
    best = callbacks.build_seed(*incumbent_cut, components, *incumbent_points);
    if (timing != nullptr)
      ++timing->component_search_seed_graph;
  } else if (timing != nullptr) {
    ++timing->component_search_seed_none;
  }

  const PbctopoComponentPairShiftTable pair_table =
      build_component_pair_shift_table(components, shifted_components,
                                       cutoffs.clash_cutoff, timing);
  const std::vector<std::size_t> atom_component_index =
      build_atom_component_index(atoms.size(), components);
  const std::vector<std::size_t> search_order = build_component_search_order(
      components, anchor_indices, atom_component_index, pair_table,
      identity_shift_pos);
  const PbctopoSearchContinuityPrefix continuity_prefix =
      build_search_continuity_prefix(continuity_pairs, atom_component_index,
                                     search_order, components.size());
  const std::vector<Int3> normalized_warm_offsets =
      (warm_component_offsets != nullptr &&
       warm_component_offsets->size() == components.size())
          ? normalize_search_warm_offsets(*warm_component_offsets, components)
          : std::vector<Int3>{};
  const std::vector<Int3> *normalized_warm_offsets_ptr =
      normalized_warm_offsets.empty() ? nullptr : &normalized_warm_offsets;

  const PbctopoShiftMask full_shift_mask = PbctopoShiftMask::full(shifts.size());
  std::vector<PbctopoShiftMask> initial_allowed_masks(components.size(),
                                                      full_shift_mask);
  if (!search_order.empty() && identity_shift_pos < shifts.size()) {
    const std::size_t root_comp_idx = search_order.front();
    initial_allowed_masks[root_comp_idx] =
        PbctopoShiftMask::only(shifts.size(), identity_shift_pos);
    selected_shift[root_comp_idx] = identity_shift_pos;
    for (std::size_t depth = 1; depth < search_order.size(); ++depth) {
      const std::size_t future_comp_idx = search_order[depth];
      initial_allowed_masks[future_comp_idx].intersect_with(
          pair_table.compat_mask(root_comp_idx, identity_shift_pos,
                                 future_comp_idx));
      if (initial_allowed_masks[future_comp_idx].empty() && timing != nullptr) {
        ++timing->component_pair_table_prunes;
        ++timing->component_search_zero_future_mask_prunes;
      }
    }
  }

  double initial_continuity_path = 0.0;
  double initial_continuity_max_d2 = 0.0;
  const double continuity_cutoff2 =
      cutoffs.continuity_cutoff * cutoffs.continuity_cutoff;
  if (!search_order.empty() &&
      !consume_search_continuity_prefix(
          placed_points, continuity_pairs,
          continuity_prefix.pairs_closing_at_depth[0], continuity_cutoff2,
          best.valid ? &best : nullptr, initial_continuity_path,
          initial_continuity_max_d2, timing)) {
    return best;
  }

  std::function<void(std::size_t, double, double, double, double, std::size_t,
                     std::size_t, std::vector<PbctopoShiftMask> &)>
      dfs = [&](std::size_t depth, double shift_norm2,
                double continuity_path_prefix,
                double continuity_max_d2_prefix, double min_pair_d2,
                std::size_t clash_atom_a, std::size_t clash_atom_b,
                std::vector<PbctopoShiftMask> &allowed_masks) {
        if (timing != nullptr &&
            depth < timing->component_search_depth_visits.size()) {
          ++timing->component_search_depth_visits[depth];
        }
        if (depth >= search_order.size()) {
          PbctopoComponentCandidate candidate = callbacks.score_assignment(
              atoms, placed_points, components, min_pair_d2, clash_atom_a,
              clash_atom_b, continuity_pairs, anchor_indices,
              cutoffs.continuity_cutoff, shift_norm2,
              best.valid ? &best : nullptr, timing, anchor_backend,
              anchor_context);
          if (callbacks.better_candidate(candidate, best)) {
            best = candidate;
            best_points_out = placed_points;
          }
          return;
        }

        const std::size_t actual_comp_idx = search_order[depth];
        const std::vector<PbctopoComponentSearchShiftOption> ordered_shifts =
            build_component_search_shift_order(
                shifts, actual_comp_idx, depth, search_order, selected_shift,
                allowed_masks, pair_table, min_pair_d2, clash_atom_a,
                clash_atom_b, normalized_warm_offsets_ptr, timing);
        for (const auto &option : ordered_shifts) {
          const std::size_t shift_idx = option.shift_idx;
          if (timing != nullptr)
            ++timing->component_shift_trials;
          const auto &component = components[actual_comp_idx];
          const auto &shifted = shifted_components[actual_comp_idx][shift_idx];
          for (std::size_t local_idx = 0; local_idx < component.size();
               ++local_idx) {
            placed_points[component[local_idx]] = shifted.points[local_idx];
          }
          selected_shift[actual_comp_idx] = shift_idx;
          std::vector<std::pair<std::size_t, PbctopoShiftMask>>
              allowed_mask_restore;
          allowed_mask_restore.reserve(search_order.size() - depth);
          if (actual_comp_idx < allowed_masks.size() &&
              shift_idx < shifts.size()) {
            allowed_mask_restore.push_back(
                {actual_comp_idx, allowed_masks[actual_comp_idx]});
            allowed_masks[actual_comp_idx] =
                PbctopoShiftMask::only(shifts.size(), shift_idx);
          }
          for (std::size_t future_depth = depth + 1;
               future_depth < search_order.size(); ++future_depth) {
            const std::size_t future_comp_idx = search_order[future_depth];
            PbctopoShiftMask updated_mask = allowed_masks[future_comp_idx];
            if (updated_mask.intersect_with(pair_table.compat_mask(
                    actual_comp_idx, shift_idx, future_comp_idx))) {
              allowed_mask_restore.push_back(
                  {future_comp_idx, allowed_masks[future_comp_idx]});
              allowed_masks[future_comp_idx] = std::move(updated_mask);
            }
          }
          double next_continuity_path = continuity_path_prefix;
          double next_continuity_max_d2 = continuity_max_d2_prefix;
          bool prefix_ok = true;
          if (depth < continuity_prefix.pairs_closing_at_depth.size()) {
            prefix_ok = consume_search_continuity_prefix(
                placed_points, continuity_pairs,
                continuity_prefix.pairs_closing_at_depth[depth],
                continuity_cutoff2, best.valid ? &best : nullptr,
                next_continuity_path, next_continuity_max_d2, timing);
          }
          if (prefix_ok) {
            dfs(depth + 1, shift_norm2 + shifted.weighted_shift_norm2,
                next_continuity_path, next_continuity_max_d2,
                option.next_min_pair_d2, option.clash_atom_a,
                option.clash_atom_b, allowed_masks);
          }
          for (auto restore_it = allowed_mask_restore.rbegin();
               restore_it != allowed_mask_restore.rend(); ++restore_it) {
            allowed_masks[restore_it->first] = restore_it->second;
          }
          selected_shift[actual_comp_idx] = shifts.size();
        }
      };

  dfs(search_order.empty() ? 0 : 1, 0.0, initial_continuity_path,
      initial_continuity_max_d2, std::numeric_limits<double>::infinity(),
      std::numeric_limits<std::size_t>::max(),
      std::numeric_limits<std::size_t>::max(), initial_allowed_masks);
  return best;
}

} // namespace titan_pbctopo
