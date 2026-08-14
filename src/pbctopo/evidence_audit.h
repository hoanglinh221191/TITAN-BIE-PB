#ifndef TITAN_PBCTOPO_EVIDENCE_AUDIT_H
#define TITAN_PBCTOPO_EVIDENCE_AUDIT_H

#include "types.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace titan_pbctopo {

inline constexpr std::string_view PBCTOPO_EOBS_PAIR_HASH_ALGORITHM =
    "sha256-128";
inline constexpr std::string_view PBCTOPO_EOBS_PAIR_HASH_DOMAIN =
    "titan-pbctopo-eobs-pairs-v1";
inline constexpr std::string_view PBCTOPO_EOBS_HYPOTHESIS_HASH_DOMAIN =
    "titan-pbctopo-eobs-shift-hypotheses-v1";
inline constexpr std::string_view PBCTOPO_EOBS_CUTOFF_POLICY =
    "wrapped-mic-heavy-closest8-per-owner-fragment-pair-v1";
inline constexpr std::uint32_t PBCTOPO_EOBS_AUDIT_POLICY_VERSION = 1;
inline constexpr std::size_t PBCTOPO_EOBS_MAX_PAIRS_PER_UNIT_PAIR = 8;
inline constexpr double PBCTOPO_EOBS_SELECTION_TIE_TOLERANCE_D2_A2 = 1.0e-9;
inline constexpr double PBCTOPO_EOBS_LOSS_TOLERANCE_D2_A2 = 1.0e-6;

struct PbctopoCanonicalAtomPair {
  std::uint64_t a = 0;
  std::uint64_t b = 0;

  auto operator<=>(const PbctopoCanonicalAtomPair &) const = default;
};

struct PbctopoSelectedEvidenceAudit {
  PbctopoEvidenceState state = PbctopoEvidenceState::NotEvaluated;
  bool evaluated = false;
  bool consistent = false;
  bool evidence_graph_connected = false;
  std::size_t contacts_checked = 0;
  std::size_t contacts_lost = 0;
  PbctopoHash128 observed_pair_hash;
  PbctopoHash128 lost_pair_hash;
  PbctopoHash128 all_hypotheses_hash;
  PbctopoHash128 selected_hypothesis_hash;
  PbctopoHash128 selected_compatible_pair_hash;
  std::size_t hypothesis_count = 0;
  std::size_t component_relation_count = 0;
  std::size_t ambiguous_component_relation_count = 0;
  std::size_t compatible_hypothesis_count = 0;
  std::size_t supported_relation_count = 0;
  std::size_t selected_compatible_contact_count = 0;
  std::size_t alternative_hypothesis_contact_count = 0;
  std::size_t selected_compatible_contact_loss_count = 0;
  std::size_t no_supported_contact_relation_count = 0;
  bool selected_relation_has_support = false;
};

[[nodiscard]] PbctopoContactHypothesisSet
build_pbctopo_contact_hypotheses(
    std::span<const PbctopoAtom> atoms,
    std::span<const PbctopoAnchorContactEdge> observed_edges,
    std::span<const std::size_t> atom_component_index,
    std::span<const Int3> spatial_lift_atom_images,
    const PbctopoImageAssignment &identity_assignment);

[[nodiscard]] std::expected<PbctopoHash128, std::string>
hash_pbctopo_contact_hypotheses(
    std::span<const PbctopoContactShiftHypothesis> hypotheses);

[[nodiscard]] std::expected<PbctopoHash128, std::string>
hash_pbctopo_canonical_atom_pairs(
    std::span<const PbctopoCanonicalAtomPair> pairs);

template <class Point>
[[nodiscard]] std::expected<PbctopoSelectedEvidenceAudit, std::string>
audit_pbctopo_selected_evidence(
    std::span<const PbctopoAtom> atoms, std::span<const Point> placed_points,
    std::span<const PbctopoAnchorContactEdge> observed_edges,
    bool hard_feasible, bool evidence_graph_connected) {
  PbctopoSelectedEvidenceAudit result;
  result.evidence_graph_connected = evidence_graph_connected;
  if (!hard_feasible)
    return result;
  if (atoms.size() != placed_points.size())
    return std::unexpected("selected evidence audit atom/point size mismatch");

  std::vector<PbctopoCanonicalAtomPair> observed_pairs;
  std::vector<PbctopoCanonicalAtomPair> lost_pairs;
  observed_pairs.reserve(observed_edges.size());
  lost_pairs.reserve(observed_edges.size());
  for (const auto &edge : observed_edges) {
    if (edge.a >= atoms.size() || edge.b >= atoms.size() || edge.a == edge.b ||
        !std::isfinite(edge.cutoff_d2) || edge.cutoff_d2 < 0.0) {
      return std::unexpected("selected evidence audit contains an invalid edge");
    }
    const std::uint64_t global_a =
        static_cast<std::uint64_t>(atoms[edge.a].global_atom_id);
    const std::uint64_t global_b =
        static_cast<std::uint64_t>(atoms[edge.b].global_atom_id);
    const PbctopoCanonicalAtomPair pair{std::min(global_a, global_b),
                                        std::max(global_a, global_b)};
    observed_pairs.push_back(pair);
    const double dx = placed_points[edge.a].x - placed_points[edge.b].x;
    const double dy = placed_points[edge.a].y - placed_points[edge.b].y;
    const double dz = placed_points[edge.a].z - placed_points[edge.b].z;
    const double d2 = dx * dx + dy * dy + dz * dz;
    if (!std::isfinite(d2))
      return std::unexpected("selected evidence audit found non-finite distance");
    if (d2 > edge.cutoff_d2 + PBCTOPO_EOBS_LOSS_TOLERANCE_D2_A2)
      lost_pairs.push_back(pair);
  }
  std::sort(observed_pairs.begin(), observed_pairs.end());
  std::sort(lost_pairs.begin(), lost_pairs.end());
  if (std::adjacent_find(observed_pairs.begin(), observed_pairs.end()) !=
      observed_pairs.end()) {
    return std::unexpected("selected evidence audit contains duplicate pairs");
  }
  const auto observed_hash =
      hash_pbctopo_canonical_atom_pairs(observed_pairs);
  if (!observed_hash)
    return std::unexpected(observed_hash.error());
  const auto lost_hash = hash_pbctopo_canonical_atom_pairs(lost_pairs);
  if (!lost_hash)
    return std::unexpected(lost_hash.error());
  result.observed_pair_hash = *observed_hash;
  result.lost_pair_hash = *lost_hash;
  result.evaluated = true;
  result.contacts_checked = observed_pairs.size();
  result.contacts_lost = lost_pairs.size();
  result.consistent = lost_pairs.empty();
  result.state = !result.consistent
                     ? PbctopoEvidenceState::Contradicted
                     : evidence_graph_connected
                           ? PbctopoEvidenceState::SupportedConnected
                           : PbctopoEvidenceState::WeakDisconnected;
  return result;
}

template <class Point>
[[nodiscard]] std::expected<PbctopoSelectedEvidenceAudit, std::string>
audit_pbctopo_selected_hypotheses(
    std::span<const PbctopoAtom> atoms, std::span<const Point> placed_points,
    const PbctopoContactHypothesisSet &hypothesis_set,
    const PbctopoImageAssignment &selected_assignment, bool hard_feasible) {
  PbctopoSelectedEvidenceAudit result;
  if (!hard_feasible)
    return result;
  if (!hypothesis_set.complete)
    return result;
  if (!selected_assignment.valid || atoms.size() != placed_points.size() ||
      selected_assignment.component_offsets.size() !=
          selected_assignment.component_ids.size()) {
    return std::unexpected("selected hypothesis audit input mismatch");
  }

  result.all_hypotheses_hash = hypothesis_set.all_hypotheses_hash;
  result.observed_pair_hash = hypothesis_set.raw_observed_pair_hash;
  result.contacts_checked = hypothesis_set.raw_observed_pair_count;
  result.hypothesis_count = hypothesis_set.hypotheses.size();
  result.component_relation_count = hypothesis_set.component_relation_count;
  result.ambiguous_component_relation_count =
      hypothesis_set.ambiguous_component_relation_count;

  auto offset_for_id = [&](std::uint64_t id) -> const Int3 * {
    for (std::size_t idx = 0; idx < selected_assignment.component_ids.size(); ++idx) {
      if (selected_assignment.component_ids[idx] == id)
        return &selected_assignment.component_offsets[idx];
    }
    return nullptr;
  };

  std::vector<PbctopoContactShiftHypothesis> compatible;
  std::vector<PbctopoCanonicalAtomPair> observed_pairs;
  std::vector<PbctopoCanonicalAtomPair> lost_pairs;
  std::vector<std::pair<std::size_t, std::size_t>> connected_edges;
  bool unsupported_relation = false;
  std::size_t relation_begin = 0;
  while (relation_begin < hypothesis_set.hypotheses.size()) {
    std::size_t relation_end = relation_begin + 1;
    const auto &first = hypothesis_set.hypotheses[relation_begin];
    while (relation_end < hypothesis_set.hypotheses.size() &&
           hypothesis_set.hypotheses[relation_end].component_id_a ==
               first.component_id_a &&
           hypothesis_set.hypotheses[relation_end].component_id_b ==
               first.component_id_b) {
      ++relation_end;
    }
    const Int3 *offset_a = offset_for_id(first.component_id_a);
    const Int3 *offset_b = offset_for_id(first.component_id_b);
    if (offset_a == nullptr || offset_b == nullptr)
      return std::unexpected("selected assignment misses an evidence component");
    Int3 selected_relative{};
    if (!checked_int3_subtract(*offset_b, *offset_a, selected_relative))
      return std::unexpected("selected evidence relative shift overflow");
    bool matched = false;
    for (std::size_t idx = relation_begin; idx < relation_end; ++idx) {
      const auto &hypothesis = hypothesis_set.hypotheses[idx];
      if (!same_delta(hypothesis.relative_shift, selected_relative))
      {
        result.alternative_hypothesis_contact_count += hypothesis.edges.size();
        continue;
      }
      matched = true;
      compatible.push_back(hypothesis);
      connected_edges.emplace_back(hypothesis.component_a,
                                   hypothesis.component_b);
      for (const auto &edge : hypothesis.edges) {
        if (edge.a >= atoms.size() || edge.b >= atoms.size() || edge.a == edge.b)
          return std::unexpected("selected hypothesis contains an invalid edge");
        const auto global_a = static_cast<std::uint64_t>(atoms[edge.a].global_atom_id);
        const auto global_b = static_cast<std::uint64_t>(atoms[edge.b].global_atom_id);
        const PbctopoCanonicalAtomPair pair{std::min(global_a, global_b),
                                            std::max(global_a, global_b)};
        observed_pairs.push_back(pair);
        const double dx = placed_points[edge.a].x - placed_points[edge.b].x;
        const double dy = placed_points[edge.a].y - placed_points[edge.b].y;
        const double dz = placed_points[edge.a].z - placed_points[edge.b].z;
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (!std::isfinite(d2))
          return std::unexpected("selected hypothesis audit found non-finite distance");
        if (d2 > edge.cutoff_d2 + PBCTOPO_EOBS_LOSS_TOLERANCE_D2_A2)
          lost_pairs.push_back(pair);
      }
    }
    unsupported_relation = unsupported_relation || !matched;
    if (matched) ++result.supported_relation_count;
    if (!matched) ++result.no_supported_contact_relation_count;
    relation_begin = relation_end;
  }

  std::sort(observed_pairs.begin(), observed_pairs.end());
  std::sort(lost_pairs.begin(), lost_pairs.end());
  if (std::adjacent_find(observed_pairs.begin(), observed_pairs.end()) !=
          observed_pairs.end() ||
      std::adjacent_find(lost_pairs.begin(), lost_pairs.end()) !=
          lost_pairs.end()) {
    return std::unexpected("selected hypothesis audit contains duplicate pairs");
  }
  auto observed_hash = hash_pbctopo_canonical_atom_pairs(observed_pairs);
  auto lost_hash = hash_pbctopo_canonical_atom_pairs(lost_pairs);
  auto compatible_hashable = compatible;
  for (auto &hypothesis : compatible_hashable) {
    for (auto &edge : hypothesis.edges) {
      edge.a = atoms[edge.a].global_atom_id;
      edge.b = atoms[edge.b].global_atom_id;
      if (edge.b < edge.a) std::swap(edge.a, edge.b);
    }
    std::sort(hypothesis.edges.begin(), hypothesis.edges.end(),
              [](const auto &a, const auto &b) {
                return std::pair{a.a, a.b} < std::pair{b.a, b.b};
              });
  }
  std::expected<PbctopoHash128, std::string> selected_hash =
      compatible_hashable.empty()
          ? std::expected<PbctopoHash128, std::string>(PbctopoHash128{})
          : hash_pbctopo_contact_hypotheses(compatible_hashable);
  if (!observed_hash) return std::unexpected(observed_hash.error());
  if (!lost_hash) return std::unexpected(lost_hash.error());
  if (!selected_hash) return std::unexpected(selected_hash.error());
  result.selected_compatible_pair_hash = compatible.empty()
                                             ? PbctopoHash128{}
                                             : *observed_hash;
  result.lost_pair_hash = *lost_hash;
  result.selected_hypothesis_hash = *selected_hash;
  result.selected_compatible_contact_count = observed_pairs.size();
  result.contacts_lost = lost_pairs.size();
  result.selected_compatible_contact_loss_count = lost_pairs.size();
  result.compatible_hypothesis_count = compatible.size();
  result.selected_relation_has_support = !compatible.empty();

  std::vector<unsigned char> seen(selected_assignment.component_ids.size(), 0);
  if (!seen.empty()) {
    std::vector<std::vector<std::size_t>> adjacency(seen.size());
    for (const auto &[a, b] : connected_edges) {
      if (a < seen.size() && b < seen.size()) {
        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
      }
    }
    std::vector<std::size_t> stack{0};
    seen[0] = 1;
    while (!stack.empty()) {
      const auto node = stack.back();
      stack.pop_back();
      for (const auto next : adjacency[node]) {
        if (seen[next] == 0) { seen[next] = 1; stack.push_back(next); }
      }
    }
  }
  result.evidence_graph_connected = !seen.empty() &&
      std::all_of(seen.begin(), seen.end(), [](unsigned char v) { return v != 0; });
  result.evaluated = true;
  result.consistent = lost_pairs.empty();
  if (!result.consistent) {
    result.state = PbctopoEvidenceState::Contradicted;
  } else if (unsupported_relation ||
             hypothesis_set.ambiguous_component_relation_count != 0) {
    result.state = PbctopoEvidenceState::WeakAmbiguous;
  } else if (result.evidence_graph_connected) {
    result.state = PbctopoEvidenceState::SupportedConnected;
  } else {
    result.state = PbctopoEvidenceState::WeakDisconnected;
  }
  return result;
}

} // namespace titan_pbctopo

#endif
