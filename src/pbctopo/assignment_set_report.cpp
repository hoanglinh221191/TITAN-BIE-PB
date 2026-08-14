#include "assignment_set_report.h"

#include "image_assignment.h"
#include "observation_binding.h"
#include "observation_schema.h"
#include "report.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

namespace titan_pbctopo {
namespace {

std::string format_double(double value) {
  if (!std::isfinite(value))
    return "na";
  std::ostringstream out;
  out.imbue(std::locale::classic());
  out << std::setprecision(std::numeric_limits<double>::max_digits10)
      << (value == 0.0 ? 0.0 : value);
  return out.str();
}

std::string observation_class(const PbctopoCutCandidate &candidate,
                               const PbctopoIdentifiability &identifiability,
                               bool selected, bool equivalent_set_complete,
                               bool credible_alternative_set_complete,
                               std::size_t spatial_lift_ambiguities,
                               std::size_t spatial_lift_cycle_residuals) {
  if (!candidate.hard_feasibility.certified())
    return "not_hard_feasible";
  if (!identifiability.identified)
    return "weak_observation";
  if (selected && candidate.search_equivalence_known &&
      candidate.equivalent_to_search_best && equivalent_set_complete &&
      credible_alternative_set_complete &&
      identifiability.objective_uniqueness_known &&
      identifiability.objective_unique_within_search_domain &&
      identifiability.evidence_uniqueness_known &&
      identifiability.evidence_unique_within_search_domain &&
      candidate.search_domain_audit_known &&
       candidate.search_domain_evidence_complete &&
       candidate.component_search_shell_radius > 0 &&
       !candidate.selected_assignment_touches_shell_boundary &&
       spatial_lift_ambiguities == 0 &&
       spatial_lift_cycle_residuals == 0) {
    if (candidate.domain_completeness ==
        PbctopoDomainCompleteness::BoundedExhaustive)
      return "strong_bounded_anchor";
  }
  return "ambiguous_anchor";
}

std::string serialize_component_assignment(
    const PbctopoImageAssignment &assignment) {
  if (!assignment.valid ||
      assignment.component_ids.size() != assignment.component_offsets.size()) {
    return "na";
  }
  struct Entry {
    std::uint64_t component_id = 0;
    Int3 offset{};
  };
  std::vector<Entry> entries;
  entries.reserve(assignment.component_offsets.size());
  for (std::size_t idx = 0; idx < assignment.component_offsets.size(); ++idx) {
    entries.push_back(
        {assignment.component_ids[idx], assignment.component_offsets[idx]});
  }
  std::sort(entries.begin(), entries.end(),
            [](const Entry &lhs, const Entry &rhs) {
              if (lhs.component_id != rhs.component_id)
                return lhs.component_id < rhs.component_id;
              if (lhs.offset.x != rhs.offset.x)
                return lhs.offset.x < rhs.offset.x;
              if (lhs.offset.y != rhs.offset.y)
                return lhs.offset.y < rhs.offset.y;
              return lhs.offset.z < rhs.offset.z;
            });
  std::ostringstream out;
  out.imbue(std::locale::classic());
  for (std::size_t idx = 0; idx < entries.size(); ++idx) {
    if (idx != 0)
      out << '|';
    out << entries[idx].component_id << ':' << entries[idx].offset.x << ':'
        << entries[idx].offset.y << ':' << entries[idx].offset.z;
  }
  return out.str();
}

bool equivalent_under_order(const PbctopoCutCandidate &lhs,
                            const PbctopoCutCandidate &rhs,
                            PbctopoCutCandidateBetterFn better_candidate) {
  return better_candidate != nullptr && !better_candidate(lhs, rhs) &&
         !better_candidate(rhs, lhs);
}

} // namespace

void write_pbctopo_assignment_set_header(std::ostream &out) {
  out << "schema_version,producer_version,coordinate_unit,time_unit,"
         "endianness,floating_format,identity_hash_algorithm,"
         "binding_hash_algorithm,payload_checksum_algorithm,binding_mode,"
         "negative_zero_policy,"
         "nonfinite_policy,observation_scope,topology_epoch_id,"
         "frame,time_ps,candidate_rank,selected_framewise,"
         "candidate_set_scope,observation_class,search_equivalence_known,"
         "equivalent_to_search_best,equivalent_under_output_order,"
         "equivalent_set_complete,retained_certified_assignment_count,"
         "search_equivalent_assignments_exported,"
         "equivalent_best_assignments_within_domain,identified,"
         "unique_within_search_domain,uniqueness_search_exhaustive,"
         "objective_uniqueness_known,"
         "objective_unique_within_search_domain,"
         "evidence_uniqueness_known,"
         "evidence_unique_within_search_domain,"
         "feasible_assignment_set_enumerated,"
         "feasible_assignment_set_semantics,"
         "feasible_assignments_within_domain,"
         "soft_score_valid_assignments_within_domain,"
         "search_domain_assignments_exported,"
         "credible_alternative_set_complete,"
         "component_search_shell_radius,search_domain_audit_known,"
         "search_domain_evidence_complete,"
         "selected_assignment_touches_shell_boundary,"
         "domain_expansion_attempted,domain_expansion_exhausted,"
         "bounded_domain_only,domain_completeness,"
         "hard_feasible,layout_signature,assignment_signature,"
         "evidence_relative_relation_signature,"
         "evidence_compatible_hypothesis_hash,"
         "evidence_compatible_pair_hash,"
         "evidence_compatible_hypothesis_count,"
         "evidence_supported_relation_count,"
         "evidence_compatible_contact_count,"
         "evidence_no_support_relation_count,"
         "component_count,global_gauge_x,global_gauge_y,global_gauge_z,"
         "component_assignment,normalized_gram_g00,normalized_gram_g01,"
         "normalized_gram_g02,normalized_gram_g11,normalized_gram_g12,"
         "normalized_gram_g22,broken_edges,continuity_max_d2,"
         "continuity_path,anchor_max_d2,anchor_mst2,min_pair_d2,"
         "shift_norm2\n";
}

void write_pbctopo_assignment_set_rows(
    std::ostream &out, std::size_t frame, double time_ps,
    std::span<const PbctopoCutCandidate> certified_candidates,
    const PbctopoCutCandidate *framewise_selected,
    const PbctopoIdentifiability &identifiability,
    const std::array<double, 6> &normalized_gram,
    PbctopoCutCandidateBetterFn better_candidate,
    std::uint64_t topology_epoch_id,
    std::size_t soft_score_valid_assignments_within_domain,
    std::size_t spatial_lift_ambiguities,
    std::size_t spatial_lift_cycle_residuals) {
  std::vector<const PbctopoCutCandidate *> ordered;
  ordered.reserve(certified_candidates.size());
  for (const auto &candidate : certified_candidates) {
    if (candidate.valid && candidate.assignment.valid &&
        candidate.hard_feasibility.certified()) {
      ordered.push_back(&candidate);
    }
  }
  if (better_candidate != nullptr) {
    std::stable_sort(ordered.begin(), ordered.end(),
                     [better_candidate](const auto *lhs, const auto *rhs) {
                       return better_candidate(*lhs, *rhs);
                     });
  }
  if (ordered.empty())
    return;

  const PbctopoCutCandidate &pool_best = *ordered.front();
  std::size_t search_equivalent_exported = 0;
  std::size_t search_domain_assignments_exported = 0;
  for (const auto *candidate : ordered) {
    if (candidate->search_equivalence_known) {
      ++search_domain_assignments_exported;
    }
    if (candidate->search_equivalence_known &&
        candidate->equivalent_to_search_best) {
      ++search_equivalent_exported;
    }
  }
  const bool equivalent_set_complete =
      identifiability.uniqueness_search_exhaustive &&
      identifiability.equivalent_best_assignments_within_domain > 0 &&
      search_equivalent_exported ==
          identifiability.equivalent_best_assignments_within_domain;
  const bool all_retained_assignments_have_search_provenance =
      search_domain_assignments_exported == ordered.size();
  const bool credible_alternative_set_complete =
      identifiability.feasible_assignment_set_enumerated &&
      identifiability.feasible_assignments_within_domain > 0 &&
      all_retained_assignments_have_search_provenance &&
      search_domain_assignments_exported ==
          identifiability.feasible_assignments_within_domain;

  for (std::size_t rank = 0; rank < ordered.size(); ++rank) {
    const PbctopoCutCandidate &candidate = *ordered[rank];
    const bool selected =
        framewise_selected != nullptr && framewise_selected->assignment.valid &&
        same_pbctopo_image_assignment(candidate.assignment,
                                       framewise_selected->assignment);
    const bool equivalent_under_output_order =
        equivalent_under_order(pool_best, candidate, better_candidate);
    out << VIBE_OBSERVATION_SCHEMA_VERSION << ','
        << vibe_observation_producer_version() << ','
        << VIBE_OBSERVATION_COORDINATE_UNIT << ','
        << VIBE_OBSERVATION_TIME_UNIT << ','
        << VIBE_OBSERVATION_ENDIANNESS << ','
        << VIBE_OBSERVATION_FLOATING_FORMAT << ','
        << VIBE_OBSERVATION_IDENTITY_HASH_ALGORITHM << ','
        << VIBE_OBSERVATION_BINDING_HASH_ALGORITHM << ','
        << VIBE_OBSERVATION_PAYLOAD_CHECKSUM_ALGORITHM << ','
        << VIBE_OBSERVATION_BINDING_MODE << ','
        << VIBE_OBSERVATION_NEGATIVE_ZERO_POLICY << ','
        << VIBE_OBSERVATION_NONFINITE_POLICY << ','
        << VIBE_OBSERVATION_SCOPE << ','
        << topology_epoch_id << ',' << frame << ',' << format_double(time_ps)
        << ',' << rank + 1 << ','
        << (selected ? 1 : 0)
        << ",retained_certified_band,"
        << observation_class(candidate, identifiability, selected,
                             equivalent_set_complete,
                             credible_alternative_set_complete,
                             spatial_lift_ambiguities,
                             spatial_lift_cycle_residuals)
        << ','
        << (candidate.search_equivalence_known ? 1 : 0) << ',';
    if (candidate.search_equivalence_known) {
      out << (candidate.equivalent_to_search_best ? 1 : 0);
    } else {
      out << "na";
    }
    out << ',' << (equivalent_under_output_order ? 1 : 0) << ','
        << (equivalent_set_complete ? 1 : 0) << ',' << ordered.size() << ','
        << search_equivalent_exported << ','
        << identifiability.equivalent_best_assignments_within_domain << ','
        << (identifiability.identified ? 1 : 0) << ','
        << (identifiability.unique_within_search_domain ? 1 : 0) << ','
        << (identifiability.uniqueness_search_exhaustive ? 1 : 0) << ','
        << (identifiability.objective_uniqueness_known ? 1 : 0) << ','
        << (identifiability.objective_unique_within_search_domain ? 1 : 0)
        << ',' << (identifiability.evidence_uniqueness_known ? 1 : 0) << ','
        << (identifiability.evidence_unique_within_search_domain ? 1 : 0)
        << ','
        << (identifiability.feasible_assignment_set_enumerated ? 1 : 0)
        << ','
        << (identifiability.feasible_assignment_set_enumerated
                ? "hard_constraints_within_search_domain"
                : "not_enumerated")
        << ',' << identifiability.feasible_assignments_within_domain << ','
        << soft_score_valid_assignments_within_domain << ','
        << search_domain_assignments_exported << ','
        << (credible_alternative_set_complete ? 1 : 0) << ','
        << candidate.component_search_shell_radius << ','
        << (candidate.search_domain_audit_known ? 1 : 0) << ','
        << (candidate.search_domain_evidence_complete ? 1 : 0) << ',';
    if (candidate.search_domain_audit_known)
      out << (candidate.selected_assignment_touches_shell_boundary ? 1 : 0);
    else
      out << "na";
    out << ',' << (candidate.domain_expansion_attempted ? 1 : 0) << ','
        << (candidate.domain_expansion_exhausted ? 1 : 0) << ','
        << (candidate.bounded_domain_only ? 1 : 0) << ','
        << pbctopo_domain_completeness_name(candidate.domain_completeness)
        << ','
        << (candidate.hard_feasibility.certified() ? 1 : 0) << ','
        << candidate.assignment.layout_signature << ','
        << candidate.assignment.signature << ','
        << candidate.evidence_relative_relation_signature << ','
        << format_pbctopo_hash128(
               candidate.evidence_compatible_hypothesis_hash) << ','
        << format_pbctopo_hash128(candidate.evidence_compatible_pair_hash)
        << ',' << candidate.evidence_compatible_hypothesis_count << ','
        << candidate.evidence_supported_relation_count << ','
        << candidate.evidence_compatible_contact_count << ','
        << candidate.evidence_no_support_relation_count << ','
        << candidate.assignment.component_offsets.size() << ','
        << candidate.assignment.global_gauge.x << ','
        << candidate.assignment.global_gauge.y << ','
        << candidate.assignment.global_gauge.z << ','
        << serialize_component_assignment(candidate.assignment);
    for (double value : normalized_gram)
      out << ',' << format_double(value);
    out << ',' << candidate.broken_edges << ','
        << format_double(candidate.continuity_max_d2) << ','
        << format_double(candidate.continuity_path) << ','
        << format_double(candidate.anchor_max_d2) << ','
        << format_double(candidate.anchor_mst2) << ','
        << format_double(candidate.min_pair_d2) << ','
        << format_double(candidate.shift_norm2) << '\n';
  }
}

} // namespace titan_pbctopo
