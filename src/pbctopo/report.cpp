#include "pbctopo/report.h"

#include <cmath>
#include <ostream>

namespace titan_pbctopo {
namespace {

std::string csv_escape(std::string value) {
  const bool needs_quotes = value.find_first_of(",\"\r\n") != std::string::npos;
  if (!needs_quotes)
    return value;
  std::string escaped = "\"";
  for (char ch : value) {
    if (ch == '"')
      escaped += "\"\"";
    else
      escaped += ch;
  }
  escaped += '"';
  return escaped;
}

} // namespace

const char *pbctopo_frame_status_name(PbctopoFrameStatus status) {
  switch (status) {
  case PbctopoFrameStatus::Certified:
    return "certified";
  case PbctopoFrameStatus::Rescued:
    return "rescued";
  case PbctopoFrameStatus::Fallback:
  default:
    return "fallback";
  }
}

const char *pbctopo_frame_mode_name(PbctopoFrameMode mode) {
  switch (mode) {
  case PbctopoFrameMode::None:
    return "none";
  case PbctopoFrameMode::Virtual:
    return "virtual";
  case PbctopoFrameMode::EarlyProjection:
    return "early_projection";
  case PbctopoFrameMode::ComponentGraph:
    return "exact";
  case PbctopoFrameMode::ComponentSearchExact:
    return "component_exact";
  case PbctopoFrameMode::ComponentSearchBeam:
    return "beam";
  case PbctopoFrameMode::Warm:
    return "warm";
  case PbctopoFrameMode::HullRescue:
    return "hull_rescue";
  case PbctopoFrameMode::DistanceRescue:
    return "distance_rescue";
  case PbctopoFrameMode::Pbci:
  default:
    return "pbci";
  }
}

const char *pbctopo_failure_reason_name(PbctopoFailureReason reason) {
  switch (reason) {
  case PbctopoFailureReason::None:
    return "none";
  case PbctopoFailureReason::NoSeed:
    return "no_seed";
  case PbctopoFailureReason::ProjectionInvalid:
    return "projection_invalid";
  case PbctopoFailureReason::FullVerifyClash:
    return "full_verify_clash";
  case PbctopoFailureReason::ContinuityReject:
    return "continuity_reject";
  case PbctopoFailureReason::AnchorReject:
    return "anchor_reject";
  case PbctopoFailureReason::SplitOwnerReject:
    return "split_owner_reject";
  case PbctopoFailureReason::ComponentGraphFailed:
    return "component_graph_failed";
  case PbctopoFailureReason::BeamEmpty:
    return "beam_empty";
  case PbctopoFailureReason::NotBetterThanIncumbent:
    return "not_better_than_incumbent";
  case PbctopoFailureReason::FallbackPbci:
  default:
    return "fallback_pbci";
  }
}

std::string format_pbctopo_report_double(double value) {
  return std::isfinite(value) ? std::to_string(value) : std::string("na");
}

std::string format_pbctopo_report_size(std::size_t value) {
  return value == std::numeric_limits<std::size_t>::max()
             ? std::string("na")
             : std::to_string(value);
}

std::string format_pbctopo_score_vector(const PbctopoScoreVector &score) {
  return "valid=" + std::string(score.valid ? "1" : "0") +
         ";ultra=" + std::string(score.ultra_clash ? "1" : "0") +
         ";broken=" + format_pbctopo_report_size(score.broken_edges) +
         ";cont_max=" +
         format_pbctopo_report_double(score.continuity_max_d2) +
         ";cont_path=" + format_pbctopo_report_double(score.continuity_path) +
         ";anchor_max=" + format_pbctopo_report_double(score.anchor_max_d2) +
         ";anchor_mst=" + format_pbctopo_report_double(score.anchor_mst2) +
         ";anchor_order_max=" +
         format_pbctopo_report_double(score.anchor_graph.order_max_d2) +
         ";anchor_order_path=" +
         format_pbctopo_report_double(score.anchor_graph.order_path_d2) +
         ";anchor_endpoint_gap=" +
         format_pbctopo_report_double(score.anchor_graph.endpoint_image_gap) +
         ";anchor_residue_violations=" +
         format_pbctopo_report_size(
             score.anchor_graph.residue_order_violations) +
         ";anchor_contact_lost=" +
         format_pbctopo_report_size(score.anchor_graph.contact_map_lost) +
         ";anchor_contact_max=" +
         format_pbctopo_report_double(score.anchor_graph.contact_map_max_d2) +
         ";anchor_ambiguity=" +
         format_pbctopo_report_size(score.anchor_graph.ambiguity_count) +
         ";rg2=" + format_pbctopo_report_double(score.principal_rg2_max) +
         "/" + format_pbctopo_report_double(score.principal_rg2_mid) + "/" +
         format_pbctopo_report_double(score.principal_rg2_min) +
         ";shape=" + format_pbctopo_report_double(score.shape_kappa2) +
         "/" + format_pbctopo_report_double(score.shape_asphericity) + "/" +
         format_pbctopo_report_double(score.shape_acylindricity) +
         ";edge=" + format_pbctopo_report_double(score.edge_stretch) +
         ";shift=" + format_pbctopo_report_double(score.shift_norm2) +
         ";gap=" + format_pbctopo_report_double(score.gap_score);
}

void write_pbctopo_frame_report_header(std::ostream &out) {
  out << "frame,time_ps,status,mode,component_count,safe_candidate_count,"
         "best_score_vector,second_score_vector,score_gap,clash_margin,"
         "continuity_margin,anchor_margin,split_owner_margin,"
         "anchor_order_max_d2,anchor_order_path_d2,"
         "anchor_endpoint_image_gap,anchor_residue_order_violations,"
         "anchor_contact_map_lost,anchor_contact_map_max_d2,"
         "anchor_ambiguity_count,fallback_reason,"
         "time_total_ms,time_graph_ms,time_virtual_ms,time_component_ms,"
         "cycle_support,cycle_mismatch,inconsistent_edge_count\n";
}

void write_pbctopo_frame_report_row(std::ostream &out,
                                    const PbctopoFrameReport &report) {
  out << report.frame << ',' << format_pbctopo_report_double(report.time_ps)
      << ',' << pbctopo_frame_status_name(report.status) << ','
      << pbctopo_frame_mode_name(report.mode) << ',' << report.component_count
      << ',' << report.safe_candidate_count << ','
      << csv_escape(format_pbctopo_score_vector(report.best_score)) << ','
      << csv_escape(format_pbctopo_score_vector(report.second_score)) << ','
      << format_pbctopo_report_double(report.score_gap) << ','
      << format_pbctopo_report_double(report.clash_margin) << ','
      << format_pbctopo_report_double(report.continuity_margin) << ','
      << format_pbctopo_report_double(report.anchor_margin) << ','
      << format_pbctopo_report_double(report.split_owner_margin) << ','
      << format_pbctopo_report_double(report.anchor_graph.order_max_d2) << ','
      << format_pbctopo_report_double(report.anchor_graph.order_path_d2)
      << ','
      << format_pbctopo_report_double(report.anchor_graph.endpoint_image_gap)
      << ','
      << format_pbctopo_report_size(
             report.anchor_graph.residue_order_violations)
      << ',' << format_pbctopo_report_size(report.anchor_graph.contact_map_lost)
      << ','
      << format_pbctopo_report_double(report.anchor_graph.contact_map_max_d2)
      << ',' << format_pbctopo_report_size(report.anchor_graph.ambiguity_count)
      << ',' << pbctopo_failure_reason_name(report.fallback_reason) << ','
      << format_pbctopo_report_double(report.time_total_ms) << ','
      << format_pbctopo_report_double(report.time_graph_ms) << ','
      << format_pbctopo_report_double(report.time_virtual_ms) << ','
      << format_pbctopo_report_double(report.time_component_ms) << ','
      << report.cycle_support << ','
      << format_pbctopo_report_double(report.cycle_mismatch) << ','
      << report.inconsistent_edge_count << '\n';
}

} // namespace titan_pbctopo
