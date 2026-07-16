#ifndef TITAN_PBCTOPO_REPORT_H
#define TITAN_PBCTOPO_REPORT_H

#include <cstddef>
#include <iosfwd>
#include <limits>
#include <string>

#include "types.h"

namespace titan_pbctopo {

enum class PbctopoFrameStatus { Certified, Rescued, Fallback };

enum class PbctopoFrameMode {
  None,
  Virtual,
  EarlyProjection,
  ComponentGraph,
  ComponentSearchExact,
  ComponentSearchBeam,
  Warm,
  HullRescue,
  DistanceRescue,
  Pbci,
  LlpsLocal,
};

enum class PbctopoFailureReason {
  None,
  NoSeed,
  ProjectionInvalid,
  FullVerifyClash,
  ContinuityReject,
  AnchorReject,
  SplitOwnerReject,
  ComponentGraphFailed,
  BeamEmpty,
  NotBetterThanIncumbent,
  CertificationReject,
  FallbackPbci,
};

struct PbctopoScoreVector {
  bool valid = false;
  bool ultra_clash = false;
  std::size_t broken_edges = 0;
  double continuity_max_d2 = std::numeric_limits<double>::quiet_NaN();
  double continuity_path = std::numeric_limits<double>::quiet_NaN();
  double anchor_max_d2 = std::numeric_limits<double>::quiet_NaN();
  double anchor_mst2 = std::numeric_limits<double>::quiet_NaN();
  PbctopoAnchorGraphScore anchor_graph;
  double principal_rg2_max = std::numeric_limits<double>::quiet_NaN();
  double principal_rg2_mid = std::numeric_limits<double>::quiet_NaN();
  double principal_rg2_min = std::numeric_limits<double>::quiet_NaN();
  double shape_kappa2 = std::numeric_limits<double>::quiet_NaN();
  double shape_asphericity = std::numeric_limits<double>::quiet_NaN();
  double shape_acylindricity = std::numeric_limits<double>::quiet_NaN();
  double edge_stretch = std::numeric_limits<double>::quiet_NaN();
  double shift_norm2 = std::numeric_limits<double>::quiet_NaN();
  double gap_score = std::numeric_limits<double>::quiet_NaN();
};

struct PbctopoFrameReport {
  std::size_t frame = 0;
  double time_ps = 0.0;
  PbctopoFrameStatus status = PbctopoFrameStatus::Fallback;
  PbctopoFrameMode mode = PbctopoFrameMode::None;
  std::size_t component_count = 0;
  std::size_t safe_candidate_count = 0;
  PbctopoScoreVector best_score;
  PbctopoScoreVector second_score;
  std::string score_gap_metric = "na";
  double score_gap = std::numeric_limits<double>::quiet_NaN();
  double clash_margin_A = std::numeric_limits<double>::quiet_NaN();
  double continuity_margin_A = std::numeric_limits<double>::quiet_NaN();
  double anchor_max_edge_A = std::numeric_limits<double>::quiet_NaN();
  double split_owner_gap_scaled = std::numeric_limits<double>::quiet_NaN();
  PbctopoAnchorGraphScore anchor_graph;
  PbctopoHardFeasibility hard_feasibility;
  std::uint64_t assignment_signature = 0;
  std::size_t assignment_component_count = 0;
  bool temporal_selected = false;
  bool temporal_changed_from_greedy = false;
  std::uint64_t temporal_assignment_signature = 0;
  bool virtual_search_exhaustive = false;
  bool virtual_branch_and_bound_used = false;
  std::uint64_t virtual_state_space = 0;
  std::size_t virtual_states_evaluated = 0;
  double virtual_omitted_seam_lower_bound =
      std::numeric_limits<double>::quiet_NaN();
  double virtual_lower_bound_gap =
      std::numeric_limits<double>::quiet_NaN();
  PbctopoFailureReason fallback_reason = PbctopoFailureReason::None;
  double time_total_ms = 0.0;
  double time_graph_ms = 0.0;
  double time_virtual_ms = 0.0;
  double time_component_ms = 0.0;
  std::size_t cycle_support = 0;
  double cycle_mismatch = 0.0;
  std::size_t inconsistent_edge_count = 0;
  std::size_t fundamental_cycle_count = 0;
  double fundamental_cycle_residual_rms = 0.0;
  double fundamental_cycle_residual_max = 0.0;
};

const char *pbctopo_frame_status_name(PbctopoFrameStatus status);
const char *pbctopo_frame_mode_name(PbctopoFrameMode mode);
const char *pbctopo_failure_reason_name(PbctopoFailureReason reason);
std::string format_pbctopo_report_double(double value);
std::string format_pbctopo_report_size(std::size_t value);
std::string format_pbctopo_score_vector(const PbctopoScoreVector &score);
void write_pbctopo_frame_report_header(std::ostream &out);
void write_pbctopo_frame_report_row(std::ostream &out,
                                    const PbctopoFrameReport &report);

} // namespace titan_pbctopo

#endif
