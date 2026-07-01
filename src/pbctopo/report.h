#ifndef TITAN_PBCTOPO_REPORT_H
#define TITAN_PBCTOPO_REPORT_H

#include <cstddef>
#include <iosfwd>
#include <limits>
#include <string>

namespace titan_pbctopo {

struct PbctopoAnchorGraphScore {
  double order_max_d2 = std::numeric_limits<double>::infinity();
  double order_path_d2 = std::numeric_limits<double>::infinity();
  double endpoint_image_gap = std::numeric_limits<double>::infinity();
  std::size_t residue_order_violations = 0;
  std::size_t contact_map_lost = 0;
  double contact_map_max_d2 = std::numeric_limits<double>::infinity();
  std::size_t ambiguity_count = 0;
};

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
  double score_gap = std::numeric_limits<double>::quiet_NaN();
  double clash_margin = std::numeric_limits<double>::quiet_NaN();
  double continuity_margin = std::numeric_limits<double>::quiet_NaN();
  double anchor_margin = std::numeric_limits<double>::quiet_NaN();
  double split_owner_margin = std::numeric_limits<double>::quiet_NaN();
  PbctopoAnchorGraphScore anchor_graph;
  PbctopoFailureReason fallback_reason = PbctopoFailureReason::None;
  double time_total_ms = 0.0;
  double time_graph_ms = 0.0;
  double time_virtual_ms = 0.0;
  double time_component_ms = 0.0;
  std::size_t cycle_support = 0;
  double cycle_mismatch = 0.0;
  std::size_t inconsistent_edge_count = 0;
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
