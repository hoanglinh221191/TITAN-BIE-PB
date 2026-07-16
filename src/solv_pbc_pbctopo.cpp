#include "Multi_MOL.h"
#include "BoxMatrix.h"
#include "aligned_allocator.h"
#include "pbctopo/config.h"
#include "pbctopo/anchor_graph.h"
#include "pbctopo/candidate_order.h"
#include "pbctopo/candidate_scoring.h"
#include "pbctopo/certification.h"
#include "pbctopo/component_graph.h"
#include "pbctopo/component_metric.h"
#include "pbctopo/component_refine.h"
#include "pbctopo/component_search.h"
#include "pbctopo/component_shift_domain.h"
#include "pbctopo/llps_local.h"
#include "pbctopo/image_assignment.h"
#include "pbctopo/lattice_math.h"
#include "pbctopo/local_graph.h"
#include "pbctopo/logging.h"
#include "pbctopo/report.h"
#include "pbctopo/steric_clash.h"
#include "pbctopo/temporal_path.h"
#include "pbctopo/temporal_replay.h"
#include "pbctopo/types.h"
#include "pbctopo/virtual_seed.h"
#include "pbctopo/virtual_discrete_search.h"
#include "pbctopo_hull.h"
#include "titan_gyration.h"
#include "titan_simd_dispatch.h"
#include "titan_executor.h"
#include "titan_ui.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

namespace {

using pbctopo_hull::cross;
using pbctopo_hull::dot;
using pbctopo_hull::norm2;
using pbctopo_hull::PbctopoHullScratch;
using pbctopo_hull::PbctopoHullSpan;
using pbctopo_hull::Vec3;
using titan_gyration::build_point_stats;
using titan_gyration::compute_bbox_volume;
using titan_gyration::compute_gyration_shape_metrics;
using titan_gyration::compute_principal_rg2;
using titan_gyration::compute_rg2;
using titan_gyration::compute_rg2_from_stats;
using titan_gyration::GyrationShapeMetrics;
using titan_gyration::PointStats;
using PbctopoAnchorBackend = titan_pbctopo::AnchorBackend;
using PbctopoAnchorMetricMode = titan_pbctopo::AnchorMetricMode;
using PbctopoRuntimeCutoffs = titan_pbctopo::RuntimeCutoffs;
using PbctopoSolveMode = titan_pbctopo::SolveMode;
using PbctopoTemporalMode = titan_pbctopo::TemporalMode;
using PbctopoTemporalCandidate = titan_pbctopo::PbctopoTemporalCandidate;
using PbctopoAnchorContactEdge = titan_pbctopo::PbctopoAnchorContactEdge;
using PbctopoAnchorGraphContext = titan_pbctopo::PbctopoAnchorGraphContext;
using PbctopoAnchorGraphScore = titan_pbctopo::PbctopoAnchorGraphScore;
using PbctopoAtom = titan_pbctopo::PbctopoAtom;
using PbctopoComponentCandidate = titan_pbctopo::PbctopoComponentCandidate;
using PbctopoComponentCandidatePrefixStage =
    titan_pbctopo::PbctopoComponentCandidatePrefixStage;
using PbctopoComponentBaseGeometry = titan_pbctopo::PbctopoComponentBaseGeometry;
using PbctopoComponentSearchCallbacks =
    titan_pbctopo::PbctopoComponentSearchCallbacks;
using PbctopoContinuityPair = titan_pbctopo::PbctopoContinuityPair;
using PbctopoCutScoreBreakdown = titan_pbctopo::PbctopoCutScoreBreakdown;
using PbctopoCutCandidate = titan_pbctopo::PbctopoCutCandidate;
using PbctopoCycleSyncMetrics = titan_pbctopo::PbctopoCycleSyncMetrics;
using PbctopoCycleSyncResult = titan_pbctopo::PbctopoCycleSyncResult;
using PbctopoEdgeScanSoA = titan_pbctopo::PbctopoEdgeScanSoA;
using PbctopoFailureReason = titan_pbctopo::PbctopoFailureReason;
using PbctopoFrameMode = titan_pbctopo::PbctopoFrameMode;
using PbctopoFrameReport = titan_pbctopo::PbctopoFrameReport;
using PbctopoFrameStatus = titan_pbctopo::PbctopoFrameStatus;
using PbctopoGraphEdge = titan_pbctopo::PbctopoGraphEdge;
using PbctopoHardFeasibility = titan_pbctopo::PbctopoHardFeasibility;
using PbctopoImageAssignment = titan_pbctopo::PbctopoImageAssignment;
using PbctopoImageEquationResiduals =
    titan_pbctopo::PbctopoImageEquationResiduals;
using PbctopoLocalEdge = titan_pbctopo::PbctopoLocalEdge;
using PbctopoLocalEdgeKind = titan_pbctopo::PbctopoLocalEdgeKind;
using PbctopoLlpsLocalUnits = titan_pbctopo::PbctopoLlpsLocalUnits;
using PbctopoPairShiftCandidate = titan_pbctopo::PbctopoPairShiftCandidate;
using PbctopoPairShiftCandidateSet = titan_pbctopo::PbctopoPairShiftCandidateSet;
using PbctopoPairShiftScore = titan_pbctopo::PbctopoPairShiftScore;
using PbctopoRefineScorePhase = titan_pbctopo::PbctopoRefineScorePhase;
using PbctopoRefineScoringCallbacks =
    titan_pbctopo::PbctopoRefineScoringCallbacks;
using PbctopoScoreVector = titan_pbctopo::PbctopoScoreVector;
using PbctopoShift = titan_pbctopo::PbctopoShift;
using PbctopoShiftedComponent = titan_pbctopo::PbctopoShiftedComponent;
using PbctopoShiftedComponentDescriptor =
    titan_pbctopo::PbctopoShiftedComponentDescriptor;
using PbctopoStage1EdgeContext = titan_pbctopo::PbctopoStage1EdgeContext;
using PbctopoRange = titan_pbctopo::PbctopoRange;
using PbctopoSystemStaticData = titan_pbctopo::PbctopoSystemStaticData;
using PbctopoTrajectoryContext = titan_pbctopo::PbctopoTrajectoryContext;
using PbctopoTimingStats = titan_pbctopo::PbctopoTimingStats;
using PbctopoVirtualSeedCenter = titan_pbctopo::PbctopoVirtualSeedCenter;
using PbctopoVirtualSeedSource = titan_pbctopo::PbctopoVirtualSeedSource;
using PbctopoLatticeMetric = titan_pbctopo::PbctopoLatticeMetric;
using PbctopoLexicographicGap = titan_pbctopo::PbctopoLexicographicGap;
using PbctopoVirtualDiscreteSearchOptions =
    titan_pbctopo::PbctopoVirtualDiscreteSearchOptions;
using PbctopoVirtualDiscreteSearchPlan =
    titan_pbctopo::PbctopoVirtualDiscreteSearchPlan;
using PbctopoVirtualDiscreteSeedHint =
    titan_pbctopo::PbctopoVirtualDiscreteSeedHint;
using PbctopoWarmState = titan_pbctopo::PbctopoWarmState;
using ClashCellKey = titan_pbctopo::ClashCellKey;
using ClashCellKeyHash = titan_pbctopo::ClashCellKeyHash;
using Int3 = titan_pbctopo::Int3;
using titan_pbctopo::is_zero_delta;
using titan_pbctopo::same_delta;
using titan_pbctopo::anchor_graph_metrics_enhanced;
using titan_pbctopo::build_anchor_graph_context;
using titan_pbctopo::build_atom_component_index;
using titan_pbctopo::build_coarse_anchor_subset;
using titan_pbctopo::build_component_graph_edge_candidates;
using titan_pbctopo::build_component_points_from_offsets;
using titan_pbctopo::build_edge_scan_soa;
using titan_pbctopo::build_stage1_edge_context;
using titan_pbctopo::build_tapered_stage1_anchor_subset;
using titan_pbctopo::build_third_anchor_subset;
using titan_pbctopo::clear_component_trial_hints;
using titan_pbctopo::compare_anchor_graph_score;
using titan_pbctopo::compare_pbctopo_metric;
using titan_pbctopo::component_trial_hints_compatible;
using titan_pbctopo::better_component_graph_edge;
using titan_pbctopo::better_component_candidate;
using titan_pbctopo::better_cut_candidate;
using titan_pbctopo::better_preclash_cut_candidate;
using titan_pbctopo::better_stage1_cut_candidate;
using titan_pbctopo::candidate_anchor_prefix_equivalent;
using titan_pbctopo::candidate_close_enough_to_skip_search;
using titan_pbctopo::candidate_has_projection_clash_margin;
using titan_pbctopo::candidate_lexicographic_gap;
using titan_pbctopo::candidate_skip_search_reject_reason;
using titan_pbctopo::component_candidate_cannot_beat_incumbent_prefix;
using titan_pbctopo::component_candidate_is_outlier;
using titan_pbctopo::component_candidate_to_cut_candidate;
using titan_pbctopo::compare_candidate_prefix_before_edge_stretch;
using titan_pbctopo::compare_pbctopo_double_metric;
using titan_pbctopo::compare_principal_rg2;
using titan_pbctopo::compare_shape_kappa;
using titan_pbctopo::compare_shape_moment;
using titan_pbctopo::cut_candidate_cannot_beat_incumbent;
using titan_pbctopo::make_score_vector;
using titan_pbctopo::preclash_cut_candidate_cannot_beat_incumbent;
using titan_pbctopo::graph_edge_has_zero_shift;
using titan_pbctopo::graph_edge_matches_delta;
using titan_pbctopo::append_atoms_for_system;
using titan_pbctopo::append_metadata_graph_edges_for_system;
using titan_pbctopo::append_topology_bond_edges_for_system;
using titan_pbctopo::build_pbctopo_llps_local_units;
using titan_pbctopo::build_pbctopo_trajectory_context;
using titan_pbctopo::improve_component_offsets_locally;
using titan_pbctopo::make_pbctopo_image_assignment;
using titan_pbctopo::build_heavy_atom_subset;
using titan_pbctopo::build_local_neighborhood_graph;
using titan_pbctopo::build_projection_anchor_subset;
using titan_pbctopo::build_projection_continuity_subset;
using titan_pbctopo::build_projection_edge_subset;
using titan_pbctopo::build_ranges_for_system;
using titan_pbctopo::merge_pbctopo_llps_unit_source;
using titan_pbctopo::pbctopo_topology_source_label;
using titan_pbctopo::pbctopo_atoms_clash;
using titan_pbctopo::pbctopo_max_clash_cutoff;
using titan_pbctopo::pbctopo_component_shift_metric;
using titan_pbctopo::prepare_pbctopo_steric_radii;
using titan_pbctopo::certify_pbctopo_candidate;
using titan_pbctopo::evaluate_pbctopo_image_equations;
using titan_pbctopo::rebase_component_offsets;
using titan_pbctopo::replay_pbctopo_assignment_delta;
using titan_pbctopo::same_component_layout;
using titan_pbctopo::same_pbctopo_image_assignment;
using titan_pbctopo::score_component_offset_solution;
using titan_pbctopo::score_component_placed_points_solution;
using titan_pbctopo::search_component_shift_candidate;
using titan_pbctopo::search_component_shift_candidate_beam;
using titan_pbctopo::kPbctopoExactAssignmentBudget;
using titan_pbctopo::pbctopo_component_assignment_upper_bound;
using titan_pbctopo::solve_component_offsets_from_ranked_edges;
using titan_pbctopo::solve_component_offsets_via_pair_graph;
using titan_pbctopo::unwrap_pbctopo_llps_local_points;
using titan_pbctopo::write_pbctopo_frame_report_header;
using titan_pbctopo::write_pbctopo_frame_report_row;
using titan_pbctopo::extrapolate_virtual_seed_center;
using titan_pbctopo::build_pbctopo_virtual_discrete_search_plan;

constexpr double kDistanceEps = 1.0e-12;
constexpr int kPackedCellBias = 1 << 20;
constexpr uint64_t kPackedCellMask = (1ull << 21) - 1ull;
constexpr size_t kPbctopoDenseClashGridMaxCells = 1u << 24;
constexpr size_t kPbctopoExactComponentLimit =
    titan_pbctopo::kDefaultExactComponentLimit;
constexpr size_t kPbctopoMergeModeMaxRawComponents = 16;
constexpr size_t kPbctopoMergedNodeMaxAtoms = 4096;
constexpr size_t kPbctopoMergedNodeMaxRawComponents = 4;
constexpr size_t kPbctopoBeamWidth = titan_pbctopo::kDefaultBeamWidth;
constexpr size_t kPbctopoExactComponentLimitMax =
    titan_pbctopo::kMaxExactComponentLimit;
constexpr size_t kPbctopoBeamWidthMax = titan_pbctopo::kMaxBeamWidth;

void reset_pbctopo_log_file() {
  titan_pbctopo::reset_log_file();
}

void append_pbctopo_log_line(const string &message) {
  titan_pbctopo::append_log_line(message);
}

ofstream open_pbctopo_frame_report() {
  namespace fs = std::filesystem;
  static bool initialized = false;
  static fs::path initialized_cwd;

  std::error_code ec;
  const fs::path cwd = fs::current_path(ec);
  const bool reset_report =
      !initialized || ec || initialized_cwd.empty() || initialized_cwd != cwd;

  ofstream out("pbctopo_frame_report.csv",
               ios::out | (reset_report ? ios::trunc : ios::app));
  if (out.is_open() && reset_report)
    write_pbctopo_frame_report_header(out);

  initialized = true;
  initialized_cwd = cwd;
  return out;
}

struct PbctopoFragment {
  size_t id = 0;
  int owner = 0;
  size_t source_index = 0;
  size_t begin_local = 0;
  size_t end_local = 0;
  size_t base_atom_begin = 0;
  size_t base_atom_end = 0;
  size_t atom_count = 0;
  bool topo_backed = false;
};

struct PbctopoAtomCopy {
  size_t global_atom_id = 0;
  size_t fragment_id = 0;
  size_t node_id = 0;
  int owner = 0;
  size_t local_atom_index = 0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct PbctopoFragmentCopyNode {
  size_t node_id = 0;
  size_t fragment_id = 0;
  size_t shift_index = 26;
  int shift_ix = 0;
  int shift_iy = 0;
  int shift_iz = 0;
  double tx = 0.0;
  double ty = 0.0;
  double tz = 0.0;
  size_t atom_begin = 0;
  size_t atom_end = 0;
  double min_x = 0.0;
  double min_y = 0.0;
  double min_z = 0.0;
  double max_x = 0.0;
  double max_y = 0.0;
  double max_z = 0.0;
};

struct PbctopoOwnerPairSplitDetail {
  int owner_a = -1;
  int owner_b = -1;
  Vec3 raw_delta_scaled{};
  Vec3 minimg_delta_scaled{};
  double raw_norm = 0.0;
  double minimg_norm = 0.0;
  bool valid = false;
};

struct PbctopoTranslationCache {
  Vec3 by_shift[5][5][5]{};
};

struct PbctopoCandidate {
  vector<size_t> node_ids;
  double hull_volume = numeric_limits<double>::infinity();
  double bbox_volume = numeric_limits<double>::infinity();
  double rg2 = numeric_limits<double>::infinity();
  double shape_kappa2 = numeric_limits<double>::infinity();
  double shape_asphericity = numeric_limits<double>::infinity();
  double shape_acylindricity = numeric_limits<double>::infinity();
  double fragment_mst2 = numeric_limits<double>::infinity();
  double fragment_rg2 = numeric_limits<double>::infinity();
  double fragment_kappa2 = numeric_limits<double>::infinity();
  double fragment_asphericity = numeric_limits<double>::infinity();
  double fragment_acylindricity = numeric_limits<double>::infinity();
  double shift_norm2 = numeric_limits<double>::infinity();
  double min_interfragment_d2 = numeric_limits<double>::infinity();
  size_t connectivity_edges = 0;
  bool hull_valid = false;
  bool hard_clash = false;
  bool valid = false;
};

struct PbctopoCutOption {
  double cut = 0.0;
  double gap = 0.0;
};

struct PbctopoSeed2D {
  double cx = 1.5;
  double cy = 1.5;
  double score = -numeric_limits<double>::infinity();
  size_t ix = 0;
  size_t iy = 0;
};

struct PbctopoSeed3D {
  double cx = 1.5;
  double cy = 1.5;
  double cz = 1.5;
  double score = -numeric_limits<double>::infinity();
  size_t ix = 0;
  size_t iy = 0;
  size_t iz = 0;
  PbctopoVirtualSeedSource source = PbctopoVirtualSeedSource::Density;
};

using PbctopoGyrationShapeMetrics = GyrationShapeMetrics;

string format_int3(const Int3 &value) {
  return "(" + to_string(value.x) + "," + to_string(value.y) + "," +
         to_string(value.z) + ")";
}

string format_vec3(const Vec3 &value) {
  return "(" + to_string(value.x) + "," + to_string(value.y) + "," +
         to_string(value.z) + ")";
}

string format_component_offsets(const vector<Int3> &offsets) {
  string formatted = "[";
  for (size_t idx = 0; idx < offsets.size(); ++idx) {
    if (idx > 0)
      formatted += " ";
    formatted += to_string(idx) + ":" + format_int3(offsets[idx]);
  }
  formatted += "]";
  return formatted;
}

string format_graph_edges(const vector<PbctopoGraphEdge> &edges) {
  string formatted = "[";
  for (size_t idx = 0; idx < edges.size(); ++idx) {
    if (idx > 0)
      formatted += " ";
    const auto &edge = edges[idx];
    formatted += to_string(edge.a) + "-" + to_string(edge.b) + ":" +
                 format_int3({edge.dx, edge.dy, edge.dz}) + "/min=" +
                 to_string(edge.min_pair_d2) + "/bbox=" +
                 to_string(edge.bbox_d2) + "/conf=" +
                 to_string(edge.confidence);
  }
  formatted += "]";
  return formatted;
}

enum class PbctopoComponentMode {
  Disabled,
  Exact,
  MergeThenGraph,
  Beam,
};

struct PbctopoMergedComponentPlan {
  vector<vector<size_t>> merged_components;
  vector<vector<size_t>> merged_to_raw;
  vector<size_t> raw_to_merged;
  vector<Int3> raw_internal_offsets;
};

struct PbctopoGraphSeedQuickScore {
  bool valid = false;
  double continuity_max_d2 = numeric_limits<double>::infinity();
  double continuity_path = numeric_limits<double>::infinity();
  double anchor_max_d2 = numeric_limits<double>::infinity();
  double anchor_mst2 = numeric_limits<double>::infinity();
  PbctopoAnchorGraphScore anchor_graph;
  double principal_rg2_max = numeric_limits<double>::infinity();
  double principal_rg2_mid = numeric_limits<double>::infinity();
  double principal_rg2_min = numeric_limits<double>::infinity();
  double shift_norm2 = numeric_limits<double>::infinity();
};

struct PbctopoCoarseSeedOption {
  uint64_t state_key = 0;
  PbctopoSeed3D seed;
};

struct PbctopoFrameScratch {
  vector<PbctopoAtom> base_atoms;
  vector<PbctopoAtom> projection_atoms;
  vector<PbctopoAtom> component_atoms;
  vector<PbctopoLocalEdge> local_edges;
  vector<PbctopoLocalEdge> component_scoring_edges;
  vector<PbctopoLocalEdge> projection_edges;
  vector<PbctopoContinuityPair> projection_continuity_pairs;
  vector<size_t> projection_anchor_indices;
  unordered_set<uint64_t> local_edge_set;
  unordered_set<uint64_t> projection_edge_set;
  vector<vector<size_t>> components;
  vector<size_t> atom_component_index;
  vector<size_t> component_sizes;
  vector<PbctopoShift> shifts;
  vector<PbctopoSeed3D> volume_seeds;
  vector<Vec3> warm_points;
  vector<Vec3> coarse_warm_points;
  vector<Vec3> placed_points;
  vector<Vec3> best_safe_points;
  vector<Vec3> verify_points;
  vector<Vec3> component_points;
  vector<PbctopoCutCandidate> coarse_candidates;
  vector<PbctopoCoarseSeedOption> coarse_seed_options;
  vector<PbctopoCutCandidate> coarse_results;
  vector<size_t> projection_index;
  titan_pbctopo::PbctopoVirtualDiscreteSearchWorkspace virtual_search_workspace;
  PbctopoVirtualDiscreteSearchPlan virtual_search_plan;

  void reserve_from_context(const PbctopoTrajectoryContext &context) {
    base_atoms.reserve(context.total_atoms_hint);
    projection_atoms.reserve(context.heavy_atoms_hint);
    component_atoms.reserve(context.total_atoms_hint);
    local_edges.reserve(context.local_edges_hint);
    component_scoring_edges.reserve(context.local_edges_hint + 32);
    projection_edges.reserve(context.projection_edges_hint);
    projection_continuity_pairs.reserve(context.continuity_pair_hint);
    projection_anchor_indices.reserve(context.continuity_pair_hint);
    local_edge_set.reserve(context.local_edges_hint * 2 + 1);
    projection_edge_set.reserve(context.projection_edges_hint * 2 + 1);
    components.reserve(context.component_hint);
    atom_component_index.reserve(context.total_atoms_hint);
    component_sizes.reserve(context.component_hint);
    shifts.reserve(27);
    volume_seeds.reserve(context.coarse_seed_hint);
    warm_points.reserve(context.total_atoms_hint);
    coarse_warm_points.reserve(context.heavy_atoms_hint);
    placed_points.reserve(context.heavy_atoms_hint);
    best_safe_points.reserve(context.total_atoms_hint);
    verify_points.reserve(context.total_atoms_hint);
    component_points.reserve(context.total_atoms_hint);
    coarse_candidates.reserve(context.coarse_candidate_hint);
    coarse_seed_options.reserve(context.coarse_seed_hint);
    coarse_results.reserve(context.coarse_seed_hint);
    projection_index.reserve(context.total_atoms_hint);
  }

  void reset_frame() {
    base_atoms.clear();
    projection_atoms.clear();
    component_atoms.clear();
    local_edges.clear();
    component_scoring_edges.clear();
    projection_edges.clear();
    projection_continuity_pairs.clear();
    projection_anchor_indices.clear();
    local_edge_set.clear();
    projection_edge_set.clear();
    components.clear();
    atom_component_index.clear();
    component_sizes.clear();
    shifts.clear();
    volume_seeds.clear();
    warm_points.clear();
    coarse_warm_points.clear();
    placed_points.clear();
    best_safe_points.clear();
    verify_points.clear();
    component_points.clear();
    coarse_candidates.clear();
    coarse_seed_options.clear();
    coarse_results.clear();
    projection_index.clear();
    virtual_search_plan.x_states.clear();
    virtual_search_plan.y_states.clear();
    virtual_search_plan.z_states.clear();
    virtual_search_plan.states.clear();
  }

  size_t point_capacity_bytes_estimate() const noexcept {
    return (warm_points.capacity() + coarse_warm_points.capacity() +
            placed_points.capacity() + best_safe_points.capacity() +
            verify_points.capacity() + component_points.capacity()) *
           sizeof(Vec3);
  }
};

struct PbctopoUltraclashScratch {
  vector<int> cell_x;
  vector<int> cell_y;
  vector<int> cell_z;
  vector<int> head;
  vector<uint32_t> head_epoch;
  vector<int> next;
  uint32_t current_epoch = 0;
};

uint64_t pack_cell_key(int x, int y, int z) noexcept {
  // Current pbctopo workloads stay comfortably inside +/- 2^20 cells.
  const uint64_t ux =
      static_cast<uint64_t>(static_cast<uint32_t>(x + kPackedCellBias)) &
      kPackedCellMask;
  const uint64_t uy =
      static_cast<uint64_t>(static_cast<uint32_t>(y + kPackedCellBias)) &
      kPackedCellMask;
  const uint64_t uz =
      static_cast<uint64_t>(static_cast<uint32_t>(z + kPackedCellBias)) &
      kPackedCellMask;
  return (ux << 42) | (uy << 21) | uz;
}

uint64_t pack_cell_key(const ClashCellKey &key) noexcept {
  return pack_cell_key(key.x, key.y, key.z);
}

struct PersistentClashGrid {
  unordered_map<ClashCellKey, vector<size_t>, ClashCellKeyHash> cells;
  vector<ClashCellKey> atom_keys;
  vector<size_t> atom_offsets;
  double inv_grid = 1.0;
};

struct NodeVectorHash {
  size_t operator()(const vector<size_t> &values) const noexcept {
    uint64_t hash = 1469598103934665603ull;
    for (size_t value : values) {
      hash ^= static_cast<uint64_t>(value) + 0x9e3779b97f4a7c15ull +
              (hash << 6) + (hash >> 2);
      hash *= 1099511628211ull;
    }
    return static_cast<size_t>(hash);
  }
};

double compute_mst2(const vector<Vec3> &points);
vector<size_t> collect_fragment_options(
    size_t fragment_id, const vector<size_t> &selected_nodes,
    const vector<PbctopoFragmentCopyNode> &nodes,
    const vector<vector<size_t>> &adjacency, size_t max_options);
vector<pair<int, int>>
build_contacting_owner_pairs(const vector<PbctopoLocalEdge> &local_edges,
                             const vector<PbctopoAtom> &atoms, box Box,
                             double contact_cutoff,
                             bool hard_required_only = false);
unordered_map<int, Vec3>
compute_owner_com_cartesian(const vector<Vec3> &placed_points,
                            const vector<PbctopoAtom> &atoms);
unordered_map<int, Vec3>
compute_owner_com_scaled(box Box,
                         const unordered_map<int, Vec3> &owner_com_cartesian);
bool candidate_has_split_contacting_owners(
    const vector<Vec3> &placed_points, const vector<PbctopoAtom> &atoms, box Box,
    const vector<pair<int, int>> &contacting_owner_pairs, size_t &split_owner_pairs,
    double &max_owner_raw_minus_minimg_delta,
    PbctopoOwnerPairSplitDetail *detail = nullptr);
bool reimage_contacting_owner_points_minimum_image(
    vector<Vec3> &placed_points, const vector<PbctopoAtom> &atoms, box Box,
    const vector<pair<int, int>> &contacting_owner_pairs);
bool rescue_split_contacting_owner_candidate(
    const vector<Vec3> &seed_points, const vector<PbctopoAtom> &atoms, box Box,
    const vector<pair<int, int>> &contacting_owner_pairs,
    const vector<PbctopoShift> &shifts,
    const PbctopoOwnerPairSplitDetail &split_detail, vector<Vec3> &rescued_points,
    Int3 &applied_shift, double &best_raw_delta_norm);
PbctopoCutCandidate score_cut_candidate(const vector<PbctopoLocalEdge> &local_edges,
                                        const vector<Vec3> &placed_points,
                                        double cut_x, double cut_y, double cut_z,
                                        double gap_score, const PointStats &stats,
                                        double shift_norm2,
                                        const PbctopoStage1EdgeContext *edge_context = nullptr,
                                        const PbctopoEdgeScanSoA *edge_scan_soa = nullptr,
                                        const PbctopoCutCandidate *edge_prune_threshold = nullptr,
                                        const PbctopoCutCandidate *pre_edge_seed = nullptr,
                                        PbctopoRefineScorePhase phase =
                                            PbctopoRefineScorePhase::Exact);
PbctopoCutCandidate score_cut_candidate_principal_only(
    const vector<PbctopoLocalEdge> &local_edges,
    const vector<Vec3> &placed_points, double cut_x, double cut_y, double cut_z,
    double gap_score, const PointStats &stats, double shift_norm2,
    PbctopoCutScoreBreakdown *breakdown = nullptr,
    const PbctopoStage1EdgeContext *edge_context = nullptr,
    const PbctopoEdgeScanSoA *edge_scan_soa = nullptr,
    const PbctopoCutCandidate *edge_prune_threshold = nullptr,
    const PbctopoCutCandidate *pre_edge_seed = nullptr);
bool edge_scan_cannot_beat_threshold(const PbctopoCutCandidate &partial,
                                     const PbctopoCutCandidate *threshold,
                                     PbctopoRefineScorePhase phase);
bool scaled_span_fits_unit_box(box Box, const vector<Vec3> &placed_points,
                               double tolerance = 1.0e-5);
bool continuity_metrics_ok(
    const vector<Vec3> &placed_points,
    const vector<PbctopoContinuityPair> &continuity_pairs, double hard_cutoff,
    double &path_sum, double &max_d2,
    const PbctopoCutCandidate *incumbent = nullptr);
void analyze_cut_candidate_ultraclash(const vector<Vec3> &placed_points,
                                       const unordered_set<uint64_t> &edge_set,
                                       PbctopoCutCandidate &candidate,
                                       double cutoff,
                                       const vector<PbctopoAtom> *atoms,
                                       const vector<size_t> *atom_component_index);
PbctopoAnchorBackend resolve_pbctopo_anchor_backend(const arguments &Arguments);
const char *pbctopo_anchor_backend_name(PbctopoAnchorBackend backend);
PbctopoSolveMode resolve_pbctopo_solve_mode(const arguments &Arguments);
const char *pbctopo_solve_mode_name(PbctopoSolveMode mode);
PbctopoCutCandidate evaluate_virtual_box_candidate(
    const vector<PbctopoAtom> &atoms, const vector<PbctopoLocalEdge> &local_edges,
    const unordered_set<uint64_t> &edge_set,
    const PbctopoTranslationCache &translation_cache, box Box, double cut_x,
    double cut_y, double cut_z, double gap_score, vector<Vec3> &placed_points,
    double clash_cutoff,
    const vector<PbctopoContinuityPair> *continuity_pairs = nullptr,
    const vector<size_t> *anchor_indices = nullptr,
    const vector<pair<int, int>> *contacting_owner_pairs = nullptr,
    double continuity_cutoff = 20.0,
    const PbctopoCutCandidate *incumbent = nullptr,
    bool skip_clash = false,
    double frame_time_ps = numeric_limits<double>::quiet_NaN(),
    const char *reject_context = nullptr,
    const PbctopoAnchorGraphContext *anchor_context = nullptr,
    PbctopoVirtualSeedSource seed_source = PbctopoVirtualSeedSource::None,
    const vector<size_t> *atom_component_index = nullptr);
void append_interface_anchor_edges(
    const vector<PbctopoAtom> &atoms, const vector<Vec3> &placed_points,
    const vector<size_t> &atom_component_index,
    const unordered_set<uint64_t> &edge_set, double contact_cutoff,
    vector<PbctopoLocalEdge> &edges);
PbctopoCutCandidate rescore_cut_candidate(
    const vector<PbctopoLocalEdge> &scoring_edges,
    const vector<Vec3> &placed_points, const PbctopoCutCandidate &source);

PbctopoRefineScoringCallbacks make_component_refine_callbacks() {
  PbctopoRefineScoringCallbacks callbacks;
  callbacks.score_cut_candidate = score_cut_candidate;
  callbacks.score_cut_candidate_principal_only =
      score_cut_candidate_principal_only;
  callbacks.continuity_metrics_ok = continuity_metrics_ok;
  callbacks.analyze_cut_candidate_ultraclash =
      analyze_cut_candidate_ultraclash;
  callbacks.better_cut_candidate = better_cut_candidate;
  callbacks.better_preclash_cut_candidate = better_preclash_cut_candidate;
  callbacks.better_stage1_cut_candidate = better_stage1_cut_candidate;
  callbacks.cut_candidate_cannot_beat_incumbent =
      cut_candidate_cannot_beat_incumbent;
  callbacks.edge_scan_cannot_beat_threshold = edge_scan_cannot_beat_threshold;
  callbacks.candidate_close_enough_to_skip_search =
      candidate_close_enough_to_skip_search;
  return callbacks;
}

[[maybe_unused]] Vec3 point_from_atom_copy(const PbctopoAtomCopy &atom) {
  return {atom.x, atom.y, atom.z};
}

PbctopoRuntimeCutoffs resolve_pbctopo_cutoffs(const arguments &Arguments) {
  return titan_pbctopo::resolve_runtime_cutoffs(Arguments);
}

PbctopoAnchorBackend resolve_pbctopo_anchor_backend(const arguments &Arguments) {
  return titan_pbctopo::resolve_anchor_backend(Arguments);
}

const char *pbctopo_anchor_backend_name(PbctopoAnchorBackend backend) {
  return titan_pbctopo::anchor_backend_name(backend);
}

PbctopoSolveMode resolve_pbctopo_solve_mode(const arguments &Arguments) {
  return titan_pbctopo::resolve_solve_mode(Arguments);
}

const char *pbctopo_solve_mode_name(PbctopoSolveMode mode) {
  return titan_pbctopo::solve_mode_name(mode);
}

uint64_t atom_pair_key(size_t a, size_t b) {
  if (a > b)
    swap(a, b);
  return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

void scaled_to_cartesian(box Box, double sx, double sy, double sz,
                         double &x, double &y, double &z) {
  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);
  x = sx * v1[0] + sy * v2[0] + sz * v3[0];
  y = sx * v1[1] + sy * v2[1] + sz * v3[1];
  z = sx * v1[2] + sy * v2[2] + sz * v3[2];
}

PbctopoTranslationCache build_virtual_box_translation_cache(box Box) {
  PbctopoTranslationCache cache;
  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);
  for (int ix = -2; ix <= 2; ++ix) {
    for (int iy = -2; iy <= 2; ++iy) {
      for (int iz = -2; iz <= 2; ++iz) {
        cache.by_shift[ix + 2][iy + 2][iz + 2] = {
            ix * v1[0] + iy * v2[0] + iz * v3[0],
            ix * v1[1] + iy * v2[1] + iz * v3[1],
            ix * v1[2] + iy * v2[2] + iz * v3[2]};
      }
    }
  }
  return cache;
}

uint64_t oriented_edge_key(int u, int v) {
  return (static_cast<uint64_t>(static_cast<uint32_t>(u)) << 32) |
         static_cast<uint32_t>(v);
}

bool partial_component_ultraclash(
    const vector<size_t> &component, const vector<Vec3> &placed_points,
    const PbctopoShiftedComponent &shifted_component,
    const PersistentClashGrid &clash_grid,
    const unordered_set<uint64_t> &edge_set, double cutoff,
    double &min_pair_d2, size_t &clash_atom_a, size_t &clash_atom_b,
    PbctopoTimingStats *timing = nullptr) {
  using clock = chrono::steady_clock;
  const auto started = clock::now();
  auto finish = [&](bool result) {
    if (timing != nullptr) {
      ++timing->component_partial_checks;
      timing->component_partial_s +=
          chrono::duration<double>(clock::now() - started).count();
      if (result)
        ++timing->component_partial_prunes;
    }
    return result;
  };
  if (component.empty())
    return finish(false);

  const double cutoff2 = cutoff * cutoff;

  for (size_t local_idx = 0; local_idx < component.size(); ++local_idx) {
    const size_t atom_idx = component[local_idx];
    const auto &point = shifted_component.points[local_idx];
    const auto &key = shifted_component.clash_keys[local_idx];
    const int cx = key.x;
    const int cy = key.y;
    const int cz = key.z;
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          const auto it = clash_grid.cells.find({cx + dx, cy + dy, cz + dz});
          if (it == clash_grid.cells.end())
            continue;
          for (size_t other_idx : it->second) {
            if (edge_set.find(atom_pair_key(atom_idx, other_idx)) != edge_set.end())
              continue;
            const auto &other = placed_points[other_idx];
            const double ddx = point.x - other.x;
            const double ddy = point.y - other.y;
            const double ddz = point.z - other.z;
            const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (d2 < min_pair_d2) {
              min_pair_d2 = d2;
              clash_atom_a = atom_idx;
              clash_atom_b = other_idx;
            }
            if (d2 + kDistanceEps < cutoff2)
              return finish(true);
          }
        }
      }
    }
  }

  return finish(false);
}

void insert_component_into_clash_grid(const vector<size_t> &component,
                                      const PbctopoShiftedComponent &shifted_component,
                                      PersistentClashGrid &clash_grid) {
  for (size_t local_idx = 0; local_idx < component.size(); ++local_idx) {
    const size_t atom_idx = component[local_idx];
    const ClashCellKey key = shifted_component.clash_keys[local_idx];
    auto &bucket = clash_grid.cells[key];
    clash_grid.atom_keys[atom_idx] = key;
    clash_grid.atom_offsets[atom_idx] = bucket.size();
    bucket.push_back(atom_idx);
  }
}

bool scaled_span_fits_unit_box(box Box, const vector<Vec3> &placed_points,
                               double tolerance) {
  if (placed_points.empty())
    return false;

  const PbctopoLatticeMetric lattice(Box);
  double min_sx = numeric_limits<double>::infinity();
  double min_sy = numeric_limits<double>::infinity();
  double min_sz = numeric_limits<double>::infinity();
  double max_sx = -numeric_limits<double>::infinity();
  double max_sy = -numeric_limits<double>::infinity();
  double max_sz = -numeric_limits<double>::infinity();

  for (const auto &point : placed_points) {
    const auto scaled = lattice.cartesian_to_scaled(point.x, point.y, point.z);
    const double sx = scaled[0];
    const double sy = scaled[1];
    const double sz = scaled[2];
    min_sx = min(min_sx, sx);
    min_sy = min(min_sy, sy);
    min_sz = min(min_sz, sz);
    max_sx = max(max_sx, sx);
    max_sy = max(max_sy, sy);
    max_sz = max(max_sz, sz);
  }

  return (max_sx - min_sx) <= (1.0 + tolerance) &&
         (max_sy - min_sy) <= (1.0 + tolerance) &&
         (max_sz - min_sz) <= (1.0 + tolerance);
}

vector<pair<int, int>>
build_contacting_owner_pairs(const vector<PbctopoLocalEdge> &local_edges,
                             const vector<PbctopoAtom> &atoms, box Box,
                             double contact_cutoff,
                             bool hard_required_only) {
  vector<pair<int, int>> owner_pairs;
  unordered_set<uint64_t> seen;
  seen.reserve(local_edges.size());
  auto register_owner_pair = [&](int owner_a, int owner_b) {
    if (owner_a == owner_b)
      return false;
    const int lo = min(owner_a, owner_b);
    const int hi = max(owner_a, owner_b);
    const uint64_t key =
        (static_cast<uint64_t>(static_cast<uint32_t>(lo)) << 32) |
        static_cast<uint32_t>(hi);
    if (!seen.insert(key).second)
      return false;
    owner_pairs.emplace_back(lo, hi);
    return true;
  };
  for (const auto &edge : local_edges) {
    if (edge.a >= atoms.size() || edge.b >= atoms.size())
      continue;
    if (hard_required_only &&
        edge.kind != PbctopoLocalEdgeKind::TopologyBond) {
      continue;
    }
    register_owner_pair(atoms[edge.a].owner, atoms[edge.b].owner);
  }

  if (hard_required_only)
    return owner_pairs;

  if (atoms.size() < 2 || contact_cutoff <= kDistanceEps)
    return owner_pairs;

  const BoxMatrix box_matrix(Box);
  const PbctopoLatticeMetric lattice(Box);
  const double cutoff2 = contact_cutoff * contact_cutoff;
  const double sx_bound =
      min(0.5, contact_cutoff * sqrt(box_matrix.inv[0][0] * box_matrix.inv[0][0] +
                                     box_matrix.inv[1][0] * box_matrix.inv[1][0] +
                                     box_matrix.inv[2][0] * box_matrix.inv[2][0]) +
                    1.0e-9);
  const double sy_bound =
      min(0.5, contact_cutoff * sqrt(box_matrix.inv[0][1] * box_matrix.inv[0][1] +
                                     box_matrix.inv[1][1] * box_matrix.inv[1][1] +
                                     box_matrix.inv[2][1] * box_matrix.inv[2][1]) +
                    1.0e-9);
  const double sz_bound =
      min(0.5, contact_cutoff * sqrt(box_matrix.inv[0][2] * box_matrix.inv[0][2] +
                                     box_matrix.inv[1][2] * box_matrix.inv[1][2] +
                                     box_matrix.inv[2][2] * box_matrix.inv[2][2]) +
                    1.0e-9);
  const int nx = max(1, static_cast<int>(floor(1.0 / max(sx_bound, 1.0e-6))));
  const int ny = max(1, static_cast<int>(floor(1.0 / max(sy_bound, 1.0e-6))));
  const int nz = max(1, static_cast<int>(floor(1.0 / max(sz_bound, 1.0e-6))));
  const double cell_sx = 1.0 / static_cast<double>(nx);
  const double cell_sy = 1.0 / static_cast<double>(ny);
  const double cell_sz = 1.0 / static_cast<double>(nz);
  auto wrap_index = [](int value, int mod) {
    int wrapped = value % mod;
    if (wrapped < 0)
      wrapped += mod;
    return wrapped;
  };

  unordered_map<ClashCellKey, vector<size_t>, ClashCellKeyHash> contact_grid;
  contact_grid.reserve(max<size_t>(64, atoms.size() / 2));
  for (size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
    const auto &atom = atoms[atom_idx];
    if (!atom.is_heavy)
      continue;
    const int cx = min(nx - 1, static_cast<int>(floor(atom.sx / cell_sx)));
    const int cy = min(ny - 1, static_cast<int>(floor(atom.sy / cell_sy)));
    const int cz = min(nz - 1, static_cast<int>(floor(atom.sz / cell_sz)));
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          const ClashCellKey key{wrap_index(cx + dx, nx), wrap_index(cy + dy, ny),
                                 wrap_index(cz + dz, nz)};
          const auto it = contact_grid.find(key);
          if (it == contact_grid.end())
            continue;
          for (size_t other_idx : it->second) {
            const auto &other = atoms[other_idx];
            if (atom.owner == other.owner)
              continue;
            const auto mic = lattice.nearest_image_cartesian(
                other.x - atom.x, other.y - atom.y, other.z - atom.z);
            const double d2 = mic.d2;
            if (d2 > cutoff2)
              continue;
            register_owner_pair(atom.owner, other.owner);
          }
        }
      }
    }
    contact_grid[{cx, cy, cz}].push_back(atom_idx);
  }
  return owner_pairs;
}

unordered_map<int, Vec3>
compute_owner_com_cartesian(const vector<Vec3> &placed_points,
                            const vector<PbctopoAtom> &atoms) {
  unordered_map<int, Vec3> owner_sums;
  unordered_map<int, size_t> owner_counts;
  owner_sums.reserve(8);
  owner_counts.reserve(8);
  const size_t atom_count = min(placed_points.size(), atoms.size());
  for (size_t atom_idx = 0; atom_idx < atom_count; ++atom_idx) {
    const int owner = atoms[atom_idx].owner;
    Vec3 &sum = owner_sums[owner];
    sum.x += placed_points[atom_idx].x;
    sum.y += placed_points[atom_idx].y;
    sum.z += placed_points[atom_idx].z;
    ++owner_counts[owner];
  }
  for (auto &[owner, sum] : owner_sums) {
    const size_t count = owner_counts[owner];
    if (count == 0)
      continue;
    const double inv_count = 1.0 / static_cast<double>(count);
    sum.x *= inv_count;
    sum.y *= inv_count;
    sum.z *= inv_count;
  }
  return owner_sums;
}

unordered_map<int, Vec3>
compute_owner_com_scaled(box Box,
                         const unordered_map<int, Vec3> &owner_com_cartesian) {
  unordered_map<int, Vec3> owner_com_scaled;
  owner_com_scaled.reserve(owner_com_cartesian.size());
  const PbctopoLatticeMetric lattice(Box);
  for (const auto &[owner, point] : owner_com_cartesian) {
    const auto scaled = lattice.cartesian_to_scaled(point.x, point.y, point.z);
    owner_com_scaled[owner] = {scaled[0], scaled[1], scaled[2]};
  }
  return owner_com_scaled;
}

bool candidate_has_split_contacting_owners(
    const vector<Vec3> &placed_points, const vector<PbctopoAtom> &atoms, box Box,
    const vector<pair<int, int>> &contacting_owner_pairs, size_t &split_owner_pairs,
    double &max_owner_raw_minus_minimg_delta,
    PbctopoOwnerPairSplitDetail *detail) {
  split_owner_pairs = 0;
  max_owner_raw_minus_minimg_delta = 0.0;
  if (detail != nullptr)
    *detail = {};
  if (contacting_owner_pairs.empty() || placed_points.empty() || atoms.empty())
    return false;

  const auto owner_com_cartesian = compute_owner_com_cartesian(placed_points, atoms);
  const auto owner_com_scaled = compute_owner_com_scaled(Box, owner_com_cartesian);
  const PbctopoLatticeMetric lattice(Box);
  const double box_scale = max(lattice.shortest_basis_length(), 1.0e-12);

  for (const auto &[owner_a, owner_b] : contacting_owner_pairs) {
    const auto it_a = owner_com_scaled.find(owner_a);
    const auto it_b = owner_com_scaled.find(owner_b);
    if (it_a == owner_com_scaled.end() || it_b == owner_com_scaled.end())
      continue;

    const Vec3 raw_delta_scaled{it_b->second.x - it_a->second.x,
                                it_b->second.y - it_a->second.y,
                                it_b->second.z - it_a->second.z};
    const auto mic = lattice.nearest_image_scaled(
        raw_delta_scaled.x, raw_delta_scaled.y, raw_delta_scaled.z);
    const Vec3 minimg_delta_scaled{mic.sx, mic.sy, mic.sz};
    const Vec3 split_delta{raw_delta_scaled.x - minimg_delta_scaled.x,
                           raw_delta_scaled.y - minimg_delta_scaled.y,
                           raw_delta_scaled.z - minimg_delta_scaled.z};
    const double axis_shift_max =
        max({fabs(split_delta.x), fabs(split_delta.y), fabs(split_delta.z)});
    const auto raw_cart = lattice.scaled_to_cartesian(
        raw_delta_scaled.x, raw_delta_scaled.y, raw_delta_scaled.z);
    const double raw_norm_A = sqrt(raw_cart[0] * raw_cart[0] +
                                   raw_cart[1] * raw_cart[1] +
                                   raw_cart[2] * raw_cart[2]);
    const double minimg_norm_A = sqrt(max(0.0, mic.d2));
    const double norm_gap = (raw_norm_A - minimg_norm_A) / box_scale;
    if (axis_shift_max <= 0.5 || norm_gap <= 0.20)
      continue;

    ++split_owner_pairs;
    if (norm_gap > max_owner_raw_minus_minimg_delta) {
      max_owner_raw_minus_minimg_delta = norm_gap;
      if (detail != nullptr) {
        detail->owner_a = owner_a;
        detail->owner_b = owner_b;
        detail->raw_delta_scaled = raw_delta_scaled;
        detail->minimg_delta_scaled = minimg_delta_scaled;
        detail->raw_norm = raw_norm_A;
        detail->minimg_norm = minimg_norm_A;
        detail->valid = true;
      }
    }
  }

  return split_owner_pairs > 0;
}

bool reimage_contacting_owner_points_minimum_image(
    vector<Vec3> &placed_points, const vector<PbctopoAtom> &atoms, box Box,
    const vector<pair<int, int>> &contacting_owner_pairs) {
  if (contacting_owner_pairs.empty() || placed_points.size() != atoms.size())
    return false;

  const auto owner_com_cartesian = compute_owner_com_cartesian(placed_points, atoms);
  const auto owner_com_scaled = compute_owner_com_scaled(Box, owner_com_cartesian);
  const PbctopoLatticeMetric lattice(Box);
  if (owner_com_scaled.size() < 2)
    return false;

  unordered_map<int, vector<int>> adjacency;
  adjacency.reserve(owner_com_scaled.size());
  for (const auto &[owner_a, owner_b] : contacting_owner_pairs) {
    if (owner_com_scaled.find(owner_a) == owner_com_scaled.end() ||
        owner_com_scaled.find(owner_b) == owner_com_scaled.end()) {
      continue;
    }
    adjacency[owner_a].push_back(owner_b);
    adjacency[owner_b].push_back(owner_a);
  }
  if (adjacency.empty())
    return false;

  vector<int> owners;
  owners.reserve(owner_com_scaled.size());
  for (const auto &[owner, _] : owner_com_scaled)
    owners.push_back(owner);
  sort(owners.begin(), owners.end());

  unordered_map<int, Int3> owner_shifts;
  owner_shifts.reserve(owner_com_scaled.size());
  queue<int> pending;
  for (int root_owner : owners) {
    if (owner_shifts.find(root_owner) != owner_shifts.end())
      continue;
    owner_shifts[root_owner] = {0, 0, 0};
    pending.push(root_owner);
    while (!pending.empty()) {
      const int owner_u = pending.front();
      pending.pop();
      const auto adj_it = adjacency.find(owner_u);
      if (adj_it == adjacency.end())
        continue;
      const Vec3 &com_u = owner_com_scaled.at(owner_u);
      const Int3 shift_u = owner_shifts.at(owner_u);
      for (int owner_v : adj_it->second) {
        const Vec3 &com_v = owner_com_scaled.at(owner_v);
        const auto mic = lattice.nearest_image_scaled(
            com_v.x - com_u.x, com_v.y - com_u.y, com_v.z - com_u.z);
        const Int3 proposed_shift{shift_u.x + mic.image_x,
                                  shift_u.y + mic.image_y,
                                  shift_u.z + mic.image_z};
        const auto [shift_it, inserted] =
            owner_shifts.emplace(owner_v, proposed_shift);
        if (inserted)
          pending.push(owner_v);
      }
    }
  }

  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);
  bool changed = false;
  for (size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
    const auto shift_it = owner_shifts.find(atoms[atom_idx].owner);
    if (shift_it == owner_shifts.end())
      continue;
    const Int3 &shift = shift_it->second;
    if (shift.x == 0 && shift.y == 0 && shift.z == 0)
      continue;
    placed_points[atom_idx].x += shift.x * v1[0] + shift.y * v2[0] + shift.z * v3[0];
    placed_points[atom_idx].y += shift.x * v1[1] + shift.y * v2[1] + shift.z * v3[1];
    placed_points[atom_idx].z += shift.x * v1[2] + shift.y * v2[2] + shift.z * v3[2];
    changed = true;
  }
  return changed;
}

void build_owner_spans(const vector<PbctopoAtom> &atoms, int owner,
                       vector<PbctopoHullSpan> &spans) {
  spans.clear();
  size_t span_begin = numeric_limits<size_t>::max();
  for (size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
    if (atoms[atom_idx].owner == owner) {
      if (span_begin == numeric_limits<size_t>::max())
        span_begin = atom_idx;
      continue;
    }
    if (span_begin != numeric_limits<size_t>::max()) {
      spans.push_back({span_begin, atom_idx});
      span_begin = numeric_limits<size_t>::max();
    }
  }
  if (span_begin != numeric_limits<size_t>::max())
    spans.push_back({span_begin, atoms.size()});
}

bool rescue_split_contacting_owner_candidate(
    const vector<Vec3> &seed_points, const vector<PbctopoAtom> &atoms, box Box,
    const vector<pair<int, int>> &contacting_owner_pairs,
    const vector<PbctopoShift> &shifts,
    const PbctopoOwnerPairSplitDetail &split_detail, vector<Vec3> &rescued_points,
    Int3 &applied_shift, double &best_raw_delta_norm) {
  rescued_points.clear();
  applied_shift = {};
  best_raw_delta_norm = numeric_limits<double>::infinity();
  if (!split_detail.valid || seed_points.size() != atoms.size() ||
      contacting_owner_pairs.empty() || shifts.empty()) {
    return false;
  }

  const auto owner_com_cartesian = compute_owner_com_cartesian(seed_points, atoms);
  const auto owner_com_scaled = compute_owner_com_scaled(Box, owner_com_cartesian);
  size_t best_shift_idx = numeric_limits<size_t>::max();
  double best_score = numeric_limits<double>::infinity();
  int owner_to_move = -1;
  vector<int> move_candidates;
  move_candidates.reserve(2);
  move_candidates.push_back(split_detail.owner_b);
  if (split_detail.owner_a != split_detail.owner_b)
    move_candidates.push_back(split_detail.owner_a);

  for (int move_owner : move_candidates) {
    const int ref_owner =
        move_owner == split_detail.owner_b ? split_detail.owner_a
                                           : split_detail.owner_b;
    const auto ref_it = owner_com_scaled.find(ref_owner);
    const auto move_it = owner_com_scaled.find(move_owner);
    if (ref_it == owner_com_scaled.end() || move_it == owner_com_scaled.end())
      continue;

    const Vec3 &ref_com = ref_it->second;
    const Vec3 &move_com = move_it->second;
    for (size_t shift_idx = 0; shift_idx < shifts.size(); ++shift_idx) {
      const auto &shift = shifts[shift_idx];
      const Vec3 shifted_move_scaled{move_com.x + static_cast<double>(shift.ix),
                                     move_com.y + static_cast<double>(shift.iy),
                                     move_com.z + static_cast<double>(shift.iz)};
      const Vec3 raw_delta_scaled{shifted_move_scaled.x - ref_com.x,
                                  shifted_move_scaled.y - ref_com.y,
                                  shifted_move_scaled.z - ref_com.z};
      const double raw_norm = sqrt(norm2(raw_delta_scaled));
      double score = raw_norm;
      for (const auto &[owner_a, owner_b] : contacting_owner_pairs) {
        if (owner_a != move_owner && owner_b != move_owner)
          continue;
        const int other_owner = (owner_a == move_owner) ? owner_b : owner_a;
        const auto other_it = owner_com_scaled.find(other_owner);
        if (other_it == owner_com_scaled.end())
          continue;
        Vec3 pair_raw_delta_scaled{};
        if (owner_a == move_owner) {
          pair_raw_delta_scaled = {
              other_it->second.x - shifted_move_scaled.x,
              other_it->second.y - shifted_move_scaled.y,
              other_it->second.z - shifted_move_scaled.z};
        } else {
          pair_raw_delta_scaled = {
              shifted_move_scaled.x - other_it->second.x,
              shifted_move_scaled.y - other_it->second.y,
              shifted_move_scaled.z - other_it->second.z};
        }
        score += 0.25 * sqrt(norm2(pair_raw_delta_scaled));
      }
      if (score + 1.0e-9 < best_score) {
        best_score = score;
        best_raw_delta_norm = raw_norm;
        best_shift_idx = shift_idx;
        applied_shift = {shift.ix, shift.iy, shift.iz};
        owner_to_move = move_owner;
      }
    }
  }

  if (best_shift_idx == numeric_limits<size_t>::max() || owner_to_move < 0)
    return false;
  if (applied_shift.x == 0 && applied_shift.y == 0 && applied_shift.z == 0)
    return false;

  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);
  const Vec3 translation{
      applied_shift.x * v1[0] + applied_shift.y * v2[0] + applied_shift.z * v3[0],
      applied_shift.x * v1[1] + applied_shift.y * v2[1] + applied_shift.z * v3[1],
      applied_shift.x * v1[2] + applied_shift.y * v2[2] + applied_shift.z * v3[2]};

  rescued_points = seed_points;
  for (size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
    if (atoms[atom_idx].owner != owner_to_move)
      continue;
    rescued_points[atom_idx].x += translation.x;
    rescued_points[atom_idx].y += translation.y;
    rescued_points[atom_idx].z += translation.z;
  }
  return true;
}

bool continuity_metrics_ok(
    const vector<Vec3> &placed_points,
    const vector<PbctopoContinuityPair> &continuity_pairs, double hard_cutoff,
    double &path_sum, double &max_d2, const PbctopoCutCandidate *incumbent) {
  path_sum = 0.0;
  max_d2 = 0.0;
  if (continuity_pairs.empty())
    return true;

  const double hard_cutoff2 = hard_cutoff * hard_cutoff;
  const bool has_incumbent =
      (incumbent != nullptr && incumbent->valid &&
       isfinite(incumbent->continuity_path) &&
       isfinite(incumbent->continuity_max_d2));

  auto consume_pair_metric = [&](double d2, double distance) -> bool {
    if (d2 > max_d2)
      max_d2 = d2;
    if (d2 > hard_cutoff2)
      return false;
    path_sum += distance;
    if (has_incumbent && max_d2 > incumbent->continuity_max_d2 + 1.0e-6 &&
        path_sum > incumbent->continuity_path + 1.0e-6) {
      return false;
    }
    return true;
  };

#if TITAN_SIMD_X86 && TITAN_RUNTIME_TARGET_ATTRS
  if (titan_simd_dispatch::runtime_has_avx2() &&
      continuity_pairs.size() >= 8) {
    alignas(32) double ax[4];
    alignas(32) double ay[4];
    alignas(32) double az[4];
    alignas(32) double bx[4];
    alignas(32) double by[4];
    alignas(32) double bz[4];
    alignas(32) double d2_lane[4];
    alignas(32) double dist_lane[4];
    const __m256d zero = _mm256_setzero_pd();
    size_t pair_idx = 0;
    for (; pair_idx + 4 <= continuity_pairs.size(); pair_idx += 4) {
      for (size_t lane = 0; lane < 4; ++lane) {
        const auto &pair = continuity_pairs[pair_idx + lane];
        const Vec3 &point_a = placed_points[pair.a];
        const Vec3 &point_b = placed_points[pair.b];
        ax[lane] = point_a.x;
        ay[lane] = point_a.y;
        az[lane] = point_a.z;
        bx[lane] = point_b.x;
        by[lane] = point_b.y;
        bz[lane] = point_b.z;
      }

      const __m256d dx =
          _mm256_sub_pd(_mm256_load_pd(ax), _mm256_load_pd(bx));
      const __m256d dy =
          _mm256_sub_pd(_mm256_load_pd(ay), _mm256_load_pd(by));
      const __m256d dz =
          _mm256_sub_pd(_mm256_load_pd(az), _mm256_load_pd(bz));
      const __m256d d2 = _mm256_add_pd(
          _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)),
          _mm256_mul_pd(dz, dz));
      const __m256d dist = _mm256_sqrt_pd(_mm256_max_pd(d2, zero));
      _mm256_store_pd(d2_lane, d2);
      _mm256_store_pd(dist_lane, dist);

      for (size_t lane = 0; lane < 4; ++lane) {
        if (!consume_pair_metric(d2_lane[lane], dist_lane[lane]))
          return false;
      }
    }
    for (; pair_idx < continuity_pairs.size(); ++pair_idx) {
      const auto &pair = continuity_pairs[pair_idx];
      const Vec3 delta = placed_points[pair.a] - placed_points[pair.b];
      const double d2 = norm2(delta);
      if (!consume_pair_metric(d2, sqrt(max(d2, 0.0))))
        return false;
    }
    return true;
  }
#endif

  for (const auto &pair : continuity_pairs) {
    const Vec3 delta = placed_points[pair.a] - placed_points[pair.b];
    const double d2 = norm2(delta);
    if (!consume_pair_metric(d2, sqrt(max(d2, 0.0))))
      return false;
  }
  return true;
}

vector<size_t>
component_size_signature(const vector<vector<size_t>> &components) {
  vector<size_t> sizes;
  sizes.reserve(components.size());
  for (const auto &component : components)
    sizes.push_back(component.size());
  return sizes;
}

vector<Int3> extract_component_offsets_from_points(
    const vector<PbctopoAtom> &atoms, const vector<vector<size_t>> &components,
    box Box, const vector<Vec3> &placed_points) {
  vector<Int3> offsets;
  if (placed_points.size() != atoms.size())
    return offsets;
  offsets.resize(components.size());
  const PbctopoLatticeMetric lattice(Box);
  for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
    if (components[comp_idx].empty())
      continue;
    const size_t atom_idx = components[comp_idx].front();
    const double dx = placed_points[atom_idx].x - atoms[atom_idx].x;
    const double dy = placed_points[atom_idx].y - atoms[atom_idx].y;
    const double dz = placed_points[atom_idx].z - atoms[atom_idx].z;
    const auto scaled = lattice.cartesian_to_scaled(dx, dy, dz);
    offsets[comp_idx] = {static_cast<int>(llround(scaled[0])),
                         static_cast<int>(llround(scaled[1])),
                         static_cast<int>(llround(scaled[2]))};
  }
  rebase_component_offsets(offsets, Box);
  return offsets;
}

const char *pbctopo_component_mode_name(PbctopoComponentMode mode) {
  switch (mode) {
  case PbctopoComponentMode::Disabled:
    return "disabled";
  case PbctopoComponentMode::Exact:
    return "small_graph";
  case PbctopoComponentMode::MergeThenGraph:
    return "merge_then_graph";
  case PbctopoComponentMode::Beam:
    return "beam";
  }
  return "unknown";
}

PbctopoComponentMode select_pbctopo_component_mode(size_t raw_component_count,
                                                   size_t exact_component_limit) {
  if (raw_component_count == 0)
    return PbctopoComponentMode::Disabled;
  const size_t bounded_exact_limit =
      min(exact_component_limit, kPbctopoExactComponentLimitMax);
  if (raw_component_count <= bounded_exact_limit)
    return PbctopoComponentMode::Exact;
  if (raw_component_count <= kPbctopoMergeModeMaxRawComponents)
    return PbctopoComponentMode::MergeThenGraph;
  return PbctopoComponentMode::Beam;
}

struct PbctopoRawComponentMetadata {
  size_t raw_index = 0;
  int owner = 0;
  uint64_t chain_key = 0;
  size_t atom_count = 0;
  size_t min_local_atom_index = numeric_limits<size_t>::max();
  size_t max_local_atom_index = 0;
  size_t min_atom_index = numeric_limits<size_t>::max();
  size_t first_atom_index = numeric_limits<size_t>::max();
  size_t last_atom_index = numeric_limits<size_t>::max();
};

PbctopoMergedComponentPlan build_chain_merged_component_plan(
    const vector<PbctopoAtom> &atoms, box Box,
    const vector<vector<size_t>> &raw_components) {
  PbctopoMergedComponentPlan plan;
  plan.raw_to_merged.assign(raw_components.size(), numeric_limits<size_t>::max());
  plan.raw_internal_offsets.assign(raw_components.size(), {});
  if (raw_components.empty())
    return plan;

  vector<PbctopoRawComponentMetadata> metadata;
  metadata.reserve(raw_components.size());
  for (size_t raw_idx = 0; raw_idx < raw_components.size(); ++raw_idx) {
    const auto &component = raw_components[raw_idx];
    if (component.empty())
      continue;
    PbctopoRawComponentMetadata meta;
    meta.raw_index = raw_idx;
    meta.atom_count = component.size();
    const auto &first_atom = atoms[component.front()];
    meta.owner = first_atom.owner;
    meta.chain_key = first_atom.chain_key;
    for (size_t atom_idx : component) {
      if (atom_idx >= atoms.size())
        continue;
      meta.min_local_atom_index =
          min(meta.min_local_atom_index, atoms[atom_idx].local_atom_index);
      meta.max_local_atom_index =
          max(meta.max_local_atom_index, atoms[atom_idx].local_atom_index);
      meta.min_atom_index = min(meta.min_atom_index, atom_idx);
      if (meta.first_atom_index == numeric_limits<size_t>::max() ||
          atoms[atom_idx].local_atom_index <
              atoms[meta.first_atom_index].local_atom_index) {
        meta.first_atom_index = atom_idx;
      }
      if (meta.last_atom_index == numeric_limits<size_t>::max() ||
          atoms[atom_idx].local_atom_index >
              atoms[meta.last_atom_index].local_atom_index) {
        meta.last_atom_index = atom_idx;
      }
    }
    metadata.push_back(meta);
  }

  stable_sort(metadata.begin(), metadata.end(),
              [](const PbctopoRawComponentMetadata &lhs,
                 const PbctopoRawComponentMetadata &rhs) {
                if (lhs.owner != rhs.owner)
                  return lhs.owner < rhs.owner;
                if (lhs.chain_key != rhs.chain_key)
                  return lhs.chain_key < rhs.chain_key;
                if (lhs.min_local_atom_index != rhs.min_local_atom_index)
                  return lhs.min_local_atom_index < rhs.min_local_atom_index;
                return lhs.raw_index < rhs.raw_index;
              });

  vector<size_t> current_raw_group;
  vector<size_t> current_atoms;
  vector<size_t> first_atom_by_raw(raw_components.size(),
                                   numeric_limits<size_t>::max());
  vector<size_t> last_atom_by_raw(raw_components.size(),
                                  numeric_limits<size_t>::max());
  for (const auto &meta : metadata) {
    first_atom_by_raw[meta.raw_index] = meta.first_atom_index;
    last_atom_by_raw[meta.raw_index] = meta.last_atom_index;
  }
  const PbctopoLatticeMetric lattice(Box);
  current_raw_group.reserve(kPbctopoMergedNodeMaxRawComponents);
  current_atoms.reserve(kPbctopoMergedNodeMaxAtoms);
  int current_owner = 0;
  uint64_t current_chain_key = 0;
  size_t current_last_local_atom = 0;
  bool has_current_group = false;

  auto flush_group = [&]() {
    if (current_raw_group.empty())
      return;
    Int3 previous_offset{};
    for (size_t pos = 1; pos < current_raw_group.size(); ++pos) {
      const size_t previous_raw = current_raw_group[pos - 1];
      const size_t current_raw = current_raw_group[pos];
      const size_t previous_atom = last_atom_by_raw[previous_raw];
      const size_t current_atom = first_atom_by_raw[current_raw];
      if (previous_atom >= atoms.size() || current_atom >= atoms.size())
        continue;
      const auto mic = lattice.nearest_image_scaled(
          atoms[current_atom].sx -
              (atoms[previous_atom].sx + previous_offset.x),
          atoms[current_atom].sy -
              (atoms[previous_atom].sy + previous_offset.y),
          atoms[current_atom].sz -
              (atoms[previous_atom].sz + previous_offset.z));
      if (!mic.valid)
        continue;
      previous_offset = {mic.image_x, mic.image_y, mic.image_z};
      plan.raw_internal_offsets[current_raw] = previous_offset;
    }
    const size_t merged_idx = plan.merged_components.size();
    sort(current_atoms.begin(), current_atoms.end());
    plan.merged_components.push_back(current_atoms);
    plan.merged_to_raw.push_back(current_raw_group);
    for (size_t raw_idx : current_raw_group) {
      if (raw_idx < plan.raw_to_merged.size())
        plan.raw_to_merged[raw_idx] = merged_idx;
    }
    current_raw_group.clear();
    current_atoms.clear();
    has_current_group = false;
  };

  for (const auto &meta : metadata) {
    const bool contiguous_with_current =
        has_current_group && meta.owner == current_owner &&
        meta.chain_key == current_chain_key &&
        meta.min_local_atom_index == current_last_local_atom + 1;
    const bool fits_current =
        contiguous_with_current &&
        current_raw_group.size() < kPbctopoMergedNodeMaxRawComponents &&
        current_atoms.size() + meta.atom_count <= kPbctopoMergedNodeMaxAtoms;
    if (!fits_current)
      flush_group();
    if (!has_current_group) {
      current_owner = meta.owner;
      current_chain_key = meta.chain_key;
      current_last_local_atom = meta.max_local_atom_index;
      has_current_group = true;
    } else {
      current_last_local_atom = meta.max_local_atom_index;
    }
    current_raw_group.push_back(meta.raw_index);
    current_atoms.insert(current_atoms.end(), raw_components[meta.raw_index].begin(),
                         raw_components[meta.raw_index].end());
  }
  flush_group();

  return plan;
}

vector<Int3> expand_merged_offsets_to_raw_components(
    const vector<Int3> &merged_offsets, box Box,
    const vector<vector<size_t>> &raw_components,
    const PbctopoMergedComponentPlan &plan) {
  vector<Int3> raw_offsets(raw_components.size(), {0, 0, 0});
  for (size_t merged_idx = 0; merged_idx < plan.merged_to_raw.size(); ++merged_idx) {
    const Int3 merged_offset =
        merged_idx < merged_offsets.size() ? merged_offsets[merged_idx] : Int3{};
    for (size_t raw_idx : plan.merged_to_raw[merged_idx]) {
      if (raw_idx < raw_offsets.size()) {
        const Int3 internal = raw_idx < plan.raw_internal_offsets.size()
                                  ? plan.raw_internal_offsets[raw_idx]
                                  : Int3{};
        raw_offsets[raw_idx] = merged_offset + internal;
      }
    }
  }
  rebase_component_offsets(raw_offsets, Box);
  return raw_offsets;
}

bool derive_merged_warm_offsets(const vector<Int3> &raw_offsets,
                                const PbctopoMergedComponentPlan &plan,
                                vector<Int3> &merged_offsets) {
  if (raw_offsets.empty() || plan.merged_to_raw.empty())
    return false;
  merged_offsets.clear();
  merged_offsets.reserve(plan.merged_to_raw.size());
  for (const auto &group : plan.merged_to_raw) {
    if (group.empty())
      return false;
    const size_t first_raw_idx = group.front();
    if (first_raw_idx >= raw_offsets.size())
      return false;
    const Int3 first_internal =
        first_raw_idx < plan.raw_internal_offsets.size()
            ? plan.raw_internal_offsets[first_raw_idx]
            : Int3{};
    const Int3 candidate = raw_offsets[first_raw_idx] - first_internal;
    for (size_t raw_idx : group) {
      if (raw_idx >= raw_offsets.size())
        return false;
      const Int3 internal = raw_idx < plan.raw_internal_offsets.size()
                                ? plan.raw_internal_offsets[raw_idx]
                                : Int3{};
      const Int3 residual = raw_offsets[raw_idx] - internal;
      if (!same_delta(residual, candidate)) {
        return false;
      }
    }
    merged_offsets.push_back(candidate);
  }
  return true;
}

vector<PbctopoAtom> apply_merged_component_internal_offsets(
    const vector<PbctopoAtom> &atoms,
    const vector<vector<size_t>> &raw_components,
    const PbctopoMergedComponentPlan &plan, box Box) {
  vector<PbctopoAtom> aligned = atoms;
  const PbctopoLatticeMetric lattice(Box);
  const size_t raw_count =
      min(raw_components.size(), plan.raw_internal_offsets.size());
  for (size_t raw_idx = 0; raw_idx < raw_count; ++raw_idx) {
    const Int3 offset = plan.raw_internal_offsets[raw_idx];
    for (size_t atom_idx : raw_components[raw_idx]) {
      if (atom_idx >= aligned.size())
        continue;
      auto &atom = aligned[atom_idx];
      atom.sx += offset.x;
      atom.sy += offset.y;
      atom.sz += offset.z;
      const auto cart = lattice.scaled_to_cartesian(atom.sx, atom.sy, atom.sz);
      atom.x = cart[0];
      atom.y = cart[1];
      atom.z = cart[2];
    }
  }
  return aligned;
}

bool graph_seed_quick_score_better(const PbctopoGraphSeedQuickScore &lhs,
                                   const PbctopoGraphSeedQuickScore &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  int metric_cmp = compare_pbctopo_double_metric(
      lhs.continuity_max_d2, rhs.continuity_max_d2, 1.0e-6);
  if (metric_cmp != 0)
    return metric_cmp < 0;
  metric_cmp = compare_pbctopo_double_metric(lhs.continuity_path,
                                             rhs.continuity_path, 1.0e-6);
  if (metric_cmp != 0)
    return metric_cmp < 0;
  int anchor_cmp = compare_anchor_graph_score(lhs.anchor_graph,
                                              rhs.anchor_graph);
  if (anchor_cmp == 0) {
    anchor_cmp = compare_pbctopo_double_metric(lhs.anchor_max_d2,
                                               rhs.anchor_max_d2, 1.0e-6);
  }
  if (anchor_cmp == 0) {
    anchor_cmp = compare_pbctopo_double_metric(lhs.anchor_mst2,
                                               rhs.anchor_mst2, 1.0e-6);
  }
  if (anchor_cmp != 0)
    return anchor_cmp < 0;
  const int principal_max_cmp =
      compare_principal_rg2(lhs.principal_rg2_max, rhs.principal_rg2_max);
  if (principal_max_cmp != 0)
    return principal_max_cmp < 0;
  const int principal_mid_cmp =
      compare_principal_rg2(lhs.principal_rg2_mid, rhs.principal_rg2_mid);
  if (principal_mid_cmp != 0)
    return principal_mid_cmp < 0;
  const int principal_min_cmp =
      compare_principal_rg2(lhs.principal_rg2_min, rhs.principal_rg2_min);
  if (principal_min_cmp != 0)
    return principal_min_cmp < 0;
  metric_cmp = compare_pbctopo_double_metric(lhs.shift_norm2, rhs.shift_norm2,
                                              1.0e-6);
  if (metric_cmp != 0)
    return metric_cmp < 0;
  return false;
}

PbctopoGraphSeedQuickScore score_component_graph_seed_quick(
    const vector<PbctopoAtom> &atoms, const vector<vector<size_t>> &components,
    const vector<Vec3> &placed_points, const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices, double continuity_cutoff,
    double shift_norm2,
    PbctopoAnchorBackend anchor_backend = PbctopoAnchorBackend::Auto,
    const PbctopoAnchorGraphContext *anchor_context = nullptr) {
  PbctopoGraphSeedQuickScore score;
  if (placed_points.size() != atoms.size())
    return score;
  if (!continuity_metrics_ok(placed_points, continuity_pairs, continuity_cutoff,
                             score.continuity_path, score.continuity_max_d2,
                             nullptr)) {
    return score;
  }
  const vector<size_t> half_anchor_indices =
      build_coarse_anchor_subset(atoms, anchor_indices);
  const vector<size_t> third_anchor_indices =
      build_third_anchor_subset(atoms, anchor_indices);
  const vector<size_t> *quick_anchor_indices = &anchor_indices;
  if (anchor_indices.size() > 256 && third_anchor_indices.size() >= 2) {
    quick_anchor_indices = &third_anchor_indices;
  } else if (half_anchor_indices.size() >= 2) {
    quick_anchor_indices = &half_anchor_indices;
  }
  if (!anchor_graph_metrics_enhanced(
          placed_points, *quick_anchor_indices, continuity_pairs,
          score.anchor_mst2, score.anchor_max_d2, score.anchor_graph, nullptr,
          PbctopoRefineScorePhase::Exact, anchor_backend, nullptr,
          anchor_context)) {
    return score;
  }
  const PointStats stats = build_point_stats(placed_points);
  if (stats.count == 0)
    return score;
  const array<double, 3> principal_rg2 = compute_principal_rg2(stats);
  score.principal_rg2_max = principal_rg2[0];
  score.principal_rg2_mid = principal_rg2[1];
  score.principal_rg2_min = principal_rg2[2];
  score.shift_norm2 = shift_norm2;
  score.valid = true;
  return score;
}

PbctopoVirtualSeedCenter
virtual_seed_center_from_candidate(const PbctopoCutCandidate &candidate) {
  return {candidate.cut_x, candidate.cut_y, candidate.cut_z};
}

bool kinematic_virtual_candidate_usable(
    const PbctopoCutCandidate &candidate) {
  return candidate.valid && !candidate.ultra_clash &&
         isfinite(candidate.cut_x) && isfinite(candidate.cut_y) &&
         isfinite(candidate.cut_z);
}

bool has_consecutive_kinematic_seed_history(const PbctopoWarmState &warm_state,
                                            size_t frame_idx) {
  return warm_state.has_kinematic_previous_virtual_candidate &&
         warm_state.has_kinematic_latest_virtual_candidate &&
         warm_state.kinematic_previous_virtual_frame + 1 ==
             warm_state.kinematic_latest_virtual_frame &&
         warm_state.kinematic_latest_virtual_frame + 1 == frame_idx &&
         kinematic_virtual_candidate_usable(
             warm_state.kinematic_previous_virtual_candidate) &&
         kinematic_virtual_candidate_usable(
             warm_state.kinematic_latest_virtual_candidate);
}

void reset_kinematic_seed_history(PbctopoWarmState &warm_state,
                                  PbctopoTimingStats &timing) {
  if (warm_state.has_kinematic_previous_virtual_candidate ||
      warm_state.has_kinematic_latest_virtual_candidate) {
    ++timing.pbctopo_kinematic_seed_history_resets;
  }
  warm_state.kinematic_previous_virtual_candidate = PbctopoCutCandidate{};
  warm_state.kinematic_latest_virtual_candidate = PbctopoCutCandidate{};
  warm_state.kinematic_previous_virtual_frame =
      numeric_limits<size_t>::max();
  warm_state.kinematic_latest_virtual_frame =
      numeric_limits<size_t>::max();
  warm_state.has_kinematic_previous_virtual_candidate = false;
  warm_state.has_kinematic_latest_virtual_candidate = false;
}

void update_kinematic_seed_history(PbctopoWarmState &warm_state,
                                   const PbctopoCutCandidate &candidate,
                                   size_t frame_idx,
                                   PbctopoTimingStats &timing) {
  if (!kinematic_virtual_candidate_usable(candidate)) {
    reset_kinematic_seed_history(warm_state, timing);
    return;
  }

  if (warm_state.has_kinematic_latest_virtual_candidate &&
      warm_state.kinematic_latest_virtual_frame + 1 == frame_idx) {
    warm_state.kinematic_previous_virtual_candidate =
        warm_state.kinematic_latest_virtual_candidate;
    warm_state.kinematic_previous_virtual_frame =
        warm_state.kinematic_latest_virtual_frame;
    warm_state.has_kinematic_previous_virtual_candidate = true;
  } else if (warm_state.has_kinematic_latest_virtual_candidate ||
             warm_state.has_kinematic_previous_virtual_candidate) {
    reset_kinematic_seed_history(warm_state, timing);
  }

  warm_state.kinematic_latest_virtual_candidate = candidate;
  warm_state.kinematic_latest_virtual_frame = frame_idx;
  warm_state.has_kinematic_latest_virtual_candidate = true;
}

struct PbctopoGraphSeedSolveResult {
  const char *seed_name = "baseline";
  vector<Int3> offsets;
  PbctopoGraphSeedQuickScore quick_score;
  bool valid = false;
};

vector<PbctopoGraphSeedSolveResult> build_large_component_graph_seed_candidates(
    const vector<PbctopoAtom> &atoms, box Box,
    const vector<vector<size_t>> &components, const vector<PbctopoShift> &shifts,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices, const PbctopoRuntimeCutoffs &cutoffs,
    const vector<Int3> *warm_component_offsets, PbctopoTimingStats *timing,
    PbctopoAnchorBackend anchor_backend = PbctopoAnchorBackend::Auto,
    const PbctopoAnchorGraphContext *anchor_context = nullptr) {
  const vector<PbctopoGraphEdge> base_edges = build_component_graph_edge_candidates(
      atoms, components, shifts, cutoffs, timing, warm_component_offsets);
  vector<PbctopoGraphSeedSolveResult> results;
  results.reserve(3);
  unordered_set<string> seen_offset_signatures;

  auto try_seed = [&](const char *seed_name,
                      const function<void(vector<PbctopoGraphEdge> &)> &sort_fn) {
    vector<PbctopoGraphEdge> ranked_edges = base_edges;
    sort_fn(ranked_edges);
    vector<Int3> offsets;
    if (!solve_component_offsets_from_ranked_edges(components.size(), Box, components,
                                                   ranked_edges, offsets)) {
      return;
    }
    const string signature = format_component_offsets(offsets);
    if (!seen_offset_signatures.insert(signature).second)
      return;
    vector<Vec3> seed_points;
    double seed_shift_norm2 = 0.0;
    build_component_points_from_offsets(atoms, components, Box, offsets, seed_points,
                                        seed_shift_norm2);
    PbctopoGraphSeedSolveResult result;
    result.seed_name = seed_name;
    result.offsets = offsets;
    result.quick_score = score_component_graph_seed_quick(
        atoms, components, seed_points, continuity_pairs, anchor_indices,
        cutoffs.continuity_cutoff, seed_shift_norm2, anchor_backend,
        anchor_context);
    result.valid = result.quick_score.valid;
    results.push_back(std::move(result));
  };

  try_seed("baseline", [&](vector<PbctopoGraphEdge> &ranked_edges) {
    sort(ranked_edges.begin(), ranked_edges.end(), better_component_graph_edge);
  });

  try_seed("nonzero", [&](vector<PbctopoGraphEdge> &ranked_edges) {
    stable_sort(ranked_edges.begin(), ranked_edges.end(),
                [](const PbctopoGraphEdge &lhs, const PbctopoGraphEdge &rhs) {
                  if (!lhs.valid)
                    return false;
                  if (!rhs.valid)
                    return true;
                  if (lhs.ultra_clash != rhs.ultra_clash)
                    return rhs.ultra_clash;
                  const bool lhs_zero = graph_edge_has_zero_shift(lhs);
                  const bool rhs_zero = graph_edge_has_zero_shift(rhs);
                  if (lhs_zero != rhs_zero)
                    return !lhs_zero;
                  int metric_cmp = compare_pbctopo_double_metric(
                      lhs.contact_fraction, rhs.contact_fraction, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp > 0;
                  metric_cmp = compare_pbctopo_double_metric(
                      lhs.confidence, rhs.confidence, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp > 0;
                  metric_cmp = compare_pbctopo_double_metric(
                      lhs.soft_penalty_mean, rhs.soft_penalty_mean, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp < 0;
                  metric_cmp = compare_pbctopo_double_metric(
                      lhs.min_pair_d2, rhs.min_pair_d2, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp > 0;
                  metric_cmp = compare_pbctopo_double_metric(
                      lhs.bbox_d2, rhs.bbox_d2, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp < 0;
                  metric_cmp = compare_pbctopo_double_metric(
                      lhs.shift_norm2, rhs.shift_norm2, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp < 0;
                  if (lhs.a != rhs.a)
                    return lhs.a < rhs.a;
                  return lhs.b < rhs.b;
                });
  });

  try_seed("warm", [&](vector<PbctopoGraphEdge> &ranked_edges) {
    if (warm_component_offsets == nullptr ||
        warm_component_offsets->size() != components.size()) {
      stable_sort(ranked_edges.begin(), ranked_edges.end(),
                  [](const PbctopoGraphEdge &lhs, const PbctopoGraphEdge &rhs) {
                    if (!lhs.valid)
                      return false;
                    if (!rhs.valid)
                      return true;
                    if (lhs.ultra_clash != rhs.ultra_clash)
                      return rhs.ultra_clash;
                    const bool lhs_zero = graph_edge_has_zero_shift(lhs);
                    const bool rhs_zero = graph_edge_has_zero_shift(rhs);
                    if (lhs_zero != rhs_zero)
                      return !lhs_zero;
                    return better_component_graph_edge(lhs, rhs);
                  });
      return;
    }
    stable_sort(ranked_edges.begin(), ranked_edges.end(),
                [&](const PbctopoGraphEdge &lhs, const PbctopoGraphEdge &rhs) {
                  if (!lhs.valid)
                    return false;
                  if (!rhs.valid)
                    return true;
                  if (lhs.ultra_clash != rhs.ultra_clash)
                    return rhs.ultra_clash;
                  const Int3 lhs_expected =
                      (*warm_component_offsets)[lhs.b] - (*warm_component_offsets)[lhs.a];
                  const Int3 rhs_expected =
                      (*warm_component_offsets)[rhs.b] - (*warm_component_offsets)[rhs.a];
                  const bool lhs_warm = graph_edge_matches_delta(lhs, lhs_expected);
                  const bool rhs_warm = graph_edge_matches_delta(rhs, rhs_expected);
                  if (lhs_warm != rhs_warm)
                    return lhs_warm;
                  const bool lhs_zero = graph_edge_has_zero_shift(lhs);
                  const bool rhs_zero = graph_edge_has_zero_shift(rhs);
                  if (lhs_zero != rhs_zero)
                    return !lhs_zero;
                  int metric_cmp = compare_pbctopo_double_metric(
                      lhs.contact_fraction, rhs.contact_fraction, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp > 0;
                  metric_cmp = compare_pbctopo_double_metric(
                      lhs.confidence, rhs.confidence, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp > 0;
                  metric_cmp = compare_pbctopo_double_metric(
                      lhs.soft_penalty_mean, rhs.soft_penalty_mean, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp < 0;
                  metric_cmp = compare_pbctopo_double_metric(
                      lhs.min_pair_d2, rhs.min_pair_d2, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp > 0;
                  metric_cmp = compare_pbctopo_double_metric(
                      lhs.bbox_d2, rhs.bbox_d2, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp < 0;
                  metric_cmp = compare_pbctopo_double_metric(
                      lhs.shift_norm2, rhs.shift_norm2, 1.0e-9);
                  if (metric_cmp != 0)
                    return metric_cmp < 0;
                  if (lhs.a != rhs.a)
                    return lhs.a < rhs.a;
                  return lhs.b < rhs.b;
                });
  });

  return results;
}

void remove_component_from_clash_grid(const vector<size_t> &component,
                                      PersistentClashGrid &clash_grid) {
  for (size_t atom_idx : component) {
    const ClashCellKey key = clash_grid.atom_keys[atom_idx];
    const auto it = clash_grid.cells.find(key);
    if (it == clash_grid.cells.end())
      continue;
    auto &bucket = it->second;
    const size_t remove_pos = clash_grid.atom_offsets[atom_idx];
    if (remove_pos != bucket.size() - 1) {
      const size_t moved_atom = bucket.back();
      bucket[remove_pos] = moved_atom;
      clash_grid.atom_offsets[moved_atom] = remove_pos;
    }
    bucket.pop_back();
    if (bucket.empty())
      clash_grid.cells.erase(it);
  }
}

vector<vector<size_t>>
build_connected_components_from_local_edges(size_t atom_count,
                                            const vector<PbctopoLocalEdge> &edges) {
  vector<vector<size_t>> adjacency(atom_count);
  adjacency.reserve(atom_count);
  for (const auto &edge : edges) {
    adjacency[edge.a].push_back(edge.b);
    adjacency[edge.b].push_back(edge.a);
  }

  vector<vector<size_t>> components;
  vector<char> seen(atom_count, 0);
  vector<size_t> stack;
  for (size_t atom = 0; atom < atom_count; ++atom) {
    if (seen[atom] != 0)
      continue;
    seen[atom] = 1;
    stack.clear();
    stack.push_back(atom);
    vector<size_t> component;
    while (!stack.empty()) {
      const size_t current = stack.back();
      stack.pop_back();
      component.push_back(current);
      for (size_t nb : adjacency[current]) {
        if (seen[nb] != 0)
          continue;
        seen[nb] = 1;
        stack.push_back(nb);
      }
    }
    sort(component.begin(), component.end());
    components.push_back(std::move(component));
  }

  sort(components.begin(), components.end(),
       [](const vector<size_t> &lhs, const vector<size_t> &rhs) {
         if (lhs.size() != rhs.size())
           return lhs.size() > rhs.size();
         return lhs.front() < rhs.front();
       });
  return components;
}

struct PbctopoEdgeScanScratch {
  AlignedDoubleVector xs;
  AlignedDoubleVector ys;
  AlignedDoubleVector zs;

  void resize(size_t point_count) {
    xs.resize(point_count);
    ys.resize(point_count);
    zs.resize(point_count);
  }
};

PbctopoEdgeScanScratch &pbctopo_edge_scan_scratch() {
  thread_local PbctopoEdgeScanScratch scratch;
  return scratch;
}

void sync_edge_scan_coords_soa(const vector<Vec3> &placed_points,
                               PbctopoEdgeScanScratch &scratch) {
  scratch.resize(placed_points.size());
  for (size_t i = 0; i < placed_points.size(); ++i) {
    scratch.xs[i] = placed_points[i].x;
    scratch.ys[i] = placed_points[i].y;
    scratch.zs[i] = placed_points[i].z;
  }
}

bool update_cut_candidate_with_edge_distance(
    PbctopoCutCandidate &candidate, double d2, double ref_d2, double inv_ref_d2,
    double break_cutoff_d2, const PbctopoCutCandidate *edge_prune_threshold,
    PbctopoRefineScorePhase phase) {
  if (d2 > break_cutoff_d2)
    ++candidate.broken_edges;
  if (d2 > ref_d2)
    candidate.edge_stretch += (d2 - ref_d2) * inv_ref_d2;
  return edge_prune_threshold == nullptr ||
         !edge_scan_cannot_beat_threshold(candidate, edge_prune_threshold,
                                          phase);
}

#if TITAN_SIMD_X86 && TITAN_RUNTIME_TARGET_ATTRS
TITAN_TARGET_AVX2
bool accumulate_edge_scan_avx2(
    PbctopoCutCandidate &candidate, const PbctopoEdgeScanSoA &edge_scan_soa,
    const vector<Vec3> &placed_points,
    const PbctopoStage1EdgeContext *edge_context,
    const PbctopoCutCandidate *edge_prune_threshold,
    PbctopoRefineScorePhase phase) {
  if (placed_points.empty())
    return true;

  auto &scratch = pbctopo_edge_scan_scratch();
  sync_edge_scan_coords_soa(placed_points, scratch);
  const double *xs = scratch.xs.data();
  const double *ys = scratch.ys.data();
  const double *zs = scratch.zs.data();

  const auto process_block = [&](const size_t *edge_indices,
                                 size_t block_count) -> bool {
    alignas(32) long long atom_a_lane[4]{};
    alignas(32) long long atom_b_lane[4]{};
    alignas(32) double ref_lane[4]{};
    alignas(32) double inv_ref_lane[4]{};
    alignas(32) double cutoff_lane[4]{};
    for (size_t lane = 0; lane < block_count; ++lane) {
      const size_t edge_idx = edge_indices[lane];
      atom_a_lane[lane] = static_cast<long long>(edge_scan_soa.a[edge_idx]);
      atom_b_lane[lane] = static_cast<long long>(edge_scan_soa.b[edge_idx]);
      ref_lane[lane] = edge_scan_soa.ref_d2[edge_idx];
      inv_ref_lane[lane] = edge_scan_soa.inv_ref_d2[edge_idx];
      cutoff_lane[lane] = edge_scan_soa.break_cutoff_d2[edge_idx];
    }

    const __m256i gather_a = _mm256_set_epi64x(
        atom_a_lane[3], atom_a_lane[2], atom_a_lane[1], atom_a_lane[0]);
    const __m256i gather_b = _mm256_set_epi64x(
        atom_b_lane[3], atom_b_lane[2], atom_b_lane[1], atom_b_lane[0]);
    const __m256d ax = _mm256_i64gather_pd(xs, gather_a, 8);
    const __m256d ay = _mm256_i64gather_pd(ys, gather_a, 8);
    const __m256d az = _mm256_i64gather_pd(zs, gather_a, 8);
    const __m256d bx = _mm256_i64gather_pd(xs, gather_b, 8);
    const __m256d by = _mm256_i64gather_pd(ys, gather_b, 8);
    const __m256d bz = _mm256_i64gather_pd(zs, gather_b, 8);
    const __m256d dx = _mm256_sub_pd(ax, bx);
    const __m256d dy = _mm256_sub_pd(ay, by);
    const __m256d dz = _mm256_sub_pd(az, bz);
    const __m256d d2 = _mm256_add_pd(
        _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)),
        _mm256_mul_pd(dz, dz));
    alignas(32) double d2_lane[4];
    _mm256_store_pd(d2_lane, d2);

    for (size_t lane = 0; lane < block_count; ++lane) {
      if (!update_cut_candidate_with_edge_distance(
              candidate, d2_lane[lane], ref_lane[lane], inv_ref_lane[lane],
              cutoff_lane[lane], edge_prune_threshold, phase)) {
        return false;
      }
    }
    return true;
  };

  if (edge_context != nullptr) {
    size_t block[4]{};
    size_t block_count = 0;
    for (size_t edge_idx : edge_context->variable_edge_indices) {
      if (edge_idx >= edge_scan_soa.a.size())
        continue;
      block[block_count++] = edge_idx;
      if (block_count == 4) {
        if (!process_block(block, block_count))
          return false;
        block_count = 0;
      }
    }
    if (block_count != 0 && !process_block(block, block_count))
      return false;
    return true;
  }

  size_t edge_idx = 0;
  size_t block[4]{};
  while (edge_idx + 4 <= edge_scan_soa.a.size()) {
    block[0] = edge_idx + 0;
    block[1] = edge_idx + 1;
    block[2] = edge_idx + 2;
    block[3] = edge_idx + 3;
    if (!process_block(block, 4))
      return false;
    edge_idx += 4;
  }
  for (; edge_idx < edge_scan_soa.a.size(); ++edge_idx) {
    const size_t atom_a = edge_scan_soa.a[edge_idx];
    const size_t atom_b = edge_scan_soa.b[edge_idx];
    const Vec3 delta = placed_points[atom_a] - placed_points[atom_b];
    const double d2 = norm2(delta);
    if (!update_cut_candidate_with_edge_distance(
            candidate, d2, edge_scan_soa.ref_d2[edge_idx],
            edge_scan_soa.inv_ref_d2[edge_idx],
            edge_scan_soa.break_cutoff_d2[edge_idx], edge_prune_threshold,
            phase)) {
      return false;
    }
  }
  return true;
}
#endif

bool accumulate_edge_scan(
    PbctopoCutCandidate &candidate, const vector<PbctopoLocalEdge> &local_edges,
    const vector<Vec3> &placed_points,
    const PbctopoStage1EdgeContext *edge_context,
    const PbctopoEdgeScanSoA *edge_scan_soa,
    const PbctopoCutCandidate *edge_prune_threshold,
    PbctopoRefineScorePhase phase) {
  const bool use_soa =
      edge_scan_soa != nullptr && edge_scan_soa->a.size() == local_edges.size() &&
      edge_scan_soa->b.size() == local_edges.size() &&
      edge_scan_soa->ref_d2.size() == local_edges.size() &&
      edge_scan_soa->inv_ref_d2.size() == local_edges.size() &&
      edge_scan_soa->break_cutoff_d2.size() == local_edges.size();

#if TITAN_SIMD_X86 && TITAN_RUNTIME_TARGET_ATTRS
  const bool simd_enabled = use_soa &&
                            titan_simd_dispatch::runtime_has_avx2() &&
                            local_edges.size() >= 8;
  if (simd_enabled) {
    return accumulate_edge_scan_avx2(candidate, *edge_scan_soa, placed_points,
                                     edge_context, edge_prune_threshold,
                                     phase);
  }
#endif

  const auto accumulate_edge = [&](size_t edge_idx) {
    size_t atom_a = 0;
    size_t atom_b = 0;
    double ref_d2 = 0.0;
    double inv_ref_d2 = 0.0;
    double break_cutoff_d2 = 0.0;
    if (use_soa) {
      atom_a = edge_scan_soa->a[edge_idx];
      atom_b = edge_scan_soa->b[edge_idx];
      ref_d2 = edge_scan_soa->ref_d2[edge_idx];
      inv_ref_d2 = edge_scan_soa->inv_ref_d2[edge_idx];
      break_cutoff_d2 = edge_scan_soa->break_cutoff_d2[edge_idx];
    } else {
      const auto &edge = local_edges[edge_idx];
      atom_a = edge.a;
      atom_b = edge.b;
      ref_d2 = edge.ref_d2;
      inv_ref_d2 = edge.inv_ref_d2;
      break_cutoff_d2 = edge.break_cutoff_d2;
    }
    const Vec3 delta = placed_points[atom_a] - placed_points[atom_b];
    const double d2 = norm2(delta);
    return update_cut_candidate_with_edge_distance(
        candidate, d2, ref_d2, inv_ref_d2, break_cutoff_d2,
        edge_prune_threshold, phase);
  };

  if (edge_context != nullptr) {
    for (size_t edge_idx : edge_context->variable_edge_indices) {
      if (edge_idx >= local_edges.size())
        continue;
      if (!accumulate_edge(edge_idx))
        return false;
    }
    return true;
  }

  for (size_t edge_idx = 0; edge_idx < local_edges.size(); ++edge_idx) {
    if (!accumulate_edge(edge_idx))
      return false;
  }
  return true;
}

bool edge_scan_cannot_beat_threshold(const PbctopoCutCandidate &partial,
                                     const PbctopoCutCandidate *threshold,
                                     PbctopoRefineScorePhase phase) {
  if (threshold == nullptr || !threshold->valid)
    return false;
  const int prefix_cmp =
      compare_candidate_prefix_before_edge_stretch(partial, *threshold, phase);
  if (prefix_cmp > 0)
    return true;
  if (prefix_cmp < 0)
    return false;
  return partial.edge_stretch > threshold->edge_stretch + 1.0e-6;
}

void append_interface_anchor_edges(
    const vector<PbctopoAtom> &atoms, const vector<Vec3> &placed_points,
    const vector<size_t> &atom_component_index,
    const unordered_set<uint64_t> &edge_set, double contact_cutoff,
    vector<PbctopoLocalEdge> &edges) {
  if (atoms.size() != placed_points.size() ||
      atom_component_index.size() != atoms.size() || atoms.size() < 2 ||
      contact_cutoff <= kDistanceEps) {
    return;
  }

  struct AnchorCandidate {
    size_t a = 0;
    size_t b = 0;
    double d2 = numeric_limits<double>::infinity();
  };

  unordered_map<ClashCellKey, vector<size_t>, ClashCellKeyHash> grid;
  grid.reserve(max<size_t>(64, atoms.size() / 4));
  unordered_map<uint64_t, vector<AnchorCandidate>> pair_candidates;
  pair_candidates.reserve(16);
  const double inv_grid = 1.0 / contact_cutoff;
  const double contact_cutoff2 = contact_cutoff * contact_cutoff;

  for (size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
    const auto &point = placed_points[atom_idx];
    const int cx = static_cast<int>(floor(point.x * inv_grid));
    const int cy = static_cast<int>(floor(point.y * inv_grid));
    const int cz = static_cast<int>(floor(point.z * inv_grid));
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          const auto it = grid.find({cx + dx, cy + dy, cz + dz});
          if (it == grid.end())
            continue;
          for (size_t other_idx : it->second) {
            if (atoms[atom_idx].owner == atoms[other_idx].owner)
              continue;
            const size_t comp_a = atom_component_index[atom_idx];
            const size_t comp_b = atom_component_index[other_idx];
            if (comp_a == numeric_limits<size_t>::max() ||
                comp_b == numeric_limits<size_t>::max() || comp_a == comp_b) {
              continue;
            }
            if (edge_set.find(atom_pair_key(atom_idx, other_idx)) != edge_set.end())
              continue;
            const auto &other = placed_points[other_idx];
            const double ddx = point.x - other.x;
            const double ddy = point.y - other.y;
            const double ddz = point.z - other.z;
            const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (d2 > contact_cutoff2)
              continue;
            const size_t lo = min(comp_a, comp_b);
            const size_t hi = max(comp_a, comp_b);
            const uint64_t pair_key = (static_cast<uint64_t>(lo) << 32) |
                                      static_cast<uint64_t>(hi);
            pair_candidates[pair_key].push_back({other_idx, atom_idx, d2});
          }
        }
      }
    }
    grid[{cx, cy, cz}].push_back(atom_idx);
  }

  constexpr size_t kMaxAnchorsPerPair = 4;
  constexpr size_t kMinContactsPerPair = 2;
  for (auto &entry : pair_candidates) {
    auto &candidates = entry.second;
    if (candidates.size() < kMinContactsPerPair)
      continue;
    sort(candidates.begin(), candidates.end(),
         [](const AnchorCandidate &lhs, const AnchorCandidate &rhs) {
           const int metric_cmp =
               compare_pbctopo_double_metric(lhs.d2, rhs.d2, 1.0e-9);
           if (metric_cmp != 0)
             return metric_cmp < 0;
           if (lhs.a != rhs.a)
             return lhs.a < rhs.a;
           return lhs.b < rhs.b;
         });
    unordered_set<size_t> used_a;
    unordered_set<size_t> used_b;
    size_t added = 0;
    for (const auto &candidate : candidates) {
      if (used_a.find(candidate.a) != used_a.end() ||
          used_b.find(candidate.b) != used_b.end()) {
        continue;
      }
      const double ref_d2 = max(candidate.d2, 1.0e-6);
      edges.push_back({candidate.a, candidate.b, ref_d2, 1.0 / ref_d2,
                       max(36.0, ref_d2 * 4.0)});
      used_a.insert(candidate.a);
      used_b.insert(candidate.b);
      ++added;
      if (added >= kMaxAnchorsPerPair)
        break;
    }
  }
}

PbctopoCutCandidate rescore_cut_candidate(
    const vector<PbctopoLocalEdge> &scoring_edges,
    const vector<Vec3> &placed_points, const PbctopoCutCandidate &source) {
  if (!source.valid)
    return source;
  const PointStats stats = build_point_stats(placed_points);
  PbctopoCutCandidate rescored = score_cut_candidate(
      scoring_edges, placed_points, source.cut_x, source.cut_y, source.cut_z,
      source.gap_score, stats, source.shift_norm2);
  rescored.min_pair_d2 = source.min_pair_d2;
  rescored.clash_atom_a = source.clash_atom_a;
  rescored.clash_atom_b = source.clash_atom_b;
  rescored.continuity_path = source.continuity_path;
  rescored.continuity_max_d2 = source.continuity_max_d2;
  rescored.anchor_mst2 = source.anchor_mst2;
  rescored.anchor_max_d2 = source.anchor_max_d2;
  rescored.anchor_graph = source.anchor_graph;
  rescored.ultra_clash = source.ultra_clash;
  rescored.valid = source.valid;
  return rescored;
}

vector<PbctopoAtom> unwrap_atoms_by_local_graph(
    const vector<PbctopoAtom> &wrapped_atoms,
    const vector<PbctopoLocalEdge> &edges, box Box) {
  vector<PbctopoAtom> unwrapped = wrapped_atoms;
  if (wrapped_atoms.empty())
    return unwrapped;
  const PbctopoLatticeMetric lattice(Box);

  vector<vector<size_t>> adjacency(wrapped_atoms.size());
  adjacency.reserve(wrapped_atoms.size());
  for (const auto &edge : edges) {
    adjacency[edge.a].push_back(edge.b);
    adjacency[edge.b].push_back(edge.a);
  }

  vector<char> seen(wrapped_atoms.size(), 0);
  vector<size_t> stack;
  for (size_t atom = 0; atom < wrapped_atoms.size(); ++atom) {
    if (seen[atom] != 0)
      continue;
    seen[atom] = 1;
    stack.clear();
    stack.push_back(atom);
    unwrapped[atom].sx = wrapped_atoms[atom].sx;
    unwrapped[atom].sy = wrapped_atoms[atom].sy;
    unwrapped[atom].sz = wrapped_atoms[atom].sz;

    while (!stack.empty()) {
      const size_t current = stack.back();
      stack.pop_back();
      for (size_t nb : adjacency[current]) {
        if (seen[nb] != 0)
          continue;
        const auto mic = lattice.nearest_image_scaled(
            wrapped_atoms[nb].sx - wrapped_atoms[current].sx,
            wrapped_atoms[nb].sy - wrapped_atoms[current].sy,
            wrapped_atoms[nb].sz - wrapped_atoms[current].sz);
        unwrapped[nb].sx = unwrapped[current].sx + mic.sx;
        unwrapped[nb].sy = unwrapped[current].sy + mic.sy;
        unwrapped[nb].sz = unwrapped[current].sz + mic.sz;
        seen[nb] = 1;
        stack.push_back(nb);
      }
    }
  }

  for (auto &atom : unwrapped) {
    scaled_to_cartesian(Box, atom.sx, atom.sy, atom.sz, atom.x, atom.y, atom.z);
  }
  return unwrapped;
}

PbctopoComponentCandidate build_component_candidate_seed(
    const PbctopoCutCandidate &source, const vector<vector<size_t>> &components,
    const vector<Vec3> &placed_points) {
  PbctopoComponentCandidate seed;
  if (!source.valid || placed_points.empty())
    return seed;
  seed.continuity_path = source.continuity_path;
  seed.continuity_max_d2 = source.continuity_max_d2;
  seed.anchor_mst2 = source.anchor_mst2;
  seed.anchor_max_d2 = source.anchor_max_d2;
  seed.anchor_graph = source.anchor_graph;
  seed.principal_rg2_max = source.principal_rg2_max;
  seed.principal_rg2_mid = source.principal_rg2_mid;
  seed.principal_rg2_min = source.principal_rg2_min;
  seed.bbox_volume = source.bbox_volume;
  seed.rg2 = source.rg2;
  seed.shape_kappa2 = source.shape_kappa2;
  seed.shape_asphericity = source.shape_asphericity;
  seed.shape_acylindricity = source.shape_acylindricity;
  seed.shift_norm2 = source.shift_norm2;
  seed.min_pair_d2 = source.min_pair_d2;
  seed.clash_atom_a = source.clash_atom_a;
  seed.clash_atom_b = source.clash_atom_b;
  seed.ultra_clash = source.ultra_clash;
  vector<Vec3> component_centers;
  component_centers.reserve(components.size());
  for (const auto &component : components) {
    if (component.empty())
      continue;
    Vec3 center{};
    for (size_t atom_idx : component)
      center = center + placed_points[atom_idx];
    component_centers.push_back(
        center * (1.0 / static_cast<double>(component.size())));
  }
  seed.component_rg2 = compute_rg2(component_centers);
  const PointStats component_stats = build_point_stats(component_centers);
  const auto component_principal_rg2 = compute_principal_rg2(component_stats);
  const auto component_shape_metrics =
      compute_gyration_shape_metrics(component_principal_rg2);
  seed.component_kappa2 = component_shape_metrics.kappa2;
  seed.component_asphericity = component_shape_metrics.asphericity;
  seed.component_acylindricity = component_shape_metrics.acylindricity;
  seed.component_mst2 = compute_mst2(component_centers);
  seed.valid = true;
  return seed;
}

PbctopoComponentCandidate score_component_assignment(
    const vector<PbctopoAtom> &atoms, const vector<Vec3> &placed_points,
    const vector<vector<size_t>> &components, double min_pair_d2,
    size_t clash_atom_a, size_t clash_atom_b,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices, double continuity_cutoff,
    double shift_norm2, const PbctopoComponentCandidate *incumbent = nullptr,
    PbctopoTimingStats *timing = nullptr,
    PbctopoAnchorBackend anchor_backend = PbctopoAnchorBackend::Auto,
    const PbctopoAnchorGraphContext *anchor_context = nullptr) {
  using clock = chrono::steady_clock;
  const auto started = clock::now();
  auto finish = [&](PbctopoComponentCandidate candidate) {
    if (timing != nullptr) {
      ++timing->component_score_calls;
      timing->component_score_s +=
          chrono::duration<double>(clock::now() - started).count();
    }
    return candidate;
  };
  PbctopoComponentCandidate candidate;
  const bool debug_small_search = (components.size() <= 2);
  if (atoms.empty() || placed_points.size() != atoms.size())
    return finish(candidate);
  candidate.shift_norm2 = shift_norm2;
  candidate.min_pair_d2 = min_pair_d2;
  candidate.clash_atom_a = clash_atom_a;
  candidate.clash_atom_b = clash_atom_b;
  PbctopoCutCandidate incumbent_cut_storage;
  const PbctopoCutCandidate *incumbent_cut = nullptr;
  if (incumbent != nullptr && incumbent->valid && !incumbent->ultra_clash) {
    incumbent_cut_storage = component_candidate_to_cut_candidate(*incumbent);
    incumbent_cut = &incumbent_cut_storage;
  }

  auto stage_started = clock::now();
  if (!continuity_metrics_ok(placed_points, continuity_pairs, continuity_cutoff,
                             candidate.continuity_path,
                             candidate.continuity_max_d2, incumbent_cut)) {
    if (timing != nullptr) {
      ++timing->component_score_continuity_rejects;
      timing->component_score_continuity_s +=
          chrono::duration<double>(clock::now() - stage_started).count();
    }
    if (debug_small_search) {
      append_pbctopo_log_line(
          "[pbctopo-debug] component leaf rejected by continuity: max_d2=" +
          to_string(candidate.continuity_max_d2) +
          ", path=" + to_string(candidate.continuity_path));
    }
    return finish(candidate);
  }
  if (timing != nullptr) {
    timing->component_score_continuity_s +=
        chrono::duration<double>(clock::now() - stage_started).count();
  }
  if (component_candidate_cannot_beat_incumbent_prefix(
          candidate, incumbent, PbctopoComponentCandidatePrefixStage::Continuity)) {
    if (timing != nullptr)
      ++timing->component_score_prefix_prunes_after_continuity;
    return finish(candidate);
  }

  stage_started = clock::now();
  if (!anchor_graph_metrics_enhanced(
          placed_points, anchor_indices, continuity_pairs, candidate.anchor_mst2,
          candidate.anchor_max_d2, candidate.anchor_graph, incumbent_cut,
          PbctopoRefineScorePhase::Exact, anchor_backend, nullptr,
          anchor_context)) {
    if (timing != nullptr) {
      ++timing->component_score_anchor_rejects;
      timing->component_score_anchor_s +=
          chrono::duration<double>(clock::now() - stage_started).count();
    }
    if (debug_small_search) {
      append_pbctopo_log_line(
          "[pbctopo-debug] component leaf rejected by anchor-graph: max_d2=" +
          to_string(candidate.anchor_max_d2) +
          ", mst2=" + to_string(candidate.anchor_mst2));
    }
    return finish(candidate);
  }
  if (timing != nullptr) {
    timing->component_score_anchor_s +=
        chrono::duration<double>(clock::now() - stage_started).count();
  }
  if (component_candidate_cannot_beat_incumbent_prefix(
          candidate, incumbent, PbctopoComponentCandidatePrefixStage::Anchor)) {
    if (timing != nullptr)
      ++timing->component_score_prefix_prunes_after_anchor;
    return finish(candidate);
  }

  stage_started = clock::now();
  const PointStats stats = build_point_stats(placed_points);
  if (stats.count == 0)
    return finish(candidate);
  const array<double, 3> principal_rg2 = compute_principal_rg2(stats);
  candidate.principal_rg2_max = principal_rg2[0];
  candidate.principal_rg2_mid = principal_rg2[1];
  candidate.principal_rg2_min = principal_rg2[2];
  const auto shape_metrics = compute_gyration_shape_metrics(principal_rg2);
  candidate.shape_kappa2 = shape_metrics.kappa2;
  candidate.shape_asphericity = shape_metrics.asphericity;
  candidate.shape_acylindricity = shape_metrics.acylindricity;
  candidate.rg2 = compute_rg2_from_stats(stats);
  candidate.bbox_volume = compute_bbox_volume(stats);
  if (timing != nullptr) {
    timing->component_score_stats_principal_s +=
        chrono::duration<double>(clock::now() - stage_started).count();
  }
  if (component_candidate_cannot_beat_incumbent_prefix(
          candidate, incumbent, PbctopoComponentCandidatePrefixStage::Principal)) {
    if (timing != nullptr)
      ++timing->component_score_prefix_prunes_after_principal;
    return finish(candidate);
  }

  stage_started = clock::now();
  vector<Vec3> component_centers;
  component_centers.reserve(components.size());
  for (const auto &component : components) {
    if (component.empty())
      continue;
    Vec3 center{};
    for (size_t atom_idx : component)
      center = center + placed_points[atom_idx];
    component_centers.push_back(
        center * (1.0 / static_cast<double>(component.size())));
  }
  candidate.component_rg2 = compute_rg2(component_centers);
  const PointStats component_stats = build_point_stats(component_centers);
  const auto component_principal_rg2 = compute_principal_rg2(component_stats);
  const auto component_shape_metrics =
      compute_gyration_shape_metrics(component_principal_rg2);
  candidate.component_kappa2 = component_shape_metrics.kappa2;
  candidate.component_asphericity = component_shape_metrics.asphericity;
  candidate.component_acylindricity = component_shape_metrics.acylindricity;
  candidate.component_mst2 = compute_mst2(component_centers);
  if (timing != nullptr) {
    timing->component_score_component_metrics_s +=
        chrono::duration<double>(clock::now() - stage_started).count();
  }
  if (component_candidate_cannot_beat_incumbent_prefix(
          candidate, incumbent,
          PbctopoComponentCandidatePrefixStage::ComponentMetrics)) {
    if (timing != nullptr)
      ++timing->component_score_prefix_prunes_after_component_metrics;
    return finish(candidate);
  }
  if (component_candidate_cannot_beat_incumbent_prefix(
          candidate, incumbent, PbctopoComponentCandidatePrefixStage::Final)) {
    if (timing != nullptr)
      ++timing->component_score_prefix_prunes_after_final;
    return finish(candidate);
  }
  candidate.valid = true;
  if (debug_small_search) {
    append_pbctopo_log_line(
        "[pbctopo-debug] component leaf valid: cont_max_d2=" +
        to_string(candidate.continuity_max_d2) +
        ", cont_path=" + to_string(candidate.continuity_path) +
        ", anchor_max_d2=" + to_string(candidate.anchor_max_d2) +
        ", anchor_mst2=" + to_string(candidate.anchor_mst2) +
        ", princ_max=" + to_string(candidate.principal_rg2_max) +
        ", comp_mst2=" + to_string(candidate.component_mst2) +
        ", kappa2=" + to_string(candidate.shape_kappa2) +
        ", asph=" + to_string(candidate.shape_asphericity));
  }
  return finish(candidate);
}

PbctopoComponentSearchCallbacks make_component_search_callbacks() {
  PbctopoComponentSearchCallbacks callbacks;
  callbacks.build_seed = build_component_candidate_seed;
  callbacks.score_assignment = score_component_assignment;
  callbacks.better_candidate = better_component_candidate;
  return callbacks;
}

bool solve_large_component_merge_mode(
    const vector<PbctopoAtom> &atoms, box Box,
    const vector<vector<size_t>> &raw_components,
    const PbctopoMergedComponentPlan &merge_plan,
    const vector<PbctopoShift> &shifts,
    const vector<PbctopoLocalEdge> &local_edges,
    const unordered_set<uint64_t> &edge_set, vector<Vec3> &placed_points,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices,
    const PbctopoRuntimeCutoffs &cutoffs, PbctopoCutCandidate &best_cut,
    double &shift_norm2, PbctopoTimingStats *timing,
    PbctopoAnchorBackend anchor_backend,
    const vector<Int3> *warm_component_offsets,
    PbctopoWarmState *warm_state,
    const PbctopoCutCandidate *warm_reference_candidate,
    bool allow_exact_search_on_merged,
    vector<Int3> *best_raw_offsets_out = nullptr,
    const PbctopoAnchorGraphContext *anchor_context = nullptr) {
  best_cut = {};
  shift_norm2 = 0.0;
  placed_points.clear();
  if (merge_plan.merged_components.empty())
    return false;

  const vector<PbctopoAtom> merged_atoms =
      apply_merged_component_internal_offsets(atoms, raw_components,
                                              merge_plan, Box);
  if (timing != nullptr) {
    ++timing->component_large_mode_internal_offsets_preserved;
    timing->component_large_mode_internal_offset_nonzero +=
        count_if(merge_plan.raw_internal_offsets.begin(),
                 merge_plan.raw_internal_offsets.end(),
                 [](const Int3 &offset) { return !is_zero_delta(offset); });
  }

  vector<Int3> merged_warm_offsets;
  const vector<Int3> *merged_warm_offsets_ptr =
      (warm_component_offsets != nullptr &&
       derive_merged_warm_offsets(*warm_component_offsets, merge_plan,
                                  merged_warm_offsets))
          ? &merged_warm_offsets
          : nullptr;

  vector<Vec3> best_points_local;
  vector<Int3> best_raw_offsets;
  auto consider_raw_offsets = [&](const vector<Int3> &raw_offsets) {
    vector<Int3> refine_offsets = raw_offsets;
    vector<Vec3> refine_points;
    PbctopoCutCandidate candidate = improve_component_offsets_locally(
        atoms, raw_components, Box, shifts, local_edges, edge_set,
        continuity_pairs, anchor_indices, refine_offsets, cutoffs.clash_cutoff,
        cutoffs.continuity_cutoff, refine_points,
        make_component_refine_callbacks(), anchor_backend, timing, warm_state,
        warm_reference_candidate, anchor_context);
    if (candidate.valid && !candidate.ultra_clash &&
        refine_points.size() == atoms.size() &&
        (!best_cut.valid || better_cut_candidate(candidate, best_cut))) {
      best_cut = candidate;
      best_points_local = refine_points;
      best_raw_offsets = refine_offsets;
    }
  };

  const vector<PbctopoGraphSeedSolveResult> graph_seed_results =
      build_large_component_graph_seed_candidates(
          merged_atoms, Box, merge_plan.merged_components, shifts,
          continuity_pairs, anchor_indices, cutoffs, merged_warm_offsets_ptr, timing,
          anchor_backend, anchor_context);
  if (timing != nullptr)
    timing->component_large_mode_graph_seed_count += graph_seed_results.size();

  const PbctopoGraphSeedSolveResult *selected_seed = nullptr;
  for (const auto &seed_result : graph_seed_results) {
    if (!seed_result.valid)
      continue;
    if (selected_seed == nullptr ||
        graph_seed_quick_score_better(seed_result.quick_score,
                                      selected_seed->quick_score)) {
      selected_seed = &seed_result;
    }
  }
  if (selected_seed == nullptr && !graph_seed_results.empty())
    selected_seed = &graph_seed_results.front();

  if (selected_seed != nullptr &&
      selected_seed->offsets.size() == merge_plan.merged_components.size()) {
    if (timing != nullptr) {
      ++timing->component_large_mode_graph_seed_refined;
      if (string(selected_seed->seed_name) != "baseline")
        ++timing->component_large_mode_graph_seed_nonbaseline_win;
    }
    consider_raw_offsets(expand_merged_offsets_to_raw_components(
        selected_seed->offsets, Box, raw_components, merge_plan));
  }

  if (allow_exact_search_on_merged) {
    vector<Vec3> merged_search_points;
    const PbctopoCutCandidate *search_incumbent =
        (best_cut.valid && !best_cut.ultra_clash) ? &best_cut : warm_reference_candidate;
    const vector<Vec3> *search_incumbent_points =
        (best_cut.valid && !best_cut.ultra_clash &&
         best_points_local.size() == atoms.size())
            ? &best_points_local
            : nullptr;
    if (timing != nullptr)
      ++timing->component_search_frames;
    PbctopoComponentCandidate merged_candidate = search_component_shift_candidate(
        merged_atoms, Box, merge_plan.merged_components, shifts, continuity_pairs,
        anchor_indices, cutoffs, merged_search_points,
        make_component_search_callbacks(), merged_warm_offsets_ptr,
        search_incumbent, search_incumbent_points, timing, anchor_backend,
        anchor_context);
    if (merged_candidate.valid && !merged_candidate.ultra_clash &&
        merged_search_points.size() == atoms.size()) {
      const vector<Int3> merged_search_offsets = extract_component_offsets_from_points(
          merged_atoms, merge_plan.merged_components, Box,
          merged_search_points);
      if (merged_search_offsets.size() == merge_plan.merged_components.size()) {
        consider_raw_offsets(expand_merged_offsets_to_raw_components(
            merged_search_offsets, Box, raw_components, merge_plan));
      }
    }
  }

  if (!best_cut.valid || best_points_local.size() != atoms.size())
    return false;
  placed_points = best_points_local;
  shift_norm2 = best_cut.shift_norm2;
  if (best_raw_offsets_out != nullptr)
    *best_raw_offsets_out = best_raw_offsets;
  return true;
}

double clamp_virtual_box_center(double center) {
  if (center < 0.5)
    return 0.5;
  if (center > 2.5)
    return 2.5;
  return center;
}

bool find_virtual_box_shift(double wrapped_coord, double center, int &shift_idx) {
  const double low = center - 0.5;
  const double high = center + 0.5;
  for (int candidate = 0; candidate <= 2; ++candidate) {
    const double coord = wrapped_coord + static_cast<double>(candidate);
    if (coord + kDistanceEps >= low && coord < high - kDistanceEps) {
      shift_idx = candidate;
      return true;
    }
  }
  const double mid = 0.5 * (low + high);
  double best_dist = numeric_limits<double>::infinity();
  int best_candidate = 0;
  for (int candidate = 0; candidate <= 2; ++candidate) {
    const double coord = wrapped_coord + static_cast<double>(candidate);
    const double dist = fabs(coord - mid);
    if (dist < best_dist) {
      best_dist = dist;
      best_candidate = candidate;
    }
  }
  shift_idx = best_candidate;
  return true;
}

void push_unique_xy_seed(vector<PbctopoSeed2D> &best, const PbctopoSeed2D &seed,
                         size_t max_keep) {
  for (const auto &existing : best) {
    if (fabs(existing.cx - seed.cx) <= 1.0e-9 &&
        fabs(existing.cy - seed.cy) <= 1.0e-9) {
      return;
    }
  }
  best.push_back(seed);
  sort(best.begin(), best.end(), [](const PbctopoSeed2D &lhs,
                                    const PbctopoSeed2D &rhs) {
    int cmp = compare_pbctopo_double_metric(lhs.score, rhs.score, 1.0e-9);
    if (cmp != 0)
      return cmp > 0;
    cmp = compare_pbctopo_double_metric(lhs.cx, rhs.cx, 1.0e-9);
    if (cmp != 0)
      return cmp < 0;
    return compare_pbctopo_double_metric(lhs.cy, rhs.cy, 1.0e-9) < 0;
  });
  if (best.size() > max_keep)
    best.resize(max_keep);
}

void push_unique_axis_seed(vector<PbctopoCutOption> &best,
                           const PbctopoCutOption &seed, size_t max_keep) {
  for (const auto &existing : best) {
    if (fabs(existing.cut - seed.cut) <= 1.0e-9)
      return;
  }
  best.push_back(seed);
  sort(best.begin(), best.end(),
       [](const PbctopoCutOption &lhs, const PbctopoCutOption &rhs) {
         const int gap_cmp =
             compare_pbctopo_double_metric(lhs.gap, rhs.gap, 1.0e-9);
         if (gap_cmp != 0)
           return gap_cmp > 0;
         return compare_pbctopo_double_metric(lhs.cut, rhs.cut, 1.0e-9) < 0;
       });
  if (best.size() > max_keep)
    best.resize(max_keep);
}

void push_unique_xyz_seed(vector<PbctopoSeed3D> &best, const PbctopoSeed3D &seed,
                          size_t max_keep) {
  for (const auto &existing : best) {
    if (fabs(existing.cx - seed.cx) <= 1.0e-9 &&
        fabs(existing.cy - seed.cy) <= 1.0e-9 &&
        fabs(existing.cz - seed.cz) <= 1.0e-9) {
      return;
    }
  }
  best.push_back(seed);
  sort(best.begin(), best.end(), [](const PbctopoSeed3D &lhs,
                                    const PbctopoSeed3D &rhs) {
    int cmp = compare_pbctopo_double_metric(lhs.score, rhs.score, 1.0e-9);
    if (cmp != 0)
      return cmp > 0;
    cmp = compare_pbctopo_double_metric(lhs.cx, rhs.cx, 1.0e-9);
    if (cmp != 0)
      return cmp < 0;
    cmp = compare_pbctopo_double_metric(lhs.cy, rhs.cy, 1.0e-9);
    if (cmp != 0)
      return cmp < 0;
    return compare_pbctopo_double_metric(lhs.cz, rhs.cz, 1.0e-9) < 0;
  });
  if (best.size() > max_keep)
    best.resize(max_keep);
}

double projection_axis_coord(const PbctopoAtom &atom, char axis) {
  if (axis == 'x')
    return atom.sx;
  if (axis == 'y')
    return atom.sy;
  return atom.sz;
}

vector<double> build_projection_density_grid(const vector<PbctopoAtom> &atoms,
                                             size_t mesh_bins, double radius,
                                             char axis_a, char axis_b) {
  vector<double> scores(mesh_bins * mesh_bins, 0.0);
  if (atoms.empty() || mesh_bins == 0)
    return scores;

  const double span = 2.0;
  const double step =
      (mesh_bins == 1) ? 0.0 : span / static_cast<double>(mesh_bins - 1);
  auto clamp_index = [mesh_bins](long long value) -> size_t {
    if (value < 0)
      return 0;
    const long long max_index = static_cast<long long>(mesh_bins) - 1;
    if (value > max_index)
      return static_cast<size_t>(max_index);
    return static_cast<size_t>(value);
  };

  for (const auto &atom : atoms) {
    const double wrapped_a = projection_axis_coord(atom, axis_a);
    const double wrapped_b = projection_axis_coord(atom, axis_b);
    for (int shift_a = 0; shift_a <= 2; ++shift_a) {
      const double shifted_a = wrapped_a + static_cast<double>(shift_a);
      if (shifted_a + radius < 0.5 || shifted_a - radius > 2.5)
        continue;
      const size_t ia_begin =
          (step <= kDistanceEps)
              ? 0
              : clamp_index(static_cast<long long>(
                    ceil((shifted_a - radius - 0.5) / step)));
      const size_t ia_end =
          (step <= kDistanceEps)
              ? 0
              : clamp_index(static_cast<long long>(
                    floor((shifted_a + radius - 0.5) / step)));
      for (int shift_b = 0; shift_b <= 2; ++shift_b) {
        const double shifted_b = wrapped_b + static_cast<double>(shift_b);
        if (shifted_b + radius < 0.5 || shifted_b - radius > 2.5)
          continue;
        const size_t ib_begin =
            (step <= kDistanceEps)
                ? 0
                : clamp_index(static_cast<long long>(
                      ceil((shifted_b - radius - 0.5) / step)));
        const size_t ib_end =
            (step <= kDistanceEps)
                ? 0
                : clamp_index(static_cast<long long>(
                      floor((shifted_b + radius - 0.5) / step)));
        for (size_t ia = ia_begin; ia <= ia_end; ++ia) {
          const double ca = 0.5 + step * static_cast<double>(ia);
          const double da = fabs(shifted_a - ca);
          if (da > radius)
            continue;
          for (size_t ib = ib_begin; ib <= ib_end; ++ib) {
            const double cb = 0.5 + step * static_cast<double>(ib);
            const double db = fabs(shifted_b - cb);
            if (db > radius)
              continue;
            scores[ia * mesh_bins + ib] += (radius - da) * (radius - db);
          }
        }
      }
    }
  }
  return scores;
}

double projection_grid_center(size_t mesh_bins, size_t idx) {
  if (mesh_bins <= 1)
    return 1.5;
  const double step = 2.0 / static_cast<double>(mesh_bins - 1);
  return 0.5 + step * static_cast<double>(idx);
}

bool projection_density_is_local_maximum(const vector<double> &scores,
                                         size_t mesh_bins, size_t ix,
                                         size_t iy) {
  const double score = scores[ix * mesh_bins + iy];
  if (!(score > 0.0))
    return false;

  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      if (dx == 0 && dy == 0)
        continue;
      const long long nx_ll = static_cast<long long>(ix) + dx;
      const long long ny_ll = static_cast<long long>(iy) + dy;
      if (nx_ll < 0 || ny_ll < 0)
        continue;
      const size_t nx = static_cast<size_t>(nx_ll);
      const size_t ny = static_cast<size_t>(ny_ll);
      if (nx >= mesh_bins || ny >= mesh_bins)
        continue;
      const double neighbor = scores[nx * mesh_bins + ny];
      if (neighbor > score + 1.0e-12)
        return false;
      if (fabs(neighbor - score) <= 1.0e-12 &&
          (nx < ix || (nx == ix && ny < iy))) {
        return false;
      }
    }
  }
  return true;
}

bool projection_seed2d_better(const PbctopoSeed2D &lhs,
                              const PbctopoSeed2D &rhs) {
  const int score_cmp =
      compare_pbctopo_double_metric(lhs.score, rhs.score, 1.0e-12);
  if (score_cmp != 0)
    return score_cmp > 0;
  if (lhs.ix != rhs.ix)
    return lhs.ix < rhs.ix;
  return lhs.iy < rhs.iy;
}

void collect_projection_plane_hotspots(const vector<double> &scores,
                                       size_t mesh_bins, size_t top_k,
                                       vector<PbctopoSeed2D> &hotspots) {
  hotspots.clear();
  if (mesh_bins == 0 || top_k == 0 || scores.size() != mesh_bins * mesh_bins)
    return;

  vector<PbctopoSeed2D> local_maxima;
  vector<PbctopoSeed2D> ranked_cells;
  local_maxima.reserve(mesh_bins * mesh_bins / 4 + 1);
  ranked_cells.reserve(mesh_bins * mesh_bins);

  for (size_t ix = 0; ix < mesh_bins; ++ix) {
    const double cx = projection_grid_center(mesh_bins, ix);
    for (size_t iy = 0; iy < mesh_bins; ++iy) {
      const double score = scores[ix * mesh_bins + iy];
      if (!(score > 0.0))
        continue;
      PbctopoSeed2D seed{cx, projection_grid_center(mesh_bins, iy), score, ix,
                         iy};
      ranked_cells.push_back(seed);
      if (projection_density_is_local_maximum(scores, mesh_bins, ix, iy))
        local_maxima.push_back(seed);
    }
  }

  sort(local_maxima.begin(), local_maxima.end(), projection_seed2d_better);
  sort(ranked_cells.begin(), ranked_cells.end(), projection_seed2d_better);

  auto append_from = [&](const vector<PbctopoSeed2D> &source,
                         long long min_bin_gap) {
    for (const auto &seed : source) {
      if (hotspots.size() >= top_k)
        break;
      bool too_close = false;
      for (const auto &existing : hotspots) {
        const long long dx =
            llabs(static_cast<long long>(existing.ix) -
                  static_cast<long long>(seed.ix));
        const long long dy =
            llabs(static_cast<long long>(existing.iy) -
                  static_cast<long long>(seed.iy));
        if (max(dx, dy) <= min_bin_gap) {
          too_close = true;
          break;
        }
      }
      if (!too_close)
        hotspots.push_back(seed);
    }
  };

  const long long nms_bin_gap = (mesh_bins >= 10) ? 1ll : 0ll;
  append_from(local_maxima, nms_bin_gap);
  if (hotspots.size() < min(top_k, static_cast<size_t>(4)))
    append_from(ranked_cells, nms_bin_gap);
  if (hotspots.size() < top_k)
    append_from(ranked_cells, 0);
}

double projection_grid_score(const vector<double> &scores, size_t mesh_bins,
                             size_t ia, size_t ib) {
  if (ia >= mesh_bins || ib >= mesh_bins || scores.size() != mesh_bins * mesh_bins)
    return 0.0;
  return scores[ia * mesh_bins + ib];
}

void append_dense_projection_volume_seeds(const vector<double> &xy_scores,
                                          const vector<double> &xz_scores,
                                          const vector<double> &yz_scores,
                                          size_t mesh_bins, size_t top_k,
                                          vector<PbctopoSeed3D> &seeds) {
  for (size_t ix = 0; ix < mesh_bins; ++ix) {
    const double cx = projection_grid_center(mesh_bins, ix);
    for (size_t iy = 0; iy < mesh_bins; ++iy) {
      const double cy = projection_grid_center(mesh_bins, iy);
      const double xy_score = projection_grid_score(xy_scores, mesh_bins, ix, iy);
      for (size_t iz = 0; iz < mesh_bins; ++iz) {
        const double score =
            xy_score + projection_grid_score(xz_scores, mesh_bins, ix, iz) +
            projection_grid_score(yz_scores, mesh_bins, iy, iz);
        push_unique_xyz_seed(
            seeds,
            {cx, cy, projection_grid_center(mesh_bins, iz), score, ix, iy, iz},
            top_k);
      }
    }
  }
}

void collect_projection_volume_seeds(const vector<PbctopoAtom> &atoms,
                                     size_t mesh_bins, size_t top_k,
                                     double radius,
                                     vector<PbctopoSeed3D> &seeds) {
  seeds.clear();
  if (atoms.empty() || mesh_bins == 0 || top_k == 0)
    return;

  const auto xy_scores =
      build_projection_density_grid(atoms, mesh_bins, radius, 'x', 'y');
  const auto xz_scores =
      build_projection_density_grid(atoms, mesh_bins, radius, 'x', 'z');
  const auto yz_scores =
      build_projection_density_grid(atoms, mesh_bins, radius, 'y', 'z');

  vector<PbctopoSeed2D> xy_hotspots;
  vector<PbctopoSeed2D> xz_hotspots;
  vector<PbctopoSeed2D> yz_hotspots;
  const size_t plane_top_k =
      min(mesh_bins * mesh_bins, max<size_t>(top_k * 2, static_cast<size_t>(12)));
  collect_projection_plane_hotspots(xy_scores, mesh_bins, plane_top_k,
                                    xy_hotspots);
  collect_projection_plane_hotspots(xz_scores, mesh_bins, plane_top_k,
                                    xz_hotspots);
  collect_projection_plane_hotspots(yz_scores, mesh_bins, plane_top_k,
                                    yz_hotspots);

  seeds.reserve(top_k);
  auto push_seed = [&](size_t ix, size_t iy, size_t iz, double score) {
    push_unique_xyz_seed(
        seeds,
        {projection_grid_center(mesh_bins, ix), projection_grid_center(mesh_bins, iy),
         projection_grid_center(mesh_bins, iz), score, ix, iy, iz},
        top_k);
  };

  for (const auto &xy_seed : xy_hotspots) {
    for (const auto &xz_seed : xz_hotspots) {
      if (xy_seed.ix != xz_seed.ix)
        continue;
      for (const auto &yz_seed : yz_hotspots) {
        if (xy_seed.iy != yz_seed.ix || xz_seed.iy != yz_seed.iy)
          continue;
        push_seed(xy_seed.ix, xy_seed.iy, xz_seed.iy,
                  xy_seed.score + xz_seed.score + yz_seed.score);
      }
    }
  }

  if (seeds.size() < top_k) {
    for (const auto &xy_seed : xy_hotspots) {
      for (const auto &xz_seed : xz_hotspots) {
        if (xy_seed.ix != xz_seed.ix)
          continue;
        push_seed(xy_seed.ix, xy_seed.iy, xz_seed.iy,
                  xy_seed.score + xz_seed.score +
                      projection_grid_score(yz_scores, mesh_bins, xy_seed.iy,
                                            xz_seed.iy));
      }
    }
  }

  if (seeds.size() < top_k) {
    for (const auto &xy_seed : xy_hotspots) {
      for (const auto &yz_seed : yz_hotspots) {
        if (xy_seed.iy != yz_seed.ix)
          continue;
        push_seed(xy_seed.ix, xy_seed.iy, yz_seed.iy,
                  xy_seed.score + yz_seed.score +
                      projection_grid_score(xz_scores, mesh_bins, xy_seed.ix,
                                            yz_seed.iy));
      }
    }
  }

  if (seeds.size() < top_k) {
    for (const auto &xz_seed : xz_hotspots) {
      for (const auto &yz_seed : yz_hotspots) {
        if (xz_seed.iy != yz_seed.iy)
          continue;
        push_seed(xz_seed.ix, yz_seed.ix, xz_seed.iy,
                  xz_seed.score + yz_seed.score +
                      projection_grid_score(xy_scores, mesh_bins, xz_seed.ix,
                                            yz_seed.ix));
      }
    }
  }

  if (seeds.size() < top_k)
    append_dense_projection_volume_seeds(xy_scores, xz_scores, yz_scores,
                                         mesh_bins, top_k, seeds);
}

[[maybe_unused]] vector<PbctopoSeed2D>
collect_xy_density_seeds(const vector<PbctopoAtom> &atoms, size_t mesh_bins,
                         size_t top_k, double radius) {
  vector<PbctopoSeed2D> seeds;
  if (atoms.empty() || mesh_bins == 0 || top_k == 0)
    return seeds;

  const double span = 2.0;
  const double step =
      (mesh_bins == 1) ? 0.0 : span / static_cast<double>(mesh_bins - 1);
  vector<double> scores(mesh_bins * mesh_bins, 0.0);
  auto clamp_index = [mesh_bins](long long value) -> size_t {
    if (value < 0)
      return 0;
    const long long max_index = static_cast<long long>(mesh_bins) - 1;
    if (value > max_index)
      return static_cast<size_t>(max_index);
    return static_cast<size_t>(value);
  };

  for (const auto &atom : atoms) {
    for (int shift_x = 0; shift_x <= 2; ++shift_x) {
      const double shifted_x = atom.sx + static_cast<double>(shift_x);
      if (shifted_x + radius < 0.5 || shifted_x - radius > 2.5)
        continue;
      const size_t ix_begin =
          (step <= kDistanceEps)
              ? 0
              : clamp_index(static_cast<long long>(
                    ceil((shifted_x - radius - 0.5) / step)));
      const size_t ix_end =
          (step <= kDistanceEps)
              ? 0
              : clamp_index(static_cast<long long>(
                    floor((shifted_x + radius - 0.5) / step)));
      for (int shift_y = 0; shift_y <= 2; ++shift_y) {
        const double shifted_y = atom.sy + static_cast<double>(shift_y);
        if (shifted_y + radius < 0.5 || shifted_y - radius > 2.5)
          continue;
        const size_t iy_begin =
            (step <= kDistanceEps)
                ? 0
                : clamp_index(static_cast<long long>(
                      ceil((shifted_y - radius - 0.5) / step)));
        const size_t iy_end =
            (step <= kDistanceEps)
                ? 0
                : clamp_index(static_cast<long long>(
                      floor((shifted_y + radius - 0.5) / step)));
        for (size_t ix = ix_begin; ix <= ix_end; ++ix) {
          const double cx = 0.5 + step * static_cast<double>(ix);
          const double dx = fabs(shifted_x - cx);
          if (dx > radius)
            continue;
          for (size_t iy = iy_begin; iy <= iy_end; ++iy) {
            const double cy = 0.5 + step * static_cast<double>(iy);
            const double dy = fabs(shifted_y - cy);
            if (dy > radius)
              continue;
            scores[ix * mesh_bins + iy] += (radius - dx) * (radius - dy);
          }
        }
      }
    }
  }

  for (size_t ix = 0; ix < mesh_bins; ++ix) {
    const double cx = 0.5 + step * static_cast<double>(ix);
    for (size_t iy = 0; iy < mesh_bins; ++iy) {
      const double cy = 0.5 + step * static_cast<double>(iy);
      push_unique_xy_seed(seeds, {cx, cy, scores[ix * mesh_bins + iy]}, top_k);
    }
  }
  return seeds;
}

[[maybe_unused]] vector<PbctopoCutOption>
collect_z_density_seeds(const vector<PbctopoAtom> &atoms, size_t mesh_bins,
                        size_t top_k, double radius) {
  vector<PbctopoCutOption> seeds;
  if (atoms.empty() || mesh_bins == 0 || top_k == 0)
    return seeds;

  const double span = 2.0;
  const double step =
      (mesh_bins == 1) ? 0.0 : span / static_cast<double>(mesh_bins - 1);
  vector<double> scores(mesh_bins, 0.0);
  auto clamp_index = [mesh_bins](long long value) -> size_t {
    if (value < 0)
      return 0;
    const long long max_index = static_cast<long long>(mesh_bins) - 1;
    if (value > max_index)
      return static_cast<size_t>(max_index);
    return static_cast<size_t>(value);
  };

  for (const auto &atom : atoms) {
    for (int shift_z = 0; shift_z <= 2; ++shift_z) {
      const double shifted_z = atom.sz + static_cast<double>(shift_z);
      if (shifted_z + radius < 0.5 || shifted_z - radius > 2.5)
        continue;
      const size_t iz_begin =
          (step <= kDistanceEps)
              ? 0
              : clamp_index(static_cast<long long>(
                    ceil((shifted_z - radius - 0.5) / step)));
      const size_t iz_end =
          (step <= kDistanceEps)
              ? 0
              : clamp_index(static_cast<long long>(
                    floor((shifted_z + radius - 0.5) / step)));
      for (size_t iz = iz_begin; iz <= iz_end; ++iz) {
        const double cz = 0.5 + step * static_cast<double>(iz);
        const double dz = fabs(shifted_z - cz);
        if (dz > radius)
          continue;
        scores[iz] += radius - dz;
      }
    }
  }

  for (size_t iz = 0; iz < mesh_bins; ++iz) {
    const double cz = 0.5 + step * static_cast<double>(iz);
    push_unique_axis_seed(seeds, {cz, scores[iz]}, top_k);
  }
  return seeds;
}

void push_top_virtual_box_candidate(vector<PbctopoCutCandidate> &top_candidates,
                                    const PbctopoCutCandidate &candidate,
                                    size_t max_keep) {
  if (!candidate.valid)
    return;
  top_candidates.push_back(candidate);
  sort(top_candidates.begin(), top_candidates.end(), better_cut_candidate);
  if (top_candidates.size() > max_keep)
    top_candidates.resize(max_keep);
}

[[maybe_unused]] bool
build_virtual_box_positions(const vector<PbctopoAtom> &atoms, box Box,
                            double center_x, double center_y, double center_z,
                            vector<Vec3> &placed_points, PointStats &stats,
                            double &shift_norm2) {
  placed_points.resize(atoms.size());
  stats = PointStats{};
  shift_norm2 = 0.0;

  const int global_shift_x = static_cast<int>(floor(center_x)) - 1;
  const int global_shift_y = static_cast<int>(floor(center_y)) - 1;
  const int global_shift_z = static_cast<int>(floor(center_z)) - 1;

  for (size_t idx = 0; idx < atoms.size(); ++idx) {
    const auto &atom = atoms[idx];
    int search_shift_x = -1;
    int search_shift_y = -1;
    int search_shift_z = -1;
    if (!find_virtual_box_shift(atom.sx, center_x, search_shift_x) ||
        !find_virtual_box_shift(atom.sy, center_y, search_shift_y) ||
        !find_virtual_box_shift(atom.sz, center_z, search_shift_z)) {
      return false;
    }

    const int actual_shift_x = search_shift_x - 1 - global_shift_x;
    const int actual_shift_y = search_shift_y - 1 - global_shift_y;
    const int actual_shift_z = search_shift_z - 1 - global_shift_z;

    Vec3 point{};
    scaled_to_cartesian(Box, atom.sx + static_cast<double>(actual_shift_x),
                        atom.sy + static_cast<double>(actual_shift_y),
                        atom.sz + static_cast<double>(actual_shift_z), point.x,
                        point.y, point.z);
    placed_points[idx] = point;
    stats.add(point.x, point.y, point.z);
    shift_norm2 += static_cast<double>(actual_shift_x * actual_shift_x +
                                       actual_shift_y * actual_shift_y +
                                       actual_shift_z * actual_shift_z);
  }

  return stats.count == atoms.size();
}

PbctopoCutCandidate evaluate_virtual_box_candidate(
    const vector<PbctopoAtom> &atoms, const vector<PbctopoLocalEdge> &local_edges,
    const unordered_set<uint64_t> &edge_set,
    const PbctopoTranslationCache &translation_cache, box Box, double cut_x,
    double cut_y, double cut_z, double gap_score, vector<Vec3> &placed_points,
    double clash_cutoff, const vector<PbctopoContinuityPair> *continuity_pairs,
    const vector<size_t> *anchor_indices,
    const vector<pair<int, int>> *contacting_owner_pairs,
    double continuity_cutoff, const PbctopoCutCandidate *incumbent,
    bool skip_clash, double frame_time_ps, const char *reject_context,
    const PbctopoAnchorGraphContext *anchor_context,
    PbctopoVirtualSeedSource seed_source,
    const vector<size_t> *atom_component_index) {
  PbctopoCutCandidate candidate;
  candidate.cut_x = clamp_virtual_box_center(cut_x);
  candidate.cut_y = clamp_virtual_box_center(cut_y);
  candidate.cut_z = clamp_virtual_box_center(cut_z);
  candidate.gap_score = gap_score;
  candidate.virtual_seed_source = seed_source;
  candidate.broken_edges = 0;
  candidate.edge_stretch = 0.0;
  candidate.shift_norm2 = 0.0;
  const PbctopoCutCandidate *active_incumbent =
      (incumbent != nullptr && incumbent->valid && !incumbent->ultra_clash)
          ? incumbent
          : nullptr;
  const PbctopoCutCandidate *metric_incumbent =
      (active_incumbent != nullptr && active_incumbent->broken_edges == 0)
          ? active_incumbent
          : nullptr;

  if (atoms.empty())
    return candidate;

  placed_points.resize(atoms.size());
  PointStats stats;
  vector<Int3> atom_image_offsets(atoms.size());

  const int global_shift_x = static_cast<int>(floor(candidate.cut_x)) - 1;
  const int global_shift_y = static_cast<int>(floor(candidate.cut_y)) - 1;
  const int global_shift_z = static_cast<int>(floor(candidate.cut_z)) - 1;

  for (size_t idx = 0; idx < atoms.size(); ++idx) {
    const auto &atom = atoms[idx];
    int search_shift_x = -1;
    int search_shift_y = -1;
    int search_shift_z = -1;
    if (!find_virtual_box_shift(atom.sx, candidate.cut_x, search_shift_x) ||
        !find_virtual_box_shift(atom.sy, candidate.cut_y, search_shift_y) ||
        !find_virtual_box_shift(atom.sz, candidate.cut_z, search_shift_z)) {
      return PbctopoCutCandidate{};
    }

    const int actual_shift_x = search_shift_x - 1 - global_shift_x;
    const int actual_shift_y = search_shift_y - 1 - global_shift_y;
    const int actual_shift_z = search_shift_z - 1 - global_shift_z;
    const Vec3 &translation =
        translation_cache.by_shift[actual_shift_x + 2][actual_shift_y + 2]
                                  [actual_shift_z + 2];
    Vec3 point{atom.x + translation.x, atom.y + translation.y,
               atom.z + translation.z};
    placed_points[idx] = point;
    stats.add(point.x, point.y, point.z);
    atom_image_offsets[idx] = {actual_shift_x, actual_shift_y, actual_shift_z};
  }

  if (stats.count != atoms.size())
    return PbctopoCutCandidate{};

  const vector<Int3> *metric_offsets = &atom_image_offsets;
  vector<Int3> component_image_offsets;
  vector<unsigned char> component_offset_seen;
  bool component_assignment_consistent =
      atom_component_index != nullptr &&
      atom_component_index->size() == atoms.size();
  if (component_assignment_consistent) {
    size_t component_count = 0;
    for (size_t comp_idx : *atom_component_index) {
      if (comp_idx >= atoms.size()) {
        component_assignment_consistent = false;
        break;
      }
      component_count = max(component_count, comp_idx + 1);
    }
    if (component_assignment_consistent && component_count > 0) {
      component_image_offsets.resize(component_count);
      component_offset_seen.assign(component_count, 0);
      for (size_t atom_idx = 0; atom_idx < atom_image_offsets.size(); ++atom_idx) {
        const size_t comp_idx = (*atom_component_index)[atom_idx];
        if (component_offset_seen[comp_idx] == 0) {
          component_image_offsets[comp_idx] = atom_image_offsets[atom_idx];
          component_offset_seen[comp_idx] = 1;
        } else if (!same_delta(component_image_offsets[comp_idx],
                               atom_image_offsets[atom_idx])) {
          component_assignment_consistent = false;
          break;
        }
      }
      if (component_assignment_consistent &&
          find(component_offset_seen.begin(), component_offset_seen.end(), 0) ==
              component_offset_seen.end()) {
        metric_offsets = &component_image_offsets;
      }
    }
  }
  candidate.shift_norm2 = pbctopo_component_shift_metric(Box, *metric_offsets);

  if (contacting_owner_pairs != nullptr) {
    PbctopoOwnerPairSplitDetail split_detail;
    if (candidate_has_split_contacting_owners(
            placed_points, atoms, Box, *contacting_owner_pairs,
            candidate.split_owner_pairs, candidate.max_owner_raw_minus_minimg_delta,
            &split_detail)) {
      if (reject_context != nullptr && split_detail.valid &&
          isfinite(frame_time_ps)) {
        append_pbctopo_log_line(
            "[pbctopo] reject split_contacting_owner_pair at frame " +
            to_string(static_cast<long long>(frame_time_ps)) + " via " +
            reject_context + ": owners=" + to_string(split_detail.owner_a) + "-" +
            to_string(split_detail.owner_b) +
            ", raw_delta_scaled=" + format_vec3(split_detail.raw_delta_scaled) +
            ", min_image_delta_scaled=" +
            format_vec3(split_detail.minimg_delta_scaled));
      }
      return candidate;
    }
  }

  if (continuity_pairs != nullptr) {
    double continuity_path = numeric_limits<double>::infinity();
    double continuity_max_d2 = numeric_limits<double>::infinity();
    if (!continuity_metrics_ok(placed_points, *continuity_pairs, continuity_cutoff,
                               continuity_path, continuity_max_d2,
                               metric_incumbent)) {
      return PbctopoCutCandidate{};
    }
    candidate.continuity_path = continuity_path;
    candidate.continuity_max_d2 = continuity_max_d2;
  }
  if (anchor_indices != nullptr) {
    double anchor_mst2 = numeric_limits<double>::infinity();
    double anchor_max_d2 = numeric_limits<double>::infinity();
    PbctopoAnchorGraphScore anchor_graph_score;
    const vector<PbctopoContinuityPair> empty_continuity_pairs;
    const vector<PbctopoContinuityPair> &anchor_continuity_pairs =
        (continuity_pairs != nullptr) ? *continuity_pairs
                                      : empty_continuity_pairs;
    if (!anchor_graph_metrics_enhanced(
            placed_points, *anchor_indices, anchor_continuity_pairs,
            anchor_mst2, anchor_max_d2, anchor_graph_score, metric_incumbent,
            PbctopoRefineScorePhase::Exact, PbctopoAnchorBackend::Prim,
            nullptr, anchor_context)) {
      return PbctopoCutCandidate{};
    }
    candidate.anchor_mst2 = anchor_mst2;
    candidate.anchor_max_d2 = anchor_max_d2;
    candidate.anchor_graph = anchor_graph_score;
  }

  candidate.bbox_volume = compute_bbox_volume(stats);
  candidate.rg2 = compute_rg2_from_stats(stats);
  const array<double, 3> principal_rg2 = compute_principal_rg2(stats);
  candidate.principal_rg2_max = principal_rg2[0];
  candidate.principal_rg2_mid = principal_rg2[1];
  candidate.principal_rg2_min = principal_rg2[2];
  const auto shape_metrics = compute_gyration_shape_metrics(principal_rg2);
  candidate.shape_kappa2 = shape_metrics.kappa2;
  candidate.shape_asphericity = shape_metrics.asphericity;
  candidate.shape_acylindricity = shape_metrics.acylindricity;

  for (const auto &edge : local_edges) {
    const Vec3 delta = placed_points[edge.a] - placed_points[edge.b];
    const double d2 = norm2(delta);
    if (d2 > edge.break_cutoff_d2)
      ++candidate.broken_edges;
    if (d2 > edge.ref_d2)
      candidate.edge_stretch += (d2 - edge.ref_d2) * edge.inv_ref_d2;
  }

  candidate.valid = true;
  if (!skip_clash &&
      preclash_cut_candidate_cannot_beat_incumbent(candidate,
                                                   active_incumbent)) {
    return PbctopoCutCandidate{};
  }
  if (!skip_clash)
    analyze_cut_candidate_ultraclash(placed_points, edge_set, candidate,
                                     clash_cutoff, &atoms,
                                     atom_component_index);
  return candidate;
}

[[maybe_unused]] vector<PbctopoCutOption>
collect_cut_options(const vector<PbctopoAtom> &atoms, char axis,
                    size_t max_options) {
  vector<double> coords;
  coords.reserve(atoms.size() * 2 + 2);
  coords.push_back(0.5);
  coords.push_back(2.5);
  for (const auto &atom : atoms) {
    const double coord =
        (axis == 'x') ? atom.sx : ((axis == 'y') ? atom.sy : atom.sz);
    coords.push_back(coord + 0.5);
    coords.push_back(coord + 1.5);
  }
  if (coords.empty())
    return {{1.5, 2.0}};

  sort(coords.begin(), coords.end());
  vector<PbctopoCutOption> options;
  options.reserve(coords.size());
  for (size_t i = 0; i + 1 < coords.size(); ++i) {
    const double current = coords[i];
    const double next = coords[i + 1];
    const double gap = next - current;
    if (gap <= 1.0e-9)
      continue;
    const double center = 0.5 * (current + next);
    if (center + 1.0e-9 < 0.5 || center - 1.0e-9 > 2.5)
      continue;
    options.push_back({center, gap});
  }

  sort(options.begin(), options.end(),
       [](const PbctopoCutOption &lhs, const PbctopoCutOption &rhs) {
          const int metric_cmp =
              compare_pbctopo_double_metric(lhs.gap, rhs.gap, 1.0e-9);
          if (metric_cmp != 0)
            return metric_cmp > 0;
         return lhs.cut < rhs.cut;
       });

  vector<PbctopoCutOption> best;
  best.reserve(min(max_options, options.size() + 1));
  auto push_unique = [&](double cut, double gap) {
    for (const auto &item : best) {
      if (fabs(item.cut - cut) <= 1.0e-6)
        return;
    }
    best.push_back({cut, gap});
  };

  for (const auto &option : options) {
    push_unique(option.cut, option.gap);
    if (best.size() >= max_options)
      break;
  }
  push_unique(1.5, 0.0);
  return best;
}

[[maybe_unused]] void
build_cut_positions(const vector<PbctopoAtom> &atoms, box Box, double cut_x,
                    double cut_y, double cut_z, vector<Vec3> &placed_points,
                    PointStats &stats, double &shift_norm2) {
  auto select_supercell_coord = [](double wrapped_coord, double center,
                                   double &selected_coord, int &selected_shift) {
    const double lo = center - 0.5;
    const double hi = center + 0.5;
    constexpr double kSelectEps = 1.0e-9;
    for (int copy_idx = 0; copy_idx < 3; ++copy_idx) {
      const double candidate = wrapped_coord + static_cast<double>(copy_idx);
      if (candidate + kSelectEps < lo)
        continue;
      if (candidate >= hi - kSelectEps)
        continue;
      selected_coord = candidate;
      selected_shift = copy_idx - 1;
      return;
    }

    const int copy_idx =
        max(0, min(2, static_cast<int>(ceil(lo - wrapped_coord - kSelectEps))));
    selected_coord = wrapped_coord + static_cast<double>(copy_idx);
    selected_shift = copy_idx - 1;
  };

  placed_points.resize(atoms.size());
  stats = PointStats{};
  shift_norm2 = 0.0;
  for (size_t idx = 0; idx < atoms.size(); ++idx) {
    const auto &atom = atoms[idx];
    double sx = 0.0;
    double sy = 0.0;
    double sz = 0.0;
    int shift_x = 0;
    int shift_y = 0;
    int shift_z = 0;
    select_supercell_coord(atom.sx, cut_x, sx, shift_x);
    select_supercell_coord(atom.sy, cut_y, sy, shift_y);
    select_supercell_coord(atom.sz, cut_z, sz, shift_z);
    Vec3 point{};
    scaled_to_cartesian(Box, sx, sy, sz, point.x, point.y, point.z);
    placed_points[idx] = point;
    stats.add(point.x, point.y, point.z);

    shift_norm2 += static_cast<double>(shift_x * shift_x + shift_y * shift_y +
                                       shift_z * shift_z);
  }
}

PbctopoCutCandidate score_cut_candidate(const vector<PbctopoLocalEdge> &local_edges,
                                        const vector<Vec3> &placed_points,
                                        double cut_x, double cut_y, double cut_z,
                                        double gap_score, const PointStats &stats,
                                        double shift_norm2,
                                        const PbctopoStage1EdgeContext *edge_context,
                                        const PbctopoEdgeScanSoA *edge_scan_soa,
                                        const PbctopoCutCandidate *edge_prune_threshold,
                                        const PbctopoCutCandidate *pre_edge_seed,
                                        PbctopoRefineScorePhase phase) {
  PbctopoCutCandidate candidate =
      (pre_edge_seed != nullptr) ? *pre_edge_seed : PbctopoCutCandidate{};
  candidate.cut_x = cut_x;
  candidate.cut_y = cut_y;
  candidate.cut_z = cut_z;
  candidate.gap_score = gap_score;
  candidate.broken_edges = 0;
  candidate.edge_stretch = 0.0;
  candidate.valid = false;

  if (stats.count == 0 || placed_points.empty())
    return candidate;

  candidate.bbox_volume = compute_bbox_volume(stats);
  candidate.rg2 = compute_rg2_from_stats(stats);
  const array<double, 3> principal_rg2 = compute_principal_rg2(stats);
  candidate.principal_rg2_max = principal_rg2[0];
  candidate.principal_rg2_mid = principal_rg2[1];
  candidate.principal_rg2_min = principal_rg2[2];
  const auto shape_metrics = compute_gyration_shape_metrics(principal_rg2);
  candidate.shape_kappa2 = shape_metrics.kappa2;
  candidate.shape_asphericity = shape_metrics.asphericity;
  candidate.shape_acylindricity = shape_metrics.acylindricity;
  candidate.shift_norm2 = shift_norm2;

  candidate.broken_edges = (edge_context != nullptr)
                               ? edge_context->fixed_broken_edges
                               : 0;
  candidate.edge_stretch = (edge_context != nullptr)
                               ? edge_context->fixed_edge_stretch
                               : 0.0;

  if (!accumulate_edge_scan(candidate, local_edges, placed_points, edge_context,
                            edge_scan_soa, edge_prune_threshold, phase)) {
    return {};
  }

  candidate.valid = true;
  return candidate;
}

PbctopoCutCandidate score_cut_candidate_principal_only(
    const vector<PbctopoLocalEdge> &local_edges,
    const vector<Vec3> &placed_points, double cut_x, double cut_y, double cut_z,
    double gap_score, const PointStats &stats, double shift_norm2,
    PbctopoCutScoreBreakdown *breakdown,
    const PbctopoStage1EdgeContext *edge_context,
    const PbctopoEdgeScanSoA *edge_scan_soa,
    const PbctopoCutCandidate *edge_prune_threshold,
    const PbctopoCutCandidate *pre_edge_seed) {
  PbctopoCutCandidate candidate =
      (pre_edge_seed != nullptr) ? *pre_edge_seed : PbctopoCutCandidate{};
  candidate.cut_x = cut_x;
  candidate.cut_y = cut_y;
  candidate.cut_z = cut_z;
  candidate.gap_score = gap_score;
  candidate.broken_edges = 0;
  candidate.edge_stretch = 0.0;
  candidate.valid = false;

  if (stats.count == 0 || placed_points.empty())
    return candidate;

  chrono::steady_clock::time_point section_started;
  if (breakdown != nullptr)
    section_started = chrono::steady_clock::now();
  const array<double, 3> principal_rg2 = compute_principal_rg2(stats);
  if (breakdown != nullptr) {
    breakdown->principal_s +=
        chrono::duration<double>(chrono::steady_clock::now() - section_started)
            .count();
    section_started = chrono::steady_clock::now();
  }
  candidate.principal_rg2_max = principal_rg2[0];
  candidate.principal_rg2_mid = principal_rg2[1];
  candidate.principal_rg2_min = principal_rg2[2];
  const auto shape_metrics = compute_gyration_shape_metrics(principal_rg2);
  candidate.shape_kappa2 = shape_metrics.kappa2;
  candidate.shape_asphericity = shape_metrics.asphericity;
  candidate.shape_acylindricity = shape_metrics.acylindricity;
  candidate.shift_norm2 = shift_norm2;

  candidate.broken_edges = (edge_context != nullptr)
                               ? edge_context->fixed_broken_edges
                               : 0;
  candidate.edge_stretch = (edge_context != nullptr)
                               ? edge_context->fixed_edge_stretch
                               : 0.0;

  if (!accumulate_edge_scan(candidate, local_edges, placed_points, edge_context,
                            edge_scan_soa, edge_prune_threshold,
                            PbctopoRefineScorePhase::CoarseStage1)) {
    if (breakdown != nullptr) {
      breakdown->edge_scan_s +=
          chrono::duration<double>(chrono::steady_clock::now() -
                                   section_started)
              .count();
    }
    return {};
  }
  if (breakdown != nullptr) {
    breakdown->edge_scan_s +=
        chrono::duration<double>(chrono::steady_clock::now() - section_started)
            .count();
  }

  candidate.valid = true;
  return candidate;
}


void analyze_cut_candidate_ultraclash(const vector<Vec3> &placed_points,
                                       const unordered_set<uint64_t> &edge_set,
                                       PbctopoCutCandidate &candidate,
                                       double cutoff,
                                       const vector<PbctopoAtom> *atoms,
                                       const vector<size_t> *atom_component_index) {
  thread_local PbctopoUltraclashScratch scratch;
  candidate.min_pair_d2 = numeric_limits<double>::infinity();
  candidate.clash_atom_a = numeric_limits<size_t>::max();
  candidate.clash_atom_b = numeric_limits<size_t>::max();
  candidate.ultra_clash = false;
  if (!candidate.valid || placed_points.size() < 2)
    return;

  const bool filter_heavy =
      atoms != nullptr && atoms->size() == placed_points.size();
  const double grid_cutoff =
      filter_heavy ? pbctopo_max_clash_cutoff(*atoms, cutoff) : cutoff;
  const double inv_grid = 1.0 / grid_cutoff;
  const double cutoff2 = cutoff * cutoff;
  const bool filter_same_component =
      atom_component_index != nullptr &&
      atom_component_index->size() == placed_points.size() &&
      candidate.broken_edges == 0;
  auto accumulate_pair = [&](size_t atom_idx, size_t other_idx) {
    if (filter_same_component &&
        (*atom_component_index)[atom_idx] ==
            (*atom_component_index)[other_idx])
      return false;
    if (edge_set.find(atom_pair_key(atom_idx, other_idx)) != edge_set.end())
      return false;
    const auto &point = placed_points[atom_idx];
    const auto &other = placed_points[other_idx];
    const double dx_atom = point.x - other.x;
    const double dy_atom = point.y - other.y;
    const double dz_atom = point.z - other.z;
    const double d2 = dx_atom * dx_atom + dy_atom * dy_atom + dz_atom * dz_atom;
    if (d2 < candidate.min_pair_d2) {
      candidate.min_pair_d2 = d2;
      candidate.clash_atom_a = atom_idx;
      candidate.clash_atom_b = other_idx;
    }
    const bool clashes =
        filter_heavy
            ? pbctopo_atoms_clash(d2, (*atoms)[atom_idx],
                                   (*atoms)[other_idx], cutoff)
            : d2 + kDistanceEps < cutoff2;
    if (clashes) {
      candidate.ultra_clash = true;
      return true;
    }
    return false;
  };

  auto analyze_with_hash_grid = [&]() {
    unordered_map<uint64_t, vector<size_t>> clash_grid;
    clash_grid.reserve(max<size_t>(64, placed_points.size()));
    for (size_t atom_idx = 0; atom_idx < placed_points.size(); ++atom_idx) {
      if (filter_heavy && !(*atoms)[atom_idx].is_heavy)
        continue;
      const auto &point = placed_points[atom_idx];
      const int cx = static_cast<int>(floor(point.x * inv_grid));
      const int cy = static_cast<int>(floor(point.y * inv_grid));
      const int cz = static_cast<int>(floor(point.z * inv_grid));

      for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dz = -1; dz <= 1; ++dz) {
            const auto it =
                clash_grid.find(pack_cell_key(cx + dx, cy + dy, cz + dz));
            if (it == clash_grid.end())
              continue;
            for (size_t other_idx : it->second) {
              if (accumulate_pair(atom_idx, other_idx))
                return;
            }
          }
        }
      }

      clash_grid[pack_cell_key(cx, cy, cz)].push_back(atom_idx);
    }
  };

  scratch.cell_x.resize(placed_points.size());
  scratch.cell_y.resize(placed_points.size());
  scratch.cell_z.resize(placed_points.size());
  int min_cx = numeric_limits<int>::max();
  int min_cy = numeric_limits<int>::max();
  int min_cz = numeric_limits<int>::max();
  int max_cx = numeric_limits<int>::min();
  int max_cy = numeric_limits<int>::min();
  int max_cz = numeric_limits<int>::min();

  for (size_t atom_idx = 0; atom_idx < placed_points.size(); ++atom_idx) {
    if (filter_heavy && !(*atoms)[atom_idx].is_heavy)
      continue;
    const auto &point = placed_points[atom_idx];
    const int cx = static_cast<int>(floor(point.x * inv_grid));
    const int cy = static_cast<int>(floor(point.y * inv_grid));
    const int cz = static_cast<int>(floor(point.z * inv_grid));
    scratch.cell_x[atom_idx] = cx;
    scratch.cell_y[atom_idx] = cy;
    scratch.cell_z[atom_idx] = cz;
    min_cx = min(min_cx, cx);
    min_cy = min(min_cy, cy);
    min_cz = min(min_cz, cz);
    max_cx = max(max_cx, cx);
    max_cy = max(max_cy, cy);
    max_cz = max(max_cz, cz);
  }
  if (min_cx == numeric_limits<int>::max())
    return;

  const size_t nx = static_cast<size_t>(max_cx - min_cx) + 1;
  const size_t ny = static_cast<size_t>(max_cy - min_cy) + 1;
  const size_t nz = static_cast<size_t>(max_cz - min_cz) + 1;
  const size_t max_cells = kPbctopoDenseClashGridMaxCells;
  if (nx == 0 || ny == 0 || nz == 0 || nx > max_cells || ny > max_cells ||
      nz > max_cells || ny > (max_cells / nx) ||
      nz > (max_cells / (nx * ny))) {
    analyze_with_hash_grid();
    return;
  }

  const size_t total_cells = nx * ny * nz;
  scratch.head.resize(total_cells);
  scratch.head_epoch.resize(total_cells, 0u);
  ++scratch.current_epoch;
  if (scratch.current_epoch == 0u) {
    scratch.current_epoch = 1u;
    fill(scratch.head_epoch.begin(), scratch.head_epoch.end(), 0u);
  }
  scratch.next.resize(placed_points.size());
  const size_t plane_stride = nx * ny;
  auto cell_index = [&](int cx, int cy, int cz) {
    return (static_cast<size_t>(cz - min_cz) * plane_stride) +
           (static_cast<size_t>(cy - min_cy) * nx) +
           static_cast<size_t>(cx - min_cx);
  };

  for (size_t atom_idx = 0; atom_idx < placed_points.size(); ++atom_idx) {
    if (filter_heavy && !(*atoms)[atom_idx].is_heavy)
      continue;
    const int cx = scratch.cell_x[atom_idx];
    const int cy = scratch.cell_y[atom_idx];
    const int cz = scratch.cell_z[atom_idx];
    for (int nz_cell = max(min_cz, cz - 1); nz_cell <= min(max_cz, cz + 1);
         ++nz_cell) {
      for (int ny_cell = max(min_cy, cy - 1); ny_cell <= min(max_cy, cy + 1);
           ++ny_cell) {
        for (int nx_cell = max(min_cx, cx - 1); nx_cell <= min(max_cx, cx + 1);
             ++nx_cell) {
          const size_t neighbor_cell_idx = cell_index(nx_cell, ny_cell, nz_cell);
          int other_atom =
              (scratch.head_epoch[neighbor_cell_idx] == scratch.current_epoch)
                  ? scratch.head[neighbor_cell_idx]
                  : -1;
          while (other_atom >= 0) {
            if (accumulate_pair(atom_idx, static_cast<size_t>(other_atom)))
              return;
            other_atom = scratch.next[static_cast<size_t>(other_atom)];
          }
        }
      }
    }

    const size_t own_cell_idx = cell_index(cx, cy, cz);
    scratch.next[atom_idx] =
        (scratch.head_epoch[own_cell_idx] == scratch.current_epoch)
            ? scratch.head[own_cell_idx]
            : -1;
    scratch.head[own_cell_idx] = static_cast<int>(atom_idx);
    scratch.head_epoch[own_cell_idx] = scratch.current_epoch;
  }
}

void apply_cut_candidate(const vector<PbctopoAtom> &atoms,
                         const vector<Vec3> &placed_points,
                         Multi_MOL &mol1_system, Multi_MOL &mol2_system,
                         size_t frame_index) {
  if (atoms.empty() || placed_points.size() != atoms.size())
    return;

  box Box = mol1_system.m_mol[frame_index].Box;
  const PbctopoLatticeMetric lattice(Box);
  double min_sx = numeric_limits<double>::infinity();
  double min_sy = numeric_limits<double>::infinity();
  double min_sz = numeric_limits<double>::infinity();
  double max_sx = -numeric_limits<double>::infinity();
  double max_sy = -numeric_limits<double>::infinity();
  double max_sz = -numeric_limits<double>::infinity();

  for (const auto &point : placed_points) {
    const auto scaled = lattice.cartesian_to_scaled(point.x, point.y, point.z);
    const double sx = scaled[0];
    const double sy = scaled[1];
    const double sz = scaled[2];
    min_sx = min(min_sx, sx);
    min_sy = min(min_sy, sy);
    min_sz = min(min_sz, sz);
    max_sx = max(max_sx, sx);
    max_sy = max(max_sy, sy);
    max_sz = max(max_sz, sz);
  }

  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);
  const Vec3 recenter_shift{
      (0.5 - 0.5 * (min_sx + max_sx)) * v1[0] +
          (0.5 - 0.5 * (min_sy + max_sy)) * v2[0] +
          (0.5 - 0.5 * (min_sz + max_sz)) * v3[0],
      (0.5 - 0.5 * (min_sx + max_sx)) * v1[1] +
          (0.5 - 0.5 * (min_sy + max_sy)) * v2[1] +
          (0.5 - 0.5 * (min_sz + max_sz)) * v3[1],
      (0.5 - 0.5 * (min_sx + max_sx)) * v1[2] +
          (0.5 - 0.5 * (min_sy + max_sy)) * v2[2] +
          (0.5 - 0.5 * (min_sz + max_sz)) * v3[2]};

  for (size_t idx = 0; idx < atoms.size(); ++idx) {
    const auto &atom = atoms[idx];
    const Vec3 point{placed_points[idx].x + recenter_shift.x,
                     placed_points[idx].y + recenter_shift.y,
                     placed_points[idx].z + recenter_shift.z};
    if (atom.owner == 0) {
      mol1_system.m_mol[frame_index].pos.x[atom.local_atom_index] = point.x;
      mol1_system.m_mol[frame_index].pos.y[atom.local_atom_index] = point.y;
      mol1_system.m_mol[frame_index].pos.z[atom.local_atom_index] = point.z;
    } else if (frame_index < mol2_system.m_mol.size()) {
      mol2_system.m_mol[frame_index].pos.x[atom.local_atom_index] = point.x;
      mol2_system.m_mol[frame_index].pos.y[atom.local_atom_index] = point.y;
      mol2_system.m_mol[frame_index].pos.z[atom.local_atom_index] = point.z;
    }
  }
}

void apply_pbctopo_points_direct(const vector<PbctopoAtom> &atoms,
                                 const vector<Vec3> &placed_points,
                                 Multi_MOL &mol1_system,
                                 Multi_MOL &mol2_system,
                                 size_t frame_index) {
  if (atoms.empty() || placed_points.size() != atoms.size())
    return;

  for (size_t idx = 0; idx < atoms.size(); ++idx) {
    const auto &atom = atoms[idx];
    const Vec3 &point = placed_points[idx];
    if (atom.owner == 0) {
      mol1_system.m_mol[frame_index].pos.x[atom.local_atom_index] = point.x;
      mol1_system.m_mol[frame_index].pos.y[atom.local_atom_index] = point.y;
      mol1_system.m_mol[frame_index].pos.z[atom.local_atom_index] = point.z;
    } else if (frame_index < mol2_system.m_mol.size()) {
      mol2_system.m_mol[frame_index].pos.x[atom.local_atom_index] = point.x;
      mol2_system.m_mol[frame_index].pos.y[atom.local_atom_index] = point.y;
      mol2_system.m_mol[frame_index].pos.z[atom.local_atom_index] = point.z;
    }
  }
}

[[maybe_unused]] void
append_fragments_for_system(const Multi_MOL &system, const MOL &mol, int owner,
                            vector<PbctopoFragment> &fragments,
                            vector<PbctopoAtom> &atoms) {
  const auto ranges = build_ranges_for_system(system, mol);
  for (const auto &range : ranges) {
    PbctopoFragment fragment;
    fragment.id = fragments.size();
    fragment.owner = owner;
    fragment.source_index = range.source_index;
    fragment.begin_local = range.begin;
    fragment.end_local = range.end;
    fragment.base_atom_begin = atoms.size();
    fragment.atom_count = range.end - range.begin;
    fragment.topo_backed = range.topo_backed;

    for (size_t local = range.begin; local < range.end; ++local) {
      double r[3]{mol.pos.x[local], mol.pos.y[local], mol.pos.z[local]};
      double s[3]{};
      mol.Box.to_scaled(r, s);
      for (double &coord : s)
        coord -= floor(coord);
      scaled_to_cartesian(mol.Box, s[0], s[1], s[2], r[0], r[1], r[2]);
      PbctopoAtom atom;
      atom.global_atom_id = atoms.size();
      atom.fragment_id = fragment.id;
      atom.owner = owner;
      atom.local_atom_index = local;
      atom.x = r[0];
      atom.y = r[1];
      atom.z = r[2];
      atom.sx = s[0];
      atom.sy = s[1];
      atom.sz = s[2];
      atoms.push_back(atom);
    }
    fragment.base_atom_end = atoms.size();
    fragments.push_back(fragment);
  }
}

vector<PbctopoShift> build_topology_shifts(box Box, int shell_radius = 1) {
  return titan_pbctopo::build_component_image_shifts(Box, shell_radius);
}

[[maybe_unused]] void
build_supercell_for_frame(const vector<PbctopoFragment> &fragments,
                          const vector<PbctopoAtom> &base_atoms,
                          const vector<PbctopoShift> &shifts,
                          vector<PbctopoFragmentCopyNode> &nodes,
                          vector<PbctopoAtomCopy> &atom_copies) {
  nodes.clear();
  atom_copies.clear();
  nodes.reserve(fragments.size() * shifts.size());
  atom_copies.reserve(base_atoms.size() * shifts.size());

  for (const auto &fragment : fragments) {
    for (const auto &shift : shifts) {
      PbctopoFragmentCopyNode node;
      node.node_id = nodes.size();
      node.fragment_id = fragment.id;
      node.shift_index = shift.shift_index;
      node.shift_ix = shift.ix;
      node.shift_iy = shift.iy;
      node.shift_iz = shift.iz;
      node.tx = shift.tx;
      node.ty = shift.ty;
      node.tz = shift.tz;
      node.atom_begin = atom_copies.size();
      node.min_x = numeric_limits<double>::infinity();
      node.min_y = numeric_limits<double>::infinity();
      node.min_z = numeric_limits<double>::infinity();
      node.max_x = -numeric_limits<double>::infinity();
      node.max_y = -numeric_limits<double>::infinity();
      node.max_z = -numeric_limits<double>::infinity();

      for (size_t atom_idx = fragment.base_atom_begin;
           atom_idx < fragment.base_atom_end; ++atom_idx) {
        const PbctopoAtom &base_atom = base_atoms[atom_idx];
        PbctopoAtomCopy atom_copy;
        atom_copy.global_atom_id = base_atom.global_atom_id;
        atom_copy.fragment_id = fragment.id;
        atom_copy.node_id = node.node_id;
        atom_copy.owner = base_atom.owner;
        atom_copy.local_atom_index = base_atom.local_atom_index;
        atom_copy.x = base_atom.x + shift.tx;
        atom_copy.y = base_atom.y + shift.ty;
        atom_copy.z = base_atom.z + shift.tz;
        node.min_x = min(node.min_x, atom_copy.x);
        node.min_y = min(node.min_y, atom_copy.y);
        node.min_z = min(node.min_z, atom_copy.z);
        node.max_x = max(node.max_x, atom_copy.x);
        node.max_y = max(node.max_y, atom_copy.y);
        node.max_z = max(node.max_z, atom_copy.z);
        atom_copies.push_back(atom_copy);
      }

      node.atom_end = atom_copies.size();
      nodes.push_back(node);
    }
  }
}

bool are_topology_neighbors(const PbctopoFragment &lhs,
                            const PbctopoFragment &rhs) {
  if (!lhs.topo_backed || !rhs.topo_backed || lhs.owner != rhs.owner)
    return false;
  if (lhs.source_index + 1 == rhs.source_index ||
      rhs.source_index + 1 == lhs.source_index)
    return true;
  return false;
}

double axis_gap(double min_a, double max_a, double min_b, double max_b) {
  if (max_a < min_b)
    return min_b - max_a;
  if (max_b < min_a)
    return min_a - max_b;
  return 0.0;
}

double bbox_distance2(const PbctopoFragmentCopyNode &lhs,
                      const PbctopoFragmentCopyNode &rhs) {
  const double dx = axis_gap(lhs.min_x, lhs.max_x, rhs.min_x, rhs.max_x);
  const double dy = axis_gap(lhs.min_y, lhs.max_y, rhs.min_y, rhs.max_y);
  const double dz = axis_gap(lhs.min_z, lhs.max_z, rhs.min_z, rhs.max_z);
  return dx * dx + dy * dy + dz * dz;
}

[[maybe_unused]] vector<vector<size_t>>
build_fragment_copy_graph(const vector<PbctopoFragment> &fragments,
                          const vector<PbctopoFragmentCopyNode> &nodes,
                          double geometry_cutoff, double topology_cutoff) {
  vector<vector<size_t>> adjacency(nodes.size());
  const double geometry_cutoff2 = geometry_cutoff * geometry_cutoff;
  const double topology_cutoff2 = topology_cutoff * topology_cutoff;
  const double max_cutoff = max(geometry_cutoff, topology_cutoff);
  const double inv_cell = 1.0 / max_cutoff;
  unordered_map<ClashCellKey, vector<size_t>, ClashCellKeyHash> cell_nodes;
  cell_nodes.reserve(max<size_t>(64, nodes.size() * 2));

  for (size_t node_id = 0; node_id < nodes.size(); ++node_id) {
    const auto &node = nodes[node_id];
    const int min_cx = static_cast<int>(floor((node.min_x - max_cutoff) * inv_cell));
    const int min_cy = static_cast<int>(floor((node.min_y - max_cutoff) * inv_cell));
    const int min_cz = static_cast<int>(floor((node.min_z - max_cutoff) * inv_cell));
    const int max_cx = static_cast<int>(floor((node.max_x + max_cutoff) * inv_cell));
    const int max_cy = static_cast<int>(floor((node.max_y + max_cutoff) * inv_cell));
    const int max_cz = static_cast<int>(floor((node.max_z + max_cutoff) * inv_cell));

    for (int cx = min_cx; cx <= max_cx; ++cx) {
      for (int cy = min_cy; cy <= max_cy; ++cy) {
        for (int cz = min_cz; cz <= max_cz; ++cz)
          cell_nodes[ClashCellKey{cx, cy, cz}].push_back(node_id);
      }
    }
  }

  unordered_set<uint64_t> seen_pairs;
  seen_pairs.reserve(max<size_t>(128, nodes.size() * 8));
  for (const auto &cell_entry : cell_nodes) {
    const auto &cell_list = cell_entry.second;
    for (size_t a = 0; a < cell_list.size(); ++a) {
      for (size_t b = a + 1; b < cell_list.size(); ++b) {
        size_t lhs_id = cell_list[a];
        size_t rhs_id = cell_list[b];
        if (lhs_id == rhs_id)
          continue;
        if (lhs_id > rhs_id)
          swap(lhs_id, rhs_id);
        const uint64_t pair_key =
            (static_cast<uint64_t>(lhs_id) << 32) | static_cast<uint64_t>(rhs_id);
        if (!seen_pairs.insert(pair_key).second)
          continue;

        const auto &lhs = nodes[lhs_id];
        const auto &rhs = nodes[rhs_id];
        if (lhs.fragment_id == rhs.fragment_id)
          continue;

        const bool topo_related =
            are_topology_neighbors(fragments[lhs.fragment_id],
                                   fragments[rhs.fragment_id]);
        const double cutoff2 = topo_related ? topology_cutoff2 : geometry_cutoff2;
        if (bbox_distance2(lhs, rhs) > cutoff2)
          continue;

        adjacency[lhs_id].push_back(rhs_id);
        adjacency[rhs_id].push_back(lhs_id);
      }
    }
  }
  return adjacency;
}

double node_shift_norm2(const PbctopoFragmentCopyNode &node) {
  return static_cast<double>(node.shift_ix * node.shift_ix +
                             node.shift_iy * node.shift_iy +
                             node.shift_iz * node.shift_iz);
}

[[maybe_unused]] void
sort_fragment_copy_graph(vector<vector<size_t>> &adjacency,
                         const vector<PbctopoFragmentCopyNode> &nodes) {
  for (size_t node_id = 0; node_id < adjacency.size(); ++node_id) {
    auto &neighbors = adjacency[node_id];
    sort(neighbors.begin(), neighbors.end(),
         [&](size_t lhs, size_t rhs) {
           const double lhs_d2 = bbox_distance2(nodes[node_id], nodes[lhs]);
           const double rhs_d2 = bbox_distance2(nodes[node_id], nodes[rhs]);
           if (fabs(lhs_d2 - rhs_d2) > 1.0e-9)
             return lhs_d2 < rhs_d2;
           const double lhs_shift = node_shift_norm2(nodes[lhs]);
           const double rhs_shift = node_shift_norm2(nodes[rhs]);
           if (fabs(lhs_shift - rhs_shift) > 1.0e-9)
             return lhs_shift < rhs_shift;
           return lhs < rhs;
         });
  }
}

double min_bbox_distance2_to_selected(
    size_t node_id, const vector<size_t> &selected_nodes,
    const vector<PbctopoFragmentCopyNode> &nodes) {
  double best_d2 = numeric_limits<double>::infinity();
  for (size_t selected_id : selected_nodes)
    best_d2 = min(best_d2, bbox_distance2(nodes[node_id], nodes[selected_id]));
  return best_d2;
}

[[maybe_unused]] bool
build_candidate_block_from_root(size_t root_node,
                                const vector<PbctopoFragment> &fragments,
                                const vector<PbctopoFragmentCopyNode> &nodes,
                                const vector<vector<size_t>> &adjacency,
                                size_t total_atom_count,
                                vector<size_t> &component_nodes) {
  component_nodes.clear();
  if (root_node >= nodes.size())
    return false;

  vector<int> fragment_assignment(fragments.size(), -1);
  const size_t root_fragment = nodes[root_node].fragment_id;
  fragment_assignment[root_fragment] = static_cast<int>(root_node);
  component_nodes.push_back(root_node);

  size_t component_atom_count = fragments[root_fragment].atom_count;
  if (component_atom_count > total_atom_count)
    return false;

  while (component_nodes.size() < fragments.size()) {
    size_t next_fragment = numeric_limits<size_t>::max();
    size_t next_node = numeric_limits<size_t>::max();
    double next_d2 = numeric_limits<double>::infinity();
    double next_shift = numeric_limits<double>::infinity();

    for (size_t fragment_id = 0; fragment_id < fragments.size(); ++fragment_id) {
      if (fragment_assignment[fragment_id] != -1)
        continue;
      vector<size_t> options =
          collect_fragment_options(fragment_id, component_nodes, nodes, adjacency, 1);
      if (options.empty())
        continue;
      const size_t candidate_node = options.front();
      const double candidate_d2 =
          min_bbox_distance2_to_selected(candidate_node, component_nodes, nodes);
      const double candidate_shift = node_shift_norm2(nodes[candidate_node]);
      if (next_fragment == numeric_limits<size_t>::max() ||
          candidate_d2 < next_d2 - 1.0e-9 ||
          (fabs(candidate_d2 - next_d2) <= 1.0e-9 &&
           candidate_shift < next_shift - 1.0e-9) ||
          (fabs(candidate_d2 - next_d2) <= 1.0e-9 &&
           fabs(candidate_shift - next_shift) <= 1.0e-9 &&
           candidate_node < next_node)) {
        next_fragment = fragment_id;
        next_node = candidate_node;
        next_d2 = candidate_d2;
        next_shift = candidate_shift;
      }
    }

    if (next_fragment == numeric_limits<size_t>::max())
      return false;

    fragment_assignment[next_fragment] = static_cast<int>(next_node);
    component_nodes.push_back(next_node);
    component_atom_count += fragments[next_fragment].atom_count;
    if (component_atom_count > total_atom_count)
      return false;
  }

  if (component_atom_count != total_atom_count ||
      component_nodes.size() != fragments.size()) {
    return false;
  }

  sort(component_nodes.begin(), component_nodes.end());
  return true;
}

vector<size_t> collect_fragment_options(
    size_t fragment_id, const vector<size_t> &selected_nodes,
    const vector<PbctopoFragmentCopyNode> &nodes,
    const vector<vector<size_t>> &adjacency, size_t max_options) {
  unordered_map<size_t, double> candidate_dist2;
  candidate_dist2.reserve(16);
  for (size_t selected_node : selected_nodes) {
    for (size_t neighbor : adjacency[selected_node]) {
      if (nodes[neighbor].fragment_id != fragment_id)
        continue;
      const double d2 = bbox_distance2(nodes[selected_node], nodes[neighbor]);
      const auto it = candidate_dist2.find(neighbor);
      if (it == candidate_dist2.end() || d2 < it->second)
        candidate_dist2[neighbor] = d2;
    }
  }

  vector<pair<size_t, double>> ranked(candidate_dist2.begin(),
                                      candidate_dist2.end());
  sort(ranked.begin(), ranked.end(),
       [&](const pair<size_t, double> &lhs, const pair<size_t, double> &rhs) {
         int metric_cmp = compare_pbctopo_double_metric(
             lhs.second, rhs.second, 1.0e-9);
         if (metric_cmp != 0)
           return metric_cmp < 0;
         const double lhs_shift = node_shift_norm2(nodes[lhs.first]);
         const double rhs_shift = node_shift_norm2(nodes[rhs.first]);
         metric_cmp = compare_pbctopo_double_metric(lhs_shift, rhs_shift,
                                                     1.0e-9);
         if (metric_cmp != 0)
           return metric_cmp < 0;
         return lhs.first < rhs.first;
       });

  vector<size_t> options;
  options.reserve(min(max_options, ranked.size()));
  for (size_t i = 0; i < ranked.size() && i < max_options; ++i)
    options.push_back(ranked[i].first);
  return options;
}

[[maybe_unused]] void enumerate_connected_block_candidates_recursive(
    const vector<PbctopoFragment> &fragments,
    const vector<PbctopoFragmentCopyNode> &nodes,
    const vector<vector<size_t>> &adjacency, size_t max_options_per_fragment,
    vector<int> &assignment, vector<size_t> &selected_nodes,
    unordered_set<vector<size_t>, NodeVectorHash> &seen_components,
    vector<vector<size_t>> &out_components) {
  if (selected_nodes.size() == fragments.size()) {
    vector<size_t> candidate_nodes = selected_nodes;
    sort(candidate_nodes.begin(), candidate_nodes.end());
    if (seen_components.insert(candidate_nodes).second)
      out_components.push_back(std::move(candidate_nodes));
    return;
  }

  size_t next_fragment = numeric_limits<size_t>::max();
  vector<size_t> next_options;
  for (size_t fragment_id = 0; fragment_id < fragments.size(); ++fragment_id) {
    if (assignment[fragment_id] != -1)
      continue;
    vector<size_t> options =
        collect_fragment_options(fragment_id, selected_nodes, nodes, adjacency,
                                 max_options_per_fragment);
    if (options.empty())
      continue;
    if (next_fragment == numeric_limits<size_t>::max() ||
        options.size() < next_options.size()) {
      next_fragment = fragment_id;
      next_options = std::move(options);
    }
  }

  if (next_fragment == numeric_limits<size_t>::max())
    return;

  for (size_t node_id : next_options) {
    assignment[next_fragment] = static_cast<int>(node_id);
    selected_nodes.push_back(node_id);
    enumerate_connected_block_candidates_recursive(
        fragments, nodes, adjacency, max_options_per_fragment, assignment,
        selected_nodes, seen_components, out_components);
    selected_nodes.pop_back();
    assignment[next_fragment] = -1;
  }
}

[[maybe_unused]] vector<vector<size_t>> enumerate_connected_block_candidates(
    const vector<PbctopoFragment> &fragments,
    const vector<PbctopoFragmentCopyNode> &nodes,
    const vector<vector<size_t>> &adjacency, size_t max_options_per_fragment) {
  vector<vector<size_t>> out_components;
  if (fragments.empty())
    return out_components;

  unordered_set<vector<size_t>, NodeVectorHash> seen_components;
  const size_t root_fragment = 0;
  for (size_t node_id = 0; node_id < nodes.size(); ++node_id) {
    if (nodes[node_id].fragment_id != root_fragment)
      continue;
    vector<int> assignment(fragments.size(), -1);
    vector<size_t> selected_nodes{node_id};
    assignment[root_fragment] = static_cast<int>(node_id);
    enumerate_connected_block_candidates_recursive(
        fragments, nodes, adjacency, max_options_per_fragment, assignment,
        selected_nodes, seen_components, out_components);
  }

  return out_components;
}

Vec3 node_bbox_center(const PbctopoFragmentCopyNode &node) {
  return {(node.min_x + node.max_x) * 0.5, (node.min_y + node.max_y) * 0.5,
          (node.min_z + node.max_z) * 0.5};
}

double compute_mst2(const vector<Vec3> &points) {
  if (points.size() < 2)
    return 0.0;

  vector<double> best(points.size(), numeric_limits<double>::infinity());
  vector<char> used(points.size(), 0);
  best[0] = 0.0;
  double mst2 = 0.0;

  for (size_t step = 0; step < points.size(); ++step) {
    size_t next = numeric_limits<size_t>::max();
    double next_cost = numeric_limits<double>::infinity();
    for (size_t i = 0; i < points.size(); ++i) {
      if (used[i] != 0)
        continue;
      if (best[i] < next_cost) {
        next = i;
        next_cost = best[i];
      }
    }
    if (next == numeric_limits<size_t>::max())
      break;
    used[next] = 1;
    mst2 += next_cost;
    for (size_t i = 0; i < points.size(); ++i) {
      if (used[i] != 0)
        continue;
      const double d2 = norm2(points[next] - points[i]);
      if (d2 < best[i])
        best[i] = d2;
    }
  }

  return mst2;
}

size_t count_component_connectivity_edges(
    const vector<size_t> &component_nodes,
    const vector<vector<size_t>> &adjacency) {
  if (component_nodes.empty() || adjacency.empty())
    return 0;
  vector<char> node_mask(adjacency.size(), 0);
  for (size_t node_id : component_nodes) {
    if (node_id < node_mask.size())
      node_mask[node_id] = 1;
  }
  size_t edge_count = 0;
  for (size_t node_id : component_nodes) {
    if (node_id >= adjacency.size())
      continue;
    for (size_t neighbor : adjacency[node_id]) {
      if (neighbor >= node_mask.size() || node_mask[neighbor] == 0 ||
          node_id >= neighbor) {
        continue;
      }
      ++edge_count;
    }
  }
  return edge_count;
}

PbctopoCandidate
score_coarse_node_combination(const vector<size_t> &component_nodes,
                              const vector<PbctopoFragmentCopyNode> &nodes,
                              const vector<vector<size_t>> &adjacency) {
  PbctopoCandidate candidate;
  candidate.node_ids = component_nodes;
  if (component_nodes.empty())
    return candidate;

  vector<Vec3> centers;
  centers.reserve(component_nodes.size());
  double min_x = numeric_limits<double>::infinity();
  double min_y = numeric_limits<double>::infinity();
  double min_z = numeric_limits<double>::infinity();
  double max_x = -numeric_limits<double>::infinity();
  double max_y = -numeric_limits<double>::infinity();
  double max_z = -numeric_limits<double>::infinity();
  candidate.shift_norm2 = 0.0;

  for (size_t node_id : component_nodes) {
    const auto &node = nodes[node_id];
    min_x = min(min_x, node.min_x);
    min_y = min(min_y, node.min_y);
    min_z = min(min_z, node.min_z);
    max_x = max(max_x, node.max_x);
    max_y = max(max_y, node.max_y);
    max_z = max(max_z, node.max_z);
    centers.push_back(node_bbox_center(node));
    candidate.shift_norm2 += node_shift_norm2(node);
  }

  candidate.bbox_volume =
      max(0.0, max_x - min_x) * max(0.0, max_y - min_y) * max(0.0, max_z - min_z);
  candidate.rg2 = compute_rg2(centers);
  candidate.fragment_rg2 = candidate.rg2;
  const PointStats center_stats = build_point_stats(centers);
  const auto center_principal_rg2 = compute_principal_rg2(center_stats);
  const auto center_shape_metrics =
      compute_gyration_shape_metrics(center_principal_rg2);
  candidate.shape_kappa2 = center_shape_metrics.kappa2;
  candidate.shape_asphericity = center_shape_metrics.asphericity;
  candidate.shape_acylindricity = center_shape_metrics.acylindricity;
  candidate.fragment_kappa2 = center_shape_metrics.kappa2;
  candidate.fragment_asphericity = center_shape_metrics.asphericity;
  candidate.fragment_acylindricity = center_shape_metrics.acylindricity;
  candidate.fragment_mst2 = compute_mst2(centers);
  candidate.connectivity_edges =
      count_component_connectivity_edges(component_nodes, adjacency);
  candidate.valid = true;
  return candidate;
}

bool better_pbctopo_coarse_candidate(const PbctopoCandidate &lhs,
                                     const PbctopoCandidate &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  if (lhs.connectivity_edges != rhs.connectivity_edges)
    return lhs.connectivity_edges > rhs.connectivity_edges;
  int metric_cmp =
      compare_pbctopo_double_metric(lhs.fragment_mst2, rhs.fragment_mst2);
  if (metric_cmp != 0)
    return metric_cmp < 0;
  const int fragment_kappa_cmp =
      compare_shape_kappa(lhs.fragment_kappa2, rhs.fragment_kappa2);
  if (fragment_kappa_cmp != 0)
    return fragment_kappa_cmp < 0;
  const int fragment_asph_cmp = compare_shape_moment(
      lhs.fragment_asphericity, rhs.fragment_asphericity);
  if (fragment_asph_cmp != 0)
    return fragment_asph_cmp < 0;
  const int fragment_acyl_cmp = compare_shape_moment(
      lhs.fragment_acylindricity, rhs.fragment_acylindricity);
  if (fragment_acyl_cmp != 0)
    return fragment_acyl_cmp < 0;
  metric_cmp = compare_pbctopo_double_metric(lhs.bbox_volume, rhs.bbox_volume);
  if (metric_cmp != 0)
    return metric_cmp < 0;
  const int shape_kappa_cmp =
      compare_shape_kappa(lhs.shape_kappa2, rhs.shape_kappa2);
  if (shape_kappa_cmp != 0)
    return shape_kappa_cmp < 0;
  const int shape_asph_cmp =
      compare_shape_moment(lhs.shape_asphericity, rhs.shape_asphericity);
  if (shape_asph_cmp != 0)
    return shape_asph_cmp < 0;
  const int shape_acyl_cmp =
      compare_shape_moment(lhs.shape_acylindricity, rhs.shape_acylindricity);
  if (shape_acyl_cmp != 0)
    return shape_acyl_cmp < 0;
  metric_cmp = compare_pbctopo_double_metric(lhs.shift_norm2, rhs.shift_norm2);
  if (metric_cmp != 0)
    return metric_cmp < 0;
  return lhs.node_ids < rhs.node_ids;
}

[[maybe_unused]] void enumerate_small_fragment_shift_combinations_recursive(
    size_t fragment_id, const vector<PbctopoFragment> &fragments,
    const vector<PbctopoFragmentCopyNode> &nodes,
    const vector<vector<size_t>> &adjacency, vector<size_t> &selected_nodes,
    vector<PbctopoCandidate> &shortlist, size_t shortlist_limit) {
  if (fragment_id == fragments.size()) {
    PbctopoCandidate candidate =
        score_coarse_node_combination(selected_nodes, nodes, adjacency);
    shortlist.push_back(std::move(candidate));
    if (shortlist.size() > shortlist_limit * 2) {
      sort(shortlist.begin(), shortlist.end(), better_pbctopo_coarse_candidate);
      shortlist.resize(shortlist_limit);
    }
    return;
  }

  const size_t base_node = fragment_id * 27;
  for (size_t shift_idx = 0; shift_idx < 27; ++shift_idx) {
    selected_nodes.push_back(base_node + shift_idx);
    enumerate_small_fragment_shift_combinations_recursive(
        fragment_id + 1, fragments, nodes, adjacency, selected_nodes,
        shortlist, shortlist_limit);
    selected_nodes.pop_back();
  }
}

[[maybe_unused]] vector<vector<size_t>>
enumerate_small_fragment_shift_combinations(
    const vector<PbctopoFragment> &fragments,
    const vector<PbctopoFragmentCopyNode> &nodes,
    const vector<vector<size_t>> &adjacency, size_t shortlist_limit) {
  vector<PbctopoCandidate> shortlist;
  vector<size_t> selected_nodes;
  selected_nodes.reserve(fragments.size());
  enumerate_small_fragment_shift_combinations_recursive(
      0, fragments, nodes, adjacency, selected_nodes, shortlist,
      shortlist_limit);
  sort(shortlist.begin(), shortlist.end(), better_pbctopo_coarse_candidate);
  if (shortlist.size() > shortlist_limit)
    shortlist.resize(shortlist_limit);

  vector<vector<size_t>> components;
  components.reserve(shortlist.size());
  for (auto &candidate : shortlist)
    components.push_back(candidate.node_ids);
  return components;
}

[[maybe_unused]] PbctopoCandidate
score_component_candidate(const vector<size_t> &component_nodes,
                          const vector<PbctopoFragmentCopyNode> &nodes,
                          const vector<PbctopoAtomCopy> &atom_copies,
                          const vector<vector<size_t>> &adjacency) {
  PbctopoCandidate candidate;
  candidate.node_ids = component_nodes;
  if (component_nodes.empty())
    return candidate;

  vector<Vec3> centers;
  PointStats stats;
  candidate.shift_norm2 = 0.0;
  centers.reserve(component_nodes.size());
  for (size_t node_id : component_nodes) {
    const auto &node = nodes[node_id];
    candidate.shift_norm2 +=
        static_cast<double>(node.shift_ix * node.shift_ix + node.shift_iy * node.shift_iy +
                            node.shift_iz * node.shift_iz) *
        static_cast<double>(node.atom_end - node.atom_begin);
    centers.push_back(node_bbox_center(node));
    for (size_t atom_idx = node.atom_begin; atom_idx < node.atom_end; ++atom_idx) {
      const auto &atom = atom_copies[atom_idx];
      stats.add(atom.x, atom.y, atom.z);
    }
  }

  candidate.hull_valid = false;
  candidate.hull_volume = numeric_limits<double>::infinity();
  if (stats.count > 0) {
    candidate.bbox_volume = compute_bbox_volume(stats);
    candidate.rg2 = compute_rg2_from_stats(stats);
    const auto principal_rg2 = compute_principal_rg2(stats);
    const auto shape_metrics = compute_gyration_shape_metrics(principal_rg2);
    candidate.shape_kappa2 = shape_metrics.kappa2;
    candidate.shape_asphericity = shape_metrics.asphericity;
    candidate.shape_acylindricity = shape_metrics.acylindricity;
  }
  candidate.fragment_rg2 = compute_rg2(centers);
  const PointStats center_stats = build_point_stats(centers);
  const auto center_principal_rg2 = compute_principal_rg2(center_stats);
  const auto center_shape_metrics =
      compute_gyration_shape_metrics(center_principal_rg2);
  candidate.fragment_kappa2 = center_shape_metrics.kappa2;
  candidate.fragment_asphericity = center_shape_metrics.asphericity;
  candidate.fragment_acylindricity = center_shape_metrics.acylindricity;
  candidate.fragment_mst2 = compute_mst2(centers);
  candidate.connectivity_edges =
      count_component_connectivity_edges(component_nodes, adjacency);
  candidate.valid = stats.count > 0;
  return candidate;
}

[[maybe_unused]] bool
analyze_candidate_hard_clash(const vector<size_t> &component_nodes,
                             const vector<PbctopoFragmentCopyNode> &nodes,
                             const vector<PbctopoAtomCopy> &atom_copies,
                             double clash_cutoff,
                             double &min_interfragment_d2) {
  min_interfragment_d2 = numeric_limits<double>::infinity();
  if (component_nodes.empty())
    return false;

  size_t total_atoms = 0;
  for (size_t node_id : component_nodes) {
    const auto &node = nodes[node_id];
    total_atoms += node.atom_end - node.atom_begin;
  }

  unordered_map<ClashCellKey, vector<size_t>, ClashCellKeyHash> clash_grid;
  clash_grid.reserve(max<size_t>(64, total_atoms / 2));
  const double inv_grid = 1.0 / clash_cutoff;
  const double clash_cutoff2 = clash_cutoff * clash_cutoff;

  for (size_t node_id : component_nodes) {
    const auto &node = nodes[node_id];
    for (size_t atom_idx = node.atom_begin; atom_idx < node.atom_end; ++atom_idx) {
      const auto &atom = atom_copies[atom_idx];
      const int cx = static_cast<int>(floor(atom.x * inv_grid));
      const int cy = static_cast<int>(floor(atom.y * inv_grid));
      const int cz = static_cast<int>(floor(atom.z * inv_grid));

      for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dz = -1; dz <= 1; ++dz) {
            const ClashCellKey neighbor{cx + dx, cy + dy, cz + dz};
            const auto it = clash_grid.find(neighbor);
            if (it == clash_grid.end())
              continue;
            for (size_t other_idx : it->second) {
              const auto &other = atom_copies[other_idx];
              if (other.fragment_id == atom.fragment_id)
                continue;
              const double dx_atom = atom.x - other.x;
              const double dy_atom = atom.y - other.y;
              const double dz_atom = atom.z - other.z;
              const double d2 =
                  dx_atom * dx_atom + dy_atom * dy_atom + dz_atom * dz_atom;
              if (d2 < min_interfragment_d2)
                min_interfragment_d2 = d2;
              if (d2 + kDistanceEps < clash_cutoff2)
                return true;
            }
          }
        }
      }

      clash_grid[ClashCellKey{cx, cy, cz}].push_back(atom_idx);
    }
  }

  return false;
}

[[maybe_unused]] bool
better_pbctopo_candidate(const PbctopoCandidate &lhs,
                         const PbctopoCandidate &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  if (lhs.hard_clash != rhs.hard_clash)
    return rhs.hard_clash;
  if (lhs.connectivity_edges != rhs.connectivity_edges)
    return lhs.connectivity_edges > rhs.connectivity_edges;
  int metric_cmp = compare_pbctopo_double_metric(lhs.fragment_mst2,
                                                  rhs.fragment_mst2);
  if (metric_cmp != 0)
    return metric_cmp < 0;
  const int fragment_kappa_cmp =
      compare_shape_kappa(lhs.fragment_kappa2, rhs.fragment_kappa2);
  if (fragment_kappa_cmp != 0)
    return fragment_kappa_cmp < 0;
  const int fragment_asph_cmp = compare_shape_moment(
      lhs.fragment_asphericity, rhs.fragment_asphericity);
  if (fragment_asph_cmp != 0)
    return fragment_asph_cmp < 0;
  const int fragment_acyl_cmp = compare_shape_moment(
      lhs.fragment_acylindricity, rhs.fragment_acylindricity);
  if (fragment_acyl_cmp != 0)
    return fragment_acyl_cmp < 0;
  const int shape_kappa_cmp =
      compare_shape_kappa(lhs.shape_kappa2, rhs.shape_kappa2);
  if (shape_kappa_cmp != 0)
    return shape_kappa_cmp < 0;
  const int shape_asph_cmp =
      compare_shape_moment(lhs.shape_asphericity, rhs.shape_asphericity);
  if (shape_asph_cmp != 0)
    return shape_asph_cmp < 0;
  const int shape_acyl_cmp =
      compare_shape_moment(lhs.shape_acylindricity, rhs.shape_acylindricity);
  if (shape_acyl_cmp != 0)
    return shape_acyl_cmp < 0;
  metric_cmp = compare_pbctopo_double_metric(lhs.bbox_volume,
                                              rhs.bbox_volume);
  if (metric_cmp != 0)
    return metric_cmp < 0;
  if (lhs.hull_valid && rhs.hull_valid &&
      (metric_cmp = compare_pbctopo_double_metric(lhs.hull_volume,
                                                   rhs.hull_volume)) != 0)
    return metric_cmp < 0;
  metric_cmp = compare_pbctopo_double_metric(lhs.min_interfragment_d2,
                                              rhs.min_interfragment_d2);
  if (metric_cmp != 0)
    return metric_cmp < 0;
  metric_cmp = compare_pbctopo_double_metric(lhs.shift_norm2, rhs.shift_norm2);
  if (metric_cmp != 0)
    return metric_cmp < 0;
  return lhs.node_ids < rhs.node_ids;
}

[[maybe_unused]] void
apply_pbctopo_candidate(Multi_MOL &mol1_system, Multi_MOL &mol2_system,
                        size_t frame_index, const PbctopoCandidate &candidate,
                        const vector<PbctopoFragmentCopyNode> &nodes,
                        const vector<PbctopoAtomCopy> &atom_copies) {
  for (size_t node_id : candidate.node_ids) {
    const auto &node = nodes[node_id];
    for (size_t atom_idx = node.atom_begin; atom_idx < node.atom_end; ++atom_idx) {
      const auto &atom_copy = atom_copies[atom_idx];
      if (atom_copy.owner == 0) {
        mol1_system.m_mol[frame_index].pos.x[atom_copy.local_atom_index] =
            atom_copy.x;
        mol1_system.m_mol[frame_index].pos.y[atom_copy.local_atom_index] =
            atom_copy.y;
        mol1_system.m_mol[frame_index].pos.z[atom_copy.local_atom_index] =
            atom_copy.z;
      } else if (frame_index < mol2_system.m_mol.size()) {
        mol2_system.m_mol[frame_index].pos.x[atom_copy.local_atom_index] =
            atom_copy.x;
        mol2_system.m_mol[frame_index].pos.y[atom_copy.local_atom_index] =
            atom_copy.y;
        mol2_system.m_mol[frame_index].pos.z[atom_copy.local_atom_index] =
            atom_copy.z;
      }
    }
  }
}

} // namespace

int Multi_MOL::solve_pbc_pbctopo(Multi_MOL &mol2, arguments &Arguments,
                                 TitanExecutor *executor,
                                 size_t frame_report_offset) {
  using clock = chrono::steady_clock;
  const size_t frame_count =
      mol2.m_mol.empty() ? this->m_mol.size()
                         : min(this->m_mol.size(), mol2.m_mol.size());
  if (frame_count == 0)
    return 0;

  cache_topology_vdw_radii_from_parameters();
  mol2.cache_topology_vdw_radii_from_parameters();

  reset_pbctopo_log_file();
  ofstream pbctopo_frame_report = open_pbctopo_frame_report();
  size_t fallback_frames = 0;
  size_t valid_candidates = 0;
  PbctopoTimingStats timing;
  PbctopoWarmState warm_state;
  const PbctopoRuntimeCutoffs runtime_cutoffs =
      resolve_pbctopo_cutoffs(Arguments);
  const PbctopoSolveMode solve_mode = resolve_pbctopo_solve_mode(Arguments);
  const PbctopoAnchorBackend anchor_backend =
      resolve_pbctopo_anchor_backend(Arguments);
  const PbctopoAnchorMetricMode anchor_metric_mode =
      titan_pbctopo::resolve_anchor_metric_mode(Arguments);
  const PbctopoTemporalMode temporal_mode =
      titan_pbctopo::resolve_temporal_mode(Arguments);
  const size_t requested_exact_component_limit =
      Arguments.pbc.pbctopo_exact_component_limit_set
          ? static_cast<size_t>(Arguments.pbc.pbctopo_exact_component_limit)
          : kPbctopoExactComponentLimit;
  const size_t exact_component_limit =
      min(requested_exact_component_limit, kPbctopoExactComponentLimitMax);
  const size_t requested_beam_width =
      Arguments.pbc.pbctopo_beam_width_set
          ? static_cast<size_t>(Arguments.pbc.pbctopo_beam_width)
          : kPbctopoBeamWidth;
  const size_t beam_width =
      max<size_t>(1, min(requested_beam_width, kPbctopoBeamWidthMax));
  append_pbctopo_log_line("[pbctopo] begin run: frames=" +
                          to_string(frame_count) + ", pbctopo_mode=" +
                          pbctopo_solve_mode_name(solve_mode) +
                          ", anchor_backend=" +
                          pbctopo_anchor_backend_name(anchor_backend) +
                          ", anchor_metrics=" +
                          titan_pbctopo::anchor_metric_mode_name(
                              anchor_metric_mode) +
                          ", temporal=" +
                          titan_pbctopo::temporal_mode_name(temporal_mode) +
                          ", vdw_clash_scale=" +
                          to_string(runtime_cutoffs.vdw_clash_scale) +
                          ", bond_cutoff=" +
                          to_string(runtime_cutoffs.bond_cutoff) +
                          ", clash_cutoff=" +
                          to_string(runtime_cutoffs.clash_cutoff) +
                          ", contact_cutoff=" +
                          to_string(runtime_cutoffs.pair_contact_cutoff) +
                          ", continuity_cutoff=" +
                          to_string(runtime_cutoffs.continuity_cutoff) +
                          ", exact_component_limit=" +
                          to_string(exact_component_limit) + ", beam_width=" +
                          to_string(beam_width));
  const PbctopoTrajectoryContext trajectory_context =
      build_pbctopo_trajectory_context(*this, mol2, frame_count);
  vector<vector<PbctopoCutCandidate>> temporal_cut_candidates(frame_count);
  vector<PbctopoCutCandidate> greedy_selected_candidates(frame_count);
  struct PbctopoTemporalComponentLayout {
    vector<uint32_t> atom_component_index;
    size_t component_count = 0;
  };
  vector<PbctopoTemporalComponentLayout> temporal_component_layouts;
  vector<size_t> temporal_frame_layout(
      frame_count, numeric_limits<size_t>::max());
  vector<PbctopoFrameReport> deferred_frame_reports(frame_count);
  vector<unsigned char> deferred_frame_report_valid(frame_count, 0);
  auto write_or_defer_frame_report = [&](size_t frame_idx,
                                         const PbctopoFrameReport &report) {
    if (temporal_mode == PbctopoTemporalMode::Viterbi) {
      deferred_frame_reports[frame_idx] = report;
      deferred_frame_report_valid[frame_idx] = 1;
    } else if (pbctopo_frame_report.is_open()) {
      write_pbctopo_frame_report_row(pbctopo_frame_report, report);
    }
  };
  const bool has_ligand_static = trajectory_context.lig_static.atom_count > 0;
  const bool static_topology_ready =
      trajectory_context.rec_static.topology_graph_available &&
      (!has_ligand_static ||
       trajectory_context.lig_static.topology_graph_available);
  const bool static_metadata_ready =
      !static_topology_ready &&
      trajectory_context.rec_static.metadata_graph_available &&
      (!has_ligand_static ||
       trajectory_context.lig_static.metadata_graph_available);
  const string static_graph_source =
      static_topology_ready
          ? "gromacs_topology"
          : (static_metadata_ready ? "metadata_chain_residue"
                                   : "geometry_cutoff");
  append_pbctopo_log_line(
      "[pbctopo] local graph prepared: source=" + static_graph_source +
      ", selected_atoms=" +
      to_string(trajectory_context.rec_static.atom_count +
                trajectory_context.lig_static.atom_count) +
      ", topology_bonds=" +
      to_string(trajectory_context.rec_static.topology_bonds_available +
                trajectory_context.lig_static.topology_bonds_available) +
      ", projected_topology_bonds=" +
      to_string(
          trajectory_context.rec_static.projected_topology_bonds.size() +
          trajectory_context.lig_static.projected_topology_bonds.size()) +
      ", metadata_edges=" +
      to_string(trajectory_context.rec_static.projected_metadata_edges.size() +
                trajectory_context.lig_static.projected_metadata_edges.size()) +
      ", vdw_topology_radii=" +
      to_string(trajectory_context.rec_static.vdw_topology_count +
                trajectory_context.lig_static.vdw_topology_count) +
      ", vdw_parameter_radii=" +
      to_string(trajectory_context.rec_static.vdw_parameter_count +
                trajectory_context.lig_static.vdw_parameter_count) +
      ", vdw_inferred_radii=" +
      to_string(trajectory_context.rec_static.vdw_inferred_count +
                trajectory_context.lig_static.vdw_inferred_count));
  PbctopoFrameScratch scratch;
  scratch.reserve_from_context(trajectory_context);
  PbctopoHullScratch hull_scratch;
  hull_scratch.reserve(trajectory_context.total_atoms_hint);
  auto &base_atoms = scratch.base_atoms;
  auto &projection_atoms = scratch.projection_atoms;
  auto &component_atoms = scratch.component_atoms;
  auto &local_edges = scratch.local_edges;
  auto &component_scoring_edges = scratch.component_scoring_edges;
  auto &projection_edges = scratch.projection_edges;
  auto &projection_continuity_pairs = scratch.projection_continuity_pairs;
  auto &projection_anchor_indices = scratch.projection_anchor_indices;
  auto &local_edge_set = scratch.local_edge_set;
  auto &projection_edge_set = scratch.projection_edge_set;
  auto &components = scratch.components;
  auto &component_sizes = scratch.component_sizes;
  auto &shifts = scratch.shifts;
  auto &volume_seeds = scratch.volume_seeds;
  auto &warm_points = scratch.warm_points;
  auto &coarse_warm_points = scratch.coarse_warm_points;
  auto &placed_points = scratch.placed_points;
  auto &best_safe_points = scratch.best_safe_points;
  auto &verify_points = scratch.verify_points;
  auto &component_points = scratch.component_points;
  auto &coarse_candidates = scratch.coarse_candidates;
  auto &coarse_seed_options = scratch.coarse_seed_options;
  auto &coarse_results = scratch.coarse_results;
  auto &projection_index = scratch.projection_index;
  auto &atom_component_index = scratch.atom_component_index;
  auto &virtual_search_workspace = scratch.virtual_search_workspace;
  auto &virtual_search_plan = scratch.virtual_search_plan;
  vector<PbctopoHullSpan> hull_owner_spans;
  hull_owner_spans.reserve(8);
  vector<Vec3> hull_rescue_points;
  hull_rescue_points.reserve(trajectory_context.total_atoms_hint);
  bool projection_clash_skip_active = false;
  string llps_local_unit_source = "none";
  auto fallback_pbci_frame = [&](size_t frame_index) {
    ++timing.pbci_fallback_frames;
    int chain_ref = Arguments.pbc.chain_ref;
    const size_t chain_count = this->Index.empty() ? 0 : this->Index.size() - 1;
    if (chain_count == 0) {
      chain_ref = 0;
    } else if (chain_ref < 0 ||
               static_cast<size_t>(chain_ref) >= chain_count) {
      chain_ref = 0;
    }

    this->find_best_image_envelope(frame_index, chain_ref, Arguments);
    if (frame_index < mol2.m_mol.size() && !mol2.m_mol[frame_index].pos.x.empty())
      this->find_best_image_of_mol2(mol2, frame_index, chain_ref, Arguments);

    vector<PbctopoAtom> fallback_atoms;
    fallback_atoms.reserve(trajectory_context.total_atoms_hint);
    append_atoms_for_system(this->m_mol[frame_index], &this->info,
                            &trajectory_context.rec_static, 0, fallback_atoms);
    if (frame_index < mol2.m_mol.size() && !mol2.m_mol[frame_index].pos.x.empty()) {
      append_atoms_for_system(mol2.m_mol[frame_index], &mol2.info,
                              &trajectory_context.lig_static, 1, fallback_atoms);
    }
    prepare_pbctopo_steric_radii(fallback_atoms,
                                 runtime_cutoffs.vdw_clash_scale);
    vector<pair<int, int>> fallback_contacting_owner_pairs =
        build_contacting_owner_pairs({}, fallback_atoms,
                                     this->m_mol[frame_index].Box,
                                     runtime_cutoffs.pair_contact_cutoff,
                                     false);
    vector<Vec3> fallback_points;
    fallback_points.reserve(fallback_atoms.size());
    for (const auto &atom : fallback_atoms)
      fallback_points.push_back({atom.x, atom.y, atom.z});
    if (!reimage_contacting_owner_points_minimum_image(
            fallback_points, fallback_atoms, this->m_mol[frame_index].Box,
            fallback_contacting_owner_pairs)) {
      return;
    }
    for (size_t atom_idx = 0; atom_idx < fallback_atoms.size(); ++atom_idx) {
      const auto &atom = fallback_atoms[atom_idx];
      const auto &point = fallback_points[atom_idx];
      if (atom.owner == 0) {
        this->m_mol[frame_index].pos.x[atom.local_atom_index] = point.x;
        this->m_mol[frame_index].pos.y[atom.local_atom_index] = point.y;
        this->m_mol[frame_index].pos.z[atom.local_atom_index] = point.z;
      } else if (frame_index < mol2.m_mol.size()) {
        mol2.m_mol[frame_index].pos.x[atom.local_atom_index] = point.x;
        mol2.m_mol[frame_index].pos.y[atom.local_atom_index] = point.y;
        mol2.m_mol[frame_index].pos.z[atom.local_atom_index] = point.z;
      }
    }
    append_pbctopo_log_line(
        "[pbctopo] fallback pbci=0 reimaged contacting owners at frame " +
        to_string(static_cast<long long>(this->m_mol[frame_index].time_of_snap)));
  };

  for (size_t frame_idx = 0; frame_idx < frame_count; ++frame_idx) {
    ++timing.frames;
    const bool kinematic_seed_history_ready =
        has_consecutive_kinematic_seed_history(warm_state, frame_idx);
    const PbctopoCutCandidate kinematic_seed_previous =
        warm_state.kinematic_previous_virtual_candidate;
    const PbctopoCutCandidate kinematic_seed_latest =
        warm_state.kinematic_latest_virtual_candidate;
    PbctopoFrameReport frame_report;
    frame_report.frame = frame_report_offset + frame_idx;
    frame_report.time_ps = this->m_mol[frame_idx].time_of_snap;
    const auto frame_report_started = clock::now();
    const double frame_graph_s_before = timing.graph_s + timing.component_graph_s;
    const double frame_virtual_s_before =
        timing.event_planes_s + timing.coarse_s + timing.refine_s +
        timing.early_projection_s;
    const double frame_component_s_before =
        timing.component_prep_s + timing.component_graph_s +
        timing.component_search_s + timing.warm_component_s;
    const bool skip_projection_candidate_clash = projection_clash_skip_active;
    if (skip_projection_candidate_clash)
      ++timing.projection_clash_skip_frames;
    scratch.reset_frame();
    auto stage_started = clock::now();
    append_atoms_for_system(this->m_mol[frame_idx], &this->info,
                            &trajectory_context.rec_static, 0, base_atoms);
    if (frame_idx < mol2.m_mol.size() && !mol2.m_mol[frame_idx].pos.x.empty()) {
      append_atoms_for_system(mol2.m_mol[frame_idx], &mol2.info,
                              &trajectory_context.lig_static, 1, base_atoms);
    }
    prepare_pbctopo_steric_radii(base_atoms,
                                 runtime_cutoffs.vdw_clash_scale);
    timing.append_s += chrono::duration<double>(clock::now() - stage_started).count();

    if (base_atoms.empty()) {
      frame_report.status = PbctopoFrameStatus::Fallback;
      frame_report.mode = PbctopoFrameMode::Pbci;
      frame_report.fallback_reason = PbctopoFailureReason::NoSeed;
      frame_report.time_total_ms =
          chrono::duration<double, milli>(clock::now() - frame_report_started)
              .count();
      write_or_defer_frame_report(frame_idx, frame_report);
      reset_kinematic_seed_history(warm_state, timing);
      continue;
    }

    stage_started = clock::now();
    build_heavy_atom_subset(base_atoms, projection_atoms);
    timing.heavy_subset_s +=
        chrono::duration<double>(clock::now() - stage_started).count();
    stage_started = clock::now();
    const PbctopoTranslationCache translation_cache =
        build_virtual_box_translation_cache(this->m_mol[frame_idx].Box);
    timing.translation_cache_s +=
        chrono::duration<double>(clock::now() - stage_started).count();

    stage_started = clock::now();
    const bool has_ligand_frame =
        frame_idx < mol2.m_mol.size() && !mol2.m_mol[frame_idx].pos.x.empty();
    const bool rec_topology_graph_ready =
        trajectory_context.rec_static.topology_graph_available;
    const bool lig_topology_graph_ready =
        !has_ligand_frame ||
        trajectory_context.lig_static.topology_graph_available;
    const bool use_static_topology_graph =
        rec_topology_graph_ready && lig_topology_graph_ready;
    const bool any_topology_graph =
        rec_topology_graph_ready ||
        (has_ligand_frame &&
         trajectory_context.lig_static.topology_graph_available);
    const bool rec_metadata_graph_ready =
        !rec_topology_graph_ready &&
        trajectory_context.rec_static.metadata_graph_available;
    const bool lig_metadata_graph_ready =
        !has_ligand_frame ||
        (!trajectory_context.lig_static.topology_graph_available &&
         trajectory_context.lig_static.metadata_graph_available);
    const bool use_static_metadata_graph =
        !any_topology_graph && rec_metadata_graph_ready &&
        lig_metadata_graph_ready;
    if (use_static_topology_graph) {
      local_edges.clear();
      local_edge_set.clear();
      const size_t projected_topology_bonds =
          trajectory_context.rec_static.projected_topology_bonds.size() +
          (has_ligand_frame
               ? trajectory_context.lig_static.projected_topology_bonds.size()
               : 0);
      local_edges.reserve(projected_topology_bonds);
      local_edge_set.reserve(projected_topology_bonds * 2 + 1);
      ++timing.pbctopo_topology_static_graph_frames;
    } else if (use_static_metadata_graph) {
      local_edges.clear();
      local_edge_set.clear();
      const size_t projected_metadata_edges =
          trajectory_context.rec_static.projected_metadata_edges.size() +
          (has_ligand_frame
               ? trajectory_context.lig_static.projected_metadata_edges.size()
               : 0);
      local_edges.reserve(projected_metadata_edges);
      local_edge_set.reserve(projected_metadata_edges * 2 + 1);
      ++timing.pbctopo_metadata_graph_frames;
    } else {
      build_local_neighborhood_graph(base_atoms, this->m_mol[frame_idx].Box,
                                     local_edges, local_edge_set,
                                     runtime_cutoffs.bond_cutoff);
      timing.pbctopo_geometry_edges_added += local_edges.size();
      ++timing.pbctopo_geometry_fallback_frames;
    }
    const auto geometry_done = clock::now();
    timing.pbctopo_local_graph_geometry_s +=
        chrono::duration<double>(geometry_done - stage_started).count();
    if (use_static_metadata_graph) {
      append_metadata_graph_edges_for_system(trajectory_context.rec_static, 0,
                                             base_atoms, local_edges,
                                             local_edge_set, timing);
      if (has_ligand_frame) {
        append_metadata_graph_edges_for_system(
            trajectory_context.lig_static,
            trajectory_context.rec_static.atom_count, base_atoms, local_edges,
            local_edge_set, timing);
      }
    } else {
      append_topology_bond_edges_for_system(trajectory_context.rec_static, 0,
                                            base_atoms,
                                            this->m_mol[frame_idx].Box,
                                            local_edges, local_edge_set,
                                            timing);
      if (has_ligand_frame) {
        append_topology_bond_edges_for_system(
            trajectory_context.lig_static,
            trajectory_context.rec_static.atom_count, base_atoms,
            this->m_mol[frame_idx].Box, local_edges, local_edge_set, timing);
      }
    }
    timing.pbctopo_topology_edge_append_s +=
        chrono::duration<double>(clock::now() - geometry_done).count();
    timing.graph_s += chrono::duration<double>(clock::now() - stage_started).count();

    if (solve_mode == PbctopoSolveMode::LlpsLocal) {
      const PbctopoLlpsLocalUnits llps_units =
          build_pbctopo_llps_local_units(*this, mol2, base_atoms, local_edges);
      llps_local_unit_source =
          merge_pbctopo_llps_unit_source(llps_local_unit_source,
                                         llps_units.source);
      vector<Vec3> llps_points = unwrap_pbctopo_llps_local_points(
          base_atoms, local_edges, this->m_mol[frame_idx].Box, llps_units,
          timing);
      stage_started = clock::now();
      apply_pbctopo_points_direct(base_atoms, llps_points, *this, mol2,
                                  frame_idx);
      timing.apply_s +=
          chrono::duration<double>(clock::now() - stage_started).count();
      ++timing.apply_calls;
      ++timing.pbctopo_llps_local_frames;
      timing.pbctopo_llps_local_units += llps_units.units.size();
      ++valid_candidates;

      frame_report.status = PbctopoFrameStatus::Certified;
      frame_report.mode = PbctopoFrameMode::LlpsLocal;
      frame_report.component_count = llps_units.units.size();
      frame_report.safe_candidate_count = 1;
      frame_report.fallback_reason = PbctopoFailureReason::None;
      frame_report.time_graph_ms =
          max(0.0, (timing.graph_s + timing.component_graph_s -
                    frame_graph_s_before) *
                       1000.0);
      frame_report.time_virtual_ms =
          max(0.0, (timing.coarse_s + timing.refine_s +
                    timing.early_projection_s - frame_virtual_s_before) *
                       1000.0);
      frame_report.time_component_ms =
          max(0.0, (timing.component_prep_s + timing.component_graph_s +
                    timing.component_search_s + timing.warm_component_s -
                    frame_component_s_before) *
                       1000.0);
      frame_report.time_total_ms =
          chrono::duration<double, milli>(clock::now() - frame_report_started)
              .count();
      write_or_defer_frame_report(frame_idx, frame_report);
      reset_kinematic_seed_history(warm_state, timing);
      continue;
    }

    stage_started = clock::now();
    build_projection_edge_subset(base_atoms, projection_atoms, local_edges,
                                 projection_edges, projection_edge_set,
                                 projection_index);
    build_projection_continuity_subset(base_atoms, projection_atoms,
                                       trajectory_context.continuity_pairs,
                                       projection_continuity_pairs,
                                       projection_index);
    build_projection_anchor_subset(base_atoms,
                                   trajectory_context.continuity_anchor_indices,
                                   projection_anchor_indices, projection_index);
    const vector<pair<int, int>> local_hard_required_owner_pairs =
        build_contacting_owner_pairs(local_edges, base_atoms,
                                     this->m_mol[frame_idx].Box,
                                     runtime_cutoffs.pair_contact_cutoff,
                                     true);
    const vector<pair<int, int>> projection_hard_required_owner_pairs =
        build_contacting_owner_pairs(projection_edges, projection_atoms,
                                     this->m_mol[frame_idx].Box,
                                     runtime_cutoffs.pair_contact_cutoff,
                                     true);
    timing.pbctopo_hard_required_owner_pairs +=
        local_hard_required_owner_pairs.size();
    timing.projection_edge_s +=
        chrono::duration<double>(clock::now() - stage_started).count();

    PbctopoAnchorGraphContext base_anchor_context =
        build_anchor_graph_context(base_atoms, this->m_mol[frame_idx].Box,
                                   trajectory_context.continuity_pairs,
                                   trajectory_context.continuity_anchor_indices,
                                   runtime_cutoffs.pair_contact_cutoff,
                                   &timing);
    PbctopoAnchorGraphContext projection_anchor_context =
        build_anchor_graph_context(projection_atoms, this->m_mol[frame_idx].Box,
                                   projection_continuity_pairs,
                                   projection_anchor_indices,
                                   runtime_cutoffs.pair_contact_cutoff,
                                   &timing);
    const bool use_anchor_order =
        anchor_metric_mode == PbctopoAnchorMetricMode::Full ||
        anchor_metric_mode == PbctopoAnchorMetricMode::NoContact;
    const bool use_anchor_contact =
        anchor_metric_mode == PbctopoAnchorMetricMode::Full ||
        anchor_metric_mode == PbctopoAnchorMetricMode::NoOrder;
    base_anchor_context.use_order_metrics = use_anchor_order;
    base_anchor_context.use_contact_metrics = use_anchor_contact;
    projection_anchor_context.use_order_metrics = use_anchor_order;
    projection_anchor_context.use_contact_metrics = use_anchor_contact;

    stage_started = clock::now();
    components = build_connected_components_from_local_edges(base_atoms.size(), local_edges);
    atom_component_index =
        build_atom_component_index(base_atoms.size(), components);
    if (temporal_mode == PbctopoTemporalMode::Viterbi &&
        !atom_component_index.empty() &&
        components.size() <= numeric_limits<uint32_t>::max()) {
      vector<uint32_t> compact_layout;
      compact_layout.reserve(atom_component_index.size());
      bool layout_valid = true;
      for (size_t comp_idx : atom_component_index) {
        if (comp_idx >= components.size() ||
            comp_idx > numeric_limits<uint32_t>::max()) {
          layout_valid = false;
          break;
        }
        compact_layout.push_back(static_cast<uint32_t>(comp_idx));
      }
      if (layout_valid) {
        size_t layout_idx = 0;
        for (; layout_idx < temporal_component_layouts.size(); ++layout_idx) {
          if (temporal_component_layouts[layout_idx].component_count ==
                  components.size() &&
              temporal_component_layouts[layout_idx].atom_component_index ==
                  compact_layout) {
            break;
          }
        }
        if (layout_idx == temporal_component_layouts.size()) {
          temporal_component_layouts.push_back(
              {std::move(compact_layout), components.size()});
        }
        temporal_frame_layout[frame_idx] = layout_idx;
      }
    }
    timing.component_split_s +=
        chrono::duration<double>(clock::now() - stage_started).count();
    const size_t raw_component_count = components.size();
    const PbctopoComponentMode component_mode =
        select_pbctopo_component_mode(raw_component_count,
                                      exact_component_limit);
    const bool can_use_component_logic =
        component_mode != PbctopoComponentMode::Disabled;
    component_sizes = component_size_signature(components);
    const bool same_component_warm_layout =
        warm_state.has_component_offsets &&
        same_component_layout(components, warm_state.previous_component_sizes);
    if (!same_component_warm_layout)
      clear_component_trial_hints(warm_state);
    bool component_context_ready = false;
    auto ensure_component_context = [&]() {
      if (component_context_ready)
        return;
      const auto component_context_started = clock::now();
      component_atoms =
          unwrap_atoms_by_local_graph(base_atoms, local_edges,
                                      this->m_mol[frame_idx].Box);
      if (can_use_component_logic)
        shifts = build_topology_shifts(this->m_mol[frame_idx].Box);
      component_context_ready = true;
      ++timing.component_context_builds;
      timing.component_context_s +=
          chrono::duration<double>(clock::now() - component_context_started).count();
    };
    auto extract_offsets_and_time = [&](const vector<PbctopoAtom> &source_atoms,
                                        const vector<Vec3> &points) {
      const auto extract_started = clock::now();
      auto offsets = extract_component_offsets_from_points(
          source_atoms, components, this->m_mol[frame_idx].Box, points);
      ++timing.component_offset_extracts;
      timing.component_offset_extract_s +=
          chrono::duration<double>(clock::now() - extract_started).count();
      return offsets;
    };
    const PbctopoLatticeMetric certification_lattice(
        this->m_mol[frame_idx].Box);
    auto certify_candidate_for_output =
        [&](PbctopoCutCandidate &candidate, const vector<Vec3> &points,
            const vector<Int3> *known_offsets = nullptr) {
          ++timing.pbctopo_certification_attempts;
          vector<Int3> extracted_offsets;
          const vector<Int3> *offsets = known_offsets;
          if (offsets == nullptr || offsets->size() != components.size()) {
            ensure_component_context();
            extracted_offsets =
                extract_offsets_and_time(component_atoms, points);
            offsets = &extracted_offsets;
          }
          candidate.assignment = make_pbctopo_image_assignment(
              *offsets, components);
          if (!candidate.assignment.valid) {
            ++timing.pbctopo_certification_assignment_rejects;
            candidate.hard_feasibility = {};
            candidate.hard_feasibility.evaluated = candidate.valid;
            return false;
          }
          const size_t cycle_mismatches =
              candidate.hard_feasibility.cycle_mismatches;
          const PbctopoImageEquationResiduals image_equations =
              evaluate_pbctopo_image_equations(
                  base_atoms, points, local_edges, certification_lattice);
          const auto required_contacts =
              titan_pbctopo::evaluate_pbctopo_required_contacts(
                  base_atoms, points, local_edges);
          timing.pbctopo_certification_topology_edges_checked +=
              image_equations.topology_edges_checked;
          timing.pbctopo_certification_topology_edge_residuals +=
              image_equations.topology_edge_residuals;
          timing.pbctopo_certification_atom_image_residuals +=
              image_equations.atom_image_residuals;
          timing.pbctopo_certification_required_contacts_checked +=
              required_contacts.contacts_checked;
          timing.pbctopo_certification_required_contacts_lost +=
              required_contacts.contacts_lost;
          if (!certify_pbctopo_candidate(candidate, image_equations,
                                         cycle_mismatches,
                                         required_contacts)) {
            ++timing.pbctopo_certification_hard_rejects;
            return false;
          }
          ++timing.pbctopo_certification_successes;
          return true;
        };
    auto record_temporal_certified_candidate =
        [&](const PbctopoCutCandidate &source_candidate,
            const vector<Vec3> &candidate_points,
            const vector<Int3> *known_offsets = nullptr) {
          if (temporal_mode != PbctopoTemporalMode::Viterbi ||
              !source_candidate.valid || source_candidate.ultra_clash ||
              candidate_points.size() != base_atoms.size()) {
            return;
          }
          PbctopoCutCandidate candidate = source_candidate;
          if (!candidate.hard_feasibility.certified() ||
              !candidate.assignment.valid) {
            if (!certify_candidate_for_output(candidate, candidate_points,
                                              known_offsets)) {
              return;
            }
          }
          auto &frame_candidates = temporal_cut_candidates[frame_idx];
          const bool duplicate = any_of(
              frame_candidates.begin(), frame_candidates.end(),
              [&](const PbctopoCutCandidate &existing) {
                return same_pbctopo_image_assignment(existing.assignment,
                                                      candidate.assignment);
              });
          if (duplicate) {
            ++timing.pbctopo_temporal_assignment_deduped;
            return;
          }
          frame_candidates.push_back(std::move(candidate));
          sort(frame_candidates.begin(), frame_candidates.end(),
               better_cut_candidate);
          const size_t keep =
              static_cast<size_t>(Arguments.pbc.pbctopo_temporal_k);
          if (frame_candidates.size() > keep)
            frame_candidates.resize(keep);
          ++timing.pbctopo_temporal_certified_candidates;
        };
    auto apply_cut_candidate_timed = [&](const vector<Vec3> &points) {
      const auto apply_started = clock::now();
      apply_cut_candidate(base_atoms, points, *this, mol2, frame_idx);
      ++timing.apply_calls;
      timing.apply_s +=
          chrono::duration<double>(clock::now() - apply_started).count();
    };

    constexpr size_t kProjectionTopSeeds = 24;
    constexpr size_t kProjectionKeepCandidates = 3;
    warm_points.reserve(base_atoms.size());
    coarse_warm_points.reserve(projection_atoms.size());
    vector<Vec3> prefetched_safe_points;
    prefetched_safe_points.reserve(base_atoms.size());
    PbctopoCutCandidate prefetched_safe_candidate;
    PbctopoCutCandidate best_verified_virtual_candidate;
    PbctopoCutCandidate best_virtual_geometry_reference;
    PbctopoCutCandidate best_coarse_virtual_geometry_seed;
    vector<Vec3> split_rescue_points;
    split_rescue_points.reserve(base_atoms.size());
    PbctopoCutCandidate split_rescue_candidate;
    PbctopoOwnerPairSplitDetail split_rescue_detail;
    string split_rescue_source;
    bool split_rescue_available = false;
    auto record_split_rescue_seed =
        [&](const PbctopoCutCandidate &candidate,
            const vector<Vec3> &candidate_points,
            const PbctopoOwnerPairSplitDetail &detail,
            const string &source_label) {
          if (!candidate.valid || candidate.ultra_clash || !detail.valid ||
              candidate_points.size() != base_atoms.size()) {
            return;
          }
          const bool better_seed =
              !split_rescue_available ||
              detail.raw_norm + 1.0e-9 < split_rescue_detail.raw_norm ||
              (fabs(detail.raw_norm - split_rescue_detail.raw_norm) <= 1.0e-9 &&
               detail.minimg_norm + 1.0e-9 < split_rescue_detail.minimg_norm);
          if (!better_seed)
            return;
          split_rescue_available = true;
          split_rescue_candidate = candidate;
          split_rescue_points = candidate_points;
          split_rescue_detail = detail;
          split_rescue_source = source_label;
        };
    auto validate_rescue_candidate_points =
        [&](const vector<Vec3> &candidate_points,
            const PbctopoCutCandidate *incumbent_cut,
            PbctopoCutCandidate &candidate_out) {
          if (candidate_points.size() != base_atoms.size())
            return false;
          size_t split_owner_pairs = 0;
          double split_gap = 0.0;
          PbctopoOwnerPairSplitDetail detail;
          if (candidate_has_split_contacting_owners(
                  candidate_points, base_atoms, this->m_mol[frame_idx].Box,
                  local_hard_required_owner_pairs, split_owner_pairs, split_gap,
                  &detail)) {
            return false;
          }

          const vector<Int3> rescue_offsets = extract_component_offsets_from_points(
              base_atoms, components, this->m_mol[frame_idx].Box,
              candidate_points);
          const double rescue_shift_norm2 = pbctopo_component_shift_metric(
              this->m_mol[frame_idx].Box, rescue_offsets);

          const PbctopoLatticeMetric rescue_lattice(
              this->m_mol[frame_idx].Box);
          double min_sx = numeric_limits<double>::infinity();
          double min_sy = numeric_limits<double>::infinity();
          double min_sz = numeric_limits<double>::infinity();
          double max_sx = -numeric_limits<double>::infinity();
          double max_sy = -numeric_limits<double>::infinity();
          double max_sz = -numeric_limits<double>::infinity();
          for (const auto &point : candidate_points) {
            const auto scaled = rescue_lattice.cartesian_to_scaled(
                point.x, point.y, point.z);
            const double sx = scaled[0];
            const double sy = scaled[1];
            const double sz = scaled[2];
            min_sx = min(min_sx, sx);
            min_sy = min(min_sy, sy);
            min_sz = min(min_sz, sz);
            max_sx = max(max_sx, sx);
            max_sy = max(max_sy, sy);
            max_sz = max(max_sz, sz);
          }
          auto normalize_rescue_center = [](double center) {
            if (!isfinite(center))
              return 1.5;
            center = fmod(center - 0.5, 1.0);
            if (center < 0.0)
              center += 1.0;
            return clamp_virtual_box_center(center + 0.5);
          };
          const double rescue_cut_x =
              normalize_rescue_center(0.5 * (min_sx + max_sx));
          const double rescue_cut_y =
              normalize_rescue_center(0.5 * (min_sy + max_sy));
          const double rescue_cut_z =
              normalize_rescue_center(0.5 * (min_sz + max_sz));

          const PointStats rescue_stats = build_point_stats(candidate_points);
          candidate_out = score_cut_candidate(
              component_scoring_edges, candidate_points, rescue_cut_x,
              rescue_cut_y, rescue_cut_z, 0.0, rescue_stats,
              rescue_shift_norm2);
          double continuity_path = numeric_limits<double>::infinity();
          double continuity_max_d2 = numeric_limits<double>::infinity();
          double anchor_mst2 = numeric_limits<double>::infinity();
          double anchor_max_d2 = numeric_limits<double>::infinity();
          PbctopoAnchorGraphScore anchor_graph_score;
          const bool rescue_ok =
              candidate_out.valid &&
              continuity_metrics_ok(candidate_points,
                                    trajectory_context.continuity_pairs,
                                    runtime_cutoffs.continuity_cutoff,
                                    continuity_path, continuity_max_d2,
                                    incumbent_cut) &&
              anchor_graph_metrics_enhanced(
                  candidate_points, trajectory_context.continuity_anchor_indices,
                  trajectory_context.continuity_pairs, anchor_mst2,
                  anchor_max_d2, anchor_graph_score, incumbent_cut,
                  PbctopoRefineScorePhase::Exact, anchor_backend, nullptr,
                  &base_anchor_context);
          if (!rescue_ok)
            return false;

          candidate_out.continuity_path = continuity_path;
          candidate_out.continuity_max_d2 = continuity_max_d2;
          candidate_out.anchor_mst2 = anchor_mst2;
          candidate_out.anchor_max_d2 = anchor_max_d2;
          candidate_out.anchor_graph = anchor_graph_score;
          analyze_cut_candidate_ultraclash(
              candidate_points, local_edge_set, candidate_out,
              runtime_cutoffs.clash_cutoff, &base_atoms,
              &atom_component_index);
          return !candidate_out.ultra_clash;
        };
    if (warm_state.has_projection_candidate) {
      const auto early_projection_started = clock::now();
      PbctopoCutCandidate early_projection_coarse = evaluate_virtual_box_candidate(
          projection_atoms, projection_edges, projection_edge_set,
          translation_cache, this->m_mol[frame_idx].Box,
          warm_state.previous_projection_candidate.cut_x,
          warm_state.previous_projection_candidate.cut_y,
          warm_state.previous_projection_candidate.cut_z,
          warm_state.previous_projection_candidate.gap_score, coarse_warm_points,
          runtime_cutoffs.clash_cutoff, &projection_continuity_pairs,
          &projection_anchor_indices, &projection_hard_required_owner_pairs,
          runtime_cutoffs.continuity_cutoff, nullptr,
          skip_projection_candidate_clash, numeric_limits<double>::quiet_NaN(),
          nullptr, &projection_anchor_context,
          PbctopoVirtualSeedSource::Previous);
      if (skip_projection_candidate_clash)
        ++timing.projection_clash_skip_evals;
      if (early_projection_coarse.valid && !early_projection_coarse.ultra_clash) {
        PbctopoCutCandidate early_projection = evaluate_virtual_box_candidate(
            base_atoms, local_edges, local_edge_set, translation_cache,
            this->m_mol[frame_idx].Box,
            warm_state.previous_projection_candidate.cut_x,
            warm_state.previous_projection_candidate.cut_y,
            warm_state.previous_projection_candidate.cut_z,
            warm_state.previous_projection_candidate.gap_score, warm_points,
            runtime_cutoffs.clash_cutoff, &trajectory_context.continuity_pairs,
            &trajectory_context.continuity_anchor_indices,
            &local_hard_required_owner_pairs, runtime_cutoffs.continuity_cutoff,
            nullptr, false, this->m_mol[frame_idx].time_of_snap,
            "warm_projection", &base_anchor_context,
            PbctopoVirtualSeedSource::Previous, &atom_component_index);
        if (early_projection.valid && !early_projection.ultra_clash &&
          candidate_close_enough_to_skip_search(
              early_projection, warm_state.previous_projection_candidate)) {
          ++timing.early_projection_used;
          warm_state.previous_projection_candidate = early_projection;
          warm_state.has_projection_candidate = true;
          if (better_cut_candidate(early_projection,
                                   best_verified_virtual_candidate)) {
            best_verified_virtual_candidate = early_projection;
          }
          if (better_cut_candidate(early_projection, prefetched_safe_candidate)) {
            prefetched_safe_candidate = early_projection;
            prefetched_safe_points = warm_points;
            frame_report.mode = PbctopoFrameMode::EarlyProjection;
          }
        }
      }
      timing.early_projection_s +=
          chrono::duration<double>(clock::now() - early_projection_started)
              .count();
    }
    if (can_use_component_logic && warm_state.has_component_offsets &&
        same_component_layout(components, warm_state.previous_component_sizes)) {
      const auto warm_component_started = clock::now();
      ensure_component_context();
      vector<Int3> warm_offsets = warm_state.previous_component_offsets;
      rebase_component_offsets(warm_offsets, this->m_mol[frame_idx].Box);
      PbctopoCutCandidate early_component = score_component_offset_solution(
          component_atoms, components, this->m_mol[frame_idx].Box, warm_offsets,
          local_edges, local_edge_set, runtime_cutoffs.clash_cutoff,
          trajectory_context.continuity_pairs,
          trajectory_context.continuity_anchor_indices,
          runtime_cutoffs.continuity_cutoff, warm_points,
          make_component_refine_callbacks(),
          &warm_state.previous_projection_candidate, nullptr, true,
          PbctopoRefineScorePhase::Exact, nullptr, nullptr, nullptr,
          anchor_backend, &base_anchor_context);
      PbctopoOwnerPairSplitDetail early_component_split_detail;
      size_t early_component_split_owner_pairs = 0;
      double early_component_split_gap = 0.0;
      const bool early_component_split =
          early_component.valid && !early_component.ultra_clash &&
          candidate_has_split_contacting_owners(
              warm_points, base_atoms, this->m_mol[frame_idx].Box,
              local_hard_required_owner_pairs, early_component_split_owner_pairs,
              early_component_split_gap, &early_component_split_detail);
      if (early_component_split && early_component_split_detail.valid) {
        record_split_rescue_seed(early_component, warm_points,
                                 early_component_split_detail,
                                 "warm_component");
        append_pbctopo_log_line(
            "[pbctopo] reject split_contacting_owner_pair at frame " +
            to_string(static_cast<long long>(this->m_mol[frame_idx].time_of_snap)) +
            " via warm_component: owners=" +
            to_string(early_component_split_detail.owner_a) + "-" +
            to_string(early_component_split_detail.owner_b) +
            ", raw_delta_scaled=" +
            format_vec3(early_component_split_detail.raw_delta_scaled) +
            ", min_image_delta_scaled=" +
            format_vec3(early_component_split_detail.minimg_delta_scaled));
      }
      if (!early_component_split && early_component.valid &&
          !early_component.ultra_clash &&
          candidate_close_enough_to_skip_search(
              early_component, warm_state.previous_projection_candidate)) {
        ++timing.component_warm_used;
        warm_state.previous_projection_candidate = early_component;
        warm_state.previous_component_offsets = warm_offsets;
        warm_state.previous_component_sizes = component_sizes;
        warm_state.has_projection_candidate = true;
        warm_state.has_component_offsets = true;
        if (better_cut_candidate(early_component, prefetched_safe_candidate)) {
          prefetched_safe_candidate = early_component;
          prefetched_safe_points = warm_points;
          frame_report.mode = PbctopoFrameMode::Warm;
        }
      }
      timing.warm_component_s +=
          chrono::duration<double>(clock::now() - warm_component_started)
              .count();
    }

    component_scoring_edges = local_edges;
    bool primary_component_graph_used = false;
    if (can_use_component_logic &&
        component_mode == PbctopoComponentMode::Exact) {
      ++timing.component_graph_primary_attempts;
      ensure_component_context();
      PbctopoCutCandidate primary_cut;
      double primary_shift_norm2 = 0.0;
      vector<Int3> primary_tree_offsets;
      vector<Int3> primary_refined_offsets;
      vector<PbctopoGraphEdge> primary_tree_edges;
      PbctopoCycleSyncMetrics primary_cycle_metrics;
      const PbctopoCutCandidate *primary_warm_reference =
          warm_state.has_projection_candidate &&
                  warm_state.previous_projection_candidate.valid &&
                  !warm_state.previous_projection_candidate.ultra_clash
              ? &warm_state.previous_projection_candidate
              : nullptr;
      auto try_primary_component_shell = [&](int shell_radius) {
        if (shell_radius > 1)
          shifts = build_topology_shifts(this->m_mol[frame_idx].Box,
                                         shell_radius);
        timing.component_image_shell_max_radius =
            max(timing.component_image_shell_max_radius,
                static_cast<size_t>(shell_radius));
        timing.component_image_shell_max_shift_count =
            max(timing.component_image_shell_max_shift_count, shifts.size());
        primary_cut = {};
        primary_shift_norm2 = 0.0;
        primary_tree_offsets.clear();
        primary_refined_offsets.clear();
        primary_tree_edges.clear();
        primary_cycle_metrics = {};
        component_points.clear();

        const auto primary_graph_started = clock::now();
        const bool graph_ok = solve_component_offsets_via_pair_graph(
            component_atoms, this->m_mol[frame_idx].Box, components, shifts,
            component_scoring_edges, local_edge_set, component_points,
            trajectory_context.continuity_pairs,
            trajectory_context.continuity_anchor_indices, runtime_cutoffs,
            primary_cut, primary_shift_norm2, beam_width, &timing,
            make_component_refine_callbacks(), anchor_backend,
            same_component_warm_layout
                ? &warm_state.previous_component_offsets
                : nullptr,
            &warm_state, primary_warm_reference, &primary_tree_offsets,
            &primary_refined_offsets, &primary_tree_edges,
            &primary_cycle_metrics, &base_anchor_context);
        timing.component_graph_s +=
            chrono::duration<double>(clock::now() - primary_graph_started)
                .count();
        frame_report.cycle_support = primary_cycle_metrics.cycle_support;
        frame_report.cycle_mismatch = primary_cycle_metrics.cycle_mismatch;
        frame_report.inconsistent_edge_count =
            primary_cycle_metrics.inconsistent_edge_count;
        frame_report.fundamental_cycle_count =
            primary_cycle_metrics.fundamental_cycle_count;
        frame_report.fundamental_cycle_residual_rms =
            primary_cycle_metrics.fundamental_cycle_residual_rms;
        frame_report.fundamental_cycle_residual_max =
            primary_cycle_metrics.fundamental_cycle_residual_max;
        primary_cut.hard_feasibility.cycle_mismatches =
            primary_cycle_metrics.inconsistent_edge_count;

        size_t split_pairs = 0;
        double split_gap = 0.0;
        const bool split =
            graph_ok && primary_cut.valid && !primary_cut.ultra_clash &&
            component_points.size() == base_atoms.size() &&
            candidate_has_split_contacting_owners(
                component_points, base_atoms, this->m_mol[frame_idx].Box,
                local_hard_required_owner_pairs, split_pairs, split_gap);
        return graph_ok && primary_cut.valid && !primary_cut.ultra_clash &&
               !split && component_points.size() == base_atoms.size() &&
               primary_refined_offsets.size() == components.size() &&
               certify_candidate_for_output(primary_cut, component_points,
                                            &primary_refined_offsets);
      };

      bool primary_certified = try_primary_component_shell(1);
      if (!primary_certified) {
        ++timing.component_image_shell_expansions;
        primary_certified = try_primary_component_shell(2);
        if (primary_certified)
          ++timing.component_image_shell_retry_successes;
      }
      if (primary_certified) {
        record_temporal_certified_candidate(primary_cut, component_points,
                                            &primary_refined_offsets);
        prefetched_safe_candidate = primary_cut;
        prefetched_safe_points = component_points;
        warm_state.previous_component_offsets = primary_refined_offsets;
        warm_state.previous_component_sizes = component_sizes;
        warm_state.has_component_offsets = true;
        frame_report.mode = PbctopoFrameMode::ComponentGraph;
        primary_component_graph_used = true;
        ++timing.component_graph_primary_certified;
        ++timing.component_graph_primary_virtual_skips;
        ++timing.component_graph_used;
      } else {
        ++timing.component_graph_primary_fallbacks;
      }
    }

    vector<PbctopoVirtualDiscreteSeedHint> virtual_seed_hints;
    if (!primary_component_graph_used) {
      stage_started = clock::now();
      collect_projection_volume_seeds(projection_atoms, 16,
                                      kProjectionTopSeeds, 0.22,
                                      volume_seeds);
    if (warm_state.has_projection_candidate) {
      PbctopoSeed3D previous_seed;
      previous_seed.cx =
          clamp_virtual_box_center(warm_state.previous_projection_candidate.cut_x);
      previous_seed.cy =
          clamp_virtual_box_center(warm_state.previous_projection_candidate.cut_y);
      previous_seed.cz =
          clamp_virtual_box_center(warm_state.previous_projection_candidate.cut_z);
      previous_seed.score =
          max(1.0, warm_state.previous_projection_candidate.gap_score + 1.0);
      previous_seed.source = PbctopoVirtualSeedSource::Previous;
      volume_seeds.push_back(previous_seed);
      ++timing.previous_seed_used;
    }
    if (kinematic_seed_history_ready) {
      const PbctopoVirtualSeedCenter predicted_center =
          extrapolate_virtual_seed_center(
              virtual_seed_center_from_candidate(kinematic_seed_previous),
              virtual_seed_center_from_candidate(kinematic_seed_latest),
              PbctopoLatticeMetric(this->m_mol[frame_idx].Box));
      PbctopoSeed3D kinematic_seed;
      kinematic_seed.cx = predicted_center.cx;
      kinematic_seed.cy = predicted_center.cy;
      kinematic_seed.cz = predicted_center.cz;
      kinematic_seed.score =
          max(1.0, kinematic_seed_latest.gap_score + 0.5);
      kinematic_seed.source = PbctopoVirtualSeedSource::Kinematic;
      volume_seeds.push_back(kinematic_seed);
      ++timing.pbctopo_kinematic_seed_candidates;
    }
    timing.seed_s += chrono::duration<double>(clock::now() - stage_started).count();

    virtual_seed_hints.reserve(volume_seeds.size());
    for (const auto &seed : volume_seeds) {
      virtual_seed_hints.push_back(
          {seed.cx, seed.cy, seed.cz, seed.score, seed.source});
    }
    stage_started = clock::now();
    PbctopoVirtualDiscreteSearchOptions discrete_options;
    discrete_options.top_k_per_axis = 4;
    discrete_options.max_exhaustive_states = 4096;
    build_pbctopo_virtual_discrete_search_plan(
        projection_atoms, projection_edges, this->m_mol[frame_idx].Box,
        virtual_seed_hints, discrete_options, virtual_search_workspace,
        virtual_search_plan);
    const double seam_sweep_s =
        chrono::duration<double>(clock::now() - stage_started).count();
    timing.event_planes_s += seam_sweep_s;
    timing.pbctopo_virtual_seam_sweep_s += seam_sweep_s;
    ++timing.pbctopo_virtual_discrete_frames;
    timing.pbctopo_virtual_axis_states_x +=
        virtual_search_plan.x_states.size();
    timing.pbctopo_virtual_axis_states_y +=
        virtual_search_plan.y_states.size();
    timing.pbctopo_virtual_axis_states_z +=
        virtual_search_plan.z_states.size();
    timing.pbctopo_virtual_sweep_events_x +=
        virtual_search_plan.event_count_x;
    timing.pbctopo_virtual_sweep_events_y +=
        virtual_search_plan.event_count_y;
    timing.pbctopo_virtual_sweep_events_z +=
        virtual_search_plan.event_count_z;
    if (virtual_search_plan.workspace_reused)
      ++timing.pbctopo_virtual_sweep_workspace_reuse_frames;
    timing.pbctopo_virtual_sweep_workspace_peak_bytes =
        max(timing.pbctopo_virtual_sweep_workspace_peak_bytes,
            virtual_search_plan.workspace_bytes_estimate);
    if (numeric_limits<uint64_t>::max() - timing.pbctopo_virtual_state_space <
        virtual_search_plan.total_state_space) {
      timing.pbctopo_virtual_state_space = numeric_limits<uint64_t>::max();
    } else {
      timing.pbctopo_virtual_state_space +=
          virtual_search_plan.total_state_space;
    }
    timing.pbctopo_virtual_states_evaluated +=
        virtual_search_plan.states.size();
    if (virtual_search_plan.search_exhaustive)
      ++timing.pbctopo_virtual_search_exhaustive_frames;
    else
      ++timing.pbctopo_virtual_search_shortlist_frames;
    if (virtual_search_plan.branch_and_bound_used)
      ++timing.pbctopo_virtual_branch_and_bound_frames;
    if (isfinite(virtual_search_plan.omitted_seam_lower_bound)) {
      timing.pbctopo_virtual_omitted_seam_lower_bound_min =
          min(timing.pbctopo_virtual_omitted_seam_lower_bound_min,
              virtual_search_plan.omitted_seam_lower_bound);
    }
    if (isfinite(virtual_search_plan.lower_bound_gap)) {
      timing.pbctopo_virtual_lower_bound_gap_min =
          min(timing.pbctopo_virtual_lower_bound_gap_min,
              virtual_search_plan.lower_bound_gap);
    }
    frame_report.virtual_search_exhaustive =
        virtual_search_plan.search_exhaustive;
    frame_report.virtual_branch_and_bound_used =
        virtual_search_plan.branch_and_bound_used;
    frame_report.virtual_state_space = virtual_search_plan.total_state_space;
    frame_report.virtual_states_evaluated = virtual_search_plan.states.size();
    frame_report.virtual_omitted_seam_lower_bound =
        virtual_search_plan.omitted_seam_lower_bound;
    frame_report.virtual_lower_bound_gap = virtual_search_plan.lower_bound_gap;
    }

    PbctopoCutCandidate best_candidate;
    PbctopoCutCandidate best_safe_candidate;
    bool best_safe_certification_rejected = false;
    size_t safe_candidates = 0;
    PbctopoCutCandidate report_best_candidate;
    PbctopoCutCandidate report_second_candidate;
    bool frame_anchor_ambiguity_counted = false;
    auto record_report_safe_candidate = [&](const PbctopoCutCandidate &candidate) {
      if (!candidate.valid)
        return;
      ++frame_report.safe_candidate_count;
      if (better_cut_candidate(candidate, report_best_candidate)) {
        report_second_candidate = report_best_candidate;
        report_best_candidate = candidate;
      } else if (better_cut_candidate(candidate, report_second_candidate)) {
        report_second_candidate = candidate;
      }
    };
    auto record_frame_failure = [&](PbctopoFailureReason reason) {
      if (frame_report.fallback_reason == PbctopoFailureReason::None)
        frame_report.fallback_reason = reason;
    };
    auto update_kinematic_history_for_frame = [&]() {
      if (best_verified_virtual_candidate.valid &&
          !best_verified_virtual_candidate.ultra_clash) {
        update_kinematic_seed_history(warm_state,
                                      best_verified_virtual_candidate,
                                      frame_idx, timing);
      } else {
        reset_kinematic_seed_history(warm_state, timing);
      }
    };
    auto finalize_frame_report =
        [&](PbctopoFrameStatus status, PbctopoFrameMode mode,
            PbctopoFailureReason reason,
            const PbctopoCutCandidate *selected_candidate) {
          if (selected_candidate != nullptr && selected_candidate->valid) {
            greedy_selected_candidates[frame_idx] = *selected_candidate;
            frame_report.hard_feasibility =
                selected_candidate->hard_feasibility;
            frame_report.assignment_signature =
                selected_candidate->assignment.signature;
            frame_report.assignment_component_count =
                selected_candidate->assignment.component_offsets.size();
          }
          frame_report.status = status;
          if (frame_report.mode == PbctopoFrameMode::None ||
              status != PbctopoFrameStatus::Certified) {
            frame_report.mode = mode;
          }
          if (reason != PbctopoFailureReason::None)
            frame_report.fallback_reason = reason;
          else if (status != PbctopoFrameStatus::Fallback)
            frame_report.fallback_reason = PbctopoFailureReason::None;
          frame_report.component_count = components.size();
          frame_report.safe_candidate_count =
              max(frame_report.safe_candidate_count, safe_candidates);

          const PbctopoCutCandidate *best_for_report =
              selected_candidate != nullptr && selected_candidate->valid
                  ? selected_candidate
                  : (report_best_candidate.valid ? &report_best_candidate
                                                 : nullptr);
          const PbctopoCutCandidate *second_for_report = nullptr;
          if (best_for_report == selected_candidate &&
              report_best_candidate.valid) {
            const bool report_best_differs =
                better_cut_candidate(*best_for_report, report_best_candidate) ||
                better_cut_candidate(report_best_candidate, *best_for_report);
            second_for_report = report_best_differs
                                    ? &report_best_candidate
                                    : (report_second_candidate.valid
                                           ? &report_second_candidate
                                           : nullptr);
          } else if (report_second_candidate.valid) {
            second_for_report = &report_second_candidate;
          }
          if (best_for_report != nullptr && best_for_report->valid) {
            frame_report.anchor_graph = best_for_report->anchor_graph;
            if (second_for_report != nullptr && second_for_report->valid &&
                candidate_anchor_prefix_equivalent(
                    *best_for_report, *second_for_report) &&
                (better_cut_candidate(*best_for_report,
                                      *second_for_report) ||
                 better_cut_candidate(*second_for_report,
                                      *best_for_report))) {
              frame_report.anchor_graph.ambiguity_count = 1;
              if (!frame_anchor_ambiguity_counted) {
                ++timing.pbctopo_anchor_ambiguous_frames;
                frame_anchor_ambiguity_counted = true;
              }
            }
            frame_report.best_score = make_score_vector(*best_for_report);
            frame_report.best_score.anchor_graph.ambiguity_count =
                frame_report.anchor_graph.ambiguity_count;
          }
          if (second_for_report != nullptr && second_for_report->valid) {
            frame_report.second_score = make_score_vector(*second_for_report);
            frame_report.second_score.anchor_graph.ambiguity_count =
                frame_report.anchor_graph.ambiguity_count;
          }
          if (best_for_report != nullptr && best_for_report->valid &&
              second_for_report != nullptr && second_for_report->valid) {
            const PbctopoLexicographicGap score_gap =
                candidate_lexicographic_gap(*best_for_report,
                                            *second_for_report);
            frame_report.score_gap_metric = score_gap.metric;
            frame_report.score_gap = score_gap.gap;
          }
          if (best_for_report != nullptr && best_for_report->valid) {
            if (isfinite(best_for_report->min_pair_d2)) {
              double active_clash_cutoff = runtime_cutoffs.clash_cutoff;
              if (best_for_report->clash_atom_a < base_atoms.size() &&
                  best_for_report->clash_atom_b < base_atoms.size()) {
                active_clash_cutoff =
                    titan_pbctopo::pbctopo_pair_clash_threshold(
                        base_atoms[best_for_report->clash_atom_a],
                        base_atoms[best_for_report->clash_atom_b],
                        runtime_cutoffs.clash_cutoff)
                        .cutoff;
              }
              frame_report.clash_margin_A =
                  sqrt(max(0.0, best_for_report->min_pair_d2)) -
                  active_clash_cutoff;
            }
            if (isfinite(best_for_report->continuity_max_d2)) {
              frame_report.continuity_margin_A =
                  runtime_cutoffs.continuity_cutoff -
                  sqrt(max(0.0, best_for_report->continuity_max_d2));
            }
            if (isfinite(best_for_report->anchor_max_d2)) {
              frame_report.anchor_max_edge_A =
                  sqrt(max(0.0, best_for_report->anchor_max_d2));
            }
            frame_report.anchor_graph.ambiguity_count =
                frame_report.best_score.anchor_graph.ambiguity_count;
            frame_report.split_owner_gap_scaled =
                best_for_report->max_owner_raw_minus_minimg_delta;
          }
          frame_report.time_graph_ms =
              max(0.0, (timing.graph_s + timing.component_graph_s -
                        frame_graph_s_before) *
                           1000.0);
          frame_report.time_virtual_ms =
              max(0.0, (timing.event_planes_s + timing.coarse_s +
                        timing.refine_s + timing.early_projection_s -
                        frame_virtual_s_before) *
                           1000.0);
          frame_report.time_component_ms =
              max(0.0, (timing.component_prep_s + timing.component_graph_s +
                        timing.component_search_s + timing.warm_component_s -
                        frame_component_s_before) *
                           1000.0);
          frame_report.time_total_ms =
              chrono::duration<double, milli>(clock::now() - frame_report_started)
                  .count();
          write_or_defer_frame_report(frame_idx, frame_report);
        };
    auto accept_rescue_candidate =
        [&](const vector<Vec3> &candidate_points,
            PbctopoCutCandidate candidate) {
          if (!certify_candidate_for_output(candidate, candidate_points)) {
            ++timing.pbctopo_diagnostic_only_candidates;
            record_frame_failure(PbctopoFailureReason::CertificationReject);
            return false;
          }
          record_temporal_certified_candidate(candidate, candidate_points);
          record_report_safe_candidate(candidate);
          best_safe_candidate = candidate;
          best_candidate = candidate;
          best_safe_points = candidate_points;
          safe_candidates = max<size_t>(safe_candidates, 1);
          projection_clash_skip_active =
              candidate_has_projection_clash_margin(best_safe_candidate,
                                                    runtime_cutoffs.clash_cutoff);
          warm_state.previous_projection_candidate = best_safe_candidate;
          warm_state.has_projection_candidate = true;
          update_kinematic_history_for_frame();
          apply_cut_candidate_timed(best_safe_points);
          return true;
        };
    placed_points.reserve(projection_atoms.size());
    best_safe_points.reserve(base_atoms.size());
    coarse_candidates.reserve(12);
    coarse_seed_options.reserve(virtual_search_plan.states.size());
    if (prefetched_safe_candidate.valid && !prefetched_safe_candidate.ultra_clash &&
        prefetched_safe_points.size() == base_atoms.size()) {
      best_candidate = prefetched_safe_candidate;
      best_safe_candidate = prefetched_safe_candidate;
      best_safe_points = prefetched_safe_points;
      safe_candidates = 1;
      record_temporal_certified_candidate(prefetched_safe_candidate,
                                          prefetched_safe_points);
      record_report_safe_candidate(prefetched_safe_candidate);
    }

    stage_started = clock::now();
    const bool has_kinematic_hint = any_of(
        virtual_seed_hints.begin(), virtual_seed_hints.end(),
        [](const PbctopoVirtualDiscreteSeedHint &hint) {
          return hint.source == PbctopoVirtualSeedSource::Kinematic;
        });
    size_t kinematic_states_evaluated = 0;
    for (size_t state_idx = 0;
         state_idx < virtual_search_plan.states.size(); ++state_idx) {
      const auto &state = virtual_search_plan.states[state_idx];
      PbctopoSeed3D seed;
      seed.cx = state.cx;
      seed.cy = state.cy;
      seed.cz = state.cz;
      // Density and temporal hints order exact states only; they are not part
      // of the physical candidate score.
      seed.score = 0.0;
      seed.ix = state.ix;
      seed.iy = state.iy;
      seed.iz = state.iz;
      seed.source = state.source;
      if (seed.source == PbctopoVirtualSeedSource::Kinematic) {
        ++timing.pbctopo_kinematic_seed_evaluated;
        ++kinematic_states_evaluated;
      }
      coarse_seed_options.push_back(
          {static_cast<uint64_t>(state_idx), seed});
    }
    if (has_kinematic_hint && kinematic_states_evaluated == 0)
      ++timing.pbctopo_kinematic_seed_deduped;
    timing.coarse_candidates += coarse_seed_options.size();
    const size_t requested_threads =
        (Arguments.nt > 1) ? Arguments.nt : static_cast<size_t>(1);
    const size_t executor_capacity =
        (executor != nullptr) ? max<size_t>(1, executor->capacity()) : requested_threads;
    const size_t coarse_threads =
        min(executor_capacity, coarse_seed_options.size());
    atomic<size_t> coarse_projection_skip_evals{0};
    coarse_results.assign(coarse_seed_options.size(), PbctopoCutCandidate{});
    if (executor != nullptr && coarse_threads > 1 && coarse_seed_options.size() >= 4 &&
        projection_atoms.size() >= 2000) {
      executor->parallel_for(
          coarse_seed_options.size(), coarse_threads, 1,
          [&](size_t, size_t begin, size_t end) {
            vector<Vec3> local_points;
            local_points.reserve(projection_atoms.size());
            for (size_t option_idx = begin; option_idx < end; ++option_idx) {
              const auto &seed = coarse_seed_options[option_idx].seed;
              coarse_results[option_idx] = evaluate_virtual_box_candidate(
                  projection_atoms, projection_edges, projection_edge_set,
                  translation_cache, this->m_mol[frame_idx].Box, seed.cx,
                  seed.cy, seed.cz, seed.score,
                  local_points, runtime_cutoffs.clash_cutoff,
                  &projection_continuity_pairs, &projection_anchor_indices,
                  &projection_hard_required_owner_pairs,
                  runtime_cutoffs.continuity_cutoff, nullptr,
                  skip_projection_candidate_clash,
                  numeric_limits<double>::quiet_NaN(), nullptr,
                  &projection_anchor_context, seed.source);
              if (skip_projection_candidate_clash)
                coarse_projection_skip_evals.fetch_add(
                    1, memory_order_relaxed);
            }
          });
    } else if (coarse_threads > 1 && coarse_seed_options.size() >= 4 &&
               projection_atoms.size() >= 2000) {
      vector<titan::jthread> workers;
      workers.reserve(coarse_threads);
      const size_t chunk =
          (coarse_seed_options.size() + coarse_threads - 1) / coarse_threads;
      for (size_t thread_idx = 0; thread_idx < coarse_threads; ++thread_idx) {
        const size_t begin = thread_idx * chunk;
        const size_t end = min(begin + chunk, coarse_seed_options.size());
        if (begin >= end)
          continue;
        workers.emplace_back([&, begin, end]() {
          vector<Vec3> local_points;
          local_points.reserve(projection_atoms.size());
          for (size_t option_idx = begin; option_idx < end; ++option_idx) {
            const auto &seed = coarse_seed_options[option_idx].seed;
            coarse_results[option_idx] = evaluate_virtual_box_candidate(
                projection_atoms, projection_edges, projection_edge_set,
                translation_cache, this->m_mol[frame_idx].Box, seed.cx, seed.cy,
                seed.cz, seed.score,
                local_points, runtime_cutoffs.clash_cutoff,
                &projection_continuity_pairs, &projection_anchor_indices,
                &projection_hard_required_owner_pairs,
                runtime_cutoffs.continuity_cutoff, nullptr,
                skip_projection_candidate_clash,
                numeric_limits<double>::quiet_NaN(), nullptr,
                &projection_anchor_context, seed.source);
            if (skip_projection_candidate_clash)
              coarse_projection_skip_evals.fetch_add(
                  1, memory_order_relaxed);
          }
        });
      }
    } else {
      for (size_t option_idx = 0; option_idx < coarse_seed_options.size();
           ++option_idx) {
        const auto &seed = coarse_seed_options[option_idx].seed;
        coarse_results[option_idx] = evaluate_virtual_box_candidate(
            projection_atoms, projection_edges, projection_edge_set,
            translation_cache, this->m_mol[frame_idx].Box, seed.cx, seed.cy,
            seed.cz, seed.score,
            placed_points, runtime_cutoffs.clash_cutoff,
            &projection_continuity_pairs, &projection_anchor_indices,
            &projection_hard_required_owner_pairs, runtime_cutoffs.continuity_cutoff,
            nullptr, skip_projection_candidate_clash,
            numeric_limits<double>::quiet_NaN(), nullptr,
            &projection_anchor_context, seed.source);
        if (skip_projection_candidate_clash)
          coarse_projection_skip_evals.fetch_add(1, memory_order_relaxed);
      }
    }
    timing.projection_clash_skip_evals +=
        coarse_projection_skip_evals.load(memory_order_relaxed);
    for (const auto &candidate : coarse_results) {
      if (better_preclash_cut_candidate(
              candidate, best_coarse_virtual_geometry_seed)) {
        best_coarse_virtual_geometry_seed = candidate;
      }
      push_top_virtual_box_candidate(coarse_candidates, candidate,
                                     kProjectionKeepCandidates);
    }
    timing.coarse_s += chrono::duration<double>(clock::now() - stage_started).count();

    stage_started = clock::now();
    verify_points.reserve(base_atoms.size());
    if (best_coarse_virtual_geometry_seed.valid) {
      ++timing.pbctopo_verify_geometry_buffer_reuses;
      ++timing.pbctopo_virtual_geometry_reference_evals;
      PbctopoCutCandidate geometry_reference = evaluate_virtual_box_candidate(
          base_atoms, local_edges, local_edge_set, translation_cache,
          this->m_mol[frame_idx].Box,
          best_coarse_virtual_geometry_seed.cut_x,
          best_coarse_virtual_geometry_seed.cut_y,
          best_coarse_virtual_geometry_seed.cut_z,
          best_coarse_virtual_geometry_seed.gap_score,
          verify_points, runtime_cutoffs.clash_cutoff,
          &trajectory_context.continuity_pairs,
          &trajectory_context.continuity_anchor_indices,
          &local_hard_required_owner_pairs, runtime_cutoffs.continuity_cutoff,
          nullptr, false, this->m_mol[frame_idx].time_of_snap,
          "virtual_geometry_reference", &base_anchor_context,
          best_coarse_virtual_geometry_seed.virtual_seed_source,
          &atom_component_index);
      if (geometry_reference.valid) {
        if (geometry_reference.ultra_clash) {
          ++timing.pbctopo_virtual_geometry_reference_clashing;
        } else {
          best_virtual_geometry_reference = geometry_reference;
          ++timing.pbctopo_virtual_geometry_reference_safe;
          ++safe_candidates;
          record_temporal_certified_candidate(geometry_reference,
                                              verify_points);
          record_report_safe_candidate(geometry_reference);
          best_verified_virtual_candidate = geometry_reference;
          if (better_cut_candidate(geometry_reference, best_safe_candidate)) {
            best_safe_candidate = geometry_reference;
            best_safe_points = verify_points;
            frame_report.mode = PbctopoFrameMode::Virtual;
          }
          if (better_cut_candidate(geometry_reference, best_candidate)) {
            best_candidate = geometry_reference;
          }
        }
      }
    }
    size_t refine_limit = coarse_candidates.size();
    if (coarse_candidates.size() >= 2) {
      const auto &top1 = coarse_candidates[0];
      const auto &top2 = coarse_candidates[1];
      const bool top1_clear_principal =
          compare_principal_rg2(top1.principal_rg2_max, top2.principal_rg2_max) <
              0 &&
          compare_principal_rg2(top1.principal_rg2_mid, top2.principal_rg2_mid) <
              0 &&
          compare_principal_rg2(top1.principal_rg2_min, top2.principal_rg2_min) <
              0;
      if (top1.broken_edges + 1 < top2.broken_edges ||
          top1.edge_stretch + 1.0e-6 < top2.edge_stretch * 0.85 ||
          top1_clear_principal) {
        refine_limit = 1;
      } else if (coarse_candidates.size() >= 3) {
        const auto &top3 = coarse_candidates[2];
        const bool top2_clear_principal =
            compare_principal_rg2(top2.principal_rg2_max,
                                  top3.principal_rg2_max) < 0 &&
            compare_principal_rg2(top2.principal_rg2_mid,
                                  top3.principal_rg2_mid) < 0 &&
            compare_principal_rg2(top2.principal_rg2_min,
                                  top3.principal_rg2_min) < 0;
        if (top2.broken_edges + 1 < top3.broken_edges ||
            top2.edge_stretch + 1.0e-6 < top3.edge_stretch * 0.85 ||
            top2_clear_principal) {
          refine_limit = 2;
        }
      }
    }
    for (size_t coarse_idx = 0; coarse_idx < refine_limit; ++coarse_idx) {
      const auto &coarse = coarse_candidates[coarse_idx];
      const PbctopoCutCandidate *verify_incumbent =
          (best_safe_candidate.valid && !best_safe_candidate.ultra_clash)
              ? &best_safe_candidate
              : nullptr;
      PbctopoCutCandidate verified_candidate = evaluate_virtual_box_candidate(
          base_atoms, local_edges, local_edge_set, translation_cache,
          this->m_mol[frame_idx].Box,
          coarse.cut_x, coarse.cut_y, coarse.cut_z, coarse.gap_score,
          verify_points,
          runtime_cutoffs.clash_cutoff, &trajectory_context.continuity_pairs,
          &trajectory_context.continuity_anchor_indices,
          &local_hard_required_owner_pairs, runtime_cutoffs.continuity_cutoff,
          verify_incumbent, false, this->m_mol[frame_idx].time_of_snap,
          "verified_candidate", &base_anchor_context,
          coarse.virtual_seed_source, &atom_component_index);
      ++timing.refine_candidates;
      if (verified_candidate.valid)
        ++valid_candidates;
      if (!verified_candidate.valid)
        continue;
      const bool virtual_geometry_outlier = component_candidate_is_outlier(
          verified_candidate,
          best_virtual_geometry_reference.valid
              ? &best_virtual_geometry_reference
              : nullptr);
      if (!verified_candidate.ultra_clash) {
        if (virtual_geometry_outlier) {
          ++timing.pbctopo_virtual_geometry_outlier_observations;
          append_pbctopo_log_line(
              "[pbctopo] observe virtual geometry deviation at frame " +
              to_string(static_cast<long long>(
                  this->m_mol[frame_idx].time_of_snap)) +
              ": principal_rg2_max=" +
              to_string(verified_candidate.principal_rg2_max) +
              ", reference_principal_rg2_max=" +
              to_string(best_virtual_geometry_reference.principal_rg2_max) +
              ", kappa2=" + to_string(verified_candidate.shape_kappa2) +
              ", reference_kappa2=" +
              to_string(best_virtual_geometry_reference.shape_kappa2) +
              ", policy=soft_observation");
        }
        if (better_preclash_cut_candidate(
                verified_candidate, best_virtual_geometry_reference)) {
          best_virtual_geometry_reference = verified_candidate;
        }
        ++safe_candidates;
        record_temporal_certified_candidate(verified_candidate,
                                            verify_points);
        record_report_safe_candidate(verified_candidate);
        if (better_cut_candidate(verified_candidate,
                                 best_verified_virtual_candidate)) {
          best_verified_virtual_candidate = verified_candidate;
        }
        if (better_cut_candidate(verified_candidate, best_safe_candidate)) {
          best_safe_candidate = verified_candidate;
          best_safe_points = verify_points;
          frame_report.mode = PbctopoFrameMode::Virtual;
        }
      }
      if (better_cut_candidate(verified_candidate, best_candidate)) {
        best_candidate = verified_candidate;
      }
    }
    timing.refine_s += chrono::duration<double>(clock::now() - stage_started).count();
    if (best_verified_virtual_candidate.valid &&
        !best_verified_virtual_candidate.ultra_clash &&
        best_verified_virtual_candidate.virtual_seed_source ==
            PbctopoVirtualSeedSource::Kinematic) {
      ++timing.pbctopo_kinematic_seed_best;
    }

    // Keep component scoring aligned with the original local graph.
    bool component_search_used = primary_component_graph_used;
    stage_started = clock::now();
    if (can_use_component_logic) {
      ++timing.component_frames;
      ensure_component_context();
      timing.component_prep_s +=
          chrono::duration<double>(clock::now() - stage_started).count();
      size_t effective_component_count = components.size();
      string component_graph_fallthrough_reason;
      const string component_mode_label =
          pbctopo_component_mode_name(component_mode);
      const PbctopoCutCandidate *component_warm_reference =
          (best_safe_candidate.valid && !best_safe_candidate.ultra_clash)
              ? &best_safe_candidate
              : ((warm_state.has_projection_candidate &&
                  warm_state.previous_projection_candidate.valid &&
                  !warm_state.previous_projection_candidate.ultra_clash)
                     ? &warm_state.previous_projection_candidate
                     : nullptr);
      auto accept_component_cut = [&](const PbctopoCutCandidate &candidate,
                                      const vector<Vec3> &candidate_points,
                                      const vector<Int3> *accepted_offsets,
                                      const string &source_label) {
        if (!candidate.valid || candidate.ultra_clash ||
            candidate_points.size() != base_atoms.size()) {
          return false;
        }
        const PbctopoCutCandidate *virtual_cloud_reference =
            best_virtual_geometry_reference.valid &&
                    !best_virtual_geometry_reference.ultra_clash
                ? &best_virtual_geometry_reference
                : nullptr;
        if (component_candidate_is_outlier(candidate,
                                           virtual_cloud_reference)) {
          if (source_label.find("search") != string::npos ||
              source_label.find("beam") != string::npos) {
            ++timing.component_search_outlier_observations;
          } else {
            ++timing.component_graph_outlier_observations;
          }
          append_pbctopo_log_line(
              "[pbctopo] observe component geometry deviation at frame " +
              to_string(static_cast<long long>(
                  this->m_mol[frame_idx].time_of_snap)) +
              " via " + source_label + ": principal_rg2_max=" +
              to_string(candidate.principal_rg2_max) +
              ", virtual_principal_rg2_max=" +
              to_string(virtual_cloud_reference->principal_rg2_max) +
              ", kappa2=" + to_string(candidate.shape_kappa2) +
              ", virtual_kappa2=" +
              to_string(virtual_cloud_reference->shape_kappa2) +
              ", asphericity=" +
              to_string(candidate.shape_asphericity) +
              ", virtual_asphericity=" +
              to_string(virtual_cloud_reference->shape_asphericity) +
              ", policy=soft_observation");
        }
        PbctopoOwnerPairSplitDetail split_detail;
        size_t split_owner_pairs = 0;
        double max_owner_raw_minus_minimg_delta = 0.0;
        if (candidate_has_split_contacting_owners(
                candidate_points, base_atoms, this->m_mol[frame_idx].Box,
                local_hard_required_owner_pairs, split_owner_pairs,
                max_owner_raw_minus_minimg_delta, &split_detail)) {
          record_frame_failure(PbctopoFailureReason::SplitOwnerReject);
          if (split_detail.valid) {
            record_split_rescue_seed(candidate, candidate_points, split_detail,
                                     source_label);
            append_pbctopo_log_line(
                "[pbctopo] reject split_contacting_owner_pair at frame " +
                to_string(static_cast<long long>(
                    this->m_mol[frame_idx].time_of_snap)) +
                " via " + source_label + ": owners=" +
                to_string(split_detail.owner_a) + "-" +
                to_string(split_detail.owner_b) +
                ", raw_delta_scaled=" +
                format_vec3(split_detail.raw_delta_scaled) +
                ", min_image_delta_scaled=" +
                format_vec3(split_detail.minimg_delta_scaled));
          }
          return false;
        }
        bool accepted = false;
        bool better_than_best = false;
        if (!best_safe_candidate.valid ||
            better_cut_candidate(candidate, best_safe_candidate)) {
          best_safe_candidate = candidate;
          if (accepted_offsets != nullptr &&
              accepted_offsets->size() == components.size()) {
            best_safe_candidate.assignment = make_pbctopo_image_assignment(
                *accepted_offsets, components);
          }
          best_safe_points = candidate_points;
          accepted = true;
          better_than_best = true;
        } else if (candidate_close_enough_to_skip_search(candidate,
                                                         best_safe_candidate)) {
          accepted = true;
        }
        if (!accepted)
          return false;
        record_temporal_certified_candidate(candidate, candidate_points,
                                            accepted_offsets);
        record_report_safe_candidate(candidate);
        if (source_label.find("beam") != string::npos) {
          frame_report.mode = PbctopoFrameMode::ComponentSearchBeam;
        } else if (source_label.find("search") != string::npos) {
          frame_report.mode = PbctopoFrameMode::ComponentSearchExact;
        } else if (source_label.find("warm") != string::npos) {
          frame_report.mode = PbctopoFrameMode::Warm;
        } else {
          frame_report.mode = PbctopoFrameMode::ComponentGraph;
        }
        component_search_used = true;
        if (accepted_offsets != nullptr &&
            accepted_offsets->size() == components.size()) {
          warm_state.previous_component_offsets = *accepted_offsets;
        } else {
          warm_state.previous_component_offsets =
              extract_offsets_and_time(component_atoms, candidate_points);
        }
        if (better_than_best) {
          append_pbctopo_log_line(
              "[pbctopo] " + source_label + ": mode=" + component_mode_label +
              ", raw_components=" + to_string(raw_component_count) +
              ", effective_components=" +
              to_string(effective_component_count) + ", min_pair_d2=" +
              to_string(candidate.min_pair_d2) + ", kappa2=" +
              to_string(candidate.shape_kappa2) + ", asph=" +
              to_string(candidate.shape_asphericity) + ", bbox=" +
              to_string(candidate.bbox_volume));
        }
        return true;
      };
      auto accept_component_candidate =
          [&](const PbctopoComponentCandidate &component_candidate,
              const vector<Vec3> &candidate_points,
              const vector<Int3> *accepted_offsets,
              const string &source_label) {
            if (!component_candidate.valid || component_candidate.ultra_clash ||
                candidate_points.size() != base_atoms.size()) {
              return false;
            }
            PointStats component_score_stats;
            for (const auto &point : candidate_points)
              component_score_stats.add(point.x, point.y, point.z);
            double continuity_path = numeric_limits<double>::infinity();
            double continuity_max_d2 = numeric_limits<double>::infinity();
            if (!continuity_metrics_ok(candidate_points,
                                       trajectory_context.continuity_pairs,
                                       runtime_cutoffs.continuity_cutoff,
                                       continuity_path, continuity_max_d2,
                                       &best_safe_candidate)) {
              record_frame_failure(PbctopoFailureReason::ContinuityReject);
              return false;
            }
            PbctopoCutCandidate component_as_cut =
                score_cut_candidate(component_scoring_edges, candidate_points,
                                    1.5, 1.5, 1.5, 0.0, component_score_stats,
                                    component_candidate.shift_norm2);
            component_as_cut.min_pair_d2 = component_candidate.min_pair_d2;
            component_as_cut.continuity_path = continuity_path;
            component_as_cut.continuity_max_d2 = continuity_max_d2;
            component_as_cut.anchor_mst2 = component_candidate.anchor_mst2;
            component_as_cut.anchor_max_d2 = component_candidate.anchor_max_d2;
            component_as_cut.anchor_graph = component_candidate.anchor_graph;
            component_as_cut.clash_atom_a = component_candidate.clash_atom_a;
            component_as_cut.clash_atom_b = component_candidate.clash_atom_b;
            component_as_cut.ultra_clash = component_candidate.ultra_clash;
            component_as_cut.valid = component_candidate.valid;
            return accept_component_cut(component_as_cut, candidate_points,
                                        accepted_offsets, source_label);
          };
      if (same_component_warm_layout) {
        vector<Int3> warm_offsets = warm_state.previous_component_offsets;
        rebase_component_offsets(warm_offsets, this->m_mol[frame_idx].Box);
        PbctopoCutCandidate warm_cut = improve_component_offsets_locally(
            component_atoms, components, this->m_mol[frame_idx].Box, shifts,
            component_scoring_edges, local_edge_set,
            trajectory_context.continuity_pairs,
            trajectory_context.continuity_anchor_indices, warm_offsets,
            runtime_cutoffs.clash_cutoff,
            runtime_cutoffs.continuity_cutoff, component_points,
            make_component_refine_callbacks(), anchor_backend, &timing,
            &warm_state, component_warm_reference, &base_anchor_context);
        PbctopoOwnerPairSplitDetail warm_cut_split_detail;
        size_t warm_cut_split_owner_pairs = 0;
        double warm_cut_split_gap = 0.0;
        const bool warm_cut_split =
            warm_cut.valid && !warm_cut.ultra_clash &&
            candidate_has_split_contacting_owners(
                component_points, base_atoms, this->m_mol[frame_idx].Box,
                local_hard_required_owner_pairs, warm_cut_split_owner_pairs,
                warm_cut_split_gap, &warm_cut_split_detail);
        if (warm_cut_split && warm_cut_split_detail.valid) {
          record_frame_failure(PbctopoFailureReason::SplitOwnerReject);
          record_split_rescue_seed(warm_cut, component_points,
                                   warm_cut_split_detail,
                                   "warm_component_local_improve");
          append_pbctopo_log_line(
              "[pbctopo] reject split_contacting_owner_pair at frame " +
              to_string(
                  static_cast<long long>(this->m_mol[frame_idx].time_of_snap)) +
              " via warm_component_local_improve: owners=" +
              to_string(warm_cut_split_detail.owner_a) + "-" +
              to_string(warm_cut_split_detail.owner_b) +
              ", raw_delta_scaled=" +
              format_vec3(warm_cut_split_detail.raw_delta_scaled) +
              ", min_image_delta_scaled=" +
              format_vec3(warm_cut_split_detail.minimg_delta_scaled));
        }
        if (!warm_cut_split && warm_cut.valid && !warm_cut.ultra_clash) {
          bool warm_accepted = false;
          if (!best_safe_candidate.valid ||
              better_cut_candidate(warm_cut, best_safe_candidate)) {
            best_safe_candidate = warm_cut;
            best_safe_points = component_points;
            warm_accepted = true;
            frame_report.mode = PbctopoFrameMode::Warm;
            record_temporal_certified_candidate(warm_cut, component_points,
                                                &warm_offsets);
            record_report_safe_candidate(warm_cut);
          } else if (candidate_close_enough_to_skip_search(
                         warm_cut, best_safe_candidate)) {
            warm_accepted = true;
          }
          if (warm_accepted) {
            ++timing.component_warm_used;
            warm_state.previous_component_offsets = warm_offsets;
          }
        }
      }
      if (component_mode == PbctopoComponentMode::Exact) {
        double component_shift_norm2 = 0.0;
        PbctopoCutCandidate graph_as_cut;
        vector<Int3> graph_tree_offsets;
        vector<Int3> graph_refined_offsets;
        vector<PbctopoGraphEdge> graph_tree_edges;
        PbctopoCycleSyncMetrics graph_cycle_metrics;
        if (!component_search_used) {
          const auto component_graph_started = clock::now();
          const bool graph_ok = solve_component_offsets_via_pair_graph(
              component_atoms, this->m_mol[frame_idx].Box, components, shifts,
              component_scoring_edges, local_edge_set, component_points,
              trajectory_context.continuity_pairs,
              trajectory_context.continuity_anchor_indices, runtime_cutoffs,
              graph_as_cut, component_shift_norm2, beam_width, &timing,
              make_component_refine_callbacks(), anchor_backend,
              same_component_warm_layout ? &warm_state.previous_component_offsets
                                         : nullptr,
              &warm_state, component_warm_reference, &graph_tree_offsets,
              &graph_refined_offsets, &graph_tree_edges, &graph_cycle_metrics,
              &base_anchor_context);
          timing.component_graph_s +=
              chrono::duration<double>(clock::now() - component_graph_started)
                  .count();
          frame_report.cycle_support = graph_cycle_metrics.cycle_support;
          frame_report.cycle_mismatch = graph_cycle_metrics.cycle_mismatch;
          frame_report.inconsistent_edge_count =
              graph_cycle_metrics.inconsistent_edge_count;
          frame_report.fundamental_cycle_count =
              graph_cycle_metrics.fundamental_cycle_count;
          frame_report.fundamental_cycle_residual_rms =
              graph_cycle_metrics.fundamental_cycle_residual_rms;
          frame_report.fundamental_cycle_residual_max =
              graph_cycle_metrics.fundamental_cycle_residual_max;
          graph_as_cut.hard_feasibility.cycle_mismatches =
              graph_cycle_metrics.inconsistent_edge_count;
          if (!graph_ok) {
            ++timing.component_graph_failed;
            component_graph_fallthrough_reason = "graph_fail";
            record_frame_failure(PbctopoFailureReason::ComponentGraphFailed);
            clear_component_trial_hints(warm_state);
          }
          if (graph_ok && component_points.size() != base_atoms.size()) {
            ++timing.component_graph_points_mismatch;
            component_graph_fallthrough_reason = "points_mismatch";
          }
          if (graph_ok && component_points.size() == base_atoms.size()) {
            if (accept_component_cut(graph_as_cut, component_points,
                                     &graph_refined_offsets,
                                     "component graph solve")) {
              ++timing.component_graph_used;
            } else if (!graph_as_cut.valid) {
              ++timing.component_graph_invalid;
              component_graph_fallthrough_reason = "invalid";
              record_frame_failure(PbctopoFailureReason::ProjectionInvalid);
            } else if (graph_as_cut.ultra_clash) {
              ++timing.component_graph_ultraclash;
              component_graph_fallthrough_reason = "ultra_clash";
              record_frame_failure(PbctopoFailureReason::FullVerifyClash);
            } else {
              ++timing.component_graph_not_better_not_close;
              record_frame_failure(PbctopoFailureReason::NotBetterThanIncumbent);
              const char *close_reject_reason =
                  candidate_skip_search_reject_reason(graph_as_cut,
                                                      best_safe_candidate);
              component_graph_fallthrough_reason =
                  string("not_better_not_close/") +
                  (close_reject_reason != nullptr ? close_reject_reason
                                                 : "unknown");
            }
          }
        }

        if (!component_search_used) {
          ++timing.component_search_frames;
          auto component_clash_atom_summary =
              [&](size_t atom_idx) -> string {
            if (atom_idx >= component_atoms.size())
              return "na";
            size_t component_idx = numeric_limits<size_t>::max();
            for (size_t idx = 0; idx < components.size(); ++idx) {
              if (binary_search(components[idx].begin(), components[idx].end(),
                                atom_idx)) {
                component_idx = idx;
                break;
              }
            }
            const PbctopoAtom &atom = component_atoms[atom_idx];
            return to_string(atom_idx) + "/global=" +
                   to_string(atom.global_atom_id) + "/owner=" +
                   to_string(atom.owner) + "/component=" +
                   (component_idx == numeric_limits<size_t>::max()
                        ? string("na")
                        : to_string(component_idx)) +
                   "/heavy=" + (atom.is_heavy ? "yes" : "no");
          };
          append_pbctopo_log_line(
              "[pbctopo] component graph fell through to search at frame " +
              to_string(static_cast<long long>(
                  this->m_mol[frame_idx].time_of_snap)) +
              ": mode=" + component_mode_label +
              ", raw_components=" + to_string(raw_component_count) +
              ", effective_components=" +
              to_string(effective_component_count) + ", reason=" +
              component_graph_fallthrough_reason + ", anchor_count=" +
              to_string(trajectory_context.continuity_anchor_indices.size()) +
              ", anchor_policy=full" +
              ", graph_valid=" + string(graph_as_cut.valid ? "yes" : "no") +
              ", graph_ultra_clash=" +
              string(graph_as_cut.ultra_clash ? "yes" : "no") +
              ", graph_anchor_max_d2=" + to_string(graph_as_cut.anchor_max_d2) +
              ", ref_anchor_max_d2=" + to_string(best_safe_candidate.anchor_max_d2) +
              ", graph_anchor_mst2=" + to_string(graph_as_cut.anchor_mst2) +
              ", ref_anchor_mst2=" + to_string(best_safe_candidate.anchor_mst2) +
              ", graph_min_pair_d2=" + to_string(graph_as_cut.min_pair_d2) +
              ", graph_clash_atom_a=" +
              component_clash_atom_summary(graph_as_cut.clash_atom_a) +
              ", graph_clash_atom_b=" +
              component_clash_atom_summary(graph_as_cut.clash_atom_b) +
              ", graph_kappa2=" + to_string(graph_as_cut.shape_kappa2) +
              ", graph_asph=" + to_string(graph_as_cut.shape_asphericity) +
              ", graph_bbox=" + to_string(graph_as_cut.bbox_volume));
          append_pbctopo_log_line(
              "[pbctopo] component graph offsets at frame " +
              to_string(static_cast<long long>(
                  this->m_mol[frame_idx].time_of_snap)) +
              ": tree_offsets=" + format_component_offsets(graph_tree_offsets) +
              ", refined_offsets=" + format_component_offsets(graph_refined_offsets) +
              ", tree_edges=" + format_graph_edges(graph_tree_edges));
          const auto component_search_started = clock::now();
          const size_t exact_state_upper_bound =
              pbctopo_component_assignment_upper_bound(
                  components.size(), shifts.size(),
                  kPbctopoExactAssignmentBudget + 1);
          timing.component_search_exact_state_upper_bound =
              max(timing.component_search_exact_state_upper_bound,
                  exact_state_upper_bound);
          const bool use_budgeted_beam =
              exact_state_upper_bound > kPbctopoExactAssignmentBudget;
          if (use_budgeted_beam) {
            ++timing.component_search_exact_budget_skips;
            append_pbctopo_log_line(
                "[pbctopo] replace exact search with beam at frame " +
                to_string(static_cast<long long>(
                    this->m_mol[frame_idx].time_of_snap)) +
                ": assignment_upper_bound=" +
                to_string(exact_state_upper_bound) +
                ", budget=" + to_string(kPbctopoExactAssignmentBudget));
          }
          const auto search_incumbent =
              (best_safe_candidate.valid && !best_safe_candidate.ultra_clash)
                  ? &best_safe_candidate
                  : nullptr;
          const auto search_incumbent_points =
              (search_incumbent != nullptr &&
               best_safe_points.size() == base_atoms.size())
                  ? &best_safe_points
                  : nullptr;
          PbctopoComponentCandidate component_candidate =
              use_budgeted_beam
                  ? search_component_shift_candidate_beam(
                        component_atoms, this->m_mol[frame_idx].Box,
                        components, shifts,
                        trajectory_context.continuity_pairs,
                        trajectory_context.continuity_anchor_indices,
                        runtime_cutoffs, component_points, beam_width,
                        make_component_search_callbacks(),
                        same_component_warm_layout
                            ? &warm_state.previous_component_offsets
                            : nullptr,
                        search_incumbent, search_incumbent_points, &timing,
                        anchor_backend, &base_anchor_context)
                  : search_component_shift_candidate(
                        component_atoms, this->m_mol[frame_idx].Box,
                        components, shifts,
                        trajectory_context.continuity_pairs,
                        trajectory_context.continuity_anchor_indices,
                        runtime_cutoffs, component_points,
                        make_component_search_callbacks(),
                        same_component_warm_layout
                            ? &warm_state.previous_component_offsets
                            : nullptr,
                        search_incumbent, search_incumbent_points, &timing,
                        anchor_backend, &base_anchor_context);
          timing.component_search_s +=
              chrono::duration<double>(clock::now() - component_search_started)
                  .count();
          if (accept_component_candidate(component_candidate, component_points,
                                         nullptr,
                                         use_budgeted_beam
                                             ? "component beam budget winner"
                                             : "component search winner")) {
            const vector<Int3> search_offsets =
                extract_offsets_and_time(component_atoms, component_points);
            warm_state.previous_component_offsets = search_offsets;
            append_pbctopo_log_line(
                "[pbctopo] component search winner offsets at frame " +
                to_string(static_cast<long long>(
                    this->m_mol[frame_idx].time_of_snap)) +
                ": offsets=" + format_component_offsets(search_offsets));
          }
          if (component_search_used || !best_safe_candidate.valid) {
            append_pbctopo_log_line(
                "[pbctopo] component-shift search: mode=" + component_mode_label +
                ", raw_components=" + to_string(raw_component_count) +
                ", effective_components=" +
                to_string(effective_component_count) + ", components=" +
                to_string(components.size()) + ", ultra_clash=" +
                string(component_candidate.ultra_clash ? "yes" : "no") +
                ", min_pair_d2=" + to_string(component_candidate.min_pair_d2) +
                ", kappa2=" + to_string(component_candidate.shape_kappa2) +
                ", asph=" + to_string(component_candidate.shape_asphericity) +
                ", bbox=" + to_string(component_candidate.bbox_volume));
          }
        }
      } else if (component_mode == PbctopoComponentMode::MergeThenGraph) {
        ++timing.component_large_mode_merge_frames;
        const PbctopoMergedComponentPlan merge_plan =
            build_chain_merged_component_plan(
                component_atoms, this->m_mol[frame_idx].Box, components);
        effective_component_count = merge_plan.merged_components.size();
        if (raw_component_count > effective_component_count)
          timing.component_large_mode_merge_count +=
              raw_component_count - effective_component_count;
        const size_t merged_exact_state_upper_bound =
            pbctopo_component_assignment_upper_bound(
                effective_component_count, shifts.size(),
                kPbctopoExactAssignmentBudget + 1);
        timing.component_search_exact_state_upper_bound =
            max(timing.component_search_exact_state_upper_bound,
                merged_exact_state_upper_bound);
        const bool merged_exact_within_budget =
            merged_exact_state_upper_bound <= kPbctopoExactAssignmentBudget;
        const bool allow_exact_search_on_merged =
            effective_component_count <= exact_component_limit &&
            merged_exact_within_budget;
        if (effective_component_count <= exact_component_limit &&
            !merged_exact_within_budget) {
          ++timing.component_search_exact_budget_skips;
          append_pbctopo_log_line(
              "[pbctopo] skip merged exact search at frame " +
              to_string(static_cast<long long>(
                  this->m_mol[frame_idx].time_of_snap)) +
              ": assignment_upper_bound=" +
              to_string(merged_exact_state_upper_bound) +
              ", budget=" + to_string(kPbctopoExactAssignmentBudget) +
              ", policy=component_graph_refine");
        }
        if (!allow_exact_search_on_merged)
          ++timing.component_large_mode_pair_only_frames;

        PbctopoCutCandidate merged_mode_cut;
        vector<Int3> merged_mode_raw_offsets;
        double component_shift_norm2 = 0.0;
        const auto merge_mode_started = clock::now();
        const bool merge_ok = solve_large_component_merge_mode(
            component_atoms, this->m_mol[frame_idx].Box, components, merge_plan,
            shifts, component_scoring_edges, local_edge_set, component_points,
            trajectory_context.continuity_pairs,
            trajectory_context.continuity_anchor_indices, runtime_cutoffs,
            merged_mode_cut, component_shift_norm2, &timing, anchor_backend,
            same_component_warm_layout ? &warm_state.previous_component_offsets
                                       : nullptr,
            &warm_state, component_warm_reference, allow_exact_search_on_merged,
            &merged_mode_raw_offsets, &base_anchor_context);
        timing.component_graph_s +=
            chrono::duration<double>(clock::now() - merge_mode_started).count();
        if (merge_ok &&
            accept_component_cut(merged_mode_cut, component_points,
                                 &merged_mode_raw_offsets,
                                 "component merge-graph solve")) {
          ++timing.component_graph_used;
        } else {
          ++timing.component_large_mode_fallback_frames;
          record_frame_failure(!merge_ok
                                   ? PbctopoFailureReason::ComponentGraphFailed
                                   : PbctopoFailureReason::NotBetterThanIncumbent);
          const char *merge_reject_reason =
              merge_ok ? candidate_skip_search_reject_reason(
                             merged_mode_cut, best_safe_candidate)
                       : nullptr;
          component_graph_fallthrough_reason = !merge_ok
                                                   ? "graph_invalid"
                                                   : string("not_better_not_close/") +
                                                          (merge_reject_reason != nullptr
                                                               ? merge_reject_reason
                                                               : "unknown");
          append_pbctopo_log_line(
              "[pbctopo] large-component merge mode fell back at frame " +
              to_string(static_cast<long long>(
                  this->m_mol[frame_idx].time_of_snap)) +
              ": mode=" + component_mode_label +
              ", raw_components=" + to_string(raw_component_count) +
              ", effective_components=" +
              to_string(effective_component_count) + ", reason=" +
              component_graph_fallthrough_reason);
        }
      } else if (component_mode == PbctopoComponentMode::Beam) {
        ++timing.component_large_mode_beam_frames;
        effective_component_count = components.size();
        const auto beam_started = clock::now();
        PbctopoComponentCandidate beam_candidate =
            search_component_shift_candidate_beam(
                component_atoms, this->m_mol[frame_idx].Box, components, shifts,
                trajectory_context.continuity_pairs,
                trajectory_context.continuity_anchor_indices, runtime_cutoffs,
                component_points, beam_width, make_component_search_callbacks(),
                same_component_warm_layout
                    ? &warm_state.previous_component_offsets
                    : nullptr,
                (best_safe_candidate.valid && !best_safe_candidate.ultra_clash)
                    ? &best_safe_candidate
                    : nullptr,
                (best_safe_candidate.valid && !best_safe_candidate.ultra_clash &&
                 best_safe_points.size() == base_atoms.size())
                    ? &best_safe_points
                    : nullptr,
                &timing, anchor_backend, &base_anchor_context);
        timing.component_search_s +=
            chrono::duration<double>(clock::now() - beam_started).count();
        if (!accept_component_candidate(beam_candidate, component_points, nullptr,
                                        "component beam solve")) {
          ++timing.component_large_mode_fallback_frames;
          record_frame_failure(beam_candidate.valid && !beam_candidate.ultra_clash
                                   ? PbctopoFailureReason::NotBetterThanIncumbent
                                   : PbctopoFailureReason::BeamEmpty);
          PbctopoCutCandidate beam_as_cut =
              component_candidate_to_cut_candidate(beam_candidate);
          const char *beam_reject_reason =
              (beam_candidate.valid && !beam_candidate.ultra_clash)
                  ? candidate_skip_search_reject_reason(beam_as_cut,
                                                        best_safe_candidate)
                  : nullptr;
          component_graph_fallthrough_reason =
              beam_candidate.valid && !beam_candidate.ultra_clash
                  ? string("not_better_not_close/") +
                        (beam_reject_reason != nullptr
                             ? beam_reject_reason
                             : "unknown")
                  : "beam_no_complete_assignment";
          append_pbctopo_log_line(
              "[pbctopo] large-component beam mode fell back at frame " +
              to_string(static_cast<long long>(
                  this->m_mol[frame_idx].time_of_snap)) +
              ": mode=" + component_mode_label +
              ", raw_components=" + to_string(raw_component_count) +
              ", effective_components=" +
              to_string(effective_component_count) + ", reason=" +
              component_graph_fallthrough_reason + ", min_pair_d2=" +
              to_string(beam_candidate.min_pair_d2) + ", kappa2=" +
              to_string(beam_candidate.shape_kappa2) + ", asph=" +
              to_string(beam_candidate.shape_asphericity) + ", bbox=" +
              to_string(beam_candidate.bbox_volume));
        }
      }
    } else {
      timing.component_prep_s +=
          chrono::duration<double>(clock::now() - stage_started).count();
    }
    warm_state.previous_component_sizes = component_sizes;
    warm_state.has_component_offsets =
        warm_state.previous_component_offsets.size() == components.size() &&
        !components.empty();
    if (!warm_state.has_component_offsets)
      clear_component_trial_hints(warm_state);

    if (best_safe_candidate.valid && !best_safe_candidate.ultra_clash &&
        best_safe_points.size() == base_atoms.size()) {
      if (certify_candidate_for_output(best_safe_candidate,
                                       best_safe_points)) {
        record_temporal_certified_candidate(best_safe_candidate,
                                            best_safe_points);
        projection_clash_skip_active =
            candidate_has_projection_clash_margin(best_safe_candidate,
                                                  runtime_cutoffs.clash_cutoff);
        warm_state.previous_projection_candidate = best_safe_candidate;
        warm_state.has_projection_candidate = true;
        update_kinematic_history_for_frame();
        apply_cut_candidate_timed(best_safe_points);
        finalize_frame_report(PbctopoFrameStatus::Certified,
                              frame_report.mode == PbctopoFrameMode::None
                                  ? PbctopoFrameMode::Virtual
                                  : frame_report.mode,
                              PbctopoFailureReason::None,
                              &best_safe_candidate);
        continue;
      }
      ++timing.pbctopo_diagnostic_only_candidates;
      best_safe_certification_rejected = true;
      record_frame_failure(PbctopoFailureReason::CertificationReject);
      append_pbctopo_log_line(
          "[pbctopo] reject output candidate at frame " +
          to_string(static_cast<long long>(
              this->m_mol[frame_idx].time_of_snap)) +
          ": assignment_available=" +
          string(best_safe_candidate.hard_feasibility.assignment_available
                     ? "yes"
                     : "no") +
          ", hard_edge_residuals=" +
          to_string(best_safe_candidate.hard_feasibility.hard_edge_residuals) +
          ", cycle_mismatches=" +
          to_string(best_safe_candidate.hard_feasibility.cycle_mismatches) +
          ", hard_clashes=" +
          to_string(best_safe_candidate.hard_feasibility.hard_clashes));
    }

    if ((!best_safe_candidate.valid || best_safe_candidate.ultra_clash ||
         best_safe_certification_rejected) &&
        split_rescue_available) {
      ensure_component_context();
      ++timing.hull_rescue_attempts;
      PbctopoCutCandidate best_hull_candidate;
      vector<Vec3> best_hull_points;
      best_hull_points.reserve(base_atoms.size());
      Int3 best_hull_shift{};
      int best_hull_owner = -1;
      double best_hull_volume = numeric_limits<double>::infinity();
      size_t hull_valid_candidates = 0;
      vector<int> rescue_owner_candidates;
      rescue_owner_candidates.reserve(2);
      rescue_owner_candidates.push_back(split_rescue_detail.owner_b);
      if (split_rescue_detail.owner_a != split_rescue_detail.owner_b)
        rescue_owner_candidates.push_back(split_rescue_detail.owner_a);
      for (int rescue_owner : rescue_owner_candidates) {
        build_owner_spans(base_atoms, rescue_owner, hull_owner_spans);
        if (hull_owner_spans.empty())
          continue;
        for (const auto &shift : shifts) {
          if (shift.ix == 0 && shift.iy == 0 && shift.iz == 0)
            continue;
          hull_rescue_points.clear();
          pbctopo_hull::shift_points_by_spans(
              split_rescue_points, hull_owner_spans,
              Vec3{shift.tx, shift.ty, shift.tz}, hull_rescue_points,
              &hull_scratch);
          PbctopoCutCandidate hull_candidate;
          const PbctopoCutCandidate *hull_incumbent =
              best_hull_candidate.valid ? &best_hull_candidate : nullptr;
          if (!validate_rescue_candidate_points(hull_rescue_points,
                                                hull_incumbent,
                                                hull_candidate)) {
            continue;
          }
          double hull_volume = numeric_limits<double>::infinity();
          if (!pbctopo_hull::quickhull_volume(hull_rescue_points, hull_volume,
                                              &hull_scratch)) {
            continue;
          }
          ++hull_valid_candidates;
          const bool better_hull =
              !best_hull_candidate.valid ||
              hull_volume + 1.0e-6 < best_hull_volume ||
              (fabs(hull_volume - best_hull_volume) <= 1.0e-6 &&
               better_cut_candidate(hull_candidate, best_hull_candidate));
          if (!better_hull)
            continue;
          best_hull_candidate = hull_candidate;
          best_hull_volume = hull_volume;
          best_hull_shift = {shift.ix, shift.iy, shift.iz};
          best_hull_owner = rescue_owner;
          best_hull_points = hull_rescue_points;
        }
      }
      if (best_hull_candidate.valid) {
        append_pbctopo_log_line(
            "[pbctopo] hull rescue candidate at frame " +
            to_string(static_cast<long long>(
                this->m_mol[frame_idx].time_of_snap)) +
            " from " + split_rescue_source + ": moved_owner=" +
            to_string(best_hull_owner) + ", applied_shift=" +
            format_int3(best_hull_shift) + ", hull_volume=" +
            to_string(best_hull_volume) + ", valid_candidates=" +
            to_string(hull_valid_candidates) + ", avx2=" +
            string(pbctopo_hull::avx2_available() ? "on" : "off"));
        if (accept_rescue_candidate(best_hull_points, best_hull_candidate)) {
          ++timing.hull_rescue_successes;
          finalize_frame_report(PbctopoFrameStatus::Rescued,
                                PbctopoFrameMode::HullRescue,
                                PbctopoFailureReason::None,
                                &best_safe_candidate);
          continue;
        }
      }

      ++timing.distance_rescue_attempts;
      vector<Vec3> rescued_points;
      rescued_points.reserve(base_atoms.size());
      Int3 rescue_shift{};
      double rescue_raw_delta_norm = numeric_limits<double>::infinity();
      if (rescue_split_contacting_owner_candidate(
              split_rescue_points, base_atoms, this->m_mol[frame_idx].Box,
              local_hard_required_owner_pairs, shifts, split_rescue_detail,
              rescued_points, rescue_shift, rescue_raw_delta_norm)) {
        PbctopoCutCandidate rescued_candidate;
        if (validate_rescue_candidate_points(rescued_points, nullptr,
                                             rescued_candidate)) {
          append_pbctopo_log_line(
              "[pbctopo] distance rescue candidate at frame " +
              to_string(static_cast<long long>(
                  this->m_mol[frame_idx].time_of_snap)) +
              " from " + split_rescue_source + ": applied_shift=" +
              format_int3(rescue_shift) + ", rescued_raw_delta_norm=" +
              to_string(rescue_raw_delta_norm));
          if (accept_rescue_candidate(rescued_points, rescued_candidate)) {
            ++timing.distance_rescue_successes;
            finalize_frame_report(PbctopoFrameStatus::Rescued,
                                  PbctopoFrameMode::DistanceRescue,
                                  PbctopoFailureReason::None,
                                  &best_safe_candidate);
            continue;
          }
        }
      }
    }

    if (best_candidate.valid)
      ++timing.pbctopo_diagnostic_only_candidates;
    const PbctopoFailureReason fallback_reason =
        frame_report.fallback_reason == PbctopoFailureReason::CertificationReject
            ? PbctopoFailureReason::CertificationReject
            : PbctopoFailureReason::FallbackPbci;
    const string fallback_message =
        "[pbctopo] frame " +
        to_string(static_cast<long long>(this->m_mol[frame_idx].time_of_snap)) +
        " fell back to pbci=0 after no certified image assignment: "
        "safe_candidates=" +
        to_string(safe_candidates) + ", best.min_pair_d2=" +
        to_string(best_candidate.min_pair_d2) + ", best.kappa2=" +
        to_string(best_candidate.shape_kappa2) + ", best.asph=" +
        to_string(best_candidate.shape_asphericity) + ", best.bbox=" +
        to_string(best_candidate.bbox_volume);
    append_pbctopo_log_line(fallback_message);
    TitanUI::print_warning(fallback_message);
    ++fallback_frames;
    projection_clash_skip_active = false;
    clear_component_trial_hints(warm_state);
    reset_kinematic_seed_history(warm_state, timing);
    fallback_pbci_frame(frame_idx);
    finalize_frame_report(PbctopoFrameStatus::Fallback,
                          PbctopoFrameMode::Pbci, fallback_reason,
                          best_candidate.valid ? &best_candidate : nullptr);
  }

  if (temporal_mode == PbctopoTemporalMode::Viterbi) {
    vector<vector<PbctopoTemporalCandidate>> temporal_frames(frame_count);
    for (size_t frame_idx = 0; frame_idx < frame_count; ++frame_idx) {
      const auto &cuts = temporal_cut_candidates[frame_idx];
      const auto &greedy = greedy_selected_candidates[frame_idx];
      if (cuts.empty() || !greedy.valid || !greedy.assignment.valid ||
          !greedy.hard_feasibility.certified())
        continue;
      ++timing.pbctopo_temporal_candidate_frames;
      timing.pbctopo_temporal_candidates += cuts.size();
      temporal_frames[frame_idx].reserve(cuts.size());
      const double denominator =
          static_cast<double>(max<size_t>(1, cuts.size() - 1));
      for (size_t rank = 0; rank < cuts.size(); ++rank) {
        PbctopoTemporalCandidate candidate;
        candidate.center = {cuts[rank].cut_x, cuts[rank].cut_y,
                            cuts[rank].cut_z};
        candidate.emission_cost = static_cast<double>(rank) / denominator;
        candidate.component_offsets =
            cuts[rank].assignment.component_offsets;
        candidate.assignment_signature = cuts[rank].assignment.signature;
        candidate.certified = cuts[rank].hard_feasibility.certified() &&
                              cuts[rank].assignment.valid;
        temporal_frames[frame_idx].push_back(std::move(candidate));
      }
    }

    const auto temporal_result = titan_pbctopo::solve_pbctopo_temporal_viterbi(
        temporal_frames, Arguments.pbc.pbctopo_temporal_penalty);
    ofstream temporal_report("pbctopo_temporal_path.csv", ios::trunc);
    if (temporal_report.is_open())
      temporal_report << "frame,time_ps,selected_rank,cut_x,cut_y,cut_z,"
                         "emission_cost,assignment_signature,"
                         "changed_from_greedy,replay_applied\n";
    const size_t invalid = numeric_limits<size_t>::max();
    if (temporal_result.valid) {
      for (size_t frame_idx = 0; frame_idx < frame_count; ++frame_idx) {
        const size_t selected = temporal_result.selected[frame_idx];
        if (selected == invalid ||
            selected >= temporal_cut_candidates[frame_idx].size())
          continue;
        const auto &chosen = temporal_cut_candidates[frame_idx][selected];
        const auto &greedy = greedy_selected_candidates[frame_idx];
        const bool changed = !same_pbctopo_image_assignment(
            chosen.assignment, greedy.assignment);
        if (changed)
          ++timing.pbctopo_temporal_proposed_changes;
        ++timing.pbctopo_temporal_viterbi_frames;
        ++timing.pbctopo_temporal_replay_frames;

        bool replay_applied = !changed;
        if (changed) {
          const size_t layout_idx = temporal_frame_layout[frame_idx];
          vector<PbctopoAtom> replay_atoms;
          replay_atoms.reserve(trajectory_context.total_atoms_hint);
          append_atoms_for_system(this->m_mol[frame_idx], &this->info,
                                  &trajectory_context.rec_static, 0,
                                  replay_atoms);
          if (frame_idx < mol2.m_mol.size() &&
              !mol2.m_mol[frame_idx].pos.x.empty()) {
            append_atoms_for_system(mol2.m_mol[frame_idx], &mol2.info,
                                    &trajectory_context.lig_static, 1,
                                    replay_atoms);
          }
          vector<Vec3> applied_points;
          applied_points.reserve(replay_atoms.size());
          for (const auto &atom : replay_atoms)
            applied_points.push_back({atom.x, atom.y, atom.z});
          vector<Vec3> replayed_points;
          if (layout_idx < temporal_component_layouts.size() &&
              replay_pbctopo_assignment_delta(
                  applied_points,
                  temporal_component_layouts[layout_idx].atom_component_index,
                  greedy.assignment, chosen.assignment,
                  this->m_mol[frame_idx].Box, replayed_points)) {
            apply_pbctopo_points_direct(replay_atoms, replayed_points, *this,
                                        mol2, frame_idx);
            replay_applied = true;
            ++timing.pbctopo_temporal_replay_applied;
            ++timing.pbctopo_temporal_path_changes;
          } else {
            ++timing.pbctopo_temporal_replay_skipped;
          }
        }
        if (deferred_frame_report_valid[frame_idx] != 0) {
          auto &report = deferred_frame_reports[frame_idx];
          report.temporal_selected = true;
          report.temporal_changed_from_greedy = changed && replay_applied;
          report.temporal_assignment_signature = chosen.assignment.signature;
          if (replay_applied) {
            report.assignment_signature = chosen.assignment.signature;
            report.assignment_component_count =
                chosen.assignment.component_offsets.size();
            report.hard_feasibility = chosen.hard_feasibility;
            report.best_score = make_score_vector(chosen);
          }
        }
        if (temporal_report.is_open()) {
          temporal_report << frame_report_offset + frame_idx << ','
                          << this->m_mol[frame_idx].time_of_snap << ','
                          << selected << ',' << chosen.cut_x << ','
                          << chosen.cut_y << ',' << chosen.cut_z << ','
                          << temporal_frames[frame_idx][selected].emission_cost
                          << ',' << chosen.assignment.signature << ','
                          << (changed ? 1 : 0) << ','
                          << (replay_applied ? 1 : 0) << '\n';
        }
      }
      append_pbctopo_log_line(
          "[pbctopo temporal] mode=viterbi, frames=" +
          to_string(timing.pbctopo_temporal_viterbi_frames) +
          ", candidates=" + to_string(timing.pbctopo_temporal_candidates) +
          ", proposed_changes=" +
          to_string(timing.pbctopo_temporal_proposed_changes) +
          ", applied_changes=" +
          to_string(timing.pbctopo_temporal_path_changes) +
          ", replay=immutable_assignment_delta" +
          ", total_cost=" + to_string(temporal_result.total_cost));
    }
    if (pbctopo_frame_report.is_open()) {
      for (size_t frame_idx = 0; frame_idx < frame_count; ++frame_idx) {
        if (deferred_frame_report_valid[frame_idx] != 0) {
          write_pbctopo_frame_report_row(pbctopo_frame_report,
                                         deferred_frame_reports[frame_idx]);
        }
      }
    }
  }

  timing.pbctopo_frame_point_buffer_peak_bytes_estimate =
      max(timing.pbctopo_frame_point_buffer_peak_bytes_estimate,
          scratch.point_capacity_bytes_estimate());
  if (fallback_frames > 0) {
    const string fallback_summary = "[pbctopo] Fell back to pbci=0 on " +
                                    to_string(fallback_frames) + " frame(s).";
    append_pbctopo_log_line(fallback_summary);
    TitanUI::print_warning(fallback_summary);
  }
  append_pbctopo_log_line(
      "[pbctopo timing] frames=" + to_string(timing.frames) +
      ", append=" + to_string(timing.append_s) + " s, heavy_subset=" +
      to_string(timing.heavy_subset_s) + " s, translation_cache=" +
      to_string(timing.translation_cache_s) + " s, graph=" +
      to_string(timing.graph_s) + " s, pbctopo_local_graph_geometry_s=" +
      to_string(timing.pbctopo_local_graph_geometry_s) +
      " s, pbctopo_topology_edge_append_s=" +
      to_string(timing.pbctopo_topology_edge_append_s) +
      " s, projection_edges=" +
      to_string(timing.projection_edge_s) + " s, component_split=" +
      to_string(timing.component_split_s) + " s, seeds=" +
      to_string(timing.seed_s) + " s, event_planes=" +
      to_string(timing.event_planes_s) +
      " s, pbctopo_density_cut_plane_build=" +
      to_string(timing.pbctopo_density_cut_plane_build_s) +
      " s, pbctopo_virtual_seam_sweep=" +
      to_string(timing.pbctopo_virtual_seam_sweep_s) + " s, coarse=" +
      to_string(timing.coarse_s) + " s, refine=" +
      to_string(timing.refine_s) + " s, early_projection=" +
      to_string(timing.early_projection_s) + " s, warm_component=" +
      to_string(timing.warm_component_s) + " s, component_prep=" +
      to_string(timing.component_prep_s) +
      " s, component_context=" + to_string(timing.component_context_s) +
      " s, component_graph=" + to_string(timing.component_graph_s) +
      " s, component_graph_precompute=" +
      to_string(timing.component_graph_precompute_s) +
      " s, component_graph_lazy_materialize=" +
      to_string(timing.component_graph_lazy_materialize_s) +
      " s, component_graph_pair_score=" +
      to_string(timing.component_graph_pair_score_s) +
      " s, component_graph_tree=" +
      to_string(timing.component_graph_tree_s) +
      " s, component_graph_refine=" +
      to_string(timing.component_graph_refine_s) +
      " s, component_graph_refine_rebase=" +
      to_string(timing.component_graph_refine_rebase_s) +
      " s, component_graph_refine_coarse=" +
      to_string(timing.component_graph_refine_coarse_s) +
      " s, component_graph_refine_coarse_stage1=" +
      to_string(timing.component_graph_refine_coarse_stage1_s) +
      " s, component_graph_refine_coarse_stage1_build_points=" +
      to_string(timing.component_graph_refine_coarse_stage1_build_points_s) +
      " s, component_graph_refine_coarse_stage1_continuity=" +
      to_string(timing.component_graph_refine_coarse_stage1_continuity_s) +
      " s, component_graph_refine_coarse_stage1_anchor=" +
      to_string(timing.component_graph_refine_coarse_stage1_anchor_s) +
      " s, component_graph_refine_coarse_stage1_cut_score=" +
      to_string(timing.component_graph_refine_coarse_stage1_cut_score_s) +
      " s, component_graph_refine_coarse_stage1_cut_score_principal=" +
      to_string(timing.component_graph_refine_coarse_stage1_cut_score_principal_s) +
      " s, component_graph_refine_coarse_stage1_cut_score_edge_scan=" +
      to_string(timing.component_graph_refine_coarse_stage1_cut_score_edge_scan_s) +
      " s, component_graph_refine_coarse_stage2=" +
      to_string(timing.component_graph_refine_coarse_stage2_s) +
      " s, component_graph_refine_coarse_stage2_build_points=" +
      to_string(timing.component_graph_refine_coarse_stage2_build_points_s) +
      " s, component_graph_refine_coarse_stage2_continuity=" +
      to_string(timing.component_graph_refine_coarse_stage2_continuity_s) +
      " s, component_graph_refine_coarse_stage2_anchor=" +
      to_string(timing.component_graph_refine_coarse_stage2_anchor_s) +
      " s, component_graph_refine_coarse_stage2_cut_score=" +
      to_string(timing.component_graph_refine_coarse_stage2_cut_score_s) +
      " s, component_graph_refine_build_points=" +
      to_string(timing.component_graph_refine_build_points_s) +
      " s, component_graph_refine_continuity=" +
      to_string(timing.component_graph_refine_continuity_s) +
      " s, component_graph_refine_anchor=" +
      to_string(timing.component_graph_refine_anchor_s) +
      " s, component_graph_refine_anchor_spatial_build=" +
      to_string(timing.component_graph_refine_anchor_spatial_build_s) +
      " s, component_graph_refine_anchor_spatial_query=" +
      to_string(timing.component_graph_refine_anchor_spatial_query_s) +
      " s, component_graph_refine_anchor_dense_fallback=" +
      to_string(timing.component_graph_refine_anchor_dense_fallback_s) +
      " s, pbctopo_anchor_context_build=" +
      to_string(timing.pbctopo_anchor_context_build_s) +
      " s, component_graph_refine_anchor_order=" +
      to_string(timing.component_graph_refine_anchor_order_s) +
      " s, component_graph_refine_anchor_contact=" +
      to_string(timing.component_graph_refine_anchor_contact_s) +
      " s, component_graph_refine_cut_score=" +
      to_string(timing.component_graph_refine_cut_score_s) +
      " s, component_graph_refine_ultraclash=" +
      to_string(timing.component_graph_refine_ultraclash_s) +
      " s, component_graph_refine_stats=" +
      to_string(timing.component_graph_refine_stats_s) +
      " s, component_cycle_sync=" +
      to_string(timing.component_cycle_sync_s) +
      " s, component_search=" + to_string(timing.component_search_s) +
      " s, offset_extract=" + to_string(timing.component_offset_extract_s) +
      " s, apply=" + to_string(timing.apply_s) + " s");
  append_pbctopo_log_line(
      "[pbctopo timing] coarse_candidates=" +
      to_string(timing.coarse_candidates) + ", refine_candidates=" +
      to_string(timing.refine_candidates) + ", projection_clash_skip_frames=" +
      to_string(timing.projection_clash_skip_frames) +
      ", projection_clash_skip_evals=" +
      to_string(timing.projection_clash_skip_evals) + ", component_frames=" +
      to_string(timing.component_frames) + ", component_shift_trials=" +
      to_string(timing.component_shift_trials) +
      ", hull_rescue_attempts=" + to_string(timing.hull_rescue_attempts) +
      ", hull_rescue_successes=" + to_string(timing.hull_rescue_successes) +
      ", distance_rescue_attempts=" +
      to_string(timing.distance_rescue_attempts) +
      ", distance_rescue_successes=" +
      to_string(timing.distance_rescue_successes) +
      ", pbci_fallback_frames=" + to_string(timing.pbci_fallback_frames) +
      ", pbctopo_llps_local_frames=" +
      to_string(timing.pbctopo_llps_local_frames) +
      ", pbctopo_llps_local_units=" +
      to_string(timing.pbctopo_llps_local_units) +
      ", pbctopo_llps_local_unit_source=" + llps_local_unit_source +
      ", pbctopo_llps_local_edges_used=" +
      to_string(timing.pbctopo_llps_local_edges_used) +
      ", pbctopo_llps_local_edges_skipped_cross_unit=" +
      to_string(timing.pbctopo_llps_local_edges_skipped_cross_unit) +
      ", pbctopo_llps_local_disconnected_units=" +
      to_string(timing.pbctopo_llps_local_disconnected_units) +
      ", early_projection_used=" + to_string(timing.early_projection_used) +
      ", previous_seed_used=" + to_string(timing.previous_seed_used) +
      ", pbctopo_kinematic_seed_candidates=" +
      to_string(timing.pbctopo_kinematic_seed_candidates) +
      ", pbctopo_kinematic_seed_evaluated=" +
      to_string(timing.pbctopo_kinematic_seed_evaluated) +
      ", pbctopo_kinematic_seed_deduped=" +
      to_string(timing.pbctopo_kinematic_seed_deduped) +
      ", pbctopo_kinematic_seed_best=" +
      to_string(timing.pbctopo_kinematic_seed_best) +
      ", pbctopo_kinematic_seed_history_resets=" +
      to_string(timing.pbctopo_kinematic_seed_history_resets) +
      ", pbctopo_density_cut_planes_frames=" +
      to_string(timing.pbctopo_density_cut_planes_frames) +
      ", pbctopo_density_cut_planes_x=" +
      to_string(timing.pbctopo_density_cut_planes_x) +
      ", pbctopo_density_cut_planes_y=" +
      to_string(timing.pbctopo_density_cut_planes_y) +
      ", pbctopo_density_cut_planes_z=" +
      to_string(timing.pbctopo_density_cut_planes_z) +
      ", pbctopo_density_cut_planes_dense_equiv_x=" +
      to_string(timing.pbctopo_density_cut_planes_dense_equiv_x) +
      ", pbctopo_density_cut_planes_dense_equiv_y=" +
      to_string(timing.pbctopo_density_cut_planes_dense_equiv_y) +
      ", pbctopo_density_cut_planes_dense_equiv_z=" +
      to_string(timing.pbctopo_density_cut_planes_dense_equiv_z) +
      ", pbctopo_density_cut_plane_collisions=" +
      to_string(timing.pbctopo_density_cut_plane_collisions) +
      ", pbctopo_density_cut_plane_dense_fallbacks=" +
      to_string(timing.pbctopo_density_cut_plane_dense_fallbacks) +
      ", pbctopo_virtual_discrete_frames=" +
      to_string(timing.pbctopo_virtual_discrete_frames) +
      ", pbctopo_virtual_axis_states_x=" +
      to_string(timing.pbctopo_virtual_axis_states_x) +
      ", pbctopo_virtual_axis_states_y=" +
      to_string(timing.pbctopo_virtual_axis_states_y) +
      ", pbctopo_virtual_axis_states_z=" +
      to_string(timing.pbctopo_virtual_axis_states_z) +
      ", pbctopo_virtual_state_space=" +
      to_string(timing.pbctopo_virtual_state_space) +
      ", pbctopo_virtual_states_evaluated=" +
      to_string(timing.pbctopo_virtual_states_evaluated) +
      ", pbctopo_virtual_search_exhaustive_frames=" +
      to_string(timing.pbctopo_virtual_search_exhaustive_frames) +
      ", pbctopo_virtual_search_shortlist_frames=" +
      to_string(timing.pbctopo_virtual_search_shortlist_frames) +
      ", pbctopo_virtual_branch_and_bound_frames=" +
      to_string(timing.pbctopo_virtual_branch_and_bound_frames) +
      ", pbctopo_virtual_sweep_workspace_reuse_frames=" +
      to_string(timing.pbctopo_virtual_sweep_workspace_reuse_frames) +
      ", pbctopo_virtual_sweep_events_x=" +
      to_string(timing.pbctopo_virtual_sweep_events_x) +
      ", pbctopo_virtual_sweep_events_y=" +
      to_string(timing.pbctopo_virtual_sweep_events_y) +
      ", pbctopo_virtual_sweep_events_z=" +
      to_string(timing.pbctopo_virtual_sweep_events_z) +
      ", pbctopo_virtual_sweep_workspace_peak_bytes=" +
      to_string(timing.pbctopo_virtual_sweep_workspace_peak_bytes) +
      ", pbctopo_verify_geometry_buffer_reuses=" +
      to_string(timing.pbctopo_verify_geometry_buffer_reuses) +
      ", pbctopo_frame_point_buffer_peak_bytes_estimate=" +
      to_string(timing.pbctopo_frame_point_buffer_peak_bytes_estimate) +
      ", pbctopo_virtual_omitted_seam_lower_bound_min=" +
      (isfinite(timing.pbctopo_virtual_omitted_seam_lower_bound_min)
           ? to_string(timing.pbctopo_virtual_omitted_seam_lower_bound_min)
           : string("na")) +
      ", pbctopo_virtual_lower_bound_gap_min=" +
      (isfinite(timing.pbctopo_virtual_lower_bound_gap_min)
           ? to_string(timing.pbctopo_virtual_lower_bound_gap_min)
           : string("na")) +
       ", component_graph_edges=" + to_string(timing.component_graph_edges) +
       ", component_graph_shift_descriptors=" +
       to_string(timing.component_graph_shift_descriptors) +
       ", component_graph_shift_materialized=" +
       to_string(timing.component_graph_shift_materialized) +
       ", component_graph_bbox_only_scores=" +
       to_string(timing.component_graph_bbox_only_scores) +
        ", component_graph_used=" + to_string(timing.component_graph_used) +
        ", component_graph_primary_attempts=" +
        to_string(timing.component_graph_primary_attempts) +
        ", component_graph_primary_certified=" +
        to_string(timing.component_graph_primary_certified) +
        ", component_graph_primary_virtual_skips=" +
        to_string(timing.component_graph_primary_virtual_skips) +
        ", component_graph_primary_fallbacks=" +
        to_string(timing.component_graph_primary_fallbacks) +
        ", component_image_shell_expansions=" +
        to_string(timing.component_image_shell_expansions) +
        ", component_image_shell_retry_successes=" +
        to_string(timing.component_image_shell_retry_successes) +
        ", component_image_shell_max_radius=" +
        to_string(timing.component_image_shell_max_radius) +
        ", component_image_shell_max_shift_count=" +
        to_string(timing.component_image_shell_max_shift_count) +
         ", component_warm_used=" + to_string(timing.component_warm_used) +
       ", component_warm_skip_search=" + to_string(timing.component_warm_skip_search) +
       ", component_search_frames=" + to_string(timing.component_search_frames) +
       ", component_graph_failed=" + to_string(timing.component_graph_failed) +
       ", component_graph_points_mismatch=" +
       to_string(timing.component_graph_points_mismatch) +
       ", component_graph_invalid=" + to_string(timing.component_graph_invalid) +
       ", component_graph_ultraclash=" +
       to_string(timing.component_graph_ultraclash) +
       ", component_graph_not_better_not_close=" +
       to_string(timing.component_graph_not_better_not_close) +
        ", component_refine_warm_shift_hits=" +
        to_string(timing.component_refine_warm_shift_hits) +
       ", component_refine_shortlist_hint_hits=" +
       to_string(timing.component_refine_shortlist_hint_hits) +
       ", component_pair_warm_shift_hits=" +
       to_string(timing.component_pair_warm_shift_hits) +
       ", component_refine_fanout_reduced=" +
       to_string(timing.component_refine_fanout_reduced) +
       ", component_graph_refine_score_calls=" +
       to_string(timing.component_graph_refine_score_calls) +
      ", component_graph_refine_coarse_score_calls=" +
      to_string(timing.component_graph_refine_coarse_score_calls) +
      ", component_graph_refine_coarse_stage1_score_calls=" +
      to_string(timing.component_graph_refine_coarse_stage1_score_calls) +
      ", component_graph_refine_coarse_stage2_score_calls=" +
      to_string(timing.component_graph_refine_coarse_stage2_score_calls) +
      ", component_graph_refine_exact_rescore_calls=" +
      to_string(timing.component_graph_refine_exact_rescore_calls) +
      ", component_graph_refine_preclash_prunes=" +
      to_string(timing.component_graph_refine_preclash_prunes) +
      ", component_graph_refine_trial_build_incremental=" +
      to_string(timing.component_graph_refine_trial_build_incremental) +
      ", component_graph_refine_trial_build_full=" +
      to_string(timing.component_graph_refine_trial_build_full) +
      ", component_graph_refine_trial_build_global_rebase=" +
      to_string(timing.component_graph_refine_trial_build_global_rebase) +
      ", component_graph_refine_stats_cached=" +
      to_string(timing.component_graph_refine_stats_cached) +
      ", component_graph_refine_stats_full=" +
      to_string(timing.component_graph_refine_stats_full) +
      ", component_graph_refine_anchor_spatial_calls=" +
      to_string(timing.component_graph_refine_anchor_spatial_calls) +
      ", component_graph_refine_anchor_spatial_fallbacks=" +
      to_string(timing.component_graph_refine_anchor_spatial_fallbacks) +
      ", component_graph_refine_anchor_spatial_disconnects=" +
      to_string(timing.component_graph_refine_anchor_spatial_disconnects) +
      ", component_cycle_sync_attempts=" +
      to_string(timing.component_cycle_sync_attempts) +
      ", component_cycle_sync_used=" +
      to_string(timing.component_cycle_sync_used) +
      ", component_cycle_sync_pair_candidates=" +
      to_string(timing.component_cycle_sync_pair_candidates) +
      ", component_cycle_sync_cycle_support=" +
      to_string(timing.component_cycle_sync_cycle_support) +
      ", component_cycle_sync_cycle_mismatch=" +
      to_string(timing.component_cycle_sync_cycle_mismatch) +
      ", component_cycle_sync_inconsistent_edges=" +
      to_string(timing.component_cycle_sync_inconsistent_edges) +
      ", pbctopo_topology_bonds_available=" +
      to_string(timing.pbctopo_topology_bonds_available) +
      ", pbctopo_topology_bond_edges_added=" +
      to_string(timing.pbctopo_topology_bond_edges_added) +
      ", pbctopo_topology_bond_edges_skipped=" +
      to_string(timing.pbctopo_topology_bond_edges_skipped) +
      ", pbctopo_geometry_edges_added=" +
      to_string(timing.pbctopo_geometry_edges_added) +
      ", pbctopo_topology_duplicate_edges=" +
      to_string(timing.pbctopo_topology_duplicate_edges) +
      ", pbctopo_topology_projected_bonds=" +
      to_string(timing.pbctopo_topology_projected_bonds) +
      ", pbctopo_topology_static_graph_frames=" +
      to_string(timing.pbctopo_topology_static_graph_frames) +
      ", pbctopo_metadata_graph_frames=" +
      to_string(timing.pbctopo_metadata_graph_frames) +
      ", pbctopo_metadata_residue_groups=" +
      to_string(timing.pbctopo_metadata_residue_groups) +
      ", pbctopo_metadata_intra_edges=" +
      to_string(timing.pbctopo_metadata_intra_edges) +
      ", pbctopo_metadata_inter_edges=" +
      to_string(timing.pbctopo_metadata_inter_edges) +
      ", pbctopo_geometry_fallback_frames=" +
      to_string(timing.pbctopo_geometry_fallback_frames) +
      ", pbctopo_topology_source=" +
      pbctopo_topology_source_label(timing) +
      ", pbctopo_topology_component_source=" +
      string((timing.frames > 0 &&
              timing.pbctopo_topology_static_graph_frames == timing.frames)
                 ? "topology"
                 : ((timing.frames > 0 &&
                     timing.pbctopo_metadata_graph_frames == timing.frames)
                        ? "metadata"
                        : ((timing.pbctopo_topology_bond_edges_added > 0 ||
                            timing.pbctopo_metadata_graph_frames > 0)
                               ? "mixed"
                               : "geometry"))) +
      ", pbctopo_anchor_order_edges=" +
      to_string(timing.pbctopo_anchor_order_edges) +
      ", pbctopo_anchor_endpoint_edges=" +
      to_string(timing.pbctopo_anchor_endpoint_edges) +
      ", pbctopo_anchor_contact_pairs=" +
      to_string(timing.pbctopo_anchor_contact_pairs) +
      ", pbctopo_soft_observed_contact_pairs=" +
      to_string(timing.pbctopo_soft_observed_contact_pairs) +
      ", pbctopo_hard_required_owner_pairs=" +
      to_string(timing.pbctopo_hard_required_owner_pairs) +
      ", pbctopo_anchor_residue_order_violations=" +
      to_string(timing.pbctopo_anchor_residue_order_violations) +
      ", pbctopo_anchor_ambiguous_frames=" +
      to_string(timing.pbctopo_anchor_ambiguous_frames) +
      ", component_context_builds=" +
      to_string(timing.component_context_builds) +
      ", component_offset_extracts=" +
      to_string(timing.component_offset_extracts) + ", apply_calls=" +
      to_string(timing.apply_calls) +
      ", component_partial_checks=" +
      to_string(timing.component_partial_checks) +
      ", component_partial_prunes=" +
      to_string(timing.component_partial_prunes) +
      ", component_bbox_skips=" + to_string(timing.component_bbox_skips) +
      ", component_score_calls=" +
      to_string(timing.component_score_calls) + ", component_partial_time=" +
      to_string(timing.component_partial_s) + " s, component_score_time=" +
      to_string(timing.component_score_s) + " s");
  append_pbctopo_log_line(
      "[pbctopo timing] component_pair_table_builds=" +
      to_string(timing.component_pair_table_builds) +
      ", component_pair_table_time=" +
      to_string(timing.component_pair_table_s) +
      " s, component_pair_table_prunes=" +
      to_string(timing.component_pair_table_prunes) +
      ", component_pair_table_component_count=" +
      to_string(timing.component_pair_table_component_count) +
      ", component_pair_table_shift_count=" +
      to_string(timing.component_pair_table_shift_count) +
      ", component_pair_table_dense_entries_equivalent=" +
      to_string(timing.component_pair_table_dense_entries_equivalent) +
      ", component_pair_table_materialized_entries=" +
      to_string(timing.component_pair_table_materialized_entries) +
      ", component_pair_table_bbox_prefiltered=" +
      to_string(timing.component_pair_table_bbox_prefiltered) +
      ", component_pair_table_sparse_rows=" +
      to_string(timing.component_pair_table_sparse_rows) +
      ", component_pair_table_peak_bytes_estimate=" +
      to_string(timing.component_pair_table_peak_bytes_estimate) +
      ", component_search_zero_current_mask_prunes=" +
      to_string(timing.component_search_zero_current_mask_prunes) +
      ", component_search_zero_future_mask_prunes=" +
      to_string(timing.component_search_zero_future_mask_prunes) +
      ", component_search_fixed_minor_frames=" +
      to_string(timing.component_search_fixed_minor_frames) +
      ", component_search_fixed_minor_components=" +
      to_string(timing.component_search_fixed_minor_components) +
      ", component_search_seed_warm=" +
      to_string(timing.component_search_seed_warm) +
      ", component_search_seed_graph=" +
      to_string(timing.component_search_seed_graph) +
      ", component_search_seed_none=" +
      to_string(timing.component_search_seed_none) +
      ", component_search_prefix_continuity_checks=" +
      to_string(timing.component_search_prefix_continuity_checks) +
      ", component_search_prefix_continuity_prunes=" +
      to_string(timing.component_search_prefix_continuity_prunes) +
      ", component_search_prefix_continuity_hard_rejects=" +
      to_string(timing.component_search_prefix_continuity_hard_rejects) +
      ", component_search_exact_budget_skips=" +
      to_string(timing.component_search_exact_budget_skips) +
      ", component_search_exact_state_upper_bound=" +
      to_string(timing.component_search_exact_state_upper_bound) +
      ", component_search_prefix_continuity_time=" +
      to_string(timing.component_search_prefix_continuity_s) +
      " s" +
      ", component_large_mode_merge_frames=" +
      to_string(timing.component_large_mode_merge_frames) +
      ", component_large_mode_pair_only_frames=" +
      to_string(timing.component_large_mode_pair_only_frames) +
      ", component_large_mode_beam_frames=" +
      to_string(timing.component_large_mode_beam_frames) +
      ", component_large_mode_merge_count=" +
      to_string(timing.component_large_mode_merge_count) +
      ", component_large_mode_internal_offsets_preserved=" +
      to_string(timing.component_large_mode_internal_offsets_preserved) +
      ", component_large_mode_internal_offset_nonzero=" +
      to_string(timing.component_large_mode_internal_offset_nonzero) +
      ", component_large_mode_graph_seed_count=" +
      to_string(timing.component_large_mode_graph_seed_count) +
      ", component_large_mode_graph_seed_refined=" +
      to_string(timing.component_large_mode_graph_seed_refined) +
      ", component_large_mode_graph_seed_nonbaseline_win=" +
      to_string(timing.component_large_mode_graph_seed_nonbaseline_win) +
      ", component_large_mode_beam_expanded=" +
      to_string(timing.component_large_mode_beam_expanded) +
      ", component_large_mode_beam_kept=" +
      to_string(timing.component_large_mode_beam_kept) +
      ", component_large_mode_beam_complete_assignments=" +
      to_string(timing.component_large_mode_beam_complete_assignments) +
      ", component_large_mode_beam_final_rescored=" +
      to_string(timing.component_large_mode_beam_final_rescored) +
      ", component_large_mode_fallback_frames=" +
      to_string(timing.component_large_mode_fallback_frames) +
      ", component_search_shift_warm_hits=" +
      to_string(timing.component_search_shift_warm_hits) +
      ", pbctopo_certification_attempts=" +
      to_string(timing.pbctopo_certification_attempts) +
      ", pbctopo_certification_successes=" +
      to_string(timing.pbctopo_certification_successes) +
      ", pbctopo_certification_assignment_rejects=" +
      to_string(timing.pbctopo_certification_assignment_rejects) +
      ", pbctopo_certification_hard_rejects=" +
      to_string(timing.pbctopo_certification_hard_rejects) +
      ", pbctopo_certification_topology_edges_checked=" +
      to_string(timing.pbctopo_certification_topology_edges_checked) +
      ", pbctopo_certification_topology_edge_residuals=" +
      to_string(timing.pbctopo_certification_topology_edge_residuals) +
      ", pbctopo_certification_atom_image_residuals=" +
      to_string(timing.pbctopo_certification_atom_image_residuals) +
      ", pbctopo_certification_required_contacts_checked=" +
      to_string(timing.pbctopo_certification_required_contacts_checked) +
      ", pbctopo_certification_required_contacts_lost=" +
      to_string(timing.pbctopo_certification_required_contacts_lost) +
      ", pbctopo_diagnostic_only_candidates=" +
      to_string(timing.pbctopo_diagnostic_only_candidates) +
      ", pbctopo_virtual_geometry_reference_evals=" +
      to_string(timing.pbctopo_virtual_geometry_reference_evals) +
      ", pbctopo_virtual_geometry_reference_safe=" +
      to_string(timing.pbctopo_virtual_geometry_reference_safe) +
      ", pbctopo_virtual_geometry_reference_clashing=" +
      to_string(timing.pbctopo_virtual_geometry_reference_clashing) +
      ", pbctopo_virtual_geometry_outlier_rejects=" +
      to_string(timing.pbctopo_virtual_geometry_outlier_rejects) +
      ", pbctopo_virtual_geometry_outlier_observations=" +
      to_string(timing.pbctopo_virtual_geometry_outlier_observations) +
      ", component_graph_outlier_rejects=" +
      to_string(timing.component_graph_outlier_rejects) +
      ", component_search_outlier_rejects=" +
      to_string(timing.component_search_outlier_rejects) +
      ", component_graph_outlier_observations=" +
      to_string(timing.component_graph_outlier_observations) +
      ", component_search_outlier_observations=" +
      to_string(timing.component_search_outlier_observations) +
      ", component_score_continuity_time=" +
      to_string(timing.component_score_continuity_s) +
      " s, component_score_anchor_time=" +
      to_string(timing.component_score_anchor_s) +
      " s, component_score_stats_principal_time=" +
      to_string(timing.component_score_stats_principal_s) +
      " s, component_score_component_metrics_time=" +
      to_string(timing.component_score_component_metrics_s) +
      " s, component_score_continuity_rejects=" +
      to_string(timing.component_score_continuity_rejects) +
      ", component_score_anchor_rejects=" +
      to_string(timing.component_score_anchor_rejects) +
      ", component_score_prefix_prunes_after_continuity=" +
      to_string(timing.component_score_prefix_prunes_after_continuity) +
      ", component_score_prefix_prunes_after_anchor=" +
      to_string(timing.component_score_prefix_prunes_after_anchor) +
      ", component_score_prefix_prunes_after_principal=" +
      to_string(timing.component_score_prefix_prunes_after_principal) +
      ", component_score_prefix_prunes_after_component_metrics=" +
      to_string(timing.component_score_prefix_prunes_after_component_metrics) +
      ", component_score_prefix_prunes_after_final=" +
      to_string(timing.component_score_prefix_prunes_after_final) +
      ", component_search_depth_visits=[d0:" +
      to_string(timing.component_search_depth_visits[0]) + ",d1:" +
      to_string(timing.component_search_depth_visits[1]) + ",d2:" +
      to_string(timing.component_search_depth_visits[2]) + ",d3:" +
      to_string(timing.component_search_depth_visits[3]) + ",d4:" +
      to_string(timing.component_search_depth_visits[4]) + ",d5:" +
      to_string(timing.component_search_depth_visits[5]) + ",d6:" +
      to_string(timing.component_search_depth_visits[6]) + ",d7:" +
      to_string(timing.component_search_depth_visits[7]) + "]");
  append_pbctopo_log_line(
      "[pbctopo math] pair_shift_score=normalized, shift_metric="
      "cartesian_component_mean_A2, gauge=closest_lattice_mean, cycle_residual="
      "fundamental_basis, anchor_metrics=" +
      string(titan_pbctopo::anchor_metric_mode_name(anchor_metric_mode)) +
      ", clash_model=vdw_radius, vdw_pair_checks=" +
      to_string(timing.pbctopo_vdw_pair_checks) +
      ", vdw_radius_pairs=" + to_string(timing.pbctopo_vdw_radius_pairs) +
      ", vdw_fallback_pairs=" +
      to_string(timing.pbctopo_vdw_radius_fallback_pairs) +
      ", temporal_mode=" +
      string(titan_pbctopo::temporal_mode_name(temporal_mode)) +
      ", temporal_candidate_frames=" +
      to_string(timing.pbctopo_temporal_candidate_frames) +
      ", temporal_candidates=" +
      to_string(timing.pbctopo_temporal_candidates) +
      ", temporal_certified_candidates=" +
      to_string(timing.pbctopo_temporal_certified_candidates) +
      ", temporal_assignment_deduped=" +
      to_string(timing.pbctopo_temporal_assignment_deduped) +
      ", temporal_diagnostic_only_frames=" +
      to_string(timing.pbctopo_temporal_diagnostic_only_frames) +
      ", temporal_proposed_changes=" +
      to_string(timing.pbctopo_temporal_proposed_changes) +
      ", temporal_path_changes=" +
      to_string(timing.pbctopo_temporal_path_changes) +
      ", temporal_replay_frames=" +
      to_string(timing.pbctopo_temporal_replay_frames) +
      ", temporal_replay_applied=" +
      to_string(timing.pbctopo_temporal_replay_applied) +
      ", temporal_replay_skipped=" +
      to_string(timing.pbctopo_temporal_replay_skipped));
  append_pbctopo_log_line("[pbctopo] Valid candidate evaluations: " +
                          to_string(valid_candidates));
  append_pbctopo_log_line("[pbctopo] end run");
  return 0;
}
