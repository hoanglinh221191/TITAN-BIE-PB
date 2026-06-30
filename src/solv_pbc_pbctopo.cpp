#include "Multi_MOL.h"
#include "BoxMatrix.h"
#include "aligned_allocator.h"
#include "pbctopo/config.h"
#include "pbctopo/logging.h"
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
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
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
using PbctopoRuntimeCutoffs = titan_pbctopo::RuntimeCutoffs;

constexpr int kTopoImageCoeff[27][3] = {
    {1, 0, 0},   {-1, 0, 0},  {0, 1, 0},   {0, -1, 0},  {0, 0, 1},
    {0, 0, -1},  {1, 1, 0},   {1, -1, 0},  {-1, 1, 0},  {-1, -1, 0},
    {1, 0, 1},   {1, 0, -1},  {-1, 0, 1},  {-1, 0, -1}, {0, 1, 1},
    {0, -1, 1},  {0, 1, -1},  {0, -1, -1}, {1, 1, 1},   {-1, 1, 1},
    {1, -1, 1},  {1, 1, -1},  {-1, -1, 1}, {-1, 1, -1}, {1, -1, -1},
    {-1, -1, -1}, {0, 0, 0}};

constexpr double kDistanceEps = 1.0e-12;
constexpr int kPackedCellBias = 1 << 20;
constexpr uint64_t kPackedCellMask = (1ull << 21) - 1ull;
constexpr size_t kPbctopoSpatialAnchorThreshold = 512;
constexpr size_t kPbctopoSpatialAnchorK = 12;
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
constexpr size_t kPbctopoBeamFinalRescore = 4;

void reset_pbctopo_log_file() {
  titan_pbctopo::reset_log_file();
}

void append_pbctopo_log_line(const string &message) {
  titan_pbctopo::append_log_line(message);
}

struct PbctopoRange {
  size_t begin = 0;
  size_t end = 0;
  size_t source_index = 0;
  bool topo_backed = false;
};

struct PbctopoAtom {
  size_t global_atom_id = 0;
  size_t fragment_id = 0;
  int owner = 0;
  size_t local_atom_index = 0;
  uint64_t chain_key = 0;
  int residue_id = numeric_limits<int>::min();
  bool has_topology_key = false;
  bool is_heavy = true;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double sx = 0.0;
  double sy = 0.0;
  double sz = 0.0;
};

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

struct PbctopoShift {
  size_t shift_index = 26;
  int ix = 0;
  int iy = 0;
  int iz = 0;
  double tx = 0.0;
  double ty = 0.0;
  double tz = 0.0;
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
};

using PbctopoGyrationShapeMetrics = GyrationShapeMetrics;

struct PbctopoCutCandidate {
  double cut_x = 0.0;
  double cut_y = 0.0;
  double cut_z = 0.0;
  double gap_score = 0.0;
  double bbox_volume = numeric_limits<double>::infinity();
  double rg2 = numeric_limits<double>::infinity();
  double shape_kappa2 = numeric_limits<double>::infinity();
  double shape_asphericity = numeric_limits<double>::infinity();
  double shape_acylindricity = numeric_limits<double>::infinity();
  double shift_norm2 = numeric_limits<double>::infinity();
  double min_pair_d2 = numeric_limits<double>::infinity();
  double continuity_path = numeric_limits<double>::infinity();
  double continuity_max_d2 = numeric_limits<double>::infinity();
  double anchor_mst2 = numeric_limits<double>::infinity();
  double anchor_max_d2 = numeric_limits<double>::infinity();
  double principal_rg2_max = numeric_limits<double>::infinity();
  double principal_rg2_mid = numeric_limits<double>::infinity();
  double principal_rg2_min = numeric_limits<double>::infinity();
  double edge_stretch = numeric_limits<double>::infinity();
  double max_owner_raw_minus_minimg_delta = 0.0;
  size_t broken_edges = numeric_limits<size_t>::max();
  size_t split_owner_pairs = 0;
  size_t clash_atom_a = numeric_limits<size_t>::max();
  size_t clash_atom_b = numeric_limits<size_t>::max();
  bool ultra_clash = false;
  bool valid = false;
};

struct PbctopoLocalEdge {
  size_t a = 0;
  size_t b = 0;
  double ref_d2 = 0.0;
  double inv_ref_d2 = 0.0;
  double break_cutoff_d2 = 0.0;
};

struct PbctopoContinuityPair {
  size_t a = 0;
  size_t b = 0;
};

struct PbctopoComponentCandidate {
  double continuity_path = numeric_limits<double>::infinity();
  double continuity_max_d2 = numeric_limits<double>::infinity();
  double anchor_mst2 = numeric_limits<double>::infinity();
  double anchor_max_d2 = numeric_limits<double>::infinity();
  double principal_rg2_max = numeric_limits<double>::infinity();
  double principal_rg2_mid = numeric_limits<double>::infinity();
  double principal_rg2_min = numeric_limits<double>::infinity();
  double shape_kappa2 = numeric_limits<double>::infinity();
  double shape_asphericity = numeric_limits<double>::infinity();
  double shape_acylindricity = numeric_limits<double>::infinity();
  double component_mst2 = numeric_limits<double>::infinity();
  double component_rg2 = numeric_limits<double>::infinity();
  double component_kappa2 = numeric_limits<double>::infinity();
  double component_asphericity = numeric_limits<double>::infinity();
  double component_acylindricity = numeric_limits<double>::infinity();
  double bbox_volume = numeric_limits<double>::infinity();
  double rg2 = numeric_limits<double>::infinity();
  double shift_norm2 = numeric_limits<double>::infinity();
  double min_pair_d2 = numeric_limits<double>::infinity();
  size_t clash_atom_a = numeric_limits<size_t>::max();
  size_t clash_atom_b = numeric_limits<size_t>::max();
  bool ultra_clash = false;
  bool valid = false;
};

struct Int3 {
  int x = 0;
  int y = 0;
  int z = 0;
};

struct PbctopoPairShiftScore {
  size_t shift_index = 0;
  int ix = 0;
  int iy = 0;
  int iz = 0;
  size_t contact_atoms = 0;
  double soft_penalty = 0.0;
  double min_pair_d2 = numeric_limits<double>::infinity();
  double bbox_d2 = numeric_limits<double>::infinity();
  double shift_norm2 = numeric_limits<double>::infinity();
  double quality_score = -numeric_limits<double>::infinity();
  bool ultra_clash = false;
  bool valid = false;
};

struct PbctopoGraphEdge {
  size_t a = 0;
  size_t b = 0;
  int dx = 0;
  int dy = 0;
  int dz = 0;
  size_t option_rank = 0;
  size_t contact_atoms = 0;
  double soft_penalty = 0.0;
  double min_pair_d2 = numeric_limits<double>::infinity();
  double bbox_d2 = numeric_limits<double>::infinity();
  double shift_norm2 = numeric_limits<double>::infinity();
  double confidence = 0.0;
  int cycle_support = 0;
  int cycle_mismatch = 0;
  bool ultra_clash = false;
  bool valid = false;
};

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

struct PbctopoTimingStats {
  size_t frames = 0;
  size_t coarse_candidates = 0;
  size_t refine_candidates = 0;
  size_t projection_clash_skip_frames = 0;
  size_t projection_clash_skip_evals = 0;
  size_t component_frames = 0;
  size_t component_shift_trials = 0;
  size_t component_graph_refine_coarse_score_calls = 0;
  size_t component_graph_refine_coarse_stage1_score_calls = 0;
  size_t component_graph_refine_coarse_stage2_score_calls = 0;
  size_t component_graph_refine_exact_rescore_calls = 0;
  size_t component_partial_checks = 0;
  size_t component_partial_prunes = 0;
  size_t component_bbox_skips = 0;
  size_t component_score_calls = 0;
  size_t component_pair_table_builds = 0;
  size_t component_pair_table_prunes = 0;
  size_t component_search_zero_current_mask_prunes = 0;
  size_t component_search_zero_future_mask_prunes = 0;
  size_t component_search_fixed_minor_components = 0;
  size_t component_search_fixed_minor_frames = 0;
  size_t component_search_shift_warm_hits = 0;
  size_t component_search_seed_warm = 0;
  size_t component_search_seed_graph = 0;
  size_t component_search_seed_none = 0;
  size_t component_large_mode_merge_frames = 0;
  size_t component_large_mode_pair_only_frames = 0;
  size_t component_large_mode_beam_frames = 0;
  size_t component_large_mode_merge_count = 0;
  size_t component_large_mode_graph_seed_count = 0;
  size_t component_large_mode_graph_seed_refined = 0;
  size_t component_large_mode_graph_seed_nonbaseline_win = 0;
  size_t component_large_mode_beam_expanded = 0;
  size_t component_large_mode_beam_kept = 0;
  size_t component_large_mode_beam_complete_assignments = 0;
  size_t component_large_mode_beam_final_rescored = 0;
  size_t component_large_mode_fallback_frames = 0;
  size_t component_graph_outlier_rejects = 0;
  size_t component_search_outlier_rejects = 0;
  size_t component_score_continuity_rejects = 0;
  size_t component_score_anchor_rejects = 0;
  size_t component_score_prefix_prunes_after_continuity = 0;
  size_t component_score_prefix_prunes_after_anchor = 0;
  size_t component_score_prefix_prunes_after_principal = 0;
  size_t component_score_prefix_prunes_after_component_metrics = 0;
  size_t component_score_prefix_prunes_after_final = 0;
  size_t component_graph_edges = 0;
  size_t component_graph_shift_descriptors = 0;
  size_t component_graph_shift_materialized = 0;
  size_t component_graph_bbox_only_scores = 0;
  size_t component_graph_used = 0;
  size_t component_warm_used = 0;
  size_t component_warm_skip_search = 0;
  size_t component_search_frames = 0;
  size_t component_graph_failed = 0;
  size_t component_graph_points_mismatch = 0;
  size_t component_graph_invalid = 0;
  size_t component_graph_ultraclash = 0;
  size_t component_graph_not_better_not_close = 0;
  size_t component_refine_warm_shift_hits = 0;
  size_t component_refine_shortlist_hint_hits = 0;
  size_t component_pair_warm_shift_hits = 0;
  size_t component_refine_fanout_reduced = 0;
  size_t hull_rescue_attempts = 0;
  size_t hull_rescue_successes = 0;
  size_t distance_rescue_attempts = 0;
  size_t distance_rescue_successes = 0;
  size_t pbci_fallback_frames = 0;
  size_t component_graph_refine_score_calls = 0;
  size_t component_graph_refine_preclash_prunes = 0;
  size_t component_graph_refine_trial_build_incremental = 0;
  size_t component_graph_refine_trial_build_full = 0;
  size_t component_graph_refine_trial_build_global_rebase = 0;
  size_t component_graph_refine_stats_cached = 0;
  size_t component_graph_refine_stats_full = 0;
  size_t component_graph_refine_anchor_spatial_calls = 0;
  size_t component_graph_refine_anchor_spatial_fallbacks = 0;
  size_t component_graph_refine_anchor_spatial_disconnects = 0;
  size_t component_context_builds = 0;
  size_t component_offset_extracts = 0;
  size_t apply_calls = 0;
  size_t early_projection_used = 0;
  size_t previous_seed_used = 0;
  array<size_t, 8> component_search_depth_visits{};
  double append_s = 0.0;
  double heavy_subset_s = 0.0;
  double translation_cache_s = 0.0;
  double graph_s = 0.0;
  double projection_edge_s = 0.0;
  double component_split_s = 0.0;
  double event_planes_s = 0.0;
  double seed_s = 0.0;
  double coarse_s = 0.0;
  double refine_s = 0.0;
  double component_prep_s = 0.0;
  double component_search_s = 0.0;
  double component_graph_s = 0.0;
  double component_graph_precompute_s = 0.0;
  double component_graph_lazy_materialize_s = 0.0;
  double component_graph_pair_score_s = 0.0;
  double component_graph_tree_s = 0.0;
  double component_graph_refine_s = 0.0;
  double component_graph_refine_rebase_s = 0.0;
  double component_graph_refine_coarse_s = 0.0;
  double component_graph_refine_coarse_stage1_s = 0.0;
  double component_graph_refine_coarse_stage1_build_points_s = 0.0;
  double component_graph_refine_coarse_stage1_continuity_s = 0.0;
  double component_graph_refine_coarse_stage1_anchor_s = 0.0;
  double component_graph_refine_coarse_stage1_cut_score_s = 0.0;
  double component_graph_refine_coarse_stage1_cut_score_principal_s = 0.0;
  double component_graph_refine_coarse_stage1_cut_score_edge_scan_s = 0.0;
  double component_graph_refine_coarse_stage2_s = 0.0;
  double component_graph_refine_coarse_stage2_build_points_s = 0.0;
  double component_graph_refine_coarse_stage2_continuity_s = 0.0;
  double component_graph_refine_coarse_stage2_anchor_s = 0.0;
  double component_graph_refine_coarse_stage2_cut_score_s = 0.0;
  double component_graph_refine_build_points_s = 0.0;
  double component_graph_refine_continuity_s = 0.0;
  double component_graph_refine_anchor_s = 0.0;
  double component_graph_refine_anchor_spatial_build_s = 0.0;
  double component_graph_refine_anchor_spatial_query_s = 0.0;
  double component_graph_refine_anchor_dense_fallback_s = 0.0;
  double component_graph_refine_cut_score_s = 0.0;
  double component_graph_refine_ultraclash_s = 0.0;
  double component_graph_refine_stats_s = 0.0;
  double early_projection_s = 0.0;
  double warm_component_s = 0.0;
  double component_partial_s = 0.0;
  double component_score_s = 0.0;
  double component_pair_table_s = 0.0;
  double component_score_continuity_s = 0.0;
  double component_score_anchor_s = 0.0;
  double component_score_stats_principal_s = 0.0;
  double component_score_component_metrics_s = 0.0;
  double component_context_s = 0.0;
  double component_offset_extract_s = 0.0;
  double apply_s = 0.0;
};

enum class PbctopoRefineScorePhase {
  Exact,
  CoarseStage1,
  CoarseStage2,
};

enum class PbctopoComponentMode {
  Disabled,
  Exact,
  MergeThenGraph,
  Beam,
};

struct PbctopoCutScoreBreakdown {
  double principal_s = 0.0;
  double edge_scan_s = 0.0;
};

struct PbctopoStage1EdgeContext {
  size_t fixed_broken_edges = 0;
  double fixed_edge_stretch = 0.0;
  vector<size_t> variable_edge_indices;
};

struct PbctopoEdgeScanSoA {
  vector<size_t> a;
  vector<size_t> b;
  vector<double> ref_d2;
  vector<double> inv_ref_d2;
  vector<double> break_cutoff_d2;
};

struct PbctopoSparseAnchorEdge {
  size_t a = 0;
  size_t b = 0;
  double d2 = numeric_limits<double>::infinity();
};

struct PbctopoComponentPairShiftEntry {
  double min_pair_d2 = numeric_limits<double>::infinity();
  double bbox_d2 = numeric_limits<double>::infinity();
  size_t clash_atom_a = numeric_limits<size_t>::max();
  size_t clash_atom_b = numeric_limits<size_t>::max();
  bool hard_clash = false;
};

struct PbctopoComponentPairShiftTable {
  size_t component_count = 0;
  size_t shift_count = 0;
  vector<PbctopoComponentPairShiftEntry> entries;
  vector<uint32_t> compat_masks;

  void reset(size_t new_component_count, size_t new_shift_count) {
    component_count = new_component_count;
    shift_count = new_shift_count;
    entries.assign(component_count * component_count * shift_count * shift_count,
                   {});
    compat_masks.assign(component_count * component_count * shift_count, 0u);
  }

  size_t index(size_t lhs_comp, size_t lhs_shift, size_t rhs_comp,
               size_t rhs_shift) const {
    return ((((lhs_comp * component_count) + rhs_comp) * shift_count) +
            lhs_shift) *
               shift_count +
           rhs_shift;
  }

  size_t compat_index(size_t lhs_comp, size_t lhs_shift, size_t rhs_comp) const {
    return (((lhs_comp * component_count) + rhs_comp) * shift_count) + lhs_shift;
  }

  PbctopoComponentPairShiftEntry &at(size_t lhs_comp, size_t lhs_shift,
                                     size_t rhs_comp, size_t rhs_shift) {
    return entries[index(lhs_comp, lhs_shift, rhs_comp, rhs_shift)];
  }

  const PbctopoComponentPairShiftEntry &
  at(size_t lhs_comp, size_t lhs_shift, size_t rhs_comp,
     size_t rhs_shift) const {
    return entries[index(lhs_comp, lhs_shift, rhs_comp, rhs_shift)];
  }

  uint32_t &compat_mask(size_t lhs_comp, size_t lhs_shift, size_t rhs_comp) {
    return compat_masks[compat_index(lhs_comp, lhs_shift, rhs_comp)];
  }

  uint32_t compat_mask(size_t lhs_comp, size_t lhs_shift, size_t rhs_comp) const {
    return compat_masks[compat_index(lhs_comp, lhs_shift, rhs_comp)];
  }
};

struct PbctopoMergedComponentPlan {
  vector<vector<size_t>> merged_components;
  vector<vector<size_t>> merged_to_raw;
  vector<size_t> raw_to_merged;
};

struct PbctopoGraphSeedQuickScore {
  bool valid = false;
  double continuity_max_d2 = numeric_limits<double>::infinity();
  double continuity_path = numeric_limits<double>::infinity();
  double anchor_max_d2 = numeric_limits<double>::infinity();
  double anchor_mst2 = numeric_limits<double>::infinity();
  double principal_rg2_max = numeric_limits<double>::infinity();
  double principal_rg2_mid = numeric_limits<double>::infinity();
  double principal_rg2_min = numeric_limits<double>::infinity();
  double shift_norm2 = numeric_limits<double>::infinity();
};

enum class PbctopoComponentCandidatePrefixStage {
  Continuity,
  Anchor,
  Principal,
  ComponentMetrics,
  Final,
};

struct PbctopoAnchorKdNode {
  size_t point_idx = 0;
  int axis = 0;
  int left = -1;
  int right = -1;
  double min_x = 0.0;
  double max_x = 0.0;
  double min_y = 0.0;
  double max_y = 0.0;
  double min_z = 0.0;
  double max_z = 0.0;
};

struct PbctopoWarmState {
  PbctopoCutCandidate previous_projection_candidate;
  vector<size_t> previous_component_sizes;
  vector<Int3> previous_component_offsets;
  vector<vector<Int3>> previous_stage1_delta_hints_by_component;
  vector<vector<Int3>> previous_stage2_delta_hints_by_component;
  bool has_projection_candidate = false;
  bool has_component_offsets = false;
  bool has_component_trial_hints = false;
};

struct PbctopoCoarseSeedOption {
  uint64_t state_key = 0;
  PbctopoSeed3D seed;
};

struct PbctopoSystemStaticData {
  size_t atom_count = 0;
  size_t heavy_count = 0;
  vector<char> heavy_mask;
  vector<PbctopoRange> ranges;
  vector<PbctopoContinuityPair> continuity_pairs_local;
  vector<size_t> continuity_anchor_local;
};

struct PbctopoTrajectoryContext {
  size_t total_atoms_hint = 0;
  size_t heavy_atoms_hint = 0;
  size_t local_edges_hint = 0;
  size_t projection_edges_hint = 0;
  size_t component_hint = 0;
  size_t plane_hint = 0;
  size_t virtual_state_hint = 512;
  size_t coarse_seed_hint = 32;
  size_t coarse_candidate_hint = 16;
  size_t continuity_pair_hint = 0;
  PbctopoSystemStaticData rec_static;
  PbctopoSystemStaticData lig_static;
  vector<PbctopoContinuityPair> continuity_pairs;
  vector<size_t> continuity_anchor_indices;
};

struct PbctopoResidueChainKey {
  int owner = 0;
  uint64_t chain_key = 0;

  bool operator==(const PbctopoResidueChainKey &other) const {
    return owner == other.owner && chain_key == other.chain_key;
  }
};

struct PbctopoResidueChainKeyHash {
  size_t operator()(const PbctopoResidueChainKey &key) const {
    const uint64_t owner_bits =
        static_cast<uint64_t>(static_cast<uint32_t>(key.owner));
    const uint64_t mixed = key.chain_key ^ (owner_bits * 0x9e3779b97f4a7c15ull);
    return static_cast<size_t>(mixed ^ (mixed >> 32));
  }
};

vector<int> build_topology_residue_sequence_indices(
    const vector<PbctopoAtom> &atoms) {
  vector<int> residue_sequence_indices(atoms.size(), numeric_limits<int>::min());
  struct ChainState {
    int last_residue_id = numeric_limits<int>::min();
    int current_rank = -1;
    bool initialized = false;
  };
  unordered_map<PbctopoResidueChainKey, ChainState, PbctopoResidueChainKeyHash>
      chain_states;
  chain_states.reserve(max<size_t>(8, atoms.size() / 64));

  for (size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
    const auto &atom = atoms[atom_idx];
    if (!atom.has_topology_key)
      continue;
    const PbctopoResidueChainKey key{atom.owner, atom.chain_key};
    auto &state = chain_states[key];
    if (!state.initialized) {
      state.initialized = true;
      state.last_residue_id = atom.residue_id;
      state.current_rank = 0;
    } else if (atom.residue_id != state.last_residue_id) {
      state.last_residue_id = atom.residue_id;
      ++state.current_rank;
    }
    residue_sequence_indices[atom_idx] = state.current_rank;
  }

  return residue_sequence_indices;
}

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
  unordered_set<uint64_t> seen_virtual_states;
  vector<Vec3> warm_points;
  vector<Vec3> coarse_warm_points;
  vector<Vec3> placed_points;
  vector<Vec3> best_points;
  vector<Vec3> best_safe_points;
  vector<Vec3> refine_points;
  vector<Vec3> verify_points;
  vector<Vec3> trial_points;
  vector<Vec3> component_points;
  vector<PbctopoCutCandidate> coarse_candidates;
  vector<PbctopoCoarseSeedOption> coarse_seed_options;
  vector<PbctopoCutCandidate> coarse_results;
  vector<double> x_planes;
  vector<double> y_planes;
  vector<double> z_planes;
  vector<size_t> projection_index;

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
    seen_virtual_states.reserve(context.virtual_state_hint);
    warm_points.reserve(context.total_atoms_hint);
    coarse_warm_points.reserve(context.heavy_atoms_hint);
    placed_points.reserve(context.heavy_atoms_hint);
    best_points.reserve(context.total_atoms_hint);
    best_safe_points.reserve(context.total_atoms_hint);
    refine_points.reserve(context.heavy_atoms_hint);
    verify_points.reserve(context.total_atoms_hint);
    trial_points.reserve(context.total_atoms_hint);
    component_points.reserve(context.total_atoms_hint);
    coarse_candidates.reserve(context.coarse_candidate_hint);
    coarse_seed_options.reserve(context.coarse_seed_hint);
    coarse_results.reserve(context.coarse_seed_hint);
    x_planes.reserve(context.plane_hint);
    y_planes.reserve(context.plane_hint);
    z_planes.reserve(context.plane_hint);
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
    seen_virtual_states.clear();
    warm_points.clear();
    coarse_warm_points.clear();
    placed_points.clear();
    best_points.clear();
    best_safe_points.clear();
    refine_points.clear();
    verify_points.clear();
    trial_points.clear();
    component_points.clear();
    coarse_candidates.clear();
    coarse_seed_options.clear();
    coarse_results.clear();
    x_planes.clear();
    y_planes.clear();
    z_planes.clear();
    projection_index.clear();
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

struct ClashCellKey {
  int x = 0;
  int y = 0;
  int z = 0;

  bool operator==(const ClashCellKey &rhs) const {
    return x == rhs.x && y == rhs.y && z == rhs.z;
  }
};

struct ClashCellKeyHash {
  size_t operator()(const ClashCellKey &key) const noexcept {
    const uint64_t ux = static_cast<uint32_t>(key.x + (1 << 20));
    const uint64_t uy = static_cast<uint32_t>(key.y + (1 << 20));
    const uint64_t uz = static_cast<uint32_t>(key.z + (1 << 20));
    return static_cast<size_t>((ux << 42) ^ (uy << 21) ^ uz);
  }
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

struct PbctopoShiftedComponent {
  vector<Vec3> points;
  vector<ClashCellKey> clash_keys;
  vector<ClashCellKey> contact_keys;
  unordered_map<uint64_t, vector<size_t>> clash_grid;
  unordered_map<uint64_t, vector<size_t>> contact_grid;
  Vec3 center{};
  double min_x = numeric_limits<double>::infinity();
  double min_y = numeric_limits<double>::infinity();
  double min_z = numeric_limits<double>::infinity();
  double max_x = -numeric_limits<double>::infinity();
  double max_y = -numeric_limits<double>::infinity();
  double max_z = -numeric_limits<double>::infinity();
  double weighted_shift_norm2 = 0.0;
};

struct PbctopoShiftedComponentDescriptor {
  Vec3 center{};
  double min_x = numeric_limits<double>::infinity();
  double min_y = numeric_limits<double>::infinity();
  double min_z = numeric_limits<double>::infinity();
  double max_x = -numeric_limits<double>::infinity();
  double max_y = -numeric_limits<double>::infinity();
  double max_z = -numeric_limits<double>::infinity();
  double weighted_shift_norm2 = 0.0;
};

struct PbctopoComponentBaseGeometry {
  vector<double> x;
  vector<double> y;
  vector<double> z;
  Vec3 center{};
  double min_x = numeric_limits<double>::infinity();
  double min_y = numeric_limits<double>::infinity();
  double min_z = numeric_limits<double>::infinity();
  double max_x = -numeric_limits<double>::infinity();
  double max_y = -numeric_limits<double>::infinity();
  double max_z = -numeric_limits<double>::infinity();
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

Int3 operator+(const Int3 &lhs, const Int3 &rhs) {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Int3 operator-(const Int3 &lhs, const Int3 &rhs) {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

bool is_zero_delta(const Int3 &delta) {
  return delta.x == 0 && delta.y == 0 && delta.z == 0;
}

bool same_delta(const Int3 &lhs, const Int3 &rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

size_t find_shift_index(const vector<PbctopoShift> &shifts, const Int3 &delta) {
  for (size_t idx = 0; idx < shifts.size(); ++idx) {
    if (shifts[idx].ix == delta.x && shifts[idx].iy == delta.y &&
        shifts[idx].iz == delta.z) {
      return idx;
    }
  }
  return numeric_limits<size_t>::max();
}

bool append_shift_index_if_present(const vector<PbctopoShift> &shifts,
                                   const Int3 &delta, bool allow_identity,
                                   vector<size_t> &ordered_indices,
                                   vector<unsigned char> &seen_indices) {
  if (!allow_identity && is_zero_delta(delta))
    return false;
  const size_t shift_idx = find_shift_index(shifts, delta);
  if (shift_idx == numeric_limits<size_t>::max() ||
      shift_idx >= seen_indices.size() || seen_indices[shift_idx] != 0) {
    return false;
  }
  seen_indices[shift_idx] = 1;
  ordered_indices.push_back(shift_idx);
  return true;
}

void clear_component_trial_hints(PbctopoWarmState &warm_state) {
  warm_state.previous_stage1_delta_hints_by_component.clear();
  warm_state.previous_stage2_delta_hints_by_component.clear();
  warm_state.has_component_trial_hints = false;
}

bool component_trial_hints_compatible(const PbctopoWarmState &warm_state,
                                      size_t component_count) {
  return warm_state.has_component_trial_hints &&
         warm_state.previous_stage1_delta_hints_by_component.size() ==
             component_count &&
         warm_state.previous_stage2_delta_hints_by_component.size() ==
             component_count;
}

vector<size_t> build_refine_shift_order(
    const vector<PbctopoShift> &shifts, const Int3 *warm_target_delta,
    const vector<Int3> *stage2_hints, const vector<Int3> *stage1_hints,
    PbctopoTimingStats *timing) {
  vector<size_t> ordered_indices;
  ordered_indices.reserve(shifts.size());
  vector<unsigned char> seen_indices(shifts.size(), 0);

  if (warm_target_delta != nullptr &&
      append_shift_index_if_present(shifts, *warm_target_delta, false,
                                    ordered_indices, seen_indices) &&
      timing != nullptr) {
    ++timing->component_refine_warm_shift_hits;
  }
  auto append_hint_deltas = [&](const vector<Int3> *hint_deltas) {
    if (hint_deltas == nullptr)
      return;
    for (const auto &delta : *hint_deltas) {
      if (append_shift_index_if_present(shifts, delta, false, ordered_indices,
                                        seen_indices) &&
          timing != nullptr) {
        ++timing->component_refine_shortlist_hint_hits;
      }
    }
  };
  append_hint_deltas(stage2_hints);
  append_hint_deltas(stage1_hints);

  for (size_t shift_idx = 0; shift_idx < shifts.size(); ++shift_idx) {
    const auto &shift = shifts[shift_idx];
    if ((shift.ix == 0 && shift.iy == 0 && shift.iz == 0) ||
        seen_indices[shift_idx] != 0) {
      continue;
    }
    seen_indices[shift_idx] = 1;
    ordered_indices.push_back(shift_idx);
  }
  return ordered_indices;
}

vector<size_t> build_pair_shift_order(const vector<PbctopoShift> &shifts,
                                      const Int3 *expected_delta,
                                      PbctopoTimingStats *timing) {
  vector<size_t> ordered_indices;
  ordered_indices.reserve(shifts.size());
  vector<unsigned char> seen_indices(shifts.size(), 0);

  if (expected_delta != nullptr) {
    if (append_shift_index_if_present(shifts, *expected_delta, true,
                                      ordered_indices, seen_indices) &&
        timing != nullptr) {
      ++timing->component_pair_warm_shift_hits;
    }
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          if (dx == 0 && dy == 0 && dz == 0)
            continue;
          const Int3 neighbor = *expected_delta + Int3{dx, dy, dz};
          if (append_shift_index_if_present(shifts, neighbor, true,
                                            ordered_indices, seen_indices) &&
              timing != nullptr) {
            ++timing->component_pair_warm_shift_hits;
          }
        }
      }
    }
  }

  for (size_t shift_idx = 0; shift_idx < shifts.size(); ++shift_idx) {
    if (seen_indices[shift_idx] != 0)
      continue;
    seen_indices[shift_idx] = 1;
    ordered_indices.push_back(shift_idx);
  }
  return ordered_indices;
}

double compute_mst2(const vector<Vec3> &points);
vector<size_t> collect_fragment_options(
    size_t fragment_id, const vector<size_t> &selected_nodes,
    const vector<PbctopoFragmentCopyNode> &nodes,
    const vector<vector<size_t>> &adjacency, size_t max_options);
bool better_cut_candidate(const PbctopoCutCandidate &lhs,
                          const PbctopoCutCandidate &rhs);
bool better_preclash_cut_candidate(const PbctopoCutCandidate &lhs,
                                   const PbctopoCutCandidate &rhs);
vector<pair<int, int>>
build_contacting_owner_pairs(const vector<PbctopoLocalEdge> &local_edges,
                             const vector<PbctopoAtom> &atoms, box Box,
                             double contact_cutoff);
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
bool preclash_cut_candidate_cannot_beat_incumbent(
    const PbctopoCutCandidate &candidate,
    const PbctopoCutCandidate *incumbent);
bool cut_candidate_cannot_beat_incumbent(const PbctopoCutCandidate &partial,
                                         const PbctopoCutCandidate *incumbent);
vector<PbctopoContinuityPair>
build_local_continuity_pairs(const Multi_MOL &system,
                             const PbctopoSystemStaticData &static_data);
void append_global_continuity_pairs(
    const vector<PbctopoContinuityPair> &local_pairs, size_t offset,
    vector<PbctopoContinuityPair> &global_pairs);
bool same_component_layout(const vector<vector<size_t>> &components,
                           const vector<size_t> &sizes);
bool candidate_close_enough_to_skip_search(const PbctopoCutCandidate &candidate,
                                           const PbctopoCutCandidate &reference);
void build_projection_continuity_subset(
    const vector<PbctopoAtom> &base_atoms,
    const vector<PbctopoAtom> &projection_atoms,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    vector<PbctopoContinuityPair> &projection_continuity_pairs,
    vector<size_t> &projection_index);
void build_anchor_index_list(size_t atom_count,
                             const vector<PbctopoContinuityPair> &continuity_pairs,
                             vector<size_t> &anchor_indices);
void append_global_anchor_indices(const vector<size_t> &local_indices, size_t offset,
                                  vector<size_t> &global_indices);
void build_projection_anchor_subset(const vector<PbctopoAtom> &base_atoms,
                                    const vector<size_t> &anchor_indices,
                                    vector<size_t> &projection_anchor_indices,
                                    const vector<size_t> &projection_index);
PbctopoStage1EdgeContext build_stage1_edge_context(
    const vector<PbctopoLocalEdge> &local_edges,
    const vector<size_t> &atom_component_index, size_t moving_component_idx,
    const vector<Vec3> &placed_points);
PbctopoEdgeScanSoA
build_edge_scan_soa(const vector<PbctopoLocalEdge> &local_edges);
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
bool better_stage1_cut_candidate(const PbctopoCutCandidate &lhs,
                                 const PbctopoCutCandidate &rhs);
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
bool anchor_graph_metrics_ok(const vector<Vec3> &placed_points,
                             const vector<size_t> &anchor_indices,
                             double &mst2, double &max_d2,
                             const PbctopoCutCandidate *incumbent = nullptr);
bool anchor_graph_metrics_spatial_or_dense(
    const vector<Vec3> &placed_points, const vector<size_t> &anchor_indices,
    const vector<PbctopoContinuityPair> &continuity_pairs, double &mst2,
    double &max_d2, const PbctopoCutCandidate *incumbent,
    PbctopoRefineScorePhase phase, PbctopoAnchorBackend backend,
    PbctopoTimingStats *timing = nullptr);
void analyze_cut_candidate_ultraclash(const vector<Vec3> &placed_points,
                                      const unordered_set<uint64_t> &edge_set,
                                      PbctopoCutCandidate &candidate,
                                      double cutoff);
vector<double> build_virtual_box_event_planes(const vector<PbctopoAtom> &atoms,
                                              char axis);
uint64_t virtual_box_state_key(const vector<double> &x_planes,
                               const vector<double> &y_planes,
                               const vector<double> &z_planes, double cut_x,
                               double cut_y, double cut_z);
PbctopoAnchorBackend resolve_pbctopo_anchor_backend(const arguments &Arguments);
const char *pbctopo_anchor_backend_name(PbctopoAnchorBackend backend);
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
    const char *reject_context = nullptr);
vector<size_t> build_atom_component_index(
    size_t atom_count, const vector<vector<size_t>> &components);
void append_interface_anchor_edges(
    const vector<PbctopoAtom> &atoms, const vector<Vec3> &placed_points,
    const vector<size_t> &atom_component_index,
    const unordered_set<uint64_t> &edge_set, double contact_cutoff,
    vector<PbctopoLocalEdge> &edges);
PbctopoCutCandidate rescore_cut_candidate(
    const vector<PbctopoLocalEdge> &scoring_edges,
    const vector<Vec3> &placed_points, const PbctopoCutCandidate &source);

[[maybe_unused]] Vec3 point_from_atom_copy(const PbctopoAtomCopy &atom) {
  return {atom.x, atom.y, atom.z};
}

double wrap_unit(double value) {
  value -= floor(value);
  if (value >= 1.0)
    value -= 1.0;
  if (value < 0.0)
    value += 1.0;
  return value;
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

uint64_t atom_pair_key(size_t a, size_t b) {
  if (a > b)
    swap(a, b);
  return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

double cart_dist_sq_from_scaled_delta(box Box, double dsx, double dsy,
                                      double dsz) {
  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);
  const double dx = v1[0] * dsx + v2[0] * dsy + v3[0] * dsz;
  const double dy = v1[1] * dsx + v2[1] * dsy + v3[1] * dsz;
  const double dz = v1[2] * dsx + v2[2] * dsy + v3[2] * dsz;
  return dx * dx + dy * dy + dz * dz;
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

vector<PbctopoRange> build_ranges_from_index(const vector<class index> &Index,
                                             size_t atom_count) {
  vector<PbctopoRange> ranges;
  ranges.reserve(Index.size());
  for (size_t k = 0; k + 1 < Index.size(); ++k) {
    const size_t begin = Index[k].get_res_pos_e(0) - 1;
    const size_t end = Index[k + 1].get_res_pos_e(0) - 1;
    if (begin >= atom_count || end > atom_count || begin >= end)
      continue;
    ranges.push_back({begin, end, ranges.size(), true});
  }
  return ranges;
}

vector<PbctopoRange> build_ranges_from_info(const vector<mol_info> &info,
                                            size_t atom_count) {
  vector<PbctopoRange> ranges;
  if (info.size() != atom_count || atom_count == 0)
    return ranges;

  size_t begin = 0;
  size_t source_index = 0;
  for (size_t i = 1; i < atom_count; ++i) {
    if (info[i].get_chainnum() != info[i - 1].get_chainnum() ||
        info[i].get_chain() != info[i - 1].get_chain()) {
      if (begin < i)
        ranges.push_back({begin, i, source_index++, true});
      begin = i;
    }
  }
  if (begin < atom_count)
    ranges.push_back({begin, atom_count, source_index, true});
  return ranges;
}

vector<PbctopoRange> build_ranges_for_system(const Multi_MOL &system,
                                             const MOL &mol) {
  vector<PbctopoRange> ranges =
      build_ranges_from_index(system.Index, mol.pos.x.size());
  if (!ranges.empty())
    return ranges;

  ranges = build_ranges_from_info(system.info, mol.pos.x.size());
  if (!ranges.empty())
    return ranges;

  if (!mol.pos.x.empty())
    ranges.push_back({0, mol.pos.x.size(), 0, false});
  return ranges;
}

char first_alpha_upper(const string &text) {
  for (unsigned char ch : text) {
    if (isalpha(ch))
      return static_cast<char>(toupper(ch));
  }
  return '\0';
}

bool detect_heavy_atom(const vector<mol_info> *info, size_t local_atom_index) {
  if (info == nullptr || local_atom_index >= info->size())
    return true;
  const auto &atom_info = (*info)[local_atom_index];
  const string atom_name = atom_info.get_atom();
  const string &atom_type = atom_info.get_atom_type();
  const char atom_symbol = first_alpha_upper(atom_name);
  const char type_symbol = first_alpha_upper(atom_type);
  return atom_symbol != 'H' && type_symbol != 'H';
}

void build_heavy_atom_subset(const vector<PbctopoAtom> &atoms,
                             vector<PbctopoAtom> &heavy_atoms) {
  heavy_atoms.clear();
  heavy_atoms.reserve(atoms.size());
  for (const auto &atom : atoms) {
    if (atom.is_heavy)
      heavy_atoms.push_back(atom);
  }
  if (heavy_atoms.empty())
    heavy_atoms = atoms;
}

PbctopoSystemStaticData
build_pbctopo_system_static_data(const Multi_MOL &system, const MOL &mol) {
  PbctopoSystemStaticData data;
  data.atom_count = mol.pos.x.size();
  data.ranges = build_ranges_for_system(system, mol);
  data.heavy_mask.assign(data.atom_count, 1);
  data.heavy_count = 0;
  for (size_t atom_idx = 0; atom_idx < data.atom_count; ++atom_idx) {
    const bool is_heavy = detect_heavy_atom(&system.info, atom_idx);
    data.heavy_mask[atom_idx] = is_heavy ? 1 : 0;
    data.heavy_count += is_heavy ? 1u : 0u;
  }
  if (data.heavy_count == 0 && data.atom_count > 0) {
    std::fill(data.heavy_mask.begin(), data.heavy_mask.end(), 1);
    data.heavy_count = data.atom_count;
  }
  PbctopoSystemStaticData continuity_data = data;
  continuity_data.ranges = build_ranges_from_info(system.info, mol.pos.x.size());
  if (continuity_data.ranges.empty())
    continuity_data.ranges = data.ranges;
  data.continuity_pairs_local =
      build_local_continuity_pairs(system, continuity_data);
  build_anchor_index_list(data.atom_count, data.continuity_pairs_local,
                          data.continuity_anchor_local);
  return data;
}

PbctopoTrajectoryContext
build_pbctopo_trajectory_context(const Multi_MOL &mol1_system,
                                 const Multi_MOL &mol2_system,
                                 size_t frame_count) {
  PbctopoTrajectoryContext context;
  if (frame_count == 0 || mol1_system.m_mol.empty())
    return context;

  const MOL &rec = mol1_system.m_mol.front();
  context.rec_static = build_pbctopo_system_static_data(mol1_system, rec);
  size_t lig_atoms = 0;
  size_t lig_heavy = 0;
  if (!mol2_system.m_mol.empty() && !mol2_system.m_mol.front().pos.x.empty()) {
    const MOL &lig = mol2_system.m_mol.front();
    context.lig_static = build_pbctopo_system_static_data(mol2_system, lig);
    lig_atoms = context.lig_static.atom_count;
    lig_heavy = context.lig_static.heavy_count;
  }

  context.total_atoms_hint =
      max<size_t>(1, context.rec_static.atom_count + lig_atoms);
  context.heavy_atoms_hint =
      max<size_t>(1, context.rec_static.heavy_count + lig_heavy);
  context.local_edges_hint = max<size_t>(1024, context.total_atoms_hint * 8);
  context.projection_edges_hint =
      max<size_t>(1024, context.heavy_atoms_hint * 8);
  append_global_continuity_pairs(context.rec_static.continuity_pairs_local, 0,
                                 context.continuity_pairs);
  append_global_continuity_pairs(context.lig_static.continuity_pairs_local,
                                 context.rec_static.atom_count,
                                 context.continuity_pairs);
  append_global_anchor_indices(context.rec_static.continuity_anchor_local, 0,
                               context.continuity_anchor_indices);
  append_global_anchor_indices(context.lig_static.continuity_anchor_local,
                               context.rec_static.atom_count,
                               context.continuity_anchor_indices);
  context.continuity_pair_hint =
      max<size_t>(64, context.continuity_pairs.size() + 1);
  context.component_hint =
      max<size_t>(8, mol1_system.Index.size() + mol2_system.Index.size());
  context.plane_hint = max<size_t>(32, context.heavy_atoms_hint * 2 + 4);
  context.virtual_state_hint = 512;
  context.coarse_seed_hint = 32;
  context.coarse_candidate_hint = 16;
  return context;
}

void build_projection_edge_subset(const vector<PbctopoAtom> &base_atoms,
                                  const vector<PbctopoAtom> &projection_atoms,
                                  const vector<PbctopoLocalEdge> &local_edges,
                                  vector<PbctopoLocalEdge> &projection_edges,
                                  unordered_set<uint64_t> &projection_edge_set,
                                  vector<size_t> &projection_index) {
  projection_edges.clear();
  projection_edge_set.clear();
  if (projection_atoms.size() == base_atoms.size()) {
    projection_edges = local_edges;
    projection_edge_set.reserve(local_edges.size() * 2 + 1);
    for (const auto &edge : projection_edges)
      projection_edge_set.insert(atom_pair_key(edge.a, edge.b));
    return;
  }

  const size_t invalid_index = numeric_limits<size_t>::max();
  projection_index.assign(base_atoms.size(), invalid_index);
  for (size_t idx = 0; idx < projection_atoms.size(); ++idx)
    projection_index[projection_atoms[idx].global_atom_id] = idx;

  projection_edges.reserve(local_edges.size());
  projection_edge_set.reserve(local_edges.size() + 1);
  for (const auto &edge : local_edges) {
    const size_t mapped_a = projection_index[base_atoms[edge.a].global_atom_id];
    const size_t mapped_b = projection_index[base_atoms[edge.b].global_atom_id];
    if (mapped_a == invalid_index || mapped_b == invalid_index)
      continue;
    if (!projection_edge_set.insert(atom_pair_key(mapped_a, mapped_b)).second)
      continue;
    projection_edges.push_back(
        {mapped_a, mapped_b, edge.ref_d2, edge.inv_ref_d2, edge.break_cutoff_d2});
  }
}

vector<PbctopoContinuityPair>
build_local_continuity_pairs(const Multi_MOL &system,
                             const PbctopoSystemStaticData &static_data) {
  vector<PbctopoContinuityPair> continuity_pairs;
  continuity_pairs.reserve(static_data.atom_count);
  if (system.info.size() == static_data.atom_count && static_data.atom_count > 0) {
    auto atom_priority = [&](size_t atom_idx) {
      const string atom_name = system.info[atom_idx].get_atom();
      if (atom_name == "CA")
        return 0;
      if (atom_name == "C4'" || atom_name == "C4*")
        return 1;
      if (atom_name == "BB")
        return 2;
      if (atom_idx < static_data.heavy_mask.size() && static_data.heavy_mask[atom_idx] != 0)
        return 3;
      return 4;
    };

    size_t prev_anchor = numeric_limits<size_t>::max();
    string prev_chain;
    size_t prev_chainnum = numeric_limits<size_t>::max();
    size_t i = 0;
    while (i < static_data.atom_count) {
      const auto &start_info = system.info[i];
      const size_t chainnum = start_info.get_chainnum();
      const string chain = start_info.get_chain();
      const size_t resnr = start_info.get_resnr();
      size_t j = i + 1;
      while (j < static_data.atom_count &&
             system.info[j].get_chainnum() == chainnum &&
             system.info[j].get_chain() == chain &&
             system.info[j].get_resnr() == resnr) {
        ++j;
      }

      size_t anchor = numeric_limits<size_t>::max();
      int best_priority = numeric_limits<int>::max();
      for (size_t atom_idx = i; atom_idx < j; ++atom_idx) {
        const int priority = atom_priority(atom_idx);
        if (priority < best_priority) {
          best_priority = priority;
          anchor = atom_idx;
        }
      }

      if (anchor != numeric_limits<size_t>::max()) {
        if (prev_anchor != numeric_limits<size_t>::max() &&
            chainnum == prev_chainnum && chain == prev_chain) {
          continuity_pairs.push_back({prev_anchor, anchor});
        }
        prev_anchor = anchor;
        prev_chainnum = chainnum;
        prev_chain = chain;
      }
      i = j;
    }

    if (!continuity_pairs.empty())
      return continuity_pairs;
  }

  for (const auto &range : static_data.ranges) {
    size_t prev = numeric_limits<size_t>::max();
    for (size_t atom_idx = range.begin; atom_idx < range.end; ++atom_idx) {
      if (atom_idx < static_data.heavy_mask.size() &&
          static_data.heavy_mask[atom_idx] == 0) {
        continue;
      }
      if (prev != numeric_limits<size_t>::max())
        continuity_pairs.push_back({prev, atom_idx});
      prev = atom_idx;
    }
  }
  return continuity_pairs;
}

void append_global_continuity_pairs(
    const vector<PbctopoContinuityPair> &local_pairs, size_t offset,
    vector<PbctopoContinuityPair> &global_pairs) {
  global_pairs.reserve(global_pairs.size() + local_pairs.size());
  for (const auto &pair : local_pairs)
    global_pairs.push_back({pair.a + offset, pair.b + offset});
}

void build_projection_continuity_subset(
    const vector<PbctopoAtom> &base_atoms, const vector<PbctopoAtom> &projection_atoms,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    vector<PbctopoContinuityPair> &projection_continuity_pairs,
    vector<size_t> &projection_index) {
  projection_continuity_pairs.clear();
  if (projection_atoms.size() == base_atoms.size()) {
    projection_continuity_pairs = continuity_pairs;
    return;
  }
  const size_t invalid_index = numeric_limits<size_t>::max();
  if (projection_index.size() != base_atoms.size())
    projection_index.assign(base_atoms.size(), invalid_index);
  for (const auto &pair : continuity_pairs) {
    const size_t mapped_a = projection_index[base_atoms[pair.a].global_atom_id];
    const size_t mapped_b = projection_index[base_atoms[pair.b].global_atom_id];
    if (mapped_a == invalid_index || mapped_b == invalid_index)
      continue;
    projection_continuity_pairs.push_back({mapped_a, mapped_b});
  }
}

void build_anchor_index_list(size_t atom_count,
                             const vector<PbctopoContinuityPair> &continuity_pairs,
                             vector<size_t> &anchor_indices) {
  anchor_indices.clear();
  if (atom_count == 0 || continuity_pairs.empty())
    return;
  vector<char> seen(atom_count, 0);
  for (const auto &pair : continuity_pairs) {
    if (pair.a < atom_count && seen[pair.a] == 0) {
      seen[pair.a] = 1;
      anchor_indices.push_back(pair.a);
    }
    if (pair.b < atom_count && seen[pair.b] == 0) {
      seen[pair.b] = 1;
      anchor_indices.push_back(pair.b);
    }
  }
}

void append_global_anchor_indices(const vector<size_t> &local_indices, size_t offset,
                                  vector<size_t> &global_indices) {
  global_indices.reserve(global_indices.size() + local_indices.size());
  for (size_t local_idx : local_indices)
    global_indices.push_back(local_idx + offset);
}

void build_projection_anchor_subset(const vector<PbctopoAtom> &base_atoms,
                                    const vector<size_t> &anchor_indices,
                                    vector<size_t> &projection_anchor_indices,
                                    const vector<size_t> &projection_index) {
  projection_anchor_indices.clear();
  if (anchor_indices.empty())
    return;
  if (projection_index.size() != base_atoms.size()) {
    projection_anchor_indices = anchor_indices;
    return;
  }
  const size_t invalid_index = numeric_limits<size_t>::max();
  projection_anchor_indices.reserve(anchor_indices.size());
  for (size_t anchor_idx : anchor_indices) {
    if (anchor_idx >= base_atoms.size())
      continue;
    const size_t mapped = projection_index[base_atoms[anchor_idx].global_atom_id];
    if (mapped == invalid_index)
      continue;
    projection_anchor_indices.push_back(mapped);
  }
}

vector<size_t> build_sampled_anchor_subset(const vector<PbctopoAtom> &atoms,
                                           const vector<size_t> &anchor_indices,
                                           size_t interior_step,
                                           size_t keep_all_cutoff) {
  if (anchor_indices.empty() || atoms.empty())
    return anchor_indices;
  if (interior_step == 0)
    return anchor_indices;

  vector<size_t> coarse_anchor_indices;
  coarse_anchor_indices.reserve(anchor_indices.size());

  auto append_chain_subset = [&](size_t begin, size_t end) {
    if (begin >= end)
      return;
    const size_t count = end - begin;
    if (count <= keep_all_cutoff) {
      coarse_anchor_indices.insert(coarse_anchor_indices.end(),
                                   anchor_indices.begin() + begin,
                                   anchor_indices.begin() + end);
      return;
    }
    coarse_anchor_indices.push_back(anchor_indices[begin]);
    for (size_t idx = begin + 1; idx + 1 < end; idx += interior_step)
      coarse_anchor_indices.push_back(anchor_indices[idx]);
    coarse_anchor_indices.push_back(anchor_indices[end - 1]);
  };

  size_t chain_begin = 0;
  while (chain_begin < anchor_indices.size()) {
    size_t chain_end = chain_begin + 1;
    const size_t anchor_idx = anchor_indices[chain_begin];
    const uint64_t chain_key =
        (anchor_idx < atoms.size()) ? atoms[anchor_idx].chain_key : 0;
    while (chain_end < anchor_indices.size()) {
      const size_t next_anchor_idx = anchor_indices[chain_end];
      const uint64_t next_chain_key =
          (next_anchor_idx < atoms.size()) ? atoms[next_anchor_idx].chain_key : 0;
      if (next_chain_key != chain_key)
        break;
      ++chain_end;
    }
    append_chain_subset(chain_begin, chain_end);
    chain_begin = chain_end;
  }

  return coarse_anchor_indices;
}

vector<size_t> build_coarse_anchor_subset(const vector<PbctopoAtom> &atoms,
                                          const vector<size_t> &anchor_indices) {
  return build_sampled_anchor_subset(atoms, anchor_indices, 2, 6);
}

vector<size_t> build_third_anchor_subset(const vector<PbctopoAtom> &atoms,
                                         const vector<size_t> &anchor_indices) {
  return build_sampled_anchor_subset(atoms, anchor_indices, 3, 9);
}

vector<size_t> build_tapered_stage1_anchor_subset(
    const vector<PbctopoAtom> &atoms, const vector<size_t> &anchor_indices) {
  if (anchor_indices.empty() || atoms.empty())
    return anchor_indices;

  vector<size_t> stage1_anchor_indices;
  const bool ultra_large_stage1 = anchor_indices.size() >= 1024;
  stage1_anchor_indices.reserve(max<size_t>(
      8, anchor_indices.size() / (ultra_large_stage1 ? 8 : 5)));

  auto append_unique = [&](size_t anchor_idx) {
    if (!stage1_anchor_indices.empty() &&
        stage1_anchor_indices.back() == anchor_idx) {
      return;
    }
    stage1_anchor_indices.push_back(anchor_idx);
  };

  auto append_chain_subset = [&](size_t begin, size_t end) {
    if (begin >= end)
      return;
    const size_t count = end - begin;
    if (count <= 18) {
      stage1_anchor_indices.insert(stage1_anchor_indices.end(),
                                   anchor_indices.begin() + begin,
                                   anchor_indices.begin() + end);
      return;
    }

    const size_t outer_span =
        max<size_t>(4, count / (ultra_large_stage1 ? 5 : 4));
    const size_t middle_begin = min(end, begin + outer_span);
    const size_t middle_end = (end > outer_span) ? max(begin, end - outer_span)
                                                 : begin;
    const size_t outer_step = ultra_large_stage1 ? 6 : 5;
    const size_t middle_step = ultra_large_stage1 ? 12 : 8;

    append_unique(anchor_indices[begin]);
    for (size_t idx = begin + 1; idx < middle_begin; idx += outer_step)
      append_unique(anchor_indices[idx]);
    for (size_t idx = middle_begin; idx < middle_end; idx += middle_step)
      append_unique(anchor_indices[idx]);
    for (size_t idx = middle_end; idx + 1 < end; idx += outer_step)
      append_unique(anchor_indices[idx]);
    append_unique(anchor_indices[end - 1]);

    const size_t quarter = begin + count / 4;
    const size_t halfway = begin + count / 2;
    const size_t three_quarter = begin + (3 * count) / 4;
    if (quarter < end)
      append_unique(anchor_indices[quarter]);
    if (halfway < end)
      append_unique(anchor_indices[halfway]);
    if (three_quarter < end)
      append_unique(anchor_indices[three_quarter]);
  };

  size_t chain_begin = 0;
  while (chain_begin < anchor_indices.size()) {
    size_t chain_end = chain_begin + 1;
    const size_t anchor_idx = anchor_indices[chain_begin];
    const uint64_t chain_key =
        (anchor_idx < atoms.size()) ? atoms[anchor_idx].chain_key : 0;
    while (chain_end < anchor_indices.size()) {
      const size_t next_anchor_idx = anchor_indices[chain_end];
      const uint64_t next_chain_key =
          (next_anchor_idx < atoms.size()) ? atoms[next_anchor_idx].chain_key
                                           : 0;
      if (next_chain_key != chain_key)
        break;
      ++chain_end;
    }
    append_chain_subset(chain_begin, chain_end);
    chain_begin = chain_end;
  }

  return stage1_anchor_indices;
}

void append_atoms_for_system(const MOL &mol, const vector<mol_info> *info,
                             const PbctopoSystemStaticData *static_data,
                             int owner, vector<PbctopoAtom> &atoms) {
  for (size_t local = 0; local < mol.pos.x.size(); ++local) {
    double r[3]{mol.pos.x[local], mol.pos.y[local], mol.pos.z[local]};
    double s[3]{};
    mol.Box.to_scaled(r, s);
    s[0] = wrap_unit(s[0]);
    s[1] = wrap_unit(s[1]);
    s[2] = wrap_unit(s[2]);
    scaled_to_cartesian(mol.Box, s[0], s[1], s[2], r[0], r[1], r[2]);

    PbctopoAtom atom;
    atom.global_atom_id = atoms.size();
    atom.fragment_id = 0;
    atom.owner = owner;
    atom.local_atom_index = local;
    if (info != nullptr && local < info->size()) {
      const auto &atom_info = (*info)[local];
      const uint64_t chainnum =
          static_cast<uint64_t>(atom_info.get_chainnum() & 0xffffffffu);
      const uint64_t chainhash =
          static_cast<uint64_t>(hash<string>{}(atom_info.get_chain()));
      atom.chain_key = (chainnum << 32) ^ (chainhash & 0xffffffffu);
      atom.residue_id = static_cast<int>(atom_info.get_resnr());
      atom.has_topology_key = true;
    }
    atom.is_heavy =
        (static_data != nullptr && local < static_data->heavy_mask.size())
            ? (static_data->heavy_mask[local] != 0)
            : detect_heavy_atom(info, local);
    atom.x = r[0];
    atom.y = r[1];
    atom.z = r[2];
    atom.sx = s[0];
    atom.sy = s[1];
    atom.sz = s[2];
    atoms.push_back(atom);
  }
}

void build_local_neighborhood_graph(
    const vector<PbctopoAtom> &atoms, box Box, vector<PbctopoLocalEdge> &edges,
    unordered_set<uint64_t> &edge_set, double bond_cutoff) {
  edges.clear();
  edge_set.clear();
  if (atoms.size() < 2)
    return;
  const vector<int> residue_sequence_indices =
      build_topology_residue_sequence_indices(atoms);

  const double kBondCutoff = bond_cutoff;
  const double kBondCutoff2 = bond_cutoff * bond_cutoff;
  edges.reserve(atoms.size() * 2);
  edge_set.reserve(atoms.size() * 2);
  const BoxMatrix box_matrix(Box);
  const double sx_bound = min(0.5, kBondCutoff * sqrt(box_matrix.inv[0][0] * box_matrix.inv[0][0] +
                                                      box_matrix.inv[0][1] * box_matrix.inv[0][1] +
                                                      box_matrix.inv[0][2] * box_matrix.inv[0][2]) +
                                     1.0e-9);
  const double sy_bound = min(0.5, kBondCutoff * sqrt(box_matrix.inv[1][0] * box_matrix.inv[1][0] +
                                                      box_matrix.inv[1][1] * box_matrix.inv[1][1] +
                                                      box_matrix.inv[1][2] * box_matrix.inv[1][2]) +
                                     1.0e-9);
  const double sz_bound = min(0.5, kBondCutoff * sqrt(box_matrix.inv[2][0] * box_matrix.inv[2][0] +
                                                      box_matrix.inv[2][1] * box_matrix.inv[2][1] +
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

  unordered_map<ClashCellKey, vector<size_t>, ClashCellKeyHash> bond_grid;
  bond_grid.reserve(max<size_t>(64, atoms.size() / 2));
  for (size_t i = 0; i < atoms.size(); ++i) {
    const auto &ai = atoms[i];
    const int cx = min(nx - 1, static_cast<int>(floor(ai.sx / cell_sx)));
    const int cy = min(ny - 1, static_cast<int>(floor(ai.sy / cell_sy)));
    const int cz = min(nz - 1, static_cast<int>(floor(ai.sz / cell_sz)));
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          const ClashCellKey key{wrap_index(cx + dx, nx), wrap_index(cy + dy, ny),
                                 wrap_index(cz + dz, nz)};
          const auto it = bond_grid.find(key);
          if (it == bond_grid.end())
            continue;
          for (size_t j : it->second) {
            const auto &aj = atoms[j];
            if (ai.owner != aj.owner)
              continue;
            if (ai.has_topology_key && aj.has_topology_key) {
              if (ai.chain_key != aj.chain_key)
                continue;
              const int ai_residue_rank = residue_sequence_indices[i];
              const int aj_residue_rank = residue_sequence_indices[j];
              if (ai_residue_rank != numeric_limits<int>::min() &&
                  aj_residue_rank != numeric_limits<int>::min()) {
                if (abs(ai_residue_rank - aj_residue_rank) > 1)
                  continue;
              } else if (abs(ai.residue_id - aj.residue_id) > 1) {
                continue;
              }
            }
            double dsx = aj.sx - ai.sx;
            double dsy = aj.sy - ai.sy;
            double dsz = aj.sz - ai.sz;
            dsx -= round(dsx);
            dsy -= round(dsy);
            dsz -= round(dsz);
            const double d2 = cart_dist_sq_from_scaled_delta(Box, dsx, dsy, dsz);
            if (d2 > kBondCutoff2)
              continue;
            const uint64_t pair_key = atom_pair_key(i, j);
            if (!edge_set.insert(pair_key).second)
              continue;
            const double ref_d2 = max(d2, 1.0e-6);
            edges.push_back(
                {j, i, d2, 1.0 / ref_d2, max(9.0, ref_d2 * 4.0)});
          }
        }
      }
    }
    bond_grid[{cx, cy, cz}].push_back(i);
  }
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

bool component_bbox_may_clash(const PbctopoShiftedComponent &lhs,
                              const PbctopoShiftedComponent &rhs,
                              double cutoff) {
  const auto gap = [](double min_a, double max_a, double min_b, double max_b) {
    if (max_a < min_b)
      return min_b - max_a;
    if (max_b < min_a)
      return min_a - max_b;
    return 0.0;
  };
  const double dx = gap(lhs.min_x, lhs.max_x, rhs.min_x, rhs.max_x);
  if (dx > cutoff)
    return false;
  const double dy = gap(lhs.min_y, lhs.max_y, rhs.min_y, rhs.max_y);
  if (dy > cutoff)
    return false;
  const double dz = gap(lhs.min_z, lhs.max_z, rhs.min_z, rhs.max_z);
  return dz <= cutoff;
}

PbctopoComponentBaseGeometry
build_component_base_geometry(const vector<PbctopoAtom> &atoms,
                              const vector<size_t> &component) {
  PbctopoComponentBaseGeometry base;
  base.x.resize(component.size());
  base.y.resize(component.size());
  base.z.resize(component.size());
  for (size_t local_idx = 0; local_idx < component.size(); ++local_idx) {
    const size_t atom_idx = component[local_idx];
    const auto &atom = atoms[atom_idx];
    base.x[local_idx] = atom.x;
    base.y[local_idx] = atom.y;
    base.z[local_idx] = atom.z;
    base.center = base.center + Vec3{atom.x, atom.y, atom.z};
    base.min_x = min(base.min_x, atom.x);
    base.min_y = min(base.min_y, atom.y);
    base.min_z = min(base.min_z, atom.z);
    base.max_x = max(base.max_x, atom.x);
    base.max_y = max(base.max_y, atom.y);
    base.max_z = max(base.max_z, atom.z);
  }
  if (!component.empty()) {
    base.center = base.center * (1.0 / static_cast<double>(component.size()));
  }
  return base;
}

PbctopoShiftedComponentDescriptor
describe_shifted_component(const PbctopoComponentBaseGeometry &base,
                           const PbctopoShift &shift) {
  PbctopoShiftedComponentDescriptor descriptor;
  descriptor.weighted_shift_norm2 =
      static_cast<double>(shift.ix * shift.ix + shift.iy * shift.iy +
                          shift.iz * shift.iz) *
      static_cast<double>(base.x.size());
  descriptor.center = {base.center.x + shift.tx, base.center.y + shift.ty,
                       base.center.z + shift.tz};
  descriptor.min_x = base.min_x + shift.tx;
  descriptor.min_y = base.min_y + shift.ty;
  descriptor.min_z = base.min_z + shift.tz;
  descriptor.max_x = base.max_x + shift.tx;
  descriptor.max_y = base.max_y + shift.ty;
  descriptor.max_z = base.max_z + shift.tz;
  return descriptor;
}

void fill_shifted_component_from_base(
    const PbctopoComponentBaseGeometry &base, const PbctopoShift &shift,
    double inv_grid, double inv_contact_grid, bool build_clash_grid,
    PbctopoShiftedComponent &shifted) {
  const PbctopoShiftedComponentDescriptor descriptor =
      describe_shifted_component(base, shift);
  const size_t atom_count = base.x.size();
  shifted.points.resize(atom_count);
  shifted.clash_keys.resize(atom_count);
  shifted.contact_keys.resize(atom_count);
  shifted.clash_grid.clear();
  if (build_clash_grid)
    shifted.clash_grid.reserve(max<size_t>(32, atom_count));
  shifted.contact_grid.clear();
  shifted.contact_grid.reserve(max<size_t>(32, atom_count));
  shifted.weighted_shift_norm2 = descriptor.weighted_shift_norm2;
  shifted.center = descriptor.center;
  shifted.min_x = descriptor.min_x;
  shifted.min_y = descriptor.min_y;
  shifted.min_z = descriptor.min_z;
  shifted.max_x = descriptor.max_x;
  shifted.max_y = descriptor.max_y;
  shifted.max_z = descriptor.max_z;
  for (size_t local_idx = 0; local_idx < atom_count; ++local_idx) {
    const Vec3 point{base.x[local_idx] + shift.tx,
                     base.y[local_idx] + shift.ty,
                     base.z[local_idx] + shift.tz};
    shifted.points[local_idx] = point;
    shifted.clash_keys[local_idx] = {
        static_cast<int>(floor(point.x * inv_grid)),
        static_cast<int>(floor(point.y * inv_grid)),
        static_cast<int>(floor(point.z * inv_grid))};
    if (build_clash_grid)
      shifted.clash_grid[pack_cell_key(shifted.clash_keys[local_idx])]
          .push_back(local_idx);
    const ClashCellKey contact_key{
        static_cast<int>(floor(point.x * inv_contact_grid)),
        static_cast<int>(floor(point.y * inv_contact_grid)),
        static_cast<int>(floor(point.z * inv_contact_grid))};
    shifted.contact_keys[local_idx] = contact_key;
    shifted.contact_grid[pack_cell_key(contact_key)].push_back(local_idx);
  }
}

class PbctopoLazyShiftedComponents {
public:
  PbctopoLazyShiftedComponents(const vector<PbctopoAtom> &atoms,
                               const vector<vector<size_t>> &components,
                               const vector<PbctopoShift> &shifts,
                               double clash_cutoff, double contact_cutoff,
                               PbctopoTimingStats *timing,
                               bool build_clash_grid = false)
      : shifts_(shifts), inv_grid_(1.0 / clash_cutoff),
        inv_contact_grid_(1.0 / contact_cutoff),
        build_clash_grid_(build_clash_grid), timing_(timing) {
    bases_.reserve(components.size());
    descriptors_.resize(components.size());
    materialized_.resize(components.size());
    materialized_flags_.resize(components.size());
    for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
      bases_.push_back(build_component_base_geometry(atoms, components[comp_idx]));
      auto &component_descriptors = descriptors_[comp_idx];
      component_descriptors.reserve(shifts_.size());
      for (const auto &shift : shifts_) {
        component_descriptors.push_back(
            describe_shifted_component(bases_.back(), shift));
      }
      materialized_[comp_idx].resize(shifts_.size());
      materialized_flags_[comp_idx].assign(shifts_.size(), 0);
    }
    if (timing_ != nullptr) {
      timing_->component_graph_shift_descriptors +=
          components.size() * shifts_.size();
    }
  }

  const PbctopoShiftedComponentDescriptor &descriptor(size_t comp_idx,
                                                      size_t shift_idx) const {
    return descriptors_[comp_idx][shift_idx];
  }

  const PbctopoShiftedComponent &materialize(size_t comp_idx, size_t shift_idx) {
    if (materialized_flags_[comp_idx][shift_idx] == 0) {
      const auto started = chrono::steady_clock::now();
      fill_shifted_component_from_base(
          bases_[comp_idx], shifts_[shift_idx], inv_grid_, inv_contact_grid_,
          build_clash_grid_, materialized_[comp_idx][shift_idx]);
      materialized_flags_[comp_idx][shift_idx] = 1;
      if (timing_ != nullptr) {
        ++timing_->component_graph_shift_materialized;
        timing_->component_graph_lazy_materialize_s +=
            chrono::duration<double>(chrono::steady_clock::now() - started)
                .count();
      }
    }
    return materialized_[comp_idx][shift_idx];
  }

private:
  const vector<PbctopoShift> &shifts_;
  vector<PbctopoComponentBaseGeometry> bases_;
  vector<vector<PbctopoShiftedComponentDescriptor>> descriptors_;
  vector<vector<PbctopoShiftedComponent>> materialized_;
  vector<vector<unsigned char>> materialized_flags_;
  double inv_grid_ = 1.0;
  double inv_contact_grid_ = 1.0;
  bool build_clash_grid_ = false;
  PbctopoTimingStats *timing_ = nullptr;
};

vector<vector<PbctopoShiftedComponent>>
precompute_shifted_components(const vector<PbctopoAtom> &atoms,
                              const vector<vector<size_t>> &components,
                              const vector<PbctopoShift> &shifts,
                              double clash_cutoff, double contact_cutoff,
                              bool build_clash_grid = false) {
  vector<vector<PbctopoShiftedComponent>> precomputed(components.size());
  const double inv_grid = 1.0 / clash_cutoff;
  const double inv_contact_grid = 1.0 / contact_cutoff;
  for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
    const PbctopoComponentBaseGeometry base =
        build_component_base_geometry(atoms, components[comp_idx]);
    auto &shifted_list = precomputed[comp_idx];
    shifted_list.resize(shifts.size());
    for (size_t shift_idx = 0; shift_idx < shifts.size(); ++shift_idx) {
      fill_shifted_component_from_base(base, shifts[shift_idx], inv_grid,
                                       inv_contact_grid, build_clash_grid,
                                       shifted_list[shift_idx]);
    }
  }
  return precomputed;
}

double shifted_component_bbox_distance2(const PbctopoShiftedComponent &lhs,
                                        const PbctopoShiftedComponent &rhs) {
  const auto gap = [](double min_a, double max_a, double min_b, double max_b) {
    if (max_a < min_b)
      return min_b - max_a;
    if (max_b < min_a)
      return min_a - max_b;
    return 0.0;
  };
  const double dx = gap(lhs.min_x, lhs.max_x, rhs.min_x, rhs.max_x);
  const double dy = gap(lhs.min_y, lhs.max_y, rhs.min_y, rhs.max_y);
  const double dz = gap(lhs.min_z, lhs.max_z, rhs.min_z, rhs.max_z);
  return dx * dx + dy * dy + dz * dz;
}

double shifted_component_bbox_distance2(
    const PbctopoShiftedComponentDescriptor &lhs,
    const PbctopoShiftedComponentDescriptor &rhs) {
  const auto gap = [](double min_a, double max_a, double min_b, double max_b) {
    if (max_a < min_b)
      return min_b - max_a;
    if (max_b < min_a)
      return min_a - max_b;
    return 0.0;
  };
  const double dx = gap(lhs.min_x, lhs.max_x, rhs.min_x, rhs.max_x);
  const double dy = gap(lhs.min_y, lhs.max_y, rhs.min_y, rhs.max_y);
  const double dz = gap(lhs.min_z, lhs.max_z, rhs.min_z, rhs.max_z);
  return dx * dx + dy * dy + dz * dz;
}

PbctopoPairShiftScore make_bbox_only_pair_shift_score(const PbctopoShift &shift,
                                                       double bbox_d2) {
  PbctopoPairShiftScore score;
  score.shift_index = shift.shift_index;
  score.ix = shift.ix;
  score.iy = shift.iy;
  score.iz = shift.iz;
  score.shift_norm2 = static_cast<double>(shift.ix * shift.ix +
                                          shift.iy * shift.iy +
                                          shift.iz * shift.iz);
  score.bbox_d2 = bbox_d2;
  score.min_pair_d2 = bbox_d2;
  score.valid = true;
  return score;
}

bool better_pair_shift_score(const PbctopoPairShiftScore &lhs,
                             const PbctopoPairShiftScore &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  if (lhs.ultra_clash != rhs.ultra_clash)
    return rhs.ultra_clash;
  if (lhs.contact_atoms != rhs.contact_atoms)
    return lhs.contact_atoms > rhs.contact_atoms;
  if (fabs(lhs.soft_penalty - rhs.soft_penalty) > 1.0e-9)
    return lhs.soft_penalty < rhs.soft_penalty;
  if (fabs(lhs.min_pair_d2 - rhs.min_pair_d2) > 1.0e-9)
    return lhs.min_pair_d2 < rhs.min_pair_d2;
  if (fabs(lhs.bbox_d2 - rhs.bbox_d2) > 1.0e-9)
    return lhs.bbox_d2 < rhs.bbox_d2;
  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-9)
    return lhs.shift_norm2 < rhs.shift_norm2;
  if (lhs.ix != rhs.ix)
    return lhs.ix < rhs.ix;
  if (lhs.iy != rhs.iy)
    return lhs.iy < rhs.iy;
  return lhs.iz < rhs.iz;
}

bool pair_shift_cannot_beat_incumbent(
    const PbctopoPairShiftScore &partial, size_t remaining_lhs_points,
    const PbctopoPairShiftScore *incumbent) {
  if (incumbent == nullptr || !incumbent->valid || incumbent->ultra_clash ||
      partial.ultra_clash) {
    return false;
  }
  const size_t max_contact_atoms = partial.contact_atoms + remaining_lhs_points;
  if (max_contact_atoms < incumbent->contact_atoms)
    return true;
  if (max_contact_atoms == incumbent->contact_atoms &&
      partial.soft_penalty > incumbent->soft_penalty + 1.0e-9) {
    return true;
  }
  return false;
}

bool better_component_graph_edge(const PbctopoGraphEdge &lhs,
                                 const PbctopoGraphEdge &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  if (lhs.ultra_clash != rhs.ultra_clash)
    return rhs.ultra_clash;
  if (lhs.contact_atoms != rhs.contact_atoms)
    return lhs.contact_atoms > rhs.contact_atoms;
  if (fabs(lhs.confidence - rhs.confidence) > 1.0e-9)
    return lhs.confidence > rhs.confidence;
  if (fabs(lhs.soft_penalty - rhs.soft_penalty) > 1.0e-9)
    return lhs.soft_penalty < rhs.soft_penalty;
  if (fabs(lhs.min_pair_d2 - rhs.min_pair_d2) > 1.0e-9)
    return lhs.min_pair_d2 < rhs.min_pair_d2;
  if (fabs(lhs.bbox_d2 - rhs.bbox_d2) > 1.0e-9)
    return lhs.bbox_d2 < rhs.bbox_d2;
  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-9)
    return lhs.shift_norm2 < rhs.shift_norm2;
  if (lhs.a != rhs.a)
    return lhs.a < rhs.a;
  return lhs.b < rhs.b;
}

PbctopoPairShiftScore score_component_pair_shift(
    const PbctopoShiftedComponent &lhs, const PbctopoShiftedComponent &rhs,
    const PbctopoShift &shift, double hard_cutoff, double contact_cutoff,
    const PbctopoPairShiftScore *incumbent = nullptr) {
  PbctopoPairShiftScore score;
  score.shift_index = shift.shift_index;
  score.ix = shift.ix;
  score.iy = shift.iy;
  score.iz = shift.iz;
  score.shift_norm2 =
      static_cast<double>(shift.ix * shift.ix + shift.iy * shift.iy + shift.iz * shift.iz);
  score.bbox_d2 = shifted_component_bbox_distance2(lhs, rhs);
  score.min_pair_d2 = score.bbox_d2;
  score.valid = true;

  const double contact_cutoff2 = contact_cutoff * contact_cutoff;
  const double hard_cutoff2 = hard_cutoff * hard_cutoff;
  constexpr double kSoftMin2 = 2.6 * 2.6;
  if (score.bbox_d2 > contact_cutoff2)
    return score;

  for (size_t lhs_idx = 0; lhs_idx < lhs.points.size(); ++lhs_idx) {
    const auto &lhs_point = lhs.points[lhs_idx];
    const ClashCellKey &contact_key = lhs.contact_keys[lhs_idx];
    const int cx = contact_key.x;
    const int cy = contact_key.y;
    const int cz = contact_key.z;
    bool atom_has_contact = false;
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dz = -1; dz <= 1; ++dz) {
            const auto it =
                rhs.contact_grid.find(pack_cell_key(cx + dx, cy + dy, cz + dz));
            if (it == rhs.contact_grid.end())
              continue;
            for (size_t rhs_idx : it->second) {
              const auto &rhs_point = rhs.points[rhs_idx];
              const double ddx = lhs_point.x - rhs_point.x;
              const double ddy = lhs_point.y - rhs_point.y;
              const double ddz = lhs_point.z - rhs_point.z;
              const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
              if (d2 < score.min_pair_d2)
                score.min_pair_d2 = d2;
              if (d2 + kDistanceEps < hard_cutoff2) {
                score.ultra_clash = true;
                score.quality_score = -numeric_limits<double>::infinity();
                return score;
              }
              if (d2 < kSoftMin2)
                score.soft_penalty += (kSoftMin2 - d2);
              if (d2 <= contact_cutoff2)
                atom_has_contact = true;
            }
          }
        }
      }
    if (atom_has_contact)
      ++score.contact_atoms;
    const size_t remaining_lhs_points = lhs.points.size() - lhs_idx - 1;
    if (pair_shift_cannot_beat_incumbent(score, remaining_lhs_points, incumbent))
      return score;
  }

  score.quality_score = static_cast<double>(score.contact_atoms) * 1000.0 -
                        score.soft_penalty * 50.0 - score.bbox_d2 * 5.0 -
                        score.shift_norm2;
  return score;
}

struct ComponentGraphDsu {
  vector<size_t> parent;
  vector<size_t> rank;

  explicit ComponentGraphDsu(size_t n) : parent(n), rank(n, 0) {
    iota(parent.begin(), parent.end(), 0);
  }

  size_t find(size_t x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  }

  bool unite(size_t a, size_t b) {
    a = find(a);
    b = find(b);
    if (a == b)
      return false;
    if (rank[a] < rank[b])
      swap(a, b);
    parent[b] = a;
    if (rank[a] == rank[b])
      ++rank[a];
    return true;
  }
};

void rebase_component_offsets(vector<Int3> &offsets,
                              const vector<vector<size_t>> &components) {
  if (offsets.empty() || components.empty())
    return;
  auto choose_axis = [&](char axis) {
    double weighted_sum = 0.0;
    double total_weight = 0.0;
    for (size_t idx = 0; idx < offsets.size(); ++idx) {
      const double weight = static_cast<double>(components[idx].size());
      const int value =
          (axis == 'x') ? offsets[idx].x : ((axis == 'y') ? offsets[idx].y : offsets[idx].z);
      weighted_sum += weight * static_cast<double>(value);
      total_weight += weight;
    }
    if (total_weight <= 0.0)
      return 0;
    return static_cast<int>(llround(weighted_sum / total_weight));
  };

  const Int3 global_shift{choose_axis('x'), choose_axis('y'), choose_axis('z')};
  for (auto &offset : offsets)
    offset = offset - global_shift;
}

void build_component_points_from_offsets(
    const vector<PbctopoAtom> &atoms, const vector<vector<size_t>> &components,
    box Box, const vector<Int3> &offsets, vector<Vec3> &placed_points,
    double &shift_norm2) {
  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);
  placed_points.resize(atoms.size());
  shift_norm2 = 0.0;
  for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
    const auto &offset = offsets[comp_idx];
    const Vec3 translation{
        offset.x * v1[0] + offset.y * v2[0] + offset.z * v3[0],
        offset.x * v1[1] + offset.y * v2[1] + offset.z * v3[1],
        offset.x * v1[2] + offset.y * v2[2] + offset.z * v3[2]};
    const double comp_shift_norm2 =
        static_cast<double>(offset.x * offset.x + offset.y * offset.y + offset.z * offset.z);
    shift_norm2 += comp_shift_norm2 * static_cast<double>(components[comp_idx].size());
    for (size_t atom_idx : components[comp_idx]) {
      placed_points[atom_idx] = {atoms[atom_idx].x + translation.x,
                                 atoms[atom_idx].y + translation.y,
                                 atoms[atom_idx].z + translation.z};
    }
  }
}

struct PbctopoBoxBasis {
  Vec3 v1{};
  Vec3 v2{};
  Vec3 v3{};
};

PbctopoBoxBasis build_box_basis(box Box) {
  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);
  return {{v1[0], v1[1], v1[2]},
          {v2[0], v2[1], v2[2]},
          {v3[0], v3[1], v3[2]}};
}

Vec3 component_translation(const PbctopoBoxBasis &basis, const Int3 &offset) {
  return {offset.x * basis.v1.x + offset.y * basis.v2.x +
              offset.z * basis.v3.x,
          offset.x * basis.v1.y + offset.y * basis.v2.y +
              offset.z * basis.v3.y,
          offset.x * basis.v1.z + offset.y * basis.v2.z +
              offset.z * basis.v3.z};
}

PointStats translate_point_stats(const PointStats &base, const Vec3 &delta) {
  if (base.count == 0)
    return {};
  PointStats translated = base;
  const double n = static_cast<double>(base.count);
  translated.sum_x = base.sum_x + n * delta.x;
  translated.sum_y = base.sum_y + n * delta.y;
  translated.sum_z = base.sum_z + n * delta.z;
  translated.sum_x2 =
      base.sum_x2 + 2.0 * delta.x * base.sum_x + n * delta.x * delta.x;
  translated.sum_y2 =
      base.sum_y2 + 2.0 * delta.y * base.sum_y + n * delta.y * delta.y;
  translated.sum_z2 =
      base.sum_z2 + 2.0 * delta.z * base.sum_z + n * delta.z * delta.z;
  translated.sum_xy = base.sum_xy + delta.x * base.sum_y +
                      delta.y * base.sum_x + n * delta.x * delta.y;
  translated.sum_xz = base.sum_xz + delta.x * base.sum_z +
                      delta.z * base.sum_x + n * delta.x * delta.z;
  translated.sum_yz = base.sum_yz + delta.y * base.sum_z +
                      delta.z * base.sum_y + n * delta.y * delta.z;
  translated.sum_sq =
      translated.sum_x2 + translated.sum_y2 + translated.sum_z2;
  translated.min_x = base.min_x + delta.x;
  translated.min_y = base.min_y + delta.y;
  translated.min_z = base.min_z + delta.z;
  translated.max_x = base.max_x + delta.x;
  translated.max_y = base.max_y + delta.y;
  translated.max_z = base.max_z + delta.z;
  return translated;
}

void merge_point_stats(PointStats &dst, const PointStats &src) {
  if (src.count == 0)
    return;
  if (dst.count == 0) {
    dst = src;
    return;
  }
  dst.count += src.count;
  dst.sum_x += src.sum_x;
  dst.sum_y += src.sum_y;
  dst.sum_z += src.sum_z;
  dst.sum_sq += src.sum_sq;
  dst.sum_x2 += src.sum_x2;
  dst.sum_y2 += src.sum_y2;
  dst.sum_z2 += src.sum_z2;
  dst.sum_xy += src.sum_xy;
  dst.sum_xz += src.sum_xz;
  dst.sum_yz += src.sum_yz;
  dst.min_x = min(dst.min_x, src.min_x);
  dst.min_y = min(dst.min_y, src.min_y);
  dst.min_z = min(dst.min_z, src.min_z);
  dst.max_x = max(dst.max_x, src.max_x);
  dst.max_y = max(dst.max_y, src.max_y);
  dst.max_z = max(dst.max_z, src.max_z);
}

double compute_shift_norm2_from_offsets(
    const vector<vector<size_t>> &components, const vector<Int3> &offsets) {
  double shift_norm2 = 0.0;
  const size_t count = min(components.size(), offsets.size());
  for (size_t comp_idx = 0; comp_idx < count; ++comp_idx) {
    const auto &offset = offsets[comp_idx];
    const double comp_shift_norm2 =
        static_cast<double>(offset.x * offset.x + offset.y * offset.y +
                            offset.z * offset.z);
    shift_norm2 +=
        comp_shift_norm2 * static_cast<double>(components[comp_idx].size());
  }
  return shift_norm2;
}

void accumulate_component_graph_refine_build_points(
    PbctopoTimingStats *timing, PbctopoRefineScorePhase phase, double delta) {
  if (timing == nullptr)
    return;
  switch (phase) {
  case PbctopoRefineScorePhase::Exact:
    timing->component_graph_refine_build_points_s += delta;
    break;
  case PbctopoRefineScorePhase::CoarseStage1:
    timing->component_graph_refine_coarse_s += delta;
    timing->component_graph_refine_coarse_stage1_s += delta;
    timing->component_graph_refine_coarse_stage1_build_points_s += delta;
    break;
  case PbctopoRefineScorePhase::CoarseStage2:
    timing->component_graph_refine_coarse_s += delta;
    timing->component_graph_refine_coarse_stage2_s += delta;
    timing->component_graph_refine_coarse_stage2_build_points_s += delta;
    break;
  }
}

struct PbctopoRefineScratchState {
  vector<Vec3> *points = nullptr;
  size_t generation = 0;
  size_t component = numeric_limits<size_t>::max();
  bool valid = false;
};

struct PbctopoRefinePlacementWorkspace {
  const vector<PbctopoAtom> *atoms = nullptr;
  const vector<vector<size_t>> *components = nullptr;
  box source_box{};
  PbctopoBoxBasis basis{};
  vector<PointStats> base_component_stats;
  vector<Int3> baseline_offsets;
  vector<Vec3> baseline_points;
  vector<PbctopoRefineScratchState> scratch_states;
  size_t generation = 0;
  bool initialized = false;

  void initialize(const vector<PbctopoAtom> &atoms_ref,
                  const vector<vector<size_t>> &components_ref, box Box) {
    atoms = &atoms_ref;
    components = &components_ref;
    source_box = Box;
    basis = build_box_basis(Box);
    base_component_stats.assign(components_ref.size(), {});
    for (size_t comp_idx = 0; comp_idx < components_ref.size(); ++comp_idx) {
      PointStats stats;
      for (size_t atom_idx : components_ref[comp_idx]) {
        if (atom_idx >= atoms_ref.size())
          continue;
        const auto &atom = atoms_ref[atom_idx];
        stats.add(atom.x, atom.y, atom.z);
      }
      base_component_stats[comp_idx] = stats;
    }
    initialized = true;
  }

  void reset_baseline(const vector<Int3> &offsets,
                      const vector<Vec3> &placed_points) {
    vector<size_t> changed_components;
    const bool can_update_scratch =
        (initialized && atoms != nullptr && components != nullptr &&
         placed_points.size() == atoms->size() &&
         baseline_points.size() == atoms->size() &&
         baseline_offsets.size() == offsets.size() &&
         offsets.size() == components->size());
    if (can_update_scratch) {
      for (size_t comp_idx = 0; comp_idx < offsets.size(); ++comp_idx) {
        if (!same_delta(offsets[comp_idx], baseline_offsets[comp_idx]))
          changed_components.push_back(comp_idx);
      }
      for (auto &state : scratch_states) {
        if (!state.valid || state.points == nullptr ||
            state.generation != generation ||
            state.points->size() != atoms->size()) {
          state.valid = false;
          continue;
        }
        if (state.component != numeric_limits<size_t>::max()) {
          if (state.component >= components->size()) {
            state.valid = false;
            continue;
          }
          for (size_t atom_idx : (*components)[state.component]) {
            if (atom_idx < baseline_points.size() &&
                atom_idx < state.points->size()) {
              (*state.points)[atom_idx] = baseline_points[atom_idx];
            }
          }
        }
        for (size_t comp_idx : changed_components) {
          for (size_t atom_idx : (*components)[comp_idx]) {
            if (atom_idx < placed_points.size() &&
                atom_idx < state.points->size()) {
              (*state.points)[atom_idx] = placed_points[atom_idx];
            }
          }
        }
        state.component = numeric_limits<size_t>::max();
      }
    } else {
      for (auto &state : scratch_states)
        state.valid = false;
    }
    baseline_offsets = offsets;
    baseline_points = placed_points;
    ++generation;
    if (can_update_scratch) {
      for (auto &state : scratch_states) {
        if (state.valid)
          state.generation = generation;
      }
    }
  }

  PbctopoRefineScratchState &scratch_state(vector<Vec3> &points) {
    for (auto &state : scratch_states) {
      if (state.points == &points)
        return state;
    }
    scratch_states.push_back({&points, 0, numeric_limits<size_t>::max(), false});
    return scratch_states.back();
  }

  bool single_component_delta(size_t comp_idx, const vector<Int3> &offsets,
                              bool &global_rebase) const {
    global_rebase = false;
    if (!initialized || components == nullptr || atoms == nullptr ||
        offsets.size() != baseline_offsets.size() ||
        baseline_offsets.size() != components->size() ||
        baseline_points.size() != atoms->size() || comp_idx >= components->size()) {
      global_rebase = true;
      return false;
    }
    size_t changed_count = 0;
    size_t changed_idx = numeric_limits<size_t>::max();
    for (size_t idx = 0; idx < offsets.size(); ++idx) {
      if (same_delta(offsets[idx], baseline_offsets[idx]))
        continue;
      ++changed_count;
      changed_idx = idx;
      if (changed_count > 1) {
        global_rebase = true;
        return false;
      }
    }
    if (changed_count == 0)
      return true;
    if (changed_idx == comp_idx)
      return true;
    global_rebase = true;
    return false;
  }

  PointStats stats_from_offsets(const vector<Int3> &offsets) const {
    PointStats merged;
    if (components == nullptr)
      return merged;
    const size_t count =
        min({components->size(), offsets.size(), base_component_stats.size()});
    for (size_t comp_idx = 0; comp_idx < count; ++comp_idx) {
      const Vec3 translation = component_translation(basis, offsets[comp_idx]);
      merge_point_stats(
          merged, translate_point_stats(base_component_stats[comp_idx],
                                        translation));
    }
    return merged;
  }

  bool prepare_trial(size_t comp_idx, const vector<Int3> &offsets,
                     vector<Vec3> &out_points, PointStats &out_stats,
                     double &shift_norm2, PbctopoTimingStats *timing,
                     PbctopoRefineScorePhase phase) {
    bool global_rebase = false;
    const bool incremental = single_component_delta(comp_idx, offsets,
                                                    global_rebase);
    auto build_started = chrono::steady_clock::now();
    if (incremental) {
      auto &state = scratch_state(out_points);
      if (!state.valid || state.generation != generation ||
          (state.component != comp_idx &&
           state.component != numeric_limits<size_t>::max())) {
        out_points = baseline_points;
        state.valid = true;
        state.generation = generation;
      }
      state.component = comp_idx;
      const Vec3 translation = component_translation(basis, offsets[comp_idx]);
      for (size_t atom_idx : (*components)[comp_idx]) {
        if (atom_idx >= atoms->size() || atom_idx >= out_points.size())
          continue;
        const auto &atom = (*atoms)[atom_idx];
        out_points[atom_idx] = {atom.x + translation.x, atom.y + translation.y,
                                atom.z + translation.z};
      }
      if (timing != nullptr) {
        ++timing->component_graph_refine_trial_build_incremental;
        accumulate_component_graph_refine_build_points(
            timing, phase,
            chrono::duration<double>(chrono::steady_clock::now() -
                                     build_started)
                .count());
      }
      auto stats_started = chrono::steady_clock::now();
      out_stats = stats_from_offsets(offsets);
      shift_norm2 = compute_shift_norm2_from_offsets(*components, offsets);
      if (timing != nullptr) {
        ++timing->component_graph_refine_stats_cached;
        timing->component_graph_refine_stats_s +=
            chrono::duration<double>(chrono::steady_clock::now() -
                                     stats_started)
                .count();
      }
      return true;
    }

    if (timing != nullptr && global_rebase)
      ++timing->component_graph_refine_trial_build_global_rebase;
    build_component_points_from_offsets(*atoms, *components, source_box, offsets,
                                        out_points, shift_norm2);
    scratch_state(out_points).valid = false;
    if (timing != nullptr) {
      ++timing->component_graph_refine_trial_build_full;
      accumulate_component_graph_refine_build_points(
          timing, phase,
          chrono::duration<double>(chrono::steady_clock::now() - build_started)
              .count());
    }
    auto stats_started = chrono::steady_clock::now();
    out_stats = build_point_stats(out_points);
    if (timing != nullptr) {
      ++timing->component_graph_refine_stats_full;
      timing->component_graph_refine_stats_s +=
          chrono::duration<double>(chrono::steady_clock::now() -
                                   stats_started)
              .count();
    }
    return false;
  }
};

bool scaled_span_fits_unit_box(box Box, const vector<Vec3> &placed_points,
                               double tolerance) {
  if (placed_points.empty())
    return false;

  const BoxMatrix box_matrix(Box);
  double min_sx = numeric_limits<double>::infinity();
  double min_sy = numeric_limits<double>::infinity();
  double min_sz = numeric_limits<double>::infinity();
  double max_sx = -numeric_limits<double>::infinity();
  double max_sy = -numeric_limits<double>::infinity();
  double max_sz = -numeric_limits<double>::infinity();

  for (const auto &point : placed_points) {
    const double sx = box_matrix.inv[0][0] * point.x +
                      box_matrix.inv[0][1] * point.y +
                      box_matrix.inv[0][2] * point.z;
    const double sy = box_matrix.inv[1][0] * point.x +
                      box_matrix.inv[1][1] * point.y +
                      box_matrix.inv[1][2] * point.z;
    const double sz = box_matrix.inv[2][0] * point.x +
                      box_matrix.inv[2][1] * point.y +
                      box_matrix.inv[2][2] * point.z;
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
                             double contact_cutoff) {
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
    register_owner_pair(atoms[edge.a].owner, atoms[edge.b].owner);
  }

  if (atoms.size() < 2 || contact_cutoff <= kDistanceEps)
    return owner_pairs;

  const BoxMatrix box_matrix(Box);
  const double cutoff2 = contact_cutoff * contact_cutoff;
  const double sx_bound =
      min(0.5, contact_cutoff * sqrt(box_matrix.inv[0][0] * box_matrix.inv[0][0] +
                                     box_matrix.inv[0][1] * box_matrix.inv[0][1] +
                                     box_matrix.inv[0][2] * box_matrix.inv[0][2]) +
                    1.0e-9);
  const double sy_bound =
      min(0.5, contact_cutoff * sqrt(box_matrix.inv[1][0] * box_matrix.inv[1][0] +
                                     box_matrix.inv[1][1] * box_matrix.inv[1][1] +
                                     box_matrix.inv[1][2] * box_matrix.inv[1][2]) +
                    1.0e-9);
  const double sz_bound =
      min(0.5, contact_cutoff * sqrt(box_matrix.inv[2][0] * box_matrix.inv[2][0] +
                                     box_matrix.inv[2][1] * box_matrix.inv[2][1] +
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
            double dx_cart = other.x - atom.x;
            double dy_cart = other.y - atom.y;
            double dz_cart = other.z - atom.z;
            box_matrix.apply_mic(dx_cart, dy_cart, dz_cart);
            const double d2 = dx_cart * dx_cart + dy_cart * dy_cart +
                              dz_cart * dz_cart;
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
  const BoxMatrix box_matrix(Box);
  for (const auto &[owner, point] : owner_com_cartesian) {
    owner_com_scaled[owner] = {
        box_matrix.inv[0][0] * point.x + box_matrix.inv[0][1] * point.y +
            box_matrix.inv[0][2] * point.z,
        box_matrix.inv[1][0] * point.x + box_matrix.inv[1][1] * point.y +
            box_matrix.inv[1][2] * point.z,
        box_matrix.inv[2][0] * point.x + box_matrix.inv[2][1] * point.y +
            box_matrix.inv[2][2] * point.z};
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

  for (const auto &[owner_a, owner_b] : contacting_owner_pairs) {
    const auto it_a = owner_com_scaled.find(owner_a);
    const auto it_b = owner_com_scaled.find(owner_b);
    if (it_a == owner_com_scaled.end() || it_b == owner_com_scaled.end())
      continue;

    const Vec3 raw_delta_scaled{it_b->second.x - it_a->second.x,
                                it_b->second.y - it_a->second.y,
                                it_b->second.z - it_a->second.z};
    const Vec3 minimg_delta_scaled{
        raw_delta_scaled.x - round(raw_delta_scaled.x),
        raw_delta_scaled.y - round(raw_delta_scaled.y),
        raw_delta_scaled.z - round(raw_delta_scaled.z)};
    const Vec3 split_delta{raw_delta_scaled.x - minimg_delta_scaled.x,
                           raw_delta_scaled.y - minimg_delta_scaled.y,
                           raw_delta_scaled.z - minimg_delta_scaled.z};
    const double axis_shift_max =
        max({fabs(split_delta.x), fabs(split_delta.y), fabs(split_delta.z)});
    const double raw_norm = sqrt(raw_delta_scaled.x * raw_delta_scaled.x +
                                 raw_delta_scaled.y * raw_delta_scaled.y +
                                 raw_delta_scaled.z * raw_delta_scaled.z);
    const double minimg_norm =
        sqrt(minimg_delta_scaled.x * minimg_delta_scaled.x +
             minimg_delta_scaled.y * minimg_delta_scaled.y +
             minimg_delta_scaled.z * minimg_delta_scaled.z);
    const double norm_gap = raw_norm - minimg_norm;
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
        detail->raw_norm = raw_norm;
        detail->minimg_norm = minimg_norm;
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
        const Int3 proposed_shift{
            shift_u.x - static_cast<int>(llround(round(com_v.x - com_u.x))),
            shift_u.y - static_cast<int>(llround(round(com_v.y - com_u.y))),
            shift_u.z - static_cast<int>(llround(round(com_v.z - com_u.z)))};        
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

inline void relax_anchor_best_scalar(double *best, const double *px,
                                     const double *py, const double *pz,
                                     size_t n, double sx, double sy,
                                     double sz) noexcept {
  for (size_t i = 0; i < n; ++i) {
    const double dx = sx - px[i];
    const double dy = sy - py[i];
    const double dz = sz - pz[i];
    const double d2 = dx * dx + dy * dy + dz * dz;
    best[i] = min(best[i], d2);
  }
}

#if TITAN_SIMD_X86
TITAN_TARGET_AVX2
inline void relax_anchor_best_avx2(double *best, const double *px,
                                   const double *py, const double *pz,
                                   size_t n, double sx, double sy,
                                   double sz) noexcept {
  const __m256d vsx = _mm256_set1_pd(sx);
  const __m256d vsy = _mm256_set1_pd(sy);
  const __m256d vsz = _mm256_set1_pd(sz);
  size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    const __m256d dx = _mm256_sub_pd(vsx, _mm256_loadu_pd(px + i));
    const __m256d dy = _mm256_sub_pd(vsy, _mm256_loadu_pd(py + i));
    const __m256d dz = _mm256_sub_pd(vsz, _mm256_loadu_pd(pz + i));
    const __m256d d2 = _mm256_add_pd(
        _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)),
        _mm256_mul_pd(dz, dz));
    const __m256d best_vec = _mm256_loadu_pd(best + i);
    _mm256_storeu_pd(best + i, _mm256_min_pd(best_vec, d2));
  }
  for (; i < n; ++i) {
    const double dx = sx - px[i];
    const double dy = sy - py[i];
    const double dz = sz - pz[i];
    const double d2 = dx * dx + dy * dy + dz * dz;
    best[i] = min(best[i], d2);
  }
}

TITAN_TARGET_AVX2_FMA
inline void relax_anchor_best_avx2_fma(double *best, const double *px,
                                       const double *py, const double *pz,
                                       size_t n, double sx, double sy,
                                       double sz) noexcept {
  const __m256d vsx = _mm256_set1_pd(sx);
  const __m256d vsy = _mm256_set1_pd(sy);
  const __m256d vsz = _mm256_set1_pd(sz);
  size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    const __m256d dx = _mm256_sub_pd(vsx, _mm256_loadu_pd(px + i));
    const __m256d dy = _mm256_sub_pd(vsy, _mm256_loadu_pd(py + i));
    const __m256d dz = _mm256_sub_pd(vsz, _mm256_loadu_pd(pz + i));
    const __m256d d2 =
        _mm256_fmadd_pd(dz, dz,
                        _mm256_fmadd_pd(dy, dy, _mm256_mul_pd(dx, dx)));
    const __m256d best_vec = _mm256_loadu_pd(best + i);
    _mm256_storeu_pd(best + i, _mm256_min_pd(best_vec, d2));
  }
  for (; i < n; ++i) {
    const double dx = sx - px[i];
    const double dy = sy - py[i];
    const double dz = sz - pz[i];
    const double d2 = dx * dx + dy * dy + dz * dz;
    best[i] = min(best[i], d2);
  }
}
#endif

inline void relax_anchor_best(double *best, const double *px, const double *py,
                              const double *pz, size_t n, double sx,
                              double sy, double sz) noexcept {
  using Kernel = void (*)(double *, const double *, const double *,
                          const double *, size_t, double, double,
                          double) noexcept;
  static const Kernel kernel = []() noexcept -> Kernel {
#if TITAN_SIMD_X86 && TITAN_RUNTIME_TARGET_ATTRS
    if (titan_simd_dispatch::runtime_has_avx2_fma())
      return &relax_anchor_best_avx2_fma;
    if (titan_simd_dispatch::runtime_has_avx2())
      return &relax_anchor_best_avx2;
#endif
    return &relax_anchor_best_scalar;
  }();
  kernel(best, px, py, pz, n, sx, sy, sz);
}

inline size_t select_anchor_next_scalar(const double *best,
                                        const unsigned char *used, size_t n,
                                        double &next_cost) noexcept {
  size_t next = numeric_limits<size_t>::max();
  next_cost = numeric_limits<double>::infinity();
  for (size_t i = 0; i < n; ++i) {
    if (used[i] != 0)
      continue;
    if (best[i] < next_cost) {
      next = i;
      next_cost = best[i];
    }
  }
  return next;
}

#if TITAN_SIMD_X86
TITAN_TARGET_AVX2
inline size_t select_anchor_next_avx2(const double *best,
                                      const unsigned char *used, size_t n,
                                      double &next_cost) noexcept {
  const double inf = numeric_limits<double>::infinity();
  size_t next = numeric_limits<size_t>::max();
  next_cost = inf;
  alignas(32) double lane_best[4];
  size_t i = 0;
  for (; i + 4 <= n; i += 4) {
    const __m256d best_vec = _mm256_loadu_pd(best + i);
    _mm256_store_pd(lane_best, best_vec);
    for (size_t lane = 0; lane < 4; ++lane) {
      if (used[i + lane] != 0)
        continue;
      if (lane_best[lane] < next_cost) {
        next = i + lane;
        next_cost = lane_best[lane];
      }
    }
  }
  for (; i < n; ++i) {
    if (used[i] != 0)
      continue;
    if (best[i] < next_cost) {
      next = i;
      next_cost = best[i];
    }
  }
  return next;
}
#endif

inline size_t select_anchor_next(const double *best, const unsigned char *used,
                                 size_t n, double &next_cost) noexcept {
  using Kernel =
      size_t (*)(const double *, const unsigned char *, size_t, double &) noexcept;
  static const Kernel kernel = []() noexcept -> Kernel {
#if TITAN_SIMD_X86 && TITAN_RUNTIME_TARGET_ATTRS
    if (titan_simd_dispatch::runtime_has_avx2())
      return &select_anchor_next_avx2;
#endif
    return &select_anchor_next_scalar;
  }();
  return kernel(best, used, n, next_cost);
}

bool anchor_graph_metrics_ok(const vector<Vec3> &placed_points,
                             const vector<size_t> &anchor_indices,
                             double &mst2, double &max_d2,
                             const PbctopoCutCandidate *incumbent) {
  mst2 = 0.0;
  max_d2 = 0.0;
  if (anchor_indices.size() < 2)
    return true;

  const size_t n = anchor_indices.size();
  thread_local vector<double> best_storage;
  thread_local vector<unsigned char> used_storage;
  thread_local vector<double> px_storage;
  thread_local vector<double> py_storage;
  thread_local vector<double> pz_storage;

  if (best_storage.size() < n)
    best_storage.resize(n);
  if (used_storage.size() < n)
    used_storage.resize(n);
  if (px_storage.size() < n)
    px_storage.resize(n);
  if (py_storage.size() < n)
    py_storage.resize(n);
  if (pz_storage.size() < n)
    pz_storage.resize(n);

  const double inf = numeric_limits<double>::infinity();
  fill(best_storage.begin(), best_storage.begin() + n, inf);
  fill(used_storage.begin(), used_storage.begin() + n, static_cast<unsigned char>(0));
  for (size_t i = 0; i < n; ++i) {
    const Vec3 &point = placed_points[anchor_indices[i]];
    px_storage[i] = point.x;
    py_storage[i] = point.y;
    pz_storage[i] = point.z;
  }

  double *best = best_storage.data();
  unsigned char *used = used_storage.data();
  const double *px = px_storage.data();
  const double *py = py_storage.data();
  const double *pz = pz_storage.data();
  best[0] = 0.0;
  const bool has_incumbent =
      (incumbent != nullptr && incumbent->valid &&
       isfinite(incumbent->anchor_mst2) && isfinite(incumbent->anchor_max_d2));
  const double incumbent_anchor_max_d2 =
      has_incumbent ? incumbent->anchor_max_d2 : inf;

  for (size_t step = 0; step < n; ++step) {
    double next_cost = inf;
    const size_t next = select_anchor_next(best, used, n, next_cost);
    if (next == numeric_limits<size_t>::max())
      return false;
    used[next] = 1;
    mst2 += next_cost;
    max_d2 = max(max_d2, next_cost);
    if (has_incumbent && max_d2 > incumbent_anchor_max_d2 + 1.0e-6) {
      return false;
    }
    const double sx = px[next];
    const double sy = py[next];
    const double sz = pz[next];
    relax_anchor_best(best, px, py, pz, n, sx, sy, sz);
  }
  return true;
}

double anchor_axis_value(const Vec3 &point, int axis) {
  switch (axis) {
  case 1:
    return point.y;
  case 2:
    return point.z;
  case 0:
  default:
    return point.x;
  }
}

double anchor_kd_bbox_distance2(const PbctopoAnchorKdNode &node,
                                const Vec3 &point) {
  const double dx = (point.x < node.min_x)
                        ? (node.min_x - point.x)
                        : ((point.x > node.max_x) ? (point.x - node.max_x) : 0.0);
  const double dy = (point.y < node.min_y)
                        ? (node.min_y - point.y)
                        : ((point.y > node.max_y) ? (point.y - node.max_y) : 0.0);
  const double dz = (point.z < node.min_z)
                        ? (node.min_z - point.z)
                        : ((point.z > node.max_z) ? (point.z - node.max_z) : 0.0);
  return dx * dx + dy * dy + dz * dz;
}

int build_anchor_kdtree_recursive(vector<size_t> &order, size_t begin, size_t end,
                                  const vector<Vec3> &points,
                                  vector<PbctopoAnchorKdNode> &nodes) {
  if (begin >= end)
    return -1;

  double min_x = numeric_limits<double>::infinity();
  double max_x = -numeric_limits<double>::infinity();
  double min_y = numeric_limits<double>::infinity();
  double max_y = -numeric_limits<double>::infinity();
  double min_z = numeric_limits<double>::infinity();
  double max_z = -numeric_limits<double>::infinity();
  for (size_t idx = begin; idx < end; ++idx) {
    const Vec3 &point = points[order[idx]];
    min_x = min(min_x, point.x);
    max_x = max(max_x, point.x);
    min_y = min(min_y, point.y);
    max_y = max(max_y, point.y);
    min_z = min(min_z, point.z);
    max_z = max(max_z, point.z);
  }

  int axis = 0;
  const double span_x = max_x - min_x;
  const double span_y = max_y - min_y;
  const double span_z = max_z - min_z;
  if (span_y > span_x && span_y >= span_z) {
    axis = 1;
  } else if (span_z > span_x && span_z > span_y) {
    axis = 2;
  }

  const size_t mid = begin + (end - begin) / 2;
  nth_element(order.begin() + begin, order.begin() + mid, order.begin() + end,
              [&](size_t lhs, size_t rhs) {
                return anchor_axis_value(points[lhs], axis) <
                       anchor_axis_value(points[rhs], axis);
              });

  const int node_idx = static_cast<int>(nodes.size());
  nodes.push_back({});
  const int left = build_anchor_kdtree_recursive(order, begin, mid, points, nodes);
  const int right =
      build_anchor_kdtree_recursive(order, mid + 1, end, points, nodes);

  PbctopoAnchorKdNode &node = nodes[node_idx];
  node.point_idx = order[mid];
  node.axis = axis;
  node.left = left;
  node.right = right;
  node.min_x = min_x;
  node.max_x = max_x;
  node.min_y = min_y;
  node.max_y = max_y;
  node.min_z = min_z;
  node.max_z = max_z;
  return node_idx;
}

void query_anchor_kdtree_knn(
    const vector<PbctopoAnchorKdNode> &nodes, int node_idx,
    const vector<Vec3> &points, const Vec3 &query_point, size_t query_idx,
    size_t knn_k, priority_queue<pair<double, size_t>> &neighbors) {
  if (node_idx < 0)
    return;

  const PbctopoAnchorKdNode &node = nodes[static_cast<size_t>(node_idx)];
  const Vec3 &point = points[node.point_idx];
  if (node.point_idx != query_idx) {
    const double dx = query_point.x - point.x;
    const double dy = query_point.y - point.y;
    const double dz = query_point.z - point.z;
    const double d2 = dx * dx + dy * dy + dz * dz;
    if (neighbors.size() < knn_k) {
      neighbors.push({d2, node.point_idx});
    } else if (d2 < neighbors.top().first) {
      neighbors.pop();
      neighbors.push({d2, node.point_idx});
    }
  }

  const double query_axis = anchor_axis_value(query_point, node.axis);
  const double split_axis = anchor_axis_value(point, node.axis);
  const int near_child = (query_axis <= split_axis) ? node.left : node.right;
  const int far_child = (query_axis <= split_axis) ? node.right : node.left;

  query_anchor_kdtree_knn(nodes, near_child, points, query_point, query_idx,
                          knn_k, neighbors);
  if (far_child < 0)
    return;

  const double worst_d2 =
      (neighbors.size() < knn_k) ? numeric_limits<double>::infinity()
                                 : neighbors.top().first;
  if (anchor_kd_bbox_distance2(nodes[static_cast<size_t>(far_child)], query_point) <=
      worst_d2 + kDistanceEps) {
    query_anchor_kdtree_knn(nodes, far_child, points, query_point, query_idx,
                            knn_k, neighbors);
  }
}

bool build_anchor_knn_edges_kdtree(
    const vector<Vec3> &placed_points, const vector<size_t> &anchor_indices,
    const vector<PbctopoContinuityPair> &continuity_pairs, size_t knn_k,
    vector<PbctopoSparseAnchorEdge> &edges, double &build_s, double &query_s) {
  edges.clear();
  build_s = 0.0;
  query_s = 0.0;
  const size_t n = anchor_indices.size();
  if (n < 2)
    return true;

  vector<Vec3> anchor_points(n);
  unordered_map<size_t, size_t> anchor_to_local;
  anchor_to_local.reserve(n * 2 + 1);
  for (size_t idx = 0; idx < n; ++idx) {
    anchor_points[idx] = placed_points[anchor_indices[idx]];
    anchor_to_local.emplace(anchor_indices[idx], idx);
  }

  const auto build_started = chrono::steady_clock::now();
  vector<size_t> order(n);
  iota(order.begin(), order.end(), 0);
  vector<PbctopoAnchorKdNode> nodes;
  nodes.reserve(n);
  const int root = build_anchor_kdtree_recursive(order, 0, n, anchor_points, nodes);
  build_s = chrono::duration<double>(chrono::steady_clock::now() - build_started)
                .count();
  if (root < 0)
    return false;

  unordered_set<uint64_t> seen_edges;
  seen_edges.reserve(n * (knn_k + 2));
  edges.reserve(n * (knn_k + 2));
  auto append_unique_edge = [&](size_t a, size_t b, double d2) {
    if (a == b)
      return;
    const size_t lo = min(a, b);
    const size_t hi = max(a, b);
    const uint64_t key =
        (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi);
    if (!seen_edges.insert(key).second)
      return;
    edges.push_back({lo, hi, d2});
  };

  for (const auto &pair : continuity_pairs) {
    const auto it_a = anchor_to_local.find(pair.a);
    if (it_a == anchor_to_local.end())
      continue;
    const auto it_b = anchor_to_local.find(pair.b);
    if (it_b == anchor_to_local.end())
      continue;
    const Vec3 &lhs_point = placed_points[pair.a];
    const Vec3 &rhs_point = placed_points[pair.b];
    const double dx = lhs_point.x - rhs_point.x;
    const double dy = lhs_point.y - rhs_point.y;
    const double dz = lhs_point.z - rhs_point.z;
    append_unique_edge(it_a->second, it_b->second, dx * dx + dy * dy + dz * dz);
  }

  const size_t k_eff = min(knn_k, n - 1);
  const auto query_started = chrono::steady_clock::now();
  for (size_t query_idx = 0; query_idx < n; ++query_idx) {
    priority_queue<pair<double, size_t>> neighbors;
    query_anchor_kdtree_knn(nodes, root, anchor_points, anchor_points[query_idx],
                            query_idx, k_eff, neighbors);
    while (!neighbors.empty()) {
      append_unique_edge(query_idx, neighbors.top().second,
                         neighbors.top().first);
      neighbors.pop();
    }
  }
  query_s = chrono::duration<double>(chrono::steady_clock::now() - query_started)
                .count();
  return true;
}

bool compute_sparse_anchor_mst(const vector<PbctopoSparseAnchorEdge> &edges,
                               size_t node_count, double &mst2,
                               double &max_d2) {
  mst2 = 0.0;
  max_d2 = 0.0;
  if (node_count < 2)
    return true;
  if (edges.size() + 1 < node_count)
    return false;

  vector<PbctopoSparseAnchorEdge> ranked_edges = edges;
  sort(ranked_edges.begin(), ranked_edges.end(),
       [](const PbctopoSparseAnchorEdge &lhs,
          const PbctopoSparseAnchorEdge &rhs) { return lhs.d2 < rhs.d2; });

  vector<size_t> parent(node_count);
  vector<unsigned char> rank(node_count, 0);
  for (size_t idx = 0; idx < node_count; ++idx)
    parent[idx] = idx;

  const auto find_root = [&](size_t node, auto &&self) -> size_t {
    if (parent[node] == node)
      return node;
    parent[node] = self(parent[node], self);
    return parent[node];
  };
  auto unite = [&](size_t lhs, size_t rhs) {
    lhs = find_root(lhs, find_root);
    rhs = find_root(rhs, find_root);
    if (lhs == rhs)
      return false;
    if (rank[lhs] < rank[rhs])
      swap(lhs, rhs);
    parent[rhs] = lhs;
    if (rank[lhs] == rank[rhs])
      ++rank[lhs];
    return true;
  };

  size_t used_edges = 0;
  for (const auto &edge : ranked_edges) {
    if (!unite(edge.a, edge.b))
      continue;
    mst2 += edge.d2;
    max_d2 = max(max_d2, edge.d2);
    ++used_edges;
    if (used_edges + 1 == node_count)
      return true;
  }
  return false;
}

bool anchor_graph_metrics_spatial_or_dense(
    const vector<Vec3> &placed_points, const vector<size_t> &anchor_indices,
    const vector<PbctopoContinuityPair> &continuity_pairs, double &mst2,
    double &max_d2, const PbctopoCutCandidate *incumbent,
    PbctopoRefineScorePhase phase, PbctopoAnchorBackend backend,
    PbctopoTimingStats *timing) {
  const bool allow_spatial =
      backend == PbctopoAnchorBackend::SpatialAll &&
      anchor_indices.size() >= kPbctopoSpatialAnchorThreshold;
  if (!allow_spatial) {
    return anchor_graph_metrics_ok(placed_points, anchor_indices, mst2, max_d2,
                                   incumbent);
  }

  if (timing != nullptr)
    ++timing->component_graph_refine_anchor_spatial_calls;
  vector<PbctopoSparseAnchorEdge> sparse_edges;
  double build_s = 0.0;
  double query_s = 0.0;
  const bool built = build_anchor_knn_edges_kdtree(
      placed_points, anchor_indices, continuity_pairs, kPbctopoSpatialAnchorK,
      sparse_edges, build_s, query_s);
  if (timing != nullptr) {
    timing->component_graph_refine_anchor_spatial_build_s += build_s;
    timing->component_graph_refine_anchor_spatial_query_s += query_s;
  }
  const bool sparse_ok =
      built && compute_sparse_anchor_mst(sparse_edges, anchor_indices.size(), mst2,
                                        max_d2);
  if (sparse_ok)
    return true;

  if (timing != nullptr) {
    ++timing->component_graph_refine_anchor_spatial_fallbacks;
    ++timing->component_graph_refine_anchor_spatial_disconnects;
  }
  const auto dense_started = chrono::steady_clock::now();
  const bool dense_ok = anchor_graph_metrics_ok(placed_points, anchor_indices, mst2,
                                                max_d2, incumbent);
  if (timing != nullptr) {
    timing->component_graph_refine_anchor_dense_fallback_s +=
        chrono::duration<double>(chrono::steady_clock::now() - dense_started)
            .count();
  }
  return dense_ok;
}

PbctopoCutCandidate score_component_placed_points_solution(
    const vector<Vec3> &placed_points, double shift_norm2,
    const vector<PbctopoLocalEdge> &local_edges,
    const unordered_set<uint64_t> &edge_set, double clash_cutoff,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices, double continuity_cutoff,
    const PbctopoCutCandidate *incumbent = nullptr,
    PbctopoTimingStats *timing = nullptr, bool run_ultraclash = true,
    PbctopoRefineScorePhase phase = PbctopoRefineScorePhase::Exact,
    const PbctopoStage1EdgeContext *stage1_edge_context = nullptr,
    const PbctopoEdgeScanSoA *edge_scan_soa = nullptr,
    const PbctopoCutCandidate *edge_prune_threshold = nullptr,
    PbctopoAnchorBackend anchor_backend = PbctopoAnchorBackend::Auto,
    const PointStats *precomputed_stats = nullptr) {
  if (timing != nullptr) {
    switch (phase) {
    case PbctopoRefineScorePhase::Exact:
      ++timing->component_graph_refine_score_calls;
      break;
    case PbctopoRefineScorePhase::CoarseStage1:
      ++timing->component_graph_refine_coarse_score_calls;
      ++timing->component_graph_refine_coarse_stage1_score_calls;
      break;
    case PbctopoRefineScorePhase::CoarseStage2:
      ++timing->component_graph_refine_coarse_score_calls;
      ++timing->component_graph_refine_coarse_stage2_score_calls;
      break;
    }
  }
  auto accumulate_continuity = [&](double delta) {
    if (timing == nullptr)
      return;
    switch (phase) {
    case PbctopoRefineScorePhase::Exact:
      timing->component_graph_refine_continuity_s += delta;
      break;
    case PbctopoRefineScorePhase::CoarseStage1:
      timing->component_graph_refine_coarse_s += delta;
      timing->component_graph_refine_coarse_stage1_s += delta;
      timing->component_graph_refine_coarse_stage1_continuity_s += delta;
      break;
    case PbctopoRefineScorePhase::CoarseStage2:
      timing->component_graph_refine_coarse_s += delta;
      timing->component_graph_refine_coarse_stage2_s += delta;
      timing->component_graph_refine_coarse_stage2_continuity_s += delta;
      break;
    }
  };
  auto accumulate_anchor = [&](double delta) {
    if (timing == nullptr)
      return;
    switch (phase) {
    case PbctopoRefineScorePhase::Exact:
      timing->component_graph_refine_anchor_s += delta;
      break;
    case PbctopoRefineScorePhase::CoarseStage1:
      timing->component_graph_refine_coarse_s += delta;
      timing->component_graph_refine_coarse_stage1_s += delta;
      timing->component_graph_refine_coarse_stage1_anchor_s += delta;
      break;
    case PbctopoRefineScorePhase::CoarseStage2:
      timing->component_graph_refine_coarse_s += delta;
      timing->component_graph_refine_coarse_stage2_s += delta;
      timing->component_graph_refine_coarse_stage2_anchor_s += delta;
      break;
    }
  };
  auto accumulate_cut_score = [&](double delta) {
    if (timing == nullptr)
      return;
    switch (phase) {
    case PbctopoRefineScorePhase::Exact:
      timing->component_graph_refine_cut_score_s += delta;
      break;
    case PbctopoRefineScorePhase::CoarseStage1:
      timing->component_graph_refine_coarse_s += delta;
      timing->component_graph_refine_coarse_stage1_s += delta;
      timing->component_graph_refine_coarse_stage1_cut_score_s += delta;
      break;
    case PbctopoRefineScorePhase::CoarseStage2:
      timing->component_graph_refine_coarse_s += delta;
      timing->component_graph_refine_coarse_stage2_s += delta;
        timing->component_graph_refine_coarse_stage2_cut_score_s += delta;
      break;
    }
  };
  chrono::steady_clock::time_point stage_started;
  if (timing != nullptr)
    stage_started = chrono::steady_clock::now();
  double continuity_path = numeric_limits<double>::infinity();
  double continuity_max_d2 = numeric_limits<double>::infinity();
  if (!continuity_metrics_ok(placed_points, continuity_pairs, continuity_cutoff,
                             continuity_path, continuity_max_d2, incumbent)) {
    if (timing != nullptr) {
      accumulate_continuity(
          chrono::duration<double>(chrono::steady_clock::now() - stage_started)
              .count());
    }
    return {};
  }
  if (timing != nullptr) {
    accumulate_continuity(
        chrono::duration<double>(chrono::steady_clock::now() - stage_started)
            .count());
    stage_started = chrono::steady_clock::now();
  }
  double anchor_mst2 = numeric_limits<double>::infinity();
  double anchor_max_d2 = numeric_limits<double>::infinity();
  if (!anchor_graph_metrics_spatial_or_dense(
          placed_points, anchor_indices, continuity_pairs, anchor_mst2,
          anchor_max_d2, incumbent, phase, anchor_backend, timing)) {
    if (timing != nullptr) {
      accumulate_anchor(
          chrono::duration<double>(chrono::steady_clock::now() - stage_started)
              .count());
    }
    return {};
  }
  if (timing != nullptr) {
    accumulate_anchor(
        chrono::duration<double>(chrono::steady_clock::now() - stage_started)
            .count());
    stage_started = chrono::steady_clock::now();
  }
  PointStats stats_storage;
  const PointStats *stats_ptr = precomputed_stats;
  if (stats_ptr == nullptr) {
    const auto stats_started = chrono::steady_clock::now();
    stats_storage = build_point_stats(placed_points);
    stats_ptr = &stats_storage;
    if (timing != nullptr) {
      ++timing->component_graph_refine_stats_full;
      timing->component_graph_refine_stats_s +=
          chrono::duration<double>(chrono::steady_clock::now() -
                                   stats_started)
              .count();
    }
  }
  const PointStats &stats = *stats_ptr;
  PbctopoCutCandidate candidate;
  PbctopoCutScoreBreakdown cut_breakdown;
  PbctopoCutCandidate pre_edge_seed;
  pre_edge_seed.continuity_path = continuity_path;
  pre_edge_seed.continuity_max_d2 = continuity_max_d2;
  pre_edge_seed.anchor_mst2 = anchor_mst2;
  pre_edge_seed.anchor_max_d2 = anchor_max_d2;
  if (phase == PbctopoRefineScorePhase::CoarseStage1) {
    candidate = score_cut_candidate_principal_only(
        local_edges, placed_points, 1.5, 1.5, 1.5, 0.0, stats, shift_norm2,
        &cut_breakdown, stage1_edge_context, edge_scan_soa,
        edge_prune_threshold, &pre_edge_seed);
  } else {
    candidate =
        score_cut_candidate(local_edges, placed_points, 1.5, 1.5, 1.5, 0.0,
                            stats, shift_norm2,
                            (phase == PbctopoRefineScorePhase::CoarseStage2)
                                ? stage1_edge_context
                                : nullptr,
                            edge_scan_soa, edge_prune_threshold, &pre_edge_seed,
                            phase);
  }
  if (!candidate.valid) {
    if (timing != nullptr) {
      accumulate_cut_score(
          chrono::duration<double>(chrono::steady_clock::now() - stage_started)
              .count());
    }
    return {};
  }
  if (timing != nullptr && phase == PbctopoRefineScorePhase::CoarseStage1) {
    timing->component_graph_refine_coarse_stage1_cut_score_principal_s +=
        cut_breakdown.principal_s;
    timing->component_graph_refine_coarse_stage1_cut_score_edge_scan_s +=
        cut_breakdown.edge_scan_s;
  }
  candidate.continuity_path = continuity_path;
  candidate.continuity_max_d2 = continuity_max_d2;
  candidate.anchor_mst2 = anchor_mst2;
  candidate.anchor_max_d2 = anchor_max_d2;
  if (cut_candidate_cannot_beat_incumbent(candidate, incumbent)) {
    if (timing != nullptr) {
      accumulate_cut_score(
          chrono::duration<double>(chrono::steady_clock::now() - stage_started)
              .count());
      if (phase == PbctopoRefineScorePhase::Exact)
        ++timing->component_graph_refine_preclash_prunes;
    }
    return {};
  }
  if (timing != nullptr) {
    accumulate_cut_score(
        chrono::duration<double>(chrono::steady_clock::now() - stage_started)
            .count());
  }
  if (run_ultraclash) {
    if (timing != nullptr)
      stage_started = chrono::steady_clock::now();
    analyze_cut_candidate_ultraclash(placed_points, edge_set, candidate,
                                     clash_cutoff);
    if (timing != nullptr) {
      timing->component_graph_refine_ultraclash_s +=
          chrono::duration<double>(chrono::steady_clock::now() - stage_started)
              .count();
    }
  }
  return candidate;
}

PbctopoCutCandidate score_component_offset_solution(
    const vector<PbctopoAtom> &atoms, const vector<vector<size_t>> &components,
    box Box, const vector<Int3> &offsets,
    const vector<PbctopoLocalEdge> &local_edges,
    const unordered_set<uint64_t> &edge_set, double clash_cutoff,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices,
    double continuity_cutoff, vector<Vec3> &placed_points,
    const PbctopoCutCandidate *incumbent = nullptr,
    PbctopoTimingStats *timing = nullptr, bool run_ultraclash = true,
    PbctopoRefineScorePhase phase = PbctopoRefineScorePhase::Exact,
    const PbctopoStage1EdgeContext *stage1_edge_context = nullptr,
    const PbctopoEdgeScanSoA *edge_scan_soa = nullptr,
    const PbctopoCutCandidate *edge_prune_threshold = nullptr,
    PbctopoAnchorBackend anchor_backend = PbctopoAnchorBackend::Auto) {
  auto accumulate_build_points = [&](double delta) {
    if (timing == nullptr)
      return;
    switch (phase) {
    case PbctopoRefineScorePhase::Exact:
      timing->component_graph_refine_build_points_s += delta;
      break;
    case PbctopoRefineScorePhase::CoarseStage1:
      timing->component_graph_refine_coarse_s += delta;
      timing->component_graph_refine_coarse_stage1_s += delta;
      timing->component_graph_refine_coarse_stage1_build_points_s += delta;
      break;
    case PbctopoRefineScorePhase::CoarseStage2:
      timing->component_graph_refine_coarse_s += delta;
      timing->component_graph_refine_coarse_stage2_s += delta;
      timing->component_graph_refine_coarse_stage2_build_points_s += delta;
      break;
    }
  };
  double shift_norm2 = 0.0;
  chrono::steady_clock::time_point stage_started;
  if (timing != nullptr)
    stage_started = chrono::steady_clock::now();
  build_component_points_from_offsets(atoms, components, Box, offsets, placed_points,
                                      shift_norm2);
  if (timing != nullptr) {
    accumulate_build_points(
        chrono::duration<double>(chrono::steady_clock::now() - stage_started)
            .count());
  }
  return score_component_placed_points_solution(
      placed_points, shift_norm2, local_edges, edge_set, clash_cutoff,
      continuity_pairs, anchor_indices, continuity_cutoff, incumbent, timing,
      run_ultraclash, phase, stage1_edge_context, edge_scan_soa,
      edge_prune_threshold, anchor_backend);
}

PbctopoCutCandidate improve_component_offsets_locally(
    const vector<PbctopoAtom> &atoms, const vector<vector<size_t>> &components,
    box Box, const vector<PbctopoShift> &shifts,
    const vector<PbctopoLocalEdge> &local_edges,
    const unordered_set<uint64_t> &edge_set,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices, vector<Int3> &offsets,
    double clash_cutoff, double continuity_cutoff, vector<Vec3> &placed_points,
    PbctopoAnchorBackend anchor_backend = PbctopoAnchorBackend::Auto,
    PbctopoTimingStats *timing = nullptr, PbctopoWarmState *warm_state = nullptr,
    const PbctopoCutCandidate *warm_reference_candidate = nullptr) {
  const vector<size_t> &exact_anchor_indices = anchor_indices;
  PbctopoCutCandidate best = score_component_offset_solution(
      atoms, components, Box, offsets, local_edges, edge_set, clash_cutoff,
      continuity_pairs, exact_anchor_indices, continuity_cutoff, placed_points,
      nullptr, timing, true, PbctopoRefineScorePhase::Exact, nullptr, nullptr,
      nullptr, anchor_backend);
  if (!best.valid || offsets.size() < 2) {
    if (warm_state != nullptr && (!best.valid || best.ultra_clash))
      clear_component_trial_hints(*warm_state);
    return best;
  }
  PbctopoRefinePlacementWorkspace refine_workspace;
  refine_workspace.initialize(atoms, components, Box);

  const vector<size_t> half_anchor_indices =
      build_coarse_anchor_subset(atoms, anchor_indices);
  const vector<size_t> third_anchor_indices =
      build_third_anchor_subset(atoms, anchor_indices);
  const vector<size_t> tapered_stage1_anchor_indices =
      build_tapered_stage1_anchor_subset(atoms, anchor_indices);
  const bool use_three_stage_anchor =
      anchor_indices.size() > 256 &&
      tapered_stage1_anchor_indices.size() >= 2 &&
      third_anchor_indices.size() >= 2 &&
      tapered_stage1_anchor_indices.size() + 4 < third_anchor_indices.size() &&
      third_anchor_indices.size() + 2 < anchor_indices.size();
  const bool use_two_stage_anchor =
      !use_three_stage_anchor &&
      anchor_indices.size() >= 8 && half_anchor_indices.size() >= 2 &&
      half_anchor_indices.size() + 2 < anchor_indices.size();
  const vector<size_t> atom_component_index =
      (use_two_stage_anchor || use_three_stage_anchor)
          ? build_atom_component_index(atoms.size(), components)
          : vector<size_t>{};
  const PbctopoEdgeScanSoA edge_scan_soa =
      (use_two_stage_anchor || use_three_stage_anchor)
          ? build_edge_scan_soa(local_edges)
          : PbctopoEdgeScanSoA{};
  const bool warm_offsets_compatible =
      (warm_state != nullptr && warm_state->has_component_offsets &&
       warm_state->previous_component_offsets.size() == offsets.size() &&
       same_component_layout(components, warm_state->previous_component_sizes));
  const bool warm_trial_hints_compatible =
      (warm_state != nullptr &&
       component_trial_hints_compatible(*warm_state, offsets.size()) &&
       same_component_layout(components, warm_state->previous_component_sizes));
  const bool strong_warm_candidate =
      (warm_reference_candidate != nullptr && warm_reference_candidate->valid &&
       !warm_reference_candidate->ultra_clash &&
       candidate_close_enough_to_skip_search(best, *warm_reference_candidate));
  size_t stage1_keep = 6;
  size_t stage2_keep = 3;
  size_t exact_keep = 3;
  if (strong_warm_candidate &&
      (use_two_stage_anchor || use_three_stage_anchor)) {
    if (use_three_stage_anchor)
      stage1_keep = 4;
    stage2_keep = 2;
    exact_keep = 2;
    if (timing != nullptr)
      ++timing->component_refine_fanout_reduced;
  }

  struct PbctopoRefineTrialCandidate {
    vector<Int3> offsets;
    PbctopoCutCandidate candidate;
  };
  auto push_top_refine_trial =
      [](vector<PbctopoRefineTrialCandidate> &top_trials,
         const vector<Int3> &trial_offsets,
         const PbctopoCutCandidate &trial_candidate, size_t max_keep,
         bool (*better_candidate)(const PbctopoCutCandidate &,
                                  const PbctopoCutCandidate &)) {
        if (!trial_candidate.valid)
          return;
        PbctopoRefineTrialCandidate entry;
        entry.offsets = trial_offsets;
        entry.candidate = trial_candidate;
        auto insert_pos = top_trials.begin();
        while (insert_pos != top_trials.end() &&
               !better_candidate(entry.candidate, insert_pos->candidate)) {
          ++insert_pos;
        }
         top_trials.insert(insert_pos, std::move(entry));
         if (top_trials.size() > max_keep)
           top_trials.resize(max_keep);
      };
  vector<vector<Int3>> stage1_absolute_hint_offsets(offsets.size());
  vector<vector<Int3>> stage2_absolute_hint_offsets(offsets.size());
  auto score_refine_trial =
      [&](size_t trial_comp_idx, const vector<Int3> &trial_offsets_ref,
          const vector<size_t> &trial_anchor_indices, vector<Vec3> &trial_points,
          const PbctopoCutCandidate *incumbent, bool run_ultraclash,
          PbctopoRefineScorePhase phase,
          const PbctopoStage1EdgeContext *stage1_edge_context,
          const PbctopoEdgeScanSoA *edge_scan_soa,
          const PbctopoCutCandidate *edge_prune_threshold) {
        PointStats trial_stats;
        double trial_shift_norm2 = 0.0;
        refine_workspace.prepare_trial(trial_comp_idx, trial_offsets_ref,
                                       trial_points, trial_stats,
                                       trial_shift_norm2, timing, phase);
        return score_component_placed_points_solution(
            trial_points, trial_shift_norm2, local_edges, edge_set,
            clash_cutoff, continuity_pairs, trial_anchor_indices,
            continuity_cutoff, incumbent, timing, run_ultraclash, phase,
            stage1_edge_context, edge_scan_soa, edge_prune_threshold,
            anchor_backend, &trial_stats);
      };
  auto capture_component_hint_offsets =
      [](const vector<PbctopoRefineTrialCandidate> &source_trials,
         size_t comp_idx, vector<vector<Int3>> &target_offsets) {
        if (comp_idx >= target_offsets.size())
          return;
        auto &dst = target_offsets[comp_idx];
        dst.clear();
        dst.reserve(source_trials.size());
        for (const auto &trial : source_trials) {
          if (comp_idx < trial.offsets.size())
            dst.push_back(trial.offsets[comp_idx]);
        }
      };
  for (int pass = 0; pass < 2; ++pass) {
    bool improved = false;
    vector<Int3> trial_offsets = offsets;
    vector<Vec3> local_points;
    local_points.reserve(atoms.size());
    vector<Vec3> coarse_points;
    coarse_points.reserve(atoms.size());
    for (size_t comp_idx = 1; comp_idx < offsets.size(); ++comp_idx) {
      refine_workspace.reset_baseline(offsets, placed_points);
      vector<Int3> best_offsets = offsets;
      PbctopoCutCandidate best_local = best;
      PbctopoStage1EdgeContext stage1_edge_context;
      if (use_two_stage_anchor || use_three_stage_anchor) {
        stage1_edge_context = build_stage1_edge_context(
            local_edges, atom_component_index, comp_idx, placed_points);
      }
      Int3 warm_target_delta;
      const Int3 *warm_target_delta_ptr = nullptr;
      if (warm_offsets_compatible &&
          comp_idx < warm_state->previous_component_offsets.size()) {
        warm_target_delta =
            warm_state->previous_component_offsets[comp_idx] - offsets[comp_idx];
        warm_target_delta_ptr = &warm_target_delta;
      }
      const vector<Int3> *stage1_hint_deltas =
          (warm_trial_hints_compatible &&
           comp_idx < warm_state->previous_stage1_delta_hints_by_component.size())
              ? &warm_state->previous_stage1_delta_hints_by_component[comp_idx]
              : nullptr;
      const vector<Int3> *stage2_hint_deltas =
          (warm_trial_hints_compatible &&
           comp_idx < warm_state->previous_stage2_delta_hints_by_component.size())
              ? &warm_state->previous_stage2_delta_hints_by_component[comp_idx]
              : nullptr;
      const vector<size_t> ordered_shift_indices = build_refine_shift_order(
          shifts, warm_target_delta_ptr, stage2_hint_deltas, stage1_hint_deltas,
          timing);
      if (!use_two_stage_anchor && !use_three_stage_anchor) {
        for (size_t shift_pos : ordered_shift_indices) {
          const auto &shift = shifts[shift_pos];
          trial_offsets = offsets;
          trial_offsets[comp_idx] =
              offsets[comp_idx] + Int3{shift.ix, shift.iy, shift.iz};
          if (timing != nullptr) {
            const auto rebase_started = chrono::steady_clock::now();
            rebase_component_offsets(trial_offsets, components);
            timing->component_graph_refine_rebase_s +=
                chrono::duration<double>(chrono::steady_clock::now() -
                                         rebase_started)
                    .count();
          } else {
            rebase_component_offsets(trial_offsets, components);
          }
          PbctopoCutCandidate trial = score_refine_trial(
              comp_idx, trial_offsets, exact_anchor_indices, local_points,
              &best_local, true, PbctopoRefineScorePhase::Exact, nullptr,
              nullptr, nullptr);
          if (better_cut_candidate(trial, best_local)) {
            best_local = trial;
            best_offsets = trial_offsets;
            placed_points = local_points;
          }
        }
      } else if (use_two_stage_anchor) {
        vector<PbctopoRefineTrialCandidate> top_trials;
        top_trials.reserve(stage2_keep);
        for (size_t shift_pos : ordered_shift_indices) {
          const auto &shift = shifts[shift_pos];
          trial_offsets = offsets;
          trial_offsets[comp_idx] =
              offsets[comp_idx] + Int3{shift.ix, shift.iy, shift.iz};
          if (timing != nullptr) {
            const auto rebase_started = chrono::steady_clock::now();
            rebase_component_offsets(trial_offsets, components);
            timing->component_graph_refine_rebase_s +=
                chrono::duration<double>(chrono::steady_clock::now() -
                                         rebase_started)
                    .count();
          } else {
            rebase_component_offsets(trial_offsets, components);
          }
          const PbctopoCutCandidate *stage2_prune_threshold =
              (top_trials.size() >= stage2_keep) ? &top_trials.back().candidate
                                                 : nullptr;
          PbctopoCutCandidate coarse_trial = score_refine_trial(
              comp_idx, trial_offsets, half_anchor_indices, coarse_points,
              nullptr, false, PbctopoRefineScorePhase::CoarseStage2,
              &stage1_edge_context, &edge_scan_soa, stage2_prune_threshold);
          push_top_refine_trial(top_trials, trial_offsets, coarse_trial,
                                stage2_keep,
                                better_preclash_cut_candidate);
        }
        capture_component_hint_offsets(top_trials, comp_idx,
                                       stage2_absolute_hint_offsets);
        for (size_t shortlist_idx = 0;
             shortlist_idx < top_trials.size() && shortlist_idx < exact_keep;
             ++shortlist_idx) {
          const auto &shortlisted = top_trials[shortlist_idx];
          if (timing != nullptr)
            ++timing->component_graph_refine_exact_rescore_calls;
          PbctopoCutCandidate trial = score_refine_trial(
              comp_idx, shortlisted.offsets, exact_anchor_indices, local_points,
              &best_local, true, PbctopoRefineScorePhase::Exact, nullptr,
              nullptr, nullptr);
          if (better_cut_candidate(trial, best_local)) {
            best_local = trial;
            best_offsets = shortlisted.offsets;
            placed_points = local_points;
          }
        }
      } else {
        vector<PbctopoRefineTrialCandidate> sparse_top_trials;
        sparse_top_trials.reserve(stage1_keep);
        for (size_t shift_pos : ordered_shift_indices) {
          const auto &shift = shifts[shift_pos];
          trial_offsets = offsets;
          trial_offsets[comp_idx] =
              offsets[comp_idx] + Int3{shift.ix, shift.iy, shift.iz};
          if (timing != nullptr) {
            const auto rebase_started = chrono::steady_clock::now();
            rebase_component_offsets(trial_offsets, components);
            timing->component_graph_refine_rebase_s +=
                chrono::duration<double>(chrono::steady_clock::now() -
                                         rebase_started)
                    .count();
          } else {
            rebase_component_offsets(trial_offsets, components);
          }
          const PbctopoCutCandidate *stage1_prune_threshold =
              (sparse_top_trials.size() >= stage1_keep)
                  ? &sparse_top_trials.back().candidate
                  : nullptr;
          PbctopoCutCandidate sparse_trial = score_refine_trial(
              comp_idx, trial_offsets, tapered_stage1_anchor_indices,
              coarse_points, nullptr, false,
              PbctopoRefineScorePhase::CoarseStage1, &stage1_edge_context,
              &edge_scan_soa, stage1_prune_threshold);
          push_top_refine_trial(sparse_top_trials, trial_offsets, sparse_trial,
                                stage1_keep,
                                better_stage1_cut_candidate);
        }
        capture_component_hint_offsets(sparse_top_trials, comp_idx,
                                       stage1_absolute_hint_offsets);

        vector<PbctopoRefineTrialCandidate> top_trials;
        top_trials.reserve(stage2_keep);
        for (const auto &shortlisted : sparse_top_trials) {
          const PbctopoCutCandidate *stage2_prune_threshold =
              (top_trials.size() >= stage2_keep)
                  ? &top_trials.back().candidate
                  : nullptr;
          PbctopoCutCandidate coarse_trial = score_refine_trial(
              comp_idx, shortlisted.offsets, third_anchor_indices,
              coarse_points, nullptr, false,
              PbctopoRefineScorePhase::CoarseStage2, &stage1_edge_context,
              &edge_scan_soa, stage2_prune_threshold);
          push_top_refine_trial(top_trials, shortlisted.offsets, coarse_trial,
                                stage2_keep,
                                better_preclash_cut_candidate);
        }
        capture_component_hint_offsets(top_trials, comp_idx,
                                       stage2_absolute_hint_offsets);

        for (size_t shortlist_idx = 0;
             shortlist_idx < top_trials.size() && shortlist_idx < exact_keep;
             ++shortlist_idx) {
          const auto &shortlisted = top_trials[shortlist_idx];
          if (timing != nullptr)
            ++timing->component_graph_refine_exact_rescore_calls;
          PbctopoCutCandidate trial = score_refine_trial(
              comp_idx, shortlisted.offsets, exact_anchor_indices, local_points,
              &best_local, true, PbctopoRefineScorePhase::Exact, nullptr,
              nullptr, nullptr);
          if (better_cut_candidate(trial, best_local)) {
            best_local = trial;
            best_offsets = shortlisted.offsets;
            placed_points = local_points;
          }
        }
      }
      if (better_cut_candidate(best_local, best)) {
        offsets = best_offsets;
        best = best_local;
        improved = true;
      }
    }
    if (!improved)
      break;
  }

  if (warm_state != nullptr) {
    if (!best.valid || best.ultra_clash) {
      clear_component_trial_hints(*warm_state);
    } else {
      warm_state->previous_stage1_delta_hints_by_component.assign(offsets.size(),
                                                                  {});
      warm_state->previous_stage2_delta_hints_by_component.assign(offsets.size(),
                                                                  {});
      auto convert_absolute_offsets_to_deltas =
          [&](const vector<vector<Int3>> &absolute_offsets,
              vector<vector<Int3>> &delta_hints) {
            for (size_t comp_idx = 0; comp_idx < absolute_offsets.size() &&
                                      comp_idx < delta_hints.size();
                 ++comp_idx) {
              for (const auto &absolute_offset : absolute_offsets[comp_idx]) {
                const Int3 delta = absolute_offset - offsets[comp_idx];
                if (is_zero_delta(delta))
                  continue;
                auto &dst = delta_hints[comp_idx];
                if (find_if(dst.begin(), dst.end(),
                            [&](const Int3 &existing) {
                              return same_delta(existing, delta);
                            }) == dst.end()) {
                  dst.push_back(delta);
                }
              }
            }
          };
      convert_absolute_offsets_to_deltas(
          stage1_absolute_hint_offsets,
          warm_state->previous_stage1_delta_hints_by_component);
      convert_absolute_offsets_to_deltas(
          stage2_absolute_hint_offsets,
          warm_state->previous_stage2_delta_hints_by_component);
      warm_state->has_component_trial_hints = false;
      for (size_t comp_idx = 0;
           comp_idx < warm_state->previous_stage1_delta_hints_by_component.size();
           ++comp_idx) {
        if (!warm_state->previous_stage1_delta_hints_by_component[comp_idx]
                 .empty() ||
            !warm_state->previous_stage2_delta_hints_by_component[comp_idx]
                 .empty()) {
          warm_state->has_component_trial_hints = true;
          break;
        }
      }
    }
  }
  return best;
}

vector<size_t>
component_size_signature(const vector<vector<size_t>> &components) {
  vector<size_t> sizes;
  sizes.reserve(components.size());
  for (const auto &component : components)
    sizes.push_back(component.size());
  return sizes;
}

bool same_component_layout(const vector<vector<size_t>> &components,
                           const vector<size_t> &sizes) {
  if (components.size() != sizes.size())
    return false;
  for (size_t idx = 0; idx < components.size(); ++idx) {
    if (components[idx].size() != sizes[idx])
      return false;
  }
  return true;
}

vector<Int3> extract_component_offsets_from_points(
    const vector<PbctopoAtom> &atoms, const vector<vector<size_t>> &components,
    box Box, const vector<Vec3> &placed_points) {
  vector<Int3> offsets;
  if (placed_points.size() != atoms.size())
    return offsets;
  offsets.resize(components.size());
  const BoxMatrix box_matrix(Box);
  for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
    if (components[comp_idx].empty())
      continue;
    const size_t atom_idx = components[comp_idx].front();
    const double dx = placed_points[atom_idx].x - atoms[atom_idx].x;
    const double dy = placed_points[atom_idx].y - atoms[atom_idx].y;
    const double dz = placed_points[atom_idx].z - atoms[atom_idx].z;
    offsets[comp_idx] = {
        static_cast<int>(llround(box_matrix.inv[0][0] * dx +
                                 box_matrix.inv[0][1] * dy +
                                 box_matrix.inv[0][2] * dz)),
        static_cast<int>(llround(box_matrix.inv[1][0] * dx +
                                 box_matrix.inv[1][1] * dy +
                                 box_matrix.inv[1][2] * dz)),
        static_cast<int>(llround(box_matrix.inv[2][0] * dx +
                                 box_matrix.inv[2][1] * dy +
                                 box_matrix.inv[2][2] * dz))};
  }
  rebase_component_offsets(offsets, components);
  return offsets;
}

const char *pbctopo_component_mode_name(PbctopoComponentMode mode) {
  switch (mode) {
  case PbctopoComponentMode::Disabled:
    return "disabled";
  case PbctopoComponentMode::Exact:
    return "exact";
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
};

PbctopoMergedComponentPlan build_chain_merged_component_plan(
    const vector<PbctopoAtom> &atoms,
    const vector<vector<size_t>> &raw_components) {
  PbctopoMergedComponentPlan plan;
  plan.raw_to_merged.assign(raw_components.size(), numeric_limits<size_t>::max());
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
  current_raw_group.reserve(kPbctopoMergedNodeMaxRawComponents);
  current_atoms.reserve(kPbctopoMergedNodeMaxAtoms);
  int current_owner = 0;
  uint64_t current_chain_key = 0;
  size_t current_last_local_atom = 0;
  bool has_current_group = false;

  auto flush_group = [&]() {
    if (current_raw_group.empty())
      return;
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
    const vector<Int3> &merged_offsets, const vector<vector<size_t>> &raw_components,
    const PbctopoMergedComponentPlan &plan) {
  vector<Int3> raw_offsets(raw_components.size(), {0, 0, 0});
  for (size_t merged_idx = 0; merged_idx < plan.merged_to_raw.size(); ++merged_idx) {
    const Int3 merged_offset =
        merged_idx < merged_offsets.size() ? merged_offsets[merged_idx] : Int3{};
    for (size_t raw_idx : plan.merged_to_raw[merged_idx]) {
      if (raw_idx < raw_offsets.size())
        raw_offsets[raw_idx] = merged_offset;
    }
  }
  rebase_component_offsets(raw_offsets, raw_components);
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
    const Int3 &candidate = raw_offsets[first_raw_idx];
    for (size_t raw_idx : group) {
      if (raw_idx >= raw_offsets.size())
        return false;
      const auto &raw_offset = raw_offsets[raw_idx];
      if (raw_offset.x != candidate.x || raw_offset.y != candidate.y ||
          raw_offset.z != candidate.z) {
        return false;
      }
    }
    merged_offsets.push_back(candidate);
  }
  return true;
}

int compare_principal_rg2(double lhs, double rhs);

bool graph_seed_quick_score_better(const PbctopoGraphSeedQuickScore &lhs,
                                   const PbctopoGraphSeedQuickScore &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  if (fabs(lhs.continuity_max_d2 - rhs.continuity_max_d2) > 1.0e-6)
    return lhs.continuity_max_d2 < rhs.continuity_max_d2;
  if (fabs(lhs.continuity_path - rhs.continuity_path) > 1.0e-6)
    return lhs.continuity_path < rhs.continuity_path;
  if (fabs(lhs.anchor_max_d2 - rhs.anchor_max_d2) > 1.0e-6)
    return lhs.anchor_max_d2 < rhs.anchor_max_d2;
  if (fabs(lhs.anchor_mst2 - rhs.anchor_mst2) > 1.0e-6)
    return lhs.anchor_mst2 < rhs.anchor_mst2;
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
  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-6)
    return lhs.shift_norm2 < rhs.shift_norm2;
  return false;
}

PbctopoGraphSeedQuickScore score_component_graph_seed_quick(
    const vector<PbctopoAtom> &atoms, const vector<vector<size_t>> &components,
    const vector<Vec3> &placed_points, const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices, double continuity_cutoff,
    double shift_norm2) {
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
  if (!anchor_graph_metrics_ok(placed_points, *quick_anchor_indices, score.anchor_mst2,
                               score.anchor_max_d2, nullptr)) {
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

constexpr double kPrincipalRgRelTol = 0.08;
constexpr double kPrincipalRgAbsTol = 5.0;
constexpr double kShapeKappaRelTol = 0.10;
constexpr double kShapeKappaAbsTol = 0.02;
constexpr double kShapeMomentRelTol = 0.08;
constexpr double kShapeMomentAbsTol = 5.0;

double principal_rg2_tolerance(double reference) {
  return max(kPrincipalRgAbsTol, max(0.0, reference) * kPrincipalRgRelTol);
}

bool principal_rg2_significantly_worse(double value, double reference) {
  if (!isfinite(value) || !isfinite(reference))
    return false;
  return value > reference + principal_rg2_tolerance(reference);
}

bool principal_rg2_significantly_better(double value, double reference) {
  if (!isfinite(value) || !isfinite(reference))
    return false;
  return value + principal_rg2_tolerance(reference) < reference;
}

int compare_principal_rg2(double lhs, double rhs) {
  if (principal_rg2_significantly_better(lhs, rhs))
    return -1;
  if (principal_rg2_significantly_better(rhs, lhs))
    return 1;
  return 0;
}

double shape_kappa_tolerance(double reference) {
  return max(kShapeKappaAbsTol, max(0.0, reference) * kShapeKappaRelTol);
}

bool shape_kappa_significantly_worse(double value, double reference) {
  if (!isfinite(value) || !isfinite(reference))
    return false;
  return value > reference + shape_kappa_tolerance(reference);
}

bool shape_kappa_significantly_better(double value, double reference) {
  if (!isfinite(value) || !isfinite(reference))
    return false;
  return value + shape_kappa_tolerance(reference) < reference;
}

int compare_shape_kappa(double lhs, double rhs) {
  if (shape_kappa_significantly_better(lhs, rhs))
    return -1;
  if (shape_kappa_significantly_better(rhs, lhs))
    return 1;
  return 0;
}

double shape_moment_tolerance(double reference) {
  return max(kShapeMomentAbsTol, max(0.0, reference) * kShapeMomentRelTol);
}

bool shape_moment_significantly_worse(double value, double reference) {
  if (!isfinite(value) || !isfinite(reference))
    return false;
  return value > reference + shape_moment_tolerance(reference);
}

bool shape_moment_significantly_better(double value, double reference) {
  if (!isfinite(value) || !isfinite(reference))
    return false;
  return value + shape_moment_tolerance(reference) < reference;
}

int compare_shape_moment(double lhs, double rhs) {
  if (shape_moment_significantly_better(lhs, rhs))
    return -1;
  if (shape_moment_significantly_better(rhs, lhs))
    return 1;
  return 0;
}

bool component_candidate_is_outlier(const PbctopoCutCandidate &candidate,
                                    const PbctopoCutCandidate *reference) {
  if (!candidate.valid || candidate.ultra_clash || reference == nullptr ||
      !reference->valid || reference->ultra_clash) {
    return false;
  }
  return principal_rg2_significantly_worse(candidate.principal_rg2_max,
                                           reference->principal_rg2_max) &&
         (shape_moment_significantly_worse(candidate.shape_asphericity,
                                           reference->shape_asphericity) ||
          shape_kappa_significantly_worse(candidate.shape_kappa2,
                                          reference->shape_kappa2));
}

const char *candidate_skip_search_reject_reason(
    const PbctopoCutCandidate &candidate, const PbctopoCutCandidate &reference) {
  if (!candidate.valid || candidate.ultra_clash)
    return "invalid_or_ultraclash";
  if (!reference.valid || reference.ultra_clash)
    return nullptr;
  if (candidate.broken_edges > reference.broken_edges)
    return "broken_edges";
  if (isfinite(candidate.continuity_max_d2) &&
      isfinite(reference.continuity_max_d2) &&
      candidate.continuity_max_d2 >
          max(reference.continuity_max_d2 * 1.05,
              reference.continuity_max_d2 + 4.0)) {
    return "continuity_max_d2";
  }
  if (isfinite(candidate.continuity_path) &&
      isfinite(reference.continuity_path) &&
      candidate.continuity_path >
          max(reference.continuity_path * 1.05,
              reference.continuity_path + 20.0)) {
    return "continuity_path";
  }
  if (isfinite(candidate.anchor_max_d2) &&
      isfinite(reference.anchor_max_d2) &&
      candidate.anchor_max_d2 >
          max(reference.anchor_max_d2 * 1.10, reference.anchor_max_d2 + 25.0)) {
    return "anchor_max_d2";
  }
  if (isfinite(candidate.anchor_mst2) && isfinite(reference.anchor_mst2) &&
      candidate.anchor_mst2 >
          max(reference.anchor_mst2 * 1.10, reference.anchor_mst2 + 100.0)) {
    return "anchor_mst2";
  }
  if (principal_rg2_significantly_worse(candidate.principal_rg2_max,
                                        reference.principal_rg2_max)) {
    return "principal_rg2_max";
  }
  if (principal_rg2_significantly_worse(candidate.principal_rg2_mid,
                                        reference.principal_rg2_mid)) {
    return "principal_rg2_mid";
  }
  if (principal_rg2_significantly_worse(candidate.principal_rg2_min,
                                        reference.principal_rg2_min)) {
    return "principal_rg2_min";
  }
  if (shape_kappa_significantly_worse(candidate.shape_kappa2,
                                      reference.shape_kappa2)) {
    return "shape_kappa2";
  }
  if (shape_moment_significantly_worse(candidate.shape_asphericity,
                                       reference.shape_asphericity)) {
    return "shape_asphericity";
  }
  if (shape_moment_significantly_worse(candidate.shape_acylindricity,
                                       reference.shape_acylindricity)) {
    return "shape_acylindricity";
  }
  if (candidate.edge_stretch >
      max(reference.edge_stretch * 1.10, reference.edge_stretch + 0.05))
    return "edge_stretch";
  if (isfinite(reference.min_pair_d2) &&
      candidate.min_pair_d2 + 1.0e-6 <
          min(reference.min_pair_d2 * 0.90, reference.min_pair_d2 - 0.05))
    return "min_pair_d2";
  return nullptr;
}

bool candidate_close_enough_to_skip_search(const PbctopoCutCandidate &candidate,
                                           const PbctopoCutCandidate &reference) {
  return candidate_skip_search_reject_reason(candidate, reference) == nullptr;
}

bool candidate_has_projection_clash_margin(const PbctopoCutCandidate &candidate,
                                           double clash_cutoff) {
  if (!candidate.valid || candidate.ultra_clash || !isfinite(candidate.min_pair_d2))
    return false;
  const double margin_distance =
      max(clash_cutoff * 1.10, clash_cutoff + 0.20);
  return candidate.min_pair_d2 >= margin_distance * margin_distance;
}

Int3 edge_delta(const PbctopoGraphEdge &edge) {
  return {edge.dx, edge.dy, edge.dz};
}

bool same_pair(const PbctopoGraphEdge &lhs, const PbctopoGraphEdge &rhs) {
  return lhs.a == rhs.a && lhs.b == rhs.b;
}

bool same_edge_option(const PbctopoGraphEdge &lhs, const PbctopoGraphEdge &rhs) {
  return same_pair(lhs, rhs) && lhs.dx == rhs.dx && lhs.dy == rhs.dy &&
         lhs.dz == rhs.dz;
}

bool build_component_tree_from_ranked_edges(
    size_t component_count, const vector<PbctopoGraphEdge> &ranked_edges,
    const vector<PbctopoGraphEdge> &forced_edges, vector<PbctopoGraphEdge> &tree_edges) {
  tree_edges.clear();
  ComponentGraphDsu dsu(component_count);
  auto add_edge = [&](const PbctopoGraphEdge &edge) {
    if (!edge.valid || edge.ultra_clash)
      return false;
    if (!dsu.unite(edge.a, edge.b))
      return false;
    tree_edges.push_back(edge);
    return true;
  };

  for (const auto &edge : forced_edges) {
    if (!add_edge(edge))
      return false;
  }

  for (const auto &edge : ranked_edges) {
    bool pair_forced = false;
    bool exact_forced = false;
    for (const auto &forced : forced_edges) {
      if (same_pair(edge, forced))
        pair_forced = true;
      if (same_edge_option(edge, forced)) {
        exact_forced = true;
        break;
      }
    }
    if (pair_forced && !exact_forced)
      continue;
    add_edge(edge);
    if (tree_edges.size() >= component_count - 1)
      break;
  }

  return tree_edges.size() + 1 == component_count;
}

bool solve_component_offsets_from_tree(
    size_t component_count, const vector<PbctopoGraphEdge> &tree_edges,
    vector<Int3> &offsets) {
  offsets.assign(component_count, {numeric_limits<int>::max(),
                                   numeric_limits<int>::max(),
                                   numeric_limits<int>::max()});
  vector<vector<pair<size_t, Int3>>> adjacency(component_count);
  for (const auto &edge : tree_edges) {
    const Int3 delta = edge_delta(edge);
    adjacency[edge.a].push_back({edge.b, delta});
    adjacency[edge.b].push_back({edge.a, {-delta.x, -delta.y, -delta.z}});
  }

  vector<char> seen(component_count, 0);
  queue<size_t> q;
  offsets[0] = {0, 0, 0};
  seen[0] = 1;
  q.push(0);
  while (!q.empty()) {
    const size_t current = q.front();
    q.pop();
    for (const auto &edge : adjacency[current]) {
      const size_t next = edge.first;
      if (seen[next] != 0)
        continue;
      offsets[next] = offsets[current] + edge.second;
      seen[next] = 1;
      q.push(next);
    }
  }
  return find(seen.begin(), seen.end(), 0) == seen.end();
}

bool graph_edge_has_zero_shift(const PbctopoGraphEdge &edge) {
  return edge.dx == 0 && edge.dy == 0 && edge.dz == 0;
}

bool graph_edge_matches_delta(const PbctopoGraphEdge &edge, const Int3 &delta) {
  return edge.dx == delta.x && edge.dy == delta.y && edge.dz == delta.z;
}

vector<PbctopoGraphEdge> build_component_graph_edge_candidates(
    const vector<PbctopoAtom> &atoms,
    const vector<vector<size_t>> &components, const vector<PbctopoShift> &shifts,
    const PbctopoRuntimeCutoffs &cutoffs, PbctopoTimingStats *timing,
    const vector<Int3> *warm_component_offsets) {
  size_t identity_shift_pos = 0;
  for (size_t shift_idx = 0; shift_idx < shifts.size(); ++shift_idx) {
    if (shifts[shift_idx].ix == 0 && shifts[shift_idx].iy == 0 &&
        shifts[shift_idx].iz == 0) {
      identity_shift_pos = shift_idx;
      break;
    }
  }

  const auto precompute_started = chrono::steady_clock::now();
  PbctopoLazyShiftedComponents shifted_components(
      atoms, components, shifts, cutoffs.clash_cutoff,
      cutoffs.pair_contact_cutoff, timing);
  if (timing != nullptr) {
    timing->component_graph_precompute_s +=
        chrono::duration<double>(chrono::steady_clock::now() - precompute_started)
            .count();
  }

  vector<PbctopoGraphEdge> ranked_edges;
  ranked_edges.reserve((components.size() * (components.size() - 1)) / 2);
  const auto pair_score_started = chrono::steady_clock::now();
  const double lazy_materialize_before =
      timing != nullptr ? timing->component_graph_lazy_materialize_s : 0.0;
  const double contact_cutoff2 =
      cutoffs.pair_contact_cutoff * cutoffs.pair_contact_cutoff;
  for (size_t lhs = 0; lhs < components.size(); ++lhs) {
    for (size_t rhs = lhs + 1; rhs < components.size(); ++rhs) {
      PbctopoPairShiftScore best_score;
      PbctopoPairShiftScore second_score;
      Int3 expected_delta;
      const Int3 *expected_delta_ptr = nullptr;
      if (warm_component_offsets != nullptr &&
          warm_component_offsets->size() == components.size()) {
        expected_delta = (*warm_component_offsets)[rhs] - (*warm_component_offsets)[lhs];
        expected_delta_ptr = &expected_delta;
      }
      const vector<size_t> ordered_shift_indices =
          build_pair_shift_order(shifts, expected_delta_ptr, timing);
      for (size_t shift_idx : ordered_shift_indices) {
        const auto &shift = shifts[shift_idx];
        const PbctopoPairShiftScore *incumbent =
            (best_score.valid && !best_score.ultra_clash) ? &best_score : nullptr;
        const double bbox_d2 = shifted_component_bbox_distance2(
            shifted_components.descriptor(lhs, identity_shift_pos),
            shifted_components.descriptor(rhs, shift_idx));
        PbctopoPairShiftScore score;
        if (bbox_d2 > contact_cutoff2) {
          score = make_bbox_only_pair_shift_score(shift, bbox_d2);
          if (timing != nullptr)
            ++timing->component_graph_bbox_only_scores;
        } else {
          score = score_component_pair_shift(
              shifted_components.materialize(lhs, identity_shift_pos),
              shifted_components.materialize(rhs, shift_idx), shift,
              cutoffs.clash_cutoff, cutoffs.pair_contact_cutoff, incumbent);
        }
        if (better_pair_shift_score(score, best_score)) {
          second_score = best_score;
          best_score = score;
        } else if (better_pair_shift_score(score, second_score)) {
          second_score = score;
        }
      }
      PbctopoGraphEdge edge;
      edge.a = lhs;
      edge.b = rhs;
      edge.dx = best_score.ix;
      edge.dy = best_score.iy;
      edge.dz = best_score.iz;
      edge.contact_atoms = best_score.contact_atoms;
      edge.soft_penalty = best_score.soft_penalty;
      edge.min_pair_d2 = best_score.min_pair_d2;
      edge.bbox_d2 = best_score.bbox_d2;
      edge.shift_norm2 = best_score.shift_norm2;
      edge.confidence = best_score.quality_score -
                        (second_score.valid ? second_score.quality_score
                                            : best_score.quality_score - 1.0);
      edge.ultra_clash = best_score.ultra_clash;
      edge.valid = best_score.valid;
      ranked_edges.push_back(edge);
    }
  }
  if (timing != nullptr) {
    const double lazy_materialize_delta =
        timing->component_graph_lazy_materialize_s - lazy_materialize_before;
    timing->component_graph_pair_score_s +=
        max(0.0, chrono::duration<double>(chrono::steady_clock::now() -
                                          pair_score_started)
                     .count() -
                     lazy_materialize_delta);
    timing->component_graph_edges += ranked_edges.size();
  }
  return ranked_edges;
}

bool solve_component_offsets_via_pair_graph(
    const vector<PbctopoAtom> &atoms, box Box,
    const vector<vector<size_t>> &components, const vector<PbctopoShift> &shifts,
    const vector<PbctopoLocalEdge> &local_edges,
    const unordered_set<uint64_t> &edge_set, vector<Vec3> &placed_points,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices,
    const PbctopoRuntimeCutoffs &cutoffs, PbctopoCutCandidate &best_cut,
    double &shift_norm2, PbctopoTimingStats *timing,
    PbctopoAnchorBackend anchor_backend = PbctopoAnchorBackend::Auto,
    const vector<Int3> *warm_component_offsets = nullptr,
    PbctopoWarmState *warm_state = nullptr,
    const PbctopoCutCandidate *warm_reference_candidate = nullptr,
    vector<Int3> *tree_offsets_out = nullptr,
    vector<Int3> *refined_offsets_out = nullptr,
    vector<PbctopoGraphEdge> *tree_edges_out = nullptr) {
  if (components.empty())
    return false;

  vector<PbctopoGraphEdge> ranked_edges = build_component_graph_edge_candidates(
      atoms, components, shifts, cutoffs, timing, warm_component_offsets);
  const auto tree_started = chrono::steady_clock::now();
  sort(ranked_edges.begin(), ranked_edges.end(), better_component_graph_edge);
  vector<PbctopoGraphEdge> tree_edges;
  if (!build_component_tree_from_ranked_edges(components.size(), ranked_edges, {},
                                              tree_edges))
    return false;
  vector<Int3> offsets;
  if (!solve_component_offsets_from_tree(components.size(), tree_edges, offsets))
    return false;
  if (timing != nullptr) {
    timing->component_graph_tree_s +=
        chrono::duration<double>(chrono::steady_clock::now() - tree_started)
            .count();
  }

  rebase_component_offsets(offsets, components);
  if (tree_offsets_out != nullptr)
    *tree_offsets_out = offsets;
  if (tree_edges_out != nullptr)
    *tree_edges_out = tree_edges;
  const auto refine_started = chrono::steady_clock::now();
  best_cut = improve_component_offsets_locally(atoms, components, Box, shifts,
                                               local_edges, edge_set,
                                               continuity_pairs, anchor_indices,
                                               offsets,
                                               cutoffs.clash_cutoff,
                                               cutoffs.continuity_cutoff,
                                               placed_points, anchor_backend,
                                               timing, warm_state,
                                               warm_reference_candidate);
  if (timing != nullptr) {
    timing->component_graph_refine_s +=
        chrono::duration<double>(chrono::steady_clock::now() - refine_started)
            .count();
  }
  if (refined_offsets_out != nullptr)
    *refined_offsets_out = offsets;
  shift_norm2 = best_cut.shift_norm2;
  return best_cut.valid && placed_points.size() == atoms.size();
}

struct PbctopoGraphSeedSolveResult {
  const char *seed_name = "baseline";
  vector<Int3> offsets;
  PbctopoGraphSeedQuickScore quick_score;
  bool valid = false;
};

bool solve_component_offsets_from_ranked_edges(
    size_t component_count, const vector<vector<size_t>> &components,
    const vector<PbctopoGraphEdge> &ranked_edges, vector<Int3> &offsets) {
  vector<PbctopoGraphEdge> tree_edges;
  if (!build_component_tree_from_ranked_edges(component_count, ranked_edges, {},
                                              tree_edges))
    return false;
  if (!solve_component_offsets_from_tree(component_count, tree_edges, offsets))
    return false;
  rebase_component_offsets(offsets, components);
  return true;
}

vector<PbctopoGraphSeedSolveResult> build_large_component_graph_seed_candidates(
    const vector<PbctopoAtom> &atoms, box Box,
    const vector<vector<size_t>> &components, const vector<PbctopoShift> &shifts,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices, const PbctopoRuntimeCutoffs &cutoffs,
    const vector<Int3> *warm_component_offsets, PbctopoTimingStats *timing) {
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
    if (!solve_component_offsets_from_ranked_edges(components.size(), components,
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
        cutoffs.continuity_cutoff, seed_shift_norm2);
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
                  if (lhs.contact_atoms != rhs.contact_atoms)
                    return lhs.contact_atoms > rhs.contact_atoms;
                  if (fabs(lhs.confidence - rhs.confidence) > 1.0e-9)
                    return lhs.confidence > rhs.confidence;
                  if (fabs(lhs.soft_penalty - rhs.soft_penalty) > 1.0e-9)
                    return lhs.soft_penalty < rhs.soft_penalty;
                  if (fabs(lhs.min_pair_d2 - rhs.min_pair_d2) > 1.0e-9)
                    return lhs.min_pair_d2 > rhs.min_pair_d2;
                  if (fabs(lhs.bbox_d2 - rhs.bbox_d2) > 1.0e-9)
                    return lhs.bbox_d2 < rhs.bbox_d2;
                  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-9)
                    return lhs.shift_norm2 < rhs.shift_norm2;
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
                  if (lhs.contact_atoms != rhs.contact_atoms)
                    return lhs.contact_atoms > rhs.contact_atoms;
                  if (fabs(lhs.confidence - rhs.confidence) > 1.0e-9)
                    return lhs.confidence > rhs.confidence;
                  if (fabs(lhs.soft_penalty - rhs.soft_penalty) > 1.0e-9)
                    return lhs.soft_penalty < rhs.soft_penalty;
                  if (fabs(lhs.min_pair_d2 - rhs.min_pair_d2) > 1.0e-9)
                    return lhs.min_pair_d2 > rhs.min_pair_d2;
                  if (fabs(lhs.bbox_d2 - rhs.bbox_d2) > 1.0e-9)
                    return lhs.bbox_d2 < rhs.bbox_d2;
                  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-9)
                    return lhs.shift_norm2 < rhs.shift_norm2;
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

vector<size_t> build_atom_component_index(
    size_t atom_count, const vector<vector<size_t>> &components) {
  vector<size_t> atom_component_index(atom_count, numeric_limits<size_t>::max());
  for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
    for (size_t atom_idx : components[comp_idx]) {
      if (atom_idx < atom_component_index.size())
        atom_component_index[atom_idx] = comp_idx;
    }
  }
  return atom_component_index;
}

PbctopoStage1EdgeContext build_stage1_edge_context(
    const vector<PbctopoLocalEdge> &local_edges,
    const vector<size_t> &atom_component_index, size_t moving_component_idx,
    const vector<Vec3> &placed_points) {
  PbctopoStage1EdgeContext context;
  if (atom_component_index.empty() || placed_points.empty())
    return context;

  context.variable_edge_indices.reserve(local_edges.size() / 4 + 1);
  for (size_t edge_idx = 0; edge_idx < local_edges.size(); ++edge_idx) {
    const auto &edge = local_edges[edge_idx];
    const size_t comp_a = (edge.a < atom_component_index.size())
                              ? atom_component_index[edge.a]
                              : numeric_limits<size_t>::max();
    const size_t comp_b = (edge.b < atom_component_index.size())
                              ? atom_component_index[edge.b]
                              : numeric_limits<size_t>::max();
    const bool variable =
        comp_a != comp_b &&
        (comp_a == moving_component_idx || comp_b == moving_component_idx);
    if (variable) {
      context.variable_edge_indices.push_back(edge_idx);
      continue;
    }

    const Vec3 delta = placed_points[edge.a] - placed_points[edge.b];
    const double d2 = norm2(delta);
    if (d2 > edge.break_cutoff_d2)
      ++context.fixed_broken_edges;
    if (d2 > edge.ref_d2)
      context.fixed_edge_stretch += (d2 - edge.ref_d2) * edge.inv_ref_d2;
  }
  return context;
}

PbctopoEdgeScanSoA
build_edge_scan_soa(const vector<PbctopoLocalEdge> &local_edges) {
  PbctopoEdgeScanSoA soa;
  soa.a.reserve(local_edges.size());
  soa.b.reserve(local_edges.size());
  soa.ref_d2.reserve(local_edges.size());
  soa.inv_ref_d2.reserve(local_edges.size());
  soa.break_cutoff_d2.reserve(local_edges.size());
  for (const auto &edge : local_edges) {
    soa.a.push_back(edge.a);
    soa.b.push_back(edge.b);
    soa.ref_d2.push_back(edge.ref_d2);
    soa.inv_ref_d2.push_back(edge.inv_ref_d2);
    soa.break_cutoff_d2.push_back(edge.break_cutoff_d2);
  }
  return soa;
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

int compare_candidate_prefix_before_edge_stretch(
    const PbctopoCutCandidate &lhs, const PbctopoCutCandidate &rhs,
    PbctopoRefineScorePhase phase) {
  if (lhs.broken_edges != rhs.broken_edges)
    return (lhs.broken_edges < rhs.broken_edges) ? -1 : 1;
  const bool lhs_has_cont =
      isfinite(lhs.continuity_max_d2) && isfinite(lhs.continuity_path);
  const bool rhs_has_cont =
      isfinite(rhs.continuity_max_d2) && isfinite(rhs.continuity_path);
  if (lhs_has_cont != rhs_has_cont)
    return lhs_has_cont ? -1 : 1;
  if (lhs_has_cont &&
      fabs(lhs.continuity_max_d2 - rhs.continuity_max_d2) > 1.0e-6) {
    return (lhs.continuity_max_d2 < rhs.continuity_max_d2) ? -1 : 1;
  }
  if (lhs_has_cont &&
      fabs(lhs.continuity_path - rhs.continuity_path) > 1.0e-6) {
    return (lhs.continuity_path < rhs.continuity_path) ? -1 : 1;
  }
  const bool lhs_has_anchor =
      isfinite(lhs.anchor_max_d2) && isfinite(lhs.anchor_mst2);
  const bool rhs_has_anchor =
      isfinite(rhs.anchor_max_d2) && isfinite(rhs.anchor_mst2);
  if (lhs_has_anchor != rhs_has_anchor)
    return lhs_has_anchor ? -1 : 1;
  if (lhs_has_anchor && fabs(lhs.anchor_max_d2 - rhs.anchor_max_d2) > 1.0e-6)
    return (lhs.anchor_max_d2 < rhs.anchor_max_d2) ? -1 : 1;
  if (lhs_has_anchor && fabs(lhs.anchor_mst2 - rhs.anchor_mst2) > 1.0e-6)
    return (lhs.anchor_mst2 < rhs.anchor_mst2) ? -1 : 1;
  const int principal_max_cmp =
      compare_principal_rg2(lhs.principal_rg2_max, rhs.principal_rg2_max);
  if (principal_max_cmp != 0)
    return principal_max_cmp;
  const int principal_mid_cmp =
      compare_principal_rg2(lhs.principal_rg2_mid, rhs.principal_rg2_mid);
  if (principal_mid_cmp != 0)
    return principal_mid_cmp;
  if (phase == PbctopoRefineScorePhase::CoarseStage1) {
    const int principal_min_cmp =
        compare_principal_rg2(lhs.principal_rg2_min, rhs.principal_rg2_min);
    if (principal_min_cmp != 0)
      return principal_min_cmp;
  }
  return 0;
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
           if (fabs(lhs.d2 - rhs.d2) > 1.0e-9)
             return lhs.d2 < rhs.d2;
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
        const double dsx = wrapped_atoms[nb].sx - wrapped_atoms[current].sx -
                           round(wrapped_atoms[nb].sx - wrapped_atoms[current].sx);
        const double dsy = wrapped_atoms[nb].sy - wrapped_atoms[current].sy -
                           round(wrapped_atoms[nb].sy - wrapped_atoms[current].sy);
        const double dsz = wrapped_atoms[nb].sz - wrapped_atoms[current].sz -
                           round(wrapped_atoms[nb].sz - wrapped_atoms[current].sz);
        unwrapped[nb].sx = unwrapped[current].sx + dsx;
        unwrapped[nb].sy = unwrapped[current].sy + dsy;
        unwrapped[nb].sz = unwrapped[current].sz + dsz;
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

PbctopoCutCandidate
component_candidate_to_cut_candidate(const PbctopoComponentCandidate &candidate) {
  PbctopoCutCandidate cut;
  cut.continuity_path = candidate.continuity_path;
  cut.continuity_max_d2 = candidate.continuity_max_d2;
  cut.anchor_mst2 = candidate.anchor_mst2;
  cut.anchor_max_d2 = candidate.anchor_max_d2;
  cut.principal_rg2_max = candidate.principal_rg2_max;
  cut.principal_rg2_mid = candidate.principal_rg2_mid;
  cut.principal_rg2_min = candidate.principal_rg2_min;
  cut.bbox_volume = candidate.bbox_volume;
  cut.rg2 = candidate.rg2;
  cut.shape_kappa2 = candidate.shape_kappa2;
  cut.shape_asphericity = candidate.shape_asphericity;
  cut.shape_acylindricity = candidate.shape_acylindricity;
  cut.shift_norm2 = candidate.shift_norm2;
  cut.min_pair_d2 = candidate.min_pair_d2;
  cut.clash_atom_a = candidate.clash_atom_a;
  cut.clash_atom_b = candidate.clash_atom_b;
  cut.ultra_clash = candidate.ultra_clash;
  cut.valid = candidate.valid;
  return cut;
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

PbctopoComponentPairShiftEntry compute_component_pair_shift_entry(
    const vector<size_t> &lhs_component, const PbctopoShiftedComponent &lhs,
    const vector<size_t> &rhs_component, const PbctopoShiftedComponent &rhs,
    double clash_cutoff) {
  PbctopoComponentPairShiftEntry entry;
  entry.bbox_d2 = shifted_component_bbox_distance2(lhs, rhs);
  if (!component_bbox_may_clash(lhs, rhs, clash_cutoff))
    return entry;

  const double cutoff2 = clash_cutoff * clash_cutoff;
  const vector<size_t> *scan_component = &lhs_component;
  const vector<size_t> *grid_component = &rhs_component;
  const PbctopoShiftedComponent *scan_shifted = &lhs;
  const PbctopoShiftedComponent *grid_shifted = &rhs;
  if (lhs_component.size() > rhs_component.size()) {
    scan_component = &rhs_component;
    grid_component = &lhs_component;
    scan_shifted = &rhs;
    grid_shifted = &lhs;
  }

  for (size_t local_idx = 0; local_idx < scan_component->size(); ++local_idx) {
    const size_t atom_idx = (*scan_component)[local_idx];
    const auto &point = scan_shifted->points[local_idx];
    const auto &key = scan_shifted->clash_keys[local_idx];
    const int cx = key.x;
    const int cy = key.y;
    const int cz = key.z;
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          const auto it =
              grid_shifted->clash_grid.find(pack_cell_key(cx + dx, cy + dy, cz + dz));
          if (it == grid_shifted->clash_grid.end())
            continue;
          for (size_t other_local_idx : it->second) {
            const size_t other_idx = (*grid_component)[other_local_idx];
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
    const vector<vector<size_t>> &components,
    const vector<vector<PbctopoShiftedComponent>> &shifted_components,
    double clash_cutoff, PbctopoTimingStats *timing = nullptr) {
  using clock = chrono::steady_clock;
  const auto started = clock::now();
  PbctopoComponentPairShiftTable table;
  const size_t component_count = components.size();
  const size_t shift_count =
      shifted_components.empty() ? 0 : shifted_components.front().size();
  table.reset(component_count, shift_count);

  for (size_t lhs_comp = 0; lhs_comp < component_count; ++lhs_comp) {
    for (size_t rhs_comp = lhs_comp + 1; rhs_comp < component_count; ++rhs_comp) {
      for (size_t lhs_shift = 0; lhs_shift < shift_count; ++lhs_shift) {
        for (size_t rhs_shift = 0; rhs_shift < shift_count; ++rhs_shift) {
          const auto entry = compute_component_pair_shift_entry(
              components[lhs_comp], shifted_components[lhs_comp][lhs_shift],
              components[rhs_comp], shifted_components[rhs_comp][rhs_shift], clash_cutoff);
          table.at(lhs_comp, lhs_shift, rhs_comp, rhs_shift) = entry;
          table.at(rhs_comp, rhs_shift, lhs_comp, lhs_shift) = entry;
          if (!entry.hard_clash) {
            table.compat_mask(lhs_comp, lhs_shift, rhs_comp) |=
                (1u << rhs_shift);
            table.compat_mask(rhs_comp, rhs_shift, lhs_comp) |=
                (1u << lhs_shift);
          }
        }
      }
    }
  }

  if (timing != nullptr) {
    ++timing->component_pair_table_builds;
    timing->component_pair_table_s +=
        chrono::duration<double>(clock::now() - started).count();
  }
  return table;
}

uint32_t component_search_full_shift_mask(size_t shift_count) {
  if (shift_count >= 32)
    return numeric_limits<uint32_t>::max();
  return shift_count == 0 ? 0u : ((1u << shift_count) - 1u);
}

size_t find_component_search_shift_index(const vector<PbctopoShift> &shifts,
                                         const Int3 &offset) {
  for (size_t shift_idx = 0; shift_idx < shifts.size(); ++shift_idx) {
    if (shifts[shift_idx].ix == offset.x && shifts[shift_idx].iy == offset.y &&
        shifts[shift_idx].iz == offset.z) {
      return shift_idx;
    }
  }
  return shifts.size();
}

vector<size_t> build_component_search_order(
    const vector<vector<size_t>> &components, const vector<size_t> &anchor_indices,
    const vector<size_t> &atom_component_index,
    const PbctopoComponentPairShiftTable &pair_table, size_t identity_shift_pos) {
  vector<size_t> order;
  if (components.empty())
    return order;
  order.reserve(components.size());
  order.push_back(0);
  if (components.size() == 1)
    return order;

  vector<size_t> anchor_counts(components.size(), 0);
  for (size_t anchor_idx : anchor_indices) {
    if (anchor_idx < atom_component_index.size()) {
      const size_t comp_idx = atom_component_index[anchor_idx];
      if (comp_idx < anchor_counts.size())
        ++anchor_counts[comp_idx];
    }
  }

  vector<size_t> remaining;
  remaining.reserve(components.size() - 1);
  for (size_t comp_idx = 1; comp_idx < components.size(); ++comp_idx)
    remaining.push_back(comp_idx);

  stable_sort(remaining.begin(), remaining.end(),
              [&](size_t lhs, size_t rhs) {
                if (components[lhs].size() != components[rhs].size())
                  return components[lhs].size() > components[rhs].size();
                if (anchor_counts[lhs] != anchor_counts[rhs])
                  return anchor_counts[lhs] > anchor_counts[rhs];

                const size_t lhs_compatible_shifts = static_cast<size_t>(popcount(
                    pair_table.compat_mask(0, identity_shift_pos, lhs)));
                const size_t rhs_compatible_shifts = static_cast<size_t>(popcount(
                    pair_table.compat_mask(0, identity_shift_pos, rhs)));
                double lhs_best_min_pair = numeric_limits<double>::infinity();
                double rhs_best_min_pair = numeric_limits<double>::infinity();
                for (size_t shift_idx = 0; shift_idx < pair_table.shift_count; ++shift_idx) {
                  const auto &lhs_entry =
                      pair_table.at(0, identity_shift_pos, lhs, shift_idx);
                  const auto &rhs_entry =
                      pair_table.at(0, identity_shift_pos, rhs, shift_idx);
                  lhs_best_min_pair = min(lhs_best_min_pair, lhs_entry.min_pair_d2);
                  rhs_best_min_pair = min(rhs_best_min_pair, rhs_entry.min_pair_d2);
                }
                if (lhs_compatible_shifts != rhs_compatible_shifts)
                  return lhs_compatible_shifts < rhs_compatible_shifts;
                if (fabs(lhs_best_min_pair - rhs_best_min_pair) > 1.0e-9)
                  return lhs_best_min_pair < rhs_best_min_pair;
                return lhs < rhs;
              });

  order.insert(order.end(), remaining.begin(), remaining.end());
  return order;
}

bool component_candidate_cannot_beat_incumbent_prefix(
    const PbctopoComponentCandidate &candidate,
    const PbctopoComponentCandidate *incumbent,
    PbctopoComponentCandidatePrefixStage stage) {
  if (incumbent == nullptr || !incumbent->valid || incumbent->ultra_clash ||
      candidate.ultra_clash) {
    return false;
  }

  const bool candidate_has_cont =
      isfinite(candidate.continuity_max_d2) && isfinite(candidate.continuity_path);
  const bool incumbent_has_cont =
      isfinite(incumbent->continuity_max_d2) && isfinite(incumbent->continuity_path);
  if (candidate_has_cont != incumbent_has_cont)
    return !candidate_has_cont;
  if (candidate_has_cont) {
    if (candidate.continuity_max_d2 > incumbent->continuity_max_d2 + 1.0e-6)
      return true;
    if (candidate.continuity_max_d2 + 1.0e-6 < incumbent->continuity_max_d2)
      return false;
    if (candidate.continuity_path > incumbent->continuity_path + 1.0e-6)
      return true;
    if (candidate.continuity_path + 1.0e-6 < incumbent->continuity_path)
      return false;
  }
  if (stage == PbctopoComponentCandidatePrefixStage::Continuity)
    return false;

  const bool candidate_has_anchor =
      isfinite(candidate.anchor_max_d2) && isfinite(candidate.anchor_mst2);
  const bool incumbent_has_anchor =
      isfinite(incumbent->anchor_max_d2) && isfinite(incumbent->anchor_mst2);
  if (candidate_has_anchor != incumbent_has_anchor)
    return !candidate_has_anchor;
  if (candidate_has_anchor) {
    if (candidate.anchor_max_d2 > incumbent->anchor_max_d2 + 1.0e-6)
      return true;
    if (candidate.anchor_max_d2 + 1.0e-6 < incumbent->anchor_max_d2)
      return false;
    if (candidate.anchor_mst2 > incumbent->anchor_mst2 + 1.0e-6)
      return true;
    if (candidate.anchor_mst2 + 1.0e-6 < incumbent->anchor_mst2)
      return false;
  }
  if (stage == PbctopoComponentCandidatePrefixStage::Anchor)
    return false;

  const int principal_max_cmp =
      compare_principal_rg2(candidate.principal_rg2_max, incumbent->principal_rg2_max);
  if (principal_max_cmp > 0)
    return true;
  if (principal_max_cmp < 0)
    return false;
  const int principal_mid_cmp =
      compare_principal_rg2(candidate.principal_rg2_mid, incumbent->principal_rg2_mid);
  if (principal_mid_cmp > 0)
    return true;
  if (principal_mid_cmp < 0)
    return false;
  const int principal_min_cmp =
      compare_principal_rg2(candidate.principal_rg2_min, incumbent->principal_rg2_min);
  if (principal_min_cmp > 0)
    return true;
  if (principal_min_cmp < 0)
    return false;
  if (stage == PbctopoComponentCandidatePrefixStage::Principal)
    return false;

  if (candidate.component_mst2 > incumbent->component_mst2 + 1.0e-6)
    return true;
  if (candidate.component_mst2 + 1.0e-6 < incumbent->component_mst2)
    return false;
  const int component_kappa_cmp =
      compare_shape_kappa(candidate.component_kappa2,
                          incumbent->component_kappa2);
  if (component_kappa_cmp > 0)
    return true;
  if (component_kappa_cmp < 0)
    return false;
  const int component_asph_cmp =
      compare_shape_moment(candidate.component_asphericity,
                           incumbent->component_asphericity);
  if (component_asph_cmp > 0)
    return true;
  if (component_asph_cmp < 0)
    return false;
  const int component_acyl_cmp =
      compare_shape_moment(candidate.component_acylindricity,
                           incumbent->component_acylindricity);
  if (component_acyl_cmp > 0)
    return true;
  if (component_acyl_cmp < 0)
    return false;
  if (stage == PbctopoComponentCandidatePrefixStage::ComponentMetrics)
    return false;

  const int shape_kappa_cmp =
      compare_shape_kappa(candidate.shape_kappa2, incumbent->shape_kappa2);
  if (shape_kappa_cmp > 0)
    return true;
  if (shape_kappa_cmp < 0)
    return false;
  const int shape_asph_cmp =
      compare_shape_moment(candidate.shape_asphericity,
                           incumbent->shape_asphericity);
  if (shape_asph_cmp > 0)
    return true;
  if (shape_asph_cmp < 0)
    return false;
  const int shape_acyl_cmp =
      compare_shape_moment(candidate.shape_acylindricity,
                           incumbent->shape_acylindricity);
  if (shape_acyl_cmp > 0)
    return true;
  if (shape_acyl_cmp < 0)
    return false;
  if (candidate.shift_norm2 > incumbent->shift_norm2 + 1.0e-6)
    return true;
  if (candidate.shift_norm2 + 1.0e-6 < incumbent->shift_norm2)
    return false;
  if (candidate.min_pair_d2 + 1.0e-9 < incumbent->min_pair_d2)
    return true;
  return false;
}

PbctopoComponentCandidate score_component_assignment(
    const vector<PbctopoAtom> &atoms, const vector<Vec3> &placed_points,
    const vector<vector<size_t>> &components, double min_pair_d2,
    size_t clash_atom_a, size_t clash_atom_b,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices, double continuity_cutoff,
    double shift_norm2, const PbctopoComponentCandidate *incumbent = nullptr,
    PbctopoTimingStats *timing = nullptr) {
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
  if (!anchor_graph_metrics_ok(placed_points, anchor_indices,
                               candidate.anchor_mst2,
                               candidate.anchor_max_d2, incumbent_cut)) {
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

bool better_component_candidate(const PbctopoComponentCandidate &lhs,
                                const PbctopoComponentCandidate &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  if (lhs.ultra_clash != rhs.ultra_clash)
    return rhs.ultra_clash;
  const bool lhs_has_cont =
      isfinite(lhs.continuity_max_d2) && isfinite(lhs.continuity_path);
  const bool rhs_has_cont =
      isfinite(rhs.continuity_max_d2) && isfinite(rhs.continuity_path);
  if (lhs_has_cont != rhs_has_cont)
    return lhs_has_cont;
  if (lhs_has_cont &&
      fabs(lhs.continuity_max_d2 - rhs.continuity_max_d2) > 1.0e-6)
    return lhs.continuity_max_d2 < rhs.continuity_max_d2;
  if (lhs_has_cont &&
      fabs(lhs.continuity_path - rhs.continuity_path) > 1.0e-6)
    return lhs.continuity_path < rhs.continuity_path;
  const bool lhs_has_anchor =
      isfinite(lhs.anchor_max_d2) && isfinite(lhs.anchor_mst2);
  const bool rhs_has_anchor =
      isfinite(rhs.anchor_max_d2) && isfinite(rhs.anchor_mst2);
  if (lhs_has_anchor != rhs_has_anchor)
    return lhs_has_anchor;
  if (lhs_has_anchor && fabs(lhs.anchor_max_d2 - rhs.anchor_max_d2) > 1.0e-6)
    return lhs.anchor_max_d2 < rhs.anchor_max_d2;
  if (lhs_has_anchor && fabs(lhs.anchor_mst2 - rhs.anchor_mst2) > 1.0e-6)
    return lhs.anchor_mst2 < rhs.anchor_mst2;
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
  if (fabs(lhs.component_mst2 - rhs.component_mst2) > 1.0e-6)
    return lhs.component_mst2 < rhs.component_mst2;
  const int component_kappa_cmp =
      compare_shape_kappa(lhs.component_kappa2, rhs.component_kappa2);
  if (component_kappa_cmp != 0)
    return component_kappa_cmp < 0;
  const int component_asph_cmp =
      compare_shape_moment(lhs.component_asphericity, rhs.component_asphericity);
  if (component_asph_cmp != 0)
    return component_asph_cmp < 0;
  const int component_acyl_cmp = compare_shape_moment(
      lhs.component_acylindricity, rhs.component_acylindricity);
  if (component_acyl_cmp != 0)
    return component_acyl_cmp < 0;
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
  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-6)
    return lhs.shift_norm2 < rhs.shift_norm2;
  if (fabs(lhs.min_pair_d2 - rhs.min_pair_d2) > 1.0e-9)
    return lhs.min_pair_d2 > rhs.min_pair_d2;
  return false;
}

struct PbctopoComponentSearchShiftOption {
  size_t shift_idx = 0;
  double next_min_pair_d2 = numeric_limits<double>::infinity();
  double ordering_min_pair_d2 = numeric_limits<double>::infinity();
  double shift_norm2 = numeric_limits<double>::infinity();
  size_t future_remaining_shifts = numeric_limits<size_t>::max();
  size_t clash_atom_a = numeric_limits<size_t>::max();
  size_t clash_atom_b = numeric_limits<size_t>::max();
  bool warm_priority = false;
};

struct PbctopoBeamState {
  vector<size_t> selected_shift;
  vector<uint32_t> allowed_masks;
  double shift_norm2 = 0.0;
  double min_pair_d2 = numeric_limits<double>::infinity();
  size_t clash_atom_a = numeric_limits<size_t>::max();
  size_t clash_atom_b = numeric_limits<size_t>::max();
  size_t future_remaining_shifts = numeric_limits<size_t>::max();
  size_t warm_hits = 0;
};

vector<Int3> normalize_search_warm_offsets(
    const vector<Int3> &warm_offsets, const vector<vector<size_t>> &components) {
  if (warm_offsets.empty())
    return {};
  vector<Int3> rebased = warm_offsets;
  rebase_component_offsets(rebased, components);
  vector<Int3> normalized;
  normalized.reserve(rebased.size());
  const Int3 root = rebased.front();
  for (const auto &offset : rebased) {
    normalized.push_back(
        {offset.x - root.x, offset.y - root.y, offset.z - root.z});
  }
  return normalized;
}

vector<PbctopoComponentSearchShiftOption> build_component_search_shift_order(
    const vector<PbctopoShift> &shifts, size_t actual_comp_idx, size_t depth,
    const vector<size_t> &search_order, const vector<size_t> &selected_shift,
    const vector<uint32_t> &allowed_masks,
    const PbctopoComponentPairShiftTable &pair_table, double current_min_pair_d2,
    size_t clash_atom_a, size_t clash_atom_b,
    const vector<Int3> *normalized_warm_offsets = nullptr,
    PbctopoTimingStats *timing = nullptr) {
  vector<PbctopoComponentSearchShiftOption> ordered;
  if (actual_comp_idx >= allowed_masks.size())
    return ordered;
  uint32_t current_mask = allowed_masks[actual_comp_idx];
  if (current_mask == 0u) {
    if (timing != nullptr) {
      ++timing->component_pair_table_prunes;
      ++timing->component_search_zero_current_mask_prunes;
    }
    return ordered;
  }
  ordered.reserve(static_cast<size_t>(popcount(current_mask)));

  Int3 warm_delta;
  const bool has_warm_delta =
      (depth > 0 && normalized_warm_offsets != nullptr &&
       actual_comp_idx < normalized_warm_offsets->size())
          ? (warm_delta = (*normalized_warm_offsets)[actual_comp_idx], true)
          : false;

  while (current_mask != 0u) {
    const size_t shift_idx = static_cast<size_t>(countr_zero(current_mask));
    current_mask &= (current_mask - 1u);
    if (shift_idx >= shifts.size())
      continue;
    PbctopoComponentSearchShiftOption option;
    option.shift_idx = shift_idx;
    option.next_min_pair_d2 = current_min_pair_d2;
    option.ordering_min_pair_d2 = numeric_limits<double>::infinity();
    option.future_remaining_shifts = 0;
    option.shift_norm2 = static_cast<double>(shifts[shift_idx].ix * shifts[shift_idx].ix +
                                             shifts[shift_idx].iy * shifts[shift_idx].iy +
                                             shifts[shift_idx].iz * shifts[shift_idx].iz);
    option.clash_atom_a = clash_atom_a;
    option.clash_atom_b = clash_atom_b;
    option.warm_priority =
        has_warm_delta && shifts[shift_idx].ix == warm_delta.x &&
        shifts[shift_idx].iy == warm_delta.y && shifts[shift_idx].iz == warm_delta.z;

    for (size_t prev_depth = 0; prev_depth < depth; ++prev_depth) {
      const size_t prev_comp_idx = search_order[prev_depth];
      const size_t prev_shift_idx = selected_shift[prev_comp_idx];
      if (prev_shift_idx >= pair_table.shift_count)
        continue;
      const auto &entry =
          pair_table.at(actual_comp_idx, shift_idx, prev_comp_idx, prev_shift_idx);
      option.ordering_min_pair_d2 =
          min(option.ordering_min_pair_d2, entry.min_pair_d2);
      if (entry.min_pair_d2 < option.next_min_pair_d2) {
        option.next_min_pair_d2 = entry.min_pair_d2;
        option.clash_atom_a = entry.clash_atom_a;
        option.clash_atom_b = entry.clash_atom_b;
      }
    }

    bool future_zero_mask = false;
    for (size_t future_depth = depth + 1; future_depth < search_order.size();
         ++future_depth) {
      const size_t future_comp_idx = search_order[future_depth];
      if (future_comp_idx >= allowed_masks.size())
        continue;
      const uint32_t next_mask =
          allowed_masks[future_comp_idx] &
          pair_table.compat_mask(actual_comp_idx, shift_idx, future_comp_idx);
      if (next_mask == 0u) {
        future_zero_mask = true;
        break;
      }
      option.future_remaining_shifts += static_cast<size_t>(popcount(next_mask));
    }
    if (future_zero_mask) {
      if (timing != nullptr)
        ++timing->component_pair_table_prunes;
      if (timing != nullptr)
        ++timing->component_search_zero_future_mask_prunes;
      continue;
    }
    ordered.push_back(option);
  }

  sort(ordered.begin(), ordered.end(),
       [](const PbctopoComponentSearchShiftOption &lhs,
          const PbctopoComponentSearchShiftOption &rhs) {
         if (lhs.warm_priority != rhs.warm_priority)
           return lhs.warm_priority;
         if (lhs.future_remaining_shifts != rhs.future_remaining_shifts)
           return lhs.future_remaining_shifts < rhs.future_remaining_shifts;
         if (fabs(lhs.ordering_min_pair_d2 - rhs.ordering_min_pair_d2) > 1.0e-9)
           return lhs.ordering_min_pair_d2 > rhs.ordering_min_pair_d2;
         if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-9)
           return lhs.shift_norm2 < rhs.shift_norm2;
         return lhs.shift_idx < rhs.shift_idx;
       });
  if (!ordered.empty() && ordered.front().warm_priority && timing != nullptr)
    ++timing->component_search_shift_warm_hits;
  return ordered;
}

bool better_component_beam_state(const PbctopoBeamState &lhs,
                                 const PbctopoBeamState &rhs) {
  const bool lhs_has_pair = isfinite(lhs.min_pair_d2);
  const bool rhs_has_pair = isfinite(rhs.min_pair_d2);
  if (lhs_has_pair != rhs_has_pair)
    return lhs_has_pair;
  if (lhs_has_pair && fabs(lhs.min_pair_d2 - rhs.min_pair_d2) > 1.0e-9)
    return lhs.min_pair_d2 > rhs.min_pair_d2;
  if (lhs.warm_hits != rhs.warm_hits)
    return lhs.warm_hits > rhs.warm_hits;
  if (lhs.future_remaining_shifts != rhs.future_remaining_shifts)
    return lhs.future_remaining_shifts < rhs.future_remaining_shifts;
  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-9)
    return lhs.shift_norm2 < rhs.shift_norm2;
  return lhs.selected_shift < rhs.selected_shift;
}

bool build_component_points_from_selected_shifts(
    const vector<PbctopoAtom> &atoms, const vector<vector<size_t>> &components,
    const vector<vector<PbctopoShiftedComponent>> &shifted_components,
    const vector<size_t> &selected_shift, vector<Vec3> &placed_points,
    double &shift_norm2) {
  if (components.size() != selected_shift.size())
    return false;
  placed_points.resize(atoms.size());
  shift_norm2 = 0.0;
  for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
    const size_t shift_idx = selected_shift[comp_idx];
    if (shift_idx >= shifted_components[comp_idx].size())
      return false;
    const auto &component = components[comp_idx];
    const auto &shifted = shifted_components[comp_idx][shift_idx];
    if (shifted.points.size() != component.size())
      return false;
    shift_norm2 += shifted.weighted_shift_norm2;
    for (size_t local_idx = 0; local_idx < component.size(); ++local_idx)
      placed_points[component[local_idx]] = shifted.points[local_idx];
  }
  return true;
}

PbctopoComponentCandidate search_component_shift_candidate_beam(
    const vector<PbctopoAtom> &atoms, box Box,
    const vector<vector<size_t>> &components, const vector<PbctopoShift> &shifts,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices,
    const PbctopoRuntimeCutoffs &cutoffs, vector<Vec3> &best_points_out,
    size_t beam_width,
    const vector<Int3> *warm_component_offsets = nullptr,
    const PbctopoCutCandidate *incumbent_cut = nullptr,
    const vector<Vec3> *incumbent_points = nullptr,
    PbctopoTimingStats *timing = nullptr) {
  (void)Box;
  beam_width = max<size_t>(1, min(beam_width, kPbctopoBeamWidthMax));
  PbctopoComponentCandidate best;
  if (components.empty())
    return best;

  vector<Vec3> placed_points(atoms.size());
  best_points_out.clear();
  best_points_out.reserve(atoms.size());

  const auto shifted_components =
      precompute_shifted_components(atoms, components, shifts,
                                    cutoffs.clash_cutoff,
                                    cutoffs.pair_contact_cutoff);

  size_t identity_shift_pos = 0;
  for (size_t shift_idx = 0; shift_idx < shifts.size(); ++shift_idx) {
    if (shifts[shift_idx].ix == 0 && shifts[shift_idx].iy == 0 &&
        shifts[shift_idx].iz == 0) {
      identity_shift_pos = shift_idx;
      break;
    }
  }

  if (incumbent_cut != nullptr && incumbent_cut->valid &&
      !incumbent_cut->ultra_clash && incumbent_points != nullptr &&
      incumbent_points->size() == atoms.size()) {
    best = build_component_candidate_seed(*incumbent_cut, components, *incumbent_points);
    if (timing != nullptr)
      ++timing->component_search_seed_graph;
  } else if (timing != nullptr) {
    ++timing->component_search_seed_none;
  }

  const PbctopoComponentPairShiftTable pair_table =
      build_component_pair_shift_table(components, shifted_components,
                                       cutoffs.clash_cutoff, timing);
  const vector<size_t> atom_component_index =
      build_atom_component_index(atoms.size(), components);
  const vector<size_t> search_order = build_component_search_order(
      components, anchor_indices, atom_component_index, pair_table, identity_shift_pos);
  const vector<Int3> normalized_warm_offsets =
      (warm_component_offsets != nullptr &&
       warm_component_offsets->size() == components.size())
          ? normalize_search_warm_offsets(*warm_component_offsets, components)
          : vector<Int3>{};
  const vector<Int3> *normalized_warm_offsets_ptr =
      normalized_warm_offsets.empty() ? nullptr : &normalized_warm_offsets;

  const uint32_t full_shift_mask = component_search_full_shift_mask(shifts.size());
  vector<uint32_t> initial_allowed_masks(components.size(), full_shift_mask);
  vector<size_t> initial_selected_shift(components.size(), shifts.size());
  if (!search_order.empty() && identity_shift_pos < 32) {
    const size_t root_comp_idx = search_order.front();
    initial_allowed_masks[root_comp_idx] = (1u << identity_shift_pos);
    initial_selected_shift[root_comp_idx] = identity_shift_pos;
    for (size_t depth = 1; depth < search_order.size(); ++depth) {
      const size_t future_comp_idx = search_order[depth];
      initial_allowed_masks[future_comp_idx] &=
          pair_table.compat_mask(root_comp_idx, identity_shift_pos, future_comp_idx);
      if (initial_allowed_masks[future_comp_idx] == 0u && timing != nullptr) {
        ++timing->component_pair_table_prunes;
        ++timing->component_search_zero_future_mask_prunes;
      }
    }
  }

  vector<PbctopoBeamState> beam_states;
  PbctopoBeamState root_state;
  root_state.selected_shift = initial_selected_shift;
  root_state.allowed_masks = initial_allowed_masks;
  root_state.future_remaining_shifts = 0;
  beam_states.push_back(std::move(root_state));
  if (timing != nullptr)
    ++timing->component_search_depth_visits[0];

  for (size_t depth = search_order.empty() ? 0 : 1; depth < search_order.size();
       ++depth) {
    vector<PbctopoBeamState> next_states;
    next_states.reserve(beam_states.size() * min<size_t>(beam_width,
                                                         shifts.size()));
    const size_t actual_comp_idx = search_order[depth];
    for (const auto &state : beam_states) {
      if (timing != nullptr &&
          depth < timing->component_search_depth_visits.size()) {
        ++timing->component_search_depth_visits[depth];
      }
      const vector<PbctopoComponentSearchShiftOption> ordered_shifts =
          build_component_search_shift_order(
              shifts, actual_comp_idx, depth, search_order, state.selected_shift,
              state.allowed_masks, pair_table, state.min_pair_d2, state.clash_atom_a,
              state.clash_atom_b, normalized_warm_offsets_ptr, timing);
      for (const auto &option : ordered_shifts) {
        if (timing != nullptr) {
          ++timing->component_shift_trials;
          ++timing->component_large_mode_beam_expanded;
        }
        PbctopoBeamState next_state = state;
        next_state.selected_shift[actual_comp_idx] = option.shift_idx;
        if (actual_comp_idx < next_state.allowed_masks.size() && option.shift_idx < 32)
          next_state.allowed_masks[actual_comp_idx] = (1u << option.shift_idx);
        for (size_t future_depth = depth + 1; future_depth < search_order.size();
             ++future_depth) {
          const size_t future_comp_idx = search_order[future_depth];
          next_state.allowed_masks[future_comp_idx] &=
              pair_table.compat_mask(actual_comp_idx, option.shift_idx, future_comp_idx);
        }
        next_state.shift_norm2 +=
            shifted_components[actual_comp_idx][option.shift_idx].weighted_shift_norm2;
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
      partial_sort(next_states.begin(),
                   next_states.begin() + static_cast<ptrdiff_t>(beam_width),
                   next_states.end(), better_component_beam_state);
      next_states.resize(beam_width);
    } else {
      sort(next_states.begin(), next_states.end(), better_component_beam_state);
    }
    if (timing != nullptr)
      timing->component_large_mode_beam_kept += next_states.size();
    beam_states = std::move(next_states);
  }

  if (timing != nullptr)
    timing->component_large_mode_beam_complete_assignments += beam_states.size();
  const size_t exact_rescore_count =
      min(kPbctopoBeamFinalRescore, beam_states.size());
  if (beam_states.size() > exact_rescore_count) {
    partial_sort(beam_states.begin(),
                 beam_states.begin() + static_cast<ptrdiff_t>(exact_rescore_count),
                 beam_states.end(), better_component_beam_state);
  } else {
    sort(beam_states.begin(), beam_states.end(), better_component_beam_state);
  }
  if (timing != nullptr)
    timing->component_large_mode_beam_final_rescored += exact_rescore_count;
  for (size_t state_idx = 0; state_idx < exact_rescore_count; ++state_idx) {
    double candidate_shift_norm2 = 0.0;
    if (!build_component_points_from_selected_shifts(
            atoms, components, shifted_components, beam_states[state_idx].selected_shift,
            placed_points, candidate_shift_norm2)) {
      continue;
    }
    PbctopoComponentCandidate candidate = score_component_assignment(
        atoms, placed_points, components, beam_states[state_idx].min_pair_d2,
        beam_states[state_idx].clash_atom_a, beam_states[state_idx].clash_atom_b,
        continuity_pairs, anchor_indices, cutoffs.continuity_cutoff,
        candidate_shift_norm2, best.valid ? &best : nullptr, timing);
    if (better_component_candidate(candidate, best)) {
      best = candidate;
      best_points_out = placed_points;
    }
  }

  return best;
}

PbctopoComponentCandidate search_component_shift_candidate(
    const vector<PbctopoAtom> &atoms, box Box,
    const vector<vector<size_t>> &components, const vector<PbctopoShift> &shifts,
    const vector<PbctopoContinuityPair> &continuity_pairs,
    const vector<size_t> &anchor_indices,
    const PbctopoRuntimeCutoffs &cutoffs, vector<Vec3> &best_points_out,
    const vector<Int3> *warm_component_offsets = nullptr,
    const PbctopoCutCandidate *incumbent_cut = nullptr,
    const vector<Vec3> *incumbent_points = nullptr,
    PbctopoTimingStats *timing = nullptr) {
  (void)Box;
  PbctopoComponentCandidate best;
  if (components.empty())
    return best;

  vector<Vec3> placed_points(atoms.size());
  best_points_out.clear();
  best_points_out.reserve(atoms.size());
  vector<size_t> selected_shift(components.size(), shifts.size());
  const auto shifted_components =
      precompute_shifted_components(atoms, components, shifts,
                                    cutoffs.clash_cutoff,
                                    cutoffs.pair_contact_cutoff);

  size_t identity_shift_pos = 0;
  for (size_t shift_idx = 0; shift_idx < shifts.size(); ++shift_idx) {
    if (shifts[shift_idx].ix == 0 && shifts[shift_idx].iy == 0 &&
        shifts[shift_idx].iz == 0) {
      identity_shift_pos = shift_idx;
      break;
    }
  }

  for (const auto &component : components) {
    for (size_t atom_idx : component) {
      placed_points[atom_idx] = {atoms[atom_idx].x, atoms[atom_idx].y,
                                 atoms[atom_idx].z};
    }
  }

  if (incumbent_cut != nullptr && incumbent_cut->valid && !incumbent_cut->ultra_clash &&
      incumbent_points != nullptr && incumbent_points->size() == atoms.size()) {
    best = build_component_candidate_seed(*incumbent_cut, components, *incumbent_points);
    if (timing != nullptr)
      ++timing->component_search_seed_graph;
  } else if (timing != nullptr) {
    ++timing->component_search_seed_none;
  }

  const PbctopoComponentPairShiftTable pair_table =
      build_component_pair_shift_table(components, shifted_components,
                                       cutoffs.clash_cutoff, timing);
  const vector<size_t> atom_component_index =
      build_atom_component_index(atoms.size(), components);
  const vector<size_t> search_order = build_component_search_order(
      components, anchor_indices, atom_component_index, pair_table, identity_shift_pos);
  const vector<Int3> normalized_warm_offsets =
      (warm_component_offsets != nullptr &&
       warm_component_offsets->size() == components.size())
          ? normalize_search_warm_offsets(*warm_component_offsets, components)
          : vector<Int3>{};
  const vector<Int3> *normalized_warm_offsets_ptr =
      normalized_warm_offsets.empty() ? nullptr : &normalized_warm_offsets;

  const uint32_t full_shift_mask = component_search_full_shift_mask(shifts.size());
  vector<uint32_t> initial_allowed_masks(components.size(), full_shift_mask);
  if (!search_order.empty() && identity_shift_pos < 32) {
    const size_t root_comp_idx = search_order.front();
    initial_allowed_masks[root_comp_idx] = (1u << identity_shift_pos);
    selected_shift[root_comp_idx] = identity_shift_pos;
    for (size_t depth = 1; depth < search_order.size(); ++depth) {
      const size_t future_comp_idx = search_order[depth];
      initial_allowed_masks[future_comp_idx] &=
          pair_table.compat_mask(root_comp_idx, identity_shift_pos, future_comp_idx);
      if (initial_allowed_masks[future_comp_idx] == 0u && timing != nullptr) {
        ++timing->component_pair_table_prunes;
        ++timing->component_search_zero_future_mask_prunes;
      }
    }
  }

  function<void(size_t, double, double, size_t, size_t, vector<uint32_t> &)> dfs =
      [&](size_t depth, double shift_norm2, double min_pair_d2,
          size_t clash_atom_a, size_t clash_atom_b,
          vector<uint32_t> &allowed_masks) {
        if (timing != nullptr &&
            depth < timing->component_search_depth_visits.size()) {
          ++timing->component_search_depth_visits[depth];
        }
        if (depth >= search_order.size()) {
          PbctopoComponentCandidate candidate = score_component_assignment(
              atoms, placed_points, components, min_pair_d2, clash_atom_a,
              clash_atom_b, continuity_pairs, anchor_indices,
              cutoffs.continuity_cutoff, shift_norm2,
              best.valid ? &best : nullptr, timing);
          if (better_component_candidate(candidate, best)) {
            best = candidate;
            best_points_out = placed_points;
          }
          return;
        }

        const size_t actual_comp_idx = search_order[depth];
        const vector<PbctopoComponentSearchShiftOption> ordered_shifts =
            build_component_search_shift_order(
                shifts, actual_comp_idx, depth, search_order, selected_shift,
                allowed_masks, pair_table, min_pair_d2, clash_atom_a, clash_atom_b,
                normalized_warm_offsets_ptr, timing);
        for (const auto &option : ordered_shifts) {
          const size_t shift_idx = option.shift_idx;
          if (timing != nullptr)
            ++timing->component_shift_trials;
          const auto &component = components[actual_comp_idx];
          const auto &shifted = shifted_components[actual_comp_idx][shift_idx];
          for (size_t local_idx = 0; local_idx < component.size(); ++local_idx) {
            placed_points[component[local_idx]] = shifted.points[local_idx];
          }
          selected_shift[actual_comp_idx] = shift_idx;
          vector<pair<size_t, uint32_t>> allowed_mask_restore;
          allowed_mask_restore.reserve(search_order.size() - depth);
          if (actual_comp_idx < allowed_masks.size() && shift_idx < 32) {
            allowed_mask_restore.push_back(
                {actual_comp_idx, allowed_masks[actual_comp_idx]});
            allowed_masks[actual_comp_idx] = (1u << shift_idx);
          }
          for (size_t future_depth = depth + 1; future_depth < search_order.size();
               ++future_depth) {
            const size_t future_comp_idx = search_order[future_depth];
            const uint32_t previous_mask = allowed_masks[future_comp_idx];
            const uint32_t updated_mask =
                previous_mask &
                pair_table.compat_mask(actual_comp_idx, shift_idx, future_comp_idx);
            if (updated_mask != previous_mask) {
              allowed_mask_restore.push_back({future_comp_idx, previous_mask});
              allowed_masks[future_comp_idx] = updated_mask;
            }
          }
          dfs(depth + 1, shift_norm2 + shifted.weighted_shift_norm2,
              option.next_min_pair_d2, option.clash_atom_a,
              option.clash_atom_b, allowed_masks);
          for (auto restore_it = allowed_mask_restore.rbegin();
               restore_it != allowed_mask_restore.rend(); ++restore_it) {
            allowed_masks[restore_it->first] = restore_it->second;
          }
          selected_shift[actual_comp_idx] = shifts.size();
        }
      };

  dfs(search_order.empty() ? 0 : 1, 0.0, numeric_limits<double>::infinity(),
      numeric_limits<size_t>::max(), numeric_limits<size_t>::max(),
      initial_allowed_masks);
  return best;
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
    vector<Int3> *best_raw_offsets_out = nullptr) {
  best_cut = {};
  shift_norm2 = 0.0;
  placed_points.clear();
  if (merge_plan.merged_components.empty())
    return false;

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
        cutoffs.continuity_cutoff, refine_points, anchor_backend, timing,
        warm_state, warm_reference_candidate);
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
          atoms, Box, merge_plan.merged_components, shifts, continuity_pairs,
          anchor_indices, cutoffs, merged_warm_offsets_ptr, timing);
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
        selected_seed->offsets, raw_components, merge_plan));
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
        atoms, Box, merge_plan.merged_components, shifts, continuity_pairs,
        anchor_indices, cutoffs, merged_search_points, merged_warm_offsets_ptr,
        search_incumbent, search_incumbent_points, timing);
    if (merged_candidate.valid && !merged_candidate.ultra_clash &&
        merged_search_points.size() == atoms.size()) {
      const vector<Int3> merged_search_offsets = extract_component_offsets_from_points(
          atoms, merge_plan.merged_components, Box, merged_search_points);
      if (merged_search_offsets.size() == merge_plan.merged_components.size()) {
        consider_raw_offsets(expand_merged_offsets_to_raw_components(
            merged_search_offsets, raw_components, merge_plan));
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

vector<double> build_virtual_box_event_planes(const vector<PbctopoAtom> &atoms,
                                              char axis) {
  vector<double> planes;
  planes.reserve(atoms.size() * 2 + 2);
  planes.push_back(0.5);
  planes.push_back(2.5);
  for (const auto &atom : atoms) {
    const double coord =
        (axis == 'x') ? atom.sx : ((axis == 'y') ? atom.sy : atom.sz);
    planes.push_back(coord + 0.5);
    planes.push_back(coord + 1.5);
  }
  sort(planes.begin(), planes.end());
  vector<double> unique_planes;
  unique_planes.reserve(planes.size());
  for (double value : planes) {
    if (value + 1.0e-9 < 0.5 || value - 1.0e-9 > 2.5)
      continue;
    if (!unique_planes.empty() &&
        fabs(unique_planes.back() - value) <= 1.0e-9) {
      continue;
    }
    unique_planes.push_back(value);
  }
  if (unique_planes.size() < 2) {
    unique_planes = {0.5, 2.5};
  } else {
    if (fabs(unique_planes.front() - 0.5) > 1.0e-9)
      unique_planes.insert(unique_planes.begin(), 0.5);
    if (fabs(unique_planes.back() - 2.5) > 1.0e-9)
      unique_planes.push_back(2.5);
  }
  return unique_planes;
}

size_t virtual_box_interval_index(const vector<double> &planes, double center) {
  if (planes.size() < 2)
    return 0;
  const double clamped = clamp_virtual_box_center(center);
  auto it = upper_bound(planes.begin(), planes.end(), clamped + 1.0e-9);
  if (it == planes.begin())
    return 0;
  const size_t idx = static_cast<size_t>(distance(planes.begin(), it) - 1);
  return min(idx, planes.size() - 2);
}

uint64_t virtual_box_state_key(const vector<double> &x_planes,
                               const vector<double> &y_planes,
                               const vector<double> &z_planes, double cut_x,
                               double cut_y, double cut_z) {
  const uint64_t ix = static_cast<uint64_t>(
      virtual_box_interval_index(x_planes, cut_x));
  const uint64_t iy = static_cast<uint64_t>(
      virtual_box_interval_index(y_planes, cut_y));
  const uint64_t iz = static_cast<uint64_t>(
      virtual_box_interval_index(z_planes, cut_z));
  return (ix << 42) ^ (iy << 21) ^ iz;
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
    if (fabs(lhs.score - rhs.score) > 1.0e-9)
      return lhs.score > rhs.score;
    if (fabs(lhs.cx - rhs.cx) > 1.0e-9)
      return lhs.cx < rhs.cx;
    return lhs.cy < rhs.cy;
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
         if (fabs(lhs.gap - rhs.gap) > 1.0e-9)
           return lhs.gap > rhs.gap;
         return lhs.cut < rhs.cut;
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
    if (fabs(lhs.score - rhs.score) > 1.0e-9)
      return lhs.score > rhs.score;
    if (fabs(lhs.cx - rhs.cx) > 1.0e-9)
      return lhs.cx < rhs.cx;
    if (fabs(lhs.cy - rhs.cy) > 1.0e-9)
      return lhs.cy < rhs.cy;
    return lhs.cz < rhs.cz;
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
  if (fabs(lhs.score - rhs.score) > 1.0e-12)
    return lhs.score > rhs.score;
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
    bool skip_clash, double frame_time_ps, const char *reject_context) {
  PbctopoCutCandidate candidate;
  candidate.cut_x = clamp_virtual_box_center(cut_x);
  candidate.cut_y = clamp_virtual_box_center(cut_y);
  candidate.cut_z = clamp_virtual_box_center(cut_z);
  candidate.gap_score = gap_score;
  candidate.broken_edges = 0;
  candidate.edge_stretch = 0.0;
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
    candidate.shift_norm2 +=
        static_cast<double>(actual_shift_x * actual_shift_x +
                            actual_shift_y * actual_shift_y +
                            actual_shift_z * actual_shift_z);
  }

  if (stats.count != atoms.size())
    return PbctopoCutCandidate{};

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
    if (!anchor_graph_metrics_ok(placed_points, *anchor_indices, anchor_mst2,
                                 anchor_max_d2, metric_incumbent)) {
      return PbctopoCutCandidate{};
    }
    candidate.anchor_mst2 = anchor_mst2;
    candidate.anchor_max_d2 = anchor_max_d2;
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
                                     clash_cutoff);
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
          if (fabs(lhs.gap - rhs.gap) > 1.0e-9)
           return lhs.gap > rhs.gap;
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

bool better_stage1_cut_candidate(const PbctopoCutCandidate &lhs,
                                 const PbctopoCutCandidate &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  if (lhs.broken_edges != rhs.broken_edges)
    return lhs.broken_edges < rhs.broken_edges;
  const bool lhs_has_cont =
      isfinite(lhs.continuity_max_d2) && isfinite(lhs.continuity_path);
  const bool rhs_has_cont =
      isfinite(rhs.continuity_max_d2) && isfinite(rhs.continuity_path);
  if (lhs_has_cont != rhs_has_cont)
    return lhs_has_cont;
  if (lhs_has_cont &&
      fabs(lhs.continuity_max_d2 - rhs.continuity_max_d2) > 1.0e-6)
    return lhs.continuity_max_d2 < rhs.continuity_max_d2;
  if (lhs_has_cont &&
      fabs(lhs.continuity_path - rhs.continuity_path) > 1.0e-6)
    return lhs.continuity_path < rhs.continuity_path;
  const bool lhs_has_anchor =
      isfinite(lhs.anchor_max_d2) && isfinite(lhs.anchor_mst2);
  const bool rhs_has_anchor =
      isfinite(rhs.anchor_max_d2) && isfinite(rhs.anchor_mst2);
  if (lhs_has_anchor != rhs_has_anchor)
    return lhs_has_anchor;
  if (lhs_has_anchor && fabs(lhs.anchor_max_d2 - rhs.anchor_max_d2) > 1.0e-6)
    return lhs.anchor_max_d2 < rhs.anchor_max_d2;
  if (lhs_has_anchor && fabs(lhs.anchor_mst2 - rhs.anchor_mst2) > 1.0e-6)
    return lhs.anchor_mst2 < rhs.anchor_mst2;
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
  const int kappa_cmp = compare_shape_kappa(lhs.shape_kappa2, rhs.shape_kappa2);
  if (kappa_cmp != 0)
    return kappa_cmp < 0;
  const int asph_cmp =
      compare_shape_moment(lhs.shape_asphericity, rhs.shape_asphericity);
  if (asph_cmp != 0)
    return asph_cmp < 0;
  const int acyl_cmp =
      compare_shape_moment(lhs.shape_acylindricity, rhs.shape_acylindricity);
  if (acyl_cmp != 0)
    return acyl_cmp < 0;
  if (fabs(lhs.edge_stretch - rhs.edge_stretch) > 1.0e-6)
    return lhs.edge_stretch < rhs.edge_stretch;
  if (fabs(lhs.gap_score - rhs.gap_score) > 1.0e-9)
    return lhs.gap_score > rhs.gap_score;
  return false;
}

bool better_preclash_cut_candidate(const PbctopoCutCandidate &lhs,
                                   const PbctopoCutCandidate &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  if (lhs.broken_edges != rhs.broken_edges)
    return lhs.broken_edges < rhs.broken_edges;
  const bool lhs_has_cont =
      isfinite(lhs.continuity_max_d2) && isfinite(lhs.continuity_path);
  const bool rhs_has_cont =
      isfinite(rhs.continuity_max_d2) && isfinite(rhs.continuity_path);
  if (lhs_has_cont != rhs_has_cont)
    return lhs_has_cont;
  if (lhs_has_cont &&
      fabs(lhs.continuity_max_d2 - rhs.continuity_max_d2) > 1.0e-6)
    return lhs.continuity_max_d2 < rhs.continuity_max_d2;
  if (lhs_has_cont &&
      fabs(lhs.continuity_path - rhs.continuity_path) > 1.0e-6)
    return lhs.continuity_path < rhs.continuity_path;
  const bool lhs_has_anchor =
      isfinite(lhs.anchor_max_d2) && isfinite(lhs.anchor_mst2);
  const bool rhs_has_anchor =
      isfinite(rhs.anchor_max_d2) && isfinite(rhs.anchor_mst2);
  if (lhs_has_anchor != rhs_has_anchor)
    return lhs_has_anchor;
  if (lhs_has_anchor && fabs(lhs.anchor_max_d2 - rhs.anchor_max_d2) > 1.0e-6)
    return lhs.anchor_max_d2 < rhs.anchor_max_d2;
  if (lhs_has_anchor && fabs(lhs.anchor_mst2 - rhs.anchor_mst2) > 1.0e-6)
    return lhs.anchor_mst2 < rhs.anchor_mst2;
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
  const int kappa_cmp = compare_shape_kappa(lhs.shape_kappa2, rhs.shape_kappa2);
  if (kappa_cmp != 0)
    return kappa_cmp < 0;
  const int asph_cmp =
      compare_shape_moment(lhs.shape_asphericity, rhs.shape_asphericity);
  if (asph_cmp != 0)
    return asph_cmp < 0;
  const int acyl_cmp =
      compare_shape_moment(lhs.shape_acylindricity, rhs.shape_acylindricity);
  if (acyl_cmp != 0)
    return acyl_cmp < 0;
  if (fabs(lhs.edge_stretch - rhs.edge_stretch) > 1.0e-6)
    return lhs.edge_stretch < rhs.edge_stretch;
  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-6)
    return lhs.shift_norm2 < rhs.shift_norm2;
  if (fabs(lhs.gap_score - rhs.gap_score) > 1.0e-9)
    return lhs.gap_score > rhs.gap_score;
  return false;
}

bool better_cut_candidate(const PbctopoCutCandidate &lhs,
                          const PbctopoCutCandidate &rhs) {
  if (!lhs.valid)
    return false;
  if (!rhs.valid)
    return true;
  if (lhs.ultra_clash != rhs.ultra_clash)
    return rhs.ultra_clash;
  if (lhs.broken_edges != rhs.broken_edges)
    return lhs.broken_edges < rhs.broken_edges;
  const bool lhs_has_cont =
      isfinite(lhs.continuity_max_d2) && isfinite(lhs.continuity_path);
  const bool rhs_has_cont =
      isfinite(rhs.continuity_max_d2) && isfinite(rhs.continuity_path);
  if (lhs_has_cont != rhs_has_cont)
    return lhs_has_cont;
  if (lhs_has_cont &&
      fabs(lhs.continuity_max_d2 - rhs.continuity_max_d2) > 1.0e-6)
    return lhs.continuity_max_d2 < rhs.continuity_max_d2;
  if (lhs_has_cont &&
      fabs(lhs.continuity_path - rhs.continuity_path) > 1.0e-6)
    return lhs.continuity_path < rhs.continuity_path;
  const bool lhs_has_anchor =
      isfinite(lhs.anchor_max_d2) && isfinite(lhs.anchor_mst2);
  const bool rhs_has_anchor =
      isfinite(rhs.anchor_max_d2) && isfinite(rhs.anchor_mst2);
  if (lhs_has_anchor != rhs_has_anchor)
    return lhs_has_anchor;
  if (lhs_has_anchor && fabs(lhs.anchor_max_d2 - rhs.anchor_max_d2) > 1.0e-6)
    return lhs.anchor_max_d2 < rhs.anchor_max_d2;
  if (lhs_has_anchor && fabs(lhs.anchor_mst2 - rhs.anchor_mst2) > 1.0e-6)
    return lhs.anchor_mst2 < rhs.anchor_mst2;
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
  const int kappa_cmp = compare_shape_kappa(lhs.shape_kappa2, rhs.shape_kappa2);
  if (kappa_cmp != 0)
    return kappa_cmp < 0;
  const int asph_cmp =
      compare_shape_moment(lhs.shape_asphericity, rhs.shape_asphericity);
  if (asph_cmp != 0)
    return asph_cmp < 0;
  const int acyl_cmp =
      compare_shape_moment(lhs.shape_acylindricity, rhs.shape_acylindricity);
  if (acyl_cmp != 0)
    return acyl_cmp < 0;
  if (fabs(lhs.edge_stretch - rhs.edge_stretch) > 1.0e-6)
    return lhs.edge_stretch < rhs.edge_stretch;
  if (fabs(lhs.min_pair_d2 - rhs.min_pair_d2) > 1.0e-9)
    return lhs.min_pair_d2 > rhs.min_pair_d2;
  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-6)
    return lhs.shift_norm2 < rhs.shift_norm2;
  if (fabs(lhs.cut_x - rhs.cut_x) > 1.0e-9)
    return lhs.cut_x < rhs.cut_x;
  if (fabs(lhs.cut_y - rhs.cut_y) > 1.0e-9)
    return lhs.cut_y < rhs.cut_y;
  return lhs.cut_z < rhs.cut_z;
}

bool preclash_cut_candidate_cannot_beat_incumbent(
    const PbctopoCutCandidate &candidate,
    const PbctopoCutCandidate *incumbent) {
  if (!candidate.valid || incumbent == nullptr || !incumbent->valid ||
      incumbent->ultra_clash)
    return false;
  return better_preclash_cut_candidate(*incumbent, candidate);
}

bool cut_candidate_cannot_beat_incumbent(const PbctopoCutCandidate &partial,
                                         const PbctopoCutCandidate *incumbent) {
  if (incumbent == nullptr || !incumbent->valid || incumbent->ultra_clash ||
      partial.ultra_clash) {
    return false;
  }
  if (partial.broken_edges > incumbent->broken_edges)
    return true;
  if (isfinite(partial.continuity_max_d2) &&
      isfinite(incumbent->continuity_max_d2) &&
      partial.continuity_max_d2 > incumbent->continuity_max_d2 + 1.0e-6 &&
      partial.continuity_path > incumbent->continuity_path + 1.0e-6) {
    return true;
  }
  if (isfinite(partial.anchor_max_d2) && isfinite(incumbent->anchor_max_d2) &&
      partial.anchor_max_d2 > incumbent->anchor_max_d2 + 1.0e-6 &&
      partial.anchor_mst2 > incumbent->anchor_mst2 + 1.0e-6) {
    return true;
  }
  if (principal_rg2_significantly_worse(partial.principal_rg2_max,
                                        incumbent->principal_rg2_max) &&
      principal_rg2_significantly_worse(partial.principal_rg2_mid,
                                        incumbent->principal_rg2_mid)) {
    return true;
  }
  if (shape_kappa_significantly_worse(partial.shape_kappa2,
                                      incumbent->shape_kappa2) &&
      shape_moment_significantly_worse(partial.shape_asphericity,
                                       incumbent->shape_asphericity)) {
    return true;
  }
  if (partial.broken_edges == incumbent->broken_edges &&
      partial.edge_stretch > incumbent->edge_stretch + 1.0e-6) {
    return true;
  }
  return false;
}

void analyze_cut_candidate_ultraclash(const vector<Vec3> &placed_points,
                                      const unordered_set<uint64_t> &edge_set,
                                      PbctopoCutCandidate &candidate,
                                      double cutoff) {
  thread_local PbctopoUltraclashScratch scratch;
  candidate.min_pair_d2 = numeric_limits<double>::infinity();
  candidate.clash_atom_a = numeric_limits<size_t>::max();
  candidate.clash_atom_b = numeric_limits<size_t>::max();
  candidate.ultra_clash = false;
  if (!candidate.valid || placed_points.size() < 2)
    return;

  const double inv_grid = 1.0 / cutoff;
  const double cutoff2 = cutoff * cutoff;
  auto accumulate_pair = [&](size_t atom_idx, size_t other_idx) {
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
    if (d2 + kDistanceEps < cutoff2) {
      candidate.ultra_clash = true;
      return true;
    }
    return false;
  };

  auto analyze_with_hash_grid = [&]() {
    unordered_map<uint64_t, vector<size_t>> clash_grid;
    clash_grid.reserve(max<size_t>(64, placed_points.size()));
    for (size_t atom_idx = 0; atom_idx < placed_points.size(); ++atom_idx) {
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
  const BoxMatrix box_matrix(Box);
  double min_sx = numeric_limits<double>::infinity();
  double min_sy = numeric_limits<double>::infinity();
  double min_sz = numeric_limits<double>::infinity();
  double max_sx = -numeric_limits<double>::infinity();
  double max_sy = -numeric_limits<double>::infinity();
  double max_sz = -numeric_limits<double>::infinity();

  for (const auto &point : placed_points) {
    const double sx = box_matrix.inv[0][0] * point.x +
                      box_matrix.inv[0][1] * point.y +
                      box_matrix.inv[0][2] * point.z;
    const double sy = box_matrix.inv[1][0] * point.x +
                      box_matrix.inv[1][1] * point.y +
                      box_matrix.inv[1][2] * point.z;
    const double sz = box_matrix.inv[2][0] * point.x +
                      box_matrix.inv[2][1] * point.y +
                      box_matrix.inv[2][2] * point.z;
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

vector<PbctopoShift> build_topology_shifts(box Box) {
  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);

  vector<PbctopoShift> shifts;
  shifts.reserve(27);
  for (size_t idx = 0; idx < 27; ++idx) {
    const int *coeff = kTopoImageCoeff[idx];
    shifts.push_back({idx,
                      coeff[0],
                      coeff[1],
                      coeff[2],
                      coeff[0] * v1[0] + coeff[1] * v2[0] + coeff[2] * v3[0],
                      coeff[0] * v1[1] + coeff[1] * v2[1] + coeff[2] * v3[1],
                      coeff[0] * v1[2] + coeff[1] * v2[2] + coeff[2] * v3[2]});
  }
  stable_sort(shifts.begin(), shifts.end(),
              [](const PbctopoShift &lhs, const PbctopoShift &rhs) {
                const int lhs_norm2 =
                    lhs.ix * lhs.ix + lhs.iy * lhs.iy + lhs.iz * lhs.iz;
                const int rhs_norm2 =
                    rhs.ix * rhs.ix + rhs.iy * rhs.iy + rhs.iz * rhs.iz;
                if (lhs_norm2 != rhs_norm2)
                  return lhs_norm2 < rhs_norm2;
                return lhs.shift_index < rhs.shift_index;
              });
  return shifts;
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
         if (fabs(lhs.second - rhs.second) > 1.0e-9)
           return lhs.second < rhs.second;
         const double lhs_shift = node_shift_norm2(nodes[lhs.first]);
         const double rhs_shift = node_shift_norm2(nodes[rhs.first]);
         if (fabs(lhs_shift - rhs_shift) > 1.0e-9)
           return lhs_shift < rhs_shift;
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
  if (fabs(lhs.fragment_mst2 - rhs.fragment_mst2) > 1.0e-6)
    return lhs.fragment_mst2 < rhs.fragment_mst2;
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
  if (fabs(lhs.bbox_volume - rhs.bbox_volume) > 1.0e-6)
    return lhs.bbox_volume < rhs.bbox_volume;
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
  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-6)
    return lhs.shift_norm2 < rhs.shift_norm2;
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
  if (fabs(lhs.fragment_mst2 - rhs.fragment_mst2) > 1.0e-6)
    return lhs.fragment_mst2 < rhs.fragment_mst2;
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
  if (fabs(lhs.bbox_volume - rhs.bbox_volume) > 1.0e-6)
    return lhs.bbox_volume < rhs.bbox_volume;
  if (lhs.hull_valid && rhs.hull_valid &&
      fabs(lhs.hull_volume - rhs.hull_volume) > 1.0e-6)
    return lhs.hull_volume < rhs.hull_volume;
  if (fabs(lhs.min_interfragment_d2 - rhs.min_interfragment_d2) > 1.0e-6)
    return lhs.min_interfragment_d2 < rhs.min_interfragment_d2;
  if (fabs(lhs.shift_norm2 - rhs.shift_norm2) > 1.0e-6)
    return lhs.shift_norm2 < rhs.shift_norm2;
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
                                 TitanExecutor *executor) {
  using clock = chrono::steady_clock;
  const size_t frame_count =
      mol2.m_mol.empty() ? this->m_mol.size()
                         : min(this->m_mol.size(), mol2.m_mol.size());
  if (frame_count == 0)
    return 0;

  reset_pbctopo_log_file();
  size_t fallback_frames = 0;
  size_t partial_components = 0;
  size_t valid_candidates = 0;
  size_t hull_fallback_candidates = 0;
  PbctopoTimingStats timing;
  PbctopoWarmState warm_state;
  const PbctopoRuntimeCutoffs runtime_cutoffs =
      resolve_pbctopo_cutoffs(Arguments);
  const PbctopoAnchorBackend anchor_backend =
      resolve_pbctopo_anchor_backend(Arguments);
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
                          to_string(frame_count) + ", anchor_backend=" +
                          pbctopo_anchor_backend_name(anchor_backend) +
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
  auto &seen_virtual_states = scratch.seen_virtual_states;
  auto &warm_points = scratch.warm_points;
  auto &coarse_warm_points = scratch.coarse_warm_points;
  auto &placed_points = scratch.placed_points;
  auto &best_points = scratch.best_points;
  auto &best_safe_points = scratch.best_safe_points;
  auto &refine_points = scratch.refine_points;
  auto &verify_points = scratch.verify_points;
  auto &trial_points = scratch.trial_points;
  auto &component_points = scratch.component_points;
  auto &coarse_candidates = scratch.coarse_candidates;
  auto &coarse_seed_options = scratch.coarse_seed_options;
  auto &coarse_results = scratch.coarse_results;
  auto &x_planes = scratch.x_planes;
  auto &y_planes = scratch.y_planes;
  auto &z_planes = scratch.z_planes;
  auto &projection_index = scratch.projection_index;
  vector<PbctopoHullSpan> hull_owner_spans;
  hull_owner_spans.reserve(8);
  vector<Vec3> hull_rescue_points;
  hull_rescue_points.reserve(trajectory_context.total_atoms_hint);
  bool projection_clash_skip_active = false;
  auto fallback_pbci_frame = [&](size_t frame_index) {
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
    vector<pair<int, int>> fallback_contacting_owner_pairs =
        build_contacting_owner_pairs({}, fallback_atoms,
                                     this->m_mol[frame_index].Box,
                                     runtime_cutoffs.pair_contact_cutoff);
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
    ++timing.pbci_fallback_frames;
  };

  for (size_t frame_idx = 0; frame_idx < frame_count; ++frame_idx) {
    ++timing.frames;
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
    timing.append_s += chrono::duration<double>(clock::now() - stage_started).count();

    if (base_atoms.empty())
      continue;

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
    build_local_neighborhood_graph(base_atoms, this->m_mol[frame_idx].Box,
                                   local_edges, local_edge_set,
                                   runtime_cutoffs.bond_cutoff);
    timing.graph_s += chrono::duration<double>(clock::now() - stage_started).count();

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
    const vector<pair<int, int>> local_contacting_owner_pairs =
        build_contacting_owner_pairs(local_edges, base_atoms,
                                     this->m_mol[frame_idx].Box,
                                     runtime_cutoffs.pair_contact_cutoff);
    const vector<pair<int, int>> projection_contacting_owner_pairs =
        build_contacting_owner_pairs(projection_edges, projection_atoms,
                                     this->m_mol[frame_idx].Box,
                                     runtime_cutoffs.pair_contact_cutoff);
    timing.projection_edge_s +=
        chrono::duration<double>(clock::now() - stage_started).count();

    stage_started = clock::now();
    components = build_connected_components_from_local_edges(base_atoms.size(), local_edges);
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
      if (!can_use_component_logic || component_context_ready)
        return;
      const auto component_context_started = clock::now();
      component_atoms =
          unwrap_atoms_by_local_graph(base_atoms, local_edges,
                                      this->m_mol[frame_idx].Box);
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
                  local_contacting_owner_pairs, split_owner_pairs, split_gap,
                  &detail)) {
            return false;
          }

          const vector<Int3> rescue_offsets = extract_component_offsets_from_points(
              base_atoms, components, this->m_mol[frame_idx].Box,
              candidate_points);
          double rescue_shift_norm2 = 0.0;
          for (size_t comp_idx = 0; comp_idx < rescue_offsets.size() &&
                                   comp_idx < components.size();
               ++comp_idx) {
            const auto &offset = rescue_offsets[comp_idx];
            const double comp_shift_norm2 = static_cast<double>(
                offset.x * offset.x + offset.y * offset.y + offset.z * offset.z);
            rescue_shift_norm2 +=
                comp_shift_norm2 * static_cast<double>(components[comp_idx].size());
          }

          const BoxMatrix rescue_box(this->m_mol[frame_idx].Box);
          double min_sx = numeric_limits<double>::infinity();
          double min_sy = numeric_limits<double>::infinity();
          double min_sz = numeric_limits<double>::infinity();
          double max_sx = -numeric_limits<double>::infinity();
          double max_sy = -numeric_limits<double>::infinity();
          double max_sz = -numeric_limits<double>::infinity();
          for (const auto &point : candidate_points) {
            const double sx = rescue_box.inv[0][0] * point.x +
                              rescue_box.inv[0][1] * point.y +
                              rescue_box.inv[0][2] * point.z;
            const double sy = rescue_box.inv[1][0] * point.x +
                              rescue_box.inv[1][1] * point.y +
                              rescue_box.inv[1][2] * point.z;
            const double sz = rescue_box.inv[2][0] * point.x +
                              rescue_box.inv[2][1] * point.y +
                              rescue_box.inv[2][2] * point.z;
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
          const bool rescue_ok =
              candidate_out.valid &&
              continuity_metrics_ok(candidate_points,
                                    trajectory_context.continuity_pairs,
                                    runtime_cutoffs.continuity_cutoff,
                                    continuity_path, continuity_max_d2,
                                    incumbent_cut) &&
              anchor_graph_metrics_ok(
                  candidate_points, trajectory_context.continuity_anchor_indices,
                  anchor_mst2, anchor_max_d2, incumbent_cut);
          if (!rescue_ok)
            return false;

          candidate_out.continuity_path = continuity_path;
          candidate_out.continuity_max_d2 = continuity_max_d2;
          candidate_out.anchor_mst2 = anchor_mst2;
          candidate_out.anchor_max_d2 = anchor_max_d2;
          analyze_cut_candidate_ultraclash(candidate_points, local_edge_set,
                                           candidate_out,
                                           runtime_cutoffs.clash_cutoff);
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
          &projection_anchor_indices, &projection_contacting_owner_pairs,
          runtime_cutoffs.continuity_cutoff, nullptr,
          skip_projection_candidate_clash, numeric_limits<double>::quiet_NaN(),
          nullptr);
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
            &local_contacting_owner_pairs, runtime_cutoffs.continuity_cutoff,
            nullptr, false, this->m_mol[frame_idx].time_of_snap,
            "warm_projection");
        if (early_projection.valid && !early_projection.ultra_clash &&
          candidate_close_enough_to_skip_search(
              early_projection, warm_state.previous_projection_candidate)) {
          ++timing.early_projection_used;
          warm_state.previous_projection_candidate = early_projection;
          warm_state.has_projection_candidate = true;
          if (better_cut_candidate(early_projection, prefetched_safe_candidate)) {
            prefetched_safe_candidate = early_projection;
            prefetched_safe_points = warm_points;
          }
          timing.early_projection_s +=
              chrono::duration<double>(clock::now() - early_projection_started)
                  .count();
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
      rebase_component_offsets(warm_offsets, components);
      PbctopoCutCandidate early_component = score_component_offset_solution(
          component_atoms, components, this->m_mol[frame_idx].Box, warm_offsets,
          local_edges, local_edge_set, runtime_cutoffs.clash_cutoff,
          trajectory_context.continuity_pairs,
          trajectory_context.continuity_anchor_indices,
          runtime_cutoffs.continuity_cutoff, warm_points,
          &warm_state.previous_projection_candidate, nullptr, true,
          PbctopoRefineScorePhase::Exact, nullptr, nullptr, nullptr,
          anchor_backend);
      PbctopoOwnerPairSplitDetail early_component_split_detail;
      size_t early_component_split_owner_pairs = 0;
      double early_component_split_gap = 0.0;
      const bool early_component_split =
          early_component.valid && !early_component.ultra_clash &&
          candidate_has_split_contacting_owners(
              warm_points, base_atoms, this->m_mol[frame_idx].Box,
              local_contacting_owner_pairs, early_component_split_owner_pairs,
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
        }
        timing.warm_component_s +=
            chrono::duration<double>(clock::now() - warm_component_started)
                .count();
      }
      timing.warm_component_s +=
          chrono::duration<double>(clock::now() - warm_component_started)
              .count();
    }

    stage_started = clock::now();
    collect_projection_volume_seeds(projection_atoms, 16, kProjectionTopSeeds, 0.22,
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
      volume_seeds.push_back(previous_seed);
      ++timing.previous_seed_used;
    }
    timing.seed_s += chrono::duration<double>(clock::now() - stage_started).count();
    const double coarse_step = 2.0 / 15.0;
    stage_started = clock::now();
    x_planes = build_virtual_box_event_planes(projection_atoms, 'x');
    y_planes = build_virtual_box_event_planes(projection_atoms, 'y');
    z_planes = build_virtual_box_event_planes(projection_atoms, 'z');
    timing.event_planes_s +=
        chrono::duration<double>(clock::now() - stage_started).count();
    seen_virtual_states.reserve(512);

    PbctopoCutCandidate best_candidate;
    PbctopoCutCandidate best_safe_candidate;
    size_t safe_candidates = 0;
    auto accept_rescue_candidate =
        [&](const vector<Vec3> &candidate_points,
            const PbctopoCutCandidate &candidate) {
          best_safe_candidate = candidate;
          best_candidate = candidate;
          best_safe_points = candidate_points;
          best_points = candidate_points;
          safe_candidates = max<size_t>(safe_candidates, 1);
          projection_clash_skip_active =
              candidate_has_projection_clash_margin(best_safe_candidate,
                                                    runtime_cutoffs.clash_cutoff);
          warm_state.previous_projection_candidate = best_safe_candidate;
          warm_state.has_projection_candidate = true;
          apply_cut_candidate_timed(best_safe_points);
        };
    placed_points.reserve(projection_atoms.size());
    best_points.reserve(base_atoms.size());
    best_safe_points.reserve(base_atoms.size());
    coarse_candidates.reserve(12);
    coarse_seed_options.reserve(volume_seeds.size());
    if (prefetched_safe_candidate.valid && !prefetched_safe_candidate.ultra_clash &&
        prefetched_safe_points.size() == base_atoms.size()) {
      best_candidate = prefetched_safe_candidate;
      best_safe_candidate = prefetched_safe_candidate;
      best_points = prefetched_safe_points;
      best_safe_points = prefetched_safe_points;
      safe_candidates = 1;
    }

    stage_started = clock::now();
    for (const auto &seed : volume_seeds) {
      const uint64_t state_key = virtual_box_state_key(x_planes, y_planes, z_planes,
                                                       seed.cx, seed.cy, seed.cz);
      if (!seen_virtual_states.insert(state_key).second)
        continue;
      coarse_seed_options.push_back({state_key, seed});
    }
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
                  &projection_contacting_owner_pairs,
                  runtime_cutoffs.continuity_cutoff, nullptr,
                  skip_projection_candidate_clash,
                  numeric_limits<double>::quiet_NaN(), nullptr);
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
                &projection_contacting_owner_pairs,
                runtime_cutoffs.continuity_cutoff, nullptr,
                skip_projection_candidate_clash,
                numeric_limits<double>::quiet_NaN(), nullptr);
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
            &projection_contacting_owner_pairs, runtime_cutoffs.continuity_cutoff,
            nullptr, skip_projection_candidate_clash,
            numeric_limits<double>::quiet_NaN(), nullptr);
        if (skip_projection_candidate_clash)
          coarse_projection_skip_evals.fetch_add(1, memory_order_relaxed);
      }
    }
    timing.projection_clash_skip_evals +=
        coarse_projection_skip_evals.load(memory_order_relaxed);
    for (const auto &candidate : coarse_results) {
      push_top_virtual_box_candidate(coarse_candidates, candidate,
                                     kProjectionKeepCandidates);
    }
    timing.coarse_s += chrono::duration<double>(clock::now() - stage_started).count();

    auto evaluate_refine_center = [&](double center_x, double center_y,
                                      double center_z, double gap_score,
                                      vector<Vec3> &scratch_points,
                                      const PbctopoCutCandidate *incumbent) {
      const double clamped_x = clamp_virtual_box_center(center_x);
      const double clamped_y = clamp_virtual_box_center(center_y);
      const double clamped_z = clamp_virtual_box_center(center_z);
      const uint64_t state_key = virtual_box_state_key(x_planes, y_planes, z_planes,
                                                       clamped_x, clamped_y,
                                                       clamped_z);
       if (!seen_virtual_states.insert(state_key).second)
         return PbctopoCutCandidate{};
        PbctopoCutCandidate candidate = evaluate_virtual_box_candidate(
            projection_atoms, projection_edges, projection_edge_set,
            translation_cache, this->m_mol[frame_idx].Box, clamped_x, clamped_y,
            clamped_z, gap_score,
            scratch_points, runtime_cutoffs.clash_cutoff,
            &projection_continuity_pairs, &projection_anchor_indices,
            &projection_contacting_owner_pairs, runtime_cutoffs.continuity_cutoff,
            incumbent, skip_projection_candidate_clash,
            numeric_limits<double>::quiet_NaN(), nullptr);
        if (skip_projection_candidate_clash)
          ++timing.projection_clash_skip_evals;
        ++timing.refine_candidates;
        ++valid_candidates;
        return candidate;
     };

    stage_started = clock::now();
    refine_points.reserve(projection_atoms.size());
    verify_points.reserve(base_atoms.size());
    trial_points.reserve(projection_atoms.size());
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
      double center_x = coarse.cut_x;
      double center_y = coarse.cut_y;
      double center_z = coarse.cut_z;
      PbctopoCutCandidate local_best = evaluate_refine_center(
          center_x, center_y, center_z, coarse.gap_score, refine_points, nullptr);
      double step = 0.5 * coarse_step;
      for (int pass = 0; pass < 4 && step >= (coarse_step / 16.0); ++pass) {
        bool improved = false;
        const array<array<double, 3>, 6> directions = {
            array<double, 3>{step, 0.0, 0.0},   array<double, 3>{-step, 0.0, 0.0},
            array<double, 3>{0.0, step, 0.0},   array<double, 3>{0.0, -step, 0.0},
            array<double, 3>{0.0, 0.0, step},   array<double, 3>{0.0, 0.0, -step}};
        for (const auto &dir : directions) {
          trial_points.clear();
          const PbctopoCutCandidate *local_incumbent =
              (local_best.valid && !local_best.ultra_clash) ? &local_best : nullptr;
          PbctopoCutCandidate trial = evaluate_refine_center(
              center_x + dir[0], center_y + dir[1], center_z + dir[2],
              coarse.gap_score, trial_points, local_incumbent);
          if (better_cut_candidate(trial, local_best)) {
            local_best = trial;
            center_x = trial.cut_x;
            center_y = trial.cut_y;
            center_z = trial.cut_z;
            refine_points = trial_points;
            improved = true;
          }
        }
        if (!improved)
          step *= 0.5;
      }
      if (!local_best.valid)
        continue;
      const PbctopoCutCandidate *verify_incumbent =
          (best_safe_candidate.valid && !best_safe_candidate.ultra_clash)
              ? &best_safe_candidate
              : nullptr;
      PbctopoCutCandidate verified_candidate = evaluate_virtual_box_candidate(
          base_atoms, local_edges, local_edge_set, translation_cache,
          this->m_mol[frame_idx].Box,
          center_x, center_y, center_z, coarse.gap_score, verify_points,
          runtime_cutoffs.clash_cutoff, &trajectory_context.continuity_pairs,
          &trajectory_context.continuity_anchor_indices,
          &local_contacting_owner_pairs, runtime_cutoffs.continuity_cutoff,
          verify_incumbent, false, this->m_mol[frame_idx].time_of_snap,
          "verified_candidate");
      if (!verified_candidate.valid)
        continue;
      if (!verified_candidate.ultra_clash) {
        ++safe_candidates;
        if (better_cut_candidate(verified_candidate, best_safe_candidate)) {
          best_safe_candidate = verified_candidate;
          best_safe_points = verify_points;
        }
      }
      if (better_cut_candidate(verified_candidate, best_candidate)) {
        best_candidate = verified_candidate;
        best_points = verify_points;
      }
    }
    timing.refine_s += chrono::duration<double>(clock::now() - stage_started).count();

    component_scoring_edges = local_edges;
    // Keep component scoring aligned with the original virtual-box local graph.
    // Experimental interface-anchor edges can overweight noisy inter-chain contacts
    // and pull component solutions away from the compact box-consistent candidate.

    bool component_search_used = false;
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
        PbctopoOwnerPairSplitDetail split_detail;
        size_t split_owner_pairs = 0;
        double max_owner_raw_minus_minimg_delta = 0.0;
        if (candidate_has_split_contacting_owners(
                candidate_points, base_atoms, this->m_mol[frame_idx].Box,
                local_contacting_owner_pairs, split_owner_pairs,
                max_owner_raw_minus_minimg_delta, &split_detail)) {
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
          best_safe_points = candidate_points;
          accepted = true;
          better_than_best = true;
        } else if (candidate_close_enough_to_skip_search(candidate,
                                                         best_safe_candidate)) {
          accepted = true;
        }
        if (!accepted)
          return false;
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
            component_as_cut.clash_atom_a = component_candidate.clash_atom_a;
            component_as_cut.clash_atom_b = component_candidate.clash_atom_b;
            component_as_cut.ultra_clash = component_candidate.ultra_clash;
            component_as_cut.valid = component_candidate.valid;
            return accept_component_cut(component_as_cut, candidate_points,
                                        accepted_offsets, source_label);
          };
      if (same_component_warm_layout) {
        vector<Int3> warm_offsets = warm_state.previous_component_offsets;
        rebase_component_offsets(warm_offsets, components);
        PbctopoCutCandidate warm_cut = improve_component_offsets_locally(
            component_atoms, components, this->m_mol[frame_idx].Box, shifts,
            component_scoring_edges, local_edge_set,
            trajectory_context.continuity_pairs,
            trajectory_context.continuity_anchor_indices, warm_offsets,
            runtime_cutoffs.clash_cutoff,
            runtime_cutoffs.continuity_cutoff, component_points, anchor_backend,
            &timing,
            &warm_state, component_warm_reference);
        PbctopoOwnerPairSplitDetail warm_cut_split_detail;
        size_t warm_cut_split_owner_pairs = 0;
        double warm_cut_split_gap = 0.0;
        const bool warm_cut_split =
            warm_cut.valid && !warm_cut.ultra_clash &&
            candidate_has_split_contacting_owners(
                component_points, base_atoms, this->m_mol[frame_idx].Box,
                local_contacting_owner_pairs, warm_cut_split_owner_pairs,
                warm_cut_split_gap, &warm_cut_split_detail);
        if (warm_cut_split && warm_cut_split_detail.valid) {
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
        if (!component_search_used) {
          const auto component_graph_started = clock::now();
          const bool graph_ok = solve_component_offsets_via_pair_graph(
              component_atoms, this->m_mol[frame_idx].Box, components, shifts,
              component_scoring_edges, local_edge_set, component_points,
              trajectory_context.continuity_pairs,
              trajectory_context.continuity_anchor_indices, runtime_cutoffs,
              graph_as_cut, component_shift_norm2, &timing, anchor_backend,
              same_component_warm_layout ? &warm_state.previous_component_offsets
                                         : nullptr,
              &warm_state, component_warm_reference, &graph_tree_offsets,
              &graph_refined_offsets, &graph_tree_edges);
          timing.component_graph_s +=
              chrono::duration<double>(clock::now() - component_graph_started)
                  .count();
          if (!graph_ok) {
            ++timing.component_graph_failed;
            component_graph_fallthrough_reason = "graph_fail";
            clear_component_trial_hints(warm_state);
          }
          if (graph_ok && component_points.size() != base_atoms.size()) {
            ++timing.component_graph_points_mismatch;
            component_graph_fallthrough_reason = "points_mismatch";
          }
          if (graph_ok && component_points.size() == base_atoms.size()) {
            if (accept_component_cut(graph_as_cut, component_points, nullptr,
                                     "component graph solve")) {
              ++timing.component_graph_used;
            } else if (!graph_as_cut.valid) {
              ++timing.component_graph_invalid;
              component_graph_fallthrough_reason = "invalid";
            } else if (graph_as_cut.ultra_clash) {
              ++timing.component_graph_ultraclash;
              component_graph_fallthrough_reason = "ultra_clash";
            } else {
              ++timing.component_graph_not_better_not_close;
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
          PbctopoComponentCandidate component_candidate =
              search_component_shift_candidate(
                  component_atoms, this->m_mol[frame_idx].Box, components, shifts,
                  trajectory_context.continuity_pairs,
                  trajectory_context.continuity_anchor_indices, runtime_cutoffs,
                  component_points,
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
                  &timing);
          timing.component_search_s +=
              chrono::duration<double>(clock::now() - component_search_started)
                  .count();
          if (accept_component_candidate(component_candidate, component_points,
                                         nullptr, "component search winner")) {
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
            build_chain_merged_component_plan(component_atoms, components);
        effective_component_count = merge_plan.merged_components.size();
        if (raw_component_count > effective_component_count)
          timing.component_large_mode_merge_count +=
              raw_component_count - effective_component_count;
        const bool allow_exact_search_on_merged =
            effective_component_count <= exact_component_limit;
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
            &merged_mode_raw_offsets);
        timing.component_graph_s +=
            chrono::duration<double>(clock::now() - merge_mode_started).count();
        if (merge_ok &&
            accept_component_cut(merged_mode_cut, component_points,
                                 &merged_mode_raw_offsets,
                                 "component merge-graph solve")) {
          ++timing.component_graph_used;
        } else {
          ++timing.component_large_mode_fallback_frames;
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
                component_points, beam_width,
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
                &timing);
        timing.component_search_s +=
            chrono::duration<double>(clock::now() - beam_started).count();
        if (!accept_component_candidate(beam_candidate, component_points, nullptr,
                                        "component beam solve")) {
          ++timing.component_large_mode_fallback_frames;
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
      projection_clash_skip_active =
          candidate_has_projection_clash_margin(best_safe_candidate,
                                                runtime_cutoffs.clash_cutoff);
      warm_state.previous_projection_candidate = best_safe_candidate;
      warm_state.has_projection_candidate = true;
      apply_cut_candidate_timed(best_safe_points);
      continue;
    }

    if ((!best_safe_candidate.valid || best_safe_candidate.ultra_clash) &&
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
        ++timing.hull_rescue_successes;
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
        accept_rescue_candidate(best_hull_points, best_hull_candidate);
        continue;
      }

      ++timing.distance_rescue_attempts;
      vector<Vec3> rescued_points;
      rescued_points.reserve(base_atoms.size());
      Int3 rescue_shift{};
      double rescue_raw_delta_norm = numeric_limits<double>::infinity();
      if (rescue_split_contacting_owner_candidate(
              split_rescue_points, base_atoms, this->m_mol[frame_idx].Box,
              local_contacting_owner_pairs, shifts, split_rescue_detail,
              rescued_points, rescue_shift, rescue_raw_delta_norm)) {
        PbctopoCutCandidate rescued_candidate;
        if (validate_rescue_candidate_points(rescued_points, nullptr,
                                             rescued_candidate)) {
          ++timing.distance_rescue_successes;
          append_pbctopo_log_line(
              "[pbctopo] distance rescue candidate at frame " +
              to_string(static_cast<long long>(
                  this->m_mol[frame_idx].time_of_snap)) +
              " from " + split_rescue_source + ": applied_shift=" +
              format_int3(rescue_shift) + ", rescued_raw_delta_norm=" +
              to_string(rescue_raw_delta_norm));
          accept_rescue_candidate(rescued_points, rescued_candidate);
          continue;
        }
      }
    }

    if (!best_candidate.valid || best_candidate.ultra_clash) {
      const string fallback_message =
          "[pbctopo] frame " +
          to_string(static_cast<long long>(this->m_mol[frame_idx].time_of_snap)) +
          " fell back to pbci=0 after unsafe virtual-box search: safe_candidates=" +
          to_string(safe_candidates) + ", best.min_pair_d2=" +
          to_string(best_candidate.min_pair_d2) + ", best.kappa2=" +
          to_string(best_candidate.shape_kappa2) + ", best.asph=" +
          to_string(best_candidate.shape_asphericity) + ", best.bbox=" +
          to_string(best_candidate.bbox_volume);
      append_pbctopo_log_line(fallback_message);
      TitanUI::print_warning(
          fallback_message);
      ++fallback_frames;
      projection_clash_skip_active = false;
      clear_component_trial_hints(warm_state);
      fallback_pbci_frame(frame_idx);
      continue;
    }

    projection_clash_skip_active =
        candidate_has_projection_clash_margin(best_candidate,
                                              runtime_cutoffs.clash_cutoff);
    warm_state.previous_projection_candidate = best_candidate;
    warm_state.has_projection_candidate = best_candidate.valid;
    apply_cut_candidate_timed(best_points);
  }

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
      to_string(timing.graph_s) + " s, projection_edges=" +
      to_string(timing.projection_edge_s) + " s, component_split=" +
      to_string(timing.component_split_s) + " s, seeds=" +
      to_string(timing.seed_s) + " s, event_planes=" +
      to_string(timing.event_planes_s) + " s, coarse=" +
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
      " s, component_graph_refine_cut_score=" +
      to_string(timing.component_graph_refine_cut_score_s) +
      " s, component_graph_refine_ultraclash=" +
      to_string(timing.component_graph_refine_ultraclash_s) +
      " s, component_graph_refine_stats=" +
      to_string(timing.component_graph_refine_stats_s) +
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
      ", early_projection_used=" + to_string(timing.early_projection_used) +
      ", previous_seed_used=" + to_string(timing.previous_seed_used) +
       ", component_graph_edges=" + to_string(timing.component_graph_edges) +
       ", component_graph_shift_descriptors=" +
       to_string(timing.component_graph_shift_descriptors) +
       ", component_graph_shift_materialized=" +
       to_string(timing.component_graph_shift_materialized) +
       ", component_graph_bbox_only_scores=" +
       to_string(timing.component_graph_bbox_only_scores) +
        ", component_graph_used=" + to_string(timing.component_graph_used) +
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
      ", component_large_mode_merge_frames=" +
      to_string(timing.component_large_mode_merge_frames) +
      ", component_large_mode_pair_only_frames=" +
      to_string(timing.component_large_mode_pair_only_frames) +
      ", component_large_mode_beam_frames=" +
      to_string(timing.component_large_mode_beam_frames) +
      ", component_large_mode_merge_count=" +
      to_string(timing.component_large_mode_merge_count) +
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
      ", component_graph_outlier_rejects=" +
      to_string(timing.component_graph_outlier_rejects) +
      ", component_search_outlier_rejects=" +
      to_string(timing.component_search_outlier_rejects) +
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
  append_pbctopo_log_line("[pbctopo] Valid N-atom block candidates: " +
                          to_string(valid_candidates) +
                          ", rejected partial/invalid components: " +
                          to_string(partial_components) +
                          ", hull fallbacks: " +
                          to_string(hull_fallback_candidates));
  append_pbctopo_log_line("[pbctopo] end run");
  return 0;
}
