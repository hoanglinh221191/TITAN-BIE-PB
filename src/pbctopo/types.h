#ifndef TITAN_PBCTOPO_TYPES_H
#define TITAN_PBCTOPO_TYPES_H

#include "../pbctopo_hull.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace titan_pbctopo {

using pbctopo_hull::Vec3;

struct PbctopoAnchorGraphScore {
  double order_max_d2 = std::numeric_limits<double>::infinity();
  double order_path_d2 = std::numeric_limits<double>::infinity();
  double endpoint_image_gap = std::numeric_limits<double>::infinity();
  std::size_t residue_order_violations = 0;
  std::size_t contact_map_lost = 0;
  double contact_map_max_d2 = std::numeric_limits<double>::infinity();
  std::size_t ambiguity_count = 0;
};

struct PbctopoAtom {
  std::size_t global_atom_id = 0;
  std::size_t fragment_id = 0;
  int owner = 0;
  std::size_t local_atom_index = 0;
  std::uint64_t chain_key = 0;
  int residue_id = std::numeric_limits<int>::min();
  bool has_topology_key = false;
  bool is_heavy = true;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double sx = 0.0;
  double sy = 0.0;
  double sz = 0.0;
};

struct PbctopoCutCandidate {
  double cut_x = 0.0;
  double cut_y = 0.0;
  double cut_z = 0.0;
  double gap_score = 0.0;
  double bbox_volume = std::numeric_limits<double>::infinity();
  double rg2 = std::numeric_limits<double>::infinity();
  double shape_kappa2 = std::numeric_limits<double>::infinity();
  double shape_asphericity = std::numeric_limits<double>::infinity();
  double shape_acylindricity = std::numeric_limits<double>::infinity();
  double shift_norm2 = std::numeric_limits<double>::infinity();
  double min_pair_d2 = std::numeric_limits<double>::infinity();
  double continuity_path = std::numeric_limits<double>::infinity();
  double continuity_max_d2 = std::numeric_limits<double>::infinity();
  double anchor_mst2 = std::numeric_limits<double>::infinity();
  double anchor_max_d2 = std::numeric_limits<double>::infinity();
  PbctopoAnchorGraphScore anchor_graph;
  double principal_rg2_max = std::numeric_limits<double>::infinity();
  double principal_rg2_mid = std::numeric_limits<double>::infinity();
  double principal_rg2_min = std::numeric_limits<double>::infinity();
  double edge_stretch = std::numeric_limits<double>::infinity();
  double max_owner_raw_minus_minimg_delta = 0.0;
  std::size_t broken_edges = std::numeric_limits<std::size_t>::max();
  std::size_t split_owner_pairs = 0;
  std::size_t clash_atom_a = std::numeric_limits<std::size_t>::max();
  std::size_t clash_atom_b = std::numeric_limits<std::size_t>::max();
  bool ultra_clash = false;
  bool valid = false;
};

struct PbctopoLocalEdge {
  std::size_t a = 0;
  std::size_t b = 0;
  double ref_d2 = 0.0;
  double inv_ref_d2 = 0.0;
  double break_cutoff_d2 = 0.0;
};

struct PbctopoContinuityPair {
  std::size_t a = 0;
  std::size_t b = 0;
};

struct PbctopoComponentCandidate {
  double continuity_path = std::numeric_limits<double>::infinity();
  double continuity_max_d2 = std::numeric_limits<double>::infinity();
  double anchor_mst2 = std::numeric_limits<double>::infinity();
  double anchor_max_d2 = std::numeric_limits<double>::infinity();
  PbctopoAnchorGraphScore anchor_graph;
  double principal_rg2_max = std::numeric_limits<double>::infinity();
  double principal_rg2_mid = std::numeric_limits<double>::infinity();
  double principal_rg2_min = std::numeric_limits<double>::infinity();
  double shape_kappa2 = std::numeric_limits<double>::infinity();
  double shape_asphericity = std::numeric_limits<double>::infinity();
  double shape_acylindricity = std::numeric_limits<double>::infinity();
  double component_mst2 = std::numeric_limits<double>::infinity();
  double component_rg2 = std::numeric_limits<double>::infinity();
  double component_kappa2 = std::numeric_limits<double>::infinity();
  double component_asphericity = std::numeric_limits<double>::infinity();
  double component_acylindricity = std::numeric_limits<double>::infinity();
  double bbox_volume = std::numeric_limits<double>::infinity();
  double rg2 = std::numeric_limits<double>::infinity();
  double shift_norm2 = std::numeric_limits<double>::infinity();
  double min_pair_d2 = std::numeric_limits<double>::infinity();
  std::size_t clash_atom_a = std::numeric_limits<std::size_t>::max();
  std::size_t clash_atom_b = std::numeric_limits<std::size_t>::max();
  bool ultra_clash = false;
  bool valid = false;
};

struct Int3 {
  int x = 0;
  int y = 0;
  int z = 0;
};

struct PbctopoTimingStats {
  std::size_t frames = 0;
  std::size_t coarse_candidates = 0;
  std::size_t refine_candidates = 0;
  std::size_t projection_clash_skip_frames = 0;
  std::size_t projection_clash_skip_evals = 0;
  std::size_t component_frames = 0;
  std::size_t component_shift_trials = 0;
  std::size_t component_graph_refine_coarse_score_calls = 0;
  std::size_t component_graph_refine_coarse_stage1_score_calls = 0;
  std::size_t component_graph_refine_coarse_stage2_score_calls = 0;
  std::size_t component_graph_refine_exact_rescore_calls = 0;
  std::size_t component_partial_checks = 0;
  std::size_t component_partial_prunes = 0;
  std::size_t component_bbox_skips = 0;
  std::size_t component_score_calls = 0;
  std::size_t component_pair_table_builds = 0;
  std::size_t component_pair_table_prunes = 0;
  std::size_t component_search_zero_current_mask_prunes = 0;
  std::size_t component_search_zero_future_mask_prunes = 0;
  std::size_t component_search_fixed_minor_components = 0;
  std::size_t component_search_fixed_minor_frames = 0;
  std::size_t component_search_shift_warm_hits = 0;
  std::size_t component_search_seed_warm = 0;
  std::size_t component_search_seed_graph = 0;
  std::size_t component_search_seed_none = 0;
  std::size_t component_large_mode_merge_frames = 0;
  std::size_t component_large_mode_pair_only_frames = 0;
  std::size_t component_large_mode_beam_frames = 0;
  std::size_t component_large_mode_merge_count = 0;
  std::size_t component_large_mode_graph_seed_count = 0;
  std::size_t component_large_mode_graph_seed_refined = 0;
  std::size_t component_large_mode_graph_seed_nonbaseline_win = 0;
  std::size_t component_large_mode_beam_expanded = 0;
  std::size_t component_large_mode_beam_kept = 0;
  std::size_t component_large_mode_beam_complete_assignments = 0;
  std::size_t component_large_mode_beam_final_rescored = 0;
  std::size_t component_large_mode_fallback_frames = 0;
  std::size_t component_graph_outlier_rejects = 0;
  std::size_t component_search_outlier_rejects = 0;
  std::size_t component_score_continuity_rejects = 0;
  std::size_t component_score_anchor_rejects = 0;
  std::size_t component_score_prefix_prunes_after_continuity = 0;
  std::size_t component_score_prefix_prunes_after_anchor = 0;
  std::size_t component_score_prefix_prunes_after_principal = 0;
  std::size_t component_score_prefix_prunes_after_component_metrics = 0;
  std::size_t component_score_prefix_prunes_after_final = 0;
  std::size_t component_graph_edges = 0;
  std::size_t component_graph_shift_descriptors = 0;
  std::size_t component_graph_shift_materialized = 0;
  std::size_t component_graph_bbox_only_scores = 0;
  std::size_t component_graph_used = 0;
  std::size_t component_warm_used = 0;
  std::size_t component_warm_skip_search = 0;
  std::size_t component_search_frames = 0;
  std::size_t component_graph_failed = 0;
  std::size_t component_graph_points_mismatch = 0;
  std::size_t component_graph_invalid = 0;
  std::size_t component_graph_ultraclash = 0;
  std::size_t component_graph_not_better_not_close = 0;
  std::size_t component_refine_warm_shift_hits = 0;
  std::size_t component_refine_shortlist_hint_hits = 0;
  std::size_t component_pair_warm_shift_hits = 0;
  std::size_t component_refine_fanout_reduced = 0;
  std::size_t hull_rescue_attempts = 0;
  std::size_t hull_rescue_successes = 0;
  std::size_t distance_rescue_attempts = 0;
  std::size_t distance_rescue_successes = 0;
  std::size_t pbci_fallback_frames = 0;
  std::size_t component_graph_refine_score_calls = 0;
  std::size_t component_graph_refine_preclash_prunes = 0;
  std::size_t component_graph_refine_trial_build_incremental = 0;
  std::size_t component_graph_refine_trial_build_full = 0;
  std::size_t component_graph_refine_trial_build_global_rebase = 0;
  std::size_t component_graph_refine_stats_cached = 0;
  std::size_t component_graph_refine_stats_full = 0;
  std::size_t component_graph_refine_anchor_spatial_calls = 0;
  std::size_t component_graph_refine_anchor_spatial_fallbacks = 0;
  std::size_t component_graph_refine_anchor_spatial_disconnects = 0;
  std::size_t pbctopo_anchor_order_edges = 0;
  std::size_t pbctopo_anchor_endpoint_edges = 0;
  std::size_t pbctopo_anchor_contact_pairs = 0;
  std::size_t pbctopo_anchor_residue_order_violations = 0;
  std::size_t pbctopo_anchor_ambiguous_frames = 0;
  std::size_t component_cycle_sync_attempts = 0;
  std::size_t component_cycle_sync_used = 0;
  std::size_t component_cycle_sync_pair_candidates = 0;
  std::size_t component_cycle_sync_cycle_support = 0;
  double component_cycle_sync_cycle_mismatch = 0.0;
  std::size_t component_cycle_sync_inconsistent_edges = 0;
  std::size_t pbctopo_topology_bonds_available = 0;
  std::size_t pbctopo_topology_bond_edges_added = 0;
  std::size_t pbctopo_topology_bond_edges_skipped = 0;
  std::size_t component_context_builds = 0;
  std::size_t component_offset_extracts = 0;
  std::size_t apply_calls = 0;
  std::size_t early_projection_used = 0;
  std::size_t previous_seed_used = 0;
  std::array<std::size_t, 8> component_search_depth_visits{};
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
  double pbctopo_anchor_context_build_s = 0.0;
  double component_graph_refine_anchor_order_s = 0.0;
  double component_graph_refine_anchor_contact_s = 0.0;
  double component_graph_refine_cut_score_s = 0.0;
  double component_graph_refine_ultraclash_s = 0.0;
  double component_graph_refine_stats_s = 0.0;
  double component_cycle_sync_s = 0.0;
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

struct PbctopoAnchorContactEdge {
  std::size_t a = 0;
  std::size_t b = 0;
  double cutoff_d2 = std::numeric_limits<double>::infinity();
};

struct PbctopoAnchorGraphContext {
  std::vector<PbctopoContinuityPair> ordered_edges;
  std::vector<PbctopoContinuityPair> endpoint_edges;
  std::vector<PbctopoAnchorContactEdge> contact_edges;
  std::size_t anchor_count = 0;
  std::size_t residue_order_violations = 0;
  double box_m[3][3]{};
  double box_inv[3][3]{};
  bool has_box = false;
};

} // namespace titan_pbctopo

#endif
