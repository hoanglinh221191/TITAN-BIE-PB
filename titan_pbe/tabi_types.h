#ifndef TITAN_PBE_TABI_TYPES_H
#define TITAN_PBE_TABI_TYPES_H

#include "../aligned_allocator.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

class TitanExecutor;

namespace TitanPBE::Tabi {

struct Atom {
  double x{0.0};
  double y{0.0};
  double z{0.0};
  double q{0.0};
  double radius{0.0};
};

struct Panel {
  double x{0.0};
  double y{0.0};
  double z{0.0};
  double nx{0.0};
  double ny{0.0};
  double nz{0.0};
  double area{0.0};
  double diagonal{1.0};
  double solid_angle{0.0};
  std::size_t atom_index{0};
  double v0x{0.0};
  double v0y{0.0};
  double v0z{0.0};
  double v1x{0.0};
  double v1y{0.0};
  double v1z{0.0};
  double v2x{0.0};
  double v2y{0.0};
  double v2z{0.0};
  std::size_t v0_index{static_cast<std::size_t>(-1)};
  std::size_t v1_index{static_cast<std::size_t>(-1)};
  std::size_t v2_index{static_cast<std::size_t>(-1)};
  bool has_triangle_geometry{false};
};

struct Vertex {
  double x{0.0};
  double y{0.0};
  double z{0.0};
  double nx{0.0};
  double ny{0.0};
  double nz{0.0};
  std::size_t atom_index{0};
};

struct Triangle {
  std::size_t v0{0};
  std::size_t v1{0};
  std::size_t v2{0};
  std::size_t atom_index{0};
  double cx{0.0};
  double cy{0.0};
  double cz{0.0};
  double nx{0.0};
  double ny{0.0};
  double nz{0.0};
  double area{0.0};
};

struct SurfaceMesh {
  std::vector<Vertex> vertices;
  std::vector<Triangle> triangles;
  std::vector<Panel> panels;
  std::size_t exposed_candidates{0};
  std::size_t buried_candidates{0};
  std::size_t solid_angle_fallbacks{0};
  std::size_t nonmanifold_edges{0};
  double avg_valence{0.0};
  double diag_min{1.0};
  double diag_max{1.0};
  std::size_t cell_count{0};
  std::size_t max_neighbor_count{0};
  double mean_neighbor_count{0.0};
  double neighbor_build_ms{0.0};
  double mesh_build_ms{0.0};
  double solid_angle_ms{0.0};
  std::string sphere_mesh_mode{"latlong"};
  double quality_total_area_a2{0.0};
  double quality_area_rel_error_vs_4piR2{0.0};
  double quality_triangle_area_min{0.0};
  double quality_triangle_area_mean{0.0};
  double quality_triangle_area_max{0.0};
  double quality_triangle_aspect_max{0.0};
  double quality_normal_radial_dot_min{0.0};
  double quality_normal_radial_dot_mean{0.0};
  double quality_centroid_radial_error_max{0.0};
};

inline constexpr std::size_t kNearBlockCount = 4;
inline constexpr std::size_t kNearPairTypeCount = 4;
inline constexpr std::size_t kNearOrientationCount = 3;
inline constexpr std::size_t kNearBlockAuditBucketCount =
    kNearBlockCount * kNearPairTypeCount * kNearOrientationCount;

struct NearBlockAuditBucket {
  std::size_t accepted{0};
  std::size_t rejected{0};
  std::size_t candidate_worse_than_centroid{0};
  std::size_t reference_nonconverged{0};
  std::size_t orientation_sign_mismatch{0};
  std::size_t strict_ref_reject{0};
  double max_candidate_ref_rel{0.0};
  double sum_candidate_ref_rel{0.0};
  double max_centroid_ref_rel{0.0};
  std::size_t k2_laplace_sign_ok{0};
  double max_k2_yukawa_smooth_ref_rel{0.0};
  double max_k2_orientation_pair_ref_rel{0.0};
  std::size_t sample_count{0};
};

using NearBlockAuditBuckets =
    std::array<NearBlockAuditBucket, kNearBlockAuditBucketCount>;

inline constexpr std::size_t kGeometricGapBucketCount = 4;
inline constexpr std::size_t kGeometricCenterBucketCount = 3;
inline constexpr std::size_t kGeometricKernelCount = 3;
inline constexpr std::size_t kGeometricNearAuditBucketCount =
    kGeometricGapBucketCount * kGeometricCenterBucketCount *
    kGeometricKernelCount;

struct GeometricNearAuditBucket {
  std::size_t pairs{0};
  double min_gap{0.0};
  double max_gap_ratio{0.0};
  double sum_gap_ratio{0.0};
  double sum_center_ratio{0.0};
};

using GeometricNearAuditBuckets =
    std::array<GeometricNearAuditBucket, kGeometricNearAuditBucketCount>;

struct Options {
  double eps_in{1.0};
  double eps_out{80.0};
  double salt_m{0.15};
  double temperature_k{298.15};
  int max_iter{10000};
  double residual_tol{1.0e-6};
  int gmres_restart{32};
  std::string sum_mode{"auto"};
  double probe_radius{0.0};
  std::string surface_mode{"vdw"};
  std::string surface_geometry{"triangle"};
  std::string diagonal_mode{"solid-angle"};
  int surface_point_count{96};
  std::string sphere_mesh_mode{"auto"};
  int tree_leaf_size{32};
  double tree_opening_theta{0.45};
  int taylor_order{3};
  std::string derivative_mode{"analytic"};
  std::string formulation{"jcp2013"};
  std::string energy_projection{"boundary"};
  std::string preconditioner{"auto"};
  std::string preconditioner_block_near{"singular"};
  std::string near_integration{"auto"};
  std::string quad_buffer{"auto"};
  std::string leaf_near_policy{"auto"};
  std::string leaf_panel_kernel{"auto"};
  std::string adaptive_level{"strict"};
  std::string cluster_policy{"off"};
  std::string interaction_mode{"auto"};
  std::string node_bounds{"moment"};
  std::string tree_builder{"grouped-morton"};
  std::string mac_safety{"off"};
  double mac_gap_factor{1.5};
  std::string moment_audit{"off"};
  std::string matvec_audit{"off"};
  int matvec_audit_sample_count{1024};
  std::string progress{"off"};
  int progress_interval_seconds{5};
  std::string true_residual_probe{"restart"};
  std::string true_residual_probe_scope{"full"};
  int true_residual_probe_interval{8};
  int true_residual_sample_count{1024};
};

struct Result {
  std::string status{"ok"};
  std::string reaction_energy_source{"tabi_boundary_integral"};
  std::string convergence_reason{"not_started"};
  std::string residual_history{};
  std::size_t panel_count{0};
  std::size_t memory_bytes{0};
  std::size_t peak_memory_bytes{0};
  int surface_point_count_requested{96};
  double surface_area_a2{0.0};
  double panel_density_per_a2{0.0};
  std::string panel_budget_status{"off"};
  double surface_build_ms{0.0};
  double surface_neighbor_ms{0.0};
  double surface_mesh_ms{0.0};
  double solid_angle_ms{0.0};
  std::string sphere_mesh_mode{"latlong"};
  double surface_quality_total_area_a2{0.0};
  double surface_quality_area_rel_error_vs_4piR2{0.0};
  double surface_quality_triangle_area_min{0.0};
  double surface_quality_triangle_area_mean{0.0};
  double surface_quality_triangle_area_max{0.0};
  double surface_quality_triangle_aspect_max{0.0};
  double surface_quality_normal_radial_dot_min{0.0};
  double surface_quality_normal_radial_dot_mean{0.0};
  double surface_quality_centroid_radial_error_max{0.0};
  std::size_t exposed_points{0};
  std::size_t buried_points{0};
  double total_charge{0.0};
  int iterations{0};
  double residual{0.0};
  int gmres_restart_requested{32};
  int gmres_restart_effective{32};
  double best_residual{0.0};
  int best_residual_iter{0};
  double final_residual{0.0};
  double best_preconditioned_residual{0.0};
  double final_preconditioned_residual{0.0};
  double true_residual_abs{0.0};
  double true_residual_rel{0.0};
  std::string true_residual_probe_mode{"restart"};
  std::string true_residual_probe_scope{"full"};
  int true_residual_probe_interval{8};
  int true_residual_sample_count{1024};
  std::size_t true_residual_probe_calls{0};
  double true_residual_probe_ms{0.0};
  double sampled_true_residual_rel{0.0};
  std::size_t sampled_true_residual_rows{0};
  std::size_t sampled_probe_calls{0};
  double sampled_probe_ms{0.0};
  std::size_t full_true_residual_confirm_calls{0};
  double full_true_residual_confirm_ms{0.0};
  bool converged{false};
  std::string sum_mode_effective{"direct"};
  double tree_build_ms{0.0};
  std::size_t tree_nodes{0};
  std::size_t tree_leaf_count{0};
  std::size_t near_interactions{0};
  std::size_t far_interactions{0};
  std::string matvec_mode{"direct"};
  double rhs_ms{0.0};
  double gmres_ms{0.0};
  double matvec_ms_total{0.0};
  double matvec_ms_avg{0.0};
  std::size_t matvec_calls{0};
  double tree_vs_direct_matvec_rel{0.0};
  std::size_t matvec_audit_rows{0};
  double matvec_audit_ms{0.0};
  double orthogonalization_ms{0.0};
  double preconditioner_ms{0.0};
  std::string surface_geometry{"triangle"};
  std::size_t mesh_vertices{0};
  std::size_t mesh_triangles{0};
  std::size_t boundary_nodes{0};
  double avg_valence{0.0};
  std::size_t nonmanifold_edges{0};
  std::size_t solid_angle_fallbacks{0};
  double diag_min{1.0};
  double diag_max{1.0};
  int taylor_order{3};
  std::string derivative_mode{"analytic"};
  double interaction_build_ms{0.0};
  double tree_eval_ms{0.0};
  double tree_moment_ms{0.0};
  double tree_traversal_ms{0.0};
  double near_eval_ms{0.0};
  double far_eval_ms{0.0};
  double singular_near_eval_ms{0.0};
  double leaf_near_eval_ms{0.0};
  double leaf_panel_eval_ms{0.0};
  double leaf_tile_eval_ms{0.0};
  double leaf_residual_eval_ms{0.0};
  double smooth_leaf_node_eval_ms{0.0};
  std::size_t singular_near_interactions{0};
  std::size_t leaf_near_interactions{0};
  std::size_t far_node_interactions{0};
  std::size_t leaf_panel_interactions{0};
  std::size_t leaf_tile_interactions{0};
  std::size_t leaf_tile_panel_interactions{0};
  std::size_t leaf_residual_panel_interactions{0};
  std::size_t smooth_leaf_node_interactions{0};
  std::size_t leaf_pruned_interactions{0};
  std::size_t adaptive_smooth_leaf_accepted{0};
  std::size_t adaptive_smooth_leaf_rejected{0};
  std::size_t adaptive_span_tile_accepted{0};
  std::size_t adaptive_span_tile_rejected{0};
  std::size_t adaptive_exact_leaf_fallbacks{0};
  std::size_t adaptive_reject_one_ring_overlap{0};
  std::size_t adaptive_reject_error_bound_failed{0};
  std::size_t adaptive_reject_p4_check_failed{0};
  std::size_t adaptive_reject_too_close{0};
  double leaf_work_reduction_ratio{0.0};
  std::size_t target_cluster_interactions{0};
  std::size_t target_cluster_targets{0};
  std::size_t target_cluster_panels_pruned{0};
  std::size_t target_cluster_rejected{0};
  double target_cluster_eval_ms{0.0};
  std::string adaptive_level{"strict"};
  std::string leaf_near_policy{"exact"};
  std::string leaf_panel_kernel{"scalar"};
  std::string interaction_mode{"cached"};
  std::string node_bounds_mode{"moment"};
  std::string tree_builder{"grouped-morton"};
  std::string mac_safety{"off"};
  double mac_gap_factor{1.5};
  std::size_t tree_depth_max{0};
  std::size_t mac_reject_radius_limited{0};
  double mac_accept_rate{0.0};
  double node_radius_max{0.0};
  double node_radius_mean{0.0};
  double leaf_panel_work_per_matvec{0.0};
  double far_node_work_per_matvec{0.0};
  std::size_t far_geometric_unsafe_count{0};
  std::size_t far_geometric_unsafe_panels{0};
  double far_geometric_min_gap{0.0};
  double far_geometric_gap_ratio_min{0.0};
  double far_geometric_gap_ratio_mean{0.0};
  std::size_t far_geometric_audit_samples{0};
  std::size_t mac_physical_gap_rejects{0};
  std::string moment_audit{"off"};
  std::size_t moment_audit_nodes{0};
  double moment_translation_rel_error_max{0.0};
  double moment_translation_rel_error_mean{0.0};
  std::string quad_buffer_auto_candidate{"off"};
  std::string leaf_tile_status{"off"};
  double leaf_panel_pairs_per_second{0.0};
  std::string near_integration{"centroid"};
  std::string quad_buffer{"off"};
  std::size_t near_self_pairs{0};
  std::size_t near_one_ring_pairs{0};
  std::size_t near_duffy_pairs{0};
  std::size_t near_centroid_pairs{0};
  std::size_t near_fallbacks{0};
  std::size_t near_geometric_buffer_pairs{0};
  std::size_t near_geometric_buffer_gauss_pairs{0};
  std::size_t near_geometric_buffer_adaptive_pairs{0};
  std::size_t near_geometric_buffer_centroid_pairs{0};
  double near_geometric_buffer_min_gap{0.0};
  double near_geometric_buffer_mean_gap_ratio{0.0};
  double near_geometric_buffer_max_edge_mean{0.0};
  double near_integration_ms{0.0};
  std::size_t duffy_subdivisions_max{0};
  double duffy_subdivisions_mean{0.0};
  std::string near_kernel_policy{"centroid"};
  std::size_t near_k1_accepted{0};
  std::size_t near_k1_rejected{0};
  std::size_t near_k2_accepted{0};
  std::size_t near_k2_rejected{0};
  std::size_t near_k3_accepted{0};
  std::size_t near_k3_rejected{0};
  std::size_t near_k4_accepted{0};
  std::size_t near_k4_rejected{0};
  std::size_t near_one_ring_k1_accepted{0};
  std::size_t near_one_ring_k1_rejected{0};
  std::size_t near_one_ring_k2_accepted{0};
  std::size_t near_one_ring_k2_rejected{0};
  std::size_t near_one_ring_k3_accepted{0};
  std::size_t near_one_ring_k3_rejected{0};
  std::size_t near_one_ring_k4_accepted{0};
  std::size_t near_one_ring_k4_rejected{0};
  std::string near_k1_policy{"centroid"};
  std::string near_k2_policy{"centroid"};
  std::string near_k3_policy{"centroid"};
  std::string near_k4_policy{"centroid"};
  NearBlockAuditBuckets near_block_audit{};
  GeometricNearAuditBuckets geometric_near_audit{};
  double energy_kcal_mol{0.0};
  double energy_kj_mol{0.0};
  double energy_kt{0.0};
  std::string formulation{"jcp2013"};
  std::string energy_projection{"boundary"};
  std::string preconditioner_requested{"auto"};
  std::string preconditioner_effective{"diagonal"};
  std::string preconditioner_block_near{"singular"};
  std::string surface_generator{"titan_sphere_patch"};
  std::string discretization{"centroid"};
  double rhs_scale{1.0};
  double dielectric_contrast{0.0};
  double operator_contrast{0.0};
  std::size_t unknown_count{0};
  std::string bie_system{"off"};
  double lhs_diag_phi{0.0};
  double lhs_diag_phin{0.0};
  std::string energy_equation{"off"};
  std::string kernel_convention{"off"};
  double source_norm_s1{0.0};
  double source_norm_s2{0.0};
  double phi_norm{0.0};
  double phin_norm{0.0};
  double kappa_inv_a{0.0};
  double coulomb_factor_kcal_mol{0.0};
  double boundary_energy_kcal_mol{0.0};
  double reference_energy_kcal_mol{0.0};
  double selected_energy_scale{1.0};
  bool fallback_used{false};
};

struct TreecodeStats {
  double tree_build_ms{0.0};
  std::size_t tree_nodes{0};
  std::size_t tree_leaf_count{0};
  std::size_t near_interactions{0};
  std::size_t far_interactions{0};
  std::string matvec_mode{"direct"};
  int taylor_order{3};
  std::string derivative_mode{"analytic"};
  double interaction_build_ms{0.0};
  double tree_eval_ms{0.0};
  double tree_moment_ms{0.0};
  double tree_traversal_ms{0.0};
  double near_eval_ms{0.0};
  double far_eval_ms{0.0};
  double singular_near_eval_ms{0.0};
  double leaf_near_eval_ms{0.0};
  double leaf_panel_eval_ms{0.0};
  double leaf_tile_eval_ms{0.0};
  double leaf_residual_eval_ms{0.0};
  double smooth_leaf_node_eval_ms{0.0};
  std::size_t singular_near_interactions{0};
  std::size_t leaf_near_interactions{0};
  std::size_t far_node_interactions{0};
  std::size_t leaf_panel_interactions{0};
  std::size_t leaf_tile_interactions{0};
  std::size_t leaf_tile_panel_interactions{0};
  std::size_t leaf_residual_panel_interactions{0};
  std::size_t smooth_leaf_node_interactions{0};
  std::size_t leaf_pruned_interactions{0};
  std::size_t adaptive_smooth_leaf_accepted{0};
  std::size_t adaptive_smooth_leaf_rejected{0};
  std::size_t adaptive_span_tile_accepted{0};
  std::size_t adaptive_span_tile_rejected{0};
  std::size_t adaptive_exact_leaf_fallbacks{0};
  std::size_t adaptive_reject_one_ring_overlap{0};
  std::size_t adaptive_reject_error_bound_failed{0};
  std::size_t adaptive_reject_p4_check_failed{0};
  std::size_t adaptive_reject_too_close{0};
  double leaf_work_reduction_ratio{0.0};
  std::size_t target_cluster_interactions{0};
  std::size_t target_cluster_targets{0};
  std::size_t target_cluster_panels_pruned{0};
  std::size_t target_cluster_rejected{0};
  double target_cluster_eval_ms{0.0};
  std::string adaptive_level{"strict"};
  std::string leaf_near_policy{"exact"};
  std::string leaf_panel_kernel{"scalar"};
  std::string interaction_mode{"cached"};
  std::string node_bounds_mode{"moment"};
  std::string tree_builder{"grouped-morton"};
  std::string mac_safety{"off"};
  double mac_gap_factor{1.5};
  std::size_t tree_depth_max{0};
  std::size_t mac_reject_radius_limited{0};
  double mac_accept_rate{0.0};
  double node_radius_max{0.0};
  double node_radius_mean{0.0};
  double leaf_panel_work_per_matvec{0.0};
  double far_node_work_per_matvec{0.0};
  std::size_t far_geometric_unsafe_count{0};
  std::size_t far_geometric_unsafe_panels{0};
  double far_geometric_min_gap{0.0};
  double far_geometric_gap_ratio_min{0.0};
  double far_geometric_gap_ratio_mean{0.0};
  std::size_t far_geometric_audit_samples{0};
  std::size_t mac_physical_gap_rejects{0};
  std::string moment_audit{"off"};
  std::size_t moment_audit_nodes{0};
  double moment_translation_rel_error_max{0.0};
  double moment_translation_rel_error_mean{0.0};
  std::string leaf_tile_status{"off"};
  double leaf_panel_pairs_per_second{0.0};
  std::size_t matvec_calls{0};
  std::string near_integration{"centroid"};
  std::string quad_buffer{"off"};
  std::size_t near_self_pairs{0};
  std::size_t near_one_ring_pairs{0};
  std::size_t near_duffy_pairs{0};
  std::size_t near_centroid_pairs{0};
  std::size_t near_fallbacks{0};
  std::size_t near_geometric_buffer_pairs{0};
  std::size_t near_geometric_buffer_gauss_pairs{0};
  std::size_t near_geometric_buffer_adaptive_pairs{0};
  std::size_t near_geometric_buffer_centroid_pairs{0};
  double near_geometric_buffer_min_gap{0.0};
  double near_geometric_buffer_mean_gap_ratio{0.0};
  double near_geometric_buffer_max_edge_mean{0.0};
  double near_integration_ms{0.0};
  std::size_t duffy_subdivisions_max{0};
  double duffy_subdivisions_mean{0.0};
  std::string near_kernel_policy{"centroid"};
  std::size_t near_k1_accepted{0};
  std::size_t near_k1_rejected{0};
  std::size_t near_k2_accepted{0};
  std::size_t near_k2_rejected{0};
  std::size_t near_k3_accepted{0};
  std::size_t near_k3_rejected{0};
  std::size_t near_k4_accepted{0};
  std::size_t near_k4_rejected{0};
  std::size_t near_one_ring_k1_accepted{0};
  std::size_t near_one_ring_k1_rejected{0};
  std::size_t near_one_ring_k2_accepted{0};
  std::size_t near_one_ring_k2_rejected{0};
  std::size_t near_one_ring_k3_accepted{0};
  std::size_t near_one_ring_k3_rejected{0};
  std::size_t near_one_ring_k4_accepted{0};
  std::size_t near_one_ring_k4_rejected{0};
  std::string near_k1_policy{"centroid"};
  std::string near_k2_policy{"centroid"};
  std::string near_k3_policy{"centroid"};
  std::string near_k4_policy{"centroid"};
  NearBlockAuditBuckets near_block_audit{};
  GeometricNearAuditBuckets geometric_near_audit{};
};

inline void merge_near_kernel_block_stats(TreecodeStats &dst,
                                          const TreecodeStats &src) noexcept {
  auto merge_geometric_means = [](std::size_t old_count,
                                  std::size_t add_count, double old_value,
                                  double add_value) noexcept {
    const std::size_t total = old_count + add_count;
    if (total == 0)
      return 0.0;
    return (old_value * static_cast<double>(old_count) +
            add_value * static_cast<double>(add_count)) /
           static_cast<double>(total);
  };
  const std::size_t old_geometric_count = dst.near_geometric_buffer_pairs;
  dst.near_geometric_buffer_pairs += src.near_geometric_buffer_pairs;
  dst.near_geometric_buffer_gauss_pairs +=
      src.near_geometric_buffer_gauss_pairs;
  dst.near_geometric_buffer_adaptive_pairs +=
      src.near_geometric_buffer_adaptive_pairs;
  dst.near_geometric_buffer_centroid_pairs +=
      src.near_geometric_buffer_centroid_pairs;
  if (src.near_geometric_buffer_pairs > 0) {
    if (old_geometric_count == 0 || dst.near_geometric_buffer_min_gap <= 0.0) {
      dst.near_geometric_buffer_min_gap =
          src.near_geometric_buffer_min_gap;
    } else if (src.near_geometric_buffer_min_gap > 0.0) {
      dst.near_geometric_buffer_min_gap =
          std::min(dst.near_geometric_buffer_min_gap,
                   src.near_geometric_buffer_min_gap);
    }
    dst.near_geometric_buffer_mean_gap_ratio = merge_geometric_means(
        old_geometric_count, src.near_geometric_buffer_pairs,
        dst.near_geometric_buffer_mean_gap_ratio,
        src.near_geometric_buffer_mean_gap_ratio);
    dst.near_geometric_buffer_max_edge_mean = merge_geometric_means(
        old_geometric_count, src.near_geometric_buffer_pairs,
        dst.near_geometric_buffer_max_edge_mean,
        src.near_geometric_buffer_max_edge_mean);
  }
  dst.near_k1_accepted += src.near_k1_accepted;
  dst.near_k1_rejected += src.near_k1_rejected;
  dst.near_k2_accepted += src.near_k2_accepted;
  dst.near_k2_rejected += src.near_k2_rejected;
  dst.near_k3_accepted += src.near_k3_accepted;
  dst.near_k3_rejected += src.near_k3_rejected;
  dst.near_k4_accepted += src.near_k4_accepted;
  dst.near_k4_rejected += src.near_k4_rejected;
  dst.near_one_ring_k1_accepted += src.near_one_ring_k1_accepted;
  dst.near_one_ring_k1_rejected += src.near_one_ring_k1_rejected;
  dst.near_one_ring_k2_accepted += src.near_one_ring_k2_accepted;
  dst.near_one_ring_k2_rejected += src.near_one_ring_k2_rejected;
  dst.near_one_ring_k3_accepted += src.near_one_ring_k3_accepted;
  dst.near_one_ring_k3_rejected += src.near_one_ring_k3_rejected;
  dst.near_one_ring_k4_accepted += src.near_one_ring_k4_accepted;
  dst.near_one_ring_k4_rejected += src.near_one_ring_k4_rejected;
  auto merge_policy = [](std::size_t accepted, std::size_t rejected,
                         std::string &dst_policy,
                         const std::string &src_policy) noexcept {
    if (accepted > 0 && rejected > 0) {
      dst_policy = "mixed_guarded";
    } else if (accepted > 0 && !src_policy.empty() &&
               src_policy != "centroid") {
      dst_policy = src_policy;
    } else if (rejected > 0 && !src_policy.empty() &&
               src_policy != "centroid") {
      dst_policy = src_policy;
    } else if (rejected > 0) {
      dst_policy = "centroid_guarded";
    }
  };
  merge_policy(dst.near_k1_accepted, dst.near_k1_rejected,
               dst.near_k1_policy, src.near_k1_policy);
  merge_policy(dst.near_k2_accepted, dst.near_k2_rejected,
               dst.near_k2_policy, src.near_k2_policy);
  merge_policy(dst.near_k3_accepted, dst.near_k3_rejected,
               dst.near_k3_policy, src.near_k3_policy);
  merge_policy(dst.near_k4_accepted, dst.near_k4_rejected,
               dst.near_k4_policy, src.near_k4_policy);
  for (std::size_t i = 0; i < dst.near_block_audit.size(); ++i) {
    auto &d = dst.near_block_audit[i];
    const auto &s = src.near_block_audit[i];
    d.accepted += s.accepted;
    d.rejected += s.rejected;
    d.candidate_worse_than_centroid += s.candidate_worse_than_centroid;
    d.reference_nonconverged += s.reference_nonconverged;
    d.orientation_sign_mismatch += s.orientation_sign_mismatch;
    d.strict_ref_reject += s.strict_ref_reject;
    if (s.max_candidate_ref_rel > d.max_candidate_ref_rel)
      d.max_candidate_ref_rel = s.max_candidate_ref_rel;
    d.sum_candidate_ref_rel += s.sum_candidate_ref_rel;
    if (s.max_centroid_ref_rel > d.max_centroid_ref_rel)
      d.max_centroid_ref_rel = s.max_centroid_ref_rel;
    d.k2_laplace_sign_ok += s.k2_laplace_sign_ok;
    if (s.max_k2_yukawa_smooth_ref_rel >
        d.max_k2_yukawa_smooth_ref_rel)
      d.max_k2_yukawa_smooth_ref_rel =
          s.max_k2_yukawa_smooth_ref_rel;
    if (s.max_k2_orientation_pair_ref_rel >
        d.max_k2_orientation_pair_ref_rel)
      d.max_k2_orientation_pair_ref_rel =
          s.max_k2_orientation_pair_ref_rel;
    d.sample_count += s.sample_count;
  }
  for (std::size_t i = 0; i < dst.geometric_near_audit.size(); ++i) {
    auto &d = dst.geometric_near_audit[i];
    const auto &s = src.geometric_near_audit[i];
    if (s.pairs == 0)
      continue;
    if (d.pairs == 0 || d.min_gap <= 0.0)
      d.min_gap = s.min_gap;
    else if (s.min_gap > 0.0)
      d.min_gap = std::min(d.min_gap, s.min_gap);
    d.max_gap_ratio = std::max(d.max_gap_ratio, s.max_gap_ratio);
    d.sum_gap_ratio += s.sum_gap_ratio;
    d.sum_center_ratio += s.sum_center_ratio;
    d.pairs += s.pairs;
  }
}

struct TreeNode {
  double cx{0.0};
  double cy{0.0};
  double cz{0.0};
  double radius{0.0};
  double area_sum{0.0};
  double min_x{0.0};
  double min_y{0.0};
  double min_z{0.0};
  double max_x{0.0};
  double max_y{0.0};
  double max_z{0.0};
  double max_edge_length{0.0};
  std::size_t first_child{0};
  std::size_t child_count{0};
  std::size_t first_panel{0};
  std::size_t panel_count{0};
  std::size_t depth{0};
};

struct TreeInteractionLists {
  std::vector<std::size_t> singular_offsets;
  std::vector<std::size_t> singular_indices;
  std::vector<std::size_t> leaf_offsets;
  std::vector<std::size_t> leaf_indices;
  std::vector<std::size_t> leaf_tile_offsets;
  std::vector<std::size_t> leaf_tile_node_indices;
  std::vector<std::size_t> leaf_tile_exclusion_offsets;
  std::vector<std::size_t> leaf_tile_exclusion_indices;
  std::vector<std::size_t> leaf_tile_span_offsets;
  std::vector<std::size_t> leaf_tile_span_firsts;
  std::vector<std::size_t> leaf_tile_span_counts;
  std::vector<std::size_t> smooth_leaf_offsets;
  std::vector<std::size_t> smooth_leaf_node_indices;
  std::vector<std::size_t> far_offsets;
  std::vector<std::size_t> far_node_indices;
  std::size_t singular_pair_count{0};
  std::size_t leaf_pair_count{0};
  std::size_t leaf_tile_node_count{0};
  std::size_t leaf_tile_panel_count{0};
  std::size_t smooth_leaf_node_count{0};
  std::size_t leaf_pruned_count{0};
  std::size_t adaptive_smooth_leaf_accepted{0};
  std::size_t adaptive_smooth_leaf_rejected{0};
  std::size_t adaptive_span_tile_accepted{0};
  std::size_t adaptive_span_tile_rejected{0};
  std::size_t adaptive_exact_leaf_fallbacks{0};
  std::size_t adaptive_reject_one_ring_overlap{0};
  std::size_t adaptive_reject_error_bound_failed{0};
  std::size_t adaptive_reject_p4_check_failed{0};
  std::size_t adaptive_reject_too_close{0};
  std::size_t target_cluster_interactions{0};
  std::size_t target_cluster_targets{0};
  std::size_t target_cluster_panels_pruned{0};
  std::size_t target_cluster_rejected{0};
  std::size_t far_node_count{0};
  std::size_t mac_reject_radius_limited{0};
  std::size_t far_geometric_unsafe_count{0};
  std::size_t far_geometric_unsafe_panels{0};
  double far_geometric_min_gap{0.0};
  double far_geometric_gap_ratio_min{0.0};
  double far_geometric_gap_ratio_sum{0.0};
  std::size_t far_geometric_audit_samples{0};
  std::size_t mac_physical_gap_rejects{0};
  double build_ms{0.0};

  [[nodiscard]] bool valid_for(std::size_t panel_count) const noexcept {
    return singular_offsets.size() == panel_count + 1 &&
           leaf_offsets.size() == panel_count + 1 &&
           leaf_tile_offsets.size() == panel_count + 1 &&
           leaf_tile_exclusion_offsets.size() ==
               leaf_tile_node_indices.size() + 1 &&
           leaf_tile_span_offsets.size() ==
               leaf_tile_node_indices.size() + 1 &&
           smooth_leaf_offsets.size() == panel_count + 1 &&
           far_offsets.size() == panel_count + 1;
  }
};

struct TreePanelCache {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> nx;
  std::vector<double> ny;
  std::vector<double> nz;
  std::vector<double> area;
  AlignedDoubleVector ordered_x;
  AlignedDoubleVector ordered_y;
  AlignedDoubleVector ordered_z;
  AlignedDoubleVector ordered_nx;
  AlignedDoubleVector ordered_ny;
  AlignedDoubleVector ordered_nz;
  AlignedDoubleVector ordered_area;
  std::vector<std::size_t> ordered_original_index;

  [[nodiscard]] bool valid_for(std::size_t panel_count) const noexcept {
    return x.size() == panel_count && y.size() == panel_count &&
           z.size() == panel_count && nx.size() == panel_count &&
           ny.size() == panel_count && nz.size() == panel_count &&
           area.size() == panel_count &&
           ordered_x.size() == panel_count &&
           ordered_y.size() == panel_count &&
           ordered_z.size() == panel_count &&
           ordered_nx.size() == panel_count &&
           ordered_ny.size() == panel_count &&
           ordered_nz.size() == panel_count &&
           ordered_area.size() == panel_count &&
           ordered_original_index.size() == panel_count;
  }
};

struct TreeWorkspace {
  std::vector<std::size_t> order;
  std::vector<std::uint64_t> morton;
  std::vector<TreeNode> nodes;
  std::vector<std::size_t> child_indices;
  std::size_t root_index{0};
  std::size_t leaf_count{0};
  std::size_t depth_max{0};
  double build_ms{0.0};
  double node_radius_max{0.0};
  double node_radius_mean{0.0};
  std::size_t leaf_size{32};
  double opening_theta{0.45};
  int taylor_order{3};
  std::string derivative_mode{"analytic"};
  std::string leaf_near_policy{"exact"};
  std::string leaf_panel_kernel{"auto"};
  std::string adaptive_level{"strict"};
  std::string cluster_policy{"off"};
  std::string interaction_mode{"cached"};
  std::string node_bounds_mode{"moment"};
  std::string tree_builder{"grouped-morton"};
  std::string quad_buffer{"off"};
  std::string mac_safety{"off"};
  double mac_gap_factor{1.5};
  std::string moment_audit{"off"};
  TreePanelCache panel_cache;
  TreeInteractionLists interactions;
  std::vector<std::pair<std::size_t, std::size_t>> moment_levels;
};

std::vector<Panel> build_surface(const std::vector<Atom> &atoms,
                                 const Options &options,
                                 Result *result = nullptr,
                                 TitanExecutor *executor = nullptr);

void apply_direct_operator(const std::vector<Panel> &panels,
                           const std::vector<double> &density,
                           double contrast, double kappa,
                           TitanExecutor &executor,
                           std::vector<double> &out,
                           TreecodeStats *stats = nullptr);

TreeWorkspace build_tree_workspace(const std::vector<Panel> &panels,
                                   TitanExecutor &executor,
                                   std::size_t leaf_size = 32,
                                   double opening_theta = 0.45,
                                   const std::string &leaf_near_policy =
                                       "auto",
                                   const std::string &leaf_panel_kernel =
                                       "auto",
                                   const std::string &adaptive_level =
                                       "strict",
                                   const std::string &quad_buffer = "off",
                                   const std::string &interaction_mode =
                                       "cached",
                                   const std::string &node_bounds_mode =
                                       "moment",
                                   const std::string &tree_builder =
                                       "grouped-morton",
                                   const std::string &mac_safety = "off",
                                   double mac_gap_factor = 1.5,
                                   const std::string &moment_audit = "off");

void apply_treecode_operator(const std::vector<Panel> &panels,
                             const std::vector<double> &density,
                             double contrast, double kappa,
                             const TreeWorkspace &tree,
                             TitanExecutor &executor,
                             std::vector<double> &out,
                             TreecodeStats *stats = nullptr);

void apply_jcp_treecode_operator(const std::vector<Panel> &panels,
                                 const std::vector<double> &unknowns,
                                 double eps_ratio, double kappa,
                                 double lhs_phi, double lhs_phin,
                                 bool coupled, const TreeWorkspace &tree,
                                 const std::string &near_integration,
                                 const std::string &quad_buffer,
                                 TitanExecutor &executor,
                                 std::vector<double> &out,
                                 TreecodeStats *stats = nullptr);

Result solve_linear_pb(const std::vector<Atom> &atoms, const Options &options,
                       TitanExecutor &executor);

} // namespace TitanPBE::Tabi

#endif // TITAN_PBE_TABI_TYPES_H
