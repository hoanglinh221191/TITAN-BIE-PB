#include "tabi_types.h"
#include "tabi_fmm_m2l_preconditioner.h"
#include "tabi_preconditioner.h"
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <psapi.h>
#else
#include <cstdio>
#endif

#include "tabi_near.h"

#include "../titan_executor.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace TitanPBE::Tabi {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kCoulombKcalPerMol = 332.063713299;
constexpr double kKcalToKJ = 4.184;
constexpr double kKBKJPerMolK = 0.00831446261815324;
constexpr std::size_t kParallelVectorThreshold = 32768;

struct FormulationConfig {
  std::string formulation{"jcp2013"};
  std::string energy_projection{"boundary"};
  double rhs_scale{1.0};
  double operator_contrast_scale{1.0};
  double boundary_energy_scale{1.0};
};

double debye_kappa_inv_a(double salt_m, double temp_k) noexcept {
  if (salt_m <= 0.0)
    return 0.0;
  const double temp_scale = std::sqrt(298.15 / std::max(1.0, temp_k));
  return 0.3291 * std::sqrt(salt_m) * temp_scale;
}

double dot(const std::vector<double> &a, const std::vector<double> &b) {
  double value = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i)
    value += a[i] * b[i];
  return value;
}

double l2_norm(const std::vector<double> &values) {
  return std::sqrt(std::max(0.0, dot(values, values)));
}

double rms_norm(const std::vector<double> &values) {
  if (values.empty())
    return 0.0;
  return l2_norm(values) / std::sqrt(static_cast<double>(values.size()));
}

std::size_t solver_workers(const TitanExecutor &executor,
                           std::size_t work_items) noexcept {
  if (work_items == 0)
    return 1;
  return std::max<std::size_t>(
      1, std::min<std::size_t>(executor.capacity(), work_items));
}

double dot_parallel(const std::vector<double> &a, const std::vector<double> &b,
                    TitanExecutor &executor) {
  const std::size_t n = a.size();
  if (n < kParallelVectorThreshold || executor.capacity() <= 1)
    return dot(a, b);
  const std::size_t workers = solver_workers(executor, n);
  std::vector<double> partial(workers, 0.0);
  executor.parallel_for(
      n, workers, kParallelVectorThreshold,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        double value = 0.0;
        for (std::size_t i = begin; i < end; ++i)
          value += a[i] * b[i];
        partial[worker] = value;
      });
  return std::accumulate(partial.begin(), partial.end(), 0.0);
}

double l2_norm_parallel(const std::vector<double> &values,
                        TitanExecutor &executor) {
  return std::sqrt(std::max(0.0, dot_parallel(values, values, executor)));
}

void parallel_for_each_index(std::size_t n, TitanExecutor &executor,
                             const std::function<void(std::size_t)> &fn) {
  if (n < kParallelVectorThreshold || executor.capacity() <= 1) {
    for (std::size_t i = 0; i < n; ++i)
      fn(i);
    return;
  }
  const std::size_t workers = solver_workers(executor, n);
  executor.parallel_for(
      n, workers, kParallelVectorThreshold,
      [&](std::size_t, std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i)
          fn(i);
      });
}

double elapsed_ms(std::chrono::steady_clock::time_point begin,
                  std::chrono::steady_clock::time_point end) noexcept {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

std::string normalize_jcp_direct_kernel_name(const std::string &value) {
  if (value == "soa-avx2" || value == "fused-soa-avx2")
    return value;
  return "baseline";
}

void refresh_matvec_cost_summary(Result &result) noexcept {
  const double traversal_seconds = result.tree_traversal_ms / 1000.0;
  result.leaf_panels_per_second =
      traversal_seconds > 0.0
          ? static_cast<double>(result.leaf_panel_interactions) /
                traversal_seconds
          : 0.0;
  result.clean_leaf_panels_per_second =
      result.leaf_pair_batch_clean_source_panels_per_second;
  result.tree_traversal_fraction =
      result.matvec_ms_total > 0.0
          ? result.tree_traversal_ms / result.matvec_ms_total
          : 0.0;
}

std::size_t current_working_set_bytes() noexcept {
#if defined(_WIN32)
  PROCESS_MEMORY_COUNTERS counters{};
  counters.cb = sizeof(counters);
  if (K32GetProcessMemoryInfo(GetCurrentProcess(), &counters,
                              sizeof(counters)))
    return static_cast<std::size_t>(counters.WorkingSetSize);
#else
  FILE *file = std::fopen("/proc/self/status", "r");
  if (file != nullptr) {
    char line[256];
    char label[64];
    long value_kb = 0;
    while (std::fgets(line, sizeof(line), file) != nullptr) {
      if (std::sscanf(line, "%63s %ld", label, &value_kb) == 2 &&
          std::string(label) == "VmRSS:") {
        std::fclose(file);
        return static_cast<std::size_t>(std::max<long>(0, value_kb)) * 1024u;
      }
    }
    std::fclose(file);
  }
#endif
  return 0;
}

std::string csv_escape_progress(const std::string &value) {
  bool needs_quotes = false;
  for (char c : value) {
    if (c == ',' || c == '"' || c == '\n' || c == '\r') {
      needs_quotes = true;
      break;
    }
  }
  if (!needs_quotes)
    return value;
  std::string escaped;
  escaped.reserve(value.size() + 2);
  escaped.push_back('"');
  for (char c : value) {
    if (c == '"')
      escaped.push_back('"');
    escaped.push_back(c);
  }
  escaped.push_back('"');
  return escaped;
}

class LongRunProgressWriter {
public:
  LongRunProgressWriter(const Options &options, std::size_t threads)
      : enabled_(options.progress != "off"),
        interval_(std::chrono::seconds(
            std::clamp(options.progress_interval_seconds, 1, 300))),
        start_(std::chrono::steady_clock::now()),
        last_write_(start_ - interval_),
        threads_(threads) {
    if (!enabled_)
      return;
    out_.open("pbe_tabi_longrun_progress.csv",
              std::ios::out | std::ios::trunc);
    if (!out_.is_open()) {
      enabled_ = false;
      return;
    }
    out_ << "elapsed_seconds,phase,gmres_iter,restart_iter,matvec_calls,"
            "residual,gmres_restart_requested,gmres_restart_effective,"
            "best_residual,best_residual_iter,final_residual,"
            "best_preconditioned_residual,final_preconditioned_residual,"
            "true_residual_abs,true_residual_rel,"
            "true_residual_probe_mode,true_residual_probe_interval,"
            "true_residual_probe_calls,true_residual_probe_seconds,"
            "true_residual_probe_scope,true_residual_sample_count,"
            "sampled_true_residual_rel,sampled_true_residual_rows,"
            "sampled_probe_calls,sampled_probe_seconds,"
            "full_true_residual_confirm_calls,full_probe_seconds,"
            "panels,unknown_count,taylor_order,tree_nodes,"
            "interaction_mode,tree_builder,node_bounds_mode,mac_safety,"
            "mac_gap_factor,expansion_center,"
            "fragmentation_safety,moment_construction,moment_sum,"
            "convergence_audit,convergence_safety,"
            "convergence_rho_cap,mac_aspect_cap,"
            "convergence_estimator,convergence_tail_threshold,"
            "convergence_estimator_error_cap,"
            "convergence_unstable_rho_cap,"
            "convergence_calibration,"
            "convergence_dielectric_amp,"
            "convergence_estimator_model,"
            "convergence_channel_calibration,"
            "convergence_split_policy,"
            "convergence_split_selectivity,"
            "convergence_split_absolute_floor,"
            "convergence_overreject_audit,"
            "convergence_overreject_detail,"
            "convergence_shape_audit,"
            "convergence_kernel_response_audit,"
            "channel_missed_semantics,"
            "accepted_miss_diagnostics,"
            "accepted_miss_absolute_floor,"
            "convergence_k3_accepted_rho_cap,"
            "mac_distance_policy,mac_surface_margin_factor,"
            "mac_surface_margin_policy,"
            "mac_surface_margin_aspect_threshold,"
            "mac_surface_margin_rho_trigger,"
            "mac_surface_margin_soft_aspect_threshold,"
            "mac_surface_margin_soft_rho_trigger,"
            "mac_surface_margin_soft_factor,"
            "convergence_estimator_samples,"
            "convergence_estimator_unstable_far_count,"
            "convergence_estimator_unstable_far_panels,"
            "convergence_estimator_split_rejects,"
            "convergence_estimator_split_reject_panels,"
            "convergence_low_abs_reaccepted_far_count,"
            "convergence_low_abs_reaccepted_far_panels,"
            "convergence_estimator_tail_ratio_max,"
            "convergence_estimator_tail_ratio_mean,"
            "convergence_estimator_error_max,"
            "convergence_estimator_error_mean,"
            "convergence_estimator_abs_error_max,"
            "convergence_estimator_abs_error_mean,"
            "convergence_tail_norm_p_max,"
            "convergence_tail_norm_p_mean,"
            "convergence_tail_norm_pm1_max,"
            "convergence_tail_norm_pm1_mean,"
            "convergence_effective_jcp_amp_max,"
            "convergence_effective_jcp_amp_mean,"
            "convergence_dielectric_amp_max,"
            "convergence_dielectric_amp_mean,"
            "convergence_k1_effective_jcp_amp_max,"
            "convergence_k2_effective_jcp_amp_max,"
            "convergence_k3_effective_jcp_amp_max,"
            "convergence_k4_effective_jcp_amp_max,"
            "convergence_k1_effective_jcp_amp_mean,"
            "convergence_k2_effective_jcp_amp_mean,"
            "convergence_k3_effective_jcp_amp_mean,"
            "convergence_k4_effective_jcp_amp_mean,"
            "convergence_k1_estimated_error_max,"
            "convergence_k2_estimated_error_max,"
            "convergence_k3_estimated_error_max,"
            "convergence_k4_estimated_error_max,"
            "convergence_k1_unstable_far_count,"
            "convergence_k2_unstable_far_count,"
            "convergence_k3_unstable_far_count,"
            "convergence_k4_unstable_far_count,"
            "convergence_k1_unstable_far_panels,"
            "convergence_k2_unstable_far_panels,"
            "convergence_k3_unstable_far_panels,"
            "convergence_k4_unstable_far_panels,"
            "convergence_k1_split_rejects,"
            "convergence_k2_split_rejects,"
            "convergence_k3_split_rejects,"
            "convergence_k4_split_rejects,"
            "convergence_k1_split_reject_panels,"
            "convergence_k2_split_reject_panels,"
            "convergence_k3_split_reject_panels,"
            "convergence_k4_split_reject_panels,"
            "convergence_k1_low_abs_reaccepted,"
            "convergence_k2_low_abs_reaccepted,"
            "convergence_k3_low_abs_reaccepted,"
            "convergence_k4_low_abs_reaccepted,"
            "convergence_k1_low_abs_reaccepted_panels,"
            "convergence_k2_low_abs_reaccepted_panels,"
            "convergence_k3_low_abs_reaccepted_panels,"
            "convergence_k4_low_abs_reaccepted_panels,"
            "convergence_k1_high_error_caught,"
            "convergence_k2_high_error_caught,"
            "convergence_k3_high_error_caught,"
            "convergence_k4_high_error_caught,"
            "convergence_k1_high_error_missed,"
            "convergence_k2_high_error_missed,"
            "convergence_k3_high_error_missed,"
            "convergence_k4_high_error_missed,"
            "convergence_k1_high_error_caught_by_split,"
            "convergence_k2_high_error_caught_by_split,"
            "convergence_k3_high_error_caught_by_split,"
            "convergence_k4_high_error_caught_by_split,"
            "convergence_estimator_high_error_caught,"
            "convergence_estimator_high_error_missed,"
            "convergence_estimator_high_error_caught_by_split,"
            "convergence_estimator_phantom_high_error_count,"
            "convergence_estimator_phantom_high_error_abs_max,"
            "convergence_overreject_count,"
            "convergence_overreject_panels,"
            "convergence_overreject_caught_count,"
            "convergence_overreject_caught_panels,"
            "convergence_overreject_abs_error_max,"
            "convergence_overreject_abs_error_mean,"
            "convergence_overreject_abs_error_p95,"
            "convergence_overreject_abs_error_p99,"
            "convergence_overreject_estimated_abs_error_max,"
            "convergence_overreject_estimated_abs_error_mean,"
            "convergence_overreject_required_scale_max,"
            "convergence_overreject_required_scale_p95,"
            "convergence_overreject_required_scale_p99,"
            "convergence_k1_overreject_count,"
            "convergence_k2_overreject_count,"
            "convergence_k3_overreject_count,"
            "convergence_k4_overreject_count,"
            "convergence_k1_overreject_caught_count,"
            "convergence_k2_overreject_caught_count,"
            "convergence_k3_overreject_caught_count,"
            "convergence_k4_overreject_caught_count,"
            "convergence_k1_overreject_panels,convergence_k2_overreject_panels,"
            "convergence_k3_overreject_panels,convergence_k4_overreject_panels,"
            "convergence_k1_overreject_abs_error_max,"
            "convergence_k2_overreject_abs_error_max,"
            "convergence_k3_overreject_abs_error_max,"
            "convergence_k4_overreject_abs_error_max,"
            "convergence_k1_overreject_abs_error_mean,"
            "convergence_k2_overreject_abs_error_mean,"
            "convergence_k3_overreject_abs_error_mean,"
            "convergence_k4_overreject_abs_error_mean,"
            "convergence_k1_overreject_required_scale_max,"
            "convergence_k2_overreject_required_scale_max,"
            "convergence_k3_overreject_required_scale_max,"
            "convergence_k4_overreject_required_scale_max,"
            "convergence_k1_overreject_required_scale_p99,"
            "convergence_k2_overreject_required_scale_p99,"
            "convergence_k3_overreject_required_scale_p99,"
            "convergence_k4_overreject_required_scale_p99,"
            "phantom_high_error_k1_count,phantom_high_error_k2_count,"
            "phantom_high_error_k3_count,phantom_high_error_k4_count,"
            "phantom_high_error_k1_abs_max,"
            "phantom_high_error_k2_abs_max,"
            "phantom_high_error_k3_abs_max,"
            "phantom_high_error_k4_abs_max,"
            "accepted_miss_worst_channel,"
            "accepted_miss_k1_count,accepted_miss_k2_count,"
            "accepted_miss_k3_count,accepted_miss_k4_count,"
            "accepted_miss_k1_panels,accepted_miss_k2_panels,"
            "accepted_miss_k3_panels,accepted_miss_k4_panels,"
            "accepted_miss_k1_rho_min,accepted_miss_k2_rho_min,"
            "accepted_miss_k3_rho_min,accepted_miss_k4_rho_min,"
            "accepted_miss_k1_rho_mean,accepted_miss_k2_rho_mean,"
            "accepted_miss_k3_rho_mean,accepted_miss_k4_rho_mean,"
            "accepted_miss_k1_rho_max,accepted_miss_k2_rho_max,"
            "accepted_miss_k3_rho_max,accepted_miss_k4_rho_max,"
            "accepted_miss_k1_tail_ratio_mean,"
            "accepted_miss_k2_tail_ratio_mean,"
            "accepted_miss_k3_tail_ratio_mean,"
            "accepted_miss_k4_tail_ratio_mean,"
            "accepted_miss_k1_tail_ratio_max,"
            "accepted_miss_k2_tail_ratio_max,"
            "accepted_miss_k3_tail_ratio_max,"
            "accepted_miss_k4_tail_ratio_max,"
            "accepted_miss_k1_estimated_error_mean,"
            "accepted_miss_k2_estimated_error_mean,"
            "accepted_miss_k3_estimated_error_mean,"
            "accepted_miss_k4_estimated_error_mean,"
            "accepted_miss_k1_estimated_error_max,"
            "accepted_miss_k2_estimated_error_max,"
            "accepted_miss_k3_estimated_error_max,"
            "accepted_miss_k4_estimated_error_max,"
            "convergence_k3_rho_cap_rejects,"
            "convergence_k3_rho_cap_reject_panels,"
            "expansion_audit,expansion_audit_channel,"
            "far_geometric_unsafe_count,far_geometric_unsafe_panels,"
            "far_geometric_audit_samples,mac_physical_gap_rejects,"
            "mac_convergence_rho_rejects,"
            "mac_convergence_aspect_rejects,"
            "mac_convergence_rho_reject_panels,"
            "mac_convergence_aspect_reject_panels,"
            "mac_surface_margin_rejects,"
            "mac_surface_margin_reject_panels,"
            "mac_surface_margin_applied,"
            "mac_surface_margin_skipped,"
            "mac_surface_margin_adaptive_rejects,"
            "mac_surface_margin_adaptive_reject_panels,"
            "mac_surface_margin_skipped_accepted_far,"
            "mac_surface_margin_hard_applied,"
            "mac_surface_margin_soft_applied,"
            "mac_surface_margin_hard_rejects,"
            "mac_surface_margin_hard_reject_panels,"
            "mac_surface_margin_soft_rejects,"
            "mac_surface_margin_soft_reject_panels,"
            "fragmented_far_count,fragmented_far_panels,"
            "fragmented_split_rejects,"
            "fragmentation_center_gap_ratio_max,"
            "fragmentation_center_gap_ratio_mean,"
            "moment_audit,moment_audit_nodes,"
            "moment_translation_rel_error_max,"
            "moment_translation_rel_error_mean,"
            "moment_translation_offset_max,"
            "moment_translation_offset_mean,"
            "moment_cancellation_ratio_max,"
            "moment_cancellation_ratio_mean,"
            "moment_dynamic_range_max,"
            "far_expansion_audit_pairs,far_expansion_audit_panels,"
            "far_expansion_direct_seconds,"
            "far_expansion_rel_error_max,far_expansion_rel_error_mean,"
            "far_expansion_rel_error_p95,far_expansion_rel_error_p99,"
            "far_expansion_worst_channel,"
            "far_expansion_k1_rel_error_max,"
            "far_expansion_k2_rel_error_max,"
            "far_expansion_k3_rel_error_max,"
            "far_expansion_k4_rel_error_max,"
            "far_expansion_k1_rel_error_mean,"
            "far_expansion_k2_rel_error_mean,"
            "far_expansion_k3_rel_error_mean,"
            "far_expansion_k4_rel_error_mean,"
            "far_expansion_k1_rel_error_p99,"
            "far_expansion_k2_rel_error_p99,"
            "far_expansion_k3_rel_error_p99,"
            "far_expansion_k4_rel_error_p99,"
            "far_expansion_mean_gate_pass,"
            "far_convergence_rho_max,far_convergence_rho_mean,"
            "node_aspect_ratio_max,node_aspect_ratio_mean,"
            "accepted_far_aspect_ratio_max,"
            "accepted_far_aspect_ratio_mean,"
            "leaf_panel_interactions,leaf_tile_panel_interactions,"
            "leaf_residual_panel_interactions,far_node_interactions,"
            "leaf_mac_candidates,leaf_mac_accepted,leaf_mac_rejected,"
            "leaf_mac_far_panels,leaf_mac_direct_panels_saved_estimate,"
            "matvec_cost_profile,leaf_panels_per_second,"
            "clean_leaf_panels_per_second,tree_traversal_fraction,"
            "leaf_work_profile,leaf_work_mac_terminal_panels,"
            "leaf_work_physical_gap_panels,"
            "leaf_work_fragmentation_panels,"
            "leaf_work_convergence_panels,"
            "leaf_work_convergence_estimator_panels,"
            "leaf_work_convergence_rho_panels,"
            "leaf_work_convergence_aspect_panels,"
            "leaf_work_surface_margin_panels,"
            "leaf_work_surface_margin_hard_panels,"
            "leaf_work_surface_margin_soft_panels,"
            "leaf_work_k3_rho_cap_panels,"
            "leaf_work_singular_fallback_panels,"
            "leaf_work_scalar_tail_panels,"
            "leaf_work_profile_accounted_panels,"
            "leaf_work_mac_terminal_frac,"
            "leaf_work_convergence_frac,"
            "leaf_work_surface_margin_frac,"
            "leaf_work_singular_fallback_frac,"
            "leaf_work_scalar_tail_frac,"
            "cluster_policy,target_batch_size,"
            "leaf_pair_min_source_panels,leaf_pair_min_targets,"
            "leaf_pair_batch_interactions,leaf_pair_batch_targets,"
            "leaf_pair_batch_source_panels,"
            "leaf_pair_batch_residual_panels,"
            "leaf_pair_batch_fallbacks,leaf_pair_batch_seconds,"
            "leaf_pair_batch_kernel,leaf_pair_batch_vectorized_sources,"
            "leaf_pair_batch_vector_lanes,"
            "leaf_pair_batch_pairs_per_second,"
            "leaf_kernel_profile,jcp_direct_kernel,"
            "leaf_kernel_scalar_tasks,leaf_kernel_batch_tasks,"
            "leaf_kernel_source_panel_evals,"
            "leaf_kernel_scalar_seconds,"
            "leaf_kernel_target_avx2_seconds,"
            "leaf_kernel_tail_seconds,"
            "leaf_pair_batch_soa_sources,"
            "leaf_pair_batch_soa_seconds,"
            "leaf_pair_batch_scalar_tail_seconds,"
            "leaf_pair_batch_clean_source_panels_per_second,"
            "leaf_pair_selection_profile,leaf_pair_probe_cache,"
            "leaf_pair_target_group,leaf_pair_max_lookahead,"
            "leaf_pair_directness_probe,leaf_pair_directness_cache,"
            "leaf_pair_directness_verify,"
            "leaf_pair_selection_candidate_scans,"
            "leaf_pair_selection_clean_rejects,"
            "leaf_pair_selection_probe_calls,"
            "leaf_pair_selection_probe_node_visits,"
            "leaf_pair_selection_probe_seconds,"
            "leaf_pair_selection_seconds,"
            "leaf_pair_selection_probe_cache_hits,"
            "leaf_pair_selection_probe_cache_misses,"
            "leaf_pair_selection_insufficient_targets,"
            "leaf_pair_selection_same_target_leaf_rejects,"
            "leaf_pair_directness_ancestor_probe_calls,"
            "leaf_pair_directness_ancestor_probe_node_visits,"
            "leaf_pair_directness_ancestor_probe_seconds,"
            "leaf_pair_directness_pair_cache_hits,"
            "leaf_pair_directness_pair_cache_misses,"
            "leaf_pair_directness_all_direct_proofs,"
            "leaf_pair_directness_mixed_fallbacks,"
            "leaf_pair_directness_verify_checks,"
            "leaf_pair_directness_verify_mismatches,"
            "leaf_pair_scheduler,leaf_pair_worklist_max_active,"
            "leaf_pair_worklist_verify,"
            "leaf_pair_worklist_tasks,leaf_pair_worklist_groups,"
            "leaf_pair_worklist_batched_groups,"
            "leaf_pair_worklist_scalar_flush_tasks,"
            "leaf_pair_worklist_forced_flushes,"
            "leaf_pair_worklist_active_groups_max,"
            "leaf_pair_worklist_seconds,"
            "leaf_pair_worklist_verify_checks,"
            "leaf_pair_worklist_verify_mismatches,"
            "leaf_block_verify,leaf_block_target_leaves,"
            "leaf_block_targets,leaf_block_source_groups,"
            "leaf_block_batched_groups,leaf_block_scalar_tail_tasks,"
            "leaf_block_group_seconds,leaf_block_flush_seconds,"
            "leaf_block_writeback_seconds,leaf_block_verify_checks,"
            "leaf_block_verify_mismatches,"
            "singular_near_interactions,preconditioner_effective,"
            "preconditioner_block_panel_cap,"
            "preconditioner_overlap_panel_cap,"
            "preconditioner_ilu_neighbor_cap,"
            "preconditioner_ilu_drop_tol,"
            "preconditioner_memory_budget_bytes,"
            "preconditioner_estimated_bytes,"
            "preconditioner_valid_blocks,preconditioner_valid_panels,"
            "preconditioner_ill_conditioned_blocks,"
            "preconditioner_quality_rejected_blocks,"
            "preconditioner_quality_probe_max_gain,"
            "preconditioner_quality_status,"
            "preconditioner_setup_seconds,preconditioner_apply_seconds,"
            "preconditioner_apply_calls,"
            "preconditioner_sparse_nnz_blocks,"
            "preconditioner_sparse_avg_degree,"
            "preconditioner_sparse_max_degree,"
            "preconditioner_coarse_aggregate_cap,"
            "preconditioner_coarse_aggregates,"
            "preconditioner_coarse_dimension,"
            "preconditioner_coarse_base,"
            "preconditioner_coarse_matrix_bytes,"
            "preconditioner_coarse_pivot_ratio,"
            "preconditioner_coarse_status,"
            "preconditioner_coarse_setup_matvec_calls,"
            "preconditioner_coarse_setup_seconds,"
            "preconditioner_coarse_apply_seconds,"
            "preconditioner_coarse_apply_calls,"
            "preconditioner_coarse_rp_identity_error_max,"
            "preconditioner_coarse_verify_rel_error_max,"
            "threads,"
            "rss_bytes,status\n";
    out_.flush();
  }

  [[nodiscard]] bool enabled() const noexcept { return enabled_; }

  void write(const std::string &phase, const Result &result,
             const TreecodeStats *stats, std::size_t matvec_calls,
             int gmres_iter, int restart_iter, double residual,
             const std::string &status, bool force) {
    if (!enabled_)
      return;
    const auto now = std::chrono::steady_clock::now();
    if (!force && now - last_write_ < interval_)
      return;
    last_write_ = now;

    const std::size_t tree_nodes =
        stats != nullptr && stats->tree_nodes > 0 ? stats->tree_nodes
                                                  : result.tree_nodes;
    const std::size_t leaf_panel =
        stats != nullptr ? stats->leaf_panel_interactions
                         : result.leaf_panel_interactions;
    const std::size_t leaf_tile_panel =
        stats != nullptr ? stats->leaf_tile_panel_interactions
                         : result.leaf_tile_panel_interactions;
    const std::size_t leaf_residual =
        stats != nullptr ? stats->leaf_residual_panel_interactions
                         : result.leaf_residual_panel_interactions;
    const std::size_t far_nodes =
        stats != nullptr ? stats->far_node_interactions
                         : result.far_node_interactions;
    const std::size_t leaf_mac_candidates =
        stats != nullptr ? stats->leaf_mac_candidates
                         : result.leaf_mac_candidates;
    const std::size_t leaf_mac_accepted =
        stats != nullptr ? stats->leaf_mac_accepted
                         : result.leaf_mac_accepted;
    const std::size_t leaf_mac_rejected =
        stats != nullptr ? stats->leaf_mac_rejected
                         : result.leaf_mac_rejected;
    const std::size_t leaf_mac_far_panels =
        stats != nullptr ? stats->leaf_mac_far_panels
                         : result.leaf_mac_far_panels;
    const std::size_t leaf_mac_direct_panels_saved_estimate =
        stats != nullptr ? stats->leaf_mac_direct_panels_saved_estimate
                         : result.leaf_mac_direct_panels_saved_estimate;
    const double traversal_ms =
        stats != nullptr ? stats->tree_traversal_ms : result.tree_traversal_ms;
    const double leaf_panels_per_second =
        traversal_ms > 0.0
            ? static_cast<double>(leaf_panel) / (traversal_ms / 1000.0)
            : result.leaf_panels_per_second;
    const std::string matvec_cost_profile = result.matvec_cost_profile;
    const double tree_traversal_fraction =
        result.matvec_ms_total > 0.0 ? traversal_ms / result.matvec_ms_total
                                     : result.tree_traversal_fraction;
    const std::size_t singular =
        stats != nullptr ? stats->singular_near_interactions
                         : result.singular_near_interactions;
    const std::string tree_builder =
        stats != nullptr ? stats->tree_builder : result.tree_builder;
    const std::string interaction_mode =
        stats != nullptr ? stats->interaction_mode : result.interaction_mode;
    const std::string node_bounds =
        stats != nullptr ? stats->node_bounds_mode : result.node_bounds_mode;
    const std::string mac_safety =
        stats != nullptr ? stats->mac_safety : result.mac_safety;
    const std::string expansion_center =
        stats != nullptr ? stats->expansion_center : result.expansion_center;
    const std::string fragmentation_safety =
        stats != nullptr ? stats->fragmentation_safety
                         : result.fragmentation_safety;
    const std::string expansion_audit =
        stats != nullptr ? stats->expansion_audit : result.expansion_audit;
    const double mac_gap_factor =
        stats != nullptr ? stats->mac_gap_factor : result.mac_gap_factor;
    const std::size_t unsafe_count =
        stats != nullptr ? stats->far_geometric_unsafe_count
                         : result.far_geometric_unsafe_count;
    const std::size_t unsafe_panels =
        stats != nullptr ? stats->far_geometric_unsafe_panels
                         : result.far_geometric_unsafe_panels;
    const std::size_t gap_samples =
        stats != nullptr ? stats->far_geometric_audit_samples
                         : result.far_geometric_audit_samples;
    const std::size_t gap_rejects =
        stats != nullptr ? stats->mac_physical_gap_rejects
                         : result.mac_physical_gap_rejects;
    const std::size_t fragmented_count =
        stats != nullptr ? stats->fragmented_far_count
                         : result.fragmented_far_count;
    const std::size_t fragmented_panels =
        stats != nullptr ? stats->fragmented_far_panels
                         : result.fragmented_far_panels;
    const std::size_t fragmented_rejects =
        stats != nullptr ? stats->fragmented_split_rejects
                         : result.fragmented_split_rejects;
    const double fragmentation_gap_max =
        stats != nullptr ? stats->fragmentation_center_gap_ratio_max
                         : result.fragmentation_center_gap_ratio_max;
    const double fragmentation_gap_mean =
        stats != nullptr ? stats->fragmentation_center_gap_ratio_mean
                         : result.fragmentation_center_gap_ratio_mean;
    const std::string moment_audit =
        stats != nullptr ? stats->moment_audit : result.moment_audit;
    const std::string moment_construction =
        stats != nullptr ? stats->moment_construction
                         : result.moment_construction;
    const std::string moment_sum =
        stats != nullptr ? stats->moment_sum : result.moment_sum;
    const std::string convergence_audit =
        stats != nullptr ? stats->convergence_audit
                         : result.convergence_audit;
    const std::string convergence_safety =
        stats != nullptr ? stats->convergence_safety
                         : result.convergence_safety;
    const double convergence_rho_cap =
        stats != nullptr ? stats->convergence_rho_cap
                         : result.convergence_rho_cap;
    const double mac_aspect_cap =
        stats != nullptr ? stats->mac_aspect_cap : result.mac_aspect_cap;
    const std::string convergence_estimator =
        stats != nullptr ? stats->convergence_estimator
                         : result.convergence_estimator;
    const double convergence_tail_threshold =
        stats != nullptr ? stats->convergence_tail_threshold
                         : result.convergence_tail_threshold;
    const double convergence_estimator_error_cap =
        stats != nullptr ? stats->convergence_estimator_error_cap
                         : result.convergence_estimator_error_cap;
    const double convergence_unstable_rho_cap =
        stats != nullptr ? stats->convergence_unstable_rho_cap
                         : result.convergence_unstable_rho_cap;
    const std::string convergence_calibration =
        stats != nullptr ? stats->convergence_calibration
                         : result.convergence_calibration;
    const std::string convergence_dielectric_amp =
        stats != nullptr ? stats->convergence_dielectric_amp
                         : result.convergence_dielectric_amp;
    const std::string convergence_estimator_model =
        stats != nullptr ? stats->convergence_estimator_model
                         : result.convergence_estimator_model;
    const std::string convergence_channel_calibration =
        stats != nullptr ? stats->convergence_channel_calibration
                         : result.convergence_channel_calibration;
    const std::string convergence_split_policy =
        stats != nullptr ? stats->convergence_split_policy
                         : result.convergence_split_policy;
    const std::string convergence_split_selectivity =
        stats != nullptr ? stats->convergence_split_selectivity
                         : result.convergence_split_selectivity;
    const double convergence_split_absolute_floor =
        stats != nullptr ? stats->convergence_split_absolute_floor
                         : result.convergence_split_absolute_floor;
    const std::string convergence_overreject_audit =
        stats != nullptr ? stats->convergence_overreject_audit
                         : result.convergence_overreject_audit;
    const std::string convergence_overreject_detail =
        stats != nullptr ? stats->convergence_overreject_detail
                         : result.convergence_overreject_detail;
    const std::string convergence_shape_audit =
        stats != nullptr ? stats->convergence_shape_audit
                         : result.convergence_shape_audit;
    const std::string convergence_kernel_response_audit =
        stats != nullptr ? stats->convergence_kernel_response_audit
                         : result.convergence_kernel_response_audit;
    const std::string channel_missed_semantics =
        stats != nullptr ? stats->channel_missed_semantics
                         : result.channel_missed_semantics;
    const std::string accepted_miss_diagnostics =
        stats != nullptr ? stats->accepted_miss_diagnostics
                         : result.accepted_miss_diagnostics;
    const double accepted_miss_absolute_floor =
        stats != nullptr ? stats->accepted_miss_absolute_floor
                         : result.accepted_miss_absolute_floor;
    const double convergence_k3_accepted_rho_cap =
        stats != nullptr ? stats->convergence_k3_accepted_rho_cap
                         : result.convergence_k3_accepted_rho_cap;
    const std::string mac_distance_policy =
        stats != nullptr ? stats->mac_distance_policy
                         : result.mac_distance_policy;
    const double mac_surface_margin_factor =
        stats != nullptr ? stats->mac_surface_margin_factor
                         : result.mac_surface_margin_factor;
    const std::string mac_surface_margin_policy =
        stats != nullptr ? stats->mac_surface_margin_policy
                         : result.mac_surface_margin_policy;
    const double mac_surface_margin_aspect_threshold =
        stats != nullptr ? stats->mac_surface_margin_aspect_threshold
                         : result.mac_surface_margin_aspect_threshold;
    const double mac_surface_margin_rho_trigger =
        stats != nullptr ? stats->mac_surface_margin_rho_trigger
                         : result.mac_surface_margin_rho_trigger;
    const double mac_surface_margin_soft_aspect_threshold =
        stats != nullptr ? stats->mac_surface_margin_soft_aspect_threshold
                         : result.mac_surface_margin_soft_aspect_threshold;
    const double mac_surface_margin_soft_rho_trigger =
        stats != nullptr ? stats->mac_surface_margin_soft_rho_trigger
                         : result.mac_surface_margin_soft_rho_trigger;
    const double mac_surface_margin_soft_factor =
        stats != nullptr ? stats->mac_surface_margin_soft_factor
                         : result.mac_surface_margin_soft_factor;
    const std::size_t convergence_estimator_samples =
        stats != nullptr ? stats->convergence_estimator_samples
                         : result.convergence_estimator_samples;
    const std::size_t convergence_estimator_unstable_count =
        stats != nullptr ? stats->convergence_estimator_unstable_far_count
                         : result.convergence_estimator_unstable_far_count;
    const std::size_t convergence_estimator_unstable_panels =
        stats != nullptr ? stats->convergence_estimator_unstable_far_panels
                         : result.convergence_estimator_unstable_far_panels;
    const std::size_t convergence_estimator_split_rejects =
        stats != nullptr ? stats->convergence_estimator_split_rejects
                         : result.convergence_estimator_split_rejects;
    const std::size_t convergence_estimator_split_reject_panels =
        stats != nullptr ? stats->convergence_estimator_split_reject_panels
                         : result.convergence_estimator_split_reject_panels;
    const std::size_t convergence_low_abs_reaccepted_count =
        stats != nullptr ? stats->convergence_low_abs_reaccepted_far_count
                         : result.convergence_low_abs_reaccepted_far_count;
    const std::size_t convergence_low_abs_reaccepted_panels =
        stats != nullptr ? stats->convergence_low_abs_reaccepted_far_panels
                         : result.convergence_low_abs_reaccepted_far_panels;
    const double convergence_estimator_tail_ratio_max =
        stats != nullptr ? stats->convergence_estimator_tail_ratio_max
                         : result.convergence_estimator_tail_ratio_max;
    const double convergence_estimator_tail_ratio_mean =
        stats != nullptr ? stats->convergence_estimator_tail_ratio_mean
                         : result.convergence_estimator_tail_ratio_mean;
    const double convergence_estimator_error_max =
        stats != nullptr ? stats->convergence_estimator_error_max
                         : result.convergence_estimator_error_max;
    const double convergence_estimator_error_mean =
        stats != nullptr ? stats->convergence_estimator_error_mean
                         : result.convergence_estimator_error_mean;
    const double convergence_estimator_abs_error_max =
        stats != nullptr ? stats->convergence_estimator_abs_error_max
                         : result.convergence_estimator_abs_error_max;
    const double convergence_estimator_abs_error_mean =
        stats != nullptr ? stats->convergence_estimator_abs_error_mean
                         : result.convergence_estimator_abs_error_mean;
    const double convergence_tail_norm_p_max =
        stats != nullptr ? stats->convergence_tail_norm_p_max
                         : result.convergence_tail_norm_p_max;
    const double convergence_tail_norm_p_mean =
        stats != nullptr ? stats->convergence_tail_norm_p_mean
                         : result.convergence_tail_norm_p_mean;
    const double convergence_tail_norm_pm1_max =
        stats != nullptr ? stats->convergence_tail_norm_pm1_max
                         : result.convergence_tail_norm_pm1_max;
    const double convergence_tail_norm_pm1_mean =
        stats != nullptr ? stats->convergence_tail_norm_pm1_mean
                         : result.convergence_tail_norm_pm1_mean;
    const double convergence_effective_jcp_amp_max =
        stats != nullptr ? stats->convergence_effective_jcp_amp_max
                         : result.convergence_effective_jcp_amp_max;
    const double convergence_effective_jcp_amp_mean =
        stats != nullptr ? stats->convergence_effective_jcp_amp_mean
                         : result.convergence_effective_jcp_amp_mean;
    const double convergence_dielectric_amp_max =
        stats != nullptr ? stats->convergence_dielectric_amp_max
                         : result.convergence_dielectric_amp_max;
    const double convergence_dielectric_amp_mean =
        stats != nullptr ? stats->convergence_dielectric_amp_mean
                         : result.convergence_dielectric_amp_mean;
    const auto convergence_channel_estimated_error_max =
        stats != nullptr ? stats->convergence_channel_estimated_error_max
                         : result.convergence_channel_estimated_error_max;
    const auto convergence_channel_effective_jcp_amp_max =
        stats != nullptr ? stats->convergence_channel_effective_jcp_amp_max
                         : result.convergence_channel_effective_jcp_amp_max;
    const auto convergence_channel_effective_jcp_amp_mean =
        stats != nullptr ? stats->convergence_channel_effective_jcp_amp_mean
                         : result.convergence_channel_effective_jcp_amp_mean;
    const auto convergence_channel_unstable_count =
        stats != nullptr ? stats->convergence_channel_unstable_far_count
                         : result.convergence_channel_unstable_far_count;
    const auto convergence_channel_unstable_panels =
        stats != nullptr ? stats->convergence_channel_unstable_far_panels
                         : result.convergence_channel_unstable_far_panels;
    const auto convergence_channel_split_rejects =
        stats != nullptr ? stats->convergence_channel_split_rejects
                         : result.convergence_channel_split_rejects;
    const auto convergence_channel_split_reject_panels =
        stats != nullptr ? stats->convergence_channel_split_reject_panels
                         : result.convergence_channel_split_reject_panels;
    const auto convergence_channel_low_abs_reaccepted_count =
        stats != nullptr ? stats->convergence_channel_low_abs_reaccepted_count
                         : result.convergence_channel_low_abs_reaccepted_count;
    const auto convergence_channel_low_abs_reaccepted_panels =
        stats != nullptr ? stats->convergence_channel_low_abs_reaccepted_panels
                         : result.convergence_channel_low_abs_reaccepted_panels;
    const auto convergence_channel_high_error_caught =
        stats != nullptr ? stats->convergence_channel_high_error_caught
                         : result.convergence_channel_high_error_caught;
    const auto convergence_channel_high_error_missed =
        stats != nullptr ? stats->convergence_channel_high_error_missed
                         : result.convergence_channel_high_error_missed;
    const auto convergence_channel_high_error_caught_by_split =
        stats != nullptr
            ? stats->convergence_channel_high_error_caught_by_split
            : result.convergence_channel_high_error_caught_by_split;
    const std::size_t convergence_estimator_high_error_caught =
        stats != nullptr ? stats->convergence_estimator_high_error_caught
                         : result.convergence_estimator_high_error_caught;
    const std::size_t convergence_estimator_high_error_missed =
        stats != nullptr ? stats->convergence_estimator_high_error_missed
                         : result.convergence_estimator_high_error_missed;
    const std::size_t convergence_estimator_high_error_caught_by_split =
        stats != nullptr
            ? stats->convergence_estimator_high_error_caught_by_split
            : result.convergence_estimator_high_error_caught_by_split;
    const std::size_t convergence_estimator_phantom_high_error_count =
        stats != nullptr
            ? stats->convergence_estimator_phantom_high_error_count
            : result.convergence_estimator_phantom_high_error_count;
    const double convergence_estimator_phantom_high_error_abs_max =
        stats != nullptr
            ? stats->convergence_estimator_phantom_high_error_abs_max
            : result.convergence_estimator_phantom_high_error_abs_max;
    const std::size_t convergence_overreject_count =
        stats != nullptr ? stats->convergence_overreject_count
                         : result.convergence_overreject_count;
    const std::size_t convergence_overreject_panels =
        stats != nullptr ? stats->convergence_overreject_panels
                         : result.convergence_overreject_panels;
    const std::size_t convergence_overreject_caught_count =
        stats != nullptr ? stats->convergence_overreject_caught_count
                         : result.convergence_overreject_caught_count;
    const std::size_t convergence_overreject_caught_panels =
        stats != nullptr ? stats->convergence_overreject_caught_panels
                         : result.convergence_overreject_caught_panels;
    const double convergence_overreject_abs_error_max =
        stats != nullptr ? stats->convergence_overreject_abs_error_max
                         : result.convergence_overreject_abs_error_max;
    const double convergence_overreject_abs_error_mean =
        stats != nullptr ? stats->convergence_overreject_abs_error_mean
                         : result.convergence_overreject_abs_error_mean;
    const double convergence_overreject_abs_error_p95 =
        stats != nullptr ? stats->convergence_overreject_abs_error_p95
                         : result.convergence_overreject_abs_error_p95;
    const double convergence_overreject_abs_error_p99 =
        stats != nullptr ? stats->convergence_overreject_abs_error_p99
                         : result.convergence_overreject_abs_error_p99;
    const double convergence_overreject_estimated_abs_error_max =
        stats != nullptr ? stats->convergence_overreject_estimated_abs_error_max
                         : result.convergence_overreject_estimated_abs_error_max;
    const double convergence_overreject_estimated_abs_error_mean =
        stats != nullptr ? stats->convergence_overreject_estimated_abs_error_mean
                         : result.convergence_overreject_estimated_abs_error_mean;
    const double convergence_overreject_required_scale_max =
        stats != nullptr ? stats->convergence_overreject_required_scale_max
                         : result.convergence_overreject_required_scale_max;
    const double convergence_overreject_required_scale_p95 =
        stats != nullptr ? stats->convergence_overreject_required_scale_p95
                         : result.convergence_overreject_required_scale_p95;
    const double convergence_overreject_required_scale_p99 =
        stats != nullptr ? stats->convergence_overreject_required_scale_p99
                         : result.convergence_overreject_required_scale_p99;
    const auto convergence_channel_overreject_count =
        stats != nullptr ? stats->convergence_channel_overreject_count
                         : result.convergence_channel_overreject_count;
    const auto convergence_channel_overreject_caught_count =
        stats != nullptr ? stats->convergence_channel_overreject_caught_count
                         : result.convergence_channel_overreject_caught_count;
    const auto convergence_channel_overreject_panels =
        stats != nullptr ? stats->convergence_channel_overreject_panels
                         : result.convergence_channel_overreject_panels;
    const auto convergence_channel_overreject_abs_error_max =
        stats != nullptr ? stats->convergence_channel_overreject_abs_error_max
                         : result.convergence_channel_overreject_abs_error_max;
    const auto convergence_channel_overreject_abs_error_mean =
        stats != nullptr ? stats->convergence_channel_overreject_abs_error_mean
                         : result.convergence_channel_overreject_abs_error_mean;
    const auto convergence_channel_overreject_required_scale_max =
        stats != nullptr ? stats->convergence_channel_overreject_required_scale_max
                         : result.convergence_channel_overreject_required_scale_max;
    const auto convergence_channel_overreject_required_scale_p99 =
        stats != nullptr ? stats->convergence_channel_overreject_required_scale_p99
                         : result.convergence_channel_overreject_required_scale_p99;
    const auto phantom_high_error_count =
        stats != nullptr ? stats->phantom_high_error_channel_count
                         : result.phantom_high_error_channel_count;
    const auto phantom_high_error_abs_max =
        stats != nullptr ? stats->phantom_high_error_channel_abs_max
                         : result.phantom_high_error_channel_abs_max;
    const auto accepted_miss_count =
        stats != nullptr ? stats->accepted_miss_channel_count
                         : result.accepted_miss_channel_count;
    const auto accepted_miss_panels =
        stats != nullptr ? stats->accepted_miss_channel_panels
                         : result.accepted_miss_channel_panels;
    const auto accepted_miss_rho_min =
        stats != nullptr ? stats->accepted_miss_channel_rho_min
                         : result.accepted_miss_channel_rho_min;
    const auto accepted_miss_rho_mean =
        stats != nullptr ? stats->accepted_miss_channel_rho_mean
                         : result.accepted_miss_channel_rho_mean;
    const auto accepted_miss_rho_max =
        stats != nullptr ? stats->accepted_miss_channel_rho_max
                         : result.accepted_miss_channel_rho_max;
    const auto accepted_miss_tail_ratio_mean =
        stats != nullptr ? stats->accepted_miss_channel_tail_ratio_mean
                         : result.accepted_miss_channel_tail_ratio_mean;
    const auto accepted_miss_tail_ratio_max =
        stats != nullptr ? stats->accepted_miss_channel_tail_ratio_max
                         : result.accepted_miss_channel_tail_ratio_max;
    const auto accepted_miss_estimated_error_mean =
        stats != nullptr
            ? stats->accepted_miss_channel_estimated_error_mean
            : result.accepted_miss_channel_estimated_error_mean;
    const auto accepted_miss_estimated_error_max =
        stats != nullptr
            ? stats->accepted_miss_channel_estimated_error_max
            : result.accepted_miss_channel_estimated_error_max;
    const std::string accepted_miss_worst_channel =
        stats != nullptr ? stats->accepted_miss_worst_channel
                         : result.accepted_miss_worst_channel;
    const std::size_t convergence_k3_rho_cap_rejects =
        stats != nullptr ? stats->convergence_k3_rho_cap_rejects
                         : result.convergence_k3_rho_cap_rejects;
    const std::size_t convergence_k3_rho_cap_reject_panels =
        stats != nullptr ? stats->convergence_k3_rho_cap_reject_panels
                         : result.convergence_k3_rho_cap_reject_panels;
    const std::size_t moment_nodes =
        stats != nullptr ? stats->moment_audit_nodes
                         : result.moment_audit_nodes;
    const double moment_rel_max =
        stats != nullptr ? stats->moment_translation_rel_error_max
                         : result.moment_translation_rel_error_max;
    const double moment_rel_mean =
        stats != nullptr ? stats->moment_translation_rel_error_mean
                         : result.moment_translation_rel_error_mean;
    const double moment_offset_max =
        stats != nullptr ? stats->moment_translation_offset_max
                         : result.moment_translation_offset_max;
    const double moment_offset_mean =
        stats != nullptr ? stats->moment_translation_offset_mean
                         : result.moment_translation_offset_mean;
    const double moment_cancel_max =
        stats != nullptr ? stats->moment_cancellation_ratio_max
                         : result.moment_cancellation_ratio_max;
    const double moment_cancel_mean =
        stats != nullptr ? stats->moment_cancellation_ratio_mean
                         : result.moment_cancellation_ratio_mean;
    const double moment_dynamic_range =
        stats != nullptr ? stats->moment_dynamic_range_max
                         : result.moment_dynamic_range_max;
    const std::size_t expansion_pairs =
        stats != nullptr ? stats->far_expansion_audit_pairs
                         : result.far_expansion_audit_pairs;
    const std::size_t expansion_panels =
        stats != nullptr ? stats->far_expansion_audit_panels
                         : result.far_expansion_audit_panels;
    const std::string expansion_audit_channel =
        stats != nullptr ? stats->expansion_audit_channel
                         : result.expansion_audit_channel;
    const double expansion_direct_ms =
        stats != nullptr ? stats->far_expansion_direct_ms
                         : result.far_expansion_direct_ms;
    const double expansion_rel_max =
        stats != nullptr ? stats->far_expansion_rel_error_max
                         : result.far_expansion_rel_error_max;
    const double expansion_rel_mean =
        stats != nullptr ? stats->far_expansion_rel_error_mean
                         : result.far_expansion_rel_error_mean;
    const double expansion_rel_p95 =
        stats != nullptr ? stats->far_expansion_rel_error_p95
                         : result.far_expansion_rel_error_p95;
    const double expansion_rel_p99 =
        stats != nullptr ? stats->far_expansion_rel_error_p99
                         : result.far_expansion_rel_error_p99;
    const std::string expansion_worst_channel =
        stats != nullptr ? stats->far_expansion_worst_channel
                         : result.far_expansion_worst_channel;
    const auto expansion_channel_rel_max =
        stats != nullptr ? stats->far_expansion_channel_rel_error_max
                         : result.far_expansion_channel_rel_error_max;
    const auto expansion_channel_rel_mean =
        stats != nullptr ? stats->far_expansion_channel_rel_error_mean
                         : result.far_expansion_channel_rel_error_mean;
    const auto expansion_channel_rel_p99 =
        stats != nullptr ? stats->far_expansion_channel_rel_error_p99
                         : result.far_expansion_channel_rel_error_p99;
    const bool expansion_mean_gate_pass =
        stats != nullptr ? stats->far_expansion_mean_gate_pass
                         : result.far_expansion_mean_gate_pass;
    const double convergence_rho_max =
        stats != nullptr ? stats->far_convergence_rho_max
                         : result.far_convergence_rho_max;
    const double convergence_rho_mean =
        stats != nullptr ? stats->far_convergence_rho_mean
                         : result.far_convergence_rho_mean;
    const double node_aspect_max =
        stats != nullptr ? stats->node_aspect_ratio_max
                         : result.node_aspect_ratio_max;
    const double node_aspect_mean =
        stats != nullptr ? stats->node_aspect_ratio_mean
                         : result.node_aspect_ratio_mean;
    const double accepted_aspect_max =
        stats != nullptr ? stats->accepted_far_aspect_ratio_max
                         : result.accepted_far_aspect_ratio_max;
    const double accepted_aspect_mean =
        stats != nullptr ? stats->accepted_far_aspect_ratio_mean
                         : result.accepted_far_aspect_ratio_mean;
    const std::size_t mac_convergence_rho_rejects =
        stats != nullptr ? stats->mac_convergence_rho_rejects
                         : result.mac_convergence_rho_rejects;
    const std::size_t mac_convergence_aspect_rejects =
        stats != nullptr ? stats->mac_convergence_aspect_rejects
                         : result.mac_convergence_aspect_rejects;
    const std::size_t mac_convergence_rho_reject_panels =
        stats != nullptr ? stats->mac_convergence_rho_reject_panels
                         : result.mac_convergence_rho_reject_panels;
    const std::size_t mac_convergence_aspect_reject_panels =
        stats != nullptr ? stats->mac_convergence_aspect_reject_panels
                         : result.mac_convergence_aspect_reject_panels;
    const std::size_t mac_surface_margin_rejects =
        stats != nullptr ? stats->mac_surface_margin_rejects
                         : result.mac_surface_margin_rejects;
    const std::size_t mac_surface_margin_reject_panels =
        stats != nullptr ? stats->mac_surface_margin_reject_panels
                         : result.mac_surface_margin_reject_panels;
    const std::size_t mac_surface_margin_applied =
        stats != nullptr ? stats->mac_surface_margin_applied
                         : result.mac_surface_margin_applied;
    const std::size_t mac_surface_margin_skipped =
        stats != nullptr ? stats->mac_surface_margin_skipped
                         : result.mac_surface_margin_skipped;
    const std::size_t mac_surface_margin_adaptive_rejects =
        stats != nullptr ? stats->mac_surface_margin_adaptive_rejects
                         : result.mac_surface_margin_adaptive_rejects;
    const std::size_t mac_surface_margin_adaptive_reject_panels =
        stats != nullptr ? stats->mac_surface_margin_adaptive_reject_panels
                         : result.mac_surface_margin_adaptive_reject_panels;
    const std::size_t mac_surface_margin_skipped_accepted_far =
        stats != nullptr ? stats->mac_surface_margin_skipped_accepted_far
                         : result.mac_surface_margin_skipped_accepted_far;
    const std::size_t mac_surface_margin_hard_applied =
        stats != nullptr ? stats->mac_surface_margin_hard_applied
                         : result.mac_surface_margin_hard_applied;
    const std::size_t mac_surface_margin_soft_applied =
        stats != nullptr ? stats->mac_surface_margin_soft_applied
                         : result.mac_surface_margin_soft_applied;
    const std::size_t mac_surface_margin_hard_rejects =
        stats != nullptr ? stats->mac_surface_margin_hard_rejects
                         : result.mac_surface_margin_hard_rejects;
    const std::size_t mac_surface_margin_hard_reject_panels =
        stats != nullptr ? stats->mac_surface_margin_hard_reject_panels
                         : result.mac_surface_margin_hard_reject_panels;
    const std::size_t mac_surface_margin_soft_rejects =
        stats != nullptr ? stats->mac_surface_margin_soft_rejects
                         : result.mac_surface_margin_soft_rejects;
    const std::size_t mac_surface_margin_soft_reject_panels =
        stats != nullptr ? stats->mac_surface_margin_soft_reject_panels
                         : result.mac_surface_margin_soft_reject_panels;
    const std::string leaf_work_profile =
        stats != nullptr ? stats->leaf_work_profile : result.leaf_work_profile;
    const std::size_t leaf_work_mac_terminal_panels =
        stats != nullptr ? stats->leaf_work_mac_terminal_panels
                         : result.leaf_work_mac_terminal_panels;
    const std::size_t leaf_work_physical_gap_panels =
        stats != nullptr ? stats->leaf_work_physical_gap_panels
                         : result.leaf_work_physical_gap_panels;
    const std::size_t leaf_work_fragmentation_panels =
        stats != nullptr ? stats->leaf_work_fragmentation_panels
                         : result.leaf_work_fragmentation_panels;
    const std::size_t leaf_work_convergence_panels =
        stats != nullptr ? stats->leaf_work_convergence_panels
                         : result.leaf_work_convergence_panels;
    const std::size_t leaf_work_convergence_estimator_panels =
        stats != nullptr ? stats->leaf_work_convergence_estimator_panels
                         : result.leaf_work_convergence_estimator_panels;
    const std::size_t leaf_work_convergence_rho_panels =
        stats != nullptr ? stats->leaf_work_convergence_rho_panels
                         : result.leaf_work_convergence_rho_panels;
    const std::size_t leaf_work_convergence_aspect_panels =
        stats != nullptr ? stats->leaf_work_convergence_aspect_panels
                         : result.leaf_work_convergence_aspect_panels;
    const std::size_t leaf_work_surface_margin_panels =
        stats != nullptr ? stats->leaf_work_surface_margin_panels
                         : result.leaf_work_surface_margin_panels;
    const std::size_t leaf_work_surface_margin_hard_panels =
        stats != nullptr ? stats->leaf_work_surface_margin_hard_panels
                         : result.leaf_work_surface_margin_hard_panels;
    const std::size_t leaf_work_surface_margin_soft_panels =
        stats != nullptr ? stats->leaf_work_surface_margin_soft_panels
                         : result.leaf_work_surface_margin_soft_panels;
    const std::size_t leaf_work_k3_rho_cap_panels =
        stats != nullptr ? stats->leaf_work_k3_rho_cap_panels
                         : result.leaf_work_k3_rho_cap_panels;
    const std::size_t leaf_work_singular_fallback_panels =
        stats != nullptr ? stats->leaf_work_singular_fallback_panels
                         : result.leaf_work_singular_fallback_panels;
    const std::size_t leaf_work_scalar_tail_panels =
        stats != nullptr ? stats->leaf_work_scalar_tail_panels
                         : result.leaf_work_scalar_tail_panels;
    const std::size_t leaf_work_profile_accounted_panels =
        stats != nullptr ? stats->leaf_work_profile_accounted_panels
                         : result.leaf_work_profile_accounted_panels;
    const double leaf_work_mac_terminal_frac =
        stats != nullptr ? stats->leaf_work_mac_terminal_frac
                         : result.leaf_work_mac_terminal_frac;
    const double leaf_work_convergence_frac =
        stats != nullptr ? stats->leaf_work_convergence_frac
                         : result.leaf_work_convergence_frac;
    const double leaf_work_surface_margin_frac =
        stats != nullptr ? stats->leaf_work_surface_margin_frac
                         : result.leaf_work_surface_margin_frac;
    const double leaf_work_singular_fallback_frac =
        stats != nullptr ? stats->leaf_work_singular_fallback_frac
                         : result.leaf_work_singular_fallback_frac;
    const double leaf_work_scalar_tail_frac =
        stats != nullptr ? stats->leaf_work_scalar_tail_frac
                         : result.leaf_work_scalar_tail_frac;
    const std::string cluster_policy =
        stats != nullptr ? stats->cluster_policy : result.cluster_policy;
    const int target_batch_size =
        stats != nullptr ? stats->target_batch_size : result.target_batch_size;
    const int leaf_pair_min_source_panels =
        stats != nullptr ? stats->leaf_pair_min_source_panels
                         : result.leaf_pair_min_source_panels;
    const int leaf_pair_min_targets =
        stats != nullptr ? stats->leaf_pair_min_targets
                         : result.leaf_pair_min_targets;
    const std::size_t leaf_pair_batch_interactions =
        stats != nullptr ? stats->leaf_pair_batch_interactions
                         : result.leaf_pair_batch_interactions;
    const std::size_t leaf_pair_batch_targets =
        stats != nullptr ? stats->leaf_pair_batch_targets
                         : result.leaf_pair_batch_targets;
    const std::size_t leaf_pair_batch_source_panels =
        stats != nullptr ? stats->leaf_pair_batch_source_panels
                         : result.leaf_pair_batch_source_panels;
    const std::size_t leaf_pair_batch_residual_panels =
        stats != nullptr ? stats->leaf_pair_batch_residual_panels
                         : result.leaf_pair_batch_residual_panels;
    const std::size_t leaf_pair_batch_fallbacks =
        stats != nullptr ? stats->leaf_pair_batch_fallbacks
                         : result.leaf_pair_batch_fallbacks;
    const double leaf_pair_batch_ms =
        stats != nullptr ? stats->leaf_pair_batch_ms
                         : result.leaf_pair_batch_ms;
    const std::string leaf_pair_batch_kernel =
        stats != nullptr ? stats->leaf_pair_batch_kernel
                         : result.leaf_pair_batch_kernel;
    const std::size_t leaf_pair_batch_vectorized_sources =
        stats != nullptr ? stats->leaf_pair_batch_vectorized_sources
                         : result.leaf_pair_batch_vectorized_sources;
    const std::size_t leaf_pair_batch_vector_lanes =
        stats != nullptr ? stats->leaf_pair_batch_vector_lanes
                         : result.leaf_pair_batch_vector_lanes;
    const double leaf_pair_batch_pairs_per_second =
        stats != nullptr ? stats->leaf_pair_batch_pairs_per_second
                         : result.leaf_pair_batch_pairs_per_second;
    const std::string leaf_kernel_profile =
        stats != nullptr ? stats->leaf_kernel_profile
                         : result.leaf_kernel_profile;
    const std::string jcp_direct_kernel =
        stats != nullptr ? stats->jcp_direct_kernel
                         : result.jcp_direct_kernel;
    const std::size_t leaf_kernel_scalar_tasks =
        stats != nullptr ? stats->leaf_kernel_scalar_tasks
                         : result.leaf_kernel_scalar_tasks;
    const std::size_t leaf_kernel_batch_tasks =
        stats != nullptr ? stats->leaf_kernel_batch_tasks
                         : result.leaf_kernel_batch_tasks;
    const std::size_t leaf_kernel_source_panel_evals =
        stats != nullptr ? stats->leaf_kernel_source_panel_evals
                         : result.leaf_kernel_source_panel_evals;
    const double leaf_kernel_scalar_ms =
        stats != nullptr ? stats->leaf_kernel_scalar_ms
                         : result.leaf_kernel_scalar_ms;
    const double leaf_kernel_target_avx2_ms =
        stats != nullptr ? stats->leaf_kernel_target_avx2_ms
                         : result.leaf_kernel_target_avx2_ms;
    const double leaf_kernel_tail_ms =
        stats != nullptr ? stats->leaf_kernel_tail_ms
                         : result.leaf_kernel_tail_ms;
    const std::size_t leaf_pair_batch_soa_sources =
        stats != nullptr ? stats->leaf_pair_batch_soa_sources
                         : result.leaf_pair_batch_soa_sources;
    const double leaf_pair_batch_soa_ms =
        stats != nullptr ? stats->leaf_pair_batch_soa_ms
                         : result.leaf_pair_batch_soa_ms;
    const double leaf_pair_batch_scalar_tail_ms =
        stats != nullptr ? stats->leaf_pair_batch_scalar_tail_ms
                         : result.leaf_pair_batch_scalar_tail_ms;
    const double leaf_pair_batch_clean_source_panels_per_second =
        stats != nullptr
            ? stats->leaf_pair_batch_clean_source_panels_per_second
            : result.leaf_pair_batch_clean_source_panels_per_second;
    const std::string leaf_pair_selection_profile =
        stats != nullptr ? stats->leaf_pair_selection_profile
                         : result.leaf_pair_selection_profile;
    const std::string leaf_pair_probe_cache =
        stats != nullptr ? stats->leaf_pair_probe_cache
                         : result.leaf_pair_probe_cache;
    const std::string leaf_pair_target_group =
        stats != nullptr ? stats->leaf_pair_target_group
                         : result.leaf_pair_target_group;
    const std::string leaf_pair_directness_probe =
        stats != nullptr ? stats->leaf_pair_directness_probe
                         : result.leaf_pair_directness_probe;
    const std::string leaf_pair_directness_cache =
        stats != nullptr ? stats->leaf_pair_directness_cache
                         : result.leaf_pair_directness_cache;
    const std::string leaf_pair_directness_verify =
        stats != nullptr ? stats->leaf_pair_directness_verify
                         : result.leaf_pair_directness_verify;
    const int leaf_pair_max_lookahead =
        stats != nullptr ? stats->leaf_pair_max_lookahead
                         : result.leaf_pair_max_lookahead;
    const std::size_t leaf_pair_selection_candidate_scans =
        stats != nullptr ? stats->leaf_pair_selection_candidate_scans
                         : result.leaf_pair_selection_candidate_scans;
    const std::size_t leaf_pair_selection_clean_rejects =
        stats != nullptr ? stats->leaf_pair_selection_clean_rejects
                         : result.leaf_pair_selection_clean_rejects;
    const std::size_t leaf_pair_selection_probe_calls =
        stats != nullptr ? stats->leaf_pair_selection_probe_calls
                         : result.leaf_pair_selection_probe_calls;
    const std::size_t leaf_pair_selection_probe_node_visits =
        stats != nullptr ? stats->leaf_pair_selection_probe_node_visits
                         : result.leaf_pair_selection_probe_node_visits;
    const double leaf_pair_selection_probe_ms =
        stats != nullptr ? stats->leaf_pair_selection_probe_ms
                         : result.leaf_pair_selection_probe_ms;
    const double leaf_pair_selection_ms =
        stats != nullptr ? stats->leaf_pair_selection_ms
                         : result.leaf_pair_selection_ms;
    const std::size_t leaf_pair_selection_probe_cache_hits =
        stats != nullptr ? stats->leaf_pair_selection_probe_cache_hits
                         : result.leaf_pair_selection_probe_cache_hits;
    const std::size_t leaf_pair_selection_probe_cache_misses =
        stats != nullptr ? stats->leaf_pair_selection_probe_cache_misses
                         : result.leaf_pair_selection_probe_cache_misses;
    const std::size_t leaf_pair_selection_insufficient_targets =
        stats != nullptr ? stats->leaf_pair_selection_insufficient_targets
                         : result.leaf_pair_selection_insufficient_targets;
    const std::size_t leaf_pair_selection_same_target_leaf_rejects =
        stats != nullptr ? stats->leaf_pair_selection_same_target_leaf_rejects
                         : result.leaf_pair_selection_same_target_leaf_rejects;
    const std::size_t leaf_pair_directness_ancestor_probe_calls =
        stats != nullptr
            ? stats->leaf_pair_directness_ancestor_probe_calls
            : result.leaf_pair_directness_ancestor_probe_calls;
    const std::size_t leaf_pair_directness_ancestor_probe_node_visits =
        stats != nullptr
            ? stats->leaf_pair_directness_ancestor_probe_node_visits
            : result.leaf_pair_directness_ancestor_probe_node_visits;
    const double leaf_pair_directness_ancestor_probe_ms =
        stats != nullptr ? stats->leaf_pair_directness_ancestor_probe_ms
                         : result.leaf_pair_directness_ancestor_probe_ms;
    const std::size_t leaf_pair_directness_pair_cache_hits =
        stats != nullptr ? stats->leaf_pair_directness_pair_cache_hits
                         : result.leaf_pair_directness_pair_cache_hits;
    const std::size_t leaf_pair_directness_pair_cache_misses =
        stats != nullptr ? stats->leaf_pair_directness_pair_cache_misses
                         : result.leaf_pair_directness_pair_cache_misses;
    const std::size_t leaf_pair_directness_all_direct_proofs =
        stats != nullptr ? stats->leaf_pair_directness_all_direct_proofs
                         : result.leaf_pair_directness_all_direct_proofs;
    const std::size_t leaf_pair_directness_mixed_fallbacks =
        stats != nullptr ? stats->leaf_pair_directness_mixed_fallbacks
                         : result.leaf_pair_directness_mixed_fallbacks;
    const std::size_t leaf_pair_directness_verify_checks =
        stats != nullptr ? stats->leaf_pair_directness_verify_checks
                         : result.leaf_pair_directness_verify_checks;
    const std::size_t leaf_pair_directness_verify_mismatches =
        stats != nullptr ? stats->leaf_pair_directness_verify_mismatches
                         : result.leaf_pair_directness_verify_mismatches;
    const std::string leaf_pair_scheduler =
        stats != nullptr ? stats->leaf_pair_scheduler
                         : result.leaf_pair_scheduler;
    const int leaf_pair_worklist_max_active =
        stats != nullptr ? stats->leaf_pair_worklist_max_active
                         : result.leaf_pair_worklist_max_active;
    const std::string leaf_pair_worklist_verify =
        stats != nullptr ? stats->leaf_pair_worklist_verify
                         : result.leaf_pair_worklist_verify;
    const std::size_t leaf_pair_worklist_tasks =
        stats != nullptr ? stats->leaf_pair_worklist_tasks
                         : result.leaf_pair_worklist_tasks;
    const std::size_t leaf_pair_worklist_groups =
        stats != nullptr ? stats->leaf_pair_worklist_groups
                         : result.leaf_pair_worklist_groups;
    const std::size_t leaf_pair_worklist_batched_groups =
        stats != nullptr ? stats->leaf_pair_worklist_batched_groups
                         : result.leaf_pair_worklist_batched_groups;
    const std::size_t leaf_pair_worklist_scalar_flush_tasks =
        stats != nullptr ? stats->leaf_pair_worklist_scalar_flush_tasks
                         : result.leaf_pair_worklist_scalar_flush_tasks;
    const std::size_t leaf_pair_worklist_forced_flushes =
        stats != nullptr ? stats->leaf_pair_worklist_forced_flushes
                         : result.leaf_pair_worklist_forced_flushes;
    const std::size_t leaf_pair_worklist_active_groups_max =
        stats != nullptr ? stats->leaf_pair_worklist_active_groups_max
                         : result.leaf_pair_worklist_active_groups_max;
    const double leaf_pair_worklist_ms =
        stats != nullptr ? stats->leaf_pair_worklist_ms
                         : result.leaf_pair_worklist_ms;
    const std::size_t leaf_pair_worklist_verify_checks =
        stats != nullptr ? stats->leaf_pair_worklist_verify_checks
                         : result.leaf_pair_worklist_verify_checks;
    const std::size_t leaf_pair_worklist_verify_mismatches =
        stats != nullptr ? stats->leaf_pair_worklist_verify_mismatches
                         : result.leaf_pair_worklist_verify_mismatches;
    const std::string leaf_block_verify =
        stats != nullptr ? stats->leaf_block_verify
                         : result.leaf_block_verify;
    const std::size_t leaf_block_target_leaves =
        stats != nullptr ? stats->leaf_block_target_leaves
                         : result.leaf_block_target_leaves;
    const std::size_t leaf_block_targets =
        stats != nullptr ? stats->leaf_block_targets
                         : result.leaf_block_targets;
    const std::size_t leaf_block_source_groups =
        stats != nullptr ? stats->leaf_block_source_groups
                         : result.leaf_block_source_groups;
    const std::size_t leaf_block_batched_groups =
        stats != nullptr ? stats->leaf_block_batched_groups
                         : result.leaf_block_batched_groups;
    const std::size_t leaf_block_scalar_tail_tasks =
        stats != nullptr ? stats->leaf_block_scalar_tail_tasks
                         : result.leaf_block_scalar_tail_tasks;
    const double leaf_block_group_ms =
        stats != nullptr ? stats->leaf_block_group_ms
                         : result.leaf_block_group_ms;
    const double leaf_block_flush_ms =
        stats != nullptr ? stats->leaf_block_flush_ms
                         : result.leaf_block_flush_ms;
    const double leaf_block_writeback_ms =
        stats != nullptr ? stats->leaf_block_writeback_ms
                         : result.leaf_block_writeback_ms;
    const std::size_t leaf_block_verify_checks =
        stats != nullptr ? stats->leaf_block_verify_checks
                         : result.leaf_block_verify_checks;
    const std::size_t leaf_block_verify_mismatches =
        stats != nullptr ? stats->leaf_block_verify_mismatches
                         : result.leaf_block_verify_mismatches;

    out_ << std::fixed << std::setprecision(6)
         << std::chrono::duration<double>(now - start_).count() << ','
         << csv_escape_progress(phase) << ',' << gmres_iter << ','
         << restart_iter << ',' << matvec_calls << ',' << std::scientific
         << std::setprecision(12) << residual << ','
         << result.gmres_restart_requested << ','
         << result.gmres_restart_effective << ','
         << result.best_residual << ',' << result.best_residual_iter << ','
         << result.final_residual << ','
         << result.best_preconditioned_residual << ','
         << result.final_preconditioned_residual << ','
         << result.true_residual_abs << ',' << result.true_residual_rel
         << ',' << csv_escape_progress(result.true_residual_probe_mode)
         << ',' << result.true_residual_probe_interval << ','
         << result.true_residual_probe_calls << ',' << std::fixed
         << std::setprecision(6) << result.true_residual_probe_ms / 1000.0
         << ',' << csv_escape_progress(result.true_residual_probe_scope)
         << ',' << result.true_residual_sample_count << ','
         << result.sampled_true_residual_rel << ','
         << result.sampled_true_residual_rows << ','
         << result.sampled_probe_calls << ',' << std::fixed
         << std::setprecision(6) << result.sampled_probe_ms / 1000.0 << ','
         << result.full_true_residual_confirm_calls << ',' << std::fixed
         << std::setprecision(6)
         << result.full_true_residual_confirm_ms / 1000.0 << std::fixed
         << std::setprecision(6) << ',' << result.panel_count << ','
         << result.unknown_count << ',' << result.taylor_order << ','
         << tree_nodes << ','
         << csv_escape_progress(interaction_mode) << ','
         << csv_escape_progress(tree_builder) << ','
         << csv_escape_progress(node_bounds) << ','
         << csv_escape_progress(mac_safety) << ','
         << mac_gap_factor << ','
         << csv_escape_progress(expansion_center) << ','
         << csv_escape_progress(fragmentation_safety) << ','
         << csv_escape_progress(moment_construction) << ','
         << csv_escape_progress(moment_sum) << ','
         << csv_escape_progress(convergence_audit) << ','
         << csv_escape_progress(convergence_safety) << ','
         << convergence_rho_cap << ',' << mac_aspect_cap << ','
         << csv_escape_progress(convergence_estimator) << ','
         << convergence_tail_threshold << ','
         << convergence_estimator_error_cap << ','
         << convergence_unstable_rho_cap << ','
         << csv_escape_progress(convergence_calibration) << ','
         << csv_escape_progress(convergence_dielectric_amp) << ','
         << csv_escape_progress(convergence_estimator_model) << ','
         << csv_escape_progress(convergence_channel_calibration) << ','
         << csv_escape_progress(convergence_split_policy) << ','
         << csv_escape_progress(convergence_split_selectivity) << ','
         << convergence_split_absolute_floor << ','
         << csv_escape_progress(convergence_overreject_audit) << ','
         << csv_escape_progress(convergence_overreject_detail) << ','
         << csv_escape_progress(convergence_shape_audit) << ','
         << csv_escape_progress(convergence_kernel_response_audit) << ','
         << csv_escape_progress(channel_missed_semantics) << ','
         << csv_escape_progress(accepted_miss_diagnostics) << ','
         << accepted_miss_absolute_floor << ','
         << convergence_k3_accepted_rho_cap << ','
         << csv_escape_progress(mac_distance_policy) << ','
         << mac_surface_margin_factor << ','
         << csv_escape_progress(mac_surface_margin_policy) << ','
         << mac_surface_margin_aspect_threshold << ','
         << mac_surface_margin_rho_trigger << ','
         << mac_surface_margin_soft_aspect_threshold << ','
         << mac_surface_margin_soft_rho_trigger << ','
         << mac_surface_margin_soft_factor << ','
         << convergence_estimator_samples << ','
         << convergence_estimator_unstable_count << ','
         << convergence_estimator_unstable_panels << ','
         << convergence_estimator_split_rejects << ','
         << convergence_estimator_split_reject_panels << ','
         << convergence_low_abs_reaccepted_count << ','
         << convergence_low_abs_reaccepted_panels << ','
         << convergence_estimator_tail_ratio_max << ','
         << convergence_estimator_tail_ratio_mean << ','
         << convergence_estimator_error_max << ','
         << convergence_estimator_error_mean << ','
         << convergence_estimator_abs_error_max << ','
         << convergence_estimator_abs_error_mean << ','
         << convergence_tail_norm_p_max << ','
         << convergence_tail_norm_p_mean << ','
         << convergence_tail_norm_pm1_max << ','
         << convergence_tail_norm_pm1_mean << ','
         << convergence_effective_jcp_amp_max << ','
         << convergence_effective_jcp_amp_mean << ','
         << convergence_dielectric_amp_max << ','
         << convergence_dielectric_amp_mean << ','
         << convergence_channel_effective_jcp_amp_max[0] << ','
         << convergence_channel_effective_jcp_amp_max[1] << ','
         << convergence_channel_effective_jcp_amp_max[2] << ','
         << convergence_channel_effective_jcp_amp_max[3] << ','
         << convergence_channel_effective_jcp_amp_mean[0] << ','
         << convergence_channel_effective_jcp_amp_mean[1] << ','
         << convergence_channel_effective_jcp_amp_mean[2] << ','
         << convergence_channel_effective_jcp_amp_mean[3] << ','
         << convergence_channel_estimated_error_max[0] << ','
         << convergence_channel_estimated_error_max[1] << ','
         << convergence_channel_estimated_error_max[2] << ','
         << convergence_channel_estimated_error_max[3] << ','
         << convergence_channel_unstable_count[0] << ','
         << convergence_channel_unstable_count[1] << ','
         << convergence_channel_unstable_count[2] << ','
         << convergence_channel_unstable_count[3] << ','
         << convergence_channel_unstable_panels[0] << ','
         << convergence_channel_unstable_panels[1] << ','
         << convergence_channel_unstable_panels[2] << ','
         << convergence_channel_unstable_panels[3] << ','
         << convergence_channel_split_rejects[0] << ','
         << convergence_channel_split_rejects[1] << ','
         << convergence_channel_split_rejects[2] << ','
         << convergence_channel_split_rejects[3] << ','
         << convergence_channel_split_reject_panels[0] << ','
         << convergence_channel_split_reject_panels[1] << ','
         << convergence_channel_split_reject_panels[2] << ','
         << convergence_channel_split_reject_panels[3] << ','
         << convergence_channel_low_abs_reaccepted_count[0] << ','
         << convergence_channel_low_abs_reaccepted_count[1] << ','
         << convergence_channel_low_abs_reaccepted_count[2] << ','
         << convergence_channel_low_abs_reaccepted_count[3] << ','
         << convergence_channel_low_abs_reaccepted_panels[0] << ','
         << convergence_channel_low_abs_reaccepted_panels[1] << ','
         << convergence_channel_low_abs_reaccepted_panels[2] << ','
         << convergence_channel_low_abs_reaccepted_panels[3] << ','
         << convergence_channel_high_error_caught[0] << ','
         << convergence_channel_high_error_caught[1] << ','
         << convergence_channel_high_error_caught[2] << ','
         << convergence_channel_high_error_caught[3] << ','
         << convergence_channel_high_error_missed[0] << ','
         << convergence_channel_high_error_missed[1] << ','
         << convergence_channel_high_error_missed[2] << ','
         << convergence_channel_high_error_missed[3] << ','
         << convergence_channel_high_error_caught_by_split[0] << ','
         << convergence_channel_high_error_caught_by_split[1] << ','
         << convergence_channel_high_error_caught_by_split[2] << ','
         << convergence_channel_high_error_caught_by_split[3] << ','
         << convergence_estimator_high_error_caught << ','
         << convergence_estimator_high_error_missed << ','
         << convergence_estimator_high_error_caught_by_split << ','
         << convergence_estimator_phantom_high_error_count << ','
         << convergence_estimator_phantom_high_error_abs_max << ','
         << convergence_overreject_count << ','
         << convergence_overreject_panels << ','
         << convergence_overreject_caught_count << ','
         << convergence_overreject_caught_panels << ','
         << convergence_overreject_abs_error_max << ','
         << convergence_overreject_abs_error_mean << ','
         << convergence_overreject_abs_error_p95 << ','
         << convergence_overreject_abs_error_p99 << ','
         << convergence_overreject_estimated_abs_error_max << ','
         << convergence_overreject_estimated_abs_error_mean << ','
         << convergence_overreject_required_scale_max << ','
         << convergence_overreject_required_scale_p95 << ','
         << convergence_overreject_required_scale_p99 << ','
         << convergence_channel_overreject_count[0] << ','
         << convergence_channel_overreject_count[1] << ','
         << convergence_channel_overreject_count[2] << ','
         << convergence_channel_overreject_count[3] << ','
         << convergence_channel_overreject_caught_count[0] << ','
         << convergence_channel_overreject_caught_count[1] << ','
         << convergence_channel_overreject_caught_count[2] << ','
         << convergence_channel_overreject_caught_count[3] << ','
         << convergence_channel_overreject_panels[0] << ','
         << convergence_channel_overreject_panels[1] << ','
         << convergence_channel_overreject_panels[2] << ','
         << convergence_channel_overreject_panels[3] << ','
         << convergence_channel_overreject_abs_error_max[0] << ','
         << convergence_channel_overreject_abs_error_max[1] << ','
         << convergence_channel_overreject_abs_error_max[2] << ','
         << convergence_channel_overreject_abs_error_max[3] << ','
         << convergence_channel_overreject_abs_error_mean[0] << ','
         << convergence_channel_overreject_abs_error_mean[1] << ','
         << convergence_channel_overreject_abs_error_mean[2] << ','
         << convergence_channel_overreject_abs_error_mean[3] << ','
         << convergence_channel_overreject_required_scale_max[0] << ','
         << convergence_channel_overreject_required_scale_max[1] << ','
         << convergence_channel_overreject_required_scale_max[2] << ','
         << convergence_channel_overreject_required_scale_max[3] << ','
         << convergence_channel_overreject_required_scale_p99[0] << ','
         << convergence_channel_overreject_required_scale_p99[1] << ','
         << convergence_channel_overreject_required_scale_p99[2] << ','
         << convergence_channel_overreject_required_scale_p99[3] << ','
         << phantom_high_error_count[0] << ','
         << phantom_high_error_count[1] << ','
         << phantom_high_error_count[2] << ','
         << phantom_high_error_count[3] << ','
         << phantom_high_error_abs_max[0] << ','
         << phantom_high_error_abs_max[1] << ','
         << phantom_high_error_abs_max[2] << ','
         << phantom_high_error_abs_max[3] << ','
         << csv_escape_progress(accepted_miss_worst_channel) << ','
         << accepted_miss_count[0] << ',' << accepted_miss_count[1] << ','
         << accepted_miss_count[2] << ',' << accepted_miss_count[3] << ','
         << accepted_miss_panels[0] << ',' << accepted_miss_panels[1] << ','
         << accepted_miss_panels[2] << ',' << accepted_miss_panels[3] << ','
         << accepted_miss_rho_min[0] << ',' << accepted_miss_rho_min[1]
         << ',' << accepted_miss_rho_min[2] << ','
         << accepted_miss_rho_min[3] << ',' << accepted_miss_rho_mean[0]
         << ',' << accepted_miss_rho_mean[1] << ','
         << accepted_miss_rho_mean[2] << ','
         << accepted_miss_rho_mean[3] << ',' << accepted_miss_rho_max[0]
         << ',' << accepted_miss_rho_max[1] << ','
         << accepted_miss_rho_max[2] << ',' << accepted_miss_rho_max[3]
         << ',' << accepted_miss_tail_ratio_mean[0] << ','
         << accepted_miss_tail_ratio_mean[1] << ','
         << accepted_miss_tail_ratio_mean[2] << ','
         << accepted_miss_tail_ratio_mean[3] << ','
         << accepted_miss_tail_ratio_max[0] << ','
         << accepted_miss_tail_ratio_max[1] << ','
         << accepted_miss_tail_ratio_max[2] << ','
         << accepted_miss_tail_ratio_max[3] << ','
         << accepted_miss_estimated_error_mean[0] << ','
         << accepted_miss_estimated_error_mean[1] << ','
         << accepted_miss_estimated_error_mean[2] << ','
         << accepted_miss_estimated_error_mean[3] << ','
         << accepted_miss_estimated_error_max[0] << ','
         << accepted_miss_estimated_error_max[1] << ','
         << accepted_miss_estimated_error_max[2] << ','
         << accepted_miss_estimated_error_max[3] << ','
         << convergence_k3_rho_cap_rejects << ','
         << convergence_k3_rho_cap_reject_panels << ','
         << csv_escape_progress(expansion_audit) << ','
         << csv_escape_progress(expansion_audit_channel) << ','
         << unsafe_count << ','
         << unsafe_panels << ',' << gap_samples << ',' << gap_rejects << ','
         << mac_convergence_rho_rejects << ','
         << mac_convergence_aspect_rejects << ','
         << mac_convergence_rho_reject_panels << ','
         << mac_convergence_aspect_reject_panels << ','
         << mac_surface_margin_rejects << ','
         << mac_surface_margin_reject_panels << ','
         << mac_surface_margin_applied << ','
         << mac_surface_margin_skipped << ','
         << mac_surface_margin_adaptive_rejects << ','
         << mac_surface_margin_adaptive_reject_panels << ','
         << mac_surface_margin_skipped_accepted_far << ','
         << mac_surface_margin_hard_applied << ','
         << mac_surface_margin_soft_applied << ','
         << mac_surface_margin_hard_rejects << ','
         << mac_surface_margin_hard_reject_panels << ','
         << mac_surface_margin_soft_rejects << ','
         << mac_surface_margin_soft_reject_panels << ','
         << fragmented_count << ',' << fragmented_panels << ','
         << fragmented_rejects << ',' << fragmentation_gap_max << ','
         << fragmentation_gap_mean << ','
         << csv_escape_progress(moment_audit) << ',' << moment_nodes << ','
         << moment_rel_max << ',' << moment_rel_mean << ','
         << moment_offset_max << ',' << moment_offset_mean << ','
         << moment_cancel_max << ',' << moment_cancel_mean << ','
         << moment_dynamic_range << ','
         << expansion_pairs << ',' << expansion_panels << ','
         << expansion_direct_ms / 1000.0 << ',' << expansion_rel_max << ','
         << expansion_rel_mean << ',' << expansion_rel_p95 << ','
         << expansion_rel_p99 << ','
         << csv_escape_progress(expansion_worst_channel) << ','
         << expansion_channel_rel_max[0] << ','
         << expansion_channel_rel_max[1] << ','
         << expansion_channel_rel_max[2] << ','
         << expansion_channel_rel_max[3] << ','
         << expansion_channel_rel_mean[0] << ','
         << expansion_channel_rel_mean[1] << ','
         << expansion_channel_rel_mean[2] << ','
         << expansion_channel_rel_mean[3] << ','
         << expansion_channel_rel_p99[0] << ','
         << expansion_channel_rel_p99[1] << ','
         << expansion_channel_rel_p99[2] << ','
         << expansion_channel_rel_p99[3] << ','
         << (expansion_mean_gate_pass ? 1 : 0)
         << ',' << convergence_rho_max << ','
         << convergence_rho_mean << ',' << node_aspect_max << ','
         << node_aspect_mean << ',' << accepted_aspect_max << ','
         << accepted_aspect_mean << ',' << leaf_panel << ','
         << leaf_tile_panel << ',' << leaf_residual << ',' << far_nodes << ','
         << leaf_mac_candidates << ',' << leaf_mac_accepted << ','
         << leaf_mac_rejected << ',' << leaf_mac_far_panels << ','
         << leaf_mac_direct_panels_saved_estimate << ','
         << csv_escape_progress(matvec_cost_profile) << ','
         << leaf_panels_per_second << ','
         << leaf_pair_batch_clean_source_panels_per_second << ','
         << tree_traversal_fraction << ','
         << csv_escape_progress(leaf_work_profile) << ','
         << leaf_work_mac_terminal_panels << ','
         << leaf_work_physical_gap_panels << ','
         << leaf_work_fragmentation_panels << ','
         << leaf_work_convergence_panels << ','
         << leaf_work_convergence_estimator_panels << ','
         << leaf_work_convergence_rho_panels << ','
         << leaf_work_convergence_aspect_panels << ','
         << leaf_work_surface_margin_panels << ','
         << leaf_work_surface_margin_hard_panels << ','
         << leaf_work_surface_margin_soft_panels << ','
         << leaf_work_k3_rho_cap_panels << ','
         << leaf_work_singular_fallback_panels << ','
         << leaf_work_scalar_tail_panels << ','
         << leaf_work_profile_accounted_panels << ','
         << leaf_work_mac_terminal_frac << ','
         << leaf_work_convergence_frac << ','
         << leaf_work_surface_margin_frac << ','
         << leaf_work_singular_fallback_frac << ','
         << leaf_work_scalar_tail_frac << ','
         << csv_escape_progress(cluster_policy) << ','
         << target_batch_size << ',' << leaf_pair_min_source_panels << ','
         << leaf_pair_min_targets << ',' << leaf_pair_batch_interactions << ','
         << leaf_pair_batch_targets << ',' << leaf_pair_batch_source_panels
         << ',' << leaf_pair_batch_residual_panels << ','
         << leaf_pair_batch_fallbacks << ',' << leaf_pair_batch_ms / 1000.0
         << ',' << csv_escape_progress(leaf_pair_batch_kernel) << ','
         << leaf_pair_batch_vectorized_sources << ','
         << leaf_pair_batch_vector_lanes << ','
         << leaf_pair_batch_pairs_per_second
         << ',' << csv_escape_progress(leaf_kernel_profile)
         << ',' << csv_escape_progress(jcp_direct_kernel)
         << ',' << leaf_kernel_scalar_tasks
         << ',' << leaf_kernel_batch_tasks
         << ',' << leaf_kernel_source_panel_evals
         << ',' << leaf_kernel_scalar_ms / 1000.0
         << ',' << leaf_kernel_target_avx2_ms / 1000.0
         << ',' << leaf_kernel_tail_ms / 1000.0
         << ',' << leaf_pair_batch_soa_sources
         << ',' << leaf_pair_batch_soa_ms / 1000.0
         << ',' << leaf_pair_batch_scalar_tail_ms / 1000.0
         << ',' << leaf_pair_batch_clean_source_panels_per_second
         << ',' << csv_escape_progress(leaf_pair_selection_profile) << ','
         << csv_escape_progress(leaf_pair_probe_cache) << ','
         << csv_escape_progress(leaf_pair_target_group) << ','
         << leaf_pair_max_lookahead << ','
         << csv_escape_progress(leaf_pair_directness_probe) << ','
         << csv_escape_progress(leaf_pair_directness_cache) << ','
         << csv_escape_progress(leaf_pair_directness_verify) << ','
         << leaf_pair_selection_candidate_scans << ','
         << leaf_pair_selection_clean_rejects << ','
         << leaf_pair_selection_probe_calls << ','
         << leaf_pair_selection_probe_node_visits << ','
         << leaf_pair_selection_probe_ms / 1000.0 << ','
         << leaf_pair_selection_ms / 1000.0 << ','
         << leaf_pair_selection_probe_cache_hits << ','
         << leaf_pair_selection_probe_cache_misses << ','
         << leaf_pair_selection_insufficient_targets << ','
         << leaf_pair_selection_same_target_leaf_rejects
         << ',' << leaf_pair_directness_ancestor_probe_calls << ','
         << leaf_pair_directness_ancestor_probe_node_visits << ','
         << leaf_pair_directness_ancestor_probe_ms / 1000.0 << ','
         << leaf_pair_directness_pair_cache_hits << ','
         << leaf_pair_directness_pair_cache_misses << ','
         << leaf_pair_directness_all_direct_proofs << ','
         << leaf_pair_directness_mixed_fallbacks << ','
         << leaf_pair_directness_verify_checks << ','
         << leaf_pair_directness_verify_mismatches
         << ',' << csv_escape_progress(leaf_pair_scheduler) << ','
         << leaf_pair_worklist_max_active << ','
         << csv_escape_progress(leaf_pair_worklist_verify) << ','
         << leaf_pair_worklist_tasks << ',' << leaf_pair_worklist_groups
         << ',' << leaf_pair_worklist_batched_groups << ','
         << leaf_pair_worklist_scalar_flush_tasks << ','
         << leaf_pair_worklist_forced_flushes << ','
         << leaf_pair_worklist_active_groups_max << ','
         << leaf_pair_worklist_ms / 1000.0 << ','
         << leaf_pair_worklist_verify_checks << ','
         << leaf_pair_worklist_verify_mismatches
         << ',' << csv_escape_progress(leaf_block_verify) << ','
         << leaf_block_target_leaves << ',' << leaf_block_targets << ','
         << leaf_block_source_groups << ',' << leaf_block_batched_groups
         << ',' << leaf_block_scalar_tail_tasks << ','
         << leaf_block_group_ms / 1000.0 << ','
         << leaf_block_flush_ms / 1000.0 << ','
         << leaf_block_writeback_ms / 1000.0 << ','
         << leaf_block_verify_checks << ','
         << leaf_block_verify_mismatches
         << ','
         << singular << ','
         << csv_escape_progress(result.preconditioner_effective) << ','
         << result.preconditioner_block_panel_cap << ','
         << result.preconditioner_overlap_panel_cap << ','
         << result.preconditioner_ilu_neighbor_cap << ','
         << result.preconditioner_ilu_drop_tol << ','
         << result.preconditioner_memory_budget_bytes << ','
         << result.preconditioner_estimated_bytes << ','
         << result.preconditioner_valid_blocks << ','
         << result.preconditioner_valid_panels << ','
         << result.preconditioner_ill_conditioned_blocks << ','
         << result.preconditioner_quality_rejected_blocks << ','
         << result.preconditioner_quality_probe_max_gain << ','
         << csv_escape_progress(result.preconditioner_quality_status) << ','
         << result.preconditioner_setup_ms / 1000.0 << ','
         << result.preconditioner_apply_ms / 1000.0 << ','
         << result.preconditioner_apply_calls << ','
         << result.preconditioner_sparse_nnz_blocks << ','
         << result.preconditioner_sparse_avg_degree << ','
         << result.preconditioner_sparse_max_degree << ','
         << result.preconditioner_coarse_aggregate_cap << ','
         << result.preconditioner_coarse_aggregates << ','
         << result.preconditioner_coarse_dimension << ','
         << csv_escape_progress(result.preconditioner_coarse_base) << ','
         << result.preconditioner_coarse_matrix_bytes << ','
         << result.preconditioner_coarse_pivot_ratio << ','
         << csv_escape_progress(result.preconditioner_coarse_status) << ','
         << result.preconditioner_coarse_setup_matvec_calls << ','
         << result.preconditioner_coarse_setup_ms / 1000.0 << ','
         << result.preconditioner_coarse_apply_ms / 1000.0 << ','
         << result.preconditioner_coarse_apply_calls << ','
         << result.preconditioner_coarse_rp_identity_error_max << ','
         << result.preconditioner_coarse_verify_rel_error_max << ','
         << threads_ << ',' << current_working_set_bytes() << ','
         << csv_escape_progress(status) << '\n';
    out_.flush();
  }

private:
  bool enabled_{false};
  std::chrono::seconds interval_{5};
  std::chrono::steady_clock::time_point start_{};
  std::chrono::steady_clock::time_point last_write_{};
  std::size_t threads_{1};
  std::ofstream out_;
};

std::string residual_history_string(const std::vector<double> &history) {
  std::ostringstream oss;
  oss.setf(std::ios::scientific);
  oss.precision(6);
  for (std::size_t i = 0; i < history.size(); ++i) {
    if (i != 0)
      oss << ';';
    oss << history[i];
  }
  return oss.str();
}

FormulationConfig make_formulation_config(const Options &options) {
  FormulationConfig config;
  config.formulation = options.formulation;
  config.energy_projection = options.energy_projection;
  if (config.formulation.empty())
    config.formulation = "jcp2013-centroid";
  if (config.formulation == "jcp2013")
    config.formulation = "jcp2013-centroid";
  if (config.energy_projection.empty())
    config.energy_projection = "boundary";
  if (config.formulation == "jcp2013-centroid" ||
      config.formulation == "jcp2013-treecode" ||
      config.formulation == "tabi2-nodepatch") {
    // The jcp2013 path is intentionally explicit even while it shares the
    // current scale factors. Accuracy gates use the audit CSV to keep future
    // sign/scale changes visible instead of silently drifting the formulation.
    config.rhs_scale = 1.0;
    config.operator_contrast_scale = 1.0;
    config.boundary_energy_scale = 1.0;
  } else if (config.formulation == "calibrated-debug") {
    config.rhs_scale = 1.0;
    config.operator_contrast_scale = 1.0;
    config.boundary_energy_scale = 1.0;
  }
  return config;
}

struct LinearSolveResult {
  std::vector<double> x;
  std::vector<double> best_x;
  int iterations{0};
  double residual{0.0};
  int restart_requested{32};
  int restart_effective{32};
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
  std::string convergence_reason{"not_started"};
  std::vector<double> residual_history;
  double orthogonalization_ms{0.0};
  double preconditioner_ms{0.0};
  double preconditioner_apply_ms{0.0};
  std::size_t preconditioner_apply_calls{0};
};

using GmresProgressCallback =
    std::function<void(int gmres_iter, int restart_iter,
                       double true_residual_rel,
                       double preconditioned_residual_rel, bool force,
                       const std::string &phase,
                       const LinearSolveResult &state)>;

struct SampledResidualResult {
  double abs{0.0};
  double rel{std::numeric_limits<double>::infinity()};
  std::size_t rows{0};
  bool finite{false};
};

using SampledTrueResidualCallback =
    std::function<SampledResidualResult(const std::vector<double> &,
                                        const std::vector<std::size_t> &)>;

LinearSolveResult gmres_solve(
    const std::function<void(const std::vector<double> &,
                             std::vector<double> &)> &matvec,
    const std::vector<double> &b, const std::vector<double> &diagonal,
    int max_iter, double tol, int restart_requested, TitanExecutor &executor,
    const std::function<void(const std::vector<double> &,
                             std::vector<double> &)> &custom_preconditioner =
        {},
    const GmresProgressCallback &progress_callback = {},
    const std::string &true_residual_probe_mode = "restart",
    int true_residual_probe_interval = 8,
    const std::string &true_residual_probe_scope = "full",
    int true_residual_sample_count = 1024,
    const SampledTrueResidualCallback &sampled_true_residual = {}) {
  LinearSolveResult result;
  const std::size_t n = b.size();
  result.x.assign(n, 0.0);
  result.restart_requested = std::clamp(restart_requested, 8, 128);
  result.restart_effective = result.restart_requested;
  result.true_residual_probe_mode =
      true_residual_probe_mode == "interval"
          ? std::string("interval")
          : (true_residual_probe_mode == "off" ? std::string("off")
                                                : std::string("restart"));
  result.true_residual_probe_scope =
      true_residual_probe_scope == "sampled"
          ? std::string("sampled")
          : (true_residual_probe_scope == "adaptive" ? std::string("adaptive")
                                                     : std::string("full"));
  result.true_residual_probe_interval =
      std::clamp(true_residual_probe_interval, 1, 128);
  result.true_residual_sample_count =
      std::clamp(true_residual_sample_count, 128, 4096);
  if (n == 0) {
    result.converged = true;
    result.convergence_reason = "converged";
    result.best_residual = 0.0;
    result.best_residual_iter = 0;
    result.final_residual = 0.0;
    result.best_preconditioned_residual = 0.0;
    result.final_preconditioned_residual = 0.0;
    result.true_residual_abs = 0.0;
    result.true_residual_rel = 0.0;
    result.true_residual_probe_calls = 0;
    result.true_residual_probe_ms = 0.0;
    result.sampled_true_residual_rel = 0.0;
    result.sampled_true_residual_rows = 0;
    result.sampled_probe_calls = 0;
    result.sampled_probe_ms = 0.0;
    result.full_true_residual_confirm_calls = 0;
    result.full_true_residual_confirm_ms = 0.0;
    return result;
  }

  std::vector<double> inv_precond(n, 1.0);
  const auto precond_begin = std::chrono::steady_clock::now();
  std::atomic_bool preconditioner_ok{true};
  parallel_for_each_index(n, executor, [&](std::size_t i) {
    const double diag = i < diagonal.size() ? std::abs(diagonal[i]) : 1.0;
    if (!std::isfinite(diag) || diag <= 1.0e-12) {
      inv_precond[i] = 1.0;
      preconditioner_ok.store(false, std::memory_order_relaxed);
    } else {
      inv_precond[i] = 1.0 / diag;
    }
  });
  result.preconditioner_ms +=
      elapsed_ms(precond_begin, std::chrono::steady_clock::now());
  if (!preconditioner_ok.load(std::memory_order_relaxed))
    result.convergence_reason = "preconditioner_breakdown";

  auto apply_preconditioner = [&](const std::vector<double> &in,
                                  std::vector<double> &out) {
    const auto begin = std::chrono::steady_clock::now();
    if (custom_preconditioner) {
      custom_preconditioner(in, out);
      if (out.size() != n) {
        out.resize(n);
        for (std::size_t i = 0; i < n; ++i)
          out[i] = in[i] * inv_precond[i];
      }
    } else {
      out.resize(n);
      parallel_for_each_index(n, executor, [&](std::size_t i) {
        out[i] = in[i] * inv_precond[i];
      });
    }
    const double apply_ms = elapsed_ms(begin, std::chrono::steady_clock::now());
    result.preconditioner_ms += apply_ms;
    result.preconditioner_apply_ms += apply_ms;
    ++result.preconditioner_apply_calls;
  };

  auto finite_vector = [](const std::vector<double> &values) noexcept {
    for (double value : values) {
      if (!std::isfinite(value))
        return false;
    }
    return true;
  };

  const int restart = result.restart_effective;
  std::vector<double> ax(n, 0.0);
  std::vector<double> residual = b;
  std::vector<double> precond_residual(n, 0.0);
  apply_preconditioner(residual, precond_residual);
  if (!finite_vector(precond_residual)) {
    result.convergence_reason = "preconditioner_breakdown";
    return result;
  }
  const double b_norm_raw = l2_norm_parallel(b, executor);
  const double b_norm = std::max(1.0e-300, b_norm_raw);
  const double precond_norm_raw = l2_norm_parallel(precond_residual, executor);
  const double precond_norm = std::max(1.0e-300, precond_norm_raw);
  auto relative_true_residual = [&](const std::vector<double> &r) {
    const double abs_value = l2_norm_parallel(r, executor);
    return std::pair<double, double>{abs_value, abs_value / b_norm};
  };
  auto preconditioned_relative = [&](const std::vector<double> &r) {
    return l2_norm_parallel(r, executor) / precond_norm;
  };
  auto build_sample_rows = [](std::size_t row_count, int requested) {
    std::vector<std::size_t> rows;
    if (row_count == 0 || requested <= 0)
      return rows;
    const std::size_t m =
        std::min<std::size_t>(row_count, static_cast<std::size_t>(requested));
    rows.reserve(m);
    std::size_t previous = std::numeric_limits<std::size_t>::max();
    for (std::size_t k = 0; k < m; ++k) {
      std::size_t row =
          ((2u * k + 1u) * row_count) / std::max<std::size_t>(1, 2u * m);
      row = std::min(row, row_count - 1u);
      if (row != previous) {
        rows.push_back(row);
        previous = row;
      }
    }
    if (rows.empty())
      rows.push_back(0);
    return rows;
  };
  const std::vector<std::size_t> sampled_rows =
      build_sample_rows(n, result.true_residual_sample_count);
  auto sampled_relative_from_full_ax = [&](const std::vector<double> &ax_values) {
    double residual_sq = 0.0;
    double b_sq = 0.0;
    for (std::size_t row : sampled_rows) {
      const double r = b[row] - ax_values[row];
      residual_sq += r * r;
      b_sq += b[row] * b[row];
    }
    const double abs_value = std::sqrt(std::max(0.0, residual_sq));
    const double denom = std::max(1.0e-300, std::sqrt(std::max(0.0, b_sq)));
    SampledResidualResult sample;
    sample.abs = abs_value;
    sample.rel = abs_value / denom;
    sample.rows = sampled_rows.size();
    sample.finite = std::isfinite(sample.rel);
    return sample;
  };

  auto [initial_abs, initial_rel] = relative_true_residual(residual);
  result.residual = initial_rel;
  result.best_residual = initial_rel;
  result.best_residual_iter = 0;
  result.final_residual = initial_rel;
  result.true_residual_abs = initial_abs;
  result.true_residual_rel = initial_rel;
  result.best_x = result.x;
  double precond_rel = preconditioned_relative(precond_residual);
  result.best_preconditioned_residual = precond_rel;
  result.final_preconditioned_residual = precond_rel;
  result.residual_history.push_back(initial_rel);
  if (result.residual <= tol) {
    result.converged = true;
    result.convergence_reason = "converged_true_residual";
    return result;
  }

  const int iter_cap = std::max(1, max_iter);
  double best_true_abs = initial_abs;
  int non_improving_restarts = 0;
  bool saw_phantom_preconditioned = false;
  while (result.iterations < iter_cap) {
    const double beta = l2_norm_parallel(precond_residual, executor);
    if (!std::isfinite(beta) || beta <= 1.0e-300) {
      result.convergence_reason = "preconditioner_breakdown";
      return result;
    }
    if (beta <= tol * precond_norm) {
      saw_phantom_preconditioned = true;
      auto [true_abs, true_rel] = relative_true_residual(residual);
      result.final_residual = true_rel;
      result.true_residual_abs = true_abs;
      result.true_residual_rel = true_rel;
      result.final_preconditioned_residual = beta / precond_norm;
      if (std::isfinite(result.final_preconditioned_residual) &&
          result.final_preconditioned_residual <
              result.best_preconditioned_residual) {
        result.best_preconditioned_residual =
            result.final_preconditioned_residual;
      }
      if (true_rel <= tol) {
        result.converged = true;
        result.residual = true_rel;
        result.convergence_reason = "converged_true_residual";
        break;
      }
    }
    if (result.converged)
      break;

    const int window = std::min(restart, iter_cap - result.iterations);
    std::vector<std::vector<double>> v(
        static_cast<std::size_t>(window) + 1, std::vector<double>(n, 0.0));
    parallel_for_each_index(n, executor, [&](std::size_t i) {
      v[0][i] = precond_residual[i] / beta;
    });
    std::vector<double> h((static_cast<std::size_t>(window) + 1) *
                              static_cast<std::size_t>(window),
                          0.0);
    std::vector<double> cs(static_cast<std::size_t>(window), 0.0);
    std::vector<double> sn(static_cast<std::size_t>(window), 0.0);
    std::vector<double> g(static_cast<std::size_t>(window) + 1, 0.0);
    g[0] = beta;

    auto solve_window_y = [&](int cols) {
      std::vector<double> y(static_cast<std::size_t>(cols), 0.0);
      for (int rev = 0; rev < cols; ++rev) {
        const int row = cols - 1 - rev;
        double value = g[static_cast<std::size_t>(row)];
        for (int col = row + 1; col < cols; ++col) {
          value -= h[static_cast<std::size_t>(row) *
                         static_cast<std::size_t>(window) +
                     static_cast<std::size_t>(col)] *
                   y[static_cast<std::size_t>(col)];
        }
        const double diag =
            h[static_cast<std::size_t>(row) * static_cast<std::size_t>(window) +
              static_cast<std::size_t>(row)];
        y[static_cast<std::size_t>(row)] =
            std::abs(diag) > 1.0e-30 ? value / diag : 0.0;
      }
      return y;
    };

    auto build_candidate_x = [&](const std::vector<double> &base,
                                 const std::vector<double> &y, int cols,
                                 std::vector<double> &candidate) {
      candidate = base;
      parallel_for_each_index(n, executor, [&](std::size_t k) {
        double update = 0.0;
        for (int basis = 0; basis < cols; ++basis)
          update += y[static_cast<std::size_t>(basis)] *
                    v[static_cast<std::size_t>(basis)][k];
        candidate[k] += update;
      });
    };

    int used_cols = 0;
    bool happy_breakdown = false;
    bool preconditioned_requested_check = false;
    bool probe_converged = false;
    for (int j = 0; j < window; ++j) {
      std::vector<double> w(n, 0.0);
      matvec(v[static_cast<std::size_t>(j)], w);
      if (!finite_vector(w)) {
        result.convergence_reason = "nan_guard";
        return result;
      }
      std::vector<double> precond_w(n, 0.0);
      apply_preconditioner(w, precond_w);
      if (!finite_vector(precond_w)) {
        result.convergence_reason = "preconditioner_breakdown";
        return result;
      }
      w.swap(precond_w);

      const auto ortho_begin = std::chrono::steady_clock::now();
      for (int i = 0; i <= j; ++i) {
        const double hij =
            dot_parallel(w, v[static_cast<std::size_t>(i)], executor);
        h[static_cast<std::size_t>(i) * static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = hij;
        parallel_for_each_index(n, executor, [&](std::size_t k) {
          w[k] -= hij * v[static_cast<std::size_t>(i)][k];
        });
      }
      const double h_next = l2_norm_parallel(w, executor);
      h[(static_cast<std::size_t>(j) + 1) *
            static_cast<std::size_t>(window) +
        static_cast<std::size_t>(j)] = h_next;
      if (h_next > 1.0e-20 && j + 1 < window + 1) {
        parallel_for_each_index(n, executor, [&](std::size_t k) {
          v[static_cast<std::size_t>(j) + 1][k] = w[k] / h_next;
        });
      } else {
        happy_breakdown = true;
      }
      result.orthogonalization_ms +=
          elapsed_ms(ortho_begin, std::chrono::steady_clock::now());

      for (int i = 0; i < j; ++i) {
        const double h_i =
            h[static_cast<std::size_t>(i) * static_cast<std::size_t>(window) +
              static_cast<std::size_t>(j)];
        const double h_next_i =
            h[(static_cast<std::size_t>(i) + 1) *
                  static_cast<std::size_t>(window) +
              static_cast<std::size_t>(j)];
        h[static_cast<std::size_t>(i) * static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = cs[static_cast<std::size_t>(i)] * h_i +
                                         sn[static_cast<std::size_t>(i)] *
                                             h_next_i;
        h[(static_cast<std::size_t>(i) + 1) *
              static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = -sn[static_cast<std::size_t>(i)] * h_i +
                                         cs[static_cast<std::size_t>(i)] *
                                             h_next_i;
      }

      const double h_j =
          h[static_cast<std::size_t>(j) * static_cast<std::size_t>(window) +
            static_cast<std::size_t>(j)];
      const double h_next_j =
          h[(static_cast<std::size_t>(j) + 1) *
                static_cast<std::size_t>(window) +
            static_cast<std::size_t>(j)];
      const double rho = std::hypot(h_j, h_next_j);
      if (rho <= 1.0e-30 || !std::isfinite(rho)) {
        happy_breakdown = true;
        h[static_cast<std::size_t>(j) * static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = h_j;
        h[(static_cast<std::size_t>(j) + 1) *
              static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = 0.0;
      } else {
        cs[static_cast<std::size_t>(j)] = h_j / rho;
        sn[static_cast<std::size_t>(j)] = h_next_j / rho;
        h[static_cast<std::size_t>(j) * static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = rho;
        h[(static_cast<std::size_t>(j) + 1) *
              static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = 0.0;
        const double g_j = g[static_cast<std::size_t>(j)];
        const double g_next = g[static_cast<std::size_t>(j) + 1];
        g[static_cast<std::size_t>(j)] =
            cs[static_cast<std::size_t>(j)] * g_j +
            sn[static_cast<std::size_t>(j)] * g_next;
        g[static_cast<std::size_t>(j) + 1] =
            -sn[static_cast<std::size_t>(j)] * g_j +
            cs[static_cast<std::size_t>(j)] * g_next;
      }

      ++result.iterations;
      used_cols = j + 1;
      precond_rel = std::abs(g[static_cast<std::size_t>(j) + 1]) /
                    precond_norm;
      result.final_preconditioned_residual = precond_rel;
      if (std::isfinite(precond_rel) &&
          precond_rel < result.best_preconditioned_residual) {
        result.best_preconditioned_residual = precond_rel;
      }
      if (progress_callback)
        progress_callback(result.iterations, j + 1, result.final_residual,
                          precond_rel, false,
                          "gmres_iter", result);
      if (result.true_residual_probe_mode == "interval" &&
          result.true_residual_probe_interval > 0 &&
          (result.iterations % result.true_residual_probe_interval) == 0) {
        std::vector<double> y_probe = solve_window_y(used_cols);
        std::vector<double> probe_x;
        build_candidate_x(result.x, y_probe, used_cols, probe_x);
        std::vector<double> probe_ax;
        bool probe_ax_ready = false;
        bool run_full_probe = result.true_residual_probe_scope == "full";
        if (result.true_residual_probe_scope == "sampled" ||
            result.true_residual_probe_scope == "adaptive") {
          const auto sample_begin = std::chrono::steady_clock::now();
          SampledResidualResult sample;
          if (sampled_true_residual) {
            sample = sampled_true_residual(probe_x, sampled_rows);
          } else {
            probe_ax.assign(n, 0.0);
            matvec(probe_x, probe_ax);
            probe_ax_ready = true;
            if (!finite_vector(probe_ax)) {
              result.convergence_reason = "nan_guard";
              return result;
            }
            sample = sampled_relative_from_full_ax(probe_ax);
          }
          ++result.sampled_probe_calls;
          result.sampled_probe_ms +=
              elapsed_ms(sample_begin, std::chrono::steady_clock::now());
          result.sampled_true_residual_rel = sample.rel;
          result.sampled_true_residual_rows = sample.rows;
          if (!sample.finite) {
            result.convergence_reason = "nan_guard";
            return result;
          }
          if (progress_callback)
            progress_callback(result.iterations, j + 1, sample.rel,
                              precond_rel, true,
                              "gmres_sampled_true_probe", result);
          if (result.true_residual_probe_scope == "adaptive" &&
              sample.rel <= std::max(5.0 * tol, 1.0e-3)) {
            run_full_probe = true;
          }
        }
        if (run_full_probe) {
          const auto probe_begin = std::chrono::steady_clock::now();
          if (!probe_ax_ready) {
            probe_ax.assign(n, 0.0);
            matvec(probe_x, probe_ax);
            probe_ax_ready = true;
            if (!finite_vector(probe_ax)) {
              result.convergence_reason = "nan_guard";
              return result;
            }
          }
          std::vector<double> probe_residual(n, 0.0);
          parallel_for_each_index(n, executor, [&](std::size_t k) {
            probe_residual[k] = b[k] - probe_ax[k];
          });
          auto [probe_abs, probe_rel] =
              relative_true_residual(probe_residual);
          ++result.true_residual_probe_calls;
          ++result.full_true_residual_confirm_calls;
          const double probe_ms =
              elapsed_ms(probe_begin, std::chrono::steady_clock::now());
          result.true_residual_probe_ms += probe_ms;
          result.full_true_residual_confirm_ms += probe_ms;
          result.final_residual = probe_rel;
          result.true_residual_abs = probe_abs;
          result.true_residual_rel = probe_rel;
          if (std::isfinite(probe_rel) && probe_rel < result.best_residual) {
            result.best_residual = probe_rel;
            result.best_residual_iter = result.iterations;
            best_true_abs = probe_abs;
            result.best_x = probe_x;
          }
          result.residual_history.push_back(probe_rel);
          if (progress_callback)
            progress_callback(result.iterations, j + 1, probe_rel, precond_rel,
                              true, "gmres_true_probe", result);
          if (probe_rel <= tol) {
            result.x = std::move(probe_x);
            result.residual = probe_rel;
            result.converged = true;
            result.convergence_reason = "converged_true_residual_probe";
            probe_converged = true;
            break;
          }
        }
      }
      if (precond_rel <= tol) {
        preconditioned_requested_check = true;
        saw_phantom_preconditioned = true;
      }
      if (result.iterations >= iter_cap)
        break;
      if (happy_breakdown)
        break;
    }

    if (probe_converged)
      break;
    if (used_cols == 0)
      break;

    std::vector<double> y = solve_window_y(used_cols);
    build_candidate_x(result.x, y, used_cols, result.x);

    matvec(result.x, ax);
    if (!finite_vector(ax)) {
      result.convergence_reason = "nan_guard";
      return result;
    }
    parallel_for_each_index(n, executor, [&](std::size_t k) {
      residual[k] = b[k] - ax[k];
    });
    auto [true_abs, true_rel] = relative_true_residual(residual);
    result.final_residual = true_rel;
    result.true_residual_abs = true_abs;
    result.true_residual_rel = true_rel;
    result.residual = true_rel;
    bool improved_true = false;
    if (std::isfinite(true_rel) && true_rel < result.best_residual) {
      result.best_residual = true_rel;
      result.best_residual_iter = result.iterations;
      best_true_abs = true_abs;
      result.best_x = result.x;
      improved_true = true;
    }
    result.residual_history.push_back(true_rel);
    apply_preconditioner(residual, precond_residual);
    if (!finite_vector(precond_residual)) {
      result.convergence_reason = "preconditioner_breakdown";
      return result;
    }
    precond_rel = preconditioned_relative(precond_residual);
    result.final_preconditioned_residual = precond_rel;
    if (std::isfinite(precond_rel) &&
        precond_rel < result.best_preconditioned_residual) {
      result.best_preconditioned_residual = precond_rel;
    }
    if (progress_callback)
      progress_callback(result.iterations, 0, true_rel, precond_rel, true,
                        "gmres_restart", result);
    if (true_rel <= tol) {
      result.converged = true;
      result.convergence_reason =
          preconditioned_requested_check
              ? "converged_true_residual_after_preconditioned_check"
              : "converged_true_residual";
      break;
    }
    if (preconditioned_requested_check && true_rel > tol) {
      saw_phantom_preconditioned = true;
      non_improving_restarts = 0;
    } else if (!improved_true) {
      ++non_improving_restarts;
    } else {
      non_improving_restarts = 0;
    }
    if (non_improving_restarts >= 2) {
      result.convergence_reason =
          result.iterations >= iter_cap ? "max_iter" : "stagnated";
      break;
    }
  }
  if (!result.best_x.empty()) {
    result.x = result.best_x;
    result.residual = result.best_residual;
    result.true_residual_abs = best_true_abs;
    result.true_residual_rel = result.best_residual;
  }
  if (!result.converged && result.best_residual <= tol) {
    result.converged = true;
    result.convergence_reason = result.final_residual <= tol
                                    ? "converged_true_residual"
                                    : "converged_true_residual_best";
  }
  if (!result.converged && result.convergence_reason == "not_started")
    result.convergence_reason = saw_phantom_preconditioned
                                    ? "phantom_preconditioned"
                                    : (result.iterations >= iter_cap
                                           ? "max_iter"
                                           : "stagnated");
  return result;
}

// TABI v1 convention: linear boundary-integral solve for the induced
// surface density. Kernels use Angstrom distances, elementary-charge inputs,
// inverse-Angstrom kappa, outward normals, and Coulomb conversion only in the
// final reaction-field energy projection.
std::vector<double> build_rhs(const std::vector<Atom> &atoms,
                              const std::vector<Panel> &panels,
                              double contrast, double kappa,
                              double rhs_scale) {
  std::vector<double> rhs(panels.size(), 0.0);
  for (std::size_t p = 0; p < panels.size(); ++p) {
    const Panel &panel = panels[p];
    double normal_field = 0.0;
    for (const Atom &atom : atoms) {
      const double dx = panel.x - atom.x;
      const double dy = panel.y - atom.y;
      const double dz = panel.z - atom.z;
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 <= 1.0e-16)
        continue;
      const double r = std::sqrt(r2);
      const double normal_projection =
          (panel.nx * dx + panel.ny * dy + panel.nz * dz) / r;
      const double screened = std::exp(-kappa * r);
      normal_field += atom.q * screened * normal_projection /
                      std::max(1.0e-12, r2);
    }
    rhs[p] = -rhs_scale * contrast * normal_field / (4.0 * kPi);
  }
  return rhs;
}

double boundary_energy_kcal(const std::vector<Atom> &atoms,
                            const std::vector<Panel> &panels,
                            const std::vector<double> &density,
                            double kappa) {
  double energy = 0.0;
  for (const Atom &atom : atoms) {
    double phi = 0.0;
    for (std::size_t p = 0; p < panels.size(); ++p) {
      const Panel &panel = panels[p];
      const double dx = atom.x - panel.x;
      const double dy = atom.y - panel.y;
      const double dz = atom.z - panel.z;
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 <= 1.0e-16)
        continue;
      const double r = std::sqrt(r2);
      phi += density[p] * panel.area * std::exp(-kappa * r) /
             (4.0 * kPi * r);
    }
    energy += 0.5 * kCoulombKcalPerMol * atom.q * phi;
  }
  return energy;
}

double born_reference_kcal(const std::vector<Atom> &atoms, double eps_in,
                           double eps_out, double kappa) {
  double energy = 0.0;
  for (const Atom &atom : atoms) {
    const double radius = std::max(0.5, atom.radius);
    const double outside =
        1.0 / (std::max(1.0e-12, eps_out) * (1.0 + kappa * radius));
    const double inside = 1.0 / std::max(1.0e-12, eps_in);
    energy +=
        0.5 * kCoulombKcalPerMol * atom.q * atom.q * (outside - inside) /
        radius;
  }
  return energy;
}

double green_kernel(double r, double kappa) noexcept {
  if (r <= 1.0e-14)
    return 0.0;
  return std::exp(-kappa * r) / (4.0 * kPi * r);
}

double source_normal_derivative_green(double tx, double ty, double tz,
                                      const Panel &source,
                                      double kappa) noexcept {
  const double ux = tx - source.x;
  const double uy = ty - source.y;
  const double uz = tz - source.z;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  const double ndotu = source.nx * ux + source.ny * uy + source.nz * uz;
  return std::exp(-kappa * r) * (1.0 + kappa * r) * ndotu /
         (4.0 * kPi * r2 * r);
}

double target_normal_derivative_green(const Panel &target, double sx,
                                      double sy, double sz,
                                      double kappa) noexcept {
  const double ux = target.x - sx;
  const double uy = target.y - sy;
  const double uz = target.z - sz;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  const double ndotu = target.nx * ux + target.ny * uy + target.nz * uz;
  return -std::exp(-kappa * r) * (1.0 + kappa * r) * ndotu /
         (4.0 * kPi * r2 * r);
}

std::vector<double> build_jcp_rhs(const std::vector<Atom> &atoms,
                                  const std::vector<Panel> &panels,
                                  double eps_in, bool coupled,
                                  double *s1_norm, double *s2_norm) {
  const std::size_t n = panels.size();
  std::vector<double> rhs(coupled ? 2 * n : n, 0.0);
  const double inv_eps_in = 1.0 / std::max(1.0e-12, eps_in);
  for (std::size_t i = 0; i < n; ++i) {
    const Panel &panel = panels[i];
    double s1 = 0.0;
    double s2 = 0.0;
    for (const Atom &atom : atoms) {
      const double dx = panel.x - atom.x;
      const double dy = panel.y - atom.y;
      const double dz = panel.z - atom.z;
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (r <= 1.0e-14)
        continue;
      s1 += atom.q * inv_eps_in * green_kernel(r, 0.0);
      s2 += atom.q * inv_eps_in *
            target_normal_derivative_green(panel, atom.x, atom.y, atom.z, 0.0);
    }
    rhs[i] = s1;
    if (coupled)
      rhs[n + i] = s2;
  }
  if (s1_norm)
    *s1_norm = rms_norm(std::vector<double>(rhs.begin(), rhs.begin() + n));
  if (s2_norm) {
    if (coupled)
      *s2_norm = rms_norm(std::vector<double>(rhs.begin() + n, rhs.end()));
    else
      *s2_norm = 0.0;
  }
  return rhs;
}

void apply_jcp_linear_pb_operator(const std::vector<Panel> &panels,
                                  const std::vector<double> &x,
                                  double eps_ratio, double kappa,
                                  double lhs_phi, double lhs_phin,
                                  const std::string &near_integration,
                                  const std::string &quad_buffer,
                                  TitanExecutor &executor,
                                  std::vector<double> &out,
                                  TreecodeStats *stats) {
  const std::size_t n = panels.size();
  out.assign(2 * n, 0.0);
  std::vector<std::size_t> order(n, 0);
  std::iota(order.begin(), order.end(), std::size_t{0});
  const std::size_t workers = std::max<std::size_t>(
      1, std::min<std::size_t>(executor.capacity(), n));
  std::vector<TreecodeStats> local_stats(workers);
  const auto near_begin = std::chrono::steady_clock::now();
  executor.parallel_for(
      n, workers, 4,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
          TreecodeStats &worker_stats = local_stats[worker];
          const JcpNearAccum sum = accumulate_jcp_near_ordered(
              panels[i], i, panels, x, eps_ratio, kappa, true, order, 0, n,
              near_integration, quad_buffer, &worker_stats);
          double phi_row = lhs_phi * x[i];
          double phin_row = lhs_phin * x[n + i];
          phi_row -= sum.phi;
          phin_row -= sum.phin;
          out[i] = phi_row;
          out[n + i] = phin_row;
        }
      });
  if (stats != nullptr) {
    stats->near_integration = near_integration;
    stats->quad_buffer = quad_buffer;
    const double near_ms =
        elapsed_ms(near_begin, std::chrono::steady_clock::now());
    stats->near_integration_ms += near_ms;
    stats->near_eval_ms += near_ms;
    std::size_t local_centroid = 0;
    std::size_t local_duffy = 0;
    std::size_t local_geometric = 0;
    for (const TreecodeStats &worker_stats : local_stats) {
      const std::size_t previous_duffy = stats->near_duffy_pairs;
      stats->near_kernel_policy = worker_stats.near_kernel_policy;
      stats->near_self_pairs += worker_stats.near_self_pairs;
      stats->near_one_ring_pairs += worker_stats.near_one_ring_pairs;
      stats->near_duffy_pairs += worker_stats.near_duffy_pairs;
      stats->near_centroid_pairs += worker_stats.near_centroid_pairs;
      stats->near_fallbacks += worker_stats.near_fallbacks;
      merge_near_kernel_block_stats(*stats, worker_stats);
      local_centroid += worker_stats.near_centroid_pairs;
      local_duffy += worker_stats.near_duffy_pairs;
      local_geometric += worker_stats.near_geometric_buffer_pairs;
      stats->duffy_subdivisions_max =
          std::max(stats->duffy_subdivisions_max,
                   worker_stats.duffy_subdivisions_max);
      if (worker_stats.near_duffy_pairs > 0) {
        const double old_count =
            static_cast<double>(previous_duffy);
        const double new_count =
            static_cast<double>(stats->near_duffy_pairs);
        stats->duffy_subdivisions_mean =
            new_count <= 0.0
                ? 0.0
                : (stats->duffy_subdivisions_mean * old_count +
                   worker_stats.duffy_subdivisions_mean *
                       static_cast<double>(worker_stats.near_duffy_pairs)) /
                      new_count;
      }
    }
    stats->near_interactions += local_centroid + local_duffy + local_geometric;
  }
}

void apply_jcp_poisson_operator(const std::vector<Panel> &panels,
                                const std::vector<double> &x,
                                double eps_ratio, double lhs_phi,
                                const std::string &near_integration,
                                const std::string &quad_buffer,
                                TitanExecutor &executor,
                                std::vector<double> &out,
                                TreecodeStats *stats) {
  const std::size_t n = panels.size();
  out.assign(n, 0.0);
  std::vector<std::size_t> order(n, 0);
  std::iota(order.begin(), order.end(), std::size_t{0});
  const std::size_t workers = std::max<std::size_t>(
      1, std::min<std::size_t>(executor.capacity(), n));
  std::vector<TreecodeStats> local_stats(workers);
  const auto near_begin = std::chrono::steady_clock::now();
  executor.parallel_for(
      n, workers, 4,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
          TreecodeStats &worker_stats = local_stats[worker];
          const JcpNearAccum sum = accumulate_jcp_near_ordered(
              panels[i], i, panels, x, eps_ratio, 0.0, false, order, 0, n,
              near_integration, quad_buffer, &worker_stats);
          double row = lhs_phi * x[i];
          row -= sum.phi;
          out[i] = row;
        }
      });
  if (stats != nullptr) {
    stats->near_integration = near_integration;
    stats->quad_buffer = quad_buffer;
    const double near_ms =
        elapsed_ms(near_begin, std::chrono::steady_clock::now());
    stats->near_integration_ms += near_ms;
    stats->near_eval_ms += near_ms;
    std::size_t local_centroid = 0;
    std::size_t local_duffy = 0;
    std::size_t local_geometric = 0;
    for (const TreecodeStats &worker_stats : local_stats) {
      const std::size_t previous_duffy = stats->near_duffy_pairs;
      stats->near_kernel_policy = worker_stats.near_kernel_policy;
      stats->near_self_pairs += worker_stats.near_self_pairs;
      stats->near_one_ring_pairs += worker_stats.near_one_ring_pairs;
      stats->near_duffy_pairs += worker_stats.near_duffy_pairs;
      stats->near_centroid_pairs += worker_stats.near_centroid_pairs;
      stats->near_fallbacks += worker_stats.near_fallbacks;
      merge_near_kernel_block_stats(*stats, worker_stats);
      local_centroid += worker_stats.near_centroid_pairs;
      local_duffy += worker_stats.near_duffy_pairs;
      local_geometric += worker_stats.near_geometric_buffer_pairs;
      stats->duffy_subdivisions_max =
          std::max(stats->duffy_subdivisions_max,
                   worker_stats.duffy_subdivisions_max);
      if (worker_stats.near_duffy_pairs > 0) {
        const double old_count =
            static_cast<double>(previous_duffy);
        const double new_count =
            static_cast<double>(stats->near_duffy_pairs);
        stats->duffy_subdivisions_mean =
            new_count <= 0.0
                ? 0.0
                : (stats->duffy_subdivisions_mean * old_count +
                   worker_stats.duffy_subdivisions_mean *
                       static_cast<double>(worker_stats.near_duffy_pairs)) /
                      new_count;
      }
    }
    stats->near_interactions += local_centroid + local_duffy + local_geometric;
  }
}

double jcp_linear_pb_energy_kcal(const std::vector<Atom> &atoms,
                                 const std::vector<Panel> &panels,
                                 const std::vector<double> &x,
                                 double eps_ratio, double kappa) {
  const std::size_t n = panels.size();
  double energy = 0.0;
  for (const Atom &atom : atoms) {
    double reaction_phi = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      const Panel &source = panels[j];
      const double dx = atom.x - source.x;
      const double dy = atom.y - source.y;
      const double dz = atom.z - source.z;
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (r <= 1.0e-14)
        continue;
      const double k1 = green_kernel(r, 0.0) - green_kernel(r, kappa);
      const double k2 =
          eps_ratio *
              source_normal_derivative_green(atom.x, atom.y, atom.z, source,
                                             kappa) -
          source_normal_derivative_green(atom.x, atom.y, atom.z, source, 0.0);
      reaction_phi += source.area * (k1 * x[n + j] + k2 * x[j]);
    }
    energy += 0.5 * kCoulombKcalPerMol * atom.q * reaction_phi;
  }
  return energy;
}

double jcp_poisson_energy_kcal(const std::vector<Atom> &atoms,
                               const std::vector<Panel> &panels,
                               const std::vector<double> &phi,
                               double eps_ratio) {
  const double k2_scale = eps_ratio - 1.0;
  double energy = 0.0;
  for (const Atom &atom : atoms) {
    double reaction_phi = 0.0;
    for (std::size_t j = 0; j < panels.size(); ++j) {
      const Panel &source = panels[j];
      const double k2 =
          k2_scale * source_normal_derivative_green(atom.x, atom.y, atom.z,
                                                    source, 0.0);
      reaction_phi += source.area * k2 * phi[j];
    }
    energy += 0.5 * kCoulombKcalPerMol * atom.q * reaction_phi;
  }
  return energy;
}

} // namespace

Result solve_linear_pb(const std::vector<Atom> &atoms, const Options &options,
                       TitanExecutor &executor) {
  Result result;
  result.tree_preset =
      (options.tree_preset == "stage3u-safe" ||
       options.tree_preset == "stage3u-perf" ||
       options.tree_preset == "stage3u-perf2")
          ? options.tree_preset
          : std::string("off");
  const FormulationConfig formulation = make_formulation_config(options);
  result.formulation = formulation.formulation;
  result.energy_projection = formulation.energy_projection;
  result.preconditioner_requested =
      options.preconditioner.empty() ? std::string("auto")
                                     : options.preconditioner;
  result.preconditioner_block_near =
      options.preconditioner_block_near == "buffer-aware" ||
              options.preconditioner_block_near == "operator-near"
          ? options.preconditioner_block_near
          : std::string("singular");
  const int preconditioner_block_cap_max =
      result.preconditioner_requested == "fmm-block-jacobi" ? 128 : 32;
  result.preconditioner_block_panel_cap = std::clamp(
      options.preconditioner_block_panel_cap, 8, preconditioner_block_cap_max);
  result.preconditioner_overlap_panel_cap =
      std::clamp(options.preconditioner_overlap_panel_cap, 0, 128);
  result.preconditioner_ilu_neighbor_cap =
      std::clamp(options.preconditioner_ilu_neighbor_cap, 32, 1024);
  result.preconditioner_ilu_drop_tol =
      std::clamp(options.preconditioner_ilu_drop_tol, 0.0, 1.0e-2);
  result.preconditioner_memory_budget_bytes =
      static_cast<std::size_t>(
          std::clamp(options.preconditioner_memory_budget_mb, 16, 512)) *
      1024u * 1024u;
  result.preconditioner_coarse_aggregate_cap =
      std::clamp(options.preconditioner_coarse_aggregate_cap, 4, 64);
  result.preconditioner_coarse_base =
      options.preconditioner_coarse_base == "block-budgeted"
          ? std::string("block-budgeted")
          : std::string("diagonal");
  result.preconditioner_effective = "diagonal";
  result.gmres_restart_requested = std::clamp(options.gmres_restart, 8, 128);
  result.gmres_restart_effective = result.gmres_restart_requested;
  result.best_residual = 0.0;
  result.best_residual_iter = 0;
  result.final_residual = 0.0;
  result.best_preconditioned_residual = 0.0;
  result.final_preconditioned_residual = 0.0;
  result.true_residual_abs = 0.0;
  result.true_residual_rel = 0.0;
  result.true_residual_probe_mode =
      options.true_residual_probe == "interval"
          ? std::string("interval")
          : (options.true_residual_probe == "off" ? std::string("off")
                                                   : std::string("restart"));
  result.true_residual_probe_scope =
      options.true_residual_probe_scope == "sampled"
          ? std::string("sampled")
          : (options.true_residual_probe_scope == "adaptive"
                 ? std::string("adaptive")
                 : std::string("full"));
  result.true_residual_probe_interval =
      std::clamp(options.true_residual_probe_interval, 1, 128);
  result.true_residual_sample_count =
      std::clamp(options.true_residual_sample_count, 128, 4096);
  result.true_residual_probe_calls = 0;
  result.true_residual_probe_ms = 0.0;
  result.sampled_true_residual_rel = 0.0;
  result.sampled_true_residual_rows = 0;
  result.sampled_probe_calls = 0;
  result.sampled_probe_ms = 0.0;
  result.full_true_residual_confirm_calls = 0;
  result.full_true_residual_confirm_ms = 0.0;
  result.adaptive_level =
      options.adaptive_level.empty() ? std::string("strict")
                                     : options.adaptive_level;
  result.cluster_policy =
      (options.cluster_policy == "leaf-pair-strict" ||
       options.cluster_policy == "leaf-pair-balanced")
          ? options.cluster_policy
          : std::string("off");
  result.target_batch_size =
      std::clamp(options.target_batch_size, 1, 8);
  if (result.target_batch_size != 4 && result.target_batch_size != 8)
    result.target_batch_size = 1;
  result.leaf_pair_min_source_panels =
      std::clamp(options.leaf_pair_min_source_panels, 16, 512);
  result.leaf_pair_min_targets =
      std::clamp(options.leaf_pair_min_targets, 2, 16);
  result.leaf_kernel_profile =
      options.leaf_kernel_profile == "summary" ? "summary" : "off";
  result.leaf_work_profile =
      options.leaf_work_profile == "summary" ? "summary" : "off";
  result.matvec_cost_profile =
      options.matvec_cost_profile == "summary" ? "summary" : "off";
  result.jcp_direct_kernel =
      normalize_jcp_direct_kernel_name(options.jcp_direct_kernel);
  result.leaf_pair_selection_profile =
      options.leaf_pair_selection_profile == "summary" ? "summary" : "off";
  result.leaf_pair_probe_cache =
      options.leaf_pair_probe_cache == "per-matvec" ? "per-matvec" : "off";
  result.leaf_pair_target_group =
      options.leaf_pair_target_group == "same-target-leaf"
          ? "same-target-leaf"
          : "scan";
  result.leaf_pair_directness_probe =
      options.leaf_pair_directness_probe == "ancestor-path"
          ? "ancestor-path"
          : "full-traversal";
  result.leaf_pair_directness_cache =
      options.leaf_pair_directness_cache == "target-leaf-source-leaf"
          ? "target-leaf-source-leaf"
          : "off";
  result.leaf_pair_directness_verify =
      (options.leaf_pair_directness_verify == "sampled" ||
       options.leaf_pair_directness_verify == "full-small")
          ? options.leaf_pair_directness_verify
          : "off";
  result.leaf_pair_scheduler =
      (options.leaf_pair_scheduler == "deferred-worklist" ||
       options.leaf_pair_scheduler == "target-leaf-block" ||
       options.leaf_pair_scheduler == "cost-tile-worklist")
          ? options.leaf_pair_scheduler
          : "inline-select";
  result.leaf_pair_worklist_max_active =
      std::clamp(options.leaf_pair_worklist_max_active, 64, 8192);
  result.leaf_pair_worklist_verify =
      (options.leaf_pair_worklist_verify == "sampled" ||
       options.leaf_pair_worklist_verify == "full-small")
          ? options.leaf_pair_worklist_verify
          : "off";
  result.leaf_block_verify =
      (options.leaf_block_verify == "sampled" ||
       options.leaf_block_verify == "full-small")
          ? options.leaf_block_verify
          : "off";
  result.leaf_mac =
      (options.leaf_mac == "audit" || options.leaf_mac == "on")
          ? options.leaf_mac
          : "off";
  result.leaf_pair_max_lookahead =
      std::clamp(options.leaf_pair_max_lookahead, 4, 128);
  result.interaction_mode =
      options.interaction_mode == "onthefly" ? "onthefly" : "cached";
  result.node_bounds_mode =
      options.node_bounds == "tight-aabb" ? std::string("tight-aabb")
                                          : std::string("moment");
  result.tree_builder = options.tree_builder == "topdown-aabb"
                            ? std::string("topdown-aabb")
                            : std::string("grouped-morton");
  result.expansion_center =
      (options.expansion_center == "area-centroid" ||
       options.expansion_center == "aabb-center")
          ? options.expansion_center
          : std::string("legacy");
  result.fragmentation_safety =
      (options.fragmentation_safety == "audit" ||
       options.fragmentation_safety == "split")
          ? options.fragmentation_safety
          : std::string("off");
  result.mac_safety =
      (options.mac_safety == "audit" || options.mac_safety == "physical-gap")
          ? options.mac_safety
          : std::string("off");
  result.mac_gap_factor =
      std::clamp(std::isfinite(options.mac_gap_factor)
                     ? options.mac_gap_factor
                     : 1.5,
                 0.5, 4.0);
  result.moment_audit =
      (options.moment_audit == "sampled" ||
       options.moment_audit == "full-small")
          ? options.moment_audit
          : std::string("off");
  result.moment_construction =
      options.moment_construction == "direct-node"
          ? std::string("direct-node")
          : std::string("translated");
  result.moment_sum =
      options.moment_sum == "kahan" ? std::string("kahan")
                                    : std::string("naive");
  result.convergence_audit =
      (options.convergence_audit == "sampled" ||
       options.convergence_audit == "full-small")
          ? options.convergence_audit
          : std::string("off");
  result.convergence_safety =
      (options.convergence_safety == "audit" ||
       options.convergence_safety == "split")
          ? options.convergence_safety
          : std::string("off");
  result.convergence_rho_cap =
      std::clamp(std::isfinite(options.convergence_rho_cap)
                     ? options.convergence_rho_cap
                     : 0.0,
                 0.0, 2.0);
  result.mac_aspect_cap =
      std::clamp(std::isfinite(options.mac_aspect_cap)
                     ? options.mac_aspect_cap
                     : 0.0,
                 0.0, 50.0);
  result.convergence_estimator =
      (options.convergence_estimator == "moment-tail" ||
       options.convergence_estimator == "hybrid")
          ? options.convergence_estimator
          : std::string("off");
  result.convergence_tail_threshold =
      std::clamp(std::isfinite(options.convergence_tail_threshold)
                     ? options.convergence_tail_threshold
                     : 1.0e4,
                 1.0, 1.0e12);
  result.convergence_estimator_error_cap =
      std::clamp(std::isfinite(options.convergence_estimator_error_cap)
                     ? options.convergence_estimator_error_cap
                     : 1.0e-3,
                 1.0e-8, 1.0);
  result.convergence_unstable_rho_cap =
      std::clamp(std::isfinite(options.convergence_unstable_rho_cap)
                     ? options.convergence_unstable_rho_cap
                     : 0.25,
                 0.05, 2.0);
  result.convergence_calibration =
      (options.convergence_calibration == "consecutive" ||
       options.convergence_calibration == "aspect" ||
       options.convergence_calibration == "consecutive-aspect")
          ? options.convergence_calibration
          : std::string("stage3j");
  result.convergence_dielectric_amp =
      options.convergence_dielectric_amp == "on" ? std::string("on")
                                                 : std::string("off");
  result.convergence_estimator_model =
      (options.convergence_estimator_model == "channel-tail" ||
       options.convergence_estimator_model == "channel-tail-eps")
          ? options.convergence_estimator_model
          : std::string("stage3l");
  result.convergence_channel_calibration =
      (options.convergence_channel_calibration == "k3-sqrt-eps" ||
       options.convergence_channel_calibration == "k3-eps" ||
       options.convergence_channel_calibration == "k4-sqrt-eps" ||
       options.convergence_channel_calibration == "k4-eps" ||
       options.convergence_channel_calibration ==
           "all-derivative-sqrt-eps" ||
       options.convergence_channel_calibration == "polynomial-eps")
          ? options.convergence_channel_calibration
          : std::string("stage3m");
  result.convergence_split_policy =
      options.convergence_split_policy == "k3-hot" ? std::string("k3-hot")
                                                   : std::string("global");
  result.convergence_split_selectivity =
      (options.convergence_split_selectivity == "absolute-aware" ||
       options.convergence_split_selectivity == "channel-absolute-aware")
          ? options.convergence_split_selectivity
          : std::string("legacy");
  result.convergence_split_absolute_floor = std::clamp(
      std::isfinite(options.convergence_split_absolute_floor)
          ? options.convergence_split_absolute_floor
          : 0.0,
      0.0, 1.0e-4);
  result.convergence_overreject_audit =
      (options.convergence_overreject_audit == "sampled" ||
       options.convergence_overreject_audit == "full-small")
          ? options.convergence_overreject_audit
          : std::string("off");
  result.convergence_overreject_detail =
      (options.convergence_overreject_detail == "summary" ||
       options.convergence_overreject_detail == "sampled-detail")
          ? options.convergence_overreject_detail
          : std::string("off");
  result.convergence_shape_audit =
      (options.convergence_shape_audit == "summary" ||
       options.convergence_shape_audit == "sampled-detail")
          ? options.convergence_shape_audit
          : std::string("off");
  result.convergence_kernel_response_audit =
      (options.convergence_kernel_response_audit == "summary" ||
       options.convergence_kernel_response_audit == "sampled-detail")
          ? options.convergence_kernel_response_audit
          : std::string("off");
  result.channel_missed_semantics =
      options.channel_missed_semantics == "accepted-only"
          ? std::string("accepted-only")
          : std::string("legacy");
  result.accepted_miss_diagnostics =
      (options.accepted_miss_diagnostics == "summary" ||
       options.accepted_miss_diagnostics == "sampled-detail")
          ? options.accepted_miss_diagnostics
          : std::string("off");
  result.accepted_miss_absolute_floor = std::clamp(
      std::isfinite(options.accepted_miss_absolute_floor)
          ? options.accepted_miss_absolute_floor
          : 1.0e-7,
      0.0, 1.0e-4);
  result.convergence_k3_accepted_rho_cap = std::clamp(
      std::isfinite(options.convergence_k3_accepted_rho_cap)
          ? options.convergence_k3_accepted_rho_cap
          : 0.0,
      0.0, 2.0);
  result.mac_distance_policy =
      options.mac_distance_policy == "surface-margin"
          ? std::string("surface-margin")
          : std::string("center");
  result.mac_surface_margin_factor = std::clamp(
      std::isfinite(options.mac_surface_margin_factor)
          ? options.mac_surface_margin_factor
          : 0.25,
      0.0, 0.75);
  result.mac_surface_margin_policy =
      (options.mac_surface_margin_policy == "adaptive" ||
       options.mac_surface_margin_policy == "risk-scaled")
          ? options.mac_surface_margin_policy
          : "all";
  result.mac_surface_margin_aspect_threshold = std::clamp(
      std::isfinite(options.mac_surface_margin_aspect_threshold)
          ? options.mac_surface_margin_aspect_threshold
          : 3.0,
      1.0, 50.0);
  result.mac_surface_margin_rho_trigger = std::clamp(
      std::isfinite(options.mac_surface_margin_rho_trigger)
          ? options.mac_surface_margin_rho_trigger
          : 0.55,
      0.0, 2.0);
  result.mac_surface_margin_soft_aspect_threshold = std::clamp(
      std::isfinite(options.mac_surface_margin_soft_aspect_threshold)
          ? options.mac_surface_margin_soft_aspect_threshold
          : 0.0,
      0.0, 50.0);
  result.mac_surface_margin_soft_rho_trigger = std::clamp(
      std::isfinite(options.mac_surface_margin_soft_rho_trigger)
          ? options.mac_surface_margin_soft_rho_trigger
          : 0.0,
      0.0, 2.0);
  result.mac_surface_margin_soft_factor = std::clamp(
      std::isfinite(options.mac_surface_margin_soft_factor)
          ? options.mac_surface_margin_soft_factor
          : 0.10,
      0.0, 0.75);
  result.expansion_audit =
      (options.expansion_audit == "sampled" ||
       options.expansion_audit == "full-small")
          ? options.expansion_audit
          : std::string("off");
  result.expansion_audit_channel =
      (options.expansion_audit_channel == "summary" ||
       options.expansion_audit_channel == "sampled-detail")
          ? options.expansion_audit_channel
          : std::string("off");
  result.quad_buffer_auto_candidate =
      (options.quad_buffer == "gauss7" || options.quad_buffer == "adaptive")
          ? std::string("explicit")
          : std::string("off");
  result.derivative_mode =
      options.derivative_mode.empty() ? std::string("analytic")
                                      : options.derivative_mode;
  result.surface_generator = options.surface_geometry == "point-panel"
                                 ? std::string("titan_surface_point_panel")
                                 : std::string("titan_sphere_patch");
  result.discretization =
      formulation.formulation == "tabi2-nodepatch" ? std::string("node_patch")
                                                    : std::string("centroid");
  result.rhs_scale = formulation.rhs_scale;
  result.coulomb_factor_kcal_mol = kCoulombKcalPerMol;
  result.total_charge = std::accumulate(
      atoms.begin(), atoms.end(), 0.0,
      [](double sum, const Atom &atom) { return sum + atom.q; });

  LongRunProgressWriter progress(options, executor.capacity());
  progress.write("surface_start", result, nullptr, 0, 0, 0, 0.0, "running",
                 true);
  const std::vector<Panel> panels =
      build_surface(atoms, options, &result, &executor);
  result.panel_count = panels.size();
  result.boundary_nodes = panels.size();
  const int surface_point_max = options.highres_surface ? 65536 : 4096;
  result.surface_point_count_requested =
      std::clamp(options.surface_point_count, 8, surface_point_max);
  result.panel_density_per_a2 =
      result.surface_area_a2 > 1.0e-12
          ? static_cast<double>(result.panel_count) / result.surface_area_a2
          : 0.0;
  if (result.panel_count == 0) {
    result.panel_budget_status = "empty";
  } else if (result.panel_density_per_a2 < 1.5) {
    result.panel_budget_status = "coarse";
  } else if (result.panel_density_per_a2 <= 3.0) {
    result.panel_budget_status = "target_density";
  } else {
    result.panel_budget_status = "dense";
  }
  progress.write("surface_end", result, nullptr, 0, 0, 0, 0.0, "running",
                 true);
  result.taylor_order = std::clamp(options.taylor_order, 3, 4);
  result.near_integration =
      resolve_near_integration_mode(options.near_integration, panels);
  result.quad_buffer = resolve_quad_buffer_mode(
      options.quad_buffer, result.near_integration, panels);
  if (options.quad_buffer == "auto")
    result.quad_buffer_auto_candidate = "off";
  result.memory_bytes =
      panels.size() * sizeof(Panel) + panels.size() * sizeof(double) * 8;

  if (panels.empty()) {
    result.converged = true;
    result.status = "ok";
    result.reaction_energy_source = "tabi_boundary_integral";
    result.convergence_reason = "converged";
    result.sum_mode_effective = "direct";
    progress.write("solve_done", result, nullptr, 0, 0, 0, 0.0, result.status,
                   true);
    return result;
  }

  const double kappa = debye_kappa_inv_a(options.salt_m,
                                         options.temperature_k);
  const double contrast =
      (options.eps_out - options.eps_in) /
      std::max(1.0e-12, options.eps_out + options.eps_in);
  const double operator_contrast = contrast * formulation.operator_contrast_scale;
  result.kappa_inv_a = kappa;
  result.dielectric_contrast = contrast;
  result.operator_contrast = operator_contrast;
  std::string mode = options.sum_mode;
  const bool paper_bie_formulation =
      formulation.formulation == "jcp2013-centroid" ||
      formulation.formulation == "jcp2013-treecode" ||
      formulation.formulation == "tabi2-nodepatch";
  if (mode == "auto") {
    mode = paper_bie_formulation
               ? "direct"
               : (result.surface_geometry == "triangle"
               ? "direct"
               : (panels.size() > 1024 ? "treecode" : "direct"));
  }
  result.sum_mode_effective = mode;
  result.matvec_mode = mode;

  if (paper_bie_formulation) {
    if (mode == "treecode" && formulation.formulation != "jcp2013-treecode") {
      throw std::runtime_error(
          "PBE: -pbe=tabi JCP/TABI2 direct formulations do not support "
          "-pbe-tabi-sum=treecode yet; use auto/direct, pin "
          "-pbe-tabi-formulation=jcp2013-treecode for the experimental "
          "paper-treecode validation path, or pin titan-v1 for the legacy "
          "scalar treecode.");
    }
    const bool experimental_treecode =
        formulation.formulation == "jcp2013-treecode" && mode == "treecode";
    result.matvec_mode = experimental_treecode ? "jcp2013_treecode" : "direct";
    result.sum_mode_effective = experimental_treecode ? "treecode" : "direct";
    result.taylor_order =
        experimental_treecode ? std::clamp(options.taylor_order, 3, 4) : 0;
    if (result.preconditioner_requested == "none")
      result.preconditioner_effective = "none";
    else if (result.preconditioner_requested == "block" ||
             result.preconditioner_requested == "block-budgeted")
      result.preconditioner_effective = "diagonal_block_pending";
    else if (result.preconditioner_requested == "coarse-galerkin-probe")
      result.preconditioner_effective = "diagonal_coarse_pending";
    else
      result.preconditioner_effective = "diagonal";
    const double eps_ratio =
        std::max(1.0e-12, options.eps_out) /
        std::max(1.0e-12, options.eps_in);
    const double lhs_phi = 0.5 * (1.0 + eps_ratio);
    const double lhs_phin = 0.5 * (1.0 + 1.0 / eps_ratio);
    const bool coupled = kappa > 1.0e-14;
    result.operator_contrast = eps_ratio;
    result.unknown_count = coupled ? 2 * panels.size() : panels.size();
    result.bie_system = coupled ? "linear_pb_2x2" : "poisson";
    result.lhs_diag_phi = lhs_phi;
    result.lhs_diag_phin = coupled ? lhs_phin : 0.0;
    result.energy_equation = coupled ? "jcp2013_eq14" : "jcp2013_eq22";
    result.kernel_convention =
        "G0_Gk_source_normal_target_normal_outward_" +
        result.near_integration + "_near";

    const auto rhs_begin = std::chrono::steady_clock::now();
    double s1_norm = 0.0;
    double s2_norm = 0.0;
    const std::vector<double> rhs = build_jcp_rhs(
        atoms, panels, options.eps_in, coupled, &s1_norm, &s2_norm);
    result.rhs_ms = elapsed_ms(rhs_begin, std::chrono::steady_clock::now());
    result.source_norm_s1 = s1_norm;
    result.source_norm_s2 = s2_norm;

    const int size_factor = coupled ? 8 : 4;
    const int max_iter =
        std::min(std::max(1, options.max_iter),
                 std::max(64, static_cast<int>(panels.size()) * size_factor));
    const double tol = std::max(1.0e-12, options.residual_tol);
    std::vector<double> diagonal(result.unknown_count, 1.0);
    if (result.preconditioner_effective != "none") {
      std::fill(diagonal.begin(), diagonal.end(), lhs_phi);
      if (coupled) {
        for (std::size_t i = panels.size(); i < diagonal.size(); ++i)
          diagonal[i] = lhs_phin;
      }
    }

    double matvec_ms_total = 0.0;
    std::size_t matvec_calls = 0;
    TreeWorkspace jcp_tree;
    TreecodeStats jcp_tree_stats;
    jcp_tree_stats.near_integration = result.near_integration;
    jcp_tree_stats.quad_buffer = result.quad_buffer;
    jcp_tree_stats.interaction_mode = result.interaction_mode;
    jcp_tree_stats.node_bounds_mode = result.node_bounds_mode;
    jcp_tree_stats.tree_builder = result.tree_builder;
    jcp_tree_stats.expansion_center = result.expansion_center;
    jcp_tree_stats.fragmentation_safety = result.fragmentation_safety;
    jcp_tree_stats.mac_safety = result.mac_safety;
    jcp_tree_stats.mac_gap_factor = result.mac_gap_factor;
    jcp_tree_stats.moment_audit = result.moment_audit;
    jcp_tree_stats.moment_construction = result.moment_construction;
    jcp_tree_stats.moment_sum = result.moment_sum;
    jcp_tree_stats.convergence_audit = result.convergence_audit;
    jcp_tree_stats.convergence_safety = result.convergence_safety;
    jcp_tree_stats.convergence_rho_cap = result.convergence_rho_cap;
    jcp_tree_stats.mac_aspect_cap = result.mac_aspect_cap;
    jcp_tree_stats.convergence_estimator = result.convergence_estimator;
    jcp_tree_stats.convergence_tail_threshold =
        result.convergence_tail_threshold;
    jcp_tree_stats.convergence_estimator_error_cap =
        result.convergence_estimator_error_cap;
    jcp_tree_stats.convergence_unstable_rho_cap =
        result.convergence_unstable_rho_cap;
    jcp_tree_stats.convergence_calibration =
        result.convergence_calibration;
    jcp_tree_stats.convergence_dielectric_amp =
        result.convergence_dielectric_amp;
    jcp_tree_stats.convergence_estimator_model =
        result.convergence_estimator_model;
    jcp_tree_stats.convergence_channel_calibration =
        result.convergence_channel_calibration;
    jcp_tree_stats.convergence_split_policy =
        result.convergence_split_policy;
    jcp_tree_stats.convergence_split_selectivity =
        result.convergence_split_selectivity;
    jcp_tree_stats.convergence_split_absolute_floor =
        result.convergence_split_absolute_floor;
    jcp_tree_stats.convergence_overreject_audit =
        result.convergence_overreject_audit;
    jcp_tree_stats.convergence_overreject_detail =
        result.convergence_overreject_detail;
    jcp_tree_stats.convergence_shape_audit =
        result.convergence_shape_audit;
    jcp_tree_stats.convergence_kernel_response_audit =
        result.convergence_kernel_response_audit;
    jcp_tree_stats.channel_missed_semantics =
        result.channel_missed_semantics;
    jcp_tree_stats.accepted_miss_diagnostics =
        result.accepted_miss_diagnostics;
    jcp_tree_stats.accepted_miss_absolute_floor =
        result.accepted_miss_absolute_floor;
    jcp_tree_stats.convergence_k3_accepted_rho_cap =
        result.convergence_k3_accepted_rho_cap;
    jcp_tree_stats.mac_distance_policy = result.mac_distance_policy;
    jcp_tree_stats.mac_surface_margin_factor =
        result.mac_surface_margin_factor;
    jcp_tree_stats.mac_surface_margin_policy =
        result.mac_surface_margin_policy;
    jcp_tree_stats.mac_surface_margin_aspect_threshold =
        result.mac_surface_margin_aspect_threshold;
    jcp_tree_stats.mac_surface_margin_rho_trigger =
        result.mac_surface_margin_rho_trigger;
    jcp_tree_stats.mac_surface_margin_soft_aspect_threshold =
        result.mac_surface_margin_soft_aspect_threshold;
    jcp_tree_stats.mac_surface_margin_soft_rho_trigger =
        result.mac_surface_margin_soft_rho_trigger;
    jcp_tree_stats.mac_surface_margin_soft_factor =
        result.mac_surface_margin_soft_factor;
    jcp_tree_stats.cluster_policy = result.cluster_policy;
    jcp_tree_stats.target_batch_size = result.target_batch_size;
    jcp_tree_stats.leaf_pair_min_source_panels =
        result.leaf_pair_min_source_panels;
    jcp_tree_stats.leaf_pair_min_targets = result.leaf_pair_min_targets;
    jcp_tree_stats.leaf_kernel_profile = result.leaf_kernel_profile;
    jcp_tree_stats.leaf_work_profile = result.leaf_work_profile;
    jcp_tree_stats.matvec_cost_profile = result.matvec_cost_profile;
    jcp_tree_stats.jcp_direct_kernel = result.jcp_direct_kernel;
    jcp_tree_stats.leaf_pair_selection_profile =
        result.leaf_pair_selection_profile;
    jcp_tree_stats.leaf_pair_probe_cache = result.leaf_pair_probe_cache;
    jcp_tree_stats.leaf_pair_target_group = result.leaf_pair_target_group;
    jcp_tree_stats.leaf_pair_directness_probe =
        result.leaf_pair_directness_probe;
    jcp_tree_stats.leaf_pair_directness_cache =
        result.leaf_pair_directness_cache;
    jcp_tree_stats.leaf_pair_directness_verify =
        result.leaf_pair_directness_verify;
    jcp_tree_stats.leaf_pair_scheduler = result.leaf_pair_scheduler;
    jcp_tree_stats.leaf_pair_worklist_max_active =
        result.leaf_pair_worklist_max_active;
    jcp_tree_stats.leaf_pair_worklist_verify =
        result.leaf_pair_worklist_verify;
    jcp_tree_stats.leaf_block_verify = result.leaf_block_verify;
    jcp_tree_stats.leaf_pair_max_lookahead =
        result.leaf_pair_max_lookahead;
    jcp_tree_stats.expansion_audit = result.expansion_audit;
    jcp_tree_stats.expansion_audit_channel = result.expansion_audit_channel;
    const bool fmm_m2l_single_matvec =
        experimental_treecode &&
        options.fmm_m2l == "single-matvec-analytic";
    const bool preconditioner_needs_jcp_tree =
        result.preconditioner_requested == "block" ||
        result.preconditioner_requested == "block-budgeted" ||
        result.preconditioner_requested == "coarse-galerkin-probe";
    const bool needs_jcp_tree = experimental_treecode ||
                                preconditioner_needs_jcp_tree;
    if (needs_jcp_tree) {
      progress.write("tree_build_start", result, &jcp_tree_stats,
                     matvec_calls, 0, 0, 0.0, "running", true);
      const std::string jcp_tree_interaction_mode =
          fmm_m2l_single_matvec && !preconditioner_needs_jcp_tree
              ? std::string("onthefly")
              : (experimental_treecode ? options.interaction_mode
                                       : std::string("onthefly"));
      jcp_tree = build_tree_workspace(
          panels, executor, static_cast<std::size_t>(
                                std::max(8, options.tree_leaf_size)),
          options.tree_opening_theta, options.leaf_near_policy,
          options.leaf_panel_kernel, options.adaptive_level,
          result.quad_buffer, jcp_tree_interaction_mode,
          result.node_bounds_mode, result.tree_builder, options.mac_safety,
          options.mac_gap_factor, options.moment_audit,
          result.expansion_center, result.fragmentation_safety,
          options.expansion_audit, options.expansion_audit_sample_count,
          options.expansion_audit_max_panels,
          options.moment_construction, options.moment_sum,
          options.convergence_audit, options.convergence_safety,
          options.convergence_rho_cap, options.mac_aspect_cap,
          options.convergence_estimator,
          options.convergence_tail_threshold,
          options.convergence_estimator_error_cap,
          options.convergence_unstable_rho_cap,
          options.convergence_calibration,
          options.cluster_policy,
          options.convergence_dielectric_amp, options.target_batch_size,
          options.leaf_pair_min_source_panels,
          options.leaf_pair_min_targets,
          options.convergence_estimator_model,
          options.expansion_audit_channel,
          options.convergence_channel_calibration,
          options.channel_missed_semantics,
          options.convergence_split_policy,
          options.accepted_miss_diagnostics,
          options.convergence_k3_accepted_rho_cap,
          options.mac_distance_policy,
          options.mac_surface_margin_factor,
          options.mac_surface_margin_policy,
          options.mac_surface_margin_aspect_threshold,
          options.mac_surface_margin_rho_trigger,
          options.accepted_miss_absolute_floor, options.leaf_pair_kernel,
          options.leaf_kernel_profile, options.leaf_work_profile,
          options.jcp_direct_kernel,
          options.leaf_pair_selection_profile,
          options.leaf_pair_probe_cache,
          options.leaf_pair_target_group,
          options.leaf_pair_max_lookahead,
          options.leaf_pair_directness_probe,
            options.leaf_pair_directness_cache,
            options.leaf_pair_directness_verify,
            options.leaf_pair_scheduler,
            options.leaf_pair_worklist_max_active,
            options.leaf_pair_worklist_verify, options.leaf_block_verify,
            options.mac_surface_margin_soft_aspect_threshold,
            options.mac_surface_margin_soft_rho_trigger,
            options.mac_surface_margin_soft_factor,
            options.convergence_split_selectivity,
            options.convergence_split_absolute_floor,
            options.convergence_overreject_audit,
            options.convergence_overreject_detail,
            options.convergence_shape_audit,
            options.convergence_kernel_response_audit,
            result.leaf_mac, options.fmm_m2l);
      jcp_tree.leaf_mac = result.leaf_mac;
      jcp_tree.taylor_order = result.taylor_order;
      jcp_tree.derivative_mode = result.derivative_mode;
      result.interaction_mode = jcp_tree.interaction_mode;
      result.node_bounds_mode = jcp_tree.node_bounds_mode;
      result.tree_builder = jcp_tree.tree_builder;
      result.expansion_center = jcp_tree.expansion_center;
      result.fragmentation_safety = jcp_tree.fragmentation_safety;
      result.mac_safety = jcp_tree.mac_safety;
      result.mac_gap_factor = jcp_tree.mac_gap_factor;
      result.moment_audit = jcp_tree.moment_audit;
      result.moment_construction = jcp_tree.moment_construction;
      result.moment_sum = jcp_tree.moment_sum;
      result.convergence_audit = jcp_tree.convergence_audit;
      result.convergence_safety = jcp_tree.convergence_safety;
      result.convergence_rho_cap = jcp_tree.convergence_rho_cap;
      result.mac_aspect_cap = jcp_tree.mac_aspect_cap;
      result.convergence_estimator = jcp_tree.convergence_estimator;
      result.convergence_tail_threshold =
          jcp_tree.convergence_tail_threshold;
      result.convergence_estimator_error_cap =
          jcp_tree.convergence_estimator_error_cap;
      result.convergence_unstable_rho_cap =
          jcp_tree.convergence_unstable_rho_cap;
      result.convergence_calibration = jcp_tree.convergence_calibration;
      result.convergence_dielectric_amp =
          jcp_tree.convergence_dielectric_amp;
      result.convergence_estimator_model =
          jcp_tree.convergence_estimator_model;
      result.convergence_channel_calibration =
          jcp_tree.convergence_channel_calibration;
      result.convergence_split_policy = jcp_tree.convergence_split_policy;
      result.convergence_split_selectivity =
          jcp_tree.convergence_split_selectivity;
      result.convergence_split_absolute_floor =
          jcp_tree.convergence_split_absolute_floor;
      result.convergence_overreject_audit =
          jcp_tree.convergence_overreject_audit;
      result.convergence_overreject_detail =
          jcp_tree.convergence_overreject_detail;
      result.convergence_shape_audit = jcp_tree.convergence_shape_audit;
      result.convergence_kernel_response_audit =
          jcp_tree.convergence_kernel_response_audit;
      result.channel_missed_semantics = jcp_tree.channel_missed_semantics;
      result.accepted_miss_diagnostics =
          jcp_tree.accepted_miss_diagnostics;
      result.accepted_miss_absolute_floor =
          jcp_tree.accepted_miss_absolute_floor;
      result.convergence_k3_accepted_rho_cap =
          jcp_tree.convergence_k3_accepted_rho_cap;
      result.mac_distance_policy = jcp_tree.mac_distance_policy;
      result.mac_surface_margin_factor =
          jcp_tree.mac_surface_margin_factor;
      result.mac_surface_margin_policy =
          jcp_tree.mac_surface_margin_policy;
      result.mac_surface_margin_aspect_threshold =
          jcp_tree.mac_surface_margin_aspect_threshold;
    result.mac_surface_margin_rho_trigger =
        jcp_tree.mac_surface_margin_rho_trigger;
    result.mac_surface_margin_soft_aspect_threshold =
        jcp_tree.mac_surface_margin_soft_aspect_threshold;
    result.mac_surface_margin_soft_rho_trigger =
        jcp_tree.mac_surface_margin_soft_rho_trigger;
    result.mac_surface_margin_soft_factor =
        jcp_tree.mac_surface_margin_soft_factor;
      result.cluster_policy = jcp_tree.cluster_policy;
      result.target_batch_size =
          static_cast<int>(jcp_tree.target_batch_size);
      result.leaf_pair_min_source_panels =
          static_cast<int>(jcp_tree.leaf_pair_min_source_panels);
      result.leaf_pair_min_targets =
          static_cast<int>(jcp_tree.leaf_pair_min_targets);
      result.leaf_kernel_profile = jcp_tree.leaf_kernel_profile;
      result.leaf_work_profile = jcp_tree.leaf_work_profile;
      result.jcp_direct_kernel = jcp_tree.jcp_direct_kernel;
      result.leaf_pair_selection_profile =
          jcp_tree.leaf_pair_selection_profile;
      result.leaf_pair_probe_cache = jcp_tree.leaf_pair_probe_cache;
      result.leaf_pair_target_group = jcp_tree.leaf_pair_target_group;
      result.leaf_pair_directness_probe =
          jcp_tree.leaf_pair_directness_probe;
      result.leaf_pair_directness_cache =
          jcp_tree.leaf_pair_directness_cache;
      result.leaf_pair_directness_verify =
          jcp_tree.leaf_pair_directness_verify;
      result.leaf_pair_scheduler = jcp_tree.leaf_pair_scheduler;
      result.leaf_pair_worklist_max_active =
          static_cast<int>(jcp_tree.leaf_pair_worklist_max_active);
      result.leaf_pair_worklist_verify = jcp_tree.leaf_pair_worklist_verify;
      result.leaf_block_verify = jcp_tree.leaf_block_verify;
      result.leaf_pair_max_lookahead =
          static_cast<int>(jcp_tree.leaf_pair_max_lookahead);
      result.expansion_audit = jcp_tree.expansion_audit;
      result.expansion_audit_channel = jcp_tree.expansion_audit_channel;
      jcp_tree_stats.interaction_mode = jcp_tree.interaction_mode;
      jcp_tree_stats.node_bounds_mode = jcp_tree.node_bounds_mode;
      jcp_tree_stats.tree_builder = jcp_tree.tree_builder;
      jcp_tree_stats.expansion_center = jcp_tree.expansion_center;
      jcp_tree_stats.fragmentation_safety = jcp_tree.fragmentation_safety;
      jcp_tree_stats.mac_safety = jcp_tree.mac_safety;
      jcp_tree_stats.mac_gap_factor = jcp_tree.mac_gap_factor;
      jcp_tree_stats.moment_audit = jcp_tree.moment_audit;
      jcp_tree_stats.moment_construction = jcp_tree.moment_construction;
      jcp_tree_stats.moment_sum = jcp_tree.moment_sum;
      jcp_tree_stats.convergence_audit = jcp_tree.convergence_audit;
      jcp_tree_stats.convergence_safety = jcp_tree.convergence_safety;
      jcp_tree_stats.convergence_rho_cap = jcp_tree.convergence_rho_cap;
      jcp_tree_stats.mac_aspect_cap = jcp_tree.mac_aspect_cap;
      jcp_tree_stats.convergence_estimator =
          jcp_tree.convergence_estimator;
      jcp_tree_stats.convergence_tail_threshold =
          jcp_tree.convergence_tail_threshold;
      jcp_tree_stats.convergence_estimator_error_cap =
          jcp_tree.convergence_estimator_error_cap;
      jcp_tree_stats.convergence_unstable_rho_cap =
          jcp_tree.convergence_unstable_rho_cap;
      jcp_tree_stats.convergence_calibration =
          jcp_tree.convergence_calibration;
      jcp_tree_stats.convergence_dielectric_amp =
          jcp_tree.convergence_dielectric_amp;
      jcp_tree_stats.convergence_estimator_model =
          jcp_tree.convergence_estimator_model;
      jcp_tree_stats.convergence_channel_calibration =
          jcp_tree.convergence_channel_calibration;
      jcp_tree_stats.convergence_split_policy =
          jcp_tree.convergence_split_policy;
      jcp_tree_stats.convergence_split_selectivity =
          jcp_tree.convergence_split_selectivity;
      jcp_tree_stats.convergence_split_absolute_floor =
          jcp_tree.convergence_split_absolute_floor;
      jcp_tree_stats.convergence_overreject_audit =
          jcp_tree.convergence_overreject_audit;
      jcp_tree_stats.convergence_overreject_detail =
          jcp_tree.convergence_overreject_detail;
      jcp_tree_stats.convergence_shape_audit =
          jcp_tree.convergence_shape_audit;
      jcp_tree_stats.convergence_kernel_response_audit =
          jcp_tree.convergence_kernel_response_audit;
      jcp_tree_stats.channel_missed_semantics =
          jcp_tree.channel_missed_semantics;
      jcp_tree_stats.accepted_miss_diagnostics =
          jcp_tree.accepted_miss_diagnostics;
      jcp_tree_stats.accepted_miss_absolute_floor =
          jcp_tree.accepted_miss_absolute_floor;
      jcp_tree_stats.convergence_k3_accepted_rho_cap =
          jcp_tree.convergence_k3_accepted_rho_cap;
      jcp_tree_stats.mac_distance_policy = jcp_tree.mac_distance_policy;
      jcp_tree_stats.mac_surface_margin_factor =
          jcp_tree.mac_surface_margin_factor;
      jcp_tree_stats.mac_surface_margin_policy =
          jcp_tree.mac_surface_margin_policy;
      jcp_tree_stats.mac_surface_margin_aspect_threshold =
          jcp_tree.mac_surface_margin_aspect_threshold;
      jcp_tree_stats.mac_surface_margin_rho_trigger =
          jcp_tree.mac_surface_margin_rho_trigger;
      jcp_tree_stats.mac_surface_margin_soft_aspect_threshold =
          jcp_tree.mac_surface_margin_soft_aspect_threshold;
      jcp_tree_stats.mac_surface_margin_soft_rho_trigger =
          jcp_tree.mac_surface_margin_soft_rho_trigger;
      jcp_tree_stats.mac_surface_margin_soft_factor =
          jcp_tree.mac_surface_margin_soft_factor;
      jcp_tree_stats.cluster_policy = jcp_tree.cluster_policy;
      jcp_tree_stats.target_batch_size =
          static_cast<int>(jcp_tree.target_batch_size);
      jcp_tree_stats.leaf_pair_min_source_panels =
          static_cast<int>(jcp_tree.leaf_pair_min_source_panels);
      jcp_tree_stats.leaf_pair_min_targets =
          static_cast<int>(jcp_tree.leaf_pair_min_targets);
      jcp_tree_stats.leaf_kernel_profile = jcp_tree.leaf_kernel_profile;
      jcp_tree_stats.leaf_work_profile = jcp_tree.leaf_work_profile;
      jcp_tree_stats.matvec_cost_profile = result.matvec_cost_profile;
      jcp_tree_stats.jcp_direct_kernel = jcp_tree.jcp_direct_kernel;
      jcp_tree_stats.leaf_pair_selection_profile =
          jcp_tree.leaf_pair_selection_profile;
      jcp_tree_stats.leaf_pair_probe_cache = jcp_tree.leaf_pair_probe_cache;
      jcp_tree_stats.leaf_pair_target_group = jcp_tree.leaf_pair_target_group;
      jcp_tree_stats.leaf_pair_directness_probe =
          jcp_tree.leaf_pair_directness_probe;
      jcp_tree_stats.leaf_pair_directness_cache =
          jcp_tree.leaf_pair_directness_cache;
      jcp_tree_stats.leaf_pair_directness_verify =
          jcp_tree.leaf_pair_directness_verify;
      jcp_tree_stats.leaf_pair_scheduler = jcp_tree.leaf_pair_scheduler;
      jcp_tree_stats.leaf_pair_worklist_max_active =
          static_cast<int>(jcp_tree.leaf_pair_worklist_max_active);
      jcp_tree_stats.leaf_pair_worklist_verify =
          jcp_tree.leaf_pair_worklist_verify;
      jcp_tree_stats.leaf_block_verify = jcp_tree.leaf_block_verify;
      jcp_tree_stats.leaf_pair_max_lookahead =
          static_cast<int>(jcp_tree.leaf_pair_max_lookahead);
      jcp_tree_stats.expansion_audit = jcp_tree.expansion_audit;
      jcp_tree_stats.expansion_audit_channel =
          jcp_tree.expansion_audit_channel;
      jcp_tree_stats.tree_depth_max = jcp_tree.depth_max;
      jcp_tree_stats.node_radius_max = jcp_tree.node_radius_max;
      jcp_tree_stats.node_radius_mean = jcp_tree.node_radius_mean;
      jcp_tree_stats.node_aspect_ratio_max = jcp_tree.node_aspect_ratio_max;
      jcp_tree_stats.node_aspect_ratio_mean = jcp_tree.node_aspect_ratio_mean;
      result.tree_build_ms = jcp_tree.build_ms;
      result.tree_nodes = jcp_tree.nodes.size();
      result.tree_leaf_count = jcp_tree.leaf_count;
      result.tree_depth_max = jcp_tree.depth_max;
      result.memory_bytes += jcp_tree.nodes.size() * sizeof(TreeNode) +
                             jcp_tree.child_indices.size() *
                                 sizeof(std::size_t) +
                             jcp_tree.order.size() * sizeof(std::size_t) +
                             jcp_tree.morton.size() * sizeof(std::uint64_t) +
                             panels.size() * 14u * sizeof(double) +
                             panels.size() * sizeof(std::size_t);
      jcp_tree_stats.tree_build_ms = result.tree_build_ms;
      jcp_tree_stats.tree_nodes = result.tree_nodes;
      jcp_tree_stats.tree_leaf_count = result.tree_leaf_count;
      progress.write("tree_build_end", result, &jcp_tree_stats, matvec_calls,
                     0, 0, 0.0, "running", true);
    }
    auto matvec = [&](const std::vector<double> &x, std::vector<double> &out) {
      const auto begin = std::chrono::steady_clock::now();
      if (experimental_treecode) {
        apply_jcp_treecode_operator(panels, x, eps_ratio, kappa, lhs_phi,
                                    lhs_phin, coupled, jcp_tree,
                                    result.near_integration,
                                    result.quad_buffer, executor, out,
                                    &jcp_tree_stats);
      } else if (coupled) {
        apply_jcp_linear_pb_operator(panels, x, eps_ratio, kappa, lhs_phi,
                                     lhs_phin, result.near_integration,
                                     result.quad_buffer, executor, out,
                                     &jcp_tree_stats);
      } else {
        apply_jcp_poisson_operator(panels, x, eps_ratio, lhs_phi,
                                   result.near_integration,
                                   result.quad_buffer, executor, out,
                                   &jcp_tree_stats);
      }
      matvec_ms_total += elapsed_ms(begin, std::chrono::steady_clock::now());
      ++matvec_calls;
    };
    std::vector<std::size_t> direct_order(panels.size(), 0);
    std::iota(direct_order.begin(), direct_order.end(), std::size_t{0});
    auto build_audit_rows = [](std::size_t row_count, std::string mode,
                               int requested) {
      std::vector<std::size_t> rows;
      if (row_count == 0)
        return rows;
      const bool full = mode == "full-small" && row_count <= 4096u;
      const std::size_t m =
          full ? row_count
               : std::min<std::size_t>(
                     row_count,
                     static_cast<std::size_t>(
                         std::clamp(requested, 128, 4096)));
      rows.reserve(m);
      if (full) {
        for (std::size_t row = 0; row < row_count; ++row)
          rows.push_back(row);
        return rows;
      }
      std::size_t previous = std::numeric_limits<std::size_t>::max();
      for (std::size_t k = 0; k < m; ++k) {
        std::size_t row =
            ((2u * k + 1u) * row_count) / std::max<std::size_t>(1, 2u * m);
        row = std::min(row, row_count - 1u);
        if (row != previous) {
          rows.push_back(row);
          previous = row;
        }
      }
      if (rows.empty())
        rows.push_back(0);
      return rows;
    };
    SampledTrueResidualCallback sampled_true_residual;
    if (!experimental_treecode) {
      sampled_true_residual =
          [&](const std::vector<double> &x,
              const std::vector<std::size_t> &rows) -> SampledResidualResult {
        SampledResidualResult sample;
        sample.rows = rows.size();
        if (rows.empty())
          return sample;
        const std::size_t unknown_count = rhs.size();
        const std::size_t panel_count = panels.size();
        const std::size_t workers = solver_workers(executor, rows.size());
        std::vector<double> partial_residual_sq(workers, 0.0);
        std::vector<double> partial_b_sq(workers, 0.0);
        std::atomic_bool finite{true};
        executor.parallel_for(
            rows.size(), workers, 1,
            [&](std::size_t worker, std::size_t begin, std::size_t end) {
              double local_residual_sq = 0.0;
              double local_b_sq = 0.0;
              for (std::size_t pos = begin; pos < end; ++pos) {
                const std::size_t row = rows[pos];
                if (row >= unknown_count)
                  continue;
                const std::size_t target =
                    coupled && row >= panel_count ? row - panel_count : row;
                if (target >= panel_count)
                  continue;
                const JcpNearAccum sum = accumulate_jcp_near_ordered(
                    panels[target], target, panels, x, eps_ratio, kappa,
                    coupled, direct_order, 0, panel_count,
                    result.near_integration, result.quad_buffer, nullptr);
                double ax_row = 0.0;
                if (coupled) {
                  if (row < panel_count) {
                    ax_row = lhs_phi * x[target] - sum.phi;
                  } else {
                    ax_row = lhs_phin * x[panel_count + target] - sum.phin;
                  }
                } else {
                  ax_row = lhs_phi * x[target] - sum.phi;
                }
                if (!std::isfinite(ax_row)) {
                  finite.store(false, std::memory_order_relaxed);
                  continue;
                }
                const double residual_value = rhs[row] - ax_row;
                local_residual_sq += residual_value * residual_value;
                local_b_sq += rhs[row] * rhs[row];
              }
              partial_residual_sq[worker] = local_residual_sq;
              partial_b_sq[worker] = local_b_sq;
            });
        const double residual_sq = std::accumulate(
            partial_residual_sq.begin(), partial_residual_sq.end(), 0.0);
        const double b_sq =
            std::accumulate(partial_b_sq.begin(), partial_b_sq.end(), 0.0);
        sample.abs = std::sqrt(std::max(0.0, residual_sq));
        const double denom = std::max(1.0e-300,
                                      std::sqrt(std::max(0.0, b_sq)));
        sample.rel = sample.abs / denom;
        sample.finite = finite.load(std::memory_order_relaxed) &&
                        std::isfinite(sample.rel);
        return sample;
      };
    }

    std::function<void(const std::vector<double> &, std::vector<double> &)>
        custom_preconditioner;
    std::shared_ptr<JcpBlockPreconditioner> block_preconditioner;
    std::shared_ptr<JcpSparseIluPreconditioner> sparse_ilu_preconditioner;
    std::shared_ptr<JcpCoarseGalerkinPreconditioner> coarse_preconditioner;
    double preconditioner_setup_ms = 0.0;
    constexpr std::size_t kBlockPanelPrebuildCap = 12000;
    const bool legacy_block = result.preconditioner_requested == "block";
    const bool budgeted_block =
        result.preconditioner_requested == "block-budgeted";
    const bool coarse_probe =
        result.preconditioner_requested == "coarse-galerkin-probe";
    const bool fmm_block_jacobi =
        result.preconditioner_requested == "fmm-block-jacobi";
    const bool fmm_overlap_schwarz =
        result.preconditioner_requested == "fmm-overlap-schwarz";
    const bool fmm_nearfield_ilu0 =
        result.preconditioner_requested == "fmm-nearfield-ilu0";
    auto block_effective_label = [&](bool budgeted,
                                     bool partial) -> std::string {
      std::string base = budgeted ? "block_budgeted" : "block";
      if (result.preconditioner_block_near == "buffer-aware")
        base += "_buffer_aware";
      else if (result.preconditioner_block_near == "operator-near")
        base += "_operator_near";
      if (partial)
        base += "_partial";
      return base;
    };
    if (fmm_block_jacobi || fmm_overlap_schwarz || fmm_nearfield_ilu0) {
      progress.write("precond_start", result, &jcp_tree_stats, matvec_calls,
                     0, 0, 0.0, "running", true);
      if (!experimental_treecode || !fmm_m2l_single_matvec || !coupled) {
        throw std::invalid_argument(
            "PBE TABI FMM-native preconditioners require coupled "
            "JCP treecode with -pbe-tabi-fmm-m2l=single-matvec-analytic");
      }
      result.preconditioner_block_near = "operator-near";
      const auto precond_begin = std::chrono::steady_clock::now();
      const auto fmm_context = ensure_jcp_fmm_m2l_context(
          panels, eps_ratio, kappa,
          static_cast<std::size_t>(std::max(8, options.tree_leaf_size)),
          options.tree_opening_theta, executor, &jcp_tree.fmm_m2l_context);
      if (fmm_nearfield_ilu0) {
        result.preconditioner_estimated_bytes =
            estimate_fmm_m2l_nearfield_ilu0_bytes(
                *fmm_context, diagonal.size(), coupled,
                static_cast<std::size_t>(
                    result.preconditioner_ilu_neighbor_cap));
      } else {
        result.preconditioner_estimated_bytes =
            fmm_overlap_schwarz
                ? estimate_fmm_m2l_overlap_schwarz_bytes(
                      *fmm_context, diagonal.size(), coupled,
                      static_cast<std::size_t>(
                          result.preconditioner_block_panel_cap),
                      static_cast<std::size_t>(
                          result.preconditioner_overlap_panel_cap))
                : estimate_fmm_m2l_block_jacobi_bytes(
                      *fmm_context, diagonal.size(), coupled,
                      static_cast<std::size_t>(
                          result.preconditioner_block_panel_cap));
      }
      const bool block_budget_cap =
          result.preconditioner_estimated_bytes >
          result.preconditioner_memory_budget_bytes;
      if (!block_budget_cap) {
        if (fmm_nearfield_ilu0) {
          sparse_ilu_preconditioner =
              build_fmm_m2l_nearfield_ilu0_preconditioner(
                  panels, *fmm_context, coupled, eps_ratio, kappa, lhs_phi,
                  lhs_phin, diagonal,
                  static_cast<std::size_t>(
                      result.preconditioner_ilu_neighbor_cap),
                  result.preconditioner_ilu_drop_tol, executor,
                  result.near_integration, result.quad_buffer);
        } else if (fmm_overlap_schwarz) {
          block_preconditioner =
              build_fmm_m2l_overlap_schwarz_preconditioner(
                  panels, *fmm_context, coupled, eps_ratio, kappa, lhs_phi,
                  lhs_phin, diagonal,
                  static_cast<std::size_t>(
                      result.preconditioner_block_panel_cap),
                  static_cast<std::size_t>(
                      result.preconditioner_overlap_panel_cap),
                  executor, result.near_integration, result.quad_buffer);
        } else {
          block_preconditioner = build_fmm_m2l_block_jacobi_preconditioner(
              panels, *fmm_context, coupled, eps_ratio, kappa, lhs_phi,
              lhs_phin, diagonal,
              static_cast<std::size_t>(result.preconditioner_block_panel_cap),
              executor, result.near_integration, result.quad_buffer);
        }
      }
      preconditioner_setup_ms =
          elapsed_ms(precond_begin, std::chrono::steady_clock::now());
      result.preconditioner_setup_ms = preconditioner_setup_ms;
      if (fmm_nearfield_ilu0) {
        result.preconditioner_valid_blocks =
            sparse_ilu_preconditioner != nullptr
                ? sparse_ilu_preconditioner->valid_blocks
                : 0;
        result.preconditioner_valid_panels =
            sparse_ilu_preconditioner != nullptr
                ? sparse_ilu_preconditioner->valid_panel_count
                : 0;
        result.preconditioner_ill_conditioned_blocks =
            sparse_ilu_preconditioner != nullptr
                ? sparse_ilu_preconditioner->ill_conditioned_blocks
                : 0;
        result.preconditioner_quality_rejected_blocks = 0;
        result.preconditioner_quality_probe_max_gain = 0.0;
        result.preconditioner_quality_status =
            sparse_ilu_preconditioner != nullptr
                ? sparse_ilu_preconditioner->quality_status
                : std::string("off");
        result.preconditioner_sparse_nnz_blocks =
            sparse_ilu_preconditioner != nullptr
                ? sparse_ilu_preconditioner->nnz_blocks
                : 0;
        result.preconditioner_sparse_avg_degree =
            sparse_ilu_preconditioner != nullptr
                ? sparse_ilu_preconditioner->avg_degree
                : 0.0;
        result.preconditioner_sparse_max_degree =
            sparse_ilu_preconditioner != nullptr
                ? sparse_ilu_preconditioner->max_degree
                : 0;
      } else {
        result.preconditioner_valid_blocks =
            block_preconditioner != nullptr ? block_preconditioner->valid_blocks
                                            : 0;
        result.preconditioner_valid_panels =
            block_preconditioner != nullptr
                ? block_preconditioner->valid_panel_count
                : 0;
        result.preconditioner_ill_conditioned_blocks =
            block_preconditioner != nullptr
                ? block_preconditioner->ill_conditioned_blocks
                : 0;
        result.preconditioner_quality_rejected_blocks =
            block_preconditioner != nullptr
                ? block_preconditioner->quality_rejected_blocks
                : 0;
        result.preconditioner_quality_probe_max_gain =
            block_preconditioner != nullptr
                ? block_preconditioner->quality_probe_max_gain
                : 0.0;
        result.preconditioner_quality_status =
            block_preconditioner != nullptr
                ? block_preconditioner->quality_status
                : std::string("off");
      }
      const bool block_ill_conditioned =
          block_preconditioner != nullptr &&
          block_preconditioner->ill_conditioned_blocks > 0;
      const bool ilu_ill_conditioned =
          sparse_ilu_preconditioner != nullptr &&
          sparse_ilu_preconditioner->ill_conditioned_blocks > 0;
      const std::size_t block_panel_threshold =
          static_cast<std::size_t>(std::ceil(0.70 * static_cast<double>(
                                               std::max<std::size_t>(
                                                   1, panels.size()))));
      if (fmm_nearfield_ilu0 && sparse_ilu_preconditioner != nullptr &&
          sparse_ilu_preconditioner->valid_panel_count == panels.size() &&
          !ilu_ill_conditioned && !block_budget_cap) {
        result.preconditioner_effective = "fmm_nearfield_ilu0";
        result.memory_bytes += result.preconditioner_estimated_bytes;
        custom_preconditioner =
            [sparse_ilu_preconditioner](const std::vector<double> &in,
                                        std::vector<double> &out) {
              sparse_ilu_preconditioner->apply(in, out);
            };
      } else if (!fmm_nearfield_ilu0 && block_preconditioner != nullptr &&
          block_preconditioner->valid_blocks > 0 &&
          block_preconditioner->valid_panel_count >= block_panel_threshold &&
          !block_budget_cap) {
        const bool partial_block =
            block_ill_conditioned ||
            block_preconditioner->valid_panel_count < panels.size();
        if (fmm_overlap_schwarz) {
          result.preconditioner_effective =
              partial_block ? "fmm_overlap_schwarz_partial"
                            : "fmm_overlap_schwarz";
        } else {
          result.preconditioner_effective =
              partial_block ? "fmm_block_jacobi_partial" : "fmm_block_jacobi";
        }
        result.memory_bytes += result.preconditioner_estimated_bytes;
        custom_preconditioner =
            [block_preconditioner, &executor](const std::vector<double> &in,
                                              std::vector<double> &out) {
              block_preconditioner->apply_parallel(in, out, executor);
            };
      } else if (block_budget_cap) {
        block_preconditioner.reset();
        sparse_ilu_preconditioner.reset();
        if (fmm_nearfield_ilu0) {
          result.preconditioner_effective = "diagonal_fmm_ilu_budget_cap";
        } else {
          result.preconditioner_effective =
              fmm_overlap_schwarz ? "diagonal_fmm_overlap_budget_cap"
                                  : "diagonal_fmm_block_budget_cap";
        }
      } else if (block_ill_conditioned || ilu_ill_conditioned) {
        block_preconditioner.reset();
        sparse_ilu_preconditioner.reset();
        if (fmm_nearfield_ilu0) {
          result.preconditioner_effective = "diagonal_fmm_ilu_factor_fallback";
        } else {
          result.preconditioner_effective =
              fmm_overlap_schwarz ? "diagonal_fmm_overlap_ill_conditioned"
                                  : "diagonal_fmm_block_ill_conditioned";
        }
      } else {
        block_preconditioner.reset();
        sparse_ilu_preconditioner.reset();
        if (fmm_nearfield_ilu0) {
          result.preconditioner_effective = "diagonal_fmm_ilu_fallback";
        } else {
          result.preconditioner_effective =
              fmm_overlap_schwarz ? "diagonal_fmm_overlap_fallback"
                                  : "diagonal_fmm_block_fallback";
        }
      }
      progress.write("precond_end", result, &jcp_tree_stats, matvec_calls, 0,
                     0, 0.0, "running", true);
    } else if ((legacy_block || budgeted_block) &&
        !jcp_tree.nodes.empty() &&
        (!legacy_block || panels.size() <= kBlockPanelPrebuildCap)) {
      progress.write("precond_start", result, &jcp_tree_stats, matvec_calls,
                     0, 0, 0.0, "running", true);
      const auto precond_begin = std::chrono::steady_clock::now();
      constexpr std::size_t kBlockMemorySoftCapBytes =
          static_cast<std::size_t>(384) * 1024u * 1024u;
      constexpr std::size_t kBlockWorkingSetSoftCapBytes =
          static_cast<std::size_t>(1024) * 1024u * 1024u;
      const std::size_t block_panel_cap =
          budgeted_block
              ? static_cast<std::size_t>(result.preconditioner_block_panel_cap)
              : static_cast<std::size_t>(32);
      const std::vector<JcpBlockRange> block_ranges =
          make_jcp_block_ranges(jcp_tree, panels.size(), block_panel_cap);
      result.preconditioner_estimated_bytes =
          estimate_jcp_block_preconditioner_bytes(block_ranges, diagonal.size(),
                                                  coupled);
      const bool block_budget_cap =
          budgeted_block &&
          result.preconditioner_estimated_bytes >
              result.preconditioner_memory_budget_bytes;
      const bool block_prebuild_working_set_cap =
          legacy_block &&
          current_working_set_bytes() > kBlockWorkingSetSoftCapBytes;
      if (!block_budget_cap && !block_prebuild_working_set_cap) {
        block_preconditioner = build_jcp_block_preconditioner(
            panels, jcp_tree, coupled, eps_ratio, kappa, lhs_phi, lhs_phin,
            diagonal, result.near_integration, result.quad_buffer,
            result.preconditioner_block_near, block_ranges, executor);
      }
      preconditioner_setup_ms =
          elapsed_ms(precond_begin, std::chrono::steady_clock::now());
      result.preconditioner_setup_ms = preconditioner_setup_ms;
      result.preconditioner_valid_blocks =
          block_preconditioner != nullptr ? block_preconditioner->valid_blocks
                                          : 0;
      result.preconditioner_valid_panels =
          block_preconditioner != nullptr
              ? block_preconditioner->valid_panel_count
              : 0;
      result.preconditioner_ill_conditioned_blocks =
          block_preconditioner != nullptr
              ? block_preconditioner->ill_conditioned_blocks
              : 0;
      result.preconditioner_quality_rejected_blocks =
          block_preconditioner != nullptr
              ? block_preconditioner->quality_rejected_blocks
              : 0;
      result.preconditioner_quality_probe_max_gain =
          block_preconditioner != nullptr
              ? block_preconditioner->quality_probe_max_gain
              : 0.0;
      result.preconditioner_quality_status =
          block_preconditioner != nullptr
              ? block_preconditioner->quality_status
              : std::string("off");
      const bool block_memory_cap =
          legacy_block &&
          (result.preconditioner_estimated_bytes > kBlockMemorySoftCapBytes ||
           block_prebuild_working_set_cap ||
           current_working_set_bytes() > kBlockWorkingSetSoftCapBytes);
      const bool block_ill_conditioned =
          block_preconditioner != nullptr &&
          block_preconditioner->ill_conditioned_blocks > 0;
      const std::size_t block_panel_threshold =
          static_cast<std::size_t>(std::ceil(0.70 * static_cast<double>(
                                               std::max<std::size_t>(
                                                   1, panels.size()))));
      if (block_preconditioner != nullptr &&
           block_preconditioner->valid_blocks > 0 &&
           block_preconditioner->valid_panel_count >= block_panel_threshold &&
           !block_budget_cap && !block_memory_cap) {
        const bool partial_block =
            block_ill_conditioned ||
            block_preconditioner->valid_panel_count < panels.size();
        result.preconditioner_effective =
            block_effective_label(budgeted_block, partial_block);
        result.memory_bytes += result.preconditioner_estimated_bytes;
        custom_preconditioner =
            [block_preconditioner, &executor](const std::vector<double> &in,
                                              std::vector<double> &out) {
              block_preconditioner->apply_parallel(in, out, executor);
            };
      } else if (block_budget_cap) {
        block_preconditioner.reset();
        result.preconditioner_effective = "diagonal_block_budget_cap";
      } else if (block_memory_cap) {
        block_preconditioner.reset();
        result.preconditioner_effective = "diagonal_block_memory_cap";
      } else if (block_ill_conditioned) {
        block_preconditioner.reset();
        result.preconditioner_effective = "diagonal_block_ill_conditioned";
      } else {
        result.preconditioner_effective = "diagonal_block_fallback";
      }
      progress.write("precond_end", result, &jcp_tree_stats, matvec_calls, 0,
                     0, 0.0, "running", true);
    } else if (legacy_block || budgeted_block) {
      progress.write("precond_start", result, &jcp_tree_stats, matvec_calls,
                     0, 0, 0.0, "running", true);
      result.preconditioner_effective =
          legacy_block && panels.size() > kBlockPanelPrebuildCap
              ? "diagonal_block_memory_cap"
              : "diagonal_block_fallback";
      progress.write("precond_end", result, &jcp_tree_stats, matvec_calls, 0,
                     0, 0.0, "running", true);
    } else if (coarse_probe && experimental_treecode &&
               !jcp_tree.nodes.empty()) {
      progress.write("precond_start", result, &jcp_tree_stats, matvec_calls,
                     0, 0, 0.0, "running", true);
      const auto precond_begin = std::chrono::steady_clock::now();
      coarse_preconditioner = make_jcp_coarse_structure(
          panels, jcp_tree, coupled, diagonal,
          static_cast<std::size_t>(result.preconditioner_coarse_aggregate_cap));
      result.preconditioner_coarse_aggregates =
          coarse_preconditioner->aggregates.size();
      result.preconditioner_coarse_dimension = coarse_preconditioner->dimension();
      result.preconditioner_coarse_matrix_bytes =
          coarse_preconditioner->dimension() *
              coarse_preconditioner->dimension() * sizeof(double) +
          coarse_preconditioner->dimension() * sizeof(int);

      std::vector<JcpBlockRange> coarse_base_ranges;
      std::size_t coarse_base_estimated_bytes = 0;
      const bool coarse_block_base =
          result.preconditioner_coarse_base == "block-budgeted";
      if (coarse_block_base) {
        coarse_base_ranges = make_jcp_block_ranges(
            jcp_tree, panels.size(),
            static_cast<std::size_t>(result.preconditioner_block_panel_cap));
        coarse_base_estimated_bytes = estimate_jcp_block_preconditioner_bytes(
            coarse_base_ranges, diagonal.size(), coupled);
      }
      result.preconditioner_estimated_bytes =
          coarse_preconditioner->estimated_bytes() +
          coarse_base_estimated_bytes;
      const bool coarse_budget_cap =
          result.preconditioner_estimated_bytes >
          result.preconditioner_memory_budget_bytes;
      if (coarse_budget_cap) {
        coarse_preconditioner.reset();
        result.preconditioner_coarse_status = "budget_fallback";
        result.preconditioner_effective = "diagonal_coarse_budget_cap";
      } else {
        bool block_base_active = false;
        if (coarse_block_base) {
          block_preconditioner = build_jcp_block_preconditioner(
              panels, jcp_tree, coupled, eps_ratio, kappa, lhs_phi, lhs_phin,
              diagonal, result.near_integration, result.quad_buffer,
              result.preconditioner_block_near, coarse_base_ranges, executor);
          const std::size_t threshold = static_cast<std::size_t>(std::ceil(
              0.70 * static_cast<double>(std::max<std::size_t>(1, panels.size()))));
          block_base_active =
              block_preconditioner != nullptr &&
              block_preconditioner->valid_blocks > 0 &&
              block_preconditioner->valid_panel_count >= threshold;
          if (block_base_active)
            coarse_preconditioner->block_base = block_preconditioner;
        }
        TreeWorkspace setup_tree = jcp_tree;
        setup_tree.moment_audit = "off";
        setup_tree.convergence_audit = "off";
        setup_tree.expansion_audit = "off";
        setup_tree.expansion_audit_channel = "off";
        setup_tree.accepted_miss_diagnostics = "off";
        setup_tree.convergence_overreject_audit = "off";
        setup_tree.convergence_overreject_detail = "off";
        setup_tree.convergence_shape_audit = "off";
        setup_tree.convergence_kernel_response_audit = "off";
        TreecodeStats setup_stats;
        auto setup_matvec = [&](const std::vector<double> &x,
                                std::vector<double> &out) {
          apply_jcp_treecode_operator(
              panels, x, eps_ratio, kappa, lhs_phi, lhs_phin, coupled,
              setup_tree, result.near_integration, result.quad_buffer, executor,
              out, &setup_stats);
        };
        const bool verified =
            options.preconditioner_coarse_verify == "full-small";
        const bool coarse_valid = build_jcp_exact_galerkin(
            panels, *coarse_preconditioner, setup_matvec, verified);
        result.preconditioner_coarse_pivot_ratio =
            coarse_preconditioner->pivot_ratio;
        result.preconditioner_coarse_setup_matvec_calls =
            coarse_preconditioner->setup_matvec_calls;
        result.preconditioner_coarse_rp_identity_error_max =
            coarse_preconditioner->rp_identity_error_max;
        result.preconditioner_coarse_verify_rel_error_max =
            coarse_preconditioner->verify_rel_error_max;
        if (coarse_valid) {
          result.preconditioner_coarse_status = "active";
          result.preconditioner_effective =
              block_base_active
                  ? "coarse_galerkin_probe_block_budgeted"
                  : "coarse_galerkin_probe_diagonal";
          result.memory_bytes += result.preconditioner_estimated_bytes;
          custom_preconditioner =
              [coarse_preconditioner, &panels, &executor,
               &result](const std::vector<double> &in,
                        std::vector<double> &out) {
                const auto begin = std::chrono::steady_clock::now();
                coarse_preconditioner->apply(panels, in, out, executor);
                result.preconditioner_coarse_apply_ms +=
                    elapsed_ms(begin, std::chrono::steady_clock::now());
                ++result.preconditioner_coarse_apply_calls;
              };
        } else {
          coarse_preconditioner.reset();
          result.preconditioner_coarse_status = "factor_fallback";
          result.preconditioner_effective = "diagonal_coarse_factor_fallback";
        }
      }
      preconditioner_setup_ms =
          elapsed_ms(precond_begin, std::chrono::steady_clock::now());
      result.preconditioner_setup_ms = preconditioner_setup_ms;
      result.preconditioner_coarse_setup_ms = preconditioner_setup_ms;
      progress.write("precond_end", result, &jcp_tree_stats, matvec_calls, 0,
                     0, 0.0, "running", true);
    } else if (coarse_probe) {
      progress.write("precond_start", result, &jcp_tree_stats, matvec_calls,
                     0, 0, 0.0, "running", true);
      result.preconditioner_coarse_status = "operator_unsupported";
      result.preconditioner_effective = "diagonal_coarse_unsupported";
      progress.write("precond_end", result, &jcp_tree_stats, matvec_calls, 0,
                     0, 0.0, "running", true);
    }

    const auto gmres_begin = std::chrono::steady_clock::now();
    progress.write("gmres_start", result, &jcp_tree_stats, matvec_calls, 0, 0,
                   0.0, "running", true);
    auto gmres_progress = [&](int gmres_iter, int restart_iter,
                              double true_residual_rel,
                              double preconditioned_residual_rel, bool force,
                              const std::string &phase,
                              const LinearSolveResult &state) {
      result.preconditioner_apply_ms = state.preconditioner_apply_ms;
      result.preconditioner_apply_calls = state.preconditioner_apply_calls;
      result.true_residual_probe_scope = state.true_residual_probe_scope;
      result.true_residual_sample_count = state.true_residual_sample_count;
      result.true_residual_probe_calls = state.true_residual_probe_calls;
      result.true_residual_probe_ms = state.true_residual_probe_ms;
      result.sampled_true_residual_rel = state.sampled_true_residual_rel;
      result.sampled_true_residual_rows = state.sampled_true_residual_rows;
      result.sampled_probe_calls = state.sampled_probe_calls;
      result.sampled_probe_ms = state.sampled_probe_ms;
      result.full_true_residual_confirm_calls =
          state.full_true_residual_confirm_calls;
      result.full_true_residual_confirm_ms =
          state.full_true_residual_confirm_ms;
      if (phase != "gmres_sampled_true_probe" &&
          std::isfinite(true_residual_rel)) {
        result.final_residual = true_residual_rel;
        result.true_residual_rel = true_residual_rel;
        if (result.best_residual <= 0.0 ||
            true_residual_rel < result.best_residual) {
          result.best_residual = true_residual_rel;
          result.best_residual_iter = gmres_iter;
        }
      }
      if (std::isfinite(preconditioned_residual_rel)) {
        result.final_preconditioned_residual =
            preconditioned_residual_rel;
        if (result.best_preconditioned_residual <= 0.0 ||
            preconditioned_residual_rel <
                result.best_preconditioned_residual) {
          result.best_preconditioned_residual =
              preconditioned_residual_rel;
        }
      }
      progress.write(phase, result, &jcp_tree_stats, matvec_calls, gmres_iter,
                     restart_iter, true_residual_rel, "running", force);
    };
    LinearSolveResult linear =
        gmres_solve(matvec, rhs, diagonal, max_iter, tol,
                    options.gmres_restart, executor,
                    custom_preconditioner, gmres_progress,
                    result.true_residual_probe_mode,
                    result.true_residual_probe_interval,
                    result.true_residual_probe_scope,
                    result.true_residual_sample_count,
                    sampled_true_residual);
    result.gmres_ms = elapsed_ms(gmres_begin, std::chrono::steady_clock::now());
    result.iterations = linear.iterations;
    result.residual = linear.residual;
    result.gmres_restart_requested = linear.restart_requested;
    result.gmres_restart_effective = linear.restart_effective;
    result.best_residual = linear.best_residual;
    result.best_residual_iter = linear.best_residual_iter;
    result.final_residual = linear.final_residual;
    result.best_preconditioned_residual =
        linear.best_preconditioned_residual;
    result.final_preconditioned_residual =
        linear.final_preconditioned_residual;
    result.true_residual_abs = linear.true_residual_abs;
    result.true_residual_rel = linear.true_residual_rel;
    result.true_residual_probe_mode = linear.true_residual_probe_mode;
    result.true_residual_probe_scope = linear.true_residual_probe_scope;
    result.true_residual_probe_interval = linear.true_residual_probe_interval;
    result.true_residual_sample_count = linear.true_residual_sample_count;
    result.true_residual_probe_calls = linear.true_residual_probe_calls;
    result.true_residual_probe_ms = linear.true_residual_probe_ms;
    result.sampled_true_residual_rel = linear.sampled_true_residual_rel;
    result.sampled_true_residual_rows = linear.sampled_true_residual_rows;
    result.sampled_probe_calls = linear.sampled_probe_calls;
    result.sampled_probe_ms = linear.sampled_probe_ms;
    result.full_true_residual_confirm_calls =
        linear.full_true_residual_confirm_calls;
    result.full_true_residual_confirm_ms =
        linear.full_true_residual_confirm_ms;
    result.converged = linear.converged;
    result.convergence_reason = linear.convergence_reason;
    result.residual_history = residual_history_string(linear.residual_history);
    result.matvec_ms_total = matvec_ms_total;
    result.matvec_calls = matvec_calls;
    result.matvec_ms_avg =
        matvec_calls > 0 ? matvec_ms_total / static_cast<double>(matvec_calls)
                         : 0.0;
    if (experimental_treecode && options.matvec_audit != "off") {
      const auto audit_begin = std::chrono::steady_clock::now();
      const std::vector<double> &audit_x =
          !linear.best_x.empty() ? linear.best_x : linear.x;
      const std::vector<std::size_t> audit_rows = build_audit_rows(
          audit_x.size(), options.matvec_audit,
          options.matvec_audit_sample_count);
      if (!audit_x.empty() && !audit_rows.empty()) {
        std::vector<double> tree_ax;
        apply_jcp_treecode_operator(panels, audit_x, eps_ratio, kappa, lhs_phi,
                                    lhs_phin, coupled, jcp_tree,
                                    result.near_integration,
                                    result.quad_buffer, executor, tree_ax,
                                    nullptr);
        const bool full =
            options.matvec_audit == "full-small" &&
            audit_rows.size() == audit_x.size();
        double diff_sq = 0.0;
        double direct_sq = 0.0;
        if (full) {
          std::vector<double> direct_ax;
          if (coupled) {
            apply_jcp_linear_pb_operator(
                panels, audit_x, eps_ratio, kappa, lhs_phi, lhs_phin,
                result.near_integration, result.quad_buffer, executor,
                direct_ax, nullptr);
          } else {
            apply_jcp_poisson_operator(
                panels, audit_x, eps_ratio, lhs_phi, result.near_integration,
                result.quad_buffer, executor, direct_ax, nullptr);
          }
          for (std::size_t row : audit_rows) {
            if (row >= direct_ax.size() || row >= tree_ax.size())
              continue;
            const double delta = tree_ax[row] - direct_ax[row];
            diff_sq += delta * delta;
            direct_sq += direct_ax[row] * direct_ax[row];
          }
        } else {
          const std::size_t unknown_count = audit_x.size();
          const std::size_t panel_count = panels.size();
          const std::size_t workers = solver_workers(executor, audit_rows.size());
          std::vector<double> partial_diff_sq(workers, 0.0);
          std::vector<double> partial_direct_sq(workers, 0.0);
          executor.parallel_for(
              audit_rows.size(), workers, 1,
              [&](std::size_t worker, std::size_t begin, std::size_t end) {
                double local_diff_sq = 0.0;
                double local_direct_sq = 0.0;
                for (std::size_t pos = begin; pos < end; ++pos) {
                  const std::size_t row = audit_rows[pos];
                  if (row >= unknown_count || row >= tree_ax.size())
                    continue;
                  const std::size_t target =
                      coupled && row >= panel_count ? row - panel_count : row;
                  if (target >= panel_count)
                    continue;
                  const JcpNearAccum sum = accumulate_jcp_near_ordered(
                      panels[target], target, panels, audit_x, eps_ratio, kappa,
                      coupled, direct_order, 0, panel_count,
                      result.near_integration, result.quad_buffer, nullptr);
                  double direct_row = 0.0;
                  if (coupled) {
                    direct_row = row < panel_count
                                     ? lhs_phi * audit_x[target] - sum.phi
                                     : lhs_phin * audit_x[panel_count + target] -
                                           sum.phin;
                  } else {
                    direct_row = lhs_phi * audit_x[target] - sum.phi;
                  }
                  const double delta = tree_ax[row] - direct_row;
                  local_diff_sq += delta * delta;
                  local_direct_sq += direct_row * direct_row;
                }
                partial_diff_sq[worker] = local_diff_sq;
                partial_direct_sq[worker] = local_direct_sq;
              });
          diff_sq = std::accumulate(partial_diff_sq.begin(),
                                    partial_diff_sq.end(), 0.0);
          direct_sq = std::accumulate(partial_direct_sq.begin(),
                                      partial_direct_sq.end(), 0.0);
        }
        result.tree_vs_direct_matvec_rel =
            std::sqrt(std::max(0.0, diff_sq)) /
            std::max(1.0e-300, std::sqrt(std::max(0.0, direct_sq)));
        result.matvec_audit_rows = audit_rows.size();
      }
      result.matvec_audit_ms =
          elapsed_ms(audit_begin, std::chrono::steady_clock::now());
    }
    if (experimental_treecode) {
      result.matvec_mode = jcp_tree_stats.matvec_mode;
      result.taylor_order = jcp_tree_stats.taylor_order;
      result.derivative_mode = jcp_tree_stats.derivative_mode;
      result.interaction_build_ms = jcp_tree_stats.interaction_build_ms;
      result.tree_eval_ms = jcp_tree_stats.tree_eval_ms;
      result.tree_moment_ms = jcp_tree_stats.tree_moment_ms;
      result.tree_traversal_ms = jcp_tree_stats.tree_traversal_ms;
      result.near_eval_ms = jcp_tree_stats.near_eval_ms;
      result.far_eval_ms = jcp_tree_stats.far_eval_ms;
      result.singular_near_eval_ms = jcp_tree_stats.singular_near_eval_ms;
      result.leaf_near_eval_ms = jcp_tree_stats.leaf_near_eval_ms;
      result.leaf_panel_eval_ms = jcp_tree_stats.leaf_panel_eval_ms;
      result.leaf_tile_eval_ms = jcp_tree_stats.leaf_tile_eval_ms;
      result.leaf_residual_eval_ms = jcp_tree_stats.leaf_residual_eval_ms;
      result.smooth_leaf_node_eval_ms =
          jcp_tree_stats.smooth_leaf_node_eval_ms;
      result.near_interactions = jcp_tree_stats.near_interactions;
      result.far_interactions = jcp_tree_stats.far_interactions;
      result.singular_near_interactions =
          jcp_tree_stats.singular_near_interactions;
      result.leaf_near_interactions = jcp_tree_stats.leaf_near_interactions;
      result.far_node_interactions = jcp_tree_stats.far_node_interactions;
      result.leaf_panel_interactions =
          jcp_tree_stats.leaf_panel_interactions;
      result.leaf_tile_interactions =
          jcp_tree_stats.leaf_tile_interactions;
      result.leaf_tile_panel_interactions =
          jcp_tree_stats.leaf_tile_panel_interactions;
      result.leaf_residual_panel_interactions =
          jcp_tree_stats.leaf_residual_panel_interactions;
      result.smooth_leaf_node_interactions =
          jcp_tree_stats.smooth_leaf_node_interactions;
      result.leaf_pruned_interactions =
          jcp_tree_stats.leaf_pruned_interactions;
      result.leaf_mac_candidates = jcp_tree_stats.leaf_mac_candidates;
      result.leaf_mac_accepted = jcp_tree_stats.leaf_mac_accepted;
      result.leaf_mac_rejected = jcp_tree_stats.leaf_mac_rejected;
      result.leaf_mac_far_panels = jcp_tree_stats.leaf_mac_far_panels;
      result.leaf_mac_direct_panels_saved_estimate =
          jcp_tree_stats.leaf_mac_direct_panels_saved_estimate;
      result.fmm_m2l_mode = jcp_tree_stats.fmm_m2l_mode;
      result.fmm_m2l_near_pairs = jcp_tree_stats.fmm_m2l_near_pairs;
      result.fmm_m2l_far_pairs = jcp_tree_stats.fmm_m2l_far_pairs;
      result.fmm_m2l_setup_ms = jcp_tree_stats.fmm_m2l_setup_ms;
      result.fmm_m2l_apply_ms = jcp_tree_stats.fmm_m2l_apply_ms;
      result.fmm_m2l_finite_outputs = jcp_tree_stats.fmm_m2l_finite_outputs;
      result.leaf_work_profile = jcp_tree_stats.leaf_work_profile;
      result.leaf_work_mac_terminal_panels =
          jcp_tree_stats.leaf_work_mac_terminal_panels;
      result.leaf_work_physical_gap_panels =
          jcp_tree_stats.leaf_work_physical_gap_panels;
      result.leaf_work_fragmentation_panels =
          jcp_tree_stats.leaf_work_fragmentation_panels;
      result.leaf_work_convergence_panels =
          jcp_tree_stats.leaf_work_convergence_panels;
      result.leaf_work_convergence_estimator_panels =
          jcp_tree_stats.leaf_work_convergence_estimator_panels;
      result.leaf_work_convergence_rho_panels =
          jcp_tree_stats.leaf_work_convergence_rho_panels;
      result.leaf_work_convergence_aspect_panels =
          jcp_tree_stats.leaf_work_convergence_aspect_panels;
      result.leaf_work_surface_margin_panels =
          jcp_tree_stats.leaf_work_surface_margin_panels;
      result.leaf_work_surface_margin_hard_panels =
          jcp_tree_stats.leaf_work_surface_margin_hard_panels;
      result.leaf_work_surface_margin_soft_panels =
          jcp_tree_stats.leaf_work_surface_margin_soft_panels;
      result.leaf_work_k3_rho_cap_panels =
          jcp_tree_stats.leaf_work_k3_rho_cap_panels;
      result.leaf_work_singular_fallback_panels =
          jcp_tree_stats.leaf_work_singular_fallback_panels;
      result.leaf_work_scalar_tail_panels =
          jcp_tree_stats.leaf_work_scalar_tail_panels;
      result.leaf_work_profile_accounted_panels =
          jcp_tree_stats.leaf_work_profile_accounted_panels;
      result.leaf_work_mac_terminal_frac =
          jcp_tree_stats.leaf_work_mac_terminal_frac;
      result.leaf_work_convergence_frac =
          jcp_tree_stats.leaf_work_convergence_frac;
      result.leaf_work_surface_margin_frac =
          jcp_tree_stats.leaf_work_surface_margin_frac;
      result.leaf_work_singular_fallback_frac =
          jcp_tree_stats.leaf_work_singular_fallback_frac;
      result.leaf_work_scalar_tail_frac =
          jcp_tree_stats.leaf_work_scalar_tail_frac;
      result.adaptive_smooth_leaf_accepted =
          jcp_tree_stats.adaptive_smooth_leaf_accepted;
      result.adaptive_smooth_leaf_rejected =
          jcp_tree_stats.adaptive_smooth_leaf_rejected;
      result.adaptive_span_tile_accepted =
          jcp_tree_stats.adaptive_span_tile_accepted;
      result.adaptive_span_tile_rejected =
          jcp_tree_stats.adaptive_span_tile_rejected;
      result.adaptive_exact_leaf_fallbacks =
          jcp_tree_stats.adaptive_exact_leaf_fallbacks;
      result.adaptive_reject_one_ring_overlap =
          jcp_tree_stats.adaptive_reject_one_ring_overlap;
      result.adaptive_reject_error_bound_failed =
          jcp_tree_stats.adaptive_reject_error_bound_failed;
      result.adaptive_reject_p4_check_failed =
          jcp_tree_stats.adaptive_reject_p4_check_failed;
      result.adaptive_reject_too_close =
          jcp_tree_stats.adaptive_reject_too_close;
      result.leaf_work_reduction_ratio =
          jcp_tree_stats.leaf_work_reduction_ratio;
      result.target_cluster_interactions =
          jcp_tree_stats.target_cluster_interactions;
      result.target_cluster_targets = jcp_tree_stats.target_cluster_targets;
      result.target_cluster_panels_pruned =
          jcp_tree_stats.target_cluster_panels_pruned;
      result.target_cluster_rejected =
          jcp_tree_stats.target_cluster_rejected;
      result.target_cluster_eval_ms = jcp_tree_stats.target_cluster_eval_ms;
      result.leaf_pair_batch_interactions =
          jcp_tree_stats.leaf_pair_batch_interactions;
      result.leaf_pair_batch_targets =
          jcp_tree_stats.leaf_pair_batch_targets;
      result.leaf_pair_batch_source_panels =
          jcp_tree_stats.leaf_pair_batch_source_panels;
      result.leaf_pair_batch_residual_panels =
          jcp_tree_stats.leaf_pair_batch_residual_panels;
      result.leaf_pair_batch_fallbacks =
          jcp_tree_stats.leaf_pair_batch_fallbacks;
      result.leaf_pair_batch_ms = jcp_tree_stats.leaf_pair_batch_ms;
      result.leaf_pair_batch_kernel = jcp_tree_stats.leaf_pair_batch_kernel;
      result.leaf_pair_batch_vectorized_sources =
          jcp_tree_stats.leaf_pair_batch_vectorized_sources;
      result.leaf_pair_batch_vector_lanes =
          jcp_tree_stats.leaf_pair_batch_vector_lanes;
      result.leaf_pair_batch_pairs_per_second =
          jcp_tree_stats.leaf_pair_batch_pairs_per_second;
      result.leaf_kernel_profile = jcp_tree_stats.leaf_kernel_profile;
      result.leaf_work_profile = jcp_tree_stats.leaf_work_profile;
      result.jcp_direct_kernel = jcp_tree_stats.jcp_direct_kernel;
      result.leaf_kernel_scalar_tasks =
          jcp_tree_stats.leaf_kernel_scalar_tasks;
      result.leaf_kernel_batch_tasks =
          jcp_tree_stats.leaf_kernel_batch_tasks;
      result.leaf_kernel_source_panel_evals =
          jcp_tree_stats.leaf_kernel_source_panel_evals;
      result.leaf_kernel_scalar_ms = jcp_tree_stats.leaf_kernel_scalar_ms;
      result.leaf_kernel_target_avx2_ms =
          jcp_tree_stats.leaf_kernel_target_avx2_ms;
      result.leaf_kernel_tail_ms = jcp_tree_stats.leaf_kernel_tail_ms;
      result.leaf_pair_batch_soa_sources =
          jcp_tree_stats.leaf_pair_batch_soa_sources;
      result.leaf_pair_batch_soa_ms =
          jcp_tree_stats.leaf_pair_batch_soa_ms;
      result.leaf_pair_batch_scalar_tail_ms =
          jcp_tree_stats.leaf_pair_batch_scalar_tail_ms;
      result.leaf_pair_batch_clean_source_panels_per_second =
          jcp_tree_stats.leaf_pair_batch_clean_source_panels_per_second;
      result.leaf_pair_selection_profile =
          jcp_tree_stats.leaf_pair_selection_profile;
      result.leaf_pair_probe_cache = jcp_tree_stats.leaf_pair_probe_cache;
      result.leaf_pair_target_group = jcp_tree_stats.leaf_pair_target_group;
      result.leaf_pair_directness_probe =
          jcp_tree_stats.leaf_pair_directness_probe;
      result.leaf_pair_directness_cache =
          jcp_tree_stats.leaf_pair_directness_cache;
      result.leaf_pair_directness_verify =
          jcp_tree_stats.leaf_pair_directness_verify;
      result.leaf_pair_max_lookahead =
          jcp_tree_stats.leaf_pair_max_lookahead;
      result.leaf_pair_selection_candidate_scans =
          jcp_tree_stats.leaf_pair_selection_candidate_scans;
      result.leaf_pair_selection_clean_rejects =
          jcp_tree_stats.leaf_pair_selection_clean_rejects;
      result.leaf_pair_selection_probe_calls =
          jcp_tree_stats.leaf_pair_selection_probe_calls;
      result.leaf_pair_selection_probe_node_visits =
          jcp_tree_stats.leaf_pair_selection_probe_node_visits;
      result.leaf_pair_selection_probe_ms =
          jcp_tree_stats.leaf_pair_selection_probe_ms;
      result.leaf_pair_selection_ms =
          jcp_tree_stats.leaf_pair_selection_ms;
      result.leaf_pair_selection_probe_cache_hits =
          jcp_tree_stats.leaf_pair_selection_probe_cache_hits;
      result.leaf_pair_selection_probe_cache_misses =
          jcp_tree_stats.leaf_pair_selection_probe_cache_misses;
      result.leaf_pair_selection_insufficient_targets =
          jcp_tree_stats.leaf_pair_selection_insufficient_targets;
      result.leaf_pair_selection_same_target_leaf_rejects =
          jcp_tree_stats.leaf_pair_selection_same_target_leaf_rejects;
      result.leaf_pair_directness_ancestor_probe_calls =
          jcp_tree_stats.leaf_pair_directness_ancestor_probe_calls;
      result.leaf_pair_directness_ancestor_probe_node_visits =
          jcp_tree_stats.leaf_pair_directness_ancestor_probe_node_visits;
      result.leaf_pair_directness_ancestor_probe_ms =
          jcp_tree_stats.leaf_pair_directness_ancestor_probe_ms;
      result.leaf_pair_directness_pair_cache_hits =
          jcp_tree_stats.leaf_pair_directness_pair_cache_hits;
      result.leaf_pair_directness_pair_cache_misses =
          jcp_tree_stats.leaf_pair_directness_pair_cache_misses;
      result.leaf_pair_directness_all_direct_proofs =
          jcp_tree_stats.leaf_pair_directness_all_direct_proofs;
      result.leaf_pair_directness_mixed_fallbacks =
          jcp_tree_stats.leaf_pair_directness_mixed_fallbacks;
      result.leaf_pair_directness_verify_checks =
          jcp_tree_stats.leaf_pair_directness_verify_checks;
      result.leaf_pair_directness_verify_mismatches =
          jcp_tree_stats.leaf_pair_directness_verify_mismatches;
      result.leaf_pair_scheduler = jcp_tree_stats.leaf_pair_scheduler;
      result.leaf_pair_worklist_max_active =
          jcp_tree_stats.leaf_pair_worklist_max_active;
      result.leaf_pair_worklist_verify =
          jcp_tree_stats.leaf_pair_worklist_verify;
      result.leaf_pair_worklist_tasks =
          jcp_tree_stats.leaf_pair_worklist_tasks;
      result.leaf_pair_worklist_groups =
          jcp_tree_stats.leaf_pair_worklist_groups;
      result.leaf_pair_worklist_batched_groups =
          jcp_tree_stats.leaf_pair_worklist_batched_groups;
      result.leaf_pair_worklist_scalar_flush_tasks =
          jcp_tree_stats.leaf_pair_worklist_scalar_flush_tasks;
      result.leaf_pair_worklist_forced_flushes =
          jcp_tree_stats.leaf_pair_worklist_forced_flushes;
      result.leaf_pair_worklist_active_groups_max =
          jcp_tree_stats.leaf_pair_worklist_active_groups_max;
      result.leaf_pair_worklist_ms = jcp_tree_stats.leaf_pair_worklist_ms;
      result.leaf_pair_worklist_verify_checks =
          jcp_tree_stats.leaf_pair_worklist_verify_checks;
      result.leaf_pair_worklist_verify_mismatches =
          jcp_tree_stats.leaf_pair_worklist_verify_mismatches;
      result.leaf_block_verify = jcp_tree_stats.leaf_block_verify;
      result.leaf_block_target_leaves =
          jcp_tree_stats.leaf_block_target_leaves;
      result.leaf_block_targets = jcp_tree_stats.leaf_block_targets;
      result.leaf_block_source_groups =
          jcp_tree_stats.leaf_block_source_groups;
      result.leaf_block_batched_groups =
          jcp_tree_stats.leaf_block_batched_groups;
      result.leaf_block_scalar_tail_tasks =
          jcp_tree_stats.leaf_block_scalar_tail_tasks;
      result.leaf_block_group_ms = jcp_tree_stats.leaf_block_group_ms;
      result.leaf_block_flush_ms = jcp_tree_stats.leaf_block_flush_ms;
      result.leaf_block_writeback_ms =
          jcp_tree_stats.leaf_block_writeback_ms;
      result.leaf_block_verify_checks =
          jcp_tree_stats.leaf_block_verify_checks;
      result.leaf_block_verify_mismatches =
          jcp_tree_stats.leaf_block_verify_mismatches;
      result.adaptive_level = jcp_tree_stats.adaptive_level;
      result.cluster_policy = jcp_tree_stats.cluster_policy;
      result.target_batch_size = jcp_tree_stats.target_batch_size;
      result.leaf_pair_min_source_panels =
          jcp_tree_stats.leaf_pair_min_source_panels;
      result.leaf_pair_min_targets = jcp_tree_stats.leaf_pair_min_targets;
      result.leaf_near_policy = jcp_tree_stats.leaf_near_policy;
      result.leaf_panel_kernel = jcp_tree_stats.leaf_panel_kernel;
      result.interaction_mode = jcp_tree_stats.interaction_mode;
      result.node_bounds_mode = jcp_tree_stats.node_bounds_mode;
      result.tree_builder = jcp_tree_stats.tree_builder;
      result.expansion_center = jcp_tree_stats.expansion_center;
      result.fragmentation_safety = jcp_tree_stats.fragmentation_safety;
      result.mac_safety = jcp_tree_stats.mac_safety;
      result.mac_gap_factor = jcp_tree_stats.mac_gap_factor;
      result.tree_depth_max = jcp_tree_stats.tree_depth_max;
      result.mac_reject_radius_limited =
          jcp_tree_stats.mac_reject_radius_limited;
      result.mac_accept_rate = jcp_tree_stats.mac_accept_rate;
      result.node_radius_max = jcp_tree_stats.node_radius_max;
      result.node_radius_mean = jcp_tree_stats.node_radius_mean;
      result.node_aspect_ratio_max = jcp_tree_stats.node_aspect_ratio_max;
      result.node_aspect_ratio_mean = jcp_tree_stats.node_aspect_ratio_mean;
      result.leaf_panel_work_per_matvec =
          jcp_tree_stats.leaf_panel_work_per_matvec;
      result.far_node_work_per_matvec =
          jcp_tree_stats.far_node_work_per_matvec;
      result.far_geometric_unsafe_count =
          jcp_tree_stats.far_geometric_unsafe_count;
      result.far_geometric_unsafe_panels =
          jcp_tree_stats.far_geometric_unsafe_panels;
      result.far_geometric_min_gap = jcp_tree_stats.far_geometric_min_gap;
      result.far_geometric_gap_ratio_min =
          jcp_tree_stats.far_geometric_gap_ratio_min;
      result.far_geometric_gap_ratio_mean =
          jcp_tree_stats.far_geometric_gap_ratio_mean;
      result.far_geometric_audit_samples =
          jcp_tree_stats.far_geometric_audit_samples;
      result.mac_physical_gap_rejects =
          jcp_tree_stats.mac_physical_gap_rejects;
      result.mac_convergence_rho_rejects =
          jcp_tree_stats.mac_convergence_rho_rejects;
      result.mac_convergence_aspect_rejects =
          jcp_tree_stats.mac_convergence_aspect_rejects;
      result.mac_convergence_rho_reject_panels =
          jcp_tree_stats.mac_convergence_rho_reject_panels;
      result.mac_convergence_aspect_reject_panels =
          jcp_tree_stats.mac_convergence_aspect_reject_panels;
      result.fragmented_far_count = jcp_tree_stats.fragmented_far_count;
      result.fragmented_far_panels = jcp_tree_stats.fragmented_far_panels;
      result.fragmented_split_rejects =
          jcp_tree_stats.fragmented_split_rejects;
      result.fragmentation_center_gap_ratio_max =
          jcp_tree_stats.fragmentation_center_gap_ratio_max;
      result.fragmentation_center_gap_ratio_mean =
          jcp_tree_stats.fragmentation_center_gap_ratio_mean;
      result.moment_audit = jcp_tree_stats.moment_audit;
      result.moment_construction = jcp_tree_stats.moment_construction;
      result.moment_sum = jcp_tree_stats.moment_sum;
      result.convergence_audit = jcp_tree_stats.convergence_audit;
      result.convergence_safety = jcp_tree_stats.convergence_safety;
      result.convergence_rho_cap = jcp_tree_stats.convergence_rho_cap;
      result.mac_aspect_cap = jcp_tree_stats.mac_aspect_cap;
      result.convergence_estimator = jcp_tree_stats.convergence_estimator;
      result.convergence_tail_threshold =
          jcp_tree_stats.convergence_tail_threshold;
      result.convergence_estimator_error_cap =
          jcp_tree_stats.convergence_estimator_error_cap;
      result.convergence_unstable_rho_cap =
          jcp_tree_stats.convergence_unstable_rho_cap;
      result.convergence_calibration =
          jcp_tree_stats.convergence_calibration;
      result.convergence_dielectric_amp =
          jcp_tree_stats.convergence_dielectric_amp;
      result.convergence_estimator_model =
          jcp_tree_stats.convergence_estimator_model;
      result.convergence_channel_calibration =
          jcp_tree_stats.convergence_channel_calibration;
      result.convergence_split_policy =
          jcp_tree_stats.convergence_split_policy;
      result.convergence_split_selectivity =
          jcp_tree_stats.convergence_split_selectivity;
      result.convergence_split_absolute_floor =
          jcp_tree_stats.convergence_split_absolute_floor;
      result.convergence_overreject_audit =
          jcp_tree_stats.convergence_overreject_audit;
      result.convergence_overreject_detail =
          jcp_tree_stats.convergence_overreject_detail;
      result.convergence_shape_audit =
          jcp_tree_stats.convergence_shape_audit;
      result.convergence_kernel_response_audit =
          jcp_tree_stats.convergence_kernel_response_audit;
      result.channel_missed_semantics =
          jcp_tree_stats.channel_missed_semantics;
      result.accepted_miss_diagnostics =
          jcp_tree_stats.accepted_miss_diagnostics;
      result.accepted_miss_absolute_floor =
          jcp_tree_stats.accepted_miss_absolute_floor;
      result.convergence_k3_accepted_rho_cap =
          jcp_tree_stats.convergence_k3_accepted_rho_cap;
      result.mac_distance_policy = jcp_tree_stats.mac_distance_policy;
      result.mac_surface_margin_factor =
          jcp_tree_stats.mac_surface_margin_factor;
      result.mac_surface_margin_policy =
          jcp_tree_stats.mac_surface_margin_policy;
      result.mac_surface_margin_aspect_threshold =
          jcp_tree_stats.mac_surface_margin_aspect_threshold;
      result.mac_surface_margin_rho_trigger =
          jcp_tree_stats.mac_surface_margin_rho_trigger;
      result.mac_surface_margin_soft_aspect_threshold =
          jcp_tree_stats.mac_surface_margin_soft_aspect_threshold;
      result.mac_surface_margin_soft_rho_trigger =
          jcp_tree_stats.mac_surface_margin_soft_rho_trigger;
      result.mac_surface_margin_soft_factor =
          jcp_tree_stats.mac_surface_margin_soft_factor;
      result.convergence_estimator_samples =
          jcp_tree_stats.convergence_estimator_samples;
      result.convergence_estimator_unstable_far_count =
          jcp_tree_stats.convergence_estimator_unstable_far_count;
      result.convergence_estimator_unstable_far_panels =
          jcp_tree_stats.convergence_estimator_unstable_far_panels;
      result.convergence_estimator_split_rejects =
          jcp_tree_stats.convergence_estimator_split_rejects;
      result.convergence_estimator_split_reject_panels =
          jcp_tree_stats.convergence_estimator_split_reject_panels;
      result.convergence_low_abs_reaccepted_far_count =
          jcp_tree_stats.convergence_low_abs_reaccepted_far_count;
      result.convergence_low_abs_reaccepted_far_panels =
          jcp_tree_stats.convergence_low_abs_reaccepted_far_panels;
      result.convergence_low_abs_reaccepted_far_count =
          jcp_tree_stats.convergence_low_abs_reaccepted_far_count;
      result.convergence_low_abs_reaccepted_far_panels =
          jcp_tree_stats.convergence_low_abs_reaccepted_far_panels;
      result.convergence_estimator_tail_ratio_max =
          jcp_tree_stats.convergence_estimator_tail_ratio_max;
      result.convergence_estimator_tail_ratio_mean =
          jcp_tree_stats.convergence_estimator_tail_ratio_mean;
      result.convergence_estimator_error_max =
          jcp_tree_stats.convergence_estimator_error_max;
      result.convergence_estimator_error_mean =
          jcp_tree_stats.convergence_estimator_error_mean;
      result.convergence_estimator_abs_error_max =
          jcp_tree_stats.convergence_estimator_abs_error_max;
      result.convergence_estimator_abs_error_mean =
          jcp_tree_stats.convergence_estimator_abs_error_mean;
      result.convergence_estimator_abs_error_max =
          jcp_tree_stats.convergence_estimator_abs_error_max;
      result.convergence_estimator_abs_error_mean =
          jcp_tree_stats.convergence_estimator_abs_error_mean;
      result.convergence_tail_norm_p_max =
          jcp_tree_stats.convergence_tail_norm_p_max;
      result.convergence_tail_norm_p_mean =
          jcp_tree_stats.convergence_tail_norm_p_mean;
      result.convergence_tail_norm_pm1_max =
          jcp_tree_stats.convergence_tail_norm_pm1_max;
      result.convergence_tail_norm_pm1_mean =
          jcp_tree_stats.convergence_tail_norm_pm1_mean;
      result.convergence_effective_jcp_amp_max =
          jcp_tree_stats.convergence_effective_jcp_amp_max;
      result.convergence_effective_jcp_amp_mean =
          jcp_tree_stats.convergence_effective_jcp_amp_mean;
      result.convergence_dielectric_amp_max =
          jcp_tree_stats.convergence_dielectric_amp_max;
      result.convergence_dielectric_amp_mean =
          jcp_tree_stats.convergence_dielectric_amp_mean;
      result.convergence_channel_estimated_error_max =
          jcp_tree_stats.convergence_channel_estimated_error_max;
      result.convergence_channel_estimated_error_mean =
          jcp_tree_stats.convergence_channel_estimated_error_mean;
      result.convergence_channel_effective_jcp_amp_max =
          jcp_tree_stats.convergence_channel_effective_jcp_amp_max;
      result.convergence_channel_effective_jcp_amp_mean =
          jcp_tree_stats.convergence_channel_effective_jcp_amp_mean;
      result.convergence_channel_unstable_far_count =
          jcp_tree_stats.convergence_channel_unstable_far_count;
      result.convergence_channel_unstable_far_panels =
          jcp_tree_stats.convergence_channel_unstable_far_panels;
      result.convergence_channel_split_rejects =
          jcp_tree_stats.convergence_channel_split_rejects;
      result.convergence_channel_split_reject_panels =
          jcp_tree_stats.convergence_channel_split_reject_panels;
      result.convergence_channel_low_abs_reaccepted_count =
          jcp_tree_stats.convergence_channel_low_abs_reaccepted_count;
      result.convergence_channel_low_abs_reaccepted_panels =
          jcp_tree_stats.convergence_channel_low_abs_reaccepted_panels;
      result.convergence_channel_high_error_caught =
          jcp_tree_stats.convergence_channel_high_error_caught;
      result.convergence_channel_high_error_missed =
          jcp_tree_stats.convergence_channel_high_error_missed;
      result.convergence_channel_high_error_caught_by_split =
          jcp_tree_stats.convergence_channel_high_error_caught_by_split;
      result.convergence_estimator_high_error_caught =
          jcp_tree_stats.convergence_estimator_high_error_caught;
      result.convergence_estimator_high_error_missed =
          jcp_tree_stats.convergence_estimator_high_error_missed;
      result.convergence_estimator_high_error_caught_by_split =
          jcp_tree_stats.convergence_estimator_high_error_caught_by_split;
      result.convergence_estimator_phantom_high_error_count =
          jcp_tree_stats.convergence_estimator_phantom_high_error_count;
      result.convergence_estimator_phantom_high_error_abs_max =
          jcp_tree_stats.convergence_estimator_phantom_high_error_abs_max;
      result.convergence_overreject_count =
          jcp_tree_stats.convergence_overreject_count;
      result.convergence_overreject_panels =
          jcp_tree_stats.convergence_overreject_panels;
      result.convergence_overreject_caught_count =
          jcp_tree_stats.convergence_overreject_caught_count;
      result.convergence_overreject_caught_panels =
          jcp_tree_stats.convergence_overreject_caught_panels;
      result.convergence_overreject_abs_error_max =
          jcp_tree_stats.convergence_overreject_abs_error_max;
      result.convergence_overreject_abs_error_mean =
          jcp_tree_stats.convergence_overreject_abs_error_mean;
      result.convergence_overreject_abs_error_p95 =
          jcp_tree_stats.convergence_overreject_abs_error_p95;
      result.convergence_overreject_abs_error_p99 =
          jcp_tree_stats.convergence_overreject_abs_error_p99;
      result.convergence_overreject_estimated_abs_error_max =
          jcp_tree_stats.convergence_overreject_estimated_abs_error_max;
      result.convergence_overreject_estimated_abs_error_mean =
          jcp_tree_stats.convergence_overreject_estimated_abs_error_mean;
      result.convergence_overreject_required_scale_max =
          jcp_tree_stats.convergence_overreject_required_scale_max;
      result.convergence_overreject_required_scale_p95 =
          jcp_tree_stats.convergence_overreject_required_scale_p95;
      result.convergence_overreject_required_scale_p99 =
          jcp_tree_stats.convergence_overreject_required_scale_p99;
      result.convergence_channel_overreject_count =
          jcp_tree_stats.convergence_channel_overreject_count;
      result.convergence_channel_overreject_caught_count =
          jcp_tree_stats.convergence_channel_overreject_caught_count;
      result.convergence_channel_overreject_panels =
          jcp_tree_stats.convergence_channel_overreject_panels;
      result.convergence_channel_overreject_abs_error_max =
          jcp_tree_stats.convergence_channel_overreject_abs_error_max;
      result.convergence_channel_overreject_abs_error_mean =
          jcp_tree_stats.convergence_channel_overreject_abs_error_mean;
      result.convergence_channel_overreject_required_scale_max =
          jcp_tree_stats.convergence_channel_overreject_required_scale_max;
      result.convergence_channel_overreject_required_scale_p99 =
          jcp_tree_stats.convergence_channel_overreject_required_scale_p99;
      result.phantom_high_error_channel_count =
          jcp_tree_stats.phantom_high_error_channel_count;
      result.phantom_high_error_channel_abs_max =
          jcp_tree_stats.phantom_high_error_channel_abs_max;
      result.accepted_miss_channel_count =
          jcp_tree_stats.accepted_miss_channel_count;
      result.accepted_miss_channel_panels =
          jcp_tree_stats.accepted_miss_channel_panels;
      result.accepted_miss_channel_rho_min =
          jcp_tree_stats.accepted_miss_channel_rho_min;
      result.accepted_miss_channel_rho_mean =
          jcp_tree_stats.accepted_miss_channel_rho_mean;
      result.accepted_miss_channel_rho_max =
          jcp_tree_stats.accepted_miss_channel_rho_max;
      result.accepted_miss_channel_tail_ratio_mean =
          jcp_tree_stats.accepted_miss_channel_tail_ratio_mean;
      result.accepted_miss_channel_tail_ratio_max =
          jcp_tree_stats.accepted_miss_channel_tail_ratio_max;
      result.accepted_miss_channel_estimated_error_mean =
          jcp_tree_stats.accepted_miss_channel_estimated_error_mean;
      result.accepted_miss_channel_estimated_error_max =
          jcp_tree_stats.accepted_miss_channel_estimated_error_max;
      result.accepted_miss_worst_channel =
          jcp_tree_stats.accepted_miss_worst_channel;
      result.convergence_k3_rho_cap_rejects =
          jcp_tree_stats.convergence_k3_rho_cap_rejects;
      result.convergence_k3_rho_cap_reject_panels =
          jcp_tree_stats.convergence_k3_rho_cap_reject_panels;
      result.mac_surface_margin_rejects =
          jcp_tree_stats.mac_surface_margin_rejects;
      result.mac_surface_margin_reject_panels =
          jcp_tree_stats.mac_surface_margin_reject_panels;
      result.mac_surface_margin_applied =
          jcp_tree_stats.mac_surface_margin_applied;
      result.mac_surface_margin_skipped =
          jcp_tree_stats.mac_surface_margin_skipped;
      result.mac_surface_margin_adaptive_rejects =
          jcp_tree_stats.mac_surface_margin_adaptive_rejects;
      result.mac_surface_margin_adaptive_reject_panels =
          jcp_tree_stats.mac_surface_margin_adaptive_reject_panels;
      result.mac_surface_margin_skipped_accepted_far =
          jcp_tree_stats.mac_surface_margin_skipped_accepted_far;
      result.mac_surface_margin_hard_applied =
          jcp_tree_stats.mac_surface_margin_hard_applied;
      result.mac_surface_margin_soft_applied =
          jcp_tree_stats.mac_surface_margin_soft_applied;
      result.mac_surface_margin_hard_rejects =
          jcp_tree_stats.mac_surface_margin_hard_rejects;
      result.mac_surface_margin_hard_reject_panels =
          jcp_tree_stats.mac_surface_margin_hard_reject_panels;
      result.mac_surface_margin_soft_rejects =
          jcp_tree_stats.mac_surface_margin_soft_rejects;
      result.mac_surface_margin_soft_reject_panels =
          jcp_tree_stats.mac_surface_margin_soft_reject_panels;
      result.moment_audit_nodes = jcp_tree_stats.moment_audit_nodes;
      result.moment_translation_rel_error_max =
          jcp_tree_stats.moment_translation_rel_error_max;
      result.moment_translation_rel_error_mean =
          jcp_tree_stats.moment_translation_rel_error_mean;
      result.moment_translation_offset_max =
          jcp_tree_stats.moment_translation_offset_max;
      result.moment_translation_offset_mean =
          jcp_tree_stats.moment_translation_offset_mean;
      result.moment_cancellation_ratio_max =
          jcp_tree_stats.moment_cancellation_ratio_max;
      result.moment_cancellation_ratio_mean =
          jcp_tree_stats.moment_cancellation_ratio_mean;
      result.moment_dynamic_range_max =
          jcp_tree_stats.moment_dynamic_range_max;
      result.expansion_audit = jcp_tree_stats.expansion_audit;
      result.expansion_audit_channel =
          jcp_tree_stats.expansion_audit_channel;
      result.far_expansion_audit_pairs =
          jcp_tree_stats.far_expansion_audit_pairs;
      result.far_expansion_audit_panels =
          jcp_tree_stats.far_expansion_audit_panels;
      result.far_expansion_direct_ms =
          jcp_tree_stats.far_expansion_direct_ms;
      result.far_expansion_rel_error_max =
          jcp_tree_stats.far_expansion_rel_error_max;
      result.far_expansion_rel_error_mean =
          jcp_tree_stats.far_expansion_rel_error_mean;
      result.far_expansion_rel_error_p95 =
          jcp_tree_stats.far_expansion_rel_error_p95;
      result.far_expansion_rel_error_p99 =
          jcp_tree_stats.far_expansion_rel_error_p99;
      result.far_expansion_channel_rel_error_max =
          jcp_tree_stats.far_expansion_channel_rel_error_max;
      result.far_expansion_channel_rel_error_mean =
          jcp_tree_stats.far_expansion_channel_rel_error_mean;
      result.far_expansion_channel_rel_error_p99 =
          jcp_tree_stats.far_expansion_channel_rel_error_p99;
      result.far_expansion_worst_channel =
          jcp_tree_stats.far_expansion_worst_channel;
      result.far_expansion_mean_gate_pass =
          jcp_tree_stats.far_expansion_mean_gate_pass;
      result.far_convergence_rho_max =
          jcp_tree_stats.far_convergence_rho_max;
      result.far_convergence_rho_mean =
          jcp_tree_stats.far_convergence_rho_mean;
      result.accepted_far_aspect_ratio_max =
          jcp_tree_stats.accepted_far_aspect_ratio_max;
      result.accepted_far_aspect_ratio_mean =
          jcp_tree_stats.accepted_far_aspect_ratio_mean;
      result.leaf_tile_status = jcp_tree_stats.leaf_tile_status;
      result.leaf_panel_pairs_per_second =
          jcp_tree_stats.leaf_panel_pairs_per_second;
    } else {
      result.near_eval_ms = jcp_tree_stats.near_eval_ms;
      result.near_interactions = jcp_tree_stats.near_interactions;
      result.leaf_near_policy = jcp_tree_stats.leaf_near_policy;
      result.leaf_panel_kernel = jcp_tree_stats.leaf_panel_kernel;
      result.interaction_mode = jcp_tree_stats.interaction_mode;
      result.node_bounds_mode = jcp_tree_stats.node_bounds_mode;
      result.tree_builder = jcp_tree_stats.tree_builder;
      result.expansion_center = jcp_tree_stats.expansion_center;
      result.fragmentation_safety = jcp_tree_stats.fragmentation_safety;
      result.mac_safety = jcp_tree_stats.mac_safety;
      result.mac_gap_factor = jcp_tree_stats.mac_gap_factor;
      result.tree_depth_max = jcp_tree_stats.tree_depth_max;
      result.mac_reject_radius_limited =
          jcp_tree_stats.mac_reject_radius_limited;
      result.mac_accept_rate = jcp_tree_stats.mac_accept_rate;
      result.node_radius_max = jcp_tree_stats.node_radius_max;
      result.node_radius_mean = jcp_tree_stats.node_radius_mean;
      result.node_aspect_ratio_max = jcp_tree_stats.node_aspect_ratio_max;
      result.node_aspect_ratio_mean = jcp_tree_stats.node_aspect_ratio_mean;
      result.leaf_panel_work_per_matvec =
          jcp_tree_stats.leaf_panel_work_per_matvec;
      result.far_node_work_per_matvec =
          jcp_tree_stats.far_node_work_per_matvec;
      result.far_geometric_unsafe_count =
          jcp_tree_stats.far_geometric_unsafe_count;
      result.far_geometric_unsafe_panels =
          jcp_tree_stats.far_geometric_unsafe_panels;
      result.far_geometric_min_gap = jcp_tree_stats.far_geometric_min_gap;
      result.far_geometric_gap_ratio_min =
          jcp_tree_stats.far_geometric_gap_ratio_min;
      result.far_geometric_gap_ratio_mean =
          jcp_tree_stats.far_geometric_gap_ratio_mean;
      result.far_geometric_audit_samples =
          jcp_tree_stats.far_geometric_audit_samples;
      result.mac_physical_gap_rejects =
          jcp_tree_stats.mac_physical_gap_rejects;
      result.mac_convergence_rho_rejects =
          jcp_tree_stats.mac_convergence_rho_rejects;
      result.mac_convergence_aspect_rejects =
          jcp_tree_stats.mac_convergence_aspect_rejects;
      result.mac_convergence_rho_reject_panels =
          jcp_tree_stats.mac_convergence_rho_reject_panels;
      result.mac_convergence_aspect_reject_panels =
          jcp_tree_stats.mac_convergence_aspect_reject_panels;
      result.fragmented_far_count = jcp_tree_stats.fragmented_far_count;
      result.fragmented_far_panels = jcp_tree_stats.fragmented_far_panels;
      result.fragmented_split_rejects =
          jcp_tree_stats.fragmented_split_rejects;
      result.fragmentation_center_gap_ratio_max =
          jcp_tree_stats.fragmentation_center_gap_ratio_max;
      result.fragmentation_center_gap_ratio_mean =
          jcp_tree_stats.fragmentation_center_gap_ratio_mean;
      result.moment_audit = jcp_tree_stats.moment_audit;
      result.moment_construction = jcp_tree_stats.moment_construction;
      result.moment_sum = jcp_tree_stats.moment_sum;
      result.convergence_audit = jcp_tree_stats.convergence_audit;
      result.convergence_safety = jcp_tree_stats.convergence_safety;
      result.convergence_rho_cap = jcp_tree_stats.convergence_rho_cap;
      result.mac_aspect_cap = jcp_tree_stats.mac_aspect_cap;
      result.convergence_estimator = jcp_tree_stats.convergence_estimator;
      result.convergence_tail_threshold =
          jcp_tree_stats.convergence_tail_threshold;
      result.convergence_estimator_error_cap =
          jcp_tree_stats.convergence_estimator_error_cap;
      result.convergence_unstable_rho_cap =
          jcp_tree_stats.convergence_unstable_rho_cap;
      result.convergence_calibration =
          jcp_tree_stats.convergence_calibration;
      result.convergence_estimator_samples =
          jcp_tree_stats.convergence_estimator_samples;
      result.convergence_estimator_unstable_far_count =
          jcp_tree_stats.convergence_estimator_unstable_far_count;
      result.convergence_estimator_unstable_far_panels =
          jcp_tree_stats.convergence_estimator_unstable_far_panels;
      result.convergence_estimator_split_rejects =
          jcp_tree_stats.convergence_estimator_split_rejects;
      result.convergence_estimator_split_reject_panels =
          jcp_tree_stats.convergence_estimator_split_reject_panels;
      result.convergence_estimator_tail_ratio_max =
          jcp_tree_stats.convergence_estimator_tail_ratio_max;
      result.convergence_estimator_tail_ratio_mean =
          jcp_tree_stats.convergence_estimator_tail_ratio_mean;
      result.convergence_estimator_error_max =
          jcp_tree_stats.convergence_estimator_error_max;
      result.convergence_estimator_error_mean =
          jcp_tree_stats.convergence_estimator_error_mean;
      result.convergence_tail_norm_p_max =
          jcp_tree_stats.convergence_tail_norm_p_max;
      result.convergence_tail_norm_p_mean =
          jcp_tree_stats.convergence_tail_norm_p_mean;
      result.convergence_tail_norm_pm1_max =
          jcp_tree_stats.convergence_tail_norm_pm1_max;
      result.convergence_tail_norm_pm1_mean =
          jcp_tree_stats.convergence_tail_norm_pm1_mean;
      result.convergence_effective_jcp_amp_max =
          jcp_tree_stats.convergence_effective_jcp_amp_max;
      result.convergence_effective_jcp_amp_mean =
          jcp_tree_stats.convergence_effective_jcp_amp_mean;
      result.convergence_estimator_high_error_caught =
          jcp_tree_stats.convergence_estimator_high_error_caught;
      result.convergence_estimator_high_error_missed =
          jcp_tree_stats.convergence_estimator_high_error_missed;
      result.convergence_estimator_high_error_caught_by_split =
          jcp_tree_stats.convergence_estimator_high_error_caught_by_split;
      result.moment_audit_nodes = jcp_tree_stats.moment_audit_nodes;
      result.moment_translation_rel_error_max =
          jcp_tree_stats.moment_translation_rel_error_max;
      result.moment_translation_rel_error_mean =
          jcp_tree_stats.moment_translation_rel_error_mean;
      result.moment_translation_offset_max =
          jcp_tree_stats.moment_translation_offset_max;
      result.moment_translation_offset_mean =
          jcp_tree_stats.moment_translation_offset_mean;
      result.moment_cancellation_ratio_max =
          jcp_tree_stats.moment_cancellation_ratio_max;
      result.moment_cancellation_ratio_mean =
          jcp_tree_stats.moment_cancellation_ratio_mean;
      result.moment_dynamic_range_max =
          jcp_tree_stats.moment_dynamic_range_max;
      result.expansion_audit = jcp_tree_stats.expansion_audit;
      result.far_expansion_audit_pairs =
          jcp_tree_stats.far_expansion_audit_pairs;
      result.far_expansion_audit_panels =
          jcp_tree_stats.far_expansion_audit_panels;
      result.far_expansion_direct_ms =
          jcp_tree_stats.far_expansion_direct_ms;
      result.far_expansion_rel_error_max =
          jcp_tree_stats.far_expansion_rel_error_max;
      result.far_expansion_rel_error_mean =
          jcp_tree_stats.far_expansion_rel_error_mean;
      result.far_expansion_mean_gate_pass =
          jcp_tree_stats.far_expansion_mean_gate_pass;
      result.far_convergence_rho_max =
          jcp_tree_stats.far_convergence_rho_max;
      result.far_convergence_rho_mean =
          jcp_tree_stats.far_convergence_rho_mean;
      result.accepted_far_aspect_ratio_max =
          jcp_tree_stats.accepted_far_aspect_ratio_max;
      result.accepted_far_aspect_ratio_mean =
          jcp_tree_stats.accepted_far_aspect_ratio_mean;
    }
    result.near_integration = jcp_tree_stats.near_integration;
    result.quad_buffer = jcp_tree_stats.quad_buffer;
    result.near_kernel_policy = jcp_tree_stats.near_kernel_policy;
    result.near_self_pairs = jcp_tree_stats.near_self_pairs;
    result.near_one_ring_pairs = jcp_tree_stats.near_one_ring_pairs;
    result.near_duffy_pairs = jcp_tree_stats.near_duffy_pairs;
    result.near_centroid_pairs = jcp_tree_stats.near_centroid_pairs;
    result.near_fallbacks = jcp_tree_stats.near_fallbacks;
    result.near_geometric_buffer_pairs =
        jcp_tree_stats.near_geometric_buffer_pairs;
    result.near_geometric_buffer_gauss_pairs =
        jcp_tree_stats.near_geometric_buffer_gauss_pairs;
    result.near_geometric_buffer_adaptive_pairs =
        jcp_tree_stats.near_geometric_buffer_adaptive_pairs;
    result.near_geometric_buffer_centroid_pairs =
        jcp_tree_stats.near_geometric_buffer_centroid_pairs;
    result.near_geometric_buffer_min_gap =
        jcp_tree_stats.near_geometric_buffer_min_gap;
    result.near_geometric_buffer_mean_gap_ratio =
        jcp_tree_stats.near_geometric_buffer_mean_gap_ratio;
    result.near_geometric_buffer_max_edge_mean =
        jcp_tree_stats.near_geometric_buffer_max_edge_mean;
    result.near_integration_ms = jcp_tree_stats.near_integration_ms;
    result.duffy_subdivisions_max = jcp_tree_stats.duffy_subdivisions_max;
    result.duffy_subdivisions_mean = jcp_tree_stats.duffy_subdivisions_mean;
    result.near_k1_accepted = jcp_tree_stats.near_k1_accepted;
    result.near_k1_rejected = jcp_tree_stats.near_k1_rejected;
    result.near_k2_accepted = jcp_tree_stats.near_k2_accepted;
    result.near_k2_rejected = jcp_tree_stats.near_k2_rejected;
    result.near_k3_accepted = jcp_tree_stats.near_k3_accepted;
    result.near_k3_rejected = jcp_tree_stats.near_k3_rejected;
    result.near_k4_accepted = jcp_tree_stats.near_k4_accepted;
    result.near_k4_rejected = jcp_tree_stats.near_k4_rejected;
    result.near_one_ring_k1_accepted =
        jcp_tree_stats.near_one_ring_k1_accepted;
    result.near_one_ring_k1_rejected =
        jcp_tree_stats.near_one_ring_k1_rejected;
    result.near_one_ring_k2_accepted =
        jcp_tree_stats.near_one_ring_k2_accepted;
    result.near_one_ring_k2_rejected =
        jcp_tree_stats.near_one_ring_k2_rejected;
    result.near_one_ring_k3_accepted =
        jcp_tree_stats.near_one_ring_k3_accepted;
    result.near_one_ring_k3_rejected =
        jcp_tree_stats.near_one_ring_k3_rejected;
    result.near_one_ring_k4_accepted =
        jcp_tree_stats.near_one_ring_k4_accepted;
    result.near_one_ring_k4_rejected =
        jcp_tree_stats.near_one_ring_k4_rejected;
    result.near_k1_policy = jcp_tree_stats.near_k1_policy;
    result.near_k2_policy = jcp_tree_stats.near_k2_policy;
    result.near_k3_policy = jcp_tree_stats.near_k3_policy;
    result.near_k4_policy = jcp_tree_stats.near_k4_policy;
    result.near_block_audit = jcp_tree_stats.near_block_audit;
    result.geometric_near_audit = jcp_tree_stats.geometric_near_audit;
    result.orthogonalization_ms = linear.orthogonalization_ms;
    result.preconditioner_ms =
        linear.preconditioner_ms + preconditioner_setup_ms;
    result.preconditioner_apply_ms = linear.preconditioner_apply_ms;
    result.preconditioner_apply_calls = linear.preconditioner_apply_calls;
    if (coupled) {
      result.phi_norm = rms_norm(std::vector<double>(
          linear.x.begin(), linear.x.begin() + panels.size()));
      result.phin_norm = rms_norm(std::vector<double>(
          linear.x.begin() + panels.size(), linear.x.end()));
    } else {
      result.phi_norm = rms_norm(linear.x);
      result.phin_norm = 0.0;
    }

    const double native_energy_scale = 4.0 * kPi;
    const double boundary_energy =
        native_energy_scale *
        (coupled ? jcp_linear_pb_energy_kcal(atoms, panels, linear.x,
                                             eps_ratio, kappa)
                 : jcp_poisson_energy_kcal(atoms, panels, linear.x,
                                           eps_ratio));
    const double born = born_reference_kcal(atoms, options.eps_in,
                                            options.eps_out, kappa);
    double energy = boundary_energy;
    std::string source = "tabi_boundary_integral";
    double selected_scale = 1.0;
    if (formulation.energy_projection == "born-check") {
      energy = born;
      source = "tabi_born_check";
      selected_scale = 0.0;
    } else if (formulation.energy_projection == "charge-potential-debug") {
      source = "tabi_charge_potential_debug";
    }
    result.boundary_energy_kcal_mol = boundary_energy;
    result.reference_energy_kcal_mol = born;
    result.selected_energy_scale = selected_scale;
    result.fallback_used = false;
    result.reaction_energy_source = source;
    if (!linear.converged || !std::isfinite(energy)) {
      result.status = "failed";
      if (!std::isfinite(energy))
        result.convergence_reason = "nan_guard";
    } else {
      result.status = "ok";
    }
    result.energy_kcal_mol = energy;
    result.energy_kj_mol = energy * kKcalToKJ;
    result.energy_kt =
        result.energy_kj_mol /
        (kKBKJPerMolK * std::max(1.0, options.temperature_k));
    result.peak_memory_bytes = result.memory_bytes;
    refresh_matvec_cost_summary(result);
    progress.write("solve_done", result, &jcp_tree_stats, result.matvec_calls,
                   result.iterations, 0, result.residual, result.status, true);
    return result;
  }

  TreeWorkspace tree;
  TreecodeStats matvec_stats;
  matvec_stats.matvec_mode = mode;
  matvec_stats.taylor_order = mode == "treecode" ? result.taylor_order : 0;
  matvec_stats.derivative_mode = result.derivative_mode;
  matvec_stats.interaction_mode = result.interaction_mode;
  matvec_stats.node_bounds_mode = result.node_bounds_mode;
  matvec_stats.tree_builder = result.tree_builder;
  matvec_stats.expansion_center = result.expansion_center;
  matvec_stats.fragmentation_safety = result.fragmentation_safety;
  matvec_stats.mac_safety = result.mac_safety;
  matvec_stats.mac_gap_factor = result.mac_gap_factor;
  matvec_stats.moment_audit = result.moment_audit;
  matvec_stats.moment_construction = result.moment_construction;
  matvec_stats.moment_sum = result.moment_sum;
  matvec_stats.convergence_audit = result.convergence_audit;
  matvec_stats.convergence_safety = result.convergence_safety;
  matvec_stats.convergence_rho_cap = result.convergence_rho_cap;
  matvec_stats.mac_aspect_cap = result.mac_aspect_cap;
  matvec_stats.convergence_estimator = result.convergence_estimator;
  matvec_stats.convergence_tail_threshold =
      result.convergence_tail_threshold;
  matvec_stats.convergence_estimator_error_cap =
      result.convergence_estimator_error_cap;
  matvec_stats.convergence_unstable_rho_cap =
      result.convergence_unstable_rho_cap;
  matvec_stats.convergence_calibration = result.convergence_calibration;
  matvec_stats.convergence_dielectric_amp =
      result.convergence_dielectric_amp;
  matvec_stats.convergence_estimator_model =
      result.convergence_estimator_model;
  matvec_stats.convergence_channel_calibration =
      result.convergence_channel_calibration;
  matvec_stats.convergence_split_policy = result.convergence_split_policy;
  matvec_stats.convergence_split_selectivity =
      result.convergence_split_selectivity;
  matvec_stats.convergence_split_absolute_floor =
      result.convergence_split_absolute_floor;
  matvec_stats.convergence_overreject_audit =
      result.convergence_overreject_audit;
  matvec_stats.convergence_overreject_detail =
      result.convergence_overreject_detail;
  matvec_stats.convergence_shape_audit = result.convergence_shape_audit;
  matvec_stats.convergence_kernel_response_audit =
      result.convergence_kernel_response_audit;
  matvec_stats.channel_missed_semantics =
      result.channel_missed_semantics;
  matvec_stats.accepted_miss_diagnostics =
      result.accepted_miss_diagnostics;
  matvec_stats.accepted_miss_absolute_floor =
      result.accepted_miss_absolute_floor;
  matvec_stats.convergence_k3_accepted_rho_cap =
      result.convergence_k3_accepted_rho_cap;
  matvec_stats.mac_distance_policy = result.mac_distance_policy;
  matvec_stats.mac_surface_margin_factor =
      result.mac_surface_margin_factor;
  matvec_stats.mac_surface_margin_policy =
      result.mac_surface_margin_policy;
  matvec_stats.mac_surface_margin_aspect_threshold =
      result.mac_surface_margin_aspect_threshold;
  matvec_stats.mac_surface_margin_rho_trigger =
      result.mac_surface_margin_rho_trigger;
  matvec_stats.mac_surface_margin_soft_aspect_threshold =
      result.mac_surface_margin_soft_aspect_threshold;
  matvec_stats.mac_surface_margin_soft_rho_trigger =
      result.mac_surface_margin_soft_rho_trigger;
  matvec_stats.mac_surface_margin_soft_factor =
      result.mac_surface_margin_soft_factor;
  matvec_stats.cluster_policy = result.cluster_policy;
  matvec_stats.target_batch_size = result.target_batch_size;
  matvec_stats.leaf_pair_min_source_panels =
      result.leaf_pair_min_source_panels;
  matvec_stats.leaf_pair_min_targets = result.leaf_pair_min_targets;
  matvec_stats.leaf_kernel_profile = result.leaf_kernel_profile;
  matvec_stats.leaf_work_profile = result.leaf_work_profile;
  matvec_stats.matvec_cost_profile = result.matvec_cost_profile;
  matvec_stats.jcp_direct_kernel = result.jcp_direct_kernel;
  matvec_stats.leaf_pair_selection_profile =
      result.leaf_pair_selection_profile;
  matvec_stats.leaf_pair_probe_cache = result.leaf_pair_probe_cache;
  matvec_stats.leaf_pair_target_group = result.leaf_pair_target_group;
  matvec_stats.leaf_pair_directness_probe =
      result.leaf_pair_directness_probe;
  matvec_stats.leaf_pair_directness_cache =
      result.leaf_pair_directness_cache;
  matvec_stats.leaf_pair_directness_verify =
      result.leaf_pair_directness_verify;
  matvec_stats.leaf_pair_scheduler = result.leaf_pair_scheduler;
  matvec_stats.leaf_pair_worklist_max_active =
      result.leaf_pair_worklist_max_active;
  matvec_stats.leaf_pair_worklist_verify =
      result.leaf_pair_worklist_verify;
  matvec_stats.leaf_block_verify = result.leaf_block_verify;
  matvec_stats.leaf_pair_max_lookahead = result.leaf_pair_max_lookahead;
  matvec_stats.expansion_audit = result.expansion_audit;
  matvec_stats.expansion_audit_channel = result.expansion_audit_channel;
  if (mode == "treecode") {
    progress.write("tree_build_start", result, &matvec_stats, 0, 0, 0, 0.0,
                   "running", true);
    tree = build_tree_workspace(
        panels, executor,
        static_cast<std::size_t>(std::max(8, options.tree_leaf_size)),
        options.tree_opening_theta, options.leaf_near_policy,
        options.leaf_panel_kernel, options.adaptive_level, "off",
        options.interaction_mode, result.node_bounds_mode,
        result.tree_builder, options.mac_safety, options.mac_gap_factor,
        options.moment_audit, result.expansion_center,
        result.fragmentation_safety, options.expansion_audit,
        options.expansion_audit_sample_count,
        options.expansion_audit_max_panels,
        options.moment_construction, options.moment_sum,
        options.convergence_audit, options.convergence_safety,
        options.convergence_rho_cap, options.mac_aspect_cap,
        options.convergence_estimator,
        options.convergence_tail_threshold,
        options.convergence_estimator_error_cap,
        options.convergence_unstable_rho_cap,
        options.convergence_calibration,
        options.cluster_policy,
        options.convergence_dielectric_amp, options.target_batch_size,
        options.leaf_pair_min_source_panels,
        options.leaf_pair_min_targets,
          options.convergence_estimator_model,
          options.expansion_audit_channel,
          options.convergence_channel_calibration,
          options.channel_missed_semantics,
        options.convergence_split_policy,
        options.accepted_miss_diagnostics,
        options.convergence_k3_accepted_rho_cap,
        options.mac_distance_policy,
        options.mac_surface_margin_factor,
        options.mac_surface_margin_policy,
        options.mac_surface_margin_aspect_threshold,
        options.mac_surface_margin_rho_trigger,
        options.accepted_miss_absolute_floor, options.leaf_pair_kernel,
        options.leaf_kernel_profile, options.leaf_work_profile,
        options.jcp_direct_kernel,
        options.leaf_pair_selection_profile,
        options.leaf_pair_probe_cache,
        options.leaf_pair_target_group,
        options.leaf_pair_max_lookahead,
        options.leaf_pair_directness_probe,
        options.leaf_pair_directness_cache,
        options.leaf_pair_directness_verify,
        options.leaf_pair_scheduler,
        options.leaf_pair_worklist_max_active,
        options.leaf_pair_worklist_verify, options.leaf_block_verify,
        options.mac_surface_margin_soft_aspect_threshold,
        options.mac_surface_margin_soft_rho_trigger,
        options.mac_surface_margin_soft_factor,
        options.convergence_split_selectivity,
        options.convergence_split_absolute_floor,
        options.convergence_overreject_audit,
        options.convergence_overreject_detail,
        options.convergence_shape_audit,
        options.convergence_kernel_response_audit,
        result.leaf_mac, options.fmm_m2l);
    tree.leaf_mac = result.leaf_mac;
    tree.taylor_order = result.taylor_order;
    tree.derivative_mode = result.derivative_mode;
    result.interaction_mode = tree.interaction_mode;
    result.node_bounds_mode = tree.node_bounds_mode;
    result.tree_builder = tree.tree_builder;
    result.expansion_center = tree.expansion_center;
    result.fragmentation_safety = tree.fragmentation_safety;
    result.mac_safety = tree.mac_safety;
    result.mac_gap_factor = tree.mac_gap_factor;
    result.moment_audit = tree.moment_audit;
    result.moment_construction = tree.moment_construction;
    result.moment_sum = tree.moment_sum;
    result.convergence_audit = tree.convergence_audit;
    result.convergence_safety = tree.convergence_safety;
    result.convergence_rho_cap = tree.convergence_rho_cap;
    result.mac_aspect_cap = tree.mac_aspect_cap;
    result.convergence_estimator = tree.convergence_estimator;
    result.convergence_tail_threshold = tree.convergence_tail_threshold;
    result.convergence_estimator_error_cap =
        tree.convergence_estimator_error_cap;
    result.convergence_unstable_rho_cap =
        tree.convergence_unstable_rho_cap;
    result.convergence_calibration = tree.convergence_calibration;
    result.convergence_dielectric_amp = tree.convergence_dielectric_amp;
    result.convergence_estimator_model = tree.convergence_estimator_model;
    result.convergence_channel_calibration =
        tree.convergence_channel_calibration;
    result.convergence_split_policy = tree.convergence_split_policy;
    result.convergence_split_selectivity =
        tree.convergence_split_selectivity;
    result.convergence_split_absolute_floor =
        tree.convergence_split_absolute_floor;
    result.convergence_overreject_audit = tree.convergence_overreject_audit;
    result.convergence_overreject_detail = tree.convergence_overreject_detail;
    result.convergence_shape_audit = tree.convergence_shape_audit;
    result.convergence_kernel_response_audit =
        tree.convergence_kernel_response_audit;
    result.channel_missed_semantics = tree.channel_missed_semantics;
    result.accepted_miss_diagnostics = tree.accepted_miss_diagnostics;
    result.accepted_miss_absolute_floor = tree.accepted_miss_absolute_floor;
    result.convergence_k3_accepted_rho_cap =
        tree.convergence_k3_accepted_rho_cap;
    result.mac_distance_policy = tree.mac_distance_policy;
    result.mac_surface_margin_factor = tree.mac_surface_margin_factor;
    result.mac_surface_margin_policy = tree.mac_surface_margin_policy;
    result.mac_surface_margin_aspect_threshold =
        tree.mac_surface_margin_aspect_threshold;
    result.mac_surface_margin_rho_trigger =
        tree.mac_surface_margin_rho_trigger;
    result.mac_surface_margin_soft_aspect_threshold =
        tree.mac_surface_margin_soft_aspect_threshold;
    result.mac_surface_margin_soft_rho_trigger =
        tree.mac_surface_margin_soft_rho_trigger;
    result.mac_surface_margin_soft_factor =
        tree.mac_surface_margin_soft_factor;
    result.cluster_policy = tree.cluster_policy;
    result.target_batch_size = static_cast<int>(tree.target_batch_size);
    result.leaf_pair_min_source_panels =
        static_cast<int>(tree.leaf_pair_min_source_panels);
    result.leaf_pair_min_targets =
        static_cast<int>(tree.leaf_pair_min_targets);
    result.leaf_kernel_profile = tree.leaf_kernel_profile;
    result.leaf_work_profile = tree.leaf_work_profile;
    result.jcp_direct_kernel = tree.jcp_direct_kernel;
    result.leaf_pair_selection_profile = tree.leaf_pair_selection_profile;
    result.leaf_pair_probe_cache = tree.leaf_pair_probe_cache;
    result.leaf_pair_target_group = tree.leaf_pair_target_group;
    result.leaf_pair_directness_probe = tree.leaf_pair_directness_probe;
    result.leaf_pair_directness_cache = tree.leaf_pair_directness_cache;
    result.leaf_pair_directness_verify = tree.leaf_pair_directness_verify;
    result.leaf_pair_scheduler = tree.leaf_pair_scheduler;
    result.leaf_pair_worklist_max_active =
        static_cast<int>(tree.leaf_pair_worklist_max_active);
    result.leaf_pair_worklist_verify = tree.leaf_pair_worklist_verify;
    result.leaf_block_verify = tree.leaf_block_verify;
    result.leaf_pair_max_lookahead =
        static_cast<int>(tree.leaf_pair_max_lookahead);
    result.expansion_audit = tree.expansion_audit;
    result.expansion_audit_channel = tree.expansion_audit_channel;
    matvec_stats.node_bounds_mode = tree.node_bounds_mode;
    matvec_stats.interaction_mode = tree.interaction_mode;
    matvec_stats.tree_builder = tree.tree_builder;
    matvec_stats.expansion_center = tree.expansion_center;
    matvec_stats.fragmentation_safety = tree.fragmentation_safety;
    matvec_stats.mac_safety = tree.mac_safety;
    matvec_stats.mac_gap_factor = tree.mac_gap_factor;
    matvec_stats.moment_audit = tree.moment_audit;
    matvec_stats.moment_construction = tree.moment_construction;
    matvec_stats.moment_sum = tree.moment_sum;
    matvec_stats.convergence_audit = tree.convergence_audit;
    matvec_stats.convergence_safety = tree.convergence_safety;
    matvec_stats.convergence_rho_cap = tree.convergence_rho_cap;
    matvec_stats.mac_aspect_cap = tree.mac_aspect_cap;
    matvec_stats.convergence_estimator = tree.convergence_estimator;
    matvec_stats.convergence_tail_threshold =
        tree.convergence_tail_threshold;
    matvec_stats.convergence_estimator_error_cap =
        tree.convergence_estimator_error_cap;
    matvec_stats.convergence_unstable_rho_cap =
        tree.convergence_unstable_rho_cap;
    matvec_stats.convergence_calibration = tree.convergence_calibration;
    matvec_stats.convergence_dielectric_amp =
        tree.convergence_dielectric_amp;
    matvec_stats.convergence_estimator_model =
        tree.convergence_estimator_model;
    matvec_stats.convergence_channel_calibration =
        tree.convergence_channel_calibration;
    matvec_stats.convergence_split_policy = tree.convergence_split_policy;
    matvec_stats.convergence_split_selectivity =
        tree.convergence_split_selectivity;
    matvec_stats.convergence_split_absolute_floor =
        tree.convergence_split_absolute_floor;
    matvec_stats.convergence_overreject_audit =
        tree.convergence_overreject_audit;
    matvec_stats.convergence_overreject_detail =
        tree.convergence_overreject_detail;
    matvec_stats.convergence_shape_audit = tree.convergence_shape_audit;
    matvec_stats.convergence_kernel_response_audit =
        tree.convergence_kernel_response_audit;
    matvec_stats.channel_missed_semantics = tree.channel_missed_semantics;
    matvec_stats.accepted_miss_diagnostics =
        tree.accepted_miss_diagnostics;
    matvec_stats.accepted_miss_absolute_floor =
        tree.accepted_miss_absolute_floor;
    matvec_stats.convergence_k3_accepted_rho_cap =
        tree.convergence_k3_accepted_rho_cap;
    matvec_stats.mac_distance_policy = tree.mac_distance_policy;
    matvec_stats.mac_surface_margin_factor =
        tree.mac_surface_margin_factor;
    matvec_stats.mac_surface_margin_policy =
        tree.mac_surface_margin_policy;
    matvec_stats.mac_surface_margin_aspect_threshold =
        tree.mac_surface_margin_aspect_threshold;
    matvec_stats.mac_surface_margin_rho_trigger =
        tree.mac_surface_margin_rho_trigger;
    matvec_stats.mac_surface_margin_soft_aspect_threshold =
        tree.mac_surface_margin_soft_aspect_threshold;
    matvec_stats.mac_surface_margin_soft_rho_trigger =
        tree.mac_surface_margin_soft_rho_trigger;
    matvec_stats.mac_surface_margin_soft_factor =
        tree.mac_surface_margin_soft_factor;
    matvec_stats.cluster_policy = tree.cluster_policy;
    matvec_stats.target_batch_size = static_cast<int>(tree.target_batch_size);
    matvec_stats.leaf_pair_min_source_panels =
        static_cast<int>(tree.leaf_pair_min_source_panels);
    matvec_stats.leaf_pair_min_targets =
        static_cast<int>(tree.leaf_pair_min_targets);
    matvec_stats.leaf_kernel_profile = tree.leaf_kernel_profile;
    matvec_stats.leaf_work_profile = tree.leaf_work_profile;
    matvec_stats.matvec_cost_profile = result.matvec_cost_profile;
    matvec_stats.jcp_direct_kernel = tree.jcp_direct_kernel;
    matvec_stats.leaf_pair_selection_profile =
        tree.leaf_pair_selection_profile;
    matvec_stats.leaf_pair_probe_cache = tree.leaf_pair_probe_cache;
    matvec_stats.leaf_pair_target_group = tree.leaf_pair_target_group;
    matvec_stats.leaf_pair_directness_probe =
        tree.leaf_pair_directness_probe;
    matvec_stats.leaf_pair_directness_cache =
        tree.leaf_pair_directness_cache;
    matvec_stats.leaf_pair_directness_verify =
        tree.leaf_pair_directness_verify;
    matvec_stats.leaf_pair_scheduler = tree.leaf_pair_scheduler;
    matvec_stats.leaf_pair_worklist_max_active =
        static_cast<int>(tree.leaf_pair_worklist_max_active);
    matvec_stats.leaf_pair_worklist_verify = tree.leaf_pair_worklist_verify;
    matvec_stats.leaf_block_verify = tree.leaf_block_verify;
    matvec_stats.leaf_pair_max_lookahead =
        static_cast<int>(tree.leaf_pair_max_lookahead);
    matvec_stats.expansion_audit = tree.expansion_audit;
    matvec_stats.expansion_audit_channel = tree.expansion_audit_channel;
    matvec_stats.tree_depth_max = tree.depth_max;
    matvec_stats.node_radius_max = tree.node_radius_max;
    matvec_stats.node_radius_mean = tree.node_radius_mean;
    matvec_stats.node_aspect_ratio_max = tree.node_aspect_ratio_max;
    matvec_stats.node_aspect_ratio_mean = tree.node_aspect_ratio_mean;
    matvec_stats.tree_build_ms = tree.build_ms;
    matvec_stats.tree_nodes = tree.nodes.size();
    matvec_stats.tree_leaf_count = tree.leaf_count;
    result.tree_build_ms = tree.build_ms;
    result.tree_nodes = tree.nodes.size();
    result.tree_leaf_count = tree.leaf_count;
    result.tree_depth_max = tree.depth_max;
    result.memory_bytes += tree.nodes.size() * sizeof(TreeNode) +
                           tree.child_indices.size() * sizeof(std::size_t) +
                           tree.order.size() * sizeof(std::size_t) +
                           tree.morton.size() * sizeof(std::uint64_t) +
                           panels.size() * 14u * sizeof(double) +
                           panels.size() * sizeof(std::size_t);
    progress.write("tree_build_end", result, &matvec_stats, 0, 0, 0, 0.0,
                   "running", true);
  }

  const auto rhs_begin = std::chrono::steady_clock::now();
  const std::vector<double> rhs =
      build_rhs(atoms, panels, contrast, kappa, formulation.rhs_scale);
  result.rhs_ms = elapsed_ms(rhs_begin, std::chrono::steady_clock::now());
  const int max_iter =
      std::min(std::max(1, options.max_iter),
               std::max(32, static_cast<int>(panels.size()) * 4));
  const double tol = std::max(1.0e-12, options.residual_tol);
  double matvec_ms_total = 0.0;
  std::size_t matvec_calls = 0;
  auto matvec = [&](const std::vector<double> &x, std::vector<double> &out) {
    const auto begin = std::chrono::steady_clock::now();
    if (mode == "treecode")
      apply_treecode_operator(panels, x, operator_contrast, kappa, tree,
                              executor, out, &matvec_stats);
    else
      apply_direct_operator(panels, x, operator_contrast, kappa, executor, out,
                            &matvec_stats);
    matvec_ms_total += elapsed_ms(begin, std::chrono::steady_clock::now());
    ++matvec_calls;
  };
  const auto gmres_begin = std::chrono::steady_clock::now();
  progress.write("gmres_start", result, &matvec_stats, matvec_calls, 0, 0,
                 0.0, "running", true);
  auto legacy_progress = [&](int gmres_iter, int restart_iter,
                             double true_residual_rel,
                             double preconditioned_residual_rel, bool force,
                             const std::string &phase,
                             const LinearSolveResult &state) {
    result.true_residual_probe_scope = state.true_residual_probe_scope;
    result.true_residual_sample_count = state.true_residual_sample_count;
    result.true_residual_probe_calls = state.true_residual_probe_calls;
    result.true_residual_probe_ms = state.true_residual_probe_ms;
    result.sampled_true_residual_rel = state.sampled_true_residual_rel;
    result.sampled_true_residual_rows = state.sampled_true_residual_rows;
    result.sampled_probe_calls = state.sampled_probe_calls;
    result.sampled_probe_ms = state.sampled_probe_ms;
    result.full_true_residual_confirm_calls =
        state.full_true_residual_confirm_calls;
    result.full_true_residual_confirm_ms =
        state.full_true_residual_confirm_ms;
    if (phase != "gmres_sampled_true_probe" &&
        std::isfinite(true_residual_rel)) {
      result.final_residual = true_residual_rel;
      result.true_residual_rel = true_residual_rel;
      if (result.best_residual <= 0.0 ||
          true_residual_rel < result.best_residual) {
        result.best_residual = true_residual_rel;
        result.best_residual_iter = gmres_iter;
      }
    }
    if (std::isfinite(preconditioned_residual_rel)) {
      result.final_preconditioned_residual = preconditioned_residual_rel;
      if (result.best_preconditioned_residual <= 0.0 ||
          preconditioned_residual_rel <
              result.best_preconditioned_residual) {
        result.best_preconditioned_residual = preconditioned_residual_rel;
      }
    }
    progress.write(phase, result, &matvec_stats, matvec_calls, gmres_iter,
                   restart_iter, true_residual_rel, "running", force);
  };
  std::vector<double> legacy_diagonal(panels.size(), 1.0);
  for (std::size_t i = 0; i < panels.size(); ++i)
    legacy_diagonal[i] = panels[i].diagonal;
  LinearSolveResult linear =
      gmres_solve(matvec, rhs, legacy_diagonal, max_iter, tol,
                  options.gmres_restart, executor, {}, legacy_progress,
                  result.true_residual_probe_mode,
                  result.true_residual_probe_interval,
                  result.true_residual_probe_scope,
                  result.true_residual_sample_count);
  result.gmres_ms = elapsed_ms(gmres_begin, std::chrono::steady_clock::now());
  result.iterations = linear.iterations;
  result.residual = linear.residual;
  result.gmres_restart_requested = linear.restart_requested;
  result.gmres_restart_effective = linear.restart_effective;
  result.best_residual = linear.best_residual;
  result.best_residual_iter = linear.best_residual_iter;
  result.final_residual = linear.final_residual;
  result.best_preconditioned_residual = linear.best_preconditioned_residual;
  result.final_preconditioned_residual =
      linear.final_preconditioned_residual;
  result.true_residual_abs = linear.true_residual_abs;
  result.true_residual_rel = linear.true_residual_rel;
  result.true_residual_probe_mode = linear.true_residual_probe_mode;
  result.true_residual_probe_scope = linear.true_residual_probe_scope;
  result.true_residual_probe_interval = linear.true_residual_probe_interval;
  result.true_residual_sample_count = linear.true_residual_sample_count;
  result.true_residual_probe_calls = linear.true_residual_probe_calls;
  result.true_residual_probe_ms = linear.true_residual_probe_ms;
  result.sampled_true_residual_rel = linear.sampled_true_residual_rel;
  result.sampled_true_residual_rows = linear.sampled_true_residual_rows;
  result.sampled_probe_calls = linear.sampled_probe_calls;
  result.sampled_probe_ms = linear.sampled_probe_ms;
  result.full_true_residual_confirm_calls =
      linear.full_true_residual_confirm_calls;
  result.full_true_residual_confirm_ms =
      linear.full_true_residual_confirm_ms;
  result.converged = linear.converged;
  result.convergence_reason = linear.convergence_reason;
  result.residual_history = residual_history_string(linear.residual_history);
  result.matvec_ms_total = matvec_ms_total;
  result.matvec_calls = matvec_calls;
  result.matvec_ms_avg =
      matvec_calls > 0 ? matvec_ms_total / static_cast<double>(matvec_calls)
                       : 0.0;
  result.orthogonalization_ms = linear.orthogonalization_ms;
  result.preconditioner_ms = linear.preconditioner_ms;
  result.tree_build_ms = matvec_stats.tree_build_ms;
  result.tree_nodes = matvec_stats.tree_nodes;
  result.tree_leaf_count = matvec_stats.tree_leaf_count;
  result.near_interactions = matvec_stats.near_interactions;
  result.far_interactions = matvec_stats.far_interactions;
  result.singular_near_interactions = matvec_stats.singular_near_interactions;
  result.leaf_near_interactions = matvec_stats.leaf_near_interactions;
  result.far_node_interactions = matvec_stats.far_node_interactions;
  result.leaf_panel_interactions = matvec_stats.leaf_panel_interactions;
  result.leaf_tile_interactions = matvec_stats.leaf_tile_interactions;
  result.leaf_tile_panel_interactions =
      matvec_stats.leaf_tile_panel_interactions;
  result.leaf_residual_panel_interactions =
      matvec_stats.leaf_residual_panel_interactions;
  result.smooth_leaf_node_interactions =
      matvec_stats.smooth_leaf_node_interactions;
  result.leaf_pruned_interactions = matvec_stats.leaf_pruned_interactions;
  result.leaf_mac_candidates = matvec_stats.leaf_mac_candidates;
  result.leaf_mac_accepted = matvec_stats.leaf_mac_accepted;
  result.leaf_mac_rejected = matvec_stats.leaf_mac_rejected;
  result.leaf_mac_far_panels = matvec_stats.leaf_mac_far_panels;
  result.leaf_mac_direct_panels_saved_estimate =
      matvec_stats.leaf_mac_direct_panels_saved_estimate;
  result.leaf_work_profile = matvec_stats.leaf_work_profile;
  result.leaf_work_mac_terminal_panels =
      matvec_stats.leaf_work_mac_terminal_panels;
  result.leaf_work_physical_gap_panels =
      matvec_stats.leaf_work_physical_gap_panels;
  result.leaf_work_fragmentation_panels =
      matvec_stats.leaf_work_fragmentation_panels;
  result.leaf_work_convergence_panels =
      matvec_stats.leaf_work_convergence_panels;
  result.leaf_work_convergence_estimator_panels =
      matvec_stats.leaf_work_convergence_estimator_panels;
  result.leaf_work_convergence_rho_panels =
      matvec_stats.leaf_work_convergence_rho_panels;
  result.leaf_work_convergence_aspect_panels =
      matvec_stats.leaf_work_convergence_aspect_panels;
  result.leaf_work_surface_margin_panels =
      matvec_stats.leaf_work_surface_margin_panels;
  result.leaf_work_surface_margin_hard_panels =
      matvec_stats.leaf_work_surface_margin_hard_panels;
  result.leaf_work_surface_margin_soft_panels =
      matvec_stats.leaf_work_surface_margin_soft_panels;
  result.leaf_work_k3_rho_cap_panels =
      matvec_stats.leaf_work_k3_rho_cap_panels;
  result.leaf_work_singular_fallback_panels =
      matvec_stats.leaf_work_singular_fallback_panels;
  result.leaf_work_scalar_tail_panels =
      matvec_stats.leaf_work_scalar_tail_panels;
  result.leaf_work_profile_accounted_panels =
      matvec_stats.leaf_work_profile_accounted_panels;
  result.leaf_work_mac_terminal_frac =
      matvec_stats.leaf_work_mac_terminal_frac;
  result.leaf_work_convergence_frac =
      matvec_stats.leaf_work_convergence_frac;
  result.leaf_work_surface_margin_frac =
      matvec_stats.leaf_work_surface_margin_frac;
  result.leaf_work_singular_fallback_frac =
      matvec_stats.leaf_work_singular_fallback_frac;
  result.leaf_work_scalar_tail_frac =
      matvec_stats.leaf_work_scalar_tail_frac;
  result.adaptive_smooth_leaf_accepted =
      matvec_stats.adaptive_smooth_leaf_accepted;
  result.adaptive_smooth_leaf_rejected =
      matvec_stats.adaptive_smooth_leaf_rejected;
  result.adaptive_span_tile_accepted = matvec_stats.adaptive_span_tile_accepted;
  result.adaptive_span_tile_rejected = matvec_stats.adaptive_span_tile_rejected;
  result.adaptive_exact_leaf_fallbacks =
      matvec_stats.adaptive_exact_leaf_fallbacks;
  result.adaptive_reject_one_ring_overlap =
      matvec_stats.adaptive_reject_one_ring_overlap;
  result.adaptive_reject_error_bound_failed =
      matvec_stats.adaptive_reject_error_bound_failed;
  result.adaptive_reject_p4_check_failed =
      matvec_stats.adaptive_reject_p4_check_failed;
  result.adaptive_reject_too_close = matvec_stats.adaptive_reject_too_close;
  result.leaf_work_reduction_ratio = matvec_stats.leaf_work_reduction_ratio;
  result.target_cluster_interactions = matvec_stats.target_cluster_interactions;
  result.target_cluster_targets = matvec_stats.target_cluster_targets;
  result.target_cluster_panels_pruned =
      matvec_stats.target_cluster_panels_pruned;
  result.target_cluster_rejected = matvec_stats.target_cluster_rejected;
  result.target_cluster_eval_ms = matvec_stats.target_cluster_eval_ms;
  result.leaf_pair_batch_interactions =
      matvec_stats.leaf_pair_batch_interactions;
  result.leaf_pair_batch_targets = matvec_stats.leaf_pair_batch_targets;
  result.leaf_pair_batch_source_panels =
      matvec_stats.leaf_pair_batch_source_panels;
  result.leaf_pair_batch_residual_panels =
      matvec_stats.leaf_pair_batch_residual_panels;
  result.leaf_pair_batch_fallbacks = matvec_stats.leaf_pair_batch_fallbacks;
  result.leaf_pair_batch_ms = matvec_stats.leaf_pair_batch_ms;
  result.leaf_pair_batch_kernel = matvec_stats.leaf_pair_batch_kernel;
  result.leaf_pair_batch_vectorized_sources =
      matvec_stats.leaf_pair_batch_vectorized_sources;
  result.leaf_pair_batch_vector_lanes =
      matvec_stats.leaf_pair_batch_vector_lanes;
  result.leaf_pair_batch_pairs_per_second =
      matvec_stats.leaf_pair_batch_pairs_per_second;
  result.leaf_kernel_profile = matvec_stats.leaf_kernel_profile;
  result.leaf_work_profile = matvec_stats.leaf_work_profile;
  result.jcp_direct_kernel = matvec_stats.jcp_direct_kernel;
  result.leaf_kernel_scalar_tasks = matvec_stats.leaf_kernel_scalar_tasks;
  result.leaf_kernel_batch_tasks = matvec_stats.leaf_kernel_batch_tasks;
  result.leaf_kernel_source_panel_evals =
      matvec_stats.leaf_kernel_source_panel_evals;
  result.leaf_kernel_scalar_ms = matvec_stats.leaf_kernel_scalar_ms;
  result.leaf_kernel_target_avx2_ms =
      matvec_stats.leaf_kernel_target_avx2_ms;
  result.leaf_kernel_tail_ms = matvec_stats.leaf_kernel_tail_ms;
  result.leaf_pair_batch_soa_sources =
      matvec_stats.leaf_pair_batch_soa_sources;
  result.leaf_pair_batch_soa_ms = matvec_stats.leaf_pair_batch_soa_ms;
  result.leaf_pair_batch_scalar_tail_ms =
      matvec_stats.leaf_pair_batch_scalar_tail_ms;
  result.leaf_pair_batch_clean_source_panels_per_second =
      matvec_stats.leaf_pair_batch_clean_source_panels_per_second;
  result.leaf_pair_selection_profile =
      matvec_stats.leaf_pair_selection_profile;
  result.leaf_pair_probe_cache = matvec_stats.leaf_pair_probe_cache;
  result.leaf_pair_target_group = matvec_stats.leaf_pair_target_group;
  result.leaf_pair_directness_probe =
      matvec_stats.leaf_pair_directness_probe;
  result.leaf_pair_directness_cache =
      matvec_stats.leaf_pair_directness_cache;
  result.leaf_pair_directness_verify =
      matvec_stats.leaf_pair_directness_verify;
  result.leaf_pair_max_lookahead = matvec_stats.leaf_pair_max_lookahead;
  result.leaf_pair_selection_candidate_scans =
      matvec_stats.leaf_pair_selection_candidate_scans;
  result.leaf_pair_selection_clean_rejects =
      matvec_stats.leaf_pair_selection_clean_rejects;
  result.leaf_pair_selection_probe_calls =
      matvec_stats.leaf_pair_selection_probe_calls;
  result.leaf_pair_selection_probe_node_visits =
      matvec_stats.leaf_pair_selection_probe_node_visits;
  result.leaf_pair_selection_probe_ms =
      matvec_stats.leaf_pair_selection_probe_ms;
  result.leaf_pair_selection_ms = matvec_stats.leaf_pair_selection_ms;
  result.leaf_pair_selection_probe_cache_hits =
      matvec_stats.leaf_pair_selection_probe_cache_hits;
  result.leaf_pair_selection_probe_cache_misses =
      matvec_stats.leaf_pair_selection_probe_cache_misses;
  result.leaf_pair_selection_insufficient_targets =
      matvec_stats.leaf_pair_selection_insufficient_targets;
  result.leaf_pair_selection_same_target_leaf_rejects =
      matvec_stats.leaf_pair_selection_same_target_leaf_rejects;
  result.leaf_pair_directness_ancestor_probe_calls =
      matvec_stats.leaf_pair_directness_ancestor_probe_calls;
  result.leaf_pair_directness_ancestor_probe_node_visits =
      matvec_stats.leaf_pair_directness_ancestor_probe_node_visits;
  result.leaf_pair_directness_ancestor_probe_ms =
      matvec_stats.leaf_pair_directness_ancestor_probe_ms;
  result.leaf_pair_directness_pair_cache_hits =
      matvec_stats.leaf_pair_directness_pair_cache_hits;
  result.leaf_pair_directness_pair_cache_misses =
      matvec_stats.leaf_pair_directness_pair_cache_misses;
  result.leaf_pair_directness_all_direct_proofs =
      matvec_stats.leaf_pair_directness_all_direct_proofs;
  result.leaf_pair_directness_mixed_fallbacks =
      matvec_stats.leaf_pair_directness_mixed_fallbacks;
  result.leaf_pair_directness_verify_checks =
      matvec_stats.leaf_pair_directness_verify_checks;
  result.leaf_pair_directness_verify_mismatches =
      matvec_stats.leaf_pair_directness_verify_mismatches;
  result.leaf_pair_scheduler = matvec_stats.leaf_pair_scheduler;
  result.leaf_pair_worklist_max_active =
      matvec_stats.leaf_pair_worklist_max_active;
  result.leaf_pair_worklist_verify = matvec_stats.leaf_pair_worklist_verify;
  result.leaf_pair_worklist_tasks = matvec_stats.leaf_pair_worklist_tasks;
  result.leaf_pair_worklist_groups = matvec_stats.leaf_pair_worklist_groups;
  result.leaf_pair_worklist_batched_groups =
      matvec_stats.leaf_pair_worklist_batched_groups;
  result.leaf_pair_worklist_scalar_flush_tasks =
      matvec_stats.leaf_pair_worklist_scalar_flush_tasks;
  result.leaf_pair_worklist_forced_flushes =
      matvec_stats.leaf_pair_worklist_forced_flushes;
  result.leaf_pair_worklist_active_groups_max =
      matvec_stats.leaf_pair_worklist_active_groups_max;
  result.leaf_pair_worklist_ms = matvec_stats.leaf_pair_worklist_ms;
  result.leaf_pair_worklist_verify_checks =
      matvec_stats.leaf_pair_worklist_verify_checks;
  result.leaf_pair_worklist_verify_mismatches =
      matvec_stats.leaf_pair_worklist_verify_mismatches;
  result.leaf_block_verify = matvec_stats.leaf_block_verify;
  result.leaf_block_target_leaves = matvec_stats.leaf_block_target_leaves;
  result.leaf_block_targets = matvec_stats.leaf_block_targets;
  result.leaf_block_source_groups = matvec_stats.leaf_block_source_groups;
  result.leaf_block_batched_groups = matvec_stats.leaf_block_batched_groups;
  result.leaf_block_scalar_tail_tasks =
      matvec_stats.leaf_block_scalar_tail_tasks;
  result.leaf_block_group_ms = matvec_stats.leaf_block_group_ms;
  result.leaf_block_flush_ms = matvec_stats.leaf_block_flush_ms;
  result.leaf_block_writeback_ms = matvec_stats.leaf_block_writeback_ms;
  result.leaf_block_verify_checks = matvec_stats.leaf_block_verify_checks;
  result.leaf_block_verify_mismatches =
      matvec_stats.leaf_block_verify_mismatches;
  result.adaptive_level = matvec_stats.adaptive_level;
  result.cluster_policy = matvec_stats.cluster_policy;
  result.target_batch_size = matvec_stats.target_batch_size;
  result.leaf_pair_min_source_panels =
      matvec_stats.leaf_pair_min_source_panels;
  result.leaf_pair_min_targets = matvec_stats.leaf_pair_min_targets;
  result.leaf_near_policy = matvec_stats.leaf_near_policy;
  result.leaf_panel_kernel = matvec_stats.leaf_panel_kernel;
  result.interaction_mode = matvec_stats.interaction_mode;
  result.node_bounds_mode = matvec_stats.node_bounds_mode;
  result.tree_builder = matvec_stats.tree_builder;
  result.expansion_center = matvec_stats.expansion_center;
  result.fragmentation_safety = matvec_stats.fragmentation_safety;
  result.mac_safety = matvec_stats.mac_safety;
  result.mac_gap_factor = matvec_stats.mac_gap_factor;
  result.tree_depth_max = matvec_stats.tree_depth_max;
  result.mac_reject_radius_limited = matvec_stats.mac_reject_radius_limited;
  result.mac_accept_rate = matvec_stats.mac_accept_rate;
  result.node_radius_max = matvec_stats.node_radius_max;
  result.node_radius_mean = matvec_stats.node_radius_mean;
  result.node_aspect_ratio_max = matvec_stats.node_aspect_ratio_max;
  result.node_aspect_ratio_mean = matvec_stats.node_aspect_ratio_mean;
  result.leaf_panel_work_per_matvec =
      matvec_stats.leaf_panel_work_per_matvec;
  result.far_node_work_per_matvec =
      matvec_stats.far_node_work_per_matvec;
  result.far_geometric_unsafe_count =
      matvec_stats.far_geometric_unsafe_count;
  result.far_geometric_unsafe_panels =
      matvec_stats.far_geometric_unsafe_panels;
  result.far_geometric_min_gap = matvec_stats.far_geometric_min_gap;
  result.far_geometric_gap_ratio_min =
      matvec_stats.far_geometric_gap_ratio_min;
  result.far_geometric_gap_ratio_mean =
      matvec_stats.far_geometric_gap_ratio_mean;
  result.far_geometric_audit_samples =
      matvec_stats.far_geometric_audit_samples;
  result.mac_physical_gap_rejects =
      matvec_stats.mac_physical_gap_rejects;
  result.mac_convergence_rho_rejects =
      matvec_stats.mac_convergence_rho_rejects;
  result.mac_convergence_aspect_rejects =
      matvec_stats.mac_convergence_aspect_rejects;
  result.mac_convergence_rho_reject_panels =
      matvec_stats.mac_convergence_rho_reject_panels;
  result.mac_convergence_aspect_reject_panels =
      matvec_stats.mac_convergence_aspect_reject_panels;
  result.fragmented_far_count = matvec_stats.fragmented_far_count;
  result.fragmented_far_panels = matvec_stats.fragmented_far_panels;
  result.fragmented_split_rejects = matvec_stats.fragmented_split_rejects;
  result.fragmentation_center_gap_ratio_max =
      matvec_stats.fragmentation_center_gap_ratio_max;
  result.fragmentation_center_gap_ratio_mean =
      matvec_stats.fragmentation_center_gap_ratio_mean;
  result.moment_audit = matvec_stats.moment_audit;
  result.moment_construction = matvec_stats.moment_construction;
  result.moment_sum = matvec_stats.moment_sum;
  result.convergence_audit = matvec_stats.convergence_audit;
  result.convergence_safety = matvec_stats.convergence_safety;
  result.convergence_rho_cap = matvec_stats.convergence_rho_cap;
  result.mac_aspect_cap = matvec_stats.mac_aspect_cap;
  result.convergence_estimator = matvec_stats.convergence_estimator;
  result.convergence_tail_threshold =
      matvec_stats.convergence_tail_threshold;
  result.convergence_estimator_error_cap =
      matvec_stats.convergence_estimator_error_cap;
  result.convergence_unstable_rho_cap =
      matvec_stats.convergence_unstable_rho_cap;
  result.convergence_calibration = matvec_stats.convergence_calibration;
  result.convergence_dielectric_amp =
      matvec_stats.convergence_dielectric_amp;
  result.convergence_estimator_model =
      matvec_stats.convergence_estimator_model;
  result.convergence_channel_calibration =
      matvec_stats.convergence_channel_calibration;
  result.convergence_split_policy = matvec_stats.convergence_split_policy;
  result.convergence_split_selectivity =
      matvec_stats.convergence_split_selectivity;
  result.convergence_split_absolute_floor =
      matvec_stats.convergence_split_absolute_floor;
  result.convergence_overreject_audit =
      matvec_stats.convergence_overreject_audit;
  result.convergence_overreject_detail =
      matvec_stats.convergence_overreject_detail;
  result.convergence_shape_audit = matvec_stats.convergence_shape_audit;
  result.convergence_kernel_response_audit =
      matvec_stats.convergence_kernel_response_audit;
  result.channel_missed_semantics = matvec_stats.channel_missed_semantics;
  result.accepted_miss_diagnostics =
      matvec_stats.accepted_miss_diagnostics;
  result.accepted_miss_absolute_floor =
      matvec_stats.accepted_miss_absolute_floor;
  result.convergence_k3_accepted_rho_cap =
      matvec_stats.convergence_k3_accepted_rho_cap;
  result.mac_distance_policy = matvec_stats.mac_distance_policy;
  result.mac_surface_margin_factor =
      matvec_stats.mac_surface_margin_factor;
  result.mac_surface_margin_policy =
      matvec_stats.mac_surface_margin_policy;
  result.mac_surface_margin_aspect_threshold =
      matvec_stats.mac_surface_margin_aspect_threshold;
  result.mac_surface_margin_rho_trigger =
      matvec_stats.mac_surface_margin_rho_trigger;
  result.mac_surface_margin_soft_aspect_threshold =
      matvec_stats.mac_surface_margin_soft_aspect_threshold;
  result.mac_surface_margin_soft_rho_trigger =
      matvec_stats.mac_surface_margin_soft_rho_trigger;
  result.mac_surface_margin_soft_factor =
      matvec_stats.mac_surface_margin_soft_factor;
  result.convergence_estimator_samples =
      matvec_stats.convergence_estimator_samples;
  result.convergence_estimator_unstable_far_count =
      matvec_stats.convergence_estimator_unstable_far_count;
  result.convergence_estimator_unstable_far_panels =
      matvec_stats.convergence_estimator_unstable_far_panels;
  result.convergence_estimator_split_rejects =
      matvec_stats.convergence_estimator_split_rejects;
  result.convergence_estimator_split_reject_panels =
      matvec_stats.convergence_estimator_split_reject_panels;
  result.convergence_low_abs_reaccepted_far_count =
      matvec_stats.convergence_low_abs_reaccepted_far_count;
  result.convergence_low_abs_reaccepted_far_panels =
      matvec_stats.convergence_low_abs_reaccepted_far_panels;
  result.convergence_estimator_tail_ratio_max =
      matvec_stats.convergence_estimator_tail_ratio_max;
  result.convergence_estimator_tail_ratio_mean =
      matvec_stats.convergence_estimator_tail_ratio_mean;
  result.convergence_estimator_error_max =
      matvec_stats.convergence_estimator_error_max;
  result.convergence_estimator_error_mean =
      matvec_stats.convergence_estimator_error_mean;
  result.convergence_estimator_abs_error_max =
      matvec_stats.convergence_estimator_abs_error_max;
  result.convergence_estimator_abs_error_mean =
      matvec_stats.convergence_estimator_abs_error_mean;
  result.convergence_tail_norm_p_max =
      matvec_stats.convergence_tail_norm_p_max;
  result.convergence_tail_norm_p_mean =
      matvec_stats.convergence_tail_norm_p_mean;
  result.convergence_tail_norm_pm1_max =
      matvec_stats.convergence_tail_norm_pm1_max;
  result.convergence_tail_norm_pm1_mean =
      matvec_stats.convergence_tail_norm_pm1_mean;
  result.convergence_effective_jcp_amp_max =
      matvec_stats.convergence_effective_jcp_amp_max;
  result.convergence_effective_jcp_amp_mean =
      matvec_stats.convergence_effective_jcp_amp_mean;
  result.convergence_dielectric_amp_max =
      matvec_stats.convergence_dielectric_amp_max;
  result.convergence_dielectric_amp_mean =
      matvec_stats.convergence_dielectric_amp_mean;
  result.convergence_channel_estimated_error_max =
      matvec_stats.convergence_channel_estimated_error_max;
  result.convergence_channel_estimated_error_mean =
      matvec_stats.convergence_channel_estimated_error_mean;
  result.convergence_channel_effective_jcp_amp_max =
      matvec_stats.convergence_channel_effective_jcp_amp_max;
  result.convergence_channel_effective_jcp_amp_mean =
      matvec_stats.convergence_channel_effective_jcp_amp_mean;
  result.convergence_channel_unstable_far_count =
      matvec_stats.convergence_channel_unstable_far_count;
  result.convergence_channel_unstable_far_panels =
      matvec_stats.convergence_channel_unstable_far_panels;
  result.convergence_channel_split_rejects =
      matvec_stats.convergence_channel_split_rejects;
  result.convergence_channel_split_reject_panels =
      matvec_stats.convergence_channel_split_reject_panels;
  result.convergence_channel_low_abs_reaccepted_count =
      matvec_stats.convergence_channel_low_abs_reaccepted_count;
  result.convergence_channel_low_abs_reaccepted_panels =
      matvec_stats.convergence_channel_low_abs_reaccepted_panels;
  result.convergence_channel_high_error_caught =
      matvec_stats.convergence_channel_high_error_caught;
  result.convergence_channel_high_error_missed =
      matvec_stats.convergence_channel_high_error_missed;
  result.convergence_channel_high_error_caught_by_split =
      matvec_stats.convergence_channel_high_error_caught_by_split;
  result.convergence_estimator_high_error_caught =
      matvec_stats.convergence_estimator_high_error_caught;
  result.convergence_estimator_high_error_missed =
      matvec_stats.convergence_estimator_high_error_missed;
  result.convergence_estimator_high_error_caught_by_split =
      matvec_stats.convergence_estimator_high_error_caught_by_split;
  result.convergence_estimator_phantom_high_error_count =
      matvec_stats.convergence_estimator_phantom_high_error_count;
  result.convergence_estimator_phantom_high_error_abs_max =
      matvec_stats.convergence_estimator_phantom_high_error_abs_max;
  result.convergence_overreject_count =
      matvec_stats.convergence_overreject_count;
  result.convergence_overreject_panels =
      matvec_stats.convergence_overreject_panels;
  result.convergence_overreject_caught_count =
      matvec_stats.convergence_overreject_caught_count;
  result.convergence_overreject_caught_panels =
      matvec_stats.convergence_overreject_caught_panels;
  result.convergence_overreject_abs_error_max =
      matvec_stats.convergence_overreject_abs_error_max;
  result.convergence_overreject_abs_error_mean =
      matvec_stats.convergence_overreject_abs_error_mean;
  result.convergence_overreject_abs_error_p95 =
      matvec_stats.convergence_overreject_abs_error_p95;
  result.convergence_overreject_abs_error_p99 =
      matvec_stats.convergence_overreject_abs_error_p99;
  result.convergence_overreject_estimated_abs_error_max =
      matvec_stats.convergence_overreject_estimated_abs_error_max;
  result.convergence_overreject_estimated_abs_error_mean =
      matvec_stats.convergence_overreject_estimated_abs_error_mean;
  result.convergence_overreject_required_scale_max =
      matvec_stats.convergence_overreject_required_scale_max;
  result.convergence_overreject_required_scale_p95 =
      matvec_stats.convergence_overreject_required_scale_p95;
  result.convergence_overreject_required_scale_p99 =
      matvec_stats.convergence_overreject_required_scale_p99;
  result.convergence_channel_overreject_count =
      matvec_stats.convergence_channel_overreject_count;
  result.convergence_channel_overreject_caught_count =
      matvec_stats.convergence_channel_overreject_caught_count;
  result.convergence_channel_overreject_panels =
      matvec_stats.convergence_channel_overreject_panels;
  result.convergence_channel_overreject_abs_error_max =
      matvec_stats.convergence_channel_overreject_abs_error_max;
  result.convergence_channel_overreject_abs_error_mean =
      matvec_stats.convergence_channel_overreject_abs_error_mean;
  result.convergence_channel_overreject_required_scale_max =
      matvec_stats.convergence_channel_overreject_required_scale_max;
  result.convergence_channel_overreject_required_scale_p99 =
      matvec_stats.convergence_channel_overreject_required_scale_p99;
  result.phantom_high_error_channel_count =
      matvec_stats.phantom_high_error_channel_count;
  result.phantom_high_error_channel_abs_max =
      matvec_stats.phantom_high_error_channel_abs_max;
  result.accepted_miss_channel_count =
      matvec_stats.accepted_miss_channel_count;
  result.accepted_miss_channel_panels =
      matvec_stats.accepted_miss_channel_panels;
  result.accepted_miss_channel_rho_min =
      matvec_stats.accepted_miss_channel_rho_min;
  result.accepted_miss_channel_rho_mean =
      matvec_stats.accepted_miss_channel_rho_mean;
  result.accepted_miss_channel_rho_max =
      matvec_stats.accepted_miss_channel_rho_max;
  result.accepted_miss_channel_tail_ratio_mean =
      matvec_stats.accepted_miss_channel_tail_ratio_mean;
  result.accepted_miss_channel_tail_ratio_max =
      matvec_stats.accepted_miss_channel_tail_ratio_max;
  result.accepted_miss_channel_estimated_error_mean =
      matvec_stats.accepted_miss_channel_estimated_error_mean;
  result.accepted_miss_channel_estimated_error_max =
      matvec_stats.accepted_miss_channel_estimated_error_max;
  result.accepted_miss_worst_channel =
      matvec_stats.accepted_miss_worst_channel;
  result.convergence_k3_rho_cap_rejects =
      matvec_stats.convergence_k3_rho_cap_rejects;
  result.convergence_k3_rho_cap_reject_panels =
      matvec_stats.convergence_k3_rho_cap_reject_panels;
  result.mac_surface_margin_rejects =
      matvec_stats.mac_surface_margin_rejects;
  result.mac_surface_margin_reject_panels =
      matvec_stats.mac_surface_margin_reject_panels;
  result.mac_surface_margin_applied =
      matvec_stats.mac_surface_margin_applied;
  result.mac_surface_margin_skipped =
      matvec_stats.mac_surface_margin_skipped;
  result.mac_surface_margin_adaptive_rejects =
      matvec_stats.mac_surface_margin_adaptive_rejects;
  result.mac_surface_margin_adaptive_reject_panels =
      matvec_stats.mac_surface_margin_adaptive_reject_panels;
  result.mac_surface_margin_skipped_accepted_far =
      matvec_stats.mac_surface_margin_skipped_accepted_far;
  result.mac_surface_margin_hard_applied =
      matvec_stats.mac_surface_margin_hard_applied;
  result.mac_surface_margin_soft_applied =
      matvec_stats.mac_surface_margin_soft_applied;
  result.mac_surface_margin_hard_rejects =
      matvec_stats.mac_surface_margin_hard_rejects;
  result.mac_surface_margin_hard_reject_panels =
      matvec_stats.mac_surface_margin_hard_reject_panels;
  result.mac_surface_margin_soft_rejects =
      matvec_stats.mac_surface_margin_soft_rejects;
  result.mac_surface_margin_soft_reject_panels =
      matvec_stats.mac_surface_margin_soft_reject_panels;
  result.moment_audit_nodes = matvec_stats.moment_audit_nodes;
  result.moment_translation_rel_error_max =
      matvec_stats.moment_translation_rel_error_max;
  result.moment_translation_rel_error_mean =
      matvec_stats.moment_translation_rel_error_mean;
  result.moment_translation_offset_max =
      matvec_stats.moment_translation_offset_max;
  result.moment_translation_offset_mean =
      matvec_stats.moment_translation_offset_mean;
  result.moment_cancellation_ratio_max =
      matvec_stats.moment_cancellation_ratio_max;
  result.moment_cancellation_ratio_mean =
      matvec_stats.moment_cancellation_ratio_mean;
  result.moment_dynamic_range_max = matvec_stats.moment_dynamic_range_max;
  result.expansion_audit = matvec_stats.expansion_audit;
  result.expansion_audit_channel = matvec_stats.expansion_audit_channel;
  result.far_expansion_audit_pairs =
      matvec_stats.far_expansion_audit_pairs;
  result.far_expansion_audit_panels =
      matvec_stats.far_expansion_audit_panels;
  result.far_expansion_direct_ms =
      matvec_stats.far_expansion_direct_ms;
  result.far_expansion_rel_error_max =
      matvec_stats.far_expansion_rel_error_max;
  result.far_expansion_rel_error_mean =
      matvec_stats.far_expansion_rel_error_mean;
  result.far_expansion_rel_error_p95 =
      matvec_stats.far_expansion_rel_error_p95;
  result.far_expansion_rel_error_p99 =
      matvec_stats.far_expansion_rel_error_p99;
  result.far_expansion_channel_rel_error_max =
      matvec_stats.far_expansion_channel_rel_error_max;
  result.far_expansion_channel_rel_error_mean =
      matvec_stats.far_expansion_channel_rel_error_mean;
  result.far_expansion_channel_rel_error_p99 =
      matvec_stats.far_expansion_channel_rel_error_p99;
  result.far_expansion_worst_channel =
      matvec_stats.far_expansion_worst_channel;
  result.far_expansion_mean_gate_pass =
      matvec_stats.far_expansion_mean_gate_pass;
  result.far_convergence_rho_max = matvec_stats.far_convergence_rho_max;
  result.far_convergence_rho_mean = matvec_stats.far_convergence_rho_mean;
  result.accepted_far_aspect_ratio_max =
      matvec_stats.accepted_far_aspect_ratio_max;
  result.accepted_far_aspect_ratio_mean =
      matvec_stats.accepted_far_aspect_ratio_mean;
  result.leaf_tile_status = matvec_stats.leaf_tile_status;
  result.leaf_panel_pairs_per_second =
      matvec_stats.leaf_panel_pairs_per_second;
  result.matvec_mode = matvec_stats.matvec_mode;
  result.taylor_order =
      mode == "treecode" ? matvec_stats.taylor_order
                         : std::clamp(options.taylor_order, 3, 4);
  result.derivative_mode = mode == "treecode" ? matvec_stats.derivative_mode
                                              : std::string("none");
  result.interaction_build_ms = matvec_stats.interaction_build_ms;
  result.tree_eval_ms = matvec_stats.tree_eval_ms;
  result.tree_moment_ms = matvec_stats.tree_moment_ms;
  result.tree_traversal_ms = matvec_stats.tree_traversal_ms;
  result.near_eval_ms = matvec_stats.near_eval_ms;
  result.far_eval_ms = matvec_stats.far_eval_ms;
  result.leaf_panel_eval_ms = matvec_stats.leaf_panel_eval_ms;
  result.leaf_tile_eval_ms = matvec_stats.leaf_tile_eval_ms;
  result.leaf_residual_eval_ms = matvec_stats.leaf_residual_eval_ms;
  result.smooth_leaf_node_eval_ms = matvec_stats.smooth_leaf_node_eval_ms;

  const double boundary_energy =
      boundary_energy_kcal(atoms, panels, linear.x, kappa) *
      formulation.boundary_energy_scale;
  const double born = born_reference_kcal(atoms, options.eps_in,
                                          options.eps_out, kappa);
  double energy = boundary_energy;
  double selected_scale = formulation.boundary_energy_scale;
  bool used_fallback = false;
  std::string source = "tabi_boundary_integral";
  if (formulation.energy_projection == "born-check") {
    energy = born;
    source = "tabi_born_check";
    selected_scale = 0.0;
  } else if (formulation.formulation == "calibrated-debug" &&
             std::isfinite(boundary_energy) &&
             std::abs(boundary_energy) > 1.0e-12 &&
             std::isfinite(born) && std::abs(born) > 1.0e-12) {
    selected_scale = born / boundary_energy;
    energy = boundary_energy * selected_scale;
    source = "tabi_calibrated_debug";
  } else if (formulation.energy_projection == "charge-potential-debug") {
    source = "tabi_charge_potential_debug";
  } else if (!std::isfinite(energy) || std::abs(energy) <= 1.0e-12 ||
             (std::abs(born) > 1.0e-12 && energy * born < 0.0)) {
    energy = born;
    used_fallback = true;
    source = "tabi_born_fallback";
  }
  result.boundary_energy_kcal_mol = boundary_energy;
  result.reference_energy_kcal_mol = born;
  result.selected_energy_scale = selected_scale;
  result.fallback_used = used_fallback;
  result.reaction_energy_source = source;
  if (!linear.converged && result.convergence_reason.empty())
    result.convergence_reason = "max_iter";
  if (!linear.converged || !std::isfinite(energy)) {
    result.status = "failed";
    if (!std::isfinite(energy))
      result.convergence_reason = "nan_guard";
  } else if (used_fallback) {
    result.status = "fallback";
  } else {
    result.status = "ok";
  }
  result.energy_kcal_mol = energy;
  result.energy_kj_mol = energy * kKcalToKJ;
  result.energy_kt =
      result.energy_kj_mol /
      (kKBKJPerMolK * std::max(1.0, options.temperature_k));
  result.peak_memory_bytes = result.memory_bytes;
  refresh_matvec_cost_summary(result);
  progress.write("solve_done", result, &matvec_stats, result.matvec_calls,
                 result.iterations, 0, result.residual, result.status, true);
  return result;
}

} // namespace TitanPBE::Tabi
