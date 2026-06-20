#include "tabi_fmm_m2l_fallback.h"

#include "tabi_near.h"
#include "tabi_treecode_internal.h"
#include "titan_executor.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace titan::tabi::fmm_m2l {
namespace {
using Clock = std::chrono::steady_clock;
using TitanPBE::Tabi::JcpNearPairKernel;
using TitanPBE::Tabi::Panel;
using TitanPBE::Tabi::detail::JcpAccum;

struct PackedTask {
  std::size_t source_first{0};
  std::size_t source_count{0};
  bool refined_fallback{false};
};

struct PackedBatch {
  std::vector<std::size_t> targets;
  std::size_t task_begin{0};
  std::size_t task_count{0};
  std::vector<std::size_t> self_targets;
};

void hash_mix(std::uint64_t &hash, std::uint64_t value) noexcept {
  hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
}

void add_kernel(JcpAccum &acc, const JcpNearPairKernel &kernel,
                double phi, double phin, bool coupled) noexcept {
  acc.phi += kernel.k2 * phi;
  if (!coupled)
    return;
  acc.phi += kernel.k1 * phin;
  acc.phin += kernel.k3 * phin + kernel.k4 * phi;
}
} // namespace

struct PackedJcpFallbackPlan::Impl {
  std::size_t panel_count{0};
  std::vector<std::size_t> source_order;
  TitanPBE::Tabi::TreePanelCache panel_cache;
  TitanPBE::Tabi::detail::TreeUnknownCache unknown_cache;
  std::vector<PackedBatch> batches;
  std::vector<PackedTask> tasks;
  PackedFallbackPlanMetrics metrics;
};

PackedJcpFallbackPlan::PackedJcpFallbackPlan()
    : impl_(std::make_unique<Impl>()) {}
PackedJcpFallbackPlan::~PackedJcpFallbackPlan() = default;
PackedJcpFallbackPlan::PackedJcpFallbackPlan(PackedJcpFallbackPlan &&) noexcept =
    default;
PackedJcpFallbackPlan &PackedJcpFallbackPlan::operator=(
    PackedJcpFallbackPlan &&) noexcept = default;

PackedFallbackPlanMetrics PackedJcpFallbackPlan::metrics() const noexcept {
  return impl_ ? impl_->metrics : PackedFallbackPlanMetrics{};
}

PackedJcpFallbackPlan build_packed_jcp_fallback_plan(
    const std::vector<Panel> &panels,
    const std::vector<std::size_t> &source_order,
    const std::vector<PackedFallbackLeafInput> &leaves,
    std::size_t target_batch_size, std::size_t source_tile_size) {
  if (panels.empty() || source_order.size() != panels.size())
    throw std::invalid_argument(
        "packed fallback requires a complete source order");
  if (target_batch_size != 4 || source_tile_size != 32)
    throw std::invalid_argument(
        "packed fallback requires target batch 4 and source tile 32");
  std::vector<bool> seen_order(panels.size(), false);
  for (const std::size_t panel : source_order) {
    if (panel >= panels.size() || seen_order[panel])
      throw std::invalid_argument(
          "packed fallback source order must be a panel permutation");
    seen_order[panel] = true;
  }

  PackedJcpFallbackPlan plan;
  auto &impl = *plan.impl_;
  impl.panel_count = panels.size();
  impl.source_order = source_order;
  impl.panel_cache = TitanPBE::Tabi::detail::build_tree_panel_cache(
      panels, source_order);
  impl.unknown_cache.phi.assign(panels.size(), 0.0);
  impl.unknown_cache.phin.assign(panels.size(), 0.0);
  impl.unknown_cache.ordered_phi.assign(panels.size(), 0.0);
  impl.unknown_cache.ordered_phin.assign(panels.size(), 0.0);
  impl.metrics.build_count = 1;
  impl.metrics.target_batch_size = target_batch_size;
  impl.metrics.source_tile_size = source_tile_size;

  for (const PackedFallbackLeafInput &leaf : leaves) {
    std::unordered_set<std::size_t> covered_sources;
    for (const PackedFallbackSourceSpan &span : leaf.source_spans) {
      if (span.first_panel + span.panel_count > source_order.size())
        throw std::invalid_argument("packed fallback source span is invalid");
      for (std::size_t pos = span.first_panel;
           pos < span.first_panel + span.panel_count; ++pos) {
        if (!covered_sources.insert(source_order[pos]).second)
          throw std::invalid_argument(
              "packed fallback source spans overlap within a target leaf");
      }
      const std::size_t pairs = leaf.target_panels.size() * span.panel_count;
      impl.metrics.total_pairs += pairs;
      if (span.refined_fallback)
        impl.metrics.refined_fallback_pairs += pairs;
      else
        impl.metrics.base_near_pairs += pairs;
    }
    for (std::size_t target_begin = 0;
         target_begin < leaf.target_panels.size();
         target_begin += target_batch_size) {
      PackedBatch batch;
      const std::size_t target_end = std::min(
          leaf.target_panels.size(), target_begin + target_batch_size);
      batch.targets.assign(leaf.target_panels.begin() + target_begin,
                           leaf.target_panels.begin() + target_end);
      for (const std::size_t target : batch.targets) {
        if (target >= panels.size())
          throw std::invalid_argument("packed fallback target is invalid");
        if (covered_sources.contains(target))
          batch.self_targets.push_back(target);
      }
      impl.metrics.self_pairs += batch.self_targets.size();
      batch.task_begin = impl.tasks.size();
      for (const PackedFallbackSourceSpan &span : leaf.source_spans) {
        for (std::size_t offset = 0; offset < span.panel_count;
             offset += source_tile_size) {
          impl.tasks.push_back(
              {span.first_panel + offset,
               std::min(source_tile_size, span.panel_count - offset),
               span.refined_fallback});
        }
      }
      batch.task_count = impl.tasks.size() - batch.task_begin;
      impl.batches.push_back(std::move(batch));
    }
  }
  impl.metrics.target_batches = impl.batches.size();
  impl.metrics.task_count = impl.tasks.size();
  impl.metrics.metadata_bytes =
      impl.tasks.size() * sizeof(PackedTask) +
      impl.batches.size() * sizeof(PackedBatch);
  for (const PackedBatch &batch : impl.batches)
    impl.metrics.metadata_bytes +=
        batch.targets.capacity() * sizeof(std::size_t) +
        batch.self_targets.capacity() * sizeof(std::size_t);

  std::uint64_t hash = 0xcbf29ce484222325ULL;
  hash_mix(hash, target_batch_size);
  hash_mix(hash, source_tile_size);
  for (const PackedBatch &batch : impl.batches) {
    for (const std::size_t target : batch.targets)
      hash_mix(hash, target);
    for (std::size_t i = 0; i < batch.task_count; ++i) {
      const PackedTask &task = impl.tasks[batch.task_begin + i];
      hash_mix(hash, task.source_first);
      hash_mix(hash, task.source_count);
      hash_mix(hash, task.refined_fallback ? 1 : 0);
    }
  }
  impl.metrics.worklist_hash = hash;
  return plan;
}

PackedFallbackApplyResult apply_packed_jcp_fallback(
    PackedJcpFallbackPlan &plan, const std::vector<Panel> &panels,
    const std::vector<double> &phi, const std::vector<double> &phin,
    double eps_ratio, double kappa, bool coupled,
    const std::string &near_integration, const std::string &quad_buffer,
    const std::string &kernel, TitanExecutor &executor) {
  if (!plan.impl_ || panels.size() != plan.impl_->panel_count ||
      phi.size() != panels.size() || phin.size() != panels.size())
    throw std::invalid_argument("packed fallback apply shape mismatch");
  if (kernel != "scalar" && kernel != "packed-avx2")
    throw std::invalid_argument(
        "packed fallback kernel must be scalar or packed-avx2");
  if (kernel == "packed-avx2" &&
      (near_integration != "analytic_singular" || quad_buffer != "off"))
    throw std::invalid_argument(
        "packed-avx2 requires analytic_singular and quad-buffer=off");
  if (kernel == "packed-avx2" &&
      !TitanPBE::Tabi::detail::packed_jcp_avx2_available())
    throw std::runtime_error("packed_avx2_unavailable");

  auto &impl = *plan.impl_;
  impl.unknown_cache.phi = phi;
  impl.unknown_cache.phin = phin;
  for (std::size_t pos = 0; pos < impl.source_order.size(); ++pos) {
    const std::size_t source = impl.source_order[pos];
    impl.unknown_cache.ordered_phi[pos] = phi[source];
    impl.unknown_cache.ordered_phin[pos] = phin[source];
  }

  PackedFallbackApplyResult result;
  result.phi.assign(panels.size(), 0.0);
  result.phin.assign(panels.size(), 0.0);
  result.effective_kernel = kernel;
  const auto begin_time = Clock::now();
  std::vector<std::size_t> vector_sources(impl.batches.size(), 0);
  std::vector<std::size_t> vector_lanes(impl.batches.size(), 0);
  std::vector<std::size_t> self_counts(impl.batches.size(), 0);
  const std::size_t workers = std::max<std::size_t>(1, executor.capacity());
  executor.parallel_for(
      impl.batches.size(), workers, 1,
      [&](std::size_t, std::size_t batch_begin, std::size_t batch_end) {
        for (std::size_t batch_index = batch_begin;
             batch_index < batch_end; ++batch_index) {
          const PackedBatch &batch = impl.batches[batch_index];
          std::vector<JcpAccum> total(batch.targets.size());
          if (kernel == "packed-avx2") {
            for (std::size_t task_offset = 0;
                 task_offset < batch.task_count; ++task_offset) {
              const PackedTask &task =
                  impl.tasks[batch.task_begin + task_offset];
              std::vector<JcpAccum> task_accum;
              std::size_t near_count = 0;
              TitanPBE::Tabi::detail::direct_jcp_packed_target_source_span_avx2(
                  panels, batch.targets, impl.panel_cache,
                  impl.unknown_cache, eps_ratio, kappa, coupled,
                  task.source_first, task.source_count, task_accum,
                  near_count, vector_sources[batch_index],
                  vector_lanes[batch_index]);
              for (std::size_t lane = 0; lane < total.size(); ++lane) {
                total[lane].phi += task_accum[lane].phi;
                total[lane].phin += task_accum[lane].phin;
              }
            }
            for (const std::size_t target : batch.self_targets) {
              const JcpNearPairKernel self =
                  TitanPBE::Tabi::evaluate_jcp_near_pair_kernel(
                      panels[target], target, panels[target], target,
                      eps_ratio, kappa, coupled, near_integration,
                      quad_buffer, nullptr);
              const auto position = std::find(batch.targets.begin(),
                                              batch.targets.end(), target);
              const std::size_t lane = static_cast<std::size_t>(
                  std::distance(batch.targets.begin(), position));
              add_kernel(total[lane], self, phi[target], phin[target], coupled);
              ++self_counts[batch_index];
            }
          } else {
            for (std::size_t lane = 0; lane < batch.targets.size(); ++lane) {
              const std::size_t target = batch.targets[lane];
              for (std::size_t task_offset = 0;
                   task_offset < batch.task_count; ++task_offset) {
                const PackedTask &task =
                    impl.tasks[batch.task_begin + task_offset];
                for (std::size_t pos = task.source_first;
                     pos < task.source_first + task.source_count; ++pos) {
                  const std::size_t source = impl.source_order[pos];
                  const JcpNearPairKernel selected =
                      TitanPBE::Tabi::evaluate_jcp_near_pair_kernel(
                          panels[target], target, panels[source], source,
                          eps_ratio, kappa, coupled, near_integration,
                          quad_buffer, nullptr);
                  add_kernel(total[lane], selected, phi[source], phin[source],
                             coupled);
                  if (source == target)
                    ++self_counts[batch_index];
                }
              }
            }
          }
          for (std::size_t lane = 0; lane < batch.targets.size(); ++lane) {
            const std::size_t target = batch.targets[lane];
            result.phi[target] = total[lane].phi;
            if (coupled)
              result.phin[target] = total[lane].phin;
          }
        }
      });
  result.vectorized_sources =
      std::accumulate(vector_sources.begin(), vector_sources.end(),
                      std::size_t{0});
  result.vectorized_lanes =
      std::accumulate(vector_lanes.begin(), vector_lanes.end(),
                      std::size_t{0});
  result.scalar_self_pairs =
      std::accumulate(self_counts.begin(), self_counts.end(), std::size_t{0});
  result.apply_seconds =
      std::chrono::duration<double>(Clock::now() - begin_time).count();
  result.finite_outputs =
      std::all_of(result.phi.begin(), result.phi.end(),
                  [](double value) { return std::isfinite(value); }) &&
      std::all_of(result.phin.begin(), result.phin.end(),
                  [](double value) { return std::isfinite(value); });
  return result;
}

} // namespace titan::tabi::fmm_m2l
