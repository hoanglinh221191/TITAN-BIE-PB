#ifndef TITAN_TABI_FMM_M2L_FALLBACK_H
#define TITAN_TABI_FMM_M2L_FALLBACK_H

#include "tabi_types.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class TitanExecutor;

namespace titan::tabi::fmm_m2l {

struct PackedFallbackSourceSpan {
  std::size_t first_panel{0};
  std::size_t panel_count{0};
  bool refined_fallback{false};
};

struct PackedFallbackLeafInput {
  std::vector<std::size_t> target_panels;
  std::vector<PackedFallbackSourceSpan> source_spans;
};

struct PackedFallbackPlanMetrics {
  std::uint64_t worklist_hash{0};
  std::size_t build_count{0};
  std::size_t target_batch_size{4};
  std::size_t source_tile_size{32};
  std::size_t target_batches{0};
  std::size_t task_count{0};
  std::size_t metadata_bytes{0};
  std::size_t total_pairs{0};
  std::size_t base_near_pairs{0};
  std::size_t refined_fallback_pairs{0};
  std::size_t self_pairs{0};
};

struct PackedFallbackApplyResult {
  std::vector<double> phi;
  std::vector<double> phin;
  std::string effective_kernel{"scalar"};
  std::size_t vectorized_sources{0};
  std::size_t vectorized_lanes{0};
  std::size_t scalar_self_pairs{0};
  double apply_seconds{0.0};
  bool finite_outputs{true};
};

class PackedJcpFallbackPlan {
public:
  PackedJcpFallbackPlan();
  ~PackedJcpFallbackPlan();
  PackedJcpFallbackPlan(PackedJcpFallbackPlan &&) noexcept;
  PackedJcpFallbackPlan &operator=(PackedJcpFallbackPlan &&) noexcept;
  PackedJcpFallbackPlan(const PackedJcpFallbackPlan &) = delete;
  PackedJcpFallbackPlan &operator=(const PackedJcpFallbackPlan &) = delete;

  [[nodiscard]] PackedFallbackPlanMetrics metrics() const noexcept;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  friend PackedJcpFallbackPlan build_packed_jcp_fallback_plan(
      const std::vector<TitanPBE::Tabi::Panel> &,
      const std::vector<std::size_t> &,
      const std::vector<PackedFallbackLeafInput> &, std::size_t,
      std::size_t);
  friend PackedFallbackApplyResult apply_packed_jcp_fallback(
      PackedJcpFallbackPlan &, const std::vector<TitanPBE::Tabi::Panel> &,
      const std::vector<double> &, const std::vector<double> &, double,
      double, bool, const std::string &, const std::string &,
      const std::string &, TitanExecutor &);
};

PackedJcpFallbackPlan build_packed_jcp_fallback_plan(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<std::size_t> &source_order,
    const std::vector<PackedFallbackLeafInput> &leaves,
    std::size_t target_batch_size = 4, std::size_t source_tile_size = 32);

PackedFallbackApplyResult apply_packed_jcp_fallback(
    PackedJcpFallbackPlan &plan,
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &phi, const std::vector<double> &phin,
    double eps_ratio, double kappa, bool coupled,
    const std::string &near_integration, const std::string &quad_buffer,
    const std::string &kernel, TitanExecutor &executor);

} // namespace titan::tabi::fmm_m2l

#endif
