#include "tabi_fmm_m2l_preconditioner.h"

#include "../titan_executor.h"
#include "tabi_near.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace TitanPBE::Tabi {

namespace {
constexpr std::size_t kMinFmmBlockJacobiPanelCap = 8;
constexpr std::size_t kMaxFmmBlockJacobiPanelCap = 128;
constexpr std::size_t kMaxFmmOverlapOwnedPanelCap = 32;
constexpr std::size_t kMaxFmmOverlapHaloPanelCap = 128;
constexpr std::size_t kMinFmmIluNeighborCap = 32;
constexpr std::size_t kMaxFmmIluNeighborCap = 1024;
constexpr double kMinIluPivotDet = 1.0e-18;

std::vector<JcpOwnedFactorPanelGroup> to_owned_factor_groups(
    const std::vector<titan::tabi::fmm_m2l::FmmM2lOverlapPanelGroup>
        &overlap_groups) {
  std::vector<JcpOwnedFactorPanelGroup> converted;
  converted.reserve(overlap_groups.size());
  for (const auto &group : overlap_groups) {
    JcpOwnedFactorPanelGroup out;
    out.owned_panels = group.owned_panels;
    out.factor_panels = group.factor_panels;
    converted.push_back(std::move(out));
  }
  return converted;
}

using Block2 = JcpSparseIluPreconditioner::Block2;

Block2 add(Block2 lhs, const Block2 &rhs) noexcept {
  lhs.a00 += rhs.a00;
  lhs.a01 += rhs.a01;
  lhs.a10 += rhs.a10;
  lhs.a11 += rhs.a11;
  return lhs;
}

Block2 sub(Block2 lhs, const Block2 &rhs) noexcept {
  lhs.a00 -= rhs.a00;
  lhs.a01 -= rhs.a01;
  lhs.a10 -= rhs.a10;
  lhs.a11 -= rhs.a11;
  return lhs;
}

Block2 mul(const Block2 &lhs, const Block2 &rhs) noexcept {
  return {lhs.a00 * rhs.a00 + lhs.a01 * rhs.a10,
          lhs.a00 * rhs.a01 + lhs.a01 * rhs.a11,
          lhs.a10 * rhs.a00 + lhs.a11 * rhs.a10,
          lhs.a10 * rhs.a01 + lhs.a11 * rhs.a11};
}

std::pair<double, double> mul_vec(const Block2 &block, double x0,
                                  double x1) noexcept {
  return {block.a00 * x0 + block.a01 * x1,
          block.a10 * x0 + block.a11 * x1};
}

double max_abs(const Block2 &block) noexcept {
  return std::max({std::abs(block.a00), std::abs(block.a01),
                   std::abs(block.a10), std::abs(block.a11)});
}

bool invert_block2(const Block2 &block, Block2 &inverse) noexcept {
  const double det = block.a00 * block.a11 - block.a01 * block.a10;
  const double scale = std::max(1.0, max_abs(block));
  if (!std::isfinite(det) || std::abs(det) <= kMinIluPivotDet * scale * scale)
    return false;
  const double inv_det = 1.0 / det;
  inverse = {block.a11 * inv_det, -block.a01 * inv_det,
             -block.a10 * inv_det, block.a00 * inv_det};
  return std::isfinite(inverse.a00) && std::isfinite(inverse.a01) &&
         std::isfinite(inverse.a10) && std::isfinite(inverse.a11);
}

std::size_t find_block(const std::vector<std::size_t> &row_ptr,
                       const std::vector<std::size_t> &col_index,
                       std::size_t row, std::size_t col) noexcept {
  if (row + 1 >= row_ptr.size())
    return static_cast<std::size_t>(-1);
  const std::size_t begin = row_ptr[row];
  const std::size_t end = row_ptr[row + 1];
  const auto it = std::lower_bound(col_index.begin() +
                                       static_cast<std::ptrdiff_t>(begin),
                                   col_index.begin() +
                                       static_cast<std::ptrdiff_t>(end),
                                   col);
  if (it == col_index.begin() + static_cast<std::ptrdiff_t>(end) ||
      *it != col)
    return static_cast<std::size_t>(-1);
  return static_cast<std::size_t>(std::distance(col_index.begin(), it));
}

Block2 make_jcp_block(const Panel &target, std::size_t target_index,
                      const Panel &source, std::size_t source_index,
                      double eps_ratio, double kappa, bool coupled,
                      const std::string &near_integration,
                      const std::string &quad_buffer) {
  const JcpNearPairKernel selected = evaluate_jcp_near_pair_kernel(
      target, target_index, source, source_index, eps_ratio, kappa, coupled,
      near_integration, quad_buffer, nullptr);
  if (!selected.valid)
    return {};
  if (coupled) {
    return {-selected.k2, -selected.k1, -selected.k4, -selected.k3};
  }
  return {-selected.k2, 0.0, 0.0, 0.0};
}
}

void JcpSparseIluPreconditioner::apply(const std::vector<double> &in,
                                       std::vector<double> &out) const {
  const std::size_t n = panel_count;
  const std::size_t expected = coupled ? 2 * n : n;
  out.assign(in.size(), 0.0);
  if (n == 0 || in.size() < expected || row_ptr.size() != n + 1 ||
      diag_inverse.size() != n) {
    for (std::size_t i = 0; i < in.size(); ++i)
      out[i] = in[i] * (i < inv_diag.size() ? inv_diag[i] : 1.0);
    return;
  }

  std::vector<double> y0(n, 0.0), y1(n, 0.0), x0(n, 0.0), x1(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    double rhs0 = i < in.size() ? in[i] : 0.0;
    double rhs1 = coupled && n + i < in.size() ? in[n + i] : 0.0;
    for (std::size_t p = row_ptr[i]; p < row_ptr[i + 1]; ++p) {
      const std::size_t col = col_index[p];
      if (col >= i)
        break;
      const auto [v0, v1] = mul_vec(values[p], y0[col], y1[col]);
      rhs0 -= v0;
      rhs1 -= v1;
    }
    y0[i] = rhs0;
    y1[i] = rhs1;
  }

  for (std::size_t rev = 0; rev < n; ++rev) {
    const std::size_t i = n - 1 - rev;
    double rhs0 = y0[i];
    double rhs1 = y1[i];
    for (std::size_t p = row_ptr[i]; p < row_ptr[i + 1]; ++p) {
      const std::size_t col = col_index[p];
      if (col <= i)
        continue;
      const auto [v0, v1] = mul_vec(values[p], x0[col], x1[col]);
      rhs0 -= v0;
      rhs1 -= v1;
    }
    const auto [s0, s1] = mul_vec(diag_inverse[i], rhs0, rhs1);
    x0[i] = s0;
    x1[i] = s1;
  }

  out.assign(expected, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = x0[i];
    if (coupled)
      out[n + i] = x1[i];
  }
}

std::size_t JcpSparseIluPreconditioner::estimated_bytes() const noexcept {
  return inv_diag.size() * sizeof(double) + row_ptr.size() * sizeof(std::size_t) +
         col_index.size() * sizeof(std::size_t) +
         values.size() * sizeof(Block2) + diag_inverse.size() * sizeof(Block2);
}

std::size_t estimate_fmm_m2l_block_jacobi_bytes(
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context,
    std::size_t diagonal_size, bool coupled, std::size_t block_panel_cap) {
  const std::size_t cap = std::clamp<std::size_t>(
      block_panel_cap, kMinFmmBlockJacobiPanelCap,
      kMaxFmmBlockJacobiPanelCap);
  const std::vector<std::vector<std::size_t>> panel_groups =
      context.block_panel_ids(cap);
  return estimate_jcp_block_preconditioner_bytes_from_panel_groups(
      panel_groups, diagonal_size, coupled);
}

std::shared_ptr<JcpBlockPreconditioner>
build_fmm_m2l_block_jacobi_preconditioner(
    const std::vector<Panel> &panels,
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context, bool coupled,
    double eps_ratio, double kappa, double lhs_phi, double lhs_phin,
    const std::vector<double> &diagonal, std::size_t block_panel_cap,
    TitanExecutor &executor, const std::string &near_integration,
    const std::string &quad_buffer) {
  if (panels.empty())
    return {};
  if (context.panel_count() != panels.size()) {
    throw std::invalid_argument(
        "FMM/M2L block-Jacobi preconditioner context does not match panels");
  }
  const std::size_t cap = std::clamp<std::size_t>(
      block_panel_cap, kMinFmmBlockJacobiPanelCap,
      kMaxFmmBlockJacobiPanelCap);
  const std::vector<std::vector<std::size_t>> panel_groups =
      context.block_panel_ids(cap);
  if (panel_groups.empty())
    return {};
  return build_jcp_block_preconditioner_from_panel_groups(
      panels, coupled, eps_ratio, kappa, lhs_phi, lhs_phin, diagonal,
      near_integration, quad_buffer, "operator-near", panel_groups, executor);
}

std::size_t estimate_fmm_m2l_overlap_schwarz_bytes(
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context,
    std::size_t diagonal_size, bool coupled, std::size_t owned_panel_cap,
    std::size_t overlap_panel_cap) {
  const std::size_t owned_cap = std::clamp<std::size_t>(
      owned_panel_cap, kMinFmmBlockJacobiPanelCap,
      kMaxFmmOverlapOwnedPanelCap);
  const std::size_t halo_cap =
      std::clamp<std::size_t>(overlap_panel_cap, 0,
                              kMaxFmmOverlapHaloPanelCap);
  const auto overlap_groups = context.overlap_panel_groups(owned_cap, halo_cap);
  const auto panel_groups = to_owned_factor_groups(overlap_groups);
  return estimate_jcp_block_preconditioner_bytes_from_owned_factor_groups(
      panel_groups, diagonal_size, coupled);
}

std::shared_ptr<JcpBlockPreconditioner>
build_fmm_m2l_overlap_schwarz_preconditioner(
    const std::vector<Panel> &panels,
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context, bool coupled,
    double eps_ratio, double kappa, double lhs_phi, double lhs_phin,
    const std::vector<double> &diagonal, std::size_t owned_panel_cap,
    std::size_t overlap_panel_cap, TitanExecutor &executor,
    const std::string &near_integration, const std::string &quad_buffer) {
  if (panels.empty())
    return {};
  if (context.panel_count() != panels.size()) {
    throw std::invalid_argument(
        "FMM/M2L overlap-Schwarz preconditioner context does not match panels");
  }
  const std::size_t owned_cap = std::clamp<std::size_t>(
      owned_panel_cap, kMinFmmBlockJacobiPanelCap,
      kMaxFmmOverlapOwnedPanelCap);
  const std::size_t halo_cap =
      std::clamp<std::size_t>(overlap_panel_cap, 0,
                              kMaxFmmOverlapHaloPanelCap);
  const auto overlap_groups = context.overlap_panel_groups(owned_cap, halo_cap);
  if (overlap_groups.empty())
    return {};
  const auto panel_groups = to_owned_factor_groups(overlap_groups);
  return build_jcp_block_preconditioner_from_owned_factor_groups(
      panels, coupled, eps_ratio, kappa, lhs_phi, lhs_phin, diagonal,
      near_integration, quad_buffer, "operator-near", panel_groups, executor);
}

std::size_t estimate_fmm_m2l_nearfield_ilu0_bytes(
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context,
    std::size_t diagonal_size, bool, std::size_t neighbor_cap) {
  const std::size_t cap = std::clamp<std::size_t>(
      neighbor_cap, kMinFmmIluNeighborCap, kMaxFmmIluNeighborCap);
  const auto graph = context.near_panel_graph(cap);
  std::size_t nnz = 0;
  for (const auto &row : graph)
    nnz += row.size();
  return diagonal_size * sizeof(double) +
         (graph.size() + 1) * sizeof(std::size_t) +
         nnz * (sizeof(std::size_t) + sizeof(Block2)) +
         graph.size() * sizeof(Block2);
}

std::shared_ptr<JcpSparseIluPreconditioner>
build_fmm_m2l_nearfield_ilu0_preconditioner(
    const std::vector<Panel> &panels,
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context, bool coupled,
    double eps_ratio, double kappa, double lhs_phi, double lhs_phin,
    const std::vector<double> &diagonal, std::size_t neighbor_cap,
    double drop_tol, TitanExecutor &, const std::string &near_integration,
    const std::string &quad_buffer) {
  if (panels.empty())
    return {};
  if (!coupled)
    return {};
  if (context.panel_count() != panels.size()) {
    throw std::invalid_argument(
        "FMM/M2L near-field ILU context does not match panels");
  }
  const std::size_t n = panels.size();
  const std::size_t cap = std::clamp<std::size_t>(
      neighbor_cap, kMinFmmIluNeighborCap, kMaxFmmIluNeighborCap);
  const double effective_drop =
      std::isfinite(drop_tol) && drop_tol > 0.0 ? drop_tol : 0.0;
  auto precond = std::make_shared<JcpSparseIluPreconditioner>();
  precond->panel_count = n;
  precond->coupled = coupled;
  precond->inv_diag.assign(diagonal.size(), 1.0);
  for (std::size_t i = 0; i < diagonal.size(); ++i) {
    const double diag = std::abs(diagonal[i]);
    precond->inv_diag[i] =
        std::isfinite(diag) && diag > 1.0e-12 ? 1.0 / diag : 1.0;
  }

  auto graph = context.near_panel_graph(cap);
  if (graph.size() != n)
    return {};
  precond->row_ptr.assign(n + 1, 0);
  std::size_t nnz = 0;
  std::size_t max_degree = 0;
  for (std::size_t i = 0; i < n; ++i) {
    auto &row = graph[i];
    row.push_back(i);
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
    if (row.size() > cap) {
      row.resize(cap);
      if (!std::binary_search(row.begin(), row.end(), i)) {
        row.back() = i;
        std::sort(row.begin(), row.end());
      }
    }
    precond->row_ptr[i] = nnz;
    nnz += row.size();
    max_degree = std::max(max_degree, row.size());
  }
  precond->row_ptr[n] = nnz;
  precond->col_index.reserve(nnz);
  precond->values.reserve(nnz);
  for (std::size_t i = 0; i < n; ++i) {
    for (const std::size_t col : graph[i]) {
      precond->col_index.push_back(col);
      Block2 block{};
      if (i == col) {
        block.a00 += lhs_phi;
        block.a11 += lhs_phin;
      }
      block = add(block,
                  make_jcp_block(panels[i], i, panels[col], col, eps_ratio,
                                 kappa, coupled, near_integration,
                                 quad_buffer));
      if (i != col && effective_drop > 0.0 && max_abs(block) < effective_drop)
        block = {};
      precond->values.push_back(block);
    }
  }
  precond->nnz_blocks = precond->values.size();
  precond->valid_blocks = precond->nnz_blocks;
  precond->valid_panel_count = n;
  precond->max_degree = max_degree;
  precond->avg_degree =
      n > 0 ? static_cast<double>(precond->nnz_blocks) / static_cast<double>(n)
            : 0.0;
  precond->diag_inverse.assign(n, {});

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t p = precond->row_ptr[i]; p < precond->row_ptr[i + 1];
         ++p) {
      const std::size_t k = precond->col_index[p];
      if (k >= i)
        break;
      const Block2 diag_inv_k = precond->diag_inverse[k];
      precond->values[p] = mul(precond->values[p], diag_inv_k);
      const Block2 lik = precond->values[p];
      for (std::size_t q = precond->row_ptr[k]; q < precond->row_ptr[k + 1];
           ++q) {
        const std::size_t j = precond->col_index[q];
        if (j <= k)
          continue;
        const std::size_t ij = find_block(precond->row_ptr,
                                          precond->col_index, i, j);
        if (ij == static_cast<std::size_t>(-1))
          continue;
        precond->values[ij] = sub(precond->values[ij],
                                  mul(lik, precond->values[q]));
      }
    }
    const std::size_t diag_pos =
        find_block(precond->row_ptr, precond->col_index, i, i);
    if (diag_pos == static_cast<std::size_t>(-1) ||
        !invert_block2(precond->values[diag_pos], precond->diag_inverse[i])) {
      ++precond->ill_conditioned_blocks;
      precond->quality_status = "factor_fallback";
      return precond;
    }
  }
  precond->quality_status = "pass";
  return precond;
}

} // namespace TitanPBE::Tabi
