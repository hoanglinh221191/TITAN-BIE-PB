#include "tabi_preconditioner.h"

#include "tabi_near.h"

#include "../titan_executor.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace TitanPBE::Tabi {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

std::size_t solver_workers(const TitanExecutor &executor,
                           std::size_t work_items) noexcept {
  if (work_items == 0)
    return 1;
  return std::max<std::size_t>(
      1, std::min<std::size_t>(executor.capacity(), work_items));
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

double mixed_normal_derivative_green(const Panel &target, const Panel &source,
                                     double kappa) noexcept {
  const double ux = target.x - source.x;
  const double uy = target.y - source.y;
  const double uz = target.z - source.z;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  const double inv_r = 1.0 / r;
  const double nxu =
      (target.nx * ux + target.ny * uy + target.nz * uz) * inv_r;
  const double nyu =
      (source.nx * ux + source.ny * uy + source.nz * uz) * inv_r;
  const double nxny = target.nx * source.nx + target.ny * source.ny +
                      target.nz * source.nz;
  const double expkr = std::exp(-kappa * r);
  const double inv_r2 = inv_r * inv_r;
  const double fprime =
      -expkr * (1.0 + kappa * r) * inv_r2 / (4.0 * kPi);
  const double fsecond =
      expkr * (kappa * kappa * inv_r + 2.0 * kappa * inv_r2 +
               2.0 * inv_r2 * inv_r) /
      (4.0 * kPi);
  return fsecond * nxu * nyu + (fprime * inv_r) * (nxny - nxu * nyu);
}

struct JcpKernels {
  double k1{0.0};
  double k2{0.0};
  double k3{0.0};
  double k4{0.0};
};

JcpKernels jcp_surface_kernels(const Panel &target, const Panel &source,
                               double eps_ratio, double kappa) noexcept {
  const double dx = target.x - source.x;
  const double dy = target.y - source.y;
  const double dz = target.z - source.z;
  const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (r <= 1.0e-14)
    return {};
  const double g0 = green_kernel(r, 0.0);
  const double gk = green_kernel(r, kappa);
  const double d0_src =
      source_normal_derivative_green(target.x, target.y, target.z, source, 0.0);
  const double dk_src = source_normal_derivative_green(
      target.x, target.y, target.z, source, kappa);
  const double d0_tgt = target_normal_derivative_green(
      target, source.x, source.y, source.z, 0.0);
  const double dk_tgt = target_normal_derivative_green(
      target, source.x, source.y, source.z, kappa);
  const double d2_0 = mixed_normal_derivative_green(target, source, 0.0);
  const double d2_k = mixed_normal_derivative_green(target, source, kappa);
  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  return {g0 - gk, eps_ratio * dk_src - d0_src,
          d0_tgt - inv_eps_ratio * dk_tgt, d2_k - d2_0};
}

int shared_vertex_count(const Panel &a, const Panel &b) noexcept {
  if (!a.has_triangle_geometry || !b.has_triangle_geometry)
    return 0;
  auto has_vertex = [](const Panel &panel, std::size_t vertex) noexcept {
    return vertex != static_cast<std::size_t>(-1) &&
           (panel.v0_index == vertex || panel.v1_index == vertex ||
            panel.v2_index == vertex);
  };
  int count = 0;
  if (has_vertex(b, a.v0_index))
    ++count;
  if (a.v1_index != a.v0_index && has_vertex(b, a.v1_index))
    ++count;
  if (a.v2_index != a.v0_index && a.v2_index != a.v1_index &&
      has_vertex(b, a.v2_index))
    ++count;
  return count;
}

bool is_singular_near_pair(const Panel &target, std::size_t target_index,
                           const Panel &source,
                           std::size_t source_index) noexcept {
  if (source_index == target_index)
    return true;
  return shared_vertex_count(target, source) > 0;
}

bool factor_dense_lu(std::vector<double> &matrix, std::vector<int> &pivots,
                     std::size_t dim, double *pivot_ratio = nullptr) {
  pivots.assign(dim, 0);
  double min_pivot = std::numeric_limits<double>::infinity();
  double max_pivot = 0.0;
  for (std::size_t k = 0; k < dim; ++k) {
    std::size_t pivot = k;
    double best = std::abs(matrix[k * dim + k]);
    for (std::size_t row = k + 1; row < dim; ++row) {
      const double value = std::abs(matrix[row * dim + k]);
      if (value > best) {
        best = value;
        pivot = row;
      }
    }
    if (!std::isfinite(best) || best <= 1.0e-14)
      return false;
    min_pivot = std::min(min_pivot, best);
    max_pivot = std::max(max_pivot, best);
    pivots[k] = static_cast<int>(pivot);
    if (pivot != k) {
      for (std::size_t col = 0; col < dim; ++col)
        std::swap(matrix[k * dim + col], matrix[pivot * dim + col]);
    }
    const double diag = matrix[k * dim + k];
    for (std::size_t row = k + 1; row < dim; ++row) {
      matrix[row * dim + k] /= diag;
      const double multiplier = matrix[row * dim + k];
      for (std::size_t col = k + 1; col < dim; ++col)
        matrix[row * dim + col] -= multiplier * matrix[k * dim + col];
    }
  }
  if (pivot_ratio != nullptr)
    *pivot_ratio = max_pivot > 0.0 && std::isfinite(min_pivot)
                       ? min_pivot / max_pivot
                       : 0.0;
  return true;
}

std::vector<double> make_equilibration_scale(const std::vector<double> &matrix,
                                             std::size_t dim) {
  std::vector<double> scale(dim, 1.0);
  for (std::size_t i = 0; i < dim; ++i) {
    const double diag = std::abs(matrix[i * dim + i]);
    if (std::isfinite(diag) && diag > 1.0e-24)
      scale[i] = 1.0 / std::sqrt(diag);
  }
  return scale;
}

void apply_symmetric_scale(std::vector<double> &matrix,
                           const std::vector<double> &scale,
                           std::size_t dim) {
  if (scale.size() != dim)
    return;
  for (std::size_t row = 0; row < dim; ++row) {
    const double sr = scale[row];
    for (std::size_t col = 0; col < dim; ++col)
      matrix[row * dim + col] *= sr * scale[col];
  }
}

void solve_dense_lu(const std::vector<double> &lu,
                    const std::vector<int> &pivots, std::size_t dim,
                    std::vector<double> &rhs) {
  for (std::size_t k = 0; k < dim; ++k) {
    const std::size_t pivot = static_cast<std::size_t>(pivots[k]);
    if (pivot != k)
      std::swap(rhs[k], rhs[pivot]);
  }
  for (std::size_t row = 0; row < dim; ++row) {
    for (std::size_t col = 0; col < row; ++col)
      rhs[row] -= lu[row * dim + col] * rhs[col];
  }
  for (std::size_t rev = 0; rev < dim; ++rev) {
    const std::size_t row = dim - 1 - rev;
    for (std::size_t col = row + 1; col < dim; ++col)
      rhs[row] -= lu[row * dim + col] * rhs[col];
    const double diag = lu[row * dim + row];
    rhs[row] = std::abs(diag) > 1.0e-14 ? rhs[row] / diag : 0.0;
  }
}

double norm2(const std::vector<double> &values) noexcept {
  double sum = 0.0;
  for (const double value : values)
    sum += value * value;
  return std::sqrt(sum);
}

void dense_matvec(const std::vector<double> &matrix,
                  const std::vector<double> &x, std::size_t dim,
                  std::vector<double> &out) {
  out.assign(dim, 0.0);
  for (std::size_t row = 0; row < dim; ++row) {
    double sum = 0.0;
    for (std::size_t col = 0; col < dim; ++col)
      sum += matrix[row * dim + col] * x[col];
    out[row] = sum;
  }
}

void solve_scaled_lu(const std::vector<double> &lu,
                     const std::vector<int> &pivots,
                     const std::vector<double> &scale, std::size_t dim,
                     const std::vector<double> &rhs,
                     std::vector<double> &solution) {
  solution.assign(dim, 0.0);
  if (lu.size() != dim * dim || pivots.size() != dim)
    return;
  solution = rhs;
  if (scale.size() == dim) {
    for (std::size_t i = 0; i < dim; ++i)
      solution[i] *= scale[i];
  }
  solve_dense_lu(lu, pivots, dim, solution);
  if (scale.size() == dim) {
    for (std::size_t i = 0; i < dim; ++i)
      solution[i] *= scale[i];
  }
}

bool passes_block_quality_gate(const std::vector<double> &original_matrix,
                               const std::vector<double> &lu,
                               const std::vector<int> &pivots,
                               const std::vector<double> &scale,
                               std::size_t dim, double *max_gain) {
  if (max_gain != nullptr)
    *max_gain = 0.0;
  if (dim == 0 || original_matrix.size() != dim * dim)
    return false;
  std::vector<std::vector<double>> probes;
  probes.push_back(std::vector<double>(dim, 1.0));
  std::vector<double> alternating(dim, 1.0);
  for (std::size_t i = 1; i < dim; i += 2)
    alternating[i] = -1.0;
  probes.push_back(std::move(alternating));
  std::vector<double> solution;
  std::vector<double> applied;
  for (const std::vector<double> &probe : probes) {
    solve_scaled_lu(lu, pivots, scale, dim, probe, solution);
    for (const double value : solution) {
      if (!std::isfinite(value))
        return false;
    }
    const double gain = norm2(solution) / std::max(1.0e-300, norm2(probe));
    if (max_gain != nullptr)
      *max_gain = std::max(*max_gain, gain);
    if (!std::isfinite(gain) || gain > 1.0e8)
      return false;
    dense_matvec(original_matrix, solution, dim, applied);
    for (std::size_t i = 0; i < dim; ++i)
      applied[i] -= probe[i];
    const double rel =
        norm2(applied) / std::max(1.0e-300, norm2(probe));
    if (!std::isfinite(rel) || rel > 1.0e-6)
      return false;
  }
  return true;
}

std::vector<JcpCoarseGalerkinPreconditioner::Aggregate>
make_jcp_coarse_aggregates(const std::vector<Panel> &panels,
                           const TreeWorkspace &tree,
                           std::size_t aggregate_cap) {
  std::vector<std::size_t> selected;
  if (tree.root_index < tree.nodes.size())
    selected.push_back(tree.root_index);
  if (selected.empty() || tree.order.size() != panels.size()) {
    JcpCoarseGalerkinPreconditioner::Aggregate aggregate;
    aggregate.panels.resize(panels.size());
    std::iota(aggregate.panels.begin(), aggregate.panels.end(), 0);
    for (const Panel &panel : panels)
      aggregate.area_sum += panel.area;
    return {std::move(aggregate)};
  }
  const std::size_t cap = std::clamp<std::size_t>(aggregate_cap, 4, 64);
  while (selected.size() < cap) {
    std::size_t split_pos = selected.size();
    std::size_t split_panels = 0;
    for (std::size_t pos = 0; pos < selected.size(); ++pos) {
      const TreeNode &node = tree.nodes[selected[pos]];
      if (node.child_count < 2 ||
          selected.size() - 1 + node.child_count > cap)
        continue;
      if (node.panel_count > split_panels) {
        split_panels = node.panel_count;
        split_pos = pos;
      }
    }
    if (split_pos >= selected.size())
      break;
    const TreeNode node = tree.nodes[selected[split_pos]];
    selected.erase(selected.begin() + static_cast<std::ptrdiff_t>(split_pos));
    for (std::size_t child = 0; child < node.child_count; ++child) {
      const std::size_t child_pos = node.first_child + child;
      if (child_pos < tree.child_indices.size() &&
          tree.child_indices[child_pos] < tree.nodes.size()) {
        selected.insert(selected.begin() + static_cast<std::ptrdiff_t>(split_pos),
                        tree.child_indices[child_pos]);
        ++split_pos;
      }
    }
  }
  std::vector<JcpCoarseGalerkinPreconditioner::Aggregate> aggregates;
  aggregates.reserve(selected.size());
  for (const std::size_t node_index : selected) {
    const TreeNode &node = tree.nodes[node_index];
    JcpCoarseGalerkinPreconditioner::Aggregate aggregate;
    aggregate.panels.reserve(node.panel_count);
    const std::size_t end =
        std::min(tree.order.size(), node.first_panel + node.panel_count);
    for (std::size_t pos = node.first_panel; pos < end; ++pos) {
      const std::size_t panel = tree.order[pos];
      if (panel < panels.size()) {
        aggregate.panels.push_back(panel);
        aggregate.area_sum += panels[panel].area;
      }
    }
    if (!aggregate.panels.empty())
      aggregates.push_back(std::move(aggregate));
  }
  return aggregates;
}

} // namespace

void JcpBlockPreconditioner::apply_block_with_scratch(
    const Block &block, const std::vector<double> &in,
    std::vector<double> &out, std::vector<double> &scratch) const {
  const std::size_t n = panel_count;
  const std::size_t m = block.panels.size();
  const std::size_t dim = coupled ? 2 * m : m;
  if (block.valid && block.lu.size() == dim * dim) {
    scratch.assign(dim, 0.0);
    for (std::size_t r = 0; r < m; ++r) {
      const std::size_t panel = block.panels[r];
      scratch[r] = panel < in.size() ? in[panel] : 0.0;
      if (coupled && n + panel < in.size())
        scratch[m + r] = in[n + panel];
    }
    if (block.scale.size() == dim) {
      for (std::size_t i = 0; i < dim; ++i)
        scratch[i] *= block.scale[i];
    }
    solve_dense_lu(block.lu, block.pivots, dim, scratch);
    if (block.scale.size() == dim) {
      for (std::size_t i = 0; i < dim; ++i)
        scratch[i] *= block.scale[i];
    }
    for (std::size_t r = 0; r < m; ++r) {
      const std::size_t panel = block.panels[r];
      if (panel < out.size())
        out[panel] = scratch[r];
      if (coupled && n + panel < out.size())
        out[n + panel] = scratch[m + r];
    }
  } else {
    for (const std::size_t panel : block.panels) {
      if (panel < out.size())
        out[panel] =
            in[panel] * (panel < inv_diag.size() ? inv_diag[panel] : 1.0);
      if (coupled && n + panel < out.size())
        out[n + panel] =
            in[n + panel] *
            (n + panel < inv_diag.size() ? inv_diag[n + panel] : 1.0);
    }
  }
}

void JcpBlockPreconditioner::apply_block(const Block &block,
                                         const std::vector<double> &in,
                                         std::vector<double> &out) const {
  std::vector<double> scratch;
  apply_block_with_scratch(block, in, out, scratch);
}

void JcpBlockPreconditioner::apply(const std::vector<double> &in,
                                   std::vector<double> &out) const {
  out.assign(in.size(), 0.0);
  std::vector<double> scratch;
  for (const Block &block : blocks)
    apply_block_with_scratch(block, in, out, scratch);
}

void JcpBlockPreconditioner::apply_parallel(const std::vector<double> &in,
                                            std::vector<double> &out,
                                            TitanExecutor &executor) const {
  out.assign(in.size(), 0.0);
  if (blocks.size() < 4 || executor.capacity() <= 1) {
    apply(in, out);
    return;
  }
  const std::size_t workers = solver_workers(executor, blocks.size());
  executor.parallel_for(
      blocks.size(), workers, 1,
      [&](std::size_t, std::size_t begin, std::size_t end) {
        std::vector<double> scratch;
        for (std::size_t i = begin; i < end; ++i)
          apply_block_with_scratch(blocks[i], in, out, scratch);
      });
}

std::size_t JcpBlockPreconditioner::estimated_bytes() const noexcept {
  std::size_t bytes = inv_diag.size() * sizeof(double);
  for (const Block &block : blocks) {
    bytes += block.panels.size() * sizeof(std::size_t);
    bytes += block.lu.size() * sizeof(double);
    bytes += block.pivots.size() * sizeof(int);
  }
  return bytes;
}

std::size_t JcpCoarseGalerkinPreconditioner::dimension() const noexcept {
  return coupled ? 2 * aggregates.size() : aggregates.size();
}

void JcpCoarseGalerkinPreconditioner::restrict_area_weighted(
    const std::vector<Panel> &panels, const std::vector<double> &fine,
    std::vector<double> &coarse) const {
  const std::size_t q = aggregates.size();
  coarse.assign(dimension(), 0.0);
  for (std::size_t a = 0; a < q; ++a) {
    const Aggregate &aggregate = aggregates[a];
    const double denom = std::max(1.0e-300, aggregate.area_sum);
    double phi = 0.0;
    double phin = 0.0;
    for (const std::size_t panel : aggregate.panels) {
      const double weight = panel < panels.size() ? panels[panel].area : 0.0;
      phi += weight * (panel < fine.size() ? fine[panel] : 0.0);
      if (coupled && panel_count + panel < fine.size())
        phin += weight * fine[panel_count + panel];
    }
    coarse[a] = phi / denom;
    if (coupled)
      coarse[q + a] = phin / denom;
  }
}

void JcpCoarseGalerkinPreconditioner::prolongate(
    const std::vector<double> &coarse, std::vector<double> &fine) const {
  const std::size_t q = aggregates.size();
  fine.assign(coupled ? 2 * panel_count : panel_count, 0.0);
  for (std::size_t a = 0; a < q; ++a) {
    for (const std::size_t panel : aggregates[a].panels) {
      if (panel < fine.size())
        fine[panel] = a < coarse.size() ? coarse[a] : 0.0;
      if (coupled && panel_count + panel < fine.size())
        fine[panel_count + panel] =
            q + a < coarse.size() ? coarse[q + a] : 0.0;
    }
  }
}

std::size_t JcpCoarseGalerkinPreconditioner::estimated_bytes() const noexcept {
  const std::size_t dim = dimension();
  std::size_t bytes = inv_diag.size() * sizeof(double) +
                      aggregates.size() * sizeof(Aggregate) +
                      dim * dim * sizeof(double) + dim * sizeof(int);
  for (const Aggregate &aggregate : aggregates)
    bytes += aggregate.panels.size() * sizeof(std::size_t);
  return bytes;
}

void JcpCoarseGalerkinPreconditioner::apply(
    const std::vector<Panel> &panels, const std::vector<double> &in,
    std::vector<double> &out, TitanExecutor &executor) const {
  std::vector<double> coarse_rhs;
  restrict_area_weighted(panels, in, coarse_rhs);
  std::vector<double> projected;
  prolongate(coarse_rhs, projected);
  std::vector<double> complement(in.size(), 0.0);
  for (std::size_t i = 0; i < in.size(); ++i)
    complement[i] = in[i] - (i < projected.size() ? projected[i] : 0.0);
  if (block_base != nullptr) {
    block_base->apply_parallel(complement, out, executor);
  } else {
    out.resize(in.size());
    for (std::size_t i = 0; i < in.size(); ++i)
      out[i] = complement[i] * (i < inv_diag.size() ? inv_diag[i] : 1.0);
  }
  if (valid && lu.size() == dimension() * dimension()) {
    solve_dense_lu(lu, pivots, dimension(), coarse_rhs);
    prolongate(coarse_rhs, projected);
    for (std::size_t i = 0; i < out.size(); ++i)
      out[i] += i < projected.size() ? projected[i] : 0.0;
  }
}

std::vector<JcpBlockRange>
make_jcp_block_ranges(const TreeWorkspace &tree, std::size_t panel_count,
                      std::size_t block_panel_cap) {
  std::vector<JcpBlockRange> ranges;
  const std::size_t cap = std::clamp<std::size_t>(block_panel_cap, 8, 32);
  auto append_range = [&](std::size_t begin, std::size_t end) {
    if (begin < end && end <= tree.order.size())
      ranges.push_back({begin, end});
  };
  bool emitted_leaf_block = false;
  for (const TreeNode &node : tree.nodes) {
    if (node.child_count != 0 || node.panel_count == 0)
      continue;
    emitted_leaf_block = true;
    const std::size_t leaf_end = node.first_panel + node.panel_count;
    for (std::size_t begin = node.first_panel; begin < leaf_end; begin += cap)
      append_range(begin, std::min(leaf_end, begin + cap));
  }
  if (!emitted_leaf_block) {
    for (std::size_t begin = 0; begin < panel_count; begin += cap)
      append_range(begin, std::min(panel_count, begin + cap));
  }
  return ranges;
}

std::size_t estimate_jcp_block_preconditioner_bytes(
    const std::vector<JcpBlockRange> &ranges, std::size_t diagonal_size,
    bool coupled) noexcept {
  std::size_t bytes = diagonal_size * sizeof(double);
  for (const JcpBlockRange &range : ranges) {
    const std::size_t m = range.second - range.first;
    const std::size_t dim = coupled ? 2 * m : m;
    bytes += m * sizeof(std::size_t);
    bytes += dim * dim * sizeof(double);
    bytes += dim * sizeof(int);
  }
  return bytes;
}

std::shared_ptr<JcpBlockPreconditioner> build_jcp_block_preconditioner(
    const std::vector<Panel> &panels, const TreeWorkspace &tree, bool coupled,
    double eps_ratio, double kappa, double lhs_phi, double lhs_phin,
    const std::vector<double> &diagonal, const std::string &near_integration,
    const std::string &quad_buffer, const std::string &block_near_mode,
    const std::vector<JcpBlockRange> &ranges, TitanExecutor &executor) {
  const std::size_t n = panels.size();
  if (n == 0 || tree.order.size() != n)
    return {};
  auto precond = std::make_shared<JcpBlockPreconditioner>();
  precond->panel_count = n;
  precond->coupled = coupled;
  precond->inv_diag.assign(diagonal.size(), 1.0);
  for (std::size_t i = 0; i < diagonal.size(); ++i) {
    const double diag = std::abs(diagonal[i]);
    precond->inv_diag[i] =
        std::isfinite(diag) && diag > 1.0e-12 ? 1.0 / diag : 1.0;
  }

  const bool buffer_aware = block_near_mode == "buffer-aware" ||
                            block_near_mode == "operator-near";
  const bool operator_near = block_near_mode == "operator-near";
  auto build_block = [&](std::size_t begin, std::size_t end) {
    JcpBlockPreconditioner::Block block;
    const std::size_t count = end - begin;
    block.panels.reserve(count);
    for (std::size_t pos = begin; pos < end; ++pos)
      block.panels.push_back(tree.order[pos]);
    const std::size_t m = block.panels.size();
    const std::size_t dim = coupled ? 2 * m : m;
    block.lu.assign(dim * dim, 0.0);
    for (std::size_t r = 0; r < m; ++r) {
      const std::size_t target_index = block.panels[r];
      block.lu[r * dim + r] += lhs_phi;
      if (coupled)
        block.lu[(m + r) * dim + (m + r)] += lhs_phin;
      for (std::size_t c = 0; c < m; ++c) {
        const std::size_t source_index = block.panels[c];
        const bool singular = is_singular_near_pair(
            panels[target_index], target_index, panels[source_index],
            source_index);
        if (!singular && !buffer_aware)
          continue;
        if (!singular) {
          const GeometricNearDecision geometric =
              classify_geometric_near_buffer_pair(
                  panels[target_index], panels[source_index], quad_buffer);
          if (!geometric.active)
            continue;
        }
        double k1 = 0.0;
        double k2 = 0.0;
        double k3 = 0.0;
        double k4 = 0.0;
        double weight = panels[source_index].area;
        if (operator_near) {
          const JcpNearPairKernel selected = evaluate_jcp_near_pair_kernel(
              panels[target_index], target_index, panels[source_index],
              source_index, eps_ratio, kappa, coupled, near_integration,
              quad_buffer, nullptr);
          if (!selected.valid)
            continue;
          k1 = selected.k1;
          k2 = selected.k2;
          k3 = selected.k3;
          k4 = selected.k4;
          weight = 1.0;
        } else {
          const JcpKernels kernels = jcp_surface_kernels(
              panels[target_index], panels[source_index], eps_ratio, kappa);
          k1 = kernels.k1;
          k2 = kernels.k2;
          k3 = kernels.k3;
          k4 = kernels.k4;
        }
        if (coupled) {
          block.lu[r * dim + c] -= weight * k2;
          block.lu[r * dim + (m + c)] -= weight * k1;
          block.lu[(m + r) * dim + c] -= weight * k4;
          block.lu[(m + r) * dim + (m + c)] -= weight * k3;
        } else {
          block.lu[r * dim + c] -= weight * k2;
        }
      }
    }
    constexpr double kMinBlockPivotRatio = 1.0e-10;
    const std::vector<double> original_matrix = block.lu;
    block.scale = make_equilibration_scale(block.lu, dim);
    apply_symmetric_scale(block.lu, block.scale, dim);
    block.valid =
        factor_dense_lu(block.lu, block.pivots, dim, &block.pivot_ratio);
    block.ill_conditioned =
        !std::isfinite(block.pivot_ratio) ||
        block.pivot_ratio < kMinBlockPivotRatio;
    if (block.ill_conditioned)
      block.valid = false;
    if (block.valid) {
      block.valid = passes_block_quality_gate(
          original_matrix, block.lu, block.pivots, block.scale, dim,
          &block.quality_gain);
      block.quality_rejected = !block.valid;
    }
    return block;
  };

  precond->blocks.resize(ranges.size());
  const std::size_t workers = solver_workers(executor, ranges.size());
  executor.parallel_for(
      ranges.size(), workers, 1,
      [&](std::size_t, std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i)
          precond->blocks[i] = build_block(ranges[i].first, ranges[i].second);
      });
  for (const auto &block : precond->blocks) {
    if (block.ill_conditioned)
      ++precond->ill_conditioned_blocks;
    if (block.quality_rejected)
      ++precond->quality_rejected_blocks;
    precond->quality_probe_max_gain =
        std::max(precond->quality_probe_max_gain, block.quality_gain);
    if (block.valid) {
      ++precond->valid_blocks;
      precond->valid_panel_count += block.panels.size();
    }
  }
  precond->quality_status =
      precond->quality_rejected_blocks == 0 ? "pass" : "partial_reject";
  if (precond->valid_blocks == 0 && !precond->blocks.empty())
    precond->quality_status = "all_rejected";
  return precond;
}

std::shared_ptr<JcpCoarseGalerkinPreconditioner>
make_jcp_coarse_structure(const std::vector<Panel> &panels,
                          const TreeWorkspace &tree, bool coupled,
                          const std::vector<double> &diagonal,
                          std::size_t aggregate_cap) {
  auto coarse = std::make_shared<JcpCoarseGalerkinPreconditioner>();
  coarse->panel_count = panels.size();
  coarse->coupled = coupled;
  coarse->aggregates =
      make_jcp_coarse_aggregates(panels, tree, aggregate_cap);
  coarse->inv_diag.assign(diagonal.size(), 1.0);
  for (std::size_t i = 0; i < diagonal.size(); ++i) {
    const double diag = std::abs(diagonal[i]);
    coarse->inv_diag[i] =
        std::isfinite(diag) && diag > 1.0e-12 ? 1.0 / diag : 1.0;
  }
  return coarse;
}

bool build_jcp_exact_galerkin(
    const std::vector<Panel> &panels, JcpCoarseGalerkinPreconditioner &coarse,
    const std::function<void(const std::vector<double> &,
                             std::vector<double> &)> &setup_matvec,
    bool verify) {
  const std::size_t dim = coarse.dimension();
  if (dim == 0 || !setup_matvec)
    return false;
  std::vector<double> matrix(dim * dim, 0.0);
  std::vector<double> basis(dim, 0.0);
  std::vector<double> fine;
  std::vector<double> applied;
  std::vector<double> restricted;
  for (std::size_t col = 0; col < dim; ++col) {
    std::fill(basis.begin(), basis.end(), 0.0);
    basis[col] = 1.0;
    coarse.prolongate(basis, fine);
    setup_matvec(fine, applied);
    ++coarse.setup_matvec_calls;
    coarse.restrict_area_weighted(panels, applied, restricted);
    for (std::size_t row = 0; row < dim; ++row)
      matrix[row * dim + col] = restricted[row];
  }
  for (std::size_t col = 0; col < dim; ++col) {
    std::fill(basis.begin(), basis.end(), 0.0);
    basis[col] = 1.0;
    coarse.prolongate(basis, fine);
    coarse.restrict_area_weighted(panels, fine, restricted);
    for (std::size_t row = 0; row < dim; ++row) {
      const double target = row == col ? 1.0 : 0.0;
      coarse.rp_identity_error_max =
          std::max(coarse.rp_identity_error_max,
                   std::abs(restricted[row] - target));
    }
  }
  if (verify) {
    double diff_sq = 0.0;
    double ref_sq = 0.0;
    for (std::size_t col = 0; col < dim; ++col) {
      std::fill(basis.begin(), basis.end(), 0.0);
      basis[col] = 1.0;
      coarse.prolongate(basis, fine);
      setup_matvec(fine, applied);
      ++coarse.setup_matvec_calls;
      coarse.restrict_area_weighted(panels, applied, restricted);
      for (std::size_t row = 0; row < dim; ++row) {
        const double expected = matrix[row * dim + col];
        const double delta = restricted[row] - expected;
        diff_sq += delta * delta;
        ref_sq += expected * expected;
      }
    }
    coarse.verify_rel_error_max =
        std::sqrt(diff_sq) / std::max(1.0e-300, std::sqrt(ref_sq));
  }
  coarse.lu = std::move(matrix);
  coarse.valid =
      factor_dense_lu(coarse.lu, coarse.pivots, dim, &coarse.pivot_ratio);
  constexpr double kMinCoarsePivotRatio = 1.0e-12;
  if (!std::isfinite(coarse.pivot_ratio) ||
      coarse.pivot_ratio < kMinCoarsePivotRatio)
    coarse.valid = false;
  return coarse.valid;
}

} // namespace TitanPBE::Tabi
