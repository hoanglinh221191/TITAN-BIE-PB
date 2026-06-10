#include "tabi/tabi_preconditioner.h"

#include "titan_executor.h"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using TitanPBE::Tabi::JcpBlockRange;
using TitanPBE::Tabi::Panel;
using TitanPBE::Tabi::TreeNode;
using TitanPBE::Tabi::TreeWorkspace;

namespace {

Panel make_panel(std::size_t i) {
  const double x = static_cast<double>(i % 4);
  const double y = static_cast<double>(i / 4);
  Panel p;
  p.x = x + 0.25;
  p.y = y + 0.25;
  p.z = 0.05 * static_cast<double>(i);
  p.nx = 0.0;
  p.ny = 0.0;
  p.nz = 1.0;
  p.area = 0.5;
  p.v0x = x;
  p.v0y = y;
  p.v0z = p.z;
  p.v1x = x + 0.5;
  p.v1y = y;
  p.v1z = p.z;
  p.v2x = x;
  p.v2y = y + 0.5;
  p.v2z = p.z;
  p.v0_index = 3 * i;
  p.v1_index = 3 * i + 1;
  p.v2_index = 3 * i + 2;
  p.has_triangle_geometry = true;
  return p;
}

TreeWorkspace make_tree(std::size_t n) {
  TreeWorkspace tree;
  tree.order.resize(n);
  std::iota(tree.order.begin(), tree.order.end(), std::size_t{0});
  TreeNode root;
  root.first_panel = 0;
  root.panel_count = n;
  root.first_child = 0;
  root.child_count = 0;
  tree.root_index = 0;
  tree.nodes.push_back(root);
  tree.leaf_count = 1;
  return tree;
}

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

} // namespace

int main() {
  constexpr std::size_t n = 8;
  std::vector<Panel> panels;
  panels.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    panels.push_back(make_panel(i));
  TreeWorkspace tree = make_tree(n);
  std::vector<double> diagonal(2 * n, 1.0);
  TitanExecutor executor(2);

  const std::vector<JcpBlockRange> ranges =
      TitanPBE::Tabi::make_jcp_block_ranges(tree, n, 4);
  require(ranges.size() == 1, "expected one clamped block range");
  const std::size_t estimated =
      TitanPBE::Tabi::estimate_jcp_block_preconditioner_bytes(
          ranges, diagonal.size(), true);
  require(estimated > diagonal.size() * sizeof(double),
          "expected positive block memory estimate");

  auto precond = TitanPBE::Tabi::build_jcp_block_preconditioner(
      panels, tree, true, 80.0, 0.125, 40.5, 0.50625, diagonal, "centroid",
      "adaptive", "buffer-aware", ranges, executor);
  require(precond != nullptr, "preconditioner build returned null");
  require(precond->valid_blocks > 0, "expected at least one valid block");
  require(precond->valid_panel_count == n, "expected full panel coverage");
  require(precond->estimated_bytes() <= estimated,
          "real allocation should not exceed estimate");

  std::vector<double> in(2 * n, 0.0);
  for (std::size_t i = 0; i < in.size(); ++i)
    in[i] = 0.1 + static_cast<double>(i);
  std::vector<double> out;
  precond->apply_parallel(in, out, executor);
  require(out.size() == in.size(), "apply output size mismatch");
  for (double value : out)
    require(std::isfinite(value), "non-finite preconditioner output");

  std::cout << "TABI preconditioner smoke PASS\n";
  return 0;
}
