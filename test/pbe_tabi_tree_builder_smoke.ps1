param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_tree_builder_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_tree_builder_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_tree_builder_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI tree-builder smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

if (-not (Test-Path -LiteralPath $Compiler)) {
    throw "Missing clang++ compiler: $Compiler"
}

New-CleanDir $WorkRoot

$source = @'
#include "titan_executor.h"
#include "tabi/tabi_types.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require_true(const std::string &label, bool value) {
  if (!value)
    throw std::runtime_error(label);
}

double rms(const std::vector<double> &values) {
  double sum = 0.0;
  for (double value : values)
    sum += value * value;
  return values.empty() ? 0.0
                        : std::sqrt(sum / static_cast<double>(values.size()));
}

double relative_rms(const std::vector<double> &reference,
                    const std::vector<double> &candidate) {
  std::vector<double> diff(reference.size(), 0.0);
  for (std::size_t i = 0; i < reference.size(); ++i)
    diff[i] = candidate[i] - reference[i];
  return rms(diff) / std::max(1.0e-12, rms(reference));
}

TitanPBE::Tabi::Panel make_panel(double x, double y, double z,
                                 std::size_t index) {
  TitanPBE::Tabi::Panel p;
  p.x = x;
  p.y = y;
  p.z = z;
  const double norm = std::sqrt(std::max(1.0e-12, x * x + y * y + z * z));
  p.nx = x / norm;
  p.ny = y / norm;
  p.nz = z / norm;
  p.area = 0.25 + 0.01 * static_cast<double>(index % 7);
  p.diagonal = 1.0;
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_jump_panels() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(512);
  std::size_t index = 0;
  for (int cluster = 0; cluster < 8; ++cluster) {
    const double ox = (cluster & 4) ? 18.0 : -18.0;
    const double oy = (cluster & 2) ? 7.0 : -7.0;
    const double oz = (cluster & 1) ? 3.0 : -3.0;
    for (int ix = 0; ix < 4; ++ix) {
      for (int iy = 0; iy < 4; ++iy) {
        for (int iz = 0; iz < 4; ++iz) {
          panels.push_back(make_panel(
              ox + 0.35 * static_cast<double>(ix),
              oy + 0.31 * static_cast<double>(iy),
              oz + 0.27 * static_cast<double>(iz), index++));
        }
      }
    }
  }
  return panels;
}

void check_tree(const TitanPBE::Tabi::TreeWorkspace &tree,
                std::size_t panel_count, const std::string &builder) {
  require_true(builder + " builder label", tree.tree_builder == builder);
  require_true(builder + " order size", tree.order.size() == panel_count);
  require_true(builder + " nodes", tree.nodes.size() >= tree.leaf_count);
  require_true(builder + " leaves", tree.leaf_count > 0);
  require_true(builder + " root", tree.root_index < tree.nodes.size());
  std::vector<std::size_t> sorted = tree.order;
  std::sort(sorted.begin(), sorted.end());
  for (std::size_t i = 0; i < sorted.size(); ++i)
    require_true(builder + " order permutation", sorted[i] == i);
  for (const auto &node : tree.nodes) {
    require_true(builder + " finite radius", std::isfinite(node.radius));
    require_true(builder + " finite center",
                 std::isfinite(node.cx) && std::isfinite(node.cy) &&
                     std::isfinite(node.cz));
    if (node.child_count == 0)
      continue;
    require_true(builder + " child offset",
                 node.first_child + node.child_count <=
                     tree.child_indices.size());
    for (std::size_t child = 0; child < node.child_count; ++child) {
      const std::size_t child_index = tree.child_indices[node.first_child + child];
      require_true(builder + " child index", child_index < tree.nodes.size());
    }
  }
}

} // namespace

int main() {
  const std::vector<TitanPBE::Tabi::Panel> panels = build_jump_panels();
  const std::size_t n = panels.size();
  TitanExecutor executor1(1);
  TitanExecutor executor4(4);

  TitanPBE::Tabi::TreeWorkspace grouped =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor1, 8, 0.7, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "grouped-morton");
  TitanPBE::Tabi::TreeWorkspace topdown =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor1, 8, 0.7, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb");
  check_tree(grouped, n, "grouped-morton");
  check_tree(topdown, n, "topdown-aabb");
  require_true("topdown depth", topdown.depth_max > 0);
  require_true("topdown child indirection",
               topdown.child_indices.size() >= topdown.nodes.size() - 1);
  require_true("radius diagnostic finite",
               std::isfinite(topdown.node_radius_mean) &&
                   std::isfinite(topdown.node_radius_max));

  std::vector<double> density(n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    density[i] = std::sin(0.011 * static_cast<double>(i)) +
                 0.2 * std::cos(0.017 * static_cast<double>(i));

  std::vector<double> topdown_out1;
  TitanPBE::Tabi::TreecodeStats topdown_stats1;
  TitanPBE::Tabi::apply_treecode_operator(panels, density, 0.42, 0.08,
                                          topdown, executor1, topdown_out1,
                                          &topdown_stats1);
  require_true("topdown stats builder",
               topdown_stats1.tree_builder == "topdown-aabb");
  require_true("topdown matvec finite",
               std::all_of(topdown_out1.begin(), topdown_out1.end(),
                           [](double v) { return std::isfinite(v); }));

  TitanPBE::Tabi::TreeWorkspace topdown4 =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor4, 8, 0.7, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb");
  std::vector<double> topdown_out4;
  TitanPBE::Tabi::TreecodeStats topdown_stats4;
  TitanPBE::Tabi::apply_treecode_operator(panels, density, 0.42, 0.08,
                                          topdown4, executor4, topdown_out4,
                                          &topdown_stats4);
  require_true("topdown nt deterministic",
               relative_rms(topdown_out1, topdown_out4) < 1.0e-14);

  std::cout << "PBE TABI tree-builder smoke PASS panels=" << n
            << " grouped_nodes=" << grouped.nodes.size()
            << " topdown_nodes=" << topdown.nodes.size()
            << " grouped_radius_mean=" << grouped.node_radius_mean
            << " topdown_radius_mean=" << topdown.node_radius_mean
            << " topdown_depth=" << topdown.depth_max << "\n";
  return 0;
}
'@

Set-Content -LiteralPath $SourcePath -Value $source -Encoding ASCII

$compileArgs = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1',
    ("-I" + (Join-Path $RepoRoot "src")),
    (Join-Path $RepoRoot 'src\tabi\tabi_near.cpp'),
    (Get-PbeTabiTreecodeSources -RepoRoot $RepoRoot),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI tree-builder smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI tree-builder smoke fixture failed."
}
