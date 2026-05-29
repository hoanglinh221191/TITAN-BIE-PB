param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_leaf_pair_target_avx2_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_leaf_pair_target_avx2_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_leaf_pair_target_avx2_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside target-avx2 smoke root: $full"
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
#include "titan_pbe/tabi_types.h"

#include <algorithm>
#include <cmath>
#include <iostream>
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
  require_true("same vector size", reference.size() == candidate.size());
  std::vector<double> diff(reference.size(), 0.0);
  for (std::size_t i = 0; i < reference.size(); ++i)
    diff[i] = candidate[i] - reference[i];
  return rms(diff) / std::max(1.0e-12, rms(reference));
}

TitanPBE::Tabi::Panel make_panel(double x, double y, double z, double nx,
                                 double ny, double nz, std::size_t index) {
  TitanPBE::Tabi::Panel p;
  p.x = x;
  p.y = y;
  p.z = z;
  const double norm =
      std::sqrt(std::max(1.0e-12, nx * nx + ny * ny + nz * nz));
  p.nx = nx / norm;
  p.ny = ny / norm;
  p.nz = nz / norm;
  p.area = 0.20 + 0.01 * static_cast<double>(index % 11);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_clean_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(128);
  std::size_t index = 0;
  const double centers[4][3] = {{-22.0, 0.0, 0.0}, {22.0, 0.0, 0.0},
                                {0.0, 24.0, 0.0}, {0.0, -24.0, 0.0}};
  for (const auto &center : centers) {
    for (int i = 0; i < 32; ++i) {
      const double a = 0.19634954084936207 * static_cast<double>(i);
      const double r = 1.0 + 0.04 * static_cast<double>(i % 7);
      panels.push_back(make_panel(center[0] + r * std::cos(a),
                                  center[1] + r * std::sin(a),
                                  center[2] + 0.08 * std::cos(3.0 * a),
                                  0.15 + 0.01 * i, 0.2, 1.0, index++));
    }
  }
  return panels;
}

std::vector<TitanPBE::Tabi::Panel> build_singular_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels = build_clean_fixture();
  for (std::size_t i = 0; i < panels.size(); ++i) {
    panels[i].has_triangle_geometry = true;
    panels[i].v0_index = i;
    panels[i].v1_index = i + 1000;
    panels[i].v2_index = i + 2000;
  }
  for (std::size_t i = 1; i < panels.size(); i += 8)
    panels[i].v0_index = panels[i - 1].v0_index;
  return panels;
}

std::vector<double> make_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    unknowns[i] = 0.7 + 0.20 * std::sin(0.19 * static_cast<double>(i));
    unknowns[n + i] = -1.1 + 0.15 * std::cos(0.27 * static_cast<double>(i));
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &cluster_policy, const std::string &leaf_pair_kernel) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 64, 0.05, "exact", "tile-avx2", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "off", 1.5, "off",
          "area-centroid", "off", "off", 512, 2048, "translated", "naive",
          "off", "off", 0.0, 0.0, "off", 1.0e4, 1.0e-3, 0.25,
          "stage3j", cluster_policy, "off", 4, 16, 2);
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  tree.leaf_pair_kernel = leaf_pair_kernel;
  return tree;
}

TitanPBE::Tabi::TreecodeStats apply_fixture(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    TitanPBE::Tabi::TreeWorkspace &tree, std::vector<double> &out) {
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.10, 0.0, 0.0, true, tree, "centroid", "off",
      executor, out, &stats);
  return stats;
}

} // namespace

int main() {
  TitanExecutor executor(1);

  const std::vector<TitanPBE::Tabi::Panel> clean_panels =
      build_clean_fixture();
  const std::vector<double> clean_unknowns = make_unknowns(clean_panels.size());

  TitanPBE::Tabi::TreeWorkspace scalar_tree =
      make_tree(clean_panels, executor, "leaf-pair-balanced", "scalar");
  TitanPBE::Tabi::TreeWorkspace avx2_tree =
      make_tree(clean_panels, executor, "leaf-pair-balanced", "target-avx2");
  std::vector<double> scalar_out;
  std::vector<double> avx2_out;
  const TitanPBE::Tabi::TreecodeStats scalar_stats =
      apply_fixture(clean_panels, clean_unknowns, executor, scalar_tree,
                    scalar_out);
  const TitanPBE::Tabi::TreecodeStats avx2_stats =
      apply_fixture(clean_panels, clean_unknowns, executor, avx2_tree,
                    avx2_out);
  const double rel = relative_rms(scalar_out, avx2_out);

  require_true("target-avx2 keeps scalar batch accuracy", rel <= 1.0e-12);
  require_true("scalar batch interactions are counted",
               scalar_stats.leaf_pair_batch_interactions > 0);
  require_true("target-avx2 batch interactions are counted",
               avx2_stats.leaf_pair_batch_interactions > 0);
  require_true("target-avx2 kernel is reported",
               avx2_stats.leaf_pair_batch_kernel == "target-avx2");
  require_true("target-avx2 vectorized sources counted",
               avx2_stats.leaf_pair_batch_vectorized_sources > 0);
  require_true("target-avx2 vector lanes counted",
               avx2_stats.leaf_pair_batch_vector_lanes > 0);

  const std::vector<TitanPBE::Tabi::Panel> singular_panels =
      build_singular_fixture();
  const std::vector<double> singular_unknowns =
      make_unknowns(singular_panels.size());
  TitanPBE::Tabi::TreeWorkspace singular_scalar_tree =
      make_tree(singular_panels, executor, "leaf-pair-strict", "scalar");
  TitanPBE::Tabi::TreeWorkspace singular_tree =
      make_tree(singular_panels, executor, "leaf-pair-strict", "target-avx2");
  std::vector<double> singular_scalar_out;
  std::vector<double> singular_out;
  const TitanPBE::Tabi::TreecodeStats singular_scalar_stats =
      apply_fixture(singular_panels, singular_unknowns, executor,
                    singular_scalar_tree, singular_scalar_out);
  const TitanPBE::Tabi::TreecodeStats singular_stats =
      apply_fixture(singular_panels, singular_unknowns, executor,
                    singular_tree, singular_out);
  const double singular_rel = relative_rms(singular_scalar_out, singular_out);
  require_true("singular/exclusion cases fall back",
               singular_stats.leaf_pair_batch_fallbacks > 0);
  require_true("singular/exclusion fallback remains accurate",
               singular_rel <= 1.0e-12);

  std::cout << "PBE TABI leaf-pair target-avx2 smoke PASS panels="
            << clean_panels.size() << " rel=" << rel << " batches="
            << avx2_stats.leaf_pair_batch_interactions << " vector_sources="
            << avx2_stats.leaf_pair_batch_vectorized_sources
            << " vector_lanes=" << avx2_stats.leaf_pair_batch_vector_lanes
            << " fallbacks=" << singular_stats.leaf_pair_batch_fallbacks
            << " singular_rel=" << singular_rel
            << "\n";
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
    "-I$RepoRoot",
    (Join-Path $RepoRoot 'titan_pbe\tabi_near.cpp'),
    (Get-PbeTabiTreecodeSources -RepoRoot $RepoRoot),
    (Join-Path $RepoRoot 'titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI leaf-pair target-avx2 smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI leaf-pair target-avx2 smoke fixture failed."
}
