param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_large_protein_policy_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_large_protein_policy_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_large_protein_policy_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI large-protein smoke root: $full"
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

std::vector<TitanPBE::Tabi::Panel> build_panels() {
  TitanPBE::Tabi::Options opt;
  opt.surface_mode = "vdw";
  opt.surface_geometry = "triangle";
  opt.diagonal_mode = "solid-angle";
  opt.surface_point_count = 24;
  TitanPBE::Tabi::Result result;
  std::vector<TitanPBE::Tabi::Atom> atoms;
  for (int i = 0; i < 14; ++i) {
    const double x = -8.4 + 1.28 * static_cast<double>(i);
    const double y = 0.45 * std::sin(0.7 * static_cast<double>(i));
    const double z = 0.38 * std::cos(0.5 * static_cast<double>(i));
    const double q = (i % 2 == 0) ? 0.35 : -0.32;
    atoms.push_back({x, y, z, q, 1.48 + 0.05 * static_cast<double>(i % 5)});
  }
  std::vector<TitanPBE::Tabi::Panel> panels =
      TitanPBE::Tabi::build_surface(atoms, opt, &result);
  require_true("triangle panels", panels.size() > 250);
  return panels;
}

} // namespace

int main() {
  std::vector<TitanPBE::Tabi::Panel> panels = build_panels();
  const std::size_t n = panels.size();
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    unknowns[i] = 0.05 * std::sin(0.011 * static_cast<double>(i));
    unknowns[n + i] = 0.02 * std::cos(0.017 * static_cast<double>(i));
  }

  TitanExecutor executor1(1);
  TitanPBE::Tabi::TreeWorkspace exact_tree =
      TitanPBE::Tabi::build_tree_workspace(panels, executor1, 32, 0.01,
                                           "exact", "tile-avx2");
  TitanPBE::Tabi::TreeWorkspace masked_tree =
      TitanPBE::Tabi::build_tree_workspace(panels, executor1, 32, 0.01,
                                           "large-protein", "tile-avx2");
  require_true("large policy normalized",
               masked_tree.leaf_near_policy == "large-protein");
  require_true("span tile emitted",
               masked_tree.interactions.leaf_tile_span_firsts.size() > 0);
  require_true("masked tile reduces flat residual",
               masked_tree.interactions.leaf_indices.size() <
                   exact_tree.interactions.leaf_indices.size());

  std::vector<double> exact_out;
  std::vector<double> masked_out;
  TitanPBE::Tabi::TreecodeStats exact_stats;
  TitanPBE::Tabi::TreecodeStats masked_stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 0.0125, 0.12, 0.50625, 40.5, true, exact_tree,
      "centroid", executor1, exact_out, &exact_stats);
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 0.0125, 0.12, 0.50625, 40.5, true, masked_tree,
      "centroid", executor1, masked_out, &masked_stats);
  require_true("masked stats policy",
               masked_stats.leaf_near_policy == "large-protein");
  require_true("masked exact agreement",
               relative_rms(exact_out, masked_out) <= 1.0e-12);
  require_true("masked residual stats reduced",
               masked_stats.leaf_residual_panel_interactions <
                   exact_stats.leaf_residual_panel_interactions);
  require_true("large-protein smooth disabled",
               masked_stats.smooth_leaf_node_interactions == 0);

  TitanExecutor executor4(4);
  TitanPBE::Tabi::TreeWorkspace large_tree4 =
      TitanPBE::Tabi::build_tree_workspace(panels, executor4, 32, 0.01,
                                           "large-protein", "tile-avx2");
  std::vector<double> large_out4;
  TitanPBE::Tabi::TreecodeStats large_stats4;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 0.0125, 0.12, 0.50625, 40.5, true, large_tree4,
      "centroid", executor4, large_out4, &large_stats4);
  require_true("large nt deterministic",
               relative_rms(masked_out, large_out4) < 1.0e-14);

  std::cout << "PBE TABI large-protein policy smoke PASS panels=" << n
            << " residual_exact="
            << exact_stats.leaf_residual_panel_interactions
            << " residual_large="
            << masked_stats.leaf_residual_panel_interactions
            << " tile_spans="
            << masked_tree.interactions.leaf_tile_span_firsts.size()
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
    ("-I" + (Join-Path $RepoRoot "src")),
    (Join-Path $RepoRoot 'src\tabi\tabi_surface.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_near.cpp'),
    (Get-PbeTabiTreecodeSources -RepoRoot $RepoRoot),
    (Join-Path $RepoRoot 'src\titan_surface.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI large-protein policy fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI large-protein policy fixture failed."
}
