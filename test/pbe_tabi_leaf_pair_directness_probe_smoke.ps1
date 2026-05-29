param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_leaf_pair_directness_probe_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_leaf_pair_directness_probe_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_leaf_pair_directness_probe_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside directness-probe smoke root: $full"
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
  p.area = 0.18 + 0.01 * static_cast<double>(index % 13);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(160);
  std::size_t index = 0;
  const double centers[5][3] = {{-24.0, 0.0, 0.0}, {24.0, 0.0, 0.0},
                                {0.0, 24.0, 0.0}, {0.0, -24.0, 0.0},
                                {0.0, 0.0, 26.0}};
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
    const std::string &directness_probe,
    const std::string &directness_verify) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 64, 0.05, "exact", "tile-avx2", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "off", 1.5, "off",
          "area-centroid", "off", "off", 512, 2048, "translated", "naive",
          "off", "off", 0.0, 0.0, "off", 1.0e4, 1.0e-3, 0.25,
          "stage3j", "leaf-pair-balanced", "off", 4, 16, 2);
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  tree.leaf_pair_kernel = "target-avx2";
  tree.leaf_pair_selection_profile = "summary";
  tree.leaf_pair_probe_cache = "off";
  tree.leaf_pair_target_group = "same-target-leaf";
  tree.leaf_pair_directness_probe = directness_probe;
  tree.leaf_pair_directness_cache = "off";
  tree.leaf_pair_directness_verify = directness_verify;
  tree.leaf_pair_max_lookahead = 32;
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
  const std::vector<TitanPBE::Tabi::Panel> panels = build_fixture();
  const std::vector<double> unknowns = make_unknowns(panels.size());

  TitanPBE::Tabi::TreeWorkspace full_tree =
      make_tree(panels, executor, "full-traversal", "off");
  TitanPBE::Tabi::TreeWorkspace ancestor_tree =
      make_tree(panels, executor, "ancestor-path", "off");
  TitanPBE::Tabi::TreeWorkspace verify_tree =
      make_tree(panels, executor, "ancestor-path", "full-small");

  std::vector<double> full_out;
  std::vector<double> ancestor_out;
  std::vector<double> verify_out;
  const TitanPBE::Tabi::TreecodeStats full_stats =
      apply_fixture(panels, unknowns, executor, full_tree, full_out);
  const TitanPBE::Tabi::TreecodeStats ancestor_stats =
      apply_fixture(panels, unknowns, executor, ancestor_tree, ancestor_out);
  const TitanPBE::Tabi::TreecodeStats verify_stats =
      apply_fixture(panels, unknowns, executor, verify_tree, verify_out);

  const double rel = relative_rms(full_out, ancestor_out);
  const double verify_rel = relative_rms(full_out, verify_out);
  require_true("ancestor-path keeps full traversal accuracy", rel <= 1.0e-12);
  require_true("ancestor-path verify row keeps full traversal accuracy",
               verify_rel <= 1.0e-12);
  require_true("full traversal probes were counted",
               full_stats.leaf_pair_selection_probe_calls > 0);
  require_true("ancestor probes were counted",
               ancestor_stats.leaf_pair_directness_ancestor_probe_calls > 0);
  require_true("ancestor path reduces node visits",
               ancestor_stats.leaf_pair_selection_probe_node_visits <
                   full_stats.leaf_pair_selection_probe_node_visits);
  require_true("full-small verification checked candidates",
               verify_stats.leaf_pair_directness_verify_checks > 0);
  require_true("ancestor-path has no directness mismatches",
               verify_stats.leaf_pair_directness_verify_mismatches == 0);

  std::cout << "PBE TABI leaf-pair directness-probe smoke PASS panels="
            << panels.size() << " rel=" << rel
            << " verify_rel=" << verify_rel
            << " full_visits="
            << full_stats.leaf_pair_selection_probe_node_visits
            << " ancestor_visits="
            << ancestor_stats.leaf_pair_selection_probe_node_visits
            << " verify_checks="
            << verify_stats.leaf_pair_directness_verify_checks << "\n";
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
    throw "Failed to compile TABI leaf-pair directness-probe smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI leaf-pair directness-probe smoke fixture failed."
}
