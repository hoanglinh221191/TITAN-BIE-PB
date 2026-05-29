param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_leaf_pair_worklist_scheduler_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_leaf_pair_worklist_scheduler_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_leaf_pair_worklist_scheduler_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside worklist scheduler smoke root: $full"
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

std::vector<TitanPBE::Tabi::Panel> build_clean_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(256);
  std::size_t index = 0;
  const double centers[8][3] = {{-30.0, 0.0, 0.0}, {30.0, 0.0, 0.0},
                                {0.0, -30.0, 0.0}, {0.0, 30.0, 0.0},
                                {0.0, 0.0, -32.0}, {0.0, 0.0, 32.0},
                                {-24.0, 24.0, 0.0}, {24.0, -24.0, 0.0}};
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

std::vector<TitanPBE::Tabi::Panel> build_exclusion_fixture() {
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
    const std::string &cluster_policy, const std::string &scheduler,
    const std::string &verify, int worklist_max_active,
    int min_source_panels) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 16, 0.05, "exact", "tile-avx2", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "off", 1.5, "off",
          "area-centroid", "off", "off", 512, 2048, "translated", "naive",
          "off", "off", 0.0, 0.0, "off", 1.0e4, 1.0e-3, 0.25,
          "stage3j", cluster_policy, "off", 4, min_source_panels, 2);
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  tree.leaf_pair_kernel = "target-avx2";
  tree.leaf_pair_selection_profile = "summary";
  tree.leaf_pair_probe_cache = "per-matvec";
  tree.leaf_pair_target_group = "same-target-leaf";
  tree.leaf_pair_directness_probe = "ancestor-path";
  tree.leaf_pair_directness_cache = "target-leaf-source-leaf";
  tree.leaf_pair_scheduler = scheduler;
  tree.leaf_pair_worklist_max_active =
      static_cast<std::size_t>(std::max(1, worklist_max_active));
  tree.leaf_pair_worklist_verify = verify;
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
      make_tree(clean_panels, executor, "off", "inline-select", "off", 2048,
                16);
  TitanPBE::Tabi::TreeWorkspace worklist_tree =
      make_tree(clean_panels, executor, "leaf-pair-balanced",
                "deferred-worklist", "full-small", 2048, 16);
  TitanPBE::Tabi::TreeWorkspace forced_tree =
      make_tree(clean_panels, executor, "leaf-pair-balanced",
                "deferred-worklist", "full-small", 1, 16);

  std::vector<double> scalar_out;
  std::vector<double> worklist_out;
  std::vector<double> forced_out;
  const TitanPBE::Tabi::TreecodeStats scalar_stats =
      apply_fixture(clean_panels, clean_unknowns, executor, scalar_tree,
                    scalar_out);
  const TitanPBE::Tabi::TreecodeStats worklist_stats =
      apply_fixture(clean_panels, clean_unknowns, executor, worklist_tree,
                    worklist_out);
  const TitanPBE::Tabi::TreecodeStats forced_stats =
      apply_fixture(clean_panels, clean_unknowns, executor, forced_tree,
                    forced_out);

  const double worklist_rel = relative_rms(scalar_out, worklist_out);
  const double forced_rel = relative_rms(scalar_out, forced_out);
  require_true("deferred worklist keeps scalar traversal accuracy",
               worklist_rel <= 1.0e-12);
  require_true("forced flush worklist keeps scalar traversal accuracy",
               forced_rel <= 1.0e-12);
  require_true("deferred scheduler mode is reported",
               worklist_stats.leaf_pair_scheduler == "deferred-worklist");
  require_true("worklist enqueued tasks",
               worklist_stats.leaf_pair_worklist_tasks > 0);
  require_true("worklist grouped tasks",
               worklist_stats.leaf_pair_worklist_groups > 0);
  require_true("worklist batched clean groups",
               worklist_stats.leaf_pair_worklist_batched_groups > 0);
  require_true("worklist vectorized source panels",
               worklist_stats.leaf_pair_batch_vectorized_sources > 0);
  require_true("worklist full-small verification checked tasks",
               worklist_stats.leaf_pair_worklist_verify_checks > 0);
  require_true("worklist verification has no mismatches",
               worklist_stats.leaf_pair_worklist_verify_mismatches == 0);
  require_true("scalar baseline does not use worklist",
               scalar_stats.leaf_pair_worklist_tasks == 0);
  require_true("small active limit exercises forced flush path",
               forced_stats.leaf_pair_worklist_forced_flushes > 0);
  require_true("forced flush verification has no mismatches",
               forced_stats.leaf_pair_worklist_verify_mismatches == 0);

  const std::vector<TitanPBE::Tabi::Panel> exclusion_panels =
      build_exclusion_fixture();
  const std::vector<double> exclusion_unknowns =
      make_unknowns(exclusion_panels.size());
  TitanPBE::Tabi::TreeWorkspace exclusion_scalar_tree =
      make_tree(exclusion_panels, executor, "off", "inline-select", "off",
                2048, 16);
  TitanPBE::Tabi::TreeWorkspace exclusion_worklist_tree =
      make_tree(exclusion_panels, executor, "leaf-pair-balanced",
                "deferred-worklist", "full-small", 2048, 16);
  std::vector<double> exclusion_scalar_out;
  std::vector<double> exclusion_worklist_out;
  const TitanPBE::Tabi::TreecodeStats exclusion_scalar_stats =
      apply_fixture(exclusion_panels, exclusion_unknowns, executor,
                    exclusion_scalar_tree, exclusion_scalar_out);
  const TitanPBE::Tabi::TreecodeStats exclusion_stats =
      apply_fixture(exclusion_panels, exclusion_unknowns, executor,
                    exclusion_worklist_tree, exclusion_worklist_out);
  const double exclusion_rel =
      relative_rms(exclusion_scalar_out, exclusion_worklist_out);
  require_true("singular/exclusion fallback remains accurate",
               exclusion_rel <= 1.0e-12);
  require_true("singular/exclusion fixture still uses direct leaf work",
               exclusion_stats.leaf_panel_interactions > 0 &&
                   exclusion_scalar_stats.leaf_panel_interactions > 0);

  std::cout << "PBE TABI leaf-pair worklist scheduler smoke PASS panels="
            << clean_panels.size() << " rel=" << worklist_rel
            << " forced_rel=" << forced_rel
            << " tasks=" << worklist_stats.leaf_pair_worklist_tasks
            << " groups=" << worklist_stats.leaf_pair_worklist_groups
            << " batched="
            << worklist_stats.leaf_pair_worklist_batched_groups
            << " forced_flushes="
            << forced_stats.leaf_pair_worklist_forced_flushes
            << " verify_checks="
            << worklist_stats.leaf_pair_worklist_verify_checks
            << " exclusion_rel=" << exclusion_rel << "\n";
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
    throw "Failed to compile TABI leaf-pair worklist scheduler smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI leaf-pair worklist scheduler smoke fixture failed."
}
