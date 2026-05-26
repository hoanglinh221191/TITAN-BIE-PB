param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_convergence_estimator_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_convergence_estimator_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_convergence_estimator_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside convergence estimator smoke root: $full"
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

TitanPBE::Tabi::Panel make_panel(double x, double y, double z,
                                 double nx, double ny, double nz,
                                 std::size_t index) {
  TitanPBE::Tabi::Panel p;
  p.x = x;
  p.y = y;
  p.z = z;
  const double norm =
      std::sqrt(std::max(1.0e-12, nx * nx + ny * ny + nz * nz));
  p.nx = nx / norm;
  p.ny = ny / norm;
  p.nz = nz / norm;
  p.area = 0.35 + 0.015 * static_cast<double>(index % 17);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_tail_heavy_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(96);
  std::size_t index = 0;

  for (int i = 0; i < 48; ++i) {
    const double x = -14.0 + 0.58 * static_cast<double>(i);
    const double y = (i % 2 == 0 ? -0.24 : 0.24) +
                     0.03 * static_cast<double>((i / 2) % 3);
    const double z = 0.05 * static_cast<double>((i % 7) - 3);
    panels.push_back(make_panel(x, y, z, 0.1, 0.2, 1.0, index++));
  }

  for (int i = 0; i < 24; ++i) {
    const double angle = 0.2617993877991494 * static_cast<double>(i);
    panels.push_back(make_panel(27.0 + 2.8 * std::cos(angle),
                                4.0 + 2.8 * std::sin(angle),
                                0.15 * std::sin(2.0 * angle),
                                -0.2, 0.1, 1.0, index++));
  }

  for (int i = 0; i < 24; ++i) {
    const double angle = 0.2617993877991494 * static_cast<double>(i);
    panels.push_back(make_panel(41.0 + 2.2 * std::cos(angle),
                                -5.0 + 2.2 * std::sin(angle),
                                0.10 * std::cos(3.0 * angle),
                                0.2, -0.1, 1.0, index++));
  }

  return panels;
}

std::vector<double> make_unknowns(std::size_t panel_count) {
  std::vector<double> unknowns(panel_count * 2, 0.0);
  for (std::size_t i = 0; i < panel_count; ++i) {
    const double alternating = (i % 2 == 0) ? 1.0 : -1.0;
    const double burst = (i % 9 == 0) ? 8.0 : 1.0;
    unknowns[2 * i] =
        alternating * burst *
        (0.9 + 0.25 * std::sin(0.37 * static_cast<double>(i)));
    unknowns[2 * i + 1] =
        -alternating * (1.4 + 0.35 * std::cos(0.23 * static_cast<double>(i)));
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &convergence_safety,
    const std::string &convergence_estimator, double tail_threshold,
    double error_cap, double unstable_rho_cap) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.70, "exact", "scalar", "strict",
          "adaptive", "onthefly", "tight-aabb", "topdown-aabb",
          "physical-gap", 1.5, "off", "area-centroid", "split", "sampled",
          2048, 8192, "translated", "naive", "sampled",
          convergence_safety, 0.0, 0.0, convergence_estimator,
          tail_threshold, error_cap, unstable_rho_cap);
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  return tree;
}

TitanPBE::Tabi::TreecodeStats apply_fixture(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    TitanPBE::Tabi::TreeWorkspace &tree, std::vector<double> &out) {
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.10, 0.0, 0.0, true, tree, "centroid",
      "off", executor, out, &stats);
  return stats;
}

} // namespace

int main() {
  const std::vector<TitanPBE::Tabi::Panel> panels =
      build_tail_heavy_fixture();
  const std::vector<double> unknowns = make_unknowns(panels.size());
  TitanExecutor executor(2);

  std::vector<double> direct;
  TitanPBE::Tabi::TreecodeStats direct_stats;
  TitanPBE::Tabi::TreeWorkspace empty_tree;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.10, 0.0, 0.0, true, empty_tree, "centroid",
      "off", executor, direct, &direct_stats);
  require_true("direct fallback mode",
               direct_stats.matvec_mode == "jcp2013_direct_fallback");

  TitanPBE::Tabi::TreeWorkspace audit_tree =
      make_tree(panels, executor, "audit", "hybrid", 1.0, 1.0e-8, 0.25);
  TitanPBE::Tabi::TreeWorkspace split_tree =
      make_tree(panels, executor, "split", "hybrid", 1.0, 1.0e-8, 0.25);
  TitanPBE::Tabi::TreeWorkspace loose_tree =
      make_tree(panels, executor, "split", "moment-tail", 1.0e12, 1.0,
                0.25);

  std::vector<double> audit_out;
  std::vector<double> split_out;
  std::vector<double> loose_out;
  TitanPBE::Tabi::TreecodeStats audit_stats =
      apply_fixture(panels, unknowns, executor, audit_tree, audit_out);
  TitanPBE::Tabi::TreecodeStats split_stats =
      apply_fixture(panels, unknowns, executor, split_tree, split_out);
  TitanPBE::Tabi::TreecodeStats loose_stats =
      apply_fixture(panels, unknowns, executor, loose_tree, loose_out);

  const double audit_rel = relative_rms(direct, audit_out);
  const double split_rel = relative_rms(direct, split_out);
  const double loose_rel = relative_rms(direct, loose_out);

  require_true("audit estimator sampled accepted far nodes",
               audit_stats.convergence_estimator_samples > 0);
  require_true("audit estimator reports unstable far nodes",
               audit_stats.convergence_estimator_unstable_far_count > 0);
  require_true("audit estimator records tail ratio",
               audit_stats.convergence_estimator_tail_ratio_max > 0.0);
  require_true("split estimator rejects unstable far nodes",
               split_stats.convergence_estimator_split_rejects > 0);
  require_true("split estimator records rejected panels",
               split_stats.convergence_estimator_split_reject_panels > 0);
  require_true("split increases leaf work",
               split_stats.leaf_panel_interactions >
                   audit_stats.leaf_panel_interactions);
  require_true("split reduces accepted far work",
               split_stats.far_node_interactions <=
                   audit_stats.far_node_interactions);
  require_true("split does not worsen tree/direct vector error",
               split_rel <= 1.05 * std::max(audit_rel, 1.0e-14));
  require_true("loose estimator leaves benign far path available",
               loose_stats.far_node_interactions > 0);
  require_true("loose estimator is less aggressive than hybrid split",
               loose_stats.convergence_estimator_split_rejects <
                   split_stats.convergence_estimator_split_rejects);
  if (audit_stats.far_expansion_audit_pairs > 0 &&
      split_stats.far_expansion_audit_pairs > 0) {
    require_true("split lowers audited far expansion error",
                 split_stats.far_expansion_rel_error_max <=
                     audit_stats.far_expansion_rel_error_max);
  }

  std::cout << "PBE TABI convergence estimator smoke PASS panels="
            << panels.size() << " audit_rel=" << audit_rel
            << " split_rel=" << split_rel << " loose_rel=" << loose_rel
            << " audit_samples="
            << audit_stats.convergence_estimator_samples
            << " audit_unstable="
            << audit_stats.convergence_estimator_unstable_far_count
            << " split_rejects="
            << split_stats.convergence_estimator_split_rejects
            << " loose_rejects="
            << loose_stats.convergence_estimator_split_rejects
            << " tail_ratio_max="
            << audit_stats.convergence_estimator_tail_ratio_max
            << " estimated_error_max="
            << audit_stats.convergence_estimator_error_max
            << " audit_far=" << audit_stats.far_node_interactions
            << " split_far=" << split_stats.far_node_interactions
            << " audit_leaf=" << audit_stats.leaf_panel_interactions
            << " split_leaf=" << split_stats.leaf_panel_interactions
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
    throw "Failed to compile TABI convergence estimator smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI convergence estimator smoke fixture failed."
}
