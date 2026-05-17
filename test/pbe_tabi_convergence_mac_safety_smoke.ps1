param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_convergence_mac_safety_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_convergence_mac_safety_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_convergence_mac_safety_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside convergence MAC safety smoke root: $full"
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
  p.area = 0.55 + 0.02 * static_cast<double>(index % 13);
  p.diagonal = 0.0;
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_high_aspect_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(72);
  std::size_t index = 0;

  for (int i = 0; i < 48; ++i) {
    const double x = -22.0 + 0.92 * static_cast<double>(i);
    const double y = 0.018 * static_cast<double>((i % 3) - 1);
    const double z = 0.006 * static_cast<double>((i % 5) - 2);
    panels.push_back(make_panel(x, y, z, 0.0, 0.15, 1.0, index++));
  }
  for (int i = 0; i < 24; ++i) {
    const double x = 46.0 + 0.20 * static_cast<double>(i % 6);
    const double y = 9.0 + 0.18 * static_cast<double>(i / 6);
    const double z = 0.35 + 0.025 * static_cast<double>(i % 4);
    panels.push_back(make_panel(x, y, z, -0.25, 0.05, 1.0, index++));
  }
  return panels;
}

std::vector<double> make_unknowns(std::size_t n) {
  std::vector<double> unknowns(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 5 == 0 || i % 11 == 0) ? -1.0 : 1.0;
    unknowns[i] =
        sign * (0.9 + 0.35 * std::sin(0.31 * static_cast<double>(i)) +
                0.20 * std::cos(0.17 * static_cast<double>(i)));
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &convergence_safety) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.70, "exact", "scalar", "strict",
          "adaptive", "onthefly", "tight-aabb", "topdown-aabb",
          "physical-gap", 1.5, "off", "area-centroid", "split", "sampled",
          2048, 8192, "translated", "naive", "sampled",
          convergence_safety, 0.25, 3.0);
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
      panels, unknowns, 80.0, 0.10, 0.0, 0.0, false, tree, "centroid",
      "off", executor, out, &stats);
  return stats;
}

} // namespace

int main() {
  const std::vector<TitanPBE::Tabi::Panel> panels =
      build_high_aspect_fixture();
  const std::vector<double> unknowns = make_unknowns(panels.size());
  TitanExecutor executor(2);

  std::vector<double> direct;
  TitanPBE::Tabi::TreecodeStats direct_stats;
  TitanPBE::Tabi::TreeWorkspace empty_tree;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.10, 0.0, 0.0, false, empty_tree, "centroid",
      "off", executor, direct, &direct_stats);
  require_true("direct fallback mode",
               direct_stats.matvec_mode == "jcp2013_direct_fallback");

  TitanPBE::Tabi::TreeWorkspace off_tree = make_tree(panels, executor, "off");
  TitanPBE::Tabi::TreeWorkspace audit_tree =
      make_tree(panels, executor, "audit");
  TitanPBE::Tabi::TreeWorkspace split_tree =
      make_tree(panels, executor, "split");
  require_true("audit convergence safety propagated",
               audit_tree.convergence_safety == "audit");
  require_true("split convergence safety propagated",
               split_tree.convergence_safety == "split");
  require_true("rho cap propagated", split_tree.convergence_rho_cap == 0.25);
  require_true("aspect cap propagated", split_tree.mac_aspect_cap == 3.0);

  std::vector<double> off_out;
  std::vector<double> audit_out;
  std::vector<double> split_out;
  TitanPBE::Tabi::TreecodeStats off_stats =
      apply_fixture(panels, unknowns, executor, off_tree, off_out);
  TitanPBE::Tabi::TreecodeStats audit_stats =
      apply_fixture(panels, unknowns, executor, audit_tree, audit_out);
  TitanPBE::Tabi::TreecodeStats split_stats =
      apply_fixture(panels, unknowns, executor, split_tree, split_out);

  const std::size_t audit_rejects =
      audit_stats.mac_convergence_rho_rejects +
      audit_stats.mac_convergence_aspect_rejects;
  const std::size_t split_rejects =
      split_stats.mac_convergence_rho_rejects +
      split_stats.mac_convergence_aspect_rejects;
  const double off_rel = relative_rms(direct, off_out);
  const double audit_rel = relative_rms(direct, audit_out);
  const double split_rel = relative_rms(direct, split_out);

  require_true("baseline uses far work", off_stats.far_node_interactions > 0);
  require_true("audit reports would-reject convergence MAC candidates",
               audit_rejects > 0);
  require_true("audit reports would-reject panels",
               audit_stats.mac_convergence_rho_reject_panels +
                       audit_stats.mac_convergence_aspect_reject_panels >
                   0);
  require_true("split rejects convergence MAC candidates", split_rejects > 0);
  require_true("split increases leaf work",
               split_stats.leaf_panel_interactions >
                   audit_stats.leaf_panel_interactions);
  require_true("split reduces or preserves far work",
               split_stats.far_node_interactions <=
                   audit_stats.far_node_interactions);
  require_true("split does not worsen tree/direct vector error",
               split_rel <= 1.05 * std::max(off_rel, audit_rel));
  if (audit_stats.far_expansion_audit_pairs > 0 &&
      split_stats.far_expansion_audit_pairs > 0) {
    require_true("split lowers audited far expansion error",
                 split_stats.far_expansion_rel_error_max <=
                     audit_stats.far_expansion_rel_error_max);
  }

  std::cout << "PBE TABI convergence MAC safety smoke PASS panels="
            << panels.size() << " off_rel=" << off_rel
            << " audit_rel=" << audit_rel << " split_rel=" << split_rel
            << " audit_rho_rejects="
            << audit_stats.mac_convergence_rho_rejects
            << " audit_aspect_rejects="
            << audit_stats.mac_convergence_aspect_rejects
            << " split_rho_rejects="
            << split_stats.mac_convergence_rho_rejects
            << " split_aspect_rejects="
            << split_stats.mac_convergence_aspect_rejects
            << " audit_far=" << audit_stats.far_node_interactions
            << " split_far=" << split_stats.far_node_interactions
            << " audit_leaf=" << audit_stats.leaf_panel_interactions
            << " split_leaf=" << split_stats.leaf_panel_interactions
            << " audit_expansion_rel_max="
            << audit_stats.far_expansion_rel_error_max
            << " split_expansion_rel_max="
            << split_stats.far_expansion_rel_error_max << "\n";
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
    (Join-Path $RepoRoot 'titan_pbe\tabi_treecode.cpp'),
    (Join-Path $RepoRoot 'titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI convergence MAC safety smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI convergence MAC safety smoke fixture failed."
}
