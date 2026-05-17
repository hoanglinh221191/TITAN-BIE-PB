param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_jcp_moment_construction_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_jcp_moment_construction_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_jcp_moment_construction_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside moment-construction smoke root: $full"
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
                                 std::size_t index) {
  TitanPBE::Tabi::Panel p;
  p.x = x;
  p.y = y;
  p.z = z;
  const double norm = std::sqrt(std::max(1.0e-12, x * x + y * y + z * z));
  p.nx = x / norm;
  p.ny = y / norm;
  p.nz = z / norm;
  p.area = 0.6 + 0.05 * static_cast<double>(index % 11);
  p.diagonal = 0.0;
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_hollow_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(96);
  std::size_t index = 0;
  auto add_patch = [&](double ox, double oy, double oz, int nx, int ny,
                       double sx, double sy) {
    for (int ix = 0; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy) {
        panels.push_back(make_panel(
            ox + sx * static_cast<double>(ix),
            oy + sy * static_cast<double>(iy),
            oz + 0.07 * static_cast<double>((ix + 3 * iy) % 5), index++));
      }
    }
  };
  add_patch(-65.0, -40.0, -8.0, 4, 4, 0.35, 0.31);
  add_patch(14.0, 7.0, 0.0, 4, 4, 0.33, 0.29);
  add_patch(14.5, 35.0, 1.0, 4, 4, 0.30, 0.34);
  add_patch(38.0, 36.0, 4.0, 4, 4, 0.32, 0.30);
  return panels;
}

std::vector<double> make_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 3 == 0 || i % 7 == 0) ? -1.0 : 1.0;
    unknowns[i] = sign * (0.7 + 0.40 * std::sin(0.37 * static_cast<double>(i)));
    unknowns[n + i] =
        (i % 2 == 0 ? 1.0 : -1.0) *
        (0.03 + 0.015 * std::cos(0.19 * static_cast<double>(i)));
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &moment_construction, const std::string &moment_sum) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 4, 2.0, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "sampled", "area-centroid", "audit", "sampled", 2048, 8192,
          moment_construction, moment_sum, "sampled");
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  return tree;
}

} // namespace

int main() {
  const std::vector<TitanPBE::Tabi::Panel> panels = build_hollow_fixture();
  const std::vector<double> unknowns = make_unknowns(panels.size());
  TitanExecutor executor(2);

  const double eps_ratio = 80.0;
  const double lhs_phi = 0.5 * (1.0 + eps_ratio);
  const double lhs_phin = 0.5 * (1.0 + 1.0 / eps_ratio);

  std::vector<double> direct;
  TitanPBE::Tabi::TreecodeStats direct_stats;
  TitanPBE::Tabi::TreeWorkspace empty_tree;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, 0.12, lhs_phi, lhs_phin, true, empty_tree,
      "centroid", "off", executor, direct, &direct_stats);
  require_true("direct fallback mode",
               direct_stats.matvec_mode == "jcp2013_direct_fallback");

  std::vector<double> translated_out;
  TitanPBE::Tabi::TreecodeStats translated_stats;
  TitanPBE::Tabi::TreeWorkspace translated_tree =
      make_tree(panels, executor, "translated", "naive");
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, 0.12, lhs_phi, lhs_phin, true,
      translated_tree, "centroid", "off", executor, translated_out,
      &translated_stats);

  std::vector<double> direct_node_out;
  TitanPBE::Tabi::TreecodeStats direct_node_stats;
  TitanPBE::Tabi::TreeWorkspace direct_node_tree =
      make_tree(panels, executor, "direct-node", "kahan");
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, 0.12, lhs_phi, lhs_phin, true,
      direct_node_tree, "centroid", "off", executor, direct_node_out,
      &direct_node_stats);

  const double translated_rel = relative_rms(direct, translated_out);
  const double direct_node_rel = relative_rms(direct, direct_node_out);

  require_true("translated mode propagated",
               translated_stats.moment_construction == "translated");
  require_true("direct-node mode propagated",
               direct_node_stats.moment_construction == "direct-node");
  require_true("kahan mode propagated", direct_node_stats.moment_sum == "kahan");
  require_true("moment audit sampled nodes",
               translated_stats.moment_audit_nodes > 0 &&
                   direct_node_stats.moment_audit_nodes > 0);
  require_true("translation offsets recorded",
               translated_stats.moment_translation_offset_max > 0.0);
  require_true("cancellation diagnostics recorded",
               direct_node_stats.moment_cancellation_ratio_max >= 1.0);
  require_true("dynamic range diagnostics recorded",
               direct_node_stats.moment_dynamic_range_max >= 1.0);
  require_true("node aspect diagnostics recorded",
               direct_node_stats.node_aspect_ratio_max >= 1.0);
  require_true("convergence rho diagnostics recorded",
               direct_node_stats.far_convergence_rho_max > 0.0);
  require_true("far expansion audit sampled pairs",
               direct_node_stats.far_expansion_audit_pairs > 0);
  require_true("direct-node vector finite and bounded",
               std::isfinite(direct_node_rel) && direct_node_rel < 10.0);
  if (translated_stats.far_expansion_rel_error_max > 1.0e-8 &&
      direct_node_stats.far_expansion_audit_pairs > 0) {
    require_true("direct-node does not inflate audited expansion error",
                 direct_node_stats.far_expansion_rel_error_max <=
                     1.25 * translated_stats.far_expansion_rel_error_max);
  }

  std::cout << "PBE TABI JCP moment-construction smoke PASS panels="
            << panels.size() << " translated_rel=" << translated_rel
            << " direct_node_rel=" << direct_node_rel
            << " moment_offset_max="
            << translated_stats.moment_translation_offset_max
            << " cancellation_max="
            << direct_node_stats.moment_cancellation_ratio_max
            << " dynamic_range_max="
            << direct_node_stats.moment_dynamic_range_max
            << " rho_max=" << direct_node_stats.far_convergence_rho_max
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
    (Join-Path $RepoRoot 'titan_pbe\tabi_treecode.cpp'),
    (Join-Path $RepoRoot 'titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI JCP moment-construction smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI JCP moment-construction smoke fixture failed."
}
