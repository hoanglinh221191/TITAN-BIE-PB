param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_far_expansion_audit_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_far_expansion_audit_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_far_expansion_audit_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside far-expansion smoke root: $full"
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
  p.area = 0.7 + 0.03 * static_cast<double>(index % 7);
  p.diagonal = 0.0;
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_hollow_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(64);
  std::size_t index = 0;
  auto add_patch = [&](double ox, double oy, double oz, int nx, int ny) {
    for (int ix = 0; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy) {
        panels.push_back(make_panel(
            ox + 0.28 * static_cast<double>(ix),
            oy + 0.31 * static_cast<double>(iy),
            oz + 0.05 * static_cast<double>((ix + 2 * iy) % 3), index++));
      }
    }
  };
  add_patch(-58.0, -45.0, -18.0, 4, 4);
  add_patch(16.0, 5.0, 0.0, 4, 4);
  add_patch(16.4, 29.0, 0.0, 4, 4);
  add_patch(21.0, 17.0, 7.0, 4, 4);
  return panels;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &fragmentation_safety) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 4, 2.0, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", fragmentation_safety, "sampled", 2048,
          8192);
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  return tree;
}

std::vector<double> make_unknowns(std::size_t n) {
  std::vector<double> unknowns(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 4 == 0 || i % 5 == 0) ? -1.0 : 1.0;
    unknowns[i] = sign * (0.8 + 0.45 * std::sin(0.41 * static_cast<double>(i)) +
                          0.25 * std::cos(0.23 * static_cast<double>(i)));
  }
  return unknowns;
}

} // namespace

int main() {
  const std::vector<TitanPBE::Tabi::Panel> panels = build_hollow_fixture();
  const std::vector<double> unknowns = make_unknowns(panels.size());
  TitanExecutor executor(2);

  std::vector<double> direct;
  TitanPBE::Tabi::TreecodeStats direct_stats;
  TitanPBE::Tabi::TreeWorkspace empty_tree;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.12, 0.0, 0.0, false, empty_tree, "centroid",
      "off", executor, direct, &direct_stats);
  require_true("direct fallback mode",
               direct_stats.matvec_mode == "jcp2013_direct_fallback");

  TitanPBE::Tabi::TreeWorkspace audit_tree = make_tree(panels, executor, "audit");
  TitanPBE::Tabi::TreeWorkspace split_tree = make_tree(panels, executor, "split");
  require_true("audit enabled", audit_tree.expansion_audit == "sampled");

  std::vector<double> audit_out;
  TitanPBE::Tabi::TreecodeStats audit_stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.12, 0.0, 0.0, false, audit_tree, "centroid",
      "off", executor, audit_out, &audit_stats);

  std::vector<double> split_out;
  TitanPBE::Tabi::TreecodeStats split_stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.12, 0.0, 0.0, false, split_tree, "centroid",
      "off", executor, split_out, &split_stats);

  const double audit_rel = relative_rms(direct, audit_out);
  const double split_rel = relative_rms(direct, split_out);
  require_true("audit sampled far expansion pairs",
               audit_stats.far_expansion_audit_pairs > 0);
  require_true("audit saw fragmented far candidates",
               audit_stats.fragmented_far_count > 0);
  require_true("Taylor far expansion error is visible",
               audit_stats.far_expansion_rel_error_max > 1.0e-5);
  require_true("split sharply reduces audited Taylor error",
               split_stats.far_expansion_audit_pairs == 0 ||
                   audit_stats.far_expansion_rel_error_max >
                       10.0 * split_stats.far_expansion_rel_error_max);
  require_true("split rejected fragmented far candidates",
               split_stats.fragmented_split_rejects > 0);
  require_true("split improves tree/direct vector error",
               split_rel < 0.8 * audit_rel);
  require_true("split lowers audited Taylor error or removes the bad node",
               split_stats.far_expansion_audit_pairs == 0 ||
                   split_stats.far_expansion_rel_error_max <
                       audit_stats.far_expansion_rel_error_max);

  std::cout << "PBE TABI far-expansion audit smoke PASS panels="
            << panels.size() << " audit_rel=" << audit_rel
            << " split_rel=" << split_rel
            << " audit_pairs=" << audit_stats.far_expansion_audit_pairs
            << " audit_panels=" << audit_stats.far_expansion_audit_panels
            << " expansion_rel_max="
            << audit_stats.far_expansion_rel_error_max
            << " expansion_rel_mean="
            << audit_stats.far_expansion_rel_error_mean
            << " split_rejects=" << split_stats.fragmented_split_rejects
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
    throw "Failed to compile TABI far-expansion audit smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI far-expansion audit smoke fixture failed."
}
