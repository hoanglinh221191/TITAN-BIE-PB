param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_channel_expansion_audit_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_channel_expansion_audit_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_channel_expansion_audit_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside channel expansion audit smoke root: $full"
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
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kK1 = 0;
constexpr std::size_t kK2 = 1;
constexpr std::size_t kK3 = 2;
constexpr std::size_t kK4 = 3;

void require_true(const std::string &label, bool value) {
  if (!value)
    throw std::runtime_error(label);
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
  p.area = 0.18 + 0.01 * static_cast<double>(index % 17);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(160);
  std::size_t index = 0;
  for (int i = 0; i < 36; ++i) {
    const double a = 0.17453292519943295 * static_cast<double>(i);
    panels.push_back(make_panel(-34.0 + 1.2 * std::cos(a),
                                1.2 * std::sin(a),
                                0.08 * std::sin(3.0 * a), 0.1, 0.2, 1.0,
                                index++));
  }
  for (int ix = -6; ix <= 6; ++ix) {
    if (ix == 0)
      continue;
    for (int iy = -4; iy <= 4; ++iy) {
      const double x = 24.0 + 0.95 * static_cast<double>(ix);
      const double y = 0.72 * static_cast<double>(iy);
      const double z = 0.05 * static_cast<double>((ix + 2 * iy) % 7);
      panels.push_back(make_panel(x, y, z, 0.2, -0.1, 1.0, index++));
    }
  }
  for (int i = 0; i < 24; ++i) {
    const double a = 0.2617993877991494 * static_cast<double>(i);
    panels.push_back(make_panel(48.0 + 2.2 * std::cos(a),
                                -5.0 + 2.0 * std::sin(a),
                                0.08 * std::cos(2.0 * a), -0.1, 0.3, 1.0,
                                index++));
  }
  return panels;
}

std::vector<double> make_unknowns(std::size_t n, bool k2_dominated) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    const double wave = 1.0 + 0.2 * std::sin(0.19 * static_cast<double>(i));
    if (k2_dominated) {
      unknowns[i] = 25.0 * sign * wave;
      unknowns[n + i] = 1.0e-4 * sign * wave;
    } else {
      unknowns[i] = 1.0e-4 * sign * wave;
      unknowns[n + i] = 25.0 * sign * wave;
    }
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.70, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "sampled", 512, 4096,
          "translated", "naive", "sampled", "audit", 0.0, 0.0, "hybrid",
          1.0e4, 1.0e-3, 0.25, "consecutive-aspect", "off", "on", 1, 64,
          4, "channel-tail-eps", "summary");
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  return tree;
}

TitanPBE::Tabi::TreecodeStats run_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor) {
  TitanPBE::Tabi::TreeWorkspace tree = make_tree(panels, executor);
  std::vector<double> out;
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.10, 0.0, 0.0, true, tree, "centroid", "off",
      executor, out, &stats);
  return stats;
}

} // namespace

int main() {
  const std::vector<TitanPBE::Tabi::Panel> panels = build_fixture();
  TitanExecutor executor(2);

  const TitanPBE::Tabi::TreecodeStats k2_stats =
      run_tree(panels, make_unknowns(panels.size(), true), executor);
  const TitanPBE::Tabi::TreecodeStats phin_stats =
      run_tree(panels, make_unknowns(panels.size(), false), executor);

  require_true("channel audit enabled",
               k2_stats.expansion_audit_channel == "summary");
  require_true("channel audit sampled pairs",
               k2_stats.far_expansion_audit_pairs > 0);
  require_true("overall p99 reported",
               k2_stats.far_expansion_rel_error_p99 >= 0.0);
  require_true("worst channel reported",
               k2_stats.far_expansion_worst_channel != "none");
  require_true("phi dominated case reports K2 channel",
               k2_stats.far_expansion_channel_rel_error_max[kK2] > 0.0);
  require_true("phin dominated case samples channels",
               phin_stats.far_expansion_audit_pairs > 0);
  require_true("phin dominated case reports phin-side channel",
               std::max(phin_stats.far_expansion_channel_rel_error_max[kK1],
                        phin_stats.far_expansion_channel_rel_error_max[kK3]) >
                   0.0);

  std::cout << "PBE TABI channel expansion audit smoke PASS panels="
            << panels.size() << " pairs="
            << k2_stats.far_expansion_audit_pairs << " p99="
            << k2_stats.far_expansion_rel_error_p99 << " worst="
            << k2_stats.far_expansion_worst_channel << " k2_max="
            << k2_stats.far_expansion_channel_rel_error_max[kK2]
            << " phin_k1_max="
            << phin_stats.far_expansion_channel_rel_error_max[kK1]
            << " phin_k3_max="
            << phin_stats.far_expansion_channel_rel_error_max[kK3] << "\n";
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
    throw "Failed to compile TABI channel expansion audit smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI channel expansion audit smoke fixture failed."
}
