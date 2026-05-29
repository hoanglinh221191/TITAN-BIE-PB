param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_risk_scaled_surface_margin_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_risk_scaled_surface_margin_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_risk_scaled_surface_margin_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside risk-scaled surface-margin smoke root: $full"
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

TitanPBE::Tabi::Panel panel(double x, double y, double z, std::size_t index) {
  TitanPBE::Tabi::Panel p;
  p.x = x;
  p.y = y;
  p.z = z;
  const double nx = 0.20 + 0.03 * static_cast<double>(index % 7);
  const double ny = ((index % 3) == 0 ? -0.45 : 0.55);
  const double nz = 0.80;
  const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  p.nx = nx / norm;
  p.ny = ny / norm;
  p.nz = nz / norm;
  p.area = 0.12 + 0.01 * static_cast<double>(index % 13);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture(double separation) {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(128);
  std::size_t index = 0;
  for (int side = 0; side < 2; ++side) {
    const double cx = side == 0 ? 0.0 : separation;
    for (int ix = -4; ix <= 4; ++ix) {
      for (int iy = -3; iy <= 3; ++iy) {
        if (std::abs(ix) <= 1 && std::abs(iy) <= 1)
          continue;
        const double x = cx + 0.82 * static_cast<double>(ix);
        const double y = 0.72 * static_cast<double>(iy);
        const double z = 0.10 * static_cast<double>((ix - 2 * iy + 19) % 7);
        panels.push_back(panel(x, y, z, index++));
      }
    }
  }
  return panels;
}

std::vector<double> unknowns(std::size_t n) {
  std::vector<double> values(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    values[i] = 1.2 * sign * (1.0 + 0.11 * std::sin(0.23 * i));
    values[n + i] = 0.7 * sign * (1.0 + 0.13 * std::cos(0.31 * i));
  }
  return values;
}

TitanPBE::Tabi::TreecodeStats run_case(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &values, TitanExecutor &executor,
    double hard_aspect, double hard_rho, double hard_factor,
    double soft_aspect, double soft_rho, double soft_factor) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.82, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "sampled", 256, 4096,
          "translated", "naive", "sampled", "split", 0.0, 0.0, "hybrid",
          1.0e12, 1.0e12, 0.25, "consecutive-aspect", "off", "on", 1, 64, 4,
          "channel-tail-eps", "summary", "polynomial-eps", "accepted-only",
          "global", "summary", 0.0, "surface-margin", hard_factor,
          "risk-scaled", hard_aspect, hard_rho, 1.0e-7, "auto", "off",
          "summary", "baseline", "off", "off", "scan", 32,
          "full-traversal", "off", "off", "inline-select", 2048, "off",
          "off", soft_aspect, soft_rho, soft_factor);
  tree.taylor_order = 4;
  tree.derivative_mode = "analytic";
  std::vector<double> out;
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, values, 80.0, 0.10, 0.0, 0.0, true, tree, "centroid", "off",
      executor, out, &stats);
  return stats;
}

std::size_t accepted_misses(const TitanPBE::Tabi::TreecodeStats &stats) {
  std::size_t total = stats.convergence_estimator_high_error_missed;
  for (std::size_t value : stats.convergence_channel_high_error_missed)
    total += value;
  return total;
}

} // namespace

int main() {
  TitanExecutor executor(2);

  const std::vector<TitanPBE::Tabi::Panel> high_risk = build_fixture(8.0);
  const std::vector<double> high_values = unknowns(high_risk.size());
  const TitanPBE::Tabi::TreecodeStats hard =
      run_case(high_risk, high_values, executor, 3.0, 0.55, 0.25, 2.0, 0.45,
               0.10);
  require_true("risk-scaled policy is reported",
               hard.mac_surface_margin_policy == "risk-scaled");
  require_true("hard-risk path applies full margin",
               hard.mac_surface_margin_hard_applied > 0);
  require_true("hard-risk path rejects some candidates",
               hard.mac_surface_margin_hard_rejects > 0 ||
                   hard.mac_surface_margin_rejects > 0);
  require_true("hard-risk accepted misses stay clean",
               accepted_misses(hard) == 0);

  const TitanPBE::Tabi::TreecodeStats soft =
      run_case(high_risk, high_values, executor, 50.0, 0.95, 0.25, 2.0, 0.45,
               0.10);
  require_true("soft-risk path applies soft margin",
               soft.mac_surface_margin_soft_applied > 0);
  require_true("soft-risk keeps accepted misses clean",
               accepted_misses(soft) == 0);

  const std::vector<TitanPBE::Tabi::Panel> benign = build_fixture(25.0);
  const std::vector<double> benign_values = unknowns(benign.size());
  const TitanPBE::Tabi::TreecodeStats benign_stats =
      run_case(benign, benign_values, executor, 50.0, 2.0, 0.25, 50.0, 2.0,
               0.10);
  require_true("benign fixture keeps accepted far work",
               benign_stats.far_node_interactions > 0);
  require_true("benign fixture skips margin",
               benign_stats.mac_surface_margin_skipped > 0);
  require_true("benign fixture avoids margin flood",
               benign_stats.mac_surface_margin_skipped >=
                   benign_stats.mac_surface_margin_hard_applied +
                       benign_stats.mac_surface_margin_soft_applied);

  std::cout << "PBE TABI risk-scaled surface-margin smoke PASS hard="
            << hard.mac_surface_margin_hard_applied
            << " soft=" << soft.mac_surface_margin_soft_applied
            << " skipped=" << benign_stats.mac_surface_margin_skipped
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
    (Join-Path $RepoRoot 'src\tabi\tabi_near.cpp'),
    (Get-PbeTabiTreecodeSources -RepoRoot $RepoRoot),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI risk-scaled surface-margin smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI risk-scaled surface-margin smoke fixture failed."
}
