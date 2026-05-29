param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_absolute_floor_miss_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_absolute_floor_miss_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_absolute_floor_miss_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside absolute-floor smoke root: $full"
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

TitanPBE::Tabi::Panel make_panel(double x, double y, double z,
                                 std::size_t index) {
  TitanPBE::Tabi::Panel p;
  p.x = x;
  p.y = y;
  p.z = z;
  const double nx = 0.25 + 0.04 * static_cast<double>(index % 7);
  const double ny = ((index % 2) == 0 ? 0.55 : -0.50);
  const double nz = 0.70;
  const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  p.nx = nx / norm;
  p.ny = ny / norm;
  p.nz = nz / norm;
  p.area = 0.12 + 0.013 * static_cast<double>(index % 11);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(112);
  std::size_t index = 0;
  for (int side = 0; side < 2; ++side) {
    const double cx = side == 0 ? 0.0 : 7.2;
    for (int ix = -3; ix <= 3; ++ix) {
      for (int iy = -3; iy <= 3; ++iy) {
        if (std::abs(ix) == 0 && std::abs(iy) == 0)
          continue;
        panels.push_back(make_panel(cx + 0.74 * static_cast<double>(ix),
                                    0.68 * static_cast<double>(iy),
                                    0.06 * static_cast<double>((ix + iy + 9) % 5),
                                    index++));
      }
    }
  }
  return panels;
}

std::vector<double> make_unknowns(std::size_t n, double scale) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    const double wave = 1.0 + 0.25 * std::sin(0.31 * static_cast<double>(i));
    unknowns[i] = scale * sign * wave;
    unknowns[n + i] = 0.45 * scale * sign *
                      (1.0 + 0.20 * std::cos(0.17 * static_cast<double>(i)));
  }
  return unknowns;
}

TitanPBE::Tabi::TreecodeStats run_case(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    TitanExecutor &executor, const std::string &safety, double scale,
    double absolute_floor) {
  const std::vector<double> unknowns = make_unknowns(panels.size(), scale);
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.72, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "sampled", 512, 4096,
          "translated", "naive", "sampled", safety, 0.0, 0.0, "hybrid",
          1.0e4, 1.0e-12, 0.25, "consecutive-aspect", "off", "on", 1, 64,
          4, "channel-tail-eps", "summary", "polynomial-eps",
          "accepted-only", "global", "sampled-detail", 0.0, "center", 0.25,
          "all", 3.0, 0.55, absolute_floor);
  tree.taylor_order = 4;
  tree.derivative_mode = "analytic";
  std::vector<double> out;
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.10, 0.0, 0.0, true, tree, "centroid", "off",
      executor, out, &stats);
  return stats;
}

std::size_t channel_sum(
    const std::array<std::size_t, TitanPBE::Tabi::kJcpExpansionChannelCount>
        &values) {
  std::size_t sum = 0;
  for (std::size_t value : values)
    sum += value;
  return sum;
}

} // namespace

int main() {
  TitanExecutor executor(2);
  const std::vector<TitanPBE::Tabi::Panel> panels = build_fixture();

  const TitanPBE::Tabi::TreecodeStats phantom =
      run_case(panels, executor, "audit", 1.0e-9, 1.0e-4);
  require_true("absolute floor is reported",
               std::abs(phantom.accepted_miss_absolute_floor - 1.0e-4) <
                   1.0e-18);
  require_true("phantom relative high-error counted",
               phantom.convergence_estimator_phantom_high_error_count > 0 ||
                   channel_sum(phantom.phantom_high_error_channel_count) > 0);
  require_true("phantom accepted high-error does not increment missed",
               phantom.convergence_estimator_high_error_missed == 0 &&
                   channel_sum(phantom.convergence_channel_high_error_missed) ==
                       0);

  const TitanPBE::Tabi::TreecodeStats real =
      run_case(panels, executor, "audit", 25.0, 1.0e-7);
  require_true("real high-error increments accepted missed",
               real.convergence_estimator_high_error_missed > 0 ||
                   channel_sum(real.convergence_channel_high_error_missed) > 0);

  const TitanPBE::Tabi::TreecodeStats split_real =
      run_case(panels, executor, "split", 25.0, 1.0e-7);
  require_true("split real high-error is caught by split",
               split_real.convergence_estimator_high_error_caught_by_split >
                       0 ||
                   channel_sum(split_real
                                   .convergence_channel_high_error_caught_by_split) >
                       0);

  const TitanPBE::Tabi::TreecodeStats split_phantom =
      run_case(panels, executor, "split", 1.0e-9, 1.0e-4);
  require_true("split phantom high-error is not caught_by_split",
               split_phantom.convergence_estimator_high_error_caught_by_split ==
                       0 &&
                   channel_sum(split_phantom
                                   .convergence_channel_high_error_caught_by_split) ==
                       0);
  require_true("split phantom still reports phantom counter",
               split_phantom.convergence_estimator_phantom_high_error_count >
                       0 ||
                   channel_sum(split_phantom.phantom_high_error_channel_count) >
                       0);

  std::cout << "PBE TABI absolute-floor miss smoke PASS panels="
            << panels.size()
            << " phantom="
            << phantom.convergence_estimator_phantom_high_error_count
            << " real_missed="
            << real.convergence_estimator_high_error_missed
            << " split_caught="
            << split_real.convergence_estimator_high_error_caught_by_split
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
    throw "Failed to compile TABI absolute-floor miss smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI absolute-floor miss smoke fixture failed."
}
