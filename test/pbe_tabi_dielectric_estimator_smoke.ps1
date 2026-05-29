param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_dielectric_estimator_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_dielectric_estimator_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_dielectric_estimator_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside dielectric estimator smoke root: $full"
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
  p.area = 0.25 + 0.02 * static_cast<double>(index % 13);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(120);
  std::size_t index = 0;
  for (int i = 0; i < 32; ++i) {
    const double angle = 0.19634954084936207 * static_cast<double>(i);
    panels.push_back(make_panel(-28.0 + 1.0 * std::cos(angle),
                                1.0 * std::sin(angle),
                                0.10 * std::sin(3.0 * angle), 0.0, 0.0, 1.0,
                                index++));
  }
  for (int ix = -5; ix <= 5; ++ix) {
    if (ix == 0)
      continue;
    for (int iy = -3; iy <= 3; ++iy) {
      const double x = 22.0 + 0.95 * static_cast<double>(ix);
      const double y = 0.82 * static_cast<double>(iy);
      const double z = 0.06 * static_cast<double>((ix + iy) % 5);
      panels.push_back(make_panel(x, y, z, 0.1, 0.2, 1.0, index++));
    }
  }
  for (int i = 0; i < 18; ++i) {
    const double angle = 0.3490658503988659 * static_cast<double>(i);
    panels.push_back(make_panel(42.0 + 2.0 * std::cos(angle),
                                -5.0 + 2.0 * std::sin(angle),
                                0.08 * std::cos(2.0 * angle), -0.2, 0.1, 1.0,
                                index++));
  }
  return panels;
}

std::vector<double> make_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    const double burst = (i % 7 == 0) ? 5.0 : 1.0;
    unknowns[i] =
        sign * burst * (0.8 + 0.20 * std::sin(0.29 * static_cast<double>(i)));
    unknowns[n + i] =
        -sign * (1.3 + 0.25 * std::cos(0.17 * static_cast<double>(i)));
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &dielectric_amp, double error_cap) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.70, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "off", 512, 2048, "translated",
          "naive", "sampled", "split", 0.0, 0.0, "hybrid", 1.0e12,
          error_cap, 0.25, "stage3j", "off", dielectric_amp, 1, 64, 4);
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  return tree;
}

TitanPBE::Tabi::TreecodeStats run_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    TitanPBE::Tabi::TreeWorkspace &tree, double eps_ratio) {
  std::vector<double> out;
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, 0.10, 0.0, 0.0, true, tree, "centroid",
      "off", executor, out, &stats);
  return stats;
}

} // namespace

int main() {
  const std::vector<TitanPBE::Tabi::Panel> panels = build_fixture();
  const std::vector<double> unknowns = make_unknowns(panels.size());
  TitanExecutor executor(2);

  TitanPBE::Tabi::TreeWorkspace off_probe =
      make_tree(panels, executor, "off", 1.0);
  TitanPBE::Tabi::TreeWorkspace on_probe =
      make_tree(panels, executor, "on", 1.0);
  const TitanPBE::Tabi::TreecodeStats off_stats =
      run_tree(panels, unknowns, executor, off_probe, 80.0);
  const TitanPBE::Tabi::TreecodeStats on_stats =
      run_tree(panels, unknowns, executor, on_probe, 80.0);

  require_true("off estimator sampled nodes",
               off_stats.convergence_estimator_samples > 0);
  require_true("dielectric amp reports baseline 30",
               off_stats.convergence_dielectric_amp_max >= 30.0 &&
                   off_stats.convergence_dielectric_amp_max < 31.0);
  require_true("dielectric amp on reaches eps_ratio",
               on_stats.convergence_dielectric_amp_max >= 80.0);
  require_true("effective amp follows dielectric amp",
               on_stats.convergence_effective_jcp_amp_max >=
                   on_stats.convergence_dielectric_amp_max);
  require_true("dielectric amp raises estimated error",
               on_stats.convergence_estimator_error_max >
                   off_stats.convergence_estimator_error_max);

  const double dynamic_cap =
      0.5 * (off_stats.convergence_estimator_error_max +
             on_stats.convergence_estimator_error_max);
  TitanPBE::Tabi::TreeWorkspace off_split =
      make_tree(panels, executor, "off", dynamic_cap);
  TitanPBE::Tabi::TreeWorkspace on_split =
      make_tree(panels, executor, "on", dynamic_cap);
  const TitanPBE::Tabi::TreecodeStats off_split_stats =
      run_tree(panels, unknowns, executor, off_split, 80.0);
  const TitanPBE::Tabi::TreecodeStats on_split_stats =
      run_tree(panels, unknowns, executor, on_split, 80.0);
  require_true("dielectric amp increases split rejects",
               on_split_stats.convergence_estimator_split_rejects >
                   off_split_stats.convergence_estimator_split_rejects);
  require_true("dielectric amp increases unstable nodes",
               on_split_stats.convergence_estimator_unstable_far_count >
                   off_split_stats.convergence_estimator_unstable_far_count);

  TitanPBE::Tabi::TreeWorkspace low_eps_tree =
      make_tree(panels, executor, "on", dynamic_cap);
  const TitanPBE::Tabi::TreecodeStats low_eps =
      run_tree(panels, unknowns, executor, low_eps_tree, 2.0);
  require_true("low eps keeps conservative baseline amp",
               low_eps.convergence_dielectric_amp_max >= 30.0 &&
                   low_eps.convergence_dielectric_amp_max < 31.0);
  require_true("low eps does not reject above high eps dielectric path",
               low_eps.convergence_estimator_split_rejects <=
                   on_split_stats.convergence_estimator_split_rejects);

  std::cout << "PBE TABI dielectric estimator smoke PASS panels="
            << panels.size() << " off_amp="
            << off_stats.convergence_dielectric_amp_max << " on_amp="
            << on_stats.convergence_dielectric_amp_max << " off_err="
            << off_stats.convergence_estimator_error_max << " on_err="
            << on_stats.convergence_estimator_error_max << " cap="
            << dynamic_cap << " off_rejects="
            << off_split_stats.convergence_estimator_split_rejects
            << " on_rejects="
            << on_split_stats.convergence_estimator_split_rejects
            << " low_eps_amp=" << low_eps.convergence_dielectric_amp_max
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
    throw "Failed to compile TABI dielectric estimator smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI dielectric estimator smoke fixture failed."
}
