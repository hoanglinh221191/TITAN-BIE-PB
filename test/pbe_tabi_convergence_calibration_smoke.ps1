param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_convergence_calibration_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_convergence_calibration_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_convergence_calibration_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside convergence calibration smoke root: $full"
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
  p.area = 0.28 + 0.01 * static_cast<double>(index % 11);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_calibration_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(112);
  std::size_t index = 0;

  for (int i = 0; i < 16; ++i) {
    const double angle = 0.39269908169872414 * static_cast<double>(i);
    panels.push_back(make_panel(-32.0 + 1.2 * std::cos(angle),
                                1.2 * std::sin(angle),
                                0.12 * std::sin(2.0 * angle), 0.0, 0.0, 1.0,
                                index++));
  }

  for (int ix = -4; ix <= 4; ++ix) {
    if (ix == 0)
      continue;
    for (int iy = -2; iy <= 2; ++iy) {
      const double x = 24.0 + 1.15 * static_cast<double>(ix);
      const double y = 0.95 * static_cast<double>(iy);
      const double z = 0.08 * static_cast<double>((ix + iy) % 3);
      panels.push_back(make_panel(x, y, z, 0.1, 0.2, 1.0, index++));
    }
  }

  for (int i = 0; i < 56; ++i) {
    const double x = 46.0 + 0.62 * static_cast<double>(i);
    const double y = (i % 2 == 0 ? -0.18 : 0.18);
    const double z = 0.04 * static_cast<double>((i % 5) - 2);
    panels.push_back(make_panel(x, y, z, -0.1, 0.3, 1.0, index++));
  }

  return panels;
}

std::vector<double> make_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double phase = static_cast<double>(i);
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    unknowns[2 * i] = sign * (1.0 + 0.30 * std::sin(0.31 * phase));
    unknowns[2 * i + 1] = -sign * (1.7 + 0.20 * std::cos(0.23 * phase));
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &safety, const std::string &estimator,
    const std::string &calibration, double threshold, double error_cap) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.70, "exact", "scalar", "strict",
          "adaptive", "onthefly", "tight-aabb", "topdown-aabb",
          "physical-gap", 1.5, "off", "area-centroid", "split", "sampled",
          2048, 8192, "translated", "naive", "sampled", safety, 0.0, 0.0,
          estimator, threshold, error_cap, 0.25, calibration);
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  return tree;
}

TitanPBE::Tabi::TreecodeStats run_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    TitanPBE::Tabi::TreeWorkspace &tree) {
  std::vector<double> out;
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.10, 0.0, 0.0, true, tree, "centroid", "off",
      executor, out, &stats);
  return stats;
}

} // namespace

int main() {
  const std::vector<TitanPBE::Tabi::Panel> panels =
      build_calibration_fixture();
  const std::vector<double> unknowns = make_unknowns(panels.size());
  TitanExecutor executor(2);

  TitanPBE::Tabi::TreeWorkspace stage3j_tree =
      make_tree(panels, executor, "split", "hybrid", "stage3j", 1.0e12,
                1.0e-3);
  TitanPBE::Tabi::TreeWorkspace consecutive_tree =
      make_tree(panels, executor, "split", "hybrid", "consecutive", 1.0e12,
                1.0e-3);
  TitanPBE::Tabi::TreeWorkspace aspect_tree =
      make_tree(panels, executor, "split", "hybrid", "aspect", 1.0e12,
                1.0e-3);
  TitanPBE::Tabi::TreeWorkspace combined_tree =
      make_tree(panels, executor, "split", "hybrid", "consecutive-aspect",
                1.0e12, 1.0e-3);
  TitanPBE::Tabi::TreeWorkspace benign_tree =
      make_tree(panels, executor, "split", "moment-tail",
                "consecutive-aspect", 1.0e12, 1.0);

  const TitanPBE::Tabi::TreecodeStats stage3j =
      run_tree(panels, unknowns, executor, stage3j_tree);
  const TitanPBE::Tabi::TreecodeStats consecutive =
      run_tree(panels, unknowns, executor, consecutive_tree);
  const TitanPBE::Tabi::TreecodeStats aspect =
      run_tree(panels, unknowns, executor, aspect_tree);
  const TitanPBE::Tabi::TreecodeStats combined =
      run_tree(panels, unknowns, executor, combined_tree);
  const TitanPBE::Tabi::TreecodeStats benign =
      run_tree(panels, unknowns, executor, benign_tree);

  require_true("stage3j sampled estimator nodes",
               stage3j.convergence_estimator_samples > 0);
  require_true("consecutive reports degree p-1 tail",
               consecutive.convergence_tail_norm_pm1_max > 0.0);
  require_true("consecutive keeps degree p tail stats",
               consecutive.convergence_tail_norm_p_max > 0.0);
  require_true("consecutive does not lower estimated error",
               consecutive.convergence_estimator_error_max + 1.0e-18 >=
                   stage3j.convergence_estimator_error_max);
  require_true("consecutive does not lower split rejects",
               consecutive.convergence_estimator_split_rejects >=
                   stage3j.convergence_estimator_split_rejects);

  require_true("aspect amplifies JCP estimator",
               aspect.convergence_effective_jcp_amp_max >
                   stage3j.convergence_effective_jcp_amp_max);
  require_true("aspect does not lower estimated error",
               aspect.convergence_estimator_error_max + 1.0e-18 >=
                   stage3j.convergence_estimator_error_max);
  require_true("combined keeps strongest effective amp",
               combined.convergence_effective_jcp_amp_max >=
                   aspect.convergence_effective_jcp_amp_max);
  require_true("combined keeps p-1 tail stats",
               combined.convergence_tail_norm_pm1_mean > 0.0);
  require_true("benign calibration leaves far path available",
               benign.far_node_interactions > 0);
  require_true("benign calibration does not reject all far nodes",
               benign.convergence_estimator_split_rejects <
                   std::max<std::size_t>(
                       combined.convergence_estimator_samples, 1));

  std::cout << "PBE TABI convergence calibration smoke PASS panels="
            << panels.size() << " stage3j_err="
            << stage3j.convergence_estimator_error_max
            << " consecutive_err="
            << consecutive.convergence_estimator_error_max
            << " aspect_amp=" << aspect.convergence_effective_jcp_amp_max
            << " combined_amp="
            << combined.convergence_effective_jcp_amp_max
            << " tail_p=" << combined.convergence_tail_norm_p_max
            << " tail_pm1=" << combined.convergence_tail_norm_pm1_max
            << " combined_rejects="
            << combined.convergence_estimator_split_rejects
            << " benign_far=" << benign.far_node_interactions << "\n";
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
    throw "Failed to compile TABI convergence calibration smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI convergence calibration smoke fixture failed."
}
