param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_selective_k3_split_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_selective_k3_split_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_selective_k3_split_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside selective K3 split smoke root: $full"
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
  p.area = 0.14 + 0.020 * static_cast<double>(index % 11);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(160);
  std::size_t index = 0;
  for (int i = 0; i < 40; ++i) {
    const double a = 0.15707963267948966 * static_cast<double>(i);
    panels.push_back(make_panel(-34.0 + 1.6 * std::cos(a),
                                1.2 * std::sin(a),
                                0.12 * std::sin(3.0 * a),
                                0.35 * std::cos(a), 0.65 * std::sin(a),
                                0.70, index++));
  }
  for (int ix = -6; ix <= 6; ++ix) {
    if (ix == 0)
      continue;
    for (int iy = -4; iy <= 4; ++iy) {
      panels.push_back(make_panel(24.0 + 0.92 * static_cast<double>(ix),
                                  0.70 * static_cast<double>(iy),
                                  0.06 * static_cast<double>((ix - iy) % 7),
                                  0.55 * static_cast<double>(ix),
                                  0.45 * static_cast<double>(iy), 1.0,
                                  index++));
    }
  }
  for (int i = 0; i < 28; ++i) {
    const double a = 0.2243994752564138 * static_cast<double>(i);
    panels.push_back(make_panel(50.0 + 1.9 * std::cos(a),
                                -4.5 + 1.4 * std::sin(a),
                                0.09 * std::cos(2.0 * a), -0.25, 0.60,
                                0.75 + 0.15 * std::sin(a), index++));
  }
  return panels;
}

std::vector<double> make_k3_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    const double wave = 1.0 + 0.30 * std::sin(0.29 * static_cast<double>(i));
    unknowns[i] = 1.0e-5 * sign * wave;
    unknowns[n + i] = 38.0 * sign * wave;
  }
  return unknowns;
}

std::vector<double> make_k4_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    const double wave = 1.0 + 0.35 * std::cos(0.21 * static_cast<double>(i));
    unknowns[i] = 42.0 * sign * wave;
    unknowns[n + i] = 0.02;
  }
  return unknowns;
}

std::vector<double> make_benign_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    unknowns[i] = 0.03 * (1.0 + 0.10 * std::sin(0.13 * i));
    unknowns[n + i] = 0.02 * (1.0 + 0.10 * std::cos(0.17 * i));
  }
  return unknowns;
}

TitanPBE::Tabi::TreecodeStats run_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    const std::string &split_policy, double error_cap,
    const std::string &channel_calibration = "polynomial-eps") {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.70, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "sampled", 512, 4096,
          "translated", "naive", "sampled", "split", 0.0, 0.0, "hybrid",
          1.0e12, error_cap, 0.25, "consecutive-aspect", "off", "on", 1, 64,
          4, "channel-tail-eps", "summary", channel_calibration,
          "accepted-only", split_policy);
  tree.taylor_order = 4;
  tree.derivative_mode = "analytic";
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
  const std::vector<double> k3_unknowns = make_k3_unknowns(panels.size());
  const std::vector<double> k4_unknowns = make_k4_unknowns(panels.size());
  const std::vector<double> benign_unknowns =
      make_benign_unknowns(panels.size());

  const TitanPBE::Tabi::TreecodeStats k3_hot =
      run_tree(panels, k3_unknowns, executor, "k3-hot", 1.0e-6);
  require_true("k3-hot policy reported",
               k3_hot.convergence_split_policy == "k3-hot");
  require_true("K3 dominated node split under k3-hot",
               k3_hot.convergence_channel_split_rejects[kK3] > 0);

  const TitanPBE::Tabi::TreecodeStats k4_global =
      run_tree(panels, k4_unknowns, executor, "global", 10.0, "k4-eps");
  const TitanPBE::Tabi::TreecodeStats k4_hot =
      run_tree(panels, k4_unknowns, executor, "k3-hot", 10.0, "k4-eps");
  require_true("global policy reports split rejects",
               k4_global.convergence_estimator_split_rejects > 0);
  require_true("k3-hot reduces K4 dominated split rejects",
               k4_hot.convergence_estimator_split_rejects <
                   k4_global.convergence_estimator_split_rejects);

  const TitanPBE::Tabi::TreecodeStats benign =
      run_tree(panels, benign_unknowns, executor, "k3-hot", 1.0);
  require_true("benign accepted-only missed stays zero",
               benign.convergence_estimator_high_error_missed == 0);
  require_true("benign path keeps accepted far work",
               benign.far_node_interactions > 0);

  std::cout << "PBE TABI selective K3 split smoke PASS panels="
            << panels.size() << " k3_hot_k3_splits="
            << k3_hot.convergence_channel_split_rejects[kK3]
            << " k4_global_splits="
            << k4_global.convergence_estimator_split_rejects
            << " k4_hot_splits="
            << k4_hot.convergence_estimator_split_rejects
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
    throw "Failed to compile TABI selective K3 split smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI selective K3 split smoke fixture failed."
}
