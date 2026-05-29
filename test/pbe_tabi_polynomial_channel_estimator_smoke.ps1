param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_polynomial_channel_estimator_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_polynomial_channel_estimator_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_polynomial_channel_estimator_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside polynomial channel estimator smoke root: $full"
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
  p.area = 0.16 + 0.018 * static_cast<double>(index % 10);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(168);
  std::size_t index = 0;
  for (int i = 0; i < 36; ++i) {
    const double a = 0.17453292519943295 * static_cast<double>(i);
    panels.push_back(make_panel(-36.0 + 1.35 * std::cos(a),
                                1.15 * std::sin(a),
                                0.14 * std::sin(3.0 * a),
                                0.45 * std::cos(a), 0.35 * std::sin(a),
                                0.80, index++));
  }
  for (int ix = -6; ix <= 6; ++ix) {
    if (ix == 0)
      continue;
    for (int iy = -4; iy <= 4; ++iy) {
      const double x = 24.0 + 0.96 * static_cast<double>(ix);
      const double y = 0.78 * static_cast<double>(iy);
      const double z = 0.07 * static_cast<double>((3 * ix + iy) % 7);
      panels.push_back(make_panel(x, y, z, 0.55 * static_cast<double>(ix),
                                  0.45 * static_cast<double>(iy), 1.0,
                                  index++));
    }
  }
  for (int i = 0; i < 24; ++i) {
    const double a = 0.2617993877991494 * static_cast<double>(i);
    panels.push_back(make_panel(49.0 + 2.1 * std::cos(a),
                                -4.8 + 1.5 * std::sin(a),
                                0.10 * std::cos(2.0 * a), -0.35, 0.55,
                                0.75 + 0.20 * std::sin(a), index++));
  }
  return panels;
}

std::vector<double> make_unknowns(std::size_t n, bool benign) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    const double wave = 1.0 + 0.40 * std::sin(0.19 * static_cast<double>(i));
    const double bend = 1.0 + 0.15 * std::cos(0.41 * static_cast<double>(i));
    if (benign) {
      unknowns[i] = 0.04 * wave;
      unknowns[n + i] = 0.03 * bend;
    } else {
      unknowns[i] = 34.0 * sign * wave;
      unknowns[n + i] = 27.0 * sign * bend;
    }
  }
  return unknowns;
}

std::vector<double> make_k3_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    const double wave = 1.0 + 0.35 * std::sin(0.23 * static_cast<double>(i));
    unknowns[i] = 1.0e-4 * sign * wave;
    unknowns[n + i] = 28.0 * sign * wave;
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &channel_calibration, const std::string &safety,
    double error_cap, int order) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.70, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "sampled", 512, 4096,
          "translated", "naive", "sampled", safety, 0.0, 0.0, "hybrid",
          1.0e12, error_cap, 0.25, "consecutive-aspect", "off", "on", 1, 64,
          4, "channel-tail-eps", "summary", channel_calibration,
          "accepted-only");
  tree.taylor_order = order;
  tree.derivative_mode = "analytic";
  return tree;
}

TitanPBE::Tabi::TreecodeStats run_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    const std::string &channel_calibration, const std::string &safety,
    double error_cap, int order) {
  TitanPBE::Tabi::TreeWorkspace tree =
      make_tree(panels, executor, channel_calibration, safety, error_cap,
                order);
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
  const std::vector<double> derivative_unknowns =
      make_unknowns(panels.size(), false);
  const std::vector<double> k3_unknowns = make_k3_unknowns(panels.size());
  const std::vector<double> benign_unknowns =
      make_unknowns(panels.size(), true);
  TitanExecutor executor(2);

  const TitanPBE::Tabi::TreecodeStats all_deriv_o3 = run_tree(
      panels, derivative_unknowns, executor, "all-derivative-sqrt-eps",
      "audit", 1.0, 3);
  const TitanPBE::Tabi::TreecodeStats poly_o3 =
      run_tree(panels, derivative_unknowns, executor, "polynomial-eps",
               "audit", 1.0, 3);
  const TitanPBE::Tabi::TreecodeStats poly_o4 =
      run_tree(panels, derivative_unknowns, executor, "polynomial-eps",
               "audit", 1.0, 4);

  require_true("polynomial calibration reported",
               poly_o3.convergence_channel_calibration == "polynomial-eps");
  require_true("polynomial estimator sampled nodes",
               poly_o3.convergence_estimator_samples > 0);
  require_true("K3 amp increases versus all-derivative",
               poly_o3.convergence_channel_effective_jcp_amp_max[kK3] >
                   all_deriv_o3
                       .convergence_channel_effective_jcp_amp_max[kK3]);
  require_true("K3 estimated error increases versus all-derivative",
               poly_o3.convergence_channel_estimated_error_max[kK3] >
                   all_deriv_o3.convergence_channel_estimated_error_max[kK3]);
  require_true("K4 amp uses polynomial second-derivative scale",
               poly_o3.convergence_channel_effective_jcp_amp_max[kK4] >
                   all_deriv_o3
                       .convergence_channel_effective_jcp_amp_max[kK4]);
  require_true("K4 estimated error increases versus all-derivative",
               poly_o3.convergence_channel_estimated_error_max[kK4] >
                   all_deriv_o3.convergence_channel_estimated_error_max[kK4]);
  require_true("K2 polynomial channel remains active",
               poly_o3.convergence_channel_estimated_error_max[kK2] > 0.0);
  require_true("order 4 increases K2 polynomial amp",
               poly_o4.convergence_channel_effective_jcp_amp_max[kK2] >
                   poly_o3.convergence_channel_effective_jcp_amp_max[kK2]);
  require_true("order 4 increases K3 polynomial amp",
               poly_o4.convergence_channel_effective_jcp_amp_max[kK3] >
                   poly_o3.convergence_channel_effective_jcp_amp_max[kK3]);
  require_true("order 4 increases K4 polynomial amp",
               poly_o4.convergence_channel_effective_jcp_amp_max[kK4] >
                   poly_o3.convergence_channel_effective_jcp_amp_max[kK4]);

  const TitanPBE::Tabi::TreecodeStats all_deriv_k3 = run_tree(
      panels, k3_unknowns, executor, "all-derivative-sqrt-eps", "audit", 1.0,
      3);
  const TitanPBE::Tabi::TreecodeStats poly_k3 =
      run_tree(panels, k3_unknowns, executor, "polynomial-eps", "audit", 1.0,
               3);
  const double k3_split_cap =
      0.5 * (all_deriv_k3.convergence_channel_estimated_error_max[kK3] +
             poly_k3.convergence_channel_estimated_error_max[kK3]);
  const TitanPBE::Tabi::TreecodeStats all_deriv_split = run_tree(
      panels, k3_unknowns, executor, "all-derivative-sqrt-eps", "split",
      k3_split_cap, 3);
  const TitanPBE::Tabi::TreecodeStats poly_split =
      run_tree(panels, k3_unknowns, executor, "polynomial-eps", "split",
               k3_split_cap, 3);
  require_true("polynomial K3 split counter reports rejects",
               poly_split.convergence_channel_split_rejects[kK3] > 0);
  require_true("polynomial split has K4 rejects",
               poly_split.convergence_channel_split_rejects[kK4] > 0);

  const TitanPBE::Tabi::TreecodeStats benign =
      run_tree(panels, benign_unknowns, executor, "polynomial-eps", "split",
               1.0, 3);
  require_true("benign path keeps far work", benign.far_node_interactions > 0);
  require_true("benign path not rejected wholesale",
               benign.convergence_estimator_split_rejects <
                   benign.convergence_estimator_samples);

  std::cout << "PBE TABI polynomial channel estimator smoke PASS panels="
            << panels.size() << " all_o3_k3_amp="
            << all_deriv_o3.convergence_channel_effective_jcp_amp_max[kK3]
            << " poly_o3_k3_amp="
            << poly_o3.convergence_channel_effective_jcp_amp_max[kK3]
            << " poly_o3_k2_amp="
            << poly_o3.convergence_channel_effective_jcp_amp_max[kK2]
            << " poly_o4_k2_amp="
            << poly_o4.convergence_channel_effective_jcp_amp_max[kK2]
            << " poly_o4_k3_amp="
            << poly_o4.convergence_channel_effective_jcp_amp_max[kK3]
            << " poly_o3_k4_amp="
            << poly_o3.convergence_channel_effective_jcp_amp_max[kK4]
            << " poly_o4_k4_amp="
            << poly_o4.convergence_channel_effective_jcp_amp_max[kK4]
            << " k3_split_cap=" << k3_split_cap
            << " all_k3_split="
            << all_deriv_split.convergence_channel_split_rejects[kK3]
            << " poly_k3_split="
            << poly_split.convergence_channel_split_rejects[kK3]
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
    throw "Failed to compile TABI polynomial channel estimator smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI polynomial channel estimator smoke fixture failed."
}
