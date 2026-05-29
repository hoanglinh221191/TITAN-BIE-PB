param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_k4_channel_estimator_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_k4_channel_estimator_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_k4_channel_estimator_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside K4 channel estimator smoke root: $full"
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
  panels.reserve(156);
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
      unknowns[n + i] = 1.0e-4 * sign * bend;
    }
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &channel_calibration, const std::string &safety,
    double error_cap) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.70, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "sampled", 512, 4096,
          "translated", "naive", "sampled", safety, 0.0, 0.0, "hybrid",
          1.0e4, error_cap, 0.25, "consecutive-aspect", "off", "on", 1, 64,
          4, "channel-tail-eps", "summary", channel_calibration,
          "accepted-only");
  tree.taylor_order = 4;
  tree.derivative_mode = "analytic";
  return tree;
}

TitanPBE::Tabi::TreecodeStats run_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    const std::string &channel_calibration, const std::string &safety,
    double error_cap) {
  TitanPBE::Tabi::TreeWorkspace tree =
      make_tree(panels, executor, channel_calibration, safety, error_cap);
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
  const std::vector<double> k4_unknowns = make_unknowns(panels.size(), false);
  const std::vector<double> benign_unknowns =
      make_unknowns(panels.size(), true);
  TitanExecutor executor(2);

  const TitanPBE::Tabi::TreecodeStats baseline =
      run_tree(panels, k4_unknowns, executor, "stage3m", "audit", 1.0);
  const TitanPBE::Tabi::TreecodeStats calibrated =
      run_tree(panels, k4_unknowns, executor, "all-derivative-sqrt-eps",
               "audit", 1.0);

  require_true("baseline channel calibration reported",
               baseline.convergence_channel_calibration == "stage3m");
  require_true("K4 calibration reported",
               calibrated.convergence_channel_calibration ==
                   "all-derivative-sqrt-eps");
  require_true("accepted-only semantics reported",
               calibrated.channel_missed_semantics == "accepted-only");
  require_true("K4 estimator sampled nodes",
               calibrated.convergence_estimator_samples > 0);
  require_true("K4 amp increases",
               calibrated.convergence_channel_effective_jcp_amp_max[kK4] >
                   baseline.convergence_channel_effective_jcp_amp_max[kK4]);
  require_true("K4 estimated error increases",
               calibrated.convergence_channel_estimated_error_max[kK4] >
                   baseline.convergence_channel_estimated_error_max[kK4]);

  const double split_cap =
      0.5 * (baseline.convergence_channel_estimated_error_max[kK4] +
             calibrated.convergence_channel_estimated_error_max[kK4]);
  const TitanPBE::Tabi::TreecodeStats baseline_split =
      run_tree(panels, k4_unknowns, executor, "stage3m", "split", split_cap);
  const TitanPBE::Tabi::TreecodeStats calibrated_split =
      run_tree(panels, k4_unknowns, executor, "all-derivative-sqrt-eps",
               "split", split_cap);
  require_true("K4 calibration increases K4 split rejects",
               calibrated_split.convergence_channel_split_rejects[kK4] >
                   baseline_split.convergence_channel_split_rejects[kK4]);
  require_true("K4 split counter reports rejects",
               calibrated_split.convergence_channel_split_rejects[kK4] > 0);

  const TitanPBE::Tabi::TreecodeStats benign =
      run_tree(panels, benign_unknowns, executor, "all-derivative-sqrt-eps",
               "split", 1.0);
  require_true("benign path keeps far work", benign.far_node_interactions > 0);
  require_true("benign path not rejected wholesale",
               benign.convergence_estimator_split_rejects <
                   benign.convergence_estimator_samples);

  std::cout << "PBE TABI K4 channel estimator smoke PASS panels="
            << panels.size() << " base_k4_amp="
            << baseline.convergence_channel_effective_jcp_amp_max[kK4]
            << " calibrated_k4_amp="
            << calibrated.convergence_channel_effective_jcp_amp_max[kK4]
            << " base_k4_est="
            << baseline.convergence_channel_estimated_error_max[kK4]
            << " calibrated_k4_est="
            << calibrated.convergence_channel_estimated_error_max[kK4]
            << " split_cap=" << split_cap
            << " base_k4_split_rejects="
            << baseline_split.convergence_channel_split_rejects[kK4]
            << " k4_split_rejects="
            << calibrated_split.convergence_channel_split_rejects[kK4]
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
    throw "Failed to compile TABI K4 channel estimator smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI K4 channel estimator smoke fixture failed."
}
