param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_channel_estimator_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_channel_estimator_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_channel_estimator_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside channel estimator smoke root: $full"
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

constexpr std::size_t kK2 = 1;

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
  p.area = 0.20 + 0.01 * static_cast<double>(index % 11);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(150);
  std::size_t index = 0;
  for (int i = 0; i < 34; ++i) {
    const double a = 0.18479956785822313 * static_cast<double>(i);
    panels.push_back(make_panel(-32.0 + 1.1 * std::cos(a),
                                1.1 * std::sin(a),
                                0.08 * std::sin(2.0 * a), 0.0, 0.2, 1.0,
                                index++));
  }
  for (int ix = -6; ix <= 6; ++ix) {
    if (ix == 0)
      continue;
    for (int iy = -4; iy <= 4; ++iy) {
      const double x = 23.0 + 0.90 * static_cast<double>(ix);
      const double y = 0.70 * static_cast<double>(iy);
      const double z = 0.05 * static_cast<double>((2 * ix + iy) % 5);
      panels.push_back(make_panel(x, y, z, 0.2, 0.1, 1.0, index++));
    }
  }
  for (int i = 0; i < 22; ++i) {
    const double a = 0.28559933214452665 * static_cast<double>(i);
    panels.push_back(make_panel(46.0 + 1.8 * std::cos(a),
                                -5.2 + 1.8 * std::sin(a),
                                0.06 * std::cos(3.0 * a), -0.2, 0.1, 1.0,
                                index++));
  }
  return panels;
}

std::vector<double> make_unknowns(std::size_t n, bool benign) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    const double wave = 1.0 + 0.25 * std::sin(0.23 * static_cast<double>(i));
    if (benign) {
      unknowns[i] = 0.05 * wave;
      unknowns[n + i] = 0.04 * wave;
    } else {
      unknowns[i] = 30.0 * sign * wave;
      unknowns[n + i] = 1.0e-4 * sign * wave;
    }
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &model, const std::string &safety, double error_cap,
    const std::string &audit_channel) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.70, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "sampled", 512, 4096,
          "translated", "naive", "sampled", safety, 0.0, 0.0, "hybrid",
          1.0e4, error_cap, 0.25, "consecutive-aspect", "off", "on", 1, 64,
          4, model, audit_channel);
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  return tree;
}

TitanPBE::Tabi::TreecodeStats run_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    const std::string &model, const std::string &safety, double error_cap,
    const std::string &audit_channel) {
  TitanPBE::Tabi::TreeWorkspace tree =
      make_tree(panels, executor, model, safety, error_cap, audit_channel);
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
  const std::vector<double> k2_unknowns = make_unknowns(panels.size(), false);
  const std::vector<double> benign_unknowns =
      make_unknowns(panels.size(), true);
  TitanExecutor executor(2);

  const TitanPBE::Tabi::TreecodeStats baseline =
      run_tree(panels, k2_unknowns, executor, "stage3l", "audit", 1.0,
               "summary");
  const TitanPBE::Tabi::TreecodeStats channel =
      run_tree(panels, k2_unknowns, executor, "channel-tail-eps", "audit",
               1.0, "summary");

  require_true("baseline model reported",
               baseline.convergence_estimator_model == "stage3l");
  require_true("channel model reported",
               channel.convergence_estimator_model == "channel-tail-eps");
  require_true("channel estimator sampled nodes",
               channel.convergence_estimator_samples > 0);
  require_true("channel K2 estimate reported",
               channel.convergence_channel_estimated_error_max[kK2] > 0.0);
  require_true("channel audit sampled pairs",
               channel.far_expansion_audit_pairs > 0);

  const double split_cap = std::max(
      1.0e-12, 0.5 * channel.convergence_channel_estimated_error_max[kK2]);
  const TitanPBE::Tabi::TreecodeStats split =
      run_tree(panels, k2_unknowns, executor, "channel-tail-eps", "split",
               split_cap, "summary");
  require_true("channel estimator split rejects",
               split.convergence_estimator_split_rejects > 0);
  require_true("channel estimator catches or reports high-error oracle",
               split.convergence_channel_high_error_caught[kK2] +
                       split.convergence_channel_high_error_missed[kK2] >
                   0 ||
                   split.far_expansion_audit_pairs > 0);

  const TitanPBE::Tabi::TreecodeStats benign =
      run_tree(panels, benign_unknowns, executor, "channel-tail-eps", "split",
               1.0, "summary");
  require_true("benign path keeps far work", benign.far_node_interactions > 0);
  require_true("benign path not rejected wholesale",
               benign.convergence_estimator_split_rejects <
                   benign.convergence_estimator_samples);

  std::cout << "PBE TABI channel estimator smoke PASS panels="
            << panels.size() << " baseline_err="
            << baseline.convergence_estimator_error_max << " k2_est="
            << channel.convergence_channel_estimated_error_max[kK2]
            << " split_cap=" << split_cap << " split_rejects="
            << split.convergence_estimator_split_rejects << " benign_far="
            << benign.far_node_interactions << "\n";
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
    throw "Failed to compile TABI channel estimator smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI channel estimator smoke fixture failed."
}
