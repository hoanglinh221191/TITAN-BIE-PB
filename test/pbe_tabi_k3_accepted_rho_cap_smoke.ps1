param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_k3_accepted_rho_cap_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_k3_accepted_rho_cap_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_k3_accepted_rho_cap_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside K3 accepted-rho cap smoke root: $full"
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
  const double nx = 0.45 + 0.03 * static_cast<double>(index % 7);
  const double ny = ((index % 2) == 0 ? 0.55 : -0.50);
  const double nz = 0.70;
  const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  p.nx = nx / norm;
  p.ny = ny / norm;
  p.nz = nz / norm;
  p.area = 0.12 + 0.015 * static_cast<double>(index % 11);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_two_cluster_fixture(
    double separation) {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(96);
  std::size_t index = 0;
  for (int side = 0; side < 2; ++side) {
    const double cx = side == 0 ? 0.0 : separation;
    for (int ix = -3; ix <= 3; ++ix) {
      for (int iy = -3; iy <= 3; ++iy) {
        const double x = cx + 0.70 * static_cast<double>(ix);
        const double y = 0.70 * static_cast<double>(iy);
        const double z =
            0.08 * static_cast<double>((ix * ix + iy + 17) % 5);
        panels.push_back(make_panel(x, y, z, index++));
      }
    }
  }
  return panels;
}

std::vector<double> make_k3_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    unknowns[i] = 1.0e-5 * sign;
    unknowns[n + i] =
        12.0 * sign * (1.0 + 0.25 * std::sin(0.31 * static_cast<double>(i)));
  }
  return unknowns;
}

TitanPBE::Tabi::TreecodeStats run_case(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    double theta, double k3_cap) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, theta, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "sampled", 256, 4096,
          "translated", "naive", "sampled", "audit", 0.0, 0.0, "hybrid",
          1.0e12, 1.0e12, 0.25, "consecutive-aspect", "off", "on", 1, 64, 4,
          "channel-tail-eps", "summary", "polynomial-eps", "accepted-only",
          "global", "summary", k3_cap, "center", 0.25);
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
  TitanExecutor executor(2);
  const std::vector<TitanPBE::Tabi::Panel> high_rho =
      build_two_cluster_fixture(7.4);
  const std::vector<double> high_unknowns = make_k3_unknowns(high_rho.size());
  const TitanPBE::Tabi::TreecodeStats baseline =
      run_case(high_rho, high_unknowns, executor, 0.90, 0.0);
  const TitanPBE::Tabi::TreecodeStats capped =
      run_case(high_rho, high_unknowns, executor, 0.90, 0.58);
  require_true("baseline has accepted far work", baseline.far_node_interactions > 0);
  require_true("baseline has no K3 cap rejects",
               baseline.convergence_k3_rho_cap_rejects == 0);
  require_true("cap rejects high-rho accepted far candidates",
               capped.convergence_k3_rho_cap_rejects > 0);
  require_true("cap increases leaf/direct work",
               capped.leaf_panel_interactions > baseline.leaf_panel_interactions);

  const std::vector<TitanPBE::Tabi::Panel> benign =
      build_two_cluster_fixture(22.0);
  const std::vector<double> benign_unknowns = make_k3_unknowns(benign.size());
  const TitanPBE::Tabi::TreecodeStats benign_capped =
      run_case(benign, benign_unknowns, executor, 0.55, 0.58);
  require_true("benign low-rho case keeps far work",
               benign_capped.far_node_interactions > 0);
  require_true("benign low-rho case is not rejected by K3 cap",
               benign_capped.convergence_k3_rho_cap_rejects == 0);

  std::cout << "PBE TABI K3 accepted-rho cap smoke PASS baseline_far="
            << baseline.far_node_interactions
            << " capped_rejects=" << capped.convergence_k3_rho_cap_rejects
            << " capped_panels="
            << capped.convergence_k3_rho_cap_reject_panels
            << " benign_far=" << benign_capped.far_node_interactions << "\n";
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
    throw "Failed to compile TABI K3 accepted-rho cap smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI K3 accepted-rho cap smoke fixture failed."
}
