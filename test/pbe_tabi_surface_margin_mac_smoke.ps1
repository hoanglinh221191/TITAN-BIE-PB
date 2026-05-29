param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_surface_margin_mac_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_surface_margin_mac_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_surface_margin_mac_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside surface-margin MAC smoke root: $full"
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
  const double nx = 0.20 + 0.04 * static_cast<double>(index % 9);
  const double ny = ((index % 3) == 0 ? -0.50 : 0.55);
  const double nz = 0.75;
  const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  p.nx = nx / norm;
  p.ny = ny / norm;
  p.nz = nz / norm;
  p.area = 0.10 + 0.012 * static_cast<double>(index % 13);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_hollow_fixture(double separation) {
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
        panels.push_back(make_panel(x, y, z, index++));
      }
    }
  }
  return panels;
}

std::vector<double> make_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    unknowns[i] = 1.8 * sign * (1.0 + 0.15 * std::sin(0.23 * i));
    unknowns[n + i] =
        0.9 * sign * (1.0 + 0.20 * std::cos(0.37 * i));
  }
  return unknowns;
}

TitanPBE::Tabi::TreecodeStats run_case(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    const std::string &distance_policy, double margin_factor) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.82, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "sampled", 256, 4096,
          "translated", "naive", "sampled", "audit", 0.0, 0.0, "hybrid",
          1.0e12, 1.0e12, 0.25, "consecutive-aspect", "off", "on", 1, 64, 4,
          "channel-tail-eps", "summary", "polynomial-eps", "accepted-only",
          "global", "summary", 0.0, distance_policy, margin_factor);
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
  const std::vector<TitanPBE::Tabi::Panel> hollow = build_hollow_fixture(8.0);
  const std::vector<double> hollow_unknowns = make_unknowns(hollow.size());
  const TitanPBE::Tabi::TreecodeStats center =
      run_case(hollow, hollow_unknowns, executor, "center", 0.25);
  const TitanPBE::Tabi::TreecodeStats surface =
      run_case(hollow, hollow_unknowns, executor, "surface-margin", 0.25);

  require_true("center policy has accepted far work",
               center.far_node_interactions > 0);
  require_true("surface-margin policy reported",
               surface.mac_distance_policy == "surface-margin");
  require_true("surface-margin adds MAC rejects",
               surface.mac_surface_margin_rejects > 0);
  require_true("surface-margin reduces accepted far work",
               surface.far_node_interactions < center.far_node_interactions);
  require_true("surface-margin increases leaf/direct work",
               surface.leaf_panel_interactions > center.leaf_panel_interactions);

  const std::vector<TitanPBE::Tabi::Panel> benign = build_hollow_fixture(25.0);
  const std::vector<double> benign_unknowns = make_unknowns(benign.size());
  const TitanPBE::Tabi::TreecodeStats benign_surface =
      run_case(benign, benign_unknowns, executor, "surface-margin", 0.25);
  require_true("benign compact/far case keeps far work",
               benign_surface.far_node_interactions > 0);
  require_true("benign compact/far case is less affected by margin",
               benign_surface.mac_surface_margin_rejects <
                   surface.mac_surface_margin_rejects);

  std::cout << "PBE TABI surface-margin MAC smoke PASS center_far="
            << center.far_node_interactions
            << " surface_far=" << surface.far_node_interactions
            << " surface_rejects=" << surface.mac_surface_margin_rejects
            << " benign_far=" << benign_surface.far_node_interactions << "\n";
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
    throw "Failed to compile TABI surface-margin MAC smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI surface-margin MAC smoke fixture failed."
}
