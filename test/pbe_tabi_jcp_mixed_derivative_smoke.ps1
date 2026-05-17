param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_jcp_mixed_derivative_smoke'),
    [double]$PoissonTolerance = 5.0e-4,
    [double]$YukawaTolerance = 5.0e-2
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_jcp_mixed_derivative_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_jcp_mixed_derivative_smoke.exe'
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$PoissonToleranceText = $PoissonTolerance.ToString('G17', $InvariantCulture)
$YukawaToleranceText = $YukawaTolerance.ToString('G17', $InvariantCulture)

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside mixed-derivative smoke root: $full"
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

$source = @"
#include "titan_executor.h"
#include "titan_pbe/tabi_types.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPoissonTolerance = $PoissonToleranceText;
constexpr double kYukawaTolerance = $YukawaToleranceText;

void require_true(const std::string &label, bool value) {
  if (!value)
    throw std::runtime_error(label);
}

double rms(const std::vector<double> &values) {
  double sum = 0.0;
  for (double value : values)
    sum += value * value;
  return values.empty() ? 0.0
                        : std::sqrt(sum / static_cast<double>(values.size()));
}

double relative_rms(const std::vector<double> &reference,
                    const std::vector<double> &candidate) {
  require_true("same vector size", reference.size() == candidate.size());
  std::vector<double> diff(reference.size(), 0.0);
  for (std::size_t i = 0; i < reference.size(); ++i)
    diff[i] = candidate[i] - reference[i];
  return rms(diff) / std::max(1.0e-12, rms(reference));
}

TitanPBE::Tabi::Panel make_panel(double x, double y, double z,
                                 double nx, double ny, double nz,
                                 std::size_t index) {
  TitanPBE::Tabi::Panel p;
  p.x = x;
  p.y = y;
  p.z = z;
  const double norm = std::sqrt(std::max(1.0e-12, nx * nx + ny * ny + nz * nz));
  p.nx = nx / norm;
  p.ny = ny / norm;
  p.nz = nz / norm;
  p.area = 0.8 + 0.02 * static_cast<double>(index % 9);
  p.diagonal = 0.0;
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(48);
  std::size_t index = 0;
  for (int ix = 0; ix < 4; ++ix) {
    for (int iy = 0; iy < 4; ++iy) {
      panels.push_back(make_panel(
          -28.0 + 0.35 * ix, -6.0 + 0.31 * iy, 0.2 * ((ix + iy) % 3),
          -0.5 + 0.1 * ix, -0.2 + 0.15 * iy, 1.0, index++));
    }
  }
  for (int ix = 0; ix < 4; ++ix) {
    for (int iy = 0; iy < 4; ++iy) {
      panels.push_back(make_panel(
          18.0 + 0.29 * ix, 7.0 + 0.33 * iy, 2.0 + 0.1 * ((ix + 2 * iy) % 4),
          0.3 + 0.08 * ix, 0.5 - 0.06 * iy, 0.9, index++));
    }
  }
  for (int ix = 0; ix < 4; ++ix) {
    for (int iy = 0; iy < 4; ++iy) {
      panels.push_back(make_panel(
          30.0 + 0.27 * ix, -10.0 + 0.28 * iy, -4.0 + 0.09 * ((ix + iy) % 4),
          -0.4 + 0.05 * ix, 0.6 + 0.02 * iy, 0.7, index++));
    }
  }
  return panels;
}

std::vector<double> make_unknowns(std::size_t n, bool coupled) {
  std::vector<double> unknowns(coupled ? 2 * n : n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    unknowns[i] = 0.9 * std::sin(0.21 * static_cast<double>(i)) +
                  0.4 * std::cos(0.47 * static_cast<double>(i));
    if (coupled) {
      unknowns[n + i] = 0.025 * std::cos(0.17 * static_cast<double>(i)) -
                        0.011 * std::sin(0.31 * static_cast<double>(i));
    }
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &derivative_mode) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 4, 2.0, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "off", 512, 2048,
          "translated", "naive", "sampled");
  tree.taylor_order = 3;
  tree.derivative_mode = derivative_mode;
  return tree;
}

double run_case(const std::vector<TitanPBE::Tabi::Panel> &panels,
                TitanExecutor &executor, bool coupled, double kappa,
                TitanPBE::Tabi::TreecodeStats &analytic_stats,
                TitanPBE::Tabi::TreecodeStats &finite_stats) {
  const double eps_ratio = 80.0;
  const double lhs_phi = coupled ? 0.5 * (1.0 + eps_ratio) : 1.0;
  const double lhs_phin = 0.5 * (1.0 + 1.0 / eps_ratio);
  const std::vector<double> unknowns = make_unknowns(panels.size(), coupled);
  std::vector<double> analytic_out;
  std::vector<double> finite_out;
  TitanPBE::Tabi::TreeWorkspace analytic_tree =
      make_tree(panels, executor, "analytic");
  TitanPBE::Tabi::TreeWorkspace finite_tree =
      make_tree(panels, executor, "finite-diff-debug");
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, kappa, lhs_phi, lhs_phin, coupled,
      analytic_tree, "centroid", "off", executor, analytic_out,
      &analytic_stats);
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, eps_ratio, kappa, lhs_phi, lhs_phin, coupled,
      finite_tree, "centroid", "off", executor, finite_out, &finite_stats);
  return relative_rms(finite_out, analytic_out);
}

} // namespace

int main() {
  const std::vector<TitanPBE::Tabi::Panel> panels = build_fixture();
  TitanExecutor executor(2);

  TitanPBE::Tabi::TreecodeStats poisson_analytic;
  TitanPBE::Tabi::TreecodeStats poisson_finite;
  const double poisson_rel =
      run_case(panels, executor, false, 0.0, poisson_analytic, poisson_finite);
  require_true("poisson analytic treecode far interactions",
               poisson_analytic.far_node_interactions > 0);
  require_true("poisson derivative check", poisson_rel <= kPoissonTolerance);

  TitanPBE::Tabi::TreecodeStats yukawa_analytic;
  TitanPBE::Tabi::TreecodeStats yukawa_finite;
  const double yukawa_rel =
      run_case(panels, executor, true, 0.12, yukawa_analytic, yukawa_finite);
  require_true("yukawa analytic treecode far interactions",
               yukawa_analytic.far_node_interactions > 0);
  require_true("yukawa mixed derivative check", yukawa_rel <= kYukawaTolerance);

  std::cout << "PBE TABI JCP mixed-derivative smoke PASS panels="
            << panels.size() << " poisson_rel=" << poisson_rel
            << " yukawa_rel=" << yukawa_rel
            << " yukawa_far=" << yukawa_analytic.far_node_interactions
            << "\n";
  return 0;
}
"@

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
    throw "Failed to compile TABI JCP mixed-derivative smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI JCP mixed-derivative smoke fixture failed."
}
