param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_channel_missed_semantics_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_channel_missed_semantics_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_channel_missed_semantics_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside channel semantics smoke root: $full"
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
  const double norm = std::sqrt(std::max(1.0e-12, nx * nx + ny * ny + nz * nz));
  p.nx = nx / norm;
  p.ny = ny / norm;
  p.nz = nz / norm;
  p.area = 0.20 + 0.012 * static_cast<double>(index % 7);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(132);
  std::size_t index = 0;
  for (int i = 0; i < 28; ++i) {
    const double a = 0.2243994752564138 * static_cast<double>(i);
    panels.push_back(make_panel(-30.0 + 1.0 * std::cos(a),
                                0.9 * std::sin(a),
                                0.08 * std::sin(2.0 * a), 0.1, 0.0, 1.0,
                                index++));
  }
  for (int ix = -5; ix <= 5; ++ix) {
    if (ix == 0)
      continue;
    for (int iy = -4; iy <= 4; ++iy) {
      panels.push_back(make_panel(22.0 + 0.88 * static_cast<double>(ix),
                                  0.70 * static_cast<double>(iy),
                                  0.05 * static_cast<double>((ix + iy) % 5),
                                  0.25, 0.15, 1.0, index++));
    }
  }
  return panels;
}

std::vector<double> make_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sign = (i % 2 == 0) ? 1.0 : -1.0;
    const double wave = 1.0 + 0.25 * std::sin(0.23 * static_cast<double>(i));
    unknowns[i] = 30.0 * sign * wave;
    unknowns[n + i] = 1.0e-4 * sign * wave;
  }
  return unknowns;
}

TitanPBE::Tabi::TreecodeStats run_case(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    const std::string &safety, double error_cap,
    const std::string &missed_semantics) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 8, 0.70, "exact", "scalar", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "sampled", 512, 4096,
          "translated", "naive", "sampled", safety, 0.0, 0.0, "hybrid",
          1.0e4, error_cap, 0.25, "consecutive-aspect", "off", "on", 1, 64,
          4, "channel-tail-eps", "summary", "stage3m", missed_semantics);
  tree.taylor_order = 3;
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
  const std::vector<double> unknowns = make_unknowns(panels.size());
  TitanExecutor executor(2);

  const TitanPBE::Tabi::TreecodeStats accepted =
      run_case(panels, unknowns, executor, "audit", 1.0e-12,
               "accepted-only");
  require_true("accepted-only mode reported",
               accepted.channel_missed_semantics == "accepted-only");
  require_true("accepted high-error sampled pairs",
               accepted.far_expansion_audit_pairs > 0);
  require_true("accepted high-error increments missed",
               accepted.convergence_channel_high_error_missed[kK2] > 0 ||
                   accepted.convergence_estimator_high_error_missed > 0);

  const TitanPBE::Tabi::TreecodeStats split =
      run_case(panels, unknowns, executor, "split", 1.0e-12,
               "accepted-only");
  const std::size_t accepted_missed =
      accepted.convergence_channel_high_error_missed[kK2] +
      accepted.convergence_estimator_high_error_missed;
  const std::size_t split_missed =
      split.convergence_channel_high_error_missed[kK2] +
      split.convergence_estimator_high_error_missed;
  const std::size_t split_caught_by_split =
      split.convergence_channel_high_error_caught_by_split[kK2] +
      split.convergence_estimator_high_error_caught_by_split;
  require_true("split rejects candidates",
               split.convergence_estimator_split_rejects > 0);
  require_true("split high-error caught by split",
               split_caught_by_split > 0);
  require_true("split-rejected high-error does not inflate missed",
               split_missed < accepted_missed);

  std::cout << "PBE TABI channel missed semantics smoke PASS panels="
            << panels.size() << " accepted_missed="
            << accepted_missed << " split_missed=" << split_missed
            << " split_caught_by_split="
            << split_caught_by_split
            << " split_rejects="
            << split.convergence_estimator_split_rejects << "\n";
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
    throw "Failed to compile TABI channel missed semantics smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI channel missed semantics smoke fixture failed."
}
