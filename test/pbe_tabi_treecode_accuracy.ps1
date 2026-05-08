param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_treecode_accuracy')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_treecode_accuracy.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_treecode_accuracy.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI treecode accuracy root: $full"
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

double rms(const std::vector<double> &values) {
  double sum = 0.0;
  for (double value : values)
    sum += value * value;
  return values.empty() ? 0.0 : std::sqrt(sum / static_cast<double>(values.size()));
}

double relative_rms(const std::vector<double> &reference,
                    const std::vector<double> &candidate) {
  std::vector<double> diff(reference.size(), 0.0);
  for (std::size_t i = 0; i < reference.size(); ++i)
    diff[i] = candidate[i] - reference[i];
  return rms(diff) / std::max(1.0e-12, rms(reference));
}

std::vector<TitanPBE::Tabi::Panel> build_panels() {
  TitanPBE::Tabi::Options opt;
  opt.surface_mode = "vdw";
  opt.surface_geometry = "triangle";
  opt.diagonal_mode = "solid-angle";
  TitanPBE::Tabi::Result result;
  std::vector<TitanPBE::Tabi::Panel> panels =
      TitanPBE::Tabi::build_surface(
          {{-3.0, -0.2, 0.0, 0.4, 1.55},
           {-1.5, 0.3, 0.1, -0.2, 1.70},
           {0.0, -0.1, -0.3, 0.3, 1.65},
           {1.5, 0.2, 0.2, -0.4, 1.70},
           {3.0, -0.2, 0.0, 0.2, 1.55}},
          opt, &result);
  require_true("triangle panels", panels.size() > 900);
  require_true("solid-angle diagonal", result.diag_min > 0.0);
  return panels;
}

} // namespace

int main() {
  std::vector<TitanPBE::Tabi::Panel> panels = build_panels();
  std::vector<double> density(panels.size(), 0.0);
  for (std::size_t i = 0; i < density.size(); ++i)
    density[i] = std::sin(0.017 * static_cast<double>(i)) +
                 0.25 * std::cos(0.031 * static_cast<double>(i));

  TitanExecutor executor(2);
  std::vector<double> direct;
  TitanPBE::Tabi::TreecodeStats direct_stats;
  TitanPBE::Tabi::apply_direct_operator(panels, density, 0.42, 0.12,
                                        executor, direct, &direct_stats);
  require_true("direct matvec", direct_stats.matvec_mode == "direct");

  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(panels, executor);
  require_true("tree nodes", tree.nodes.size() > 1);

  std::vector<double> tree3;
  tree.taylor_order = 3;
  TitanPBE::Tabi::TreecodeStats stats3;
  TitanPBE::Tabi::apply_treecode_operator(panels, density, 0.42, 0.12, tree,
                                          executor, tree3, &stats3);
  const double rel3 = relative_rms(direct, tree3);
  require_true("treecode p3 mode", stats3.matvec_mode == "treecode");
  require_true("treecode p3 far interactions", stats3.far_interactions > 0);
  require_true("treecode p3 relative error", rel3 < 2.5e-3);

  std::vector<double> tree4;
  tree.taylor_order = 4;
  TitanPBE::Tabi::TreecodeStats stats4;
  TitanPBE::Tabi::apply_treecode_operator(panels, density, 0.42, 0.12, tree,
                                          executor, tree4, &stats4);
  const double rel4 = relative_rms(direct, tree4);
  require_true("treecode p4 order", stats4.taylor_order == 4);
  require_true("treecode p4 not worse", rel4 <= rel3 * 1.50 + 5.0e-5);

  std::cout << "PBE TABI treecode accuracy PASS panels=" << panels.size()
            << " rel3=" << rel3 << " rel4=" << rel4
            << " far=" << stats3.far_interactions << "\n";
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
    (Join-Path $RepoRoot 'titan_pbe\tabi_surface.cpp'),
    (Join-Path $RepoRoot 'titan_pbe\tabi_treecode.cpp'),
    (Join-Path $RepoRoot 'titan_surface.cpp'),
    (Join-Path $RepoRoot 'titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI treecode accuracy fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI treecode accuracy fixture failed."
}
