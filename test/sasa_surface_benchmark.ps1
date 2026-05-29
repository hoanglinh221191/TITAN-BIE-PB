param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\sasa_surface_benchmark')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'sasa_surface_benchmark.cpp'
$ExePath = Join-Path $WorkRoot 'sasa_surface_benchmark.exe'
$CsvPath = Join-Path $WorkRoot 'sasa_surface_benchmark.csv'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing path outside SASA surface benchmark root: $full"
    }
}

Assert-UnderPath -Path $WorkRoot -Root $WorkRoot
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
    throw "Missing clang++ compiler: $Compiler"
}

$source = @'
#include "titan_executor.h"
#include "titan_surface.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

std::vector<TitanSurface::Atom> deterministic_atoms(std::size_t count) {
  std::vector<TitanSurface::Atom> atoms;
  atoms.reserve(count);
  std::uint64_t state = 0x6a09e667f3bcc909ull + count;
  auto next_unit = [&]() {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return static_cast<double>((state >> 11) & ((1ull << 53) - 1ull)) /
           static_cast<double>(1ull << 53);
  };
  const double span = std::max(18.0, std::cbrt(static_cast<double>(count)) * 6.0);
  for (std::size_t i = 0; i < count; ++i) {
    atoms.push_back({span * next_unit() - 0.5 * span,
                     span * next_unit() - 0.5 * span,
                     span * next_unit() - 0.5 * span,
                     1.2 + 0.7 * next_unit()});
  }
  return atoms;
}

double wall_ms(std::chrono::steady_clock::time_point begin,
               std::chrono::steady_clock::time_point end) {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

TitanSurface::Workspace run_case(const std::vector<TitanSurface::Atom> &atoms,
                                 const std::string &mode,
                                 TitanExecutor &executor, double &wall) {
  TitanSurface::Options options;
  options.acceleration = mode;
  options.emit_panels = false;
  options.emit_neighbor_workspace = false;
  const auto begin = std::chrono::steady_clock::now();
  TitanSurface::Workspace workspace =
      TitanSurface::build_workspace(atoms, options, {}, &executor);
  wall = wall_ms(begin, std::chrono::steady_clock::now());
  return workspace;
}

} // namespace

int main() {
  TitanExecutor executor(4);
  std::cout
      << "case,mode,atoms,targets,point_count,cell_count,neighbor_mean,"
         "neighbor_max,padded_slots,build_ms,neighbor_build_ms,exposure_ms,"
         "total_ms,speedup_vs_legacy\n";
  for (std::size_t count : {64u, 256u, 1024u}) {
    const std::vector<TitanSurface::Atom> atoms = deterministic_atoms(count);
    double legacy_wall = 0.0;
    TitanSurface::Workspace legacy =
        run_case(atoms, "legacy", executor, legacy_wall);
    for (const std::string mode : {std::string("legacy"),
                                   std::string("reordered_padded")}) {
      double wall = 0.0;
      TitanSurface::Workspace workspace =
          mode == "legacy" ? legacy : run_case(atoms, mode, executor, wall);
      if (mode == "legacy")
        wall = legacy_wall;
      const double speedup = wall > 0.0 ? legacy_wall / wall : 0.0;
      const auto &d = workspace.diagnostics;
      std::cout << "atoms_" << count << ',' << mode << ',' << count << ','
                << d.target_count << ',' << d.point_count << ','
                << d.cell_count << ',' << std::setprecision(10)
                << d.mean_neighbor_count << ',' << d.max_neighbor_count << ','
                << d.padded_neighbor_slots << ',' << d.build_ms << ','
                << d.neighbor_build_ms << ',' << d.exposure_ms << ',' << wall
                << ',' << speedup << '\n';
    }
  }
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
    (Join-Path $RepoRoot 'src\titan_surface.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile SASA surface benchmark fixture."
}

$csv = & $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "SASA surface benchmark fixture failed."
}
Set-Content -LiteralPath $CsvPath -Value $csv -Encoding ASCII
Write-Host "SASA surface benchmark PASS: $CsvPath"
