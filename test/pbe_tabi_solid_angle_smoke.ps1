param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_solid_angle_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_solid_angle_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_solid_angle_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI solid-angle smoke root: $full"
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
#include "tabi/tabi_types.h"

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

TitanPBE::Tabi::Result run_mode(const std::string &mode) {
  TitanPBE::Tabi::Options opt;
  opt.surface_mode = "vdw";
  opt.surface_geometry = "triangle";
  opt.diagonal_mode = mode;
  TitanPBE::Tabi::Result result;
  std::vector<TitanPBE::Tabi::Panel> panels =
      TitanPBE::Tabi::build_surface(
          {{0.0, 0.0, 0.0, 1.0, 1.6},
           {1.9, 0.1, 0.0, -1.0, 1.6}},
          opt, &result);
  require_true(mode + " panels", !panels.empty());
  return result;
}

} // namespace

int main() {
  const auto solid = run_mode("solid-angle");
  require_true("solid geometry", solid.surface_geometry == "triangle");
  require_true("solid diag min", solid.diag_min > 0.25);
  require_true("solid diag max", solid.diag_max <= 1.75);
  require_true("solid fallbacks", solid.solid_angle_fallbacks < 4);
  require_true("solid has curvature signal", solid.diag_min < 0.999);

  const auto flat = run_mode("flat");
  require_true("flat diag min", std::abs(flat.diag_min - 1.0) < 1.0e-12);
  require_true("flat diag max", std::abs(flat.diag_max - 1.0) < 1.0e-12);

  const auto off = run_mode("off");
  require_true("off diag min", std::abs(off.diag_min) < 1.0e-12);
  require_true("off diag max", std::abs(off.diag_max) < 1.0e-12);

  std::cout << "PBE TABI solid-angle smoke PASS diag=["
            << solid.diag_min << "," << solid.diag_max
            << "] fallbacks=" << solid.solid_angle_fallbacks << "\n";
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
    (Join-Path $RepoRoot 'src\tabi\tabi_surface.cpp'),
    (Join-Path $RepoRoot 'src\titan_surface.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI solid-angle smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI solid-angle smoke fixture failed."
}
