param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_block_precond_buffer_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_block_precond_buffer_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_block_precond_buffer_smoke.exe'
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI block preconditioner smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

if (-not (Test-Path -LiteralPath $Compiler)) { throw "Missing compiler: $Compiler" }
if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot

$source = @'
#include "tabi/tabi_near.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using TitanPBE::Tabi::Panel;

Panel make_panel(double z, std::size_t offset) {
  Panel p;
  p.v0x = 0.0; p.v0y = 0.0; p.v0z = z;
  p.v1x = 1.0; p.v1y = 0.0; p.v1z = z;
  p.v2x = 0.0; p.v2y = 1.0; p.v2z = z;
  p.v0_index = offset + 0;
  p.v1_index = offset + 1;
  p.v2_index = offset + 2;
  p.x = (p.v0x + p.v1x + p.v2x) / 3.0;
  p.y = (p.v0y + p.v1y + p.v2y) / 3.0;
  p.z = z;
  p.nx = 0.0; p.ny = 0.0; p.nz = 1.0;
  p.area = 0.5;
  p.has_triangle_geometry = true;
  return p;
}

int main() {
  const Panel target = make_panel(0.0, 0);
  const Panel source = make_panel(0.3, 10);
  auto decision = TitanPBE::Tabi::classify_geometric_near_buffer_pair(
      target, source, "adaptive");
  if (!decision.active)
    throw std::runtime_error("expected geometric buffer pair");
  TitanPBE::Tabi::TreecodeStats stats;
  auto kernel = TitanPBE::Tabi::evaluate_jcp_near_pair_kernel(
      target, 0, source, 1, 80.0, 0.125, true, "duffy_adaptive",
      "adaptive", &stats);
  if (!kernel.valid || !kernel.geometric_buffer)
    throw std::runtime_error("expected valid geometric-buffer kernel");
  if (!std::isfinite(kernel.k1) || !std::isfinite(kernel.k2) ||
      !std::isfinite(kernel.k3) || !std::isfinite(kernel.k4))
    throw std::runtime_error("non-finite geometric-buffer kernel");
  if (stats.near_geometric_buffer_pairs == 0 ||
      stats.near_geometric_buffer_adaptive_pairs == 0)
    throw std::runtime_error("missing geometric-buffer audit counters");
  std::cout << "geometric buffer kernel ok\n";
  return 0;
}
'@
Set-Content -LiteralPath $SourcePath -Value $source -Encoding ASCII
& $Compiler -std=c++23 $("-I" + (Join-Path $RepoRoot "src")) $SourcePath (Join-Path $RepoRoot 'src\tabi\tabi_near.cpp') -o $ExePath
if ($LASTEXITCODE -ne 0) { throw "Failed to compile geometric-buffer near kernel smoke." }
& $ExePath | Out-Null
if ($LASTEXITCODE -ne 0) { throw "Geometric-buffer near kernel smoke failed." }

$inputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'
$pqrPath = Join-Path $inputRoot 'sphere.pqr'
if (-not (Test-Path -LiteralPath $pqrPath)) { throw "Missing sphere PQR fixture: $pqrPath" }
$pdbPath = Join-Path $WorkRoot 'driver.pdb'
Write-DummyPdb $pdbPath
$libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
$argPath = Join-Path $WorkRoot 'arg.in'
Set-Content -LiteralPath $argPath -Value @(
    "-f=$(Convert-ToTitanPath $pdbPath)",
    '-1pdb',
    '-ff=AMBER99SB',
    '-g=1&1',
    "-path-a=$libPath",
    "-path-n=$libPath",
    '-pbe=tabi',
    "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
    '-pbe-surface=vdw',
    '-pbe-salt=0.15',
    '-pbe-max-iter=80',
    '-pbe-rms=1e-4',
    '-pbe-dphi=1e-4',
    '-pbe-timing=detail',
    '-pbe-tabi-sum=direct',
    '-pbe-tabi-formulation=jcp2013',
    '-pbe-tabi-precond=block',
    '-pbe-tabi-precond-block-near=buffer-aware',
    '-pbe-tabi-near-integration=duffy_adaptive',
    '-pbe-tabi-quad-buffer=adaptive',
    '-pbe-tabi-surface-points=24',
    '-pbe-tabi-gmres-restart=8',
    '-nt=2'
) -Encoding ASCII

Push-Location $WorkRoot
try {
    & $TitanExe -a $argPath | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE" }
} finally {
    Pop-Location
}

$timing = @(Import-Csv -LiteralPath (Join-Path $WorkRoot 'pbe_timing.csv'))[0]
foreach ($col in @('tabi_preconditioner_block_near','tabi_preconditioner_effective')) {
    if (-not ($timing.PSObject.Properties.Name -contains $col)) {
        throw "Missing timing column $col"
    }
}
if ([string]$timing.tabi_preconditioner_block_near -ne 'buffer-aware') {
    throw "Expected buffer-aware block-near timing column."
}
Write-Host "PBE TABI block preconditioner buffer smoke PASS: $WorkRoot"
