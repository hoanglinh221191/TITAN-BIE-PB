param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_matvec_cost_profile_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside matvec-cost smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    $lines = @(
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Require-Column([object]$Row, [string]$Column) {
    if (-not ($Row.PSObject.Properties.Name -contains $Column)) {
        throw "Missing column '$Column'"
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$runDir = Join-Path $WorkRoot 'sphere_matvec_cost'
New-CleanDir $runDir

$pdbPath = Join-Path $runDir 'driver.pdb'
Write-DummyPdb $pdbPath
$pqrPath = Join-Path $InputRoot 'sphere.pqr'
if (-not (Test-Path -LiteralPath $pqrPath)) {
    throw "Missing sphere PQR fixture: $pqrPath"
}
$libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
$argPath = Join-Path $runDir 'arg.in'
$argLines = @(
    "-f=$(Convert-ToTitanPath $pdbPath)",
    '-1pdb',
    '-ff=AMBER99SB',
    '-g=1&1',
    "-path-a=$libPath",
    "-path-n=$libPath",
    '-pbe=tabi',
    "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
    '-pbe-surface=vdw',
    '-pbe-scale=0.35',
    '-pbe-margin=2.0',
    '-pbe-max-grid=300000',
    '-pbe-max-iter=120',
    '-pbe-rms=1e-4',
    '-pbe-dphi=1e-4',
    '-pbe-timing=detail',
    '-pbe-tabi-tree-preset=stage3u-perf2',
    '-pbe-tabi-surface-points=8',
    '-pbe-tabi-matvec-cost-profile=summary',
    '-pbe-tabi-leaf-kernel-profile=summary',
    '-nt=2'
)
Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII

Push-Location $runDir
try {
    & $TitanExe -a $argPath | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for matvec-cost profile smoke"
    }
} finally {
    Pop-Location
}

$energyPath = Join-Path $runDir 'pbe_energy.csv'
$timingPath = Join-Path $runDir 'pbe_timing.csv'
$profilePath = Join-Path $runDir 'pbe_tabi_profile.csv'
foreach ($required in @($energyPath, $timingPath, $profilePath)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing matvec-cost smoke output: $required"
    }
}

$energy = @(Import-Csv -LiteralPath $energyPath)[0]
if ($energy.status -ne 'ok') {
    throw "Matvec-cost smoke did not finish ok: $($energy.status)"
}

$timing = @(Import-Csv -LiteralPath $timingPath)[0]
$profile = @(Import-Csv -LiteralPath $profilePath)[0]
foreach ($column in @(
    'tabi_matvec_cost_profile',
    'tabi_leaf_panels_per_second',
    'tabi_clean_leaf_panels_per_second',
    'tabi_tree_traversal_fraction',
    'tabi_matvec_calls',
    'tabi_tree_traversal_seconds',
    'tabi_preconditioner_seconds',
    'tabi_leaf_panel_interactions',
    'tabi_leaf_pair_batch_source_panels')) {
    Require-Column $timing $column
}
foreach ($column in @(
    'matvec_cost_profile',
    'leaf_panels_per_second',
    'clean_leaf_panels_per_second',
    'tree_traversal_fraction',
    'matvec_calls',
    'tree_traversal_seconds',
    'preconditioner_seconds',
    'leaf_panel_interactions',
    'leaf_pair_batch_source_panels')) {
    Require-Column $profile $column
}

if ([string]$timing.tabi_matvec_cost_profile -ne 'summary') {
    throw "Timing did not report matvec-cost profile summary: $($timing.tabi_matvec_cost_profile)"
}
if ([string]$profile.matvec_cost_profile -ne 'summary') {
    throw "Profile did not report matvec-cost profile summary: $($profile.matvec_cost_profile)"
}
if ([double]$timing.tabi_matvec_calls -le 0) {
    throw "Matvec-cost smoke did not record matvec calls."
}
if ([double]$timing.tabi_tree_traversal_fraction -lt 0.0 -or
    [double]$timing.tabi_tree_traversal_fraction -gt 2.0) {
    throw "Unexpected tree traversal fraction: $($timing.tabi_tree_traversal_fraction)"
}

Write-Host "PBE TABI matvec-cost profile smoke PASS matvecs=$($timing.tabi_matvec_calls) leaf_panels_per_s=$($timing.tabi_leaf_panels_per_second) traversal_fraction=$($timing.tabi_tree_traversal_fraction)"
