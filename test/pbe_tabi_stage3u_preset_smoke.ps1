param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_stage3u_preset_smoke')
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
        throw "Refusing to clean path outside Stage 3U preset smoke root: $full"
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

function Require-Value([object]$Row, [string]$Column, [string]$Expected) {
    if (-not ($Row.PSObject.Properties.Name -contains $Column)) {
        throw "Missing column '$Column'"
    }
    if ([string]$Row.$Column -ne $Expected) {
        throw "Expected ${Column}=$Expected, observed '$($Row.$Column)'"
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}

$runDir = Join-Path $WorkRoot 'sphere_stage3u_safe'
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
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
    '-pbe-tabi-tree-preset=stage3u-safe',
    '-pbe-tabi-surface-points=8',
    '-pbe-tabi-expansion-audit-sample-count=64',
    '-pbe-tabi-expansion-audit-max-panels=256',
    '-pbe-tabi-matvec-audit=off',
    '-nt=2'
)
Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII

Push-Location $runDir
try {
    & $TitanExe -a $argPath | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for Stage 3U preset smoke"
    }
} finally {
    Pop-Location
}

$energyPath = Join-Path $runDir 'pbe_energy.csv'
$timingPath = Join-Path $runDir 'pbe_timing.csv'
$profilePath = Join-Path $runDir 'pbe_tabi_profile.csv'
foreach ($required in @($energyPath, $timingPath, $profilePath)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing Stage 3U preset smoke output: $required"
    }
}

$energy = @(Import-Csv -LiteralPath $energyPath)[0]
if ($energy.status -ne 'ok') {
    throw "Stage 3U preset smoke did not finish ok: $($energy.status)"
}

$timing = @(Import-Csv -LiteralPath $timingPath)[0]
$profile = @(Import-Csv -LiteralPath $profilePath)[0]

Require-Value $timing 'tabi_tree_preset' 'stage3u-safe'
Require-Value $timing 'tabi_tree_builder' 'topdown-aabb'
Require-Value $timing 'tabi_node_bounds_mode' 'tight-aabb'
Require-Value $timing 'tabi_mac_safety' 'physical-gap'
Require-Value $timing 'tabi_expansion_center' 'area-centroid'
Require-Value $timing 'tabi_fragmentation_safety' 'split'
Require-Value $timing 'tabi_mac_distance_policy' 'surface-margin'
Require-Value $timing 'tabi_convergence_estimator_model' 'channel-tail-eps'
Require-Value $timing 'tabi_convergence_channel_calibration' 'polynomial-eps'
Require-Value $timing 'tabi_channel_missed_semantics' 'accepted-only'
Require-Value $timing 'tabi_cluster_policy' 'off'
Require-Value $profile 'tree_preset' 'stage3u-safe'

$perfRunDir = Join-Path $WorkRoot 'sphere_stage3u_perf'
New-CleanDir $perfRunDir
$perfPdbPath = Join-Path $perfRunDir 'driver.pdb'
Write-DummyPdb $perfPdbPath
$perfArgPath = Join-Path $perfRunDir 'arg.in'
$perfArgLines = @(
    "-f=$(Convert-ToTitanPath $perfPdbPath)",
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
    '-pbe-tabi-tree-preset=stage3u-perf',
    '-pbe-tabi-surface-points=8',
    '-nt=2'
)
Set-Content -LiteralPath $perfArgPath -Value $perfArgLines -Encoding ASCII

Push-Location $perfRunDir
try {
    & $TitanExe -a $perfArgPath | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for Stage 3U perf preset smoke"
    }
} finally {
    Pop-Location
}

$perfTimingPath = Join-Path $perfRunDir 'pbe_timing.csv'
$perfProfilePath = Join-Path $perfRunDir 'pbe_tabi_profile.csv'
foreach ($required in @($perfTimingPath, $perfProfilePath)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing Stage 3U perf preset smoke output: $required"
    }
}
$perfTiming = @(Import-Csv -LiteralPath $perfTimingPath)[0]
$perfProfile = @(Import-Csv -LiteralPath $perfProfilePath)[0]
Require-Value $perfTiming 'tabi_tree_preset' 'stage3u-perf'
Require-Value $perfTiming 'tabi_cluster_policy' 'leaf-pair-balanced'
Require-Value $perfTiming 'tabi_leaf_panel_kernel' 'tile-avx2'
Require-Value $perfTiming 'tabi_leaf_pair_batch_kernel' 'target-avx2'
Require-Value $perfTiming 'tabi_target_batch_size' '4'
Require-Value $perfProfile 'tree_preset' 'stage3u-perf'

$perf2RunDir = Join-Path $WorkRoot 'sphere_stage3u_perf2'
New-CleanDir $perf2RunDir
$perf2PdbPath = Join-Path $perf2RunDir 'driver.pdb'
Write-DummyPdb $perf2PdbPath
$perf2ArgPath = Join-Path $perf2RunDir 'arg.in'
$perf2ArgLines = @(
    "-f=$(Convert-ToTitanPath $perf2PdbPath)",
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
    '-nt=2'
)
Set-Content -LiteralPath $perf2ArgPath -Value $perf2ArgLines -Encoding ASCII

Push-Location $perf2RunDir
try {
    & $TitanExe -a $perf2ArgPath | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for Stage 3U perf2 preset smoke"
    }
} finally {
    Pop-Location
}

$perf2TimingPath = Join-Path $perf2RunDir 'pbe_timing.csv'
$perf2ProfilePath = Join-Path $perf2RunDir 'pbe_tabi_profile.csv'
foreach ($required in @($perf2TimingPath, $perf2ProfilePath)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing Stage 3U perf2 preset smoke output: $required"
    }
}
$perf2Timing = @(Import-Csv -LiteralPath $perf2TimingPath)[0]
$perf2Profile = @(Import-Csv -LiteralPath $perf2ProfilePath)[0]
Require-Value $perf2Timing 'tabi_tree_preset' 'stage3u-perf2'
Require-Value $perf2Timing 'tabi_cluster_policy' 'leaf-pair-balanced'
Require-Value $perf2Timing 'tabi_leaf_panel_kernel' 'tile-avx2'
Require-Value $perf2Timing 'tabi_leaf_pair_batch_kernel' 'target-avx2'
Require-Value $perf2Timing 'tabi_target_batch_size' '4'
Require-Value $perf2Timing 'tabi_leaf_pair_probe_cache' 'per-matvec'
Require-Value $perf2Timing 'tabi_leaf_pair_target_group' 'same-target-leaf'
Require-Value $perf2Timing 'tabi_leaf_pair_max_lookahead' '32'
Require-Value $perf2Profile 'tree_preset' 'stage3u-perf2'

Write-Host ("PBE TABI Stage 3U preset smoke PASS: safe={0} perf={1} perf2={2} perf2_cache={3} perf2_group={4} surface_margin={5}" -f
    $timing.tabi_tree_preset, $perfTiming.tabi_tree_preset,
    $perf2Timing.tabi_tree_preset,
    $perf2Timing.tabi_leaf_pair_probe_cache,
    $perf2Timing.tabi_leaf_pair_target_group,
    $timing.tabi_mac_surface_margin_factor)
