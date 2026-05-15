param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_longrun_progress_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
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
        throw "Refusing to clean path outside TABI long-run progress smoke root: $full"
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

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}

$runDir = Join-Path $WorkRoot 'mini_progress'
New-CleanDir $runDir

$pdbPath = Join-Path $runDir 'driver.pdb'
Write-DummyPdb $pdbPath
$pqrPath = Join-Path $InputRoot 'mini-protein.pqr'
if (-not (Test-Path -LiteralPath $pqrPath)) {
    throw "Missing mini-protein PQR fixture: $pqrPath"
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
    '-pbe-max-iter=800',
    '-pbe-rms=1e-4',
    '-pbe-dphi=1e-4',
    '-pbe-timing=detail',
    '-pbe-tabi-progress=auto',
    '-pbe-tabi-progress-interval=1',
    '-pbe-tabi-sum=direct',
    '-pbe-tabi-formulation=jcp2013',
    '-pbe-tabi-near-integration=centroid',
    '-nt=1'
)
Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII

Push-Location $runDir
try {
    & $TitanExe -a $argPath | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for TABI long-run progress smoke"
    }
} finally {
    Pop-Location
}

$progressPath = Join-Path $runDir 'pbe_tabi_longrun_progress.csv'
if (-not (Test-Path -LiteralPath $progressPath)) {
    throw "Missing TABI long-run progress CSV: $progressPath"
}
$rows = @(Import-Csv -LiteralPath $progressPath)
if ($rows.Count -lt 4) {
    throw "Expected at least four TABI progress rows, found $($rows.Count)"
}
$phases = @($rows | ForEach-Object { $_.phase })
foreach ($phase in @('surface_start','surface_end','gmres_start','solve_done')) {
    if ($phases -notcontains $phase) {
        throw "Missing TABI progress phase: $phase"
    }
}
$previous = -1.0
foreach ($row in $rows) {
    $elapsed = [double]::Parse($row.elapsed_seconds, $InvariantCulture)
    if ($elapsed -lt $previous) {
        throw "TABI progress elapsed time is not monotonic."
    }
    $previous = $elapsed
}
$last = $rows[-1]
if ($last.phase -ne 'solve_done' -or $last.status -ne 'ok') {
    throw "Unexpected final TABI progress row: phase=$($last.phase) status=$($last.status)"
}

Write-Host "PBE TABI long-run progress smoke PASS: $progressPath"
