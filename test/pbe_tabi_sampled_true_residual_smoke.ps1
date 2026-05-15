param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_sampled_true_residual_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI sampled smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
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

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot

$pqrPath = Join-Path $InputRoot 'sphere.pqr'
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
    '-pbe-max-iter=160',
    '-pbe-rms=1e-8',
    '-pbe-dphi=1e-8',
    '-pbe-timing=detail',
    '-pbe-tabi-progress=auto',
    '-pbe-tabi-progress-interval=1',
    '-pbe-tabi-sum=direct',
    '-pbe-tabi-formulation=jcp2013',
    '-pbe-tabi-precond=diagonal',
    '-pbe-tabi-near-integration=centroid',
    '-pbe-tabi-surface-points=24',
    '-pbe-tabi-gmres-restart=8',
    '-pbe-tabi-true-residual-probe=interval',
    '-pbe-tabi-true-residual-probe-scope=sampled',
    '-pbe-tabi-true-residual-probe-interval=1',
    '-pbe-tabi-true-residual-sample-count=128',
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
$profile = @(Import-Csv -LiteralPath (Join-Path $WorkRoot 'pbe_tabi_profile.csv'))[0]
$progress = @(Import-Csv -LiteralPath (Join-Path $WorkRoot 'pbe_tabi_longrun_progress.csv'))

foreach ($col in @('tabi_true_residual_probe_scope',
                   'tabi_true_residual_sample_count',
                   'tabi_sampled_true_residual_rel',
                   'tabi_sampled_true_residual_rows',
                   'tabi_sampled_probe_calls',
                   'tabi_sampled_probe_seconds',
                   'tabi_full_true_residual_confirm_calls',
                   'tabi_full_probe_seconds')) {
    if (-not ($timing.PSObject.Properties.Name -contains $col)) {
        throw "Missing timing sampled-probe column $col"
    }
}
foreach ($col in @('true_residual_probe_scope',
                   'true_residual_sample_count',
                   'sampled_true_residual_rel',
                   'sampled_true_residual_rows',
                   'sampled_probe_calls',
                   'sampled_probe_seconds',
                   'full_true_residual_confirm_calls',
                   'full_probe_seconds')) {
    if (-not ($profile.PSObject.Properties.Name -contains $col)) {
        throw "Missing profile sampled-probe column $col"
    }
    if (-not ($progress[0].PSObject.Properties.Name -contains $col)) {
        throw "Missing progress sampled-probe column $col"
    }
}
if ([string]$timing.tabi_true_residual_probe_mode -ne 'interval') {
    throw "Expected interval probe mode, got $($timing.tabi_true_residual_probe_mode)"
}
if ([string]$timing.tabi_true_residual_probe_scope -ne 'sampled') {
    throw "Expected sampled probe scope, got $($timing.tabi_true_residual_probe_scope)"
}
if ([int]$timing.tabi_true_residual_sample_count -ne 128) {
    throw "Expected sample count 128, got $($timing.tabi_true_residual_sample_count)"
}
$sampleCalls = [int64]$timing.tabi_sampled_probe_calls
if ($sampleCalls -le 0) {
    throw "Expected at least one sampled true-residual probe call."
}
$sampleRows = [int64]$timing.tabi_sampled_true_residual_rows
if ($sampleRows -le 0) {
    throw "Expected sampled residual rows."
}
$sampleSeconds = Convert-ToDouble $timing.tabi_sampled_probe_seconds
if ([double]::IsNaN($sampleSeconds) -or [double]::IsInfinity($sampleSeconds) -or $sampleSeconds -lt 0.0) {
    throw "Sampled probe seconds must be finite and non-negative."
}
$sampleRowsProgress = @($progress | Where-Object { $_.phase -eq 'gmres_sampled_true_probe' })
if ($sampleRowsProgress.Count -le 0) {
    throw "Progress CSV did not record gmres_sampled_true_probe rows."
}
if ($progress[-1].phase -ne 'solve_done') {
    throw "Progress CSV did not finish with solve_done."
}

Write-Host "PBE TABI sampled true-residual smoke PASS: $WorkRoot"
