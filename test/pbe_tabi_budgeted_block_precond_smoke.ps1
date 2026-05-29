param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_budgeted_block_precond_smoke'),
    [switch]$IncludeBudgetFallback
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI budgeted preconditioner smoke root: $full"
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

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot

$pqrPath = Join-Path $PSScriptRoot 'pbe_inputs_local\sphere.pqr'
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
    '-pbe-tabi-precond=block-budgeted',
    '-pbe-tabi-precond-block-near=buffer-aware',
    '-pbe-tabi-precond-block-panel-cap=32',
    '-pbe-tabi-precond-memory-budget-mb=256',
    '-pbe-tabi-near-integration=duffy_adaptive',
    '-pbe-tabi-quad-buffer=adaptive',
    '-pbe-tabi-surface-points=24',
    '-pbe-tabi-gmres-restart=8',
    '-pbe-tabi-progress=detail',
    '-pbe-tabi-progress-interval=1',
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
foreach ($col in @(
    'tabi_preconditioner_requested',
    'tabi_preconditioner_effective',
    'tabi_preconditioner_block_near',
    'tabi_preconditioner_block_panel_cap',
    'tabi_preconditioner_memory_budget_bytes',
    'tabi_preconditioner_estimated_bytes',
    'tabi_preconditioner_valid_blocks',
    'tabi_preconditioner_valid_panels',
    'tabi_preconditioner_ill_conditioned_blocks',
    'tabi_preconditioner_setup_seconds',
    'tabi_preconditioner_apply_seconds',
    'tabi_preconditioner_apply_calls'
)) {
    if (-not ($timing.PSObject.Properties.Name -contains $col)) {
        throw "Missing timing column $col"
    }
}
if ([string]$timing.tabi_preconditioner_requested -ne 'block-budgeted') {
    throw "Expected budgeted block preconditioner request."
}
if ([string]$timing.tabi_preconditioner_effective -ne 'block_budgeted_buffer_aware') {
    throw "Expected active budgeted buffer-aware block preconditioner."
}
if ([int]$timing.tabi_preconditioner_block_panel_cap -ne 32) {
    throw "Expected configured block panel cap."
}
if ([double]$timing.tabi_preconditioner_memory_budget_bytes -lt (256 * 1024 * 1024)) {
    throw "Expected configured block memory budget."
}
if ([double]$timing.tabi_preconditioner_estimated_bytes -le 0) {
    throw "Expected estimated preconditioner bytes."
}
if ([double]$timing.tabi_preconditioner_apply_calls -le 0) {
    throw "Expected preconditioner applications."
}

$progressPath = Join-Path $WorkRoot 'pbe_tabi_longrun_progress.csv'
if (-not (Test-Path -LiteralPath $progressPath)) {
    throw "Expected timeout-safe progress output: $progressPath"
}
$progress = @(Import-Csv -LiteralPath $progressPath)[-1]
foreach ($col in @(
    'preconditioner_effective',
    'preconditioner_block_panel_cap',
    'preconditioner_memory_budget_bytes',
    'preconditioner_estimated_bytes',
    'preconditioner_valid_blocks',
    'preconditioner_valid_panels',
    'preconditioner_setup_seconds',
    'preconditioner_apply_seconds',
    'preconditioner_apply_calls'
)) {
    if (-not ($progress.PSObject.Properties.Name -contains $col)) {
        throw "Missing progress column $col"
    }
}
if ([string]$progress.preconditioner_effective -ne 'block_budgeted_buffer_aware') {
    throw "Expected budgeted preconditioner in timeout-safe progress output."
}

if ($IncludeBudgetFallback) {
    $fallbackDir = Join-Path $WorkRoot 'budget_fallback'
    New-CleanDir $fallbackDir
    $fallbackPdb = Join-Path $fallbackDir 'driver.pdb'
    Write-DummyPdb $fallbackPdb
    $fallbackPqr = Join-Path $PSScriptRoot 'pbe_inputs_local\compact-mini-protein.pqr'
    if (-not (Test-Path -LiteralPath $fallbackPqr)) { throw "Missing compact PQR fixture: $fallbackPqr" }
    $fallbackArg = Join-Path $fallbackDir 'arg.in'
    Set-Content -LiteralPath $fallbackArg -Value @(
        "-f=$(Convert-ToTitanPath $fallbackPdb)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-ToTitanPath $fallbackPqr)",
        '-pbe-surface=vdw',
        '-pbe-salt=0.15',
        '-pbe-max-iter=1',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-timing=detail',
        '-pbe-tabi-sum=direct',
        '-pbe-tabi-formulation=jcp2013',
        '-pbe-tabi-precond=block-budgeted',
        '-pbe-tabi-precond-block-near=buffer-aware',
        '-pbe-tabi-precond-block-panel-cap=32',
        '-pbe-tabi-precond-memory-budget-mb=16',
        '-pbe-tabi-near-integration=duffy_adaptive',
        '-pbe-tabi-quad-buffer=adaptive',
        '-pbe-tabi-surface-points=256',
        '-pbe-tabi-gmres-restart=8',
        '-nt=2'
    ) -Encoding ASCII
    Push-Location $fallbackDir
    try {
        & $TitanExe -a $fallbackArg | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "Budget fallback run exited with code $LASTEXITCODE" }
    } finally {
        Pop-Location
    }
    $fallbackTiming = @(Import-Csv -LiteralPath (Join-Path $fallbackDir 'pbe_timing.csv'))[0]
    if ([string]$fallbackTiming.tabi_preconditioner_effective -ne 'diagonal_block_budget_cap') {
        throw "Expected diagonal_block_budget_cap fallback."
    }
    if ([double]$fallbackTiming.tabi_preconditioner_estimated_bytes -le [double]$fallbackTiming.tabi_preconditioner_memory_budget_bytes) {
        throw "Expected estimate to exceed configured budget."
    }
    if ([double]$fallbackTiming.tabi_preconditioner_valid_blocks -ne 0) {
        throw "Budget fallback must occur before LU block allocation."
    }
}

Write-Host "PBE TABI budgeted block preconditioner smoke PASS: $WorkRoot"
