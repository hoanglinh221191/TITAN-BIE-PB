param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_coarse_galerkin_precond_smoke'),
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
        throw "Refusing to clean path outside TABI coarse Galerkin smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

function Invoke-CoarseCase(
    [string]$RunDir,
    [string]$PqrPath,
    [int]$SurfacePoints,
    [string]$CoarseBase,
    [int]$MemoryBudgetMb,
    [int]$TreeLeafSize = 8,
    [int]$MaxIter = 30,
    [double]$RmsTolerance = 1.0e-4
) {
    New-CleanDir $RunDir
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-ToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-salt=0.15',
        "-pbe-max-iter=$MaxIter",
        "-pbe-rms=$($RmsTolerance.ToString('G17', [System.Globalization.CultureInfo]::InvariantCulture))",
        '-pbe-dphi=1e-4',
        '-pbe-timing=detail',
        '-pbe-tabi-tree-preset=stage3u-safe',
        "-pbe-tabi-surface-points=$SurfacePoints",
        "-pbe-tabi-tree-leaf-size=$TreeLeafSize",
        '-pbe-tabi-precond=coarse-galerkin-probe',
        '-pbe-tabi-precond-coarse-aggregate-cap=8',
        "-pbe-tabi-precond-coarse-base=$CoarseBase",
        '-pbe-tabi-precond-coarse-verify=full-small',
        "-pbe-tabi-precond-memory-budget-mb=$MemoryBudgetMb",
        '-pbe-tabi-precond-block-near=buffer-aware',
        '-pbe-tabi-precond-block-panel-cap=32',
        '-pbe-tabi-gmres-restart=8',
        '-pbe-tabi-progress=detail',
        '-pbe-tabi-progress-interval=1',
        '-nt=2'
    ) -Encoding ASCII
    Push-Location $RunDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE" }
    } finally {
        Pop-Location
    }
    return @(Import-Csv -LiteralPath (Join-Path $RunDir 'pbe_timing.csv'))[0]
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot
$spherePqr = Join-Path $PSScriptRoot 'pbe_inputs_local\sphere.pqr'
if (-not (Test-Path -LiteralPath $spherePqr)) { throw "Missing sphere PQR fixture: $spherePqr" }

$activeDir = Join-Path $WorkRoot 'active'
$timing = Invoke-CoarseCase -RunDir $activeDir -PqrPath $spherePqr -SurfacePoints 48 -CoarseBase 'diagonal' -MemoryBudgetMb 256
foreach ($col in @(
    'tabi_preconditioner_coarse_aggregate_cap',
    'tabi_preconditioner_coarse_aggregates',
    'tabi_preconditioner_coarse_dimension',
    'tabi_preconditioner_coarse_base',
    'tabi_preconditioner_coarse_matrix_bytes',
    'tabi_preconditioner_coarse_pivot_ratio',
    'tabi_preconditioner_coarse_status',
    'tabi_preconditioner_coarse_setup_matvec_calls',
    'tabi_preconditioner_coarse_setup_seconds',
    'tabi_preconditioner_coarse_apply_seconds',
    'tabi_preconditioner_coarse_apply_calls',
    'tabi_preconditioner_coarse_rp_identity_error_max',
    'tabi_preconditioner_coarse_verify_rel_error_max'
)) {
    if (-not ($timing.PSObject.Properties.Name -contains $col)) {
        throw "Missing timing column $col"
    }
}
if ([string]$timing.tabi_preconditioner_requested -ne 'coarse-galerkin-probe') {
    throw 'Expected coarse Galerkin preconditioner request.'
}
if ([string]$timing.tabi_preconditioner_effective -ne 'coarse_galerkin_probe_diagonal') {
    throw 'Expected active diagonal-base coarse Galerkin preconditioner.'
}
if ([string]$timing.tabi_preconditioner_coarse_status -ne 'active') {
    throw 'Expected active coarse status.'
}
if ([int]$timing.tabi_preconditioner_coarse_aggregates -le 0 -or
    [int]$timing.tabi_preconditioner_coarse_dimension -le 0) {
    throw 'Expected a non-empty coarse space.'
}
if ([int]$timing.tabi_preconditioner_coarse_dimension -gt 16) {
    throw 'Coarse dimension exceeded the configured coupled aggregate cap.'
}
if ([double]$timing.tabi_preconditioner_coarse_matrix_bytes -le 0 -or
    [double]$timing.tabi_preconditioner_coarse_setup_matvec_calls -lt [int]$timing.tabi_preconditioner_coarse_dimension) {
    throw 'Expected exact Galerkin setup work.'
}
if ([double]$timing.tabi_preconditioner_coarse_rp_identity_error_max -gt 1.0e-12 -or
    [double]$timing.tabi_preconditioner_coarse_verify_rel_error_max -gt 1.0e-12) {
    throw 'Expected exact coarse identity/Galerkin verification.'
}
if ([double]$timing.tabi_preconditioner_coarse_apply_calls -le 0) {
    throw 'Expected coarse preconditioner applications.'
}

$progress = @(Import-Csv -LiteralPath (Join-Path $activeDir 'pbe_tabi_longrun_progress.csv'))[-1]
foreach ($col in @(
    'preconditioner_coarse_status',
    'preconditioner_coarse_dimension',
    'preconditioner_coarse_setup_matvec_calls',
    'preconditioner_coarse_apply_calls'
)) {
    if (-not ($progress.PSObject.Properties.Name -contains $col)) {
        throw "Missing progress column $col"
    }
}

if ($IncludeBudgetFallback) {
    $compactPqr = Join-Path $PSScriptRoot 'pbe_inputs_local\compact-mini-protein.pqr'
    if (-not (Test-Path -LiteralPath $compactPqr)) { throw "Missing compact PQR fixture: $compactPqr" }
    $fallback = Invoke-CoarseCase -RunDir (Join-Path $WorkRoot 'budget_fallback') -PqrPath $compactPqr -SurfacePoints 256 -CoarseBase 'block-budgeted' -MemoryBudgetMb 16 -TreeLeafSize 128 -MaxIter 1 -RmsTolerance 2.0
    if ([string]$fallback.tabi_preconditioner_coarse_status -ne 'budget_fallback' -or
        [string]$fallback.tabi_preconditioner_effective -ne 'diagonal_coarse_budget_cap') {
        throw 'Expected coarse plus block base to fallback before allocation when over budget.'
    }
}

Write-Host "PBE TABI coarse Galerkin preconditioner smoke PASS: $WorkRoot"
