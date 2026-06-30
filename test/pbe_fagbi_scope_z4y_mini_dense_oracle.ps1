param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4y_mini_dense'),
    [string[]]$Fixtures = @('1vii-mini','1ubq-mini','2ci2-mini'),
    [int[]]$SurfacePoints = @(8,12),
    [int]$MiniAtomCap = 64,
    [int]$CaseTimeoutSeconds = 1800,
    [int]$Threads = 4,
    [switch]$PlanOnly,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$matrixScript = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4y_reviewer_hardening_matrix.ps1'
if (-not (Test-Path -LiteralPath $matrixScript)) {
    throw "Missing Z4Y reviewer-hardening harness: $matrixScript"
}

$matrixArgs = @{
    TitanExe = $TitanExe
    RepoRoot = $RepoRoot
    WorkRoot = $WorkRoot
    Stages = @('dense_oracle')
    Fixtures = $Fixtures
    SurfacePoints = $SurfacePoints
    MiniAtomCap = $MiniAtomCap
    DenseTimeoutSeconds = $CaseTimeoutSeconds
    FmmTimeoutSeconds = $CaseTimeoutSeconds
    Threads = $Threads
    DenseOracleProteinOnly = $true
}
if ($PlanOnly) { $matrixArgs.PlanOnly = $true }
if ($Force) { $matrixArgs.Force = $true }
if ($RerunFailed) { $matrixArgs.RerunFailed = $true }
if ($AllowFailedDecision) { $matrixArgs.AllowFailedDecision = $true }

& $matrixScript @matrixArgs
