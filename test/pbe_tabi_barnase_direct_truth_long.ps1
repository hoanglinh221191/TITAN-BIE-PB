param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_barnase_direct_truth_long'),
    [int]$CaseTimeoutSeconds = 300,
    [string[]]$DirectPreconditioners = @('none','diagonal'),
    [int[]]$DirectRestarts = @(64,100),
    [switch]$ProfileMemory
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$truthScript = Join-Path $PSScriptRoot 'pbe_tabi_barnase_truth_matrix.ps1'
if (-not (Test-Path -LiteralPath $truthScript)) {
    throw "Missing Barnase truth matrix harness: $truthScript"
}

$common = @{
    TitanExe = $TitanExe
    WorkRoot = $WorkRoot
    CaseTimeoutSeconds = $CaseTimeoutSeconds
    NearIntegration = 'duffy_adaptive'
    QuadBuffer = 'adaptive'
    TrueResidualProbe = 'interval'
    TrueResidualProbeScope = 'adaptive'
    TrueResidualProbeInterval = 4
    TrueResidualSampleCount = 1024
    DirectPreconditioners = $DirectPreconditioners
    DirectRestarts = $DirectRestarts
    PrecondOnly = $true
}
if ($ProfileMemory) {
    $common.ProfileMemory = $true
}

& $truthScript @common
if (-not $?) {
    throw "Barnase direct truth long matrix failed."
}

$src = Join-Path $WorkRoot 'pbe_tabi_barnase_truth_matrix.csv'
$dst = Join-Path $WorkRoot 'pbe_tabi_barnase_direct_truth_long.csv'
if (-not (Test-Path -LiteralPath $src)) {
    throw "Missing truth matrix CSV: $src"
}
Copy-Item -LiteralPath $src -Destination $dst -Force
Write-Host "PBE TABI Barnase direct truth long COMPLETE: $dst"
