param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_barnase_precond_truth'),
    [int]$CaseTimeoutSeconds = 180,
    [switch]$ProfileMemory
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$truthScript = Join-Path $PSScriptRoot 'pbe_tabi_barnase_truth_matrix.ps1'
if (-not (Test-Path -LiteralPath $truthScript)) {
    throw "Missing Barnase truth matrix harness: $truthScript"
}

if ($ProfileMemory) {
    & $truthScript -TitanExe $TitanExe -WorkRoot $WorkRoot `
        -CaseTimeoutSeconds $CaseTimeoutSeconds -NearIntegration duffy_adaptive `
        -QuadBuffer adaptive -TrueResidualProbe interval `
        -TrueResidualProbeInterval 4 -PrecondOnly -ProfileMemory
} else {
    & $truthScript -TitanExe $TitanExe -WorkRoot $WorkRoot `
        -CaseTimeoutSeconds $CaseTimeoutSeconds -NearIntegration duffy_adaptive `
        -QuadBuffer adaptive -TrueResidualProbe interval `
        -TrueResidualProbeInterval 4 -PrecondOnly
}
if (-not $?) {
    throw "Barnase preconditioner truth matrix failed."
}

$src = Join-Path $WorkRoot 'pbe_tabi_barnase_truth_matrix.csv'
$dst = Join-Path $WorkRoot 'pbe_tabi_barnase_precond_truth.csv'
if (-not (Test-Path -LiteralPath $src)) {
    throw "Missing truth matrix CSV: $src"
}
Copy-Item -LiteralPath $src -Destination $dst -Force
Write-Host "PBE TABI Barnase preconditioner truth COMPLETE: $dst"
