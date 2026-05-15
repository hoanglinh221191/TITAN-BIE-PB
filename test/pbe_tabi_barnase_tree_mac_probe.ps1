param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_barnase_tree_mac_probe'),
    [int]$CaseTimeoutSeconds = 300,
    [string]$TrueResidualProbeScope = 'adaptive',
    [int]$TrueResidualProbeInterval = 4,
    [int]$TrueResidualSampleCount = 1024,
    [string[]]$TreeInteractionModes = @('onthefly'),
    [int[]]$TreeLeafSizes = @(64,128,256),
    [double[]]$TreeOpeningThetas = @(0.5,0.7,1.0,1.2),
    [string[]]$NodeBoundsModes = @('moment','tight-aabb'),
    [string[]]$TreeBuilders = @('grouped-morton','topdown-aabb'),
    [string[]]$MacSafety = @('off'),
    [double]$MacGapFactor = 1.5,
    [string]$MomentAudit = 'off',
    [switch]$ProfileMemory
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$MacSafety = @($MacSafety | ForEach-Object {
    ([string]$_ -split ',') | Where-Object { -not [string]::IsNullOrWhiteSpace($_) }
})

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
    TrueResidualProbeScope = $TrueResidualProbeScope
    TrueResidualProbeInterval = $TrueResidualProbeInterval
    TrueResidualSampleCount = $TrueResidualSampleCount
    TreeProbeOnly = $true
    TreeLeafSizes = $TreeLeafSizes
    TreeOpeningThetas = $TreeOpeningThetas
    NodeBoundsModes = $NodeBoundsModes
    TreeBuilders = $TreeBuilders
    TreeInteractionModes = $TreeInteractionModes
    MacSafety = $MacSafety
    MacGapFactor = $MacGapFactor
    MomentAudit = $MomentAudit
}
if ($ProfileMemory) {
    $common.ProfileMemory = $true
}

& $truthScript @common
if (-not $?) {
    throw "Barnase tree MAC probe failed."
}

$src = Join-Path $WorkRoot 'pbe_tabi_barnase_truth_matrix.csv'
$dst = Join-Path $WorkRoot 'pbe_tabi_barnase_tree_mac_probe.csv'
if (-not (Test-Path -LiteralPath $src)) {
    throw "Missing tree MAC CSV: $src"
}
Copy-Item -LiteralPath $src -Destination $dst -Force
Write-Host "PBE TABI Barnase tree MAC probe COMPLETE: $dst"
