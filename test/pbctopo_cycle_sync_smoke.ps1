param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$SourceRoot = 'D:\data_MD\vibe\Npro\tetramer_WT',
    [string]$WorkRoot = 'D:\data_MD\vibe\Npro\tetramer_WT\test_pbctopo_cycle_sync_smoke'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$env:TITAN_NO_PAUSE = '1'

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
    $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

if (-not (Test-Path -LiteralPath $SourceRoot)) {
    throw "Missing source root: $SourceRoot"
}
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$arg = Join-Path $WorkRoot 'arg.in'
$argsText = @"
-ff=charmm36
-s
$SourceRoot\md1_snap0.pdb
-traj
$SourceRoot\md1_1.xtc
-pbctopo
-pbctopo-exact-component-limit=6
-pbctopo-beam-width=128
-pbctopo-anchor-backend=auto
-g=0-0
-nt=4
-eu
kCal/mol
-1pdb
-nocal
-dt
400
-b
0
-e
0
-n
$SourceRoot\index.ndx
-common_topo
$SourceRoot\topol.top
-noprint-side_chain_contact
-no-dssp
-opdb=pbctopo_cycle_sync_smoke.pdb
"@
[System.IO.File]::WriteAllText($arg, $argsText, [System.Text.Encoding]::ASCII)

Push-Location $WorkRoot
try {
    & $titanPath -a .\arg.in
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}

$debugLog = Join-Path $WorkRoot 'pbctopo_debug.log'
$report = Join-Path $WorkRoot 'pbctopo_frame_report.csv'
if (-not (Test-Path -LiteralPath $debugLog)) {
    throw "Missing pbctopo debug log"
}
if (-not (Test-Path -LiteralPath $report)) {
    throw "Missing pbctopo frame report"
}

$logText = Get-Content -Raw -LiteralPath $debugLog
foreach ($pattern in @(
        'component_cycle_sync_attempts=\d+',
        'component_cycle_sync_used=\d+',
        'component_cycle_sync_pair_candidates=\d+',
        'component_cycle_sync_cycle_support=\d+',
        'component_cycle_sync_cycle_mismatch=([0-9.]+|nan|na)',
        'component_cycle_sync_inconsistent_edges=\d+',
        'component_cycle_sync=([0-9.]+) s'
    )) {
    if ($logText -notmatch $pattern) {
        throw "Missing expected cycle-sync log pattern: $pattern"
    }
}

$rows = @(Import-Csv -LiteralPath $report)
if ($rows.Count -lt 1) {
    throw "Expected at least one frame report row"
}
foreach ($column in @('cycle_support', 'cycle_mismatch', 'inconsistent_edge_count')) {
    if (-not ($rows[0].PSObject.Properties.Name -contains $column)) {
        throw "Missing cycle column in frame report: $column"
    }
}

Write-Host 'PASS pbctopo_cycle_sync_smoke'
