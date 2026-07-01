param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$SourceRoot = 'D:\Data_MD\Npro\tetramer_conf0_unphos',
    [string]$WorkRoot = '',
    [string]$BaselineSha256 = '785E6BD8AE81395BA77E5A70653753771AFBAC94E916DBE3EADFAC7B8C69DACD'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$requiredInputs = @('md1_snap0.pdb', 'md1_1.xtc', 'index.ndx', 'topol.top')
foreach ($inputName in $requiredInputs) {
    $inputPath = Join-Path $SourceRoot $inputName
    if (-not (Test-Path -LiteralPath $inputPath)) {
        throw "Missing Npro pbctopo smoke input: $inputPath"
    }
}

if ([string]::IsNullOrWhiteSpace($WorkRoot)) {
    $WorkRoot = Join-Path $SourceRoot 'test_pbctopo_component_graph_precompute_lazy_smoke'
}
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$argPath = Join-Path $WorkRoot 'arg.in'
$outputName = 'tetramer_pbctopo_component_graph_precompute_lazy_smoke.pdb'
$outputPath = Join-Path $WorkRoot $outputName

$argsText = @"
-ff=charmm36
-s
$SourceRoot\md1_snap0.pdb
-traj
$SourceRoot\md1_1.xtc
-pbctopo
-pbctopo-exact-component-limit=16
-pbctopo-beam-width=128
-pbctopo-anchor-backend=auto
-g=7-11
-nt=8
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
-opdb=$outputName
"@
[System.IO.File]::WriteAllText($argPath, $argsText, [System.Text.Encoding]::ASCII)

Push-Location $WorkRoot
try {
    & $TitanExe -a $argPath
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}

if (-not (Test-Path -LiteralPath $outputPath)) {
    throw "Missing pbctopo output PDB: $outputPath"
}
$sha = (Get-FileHash -LiteralPath $outputPath -Algorithm SHA256).Hash.ToUpperInvariant()
if ($sha -ne $BaselineSha256.ToUpperInvariant()) {
    throw "PDB SHA256 changed: got $sha expected $BaselineSha256"
}

$debugLog = Join-Path $WorkRoot 'pbctopo_debug.log'
if (-not (Test-Path -LiteralPath $debugLog)) {
    throw "Missing pbctopo debug log: $debugLog"
}
$logText = Get-Content -Raw -LiteralPath $debugLog
foreach ($pattern in @(
        'mode=exact, raw_components=4, effective_components=4',
        'pbci_fallback_frames=0',
        'component_graph_shift_descriptors=\d+',
        'component_graph_shift_materialized=\d+',
        'component_graph_bbox_only_scores=\d+',
        'component_graph_lazy_materialize=([0-9.]+) s',
        'component_graph_refine_trial_build_incremental=\d+',
        'component_graph_refine_trial_build_full=\d+',
        'component_graph_refine_trial_build_global_rebase=\d+',
        'component_graph_refine_stats_cached=\d+',
        'component_graph_refine_stats_full=\d+',
        'component_graph_refine_stats=([0-9.]+) s',
        'component_cycle_sync_attempts=\d+',
        'component_cycle_sync_used=\d+',
        'component_cycle_sync_pair_candidates=\d+',
        'component_cycle_sync=([0-9.]+) s',
        'pbctopo_topology_bonds_available=\d+',
        'pbctopo_topology_bond_edges_added=\d+',
        'pbctopo_topology_bond_edges_skipped=\d+',
        'pbctopo_topology_source=(gromacs_topology\+geometry_cutoff|geometry_cutoff)',
        'pbctopo_topology_component_source=(topology\+geometry|geometry)'
    )) {
    if ($logText -notmatch $pattern) {
        throw "Missing expected pbctopo lazy-graph log pattern: $pattern"
    }
}

$incrementalMatch = [regex]::Match(
    $logText, 'component_graph_refine_trial_build_incremental=(\d+)')
if (-not $incrementalMatch.Success -or
    [int]$incrementalMatch.Groups[1].Value -le 0) {
    throw 'component_graph_refine_trial_build_incremental did not record any incremental refine trials'
}

Write-Host 'PASS pbctopo_component_graph_lazy_smoke'
