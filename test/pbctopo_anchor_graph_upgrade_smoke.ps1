param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_anchor_graph_upgrade_smoke')
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

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$pdb = Join-Path $WorkRoot 'anchor_graph_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$arg = Join-Path $WorkRoot 'arg.in'

@'
REMARK    PBCTOPO anchor graph upgrade smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  CA  ALA A   1       8.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  GLY A   2      10.000  10.000  10.000  1.00  0.00           C
ATOM      3  CA  SER A   3      12.000  10.000  10.000  1.00  0.00           C
ATOM      4  CA  ALA B   1       8.000  13.000  10.000  1.00  0.00           C
ATOM      5  CA  GLY B   2      10.000  13.000  10.000  1.00  0.00           C
ATOM      6  CA  SER B   3      12.000  13.000  10.000  1.00  0.00           C
ENDMDL
TITLE     Title t= 1.00000 step= 1
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
MODEL        2
ATOM      1  CA  ALA A   1       8.000  10.000  10.000  1.00  0.00           C
ATOM      2  CA  GLY A   2      10.000  10.000  10.000  1.00  0.00           C
ATOM      3  CA  SER A   3      12.000  10.000  10.000  1.00  0.00           C
ATOM      4  CA  ALA B   1       8.000  13.000  10.000  1.00  0.00           C
ATOM      5  CA  GLY B   2      10.000  13.000  10.000  1.00  0.00           C
ATOM      6  CA  SER B   3      12.000  13.000  10.000  1.00  0.00           C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ OWNER_A ]
1 2 3

[ OWNER_B ]
4 5 6
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
[ moleculetype ]
OWNER_A 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0
3 CT 3 SER CA 1 0.0 12.0

[ bonds ]
1 2 1
2 3 1

[ moleculetype ]
OWNER_B 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0
3 CT 3 SER CA 1 0.0 12.0

[ bonds ]
1 2 1
2 3 1

[ system ]
pbctopo anchor graph upgrade smoke

[ molecules ]
OWNER_A 1
OWNER_B 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

$argsText = @"
-ff=charmm36
-f
$pdb
-pbctopo
-pbctopo-bond-cutoff=2.50
-pbctopo-contact-cutoff=4.50
-pbctopo-exact-component-limit=6
-pbctopo-beam-width=64
-pbctopo-anchor-backend=auto
-g=0-1
-nt=1
-1pdb
-nocal
-n
$ndx
-common_topo
$top
-no-dssp
-noprint-side_chain_contact
-opdb=anchor_graph_upgrade_out.pdb
"@
[System.IO.File]::WriteAllText($arg, $argsText, [System.Text.Encoding]::ASCII)

Push-Location $WorkRoot
try {
    & $titanPath -a $arg
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
        'pbctopo_anchor_context_build=([0-9.]+) s',
        'pbctopo_anchor_order_edges=\d+',
        'pbctopo_anchor_endpoint_edges=\d+',
        'pbctopo_anchor_contact_pairs=([1-9]\d*)',
        'pbctopo_anchor_residue_order_violations=\d+',
        'pbctopo_anchor_ambiguous_frames=\d+',
        'component_graph_refine_anchor_order=([0-9.]+) s',
        'component_graph_refine_anchor_contact=([0-9.]+) s'
    )) {
    if ($logText -notmatch $pattern) {
        throw "Missing expected anchor graph log pattern: $pattern"
    }
}

$rows = @(Import-Csv -LiteralPath $report)
if ($rows.Count -ne 1) {
    throw "Expected exactly one report row, found $($rows.Count)"
}
$requiredColumns = @(
    'anchor_order_max_d2',
    'anchor_order_path_d2',
    'anchor_endpoint_image_gap',
    'anchor_residue_order_violations',
    'anchor_contact_map_lost',
    'anchor_contact_map_max_d2',
    'anchor_ambiguity_count'
)
foreach ($column in $requiredColumns) {
    if (-not ($rows[0].PSObject.Properties.Name -contains $column)) {
        throw "Missing frame report column: $column"
    }
}

if ($rows[0].status -ne 'certified') {
    throw "Expected certified frame, got $($rows[0].status)"
}
if ([int]$rows[0].anchor_contact_map_lost -ne 0) {
    throw "Expected intact owner contact map, got lost=$($rows[0].anchor_contact_map_lost)"
}
if ([double]$rows[0].anchor_contact_map_max_d2 -le 0.0) {
    throw "Expected positive contact-map distance metric"
}
if ($rows[0].best_score_vector -notmatch 'anchor_contact_lost=0') {
    throw "best_score_vector does not expose anchor contact-map score"
}

Write-Host 'PASS pbctopo_anchor_graph_upgrade_smoke'
