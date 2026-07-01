param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_confidence_report_smoke')
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

$pdb = Join-Path $WorkRoot 'confidence_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$arg = Join-Path $WorkRoot 'arg.in'

@'
REMARK    PBCTOPO confidence report smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  CA  ALA     1      10.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      2  CA  GLY     2      12.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      3  CA  ALA     1      50.000  50.000  50.000  1.00  0.00      MOLB C
ATOM      4  CA  GLY     2      52.000  50.000  50.000  1.00  0.00      MOLB C
ATOM      5  CA  DUM C   1      25.000  25.000  25.000  1.00  0.00      DUM  C
ENDMDL
TITLE     Title t= 1.00000 step= 1
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1
MODEL        2
ATOM      1  CA  ALA     1      10.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      2  CA  GLY     2      12.000  10.000  10.000  1.00  0.00      MOLA C
ATOM      3  CA  ALA     1      50.000  50.000  50.000  1.00  0.00      MOLB C
ATOM      4  CA  GLY     2      52.000  50.000  50.000  1.00  0.00      MOLB C
ATOM      5  CA  DUM C   1      25.000  25.000  25.000  1.00  0.00      DUM  C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ PROT ]
1 2 3 4

[ DUMMY ]
5
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
[ moleculetype ]
MOLA 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0

[ bonds ]
1 2 1

[ moleculetype ]
MOLB 3

[ atoms ]
1 CT 1 ALA CA 1 0.0 12.0
2 CT 2 GLY CA 1 0.0 12.0

[ bonds ]
1 2 1

[ moleculetype ]
DUM 3

[ atoms ]
1 CT 1 DUM CA 1 0.0 12.0

[ system ]
pbctopo confidence report smoke

[ molecules ]
MOLA 1
MOLB 1
DUM 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

$argsText = @"
-ff=charmm36
-f
$pdb
-pbctopo
-pbctopo-bond-cutoff=1.10
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
-opdb=confidence_report_out.pdb
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

$report = Join-Path $WorkRoot 'pbctopo_frame_report.csv'
if (-not (Test-Path -LiteralPath $report)) {
    throw "Missing frame report: $report"
}

$rows = @(Import-Csv -LiteralPath $report)
if ($rows.Count -ne 1) {
    throw "Expected exactly 1 frame report row, found $($rows.Count)"
}

$requiredColumns = @(
    'frame',
    'time_ps',
    'status',
    'mode',
    'component_count',
    'safe_candidate_count',
    'best_score_vector',
    'second_score_vector',
    'score_gap',
    'clash_margin',
    'continuity_margin',
    'anchor_margin',
    'split_owner_margin',
    'fallback_reason',
    'time_total_ms',
    'time_graph_ms',
    'time_virtual_ms',
    'time_component_ms',
    'cycle_support',
    'cycle_mismatch',
    'inconsistent_edge_count'
)

$header = (Get-Content -LiteralPath $report -First 1).Split(',')
foreach ($column in $requiredColumns) {
    if ($header -notcontains $column) {
        throw "Missing frame report column: $column"
    }
}

if ($rows[0].status -notmatch '^(certified|rescued|fallback)$') {
    throw "Unexpected status: $($rows[0].status)"
}
if ($rows[0].fallback_reason -notmatch '^(none|no_seed|projection_invalid|full_verify_clash|continuity_reject|anchor_reject|split_owner_reject|component_graph_failed|beam_empty|not_better_than_incumbent|fallback_pbci)$') {
    throw "Unexpected fallback reason: $($rows[0].fallback_reason)"
}
if ([double]$rows[0].time_total_ms -lt 0.0) {
    throw "Negative total timing"
}

Write-Host 'PASS pbctopo_confidence_report_smoke'
