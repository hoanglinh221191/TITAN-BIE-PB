param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_topology_continuity_ranges_smoke')
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

$pdb = Join-Path $WorkRoot 'topology_continuity_ranges_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$arg = Join-Path $WorkRoot 'arg.in'
$outName = 'topology_continuity_ranges_out.pdb'

@'
REMARK    PBCTOPO topology continuity range smoke
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
pbctopo topology continuity range smoke

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
-opdb=$outName
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
if (-not (Test-Path -LiteralPath $debugLog)) {
    throw "Missing pbctopo debug log: $debugLog"
}
$logText = Get-Content -Raw -LiteralPath $debugLog

foreach ($pattern in @(
        'raw_components=3, effective_components=3',
        'pbctopo_topology_bonds_available=2',
        'pbctopo_topology_bond_edges_added=2',
        'pbctopo_topology_source=gromacs_topology\+geometry_cutoff',
        'pbctopo_topology_component_source=topology\+geometry',
        'component_graph_used=1',
        'component_score_continuity_rejects=0'
    )) {
    if ($logText -notmatch $pattern) {
        throw "Missing expected topology continuity-range pbctopo log pattern: $pattern"
    }
}

if ($logText -match 'fell back to pbci=0') {
    throw 'Topology continuity range fixture unexpectedly fell back to pbci=0'
}

Write-Host 'PASS pbctopo_topology_continuity_ranges_smoke'
