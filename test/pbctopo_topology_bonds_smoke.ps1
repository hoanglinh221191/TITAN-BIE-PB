param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_topology_bonds_smoke')
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

$pdb = Join-Path $WorkRoot 'topology_bonds_fixture.pdb'
$ndx = Join-Path $WorkRoot 'index.ndx'
$top = Join-Path $WorkRoot 'topol.top'
$arg = Join-Path $WorkRoot 'arg.in'
$outName = 'topology_bonds_out.pdb'

@'
REMARK    PBCTOPO topology-bond smoke
TITLE     Title t= 0.00000 step= 0
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1 A1   RES A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2 A2   RES A   1       3.000   0.000   0.000  1.00  0.00           C
ATOM      3 A3   RES A   1       6.000   0.000   0.000  1.00  0.00           C
ATOM      4 B1   OTH B   2      12.000   0.000   0.000  1.00  0.00           C
ENDMDL
TITLE     Title t= 1.00000 step= 1
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
MODEL        2
ATOM      1 A1   RES A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2 A2   RES A   1       3.000   0.000   0.000  1.00  0.00           C
ATOM      3 A3   RES A   1       6.000   0.000   0.000  1.00  0.00           C
ATOM      4 B1   OTH B   2      12.000   0.000   0.000  1.00  0.00           C
ENDMDL
END
'@ | Set-Content -LiteralPath $pdb -Encoding ASCII

@'
[ TEST ]
1 2 3

[ OTHER ]
4
'@ | Set-Content -LiteralPath $ndx -Encoding ASCII

@'
[ moleculetype ]
TEST 3

[ atoms ]
1 CT 1 RES A1 1 0.0 12.0
2 CT 1 RES A2 1 0.0 12.0
3 CT 1 RES A3 1 0.0 12.0

[ bonds ]
1 2 1
2 3 1

[ moleculetype ]
OTHER 3

[ atoms ]
1 CT 1 OTH B1 1 0.0 12.0

[ system ]
pbctopo topology bonds smoke

[ molecules ]
TEST 1
OTHER 1
'@ | Set-Content -LiteralPath $top -Encoding ASCII

$argsText = @"
-ff=charmm36
-f
$pdb
-pbctopo
-pbctopo-bond-cutoff=2.15
-pbctopo-exact-component-limit=16
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
        'raw_components=2, effective_components=2',
        'pbctopo_topology_bonds_available=2',
        'pbctopo_topology_bond_edges_added=2',
        'pbctopo_topology_source=gromacs_topology\+geometry_cutoff',
        'pbctopo_topology_component_source=topology\+geometry'
    )) {
    if ($logText -notmatch $pattern) {
        throw "Missing expected topology-bond pbctopo log pattern: $pattern"
    }
}

Write-Host 'PASS pbctopo_topology_bonds_smoke'
