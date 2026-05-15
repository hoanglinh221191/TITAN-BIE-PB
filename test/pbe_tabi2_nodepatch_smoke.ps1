param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi2_nodepatch_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }
function New-CleanDir([string]$Path) {
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}
function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Encoding ASCII -Value @(
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    )
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$runDir = Join-Path $WorkRoot 'mini_protein_tabi2_nodepatch'
New-CleanDir $runDir
$pdbPath = Join-Path $runDir 'driver.pdb'
Write-DummyPdb $pdbPath
$pqrPath = Join-Path $InputRoot 'mini-protein.pqr'
$libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
$argPath = Join-Path $runDir 'arg.in'
Set-Content -LiteralPath $argPath -Encoding ASCII -Value @(
    "-f=$(Convert-ToTitanPath $pdbPath)",
    '-1pdb',
    '-ff=AMBER99SB',
    '-g=1&1',
    "-path-a=$libPath",
    "-path-n=$libPath",
    '-pbe=tabi',
    "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
    '-pbe-surface=vdw',
    '-pbe-timing=detail',
    '-pbe-max-iter=800',
    '-pbe-rms=1e-4',
    '-pbe-tabi-sum=direct',
    '-pbe-tabi-formulation=tabi2-nodepatch',
    '-pbe-tabi-precond=block',
    '-nt=1'
)
Push-Location $runDir
try {
    & $TitanExe -a $argPath | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE" }
} finally {
    Pop-Location
}

$energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
$audit = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_formulation_audit.csv'))[0]
$profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
if ($energy.status -ne 'ok') { throw "TABI2 node-patch smoke failed status=$($energy.status)" }
if ($audit.formulation -ne 'tabi2-nodepatch' -or $audit.discretization -ne 'node_patch') {
    throw "TABI2 audit did not report node-patch mode."
}
if ($audit.preconditioner_effective -ne 'diagonal_block_pending') {
    throw "TABI2 block preconditioner audit did not report pending block path: $($audit.preconditioner_effective)"
}
if ([int64]$profile.unknown_count -le 0) { throw 'TABI2 profile missing unknown count.' }
Write-Host ("PBE TABI2 node-patch smoke PASS: panels={0} unknowns={1} E_kJ={2}" -f
    $audit.panels, $profile.unknown_count, $energy.E_rxn_kJmol)
