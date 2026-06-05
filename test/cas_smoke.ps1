param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\cas_smoke'),
    [string]$LibRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) 'lib_ff')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$env:TITAN_NO_PAUSE = '1'

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Write-CasFixture([string]$Path) {
    @(
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'MODEL        1',
        'ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  LYS A   1       1.450   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   LYS A   1       2.050   1.350   0.000  1.00  0.00           C',
        'ATOM      4  O   LYS A   1       1.450   2.350   0.000  1.00  0.00           O',
        'ATOM      5  CB  LYS A   1       1.950  -0.750   1.250  1.00  0.00           C',
        'ATOM      6  CG  LYS A   1       3.350  -0.750   1.250  1.00  0.00           C',
        'ATOM      7  CD  LYS A   1       3.900  -1.950   1.900  1.00  0.00           C',
        'ATOM      8  CE  LYS A   1       5.200  -1.950   1.900  1.00  0.00           C',
        'ATOM      9  NZ  LYS A   1       5.800  -3.100   2.450  1.00  0.00           N',
        'ATOM     10  N   ALA B   2       4.500   0.000   0.000  1.00  0.00           N',
        'ATOM     11  CA  ALA B   2       5.950   0.000   0.000  1.00  0.00           C',
        'ATOM     12  C   ALA B   2       6.550   1.350   0.000  1.00  0.00           C',
        'ATOM     13  O   ALA B   2       5.950   2.350   0.000  1.00  0.00           O',
        'ATOM     14  CB  ALA B   2       6.450  -0.750   1.250  1.00  0.00           C',
        'ENDMDL',
        'MODEL        2',
        'ATOM      1  N   LYS A   1       0.050   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  LYS A   1       1.500   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   LYS A   1       2.100   1.350   0.000  1.00  0.00           C',
        'ATOM      4  O   LYS A   1       1.500   2.350   0.000  1.00  0.00           O',
        'ATOM      5  CB  LYS A   1       2.000  -0.750   1.250  1.00  0.00           C',
        'ATOM      6  CG  LYS A   1       3.400  -0.750   1.250  1.00  0.00           C',
        'ATOM      7  CD  LYS A   1       3.950  -1.950   1.900  1.00  0.00           C',
        'ATOM      8  CE  LYS A   1       5.250  -1.950   1.900  1.00  0.00           C',
        'ATOM      9  NZ  LYS A   1       5.850  -3.100   2.450  1.00  0.00           N',
        'ATOM     10  N   ALA B   2       4.550   0.000   0.000  1.00  0.00           N',
        'ATOM     11  CA  ALA B   2       6.000   0.000   0.000  1.00  0.00           C',
        'ATOM     12  C   ALA B   2       6.600   1.350   0.000  1.00  0.00           C',
        'ATOM     13  O   ALA B   2       6.000   2.350   0.000  1.00  0.00           O',
        'ATOM     14  CB  ALA B   2       6.500  -0.750   1.250  1.00  0.00           C',
        'ENDMDL',
        'END'
    ) | Set-Content -LiteralPath $Path -Encoding ASCII
}

function Write-ArgFile([string]$Path, [string]$Pdb, [string]$Mol2Type, [string[]]$ExtraLines) {
    $lib = Convert-ToTitanPath ((Resolve-Path -LiteralPath $LibRoot).Path.TrimEnd('\') + '\')
    @(
        "-f=$(Convert-ToTitanPath $Pdb)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        '-mol-1-type=protein',
        "-mol-2-type=$Mol2Type",
        "-path-a=$lib",
        "-path-n=$lib",
        '-nt=2'
    ) + $ExtraLines | Set-Content -LiteralPath $Path -Encoding ASCII
}

function Invoke-Titan([string]$CaseDir, [string]$ArgPath, [bool]$ExpectSuccess) {
    Push-Location $CaseDir
    $oldErrorActionPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = 'Continue'
        $output = & $TitanExe -a $ArgPath 2>&1
        $code = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $oldErrorActionPreference
        Pop-Location
    }
    if ($ExpectSuccess -and $code -ne 0) {
        throw "TITAN exited with code $code. Output:`n$output"
    }
    if ((-not $ExpectSuccess) -and $code -eq 0) {
        throw "TITAN was expected to fail but exited successfully. Output:`n$output"
    }
    return ($output -join "`n")
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $LibRoot)) {
    throw "Missing TITAN force-field library directory: $LibRoot"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$caseOn = Join-Path $WorkRoot 'cas_on'
New-Item -ItemType Directory -Path $caseOn -Force | Out-Null
$pdbOn = Join-Path $caseOn 'cas_fixture.pdb'
Write-CasFixture $pdbOn
$argOn = Join-Path $caseOn 'arg.in'
Write-ArgFile $argOn $pdbOn 'protein' @('-cas=true')
$onOutput = Invoke-Titan $caseOn $argOn $true

$required = @(
    'cas_receptor.dat',
    'cas_ligand.dat',
    'cas_hotspots_b_factor.pdb',
    'cas_summary.dat',
    'cas_performance.dat',
    'gbsa_summary.dat'
)
foreach ($file in $required) {
    if (-not (Test-Path -LiteralPath (Join-Path $caseOn $file))) {
        throw "cas=true: missing expected output $file"
    }
}
if ($onOutput -notmatch 'CAS auto-enabled GBSA') {
    throw 'cas=true: parser did not log GBSA auto-enable message'
}
if (-not (Select-String -Path (Join-Path $caseOn 'cas_receptor.dat') -Pattern 'role,chainnum,chain,resnr,resname,mutation_status' -Quiet)) {
    throw 'cas=true: receptor output lacks stable residue identity columns'
}
if (-not (Select-String -Path (Join-Path $caseOn 'cas_receptor.dat') -Pattern 'valid_mutation' -Quiet)) {
    throw 'cas=true: receptor output lacks valid mutation rows'
}
if (-not (Select-String -Path (Join-Path $caseOn 'cas_receptor.dat') -Pattern 'noop_alanine' -Quiet)) {
    throw 'cas=true: receptor output lacks ALA noop policy rows'
}
if (-not (Select-String -Path (Join-Path $caseOn 'cas_summary.dat') -Pattern 'GBSA alanine-like side-chain neutralization scan' -Quiet)) {
    throw 'cas=true: summary does not document CAS method'
}
if (-not (Select-String -Path (Join-Path $caseOn 'cas_performance.dat') -Pattern 'mutation_evaluations' -Quiet)) {
    throw 'cas=true: performance output lacks mutation counter'
}

$caseOff = Join-Path $WorkRoot 'cas_off'
New-Item -ItemType Directory -Path $caseOff -Force | Out-Null
$pdbOff = Join-Path $caseOff 'cas_fixture.pdb'
Copy-Item -LiteralPath $pdbOn -Destination $pdbOff
$argOff = Join-Path $caseOff 'arg.in'
Write-ArgFile $argOff $pdbOff 'protein' @('-cas=false')
Invoke-Titan $caseOff $argOff $true | Out-Null
foreach ($file in @('cas_receptor.dat', 'cas_summary.dat', 'cas_performance.dat')) {
    if (Test-Path -LiteralPath (Join-Path $caseOff $file)) {
        throw "cas=false: unexpected CAS output $file"
    }
}

$caseNocal = Join-Path $WorkRoot 'cas_nocal'
New-Item -ItemType Directory -Path $caseNocal -Force | Out-Null
$pdbNocal = Join-Path $caseNocal 'cas_fixture.pdb'
Copy-Item -LiteralPath $pdbOn -Destination $pdbNocal
$argNocal = Join-Path $caseNocal 'arg.in'
Write-ArgFile $argNocal $pdbNocal 'protein' @('-nocal=true', '-cas=true')
$failOutput = Invoke-Titan $caseNocal $argNocal $false
if ($failOutput -notmatch 'CAS requires GBSA energy evaluation') {
    throw 'nocal+cas: rejection message does not explain CAS/GBSA requirement'
}

$caseSmallLig = Join-Path $WorkRoot 'cas_small_ligand'
New-Item -ItemType Directory -Path $caseSmallLig -Force | Out-Null
$pdbSmallLig = Join-Path $caseSmallLig 'cas_fixture.pdb'
Copy-Item -LiteralPath $pdbOn -Destination $pdbSmallLig
$argSmallLig = Join-Path $caseSmallLig 'arg.in'
Write-ArgFile $argSmallLig $pdbSmallLig 'lig' @('-cas=true')
Invoke-Titan $caseSmallLig $argSmallLig $true | Out-Null
if (-not (Test-Path -LiteralPath (Join-Path $caseSmallLig 'cas_receptor.dat'))) {
    throw 'small-ligand: receptor CAS output is missing'
}
if (Test-Path -LiteralPath (Join-Path $caseSmallLig 'cas_ligand.dat')) {
    throw 'small-ligand: ligand CAS output should be skipped for non-protein mol2'
}
if (-not (Select-String -Path (Join-Path $caseSmallLig 'cas_summary.dat') -Pattern 'ligand_note: mol2 is not protein-like' -Quiet)) {
    throw 'small-ligand: summary does not document skipped ligand scan'
}

Write-Host 'CAS smoke PASS'
