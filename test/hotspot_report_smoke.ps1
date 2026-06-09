param(
  [string]$TitanExe = ".\titan.exe",
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\hotspot_report_smoke")
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$TitanExe = if ([System.IO.Path]::IsPathRooted($TitanExe)) {
  (Resolve-Path -LiteralPath $TitanExe).Path
} else {
  (Resolve-Path -LiteralPath (Join-Path $repo $TitanExe)).Path
}

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$unitExe = Join-Path $WorkRoot "titan_hotspot_report_smoke.exe"
& $Compiler -std=c++23 "-I$(Join-Path $repo 'src')" `
  (Join-Path $repo "test\titan_hotspot_report_smoke.cpp") `
  (Join-Path $repo "src\titan_hotspot_report.cpp") `
  -o $unitExe
if ($LASTEXITCODE -ne 0) {
  throw "Failed to compile hotspot report smoke"
}
& $unitExe
if ($LASTEXITCODE -ne 0) {
  throw "Hotspot report C++ smoke failed"
}

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

function Invoke-Titan([string]$CaseDir, [string]$ArgText, [bool]$ExpectSuccess = $true) {
  $arg = Join-Path $CaseDir "arg.in"
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Push-Location $CaseDir
  try {
    $output = & $TitanExe -a $arg 2>&1
    $code = $LASTEXITCODE
  } finally {
    Pop-Location
  }
  if ($ExpectSuccess -and $code -ne 0) {
    throw "TITAN failed with code $code. Output:`n$output"
  }
  if ((-not $ExpectSuccess) -and $code -eq 0) {
    throw "TITAN was expected to fail. Output:`n$output"
  }
  return ($output -join "`n")
}

if (Test-Path -LiteralPath $WorkRoot) {
  Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$fixture = Join-Path $WorkRoot "hotspot_fixture.pdb"
Write-CasFixture $fixture
$lib = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $repo "lib_ff")).Path.TrimEnd('\') + '\')
$base = @"
-f=$(Convert-ToTitanPath $fixture)
-1pdb=true
-ff=AMBER99SB
-g=1&1
-mol-1-type=protein
-mol-2-type=protein
-path-a=$lib
-path-n=$lib
-nt=2
"@

$caseBoth = Join-Path $WorkRoot "both"
New-Item -ItemType Directory -Force -Path $caseBoth | Out-Null
Copy-Item -LiteralPath $fixture -Destination (Join-Path $caseBoth "hotspot_fixture.pdb")
Invoke-Titan $caseBoth ($base + "`n-hotspot-report=true`n") | Out-Null
foreach ($file in @("hotspot_residue_ranking.csv", "hotspot_summary.dat", "hotspot_bfactor.pdb", "cas_receptor.dat", "edecomp_mm_residue_ranking.csv")) {
  if (-not (Test-Path -LiteralPath (Join-Path $caseBoth $file))) {
    throw "hotspot both: missing $file"
  }
}
if (-not (Select-String -Path (Join-Path $caseBoth "hotspot_summary.dat") -Pattern "CAS method: GBSA alanine-like side-chain neutralization" -Quiet)) {
  throw "hotspot both: summary missing CAS method note"
}
if (-not (Select-String -Path (Join-Path $caseBoth "hotspot_summary.dat") -Pattern "EDecomp note: endpoint entity decomposition" -Quiet)) {
  throw "hotspot both: summary missing EDecomp note"
}

$caseEdecomp = Join-Path $WorkRoot "edecomp_only"
New-Item -ItemType Directory -Force -Path $caseEdecomp | Out-Null
Copy-Item -LiteralPath $fixture -Destination (Join-Path $caseEdecomp "hotspot_fixture.pdb")
Invoke-Titan $caseEdecomp ($base + "`n-hotspot-report=true`n-hotspot-source=edecomp`n") | Out-Null
if (-not (Test-Path -LiteralPath (Join-Path $caseEdecomp "hotspot_residue_ranking.csv"))) {
  throw "hotspot edecomp-only: ranking missing"
}
if (Test-Path -LiteralPath (Join-Path $caseEdecomp "cas_receptor.dat")) {
  throw "hotspot edecomp-only: CAS output should not be generated"
}
if (-not (Select-String -Path (Join-Path $caseEdecomp "hotspot_summary.dat") -Pattern "bfactor_status:" -Quiet)) {
  throw "hotspot edecomp-only: summary missing B-factor status"
}

$caseOff = Join-Path $WorkRoot "off"
New-Item -ItemType Directory -Force -Path $caseOff | Out-Null
Copy-Item -LiteralPath $fixture -Destination (Join-Path $caseOff "hotspot_fixture.pdb")
Invoke-Titan $caseOff $base | Out-Null
if (Test-Path -LiteralPath (Join-Path $caseOff "hotspot_residue_ranking.csv")) {
  throw "hotspot off: unexpected ranking output"
}

Write-Host "PASS hotspot_report_smoke"
