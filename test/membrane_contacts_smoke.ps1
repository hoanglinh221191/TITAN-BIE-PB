param(
  [string]$TitanExe = ".\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\membrane_contacts_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null
$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$outputs = @(
  "membrane_protein_lipid_contacts.csv",
  "membrane_lipid_enrichment.csv",
  "membrane_curvature_map.csv",
  "membrane_proximity_coupling.csv",
  "membrane_coupling_summary.dat",
  "membrane_lipid_enrichment_average.csv",
  "membrane_proximity_coupling_average.csv",
  "membrane_coupling_summary_average.dat",
  "membrane_curvature_map_average.csv"
)

function Clear-Outputs {
  foreach ($file in $outputs) {
    $path = Join-Path $repo $file
    if (Test-Path -LiteralPath $path) {
      Remove-Item -LiteralPath $path -Force
    }
  }
  Get-ChildItem -LiteralPath $repo -Filter "membrane_protein_lipid_contacts_[0-9]*.csv" -ErrorAction SilentlyContinue |
    Remove-Item -Force
  Get-ChildItem -LiteralPath $repo -Filter "membrane_lipid_enrichment_[0-9]*.csv" -ErrorAction SilentlyContinue |
    Remove-Item -Force
  Get-ChildItem -LiteralPath $repo -Filter "membrane_coupling_summary_[0-9]*.dat" -ErrorAction SilentlyContinue |
    Remove-Item -Force
  Get-ChildItem -LiteralPath $repo -Filter "membrane_curvature_map_[0-9]*.csv" -ErrorAction SilentlyContinue |
    Remove-Item -Force
  Get-ChildItem -LiteralPath $repo -Filter "membrane_proximity_coupling_[0-9]*.csv" -ErrorAction SilentlyContinue |
    Remove-Item -Force
}

function Write-Fixture([string]$Path, [string]$NdxPath) {
  $pdb = @"
CRYST1   30.000   30.000   40.000  90.00  90.00  90.00 P 1           1
MODEL        1
TITLE     t= 0.00000 step= 0
ATOM      1 N    LEU A   1       5.000   5.000   6.000  1.00  0.00           N
ATOM      2 CA   LEU A   1       5.500   5.000   6.000  1.00  0.00           C
ATOM      3 CB   LEU A   1       5.500   5.500   5.800  1.00  0.00           C
ATOM      4 CG   LEU A   1       5.400   5.400   5.400  1.00  0.00           C
ATOM      5 N    LYS B   1      22.000  22.000   6.000  1.00  0.00           N
ATOM      6 CA   LYS B   1      22.500  22.000   6.000  1.00  0.00           C
ATOM      7 CB   LYS B   1      22.500  22.500   5.800  1.00  0.00           C
ATOM      8 CG   LYS B   1      22.400  22.400   5.400  1.00  0.00           C
ATOM      9 PO4  DOPCL   2       6.000   5.500   5.000  1.00  0.00           P
ATOM     10 C1A  DOPCL   2       6.000   5.500   4.000  1.00  0.00           C
ATOM     11 C2A  DOPCL   2       6.000   5.500   3.000  1.00  0.00           C
ATOM     12 PO4  DOPSL   3      22.000  22.000   5.000  1.00  0.00           P
ATOM     13 C1A  DOPSL   3      22.000  22.000   4.000  1.00  0.00           C
ATOM     14 C2A  DOPSL   3      22.000  22.000   3.000  1.00  0.00           C
ATOM     15 PO4  DOPCL   4      15.000  15.000  -5.000  1.00  0.00           P
ATOM     16 C1A  DOPCL   4      15.000  15.000  -4.000  1.00  0.00           C
ATOM     17 C2A  DOPCL   4      15.000  15.000  -3.000  1.00  0.00           C
ATOM     18 PO4  DOPSL   5      26.000  26.000  -5.000  1.00  0.00           P
ATOM     19 C1A  DOPSL   5      26.000  26.000  -4.000  1.00  0.00           C
ATOM     20 C2A  DOPSL   5      26.000  26.000  -3.000  1.00  0.00           C
ENDMDL
MODEL        2
TITLE     t= 1.00000 step= 1
ATOM      1 N    LEU A   1       5.000   5.000   6.000  1.00  0.00           N
ATOM      2 CA   LEU A   1       5.500   5.000   6.000  1.00  0.00           C
ATOM      3 CB   LEU A   1       5.500   5.500   5.800  1.00  0.00           C
ATOM      4 CG   LEU A   1       5.400   5.400   5.400  1.00  0.00           C
ATOM      5 N    LYS B   1      22.000  22.000   6.000  1.00  0.00           N
ATOM      6 CA   LYS B   1      22.500  22.000   6.000  1.00  0.00           C
ATOM      7 CB   LYS B   1      22.500  22.500   5.800  1.00  0.00           C
ATOM      8 CG   LYS B   1      22.400  22.400   5.400  1.00  0.00           C
ATOM      9 PO4  DOPCL   2       6.200   5.500   5.000  1.00  0.00           P
ATOM     10 C1A  DOPCL   2       6.200   5.500   4.000  1.00  0.00           C
ATOM     11 C2A  DOPCL   2       6.200   5.500   3.000  1.00  0.00           C
ATOM     12 PO4  DOPSL   3      22.000  22.000   5.000  1.00  0.00           P
ATOM     13 C1A  DOPSL   3      22.000  22.000   4.000  1.00  0.00           C
ATOM     14 C2A  DOPSL   3      22.000  22.000   3.000  1.00  0.00           C
ATOM     15 PO4  DOPCL   4      15.000  15.000  -5.000  1.00  0.00           P
ATOM     16 C1A  DOPCL   4      15.000  15.000  -4.000  1.00  0.00           C
ATOM     17 C2A  DOPCL   4      15.000  15.000  -3.000  1.00  0.00           C
ATOM     18 PO4  DOPSL   5      26.000  26.000  -5.000  1.00  0.00           P
ATOM     19 C1A  DOPSL   5      26.000  26.000  -4.000  1.00  0.00           C
ATOM     20 C2A  DOPSL   5      26.000  26.000  -3.000  1.00  0.00           C
ENDMDL
"@
  [System.IO.File]::WriteAllText($Path, $pdb, [System.Text.Encoding]::ASCII)
  [System.IO.File]::WriteAllText(
    $NdxPath,
    "[ protein ]`r`n1 2 3 4 5 6 7 8`r`n`r`n[ lipid ]`r`n9 10 11 12 13 14 15 16 17 18 19 20`r`n",
    [System.Text.Encoding]::ASCII)
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-Outputs
  Push-Location $repo
  try {
    & $titanPath -a $arg
    if ($LASTEXITCODE -ne 0) {
      throw "${Name}: TITAN exited with code $LASTEXITCODE"
    }
  } finally {
    Pop-Location
  }
}

function Assert-Exists([string]$File) {
  if (-not (Test-Path -LiteralPath (Join-Path $repo $File))) {
    throw "missing expected output $File"
  }
}

function Assert-NotExists([string]$File) {
  if (Test-Path -LiteralPath (Join-Path $repo $File)) {
    throw "unexpected output $File"
  }
}

function Assert-Contains([string]$File, [string]$Pattern) {
  if (-not (Select-String -Path (Join-Path $repo $File) -Pattern $Pattern -Quiet)) {
    throw "$File does not contain pattern $Pattern"
  }
}

$fixture = Join-Path $work "protein_lipid_synthetic.pdb"
$index = Join-Path $work "protein_lipid_synthetic.ndx"
Write-Fixture $fixture $index

$base = @"
-f=$fixture
-n=$index
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
-membrane-group=auto
-membrane-normal=z
-membrane-grid=2,2
"@

Run-TitanCase "contacts" ($base + "`r`n" + @"
-membrane-protein-contacts=true
"@)
Assert-Exists "membrane_protein_lipid_contacts.csv"
Assert-Exists "membrane_proximity_coupling.csv"
Assert-Exists "membrane_coupling_summary.dat"
Assert-Contains "membrane_protein_lipid_contacts.csv" "mol1\|A:1:LEU"
Assert-Contains "membrane_protein_lipid_contacts.csv" "mol1\|B:1:LYS"
Assert-Contains "membrane_proximity_coupling.csv" "apl_A2"
Assert-Contains "membrane_proximity_coupling.csv" "contact_probability"
Assert-Contains "membrane_coupling_summary.dat" "protein_lipid_contact_cut_A: 6.5"

Run-TitanCase "enrichment" ($base + "`r`n" + @"
-membrane-enrichment=true
"@)
Assert-Exists "membrane_lipid_enrichment.csv"
Assert-Exists "membrane_protein_lipid_contacts.csv"
Assert-Exists "membrane_proximity_coupling.csv"
Assert-Contains "membrane_lipid_enrichment.csv" "DOPC"
Assert-Contains "membrane_lipid_enrichment.csv" "DOPS"
Assert-Contains "membrane_proximity_coupling.csv" "lipid_enrichment"

Run-TitanCase "curvature" ($base + "`r`n" + @"
-membrane-enrichment=true
-membrane-curvature=grid
"@)
Assert-Exists "membrane_curvature_map.csv"
Assert-Exists "membrane_proximity_coupling.csv"
Assert-Contains "membrane_curvature_map.csv" "curvature_proxy_Ainv"
Assert-Contains "membrane_proximity_coupling.csv" "thickness_A"
Assert-Contains "membrane_proximity_coupling.csv" "order_S"
Assert-Contains "membrane_proximity_coupling.csv" "curvature_proxy_Ainv"

Run-TitanCase "off" $base
Assert-NotExists "membrane_protein_lipid_contacts.csv"
Assert-NotExists "membrane_lipid_enrichment.csv"
Assert-NotExists "membrane_curvature_map.csv"
Assert-NotExists "membrane_proximity_coupling.csv"
Assert-NotExists "membrane_coupling_summary.dat"

$listFile = Join-Path $work "same_pdb_twice.txt"
[System.IO.File]::WriteAllText($listFile, $fixture + "`r`n" + $fixture + "`r`n", [System.Text.Encoding]::ASCII)
Run-TitanCase "multi" (@"
-f=$fixture
-i=$listFile
-n=$index
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
-membrane-protein-contacts=true
-membrane-enrichment=true
-membrane-curvature=grid
-membrane-group=auto
-membrane-normal=z
-membrane-grid=2,2
"@)
foreach ($file in @(
    "membrane_protein_lipid_contacts_0.csv",
    "membrane_protein_lipid_contacts_1.csv",
    "membrane_lipid_enrichment_0.csv",
    "membrane_lipid_enrichment_1.csv",
    "membrane_coupling_summary_0.dat",
    "membrane_coupling_summary_1.dat",
    "membrane_lipid_enrichment_average.csv",
    "membrane_coupling_summary_average.dat",
    "membrane_curvature_map_average.csv")) {
  Assert-Exists $file
}
foreach ($file in @(
    "membrane_proximity_coupling_0.csv",
    "membrane_proximity_coupling_1.csv",
    "membrane_proximity_coupling_average.csv")) {
  Assert-Exists $file
}
Assert-Contains "membrane_proximity_coupling_average.csv" "delta_near_minus_far"
Assert-Contains "membrane_coupling_summary_average.dat" "replicates: 2"

Write-Host "PASS membrane_contacts_smoke"
