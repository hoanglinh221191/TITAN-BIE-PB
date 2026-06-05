param(
  [string]$TitanExe = ".\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\edecomp_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null

$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$outputs = @(
  "edecomp_mm_residue_matrix_elec.csv",
  "edecomp_mm_residue_matrix_vdw.csv",
  "edecomp_mm_residue_matrix_total.csv",
  "edecomp_mm_residue_edges.csv",
  "edecomp_mm_residue_ranking.csv",
  "edecomp_mm_atom_edges.csv",
  "edecomp_mm_atom_ranking.csv",
  "edecomp_summary.dat"
)

$gbsaOutputs = @(
  "edecomp_mmgbsa_residue_map.csv",
  "edecomp_mmgbsa_residue_ranking.csv",
  "edecomp_mmgbsa_atom_map.csv",
  "edecomp_mmgbsa_atom_ranking.csv"
)

$pbsaOutputs = @(
  "edecomp_mmpbsa_residue_map.csv",
  "edecomp_mmpbsa_residue_ranking.csv",
  "edecomp_mmpbsa_atom_map.csv",
  "edecomp_mmpbsa_atom_ranking.csv"
)

$legacyPatterns = @(
  "average_emap*",
  "gbsa_average_emap*",
  "gbsa_total_energy_maptype*",
  "out_*_stat_energy_res.dat",
  "out_*_energy_stat_res.tsv"
)

function Clear-EDecompOutputs {
  foreach ($file in ($outputs + $gbsaOutputs + $pbsaOutputs)) {
    $path = Join-Path $repo $file
    if (Test-Path -LiteralPath $path) {
      Remove-Item -LiteralPath $path -Force
    }
  }
  Get-ChildItem -LiteralPath $repo -Filter "edecomp_*_*.csv" -ErrorAction SilentlyContinue |
    Remove-Item -Force
  Get-ChildItem -LiteralPath $repo -Filter "edecomp_summary*.dat" -ErrorAction SilentlyContinue |
    Remove-Item -Force
  foreach ($pattern in $legacyPatterns) {
    Get-ChildItem -LiteralPath $repo -Filter $pattern -ErrorAction SilentlyContinue |
      Remove-Item -Force
  }
}

function Assert-NoLegacyEnergyOutputs([string]$Name) {
  foreach ($pattern in $legacyPatterns) {
    $matches = @(Get-ChildItem -LiteralPath $repo -Filter $pattern -ErrorAction SilentlyContinue)
    if ($matches.Count -gt 0) {
      throw "${Name}: legacy energy-map output was generated: $($matches[0].Name)"
    }
  }
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-EDecompOutputs
  Push-Location $repo
  try {
    & $titanPath -a $arg
    if ($LASTEXITCODE -ne 0) {
      throw "${Name}: TITAN exited with code $LASTEXITCODE"
    }
  }
  finally {
    Pop-Location
  }
}

$fixture = Join-Path $repo "example\2\md1_short.pdb"
$ffRoot = Join-Path $repo "lib_ff"

$base = @"
-f=$fixture
-1pdb=true
-path-a=$ffRoot
-path-n=$ffRoot
-g=1-2
-edecomp=true
-edecomp-source=mm
-edecomp-level=both
-edecomp-edge-top=5
-edecomp-rank-by=abs-mean
"@

Run-TitanCase "mm" $base

Push-Location $repo
try {
  foreach ($file in $outputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "mm: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "edecomp_summary.dat" -Pattern "map_semantics: MM matrices are direct pairwise" -Quiet)) {
    throw "mm: summary does not document map semantics"
  }
  if (-not (Select-String -Path "edecomp_mm_atom_ranking.csv" -Pattern "abs_mean,mean_abs" -Quiet)) {
    throw "mm: atom ranking does not export absolute-value ranking fields"
  }
  if (Test-Path -LiteralPath "edecomp_mm_atom_matrix.csv") {
    throw "mm: dense atom matrix should not be exported"
  }
  Assert-NoLegacyEnergyOutputs "mm"

  foreach ($file in $gbsaOutputs) {
    if (Test-Path -LiteralPath $file) {
      throw "mm: GBSA output generated for pure MM source: $file"
    }
  }
  foreach ($file in $pbsaOutputs) {
    if (Test-Path -LiteralPath $file) {
      throw "mm: PBSA output generated for pure MM source: $file"
    }
  }
}
finally {
  Pop-Location
}

$gbsa = $base.Replace("-edecomp-source=mm", "-edecomp-source=gbsa") + "`r`n" + @"
-gbsa=true
"@
Run-TitanCase "gbsa" $gbsa

Push-Location $repo
try {
  foreach ($file in $gbsaOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "gbsa: missing expected output $file"
    }
  }
  if (Get-ChildItem -LiteralPath $repo -Filter "edecomp_mmgbsa_*matrix*.csv" -ErrorAction SilentlyContinue) {
    throw "gbsa: MM/GBSA generated a matrix output"
  }
  if (-not (Select-String -Path "edecomp_summary.dat" -Pattern "MM/GBSA maps are endpoint entity decompositions" -Quiet)) {
    throw "gbsa: summary does not document endpoint map semantics"
  }
  Assert-NoLegacyEnergyOutputs "gbsa"
}
finally {
  Pop-Location
}

$pbsaBase = $base.Replace("-edecomp-source=mm", "-edecomp-source=pbsa").Replace("-g=1-2", "-g=1&1") + "`r`n" + @"
-pbsa=true
-pbsa-components=complex,receptor,ligand
-pbsa-ie-min-samples=1
-pbsa-ie-block=1
-pbe-scale=0.35
-pbe-margin=2.0
-pbe-max-iter=80
-pbe-rms=1e-4
-pbe-dphi=1e-3
-pbe-max-grid=200000
-pbe-map=off
-nt=1
"@
Run-TitanCase "pbsa" $pbsaBase

Push-Location $repo
try {
  foreach ($file in $pbsaOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "pbsa: missing expected output $file"
    }
  }
  if (Get-ChildItem -LiteralPath $repo -Filter "edecomp_mmpbsa_*matrix*.csv" -ErrorAction SilentlyContinue) {
    throw "pbsa: MM/PBSA generated a matrix output"
  }
  if (-not (Select-String -Path "edecomp_summary.dat" -Pattern "MM/PBSA maps are endpoint entity decompositions" -Quiet)) {
    throw "pbsa: summary does not document endpoint map semantics"
  }
  if (-not (Select-String -Path "edecomp_mmpbsa_residue_map.csv" -Pattern "DeltaDeltaG_solv_PB" -Quiet)) {
    throw "pbsa: map does not use PBSA component labels"
  }
  Assert-NoLegacyEnergyOutputs "pbsa"
}
finally {
  Pop-Location
}

$all = $base.Replace("-edecomp-source=mm", "-edecomp-source=all") + "`r`n" + @"
-gbsa=true
"@
Run-TitanCase "all" $all

Push-Location $repo
try {
  foreach ($file in ($outputs + $gbsaOutputs)) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "all: missing expected output $file"
    }
  }
  foreach ($file in $pbsaOutputs) {
    if (Test-Path -LiteralPath $file) {
      throw "all: PBSA output generated without -pbsa=true: $file"
    }
  }
  Assert-NoLegacyEnergyOutputs "all"
}
finally {
  Pop-Location
}

$allPbsa = $base.Replace("-edecomp-source=mm", "-edecomp-source=all").Replace("-g=1-2", "-g=1&1") + "`r`n" + @"
-gbsa=true
-pbsa=true
-pbsa-components=complex,receptor,ligand
-pbsa-decomp=residue
-pbsa-ie-min-samples=1
-pbsa-ie-block=1
-pbe-scale=0.35
-pbe-margin=2.0
-pbe-max-iter=80
-pbe-rms=1e-4
-pbe-dphi=1e-3
-pbe-max-grid=200000
-pbe-map=off
-nt=1
"@
Run-TitanCase "all_pbsa" $allPbsa

Push-Location $repo
try {
  foreach ($file in ($outputs + $gbsaOutputs + $pbsaOutputs)) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "all_pbsa: missing expected output $file"
    }
  }
  Assert-NoLegacyEnergyOutputs "all_pbsa"
  Write-Host "PASS edecomp_smoke"
}
finally {
  Pop-Location
}
