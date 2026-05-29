param(
  [string]$TitanExe = ".\test\bin\native\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\cmap_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null
$singleFramePdb = Join-Path $repo "example\2\md1_short.pdb"
$multiFramePdb = Join-Path $work "md1_short_two_frames.pdb"
$frameOne = [System.IO.File]::ReadAllText($singleFramePdb)
$frameTwo = $frameOne.Replace("MODEL        1", "MODEL        2").Replace(
  "t= 100000.00000", "t= 100001.00000")
[System.IO.File]::WriteAllText(
  $multiFramePdb, $frameOne + "`r`n" + $frameTwo, [System.Text.Encoding]::ASCII)

$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$newOutputs = @(
  "contact_probability_matrix.csv",
  "contact_probability_edges.csv",
  "contact_probability_summary.dat",
  "contact_count_timeseries.csv",
  "contact_residue_occupancy.csv"
)
$legacyOutputs = @(
  "average_contact_map.dat",
  "out_sc.txt",
  "out_sc_time.dat"
)
$averageOutputs = @(
  "contact_probability_matrix_average.csv",
  "contact_probability_edges_average.csv",
  "contact_probability_summary_average.dat",
  "contact_residue_occupancy_average.csv"
)

function Clear-ContactOutputs {
  foreach ($file in ($newOutputs + $averageOutputs + $legacyOutputs)) {
    $path = Join-Path $repo $file
    if (Test-Path -LiteralPath $path) {
      Remove-Item -LiteralPath $path -Force
    }
  }
  Get-ChildItem -LiteralPath $repo -Filter "contact_probability_*_[0-9]*.*" -ErrorAction SilentlyContinue |
    Remove-Item -Force
  Get-ChildItem -LiteralPath $repo -Filter "contact_map_traj_*.dat" -ErrorAction SilentlyContinue |
    Remove-Item -Force
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-ContactOutputs
  & $titanPath -a $arg
  if ($LASTEXITCODE -ne 0) {
    throw "${Name}: TITAN exited with code $LASTEXITCODE"
  }
}

function Expect-TitanReject([string]$Name, [string]$ArgText, [string]$ExpectedMessage) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  $stdout = Join-Path $work ("stdout_" + $Name + ".txt")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-ContactOutputs
  $savedErrorActionPreference = $ErrorActionPreference
  try {
    $ErrorActionPreference = "Continue"
    & $titanPath -a $arg *> $stdout
    $exitCode = $LASTEXITCODE
  }
  finally {
    $ErrorActionPreference = $savedErrorActionPreference
  }
  if ($exitCode -eq 0) {
    throw "${Name}: TITAN unexpectedly accepted invalid arguments"
  }
  if (-not (Select-String -Path $stdout -Pattern $ExpectedMessage -Quiet)) {
    throw "${Name}: expected validation message was not emitted"
  }
}

$base = @"
-f=..\..\example\2\md1_short.pdb
-1pdb
-nocal
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
"@

Push-Location $repo
try {
  Run-TitanCase "default" ($base + "`r`n" + @"
-cmap
-cmap-min-seqsep=0
"@)
  foreach ($file in $newOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "default: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "contact_probability_summary.dat" -Pattern "mode: sidechain-com" -Quiet) -or
      -not (Select-String -Path "contact_probability_summary.dat" -Pattern "cutoff_A: 6.500000" -Quiet)) {
    throw "default: summary does not resolve sidechain-com at 6.5 A"
  }

  $interface = @"
-f=$multiFramePdb
-1pdb
-nocal
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
-n=..\..\example\2\index.ndx
-g=3-4
-cmap
-cmap-scope=interface
"@
  Run-TitanCase "interface" $interface
  if (-not (Select-String -Path "contact_probability_summary.dat" -Pattern "scope: interface" -Quiet) -or
      -not (Select-String -Path "contact_probability_summary.dat" -Pattern "matrix_layout: rectangular" -Quiet) -or
      -not (Select-String -Path "contact_probability_summary.dat" -Pattern "frames_seen: 1" -Quiet) -or
      -not (Select-String -Path "contact_probability_summary.dat" -Pattern "valid_frames: 1" -Quiet) -or
      -not (Select-String -Path "contact_probability_summary.dat" -Pattern "skipped_incomplete_frames: 0" -Quiet)) {
    throw "interface: summary does not describe a populated rectangular output"
  }
  if (-not (Select-String -Path "contact_probability_matrix.csv" -Pattern "mol2\|" -Quiet)) {
    throw "interface: matrix does not contain molecule 2 labels"
  }

  $legacy = $base + "`r`n" + @"
-n=..\..\example\2\index.ndx
-g=3-4
-contact-map=yes
-dsc=7.0
"@
  Run-TitanCase "legacy" $legacy
  if (-not (Select-String -Path "contact_probability_summary.dat" -Pattern "mode: sidechain-com" -Quiet) -or
      -not (Select-String -Path "contact_probability_summary.dat" -Pattern "scope: interface" -Quiet) -or
      -not (Select-String -Path "contact_probability_summary.dat" -Pattern "cutoff_A: 7.000000" -Quiet)) {
    throw "legacy: deprecated alias did not map to interface COM CMap"
  }
  foreach ($file in $legacyOutputs) {
    if (Test-Path -LiteralPath $file) {
      throw "legacy: old contact output still generated: $file"
    }
  }

  Expect-TitanReject "interface_domain_rejected" ($interface + "`r`n-cmap-domain=both") "CMap: -cmap-domain and -cmap-min-seqsep apply only to -cmap-scope=mol1."
  Expect-TitanReject "mixed_alias_rejected" ($legacy + "`r`n-cmap") "CMap: do not mix deprecated -contact-map with -cmap options"

  $listFile = Join-Path $work "same_pdb_twice.txt"
  [System.IO.File]::WriteAllText(
    $listFile, $multiFramePdb + "`r`n" + $multiFramePdb + "`r`n",
    [System.Text.Encoding]::ASCII)
  $pooled = @"
-f=$multiFramePdb
-i=$listFile
-1pdb
-nocal
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
-cmap
-cmap-min-seqsep=0
"@
  Run-TitanCase "pooled" $pooled
  foreach ($file in @("contact_probability_matrix_0.csv", "contact_probability_matrix_1.csv") + $averageOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "pooled: missing expected output $file"
    }
  }
  $matrix0 = ([System.IO.File]::ReadAllText(
      (Join-Path $repo "contact_probability_matrix_0.csv"))) -replace "mol1\|[^,`r`n]+", "mol1|residue"
  $matrix1 = ([System.IO.File]::ReadAllText(
      (Join-Path $repo "contact_probability_matrix_1.csv"))) -replace "mol1\|[^,`r`n]+", "mol1|residue"
  if ($matrix0 -ne $matrix1) {
    throw "pooled: identical PDB inputs did not preserve CMap semantics across reload"
  }
  $matrixAverage = ([System.IO.File]::ReadAllText(
      (Join-Path $repo "contact_probability_matrix_average.csv"))) -replace "mol1\|[^,`r`n]+", "mol1|residue"
  if ($matrix0 -ne $matrixAverage) {
    throw "pooled: identical PDB inputs did not aggregate to their shared CMap"
  }

  Run-TitanCase "disabled" $base
  foreach ($file in $newOutputs) {
    if (Test-Path -LiteralPath $file) {
      throw "disabled: CMap output generated without -cmap or alias: $file"
    }
  }

  Write-Host "PASS cmap_smoke"
}
finally {
  Pop-Location
}
