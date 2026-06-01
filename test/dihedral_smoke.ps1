param(
  [string]$TitanExe = ".\test\bin\native\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\dihedral_inputs_local"
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

$outputs = @(
  "dihedral_backbone_timeseries.csv",
  "dihedral_sidechain_timeseries.csv",
  "dihedral_ramachandran.csv",
  "dihedral_rotamer_populations.csv",
  "dihedral_rotamer_transitions.csv",
  "dihedral_summary.dat",
  "dihedral_feature_matrix.csv"
)

function Clear-DihedralOutputs {
  foreach ($file in $outputs) {
    $path = Join-Path $repo $file
    if (Test-Path -LiteralPath $path) {
      Remove-Item -LiteralPath $path -Force
    }
  }
  Get-ChildItem -LiteralPath $repo -Filter "dihedral_*_[0-9]*.*" -ErrorAction SilentlyContinue |
    Remove-Item -Force
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-DihedralOutputs
  & $titanPath -a $arg
  if ($LASTEXITCODE -ne 0) {
    throw "${Name}: TITAN exited with code $LASTEXITCODE"
  }
}

function Expect-TitanReject([string]$Name, [string]$ArgText, [string]$ExpectedMessage) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  $stdout = Join-Path $work ("stdout_" + $Name + ".txt")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-DihedralOutputs
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
  Run-TitanCase "backbone" ($base + "`r`n" + @"
-dihedral
-dihedral-kind=backbone
"@)
  foreach ($file in @("dihedral_backbone_timeseries.csv", "dihedral_ramachandran.csv", "dihedral_summary.dat")) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "backbone: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "dihedral_summary.dat" -Pattern "kind: backbone" -Quiet) -or
      -not (Select-String -Path "dihedral_backbone_timeseries.csv" -Pattern "phi_deg,psi_deg,omega_deg" -Quiet)) {
    throw "backbone: dihedral output does not match backbone contract"
  }
  if (Test-Path -LiteralPath "dihedral_feature_matrix.csv") {
    throw "backbone: feature matrix generated without -dihedral-export-tica"
  }

  Run-TitanCase "sidechain" ($base + "`r`n" + @"
-dihedral
-dihedral-kind=sidechain
"@)
  foreach ($file in @("dihedral_sidechain_timeseries.csv", "dihedral_rotamer_populations.csv", "dihedral_rotamer_transitions.csv", "dihedral_summary.dat")) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "sidechain: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "dihedral_sidechain_timeseries.csv" -Pattern "chi1_deg,chi2_deg" -Quiet) -or
      -not (Select-String -Path "dihedral_summary.dat" -Pattern "kind: sidechain" -Quiet)) {
    throw "sidechain: output does not match chi1/chi2 contract"
  }
  if (Test-Path -LiteralPath "dihedral_backbone_timeseries.csv") {
    throw "sidechain: backbone output generated for sidechain-only run"
  }

  Run-TitanCase "both_feature_matrix" ($base + "`r`n" + @"
-dihedral
-dihedral-kind=both
-dihedral-residues=all
-dihedral-stride=1
-dihedral-export-tica
"@)
  if (-not (Test-Path -LiteralPath "dihedral_feature_matrix.csv")) {
    throw "both_feature_matrix: missing dihedral_feature_matrix.csv"
  }
  if (-not (Select-String -Path "dihedral_feature_matrix.csv" -Pattern "_phi_sin" -Quiet) -or
      -not (Select-String -Path "dihedral_feature_matrix.csv" -Pattern "_omega_cos" -Quiet) -or
      -not (Select-String -Path "dihedral_feature_matrix.csv" -Pattern "_chi1_sin" -Quiet)) {
    throw "both_feature_matrix: missing sin/cos torsion columns"
  }

  Expect-TitanReject "bad_kind_rejected" ($base + "`r`n-dihedral`r`n-dihedral-kind=chi3") "Dihedral: -dihedral-kind must be backbone, sidechain, or both."
  Expect-TitanReject "bad_stride_rejected" ($base + "`r`n-dihedral`r`n-dihedral-stride=0") "Dihedral: -dihedral-stride must be >= 1."

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
-dihedral
-dihedral-export-tica
"@
  Run-TitanCase "multi" $pooled
  foreach ($file in @("dihedral_backbone_timeseries_0.csv", "dihedral_backbone_timeseries_1.csv", "dihedral_summary_0.dat", "dihedral_summary_1.dat")) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "multi: missing indexed output $file"
    }
  }

  Run-TitanCase "disabled" $base
  foreach ($file in @("dihedral_backbone_timeseries.csv", "dihedral_sidechain_timeseries.csv", "dihedral_ramachandran.csv", "dihedral_rotamer_populations.csv", "dihedral_rotamer_transitions.csv", "dihedral_summary.dat", "dihedral_feature_matrix.csv")) {
    if (Test-Path -LiteralPath $file) {
      throw "disabled: dihedral output generated without -dihedral: $file"
    }
  }

  Write-Host "PASS dihedral_smoke"
}
finally {
  Pop-Location
}
