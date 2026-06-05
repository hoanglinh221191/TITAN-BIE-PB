param(
  [string]$TitanExe = ".\titan.exe",
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\feature_export_smoke")
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\feature_export_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

$unitExe = Join-Path $WorkRoot "titan_observable_features_smoke.exe"
& $Compiler -std=c++23 -O2 -Isrc -Isrc\eigen `
  test\titan_observable_features_smoke.cpp src\titan_observable_features.cpp `
  -o $unitExe
if ($LASTEXITCODE -ne 0) {
  throw "feature_export_smoke: failed to compile router smoke"
}
& $unitExe
if ($LASTEXITCODE -ne 0) {
  throw "feature_export_smoke: router smoke failed"
}

$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$outputs = @(
  "features_by_frame.csv",
  "features_metadata.csv",
  "feature_export_summary.dat",
  "features_by_frame_0.csv",
  "features_by_frame_1.csv",
  "features_metadata_0.csv",
  "features_metadata_1.csv",
  "feature_export_summary_0.dat",
  "feature_export_summary_1.dat",
  "dihedral_backbone_timeseries.csv",
  "dihedral_summary.dat"
)

function Clear-FeatureOutputs {
  foreach ($file in $outputs) {
    $path = Join-Path $repo $file
    if (Test-Path -LiteralPath $path) {
      Remove-Item -LiteralPath $path -Force
    }
  }
  Get-ChildItem -LiteralPath $repo -Filter "features_*_[0-9]*.csv" -ErrorAction SilentlyContinue |
    Remove-Item -Force
  Get-ChildItem -LiteralPath $repo -Filter "feature_export_summary*.dat" -ErrorAction SilentlyContinue |
    Remove-Item -Force
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-FeatureOutputs
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

function Expect-TitanReject([string]$Name, [string]$ArgText, [string]$ExpectedMessage) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  $stdout = Join-Path $work ("stdout_" + $Name + ".txt")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-FeatureOutputs
  Push-Location $repo
  $savedErrorActionPreference = $ErrorActionPreference
  try {
    $ErrorActionPreference = "Continue"
    & $titanPath -a $arg *> $stdout
    $exitCode = $LASTEXITCODE
  }
  finally {
    $ErrorActionPreference = $savedErrorActionPreference
    Pop-Location
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
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
"@

Run-TitanCase "dihedral_export" ($base + "`r`n" + @"
-feature-export=true
-feature-source=dihedral
"@)
Push-Location $repo
try {
  foreach ($file in @("features_by_frame.csv", "features_metadata.csv", "feature_export_summary.dat")) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "dihedral_export: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "feature_export_summary.dat" -Pattern "requested_source: dihedral" -Quiet) -or
      -not (Select-String -Path "feature_export_summary.dat" -Pattern "resolved_source: dihedral" -Quiet) -or
      -not (Select-String -Path "features_metadata.csv" -Pattern "^0,dihedral," -Quiet)) {
    throw "dihedral_export: feature outputs do not identify dihedral source"
  }
}
finally {
  Pop-Location
}

Run-TitanCase "qnative_export" ($base + "`r`n" + @"
-feature-export=true
-feature-source=qnative
-qnative-min-seqsep=0
"@)
Push-Location $repo
try {
  if (-not (Select-String -Path "features_metadata.csv" -Pattern "qnative,Q,fraction" -Quiet) -or
      -not (Select-String -Path "feature_export_summary.dat" -Pattern "feature_dim: 1" -Quiet)) {
    throw "qnative_export: qnative scalar feature was not exported"
  }
}
finally {
  Pop-Location
}

Expect-TitanReject "bad_source" ($base + "`r`n-feature-export=true`r`n-feature-source=bad") "Feature export: -feature-source must be one of"

Write-Host "PASS feature_export_smoke"
