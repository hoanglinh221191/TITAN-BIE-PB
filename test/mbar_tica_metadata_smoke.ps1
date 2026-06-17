param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\mbar_tica_metadata_smoke")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

$ticaExe = Join-Path $WorkRoot "titan_tica_estimator_smoke.exe"
& $Compiler -std=c++23 -O2 -Isrc `
  test\titan_tica_estimator_smoke.cpp src\titan_tica_estimator.cpp `
  -o $ticaExe
if ($LASTEXITCODE -ne 0) {
  throw "mbar_tica_metadata_smoke: failed to compile TICA estimator smoke"
}
& $ticaExe
if ($LASTEXITCODE -ne 0) {
  throw "mbar_tica_metadata_smoke: TICA estimator smoke failed"
}

$mbarExe = Join-Path $WorkRoot "titan_mbar_weights_smoke.exe"
& $Compiler -std=c++23 -O2 -Isrc `
  test\titan_mbar_weights_smoke.cpp src\titan_mbar_weights.cpp `
  -o $mbarExe
if ($LASTEXITCODE -ne 0) {
  throw "mbar_tica_metadata_smoke: failed to compile MBAR weights smoke"
}
& $mbarExe
if ($LASTEXITCODE -ne 0) {
  throw "mbar_tica_metadata_smoke: MBAR weights smoke failed"
}

Write-Host "PASS mbar_tica_metadata_smoke"
