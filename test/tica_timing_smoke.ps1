param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\tica_timing_smoke")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

$ticaTimingExe = Join-Path $WorkRoot "titan_tica_timing_smoke.exe"
& $Compiler -std=c++23 -O2 -Isrc -Isrc\eigen `
  test\titan_tica_timing_smoke.cpp src\titan_kinetics.cpp `
  -o $ticaTimingExe
if ($LASTEXITCODE -ne 0) {
  throw "tica_timing_smoke: failed to compile TICA timing smoke"
}
& $ticaTimingExe
if ($LASTEXITCODE -ne 0) {
  throw "tica_timing_smoke: TICA timing smoke failed"
}

Write-Host "PASS tica_timing_smoke"
