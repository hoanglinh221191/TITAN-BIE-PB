param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\hbond_kinetics_smoke")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

$hbkExe = Join-Path $WorkRoot "titan_hbond_kinetics_smoke.exe"
& $Compiler -std=c++23 -O2 -Isrc -Isrc\eigen `
  test\titan_hbond_kinetics_smoke.cpp `
  src\titan_hbond_kinetics.cpp `
  src\titan_interaction_primitives.cpp `
  src\titan_executor.cpp `
  -o $hbkExe
if ($LASTEXITCODE -ne 0) {
  throw "hbond_kinetics_smoke: failed to compile HBond kinetics smoke"
}
& $hbkExe
if ($LASTEXITCODE -ne 0) {
  throw "hbond_kinetics_smoke: HBond kinetics smoke failed"
}

Write-Host "PASS hbond_kinetics_smoke"
