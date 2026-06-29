param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\algebra_smoke")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

$exe = Join-Path $WorkRoot "titan_algebra_smoke.exe"
& $Compiler -std=c++23 -O2 -mavx2 -Isrc `
  test\titan_algebra_smoke.cpp `
  -o $exe
if ($LASTEXITCODE -ne 0) {
  throw "algebra_smoke: failed to compile"
}

& $exe
if ($LASTEXITCODE -ne 0) {
  throw "algebra_smoke: smoke executable failed"
}

Write-Host "PASS algebra_smoke"
