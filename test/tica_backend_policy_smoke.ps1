param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\tica_backend_policy_smoke")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

$policyExe = Join-Path $WorkRoot "titan_tica_backend_policy_smoke.exe"
& $Compiler -std=c++23 -O2 -Isrc -Isrc\eigen `
  test\titan_tica_backend_policy_smoke.cpp `
  -o $policyExe
if ($LASTEXITCODE -ne 0) {
  throw "tica_backend_policy_smoke: failed to compile TICA backend policy smoke"
}
& $policyExe
if ($LASTEXITCODE -ne 0) {
  throw "tica_backend_policy_smoke: TICA backend policy smoke failed"
}

Write-Host "PASS tica_backend_policy_smoke"
