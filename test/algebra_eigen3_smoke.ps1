param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\algebra_eigen3_smoke")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$srcDir = Join-Path $repo "src"
$smokeSource = Join-Path $repo "test\titan_algebra_eigen3_smoke.cpp"
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

$exe = Join-Path $WorkRoot "titan_algebra_eigen3_smoke.exe"
& $Compiler -std=c++23 -O2 -mavx2 "-I$srcDir" $smokeSource -o $exe
if ($LASTEXITCODE -ne 0) {
  throw "algebra_eigen3_smoke: failed to compile"
}

& $exe
if ($LASTEXITCODE -ne 0) {
  throw "algebra_eigen3_smoke: smoke executable failed"
}

Write-Host "PASS algebra_eigen3_smoke"
