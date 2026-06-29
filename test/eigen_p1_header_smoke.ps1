param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\eigen_p1_header_smoke")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$srcDir = Join-Path $repo "src"
$smokeSource = Join-Path $repo "test\titan_eigen_p1_header_smoke.cpp"
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

foreach ($header in @("src\titan_rms.h", "src\titan_s2.h")) {
  $path = Join-Path $repo $header
  $eigenLeak = Select-String -LiteralPath $path -Pattern "Eigen" -SimpleMatch
  if ($eigenLeak) {
    throw "eigen_p1_header_smoke: $header leaked Eigen token"
  }
}

$obj = Join-Path $WorkRoot "titan_eigen_p1_header_smoke.o"
& $Compiler -std=c++23 -O2 -mavx2 "-I$srcDir" -c $smokeSource -o $obj
if ($LASTEXITCODE -ne 0) {
  throw "eigen_p1_header_smoke: header-only compile failed"
}

Write-Host "PASS eigen_p1_header_smoke"
