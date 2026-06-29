param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\alignment_firewall_smoke")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$srcDir = Join-Path $repo "src"
$eigenDir = Join-Path $repo "src\eigen"
$smokeSource = Join-Path $repo "test\titan_alignment_firewall_smoke.cpp"
$adapterSource = Join-Path $repo "test\titan_alignment_adapter_smoke.cpp"
$alignmentSource = Join-Path $repo "src\titan_alignment.cpp"
$alignmentHeader = Join-Path $repo "src\titan_alignment.h"
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

$eigenLeak = Select-String -LiteralPath $alignmentHeader -Pattern "Eigen" -SimpleMatch
if ($eigenLeak) {
  throw "alignment_firewall_smoke: titan_alignment.h leaked Eigen token"
}

$exe = Join-Path $WorkRoot "titan_alignment_firewall_smoke.exe"

$headerProbeObj = Join-Path $WorkRoot "titan_alignment_header_probe.o"
& $Compiler -std=c++23 -O2 -mavx2 "-I$srcDir" -c $smokeSource -o $headerProbeObj
if ($LASTEXITCODE -ne 0) {
  throw "alignment_firewall_smoke: header firewall compile failed"
}

& $Compiler -std=c++23 -O2 -mavx2 "-I$srcDir" "-I$eigenDir" `
  $smokeSource `
  $alignmentSource `
  -o $exe
if ($LASTEXITCODE -ne 0) {
  throw "alignment_firewall_smoke: failed to compile"
}

& $exe
if ($LASTEXITCODE -ne 0) {
  throw "alignment_firewall_smoke: smoke executable failed"
}

$adapterExe = Join-Path $WorkRoot "titan_alignment_adapter_smoke.exe"
& $Compiler -std=c++23 -O2 -mavx2 "-I$srcDir" "-I$eigenDir" `
  $adapterSource `
  $alignmentSource `
  -o $adapterExe
if ($LASTEXITCODE -ne 0) {
  throw "alignment_firewall_smoke: adapter smoke failed to compile"
}

& $adapterExe
if ($LASTEXITCODE -ne 0) {
  throw "alignment_firewall_smoke: adapter smoke executable failed"
}

Write-Host "PASS alignment_firewall_smoke"
