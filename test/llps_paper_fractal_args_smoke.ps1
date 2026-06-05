param(
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\llps_paper_fractal_args_smoke")
)

$ErrorActionPreference = "Stop"
$repo = Split-Path $PSScriptRoot -Parent
$full = [System.IO.Path]::GetFullPath($WorkRoot)
$allowed = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "tmp"))
if (-not $full.StartsWith($allowed, [System.StringComparison]::OrdinalIgnoreCase)) {
  throw "Refusing path outside test tmp root: $full"
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$exe = Join-Path $WorkRoot "llps_paper_fractal_args_smoke.exe"
$Compiler = if ($env:CXX) { $env:CXX } else { "g++" }
& $Compiler -std=c++23 -O2 -I "$repo\src" `
  "$repo\test\llps_paper_fractal_args_smoke.cpp" `
  "$repo\src\parse_arguments.cpp" `
  -o $exe
if ($LASTEXITCODE -ne 0) { throw "Failed to compile fractal args smoke" }
& $exe
if ($LASTEXITCODE -ne 0) { throw "Fractal args smoke failed" }
Write-Host "PASS llps_paper_fractal_args_smoke"
