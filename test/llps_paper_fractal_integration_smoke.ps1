param(
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\llps_paper_fractal_integration_smoke")
)

$ErrorActionPreference = "Stop"
$repo = Split-Path $PSScriptRoot -Parent
$full = [System.IO.Path]::GetFullPath($WorkRoot)
$allowed = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "tmp"))
if (-not $full.StartsWith($allowed, [System.StringComparison]::OrdinalIgnoreCase)) {
  throw "Refusing path outside test tmp root: $full"
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$exe = Join-Path $WorkRoot "llps_paper_fractal_integration_smoke.exe"
$Compiler = if ($env:CXX) { $env:CXX } else { "g++" }
& $Compiler -std=c++23 -O2 -DTITAN_LLPS_TESTING -ffunction-sections -fdata-sections `
  -I "$repo\src" `
  "$repo\test\llps_paper_fractal_integration_smoke.cpp" `
  "$repo\src\titan_llps_paper.cpp" `
  "-Wl,--gc-sections" -o $exe
if ($LASTEXITCODE -ne 0) { throw "Failed to compile fractal integration smoke" }
Push-Location $WorkRoot
try {
  & $exe
  if ($LASTEXITCODE -ne 0) { throw "Fractal integration smoke failed" }
} finally {
  Pop-Location
}
