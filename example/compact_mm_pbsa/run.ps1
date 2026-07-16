param(
    [Parameter(Mandatory=$true)] [string]$TitanExe,
    [string]$ResultsDir = (Join-Path $PSScriptRoot 'results')
)
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan.exe: $TitanExe" }
if (Test-Path -LiteralPath $ResultsDir) { throw "Results directory already exists: $ResultsDir" }
New-Item -ItemType Directory -Path $ResultsDir | Out-Null
Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'arg.in') -Destination $ResultsDir
Push-Location $ResultsDir
try {
    & (Resolve-Path -LiteralPath $TitanExe).Path -a .\arg.in
    if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE" }
} finally {
    Pop-Location
}
foreach ($required in @('mm_pbsa_energy.csv','mm_pbsa_summary.dat','pbe_tabi_profile.csv')) {
    if (-not (Test-Path -LiteralPath (Join-Path $ResultsDir $required))) {
        throw "Missing expected output: $required"
    }
}
Write-Host "Compact MM/PBSA example completed: $ResultsDir"
