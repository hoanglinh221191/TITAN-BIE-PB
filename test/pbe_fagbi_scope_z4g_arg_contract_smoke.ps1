param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4g_arg_contract')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

if (Test-Path -LiteralPath $WorkRoot) { Remove-Item -LiteralPath $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$pqr = Join-Path $WorkRoot 'full.pqr'
Set-Content -LiteralPath $pqr -Encoding ASCII -Value @(
    'ATOM      1 N    ALA A     1      0.000   0.000   0.000  -0.300 1.500',
    'END'
)
$argPath = Write-FmmCliArgFile -RunDir $WorkRoot -RepoRoot $RepoRoot -PqrPath $pqr `
    -Backend tabi -Mode treecode -FmmMode single-matvec-analytic -SurfacePoints 8 `
    -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
    -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
    -OpeningTheta 0.6 -AuditSampleCount 512 -Preconditioner fmm-nearfield-ilu0 `
    -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
    -NearIntegration analytic_singular -TaylorOrder 4
$lines = @(Get-Content -LiteralPath $argPath)
foreach ($expected in @(
        '-pbe-tabi-near-integration=analytic_singular',
        '-pbe-tabi-taylor-order=4',
        '-pbe-tabi-precond=fmm-nearfield-ilu0',
        '-pbe-tabi-precond-ilu-neighbor-cap=64')) {
    if ($lines -notcontains $expected) { throw "Missing Z4G arg policy: $expected" }
}
Write-Host 'PASS pbe_fagbi_scope_z4g_arg_contract_smoke'
