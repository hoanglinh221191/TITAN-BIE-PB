param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4y_contract')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$script = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4y_reviewer_hardening_matrix.ps1'
$output = & powershell -NoProfile -ExecutionPolicy Bypass -File $script `
    -TitanExe $TitanExe `
    -WorkRoot $WorkRoot `
    -Quick `
    -Stages dense_oracle `
    -PlanOnly 2>&1 | Out-String

if ($LASTEXITCODE -ne 0) {
    throw "Z4Y plan-only contract failed:`n$output"
}
if ($output -notmatch 'dense_oracle_1vii-mini_s8_tabi_direct_dense') {
    throw "Z4Y dense oracle must use 1vii-mini direct dense calibration:`n$output"
}
if ($output -match 'dense_oracle_1vii_s8_tabi_direct_dense') {
    throw "Z4Y dense oracle must not plan full 1VII direct dense in quick mode:`n$output"
}
if ($output -notmatch 'dense_oracle_sphere_coarse_s8_tabi_direct_dense' -or
    $output -notmatch 'dense_oracle_kirkwood_coarse_s8_tabi_direct_dense') {
    throw "Z4Y dense oracle must keep analytic dense calibration rows:`n$output"
}

Write-Host 'PASS pbe_fagbi_scope_z4y_contract_smoke'
