param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4y_mini_contract')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$script = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4y_mini_dense_oracle.ps1'
$output = & powershell -NoProfile -ExecutionPolicy Bypass -File $script `
    -TitanExe $TitanExe `
    -WorkRoot $WorkRoot `
    -PlanOnly 2>&1 | Out-String

if ($LASTEXITCODE -ne 0) {
    throw "Z4Y mini dense plan-only contract failed:`n$output"
}

foreach ($fixture in @('1vii-mini','1ubq-mini','2ci2-mini')) {
    foreach ($surface in @(8,12)) {
        $direct = 'dense_oracle_{0}_s{1}_tabi_direct_dense' -f $fixture,$surface
        $p4 = 'dense_oracle_{0}_s{1}_fmm_p4_guarded' -f $fixture,$surface
        if ($output -notmatch [regex]::Escape($direct)) {
            throw "Z4Y mini dense oracle must plan $direct`n$output"
        }
        if ($output -notmatch [regex]::Escape($p4)) {
            throw "Z4Y mini dense oracle must plan $p4`n$output"
        }
    }
}

foreach ($forbidden in @(
    'dense_oracle_1vii_s8_tabi_direct_dense',
    'dense_oracle_1ubq_s8_tabi_direct_dense',
    'dense_oracle_2ci2_s8_tabi_direct_dense',
    'dense_oracle_sphere_coarse_s8_tabi_direct_dense',
    'dense_oracle_kirkwood_coarse_s8_tabi_direct_dense'
)) {
    if ($output -match [regex]::Escape($forbidden)) {
        throw "Z4Y mini dense wrapper must not plan non-mini or analytic rows: $forbidden`n$output"
    }
}

Write-Host 'PASS pbe_fagbi_scope_z4y_mini_dense_contract_smoke'
