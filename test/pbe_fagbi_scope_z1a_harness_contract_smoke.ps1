param(
    [string]$ScriptPath = (Join-Path $PSScriptRoot 'pbe_fagbi_scope_z1a_parity_matrix.ps1')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if (-not (Test-Path -LiteralPath $ScriptPath)) {
    throw "Missing Scope Z1A parity harness: $ScriptPath"
}

$source = Get-Content -LiteralPath $ScriptPath -Raw
foreach ($required in @(
    'fagbi_scope_z1a_matrix.csv',
    'fagbi_scope_z1a_decision.csv',
    'fagbi_scope_z1a_summary.dat',
    'fagbi_scope_z1a_accuracy_plot.svg',
    'TITAN_FAGBI_SCOPE_Z1A_PARITY_VALIDATION_REPORT.txt',
    'z1a_parity_pass',
    'z1a_accuracy_gap',
    'z1a_runtime_gap',
    'blocked_profile_schema',
    'skipped_missing_fixture',
    'pbe_fmm_m2l_cli_helpers.ps1'
)) {
    if ($source -notlike "*$required*") {
        throw "Scope Z1A harness is missing contract token: $required"
    }
}

Write-Host 'pbe_fagbi_scope_z1a_harness_contract_smoke ok'
