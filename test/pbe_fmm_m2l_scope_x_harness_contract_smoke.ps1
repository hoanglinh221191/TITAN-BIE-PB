param(
    [string]$ScriptPath = (Join-Path $PSScriptRoot 'pbe_fmm_m2l_scope_x_sparse_ilu_probe.ps1')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if (-not (Test-Path -LiteralPath $ScriptPath)) {
    throw "Missing Scope X harness: $ScriptPath"
}

$source = Get-Content -Raw -LiteralPath $ScriptPath
$requiredMarkers = @(
    'Select-ScopeXLongCandidates',
    'Invoke-ScopeXDirectVerification',
    'Get-ScopeXRecoveredRunResult',
    '[switch]$Resume',
    '@($candidateRows | Where-Object',
    "`$_.stage -eq 'long'",
    'residual_improvement_vs_baseline',
    'energy_delta_abs',
    'speedup_vs_direct',
    'memory_ratio_vs_baseline',
    'direct_verification_status',
    'strict_verified',
    'strict_fmm_pending_direct',
    'useful_not_strict',
    'memory_limited',
    'factor_limited',
    'no_useful_improvement'
)

foreach ($marker in $requiredMarkers) {
    if (-not $source.Contains($marker)) {
        throw "Scope X.1 harness contract missing marker: $marker"
    }
}

Write-Host 'pbe_fmm_m2l_scope_x_harness_contract_smoke ok'
