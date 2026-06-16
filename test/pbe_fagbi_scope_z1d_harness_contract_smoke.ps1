$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

function Require-Equal([string]$Label, [object]$Actual, [object]$Expected) {
    if ($Actual -ne $Expected) {
        throw "$Label expected '$Expected', got '$Actual'"
    }
}

$mapping = @(
    @{ target=5120; request=1281; actual=5120 },
    @{ target=20480; request=5121; actual=20480 },
    @{ target=81920; request=20481; actual=81920 }
)
foreach ($row in $mapping) {
    $resolved = Get-FagbiHighresSurfaceResolution -TargetPanels $row.target
    Require-Equal "Z1D target $($row.target) request" `
        $resolved.requested_points $row.request
    Require-Equal "Z1D target $($row.target) actual" `
        $resolved.expected_panels $row.actual
}

$quickCases = @(Get-FagbiZ1dCaseSpecs -Quick)
$fullCases = @(Get-FagbiZ1dCaseSpecs)
Require-Equal 'quick case count' $quickCases.Count 2
Require-Equal 'full case count' $fullCases.Count 8
Require-Equal 'quick fixtures' (@($quickCases.fixture | Sort-Object -Unique) -join ',') 'born_centered'
Require-Equal 'quick salts' (@($quickCases.salt_M | Sort-Object -Unique) -join ',') '0.15'
Require-Equal 'quick targets' (@($quickCases.target_panels | Sort-Object -Unique) -join ',') '5120,20480'
Require-Equal 'full fixtures' (@($fullCases.fixture | Sort-Object -Unique) -join ',') 'born_centered,kirkwood_offcenter'
Require-Equal 'full salts' (@($fullCases.salt_M | Sort-Object -Unique) -join ',') '0,0.15'
Require-Equal 'full targets' (@($fullCases.target_panels | Sort-Object -Unique) -join ',') '5120,20480'

$scriptPath = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z1d_high_resolution_surface.ps1'
if (-not (Test-Path -LiteralPath $scriptPath)) {
    throw "Missing Scope Z1D harness: $scriptPath"
}
$scriptText = Get-Content -LiteralPath $scriptPath -Raw
foreach ($token in @(
    'fagbi_scope_z1d_matrix.csv',
    'fagbi_scope_z1d_convergence.csv',
    'fagbi_scope_z1d_decision.csv',
    'fagbi_scope_z1d_summary.dat',
    'fagbi_scope_z1d_highres_plot.svg',
    'TITAN_FMM_M2L_SCOPE_Z1D_HIGH_RESOLUTION_REPORT.txt',
    'z1d_resolution_converged',
    'z1d_resolution_limited',
    'z1d_operator_plateau',
    'z1d_runtime_limited',
    '-pbe-tabi-highres-surface=true',
    'analytic_singular'
)) {
    if ($scriptText -notmatch [regex]::Escape($token)) {
        throw "Scope Z1D harness is missing required token: $token"
    }
}

Write-Host 'pbe_fagbi_scope_z1d_harness_contract_smoke ok'
