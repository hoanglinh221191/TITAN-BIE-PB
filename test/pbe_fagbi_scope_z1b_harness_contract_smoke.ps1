param(
    [string]$ScriptPath = (Join-Path $PSScriptRoot 'pbe_fagbi_scope_z1b_near_field_attribution.ps1'),
    [string]$HelperPath = (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

if (-not (Test-Path -LiteralPath $HelperPath)) {
    throw "Missing shared FAGBI validation helper: $HelperPath"
}
if (-not (Test-Path -LiteralPath $ScriptPath)) {
    throw "Missing Scope Z1B attribution harness: $ScriptPath"
}

$helperSource = Get-Content -LiteralPath $HelperPath -Raw
foreach ($required in @(
    'Get-FagbiBornEnergyKJ',
    'Get-FagbiKirkwoodEnergyKJ',
    'Write-FagbiBornPqr',
    'Write-FagbiKirkwoodPqr',
    'Get-FagbiIcosphereTargetPoints'
)) {
    if ($helperSource -notlike "*$required*") {
        throw "Shared FAGBI helper is missing contract token: $required"
    }
}

. $HelperPath
$targets = @(0..3 | ForEach-Object { Get-FagbiIcosphereTargetPoints -Level $_ })
if (($targets -join ',') -ne '20,80,320,1280') {
    throw "Unexpected icosphere target sequence: $($targets -join ',')"
}
$born = Get-FagbiBornEnergyKJ -Charge 1.0 -RadiusA 3.0 -SaltM 0.15
$kirkwood = Get-FagbiKirkwoodEnergyKJ -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0 -SaltM 0.15
if ([double]::IsNaN($born) -or [double]::IsInfinity($born) -or
    [double]::IsNaN($kirkwood) -or [double]::IsInfinity($kirkwood)) {
    throw 'Shared FAGBI helper returned non-finite oracle energy.'
}

$source = Get-Content -LiteralPath $ScriptPath -Raw
foreach ($required in @(
    'fagbi_scope_z1b_matrix.csv',
    'fagbi_scope_z1b_convergence.csv',
    'fagbi_scope_z1b_decision.csv',
    'fagbi_scope_z1b_summary.dat',
    'fagbi_scope_z1b_near_field_plot.svg',
    'TITAN_FMM_M2L_SCOPE_Z1B_NEAR_FIELD_ATTRIBUTION_REPORT.txt',
    'z1b_quadrature_dominant',
    'z1b_mesh_dominant',
    'z1b_mixed_quadrature_mesh',
    'z1b_near_method_regression',
    'z1b_attribution_inconclusive',
    'analytic_singular',
    'duffy_adaptive',
    'quad-buffer=adaptive',
    'actual_panels',
    'near_one_ring_k4_rejected'
)) {
    if ($source -notlike "*$required*") {
        throw "Scope Z1B harness is missing contract token: $required"
    }
}

Write-Host 'pbe_fagbi_scope_z1b_harness_contract_smoke ok'
