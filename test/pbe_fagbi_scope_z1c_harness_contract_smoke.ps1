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
    @{ mesh='icosphere'; target=80; request=80; actual=80 },
    @{ mesh='icosphere'; target=320; request=320; actual=320 },
    @{ mesh='icosphere'; target=1280; request=1280; actual=1280 },
    @{ mesh='icosphere'; target=5120; request=1281; actual=5120 },
    @{ mesh='latlong'; target=80; request=16; actual=80 },
    @{ mesh='latlong'; target=320; request=81; actual=360 },
    @{ mesh='latlong'; target=1280; request=289; actual=1224 },
    @{ mesh='latlong'; target=5120; request=1225; actual=5040 }
)
foreach ($row in $mapping) {
    $resolved = Get-FagbiSurfaceResolution `
        -MeshMode $row.mesh `
        -TargetPanels $row.target
    Require-Equal "$($row.mesh) target $($row.target) request" `
        $resolved.requested_points $row.request
    Require-Equal "$($row.mesh) target $($row.target) actual" `
        $resolved.expected_panels $row.actual
}

$quickCases = @(Get-FagbiZ1cCaseSpecs -Quick)
$fullCases = @(Get-FagbiZ1cCaseSpecs)
Require-Equal 'quick case count' $quickCases.Count 8
Require-Equal 'full case count' $fullCases.Count 32
Require-Equal 'quick salts' (@($quickCases.salt_M | Sort-Object -Unique) -join ',') '0.15'
Require-Equal 'quick targets' (@($quickCases.target_panels | Sort-Object -Unique) -join ',') '320,1280'
Require-Equal 'full salts' (@($fullCases.salt_M | Sort-Object -Unique) -join ',') '0,0.15'
Require-Equal 'full targets' (@($fullCases.target_panels | Sort-Object -Unique) -join ',') '80,320,1280,5120'

$scriptPath = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z1c_surface_discretization.ps1'
if (-not (Test-Path -LiteralPath $scriptPath)) {
    throw "Missing Scope Z1C harness: $scriptPath"
}
$scriptText = Get-Content -LiteralPath $scriptPath -Raw
foreach ($token in @(
    'fagbi_scope_z1c_matrix.csv',
    'fagbi_scope_z1c_convergence.csv',
    'fagbi_scope_z1c_decision.csv',
    'fagbi_scope_z1c_summary.dat',
    'fagbi_scope_z1c_geometry_energy_plot.svg',
    'TITAN_FMM_M2L_SCOPE_Z1C_SURFACE_DISCRETIZATION_REPORT.txt',
    'resolution_converged',
    'resolution_limited',
    'operator_discretization_plateau',
    'mesh_quality_dominant',
    'analytic_singular'
)) {
    if ($scriptText -notmatch [regex]::Escape($token)) {
        throw "Scope Z1C harness is missing required token: $token"
    }
}
if ($scriptText -notmatch 'TITAN_LIB_FF') {
    throw 'Scope Z1C harness must support TITAN_LIB_FF for linked worktrees.'
}

$writerText = Get-Content -LiteralPath (Join-Path $RepoRoot 'src\tabi\analyzer_impl.cpp') -Raw
foreach ($column in @(
    'planar_area',
    'planar_area_rel_error_vs_4piR2',
    'enclosed_volume',
    'enclosed_volume_rel_error_vs_4piR3_over3',
    'triangle_area_cv'
)) {
    if ($writerText -notmatch [regex]::Escape($column)) {
        throw "Surface-quality CSV schema is missing column: $column"
    }
}

Write-Host 'pbe_fagbi_scope_z1c_harness_contract_smoke ok'
