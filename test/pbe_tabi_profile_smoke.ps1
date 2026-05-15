param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_profile_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootFull = [System.IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    if (-not $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside TABI profile smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    $lines = @(
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}

$runDir = Join-Path $WorkRoot 'mini_protein_detail'
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
New-CleanDir $runDir

$pdbPath = Join-Path $runDir 'driver.pdb'
Write-DummyPdb $pdbPath
$pqrPath = Join-Path $InputRoot 'sphere.pqr'
if (-not (Test-Path -LiteralPath $pqrPath)) {
    throw "Missing sphere PQR fixture: $pqrPath"
}
$libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
$argPath = Join-Path $runDir 'arg.in'
$argLines = @(
    "-f=$(Convert-ToTitanPath $pdbPath)",
    '-1pdb',
    '-ff=AMBER99SB',
    '-g=1&1',
    "-path-a=$libPath",
    "-path-n=$libPath",
    '-pbe=tabi',
    "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
    '-pbe-surface=vdw',
    '-pbe-max-iter=160',
    '-pbe-rms=1e-4',
    '-pbe-dphi=1e-4',
    '-pbe-timing=detail',
    '-pbe-tabi-sum=direct',
    '-pbe-tabi-formulation=jcp2013',
    '-pbe-tabi-near-integration=duffy_adaptive',
    '-pbe-tabi-surface-points=24',
    '-pbe-tabi-gmres-restart=8',
    '-nt=1'
)
Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII

Push-Location $runDir
try {
    & $TitanExe -a $argPath | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE for TABI profile smoke"
    }
} finally {
    Pop-Location
}

$energyPath = Join-Path $runDir 'pbe_energy.csv'
$timingPath = Join-Path $runDir 'pbe_timing.csv'
$profilePath = Join-Path $runDir 'pbe_tabi_profile.csv'
$formulationAuditPath = Join-Path $runDir 'pbe_tabi_formulation_audit.csv'
$nearBlockAuditPath = Join-Path $runDir 'pbe_tabi_near_block_audit.csv'
foreach ($required in @($energyPath, $timingPath, $profilePath, $formulationAuditPath, $nearBlockAuditPath)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing TABI profile output: $required"
    }
}

$energy = @(Import-Csv -LiteralPath $energyPath)[0]
if ($energy.status -ne 'ok') {
    throw "TABI profile smoke did not finish ok: $($energy.status)"
}
if ($energy.reaction_energy_source -ne 'tabi_boundary_integral') {
    throw "TABI profile smoke used fallback: $($energy.reaction_energy_source)"
}

$timing = @(Import-Csv -LiteralPath $timingPath)[0]
foreach ($column in @(
    'tabi_rhs_seconds','tabi_gmres_seconds','tabi_matvec_seconds_total',
    'tabi_matvec_seconds_avg','tabi_matvec_calls',
    'tabi_orthogonalization_seconds','tabi_preconditioner_seconds',
    'tabi_surface_neighbor_seconds','tabi_surface_mesh_seconds',
    'tabi_solid_angle_seconds','tabi_tree_moment_seconds',
    'tabi_tree_traversal_seconds','tabi_near_eval_seconds',
    'tabi_singular_near_eval_seconds','tabi_leaf_near_eval_seconds',
    'tabi_leaf_panel_eval_seconds','tabi_leaf_tile_eval_seconds',
    'tabi_leaf_residual_eval_seconds','tabi_smooth_leaf_node_eval_seconds',
    'tabi_leaf_panel_kernel','tabi_leaf_tile_status',
    'tabi_leaf_panel_pairs_per_second',
    'tabi_singular_near_interactions','tabi_leaf_near_interactions',
    'tabi_far_node_interactions','tabi_leaf_panel_interactions',
    'tabi_leaf_tile_interactions','tabi_leaf_tile_panel_interactions',
    'tabi_leaf_residual_panel_interactions',
    'tabi_smooth_leaf_node_interactions','tabi_leaf_pruned_interactions',
    'tabi_leaf_near_policy',
    'tabi_near_integration','tabi_near_kernel_policy',
    'tabi_near_k1_accepted','tabi_near_k1_rejected',
    'tabi_near_k2_accepted','tabi_near_k2_rejected',
    'tabi_near_k3_accepted','tabi_near_k3_rejected',
    'tabi_near_k4_accepted','tabi_near_k4_rejected',
    'tabi_near_one_ring_k1_accepted','tabi_near_one_ring_k1_rejected',
    'tabi_near_one_ring_k2_accepted','tabi_near_one_ring_k2_rejected',
    'tabi_near_one_ring_k3_accepted','tabi_near_one_ring_k3_rejected',
    'tabi_near_one_ring_k4_accepted','tabi_near_one_ring_k4_rejected',
    'tabi_near_k1_policy','tabi_near_k2_policy',
    'tabi_near_k3_policy','tabi_near_k4_policy',
    'tabi_near_self_pairs',
    'tabi_near_one_ring_pairs','tabi_near_duffy_pairs',
    'tabi_near_centroid_pairs','tabi_near_fallbacks',
    'tabi_near_integration_seconds',
    'tabi_far_eval_seconds','tabi_status','tabi_convergence_reason',
    'tabi_reaction_energy_source','tabi_formulation',
    'tabi_energy_projection','tabi_fallback_used',
    'tabi_best_preconditioned_residual',
    'tabi_final_preconditioned_residual',
    'tabi_true_residual_abs','tabi_true_residual_rel',
    'tabi_true_residual_probe_mode',
    'tabi_true_residual_probe_scope',
    'tabi_true_residual_probe_interval',
    'tabi_true_residual_sample_count',
    'tabi_true_residual_probe_calls',
    'tabi_true_residual_probe_seconds',
    'tabi_sampled_true_residual_rel',
    'tabi_sampled_true_residual_rows',
    'tabi_sampled_probe_calls',
    'tabi_sampled_probe_seconds',
    'tabi_full_true_residual_confirm_calls',
    'tabi_full_probe_seconds',
    'tabi_interaction_mode','tabi_quad_buffer_auto_candidate',
    'tabi_unknown_count','tabi_bie_system','tabi_lhs_diag_phi',
    'tabi_lhs_diag_phin','tabi_energy_equation',
    'tabi_kernel_convention','tabi_source_norm_s1',
    'tabi_source_norm_s2','tabi_phi_norm','tabi_phin_norm',
    'tabi_boundary_energy_kJmol','tabi_reference_energy_kJmol')) {
    if (-not ($timing.PSObject.Properties.Name -contains $column)) {
        throw "pbe_timing.csv missing TABI detail column: $column"
    }
}

$profile = @(Import-Csv -LiteralPath $profilePath)
if ($profile.Count -ne 1) {
    throw "Expected one pbe_tabi_profile.csv row, found $($profile.Count)"
}
$row = $profile[0]
if ($row.status -ne 'ok' -or $row.reaction_energy_source -ne 'tabi_boundary_integral') {
    throw "Unexpected TABI profile status/source: status=$($row.status) source=$($row.reaction_energy_source)"
}
if ([int64]$row.matvec_calls -le 0) {
    throw "TABI profile did not record matvec calls."
}
foreach ($column in @('formulation','energy_projection','fallback_used',
        'near_integration','near_kernel_policy',
        'near_k1_accepted','near_k1_rejected',
        'near_k2_accepted','near_k2_rejected',
        'near_k3_accepted','near_k3_rejected',
        'near_k4_accepted','near_k4_rejected',
        'near_one_ring_k1_accepted','near_one_ring_k1_rejected',
        'near_one_ring_k2_accepted','near_one_ring_k2_rejected',
        'near_one_ring_k3_accepted','near_one_ring_k3_rejected',
        'near_one_ring_k4_accepted','near_one_ring_k4_rejected',
        'near_k1_policy','near_k2_policy',
        'near_k3_policy','near_k4_policy',
        'near_self_pairs','near_one_ring_pairs',
        'near_duffy_pairs','near_centroid_pairs','near_fallbacks',
        'near_integration_seconds',
        'best_preconditioned_residual','final_preconditioned_residual',
        'true_residual_abs','true_residual_rel',
        'true_residual_probe_mode','true_residual_probe_interval',
        'true_residual_probe_scope','true_residual_sample_count',
        'true_residual_probe_calls','true_residual_probe_seconds',
        'sampled_true_residual_rel','sampled_true_residual_rows',
        'sampled_probe_calls','sampled_probe_seconds',
        'full_true_residual_confirm_calls','full_probe_seconds',
        'interaction_mode','quad_buffer_auto_candidate',
        'singular_near_eval_seconds','leaf_near_eval_seconds',
        'leaf_panel_eval_seconds','leaf_tile_eval_seconds',
        'leaf_residual_eval_seconds','smooth_leaf_node_eval_seconds',
        'leaf_panel_kernel','leaf_tile_status','leaf_panel_pairs_per_second',
        'singular_near_interactions','leaf_near_interactions',
        'far_node_interactions','leaf_panel_interactions',
        'leaf_tile_interactions','leaf_tile_panel_interactions',
        'leaf_residual_panel_interactions',
        'smooth_leaf_node_interactions','leaf_pruned_interactions',
        'leaf_near_policy',
        'unknown_count','bie_system','lhs_diag_phi','lhs_diag_phin',
        'energy_equation','kernel_convention','source_norm_s1',
        'source_norm_s2','phi_norm','phin_norm',
        'boundary_energy_kJmol','reference_energy_kJmol',
        'selected_energy_scale')) {
    if (-not ($row.PSObject.Properties.Name -contains $column)) {
        throw "pbe_tabi_profile.csv missing formulation column: $column"
    }
}
$audit = @(Import-Csv -LiteralPath $formulationAuditPath)
if ($audit.Count -ne 1 -or $audit[0].reaction_energy_source -ne 'tabi_boundary_integral') {
    throw "Unexpected TABI formulation audit row."
}
$nearBlockAudit = @(Import-Csv -LiteralPath $nearBlockAuditPath)
if ($nearBlockAudit.Count -lt 4) {
    throw "Expected TABI near block audit rows, found $($nearBlockAudit.Count)"
}
foreach ($column in @('block','pair_type','orientation','accepted','rejected',
        'candidate_worse_than_centroid','reference_nonconverged',
        'orientation_sign_mismatch','strict_ref_reject',
        'max_candidate_ref_rel','mean_candidate_ref_rel',
        'max_centroid_ref_rel','k2_laplace_sign_ok',
        'k2_yukawa_smooth_ref_rel','k2_orientation_pair_ref_rel')) {
    if (-not ($nearBlockAudit[0].PSObject.Properties.Name -contains $column)) {
        throw "pbe_tabi_near_block_audit.csv missing column: $column"
    }
}
if ((Convert-ToDouble $row.gmres_seconds) -le 0.0) {
    throw "TABI profile did not record GMRES wall time."
}
if ([string]::IsNullOrWhiteSpace($row.residual_history)) {
    throw "TABI profile did not record residual history."
}

Write-Host ("PBE TABI profile smoke PASS: panels={0} matvec_calls={1} E_kJ={2}" -f
    $row.panels, $row.matvec_calls, $row.energy_kJmol)
