param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_topdown_physical_gap_matrix'),
    [string[]]$Fixtures = @('mini-protein','charged-peptide','compact-mini-protein'),
    [int]$SurfacePoints = 48,
    [int]$TreeLeafSize = 128,
    [double]$OpeningTheta = 0.7,
    [int]$Threads = 8,
    [int]$CaseTimeoutSeconds = 0,
    [int]$MatvecAuditSampleCount = 1024,
    [switch]$IncludeExpansionAuditRows,
    [switch]$IncludeStage3HRows,
    [switch]$IncludeStage3IRows,
    [switch]$IncludeStage3JRows,
    [switch]$IncludeStage3KRows,
    [switch]$IncludeStage3LRows,
    [switch]$IncludeStage3MRows,
    [switch]$IncludeStage3NRows,
    [switch]$IncludeStage3ORows,
    [switch]$IncludeStage3PRows,
    [switch]$IncludeStage3QRows,
    [switch]$IncludeStage3RRows,
    [switch]$IncludeStage3SRows,
    [switch]$IncludeStage3TRows,
    [switch]$IncludeStage3URows,
    [switch]$IncludeStage3WRows,
    [switch]$IncludeStage3XRows,
    [switch]$IncludeStage3YRows,
    [switch]$IncludeStage3ZRows,
    [switch]$IncludeStage3AARows,
    [switch]$IncludeStage3ABRows,
    [switch]$IncludeStage3AIRows,
    [switch]$IncludeStage3AJRows,
    [switch]$IncludeStage3AKRows,
    [switch]$IncludeStage3ACRows,
    [switch]$IncludeStage3ADRows,
    [switch]$IncludeStage3AERows,
    [switch]$IncludeStage3AGRows,
    [switch]$IncludeStage3AHRows,
    [switch]$IncludeStage3ALRows,
    [int]$ExpansionAuditSampleCount = 512,
    [int]$ExpansionAuditMaxPanels = 2048,
    [switch]$Quick,
    [switch]$AllowFailedStatus
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

$FixtureMap = @{
    'mini-protein' = 'mini-protein.pqr'
    'charged-peptide' = 'charged-peptide.pqr'
    'compact-mini-protein' = 'compact-mini-protein.pqr'
}

$TreeRows = @(
    @{ CaseId = 'grouped_morton_off'; TreeBuilder = 'grouped-morton'; MacSafety = 'off'; ExpansionCenter = 'legacy'; FragmentationSafety = 'off'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ExpansionAudit = 'off' },
    @{ CaseId = 'topdown_aabb_off'; TreeBuilder = 'topdown-aabb'; MacSafety = 'off'; ExpansionCenter = 'legacy'; FragmentationSafety = 'off'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ExpansionAudit = 'off' },
    @{ CaseId = 'topdown_aabb_audit'; TreeBuilder = 'topdown-aabb'; MacSafety = 'audit'; ExpansionCenter = 'legacy'; FragmentationSafety = 'off'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ExpansionAudit = 'off' },
    @{ CaseId = 'topdown_aabb_physical_gap_legacy_off'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'legacy'; FragmentationSafety = 'off'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ExpansionAudit = 'off' },
    @{ CaseId = 'topdown_aabb_physical_gap_area_centroid_audit'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'audit'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ExpansionAudit = 'off' },
    @{ CaseId = 'topdown_aabb_physical_gap_area_centroid_split'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ExpansionAudit = 'off' }
)
if ($IncludeExpansionAuditRows) {
    $TreeRows += @{ CaseId = 'topdown_aabb_physical_gap_area_centroid_split_expansion_audit'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ExpansionAudit = 'sampled' }
}
if ($IncludeStage3HRows) {
    $TreeRows += @{ CaseId = 'stage3h_area_split_translated_naive'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ExpansionAudit = 'sampled' }
    $TreeRows += @{ CaseId = 'stage3h_area_split_translated_kahan'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'kahan'; ConvergenceAudit = 'sampled'; ExpansionAudit = 'sampled' }
    $TreeRows += @{ CaseId = 'stage3h_area_split_direct_node_kahan'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'direct-node'; MomentSum = 'kahan'; ConvergenceAudit = 'sampled'; ExpansionAudit = 'sampled' }
}
if ($IncludeStage3IRows) {
    $TreeRows += @{ CaseId = 'stage3i_o3_theta0p7_baseline'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'off'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.7 }
    $TreeRows += @{ CaseId = 'stage3i_o3_theta0p5'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'off'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.5 }
    $TreeRows += @{ CaseId = 'stage3i_o3_theta0p4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'off'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.4 }
    $TreeRows += @{ CaseId = 'stage3i_o4_theta0p5'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'off'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.5 }
    $TreeRows += @{ CaseId = 'stage3i_o4_theta0p6'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'off'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3i_o3_theta0p7_rho_cap0p5'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.5; MacAspectCap = 0.0; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.7 }
    $TreeRows += @{ CaseId = 'stage3i_o3_theta0p7_aspect_cap3'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 3.0; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.7 }
    $TreeRows += @{ CaseId = 'stage3i_o3_theta0p7_rho0p5_aspect3'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.5; MacAspectCap = 3.0; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.7 }
}
if ($IncludeStage3JRows) {
    $TreeRows += @{ CaseId = 'stage3j_center_area_o3_theta0p6'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'off'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'off'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3j_center_aabb_o3_theta0p6'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'aabb-center'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'off'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'off'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3j_o3_theta0p6_rho_cap0p35'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.35; MacAspectCap = 0.0; ConvergenceEstimator = 'off'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3j_o3_theta0p6_rho_cap0p30'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.30; MacAspectCap = 0.0; ConvergenceEstimator = 'off'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3j_o3_theta0p6_rho_cap0p25'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.25; MacAspectCap = 0.0; ConvergenceEstimator = 'off'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3j_o3_theta0p6_hybrid_t1e3'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceTailThreshold = 1.0e3; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3j_o3_theta0p6_hybrid_t1e4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3j_o4_theta0p6_hybrid_t1e4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.6 }
}
if ($IncludeStage3KRows) {
    $TreeRows += @{ CaseId = 'stage3k_o3_theta0p6_hybrid_stage3j_t1e4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'stage3j'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3k_o3_theta0p6_hybrid_consecutive_t1e4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3k_o3_theta0p6_hybrid_aspect_t1e4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'aspect'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3k_o3_theta0p6_hybrid_consecutive_aspect_t1e4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3k_o3_theta0p6_hybrid_consecutive_aspect_cap3e3'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 3.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3k_o4_theta0p6_hybrid_consecutive_aspect_t1e4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.6 }
}
if ($IncludeStage3LRows) {
    $TreeRows += @{ CaseId = 'stage3l_o3_theta0p6_hybrid_consecutive_aspect_dielectric_off_cluster_off'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'off'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3l_o3_theta0p6_hybrid_consecutive_aspect_dielectric_on_cluster_off'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3l_o3_theta0p5_hybrid_consecutive_aspect_dielectric_on_cluster_off'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.5 }
    $TreeRows += @{ CaseId = 'stage3l_o3_theta0p6_hybrid_consecutive_aspect_dielectric_on_leaf_pair_strict_b4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ClusterPolicy = 'leaf-pair-strict'; TargetBatchSize = 4; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3l_o3_theta0p6_hybrid_consecutive_aspect_dielectric_on_leaf_pair_balanced_b4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ClusterPolicy = 'leaf-pair-balanced'; TargetBatchSize = 4; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3l_o4_theta0p6_hybrid_consecutive_aspect_dielectric_on_cluster_off'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.6 }
}
if ($IncludeStage3MRows) {
    $TreeRows += @{ CaseId = 'stage3m_stage3l_o3_theta0p6_baseline'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'stage3l'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3m_o3_theta0p6_channel_audit'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'audit'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3m_o3_theta0p6_channel_tail'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3m_o3_theta0p6_channel_tail_eps'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3m_o3_theta0p5_channel_tail_eps'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.5 }
    $TreeRows += @{ CaseId = 'stage3m_o4_theta0p6_channel_tail_eps_fallback'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.6 }
}
if ($IncludeStage3NRows) {
    $TreeRows += @{ CaseId = 'stage3n_o3_theta0p6_k3_sqrt_eps'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'k3-sqrt-eps'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3n_o3_theta0p5_k3_sqrt_eps'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'k3-sqrt-eps'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.5 }
    $TreeRows += @{ CaseId = 'stage3n_o3_theta0p6_k3_eps_strict'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'k3-eps'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3n_o3_theta0p6_all_derivative_sqrt_eps'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'all-derivative-sqrt-eps'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3n_o4_theta0p6_k3_sqrt_eps_fallback'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'k3-sqrt-eps'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.6 }
}
if ($IncludeStage3ORows) {
    $TreeRows += @{ CaseId = 'stage3o_baseline_o4_theta0p6_k3_sqrt_eps_legacy'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'k3-sqrt-eps'; ChannelMissedSemantics = 'legacy'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3o_o4_theta0p6_all_derivative_sqrt_eps_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'all-derivative-sqrt-eps'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3o_o4_theta0p5_all_derivative_sqrt_eps_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'all-derivative-sqrt-eps'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.5 }
    $TreeRows += @{ CaseId = 'stage3o_o3_theta0p6_all_derivative_sqrt_eps_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'all-derivative-sqrt-eps'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.6 }
    $TreeRows += @{ CaseId = 'stage3o_o4_theta0p6_k4_eps_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'k4-eps'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.6 }
}
if ($IncludeStage3PRows) {
    $TreeRows += @{ CaseId = 'stage3p_o3_theta0p50_polynomial_eps_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.50 }
    $TreeRows += @{ CaseId = 'stage3p_o3_theta0p45_polynomial_eps_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.45 }
    $TreeRows += @{ CaseId = 'stage3p_o3_theta0p40_polynomial_eps_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.40 }
    $TreeRows += @{ CaseId = 'stage3p_o4_theta0p60_polynomial_eps_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
}
if ($IncludeStage3QRows) {
    $TreeRows += @{ CaseId = 'stage3q_o3_theta0p50_polynomial_eps_global_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.50 }
    $TreeRows += @{ CaseId = 'stage3q_o3_theta0p45_polynomial_eps_global_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.45 }
    $TreeRows += @{ CaseId = 'stage3q_o3_theta0p40_polynomial_eps_global_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.40 }
    $TreeRows += @{ CaseId = 'stage3q_o4_theta0p60_polynomial_eps_global_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3q_o4_theta0p60_polynomial_eps_k3_hot_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'k3-hot'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3q_o3_theta0p45_polynomial_eps_k3_hot_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'k3-hot'; ChannelMissedSemantics = 'accepted-only'; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 3; Theta = 0.45 }
}
if ($IncludeStage3RRows) {
    $TreeRows += @{ CaseId = 'stage3r_baseline_o4_theta0p60_polynomial_eps_global_accepted_only'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'center'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3r_o4_theta0p60_k3_rho_cap0p58'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; K3AcceptedRhoCap = 0.58; MacDistancePolicy = 'center'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3r_o4_theta0p60_k3_rho_cap0p57'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; K3AcceptedRhoCap = 0.57; MacDistancePolicy = 'center'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3r_o4_theta0p60_k3_rho_cap0p55'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; K3AcceptedRhoCap = 0.55; MacDistancePolicy = 'center'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3r_o4_theta0p60_k3_rho_cap0p52'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; K3AcceptedRhoCap = 0.52; MacDistancePolicy = 'center'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3r_o4_theta0p60_surface_margin0p25'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3r_o4_theta0p60_k3_rho_cap0p55_surface_margin0p25'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; K3AcceptedRhoCap = 0.55; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
}
if ($IncludeStage3SRows) {
    $TreeRows += @{ CaseId = 'stage3s_candidate_o4_theta0p60_surface_margin0p25_poly_abs1e7'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3s_candidate_o4_theta0p60_surface_margin0p25_poly_abs1e8'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-8; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3s_baseline_o4_theta0p60_center_poly_abs1e7'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'center'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
}
if ($IncludeStage3TRows) {
    $TreeRows += @{ CaseId = 'stage3t_promotion_readiness_o4_theta0p60_surface_margin0p25_poly_abs1e7'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
}
if ($IncludeStage3URows) {
    $TreeRows += @{ CaseId = 'stage3u_promotion_freeze_o4_theta0p60_surface_margin0p25_poly_abs1e7_cluster_off'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3u_timing_o4_theta0p60_surface_margin0p25_poly_abs1e7_leaf_pair_balanced_b4'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'leaf-pair-balanced'; TargetBatchSize = 4; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3u_floor_sensitivity_o4_theta0p60_surface_margin0p25_poly_abs1e8_cluster_off'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-8; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; ClusterPolicy = 'off'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
}
if ($IncludeStage3WRows) {
    $TreeRows += @{ CaseId = 'stage3w_safe_baseline_cluster_off'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; LeafPanelKernel = 'auto'; ClusterPolicy = 'off'; LeafPairKernel = 'auto'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3w_perf_isolate_cluster_off_tile_avx2'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'off'; LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'off'; LeafPairKernel = 'target-avx2'; TargetBatchSize = 4; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'off'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3w_perf_candidate_leaf_pair_balanced_b4_min64'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'off'; LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2'; TargetBatchSize = 4; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'off'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3w_perf_validation_leaf_pair_balanced_b4_min64_audit'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'sampled'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'sampled'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'summary'; LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2'; TargetBatchSize = 4; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'sampled'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3w_perf_candidate_leaf_pair_balanced_b4_min32'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'off'; LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2'; TargetBatchSize = 4; LeafPairMinSourcePanels = 32; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'off'; TaylorOrder = 4; Theta = 0.60 }
    $TreeRows += @{ CaseId = 'stage3w_perf_candidate_leaf_pair_balanced_b4_min128'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = 'off'; LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2'; TargetBatchSize = 4; LeafPairMinSourcePanels = 128; LeafPairMinTargets = 4; ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25; ExpansionAudit = 'off'; TaylorOrder = 4; Theta = 0.60 }
}
if ($IncludeStage3XRows) {
    function Add-Stage3XRow([string]$CaseId, [string]$ProbeCache, [string]$TargetGroup, [int]$Lookahead, [string]$ExpansionAudit = 'off', [string]$ExpansionAuditChannel = 'off', [string]$MomentAudit = 'off', [string]$ConvergenceAudit = 'off') {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split';
            MomentAudit = $MomentAudit; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = $ConvergenceAudit;
            ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid';
            ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only';
            AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0;
            MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = $ExpansionAuditChannel;
            LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2';
            TargetBatchSize = 4; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; LeafPairSelectionProfile = 'summary';
            LeafPairProbeCache = $ProbeCache; LeafPairTargetGroup = $TargetGroup; LeafPairMaxLookahead = $Lookahead;
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25;
            ExpansionAudit = $ExpansionAudit; TaylorOrder = 4; Theta = 0.60
        }
    }
    Add-Stage3XRow 'stage3x_perf_baseline_scan_no_cache_lh32' 'off' 'scan' 32
    Add-Stage3XRow 'stage3x_probe_cache_per_matvec_lh32' 'per-matvec' 'scan' 32
    Add-Stage3XRow 'stage3x_same_target_leaf_lh32' 'off' 'same-target-leaf' 32
    Add-Stage3XRow 'stage3x_cache_same_leaf_lh32' 'per-matvec' 'same-target-leaf' 32
    Add-Stage3XRow 'stage3x_cache_same_leaf_lh16' 'per-matvec' 'same-target-leaf' 16
    Add-Stage3XRow 'stage3x_cache_same_leaf_lh64' 'per-matvec' 'same-target-leaf' 64
    Add-Stage3XRow 'stage3x_validation_cache_same_leaf_lh32_audit' 'per-matvec' 'same-target-leaf' 32 'sampled' 'summary' 'sampled' 'sampled'
}
if ($IncludeStage3YRows) {
    function Add-Stage3YRow([string]$CaseId, [string]$DirectnessProbe, [string]$DirectnessCache, [string]$DirectnessVerify, [string]$ExpansionAudit = 'off', [string]$ExpansionAuditChannel = 'off', [string]$MomentAudit = 'off', [string]$ConvergenceAudit = 'off') {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split';
            MomentAudit = $MomentAudit; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = $ConvergenceAudit;
            ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid';
            ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only';
            AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0;
            MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = $ExpansionAuditChannel;
            LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2';
            TargetBatchSize = 4; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; LeafPairSelectionProfile = 'summary';
            LeafPairProbeCache = 'per-matvec'; LeafPairTargetGroup = 'same-target-leaf'; LeafPairMaxLookahead = 32;
            LeafPairDirectnessProbe = $DirectnessProbe; LeafPairDirectnessCache = $DirectnessCache; LeafPairDirectnessVerify = $DirectnessVerify;
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25;
            ExpansionAudit = $ExpansionAudit; TaylorOrder = 4; Theta = 0.60
        }
    }
    Add-Stage3YRow 'stage3y_perf2_full_traversal_baseline' 'full-traversal' 'off' 'off'
    Add-Stage3YRow 'stage3y_ancestor_path_only' 'ancestor-path' 'off' 'off'
    Add-Stage3YRow 'stage3y_ancestor_path_pair_cache' 'ancestor-path' 'target-leaf-source-leaf' 'off'
    Add-Stage3YRow 'stage3y_ancestor_path_pair_cache_verify_sampled' 'ancestor-path' 'target-leaf-source-leaf' 'sampled'
    Add-Stage3YRow 'stage3y_validation_ancestor_pair_cache_audit' 'ancestor-path' 'target-leaf-source-leaf' 'sampled' 'sampled' 'summary' 'sampled' 'sampled'
}
if ($IncludeStage3ZRows) {
    function Add-Stage3ZRow(
        [string]$CaseId,
        [string]$Scheduler,
        [int]$MaxActive,
        [int]$MinSourcePanels,
        [string]$WorklistVerify = 'off',
        [string]$ExpansionAudit = 'off',
        [string]$ExpansionAuditChannel = 'off',
        [string]$MomentAudit = 'off',
        [string]$ConvergenceAudit = 'off'
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split';
            MomentAudit = $MomentAudit; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = $ConvergenceAudit;
            ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid';
            ConvergenceCalibration = 'consecutive-aspect'; ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only';
            AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0;
            MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; ExpansionAuditChannel = $ExpansionAuditChannel;
            LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2';
            TargetBatchSize = 4; LeafPairMinSourcePanels = $MinSourcePanels; LeafPairMinTargets = 4; LeafPairSelectionProfile = 'summary';
            LeafPairProbeCache = 'per-matvec'; LeafPairTargetGroup = 'same-target-leaf'; LeafPairMaxLookahead = 32;
            LeafPairDirectnessProbe = 'ancestor-path'; LeafPairDirectnessCache = 'target-leaf-source-leaf'; LeafPairDirectnessVerify = 'off';
            LeafPairScheduler = $Scheduler; LeafPairWorklistMaxActive = $MaxActive; LeafPairWorklistVerify = $WorklistVerify;
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25;
            ExpansionAudit = $ExpansionAudit; TaylorOrder = 4; Theta = 0.60
        }
    }
    Add-Stage3ZRow 'stage3z_perf2_inline_select_baseline' 'inline-select' 2048 64
    Add-Stage3ZRow 'stage3z_deferred_worklist_max2048_min64' 'deferred-worklist' 2048 64
    Add-Stage3ZRow 'stage3z_deferred_worklist_max512_min64' 'deferred-worklist' 512 64
    Add-Stage3ZRow 'stage3z_deferred_worklist_max4096_min64' 'deferred-worklist' 4096 64
    Add-Stage3ZRow 'stage3z_deferred_worklist_max2048_min32' 'deferred-worklist' 2048 32
    Add-Stage3ZRow 'stage3z_validation_deferred_worklist_audit' 'deferred-worklist' 2048 64 'sampled' 'sampled' 'summary' 'sampled' 'sampled'
}
if ($IncludeStage3AARows) {
    function Add-Stage3AARow(
        [string]$CaseId,
        [string]$Scheduler,
        [int]$MinSourcePanels,
        [string]$LeafBlockVerify = 'off',
        [string]$ExpansionAudit = 'off',
        [string]$ExpansionAuditChannel = 'off',
        [string]$ConvergenceAudit = 'off'
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid';
            FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive';
            ConvergenceAudit = $ConvergenceAudit; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0;
            ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceCalibration = 'consecutive-aspect';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global';
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7;
            K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25;
            LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2';
            TargetBatchSize = 4; LeafPairMinSourcePanels = $MinSourcePanels; LeafPairMinTargets = 4;
            LeafPairSelectionProfile = 'summary'; LeafPairProbeCache = 'per-matvec'; LeafPairTargetGroup = 'same-target-leaf';
            LeafPairMaxLookahead = 32; LeafPairDirectnessProbe = 'ancestor-path'; LeafPairDirectnessCache = 'target-leaf-source-leaf';
            LeafPairDirectnessVerify = 'off'; LeafPairScheduler = $Scheduler; LeafPairWorklistMaxActive = 2048;
            LeafPairWorklistVerify = 'off'; LeafBlockVerify = $LeafBlockVerify;
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25;
            ExpansionAudit = $ExpansionAudit; ExpansionAuditChannel = $ExpansionAuditChannel; TaylorOrder = 4; Theta = 0.60
        }
    }
    $TreeRows += @{ CaseId = 'stage3aa_safe_baseline'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'off'; LeafPairKernel = 'auto'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ExpansionAudit = 'off'; TaylorOrder = 4; Theta = 0.60 }
    Add-Stage3AARow 'stage3aa_deferred_worklist_diagnostic_min64' 'deferred-worklist' 64
    Add-Stage3AARow 'stage3aa_target_leaf_block_min64' 'target-leaf-block' 64
    Add-Stage3AARow 'stage3aa_target_leaf_block_min32' 'target-leaf-block' 32
    Add-Stage3AARow 'stage3aa_target_leaf_block_min128' 'target-leaf-block' 128
    Add-Stage3AARow 'stage3aa_validation_target_leaf_block_audit' 'target-leaf-block' 64 'sampled' 'sampled' 'summary' 'sampled'
}
if ($IncludeStage3ABRows) {
    function Add-Stage3ABRow(
        [string]$CaseId,
        [string]$Scheduler,
        [int]$MinSourcePanels,
        [string]$DirectKernel,
        [string]$LeafBlockVerify = 'off',
        [string]$ExpansionAudit = 'off',
        [string]$ExpansionAuditChannel = 'off',
        [string]$ConvergenceAudit = 'off'
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid';
            FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive';
            ConvergenceAudit = $ConvergenceAudit; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0;
            ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceCalibration = 'consecutive-aspect';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global';
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7;
            K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25;
            LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2';
            LeafKernelProfile = 'summary'; JcpDirectKernel = $DirectKernel;
            TargetBatchSize = 4; LeafPairMinSourcePanels = $MinSourcePanels; LeafPairMinTargets = 4;
            LeafPairSelectionProfile = 'summary'; LeafPairProbeCache = 'per-matvec'; LeafPairTargetGroup = 'same-target-leaf';
            LeafPairMaxLookahead = 32; LeafPairDirectnessProbe = 'ancestor-path'; LeafPairDirectnessCache = 'target-leaf-source-leaf';
            LeafPairDirectnessVerify = 'off'; LeafPairScheduler = $Scheduler; LeafPairWorklistMaxActive = 2048;
            LeafPairWorklistVerify = 'off'; LeafBlockVerify = $LeafBlockVerify;
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25;
            ExpansionAudit = $ExpansionAudit; ExpansionAuditChannel = $ExpansionAuditChannel; TaylorOrder = 4; Theta = 0.60
        }
    }
    $TreeRows += @{ CaseId = 'stage3ab_safe_baseline'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'off'; LeafPairKernel = 'auto'; LeafKernelProfile = 'summary'; JcpDirectKernel = 'baseline'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ExpansionAudit = 'off'; TaylorOrder = 4; Theta = 0.60 }
    Add-Stage3ABRow 'stage3ab_target_leaf_block_min64_baseline' 'target-leaf-block' 64 'baseline'
    Add-Stage3ABRow 'stage3ab_target_leaf_block_min64_soa_avx2' 'target-leaf-block' 64 'soa-avx2'
    Add-Stage3ABRow 'stage3ab_target_leaf_block_min32_soa_avx2' 'target-leaf-block' 32 'soa-avx2'
    Add-Stage3ABRow 'stage3ab_target_leaf_block_min128_soa_avx2' 'target-leaf-block' 128 'soa-avx2'
    Add-Stage3ABRow 'stage3ab_validation_target_leaf_block_soa_audit' 'target-leaf-block' 64 'soa-avx2' 'sampled' 'sampled' 'summary' 'sampled'
}
if ($IncludeStage3AIRows) {
    function Add-Stage3AIRow(
        [string]$CaseId,
        [string]$Scheduler,
        [int]$MinSourcePanels,
        [string]$DirectKernel,
        [string]$ExpansionAudit = 'off',
        [string]$ExpansionAuditChannel = 'off',
        [string]$ConvergenceAudit = 'off'
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid';
            FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive';
            ConvergenceAudit = $ConvergenceAudit; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0;
            ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceCalibration = 'consecutive-aspect';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global';
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7;
            K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25;
            LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2';
            LeafKernelProfile = 'summary'; MatvecCostProfile = 'summary'; JcpDirectKernel = $DirectKernel;
            TargetBatchSize = 4; LeafPairMinSourcePanels = $MinSourcePanels; LeafPairMinTargets = 4;
            LeafPairSelectionProfile = 'summary'; LeafPairProbeCache = 'per-matvec'; LeafPairTargetGroup = 'same-target-leaf';
            LeafPairMaxLookahead = 32; LeafPairDirectnessProbe = 'ancestor-path'; LeafPairDirectnessCache = 'target-leaf-source-leaf';
            LeafPairDirectnessVerify = 'off'; LeafPairScheduler = $Scheduler; LeafPairWorklistMaxActive = 2048;
            LeafPairWorklistVerify = 'off'; LeafBlockVerify = 'off';
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25;
            ExpansionAudit = $ExpansionAudit; ExpansionAuditChannel = $ExpansionAuditChannel; TaylorOrder = 4; Theta = 0.60
        }
    }
    $TreeRows += @{ CaseId = 'stage3ai_stage3u_safe_baseline'; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = 'off'; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceCalibration = 'consecutive-aspect'; ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global'; ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25; LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'off'; LeafPairKernel = 'auto'; LeafKernelProfile = 'summary'; MatvecCostProfile = 'summary'; JcpDirectKernel = 'baseline'; TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; ExpansionAudit = 'off'; TaylorOrder = 4; Theta = 0.60 }
    Add-Stage3AIRow 'stage3ai_current_best_perf_min64_baseline' 'target-leaf-block' 64 'baseline'
    Add-Stage3AIRow 'stage3ai_soa_avx2_min64' 'target-leaf-block' 64 'soa-avx2'
    Add-Stage3AIRow 'stage3ai_fused_soa_avx2_min64' 'target-leaf-block' 64 'fused-soa-avx2'
    Add-Stage3AIRow 'stage3ai_fused_soa_avx2_min32' 'target-leaf-block' 32 'fused-soa-avx2'
    Add-Stage3AIRow 'stage3ai_fused_soa_avx2_min128' 'target-leaf-block' 128 'fused-soa-avx2'
    Add-Stage3AIRow 'stage3ai_validation_fused_soa_avx2_audit' 'target-leaf-block' 64 'fused-soa-avx2' 'sampled' 'summary' 'sampled'
}
if ($IncludeStage3AJRows) {
    function Add-Stage3AJRow(
        [string]$CaseId,
        [string]$Preconditioner,
        [int]$BlockPanelCap,
        [int]$MemoryBudgetMb,
        [int]$Restart,
        [string]$ExpansionAudit = 'off',
        [string]$ExpansionAuditChannel = 'off',
        [string]$ConvergenceAudit = 'off'
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid';
            FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive';
            ConvergenceAudit = $ConvergenceAudit; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0;
            ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceCalibration = 'consecutive-aspect';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global';
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7;
            K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25;
            LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2';
            LeafKernelProfile = 'summary'; MatvecCostProfile = 'summary'; JcpDirectKernel = 'fused-soa-avx2';
            TargetBatchSize = 4; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4;
            LeafPairSelectionProfile = 'summary'; LeafPairProbeCache = 'per-matvec'; LeafPairTargetGroup = 'same-target-leaf';
            LeafPairMaxLookahead = 32; LeafPairDirectnessProbe = 'ancestor-path'; LeafPairDirectnessCache = 'target-leaf-source-leaf';
            LeafPairDirectnessVerify = 'off'; LeafPairScheduler = 'target-leaf-block'; LeafPairWorklistMaxActive = 2048;
            LeafPairWorklistVerify = 'off'; LeafBlockVerify = 'off';
            Preconditioner = $Preconditioner; PreconditionerBlockNear = 'buffer-aware';
            PreconditionerBlockPanelCap = $BlockPanelCap; PreconditionerMemoryBudgetMb = $MemoryBudgetMb; GmresRestart = $Restart;
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25;
            ExpansionAudit = $ExpansionAudit; ExpansionAuditChannel = $ExpansionAuditChannel; TaylorOrder = 4; Theta = 0.60
        }
    }
    Add-Stage3AJRow 'stage3aj_fused_diagonal_r64_control' 'diagonal' 32 256 64
    Add-Stage3AJRow 'stage3aj_budgeted_buffer_cap32_r64' 'block-budgeted' 32 256 64
    Add-Stage3AJRow 'stage3aj_budgeted_buffer_cap16_r64' 'block-budgeted' 16 256 64
    Add-Stage3AJRow 'stage3aj_budgeted_buffer_cap32_r32' 'block-budgeted' 32 256 32
    Add-Stage3AJRow 'stage3aj_budgeted_buffer_cap32_r100' 'block-budgeted' 32 256 100
    Add-Stage3AJRow 'stage3aj_validation_budgeted_cap32_r64_audit' 'block-budgeted' 32 256 64 'sampled' 'summary' 'sampled'
}
if ($IncludeStage3AKRows) {
    function Add-Stage3AKRow(
        [string]$CaseId,
        [string]$Preconditioner,
        [int]$CoarseAggregateCap,
        [string]$CoarseBase,
        [int]$TrueResidualProbeInterval,
        [string]$ExpansionAudit = 'off',
        [string]$ExpansionAuditChannel = 'off',
        [string]$ConvergenceAudit = 'off'
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid';
            FragmentationSafety = 'split'; MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive';
            ConvergenceAudit = $ConvergenceAudit; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0;
            ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceCalibration = 'consecutive-aspect';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global';
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7;
            K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25;
            LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'leaf-pair-balanced'; LeafPairKernel = 'target-avx2';
            LeafKernelProfile = 'summary'; MatvecCostProfile = 'summary'; JcpDirectKernel = 'fused-soa-avx2';
            TargetBatchSize = 4; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4;
            LeafPairSelectionProfile = 'summary'; LeafPairProbeCache = 'per-matvec'; LeafPairTargetGroup = 'same-target-leaf';
            LeafPairMaxLookahead = 32; LeafPairDirectnessProbe = 'ancestor-path'; LeafPairDirectnessCache = 'target-leaf-source-leaf';
            LeafPairDirectnessVerify = 'off'; LeafPairScheduler = 'target-leaf-block'; LeafPairWorklistMaxActive = 2048;
            LeafPairWorklistVerify = 'off'; LeafBlockVerify = 'off';
            Preconditioner = $Preconditioner; PreconditionerBlockNear = 'buffer-aware';
            PreconditionerBlockPanelCap = 32; PreconditionerMemoryBudgetMb = 256; GmresRestart = 64;
            PreconditionerCoarseAggregateCap = $CoarseAggregateCap; PreconditionerCoarseBase = $CoarseBase;
            PreconditionerCoarseVerify = 'off'; TrueResidualProbeInterval = $TrueResidualProbeInterval;
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25;
            ExpansionAudit = $ExpansionAudit; ExpansionAuditChannel = $ExpansionAuditChannel; TaylorOrder = 4; Theta = 0.60
        }
    }
    Add-Stage3AKRow 'stage3ak_fused_diagonal_interval1_control' 'diagonal' 16 'diagonal' 1
    Add-Stage3AKRow 'stage3ak_budgeted_interval1_control' 'block-budgeted' 16 'diagonal' 1
    Add-Stage3AKRow 'stage3ak_coarse_diag_cap8_interval1' 'coarse-galerkin-probe' 8 'diagonal' 1
    Add-Stage3AKRow 'stage3ak_coarse_diag_cap16_interval1' 'coarse-galerkin-probe' 16 'diagonal' 1
    Add-Stage3AKRow 'stage3ak_coarse_diag_cap32_interval1' 'coarse-galerkin-probe' 32 'diagonal' 1
    Add-Stage3AKRow 'stage3ak_coarse_block_cap16_interval1' 'coarse-galerkin-probe' 16 'block-budgeted' 1
    Add-Stage3AKRow 'stage3ak_validation_coarse_diag_cap16_audit' 'coarse-galerkin-probe' 16 'diagonal' 4 'sampled' 'summary' 'sampled'
}
if ($IncludeStage3ACRows) {
    function Add-Stage3ACRow(
        [string]$CaseId,
        [int]$LeafSize,
        [string]$MarginPolicy,
        [double]$AspectThreshold,
        [double]$RhoTrigger,
        [string]$ClusterPolicy = 'off',
        [string]$Scheduler = 'inline-select',
        [string]$DirectKernel = 'baseline',
        [string]$LeafPairKernel = 'auto',
        [string]$ExpansionAudit = 'off',
        [string]$ExpansionAuditChannel = 'off',
        [string]$MomentAudit = 'off',
        [string]$ConvergenceAudit = 'off'
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid';
            FragmentationSafety = 'split'; MomentAudit = $MomentAudit; MomentConstruction = 'translated'; MomentSum = 'naive';
            ConvergenceAudit = $ConvergenceAudit; ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0;
            ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps'; ConvergenceCalibration = 'consecutive-aspect';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global';
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7;
            K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25;
            MacSurfaceMarginPolicy = $MarginPolicy; MacSurfaceMarginAspectThreshold = $AspectThreshold; MacSurfaceMarginRhoTrigger = $RhoTrigger;
            LeafWorkProfile = 'summary'; LeafPanelKernel = 'tile-avx2'; ClusterPolicy = $ClusterPolicy; LeafPairKernel = $LeafPairKernel;
            LeafKernelProfile = 'summary'; JcpDirectKernel = $DirectKernel; TargetBatchSize = $(if ($ClusterPolicy -eq 'off') { 1 } else { 4 });
            LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; LeafPairSelectionProfile = 'summary';
            LeafPairProbeCache = 'per-matvec'; LeafPairTargetGroup = 'same-target-leaf'; LeafPairMaxLookahead = 32;
            LeafPairDirectnessProbe = 'ancestor-path'; LeafPairDirectnessCache = 'target-leaf-source-leaf'; LeafPairDirectnessVerify = 'off';
            LeafPairScheduler = $Scheduler; LeafPairWorklistMaxActive = 2048; LeafPairWorklistVerify = 'off'; LeafBlockVerify = 'off';
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25;
            ExpansionAudit = $ExpansionAudit; ExpansionAuditChannel = $ExpansionAuditChannel; TaylorOrder = 4; Theta = 0.60; LeafSize = $LeafSize
        }
    }
    Add-Stage3ACRow 'stage3ac_safe_baseline_leaf128_margin_all' 128 'all' 3.0 0.55
    Add-Stage3ACRow 'stage3ac_leaf64_margin_all' 64 'all' 3.0 0.55
    Add-Stage3ACRow 'stage3ac_leaf96_margin_all' 96 'all' 3.0 0.55
    Add-Stage3ACRow 'stage3ac_leaf128_margin_all' 128 'all' 3.0 0.55
    Add-Stage3ACRow 'stage3ac_leaf192_margin_all' 192 'all' 3.0 0.55
    Add-Stage3ACRow 'stage3ac_adaptive_aspect3p0_rho0p55' 128 'adaptive' 3.0 0.55
    Add-Stage3ACRow 'stage3ac_adaptive_aspect2p5_rho0p55' 128 'adaptive' 2.5 0.55
    Add-Stage3ACRow 'stage3ac_adaptive_aspect3p0_rho0p50' 128 'adaptive' 3.0 0.50
    Add-Stage3ACRow 'stage3ac_validation_adaptive_aspect3p0_rho0p55_audit' 128 'adaptive' 3.0 0.55 'off' 'inline-select' 'baseline' 'auto' 'sampled' 'summary' 'sampled' 'sampled'
    Add-Stage3ACRow 'stage3ac_timing_adaptive_aspect3p0_rho0p55_block_soa' 128 'adaptive' 3.0 0.55 'leaf-pair-balanced' 'target-leaf-block' 'soa-avx2' 'target-avx2'
}
if ($IncludeStage3ADRows) {
    function Add-Stage3ADRow(
        [string]$CaseId,
        [int]$LeafSize,
        [string]$MarginPolicy,
        [double]$HardAspectThreshold,
        [double]$HardRhoTrigger,
        [double]$HardFactor = 0.25,
        [double]$SoftAspectThreshold = 0.0,
        [double]$SoftRhoTrigger = 0.0,
        [double]$SoftFactor = 0.10,
        [string]$ClusterPolicy = 'off',
        [string]$Scheduler = 'inline-select',
        [string]$DirectKernel = 'baseline',
        [string]$LeafPairKernel = 'auto',
        [string]$ExpansionAudit = 'sampled',
        [string]$ExpansionAuditChannel = 'summary',
        [string]$MomentAudit = 'sampled',
        [string]$ConvergenceAudit = 'sampled'
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap'; ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split';
            MomentAudit = $MomentAudit; MomentConstruction = 'translated'; MomentSum = 'naive'; ConvergenceAudit = $ConvergenceAudit;
            ConvergenceSafety = 'split'; ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid';
            ConvergenceCalibration = 'consecutive-aspect'; ConvergenceEstimatorModel = 'channel-tail-eps';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global';
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary'; AcceptedMissAbsoluteFloor = 1.0e-7;
            K3AcceptedRhoCap = 0.0; MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = $HardFactor;
            MacSurfaceMarginPolicy = $MarginPolicy; MacSurfaceMarginAspectThreshold = $HardAspectThreshold; MacSurfaceMarginRhoTrigger = $HardRhoTrigger;
            MacSurfaceMarginSoftAspectThreshold = $SoftAspectThreshold; MacSurfaceMarginSoftRhoTrigger = $SoftRhoTrigger; MacSurfaceMarginSoftFactor = $SoftFactor;
            LeafWorkProfile = 'summary'; LeafPanelKernel = 'tile-avx2'; ClusterPolicy = $ClusterPolicy; LeafPairKernel = $LeafPairKernel;
            LeafKernelProfile = 'summary'; JcpDirectKernel = $DirectKernel; TargetBatchSize = $(if ($ClusterPolicy -eq 'off') { 1 } else { 4 });
            LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4; LeafPairSelectionProfile = 'summary';
            LeafPairProbeCache = 'per-matvec'; LeafPairTargetGroup = 'same-target-leaf'; LeafPairMaxLookahead = 32;
            LeafPairDirectnessProbe = 'ancestor-path'; LeafPairDirectnessCache = 'target-leaf-source-leaf'; LeafPairDirectnessVerify = 'off';
            LeafPairScheduler = $Scheduler; LeafPairWorklistMaxActive = 2048; LeafPairWorklistVerify = 'off'; LeafBlockVerify = 'off';
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3; ConvergenceUnstableRhoCap = 0.25;
            ExpansionAudit = $ExpansionAudit; ExpansionAuditChannel = $ExpansionAuditChannel; TaylorOrder = 4; Theta = 0.60; LeafSize = $LeafSize
        }
    }
    Add-Stage3ADRow 'stage3ad_safe_baseline_leaf128_margin_all' 128 'all' 3.0 0.55 0.25 0.0 0.0 0.10
    Add-Stage3ADRow 'stage3ad_validation_adaptive_aspect3p0_rho0p55_audit' 128 'adaptive' 3.0 0.55 0.25 0.0 0.0 0.10
    Add-Stage3ADRow 'stage3ad_validation_adaptive_aspect2p5_rho0p55_audit' 128 'adaptive' 2.5 0.55 0.25 0.0 0.0 0.10
    Add-Stage3ADRow 'stage3ad_validation_adaptive_aspect3p0_rho0p50_audit' 128 'adaptive' 3.0 0.50 0.25 0.0 0.0 0.10
    Add-Stage3ADRow 'stage3ad_adaptive_aspect3p0_rho0p55_factor0p20' 128 'adaptive' 3.0 0.55 0.20 0.0 0.0 0.10
    Add-Stage3ADRow 'stage3ad_adaptive_aspect3p0_rho0p55_factor0p15' 128 'adaptive' 3.0 0.55 0.15 0.0 0.0 0.10
    Add-Stage3ADRow 'stage3ad_risk_scaled_hard3p0_rho0p55_soft2p0_rho0p45' 128 'risk-scaled' 3.0 0.55 0.25 2.0 0.45 0.10
    Add-Stage3ADRow 'stage3ad_risk_scaled_hard3p0_rho0p55_soft2p0_rho0p45_factor0p20' 128 'risk-scaled' 3.0 0.55 0.20 2.0 0.45 0.10
    Add-Stage3ADRow 'stage3ad_risk_scaled_hard2p5_rho0p55_soft2p0_rho0p45_factor0p20' 128 'risk-scaled' 2.5 0.55 0.20 2.0 0.45 0.10
    Add-Stage3ADRow 'stage3ad_risk_scaled_leaf96' 96 'risk-scaled' 3.0 0.55 0.20 2.0 0.45 0.10
    Add-Stage3ADRow 'stage3ad_risk_scaled_leaf192' 192 'risk-scaled' 3.0 0.55 0.20 2.0 0.45 0.10
    Add-Stage3ADRow 'stage3ad_timing_risk_scaled_block_soa' 128 'risk-scaled' 3.0 0.55 0.20 2.0 0.45 0.10 'leaf-pair-balanced' 'target-leaf-block' 'soa-avx2' 'target-avx2' 'off' 'off' 'off' 'off'
}
if ($IncludeStage3AERows) {
    function Add-Stage3AERow(
        [string]$CaseId,
        [string]$SplitSelectivity,
        [double]$SplitAbsoluteFloor,
        [string]$OverrejectAudit,
        [string]$ExpansionAudit = 'sampled',
        [string]$ExpansionAuditChannel = 'summary',
        [string]$ConvergenceAudit = 'sampled',
        [string]$ClusterPolicy = 'off',
        [string]$LeafPairScheduler = 'inline-select',
        [string]$JcpDirectKernel = 'baseline',
        [string]$LeafPairKernel = 'auto',
        [int]$TargetBatchSize = 1
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap';
            ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split';
            MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive';
            ConvergenceAudit = $ConvergenceAudit; ConvergenceSafety = 'split';
            ConvergenceRhoCap = 0.0; MacAspectCap = 0.0;
            ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps';
            ConvergenceCalibration = 'consecutive-aspect'; ConvergenceChannelCalibration = 'polynomial-eps';
            ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global';
            ConvergenceSplitSelectivity = $SplitSelectivity;
            ConvergenceSplitAbsoluteFloor = $SplitAbsoluteFloor;
            ConvergenceOverrejectAudit = $OverrejectAudit;
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary';
            AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0;
            MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25;
            MacSurfaceMarginPolicy = 'all'; MacSurfaceMarginAspectThreshold = 3.0;
            MacSurfaceMarginRhoTrigger = 0.55; MacSurfaceMarginSoftAspectThreshold = 0.0;
            MacSurfaceMarginSoftRhoTrigger = 0.0; MacSurfaceMarginSoftFactor = 0.10;
            LeafWorkProfile = 'summary'; LeafPanelKernel = 'tile-avx2';
            ClusterPolicy = $ClusterPolicy; LeafPairKernel = $LeafPairKernel;
            JcpDirectKernel = $JcpDirectKernel; TargetBatchSize = $TargetBatchSize;
            LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4;
            LeafPairScheduler = $LeafPairScheduler;
            ExpansionAudit = $ExpansionAudit; ExpansionAuditChannel = $ExpansionAuditChannel;
            TaylorOrder = 4; Theta = 0.60; LeafSize = 128
        }
    }
    Add-Stage3AERow 'stage3ae_safe_baseline_legacy_split' 'legacy' 0.0 'off' 'sampled' 'summary' 'sampled'
    Add-Stage3AERow 'stage3ae_overreject_audit_legacy_split' 'legacy' 0.0 'sampled' 'off' 'off' 'off'
    Add-Stage3AERow 'stage3ae_absolute_aware_floor1e8' 'absolute-aware' 1.0e-8 'sampled'
    Add-Stage3AERow 'stage3ae_absolute_aware_floor3e8' 'absolute-aware' 3.0e-8 'sampled'
    Add-Stage3AERow 'stage3ae_absolute_aware_floor1e7' 'absolute-aware' 1.0e-7 'sampled'
    Add-Stage3AERow 'stage3ae_channel_absolute_aware_floor1e8' 'channel-absolute-aware' 1.0e-8 'sampled'
    Add-Stage3AERow 'stage3ae_channel_absolute_aware_floor3e8' 'channel-absolute-aware' 3.0e-8 'sampled'
    Add-Stage3AERow 'stage3ae_channel_absolute_aware_floor1e7' 'channel-absolute-aware' 1.0e-7 'sampled'
    Add-Stage3AERow 'stage3ae_timing_best_absaware_block_soa' 'channel-absolute-aware' 3.0e-8 'off' 'off' 'off' 'off' 'leaf-pair-balanced' 'target-leaf-block' 'soa-avx2' 'target-avx2' 4
}
if ($IncludeStage3AGRows) {
    function Add-Stage3AGRow(
        [string]$CaseId,
        [string]$OverrejectAudit,
        [string]$OverrejectDetail
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap';
            ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split';
            MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive';
            ConvergenceAudit = 'off'; ConvergenceSafety = 'split';
            ConvergenceRhoCap = 0.0; MacAspectCap = 0.0;
            ConvergenceEstimator = 'hybrid'; ConvergenceEstimatorModel = 'channel-tail-eps';
            ConvergenceCalibration = 'consecutive-aspect'; ConvergenceChannelCalibration = 'polynomial-eps';
            ConvergenceDielectricAmp = 'on'; ConvergenceSplitPolicy = 'global';
            ConvergenceSplitSelectivity = 'legacy'; ConvergenceSplitAbsoluteFloor = 0.0;
            ConvergenceOverrejectAudit = $OverrejectAudit; ConvergenceOverrejectDetail = $OverrejectDetail;
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary';
            AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0;
            MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25;
            MacSurfaceMarginPolicy = 'all'; MacSurfaceMarginAspectThreshold = 3.0;
            MacSurfaceMarginRhoTrigger = 0.55; MacSurfaceMarginSoftAspectThreshold = 0.0;
            MacSurfaceMarginSoftRhoTrigger = 0.0; MacSurfaceMarginSoftFactor = 0.10;
            LeafWorkProfile = 'summary'; LeafPanelKernel = 'tile-avx2';
            ClusterPolicy = 'off'; LeafPairKernel = 'auto'; JcpDirectKernel = 'baseline';
            TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4;
            LeafPairScheduler = 'inline-select'; ExpansionAudit = 'off'; ExpansionAuditChannel = 'off';
            TaylorOrder = 4; Theta = 0.60; LeafSize = 128
        }
    }
    Add-Stage3AGRow 'stage3ag_safe_baseline_audit_off' 'off' 'off'
    Add-Stage3AGRow 'stage3ag_overreject_calibration_summary' 'sampled' 'summary'
    Add-Stage3AGRow 'stage3ag_overreject_calibration_detail' 'sampled' 'sampled-detail'
}
if ($IncludeStage3AHRows) {
    function Add-Stage3AHRow(
        [string]$CaseId,
        [string]$ShapeAudit,
        [string]$OverrejectDetail
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap';
            ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split';
            MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive';
            ConvergenceAudit = 'off'; ExpansionAudit = 'off'; ExpansionAuditChannel = 'off';
            TaylorOrder = 4; Theta = 0.60; ConvergenceSafety = 'split';
            ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid';
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3;
            ConvergenceUnstableRhoCap = 0.25; ConvergenceCalibration = 'consecutive-aspect';
            ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global';
            ConvergenceSplitSelectivity = 'legacy'; ConvergenceSplitAbsoluteFloor = 0.0;
            ConvergenceOverrejectAudit = 'sampled'; ConvergenceOverrejectDetail = $OverrejectDetail;
            ConvergenceShapeAudit = $ShapeAudit;
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary';
            AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0;
            MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25;
            MacSurfaceMarginPolicy = 'all'; MacSurfaceMarginAspectThreshold = 3.0;
            MacSurfaceMarginRhoTrigger = 0.55; MacSurfaceMarginSoftAspectThreshold = 0.0;
            MacSurfaceMarginSoftRhoTrigger = 0.0; MacSurfaceMarginSoftFactor = 0.10;
            LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'off'; LeafPairKernel = 'auto';
            TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4;
            LeafWorkProfile = 'summary'; JcpDirectKernel = 'baseline';
        }
    }
    Add-Stage3AHRow 'stage3ah_stage3ag_detail_baseline' 'off' 'sampled-detail'
    Add-Stage3AHRow 'stage3ah_shape_audit_summary' 'summary' 'off'
    Add-Stage3AHRow 'stage3ah_shape_audit_detail' 'sampled-detail' 'off'
}
if ($IncludeStage3ALRows) {
    function Add-Stage3ALRow(
        [string]$CaseId,
        [string]$KernelResponseAudit
    ) {
        $script:TreeRows += @{
            CaseId = $CaseId; TreeBuilder = 'topdown-aabb'; MacSafety = 'physical-gap';
            ExpansionCenter = 'area-centroid'; FragmentationSafety = 'split';
            MomentAudit = 'off'; MomentConstruction = 'translated'; MomentSum = 'naive';
            ConvergenceAudit = 'off'; ExpansionAudit = 'off'; ExpansionAuditChannel = 'off';
            TaylorOrder = 4; Theta = 0.60; ConvergenceSafety = 'split';
            ConvergenceRhoCap = 0.0; MacAspectCap = 0.0; ConvergenceEstimator = 'hybrid';
            ConvergenceTailThreshold = 1.0e4; ConvergenceEstimatorErrorCap = 1.0e-3;
            ConvergenceUnstableRhoCap = 0.25; ConvergenceCalibration = 'consecutive-aspect';
            ConvergenceDielectricAmp = 'on'; ConvergenceEstimatorModel = 'channel-tail-eps';
            ConvergenceChannelCalibration = 'polynomial-eps'; ConvergenceSplitPolicy = 'global';
            ConvergenceSplitSelectivity = 'legacy'; ConvergenceSplitAbsoluteFloor = 0.0;
            ConvergenceOverrejectAudit = if ($KernelResponseAudit -eq 'off') { 'off' } else { 'sampled' };
            ConvergenceOverrejectDetail = 'off'; ConvergenceShapeAudit = 'off';
            ConvergenceKernelResponseAudit = $KernelResponseAudit;
            ChannelMissedSemantics = 'accepted-only'; AcceptedMissDiagnostics = 'summary';
            AcceptedMissAbsoluteFloor = 1.0e-7; K3AcceptedRhoCap = 0.0;
            MacDistancePolicy = 'surface-margin'; MacSurfaceMarginFactor = 0.25;
            MacSurfaceMarginPolicy = 'all'; MacSurfaceMarginAspectThreshold = 3.0;
            MacSurfaceMarginRhoTrigger = 0.55; MacSurfaceMarginSoftAspectThreshold = 0.0;
            MacSurfaceMarginSoftRhoTrigger = 0.0; MacSurfaceMarginSoftFactor = 0.10;
            LeafPanelKernel = 'tile-avx2'; ClusterPolicy = 'off'; LeafPairKernel = 'auto';
            TargetBatchSize = 1; LeafPairMinSourcePanels = 64; LeafPairMinTargets = 4;
            LeafWorkProfile = 'summary'; JcpDirectKernel = 'baseline';
        }
    }
    Add-Stage3ALRow 'stage3al_safe_baseline_audit_off' 'off'
    Add-Stage3ALRow 'stage3al_kernel_response_summary' 'summary'
    Add-Stage3ALRow 'stage3al_kernel_response_detail' 'sampled-detail'
}

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDoubleOrNull([object]$Value) {
    if ($Value -is [double] -or $Value -is [single] -or $Value -is [decimal] -or
        $Value -is [int] -or $Value -is [long]) {
        return [double]$Value
    }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $null }
    return [double]::Parse($text, $InvariantCulture)
}

function Format-Number([object]$Value) {
    if ($null -eq $Value) { return 'NA' }
    if ($Value -is [string]) { return $Value }
    return ([double]$Value).ToString('G17', $InvariantCulture)
}

function Get-DerivedProduct($Left, $Right) {
    $leftValue = Convert-ToDoubleOrNull $Left
    $rightValue = Convert-ToDoubleOrNull $Right
    if ($null -eq $leftValue -or $null -eq $rightValue) { return $null }
    return $leftValue * $rightValue
}

function Get-DerivedRatio($Numerator, $Denominator) {
    $numValue = Convert-ToDoubleOrNull $Numerator
    $denValue = Convert-ToDoubleOrNull $Denominator
    if ($null -eq $numValue -or $null -eq $denValue -or [Math]::Abs($denValue) -le 1.0e-300) {
        return $null
    }
    return $numValue / $denValue
}

function Get-RowValue([object]$Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI Stage 3F root: $full"
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
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

function Resolve-PqrFixture([string]$Fixture) {
    if (-not $FixtureMap.ContainsKey($Fixture)) {
        throw "Unknown Stage 3F fixture: $Fixture"
    }
    $path = Join-Path $InputRoot $FixtureMap[$Fixture]
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Missing PQR fixture for ${Fixture}: $path"
    }
    return $path
}

function New-ArgFile(
    [string]$RunDir,
    [string]$Fixture,
    [string]$Mode,
    [string]$TreeBuilder,
    [string]$MacSafety,
    [string]$ExpansionCenter = 'legacy',
    [string]$FragmentationSafety = 'off',
    [string]$MomentAudit = 'off',
    [string]$MomentConstruction = 'translated',
    [string]$MomentSum = 'naive',
    [string]$ConvergenceAudit = 'off',
    [string]$ExpansionAudit = 'off',
    [int]$TaylorOrder = 3,
    [double]$Theta = $OpeningTheta,
    [string]$ConvergenceSafety = 'off',
    [double]$ConvergenceRhoCap = 0.0,
    [double]$MacAspectCap = 0.0,
    [string]$ConvergenceEstimator = 'off',
    [double]$ConvergenceTailThreshold = 1.0e4,
    [double]$ConvergenceEstimatorErrorCap = 1.0e-3,
    [double]$ConvergenceUnstableRhoCap = 0.25,
    [string]$ConvergenceCalibration = 'stage3j',
    [string]$ConvergenceDielectricAmp = 'off',
    [string]$ConvergenceEstimatorModel = 'stage3l',
    [string]$ConvergenceChannelCalibration = 'stage3m',
    [string]$ConvergenceSplitPolicy = 'global',
    [string]$ConvergenceSplitSelectivity = 'legacy',
    [double]$ConvergenceSplitAbsoluteFloor = 0.0,
    [string]$ConvergenceOverrejectAudit = 'off',
    [string]$ConvergenceOverrejectDetail = 'off',
    [string]$ConvergenceShapeAudit = 'off',
    [string]$ConvergenceKernelResponseAudit = 'off',
    [string]$ChannelMissedSemantics = 'legacy',
    [string]$AcceptedMissDiagnostics = 'off',
    [double]$AcceptedMissAbsoluteFloor = 1.0e-7,
    [double]$K3AcceptedRhoCap = 0.0,
    [string]$MacDistancePolicy = 'center',
    [double]$MacSurfaceMarginFactor = 0.25,
    [string]$MacSurfaceMarginPolicy = 'all',
    [double]$MacSurfaceMarginAspectThreshold = 3.0,
    [double]$MacSurfaceMarginRhoTrigger = 0.55,
    [double]$MacSurfaceMarginSoftAspectThreshold = 0.0,
    [double]$MacSurfaceMarginSoftRhoTrigger = 0.0,
    [double]$MacSurfaceMarginSoftFactor = 0.10,
    [string]$ExpansionAuditChannel = 'off',
    [string]$LeafPanelKernel = 'auto',
    [string]$ClusterPolicy = 'off',
    [string]$LeafPairKernel = 'auto',
    [string]$LeafKernelProfile = 'off',
    [string]$LeafWorkProfile = 'off',
    [string]$MatvecCostProfile = 'off',
    [string]$JcpDirectKernel = 'baseline',
    [int]$TargetBatchSize = 1,
    [int]$LeafPairMinSourcePanels = 64,
    [int]$LeafPairMinTargets = 4,
    [string]$LeafPairSelectionProfile = 'off',
    [string]$LeafPairProbeCache = 'off',
    [string]$LeafPairTargetGroup = 'scan',
    [int]$LeafPairMaxLookahead = 32,
    [string]$LeafPairDirectnessProbe = 'full-traversal',
    [string]$LeafPairDirectnessCache = 'off',
    [string]$LeafPairDirectnessVerify = 'off',
    [string]$LeafPairScheduler = 'inline-select',
    [int]$LeafPairWorklistMaxActive = 2048,
    [string]$LeafPairWorklistVerify = 'off',
    [string]$LeafBlockVerify = 'off',
    [int]$TreeLeafSizeOverride = $TreeLeafSize,
    [string]$Preconditioner = 'diagonal',
    [string]$PreconditionerBlockNear = 'singular',
    [int]$PreconditionerBlockPanelCap = 32,
    [int]$PreconditionerMemoryBudgetMb = 256,
    [int]$PreconditionerCoarseAggregateCap = 16,
    [string]$PreconditionerCoarseBase = 'diagonal',
    [string]$PreconditionerCoarseVerify = 'off',
    [int]$GmresRestart = 32,
    [int]$TrueResidualProbeInterval = 4
) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $pqrPath = Resolve-PqrFixture $Fixture
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $thetaText = $Theta.ToString('G17', $InvariantCulture)
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
        '-pbe-probe=1.4',
        '-pbe-scale=0.35',
        '-pbe-margin=2.0',
        '-pbe-max-grid=300000',
        '-pbe-map=off',
        '-pbe-max-iter=800',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-timing=detail',
        '-pbe-salt=0.15',
        '-pbe-temp=298.15',
        "-nt=$Threads",
        '-pbe-tabi-derivatives=analytic',
        '-pbe-tabi-energy-projection=boundary',
        "-pbe-tabi-precond=$Preconditioner",
        "-pbe-tabi-precond-block-near=$PreconditionerBlockNear",
        "-pbe-tabi-precond-block-panel-cap=$PreconditionerBlockPanelCap",
        "-pbe-tabi-precond-memory-budget-mb=$PreconditionerMemoryBudgetMb",
        "-pbe-tabi-precond-coarse-aggregate-cap=$PreconditionerCoarseAggregateCap",
        "-pbe-tabi-precond-coarse-base=$PreconditionerCoarseBase",
        "-pbe-tabi-precond-coarse-verify=$PreconditionerCoarseVerify",
        "-pbe-tabi-surface-points=$SurfacePoints",
        "-pbe-tabi-taylor-order=$TaylorOrder",
        "-pbe-tabi-tree-leaf-size=$TreeLeafSizeOverride",
        "-pbe-tabi-opening-theta=$thetaText",
        '-pbe-tabi-node-bounds=tight-aabb',
        '-pbe-tabi-interaction-mode=onthefly',
        '-pbe-tabi-near-integration=duffy_adaptive',
        '-pbe-tabi-quad-buffer=adaptive',
        "-pbe-tabi-gmres-restart=$GmresRestart",
        '-pbe-tabi-true-residual-probe=interval',
        '-pbe-tabi-true-residual-probe-scope=adaptive',
        "-pbe-tabi-true-residual-probe-interval=$TrueResidualProbeInterval",
        '-pbe-tabi-true-residual-sample-count=1024',
        '-pbe-tabi-progress=auto',
        '-pbe-tabi-progress-interval=5',
        '-pbe-tabi-matvec-audit-sample-count=' + $MatvecAuditSampleCount
    )
    if ($Mode -eq 'direct') {
        $argLines += '-pbe-tabi-sum=direct'
        $argLines += '-pbe-tabi-formulation=jcp2013'
        $argLines += '-pbe-tabi-matvec-audit=off'
    } else {
        $argLines += '-pbe-tabi-sum=treecode'
        $argLines += '-pbe-tabi-formulation=jcp2013-treecode'
        $argLines += "-pbe-tabi-tree-builder=$TreeBuilder"
        $argLines += "-pbe-tabi-mac-safety=$MacSafety"
        $argLines += "-pbe-tabi-expansion-center=$ExpansionCenter"
        $argLines += "-pbe-tabi-fragmentation-safety=$FragmentationSafety"
        $argLines += "-pbe-tabi-moment-audit=$MomentAudit"
        $argLines += "-pbe-tabi-moment-construction=$MomentConstruction"
        $argLines += "-pbe-tabi-moment-sum=$MomentSum"
        $argLines += "-pbe-tabi-convergence-audit=$ConvergenceAudit"
        $argLines += "-pbe-tabi-convergence-safety=$ConvergenceSafety"
        $argLines += "-pbe-tabi-convergence-rho-cap=$($ConvergenceRhoCap.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-mac-aspect-cap=$($MacAspectCap.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-convergence-estimator=$ConvergenceEstimator"
        $argLines += "-pbe-tabi-convergence-tail-threshold=$($ConvergenceTailThreshold.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-convergence-estimator-error-cap=$($ConvergenceEstimatorErrorCap.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-convergence-unstable-rho-cap=$($ConvergenceUnstableRhoCap.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-convergence-calibration=$ConvergenceCalibration"
        $argLines += "-pbe-tabi-convergence-dielectric-amp=$ConvergenceDielectricAmp"
        $argLines += "-pbe-tabi-convergence-estimator-model=$ConvergenceEstimatorModel"
        $argLines += "-pbe-tabi-convergence-channel-calibration=$ConvergenceChannelCalibration"
        $argLines += "-pbe-tabi-convergence-split-policy=$ConvergenceSplitPolicy"
        $argLines += "-pbe-tabi-convergence-split-selectivity=$ConvergenceSplitSelectivity"
        $argLines += "-pbe-tabi-convergence-split-absolute-floor=$($ConvergenceSplitAbsoluteFloor.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-convergence-overreject-audit=$ConvergenceOverrejectAudit"
        $argLines += "-pbe-tabi-convergence-overreject-detail=$ConvergenceOverrejectDetail"
        $argLines += "-pbe-tabi-convergence-shape-audit=$ConvergenceShapeAudit"
        $argLines += "-pbe-tabi-convergence-kernel-response-audit=$ConvergenceKernelResponseAudit"
        $argLines += "-pbe-tabi-channel-missed-semantics=$ChannelMissedSemantics"
        $argLines += "-pbe-tabi-accepted-miss-diagnostics=$AcceptedMissDiagnostics"
        $argLines += "-pbe-tabi-accepted-miss-absolute-floor=$($AcceptedMissAbsoluteFloor.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-convergence-k3-accepted-rho-cap=$($K3AcceptedRhoCap.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-mac-distance-policy=$MacDistancePolicy"
        $argLines += "-pbe-tabi-mac-surface-margin-factor=$($MacSurfaceMarginFactor.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-mac-surface-margin-policy=$MacSurfaceMarginPolicy"
        $argLines += "-pbe-tabi-mac-surface-margin-aspect-threshold=$($MacSurfaceMarginAspectThreshold.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-mac-surface-margin-rho-trigger=$($MacSurfaceMarginRhoTrigger.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-mac-surface-margin-soft-aspect-threshold=$($MacSurfaceMarginSoftAspectThreshold.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-mac-surface-margin-soft-rho-trigger=$($MacSurfaceMarginSoftRhoTrigger.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-mac-surface-margin-soft-factor=$($MacSurfaceMarginSoftFactor.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-leaf-panel-kernel=$LeafPanelKernel"
        $argLines += "-pbe-tabi-cluster-policy=$ClusterPolicy"
        $argLines += "-pbe-tabi-leaf-pair-kernel=$LeafPairKernel"
        $argLines += "-pbe-tabi-leaf-kernel-profile=$LeafKernelProfile"
        $argLines += "-pbe-tabi-leaf-work-profile=$LeafWorkProfile"
        $argLines += "-pbe-tabi-matvec-cost-profile=$MatvecCostProfile"
        $argLines += "-pbe-tabi-jcp-direct-kernel=$JcpDirectKernel"
        $argLines += "-pbe-tabi-target-batch-size=$TargetBatchSize"
        $argLines += "-pbe-tabi-leaf-pair-min-source-panels=$LeafPairMinSourcePanels"
        $argLines += "-pbe-tabi-leaf-pair-min-targets=$LeafPairMinTargets"
        $argLines += "-pbe-tabi-leaf-pair-selection-profile=$LeafPairSelectionProfile"
        $argLines += "-pbe-tabi-leaf-pair-probe-cache=$LeafPairProbeCache"
        $argLines += "-pbe-tabi-leaf-pair-target-group=$LeafPairTargetGroup"
        $argLines += "-pbe-tabi-leaf-pair-max-lookahead=$LeafPairMaxLookahead"
        $argLines += "-pbe-tabi-leaf-pair-directness-probe=$LeafPairDirectnessProbe"
        $argLines += "-pbe-tabi-leaf-pair-directness-cache=$LeafPairDirectnessCache"
        $argLines += "-pbe-tabi-leaf-pair-directness-verify=$LeafPairDirectnessVerify"
        $argLines += "-pbe-tabi-leaf-pair-scheduler=$LeafPairScheduler"
        $argLines += "-pbe-tabi-leaf-pair-worklist-max-active=$LeafPairWorklistMaxActive"
        $argLines += "-pbe-tabi-leaf-pair-worklist-verify=$LeafPairWorklistVerify"
        $argLines += "-pbe-tabi-leaf-block-verify=$LeafBlockVerify"
        $argLines += "-pbe-tabi-expansion-audit=$ExpansionAudit"
        $argLines += "-pbe-tabi-expansion-audit-channel=$ExpansionAuditChannel"
        $argLines += "-pbe-tabi-expansion-audit-sample-count=$ExpansionAuditSampleCount"
        $argLines += "-pbe-tabi-expansion-audit-max-panels=$ExpansionAuditMaxPanels"
        $argLines += '-pbe-tabi-mac-gap-factor=1.5'
        $argLines += '-pbe-tabi-matvec-audit=full-small'
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Invoke-TitanWithTimeout([string]$ArgPath, [string]$RunDir) {
    $stdoutPath = Join-Path $RunDir 'stdout.log'
    $stderrPath = Join-Path $RunDir 'stderr.log'
    $psi = [System.Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = $TitanExe
    $psi.WorkingDirectory = $RunDir
    $psi.UseShellExecute = $false
    $psi.CreateNoWindow = $true
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.Arguments = ('-a "{0}"' -f ($ArgPath -replace '"','\"'))

    $proc = [System.Diagnostics.Process]::new()
    $proc.StartInfo = $psi
    $start = Get-Date
    [void]$proc.Start()
    $stdoutTask = $proc.StandardOutput.ReadToEndAsync()
    $stderrTask = $proc.StandardError.ReadToEndAsync()
    $timedOut = $false
    if ($CaseTimeoutSeconds -gt 0) {
        if (-not $proc.WaitForExit($CaseTimeoutSeconds * 1000)) {
            $timedOut = $true
            $proc.Kill()
            $proc.WaitForExit()
        }
    } else {
        $proc.WaitForExit()
    }
    $stdoutTask.Wait()
    $stderrTask.Wait()
    [System.IO.File]::WriteAllText($stdoutPath, $stdoutTask.Result)
    [System.IO.File]::WriteAllText($stderrPath, $stderrTask.Result)
    $elapsed = ((Get-Date) - $start).TotalSeconds
    $exitCode = if ($timedOut) { 'timeout' } else { [string]$proc.ExitCode }
    return [PSCustomObject]@{
        timed_out = $timedOut
        exit_code = $exitCode
        seconds = $elapsed
        stdout = $stdoutPath
        stderr = $stderrPath
    }
}

function New-ResultRow(
    [string]$Fixture,
    [string]$CaseId,
    [string]$Mode,
    [string]$TreeBuilder,
    [string]$MacSafety,
    [string]$ExpansionCenter,
    [string]$FragmentationSafety,
    [string]$MomentAudit,
    [string]$MomentConstruction,
    [string]$MomentSum,
    [string]$ConvergenceAudit,
    [string]$ExpansionAudit,
    [int]$TaylorOrder,
    [double]$Theta,
    [string]$ConvergenceSafety,
    [double]$ConvergenceRhoCap,
    [double]$MacAspectCap,
    [string]$ConvergenceEstimator,
    [double]$ConvergenceTailThreshold,
    [double]$ConvergenceEstimatorErrorCap,
    [double]$ConvergenceUnstableRhoCap,
    [string]$ConvergenceCalibration = 'stage3j',
    [string]$ConvergenceDielectricAmp = 'off',
    [string]$ConvergenceEstimatorModel = 'stage3l',
    [string]$ConvergenceChannelCalibration = 'stage3m',
    [string]$ConvergenceSplitPolicy = 'global',
    [string]$ConvergenceSplitSelectivity = 'legacy',
    [double]$ConvergenceSplitAbsoluteFloor = 0.0,
    [string]$ConvergenceOverrejectAudit = 'off',
    [string]$ConvergenceOverrejectDetail = 'off',
    [string]$ConvergenceShapeAudit = 'off',
    [string]$ConvergenceKernelResponseAudit = 'off',
    [string]$ChannelMissedSemantics = 'legacy',
    [string]$AcceptedMissDiagnostics = 'off',
    [double]$AcceptedMissAbsoluteFloor = 1.0e-7,
    [double]$K3AcceptedRhoCap = 0.0,
    [string]$MacDistancePolicy = 'center',
    [double]$MacSurfaceMarginFactor = 0.25,
    [string]$MacSurfaceMarginPolicy = 'all',
    [double]$MacSurfaceMarginAspectThreshold = 3.0,
    [double]$MacSurfaceMarginRhoTrigger = 0.55,
    [double]$MacSurfaceMarginSoftAspectThreshold = 0.0,
    [double]$MacSurfaceMarginSoftRhoTrigger = 0.0,
    [double]$MacSurfaceMarginSoftFactor = 0.10,
    [string]$ExpansionAuditChannel = 'off',
    [string]$LeafPanelKernel = 'auto',
    [string]$ClusterPolicy = 'off',
    [string]$LeafPairKernel = 'auto',
    [string]$LeafKernelProfile = 'off',
    [string]$LeafWorkProfile = 'off',
    [string]$MatvecCostProfile = 'off',
    [string]$JcpDirectKernel = 'baseline',
    [int]$TargetBatchSize = 1,
    [int]$LeafPairMinSourcePanels = 64,
    [int]$LeafPairMinTargets = 4,
    [string]$LeafPairSelectionProfile = 'off',
    [string]$LeafPairProbeCache = 'off',
    [string]$LeafPairTargetGroup = 'scan',
    [int]$LeafPairMaxLookahead = 32,
    [string]$LeafPairDirectnessProbe = 'full-traversal',
    [string]$LeafPairDirectnessCache = 'off',
    [string]$LeafPairDirectnessVerify = 'off',
    [string]$LeafPairScheduler = 'inline-select',
    [int]$LeafPairWorklistMaxActive = 2048,
    [string]$LeafPairWorklistVerify = 'off',
    [string]$LeafBlockVerify = 'off',
    [int]$TreeLeafSizeOverride = $TreeLeafSize,
    [string]$Preconditioner = 'diagonal',
    [string]$PreconditionerBlockNear = 'singular',
    [int]$PreconditionerBlockPanelCap = 32,
    [int]$PreconditionerMemoryBudgetMb = 256,
    [int]$PreconditionerCoarseAggregateCap = 16,
    [string]$PreconditionerCoarseBase = 'diagonal',
    [string]$PreconditionerCoarseVerify = 'off',
    [int]$GmresRestart = 32,
    [int]$TrueResidualProbeInterval = 4,
    [string]$RunDir,
    [object]$RunInfo,
    [object]$DirectEnergy
) {
    $timingPath = Join-Path $RunDir 'pbe_timing.csv'
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $memoryPath = Join-Path $RunDir 'pbe_memory_profile.csv'
    $timing = if (Test-Path -LiteralPath $timingPath) { @(Import-Csv -LiteralPath $timingPath)[0] } else { $null }
    $energy = if (Test-Path -LiteralPath $energyPath) { @(Import-Csv -LiteralPath $energyPath)[0] } else { $null }
    $memory = if (Test-Path -LiteralPath $memoryPath) { @(Import-Csv -LiteralPath $memoryPath)[0] } else { $null }
    $energyValue = Convert-ToDoubleOrNull (Get-RowValue $energy 'E_rxn_kJmol')
    $status = if ($RunInfo.timed_out) {
        'timeout'
    } elseif ($RunInfo.exit_code -ne '0') {
        'failed'
    } elseif ($null -eq $timing -or $null -eq $energy) {
        'missing_output'
    } else {
        Get-RowValue $energy 'status' 'ok'
    }

    $driftAbs = $null
    $driftRel = $null
    $directEnergyValue = Convert-ToDoubleOrNull $DirectEnergy
    if ($null -ne $energyValue -and $null -ne $directEnergyValue) {
        $driftAbs = [math]::Abs($energyValue - $directEnergyValue)
        $driftRel = $driftAbs / [math]::Max(1.0e-12, [math]::Abs($directEnergyValue))
    }

    $k2AcceptedMissed = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_convergence_k2_high_error_missed')
    $k3AcceptedMissed = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_convergence_k3_high_error_missed')
    $k2SplitRejects = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_convergence_k2_split_rejects')
    $k3SplitRejects = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_convergence_k3_split_rejects')
    $k2EstimatedError = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_convergence_k2_estimated_error_max')
    $k3EstimatedError = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_convergence_k3_estimated_error_max')
    $k2k3AcceptedMissedDelta = if ($null -ne $k2AcceptedMissed -and $null -ne $k3AcceptedMissed) { $k2AcceptedMissed - $k3AcceptedMissed } else { $null }
    $k2k3SplitRejectDelta = if ($null -ne $k2SplitRejects -and $null -ne $k3SplitRejects) { $k2SplitRejects - $k3SplitRejects } else { $null }
    $k2k3EstimatedErrorRatio = if ($null -ne $k2EstimatedError -and $null -ne $k3EstimatedError -and [math]::Abs($k3EstimatedError) -gt 1.0e-300) { $k2EstimatedError / $k3EstimatedError } else { $null }
    $phiNorm = Get-RowValue $timing 'tabi_phi_norm'
    $phinNorm = Get-RowValue $timing 'tabi_phin_norm'
    $sourceNormS1 = Get-RowValue $timing 'tabi_source_norm_s1'
    $sourceNormS2 = Get-RowValue $timing 'tabi_source_norm_s2'
    $phantomGlobalCount = Get-RowValue $timing 'tabi_convergence_estimator_phantom_high_error_count'
    $phantomGlobalAbsMax = Get-RowValue $timing 'tabi_convergence_estimator_phantom_high_error_abs_max'
    $phantomK1Count = Get-RowValue $timing 'tabi_phantom_high_error_k1_count'
    $phantomK2Count = Get-RowValue $timing 'tabi_phantom_high_error_k2_count'
    $phantomK3Count = Get-RowValue $timing 'tabi_phantom_high_error_k3_count'
    $phantomK4Count = Get-RowValue $timing 'tabi_phantom_high_error_k4_count'
    $phantomK1AbsMax = Get-RowValue $timing 'tabi_phantom_high_error_k1_abs_max'
    $phantomK2AbsMax = Get-RowValue $timing 'tabi_phantom_high_error_k2_abs_max'
    $phantomK3AbsMax = Get-RowValue $timing 'tabi_phantom_high_error_k3_abs_max'
    $phantomK4AbsMax = Get-RowValue $timing 'tabi_phantom_high_error_k4_abs_max'

    return [PSCustomObject]@{
        fixture = $Fixture
        case_id = $CaseId
        mode = $Mode
        tree_builder = $TreeBuilder
        mac_safety = $MacSafety
        expansion_center = $ExpansionCenter
        fragmentation_safety = $FragmentationSafety
        moment_audit = $MomentAudit
        moment_construction = $MomentConstruction
        moment_sum = $MomentSum
        convergence_audit = $ConvergenceAudit
        expansion_audit = $ExpansionAudit
        requested_taylor_order = $TaylorOrder
        taylor_order = Get-RowValue $timing 'tabi_taylor_order' ([string]$TaylorOrder)
        requested_convergence_safety = $ConvergenceSafety
        requested_convergence_rho_cap = $ConvergenceRhoCap.ToString('G17', $InvariantCulture)
        requested_mac_aspect_cap = $MacAspectCap.ToString('G17', $InvariantCulture)
        requested_convergence_estimator = $ConvergenceEstimator
        requested_convergence_tail_threshold = $ConvergenceTailThreshold.ToString('G17', $InvariantCulture)
        requested_convergence_estimator_error_cap = $ConvergenceEstimatorErrorCap.ToString('G17', $InvariantCulture)
        requested_convergence_unstable_rho_cap = $ConvergenceUnstableRhoCap.ToString('G17', $InvariantCulture)
        requested_convergence_calibration = $ConvergenceCalibration
        requested_convergence_dielectric_amp = $ConvergenceDielectricAmp
        requested_convergence_estimator_model = $ConvergenceEstimatorModel
        requested_convergence_channel_calibration = $ConvergenceChannelCalibration
        requested_convergence_split_policy = $ConvergenceSplitPolicy
        requested_convergence_split_selectivity = $ConvergenceSplitSelectivity
        requested_convergence_split_absolute_floor = $ConvergenceSplitAbsoluteFloor.ToString('G17', $InvariantCulture)
        requested_convergence_overreject_audit = $ConvergenceOverrejectAudit
        requested_convergence_overreject_detail = $ConvergenceOverrejectDetail
        requested_convergence_shape_audit = $ConvergenceShapeAudit
        requested_convergence_kernel_response_audit = $ConvergenceKernelResponseAudit
        requested_channel_missed_semantics = $ChannelMissedSemantics
        requested_accepted_miss_diagnostics = $AcceptedMissDiagnostics
        requested_accepted_miss_absolute_floor = $AcceptedMissAbsoluteFloor.ToString('G17', $InvariantCulture)
        requested_k3_accepted_rho_cap = $K3AcceptedRhoCap.ToString('G17', $InvariantCulture)
        requested_mac_distance_policy = $MacDistancePolicy
        requested_mac_surface_margin_factor = $MacSurfaceMarginFactor.ToString('G17', $InvariantCulture)
        requested_mac_surface_margin_policy = $MacSurfaceMarginPolicy
        requested_mac_surface_margin_aspect_threshold = $MacSurfaceMarginAspectThreshold.ToString('G17', $InvariantCulture)
        requested_mac_surface_margin_rho_trigger = $MacSurfaceMarginRhoTrigger.ToString('G17', $InvariantCulture)
        requested_mac_surface_margin_soft_aspect_threshold = $MacSurfaceMarginSoftAspectThreshold.ToString('G17', $InvariantCulture)
        requested_mac_surface_margin_soft_rho_trigger = $MacSurfaceMarginSoftRhoTrigger.ToString('G17', $InvariantCulture)
        requested_mac_surface_margin_soft_factor = $MacSurfaceMarginSoftFactor.ToString('G17', $InvariantCulture)
        requested_expansion_audit_channel = $ExpansionAuditChannel
        requested_leaf_panel_kernel = $LeafPanelKernel
        requested_cluster_policy = $ClusterPolicy
        requested_leaf_pair_kernel = $LeafPairKernel
        requested_leaf_kernel_profile = $LeafKernelProfile
        requested_leaf_work_profile = $LeafWorkProfile
        requested_matvec_cost_profile = $MatvecCostProfile
        requested_jcp_direct_kernel = $JcpDirectKernel
        requested_target_batch_size = $TargetBatchSize
        requested_leaf_pair_min_source_panels = $LeafPairMinSourcePanels
        requested_leaf_pair_min_targets = $LeafPairMinTargets
        requested_leaf_pair_selection_profile = $LeafPairSelectionProfile
        requested_leaf_pair_probe_cache = $LeafPairProbeCache
        requested_leaf_pair_target_group = $LeafPairTargetGroup
        requested_leaf_pair_max_lookahead = $LeafPairMaxLookahead
        requested_leaf_pair_directness_probe = $LeafPairDirectnessProbe
        requested_leaf_pair_directness_cache = $LeafPairDirectnessCache
        requested_leaf_pair_directness_verify = $LeafPairDirectnessVerify
        requested_leaf_pair_scheduler = $LeafPairScheduler
        requested_leaf_pair_worklist_max_active = $LeafPairWorklistMaxActive
        requested_leaf_pair_worklist_verify = $LeafPairWorklistVerify
        requested_leaf_block_verify = $LeafBlockVerify
        requested_preconditioner = $Preconditioner
        requested_preconditioner_block_near = $PreconditionerBlockNear
        requested_preconditioner_block_panel_cap = $PreconditionerBlockPanelCap
        requested_preconditioner_memory_budget_mb = $PreconditionerMemoryBudgetMb
        requested_preconditioner_coarse_aggregate_cap = $PreconditionerCoarseAggregateCap
        requested_preconditioner_coarse_base = $PreconditionerCoarseBase
        requested_preconditioner_coarse_verify = $PreconditionerCoarseVerify
        requested_gmres_restart = $GmresRestart
        requested_true_residual_probe_interval = $TrueResidualProbeInterval
        convergence_safety = Get-RowValue $timing 'tabi_convergence_safety' $ConvergenceSafety
        convergence_rho_cap = Get-RowValue $timing 'tabi_convergence_rho_cap' ($ConvergenceRhoCap.ToString('G17', $InvariantCulture))
        mac_aspect_cap = Get-RowValue $timing 'tabi_mac_aspect_cap' ($MacAspectCap.ToString('G17', $InvariantCulture))
        convergence_estimator = Get-RowValue $timing 'tabi_convergence_estimator' $ConvergenceEstimator
        convergence_tail_threshold = Get-RowValue $timing 'tabi_convergence_tail_threshold' ($ConvergenceTailThreshold.ToString('G17', $InvariantCulture))
        convergence_estimator_error_cap = Get-RowValue $timing 'tabi_convergence_estimator_error_cap' ($ConvergenceEstimatorErrorCap.ToString('G17', $InvariantCulture))
        convergence_unstable_rho_cap = Get-RowValue $timing 'tabi_convergence_unstable_rho_cap' ($ConvergenceUnstableRhoCap.ToString('G17', $InvariantCulture))
        convergence_calibration = Get-RowValue $timing 'tabi_convergence_calibration' $ConvergenceCalibration
        convergence_dielectric_amp = Get-RowValue $timing 'tabi_convergence_dielectric_amp' $ConvergenceDielectricAmp
        convergence_estimator_model = Get-RowValue $timing 'tabi_convergence_estimator_model' $ConvergenceEstimatorModel
        convergence_channel_calibration = Get-RowValue $timing 'tabi_convergence_channel_calibration' $ConvergenceChannelCalibration
        convergence_split_policy = Get-RowValue $timing 'tabi_convergence_split_policy' $ConvergenceSplitPolicy
        convergence_split_selectivity = Get-RowValue $timing 'tabi_convergence_split_selectivity' $ConvergenceSplitSelectivity
        convergence_split_absolute_floor = Get-RowValue $timing 'tabi_convergence_split_absolute_floor' ($ConvergenceSplitAbsoluteFloor.ToString('G17', $InvariantCulture))
        convergence_overreject_audit = Get-RowValue $timing 'tabi_convergence_overreject_audit' $ConvergenceOverrejectAudit
        convergence_overreject_detail = Get-RowValue $timing 'tabi_convergence_overreject_detail' $ConvergenceOverrejectDetail
        convergence_shape_audit = Get-RowValue $timing 'tabi_convergence_shape_audit' $ConvergenceShapeAudit
        convergence_kernel_response_audit = Get-RowValue $timing 'tabi_convergence_kernel_response_audit' $ConvergenceKernelResponseAudit
        channel_missed_semantics = Get-RowValue $timing 'tabi_channel_missed_semantics' $ChannelMissedSemantics
        accepted_miss_diagnostics = Get-RowValue $timing 'tabi_accepted_miss_diagnostics' $AcceptedMissDiagnostics
        accepted_miss_absolute_floor = Get-RowValue $timing 'tabi_accepted_miss_absolute_floor' ($AcceptedMissAbsoluteFloor.ToString('G17', $InvariantCulture))
        k3_accepted_rho_cap = Get-RowValue $timing 'tabi_convergence_k3_accepted_rho_cap' ($K3AcceptedRhoCap.ToString('G17', $InvariantCulture))
        mac_distance_policy = Get-RowValue $timing 'tabi_mac_distance_policy' $MacDistancePolicy
        mac_surface_margin_factor = Get-RowValue $timing 'tabi_mac_surface_margin_factor' ($MacSurfaceMarginFactor.ToString('G17', $InvariantCulture))
        mac_surface_margin_policy = Get-RowValue $timing 'tabi_mac_surface_margin_policy' $MacSurfaceMarginPolicy
        mac_surface_margin_aspect_threshold = Get-RowValue $timing 'tabi_mac_surface_margin_aspect_threshold' ($MacSurfaceMarginAspectThreshold.ToString('G17', $InvariantCulture))
        mac_surface_margin_rho_trigger = Get-RowValue $timing 'tabi_mac_surface_margin_rho_trigger' ($MacSurfaceMarginRhoTrigger.ToString('G17', $InvariantCulture))
        mac_surface_margin_soft_aspect_threshold = Get-RowValue $timing 'tabi_mac_surface_margin_soft_aspect_threshold' ($MacSurfaceMarginSoftAspectThreshold.ToString('G17', $InvariantCulture))
        mac_surface_margin_soft_rho_trigger = Get-RowValue $timing 'tabi_mac_surface_margin_soft_rho_trigger' ($MacSurfaceMarginSoftRhoTrigger.ToString('G17', $InvariantCulture))
        mac_surface_margin_soft_factor = Get-RowValue $timing 'tabi_mac_surface_margin_soft_factor' ($MacSurfaceMarginSoftFactor.ToString('G17', $InvariantCulture))
        expansion_audit_channel = Get-RowValue $timing 'tabi_expansion_audit_channel' $ExpansionAuditChannel
        leaf_panel_kernel = Get-RowValue $timing 'tabi_leaf_panel_kernel' $LeafPanelKernel
        cluster_policy = Get-RowValue $timing 'tabi_cluster_policy' $ClusterPolicy
        leaf_pair_batch_kernel = Get-RowValue $timing 'tabi_leaf_pair_batch_kernel' $LeafPairKernel
        leaf_kernel_profile = Get-RowValue $timing 'tabi_leaf_kernel_profile' $LeafKernelProfile
        leaf_work_profile = Get-RowValue $timing 'tabi_leaf_work_profile' $LeafWorkProfile
        matvec_cost_profile = Get-RowValue $timing 'tabi_matvec_cost_profile' $MatvecCostProfile
        leaf_panels_per_second = Get-RowValue $timing 'tabi_leaf_panels_per_second'
        clean_leaf_panels_per_second = Get-RowValue $timing 'tabi_clean_leaf_panels_per_second'
        tree_traversal_fraction = Get-RowValue $timing 'tabi_tree_traversal_fraction'
        jcp_direct_kernel = Get-RowValue $timing 'tabi_jcp_direct_kernel' $JcpDirectKernel
        leaf_kernel_scalar_tasks = Get-RowValue $timing 'tabi_leaf_kernel_scalar_tasks'
        leaf_kernel_batch_tasks = Get-RowValue $timing 'tabi_leaf_kernel_batch_tasks'
        leaf_kernel_source_panel_evals = Get-RowValue $timing 'tabi_leaf_kernel_source_panel_evals'
        leaf_kernel_scalar_seconds = Get-RowValue $timing 'tabi_leaf_kernel_scalar_seconds'
        leaf_kernel_target_avx2_seconds = Get-RowValue $timing 'tabi_leaf_kernel_target_avx2_seconds'
        leaf_kernel_tail_seconds = Get-RowValue $timing 'tabi_leaf_kernel_tail_seconds'
        leaf_pair_batch_soa_sources = Get-RowValue $timing 'tabi_leaf_pair_batch_soa_sources'
        leaf_pair_batch_soa_seconds = Get-RowValue $timing 'tabi_leaf_pair_batch_soa_seconds'
        leaf_pair_batch_scalar_tail_seconds = Get-RowValue $timing 'tabi_leaf_pair_batch_scalar_tail_seconds'
        leaf_pair_batch_clean_source_panels_per_second = Get-RowValue $timing 'tabi_leaf_pair_batch_clean_source_panels_per_second'
        target_batch_size = Get-RowValue $timing 'tabi_target_batch_size' ([string]$TargetBatchSize)
        leaf_pair_min_source_panels = Get-RowValue $timing 'tabi_leaf_pair_min_source_panels' ([string]$LeafPairMinSourcePanels)
        leaf_pair_min_targets = Get-RowValue $timing 'tabi_leaf_pair_min_targets' ([string]$LeafPairMinTargets)
        leaf_pair_selection_profile = Get-RowValue $timing 'tabi_leaf_pair_selection_profile' $LeafPairSelectionProfile
        leaf_pair_probe_cache = Get-RowValue $timing 'tabi_leaf_pair_probe_cache' $LeafPairProbeCache
        leaf_pair_target_group = Get-RowValue $timing 'tabi_leaf_pair_target_group' $LeafPairTargetGroup
        leaf_pair_max_lookahead = Get-RowValue $timing 'tabi_leaf_pair_max_lookahead' ([string]$LeafPairMaxLookahead)
        leaf_pair_directness_probe = Get-RowValue $timing 'tabi_leaf_pair_directness_probe' $LeafPairDirectnessProbe
        leaf_pair_directness_cache = Get-RowValue $timing 'tabi_leaf_pair_directness_cache' $LeafPairDirectnessCache
        leaf_pair_directness_verify = Get-RowValue $timing 'tabi_leaf_pair_directness_verify' $LeafPairDirectnessVerify
        leaf_pair_scheduler = Get-RowValue $timing 'tabi_leaf_pair_scheduler' $LeafPairScheduler
        leaf_pair_worklist_max_active = Get-RowValue $timing 'tabi_leaf_pair_worklist_max_active' ([string]$LeafPairWorklistMaxActive)
        leaf_pair_worklist_verify = Get-RowValue $timing 'tabi_leaf_pair_worklist_verify' $LeafPairWorklistVerify
        leaf_block_verify = Get-RowValue $timing 'tabi_leaf_block_verify' $LeafBlockVerify
        node_bounds = 'tight-aabb'
        interaction_mode = 'onthefly'
        near_integration = 'duffy_adaptive'
        quad_buffer = 'adaptive'
        preconditioner = Get-RowValue $timing 'tabi_preconditioner_requested' $Preconditioner
        preconditioner_effective = Get-RowValue $timing 'tabi_preconditioner_effective'
        preconditioner_block_near = Get-RowValue $timing 'tabi_preconditioner_block_near' $PreconditionerBlockNear
        preconditioner_block_panel_cap = Get-RowValue $timing 'tabi_preconditioner_block_panel_cap' ([string]$PreconditionerBlockPanelCap)
        preconditioner_memory_budget_bytes = Get-RowValue $timing 'tabi_preconditioner_memory_budget_bytes'
        preconditioner_estimated_bytes = Get-RowValue $timing 'tabi_preconditioner_estimated_bytes'
        preconditioner_valid_blocks = Get-RowValue $timing 'tabi_preconditioner_valid_blocks'
        preconditioner_valid_panels = Get-RowValue $timing 'tabi_preconditioner_valid_panels'
        preconditioner_ill_conditioned_blocks = Get-RowValue $timing 'tabi_preconditioner_ill_conditioned_blocks'
        preconditioner_setup_seconds = Get-RowValue $timing 'tabi_preconditioner_setup_seconds'
        preconditioner_apply_seconds = Get-RowValue $timing 'tabi_preconditioner_apply_seconds'
        preconditioner_apply_calls = Get-RowValue $timing 'tabi_preconditioner_apply_calls'
        preconditioner_coarse_aggregate_cap = Get-RowValue $timing 'tabi_preconditioner_coarse_aggregate_cap' ([string]$PreconditionerCoarseAggregateCap)
        preconditioner_coarse_aggregates = Get-RowValue $timing 'tabi_preconditioner_coarse_aggregates'
        preconditioner_coarse_dimension = Get-RowValue $timing 'tabi_preconditioner_coarse_dimension'
        preconditioner_coarse_base = Get-RowValue $timing 'tabi_preconditioner_coarse_base' $PreconditionerCoarseBase
        preconditioner_coarse_matrix_bytes = Get-RowValue $timing 'tabi_preconditioner_coarse_matrix_bytes'
        preconditioner_coarse_pivot_ratio = Get-RowValue $timing 'tabi_preconditioner_coarse_pivot_ratio'
        preconditioner_coarse_status = Get-RowValue $timing 'tabi_preconditioner_coarse_status'
        preconditioner_coarse_setup_matvec_calls = Get-RowValue $timing 'tabi_preconditioner_coarse_setup_matvec_calls'
        preconditioner_coarse_setup_seconds = Get-RowValue $timing 'tabi_preconditioner_coarse_setup_seconds'
        preconditioner_coarse_apply_seconds = Get-RowValue $timing 'tabi_preconditioner_coarse_apply_seconds'
        preconditioner_coarse_apply_calls = Get-RowValue $timing 'tabi_preconditioner_coarse_apply_calls'
        preconditioner_coarse_rp_identity_error_max = Get-RowValue $timing 'tabi_preconditioner_coarse_rp_identity_error_max'
        preconditioner_coarse_verify_rel_error_max = Get-RowValue $timing 'tabi_preconditioner_coarse_verify_rel_error_max'
        gmres_restart = $GmresRestart
        surface_points = $SurfacePoints
        tree_leaf_size = $TreeLeafSizeOverride
        opening_theta = $Theta.ToString('G17', $InvariantCulture)
        threads = $Threads
        status = $status
        titan_exit_code = $RunInfo.exit_code
        total_seconds = Format-Number $RunInfo.seconds
        energy_kJmol = Format-Number $energyValue
        direct_energy_kJmol = Format-Number $DirectEnergy
        energy_drift_abs_kJmol = Format-Number $driftAbs
        energy_drift_rel = Format-Number $driftRel
        tree_vs_direct_matvec_rel = Get-RowValue $timing 'tabi_tree_vs_direct_matvec_rel'
        matvec_audit_rows = Get-RowValue $timing 'tabi_matvec_audit_rows'
        matvec_audit_seconds = Get-RowValue $timing 'tabi_matvec_audit_seconds'
        true_residual_rel = Get-RowValue $timing 'tabi_true_residual_rel'
        sampled_true_residual_rel = Get-RowValue $timing 'tabi_sampled_true_residual_rel'
        sampled_true_residual_rows = Get-RowValue $timing 'tabi_sampled_true_residual_rows'
        full_true_residual_confirm_calls = Get-RowValue $timing 'tabi_full_true_residual_confirm_calls'
        phi_norm = $phiNorm
        phin_norm = $phinNorm
        source_norm_s1 = $sourceNormS1
        source_norm_s2 = $sourceNormS2
        far_geometric_unsafe_count = Get-RowValue $timing 'tabi_far_geometric_unsafe_count'
        far_geometric_unsafe_panels = Get-RowValue $timing 'tabi_far_geometric_unsafe_panels'
        far_geometric_gap_ratio_min = Get-RowValue $timing 'tabi_far_geometric_gap_ratio_min'
        far_geometric_gap_ratio_mean = Get-RowValue $timing 'tabi_far_geometric_gap_ratio_mean'
        mac_physical_gap_rejects = Get-RowValue $timing 'tabi_mac_physical_gap_rejects'
        mac_convergence_rho_rejects = Get-RowValue $timing 'tabi_mac_convergence_rho_rejects'
        mac_convergence_aspect_rejects = Get-RowValue $timing 'tabi_mac_convergence_aspect_rejects'
        mac_convergence_rho_reject_panels = Get-RowValue $timing 'tabi_mac_convergence_rho_reject_panels'
        mac_convergence_aspect_reject_panels = Get-RowValue $timing 'tabi_mac_convergence_aspect_reject_panels'
        mac_surface_margin_rejects = Get-RowValue $timing 'tabi_mac_surface_margin_rejects'
        mac_surface_margin_reject_panels = Get-RowValue $timing 'tabi_mac_surface_margin_reject_panels'
        mac_surface_margin_applied = Get-RowValue $timing 'tabi_mac_surface_margin_applied'
        mac_surface_margin_skipped = Get-RowValue $timing 'tabi_mac_surface_margin_skipped'
        mac_surface_margin_adaptive_rejects = Get-RowValue $timing 'tabi_mac_surface_margin_adaptive_rejects'
        mac_surface_margin_adaptive_reject_panels = Get-RowValue $timing 'tabi_mac_surface_margin_adaptive_reject_panels'
        mac_surface_margin_skipped_accepted_far = Get-RowValue $timing 'tabi_mac_surface_margin_skipped_accepted_far'
        mac_surface_margin_hard_applied = Get-RowValue $timing 'tabi_mac_surface_margin_hard_applied'
        mac_surface_margin_soft_applied = Get-RowValue $timing 'tabi_mac_surface_margin_soft_applied'
        mac_surface_margin_hard_rejects = Get-RowValue $timing 'tabi_mac_surface_margin_hard_rejects'
        mac_surface_margin_hard_reject_panels = Get-RowValue $timing 'tabi_mac_surface_margin_hard_reject_panels'
        mac_surface_margin_soft_rejects = Get-RowValue $timing 'tabi_mac_surface_margin_soft_rejects'
        mac_surface_margin_soft_reject_panels = Get-RowValue $timing 'tabi_mac_surface_margin_soft_reject_panels'
        convergence_k3_rho_cap_rejects = Get-RowValue $timing 'tabi_convergence_k3_rho_cap_rejects'
        convergence_k3_rho_cap_reject_panels = Get-RowValue $timing 'tabi_convergence_k3_rho_cap_reject_panels'
        convergence_estimator_samples = Get-RowValue $timing 'tabi_convergence_estimator_samples'
        convergence_estimator_unstable_far_count = Get-RowValue $timing 'tabi_convergence_estimator_unstable_far_count'
        convergence_estimator_unstable_far_panels = Get-RowValue $timing 'tabi_convergence_estimator_unstable_far_panels'
        convergence_estimator_split_rejects = Get-RowValue $timing 'tabi_convergence_estimator_split_rejects'
        convergence_estimator_split_reject_panels = Get-RowValue $timing 'tabi_convergence_estimator_split_reject_panels'
        convergence_low_abs_reaccepted_far_count = Get-RowValue $timing 'tabi_convergence_low_abs_reaccepted_far_count'
        convergence_low_abs_reaccepted_far_panels = Get-RowValue $timing 'tabi_convergence_low_abs_reaccepted_far_panels'
        convergence_estimator_tail_ratio_max = Get-RowValue $timing 'tabi_convergence_estimator_tail_ratio_max'
        convergence_estimator_tail_ratio_mean = Get-RowValue $timing 'tabi_convergence_estimator_tail_ratio_mean'
        convergence_estimator_error_max = Get-RowValue $timing 'tabi_convergence_estimator_error_max'
        convergence_estimator_error_mean = Get-RowValue $timing 'tabi_convergence_estimator_error_mean'
        convergence_estimator_abs_error_max = Get-RowValue $timing 'tabi_convergence_estimator_abs_error_max'
        convergence_estimator_abs_error_mean = Get-RowValue $timing 'tabi_convergence_estimator_abs_error_mean'
        convergence_tail_norm_p_max = Get-RowValue $timing 'tabi_convergence_tail_norm_p_max'
        convergence_tail_norm_p_mean = Get-RowValue $timing 'tabi_convergence_tail_norm_p_mean'
        convergence_tail_norm_pm1_max = Get-RowValue $timing 'tabi_convergence_tail_norm_pm1_max'
        convergence_tail_norm_pm1_mean = Get-RowValue $timing 'tabi_convergence_tail_norm_pm1_mean'
        convergence_effective_jcp_amp_max = Get-RowValue $timing 'tabi_convergence_effective_jcp_amp_max'
        convergence_effective_jcp_amp_mean = Get-RowValue $timing 'tabi_convergence_effective_jcp_amp_mean'
        convergence_dielectric_amp_max = Get-RowValue $timing 'tabi_convergence_dielectric_amp_max'
        convergence_dielectric_amp_mean = Get-RowValue $timing 'tabi_convergence_dielectric_amp_mean'
        convergence_k1_effective_jcp_amp_max = Get-RowValue $timing 'tabi_convergence_k1_effective_jcp_amp_max'
        convergence_k2_effective_jcp_amp_max = Get-RowValue $timing 'tabi_convergence_k2_effective_jcp_amp_max'
        convergence_k3_effective_jcp_amp_max = Get-RowValue $timing 'tabi_convergence_k3_effective_jcp_amp_max'
        convergence_k4_effective_jcp_amp_max = Get-RowValue $timing 'tabi_convergence_k4_effective_jcp_amp_max'
        convergence_k1_unstable_far_count = Get-RowValue $timing 'tabi_convergence_k1_unstable_far_count'
        convergence_k2_unstable_far_count = Get-RowValue $timing 'tabi_convergence_k2_unstable_far_count'
        convergence_k3_unstable_far_count = Get-RowValue $timing 'tabi_convergence_k3_unstable_far_count'
        convergence_k4_unstable_far_count = Get-RowValue $timing 'tabi_convergence_k4_unstable_far_count'
        convergence_k1_split_rejects = Get-RowValue $timing 'tabi_convergence_k1_split_rejects'
        convergence_k2_split_rejects = Get-RowValue $timing 'tabi_convergence_k2_split_rejects'
        convergence_k3_split_rejects = Get-RowValue $timing 'tabi_convergence_k3_split_rejects'
        convergence_k4_split_rejects = Get-RowValue $timing 'tabi_convergence_k4_split_rejects'
        convergence_k1_low_abs_reaccepted = Get-RowValue $timing 'tabi_convergence_k1_low_abs_reaccepted'
        convergence_k2_low_abs_reaccepted = Get-RowValue $timing 'tabi_convergence_k2_low_abs_reaccepted'
        convergence_k3_low_abs_reaccepted = Get-RowValue $timing 'tabi_convergence_k3_low_abs_reaccepted'
        convergence_k4_low_abs_reaccepted = Get-RowValue $timing 'tabi_convergence_k4_low_abs_reaccepted'
        convergence_k1_estimated_error_max = Get-RowValue $timing 'tabi_convergence_k1_estimated_error_max'
        convergence_k2_estimated_error_max = Get-RowValue $timing 'tabi_convergence_k2_estimated_error_max'
        convergence_k3_estimated_error_max = Get-RowValue $timing 'tabi_convergence_k3_estimated_error_max'
        convergence_k4_estimated_error_max = Get-RowValue $timing 'tabi_convergence_k4_estimated_error_max'
        k2_k3_accepted_missed_delta = Format-Number $k2k3AcceptedMissedDelta
        k2_k3_split_reject_delta = Format-Number $k2k3SplitRejectDelta
        k2_k3_estimated_error_ratio = Format-Number $k2k3EstimatedErrorRatio
        convergence_k1_high_error_caught = Get-RowValue $timing 'tabi_convergence_k1_high_error_caught'
        convergence_k2_high_error_caught = Get-RowValue $timing 'tabi_convergence_k2_high_error_caught'
        convergence_k3_high_error_caught = Get-RowValue $timing 'tabi_convergence_k3_high_error_caught'
        convergence_k4_high_error_caught = Get-RowValue $timing 'tabi_convergence_k4_high_error_caught'
        convergence_k1_high_error_missed = Get-RowValue $timing 'tabi_convergence_k1_high_error_missed'
        convergence_k2_high_error_missed = Get-RowValue $timing 'tabi_convergence_k2_high_error_missed'
        convergence_k3_high_error_missed = Get-RowValue $timing 'tabi_convergence_k3_high_error_missed'
        convergence_k4_high_error_missed = Get-RowValue $timing 'tabi_convergence_k4_high_error_missed'
        convergence_k1_high_error_caught_by_split = Get-RowValue $timing 'tabi_convergence_k1_high_error_caught_by_split'
        convergence_k2_high_error_caught_by_split = Get-RowValue $timing 'tabi_convergence_k2_high_error_caught_by_split'
        convergence_k3_high_error_caught_by_split = Get-RowValue $timing 'tabi_convergence_k3_high_error_caught_by_split'
        convergence_k4_high_error_caught_by_split = Get-RowValue $timing 'tabi_convergence_k4_high_error_caught_by_split'
        convergence_estimator_high_error_caught = Get-RowValue $timing 'tabi_convergence_estimator_high_error_caught'
        convergence_estimator_high_error_missed = Get-RowValue $timing 'tabi_convergence_estimator_high_error_missed'
        convergence_estimator_high_error_caught_by_split = Get-RowValue $timing 'tabi_convergence_estimator_high_error_caught_by_split'
        convergence_overreject_count = Get-RowValue $timing 'tabi_convergence_overreject_count'
        convergence_overreject_panels = Get-RowValue $timing 'tabi_convergence_overreject_panels'
        convergence_overreject_caught_count = Get-RowValue $timing 'tabi_convergence_overreject_caught_count'
        convergence_overreject_caught_panels = Get-RowValue $timing 'tabi_convergence_overreject_caught_panels'
        convergence_overreject_abs_error_max = Get-RowValue $timing 'tabi_convergence_overreject_abs_error_max'
        convergence_overreject_abs_error_mean = Get-RowValue $timing 'tabi_convergence_overreject_abs_error_mean'
        convergence_overreject_abs_error_p95 = Get-RowValue $timing 'tabi_convergence_overreject_abs_error_p95'
        convergence_overreject_abs_error_p99 = Get-RowValue $timing 'tabi_convergence_overreject_abs_error_p99'
        convergence_overreject_estimated_abs_error_max = Get-RowValue $timing 'tabi_convergence_overreject_estimated_abs_error_max'
        convergence_overreject_estimated_abs_error_mean = Get-RowValue $timing 'tabi_convergence_overreject_estimated_abs_error_mean'
        convergence_overreject_required_scale_max = Get-RowValue $timing 'tabi_convergence_overreject_required_scale_max'
        convergence_overreject_required_scale_p95 = Get-RowValue $timing 'tabi_convergence_overreject_required_scale_p95'
        convergence_overreject_required_scale_p99 = Get-RowValue $timing 'tabi_convergence_overreject_required_scale_p99'
        convergence_k1_overreject_count = Get-RowValue $timing 'tabi_convergence_k1_overreject_count'
        convergence_k2_overreject_count = Get-RowValue $timing 'tabi_convergence_k2_overreject_count'
        convergence_k3_overreject_count = Get-RowValue $timing 'tabi_convergence_k3_overreject_count'
        convergence_k4_overreject_count = Get-RowValue $timing 'tabi_convergence_k4_overreject_count'
        convergence_k1_overreject_caught_count = Get-RowValue $timing 'tabi_convergence_k1_overreject_caught_count'
        convergence_k2_overreject_caught_count = Get-RowValue $timing 'tabi_convergence_k2_overreject_caught_count'
        convergence_k3_overreject_caught_count = Get-RowValue $timing 'tabi_convergence_k3_overreject_caught_count'
        convergence_k4_overreject_caught_count = Get-RowValue $timing 'tabi_convergence_k4_overreject_caught_count'
        convergence_k1_overreject_panels = Get-RowValue $timing 'tabi_convergence_k1_overreject_panels'
        convergence_k2_overreject_panels = Get-RowValue $timing 'tabi_convergence_k2_overreject_panels'
        convergence_k3_overreject_panels = Get-RowValue $timing 'tabi_convergence_k3_overreject_panels'
        convergence_k4_overreject_panels = Get-RowValue $timing 'tabi_convergence_k4_overreject_panels'
        convergence_k1_overreject_abs_error_max = Get-RowValue $timing 'tabi_convergence_k1_overreject_abs_error_max'
        convergence_k2_overreject_abs_error_max = Get-RowValue $timing 'tabi_convergence_k2_overreject_abs_error_max'
        convergence_k3_overreject_abs_error_max = Get-RowValue $timing 'tabi_convergence_k3_overreject_abs_error_max'
        convergence_k4_overreject_abs_error_max = Get-RowValue $timing 'tabi_convergence_k4_overreject_abs_error_max'
        convergence_k1_overreject_abs_error_mean = Get-RowValue $timing 'tabi_convergence_k1_overreject_abs_error_mean'
        convergence_k2_overreject_abs_error_mean = Get-RowValue $timing 'tabi_convergence_k2_overreject_abs_error_mean'
        convergence_k3_overreject_abs_error_mean = Get-RowValue $timing 'tabi_convergence_k3_overreject_abs_error_mean'
        convergence_k4_overreject_abs_error_mean = Get-RowValue $timing 'tabi_convergence_k4_overreject_abs_error_mean'
        convergence_k1_overreject_required_scale_max = Get-RowValue $timing 'tabi_convergence_k1_overreject_required_scale_max'
        convergence_k2_overreject_required_scale_max = Get-RowValue $timing 'tabi_convergence_k2_overreject_required_scale_max'
        convergence_k3_overreject_required_scale_max = Get-RowValue $timing 'tabi_convergence_k3_overreject_required_scale_max'
        convergence_k4_overreject_required_scale_max = Get-RowValue $timing 'tabi_convergence_k4_overreject_required_scale_max'
        convergence_k1_overreject_required_scale_p99 = Get-RowValue $timing 'tabi_convergence_k1_overreject_required_scale_p99'
        convergence_k2_overreject_required_scale_p99 = Get-RowValue $timing 'tabi_convergence_k2_overreject_required_scale_p99'
        convergence_k3_overreject_required_scale_p99 = Get-RowValue $timing 'tabi_convergence_k3_overreject_required_scale_p99'
        convergence_k4_overreject_required_scale_p99 = Get-RowValue $timing 'tabi_convergence_k4_overreject_required_scale_p99'
        convergence_estimator_phantom_high_error_count = $phantomGlobalCount
        convergence_estimator_phantom_high_error_abs_max = $phantomGlobalAbsMax
        phantom_high_error_k1_count = $phantomK1Count
        phantom_high_error_k2_count = $phantomK2Count
        phantom_high_error_k3_count = $phantomK3Count
        phantom_high_error_k4_count = $phantomK4Count
        phantom_high_error_k1_abs_max = $phantomK1AbsMax
        phantom_high_error_k2_abs_max = $phantomK2AbsMax
        phantom_high_error_k3_abs_max = $phantomK3AbsMax
        phantom_high_error_k4_abs_max = $phantomK4AbsMax
        phantom_global_abs_bound = Format-Number (Get-DerivedProduct $phantomGlobalCount $phantomGlobalAbsMax)
        phantom_k1_abs_bound = Format-Number (Get-DerivedProduct $phantomK1Count $phantomK1AbsMax)
        phantom_k2_abs_bound = Format-Number (Get-DerivedProduct $phantomK2Count $phantomK2AbsMax)
        phantom_k3_abs_bound = Format-Number (Get-DerivedProduct $phantomK3Count $phantomK3AbsMax)
        phantom_k4_abs_bound = Format-Number (Get-DerivedProduct $phantomK4Count $phantomK4AbsMax)
        phantom_k3_abs_over_phi_norm = Format-Number (Get-DerivedRatio $phantomK3AbsMax $phiNorm)
        phantom_k3_abs_over_phin_norm = Format-Number (Get-DerivedRatio $phantomK3AbsMax $phinNorm)
        accepted_miss_worst_channel = Get-RowValue $timing 'tabi_accepted_miss_worst_channel'
        accepted_miss_k1_count = Get-RowValue $timing 'tabi_accepted_miss_k1_count'
        accepted_miss_k2_count = Get-RowValue $timing 'tabi_accepted_miss_k2_count'
        accepted_miss_k3_count = Get-RowValue $timing 'tabi_accepted_miss_k3_count'
        accepted_miss_k4_count = Get-RowValue $timing 'tabi_accepted_miss_k4_count'
        accepted_miss_k1_rho_min = Get-RowValue $timing 'tabi_accepted_miss_k1_rho_min'
        accepted_miss_k2_rho_min = Get-RowValue $timing 'tabi_accepted_miss_k2_rho_min'
        accepted_miss_k3_rho_min = Get-RowValue $timing 'tabi_accepted_miss_k3_rho_min'
        accepted_miss_k4_rho_min = Get-RowValue $timing 'tabi_accepted_miss_k4_rho_min'
        accepted_miss_k1_rho_mean = Get-RowValue $timing 'tabi_accepted_miss_k1_rho_mean'
        accepted_miss_k2_rho_mean = Get-RowValue $timing 'tabi_accepted_miss_k2_rho_mean'
        accepted_miss_k3_rho_mean = Get-RowValue $timing 'tabi_accepted_miss_k3_rho_mean'
        accepted_miss_k4_rho_mean = Get-RowValue $timing 'tabi_accepted_miss_k4_rho_mean'
        accepted_miss_k1_rho_max = Get-RowValue $timing 'tabi_accepted_miss_k1_rho_max'
        accepted_miss_k2_rho_max = Get-RowValue $timing 'tabi_accepted_miss_k2_rho_max'
        accepted_miss_k3_rho_max = Get-RowValue $timing 'tabi_accepted_miss_k3_rho_max'
        accepted_miss_k4_rho_max = Get-RowValue $timing 'tabi_accepted_miss_k4_rho_max'
        accepted_miss_k1_tail_ratio_mean = Get-RowValue $timing 'tabi_accepted_miss_k1_tail_ratio_mean'
        accepted_miss_k2_tail_ratio_mean = Get-RowValue $timing 'tabi_accepted_miss_k2_tail_ratio_mean'
        accepted_miss_k3_tail_ratio_mean = Get-RowValue $timing 'tabi_accepted_miss_k3_tail_ratio_mean'
        accepted_miss_k4_tail_ratio_mean = Get-RowValue $timing 'tabi_accepted_miss_k4_tail_ratio_mean'
        accepted_miss_k1_tail_ratio_max = Get-RowValue $timing 'tabi_accepted_miss_k1_tail_ratio_max'
        accepted_miss_k2_tail_ratio_max = Get-RowValue $timing 'tabi_accepted_miss_k2_tail_ratio_max'
        accepted_miss_k3_tail_ratio_max = Get-RowValue $timing 'tabi_accepted_miss_k3_tail_ratio_max'
        accepted_miss_k4_tail_ratio_max = Get-RowValue $timing 'tabi_accepted_miss_k4_tail_ratio_max'
        accepted_miss_k1_estimated_error_mean = Get-RowValue $timing 'tabi_accepted_miss_k1_estimated_error_mean'
        accepted_miss_k2_estimated_error_mean = Get-RowValue $timing 'tabi_accepted_miss_k2_estimated_error_mean'
        accepted_miss_k3_estimated_error_mean = Get-RowValue $timing 'tabi_accepted_miss_k3_estimated_error_mean'
        accepted_miss_k4_estimated_error_mean = Get-RowValue $timing 'tabi_accepted_miss_k4_estimated_error_mean'
        accepted_miss_k1_estimated_error_max = Get-RowValue $timing 'tabi_accepted_miss_k1_estimated_error_max'
        accepted_miss_k2_estimated_error_max = Get-RowValue $timing 'tabi_accepted_miss_k2_estimated_error_max'
        accepted_miss_k3_estimated_error_max = Get-RowValue $timing 'tabi_accepted_miss_k3_estimated_error_max'
        accepted_miss_k4_estimated_error_max = Get-RowValue $timing 'tabi_accepted_miss_k4_estimated_error_max'
        fragmented_far_count = Get-RowValue $timing 'tabi_fragmented_far_count'
        fragmented_far_panels = Get-RowValue $timing 'tabi_fragmented_far_panels'
        fragmented_split_rejects = Get-RowValue $timing 'tabi_fragmented_split_rejects'
        fragmentation_center_gap_ratio_max = Get-RowValue $timing 'tabi_fragmentation_center_gap_ratio_max'
        fragmentation_center_gap_ratio_mean = Get-RowValue $timing 'tabi_fragmentation_center_gap_ratio_mean'
        moment_audit_nodes = Get-RowValue $timing 'tabi_moment_audit_nodes'
        moment_translation_rel_error_max = Get-RowValue $timing 'tabi_moment_translation_rel_error_max'
        moment_translation_rel_error_mean = Get-RowValue $timing 'tabi_moment_translation_rel_error_mean'
        moment_translation_offset_max = Get-RowValue $timing 'tabi_moment_translation_offset_max'
        moment_translation_offset_mean = Get-RowValue $timing 'tabi_moment_translation_offset_mean'
        moment_cancellation_ratio_max = Get-RowValue $timing 'tabi_moment_cancellation_ratio_max'
        moment_cancellation_ratio_mean = Get-RowValue $timing 'tabi_moment_cancellation_ratio_mean'
        moment_dynamic_range_max = Get-RowValue $timing 'tabi_moment_dynamic_range_max'
        far_expansion_audit_pairs = Get-RowValue $timing 'tabi_far_expansion_audit_pairs'
        far_expansion_audit_panels = Get-RowValue $timing 'tabi_far_expansion_audit_panels'
        far_expansion_direct_seconds = Get-RowValue $timing 'tabi_far_expansion_direct_seconds'
        far_expansion_rel_error_max = Get-RowValue $timing 'tabi_far_expansion_rel_error_max'
        far_expansion_rel_error_mean = Get-RowValue $timing 'tabi_far_expansion_rel_error_mean'
        far_expansion_rel_error_p95 = Get-RowValue $timing 'tabi_far_expansion_rel_error_p95'
        far_expansion_rel_error_p99 = Get-RowValue $timing 'tabi_far_expansion_rel_error_p99'
        far_expansion_worst_channel = Get-RowValue $timing 'tabi_far_expansion_worst_channel'
        far_expansion_k1_rel_error_max = Get-RowValue $timing 'tabi_far_expansion_k1_rel_error_max'
        far_expansion_k2_rel_error_max = Get-RowValue $timing 'tabi_far_expansion_k2_rel_error_max'
        far_expansion_k3_rel_error_max = Get-RowValue $timing 'tabi_far_expansion_k3_rel_error_max'
        far_expansion_k4_rel_error_max = Get-RowValue $timing 'tabi_far_expansion_k4_rel_error_max'
        far_expansion_k1_rel_error_mean = Get-RowValue $timing 'tabi_far_expansion_k1_rel_error_mean'
        far_expansion_k2_rel_error_mean = Get-RowValue $timing 'tabi_far_expansion_k2_rel_error_mean'
        far_expansion_k3_rel_error_mean = Get-RowValue $timing 'tabi_far_expansion_k3_rel_error_mean'
        far_expansion_k4_rel_error_mean = Get-RowValue $timing 'tabi_far_expansion_k4_rel_error_mean'
        far_expansion_k1_rel_error_p99 = Get-RowValue $timing 'tabi_far_expansion_k1_rel_error_p99'
        far_expansion_k2_rel_error_p99 = Get-RowValue $timing 'tabi_far_expansion_k2_rel_error_p99'
        far_expansion_k3_rel_error_p99 = Get-RowValue $timing 'tabi_far_expansion_k3_rel_error_p99'
        far_expansion_k4_rel_error_p99 = Get-RowValue $timing 'tabi_far_expansion_k4_rel_error_p99'
        far_expansion_mean_gate_pass = Get-RowValue $timing 'tabi_far_expansion_mean_gate_pass'
        far_convergence_rho_max = Get-RowValue $timing 'tabi_far_convergence_rho_max'
        far_convergence_rho_mean = Get-RowValue $timing 'tabi_far_convergence_rho_mean'
        node_aspect_ratio_max = Get-RowValue $timing 'tabi_node_aspect_ratio_max'
        node_aspect_ratio_mean = Get-RowValue $timing 'tabi_node_aspect_ratio_mean'
        accepted_far_aspect_ratio_max = Get-RowValue $timing 'tabi_accepted_far_aspect_ratio_max'
        accepted_far_aspect_ratio_mean = Get-RowValue $timing 'tabi_accepted_far_aspect_ratio_mean'
        leaf_panel_work_per_matvec = Get-RowValue $timing 'tabi_leaf_panel_work_per_matvec'
        far_node_work_per_matvec = Get-RowValue $timing 'tabi_far_node_work_per_matvec'
        leaf_panel_interactions = Get-RowValue $timing 'tabi_leaf_panel_interactions'
        leaf_work_mac_terminal_panels = Get-RowValue $timing 'tabi_leaf_work_mac_terminal_panels'
        leaf_work_physical_gap_panels = Get-RowValue $timing 'tabi_leaf_work_physical_gap_panels'
        leaf_work_fragmentation_panels = Get-RowValue $timing 'tabi_leaf_work_fragmentation_panels'
        leaf_work_convergence_panels = Get-RowValue $timing 'tabi_leaf_work_convergence_panels'
        leaf_work_convergence_estimator_panels = Get-RowValue $timing 'tabi_leaf_work_convergence_estimator_panels'
        leaf_work_convergence_rho_panels = Get-RowValue $timing 'tabi_leaf_work_convergence_rho_panels'
        leaf_work_convergence_aspect_panels = Get-RowValue $timing 'tabi_leaf_work_convergence_aspect_panels'
        leaf_work_surface_margin_panels = Get-RowValue $timing 'tabi_leaf_work_surface_margin_panels'
        leaf_work_surface_margin_hard_panels = Get-RowValue $timing 'tabi_leaf_work_surface_margin_hard_panels'
        leaf_work_surface_margin_soft_panels = Get-RowValue $timing 'tabi_leaf_work_surface_margin_soft_panels'
        leaf_work_k3_rho_cap_panels = Get-RowValue $timing 'tabi_leaf_work_k3_rho_cap_panels'
        leaf_work_singular_fallback_panels = Get-RowValue $timing 'tabi_leaf_work_singular_fallback_panels'
        leaf_work_scalar_tail_panels = Get-RowValue $timing 'tabi_leaf_work_scalar_tail_panels'
        leaf_work_profile_accounted_panels = Get-RowValue $timing 'tabi_leaf_work_profile_accounted_panels'
        leaf_work_mac_terminal_frac = Get-RowValue $timing 'tabi_leaf_work_mac_terminal_frac'
        leaf_work_convergence_frac = Get-RowValue $timing 'tabi_leaf_work_convergence_frac'
        leaf_work_surface_margin_frac = Get-RowValue $timing 'tabi_leaf_work_surface_margin_frac'
        leaf_work_singular_fallback_frac = Get-RowValue $timing 'tabi_leaf_work_singular_fallback_frac'
        leaf_work_scalar_tail_frac = Get-RowValue $timing 'tabi_leaf_work_scalar_tail_frac'
        leaf_pair_batch_interactions = Get-RowValue $timing 'tabi_leaf_pair_batch_interactions'
        leaf_pair_batch_targets = Get-RowValue $timing 'tabi_leaf_pair_batch_targets'
        leaf_pair_batch_source_panels = Get-RowValue $timing 'tabi_leaf_pair_batch_source_panels'
        leaf_pair_batch_residual_panels = Get-RowValue $timing 'tabi_leaf_pair_batch_residual_panels'
        leaf_pair_batch_fallbacks = Get-RowValue $timing 'tabi_leaf_pair_batch_fallbacks'
        leaf_pair_batch_seconds = Get-RowValue $timing 'tabi_leaf_pair_batch_seconds'
        leaf_pair_batch_vectorized_sources = Get-RowValue $timing 'tabi_leaf_pair_batch_vectorized_sources'
        leaf_pair_batch_vector_lanes = Get-RowValue $timing 'tabi_leaf_pair_batch_vector_lanes'
        leaf_pair_batch_pairs_per_second = Get-RowValue $timing 'tabi_leaf_pair_batch_pairs_per_second'
        leaf_pair_selection_candidate_scans = Get-RowValue $timing 'tabi_leaf_pair_selection_candidate_scans'
        leaf_pair_selection_clean_rejects = Get-RowValue $timing 'tabi_leaf_pair_selection_clean_rejects'
        leaf_pair_selection_probe_calls = Get-RowValue $timing 'tabi_leaf_pair_selection_probe_calls'
        leaf_pair_selection_probe_node_visits = Get-RowValue $timing 'tabi_leaf_pair_selection_probe_node_visits'
        leaf_pair_selection_probe_seconds = Get-RowValue $timing 'tabi_leaf_pair_selection_probe_seconds'
        leaf_pair_selection_seconds = Get-RowValue $timing 'tabi_leaf_pair_selection_seconds'
        leaf_pair_selection_probe_cache_hits = Get-RowValue $timing 'tabi_leaf_pair_selection_probe_cache_hits'
        leaf_pair_selection_probe_cache_misses = Get-RowValue $timing 'tabi_leaf_pair_selection_probe_cache_misses'
        leaf_pair_selection_insufficient_targets = Get-RowValue $timing 'tabi_leaf_pair_selection_insufficient_targets'
        leaf_pair_selection_same_target_leaf_rejects = Get-RowValue $timing 'tabi_leaf_pair_selection_same_target_leaf_rejects'
        leaf_pair_directness_ancestor_probe_calls = Get-RowValue $timing 'tabi_leaf_pair_directness_ancestor_probe_calls'
        leaf_pair_directness_ancestor_probe_node_visits = Get-RowValue $timing 'tabi_leaf_pair_directness_ancestor_probe_node_visits'
        leaf_pair_directness_ancestor_probe_seconds = Get-RowValue $timing 'tabi_leaf_pair_directness_ancestor_probe_seconds'
        leaf_pair_directness_pair_cache_hits = Get-RowValue $timing 'tabi_leaf_pair_directness_pair_cache_hits'
        leaf_pair_directness_pair_cache_misses = Get-RowValue $timing 'tabi_leaf_pair_directness_pair_cache_misses'
        leaf_pair_directness_all_direct_proofs = Get-RowValue $timing 'tabi_leaf_pair_directness_all_direct_proofs'
        leaf_pair_directness_mixed_fallbacks = Get-RowValue $timing 'tabi_leaf_pair_directness_mixed_fallbacks'
        leaf_pair_directness_verify_checks = Get-RowValue $timing 'tabi_leaf_pair_directness_verify_checks'
        leaf_pair_directness_verify_mismatches = Get-RowValue $timing 'tabi_leaf_pair_directness_verify_mismatches'
        leaf_pair_worklist_tasks = Get-RowValue $timing 'tabi_leaf_pair_worklist_tasks'
        leaf_pair_worklist_groups = Get-RowValue $timing 'tabi_leaf_pair_worklist_groups'
        leaf_pair_worklist_batched_groups = Get-RowValue $timing 'tabi_leaf_pair_worklist_batched_groups'
        leaf_pair_worklist_scalar_flush_tasks = Get-RowValue $timing 'tabi_leaf_pair_worklist_scalar_flush_tasks'
        leaf_pair_worklist_forced_flushes = Get-RowValue $timing 'tabi_leaf_pair_worklist_forced_flushes'
        leaf_pair_worklist_active_groups_max = Get-RowValue $timing 'tabi_leaf_pair_worklist_active_groups_max'
        leaf_pair_worklist_seconds = Get-RowValue $timing 'tabi_leaf_pair_worklist_seconds'
        leaf_pair_worklist_verify_checks = Get-RowValue $timing 'tabi_leaf_pair_worklist_verify_checks'
        leaf_pair_worklist_verify_mismatches = Get-RowValue $timing 'tabi_leaf_pair_worklist_verify_mismatches'
        leaf_block_target_leaves = Get-RowValue $timing 'tabi_leaf_block_target_leaves'
        leaf_block_targets = Get-RowValue $timing 'tabi_leaf_block_targets'
        leaf_block_source_groups = Get-RowValue $timing 'tabi_leaf_block_source_groups'
        leaf_block_batched_groups = Get-RowValue $timing 'tabi_leaf_block_batched_groups'
        leaf_block_scalar_tail_tasks = Get-RowValue $timing 'tabi_leaf_block_scalar_tail_tasks'
        leaf_block_group_seconds = Get-RowValue $timing 'tabi_leaf_block_group_seconds'
        leaf_block_flush_seconds = Get-RowValue $timing 'tabi_leaf_block_flush_seconds'
        leaf_block_writeback_seconds = Get-RowValue $timing 'tabi_leaf_block_writeback_seconds'
        leaf_block_verify_checks = Get-RowValue $timing 'tabi_leaf_block_verify_checks'
        leaf_block_verify_mismatches = Get-RowValue $timing 'tabi_leaf_block_verify_mismatches'
        far_node_interactions = Get-RowValue $timing 'tabi_far_node_interactions'
        matvec_calls = Get-RowValue $timing 'tabi_matvec_calls'
        rss_peak_bytes = Get-RowValue $memory 'peak_working_set_bytes'
        estimated_bytes_total = Get-RowValue $memory 'estimated_bytes_total'
        matvec_reduction_vs_diagonal = 'NA'
        gmres_matvec_reduction_excluding_setup = 'NA'
        setup_amortization_matvec_equivalent = 'NA'
        time_speedup_vs_fused_diagonal = 'NA'
        peak_memory_delta_vs_diagonal = 'NA'
        performance_classification = 'pending'
        run_dir = $RunDir
    }
}

function Update-Classifications([object[]]$Rows) {
    foreach ($fixture in @($Rows | Select-Object -ExpandProperty fixture -Unique)) {
        $fixtureRows = @($Rows | Where-Object { $_.fixture -eq $fixture })
        $grouped = @($fixtureRows | Where-Object {
            $_.tree_builder -eq 'grouped-morton' -and $_.mac_safety -eq 'off'
        } | Select-Object -First 1)
        $groupedSampled = if ($grouped.Count -gt 0) {
            Convert-ToDoubleOrNull $grouped[0].sampled_true_residual_rel
        } else {
            $null
        }
        $stage3ajDiagonal = @($fixtureRows | Where-Object {
            $_.case_id -eq 'stage3aj_fused_diagonal_r64_control' -and $_.status -eq 'ok'
        } | Select-Object -First 1)
        $stage3akDiagonal = @($fixtureRows | Where-Object {
            $_.case_id -eq 'stage3ak_fused_diagonal_interval1_control' -and $_.status -eq 'ok'
        } | Select-Object -First 1)
        foreach ($row in $fixtureRows) {
            if ($stage3ajDiagonal.Count -gt 0 -and $row.case_id -like 'stage3aj_*') {
                $baselineMatvec = Convert-ToDoubleOrNull $stage3ajDiagonal[0].matvec_calls
                $rowMatvec = Convert-ToDoubleOrNull $row.matvec_calls
                $baselineSeconds = Convert-ToDoubleOrNull $stage3ajDiagonal[0].total_seconds
                $rowSeconds = Convert-ToDoubleOrNull $row.total_seconds
                $baselinePeak = Convert-ToDoubleOrNull $stage3ajDiagonal[0].rss_peak_bytes
                $rowPeak = Convert-ToDoubleOrNull $row.rss_peak_bytes
                if ($null -ne $baselineMatvec -and $baselineMatvec -gt 0.0 -and $null -ne $rowMatvec) {
                    $row.matvec_reduction_vs_diagonal = Format-Number (1.0 - ($rowMatvec / $baselineMatvec))
                }
                if ($null -ne $baselineSeconds -and $null -ne $rowSeconds -and $rowSeconds -gt 0.0) {
                    $row.time_speedup_vs_fused_diagonal = Format-Number ($baselineSeconds / $rowSeconds)
                }
                if ($null -ne $baselinePeak -and $null -ne $rowPeak) {
                    $row.peak_memory_delta_vs_diagonal = Format-Number ($rowPeak - $baselinePeak)
                }
            }
            if ($stage3akDiagonal.Count -gt 0 -and $row.case_id -like 'stage3ak_*') {
                $baselineMatvec = Convert-ToDoubleOrNull $stage3akDiagonal[0].matvec_calls
                $rowMatvec = Convert-ToDoubleOrNull $row.matvec_calls
                $baselineSeconds = Convert-ToDoubleOrNull $stage3akDiagonal[0].total_seconds
                $rowSeconds = Convert-ToDoubleOrNull $row.total_seconds
                $setupSeconds = Convert-ToDoubleOrNull $row.preconditioner_coarse_setup_seconds
                if ($null -ne $baselineMatvec -and $baselineMatvec -gt 0.0 -and $null -ne $rowMatvec) {
                    $row.gmres_matvec_reduction_excluding_setup = Format-Number (1.0 - ($rowMatvec / $baselineMatvec))
                    $row.matvec_reduction_vs_diagonal = $row.gmres_matvec_reduction_excluding_setup
                }
                if ($null -ne $baselineSeconds -and $null -ne $baselineMatvec -and $baselineMatvec -gt 0.0 -and
                    $null -ne $setupSeconds -and $baselineSeconds -gt 0.0) {
                    $row.setup_amortization_matvec_equivalent =
                        Format-Number ($setupSeconds / ($baselineSeconds / $baselineMatvec))
                }
                if ($null -ne $baselineSeconds -and $null -ne $rowSeconds -and $rowSeconds -gt 0.0) {
                    $row.time_speedup_vs_fused_diagonal = Format-Number ($baselineSeconds / $rowSeconds)
                }
            }
            if ($row.mode -eq 'direct') {
                $row.performance_classification = 'direct_baseline'
                continue
            }
            if ($row.status -ne 'ok') {
                $row.performance_classification = $row.status
                continue
            }
            $matvecRel = Convert-ToDoubleOrNull $row.tree_vs_direct_matvec_rel
            $driftRel = Convert-ToDoubleOrNull $row.energy_drift_rel
            $sampledRel = Convert-ToDoubleOrNull $row.sampled_true_residual_rel
            $farMean = Convert-ToDoubleOrNull $row.far_expansion_rel_error_mean
            $farP99 = Convert-ToDoubleOrNull $row.far_expansion_rel_error_p99
            $missed = Convert-ToDoubleOrNull $row.convergence_estimator_high_error_missed
            $channelMissed = 0.0
            foreach ($channelMissName in @('convergence_k1_high_error_missed','convergence_k2_high_error_missed','convergence_k3_high_error_missed','convergence_k4_high_error_missed')) {
                $value = Convert-ToDoubleOrNull $row.$channelMissName
                if ($null -ne $value) { $channelMissed += $value }
            }
            $unsafe = Convert-ToDoubleOrNull $row.far_geometric_unsafe_count
            $rejects = Convert-ToDoubleOrNull $row.mac_physical_gap_rejects
            $matvecPass = $null -ne $matvecRel -and $matvecRel -le 2.5e-3
            $energyPass = $null -ne $driftRel -and $driftRel -le 5.0e-3
            $meanGateFlag = ([string]$row.far_expansion_mean_gate_pass) -eq '1'
            $meanLimit = if ([string]$row.requested_convergence_dielectric_amp -eq 'on') { 5.0e-5 } else { 5.0e-4 }
            $p99GatePass = if ([string]$row.requested_expansion_audit_channel -eq 'off') {
                $true
            } elseif ([string]$row.requested_channel_missed_semantics -eq 'accepted-only') {
                $true
            } else {
                $null -ne $farP99 -and $farP99 -lt 5.0e-4
            }
            $meanGatePass = $meanGateFlag -and $null -ne $farMean -and $farMean -lt $meanLimit -and
                (($null -eq $missed) -or $missed -eq 0.0) -and $channelMissed -eq 0.0 -and $p99GatePass
            $residualPass = $false
            if ($null -ne $sampledRel -and $null -ne $groupedSampled) {
                $residualPass = $sampledRel -le (2.0 * [math]::Max(1.0e-16, $groupedSampled))
            }
            if ($row.requested_convergence_estimator -ne 'off' -and $row.requested_convergence_safety -eq 'split' -and $matvecPass -and $energyPass -and $residualPass -and $meanGatePass) {
                $row.performance_classification = 'convergence_mean_gate_medium_pass'
            } elseif ($row.requested_convergence_estimator -ne 'off' -and $row.requested_convergence_safety -eq 'split' -and $matvecPass -and $energyPass -and $residualPass) {
                if ([string]$row.requested_expansion_audit_channel -ne 'off') {
                    $row.performance_classification = 'convergence_mean_gate_followup'
                } else {
                    $row.performance_classification = 'convergence_estimator_medium_gate_pass'
                }
            } elseif ($row.requested_convergence_estimator -ne 'off' -and $row.requested_convergence_safety -eq 'split' -and $matvecPass -and $energyPass) {
                $row.performance_classification = 'convergence_estimator_residual_followup'
            } elseif ($row.requested_convergence_safety -eq 'split' -and $matvecPass -and $energyPass -and $residualPass) {
                $row.performance_classification = 'convergence_split_medium_gate_pass'
            } elseif ($row.fragmentation_safety -eq 'split' -and $matvecPass -and $energyPass -and $residualPass) {
                $row.performance_classification = 'fragmentation_split_medium_gate_pass'
            } elseif ($row.requested_convergence_safety -eq 'split' -and $matvecPass -and $energyPass) {
                $row.performance_classification = 'convergence_split_residual_followup'
            } elseif ($row.mac_safety -eq 'physical-gap' -and $matvecPass -and $energyPass -and $residualPass) {
                $row.performance_classification = 'physical_gap_medium_gate_pass'
            } elseif ($row.fragmentation_safety -eq 'split' -and $matvecPass -and $energyPass) {
                $row.performance_classification = 'fragmentation_split_residual_followup'
            } elseif ($row.requested_convergence_safety -eq 'split') {
                $row.performance_classification = 'convergence_split_accuracy_followup'
            } elseif ($row.mac_safety -eq 'physical-gap' -and $matvecPass -and $energyPass) {
                $row.performance_classification = 'physical_gap_residual_followup'
            } elseif ($row.fragmentation_safety -eq 'split') {
                $row.performance_classification = 'fragmentation_split_accuracy_followup'
            } elseif ($row.mac_safety -eq 'physical-gap' -and
                      (($null -ne $rejects -and $rejects -gt 0.0) -or
                       ($null -ne $unsafe -and $unsafe -eq 0.0))) {
                $row.performance_classification = 'moment_or_mac_followup'
            } else {
                $row.performance_classification = 'tree_accuracy_regression'
            }
        }
    }
}

if ($Quick -and -not $PSBoundParameters.ContainsKey('Fixtures')) {
    $Fixtures = @('mini-protein')
}
if ($Quick -and $CaseTimeoutSeconds -le 0) {
    $CaseTimeoutSeconds = 90
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}

New-CleanDir $WorkRoot
$oldNoPause = $env:TITAN_NO_PAUSE
$env:TITAN_NO_PAUSE = '1'
$rows = @()
try {
    foreach ($fixture in $Fixtures) {
        $fixtureName = $fixture.Trim().ToLowerInvariant()
        if ([string]::IsNullOrWhiteSpace($fixtureName)) { continue }
        $directRunDir = Join-Path $WorkRoot (Join-Path $fixtureName 'direct_jcp2013')
        New-CleanDir $directRunDir
        $directArg = New-ArgFile -RunDir $directRunDir -Fixture $fixtureName `
            -Mode 'direct' -TreeBuilder 'direct' -MacSafety 'off' `
            -ExpansionCenter 'legacy' -FragmentationSafety 'off' `
            -MomentAudit 'off' -MomentConstruction 'translated' `
            -MomentSum 'naive' -ConvergenceAudit 'off' `
            -ExpansionAudit 'off' -TaylorOrder 3 -Theta $OpeningTheta `
            -ConvergenceSafety 'off' -ConvergenceRhoCap 0.0 -MacAspectCap 0.0
        $directInfo = Invoke-TitanWithTimeout -ArgPath $directArg -RunDir $directRunDir
        $directRow = New-ResultRow -Fixture $fixtureName -CaseId 'direct_jcp2013' `
            -Mode 'direct' -TreeBuilder 'direct' -MacSafety 'off' `
            -ExpansionCenter 'legacy' -FragmentationSafety 'off' `
            -MomentAudit 'off' -MomentConstruction 'translated' `
            -MomentSum 'naive' -ConvergenceAudit 'off' `
            -ExpansionAudit 'off' -TaylorOrder 3 -Theta $OpeningTheta `
            -ConvergenceSafety 'off' -ConvergenceRhoCap 0.0 -MacAspectCap 0.0 `
            -ConvergenceEstimator 'off' -ConvergenceTailThreshold 1.0e4 `
            -ConvergenceEstimatorErrorCap 1.0e-3 `
            -ConvergenceUnstableRhoCap 0.25 `
            -RunDir $directRunDir -RunInfo $directInfo -DirectEnergy $null
        $directEnergy = Convert-ToDoubleOrNull $directRow.energy_kJmol
        if ($null -ne $directEnergy) {
            $directRow.direct_energy_kJmol = Format-Number $directEnergy
            $directRow.energy_drift_abs_kJmol = '0'
            $directRow.energy_drift_rel = '0'
        }
        $rows += $directRow
        if ($directRow.status -ne 'ok' -and -not $AllowFailedStatus) {
            throw "Direct baseline failed for ${fixtureName}: $($directRow.status)"
        }

        foreach ($tree in $TreeRows) {
            $caseId = [string]$tree.CaseId
            $rowTheta = if ($tree.ContainsKey('Theta')) { [double]$tree.Theta } else { $OpeningTheta }
            $rowTaylorOrder = if ($tree.ContainsKey('TaylorOrder')) { [int]$tree.TaylorOrder } else { 3 }
            $rowConvergenceSafety = if ($tree.ContainsKey('ConvergenceSafety')) { [string]$tree.ConvergenceSafety } else { 'off' }
            $rowConvergenceRhoCap = if ($tree.ContainsKey('ConvergenceRhoCap')) { [double]$tree.ConvergenceRhoCap } else { 0.0 }
            $rowMacAspectCap = if ($tree.ContainsKey('MacAspectCap')) { [double]$tree.MacAspectCap } else { 0.0 }
            $rowConvergenceEstimator = if ($tree.ContainsKey('ConvergenceEstimator')) { [string]$tree.ConvergenceEstimator } else { 'off' }
            $rowConvergenceTailThreshold = if ($tree.ContainsKey('ConvergenceTailThreshold')) { [double]$tree.ConvergenceTailThreshold } else { 1.0e4 }
            $rowConvergenceEstimatorErrorCap = if ($tree.ContainsKey('ConvergenceEstimatorErrorCap')) { [double]$tree.ConvergenceEstimatorErrorCap } else { 1.0e-3 }
            $rowConvergenceUnstableRhoCap = if ($tree.ContainsKey('ConvergenceUnstableRhoCap')) { [double]$tree.ConvergenceUnstableRhoCap } else { 0.25 }
            $rowConvergenceCalibration = if ($tree.ContainsKey('ConvergenceCalibration')) { [string]$tree.ConvergenceCalibration } else { 'stage3j' }
            $rowConvergenceDielectricAmp = if ($tree.ContainsKey('ConvergenceDielectricAmp')) { [string]$tree.ConvergenceDielectricAmp } else { 'off' }
            $rowConvergenceEstimatorModel = if ($tree.ContainsKey('ConvergenceEstimatorModel')) { [string]$tree.ConvergenceEstimatorModel } else { 'stage3l' }
            $rowConvergenceChannelCalibration = if ($tree.ContainsKey('ConvergenceChannelCalibration')) { [string]$tree.ConvergenceChannelCalibration } else { 'stage3m' }
            $rowConvergenceSplitPolicy = if ($tree.ContainsKey('ConvergenceSplitPolicy')) { [string]$tree.ConvergenceSplitPolicy } else { 'global' }
            $rowConvergenceSplitSelectivity = if ($tree.ContainsKey('ConvergenceSplitSelectivity')) { [string]$tree.ConvergenceSplitSelectivity } else { 'legacy' }
            $rowConvergenceSplitAbsoluteFloor = if ($tree.ContainsKey('ConvergenceSplitAbsoluteFloor')) { [double]$tree.ConvergenceSplitAbsoluteFloor } else { 0.0 }
            $rowConvergenceOverrejectAudit = if ($tree.ContainsKey('ConvergenceOverrejectAudit')) { [string]$tree.ConvergenceOverrejectAudit } else { 'off' }
            $rowConvergenceOverrejectDetail = if ($tree.ContainsKey('ConvergenceOverrejectDetail')) { [string]$tree.ConvergenceOverrejectDetail } else { 'off' }
            $rowConvergenceShapeAudit = if ($tree.ContainsKey('ConvergenceShapeAudit')) { [string]$tree.ConvergenceShapeAudit } else { 'off' }
            $rowConvergenceKernelResponseAudit = if ($tree.ContainsKey('ConvergenceKernelResponseAudit')) { [string]$tree.ConvergenceKernelResponseAudit } else { 'off' }
            $rowChannelMissedSemantics = if ($tree.ContainsKey('ChannelMissedSemantics')) { [string]$tree.ChannelMissedSemantics } else { 'legacy' }
            $rowAcceptedMissDiagnostics = if ($tree.ContainsKey('AcceptedMissDiagnostics')) { [string]$tree.AcceptedMissDiagnostics } else { 'off' }
            $rowAcceptedMissAbsoluteFloor = if ($tree.ContainsKey('AcceptedMissAbsoluteFloor')) { [double]$tree.AcceptedMissAbsoluteFloor } else { 1.0e-7 }
            $rowK3AcceptedRhoCap = if ($tree.ContainsKey('K3AcceptedRhoCap')) { [double]$tree.K3AcceptedRhoCap } else { 0.0 }
            $rowMacDistancePolicy = if ($tree.ContainsKey('MacDistancePolicy')) { [string]$tree.MacDistancePolicy } else { 'center' }
            $rowMacSurfaceMarginFactor = if ($tree.ContainsKey('MacSurfaceMarginFactor')) { [double]$tree.MacSurfaceMarginFactor } else { 0.25 }
            $rowMacSurfaceMarginPolicy = if ($tree.ContainsKey('MacSurfaceMarginPolicy')) { [string]$tree.MacSurfaceMarginPolicy } else { 'all' }
            $rowMacSurfaceMarginAspectThreshold = if ($tree.ContainsKey('MacSurfaceMarginAspectThreshold')) { [double]$tree.MacSurfaceMarginAspectThreshold } else { 3.0 }
            $rowMacSurfaceMarginRhoTrigger = if ($tree.ContainsKey('MacSurfaceMarginRhoTrigger')) { [double]$tree.MacSurfaceMarginRhoTrigger } else { 0.55 }
            $rowMacSurfaceMarginSoftAspectThreshold = if ($tree.ContainsKey('MacSurfaceMarginSoftAspectThreshold')) { [double]$tree.MacSurfaceMarginSoftAspectThreshold } else { 0.0 }
            $rowMacSurfaceMarginSoftRhoTrigger = if ($tree.ContainsKey('MacSurfaceMarginSoftRhoTrigger')) { [double]$tree.MacSurfaceMarginSoftRhoTrigger } else { 0.0 }
            $rowMacSurfaceMarginSoftFactor = if ($tree.ContainsKey('MacSurfaceMarginSoftFactor')) { [double]$tree.MacSurfaceMarginSoftFactor } else { 0.10 }
            $rowExpansionAuditChannel = if ($tree.ContainsKey('ExpansionAuditChannel')) { [string]$tree.ExpansionAuditChannel } else { 'off' }
            $rowLeafPanelKernel = if ($tree.ContainsKey('LeafPanelKernel')) { [string]$tree.LeafPanelKernel } else { 'auto' }
            $rowClusterPolicy = if ($tree.ContainsKey('ClusterPolicy')) { [string]$tree.ClusterPolicy } else { 'off' }
            $rowLeafPairKernel = if ($tree.ContainsKey('LeafPairKernel')) { [string]$tree.LeafPairKernel } else { 'auto' }
            $rowLeafKernelProfile = if ($tree.ContainsKey('LeafKernelProfile')) { [string]$tree.LeafKernelProfile } else { 'off' }
            $rowLeafWorkProfile = if ($tree.ContainsKey('LeafWorkProfile')) { [string]$tree.LeafWorkProfile } else { 'off' }
            $rowMatvecCostProfile = if ($tree.ContainsKey('MatvecCostProfile')) { [string]$tree.MatvecCostProfile } else { 'off' }
            $rowJcpDirectKernel = if ($tree.ContainsKey('JcpDirectKernel')) { [string]$tree.JcpDirectKernel } else { 'baseline' }
            $rowTargetBatchSize = if ($tree.ContainsKey('TargetBatchSize')) { [int]$tree.TargetBatchSize } else { 1 }
            $rowLeafPairMinSourcePanels = if ($tree.ContainsKey('LeafPairMinSourcePanels')) { [int]$tree.LeafPairMinSourcePanels } else { 64 }
            $rowLeafPairMinTargets = if ($tree.ContainsKey('LeafPairMinTargets')) { [int]$tree.LeafPairMinTargets } else { 4 }
            $rowLeafPairSelectionProfile = if ($tree.ContainsKey('LeafPairSelectionProfile')) { [string]$tree.LeafPairSelectionProfile } else { 'off' }
            $rowLeafPairProbeCache = if ($tree.ContainsKey('LeafPairProbeCache')) { [string]$tree.LeafPairProbeCache } else { 'off' }
            $rowLeafPairTargetGroup = if ($tree.ContainsKey('LeafPairTargetGroup')) { [string]$tree.LeafPairTargetGroup } else { 'scan' }
            $rowLeafPairMaxLookahead = if ($tree.ContainsKey('LeafPairMaxLookahead')) { [int]$tree.LeafPairMaxLookahead } else { 32 }
            $rowLeafPairDirectnessProbe = if ($tree.ContainsKey('LeafPairDirectnessProbe')) { [string]$tree.LeafPairDirectnessProbe } else { 'full-traversal' }
            $rowLeafPairDirectnessCache = if ($tree.ContainsKey('LeafPairDirectnessCache')) { [string]$tree.LeafPairDirectnessCache } else { 'off' }
            $rowLeafPairDirectnessVerify = if ($tree.ContainsKey('LeafPairDirectnessVerify')) { [string]$tree.LeafPairDirectnessVerify } else { 'off' }
            $rowLeafPairScheduler = if ($tree.ContainsKey('LeafPairScheduler')) { [string]$tree.LeafPairScheduler } else { 'inline-select' }
            $rowLeafPairWorklistMaxActive = if ($tree.ContainsKey('LeafPairWorklistMaxActive')) { [int]$tree.LeafPairWorklistMaxActive } else { 2048 }
            $rowLeafPairWorklistVerify = if ($tree.ContainsKey('LeafPairWorklistVerify')) { [string]$tree.LeafPairWorklistVerify } else { 'off' }
            $rowLeafBlockVerify = if ($tree.ContainsKey('LeafBlockVerify')) { [string]$tree.LeafBlockVerify } else { 'off' }
            $rowTreeLeafSize = if ($tree.ContainsKey('LeafSize')) { [int]$tree.LeafSize } else { $TreeLeafSize }
            $rowPreconditioner = if ($tree.ContainsKey('Preconditioner')) { [string]$tree.Preconditioner } else { 'diagonal' }
            $rowPreconditionerBlockNear = if ($tree.ContainsKey('PreconditionerBlockNear')) { [string]$tree.PreconditionerBlockNear } else { 'singular' }
            $rowPreconditionerBlockPanelCap = if ($tree.ContainsKey('PreconditionerBlockPanelCap')) { [int]$tree.PreconditionerBlockPanelCap } else { 32 }
            $rowPreconditionerMemoryBudgetMb = if ($tree.ContainsKey('PreconditionerMemoryBudgetMb')) { [int]$tree.PreconditionerMemoryBudgetMb } else { 256 }
            $rowPreconditionerCoarseAggregateCap = if ($tree.ContainsKey('PreconditionerCoarseAggregateCap')) { [int]$tree.PreconditionerCoarseAggregateCap } else { 16 }
            $rowPreconditionerCoarseBase = if ($tree.ContainsKey('PreconditionerCoarseBase')) { [string]$tree.PreconditionerCoarseBase } else { 'diagonal' }
            $rowPreconditionerCoarseVerify = if ($tree.ContainsKey('PreconditionerCoarseVerify')) { [string]$tree.PreconditionerCoarseVerify } else { 'off' }
            $rowGmresRestart = if ($tree.ContainsKey('GmresRestart')) { [int]$tree.GmresRestart } else { 32 }
            $rowTrueResidualProbeInterval = if ($tree.ContainsKey('TrueResidualProbeInterval')) { [int]$tree.TrueResidualProbeInterval } else { 4 }
            $runDir = Join-Path $WorkRoot (Join-Path $fixtureName $caseId)
            New-CleanDir $runDir
            $arg = New-ArgFile -RunDir $runDir -Fixture $fixtureName `
                -Mode 'treecode' -TreeBuilder $tree.TreeBuilder `
                -MacSafety $tree.MacSafety `
                -ExpansionCenter $tree.ExpansionCenter `
                -FragmentationSafety $tree.FragmentationSafety `
                -MomentAudit $tree.MomentAudit `
                -MomentConstruction $tree.MomentConstruction `
                -MomentSum $tree.MomentSum `
                -ConvergenceAudit $tree.ConvergenceAudit `
                -ExpansionAudit $tree.ExpansionAudit `
                -TaylorOrder $rowTaylorOrder -Theta $rowTheta `
                -ConvergenceSafety $rowConvergenceSafety `
                -ConvergenceRhoCap $rowConvergenceRhoCap `
                -MacAspectCap $rowMacAspectCap `
                -ConvergenceEstimator $rowConvergenceEstimator `
                -ConvergenceTailThreshold $rowConvergenceTailThreshold `
                -ConvergenceEstimatorErrorCap $rowConvergenceEstimatorErrorCap `
                -ConvergenceUnstableRhoCap $rowConvergenceUnstableRhoCap `
                -ConvergenceCalibration $rowConvergenceCalibration `
                -ConvergenceDielectricAmp $rowConvergenceDielectricAmp `
                -ConvergenceEstimatorModel $rowConvergenceEstimatorModel `
                -ConvergenceChannelCalibration $rowConvergenceChannelCalibration `
                -ConvergenceSplitPolicy $rowConvergenceSplitPolicy `
                -ConvergenceSplitSelectivity $rowConvergenceSplitSelectivity `
                -ConvergenceSplitAbsoluteFloor $rowConvergenceSplitAbsoluteFloor `
                -ConvergenceOverrejectAudit $rowConvergenceOverrejectAudit `
                -ConvergenceOverrejectDetail $rowConvergenceOverrejectDetail `
                -ConvergenceShapeAudit $rowConvergenceShapeAudit `
                -ConvergenceKernelResponseAudit $rowConvergenceKernelResponseAudit `
                -ChannelMissedSemantics $rowChannelMissedSemantics `
                -AcceptedMissDiagnostics $rowAcceptedMissDiagnostics `
                -AcceptedMissAbsoluteFloor $rowAcceptedMissAbsoluteFloor `
                -K3AcceptedRhoCap $rowK3AcceptedRhoCap `
                -MacDistancePolicy $rowMacDistancePolicy `
                -MacSurfaceMarginFactor $rowMacSurfaceMarginFactor `
                -MacSurfaceMarginPolicy $rowMacSurfaceMarginPolicy `
                -MacSurfaceMarginAspectThreshold $rowMacSurfaceMarginAspectThreshold `
                -MacSurfaceMarginRhoTrigger $rowMacSurfaceMarginRhoTrigger `
                -MacSurfaceMarginSoftAspectThreshold $rowMacSurfaceMarginSoftAspectThreshold `
                -MacSurfaceMarginSoftRhoTrigger $rowMacSurfaceMarginSoftRhoTrigger `
                -MacSurfaceMarginSoftFactor $rowMacSurfaceMarginSoftFactor `
                -ExpansionAuditChannel $rowExpansionAuditChannel `
                -LeafPanelKernel $rowLeafPanelKernel `
                -ClusterPolicy $rowClusterPolicy `
                -LeafPairKernel $rowLeafPairKernel `
                -LeafKernelProfile $rowLeafKernelProfile `
                -LeafWorkProfile $rowLeafWorkProfile `
                -MatvecCostProfile $rowMatvecCostProfile `
                -JcpDirectKernel $rowJcpDirectKernel `
                -TargetBatchSize $rowTargetBatchSize `
                -LeafPairMinSourcePanels $rowLeafPairMinSourcePanels `
                -LeafPairMinTargets $rowLeafPairMinTargets `
                -LeafPairSelectionProfile $rowLeafPairSelectionProfile `
                -LeafPairProbeCache $rowLeafPairProbeCache `
                -LeafPairTargetGroup $rowLeafPairTargetGroup `
                -LeafPairMaxLookahead $rowLeafPairMaxLookahead `
                -LeafPairDirectnessProbe $rowLeafPairDirectnessProbe `
                -LeafPairDirectnessCache $rowLeafPairDirectnessCache `
                -LeafPairDirectnessVerify $rowLeafPairDirectnessVerify `
                -LeafPairScheduler $rowLeafPairScheduler `
                -LeafPairWorklistMaxActive $rowLeafPairWorklistMaxActive `
                -LeafPairWorklistVerify $rowLeafPairWorklistVerify `
                -LeafBlockVerify $rowLeafBlockVerify `
                -TreeLeafSizeOverride $rowTreeLeafSize `
                -Preconditioner $rowPreconditioner `
                -PreconditionerBlockNear $rowPreconditionerBlockNear `
                -PreconditionerBlockPanelCap $rowPreconditionerBlockPanelCap `
                -PreconditionerMemoryBudgetMb $rowPreconditionerMemoryBudgetMb `
                -PreconditionerCoarseAggregateCap $rowPreconditionerCoarseAggregateCap `
                -PreconditionerCoarseBase $rowPreconditionerCoarseBase `
                -PreconditionerCoarseVerify $rowPreconditionerCoarseVerify `
                -GmresRestart $rowGmresRestart `
                -TrueResidualProbeInterval $rowTrueResidualProbeInterval
            $info = Invoke-TitanWithTimeout -ArgPath $arg -RunDir $runDir
            $rows += New-ResultRow -Fixture $fixtureName -CaseId $caseId `
                -Mode 'treecode' -TreeBuilder $tree.TreeBuilder `
                -MacSafety $tree.MacSafety `
                -ExpansionCenter $tree.ExpansionCenter `
                -FragmentationSafety $tree.FragmentationSafety `
                -MomentAudit $tree.MomentAudit `
                -MomentConstruction $tree.MomentConstruction `
                -MomentSum $tree.MomentSum `
                -ConvergenceAudit $tree.ConvergenceAudit `
                -ExpansionAudit $tree.ExpansionAudit `
                -TaylorOrder $rowTaylorOrder -Theta $rowTheta `
                -ConvergenceSafety $rowConvergenceSafety `
                -ConvergenceRhoCap $rowConvergenceRhoCap `
                -MacAspectCap $rowMacAspectCap `
                -ConvergenceEstimator $rowConvergenceEstimator `
                -ConvergenceTailThreshold $rowConvergenceTailThreshold `
                -ConvergenceEstimatorErrorCap $rowConvergenceEstimatorErrorCap `
                -ConvergenceUnstableRhoCap $rowConvergenceUnstableRhoCap `
                -ConvergenceCalibration $rowConvergenceCalibration `
                -ConvergenceDielectricAmp $rowConvergenceDielectricAmp `
                -ConvergenceEstimatorModel $rowConvergenceEstimatorModel `
                -ConvergenceChannelCalibration $rowConvergenceChannelCalibration `
                -ConvergenceSplitPolicy $rowConvergenceSplitPolicy `
                -ConvergenceSplitSelectivity $rowConvergenceSplitSelectivity `
                -ConvergenceSplitAbsoluteFloor $rowConvergenceSplitAbsoluteFloor `
                -ConvergenceOverrejectAudit $rowConvergenceOverrejectAudit `
                -ConvergenceOverrejectDetail $rowConvergenceOverrejectDetail `
                -ConvergenceShapeAudit $rowConvergenceShapeAudit `
                -ConvergenceKernelResponseAudit $rowConvergenceKernelResponseAudit `
                -ChannelMissedSemantics $rowChannelMissedSemantics `
                -AcceptedMissDiagnostics $rowAcceptedMissDiagnostics `
                -AcceptedMissAbsoluteFloor $rowAcceptedMissAbsoluteFloor `
                -K3AcceptedRhoCap $rowK3AcceptedRhoCap `
                -MacDistancePolicy $rowMacDistancePolicy `
                -MacSurfaceMarginFactor $rowMacSurfaceMarginFactor `
                -MacSurfaceMarginPolicy $rowMacSurfaceMarginPolicy `
                -MacSurfaceMarginAspectThreshold $rowMacSurfaceMarginAspectThreshold `
                -MacSurfaceMarginRhoTrigger $rowMacSurfaceMarginRhoTrigger `
                -MacSurfaceMarginSoftAspectThreshold $rowMacSurfaceMarginSoftAspectThreshold `
                -MacSurfaceMarginSoftRhoTrigger $rowMacSurfaceMarginSoftRhoTrigger `
                -MacSurfaceMarginSoftFactor $rowMacSurfaceMarginSoftFactor `
                -ExpansionAuditChannel $rowExpansionAuditChannel `
                -LeafPanelKernel $rowLeafPanelKernel `
                -ClusterPolicy $rowClusterPolicy `
                -LeafPairKernel $rowLeafPairKernel `
                -LeafKernelProfile $rowLeafKernelProfile `
                -LeafWorkProfile $rowLeafWorkProfile `
                -MatvecCostProfile $rowMatvecCostProfile `
                -JcpDirectKernel $rowJcpDirectKernel `
                -TargetBatchSize $rowTargetBatchSize `
                -LeafPairMinSourcePanels $rowLeafPairMinSourcePanels `
                -LeafPairMinTargets $rowLeafPairMinTargets `
                -LeafPairSelectionProfile $rowLeafPairSelectionProfile `
                -LeafPairProbeCache $rowLeafPairProbeCache `
                -LeafPairTargetGroup $rowLeafPairTargetGroup `
                -LeafPairMaxLookahead $rowLeafPairMaxLookahead `
                -LeafPairDirectnessProbe $rowLeafPairDirectnessProbe `
                -LeafPairDirectnessCache $rowLeafPairDirectnessCache `
                -LeafPairDirectnessVerify $rowLeafPairDirectnessVerify `
                -LeafPairScheduler $rowLeafPairScheduler `
                -LeafPairWorklistMaxActive $rowLeafPairWorklistMaxActive `
                -LeafPairWorklistVerify $rowLeafPairWorklistVerify `
                -LeafBlockVerify $rowLeafBlockVerify `
                -TreeLeafSizeOverride $rowTreeLeafSize `
                -Preconditioner $rowPreconditioner `
                -PreconditionerBlockNear $rowPreconditionerBlockNear `
                -PreconditionerBlockPanelCap $rowPreconditionerBlockPanelCap `
                -PreconditionerMemoryBudgetMb $rowPreconditionerMemoryBudgetMb `
                -PreconditionerCoarseAggregateCap $rowPreconditionerCoarseAggregateCap `
                -PreconditionerCoarseBase $rowPreconditionerCoarseBase `
                -PreconditionerCoarseVerify $rowPreconditionerCoarseVerify `
                -GmresRestart $rowGmresRestart `
                -TrueResidualProbeInterval $rowTrueResidualProbeInterval `
                -RunDir $runDir -RunInfo $info `
                -DirectEnergy $directEnergy
        }
    }
} finally {
    $env:TITAN_NO_PAUSE = $oldNoPause
}

Update-Classifications $rows
$outCsv = Join-Path $WorkRoot 'pbe_tabi_topdown_physical_gap_matrix.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII

$failed = @($rows | Where-Object { $_.status -ne 'ok' })
if ($failed.Count -gt 0 -and -not $AllowFailedStatus) {
    $messages = ($failed | ForEach-Object { "$($_.fixture)/$($_.case_id)=$($_.status)" }) -join '; '
    throw "PBE TABI topdown physical-gap matrix failed: $messages"
}

Write-Host "PBE TABI topdown physical-gap/fragmentation matrix COMPLETE: $outCsv"
