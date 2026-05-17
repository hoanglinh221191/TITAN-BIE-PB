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
    [string]$ConvergenceCalibration = 'stage3j'
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
        '-pbe-tabi-precond=diagonal',
        "-pbe-tabi-surface-points=$SurfacePoints",
        "-pbe-tabi-taylor-order=$TaylorOrder",
        "-pbe-tabi-tree-leaf-size=$TreeLeafSize",
        "-pbe-tabi-opening-theta=$thetaText",
        '-pbe-tabi-node-bounds=tight-aabb',
        '-pbe-tabi-interaction-mode=onthefly',
        '-pbe-tabi-near-integration=duffy_adaptive',
        '-pbe-tabi-quad-buffer=adaptive',
        '-pbe-tabi-gmres-restart=32',
        '-pbe-tabi-true-residual-probe=interval',
        '-pbe-tabi-true-residual-probe-scope=adaptive',
        '-pbe-tabi-true-residual-probe-interval=4',
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
        $argLines += "-pbe-tabi-expansion-audit=$ExpansionAudit"
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
        convergence_safety = Get-RowValue $timing 'tabi_convergence_safety' $ConvergenceSafety
        convergence_rho_cap = Get-RowValue $timing 'tabi_convergence_rho_cap' ($ConvergenceRhoCap.ToString('G17', $InvariantCulture))
        mac_aspect_cap = Get-RowValue $timing 'tabi_mac_aspect_cap' ($MacAspectCap.ToString('G17', $InvariantCulture))
        convergence_estimator = Get-RowValue $timing 'tabi_convergence_estimator' $ConvergenceEstimator
        convergence_tail_threshold = Get-RowValue $timing 'tabi_convergence_tail_threshold' ($ConvergenceTailThreshold.ToString('G17', $InvariantCulture))
        convergence_estimator_error_cap = Get-RowValue $timing 'tabi_convergence_estimator_error_cap' ($ConvergenceEstimatorErrorCap.ToString('G17', $InvariantCulture))
        convergence_unstable_rho_cap = Get-RowValue $timing 'tabi_convergence_unstable_rho_cap' ($ConvergenceUnstableRhoCap.ToString('G17', $InvariantCulture))
        convergence_calibration = Get-RowValue $timing 'tabi_convergence_calibration' $ConvergenceCalibration
        node_bounds = 'tight-aabb'
        interaction_mode = 'onthefly'
        near_integration = 'duffy_adaptive'
        quad_buffer = 'adaptive'
        preconditioner = 'diagonal'
        surface_points = $SurfacePoints
        tree_leaf_size = $TreeLeafSize
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
        far_geometric_unsafe_count = Get-RowValue $timing 'tabi_far_geometric_unsafe_count'
        far_geometric_unsafe_panels = Get-RowValue $timing 'tabi_far_geometric_unsafe_panels'
        far_geometric_gap_ratio_min = Get-RowValue $timing 'tabi_far_geometric_gap_ratio_min'
        far_geometric_gap_ratio_mean = Get-RowValue $timing 'tabi_far_geometric_gap_ratio_mean'
        mac_physical_gap_rejects = Get-RowValue $timing 'tabi_mac_physical_gap_rejects'
        mac_convergence_rho_rejects = Get-RowValue $timing 'tabi_mac_convergence_rho_rejects'
        mac_convergence_aspect_rejects = Get-RowValue $timing 'tabi_mac_convergence_aspect_rejects'
        mac_convergence_rho_reject_panels = Get-RowValue $timing 'tabi_mac_convergence_rho_reject_panels'
        mac_convergence_aspect_reject_panels = Get-RowValue $timing 'tabi_mac_convergence_aspect_reject_panels'
        convergence_estimator_samples = Get-RowValue $timing 'tabi_convergence_estimator_samples'
        convergence_estimator_unstable_far_count = Get-RowValue $timing 'tabi_convergence_estimator_unstable_far_count'
        convergence_estimator_unstable_far_panels = Get-RowValue $timing 'tabi_convergence_estimator_unstable_far_panels'
        convergence_estimator_split_rejects = Get-RowValue $timing 'tabi_convergence_estimator_split_rejects'
        convergence_estimator_split_reject_panels = Get-RowValue $timing 'tabi_convergence_estimator_split_reject_panels'
        convergence_estimator_tail_ratio_max = Get-RowValue $timing 'tabi_convergence_estimator_tail_ratio_max'
        convergence_estimator_tail_ratio_mean = Get-RowValue $timing 'tabi_convergence_estimator_tail_ratio_mean'
        convergence_estimator_error_max = Get-RowValue $timing 'tabi_convergence_estimator_error_max'
        convergence_estimator_error_mean = Get-RowValue $timing 'tabi_convergence_estimator_error_mean'
        convergence_tail_norm_p_max = Get-RowValue $timing 'tabi_convergence_tail_norm_p_max'
        convergence_tail_norm_p_mean = Get-RowValue $timing 'tabi_convergence_tail_norm_p_mean'
        convergence_tail_norm_pm1_max = Get-RowValue $timing 'tabi_convergence_tail_norm_pm1_max'
        convergence_tail_norm_pm1_mean = Get-RowValue $timing 'tabi_convergence_tail_norm_pm1_mean'
        convergence_effective_jcp_amp_max = Get-RowValue $timing 'tabi_convergence_effective_jcp_amp_max'
        convergence_effective_jcp_amp_mean = Get-RowValue $timing 'tabi_convergence_effective_jcp_amp_mean'
        convergence_estimator_high_error_caught = Get-RowValue $timing 'tabi_convergence_estimator_high_error_caught'
        convergence_estimator_high_error_missed = Get-RowValue $timing 'tabi_convergence_estimator_high_error_missed'
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
        far_node_interactions = Get-RowValue $timing 'tabi_far_node_interactions'
        matvec_calls = Get-RowValue $timing 'tabi_matvec_calls'
        rss_peak_bytes = Get-RowValue $memory 'peak_working_set_bytes'
        estimated_bytes_total = Get-RowValue $memory 'estimated_bytes_total'
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
        foreach ($row in $fixtureRows) {
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
            $missed = Convert-ToDoubleOrNull $row.convergence_estimator_high_error_missed
            $unsafe = Convert-ToDoubleOrNull $row.far_geometric_unsafe_count
            $rejects = Convert-ToDoubleOrNull $row.mac_physical_gap_rejects
            $matvecPass = $null -ne $matvecRel -and $matvecRel -le 2.5e-3
            $energyPass = $null -ne $driftRel -and $driftRel -le 5.0e-3
            $meanGateFlag = ([string]$row.far_expansion_mean_gate_pass) -eq '1'
            $meanGatePass = $meanGateFlag -and $null -ne $farMean -and $farMean -lt 5.0e-4 -and
                (($null -eq $missed) -or $missed -eq 0.0)
            $residualPass = $false
            if ($null -ne $sampledRel -and $null -ne $groupedSampled) {
                $residualPass = $sampledRel -le (2.0 * [math]::Max(1.0e-16, $groupedSampled))
            }
            if ($row.requested_convergence_estimator -ne 'off' -and $row.requested_convergence_safety -eq 'split' -and $matvecPass -and $energyPass -and $residualPass -and $meanGatePass) {
                $row.performance_classification = 'convergence_mean_gate_medium_pass'
            } elseif ($row.requested_convergence_estimator -ne 'off' -and $row.requested_convergence_safety -eq 'split' -and $matvecPass -and $energyPass -and $residualPass) {
                $row.performance_classification = 'convergence_estimator_medium_gate_pass'
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
                -ConvergenceCalibration $rowConvergenceCalibration
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
