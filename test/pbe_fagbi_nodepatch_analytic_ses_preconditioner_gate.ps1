[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = 'test\tmp\fagbi_nodepatch_analytic_ses_preconditioner',
    [string[]]$Fixtures = @('barnase','1a63','1sh1'),
    [ValidateRange(0.25,4.0)][double]$SesGridSpacing = 1.25,
    [ValidateRange(30,86400)][int]$CaseTimeoutSeconds = 1800,
    [ValidateRange(256,8192)][int]$RasMemoryBudgetMb = 1024,
    [switch]$Quick,
    [switch]$Force,
    [switch]$RerunFailed,
    [string[]]$RerunStages = @(),
    [switch]$RecoverVerifiedProfiles,
    [switch]$IncludeNodePatchControl,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot))
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$analyticRoot = Join-Path $WorkRoot 'analytic_ses_ablation'
$nodePatchRoot = Join-Path $WorkRoot 'nodepatch_control'
$ablationScript = Join-Path $PSScriptRoot 'pbe_fagbi_graph_preconditioner_ablation.ps1'
$nodePatchScript = Join-Path $PSScriptRoot 'pbe_tabi2_nodepatch_analytic_ses_control.ps1'
$ablationReport = Join-Path $RepoRoot 'TITAN_FAGBI_ANALYTIC_SES_GRAPH_PRECONDITIONER_REPORT.txt'
$rootReport = Join-Path $RepoRoot 'TITAN_FAGBI_NODEPATCH_ANALYTIC_SES_PRECONDITIONER_REPORT.txt'

$ablationArgs = @{
    TitanExe = $TitanExe
    WorkRoot = $analyticRoot
    ReportPath = $ablationReport
    Fixtures = $Fixtures
    PrimarySurface = 'analytic-ses'
    SesGridSpacing = $SesGridSpacing
    CaseTimeoutSeconds = $CaseTimeoutSeconds
    RasMemoryBudgetMb = $RasMemoryBudgetMb
    AllowFailedDecision = $true
}
if ($Quick) { $ablationArgs.Quick = $true }
if ($Force) { $ablationArgs.Force = $true }
if ($RerunFailed) { $ablationArgs.RerunFailed = $true }
if ($RerunStages.Count -gt 0) { $ablationArgs.RerunStages = $RerunStages }
if ($RecoverVerifiedProfiles) { $ablationArgs.RecoverVerifiedProfiles = $true }
& $ablationScript @ablationArgs

if ($IncludeNodePatchControl) {
    & $nodePatchScript -TitanExe $TitanExe -WorkRoot $nodePatchRoot `
        -CaseTimeoutSeconds ([Math]::Min(1800,$CaseTimeoutSeconds))
}

$matrixPath = Join-Path $analyticRoot 'fagbi_graph_preconditioner_matrix.csv'
$nodePatchPath = Join-Path $nodePatchRoot 'nodepatch_analytic_ses_control.csv'
if (-not (Test-Path -LiteralPath $matrixPath)) {
    throw 'Missing analytic-SES ablation artifact.'
}
$matrix = @(Import-Csv -LiteralPath $matrixPath | Where-Object surface_protocol -eq 'analytic-ses')
$nodePatch = if ($IncludeNodePatchControl -and (Test-Path -LiteralPath $nodePatchPath)) {
    @(Import-Csv -LiteralPath $nodePatchPath)[-1]
} else { $null }
$fixtureNames = if ($Quick) { @('1sh1') } else {
    @($Fixtures | ForEach-Object { ([string]$_).ToLowerInvariant() })
}
$ranking = [Collections.Generic.List[object]]::new()
$allIluBest = $true
$allIluPass = $true
$runtimeLimited = $false
foreach ($fixture in $fixtureNames) {
    $rows = @($matrix | Where-Object fixture -eq $fixture)
    if (@($rows | Where-Object failure_class -eq 'timeout').Count -gt 0) {
        $runtimeLimited = $true
    }
    $ilu = @($rows | Where-Object preconditioner_label -in @(
        'distance_knn_ilu0','fmm_near_graph_ilu0'))
    $nonfactorized = @($rows | Where-Object preconditioner_label -notin @(
        'distance_knn_ilu0','fmm_near_graph_ilu0'))
    $iluPass = $ilu.Count -eq 2 -and @($ilu | Where-Object {
        $_.failure_class -ne 'pass' -or [double]$_.true_residual -gt 1e-4 -or
        [double]$_.audit_rms -gt 1e-2
    }).Count -eq 0
    $iluMin = if ($ilu.Count -gt 0) {
        ($ilu | Measure-Object -Property iterations -Minimum).Minimum
    } else { [double]::PositiveInfinity }
    $nonfactorizedMin = if ($nonfactorized.Count -gt 0) {
        ($nonfactorized | Measure-Object -Property iterations -Minimum).Minimum
    } else { [double]::PositiveInfinity }
    $iluBest = $iluPass -and [double]$iluMin -le [double]$nonfactorizedMin
    $allIluPass = $allIluPass -and $iluPass
    $allIluBest = $allIluBest -and $iluBest
    $ranking.Add([pscustomobject]@{
        fixture = $fixture
        panels = if ($rows.Count -gt 0) { [int64]$rows[0].panels } else { 0 }
        ilu_pass = $iluPass
        best_ilu_iterations = $iluMin
        best_nonfactorized_iterations = $nonfactorizedMin
        ilu_iteration_advantage = if ([double]::IsInfinity([double]$nonfactorizedMin)) {
            [double]::NaN
        } else { [double]$nonfactorizedMin - [double]$iluMin }
        ilu_best = $iluBest
    })
}
$decision = if ($runtimeLimited) {
    'analytic_ses_preconditioner_runtime_limited'
} elseif (-not $allIluPass) {
    'analytic_ses_ilu_convergence_or_accuracy_gap'
} elseif (-not $allIluBest) {
    'analytic_ses_ilu_not_best'
} elseif ($IncludeNodePatchControl -and
    ($null -eq $nodePatch -or [string]$nodePatch.contract_pass -ne 'True')) {
    'analytic_ses_ilu_best_nodepatch_control_gap'
} elseif ($IncludeNodePatchControl) {
    'analytic_ses_ilu_best_nodepatch_fmm_pending'
} else {
    'analytic_ses_ilu_best'
}
$rankingPath = Join-Path $WorkRoot 'fagbi_nodepatch_analytic_ses_ranking.csv'
$decisionPath = Join-Path $WorkRoot 'fagbi_nodepatch_analytic_ses_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_nodepatch_analytic_ses_summary.dat'
$ranking.ToArray() | Export-Csv -LiteralPath $rankingPath -NoTypeInformation -Encoding UTF8
@([pscustomobject]@{
    decision = $decision
    analytic_ses_ilu_pass = $allIluPass
    analytic_ses_ilu_best = $allIluBest
    nodepatch_direct_control_pass = if ($null -ne $nodePatch) { [string]$nodePatch.contract_pass } else { 'not_run' }
    nodepatch_ilu_comparison_available = $false
    nodepatch_limitation = if ($null -ne $nodePatch) {
        [string]$nodePatch.implementation_scope
    } else { 'not_run_in_analytic_ses_slice' }
    matrix = $matrixPath
}) | Export-Csv -LiteralPath $decisionPath -NoTypeInformation -Encoding UTF8
@(
    "decision $decision",
    "analytic_ses_ilu_pass $allIluPass",
    "analytic_ses_ilu_best $allIluBest",
    'nodepatch_ilu_comparison_available false',
    "matrix $matrixPath",
    "ranking $rankingPath"
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII
$report = @(
    'TITAN NODE-PATCH / ANALYTIC-SES PRECONDITIONER GATE',
    '==================================================',
    '',
    "Decision: $decision",
    "Analytic-SES spacing: $SesGridSpacing A",
    "Fixtures: $($fixtureNames -join ', ')",
    "Analytic-SES ILU rows pass: $allIluPass",
    "ILU has the lowest iteration count on every completed fixture: $allIluBest",
    '',
    'Node-patch scope:',
    'The current tabi2-nodepatch path is a direct report/control path. It does not yet implement a vertex-patch FMM operator or FMM-native ILU/RAS assembly. The node-patch control therefore verifies analytic-SES compatibility only; it is not included in the six-preconditioner ranking.',
    '',
    'Scientific interpretation:',
    'The analytic-SES matrix can test whether sparse ILU(0) remains the best algebraic family after changing the molecular surface. It cannot yet establish the same result for a genuine node-patch discretization. Doing so requires a mass-lumped vertex-patch operator and matching near-field block assembly rather than relabeling the centroid operator.',
    '',
    "Matrix: $matrixPath",
    "Ranking: $rankingPath",
    "Decision: $decisionPath",
    "Figure: $(Join-Path $analyticRoot 'fagbi_graph_preconditioner_convergence.svg')"
)
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII
Write-Host "Node-patch / analytic-SES decision: $decision"
Write-Host "Matrix: $matrixPath"
Write-Host "Figure: $(Join-Path $analyticRoot 'fagbi_graph_preconditioner_convergence.svg')"
if (-not $AllowFailedDecision -and $decision -notin @(
        'analytic_ses_ilu_best','analytic_ses_ilu_best_nodepatch_fmm_pending')) {
    throw "Node-patch / analytic-SES decision: $decision"
}
