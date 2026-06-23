param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4n2'),
    [int]$CaseTimeoutSeconds = 3600,
    [switch]$Quick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
if (-not [IO.Path]::IsPathRooted($FixtureRoot)) {
    $FixtureRoot = Join-Path $RepoRoot $FixtureRoot
}
$FixtureRoot = [IO.Path]::GetFullPath($FixtureRoot)
if (-not [IO.Path]::IsPathRooted($WorkRoot)) {
    $WorkRoot = Join-Path $RepoRoot $WorkRoot
}
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$RuntimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
        $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
    }
    $RuntimeRepoRoot = Split-Path -Parent $commonGitDir
}

$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4n2_matrix.csv'
$outlierCsv = Join-Path $WorkRoot 'fagbi_scope_z4n2_outliers.csv'
$geometryCsv = Join-Path $WorkRoot 'fagbi_scope_z4n2_geometry.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4n2_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4n2_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4n2_refinement_plot.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4N2_OUTLIER_REFINEMENT_REPORT.txt'
$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$manifest = @(Import-Csv -LiteralPath $manifestPath)

function Get-Z4n2Fixture([string]$Name) {
    $row = @($manifest | Where-Object fixture -eq $Name)
    if ($row.Count -ne 1) { throw "Expected one public fixture row for $Name" }
    return $row[0]
}

function Get-Z4n2Double([object]$Value,[double]$Default=[double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Invoke-Z4n2Case(
    [pscustomobject]$Fixture,
    [string]$Condition,
    [string]$RefineMode,
    [string]$GeometryGuard,
    [double]$GeometryTheta
) {
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_$Condition"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $pqr = Join-Path $FixtureRoot $Fixture.pqr
    $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RuntimeRepoRoot `
        -PqrPath $pqr -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -SurfacePoints 8 -MaxIter 720 `
        -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount 256 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder 4 -FmmCartesianAudit off -FmmLocal gmres-p4 `
        -FmmLocalRefine $RefineMode -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache green-jet-compressed -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback packed-avx2 -FmmLocalPlanBuilder metric-cache `
        -FmmLocalRefineProbeCount 4 -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard $GeometryGuard -FmmLocalGeometryTheta $GeometryTheta
    $bare = @(Get-Content -LiteralPath $arg | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) { throw "Generated arg.in contains bare options: $($bare -join ', ')" }
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) {
        return [pscustomobject]@{
            fixture=$Fixture.fixture;condition=$Condition;status='timeout'
            refine_mode=$RefineMode;geometry_guard=$GeometryGuard
            geometry_theta=$GeometryTheta;panels=0;iterations=0
            audit_rms=[double]::NaN;sampled_direct_true_residual=[double]::NaN
            linf_max_abs_tail=[double]::NaN;linf_max_relative_tail=[double]::NaN
            linf_outlier_panels=0;geometry_rejected_cells=0
            geometry_direct_fallback_pairs=0;geometry_gap_min=[double]::NaN
            geometry_ratio_max=[double]::NaN;direct_fallback_pairs=0
            far_pairs_after=0;wall_seconds=$run.wall_seconds
            peak_rss_mb=$run.peak_working_set_mb;run_dir=$runDir
        }
    }
    if ($run.exit_code -ne 0) {
        $log = (Get-Content -Raw (Join-Path $runDir 'stdout.log')) +
            (Get-Content -Raw (Join-Path $runDir 'stderr.log'))
        return [pscustomobject]@{
            fixture=$Fixture.fixture;condition=$Condition;status='failed'
            refine_mode=$RefineMode;geometry_guard=$GeometryGuard
            geometry_theta=$GeometryTheta;panels=0;iterations=0
            audit_rms=[double]::NaN;sampled_direct_true_residual=[double]::NaN
            linf_max_abs_tail=[double]::NaN;linf_max_relative_tail=[double]::NaN
            linf_outlier_panels=0;geometry_rejected_cells=0
            geometry_direct_fallback_pairs=0;geometry_gap_min=[double]::NaN
            geometry_ratio_max=[double]::NaN;direct_fallback_pairs=0
            far_pairs_after=0;wall_seconds=$run.wall_seconds
            peak_rss_mb=$run.peak_working_set_mb;failure_log=$log;run_dir=$runDir
        }
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
    $summaryFile = Join-Path $runDir 'pbe_tabi_operator_attribution_summary.csv'
    $summary = @(if (Test-Path -LiteralPath $summaryFile) {
        @((Import-Csv -LiteralPath $summaryFile) | Where-Object {
            $_.approximation_mode -eq 'local_m2l_p4_static_green_cached' -and
            $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
            $_.region -eq 'total' -and $_.channel -eq 'aggregate'
        } | Select-Object -First 1)
    } else { @() })
    $summaryRow = if ($summary.Count) { $summary[0] } else { $null }
    return [pscustomobject]@{
        fixture=$Fixture.fixture;condition=$Condition;status=$case.status
        refine_mode=$RefineMode;geometry_guard=$GeometryGuard
        geometry_theta=$GeometryTheta;panels=$case.panels
        iterations=$case.iterations;audit_rms=$case.tree_vs_direct_matvec_rel
        sampled_direct_true_residual=Get-Z4n2Double `
            $profile.fmm_local_sampled_direct_true_residual 0.0
        linf_max_abs_tail=Get-Z4n2Double $profile.fmm_local_linf_max_abs_tail 0.0
        linf_max_relative_tail=Get-Z4n2Double `
            $profile.fmm_local_linf_max_relative_tail 0.0
        linf_outlier_panels=[int64]$profile.fmm_local_linf_outlier_panels
        geometry_rejected_cells=[int64]$profile.fmm_local_geometry_guard_rejected_cells
        geometry_direct_fallback_pairs=[int64]$profile.fmm_local_geometry_guard_direct_fallback_pairs
        geometry_gap_min=Get-Z4n2Double $profile.fmm_local_geometry_gap_min 0.0
        geometry_ratio_max=Get-Z4n2Double $profile.fmm_local_geometry_ratio_max 0.0
        direct_fallback_pairs=if ($summaryRow) { Convert-FmmCliInt64 $summaryRow.fmm_local_direct_fallback_pairs } else { 0 }
        far_pairs_after=if ($summaryRow) { Convert-FmmCliInt64 $summaryRow.fmm_local_far_pairs_after } else { 0 }
        candidate_audit_rms=if ($summaryRow) { Convert-FmmCliDouble $summaryRow.aggregate_audit_rel } else { [double]::NaN }
        candidate_recomposition_error=if ($summaryRow) { Convert-FmmCliDouble $summaryRow.recomposition_rel_error } else { [double]::NaN }
        wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb
        run_dir=$runDir
    }
}

$conditions = @(
    @{name='baseline_rms';refine='operator-tail-static';guard='off';theta=0.0},
    @{name='linf_only';refine='operator-tail-linf';guard='off';theta=0.0},
    @{name='strict_geometry_only';refine='operator-tail-static';guard='strict-geometry-mac';theta=0.35},
    @{name='geometry_linf';refine='operator-tail-linf';guard='strict-geometry-mac';theta=0.35}
)

$fixtureNames = if ($Quick) { @('1a63') } else { @('1a63') }
$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    $fixture = Get-Z4n2Fixture $fixtureName
    foreach ($condition in $conditions) {
        $rows.Add((Invoke-Z4n2Case $fixture $condition.name `
            $condition.refine $condition.guard $condition.theta)) | Out-Null
    }
}

$matrix = @($rows.ToArray())
$baseline = @($matrix | Where-Object condition -eq 'baseline_rms' | Select-Object -First 1)
$best = @($matrix | Where-Object status -eq 'ok' |
    Sort-Object @{Expression='candidate_audit_rms';Ascending=$true} |
    Select-Object -First 1)
$baselineAudit = if ($baseline.Count) { [double]$baseline[0].candidate_audit_rms } else { [double]::NaN }
$bestAudit = if ($best.Count) { [double]$best[0].candidate_audit_rms } else { [double]::NaN }
$improvement = if ((Test-FmmCliFinite $baselineAudit) -and $bestAudit -gt 0) {
    $baselineAudit / $bestAudit
} else { [double]::NaN }
$fallbackRatio = if ($baseline.Count -and $best.Count -and [double]$baseline[0].direct_fallback_pairs -gt 0) {
    [double]$best[0].direct_fallback_pairs / [double]$baseline[0].direct_fallback_pairs
} else { [double]::NaN }

$decision = 'z4n2_estimator_gap'
if (@($matrix | Where-Object status -eq 'timeout').Count) {
    $decision = 'z4n2_runtime_limited'
} elseif (@($matrix | Where-Object status -eq 'failed').Count) {
    $decision = 'z4n2_infrastructure_gap'
} elseif ($best.Count -and
    ((([double]$best[0].candidate_audit_rms -le 8e-3) -and
      ([double]$best[0].sampled_direct_true_residual -le 8e-3)) -or
     ((Test-FmmCliFinite $improvement) -and $improvement -ge 4.0)) -and
    ((-not (Test-FmmCliFinite $fallbackRatio)) -or $fallbackRatio -le 3.0) -and
    [double]$best[0].peak_rss_mb -lt 512.0) {
    $decision = 'z4n2_outlier_refinement_promising'
} elseif ($best.Count -and (Test-FmmCliFinite $fallbackRatio) -and $fallbackRatio -gt 3.0) {
    $decision = 'z4n2_refinement_cost_excessive'
}

$matrix | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$matrix | Select-Object fixture,condition,refine_mode,linf_max_abs_tail,
    linf_max_relative_tail,linf_outlier_panels,direct_fallback_pairs,
    far_pairs_after,candidate_audit_rms,sampled_direct_true_residual |
    Export-Csv -LiteralPath $outlierCsv -NoTypeInformation
$matrix | Select-Object fixture,condition,geometry_guard,geometry_theta,
    geometry_rejected_cells,geometry_direct_fallback_pairs,geometry_gap_min,
    geometry_ratio_max,direct_fallback_pairs,candidate_audit_rms |
    Export-Csv -LiteralPath $geometryCsv -NoTypeInformation
[pscustomobject]@{decision=$decision;quick=[bool]$Quick;rows=$matrix.Count
    baseline_audit=$baselineAudit;best_condition=if ($best.Count) { $best[0].condition } else { '' }
    best_audit=$bestAudit;improvement=$improvement;fallback_ratio=$fallbackRatio} |
    Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
@('scope: Z4N.2 outlier refinement diagnostic',
  "decision: $decision","quick: $([bool]$Quick)","rows: $($matrix.Count)",
  "baseline_candidate_audit: $baselineAudit","best_candidate_audit: $bestAudit",
  "improvement: $improvement","fallback_ratio: $fallbackRatio",
  'production_default: grid','gmres_p4: explicit diagnostic only') |
    Set-Content -LiteralPath $summaryPath -Encoding ASCII

$svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="960" height="520" viewBox="0 0 960 520">',
    '<rect width="960" height="520" fill="white"/>',
    '<text x="32" y="38" font-family="Arial" font-size="20">Z4N.2 1A63 outlier refinement diagnostic</text>',
    '<line x1="90" y1="430" x2="900" y2="430" stroke="#222"/>',
    '<line x1="90" y1="80" x2="90" y2="430" stroke="#222"/>')
$i = 0
foreach ($row in $matrix) {
    $x = 140 + 180 * $i
    $audit = if (Test-FmmCliFinite ([double]$row.candidate_audit_rms)) {
        [double]$row.candidate_audit_rms
    } else { 1.0 }
    $bar = [Math]::Min(320.0, 70.0 * (-[Math]::Log10([Math]::Max(1e-6, $audit))))
    $svg += "<rect x='$($x)' y='$([int](430-$bar))' width='72' height='$([int]$bar)' fill='#3b82f6'/>"
    $svg += "<text x='$($x-18)' y='462' font-family='Arial' font-size='12'>$($row.condition)</text>"
    $svg += "<text x='$($x-8)' y='$([int](420-$bar))' font-family='Arial' font-size='12'>$('{0:G3}' -f $audit)</text>"
    ++$i
}
$svg += '<text x="96" y="76" font-family="Arial" font-size="12">-log10 audit RMS</text>'
$svg += '</svg>'
$svg | Set-Content -LiteralPath $plotPath -Encoding ASCII

@(
    'TITAN FMM/M2L Scope Z4N.2 Outlier Refinement Diagnostic',
    "Decision: $decision",
    "Rows: $($matrix.Count)",
    "Baseline candidate audit: $baselineAudit",
    "Best candidate audit: $bestAudit",
    "Improvement: $improvement",
    "Fallback ratio: $fallbackRatio",
    '',
    'Artifacts:',
    $matrixCsv,
    $outlierCsv,
    $geometryCsv,
    $decisionCsv,
    $summaryPath,
    $plotPath
) | Set-Content -LiteralPath $rootReport -Encoding ASCII

if (-not $AllowFailedDecision -and
    $decision -notin @('z4n2_outlier_refinement_promising','z4n2_estimator_gap',
                       'z4n2_refinement_cost_excessive')) {
    throw "Z4N.2 diagnostic failed with decision $decision"
}

Write-Host "Z4N.2 diagnostic decision: $decision"
