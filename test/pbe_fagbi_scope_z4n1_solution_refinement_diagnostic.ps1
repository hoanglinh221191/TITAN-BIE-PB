param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4n1'),
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

$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4n1_matrix.csv'
$planCsv = Join-Path $WorkRoot 'fagbi_scope_z4n1_plan.csv'
$auditCsv = Join-Path $WorkRoot 'fagbi_scope_z4n1_audit.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4n1_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4n1_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4n1_refinement_plot.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4N1_SOLUTION_REFINEMENT_REPORT.txt'
$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$manifest = @(Import-Csv -LiteralPath $manifestPath)

function Get-Z4n1Fixture([string]$Name) {
    $row = @($manifest | Where-Object fixture -eq $Name)
    if ($row.Count -ne 1) { throw "Expected one public fixture row for $Name" }
    return $row[0]
}

function Get-Z4n1Double([object]$Value,[double]$Default=[double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Invoke-Z4n1Case(
    [pscustomobject]$Fixture,
    [string]$Condition,
    [int]$ProbeCount,
    [string]$ProbeMode
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
        -FmmLocalRefine operator-tail-static -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache green-jet-compressed -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback packed-avx2 -FmmLocalPlanBuilder metric-cache `
        -FmmLocalRefineProbeCount $ProbeCount `
        -FmmLocalRefineProbeMode $ProbeMode
    $bare = @(Get-Content -LiteralPath $arg | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) { throw "Generated arg.in contains bare options: $($bare -join ', ')" }
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) {
        return [pscustomobject]@{fixture=$Fixture.fixture;condition=$Condition
            probe_count=$ProbeCount;probe_mode=$ProbeMode;status='timeout'
            panels=0;iterations=0;final_residual=[double]::NaN
            true_residual=[double]::NaN;gmres_audit_rms=[double]::NaN
            candidate_audit_rms=[double]::NaN;candidate_relative_rms=[double]::NaN
            candidate_reference_gap=[double]::NaN;candidate_recomposition_error=[double]::NaN
            sampled_direct_true_residual=[double]::NaN;matvec_seconds_avg=[double]::NaN
            wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb
            plan_hash='';plan_build_count=0;cache_build_count=0
            far_pairs_after=0;direct_fallback_pairs=0;refinement_seconds=0.0
            solution_probe_used='false';operator_finite='false';run_dir=$runDir}
    }
    if ($run.exit_code -ne 0) {
        $log = (Get-Content -Raw (Join-Path $runDir 'stdout.log')) +
            (Get-Content -Raw (Join-Path $runDir 'stderr.log'))
        return [pscustomobject]@{fixture=$Fixture.fixture;condition=$Condition
            probe_count=$ProbeCount;probe_mode=$ProbeMode;status='failed'
            panels=0;iterations=0;final_residual=[double]::NaN
            true_residual=[double]::NaN;gmres_audit_rms=[double]::NaN
            candidate_audit_rms=[double]::NaN;candidate_relative_rms=[double]::NaN
            candidate_reference_gap=[double]::NaN;candidate_recomposition_error=[double]::NaN
            sampled_direct_true_residual=[double]::NaN;matvec_seconds_avg=[double]::NaN
            wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb
            plan_hash='';plan_build_count=0;cache_build_count=0
            far_pairs_after=0;direct_fallback_pairs=0;refinement_seconds=0.0
            solution_probe_used='false';operator_finite='false';failure_log=$log
            run_dir=$runDir}
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
    $summaryPath = Join-Path $runDir 'pbe_tabi_operator_attribution_summary.csv'
    $summary = @(if (Test-Path -LiteralPath $summaryPath) {
        @((Import-Csv -LiteralPath $summaryPath) | Where-Object {
            $_.approximation_mode -eq 'local_m2l_p4_static_green_cached' -and
            $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
            $_.region -eq 'total' -and $_.channel -eq 'aggregate'
        } | Select-Object -First 1)
    } else { @() })
    $summaryRow = if ($summary.Count) { $summary[0] } else { $null }
    return [pscustomobject]@{
        fixture=$Fixture.fixture;condition=$Condition
        probe_count=$ProbeCount;probe_mode=$ProbeMode;status=$case.status
        panels=$case.panels;iterations=$case.iterations
        final_residual=$case.final_residual;true_residual=$case.true_residual
        gmres_audit_rms=$case.tree_vs_direct_matvec_rel
        candidate_audit_rms=if ($summaryRow) { Convert-FmmCliDouble $summaryRow.aggregate_audit_rel } else { [double]::NaN }
        candidate_relative_rms=if ($summaryRow) { Convert-FmmCliDouble $summaryRow.relative_rms } else { [double]::NaN }
        candidate_reference_gap=if ($summaryRow) { Convert-FmmCliDouble $summaryRow.reference_gap 0.0 } else { [double]::NaN }
        candidate_recomposition_error=if ($summaryRow) { Convert-FmmCliDouble $summaryRow.recomposition_rel_error } else { [double]::NaN }
        sampled_direct_true_residual=Get-Z4n1Double `
            $profile.fmm_local_sampled_direct_true_residual 0.0
        matvec_seconds_avg=$case.matvec_seconds_avg
        wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb
        plan_hash=if ($summaryRow) { [string]$summaryRow.fmm_local_plan_hash } else { '' }
        plan_build_count=[int64]$profile.fmm_local_plan_build_count
        cache_build_count=[int64]$profile.fmm_local_cache_build_count
        far_pairs_after=if ($summaryRow) { Convert-FmmCliInt64 $summaryRow.fmm_local_far_pairs_after } else { 0 }
        direct_fallback_pairs=if ($summaryRow) { Convert-FmmCliInt64 $summaryRow.fmm_local_direct_fallback_pairs } else { 0 }
        refinement_seconds=if ($summaryRow) { Convert-FmmCliDouble $summaryRow.fmm_local_plan_refinement_seconds 0.0 } else { 0.0 }
        solution_probe_used=[string]$profile.fmm_local_solution_probe_used
        operator_finite=[string]$profile.fmm_local_operator_finite
        run_dir=$runDir
    }
}

$quickConditions = @(
    @{name='r4';count=4;mode='rademacher'},
    @{name='r8';count=8;mode='rademacher'},
    @{name='r16';count=16;mode='rademacher'},
    @{name='solution';count=4;mode='solution'},
    @{name='mixed';count=8;mode='mixed'}
)
$fixtureNames = if ($Quick) { @('1a63') } else { @('1a63') }
$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    $fixture = Get-Z4n1Fixture $fixtureName
    foreach ($condition in $quickConditions) {
        $rows.Add((Invoke-Z4n1Case $fixture $condition.name `
            $condition.count $condition.mode)) | Out-Null
    }
}

$matrix = @($rows.ToArray())
$baseline = @($matrix | Where-Object condition -eq 'r4' | Select-Object -First 1)
$best = @($matrix | Where-Object status -eq 'ok' | Sort-Object candidate_audit_rms | Select-Object -First 1)
$baselineAudit = if ($baseline.Count) { [double]$baseline[0].candidate_audit_rms } else { [double]::NaN }
$bestAudit = if ($best.Count) { [double]$best[0].candidate_audit_rms } else { [double]::NaN }
$improvement = if ((Test-FmmCliFinite $baselineAudit) -and $bestAudit -gt 0) {
    $baselineAudit / $bestAudit
} else { [double]::NaN }

$decision = 'z4n1_estimator_gap'
if (@($matrix | Where-Object status -eq 'timeout').Count) {
    $decision = 'z4n1_runtime_limited'
} elseif ($best.Count -and [double]$best[0].candidate_audit_rms -le 8e-3 -and
    [double]$best[0].candidate_recomposition_error -le 1e-10) {
    if ($best[0].probe_mode -eq 'rademacher' -and [int]$best[0].probe_count -gt 4) {
        $decision = 'z4n1_probe_count_sufficient'
    } elseif ($best[0].probe_mode -ne 'rademacher') {
        $decision = 'z4n1_solution_probe_required'
    } else {
        $decision = 'z4n1_accuracy_recovered'
    }
} elseif ((Test-FmmCliFinite $improvement) -and $improvement -lt 4.0) {
    $decision = 'z4n1_estimator_gap'
}

$matrix | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$matrix | Select-Object fixture,condition,probe_count,probe_mode,plan_hash,
    far_pairs_after,direct_fallback_pairs,refinement_seconds,solution_probe_used |
    Export-Csv -LiteralPath $planCsv -NoTypeInformation
$matrix | Select-Object fixture,condition,gmres_audit_rms,candidate_audit_rms,
    candidate_relative_rms,candidate_reference_gap,candidate_recomposition_error,
    sampled_direct_true_residual,true_residual,operator_finite |
    Export-Csv -LiteralPath $auditCsv -NoTypeInformation
[pscustomobject]@{decision=$decision;quick=[bool]$Quick;rows=$matrix.Count
    baseline_audit=$baselineAudit;best_audit=$bestAudit;improvement=$improvement} |
    Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
@('scope: Z4N.1 1A63 solution-subspace refinement diagnostic',
  "decision: $decision","quick: $([bool]$Quick)","rows: $($matrix.Count)",
    "baseline_candidate_audit: $baselineAudit","best_candidate_audit: $bestAudit",
  "improvement: $improvement",'production_default: grid',
  'gmres_p4: explicit diagnostic only') |
    Set-Content -LiteralPath $summaryPath -Encoding ASCII

$svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="960" height="520" viewBox="0 0 960 520">',
    '<rect width="960" height="520" fill="white"/>',
    '<text x="32" y="38" font-family="Arial" font-size="20">Z4N.1 1A63 refinement probe diagnostic</text>')
$y = 78
foreach ($row in $matrix) {
    $svg += ('<text x="32" y="{0}" font-family="Arial" font-size="13">{1}: mode={2}, probes={3}, candidate_audit={4}, gmres_true={5}, fallback_pairs={6}, RSS={7} MB</text>' -f
        $y,$row.condition,$row.probe_mode,$row.probe_count,$row.candidate_audit_rms,
        $row.sampled_direct_true_residual,$row.direct_fallback_pairs,
        $row.peak_rss_mb)
    $y += 32
}
$svg += '</svg>'
$svg | Set-Content -LiteralPath $plotPath -Encoding ASCII

$report = @('# TITAN FMM/M2L Scope Z4N.1 Solution-Subspace Refinement Report','',
    "Decision: $decision",'',
    'This diagnostic compares post-GMRES static refinement probe plans on full 1A63. Candidate audit values come from sampled-detail attribution on the converged vector; sampled direct true residual remains the baseline GMRES solve residual. It does not change production defaults or promote gmres-p4.','')
foreach ($row in $matrix) {
    $report += ('- {0}: mode={1}, probes={2}, status={3}, candidate_audit={4}, gmres_audit={5}, sampled_direct_residual={6}, plan_hash={7}, fallback_pairs={8}, wall={9}s, RSS={10} MB' -f
        $row.condition,$row.probe_mode,$row.probe_count,$row.status,
        $row.candidate_audit_rms,$row.gmres_audit_rms,
        $row.sampled_direct_true_residual,$row.plan_hash,
        $row.direct_fallback_pairs,$row.wall_seconds,$row.peak_rss_mb)
}
$report += @('','Grid remains the production PBE default. gmres-p4 remains explicit-only.')
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

if ($decision -eq 'z4n1_runtime_limited' -and -not $AllowFailedDecision) {
    throw "Scope Z4N.1 diagnostic did not complete: $decision"
}
Write-Host "PASS pbe_fagbi_scope_z4n1_solution_refinement_diagnostic decision=$decision"
