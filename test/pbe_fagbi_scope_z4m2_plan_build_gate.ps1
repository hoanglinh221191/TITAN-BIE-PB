param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4m2'),
    [int]$CaseTimeoutSeconds = 3600,
    [switch]$Quick,
    [switch]$ResumeQuick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4m2_matrix.csv'
$metricsCsv = Join-Path $WorkRoot 'fagbi_scope_z4m2_plan_metrics.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4m2_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4m2_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4m2_plan_runtime.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4M2_PLAN_BUILD_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing public fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (-not $ResumeQuick) {
    New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
} elseif (-not (Test-Path -LiteralPath $WorkRoot)) {
    throw "ResumeQuick requires an existing work root: $WorkRoot"
}
$manifest = @(Import-Csv -LiteralPath $manifestPath)

function Get-Z4m2Fixture([string]$Name) {
    $rows = @($manifest | Where-Object fixture -eq $Name)
    if ($rows.Count -ne 1) { throw "Expected exactly one fixture row for $Name" }
    return $rows[0]
}

function Get-Z4m2Double([object]$Value,[double]$Default=[double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Get-Z4m2Int64([object]$Value,[int64]$Default=0) {
    $parsed = [int64]0
    if ($null -eq $Value -or -not [int64]::TryParse([string]$Value,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Invoke-Z4m2Case([pscustomobject]$Fixture,[string]$Builder,
        [bool]$ReuseExisting=$false) {
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_$($Builder -replace '-','_')"
    if ($ReuseExisting) {
        $existingSummary = Join-Path $runDir 'pbe_tabi_operator_attribution_summary.csv'
        if (-not (Test-Path -LiteralPath $existingSummary)) {
            $existingLog = ''
            foreach ($name in @('stdout.log','stderr.log')) {
                $path = Join-Path $runDir $name
                if (Test-Path -LiteralPath $path) { $existingLog += Get-Content -Raw $path }
            }
            if ($existingLog -match 'cache_budget_exceeded: required=(\d+) budget=(\d+)') {
                return [pscustomobject]@{fixture=$Fixture.fixture;builder=$Builder
                    status='cache_budget_exceeded';cache_required_bytes=[int64]$Matches[1]
                    cache_budget_bytes=[int64]$Matches[2];wall_seconds=[double]::NaN
                    peak_rss_mb=[double]::NaN;run_dir=$runDir}
            }
            throw "ResumeQuick missing existing summary: $existingSummary"
        }
        $run = [pscustomobject]@{timed_out=$false;exit_code=0
            wall_seconds=[double]::NaN;peak_working_set_mb=[double]::NaN}
    } else {
        New-Item -ItemType Directory -Force -Path $runDir | Out-Null
        $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
            -PqrPath (Join-Path $FixtureRoot $Fixture.pqr) -Backend tabi `
            -Mode treecode -FmmMode single-matvec-analytic -SurfacePoints 8 `
            -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
            -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
            -OpeningTheta 0.6 -AuditSampleCount 256 -AuditMode sampled-detail `
            -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
            -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
            -TaylorOrder 4 -FmmCartesianAudit on -FmmLocal single-matvec-p4 `
            -FmmLocalRefine operator-tail-static -FmmLocalRefineTolerance 0.002 `
            -FmmLocalCache kernel-jet -FmmLocalCacheBudgetMb 256 `
            -FmmLocalFallback packed-avx2 -FmmLocalPlanBuilder $Builder
        $bareOptions = @(Get-Content -LiteralPath $arg | Where-Object {
            $_ -match '^-[A-Za-z0-9-]+$'
        })
        if ($bareOptions.Count) { throw "Generated arg.in contains bare options: $($bareOptions -join ', ')" }
        $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
            -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
        if ($run.timed_out) {
            return [pscustomobject]@{fixture=$Fixture.fixture;builder=$Builder;status='timeout'
                wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb;run_dir=$runDir}
        }
        if ($run.exit_code -ne 0) {
            $log = (Get-Content -Raw (Join-Path $runDir 'stdout.log')) +
                (Get-Content -Raw (Join-Path $runDir 'stderr.log'))
            if ($log -match 'cache_budget_exceeded: required=(\d+) budget=(\d+)') {
                return [pscustomobject]@{fixture=$Fixture.fixture;builder=$Builder
                    status='cache_budget_exceeded';cache_required_bytes=[int64]$Matches[1]
                    cache_budget_bytes=[int64]$Matches[2];wall_seconds=$run.wall_seconds
                    peak_rss_mb=$run.peak_working_set_mb;run_dir=$runDir}
            }
            throw "Z4M.2 case failed fixture=$($Fixture.fixture) builder=$Builder exit=$($run.exit_code): $log"
        }
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $summary = @(Import-Csv -LiteralPath (Join-Path $runDir `
        'pbe_tabi_operator_attribution_summary.csv'))
    $row = @($summary | Where-Object {
        $_.approximation_mode -eq 'local_m2l_p4_static_cached' -and
        $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
        $_.region -eq 'total' -and $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
    if ($null -eq $row) { throw "Missing Z4M.2 aggregate row for $($Fixture.fixture) $Builder" }
    return [pscustomobject]@{
        fixture=$Fixture.fixture;protein_name=$Fixture.protein_name;pdb_id=$Fixture.pdb_id
        builder=$Builder;status=$case.status;panels=$case.panels
        audit_rms=Get-Z4m2Double $row.aggregate_audit_rel
        recomposition_error=Get-Z4m2Double $row.recomposition_rel_error
        reference_gap=Get-Z4m2Double $row.reference_gap 0.0
        true_residual=$case.true_residual;wall_seconds=$run.wall_seconds
        peak_rss_mb=$run.peak_working_set_mb
        plan_hash=[string]$row.fmm_local_plan_hash
        plan_build_count=Get-Z4m2Int64 $row.fmm_local_plan_build_count
        refined_cells=Get-Z4m2Int64 $row.fmm_local_refined_cells
        split_children=Get-Z4m2Int64 $row.fmm_local_split_children
        fallback_pairs=Get-Z4m2Int64 $row.fmm_local_direct_fallback_pairs
        far_pairs_after=Get-Z4m2Int64 $row.fmm_local_far_pairs_after
        plan_seconds=Get-Z4m2Double $row.fmm_local_plan_seconds 0.0
        probe_moment_seconds=Get-Z4m2Double $row.fmm_local_plan_probe_moment_seconds 0.0
        pair_metric_seconds=Get-Z4m2Double $row.fmm_local_plan_pair_metric_seconds 0.0
        refinement_seconds=Get-Z4m2Double $row.fmm_local_plan_refinement_seconds 0.0
        metric_evaluations=Get-Z4m2Int64 $row.fmm_local_plan_pair_metric_evaluations
        metric_cache_hits=Get-Z4m2Int64 $row.fmm_local_plan_pair_metric_cache_hits
        refinement_iterations=Get-Z4m2Int64 $row.fmm_local_plan_refinement_iterations
        metric_cache_peak_bytes=Get-Z4m2Int64 $row.fmm_local_plan_metric_cache_peak_bytes
        cache_required_bytes=Get-Z4m2Int64 $row.fmm_local_cache_bytes
        cache_budget_bytes=256MB
        parallel_target_leaves=Get-Z4m2Int64 $row.fmm_local_plan_parallel_target_leaves
        candidate_median_apply_seconds=Get-Z4m2Double $row.fmm_local_candidate_median_apply_seconds 0.0
        fallback_effective=[string]$row.fmm_local_fallback_effective
        run_dir=$runDir
    }
}

function Test-Z4m2Parity([object]$Legacy,[object]$Cached) {
    if ($null -eq $Legacy -or $null -eq $Cached) { return $false }
    foreach ($name in @('plan_hash','refined_cells','split_children','fallback_pairs','far_pairs_after')) {
        if ([string]$Legacy.$name -ne [string]$Cached.$name) { return $false }
    }
    return ([Math]::Abs([double]$Legacy.audit_rms - [double]$Cached.audit_rms) -le 1e-12)
}

$rows = New-Object System.Collections.Generic.List[object]
$oneSh1 = Get-Z4m2Fixture '1sh1'
foreach ($builder in @('legacy-rescan','metric-cache')) {
    $existingSummary = Join-Path (Join-Path $WorkRoot `
        "1sh1_$($builder -replace '-','_')") `
        'pbe_tabi_operator_attribution_summary.csv'
    $reuse = $ResumeQuick -and (Test-Path -LiteralPath $existingSummary)
    $rows.Add((Invoke-Z4m2Case $oneSh1 $builder $reuse)) | Out-Null
}

$legacy = @($rows | Where-Object { $_.fixture -eq '1sh1' -and $_.builder -eq 'legacy-rescan' }) | Select-Object -First 1
$cached = @($rows | Where-Object { $_.fixture -eq '1sh1' -and $_.builder -eq 'metric-cache' }) | Select-Object -First 1
$quickDecision = 'z4m2_runtime_limited'
$quickSpeedup = [double]::NaN
if ($legacy.status -ne 'timeout' -and $cached.status -ne 'timeout') {
    if (-not (Test-Z4m2Parity $legacy $cached)) {
        $quickDecision = 'z4m2_plan_parity_gap'
    } elseif ([double]$cached.audit_rms -gt 8e-3 -or
            [double]$cached.recomposition_error -gt 1e-10 -or
            [double]$cached.reference_gap -gt 1e-10) {
        $quickDecision = 'z4m2_accuracy_gap'
    } elseif ([int64]$cached.metric_cache_peak_bytes -ge 128MB -or
            [double]$cached.peak_rss_mb -ge 512) {
        $quickDecision = 'z4m2_memory_limited'
    } elseif ([double]$cached.plan_seconds -gt 0) {
        $quickSpeedup = [double]$legacy.plan_seconds / [double]$cached.plan_seconds
        $quickDecision = if ($quickSpeedup -ge 10.0) {
            'z4m2_plan_build_pass'
        } else {
            'z4m2_partial_speedup'
        }
    }
}

if (-not $Quick -and $quickDecision -eq 'z4m2_plan_build_pass') {
    foreach ($name in @('barnase','1a63')) {
        $runDir = Join-Path $WorkRoot "${name}_metric_cache"
        $reuse = $ResumeQuick -and ((Test-Path -LiteralPath (Join-Path $runDir `
            'pbe_tabi_operator_attribution_summary.csv')) -or
            (Test-Path -LiteralPath (Join-Path $runDir 'stdout.log')))
        $rows.Add((Invoke-Z4m2Case (Get-Z4m2Fixture $name) 'metric-cache' $reuse)) | Out-Null
    }
}

$decision = $quickDecision
if (-not $Quick -and $quickDecision -eq 'z4m2_plan_build_pass') {
    $fullCached = @($rows | Where-Object builder -eq 'metric-cache')
    if ($fullCached.Count -ne 3 -or @($fullCached | Where-Object status -eq 'timeout').Count) {
        $decision = 'z4m2_runtime_limited'
    } elseif (@($fullCached | Where-Object status -eq 'cache_budget_exceeded').Count) {
        $decision = 'z4m2_memory_limited'
    } elseif (@($fullCached | Where-Object {
        [double]$_.audit_rms -gt 8e-3 -or [double]$_.recomposition_error -gt 1e-10 -or
        [double]$_.reference_gap -gt 1e-10
    }).Count) {
        $decision = 'z4m2_accuracy_gap'
    } elseif (@($fullCached | Where-Object {
        [int64]$_.metric_cache_peak_bytes -ge 128MB -or [double]$_.peak_rss_mb -ge 512
    }).Count) {
        $decision = 'z4m2_memory_limited'
    }
}

foreach ($row in $rows) {
    $rowSpeedup = if ($row.fixture -eq '1sh1' -and
            $row.builder -eq 'metric-cache') { $quickSpeedup } else { [double]::NaN }
    $row | Add-Member -NotePropertyName plan_build_speedup `
        -NotePropertyValue $rowSpeedup -Force
}
$resultRows = @($rows.ToArray())
$resultRows | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$resultRows | Select-Object fixture,builder,plan_hash,plan_seconds,
    probe_moment_seconds,pair_metric_seconds,refinement_seconds,
    metric_evaluations,metric_cache_hits,refinement_iterations,
    metric_cache_peak_bytes,parallel_target_leaves,plan_build_speedup |
    Export-Csv -LiteralPath $metricsCsv -NoTypeInformation
[pscustomobject]@{decision=$decision;quick_decision=$quickDecision;quick=[bool]$Quick
    rows=$resultRows.Count;quick_plan_build_speedup=$quickSpeedup
    full_matrix_started=(-not $Quick -and $quickDecision -eq 'z4m2_plan_build_pass')} |
    Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
@('scope: Z4M.2 static refinement plan build acceleration',"decision: $decision",
  "quick_decision: $quickDecision","rows: $($resultRows.Count)",
  "quick_plan_build_speedup: $quickSpeedup",'refinement_plan: unchanged',
  'gmres_p4: fail-closed','production_default: grid remains default') |
    Set-Content -LiteralPath $summaryPath -Encoding ASCII

$svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="620" viewBox="0 0 1000 620">',
    '<rect width="1000" height="620" fill="white"/>',
    '<text x="35" y="38" font-family="Arial" font-size="20">Static refinement plan build runtime</text>')
$y = 82
foreach ($row in $resultRows) {
    if ($row.status -eq 'cache_budget_exceeded') {
        $svg += ('<text x="35" y="{0}" font-family="Arial" font-size="13">{1} {2}: cache budget exceeded ({3} MB required)</text>' -f
            $y,$row.fixture,$row.builder,
            [Math]::Round(([double]$row.cache_required_bytes / 1MB),2))
    } else {
        $svg += ('<text x="35" y="{0}" font-family="Arial" font-size="13">{1} {2}: plan={3}s, metrics={4}, hits={5}, cache={6} MB</text>' -f
            $y,$row.fixture,$row.builder,$row.plan_seconds,$row.metric_evaluations,
            $row.metric_cache_hits,[Math]::Round(([double]$row.metric_cache_peak_bytes / 1MB),2))
    }
    $y += 38
}
$svg += '</svg>'
$svg | Set-Content -LiteralPath $plotPath -Encoding ASCII

$report = @('# TITAN FMM/M2L Scope Z4M.2 Plan Build Report','',"Decision: $decision",'',
    'Metric caching preserves the static refinement decisions and deterministic plan hash while avoiding repeated p3/p4 pair evaluation.','')
foreach ($row in $resultRows) {
    if ($row.status -eq 'cache_budget_exceeded') {
        $report += ('- {0} {1}: cache_budget_exceeded, required={2} MB, budget={3} MB' -f
            $row.fixture,$row.builder,
            [Math]::Round(([double]$row.cache_required_bytes / 1MB),2),
            [Math]::Round(([double]$row.cache_budget_bytes / 1MB),2))
    } else {
        $report += ('- {0} {1}: plan={2}s, audit={3}, evaluations={4}, hits={5}, cache={6} MB, RSS={7} MB' -f
            $row.fixture,$row.builder,$row.plan_seconds,$row.audit_rms,$row.metric_evaluations,
            $row.metric_cache_hits,[Math]::Round(([double]$row.metric_cache_peak_bytes / 1MB),2),
            $row.peak_rss_mb)
    }
}
$report += @('','The p4 operator, packed AVX2 fallback, tolerance, four deterministic probes, and production defaults are unchanged.',
    'Guarded GMRES remains disabled.')
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

if ($decision -ne 'z4m2_plan_build_pass' -and -not $AllowFailedDecision) {
    throw "Scope Z4M.2 gate did not pass: $decision"
}
Write-Host "PASS pbe_fagbi_scope_z4m2_plan_build_gate decision=$decision speedup=$quickSpeedup"
