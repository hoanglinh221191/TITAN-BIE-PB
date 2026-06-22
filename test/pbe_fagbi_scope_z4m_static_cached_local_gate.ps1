param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4m'),
    [int]$CaseTimeoutSeconds = 3600,
    [switch]$Quick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4m_matrix.csv'
$planCsv = Join-Path $WorkRoot 'fagbi_scope_z4m_plan.csv'
$cacheCsv = Join-Path $WorkRoot 'fagbi_scope_z4m_cache.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4m_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4m_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4m_accuracy_runtime.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4M_STATIC_CACHED_LOCAL_REPORT.txt'

if (-not (Test-Path $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path $manifestPath)) { throw "Missing public fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$manifest = @(Import-Csv $manifestPath)
$auditTargets = if ($Quick) { 128 } else { 256 }

function Get-Z4mFixture([string]$Name) {
    $rows = @($manifest | Where-Object fixture -eq $Name)
    if ($rows.Count -ne 1) { throw "Expected exactly one fixture row for $Name" }
    return $rows[0]
}

function Get-Z4mDouble([object]$Value,[double]$Default=[double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Invoke-Z4mCase([pscustomobject]$Fixture,[string]$Variant) {
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_$Variant"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $localMode = if ($Variant -eq 'static_cached') { 'single-matvec-p4' } else { 'audit-p4' }
    $refineMode = if ($Variant -eq 'dynamic') { 'p3-p4-tail' } else { 'operator-tail-static' }
    $cacheMode = if ($Variant -eq 'static_cached') { 'kernel-jet' } else { 'off' }
    $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
        -PqrPath (Join-Path $FixtureRoot $Fixture.pqr) -Backend tabi `
        -Mode treecode -FmmMode single-matvec-analytic -SurfacePoints 8 `
        -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount $auditTargets `
        -AuditMode sampled-detail -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular -TaylorOrder 4 `
        -FmmCartesianAudit on -FmmLocal $localMode `
        -FmmLocalRefine $refineMode -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache $cacheMode -FmmLocalCacheBudgetMb 256
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) {
        return [pscustomobject]@{fixture=$Fixture.fixture;variant=$Variant
            status='timeout';wall_seconds=$run.wall_seconds;run_dir=$runDir}
    }
    if ($run.exit_code -ne 0) {
        $log = (Get-Content -Raw (Join-Path $runDir 'stdout.log')) +
            (Get-Content -Raw (Join-Path $runDir 'stderr.log'))
        if ($log -match 'cache_budget_exceeded') {
            return [pscustomobject]@{fixture=$Fixture.fixture;variant=$Variant
                status='cache_budget_exceeded';wall_seconds=$run.wall_seconds;run_dir=$runDir}
        }
        throw "Z4M case failed fixture=$($Fixture.fixture) variant=$Variant exit=$($run.exit_code): $log"
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $summary = @(Import-Csv (Join-Path $runDir 'pbe_tabi_operator_attribution_summary.csv'))
    $mode = if ($Variant -eq 'dynamic') { 'local_m2l_p4' } `
        elseif ($Variant -eq 'static_cached') { 'local_m2l_p4_static_cached' } `
        else { 'local_m2l_p4_static' }
    $row = @($summary | Where-Object {
        $_.approximation_mode -eq $mode -and $_.probe -eq 'solution' -and
        $_.row_component -eq 'both' -and $_.region -eq 'total' -and
        $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
    if ($null -eq $row) { throw "Missing Z4M aggregate row mode=$mode" }
    return [pscustomobject]@{
        fixture=$Fixture.fixture;protein_name=$Fixture.protein_name;pdb_id=$Fixture.pdb_id
        variant=$Variant;status=$case.status;panels=$case.panels
        audit_rms=Get-Z4mDouble $row.aggregate_audit_rel
        recomposition_error=Get-Z4mDouble $row.recomposition_rel_error
        reference_gap=Get-Z4mDouble $row.reference_gap 0.0
        true_residual=$case.true_residual;wall_seconds=$run.wall_seconds
        peak_rss_mb=$run.peak_working_set_mb
        plan_hash=[string]$row.fmm_local_plan_hash
        plan_build_count=[int64]$row.fmm_local_plan_build_count
        static_probe_count=[int64]$row.fmm_local_static_probe_count
        refined_cells=[int64]$row.fmm_local_refined_cells
        fallback_pairs=[int64]$row.fmm_local_direct_fallback_pairs
        cache_mode=[string]$row.fmm_local_cache_mode
        cache_operators=[int64]$row.fmm_local_cache_operators
        cache_bytes=[int64]$row.fmm_local_cache_bytes
        cache_build_count=[int64]$row.fmm_local_cache_build_count
        cache_hits=[int64]$row.fmm_local_cache_hits
        plan_seconds=Get-Z4mDouble $row.fmm_local_plan_seconds 0.0
        cache_seconds=Get-Z4mDouble $row.fmm_local_cache_seconds 0.0
        candidate_apply_seconds=Get-Z4mDouble $row.fmm_local_candidate_apply_seconds 0.0
        candidate_repeat_count=[int64]$row.fmm_local_candidate_repeat_count
        candidate_median_apply_seconds=Get-Z4mDouble $row.fmm_local_candidate_median_apply_seconds 0.0
        direct_reference_seconds=Get-Z4mDouble $row.fmm_local_candidate_direct_reference_seconds 0.0
        run_dir=$runDir
    }
}

$matrix = New-Object System.Collections.Generic.List[object]
if ($Quick) {
    $fixture = Get-Z4mFixture '1sh1'
    foreach ($variant in @('static_uncached','static_cached')) {
        $matrix.Add((Invoke-Z4mCase $fixture $variant)) | Out-Null
    }
} else {
    $oneA63 = Get-Z4mFixture '1a63'
    foreach ($variant in @('dynamic','static_uncached','static_cached')) {
        $matrix.Add((Invoke-Z4mCase $oneA63 $variant)) | Out-Null
    }
    foreach ($name in @('barnase','1sh1')) {
        $matrix.Add((Invoke-Z4mCase (Get-Z4mFixture $name) 'static_cached')) | Out-Null
    }
}
$rows = $matrix.ToArray()
$rows | Export-Csv $matrixCsv -NoTypeInformation
$rows | Select-Object fixture,variant,plan_hash,plan_build_count,
    static_probe_count,refined_cells,fallback_pairs,plan_seconds |
    Export-Csv $planCsv -NoTypeInformation
$rows | Select-Object fixture,variant,cache_mode,cache_operators,cache_bytes,
    cache_build_count,cache_hits,cache_seconds,candidate_repeat_count,
    candidate_apply_seconds,candidate_median_apply_seconds,direct_reference_seconds |
    Export-Csv $cacheCsv -NoTypeInformation

$required = @($rows | Where-Object variant -eq 'static_cached')
$decision = 'z4m_accuracy_speed_pass'
if (@($rows | Where-Object status -eq 'timeout').Count) {
    $decision = 'z4m_runtime_limited'
} elseif (@($rows | Where-Object status -eq 'cache_budget_exceeded').Count) {
    $decision = 'z4m_cache_budget_exceeded'
} elseif (@($required | Where-Object {
        [double]$_.reference_gap -gt 1e-10 -or
        [double]$_.recomposition_error -gt 1e-10
    }).Count) {
    $decision = 'z4m_reference_gap'
} elseif (@($required | Where-Object {
        [double]$_.audit_rms -gt 8e-3
    }).Count) {
    $decision = 'z4m_static_plan_accuracy_gap'
} elseif (@($required | Where-Object {
        [double]$_.peak_rss_mb -ge 512 -or [int64]$_.cache_bytes -gt 256MB
    }).Count) {
    $decision = 'z4m_accuracy_pass_cache_limited'
} else {
    foreach ($cached in $required) {
        $uncached = @($rows | Where-Object {
            $_.fixture -eq $cached.fixture -and $_.variant -eq 'static_uncached'
        }) | Select-Object -First 1
        if ($null -ne $uncached -and [double]$cached.candidate_median_apply_seconds -gt 0) {
            $speedup = [double]$uncached.candidate_apply_seconds /
                [double]$cached.candidate_median_apply_seconds
            $cached | Add-Member -NotePropertyName speedup_vs_uncached `
                -NotePropertyValue $speedup -Force
            if ($speedup -lt 3.0) { $decision = 'z4m_accuracy_pass_cache_limited' }
        } else {
            $cached | Add-Member -NotePropertyName speedup_vs_uncached `
                -NotePropertyValue ([double]::NaN) -Force
        }
    }
}
if (-not $Quick -and $required.Count -ne 3) { $decision = 'z4m_runtime_limited' }

$rows | Export-Csv $matrixCsv -NoTypeInformation
$rows | Select-Object fixture,variant,cache_mode,cache_operators,cache_bytes,
    cache_build_count,cache_hits,cache_seconds,candidate_repeat_count,
    candidate_apply_seconds,candidate_median_apply_seconds,
    speedup_vs_uncached,direct_reference_seconds |
    Export-Csv $cacheCsv -NoTypeInformation

[pscustomobject]@{decision=$decision;quick=[bool]$Quick;rows=$rows.Count
    required_cached_rows=$required.Count;tolerance=0.002;cache_budget_mb=256} |
    Export-Csv $decisionCsv -NoTypeInformation
@('scope: Z4M static refinement and cached p4 local gate',"decision: $decision",
  "quick: $([bool]$Quick)","rows: $($rows.Count)",'theta: 0.6',
  'static_operator_probes: 4','candidate_repeats: 5','cache_budget_mb: 256',
  'production_degree1_fmm: unchanged','gmres_p4: fail-closed',
  'production_default: grid remains default') | Set-Content $summaryPath -Encoding ASCII

$svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="620" viewBox="0 0 1000 620">',
    '<rect width="1000" height="620" fill="white"/>',
    '<text x="35" y="38" font-family="Arial" font-size="20">Static p4 plan: accuracy and cached apply time</text>')
$y = 82
foreach ($row in $rows) {
    $svg += ('<text x="35" y="{0}" font-family="Arial" font-size="13">{1} {2}: audit={3}, apply={4}s, median={5}s, cache={6} MB</text>' -f
        $y,$row.fixture,$row.variant,$row.audit_rms,$row.candidate_apply_seconds,
        $row.candidate_median_apply_seconds,([Math]::Round([double]$row.cache_bytes/1MB,2)))
    $y += 38
}
$svg += '</svg>'
$svg | Set-Content $plotPath -Encoding ASCII

$report = @('# TITAN FMM/M2L Scope Z4M Static Cached Local Report','',
    "Decision: $decision",'',
    'Z4M builds a density-independent refinement plan from four deterministic operator probes and optionally caches compact order-8 kernel jets. The p4 candidate remains diagnostic-only.','')
foreach ($row in $rows) {
    $speedupText = if ($row.PSObject.Properties.Name -contains 'speedup_vs_uncached') {
        [string]$row.speedup_vs_uncached
    } else { 'n/a' }
    $report += ('- {0} {1}: audit={2}, plan_hash={3}, plan_builds={4}, cache={5} MB, hits={6}, median_apply={7}s, RSS={8} MB' -f
        $row.fixture,$row.variant,$row.audit_rms,$row.plan_hash,$row.plan_build_count,
        ([Math]::Round([double]$row.cache_bytes/1MB,2)),$row.cache_hits,
        $row.candidate_median_apply_seconds,$row.peak_rss_mb)
    if ($speedupText -ne 'n/a') { $report += "  speedup_vs_uncached=$speedupText" }
}
$report += @('','GMRES-p4 remains fail-closed.','Grid remains the production Poisson-Boltzmann default.')
$report | Set-Content $rootReport -Encoding ASCII

$bad = @(Get-ChildItem $WorkRoot -Recurse -Filter arg.in | ForEach-Object {
    $path = $_.FullName
    Get-Content $path | Where-Object {
        $_ -match '^\s*-' -and $_ -notmatch '^\s*-[^=\s]+='
    } | ForEach-Object { "$path :: $_" }
})
if ($bad.Count) { throw "Generated arg.in contains bare options: $($bad -join '; ')" }
if (-not $AllowFailedDecision -and $decision -ne 'z4m_accuracy_speed_pass') {
    throw "Scope Z4M decision: $decision"
}
Write-Host "PASS Scope Z4M decision=$decision rows=$($rows.Count)"
