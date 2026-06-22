param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4m3'),
    [int]$CaseTimeoutSeconds = 3600,
    [switch]$Quick,
    [switch]$ResumeQuick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4m3_matrix.csv'
$cacheCsv = Join-Path $WorkRoot 'fagbi_scope_z4m3_cache.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4m3_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4m3_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4m3_memory_runtime.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4M3_COMPRESSED_CACHE_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing public fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$existingRows = @()
if ($ResumeQuick) {
    if (-not (Test-Path -LiteralPath $matrixCsv)) {
        throw "ResumeQuick requires an existing matrix: $matrixCsv"
    }
    $existingRows = @(Import-Csv -LiteralPath $matrixCsv)
} else {
    New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
}
$manifest = @(Import-Csv -LiteralPath $manifestPath)

function Get-Z4m3Fixture([string]$Name) {
    $rows = @($manifest | Where-Object fixture -eq $Name)
    if ($rows.Count -ne 1) { throw "Expected exactly one fixture row for $Name" }
    return $rows[0]
}

function Get-Z4m3Double([object]$Value,[double]$Default=[double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Get-Z4m3Int64([object]$Value,[int64]$Default=0) {
    $parsed = [int64]0
    if ($null -eq $Value -or -not [int64]::TryParse([string]$Value,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Invoke-Z4m3Case([pscustomobject]$Fixture,[string]$CacheMode) {
    $label = switch ($CacheMode) {
        'off' { 'uncached' }
        'kernel-jet' { 'legacy_kernel_jet' }
        default { 'compressed_green_jet' }
    }
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_$label"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $localMode = if ($CacheMode -eq 'off') { 'audit-p4' } else { 'single-matvec-p4' }
    $fallbackMode = if ($CacheMode -eq 'off') { 'scalar' } else { 'packed-avx2' }
    $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
        -PqrPath (Join-Path $FixtureRoot $Fixture.pqr) -Backend tabi `
        -Mode treecode -FmmMode single-matvec-analytic -SurfacePoints 8 `
        -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount 256 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder 4 -FmmCartesianAudit on -FmmLocal $localMode `
        -FmmLocalRefine operator-tail-static -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache $CacheMode -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback $fallbackMode -FmmLocalPlanBuilder metric-cache
    $bareOptions = @(Get-Content -LiteralPath $arg | Where-Object {
        $_ -match '^-[A-Za-z0-9-]+$'
    })
    if ($bareOptions.Count) { throw "Generated arg.in contains bare options: $($bareOptions -join ', ')" }
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) {
        return [pscustomobject]@{fixture=$Fixture.fixture;cache_mode=$CacheMode
            cache_label=$label;status='timeout';wall_seconds=$run.wall_seconds
            peak_rss_mb=$run.peak_working_set_mb;run_dir=$runDir}
    }
    if ($run.exit_code -ne 0) {
        $log = (Get-Content -Raw (Join-Path $runDir 'stdout.log')) +
            (Get-Content -Raw (Join-Path $runDir 'stderr.log'))
        if ($log -match 'cache_budget_exceeded: required=(\d+) budget=(\d+)') {
            return [pscustomobject]@{fixture=$Fixture.fixture;cache_mode=$CacheMode
                cache_label=$label;status='cache_budget_exceeded'
                cache_bytes=[int64]$Matches[1];cache_budget_bytes=[int64]$Matches[2]
                wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb
                run_dir=$runDir}
        }
        throw "Z4M.3 case failed fixture=$($Fixture.fixture) cache=$CacheMode exit=$($run.exit_code): $log"
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $summary = @(Import-Csv -LiteralPath (Join-Path $runDir `
        'pbe_tabi_operator_attribution_summary.csv'))
    $expectedApproximation = switch ($CacheMode) {
        'off' { 'local_m2l_p4_static' }
        'kernel-jet' { 'local_m2l_p4_static_cached' }
        default { 'local_m2l_p4_static_green_cached' }
    }
    $row = @($summary | Where-Object {
        $_.approximation_mode -eq $expectedApproximation -and
        $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
        $_.region -eq 'total' -and $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
    if ($null -eq $row) { throw "Missing Z4M.3 aggregate row for $($Fixture.fixture) $CacheMode" }
    return [pscustomobject]@{
        fixture=$Fixture.fixture;protein_name=$Fixture.protein_name;pdb_id=$Fixture.pdb_id
        cache_mode=$CacheMode;cache_label=$label;status=$case.status;panels=$case.panels
        audit_rms=Get-Z4m3Double $row.aggregate_audit_rel
        recomposition_error=Get-Z4m3Double $row.recomposition_rel_error
        reference_gap=Get-Z4m3Double $row.reference_gap 0.0
        true_residual=$case.true_residual;wall_seconds=$run.wall_seconds
        peak_rss_mb=$run.peak_working_set_mb
        plan_hash=[string]$row.fmm_local_plan_hash
        refined_cells=Get-Z4m3Int64 $row.fmm_local_refined_cells
        split_children=Get-Z4m3Int64 $row.fmm_local_split_children
        fallback_pairs=Get-Z4m3Int64 $row.fmm_local_direct_fallback_pairs
        far_pairs=Get-Z4m3Int64 $row.fmm_local_far_pairs_after
        cache_layout=[string]$row.fmm_local_cache_layout
        cache_operator_bytes=Get-Z4m3Int64 $row.fmm_local_cache_operator_bytes
        cache_bytes=Get-Z4m3Int64 $row.fmm_local_cache_bytes
        uncompressed_estimated_bytes=Get-Z4m3Int64 `
            $row.fmm_local_cache_uncompressed_estimated_bytes
        compression_ratio=Get-Z4m3Double $row.fmm_local_cache_compression_ratio 1.0
        far_pair_offsets=Get-Z4m3Int64 $row.fmm_local_cache_far_pair_offsets
        cache_build_seconds=Get-Z4m3Double $row.fmm_local_cache_build_seconds 0.0
        cache_hits=Get-Z4m3Int64 $row.fmm_local_cache_hits
        candidate_median_apply_seconds=Get-Z4m3Double `
            $row.fmm_local_candidate_median_apply_seconds 0.0
        run_dir=$runDir
    }
}

$rows = New-Object System.Collections.Generic.List[object]
$oneSh1 = Get-Z4m3Fixture '1sh1'
foreach ($mode in @('off','kernel-jet','green-jet-compressed')) {
    $existing = @($existingRows | Where-Object {
        $_.fixture -eq '1sh1' -and $_.cache_mode -eq $mode
    }) | Select-Object -First 1
    $rows.Add($(if ($null -ne $existing) { $existing } else {
        Invoke-Z4m3Case $oneSh1 $mode
    })) | Out-Null
}
$existing1a63 = @($existingRows | Where-Object {
    $_.fixture -eq '1a63' -and $_.cache_mode -eq 'green-jet-compressed'
}) | Select-Object -First 1
$rows.Add($(if ($null -ne $existing1a63) { $existing1a63 } else {
    Invoke-Z4m3Case (Get-Z4m3Fixture '1a63') 'green-jet-compressed'
})) | Out-Null
if (-not $Quick) {
    $existingBarnase = @($existingRows | Where-Object {
        $_.fixture -eq 'barnase' -and
        $_.cache_mode -eq 'green-jet-compressed'
    }) | Select-Object -First 1
    $rows.Add($(if ($null -ne $existingBarnase) { $existingBarnase } else {
        Invoke-Z4m3Case (Get-Z4m3Fixture 'barnase') 'green-jet-compressed'
    })) | Out-Null
}

$resultRows = @($rows.ToArray())
$legacy = @($resultRows | Where-Object {
    $_.fixture -eq '1sh1' -and $_.cache_mode -eq 'kernel-jet'
}) | Select-Object -First 1
$compressed = @($resultRows | Where-Object {
    $_.fixture -eq '1sh1' -and $_.cache_mode -eq 'green-jet-compressed'
}) | Select-Object -First 1
$oneA63 = @($resultRows | Where-Object {
    $_.fixture -eq '1a63' -and $_.cache_mode -eq 'green-jet-compressed'
}) | Select-Object -First 1

$decision = 'z4m3_runtime_limited'
$required = @($resultRows | Where-Object cache_mode -eq 'green-jet-compressed')
if (@($resultRows | Where-Object status -eq 'timeout').Count) {
    $decision = 'z4m3_runtime_limited'
} elseif (@($resultRows | Where-Object status -eq 'cache_budget_exceeded').Count -or
          @($resultRows | Where-Object {
              [double]$_.peak_rss_mb -ge 512
          }).Count) {
    $decision = 'z4m3_memory_limited'
} elseif (@($required | Where-Object {
        [double]$_.audit_rms -gt 8e-3 -or
        [double]$_.recomposition_error -gt 1e-10 -or
        [double]$_.reference_gap -gt 1e-10
    }).Count) {
    $decision = 'z4m3_accuracy_gap'
} elseif ($compressed.plan_hash -ne $legacy.plan_hash -or
          $compressed.refined_cells -ne $legacy.refined_cells -or
          $compressed.split_children -ne $legacy.split_children -or
          $compressed.fallback_pairs -ne $legacy.fallback_pairs -or
          $compressed.far_pairs -ne $legacy.far_pairs) {
    $decision = 'z4m3_derivative_parity_gap'
} elseif ([double]$compressed.compression_ratio -lt 2.8 -or
          [double]$oneA63.cache_bytes -gt 160MB) {
    $decision = 'z4m3_compression_insufficient'
} elseif ([double]$legacy.candidate_median_apply_seconds -gt 0 -and
          [double]$compressed.candidate_median_apply_seconds -gt
              1.2 * [double]$legacy.candidate_median_apply_seconds) {
    $decision = 'z4m3_apply_runtime_regression'
} else {
    $decision = 'z4m3_compressed_cache_pass'
}

$resultRows | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$resultRows | Select-Object fixture,cache_mode,cache_layout,cache_operator_bytes,
    cache_bytes,uncompressed_estimated_bytes,compression_ratio,far_pair_offsets,
    cache_build_seconds,cache_hits,candidate_median_apply_seconds,peak_rss_mb |
    Export-Csv -LiteralPath $cacheCsv -NoTypeInformation
[pscustomobject]@{decision=$decision;quick=[bool]$Quick;rows=$resultRows.Count
    one_a63_cache_mb=if ($null -ne $oneA63.cache_bytes) {
        [Math]::Round(([double]$oneA63.cache_bytes / 1MB),3) } else { [double]::NaN }
    compression_ratio=$compressed.compression_ratio} |
    Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
@('scope: Z4M.3 lossless compressed Green-jet cache',"decision: $decision",
  "rows: $($resultRows.Count)","1a63_cache_bytes: $($oneA63.cache_bytes)",
  "compression_ratio: $($compressed.compression_ratio)",
  'precision: double','partial_cache: disabled','gmres_p4: fail-closed',
  'production_default: grid remains default') |
    Set-Content -LiteralPath $summaryPath -Encoding ASCII

$svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="620" viewBox="0 0 1050 620">',
    '<rect width="1050" height="620" fill="white"/>',
    '<text x="35" y="38" font-family="Arial" font-size="20">Lossless p4 cache memory and apply runtime</text>')
$y = 82
foreach ($row in $resultRows) {
    $cacheMb = if ($row.PSObject.Properties.Name -contains 'cache_bytes') {
        [Math]::Round(([double]$row.cache_bytes / 1MB),2) } else { [double]::NaN }
    $svg += ('<text x="35" y="{0}" font-family="Arial" font-size="13">{1} {2}: cache={3} MB, apply={4}s, RSS={5} MB, status={6}</text>' -f
        $y,$row.fixture,$row.cache_label,$cacheMb,
        $row.candidate_median_apply_seconds,$row.peak_rss_mb,$row.status)
    $y += 38
}
$svg += '</svg>'
$svg | Set-Content -LiteralPath $plotPath -Encoding ASCII

$report = @('# TITAN FMM/M2L Scope Z4M.3 Compressed Cache Report','',
    "Decision: $decision",'',
    'The compressed cache stores two double-precision degree-9 Green jets per accepted far pair. Source-normal coefficients are derived during contraction; the static refinement plan and packed fallback are unchanged.','')
foreach ($row in $resultRows) {
    $cacheMb = if ($row.PSObject.Properties.Name -contains 'cache_bytes') {
        [Math]::Round(([double]$row.cache_bytes / 1MB),2) } else { [double]::NaN }
    $report += ('- {0} {1}: status={2}, cache={3} MB, ratio={4}, audit={5}, apply={6}s, RSS={7} MB' -f
        $row.fixture,$row.cache_label,$row.status,$cacheMb,$row.compression_ratio,
        $row.audit_rms,$row.candidate_median_apply_seconds,$row.peak_rss_mb)
}
$report += @('','No FP32, quantization, partial cache, disk cache, or tiled recomputation is used.',
    'Guarded GMRES remains disabled and grid remains the production PB default.')
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

if ($decision -ne 'z4m3_compressed_cache_pass' -and -not $AllowFailedDecision) {
    throw "Scope Z4M.3 gate did not pass: $decision"
}
Write-Host "PASS pbe_fagbi_scope_z4m3_compressed_cache_gate decision=$decision"
