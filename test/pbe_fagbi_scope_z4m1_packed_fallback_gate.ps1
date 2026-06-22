param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4m1'),
    [int]$CaseTimeoutSeconds = 3600,
    [switch]$Quick,
    [switch]$ResumeQuick,
    [switch]$PackageOnly,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4m1_matrix.csv'
$worklistCsv = Join-Path $WorkRoot 'fagbi_scope_z4m1_worklist.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4m1_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4m1_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4m1_accuracy_runtime.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4M1_PACKED_FALLBACK_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing public fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if ($PackageOnly) {
    if (-not (Test-Path -LiteralPath $matrixCsv)) {
        throw "PackageOnly requires an existing matrix: $matrixCsv"
    }
} elseif (-not $ResumeQuick) {
    New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
} elseif (-not (Test-Path -LiteralPath $matrixCsv)) {
    throw "ResumeQuick requires an existing quick matrix: $matrixCsv"
}
$manifest = @(Import-Csv -LiteralPath $manifestPath)
$auditTargets = if ($Quick) { 128 } else { 256 }

function Get-Z4m1Fixture([string]$Name) {
    $rows = @($manifest | Where-Object fixture -eq $Name)
    if ($rows.Count -ne 1) { throw "Expected exactly one fixture row for $Name" }
    return $rows[0]
}

function Get-Z4m1Double([object]$Value,[double]$Default=[double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Invoke-Z4m1Case([pscustomobject]$Fixture,[string]$Variant) {
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_$Variant"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $cached = $Variant -ne 'static_uncached'
    $fallback = if ($Variant -eq 'packed_cached') { 'packed-avx2' } else { 'scalar' }
    $localMode = if ($cached) { 'single-matvec-p4' } else { 'audit-p4' }
    $cacheMode = if ($cached) { 'kernel-jet' } else { 'off' }
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
        -FmmLocalRefine operator-tail-static -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache $cacheMode -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback $fallback
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) {
        return [pscustomobject]@{fixture=$Fixture.fixture;variant=$Variant
            status='timeout';wall_seconds=$run.wall_seconds;run_dir=$runDir}
    }
    if ($run.exit_code -ne 0) {
        $log = (Get-Content -Raw (Join-Path $runDir 'stdout.log')) +
            (Get-Content -Raw (Join-Path $runDir 'stderr.log'))
        if ($log -match 'packed_avx2_unavailable') {
            return [pscustomobject]@{fixture=$Fixture.fixture;variant=$Variant
                status='avx2_unavailable';wall_seconds=$run.wall_seconds;run_dir=$runDir}
        }
        throw "Z4M.1 case failed fixture=$($Fixture.fixture) variant=$Variant exit=$($run.exit_code): $log"
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $summary = @(Import-Csv -LiteralPath (Join-Path $runDir `
        'pbe_tabi_operator_attribution_summary.csv'))
    $mode = if ($cached) { 'local_m2l_p4_static_cached' } else { 'local_m2l_p4_static' }
    $row = @($summary | Where-Object {
        $_.approximation_mode -eq $mode -and $_.probe -eq 'solution' -and
        $_.row_component -eq 'both' -and $_.region -eq 'total' -and
        $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
    if ($null -eq $row) { throw "Missing Z4M.1 aggregate row mode=$mode" }
    return [pscustomobject]@{
        fixture=$Fixture.fixture;protein_name=$Fixture.protein_name;pdb_id=$Fixture.pdb_id
        variant=$Variant;status=$case.status;panels=$case.panels
        audit_rms=Get-Z4m1Double $row.aggregate_audit_rel
        recomposition_error=Get-Z4m1Double $row.recomposition_rel_error
        reference_gap=Get-Z4m1Double $row.reference_gap 0.0
        true_residual=$case.true_residual;wall_seconds=$run.wall_seconds
        peak_rss_mb=$run.peak_working_set_mb
        candidate_apply_seconds=Get-Z4m1Double $row.fmm_local_candidate_apply_seconds 0.0
        candidate_median_apply_seconds=Get-Z4m1Double $row.fmm_local_candidate_median_apply_seconds 0.0
        fallback_requested=[string]$row.fmm_local_fallback_requested
        fallback_effective=[string]$row.fmm_local_fallback_effective
        fallback_seconds=Get-Z4m1Double $row.fmm_local_fallback_seconds 0.0
        worklist_hash=[string]$row.fmm_local_fallback_worklist_hash
        worklist_build_count=[int64]$row.fmm_local_fallback_worklist_build_count
        worklist_tasks=[int64]$row.fmm_local_fallback_tasks
        worklist_metadata_bytes=[int64]$row.fmm_local_fallback_metadata_bytes
        base_near_pairs=[int64]$row.fmm_local_fallback_base_near_pairs
        refined_pairs=[int64]$row.fmm_local_fallback_refined_pairs
        vectorized_sources=[int64]$row.fmm_local_fallback_vectorized_sources
        vectorized_lanes=[int64]$row.fmm_local_fallback_vectorized_lanes
        scalar_self_pairs=[int64]$row.fmm_local_fallback_scalar_self_pairs
        run_dir=$runDir
    }
}

function Get-Z4m1Decision([object[]]$Rows,[bool]$RequireAllFixtures) {
    if (@($Rows | Where-Object status -eq 'timeout').Count) { return 'z4m1_runtime_limited' }
    if (@($Rows | Where-Object status -eq 'avx2_unavailable').Count) { return 'z4m1_avx2_unavailable' }
    $packed = @($Rows | Where-Object variant -eq 'packed_cached')
    if ($RequireAllFixtures -and $packed.Count -ne 3) { return 'z4m1_runtime_limited' }
    if (@($packed | Where-Object {
        [double]$_.audit_rms -gt 8e-3 -or
        [double]$_.recomposition_error -gt 1e-10 -or
        [double]$_.reference_gap -gt 1e-10
    }).Count) { return 'z4m1_kernel_parity_gap' }
    if (@($packed | Where-Object {
        [double]$_.peak_rss_mb -ge 512 -or
        [int64]$_.worklist_metadata_bytes -ge 64MB
    }).Count) { return 'z4m1_memory_limited' }
    $partial = $false
    foreach ($packedRow in $packed) {
        $scalar = @($Rows | Where-Object {
            $_.fixture -eq $packedRow.fixture -and $_.variant -eq 'scalar_cached'
        }) | Select-Object -First 1
        $uncached = @($Rows | Where-Object {
            $_.fixture -eq $packedRow.fixture -and $_.variant -eq 'static_uncached'
        }) | Select-Object -First 1
        if ($null -eq $scalar -or $null -eq $uncached -or
            [double]$packedRow.fallback_seconds -le 0 -or
            [double]$packedRow.candidate_median_apply_seconds -le 0) {
            return 'z4m1_runtime_limited'
        }
        $fallbackSpeedup = [double]$scalar.fallback_seconds /
            [double]$packedRow.fallback_seconds
        $totalSpeedup = [double]$uncached.candidate_apply_seconds /
            [double]$packedRow.candidate_median_apply_seconds
        $auditDelta = [Math]::Abs([double]$scalar.audit_rms - [double]$packedRow.audit_rms)
        $packedRow | Add-Member fallback_speedup $fallbackSpeedup -Force
        $packedRow | Add-Member total_speedup $totalSpeedup -Force
        $packedRow | Add-Member scalar_audit_delta $auditDelta -Force
        if ($auditDelta -gt 1e-12) { return 'z4m1_kernel_parity_gap' }
        if ($fallbackSpeedup -lt 4.0 -or $totalSpeedup -lt 3.0) { $partial = $true }
    }
    if ($partial) { return 'z4m1_accuracy_pass_partial_speedup' }
    return 'z4m1_accuracy_speed_pass'
}

$matrix = New-Object System.Collections.Generic.List[object]
if ($PackageOnly) {
    foreach ($row in @(Import-Csv -LiteralPath $matrixCsv)) {
        $matrix.Add($row) | Out-Null
    }
} elseif ($ResumeQuick) {
    foreach ($row in @(Import-Csv -LiteralPath $matrixCsv | Where-Object fixture -eq '1sh1')) {
        $matrix.Add($row) | Out-Null
    }
    if ($matrix.Count -ne 3) {
        throw "ResumeQuick expected three 1SH1 rows, found $($matrix.Count)"
    }
} else {
    $quickFixture = Get-Z4m1Fixture '1sh1'
    foreach ($variant in @('static_uncached','scalar_cached','packed_cached')) {
        $matrix.Add((Invoke-Z4m1Case $quickFixture $variant)) | Out-Null
    }
}
$quickRows = @($matrix.ToArray() | Where-Object fixture -eq '1sh1')
$quickDecision = Get-Z4m1Decision $quickRows $false

if (-not $Quick -and -not $PackageOnly -and
        $quickDecision -eq 'z4m1_accuracy_speed_pass') {
    foreach ($name in @('barnase','1a63')) {
        $fixture = Get-Z4m1Fixture $name
        foreach ($variant in @('static_uncached','scalar_cached','packed_cached')) {
            $matrix.Add((Invoke-Z4m1Case $fixture $variant)) | Out-Null
        }
    }
}

$rows = $matrix.ToArray()
$decision = if ($Quick) { $quickDecision } elseif ($quickDecision -ne 'z4m1_accuracy_speed_pass') {
    $quickDecision
} else {
    Get-Z4m1Decision $rows $true
}
foreach ($row in $rows) {
    foreach ($property in @('fallback_speedup','total_speedup','scalar_audit_delta')) {
        if (-not ($row.PSObject.Properties.Name -contains $property)) {
            $row | Add-Member -NotePropertyName $property `
                -NotePropertyValue ([double]::NaN) -Force
        }
    }
}
if (-not $PackageOnly) {
    $rows | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
}
if (-not $PackageOnly) {
    $rows | Select-Object fixture,variant,fallback_requested,fallback_effective,
        worklist_hash,worklist_build_count,worklist_tasks,worklist_metadata_bytes,
        base_near_pairs,refined_pairs,vectorized_sources,vectorized_lanes,
        scalar_self_pairs,fallback_seconds,fallback_speedup,total_speedup |
        Export-Csv -LiteralPath $worklistCsv -NoTypeInformation
    [pscustomobject]@{decision=$decision;quick_decision=$quickDecision;quick=[bool]$Quick
        rows=$rows.Count;full_matrix_started=(-not $Quick -and $quickDecision -eq 'z4m1_accuracy_speed_pass')} |
        Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
    @('scope: Z4M.1 packed AVX2 direct-fallback acceleration',"decision: $decision",
      "quick_decision: $quickDecision","rows: $($rows.Count)",
      'target_batch: 4','source_tile: 32','static_refinement_plan: unchanged',
      'gmres_p4: fail-closed','production_default: grid remains default') |
        Set-Content -LiteralPath $summaryPath -Encoding ASCII
}

$svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="620" viewBox="0 0 1000 620">',
    '<rect width="1000" height="620" fill="white"/>',
    '<text x="35" y="38" font-family="Arial" font-size="20">Packed AVX2 fallback accuracy and runtime</text>')
$y = 82
foreach ($row in $rows) {
    $svg += ('<text x="35" y="{0}" font-family="Arial" font-size="13">{1} {2}: audit={3}, fallback={4}s, median={5}s, RSS={6} MB</text>' -f
        $y,$row.fixture,$row.variant,$row.audit_rms,$row.fallback_seconds,
        $row.candidate_median_apply_seconds,$row.peak_rss_mb)
    $y += 38
}
$svg += '</svg>'
if (-not $PackageOnly) {
    $svg | Set-Content -LiteralPath $plotPath -Encoding ASCII
}

$report = @('# TITAN FMM/M2L Scope Z4M.1 Packed Fallback Report','',
    "Decision: $decision",'',
    'The static refinement plan and p4 far operator are unchanged. This scope packs direct fallback work into deterministic 4-target by 32-source tasks and reuses the legacy AVX2 centroid primitive.','')
foreach ($row in $rows) {
    $speedupText = if ($row.PSObject.Properties.Name -contains 'fallback_speedup') {
        ', fallback_speedup={0}x, total_speedup={1}x' -f
            ([Math]::Round([double]$row.fallback_speedup,2)),
            ([Math]::Round([double]$row.total_speedup,2))
    } else { '' }
    $report += ('- {0} {1}: audit={2}, fallback={3}s, median_apply={4}s, tasks={5}, metadata={6} MB, RSS={7} MB' -f
        $row.fixture,$row.variant,$row.audit_rms,$row.fallback_seconds,
        $row.candidate_median_apply_seconds,$row.worklist_tasks,
        ([Math]::Round([double]$row.worklist_metadata_bytes/1MB,2)),$row.peak_rss_mb) +
        $speedupText
}
$report += @('','GMRES-p4 remains fail-closed.','Grid remains the production Poisson-Boltzmann default.')
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

$bad = @(Get-ChildItem -LiteralPath $WorkRoot -Recurse -Filter arg.in | ForEach-Object {
    $path = $_.FullName
    Get-Content -LiteralPath $path | Where-Object {
        $_ -match '^\s*-' -and $_ -notmatch '^\s*-[^=\s]+='
    } | ForEach-Object { "$path :: $_" }
})
if ($bad.Count) { throw "Generated arg.in contains bare options: $($bad -join '; ')" }
if (-not $AllowFailedDecision -and $decision -ne 'z4m1_accuracy_speed_pass') {
    throw "Scope Z4M.1 decision: $decision"
}
Write-Host "PASS Scope Z4M.1 decision=$decision rows=$($rows.Count)"
