param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4n4'),
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

$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4n4_matrix.csv'
$convergenceCsv = Join-Path $WorkRoot 'fagbi_scope_z4n4_convergence.csv'
$auditCsv = Join-Path $WorkRoot 'fagbi_scope_z4n4_operator_audit.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4n4_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4n4_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4n4_policy_plot.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4N4_GUARDED_P4_POLICY_REPORT.txt'
$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)

$RuntimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
        $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
    }
    $RuntimeRepoRoot = Split-Path -Parent $commonGitDir
}

$generatedPqr = Join-Path $WorkRoot 'generated_direct_reference.pqr'
@(
    'ATOM      1  C   GLY A   1      -6.000   0.000   0.000  1.0000  1.7000',
    'ATOM      2  O   GLY A   1      -2.000   0.000   0.000 -1.0000  1.5000',
    'ATOM      3  N   GLY A   2       2.000   0.000   0.000  1.0000  1.5500',
    'ATOM      4  O   GLY A   2       6.000   0.000   0.000 -1.0000  1.5000'
) | Set-Content -LiteralPath $generatedPqr -Encoding ASCII

$manifest = @(Import-Csv -LiteralPath $manifestPath)

function Get-Z4n4Fixture([string]$Name) {
    if ($Name -eq 'generated') {
        return [pscustomobject]@{
            fixture='generated'; protein_name='generated direct reference'
            pdb_id='generated'; pqr=$generatedPqr
        }
    }
    $row = @($manifest | Where-Object fixture -eq $Name)
    if ($row.Count -ne 1) { throw "Expected one public fixture row for $Name" }
    return $row[0]
}

function Get-Z4n4Double([object]$Value,[double]$Default=[double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Get-Z4n4Int64([object]$Value,[int64]$Default=0) {
    $parsed = [int64]0
    if ($null -eq $Value -or -not [int64]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Get-Z4n4ProfileValue([object]$Profile,[string]$Name,[object]$Default='') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function New-Z4n4FailureRow(
    [pscustomobject]$Fixture,
    [string]$Variant,
    [string]$Status,
    [object]$Run,
    [string]$RunDir,
    [string]$FailureLog = ''
) {
    return [pscustomobject]@{
        fixture=$Fixture.fixture; protein_name=$Fixture.protein_name
        pdb_id=$Fixture.pdb_id; variant=$Variant; status=$Status
        panels=0; iterations=0; internal_true_residual=[double]::NaN
        final_residual=[double]::NaN; energy_kj_mol=[double]::NaN
        audit_rms=[double]::NaN; matvec_seconds_avg=[double]::NaN
        wall_seconds=$Run.wall_seconds; peak_rss_mb=$Run.peak_working_set_mb
        fmm_local_gmres_active='false'; prepare_seconds=0.0
        plan_build_count=0; cache_build_count=0; gmres_matvec_calls=0
        gmres_apply_seconds_total=0.0; gmres_apply_seconds_mean=0.0
        moment_seconds_total=0.0; m2l_seconds_total=0.0
        fallback_seconds_total=0.0; l2p_seconds_total=0.0
        sampled_direct_true_residual=[double]::NaN; operator_finite='false'
        geometry_guard=''; geometry_theta=[double]::NaN
        residual_history=''; p4_vs_degree1_matvec_ratio=[double]::NaN
        failure_log=$FailureLog; run_dir=$RunDir
    }
}

function Invoke-Z4n4Case([pscustomobject]$Fixture,[string]$Variant) {
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_$Variant"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $isP4 = $Variant -eq 'p4_guarded'
    $pqr = if ($Fixture.fixture -eq 'generated') { $Fixture.pqr } else {
        Join-Path $FixtureRoot $Fixture.pqr
    }
    $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RuntimeRepoRoot `
        -PqrPath $pqr -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -SurfacePoints 8 -MaxIter 720 `
        -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount 256 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder 4 -FmmCartesianAudit off `
        -FmmLocal $(if ($isP4) { 'gmres-p4' } else { 'off' }) `
        -FmmLocalRefine $(if ($isP4) { 'operator-tail-static' } else { 'off' }) `
        -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache $(if ($isP4) { 'green-jet-compressed' } else { 'off' }) `
        -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback $(if ($isP4) { 'packed-avx2' } else { 'scalar' }) `
        -FmmLocalPlanBuilder $(if ($isP4) { 'metric-cache' } else { 'legacy-rescan' }) `
        -FmmLocalRefineProbeCount 4 -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard $(if ($isP4) { 'strict-geometry-mac' } else { 'off' }) `
        -FmmLocalGeometryTheta $(if ($isP4) { 0.5 } else { 0.0 })
    $bare = @(Get-Content -LiteralPath $arg | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) { throw "Generated arg.in contains bare options: $($bare -join ', ')" }
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) {
        return New-Z4n4FailureRow $Fixture $Variant 'timeout' $run $runDir
    }
    if ($run.exit_code -ne 0) {
        $stdoutPath = Join-Path $runDir 'stdout.log'
        $stderrPath = Join-Path $runDir 'stderr.log'
        $log = ((if (Test-Path $stdoutPath) { Get-Content -Raw $stdoutPath } else { '' }) +
            (if (Test-Path $stderrPath) { Get-Content -Raw $stderrPath } else { '' }))
        $status = if ($log -match 'cache_budget_exceeded') {
            'cache_budget_exceeded'
        } elseif ($log -match 'non_finite') { 'non_finite' } else { 'failed' }
        return New-Z4n4FailureRow $Fixture $Variant $status $run $runDir $log
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
    return [pscustomobject]@{
        fixture=$Fixture.fixture; protein_name=$Fixture.protein_name
        pdb_id=$Fixture.pdb_id; variant=$Variant; status=$case.status
        panels=$case.panels; iterations=$case.iterations
        internal_true_residual=$case.true_residual
        final_residual=$case.final_residual
        energy_kj_mol=$case.energy_kj_mol
        audit_rms=$case.tree_vs_direct_matvec_rel
        matvec_seconds_avg=$case.matvec_seconds_avg
        wall_seconds=$run.wall_seconds
        peak_rss_mb=$run.peak_working_set_mb
        fmm_local_gmres_active=[string](Get-Z4n4ProfileValue $profile 'fmm_local_gmres_active' 'false')
        prepare_seconds=Get-Z4n4Double (Get-Z4n4ProfileValue $profile 'fmm_local_prepare_seconds' 0.0) 0.0
        plan_build_count=Get-Z4n4Int64 (Get-Z4n4ProfileValue $profile 'fmm_local_plan_build_count' 0)
        cache_build_count=Get-Z4n4Int64 (Get-Z4n4ProfileValue $profile 'fmm_local_cache_build_count' 0)
        gmres_matvec_calls=Get-Z4n4Int64 (Get-Z4n4ProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
        gmres_apply_seconds_total=Get-Z4n4Double (Get-Z4n4ProfileValue $profile 'fmm_local_gmres_apply_seconds_total' 0.0) 0.0
        gmres_apply_seconds_mean=Get-Z4n4Double (Get-Z4n4ProfileValue $profile 'fmm_local_gmres_apply_seconds_mean' 0.0) 0.0
        moment_seconds_total=Get-Z4n4Double (Get-Z4n4ProfileValue $profile 'fmm_local_moment_seconds_total' 0.0) 0.0
        m2l_seconds_total=Get-Z4n4Double (Get-Z4n4ProfileValue $profile 'fmm_local_m2l_seconds_total' 0.0) 0.0
        fallback_seconds_total=Get-Z4n4Double (Get-Z4n4ProfileValue $profile 'fmm_local_fallback_seconds_total' 0.0) 0.0
        l2p_seconds_total=Get-Z4n4Double (Get-Z4n4ProfileValue $profile 'fmm_local_l2p_seconds_total' 0.0) 0.0
        sampled_direct_true_residual=Get-Z4n4Double (Get-Z4n4ProfileValue $profile 'fmm_local_sampled_direct_true_residual' 0.0) 0.0
        operator_finite=[string](Get-Z4n4ProfileValue $profile 'fmm_local_operator_finite' 'false')
        geometry_guard=[string](Get-Z4n4ProfileValue $profile 'fmm_local_geometry_guard' '')
        geometry_theta=Get-Z4n4Double (Get-Z4n4ProfileValue $profile 'fmm_local_geometry_theta' 0.0) 0.0
        residual_history=$case.residual_history
        p4_vs_degree1_matvec_ratio=[double]::NaN
        failure_log=''; run_dir=$runDir
    }
}

function Test-Z4n4P4Accuracy([object]$Row) {
    return $Row.status -eq 'ok' -and
        (Test-FmmCliFinite ([double]$Row.internal_true_residual)) -and
        (Test-FmmCliFinite ([double]$Row.audit_rms)) -and
        (Test-FmmCliFinite ([double]$Row.sampled_direct_true_residual)) -and
        (Test-FmmCliFinite ([double]$Row.energy_kj_mol)) -and
        [double]$Row.internal_true_residual -le 1.0e-6 -and
        [double]$Row.audit_rms -le 8.0e-3 -and
        [double]$Row.sampled_direct_true_residual -le 8.0e-3 -and
        [double]$Row.peak_rss_mb -lt 512.0 -and
        [int64]$Row.plan_build_count -eq 1 -and
        [int64]$Row.cache_build_count -eq 1 -and
        [int64]$Row.gmres_matvec_calls -gt 0 -and
        $Row.fmm_local_gmres_active -eq 'true' -and
        $Row.operator_finite -eq 'true'
}

function Write-Z4n4Svg([object[]]$Rows,[string]$Path) {
    $svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="640" viewBox="0 0 1120 640">',
        '<rect width="1120" height="640" fill="white"/>',
        '<text x="34" y="38" font-family="Arial" font-size="20">Z4N.4 guarded p4 GMRES policy gate</text>',
        '<text x="34" y="66" font-family="Arial" font-size="13">p4 policy: operator-tail-static tol=0.002, strict-geometry-mac theta=0.5</text>')
    $y = 104
    foreach ($row in $Rows) {
        $ratio = if (Test-FmmCliFinite ([double]$row.p4_vs_degree1_matvec_ratio)) {
            '{0:G3}x' -f [double]$row.p4_vs_degree1_matvec_ratio
        } else { 'n/a' }
        $svg += ('<text x="34" y="{0}" font-family="Arial" font-size="13">{1} {2}: status={3}, residual={4}, audit={5}, direct-res={6}, matvec={7}s, ratio={8}, RSS={9} MB</text>' -f
            $y,$row.fixture,$row.variant,$row.status,$row.internal_true_residual,
            $row.audit_rms,$row.sampled_direct_true_residual,$row.matvec_seconds_avg,
            $ratio,$row.peak_rss_mb)
        $y += 32
    }
    $svg += '</svg>'
    $svg | Set-Content -LiteralPath $Path -Encoding ASCII
}

$fixtureNames = if ($Quick) { @('generated','1sh1') } else { @('barnase','1a63','1sh1') }
$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    $fixture = Get-Z4n4Fixture $fixtureName
    $rows.Add((Invoke-Z4n4Case $fixture 'degree1_baseline')) | Out-Null
    $rows.Add((Invoke-Z4n4Case $fixture 'p4_guarded')) | Out-Null
}
$matrix = @($rows.ToArray())

foreach ($fixtureName in $fixtureNames) {
    $baseline = @($matrix | Where-Object {
        $_.fixture -eq $fixtureName -and $_.variant -eq 'degree1_baseline'
    } | Select-Object -First 1)
    $p4 = @($matrix | Where-Object {
        $_.fixture -eq $fixtureName -and $_.variant -eq 'p4_guarded'
    } | Select-Object -First 1)
    if ($baseline.Count -and $p4.Count -and
        (Test-FmmCliFinite ([double]$baseline[0].matvec_seconds_avg)) -and
        [double]$baseline[0].matvec_seconds_avg -gt 0) {
        $p4[0].p4_vs_degree1_matvec_ratio =
            [double]$p4[0].matvec_seconds_avg / [double]$baseline[0].matvec_seconds_avg
    }
}

$trace = New-Object System.Collections.Generic.List[object]
foreach ($row in $matrix) {
    $values = @(([string]$row.residual_history) -split ';' | Where-Object { $_ })
    for ($iteration = 0; $iteration -lt $values.Count; ++$iteration) {
        $trace.Add([pscustomobject]@{
            fixture=$row.fixture; variant=$row.variant
            iteration=$iteration; residual=Get-Z4n4Double $values[$iteration]
        }) | Out-Null
    }
}

$requiredP4 = @($matrix | Where-Object variant -eq 'p4_guarded')
$timeouts = @($requiredP4 | Where-Object status -eq 'timeout')
$memoryRows = @($requiredP4 | Where-Object {
    $_.status -eq 'cache_budget_exceeded' -or
    ((Test-FmmCliFinite ([double]$_.peak_rss_mb)) -and [double]$_.peak_rss_mb -ge 512.0)
})
$accuracyRows = @($requiredP4 | Where-Object { Test-Z4n4P4Accuracy $_ })
$convergenceGaps = @($requiredP4 | Where-Object {
    $_.status -ne 'ok' -or
    -not (Test-FmmCliFinite ([double]$_.internal_true_residual)) -or
    [double]$_.internal_true_residual -gt 1.0e-6 -or
    -not (Test-FmmCliFinite ([double]$_.energy_kj_mol))
})
$operatorGaps = @($requiredP4 | Where-Object {
    $_.status -eq 'ok' -and
    ((Test-FmmCliFinite ([double]$_.audit_rms)) -and
     (Test-FmmCliFinite ([double]$_.sampled_direct_true_residual))) -and
    ([double]$_.audit_rms -gt 8.0e-3 -or
     [double]$_.sampled_direct_true_residual -gt 8.0e-3 -or
     $_.operator_finite -ne 'true')
})
$lifecycleGaps = @($requiredP4 | Where-Object {
    $_.status -eq 'ok' -and
    ([int64]$_.plan_build_count -ne 1 -or
     [int64]$_.cache_build_count -ne 1 -or
     [int64]$_.gmres_matvec_calls -le 0 -or
     $_.fmm_local_gmres_active -ne 'true')
})
$runtimeGaps = @($requiredP4 | Where-Object {
    $_.status -eq 'ok' -and
    (Test-FmmCliFinite ([double]$_.p4_vs_degree1_matvec_ratio)) -and
    [double]$_.p4_vs_degree1_matvec_ratio -gt 2.0
})

$decision = 'z4n4_guarded_p4_policy_pass'
if ($timeouts.Count) {
    $decision = 'z4n4_runtime_limited'
} elseif ($memoryRows.Count) {
    $decision = 'z4n4_memory_limited'
} elseif ($convergenceGaps.Count) {
    $decision = 'z4n4_gmres_convergence_gap'
} elseif ($operatorGaps.Count) {
    $decision = 'z4n4_operator_accuracy_gap'
} elseif ($lifecycleGaps.Count) {
    $decision = 'z4n4_cache_lifecycle_gap'
} elseif ($accuracyRows.Count -eq $requiredP4.Count -and $runtimeGaps.Count) {
    $decision = 'z4n4_accuracy_pass_runtime_gap'
}

$matrix | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$trace.ToArray() | Export-Csv -LiteralPath $convergenceCsv -NoTypeInformation
$requiredP4 | Select-Object fixture,audit_rms,sampled_direct_true_residual,
    internal_true_residual,operator_finite,plan_build_count,cache_build_count,
    gmres_matvec_calls,p4_vs_degree1_matvec_ratio |
    Export-Csv -LiteralPath $auditCsv -NoTypeInformation
[pscustomobject]@{
    decision=$decision; quick=[bool]$Quick; rows=$matrix.Count
    required_p4_rows=$requiredP4.Count; pass_rows=$accuracyRows.Count
} | Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
@(
    'scope: Z4N.4 guarded p4 GMRES policy gate',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "rows: $($matrix.Count)",
    'p4_policy: operator-tail-static tol=0.002 strict-geometry-mac theta=0.5',
    'gmres_p4: hidden explicit only',
    'tabi_treecode_default: degree1',
    'production_default: grid'
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII

Write-Z4n4Svg $matrix $plotPath

$report = @(
    '# TITAN FMM/M2L Scope Z4N.4 Guarded p4 Policy Gate',
    '',
    "Decision: $decision",
    '',
    'The p4 operator remains explicit-only. This gate validates the Z4N.3 strict-geometry budget policy inside full GMRES.',
    ''
)
foreach ($row in $matrix) {
    $report += ('- {0} {1}: status={2}, internal_residual={3}, direct_sampled_residual={4}, audit={5}, matvec={6}s, ratio={7}, wall={8}s, RSS={9} MB' -f
        $row.fixture,$row.variant,$row.status,$row.internal_true_residual,
        $row.sampled_direct_true_residual,$row.audit_rms,$row.matvec_seconds_avg,
        $row.p4_vs_degree1_matvec_ratio,$row.wall_seconds,$row.peak_rss_mb)
}
$report += @(
    '',
    'Grid remains the production PBE default. TABI treecode remains degree-1 unless gmres-p4 is explicitly selected.'
)
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

if ($decision -ne 'z4n4_guarded_p4_policy_pass' -and -not $AllowFailedDecision) {
    throw "Scope Z4N.4 gate did not pass: $decision"
}
Write-Host "Z4N.4 guarded p4 policy decision: $decision"
