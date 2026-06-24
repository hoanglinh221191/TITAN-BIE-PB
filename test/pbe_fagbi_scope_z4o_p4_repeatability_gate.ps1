param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4o'),
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

$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4o_matrix.csv'
$repeatCsv = Join-Path $WorkRoot 'fagbi_scope_z4o_repeatability.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4o_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4o_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4o_repeatability_plot.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4O_P4_REPEATABILITY_REPORT.txt'
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

$manifest = @(Import-Csv -LiteralPath $manifestPath)

function Get-Z4oFixture([string]$Name) {
    $row = @($manifest | Where-Object fixture -eq $Name)
    if ($row.Count -ne 1) { throw "Expected one public fixture row for $Name" }
    return $row[0]
}

function Get-Z4oDouble([object]$Value,[double]$Default=[double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Get-Z4oInt64([object]$Value,[int64]$Default=0) {
    $parsed = [int64]0
    if ($null -eq $Value -or -not [int64]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Get-Z4oProfileValue([object]$Profile,[string]$Name,[object]$Default='') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function New-Z4oFailureRow(
    [pscustomobject]$Fixture,
    [string]$Variant,
    [int]$Repeat,
    [string]$Status,
    [object]$Run,
    [string]$RunDir,
    [string]$FailureLog = ''
) {
    return [pscustomobject]@{
        fixture=$Fixture.fixture; protein_name=$Fixture.protein_name
        pdb_id=$Fixture.pdb_id; variant=$Variant; repeat=$Repeat; status=$Status
        panels=0; iterations=0; internal_true_residual=[double]::NaN
        final_residual=[double]::NaN; energy_kj_mol=[double]::NaN
        audit_rms=[double]::NaN; matvec_seconds_avg=[double]::NaN
        wall_seconds=$Run.wall_seconds; peak_rss_mb=$Run.peak_working_set_mb
        fmm_local_gmres_active='false'; plan_build_count=0; cache_build_count=0
        gmres_matvec_calls=0; sampled_direct_true_residual=[double]::NaN
        operator_finite='false'; p4_vs_degree1_matvec_ratio=[double]::NaN
        failure_log=$FailureLog; run_dir=$RunDir
    }
}

function Invoke-Z4oCase([pscustomobject]$Fixture,[string]$Variant,[int]$Repeat) {
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_${Variant}_r$Repeat"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $isP4 = $Variant -eq 'p4_guarded'
    $pqr = Join-Path $FixtureRoot $Fixture.pqr
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
        return New-Z4oFailureRow $Fixture $Variant $Repeat 'timeout' $run $runDir
    }
    if ($run.exit_code -ne 0) {
        $stdoutPath = Join-Path $runDir 'stdout.log'
        $stderrPath = Join-Path $runDir 'stderr.log'
        $log = ((if (Test-Path $stdoutPath) { Get-Content -Raw $stdoutPath } else { '' }) +
            (if (Test-Path $stderrPath) { Get-Content -Raw $stderrPath } else { '' }))
        $status = if ($log -match 'cache_budget_exceeded') {
            'cache_budget_exceeded'
        } elseif ($log -match 'non_finite') { 'non_finite' } else { 'failed' }
        return New-Z4oFailureRow $Fixture $Variant $Repeat $status $run $runDir $log
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
    return [pscustomobject]@{
        fixture=$Fixture.fixture; protein_name=$Fixture.protein_name
        pdb_id=$Fixture.pdb_id; variant=$Variant; repeat=$Repeat; status=$case.status
        panels=$case.panels; iterations=$case.iterations
        internal_true_residual=$case.true_residual
        final_residual=$case.final_residual
        energy_kj_mol=$case.energy_kj_mol
        audit_rms=$case.tree_vs_direct_matvec_rel
        matvec_seconds_avg=$case.matvec_seconds_avg
        wall_seconds=$run.wall_seconds
        peak_rss_mb=$run.peak_working_set_mb
        fmm_local_gmres_active=[string](Get-Z4oProfileValue $profile 'fmm_local_gmres_active' 'false')
        plan_build_count=Get-Z4oInt64 (Get-Z4oProfileValue $profile 'fmm_local_plan_build_count' 0)
        cache_build_count=Get-Z4oInt64 (Get-Z4oProfileValue $profile 'fmm_local_cache_build_count' 0)
        gmres_matvec_calls=Get-Z4oInt64 (Get-Z4oProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
        sampled_direct_true_residual=Get-Z4oDouble (Get-Z4oProfileValue $profile 'fmm_local_sampled_direct_true_residual' 0.0) 0.0
        operator_finite=[string](Get-Z4oProfileValue $profile 'fmm_local_operator_finite' 'false')
        p4_vs_degree1_matvec_ratio=[double]::NaN
        failure_log=''; run_dir=$RunDir
    }
}

function Test-Z4oP4Row([object]$Row) {
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

function Get-Z4oMedian([double[]]$Values) {
    $valid = @($Values | Where-Object { Test-FmmCliFinite $_ } | Sort-Object)
    if (-not $valid.Count) { return [double]::NaN }
    $mid = [int]($valid.Count / 2)
    if ($valid.Count % 2) { return [double]$valid[$mid] }
    return ([double]$valid[$mid - 1] + [double]$valid[$mid]) / 2.0
}

function Write-Z4oSvg([object[]]$Rows,[string]$Path) {
    $svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="620" viewBox="0 0 1120 620">',
        '<rect width="1120" height="620" fill="white"/>',
        '<text x="34" y="38" font-family="Arial" font-size="20">Z4O guarded p4 repeatability gate</text>',
        '<text x="34" y="64" font-family="Arial" font-size="13">bars show p4 audit RMS; labels show median matvec ratio</text>',
        '<line x1="90" y1="500" x2="1040" y2="500" stroke="#222"/>',
        '<line x1="90" y1="96" x2="90" y2="500" stroke="#222"/>')
    $i = 0
    foreach ($row in $Rows) {
        $x = 130 + 150 * $i
        $audit = if (Test-FmmCliFinite ([double]$row.audit_max)) { [double]$row.audit_max } else { 1.0 }
        $bar = [Math]::Min(360.0, 70.0 * (-[Math]::Log10([Math]::Max(1e-8, $audit))))
        $color = if ($row.repeatability_pass -eq 'true') { '#16a34a' } else { '#f59e0b' }
        $svg += "<rect x='$x' y='$([int](500-$bar))' width='62' height='$([int]$bar)' fill='$color'/>"
        $svg += "<text x='$($x-4)' y='526' font-family='Arial' font-size='12'>$($row.fixture)</text>"
        $svg += "<text x='$($x-8)' y='$([int](488-$bar))' font-family='Arial' font-size='11'>$('{0:G3}' -f [double]$row.median_matvec_ratio)x</text>"
        ++$i
    }
    $svg += '</svg>'
    $svg | Set-Content -LiteralPath $Path -Encoding ASCII
}

$fixtureNames = if ($Quick) { @('1sh1') } else { @('barnase','1a63','1sh1') }
$repeatCount = if ($Quick) { 2 } else { 3 }
$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    $fixture = Get-Z4oFixture $fixtureName
    $baseline = Invoke-Z4oCase $fixture 'degree1_baseline' 0
    $rows.Add($baseline) | Out-Null
    for ($repeat = 1; $repeat -le $repeatCount; ++$repeat) {
        $rows.Add((Invoke-Z4oCase $fixture 'p4_guarded' $repeat)) | Out-Null
    }
}
$matrix = @($rows.ToArray())

foreach ($fixtureName in $fixtureNames) {
    $baseline = @($matrix | Where-Object {
        $_.fixture -eq $fixtureName -and $_.variant -eq 'degree1_baseline'
    } | Select-Object -First 1)
    if (-not $baseline.Count -or
        -not (Test-FmmCliFinite ([double]$baseline[0].matvec_seconds_avg)) -or
        [double]$baseline[0].matvec_seconds_avg -le 0.0) { continue }
    foreach ($row in @($matrix | Where-Object {
            $_.fixture -eq $fixtureName -and $_.variant -eq 'p4_guarded'
        })) {
        $row.p4_vs_degree1_matvec_ratio =
            [double]$row.matvec_seconds_avg / [double]$baseline[0].matvec_seconds_avg
    }
}

$repeatRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    $p4Rows = @($matrix | Where-Object { $_.fixture -eq $fixtureName -and $_.variant -eq 'p4_guarded' })
    $energies = @($p4Rows | ForEach-Object { [double]$_.energy_kj_mol })
    $audits = @($p4Rows | ForEach-Object { [double]$_.audit_rms })
    $directResiduals = @($p4Rows | ForEach-Object { [double]$_.sampled_direct_true_residual })
    $ratios = @($p4Rows | ForEach-Object { [double]$_.p4_vs_degree1_matvec_ratio })
    $meanEnergy = if ($energies.Count) { ($energies | Measure-Object -Average).Average } else { [double]::NaN }
    $energySpread = if ($energies.Count) {
        ([double]($energies | Measure-Object -Maximum).Maximum -
         [double]($energies | Measure-Object -Minimum).Minimum)
    } else { [double]::NaN }
    $energyGate = if (Test-FmmCliFinite $meanEnergy) {
        [Math]::Max(0.5, 1.0e-4 * [Math]::Abs($meanEnergy))
    } else { [double]::NaN }
    $auditMax = if ($audits.Count) { [double]($audits | Measure-Object -Maximum).Maximum } else { [double]::NaN }
    $directResidualMax = if ($directResiduals.Count) { [double]($directResiduals | Measure-Object -Maximum).Maximum } else { [double]::NaN }
    $ratioMax = if ($ratios.Count) { [double]($ratios | Measure-Object -Maximum).Maximum } else { [double]::NaN }
    $ratioMedian = Get-Z4oMedian $ratios
    $rowPasses = @($p4Rows | Where-Object { Test-Z4oP4Row $_ }).Count
    $repeatPass = $p4Rows.Count -eq $repeatCount -and $rowPasses -eq $repeatCount -and
        (Test-FmmCliFinite $energySpread) -and (Test-FmmCliFinite $energyGate) -and
        $energySpread -le $energyGate -and
        (Test-FmmCliFinite $auditMax) -and $auditMax -le 8.0e-3 -and
        (Test-FmmCliFinite $directResidualMax) -and $directResidualMax -le 8.0e-3 -and
        (Test-FmmCliFinite $ratioMedian) -and $ratioMedian -le 1.0 -and
        (Test-FmmCliFinite $ratioMax) -and $ratioMax -le 1.25
    $repeatRows.Add([pscustomobject]@{
        fixture=$fixtureName; repeats=$p4Rows.Count; row_passes=$rowPasses
        mean_energy=$meanEnergy; energy_spread=$energySpread; energy_gate=$energyGate
        audit_max=$auditMax; sampled_direct_true_residual_max=$directResidualMax
        median_matvec_ratio=$ratioMedian; max_matvec_ratio=$ratioMax
        peak_rss_mb_max=if ($p4Rows.Count) { [double]($p4Rows.peak_rss_mb | Measure-Object -Maximum).Maximum } else { [double]::NaN }
        repeatability_pass=[string]$repeatPass
    }) | Out-Null
}
$repeatability = @($repeatRows.ToArray())

$p4All = @($matrix | Where-Object variant -eq 'p4_guarded')
$timeouts = @($p4All | Where-Object status -eq 'timeout')
$memoryRows = @($p4All | Where-Object {
    $_.status -eq 'cache_budget_exceeded' -or
    ((Test-FmmCliFinite ([double]$_.peak_rss_mb)) -and [double]$_.peak_rss_mb -ge 512.0)
})
$lifecycleRows = @($p4All | Where-Object {
    $_.status -eq 'ok' -and
    ([int64]$_.plan_build_count -ne 1 -or
     [int64]$_.cache_build_count -ne 1 -or
     [int64]$_.gmres_matvec_calls -le 0 -or
     $_.fmm_local_gmres_active -ne 'true')
})
$accuracyRows = @($p4All | Where-Object { -not (Test-Z4oP4Row $_) })
$runtimeRepeatRows = @($repeatability | Where-Object {
    $_.repeatability_pass -ne 'true' -and
    (Test-FmmCliFinite ([double]$_.audit_max)) -and [double]$_.audit_max -le 8.0e-3 -and
    (Test-FmmCliFinite ([double]$_.sampled_direct_true_residual_max)) -and
    [double]$_.sampled_direct_true_residual_max -le 8.0e-3
})

$decision = 'z4o_repeatability_pass'
if ($timeouts.Count) {
    $decision = 'z4o_runtime_limited'
} elseif ($memoryRows.Count) {
    $decision = 'z4o_memory_limited'
} elseif ($lifecycleRows.Count) {
    $decision = 'z4o_lifecycle_gap'
} elseif ($accuracyRows.Count) {
    $decision = 'z4o_repeatability_accuracy_gap'
} elseif (@($repeatability | Where-Object repeatability_pass -ne 'true').Count -eq 0) {
    $decision = 'z4o_repeatability_pass'
} elseif ($runtimeRepeatRows.Count) {
    $decision = 'z4o_repeatability_runtime_gap'
} else {
    $decision = 'z4o_repeatability_accuracy_gap'
}

$matrix | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$repeatability | Export-Csv -LiteralPath $repeatCsv -NoTypeInformation
[pscustomobject]@{
    decision=$decision; quick=[bool]$Quick; rows=$matrix.Count
    repeat_fixtures=$repeatability.Count; p4_rows=$p4All.Count
} | Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
@(
    'scope: Z4O guarded p4 repeatability gate',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "rows: $($matrix.Count)",
    "repeat_count: $repeatCount",
    'p4_policy: operator-tail-static tol=0.002 strict-geometry-mac theta=0.5',
    'gmres_p4: explicit advanced diagnostic only',
    'tabi_treecode_default: degree1',
    'production_default: grid'
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII
Write-Z4oSvg $repeatability $plotPath

$report = @(
    '# TITAN FMM/M2L Scope Z4O p4 Repeatability Report',
    '',
    "Decision: $decision",
    '',
    'This gate measures repeatability for the explicit guarded p4 path. It does not change grid or TABI treecode defaults.',
    ''
)
foreach ($row in $repeatability) {
    $report += ('- {0}: repeats={1}, audit_max={2}, direct_residual_max={3}, energy_spread={4}, median_ratio={5}, max_ratio={6}, peak_RSS={7} MB, pass={8}' -f
        $row.fixture,$row.repeats,$row.audit_max,
        $row.sampled_direct_true_residual_max,$row.energy_spread,
        $row.median_matvec_ratio,$row.max_matvec_ratio,
        $row.peak_rss_mb_max,$row.repeatability_pass)
}
$report += @(
    '',
    'Grid remains the production PBE default. TABI treecode remains degree-1 unless gmres-p4 is explicitly selected.'
)
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

if ($decision -ne 'z4o_repeatability_pass' -and -not $AllowFailedDecision) {
    throw "Scope Z4O gate did not pass: $decision"
}
Write-Host "Z4O repeatability decision: $decision"
