param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4n3'),
    [int]$CaseTimeoutSeconds = 3600,
    [switch]$Quick,
    [switch]$ConfirmCandidate,
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

$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4n3_matrix.csv'
$budgetCsv = Join-Path $WorkRoot 'fagbi_scope_z4n3_budget.csv'
$confirmMatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4n3_confirm_matrix.csv'
$confirmDecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4n3_confirm_decision.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4n3_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4n3_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4n3_budget_plot.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4N3_OUTLIER_BUDGET_REPORT.txt'
$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if ($ConfirmCandidate) {
    New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
} else {
    New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
}
$manifest = @(Import-Csv -LiteralPath $manifestPath)

function Get-Z4n3Fixture([string]$Name) {
    $row = @($manifest | Where-Object fixture -eq $Name)
    if ($row.Count -ne 1) { throw "Expected one public fixture row for $Name" }
    return $row[0]
}

function Get-Z4n3Double([object]$Value,[double]$Default=[double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Get-Z4n3Int64([object]$Value,[int64]$Default=0) {
    $parsed = [int64]0
    if ($null -eq $Value -or -not [int64]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Get-Z4n3Prop([object]$Object,[string]$Name,[object]$Default='') {
    if ($null -eq $Object) { return $Default }
    $prop = $Object.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Convert-Z4n3ConditionName([double]$Value) {
    return ('{0:G4}' -f $Value).Replace('.','p').Replace('-','m')
}

function Invoke-Z4n3Case(
    [pscustomobject]$Fixture,
    [string]$Condition,
    [string]$RefineMode,
    [double]$RefineTolerance,
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
        -FmmLocalRefine $RefineMode -FmmLocalRefineTolerance $RefineTolerance `
        -FmmLocalCache green-jet-compressed -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback packed-avx2 -FmmLocalPlanBuilder metric-cache `
        -FmmLocalRefineProbeCount 4 -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard $GeometryGuard -FmmLocalGeometryTheta $GeometryTheta
    $bare = @(Get-Content -LiteralPath $arg | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) { throw "Generated arg.in contains bare options: $($bare -join ', ')" }

    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $base = [ordered]@{
        fixture=$Fixture.fixture;condition=$Condition;status='timeout'
        refine_mode=$RefineMode;refine_tolerance=$RefineTolerance
        geometry_guard=$GeometryGuard;geometry_theta=$GeometryTheta
        panels=0;iterations=0;candidate_audit_rms=[double]::NaN
        sampled_direct_true_residual=[double]::NaN
        baseline_fallback_pairs=0;fallback_ratio=[double]::NaN
        direct_fallback_pairs=0;far_pairs_after=0;peak_rss_mb=$run.peak_working_set_mb
        wall_seconds=$run.wall_seconds;plan_build_count=0;cache_build_count=0
        plan_hash='';linf_max_relative_tail=[double]::NaN
        linf_outlier_panels=0;geometry_rejected_cells=0;run_dir=$runDir
    }
    if ($run.timed_out) { return [pscustomobject]$base }
    if ($run.exit_code -ne 0) {
        $base.status = 'failed'
        $stdoutPath = Join-Path $runDir 'stdout.log'
        $stderrPath = Join-Path $runDir 'stderr.log'
        $base.failure_log = ((if (Test-Path $stdoutPath) { Get-Content -Raw $stdoutPath } else { '' }) +
            (if (Test-Path $stderrPath) { Get-Content -Raw $stderrPath } else { '' }))
        return [pscustomobject]$base
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

    $base.status = $case.status
    $base.panels = $case.panels
    $base.iterations = $case.iterations
    $base.candidate_audit_rms = if ($summaryRow) {
        Convert-FmmCliDouble $summaryRow.aggregate_audit_rel
    } else { [double]$case.tree_vs_direct_matvec_rel }
    $base.sampled_direct_true_residual =
        Get-Z4n3Double (Get-Z4n3Prop $profile 'fmm_local_sampled_direct_true_residual' 0.0) 0.0
    $base.direct_fallback_pairs = if ($summaryRow) {
        Convert-FmmCliInt64 $summaryRow.fmm_local_direct_fallback_pairs
    } else { 0 }
    $base.far_pairs_after = if ($summaryRow) {
        Convert-FmmCliInt64 $summaryRow.fmm_local_far_pairs_after
    } else { 0 }
    $base.plan_build_count =
        Get-Z4n3Int64 (Get-Z4n3Prop $profile 'fmm_local_plan_build_count' 0)
    $base.cache_build_count =
        Get-Z4n3Int64 (Get-Z4n3Prop $profile 'fmm_local_cache_build_count' 0)
    $base.plan_hash = [string](Get-Z4n3Prop $profile 'fmm_local_plan_hash' '')
    $base.linf_max_relative_tail =
        Get-Z4n3Double (Get-Z4n3Prop $profile 'fmm_local_linf_max_relative_tail' 0.0) 0.0
    $base.linf_outlier_panels =
        Get-Z4n3Int64 (Get-Z4n3Prop $profile 'fmm_local_linf_outlier_panels' 0)
    $base.geometry_rejected_cells =
        Get-Z4n3Int64 (Get-Z4n3Prop $profile 'fmm_local_geometry_guard_rejected_cells' 0)
    return [pscustomobject]$base
}

function Add-Z4n3Ratios([object[]]$Rows,[int64]$BaselineFallbackPairs) {
    foreach ($row in $Rows) {
        $ratio = if ($BaselineFallbackPairs -gt 0) {
            [double]$row.direct_fallback_pairs / [double]$BaselineFallbackPairs
        } else { [double]::NaN }
        $row.baseline_fallback_pairs = $BaselineFallbackPairs
        $row.fallback_ratio = $ratio
    }
    return $Rows
}

function Test-Z4n3BudgetCandidate([object]$Row) {
    return $Row.status -eq 'ok' -and
        (Test-FmmCliFinite ([double]$Row.candidate_audit_rms)) -and
        (Test-FmmCliFinite ([double]$Row.sampled_direct_true_residual)) -and
        [double]$Row.candidate_audit_rms -le 8e-3 -and
        [double]$Row.sampled_direct_true_residual -le 8e-3 -and
        (Test-FmmCliFinite ([double]$Row.fallback_ratio)) -and
        [double]$Row.fallback_ratio -le 3.0 -and
        [double]$Row.peak_rss_mb -lt 512.0 -and
        [int64]$Row.plan_build_count -eq 1 -and
        [int64]$Row.cache_build_count -eq 1
}

function Test-Z4n3ConfirmCandidate([object]$Row) {
    return $Row.status -eq 'ok' -and
        (Test-FmmCliFinite ([double]$Row.candidate_audit_rms)) -and
        (Test-FmmCliFinite ([double]$Row.sampled_direct_true_residual)) -and
        (Test-FmmCliFinite ([double]$Row.fallback_ratio)) -and
        [double]$Row.candidate_audit_rms -le 8e-3 -and
        [double]$Row.sampled_direct_true_residual -le 8e-3 -and
        [double]$Row.fallback_ratio -le 3.0 -and
        [double]$Row.peak_rss_mb -lt 512.0 -and
        [int64]$Row.plan_build_count -eq 1 -and
        [int64]$Row.cache_build_count -eq 1
}

function Write-Z4n3Svg([object[]]$Rows,[string]$Path) {
    $svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="1080" height="560" viewBox="0 0 1080 560">',
        '<rect width="1080" height="560" fill="white"/>',
        '<text x="32" y="38" font-family="Arial" font-size="20">Z4N.3 outlier refinement work-budget sweep</text>',
        '<line x1="86" y1="450" x2="1020" y2="450" stroke="#222"/>',
        '<line x1="86" y1="90" x2="86" y2="450" stroke="#222"/>',
        '<text x="92" y="82" font-family="Arial" font-size="12">audit (-log10), label shows fallback ratio</text>')
    $i = 0
    foreach ($row in $Rows) {
        $x = 116 + 70 * $i
        $audit = if (Test-FmmCliFinite ([double]$row.candidate_audit_rms)) {
            [double]$row.candidate_audit_rms
        } else { 1.0 }
        $bar = [Math]::Min(330.0, 62.0 * (-[Math]::Log10([Math]::Max(1e-8, $audit))))
        $color = if (Test-Z4n3BudgetCandidate $row) { '#16a34a' } elseif ([double]$row.fallback_ratio -le 3.0) { '#f59e0b' } else { '#3b82f6' }
        $svg += "<rect x='$x' y='$([int](450-$bar))' width='38' height='$([int]$bar)' fill='$color'/>"
        $svg += "<text x='$($x-8)' y='470' font-family='Arial' font-size='10' transform='rotate(45 $($x-8),470)'>$($row.condition)</text>"
        $label = if (Test-FmmCliFinite ([double]$row.fallback_ratio)) { '{0:G3}x' -f [double]$row.fallback_ratio } else { 'n/a' }
        $svg += "<text x='$($x-6)' y='$([int](440-$bar))' font-family='Arial' font-size='10'>$label</text>"
        ++$i
    }
    $svg += '</svg>'
    $svg | Set-Content -LiteralPath $Path -Encoding ASCII
}

if ($ConfirmCandidate) {
    if (-not (Test-Path -LiteralPath $budgetCsv)) {
        throw "Missing quick budget matrix: $budgetCsv. Rerun: test\pbe_fagbi_scope_z4n3_outlier_budget_sweep.ps1 -TitanExe .\titan.exe -Quick -CaseTimeoutSeconds 900 -AllowFailedDecision"
    }
    $quickRows = @(Import-Csv -LiteralPath $budgetCsv)
    $quickBaseline = @($quickRows | Where-Object {
        $_.fixture -eq '1a63' -and $_.condition -eq 'baseline_rms'
    } | Select-Object -First 1)
    $quickCandidate = @($quickRows | Where-Object {
        $_.fixture -eq '1a63' -and $_.condition -eq 'geom_theta0p5' -and
        $_.refine_mode -eq 'operator-tail-static' -and
        $_.geometry_guard -eq 'strict-geometry-mac'
    } | Select-Object -First 1)
    if (-not $quickBaseline.Count -or -not $quickCandidate.Count) {
        throw "Quick budget matrix must contain baseline_rms and geom_theta0p5 rows before -ConfirmCandidate. Rerun quick Z4N.3 first."
    }
    $baselineFallbackPairs = [int64]$quickBaseline[0].direct_fallback_pairs
    if ($baselineFallbackPairs -le 0) {
        throw "Quick baseline_rms has invalid direct_fallback_pairs=$baselineFallbackPairs."
    }

    $confirmRows = New-Object System.Collections.Generic.List[object]
    foreach ($fixtureName in @('1a63','barnase','1sh1')) {
        $fixture = Get-Z4n3Fixture $fixtureName
        $confirmRows.Add((Invoke-Z4n3Case $fixture 'confirm_geom_theta0p5' `
            'operator-tail-static' 0.002 'strict-geometry-mac' 0.5)) | Out-Null
    }
    $confirmMatrix = @(Add-Z4n3Ratios @($confirmRows.ToArray()) $baselineFallbackPairs)
    $timeouts = @($confirmMatrix | Where-Object status -eq 'timeout')
    $failed = @($confirmMatrix | Where-Object status -eq 'failed')
    $lifecycleGaps = @($confirmMatrix | Where-Object {
        $_.status -eq 'ok' -and ([int64]$_.plan_build_count -ne 1 -or [int64]$_.cache_build_count -ne 1)
    })
    $passRows = @($confirmMatrix | Where-Object { Test-Z4n3ConfirmCandidate $_ })
    $oneA63Pass = @($passRows | Where-Object fixture -eq '1a63').Count -eq 1
    $accuracyGaps = @($confirmMatrix | Where-Object {
        $_.status -eq 'ok' -and
        ((Test-FmmCliFinite ([double]$_.candidate_audit_rms)) -and
         (Test-FmmCliFinite ([double]$_.sampled_direct_true_residual))) -and
        ([double]$_.candidate_audit_rms -gt 8e-3 -or
         [double]$_.sampled_direct_true_residual -gt 8e-3)
    })
    $budgetGaps = @($confirmMatrix | Where-Object {
        $_.status -eq 'ok' -and
        [double]$_.candidate_audit_rms -le 8e-3 -and
        [double]$_.sampled_direct_true_residual -le 8e-3 -and
        ((Test-FmmCliFinite ([double]$_.fallback_ratio)) -and
         ([double]$_.fallback_ratio -gt 3.0 -or [double]$_.peak_rss_mb -ge 512.0))
    })

    $confirmDecision = 'z4n3_cross_fixture_accuracy_gap'
    if ($timeouts.Count) {
        $confirmDecision = 'z4n3_runtime_limited'
    } elseif ($failed.Count -or $lifecycleGaps.Count) {
        $confirmDecision = 'z4n3_parser_or_lifecycle_gap'
    } elseif ($passRows.Count -eq 3) {
        $confirmDecision = 'z4n3_confirmed_budget_policy'
    } elseif ($oneA63Pass) {
        $confirmDecision = 'z4n3_1a63_only_budget_policy'
    } elseif ($budgetGaps.Count) {
        $confirmDecision = 'z4n3_cross_fixture_budget_gap'
    } elseif ($accuracyGaps.Count) {
        $confirmDecision = 'z4n3_cross_fixture_accuracy_gap'
    }

    $confirmMatrix | Export-Csv -LiteralPath $confirmMatrixCsv -NoTypeInformation
    [pscustomobject]@{
        decision=$confirmDecision;rows=$confirmMatrix.Count
        baseline_fallback_pairs=$baselineFallbackPairs
        passed_rows=$passRows.Count
        one_a63_pass=$oneA63Pass
    } | Export-Csv -LiteralPath $confirmDecisionCsv -NoTypeInformation
    @(
        'TITAN FMM/M2L Scope Z4N.3 Outlier Budget Sweep',
        "Decision: $confirmDecision",
        "Mode: ConfirmCandidate",
        "Rows: $($confirmMatrix.Count)",
        "Candidate: operator-tail-static tol=0.002 strict-geometry-mac theta=0.5",
        "Baseline fallback pairs: $baselineFallbackPairs",
        '',
        'Artifacts:',
        $matrixCsv,
        $budgetCsv,
        $confirmMatrixCsv,
        $confirmDecisionCsv,
        $decisionCsv,
        $summaryPath,
        $plotPath
    ) | Set-Content -LiteralPath $rootReport -Encoding ASCII
    Write-Host "Z4N.3 confirm decision: $confirmDecision"
    if (-not $AllowFailedDecision -and
        $confirmDecision -notin @('z4n3_confirmed_budget_policy',
            'z4n3_1a63_only_budget_policy','z4n3_cross_fixture_accuracy_gap',
            'z4n3_cross_fixture_budget_gap')) {
        throw "Z4N.3 confirm failed with decision $confirmDecision"
    }
    return
}

$rows = New-Object System.Collections.Generic.List[object]
$oneA63 = Get-Z4n3Fixture '1a63'

$staged = @(
    @{name='baseline_rms';refine='operator-tail-static';tol=0.002;guard='off';theta=0.0}
)
foreach ($tol in @(0.004,0.008,0.016,0.032)) {
    $staged += @{name="linf_tol$(Convert-Z4n3ConditionName $tol)";refine='operator-tail-linf';tol=$tol;guard='off';theta=0.0}
}
foreach ($theta in @(0.5,0.75,1.0)) {
    $staged += @{name="geom_theta$(Convert-Z4n3ConditionName $theta)";refine='operator-tail-static';tol=0.002;guard='strict-geometry-mac';theta=$theta}
}
foreach ($condition in $staged) {
    $rows.Add((Invoke-Z4n3Case $oneA63 $condition.name $condition.refine `
        $condition.tol $condition.guard $condition.theta)) | Out-Null
}

$linfCompleted = @($rows.ToArray() | Where-Object {
    $_.fixture -eq '1a63' -and $_.refine_mode -eq 'operator-tail-linf' -and
    $_.geometry_guard -eq 'off' -and $_.status -eq 'ok'
} | Sort-Object @{Expression='candidate_audit_rms';Ascending=$true},
    @{Expression='direct_fallback_pairs';Ascending=$true} | Select-Object -First 2)
foreach ($linf in $linfCompleted) {
    foreach ($theta in @(0.75,1.0)) {
        $tol = [double]$linf.refine_tolerance
        $name = "combined_linf$(Convert-Z4n3ConditionName $tol)_geom$(Convert-Z4n3ConditionName $theta)"
        $rows.Add((Invoke-Z4n3Case $oneA63 $name 'operator-tail-linf' `
            $tol 'strict-geometry-mac' $theta)) | Out-Null
    }
}

$matrix = @($rows.ToArray())
$baseline = @($matrix | Where-Object {
    $_.fixture -eq '1a63' -and $_.condition -eq 'baseline_rms'
} | Select-Object -First 1)
$baselineFallbackPairs = if ($baseline.Count) { [int64]$baseline[0].direct_fallback_pairs } else { 0 }
$matrix = @(Add-Z4n3Ratios $matrix $baselineFallbackPairs)

$oneA63Candidates = @($matrix | Where-Object { $_.fixture -eq '1a63' -and (Test-Z4n3BudgetCandidate $_) } |
    Sort-Object @{Expression='candidate_audit_rms';Ascending=$true},
        @{Expression='fallback_ratio';Ascending=$true})
$bestBudget = if ($oneA63Candidates.Count) { $oneA63Candidates[0] } else { $null }

if (-not $Quick -and $null -ne $bestBudget) {
    foreach ($fixtureName in @('barnase','1sh1')) {
        $fixture = Get-Z4n3Fixture $fixtureName
        $bestTol = [double]$bestBudget.refine_tolerance
        $bestTheta = [double]$bestBudget.geometry_theta
        $rows.Add((Invoke-Z4n3Case $fixture "best_$($bestBudget.condition)" `
            $bestBudget.refine_mode $bestTol $bestBudget.geometry_guard $bestTheta)) | Out-Null
    }
    $matrix = @(Add-Z4n3Ratios @($rows.ToArray()) $baselineFallbackPairs)
}

$timeouts = @($matrix | Where-Object status -eq 'timeout')
$failed = @($matrix | Where-Object status -eq 'failed')
$lifecycleGaps = @($matrix | Where-Object {
    $_.status -eq 'ok' -and ([int64]$_.plan_build_count -ne 1 -or [int64]$_.cache_build_count -ne 1)
})
$accurateRows = @($matrix | Where-Object {
    $_.fixture -eq '1a63' -and $_.status -eq 'ok' -and
    (Test-FmmCliFinite ([double]$_.candidate_audit_rms)) -and
    (Test-FmmCliFinite ([double]$_.sampled_direct_true_residual)) -and
    [double]$_.candidate_audit_rms -le 8e-3 -and
    [double]$_.sampled_direct_true_residual -le 8e-3
})
$budgetRows = @($matrix | Where-Object {
    $_.fixture -eq '1a63' -and $_.status -eq 'ok' -and
    (Test-FmmCliFinite ([double]$_.fallback_ratio)) -and [double]$_.fallback_ratio -le 3.0
})
$best = @($matrix | Where-Object { $_.fixture -eq '1a63' -and $_.status -eq 'ok' } |
    Sort-Object @{Expression='candidate_audit_rms';Ascending=$true},
        @{Expression='fallback_ratio';Ascending=$true} | Select-Object -First 1)

$decision = 'z4n3_no_accuracy_budget'
if ($timeouts.Count) {
    $decision = 'z4n3_runtime_limited'
} elseif ($failed.Count -or $lifecycleGaps.Count) {
    $decision = 'z4n3_parser_or_lifecycle_gap'
} elseif ($oneA63Candidates.Count) {
    $decision = 'z4n3_budget_candidate_found'
} elseif ($accurateRows.Count) {
    $decision = 'z4n3_accuracy_only_costly'
} elseif ($budgetRows.Count) {
    $decision = 'z4n3_no_accuracy_budget'
}

$matrix | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$matrix | Select-Object fixture,condition,refine_mode,refine_tolerance,
    geometry_guard,geometry_theta,candidate_audit_rms,
    sampled_direct_true_residual,direct_fallback_pairs,
    baseline_fallback_pairs,fallback_ratio,far_pairs_after,peak_rss_mb,
    plan_build_count,cache_build_count,plan_hash,linf_max_relative_tail,
    linf_outlier_panels,geometry_rejected_cells |
    Export-Csv -LiteralPath $budgetCsv -NoTypeInformation
[pscustomobject]@{
    decision=$decision;quick=[bool]$Quick;rows=$matrix.Count
    best_condition=if ($best.Count) { $best[0].condition } else { '' }
    best_audit=if ($best.Count) { $best[0].candidate_audit_rms } else { [double]::NaN }
    best_residual=if ($best.Count) { $best[0].sampled_direct_true_residual } else { [double]::NaN }
    best_fallback_ratio=if ($best.Count) { $best[0].fallback_ratio } else { [double]::NaN }
    budget_candidate=if ($null -ne $bestBudget) { $bestBudget.condition } else { '' }
} | Export-Csv -LiteralPath $decisionCsv -NoTypeInformation

@('scope: Z4N.3 outlier refinement work-budget sweep',
  "decision: $decision","quick: $([bool]$Quick)","rows: $($matrix.Count)",
  "best_condition: $(if ($best.Count) { $best[0].condition } else { '' })",
  "best_audit: $(if ($best.Count) { $best[0].candidate_audit_rms } else { [double]::NaN })",
  "best_fallback_ratio: $(if ($best.Count) { $best[0].fallback_ratio } else { [double]::NaN })",
  'production_default: grid','gmres_p4: explicit diagnostic only') |
    Set-Content -LiteralPath $summaryPath -Encoding ASCII

Write-Z4n3Svg $matrix $plotPath

@(
    'TITAN FMM/M2L Scope Z4N.3 Outlier Budget Sweep',
    "Decision: $decision",
    "Rows: $($matrix.Count)",
    "Best condition: $(if ($best.Count) { $best[0].condition } else { '' })",
    "Budget candidate: $(if ($null -ne $bestBudget) { $bestBudget.condition } else { '' })",
    '',
    'Artifacts:',
    $matrixCsv,
    $budgetCsv,
    $decisionCsv,
    $summaryPath,
    $plotPath
) | Set-Content -LiteralPath $rootReport -Encoding ASCII

if (-not $AllowFailedDecision -and $decision -notin @('z4n3_budget_candidate_found',
        'z4n3_accuracy_only_costly','z4n3_no_accuracy_budget')) {
    throw "Z4N.3 budget sweep failed with decision $decision"
}

Write-Host "Z4N.3 budget sweep decision: $decision"
