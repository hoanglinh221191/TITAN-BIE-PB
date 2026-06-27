param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$PrimaryFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$ExpansionFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public_z4q'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4w'),
    [string[]]$Fixtures = @(),
    [string[]]$Stages = @(),
    [int]$CaseTimeoutSeconds = 1800,
    [int]$Threads = 4,
    [switch]$Quick,
    [switch]$PlanOnly,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$PrimaryFixtureRoot = [IO.Path]::GetFullPath($PrimaryFixtureRoot)
$ExpansionFixtureRoot = [IO.Path]::GetFullPath($ExpansionFixtureRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4w_matrix.csv'
$TraceCsv = Join-Path $WorkRoot 'fagbi_scope_z4w_residual_trace.csv'
$SensitivityCsv = Join-Path $WorkRoot 'fagbi_scope_z4w_sensitivity_summary.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4w_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4w_summary.dat'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4W_ABLATION_SENSITIVITY_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$RuntimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
        $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
    }
    $RuntimeRepoRoot = Split-Path -Parent $commonGitDir
}
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    throw "Missing lib_ff under runtime repo root: $RuntimeRepoRoot"
}

function Get-Z4wDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or $text -eq 'nan') { return $Default }
    $parsed = 0.0
    if ([double]::TryParse($text, [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) { return $parsed }
    return $Default
}

function Get-Z4wInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $parsed = [int64]0
    if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) { return $parsed }
    return $Default
}

function Get-Z4wProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Import-Z4wManifestRows([string]$Root, [string]$ManifestName) {
    $path = Join-Path $Root $ManifestName
    if (-not (Test-Path -LiteralPath $path)) { return @() }
    return @(Import-Csv -LiteralPath $path | ForEach-Object {
        $_ | Add-Member -NotePropertyName fixture_root -NotePropertyValue $Root -Force
        $_
    })
}

function Get-Z4wFixtureRows {
    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($row in (Import-Z4wManifestRows $PrimaryFixtureRoot 'public_fixture_manifest.csv')) { $rows.Add($row) | Out-Null }
    foreach ($row in (Import-Z4wManifestRows $ExpansionFixtureRoot 'public_fixture_expansion_manifest.csv')) { $rows.Add($row) | Out-Null }
    return $rows.ToArray()
}

function Get-Z4wFixture([string]$Name, [object[]]$Rows) {
    $match = @($Rows | Where-Object { $_.fixture -eq $Name })
    if ($match.Count -ne 1) { throw "Expected exactly one cached public fixture row for '$Name', found $($match.Count)." }
    $row = $match[0]
    $pqrPath = Join-Path $row.fixture_root $row.pqr
    if (-not (Test-Path -LiteralPath $pqrPath)) { throw "Missing cached PQR for $Name`: $pqrPath" }
    $row | Add-Member -NotePropertyName pqr_path -NotePropertyValue $pqrPath -Force
    return $row
}

function New-Z4wStage(
    [int]$Order,
    [string]$Stage,
    [string]$Axis,
    [string]$Value,
    [int]$TaylorOrder = 4,
    [double]$OpeningTheta = 0.6,
    [string]$Refine = 'operator-tail-static',
    [double]$RefineTol = 0.002,
    [string]$GeometryGuard = 'strict-geometry-mac',
    [double]$GeometryTheta = 0.5,
    [int]$IluCap = 64,
    [string]$UnsupportedReason = ''
) {
    return [pscustomobject]@{
        order=$Order; stage=$Stage; axis=$Axis; value=$Value;
        taylor_order=$TaylorOrder; opening_theta=$OpeningTheta;
        refine=$Refine; refine_tol=$RefineTol;
        geometry_guard=$GeometryGuard; geometry_theta=$GeometryTheta;
        ilu_cap=$IluCap; unsupported_reason=$UnsupportedReason;
    }
}

function Get-Z4wStages {
    return @(
        (New-Z4wStage 1 'baseline_guarded_p4' 'baseline' 'p4 theta0.6 tail-on geometry-on ilu64'),
        (New-Z4wStage 2 'p_order2' 'expansion_order' '2' -TaylorOrder 2),
        (New-Z4wStage 3 'p_order3' 'expansion_order' '3' -TaylorOrder 3),
        (New-Z4wStage 4 'p_order4' 'expansion_order' '4' -TaylorOrder 4),
        (New-Z4wStage 5 'p_order5' 'expansion_order' '5' -TaylorOrder 5 -UnsupportedReason 'taylor-order=5 is reserved for future guarded local FMM/M2L work'),
        (New-Z4wStage 6 'theta0p4' 'opening_theta' '0.4' -OpeningTheta 0.4),
        (New-Z4wStage 7 'theta0p5' 'opening_theta' '0.5' -OpeningTheta 0.5),
        (New-Z4wStage 8 'theta0p6' 'opening_theta' '0.6' -OpeningTheta 0.6),
        (New-Z4wStage 9 'theta0p7' 'opening_theta' '0.7' -OpeningTheta 0.7),
        (New-Z4wStage 10 'tail_off' 'tail_indicator' 'off' -Refine off),
        (New-Z4wStage 11 'tail_on' 'tail_indicator' 'operator-tail-static tol0.002' -Refine operator-tail-static -RefineTol 0.002),
        (New-Z4wStage 12 'geometry_off' 'geometry_guard' 'off' -GeometryGuard off -GeometryTheta 0.0),
        (New-Z4wStage 13 'geometry_on' 'geometry_guard' 'strict theta0.5' -GeometryGuard strict-geometry-mac -GeometryTheta 0.5),
        (New-Z4wStage 14 'ilu_cap32' 'ilu_neighbor_cap' '32' -IluCap 32),
        (New-Z4wStage 15 'ilu_cap64' 'ilu_neighbor_cap' '64' -IluCap 64),
        (New-Z4wStage 16 'ilu_cap128' 'ilu_neighbor_cap' '128' -IluCap 128)
    )
}

function New-Z4wRow([object]$Fixture, [object]$Stage, [string]$RunDir) {
    return [pscustomobject]@{
        fixture=$Fixture.fixture; pdb_id=$Fixture.pdb_id; protein_name=$Fixture.protein_name;
        atoms=(Get-Z4wInt64 $Fixture.atom_count); stage=$Stage.stage; stage_order=$Stage.order;
        axis=$Stage.axis; value=$Stage.value; taylor_order=$Stage.taylor_order;
        opening_theta=$Stage.opening_theta; refine=$Stage.refine; refine_tol=$Stage.refine_tol;
        geometry_guard=$Stage.geometry_guard; geometry_theta=$Stage.geometry_theta; ilu_neighbor_cap=$Stage.ilu_cap;
        exit_code=0; timed_out='False'; status='planned'; failure_class='planned';
        panels=0; iterations=0; energy_kj_mol=[double]::NaN;
        final_residual=[double]::NaN; true_residual=[double]::NaN; sampled_audit_rms=[double]::NaN;
        sampled_direct_true_residual=[double]::NaN; matvec_seconds_avg=[double]::NaN; wall_seconds=[double]::NaN;
        peak_rss_mb=[double]::NaN; preconditioner_setup_seconds=[double]::NaN; preconditioner_apply_seconds=[double]::NaN;
        preconditioner_sparse_nnz_blocks=0; fmm_local_prepare_seconds=[double]::NaN;
        fmm_local_gmres_apply_seconds_mean=[double]::NaN; fmm_local_moment_seconds_total=[double]::NaN;
        fmm_local_m2l_seconds_total=[double]::NaN; fmm_local_fallback_seconds_total=[double]::NaN;
        fmm_local_l2p_seconds_total=[double]::NaN; matvec_audit_seconds=[double]::NaN;
        plan_build_count=0; cache_build_count=0; gmres_matvec_calls=0;
        geometry_guard_direct_fallback_pairs=0; geometry_guard_rejected_cells=0;
        linf_max_abs_tail=[double]::NaN; linf_outlier_panels=0;
        run_dir=$RunDir
    }
}

function Write-Z4wArg([object]$Fixture, [object]$Stage, [string]$RunDir) {
    return Write-FmmCliArgFile -RunDir $RunDir -RepoRoot $RuntimeRepoRoot `
        -PqrPath $Fixture.pqr_path -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -SurfacePoints 8 -MaxIter 720 -Restart 128 `
        -Tolerance 1e-6 -Threads $Threads -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta $Stage.opening_theta -AuditSampleCount 512 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerBlockPanelCap 32 `
        -PreconditionerOverlapPanelCap 64 -PreconditionerIluNeighborCap $Stage.ilu_cap `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder $Stage.taylor_order -FmmCartesianAudit off -FmmLocal gmres-p4 `
        -FmmLocalRefine $Stage.refine -FmmLocalRefineTolerance $Stage.refine_tol `
        -FmmLocalCache green-jet-compressed -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback packed-avx2 -FmmLocalPlanBuilder metric-cache `
        -FmmLocalRefineProbeCount 4 -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard $Stage.geometry_guard -FmmLocalGeometryTheta $Stage.geometry_theta
}

function Read-Z4wOutputs([string]$RunDir, [object]$Row) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (Test-Path -LiteralPath $energyPath) {
        $energy = @(Import-Csv -LiteralPath $energyPath)[0]
        $Row.energy_kj_mol = Get-Z4wDouble $energy.E_rxn_kJmol
        if ([string]$energy.status -ne 'ok') { $Row.status = [string]$energy.status; $Row.failure_class = 'energy_status_gap' }
    }
    if (-not (Test-Path -LiteralPath $profilePath)) {
        $Row.failure_class = 'missing_profile'
        return $Row
    }
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    if ($Row.status -eq 'ok' -and [string]$profile.status -ne 'ok') {
        $Row.status = "profile=$($profile.status)"
        $Row.failure_class = 'profile_status_gap'
    }
    $Row.panels = Get-Z4wInt64 $profile.panels
    $Row.iterations = Get-Z4wInt64 $profile.iterations
    $Row.final_residual = Get-Z4wDouble $profile.final_residual
    $Row.true_residual = Get-Z4wDouble $profile.true_residual_rel
    $Row.sampled_audit_rms = Get-Z4wDouble $profile.tree_vs_direct_matvec_rel
    $Row.sampled_direct_true_residual = Get-Z4wDouble (Get-Z4wProfileValue $profile 'fmm_local_sampled_direct_true_residual' '')
    $Row.matvec_seconds_avg = Get-Z4wDouble $profile.matvec_seconds_avg
    $Row.matvec_audit_seconds = Get-Z4wDouble $profile.matvec_audit_seconds
    $Row.preconditioner_setup_seconds = Get-Z4wDouble (Get-Z4wProfileValue $profile 'preconditioner_setup_seconds' '')
    $Row.preconditioner_apply_seconds = Get-Z4wDouble (Get-Z4wProfileValue $profile 'preconditioner_apply_seconds' '')
    $Row.preconditioner_sparse_nnz_blocks = Get-Z4wInt64 (Get-Z4wProfileValue $profile 'preconditioner_sparse_nnz_blocks' 0)
    $Row.fmm_local_prepare_seconds = Get-Z4wDouble (Get-Z4wProfileValue $profile 'fmm_local_prepare_seconds' '')
    $Row.fmm_local_gmres_apply_seconds_mean = Get-Z4wDouble (Get-Z4wProfileValue $profile 'fmm_local_gmres_apply_seconds_mean' '')
    $Row.fmm_local_moment_seconds_total = Get-Z4wDouble (Get-Z4wProfileValue $profile 'fmm_local_moment_seconds_total' '')
    $Row.fmm_local_m2l_seconds_total = Get-Z4wDouble (Get-Z4wProfileValue $profile 'fmm_local_m2l_seconds_total' '')
    $Row.fmm_local_fallback_seconds_total = Get-Z4wDouble (Get-Z4wProfileValue $profile 'fmm_local_fallback_seconds_total' '')
    $Row.fmm_local_l2p_seconds_total = Get-Z4wDouble (Get-Z4wProfileValue $profile 'fmm_local_l2p_seconds_total' '')
    $Row.plan_build_count = Get-Z4wInt64 (Get-Z4wProfileValue $profile 'fmm_local_plan_build_count' 0)
    $Row.cache_build_count = Get-Z4wInt64 (Get-Z4wProfileValue $profile 'fmm_local_cache_build_count' 0)
    $Row.gmres_matvec_calls = Get-Z4wInt64 (Get-Z4wProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
    $Row.geometry_guard_direct_fallback_pairs = Get-Z4wInt64 (Get-Z4wProfileValue $profile 'fmm_local_geometry_guard_direct_fallback_pairs' 0)
    $Row.geometry_guard_rejected_cells = Get-Z4wInt64 (Get-Z4wProfileValue $profile 'fmm_local_geometry_guard_rejected_cells' 0)
    $Row.linf_max_abs_tail = Get-Z4wDouble (Get-Z4wProfileValue $profile 'fmm_local_linf_max_abs_tail' '')
    $Row.linf_outlier_panels = Get-Z4wInt64 (Get-Z4wProfileValue $profile 'fmm_local_linf_outlier_panels' 0)
    if ($Row.status -eq 'ok') {
        if (-not (Test-FmmCliFinite ([double]$Row.sampled_audit_rms))) {
            $Row.failure_class = 'missing_sampled_audit'
        } elseif (-not (Test-FmmCliFinite ([double]$Row.sampled_direct_true_residual))) {
            $Row.failure_class = 'missing_sampled_direct_true_residual'
        } elseif ([double]$Row.sampled_audit_rms -le 8.0e-3 -and [double]$Row.sampled_direct_true_residual -le 8.0e-3) {
            $Row.failure_class = 'pass'
        } else {
            $Row.failure_class = 'operator_accuracy_gap'
        }
    }
    return $Row
}

function Test-Z4wRowFinal([object]$Row) {
    return $Row.status -eq 'ok' -and $Row.failure_class -eq 'pass'
}

function Test-Z4wShouldRun([string]$RowCsv) {
    if ($Force) { return $true }
    if (-not (Test-Path -LiteralPath $RowCsv)) { return $true }
    try {
        $row = @(Import-Csv -LiteralPath $RowCsv)[0]
        if ($row.status -eq 'planned') { return $true }
        if (-not $RerunFailed) { return $false }
        return -not (Test-Z4wRowFinal $row)
    } catch {
        return $true
    }
}

function Invoke-Z4wCase([object]$Fixture, [object]$Stage) {
    $runName = '{0}_{1}' -f $Fixture.fixture, $Stage.stage
    $runDir = Join-Path $WorkRoot $runName
    $rowCsv = Join-Path $runDir 'fagbi_scope_z4w_row.csv'
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null

    if (-not (Test-Z4wShouldRun $rowCsv)) {
        Write-Host "resume: $runName"
        return @(Import-Csv -LiteralPath $rowCsv)[0]
    }

    $row = New-Z4wRow -Fixture $Fixture -Stage $Stage -RunDir $runDir
    if (-not [string]::IsNullOrWhiteSpace($Stage.unsupported_reason)) {
        Write-Host "skip: $runName ($($Stage.unsupported_reason))"
        $row.status = 'unsupported'
        $row.failure_class = 'unsupported_by_current_parser'
        $row.exit_code = 0
        $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
        return $row
    }

    Write-Host "run: $runName (timeout ${CaseTimeoutSeconds}s)"
    $argPath = Write-Z4wArg -Fixture $Fixture -Stage $Stage -RunDir $runDir
    $bare = @(Get-Content -LiteralPath $argPath | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) { throw "Generated arg.in contains bare options: $($bare -join ', ')" }
    if ($PlanOnly) {
        $row.status = 'planned'
        $row.failure_class = 'plan_only'
        $row.wall_seconds = 0.0
        $row.peak_rss_mb = 0.0
        $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
        return $row
    }

    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $row.exit_code = $run.exit_code
    $row.timed_out = [string]$run.timed_out
    $row.wall_seconds = $run.wall_seconds
    $row.peak_rss_mb = $run.peak_working_set_mb
    if ($run.timed_out) {
        $row.status = 'timeout'
        $row.failure_class = 'runtime_limited'
    } elseif ($run.exit_code -ne 0) {
        $row.status = 'failed'
        $row.failure_class = "exit_code_$($run.exit_code)"
    } else {
        $row.status = 'ok'
        $row.failure_class = 'pass'
        $row = Read-Z4wOutputs -RunDir $runDir -Row $row
    }
    $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
    return $row
}

function Export-Z4wResidualTrace([object[]]$Rows, [string]$Path) {
    $trace = New-Object System.Collections.Generic.List[object]
    foreach ($row in $Rows) {
        $profilePath = Join-Path $row.run_dir 'pbe_tabi_profile.csv'
        if (-not (Test-Path -LiteralPath $profilePath)) { continue }
        $profile = @(Import-Csv -LiteralPath $profilePath)[0]
        $hist = [string](Get-Z4wProfileValue $profile 'residual_history' '')
        if ([string]::IsNullOrWhiteSpace($hist)) { continue }
        $values = @($hist -split ';' | Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
        for ($i = 0; $i -lt $values.Count; ++$i) {
            $res = Get-Z4wDouble $values[$i]
            $trace.Add([pscustomobject]@{
                fixture=$row.fixture; pdb_id=$row.pdb_id; stage=$row.stage; axis=$row.axis;
                iteration_index=$i; residual=$res;
                log10_residual=if (Test-FmmCliFinite $res -and $res -gt 0) { [Math]::Log10($res) } else { [double]::NaN }
            }) | Out-Null
        }
    }
    $trace.ToArray() | Export-Csv -LiteralPath $Path -NoTypeInformation
}

function Write-Z4wSensitivitySummary([object[]]$Rows, [string]$Path) {
    $summary = New-Object System.Collections.Generic.List[object]
    foreach ($axis in @($Rows.axis | Sort-Object -Unique)) {
        $set = @($Rows | Where-Object { $_.axis -eq $axis })
        $pass = @($set | Where-Object { $_.failure_class -eq 'pass' })
        $summary.Add([pscustomobject]@{
            axis=$axis; rows=$set.Count; pass_rows=$pass.Count;
            best_audit=if ($pass.Count) { ($pass | ForEach-Object { Get-Z4wDouble $_.sampled_audit_rms } | Measure-Object -Minimum).Minimum } else { [double]::NaN };
            min_wall_seconds=if ($pass.Count) { ($pass | ForEach-Object { Get-Z4wDouble $_.wall_seconds } | Measure-Object -Minimum).Minimum } else { [double]::NaN };
            interpretation=if ($axis -eq 'expansion_order') { 'p=2 and p=5 are recorded as unsupported unless the parser is extended' } else { 'measured sensitivity axis' }
        }) | Out-Null
    }
    $summary.ToArray() | Export-Csv -LiteralPath $Path -NoTypeInformation
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$fixtureRows = Get-Z4wFixtureRows
$defaultFixtures = if ($Quick) {
    @('1a63')
} elseif ($Fixtures.Count -gt 0) {
    $Fixtures
} else {
    @('barnase','1a63','1sh1','1ubq','2ci2','1vii')
}
$allStages = Get-Z4wStages
$selectedStages = if ($Stages.Count -gt 0) {
    @($allStages | Where-Object { $Stages -contains $_.stage })
} elseif ($Quick) {
    @($allStages | Where-Object { @('baseline_guarded_p4','theta0p5','theta0p6','geometry_off','geometry_on','ilu_cap64') -contains $_.stage })
} else {
    $allStages
}
if ($selectedStages.Count -eq 0) { throw "No matching Z4W stages selected. Valid stages: $(@($allStages.stage) -join ', ')" }

$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $defaultFixtures) {
    $fixture = Get-Z4wFixture -Name $fixtureName -Rows $fixtureRows
    foreach ($stage in $selectedStages) {
        $rows.Add((Invoke-Z4wCase -Fixture $fixture -Stage $stage)) | Out-Null
    }
}

$allRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $defaultFixtures) {
    foreach ($stage in $allStages) {
        $rowCsv = Join-Path (Join-Path $WorkRoot ('{0}_{1}' -f $fixtureName, $stage.stage)) 'fagbi_scope_z4w_row.csv'
        if (Test-Path -LiteralPath $rowCsv) { $allRows.Add((@(Import-Csv -LiteralPath $rowCsv)[0])) | Out-Null }
    }
}
$matrix = @($allRows.ToArray() | Sort-Object fixture,stage_order)
$matrix | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
Export-Z4wResidualTrace -Rows $matrix -Path $TraceCsv
Write-Z4wSensitivitySummary -Rows $matrix -Path $SensitivityCsv

$measured = @($matrix | Where-Object { $_.status -ne 'unsupported' -and $_.status -ne 'planned' })
$passes = @($matrix | Where-Object failure_class -eq 'pass')
$timeouts = @($matrix | Where-Object status -eq 'timeout')
$accuracyGaps = @($matrix | Where-Object failure_class -eq 'operator_accuracy_gap')
$decision = if ($PlanOnly) {
    'z4w_plan_only'
} elseif ($timeouts.Count -gt 0) {
    'z4w_runtime_limited'
} elseif ($accuracyGaps.Count -gt 0 -and $passes.Count -gt 0) {
    'z4w_sensitivity_ready_with_accuracy_gaps'
} elseif ($measured.Count -gt 0 -and $passes.Count -eq $measured.Count) {
    'z4w_sensitivity_pass'
} elseif ($measured.Count -gt 0) {
    'z4w_sensitivity_ready_with_gaps'
} else {
    'z4w_no_measured_rows'
}

[pscustomobject]@{
    decision=$decision; quick=[string][bool]$Quick; fixtures=($defaultFixtures -join ';');
    selected_stages=(@($selectedStages.stage) -join ';'); rows=$matrix.Count;
    measured_rows=$measured.Count; pass_rows=$passes.Count; timeout_rows=$timeouts.Count;
    accuracy_gap_rows=$accuracyGaps.Count; case_timeout_seconds=$CaseTimeoutSeconds
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
    'scope: Z4W ablation and sensitivity matrix',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "fixtures: $($defaultFixtures -join ', ')",
    "selected_stages: $(@($selectedStages.stage) -join ', ')",
    "rows: $($matrix.Count)",
    "case_timeout_seconds: $CaseTimeoutSeconds",
    'resume_policy: row CSV files under each run directory are reused unless -Force or -RerunFailed is passed',
    'outputs: matrix, residual trace, sensitivity summary, decision, root report',
    'production_default: grid',
    'tabi_treecode_default: degree-one; gmres-p4 explicit only'
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z4W Ablation Sensitivity Report',
    '',
    "Decision: $decision",
    '',
    "Matrix: $MatrixCsv",
    "Residual trace: $TraceCsv",
    "Sensitivity summary: $SensitivityCsv",
    '',
    'Resume: rerun the same command to continue missing rows; use -RerunFailed for failed/timeout rows or -Force for all rows.',
    '',
    'Unsupported p=2 and p=5 rows are recorded explicitly because the current CLI accepts Taylor order 3 or 4 only.'
)
$report | Set-Content -LiteralPath $RootReport -Encoding ASCII

Write-Host "Z4W decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Resume: rerun the same command to continue missing rows; use -RerunFailed for failed/timeout rows or -Force for all rows."

if (-not $AllowFailedDecision -and $decision -notin @('z4w_sensitivity_pass','z4w_plan_only')) {
    throw "Z4W decision did not pass: $decision"
}
