param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 2400,
    [int]$Threads = 4,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_t')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fmm_m2l_scope_t_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_t_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_t_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_T_CAP256_CONVERGENCE_REPORT.txt'

$AtomCap = 256
$SurfacePoints = 8
$OpeningTheta = 0.6
$AuditSampleCount = if ($Quick) { 128 } else { 512 }
$TemperatureK = 298.15
$AuditGate = 1.0e-2

if ($Quick) {
    $CaseTimeoutSeconds = [Math]::Min($CaseTimeoutSeconds, 600)
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope T work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function New-ScopeTEmptyCase([string]$Status) {
    return [pscustomobject]@{
        status=$Status; panels=0; iterations=0; final_residual=[double]::NaN
        true_residual=[double]::NaN; energy_kj_mol=[double]::NaN
        matvec_seconds_avg=[double]::NaN; total_seconds=[double]::NaN
        fmm_m2l_mode=''; fmm_near_pairs=0; fmm_far_pairs=0
        fmm_setup_seconds=0.0; fmm_apply_seconds=0.0
        fmm_finite_outputs='false'; tree_vs_direct_matvec_rel=[double]::NaN
        matvec_audit_rows=0; matvec_audit_seconds=0.0
        reaction_energy_source=''; timing_backend=''; memory_bytes=[double]::NaN
    }
}

function Test-ScopeTFinite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function Get-ScopeTConditions {
    if ($Quick) {
        return @(
            [pscustomobject]@{ name='baseline_240_r64_1e-6'; max_iter=240; restart=64; tolerance=1.0e-6 },
            [pscustomobject]@{ name='practical_360_r96_1e-4'; max_iter=360; restart=96; tolerance=1.0e-4 }
        )
    }
    return @(
        [pscustomobject]@{ name='baseline_240_r64_1e-6'; max_iter=240; restart=64; tolerance=1.0e-6 },
        [pscustomobject]@{ name='iter_360_r64_1e-6'; max_iter=360; restart=64; tolerance=1.0e-6 },
        [pscustomobject]@{ name='iter_480_r64_1e-6'; max_iter=480; restart=64; tolerance=1.0e-6 },
        [pscustomobject]@{ name='restart_480_r96_1e-6'; max_iter=480; restart=96; tolerance=1.0e-6 },
        [pscustomobject]@{ name='restart_480_r128_1e-6'; max_iter=480; restart=128; tolerance=1.0e-6 },
        [pscustomobject]@{ name='long_720_r128_1e-6'; max_iter=720; restart=128; tolerance=1.0e-6 },
        [pscustomobject]@{ name='practical_360_r96_1e-4'; max_iter=360; restart=96; tolerance=1.0e-4 },
        [pscustomobject]@{ name='practical_480_r128_1e-4'; max_iter=480; restart=128; tolerance=1.0e-4 }
    )
}

function Invoke-ScopeTCase(
    [string]$CaseRoot,
    [string]$PqrPath,
    [ValidateSet('direct','fmm')][string]$Backend,
    [int]$MaxIter,
    [int]$Restart,
    [double]$Tolerance,
    [int]$AuditRows
) {
    New-FmmCliCleanDir -Path $CaseRoot -WorkRoot $WorkRoot
    $mode = if ($Backend -eq 'fmm') { 'treecode' } else { 'direct' }
    $fmmMode = if ($Backend -eq 'fmm') { 'single-matvec-analytic' } else { 'off' }
    $arg = Write-FmmCliArgFile -RunDir $CaseRoot -RepoRoot $RepoRoot `
        -PqrPath $PqrPath -Backend tabi -Mode $mode -FmmMode $fmmMode `
        -SurfacePoints $SurfacePoints -MaxIter $MaxIter -Restart $Restart `
        -Tolerance $Tolerance -Threads $Threads -SaltM 0.15 `
        -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK $TemperatureK -OpeningTheta $OpeningTheta `
        -AuditSampleCount $AuditRows
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $CaseRoot `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        return [pscustomobject]@{
            run=$run
            result=(New-ScopeTEmptyCase "exit_$($run.exit_code)")
        }
    }
    try {
        $result = Read-FmmCliCaseResult -RunDir $CaseRoot -ExpectedFmmMode $fmmMode `
            -RequireAudit ($Backend -eq 'fmm')
    } catch {
        $result = New-ScopeTEmptyCase "parse_failed:$($_.Exception.Message)"
    }
    return [pscustomobject]@{ run=$run; result=$result }
}

function Assert-ScopeTGeneratedArgEquals([string]$Root) {
    $bad = @(Get-ChildItem -LiteralPath $Root -Recurse -Filter arg*.in |
        ForEach-Object {
            $path = $_.FullName
            Get-Content -LiteralPath $path | Where-Object {
                -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
            } | ForEach-Object { "$path :: $_" }
        })
    if ($bad.Count -ne 0) {
        throw "Generated arg.in contains non equals-syntax options: $($bad -join '; ')"
    }
}

function Get-ScopeTEnergyGate([object]$Direct) {
    $magnitude = [Math]::Max(1.0, [Math]::Abs($Direct.energy_kj_mol))
    return [Math]::Max(5.0, 0.05 * $magnitude)
}

function Get-ScopeTSpeedup([object]$Direct, [object]$Fmm) {
    if (-not (Test-ScopeTFinite $Direct.matvec_seconds_avg) -or
        -not (Test-ScopeTFinite $Fmm.matvec_seconds_avg) -or
        $Fmm.matvec_seconds_avg -le 0.0) {
        return [double]::NaN
    }
    return $Direct.matvec_seconds_avg / $Fmm.matvec_seconds_avg
}

function Test-ScopeTFmmBasic([object]$Fmm, [int]$FmmExit) {
    if ($FmmExit -ne 0) { return $false }
    if ($Fmm.fmm_m2l_mode -ne 'single-matvec-analytic' -or
        $Fmm.fmm_far_pairs -le 0 -or $Fmm.fmm_finite_outputs -ne 'true') {
        return $false
    }
    foreach ($value in @(
        $Fmm.energy_kj_mol, $Fmm.true_residual,
        $Fmm.matvec_seconds_avg, $Fmm.tree_vs_direct_matvec_rel)) {
        if (-not (Test-ScopeTFinite $value)) { return $false }
    }
    if ($Fmm.matvec_audit_rows -le 0 -or
        $Fmm.tree_vs_direct_matvec_rel -gt $AuditGate) {
        return $false
    }
    return $true
}

function Get-ScopeTFmmOnlyClass(
    [object]$FmmCase,
    [double]$Tolerance
) {
    if ($FmmCase.run.exit_code -eq 124) { return 'timeout_not_verified' }
    if ($FmmCase.run.exit_code -ne 0) { return 'status_fail' }
    $fmm = $FmmCase.result
    if ($fmm.fmm_m2l_mode -ne 'single-matvec-analytic' -or
        $fmm.fmm_far_pairs -le 0 -or $fmm.fmm_finite_outputs -ne 'true') {
        return 'telemetry_fail'
    }
    foreach ($value in @(
        $fmm.energy_kj_mol, $fmm.true_residual,
        $fmm.matvec_seconds_avg, $fmm.tree_vs_direct_matvec_rel)) {
        if (-not (Test-ScopeTFinite $value)) { return 'nonfinite_fail' }
    }
    if ($fmm.matvec_audit_rows -le 0 -or
        $fmm.tree_vs_direct_matvec_rel -gt $AuditGate) {
        return 'accuracy_limited'
    }
    if ($fmm.status -eq 'ok') {
        if ($Tolerance -le 1.0000001e-6) {
            return 'strict_fmm_converged_pending_direct'
        }
        return 'practical_fmm_converged_pending_direct'
    }
    return 'convergence_pending'
}

function Get-ScopeTVerifiedClass(
    [object]$DirectCase,
    [object]$FmmCase,
    [double]$Tolerance
) {
    if ($DirectCase.run.exit_code -eq 124 -or $FmmCase.run.exit_code -eq 124) {
        return 'timeout_not_verified'
    }
    if ($DirectCase.run.exit_code -ne 0 -or $FmmCase.run.exit_code -ne 0) {
        return 'status_fail'
    }
    $direct = $DirectCase.result
    $fmm = $FmmCase.result
    if (-not (Test-ScopeTFmmBasic -Fmm $fmm -FmmExit $FmmCase.run.exit_code)) {
        return (Get-ScopeTFmmOnlyClass -FmmCase $FmmCase -Tolerance $Tolerance)
    }
    foreach ($value in @(
        $direct.energy_kj_mol, $fmm.energy_kj_mol,
        $direct.true_residual, $fmm.true_residual,
        $direct.matvec_seconds_avg, $fmm.matvec_seconds_avg)) {
        if (-not (Test-ScopeTFinite $value)) { return 'nonfinite_fail' }
    }
    if ([Math]::Abs($fmm.energy_kj_mol - $direct.energy_kj_mol) -gt
        (Get-ScopeTEnergyGate -Direct $direct)) {
        return 'accuracy_limited'
    }
    $speedup = Get-ScopeTSpeedup -Direct $direct -Fmm $fmm
    if (-not (Test-ScopeTFinite $speedup) -or $speedup -lt 1.2) {
        return 'speed_limited'
    }
    if ($direct.status -eq 'ok' -and $fmm.status -eq 'ok') {
        if ($Tolerance -le 1.0000001e-6) {
            return 'strict_converged_verified'
        }
        return 'practical_converged'
    }
    return 'convergence_pending'
}

function Get-ScopeTBestCandidateCondition([object[]]$FmmRuns) {
    $basic = @($FmmRuns | Where-Object {
        Test-ScopeTFmmBasic -Fmm $_.case.result -FmmExit $_.case.run.exit_code
    })
    $strictOk = @($basic | Where-Object {
        $_.condition.tolerance -le 1.0000001e-6 -and $_.case.result.status -eq 'ok'
    } | Sort-Object { [double]$_.case.result.true_residual })
    if ($strictOk.Count -gt 0) { return $strictOk[0].condition.name }

    $practicalOk = @($basic | Where-Object {
        $_.condition.tolerance -gt 1.0000001e-6 -and $_.case.result.status -eq 'ok'
    } | Sort-Object { [double]$_.case.result.true_residual })
    if ($practicalOk.Count -gt 0) { return $practicalOk[0].condition.name }

    $baseline = @($basic | Where-Object {
        $_.condition.name -eq 'baseline_240_r64_1e-6'
    } | Select-Object -First 1)
    if ($baseline.Count -eq 0) { return '' }
    $baselineResidual = [double]$baseline[0].case.result.true_residual
    $improved = @($basic | Where-Object {
        [double]$_.case.result.true_residual -le 0.25 * $baselineResidual
    } | Sort-Object { [double]$_.case.result.true_residual })
    if ($improved.Count -gt 0) { return $improved[0].condition.name }
    return ''
}

function Get-ScopeTFinalClass([object[]]$Rows) {
    if (@($Rows | Where-Object { $_.scope_t_class -eq 'strict_converged_verified' }).Count -gt 0) {
        return 'strict_converged_verified'
    }
    if (@($Rows | Where-Object { $_.scope_t_class -eq 'strict_fmm_converged_pending_direct' }).Count -gt 0) {
        return 'strict_fmm_converged_pending_direct'
    }
    if (@($Rows | Where-Object { $_.scope_t_class -eq 'practical_converged' }).Count -gt 0) {
        return 'practical_converged'
    }
    if (@($Rows | Where-Object { $_.scope_t_class -eq 'accuracy_limited' }).Count -gt 0) {
        return 'accuracy_limited'
    }
    $fmmRows = @($Rows | Where-Object {
        $_.backend -eq 'fmm' -and (Test-ScopeTFinite ([double]$_.true_residual))
    })
    if ($fmmRows.Count -eq 0) {
        if (@($Rows | Where-Object { $_.scope_t_class -eq 'timeout_not_verified' }).Count -gt 0) {
            return 'timeout_not_verified'
        }
        return 'preconditioner_needed'
    }
    $baseline = @($fmmRows | Where-Object {
        $_.condition -eq 'baseline_240_r64_1e-6'
    } | Select-Object -First 1)
    if ($baseline.Count -eq 0) { return 'preconditioner_needed' }
    $baselineResidual = [double]$baseline[0].true_residual
    $bestResidual = ($fmmRows | ForEach-Object {
        [double]$_.true_residual
    } | Measure-Object -Minimum).Minimum
    $iter480 = @($fmmRows | Where-Object {
        $_.condition -eq 'iter_480_r64_1e-6'
    } | Select-Object -First 1)
    $restartRows = @($fmmRows | Where-Object {
        $_.condition -match '^restart_480_r(96|128)_1e-6$'
    })
    if ($iter480.Count -gt 0 -and $restartRows.Count -gt 0) {
        $bestRestartResidual = ($restartRows | ForEach-Object {
            [double]$_.true_residual
        } | Measure-Object -Minimum).Minimum
        if ($bestRestartResidual -le 0.5 * [double]$iter480[0].true_residual) {
            return 'restart_limited'
        }
    }
    if ($bestResidual -le 0.25 * $baselineResidual) {
        return 'iteration_limited'
    }
    if ($bestResidual -gt 0.5 * $baselineResidual) {
        return 'preconditioner_needed'
    }
    if (@($Rows | Where-Object { $_.scope_t_class -eq 'timeout_not_verified' }).Count -gt 0) {
        return 'timeout_not_verified'
    }
    return 'iteration_limited'
}

function New-ScopeTRow(
    [string]$ConditionName,
    [object]$Condition,
    [object]$Sample,
    [ValidateSet('direct','fmm')][string]$Backend,
    [object]$Case,
    [string]$DirectVerificationStatus,
    [AllowNull()][object]$DirectCase,
    [string]$Class
) {
    $directResult = if ($null -ne $DirectCase) { $DirectCase.result } else { $null }
    $directRun = if ($null -ne $DirectCase) { $DirectCase.run } else { $null }
    $energyDelta = if ($null -ne $directResult -and
        (Test-ScopeTFinite $directResult.energy_kj_mol) -and
        (Test-ScopeTFinite $Case.result.energy_kj_mol)) {
        [Math]::Abs($Case.result.energy_kj_mol - $directResult.energy_kj_mol)
    } else { [double]::NaN }
    $speedup = if ($null -ne $directResult) {
        Get-ScopeTSpeedup -Direct $directResult -Fmm $Case.result
    } else { [double]::NaN }
    return [pscustomobject]@{
        condition=$ConditionName; fixture='barnase'; atom_cap=$AtomCap
        backend=$Backend; input_atoms=$Sample.input_atoms
        sampled_atoms=$Sample.sampled_atoms; surface_points=$SurfacePoints
        opening_theta=$OpeningTheta; audit_sample_count=$AuditSampleCount
        max_iter=$Condition.max_iter; restart=$Condition.restart
        tolerance=$Condition.tolerance; exit_code=$Case.run.exit_code
        status=$Case.result.status; iterations=$Case.result.iterations
        final_residual=$Case.result.final_residual
        true_residual=$Case.result.true_residual
        energy_kj_mol=$Case.result.energy_kj_mol
        tree_vs_direct_matvec_rel=$Case.result.tree_vs_direct_matvec_rel
        matvec_audit_rows=$Case.result.matvec_audit_rows
        matvec_seconds_avg=$Case.result.matvec_seconds_avg
        wall_seconds=$Case.run.wall_seconds
        peak_working_set_bytes=$Case.run.peak_working_set_bytes
        peak_working_set_mb=$Case.run.peak_working_set_mb
        solver_memory_bytes=$Case.result.memory_bytes
        fmm_m2l_mode=$Case.result.fmm_m2l_mode
        fmm_far_pairs=$Case.result.fmm_far_pairs
        fmm_setup_seconds=$Case.result.fmm_setup_seconds
        fmm_apply_seconds=$Case.result.fmm_apply_seconds
        direct_verification_status=$DirectVerificationStatus
        direct_exit_code=if ($null -ne $directRun) { $directRun.exit_code } else { 0 }
        direct_status=if ($null -ne $directResult) { $directResult.status } else { '' }
        direct_true_residual=if ($null -ne $directResult) { $directResult.true_residual } else { [double]::NaN }
        direct_energy_kj_mol=if ($null -ne $directResult) { $directResult.energy_kj_mol } else { [double]::NaN }
        energy_delta_abs_kj_mol=$energyDelta
        matvec_speedup=$speedup
        scope_t_class=$Class
    }
}

$SamplesRoot = Join-Path $WorkRoot 'samples'
New-FmmCliCleanDir -Path $SamplesRoot -WorkRoot $WorkRoot
$sourcePqr = Resolve-FmmCliFixture -Fixture barnase -RepoRoot $RepoRoot `
    -RunDir $SamplesRoot

$Rows = New-Object System.Collections.Generic.List[object]
$fmmRuns = New-Object System.Collections.Generic.List[object]
$directRuns = @{}
$conditions = @(Get-ScopeTConditions)

if ($null -eq $sourcePqr -or -not (Test-Path -LiteralPath $sourcePqr)) {
    $emptyRun = [pscustomobject]@{
        exit_code=0; timed_out=$false; wall_seconds=[double]::NaN
        peak_working_set_bytes=0; peak_working_set_mb=0.0
    }
    $emptyCase = [pscustomobject]@{
        run=$emptyRun
        result=(New-ScopeTEmptyCase 'skipped_missing_fixture')
    }
    $emptySample = [pscustomobject]@{ input_atoms=0; sampled_atoms=0 }
    $emptyCondition = [pscustomobject]@{
        name='fixture_missing'; max_iter=0; restart=0; tolerance=0.0
    }
    $Rows.Add((New-ScopeTRow -ConditionName 'fixture_missing' `
        -Condition $emptyCondition -Sample $emptySample -Backend fmm `
        -Case $emptyCase -DirectVerificationStatus 'skipped_missing_fixture' `
        -DirectCase $null -Class 'skipped_missing_fixture'))
} else {
    $samplePqr = Join-Path $SamplesRoot "barnase_cap$AtomCap.pqr"
    $sample = Write-FmmCliSampledPqr -InputPath $sourcePqr -OutputPath $samplePqr `
        -Cap $AtomCap

    foreach ($condition in $conditions) {
        $conditionSafe = Get-FmmCliSafeName $condition.name
        $conditionRoot = Join-Path $WorkRoot $conditionSafe
        New-FmmCliCleanDir -Path $conditionRoot -WorkRoot $WorkRoot
        Write-Host "Scope T Barnase cap-256 FMM: $($condition.name)"
        $fmmCase = Invoke-ScopeTCase -CaseRoot (Join-Path $conditionRoot 'fmm') `
            -PqrPath $samplePqr -Backend fmm -MaxIter $condition.max_iter `
            -Restart $condition.restart -Tolerance $condition.tolerance `
            -AuditRows $AuditSampleCount
        $fmmRuns.Add([pscustomobject]@{
            condition=$condition
            condition_root=$conditionRoot
            case=$fmmCase
        })
    }

    $candidateCondition = if ($Quick) {
        ''
    } else {
        Get-ScopeTBestCandidateCondition -FmmRuns @($fmmRuns.ToArray())
    }
    $directConditionNames = New-Object System.Collections.Generic.HashSet[string]
    [void]$directConditionNames.Add('baseline_240_r64_1e-6')
    if (-not [string]::IsNullOrWhiteSpace($candidateCondition)) {
        [void]$directConditionNames.Add($candidateCondition)
    }

    if (-not $Quick) {
        foreach ($condition in $conditions) {
            if (-not $directConditionNames.Contains($condition.name)) { continue }
            $conditionSafe = Get-FmmCliSafeName $condition.name
            $conditionRoot = Join-Path $WorkRoot $conditionSafe
            Write-Host "Scope T Barnase cap-256 direct verification: $($condition.name)"
            $directCase = Invoke-ScopeTCase -CaseRoot (Join-Path $conditionRoot 'direct') `
                -PqrPath $samplePqr -Backend direct -MaxIter $condition.max_iter `
                -Restart $condition.restart -Tolerance $condition.tolerance `
                -AuditRows 0
            $directRuns[$condition.name] = $directCase
            $Rows.Add((New-ScopeTRow -ConditionName $condition.name `
                -Condition $condition -Sample $sample -Backend direct `
                -Case $directCase -DirectVerificationStatus 'not_applicable' `
                -DirectCase $null -Class 'direct_reference'))
        }
    }

    foreach ($entry in @($fmmRuns.ToArray())) {
        $condition = $entry.condition
        $directCase = if ($directRuns.ContainsKey($condition.name)) {
            $directRuns[$condition.name]
        } else {
            $null
        }
        $directStatus = if ($Quick) {
            'quick_fmm_only'
        } elseif ($null -ne $directCase) {
            'verified'
        } else {
            'skipped_not_selected'
        }
        $rowClass = if ($null -ne $directCase) {
            Get-ScopeTVerifiedClass -DirectCase $directCase `
                -FmmCase $entry.case -Tolerance $condition.tolerance
        } else {
            Get-ScopeTFmmOnlyClass -FmmCase $entry.case `
                -Tolerance $condition.tolerance
        }
        $Rows.Add((New-ScopeTRow -ConditionName $condition.name `
            -Condition $condition -Sample $sample -Backend fmm `
            -Case $entry.case -DirectVerificationStatus $directStatus `
            -DirectCase $directCase -Class $rowClass))
    }
}

$Rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation -Encoding ASCII
Assert-ScopeTGeneratedArgEquals -Root $WorkRoot

$rowsArray = @($Rows.ToArray())
$finalClass = Get-ScopeTFinalClass -Rows $rowsArray
$fmmRows = @($rowsArray | Where-Object { $_.backend -eq 'fmm' })
$directRows = @($rowsArray | Where-Object { $_.backend -eq 'direct' })
$timeoutRows = @($rowsArray | Where-Object { $_.scope_t_class -eq 'timeout_not_verified' })
$bestFinite = @($fmmRows | Where-Object {
    Test-ScopeTFinite ([double]$_.true_residual)
} | Sort-Object { [double]$_.true_residual } | Select-Object -First 1)

$nextScope = switch ($finalClass) {
    'strict_converged_verified' { 'plan_broader_cap256_hidden_matrix' }
    'strict_fmm_converged_pending_direct' { 'run_direct_verification_or_longer_timeout' }
    'practical_converged' { 'plan_tolerance_policy_discussion_keep_hidden' }
    'iteration_limited' { 'plan_cap256_gmres_iteration_budget_gate' }
    'restart_limited' { 'plan_cap256_gmres_restart_preset_gate' }
    'preconditioner_needed' { 'plan_cap256_preconditioner_work' }
    'accuracy_limited' { 'pause_promotion_fix_accuracy_or_energy_gate' }
    'timeout_not_verified' { 'run_longer_or_reduce_verifier_cost' }
    default { 'hold_cap256_promotion' }
}

$decision = [pscustomobject]@{
    scope_t_class=$finalClass
    quick_mode=$Quick.IsPresent
    rows=$rowsArray.Count
    fmm_rows=$fmmRows.Count
    direct_rows=$directRows.Count
    timeout_rows=$timeoutRows.Count
    best_condition=if ($bestFinite.Count -gt 0) { $bestFinite[0].condition } else { '' }
    best_fmm_true_residual=if ($bestFinite.Count -gt 0) { $bestFinite[0].true_residual } else { [double]::NaN }
    production_default='grid'
    fmm_status='hidden_experimental'
    next_scope=$nextScope
}
$decision | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation -Encoding ASCII

$Summary = New-Object System.Collections.Generic.List[string]
$Summary.Add('TITAN FMM/M2L Scope T Barnase Cap-256 Convergence Sweep')
$Summary.Add('=========================================================')
$Summary.Add('')
$Summary.Add("quick_mode: $($Quick.IsPresent)")
$Summary.Add('scope: hidden_cli_bounded_barnase_cap256_convergence')
$Summary.Add('hidden_option: -pbe-tabi-fmm-m2l=single-matvec-analytic')
$Summary.Add('production_default: grid')
$Summary.Add('fmm_status: hidden_experimental')
$Summary.Add("atom_cap: $AtomCap")
$Summary.Add("surface_points: $SurfacePoints")
$Summary.Add("opening_theta: $OpeningTheta")
$Summary.Add("audit_sample_count: $AuditSampleCount")
$Summary.Add("case_timeout_seconds: $CaseTimeoutSeconds")
$Summary.Add("scope_t_class: $finalClass")
$Summary.Add("next_scope: $nextScope")
$Summary.Add('')
$Summary.Add('Rows:')
foreach ($row in $rowsArray) {
    $Summary.Add("row: condition=$($row.condition) backend=$($row.backend) max_iter=$($row.max_iter) restart=$($row.restart) tolerance=$($row.tolerance) status=$($row.status) iter=$($row.iterations) true_residual=$($row.true_residual) energy=$($row.energy_kj_mol) audit_rel=$($row.tree_vs_direct_matvec_rel) direct_verification=$($row.direct_verification_status) direct_status=$($row.direct_status) energy_delta_abs=$($row.energy_delta_abs_kj_mol) speedup=$($row.matvec_speedup) wall=$($row.wall_seconds) peak_mb=$($row.peak_working_set_mb) class=$($row.scope_t_class)")
}
$Summary.Add('')
$Summary.Add('classification_note: strict_converged_verified requires both direct and FMM status ok at 1e-6; practical_converged applies only to 1e-4 rows that pass verification gates.')
$Summary.Add('direct_policy: full mode verifies direct only for the baseline and the best/promising FMM condition; quick mode is FMM-only.')
$Summary.Add('timeout_policy: timeout_not_verified blocks promotion planning but is not classified as solver failure.')
$Summary.Add('barnase_note: Scope T diagnoses the bounded Barnase cap-256 convergence/status gap from refreshed Scope R; it does not retry legacy TABI leaf-kernel work.')
$Summary.Add("matrix_csv: $MatrixCsv")
$Summary.Add("decision_csv: $DecisionCsv")
Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

Write-Host "FMM/M2L Scope T Barnase cap-256 convergence sweep complete: $DecisionCsv"
