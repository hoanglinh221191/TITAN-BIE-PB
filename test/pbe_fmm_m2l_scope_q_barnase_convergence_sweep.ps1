param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 4,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_q')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fmm_m2l_scope_q_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_q_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_q_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Q_BARNASE_CONVERGENCE_REPORT.txt'

$AtomCap = 128
$SurfacePoints = 8
$OpeningTheta = 0.6
$AuditSampleCount = if ($Quick) { 128 } else { 512 }
$TemperatureK = 298.15

if ($Quick) {
    $CaseTimeoutSeconds = [Math]::Min($CaseTimeoutSeconds, 420)
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope Q work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function New-ScopeQEmptyCase([string]$Status) {
    return [pscustomobject]@{
        status=$Status; panels=0; iterations=0; final_residual=[double]::NaN
        true_residual=[double]::NaN; energy_kj_mol=[double]::NaN
        matvec_seconds_avg=[double]::NaN; total_seconds=[double]::NaN
        fmm_m2l_mode=''; fmm_near_pairs=0; fmm_far_pairs=0
        fmm_setup_seconds=0.0; fmm_apply_seconds=0.0
        fmm_finite_outputs='false'; tree_vs_direct_matvec_rel=[double]::NaN
        matvec_audit_rows=0; matvec_audit_seconds=0.0
        reaction_energy_source=''; timing_backend=''
    }
}

function Test-ScopeQFinite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function Get-ScopeQEnergyGate([object]$Direct) {
    $magnitude = [Math]::Max(1.0, [Math]::Abs($Direct.energy_kj_mol))
    return [Math]::Max(2.0, 0.05 * $magnitude)
}

function Get-ScopeQSpeedup([object]$Direct, [object]$Fmm) {
    if (-not (Test-ScopeQFinite $Direct.matvec_seconds_avg) -or
        -not (Test-ScopeQFinite $Fmm.matvec_seconds_avg) -or
        $Fmm.matvec_seconds_avg -le 0.0) {
        return [double]::NaN
    }
    return $Direct.matvec_seconds_avg / $Fmm.matvec_seconds_avg
}

function Get-ScopeQRowClass(
    [int]$DirectExit,
    [int]$FmmExit,
    [object]$Direct,
    [object]$Fmm,
    [double]$Tolerance
) {
    if ($DirectExit -eq 124 -or $FmmExit -eq 124) {
        return 'timeout_not_verified'
    }
    if ($DirectExit -ne 0 -or $FmmExit -ne 0) {
        return 'status_fail'
    }
    if ($Fmm.fmm_m2l_mode -ne 'single-matvec-analytic' -or
        $Fmm.fmm_far_pairs -le 0 -or $Fmm.fmm_finite_outputs -ne 'true') {
        return 'telemetry_fail'
    }
    foreach ($value in @(
        $Direct.energy_kj_mol, $Fmm.energy_kj_mol,
        $Direct.true_residual, $Fmm.true_residual,
        $Direct.matvec_seconds_avg, $Fmm.matvec_seconds_avg)) {
        if (-not (Test-ScopeQFinite $value)) {
            return 'nonfinite_fail'
        }
    }
    if ($Fmm.matvec_audit_rows -le 0 -or
        -not (Test-ScopeQFinite $Fmm.tree_vs_direct_matvec_rel) -or
        $Fmm.tree_vs_direct_matvec_rel -gt 1.0e-2) {
        return 'accuracy_limited'
    }
    if ([Math]::Abs($Fmm.energy_kj_mol - $Direct.energy_kj_mol) -gt
        (Get-ScopeQEnergyGate -Direct $Direct)) {
        return 'accuracy_limited'
    }
    $speedup = Get-ScopeQSpeedup -Direct $Direct -Fmm $Fmm
    if (-not (Test-ScopeQFinite $speedup) -or $speedup -lt 1.2) {
        return 'speed_limited'
    }
    if ($Direct.status -eq 'ok' -and $Fmm.status -eq 'ok') {
        if ($Tolerance -le 1.0000001e-6) {
            return 'strict_converged'
        }
        return 'practical_converged'
    }
    return 'convergence_pending'
}

function Assert-ScopeQGeneratedArgEquals([string]$Root) {
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

function Invoke-ScopeQCase(
    [string]$CaseRoot,
    [string]$PqrPath,
    [string]$Mode,
    [string]$FmmMode,
    [int]$MaxIter,
    [int]$Restart,
    [double]$Tolerance,
    [int]$AuditRows
) {
    New-FmmCliCleanDir -Path $CaseRoot -WorkRoot $WorkRoot
    $arg = Write-FmmCliArgFile -RunDir $CaseRoot -RepoRoot $RepoRoot `
        -PqrPath $PqrPath -Mode $Mode -FmmMode $FmmMode `
        -SurfacePoints $SurfacePoints -MaxIter $MaxIter -Restart $Restart `
        -Tolerance $Tolerance -Threads $Threads -SaltM 0.15 `
        -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK $TemperatureK -OpeningTheta $OpeningTheta `
        -AuditSampleCount $AuditRows
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $CaseRoot `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $result = if ($run.exit_code -eq 0) {
        Read-FmmCliCaseResult -RunDir $CaseRoot -ExpectedFmmMode $FmmMode `
            -RequireAudit ($FmmMode -eq 'single-matvec-analytic')
    } else {
        New-ScopeQEmptyCase "exit_$($run.exit_code)"
    }
    return [pscustomobject]@{ run=$run; result=$result }
}

function Get-ScopeQConditions {
    if ($Quick) {
        return @(
            [pscustomobject]@{ name='baseline_120_r24_1e-6'; max_iter=120; restart=24; tolerance=1.0e-6 },
            [pscustomobject]@{ name='practical_240_r48_1e-4'; max_iter=240; restart=48; tolerance=1.0e-4 }
        )
    }
    return @(
        [pscustomobject]@{ name='baseline_120_r24_1e-6'; max_iter=120; restart=24; tolerance=1.0e-6 },
        [pscustomobject]@{ name='iter_240_r24_1e-6'; max_iter=240; restart=24; tolerance=1.0e-6 },
        [pscustomobject]@{ name='restart_240_r48_1e-6'; max_iter=240; restart=48; tolerance=1.0e-6 },
        [pscustomobject]@{ name='restart_240_r64_1e-6'; max_iter=240; restart=64; tolerance=1.0e-6 },
        [pscustomobject]@{ name='practical_240_r48_1e-4'; max_iter=240; restart=48; tolerance=1.0e-4 },
        [pscustomobject]@{ name='long_360_r64_1e-6'; max_iter=360; restart=64; tolerance=1.0e-6 }
    )
}

function Get-ScopeQFinalClass([object[]]$Rows) {
    if (@($Rows | Where-Object { $_.scope_q_class -eq 'strict_converged' }).Count -gt 0) {
        return 'strict_converged'
    }
    if (@($Rows | Where-Object { $_.scope_q_class -eq 'practical_converged' }).Count -gt 0) {
        return 'practical_converged'
    }
    if (@($Rows | Where-Object { $_.scope_q_class -eq 'accuracy_limited' }).Count -gt 0) {
        return 'accuracy_limited'
    }
    $finiteRows = @($Rows | Where-Object {
        Test-ScopeQFinite ([double]$_.fmm_true_residual)
    })
    if ($finiteRows.Count -eq 0) {
        if (@($Rows | Where-Object { $_.scope_q_class -eq 'timeout_not_verified' }).Count -gt 0) {
            return 'timeout_not_verified'
        }
        return 'not_converging'
    }
    $baseline = @($finiteRows | Where-Object { $_.condition -eq 'baseline_120_r24_1e-6' } | Select-Object -First 1)
    $iter = @($finiteRows | Where-Object { $_.condition -eq 'iter_240_r24_1e-6' } | Select-Object -First 1)
    if ($baseline.Count -gt 0 -and $iter.Count -gt 0 -and
        [double]$iter[0].fmm_true_residual -lt 0.5 * [double]$baseline[0].fmm_true_residual) {
        return 'iteration_limited'
    }
    $restartRows = @($finiteRows | Where-Object {
        $_.condition -match '^restart_240_r(48|64)_1e-6$'
    })
    if ($baseline.Count -gt 0 -and $restartRows.Count -gt 0) {
        $bestRestart = ($restartRows | ForEach-Object {
            [double]$_.fmm_true_residual
        } | Measure-Object -Minimum).Minimum
        if ($bestRestart -lt 0.5 * [double]$baseline[0].fmm_true_residual) {
            return 'restart_limited'
        }
    }
    if (@($Rows | Where-Object { $_.scope_q_class -eq 'timeout_not_verified' }).Count -gt 0) {
        return 'timeout_not_verified'
    }
    return 'not_converging'
}

$SamplesRoot = Join-Path $WorkRoot 'samples'
New-FmmCliCleanDir -Path $SamplesRoot -WorkRoot $WorkRoot
$sourcePqr = Resolve-FmmCliFixture -Fixture barnase -RepoRoot $RepoRoot `
    -RunDir $SamplesRoot

$Rows = New-Object System.Collections.Generic.List[object]

if ($null -eq $sourcePqr -or -not (Test-Path -LiteralPath $sourcePqr)) {
    $Rows.Add([pscustomobject]@{
        condition='fixture_missing'; fixture='barnase'; atom_cap=$AtomCap
        input_atoms=0; sampled_atoms=0; surface_points=$SurfacePoints
        opening_theta=$OpeningTheta; audit_sample_count=$AuditSampleCount
        max_iter=0; restart=0; tolerance=0.0
        direct_exit_code=0; fmm_exit_code=0
        direct_status='skipped_missing_fixture'; fmm_status='skipped_missing_fixture'
        direct_iterations=0; fmm_iterations=0
        direct_final_residual=[double]::NaN; fmm_final_residual=[double]::NaN
        direct_true_residual=[double]::NaN; fmm_true_residual=[double]::NaN
        direct_energy_kj_mol=[double]::NaN; fmm_energy_kj_mol=[double]::NaN
        energy_delta_abs_kj_mol=[double]::NaN
        tree_vs_direct_matvec_rel=[double]::NaN; matvec_audit_rows=0
        direct_matvec_seconds_avg=[double]::NaN
        fmm_matvec_seconds_avg=[double]::NaN; matvec_speedup=[double]::NaN
        fmm_m2l_mode=''; fmm_far_pairs=0; fmm_setup_seconds=0.0
        fmm_apply_seconds=0.0; direct_wall_seconds=[double]::NaN
        fmm_wall_seconds=[double]::NaN; scope_q_class='skipped_missing_fixture'
    })
} else {
    $samplePqr = Join-Path $SamplesRoot "barnase_cap$AtomCap.pqr"
    $sample = Write-FmmCliSampledPqr -InputPath $sourcePqr -OutputPath $samplePqr `
        -Cap $AtomCap
    foreach ($condition in Get-ScopeQConditions) {
        $conditionSafe = Get-FmmCliSafeName $condition.name
        $conditionRoot = Join-Path $WorkRoot $conditionSafe
        New-FmmCliCleanDir -Path $conditionRoot -WorkRoot $WorkRoot
        Write-Host "Scope Q Barnase direct: $($condition.name)"
        $directCase = Invoke-ScopeQCase -CaseRoot (Join-Path $conditionRoot 'direct') `
            -PqrPath $samplePqr -Mode direct -FmmMode off `
            -MaxIter $condition.max_iter -Restart $condition.restart `
            -Tolerance $condition.tolerance -AuditRows 0
        Write-Host "Scope Q Barnase FMM: $($condition.name)"
        $fmmCase = Invoke-ScopeQCase -CaseRoot (Join-Path $conditionRoot 'fmm') `
            -PqrPath $samplePqr -Mode treecode `
            -FmmMode single-matvec-analytic -MaxIter $condition.max_iter `
            -Restart $condition.restart -Tolerance $condition.tolerance `
            -AuditRows $AuditSampleCount
        $rowClass = Get-ScopeQRowClass -DirectExit $directCase.run.exit_code `
            -FmmExit $fmmCase.run.exit_code -Direct $directCase.result `
            -Fmm $fmmCase.result -Tolerance $condition.tolerance
        $delta = $fmmCase.result.energy_kj_mol - $directCase.result.energy_kj_mol
        $Rows.Add([pscustomobject]@{
            condition=$condition.name; fixture='barnase'; atom_cap=$AtomCap
            input_atoms=$sample.input_atoms; sampled_atoms=$sample.sampled_atoms
            surface_points=$SurfacePoints; opening_theta=$OpeningTheta
            audit_sample_count=$AuditSampleCount; max_iter=$condition.max_iter
            restart=$condition.restart; tolerance=$condition.tolerance
            direct_exit_code=$directCase.run.exit_code
            fmm_exit_code=$fmmCase.run.exit_code
            direct_status=$directCase.result.status
            fmm_status=$fmmCase.result.status
            direct_iterations=$directCase.result.iterations
            fmm_iterations=$fmmCase.result.iterations
            direct_final_residual=$directCase.result.final_residual
            fmm_final_residual=$fmmCase.result.final_residual
            direct_true_residual=$directCase.result.true_residual
            fmm_true_residual=$fmmCase.result.true_residual
            direct_energy_kj_mol=$directCase.result.energy_kj_mol
            fmm_energy_kj_mol=$fmmCase.result.energy_kj_mol
            energy_delta_abs_kj_mol=[Math]::Abs($delta)
            tree_vs_direct_matvec_rel=$fmmCase.result.tree_vs_direct_matvec_rel
            matvec_audit_rows=$fmmCase.result.matvec_audit_rows
            direct_matvec_seconds_avg=$directCase.result.matvec_seconds_avg
            fmm_matvec_seconds_avg=$fmmCase.result.matvec_seconds_avg
            matvec_speedup=(Get-ScopeQSpeedup -Direct $directCase.result -Fmm $fmmCase.result)
            fmm_m2l_mode=$fmmCase.result.fmm_m2l_mode
            fmm_far_pairs=$fmmCase.result.fmm_far_pairs
            fmm_setup_seconds=$fmmCase.result.fmm_setup_seconds
            fmm_apply_seconds=$fmmCase.result.fmm_apply_seconds
            direct_wall_seconds=$directCase.run.wall_seconds
            fmm_wall_seconds=$fmmCase.run.wall_seconds
            scope_q_class=$rowClass
        })
    }
}

$Rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation -Encoding ASCII
Assert-ScopeQGeneratedArgEquals -Root $WorkRoot

$rowsArray = @($Rows.ToArray())
$finalClass = Get-ScopeQFinalClass $rowsArray
$strictRows = @($rowsArray | Where-Object { $_.scope_q_class -eq 'strict_converged' })
$practicalRows = @($rowsArray | Where-Object { $_.scope_q_class -eq 'practical_converged' })
$timeoutRows = @($rowsArray | Where-Object { $_.scope_q_class -eq 'timeout_not_verified' })
$bestFinite = @($rowsArray | Where-Object {
    Test-ScopeQFinite ([double]$_.fmm_true_residual)
} | Sort-Object {[double]$_.fmm_true_residual} | Select-Object -First 1)

$nextScope = switch ($finalClass) {
    'strict_converged' { 'consider_scope_r_hidden_barnase_cap256_or_broader_regression' }
    'practical_converged' { 'plan_tolerance_policy_discussion_keep_hidden' }
    'iteration_limited' { 'plan_gmres_iteration_budget_gate' }
    'restart_limited' { 'plan_gmres_restart_tuning_gate' }
    'accuracy_limited' { 'pause_promotion_fix_accuracy_or_energy_gate' }
    'timeout_not_verified' { 'run_longer_or_reduce_verifier_cost' }
    default { 'hold_barnase_promotion_scope_n_o_only' }
}

$decision = [pscustomobject]@{
    scope_q_class=$finalClass
    quick_mode=$Quick.IsPresent
    rows=$rowsArray.Count
    strict_rows=$strictRows.Count
    practical_rows=$practicalRows.Count
    timeout_rows=$timeoutRows.Count
    best_condition=if ($bestFinite.Count -gt 0) { $bestFinite[0].condition } else { '' }
    best_fmm_true_residual=if ($bestFinite.Count -gt 0) { $bestFinite[0].fmm_true_residual } else { [double]::NaN }
    production_default='grid'
    fmm_status='hidden_experimental'
    next_scope=$nextScope
}
$decision | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation -Encoding ASCII

$Summary = New-Object System.Collections.Generic.List[string]
$Summary.Add('TITAN FMM/M2L Scope Q Barnase Convergence Sweep')
$Summary.Add('=================================================')
$Summary.Add('')
$Summary.Add("quick_mode: $($Quick.IsPresent)")
$Summary.Add('scope: hidden_cli_bounded_barnase_convergence')
$Summary.Add('hidden_option: -pbe-tabi-fmm-m2l=single-matvec-analytic')
$Summary.Add('production_default: grid')
$Summary.Add('fmm_status: hidden_experimental')
$Summary.Add("atom_cap: $AtomCap")
$Summary.Add("surface_points: $SurfacePoints")
$Summary.Add("opening_theta: $OpeningTheta")
$Summary.Add("audit_sample_count: $AuditSampleCount")
$Summary.Add("case_timeout_seconds: $CaseTimeoutSeconds")
$Summary.Add("scope_q_class: $finalClass")
$Summary.Add("next_scope: $nextScope")
$Summary.Add('')
$Summary.Add('Rows:')
foreach ($row in $rowsArray) {
    $Summary.Add("row: condition=$($row.condition) max_iter=$($row.max_iter) restart=$($row.restart) tolerance=$($row.tolerance) direct_status=$($row.direct_status) fmm_status=$($row.fmm_status) direct_iter=$($row.direct_iterations) fmm_iter=$($row.fmm_iterations) direct_true_residual=$($row.direct_true_residual) fmm_true_residual=$($row.fmm_true_residual) energy_delta_abs=$($row.energy_delta_abs_kj_mol) audit_rel=$($row.tree_vs_direct_matvec_rel) speedup=$($row.matvec_speedup) class=$($row.scope_q_class)")
}
$Summary.Add('')
$Summary.Add('classification_note: strict_converged requires both direct and FMM status ok at 1e-6; practical_converged requires the same gates at a looser 1e-4 row.')
$Summary.Add('timeout_policy: timeout_not_verified blocks promotion planning but is not classified as solver failure.')
$Summary.Add('barnase_note: Scope Q diagnoses the bounded Barnase cap-128 convergence gap from Scope P; it does not retry legacy TABI leaf-kernel work.')
$Summary.Add("matrix_csv: $MatrixCsv")
$Summary.Add("decision_csv: $DecisionCsv")
Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

Write-Host "FMM/M2L Scope Q Barnase convergence sweep complete: $DecisionCsv"
