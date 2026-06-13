param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 600,
    [int]$MaxIter = 120,
    [int]$Restart = 24,
    [double]$Tolerance = 1.0e-6,
    [int]$Threads = 4,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_o')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$RepeatCsv = Join-Path $WorkRoot 'fmm_m2l_scope_o_repeatability.csv'
$StressCsv = Join-Path $WorkRoot 'fmm_m2l_scope_o_stress.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_o_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_o_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_O_REPORT.txt'
$TemperatureK = 298.15
$OpeningTheta = 0.6
$AuditSampleCount = if ($Quick) { 128 } else { 512 }
$RepeatCount = if ($Quick) { 2 } else { 3 }
if ($Quick) {
    $MaxIter = [Math]::Min($MaxIter, 80)
    $CaseTimeoutSeconds = [Math]::Min($CaseTimeoutSeconds, 240)
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope O work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function New-EmptyFmmCase([string]$Status) {
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

function Get-MedianDouble([double[]]$Values) {
    $finite = @($Values | Where-Object { Test-FmmCliFinite $_ } | Sort-Object)
    if ($finite.Count -eq 0) { return [double]::NaN }
    $mid = [int]($finite.Count / 2)
    if (($finite.Count % 2) -eq 0) {
        return 0.5 * ($finite[$mid - 1] + $finite[$mid])
    }
    return $finite[$mid]
}

function Get-ItemCount($Value) {
    if ($null -eq $Value) {
        return 0
    }
    if ($Value -is [System.Collections.ICollection]) {
        return $Value.Count
    }
    return @($Value).Count
}

function Get-ScopeOEnergyGate([object]$Direct, [string]$Fixture,
                              [bool]$Stress) {
    $magnitude = [Math]::Max(1.0, [Math]::Abs($Direct.energy_kj_mol))
    if ($Stress -and $Fixture -eq 'research-large') {
        return [Math]::Max(10.0, 0.08 * $magnitude)
    }
    if ($Stress) {
        return [Math]::Max(5.0, 0.05 * $magnitude)
    }
    if ($Fixture -eq 'research-large') {
        return [Math]::Max(5.0, 0.08 * $magnitude)
    }
    return [Math]::Max(2.0, 0.05 * $magnitude)
}

function Get-ScopeOCaseGate(
    [int]$DirectExit,
    [int]$FmmExit,
    [AllowNull()][object]$Direct,
    [AllowNull()][object]$Fmm,
    [string]$Fixture,
    [bool]$Stress,
    [bool]$RequireSpeed
) {
    if ($DirectExit -eq 124 -or $FmmExit -eq 124) {
        return 'timeout_not_verified'
    }
    if ($DirectExit -ne 0 -or $FmmExit -ne 0 -or
        $null -eq $Direct -or $null -eq $Fmm) {
        return 'status_fail'
    }
    if ($Direct.status -ne 'ok' -or $Fmm.status -ne 'ok') {
        return 'status_fail'
    }
    if ($Fmm.fmm_m2l_mode -ne 'single-matvec-analytic' -or
        $Fmm.fmm_far_pairs -le 0 -or $Fmm.fmm_finite_outputs -ne 'true') {
        return 'telemetry_fail'
    }
    foreach ($value in @(
        $Direct.energy_kj_mol, $Fmm.energy_kj_mol,
        $Direct.true_residual, $Fmm.true_residual,
        $Direct.matvec_seconds_avg, $Fmm.matvec_seconds_avg,
        $Fmm.tree_vs_direct_matvec_rel)) {
        if (-not (Test-FmmCliFinite $value)) { return 'nonfinite_fail' }
    }
    $residualGate = [Math]::Max(
        1.0e-4, 10.0 * [Math]::Max(1.0e-12, [Math]::Abs($Direct.true_residual)))
    if ([Math]::Abs($Fmm.true_residual) -gt $residualGate) {
        return 'convergence_fail'
    }
    if ($Fmm.matvec_audit_rows -le 0 -or
        $Fmm.tree_vs_direct_matvec_rel -gt 1.0e-2) {
        return 'accuracy_fail'
    }
    if ([Math]::Abs($Fmm.energy_kj_mol - $Direct.energy_kj_mol) -gt
        (Get-ScopeOEnergyGate -Direct $Direct -Fixture $Fixture -Stress $Stress)) {
        return 'energy_fail'
    }
    $speedup = $Direct.matvec_seconds_avg / $Fmm.matvec_seconds_avg
    if ($RequireSpeed -and $speedup -lt 1.2) {
        return 'speed_fail'
    }
    return 'pass'
}

function Assert-GeneratedArgEquals([string]$Root) {
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

function Invoke-ScopeOCase(
    [string]$CaseRoot,
    [string]$PqrPath,
    [string]$Mode,
    [string]$FmmMode,
    [int]$AtomCap,
    [int]$SurfacePoints,
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
        New-EmptyFmmCase "exit_$($run.exit_code)"
    }
    return [pscustomobject]@{
        run = $run
        result = $result
    }
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$SamplesRoot = Join-Path $WorkRoot 'samples'
New-FmmCliCleanDir -Path $SamplesRoot -WorkRoot $WorkRoot
$RepeatFixtures = if ($Quick) {
    @([pscustomobject]@{ fixture='research-medium'; atom_cap=64; surface_points=8 })
} else {
    @(
        [pscustomobject]@{ fixture='research-medium'; atom_cap=128; surface_points=8 },
        [pscustomobject]@{ fixture='research-large'; atom_cap=128; surface_points=8 },
        [pscustomobject]@{ fixture='research-alt'; atom_cap=128; surface_points=8 }
    )
}
$StressFixtures = if ($Quick) {
    @()
} else {
    @(
        [pscustomobject]@{ fixture='research-medium'; atom_cap=256; surface_points=8 },
        [pscustomobject]@{ fixture='research-large'; atom_cap=256; surface_points=8 }
    )
}

$RepeatRows = New-Object System.Collections.Generic.List[object]
$RepeatAggregateRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureSpec in $RepeatFixtures) {
    $fixture = $fixtureSpec.fixture
    $sampleKey = "$fixture`_cap$($fixtureSpec.atom_cap)"
    $sourcePqr = Resolve-FmmCliFixture -Fixture $fixture -RepoRoot $RepoRoot
    $samplePqr = Join-Path $SamplesRoot "$sampleKey.pqr"
    $sample = Write-FmmCliSampledPqr -InputPath $sourcePqr -OutputPath $samplePqr `
        -Cap $fixtureSpec.atom_cap
    $fixtureRoot = Join-Path $WorkRoot "repeat_$sampleKey"
    New-FmmCliCleanDir -Path $fixtureRoot -WorkRoot $WorkRoot
    Write-Host "Scope O repeat direct: $fixture cap=$($fixtureSpec.atom_cap)"
    $directCase = Invoke-ScopeOCase -CaseRoot (Join-Path $fixtureRoot 'direct') `
        -PqrPath $samplePqr -Mode direct -FmmMode off `
        -AtomCap $fixtureSpec.atom_cap -SurfacePoints $fixtureSpec.surface_points `
        -AuditRows 0
    for ($repeat = 1; $repeat -le $RepeatCount; ++$repeat) {
        Write-Host "Scope O repeat FMM: $fixture cap=$($fixtureSpec.atom_cap) repeat=$repeat"
        $fmmCase = Invoke-ScopeOCase -CaseRoot (Join-Path $fixtureRoot "fmm_$repeat") `
            -PqrPath $samplePqr -Mode treecode `
            -FmmMode single-matvec-analytic -AtomCap $fixtureSpec.atom_cap `
            -SurfacePoints $fixtureSpec.surface_points -AuditRows $AuditSampleCount
        $gate = Get-ScopeOCaseGate -DirectExit $directCase.run.exit_code `
            -FmmExit $fmmCase.run.exit_code -Direct $directCase.result `
            -Fmm $fmmCase.result -Fixture $fixture -Stress $false `
            -RequireSpeed $false
        $delta = $fmmCase.result.energy_kj_mol - $directCase.result.energy_kj_mol
        $speedup = $directCase.result.matvec_seconds_avg /
            $fmmCase.result.matvec_seconds_avg
        $RepeatRows.Add([pscustomobject]@{
            tier='repeatability'; fixture=$fixture; repeat=$repeat
            atom_cap=$fixtureSpec.atom_cap; input_atoms=$sample.input_atoms
            sampled_atoms=$sample.sampled_atoms
            surface_points=$fixtureSpec.surface_points
            opening_theta=$OpeningTheta; audit_sample_count=$AuditSampleCount
            direct_exit_code=$directCase.run.exit_code
            fmm_exit_code=$fmmCase.run.exit_code
            direct_status=$directCase.result.status
            fmm_status=$fmmCase.result.status
            direct_panels=$directCase.result.panels
            fmm_panels=$fmmCase.result.panels
            direct_energy_kj_mol=$directCase.result.energy_kj_mol
            fmm_energy_kj_mol=$fmmCase.result.energy_kj_mol
            energy_delta_abs_kj_mol=[Math]::Abs($delta)
            direct_true_residual=$directCase.result.true_residual
            fmm_true_residual=$fmmCase.result.true_residual
            tree_vs_direct_matvec_rel=$fmmCase.result.tree_vs_direct_matvec_rel
            matvec_audit_rows=$fmmCase.result.matvec_audit_rows
            direct_matvec_seconds_avg=$directCase.result.matvec_seconds_avg
            fmm_matvec_seconds_avg=$fmmCase.result.matvec_seconds_avg
            matvec_speedup=$speedup
            fmm_m2l_mode=$fmmCase.result.fmm_m2l_mode
            fmm_far_pairs=$fmmCase.result.fmm_far_pairs
            fmm_wall_seconds=$fmmCase.run.wall_seconds
            repeat_gate=$gate
        })
    }
    $fixtureRows = @($RepeatRows | Where-Object { $_.fixture -eq $fixture })
    $energies = [double[]]@($fixtureRows | ForEach-Object { [double]$_.fmm_energy_kj_mol })
    $speedups = [double[]]@($fixtureRows | ForEach-Object { [double]$_.matvec_speedup })
    $audits = [double[]]@($fixtureRows | ForEach-Object { [double]$_.tree_vs_direct_matvec_rel })
    $energyMean = if ($energies.Count -gt 0) {
        ($energies | Measure-Object -Average).Average
    } else { [double]::NaN }
    $energySpread = if ($energies.Count -gt 0) {
        ($energies | Measure-Object -Maximum).Maximum -
            ($energies | Measure-Object -Minimum).Minimum
    } else { [double]::NaN }
    $energySpreadGate = [Math]::Max(0.5, 1.0e-4 * [Math]::Abs($energyMean))
    $auditMax = if ($audits.Count -gt 0) {
        ($audits | Measure-Object -Maximum).Maximum
    } else { [double]::NaN }
    $speedupMedian = Get-MedianDouble $speedups
    $speedupMin = if ($speedups.Count -gt 0) {
        ($speedups | Measure-Object -Minimum).Minimum
    } else { [double]::NaN }
    $rowFailures = @($fixtureRows | Where-Object { $_.repeat_gate -ne 'pass' })
    $aggregateGate = if ($rowFailures.Count -gt 0) {
        'row_fail'
    } elseif ($energySpread -gt $energySpreadGate) {
        'energy_repeatability_fail'
    } elseif ($auditMax -gt 1.0e-2) {
        'audit_repeatability_fail'
    } elseif ($speedupMedian -lt 1.2 -or $speedupMin -lt 1.0) {
        'speed_repeatability_fail'
    } else {
        'pass'
    }
    $RepeatAggregateRows.Add([pscustomobject]@{
        fixture=$fixture; atom_cap=$fixtureSpec.atom_cap
        repeats=$RepeatCount; fmm_energy_mean_kj_mol=$energyMean
        fmm_energy_spread_kj_mol=$energySpread
        energy_spread_gate_kj_mol=$energySpreadGate
        audit_max=$auditMax; speedup_median=$speedupMedian
        speedup_min=$speedupMin; repeatability_gate=$aggregateGate
    })
}

$StressRows = New-Object System.Collections.Generic.List[object]
$MediumStressPassed = $Quick
foreach ($fixtureSpec in $StressFixtures) {
    if ($fixtureSpec.fixture -eq 'research-large' -and -not $MediumStressPassed) {
        $StressRows.Add([pscustomobject]@{
            tier='stress'; fixture=$fixtureSpec.fixture; atom_cap=$fixtureSpec.atom_cap
            sampled_atoms=0; surface_points=$fixtureSpec.surface_points
            direct_exit_code=0; fmm_exit_code=0; direct_status='skipped'
            fmm_status='skipped_medium_stress_not_passed'; direct_panels=0
            fmm_panels=0; direct_energy_kj_mol=[double]::NaN
            fmm_energy_kj_mol=[double]::NaN; energy_delta_abs_kj_mol=[double]::NaN
            direct_true_residual=[double]::NaN; fmm_true_residual=[double]::NaN
            tree_vs_direct_matvec_rel=[double]::NaN; matvec_audit_rows=0
            direct_matvec_seconds_avg=[double]::NaN
            fmm_matvec_seconds_avg=[double]::NaN; matvec_speedup=[double]::NaN
            fmm_m2l_mode=''; fmm_far_pairs=0; stress_gate='skipped'
        })
        continue
    }
    $fixture = $fixtureSpec.fixture
    $sampleKey = "$fixture`_cap$($fixtureSpec.atom_cap)"
    $sourcePqr = Resolve-FmmCliFixture -Fixture $fixture -RepoRoot $RepoRoot
    $samplePqr = Join-Path $SamplesRoot "$sampleKey.pqr"
    $sample = Write-FmmCliSampledPqr -InputPath $sourcePqr -OutputPath $samplePqr `
        -Cap $fixtureSpec.atom_cap
    $fixtureRoot = Join-Path $WorkRoot "stress_$sampleKey"
    New-FmmCliCleanDir -Path $fixtureRoot -WorkRoot $WorkRoot
    Write-Host "Scope O stress direct/FMM: $fixture cap=$($fixtureSpec.atom_cap)"
    $directCase = Invoke-ScopeOCase -CaseRoot (Join-Path $fixtureRoot 'direct') `
        -PqrPath $samplePqr -Mode direct -FmmMode off `
        -AtomCap $fixtureSpec.atom_cap -SurfacePoints $fixtureSpec.surface_points `
        -AuditRows 0
    $fmmCase = Invoke-ScopeOCase -CaseRoot (Join-Path $fixtureRoot 'fmm') `
        -PqrPath $samplePqr -Mode treecode `
        -FmmMode single-matvec-analytic -AtomCap $fixtureSpec.atom_cap `
        -SurfacePoints $fixtureSpec.surface_points -AuditRows $AuditSampleCount
    $gate = Get-ScopeOCaseGate -DirectExit $directCase.run.exit_code `
        -FmmExit $fmmCase.run.exit_code -Direct $directCase.result `
        -Fmm $fmmCase.result -Fixture $fixture -Stress $true `
        -RequireSpeed $true
    if ($fixture -eq 'research-medium') { $MediumStressPassed = $gate -eq 'pass' }
    $delta = $fmmCase.result.energy_kj_mol - $directCase.result.energy_kj_mol
    $speedup = $directCase.result.matvec_seconds_avg /
        $fmmCase.result.matvec_seconds_avg
    $StressRows.Add([pscustomobject]@{
        tier='stress'; fixture=$fixture; atom_cap=$fixtureSpec.atom_cap
        sampled_atoms=$sample.sampled_atoms; surface_points=$fixtureSpec.surface_points
        direct_exit_code=$directCase.run.exit_code
        fmm_exit_code=$fmmCase.run.exit_code
        direct_status=$directCase.result.status
        fmm_status=$fmmCase.result.status
        direct_panels=$directCase.result.panels
        fmm_panels=$fmmCase.result.panels
        direct_energy_kj_mol=$directCase.result.energy_kj_mol
        fmm_energy_kj_mol=$fmmCase.result.energy_kj_mol
        energy_delta_abs_kj_mol=[Math]::Abs($delta)
        direct_true_residual=$directCase.result.true_residual
        fmm_true_residual=$fmmCase.result.true_residual
        tree_vs_direct_matvec_rel=$fmmCase.result.tree_vs_direct_matvec_rel
        matvec_audit_rows=$fmmCase.result.matvec_audit_rows
        direct_matvec_seconds_avg=$directCase.result.matvec_seconds_avg
        fmm_matvec_seconds_avg=$fmmCase.result.matvec_seconds_avg
        matvec_speedup=$speedup
        fmm_m2l_mode=$fmmCase.result.fmm_m2l_mode
        fmm_far_pairs=$fmmCase.result.fmm_far_pairs
        stress_gate=$gate
    })
}

$RepeatRows | Export-Csv -LiteralPath $RepeatCsv -NoTypeInformation -Encoding ASCII
$StressRows | Export-Csv -LiteralPath $StressCsv -NoTypeInformation -Encoding ASCII

Assert-GeneratedArgEquals -Root $WorkRoot

$repeatFailures = @($RepeatRows | Where-Object { $_.repeat_gate -ne 'pass' })
$aggregateFailures = @($RepeatAggregateRows | Where-Object {
    $_.repeatability_gate -ne 'pass'
})
$requiredStressRows = if ($Quick) { @() } else { @($StressRows | Where-Object {
    $_.stress_gate -ne 'skipped'
}) }
$stressFailures = @($requiredStressRows | Where-Object { $_.stress_gate -ne 'pass' })
$timeoutRows = @($RepeatRows | Where-Object { $_.repeat_gate -eq 'timeout_not_verified' }) +
    @($StressRows | Where-Object { $_.stress_gate -eq 'timeout_not_verified' })
$overallGate = if ($repeatFailures.Count -eq 0 -and
    $aggregateFailures.Count -eq 0 -and $stressFailures.Count -eq 0 -and
    ($Quick -or $requiredStressRows.Count -eq 2)) {
    'pass'
} elseif ($timeoutRows.Count -gt 0 -and
    ($timeoutRows.Count -eq ($repeatFailures.Count + $stressFailures.Count))) {
    'timeout_not_verified'
} else {
    'fail'
}
$failureClasses = @(
    $repeatFailures | ForEach-Object { $_.repeat_gate }
    $aggregateFailures | ForEach-Object { $_.repeatability_gate }
    $stressFailures | ForEach-Object { $_.stress_gate }
) | Sort-Object -Unique
$decision = [pscustomobject]@{
    scope_o_gate=$overallGate
    quick_mode=$Quick.IsPresent
    repeat_rows=(Get-ItemCount $RepeatRows)
    repeatability_fixtures=(Get-ItemCount $RepeatAggregateRows)
    stress_rows=(Get-ItemCount $requiredStressRows)
    timeout_rows=(Get-ItemCount $timeoutRows)
    failure_classes=($failureClasses -join ';')
    production_default='grid'
    fmm_status='hidden_experimental'
    next_scope=if ($overallGate -eq 'pass' -and $Quick) {
        'run_full_scope_o_gate'
    } elseif ($overallGate -eq 'pass') {
        'scope_p_hidden_longer_integration_regression_gate'
    } else {
        'hold_and_diagnose_scope_o_failures'
    }
}
$decision | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation -Encoding ASCII

$Summary = New-Object System.Collections.Generic.List[string]
$Summary.Add('TITAN FMM/M2L Scope O Repeatability And Cap-256 Stress Gate')
$Summary.Add('==========================================================')
$Summary.Add('')
$Summary.Add("quick_mode: $($Quick.IsPresent)")
$Summary.Add('scope: hidden_cli_repeatability_cap256_stress')
$Summary.Add('hidden_option: -pbe-tabi-fmm-m2l=single-matvec-analytic')
$Summary.Add('production_default: grid')
$Summary.Add('fmm_status: hidden_experimental')
$Summary.Add("opening_theta: $OpeningTheta")
$Summary.Add("audit_sample_count: $AuditSampleCount")
$Summary.Add("repeat_count: $RepeatCount")
$Summary.Add("case_timeout_seconds: $CaseTimeoutSeconds")
$Summary.Add("scope_o_gate: $overallGate")
$Summary.Add("failure_classes: $($decision.failure_classes)")
$Summary.Add('')
$Summary.Add('Repeatability aggregates:')
foreach ($row in $RepeatAggregateRows) {
    $Summary.Add("repeatability: fixture=$($row.fixture) atom_cap=$($row.atom_cap) repeats=$($row.repeats) energy_spread=$($row.fmm_energy_spread_kj_mol) energy_gate=$($row.energy_spread_gate_kj_mol) audit_max=$($row.audit_max) speedup_median=$($row.speedup_median) speedup_min=$($row.speedup_min) gate=$($row.repeatability_gate)")
}
$Summary.Add('')
$Summary.Add('Stress rows:')
foreach ($row in $StressRows) {
    $Summary.Add("stress: fixture=$($row.fixture) atom_cap=$($row.atom_cap) sampled_atoms=$($row.sampled_atoms) panels=$($row.fmm_panels) energy_delta_abs=$($row.energy_delta_abs_kj_mol) fmm_residual=$($row.fmm_true_residual) audit_rel=$($row.tree_vs_direct_matvec_rel) speedup=$($row.matvec_speedup) gate=$($row.stress_gate)")
}
$Summary.Add('')
$Summary.Add('acceptance: repeat rows must pass Scope N-style gates; per-fixture energy spread, audit max, median speedup, and min speedup must pass; full mode also requires cap-256 medium and gated large stress rows to pass.')
$Summary.Add('timeout_policy: timeout_not_verified blocks promotion planning but is not classified as a solver failure.')
$Summary.Add('next_gate: Scope P hidden longer integration/regression only after full Scope O passes.')
$Summary.Add("repeatability_csv: $RepeatCsv")
$Summary.Add("stress_csv: $StressCsv")
$Summary.Add("decision_csv: $DecisionCsv")
Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

if ($overallGate -ne 'pass') {
    throw "Scope O gate did not pass: $overallGate ($($decision.failure_classes))"
}
Write-Host "FMM/M2L Scope O repeatability probe complete: $DecisionCsv"
