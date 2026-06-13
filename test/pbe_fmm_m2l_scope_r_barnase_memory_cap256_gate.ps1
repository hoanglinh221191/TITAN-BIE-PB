param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 4,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_r')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fmm_m2l_scope_r_matrix.csv'
$MemoryCsv = Join-Path $WorkRoot 'fmm_m2l_scope_r_memory.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_r_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_r_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_R_BARNASE_MEMORY_CAP256_REPORT.txt'

$SurfacePoints = 8
$OpeningTheta = 0.6
$MaxIter = 240
$Restart = 64
$Tolerance = 1.0e-6
$TemperatureK = 298.15
$AuditSampleCount = if ($Quick) { 128 } else { 512 }
$RepeatCount = if ($Quick) { 1 } else { 3 }

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
    throw "Refusing to clean Scope R work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function New-ScopeREmptyCase([string]$Status) {
    return [pscustomobject]@{
        status=$Status; panels=0; iterations=0; final_residual=[double]::NaN
        true_residual=[double]::NaN; energy_kj_mol=[double]::NaN
        matvec_seconds_avg=[double]::NaN; total_seconds=[double]::NaN
        fmm_m2l_mode=''; fmm_near_pairs=0; fmm_far_pairs=0
        fmm_setup_seconds=0.0; fmm_apply_seconds=0.0
        fmm_finite_outputs='false'; tree_vs_direct_matvec_rel=[double]::NaN
        matvec_audit_rows=0; matvec_audit_seconds=0.0
        reaction_energy_source=''; timing_backend=''
        memory_bytes=[double]::NaN
    }
}

function Test-ScopeRFinite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function Get-ScopeREnergyGate([object]$Direct, [bool]$Stress) {
    $magnitude = [Math]::Max(1.0, [Math]::Abs($Direct.energy_kj_mol))
    if ($Stress) { return [Math]::Max(5.0, 0.05 * $magnitude) }
    return [Math]::Max(2.0, 0.05 * $magnitude)
}

function Get-ScopeRSpeedup([object]$Direct, [object]$Fmm) {
    if (-not (Test-ScopeRFinite $Direct.matvec_seconds_avg) -or
        -not (Test-ScopeRFinite $Fmm.matvec_seconds_avg) -or
        $Fmm.matvec_seconds_avg -le 0.0) {
        return [double]::NaN
    }
    return $Direct.matvec_seconds_avg / $Fmm.matvec_seconds_avg
}

function Read-ScopeRGridResult([string]$RunDir) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $timingPath = Join-Path $RunDir 'pbe_timing.csv'
    foreach ($required in @($energyPath, $timingPath)) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Missing grid CLI output: $required"
        }
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $timing = @(Import-Csv -LiteralPath $timingPath)[0]
    return [pscustomobject]@{
        status=[string]$energy.status
        panels=0
        iterations=Convert-FmmCliInt64 $energy.iterations_solvent_final
        final_residual=Convert-FmmCliDouble $energy.rms_solvent
        true_residual=[double]::NaN
        energy_kj_mol=Convert-FmmCliDouble $energy.E_rxn_kJmol
        matvec_seconds_avg=[double]::NaN
        total_seconds=Convert-FmmCliDouble $timing.total_seconds
        fmm_m2l_mode='off'
        fmm_near_pairs=0
        fmm_far_pairs=0
        fmm_setup_seconds=0.0
        fmm_apply_seconds=0.0
        fmm_finite_outputs='true'
        tree_vs_direct_matvec_rel=[double]::NaN
        matvec_audit_rows=0
        matvec_audit_seconds=0.0
        reaction_energy_source=[string]$energy.reaction_energy_source
        timing_backend=[string]$timing.pbe_backend
        memory_bytes=Convert-FmmCliDouble $energy.memory_bytes
    }
}

function Invoke-ScopeRCase(
    [string]$CaseRoot,
    [string]$PqrPath,
    [ValidateSet('grid','direct','fmm')][string]$Backend
) {
    New-FmmCliCleanDir -Path $CaseRoot -WorkRoot $WorkRoot
    if ($Backend -eq 'grid') {
        $arg = Write-FmmCliArgFile -RunDir $CaseRoot -RepoRoot $RepoRoot `
            -PqrPath $PqrPath -Backend grid -Mode direct -FmmMode off `
            -SurfacePoints $SurfacePoints -MaxIter $MaxIter -Restart $Restart `
            -Tolerance $Tolerance -Threads $Threads -SaltM 0.15 `
            -InternalDielectric 1.0 -ExternalDielectric 80.0 `
            -TemperatureK $TemperatureK -OpeningTheta $OpeningTheta `
            -AuditSampleCount 0
    } else {
        $mode = if ($Backend -eq 'fmm') { 'treecode' } else { 'direct' }
        $fmmMode = if ($Backend -eq 'fmm') { 'single-matvec-analytic' } else { 'off' }
        $auditRows = if ($Backend -eq 'fmm') { $AuditSampleCount } else { 0 }
        $arg = Write-FmmCliArgFile -RunDir $CaseRoot -RepoRoot $RepoRoot `
            -PqrPath $PqrPath -Backend tabi -Mode $mode -FmmMode $fmmMode `
            -SurfacePoints $SurfacePoints -MaxIter $MaxIter -Restart $Restart `
            -Tolerance $Tolerance -Threads $Threads -SaltM 0.15 `
            -InternalDielectric 1.0 -ExternalDielectric 80.0 `
            -TemperatureK $TemperatureK -OpeningTheta $OpeningTheta `
            -AuditSampleCount $auditRows
    }
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $CaseRoot `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        return [pscustomobject]@{
            run=$run
            result=(New-ScopeREmptyCase "exit_$($run.exit_code)")
        }
    }
    $result = if ($Backend -eq 'grid') {
        Read-ScopeRGridResult -RunDir $CaseRoot
    } else {
        $expected = if ($Backend -eq 'fmm') { 'single-matvec-analytic' } else { 'off' }
        Read-FmmCliCaseResult -RunDir $CaseRoot -ExpectedFmmMode $expected `
            -RequireAudit ($Backend -eq 'fmm')
    }
    return [pscustomobject]@{ run=$run; result=$result }
}

function Assert-ScopeRGeneratedArgEquals([string]$Root) {
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

function Get-ScopeRFmmGate(
    [object]$Grid,
    [object]$Direct,
    [object]$Fmm,
    [int]$DirectExit,
    [int]$FmmExit,
    [bool]$Stress
) {
    if ($DirectExit -eq 124 -or $FmmExit -eq 124) { return 'timeout_not_verified' }
    if ($DirectExit -ne 0 -or $FmmExit -ne 0) { return 'not_promotion_ready' }
    $fmmPeak = [double]$Fmm.peak_working_set_bytes
    $directPeak = [double]$Direct.peak_working_set_bytes
    $gridPeak = if ($null -ne $Grid) { [double]$Grid.peak_working_set_bytes } else { 0.0 }
    $memoryBlowup = (($directPeak -gt 0 -and $fmmPeak -gt 2.0 * $directPeak) -or
        ($gridPeak -gt 0 -and $fmmPeak -gt 4.0 * $gridPeak))
    if ($Direct.status -ne 'ok' -or $Fmm.status -ne 'ok') {
        if ($memoryBlowup) { return 'memory_regression' }
        return 'not_promotion_ready'
    }
    if ($Fmm.fmm_m2l_mode -ne 'single-matvec-analytic' -or
        $Fmm.fmm_far_pairs -le 0 -or $Fmm.fmm_finite_outputs -ne 'true') {
        return 'not_promotion_ready'
    }
    foreach ($value in @(
        $Direct.energy_kj_mol, $Fmm.energy_kj_mol,
        $Direct.true_residual, $Fmm.true_residual,
        $Direct.matvec_seconds_avg, $Fmm.matvec_seconds_avg,
        $Fmm.tree_vs_direct_matvec_rel)) {
        if (-not (Test-ScopeRFinite $value)) { return 'not_promotion_ready' }
    }
    if ($Fmm.tree_vs_direct_matvec_rel -gt 1.0e-2) { return 'not_promotion_ready' }
    if ([Math]::Abs($Fmm.energy_kj_mol - $Direct.energy_kj_mol) -gt
        (Get-ScopeREnergyGate -Direct $Direct -Stress $Stress)) {
        return 'not_promotion_ready'
    }
    $speedup = Get-ScopeRSpeedup -Direct $Direct -Fmm $Fmm
    if (-not (Test-ScopeRFinite $speedup) -or $speedup -lt 1.2) {
        return 'not_promotion_ready'
    }
    if ($memoryBlowup -and $speedup -lt 1.2) {
        return 'memory_regression'
    }
    if ($Stress) { return 'cap256_pass' }
    return 'cap128_ready'
}

function New-ScopeRRow(
    [string]$Tier,
    [int]$AtomCap,
    [int]$Repeat,
    [string]$Backend,
    [object]$Case,
    [object]$Sample,
    [string]$Gate,
    [double]$Speedup = [double]::NaN,
    [double]$EnergyDeltaAbs = [double]::NaN,
    [double]$MemoryRatioVsGrid = [double]::NaN,
    [double]$MemoryRatioVsDirect = [double]::NaN
) {
    return [pscustomobject]@{
        tier=$Tier; fixture='barnase'; atom_cap=$AtomCap; repeat=$Repeat
        backend=$Backend; input_atoms=$Sample.input_atoms
        sampled_atoms=$Sample.sampled_atoms; surface_points=$SurfacePoints
        opening_theta=$OpeningTheta; max_iter=$MaxIter; restart=$Restart
        tolerance=$Tolerance; audit_sample_count=$AuditSampleCount
        exit_code=$Case.run.exit_code; status=$Case.result.status
        iterations=$Case.result.iterations; final_residual=$Case.result.final_residual
        true_residual=$Case.result.true_residual
        energy_kj_mol=$Case.result.energy_kj_mol
        energy_delta_abs_kj_mol=$EnergyDeltaAbs
        tree_vs_direct_matvec_rel=$Case.result.tree_vs_direct_matvec_rel
        matvec_audit_rows=$Case.result.matvec_audit_rows
        matvec_seconds_avg=$Case.result.matvec_seconds_avg
        matvec_speedup=$Speedup
        wall_seconds=$Case.run.wall_seconds
        peak_working_set_bytes=$Case.run.peak_working_set_bytes
        peak_working_set_mb=$Case.run.peak_working_set_mb
        solver_memory_bytes=$Case.result.memory_bytes
        fmm_m2l_mode=$Case.result.fmm_m2l_mode
        fmm_far_pairs=$Case.result.fmm_far_pairs
        fmm_setup_seconds=$Case.result.fmm_setup_seconds
        fmm_apply_seconds=$Case.result.fmm_apply_seconds
        memory_ratio_vs_grid=$MemoryRatioVsGrid
        memory_ratio_vs_direct=$MemoryRatioVsDirect
        scope_r_gate=$Gate
    }
}

$SamplesRoot = Join-Path $WorkRoot 'samples'
New-FmmCliCleanDir -Path $SamplesRoot -WorkRoot $WorkRoot
$sourcePqr = Resolve-FmmCliFixture -Fixture barnase -RepoRoot $RepoRoot `
    -RunDir $SamplesRoot

$Rows = New-Object System.Collections.Generic.List[object]

if ($null -eq $sourcePqr -or -not (Test-Path -LiteralPath $sourcePqr)) {
    $Rows.Add([pscustomobject]@{
        tier='fixture'; fixture='barnase'; atom_cap=0; repeat=0; backend='none'
        input_atoms=0; sampled_atoms=0; surface_points=0; opening_theta=$OpeningTheta
        max_iter=$MaxIter; restart=$Restart; tolerance=$Tolerance
        audit_sample_count=$AuditSampleCount; exit_code=0
        status='skipped_missing_fixture'; iterations=0
        final_residual=[double]::NaN; true_residual=[double]::NaN
        energy_kj_mol=[double]::NaN; energy_delta_abs_kj_mol=[double]::NaN
        tree_vs_direct_matvec_rel=[double]::NaN; matvec_audit_rows=0
        matvec_seconds_avg=[double]::NaN; matvec_speedup=[double]::NaN
        wall_seconds=[double]::NaN; peak_working_set_bytes=0
        peak_working_set_mb=0.0; solver_memory_bytes=[double]::NaN
        fmm_m2l_mode=''; fmm_far_pairs=0; fmm_setup_seconds=0.0
        fmm_apply_seconds=0.0; memory_ratio_vs_grid=[double]::NaN
        memory_ratio_vs_direct=[double]::NaN; scope_r_gate='skipped_missing_fixture'
    })
} else {
    $cap128Pqr = Join-Path $SamplesRoot 'barnase_cap128.pqr'
    $cap128Sample = Write-FmmCliSampledPqr -InputPath $sourcePqr `
        -OutputPath $cap128Pqr -Cap 128
    for ($repeat = 1; $repeat -le $RepeatCount; ++$repeat) {
        $repeatRoot = Join-Path $WorkRoot "cap128_repeat$repeat"
        New-FmmCliCleanDir -Path $repeatRoot -WorkRoot $WorkRoot
        Write-Host "Scope R cap128 grid/direct/FMM repeat=$repeat"
        $grid = Invoke-ScopeRCase -CaseRoot (Join-Path $repeatRoot 'grid') `
            -PqrPath $cap128Pqr -Backend grid
        $direct = Invoke-ScopeRCase -CaseRoot (Join-Path $repeatRoot 'direct') `
            -PqrPath $cap128Pqr -Backend direct
        $fmm = Invoke-ScopeRCase -CaseRoot (Join-Path $repeatRoot 'fmm') `
            -PqrPath $cap128Pqr -Backend fmm
        $speedup = Get-ScopeRSpeedup -Direct $direct.result -Fmm $fmm.result
        $energyDelta = [Math]::Abs($fmm.result.energy_kj_mol - $direct.result.energy_kj_mol)
        $ratioGrid = if ($grid.run.peak_working_set_bytes -gt 0) {
            $fmm.run.peak_working_set_bytes / [double]$grid.run.peak_working_set_bytes
        } else { [double]::NaN }
        $ratioDirect = if ($direct.run.peak_working_set_bytes -gt 0) {
            $fmm.run.peak_working_set_bytes / [double]$direct.run.peak_working_set_bytes
        } else { [double]::NaN }
        $fmmGate = Get-ScopeRFmmGate -Grid $grid.run -Direct ([pscustomobject]@{
                status=$direct.result.status
                energy_kj_mol=$direct.result.energy_kj_mol
                true_residual=$direct.result.true_residual
                matvec_seconds_avg=$direct.result.matvec_seconds_avg
                peak_working_set_bytes=$direct.run.peak_working_set_bytes
            }) `
            -Fmm ([pscustomobject]@{
                status=$fmm.result.status
                fmm_m2l_mode=$fmm.result.fmm_m2l_mode
                fmm_far_pairs=$fmm.result.fmm_far_pairs
                fmm_finite_outputs=$fmm.result.fmm_finite_outputs
                energy_kj_mol=$fmm.result.energy_kj_mol
                true_residual=$fmm.result.true_residual
                matvec_seconds_avg=$fmm.result.matvec_seconds_avg
                tree_vs_direct_matvec_rel=$fmm.result.tree_vs_direct_matvec_rel
                matvec_audit_rows=$fmm.result.matvec_audit_rows
                peak_working_set_bytes=$fmm.run.peak_working_set_bytes
            }) -DirectExit $direct.run.exit_code -FmmExit $fmm.run.exit_code `
            -Stress $false
        $Rows.Add((New-ScopeRRow -Tier 'cap128' -AtomCap 128 -Repeat $repeat `
            -Backend 'grid' -Case $grid -Sample $cap128Sample -Gate 'grid_reference'))
        $Rows.Add((New-ScopeRRow -Tier 'cap128' -AtomCap 128 -Repeat $repeat `
            -Backend 'direct' -Case $direct -Sample $cap128Sample -Gate 'direct_reference'))
        $Rows.Add((New-ScopeRRow -Tier 'cap128' -AtomCap 128 -Repeat $repeat `
            -Backend 'fmm' -Case $fmm -Sample $cap128Sample -Gate $fmmGate `
            -Speedup $speedup -EnergyDeltaAbs $energyDelta `
            -MemoryRatioVsGrid $ratioGrid -MemoryRatioVsDirect $ratioDirect))
    }

    $cap128FmmRows = @($Rows | Where-Object { $_.tier -eq 'cap128' -and $_.backend -eq 'fmm' })
    $cap128Ready = @($cap128FmmRows | Where-Object { $_.scope_r_gate -ne 'cap128_ready' }).Count -eq 0
    if (-not $Quick -and $cap128Ready) {
        $cap256Pqr = Join-Path $SamplesRoot 'barnase_cap256.pqr'
        $cap256Sample = Write-FmmCliSampledPqr -InputPath $sourcePqr `
            -OutputPath $cap256Pqr -Cap 256
        $stressRoot = Join-Path $WorkRoot 'cap256_stress'
        New-FmmCliCleanDir -Path $stressRoot -WorkRoot $WorkRoot
        Write-Host 'Scope R cap256 direct/FMM stress'
        $direct256 = Invoke-ScopeRCase -CaseRoot (Join-Path $stressRoot 'direct') `
            -PqrPath $cap256Pqr -Backend direct
        $fmm256 = Invoke-ScopeRCase -CaseRoot (Join-Path $stressRoot 'fmm') `
            -PqrPath $cap256Pqr -Backend fmm
        $speedup256 = Get-ScopeRSpeedup -Direct $direct256.result -Fmm $fmm256.result
        $energyDelta256 = [Math]::Abs($fmm256.result.energy_kj_mol - $direct256.result.energy_kj_mol)
        $ratioDirect256 = if ($direct256.run.peak_working_set_bytes -gt 0) {
            $fmm256.run.peak_working_set_bytes / [double]$direct256.run.peak_working_set_bytes
        } else { [double]::NaN }
        $fmm256Gate = Get-ScopeRFmmGate -Grid $null -Direct ([pscustomobject]@{
                status=$direct256.result.status
                energy_kj_mol=$direct256.result.energy_kj_mol
                true_residual=$direct256.result.true_residual
                matvec_seconds_avg=$direct256.result.matvec_seconds_avg
                peak_working_set_bytes=$direct256.run.peak_working_set_bytes
            }) `
            -Fmm ([pscustomobject]@{
                status=$fmm256.result.status
                fmm_m2l_mode=$fmm256.result.fmm_m2l_mode
                fmm_far_pairs=$fmm256.result.fmm_far_pairs
                fmm_finite_outputs=$fmm256.result.fmm_finite_outputs
                energy_kj_mol=$fmm256.result.energy_kj_mol
                true_residual=$fmm256.result.true_residual
                matvec_seconds_avg=$fmm256.result.matvec_seconds_avg
                tree_vs_direct_matvec_rel=$fmm256.result.tree_vs_direct_matvec_rel
                matvec_audit_rows=$fmm256.result.matvec_audit_rows
                peak_working_set_bytes=$fmm256.run.peak_working_set_bytes
            }) -DirectExit $direct256.run.exit_code -FmmExit $fmm256.run.exit_code `
            -Stress $true
        $Rows.Add((New-ScopeRRow -Tier 'cap256' -AtomCap 256 -Repeat 1 `
            -Backend 'direct' -Case $direct256 -Sample $cap256Sample -Gate 'direct_reference'))
        $Rows.Add((New-ScopeRRow -Tier 'cap256' -AtomCap 256 -Repeat 1 `
            -Backend 'fmm' -Case $fmm256 -Sample $cap256Sample -Gate $fmm256Gate `
            -Speedup $speedup256 -EnergyDeltaAbs $energyDelta256 `
            -MemoryRatioVsDirect $ratioDirect256))
    } elseif (-not $Quick) {
        $emptyRun = [pscustomobject]@{
            exit_code=0; timed_out=$false; wall_seconds=[double]::NaN
            peak_working_set_bytes=0; peak_working_set_mb=0.0
        }
        $emptyCase = [pscustomobject]@{
            run=$emptyRun
            result=(New-ScopeREmptyCase 'skipped_cap128_not_ready')
        }
        $emptySample = [pscustomobject]@{ input_atoms=0; sampled_atoms=0 }
        $Rows.Add((New-ScopeRRow -Tier 'cap256' -AtomCap 256 -Repeat 0 `
            -Backend 'fmm' -Case $emptyCase -Sample $emptySample `
            -Gate 'skipped_cap128_not_ready'))
    }
}

$Rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation -Encoding ASCII
$Rows | Select-Object tier,fixture,atom_cap,repeat,backend,status,wall_seconds,
    peak_working_set_bytes,peak_working_set_mb,solver_memory_bytes,
    memory_ratio_vs_grid,memory_ratio_vs_direct,scope_r_gate |
    Export-Csv -LiteralPath $MemoryCsv -NoTypeInformation -Encoding ASCII
Assert-ScopeRGeneratedArgEquals -Root $WorkRoot

$rowsArray = @($Rows.ToArray())
$cap128Fmm = @($rowsArray | Where-Object { $_.tier -eq 'cap128' -and $_.backend -eq 'fmm' })
$cap256Fmm = @($rowsArray | Where-Object { $_.tier -eq 'cap256' -and $_.backend -eq 'fmm' })
$cap128Ready = $cap128Fmm.Count -gt 0 -and
    @($cap128Fmm | Where-Object { $_.scope_r_gate -ne 'cap128_ready' }).Count -eq 0
$cap256Pass = $cap256Fmm.Count -gt 0 -and
    @($cap256Fmm | Where-Object { $_.scope_r_gate -eq 'cap256_pass' }).Count -gt 0
$timeoutRows = @($rowsArray | Where-Object { $_.scope_r_gate -eq 'timeout_not_verified' })
$memoryRegressionRows = @($rowsArray | Where-Object { $_.scope_r_gate -eq 'memory_regression' })
$overallGate = if (@($rowsArray | Where-Object { $_.scope_r_gate -eq 'skipped_missing_fixture' }).Count -gt 0) {
    'skipped_missing_fixture'
} elseif ($memoryRegressionRows.Count -gt 0) {
    'memory_regression'
} elseif ($cap128Ready -and $cap256Pass) {
    'pass'
} elseif ($cap128Ready -and $timeoutRows.Count -gt 0) {
    'cap128_ready_cap256_timeout_not_verified'
} elseif ($cap128Ready -and $cap256Fmm.Count -gt 0) {
    'cap128_ready_cap256_not_promotion_ready'
} elseif ($cap128Ready -and $Quick) {
    'cap128_ready_quick'
} else {
    'not_promotion_ready'
}

$fmmCap128Energy = [double[]]@($cap128Fmm | ForEach-Object { [double]$_.energy_kj_mol })
$energySpread = if ($fmmCap128Energy.Count -gt 0) {
    ($fmmCap128Energy | Measure-Object -Maximum).Maximum -
        ($fmmCap128Energy | Measure-Object -Minimum).Minimum
} else { [double]::NaN }
$decision = [pscustomobject]@{
    scope_r_gate=$overallGate
    quick_mode=$Quick.IsPresent
    rows=$rowsArray.Count
    cap128_fmm_rows=$cap128Fmm.Count
    cap256_fmm_rows=$cap256Fmm.Count
    timeout_rows=$timeoutRows.Count
    cap128_energy_spread_kj_mol=$energySpread
    production_default='grid'
    fmm_status='hidden_experimental'
    next_scope=if ($overallGate -eq 'pass') {
        'consider_broader_hidden_regression_not_public_exposure'
    } elseif ($overallGate -eq 'cap128_ready_cap256_timeout_not_verified') {
        'reduce_cap256_direct_verifier_or_extend_timeout'
    } elseif ($overallGate -eq 'cap128_ready_quick') {
        'run_full_scope_r'
    } else {
        'hold_and_diagnose_scope_r'
    }
}
$decision | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation -Encoding ASCII

$Summary = New-Object System.Collections.Generic.List[string]
$Summary.Add('TITAN FMM/M2L Scope R Barnase Fair-Memory + Cap-256 Gate')
$Summary.Add('==========================================================')
$Summary.Add('')
$Summary.Add("quick_mode: $($Quick.IsPresent)")
$Summary.Add('scope: hidden_cli_bounded_barnase_memory_cap256')
$Summary.Add('hidden_option: -pbe-tabi-fmm-m2l=single-matvec-analytic')
$Summary.Add('production_default: grid')
$Summary.Add('fmm_status: hidden_experimental')
$Summary.Add("max_iter: $MaxIter")
$Summary.Add("restart: $Restart")
$Summary.Add("tolerance: $Tolerance")
$Summary.Add("surface_points: $SurfacePoints")
$Summary.Add("opening_theta: $OpeningTheta")
$Summary.Add("audit_sample_count: $AuditSampleCount")
$Summary.Add("case_timeout_seconds: $CaseTimeoutSeconds")
$Summary.Add("scope_r_gate: $overallGate")
$Summary.Add("cap128_energy_spread_kj_mol: $energySpread")
$Summary.Add('')
$Summary.Add('Rows:')
foreach ($row in $rowsArray) {
    $Summary.Add("row: tier=$($row.tier) cap=$($row.atom_cap) repeat=$($row.repeat) backend=$($row.backend) status=$($row.status) iter=$($row.iterations) residual=$($row.true_residual) energy=$($row.energy_kj_mol) energy_delta=$($row.energy_delta_abs_kj_mol) audit=$($row.tree_vs_direct_matvec_rel) speedup=$($row.matvec_speedup) wall=$($row.wall_seconds) peak_mb=$($row.peak_working_set_mb) ratio_grid=$($row.memory_ratio_vs_grid) ratio_direct=$($row.memory_ratio_vs_direct) gate=$($row.scope_r_gate)")
}
$Summary.Add('')
$Summary.Add('memory_note: peak_working_set_bytes is process-level memory from the launched titan.exe process; solver_memory_bytes remains the internal solver estimate.')
$Summary.Add('timeout_policy: timeout_not_verified blocks promotion planning but is not classified as solver failure.')
$Summary.Add('barnase_note: Scope R uses the Scope Q converged preset for bounded Barnase and does not retry legacy TABI leaf-kernel paths.')
$Summary.Add("matrix_csv: $MatrixCsv")
$Summary.Add("memory_csv: $MemoryCsv")
$Summary.Add("decision_csv: $DecisionCsv")
Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

Write-Host "FMM/M2L Scope R Barnase memory/cap256 gate complete: $DecisionCsv"
