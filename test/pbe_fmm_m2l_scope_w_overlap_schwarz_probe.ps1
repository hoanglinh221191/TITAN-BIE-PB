param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 2400,
    [int]$Threads = 4,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_w')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fmm_m2l_scope_w_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_w_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_w_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_W_OVERLAP_SCHWARZ_REPORT.txt'

$SurfacePoints = 8
$OpeningTheta = 0.6
$TemperatureK = 298.15
$AuditSampleCount = if ($Quick) { 128 } else { 512 }

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
    throw "Refusing to clean Scope W work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Test-ScopeWFinite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function New-ScopeWEmptyResult([string]$Status) {
    return [pscustomobject]@{
        status=$Status; panels=0; iterations=0
        true_residual=[double]::NaN; energy_kj_mol=[double]::NaN
        matvec_seconds_avg=[double]::NaN; tree_vs_direct_matvec_rel=[double]::NaN
        fmm_m2l_mode=''; fmm_far_pairs=0; fmm_finite_outputs='false'
    }
}

function Read-ScopeWPrecondProfile([string]$RunDir) {
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $profilePath)) {
        return [pscustomobject]@{
            requested='NA'; effective='NA'; block_cap=0; overlap_cap=0
            estimated_bytes=[double]::NaN; valid_blocks=0; valid_panels=0
            quality_status='NA'; setup_seconds=[double]::NaN
            apply_seconds=[double]::NaN; apply_calls=0
        }
    }
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    $overlapCap = 0
    if ($profile.PSObject.Properties.Name -contains 'preconditioner_overlap_panel_cap') {
        $overlapCap = Convert-FmmCliInt64 $profile.preconditioner_overlap_panel_cap
    }
    return [pscustomobject]@{
        requested=[string]$profile.preconditioner_requested
        effective=[string]$profile.preconditioner_effective
        block_cap=Convert-FmmCliInt64 $profile.preconditioner_block_panel_cap
        overlap_cap=$overlapCap
        estimated_bytes=Convert-FmmCliDouble $profile.preconditioner_estimated_bytes
        valid_blocks=Convert-FmmCliInt64 $profile.preconditioner_valid_blocks
        valid_panels=Convert-FmmCliInt64 $profile.preconditioner_valid_panels
        quality_status=[string]$profile.preconditioner_quality_status
        setup_seconds=Convert-FmmCliDouble $profile.preconditioner_setup_seconds
        apply_seconds=Convert-FmmCliDouble $profile.preconditioner_apply_seconds
        apply_calls=Convert-FmmCliInt64 $profile.preconditioner_apply_calls
    }
}

function Invoke-ScopeWCase(
    [string]$CaseRoot,
    [string]$PqrPath,
    [ValidateSet('direct','fmm-block','fmm-overlap')][string]$Kind,
    [int]$MaxIter,
    [int]$Restart,
    [double]$Tolerance,
    [int]$BlockCap,
    [int]$OverlapCap
) {
    New-FmmCliCleanDir -Path $CaseRoot -WorkRoot $WorkRoot
    $mode = if ($Kind -eq 'direct') { 'direct' } else { 'treecode' }
    $fmmMode = if ($Kind -eq 'direct') { 'off' } else { 'single-matvec-analytic' }
    $precond = switch ($Kind) {
        'fmm-block' { 'fmm-block-jacobi' }
        'fmm-overlap' { 'fmm-overlap-schwarz' }
        default { 'diagonal' }
    }
    $auditRows = if ($Kind -eq 'direct') { 0 } else { $AuditSampleCount }
    $arg = Write-FmmCliArgFile -RunDir $CaseRoot -RepoRoot $RepoRoot `
        -PqrPath $PqrPath -Backend tabi -Mode $mode -FmmMode $fmmMode `
        -SurfacePoints $SurfacePoints -MaxIter $MaxIter -Restart $Restart `
        -Tolerance $Tolerance -Threads $Threads -SaltM 0.15 `
        -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK $TemperatureK -OpeningTheta $OpeningTheta `
        -AuditSampleCount $auditRows -Preconditioner $precond `
        -PreconditionerBlockPanelCap $BlockCap `
        -PreconditionerOverlapPanelCap $OverlapCap
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $CaseRoot `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $result = if ($run.exit_code -eq 0) {
        $expected = if ($Kind -eq 'direct') { 'off' } else { 'single-matvec-analytic' }
        Read-FmmCliCaseResult -RunDir $CaseRoot -ExpectedFmmMode $expected `
            -RequireAudit ($Kind -ne 'direct')
    } else {
        New-ScopeWEmptyResult "exit_$($run.exit_code)"
    }
    return [pscustomobject]@{
        run=$run; result=$result; preconditioner=(Read-ScopeWPrecondProfile -RunDir $CaseRoot)
    }
}

function Invoke-ScopeWParserProbe(
    [string]$Name,
    [string]$PqrPath,
    [string]$Preconditioner,
    [string]$FmmMode,
    [int]$OverlapCap,
    [bool]$ExpectSuccess
) {
    $runDir = Join-Path $workRootFull "parser_$Name"
    New-FmmCliCleanDir -Path $runDir -WorkRoot $WorkRoot
    $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
        -PqrPath $PqrPath -Backend tabi -Mode treecode -FmmMode $FmmMode `
        -SurfacePoints $SurfacePoints -MaxIter 1 -Restart 8 `
        -Tolerance 1.0e-4 -Threads 1 -SaltM 0.15 `
        -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK $TemperatureK -OpeningTheta $OpeningTheta `
        -AuditSampleCount 0 -Preconditioner $Preconditioner `
        -PreconditionerBlockPanelCap 32 `
        -PreconditionerOverlapPanelCap $OverlapCap
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds 120
    if ($ExpectSuccess -and $run.exit_code -ne 0) {
        throw "Scope W parser probe $Name expected success, got exit $($run.exit_code)"
    }
    if ((-not $ExpectSuccess) -and $run.exit_code -eq 0) {
        throw "Scope W parser probe $Name expected rejection but succeeded"
    }
}

function Get-ScopeWSpeedup([object]$Direct, [object]$Candidate) {
    if ($null -eq $Direct) { return [double]::NaN }
    if (-not (Test-ScopeWFinite $Direct.data.result.matvec_seconds_avg) -or
        -not (Test-ScopeWFinite $Candidate.data.result.matvec_seconds_avg) -or
        $Candidate.data.result.matvec_seconds_avg -le 0.0) {
        return [double]::NaN
    }
    return $Direct.data.result.matvec_seconds_avg / $Candidate.data.result.matvec_seconds_avg
}

function Get-ScopeWEnergyDelta([object]$Direct, [object]$Candidate) {
    if ($null -eq $Direct) { return [double]::NaN }
    if (-not (Test-ScopeWFinite $Direct.data.result.energy_kj_mol) -or
        -not (Test-ScopeWFinite $Candidate.data.result.energy_kj_mol)) {
        return [double]::NaN
    }
    return [Math]::Abs($Candidate.data.result.energy_kj_mol -
        $Direct.data.result.energy_kj_mol)
}

function Get-ScopeWClass([object]$Candidate, [object]$Baseline, [object]$Direct) {
    if ($Candidate.data.run.exit_code -eq 124) { return 'timeout_not_verified' }
    if ($Candidate.data.run.exit_code -ne 0) { return 'status_fail' }
    if ($Candidate.backend -eq 'fmm_overlap_schwarz') {
        if ($Candidate.data.preconditioner.effective -notin @('fmm_overlap_schwarz','fmm_overlap_schwarz_partial')) {
            return 'preconditioner_fail'
        }
        if ($Candidate.data.preconditioner.overlap_cap -ne $Candidate.overlap_cap) {
            return 'overlap_telemetry_fail'
        }
    }
    if ($Candidate.backend -ne 'direct' -and
        ($Candidate.data.result.fmm_m2l_mode -ne 'single-matvec-analytic' -or
         $Candidate.data.result.fmm_far_pairs -le 0 -or
         $Candidate.data.result.fmm_finite_outputs -ne 'true')) {
        return 'telemetry_fail'
    }
    if ($Candidate.backend -ne 'direct' -and
        (-not (Test-ScopeWFinite $Candidate.data.result.tree_vs_direct_matvec_rel) -or
         $Candidate.data.result.tree_vs_direct_matvec_rel -gt 1.0e-2)) {
        return 'accuracy_limited'
    }
    if ($Candidate.data.result.status -eq 'ok' -and
        (Test-ScopeWFinite $Candidate.data.result.true_residual) -and
        $Candidate.data.result.true_residual -le 1.0e-6) {
        if ($null -ne $Direct) {
            $speedup = Get-ScopeWSpeedup -Direct $Direct -Candidate $Candidate
            $energyDelta = Get-ScopeWEnergyDelta -Direct $Direct -Candidate $Candidate
            $energyGate = [Math]::Max(5.0, 0.05 * [Math]::Abs($Direct.data.result.energy_kj_mol))
            if ((Test-ScopeWFinite $speedup) -and $speedup -ge 1.2 -and
                (Test-ScopeWFinite $energyDelta) -and $energyDelta -le $energyGate) {
                return 'strict_converged_verified'
            }
            return 'strict_converged_direct_gate_fail'
        }
        return 'strict_fmm_pending_direct'
    }
    if ($null -ne $Baseline -and
        (Test-ScopeWFinite $Baseline.data.result.true_residual) -and
        (Test-ScopeWFinite $Candidate.data.result.true_residual) -and
        $Candidate.data.result.true_residual -gt 0.0) {
        $improvement = $Baseline.data.result.true_residual /
            $Candidate.data.result.true_residual
        if ($improvement -ge 4.0) { return 'residual_improved_4x' }
        if ($improvement -gt 1.0) { return 'residual_improved' }
    }
    return 'not_improved'
}

function New-ScopeWRecord(
    [string]$Condition,
    [int]$AtomCap,
    [string]$Backend,
    [int]$BlockCap,
    [int]$OverlapCap,
    [int]$MaxIter,
    [int]$Restart,
    [double]$Tolerance,
    [object]$Data
) {
    return [pscustomobject]@{
        condition=$Condition; atom_cap=$AtomCap; backend=$Backend
        block_cap=$BlockCap; overlap_cap=$OverlapCap; max_iter=$MaxIter
        restart=$Restart; tolerance=$Tolerance; data=$Data
    }
}

$fixtureRoot = Join-Path $workRootFull 'fixture'
New-Item -ItemType Directory -Path $fixtureRoot -Force | Out-Null
$fixture = Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot -RunDir $fixtureRoot
if ([string]::IsNullOrWhiteSpace($fixture) -or -not (Test-Path -LiteralPath $fixture)) {
    throw 'Barnase fixture could not be resolved for Scope W.'
}

$sampleCache = @{}
function Get-ScopeWSampledPqr([int]$Cap) {
    if ($sampleCache.ContainsKey($Cap)) { return $sampleCache[$Cap].path }
    $sampleDir = Join-Path $workRootFull "sample_cap$Cap"
    New-Item -ItemType Directory -Path $sampleDir -Force | Out-Null
    $sampledPqr = Join-Path $sampleDir 'barnase_sampled.pqr'
    $sample = Write-FmmCliSampledPqr -InputPath $fixture -OutputPath $sampledPqr -Cap $Cap
    $sampleCache[$Cap] = [pscustomobject]@{ path=$sampledPqr; sampled_atoms=$sample.sampled_atoms }
    return $sampledPqr
}

$parserPqr = Get-ScopeWSampledPqr -Cap 64
Invoke-ScopeWParserProbe -Name 'overlap_accept' -PqrPath $parserPqr `
    -Preconditioner 'fmm-overlap-schwarz' -FmmMode 'single-matvec-analytic' `
    -OverlapCap 64 -ExpectSuccess $true
Invoke-ScopeWParserProbe -Name 'overlap_without_fmm_reject' -PqrPath $parserPqr `
    -Preconditioner 'fmm-overlap-schwarz' -FmmMode 'off' `
    -OverlapCap 64 -ExpectSuccess $false
Invoke-ScopeWParserProbe -Name 'overlap_cap129_reject' -PqrPath $parserPqr `
    -Preconditioner 'fmm-overlap-schwarz' -FmmMode 'single-matvec-analytic' `
    -OverlapCap 129 -ExpectSuccess $false

$records = New-Object System.Collections.Generic.List[object]
$directByCondition = @{}
$baselineByCondition = @{}

if ($Quick) {
    $condition = [pscustomobject]@{ name='cap128_240_r64_1e-6'; cap=128; max_iter=240; restart=64; tol=1.0e-6 }
    $pqr = Get-ScopeWSampledPqr -Cap $condition.cap
    $direct = Invoke-ScopeWCase -CaseRoot (Join-Path $workRootFull "$($condition.name)_direct") `
        -PqrPath $pqr -Kind direct -MaxIter $condition.max_iter `
        -Restart $condition.restart -Tolerance $condition.tol -BlockCap 32 -OverlapCap 0
    $directByCondition[$condition.name] = New-ScopeWRecord -Condition $condition.name `
        -AtomCap $condition.cap -Backend 'direct' -BlockCap 0 -OverlapCap 0 `
        -MaxIter $condition.max_iter -Restart $condition.restart `
        -Tolerance $condition.tol -Data $direct
    foreach ($overlap in @(0, 32, 64)) {
        $case = Invoke-ScopeWCase -CaseRoot (Join-Path $workRootFull "$($condition.name)_fmm_overlap$overlap") `
            -PqrPath $pqr -Kind fmm-overlap -MaxIter $condition.max_iter `
            -Restart $condition.restart -Tolerance $condition.tol `
            -BlockCap 32 -OverlapCap $overlap
        $record = New-ScopeWRecord -Condition $condition.name -AtomCap $condition.cap `
            -Backend 'fmm_overlap_schwarz' -BlockCap 32 -OverlapCap $overlap `
            -MaxIter $condition.max_iter -Restart $condition.restart `
            -Tolerance $condition.tol -Data $case
        $records.Add($record) | Out-Null
        if ($overlap -eq 0) { $baselineByCondition[$condition.name] = $record }
    }
} else {
    $pqr256 = Get-ScopeWSampledPqr -Cap 256
    $conditions = @(
        [pscustomobject]@{ name='cap256_480_r128_1e-6'; cap=256; max_iter=480; restart=128; tol=1.0e-6; overlaps=@(32,64,128) },
        [pscustomobject]@{ name='cap256_720_r128_1e-6'; cap=256; max_iter=720; restart=128; tol=1.0e-6; overlaps=@(32,64,128) }
    )
    foreach ($condition in $conditions) {
        $baseline = Invoke-ScopeWCase -CaseRoot (Join-Path $workRootFull "$($condition.name)_fmm_block_cap32") `
            -PqrPath $pqr256 -Kind fmm-block -MaxIter $condition.max_iter `
            -Restart $condition.restart -Tolerance $condition.tol `
            -BlockCap 32 -OverlapCap 0
        $baselineRecord = New-ScopeWRecord -Condition $condition.name -AtomCap $condition.cap `
            -Backend 'fmm_block_jacobi' -BlockCap 32 -OverlapCap 0 `
            -MaxIter $condition.max_iter -Restart $condition.restart `
            -Tolerance $condition.tol -Data $baseline
        $records.Add($baselineRecord) | Out-Null
        $baselineByCondition[$condition.name] = $baselineRecord
        foreach ($overlap in $condition.overlaps) {
            $case = Invoke-ScopeWCase -CaseRoot (Join-Path $workRootFull "$($condition.name)_fmm_overlap$overlap") `
                -PqrPath $pqr256 -Kind fmm-overlap -MaxIter $condition.max_iter `
                -Restart $condition.restart -Tolerance $condition.tol `
                -BlockCap 32 -OverlapCap $overlap
            $records.Add((New-ScopeWRecord -Condition $condition.name -AtomCap $condition.cap `
                -Backend 'fmm_overlap_schwarz' -BlockCap 32 -OverlapCap $overlap `
                -MaxIter $condition.max_iter -Restart $condition.restart `
                -Tolerance $condition.tol -Data $case)) | Out-Null
        }
    }

    $candidateRows = @($records | Where-Object { $_.backend -eq 'fmm_overlap_schwarz' })
    $best = @($candidateRows | Sort-Object `
        @{ Expression = {
                if ($_.data.result.status -eq 'ok' -and
                    (Test-ScopeWFinite $_.data.result.true_residual) -and
                    $_.data.result.true_residual -le 1.0e-6) { 0 } else { 1 }
            }; Ascending = $true },
        @{ Expression = {
                if (Test-ScopeWFinite $_.data.result.true_residual) {
                    $_.data.result.true_residual
                } else {
                    [double]::PositiveInfinity
                }
            }; Ascending = $true },
        @{ Expression = { -1 * $_.overlap_cap }; Ascending = $true })[0]

    $needsDirectVerification = $false
    foreach ($candidate in $candidateRows) {
        $baseline = if ($baselineByCondition.ContainsKey($candidate.condition)) {
            $baselineByCondition[$candidate.condition]
        } else {
            $null
        }
        if ($candidate.data.result.status -eq 'ok' -and
            (Test-ScopeWFinite $candidate.data.result.true_residual) -and
            $candidate.data.result.true_residual -le 1.0e-6) {
            $needsDirectVerification = $true
        }
        if ($null -ne $baseline -and
            (Test-ScopeWFinite $baseline.data.result.true_residual) -and
            (Test-ScopeWFinite $candidate.data.result.true_residual) -and
            $candidate.data.result.true_residual -gt 0.0) {
            $improvement = $baseline.data.result.true_residual /
                $candidate.data.result.true_residual
            if ($improvement -ge 4.0) {
                $needsDirectVerification = $true
            }
        }
    }

    if ($needsDirectVerification) {
        foreach ($directCondition in @(
                [pscustomobject]@{ name='cap256_480_r128_1e-6'; cap=256; max_iter=480; restart=128; tol=1.0e-6 },
                [pscustomobject]@{ name=$best.condition; cap=$best.atom_cap; max_iter=$best.max_iter; restart=$best.restart; tol=$best.tolerance }
            )) {
            if (-not $directByCondition.ContainsKey($directCondition.name)) {
                $direct = Invoke-ScopeWCase -CaseRoot (Join-Path $workRootFull "$($directCondition.name)_direct") `
                    -PqrPath $pqr256 -Kind direct -MaxIter $directCondition.max_iter `
                    -Restart $directCondition.restart -Tolerance $directCondition.tol `
                    -BlockCap 32 -OverlapCap 0
                $directByCondition[$directCondition.name] = New-ScopeWRecord -Condition $directCondition.name `
                    -AtomCap $directCondition.cap -Backend 'direct' -BlockCap 0 -OverlapCap 0 `
                    -MaxIter $directCondition.max_iter -Restart $directCondition.restart `
                    -Tolerance $directCondition.tol -Data $direct
            }
        }
    }
}

$recordsForOutput = New-Object System.Collections.Generic.List[object]
foreach ($directRecord in $directByCondition.Values) { $recordsForOutput.Add($directRecord) | Out-Null }
foreach ($record in $records) { $recordsForOutput.Add($record) | Out-Null }

$rows = New-Object System.Collections.Generic.List[object]
foreach ($record in $recordsForOutput) {
    $direct = if ($directByCondition.ContainsKey($record.condition)) { $directByCondition[$record.condition] } else { $null }
    $baseline = if ($baselineByCondition.ContainsKey($record.condition)) { $baselineByCondition[$record.condition] } else { $null }
    $class = switch ($record.backend) {
        'direct' { 'direct_reference' }
        'fmm_block_jacobi' { 'block_jacobi_reference' }
        default { Get-ScopeWClass -Candidate $record -Baseline $baseline -Direct $direct }
    }
    $speedup = Get-ScopeWSpeedup -Direct $direct -Candidate $record
    $energyDelta = Get-ScopeWEnergyDelta -Direct $direct -Candidate $record
    $residualImprovement = [double]::NaN
    $memoryRatio = [double]::NaN
    if ($null -ne $baseline -and $record.backend -eq 'fmm_overlap_schwarz' -and
        (Test-ScopeWFinite $baseline.data.result.true_residual) -and
        (Test-ScopeWFinite $record.data.result.true_residual) -and
        $record.data.result.true_residual -gt 0.0) {
        $residualImprovement = $baseline.data.result.true_residual / $record.data.result.true_residual
    }
    if ($null -ne $baseline -and $record.backend -eq 'fmm_overlap_schwarz' -and
        $baseline.data.run.peak_working_set_mb -gt 0.0) {
        $memoryRatio = $record.data.run.peak_working_set_mb / $baseline.data.run.peak_working_set_mb
    }
    $rows.Add([pscustomobject]@{
        condition=$record.condition
        atom_cap=$record.atom_cap
        backend=$record.backend
        block_cap=$record.block_cap
        overlap_cap=$record.overlap_cap
        max_iter=$record.max_iter
        restart=$record.restart
        tolerance=$record.tolerance
        exit_code=$record.data.run.exit_code
        status=$record.data.result.status
        panels=$record.data.result.panels
        iterations=$record.data.result.iterations
        true_residual=$record.data.result.true_residual
        energy_kj_mol=$record.data.result.energy_kj_mol
        energy_delta_abs=$energyDelta
        audit_rel=$record.data.result.tree_vs_direct_matvec_rel
        matvec_seconds_avg=$record.data.result.matvec_seconds_avg
        speedup_vs_direct=$speedup
        residual_improvement_vs_baseline=$residualImprovement
        memory_ratio_vs_baseline=$memoryRatio
        wall_seconds=$record.data.run.wall_seconds
        peak_working_set_mb=$record.data.run.peak_working_set_mb
        fmm_m2l_mode=$record.data.result.fmm_m2l_mode
        fmm_far_pairs=$record.data.result.fmm_far_pairs
        preconditioner_requested=$record.data.preconditioner.requested
        preconditioner_effective=$record.data.preconditioner.effective
        preconditioner_block_panel_cap=$record.data.preconditioner.block_cap
        preconditioner_overlap_panel_cap=$record.data.preconditioner.overlap_cap
        preconditioner_estimated_bytes=$record.data.preconditioner.estimated_bytes
        preconditioner_valid_blocks=$record.data.preconditioner.valid_blocks
        preconditioner_valid_panels=$record.data.preconditioner.valid_panels
        preconditioner_quality_status=$record.data.preconditioner.quality_status
        preconditioner_setup_seconds=$record.data.preconditioner.setup_seconds
        preconditioner_apply_seconds=$record.data.preconditioner.apply_seconds
        preconditioner_apply_calls=$record.data.preconditioner.apply_calls
        scope_w_class=$class
    }) | Out-Null
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation

$overlapRows = @($rows | Where-Object { $_.backend -eq 'fmm_overlap_schwarz' })
$strictVerified = @($overlapRows | Where-Object { $_.scope_w_class -eq 'strict_converged_verified' })
$strictPending = @($overlapRows | Where-Object { $_.scope_w_class -eq 'strict_fmm_pending_direct' })
$improved4x = @($overlapRows | Where-Object { $_.scope_w_class -eq 'residual_improved_4x' })
$improved = @($overlapRows | Where-Object { $_.scope_w_class -eq 'residual_improved' })
$decision = if ($strictVerified.Count -gt 0) {
    'strict_converged_verified'
} elseif ($strictPending.Count -gt 0) {
    'strict_fmm_pending_direct'
} elseif ($improved4x.Count -gt 0) {
    'residual_improved_4x'
} elseif ($improved.Count -gt 0) {
    'residual_improved'
} else {
    'no_win'
}

[pscustomobject]@{
    fixture='barnase'
    mode=if ($Quick) { 'quick' } else { 'full' }
    rows=$rows.Count
    overlap_rows=$overlapRows.Count
    strict_verified_rows=$strictVerified.Count
    strict_pending_rows=$strictPending.Count
    improved4x_rows=$improved4x.Count
    improved_rows=$improved.Count
    scope_w_decision=$decision
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'Scope W FMM Restricted Additive-Schwarz Probe',
    "mode: $(if ($Quick) { 'quick' } else { 'full' })",
    "rows: $($rows.Count)",
    "decision: $decision",
    'production_default: grid_remains_default',
    'fmm_m2l_status: hidden_experimental',
    'parser_probes: overlap_accept, overlap_without_fmm_reject, overlap_cap129_reject'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = New-Object System.Collections.Generic.List[string]
$report.Add('# TITAN FMM/M2L Scope W Overlap-Schwarz Report')
$report.Add('')
$report.Add("mode: $(if ($Quick) { 'quick' } else { 'full' })")
$report.Add("decision: $decision")
$report.Add('')
$report.Add('Rows:')
foreach ($row in $rows) {
    $report.Add(('{0} {1} owned={2} overlap={3} status={4} residual={5} audit={6} speedup={7} mem_ratio={8} precond={9} class={10}' -f
        $row.condition, $row.backend, $row.block_cap, $row.overlap_cap,
        $row.status, $row.true_residual, $row.audit_rel,
        $row.speedup_vs_direct, $row.memory_ratio_vs_baseline,
        $row.preconditioner_effective, $row.scope_w_class))
}
$report.Add('')
$report.Add('Notes:')
$report.Add('- FMM/M2L and fmm-overlap-schwarz remain hidden and experimental.')
$report.Add('- Restricted Additive-Schwarz factors owned plus halo panels but writes back owned panels only.')
$report.Add('- Halo panels are selected from FMM near-source leaf work and sorted deterministically by owned-block centroid distance and panel id.')
$report.Add('- If overlap gives no useful improvement, the next proposed direction is sparse near-field ILU rather than more leaf-local cap tuning.')
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "Scope W report: $ReportPath"
Write-Host "Scope W decision: $decision"
