param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [switch]$Resume,
    [int]$CaseTimeoutSeconds = 2400,
    [int]$Threads = 4,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_x')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fmm_m2l_scope_x_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_x_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_x_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_X_SPARSE_ILU_REPORT.txt'

$SurfacePoints = 8
$OpeningTheta = 0.6
$TemperatureK = 298.15
$AuditSampleCount = if ($Quick) { 128 } else { 512 }
$AtomCap = if ($Quick) { 128 } else { 256 }

if ($Quick) { $CaseTimeoutSeconds = [Math]::Min($CaseTimeoutSeconds, 600) }
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope X work root outside test/tmp: $workRootFull"
}
if ((Test-Path -LiteralPath $workRootFull) -and -not $Resume) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function New-ScopeXEmptyResult([string]$Status) {
    return [pscustomobject]@{
        status=$Status; panels=0; iterations=0
        true_residual=[double]::NaN; energy_kj_mol=[double]::NaN
        matvec_seconds_avg=[double]::NaN; tree_vs_direct_matvec_rel=[double]::NaN
        fmm_m2l_mode=''; fmm_far_pairs=0; fmm_finite_outputs='false'
    }
}

function Get-ScopeXRecoveredRunResult([string]$RunDir) {
    $progressPath = Join-Path $RunDir 'pbe_tabi_longrun_progress.csv'
    $wallSeconds = [double]::NaN
    $peakBytes = 0L
    if (Test-Path -LiteralPath $progressPath) {
        $progress = @(Import-Csv -LiteralPath $progressPath)
        foreach ($row in $progress) {
            $elapsed = Convert-FmmCliDouble $row.elapsed_seconds
            $rss = Convert-FmmCliInt64 $row.rss_bytes
            if (Test-FmmCliFinite $elapsed) {
                if (-not (Test-FmmCliFinite $wallSeconds) -or $elapsed -gt $wallSeconds) {
                    $wallSeconds = $elapsed
                }
            }
            if ($rss -gt $peakBytes) { $peakBytes = $rss }
        }
    }
    return [pscustomobject]@{
        exit_code=0; timed_out=$false; wall_seconds=$wallSeconds
        peak_working_set_bytes=$peakBytes
        peak_working_set_mb=$peakBytes / 1MB
    }
}

function Read-ScopeXPrecondProfile([string]$RunDir) {
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $profilePath)) {
        return [pscustomobject]@{
            requested='NA'; effective='NA'; ilu_neighbor_cap=0; ilu_drop_tol=[double]::NaN
            estimated_bytes=[double]::NaN; valid_blocks=0; valid_panels=0
            quality_status='NA'; setup_seconds=[double]::NaN; apply_seconds=[double]::NaN
            apply_calls=0; sparse_nnz_blocks=0; sparse_avg_degree=[double]::NaN
            sparse_max_degree=0
        }
    }
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    return [pscustomobject]@{
        requested=[string]$profile.preconditioner_requested
        effective=[string]$profile.preconditioner_effective
        ilu_neighbor_cap=Convert-FmmCliInt64 $profile.preconditioner_ilu_neighbor_cap
        ilu_drop_tol=Convert-FmmCliDouble $profile.preconditioner_ilu_drop_tol
        estimated_bytes=Convert-FmmCliDouble $profile.preconditioner_estimated_bytes
        valid_blocks=Convert-FmmCliInt64 $profile.preconditioner_valid_blocks
        valid_panels=Convert-FmmCliInt64 $profile.preconditioner_valid_panels
        quality_status=[string]$profile.preconditioner_quality_status
        setup_seconds=Convert-FmmCliDouble $profile.preconditioner_setup_seconds
        apply_seconds=Convert-FmmCliDouble $profile.preconditioner_apply_seconds
        apply_calls=Convert-FmmCliInt64 $profile.preconditioner_apply_calls
        sparse_nnz_blocks=Convert-FmmCliInt64 $profile.preconditioner_sparse_nnz_blocks
        sparse_avg_degree=Convert-FmmCliDouble $profile.preconditioner_sparse_avg_degree
        sparse_max_degree=Convert-FmmCliInt64 $profile.preconditioner_sparse_max_degree
    }
}

function Invoke-ScopeXCase(
    [string]$CaseRoot,
    [string]$PqrPath,
    [ValidateSet('direct','fmm-diagonal','fmm-block','fmm-ilu')][string]$Kind,
    [int]$MaxIter,
    [int]$Restart,
    [double]$Tolerance,
    [int]$IluNeighborCap
) {
    $profilePath = Join-Path $CaseRoot 'pbe_tabi_profile.csv'
    $energyPath = Join-Path $CaseRoot 'pbe_energy.csv'
    $timingPath = Join-Path $CaseRoot 'pbe_timing.csv'
    if ($Resume -and (Test-Path -LiteralPath $profilePath) -and
        (Test-Path -LiteralPath $energyPath) -and
        (Test-Path -LiteralPath $timingPath)) {
        $expected = if ($Kind -eq 'direct') { 'off' } else { 'single-matvec-analytic' }
        return [pscustomobject]@{
            run=Get-ScopeXRecoveredRunResult -RunDir $CaseRoot
            result=Read-FmmCliCaseResult -RunDir $CaseRoot -ExpectedFmmMode $expected `
                -RequireAudit ($Kind -ne 'direct')
            preconditioner=Read-ScopeXPrecondProfile -RunDir $CaseRoot
        }
    }
    New-FmmCliCleanDir -Path $CaseRoot -WorkRoot $WorkRoot
    $mode = if ($Kind -eq 'direct') { 'direct' } else { 'treecode' }
    $fmmMode = if ($Kind -eq 'direct') { 'off' } else { 'single-matvec-analytic' }
    $precond = switch ($Kind) {
        'fmm-block' { 'fmm-block-jacobi' }
        'fmm-ilu' { 'fmm-nearfield-ilu0' }
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
        -PreconditionerBlockPanelCap 32 `
        -PreconditionerIluNeighborCap $IluNeighborCap `
        -PreconditionerIluDropTol 0.0
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $CaseRoot `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $result = if ($run.exit_code -eq 0) {
        $expected = if ($Kind -eq 'direct') { 'off' } else { 'single-matvec-analytic' }
        Read-FmmCliCaseResult -RunDir $CaseRoot -ExpectedFmmMode $expected `
            -RequireAudit ($Kind -ne 'direct')
    } else {
        New-ScopeXEmptyResult "exit_$($run.exit_code)"
    }
    $precondProfile = Read-ScopeXPrecondProfile -RunDir $CaseRoot
    return [pscustomobject]@{
        run=$run; result=$result; preconditioner=$precondProfile
    }
}

function Get-ScopeXResidualImprovement([object]$Candidate, [object]$Baseline) {
    if ($null -eq $Baseline -or
        -not (Test-FmmCliFinite $Baseline.data.result.true_residual) -or
        -not (Test-FmmCliFinite $Candidate.data.result.true_residual) -or
        $Candidate.data.result.true_residual -le 0.0) {
        return [double]::NaN
    }
    return $Baseline.data.result.true_residual /
        $Candidate.data.result.true_residual
}

function Test-ScopeXIluCandidate([object]$Record) {
    return $Record.kind -eq 'fmm-ilu' -and
        $Record.data.run.exit_code -eq 0 -and
        $Record.data.preconditioner.effective -eq 'fmm_nearfield_ilu0' -and
        $Record.data.preconditioner.quality_status -eq 'pass' -and
        $Record.data.preconditioner.sparse_nnz_blocks -gt 0 -and
        $Record.data.result.fmm_m2l_mode -eq 'single-matvec-analytic' -and
        $Record.data.result.fmm_far_pairs -gt 0 -and
        $Record.data.result.fmm_finite_outputs -eq 'true' -and
        (Test-FmmCliFinite $Record.data.result.true_residual) -and
        (Test-FmmCliFinite $Record.data.result.tree_vs_direct_matvec_rel) -and
        $Record.data.result.tree_vs_direct_matvec_rel -le 1.0e-2
}

function Select-ScopeXLongCandidates([object[]]$Records) {
    return @($Records | Where-Object { Test-ScopeXIluCandidate $_ } |
        Sort-Object `
            @{ Expression = { [double]$_.data.result.true_residual }; Ascending = $true }, `
            @{ Expression = { [int]$_.ilu_neighbor_cap }; Ascending = $true } |
        Select-Object -First 2)
}

function Invoke-ScopeXDirectVerification(
    [object]$Candidate,
    [string]$PqrPath,
    [string]$Root
) {
    $name = 'direct_verify_{0}_r{1}' -f $Candidate.max_iter, $Candidate.restart
    $data = Invoke-ScopeXCase -CaseRoot (Join-Path $Root $name) `
        -PqrPath $PqrPath -Kind direct -MaxIter $Candidate.max_iter `
        -Restart $Candidate.restart -Tolerance $Candidate.tolerance `
        -IluNeighborCap 192
    return [pscustomobject]@{
        condition=$name; stage='direct-verification'; kind='direct'
        max_iter=$Candidate.max_iter; restart=$Candidate.restart
        tolerance=$Candidate.tolerance; ilu_neighbor_cap=0; data=$data
    }
}

function Get-ScopeXMetrics(
    [object]$Record,
    [object]$Baseline,
    [object]$Direct
) {
    $improvement = Get-ScopeXResidualImprovement -Candidate $Record -Baseline $Baseline
    $memoryRatio = [double]::NaN
    if ($null -ne $Baseline -and $Baseline.data.run.peak_working_set_mb -gt 0.0) {
        $memoryRatio = $Record.data.run.peak_working_set_mb /
            $Baseline.data.run.peak_working_set_mb
    }

    $directStatus = if ($Record.kind -eq 'direct') { 'self' } elseif ($null -eq $Direct) { 'not_run' } `
        elseif ($Direct.data.run.exit_code -eq 124) { 'timeout_not_verified' } `
        elseif ($Direct.data.run.exit_code -ne 0) { "exit_$($Direct.data.run.exit_code)" } `
        elseif ($Direct.data.result.status -ne 'ok') { 'completed_nonconverged' } `
        else { 'completed' }
    $energyDelta = [double]::NaN
    $speedup = [double]::NaN
    $residualRatio = [double]::NaN
    if ($null -ne $Direct) {
        if ((Test-FmmCliFinite $Direct.data.result.energy_kj_mol) -and
            (Test-FmmCliFinite $Record.data.result.energy_kj_mol)) {
            $energyDelta = [Math]::Abs($Record.data.result.energy_kj_mol -
                $Direct.data.result.energy_kj_mol)
        }
        if ((Test-FmmCliFinite $Direct.data.result.matvec_seconds_avg) -and
            (Test-FmmCliFinite $Record.data.result.matvec_seconds_avg) -and
            $Record.data.result.matvec_seconds_avg -gt 0.0) {
            $speedup = $Direct.data.result.matvec_seconds_avg /
                $Record.data.result.matvec_seconds_avg
        }
        if ((Test-FmmCliFinite $Direct.data.result.true_residual) -and
            $Direct.data.result.true_residual -gt 0.0 -and
            (Test-FmmCliFinite $Record.data.result.true_residual)) {
            $residualRatio = $Record.data.result.true_residual /
                $Direct.data.result.true_residual
        }
    }
    return [pscustomobject]@{
        residual_improvement_vs_baseline=$improvement
        memory_ratio_vs_baseline=$memoryRatio
        direct_verification_status=$directStatus
        energy_delta_abs=$energyDelta
        speedup_vs_direct=$speedup
        residual_ratio_vs_direct=$residualRatio
    }
}

function Get-ScopeXClass(
    [object]$Record,
    [object]$Baseline,
    [object]$Direct,
    [object]$Metrics
) {
    if ($Record.kind -eq 'direct') { return 'direct_reference' }
    if ($Record.kind -eq 'fmm-diagonal') { return 'diagonal_reference' }
    if ($Record.kind -eq 'fmm-block') { return 'block_reference' }
    if ($Record.data.run.exit_code -eq 124) { return 'timeout_not_verified' }
    if ($Record.data.run.exit_code -ne 0) { return 'factor_limited' }
    if ($Record.data.preconditioner.effective -ne 'fmm_nearfield_ilu0' -or
        $Record.data.preconditioner.quality_status -ne 'pass' -or
        $Record.data.preconditioner.sparse_nnz_blocks -le 0) {
        return 'factor_limited'
    }
    if ($Record.data.preconditioner.ilu_neighbor_cap -ne $Record.ilu_neighbor_cap -or
        $Record.data.result.fmm_m2l_mode -ne 'single-matvec-analytic' -or
        $Record.data.result.fmm_far_pairs -le 0 -or
        $Record.data.result.fmm_finite_outputs -ne 'true') {
        return 'factor_limited'
    }
    if (-not (Test-FmmCliFinite $Record.data.result.tree_vs_direct_matvec_rel) -or
        $Record.data.result.tree_vs_direct_matvec_rel -gt 1.0e-2) {
        return 'accuracy_limited'
    }
    if ($Record.data.run.peak_working_set_mb -ge 256.0 -or
        ((Test-FmmCliFinite $Metrics.memory_ratio_vs_baseline) -and
            $Metrics.memory_ratio_vs_baseline -ge 4.0)) {
        return 'memory_limited'
    }

    $strict = $Record.data.result.status -eq 'ok' -and
        (Test-FmmCliFinite $Record.data.result.true_residual) -and
        $Record.data.result.true_residual -le 1.0e-6
    $useful = (Test-FmmCliFinite $Metrics.residual_improvement_vs_baseline) -and
        $Metrics.residual_improvement_vs_baseline -ge 4.0
    if ($strict) {
        if ($Metrics.direct_verification_status -ne 'completed') {
            return 'strict_fmm_pending_direct'
        }
        if ($Metrics.direct_verification_status -eq 'completed') {
            $energyGate = [Math]::Max(5.0,
                0.05 * [Math]::Abs($Direct.data.result.energy_kj_mol))
            if ((Test-FmmCliFinite $Metrics.energy_delta_abs) -and
                $Metrics.energy_delta_abs -le $energyGate -and
                (Test-FmmCliFinite $Metrics.speedup_vs_direct) -and
                $Metrics.speedup_vs_direct -ge 1.2) {
                return 'strict_verified'
            }
        }
        return 'verification_failed'
    }
    if ($useful) { return 'useful_not_strict' }
    return 'no_useful_improvement'
}

$fixtureRoot = Join-Path $workRootFull 'fixtures'
New-Item -ItemType Directory -Path $fixtureRoot -Force | Out-Null
$barnase = Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot -RunDir $fixtureRoot
if ([string]::IsNullOrWhiteSpace($barnase) -or -not (Test-Path -LiteralPath $barnase)) {
    $skip = [pscustomobject]@{ scope_x_decision='skipped_missing_fixture'; reason='barnase fixture not resolved' }
    $skip | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation
    Set-Content -LiteralPath $SummaryPath -Value @('scope_x_decision=skipped_missing_fixture') -Encoding ASCII
    Set-Content -LiteralPath $ReportPath -Value @('Scope X skipped: Barnase fixture not resolved.') -Encoding ASCII
    return
}
$sampledPqr = Join-Path $fixtureRoot "barnase_cap$AtomCap.pqr"
$sampleInfo = Write-FmmCliSampledPqr -InputPath $barnase -OutputPath $sampledPqr -Cap $AtomCap

$records = New-Object System.Collections.Generic.List[object]
$initialCases = if ($Quick) {
    @(
        @{ name='direct_120_r64'; stage='quick'; kind='direct'; max_iter=120; restart=64; tol=1.0e-6; cap=192 },
        @{ name='diag_120_r64'; stage='quick'; kind='fmm-diagonal'; max_iter=120; restart=64; tol=1.0e-6; cap=192 },
        @{ name='block_120_r64'; stage='quick'; kind='fmm-block'; max_iter=120; restart=64; tol=1.0e-6; cap=192 },
        @{ name='ilu64_120_r64'; stage='quick'; kind='fmm-ilu'; max_iter=120; restart=64; tol=1.0e-6; cap=64 },
        @{ name='ilu128_120_r64'; stage='quick'; kind='fmm-ilu'; max_iter=120; restart=64; tol=1.0e-6; cap=128 },
        @{ name='ilu192_120_r64'; stage='quick'; kind='fmm-ilu'; max_iter=120; restart=64; tol=1.0e-6; cap=192 }
    )
} else {
    @(
        @{ name='block_720_r128'; stage='baseline'; kind='fmm-block'; max_iter=720; restart=128; tol=1.0e-6; cap=192 },
        @{ name='ilu64_480_r128'; stage='screen'; kind='fmm-ilu'; max_iter=480; restart=128; tol=1.0e-6; cap=64 },
        @{ name='ilu128_480_r128'; stage='screen'; kind='fmm-ilu'; max_iter=480; restart=128; tol=1.0e-6; cap=128 },
        @{ name='ilu192_480_r128'; stage='screen'; kind='fmm-ilu'; max_iter=480; restart=128; tol=1.0e-6; cap=192 },
        @{ name='ilu256_480_r128'; stage='screen'; kind='fmm-ilu'; max_iter=480; restart=128; tol=1.0e-6; cap=256 }
    )
}

foreach ($case in $initialCases) {
    $runDir = Join-Path $workRootFull $case.name
    $data = Invoke-ScopeXCase -CaseRoot $runDir -PqrPath $sampledPqr `
        -Kind $case.kind -MaxIter $case.max_iter -Restart $case.restart `
        -Tolerance $case.tol -IluNeighborCap $case.cap
    $records.Add([pscustomobject]@{
        condition=$case.name; stage=$case.stage; kind=$case.kind
        max_iter=$case.max_iter; restart=$case.restart; tolerance=$case.tol
        ilu_neighbor_cap=$case.cap; data=$data
    }) | Out-Null
}

$baseline = @($records | Where-Object { $_.kind -eq 'fmm-block' } | Select-Object -First 1)
$baseline = if ($baseline.Count -gt 0) { $baseline[0] } else { $null }

if (-not $Quick) {
    $screenRows = @($records | Where-Object { $_.stage -eq 'screen' })
    $longCandidates = @(Select-ScopeXLongCandidates -Records $screenRows)
    foreach ($candidate in $longCandidates) {
        $name = 'ilu{0}_720_r128' -f $candidate.ilu_neighbor_cap
        $data = Invoke-ScopeXCase -CaseRoot (Join-Path $workRootFull $name) `
            -PqrPath $sampledPqr -Kind fmm-ilu -MaxIter 720 -Restart 128 `
            -Tolerance 1.0e-6 -IluNeighborCap $candidate.ilu_neighbor_cap
        $records.Add([pscustomobject]@{
            condition=$name; stage='long'; kind='fmm-ilu'
            max_iter=720; restart=128; tolerance=1.0e-6
            ilu_neighbor_cap=$candidate.ilu_neighbor_cap; data=$data
        }) | Out-Null
    }
}

$directByBudget = @{}
$quickDirect = @($records | Where-Object { $_.kind -eq 'direct' } | Select-Object -First 1)
if ($quickDirect.Count -gt 0) {
    $directByBudget['120|64|1E-06'] = $quickDirect[0]
}

if (-not $Quick) {
    $longVerificationPool = @($records | Where-Object { $_.stage -eq 'long' })
    $verificationPool = if ($longVerificationPool.Count -gt 0) {
        $longVerificationPool
    } else {
        @($records | Where-Object { $_.stage -eq 'screen' })
    }
    $verificationCandidates = @($verificationPool | Where-Object {
        if (-not (Test-ScopeXIluCandidate $_)) { return $false }
        $improvement = Get-ScopeXResidualImprovement -Candidate $_ -Baseline $baseline
        $strict = $_.data.result.status -eq 'ok' -and
            $_.data.result.true_residual -le 1.0e-6
        return $strict -or ((Test-FmmCliFinite $improvement) -and $improvement -ge 4.0)
    } | Sort-Object `
        @{ Expression = {
                if ($_.data.result.status -eq 'ok' -and
                    $_.data.result.true_residual -le 1.0e-6) { 0 } else { 1 }
            }; Ascending = $true }, `
        @{ Expression = { [int]$_.ilu_neighbor_cap }; Ascending = $true }, `
        @{ Expression = { [double]$_.data.result.true_residual }; Ascending = $true })
    if ($verificationCandidates.Count -gt 0) {
        $candidate = $verificationCandidates[0]
        $direct = Invoke-ScopeXDirectVerification -Candidate $candidate `
            -PqrPath $sampledPqr -Root $workRootFull
        $key = '{0}|{1}|{2}' -f $candidate.max_iter, $candidate.restart,
            ([string]$candidate.tolerance)
        $directByBudget[$key] = $direct
        $records.Add($direct) | Out-Null
    }
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($record in $records) {
    $key = '{0}|{1}|{2}' -f $record.max_iter, $record.restart,
        ([string]$record.tolerance)
    $direct = if ($record.kind -eq 'direct') { $record } `
        elseif ($directByBudget.ContainsKey($key)) { $directByBudget[$key] } `
        else { $null }
    $metrics = Get-ScopeXMetrics -Record $record -Baseline $baseline -Direct $direct
    $class = Get-ScopeXClass -Record $record -Baseline $baseline `
        -Direct $direct -Metrics $metrics
    $rows.Add([pscustomobject]@{
        condition=$record.condition; stage=$record.stage; fixture='barnase'
        atom_cap=$AtomCap; sampled_atoms=$sampleInfo.sampled_atoms
        kind=$record.kind; max_iter=$record.max_iter; restart=$record.restart
        tolerance=$record.tolerance; exit_code=$record.data.run.exit_code
        timed_out=$record.data.run.timed_out; wall_seconds=$record.data.run.wall_seconds
        peak_working_set_mb=$record.data.run.peak_working_set_mb
        status=$record.data.result.status; iterations=$record.data.result.iterations
        true_residual=$record.data.result.true_residual
        energy_kj_mol=$record.data.result.energy_kj_mol
        matvec_seconds_avg=$record.data.result.matvec_seconds_avg
        audit_rel=$record.data.result.tree_vs_direct_matvec_rel
        fmm_m2l_mode=$record.data.result.fmm_m2l_mode
        fmm_far_pairs=$record.data.result.fmm_far_pairs
        preconditioner_requested=$record.data.preconditioner.requested
        preconditioner_effective=$record.data.preconditioner.effective
        ilu_neighbor_cap=$record.data.preconditioner.ilu_neighbor_cap
        ilu_drop_tol=$record.data.preconditioner.ilu_drop_tol
        estimated_bytes=$record.data.preconditioner.estimated_bytes
        valid_blocks=$record.data.preconditioner.valid_blocks
        valid_panels=$record.data.preconditioner.valid_panels
        preconditioner_setup_seconds=$record.data.preconditioner.setup_seconds
        preconditioner_apply_seconds=$record.data.preconditioner.apply_seconds
        preconditioner_apply_calls=$record.data.preconditioner.apply_calls
        sparse_nnz_blocks=$record.data.preconditioner.sparse_nnz_blocks
        sparse_avg_degree=$record.data.preconditioner.sparse_avg_degree
        sparse_max_degree=$record.data.preconditioner.sparse_max_degree
        quality_status=$record.data.preconditioner.quality_status
        residual_improvement_vs_baseline=$metrics.residual_improvement_vs_baseline
        memory_ratio_vs_baseline=$metrics.memory_ratio_vs_baseline
        direct_verification_status=$metrics.direct_verification_status
        direct_true_residual=if ($null -ne $direct) { $direct.data.result.true_residual } else { [double]::NaN }
        direct_energy_kj_mol=if ($null -ne $direct) { $direct.data.result.energy_kj_mol } else { [double]::NaN }
        energy_delta_abs=$metrics.energy_delta_abs
        speedup_vs_direct=$metrics.speedup_vs_direct
        residual_ratio_vs_direct=$metrics.residual_ratio_vs_direct
        scope_x_class=$class
    }) | Out-Null
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$candidateRows = @($rows | Where-Object { $_.kind -eq 'fmm-ilu' })
$best = @($candidateRows | Sort-Object `
    @{ Expression = {
            switch ($_.scope_x_class) {
                'strict_verified' { 0 }
                'strict_fmm_pending_direct' { 1 }
                'useful_not_strict' { 2 }
                'memory_limited' { 3 }
                'factor_limited' { 4 }
                default { 5 }
            }
        }; Ascending = $true }, `
    @{ Expression = { [int]$_.ilu_neighbor_cap }; Ascending = $true }, `
    @{ Expression = {
            if (Test-FmmCliFinite $_.true_residual) { [double]$_.true_residual }
            else { [double]::PositiveInfinity }
        }; Ascending = $true } |
    Select-Object -First 1)
$best = if ($best.Count -gt 0) { $best[0] } else { $null }
$finalClass = if (@($candidateRows | Where-Object { $_.scope_x_class -eq 'strict_verified' }).Count -gt 0) {
    'strict_verified'
} elseif (@($candidateRows | Where-Object { $_.scope_x_class -eq 'strict_fmm_pending_direct' }).Count -gt 0) {
    'strict_fmm_pending_direct'
} elseif (@($candidateRows | Where-Object { $_.scope_x_class -eq 'useful_not_strict' }).Count -gt 0) {
    'useful_not_strict'
} elseif (@($candidateRows | Where-Object { $_.scope_x_class -eq 'memory_limited' }).Count -gt 0) {
    'memory_limited'
} elseif (@($candidateRows | Where-Object { $_.scope_x_class -eq 'factor_limited' }).Count -gt 0) {
    'factor_limited'
} elseif (@($candidateRows | Where-Object { $_.scope_x_class -eq 'timeout_not_verified' }).Count -gt 0) {
    'timeout_not_verified'
} else {
    'no_useful_improvement'
}
[pscustomobject]@{
    scope_x_decision=$finalClass
    rows=$rows.Count
    atom_cap=$AtomCap
    best_condition=if ($null -ne $best) { $best.condition } else { '' }
    best_neighbor_cap=if ($null -ne $best) { $best.ilu_neighbor_cap } else { 0 }
    direct_verification_status=if ($null -ne $best) { $best.direct_verification_status } else { 'not_run' }
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

Set-Content -LiteralPath $SummaryPath -Encoding ASCII -Value @(
    "scope_x_decision=$finalClass",
    "fixture=barnase",
    "atom_cap=$AtomCap",
    "sampled_atoms=$($sampleInfo.sampled_atoms)",
    "best_condition=$(if ($null -ne $best) { $best.condition } else { '' })",
    "direct_verification_status=$(if ($null -ne $best) { $best.direct_verification_status } else { 'not_run' })",
    "matrix=$MatrixCsv",
    "note=FMM/M2L remains hidden and grid remains production default."
)

$lines = New-Object System.Collections.Generic.List[string]
$lines.Add('TITAN FMM/M2L Scope X Sparse Near-Field ILU Report')
$lines.Add("Decision: $finalClass")
$lines.Add("Fixture: barnase sampled cap $AtomCap ($($sampleInfo.sampled_atoms) atoms)")
$lines.Add('')
$lines.Add('Rows:')
foreach ($row in $rows) {
    $lines.Add(('{0}: stage={1} kind={2} effective={3} residual={4:G6} audit={5:G6} nnz={6} improve={7:G5} mem_ratio={8:G5} direct={9} dE={10:G5} speedup={11:G5} class={12}' -f
        $row.condition, $row.stage, $row.kind, $row.preconditioner_effective,
        [double]$row.true_residual, [double]$row.audit_rel,
        [int64]$row.sparse_nnz_blocks,
        [double]$row.residual_improvement_vs_baseline,
        [double]$row.memory_ratio_vs_baseline,
        $row.direct_verification_status, [double]$row.energy_delta_abs,
        [double]$row.speedup_vs_direct, $row.scope_x_class))
}
$lines.Add('')
$lines.Add('Production default unchanged: -pbe=grid. FMM/M2L sparse ILU is hidden/experimental.')
Set-Content -LiteralPath $ReportPath -Value $lines -Encoding ASCII

Write-Host "Scope X report: $ReportPath"
