param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 4,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_u')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fmm_m2l_scope_u_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_u_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_u_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_U_FMM_BLOCK_PRECONDITIONER_REPORT.txt'

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
    throw "Refusing to clean Scope U work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Test-ScopeUFinite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and
        (-not [double]::IsInfinity($Value))
}

function New-ScopeUEmptyResult([string]$Status) {
    return [pscustomobject]@{
        status=$Status; panels=0; iterations=0
        true_residual=[double]::NaN; energy_kj_mol=[double]::NaN
        matvec_seconds_avg=[double]::NaN; tree_vs_direct_matvec_rel=[double]::NaN
        fmm_m2l_mode=''; fmm_far_pairs=0; fmm_finite_outputs='false'
    }
}

function Read-ScopeUPrecondProfile([string]$RunDir) {
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $profilePath)) {
        return [pscustomobject]@{
            requested='NA'; effective='NA'; estimated_bytes=[double]::NaN
            valid_blocks=0; valid_panels=0; quality_status='NA'
            setup_seconds=[double]::NaN; apply_seconds=[double]::NaN
            apply_calls=0
        }
    }
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    return [pscustomobject]@{
        requested=[string]$profile.preconditioner_requested
        effective=[string]$profile.preconditioner_effective
        estimated_bytes=Convert-FmmCliDouble $profile.preconditioner_estimated_bytes
        valid_blocks=Convert-FmmCliInt64 $profile.preconditioner_valid_blocks
        valid_panels=Convert-FmmCliInt64 $profile.preconditioner_valid_panels
        quality_status=[string]$profile.preconditioner_quality_status
        setup_seconds=Convert-FmmCliDouble $profile.preconditioner_setup_seconds
        apply_seconds=Convert-FmmCliDouble $profile.preconditioner_apply_seconds
        apply_calls=Convert-FmmCliInt64 $profile.preconditioner_apply_calls
    }
}

function Invoke-ScopeUCase(
    [string]$CaseRoot,
    [string]$PqrPath,
    [ValidateSet('direct','fmm-diagonal','fmm-block')][string]$Kind,
    [int]$MaxIter,
    [int]$Restart,
    [double]$Tolerance
) {
    New-FmmCliCleanDir -Path $CaseRoot -WorkRoot $WorkRoot
    $mode = if ($Kind -eq 'direct') { 'direct' } else { 'treecode' }
    $fmmMode = if ($Kind -eq 'direct') { 'off' } else { 'single-matvec-analytic' }
    $precond = if ($Kind -eq 'fmm-block') { 'fmm-block-jacobi' } else { 'diagonal' }
    $auditRows = if ($Kind -eq 'direct') { 0 } else { $AuditSampleCount }
    $arg = Write-FmmCliArgFile -RunDir $CaseRoot -RepoRoot $RepoRoot `
        -PqrPath $PqrPath -Backend tabi -Mode $mode -FmmMode $fmmMode `
        -SurfacePoints $SurfacePoints -MaxIter $MaxIter -Restart $Restart `
        -Tolerance $Tolerance -Threads $Threads -SaltM 0.15 `
        -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK $TemperatureK -OpeningTheta $OpeningTheta `
        -AuditSampleCount $auditRows -Preconditioner $precond
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $CaseRoot `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $result = if ($run.exit_code -eq 0) {
        $expected = if ($Kind -eq 'direct') { 'off' } else { 'single-matvec-analytic' }
        Read-FmmCliCaseResult -RunDir $CaseRoot -ExpectedFmmMode $expected `
            -RequireAudit ($Kind -ne 'direct')
    } else {
        New-ScopeUEmptyResult "exit_$($run.exit_code)"
    }
    $precondProfile = Read-ScopeUPrecondProfile -RunDir $CaseRoot
    return [pscustomobject]@{
        run=$run; result=$result; preconditioner=$precondProfile
    }
}

function Get-ScopeUSpeedup([object]$Direct, [object]$Candidate) {
    if (-not (Test-ScopeUFinite $Direct.matvec_seconds_avg) -or
        -not (Test-ScopeUFinite $Candidate.matvec_seconds_avg) -or
        $Candidate.matvec_seconds_avg -le 0.0) {
        return [double]::NaN
    }
    return $Direct.matvec_seconds_avg / $Candidate.matvec_seconds_avg
}

function Get-ScopeUEnergyDelta([object]$Direct, [object]$Candidate) {
    if (-not (Test-ScopeUFinite $Direct.energy_kj_mol) -or
        -not (Test-ScopeUFinite $Candidate.energy_kj_mol)) {
        return [double]::NaN
    }
    return [Math]::Abs($Candidate.energy_kj_mol - $Direct.energy_kj_mol)
}

function Get-ScopeUClass(
    [object]$Direct,
    [object]$Diagonal,
    [object]$Candidate,
    [object]$CandidateRun
) {
    if ($CandidateRun.exit_code -eq 124) { return 'timeout_not_verified' }
    if ($CandidateRun.exit_code -ne 0 -or $Candidate.status -eq $null) {
        return 'status_fail'
    }
    if ($Candidate.fmm_m2l_mode -ne 'single-matvec-analytic' -or
        $Candidate.fmm_far_pairs -le 0 -or
        $Candidate.fmm_finite_outputs -ne 'true') {
        return 'telemetry_fail'
    }
    if (-not (Test-ScopeUFinite $Candidate.tree_vs_direct_matvec_rel) -or
        $Candidate.tree_vs_direct_matvec_rel -gt 1.0e-2) {
        return 'accuracy_limited'
    }
    $speedup = Get-ScopeUSpeedup -Direct $Direct -Candidate $Candidate
    if (-not (Test-ScopeUFinite $speedup) -or $speedup -lt 1.2) {
        return 'speed_limited'
    }
    if ($Candidate.status -eq 'ok' -and
        (Test-ScopeUFinite $Candidate.true_residual) -and
        $Candidate.true_residual -le 1.0e-6) {
        return 'strict_converged'
    }
    if ((Test-ScopeUFinite $Diagonal.true_residual) -and
        (Test-ScopeUFinite $Candidate.true_residual) -and
        $Candidate.true_residual -gt 0.0) {
        $improvement = $Diagonal.true_residual / $Candidate.true_residual
        if ($improvement -ge 4.0) { return 'residual_improved_4x' }
    }
    return 'not_improved'
}

$fixtureRoot = Join-Path $workRootFull 'fixture'
New-Item -ItemType Directory -Path $fixtureRoot -Force | Out-Null
$fixture = Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot -RunDir $fixtureRoot
if ([string]::IsNullOrWhiteSpace($fixture) -or -not (Test-Path -LiteralPath $fixture)) {
    throw 'Barnase fixture could not be resolved for Scope U.'
}

$rows = New-Object System.Collections.Generic.List[object]
$conditions = @(
    [pscustomobject]@{ name='cap128_240_r64_1e-6'; cap=128; max_iter=240; restart=64; tol=1.0e-6 }
)
if (-not $Quick) {
    $conditions += @(
        [pscustomobject]@{ name='cap256_240_r64_1e-6'; cap=256; max_iter=240; restart=64; tol=1.0e-6 },
        [pscustomobject]@{ name='cap256_360_r96_1e-6'; cap=256; max_iter=360; restart=96; tol=1.0e-6 },
        [pscustomobject]@{ name='cap256_480_r128_1e-6'; cap=256; max_iter=480; restart=128; tol=1.0e-6 }
    )
}

foreach ($condition in $conditions) {
    $sampleDir = Join-Path $workRootFull "sample_$($condition.name)"
    New-Item -ItemType Directory -Path $sampleDir -Force | Out-Null
    $sampledPqr = Join-Path $sampleDir 'barnase_sampled.pqr'
    $sample = Write-FmmCliSampledPqr -InputPath $fixture -OutputPath $sampledPqr `
        -Cap $condition.cap

    $direct = Invoke-ScopeUCase -CaseRoot (Join-Path $workRootFull "$($condition.name)_direct") `
        -PqrPath $sampledPqr -Kind direct -MaxIter $condition.max_iter `
        -Restart $condition.restart -Tolerance $condition.tol
    $diagonal = Invoke-ScopeUCase -CaseRoot (Join-Path $workRootFull "$($condition.name)_fmm_diagonal") `
        -PqrPath $sampledPqr -Kind fmm-diagonal -MaxIter $condition.max_iter `
        -Restart $condition.restart -Tolerance $condition.tol
    $block = Invoke-ScopeUCase -CaseRoot (Join-Path $workRootFull "$($condition.name)_fmm_block") `
        -PqrPath $sampledPqr -Kind fmm-block -MaxIter $condition.max_iter `
        -Restart $condition.restart -Tolerance $condition.tol

    foreach ($case in @(
            [pscustomobject]@{ label='direct'; data=$direct; class='direct_reference' },
            [pscustomobject]@{ label='fmm_diagonal'; data=$diagonal; class='diagonal_reference' },
            [pscustomobject]@{
                label='fmm_block_jacobi'; data=$block
                class=(Get-ScopeUClass -Direct $direct.result `
                    -Diagonal $diagonal.result -Candidate $block.result `
                    -CandidateRun $block.run)
            })) {
        $speedup = Get-ScopeUSpeedup -Direct $direct.result `
            -Candidate $case.data.result
        $energyDelta = Get-ScopeUEnergyDelta -Direct $direct.result `
            -Candidate $case.data.result
        $rows.Add([pscustomobject]@{
            condition=$condition.name
            atom_cap=$condition.cap
            sampled_atoms=$sample.sampled_atoms
            backend=$case.label
            max_iter=$condition.max_iter
            restart=$condition.restart
            tolerance=$condition.tol
            exit_code=$case.data.run.exit_code
            status=$case.data.result.status
            panels=$case.data.result.panels
            iterations=$case.data.result.iterations
            true_residual=$case.data.result.true_residual
            energy_kj_mol=$case.data.result.energy_kj_mol
            energy_delta_abs=$energyDelta
            audit_rel=$case.data.result.tree_vs_direct_matvec_rel
            matvec_seconds_avg=$case.data.result.matvec_seconds_avg
            speedup_vs_direct=$speedup
            wall_seconds=$case.data.run.wall_seconds
            peak_working_set_mb=$case.data.run.peak_working_set_mb
            fmm_m2l_mode=$case.data.result.fmm_m2l_mode
            fmm_far_pairs=$case.data.result.fmm_far_pairs
            preconditioner_requested=$case.data.preconditioner.requested
            preconditioner_effective=$case.data.preconditioner.effective
            preconditioner_estimated_bytes=$case.data.preconditioner.estimated_bytes
            preconditioner_valid_blocks=$case.data.preconditioner.valid_blocks
            preconditioner_valid_panels=$case.data.preconditioner.valid_panels
            preconditioner_quality_status=$case.data.preconditioner.quality_status
            preconditioner_setup_seconds=$case.data.preconditioner.setup_seconds
            preconditioner_apply_seconds=$case.data.preconditioner.apply_seconds
            preconditioner_apply_calls=$case.data.preconditioner.apply_calls
            scope_u_class=$case.class
        }) | Out-Null
    }
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation

$blockRows = @($rows | Where-Object { $_.backend -eq 'fmm_block_jacobi' })
$strictRows = @($blockRows | Where-Object { $_.scope_u_class -eq 'strict_converged' })
$improvedRows = @($blockRows | Where-Object {
        $_.scope_u_class -eq 'strict_converged' -or
        $_.scope_u_class -eq 'residual_improved_4x'
    })
$finalClass = if ($strictRows.Count -gt 0) {
    'strict_converged'
} elseif ($improvedRows.Count -gt 0) {
    'residual_improved'
} elseif (($blockRows | Where-Object { $_.scope_u_class -eq 'timeout_not_verified' }).Count -gt 0) {
    'timeout_not_verified'
} else {
    'not_improved'
}

[pscustomobject]@{
    fixture='barnase'
    mode=if ($Quick) { 'quick' } else { 'full' }
    rows=$rows.Count
    block_rows=$blockRows.Count
    strict_rows=$strictRows.Count
    improved_rows=$improvedRows.Count
    scope_u_decision=$finalClass
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'Scope U FMM-native Block-Jacobi Preconditioner Probe',
    "mode: $(if ($Quick) { 'quick' } else { 'full' })",
    "rows: $($rows.Count)",
    "decision: $finalClass",
    'production_default: grid_remains_default',
    'fmm_m2l_status: hidden_experimental'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = New-Object System.Collections.Generic.List[string]
$report.Add('# TITAN FMM/M2L Scope U FMM Block Preconditioner Report')
$report.Add('')
$report.Add("mode: $(if ($Quick) { 'quick' } else { 'full' })")
$report.Add("decision: $finalClass")
$report.Add('')
$report.Add('Rows:')
foreach ($row in $rows) {
    $report.Add(('{0} {1} status={2} residual={3} audit={4} speedup={5} precond={6} class={7}' -f
            $row.condition, $row.backend, $row.status, $row.true_residual,
            $row.audit_rel, $row.speedup_vs_direct,
            $row.preconditioner_effective, $row.scope_u_class))
}
$report.Add('')
$report.Add('Notes:')
$report.Add('- FMM/M2L remains hidden and grid remains the production default.')
$report.Add('- fmm-block-jacobi is an experimental native preconditioner built from FMM context leaf blocks.')
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

$badArgs = @(Get-ChildItem -LiteralPath $workRootFull -Recurse -Filter arg.in |
    ForEach-Object {
        $path = $_.FullName
        Get-Content -LiteralPath $path | Where-Object {
            -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
        } | ForEach-Object { "$path :: $_" }
    })
if ($badArgs.Count -ne 0) {
    throw "Generated arg.in contains non equals-syntax options: $($badArgs -join '; ')"
}

Write-Host "Scope U report: $ReportPath"
Write-Host "Scope U decision: $finalClass"
