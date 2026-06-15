param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_z5'),
    [switch]$Quick,
    [switch]$IncludeContextRows,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 4
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fmm_m2l_scope_z5_matrix.csv'
$TraceCsv = Join-Path $WorkRoot 'fmm_m2l_scope_z5_residual_trace.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_z5_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_z5_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fmm_m2l_scope_z5_money_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z5_PRECONDITIONER_MONEY_PLOT_REPORT.txt'

if ($Quick) { $CaseTimeoutSeconds = [Math]::Min($CaseTimeoutSeconds, 900) }
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope Z5 work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function New-Z5EmptyResult([string]$Status) {
    return [pscustomobject]@{
        status=$Status; panels=0; iterations=0
        final_residual=[double]::NaN; true_residual=[double]::NaN
        residual_history=''; energy_kj_mol=[double]::NaN
        matvec_seconds_avg=[double]::NaN; total_seconds=[double]::NaN
        fmm_m2l_mode=''; fmm_near_pairs=0; fmm_far_pairs=0
        fmm_setup_seconds=[double]::NaN; fmm_apply_seconds=[double]::NaN
        fmm_finite_outputs='false'; tree_vs_direct_matvec_rel=[double]::NaN
        matvec_audit_rows=0; matvec_audit_seconds=[double]::NaN
        reaction_energy_source=''; timing_backend=''; memory_bytes=[double]::NaN
    }
}

function Read-Z5PrecondProfile([string]$RunDir) {
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $profilePath)) {
        return [pscustomobject]@{
            requested='NA'; effective='NA'; ilu_neighbor_cap=0
            estimated_bytes=[double]::NaN; valid_blocks=0; valid_panels=0
            setup_seconds=[double]::NaN; apply_seconds=[double]::NaN
            apply_calls=0; sparse_nnz_blocks=0; sparse_avg_degree=[double]::NaN
            sparse_max_degree=0; quality_status='NA'
        }
    }
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    return [pscustomobject]@{
        requested=[string]$profile.preconditioner_requested
        effective=[string]$profile.preconditioner_effective
        ilu_neighbor_cap=Convert-FmmCliInt64 $profile.preconditioner_ilu_neighbor_cap
        estimated_bytes=Convert-FmmCliDouble $profile.preconditioner_estimated_bytes
        valid_blocks=Convert-FmmCliInt64 $profile.preconditioner_valid_blocks
        valid_panels=Convert-FmmCliInt64 $profile.preconditioner_valid_panels
        setup_seconds=Convert-FmmCliDouble $profile.preconditioner_setup_seconds
        apply_seconds=Convert-FmmCliDouble $profile.preconditioner_apply_seconds
        apply_calls=Convert-FmmCliInt64 $profile.preconditioner_apply_calls
        sparse_nnz_blocks=Convert-FmmCliInt64 $profile.preconditioner_sparse_nnz_blocks
        sparse_avg_degree=Convert-FmmCliDouble $profile.preconditioner_sparse_avg_degree
        sparse_max_degree=Convert-FmmCliInt64 $profile.preconditioner_sparse_max_degree
        quality_status=[string]$profile.preconditioner_quality_status
    }
}

function Write-Z5BlockedOutputs([string]$Decision, [string]$Reason) {
    [pscustomobject]@{
        scope_z5_decision=$Decision
        reason=$Reason
        matrix=$MatrixCsv
        trace=$TraceCsv
        plot=$PlotPath
    } | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation
    Set-Content -LiteralPath $SummaryPath -Encoding ASCII -Value @(
        "scope_z5_decision=$Decision",
        "reason=$Reason"
    )
    Set-Content -LiteralPath $ReportPath -Encoding ASCII -Value @(
        'TITAN FMM/M2L Scope Z5 Preconditioner Money Plot Report',
        "Decision: $Decision",
        "Reason: $Reason"
    )
}

$fixtureSourceDir = Join-Path $workRootFull 'fixture_source'
New-Item -ItemType Directory -Path $fixtureSourceDir -Force | Out-Null
$fixturePath = Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot `
    -RunDir $fixtureSourceDir
if ([string]::IsNullOrWhiteSpace($fixturePath) -or
    -not (Test-Path -LiteralPath $fixturePath)) {
    Write-Z5BlockedOutputs -Decision 'blocked_missing_fixture' `
        -Reason 'Barnase fixture could not be resolved.'
    throw 'Scope Z5 blocked: Barnase fixture could not be resolved.'
}

$atomCap = if ($Quick) { 128 } else { 256 }
$fixtureDir = Join-Path $workRootFull 'fixture'
New-Item -ItemType Directory -Path $fixtureDir -Force | Out-Null
$samplePqr = Join-Path $fixtureDir 'barnase_sampled.pqr'
$sampleInfo = Write-FmmCliSampledPqr -InputPath $fixturePath `
    -OutputPath $samplePqr -Cap $atomCap

$requiredCases = @(
    [pscustomobject]@{
        label='block_diagonal'; precond='fmm-block-jacobi'; ilu_cap=64
        block_cap=32; overlap_cap=64; max_iter=720; restart=128
        tolerance=1.0e-6
    },
    [pscustomobject]@{
        label='sparse_ilu0'; precond='fmm-nearfield-ilu0'; ilu_cap=64
        block_cap=32; overlap_cap=64; max_iter=720; restart=128
        tolerance=1.0e-6
    }
)
$contextCases = @(
    [pscustomobject]@{
        label='diagonal'; precond='diagonal'; ilu_cap=64; block_cap=32
        overlap_cap=64; max_iter=720; restart=128; tolerance=1.0e-6
    },
    [pscustomobject]@{
        label='overlap_schwarz'; precond='fmm-overlap-schwarz'; ilu_cap=64
        block_cap=32; overlap_cap=64; max_iter=720; restart=128
        tolerance=1.0e-6
    }
)
$cases = @($requiredCases)
if ((-not $Quick) -and $IncludeContextRows) {
    $cases += $contextCases
}

function Invoke-Z5Case([object]$Case) {
    $runDir = Join-Path $workRootFull $Case.label
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $argPath = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
        -PqrPath $samplePqr -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -SurfacePoints 8 `
        -MaxIter $Case.max_iter -Restart $Case.restart `
        -Tolerance $Case.tolerance -Threads $Threads -SaltM 0.15 `
        -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK 298.15 -OpeningTheta 0.6 -AuditSampleCount 512 `
        -Preconditioner $Case.precond `
        -PreconditionerBlockPanelCap $Case.block_cap `
        -PreconditionerOverlapPanelCap $Case.overlap_cap `
        -PreconditionerIluNeighborCap $Case.ilu_cap `
        -PreconditionerIluDropTol 0.0
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $readError = ''
    $result = $null
    if ($run.exit_code -eq 0 -and -not $run.timed_out) {
        try {
            $result = Read-FmmCliCaseResult -RunDir $runDir `
                -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
        } catch {
            $readError = $_.Exception.Message
        }
    }
    if ($null -eq $result) {
        $status = if ($run.timed_out) { 'timeout' }
            elseif (-not [string]::IsNullOrWhiteSpace($readError)) { 'read_error' }
            else { "exit_$($run.exit_code)" }
        $result = New-Z5EmptyResult $status
    }
    return [pscustomobject]@{
        case=$Case
        run=$run
        result=$result
        preconditioner=Read-Z5PrecondProfile -RunDir $runDir
        read_error=$readError
        run_dir=$runDir
    }
}

function Get-Z5Status([object]$Record) {
    if ($Record.run.timed_out) { return 'timeout' }
    if (-not [string]::IsNullOrWhiteSpace($Record.read_error)) {
        return 'read_error'
    }
    return $Record.result.status
}

$records = New-Object System.Collections.Generic.List[object]
foreach ($case in $cases) {
    $records.Add((Invoke-Z5Case $case)) | Out-Null
}

$baseRows = foreach ($record in $records) {
    [pscustomobject]@{
        fixture='barnase'
        atom_cap=$atomCap
        sampled_atoms=$sampleInfo.sampled_atoms
        preconditioner_label=$record.case.label
        preconditioner_requested=$record.case.precond
        preconditioner_effective=$record.preconditioner.effective
        status=Get-Z5Status $record
        read_error=$record.read_error
        exit_code=$record.run.exit_code
        timed_out=$record.run.timed_out
        wall_seconds=$record.run.wall_seconds
        peak_working_set_mb=$record.run.peak_working_set_mb
        panels=$record.result.panels
        iterations=$record.result.iterations
        true_residual=$record.result.true_residual
        final_residual=$record.result.final_residual
        energy_kj_mol=$record.result.energy_kj_mol
        matvec_seconds_avg=$record.result.matvec_seconds_avg
        audit_rel=$record.result.tree_vs_direct_matvec_rel
        matvec_audit_rows=$record.result.matvec_audit_rows
        fmm_m2l_mode=$record.result.fmm_m2l_mode
        fmm_near_pairs=$record.result.fmm_near_pairs
        fmm_far_pairs=$record.result.fmm_far_pairs
        preconditioner_setup_seconds=$record.preconditioner.setup_seconds
        preconditioner_apply_seconds=$record.preconditioner.apply_seconds
        preconditioner_apply_calls=$record.preconditioner.apply_calls
        sparse_nnz_blocks=$record.preconditioner.sparse_nnz_blocks
        sparse_avg_degree=$record.preconditioner.sparse_avg_degree
        sparse_max_degree=$record.preconditioner.sparse_max_degree
        quality_status=$record.preconditioner.quality_status
        residual_history=$record.result.residual_history
    }
}
$blockRowForRatio = @($baseRows | Where-Object {
        $_.preconditioner_label -eq 'block_diagonal'
    }) | Select-Object -First 1
$diagRowForRatio = @($baseRows | Where-Object {
        $_.preconditioner_label -eq 'diagonal'
    }) | Select-Object -First 1
$blockResidualForRatio = if ($null -ne $blockRowForRatio) {
    [double]$blockRowForRatio.true_residual
} else {
    [double]::NaN
}
$diagResidualForRatio = if ($null -ne $diagRowForRatio) {
    [double]$diagRowForRatio.true_residual
} else {
    [double]::NaN
}
$rows = foreach ($row in $baseRows) {
    $role = switch ($row.preconditioner_label) {
        'block_diagonal' { 'required_block_baseline' }
        'sparse_ilu0' { 'required_sparse_ilu0' }
        default { 'context' }
    }
    $ratioVsBlock = [double]::NaN
    if ((Test-FmmCliFinite $blockResidualForRatio) -and
        (Test-FmmCliFinite ([double]$row.true_residual)) -and
        ([double]$row.true_residual -gt 0.0)) {
        $ratioVsBlock = $blockResidualForRatio / [double]$row.true_residual
    }
    $ratioVsDiagonal = [double]::NaN
    if ((Test-FmmCliFinite $diagResidualForRatio) -and
        (Test-FmmCliFinite ([double]$row.true_residual)) -and
        ([double]$row.true_residual -gt 0.0)) {
        $ratioVsDiagonal = $diagResidualForRatio / [double]$row.true_residual
    }
    [pscustomobject]@{
        fixture=$row.fixture
        atom_cap=$row.atom_cap
        sampled_atoms=$row.sampled_atoms
        preconditioner_label=$row.preconditioner_label
        money_plot_role=$role
        preconditioner_requested=$row.preconditioner_requested
        preconditioner_effective=$row.preconditioner_effective
        status=$row.status
        read_error=$row.read_error
        exit_code=$row.exit_code
        timed_out=$row.timed_out
        wall_seconds=$row.wall_seconds
        peak_working_set_mb=$row.peak_working_set_mb
        panels=$row.panels
        iterations=$row.iterations
        true_residual=$row.true_residual
        final_residual=$row.final_residual
        residual_ratio_vs_block=$ratioVsBlock
        residual_ratio_vs_diagonal=$ratioVsDiagonal
        energy_kj_mol=$row.energy_kj_mol
        matvec_seconds_avg=$row.matvec_seconds_avg
        audit_rel=$row.audit_rel
        matvec_audit_rows=$row.matvec_audit_rows
        fmm_m2l_mode=$row.fmm_m2l_mode
        fmm_near_pairs=$row.fmm_near_pairs
        fmm_far_pairs=$row.fmm_far_pairs
        preconditioner_setup_seconds=$row.preconditioner_setup_seconds
        preconditioner_apply_seconds=$row.preconditioner_apply_seconds
        preconditioner_apply_calls=$row.preconditioner_apply_calls
        sparse_nnz_blocks=$row.sparse_nnz_blocks
        sparse_avg_degree=$row.sparse_avg_degree
        sparse_max_degree=$row.sparse_max_degree
        quality_status=$row.quality_status
        residual_history=$row.residual_history
    }
}
$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation

$traceRows = New-Object System.Collections.Generic.List[object]
foreach ($row in $rows) {
    foreach ($point in (Convert-FmmCliResidualHistory ([string]$row.residual_history))) {
        $logValue = if ($point.residual -gt 0.0) {
            [Math]::Log10($point.residual)
        } else {
            [double]::NaN
        }
        $traceRows.Add([pscustomobject]@{
            fixture=$row.fixture
            atom_cap=$row.atom_cap
            preconditioner_label=$row.preconditioner_label
            iteration=$point.iteration
            residual=$point.residual
            log10_residual=$logValue
        }) | Out-Null
    }
}
$traceRows | Export-Csv -LiteralPath $TraceCsv -NoTypeInformation

function Write-Z5SvgPlot([object[]]$TraceRows, [string]$Path) {
    $finite = @($TraceRows | Where-Object {
        Test-FmmCliFinite ([double]$_.log10_residual)
    })
    if ($finite.Count -eq 0) {
        Set-Content -LiteralPath $Path -Encoding ASCII -Value '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="520"><text x="20" y="40">No finite residual trace</text></svg>'
        return
    }
    $width = 900; $height = 520
    $left = 80; $right = 30; $top = 30; $bottom = 70
    $maxIter = [double](@($finite | Measure-Object -Property iteration -Maximum)[0].Maximum)
    $minLog = [double](@($finite | Measure-Object -Property log10_residual -Minimum)[0].Minimum)
    $maxLog = [double](@($finite | Measure-Object -Property log10_residual -Maximum)[0].Maximum)
    if ($maxIter -le 0.0) { $maxIter = 1.0 }
    if ([Math]::Abs($maxLog - $minLog) -lt 1.0e-12) {
        $maxLog += 1.0
        $minLog -= 1.0
    }
    $colors = @{
        diagonal = '#666666'
        block_diagonal = '#d55e00'
        overlap_schwarz = '#0072b2'
        sparse_ilu0 = '#009e73'
    }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add("<svg xmlns=`"http://www.w3.org/2000/svg`" width=`"$width`" height=`"$height`" viewBox=`"0 0 $width $height`">")
    $svg.Add('<rect width="100%" height="100%" fill="white"/>')
    $svg.Add("<line x1=`"$left`" y1=`"$($height-$bottom)`" x2=`"$($width-$right)`" y2=`"$($height-$bottom)`" stroke=`"black`"/>")
    $svg.Add("<line x1=`"$left`" y1=`"$top`" x2=`"$left`" y2=`"$($height-$bottom)`" stroke=`"black`"/>")
    foreach ($label in @('diagonal','block_diagonal','overlap_schwarz','sparse_ilu0')) {
        $points = @($finite | Where-Object { $_.preconditioner_label -eq $label } |
            Sort-Object {[int]$_.iteration})
        if ($points.Count -eq 0) { continue }
        $coords = foreach ($p in $points) {
            $x = $left + (($width - $left - $right) *
                ([double]$p.iteration / $maxIter))
            $y = $top + (($height - $top - $bottom) *
                (($maxLog - [double]$p.log10_residual) / ($maxLog - $minLog)))
            ('{0:F2},{1:F2}' -f $x,$y)
        }
        $color = $colors[$label]
        $svg.Add("<polyline fill=`"none`" stroke=`"$color`" stroke-width=`"2.4`" points=`"$($coords -join ' ')`"/>")
    }
    $legendY = 60
    foreach ($label in @('diagonal','block_diagonal','overlap_schwarz','sparse_ilu0')) {
        $color = $colors[$label]
        $svg.Add("<line x1=`"610`" y1=`"$legendY`" x2=`"650`" y2=`"$legendY`" stroke=`"$color`" stroke-width=`"3`"/>")
        $svg.Add("<text x=`"660`" y=`"$($legendY+5)`" font-family=`"Arial`" font-size=`"14`">$label</text>")
        $legendY += 24
    }
    $svg.Add('<text x="390" y="500" font-family="Arial" font-size="16">GMRES residual-history index</text>')
    $svg.Add('<text x="20" y="280" font-family="Arial" font-size="16" transform="rotate(-90 20,280)">log10 residual</text>')
    $svg.Add('<text x="80" y="24" font-family="Arial" font-size="18" font-weight="bold">Scope Z5 preconditioner convergence</text>')
    $svg.Add('</svg>')
    Set-Content -LiteralPath $Path -Encoding ASCII -Value $svg
}

Write-Z5SvgPlot -TraceRows $traceRows.ToArray() -Path $PlotPath

$ilu = @($rows | Where-Object { $_.preconditioner_label -eq 'sparse_ilu0' }) |
    Select-Object -First 1
$block = @($rows | Where-Object { $_.preconditioner_label -eq 'block_diagonal' }) |
    Select-Object -First 1
$z5Decision = 'no_decision'
$iluStrict = ($null -ne $ilu) -and
    ($ilu.status -eq 'ok') -and
    (Test-FmmCliFinite ([double]$ilu.true_residual)) -and
    ([double]$ilu.true_residual -le 1.0e-6) -and
    (Test-FmmCliFinite ([double]$ilu.audit_rel)) -and
    ([double]$ilu.audit_rel -le 1.0e-2)
$blockStrict = ($null -ne $block) -and
    ($block.status -eq 'ok') -and
    (Test-FmmCliFinite ([double]$block.true_residual)) -and
    ([double]$block.true_residual -le 1.0e-6)
if (@($rows | Where-Object { $_.status -eq 'read_error' }).Count -gt 0 -or
    $traceRows.Count -eq 0) {
    $z5Decision = 'blocked_profile_schema'
} elseif ($null -eq $ilu -or $null -eq $block) {
    $z5Decision = 'missing_required_rows'
} elseif ($ilu.status -eq 'timeout' -or $block.status -eq 'timeout') {
    $z5Decision = 'blocked_runtime'
} elseif ($ilu.status -ne 'ok') {
    $z5Decision = 'ilu_not_ok'
} elseif ($iluStrict -and (-not $blockStrict)) {
    $z5Decision = if ($Quick) {
        'money_plot_harness_ready'
    } else {
        'money_plot_pass'
    }
} elseif ($iluStrict) {
    $z5Decision = if ($Quick) {
        'harness_ready_but_no_contrast'
    } else {
        'no_contrast'
    }
} else {
    $z5Decision = 'ilu_not_strict'
}

[pscustomobject]@{
    scope_z5_decision=$z5Decision
    atom_cap=$atomCap
    sampled_atoms=$sampleInfo.sampled_atoms
    rows=$rows.Count
    trace_rows=$traceRows.Count
    required_rows='block_diagonal,sparse_ilu0'
    context_rows=if ($IncludeContextRows -and (-not $Quick)) {
        'diagonal,overlap_schwarz'
    } else {
        'not_run'
    }
    matrix=$MatrixCsv
    trace=$TraceCsv
    plot=$PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$lines = New-Object System.Collections.Generic.List[string]
$lines.Add('TITAN FMM/M2L Scope Z5 Preconditioner Money Plot Report')
$lines.Add("Decision: $z5Decision")
$lines.Add("Fixture: Barnase sampled cap $atomCap ($($sampleInfo.sampled_atoms) atoms)")
$lines.Add("Required rows: block_diagonal, sparse_ilu0")
$lines.Add("Context rows: $(if ($IncludeContextRows -and (-not $Quick)) { 'diagonal, overlap_schwarz' } else { 'not run' })")
$lines.Add('Direct verification source: Scope X.1 existing cap-256 direct-verification evidence; Z5B compares preconditioners on the same FMM/M2L operator.')
$lines.Add('')
$lines.Add('Rows:')
foreach ($row in $rows) {
    $lines.Add(('{0}: role={1} status={2} effective={3} iter={4} residual={5:G6} ratio_vs_block={6:G6} audit={7:G6} wall={8:G5}s peak={9:G5}MB' -f
        $row.preconditioner_label, $row.money_plot_role, $row.status,
        $row.preconditioner_effective,
        [int]$row.iterations, [double]$row.true_residual,
        [double]$row.residual_ratio_vs_block,
        [double]$row.audit_rel, [double]$row.wall_seconds,
        [double]$row.peak_working_set_mb))
}
$lines.Add('')
$lines.Add('Publication framing:')
$lines.Add('Purely local block-diagonal preconditioning is expected to stall on dense boundary-integral PB systems, while TITAN sparse near-field ILU(0) captures cross-leaf near-field coupling and restores Krylov convergence.')
$lines.Add('')
$lines.Add('Z5B note: full cap-256 decision uses the required Block-Diagonal vs sparse ILU(0) rows first; diagonal and overlap rows are context only.')
$lines.Add("Matrix: $MatrixCsv")
$lines.Add("Trace: $TraceCsv")
$lines.Add("Plot: $PlotPath")
Set-Content -LiteralPath $ReportPath -Encoding ASCII -Value $lines
Set-Content -LiteralPath $SummaryPath -Encoding ASCII -Value @(
    "scope_z5_decision=$z5Decision",
    "fixture=barnase",
    "atom_cap=$atomCap",
    "sampled_atoms=$($sampleInfo.sampled_atoms)",
    "required_rows=block_diagonal,sparse_ilu0",
    "context_rows=$(if ($IncludeContextRows -and (-not $Quick)) { 'diagonal,overlap_schwarz' } else { 'not_run' })",
    "matrix=$MatrixCsv",
    "trace=$TraceCsv",
    "plot=$PlotPath",
    "report=$ReportPath",
    'note=Z5B uses required preconditioner contrast rows first; direct verification comes from Scope X.1 evidence.'
)

Write-Host "Scope Z5 report: $ReportPath"
