param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_z5c'),
    [switch]$Quick,
    [switch]$IncludeContextRows,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 4,
    [int]$BarnaseRepeats = 3
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fmm_m2l_scope_z5c_matrix.csv'
$TraceCsv = Join-Path $WorkRoot 'fmm_m2l_scope_z5c_residual_trace.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_z5c_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_z5c_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fmm_m2l_scope_z5c_money_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z5C_MULTI_FIXTURE_MONEY_PLOT_REPORT.txt'

if ($Quick) {
    $CaseTimeoutSeconds = [Math]::Min($CaseTimeoutSeconds, 900)
    $BarnaseRepeats = 1
}
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope Z5C work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function New-Z5cEmptyResult([string]$Status) {
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

function Read-Z5cPrecondProfile([string]$RunDir) {
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

function Get-Z5cStatus([object]$Record) {
    if ($Record.run.timed_out) { return 'timeout' }
    if (-not [string]::IsNullOrWhiteSpace($Record.read_error)) {
        return 'read_error'
    }
    return $Record.result.status
}

function Get-Z5cMedian([double[]]$Values) {
    $finite = @($Values | Where-Object { Test-FmmCliFinite ([double]$_) } |
        Sort-Object)
    if ($finite.Count -eq 0) { return [double]::NaN }
    $mid = [int][Math]::Floor($finite.Count / 2)
    if (($finite.Count % 2) -eq 1) { return [double]$finite[$mid] }
    return (([double]$finite[$mid - 1]) + ([double]$finite[$mid])) / 2.0
}

function Get-Z5cMin([double[]]$Values) {
    $finite = @($Values | Where-Object { Test-FmmCliFinite ([double]$_) })
    if ($finite.Count -eq 0) { return [double]::NaN }
    return [double](@($finite | Measure-Object -Minimum)[0].Minimum)
}

function Get-Z5cMax([double[]]$Values) {
    $finite = @($Values | Where-Object { Test-FmmCliFinite ([double]$_) })
    if ($finite.Count -eq 0) { return [double]::NaN }
    return [double](@($finite | Measure-Object -Maximum)[0].Maximum)
}

function Write-Z5cBlockedOutputs([string]$Decision, [string]$Reason) {
    [pscustomobject]@{
        scope_z5c_decision=$Decision
        reason=$Reason
        matrix=$MatrixCsv
        trace=$TraceCsv
        plot=$PlotPath
    } | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation
    Set-Content -LiteralPath $SummaryPath -Encoding ASCII -Value @(
        "scope_z5c_decision=$Decision",
        "reason=$Reason"
    )
    Set-Content -LiteralPath $ReportPath -Encoding ASCII -Value @(
        'TITAN FMM/M2L Scope Z5C Multi-Fixture Money Plot Report',
        "Decision: $Decision",
        "Reason: $Reason"
    )
}

function New-Z5cFixtureSpec(
    [string]$Name,
    [int]$AtomCap,
    [int]$Repeats,
    [bool]$Gated = $false
) {
    return [pscustomobject]@{
        name=$Name
        atom_cap=$AtomCap
        repeats=$Repeats
        gated=$Gated
    }
}

$fixtureSpecs = if ($Quick) {
    @((New-Z5cFixtureSpec -Name 'barnase' -AtomCap 128 -Repeats 1))
} else {
    @(
        (New-Z5cFixtureSpec -Name 'barnase' -AtomCap 256 -Repeats $BarnaseRepeats),
        (New-Z5cFixtureSpec -Name 'research-medium' -AtomCap 256 -Repeats 1),
        (New-Z5cFixtureSpec -Name 'research-alt' -AtomCap 256 -Repeats 1),
        (New-Z5cFixtureSpec -Name 'research-large' -AtomCap 256 -Repeats 1 -Gated $true)
    )
}

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
if ((-not $Quick) -and $IncludeContextRows) { $cases += $contextCases }

function Invoke-Z5cCase([object]$Fixture, [object]$Case, [int]$Repeat, [string]$SamplePqr) {
    $runDir = Join-Path $workRootFull (Join-Path $Fixture.name ("r{0}_{1}" -f $Repeat,$Case.label))
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $argPath = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
        -PqrPath $SamplePqr -Backend tabi -Mode treecode `
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
        $result = New-Z5cEmptyResult $status
    }
    return [pscustomobject]@{
        fixture=$Fixture
        case=$Case
        repeat=$Repeat
        run=$run
        result=$result
        preconditioner=Read-Z5cPrecondProfile -RunDir $runDir
        read_error=$readError
        run_dir=$runDir
    }
}

function New-Z5cSkippedRecord([object]$Fixture, [object]$Case, [int]$Repeat, [string]$Status) {
    return [pscustomobject]@{
        fixture=$Fixture
        case=$Case
        repeat=$Repeat
        run=[pscustomobject]@{
            exit_code=0; timed_out=$false; wall_seconds=0.0
            peak_working_set_bytes=0; peak_working_set_mb=[double]::NaN
        }
        result=(New-Z5cEmptyResult $Status)
        preconditioner=[pscustomobject]@{
            requested=$Case.precond; effective='NA'; ilu_neighbor_cap=0
            estimated_bytes=[double]::NaN; valid_blocks=0; valid_panels=0
            setup_seconds=[double]::NaN; apply_seconds=[double]::NaN
            apply_calls=0; sparse_nnz_blocks=0; sparse_avg_degree=[double]::NaN
            sparse_max_degree=0; quality_status='NA'
        }
        read_error=''
        run_dir=''
    }
}

$records = New-Object System.Collections.Generic.List[object]
$sampleInfoByFixture = @{}
$nonGatedTimeout = $false
foreach ($fixture in $fixtureSpecs) {
    if ($fixture.gated -and $nonGatedTimeout) {
        foreach ($case in $cases) {
            $records.Add((New-Z5cSkippedRecord -Fixture $fixture -Case $case `
                    -Repeat 1 -Status 'skipped_gated_timeout')) | Out-Null
        }
        continue
    }
    $fixtureSourceDir = Join-Path $workRootFull (Join-Path $fixture.name 'fixture_source')
    New-Item -ItemType Directory -Path $fixtureSourceDir -Force | Out-Null
    $fixturePath = Resolve-FmmCliFixture -Fixture $fixture.name -RepoRoot $RepoRoot `
        -RunDir $fixtureSourceDir
    if ([string]::IsNullOrWhiteSpace($fixturePath) -or
        -not (Test-Path -LiteralPath $fixturePath)) {
        Write-Z5cBlockedOutputs -Decision 'blocked_missing_fixture' `
            -Reason "Fixture could not be resolved: $($fixture.name)"
        throw "Scope Z5C blocked: fixture could not be resolved: $($fixture.name)"
    }
    $fixtureDir = Join-Path $workRootFull (Join-Path $fixture.name 'fixture')
    New-Item -ItemType Directory -Path $fixtureDir -Force | Out-Null
    $samplePqr = Join-Path $fixtureDir ("{0}_sampled.pqr" -f $fixture.name)
    $sampleInfo = Write-FmmCliSampledPqr -InputPath $fixturePath `
        -OutputPath $samplePqr -Cap $fixture.atom_cap
    $sampleInfoByFixture[$fixture.name] = $sampleInfo
    for ($repeat = 1; $repeat -le $fixture.repeats; ++$repeat) {
        foreach ($case in $cases) {
            $record = Invoke-Z5cCase -Fixture $fixture -Case $case `
                -Repeat $repeat -SamplePqr $samplePqr
            $records.Add($record) | Out-Null
            if ((-not $fixture.gated) -and $record.run.timed_out) {
                $nonGatedTimeout = $true
            }
        }
    }
}

$baseRows = foreach ($record in $records) {
    $sampleInfo = $sampleInfoByFixture[$record.fixture.name]
    $sampledAtoms = if ($null -ne $sampleInfo) { $sampleInfo.sampled_atoms } else { 0 }
    [pscustomobject]@{
        fixture=$record.fixture.name
        atom_cap=$record.fixture.atom_cap
        repeat=$record.repeat
        gated=$record.fixture.gated
        sampled_atoms=$sampledAtoms
        preconditioner_label=$record.case.label
        preconditioner_requested=$record.case.precond
        preconditioner_effective=$record.preconditioner.effective
        status=Get-Z5cStatus $record
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

$ratioByFixtureRepeat = @{}
foreach ($group in ($baseRows | Group-Object fixture,repeat)) {
    $block = @($group.Group | Where-Object {
            $_.preconditioner_label -eq 'block_diagonal'
        }) | Select-Object -First 1
    $ilu = @($group.Group | Where-Object {
            $_.preconditioner_label -eq 'sparse_ilu0'
        }) | Select-Object -First 1
    $ratio = [double]::NaN
    if ($null -ne $block -and $null -ne $ilu -and
        (Test-FmmCliFinite ([double]$block.true_residual)) -and
        (Test-FmmCliFinite ([double]$ilu.true_residual)) -and
        ([double]$ilu.true_residual -gt 0.0)) {
        $ratio = [double]$block.true_residual / [double]$ilu.true_residual
    }
    $ratioByFixtureRepeat[$group.Name] = $ratio
}

$statsByFixturePrecond = @{}
foreach ($group in ($baseRows | Group-Object fixture,preconditioner_label)) {
    $residuals = [double[]]@($group.Group | ForEach-Object {
            [double]$_.true_residual
        })
    $energies = [double[]]@($group.Group | ForEach-Object {
            [double]$_.energy_kj_mol
        })
    $energyMin = Get-Z5cMin $energies
    $energyMax = Get-Z5cMax $energies
    $statsByFixturePrecond[$group.Name] = [pscustomobject]@{
        residual_median=Get-Z5cMedian $residuals
        residual_min=Get-Z5cMin $residuals
        residual_max=Get-Z5cMax $residuals
        energy_spread=if ((Test-FmmCliFinite $energyMin) -and
            (Test-FmmCliFinite $energyMax)) { $energyMax - $energyMin } else { [double]::NaN }
        audit_max=Get-Z5cMax ([double[]]@($group.Group | ForEach-Object {
                    [double]$_.audit_rel
                }))
        peak_memory_max=Get-Z5cMax ([double[]]@($group.Group | ForEach-Object {
                    [double]$_.peak_working_set_mb
                }))
    }
}

$rows = foreach ($row in $baseRows) {
    $role = switch ($row.preconditioner_label) {
        'block_diagonal' { 'required_block_baseline' }
        'sparse_ilu0' { 'required_sparse_ilu0' }
        default { 'context' }
    }
    $ratioKey = "$($row.fixture), $($row.repeat)"
    $statsKey = "$($row.fixture), $($row.preconditioner_label)"
    $ratio = if ($ratioByFixtureRepeat.ContainsKey($ratioKey)) {
        $ratioByFixtureRepeat[$ratioKey]
    } else {
        [double]::NaN
    }
    $stats = if ($statsByFixturePrecond.ContainsKey($statsKey)) {
        $statsByFixturePrecond[$statsKey]
    } else {
        [pscustomobject]@{
            residual_median=[double]::NaN; residual_min=[double]::NaN
            residual_max=[double]::NaN; energy_spread=[double]::NaN
            audit_max=[double]::NaN; peak_memory_max=[double]::NaN
        }
    }
    [pscustomobject]@{
        fixture=$row.fixture
        atom_cap=$row.atom_cap
        repeat=$row.repeat
        gated=$row.gated
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
        residual_median=$stats.residual_median
        residual_min=$stats.residual_min
        residual_max=$stats.residual_max
        energy_spread=$stats.energy_spread
        audit_max=$stats.audit_max
        peak_memory_max=$stats.peak_memory_max
        ilu_vs_block_residual_ratio=$ratio
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
            repeat=$row.repeat
            preconditioner_label=$row.preconditioner_label
            iteration=$point.iteration
            residual=$point.residual
            log10_residual=$logValue
        }) | Out-Null
    }
}
$traceRows | Export-Csv -LiteralPath $TraceCsv -NoTypeInformation

function Write-Z5cSvgPlot([object[]]$TraceRows, [string]$Path) {
    $finite = @($TraceRows | Where-Object {
        Test-FmmCliFinite ([double]$_.log10_residual)
    })
    if ($finite.Count -eq 0) {
        Set-Content -LiteralPath $Path -Encoding ASCII -Value '<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="620"><text x="20" y="40">No finite residual trace</text></svg>'
        return
    }
    $width = 1000; $height = 620
    $left = 80; $right = 280; $top = 30; $bottom = 75
    $maxIter = [double](@($finite | Measure-Object -Property iteration -Maximum)[0].Maximum)
    $minLog = [double](@($finite | Measure-Object -Property log10_residual -Minimum)[0].Minimum)
    $maxLog = [double](@($finite | Measure-Object -Property log10_residual -Maximum)[0].Maximum)
    if ($maxIter -le 0.0) { $maxIter = 1.0 }
    if ([Math]::Abs($maxLog - $minLog) -lt 1.0e-12) {
        $maxLog += 1.0
        $minLog -= 1.0
    }
    $colors = @{
        block_diagonal = '#d55e00'
        sparse_ilu0 = '#009e73'
        diagonal = '#666666'
        overlap_schwarz = '#0072b2'
    }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add("<svg xmlns=`"http://www.w3.org/2000/svg`" width=`"$width`" height=`"$height`" viewBox=`"0 0 $width $height`">")
    $svg.Add('<rect width="100%" height="100%" fill="white"/>')
    $svg.Add("<line x1=`"$left`" y1=`"$($height-$bottom)`" x2=`"$($width-$right)`" y2=`"$($height-$bottom)`" stroke=`"black`"/>")
    $svg.Add("<line x1=`"$left`" y1=`"$top`" x2=`"$left`" y2=`"$($height-$bottom)`" stroke=`"black`"/>")
    $series = @($finite | Group-Object fixture,repeat,preconditioner_label)
    $legendY = 58
    foreach ($seriesGroup in $series) {
        $points = @($seriesGroup.Group | Sort-Object {[int]$_.iteration})
        if ($points.Count -eq 0) { continue }
        $label = [string]$points[0].preconditioner_label
        $fixtureLabel = [string]$points[0].fixture
        $repeatLabel = [int]$points[0].repeat
        $coords = foreach ($p in $points) {
            $x = $left + (($width - $left - $right) *
                ([double]$p.iteration / $maxIter))
            $y = $top + (($height - $top - $bottom) *
                (($maxLog - [double]$p.log10_residual) / ($maxLog - $minLog)))
            ('{0:F2},{1:F2}' -f $x,$y)
        }
        $color = if ($colors.ContainsKey($label)) { $colors[$label] } else { '#333333' }
        $opacity = if ($repeatLabel -gt 1) { '0.55' } else { '1.0' }
        $svg.Add("<polyline fill=`"none`" stroke=`"$color`" stroke-opacity=`"$opacity`" stroke-width=`"2.2`" points=`"$($coords -join ' ')`"/>")
        if ($legendY -lt ($height - 25)) {
            $legendText = "$fixtureLabel r$repeatLabel $label"
            $svg.Add("<line x1=`"735`" y1=`"$legendY`" x2=`"765`" y2=`"$legendY`" stroke=`"$color`" stroke-opacity=`"$opacity`" stroke-width=`"3`"/>")
            $svg.Add("<text x=`"773`" y=`"$($legendY+5)`" font-family=`"Arial`" font-size=`"12`">$legendText</text>")
            $legendY += 18
        }
    }
    $svg.Add('<text x="300" y="598" font-family="Arial" font-size="16">GMRES residual-history index</text>')
    $svg.Add('<text x="20" y="320" font-family="Arial" font-size="16" transform="rotate(-90 20,320)">log10 residual</text>')
    $svg.Add('<text x="80" y="24" font-family="Arial" font-size="18" font-weight="bold">Scope Z5C multi-fixture preconditioner convergence</text>')
    $svg.Add('</svg>')
    Set-Content -LiteralPath $Path -Encoding ASCII -Value $svg
}

Write-Z5cSvgPlot -TraceRows $traceRows.ToArray() -Path $PlotPath

function Test-Z5cFixtureContrast([object[]]$Rows, [string]$FixtureName) {
    $fixtureRows = @($Rows | Where-Object {
            $_.fixture -eq $FixtureName -and -not [bool]$_.gated
        })
    if ($fixtureRows.Count -eq 0) { return $false }
    foreach ($repeatGroup in ($fixtureRows | Group-Object repeat)) {
        $block = @($repeatGroup.Group | Where-Object {
                $_.preconditioner_label -eq 'block_diagonal'
            }) | Select-Object -First 1
        $ilu = @($repeatGroup.Group | Where-Object {
                $_.preconditioner_label -eq 'sparse_ilu0'
            }) | Select-Object -First 1
        if ($null -eq $block -or $null -eq $ilu) { return $false }
        $iluStrict = ($ilu.status -eq 'ok') -and
            (Test-FmmCliFinite ([double]$ilu.true_residual)) -and
            ([double]$ilu.true_residual -le 1.0e-6) -and
            (Test-FmmCliFinite ([double]$ilu.audit_rel)) -and
            ([double]$ilu.audit_rel -le 1.0e-2)
        if (-not $iluStrict) { return $false }
        $blockNonStrict = (-not (($block.status -eq 'ok') -and
            (Test-FmmCliFinite ([double]$block.true_residual)) -and
            ([double]$block.true_residual -le 1.0e-6)))
        $ratio = [double]$ilu.ilu_vs_block_residual_ratio
        $ratioStrong = (Test-FmmCliFinite $ratio) -and ($ratio -ge 100.0)
        if (-not ($blockNonStrict -or $ratioStrong)) { return $false }
    }
    return $true
}

$hasReadError = @($rows | Where-Object { $_.status -eq 'read_error' }).Count -gt 0
$hasPrimaryTimeout = @($rows | Where-Object {
        (-not [bool]$_.gated) -and $_.status -eq 'timeout'
    }).Count -gt 0
$barnasePass = Test-Z5cFixtureContrast -Rows $rows -FixtureName 'barnase'
$mediumPass = Test-Z5cFixtureContrast -Rows $rows -FixtureName 'research-medium'
$altPass = Test-Z5cFixtureContrast -Rows $rows -FixtureName 'research-alt'
$nonBarnasePassCount = @($mediumPass, $altPass | Where-Object { $_ }).Count

$z5cDecision = if ($hasReadError -or $traceRows.Count -eq 0) {
    'blocked_profile_schema'
} elseif ($Quick -and $barnasePass) {
    'money_plot_harness_ready'
} elseif ($Quick) {
    'harness_ready_but_no_contrast'
} elseif ($hasPrimaryTimeout -and (-not $barnasePass)) {
    'blocked_runtime'
} elseif ($barnasePass -and $mediumPass -and $altPass) {
    'multi_fixture_money_plot_pass'
} elseif ($barnasePass -and $nonBarnasePassCount -gt 0) {
    'fixture_specific_contrast'
} elseif ($barnasePass) {
    'barnase_only_money_plot'
} elseif ($hasPrimaryTimeout) {
    'blocked_runtime'
} else {
    'blocked_profile_schema'
}

[pscustomobject]@{
    scope_z5c_decision=$z5cDecision
    fixtures=($fixtureSpecs | ForEach-Object { $_.name }) -join ','
    rows=$rows.Count
    trace_rows=$traceRows.Count
    barnase_pass=$barnasePass
    research_medium_pass=$mediumPass
    research_alt_pass=$altPass
    matrix=$MatrixCsv
    trace=$TraceCsv
    plot=$PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$lines = New-Object System.Collections.Generic.List[string]
$lines.Add('TITAN FMM/M2L Scope Z5C Multi-Fixture Money Plot Report')
$lines.Add("Decision: $z5cDecision")
$lines.Add('Question: does the Barnase Z5B Block-Diagonal collapse vs sparse ILU(0) generalize across sampled protein fixtures?')
$lines.Add('')
$lines.Add('Rows:')
foreach ($row in $rows) {
    $lines.Add(('{0} r{1} {2}: status={3} effective={4} iter={5} residual={6:G6} ratio={7:G6} audit={8:G6} wall={9:G5}s peak={10:G5}MB' -f
        $row.fixture, [int]$row.repeat, $row.preconditioner_label,
        $row.status, $row.preconditioner_effective, [int]$row.iterations,
        [double]$row.true_residual, [double]$row.ilu_vs_block_residual_ratio,
        [double]$row.audit_rel, [double]$row.wall_seconds,
        [double]$row.peak_working_set_mb))
}
$lines.Add('')
$lines.Add('Interpretation:')
$lines.Add('- Z5B Barnase cap-256 money plot remains established independently.')
$lines.Add('- Z5C is a generalization gate; non-Barnase failures narrow the paper claim rather than invalidating Barnase.')
$lines.Add('- Direct verification is not rerun for every Z5C row; sampled audit plus Scope X.1 direct evidence remain the verification source.')
$lines.Add('')
$lines.Add("Matrix: $MatrixCsv")
$lines.Add("Trace: $TraceCsv")
$lines.Add("Plot: $PlotPath")
Set-Content -LiteralPath $ReportPath -Encoding ASCII -Value $lines
Set-Content -LiteralPath $SummaryPath -Encoding ASCII -Value @(
    "scope_z5c_decision=$z5cDecision",
    "rows=$($rows.Count)",
    "trace_rows=$($traceRows.Count)",
    "barnase_pass=$barnasePass",
    "research_medium_pass=$mediumPass",
    "research_alt_pass=$altPass",
    "matrix=$MatrixCsv",
    "trace=$TraceCsv",
    "plot=$PlotPath",
    "report=$ReportPath"
)

Write-Host "Scope Z5C report: $ReportPath"
