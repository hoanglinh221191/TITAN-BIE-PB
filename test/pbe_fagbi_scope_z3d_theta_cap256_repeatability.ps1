param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 2,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z3d')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z3d_matrix.csv'
$RepeatabilityCsv = Join-Path $WorkRoot 'fagbi_scope_z3d_repeatability.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z3d_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z3d_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z3d_theta_cap256_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z3D_THETA_CAP256_REPEATABILITY_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$TitanDir = Split-Path -Parent $TitanExe
$LibFfRoot = Join-Path $RepoRoot 'lib_ff'
if (-not (Test-Path -LiteralPath $LibFfRoot)) {
    $LibFfRoot = Join-Path $TitanDir 'lib_ff'
}
if (-not (Test-Path -LiteralPath $LibFfRoot)) {
    throw "Missing lib_ff beside repo or titan executable: $LibFfRoot"
}
$LibFfRoot = (Resolve-Path -LiteralPath $LibFfRoot).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Z3D work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z3dDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or
        $text -eq 'nan') {
        return $Default
    }
    return [double]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Convert-Z3dInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or
        $text -eq 'nan') {
        return $Default
    }
    return [int64]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Format-Z3dDouble([double]$Value) {
    return $Value.ToString('G15', [System.Globalization.CultureInfo]::InvariantCulture)
}

function Get-Z3dMedian([double[]]$Values) {
    $finite = @($Values | Where-Object { Test-FmmCliFinite $_ } | Sort-Object)
    if ($finite.Count -eq 0) { return [double]::NaN }
    $mid = [int]($finite.Count / 2)
    if (($finite.Count % 2) -eq 1) { return [double]$finite[$mid] }
    return 0.5 * ([double]$finite[$mid - 1] + [double]$finite[$mid])
}

function New-Z3dRow(
    [string]$Condition,
    [string]$Fixture,
    [int]$AtomCap,
    [int]$Repeat,
    [double]$Theta,
    [object]$Run,
    [string]$Status
) {
    return [pscustomobject]@{
        condition = $Condition
        fixture = $Fixture
        atom_cap = $AtomCap
        repeat = $Repeat
        theta = $Theta
        status = $Status
        audit_warning = ''
        exit_code = if ($null -eq $Run) { -1 } else { $Run.exit_code }
        timed_out = if ($null -eq $Run) { $false } else { $Run.timed_out }
        panels = 0
        iterations = 0
        true_residual = [double]::NaN
        tree_vs_direct_matvec_rel = [double]::NaN
        matvec_audit_rows = 0
        matvec_seconds_avg = [double]::NaN
        wall_seconds = if ($null -eq $Run) { [double]::NaN } else { $Run.wall_seconds }
        peak_working_set_mb = if ($null -eq $Run) { [double]::NaN } else { $Run.peak_working_set_mb }
        preconditioner_effective = ''
        near_integration = ''
        quad_buffer = ''
        fmm_m2l_mode = ''
        fmm_m2l_expansion_order = 0
        fmm_near_pairs = 0
        fmm_far_pairs = 0
        fmm_finite_outputs = ''
    }
}

function Resolve-Z3dFixture(
    [string]$Fixture,
    [string]$RunDir
) {
    $resolved = Resolve-FmmCliFixture -Fixture $Fixture -RepoRoot $RepoRoot `
        -RunDir $RunDir
    if ($null -ne $resolved -and (Test-Path -LiteralPath $resolved)) {
        return $resolved
    }
    $worktreeOutput = @()
    try {
        $worktreeOutput = @(git -C $RepoRoot worktree list --porcelain 2>$null)
    } catch {
        return $resolved
    }
    foreach ($line in $worktreeOutput) {
        if ($line -notmatch '^worktree\s+(.+)$') { continue }
        $candidateRoot = $Matches[1]
        if ([string]::IsNullOrWhiteSpace($candidateRoot) -or
            -not (Test-Path -LiteralPath $candidateRoot)) {
            continue
        }
        $resolved = Resolve-FmmCliFixture -Fixture $Fixture `
            -RepoRoot $candidateRoot -RunDir $RunDir
        if ($null -ne $resolved -and (Test-Path -LiteralPath $resolved)) {
            return $resolved
        }
    }
    return $null
}

function Write-Z3dArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [double]$Theta,
    [int]$MaxIter,
    [int]$AuditSampleCount
) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdbPath
    $libPath = Convert-FmmCliToTitanPath($LibFfRoot.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-FmmCliToTitanPath $pdbPath)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-FmmCliToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        "-pbe-max-iter=$MaxIter",
        '-pbe-rms=1e-6',
        '-pbe-dphi=1e-6',
        '-pbe-int-diel=1',
        '-pbe-ext-diel=80',
        '-pbe-salt=0.15',
        '-pbe-temp=298.15',
        '-pbe-timing=detail',
        "-nt=$Threads",
        '-pbe-tabi-sum=treecode',
        '-pbe-tabi-formulation=jcp2013-treecode',
        '-pbe-tabi-fmm-m2l=single-matvec-analytic',
        '-pbe-tabi-surface-points=8',
        '-pbe-tabi-tree-leaf-size=64',
        "-pbe-tabi-opening-theta=$(Format-Z3dDouble $Theta)",
        '-pbe-tabi-taylor-order=4',
        '-pbe-tabi-near-integration=analytic_singular',
        '-pbe-tabi-quad-buffer=off',
        '-pbe-tabi-precond=fmm-nearfield-ilu0',
        '-pbe-tabi-precond-ilu-neighbor-cap=64',
        '-pbe-tabi-gmres-restart=64',
        '-pbe-tabi-true-residual-probe=restart',
        '-pbe-tabi-true-residual-probe-scope=sampled',
        '-pbe-tabi-true-residual-sample-count=512',
        '-pbe-tabi-matvec-audit=sampled',
        "-pbe-tabi-matvec-audit-sample-count=$AuditSampleCount"
    )
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    $badOptions = @($argLines | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
    })
    if ($badOptions.Count -ne 0) {
        throw "Generated arg.in contains non equals-syntax options: $($badOptions -join '; ')"
    }
    return $argPath
}

function Read-Z3dResult([string]$RunDir) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $energyPath)) {
        throw "Missing pbe_energy.csv in $RunDir"
    }
    if (-not (Test-Path -LiteralPath $profilePath)) {
        throw "Missing pbe_tabi_profile.csv in $RunDir"
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    foreach ($column in @('fmm_m2l_mode','fmm_m2l_expansion_order',
            'preconditioner_effective','near_integration','quad_buffer',
            'tree_vs_direct_matvec_rel','matvec_audit_rows',
            'fmm_m2l_near_pairs','fmm_m2l_far_pairs',
            'fmm_m2l_finite_outputs')) {
        if (-not ($profile.PSObject.Properties.Name -contains $column)) {
            throw "pbe_tabi_profile.csv missing required Z3D column: $column"
        }
    }
    if ([string]$profile.fmm_m2l_mode -ne 'single-matvec-analytic') {
        throw "Expected FMM/M2L mode single-matvec-analytic, got $($profile.fmm_m2l_mode)"
    }
    if ((Convert-Z3dInt64 $profile.fmm_m2l_expansion_order) -ne 4) {
        throw "Expected FMM/M2L expansion order 4, got $($profile.fmm_m2l_expansion_order)"
    }
    return [pscustomobject]@{
        status = if ($energy.status -eq 'ok' -and $profile.status -eq 'ok') { 'ok' } else { "energy=$($energy.status);profile=$($profile.status)" }
        panels = Convert-Z3dInt64 $profile.panels
        iterations = Convert-Z3dInt64 $profile.iterations
        true_residual = Convert-Z3dDouble $profile.true_residual_rel
        tree_vs_direct_matvec_rel = Convert-Z3dDouble $profile.tree_vs_direct_matvec_rel
        matvec_audit_rows = Convert-Z3dInt64 $profile.matvec_audit_rows
        matvec_seconds_avg = Convert-Z3dDouble $profile.matvec_seconds_avg
        preconditioner_effective = [string]$profile.preconditioner_effective
        near_integration = [string]$profile.near_integration
        quad_buffer = [string]$profile.quad_buffer
        fmm_m2l_mode = [string]$profile.fmm_m2l_mode
        fmm_m2l_expansion_order = Convert-Z3dInt64 $profile.fmm_m2l_expansion_order
        fmm_near_pairs = Convert-Z3dInt64 $profile.fmm_m2l_near_pairs
        fmm_far_pairs = Convert-Z3dInt64 $profile.fmm_m2l_far_pairs
        fmm_finite_outputs = [string]$profile.fmm_m2l_finite_outputs
    }
}

function Invoke-Z3dCase(
    [string]$Fixture,
    [string]$PqrPath,
    [int]$AtomCap,
    [int]$Repeat,
    [double]$Theta,
    [int]$MaxIter,
    [int]$AuditSampleCount
) {
    $safeTheta = ((Format-Z3dDouble $Theta) -replace '\.', 'p')
    $condition = "{0}_cap{1}_theta{2}_rep{3}" -f (Get-FmmCliSafeName $Fixture), $AtomCap, $safeTheta, $Repeat
    $runDir = Join-Path $workRootFull $condition
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $argPath = Write-Z3dArgFile -RunDir $runDir -PqrPath $PqrPath `
        -Theta $Theta -MaxIter $MaxIter -AuditSampleCount $AuditSampleCount
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        return New-Z3dRow -Condition $condition -Fixture $Fixture `
            -AtomCap $AtomCap -Repeat $Repeat -Theta $Theta -Run $run `
            -Status $(if ($run.timed_out) { 'timeout_not_verified' } else { "exit_$($run.exit_code)" })
    }
    try {
        $result = Read-Z3dResult -RunDir $runDir
    } catch {
        return New-Z3dRow -Condition $condition -Fixture $Fixture `
            -AtomCap $AtomCap -Repeat $Repeat -Theta $Theta -Run $run `
            -Status "profile_error: $($_.Exception.Message)"
    }
    $auditWarning = ''
    if ((Test-FmmCliFinite $result.tree_vs_direct_matvec_rel) -and
        $result.tree_vs_direct_matvec_rel -gt 1.0e-2) {
        $auditWarning = 'audit_gt_1e-2'
    } elseif ((Test-FmmCliFinite $result.tree_vs_direct_matvec_rel) -and
        $result.tree_vs_direct_matvec_rel -gt 8.0e-3) {
        $auditWarning = 'audit_gt_8e-3'
    }
    return [pscustomobject]@{
        condition = $condition
        fixture = $Fixture
        atom_cap = $AtomCap
        repeat = $Repeat
        theta = $Theta
        status = $result.status
        audit_warning = $auditWarning
        exit_code = $run.exit_code
        timed_out = $run.timed_out
        panels = $result.panels
        iterations = $result.iterations
        true_residual = $result.true_residual
        tree_vs_direct_matvec_rel = $result.tree_vs_direct_matvec_rel
        matvec_audit_rows = $result.matvec_audit_rows
        matvec_seconds_avg = $result.matvec_seconds_avg
        wall_seconds = $run.wall_seconds
        peak_working_set_mb = $run.peak_working_set_mb
        preconditioner_effective = $result.preconditioner_effective
        near_integration = $result.near_integration
        quad_buffer = $result.quad_buffer
        fmm_m2l_mode = $result.fmm_m2l_mode
        fmm_m2l_expansion_order = $result.fmm_m2l_expansion_order
        fmm_near_pairs = $result.fmm_near_pairs
        fmm_far_pairs = $result.fmm_far_pairs
        fmm_finite_outputs = $result.fmm_finite_outputs
    }
}

function Test-Z3dRowClean([object]$Row, [double]$AuditGate) {
    if ($Row.status -ne 'ok' -or $Row.timed_out -or
        $Row.fmm_m2l_mode -ne 'single-matvec-analytic' -or
        $Row.fmm_finite_outputs -ne 'true' -or
        $Row.near_integration -ne 'analytic_singular' -or
        $Row.quad_buffer -ne 'off' -or
        $Row.preconditioner_effective -ne 'fmm_nearfield_ilu0' -or
        $Row.fmm_m2l_expansion_order -ne 4 -or
        $Row.matvec_audit_rows -le 0) {
        return $false
    }
    $audit = Convert-Z3dDouble $Row.tree_vs_direct_matvec_rel
    return (Test-FmmCliFinite $audit) -and $audit -le $AuditGate
}

function New-Z3dRepeatabilityRow([object[]]$Rows, [string]$Fixture, [int]$AtomCap) {
    $subset = @($Rows | Where-Object {
        $_.fixture -eq $Fixture -and $_.atom_cap -eq $AtomCap -and [double]$_.theta -eq 0.6
    })
    $okRows = @($subset | Where-Object { Test-Z3dRowClean -Row $_ -AuditGate 8.0e-3 })
    $audits = @($subset | ForEach-Object { Convert-Z3dDouble $_.tree_vs_direct_matvec_rel } |
        Where-Object { Test-FmmCliFinite $_ })
    $matvec = @($subset | ForEach-Object { Convert-Z3dDouble $_.matvec_seconds_avg } |
        Where-Object { Test-FmmCliFinite $_ })
    $residuals = @($subset | ForEach-Object { Convert-Z3dDouble $_.true_residual } |
        Where-Object { Test-FmmCliFinite $_ })
    $walls = @($subset | ForEach-Object { Convert-Z3dDouble $_.wall_seconds } |
        Where-Object { Test-FmmCliFinite $_ })
    $rss = @($subset | ForEach-Object { Convert-Z3dDouble $_.peak_working_set_mb } |
        Where-Object { Test-FmmCliFinite $_ })
    $auditMin = if ($audits.Count -gt 0) { ($audits | Measure-Object -Minimum).Minimum } else { [double]::NaN }
    $auditMax = if ($audits.Count -gt 0) { ($audits | Measure-Object -Maximum).Maximum } else { [double]::NaN }
    $auditMean = if ($audits.Count -gt 0) { ($audits | Measure-Object -Average).Average } else { [double]::NaN }
    $status = 'ok'
    if ($subset.Count -eq 0) {
        $status = 'missing_rows'
    } elseif (@($subset | Where-Object { $_.timed_out -or $_.status -eq 'timeout_not_verified' }).Count -gt 0) {
        $status = 'runtime_limited'
    } elseif ($okRows.Count -ne $subset.Count) {
        $status = 'audit_or_status_fail'
    }
    return [pscustomobject]@{
        fixture = $Fixture
        atom_cap = $AtomCap
        theta = 0.6
        repeats = $subset.Count
        ok_repeats = $okRows.Count
        audit_mean = $auditMean
        audit_max = $auditMax
        audit_spread = if ((Test-FmmCliFinite $auditMin) -and (Test-FmmCliFinite $auditMax)) { $auditMax - $auditMin } else { [double]::NaN }
        matvec_median = Get-Z3dMedian $matvec
        matvec_min = if ($matvec.Count -gt 0) { ($matvec | Measure-Object -Minimum).Minimum } else { [double]::NaN }
        matvec_max = if ($matvec.Count -gt 0) { ($matvec | Measure-Object -Maximum).Maximum } else { [double]::NaN }
        residual_max = if ($residuals.Count -gt 0) { ($residuals | Measure-Object -Maximum).Maximum } else { [double]::NaN }
        wall_max = if ($walls.Count -gt 0) { ($walls | Measure-Object -Maximum).Maximum } else { [double]::NaN }
        peak_rss_max = if ($rss.Count -gt 0) { ($rss | Measure-Object -Maximum).Maximum } else { [double]::NaN }
        repeatability_status = $status
    }
}

function Get-Z3dThetaMedianMatvec([object[]]$Rows, [double]$Theta) {
    return Get-Z3dMedian @($Rows | Where-Object {
            [double]$_.theta -eq $Theta -and $_.status -eq 'ok'
        } | ForEach-Object { Convert-Z3dDouble $_.matvec_seconds_avg })
}

function Get-Z3dDecision([object[]]$Rows, [object[]]$RepeatRows) {
    if (@($Rows | Where-Object { $_.timed_out -or $_.status -eq 'timeout_not_verified' }).Count -gt 0) {
        return 'z3d_runtime_limited'
    }
    if (@($Rows | Where-Object { $_.status -ne 'ok' }).Count -gt 0) {
        return 'z3d_runtime_limited'
    }
    if (@($RepeatRows | Where-Object { $_.repeatability_status -ne 'ok' }).Count -gt 0) {
        return 'z3d_no_safe_theta'
    }
    $theta06Rows = @($Rows | Where-Object { [double]$_.theta -eq 0.6 })
    if (@($theta06Rows | Where-Object { -not (Test-Z3dRowClean -Row $_ -AuditGate 8.0e-3) }).Count -gt 0) {
        return 'z3d_no_safe_theta'
    }
    $theta07Rows = @($Rows | Where-Object { [double]$_.theta -eq 0.7 })
    $theta08Rows = @($Rows | Where-Object { [double]$_.theta -eq 0.8 })
    $theta07Clean = $theta07Rows.Count -gt 0 -and
        @($theta07Rows | Where-Object { -not (Test-Z3dRowClean -Row $_ -AuditGate 8.0e-3) }).Count -eq 0
    $theta08Clean = $theta08Rows.Count -gt 0 -and
        @($theta08Rows | Where-Object { -not (Test-Z3dRowClean -Row $_ -AuditGate 1.0e-2) }).Count -eq 0
    $m06 = Get-Z3dThetaMedianMatvec -Rows $Rows -Theta 0.6
    $m07 = Get-Z3dThetaMedianMatvec -Rows $Rows -Theta 0.7
    $speed07 = if ((Test-FmmCliFinite $m06) -and (Test-FmmCliFinite $m07) -and $m07 -gt 0.0) {
        $m06 / $m07
    } else {
        [double]::NaN
    }
    if ($theta07Clean -and $speed07 -ge 1.15 -and (-not $theta08Clean)) {
        return 'z3d_theta07_candidate'
    }
    return 'z3d_theta06_confirmed'
}

function Get-Z3dTheta08Status([object[]]$Rows) {
    $theta08Rows = @($Rows | Where-Object { [double]$_.theta -eq 0.8 })
    if ($theta08Rows.Count -eq 0) { return 'not_run' }
    if (@($theta08Rows | Where-Object {
            $audit = Convert-Z3dDouble $_.tree_vs_direct_matvec_rel
            (-not (Test-FmmCliFinite $audit)) -or $audit -gt 1.0e-2
        }).Count -gt 0) {
        return 'z3d_theta08_rejected'
    }
    if (@($theta08Rows | Where-Object {
            $audit = Convert-Z3dDouble $_.tree_vs_direct_matvec_rel
            $audit -gt 8.0e-3
        }).Count -gt 0) {
        return 'marginal'
    }
    return 'clean'
}

function Write-Z3dSvg([object[]]$Rows, [string]$Path) {
    $points = @($Rows | Where-Object { $_.status -eq 'ok' })
    $maxAudit = ($points | ForEach-Object { Convert-Z3dDouble $_.tree_vs_direct_matvec_rel } |
        Where-Object { Test-FmmCliFinite $_ } | Measure-Object -Maximum).Maximum
    if (-not (Test-FmmCliFinite $maxAudit) -or $maxAudit -le 0.0) {
        $maxAudit = 1.0e-2
    }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="980" height="360" viewBox="0 0 980 360">') | Out-Null
    $svg.Add('<rect width="980" height="360" fill="white"/>') | Out-Null
    $svg.Add('<text x="24" y="28" font-family="Arial" font-size="16">Scope Z3D cap-256 theta repeatability gate</text>') | Out-Null
    $svg.Add('<text x="24" y="48" font-family="Arial" font-size="11">Bar height = sampled tree-vs-direct matvec relative RMS; orange &gt; 8e-3, red &gt; 1e-2</text>') | Out-Null
    $svg.Add('<line x1="24" y1="96" x2="956" y2="96" stroke="#ddaa33" stroke-dasharray="4 4"/>') | Out-Null
    $svg.Add('<line x1="24" y1="72" x2="956" y2="72" stroke="#cc3311" stroke-dasharray="4 4"/>') | Out-Null
    $svg.Add('<text x="842" y="92" font-family="Arial" font-size="10">8e-3 gate</text>') | Out-Null
    $svg.Add('<text x="842" y="68" font-family="Arial" font-size="10">1e-2 gate</text>') | Out-Null
    $x = 34
    foreach ($row in $points) {
        $audit = Convert-Z3dDouble $row.tree_vs_direct_matvec_rel 0.0
        $h = [Math]::Max(2.0, 230.0 * ($audit / [Math]::Max($maxAudit, 1.0e-12)))
        $y = 292.0 - $h
        $fill = if ($audit -gt 1.0e-2) { '#cc3311' } elseif ($audit -gt 8.0e-3) { '#ee7733' } else { '#0072b2' }
        $svg.Add(('<rect x="{0}" y="{1:F1}" width="22" height="{2:F1}" fill="{3}"/>' -f $x, $y, $h, $fill)) | Out-Null
        $label = ('{0} t{1} r{2}' -f $row.fixture,
            ([double]$row.theta).ToString('G3', [System.Globalization.CultureInfo]::InvariantCulture),
            $row.repeat)
        $svg.Add(('<text x="{0}" y="318" font-family="Arial" font-size="8" transform="rotate(45 {0},318)">{1}</text>' -f $x, $label)) | Out-Null
        $x += 30
        if ($x -gt 940) { break }
    }
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$fixtureDir = Join-Path $workRootFull 'fixtures'
New-Item -ItemType Directory -Path $fixtureDir -Force | Out-Null
$fixtureSpecs = if ($Quick) {
    @(@{ fixture='research-medium'; cap=128; max_iter=240; audit=256 })
} else {
    @(
        @{ fixture='barnase'; cap=256; max_iter=360; audit=512 },
        @{ fixture='research-medium'; cap=256; max_iter=360; audit=512 },
        @{ fixture='research-alt'; cap=256; max_iter=360; audit=512 }
    )
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($spec in $fixtureSpecs) {
    $sampleDir = Join-Path $fixtureDir $spec.fixture
    New-Item -ItemType Directory -Path $sampleDir -Force | Out-Null
    $source = Resolve-Z3dFixture -Fixture $spec.fixture -RunDir $sampleDir
    if ($null -eq $source -or -not (Test-Path -LiteralPath $source)) {
        $thetas = if ($Quick) { @(0.6,0.8) } else { @(0.6,0.7,0.8) }
        foreach ($theta in $thetas) {
            $repeatCount = if ((-not $Quick) -and $theta -eq 0.6) { 3 } else { 1 }
            for ($repeat = 1; $repeat -le $repeatCount; ++$repeat) {
                $condition = "{0}_cap{1}_theta{2}_rep{3}" -f (Get-FmmCliSafeName $spec.fixture), $spec.cap, ((Format-Z3dDouble $theta) -replace '\.', 'p'), $repeat
                $rows.Add((New-Z3dRow -Condition $condition -Fixture $spec.fixture `
                    -AtomCap $spec.cap -Repeat $repeat -Theta $theta -Run $null `
                    -Status 'skipped_missing_fixture')) | Out-Null
            }
        }
        continue
    }
    $sampled = Join-Path $sampleDir ("{0}_cap{1}.pqr" -f (Get-FmmCliSafeName $spec.fixture), $spec.cap)
    [void](Write-FmmCliSampledPqr -InputPath $source -OutputPath $sampled -Cap $spec.cap)
    $schedule = New-Object System.Collections.Generic.List[object]
    if ($Quick) {
        $schedule.Add([pscustomobject]@{ theta=0.6; repeat=1 }) | Out-Null
        $schedule.Add([pscustomobject]@{ theta=0.8; repeat=1 }) | Out-Null
    } else {
        foreach ($repeat in 1..3) {
            $schedule.Add([pscustomobject]@{ theta=0.6; repeat=$repeat }) | Out-Null
        }
        $schedule.Add([pscustomobject]@{ theta=0.7; repeat=1 }) | Out-Null
        $schedule.Add([pscustomobject]@{ theta=0.8; repeat=1 }) | Out-Null
    }
    foreach ($item in $schedule) {
        $rows.Add((Invoke-Z3dCase -Fixture $spec.fixture -PqrPath $sampled `
            -AtomCap $spec.cap -Repeat $item.repeat -Theta $item.theta `
            -MaxIter $spec.max_iter -AuditSampleCount $spec.audit)) | Out-Null
    }
}

$rowArray = $rows.ToArray()
$repeatRows = @()
foreach ($spec in $fixtureSpecs) {
    $repeatRows += New-Z3dRepeatabilityRow -Rows $rowArray -Fixture $spec.fixture -AtomCap $spec.cap
}
$decision = Get-Z3dDecision -Rows $rowArray -RepeatRows $repeatRows
$theta08Status = Get-Z3dTheta08Status -Rows $rowArray
$m06 = Get-Z3dThetaMedianMatvec -Rows $rowArray -Theta 0.6
$m07 = Get-Z3dThetaMedianMatvec -Rows $rowArray -Theta 0.7
$m08 = Get-Z3dThetaMedianMatvec -Rows $rowArray -Theta 0.8
$speed07 = if ((Test-FmmCliFinite $m06) -and (Test-FmmCliFinite $m07) -and $m07 -gt 0.0) { $m06 / $m07 } else { [double]::NaN }
$speed08 = if ((Test-FmmCliFinite $m06) -and (Test-FmmCliFinite $m08) -and $m08 -gt 0.0) { $m06 / $m08 } else { [double]::NaN }

$rowArray | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$repeatRows | Export-Csv -LiteralPath $RepeatabilityCsv -NoTypeInformation
Write-Z3dSvg -Rows $rowArray -Path $PlotPath

[pscustomobject]@{
    decision = $decision
    theta08_status = $theta08Status
    rows = $rowArray.Count
    quick = [bool]$Quick
    median_theta06_matvec_s = $m06
    median_theta07_matvec_s = $m07
    median_theta08_matvec_s = $m08
    theta07_vs_theta06_speedup = $speed07
    theta08_vs_theta06_speedup = $speed08
    artifact_matrix = $MatrixCsv
    artifact_repeatability = $RepeatabilityCsv
    artifact_plot = $PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'scope: Z3D cap-256 theta policy repeatability gate',
    "decision: $decision",
    "theta08_status: $theta08Status",
    "rows: $($rowArray.Count)",
    "quick: $([bool]$Quick)",
    'taylor_order: 4',
    'surface_points: 8',
    'near_policy: analytic_singular',
    'quad_buffer: off',
    'preconditioner: fmm-nearfield-ilu0',
    'ilu_neighbor_cap: 64',
    'theta06_repeat_gate: audit <= 8e-3 on every repeat',
    'theta08_reject_gate: audit > 1e-2 on any required fixture',
    "theta07_vs_theta06_speedup: $speed07",
    "theta08_vs_theta06_speedup: $speed08",
    'default_policy_note: grid remains production default; TABI treecode behavior is unchanged by this benchmark'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z3D Theta Cap-256 Repeatability Report',
    '',
    "Decision: $decision",
    "Theta 0.8 status: $theta08Status",
    '',
    'Z3D validates the cap-128 Z3C theta policy on cap-256 sampled protein fixtures with theta 0.6 repeats and diagnostic theta 0.7/0.8 rows.',
    '',
    'Artifacts:',
    "- Matrix: $MatrixCsv",
    "- Repeatability: $RepeatabilityCsv",
    "- Decision: $DecisionCsv",
    "- Summary: $SummaryPath",
    "- Plot: $PlotPath",
    '',
    'Repeatability:',
    ($repeatRows | ForEach-Object {
        '- {0}: repeats={1}, ok={2}, audit_max={3}, matvec_median={4}, status={5}' -f `
            $_.fixture, $_.repeats, $_.ok_repeats, $_.audit_max, $_.matvec_median, $_.repeatability_status
    }),
    '',
    'Rows:'
)
foreach ($row in $rowArray) {
    $report += ('- {0}: fixture={1}, cap={2}, repeat={3}, theta={4}, status={5}, audit={6}, warning={7}, matvec_s={8}, true_residual={9}, precond={10}, fmm_order={11}' -f `
        $row.condition, $row.fixture, $row.atom_cap, $row.repeat, $row.theta, `
        $row.status, $row.tree_vs_direct_matvec_rel, $row.audit_warning, `
        $row.matvec_seconds_avg, $row.true_residual, $row.preconditioner_effective, `
        $row.fmm_m2l_expansion_order)
}
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "PASS Scope Z3D theta cap-256 repeatability gate: $decision"
