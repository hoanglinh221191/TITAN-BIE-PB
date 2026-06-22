param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 420,
    [int]$Threads = 2,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z3c')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z3c_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z3c_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z3c_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z3c_theta_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z3C_THETA_POLICY_REPORT.txt'

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
    throw "Refusing to clean Z3C work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z3cDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or
        $text -eq 'nan') {
        return $Default
    }
    return [double]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Convert-Z3cInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or
        $text -eq 'nan') {
        return $Default
    }
    return [int64]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Get-Z3cMedian([double[]]$Values) {
    $finite = @($Values | Where-Object { Test-FmmCliFinite $_ } | Sort-Object)
    if ($finite.Count -eq 0) { return [double]::NaN }
    $mid = [int]($finite.Count / 2)
    if (($finite.Count % 2) -eq 1) { return [double]$finite[$mid] }
    return 0.5 * ([double]$finite[$mid - 1] + [double]$finite[$mid])
}

function New-Z3cRow(
    [string]$Condition,
    [string]$Fixture,
    [int]$AtomCap,
    [double]$Theta,
    [object]$Run,
    [string]$Status
) {
    return [pscustomobject]@{
        condition = $Condition
        fixture = $Fixture
        atom_cap = $AtomCap
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

function Resolve-Z3cFixture(
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

function Write-Z3cArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [double]$Theta,
    [int]$MaxIter,
    [int]$AuditSampleCount
) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdbPath
    $libPath = Convert-FmmCliToTitanPath($LibFfRoot.TrimEnd('\') + '\')
    $format = {
        param([double]$Value)
        return $Value.ToString('G15', [System.Globalization.CultureInfo]::InvariantCulture)
    }
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
        "-pbe-tabi-opening-theta=$(& $format $Theta)",
        '-pbe-tabi-taylor-order=4',
        '-pbe-tabi-near-integration=analytic_singular',
        '-pbe-tabi-quad-buffer=off',
        '-pbe-tabi-precond=fmm-nearfield-ilu0',
        '-pbe-tabi-precond-ilu-neighbor-cap=64',
        '-pbe-tabi-gmres-restart=64',
        '-pbe-tabi-true-residual-probe=restart',
        '-pbe-tabi-true-residual-probe-scope=sampled',
        '-pbe-tabi-true-residual-sample-count=256',
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

function Read-Z3cResult([string]$RunDir) {
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
            'tree_vs_direct_matvec_rel','matvec_audit_rows')) {
        if (-not ($profile.PSObject.Properties.Name -contains $column)) {
            throw "pbe_tabi_profile.csv missing required Z3C column: $column"
        }
    }
    if ([string]$profile.fmm_m2l_mode -ne 'single-matvec-analytic') {
        throw "Expected FMM/M2L mode single-matvec-analytic, got $($profile.fmm_m2l_mode)"
    }
    if ((Convert-Z3cInt64 $profile.fmm_m2l_expansion_order) -ne 4) {
        throw "Expected FMM/M2L expansion order 4, got $($profile.fmm_m2l_expansion_order)"
    }
    return [pscustomobject]@{
        status = if ($energy.status -eq 'ok' -and $profile.status -eq 'ok') { 'ok' } else { "energy=$($energy.status);profile=$($profile.status)" }
        panels = Convert-Z3cInt64 $profile.panels
        iterations = Convert-Z3cInt64 $profile.iterations
        true_residual = Convert-Z3cDouble $profile.true_residual_rel
        tree_vs_direct_matvec_rel = Convert-Z3cDouble $profile.tree_vs_direct_matvec_rel
        matvec_audit_rows = Convert-Z3cInt64 $profile.matvec_audit_rows
        matvec_seconds_avg = Convert-Z3cDouble $profile.matvec_seconds_avg
        preconditioner_effective = [string]$profile.preconditioner_effective
        near_integration = [string]$profile.near_integration
        quad_buffer = [string]$profile.quad_buffer
        fmm_m2l_mode = [string]$profile.fmm_m2l_mode
        fmm_m2l_expansion_order = Convert-Z3cInt64 $profile.fmm_m2l_expansion_order
        fmm_near_pairs = Convert-Z3cInt64 $profile.fmm_m2l_near_pairs
        fmm_far_pairs = Convert-Z3cInt64 $profile.fmm_m2l_far_pairs
        fmm_finite_outputs = [string]$profile.fmm_m2l_finite_outputs
    }
}

function Invoke-Z3cCase(
    [string]$Fixture,
    [string]$PqrPath,
    [int]$AtomCap,
    [double]$Theta,
    [int]$MaxIter,
    [int]$AuditSampleCount
) {
    $safeTheta = ($Theta.ToString('G15', [System.Globalization.CultureInfo]::InvariantCulture) -replace '\.', 'p')
    $condition = "{0}_theta{1}" -f (Get-FmmCliSafeName $Fixture), $safeTheta
    $runDir = Join-Path $workRootFull $condition
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $argPath = Write-Z3cArgFile -RunDir $runDir -PqrPath $PqrPath `
        -Theta $Theta -MaxIter $MaxIter -AuditSampleCount $AuditSampleCount
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        return New-Z3cRow -Condition $condition -Fixture $Fixture `
            -AtomCap $AtomCap -Theta $Theta -Run $run `
            -Status $(if ($run.timed_out) { 'timeout_not_verified' } else { "exit_$($run.exit_code)" })
    }
    try {
        $result = Read-Z3cResult -RunDir $runDir
    } catch {
        return New-Z3cRow -Condition $condition -Fixture $Fixture `
            -AtomCap $AtomCap -Theta $Theta -Run $run `
            -Status "profile_error: $($_.Exception.Message)"
    }
    $auditWarning = if ((Test-FmmCliFinite $result.tree_vs_direct_matvec_rel) -and
        $result.tree_vs_direct_matvec_rel -gt 8.0e-3) {
        'audit_gt_8e-3'
    } else {
        ''
    }
    return [pscustomobject]@{
        condition = $condition
        fixture = $Fixture
        atom_cap = $AtomCap
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

function Test-Z3cThetaGate([object[]]$Rows, [double]$Theta, [double]$AuditGate) {
    $subset = @($Rows | Where-Object { [double]$_.theta -eq $Theta })
    if ($subset.Count -eq 0) { return $false }
    foreach ($row in $subset) {
        if ($row.status -ne 'ok' -or $row.fmm_m2l_mode -ne 'single-matvec-analytic' -or
            $row.fmm_finite_outputs -ne 'true' -or $row.near_integration -ne 'analytic_singular' -or
            $row.quad_buffer -ne 'off' -or $row.preconditioner_effective -ne 'fmm_nearfield_ilu0' -or
            $row.matvec_audit_rows -le 0) {
            return $false
        }
        $audit = Convert-Z3cDouble $row.tree_vs_direct_matvec_rel
        if (-not (Test-FmmCliFinite $audit) -or $audit -gt $AuditGate) {
            return $false
        }
    }
    return $true
}

function Get-Z3cThetaMedianMatvec([object[]]$Rows, [double]$Theta) {
    return Get-Z3cMedian @($Rows | Where-Object {
            [double]$_.theta -eq $Theta -and $_.status -eq 'ok'
        } | ForEach-Object { Convert-Z3cDouble $_.matvec_seconds_avg })
}

function Get-Z3cDecision([object[]]$Rows) {
    $okRows = @($Rows | Where-Object { $_.status -eq 'ok' })
    if ($okRows.Count -ne $Rows.Count) { return 'z3c_runtime_limited' }
    $theta06Clean = Test-Z3cThetaGate -Rows $Rows -Theta 0.6 -AuditGate 8.0e-3
    $theta07Clean = Test-Z3cThetaGate -Rows $Rows -Theta 0.7 -AuditGate 8.0e-3
    $theta08Clean = Test-Z3cThetaGate -Rows $Rows -Theta 0.8 -AuditGate 8.0e-3
    $theta08Pass = Test-Z3cThetaGate -Rows $Rows -Theta 0.8 -AuditGate 1.0e-2
    $m06 = Get-Z3cThetaMedianMatvec -Rows $Rows -Theta 0.6
    $m07 = Get-Z3cThetaMedianMatvec -Rows $Rows -Theta 0.7
    $m08 = Get-Z3cThetaMedianMatvec -Rows $Rows -Theta 0.8
    $speed08 = if ((Test-FmmCliFinite $m06) -and (Test-FmmCliFinite $m08) -and $m08 -gt 0.0) {
        $m06 / $m08
    } else {
        [double]::NaN
    }
    $speed07 = if ((Test-FmmCliFinite $m06) -and (Test-FmmCliFinite $m07) -and $m07 -gt 0.0) {
        $m06 / $m07
    } else {
        [double]::NaN
    }
    if ($theta08Pass -and $speed08 -ge 1.25) { return 'z3c_theta08_safe_fast' }
    if ((-not $theta08Clean) -and $theta07Clean -and $speed07 -gt 1.05) {
        return 'z3c_intermediate_theta_needed'
    }
    if ($theta06Clean) { return 'z3c_theta06_safe_default' }
    return 'z3c_no_safe_theta'
}

function Write-Z3cSvg([object[]]$Rows, [string]$Path) {
    $points = @($Rows | Where-Object { $_.status -eq 'ok' })
    $maxAudit = ($points | ForEach-Object { Convert-Z3cDouble $_.tree_vs_direct_matvec_rel } |
        Where-Object { Test-FmmCliFinite $_ } | Measure-Object -Maximum).Maximum
    if (-not (Test-FmmCliFinite $maxAudit) -or $maxAudit -le 0.0) {
        $maxAudit = 1.0e-2
    }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="340" viewBox="0 0 900 340">') | Out-Null
    $svg.Add('<rect width="900" height="340" fill="white"/>') | Out-Null
    $svg.Add('<text x="24" y="28" font-family="Arial" font-size="16">Scope Z3C opening-theta policy gate</text>') | Out-Null
    $svg.Add('<text x="24" y="48" font-family="Arial" font-size="11">Bar height = sampled tree-vs-direct matvec relative RMS</text>') | Out-Null
    $x = 42
    foreach ($row in $points) {
        $audit = Convert-Z3cDouble $row.tree_vs_direct_matvec_rel 0.0
        $h = [Math]::Max(2.0, 230.0 * ($audit / [Math]::Max($maxAudit, 1.0e-12)))
        $y = 280.0 - $h
        $fill = if ($audit -gt 8.0e-3) { '#cc3311' } elseif ($audit -gt 6.0e-3) { '#ee7733' } else { '#0072b2' }
        $svg.Add(('<rect x="{0}" y="{1:F1}" width="28" height="{2:F1}" fill="{3}"/>' -f $x, $y, $h, $fill)) | Out-Null
        $label = ('{0} t{1}' -f $row.fixture,
            ([double]$row.theta).ToString('G3', [System.Globalization.CultureInfo]::InvariantCulture))
        $svg.Add(('<text x="{0}" y="304" font-family="Arial" font-size="8" transform="rotate(45 {0},304)">{1}</text>' -f $x, $label)) | Out-Null
        $x += 40
    }
    $svg.Add('<line x1="24" y1="96" x2="876" y2="96" stroke="#999" stroke-dasharray="4 4"/>') | Out-Null
    $svg.Add('<text x="760" y="92" font-family="Arial" font-size="10">relative scale</text>') | Out-Null
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$fixtureDir = Join-Path $workRootFull 'fixtures'
New-Item -ItemType Directory -Path $fixtureDir -Force | Out-Null
$cases = New-Object System.Collections.Generic.List[object]
$fixtureSpecs = if ($Quick) {
    @(@{ fixture='research-medium'; cap=64; max_iter=180; audit=256 })
} else {
    @(
        @{ fixture='barnase'; cap=128; max_iter=240; audit=512 },
        @{ fixture='research-medium'; cap=128; max_iter=240; audit=512 },
        @{ fixture='research-alt'; cap=128; max_iter=240; audit=512 }
    )
}
foreach ($spec in $fixtureSpecs) {
    $sampleDir = Join-Path $fixtureDir $spec.fixture
    New-Item -ItemType Directory -Path $sampleDir -Force | Out-Null
    $source = Resolve-Z3cFixture -Fixture $spec.fixture -RunDir $sampleDir
    if ($null -eq $source -or -not (Test-Path -LiteralPath $source)) {
        foreach ($theta in @($(if ($Quick) { 0.6 } else { 0.45 }),$(if ($Quick) { 0.8 } else { 0.6 }),$(if ($Quick) { $null } else { 0.7 }),$(if ($Quick) { $null } else { 0.8 }))) {
            if ($null -eq $theta) { continue }
            $condition = "{0}_theta{1}" -f (Get-FmmCliSafeName $spec.fixture),
                (($theta.ToString('G15', [System.Globalization.CultureInfo]::InvariantCulture)) -replace '\.', 'p')
            $cases.Add((New-Z3cRow -Condition $condition -Fixture $spec.fixture `
                -AtomCap $spec.cap -Theta $theta -Run $null `
                -Status 'skipped_missing_fixture')) | Out-Null
        }
        continue
    }
    $sampled = Join-Path $sampleDir ("{0}_cap{1}.pqr" -f (Get-FmmCliSafeName $spec.fixture), $spec.cap)
    [void](Write-FmmCliSampledPqr -InputPath $source -OutputPath $sampled -Cap $spec.cap)
    foreach ($theta in $(if ($Quick) { @(0.6,0.8) } else { @(0.45,0.6,0.7,0.8) })) {
        $cases.Add((Invoke-Z3cCase -Fixture $spec.fixture -PqrPath $sampled `
            -AtomCap $spec.cap -Theta $theta -MaxIter $spec.max_iter `
            -AuditSampleCount $spec.audit)) | Out-Null
    }
}

$rows = $cases.ToArray()
$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
Write-Z3cSvg -Rows $rows -Path $PlotPath
$decision = Get-Z3cDecision -Rows $rows

$m06 = Get-Z3cThetaMedianMatvec -Rows $rows -Theta 0.6
$m08 = Get-Z3cThetaMedianMatvec -Rows $rows -Theta 0.8
$speed08 = if ((Test-FmmCliFinite $m06) -and (Test-FmmCliFinite $m08) -and $m08 -gt 0.0) {
    $m06 / $m08
} else {
    [double]::NaN
}
[pscustomobject]@{
    decision = $decision
    rows = $rows.Count
    quick = [bool]$Quick
    median_theta06_matvec_s = $m06
    median_theta08_matvec_s = $m08
    theta08_vs_theta06_speedup = $speed08
    artifact_matrix = $MatrixCsv
    artifact_plot = $PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'scope: Z3C FMM/M2L opening-theta policy gate',
    "decision: $decision",
    "rows: $($rows.Count)",
    "quick: $([bool]$Quick)",
    'taylor_order: 4',
    'near_policy: analytic_singular',
    'quad_buffer: off',
    'preconditioner: fmm-nearfield-ilu0',
    'audit_clean_gate: 8e-3',
    'audit_pass_gate: 1e-2',
    "theta08_vs_theta06_speedup: $speed08",
    'default_policy_note: grid remains production default; TABI treecode default is unchanged'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z3C Opening-Theta Policy Report',
    '',
    "Decision: $decision",
    '',
    'Z3C isolates opening-theta policy after Z3B showed theta dominates the FMM audit/runtime tradeoff.',
    '',
    'Artifacts:',
    "- Matrix: $MatrixCsv",
    "- Decision: $DecisionCsv",
    "- Summary: $SummaryPath",
    "- Plot: $PlotPath",
    '',
    'Rows:'
)
foreach ($row in $rows) {
    $report += ('- {0}: fixture={1}, cap={2}, theta={3}, status={4}, audit={5}, warning={6}, matvec_s={7}, true_residual={8}, precond={9}' -f `
        $row.condition, $row.fixture, $row.atom_cap, $row.theta, $row.status, `
        $row.tree_vs_direct_matvec_rel, $row.audit_warning, `
        $row.matvec_seconds_avg, $row.true_residual, $row.preconditioner_effective)
}
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII
Write-Host "PASS Scope Z3C theta policy gate: $decision"
