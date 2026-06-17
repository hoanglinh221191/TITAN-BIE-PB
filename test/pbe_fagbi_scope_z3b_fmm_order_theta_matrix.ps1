param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [switch]$IncludeResearchAlt,
    [int]$CaseTimeoutSeconds = 420,
    [int]$Threads = 2,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z3b')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z3b_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z3b_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z3b_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z3b_order_theta_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z3B_FMM_ORDER_THETA_MATRIX_REPORT.txt'

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
    throw "Refusing to clean Z3B work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z3bDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or
        $text -eq 'nan') {
        return $Default
    }
    return [double]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Convert-Z3bInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or
        $text -eq 'nan') {
        return $Default
    }
    return [int64]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Get-Z3bProp([object]$Row, [string]$Name, [object]$Default = $null) {
    if ($null -eq $Row) { return $Default }
    $prop = $Row.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function New-Z3bNaNRow(
    [string]$Condition,
    [string]$Fixture,
    [string]$FixtureClass,
    [int]$AtomCap,
    [int]$SurfacePoints,
    [int]$TaylorOrder,
    [double]$Theta,
    [double]$OracleEnergy,
    [object]$Run,
    [string]$Status
) {
    return [pscustomobject]@{
        condition = $Condition
        fixture = $Fixture
        fixture_class = $FixtureClass
        atom_cap = $AtomCap
        surface_points = $SurfacePoints
        theta = $Theta
        requested_taylor_order = $TaylorOrder
        fmm_m2l_expansion_order = 0
        status = $Status
        exit_code = if ($null -eq $Run) { -1 } else { $Run.exit_code }
        timed_out = if ($null -eq $Run) { $false } else { $Run.timed_out }
        panels = 0
        iterations = 0
        energy_kj_mol = [double]::NaN
        reference_energy_kj_mol = $OracleEnergy
        reference_kind = if ([double]::IsNaN($OracleEnergy)) { 'sampled_matvec_audit' } else { 'analytic_oracle' }
        relative_error = [double]::NaN
        energy_delta_kj_mol = [double]::NaN
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
        fmm_near_pairs = 0
        fmm_far_pairs = 0
        fmm_finite_outputs = ''
    }
}

function Write-Z3bArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [int]$TaylorOrder,
    [double]$Theta,
    [int]$SurfacePoints,
    [int]$MaxIter,
    [int]$Restart,
    [double]$SaltM,
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
        "-pbe-salt=$(& $format $SaltM)",
        '-pbe-temp=298.15',
        '-pbe-timing=detail',
        "-nt=$Threads",
        '-pbe-tabi-sum=treecode',
        '-pbe-tabi-formulation=jcp2013-treecode',
        '-pbe-tabi-fmm-m2l=single-matvec-analytic',
        "-pbe-tabi-surface-points=$SurfacePoints",
        '-pbe-tabi-tree-leaf-size=64',
        "-pbe-tabi-opening-theta=$(& $format $Theta)",
        "-pbe-tabi-taylor-order=$TaylorOrder",
        '-pbe-tabi-near-integration=analytic_singular',
        '-pbe-tabi-quad-buffer=off',
        '-pbe-tabi-precond=fmm-nearfield-ilu0',
        '-pbe-tabi-precond-ilu-neighbor-cap=64',
        "-pbe-tabi-gmres-restart=$Restart",
        '-pbe-tabi-true-residual-probe=restart',
        '-pbe-tabi-true-residual-probe-scope=sampled',
        '-pbe-tabi-true-residual-sample-count=256'
    )
    if ($AuditSampleCount -gt 0) {
        $argLines += '-pbe-tabi-matvec-audit=sampled'
        $argLines += "-pbe-tabi-matvec-audit-sample-count=$AuditSampleCount"
    } else {
        $argLines += '-pbe-tabi-matvec-audit=off'
    }
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

function Resolve-Z3bFixture(
    [string]$Fixture,
    [string]$RunDir
) {
    $resolved = Resolve-FmmCliFixture -Fixture $Fixture -RepoRoot $RepoRoot `
        -RunDir $RunDir
    if ($null -ne $resolved -and (Test-Path -LiteralPath $resolved)) {
        return $resolved
    }
    if ($Fixture -ne 'barnase') { return $resolved }

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

function Read-Z3bResult([string]$RunDir, [int]$ExpectedOrder) {
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
    foreach ($column in @('fmm_m2l_expansion_order','fmm_m2l_mode',
            'preconditioner_effective','near_integration','quad_buffer')) {
        if (-not ($profile.PSObject.Properties.Name -contains $column)) {
            throw "pbe_tabi_profile.csv missing required Z3B column: $column"
        }
    }
    $fmmOrder = Convert-Z3bInt64 $profile.fmm_m2l_expansion_order
    if ($fmmOrder -ne $ExpectedOrder) {
        throw "Expected FMM/M2L expansion order $ExpectedOrder, got $fmmOrder"
    }
    if ([string]$profile.fmm_m2l_mode -ne 'single-matvec-analytic') {
        throw "Expected FMM/M2L mode single-matvec-analytic, got $($profile.fmm_m2l_mode)"
    }
    return [pscustomobject]@{
        status = if ($energy.status -eq 'ok' -and $profile.status -eq 'ok') { 'ok' } else { "energy=$($energy.status);profile=$($profile.status)" }
        energy_kj_mol = Convert-Z3bDouble $energy.E_rxn_kJmol
        panels = Convert-Z3bInt64 $profile.panels
        iterations = Convert-Z3bInt64 $profile.iterations
        true_residual = Convert-Z3bDouble $profile.true_residual_rel
        fmm_m2l_expansion_order = $fmmOrder
        near_integration = [string]$profile.near_integration
        quad_buffer = [string]$profile.quad_buffer
        fmm_m2l_mode = [string]$profile.fmm_m2l_mode
        fmm_near_pairs = Convert-Z3bInt64 $profile.fmm_m2l_near_pairs
        fmm_far_pairs = Convert-Z3bInt64 $profile.fmm_m2l_far_pairs
        fmm_finite_outputs = [string]$profile.fmm_m2l_finite_outputs
        tree_vs_direct_matvec_rel = Convert-Z3bDouble $profile.tree_vs_direct_matvec_rel
        matvec_audit_rows = Convert-Z3bInt64 $profile.matvec_audit_rows
        matvec_seconds_avg = Convert-Z3bDouble $profile.matvec_seconds_avg
        preconditioner_effective = [string]$profile.preconditioner_effective
    }
}

function Invoke-Z3bCase(
    [string]$Name,
    [string]$Fixture,
    [string]$PqrPath,
    [string]$FixtureClass,
    [int]$AtomCap,
    [double]$OracleEnergy,
    [int]$TaylorOrder,
    [double]$Theta,
    [int]$SurfacePoints,
    [int]$MaxIter,
    [int]$Restart,
    [double]$SaltM,
    [int]$AuditSampleCount
) {
    $runDir = Join-Path $workRootFull $Name
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $argPath = Write-Z3bArgFile -RunDir $runDir -PqrPath $PqrPath `
        -TaylorOrder $TaylorOrder -Theta $Theta -SurfacePoints $SurfacePoints `
        -MaxIter $MaxIter -Restart $Restart -SaltM $SaltM `
        -AuditSampleCount $AuditSampleCount
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        return New-Z3bNaNRow -Condition $Name -Fixture $Fixture `
            -FixtureClass $FixtureClass -AtomCap $AtomCap `
            -SurfacePoints $SurfacePoints -TaylorOrder $TaylorOrder `
            -Theta $Theta -OracleEnergy $OracleEnergy -Run $run `
            -Status $(if ($run.timed_out) { 'timeout_not_verified' } else { "exit_$($run.exit_code)" })
    }
    try {
        $result = Read-Z3bResult -RunDir $runDir -ExpectedOrder $TaylorOrder
    } catch {
        return New-Z3bNaNRow -Condition $Name -Fixture $Fixture `
            -FixtureClass $FixtureClass -AtomCap $AtomCap `
            -SurfacePoints $SurfacePoints -TaylorOrder $TaylorOrder `
            -Theta $Theta -OracleEnergy $OracleEnergy -Run $run `
            -Status "profile_error: $($_.Exception.Message)"
    }
    $relError = [double]::NaN
    $energyDelta = [double]::NaN
    if (-not [double]::IsNaN($OracleEnergy) -and [math]::Abs($OracleEnergy) -gt 1.0e-12) {
        $energyDelta = $result.energy_kj_mol - $OracleEnergy
        $relError = [math]::Abs($energyDelta) / [math]::Abs($OracleEnergy)
    }
    return [pscustomobject]@{
        condition = $Name
        fixture = $Fixture
        fixture_class = $FixtureClass
        atom_cap = $AtomCap
        surface_points = $SurfacePoints
        theta = $Theta
        requested_taylor_order = $TaylorOrder
        fmm_m2l_expansion_order = $result.fmm_m2l_expansion_order
        status = $result.status
        exit_code = $run.exit_code
        timed_out = $run.timed_out
        panels = $result.panels
        iterations = $result.iterations
        energy_kj_mol = $result.energy_kj_mol
        reference_energy_kj_mol = $OracleEnergy
        reference_kind = if ([double]::IsNaN($OracleEnergy)) { 'sampled_matvec_audit' } else { 'analytic_oracle' }
        relative_error = $relError
        energy_delta_kj_mol = $energyDelta
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
        fmm_near_pairs = $result.fmm_near_pairs
        fmm_far_pairs = $result.fmm_far_pairs
        fmm_finite_outputs = $result.fmm_finite_outputs
    }
}

function Get-Z3bMetric([object]$Row) {
    if ($null -eq $Row) { return [double]::NaN }
    $audit = Convert-Z3bDouble $Row.tree_vs_direct_matvec_rel
    if (Test-FmmCliFinite $audit -and $audit -gt 0.0) { return $audit }
    $err = Convert-Z3bDouble $Row.relative_error
    if (Test-FmmCliFinite $err -and $err -gt 0.0) { return $err }
    return [double]::NaN
}

function Write-Z3bSvg([object[]]$Rows, [string]$Path) {
    $points = @($Rows | Where-Object { $_.status -eq 'ok' })
    $metrics = @($points | ForEach-Object { Get-Z3bMetric $_ } |
        Where-Object { Test-FmmCliFinite $_ })
    $maxY = if ($metrics.Count -gt 0) {
        [Math]::Max(1.0e-6, ($metrics | Measure-Object -Maximum).Maximum)
    } else {
        1.0
    }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="980" height="340" viewBox="0 0 980 340">') | Out-Null
    $svg.Add('<rect width="980" height="340" fill="white"/>') | Out-Null
    $svg.Add('<text x="24" y="28" font-family="Arial" font-size="16">Scope Z3B FMM order/theta matrix</text>') | Out-Null
    $svg.Add('<text x="24" y="48" font-family="Arial" font-size="11">Bar height = sampled audit RMS when available, otherwise analytic oracle relative error</text>') | Out-Null
    $x = 36
    foreach ($row in $points) {
        $metric = Get-Z3bMetric $row
        if (-not (Test-FmmCliFinite $metric)) { $metric = 0.0 }
        $h = [Math]::Max(2.0, 230.0 * ($metric / $maxY))
        $y = 280.0 - $h
        $fill = if ([int]$row.requested_taylor_order -eq 4) { '#0072b2' } else { '#d55e00' }
        $svg.Add(('<rect x="{0}" y="{1:F1}" width="22" height="{2:F1}" fill="{3}"/>' -f $x, $y, $h, $fill)) | Out-Null
        $label = ('{0} p{1} t{2}' -f $row.fixture, $row.requested_taylor_order,
            ([double]$row.theta).ToString('G3', [System.Globalization.CultureInfo]::InvariantCulture))
        $svg.Add(('<text x="{0}" y="302" font-family="Arial" font-size="8" transform="rotate(45 {0},302)">{1}</text>' -f $x, $label)) | Out-Null
        $x += 32
        if ($x -gt 930) { break }
    }
    $svg.Add('<rect x="800" y="20" width="14" height="14" fill="#d55e00"/><text x="820" y="32" font-family="Arial" font-size="11">order 3</text>') | Out-Null
    $svg.Add('<rect x="800" y="40" width="14" height="14" fill="#0072b2"/><text x="820" y="52" font-family="Arial" font-size="11">order 4</text>') | Out-Null
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

function Get-Z3bDecision([object[]]$Rows) {
    $okRows = @($Rows | Where-Object {
        $_.status -eq 'ok' -and
        $_.fmm_m2l_mode -eq 'single-matvec-analytic' -and
        $_.fmm_finite_outputs -eq 'true' -and
        $_.fmm_m2l_expansion_order -eq $_.requested_taylor_order -and
        $_.near_integration -eq 'analytic_singular' -and
        $_.quad_buffer -eq 'off'
    })
    if ($okRows.Count -ne $Rows.Count) { return 'z3b_runtime_limited' }

    $requiredAuditFailures = @($okRows | Where-Object {
        $_.matvec_audit_rows -gt 0 -and
        (Test-FmmCliFinite ([double]$_.tree_vs_direct_matvec_rel)) -and
        [double]$_.tree_vs_direct_matvec_rel -gt 1.0e-2
    })
    # Z3B isolates FMM order/theta. Low-panel analytic oracle error is already
    # covered by Z1C/Z1D surface convergence, so use a loose diagnostic gate.
    $analyticFailures = @($okRows | Where-Object {
        $_.fixture_class -eq 'analytic' -and
        (Test-FmmCliFinite ([double]$_.relative_error)) -and
        [double]$_.relative_error -gt 1.0e-1
    })
    if ($requiredAuditFailures.Count -gt 0 -or $analyticFailures.Count -gt 0) {
        return 'z3b_accuracy_gap'
    }

    $order4Wins = 0
    $order3Sufficient = 0
    $orderPairs = 0
    foreach ($group in $okRows | Group-Object {
            '{0}|{1}|{2}' -f $_.fixture, $_.fixture_class, $_.theta
        }) {
        $p3 = @($group.Group | Where-Object { [int]$_.requested_taylor_order -eq 3 })[0]
        $p4 = @($group.Group | Where-Object { [int]$_.requested_taylor_order -eq 4 })[0]
        if ($null -eq $p3 -or $null -eq $p4) { continue }
        $m3 = Get-Z3bMetric $p3
        $m4 = Get-Z3bMetric $p4
        if (-not (Test-FmmCliFinite $m3) -or -not (Test-FmmCliFinite $m4)) { continue }
        ++$orderPairs
        if ($m4 -le 0.75 * [Math]::Max($m3, 1.0e-16) -and
            [double]$p4.matvec_seconds_avg -le 1.35 * [double]$p3.matvec_seconds_avg) {
            ++$order4Wins
        }
        if ($m3 -le 1.25 * [Math]::Max($m4, 1.0e-16) -and
            [double]$p3.matvec_seconds_avg -le 0.90 * [double]$p4.matvec_seconds_avg) {
            ++$order3Sufficient
        }
    }

    $thetaDominant = 0
    $thetaGroups = 0
    foreach ($group in $okRows | Group-Object {
            '{0}|{1}|{2}' -f $_.fixture, $_.fixture_class, $_.requested_taylor_order
        }) {
        $metrics = @($group.Group | ForEach-Object { Get-Z3bMetric $_ } |
            Where-Object { Test-FmmCliFinite $_ })
        if ($metrics.Count -lt 2) { continue }
        ++$thetaGroups
        $minMetric = [Math]::Max(1.0e-16, ($metrics | Measure-Object -Minimum).Minimum)
        $maxMetric = ($metrics | Measure-Object -Maximum).Maximum
        if (($maxMetric / $minMetric) -ge 1.5) { ++$thetaDominant }
    }

    if ($thetaGroups -gt 0 -and $thetaDominant -ge [Math]::Ceiling($thetaGroups / 2.0)) {
        if ($order4Wins -gt 0 -and $order3Sufficient -gt 0) {
            return 'z3b_adaptive_order_justified'
        }
        return 'z3b_theta_dominant'
    }
    if ($orderPairs -gt 0 -and $order3Sufficient -ge [Math]::Ceiling($orderPairs / 2.0)) {
        return 'z3b_order3_sufficient'
    }
    if ($order4Wins -gt 0) { return 'z3b_order4_confirmed' }
    return 'z3b_adaptive_order_justified'
}

$fixtureDir = Join-Path $workRootFull 'fixtures'
New-Item -ItemType Directory -Path $fixtureDir -Force | Out-Null
$bornPqr = Join-Path $fixtureDir 'born_centered.pqr'
$kirkwoodPqr = Join-Path $fixtureDir 'kirkwood_offcenter.pqr'
Write-FagbiBornPqr -Path $bornPqr -Charge 1.0 -RadiusA 3.0
Write-FagbiKirkwoodPqr -Path $kirkwoodPqr -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0

$cases = New-Object System.Collections.Generic.List[object]
$cases.Add([pscustomobject]@{
    fixture='born'; fixture_class='analytic'; atom_cap=0; pqr=$bornPqr
    max_iter=100; surface_points=8; audit=128
    oracle=Get-FagbiBornEnergyKJ -Charge 1.0 -RadiusA 3.0 -SaltM 0.15
}) | Out-Null
$cases.Add([pscustomobject]@{
    fixture='kirkwood'; fixture_class='analytic'; atom_cap=0; pqr=$kirkwoodPqr
    max_iter=140; surface_points=8; audit=128
    oracle=Get-FagbiKirkwoodEnergyKJ -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0 -SaltM 0.15
}) | Out-Null

if (-not $Quick) {
    foreach ($fixtureName in @('barnase','research-medium')) {
        $sampleDir = Join-Path $fixtureDir $fixtureName
        New-Item -ItemType Directory -Path $sampleDir -Force | Out-Null
        $source = Resolve-Z3bFixture -Fixture $fixtureName -RunDir $sampleDir
        if ($null -eq $source -or -not (Test-Path -LiteralPath $source)) {
            $cases.Add([pscustomobject]@{
                fixture=$fixtureName; fixture_class='protein'; atom_cap=128
                pqr=''; max_iter=0; surface_points=8; audit=512
                oracle=[double]::NaN; missing=$true
            }) | Out-Null
            continue
        }
        $sampled = Join-Path $sampleDir ("{0}_cap128.pqr" -f (Get-FmmCliSafeName $fixtureName))
        [void](Write-FmmCliSampledPqr -InputPath $source -OutputPath $sampled -Cap 128)
        $cases.Add([pscustomobject]@{
            fixture=$fixtureName; fixture_class='protein'; atom_cap=128
            pqr=$sampled; max_iter=240; surface_points=8; audit=512
            oracle=[double]::NaN; missing=$false
        }) | Out-Null
    }
    if ($IncludeResearchAlt) {
        $sampleDir = Join-Path $fixtureDir 'research-alt'
        New-Item -ItemType Directory -Path $sampleDir -Force | Out-Null
        $source = Resolve-FmmCliFixture -Fixture 'research-alt' -RepoRoot $RepoRoot `
            -RunDir $sampleDir
        $sampled = Join-Path $sampleDir 'research_alt_cap128.pqr'
        [void](Write-FmmCliSampledPqr -InputPath $source -OutputPath $sampled -Cap 128)
        $cases.Add([pscustomobject]@{
            fixture='research-alt'; fixture_class='protein'; atom_cap=128
            pqr=$sampled; max_iter=240; surface_points=8; audit=512
            oracle=[double]::NaN; missing=$false
        }) | Out-Null
    }
}

$thetas = if ($Quick) { @(0.6) } else { @(0.45,0.6,0.8) }
$rows = New-Object System.Collections.Generic.List[object]
foreach ($case in $cases) {
    foreach ($theta in $thetas) {
        foreach ($order in @(3,4)) {
            $safeTheta = ($theta.ToString('G15', [System.Globalization.CultureInfo]::InvariantCulture) -replace '\.', 'p')
            $name = "{0}_p{1}_theta{2}" -f (Get-FmmCliSafeName $case.fixture), $order, $safeTheta
            if ((Get-Z3bProp $case 'missing' $false)) {
                $rows.Add((New-Z3bNaNRow -Condition $name -Fixture $case.fixture `
                    -FixtureClass $case.fixture_class -AtomCap $case.atom_cap `
                    -SurfacePoints $case.surface_points -TaylorOrder $order `
                    -Theta $theta -OracleEnergy $case.oracle -Run $null `
                    -Status 'skipped_missing_fixture')) | Out-Null
                continue
            }
            $rows.Add((Invoke-Z3bCase -Name $name -Fixture $case.fixture `
                -PqrPath $case.pqr -FixtureClass $case.fixture_class `
                -AtomCap $case.atom_cap -OracleEnergy $case.oracle `
                -TaylorOrder $order -Theta $theta `
                -SurfacePoints $case.surface_points -MaxIter $case.max_iter `
                -Restart 64 -SaltM 0.15 -AuditSampleCount $case.audit)) | Out-Null
        }
    }
}

$rowArray = $rows.ToArray()
$rowArray | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
Write-Z3bSvg -Rows $rowArray -Path $PlotPath

$decision = Get-Z3bDecision -Rows $rowArray
[pscustomobject]@{
    decision = $decision
    rows = $rowArray.Count
    quick = [bool]$Quick
    include_research_alt = [bool]$IncludeResearchAlt
    artifact_matrix = $MatrixCsv
    artifact_plot = $PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'scope: Z3B FMM/M2L fixed-order/theta full matrix',
    "decision: $decision",
    "rows: $($rowArray.Count)",
    "quick: $([bool]$Quick)",
    "include_research_alt: $([bool]$IncludeResearchAlt)",
    'near_policy: analytic_singular',
    'quad_buffer: off',
    'preconditioner: fmm-nearfield-ilu0',
    'analytic_oracle_gate: 0.10',
    'default_policy_note: grid remains production default; TABI treecode default is unchanged'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z3B Order/Theta Matrix Report',
    '',
    "Decision: $decision",
    '',
    'Z3B measures fixed FMM/M2L expansion order and opening-theta behavior before any adaptive-order implementation.',
    '',
    'Artifacts:',
    "- Matrix: $MatrixCsv",
    "- Decision: $DecisionCsv",
    "- Summary: $SummaryPath",
    "- Plot: $PlotPath",
    '',
    'Rows:'
)
foreach ($row in $rowArray) {
    $report += ('- {0}: fixture={1}, p={2}, theta={3}, status={4}, audit={5}, rel_error={6}, matvec_s={7}, precond={8}' -f `
        $row.condition, $row.fixture, $row.requested_taylor_order, $row.theta, `
        $row.status, $row.tree_vs_direct_matvec_rel, $row.relative_error, `
        $row.matvec_seconds_avg, $row.preconditioner_effective)
}
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "PASS Scope Z3B FMM order/theta matrix: $decision"
