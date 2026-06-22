param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 300,
    [int]$Threads = 2,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z3a')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z3a_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z3a_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z3a_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z3a_order_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z3A_FMM_ORDER_MATRIX_REPORT.txt'

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
    throw "Refusing to clean Z3A work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z3aDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or $text -eq 'nan') {
        return $Default
    }
    return [double]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Convert-Z3aInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $Default }
    return [int64]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Get-Z3aProp([object]$Row, [string]$Name, [object]$Default = $null) {
    if ($null -eq $Row) { return $Default }
    $prop = $Row.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Write-Z3aArgFile(
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
        '-pbe-tabi-gmres-restart=64',
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

function Read-Z3aResult([string]$RunDir, [int]$ExpectedOrder) {
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
    if (-not ($profile.PSObject.Properties.Name -contains 'fmm_m2l_expansion_order')) {
        throw 'pbe_tabi_profile.csv missing FMM/M2L column: fmm_m2l_expansion_order'
    }
    $fmmOrder = Convert-Z3aInt64 $profile.fmm_m2l_expansion_order
    if ($fmmOrder -ne $ExpectedOrder) {
        throw "Expected FMM/M2L expansion order $ExpectedOrder, got $fmmOrder"
    }
    if ([string]$profile.fmm_m2l_mode -ne 'single-matvec-analytic') {
        throw "Expected FMM/M2L mode single-matvec-analytic, got $($profile.fmm_m2l_mode)"
    }
    return [pscustomobject]@{
        status = if ($energy.status -eq 'ok' -and $profile.status -eq 'ok') { 'ok' } else { "energy=$($energy.status);profile=$($profile.status)" }
        energy_kj_mol = Convert-Z3aDouble $energy.E_rxn_kJmol
        panels = Convert-Z3aInt64 $profile.panels
        iterations = Convert-Z3aInt64 $profile.iterations
        true_residual = Convert-Z3aDouble $profile.true_residual_rel
        taylor_order = Convert-Z3aInt64 (Get-Z3aProp $profile 'taylor_order' '0')
        fmm_m2l_expansion_order = $fmmOrder
        near_integration = [string]$profile.near_integration
        fmm_m2l_mode = [string]$profile.fmm_m2l_mode
        fmm_near_pairs = Convert-Z3aInt64 $profile.fmm_m2l_near_pairs
        fmm_far_pairs = Convert-Z3aInt64 $profile.fmm_m2l_far_pairs
        fmm_finite_outputs = [string]$profile.fmm_m2l_finite_outputs
        tree_vs_direct_matvec_rel = Convert-Z3aDouble $profile.tree_vs_direct_matvec_rel
        matvec_audit_rows = Convert-Z3aInt64 $profile.matvec_audit_rows
        matvec_seconds_avg = Convert-Z3aDouble $profile.matvec_seconds_avg
    }
}

function Invoke-Z3aCase(
    [string]$Name,
    [string]$PqrPath,
    [string]$FixtureClass,
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
    $argPath = Write-Z3aArgFile -RunDir $runDir -PqrPath $PqrPath `
        -TaylorOrder $TaylorOrder -Theta $Theta -SurfacePoints $SurfacePoints `
        -MaxIter $MaxIter -Restart $Restart -SaltM $SaltM `
        -AuditSampleCount $AuditSampleCount
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        throw "Z3A case $Name failed with exit code $($run.exit_code). See $runDir"
    }
    $result = Read-Z3aResult -RunDir $runDir -ExpectedOrder $TaylorOrder
    $relError = [double]::NaN
    if (-not [double]::IsNaN($OracleEnergy) -and [math]::Abs($OracleEnergy) -gt 1.0e-12) {
        $relError = [math]::Abs($result.energy_kj_mol - $OracleEnergy) /
            [math]::Abs($OracleEnergy)
    }
    return [pscustomobject]@{
        condition = $Name
        fixture_class = $FixtureClass
        exit_code = $run.exit_code
        timed_out = $run.timed_out
        status = $result.status
        requested_taylor_order = $TaylorOrder
        taylor_order = $result.taylor_order
        fmm_m2l_expansion_order = $result.fmm_m2l_expansion_order
        theta = $Theta
        near_integration = $result.near_integration
        fmm_m2l_mode = $result.fmm_m2l_mode
        panels = $result.panels
        iterations = $result.iterations
        energy_kj_mol = $result.energy_kj_mol
        oracle_energy_kj_mol = $OracleEnergy
        relative_error = $relError
        true_residual = $result.true_residual
        tree_vs_direct_matvec_rel = $result.tree_vs_direct_matvec_rel
        matvec_audit_rows = $result.matvec_audit_rows
        matvec_seconds_avg = $result.matvec_seconds_avg
        fmm_near_pairs = $result.fmm_near_pairs
        fmm_far_pairs = $result.fmm_far_pairs
        fmm_finite_outputs = $result.fmm_finite_outputs
        wall_seconds = $run.wall_seconds
        peak_working_set_mb = $run.peak_working_set_mb
    }
}

function Write-Z3aSvg([object[]]$Rows, [string]$Path) {
    $points = @($Rows | Where-Object { $_.fixture_class -eq 'analytic' })
    $maxY = [Math]::Max(1.0e-6, ($points | Measure-Object -Property tree_vs_direct_matvec_rel -Maximum).Maximum)
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="760" height="300" viewBox="0 0 760 300">') | Out-Null
    $svg.Add('<rect width="760" height="300" fill="white"/>') | Out-Null
    $svg.Add('<text x="24" y="28" font-family="Arial" font-size="16">Scope Z3A FMM order matrix</text>') | Out-Null
    $x = 60
    foreach ($row in $points) {
        $h = [Math]::Max(2.0, 210.0 * ([double]$row.tree_vs_direct_matvec_rel / $maxY))
        $y = 250.0 - $h
        $fill = if ([int]$row.requested_taylor_order -eq 4) { '#0072b2' } else { '#d55e00' }
        $svg.Add(('<rect x="{0}" y="{1:F1}" width="46" height="{2:F1}" fill="{3}"/>' -f $x, $y, $h, $fill)) | Out-Null
        $svg.Add(('<text x="{0}" y="270" font-family="Arial" font-size="10" transform="rotate(30 {0},270)">p{1} {2}</text>' -f $x, $row.requested_taylor_order, $row.condition)) | Out-Null
        $x += 62
    }
    $svg.Add('<text x="24" y="292" font-family="Arial" font-size="11">Bar height = sampled tree-vs-direct matvec relative RMS</text>') | Out-Null
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$fixtureDir = Join-Path $workRootFull 'fixtures'
New-Item -ItemType Directory -Path $fixtureDir -Force | Out-Null
$bornPqr = Join-Path $fixtureDir 'born_centered.pqr'
$kirkwoodPqr = Join-Path $fixtureDir 'kirkwood_offcenter.pqr'
Write-FagbiBornPqr -Path $bornPqr -Charge 1.0 -RadiusA 3.0
Write-FagbiKirkwoodPqr -Path $kirkwoodPqr -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0

$cases = @(
    [pscustomobject]@{
        name='born'; fixture_class='analytic'; pqr=$bornPqr; max_iter=100
        oracle=Get-FagbiBornEnergyKJ -Charge 1.0 -RadiusA 3.0 -SaltM 0.15
    },
    [pscustomobject]@{
        name='kirkwood'; fixture_class='analytic'; pqr=$kirkwoodPqr; max_iter=140
        oracle=Get-FagbiKirkwoodEnergyKJ -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0 -SaltM 0.15
    }
)

if (-not $Quick) {
    $proteinDir = Join-Path $fixtureDir 'research_medium'
    New-Item -ItemType Directory -Path $proteinDir -Force | Out-Null
    $proteinSource = Resolve-FmmCliFixture -Fixture 'research-medium' -RepoRoot $RepoRoot -RunDir $proteinDir
    if ($null -ne $proteinSource -and (Test-Path -LiteralPath $proteinSource)) {
        $sampled = Join-Path $proteinDir 'research_medium_cap128.pqr'
        [void](Write-FmmCliSampledPqr -InputPath $proteinSource -OutputPath $sampled -Cap 128)
        $cases += [pscustomobject]@{
            name='research_medium_cap128'; fixture_class='protein'; pqr=$sampled
            max_iter=240; oracle=[double]::NaN
        }
    }
}

$thetas = if ($Quick) { @(0.6) } else { @(0.6,0.8) }
$rows = @()
foreach ($case in $cases) {
    foreach ($theta in $thetas) {
        foreach ($order in @(3,4)) {
            $safeTheta = ($theta.ToString('G15', [System.Globalization.CultureInfo]::InvariantCulture) -replace '\.', 'p')
            $name = "{0}_p{1}_theta{2}" -f $case.name, $order, $safeTheta
            $rows += Invoke-Z3aCase -Name $name -PqrPath $case.pqr `
                -FixtureClass $case.fixture_class -OracleEnergy $case.oracle `
                -TaylorOrder $order -Theta $theta -SurfacePoints 8 `
                -MaxIter $case.max_iter -Restart 64 -SaltM 0.15 `
                -AuditSampleCount 128
        }
    }
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
Write-Z3aSvg -Rows $rows -Path $PlotPath

$requiredOk = @($rows | Where-Object {
    $_.status -eq 'ok' -and $_.fmm_m2l_mode -eq 'single-matvec-analytic' -and
    $_.fmm_finite_outputs -eq 'true' -and $_.fmm_m2l_expansion_order -eq $_.requested_taylor_order
})
$decision = 'z3a_runtime_limited'
if ($requiredOk.Count -eq $rows.Count) {
    $order4Wins = 0
    $order3Sufficient = 0
    $thetaDominant = $false
    foreach ($group in $rows | Group-Object {
            ($_.condition -replace '_p[34]_theta[0-9p]+$', '')
        }) {
        $p3 = @($group.Group | Where-Object { [int]$_.requested_taylor_order -eq 3 })[0]
        $p4 = @($group.Group | Where-Object { [int]$_.requested_taylor_order -eq 4 })[0]
        if ($null -ne $p3 -and $null -ne $p4) {
            $m3 = if (Test-FmmCliFinite ([double]$p3.tree_vs_direct_matvec_rel)) {
                [double]$p3.tree_vs_direct_matvec_rel
            } else {
                [double]$p3.relative_error
            }
            $m4 = if (Test-FmmCliFinite ([double]$p4.tree_vs_direct_matvec_rel)) {
                [double]$p4.tree_vs_direct_matvec_rel
            } else {
                [double]$p4.relative_error
            }
            if ($m4 -lt 0.75 * [Math]::Max($m3, 1.0e-16)) { ++$order4Wins }
            if ($m3 -le 1.25 * [Math]::Max($m4, 1.0e-16) -and
                [double]$p3.matvec_seconds_avg -lt 0.9 * [double]$p4.matvec_seconds_avg) {
                ++$order3Sufficient
            }
        }
    }
    if (-not $Quick) {
        foreach ($fixtureGroup in $rows | Group-Object {
                ($_.condition -replace '_p[34]_theta[0-9p]+$', '') -replace '_theta.*$', ''
            }) {
            $theta06 = @($fixtureGroup.Group | Where-Object { [double]$_.theta -eq 0.6 })
            $theta08 = @($fixtureGroup.Group | Where-Object { [double]$_.theta -eq 0.8 })
            if ($theta06.Count -gt 0 -and $theta08.Count -gt 0) {
                $a06 = ($theta06 | Measure-Object -Property tree_vs_direct_matvec_rel -Average).Average
                $a08 = ($theta08 | Measure-Object -Property tree_vs_direct_matvec_rel -Average).Average
                if ([Math]::Abs($a06 - $a08) -gt 0.5 * [Math]::Max([Math]::Min($a06, $a08), 1.0e-16)) {
                    $thetaDominant = $true
                }
            }
        }
    }
    if ($thetaDominant) {
        $decision = 'z3a_theta_dominant'
    } elseif ($order3Sufficient -ge [Math]::Ceiling($rows.Count / 4.0)) {
        $decision = 'z3a_order3_sufficient'
    } elseif ($order4Wins -gt 0) {
        $decision = 'z3a_order4_confirmed'
    } elseif (@($rows | Where-Object { [double]$_.tree_vs_direct_matvec_rel -gt 1.0e-2 }).Count -gt 0) {
        $decision = 'z3a_accuracy_gap'
    } else {
        $decision = 'z3a_order4_confirmed'
    }
}

$decisionRow = [pscustomobject]@{
    decision = $decision
    rows = $rows.Count
    quick = [bool]$Quick
    artifact_matrix = $MatrixCsv
    artifact_plot = $PlotPath
}
$decisionRow | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'scope: Z3A FMM/M2L fixed-order sweep',
    "decision: $decision",
    "rows: $($rows.Count)",
    "quick: $([bool]$Quick)",
    'near_policy: analytic_singular',
    'default_policy_note: grid remains production default; TABI treecode default is unchanged'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z3A Fixed-Order Matrix Report',
    '',
    "Decision: $decision",
    '',
    'Z3A verifies that `-pbe-tabi-taylor-order=3|4` controls the FMM/M2L expansion order.',
    'The full adaptive-by-level policy remains a future scope.',
    '',
    'Artifacts:',
    "- Matrix: $MatrixCsv",
    "- Decision: $DecisionCsv",
    "- Plot: $PlotPath"
)
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII
Write-Host "PASS Scope Z3A FMM order matrix: $decision"
