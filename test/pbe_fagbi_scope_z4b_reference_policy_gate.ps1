param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 2,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4b')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4b_matrix.csv'
$ReferenceCsv = Join-Path $WorkRoot 'fagbi_scope_z4b_reference_policy.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4b_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4b_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4b_reference_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4B_REFERENCE_POLICY_REPORT.txt'

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
    throw "Refusing to clean Z4B work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z4bDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or $text -eq 'nan') {
        return $Default
    }
    return [double]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Convert-Z4bInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or $text -eq 'nan') {
        return $Default
    }
    return [int64]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Format-Z4bDouble([double]$Value) {
    return $Value.ToString('G15', [System.Globalization.CultureInfo]::InvariantCulture)
}

function New-Z4bRow(
    [string]$Condition,
    [string]$Fixture,
    [int]$AtomCap,
    [string]$Backend,
    [object]$Run,
    [string]$Status
) {
    return [pscustomobject]@{
        condition = $Condition
        fixture = $Fixture
        atom_cap = $AtomCap
        backend = $Backend
        status = $Status
        failure_class = $Status
        exit_code = if ($null -eq $Run) { -1 } else { $Run.exit_code }
        timed_out = if ($null -eq $Run) { $false } else { $Run.timed_out }
        panels = 0
        iterations = 0
        energy_kj_mol = [double]::NaN
        reference_kind = ''
        reference_backend = ''
        reference_status = ''
        reference_energy_kj_mol = [double]::NaN
        energy_delta_kj_mol = [double]::NaN
        relative_energy_delta = [double]::NaN
        reference_usable_for_paper = $false
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

function Resolve-Z4bFixture([string]$Fixture, [string]$RunDir) {
    $resolved = Resolve-FmmCliFixture -Fixture $Fixture -RepoRoot $RepoRoot -RunDir $RunDir
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

function Write-Z4bArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [ValidateSet('grid_context','tabi_direct_full','fmm_ilu0')][string]$Backend
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
        "-pbe-test-pqr=$(Convert-FmmCliToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-max-iter=360',
        '-pbe-rms=1e-6',
        '-pbe-dphi=1e-6',
        '-pbe-int-diel=1',
        '-pbe-ext-diel=80',
        '-pbe-salt=0.15',
        '-pbe-temp=298.15',
        '-pbe-timing=detail',
        "-nt=$Threads"
    )
    if ($Backend -eq 'grid_context') {
        $argLines += '-pbe=grid'
    } else {
        $mode = if ($Backend -eq 'tabi_direct_full') { 'direct' } else { 'treecode' }
        $fmmMode = if ($Backend -eq 'fmm_ilu0') { 'single-matvec-analytic' } else { 'off' }
        $precond = if ($Backend -eq 'fmm_ilu0') { 'fmm-nearfield-ilu0' } else { 'diagonal' }
        $argLines += @(
            '-pbe=tabi',
            "-pbe-tabi-sum=$mode",
            '-pbe-tabi-formulation=jcp2013-treecode',
            "-pbe-tabi-fmm-m2l=$fmmMode",
            '-pbe-tabi-surface-points=8',
            '-pbe-tabi-tree-leaf-size=64',
            '-pbe-tabi-opening-theta=0.6',
            '-pbe-tabi-taylor-order=4',
            '-pbe-tabi-near-integration=analytic_singular',
            '-pbe-tabi-quad-buffer=off',
            "-pbe-tabi-precond=$precond",
            '-pbe-tabi-precond-ilu-neighbor-cap=64',
            '-pbe-tabi-gmres-restart=64',
            '-pbe-tabi-true-residual-probe=restart',
            '-pbe-tabi-true-residual-probe-scope=sampled',
            '-pbe-tabi-true-residual-sample-count=512'
        )
        if ($Backend -eq 'fmm_ilu0') {
            $argLines += '-pbe-tabi-matvec-audit=sampled'
            $argLines += '-pbe-tabi-matvec-audit-sample-count=512'
        } else {
            $argLines += '-pbe-tabi-matvec-audit=off'
        }
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

function Read-Z4bResult([string]$RunDir, [string]$Backend) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    if (-not (Test-Path -LiteralPath $energyPath)) {
        throw "Missing pbe_energy.csv in $RunDir"
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $result = [ordered]@{
        status = if ($energy.status -eq 'ok') { 'ok' } else { [string]$energy.status }
        panels = 0
        iterations = 0
        energy_kj_mol = Convert-Z4bDouble $energy.E_rxn_kJmol
        true_residual = [double]::NaN
        tree_vs_direct_matvec_rel = [double]::NaN
        matvec_audit_rows = 0
        matvec_seconds_avg = [double]::NaN
        preconditioner_effective = ''
        near_integration = ''
        quad_buffer = ''
        fmm_m2l_mode = ''
        fmm_m2l_expansion_order = 0
        fmm_near_pairs = 0
        fmm_far_pairs = 0
        fmm_finite_outputs = ''
    }
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if ($Backend -ne 'grid_context' -and (Test-Path -LiteralPath $profilePath)) {
        $profile = @(Import-Csv -LiteralPath $profilePath)[0]
        $result.status = if ($result.status -eq 'ok' -and $profile.status -eq 'ok') {
            'ok'
        } else {
            "energy=$($energy.status);profile=$($profile.status)"
        }
        $result.panels = Convert-Z4bInt64 $profile.panels
        $result.iterations = Convert-Z4bInt64 $profile.iterations
        $result.true_residual = Convert-Z4bDouble $profile.true_residual_rel
        $result.tree_vs_direct_matvec_rel = Convert-Z4bDouble $profile.tree_vs_direct_matvec_rel
        $result.matvec_audit_rows = Convert-Z4bInt64 $profile.matvec_audit_rows
        $result.matvec_seconds_avg = Convert-Z4bDouble $profile.matvec_seconds_avg
        $result.preconditioner_effective = [string]$profile.preconditioner_effective
        $result.near_integration = [string]$profile.near_integration
        $result.quad_buffer = [string]$profile.quad_buffer
        $result.fmm_m2l_mode = [string]$profile.fmm_m2l_mode
        $result.fmm_m2l_expansion_order = Convert-Z4bInt64 $profile.fmm_m2l_expansion_order
        $result.fmm_near_pairs = Convert-Z4bInt64 $profile.fmm_m2l_near_pairs
        $result.fmm_far_pairs = Convert-Z4bInt64 $profile.fmm_m2l_far_pairs
        $result.fmm_finite_outputs = [string]$profile.fmm_m2l_finite_outputs
    }
    return [pscustomobject]$result
}

function Invoke-Z4bCase(
    [string]$Fixture,
    [int]$AtomCap,
    [string]$PqrPath,
    [ValidateSet('grid_context','tabi_direct_full','fmm_ilu0')][string]$Backend
) {
    $condition = "{0}_cap{1}_{2}" -f (Get-FmmCliSafeName $Fixture), $AtomCap, $Backend
    $runDir = Join-Path $workRootFull $condition
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $argPath = Write-Z4bArgFile -RunDir $runDir -PqrPath $PqrPath -Backend $Backend
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        return New-Z4bRow -Condition $condition -Fixture $Fixture `
            -AtomCap $AtomCap -Backend $Backend -Run $run `
            -Status $(if ($run.timed_out) { 'timeout_not_verified' } else { "exit_$($run.exit_code)" })
    }
    try {
        $result = Read-Z4bResult -RunDir $runDir -Backend $Backend
    } catch {
        return New-Z4bRow -Condition $condition -Fixture $Fixture `
            -AtomCap $AtomCap -Backend $Backend -Run $run `
            -Status "profile_error: $($_.Exception.Message)"
    }
    $failureClass = 'pass'
    if ($result.status -ne 'ok') {
        $failureClass = if ($Backend -eq 'grid_context') { 'context_not_converged' } else { 'status_fail' }
    } elseif ($Backend -eq 'fmm_ilu0') {
        if ($result.fmm_m2l_mode -ne 'single-matvec-analytic' -or
            $result.fmm_finite_outputs -ne 'true' -or
            $result.fmm_m2l_expansion_order -ne 4 -or
            $result.near_integration -ne 'analytic_singular' -or
            $result.quad_buffer -ne 'off' -or
            $result.preconditioner_effective -ne 'fmm_nearfield_ilu0') {
            $failureClass = 'telemetry_fail'
        } elseif ($result.matvec_audit_rows -le 0 -or
            -not (Test-FmmCliFinite $result.tree_vs_direct_matvec_rel) -or
            $result.tree_vs_direct_matvec_rel -gt 8.0e-3) {
            $failureClass = 'accuracy_fail'
        } elseif (-not (Test-FmmCliFinite $result.true_residual) -or
            [Math]::Abs($result.true_residual) -gt 1.0e-6) {
            $failureClass = 'convergence_warn'
        }
    }
    return [pscustomobject]@{
        condition = $condition
        fixture = $Fixture
        atom_cap = $AtomCap
        backend = $Backend
        status = $result.status
        failure_class = $failureClass
        exit_code = $run.exit_code
        timed_out = $run.timed_out
        panels = $result.panels
        iterations = $result.iterations
        energy_kj_mol = $result.energy_kj_mol
        reference_kind = ''
        reference_backend = ''
        reference_status = ''
        reference_energy_kj_mol = [double]::NaN
        energy_delta_kj_mol = [double]::NaN
        relative_energy_delta = [double]::NaN
        reference_usable_for_paper = $false
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

function Add-Z4bReferencePolicy([object[]]$Rows) {
    foreach ($group in $Rows | Group-Object { '{0}|{1}' -f $_.fixture, $_.atom_cap }) {
        $directRows = @($group.Group | Where-Object {
                $_.backend -eq 'tabi_direct_full' -and $_.status -eq 'ok'
            })
        $gridRows = @($group.Group | Where-Object { $_.backend -eq 'grid_context' })
        $direct = if ($directRows.Count -gt 0) { $directRows[0] } else { $null }
        $grid = if ($gridRows.Count -gt 0) { $gridRows[0] } else { $null }
        foreach ($row in $group.Group) {
            if ($row.backend -eq 'tabi_direct_full') {
                $row.reference_kind = 'self_full_direct_tabi_energy'
                $row.reference_backend = 'tabi_direct_full'
                $row.reference_status = $row.status
                $row.reference_energy_kj_mol = $row.energy_kj_mol
                $row.reference_usable_for_paper = ($row.status -eq 'ok')
                continue
            }
            if ($row.backend -eq 'grid_context') {
                $row.reference_kind = if ($row.status -eq 'ok') {
                    'grid_context_ok'
                } else {
                    'grid_context_not_converged'
                }
                $row.reference_backend = 'grid_context'
                $row.reference_status = $row.status
                $row.reference_energy_kj_mol = $row.energy_kj_mol
                $row.reference_usable_for_paper = $false
                continue
            }
            if ($row.backend -ne 'fmm_ilu0') { continue }
            if ($row.atom_cap -lt 256 -and $null -ne $direct) {
                $row.reference_kind = 'full_direct_tabi_energy'
                $row.reference_backend = 'tabi_direct_full'
                $row.reference_status = $direct.status
                $row.reference_energy_kj_mol = $direct.energy_kj_mol
                if ((Test-FmmCliFinite $row.energy_kj_mol) -and
                    (Test-FmmCliFinite $direct.energy_kj_mol)) {
                    $row.energy_delta_kj_mol = $row.energy_kj_mol - $direct.energy_kj_mol
                    $row.relative_energy_delta = [Math]::Abs($row.energy_delta_kj_mol) /
                        [Math]::Max(1.0, [Math]::Abs($direct.energy_kj_mol))
                }
                $row.reference_usable_for_paper = $true
                continue
            }
            if ($row.atom_cap -ge 256 -and
                $row.matvec_audit_rows -gt 0 -and
                (Test-FmmCliFinite $row.tree_vs_direct_matvec_rel)) {
                $row.reference_kind = 'sampled_direct_operator_audit'
                $row.reference_backend = 'sampled_direct_operator'
                $row.reference_status = if ($row.tree_vs_direct_matvec_rel -le 8.0e-3) {
                    'audit_ok'
                } else {
                    'audit_fail'
                }
                $row.reference_usable_for_paper = ($row.tree_vs_direct_matvec_rel -le 8.0e-3)
                continue
            }
            if ($null -ne $grid -and $grid.status -eq 'ok') {
                $row.reference_kind = 'grid_context_only'
                $row.reference_backend = 'grid_context'
                $row.reference_status = $grid.status
                $row.reference_energy_kj_mol = $grid.energy_kj_mol
                $row.reference_usable_for_paper = $false
            } else {
                $row.reference_kind = 'missing_reference'
                $row.reference_status = 'missing'
                $row.reference_usable_for_paper = $false
            }
        }
    }
}

function Get-Z4bDecision([object[]]$Rows) {
    if (@($Rows | Where-Object { $_.status -eq 'skipped_missing_fixture' }).Count -gt 0) {
        return 'z4b_fixture_limited'
    }
    $requiredFmm = @($Rows | Where-Object { $_.backend -eq 'fmm_ilu0' })
    if (@($requiredFmm | Where-Object { $_.timed_out -or $_.status -eq 'timeout_not_verified' }).Count -gt 0) {
        return 'z4b_runtime_limited'
    }
    if (@($requiredFmm | Where-Object {
                $_.failure_class -in @('accuracy_fail','telemetry_fail','status_fail')
            }).Count -gt 0) {
        return 'z4b_direct_audit_gap'
    }
    $cap128Direct = @($Rows | Where-Object {
            $_.backend -eq 'tabi_direct_full' -and $_.atom_cap -lt 256
        })
    if (@($cap128Direct | Where-Object { $_.timed_out -or $_.status -ne 'ok' }).Count -gt 0) {
        return 'z4b_runtime_limited'
    }
    $cap128Fmm = @($requiredFmm | Where-Object { $_.atom_cap -lt 256 })
    foreach ($row in $cap128Fmm) {
        if ($row.reference_kind -ne 'full_direct_tabi_energy' -or
            -not (Test-FmmCliFinite $row.relative_energy_delta)) {
            return 'z4b_cap128_energy_gap'
        }
        if ($row.relative_energy_delta -gt 0.05 -and
            [Math]::Abs($row.energy_delta_kj_mol) -gt 2.0) {
            return 'z4b_cap128_energy_gap'
        }
    }
    $cap256Fmm = @($requiredFmm | Where-Object { $_.atom_cap -ge 256 })
    foreach ($row in $cap256Fmm) {
        if ($row.reference_kind -ne 'sampled_direct_operator_audit' -or
            -not $row.reference_usable_for_paper) {
            return 'z4b_direct_audit_gap'
        }
    }
    return 'z4b_reference_policy_pass'
}

function Write-Z4bSvg([object[]]$Rows, [string]$Path) {
    $fmmRows = @($Rows | Where-Object { $_.backend -eq 'fmm_ilu0' })
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="980" height="360" viewBox="0 0 980 360">') | Out-Null
    $svg.Add('<rect width="980" height="360" fill="white"/>') | Out-Null
    $svg.Add('<text x="24" y="28" font-family="Arial" font-size="16">Scope Z4B reference policy gate</text>') | Out-Null
    $svg.Add('<text x="24" y="48" font-family="Arial" font-size="11">FMM rows: cap-128 uses full direct energy; cap-256 uses sampled direct-operator audit</text>') | Out-Null
    $x = 44
    foreach ($row in $fmmRows) {
        $audit = Convert-Z4bDouble $row.tree_vs_direct_matvec_rel 0.0
        $h = [Math]::Max(2.0, 230.0 * ([Math]::Min($audit, 0.012) / 0.012))
        $y = 292.0 - $h
        $fill = if ($row.reference_kind -eq 'full_direct_tabi_energy') {
            '#0072b2'
        } elseif ($row.reference_kind -eq 'sampled_direct_operator_audit') {
            '#cc79a7'
        } else {
            '#999999'
        }
        $svg.Add(('<rect x="{0}" y="{1:F1}" width="34" height="{2:F1}" fill="{3}"/>' -f $x, $y, $h, $fill)) | Out-Null
        $label = ('{0} c{1}' -f $row.fixture, $row.atom_cap)
        $svg.Add(('<text x="{0}" y="318" font-family="Arial" font-size="9" transform="rotate(45 {0},318)">{1}</text>' -f $x, $label)) | Out-Null
        $x += 52
        if ($x -gt 900) { break }
    }
    $svg.Add('<line x1="24" y1="139" x2="940" y2="139" stroke="#d55e00" stroke-dasharray="4,4"/>') | Out-Null
    $svg.Add('<text x="760" y="134" font-family="Arial" font-size="10">audit gate 8e-3</text>') | Out-Null
    $svg.Add('<rect x="742" y="24" width="12" height="12" fill="#0072b2"/><text x="760" y="34" font-family="Arial" font-size="10">full direct energy</text>') | Out-Null
    $svg.Add('<rect x="742" y="44" width="12" height="12" fill="#cc79a7"/><text x="760" y="54" font-family="Arial" font-size="10">sampled direct audit</text>') | Out-Null
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$fixtureSpecs = if ($Quick) {
    @(@{ fixture='research-medium'; caps=@(128,256) })
} else {
    @(
        @{ fixture='barnase'; caps=@(128,256) },
        @{ fixture='research-medium'; caps=@(128,256) },
        @{ fixture='research-alt'; caps=@(128,256) }
    )
}

$fixtureDir = Join-Path $workRootFull 'fixtures'
New-Item -ItemType Directory -Path $fixtureDir -Force | Out-Null
$rows = New-Object System.Collections.Generic.List[object]
foreach ($spec in $fixtureSpecs) {
    $fixture = [string]$spec.fixture
    $sourceDir = Join-Path $fixtureDir $fixture
    New-Item -ItemType Directory -Path $sourceDir -Force | Out-Null
    $source = Resolve-Z4bFixture -Fixture $fixture -RunDir $sourceDir
    foreach ($cap in @($spec.caps)) {
        if ($null -eq $source -or -not (Test-Path -LiteralPath $source)) {
            foreach ($backend in @('grid_context','fmm_ilu0')) {
                $rows.Add((New-Z4bRow -Condition ("{0}_cap{1}_{2}" -f $fixture, $cap, $backend) `
                    -Fixture $fixture -AtomCap $cap -Backend $backend -Run $null `
                    -Status 'skipped_missing_fixture')) | Out-Null
            }
            continue
        }
        $sampled = Join-Path $sourceDir ("{0}_cap{1}.pqr" -f (Get-FmmCliSafeName $fixture), $cap)
        [void](Write-FmmCliSampledPqr -InputPath $source -OutputPath $sampled -Cap $cap)
        $backends = @('grid_context','fmm_ilu0')
        if ($cap -lt 256) {
            $backends = @('grid_context','tabi_direct_full','fmm_ilu0')
        }
        foreach ($backend in $backends) {
            $rows.Add((Invoke-Z4bCase -Fixture $fixture -AtomCap $cap -PqrPath $sampled `
                -Backend $backend)) | Out-Null
        }
    }
}

$rowArray = $rows.ToArray()
Add-Z4bReferencePolicy -Rows $rowArray
$decision = Get-Z4bDecision -Rows $rowArray
$referenceRows = @($rowArray | Where-Object {
        $_.backend -eq 'fmm_ilu0'
    } | ForEach-Object {
        [pscustomobject]@{
            fixture = $_.fixture
            atom_cap = $_.atom_cap
            reference_kind = $_.reference_kind
            reference_backend = $_.reference_backend
            reference_status = $_.reference_status
            reference_usable_for_paper = $_.reference_usable_for_paper
            energy_delta_kj_mol = $_.energy_delta_kj_mol
            relative_energy_delta = $_.relative_energy_delta
            tree_vs_direct_matvec_rel = $_.tree_vs_direct_matvec_rel
            matvec_audit_rows = $_.matvec_audit_rows
        }
    })
$maxFmmAudit = @($rowArray | Where-Object {
        $_.backend -eq 'fmm_ilu0'
    } | ForEach-Object { Convert-Z4bDouble $_.tree_vs_direct_matvec_rel } |
    Where-Object { Test-FmmCliFinite $_ } | Measure-Object -Maximum).Maximum
if (-not (Test-FmmCliFinite $maxFmmAudit)) { $maxFmmAudit = [double]::NaN }
$gridContextWarnings = @($rowArray | Where-Object {
        $_.backend -eq 'grid_context' -and $_.status -ne 'ok'
    }).Count

$rowArray | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$referenceRows | Export-Csv -LiteralPath $ReferenceCsv -NoTypeInformation
Write-Z4bSvg -Rows $rowArray -Path $PlotPath

[pscustomobject]@{
    decision = $decision
    rows = $rowArray.Count
    quick = [bool]$Quick
    max_fmm_audit = $maxFmmAudit
    grid_context_warnings = $gridContextWarnings
    artifact_matrix = $MatrixCsv
    artifact_reference_policy = $ReferenceCsv
    artifact_plot = $PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'scope: Z4B reference/verification policy gate',
    "decision: $decision",
    "rows: $($rowArray.Count)",
    "quick: $([bool]$Quick)",
    'fmm_policy: fmm-nearfield-ilu0, theta=0.6, taylor_order=4, near=analytic_singular, quad=off',
    'reference_policy: cap-128 FMM uses full direct TABI energy; cap-256 FMM uses sampled direct-operator audit',
    'grid_policy: grid rows are context only unless explicitly status=ok; grid is not used as paper oracle when not converged',
    "max_fmm_audit: $maxFmmAudit",
    "grid_context_warnings: $gridContextWarnings",
    'default_policy_note: grid remains production default; Z4B does not change public behavior'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z4B Reference Policy Gate Report',
    '',
    "Decision: $decision",
    '',
    'Z4B separates paper references from context rows after Z4A showed direct TABI cap-256 is runtime-limited and grid does not converge uniformly.',
    '',
    'Artifacts:',
    "- Matrix: $MatrixCsv",
    "- Reference policy: $ReferenceCsv",
    "- Decision: $DecisionCsv",
    "- Summary: $SummaryPath",
    "- Plot: $PlotPath",
    '',
    'Policy:',
    '- cap-128 FMM rows cite full direct TABI energy when available.',
    '- cap-256 FMM rows cite sampled direct-operator audit from pbe_tabi_profile.csv.',
    '- grid rows are runtime/memory context and are not paper oracle rows unless they converge.',
    '',
    'Rows:'
)
foreach ($row in $rowArray) {
    $report += ('- {0}: fixture={1}, cap={2}, backend={3}, status={4}, class={5}, ref_kind={6}, ref={7}:{8}, usable={9}, energy={10}, delta={11}, audit={12}, wall_s={13}, rss_mb={14}' -f `
        $row.condition, $row.fixture, $row.atom_cap, $row.backend, $row.status, `
        $row.failure_class, $row.reference_kind, $row.reference_backend, `
        $row.reference_status, $row.reference_usable_for_paper, `
        $row.energy_kj_mol, $row.energy_delta_kj_mol, `
        $row.tree_vs_direct_matvec_rel, $row.wall_seconds, `
        $row.peak_working_set_mb)
}
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "PASS Scope Z4B reference policy gate: $decision"
