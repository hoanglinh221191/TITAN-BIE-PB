param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [switch]$IncludeResearchLarge,
    [switch]$IncludeDirectCap256,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 2,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4a')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4a_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4a_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4a_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4a_runtime_accuracy_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4A_PROTEIN_PANEL_BENCHMARK_REPORT.txt'

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
    throw "Refusing to clean Z4A work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z4aDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or
        $text -eq 'nan') {
        return $Default
    }
    return [double]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Convert-Z4aInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or
        $text -eq 'nan') {
        return $Default
    }
    return [int64]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Format-Z4aDouble([double]$Value) {
    return $Value.ToString('G15', [System.Globalization.CultureInfo]::InvariantCulture)
}

function Get-Z4aMedian([double[]]$Values) {
    $finite = @($Values | Where-Object { Test-FmmCliFinite $_ } | Sort-Object)
    if ($finite.Count -eq 0) { return [double]::NaN }
    $mid = [int]($finite.Count / 2)
    if (($finite.Count % 2) -eq 1) { return [double]$finite[$mid] }
    return 0.5 * ([double]$finite[$mid - 1] + [double]$finite[$mid])
}

function New-Z4aRow(
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
        reference_backend = ''
        reference_energy_kj_mol = [double]::NaN
        energy_delta_kj_mol = [double]::NaN
        relative_energy_delta = [double]::NaN
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

function Resolve-Z4aFixture([string]$Fixture, [string]$RunDir) {
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

function Write-Z4aArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [ValidateSet('grid','tabi_direct','fmm_ilu0')][string]$Backend
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
    if ($Backend -eq 'grid') {
        $argLines += '-pbe=grid'
    } else {
        $mode = if ($Backend -eq 'tabi_direct') { 'direct' } else { 'treecode' }
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

function Read-Z4aResult([string]$RunDir, [string]$Backend) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    if (-not (Test-Path -LiteralPath $energyPath)) {
        throw "Missing pbe_energy.csv in $RunDir"
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $result = [ordered]@{
        status = if ($energy.status -eq 'ok') { 'ok' } else { [string]$energy.status }
        panels = 0
        iterations = 0
        energy_kj_mol = Convert-Z4aDouble $energy.E_rxn_kJmol
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
    if ($Backend -ne 'grid' -and (Test-Path -LiteralPath $profilePath)) {
        $profile = @(Import-Csv -LiteralPath $profilePath)[0]
        $result.status = if ($result.status -eq 'ok' -and $profile.status -eq 'ok') { 'ok' } else { "energy=$($energy.status);profile=$($profile.status)" }
        $result.panels = Convert-Z4aInt64 $profile.panels
        $result.iterations = Convert-Z4aInt64 $profile.iterations
        $result.true_residual = Convert-Z4aDouble $profile.true_residual_rel
        $result.tree_vs_direct_matvec_rel = Convert-Z4aDouble $profile.tree_vs_direct_matvec_rel
        $result.matvec_audit_rows = Convert-Z4aInt64 $profile.matvec_audit_rows
        $result.matvec_seconds_avg = Convert-Z4aDouble $profile.matvec_seconds_avg
        $result.preconditioner_effective = [string]$profile.preconditioner_effective
        $result.near_integration = [string]$profile.near_integration
        $result.quad_buffer = [string]$profile.quad_buffer
        $result.fmm_m2l_mode = [string]$profile.fmm_m2l_mode
        $result.fmm_m2l_expansion_order = Convert-Z4aInt64 $profile.fmm_m2l_expansion_order
        $result.fmm_near_pairs = Convert-Z4aInt64 $profile.fmm_m2l_near_pairs
        $result.fmm_far_pairs = Convert-Z4aInt64 $profile.fmm_m2l_far_pairs
        $result.fmm_finite_outputs = [string]$profile.fmm_m2l_finite_outputs
    }
    return [pscustomobject]$result
}

function Invoke-Z4aCase(
    [string]$Fixture,
    [int]$AtomCap,
    [string]$PqrPath,
    [ValidateSet('grid','tabi_direct','fmm_ilu0')][string]$Backend
) {
    $condition = "{0}_cap{1}_{2}" -f (Get-FmmCliSafeName $Fixture), $AtomCap, $Backend
    $runDir = Join-Path $workRootFull $condition
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $argPath = Write-Z4aArgFile -RunDir $runDir -PqrPath $PqrPath -Backend $Backend
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        return New-Z4aRow -Condition $condition -Fixture $Fixture `
            -AtomCap $AtomCap -Backend $Backend -Run $run `
            -Status $(if ($run.timed_out) { 'timeout_not_verified' } else { "exit_$($run.exit_code)" })
    }
    try {
        $result = Read-Z4aResult -RunDir $runDir -Backend $Backend
    } catch {
        return New-Z4aRow -Condition $condition -Fixture $Fixture `
            -AtomCap $AtomCap -Backend $Backend -Run $run `
            -Status "profile_error: $($_.Exception.Message)"
    }
    $failureClass = 'pass'
    if ($result.status -ne 'ok') {
        $failureClass = 'status_fail'
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
        reference_backend = ''
        reference_energy_kj_mol = [double]::NaN
        energy_delta_kj_mol = [double]::NaN
        relative_energy_delta = [double]::NaN
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

function Add-Z4aReferences([object[]]$Rows) {
    foreach ($group in $Rows | Group-Object { '{0}|{1}' -f $_.fixture, $_.atom_cap }) {
        $directRows = @($group.Group | Where-Object { $_.backend -eq 'tabi_direct' -and $_.status -eq 'ok' })
        $gridRows = @($group.Group | Where-Object { $_.backend -eq 'grid' -and $_.status -eq 'ok' })
        $direct = if ($directRows.Count -gt 0) { $directRows[0] } else { $null }
        $grid = if ($gridRows.Count -gt 0) { $gridRows[0] } else { $null }
        foreach ($row in $group.Group) {
            $ref = if ($null -ne $direct) { $direct } elseif ($null -ne $grid) { $grid } else { $null }
            if ($null -eq $ref -or $row.backend -eq $ref.backend) { continue }
            $row.reference_backend = $ref.backend
            $row.reference_energy_kj_mol = $ref.energy_kj_mol
            if ((Test-FmmCliFinite $row.energy_kj_mol) -and (Test-FmmCliFinite $ref.energy_kj_mol)) {
                $row.energy_delta_kj_mol = $row.energy_kj_mol - $ref.energy_kj_mol
                $row.relative_energy_delta = [Math]::Abs($row.energy_delta_kj_mol) /
                    [Math]::Max(1.0, [Math]::Abs($ref.energy_kj_mol))
            }
        }
    }
}

function Get-Z4aDecision([object[]]$Rows) {
    $requiredFmm = @($Rows | Where-Object { $_.backend -eq 'fmm_ilu0' })
    if (@($Rows | Where-Object { $_.status -eq 'skipped_missing_fixture' }).Count -gt 0) {
        return 'z4a_fixture_limited'
    }
    if (@($requiredFmm | Where-Object { $_.timed_out -or $_.status -eq 'timeout_not_verified' }).Count -gt 0) {
        return 'z4a_fixture_limited'
    }
    if (@($requiredFmm | Where-Object { $_.failure_class -in @('accuracy_fail','telemetry_fail','status_fail') }).Count -gt 0) {
        return 'z4a_accuracy_gap'
    }
    $directRows = @($Rows | Where-Object { $_.backend -eq 'tabi_direct' })
    if (@($directRows | Where-Object { $_.status -eq 'direct_skipped_cap256_budget' -or $_.status -eq 'timeout_not_verified' }).Count -gt 0) {
        return 'z4a_direct_limited_fmm_pass'
    }
    if (@($Rows | Where-Object { $_.backend -eq 'fmm_ilu0' -and $_.reference_backend -eq 'tabi_direct' -and (Test-FmmCliFinite $_.relative_energy_delta) -and $_.relative_energy_delta -gt 0.05 }).Count -gt 0) {
        return 'z4a_accuracy_gap'
    }
    return 'z4a_protein_panel_pass'
}

function Write-Z4aSvg([object[]]$Rows, [string]$Path) {
    $points = @($Rows | Where-Object { $_.status -eq 'ok' })
    $maxWall = ($points | ForEach-Object { Convert-Z4aDouble $_.wall_seconds } |
        Where-Object { Test-FmmCliFinite $_ } | Measure-Object -Maximum).Maximum
    if (-not (Test-FmmCliFinite $maxWall) -or $maxWall -le 0.0) {
        $maxWall = 1.0
    }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="980" height="360" viewBox="0 0 980 360">') | Out-Null
    $svg.Add('<rect width="980" height="360" fill="white"/>') | Out-Null
    $svg.Add('<text x="24" y="28" font-family="Arial" font-size="16">Scope Z4A protein panel benchmark</text>') | Out-Null
    $svg.Add('<text x="24" y="48" font-family="Arial" font-size="11">Bar height = wall seconds; color = backend</text>') | Out-Null
    $x = 34
    foreach ($row in $points) {
        $wall = Convert-Z4aDouble $row.wall_seconds 0.0
        $h = [Math]::Max(2.0, 230.0 * ($wall / [Math]::Max($maxWall, 1.0e-12)))
        $y = 292.0 - $h
        $fill = switch ($row.backend) {
            'grid' { '#009e73' }
            'tabi_direct' { '#d55e00' }
            default { '#0072b2' }
        }
        $svg.Add(('<rect x="{0}" y="{1:F1}" width="22" height="{2:F1}" fill="{3}"/>' -f $x, $y, $h, $fill)) | Out-Null
        $label = ('{0} c{1} {2}' -f $row.fixture, $row.atom_cap, $row.backend)
        $svg.Add(('<text x="{0}" y="318" font-family="Arial" font-size="8" transform="rotate(45 {0},318)">{1}</text>' -f $x, $label)) | Out-Null
        $x += 30
        if ($x -gt 940) { break }
    }
    $svg.Add('<rect x="800" y="20" width="12" height="12" fill="#009e73"/><text x="818" y="30" font-family="Arial" font-size="10">grid</text>') | Out-Null
    $svg.Add('<rect x="800" y="38" width="12" height="12" fill="#d55e00"/><text x="818" y="48" font-family="Arial" font-size="10">tabi direct</text>') | Out-Null
    $svg.Add('<rect x="800" y="56" width="12" height="12" fill="#0072b2"/><text x="818" y="66" font-family="Arial" font-size="10">FMM ILU0</text>') | Out-Null
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$fixtureSpecs = if ($Quick) {
    @(@{ fixture='research-medium'; caps=@(128) })
} else {
    $list = @(
        @{ fixture='barnase'; caps=@(128,256) },
        @{ fixture='research-medium'; caps=@(128,256) },
        @{ fixture='research-alt'; caps=@(128,256) }
    )
    if ($IncludeResearchLarge) {
        $list += @{ fixture='research-large'; caps=@(128,256) }
    }
    $list
}

$fixtureDir = Join-Path $workRootFull 'fixtures'
New-Item -ItemType Directory -Path $fixtureDir -Force | Out-Null
$rows = New-Object System.Collections.Generic.List[object]
foreach ($spec in $fixtureSpecs) {
    $fixture = [string]$spec.fixture
    $sourceDir = Join-Path $fixtureDir $fixture
    New-Item -ItemType Directory -Path $sourceDir -Force | Out-Null
    $source = Resolve-Z4aFixture -Fixture $fixture -RunDir $sourceDir
    foreach ($cap in @($spec.caps)) {
        if ($null -eq $source -or -not (Test-Path -LiteralPath $source)) {
            foreach ($backend in @('grid','fmm_ilu0')) {
                $rows.Add((New-Z4aRow -Condition ("{0}_cap{1}_{2}" -f $fixture, $cap, $backend) `
                    -Fixture $fixture -AtomCap $cap -Backend $backend -Run $null `
                    -Status 'skipped_missing_fixture')) | Out-Null
            }
            continue
        }
        $sampled = Join-Path $sourceDir ("{0}_cap{1}.pqr" -f (Get-FmmCliSafeName $fixture), $cap)
        [void](Write-FmmCliSampledPqr -InputPath $source -OutputPath $sampled -Cap $cap)
        $backends = if ($Quick) { @('grid','fmm_ilu0') } else { @('grid','tabi_direct','fmm_ilu0') }
        foreach ($backend in $backends) {
            if ($backend -eq 'tabi_direct' -and $cap -ge 256 -and -not $IncludeDirectCap256) {
                $rows.Add((New-Z4aRow -Condition ("{0}_cap{1}_{2}" -f (Get-FmmCliSafeName $fixture), $cap, $backend) `
                    -Fixture $fixture -AtomCap $cap -Backend $backend -Run $null `
                    -Status 'direct_skipped_cap256_budget')) | Out-Null
                continue
            }
            $rows.Add((Invoke-Z4aCase -Fixture $fixture -AtomCap $cap -PqrPath $sampled `
                -Backend $backend)) | Out-Null
        }
    }
}

$rowArray = $rows.ToArray()
Add-Z4aReferences -Rows $rowArray
$decision = Get-Z4aDecision -Rows $rowArray
$medianFmmWall = Get-Z4aMedian @($rowArray | Where-Object {
        $_.backend -eq 'fmm_ilu0' -and $_.status -eq 'ok'
    } | ForEach-Object { Convert-Z4aDouble $_.wall_seconds })
$maxFmmAudit = @($rowArray | Where-Object {
        $_.backend -eq 'fmm_ilu0'
    } | ForEach-Object { Convert-Z4aDouble $_.tree_vs_direct_matvec_rel } |
    Where-Object { Test-FmmCliFinite $_ } | Measure-Object -Maximum).Maximum
if (-not (Test-FmmCliFinite $maxFmmAudit)) { $maxFmmAudit = [double]::NaN }

$rowArray | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
Write-Z4aSvg -Rows $rowArray -Path $PlotPath

[pscustomobject]@{
    decision = $decision
    rows = $rowArray.Count
    quick = [bool]$Quick
    include_research_large = [bool]$IncludeResearchLarge
    include_direct_cap256 = [bool]$IncludeDirectCap256
    median_fmm_wall_seconds = $medianFmmWall
    max_fmm_audit = $maxFmmAudit
    artifact_matrix = $MatrixCsv
    artifact_plot = $PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'scope: Z4A protein panel benchmark',
    "decision: $decision",
    "rows: $($rowArray.Count)",
    "quick: $([bool]$Quick)",
    "include_research_large: $([bool]$IncludeResearchLarge)",
    "include_direct_cap256: $([bool]$IncludeDirectCap256)",
    'fmm_policy: fmm-nearfield-ilu0, theta=0.6, taylor_order=4, near=analytic_singular, quad=off',
    'direct_policy: cap-256 direct rows are classified as skipped unless -IncludeDirectCap256 is set',
    "median_fmm_wall_seconds: $medianFmmWall",
    "max_fmm_audit: $maxFmmAudit",
    'default_policy_note: grid remains production default; Z4A does not change public behavior'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z4A Protein Panel Benchmark Report',
    '',
    "Decision: $decision",
    '',
    'Z4A runs a paper-style sampled protein panel using the Z3D-confirmed FMM/M2L policy: ILU0, theta 0.6, Taylor order 4, analytic-singular near integration, and quad buffer off.',
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
    $report += ('- {0}: fixture={1}, cap={2}, backend={3}, status={4}, class={5}, energy={6}, ref={7}:{8}, delta={9}, audit={10}, matvec_s={11}, wall_s={12}, rss_mb={13}' -f `
        $row.condition, $row.fixture, $row.atom_cap, $row.backend, $row.status, `
        $row.failure_class, $row.energy_kj_mol, $row.reference_backend, `
        $row.reference_energy_kj_mol, $row.energy_delta_kj_mol, `
        $row.tree_vs_direct_matvec_rel, $row.matvec_seconds_avg, `
        $row.wall_seconds, $row.peak_working_set_mb)
}
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "PASS Scope Z4A protein panel benchmark: $decision"
