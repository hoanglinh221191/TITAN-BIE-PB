param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$PrimaryFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$ExpansionFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public_z4q'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4z_apbs'),
    [string[]]$Fixtures = @(),
    [string[]]$Dimes = @(),
    [string[]]$FineScales = @(),
    [string]$ApbsExe = '',
    [int]$ApbsTimeoutSeconds = 900,
    [int]$FmmTimeoutSeconds = 1800,
    [int]$Threads = 4,
    [double]$ApbsStabilityRelativeThreshold = 0.05,
    [switch]$SkipTitan,
    [switch]$Quick,
    [switch]$PlanOnly,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$PrimaryFixtureRoot = [IO.Path]::GetFullPath($PrimaryFixtureRoot)
$ExpansionFixtureRoot = [IO.Path]::GetFullPath($ExpansionFixtureRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4z_apbs_matrix.csv'
$StabilityCsv = Join-Path $WorkRoot 'fagbi_scope_z4z_apbs_stability.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4z_apbs_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4z_apbs_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4z_apbs_grid_plot.svg'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4Z_APBS_GRID_SWEEP_REPORT.txt'

$TemperatureK = 298.15
$InternalDielectric = 1.0
$ExternalDielectric = 80.0
$SaltM = 0.15

if ($Dimes.Count -eq 0) {
    $Dimes = if ($Quick) { @(97,129) } else { @(97,129,161,193) }
}
if ($FineScales.Count -eq 0) {
    $FineScales = @(0.55)
}
if ($Fixtures.Count -eq 0) {
    $Fixtures = if ($Quick) { @('1vii') } else { @('1vii','1ubq','2ci2') }
}
$Fixtures = @($Fixtures | ForEach-Object {
    ([string]$_ -split ',') | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne '' }
})
$parsedDimes = New-Object System.Collections.Generic.List[int]
foreach ($token in @($Dimes | ForEach-Object { ([string]$_ -split ',') })) {
    $text = ([string]$token).Trim()
    if ($text -eq '') { continue }
    $value = 0
    if (-not [int]::TryParse($text, [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$value)) {
        throw "Invalid APBS dime '$text'. Use values such as -Dimes 225,257 or -Dimes 225 257."
    }
    $parsedDimes.Add($value) | Out-Null
}
$Dimes = @($parsedDimes.ToArray())
$parsedFineScales = New-Object System.Collections.Generic.List[double]
foreach ($token in @($FineScales | ForEach-Object { ([string]$_ -split ',') })) {
    $text = ([string]$token).Trim()
    if ($text -eq '') { continue }
    $value = 0.0
    if (-not [double]::TryParse($text, [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$value)) {
        throw "Invalid APBS fine-scale '$text'. Use values such as -FineScales 0.55,0.7."
    }
    $parsedFineScales.Add($value) | Out-Null
}
$FineScales = @($parsedFineScales.ToArray())

if (-not $SkipTitan) {
    if (-not (Test-Path -LiteralPath $TitanExe)) {
        throw "Missing titan executable: $TitanExe"
    }
    $TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
}

$RuntimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
        $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
    }
    $RuntimeRepoRoot = Split-Path -Parent $commonGitDir
}
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    throw "Missing lib_ff under runtime repo root: $RuntimeRepoRoot"
}

New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

function Get-Z4zDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or $text -eq 'nan') { return $Default }
    $parsed = 0.0
    if ([double]::TryParse($text, [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Format-Z4zDouble([double]$Value, [string]$Format = 'G6') {
    if (-not (Test-FmmCliFinite $Value)) { return 'NA' }
    return $Value.ToString($Format, [Globalization.CultureInfo]::InvariantCulture)
}

function Get-Z4zObjectProperty([object]$Object, [string]$Name, [object]$Default = '') {
    if ($null -eq $Object) { return $Default }
    $prop = $Object.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Import-Z4zFixtureRows([string]$Root, [string]$ManifestName) {
    $path = Join-Path $Root $ManifestName
    if (-not (Test-Path -LiteralPath $path)) { return @() }
    return @(Import-Csv -LiteralPath $path | ForEach-Object {
        $_ | Add-Member -NotePropertyName fixture_root -NotePropertyValue $Root -Force
        $_
    })
}

function Get-Z4zFixtureRows {
    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($row in (Import-Z4zFixtureRows $PrimaryFixtureRoot 'public_fixture_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    foreach ($row in (Import-Z4zFixtureRows $ExpansionFixtureRoot 'public_fixture_expansion_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    return $rows.ToArray()
}

function Get-Z4zFixture([string]$Name, [object[]]$Rows) {
    $match = @($Rows | Where-Object { $_.fixture -eq $Name })
    if ($match.Count -ne 1) {
        throw "Expected exactly one public fixture row for '$Name', found $($match.Count)."
    }
    $row = $match[0]
    $pqrPath = Join-Path $row.fixture_root $row.pqr
    if (-not (Test-Path -LiteralPath $pqrPath)) {
        throw "Missing cached PQR for $Name`: $pqrPath"
    }
    $row | Add-Member -NotePropertyName pqr_path -NotePropertyValue $pqrPath -Force
    return $row
}

function Resolve-Z4zApbsExe([string]$Requested) {
    if (-not [string]::IsNullOrWhiteSpace($Requested)) {
        if (Test-Path -LiteralPath $Requested -PathType Leaf) {
            return (Resolve-Path -LiteralPath $Requested).Path
        }
        $cmd = Get-Command $Requested -ErrorAction SilentlyContinue
        if ($null -ne $cmd) { return $cmd.Source }
        return ''
    }
    foreach ($name in @('apbs.exe','apbs')) {
        $cmd = Get-Command $name -ErrorAction SilentlyContinue
        if ($null -ne $cmd) { return $cmd.Source }
    }
    $known = @(
        (Join-Path $env:USERPROFILE '.conda\envs\titan_apbs\Library\bin\apbs.exe'),
        'C:\ProgramData\miniconda3\envs\titan_apbs\Library\bin\apbs.exe'
    )
    foreach ($candidate in $known) {
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }
    return ''
}

function Get-Z4zPqrAtoms([string]$PqrPath) {
    $atoms = New-Object System.Collections.Generic.List[object]
    foreach ($line in Get-Content -LiteralPath $PqrPath) {
        if ($line -notmatch '^(ATOM|HETATM)') { continue }
        $parts = @($line -split '\s+' | Where-Object { $_ -ne '' })
        if ($parts.Count -lt 10) { continue }
        $atoms.Add([pscustomobject]@{
            x = Get-Z4zDouble $parts[5]
            y = Get-Z4zDouble $parts[6]
            z = Get-Z4zDouble $parts[7]
            charge = Get-Z4zDouble $parts[8] 0.0
            radius = Get-Z4zDouble $parts[9] 1.5
        }) | Out-Null
    }
    return $atoms.ToArray()
}

function Get-Z4zPqrBox([string]$PqrPath, [double]$Padding) {
    $atoms = @(Get-Z4zPqrAtoms -PqrPath $PqrPath)
    if ($atoms.Count -eq 0) { throw "Cannot derive APBS box from PQR: $PqrPath" }
    $xmin = ($atoms | Measure-Object x -Minimum).Minimum
    $xmax = ($atoms | Measure-Object x -Maximum).Maximum
    $ymin = ($atoms | Measure-Object y -Minimum).Minimum
    $ymax = ($atoms | Measure-Object y -Maximum).Maximum
    $zmin = ($atoms | Measure-Object z -Minimum).Minimum
    $zmax = ($atoms | Measure-Object z -Maximum).Maximum
    return [pscustomobject]@{
        rx = [Math]::Max(20.0, ($xmax - $xmin) + $Padding)
        ry = [Math]::Max(20.0, ($ymax - $ymin) + $Padding)
        rz = [Math]::Max(20.0, ($zmax - $zmin) + $Padding)
        atom_count = $atoms.Count
    }
}

function Write-Z4zApbsInput(
    [string]$PqrPath,
    [string]$InputPath,
    [int]$Dime,
    [double]$FineScale,
    [double]$CoarsePadding
) {
    $box = Get-Z4zPqrBox -PqrPath $PqrPath -Padding $CoarsePadding
    $fmt = { param([double]$v) $v.ToString('G12', [Globalization.CultureInfo]::InvariantCulture) }
    $pqrName = Split-Path -Leaf $PqrPath
    $coarseLen = "    cglen $(& $fmt $box.rx) $(& $fmt $box.ry) $(& $fmt $box.rz)"
    $fineLen = "    fglen $(& $fmt ([Math]::Max(12.0, $box.rx * $FineScale))) $(& $fmt ([Math]::Max(12.0, $box.ry * $FineScale))) $(& $fmt ([Math]::Max(12.0, $box.rz * $FineScale)))"
    $dimeLine = "    dime $Dime $Dime $Dime"
    $commonSolvated = @(
        '    mg-auto',
        $dimeLine,
        $coarseLen,
        $fineLen,
        '    cgcent mol 1',
        '    fgcent mol 1',
        '    mol 1',
        '    lpbe',
        '    bcfl sdh',
        '    ion charge 1 conc 0.15 radius 2.0',
        '    ion charge -1 conc 0.15 radius 2.0',
        '    pdie 1.0',
        '    sdie 80.0',
        '    srfm smol',
        '    chgm spl2',
        '    sdens 10.0',
        '    srad 1.4',
        '    swin 0.3',
        '    temp 298.15',
        '    calcenergy total',
        '    calcforce no'
    )
    $commonReference = @(
        '    mg-auto',
        $dimeLine,
        $coarseLen,
        $fineLen,
        '    cgcent mol 1',
        '    fgcent mol 1',
        '    mol 1',
        '    lpbe',
        '    bcfl sdh',
        '    pdie 1.0',
        '    sdie 1.0',
        '    srfm smol',
        '    chgm spl2',
        '    sdens 10.0',
        '    srad 1.4',
        '    swin 0.3',
        '    temp 298.15',
        '    calcenergy total',
        '    calcforce no'
    )
    $lines = @(
        'read',
        "    mol pqr $pqrName",
        'end',
        'elec name solvated'
    ) + $commonSolvated + @(
        'end',
        'elec name reference'
    ) + $commonReference + @(
        'end',
        'print elecEnergy solvated - reference end',
        'quit'
    )
    Set-Content -LiteralPath $InputPath -Value $lines -Encoding ASCII
    return $box
}

function Invoke-Z4zExternalCase([string]$Exe, [string]$RunDir, [string]$Arguments, [int]$TimeoutSeconds) {
    $stdoutPath = Join-Path $RunDir 'stdout.log'
    $stderrPath = Join-Path $RunDir 'stderr.log'
    $psi = New-Object System.Diagnostics.ProcessStartInfo
    $psi.FileName = $Exe
    $psi.WorkingDirectory = $RunDir
    $psi.UseShellExecute = $false
    $psi.CreateNoWindow = $true
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.Arguments = $Arguments
    $process = New-Object System.Diagnostics.Process
    $process.StartInfo = $psi
    $stopwatch = [Diagnostics.Stopwatch]::StartNew()
    [void]$process.Start()
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    $deadline = if ($TimeoutSeconds -gt 0) { [DateTime]::UtcNow.AddSeconds($TimeoutSeconds) } else { [DateTime]::MaxValue }
    $peak = 0L
    while (-not $process.HasExited -and [DateTime]::UtcNow -lt $deadline) {
        try {
            $process.Refresh()
            $peak = [Math]::Max($peak, [int64]$process.WorkingSet64)
            $peak = [Math]::Max($peak, [int64]$process.PeakWorkingSet64)
        } catch {}
        Start-Sleep -Milliseconds 100
    }
    $timedOut = -not $process.HasExited
    if ($timedOut) {
        try { $process.Kill() } catch {}
        $process.WaitForExit()
    }
    try {
        $process.Refresh()
        $peak = [Math]::Max($peak, [int64]$process.PeakWorkingSet64)
    } catch {}
    $stdoutTask.Wait()
    $stderrTask.Wait()
    $stopwatch.Stop()
    Set-Content -LiteralPath $stdoutPath -Value $stdoutTask.Result -Encoding UTF8
    Set-Content -LiteralPath $stderrPath -Value $stderrTask.Result -Encoding UTF8
    return [pscustomobject]@{
        exit_code = if ($timedOut) { 124 } else { $process.ExitCode }
        timed_out = $timedOut
        wall_seconds = $stopwatch.Elapsed.TotalSeconds
        peak_working_set_mb = $peak / 1MB
    }
}

function Read-Z4zApbsEnergy([string]$RunDir) {
    $texts = @()
    foreach ($name in @('stdout.log','stderr.log')) {
        $path = Join-Path $RunDir $name
        if (Test-Path -LiteralPath $path) { $texts += (Get-Content -LiteralPath $path -Raw) }
    }
    $text = $texts -join "`n"
    $patterns = @(
        'Global net ELEC energy\s*=\s*([-+0-9.eE]+)',
        'Total electrostatic energy\s*=\s*([-+0-9.eE]+)',
        'elecEnergy\s*=\s*([-+0-9.eE]+)'
    )
    foreach ($pattern in $patterns) {
        $m = [regex]::Match($text, $pattern, [Text.RegularExpressions.RegexOptions]::IgnoreCase)
        if ($m.Success) {
            return [double]::Parse($m.Groups[1].Value, [Globalization.CultureInfo]::InvariantCulture)
        }
    }
    return [double]::NaN
}

function New-Z4zRow([object]$Case, [string]$RunDir) {
    return [pscustomobject]@{
        case_id = $Case.case_id
        fixture = $Case.fixture
        protein_name = $Case.protein_name
        pdb_id = $Case.pdb_id
        backend = $Case.backend
        dime = $Case.dime
        fine_scale = $Case.fine_scale
        surface_points = $Case.surface_points
        timeout_seconds = $Case.timeout_seconds
        status = 'planned'
        failure_class = 'planned'
        exit_code = ''
        timed_out = ''
        atom_count = $Case.atom_count
        panels = ''
        iterations = ''
        energy_kj_mol = ''
        titan_energy_kj_mol = ''
        apbs_energy_kj_mol = ''
        titan_minus_apbs = ''
        sampled_audit_rms = ''
        sampled_direct_true_residual = ''
        wall_seconds = ''
        peak_rss_mb = ''
        cglen_x = ''
        cglen_y = ''
        cglen_z = ''
        fglen_x = ''
        fglen_y = ''
        fglen_z = ''
        apbs_exe = ''
        run_dir = $RunDir
    }
}

function New-Z4zTitanCase([object]$Fixture) {
    return [pscustomobject]@{
        case_id = 'titan_guarded_p4_{0}' -f $Fixture.fixture
        fixture = $Fixture.fixture
        protein_name = $Fixture.protein_name
        pdb_id = $Fixture.pdb_id
        backend = 'titan_guarded_p4'
        dime = 0
        fine_scale = 0
        surface_points = 8
        timeout_seconds = $FmmTimeoutSeconds
        pqr_source = $Fixture.pqr_path
        atom_count = [int](Get-Z4zObjectProperty -Object $Fixture -Name 'atom_count' -Default 0)
    }
}

function New-Z4zApbsCase([object]$Fixture, [int]$Dime, [double]$FineScale) {
    $scaleText = $FineScale.ToString('0.###', [Globalization.CultureInfo]::InvariantCulture).Replace('.','p')
    return [pscustomobject]@{
        case_id = 'apbs_{0}_d{1}_fg{2}' -f $Fixture.fixture,$Dime,$scaleText
        fixture = $Fixture.fixture
        protein_name = $Fixture.protein_name
        pdb_id = $Fixture.pdb_id
        backend = 'apbs'
        dime = $Dime
        fine_scale = $FineScale
        surface_points = 0
        timeout_seconds = $ApbsTimeoutSeconds
        pqr_source = $Fixture.pqr_path
        atom_count = [int](Get-Z4zObjectProperty -Object $Fixture -Name 'atom_count' -Default 0)
    }
}

function Upsert-Z4zRow([System.Collections.Generic.List[object]]$Rows, [object]$Row) {
    for ($i = 0; $i -lt $Rows.Count; ++$i) {
        if ($Rows[$i].case_id -eq $Row.case_id) {
            $Rows[$i] = $Row
            return
        }
    }
    $Rows.Add($Row) | Out-Null
}

function Save-Z4zRows([object[]]$Rows) {
    $Rows | Sort-Object fixture,backend,dime,fine_scale | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
}

function Should-RunZ4zCase([object]$Case, [object[]]$Rows) {
    if ($Force) { return $true }
    $old = @($Rows | Where-Object { $_.case_id -eq $Case.case_id })
    if ($old.Count -eq 0) { return $true }
    if ($RerunFailed) {
        $status = [string]$old[-1].status
        $failure = [string]$old[-1].failure_class
        return -not ($status -eq 'ok' -or $failure -eq 'pass' -or $failure -like 'skipped_*')
    }
    return $false
}

function Copy-Z4zPqr([object]$Case, [string]$RunDir) {
    $pqr = Join-Path $RunDir (Split-Path -Leaf $Case.pqr_source)
    Copy-Item -LiteralPath $Case.pqr_source -Destination $pqr -Force
    return $pqr
}

function Get-Z4zProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Read-Z4zTitanOutputs([string]$RunDir, [object]$Row) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    if (-not (Test-Path -LiteralPath $energyPath)) {
        $Row.status = 'missing_output'
        $Row.failure_class = 'missing_energy_output'
        return $Row
    }
    $energyRows = @(Import-Csv -LiteralPath $energyPath)
    if ($energyRows.Count -eq 0) {
        $Row.status = 'missing_output'
        $Row.failure_class = 'missing_energy_output'
        return $Row
    }
    $energy = $energyRows[0]
    $Row.energy_kj_mol = Get-Z4zDouble $energy.E_rxn_kJmol
    $Row.titan_energy_kj_mol = $Row.energy_kj_mol
    $Row.status = [string]$energy.status

    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (Test-Path -LiteralPath $profilePath) {
        $profileRows = @(Import-Csv -LiteralPath $profilePath)
        if ($profileRows.Count -gt 0) {
            $profile = $profileRows[0]
            if ([string]$profile.status -ne 'ok' -and $Row.status -eq 'ok') {
                $Row.status = [string]$profile.status
            }
            $Row.panels = Get-Z4zProfileValue $profile 'panels' ''
            $Row.iterations = Get-Z4zProfileValue $profile 'iterations' ''
            $Row.sampled_audit_rms = Get-Z4zProfileValue $profile 'tree_vs_direct_matvec_rel' ''
            $Row.sampled_direct_true_residual = Get-Z4zProfileValue $profile 'fmm_local_sampled_direct_true_residual' (Get-Z4zProfileValue $profile 'sampled_direct_true_residual' '')
        }
    }
    if ($Row.status -eq 'ok') {
        $Row.failure_class = 'pass'
    } else {
        $Row.failure_class = 'solver_status_gap'
    }
    return $Row
}

function Invoke-Z4zTitanCase([object]$Case, [System.Collections.Generic.List[object]]$Rows) {
    $runDir = Join-Path $WorkRoot $Case.case_id
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $row = New-Z4zRow -Case $Case -RunDir $runDir
    $pqr = Copy-Z4zPqr -Case $Case -RunDir $runDir
    $arg = Write-FmmCliArgFile `
        -RunDir $runDir `
        -RepoRoot $RuntimeRepoRoot `
        -PqrPath $pqr `
        -Backend 'tabi' `
        -Mode 'treecode' `
        -FmmMode 'single-matvec-analytic' `
        -SurfacePoints $Case.surface_points `
        -MaxIter 720 `
        -Restart 128 `
        -Tolerance 1.0e-6 `
        -Threads $Threads `
        -SaltM $SaltM `
        -InternalDielectric $InternalDielectric `
        -ExternalDielectric $ExternalDielectric `
        -TemperatureK $TemperatureK `
        -OpeningTheta 0.6 `
        -AuditSampleCount 512 `
        -AuditMode 'sampled-detail' `
        -Preconditioner 'fmm-nearfield-ilu0' `
        -PreconditionerIluNeighborCap 64 `
        -NearIntegration 'analytic_singular' `
        -TaylorOrder 4 `
        -FmmLocal 'gmres-p4' `
        -FmmLocalRefine 'operator-tail-static' `
        -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache 'green-jet-compressed' `
        -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback 'packed-avx2' `
        -FmmLocalPlanBuilder 'metric-cache' `
        -FmmLocalGeometryGuard 'strict-geometry-mac' `
        -FmmLocalGeometryTheta 0.5
    $result = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $arg -CaseTimeoutSeconds $Case.timeout_seconds
    $row.exit_code = $result.exit_code
    $row.timed_out = [string]$result.timed_out
    $row.wall_seconds = $result.wall_seconds
    $row.peak_rss_mb = $result.peak_working_set_mb
    if ($result.timed_out) {
        $row.status = 'timeout'
        $row.failure_class = 'timeout'
    } elseif ($result.exit_code -ne 0) {
        $row.status = 'exit_nonzero'
        $row.failure_class = 'exit_nonzero'
    } else {
        $row = Read-Z4zTitanOutputs -RunDir $runDir -Row $row
        $energyValue = Get-Z4zDouble $row.energy_kj_mol
        $auditValue = Get-Z4zDouble $row.sampled_audit_rms
        if ($row.status -eq 'ok' -and (Test-FmmCliFinite $energyValue) -and ((-not (Test-FmmCliFinite $auditValue)) -or $auditValue -le 8.0e-3)) {
            $row.status = 'ok'
            $row.failure_class = 'pass'
        } elseif ($row.failure_class -eq 'pass') {
            $row.status = 'accuracy_gap'
            $row.failure_class = 'titan_accuracy_gap'
        }
    }
    Upsert-Z4zRow -Rows $Rows -Row $row
    Save-Z4zRows -Rows $Rows.ToArray()
    return $row
}

function Invoke-Z4zApbsCase([object]$Case, [System.Collections.Generic.List[object]]$Rows, [string]$ResolvedApbs) {
    $runDir = Join-Path $WorkRoot $Case.case_id
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $row = New-Z4zRow -Case $Case -RunDir $runDir
    if ([string]::IsNullOrWhiteSpace($ResolvedApbs)) {
        $row.status = 'skipped'
        $row.failure_class = 'skipped_apbs_missing'
        Upsert-Z4zRow -Rows $Rows -Row $row
        Save-Z4zRows -Rows $Rows.ToArray()
        return $row
    }
    $pqr = Copy-Z4zPqr -Case $Case -RunDir $runDir
    $input = Join-Path $runDir 'apbs.in'
    $box = Write-Z4zApbsInput -PqrPath $pqr -InputPath $input -Dime $Case.dime -FineScale $Case.fine_scale -CoarsePadding 32.0
    $row.cglen_x = $box.rx
    $row.cglen_y = $box.ry
    $row.cglen_z = $box.rz
    $row.fglen_x = [Math]::Max(12.0, $box.rx * [double]$Case.fine_scale)
    $row.fglen_y = [Math]::Max(12.0, $box.ry * [double]$Case.fine_scale)
    $row.fglen_z = [Math]::Max(12.0, $box.rz * [double]$Case.fine_scale)
    $row.atom_count = $box.atom_count
    $result = Invoke-Z4zExternalCase -Exe $ResolvedApbs -RunDir $runDir -Arguments 'apbs.in' -TimeoutSeconds $Case.timeout_seconds
    $row.apbs_exe = $ResolvedApbs
    $row.exit_code = $result.exit_code
    $row.timed_out = [string]$result.timed_out
    $row.wall_seconds = $result.wall_seconds
    $row.peak_rss_mb = $result.peak_working_set_mb
    if ($result.timed_out) {
        $row.status = 'timeout'
        $row.failure_class = 'timeout'
    } elseif ($result.exit_code -ne 0) {
        $row.status = 'exit_nonzero'
        $row.failure_class = 'exit_nonzero'
    } else {
        $energy = Read-Z4zApbsEnergy -RunDir $runDir
        $row.energy_kj_mol = $energy
        $row.apbs_energy_kj_mol = $energy
        if (Test-FmmCliFinite $energy) {
            $row.status = 'ok'
            $row.failure_class = 'pass'
        } else {
            $row.status = 'parse_gap'
            $row.failure_class = 'apbs_energy_parse_gap'
        }
    }
    Upsert-Z4zRow -Rows $Rows -Row $row
    Save-Z4zRows -Rows $Rows.ToArray()
    return $row
}

function Build-Z4zCases {
    $fixtureRows = @(Get-Z4zFixtureRows)
    $cases = New-Object System.Collections.Generic.List[object]
    foreach ($name in $Fixtures) {
        $fixture = Get-Z4zFixture -Name $name -Rows $fixtureRows
        if (-not $SkipTitan) {
            $cases.Add((New-Z4zTitanCase -Fixture $fixture)) | Out-Null
        }
        foreach ($scale in $FineScales) {
            foreach ($dime in $Dimes) {
                $cases.Add((New-Z4zApbsCase -Fixture $fixture -Dime $dime -FineScale $scale)) | Out-Null
            }
        }
    }
    return $cases.ToArray()
}

function Write-Z4zDerivedTables([object[]]$Rows) {
    $stability = @()
    foreach ($group in ($Rows | Where-Object { $_.backend -eq 'apbs' } | Group-Object fixture,fine_scale)) {
        $items = @($group.Group | Where-Object { $_.status -eq 'ok' } | Sort-Object { [int]$_.dime })
        if ($items.Count -eq 0) { continue }
        $finest = $items[-1]
        $finestEnergy = Get-Z4zDouble $finest.apbs_energy_kj_mol
        foreach ($item in $items) {
            $energy = Get-Z4zDouble $item.apbs_energy_kj_mol
            $delta = if ((Test-FmmCliFinite $energy) -and (Test-FmmCliFinite $finestEnergy)) { $energy - $finestEnergy } else { [double]::NaN }
            $rel = if ((Test-FmmCliFinite $delta) -and [Math]::Abs($finestEnergy) -gt 1e-12) { [Math]::Abs($delta) / [Math]::Abs($finestEnergy) } else { [double]::NaN }
            $titan = @($Rows | Where-Object { $_.fixture -eq $item.fixture -and $_.backend -eq 'titan_guarded_p4' -and $_.status -eq 'ok' }) | Select-Object -Last 1
            $titanEnergy = if ($null -ne $titan) { Get-Z4zDouble $titan.titan_energy_kj_mol } else { [double]::NaN }
            $titanDelta = if ((Test-FmmCliFinite $titanEnergy) -and (Test-FmmCliFinite $energy)) { $titanEnergy - $energy } else { [double]::NaN }
            $stability += [pscustomobject]@{
                fixture = $item.fixture
                protein_name = $item.protein_name
                pdb_id = $item.pdb_id
                fine_scale = $item.fine_scale
                dime = $item.dime
                apbs_energy_kj_mol = $energy
                finest_dime = $finest.dime
                finest_apbs_energy_kj_mol = $finestEnergy
                delta_to_finest_kj_mol = $delta
                relative_delta_to_finest = $rel
                titan_guarded_p4_energy_kj_mol = $titanEnergy
                titan_minus_apbs_kj_mol = $titanDelta
                apbs_wall_seconds = $item.wall_seconds
                apbs_peak_rss_mb = $item.peak_rss_mb
            }
        }
    }
    $stability | Export-Csv -LiteralPath $StabilityCsv -NoTypeInformation
}

function Write-Z4zSvg([object[]]$Rows) {
    $data = @($Rows | Where-Object { $_.backend -eq 'apbs' -and $_.status -eq 'ok' })
    if ($data.Count -eq 0) {
        Set-Content -LiteralPath $PlotPath -Value '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="520"><text x="40" y="80">No APBS rows available.</text></svg>' -Encoding ASCII
        return
    }
    $width = 920
    $height = 520
    $left = 80
    $right = 40
    $top = 40
    $bottom = 80
    $plotW = $width - $left - $right
    $plotH = $height - $top - $bottom
    $energies = @($data | ForEach-Object { Get-Z4zDouble $_.apbs_energy_kj_mol } | Where-Object { Test-FmmCliFinite $_ })
    $minE = ($energies | Measure-Object -Minimum).Minimum
    $maxE = ($energies | Measure-Object -Maximum).Maximum
    if ([Math]::Abs($maxE - $minE) -lt 1e-12) { $maxE = $minE + 1.0 }
    $minD = ($Dimes | Measure-Object -Minimum).Minimum
    $maxD = ($Dimes | Measure-Object -Maximum).Maximum
    if ($maxD -eq $minD) { $maxD = $minD + 1 }
    function X([double]$d) { return $left + (($d - $minD) / ($maxD - $minD)) * $plotW }
    function Y([double]$e) { return $top + (($maxE - $e) / ($maxE - $minE)) * $plotH }
    $palette = @('#1f77b4','#d62728','#2ca02c','#9467bd','#ff7f0e','#17becf')
    $lines = New-Object System.Collections.Generic.List[string]
    $lines.Add("<svg xmlns=`"http://www.w3.org/2000/svg`" width=`"$width`" height=`"$height`" viewBox=`"0 0 $width $height`">") | Out-Null
    $lines.Add('<rect width="100%" height="100%" fill="white"/>') | Out-Null
    $lines.Add("<line x1=`"$left`" y1=`"$($top + $plotH)`" x2=`"$($left + $plotW)`" y2=`"$($top + $plotH)`" stroke=`"black`"/>") | Out-Null
    $lines.Add("<line x1=`"$left`" y1=`"$top`" x2=`"$left`" y2=`"$($top + $plotH)`" stroke=`"black`"/>") | Out-Null
    foreach ($d in $Dimes) {
        $x = X $d
        $lines.Add("<line x1=`"$x`" y1=`"$($top + $plotH)`" x2=`"$x`" y2=`"$($top + $plotH + 6)`" stroke=`"black`"/>") | Out-Null
        $lines.Add("<text x=`"$x`" y=`"$($top + $plotH + 28)`" font-size=`"18`" text-anchor=`"middle`">$d</text>") | Out-Null
    }
    for ($i = 0; $i -le 4; ++$i) {
        $e = $minE + ($maxE - $minE) * $i / 4.0
        $y = Y $e
        $label = $e.ToString('G5', [Globalization.CultureInfo]::InvariantCulture)
        $lines.Add("<line x1=`"$($left - 6)`" y1=`"$y`" x2=`"$left`" y2=`"$y`" stroke=`"black`"/>") | Out-Null
        $lines.Add("<text x=`"$($left - 10)`" y=`"$($y + 6)`" font-size=`"16`" text-anchor=`"end`">$label</text>") | Out-Null
    }
    $lines.Add("<text x=`"$($left + $plotW / 2)`" y=`"$($height - 25)`" font-size=`"20`" text-anchor=`"middle`">APBS grid dimension</text>") | Out-Null
    $lines.Add("<text transform=`"translate(24 $($top + $plotH / 2)) rotate(-90)`" font-size=`"20`" text-anchor=`"middle`">Reaction energy</text>") | Out-Null
    $groups = @($data | Group-Object fixture | Sort-Object Name)
    for ($g = 0; $g -lt $groups.Count; ++$g) {
        $color = $palette[$g % $palette.Count]
        $points = @()
        foreach ($item in @($groups[$g].Group | Sort-Object { [int]$_.dime })) {
            $x = X ([double]$item.dime)
            $y = Y (Get-Z4zDouble $item.apbs_energy_kj_mol)
            $points += "$x,$y"
            $lines.Add("<circle cx=`"$x`" cy=`"$y`" r=`"5`" fill=`"$color`"/>") | Out-Null
        }
        if ($points.Count -gt 1) {
            $lines.Add("<polyline fill=`"none`" stroke=`"$color`" stroke-width=`"3`" points=`"$($points -join ' ')`"/>") | Out-Null
        }
        $legendY = $top + 24 + 24 * $g
        $lines.Add("<rect x=`"$($left + $plotW - 160)`" y=`"$($legendY - 12)`" width=`"16`" height=`"16`" fill=`"$color`"/>") | Out-Null
        $lines.Add("<text x=`"$($left + $plotW - 136)`" y=`"$($legendY + 2)`" font-size=`"16`">$($groups[$g].Name.ToUpperInvariant())</text>") | Out-Null
    }
    $lines.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $PlotPath -Value $lines -Encoding ASCII
}

function Write-Z4zDecision([object[]]$Rows) {
    Write-Z4zDerivedTables -Rows $Rows
    Write-Z4zSvg -Rows $Rows
    $apbsRows = @($Rows | Where-Object { $_.backend -eq 'apbs' })
    $titanRows = @($Rows | Where-Object { $_.backend -eq 'titan_guarded_p4' })
    $timeouts = @($Rows | Where-Object { $_.failure_class -eq 'timeout' })
    $apbsSkipped = ($apbsRows.Count -gt 0 -and (@($apbsRows | Where-Object { $_.failure_class -eq 'skipped_apbs_missing' }).Count -eq $apbsRows.Count))
    $nonPass = @($Rows | Where-Object { $_.failure_class -notin @('pass','skipped_apbs_missing') -and $_.status -notin @('ok','skipped') })
    $stabilityRows = @(Import-Csv -LiteralPath $StabilityCsv -ErrorAction SilentlyContinue)
    $unstable = @($stabilityRows | Where-Object {
        $rel = Get-Z4zDouble $_.relative_delta_to_finest
        (Test-FmmCliFinite $rel) -and $rel -gt $ApbsStabilityRelativeThreshold
    })
    $decision = if ($apbsSkipped) {
        'z4z_apbs_missing'
    } elseif ($timeouts.Count -gt 0) {
        'z4z_runtime_limited'
    } elseif ($nonPass.Count -gt 0) {
        'z4z_apbs_sweep_gap'
    } elseif ($unstable.Count -gt 0) {
        'z4z_apbs_grid_stability_gap'
    } else {
        'z4z_apbs_sweep_ready'
    }
    [pscustomobject]@{
        decision = $decision
        rows = $Rows.Count
        apbs_rows = $apbsRows.Count
        titan_rows = $titanRows.Count
        timeouts = $timeouts.Count
        non_pass = $nonPass.Count
        unstable_rows = $unstable.Count
        stability_relative_threshold = $ApbsStabilityRelativeThreshold
    } | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation
    @(
        "decision=$decision",
        "rows=$($Rows.Count)",
        "apbs_rows=$($apbsRows.Count)",
        "titan_rows=$($titanRows.Count)",
        "timeouts=$($timeouts.Count)",
        "matrix_csv=$MatrixCsv",
        "stability_csv=$StabilityCsv",
        "plot_svg=$PlotPath"
    ) | Set-Content -LiteralPath $SummaryPath -Encoding ascii
    @(
        "TITAN FAGBI Z4Z APBS grid sweep",
        "================================",
        "",
        "Decision: $decision",
        "Rows: $($Rows.Count)",
        "APBS rows: $($apbsRows.Count)",
        "TITAN rows: $($titanRows.Count)",
        "Timeouts: $($timeouts.Count)",
        "",
        "Matrix: $MatrixCsv",
        "Stability table: $StabilityCsv",
        "Plot: $PlotPath",
        "",
        "Notes:",
        "- APBS is used as an independent grid-solver sanity/context check, not as the dense boundary-integral oracle.",
        "- The same PQR input is used for APBS and TITAN, but their discretizations and energies are not expected to match exactly.",
        "- Rerun the same command to resume missing rows; add -RerunFailed for failed/timeout rows or -Force for all rows."
    ) | Set-Content -LiteralPath $RootReport -Encoding ascii
    return $decision
}

$rows = [System.Collections.Generic.List[object]]::new()
if (Test-Path -LiteralPath $MatrixCsv) {
    foreach ($row in Import-Csv -LiteralPath $MatrixCsv) { $rows.Add($row) | Out-Null }
}

$cases = @(Build-Z4zCases)

if ($PlanOnly) {
    $cases | Select-Object case_id,fixture,backend,dime,fine_scale,surface_points,timeout_seconds | Format-Table -AutoSize
    exit 0
}

$resolvedApbs = Resolve-Z4zApbsExe -Requested $ApbsExe
foreach ($case in $cases) {
    if (-not (Should-RunZ4zCase -Case $case -Rows $rows.ToArray())) {
        Write-Host "resume: $($case.case_id)"
        continue
    }
    Write-Host "run: $($case.case_id) (timeout $($case.timeout_seconds)s)"
    if ($case.backend -eq 'apbs') {
        [void](Invoke-Z4zApbsCase -Case $case -Rows $rows -ResolvedApbs $resolvedApbs)
    } else {
        [void](Invoke-Z4zTitanCase -Case $case -Rows $rows)
    }
}

$decision = Write-Z4zDecision -Rows $rows.ToArray()
Write-Host "Z4Z APBS decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Stability: $StabilityCsv"
Write-Host "Plot: $PlotPath"
Write-Host "Resume: rerun the same command to continue missing rows; use -RerunFailed for failed/timeout rows or -Force for all rows."
if (-not $AllowFailedDecision -and $decision -ne 'z4z_apbs_sweep_ready') {
    exit 1
}
