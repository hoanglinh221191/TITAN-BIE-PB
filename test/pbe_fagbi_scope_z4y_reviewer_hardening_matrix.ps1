param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$PrimaryFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$ExpansionFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public_z4q'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4y'),
    [ValidateSet('dense_oracle','mesh_sensitivity','apbs_sanity')]
    [string[]]$Stages = @('dense_oracle','mesh_sensitivity','apbs_sanity'),
    [string[]]$Fixtures = @(),
    [int[]]$SurfacePoints = @(),
    [string]$ApbsExe = '',
    [int]$DenseTimeoutSeconds = 1800,
    [int]$FmmTimeoutSeconds = 1800,
    [int]$ApbsTimeoutSeconds = 900,
    [int]$Threads = 4,
    [int]$MiniAtomCap = 64,
    [switch]$DenseOracleProteinOnly,
    [switch]$Quick,
    [switch]$PlanOnly,
    [switch]$Force,
    [switch]$RerunFailed,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$PrimaryFixtureRoot = [IO.Path]::GetFullPath($PrimaryFixtureRoot)
$ExpansionFixtureRoot = [IO.Path]::GetFullPath($ExpansionFixtureRoot)
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4y_matrix.csv'
$DenseOracleCsv = Join-Path $WorkRoot 'fagbi_scope_z4y_dense_oracle.csv'
$MeshCsv = Join-Path $WorkRoot 'fagbi_scope_z4y_mesh_sensitivity.csv'
$ApbsCsv = Join-Path $WorkRoot 'fagbi_scope_z4y_apbs_sanity.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4y_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4y_summary.dat'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4Y_REVIEWER_HARDENING_REPORT.txt'

$TemperatureK = 298.15
$InternalDielectric = 1.0
$ExternalDielectric = 80.0
$SaltM = 0.15

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

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

function Get-Z4yDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or $text -eq 'nan') {
        return $Default
    }
    $parsed = 0.0
    if ([double]::TryParse($text, [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-Z4yInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $parsed = [int64]0
    if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Format-Z4yDouble([double]$Value, [string]$Format = 'G6') {
    if (-not (Test-FmmCliFinite $Value)) { return 'NA' }
    return $Value.ToString($Format, [Globalization.CultureInfo]::InvariantCulture)
}

function Get-Z4yProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Get-Z4yObjectProperty([object]$Object, [string]$Name, [object]$Default = '') {
    if ($null -eq $Object) { return $Default }
    $prop = $Object.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Import-Z4yFixtureRows([string]$Root, [string]$ManifestName) {
    $path = Join-Path $Root $ManifestName
    if (-not (Test-Path -LiteralPath $path)) { return @() }
    return @(Import-Csv -LiteralPath $path | ForEach-Object {
        $_ | Add-Member -NotePropertyName fixture_root -NotePropertyValue $Root -Force
        $_
    })
}

function Get-Z4yFixtureRows {
    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($row in (Import-Z4yFixtureRows $PrimaryFixtureRoot 'public_fixture_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    foreach ($row in (Import-Z4yFixtureRows $ExpansionFixtureRoot 'public_fixture_expansion_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    return $rows.ToArray()
}

function Get-Z4yFixture([string]$Name, [object[]]$Rows) {
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

function New-Z4yMiniFixture([object]$BaseFixture, [string]$Name, [int]$AtomCap) {
    return [pscustomobject]@{
        fixture = $Name
        protein_name = "$($BaseFixture.protein_name) mini calibration subset"
        pdb_id = $BaseFixture.pdb_id
        fixture_class = 'public_protein_mini'
        atom_count = $AtomCap
        pqr_path = $BaseFixture.pqr_path
        pqr_sample_cap = $AtomCap
    }
}

function Get-Z4yDenseCalibrationFixture([string]$Name, [object[]]$Rows) {
    if ($Name.EndsWith('-mini', [StringComparison]::OrdinalIgnoreCase)) {
        $baseName = $Name.Substring(0, $Name.Length - 5).ToLowerInvariant()
        return New-Z4yMiniFixture `
            -BaseFixture (Get-Z4yFixture -Name $baseName -Rows $Rows) `
            -Name $Name.ToLowerInvariant() `
            -AtomCap $MiniAtomCap
    }
    return Get-Z4yFixture -Name $Name -Rows $Rows
}

function New-Z4yAnalyticFixture([string]$Name, [int]$Surface) {
    if ($Name -eq 'sphere_coarse') {
        return [pscustomobject]@{
            fixture = 'sphere_coarse'
            protein_name = 'Centered Born sphere'
            pdb_id = 'analytic'
            fixture_class = 'analytic_born'
            atom_count = 1
            surface_points = $Surface
            oracle_energy_kj_mol = Get-FagbiBornEnergyKJ `
                -Charge 1.0 `
                -RadiusA 3.0 `
                -SaltM $SaltM `
                -TemperatureK $TemperatureK `
                -InternalDielectric $InternalDielectric `
                -ExternalDielectric $ExternalDielectric
        }
    }
    return [pscustomobject]@{
        fixture = 'kirkwood_coarse'
        protein_name = 'Off-center Kirkwood sphere'
        pdb_id = 'analytic'
        fixture_class = 'analytic_kirkwood'
        atom_count = 2
        surface_points = $Surface
        oracle_energy_kj_mol = Get-FagbiKirkwoodEnergyKJ `
            -SphereRadiusA 10.0 `
            -Charge 1.0 `
            -OffsetA 3.0 `
            -SaltM $SaltM `
            -TemperatureK $TemperatureK `
            -InternalDielectric $InternalDielectric `
            -ExternalDielectric $ExternalDielectric `
            -LMax 160
    }
}

function Write-Z4yAnalyticPqr([object]$Fixture, [string]$Path) {
    if ($Fixture.fixture -eq 'sphere_coarse') {
        Write-FagbiBornPqr -Path $Path -Charge 1.0 -RadiusA 3.0
    } elseif ($Fixture.fixture -eq 'kirkwood_coarse') {
        Write-FagbiKirkwoodPqr -Path $Path -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0
    } else {
        throw "Not an analytic fixture: $($Fixture.fixture)"
    }
}

function Copy-Z4yPqrForCase([object]$Fixture, [string]$RunDir) {
    $pqr = Join-Path $RunDir 'input.pqr'
    if ($Fixture.fixture_class -like 'analytic_*') {
        Write-Z4yAnalyticPqr -Fixture $Fixture -Path $pqr
    } else {
        $sampleCap = [int](Get-Z4yObjectProperty -Object $Fixture -Name 'pqr_sample_cap' -Default 0)
        if ($sampleCap -gt 0) {
            Write-FmmCliSampledPqr -InputPath $Fixture.pqr_path -OutputPath $pqr -Cap $sampleCap
        } else {
            Copy-Item -LiteralPath $Fixture.pqr_path -Destination $pqr -Force
        }
    }
    return $pqr
}

function Resolve-Z4yApbsExe([string]$Requested) {
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
    $knownCondaApbs = @(
        (Join-Path $env:USERPROFILE '.conda\envs\titan_apbs\Library\bin\apbs.exe'),
        'C:\ProgramData\miniconda3\envs\titan_apbs\Library\bin\apbs.exe'
    )
    foreach ($candidate in $knownCondaApbs) {
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }
    return ''
}

function Get-Z4yPqrAtoms([string]$PqrPath) {
    $atoms = New-Object System.Collections.Generic.List[object]
    foreach ($line in Get-Content -LiteralPath $PqrPath) {
        if ($line -notmatch '^(ATOM|HETATM)') { continue }
        $parts = @($line -split '\s+' | Where-Object { $_ -ne '' })
        if ($parts.Count -lt 10) { continue }
        $x = Get-Z4yDouble $parts[5]
        $y = Get-Z4yDouble $parts[6]
        $z = Get-Z4yDouble $parts[7]
        $charge = Get-Z4yDouble $parts[8] 0.0
        $radius = Get-Z4yDouble $parts[9] 1.5
        $atoms.Add([pscustomobject]@{ x=$x; y=$y; z=$z; charge=$charge; radius=$radius }) | Out-Null
    }
    return $atoms.ToArray()
}

function Write-Z4yApbsInput([string]$PqrPath, [string]$InputPath) {
    $atoms = @(Get-Z4yPqrAtoms -PqrPath $PqrPath)
    if ($atoms.Count -eq 0) { throw "Cannot derive APBS box from PQR: $PqrPath" }
    $xmin = ($atoms | Measure-Object x -Minimum).Minimum
    $xmax = ($atoms | Measure-Object x -Maximum).Maximum
    $ymin = ($atoms | Measure-Object y -Minimum).Minimum
    $ymax = ($atoms | Measure-Object y -Maximum).Maximum
    $zmin = ($atoms | Measure-Object z -Minimum).Minimum
    $zmax = ($atoms | Measure-Object z -Maximum).Maximum
    $rx = [Math]::Max(20.0, ($xmax - $xmin) + 32.0)
    $ry = [Math]::Max(20.0, ($ymax - $ymin) + 32.0)
    $rz = [Math]::Max(20.0, ($zmax - $zmin) + 32.0)
    $fmt = { param([double]$v) $v.ToString('G12', [Globalization.CultureInfo]::InvariantCulture) }
    $pqrName = Split-Path -Leaf $PqrPath
    $coarseLen = "    cglen $(& $fmt $rx) $(& $fmt $ry) $(& $fmt $rz)"
    $fineLen = "    fglen $(& $fmt ([Math]::Max(12.0, $rx * 0.55))) $(& $fmt ([Math]::Max(12.0, $ry * 0.55))) $(& $fmt ([Math]::Max(12.0, $rz * 0.55)))"
    $commonSolvated = @(
        '    mg-auto',
        '    dime 129 129 129',
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
        '    dime 129 129 129',
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
}

function Invoke-Z4yExternalCase([string]$Exe, [string]$RunDir, [string]$Arguments, [int]$TimeoutSeconds) {
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

function Read-Z4yApbsEnergy([string]$RunDir) {
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
        if ($m.Success) { return [double]::Parse($m.Groups[1].Value, [Globalization.CultureInfo]::InvariantCulture) }
    }
    return [double]::NaN
}

function New-Z4yRow([object]$Case, [string]$RunDir) {
    return [pscustomobject]@{
        case_id = $Case.case_id
        stage = $Case.stage
        fixture = $Case.fixture
        protein_name = $Case.protein_name
        pdb_id = $Case.pdb_id
        fixture_class = $Case.fixture_class
        backend = $Case.backend
        surface_points = $Case.surface_points
        timeout_seconds = $Case.timeout_seconds
        status = 'planned'
        failure_class = 'planned'
        exit_code = ''
        timed_out = ''
        panels = ''
        iterations = ''
        energy_kj_mol = ''
        oracle_energy_kj_mol = $Case.oracle_energy_kj_mol
        direct_dense_energy_kj_mol = ''
        energy_delta_kj_mol = ''
        energy_relative_delta = ''
        sampled_audit_rms = ''
        sampled_direct_true_residual = ''
        true_residual = ''
        final_residual = ''
        matvec_seconds_avg = ''
        wall_seconds = ''
        peak_rss_mb = ''
        fmm_near_pairs = ''
        fmm_far_pairs = ''
        accepted_far_pairs = ''
        direct_fallback_pairs = ''
        total_evaluated_interactions = ''
        plan_build_count = ''
        cache_build_count = ''
        apbs_energy_kj_mol = ''
        apbs_energy_raw = ''
        apbs_delta_kj_mol = ''
        apbs_exe = ''
        run_dir = $RunDir
    }
}

function Read-Z4yTitanOutputs([string]$RunDir, [object]$Row) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    if (-not (Test-Path -LiteralPath $energyPath)) {
        $Row.failure_class = 'missing_energy_output'
        return $Row
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $Row.energy_kj_mol = Get-Z4yDouble $energy.E_rxn_kJmol
    $Row.status = [string]$energy.status

    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (Test-Path -LiteralPath $profilePath) {
        $profile = @(Import-Csv -LiteralPath $profilePath)[0]
        if ([string]$profile.status -ne 'ok' -and $Row.status -eq 'ok') {
            $Row.status = [string]$profile.status
        }
        $Row.panels = Get-Z4yInt64 $profile.panels
        $Row.iterations = Get-Z4yInt64 $profile.iterations
        $Row.final_residual = Get-Z4yDouble $profile.final_residual
        $Row.true_residual = Get-Z4yDouble $profile.true_residual_rel
        $Row.sampled_audit_rms = Get-Z4yDouble $profile.tree_vs_direct_matvec_rel
        $Row.sampled_direct_true_residual = Get-Z4yDouble (Get-Z4yProfileValue $profile 'fmm_local_sampled_direct_true_residual' '')
        $Row.matvec_seconds_avg = Get-Z4yDouble $profile.matvec_seconds_avg
        $Row.fmm_near_pairs = Get-Z4yInt64 (Get-Z4yProfileValue $profile 'fmm_m2l_near_pairs' '')
        $Row.fmm_far_pairs = Get-Z4yInt64 (Get-Z4yProfileValue $profile 'fmm_m2l_far_pairs' '')
        $Row.accepted_far_pairs = Get-Z4yInt64 (Get-Z4yProfileValue $profile 'fmm_local_far_pairs_after' '')
        $Row.direct_fallback_pairs = Get-Z4yInt64 (Get-Z4yProfileValue $profile 'fmm_local_direct_fallback_panels' '')
        $Row.total_evaluated_interactions = Get-Z4yInt64 (Get-Z4yProfileValue $profile 'fmm_local_total_evaluated_interactions' '')
        $Row.plan_build_count = Get-Z4yInt64 (Get-Z4yProfileValue $profile 'fmm_local_plan_build_count' '')
        $Row.cache_build_count = Get-Z4yInt64 (Get-Z4yProfileValue $profile 'fmm_local_cache_build_count' '')
    }
    if ($Row.status -eq 'ok') {
        $Row.failure_class = 'pass'
    } else {
        $Row.failure_class = 'solver_status_gap'
    }
    return $Row
}

function Should-RunZ4yCase([object]$Case, [object[]]$Rows) {
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

function Upsert-Z4yRow([System.Collections.Generic.List[object]]$Rows, [object]$Row) {
    for ($i = 0; $i -lt $Rows.Count; ++$i) {
        if ($Rows[$i].case_id -eq $Row.case_id) {
            $Rows[$i] = $Row
            return
        }
    }
    $Rows.Add($Row) | Out-Null
}

function Save-Z4yRows([object[]]$Rows) {
    $Rows | Sort-Object stage, fixture, surface_points, backend | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
}

function New-Z4yTitanCase([string]$Stage, [object]$Fixture, [int]$Surface, [string]$Backend) {
    $id = '{0}_{1}_s{2}_{3}' -f $Stage,$Fixture.fixture,$Surface,$Backend
    $fixtureClass = Get-Z4yObjectProperty -Object $Fixture -Name 'fixture_class' -Default 'public_protein'
    return [pscustomobject]@{
        case_id = $id
        stage = $Stage
        fixture = $Fixture.fixture
        protein_name = $Fixture.protein_name
        pdb_id = $Fixture.pdb_id
        fixture_class = $fixtureClass
        surface_points = $Surface
        backend = $Backend
        timeout_seconds = if ($Backend -eq 'tabi_direct_dense') { $DenseTimeoutSeconds } else { $FmmTimeoutSeconds }
        oracle_energy_kj_mol = Get-Z4yObjectProperty -Object $Fixture -Name 'oracle_energy_kj_mol' -Default ([double]::NaN)
        pqr_source = Get-Z4yObjectProperty -Object $Fixture -Name 'pqr_path' -Default ''
        pqr_sample_cap = Get-Z4yObjectProperty -Object $Fixture -Name 'pqr_sample_cap' -Default 0
    }
}

function New-Z4yApbsCase([object]$Fixture) {
    return [pscustomobject]@{
        case_id = 'apbs_sanity_{0}_apbs' -f $Fixture.fixture
        stage = 'apbs_sanity'
        fixture = $Fixture.fixture
        protein_name = $Fixture.protein_name
        pdb_id = $Fixture.pdb_id
        fixture_class = 'public_protein'
        surface_points = 0
        backend = 'apbs'
        timeout_seconds = $ApbsTimeoutSeconds
        oracle_energy_kj_mol = [double]::NaN
        pqr_source = $Fixture.pqr_path
        pqr_sample_cap = Get-Z4yObjectProperty -Object $Fixture -Name 'pqr_sample_cap' -Default 0
    }
}

function Invoke-Z4yTitanCase([object]$Case, [System.Collections.Generic.List[object]]$Rows) {
    $runDir = Join-Path $WorkRoot $Case.case_id
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $row = New-Z4yRow -Case $Case -RunDir $runDir
    $pqr = Copy-Z4yPqrForCase -Fixture ([pscustomobject]@{
        fixture = $Case.fixture
        fixture_class = $Case.fixture_class
        pqr_path = $Case.pqr_source
        pqr_sample_cap = Get-Z4yObjectProperty -Object $Case -Name 'pqr_sample_cap' -Default 0
    }) -RunDir $runDir

    if ($Case.backend -eq 'tabi_direct_dense') {
        $isPublicProtein = ([string]$Case.fixture_class -like 'public_protein*')
        $directPreconditioner = if ($isPublicProtein) { 'block-budgeted' } else { 'diagonal' }
        $directMaxIter = if ($isPublicProtein) { 1440 } else { 720 }
        $arg = Write-FmmCliArgFile `
            -RunDir $runDir `
            -RepoRoot $RuntimeRepoRoot `
            -PqrPath $pqr `
            -Backend 'tabi' `
            -Mode 'direct' `
            -FmmMode 'off' `
            -SurfacePoints $Case.surface_points `
            -MaxIter $directMaxIter `
            -Restart 128 `
            -Tolerance 1.0e-6 `
            -Threads $Threads `
            -SaltM $SaltM `
            -InternalDielectric $InternalDielectric `
            -ExternalDielectric $ExternalDielectric `
            -TemperatureK $TemperatureK `
            -OpeningTheta 0.6 `
            -AuditSampleCount 0 `
            -AuditMode 'off' `
            -Preconditioner $directPreconditioner `
            -NearIntegration 'analytic_singular' `
            -TaylorOrder 4
        if ($isPublicProtein) {
            Add-Content -LiteralPath $arg -Value @(
                '-pbe-tabi-precond-block-near=buffer-aware',
                '-pbe-tabi-precond-memory-budget-mb=256'
            ) -Encoding ASCII
        }
    } else {
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
    }
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
        $row = Read-Z4yTitanOutputs -RunDir $runDir -Row $row
    }
    Upsert-Z4yRow -Rows $Rows -Row $row
    Save-Z4yRows -Rows $Rows.ToArray()
    return $row
}

function Invoke-Z4yApbsCase([object]$Case, [System.Collections.Generic.List[object]]$Rows, [string]$ResolvedApbs) {
    $runDir = Join-Path $WorkRoot $Case.case_id
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $row = New-Z4yRow -Case $Case -RunDir $runDir
    if ([string]::IsNullOrWhiteSpace($ResolvedApbs)) {
        $row.status = 'skipped'
        $row.failure_class = 'skipped_apbs_missing'
        Upsert-Z4yRow -Rows $Rows -Row $row
        Save-Z4yRows -Rows $Rows.ToArray()
        return $row
    }
    $pqr = Copy-Z4yPqrForCase -Fixture ([pscustomobject]@{
        fixture = $Case.fixture
        fixture_class = $Case.fixture_class
        pqr_path = $Case.pqr_source
        pqr_sample_cap = Get-Z4yObjectProperty -Object $Case -Name 'pqr_sample_cap' -Default 0
    }) -RunDir $runDir
    $input = Join-Path $runDir 'apbs.in'
    Write-Z4yApbsInput -PqrPath $pqr -InputPath $input
    $result = Invoke-Z4yExternalCase -Exe $ResolvedApbs -RunDir $runDir -Arguments 'apbs.in' -TimeoutSeconds $Case.timeout_seconds
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
        $raw = Read-Z4yApbsEnergy -RunDir $runDir
        $row.apbs_energy_raw = $raw
        $row.apbs_energy_kj_mol = $raw
        if (Test-FmmCliFinite $raw) {
            $row.status = 'ok'
            $row.failure_class = 'pass'
        } else {
            $row.status = 'parse_gap'
            $row.failure_class = 'apbs_energy_parse_gap'
        }
    }
    Upsert-Z4yRow -Rows $Rows -Row $row
    Save-Z4yRows -Rows $Rows.ToArray()
    return $row
}

function Build-Z4yCases {
    $fixtureRows = @(Get-Z4yFixtureRows)
    $cases = New-Object System.Collections.Generic.List[object]

    if ($Stages -contains 'dense_oracle') {
        $denseSurfaces = if ($Quick) { @(8) } else { @(8,12,24) }
        if ($SurfacePoints.Count -gt 0) { $denseSurfaces = $SurfacePoints }
        if (-not $DenseOracleProteinOnly) {
            foreach ($surface in $denseSurfaces) {
                foreach ($analytic in @('sphere_coarse','kirkwood_coarse')) {
                    $fx = New-Z4yAnalyticFixture -Name $analytic -Surface $surface
                    foreach ($backend in @('tabi_direct_dense','fmm_p4_guarded')) {
                        $cases.Add((New-Z4yTitanCase -Stage 'dense_oracle' -Fixture $fx -Surface $surface -Backend $backend)) | Out-Null
                    }
                }
            }
        }
        $smallFixtures = if ($Quick) { @('1vii-mini') } else { @('1vii-mini') }
        if ($Fixtures.Count -gt 0) { $smallFixtures = @($Fixtures | Where-Object { $_ -notin @('sphere_coarse','kirkwood_coarse') }) }
        $smallSurfaces = if ($SurfacePoints.Count -gt 0) { $SurfacePoints } else { @(8) }
        foreach ($name in $smallFixtures) {
            $fx = Get-Z4yDenseCalibrationFixture -Name $name -Rows $fixtureRows
            foreach ($surface in $smallSurfaces) {
                foreach ($backend in @('tabi_direct_dense','fmm_p4_guarded')) {
                    $cases.Add((New-Z4yTitanCase -Stage 'dense_oracle' -Fixture $fx -Surface $surface -Backend $backend)) | Out-Null
                }
            }
        }
    }

    if ($Stages -contains 'mesh_sensitivity') {
        $meshFixtures = if ($Quick) { @('1vii') } else { @('1vii','1ubq') }
        if ($Fixtures.Count -gt 0) { $meshFixtures = @($Fixtures | Where-Object { $_ -notin @('sphere_coarse','kirkwood_coarse') }) }
        $meshSurfaces = if ($SurfacePoints.Count -gt 0) { $SurfacePoints } elseif ($Quick) { @(8,24) } else { @(8,24,96) }
        foreach ($name in $meshFixtures) {
            $fx = Get-Z4yFixture -Name $name -Rows $fixtureRows
            foreach ($surface in $meshSurfaces) {
                $cases.Add((New-Z4yTitanCase -Stage 'mesh_sensitivity' -Fixture $fx -Surface $surface -Backend 'fmm_p4_guarded')) | Out-Null
            }
        }
    }

    if ($Stages -contains 'apbs_sanity') {
        $apbsFixtures = if ($Quick) { @('1vii') } else { @('1vii','1ubq') }
        if ($Fixtures.Count -gt 0) { $apbsFixtures = @($Fixtures | Where-Object { $_ -notin @('sphere_coarse','kirkwood_coarse') }) }
        foreach ($name in $apbsFixtures) {
            $fx = Get-Z4yFixture -Name $name -Rows $fixtureRows
            $cases.Add((New-Z4yTitanCase -Stage 'apbs_sanity' -Fixture $fx -Surface 8 -Backend 'fmm_p4_guarded')) | Out-Null
            $cases.Add((New-Z4yApbsCase -Fixture $fx)) | Out-Null
        }
    }

    return $cases.ToArray()
}

function Write-Z4yDerivedTables([object[]]$Rows) {
    $denseRows = @()
    foreach ($group in ($Rows | Where-Object { $_.stage -eq 'dense_oracle' } | Group-Object fixture,surface_points)) {
        $items = @($group.Group)
        $direct = @($items | Where-Object { $_.backend -eq 'tabi_direct_dense' }) | Select-Object -Last 1
        $p4 = @($items | Where-Object { $_.backend -eq 'fmm_p4_guarded' }) | Select-Object -Last 1
        if ($null -eq $direct -or $null -eq $p4) { continue }
        $directEnergy = Get-Z4yDouble $direct.energy_kj_mol
        $p4Energy = Get-Z4yDouble $p4.energy_kj_mol
        $delta = if ((Test-FmmCliFinite $directEnergy) -and (Test-FmmCliFinite $p4Energy)) { $p4Energy - $directEnergy } else { [double]::NaN }
        $rel = if ((Test-FmmCliFinite $delta) -and [Math]::Abs($directEnergy) -gt 1e-12) { [Math]::Abs($delta) / [Math]::Abs($directEnergy) } else { [double]::NaN }
        $denseRows += [pscustomobject]@{
            case = $direct.fixture
            panels = $p4.panels
            surface_points = $p4.surface_points
            direct_dense_status = $direct.status
            guarded_p4_status = $p4.status
            direct_dense_energy_kj_mol = $directEnergy
            guarded_p4_energy_kj_mol = $p4Energy
            energy_error_kj_mol = $delta
            energy_relative_error = $rel
            operator_residual = $p4.sampled_audit_rms
            sampled_direct_true_residual = $p4.sampled_direct_true_residual
            wall_direct_s = $direct.wall_seconds
            wall_p4_s = $p4.wall_seconds
        }
    }
    $denseRows | Export-Csv -LiteralPath $DenseOracleCsv -NoTypeInformation

    $Rows |
        Where-Object { $_.stage -eq 'mesh_sensitivity' } |
        Select-Object fixture,protein_name,pdb_id,surface_points,panels,energy_kj_mol,sampled_audit_rms,sampled_direct_true_residual,iterations,wall_seconds,peak_rss_mb,status,failure_class |
        Export-Csv -LiteralPath $MeshCsv -NoTypeInformation

    $apbsRows = @()
    foreach ($group in ($Rows | Where-Object { $_.stage -eq 'apbs_sanity' } | Group-Object fixture)) {
        $items = @($group.Group)
        $titan = @($items | Where-Object { $_.backend -eq 'fmm_p4_guarded' }) | Select-Object -Last 1
        $apbs = @($items | Where-Object { $_.backend -eq 'apbs' }) | Select-Object -Last 1
        if ($null -eq $titan -and $null -eq $apbs) { continue }
        $titanEnergy = if ($null -ne $titan) { Get-Z4yDouble $titan.energy_kj_mol } else { [double]::NaN }
        $apbsEnergy = if ($null -ne $apbs) { Get-Z4yDouble $apbs.apbs_energy_kj_mol } else { [double]::NaN }
        $delta = if ((Test-FmmCliFinite $titanEnergy) -and (Test-FmmCliFinite $apbsEnergy)) { $titanEnergy - $apbsEnergy } else { [double]::NaN }
        $apbsRows += [pscustomobject]@{
            protein = if ($null -ne $titan) { $titan.protein_name } else { $apbs.protein_name }
            pdb = if ($null -ne $titan) { $titan.pdb_id } else { $apbs.pdb_id }
            titan_status = if ($null -ne $titan) { $titan.status } else { 'missing' }
            apbs_status = if ($null -ne $apbs) { $apbs.status } else { 'missing' }
            titan_energy_kj_mol = $titanEnergy
            apbs_energy_reported = $apbsEnergy
            titan_minus_apbs = $delta
            notes = 'APBS grid solve uses the same PQR input but not the same boundary-integral surface discretization.'
        }
    }
    $apbsRows | Export-Csv -LiteralPath $ApbsCsv -NoTypeInformation
}

function Write-Z4yDecision([object[]]$Rows, [string[]]$RequiredStages) {
    $needsDense = $RequiredStages -contains 'dense_oracle'
    $needsMesh = $RequiredStages -contains 'mesh_sensitivity'
    $needsApbs = $RequiredStages -contains 'apbs_sanity'
    $required = @($Rows | Where-Object { $_.backend -ne 'apbs' -or $_.failure_class -ne 'skipped_apbs_missing' })
    $timeouts = @($required | Where-Object { $_.failure_class -eq 'timeout' })
    $nonPass = @($required | Where-Object {
        $_.failure_class -notin @('pass','skipped_apbs_missing') -and $_.status -notin @('ok','skipped')
    })
    $densePairs = @(Import-Csv -LiteralPath $DenseOracleCsv -ErrorAction SilentlyContinue)
    $densePass = -not $needsDense
    if ($needsDense -and $densePairs.Count -gt 0) {
        $densePass = $true
        foreach ($row in $densePairs) {
            $rel = Get-Z4yDouble $row.energy_relative_error
            $op = Get-Z4yDouble $row.operator_residual
            if ($row.direct_dense_status -ne 'ok' -or $row.guarded_p4_status -ne 'ok' -or
                -not (Test-FmmCliFinite $rel) -or $rel -gt 0.05 -or
                -not (Test-FmmCliFinite $op) -or $op -gt 8.0e-3) {
                $densePass = $false
            }
        }
    }
    $meshRows = @($Rows | Where-Object { $_.stage -eq 'mesh_sensitivity' })
    $meshPass = ((-not $needsMesh) -or ($meshRows.Count -gt 0 -and (@($meshRows | Where-Object { $_.status -ne 'ok' }).Count -eq 0)))
    $apbsRows = @($Rows | Where-Object { $_.stage -eq 'apbs_sanity' -and $_.backend -eq 'apbs' })
    $apbsPass = ((-not $needsApbs) -or ($apbsRows.Count -gt 0 -and (@($apbsRows | Where-Object { $_.status -eq 'ok' }).Count -gt 0)))
    $apbsSkipped = ($needsApbs -and $apbsRows.Count -gt 0 -and (@($apbsRows | Where-Object { $_.failure_class -eq 'skipped_apbs_missing' }).Count -eq $apbsRows.Count))

    $decision = if ($timeouts.Count -gt 0) {
        'z4y_runtime_limited'
    } elseif (-not $densePass) {
        'z4y_dense_oracle_gap'
    } elseif (-not $meshPass) {
        'z4y_mesh_sensitivity_gap'
    } elseif ($apbsSkipped) {
        'z4y_ready_without_apbs'
    } elseif (-not $apbsPass) {
        'z4y_apbs_sanity_gap'
    } else {
        'z4y_reviewer_hardening_pass'
    }
    [pscustomobject]@{
        decision = $decision
        rows = $Rows.Count
        timeouts = $timeouts.Count
        non_pass = $nonPass.Count
        dense_pairs = $densePairs.Count
        mesh_rows = $meshRows.Count
        apbs_rows = $apbsRows.Count
    } | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation
    @(
        "decision=$decision",
        "rows=$($Rows.Count)",
        "timeouts=$($timeouts.Count)",
        "dense_oracle_csv=$DenseOracleCsv",
        "mesh_sensitivity_csv=$MeshCsv",
        "apbs_sanity_csv=$ApbsCsv"
    ) | Set-Content -LiteralPath $SummaryPath -Encoding ascii
    @(
        "TITAN FAGBI Z4Y reviewer hardening matrix",
        "===========================================",
        "",
        "Decision: $decision",
        "Rows: $($Rows.Count)",
        "Timeouts: $($timeouts.Count)",
        "",
        "Dense small oracle: $DenseOracleCsv",
        "Mesh sensitivity: $MeshCsv",
        "APBS sanity: $ApbsCsv",
        "",
        "Notes:",
        "- Full-protein dense O(N^2) reaction-energy oracles are intentionally not run here.",
        "- APBS rows are skipped if APBS is not installed or -ApbsExe is not supplied.",
        "- Rerun the same command to continue missing rows; add -RerunFailed for failed/timeout rows or -Force for all rows."
    ) | Set-Content -LiteralPath $RootReport -Encoding ascii
    return $decision
}

$rows = [System.Collections.Generic.List[object]]::new()
if (Test-Path -LiteralPath $MatrixCsv) {
    foreach ($row in Import-Csv -LiteralPath $MatrixCsv) { $rows.Add($row) | Out-Null }
}

$cases = @(Build-Z4yCases)
$plannedCaseIds = New-Object 'System.Collections.Generic.HashSet[string]'
foreach ($case in $cases) { [void]$plannedCaseIds.Add([string]$case.case_id) }
$activeStages = New-Object 'System.Collections.Generic.HashSet[string]'
foreach ($stage in $Stages) { [void]$activeStages.Add([string]$stage) }
$filteredRows = [System.Collections.Generic.List[object]]::new()
foreach ($row in $rows.ToArray()) {
    if (-not $activeStages.Contains([string]$row.stage)) {
        continue
    }
    if (-not $plannedCaseIds.Contains([string]$row.case_id)) {
        continue
    }
    $filteredRows.Add($row) | Out-Null
}
$rows = $filteredRows

if ($PlanOnly) {
    $cases | Select-Object case_id,stage,fixture,surface_points,backend,timeout_seconds | Format-Table -AutoSize
    exit 0
}

$resolvedApbs = Resolve-Z4yApbsExe -Requested $ApbsExe
foreach ($case in $cases) {
    if (-not (Should-RunZ4yCase -Case $case -Rows $rows.ToArray())) {
        Write-Host "resume: $($case.case_id)"
        continue
    }
    Write-Host "run: $($case.case_id) (timeout $($case.timeout_seconds)s)"
    if ($case.backend -eq 'apbs') {
        [void](Invoke-Z4yApbsCase -Case $case -Rows $rows -ResolvedApbs $resolvedApbs)
    } else {
        [void](Invoke-Z4yTitanCase -Case $case -Rows $rows)
    }
}

Write-Z4yDerivedTables -Rows $rows.ToArray()
$decision = Write-Z4yDecision -Rows $rows.ToArray() -RequiredStages $Stages
Write-Host "Z4Y decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Dense oracle: $DenseOracleCsv"
Write-Host "Mesh sensitivity: $MeshCsv"
Write-Host "APBS sanity: $ApbsCsv"
Write-Host "Resume: rerun the same command to continue missing rows; use -RerunFailed for failed/timeout rows or -Force for all rows."
if (-not $AllowFailedDecision -and $decision -notin @('z4y_reviewer_hardening_pass','z4y_ready_without_apbs')) {
    exit 1
}
