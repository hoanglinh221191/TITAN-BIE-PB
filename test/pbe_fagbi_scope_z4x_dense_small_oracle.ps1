param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4x'),
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 4,
    [switch]$Quick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
if (-not [IO.Path]::IsPathRooted($WorkRoot)) {
    $WorkRoot = Join-Path $RepoRoot $WorkRoot
}
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
$TmpRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp'))
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4x_dense_oracle_matrix.csv'
$EnergyParityCsv = Join-Path $WorkRoot 'fagbi_scope_z4x_energy_parity.csv'
$OperatorAuditCsv = Join-Path $WorkRoot 'fagbi_scope_z4x_operator_audit.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4x_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4x_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4x_dense_oracle_plot.svg'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4X_DENSE_SMALL_ORACLE_REPORT.txt'

$TemperatureK = 298.15
$InternalDielectric = 1.0
$ExternalDielectric = 80.0
$SaltM = 0.15

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

if (-not $WorkRoot.StartsWith($TmpRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Z4X work root outside test/tmp: $WorkRoot"
}
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot $TmpRoot

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

function Get-Z4xDouble([object]$Value, [double]$Default = [double]::NaN) {
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

function Get-Z4xInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $parsed = [int64]0
    if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-Z4xProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Get-Z4xAnalyticThreshold([int]$SurfacePoints) {
    if ($SurfacePoints -ge 24) { return 1.0e-2 }
    if ($SurfacePoints -ge 12) { return 2.0e-2 }
    return 5.0e-2
}

function New-Z4xFixture([string]$Fixture, [int]$SurfacePoints) {
    if ($Fixture -eq 'born_centered') {
        return [pscustomobject]@{
            fixture = 'born_centered'
            fixture_class = 'analytic_born'
            charge = 1.0
            radius_A = 3.0
            offset_A = 0.0
            surface_points = $SurfacePoints
            oracle_energy_kj_mol = Get-FagbiBornEnergyKJ `
                -Charge 1.0 `
                -RadiusA 3.0 `
                -SaltM $SaltM `
                -TemperatureK $TemperatureK `
                -InternalDielectric $InternalDielectric `
                -ExternalDielectric $ExternalDielectric
            oracle_threshold = Get-Z4xAnalyticThreshold -SurfacePoints $SurfacePoints
        }
    }
    return [pscustomobject]@{
        fixture = 'kirkwood_offcenter'
        fixture_class = 'analytic_kirkwood'
        charge = 1.0
        radius_A = 10.0
        offset_A = 3.0
        surface_points = $SurfacePoints
        oracle_energy_kj_mol = Get-FagbiKirkwoodEnergyKJ `
            -SphereRadiusA 10.0 `
            -Charge 1.0 `
            -OffsetA 3.0 `
            -SaltM $SaltM `
            -TemperatureK $TemperatureK `
            -InternalDielectric $InternalDielectric `
            -ExternalDielectric $ExternalDielectric `
            -LMax 160
        oracle_threshold = Get-Z4xAnalyticThreshold -SurfacePoints $SurfacePoints
    }
}

function Write-Z4xPqr([object]$Fixture, [string]$Path) {
    if ($Fixture.fixture -eq 'born_centered') {
        Write-FagbiBornPqr -Path $Path -Charge $Fixture.charge -RadiusA $Fixture.radius_A
    } else {
        Write-FagbiKirkwoodPqr `
            -Path $Path `
            -SphereRadiusA $Fixture.radius_A `
            -Charge $Fixture.charge `
            -OffsetA $Fixture.offset_A
    }
}

function New-Z4xBackend([string]$Name) {
    switch ($Name) {
        'tabi_direct_dense' {
            return [pscustomobject]@{
                backend = 'tabi_direct_dense'
                mode = 'direct'
                fmm_local = 'off'
                fmm_refine = 'off'
                fmm_cache = 'off'
                fmm_fallback = 'scalar'
                fmm_plan_builder = 'legacy-rescan'
                geometry_guard = 'off'
                preconditioner = 'diagonal'
                require_audit = $false
            }
        }
        'fmm_degree1_ilu0' {
            return [pscustomobject]@{
                backend = 'fmm_degree1_ilu0'
                mode = 'treecode'
                fmm_local = 'off'
                fmm_refine = 'off'
                fmm_cache = 'off'
                fmm_fallback = 'scalar'
                fmm_plan_builder = 'legacy-rescan'
                geometry_guard = 'off'
                preconditioner = 'fmm-nearfield-ilu0'
                require_audit = $true
            }
        }
        default {
            return [pscustomobject]@{
                backend = 'fmm_p4_guarded'
                mode = 'treecode'
                fmm_local = 'gmres-p4'
                fmm_refine = 'operator-tail-static'
                fmm_cache = 'green-jet-compressed'
                fmm_fallback = 'packed-avx2'
                fmm_plan_builder = 'metric-cache'
                geometry_guard = 'strict-geometry-mac'
                preconditioner = 'fmm-nearfield-ilu0'
                require_audit = $true
            }
        }
    }
}

function New-Z4xEmptyRow([object]$Fixture, [object]$Backend, [string]$RunDir) {
    return [pscustomobject]@{
        fixture = $Fixture.fixture
        fixture_class = $Fixture.fixture_class
        surface_points = $Fixture.surface_points
        backend = $Backend.backend
        status = 'planned'
        failure_class = 'planned'
        exit_code = 0
        timed_out = 'False'
        panels = 0
        iterations = 0
        energy_kj_mol = [double]::NaN
        oracle_energy_kj_mol = $Fixture.oracle_energy_kj_mol
        direct_dense_energy_kj_mol = [double]::NaN
        oracle_relative_error = [double]::NaN
        p4_vs_direct_delta_kj_mol = [double]::NaN
        p4_vs_direct_relative_error = [double]::NaN
        direct_oracle_threshold = $Fixture.oracle_threshold
        p4_energy_gate_kj_mol = [double]::NaN
        final_residual = [double]::NaN
        true_residual = [double]::NaN
        sampled_audit_rms = [double]::NaN
        sampled_direct_true_residual = [double]::NaN
        matvec_audit_rows = 0
        matvec_seconds_avg = [double]::NaN
        wall_seconds = [double]::NaN
        peak_rss_mb = [double]::NaN
        fmm_m2l_mode = ''
        fmm_local_gmres_active = 'false'
        fmm_local_operator_finite = ''
        fmm_local_plan_build_count = 0
        fmm_local_cache_build_count = 0
        fmm_local_gmres_matvec_calls = 0
        preconditioner_effective = ''
        reaction_energy_source = ''
        run_dir = $RunDir
    }
}

function Write-Z4xArgFile([object]$Fixture, [object]$Backend, [string]$RunDir, [string]$PqrPath) {
    $isP4 = $Backend.backend -eq 'fmm_p4_guarded'
    $isTree = $Backend.mode -eq 'treecode'
    $argPath = Write-FmmCliArgFile `
        -RunDir $RunDir `
        -RepoRoot $RuntimeRepoRoot `
        -PqrPath $PqrPath `
        -Backend tabi `
        -Mode $Backend.mode `
        -FmmMode $(if ($isTree) { 'single-matvec-analytic' } else { 'off' }) `
        -SurfacePoints $Fixture.surface_points `
        -MaxIter 720 `
        -Restart 128 `
        -Tolerance 1.0e-6 `
        -Threads $Threads `
        -SaltM $SaltM `
        -InternalDielectric $InternalDielectric `
        -ExternalDielectric $ExternalDielectric `
        -TemperatureK $TemperatureK `
        -OpeningTheta 0.6 `
        -AuditSampleCount $(if ($Backend.require_audit) { 512 } else { 0 }) `
        -AuditMode $(if ($Backend.require_audit) { 'sampled-detail' } else { 'off' }) `
        -Preconditioner $Backend.preconditioner `
        -PreconditionerBlockPanelCap 32 `
        -PreconditionerOverlapPanelCap 64 `
        -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular `
        -TaylorOrder 4 `
        -FmmCartesianAudit off `
        -FmmLocal $Backend.fmm_local `
        -FmmLocalRefine $Backend.fmm_refine `
        -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache $Backend.fmm_cache `
        -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback $Backend.fmm_fallback `
        -FmmLocalPlanBuilder $Backend.fmm_plan_builder `
        -FmmLocalRefineProbeCount 4 `
        -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard $Backend.geometry_guard `
        -FmmLocalGeometryTheta $(if ($isP4) { 0.5 } else { 0.0 })
    $bare = @(Get-Content -LiteralPath $argPath | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count -ne 0) {
        throw "Generated arg.in contains bare options: $($bare -join ', ')"
    }
    return $argPath
}

function Read-Z4xResult([object]$Fixture, [object]$Backend, [string]$RunDir, [object]$Run) {
    $row = New-Z4xEmptyRow -Fixture $Fixture -Backend $Backend -RunDir $RunDir
    $row.exit_code = $Run.exit_code
    $row.timed_out = [string]$Run.timed_out
    $row.wall_seconds = $Run.wall_seconds
    $row.peak_rss_mb = $Run.peak_working_set_mb
    if ($Run.timed_out) {
        $row.status = 'timeout'
        $row.failure_class = 'runtime_limited'
        return $row
    }
    if ($Run.exit_code -ne 0) {
        $row.status = 'failed'
        $row.failure_class = "exit_code_$($Run.exit_code)"
        return $row
    }
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $energyPath) -or
        -not (Test-Path -LiteralPath $profilePath)) {
        $row.status = 'read_error'
        $row.failure_class = 'blocked_profile_schema'
        return $row
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    $energyStatus = [string](Get-Z4xProfileValue $energy 'status' 'unknown')
    $profileStatus = [string](Get-Z4xProfileValue $profile 'status' 'unknown')
    $row.status = if ($energyStatus -eq 'ok' -and $profileStatus -eq 'ok') {
        'ok'
    } else {
        "energy=$energyStatus;profile=$profileStatus"
    }
    $row.energy_kj_mol = Get-Z4xDouble (Get-Z4xProfileValue $energy 'E_rxn_kJmol' '')
    $row.reaction_energy_source = [string](Get-Z4xProfileValue $energy 'reaction_energy_source' '')
    $row.panels = Get-Z4xInt64 (Get-Z4xProfileValue $profile 'panels' 0)
    $row.iterations = Get-Z4xInt64 (Get-Z4xProfileValue $profile 'iterations' 0)
    $row.final_residual = Get-Z4xDouble (Get-Z4xProfileValue $profile 'final_residual' '')
    $row.true_residual = Get-Z4xDouble (Get-Z4xProfileValue $profile 'true_residual_rel' '')
    $row.sampled_audit_rms = Get-Z4xDouble (Get-Z4xProfileValue $profile 'tree_vs_direct_matvec_rel' '')
    $row.sampled_direct_true_residual = Get-Z4xDouble (Get-Z4xProfileValue $profile 'fmm_local_sampled_direct_true_residual' '')
    $row.matvec_audit_rows = Get-Z4xInt64 (Get-Z4xProfileValue $profile 'matvec_audit_rows' 0)
    $row.matvec_seconds_avg = Get-Z4xDouble (Get-Z4xProfileValue $profile 'matvec_seconds_avg' '')
    $row.fmm_m2l_mode = [string](Get-Z4xProfileValue $profile 'fmm_m2l_mode' '')
    $row.fmm_local_gmres_active = [string](Get-Z4xProfileValue $profile 'fmm_local_gmres_active' 'false')
    $row.fmm_local_operator_finite = [string](Get-Z4xProfileValue $profile 'fmm_local_operator_finite' '')
    $row.fmm_local_plan_build_count = Get-Z4xInt64 (Get-Z4xProfileValue $profile 'fmm_local_plan_build_count' 0)
    $row.fmm_local_cache_build_count = Get-Z4xInt64 (Get-Z4xProfileValue $profile 'fmm_local_cache_build_count' 0)
    $row.fmm_local_gmres_matvec_calls = Get-Z4xInt64 (Get-Z4xProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
    $row.preconditioner_effective = [string](Get-Z4xProfileValue $profile 'preconditioner_effective' '')
    if ((Test-FmmCliFinite ([double]$row.energy_kj_mol)) -and
        (Test-FmmCliFinite ([double]$row.oracle_energy_kj_mol))) {
        $row.oracle_relative_error =
            [Math]::Abs([double]$row.energy_kj_mol - [double]$row.oracle_energy_kj_mol) /
            [Math]::Max(1.0e-12, [Math]::Abs([double]$row.oracle_energy_kj_mol))
    }
    $row.failure_class = 'pending'
    return $row
}

function Invoke-Z4xCase([object]$Fixture, [object]$Backend) {
    $caseName = Get-FmmCliSafeName "$($Fixture.fixture)_sp$($Fixture.surface_points)_$($Backend.backend)"
    $runDir = Join-Path $WorkRoot $caseName
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    $pqrPath = Join-Path $runDir 'input.pqr'
    Write-Z4xPqr -Fixture $Fixture -Path $pqrPath
    $argPath = Write-Z4xArgFile -Fixture $Fixture -Backend $Backend -RunDir $runDir -PqrPath $pqrPath
    Write-Host "run: $caseName (timeout ${CaseTimeoutSeconds}s)"
    $run = Invoke-FmmCliTitanCase `
        -TitanExe $TitanExe `
        -RunDir $runDir `
        -ArgPath $argPath `
        -CaseTimeoutSeconds $CaseTimeoutSeconds
    return Read-Z4xResult -Fixture $Fixture -Backend $Backend -RunDir $runDir -Run $run
}

function Set-Z4xFailureClass([object]$Row) {
    if ($Row.status -eq 'timeout') { $Row.failure_class = 'runtime_limited'; return $Row }
    if ($Row.status -eq 'read_error') { $Row.failure_class = 'blocked_profile_schema'; return $Row }
    if ($Row.status -ne 'ok') { $Row.failure_class = 'status_fail'; return $Row }
    if (-not (Test-FmmCliFinite ([double]$Row.energy_kj_mol))) {
        $Row.failure_class = 'nonfinite_energy'
        return $Row
    }
    if ($Row.backend -eq 'tabi_direct_dense') {
        if (-not (Test-FmmCliFinite ([double]$Row.oracle_relative_error)) -or
            [double]$Row.oracle_relative_error -gt [double]$Row.direct_oracle_threshold) {
            $Row.failure_class = 'direct_oracle_fail'
            return $Row
        }
        $Row.failure_class = 'pass'
        return $Row
    }
    if ($Row.backend -eq 'fmm_p4_guarded') {
        if (-not (Test-FmmCliFinite ([double]$Row.direct_dense_energy_kj_mol))) {
            $Row.failure_class = 'direct_reference_missing'
            return $Row
        }
        if (-not (Test-FmmCliFinite ([double]$Row.sampled_audit_rms)) -or
            [double]$Row.sampled_audit_rms -gt 8.0e-3) {
            $Row.failure_class = 'sampled_audit_fail'
            return $Row
        }
        if (-not (Test-FmmCliFinite ([double]$Row.sampled_direct_true_residual)) -or
            [double]$Row.sampled_direct_true_residual -gt 8.0e-3) {
            $Row.failure_class = 'sampled_direct_residual_fail'
            return $Row
        }
        if ([Math]::Abs([double]$Row.p4_vs_direct_delta_kj_mol) -gt
            [double]$Row.p4_energy_gate_kj_mol) {
            $Row.failure_class = 'p4_energy_delta_fail'
            return $Row
        }
        if ([int64]$Row.fmm_local_plan_build_count -ne 1 -or
            [int64]$Row.fmm_local_cache_build_count -ne 1 -or
            [int64]$Row.fmm_local_gmres_matvec_calls -le 0 -or
            $Row.fmm_local_gmres_active -ne 'true') {
            $Row.failure_class = 'cache_lifecycle_fail'
            return $Row
        }
        if ($Row.fmm_local_operator_finite -ne 'true') {
            $Row.failure_class = 'operator_finite_fail'
            return $Row
        }
        $Row.failure_class = 'pass'
        return $Row
    }
    if (-not (Test-FmmCliFinite ([double]$Row.sampled_audit_rms)) -or
        [double]$Row.sampled_audit_rms -gt 8.0e-3) {
        $Row.failure_class = 'context_audit_gap'
    } else {
        $Row.failure_class = 'context_only'
    }
    return $Row
}

function Write-Z4xSvg([object[]]$Rows, [string]$Path) {
    $points = @($Rows | Where-Object {
        $_.backend -in @('tabi_direct_dense','fmm_p4_guarded') -and
        (Test-FmmCliFinite ([double]$_.oracle_relative_error)) -and
        [double]$_.oracle_relative_error -gt 0.0
    })
    if ($points.Count -eq 0) {
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="460"><text x="24" y="42">No finite Z4X oracle rows</text></svg>' |
            Set-Content -LiteralPath $Path -Encoding ASCII
        return
    }
    $maxY = [Math]::Max(1.0e-6, ($points | ForEach-Object {
                [Math]::Log10([double]$_.oracle_relative_error)
            } | Measure-Object -Maximum).Maximum)
    $minY = ($points | ForEach-Object {
            [Math]::Log10([double]$_.oracle_relative_error)
        } | Measure-Object -Minimum).Minimum
    if ([Math]::Abs($maxY - $minY) -lt 1.0e-9) { $maxY += 1.0; $minY -= 1.0 }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="980" height="520" viewBox="0 0 980 520">')
    $svg.Add('<rect width="980" height="520" fill="white"/>')
    $svg.Add('<line x1="86" y1="420" x2="930" y2="420" stroke="#222"/>')
    $svg.Add('<line x1="86" y1="58" x2="86" y2="420" stroke="#222"/>')
    $svg.Add('<text x="22" y="282" font-family="Arial" font-size="14" transform="rotate(-90 22,282)">log10 relative error vs analytic oracle</text>')
    $barW = [Math]::Max(9.0, 760.0 / [Math]::Max(1, $points.Count) - 5.0)
    for ($i = 0; $i -lt $points.Count; ++$i) {
        $row = $points[$i]
        $log = [Math]::Log10([double]$row.oracle_relative_error)
        $x = 96 + $i * ($barW + 5.0)
        $y = 58 + (362.0 * (($maxY - $log) / ($maxY - $minY)))
        $h = 420.0 - $y
        $color = if ($row.backend -eq 'tabi_direct_dense') { '#0072B2' } else { '#D55E00' }
        $svg.Add(('<rect x="{0:F1}" y="{1:F1}" width="{2:F1}" height="{3:F1}" fill="{4}"/>' -f $x,$y,$barW,$h,$color))
        $label = "$($row.fixture):$($row.surface_points):$($row.backend -replace 'tabi_direct_dense','direct' -replace 'fmm_p4_guarded','p4')"
        $svg.Add(('<text x="{0:F1}" y="438" font-family="Arial" font-size="9" transform="rotate(50 {0:F1},438)">{1}</text>' -f $x,$label))
    }
    $svg.Add('<text x="650" y="38" font-family="Arial" font-size="12" fill="#0072B2">blue: dense direct TABI</text>')
    $svg.Add('<text x="650" y="56" font-family="Arial" font-size="12" fill="#D55E00">orange: guarded p4 FMM/M2L</text>')
    $svg.Add('</svg>')
    $svg | Set-Content -LiteralPath $Path -Encoding ASCII
}

$surfacePoints = if ($Quick) { @(8) } else { @(8,12,24) }
$backends = if ($Quick) {
    @('tabi_direct_dense','fmm_p4_guarded')
} else {
    @('tabi_direct_dense','fmm_p4_guarded','fmm_degree1_ilu0')
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($sp in $surfacePoints) {
    foreach ($fixtureName in @('born_centered','kirkwood_offcenter')) {
        $fixture = New-Z4xFixture -Fixture $fixtureName -SurfacePoints $sp
        foreach ($backendName in $backends) {
            $backend = New-Z4xBackend -Name $backendName
            $rows.Add((Invoke-Z4xCase -Fixture $fixture -Backend $backend)) | Out-Null
        }
    }
}

$allRows = @($rows.ToArray())
foreach ($row in $allRows) {
    $direct = @($allRows | Where-Object {
            $_.fixture -eq $row.fixture -and
            $_.surface_points -eq $row.surface_points -and
            $_.backend -eq 'tabi_direct_dense'
        } | Select-Object -First 1)
    if ($direct.Count -eq 1 -and
        (Test-FmmCliFinite ([double]$direct[0].energy_kj_mol)) -and
        (Test-FmmCliFinite ([double]$row.energy_kj_mol))) {
        $row.direct_dense_energy_kj_mol = [double]$direct[0].energy_kj_mol
        $row.p4_vs_direct_delta_kj_mol = [double]$row.energy_kj_mol - [double]$direct[0].energy_kj_mol
        $row.p4_vs_direct_relative_error =
            [Math]::Abs([double]$row.p4_vs_direct_delta_kj_mol) /
            [Math]::Max(1.0e-12, [Math]::Abs([double]$direct[0].energy_kj_mol))
        $row.p4_energy_gate_kj_mol =
            [Math]::Max(0.5, 0.01 * [Math]::Max(1.0, [Math]::Abs([double]$direct[0].energy_kj_mol)))
    }
}

$finalRows = @($allRows | ForEach-Object { Set-Z4xFailureClass -Row $_ })
$finalRows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation

$finalRows |
    Where-Object { $_.backend -in @('tabi_direct_dense','fmm_p4_guarded') } |
    Select-Object fixture,fixture_class,surface_points,backend,status,failure_class,
        energy_kj_mol,oracle_energy_kj_mol,direct_dense_energy_kj_mol,
        oracle_relative_error,p4_vs_direct_delta_kj_mol,
        p4_vs_direct_relative_error,p4_energy_gate_kj_mol |
    Export-Csv -LiteralPath $EnergyParityCsv -NoTypeInformation

$finalRows |
    Where-Object { $_.backend -ne 'tabi_direct_dense' } |
    Select-Object fixture,surface_points,backend,status,failure_class,
        sampled_audit_rms,sampled_direct_true_residual,matvec_audit_rows,
        fmm_local_plan_build_count,fmm_local_cache_build_count,
        fmm_local_gmres_matvec_calls,fmm_local_operator_finite |
    Export-Csv -LiteralPath $OperatorAuditCsv -NoTypeInformation

Write-Z4xSvg -Rows $finalRows -Path $PlotPath

$requiredRows = @($finalRows | Where-Object { $_.backend -in @('tabi_direct_dense','fmm_p4_guarded') })
$directRows = @($requiredRows | Where-Object { $_.backend -eq 'tabi_direct_dense' })
$p4Rows = @($requiredRows | Where-Object { $_.backend -eq 'fmm_p4_guarded' })
$timeouts = @($requiredRows | Where-Object { $_.failure_class -eq 'runtime_limited' })
$harnessGaps = @($requiredRows | Where-Object { $_.failure_class -in @('blocked_profile_schema','status_fail','nonfinite_energy') })
$directOracleGaps = @($directRows | Where-Object { $_.failure_class -eq 'direct_oracle_fail' })
$p4Gaps = @($p4Rows | Where-Object { $_.failure_class -notin @('pass') })

$decision = if ($timeouts.Count -gt 0) {
    'z4x_runtime_limited'
} elseif ($harnessGaps.Count -gt 0) {
    'z4x_harness_gap'
} elseif ($directOracleGaps.Count -gt 0) {
    'z4x_surface_oracle_gap'
} elseif ($p4Gaps.Count -gt 0) {
    'z4x_p4_dense_energy_gap'
} else {
    'z4x_dense_oracle_pass'
}

[pscustomobject]@{
    decision = $decision
    quick = [bool]$Quick
    rows = $finalRows.Count
    required_rows = $requiredRows.Count
    direct_rows = $directRows.Count
    p4_rows = $p4Rows.Count
    timeout_rows = $timeouts.Count
    direct_oracle_gap_rows = $directOracleGaps.Count
    p4_gap_rows = $p4Gaps.Count
    case_timeout_seconds = $CaseTimeoutSeconds
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
    'scope: Z4X-A dense small-system oracle gate',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "surface_points: $($surfacePoints -join ',')",
    "backends: $($backends -join ',')",
    "rows: $($finalRows.Count)",
    'fixtures: centered Born ion, off-center Kirkwood sphere',
    'purpose: calibrate guarded p4 sampled/full-protein audit policy against feasible dense direct TABI oracle rows',
    'production_default: grid',
    'tabi_treecode_default: degree-one',
    'gmres_p4: explicit advanced path only'
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII

$report = New-Object System.Collections.Generic.List[string]
$report.Add('# TITAN FMM/M2L Scope Z4X-A Dense Small-System Oracle Report')
$report.Add('')
$report.Add("Decision: $decision")
$report.Add('')
$report.Add('Purpose')
$report.Add('Calibrate the guarded p4 FMM/M2L path against full dense TABI direct solves on analytic small systems where a full dense oracle is feasible.')
$report.Add('')
$report.Add('Decision counts')
$report.Add("required_rows: $($requiredRows.Count)")
$report.Add("timeout_rows: $($timeouts.Count)")
$report.Add("direct_oracle_gap_rows: $($directOracleGaps.Count)")
$report.Add("p4_gap_rows: $($p4Gaps.Count)")
$report.Add('')
$report.Add('Required rows')
foreach ($row in $requiredRows) {
    $report.Add(('- {0} sp={1} {2}: status={3}, failure={4}, oracle_rel={5:G6}, p4_delta={6:G6}, audit={7:G6}, direct_residual={8:G6}, plan/cache={9}/{10}' -f
            $row.fixture,
            $row.surface_points,
            $row.backend,
            $row.status,
            $row.failure_class,
            [double]$row.oracle_relative_error,
            [double]$row.p4_vs_direct_delta_kj_mol,
            [double]$row.sampled_audit_rms,
            [double]$row.sampled_direct_true_residual,
            [int64]$row.fmm_local_plan_build_count,
            [int64]$row.fmm_local_cache_build_count))
}
$report.Add('')
$report.Add('Artifacts')
$report.Add("Matrix: $MatrixCsv")
$report.Add("Energy parity: $EnergyParityCsv")
$report.Add("Operator audit: $OperatorAuditCsv")
$report.Add("Decision: $DecisionCsv")
$report.Add("Plot: $PlotPath")
$report.Add('')
$report.Add('No production defaults are changed by this gate. Grid remains the production PBE default; guarded p4 remains explicit.')
$report | Set-Content -LiteralPath $RootReport -Encoding ASCII

Write-Host "Z4X-A decision: $decision"
Write-Host "Matrix: $MatrixCsv"

if ($decision -ne 'z4x_dense_oracle_pass' -and -not $AllowFailedDecision) {
    throw "Z4X-A dense small-system oracle gate did not pass: $decision"
}
