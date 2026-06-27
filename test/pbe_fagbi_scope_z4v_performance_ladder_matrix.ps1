param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$PrimaryFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$ExpansionFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public_z4q'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4v'),
    [string[]]$Fixtures = @(),
    [Alias('Stages')]
    [string[]]$StageFilter = @(),
    [int]$DirectTimeoutSeconds = 3600,
    [int]$FmmTimeoutSeconds = 1800,
    [int]$Threads = 4,
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
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4v_matrix.csv'
$ProteinCsv = Join-Path $WorkRoot 'fagbi_scope_z4v_summary_by_protein.csv'
$TimeoutCsv = Join-Path $WorkRoot 'fagbi_scope_z4v_timeout_policy.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4v_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4v_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4v_performance_ladder.svg'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4V_PERFORMANCE_LADDER_REPORT.txt'

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

function Get-Z4vDouble([object]$Value, [double]$Default = [double]::NaN) {
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

function Get-Z4vInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $parsed = [int64]0
    if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Format-Z4vDouble([double]$Value, [string]$Format = 'G4') {
    if (-not (Test-FmmCliFinite $Value)) { return 'NA' }
    return $Value.ToString($Format, [Globalization.CultureInfo]::InvariantCulture)
}

function Import-Z4vManifestRows([string]$Root, [string]$ManifestName) {
    $path = Join-Path $Root $ManifestName
    if (-not (Test-Path -LiteralPath $path)) { return @() }
    return @(Import-Csv -LiteralPath $path | ForEach-Object {
        $_ | Add-Member -NotePropertyName fixture_root -NotePropertyValue $Root -Force
        $_
    })
}

function Get-Z4vFixtureRows {
    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($row in (Import-Z4vManifestRows $PrimaryFixtureRoot 'public_fixture_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    foreach ($row in (Import-Z4vManifestRows $ExpansionFixtureRoot 'public_fixture_expansion_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    return $rows.ToArray()
}

function Get-Z4vFixture([string]$Name, [object[]]$Rows) {
    $match = @($Rows | Where-Object { $_.fixture -eq $Name })
    if ($match.Count -ne 1) {
        throw "Expected exactly one cached public fixture row for '$Name', found $($match.Count)."
    }
    $row = $match[0]
    $pqrPath = Join-Path $row.fixture_root $row.pqr
    if (-not (Test-Path -LiteralPath $pqrPath)) {
        throw "Missing cached PQR for $Name`: $pqrPath"
    }
    $row | Add-Member -NotePropertyName pqr_path -NotePropertyValue $pqrPath -Force
    return $row
}

function Get-Z4vProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function New-Z4vRow(
    [object]$Fixture,
    [object]$Stage,
    [int]$ExitCode,
    [bool]$TimedOut,
    [double]$WallSeconds,
    [double]$PeakRssMb,
    [string]$Status,
    [string]$FailureClass,
    [string]$RunDir
) {
    return [pscustomobject]@{
        fixture = $Fixture.fixture
        protein_name = $Fixture.protein_name
        pdb_id = $Fixture.pdb_id
        atoms = Get-Z4vInt64 $Fixture.atom_count
        stage = $Stage.stage
        stage_order = $Stage.order
        backend = $Stage.backend
        fmm_mode = $Stage.fmm_mode
        preconditioner = $Stage.preconditioner
        p4_local = $Stage.p4_local
        timeout_seconds = $Stage.timeout
        exit_code = $ExitCode
        timed_out = [string]$TimedOut
        status = $Status
        failure_class = $FailureClass
        panels = 0
        iterations = 0
        energy_kj_mol = [double]::NaN
        final_residual = [double]::NaN
        true_residual = [double]::NaN
        sampled_audit_rms = [double]::NaN
        sampled_direct_true_residual = [double]::NaN
        matvec_seconds_avg = [double]::NaN
        wall_seconds = $WallSeconds
        peak_rss_mb = $PeakRssMb
        fmm_near_pairs = 0
        fmm_far_pairs = 0
        fmm_finite_outputs = ''
        preconditioner_effective = ''
        plan_build_count = 0
        cache_build_count = 0
        gmres_matvec_calls = 0
        operator_finite = ''
        speedup_vs_tabi_direct = [double]::NaN
        speedup_vs_degree1_diagonal = [double]::NaN
        speedup_vs_fmm_diagonal = [double]::NaN
        paper_usable = 'false'
        run_dir = $RunDir
    }
}

function Read-Z4vCaseOutputs([string]$RunDir, [object]$Row) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    if (-not (Test-Path -LiteralPath $energyPath)) {
        $Row.failure_class = 'missing_energy_output'
        return $Row
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $Row.energy_kj_mol = Get-Z4vDouble $energy.E_rxn_kJmol
    if ([string]$energy.status -ne 'ok') {
        $Row.status = [string]$energy.status
        $Row.failure_class = 'energy_status_gap'
    }

    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (Test-Path -LiteralPath $profilePath) {
        $profile = @(Import-Csv -LiteralPath $profilePath)[0]
        $Row.status = if ($Row.status -eq 'ok' -and [string]$profile.status -eq 'ok') {
            'ok'
        } elseif ($Row.status -eq 'ok') {
            "profile=$($profile.status)"
        } else {
            $Row.status
        }
        $Row.panels = Get-Z4vInt64 $profile.panels
        $Row.iterations = Get-Z4vInt64 $profile.iterations
        $Row.final_residual = Get-Z4vDouble $profile.final_residual
        $Row.true_residual = Get-Z4vDouble $profile.true_residual_rel
        $Row.sampled_audit_rms = Get-Z4vDouble $profile.tree_vs_direct_matvec_rel
        $Row.matvec_seconds_avg = Get-Z4vDouble $profile.matvec_seconds_avg
        $Row.fmm_near_pairs = Get-Z4vInt64 (Get-Z4vProfileValue $profile 'fmm_m2l_near_pairs' 0)
        $Row.fmm_far_pairs = Get-Z4vInt64 (Get-Z4vProfileValue $profile 'fmm_m2l_far_pairs' 0)
        $Row.fmm_finite_outputs = [string](Get-Z4vProfileValue $profile 'fmm_m2l_finite_outputs' '')
        $Row.preconditioner_effective = [string](Get-Z4vProfileValue $profile 'preconditioner_effective' '')
        $Row.plan_build_count = Get-Z4vInt64 (Get-Z4vProfileValue $profile 'fmm_local_plan_build_count' 0)
        $Row.cache_build_count = Get-Z4vInt64 (Get-Z4vProfileValue $profile 'fmm_local_cache_build_count' 0)
        $Row.gmres_matvec_calls = Get-Z4vInt64 (Get-Z4vProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
        $Row.operator_finite = [string](Get-Z4vProfileValue $profile 'fmm_local_operator_finite' '')
        $Row.sampled_direct_true_residual = Get-Z4vDouble (Get-Z4vProfileValue $profile 'fmm_local_sampled_direct_true_residual' '')
    } elseif ($Row.backend -ne 'grid_context') {
        $Row.failure_class = 'missing_tabi_profile'
    }

    if ($Row.status -eq 'ok' -and $Row.failure_class -eq 'pass') {
        if ($Row.backend -like 'fmm_*' -and
            (-not (Test-FmmCliFinite ([double]$Row.sampled_audit_rms)))) {
            $Row.failure_class = 'missing_sampled_audit'
        } elseif ($Row.p4_local -eq 'gmres-p4' -and
            (-not (Test-FmmCliFinite ([double]$Row.sampled_direct_true_residual)))) {
            $Row.failure_class = 'missing_sampled_direct_true_residual'
        } else {
            $Row.paper_usable = 'true'
        }
    }
    return $Row
}

function Get-Z4vStages([int]$DirectTimeout, [int]$FmmTimeout) {
    return @(
        [pscustomobject]@{
            order=1; stage='tabi_direct'; backend='tabi_direct'; mode='direct';
            fmm_mode='off'; preconditioner='diagonal'; p4_local='off';
            audit='off'; audit_count=0; timeout=$DirectTimeout
        },
        [pscustomobject]@{
            order=2; stage='fmm_degree1_diagonal'; backend='fmm_degree1_diagonal'; mode='treecode';
            fmm_mode='single-matvec-analytic'; preconditioner='diagonal'; p4_local='off';
            audit='sampled'; audit_count=512; timeout=$FmmTimeout
        },
        [pscustomobject]@{
            order=3; stage='fmm_degree1_block'; backend='fmm_degree1_block'; mode='treecode';
            fmm_mode='single-matvec-analytic'; preconditioner='fmm-block-jacobi'; p4_local='off';
            audit='sampled'; audit_count=512; timeout=$FmmTimeout
        },
        [pscustomobject]@{
            order=4; stage='fmm_degree1_ilu0'; backend='fmm_degree1_ilu0'; mode='treecode';
            fmm_mode='single-matvec-analytic'; preconditioner='fmm-nearfield-ilu0'; p4_local='off';
            audit='sampled'; audit_count=512; timeout=$FmmTimeout
        },
        [pscustomobject]@{
            order=5; stage='fmm_p4_ilu0'; backend='fmm_p4_ilu0'; mode='treecode';
            fmm_mode='single-matvec-analytic'; preconditioner='fmm-nearfield-ilu0'; p4_local='gmres-p4';
            audit='sampled-detail'; audit_count=512; timeout=$FmmTimeout
        }
    )
}

function Write-Z4vArg([object]$Fixture, [object]$Stage, [string]$RunDir) {
    $isP4 = $Stage.p4_local -eq 'gmres-p4'
    return Write-FmmCliArgFile -RunDir $RunDir -RepoRoot $RuntimeRepoRoot `
        -PqrPath $Fixture.pqr_path -Backend tabi -Mode $Stage.mode `
        -FmmMode $Stage.fmm_mode -SurfacePoints 8 -MaxIter 720 -Restart 128 `
        -Tolerance 1e-6 -Threads $Threads -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount $Stage.audit_count -AuditMode $Stage.audit `
        -Preconditioner $Stage.preconditioner -PreconditionerBlockPanelCap 32 `
        -PreconditionerOverlapPanelCap 64 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder 4 -FmmCartesianAudit off -FmmLocal $Stage.p4_local `
        -FmmLocalRefine $(if ($isP4) { 'operator-tail-static' } else { 'off' }) `
        -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache $(if ($isP4) { 'green-jet-compressed' } else { 'off' }) `
        -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback $(if ($isP4) { 'packed-avx2' } else { 'scalar' }) `
        -FmmLocalPlanBuilder $(if ($isP4) { 'metric-cache' } else { 'legacy-rescan' }) `
        -FmmLocalRefineProbeCount 4 -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard $(if ($isP4) { 'strict-geometry-mac' } else { 'off' }) `
        -FmmLocalGeometryTheta $(if ($isP4) { 0.5 } else { 0.0 })
}

function Test-Z4vRowFinal([object]$Row) {
    return $Row.status -eq 'ok' -and $Row.failure_class -eq 'pass'
}

function Test-Z4vShouldRun([string]$RowCsv) {
    if ($Force) { return $true }
    if (-not (Test-Path -LiteralPath $RowCsv)) { return $true }
    try {
        $row = @(Import-Csv -LiteralPath $RowCsv)[0]
        if ($row.status -eq 'planned' -or $row.failure_class -eq 'plan_only') {
            return $true
        }
        if (-not $RerunFailed) { return $false }
        return -not (Test-Z4vRowFinal $row)
    } catch {
        return $true
    }
}

function Invoke-Z4vCase([object]$Fixture, [object]$Stage) {
    $runName = '{0}_{1}' -f $Fixture.fixture, $Stage.stage
    $runDir = Join-Path $WorkRoot $runName
    $rowCsv = Join-Path $runDir 'fagbi_scope_z4v_row.csv'
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null

    if (-not (Test-Z4vShouldRun $rowCsv)) {
        Write-Host "resume: $runName"
        return @(Import-Csv -LiteralPath $rowCsv)[0]
    }

    Write-Host "run: $runName (timeout $($Stage.timeout)s)"
    $argPath = Write-Z4vArg -Fixture $Fixture -Stage $Stage -RunDir $runDir
    $bare = @(Get-Content -LiteralPath $argPath | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) {
        throw "Generated arg.in contains bare options: $($bare -join ', ')"
    }
    if ($PlanOnly) {
        $row = New-Z4vRow -Fixture $Fixture -Stage $Stage -ExitCode 0 `
            -TimedOut $false -WallSeconds 0.0 -PeakRssMb 0.0 `
            -Status 'planned' -FailureClass 'plan_only' -RunDir $runDir
        $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
        return $row
    }

    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $Stage.timeout
    $row = New-Z4vRow -Fixture $Fixture -Stage $Stage -ExitCode $run.exit_code `
        -TimedOut $run.timed_out -WallSeconds $run.wall_seconds `
        -PeakRssMb $run.peak_working_set_mb -Status 'failed' `
        -FailureClass 'pass' -RunDir $runDir
    if ($run.timed_out) {
        $row.status = 'timeout'
        $row.failure_class = if ($Stage.stage -eq 'tabi_direct') {
            'direct_timeout_not_verified'
        } else {
            'runtime_limited'
        }
    } elseif ($run.exit_code -ne 0) {
        $stdoutPath = Join-Path $runDir 'stdout.log'
        $stderrPath = Join-Path $runDir 'stderr.log'
        $log = ''
        if (Test-Path -LiteralPath $stdoutPath) { $log += Get-Content -Raw -LiteralPath $stdoutPath }
        if (Test-Path -LiteralPath $stderrPath) { $log += Get-Content -Raw -LiteralPath $stderrPath }
        if ($log -match 'gmres-p4 is fail-closed') {
            $row.status = 'unsupported'
            $row.failure_class = 'unsupported_current_binary'
        } else {
            $row.status = 'failed'
            $row.failure_class = "exit_code_$($run.exit_code)"
        }
    } else {
        $row.status = 'ok'
        $row = Read-Z4vCaseOutputs -RunDir $runDir -Row $row
    }
    $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
    return $row
}

function Add-Z4vSpeedups([object[]]$Rows) {
    foreach ($fixture in @($Rows.fixture | Sort-Object -Unique)) {
        $direct = @($Rows | Where-Object { $_.fixture -eq $fixture -and $_.stage -eq 'tabi_direct' } | Select-Object -First 1)
        $diag = @($Rows | Where-Object { $_.fixture -eq $fixture -and $_.stage -eq 'fmm_degree1_diagonal' } | Select-Object -First 1)
        $directMatvec = if ($direct.Count) { Get-Z4vDouble $direct[0].matvec_seconds_avg } else { [double]::NaN }
        $diagMatvec = if ($diag.Count) { Get-Z4vDouble $diag[0].matvec_seconds_avg } else { [double]::NaN }
        foreach ($row in @($Rows | Where-Object { $_.fixture -eq $fixture })) {
            $matvec = Get-Z4vDouble $row.matvec_seconds_avg
            if ((Test-FmmCliFinite $directMatvec) -and $directMatvec -gt 0 -and
                (Test-FmmCliFinite $matvec) -and $matvec -gt 0) {
                $row.speedup_vs_tabi_direct = $directMatvec / $matvec
            }
            if ((Test-FmmCliFinite $diagMatvec) -and $diagMatvec -gt 0 -and
                (Test-FmmCliFinite $matvec) -and $matvec -gt 0) {
                $row.speedup_vs_degree1_diagonal = $diagMatvec / $matvec
            }
        }
    }
    return $Rows
}

function Write-Z4vSvg([object[]]$Rows, [string]$Path) {
    $fixtures = @('barnase','1a63','1sh1','1ubq','2ci2','1vii') | Where-Object {
        $fx = $_
        @($Rows | Where-Object { $_.fixture -eq $fx }).Count -gt 0
    }
    $stageNames = @('tabi_direct','fmm_degree1_diagonal','fmm_degree1_block','fmm_degree1_ilu0','fmm_p4_ilu0')
    $colors = @{
        tabi_direct = '#6b7280'; fmm_degree1_diagonal = '#60a5fa';
        fmm_degree1_block = '#f59e0b'; fmm_degree1_ilu0 = '#16a34a';
        fmm_p4_ilu0 = '#7c3aed'
    }
    $stageLabels = @{
        tabi_direct = 'Direct TABI'
        fmm_degree1_diagonal = 'FMM degree-one + diagonal'
        fmm_degree1_block = 'FMM degree-one + local block'
        fmm_degree1_ilu0 = 'FMM degree-one + sparse ILU(0)'
        fmm_p4_ilu0 = 'Guarded p4 + sparse ILU(0)'
    }
    function Escape-Z4vSvgText([string]$Text) {
        if ($null -eq $Text) { return '' }
        return ($Text -replace '&','&amp;' -replace '<','&lt;' -replace '>','&gt;')
    }
    function Format-Z4vTime([double]$Seconds) {
        if (-not (Test-FmmCliFinite $Seconds)) { return 'NA' }
        if ($Seconds -ge 1000.0) { return ('{0:F1} min' -f ($Seconds / 60.0)) }
        return ('{0:F1} s' -f $Seconds)
    }
    function Get-Z4vFixtureLabel([string]$Fixture) {
        $row = @($Rows | Where-Object { $_.fixture -eq $Fixture -and $_.stage -eq 'fmm_p4_ilu0' } | Select-Object -First 1)
        if ($row.Count -eq 0) { return $Fixture.ToUpperInvariant() }
        $pdb = [string]$row[0].pdb_id
        if ($pdb -eq 'PyDelPhi benchmark') { return 'Barnase' }
        return $pdb.ToUpperInvariant()
    }
    function Add-Z4vText([System.Collections.Generic.List[string]]$Svg, [int]$X, [int]$Y, [string]$Text, [int]$Size, [string]$Color = '#111827', [string]$Weight = 'normal') {
        $Svg.Add("<text x='$X' y='$Y' font-family='Arial' font-size='$Size' font-weight='$Weight' fill='$Color'>$(Escape-Z4vSvgText $Text)</text>") | Out-Null
    }

    $width = 1400
    $height = 210 + 72 * $fixtures.Count
    $left = 36
    $systemW = 145
    $cellW = 225
    $cellH = 56
    $rowGap = 16
    $top = 160
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add("<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>") | Out-Null
    $svg.Add("<rect width='$width' height='$height' fill='white'/>") | Out-Null
    Add-Z4vText $svg 34 34 'Full-protein TABI performance ladder' 24 '#111827' 'bold'
    Add-Z4vText $svg 34 60 'Each colored cell is one backend stage on the same full-protein fixture; text reports status and wall time.' 15 '#374151'
    $legendX = 34
    $legendY = 86
    foreach ($stage in $stageNames) {
        $svg.Add("<rect x='$legendX' y='$legendY' width='18' height='18' rx='3' fill='$($colors[$stage])'/>") | Out-Null
        Add-Z4vText $svg ($legendX + 26) ($legendY + 14) $stageLabels[$stage] 13 '#111827'
        $legendX += 250
    }
    Add-Z4vText $svg $left 135 'System' 14 '#111827' 'bold'
    $headerX = $left + $systemW
    foreach ($stage in $stageNames) {
        Add-Z4vText $svg ($headerX + 6) 135 $stageLabels[$stage] 12 '#111827' 'bold'
        $headerX += $cellW
    }
    $y = $top
    foreach ($fixture in $fixtures) {
        Add-Z4vText $svg $left ($y + 22) (Get-Z4vFixtureLabel $fixture) 15 '#111827' 'bold'
        $systemRow = @($Rows | Where-Object { $_.fixture -eq $fixture -and $_.stage -eq 'fmm_p4_ilu0' } | Select-Object -First 1)
        if ($systemRow.Count) {
            Add-Z4vText $svg $left ($y + 43) ("{0} panels" -f $systemRow[0].panels) 12 '#4b5563'
        }
        $x = $left + $systemW
        foreach ($stage in $stageNames) {
            $row = @($Rows | Where-Object { $_.fixture -eq $fixture -and $_.stage -eq $stage } | Select-Object -First 1)
            $wall = if ($row.Count) { Get-Z4vDouble $row[0].wall_seconds } else { [double]::NaN }
            $color = $colors[$stage]
            $svg.Add("<rect x='$x' y='$y' width='$($cellW - 12)' height='$cellH' rx='7' fill='$color' opacity='0.92'/>") | Out-Null
            if ($row.Count) {
                $status = [string]$row[0].status
                $firstLine = if ($status -eq 'ok') { 'ok' } elseif ($status -eq 'timeout') { 'timeout' } else { 'not converged' }
                $secondLine = Format-Z4vTime $wall
                if ($stage -eq 'fmm_p4_ilu0') {
                    $audit = Get-Z4vDouble $row[0].sampled_audit_rms
                    $thirdLine = if (Test-FmmCliFinite $audit) { 'audit ' + (Format-Z4vDouble $audit '0.##E+0') } else { 'audit NA' }
                } elseif ($stage -eq 'fmm_degree1_ilu0') {
                    $audit = Get-Z4vDouble $row[0].sampled_audit_rms
                    $thirdLine = if (Test-FmmCliFinite $audit) { 'audit ' + (Format-Z4vDouble $audit '0.##E+0') } else { '' }
                } else {
                    $iters = Get-Z4vDouble $row[0].iterations
                    $thirdLine = if (Test-FmmCliFinite $iters -and $iters -gt 0) { ('iter {0}' -f ([int]$iters)) } else { '' }
                }
                Add-Z4vText $svg ($x + 10) ($y + 18) $firstLine 14 'white' 'bold'
                Add-Z4vText $svg ($x + 10) ($y + 36) $secondLine 13 'white'
                if ($thirdLine) { Add-Z4vText $svg ($x + 10) ($y + 52) $thirdLine 12 'white' }
            } else {
                Add-Z4vText $svg ($x + 10) ($y + 31) 'missing' 13 'white' 'bold'
            }
            $x += $cellW
        }
        $y += ($cellH + $rowGap)
    }
    $svg.Add('</svg>') | Out-Null
    $svg | Set-Content -LiteralPath $Path -Encoding ASCII
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$fixtureRows = Get-Z4vFixtureRows
$defaultFixtures = @(if ($Quick) {
    @('1sh1')
} elseif ($Fixtures.Count -gt 0) {
    $Fixtures
} else {
    @('barnase','1a63','1sh1','1ubq','2ci2','1vii')
})
$allStages = Get-Z4vStages -DirectTimeout $DirectTimeoutSeconds -FmmTimeout $FmmTimeoutSeconds
$selectedStages = @(if ($StageFilter.Count -gt 0) {
    @($allStages | Where-Object { $StageFilter -contains $_.stage })
} else {
    $allStages
})
if ($selectedStages.Count -eq 0) {
    throw "No matching Z4V stages selected. Valid stages: $(@($allStages.stage) -join ', ')"
}
$rows = New-Object System.Collections.Generic.List[object]

foreach ($fixtureName in $defaultFixtures) {
    $fixture = Get-Z4vFixture -Name $fixtureName -Rows $fixtureRows
    foreach ($stage in $selectedStages) {
        $rows.Add((Invoke-Z4vCase -Fixture $fixture -Stage $stage)) | Out-Null
    }
}

$allRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $defaultFixtures) {
    foreach ($stage in $allStages) {
        $rowCsv = Join-Path (Join-Path $WorkRoot ('{0}_{1}' -f $fixtureName, $stage.stage)) 'fagbi_scope_z4v_row.csv'
        if (Test-Path -LiteralPath $rowCsv) {
            $allRows.Add((@(Import-Csv -LiteralPath $rowCsv)[0])) | Out-Null
        }
    }
}
$matrix = @(Add-Z4vSpeedups @($allRows.ToArray()))
$matrix | Sort-Object fixture,stage_order | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation

$proteinRows = New-Object System.Collections.Generic.List[object]
foreach ($fixture in @($matrix.fixture | Sort-Object -Unique)) {
    $set = @($matrix | Where-Object fixture -eq $fixture)
    $p4 = @($set | Where-Object stage -eq 'fmm_p4_ilu0' | Select-Object -First 1)
    $direct = @($set | Where-Object stage -eq 'tabi_direct' | Select-Object -First 1)
    $fmmIlu = @($set | Where-Object stage -eq 'fmm_degree1_ilu0' | Select-Object -First 1)
    $proteinRows.Add([pscustomobject]@{
        fixture = $fixture
        protein_name = ($set | Select-Object -First 1).protein_name
        pdb_id = ($set | Select-Object -First 1).pdb_id
        atoms = ($set | Select-Object -First 1).atoms
        panels = if ($p4.Count) { $p4[0].panels } else { ($set | Measure-Object panels -Maximum).Maximum }
        direct_status = if ($direct.Count) { $direct[0].status } else { 'missing' }
        direct_wall_seconds = if ($direct.Count) { $direct[0].wall_seconds } else { [double]::NaN }
        fmm_ilu_status = if ($fmmIlu.Count) { $fmmIlu[0].status } else { 'missing' }
        fmm_ilu_wall_seconds = if ($fmmIlu.Count) { $fmmIlu[0].wall_seconds } else { [double]::NaN }
        p4_status = if ($p4.Count) { $p4[0].status } else { 'missing' }
        p4_audit = if ($p4.Count) { $p4[0].sampled_audit_rms } else { [double]::NaN }
        p4_direct_true_residual = if ($p4.Count) { $p4[0].sampled_direct_true_residual } else { [double]::NaN }
        p4_wall_seconds = if ($p4.Count) { $p4[0].wall_seconds } else { [double]::NaN }
        p4_peak_rss_mb = if ($p4.Count) { $p4[0].peak_rss_mb } else { [double]::NaN }
    }) | Out-Null
}
$proteinRows.ToArray() | Export-Csv -LiteralPath $ProteinCsv -NoTypeInformation

@(
    [pscustomobject]@{ backend='tabi_direct'; timeout_seconds=$DirectTimeoutSeconds; interpretation='timeout means direct full reference is runtime-limited, not a solver science failure' },
    [pscustomobject]@{ backend='fmm_degree1_diagonal'; timeout_seconds=$FmmTimeoutSeconds; interpretation='FMM without ILU context row' },
    [pscustomobject]@{ backend='fmm_degree1_block'; timeout_seconds=$FmmTimeoutSeconds; interpretation='paper-style/local block preconditioner context row' },
    [pscustomobject]@{ backend='fmm_degree1_ilu0'; timeout_seconds=$FmmTimeoutSeconds; interpretation='isolates sparse ILU0 convergence benefit on degree-one FMM' },
    [pscustomobject]@{ backend='fmm_p4_ilu0'; timeout_seconds=$FmmTimeoutSeconds; interpretation='final guarded p4 advanced explicit path' }
) | Export-Csv -LiteralPath $TimeoutCsv -NoTypeInformation

$requiredP4 = @($matrix | Where-Object stage -eq 'fmm_p4_ilu0')
$p4Pass = @($requiredP4 | Where-Object {
    $_.status -eq 'ok' -and
    (Test-FmmCliFinite (Get-Z4vDouble $_.sampled_audit_rms)) -and
    (Get-Z4vDouble $_.sampled_audit_rms) -le 8.0e-3 -and
    (Test-FmmCliFinite (Get-Z4vDouble $_.sampled_direct_true_residual)) -and
    (Get-Z4vDouble $_.sampled_direct_true_residual) -le 8.0e-3 -and
    (Test-FmmCliFinite (Get-Z4vDouble $_.energy_kj_mol)) -and
    (Get-Z4vDouble $_.peak_rss_mb) -lt 512.0 -and
    (Get-Z4vInt64 $_.plan_build_count) -eq 1 -and
    (Get-Z4vInt64 $_.cache_build_count) -eq 1
})
$timeouts = @($matrix | Where-Object status -eq 'timeout')
$usableStages = @($matrix | Where-Object { $_.paper_usable -eq 'true' })
$decision = if ($requiredP4.Count -eq $defaultFixtures.Count -and $p4Pass.Count -eq $requiredP4.Count) {
    if ($timeouts.Count) { 'z4v_p4_pass_with_context_timeouts' } else { 'z4v_performance_ladder_pass' }
} elseif ($PlanOnly) {
    'z4v_plan_only'
} elseif ($p4Pass.Count -gt 0) {
    'z4v_partial_p4_pass'
} elseif ($timeouts.Count) {
    'z4v_runtime_limited'
} else {
    'z4v_accuracy_or_status_gap'
}

[pscustomobject]@{
    decision = $decision
    quick = [string][bool]$Quick
    fixtures = ($defaultFixtures -join ';')
    rows = $matrix.Count
    paper_usable_rows = $usableStages.Count
    p4_pass_rows = $p4Pass.Count
    timeout_rows = $timeouts.Count
    direct_timeout_seconds = $DirectTimeoutSeconds
    fmm_timeout_seconds = $FmmTimeoutSeconds
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

Write-Z4vSvg -Rows $matrix -Path $PlotPath

@(
    'scope: Z4V full-protein TABI performance ladder matrix',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "plan_only: $([bool]$PlanOnly)",
    "fixtures: $($defaultFixtures -join ', ')",
    "selected_stages: $(@($selectedStages.stage) -join ', ')",
    "rows: $($matrix.Count)",
    "direct_timeout_seconds: $DirectTimeoutSeconds",
    "fmm_timeout_seconds: $FmmTimeoutSeconds",
    'resume_policy: row CSV files under each run directory are reused unless -Force or -RerunFailed is passed',
    'production_default: grid',
    'tabi_treecode_default: degree-one; gmres-p4 explicit only'
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z4V Performance Ladder Report',
    '',
    "Decision: $decision",
    '',
    'Rows compare classical direct TABI, degree-one FMM without strong preconditioning, local block preconditioning, sparse near-field ILU(0), and final guarded p4+ILU0 on matched full-protein fixtures.',
    '',
    "Matrix: $MatrixCsv",
    "Protein summary: $ProteinCsv",
    "Timeout policy: $TimeoutCsv",
    '',
    'Timeout rows are retained as performance evidence and classified separately from scientific accuracy failures.'
)
$report | Set-Content -LiteralPath $RootReport -Encoding ASCII

Write-Host "Z4V decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Resume: rerun the same command to continue missing rows; use -RerunFailed for failed/timeout rows or -Force for all rows."

if ($decision -notin @('z4v_performance_ladder_pass','z4v_p4_pass_with_context_timeouts') -and
    -not $AllowFailedDecision -and -not $PlanOnly) {
    throw "Scope Z4V did not pass: $decision"
}
