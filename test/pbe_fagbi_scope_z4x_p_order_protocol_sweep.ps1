param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$PrimaryFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$ExpansionFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public_z4q'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4xb'),
    [string[]]$Fixtures = @(),
    [int]$CaseTimeoutSeconds = 1800,
    [int]$Threads = 4,
    [switch]$Quick,
    [switch]$PlanOnly,
    [switch]$Force,
    [switch]$RerunFailed,
    [string[]]$RerunStages = @(),
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$PrimaryFixtureRoot = [IO.Path]::GetFullPath($PrimaryFixtureRoot)
$ExpansionFixtureRoot = [IO.Path]::GetFullPath($ExpansionFixtureRoot)
if (-not [IO.Path]::IsPathRooted($WorkRoot)) {
    $WorkRoot = Join-Path $RepoRoot $WorkRoot
}
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)

$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4xb_matrix.csv'
$OrderSummaryCsv = Join-Path $WorkRoot 'fagbi_scope_z4xb_order_summary.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4xb_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4xb_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4xb_p_order_plot.svg'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4XB_P_ORDER_SWEEP_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

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

function Get-Z4xbDouble([object]$Value, [double]$Default = [double]::NaN) {
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

function Get-Z4xbInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $parsed = [int64]0
    if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-Z4xbProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Import-Z4xbManifestRows([string]$Root, [string]$ManifestName) {
    $path = Join-Path $Root $ManifestName
    if (-not (Test-Path -LiteralPath $path)) { return @() }
    return @(Import-Csv -LiteralPath $path | ForEach-Object {
        $_ | Add-Member -NotePropertyName fixture_root -NotePropertyValue $Root -Force
        $_
    })
}

function Get-Z4xbFixtureRows {
    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($row in (Import-Z4xbManifestRows $PrimaryFixtureRoot 'public_fixture_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    foreach ($row in (Import-Z4xbManifestRows $ExpansionFixtureRoot 'public_fixture_expansion_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    return $rows.ToArray()
}

function Get-Z4xbFixture([string]$Name, [object[]]$Rows) {
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

function New-Z4xbStage(
    [int]$Order,
    [string]$Stage,
    [int]$TaylorOrder,
    [string]$UnsupportedReason = ''
) {
    return [pscustomobject]@{
        order = $Order
        stage = $Stage
        taylor_order = $TaylorOrder
        unsupported_reason = $UnsupportedReason
    }
}

function Get-Z4xbStages {
    return @(
        (New-Z4xbStage 1 'p_order2' 2 ''),
        (New-Z4xbStage 2 'p_order3' 3 ''),
        (New-Z4xbStage 3 'p_order4' 4 ''),
        (New-Z4xbStage 4 'p_order5' 5 '')
    )
}

function New-Z4xbRow([object]$Fixture, [object]$Stage, [string]$RunDir) {
    return [pscustomobject]@{
        fixture = $Fixture.fixture
        pdb_id = $Fixture.pdb_id
        protein_name = $Fixture.protein_name
        atoms = Get-Z4xbInt64 $Fixture.atom_count
        stage = $Stage.stage
        stage_order = $Stage.order
        requested_taylor_order = $Stage.taylor_order
        effective_taylor_order = 0
        status = 'planned'
        failure_class = 'planned'
        unsupported_reason = $Stage.unsupported_reason
        exit_code = 0
        timed_out = 'False'
        panels = 0
        iterations = 0
        energy_kj_mol = [double]::NaN
        final_residual = [double]::NaN
        true_residual = [double]::NaN
        sampled_audit_rms = [double]::NaN
        sampled_direct_true_residual = [double]::NaN
        matvec_seconds_avg = [double]::NaN
        wall_seconds = [double]::NaN
        peak_rss_mb = [double]::NaN
        fmm_local_gmres_active = 'false'
        fmm_local_operator_finite = ''
        plan_build_count = 0
        cache_build_count = 0
        gmres_matvec_calls = 0
        geometry_guard = ''
        geometry_theta = [double]::NaN
        run_dir = $RunDir
    }
}

function Write-Z4xbArg([object]$Fixture, [object]$Stage, [string]$RunDir) {
    return Write-FmmCliArgFile `
        -RunDir $RunDir `
        -RepoRoot $RuntimeRepoRoot `
        -PqrPath $Fixture.pqr_path `
        -Backend tabi `
        -Mode treecode `
        -FmmMode single-matvec-analytic `
        -SurfacePoints 8 `
        -MaxIter 720 `
        -Restart 128 `
        -Tolerance 1.0e-6 `
        -Threads $Threads `
        -SaltM 0.15 `
        -InternalDielectric 1 `
        -ExternalDielectric 80 `
        -TemperatureK 298.15 `
        -OpeningTheta 0.6 `
        -AuditSampleCount 512 `
        -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerBlockPanelCap 32 `
        -PreconditionerOverlapPanelCap 64 `
        -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular `
        -TaylorOrder $Stage.taylor_order `
        -FmmCartesianAudit off `
        -FmmLocal gmres-p4 `
        -FmmLocalRefine operator-tail-static `
        -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache green-jet-compressed `
        -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback packed-avx2 `
        -FmmLocalPlanBuilder metric-cache `
        -FmmLocalRefineProbeCount 4 `
        -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard strict-geometry-mac `
        -FmmLocalGeometryTheta 0.5
}

function Read-Z4xbOutputs([string]$RunDir, [object]$Row) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (Test-Path -LiteralPath $energyPath) {
        $energy = @(Import-Csv -LiteralPath $energyPath)[0]
        $Row.energy_kj_mol = Get-Z4xbDouble $energy.E_rxn_kJmol
        if ([string]$energy.status -ne 'ok') {
            $Row.status = [string]$energy.status
            $Row.failure_class = 'energy_status_gap'
        }
    }
    if (-not (Test-Path -LiteralPath $profilePath)) {
        $Row.failure_class = 'missing_profile'
        return $Row
    }
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    if ($Row.status -eq 'ok' -and [string]$profile.status -ne 'ok') {
        $Row.status = "profile=$($profile.status)"
        $Row.failure_class = 'profile_status_gap'
    }
    $Row.panels = Get-Z4xbInt64 $profile.panels
    $Row.iterations = Get-Z4xbInt64 $profile.iterations
    $Row.final_residual = Get-Z4xbDouble $profile.final_residual
    $Row.true_residual = Get-Z4xbDouble $profile.true_residual_rel
    $Row.sampled_audit_rms = Get-Z4xbDouble $profile.tree_vs_direct_matvec_rel
    $Row.sampled_direct_true_residual = Get-Z4xbDouble (Get-Z4xbProfileValue $profile 'fmm_local_sampled_direct_true_residual' '')
    $Row.matvec_seconds_avg = Get-Z4xbDouble $profile.matvec_seconds_avg
    $Row.effective_taylor_order = Get-Z4xbInt64 (Get-Z4xbProfileValue $profile 'fmm_m2l_expansion_order' $Row.requested_taylor_order)
    $Row.fmm_local_gmres_active = [string](Get-Z4xbProfileValue $profile 'fmm_local_gmres_active' 'false')
    $Row.fmm_local_operator_finite = [string](Get-Z4xbProfileValue $profile 'fmm_local_operator_finite' '')
    $Row.plan_build_count = Get-Z4xbInt64 (Get-Z4xbProfileValue $profile 'fmm_local_plan_build_count' 0)
    $Row.cache_build_count = Get-Z4xbInt64 (Get-Z4xbProfileValue $profile 'fmm_local_cache_build_count' 0)
    $Row.gmres_matvec_calls = Get-Z4xbInt64 (Get-Z4xbProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
    $Row.geometry_guard = [string](Get-Z4xbProfileValue $profile 'fmm_local_geometry_guard' '')
    $Row.geometry_theta = Get-Z4xbDouble (Get-Z4xbProfileValue $profile 'fmm_local_geometry_theta' '')
    if ($Row.status -eq 'ok') {
        if (-not (Test-FmmCliFinite ([double]$Row.sampled_audit_rms)) -or
            -not (Test-FmmCliFinite ([double]$Row.sampled_direct_true_residual))) {
            $Row.failure_class = 'missing_operator_audit'
        } elseif ([double]$Row.true_residual -gt 1.0e-6) {
            $Row.failure_class = 'convergence_gap'
        } elseif ([double]$Row.sampled_audit_rms -gt 8.0e-3 -or
            [double]$Row.sampled_direct_true_residual -gt 8.0e-3) {
            $Row.failure_class = 'operator_accuracy_gap'
        } elseif ([int64]$Row.effective_taylor_order -ne [int64]$Row.requested_taylor_order) {
            $Row.failure_class = 'effective_order_mismatch'
        } elseif ([int64]$Row.plan_build_count -ne 1 -or
            [int64]$Row.cache_build_count -ne 1 -or
            [int64]$Row.gmres_matvec_calls -le 0 -or
            $Row.fmm_local_gmres_active -ne 'true' -or
            $Row.fmm_local_operator_finite -ne 'true') {
            $Row.failure_class = 'lifecycle_gap'
        } else {
            $Row.failure_class = 'pass'
        }
    }
    return $Row
}

function Test-Z4xbRowFinal([object]$Row) {
    return $Row.failure_class -in @('pass','unsupported_current_implementation')
}

function Test-Z4xbShouldRun([string]$RowCsv, [string]$StageName) {
    if ($Force) { return $true }
    if ($RerunStages -contains $StageName) { return $true }
    if (-not (Test-Path -LiteralPath $RowCsv)) { return $true }
    try {
        $row = @(Import-Csv -LiteralPath $RowCsv)[0]
        if ($row.status -eq 'planned') { return $true }
        if ($row.status -ne 'unsupported' -and
            $row.PSObject.Properties.Name -contains 'requested_taylor_order' -and
            $row.PSObject.Properties.Name -contains 'effective_taylor_order' -and
            [int64]$row.requested_taylor_order -ne [int64]$row.effective_taylor_order) {
            return $true
        }
        if (-not $RerunFailed) { return $false }
        return -not (Test-Z4xbRowFinal $row)
    } catch {
        return $true
    }
}

function Invoke-Z4xbCase([object]$Fixture, [object]$Stage) {
    $runName = '{0}_{1}' -f $Fixture.fixture, $Stage.stage
    $runDir = Join-Path $WorkRoot $runName
    $rowCsv = Join-Path $runDir 'fagbi_scope_z4xb_row.csv'
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null

    if (-not (Test-Z4xbShouldRun -RowCsv $rowCsv -StageName $Stage.stage)) {
        Write-Host "resume: $runName"
        return @(Import-Csv -LiteralPath $rowCsv)[0]
    }

    $row = New-Z4xbRow -Fixture $Fixture -Stage $Stage -RunDir $runDir
    if (-not [string]::IsNullOrWhiteSpace($Stage.unsupported_reason)) {
        Write-Host "skip: $runName ($($Stage.unsupported_reason))"
        $row.status = 'unsupported'
        $row.failure_class = 'unsupported_current_implementation'
        $row.exit_code = 0
        $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
        return $row
    }

    Write-Host "run: $runName (timeout ${CaseTimeoutSeconds}s)"
    $argPath = Write-Z4xbArg -Fixture $Fixture -Stage $Stage -RunDir $runDir
    $bare = @(Get-Content -LiteralPath $argPath | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) { throw "Generated arg.in contains bare options: $($bare -join ', ')" }
    if ($PlanOnly) {
        $row.status = 'planned'
        $row.failure_class = 'plan_only'
        $row.wall_seconds = 0.0
        $row.peak_rss_mb = 0.0
        $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
        return $row
    }

    $run = Invoke-FmmCliTitanCase `
        -TitanExe $TitanExe `
        -RunDir $runDir `
        -ArgPath $argPath `
        -CaseTimeoutSeconds $CaseTimeoutSeconds
    $row.exit_code = $run.exit_code
    $row.timed_out = [string]$run.timed_out
    $row.wall_seconds = $run.wall_seconds
    $row.peak_rss_mb = $run.peak_working_set_mb
    if ($run.timed_out) {
        $row.status = 'timeout'
        $row.failure_class = 'runtime_limited'
    } elseif ($run.exit_code -ne 0) {
        $row.status = 'failed'
        $row.failure_class = "exit_code_$($run.exit_code)"
    } else {
        $row.status = 'ok'
        $row.failure_class = 'pass'
        $row = Read-Z4xbOutputs -RunDir $runDir -Row $row
    }
    $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
    return $row
}

function Write-Z4xbOrderSummary([object[]]$Rows, [string]$Path) {
    $summary = New-Object System.Collections.Generic.List[object]
    foreach ($stage in @($Rows.stage | Sort-Object -Unique)) {
        $set = @($Rows | Where-Object { $_.stage -eq $stage })
        $pass = @($set | Where-Object { $_.failure_class -eq 'pass' })
        $unsupported = @($set | Where-Object { $_.failure_class -eq 'unsupported_current_implementation' })
        $summary.Add([pscustomobject]@{
            stage = $stage
            requested_taylor_order = (@($set | Select-Object -First 1).requested_taylor_order)
            rows = $set.Count
            pass_rows = $pass.Count
            unsupported_rows = $unsupported.Count
            max_audit = if ($pass.Count) { ($pass | ForEach-Object { Get-Z4xbDouble $_.sampled_audit_rms } | Measure-Object -Maximum).Maximum } else { [double]::NaN }
            max_sampled_direct_true_residual = if ($pass.Count) { ($pass | ForEach-Object { Get-Z4xbDouble $_.sampled_direct_true_residual } | Measure-Object -Maximum).Maximum } else { [double]::NaN }
            median_matvec_seconds = if ($pass.Count) {
                $vals = @($pass | ForEach-Object { Get-Z4xbDouble $_.matvec_seconds_avg } | Sort-Object)
                $vals[[int][Math]::Floor($vals.Count / 2)]
            } else { [double]::NaN }
            interpretation = if ($unsupported.Count -eq $set.Count) { 'not implemented under guarded gmres-p4 protocol' } else { 'measured guarded gmres-p4 protocol' }
        }) | Out-Null
    }
    $summary.ToArray() | Sort-Object requested_taylor_order | Export-Csv -LiteralPath $Path -NoTypeInformation
}

function Write-Z4xbSvg([object[]]$Rows, [string]$Path) {
    $points = @($Rows | Where-Object { $_.failure_class -eq 'pass' })
    if ($points.Count -eq 0) {
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="420"><text x="24" y="42">No measured passing p-order rows</text></svg>' |
            Set-Content -LiteralPath $Path -Encoding ASCII
        return
    }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="980" height="520" viewBox="0 0 980 520">')
    $svg.Add('<rect width="980" height="520" fill="white"/>')
    $svg.Add('<line x1="88" y1="420" x2="930" y2="420" stroke="#222"/>')
    $svg.Add('<line x1="88" y1="58" x2="88" y2="420" stroke="#222"/>')
    $svg.Add('<text x="22" y="282" font-family="Arial" font-size="14" transform="rotate(-90 22,282)">sampled direct true residual</text>')
    $maxY = [Math]::Max(1.0e-6, ($points | ForEach-Object { Get-Z4xbDouble $_.sampled_direct_true_residual } | Measure-Object -Maximum).Maximum)
    $barW = [Math]::Max(9.0, 760.0 / [Math]::Max(1, $points.Count) - 5.0)
    for ($i = 0; $i -lt $points.Count; ++$i) {
        $row = $points[$i]
        $value = [Math]::Max(0.0, [double]$row.sampled_direct_true_residual)
        $h = [Math]::Max(2.0, 350.0 * ($value / $maxY))
        $x = 98 + $i * ($barW + 5.0)
        $y = 420.0 - $h
        $svg.Add(('<rect x="{0:F1}" y="{1:F1}" width="{2:F1}" height="{3:F1}" fill="#0072B2"/>' -f $x,$y,$barW,$h))
        $label = "$($row.fixture):p$($row.requested_taylor_order)"
        $svg.Add(('<text x="{0:F1}" y="438" font-family="Arial" font-size="9" transform="rotate(50 {0:F1},438)">{1}</text>' -f $x,$label))
    }
    $svg.Add('</svg>')
    $svg | Set-Content -LiteralPath $Path -Encoding ASCII
}

$fixtureRows = Get-Z4xbFixtureRows
$fixtureNames = if ($Quick) {
    @('1sh1')
} elseif ($Fixtures.Count -gt 0) {
    $Fixtures
} else {
    @('barnase','1a63','1sh1','1ubq','2ci2','1vii')
}
$stages = Get-Z4xbStages
if ($RerunStages.Count -gt 0) {
    $knownStages = @($stages | ForEach-Object { $_.stage })
    $unknownStages = @($RerunStages | Where-Object { $knownStages -notcontains $_ })
    if ($unknownStages.Count -gt 0) {
        throw "Unknown -RerunStages value(s): $($unknownStages -join ', '). Known stages: $($knownStages -join ', ')"
    }
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    $fixture = Get-Z4xbFixture -Name $fixtureName -Rows $fixtureRows
    foreach ($stage in $stages) {
        $rows.Add((Invoke-Z4xbCase -Fixture $fixture -Stage $stage)) | Out-Null
    }
}

$allRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    foreach ($stage in $stages) {
        $rowCsv = Join-Path (Join-Path $WorkRoot ('{0}_{1}' -f $fixtureName, $stage.stage)) 'fagbi_scope_z4xb_row.csv'
        if (Test-Path -LiteralPath $rowCsv) {
            $allRows.Add((@(Import-Csv -LiteralPath $rowCsv)[0])) | Out-Null
        }
    }
}
$matrix = @($allRows.ToArray() | Sort-Object fixture,stage_order)
$matrix | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
Write-Z4xbOrderSummary -Rows $matrix -Path $OrderSummaryCsv
Write-Z4xbSvg -Rows $matrix -Path $PlotPath

$measured = @($matrix | Where-Object { $_.status -notin @('unsupported','planned') })
$p4Rows = @($matrix | Where-Object { $_.stage -eq 'p_order4' })
$p4Pass = @($p4Rows | Where-Object { $_.failure_class -eq 'pass' })
$lowerOrderRows = @($matrix | Where-Object { $_.stage -in @('p_order2','p_order3') })
$lowerOrderMeasured = @($lowerOrderRows | Where-Object { $_.status -notin @('unsupported','planned') })
$timeouts = @($matrix | Where-Object { $_.failure_class -eq 'runtime_limited' })
$accuracyGaps = @($measured | Where-Object { $_.failure_class -in @('operator_accuracy_gap','convergence_gap','missing_operator_audit','lifecycle_gap') })
$protocolGaps = @($measured | Where-Object { $_.failure_class -eq 'effective_order_mismatch' })
$unsupported = @($matrix | Where-Object { $_.failure_class -eq 'unsupported_current_implementation' })

$decision = if ($PlanOnly) {
    'z4xb_plan_only'
} elseif ($timeouts.Count -gt 0) {
    'z4xb_runtime_limited'
} elseif ($protocolGaps.Count -gt 0) {
    'z4xb_order_protocol_gap'
} elseif ($p4Rows.Count -gt 0 -and $p4Pass.Count -eq $p4Rows.Count -and
    $lowerOrderRows.Count -gt 0 -and $lowerOrderMeasured.Count -eq $lowerOrderRows.Count) {
    if ($accuracyGaps.Count -gt 0) { 'z4xb_p4_confirmed_lower_order_gaps' } else { 'z4xb_lower_orders_also_pass' }
} elseif ($accuracyGaps.Count -gt 0) {
    'z4xb_order_accuracy_gap'
} elseif ($p4Rows.Count -gt 0 -and $p4Pass.Count -eq $p4Rows.Count -and $unsupported.Count -gt 0) {
    'z4xb_p4_only_current_implementation'
} elseif ($p4Rows.Count -gt 0 -and $p4Pass.Count -eq $p4Rows.Count) {
    'z4xb_p4_supported_policy'
} else {
    'z4xb_harness_gap'
}

[pscustomobject]@{
    decision = $decision
    quick = [bool]$Quick
    fixtures = ($fixtureNames -join ';')
    rows = $matrix.Count
    measured_rows = $measured.Count
    lower_order_rows = $lowerOrderRows.Count
    lower_order_measured_rows = $lowerOrderMeasured.Count
    p4_rows = $p4Rows.Count
    p4_pass_rows = $p4Pass.Count
    unsupported_rows = $unsupported.Count
    timeout_rows = $timeouts.Count
    accuracy_gap_rows = $accuracyGaps.Count
    protocol_gap_rows = $protocolGaps.Count
    case_timeout_seconds = $CaseTimeoutSeconds
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
    'scope: Z4X-B common-protocol p-order sweep',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "fixtures: $($fixtureNames -join ', ')",
    'protocol: full protein, guarded gmres-p4, theta 0.6, strict geometry theta 0.5, operator-tail-static tol 0.002, green-jet-compressed, packed-avx2, ILU0 cap 64',
    'orders: p2 measured, p3 measured, p4 measured, p5 measured when rerun under the p5-enabled binary',
    'production_default: grid',
    'tabi_treecode_default: degree-one',
    'gmres_p4: explicit advanced path only',
    "matrix: $MatrixCsv",
    "summary: $OrderSummaryCsv"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII

$report = New-Object System.Collections.Generic.List[string]
$report.Add('# TITAN FMM/M2L Scope Z4X-B p-Order Protocol Sweep Report')
$report.Add('')
$report.Add("Decision: $decision")
$report.Add('')
$report.Add('Purpose')
$report.Add('Run a common-protocol Taylor-order sweep for the guarded local advanced TABI path. Orders p=2,3,4,5 are measured under the same static refinement, compressed Green-jet cache, packed AVX2 fallback, and ILU(0) protocol when supported by the current binary.')
$report.Add('')
$report.Add('Protocol')
$report.Add('- Full public protein PQR fixtures; no atom sampling.')
$report.Add('- Analytical FMM/M2L, gmres-p4, static tail tolerance 0.002, strict geometry theta 0.5.')
$report.Add('- Compressed Green-jet cache, packed AVX2 fallback, sparse near-field ILU(0) neighbor cap 64.')
$report.Add('- Grid remains production default; TABI treecode default remains degree-one.')
$report.Add('')
$report.Add('Rows')
foreach ($row in $matrix) {
    $report.Add(('- {0} {1}: status={2}, failure={3}, requested_p={4}, residual={5}, audit={6}, direct_residual={7}, wall={8}s, note={9}' -f
            $row.fixture,
            $row.stage,
            $row.status,
            $row.failure_class,
            $row.requested_taylor_order,
            $row.true_residual,
            $row.sampled_audit_rms,
            $row.sampled_direct_true_residual,
            $row.wall_seconds,
            $row.unsupported_reason))
}
$report.Add('')
$report.Add('Artifacts')
$report.Add("Matrix: $MatrixCsv")
$report.Add("Order summary: $OrderSummaryCsv")
$report.Add("Decision: $DecisionCsv")
$report.Add("Plot: $PlotPath")
$report | Set-Content -LiteralPath $RootReport -Encoding ASCII

Write-Host "Z4X-B decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Resume: rerun the same command to continue missing/mismatched rows; use -RerunFailed for failed/timeout rows, -RerunStages <stage>, or -Force for all rows."

if ($decision -notin @('z4xb_p4_supported_policy','z4xb_p4_only_current_implementation','z4xb_p4_confirmed_lower_order_gaps','z4xb_lower_orders_also_pass','z4xb_plan_only') -and
    -not $AllowFailedDecision) {
    throw "Z4X-B p-order protocol sweep did not pass: $decision"
}
