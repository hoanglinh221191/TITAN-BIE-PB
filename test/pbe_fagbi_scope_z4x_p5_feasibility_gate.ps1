param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$PrimaryFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$ExpansionFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public_z4q'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4xc'),
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

$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4xc_matrix.csv'
$FallbackCsv = Join-Path $WorkRoot 'fagbi_scope_z4xc_fallback_work.csv'
$CacheCsv = Join-Path $WorkRoot 'fagbi_scope_z4xc_cache.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4xc_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4xc_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4xc_p5_feasibility.svg'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4XC_P5_FEASIBILITY_REPORT.txt'

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

function Get-Z4xcDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $parsed = 0.0
    if ([double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-Z4xcInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $parsed = [int64]0
    if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-Z4xcProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Import-Z4xcManifestRows([string]$Root, [string]$ManifestName) {
    $path = Join-Path $Root $ManifestName
    if (-not (Test-Path -LiteralPath $path)) { return @() }
    return @(Import-Csv -LiteralPath $path | ForEach-Object {
        $_ | Add-Member -NotePropertyName fixture_root -NotePropertyValue $Root -Force
        $_
    })
}

function Get-Z4xcFixtureRows {
    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($row in (Import-Z4xcManifestRows $PrimaryFixtureRoot 'public_fixture_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    foreach ($row in (Import-Z4xcManifestRows $ExpansionFixtureRoot 'public_fixture_expansion_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    return $rows.ToArray()
}

function Get-Z4xcFixture([string]$Name, [object[]]$Rows) {
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

function New-Z4xcStage([int]$Order, [string]$Stage, [int]$TaylorOrder, [string]$LocalMode) {
    return [pscustomobject]@{
        order = $Order
        stage = $Stage
        taylor_order = $TaylorOrder
        local_mode = $LocalMode
        family = if ($LocalMode -eq 'gmres-p4') { 'gmres' } else { 'single_matvec' }
    }
}

function Get-Z4xcStages {
    return @(
        (New-Z4xcStage 1 'single_p4' 4 'single-matvec-p4'),
        (New-Z4xcStage 2 'single_p5' 5 'single-matvec-p4'),
        (New-Z4xcStage 3 'gmres_p4' 4 'gmres-p4'),
        (New-Z4xcStage 4 'gmres_p5' 5 'gmres-p4')
    )
}

function New-Z4xcRow([object]$Fixture, [object]$Stage, [string]$RunDir) {
    return [pscustomobject]@{
        fixture = $Fixture.fixture
        pdb_id = $Fixture.pdb_id
        protein_name = $Fixture.protein_name
        atoms = Get-Z4xcInt64 $Fixture.atom_count
        stage = $Stage.stage
        stage_order = $Stage.order
        family = $Stage.family
        fmm_local_mode = $Stage.local_mode
        requested_taylor_order = $Stage.taylor_order
        effective_taylor_order = 0
        status = 'planned'
        failure_class = 'planned'
        exit_code = 0
        timed_out = 'False'
        panels = 0
        iterations = 0
        energy_kj_mol = [double]::NaN
        final_residual = [double]::NaN
        true_residual = [double]::NaN
        sampled_audit_rms = [double]::NaN
        sampled_direct_true_residual = [double]::NaN
        local_audit_rms = [double]::NaN
        reference_gap = [double]::NaN
        matvec_seconds_avg = [double]::NaN
        candidate_median_apply_seconds = [double]::NaN
        wall_seconds = [double]::NaN
        peak_rss_mb = [double]::NaN
        fmm_local_operator_finite = ''
        plan_build_count = 0
        cache_build_count = 0
        gmres_matvec_calls = 0
        fallback_pairs = 0
        far_pairs = 0
        cache_bytes = 0
        cache_operator_bytes = 0
        cache_layout = ''
        fallback_speed_seconds = [double]::NaN
        run_dir = $RunDir
    }
}

function Write-Z4xcArg([object]$Fixture, [object]$Stage, [string]$RunDir) {
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
        -FmmLocal $Stage.local_mode `
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

function Read-Z4xcAttribution([string]$RunDir, [object]$Row) {
    $summaryPath = Join-Path $RunDir 'pbe_tabi_operator_attribution_summary.csv'
    if (-not (Test-Path -LiteralPath $summaryPath)) { return $Row }
    $summary = @(Import-Csv -LiteralPath $summaryPath)
    $local = @($summary | Where-Object {
        $_.approximation_mode -like 'local_m2l*green_cached' -and
        $_.probe -eq 'solution' -and
        $_.row_component -eq 'both' -and
        $_.region -eq 'total' -and
        $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
    if ($null -eq $local) { return $Row }
    $Row.local_audit_rms = Get-Z4xcDouble $local.aggregate_audit_rel
    $Row.reference_gap = Get-Z4xcDouble $local.reference_gap 0.0
    $Row.effective_taylor_order = Get-Z4xcInt64 $local.fmm_local_source_order $Row.effective_taylor_order
    $Row.fmm_local_operator_finite = [string]$local.fmm_local_finite_outputs
    $Row.plan_build_count = Get-Z4xcInt64 $local.fmm_local_plan_build_count $Row.plan_build_count
    $Row.cache_build_count = Get-Z4xcInt64 $local.fmm_local_cache_build_count $Row.cache_build_count
    $Row.fallback_pairs = Get-Z4xcInt64 $local.fmm_local_direct_fallback_pairs
    $Row.far_pairs = Get-Z4xcInt64 $local.fmm_local_far_pairs_after
    $Row.cache_bytes = Get-Z4xcInt64 $local.fmm_local_cache_bytes
    $Row.cache_operator_bytes = Get-Z4xcInt64 $local.fmm_local_cache_operator_bytes
    $Row.cache_layout = [string]$local.fmm_local_cache_layout
    $Row.candidate_median_apply_seconds = Get-Z4xcDouble $local.fmm_local_candidate_median_apply_seconds
    $Row.fallback_speed_seconds = Get-Z4xcDouble $local.fmm_local_fallback_seconds
    return $Row
}

function Read-Z4xcOutputs([string]$RunDir, [object]$Row) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (Test-Path -LiteralPath $energyPath) {
        $energy = @(Import-Csv -LiteralPath $energyPath)[0]
        $Row.energy_kj_mol = Get-Z4xcDouble $energy.E_rxn_kJmol
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
    $Row.panels = Get-Z4xcInt64 $profile.panels
    $Row.iterations = Get-Z4xcInt64 $profile.iterations
    $Row.final_residual = Get-Z4xcDouble $profile.final_residual
    $Row.true_residual = Get-Z4xcDouble $profile.true_residual_rel
    $Row.sampled_audit_rms = Get-Z4xcDouble $profile.tree_vs_direct_matvec_rel
    $Row.sampled_direct_true_residual = Get-Z4xcDouble (Get-Z4xcProfileValue $profile 'fmm_local_sampled_direct_true_residual' '')
    $Row.matvec_seconds_avg = Get-Z4xcDouble $profile.matvec_seconds_avg
    $Row.effective_taylor_order = Get-Z4xcInt64 (Get-Z4xcProfileValue $profile 'fmm_m2l_expansion_order' $Row.requested_taylor_order)
    $Row.fmm_local_operator_finite = [string](Get-Z4xcProfileValue $profile 'fmm_local_operator_finite' $Row.fmm_local_operator_finite)
    $Row.plan_build_count = Get-Z4xcInt64 (Get-Z4xcProfileValue $profile 'fmm_local_plan_build_count' $Row.plan_build_count)
    $Row.cache_build_count = Get-Z4xcInt64 (Get-Z4xcProfileValue $profile 'fmm_local_cache_build_count' $Row.cache_build_count)
    $Row.gmres_matvec_calls = Get-Z4xcInt64 (Get-Z4xcProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
    $Row = Read-Z4xcAttribution -RunDir $RunDir -Row $Row
    if ($Row.status -eq 'ok') {
        if ([int64]$Row.effective_taylor_order -ne [int64]$Row.requested_taylor_order) {
            $Row.failure_class = 'effective_order_mismatch'
        } elseif ($Row.family -eq 'gmres' -and ([double]$Row.true_residual -gt 1.0e-6)) {
            $Row.failure_class = 'convergence_gap'
        } elseif ($Row.family -eq 'gmres' -and
            ((-not (Test-FmmCliFinite ([double]$Row.sampled_audit_rms))) -or
             (-not (Test-FmmCliFinite ([double]$Row.sampled_direct_true_residual))) -or
             [double]$Row.sampled_audit_rms -gt 8.0e-3 -or
             [double]$Row.sampled_direct_true_residual -gt 8.0e-3)) {
            $Row.failure_class = 'operator_accuracy_gap'
        } elseif ($Row.family -eq 'single_matvec' -and
            ((-not (Test-FmmCliFinite ([double]$Row.local_audit_rms))) -or
             [double]$Row.local_audit_rms -gt 8.0e-3 -or
             [double]$Row.reference_gap -gt 1.0e-10)) {
            $Row.failure_class = 'single_matvec_accuracy_gap'
        } elseif ([int64]$Row.plan_build_count -ne 1 -or
            [int64]$Row.cache_build_count -ne 1 -or
            $Row.fmm_local_operator_finite -ne 'true') {
            $Row.failure_class = 'lifecycle_gap'
        } elseif ([double]$Row.peak_rss_mb -ge 512.0) {
            $Row.failure_class = 'memory_limited'
        } else {
            $Row.failure_class = 'pass'
        }
    }
    return $Row
}

function Test-Z4xcRowFinal([object]$Row) {
    return $Row.failure_class -eq 'pass'
}

function Update-Z4xcFailureFromLogs([string]$RunDir, [object]$Row) {
    $stderrPath = Join-Path $RunDir 'stderr.log'
    if (Test-Path -LiteralPath $stderrPath) {
        $stderr = Get-Content -LiteralPath $stderrPath -Raw
        if ($stderr -match 'cache_budget_exceeded') {
            $Row.status = 'failed'
            $Row.failure_class = 'memory_limited'
            $required = [regex]::Match($stderr, 'required=([0-9]+)')
            $budget = [regex]::Match($stderr, 'budget=([0-9]+)')
            if ($required.Success) { $Row.cache_bytes = [int64]$required.Groups[1].Value }
            if ($budget.Success -and [int64]$Row.cache_operator_bytes -eq 0) {
                $Row.cache_layout = 'green-jet-compressed'
            }
            if ([int64]$Row.effective_taylor_order -eq 0) {
                $Row.effective_taylor_order = $Row.requested_taylor_order
            }
        }
    }
    return $Row
}

function Test-Z4xcShouldRun([string]$RowCsv, [string]$StageName) {
    if ($Force) { return $true }
    if ($RerunStages -contains $StageName) { return $true }
    if (-not (Test-Path -LiteralPath $RowCsv)) { return $true }
    try {
        $row = @(Import-Csv -LiteralPath $RowCsv)[0]
        $row = Update-Z4xcFailureFromLogs -RunDir (Split-Path -Parent $RowCsv) -Row $row
        $row | Export-Csv -LiteralPath $RowCsv -NoTypeInformation
        if ($row.status -eq 'planned') { return $true }
        if ($row.failure_class -eq 'memory_limited' -and -not $RerunFailed) {
            return $false
        }
        if ($row.PSObject.Properties.Name -contains 'requested_taylor_order' -and
            $row.PSObject.Properties.Name -contains 'effective_taylor_order' -and
            [int64]$row.requested_taylor_order -ne [int64]$row.effective_taylor_order) {
            return $true
        }
        if (-not $RerunFailed) { return $false }
        return -not (Test-Z4xcRowFinal $row)
    } catch {
        return $true
    }
}

function Invoke-Z4xcCase([object]$Fixture, [object]$Stage) {
    $runName = '{0}_{1}' -f $Fixture.fixture, $Stage.stage
    $runDir = Join-Path $WorkRoot $runName
    $rowCsv = Join-Path $runDir 'fagbi_scope_z4xc_row.csv'
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    if (-not (Test-Z4xcShouldRun -RowCsv $rowCsv -StageName $Stage.stage)) {
        Write-Host "resume: $runName"
        return @(Import-Csv -LiteralPath $rowCsv)[0]
    }
    Write-Host "run: $runName (timeout ${CaseTimeoutSeconds}s)"
    $row = New-Z4xcRow -Fixture $Fixture -Stage $Stage -RunDir $runDir
    $argPath = Write-Z4xcArg -Fixture $Fixture -Stage $Stage -RunDir $runDir
    $bare = @(Get-Content -LiteralPath $argPath | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) { throw "Generated arg.in contains bare options: $($bare -join ', ')" }
    if ($PlanOnly) {
        $row.status = 'planned'
        $row.failure_class = 'plan_only'
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
        $row = Update-Z4xcFailureFromLogs -RunDir $runDir -Row $row
    } else {
        $row.status = 'ok'
        $row.failure_class = 'pass'
        $row = Read-Z4xcOutputs -RunDir $runDir -Row $row
    }
    $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
    return $row
}

function Write-Z4xcSvg([object[]]$Rows, [string]$Path) {
    $gmres = @($Rows | Where-Object { $_.family -eq 'gmres' -and $_.failure_class -eq 'pass' })
    if ($gmres.Count -eq 0) {
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="420"><text x="24" y="42">No passing GMRES rows</text></svg>' |
            Set-Content -LiteralPath $Path -Encoding ASCII
        return
    }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="980" height="520" viewBox="0 0 980 520">')
    $svg.Add('<rect width="980" height="520" fill="white"/>')
    $svg.Add('<line x1="90" y1="420" x2="930" y2="420" stroke="#222"/>')
    $svg.Add('<line x1="90" y1="62" x2="90" y2="420" stroke="#222"/>')
    $svg.Add('<text x="24" y="290" font-family="Arial" font-size="14" transform="rotate(-90 24,290)">matvec seconds</text>')
    $maxY = [Math]::Max(1.0e-6, ($gmres | ForEach-Object { Get-Z4xcDouble $_.matvec_seconds_avg } | Measure-Object -Maximum).Maximum)
    $barW = [Math]::Max(10.0, 760.0 / [Math]::Max(1, $gmres.Count) - 5.0)
    for ($i = 0; $i -lt $gmres.Count; ++$i) {
        $row = $gmres[$i]
        $value = [Math]::Max(0.0, [double]$row.matvec_seconds_avg)
        $h = [Math]::Max(2.0, 350.0 * ($value / $maxY))
        $x = 100 + $i * ($barW + 5.0)
        $y = 420.0 - $h
        $fill = if ([int]$row.requested_taylor_order -eq 5) { '#D55E00' } else { '#0072B2' }
        $svg.Add(('<rect x="{0:F1}" y="{1:F1}" width="{2:F1}" height="{3:F1}" fill="{4}"/>' -f $x,$y,$barW,$h,$fill))
        $label = "$($row.fixture):p$($row.requested_taylor_order)"
        $svg.Add(('<text x="{0:F1}" y="438" font-family="Arial" font-size="9" transform="rotate(50 {0:F1},438)">{1}</text>' -f $x,$label))
    }
    $svg.Add('</svg>')
    $svg | Set-Content -LiteralPath $Path -Encoding ASCII
}

$fixtureRows = Get-Z4xcFixtureRows
$fixtureNames = if ($Quick) {
    @('1sh1','1a63')
} elseif ($Fixtures.Count -gt 0) {
    $Fixtures
} else {
    @('barnase','1a63','1sh1','1ubq','2ci2','1vii')
}
$stages = Get-Z4xcStages
if ($RerunStages.Count -gt 0) {
    $knownStages = @($stages | ForEach-Object { $_.stage })
    $unknownStages = @($RerunStages | Where-Object { $knownStages -notcontains $_ })
    if ($unknownStages.Count -gt 0) {
        throw "Unknown -RerunStages value(s): $($unknownStages -join ', '). Known stages: $($knownStages -join ', ')"
    }
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    $fixture = Get-Z4xcFixture -Name $fixtureName -Rows $fixtureRows
    foreach ($stage in $stages) {
        $rows.Add((Invoke-Z4xcCase -Fixture $fixture -Stage $stage)) | Out-Null
    }
}

$allRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    foreach ($stage in $stages) {
        $rowCsv = Join-Path (Join-Path $WorkRoot ('{0}_{1}' -f $fixtureName, $stage.stage)) 'fagbi_scope_z4xc_row.csv'
        if (Test-Path -LiteralPath $rowCsv) {
            $allRows.Add((@(Import-Csv -LiteralPath $rowCsv)[0])) | Out-Null
        }
    }
}
$matrix = @($allRows.ToArray() | Sort-Object fixture,stage_order)
$matrix | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation

$fallbackRows = @()
$cacheRows = @()
foreach ($fixtureName in $fixtureNames) {
    $p4 = @($matrix | Where-Object { $_.fixture -eq $fixtureName -and $_.stage -eq 'gmres_p4' }) | Select-Object -First 1
    $p5 = @($matrix | Where-Object { $_.fixture -eq $fixtureName -and $_.stage -eq 'gmres_p5' }) | Select-Object -First 1
    if ($p4 -and $p5) {
        $fallbackRatio = if ([double]$p4.fallback_pairs -gt 0) { [double]$p5.fallback_pairs / [double]$p4.fallback_pairs } else { [double]::NaN }
        $matvecRatio = if ([double]$p4.matvec_seconds_avg -gt 0) { [double]$p5.matvec_seconds_avg / [double]$p4.matvec_seconds_avg } else { [double]::NaN }
        $fallbackRows += [pscustomobject]@{
            fixture = $fixtureName
            p4_fallback_pairs = $p4.fallback_pairs
            p5_fallback_pairs = $p5.fallback_pairs
            p5_vs_p4_fallback_ratio = $fallbackRatio
            p4_matvec_seconds = $p4.matvec_seconds_avg
            p5_matvec_seconds = $p5.matvec_seconds_avg
            p5_vs_p4_matvec_ratio = $matvecRatio
        }
    }
    foreach ($row in @($p4,$p5)) {
        if ($null -ne $row) {
            $cacheRows += [pscustomobject]@{
                fixture = $fixtureName
                stage = $row.stage
                requested_taylor_order = $row.requested_taylor_order
                cache_layout = $row.cache_layout
                cache_operator_bytes = $row.cache_operator_bytes
                cache_bytes = $row.cache_bytes
                peak_rss_mb = $row.peak_rss_mb
            }
        }
    }
}
$fallbackRows | Export-Csv -LiteralPath $FallbackCsv -NoTypeInformation
$cacheRows | Export-Csv -LiteralPath $CacheCsv -NoTypeInformation
Write-Z4xcSvg -Rows $matrix -Path $PlotPath

$gmresP4 = @($matrix | Where-Object { $_.stage -eq 'gmres_p4' })
$gmresP5 = @($matrix | Where-Object { $_.stage -eq 'gmres_p5' })
$requiredP5 = @($gmresP5 | Where-Object { $_.status -ne 'planned' })
$timeouts = @($matrix | Where-Object { $_.failure_class -eq 'runtime_limited' })
$memory = @($matrix | Where-Object { $_.failure_class -eq 'memory_limited' })
$accuracy = @($matrix | Where-Object { $_.failure_class -in @('operator_accuracy_gap','single_matvec_accuracy_gap','convergence_gap','effective_order_mismatch','lifecycle_gap') })
$p5Pass = @($requiredP5 | Where-Object { $_.failure_class -eq 'pass' })
$fallbackGood = @($fallbackRows | Where-Object { (Get-Z4xcDouble $_.p5_vs_p4_fallback_ratio 1.0) -le 0.75 })
$matvecGood = @($fallbackRows | Where-Object { (Get-Z4xcDouble $_.p5_vs_p4_matvec_ratio 1.0) -lt 1.0 })

$decision = if ($PlanOnly) {
    'z4xc_plan_only'
} elseif ($timeouts.Count -gt 0) {
    'z4xc_runtime_limited'
} elseif ($memory.Count -gt 0) {
    'z4xc_p5_memory_limited'
} elseif ($accuracy.Count -gt 0) {
    'z4xc_p5_accuracy_gap'
} elseif ($p5Pass.Count -ne $gmresP5.Count) {
    'z4xc_harness_gap'
} elseif ($fallbackGood.Count -lt $fallbackRows.Count) {
    'z4xc_p5_no_fallback_gain'
} elseif ($matvecGood.Count -eq $fallbackRows.Count) {
    'z4xc_p5_speed_pass'
} else {
    'z4xc_p5_accuracy_pass_runtime_gap'
}

[pscustomobject]@{
    decision = $decision
    quick = [bool]$Quick
    fixtures = ($fixtureNames -join ';')
    rows = $matrix.Count
    gmres_p4_rows = $gmresP4.Count
    gmres_p5_rows = $gmresP5.Count
    gmres_p5_pass_rows = $p5Pass.Count
    timeout_rows = $timeouts.Count
    memory_rows = $memory.Count
    accuracy_gap_rows = $accuracy.Count
    fallback_good_rows = $fallbackGood.Count
    matvec_good_rows = $matvecGood.Count
    case_timeout_seconds = $CaseTimeoutSeconds
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
    'scope: Z4X-C p5 guarded local FMM/M2L feasibility',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "fixtures: $($fixtureNames -join ', ')",
    'protocol: full protein, theta 0.6, strict geometry theta 0.5, operator-tail-static tol 0.002, green-jet-compressed, packed-avx2, ILU0 cap 64',
    'production_default: grid',
    'tabi_treecode_default: degree-one',
    'gmres_p4_name: retained for compatibility; effective order is reported separately',
    "matrix: $MatrixCsv",
    "fallback: $FallbackCsv",
    "cache: $CacheCsv"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII

$report = New-Object System.Collections.Generic.List[string]
$report.Add('# TITAN FMM/M2L Scope Z4X-C p5 Feasibility Report')
$report.Add('')
$report.Add("Decision: $decision")
$report.Add('')
$report.Add('Rows')
foreach ($row in $matrix) {
    $report.Add(('- {0} {1}: status={2}, failure={3}, p={4}/{5}, residual={6}, audit={7}, direct_residual={8}, local_audit={9}, fallback_pairs={10}, matvec={11}s, wall={12}s, rss={13} MB' -f
            $row.fixture,
            $row.stage,
            $row.status,
            $row.failure_class,
            $row.requested_taylor_order,
            $row.effective_taylor_order,
            $row.true_residual,
            $row.sampled_audit_rms,
            $row.sampled_direct_true_residual,
            $row.local_audit_rms,
            $row.fallback_pairs,
            $row.matvec_seconds_avg,
            $row.wall_seconds,
            $row.peak_rss_mb))
}
$report.Add('')
$report.Add('Artifacts')
$report.Add("Matrix: $MatrixCsv")
$report.Add("Fallback summary: $FallbackCsv")
$report.Add("Cache summary: $CacheCsv")
$report.Add("Decision: $DecisionCsv")
$report.Add("Plot: $PlotPath")
$report | Set-Content -LiteralPath $RootReport -Encoding ASCII

Write-Host "Z4X-C decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Resume: rerun the same command to continue missing/mismatched rows; use -RerunFailed for failed/timeout rows, -RerunStages <stage>, or -Force for all rows."

if ($decision -notin @('z4xc_p5_speed_pass','z4xc_p5_accuracy_pass_runtime_gap','z4xc_p5_no_fallback_gain','z4xc_plan_only') -and
    -not $AllowFailedDecision) {
    throw "Z4X-C p5 feasibility gate did not pass: $decision"
}
