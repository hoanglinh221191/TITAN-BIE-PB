param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$PrimaryFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$ExpansionFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public_z4q'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4xf'),
    [string[]]$Fixtures = @(),
    [int[]]$IluCaps = @(32,64,128),
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

$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4xf_matrix.csv'
$CapSummaryCsv = Join-Path $WorkRoot 'fagbi_scope_z4xf_neighbor_cap_summary.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4xf_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4xf_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4xf_neighbor_cap_plot.svg'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4XF_ILU_NEIGHBOR_CAP_SWEEP_REPORT.txt'

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

function Get-Z4xfDouble([object]$Value, [double]$Default = [double]::NaN) {
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

function Get-Z4xfInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $parsed = [int64]0
    if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Test-Z4xfFinite([object]$Value) {
    $x = Get-Z4xfDouble $Value
    return -not ([double]::IsNaN($x) -or [double]::IsInfinity($x))
}

function Get-Z4xfProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Import-Z4xfManifestRows([string]$Root, [string]$ManifestName) {
    $path = Join-Path $Root $ManifestName
    if (-not (Test-Path -LiteralPath $path)) { return @() }
    return @(Import-Csv -LiteralPath $path | ForEach-Object {
        $_ | Add-Member -NotePropertyName fixture_root -NotePropertyValue $Root -Force
        $_
    })
}

function Get-Z4xfFixtureRows {
    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($row in (Import-Z4xfManifestRows $PrimaryFixtureRoot 'public_fixture_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    foreach ($row in (Import-Z4xfManifestRows $ExpansionFixtureRoot 'public_fixture_expansion_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    return $rows.ToArray()
}

function Get-Z4xfFixture([string]$Name, [object[]]$Rows) {
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

function New-Z4xfStage([int]$Cap) {
    return [pscustomobject]@{
        stage = "ilu_cap$Cap"
        stage_order = $Cap
        ilu_neighbor_cap = $Cap
    }
}

function New-Z4xfRow([object]$Fixture, [object]$Stage, [string]$RunDir) {
    return [pscustomobject]@{
        fixture = $Fixture.fixture
        pdb_id = $Fixture.pdb_id
        protein_name = $Fixture.protein_name
        atoms = Get-Z4xfInt64 $Fixture.atom_count
        stage = $Stage.stage
        ilu_neighbor_cap = $Stage.ilu_neighbor_cap
        requested_taylor_order = 4
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
        matvec_seconds_avg = [double]::NaN
        wall_seconds = [double]::NaN
        peak_rss_mb = [double]::NaN
        preconditioner_effective = ''
        preconditioner_setup_seconds = [double]::NaN
        preconditioner_apply_seconds = [double]::NaN
        preconditioner_sparse_nnz_blocks = 0
        preconditioner_sparse_avg_degree = [double]::NaN
        preconditioner_sparse_max_degree = 0
        fmm_local_gmres_active = ''
        fmm_local_operator_finite = ''
        plan_build_count = 0
        cache_build_count = 0
        gmres_matvec_calls = 0
        fmm_local_prepare_seconds = [double]::NaN
        fmm_local_gmres_apply_seconds_mean = [double]::NaN
        fmm_local_fallback_pairs = 0
        fmm_local_far_pairs_after = 0
        run_dir = $RunDir
    }
}

function Write-Z4xfArg([object]$Fixture, [object]$Stage, [string]$RunDir) {
    return Write-FmmCliArgFile -RunDir $RunDir -RepoRoot $RuntimeRepoRoot `
        -PqrPath $Fixture.pqr_path -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -SurfacePoints 8 -MaxIter 720 -Restart 128 `
        -Tolerance 1e-6 -Threads $Threads -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount 512 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerBlockPanelCap 32 `
        -PreconditionerOverlapPanelCap 64 -PreconditionerIluNeighborCap $Stage.ilu_neighbor_cap `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder 4 -FmmCartesianAudit off -FmmLocal gmres-p4 `
        -FmmLocalRefine operator-tail-static -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache green-jet-compressed -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback packed-avx2 -FmmLocalPlanBuilder metric-cache `
        -FmmLocalRefineProbeCount 4 -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard strict-geometry-mac -FmmLocalGeometryTheta 0.5 `
        -FmmLocalIluCapSweep $true
}

function Read-Z4xfOutputs([string]$RunDir, [object]$Row) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (Test-Path -LiteralPath $energyPath) {
        $energy = @(Import-Csv -LiteralPath $energyPath)[0]
        $Row.energy_kj_mol = Get-Z4xfDouble $energy.E_rxn_kJmol
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
    $Row.panels = Get-Z4xfInt64 $profile.panels
    $Row.iterations = Get-Z4xfInt64 $profile.iterations
    $Row.final_residual = Get-Z4xfDouble $profile.final_residual
    $Row.true_residual = Get-Z4xfDouble $profile.true_residual_rel
    $Row.sampled_audit_rms = Get-Z4xfDouble $profile.tree_vs_direct_matvec_rel
    $Row.sampled_direct_true_residual =
        Get-Z4xfDouble (Get-Z4xfProfileValue $profile 'fmm_local_sampled_direct_true_residual' '')
    $Row.matvec_seconds_avg = Get-Z4xfDouble $profile.matvec_seconds_avg
    $Row.effective_taylor_order =
        Get-Z4xfInt64 (Get-Z4xfProfileValue $profile 'fmm_m2l_expansion_order' 0)
    $Row.preconditioner_effective =
        [string](Get-Z4xfProfileValue $profile 'preconditioner_effective' '')
    $Row.preconditioner_setup_seconds =
        Get-Z4xfDouble (Get-Z4xfProfileValue $profile 'preconditioner_setup_seconds' '')
    $Row.preconditioner_apply_seconds =
        Get-Z4xfDouble (Get-Z4xfProfileValue $profile 'preconditioner_apply_seconds' '')
    $Row.preconditioner_sparse_nnz_blocks =
        Get-Z4xfInt64 (Get-Z4xfProfileValue $profile 'preconditioner_sparse_nnz_blocks' 0)
    $Row.preconditioner_sparse_avg_degree =
        Get-Z4xfDouble (Get-Z4xfProfileValue $profile 'preconditioner_sparse_avg_degree' '')
    $Row.preconditioner_sparse_max_degree =
        Get-Z4xfInt64 (Get-Z4xfProfileValue $profile 'preconditioner_sparse_max_degree' 0)
    $Row.fmm_local_gmres_active =
        [string](Get-Z4xfProfileValue $profile 'fmm_local_gmres_active' '')
    $Row.fmm_local_operator_finite =
        [string](Get-Z4xfProfileValue $profile 'fmm_local_operator_finite' '')
    $Row.plan_build_count =
        Get-Z4xfInt64 (Get-Z4xfProfileValue $profile 'fmm_local_plan_build_count' 0)
    $Row.cache_build_count =
        Get-Z4xfInt64 (Get-Z4xfProfileValue $profile 'fmm_local_cache_build_count' 0)
    $Row.gmres_matvec_calls =
        Get-Z4xfInt64 (Get-Z4xfProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
    $Row.fmm_local_prepare_seconds =
        Get-Z4xfDouble (Get-Z4xfProfileValue $profile 'fmm_local_prepare_seconds' '')
    $Row.fmm_local_gmres_apply_seconds_mean =
        Get-Z4xfDouble (Get-Z4xfProfileValue $profile 'fmm_local_gmres_apply_seconds_mean' '')
    $Row.fmm_local_fallback_pairs =
        Get-Z4xfInt64 (Get-Z4xfProfileValue $profile 'fmm_local_direct_fallback_pairs' 0)
    $Row.fmm_local_far_pairs_after =
        Get-Z4xfInt64 (Get-Z4xfProfileValue $profile 'fmm_local_far_pairs_after' 0)

    $accuracyOk =
        (Test-Z4xfFinite $Row.energy_kj_mol) -and
        (Get-Z4xfDouble $Row.true_residual 1.0) -le 1.0e-6 -and
        (Get-Z4xfDouble $Row.sampled_audit_rms 1.0) -le 8.0e-3 -and
        (Get-Z4xfDouble $Row.sampled_direct_true_residual 1.0) -le 8.0e-3
    $lifecycleOk =
        $Row.fmm_local_gmres_active -eq 'true' -and
        $Row.fmm_local_operator_finite -eq 'true' -and
        [int64]$Row.plan_build_count -eq 1 -and
        [int64]$Row.cache_build_count -eq 1
    if ($Row.failure_class -eq 'pass' -and -not $accuracyOk) {
        $Row.failure_class = 'operator_accuracy_gap'
    } elseif ($Row.failure_class -eq 'pass' -and -not $lifecycleOk) {
        $Row.failure_class = 'cache_lifecycle_gap'
    } elseif ($Row.failure_class -eq 'pass' -and (Get-Z4xfDouble $Row.peak_rss_mb 999999.0) -ge 512.0) {
        $Row.failure_class = 'memory_limited'
    }
    return $Row
}

function Test-Z4xfShouldRun([string]$RowCsv, [string]$Stage) {
    if ($Force) { return $true }
    if ($RerunStages.Count -gt 0 -and
        ($RerunStages -contains $Stage -or $RerunStages -contains ($Stage -replace '^ilu_cap',''))) {
        return $true
    }
    if (-not (Test-Path -LiteralPath $RowCsv)) { return $true }
    try {
        $old = @(Import-Csv -LiteralPath $RowCsv)[0]
        if ($old.status -eq 'planned' -or $old.failure_class -eq 'plan_only') {
            return $true
        }
    } catch {
        return $true
    }
    if (-not $RerunFailed) { return $false }
    try {
        $old = @(Import-Csv -LiteralPath $RowCsv)[0]
        return ($old.failure_class -ne 'pass' -and $old.status -ne 'planned')
    } catch {
        return $true
    }
}

function Invoke-Z4xfCase([object]$Fixture, [object]$Stage) {
    $runName = '{0}_{1}' -f $Fixture.fixture, $Stage.stage
    $runDir = Join-Path $WorkRoot $runName
    $rowCsv = Join-Path $runDir 'fagbi_scope_z4xf_row.csv'
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null

    if (-not (Test-Z4xfShouldRun $rowCsv $Stage.stage)) {
        Write-Host "resume: $runName"
        return @(Import-Csv -LiteralPath $rowCsv)[0]
    }

    Write-Host "run: $runName (timeout ${CaseTimeoutSeconds}s)"
    $row = New-Z4xfRow -Fixture $Fixture -Stage $Stage -RunDir $runDir
    $argPath = Write-Z4xfArg -Fixture $Fixture -Stage $Stage -RunDir $runDir
    $bare = @(Get-Content -LiteralPath $argPath | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) {
        throw "Generated arg.in contains bare options: $($bare -join ', ')"
    }
    if ($PlanOnly) {
        $row.status = 'planned'
        $row.failure_class = 'plan_only'
        $row.wall_seconds = 0.0
        $row.peak_rss_mb = 0.0
        $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
        return $row
    }

    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
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
        $row = Read-Z4xfOutputs -RunDir $runDir -Row $row
    }
    $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
    return $row
}

function Write-Z4xfCapSummary([object[]]$Rows, [string]$Path) {
    $summary = New-Object System.Collections.Generic.List[object]
    foreach ($cap in @($Rows.ilu_neighbor_cap | Sort-Object -Unique)) {
        $set = @($Rows | Where-Object { [int]$_.ilu_neighbor_cap -eq [int]$cap })
        $pass = @($set | Where-Object { $_.failure_class -eq 'pass' })
        $summary.Add([pscustomobject]@{
            ilu_neighbor_cap = $cap
            rows = $set.Count
            pass_rows = $pass.Count
            timeout_rows = @($set | Where-Object { $_.status -eq 'timeout' }).Count
            max_audit = if ($pass.Count) { ($pass | ForEach-Object { Get-Z4xfDouble $_.sampled_audit_rms } | Measure-Object -Maximum).Maximum } else { [double]::NaN }
            max_sampled_direct_true_residual = if ($pass.Count) { ($pass | ForEach-Object { Get-Z4xfDouble $_.sampled_direct_true_residual } | Measure-Object -Maximum).Maximum } else { [double]::NaN }
            median_wall_seconds = if ($pass.Count) { @($pass | ForEach-Object { Get-Z4xfDouble $_.wall_seconds } | Sort-Object)[[int][Math]::Floor(($pass.Count - 1) / 2)] } else { [double]::NaN }
            median_matvec_seconds = if ($pass.Count) { @($pass | ForEach-Object { Get-Z4xfDouble $_.matvec_seconds_avg } | Sort-Object)[[int][Math]::Floor(($pass.Count - 1) / 2)] } else { [double]::NaN }
            max_rss_mb = if ($pass.Count) { ($pass | ForEach-Object { Get-Z4xfDouble $_.peak_rss_mb } | Measure-Object -Maximum).Maximum } else { [double]::NaN }
            median_nnz_blocks = if ($pass.Count) { @($pass | ForEach-Object { Get-Z4xfInt64 $_.preconditioner_sparse_nnz_blocks } | Sort-Object)[[int][Math]::Floor(($pass.Count - 1) / 2)] } else { 0 }
        }) | Out-Null
    }
    $summary.ToArray() | Export-Csv -LiteralPath $Path -NoTypeInformation
}

function Write-Z4xfSvg([object[]]$Rows, [string]$Path) {
    $passRows = @($Rows | Where-Object { $_.failure_class -eq 'pass' })
    $caps = @(32,64,128)
    $width = 780
    $height = 340
    $left = 70
    $bottom = 285
    $plotW = 660
    $plotH = 220
    $metrics = @{}
    foreach ($cap in $caps) {
        $set = @($passRows | Where-Object { [int]$_.ilu_neighbor_cap -eq $cap })
        $metrics[$cap] = if ($set.Count) {
            ($set | ForEach-Object { Get-Z4xfDouble $_.wall_seconds } | Measure-Object -Average).Average
        } else { 0.0 }
    }
    $maxY = [Math]::Max(1.0, ($metrics.Values | Measure-Object -Maximum).Maximum)
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add("<svg xmlns=`"http://www.w3.org/2000/svg`" width=`"$width`" height=`"$height`" viewBox=`"0 0 $width $height`">") | Out-Null
    $svg.Add('<rect width="100%" height="100%" fill="white"/>') | Out-Null
    $svg.Add("<line x1=`"$left`" y1=`"$bottom`" x2=`"$($left+$plotW)`" y2=`"$bottom`" stroke=`"#222`" stroke-width=`"1.2`"/>") | Out-Null
    $svg.Add("<line x1=`"$left`" y1=`"$($bottom-$plotH)`" x2=`"$left`" y2=`"$bottom`" stroke=`"#222`" stroke-width=`"1.2`"/>") | Out-Null
    $svg.Add("<text x=`"22`" y=`"54`" font-family=`"Arial`" font-size=`"13`" transform=`"rotate(-90 22,54)`">mean wall time (s)</text>") | Out-Null
    for ($i = 0; $i -lt $caps.Count; ++$i) {
        $cap = $caps[$i]
        $barW = 82
        $x = $left + 95 + $i * 190
        $barH = [Math]::Round(($metrics[$cap] / $maxY) * ($plotH - 20), 1)
        $y = $bottom - $barH
        $svg.Add("<rect x=`"$x`" y=`"$y`" width=`"$barW`" height=`"$barH`" fill=`"#4477AA`"/>") | Out-Null
        $svg.Add("<text x=`"$($x+$barW/2)`" y=`"$($bottom+24)`" text-anchor=`"middle`" font-family=`"Arial`" font-size=`"14`">cap $cap</text>") | Out-Null
        $svg.Add("<text x=`"$($x+$barW/2)`" y=`"$($y-8)`" text-anchor=`"middle`" font-family=`"Arial`" font-size=`"12`">$([Math]::Round($metrics[$cap],1))</text>") | Out-Null
    }
    $svg.Add("<text x=`"$($left+$plotW/2)`" y=`"326`" text-anchor=`"middle`" font-family=`"Arial`" font-size=`"13`">ILU(0) near-field neighbor cap</text>") | Out-Null
    $svg.Add('</svg>') | Out-Null
    $svg.ToArray() | Set-Content -LiteralPath $Path -Encoding ASCII
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$fixtureRows = Get-Z4xfFixtureRows
$defaultFixtures = if ($Quick) {
    @('1sh1')
} elseif ($Fixtures.Count -gt 0) {
    $Fixtures
} else {
    @('barnase','1a63','1sh1')
}
$allStages = @($IluCaps | Sort-Object -Unique | ForEach-Object { New-Z4xfStage ([int]$_) })
if ($allStages.Count -eq 0) { throw 'No ILU caps selected.' }

$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $defaultFixtures) {
    $fixture = Get-Z4xfFixture -Name $fixtureName -Rows $fixtureRows
    foreach ($stage in $allStages) {
        $rows.Add((Invoke-Z4xfCase -Fixture $fixture -Stage $stage)) | Out-Null
    }
}

$allRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $defaultFixtures) {
    foreach ($stage in $allStages) {
        $rowCsv = Join-Path (Join-Path $WorkRoot ('{0}_{1}' -f $fixtureName, $stage.stage)) 'fagbi_scope_z4xf_row.csv'
        if (Test-Path -LiteralPath $rowCsv) {
            $allRows.Add((@(Import-Csv -LiteralPath $rowCsv)[0])) | Out-Null
        }
    }
}
$matrix = @($allRows.ToArray() | Sort-Object fixture,ilu_neighbor_cap)
$matrix | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
Write-Z4xfCapSummary -Rows $matrix -Path $CapSummaryCsv
Write-Z4xfSvg -Rows $matrix -Path $PlotPath

$measured = @($matrix | Where-Object { $_.status -ne 'planned' })
$passes = @($matrix | Where-Object { $_.failure_class -eq 'pass' })
$timeouts = @($matrix | Where-Object { $_.status -eq 'timeout' })
$cap64Rows = @($matrix | Where-Object { [int]$_.ilu_neighbor_cap -eq 64 })
$cap64Passes = @($cap64Rows | Where-Object { $_.failure_class -eq 'pass' })
$alternativePasses = @($matrix | Where-Object {
    $_.failure_class -eq 'pass' -and [int]$_.ilu_neighbor_cap -ne 64
})
$accuracyGaps = @($matrix | Where-Object { $_.failure_class -eq 'operator_accuracy_gap' })

$decision = if ($PlanOnly) {
    'z4xf_plan_only'
} elseif ($timeouts.Count -gt 0) {
    'z4xf_runtime_limited'
} elseif ($cap64Rows.Count -gt 0 -and $cap64Passes.Count -ne $cap64Rows.Count) {
    'z4xf_cap64_gap'
} elseif ($alternativePasses.Count -gt 0 -and $cap64Passes.Count -eq $cap64Rows.Count) {
    'z4xf_cap_alternative_candidate'
} elseif ($accuracyGaps.Count -gt 0 -and $cap64Passes.Count -eq $cap64Rows.Count) {
    'z4xf_cap64_supported_with_alternative_gaps'
} elseif ($measured.Count -gt 0 -and $passes.Count -eq $measured.Count) {
    'z4xf_neighbor_cap_sweep_pass'
} elseif ($cap64Passes.Count -eq $cap64Rows.Count -and $cap64Rows.Count -gt 0) {
    'z4xf_cap64_supported'
} elseif ($measured.Count -gt 0) {
    'z4xf_neighbor_cap_sweep_ready_with_gaps'
} else {
    'z4xf_no_measured_rows'
}

[pscustomobject]@{
    decision = $decision
    quick = [string][bool]$Quick
    fixtures = ($defaultFixtures -join ';')
    ilu_caps = (@($allStages.ilu_neighbor_cap) -join ';')
    rows = $matrix.Count
    measured_rows = $measured.Count
    pass_rows = $passes.Count
    timeout_rows = $timeouts.Count
    accuracy_gap_rows = $accuracyGaps.Count
    case_timeout_seconds = $CaseTimeoutSeconds
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
    'scope: Z4X-F ILU neighbor-cap real sweep',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "fixtures: $($defaultFixtures -join ', ')",
    "ilu_caps: $(@($allStages.ilu_neighbor_cap) -join ', ')",
    "rows: $($matrix.Count)",
    "case_timeout_seconds: $CaseTimeoutSeconds",
    'resume_policy: per-row CSV files are reused unless -Force, -RerunFailed, or -RerunStages is passed',
    'production_default: grid',
    'tabi_treecode_default: degree-one; gmres-p4 explicit only',
    'diagnostic_flag: -pbe-tabi-fmm-local-ilu-cap-sweep=true'
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z4X-F ILU Neighbor-Cap Sweep Report',
    '',
    "Decision: $decision",
    '',
    "Matrix: $MatrixCsv",
    "Neighbor-cap summary: $CapSummaryCsv",
    "Plot: $PlotPath",
    '',
    "Fixtures: $($defaultFixtures -join ', ')",
    "ILU caps: $(@($allStages.ilu_neighbor_cap) -join ', ')",
    '',
    'This gate uses an explicit diagnostic flag so gmres-p4 remains guarded at ILU cap 64 outside the reviewer-facing neighbor-cap sweep.',
    '',
    'Resume: rerun the same command to continue missing rows; use -RerunFailed for failed/timeout rows, -RerunStages ilu_cap32, or -Force for all rows.'
)
$report | Set-Content -LiteralPath $RootReport -Encoding ASCII

Write-Host "Z4X-F decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host 'Resume: rerun the same command to continue missing rows; use -RerunFailed for failed/timeout rows, -RerunStages ilu_cap32, or -Force for all rows.'
if (-not $AllowFailedDecision -and $decision -notin @(
        'z4xf_neighbor_cap_sweep_pass',
        'z4xf_cap_alternative_candidate',
        'z4xf_cap64_supported',
        'z4xf_cap64_supported_with_alternative_gaps',
        'z4xf_plan_only')) {
    throw "Z4X-F decision was $decision"
}
