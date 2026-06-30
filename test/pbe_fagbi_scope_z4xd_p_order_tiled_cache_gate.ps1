param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$PrimaryFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$ExpansionFixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public_z4q'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4xd'),
    [string[]]$Fixtures = @(),
    [int]$CaseTimeoutSeconds = 1800,
    [int]$Threads = 4,
    [ValidateRange(64,8192)][int]$P5GlobalCacheBudgetMb = 4096,
    [switch]$Quick,
    [switch]$PlanOnly,
    [switch]$RefreshArtifactsOnly,
    [switch]$Force,
    [switch]$RerunFailed,
    [string[]]$RerunOrders = @(),
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

$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4xd_matrix.csv'
$OrderSummaryCsv = Join-Path $WorkRoot 'fagbi_scope_z4xd_order_summary.csv'
$BudgetCsv = Join-Path $WorkRoot 'fagbi_scope_z4xd_budget.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4xd_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4xd_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4xd_p_order_tiled_cache_plot.svg'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4XD_P5_TILED_CACHE_REPORT.txt'

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

function Get-Z4xdDouble([object]$Value, [double]$Default = [double]::NaN) {
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

function Get-Z4xdInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $parsed = [int64]0
    if ([int64]::TryParse([string]$Value, [Globalization.NumberStyles]::Integer,
            [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $Default
}

function Get-Z4xdProfileValue([object]$Profile, [string]$Name, [object]$Default = '') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Import-Z4xdManifestRows([string]$Root, [string]$ManifestName) {
    $path = Join-Path $Root $ManifestName
    if (-not (Test-Path -LiteralPath $path)) { return @() }
    return @(Import-Csv -LiteralPath $path | ForEach-Object {
        $_ | Add-Member -NotePropertyName fixture_root -NotePropertyValue $Root -Force
        $_
    })
}

function Get-Z4xdFixtureRows {
    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($row in (Import-Z4xdManifestRows $PrimaryFixtureRoot 'public_fixture_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    foreach ($row in (Import-Z4xdManifestRows $ExpansionFixtureRoot 'public_fixture_expansion_manifest.csv')) {
        $rows.Add($row) | Out-Null
    }
    return $rows.ToArray()
}

function Get-Z4xdFixture([string]$Name, [object[]]$Rows) {
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

function New-Z4xdStage(
    [int]$StageOrder,
    [string]$Stage,
    [int]$TaylorOrder,
    [string]$CacheLayout,
    [string]$CacheVariant,
    [bool]$ContextOnly = $false,
    [int]$CacheBudgetMb = 256
) {
    return [pscustomobject]@{
        order = $StageOrder
        stage = $Stage
        taylor_order = $TaylorOrder
        cache_layout = $CacheLayout
        cache_variant = $CacheVariant
        context_only = $ContextOnly
        cache_budget_mb = $CacheBudgetMb
    }
}

function Get-Z4xdStages {
    return @(
        (New-Z4xdStage 1 'p_order2' 2 'green-jet-compressed' 'global'),
        (New-Z4xdStage 2 'p_order3' 3 'green-jet-compressed' 'global'),
        (New-Z4xdStage 3 'p_order4' 4 'green-jet-compressed' 'global'),
        (New-Z4xdStage 4 'p_order5_global' 5 'green-jet-compressed' 'global' $true $P5GlobalCacheBudgetMb),
        (New-Z4xdStage 5 'p_order5' 5 'green-jet-tiled' 'tiled')
    )
}

function Convert-Z4xdOrderSelector([string]$Value) {
    $text = ([string]$Value).Trim().ToLowerInvariant()
    if ($text -in @('p_order5_global','p5_global','p5-global','order5_global','order5-global','5-global','global5')) {
        return 'p_order5_global'
    }
    if ($text -in @('p_order5_tiled','p5_tiled','p5-tiled','order5_tiled','order5-tiled','5-tiled','tiled5')) {
        return 'p_order5'
    }
    if ($text -match '^p?_?order([2-5])$') { return "p_order$($Matches[1])" }
    if ($text -match '^p([2-5])$') { return "p_order$($Matches[1])" }
    if ($text -match '^[2-5]$') { return "p_order$text" }
    return $text
}

function Set-Z4xdPropertyDefault([object]$Object, [string]$Name, [object]$Value) {
    if ($null -eq $Object.PSObject.Properties[$Name]) {
        $Object | Add-Member -NotePropertyName $Name -NotePropertyValue $Value -Force
    }
}

function Update-Z4xdRowStageMetadata([object]$Row, [object]$Stage, [string]$RunDir) {
    Set-Z4xdPropertyDefault -Object $Row -Name 'cache_variant' -Value $Stage.cache_variant
    Set-Z4xdPropertyDefault -Object $Row -Name 'context_only' -Value ([bool]$Stage.context_only)
    Set-Z4xdPropertyDefault -Object $Row -Name 'stage_order' -Value $Stage.order
    Set-Z4xdPropertyDefault -Object $Row -Name 'requested_cache_layout' -Value $Stage.cache_layout
    Set-Z4xdPropertyDefault -Object $Row -Name 'requested_cache_budget_mb' -Value $Stage.cache_budget_mb
    Set-Z4xdPropertyDefault -Object $Row -Name 'requested_taylor_order' -Value $Stage.taylor_order
    Set-Z4xdPropertyDefault -Object $Row -Name 'run_dir' -Value $RunDir
    $Row.cache_variant = $Stage.cache_variant
    $Row.context_only = [bool]$Stage.context_only
    $Row.stage_order = $Stage.order
    $Row.requested_cache_layout = $Stage.cache_layout
    $Row.requested_cache_budget_mb = $Stage.cache_budget_mb
    $Row.requested_taylor_order = $Stage.taylor_order
    $Row.run_dir = $RunDir
    return $Row
}

function Test-Z4xdContextOnly([object]$Row) {
    $prop = $Row.PSObject.Properties['context_only']
    return ($null -ne $prop -and [string]$prop.Value -match '^(True|true|1)$')
}

function New-Z4xdRow([object]$Fixture, [object]$Stage, [string]$RunDir) {
    return [pscustomobject]@{
        fixture = $Fixture.fixture
        pdb_id = $Fixture.pdb_id
        protein_name = $Fixture.protein_name
        atoms = Get-Z4xdInt64 $Fixture.atom_count
        order = $Stage.taylor_order
        stage = $Stage.stage
        stage_order = $Stage.order
        cache_variant = $Stage.cache_variant
        context_only = [bool]$Stage.context_only
        requested_taylor_order = $Stage.taylor_order
        effective_taylor_order = 0
        requested_cache_layout = $Stage.cache_layout
        requested_cache_budget_mb = $Stage.cache_budget_mb
        cache_layout = ''
        status = 'planned'
        failure_class = 'planned'
        exit_code = 0
        timed_out = 'False'
        panels = 0
        iterations = 0
        energy_kj_mol = [double]::NaN
        final_residual = [double]::NaN
        true_residual = [double]::NaN
        audit_rms = [double]::NaN
        sampled_direct_true_residual = [double]::NaN
        matvec_seconds_avg = [double]::NaN
        wall_seconds = [double]::NaN
        peak_rss_mb = [double]::NaN
        fmm_local_gmres_active = 'false'
        fmm_local_operator_finite = ''
        plan_build_count = 0
        cache_build_count = 0
        gmres_matvec_calls = 0
        cache_bytes = 0
        cache_operator_bytes = 0
        cache_tile_bytes = 0
        cache_tile_count = 0
        cache_tile_builds = 0
        cache_tile_seconds = [double]::NaN
        fallback_pairs = 0
        far_pairs = 0
        geometry_guard = ''
        geometry_theta = [double]::NaN
        run_dir = $RunDir
    }
}

function Write-Z4xdArg([object]$Fixture, [object]$Stage, [string]$RunDir) {
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
        -FmmLocalCache $Stage.cache_layout `
        -FmmLocalCacheBudgetMb $Stage.cache_budget_mb `
        -FmmLocalFallback packed-avx2 `
        -FmmLocalPlanBuilder metric-cache `
        -FmmLocalRefineProbeCount 4 `
        -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard strict-geometry-mac `
        -FmmLocalGeometryTheta 0.5
}

function Read-Z4xdAttribution([string]$RunDir, [object]$Row) {
    $summaryPath = Join-Path $RunDir 'pbe_tabi_operator_attribution_summary.csv'
    if (-not (Test-Path -LiteralPath $summaryPath)) { return $Row }
    $summary = @(Import-Csv -LiteralPath $summaryPath)
    $local = @($summary | Where-Object {
        $_.approximation_mode -like 'local_m2l*' -and
        $_.probe -eq 'solution' -and
        $_.row_component -eq 'both' -and
        $_.region -eq 'total' -and
        $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
    if ($null -eq $local) { return $Row }
    $Row.effective_taylor_order = Get-Z4xdInt64 $local.fmm_local_source_order $Row.effective_taylor_order
    $Row.fmm_local_operator_finite = [string]$local.fmm_local_finite_outputs
    $Row.plan_build_count = Get-Z4xdInt64 $local.fmm_local_plan_build_count $Row.plan_build_count
    $Row.cache_build_count = Get-Z4xdInt64 $local.fmm_local_cache_build_count $Row.cache_build_count
    $Row.cache_layout = [string]$local.fmm_local_cache_layout
    $Row.cache_bytes = Get-Z4xdInt64 $local.fmm_local_cache_bytes $Row.cache_bytes
    $Row.cache_operator_bytes = Get-Z4xdInt64 $local.fmm_local_cache_operator_bytes $Row.cache_operator_bytes
    $Row.cache_tile_bytes = Get-Z4xdInt64 $local.fmm_local_cache_tile_bytes $Row.cache_tile_bytes
    $Row.cache_tile_count = Get-Z4xdInt64 $local.fmm_local_cache_tile_count $Row.cache_tile_count
    $Row.cache_tile_builds = Get-Z4xdInt64 $local.fmm_local_cache_tile_builds $Row.cache_tile_builds
    $Row.cache_tile_seconds = Get-Z4xdDouble $local.fmm_local_cache_tile_seconds $Row.cache_tile_seconds
    $Row.fallback_pairs = Get-Z4xdInt64 $local.fmm_local_direct_fallback_pairs $Row.fallback_pairs
    $Row.far_pairs = Get-Z4xdInt64 $local.fmm_local_far_pairs_after $Row.far_pairs
    if (-not (Test-FmmCliFinite ([double]$Row.matvec_seconds_avg))) {
        $Row.matvec_seconds_avg = Get-Z4xdDouble $local.fmm_local_candidate_median_apply_seconds $Row.matvec_seconds_avg
    }
    return $Row
}

function Read-Z4xdMemoryProfile([string]$RunDir, [object]$Row) {
    $memoryPath = Join-Path $RunDir 'pbe_memory_profile.csv'
    if (-not (Test-Path -LiteralPath $memoryPath)) { return $Row }
    $mem = @(Import-Csv -LiteralPath $memoryPath)[0]
    $peakBytes = Get-Z4xdInt64 (Get-Z4xdProfileValue $mem 'peak_working_set_bytes' 0)
    if ($peakBytes -gt 0) {
        $Row.peak_rss_mb = [double]$peakBytes / 1048576.0
    }
    $totalSeconds = Get-Z4xdDouble (Get-Z4xdProfileValue $mem 'total_seconds' '')
    if (Test-FmmCliFinite $totalSeconds) {
        $Row.wall_seconds = $totalSeconds
    }
    return $Row
}

function Read-Z4xdOutputs([string]$RunDir, [object]$Row) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (Test-Path -LiteralPath $energyPath) {
        $energy = @(Import-Csv -LiteralPath $energyPath)[0]
        $Row.energy_kj_mol = Get-Z4xdDouble $energy.E_rxn_kJmol
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
    $Row.panels = Get-Z4xdInt64 $profile.panels
    $Row.iterations = Get-Z4xdInt64 $profile.iterations
    $Row.final_residual = Get-Z4xdDouble $profile.final_residual
    $Row.true_residual = Get-Z4xdDouble $profile.true_residual_rel
    $Row.audit_rms = Get-Z4xdDouble $profile.tree_vs_direct_matvec_rel
    $Row.sampled_direct_true_residual = Get-Z4xdDouble (Get-Z4xdProfileValue $profile 'fmm_local_sampled_direct_true_residual' '')
    $Row.matvec_seconds_avg = Get-Z4xdDouble $profile.matvec_seconds_avg
    $Row.effective_taylor_order = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_m2l_expansion_order' $Row.requested_taylor_order)
    $Row.fmm_local_gmres_active = [string](Get-Z4xdProfileValue $profile 'fmm_local_gmres_active' 'false')
    $Row.fmm_local_operator_finite = [string](Get-Z4xdProfileValue $profile 'fmm_local_operator_finite' '')
    $Row.plan_build_count = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_local_plan_build_count' 0)
    $Row.cache_build_count = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_local_cache_build_count' 0)
    $Row.gmres_matvec_calls = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_local_gmres_matvec_calls' 0)
    $Row.cache_layout = [string](Get-Z4xdProfileValue $profile 'fmm_local_cache_layout' '')
    $Row.cache_bytes = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_local_cache_bytes' 0)
    $Row.cache_operator_bytes = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_local_cache_operator_bytes' 0)
    $Row.cache_tile_bytes = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_local_cache_tile_bytes' 0)
    $Row.cache_tile_count = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_local_cache_tile_count' 0)
    $Row.cache_tile_builds = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_local_cache_tile_builds' 0)
    $Row.cache_tile_seconds = Get-Z4xdDouble (Get-Z4xdProfileValue $profile 'fmm_local_cache_tile_seconds' '')
    $Row.fallback_pairs = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_local_direct_fallback_pairs' 0)
    $Row.far_pairs = Get-Z4xdInt64 (Get-Z4xdProfileValue $profile 'fmm_local_far_pairs_after' 0)
    $Row.geometry_guard = [string](Get-Z4xdProfileValue $profile 'fmm_local_geometry_guard' '')
    $Row.geometry_theta = Get-Z4xdDouble (Get-Z4xdProfileValue $profile 'fmm_local_geometry_theta' '')
    $Row = Read-Z4xdAttribution -RunDir $RunDir -Row $Row
    $Row = Read-Z4xdMemoryProfile -RunDir $RunDir -Row $Row
    if ($Row.status -eq 'ok') {
        if (-not (Test-FmmCliFinite ([double]$Row.audit_rms)) -or
            -not (Test-FmmCliFinite ([double]$Row.sampled_direct_true_residual))) {
            $Row.failure_class = 'missing_operator_audit'
        } elseif ([double]$Row.true_residual -gt 1.0e-6) {
            $Row.failure_class = 'convergence_gap'
        } elseif ([double]$Row.audit_rms -gt 8.0e-3 -or
            [double]$Row.sampled_direct_true_residual -gt 8.0e-3) {
            $Row.failure_class = 'operator_accuracy_gap'
        } elseif ([int64]$Row.effective_taylor_order -ne [int64]$Row.requested_taylor_order) {
            $Row.failure_class = 'effective_order_mismatch'
        } elseif ([string]$Row.cache_layout -ne [string]$Row.requested_cache_layout) {
            $Row.failure_class = 'cache_layout_mismatch'
        } elseif ([int64]$Row.plan_build_count -ne 1 -or
            [int64]$Row.cache_build_count -ne 1 -or
            [int64]$Row.gmres_matvec_calls -le 0 -or
            $Row.fmm_local_gmres_active -ne 'true' -or
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

function Update-Z4xdFailureFromLogs([string]$RunDir, [object]$Row) {
    if ($Row.PSObject.Properties.Name -contains 'exit_code' -and
        [string]$Row.exit_code -eq '-1') {
        $Row.status = 'interrupted'
        $Row.failure_class = 'runtime_limited'
    }
    $stderrPath = Join-Path $RunDir 'stderr.log'
    if (Test-Path -LiteralPath $stderrPath) {
        $stderr = Get-Content -LiteralPath $stderrPath -Raw
        if ($stderr -match 'cache_budget_exceeded') {
            $Row.status = 'failed'
            $Row.failure_class = 'memory_limited'
            $required = [regex]::Match($stderr, 'required=([0-9]+)')
            if ($required.Success) { $Row.cache_bytes = [int64]$required.Groups[1].Value }
            if ([int64]$Row.effective_taylor_order -eq 0) { $Row.effective_taylor_order = $Row.requested_taylor_order }
        }
    }
    return $Row
}

function Test-Z4xdRowFinal([object]$Row) {
    return $Row.failure_class -eq 'pass'
}

function Test-Z4xdShouldRun([string]$RowCsv, [object]$Stage) {
    if ($Force) { return $true }
    $selectedStages = @($RerunStages | ForEach-Object { Convert-Z4xdOrderSelector $_ })
    $selectedOrders = @($RerunOrders | ForEach-Object { Convert-Z4xdOrderSelector $_ })
    if (-not (Test-Path -LiteralPath $RowCsv)) { return $true }
    try {
        $row = @(Import-Csv -LiteralPath $RowCsv)[0]
        $row = Update-Z4xdFailureFromLogs -RunDir (Split-Path -Parent $RowCsv) -Row $row
        $row | Export-Csv -LiteralPath $RowCsv -NoTypeInformation
        if (Test-Path -LiteralPath (Join-Path (Split-Path -Parent $RowCsv) 'pbe_tabi_profile.csv')) {
            return $false
        }
        if ($row.status -eq 'planned') { return $true }
        if ($RerunFailed) {
            return -not (Test-Z4xdRowFinal $row)
        }
        return $false
    } catch {
        return $true
    }
}

function Invoke-Z4xdCase([object]$Fixture, [object]$Stage) {
    $runName = '{0}_{1}' -f $Fixture.fixture, $Stage.stage
    $runDir = Join-Path $WorkRoot $runName
    $rowCsv = Join-Path $runDir 'fagbi_scope_z4xd_row.csv'
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null

    if (-not (Test-Z4xdShouldRun -RowCsv $rowCsv -Stage $Stage)) {
        Write-Host "resume: $runName"
        $row = @(Import-Csv -LiteralPath $rowCsv)[0]
        $row = Update-Z4xdRowStageMetadata -Row $row -Stage $Stage -RunDir $runDir
        if (Test-Path -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv')) {
            $row.status = 'ok'
            $row.failure_class = 'pass'
            $row = Read-Z4xdOutputs -RunDir $runDir -Row $row
            $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
        }
        return $row
    }

    Write-Host "run: $runName (timeout ${CaseTimeoutSeconds}s)"
    $row = New-Z4xdRow -Fixture $Fixture -Stage $Stage -RunDir $runDir
    $argPath = Write-Z4xdArg -Fixture $Fixture -Stage $Stage -RunDir $runDir
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
        $row.failure_class = if ($run.exit_code -eq -1) { 'runtime_limited' } else { "exit_code_$($run.exit_code)" }
        $row = Update-Z4xdFailureFromLogs -RunDir $runDir -Row $row
    } else {
        $row.status = 'ok'
        $row.failure_class = 'pass'
        $row = Read-Z4xdOutputs -RunDir $runDir -Row $row
    }
    $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
    return $row
}

function Get-Z4xdMedian([double[]]$Values) {
    $vals = @($Values | Where-Object { Test-FmmCliFinite $_ } | Sort-Object)
    if ($vals.Count -eq 0) { return [double]::NaN }
    return [double]$vals[[int][Math]::Floor($vals.Count / 2)]
}

function Write-Z4xdOrderSummary([object[]]$Rows, [string]$Path) {
    $summary = New-Object System.Collections.Generic.List[object]
    foreach ($stage in @($Rows.stage | Sort-Object -Unique)) {
        $set = @($Rows | Where-Object { $_.stage -eq $stage })
        $pass = @($set | Where-Object { $_.failure_class -eq 'pass' })
        $summary.Add([pscustomobject]@{
            stage = $stage
            requested_taylor_order = (@($set | Select-Object -First 1).requested_taylor_order)
            requested_cache_layout = (@($set | Select-Object -First 1).requested_cache_layout)
            requested_cache_budget_mb = (@($set | Select-Object -First 1).requested_cache_budget_mb)
            cache_variant = (@($set | Select-Object -First 1).cache_variant)
            context_only = (@($set | Select-Object -First 1).context_only)
            rows = $set.Count
            pass_rows = $pass.Count
            max_audit = if ($pass.Count) { ($pass | ForEach-Object { Get-Z4xdDouble $_.audit_rms } | Measure-Object -Maximum).Maximum } else { [double]::NaN }
            max_sampled_direct_true_residual = if ($pass.Count) { ($pass | ForEach-Object { Get-Z4xdDouble $_.sampled_direct_true_residual } | Measure-Object -Maximum).Maximum } else { [double]::NaN }
            median_matvec_seconds = Get-Z4xdMedian @($pass | ForEach-Object { Get-Z4xdDouble $_.matvec_seconds_avg })
            median_cache_bytes = Get-Z4xdMedian @($pass | ForEach-Object { [double](Get-Z4xdInt64 $_.cache_bytes) })
            median_peak_rss_mb = Get-Z4xdMedian @($pass | ForEach-Object { Get-Z4xdDouble $_.peak_rss_mb })
            max_fallback_pairs = if ($pass.Count) { ($pass | ForEach-Object { Get-Z4xdInt64 $_.fallback_pairs } | Measure-Object -Maximum).Maximum } else { 0 }
        }) | Out-Null
    }
    $summary.ToArray() | Sort-Object requested_taylor_order | Export-Csv -LiteralPath $Path -NoTypeInformation
}

function Write-Z4xdSvg([object[]]$Rows, [string]$Path) {
    $categories = @(
        [pscustomobject]@{ stage = 'p_order2'; label = 'p2'; color = '#56B4E9' },
        [pscustomobject]@{ stage = 'p_order3'; label = 'p3'; color = '#009E73' },
        [pscustomobject]@{ stage = 'p_order4'; label = 'p4'; color = '#0072B2' },
        [pscustomobject]@{ stage = 'p_order5_global'; label = 'p5G'; color = '#E69F00' },
        [pscustomobject]@{ stage = 'p_order5'; label = 'p5T'; color = '#D55E00' }
    )
    $stats = New-Object System.Collections.Generic.List[object]
    foreach ($cat in $categories) {
        $set = @($Rows | Where-Object { $_.stage -eq $cat.stage })
        $pass = @($set | Where-Object { $_.failure_class -eq 'pass' })
        $memValues = @()
        foreach ($row in $set) {
            $rss = Get-Z4xdDouble $row.peak_rss_mb
            $cacheMb = [double](Get-Z4xdInt64 $row.cache_bytes) / 1048576.0
            if ($row.failure_class -eq 'memory_limited' -and $cacheMb -gt 0.0) {
                $memValues += $cacheMb
            } elseif (Test-FmmCliFinite $rss) {
                $memValues += $rss
            } elseif ($cacheMb -gt 0.0) {
                $memValues += $cacheMb
            }
        }
        $maxAudit = if ($pass.Count) {
            ($pass | ForEach-Object { Get-Z4xdDouble $_.audit_rms } | Measure-Object -Maximum).Maximum
        } else { [double]::NaN }
        $stats.Add([pscustomobject]@{
            stage = $cat.stage
            label = $cat.label
            color = $cat.color
            pass_rows = $pass.Count
            rows = $set.Count
            max_audit = $maxAudit
            accuracy_score = if (Test-FmmCliFinite $maxAudit) { -[Math]::Log10([Math]::Max(1.0e-15, [double]$maxAudit)) } else { [double]::NaN }
            median_matvec = Get-Z4xdMedian @($pass | ForEach-Object { Get-Z4xdDouble $_.matvec_seconds_avg })
            max_memory_mb = if ($memValues.Count) { ($memValues | Measure-Object -Maximum).Maximum } else { [double]::NaN }
        }) | Out-Null
    }
    if ((@($stats.ToArray() | Where-Object {
                (Test-FmmCliFinite (Get-Z4xdDouble $_.accuracy_score)) -or
                (Test-FmmCliFinite (Get-Z4xdDouble $_.median_matvec)) -or
                (Test-FmmCliFinite (Get-Z4xdDouble $_.max_memory_mb))
            })).Count -eq 0) {
        '<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="620"><text x="24" y="42">No p-order data</text></svg>' |
            Set-Content -LiteralPath $Path -Encoding ASCII
        return
    }
    $maxAccuracy = [Math]::Max(1.0, ($stats.ToArray() | ForEach-Object { Get-Z4xdDouble $_.accuracy_score } | Measure-Object -Maximum).Maximum)
    $maxRuntime = [Math]::Max(1.0e-6, ($stats.ToArray() | ForEach-Object { Get-Z4xdDouble $_.median_matvec } | Measure-Object -Maximum).Maximum)
    $maxMemory = [Math]::Max(1.0, ($stats.ToArray() | ForEach-Object { Get-Z4xdDouble $_.max_memory_mb } | Measure-Object -Maximum).Maximum)
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="620" viewBox="0 0 1200 620">')
    $svg.Add('<rect width="1200" height="620" fill="white"/>')
    $panelX = @(62, 442, 822)
    $panelLetters = @('A','B','C')
    $axisLabels = @('-log10 audit RMS','median matvec (s)','RSS / cache MB')
    $maxima = @($maxAccuracy, $maxRuntime, $maxMemory)
    $metrics = @('accuracy_score','median_matvec','max_memory_mb')
    for ($p = 0; $p -lt 3; ++$p) {
        $left = $panelX[$p]
        $top = 62
        $height = 390
        $width = 306
        $svg.Add(('<text x="{0}" y="34" font-family="Arial" font-size="24" font-weight="700">{1}</text>' -f $left,$panelLetters[$p]))
        $svg.Add(('<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="#222" stroke-width="1.4"/>' -f $left,($top + $height),($left + $width)))
        $svg.Add(('<line x1="{0}" y1="{1}" x2="{0}" y2="{2}" stroke="#222" stroke-width="1.4"/>' -f $left,$top,($top + $height)))
        $svg.Add(('<text x="{0}" y="{1}" font-family="Arial" font-size="13" transform="rotate(-90 {0},{1})">{2}</text>' -f ($left - 38),($top + 250),$axisLabels[$p]))
        $svg.Add(('<text x="{0}" y="{1}" font-family="Arial" font-size="11" fill="#555">{2:F2}</text>' -f ($left + 4),($top + 12),[double]$maxima[$p]))
        $step = 57
        for ($i = 0; $i -lt $stats.Count; ++$i) {
            $stat = $stats[$i]
            $value = Get-Z4xdDouble $stat.PSObject.Properties[$metrics[$p]].Value
            $x = $left + 34 + $i * $step
            if (Test-FmmCliFinite $value) {
                $h = [Math]::Max(2.0, ($height - 30.0) * ($value / [Math]::Max(1.0e-12, [double]$maxima[$p])))
                $y = $top + $height - $h
                $svg.Add(('<rect x="{0:F1}" y="{1:F1}" width="34" height="{2:F1}" fill="{3}"/>' -f $x,$y,$h,$stat.color))
                $label = if ($p -eq 0) { ('{0:F1}' -f $value) } elseif ($p -eq 1) { ('{0:F2}' -f $value) } else { ('{0:F0}' -f $value) }
                $svg.Add(('<text x="{0:F1}" y="{1:F1}" font-family="Arial" font-size="9" text-anchor="middle">{2}</text>' -f ($x + 17),($y - 5),$label))
            } else {
                $svg.Add(('<text x="{0:F1}" y="{1:F1}" font-family="Arial" font-size="10" fill="#888" text-anchor="middle">NA</text>' -f ($x + 17),($top + $height - 12)))
            }
            $svg.Add(('<text x="{0:F1}" y="{1:F1}" font-family="Arial" font-size="11" text-anchor="middle">{2}</text>' -f ($x + 17),($top + $height + 22),$stat.label))
        }
    }
    $svg.Add('</svg>')
    $svg | Set-Content -LiteralPath $Path -Encoding ASCII
}

$fixtureRows = Get-Z4xdFixtureRows
$fixtureNames = if ($Fixtures.Count -gt 0) {
    $Fixtures
} elseif ($Quick) {
    @('1sh1','1a63')
} else {
    @('barnase','1a63','1sh1','1ubq','2ci2','1vii')
}
$fixtureNames = @($fixtureNames)
$allStages = Get-Z4xdStages
$stagesToRun = $allStages
if ($RerunStages.Count -gt 0 -or $RerunOrders.Count -gt 0) {
    $knownStages = @($allStages | ForEach-Object { $_.stage })
    $selected = @($RerunStages + $RerunOrders | ForEach-Object { Convert-Z4xdOrderSelector $_ })
    $unknownStages = @($selected | Where-Object { $knownStages -notcontains $_ })
    if ($unknownStages.Count -gt 0) {
        throw "Unknown rerun order/stage value(s): $($unknownStages -join ', '). Known values: $($knownStages -join ', '), 2, 3, 4, 5, p2, p3, p4, p5, p5-global, p5-tiled."
    }
    $stagesToRun = @($allStages | Where-Object { $selected -contains $_.stage })
}

$rows = New-Object System.Collections.Generic.List[object]
if (-not $RefreshArtifactsOnly) {
    foreach ($fixtureName in $fixtureNames) {
        $fixture = Get-Z4xdFixture -Name $fixtureName -Rows $fixtureRows
        foreach ($stage in $stagesToRun) {
            $rows.Add((Invoke-Z4xdCase -Fixture $fixture -Stage $stage)) | Out-Null
        }
    }
}

$allRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    foreach ($stage in $allStages) {
        $rowCsv = Join-Path (Join-Path $WorkRoot ('{0}_{1}' -f $fixtureName, $stage.stage)) 'fagbi_scope_z4xd_row.csv'
        if (Test-Path -LiteralPath $rowCsv) {
            $runDir = Split-Path -Parent $rowCsv
            $row = Update-Z4xdRowStageMetadata -Row (@(Import-Csv -LiteralPath $rowCsv)[0]) -Stage $stage -RunDir $runDir
            if (Test-Path -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv')) {
                $row.status = 'ok'
                $row.failure_class = 'pass'
                $row = Read-Z4xdOutputs -RunDir $runDir -Row $row
                $row | Export-Csv -LiteralPath $rowCsv -NoTypeInformation
            }
            $allRows.Add($row) | Out-Null
        }
    }
}
$matrix = @($allRows.ToArray() | Sort-Object fixture,stage_order)
$matrix | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
Write-Z4xdOrderSummary -Rows $matrix -Path $OrderSummaryCsv

$budgetRows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    $p4 = @($matrix | Where-Object { $_.fixture -eq $fixtureName -and $_.stage -eq 'p_order4' }) | Select-Object -First 1
    $p5 = @($matrix | Where-Object { $_.fixture -eq $fixtureName -and $_.stage -eq 'p_order5' }) | Select-Object -First 1
    $p5Global = @($matrix | Where-Object { $_.fixture -eq $fixtureName -and $_.stage -eq 'p_order5_global' }) | Select-Object -First 1
    if ($p4 -and $p5) {
        $budgetRows.Add([pscustomobject]@{
            fixture = $fixtureName
            p4_status = $p4.status
            p5_status = $p5.status
            p5_global_status = if ($p5Global) { $p5Global.status } else { 'missing' }
            p4_failure = $p4.failure_class
            p5_failure = $p5.failure_class
            p5_global_failure = if ($p5Global) { $p5Global.failure_class } else { 'missing' }
            p4_matvec_seconds = $p4.matvec_seconds_avg
            p5_matvec_seconds = $p5.matvec_seconds_avg
            p5_global_matvec_seconds = if ($p5Global) { $p5Global.matvec_seconds_avg } else { [double]::NaN }
            p5_vs_p4_matvec_ratio = if ((Get-Z4xdDouble $p4.matvec_seconds_avg) -gt 0.0) { (Get-Z4xdDouble $p5.matvec_seconds_avg) / (Get-Z4xdDouble $p4.matvec_seconds_avg) } else { [double]::NaN }
            p4_cache_bytes = $p4.cache_bytes
            p5_cache_bytes = $p5.cache_bytes
            p5_global_cache_bytes = if ($p5Global) { $p5Global.cache_bytes } else { 0 }
            p4_cache_budget_mb = $p4.requested_cache_budget_mb
            p5_cache_budget_mb = $p5.requested_cache_budget_mb
            p5_global_cache_budget_mb = if ($p5Global) { $p5Global.requested_cache_budget_mb } else { 0 }
            p5_vs_p4_cache_ratio = if ((Get-Z4xdDouble $p4.cache_bytes) -gt 0.0) { (Get-Z4xdDouble $p5.cache_bytes) / (Get-Z4xdDouble $p4.cache_bytes) } else { [double]::NaN }
            p5_global_vs_p4_cache_ratio = if ($p5Global -and (Get-Z4xdDouble $p4.cache_bytes) -gt 0.0) { (Get-Z4xdDouble $p5Global.cache_bytes) / (Get-Z4xdDouble $p4.cache_bytes) } else { [double]::NaN }
            p4_fallback_pairs = $p4.fallback_pairs
            p5_fallback_pairs = $p5.fallback_pairs
            p5_global_fallback_pairs = if ($p5Global) { $p5Global.fallback_pairs } else { 0 }
            p5_vs_p4_fallback_ratio = if ((Get-Z4xdDouble $p4.fallback_pairs) -gt 0.0) { (Get-Z4xdDouble $p5.fallback_pairs) / (Get-Z4xdDouble $p4.fallback_pairs) } else { [double]::NaN }
        }) | Out-Null
    }
}
$budgetRows.ToArray() | Export-Csv -LiteralPath $BudgetCsv -NoTypeInformation
Write-Z4xdSvg -Rows $matrix -Path $PlotPath

$candidateMatrix = @($matrix | Where-Object { -not (Test-Z4xdContextOnly $_) })
$contextMatrix = @($matrix | Where-Object { Test-Z4xdContextOnly $_ })
$measured = @($matrix | Where-Object { $_.status -notin @('planned') })
$timeouts = @($candidateMatrix | Where-Object { $_.failure_class -eq 'runtime_limited' })
$memoryRows = @($candidateMatrix | Where-Object { $_.failure_class -eq 'memory_limited' })
$accuracyGapRows = @($candidateMatrix | Where-Object { $_.failure_class -in @('operator_accuracy_gap','convergence_gap','missing_operator_audit','effective_order_mismatch','cache_layout_mismatch','lifecycle_gap') })
$p4Pass = @($matrix | Where-Object { $_.stage -eq 'p_order4' -and $_.failure_class -eq 'pass' })
$p5Pass = @($matrix | Where-Object { $_.stage -eq 'p_order5' -and $_.failure_class -eq 'pass' })
$p4MedianMatvec = Get-Z4xdMedian @($p4Pass | ForEach-Object { Get-Z4xdDouble $_.matvec_seconds_avg })
$p5MedianMatvec = Get-Z4xdMedian @($p5Pass | ForEach-Object { Get-Z4xdDouble $_.matvec_seconds_avg })
$p5MaxAudit = if ($p5Pass.Count) { ($p5Pass | ForEach-Object { Get-Z4xdDouble $_.audit_rms } | Measure-Object -Maximum).Maximum } else { [double]::NaN }
$requiredFixtureCount = $fixtureNames.Count

$decision = if ($PlanOnly) {
    'z4xd_plan_only'
} elseif ($RefreshArtifactsOnly) {
    'z4xd_artifacts_refreshed'
} elseif ($timeouts.Count -gt 0) {
    'z4xd_runtime_limited'
} elseif ($memoryRows.Count -gt 0 -and $p5Pass.Count -lt $requiredFixtureCount) {
    'z4xd_p5_memory_limited'
} elseif ($accuracyGapRows.Count -gt 0) {
    'z4xd_accuracy_gap'
} elseif ($p5Pass.Count -eq $requiredFixtureCount -and
    (Test-FmmCliFinite $p5MedianMatvec) -and
    (Test-FmmCliFinite $p4MedianMatvec) -and
    $p5MedianMatvec -lt $p4MedianMatvec -and
    $p5MaxAudit -le 8.0e-3) {
    'z4xd_p5_speed_pass'
} elseif ($p5Pass.Count -eq $requiredFixtureCount -and
    $p5MaxAudit -le 8.0e-3 -and
    (Test-FmmCliFinite $p5MedianMatvec) -and
    (Test-FmmCliFinite $p4MedianMatvec) -and
    $p5MedianMatvec -ge $p4MedianMatvec) {
    'z4xd_p5_accuracy_pass_runtime_gap'
} elseif ($p4Pass.Count -eq $requiredFixtureCount) {
    'z4xd_p4_remains_best'
} else {
    'z4xd_harness_gap'
}

[pscustomobject]@{
    decision = $decision
    quick = [bool]$Quick
    refresh_artifacts_only = [bool]$RefreshArtifactsOnly
    fixtures = ($fixtureNames -join ';')
    rows = $matrix.Count
    measured_rows = $measured.Count
    p4_pass_rows = $p4Pass.Count
    p5_pass_rows = $p5Pass.Count
    context_rows = $contextMatrix.Count
    timeout_rows = $timeouts.Count
    memory_rows = $memoryRows.Count
    accuracy_gap_rows = $accuracyGapRows.Count
    p4_median_matvec_seconds = $p4MedianMatvec
    p5_median_matvec_seconds = $p5MedianMatvec
    p5_max_audit = $p5MaxAudit
    case_timeout_seconds = $CaseTimeoutSeconds
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

@(
    'scope: Z4X-D p5 tiled Green-jet cache p-order gate',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "refresh_artifacts_only: $([bool]$RefreshArtifactsOnly)",
    "fixtures: $($fixtureNames -join ', ')",
    'protocol: full protein, guarded gmres-p4, theta 0.6, strict geometry theta 0.5, operator-tail-static tol 0.002, metric-cache plan builder, packed-avx2 fallback, ILU0 cap 64',
    "cache_policy: p2/p3/p4 use green-jet-compressed with 256 MB budget; p5 tiled is the 256 MB candidate; p5 global compressed is context-only RAM evidence with $P5GlobalCacheBudgetMb MB budget",
    'production_default: grid',
    'tabi_treecode_default: degree-one',
    'gmres_p4_name: retained for compatibility; effective Taylor order is reported separately',
    "matrix: $MatrixCsv",
    "order_summary: $OrderSummaryCsv",
    "budget: $BudgetCsv"
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII

$report = New-Object System.Collections.Generic.List[string]
$report.Add('# TITAN FMM/M2L Scope Z4X-D p5 Tiled Cache Report')
$report.Add('')
$report.Add("Decision: $decision")
$report.Add('')
$report.Add('Purpose')
$report.Add('Measure whether a tiled/streamed p5 Green-jet cache removes the p5 memory bottleneck and whether p5 improves runtime relative to p2/p3/p4 under one guarded local FMM/M2L protocol.')
$report.Add('')
$report.Add('Protocol')
$report.Add('- Full public protein PQR fixtures; no atom sampling.')
$report.Add('- Analytical FMM/M2L, explicit gmres-p4 path, Taylor order varied from p2 through p5.')
$report.Add("- p2/p3/p4 use the existing global compressed Green-jet cache with 256 MB budget; p5 is measured both as a context-only global compressed cache row with $P5GlobalCacheBudgetMb MB budget and as the 256 MB tiled-cache candidate.")
$report.Add('- Static tail tolerance 0.002, strict geometry theta 0.5, packed AVX2 fallback, sparse near-field ILU(0) neighbor cap 64.')
$report.Add('- Grid remains production default; TABI treecode default remains degree-one.')
$report.Add('')
$report.Add('Rows')
foreach ($row in $matrix) {
    $report.Add(('- {0} p{1} {2}: status={3}, failure={4}, context_only={5}, cache={6}/{7}, budget={8} MB, residual={9}, audit={10}, direct_residual={11}, fallback_pairs={12}, far_pairs={13}, matvec={14}s, cache_bytes={15}, rss={16} MB' -f
            $row.fixture,
            $row.requested_taylor_order,
            $row.cache_variant,
            $row.status,
            $row.failure_class,
            $row.context_only,
            $row.requested_cache_layout,
            $row.cache_layout,
            $row.requested_cache_budget_mb,
            $row.true_residual,
            $row.audit_rms,
            $row.sampled_direct_true_residual,
            $row.fallback_pairs,
            $row.far_pairs,
            $row.matvec_seconds_avg,
            $row.cache_bytes,
            $row.peak_rss_mb))
}
$report.Add('')
$report.Add('Artifacts')
$report.Add("Matrix: $MatrixCsv")
$report.Add("Order summary: $OrderSummaryCsv")
$report.Add("Budget summary: $BudgetCsv")
$report.Add("Decision: $DecisionCsv")
$report.Add("Plot: $PlotPath")
$report | Set-Content -LiteralPath $RootReport -Encoding ASCII

Write-Host "Z4X-D decision: $decision"
Write-Host "Matrix: $MatrixCsv"
Write-Host "Resume: rerun the same command to continue missing/mismatched rows; use -RerunFailed for failed/timeout rows, -RerunOrders 5, or -Force for all rows."

if ($decision -notin @('z4xd_p5_speed_pass','z4xd_p5_accuracy_pass_runtime_gap','z4xd_p4_remains_best','z4xd_plan_only','z4xd_artifacts_refreshed') -and
    -not $AllowFailedDecision) {
    throw "Z4X-D p5 tiled cache gate did not pass: $decision"
}
