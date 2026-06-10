param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_phase7_gate'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 120,
    [switch]$AllowFailedStatus
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$PerfMatrixScript = Join-Path $PSScriptRoot 'pbe_tabi_perf_matrix.ps1'
$TriageHelperScript = Join-Path $PSScriptRoot 'pbe_tabi_triage_helpers.ps1'

function Convert-ToDoubleOrNull([object]$Value) {
    if ($null -eq $Value) { return $null }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $null }
    $parsed = 0.0
    if ([double]::TryParse($text, [System.Globalization.NumberStyles]::Float,
            $InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $null
}

function Format-Number([object]$Value) {
    if ($null -eq $Value) { return 'NA' }
    if ($Value -is [string]) { return $Value }
    return ([double]$Value).ToString('G17', $InvariantCulture)
}

function Get-RowValue([object]$Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside Phase 7 gate root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Get-EnergyStatus($Row, $GridByFixture) {
    if ((Get-RowValue $Row 'status') -ne 'ok') { return 'failed_status' }
    if ((Get-RowValue $Row 'mode') -eq 'grid') { return 'grid_reference' }
    $fixture = Get-RowValue $Row 'fixture'
    if (-not $GridByFixture.ContainsKey($fixture)) { return 'missing_grid_reference' }
    $gridEnergy = Convert-ToDoubleOrNull (Get-RowValue $GridByFixture[$fixture] 'energy_kJmol')
    $energy = Convert-ToDoubleOrNull (Get-RowValue $Row 'energy_kJmol')
    if ($null -eq $gridEnergy -or $null -eq $energy) { return 'missing_energy' }
    $denom = [Math]::Max(1.0e-12, [Math]::Abs($gridEnergy))
    $rel = [Math]::Abs($energy - $gridEnergy) / $denom
    if ($rel -le 5.0e-3) { return 'ok' }
    return 'energy_drift_gt_5e-3'
}

function Get-EnergyRelVsGrid($Row, $GridByFixture) {
    $fixture = Get-RowValue $Row 'fixture'
    if (-not $GridByFixture.ContainsKey($fixture)) { return $null }
    $gridEnergy = Convert-ToDoubleOrNull (Get-RowValue $GridByFixture[$fixture] 'energy_kJmol')
    $energy = Convert-ToDoubleOrNull (Get-RowValue $Row 'energy_kJmol')
    if ($null -eq $gridEnergy -or $null -eq $energy) { return $null }
    return [Math]::Abs($energy - $gridEnergy) / [Math]::Max(1.0e-12, [Math]::Abs($gridEnergy))
}

function Get-AcceptedMissStatus($Row) {
    $names = @('accepted_missed_total',
               'convergence_estimator_high_error_missed',
               'convergence_k1_high_error_missed',
               'convergence_k2_high_error_missed',
               'convergence_k3_high_error_missed',
               'convergence_k4_high_error_missed')
    $seen = $false
    $total = 0.0
    foreach ($name in $names) {
        $value = Convert-ToDoubleOrNull (Get-RowValue $Row $name)
        if ($null -ne $value) {
            $seen = $true
            $total += $value
        }
    }
    if (-not $seen) { return 'unknown' }
    if ($total -eq 0.0) { return 'ok' }
    return 'accepted_miss_detected'
}

function Get-FailureMode($Row, [string]$EnergyStatus, [string]$AcceptedMissStatus) {
    $status = Get-RowValue $Row 'status'
    if ($status -ne 'ok') { return $status }
    if ($EnergyStatus -notin @('ok','grid_reference')) { return $EnergyStatus }
    if ($AcceptedMissStatus -eq 'accepted_miss_detected') { return $AcceptedMissStatus }
    $classification = Get-RowValue $Row 'performance_classification' 'ok'
    if ($classification -eq 'grid_reference' -or $classification -eq 'direct_reference') { return 'none' }
    return $classification
}

function New-NormalizedRow($Row, $GridByFixture) {
    $mode = Get-RowValue $Row 'mode'
    $fixture = Get-RowValue $Row 'fixture'
    $energyStatus = Get-EnergyStatus $Row $GridByFixture
    $acceptedMissStatus = if ($mode -eq 'tabi_treecode') { Get-AcceptedMissStatus $Row } else { 'NA' }
    $ramBytes = Convert-ToDoubleOrNull (Get-RowValue $Row 'peak_ram_bytes')
    $energyRel = Get-EnergyRelVsGrid $Row $GridByFixture
    $triage = New-TabiTriageRecord -Row $Row -RunDir (Get-RowValue $Row 'run_dir') -ExistingClassification (Get-RowValue $Row 'performance_classification')
    [PSCustomObject]@{
        fixture = $fixture
        mode = $mode
        status = Get-RowValue $Row 'status'
        wall_seconds = Get-RowValue $Row 'total_seconds'
        peak_rss_mb = if ($null -ne $ramBytes) { Format-Number ($ramBytes / 1048576.0) } else { 'NA' }
        peak_ram_bytes = Get-RowValue $Row 'peak_ram_bytes'
        energy_status = $energyStatus
        energy_rel_vs_grid = Format-Number $energyRel
        gmres_iterations = Get-RowValue $Row 'gmres_iterations'
        final_residual = Get-RowValue $Row 'final_residual'
        accepted_miss_status = $acceptedMissStatus
        time_ratio_vs_grid = Get-RowValue $Row 'time_ratio_vs_grid'
        ram_ratio_vs_grid = Get-RowValue $Row 'ram_ratio_vs_grid'
        failure_mode = Get-FailureMode $Row $energyStatus $acceptedMissStatus
        triage_classification = $triage.triage_classification
        seconds_per_matvec = $triage.seconds_per_matvec
        leaf_panel_interactions_per_matvec = $triage.leaf_panel_interactions_per_matvec
        leaf_vs_n2_ratio = $triage.leaf_vs_n2_ratio
        residual_delta_per_matvec = $triage.residual_delta_per_matvec
        preconditioner_requested = $triage.preconditioner_requested
        preconditioner_effective = $triage.preconditioner_effective
        preconditioner_downgrade = $triage.preconditioner_downgrade
        prepare_progress_phase = $triage.prepare_progress_phase
        prepare_progress_status = $triage.prepare_progress_status
        prepare_progress_rss_bytes = $triage.prepare_progress_rss_bytes
        triage_source = $triage.triage_source
        decision_reason = 'input_row'
        surface_points = Get-RowValue $Row 'surface_points'
        tree_leaf_size = Get-RowValue $Row 'tree_leaf_size'
        opening_theta = Get-RowValue $Row 'opening_theta'
        performance_classification = Get-RowValue $Row 'performance_classification'
        run_dir = Get-RowValue $Row 'run_dir'
    }
}

function Select-CandidateRows($Rows) {
    $selected = New-Object System.Collections.Generic.List[object]
    foreach ($group in ($Rows | Where-Object { $_.mode -eq 'tabi_treecode' } | Group-Object fixture)) {
        $candidate = @($group.Group | Sort-Object `
            @{ Expression = { if ($_.status -eq 'ok') { 0 } else { 1 } } }, `
            @{ Expression = { $v = Convert-ToDoubleOrNull $_.energy_rel_vs_grid; if ($null -eq $v) { [double]::PositiveInfinity } else { $v } } }, `
            @{ Expression = { $v = Convert-ToDoubleOrNull $_.time_ratio_vs_grid; if ($null -eq $v) { [double]::PositiveInfinity } else { $v } } }, `
            @{ Expression = { $sp = Convert-ToDoubleOrNull $_.surface_points; if ($null -eq $sp) { 0.0 } else { -1.0 * $sp } } } |
            Select-Object -First 1)
        if ($candidate.Count -gt 0) { $selected.Add($candidate[0]) }
    }
    return $selected.ToArray()
}

function New-Decision($Rows, [bool]$IsQuick) {
    $fixtures = @($Rows | Select-Object -ExpandProperty fixture -Unique)
    $candidateRows = @(Select-CandidateRows $Rows)
    $treePassRows = @($candidateRows | Where-Object {
        $_.status -eq 'ok' -and
        $_.energy_status -eq 'ok' -and
        $_.accepted_miss_status -eq 'ok'
    })
    $memoryAdvantage = @($treePassRows | Where-Object {
        $ratio = Convert-ToDoubleOrNull $_.ram_ratio_vs_grid
        $null -ne $ratio -and $ratio -le 0.70
    })
    $runtimeCompetitive = @($treePassRows | Where-Object {
        $ratio = Convert-ToDoubleOrNull $_.time_ratio_vs_grid
        $null -ne $ratio -and $ratio -le 1.25
    })
    $gridFailTabiPass = 0
    foreach ($fixture in $fixtures) {
        $grid = @($Rows | Where-Object { $_.fixture -eq $fixture -and $_.mode -eq 'grid' } | Select-Object -First 1)
        $tabi = @($treePassRows | Where-Object { $_.fixture -eq $fixture } | Select-Object -First 1)
        if ($grid.Count -gt 0 -and $tabi.Count -gt 0 -and $grid[0].status -ne 'ok') {
            $gridFailTabiPass += 1
        }
    }

    $missingAcceptedEvidence = @($candidateRows | Where-Object { $_.accepted_miss_status -eq 'unknown' }).Count
    $accuracyFailures = @($candidateRows | Where-Object {
        $_.status -ne 'ok' -or $_.energy_status -notin @('ok','grid_reference')
    }).Count
    $nonSkippedRows = @($Rows | Where-Object { $_.status -ne 'skipped' }).Count
    $promote = (-not $IsQuick) -and
        $candidateRows.Count -gt 0 -and
        $treePassRows.Count -eq $candidateRows.Count -and
        (($memoryAdvantage.Count -ge 2) -or
         ($runtimeCompetitive.Count -ge 2) -or
         ($gridFailTabiPass -ge 1))

    $decision = 'pause_new_solver_work'
    $reason = 'no_medium_memory_runtime_or_science_advantage'
    if ($IsQuick) {
        $decision = 'quick_smoke_not_authoritative'
        $reason = 'quick_mode_validates_harness_only'
    } elseif ($nonSkippedRows -eq 0) {
        $decision = 'keep_experimental'
        $reason = 'fixture_evidence_missing'
    } elseif ($promote) {
        $decision = 'promote_candidate'
        $reason = 'treecode_passed_accuracy_and_advantage_gate'
    } elseif ($accuracyFailures -gt 0 -or $missingAcceptedEvidence -gt 0) {
        $decision = 'keep_experimental'
        $reason = if ($accuracyFailures -gt 0) {
            'accuracy_or_runtime_failure_present'
        } else {
            'accepted_miss_evidence_missing_or_incomplete'
        }
    }

    [PSCustomObject]@{
        decision = $decision
        promote_candidate = [string]$promote
        quick_mode = [string]$IsQuick
        fixture_count = $fixtures.Count
        treecode_candidate_count = $candidateRows.Count
        treecode_pass_count = $treePassRows.Count
        memory_advantage_count = $memoryAdvantage.Count
        runtime_competitive_count = $runtimeCompetitive.Count
        grid_fail_tabi_pass_count = $gridFailTabiPass
        missing_accepted_miss_evidence_count = $missingAcceptedEvidence
        accuracy_failure_count = $accuracyFailures
        production_default = 'grid'
        tabi_public_status = if ($promote) { 'promotion_candidate_only' } else { 'hidden_or_experimental_opt_in' }
        decision_reason = $reason
    }
}

if (-not (Test-Path -LiteralPath $PerfMatrixScript)) {
    throw "Missing perf matrix script: $PerfMatrixScript"
}
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $TriageHelperScript)) {
    throw "Missing TABI triage helper script: $TriageHelperScript"
}
. $TriageHelperScript

New-CleanDir $WorkRoot
$matrixWork = Join-Path $WorkRoot 'matrix_run'
$fixtures = if ($Quick) { @('research-medium') } else { @('research-medium','research-large') }
$surfacePoints = @(24)
$timeout = if ($Quick -and -not $PSBoundParameters.ContainsKey('CaseTimeoutSeconds')) { 60 } else { $CaseTimeoutSeconds }

$perfParams = @{
    TitanExe = $TitanExe
    WorkRoot = $matrixWork
    Fixtures = $fixtures
    Modes = @('grid','tabi_direct','tabi_treecode')
    SurfacePointCounts = $surfacePoints
    MaxDirectSurfacePoints = 48
    TimingMode = 'detail'
    ProfileMemory = $true
    CaseTimeoutSeconds = $timeout
    AllowFailedStatus = $true
}
& $PerfMatrixScript @perfParams

$sourceMatrix = Join-Path $matrixWork 'pbe_tabi_perf_matrix.csv'
if (-not (Test-Path -LiteralPath $sourceMatrix)) {
    throw "Missing perf matrix output: $sourceMatrix"
}

$rawRows = @(Import-Csv -LiteralPath $sourceMatrix)
$gridByFixture = @{}
foreach ($grid in @($rawRows | Where-Object { $_.mode -eq 'grid' })) {
    if (-not $gridByFixture.ContainsKey($grid.fixture)) { $gridByFixture[$grid.fixture] = $grid }
}

$normalized = @($rawRows | ForEach-Object { New-NormalizedRow $_ $gridByFixture })
$matrixCsv = Join-Path $WorkRoot 'pbe_tabi_phase7_matrix.csv'
$normalized | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation -Encoding ASCII

$decision = New-Decision $normalized ([bool]$Quick)
$decisionCsv = Join-Path $WorkRoot 'pbe_tabi_phase7_decision.csv'
@($decision) | Export-Csv -LiteralPath $decisionCsv -NoTypeInformation -Encoding ASCII

$triageClassifications = (@($normalized | Group-Object triage_classification | Sort-Object Name | ForEach-Object {
    "$($_.Name)=$($_.Count)"
}) -join ';')
$summary = @(
    'phase: 7',
    'title: Benchmark-Gated PBE/TABI Direction',
    "quick_mode: $($decision.quick_mode)",
    "fixtures: $($fixtures -join ',')",
    'modes: grid,tabi_direct,tabi_treecode',
    "decision: $($decision.decision)",
    "promote_candidate: $($decision.promote_candidate)",
    "production_default: $($decision.production_default)",
    "tabi_public_status: $($decision.tabi_public_status)",
    "treecode_candidate_count: $($decision.treecode_candidate_count)",
    "treecode_pass_count: $($decision.treecode_pass_count)",
    "memory_advantage_count: $($decision.memory_advantage_count)",
    "runtime_competitive_count: $($decision.runtime_competitive_count)",
    "grid_fail_tabi_pass_count: $($decision.grid_fail_tabi_pass_count)",
    "missing_accepted_miss_evidence_count: $($decision.missing_accepted_miss_evidence_count)",
    "accuracy_failure_count: $($decision.accuracy_failure_count)",
    "decision_reason: $($decision.decision_reason)",
    "triage_classifications: $triageClassifications",
    'gate_policy: promote only if treecode passes accuracy/residual/accepted-miss evidence and memory/runtime/science advantage on medium fixtures',
    'note: quick mode validates harness only and must not change production defaults',
    "matrix_csv: $matrixCsv",
    "decision_csv: $decisionCsv"
)
$summaryPath = Join-Path $WorkRoot 'pbe_tabi_phase7_summary.dat'
Set-Content -LiteralPath $summaryPath -Value $summary -Encoding ASCII

Write-Host "PBE TABI Phase 7 gate COMPLETE: $summaryPath"
