param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_large_protein_triage'),
    [string[]]$Fixtures = @('barnase','research-medium','research-large'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 120,
    [int[]]$SurfacePointCounts = @(24),
    [int[]]$TreeLeafSizes = @(128),
    [double[]]$OpeningThetas = @(0.8),
    [string[]]$Threads = @('8'),
    [ValidateSet('diagonal','block','block-budgeted')]
    [string]$Preconditioner = 'block-budgeted',
    [ValidateSet('singular','buffer-aware','operator-near')]
    [string]$PreconditionerBlockNear = 'operator-near',
    [int]$PreconditionerBlockPanelCap = 32,
    [int]$PreconditionerMemoryBudgetMb = 256,
    [switch]$ProfileMemory
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
$RepoRoot = Split-Path -Parent $PSScriptRoot
$PerfMatrixScript = Join-Path $PSScriptRoot 'pbe_tabi_perf_matrix.ps1'
$TriageHelperScript = Join-Path $PSScriptRoot 'pbe_tabi_triage_helpers.ps1'
$OutCsv = Join-Path $WorkRoot 'pbe_tabi_large_protein_triage.csv'
$SummaryPath = Join-Path $WorkRoot 'pbe_tabi_large_protein_triage_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_TABI_OPERATOR_NEAR_LARGE_PROTEIN_TRIAGE_REPORT.txt'
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI large-protein triage root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Get-RowValue($Object, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Object -and $Object.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Object.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}

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

function Get-CaseRows([string]$CaseName, [hashtable]$MatrixArgs) {
    $caseRoot = Join-Path $WorkRoot $CaseName
    $MatrixArgs.WorkRoot = $caseRoot
    & $PerfMatrixScript @MatrixArgs
    $matrixCsv = Join-Path $caseRoot 'pbe_tabi_perf_matrix.csv'
    if (-not (Test-Path -LiteralPath $matrixCsv)) {
        throw "Missing TABI perf matrix output: $matrixCsv"
    }
    foreach ($row in @(Import-Csv -LiteralPath $matrixCsv)) {
        $triage = New-TabiTriageRecord -Row $row -RunDir (Get-RowValue $row 'run_dir') -ExistingClassification (Get-RowValue $row 'performance_classification')
        [PSCustomObject]@{
            case_name = $CaseName
            fixture = Get-RowValue $row 'fixture'
            mode = Get-RowValue $row 'mode'
            status = Get-RowValue $row 'status'
            triage_classification = $triage.triage_classification
            operator_near_decision = 'pending'
            total_seconds = Get-RowValue $row 'total_seconds'
            seconds_per_matvec = $triage.seconds_per_matvec
            preconditioner_setup_fraction = 'NA'
            leaf_panel_interactions = Get-RowValue $row 'leaf_panel_interactions'
            leaf_panel_interactions_per_matvec = $triage.leaf_panel_interactions_per_matvec
            leaf_vs_n2_ratio = $triage.leaf_vs_n2_ratio
            residual_delta_per_matvec = $triage.residual_delta_per_matvec
            preconditioner_requested = $triage.preconditioner_requested
            preconditioner_effective = $triage.preconditioner_effective
            preconditioner_downgrade = $triage.preconditioner_downgrade
            preconditioner_block_panel_cap = $triage.preconditioner_block_panel_cap
            preconditioner_memory_budget_bytes = $triage.preconditioner_memory_budget_bytes
            preconditioner_estimated_bytes = $triage.preconditioner_estimated_bytes
            preconditioner_valid_blocks = $triage.preconditioner_valid_blocks
            preconditioner_valid_panels = $triage.preconditioner_valid_panels
            preconditioner_ill_conditioned_blocks = $triage.preconditioner_ill_conditioned_blocks
            preconditioner_quality_rejected_blocks = $triage.preconditioner_quality_rejected_blocks
            preconditioner_quality_probe_max_gain = $triage.preconditioner_quality_probe_max_gain
            preconditioner_quality_status = $triage.preconditioner_quality_status
            preconditioner_setup_seconds = $triage.preconditioner_setup_seconds
            preconditioner_apply_seconds = $triage.preconditioner_apply_seconds
            preconditioner_apply_calls = $triage.preconditioner_apply_calls
            prepare_progress_phase = $triage.prepare_progress_phase
            prepare_progress_status = $triage.prepare_progress_status
            prepare_progress_rss_bytes = $triage.prepare_progress_rss_bytes
            panels = Get-RowValue $row 'panels'
            matvec_calls = if ((Get-RowValue $row 'matvec_calls') -ne 'NA') { Get-RowValue $row 'matvec_calls' } else { $triage.matvec_calls }
            final_residual = Get-RowValue $row 'final_residual'
            performance_classification = Get-RowValue $row 'performance_classification'
            surface_points = Get-RowValue $row 'surface_points'
            tree_leaf_size = Get-RowValue $row 'tree_leaf_size'
            opening_theta = Get-RowValue $row 'opening_theta'
            threads = Get-RowValue $row 'threads'
            run_dir = Get-RowValue $row 'run_dir'
            triage_source = $triage.triage_source
            source_matrix_csv = $matrixCsv
        }
    }
}

function Set-SetupFraction($Rows) {
    foreach ($row in $Rows) {
        $setup = Convert-ToDoubleOrNull $row.preconditioner_setup_seconds
        $total = Convert-ToDoubleOrNull $row.total_seconds
        if ($null -ne $setup -and $null -ne $total -and $total -gt 0.0) {
            $row.preconditioner_setup_fraction = Format-Number ($setup / $total)
        }
    }
}

function Set-OperatorNearDecision($Rows) {
    $diagonalByFixture = @{}
    foreach ($row in @($Rows | Where-Object { $_.case_name -eq 'diagonal' })) {
        $diagonalByFixture[$row.fixture] = $row
    }
    foreach ($row in $Rows) {
        if ($row.case_name -eq 'grid') {
            $row.operator_near_decision = 'grid_reference'
            continue
        }
        if ($row.case_name -eq 'diagonal') {
            $row.operator_near_decision = 'diagonal_reference'
            continue
        }
        $effective = [string]$row.preconditioner_effective
        $matvec = Convert-ToDoubleOrNull $row.matvec_calls
        $setup = Convert-ToDoubleOrNull $row.preconditioner_setup_seconds
        $setupFraction = Convert-ToDoubleOrNull $row.preconditioner_setup_fraction
        $secondsPerMatvec = Convert-ToDoubleOrNull $row.seconds_per_matvec
        if ($effective -ne 'block_budgeted_operator_near' -or $null -eq $matvec -or $matvec -le 0.0) {
            $row.operator_near_decision = 'preconditioner_failed'
            continue
        }
        if (($null -ne $setupFraction -and $setupFraction -gt 0.25) -or
            ($null -ne $setup -and $setup -gt 30.0)) {
            $row.operator_near_decision = 'setup_limited'
            continue
        }
        $diag = if ($diagonalByFixture.ContainsKey($row.fixture)) { $diagonalByFixture[$row.fixture] } else { $null }
        $diagMatvec = if ($null -ne $diag) { Convert-ToDoubleOrNull $diag.matvec_calls } else { $null }
        if ($row.status -eq 'ok' -and $null -ne $diagMatvec -and $matvec -lt $diagMatvec) {
            $row.operator_near_decision = 'promising'
            continue
        }
        if (($matvec -le 25.0 -or $row.triage_classification -eq 'leaf_work_super_n2') -and
            ($row.status -eq 'timeout' -or
             ($null -ne $secondsPerMatvec -and $secondsPerMatvec -gt 2.0))) {
            $row.operator_near_decision = 'leaf_kernel_limited'
            continue
        }
        $row.operator_near_decision = 'needs_review'
    }
}

function Write-TriageReport($Rows, [object[]]$FollowRows) {
    $lines = New-Object System.Collections.Generic.List[string]
    $lines.Add('TITAN TABI Operator-Near Large-Protein Triage Report')
    $lines.Add('====================================================')
    $lines.Add('')
    $lines.Add("fixtures: $($runFixtures -join ',')")
    $lines.Add("quick_mode: $([string][bool]$Quick)")
    $lines.Add("case_timeout_seconds: $timeout")
    $lines.Add("operator_near_preconditioner: $Preconditioner")
    $lines.Add("operator_near_block_near: $PreconditionerBlockNear")
    $lines.Add("block_panel_cap: $PreconditionerBlockPanelCap")
    $lines.Add("memory_budget_mb: $PreconditionerMemoryBudgetMb")
    $lines.Add('')
    $lines.Add('Main Results')
    $lines.Add('------------')
    $lines.Add('| Fixture | Case | Status | Effective | Matvecs | sec/matvec | setup_frac | quality | decision |')
    $lines.Add('| --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |')
    foreach ($row in @($Rows | Sort-Object fixture, case_name)) {
        $lines.Add(('| {0} | {1} | {2} | {3} | {4} | {5} | {6} | {7} | {8} |' -f
            $row.fixture, $row.case_name, $row.status, $row.preconditioner_effective,
            $row.matvec_calls, $row.seconds_per_matvec, $row.preconditioner_setup_fraction,
            $row.preconditioner_quality_status, $row.operator_near_decision))
    }
    if ($FollowRows.Count -gt 0) {
        $lines.Add('')
        $lines.Add('Leaf-Kernel Follow-Up')
        $lines.Add('---------------------')
        $lines.Add('| Fixture | leaf_size | theta | Status | Effective | Matvecs | sec/matvec | decision |')
        $lines.Add('| --- | ---: | ---: | --- | --- | ---: | ---: | --- |')
        foreach ($row in @($FollowRows | Sort-Object fixture, tree_leaf_size, opening_theta)) {
            $lines.Add(('| {0} | {1} | {2} | {3} | {4} | {5} | {6} | {7} |' -f
                $row.fixture, $row.tree_leaf_size, $row.opening_theta, $row.status,
                $row.preconditioner_effective, $row.matvec_calls,
                $row.seconds_per_matvec, $row.operator_near_decision))
        }
    }
    $decisionCounts = @($Rows | Group-Object operator_near_decision | Sort-Object Name | ForEach-Object {
        "$($_.Name)=$($_.Count)"
    })
    $lines.Add('')
    $lines.Add('Final Recommendation')
    $lines.Add('--------------------')
    $lines.Add("decision_counts: $($decisionCounts -join ', ')")
    if (@($Rows | Where-Object { $_.operator_near_decision -eq 'leaf_kernel_limited' }).Count -gt 0) {
        $lines.Add('recommendation: operator-near reached the leaf-kernel-limited regime on at least one fixture; prioritize leaf/tree kernel tuning for those fixtures.')
    } elseif (@($Rows | Where-Object { $_.operator_near_decision -eq 'preconditioner_failed' }).Count -gt 0) {
        $lines.Add('recommendation: operator-near still fails before useful GMRES progress on at least one fixture; keep improving preconditioner robustness before leaf tuning.')
    } elseif (@($Rows | Where-Object { $_.operator_near_decision -eq 'setup_limited' }).Count -gt 0) {
        $lines.Add('recommendation: operator-near setup overhead is the immediate bottleneck; parallelize/cache block near integration before leaf tuning.')
    } else {
        $lines.Add('recommendation: operator-near is promising on this run; expand the benchmark gate before changing production defaults.')
    }
    $lines.Add('')
    $lines.Add("triage_csv: $OutCsv")
    $lines.Add("summary_dat: $SummaryPath")
    Set-Content -LiteralPath $ReportPath -Value $lines -Encoding ASCII
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $PerfMatrixScript)) { throw "Missing TABI perf matrix script: $PerfMatrixScript" }
if (-not (Test-Path -LiteralPath $TriageHelperScript)) { throw "Missing TABI triage helper script: $TriageHelperScript" }
. $TriageHelperScript

New-CleanDir $WorkRoot

$runFixtures = if ($Quick) { @('barnase') } else { $Fixtures }
$timeout = if ($Quick -and -not $PSBoundParameters.ContainsKey('CaseTimeoutSeconds')) { 60 } else { $CaseTimeoutSeconds }
$runSurfacePoints = if ($Quick) { @(24) } else { $SurfacePointCounts }
$runLeafSizes = if ($Quick) { @(128) } else { $TreeLeafSizes }
$runThetas = if ($Quick) { @(0.8) } else { $OpeningThetas }

$baseArgs = @{
    TitanExe = $TitanExe
    Fixtures = $runFixtures
    SurfacePointCounts = $runSurfacePoints
    TreeLeafSizes = $runLeafSizes
    OpeningThetas = $runThetas
    LeafNearPolicies = @('exact')
    LeafPanelKernels = @('tile-avx2')
    AdaptiveLevels = @('strict')
    ClusterPolicies = @('off')
    Threads = $Threads
    TimingMode = 'detail'
    TabiProgressIntervalSeconds = 5
    CaseTimeoutSeconds = $timeout
    AllowFailedStatus = $true
}
if ($ProfileMemory) { $baseArgs.ProfileMemory = $true }

$rows = New-Object System.Collections.Generic.List[object]
$gridArgs = $baseArgs.Clone()
$gridArgs.Modes = @('grid')
foreach ($row in Get-CaseRows -CaseName 'grid' -MatrixArgs $gridArgs) { $rows.Add($row) }

$diagonalArgs = $baseArgs.Clone()
$diagonalArgs.Modes = @('tabi_treecode')
$diagonalArgs.Preconditioner = 'diagonal'
foreach ($row in Get-CaseRows -CaseName 'diagonal' -MatrixArgs $diagonalArgs) { $rows.Add($row) }

$operatorArgs = $baseArgs.Clone()
$operatorArgs.Modes = @('tabi_treecode')
$operatorArgs.Preconditioner = $Preconditioner
$operatorArgs.PreconditionerBlockNear = $PreconditionerBlockNear
$operatorArgs.PreconditionerBlockPanelCap = $PreconditionerBlockPanelCap
$operatorArgs.PreconditionerMemoryBudgetMb = $PreconditionerMemoryBudgetMb
foreach ($row in Get-CaseRows -CaseName 'operator_near' -MatrixArgs $operatorArgs) { $rows.Add($row) }

Set-SetupFraction $rows
Set-OperatorNearDecision $rows

$followRows = New-Object System.Collections.Generic.List[object]
$followFixtures = if ($Quick) {
    @()
} else {
    @($rows | Where-Object {
        $_.case_name -eq 'operator_near' -and
        $_.operator_near_decision -eq 'leaf_kernel_limited' -and
        ($_.fixture -eq 'barnase' -or $_.fixture -eq 'research-medium')
    } | Select-Object -ExpandProperty fixture -Unique)
}
if (@($followFixtures).Count -gt 0) {
    $followArgs = $baseArgs.Clone()
    $followArgs.Fixtures = $followFixtures
    $followArgs.Modes = @('tabi_treecode')
    $followArgs.TreeLeafSizes = @(32,64,128)
    $followArgs.OpeningThetas = @(0.45,0.8)
    $followArgs.Preconditioner = $Preconditioner
    $followArgs.PreconditionerBlockNear = $PreconditionerBlockNear
    $followArgs.PreconditionerBlockPanelCap = $PreconditionerBlockPanelCap
    $followArgs.PreconditionerMemoryBudgetMb = $PreconditionerMemoryBudgetMb
    foreach ($row in Get-CaseRows -CaseName 'leaf_followup' -MatrixArgs $followArgs) {
        $followRows.Add($row)
    }
    Set-SetupFraction $followRows
    Set-OperatorNearDecision $followRows
}

$allRows = @()
$allRows += @($rows.ToArray())
$allRows += @($followRows.ToArray())
$allRows | Export-Csv -LiteralPath $OutCsv -NoTypeInformation -Encoding ASCII

$classes = @($rows | Group-Object operator_near_decision | Sort-Object Name | ForEach-Object {
    "decision_$($_.Name): $($_.Count)"
})
$summary = @(
    'title: TABI Operator-Near Large-Protein Triage',
    "quick_mode: $([string][bool]$Quick)",
    "fixtures: $($runFixtures -join ',')",
    "case_timeout_seconds: $timeout",
    "threads: $($Threads -join ',')",
    "operator_near_preconditioner: $Preconditioner",
    "operator_near_block_near: $PreconditionerBlockNear",
    "operator_near_block_panel_cap: $PreconditionerBlockPanelCap",
    "operator_near_memory_budget_mb: $PreconditionerMemoryBudgetMb",
    'note: diagnostics are report-only and do not change production solver defaults',
    "rows: $($allRows.Count)"
) + $classes + @(
    "triage_csv: $OutCsv",
    "report_txt: $ReportPath"
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII
Write-TriageReport -Rows $rows -FollowRows $followRows

Write-Host "PBE TABI operator-near large-protein triage COMPLETE: $OutCsv"
