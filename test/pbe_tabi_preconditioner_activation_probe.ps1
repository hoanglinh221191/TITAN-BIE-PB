param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_preconditioner_activation_probe'),
    [string[]]$Fixtures = @('barnase','research-medium','research-large'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 30,
    [int[]]$SurfacePointCounts = @(24),
    [int[]]$TreeLeafSizes = @(128),
    [double[]]$OpeningThetas = @(0.8),
    [int[]]$PreconditionerMemoryBudgetMb = @(128,256,512),
    [int[]]$PreconditionerBlockPanelCap = @(16,24,32),
    [string[]]$Threads = @('8'),
    [switch]$ProfileMemory
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
$PerfMatrixScript = Join-Path $PSScriptRoot 'pbe_tabi_perf_matrix.ps1'
$TriageHelperScript = Join-Path $PSScriptRoot 'pbe_tabi_triage_helpers.ps1'
$OutCsv = Join-Path $WorkRoot 'pbe_tabi_preconditioner_activation_probe.csv'
$SummaryPath = Join-Path $WorkRoot 'pbe_tabi_preconditioner_activation_probe_summary.dat'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI preconditioner activation root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
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
            [System.Globalization.CultureInfo]::InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $null
}

function Get-ActivationClassification($Row) {
    $mode = Get-RowValue $Row 'preconditioner_mode'
    if ((Get-RowValue $Row 'mode') -eq 'grid') { return 'grid_reference' }
    if ($mode -eq 'diagonal') { return 'diagonal_only' }
    $effective = (Get-RowValue $Row 'preconditioner_effective').ToLowerInvariant()
    if ($effective -eq 'diagonal_block_budget_cap' -or
        $effective -eq 'diagonal_coarse_budget_cap') {
        return 'budget_capped'
    }
    if ($effective -eq 'diagonal_block_memory_cap') {
        return 'memory_capped'
    }
    if ($effective -eq 'diagonal_block_ill_conditioned' -or
        $effective -eq 'diagonal_coarse_factor_fallback') {
        return 'ill_conditioned'
    }
    if ($effective -like 'block*' -or $effective -like 'coarse*') {
        $matvec = Convert-ToDoubleOrNull (Get-RowValue $Row 'matvec_calls')
        $residualDelta = Convert-ToDoubleOrNull (Get-RowValue $Row 'residual_delta_per_matvec')
        if (($null -ne $matvec -and $matvec -gt 0.0) -or
            ($null -ne $residualDelta -and [Math]::Abs($residualDelta) -gt 0.0)) {
            return 'active_with_progress'
        }
        return 'active_no_progress'
    }
    return 'diagonal_only'
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $PerfMatrixScript)) { throw "Missing perf matrix script: $PerfMatrixScript" }
if (-not (Test-Path -LiteralPath $TriageHelperScript)) { throw "Missing triage helper script: $TriageHelperScript" }
. $TriageHelperScript

New-CleanDir $WorkRoot

$runFixtures = if ($Quick) { @('barnase') } else { $Fixtures }
$runBudgets = if ($Quick) { @(256) } else { $PreconditionerMemoryBudgetMb }
$runCaps = if ($Quick) { @(32) } else { $PreconditionerBlockPanelCap }

$cases = New-Object System.Collections.Generic.List[object]
$cases.Add([PSCustomObject]@{ Name='grid'; Preconditioner='NA'; BlockNear=''; Cap=0; Budget=0; Modes=@('grid') })
$cases.Add([PSCustomObject]@{ Name='diagonal'; Preconditioner='diagonal'; BlockNear=''; Cap=0; Budget=0; Modes=@('tabi_treecode') })
$cases.Add([PSCustomObject]@{ Name='block'; Preconditioner='block'; BlockNear=''; Cap=0; Budget=0; Modes=@('tabi_treecode') })
foreach ($budget in $runBudgets) {
    foreach ($cap in $runCaps) {
        $cases.Add([PSCustomObject]@{
            Name = "block-budgeted_${budget}mb_cap${cap}"
            Preconditioner = 'block-budgeted'
            BlockNear = 'buffer-aware'
            Cap = $cap
            Budget = $budget
            Modes = @('tabi_treecode')
        })
    }
}

$allRows = New-Object System.Collections.Generic.List[object]
foreach ($case in $cases) {
    $caseRoot = Join-Path $WorkRoot $case.Name
    $matrixArgs = @{
        TitanExe = $TitanExe
        WorkRoot = $caseRoot
        Fixtures = $runFixtures
        Modes = $case.Modes
        SurfacePointCounts = $SurfacePointCounts
        TreeLeafSizes = $TreeLeafSizes
        OpeningThetas = $OpeningThetas
        LeafNearPolicies = @('exact')
        LeafPanelKernels = @('tile-avx2')
        AdaptiveLevels = @('strict')
        ClusterPolicies = @('off')
        Threads = $Threads
        TimingMode = 'detail'
        TabiProgressIntervalSeconds = 5
        CaseTimeoutSeconds = $CaseTimeoutSeconds
        AllowFailedStatus = $true
    }
    if ($case.Preconditioner -ne 'NA') { $matrixArgs.Preconditioner = $case.Preconditioner }
    if (-not [string]::IsNullOrWhiteSpace($case.BlockNear)) { $matrixArgs.PreconditionerBlockNear = $case.BlockNear }
    if ($case.Cap -gt 0) { $matrixArgs.PreconditionerBlockPanelCap = $case.Cap }
    if ($case.Budget -gt 0) { $matrixArgs.PreconditionerMemoryBudgetMb = $case.Budget }
    if ($ProfileMemory) { $matrixArgs.ProfileMemory = $true }
    & $PerfMatrixScript @matrixArgs

    $matrixCsv = Join-Path $caseRoot 'pbe_tabi_perf_matrix.csv'
    if (-not (Test-Path -LiteralPath $matrixCsv)) { throw "Missing matrix CSV: $matrixCsv" }
    foreach ($row in @(Import-Csv -LiteralPath $matrixCsv)) {
        $triage = New-TabiTriageRecord -Row $row -RunDir (Get-RowValue $row 'run_dir') -ExistingClassification (Get-RowValue $row 'performance_classification')
        $out = [PSCustomObject]@{
            fixture = Get-RowValue $row 'fixture'
            mode = Get-RowValue $row 'mode'
            status = Get-RowValue $row 'status'
            preconditioner_mode = $case.Preconditioner
            preconditioner_block_near = if ($case.BlockNear) { $case.BlockNear } else { Get-RowValue $row 'preconditioner_block_near' }
            preconditioner_block_panel_cap = $triage.preconditioner_block_panel_cap
            preconditioner_memory_budget_bytes = $triage.preconditioner_memory_budget_bytes
            preconditioner_estimated_bytes = $triage.preconditioner_estimated_bytes
            preconditioner_valid_blocks = $triage.preconditioner_valid_blocks
            preconditioner_valid_panels = $triage.preconditioner_valid_panels
            preconditioner_ill_conditioned_blocks = $triage.preconditioner_ill_conditioned_blocks
            preconditioner_setup_seconds = $triage.preconditioner_setup_seconds
            preconditioner_apply_seconds = $triage.preconditioner_apply_seconds
            preconditioner_apply_calls = $triage.preconditioner_apply_calls
            preconditioner_effective = $triage.preconditioner_effective
            activation_classification = 'pending'
            triage_classification = $triage.triage_classification
            total_seconds = Get-RowValue $row 'total_seconds'
            seconds_per_matvec = $triage.seconds_per_matvec
            leaf_panel_interactions_per_matvec = $triage.leaf_panel_interactions_per_matvec
            leaf_vs_n2_ratio = $triage.leaf_vs_n2_ratio
            residual_delta_per_matvec = $triage.residual_delta_per_matvec
            matvec_calls = Get-RowValue $row 'matvec_calls'
            final_residual = Get-RowValue $row 'final_residual'
            prepare_progress_phase = $triage.prepare_progress_phase
            triage_source = $triage.triage_source
            surface_points = Get-RowValue $row 'surface_points'
            tree_leaf_size = Get-RowValue $row 'tree_leaf_size'
            opening_theta = Get-RowValue $row 'opening_theta'
            threads = Get-RowValue $row 'threads'
            run_dir = Get-RowValue $row 'run_dir'
            source_matrix_csv = $matrixCsv
        }
        $out.activation_classification = Get-ActivationClassification $out
        $allRows.Add($out)
    }
}

$allRows | Export-Csv -LiteralPath $OutCsv -NoTypeInformation -Encoding ASCII

$classes = @($allRows | Group-Object activation_classification | Sort-Object Name | ForEach-Object {
    "activation_$($_.Name): $($_.Count)"
})
$summary = @(
    'title: TABI Preconditioner Activation Probe',
    "quick_mode: $([string][bool]$Quick)",
    "fixtures: $($runFixtures -join ',')",
    "case_timeout_seconds: $CaseTimeoutSeconds",
    "threads: $($Threads -join ',')",
    'note: diagnostics are report-only and do not change solver defaults',
    "rows: $($allRows.Count)"
) + $classes + @(
    "activation_csv: $OutCsv"
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

Write-Host "PBE TABI preconditioner activation probe COMPLETE: $OutCsv"
