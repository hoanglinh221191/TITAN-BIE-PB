param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_barnase_timing_profile'),
    [int]$CaseTimeoutSeconds = 900,
    [string[]]$Threads = @('8'),
    [switch]$ProfileMemory
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
$PerfMatrixScript = Join-Path $PSScriptRoot 'pbe_tabi_perf_matrix.ps1'
$TriageHelperScript = Join-Path $PSScriptRoot 'pbe_tabi_triage_helpers.ps1'
$PerfWorkRoot = Join-Path $WorkRoot 'perf_matrix'
$OutCsv = Join-Path $WorkRoot 'pbe_tabi_barnase_timing_profile.csv'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI Barnase timing profile root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Get-ValueOrNA($Object, [string]$Name) {
    if ($null -eq $Object) { return 'NA' }
    if ($Object.PSObject.Properties.Name -contains $Name) {
        $value = $Object.$Name
        if ($null -ne $value -and [string]$value -ne '') { return $value }
    }
    return 'NA'
}

function Get-LastProgressRow([string]$RunDir) {
    if ([string]::IsNullOrWhiteSpace($RunDir) -or $RunDir -eq 'NA') {
        return $null
    }
    $path = Join-Path $RunDir 'pbe_tabi_longrun_progress.csv'
    if (-not (Test-Path -LiteralPath $path)) {
        return $null
    }
    $rows = @(Import-Csv -LiteralPath $path)
    if ($rows.Count -eq 0) {
        return $null
    }
    return $rows[-1]
}

function Get-ProgressClassification($PerfRow, $ProgressRow) {
    if ($PerfRow.status -ne 'timeout') {
        return $PerfRow.performance_classification
    }
    if ($null -eq $ProgressRow) {
        return 'no_progress_file'
    }
    switch -Regex ($ProgressRow.phase) {
        '^surface_' { return 'surface_limited' }
        '^tree_build_start$' { return 'tree_build_limited' }
        '^tree_build_end$' { return 'precond_limited' }
        '^precond_start$' { return 'precond_limited' }
        '^precond_end$' { return 'gmres_iteration_limited' }
        '^gmres_start$' { return 'gmres_iteration_limited' }
        '^gmres_' {
            [int64]$matvec = 0
            [void][int64]::TryParse([string]$ProgressRow.matvec_calls, [ref]$matvec)
            if ($matvec -gt 0) { return 'matvec_limited' }
            return 'gmres_iteration_limited'
        }
        default { return 'timeout_operator_limited' }
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $PerfMatrixScript)) {
    throw "Missing TABI perf matrix script: $PerfMatrixScript"
}
if (-not (Test-Path -LiteralPath $TriageHelperScript)) {
    throw "Missing TABI triage helper script: $TriageHelperScript"
}
. $TriageHelperScript

New-CleanDir $WorkRoot

$matrixArgs = @{
    TitanExe = $TitanExe
    WorkRoot = $PerfWorkRoot
    Fixtures = @('barnase')
    Modes = @('grid','tabi_treecode')
    SurfacePointCounts = @(24)
    TreeLeafSizes = @(128)
    OpeningThetas = @(0.8)
    LeafNearPolicies = @('exact')
    LeafPanelKernels = @('tile-avx2')
    AdaptiveLevels = @('strict')
    ClusterPolicies = @('off')
    Threads = $Threads
    Preconditioner = 'block'
    TimingMode = 'detail'
    TabiProgressIntervalSeconds = 5
    CaseTimeoutSeconds = $CaseTimeoutSeconds
}
if ($ProfileMemory) {
    $matrixArgs.ProfileMemory = $true
}

& $PerfMatrixScript @matrixArgs

$matrixCsv = Join-Path $PerfWorkRoot 'pbe_tabi_perf_matrix.csv'
if (-not (Test-Path -LiteralPath $matrixCsv)) {
    throw "Missing TABI Barnase perf matrix output: $matrixCsv"
}
$rows = @(Import-Csv -LiteralPath $matrixCsv)
if ($rows.Count -lt 2) {
    throw "Expected grid and TABI Barnase rows, found $($rows.Count)"
}

$summaryRows = foreach ($row in $rows) {
    $progress = if ($row.mode -eq 'tabi_treecode') {
        Get-LastProgressRow -RunDir $row.run_dir
    } else {
        $null
    }
    $classification = Get-ProgressClassification -PerfRow $row -ProgressRow $progress
    $triage = New-TabiTriageRecord -Row $row -RunDir $row.run_dir -ExistingClassification $classification
    [PSCustomObject]@{
        fixture = $row.fixture
        mode = $row.mode
        status = $row.status
        classification = $triage.triage_classification
        progress_classification = $classification
        triage_classification = $triage.triage_classification
        total_seconds = $row.total_seconds
        energy_kJmol = $row.energy_kJmol
        peak_ram_bytes = $row.peak_ram_bytes
        run_dir = $row.run_dir
        progress_phase = Get-ValueOrNA $progress 'phase'
        progress_elapsed_seconds = Get-ValueOrNA $progress 'elapsed_seconds'
        gmres_iter = Get-ValueOrNA $progress 'gmres_iter'
        restart_iter = Get-ValueOrNA $progress 'restart_iter'
        matvec_calls = Get-ValueOrNA $progress 'matvec_calls'
        residual = Get-ValueOrNA $progress 'residual'
        panels = Get-ValueOrNA $progress 'panels'
        unknown_count = Get-ValueOrNA $progress 'unknown_count'
        tree_nodes = Get-ValueOrNA $progress 'tree_nodes'
        leaf_panel_interactions = Get-ValueOrNA $progress 'leaf_panel_interactions'
        leaf_tile_panel_interactions = Get-ValueOrNA $progress 'leaf_tile_panel_interactions'
        leaf_residual_panel_interactions = Get-ValueOrNA $progress 'leaf_residual_panel_interactions'
        far_node_interactions = Get-ValueOrNA $progress 'far_node_interactions'
        singular_near_interactions = Get-ValueOrNA $progress 'singular_near_interactions'
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
        threads = if ($row.threads) { $row.threads } else { Get-ValueOrNA $progress 'threads' }
        rss_bytes = Get-ValueOrNA $progress 'rss_bytes'
        source_matrix_csv = $matrixCsv
    }
}

$summaryRows | Export-Csv -LiteralPath $OutCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE TABI Barnase timing profile COMPLETE: $OutCsv"
