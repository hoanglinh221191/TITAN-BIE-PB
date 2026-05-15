param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_barnase_surface_budget_probe'),
    [int[]]$SurfacePointCounts = @(8,12,16,24),
    [int[]]$TreeLeafSizes = @(32,64,128),
    [double[]]$OpeningThetas = @(0.5,0.8,1.0),
    [string[]]$Preconditioners = @('block'),
    [string[]]$Threads = @('8'),
    [int]$CaseTimeoutSeconds = 300,
    [switch]$ProfileMemory
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$PerfMatrixScript = Join-Path $PSScriptRoot 'pbe_tabi_perf_matrix.ps1'
$OutCsv = Join-Path $WorkRoot 'pbe_tabi_barnase_surface_budget_probe.csv'
$BaselineLeafPanelPerMatvec = 843545779.0

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI Barnase budget probe root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Convert-ToDoubleOrNull($Value) {
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

function Convert-ToIntOrNull($Value) {
    if ($null -eq $Value) { return $null }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $null }
    $parsed = [int64]0
    if ([int64]::TryParse($text, [ref]$parsed)) { return $parsed }
    return $null
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
    if ($rows.Count -eq 0) { return $null }
    return $rows[-1]
}

function Classify-Budget($PerfRow, $ProgressRow) {
    if ($PerfRow.mode -eq 'grid') { return 'grid_reference' }
    $budgetStatus = Get-ValueOrNA $PerfRow 'panel_budget_status'
    if ($budgetStatus -eq 'target_density') { return 'panel_budget_candidate' }
    $panels = Convert-ToIntOrNull (Get-ValueOrNA $ProgressRow 'panels')
    if ($null -eq $panels) { $panels = Convert-ToIntOrNull $PerfRow.panels }
    if ($null -eq $panels) { return 'no_panel_profile' }
    if ($panels -lt 8000) { return 'surface_too_coarse' }
    if ($panels -ge 10000 -and $panels -le 15000) { return 'panel_budget_candidate' }
    if ($panels -gt 15000 -and $panels -le 15500) { return 'panel_budget_near_target' }
    if ($panels -gt 15000) { return 'surface_still_dense' }
    return 'panel_budget_low'
}

function Classify-Work($PerfRow, $ProgressRow) {
    if ($PerfRow.mode -eq 'grid') { return 'grid_reference' }
    $leaf = Convert-ToDoubleOrNull (Get-ValueOrNA $ProgressRow 'leaf_panel_interactions')
    $matvec = Convert-ToDoubleOrNull (Get-ValueOrNA $ProgressRow 'matvec_calls')
    if ($null -eq $leaf) { $leaf = Convert-ToDoubleOrNull $PerfRow.leaf_panel_interactions }
    if ($null -eq $matvec) { $matvec = Convert-ToDoubleOrNull $PerfRow.matvec_calls }
    if ($null -eq $leaf -or $null -eq $matvec -or $matvec -le 0.0) {
        return 'no_work_profile'
    }
    $perMatvec = $leaf / $matvec
    if ($perMatvec -le ($BaselineLeafPanelPerMatvec / 5.0)) {
        return 'leaf_work_reduced_5x'
    }
    if ($perMatvec -lt $BaselineLeafPanelPerMatvec) {
        return 'leaf_work_reduced_partial'
    }
    return 'leaf_work_not_reduced'
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $PerfMatrixScript)) {
    throw "Missing TABI perf matrix script: $PerfMatrixScript"
}

New-CleanDir $WorkRoot

$allRows = @()
foreach ($precond in $Preconditioners) {
    $perfRoot = Join-Path $WorkRoot ("perf_{0}" -f ($precond -replace '[^A-Za-z0-9_-]', '_'))
    $matrixArgs = @{
        TitanExe = $TitanExe
        WorkRoot = $perfRoot
        Fixtures = @('barnase')
        Modes = @('grid','tabi_treecode')
        SurfacePointCounts = $SurfacePointCounts
        TreeLeafSizes = $TreeLeafSizes
        OpeningThetas = $OpeningThetas
        LeafNearPolicies = @('exact')
        LeafPanelKernels = @('tile-avx2')
        AdaptiveLevels = @('strict')
        ClusterPolicies = @('off')
        Threads = $Threads
        Preconditioner = $precond
        TimingMode = 'detail'
        TabiProgressIntervalSeconds = 5
        CaseTimeoutSeconds = $CaseTimeoutSeconds
        AllowFailedStatus = $true
    }
    if ($ProfileMemory) { $matrixArgs.ProfileMemory = $true }
    & $PerfMatrixScript @matrixArgs
    $matrixCsv = Join-Path $perfRoot 'pbe_tabi_perf_matrix.csv'
    if (-not (Test-Path -LiteralPath $matrixCsv)) {
        throw "Missing perf matrix output for preconditioner=$precond"
    }
    $rows = @(Import-Csv -LiteralPath $matrixCsv)
    foreach ($row in $rows) {
        $progress = if ($row.mode -eq 'tabi_treecode') {
            Get-LastProgressRow -RunDir $row.run_dir
        } else {
            $null
        }
        $matvec = Convert-ToDoubleOrNull (Get-ValueOrNA $progress 'matvec_calls')
        $leaf = Convert-ToDoubleOrNull (Get-ValueOrNA $progress 'leaf_panel_interactions')
        $leafPerMatvec =
            if ($null -ne $matvec -and $matvec -gt 0.0 -and $null -ne $leaf) {
                $leaf / $matvec
            } else {
                'NA'
            }
        $panels = Get-ValueOrNA $progress 'panels'
        if ($panels -eq 'NA') { $panels = $row.panels }
        $allRows += [PSCustomObject]@{
            fixture = $row.fixture
            mode = $row.mode
            status = $row.status
            preconditioner = $precond
            surface_points = $row.surface_points
            tree_leaf_size = $row.tree_leaf_size
            opening_theta = $row.opening_theta
            threads = $row.threads
            panels = $panels
            unknown_count = Get-ValueOrNA $progress 'unknown_count'
            panel_density_per_A2 = $row.panel_density_per_A2
            panel_budget_status = $row.panel_budget_status
            total_seconds = $row.total_seconds
            energy_kJmol = $row.energy_kJmol
            peak_ram_bytes = $row.peak_ram_bytes
            progress_phase = Get-ValueOrNA $progress 'phase'
            progress_elapsed_seconds = Get-ValueOrNA $progress 'elapsed_seconds'
            gmres_iter = Get-ValueOrNA $progress 'gmres_iter'
            matvec_calls = Get-ValueOrNA $progress 'matvec_calls'
            residual = Get-ValueOrNA $progress 'residual'
            leaf_panel_interactions = Get-ValueOrNA $progress 'leaf_panel_interactions'
            leaf_panel_per_matvec = if ($leafPerMatvec -eq 'NA') { 'NA' } else { $leafPerMatvec.ToString('G17', $InvariantCulture) }
            leaf_tile_panel_interactions = Get-ValueOrNA $progress 'leaf_tile_panel_interactions'
            leaf_residual_panel_interactions = Get-ValueOrNA $progress 'leaf_residual_panel_interactions'
            far_node_interactions = Get-ValueOrNA $progress 'far_node_interactions'
            singular_near_interactions = Get-ValueOrNA $progress 'singular_near_interactions'
            rss_bytes = Get-ValueOrNA $progress 'rss_bytes'
            budget_classification = Classify-Budget -PerfRow $row -ProgressRow $progress
            work_classification = Classify-Work -PerfRow $row -ProgressRow $progress
            perf_classification = $row.performance_classification
            run_dir = $row.run_dir
        }
    }
}

$allRows | Export-Csv -LiteralPath $OutCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE TABI Barnase surface-budget probe COMPLETE: $OutCsv"
