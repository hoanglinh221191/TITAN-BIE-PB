param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_timeout_progress_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$PerfMatrixScript = Join-Path $PSScriptRoot 'pbe_tabi_perf_matrix.ps1'
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $PerfMatrixScript)) {
    throw "Missing TABI perf matrix script: $PerfMatrixScript"
}

& $PerfMatrixScript -TitanExe $TitanExe -WorkRoot $WorkRoot `
    -Fixtures @('compact-mini-protein') `
    -Modes @('tabi_treecode') `
    -SurfacePointCounts @(24) `
    -TreeLeafSizes @(128) `
    -OpeningThetas @(0.8) `
    -LeafNearPolicies @('exact') `
    -LeafPanelKernels @('tile-avx2') `
    -AdaptiveLevels @('strict') `
    -ClusterPolicies @('off') `
    -Threads @('1') `
    -Preconditioner block `
    -TimingMode detail `
    -TabiProgressIntervalSeconds 1 `
    -CaseTimeoutSeconds 1

$matrixCsv = Join-Path $WorkRoot 'pbe_tabi_perf_matrix.csv'
if (-not (Test-Path -LiteralPath $matrixCsv)) {
    throw "Missing timeout progress matrix CSV: $matrixCsv"
}
$rows = @(Import-Csv -LiteralPath $matrixCsv)
if ($rows.Count -ne 1) {
    throw "Expected one timeout progress row, found $($rows.Count)"
}
$row = $rows[0]
if ($row.status -ne 'timeout') {
    throw "Expected timeout row, found status=$($row.status)"
}
$progressPath = Join-Path $row.run_dir 'pbe_tabi_longrun_progress.csv'
if (-not (Test-Path -LiteralPath $progressPath)) {
    throw "Timeout row did not preserve TABI progress CSV: $progressPath"
}
$progress = @(Import-Csv -LiteralPath $progressPath)
if ($progress.Count -lt 1) {
    throw "Timeout progress CSV is empty: $progressPath"
}
if ($progress[0].phase -ne 'surface_start') {
    throw "Unexpected first timeout progress phase: $($progress[0].phase)"
}

Write-Host "PBE TABI timeout progress smoke PASS: $progressPath"
