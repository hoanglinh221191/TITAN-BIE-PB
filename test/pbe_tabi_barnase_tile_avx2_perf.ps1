param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_barnase_tile_avx2_perf'),
    [int]$CaseTimeoutSeconds = 900
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$PerfMatrixScript = Join-Path $PSScriptRoot 'pbe_tabi_perf_matrix.ps1'
$PerfWorkRoot = Join-Path $WorkRoot 'perf_matrix'
$OutCsv = Join-Path $WorkRoot 'pbe_tabi_barnase_tile_avx2_perf.csv'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI barnase tile root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $PerfMatrixScript)) {
    throw "Missing TABI perf matrix script: $PerfMatrixScript"
}

New-CleanDir $WorkRoot

& $PerfMatrixScript `
    -TitanExe $TitanExe `
    -WorkRoot $PerfWorkRoot `
    -Fixtures @('barnase') `
    -Modes @('grid','tabi_treecode') `
    -SurfacePointCounts @(24) `
    -TreeLeafSizes @(128) `
    -OpeningThetas @(0.8) `
    -LeafNearPolicies @('exact') `
    -LeafPanelKernels @('tile-avx2') `
    -ProfileMemory `
    -CaseTimeoutSeconds $CaseTimeoutSeconds

$matrixCsv = Join-Path $PerfWorkRoot 'pbe_tabi_perf_matrix.csv'
if (-not (Test-Path -LiteralPath $matrixCsv)) {
    throw "Missing TABI barnase tile perf matrix output: $matrixCsv"
}

Copy-Item -LiteralPath $matrixCsv -Destination $OutCsv -Force
$rows = @(Import-Csv -LiteralPath $OutCsv)
if ($rows.Count -eq 0) {
    throw "Empty TABI barnase tile perf output: $OutCsv"
}

$tabiRows = @($rows | Where-Object { $_.mode -eq 'tabi_treecode' })
if ($tabiRows.Count -eq 0) {
    throw "Missing TABI barnase tile-avx2 row in $OutCsv"
}

$statusSummary = ($rows | ForEach-Object {
    "{0}:{1}:{2}" -f $_.mode, $_.status, $_.performance_classification
}) -join '; '
Write-Host "PBE TABI barnase tile-avx2 perf COMPLETE: $OutCsv"
Write-Host "PBE TABI barnase tile-avx2 statuses: $statusSummary"
