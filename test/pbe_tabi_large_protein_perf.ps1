param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_large_protein_perf'),
    [string[]]$Fixtures = @('charged-peptide','compact-mini-protein'),
    [string[]]$Threads = @('1'),
    [int]$CaseTimeoutSeconds = 0,
    [switch]$ProfileMemory
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$PerfMatrixScript = Join-Path $PSScriptRoot 'pbe_tabi_perf_matrix.ps1'
$PerfWorkRoot = Join-Path $WorkRoot 'perf_matrix'
$OutCsv = Join-Path $WorkRoot 'pbe_tabi_large_protein_perf.csv'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI large-protein perf root: $full"
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

$normalizedFixtures = @($Fixtures | ForEach-Object { ([string]$_).Trim().ToLowerInvariant() } | Where-Object { $_ })
$leafNearPolicies =
    if ($normalizedFixtures.Count -eq 1 -and $normalizedFixtures[0] -eq 'barnase') {
        @('exact')
    } else {
        @('exact','large-protein','adaptive-large')
    }
$adaptiveLevels =
    if ($normalizedFixtures.Count -eq 1 -and $normalizedFixtures[0] -eq 'barnase') {
        @('strict')
    } else {
        @('strict','balanced','aggressive')
    }

$matrixArgs = @{
    TitanExe = $TitanExe
    WorkRoot = $PerfWorkRoot
    Fixtures = $Fixtures
    Modes = @('grid','tabi_treecode')
    SurfacePointCounts = @(24)
    TreeLeafSizes = @(128)
    OpeningThetas = @(0.8)
    LeafNearPolicies = $leafNearPolicies
    LeafPanelKernels = @('tile-avx2')
    AdaptiveLevels = $adaptiveLevels
    ClusterPolicies = @('off')
    Threads = $Threads
    Preconditioner = 'block'
}
if ($ProfileMemory) {
    $matrixArgs.ProfileMemory = $true
}
if ($CaseTimeoutSeconds -gt 0) {
    $matrixArgs.CaseTimeoutSeconds = $CaseTimeoutSeconds
}

& $PerfMatrixScript @matrixArgs

$matrixCsv = Join-Path $PerfWorkRoot 'pbe_tabi_perf_matrix.csv'
if (-not (Test-Path -LiteralPath $matrixCsv)) {
    throw "Missing TABI large-protein perf matrix output: $matrixCsv"
}

Copy-Item -LiteralPath $matrixCsv -Destination $OutCsv -Force
$rows = @(Import-Csv -LiteralPath $OutCsv)
if ($rows.Count -eq 0) {
    throw "Empty TABI large-protein perf output: $OutCsv"
}

$tabiRows = @($rows | Where-Object { $_.mode -eq 'tabi_treecode' })
if ($tabiRows.Count -eq 0) {
    throw "Missing TABI large-protein rows in $OutCsv"
}

$statusSummary = ($rows | ForEach-Object {
    "{0}:{1}:{2}:{3}" -f $_.fixture, $_.mode, $_.status, $_.performance_classification
}) -join '; '
Write-Host "PBE TABI large-protein perf COMPLETE: $OutCsv"
Write-Host "PBE TABI large-protein statuses: $statusSummary"
