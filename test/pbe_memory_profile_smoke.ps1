param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_memory_profile_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

& (Join-Path $PSScriptRoot 'pbe_tabi_vs_grid_compare.ps1') -TitanExe $TitanExe `
    -WorkRoot $WorkRoot -Mode accuracy -RepeatCount 1 `
    -Fixtures @('sphere') -ProfileMemory
if ($LASTEXITCODE -ne 0) {
    throw 'pbe_tabi_vs_grid_compare memory profile run failed'
}

$csv = Join-Path $WorkRoot 'pbe_tabi_vs_grid_compare.csv'
if (-not (Test-Path -LiteralPath $csv)) {
    throw "Missing memory compare CSV: $csv"
}
$rows = @(Import-Csv -LiteralPath $csv)
if ($rows.Count -lt 2) {
    throw "Expected grid/TABI memory rows, found $($rows.Count)"
}
foreach ($row in $rows) {
    foreach ($column in @('memory_estimated_bytes_total','memory_estimated_bytes_solver',
            'memory_rss_after_bytes','memory_rss_delta_signed_bytes',
            'memory_peak_working_set_bytes')) {
        if (-not ($row.PSObject.Properties.Name -contains $column)) {
            throw "Memory compare missing column: $column"
        }
    }
    $estimated = Convert-ToDouble $row.memory_estimated_bytes_total
    if ($estimated -le 0.0) {
        throw "Non-positive estimated memory for $($row.fixture)/$($row.backend): $estimated"
    }
    if ($row.backend -eq 'tabi') {
        if ([int64]$row.panels -le 0 -or [int64]$row.memory_estimated_bytes_surface -le 0) {
            throw "TABI memory row missing surface/panel estimate for $($row.fixture)"
        }
    }
    if ($row.backend -eq 'grid' -and [int64]$row.grid_points -le 0) {
        throw "Grid memory row missing grid point count for $($row.fixture)"
    }
}

Write-Host "PBE memory profile smoke PASS: $csv"
