param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_vs_grid_barnase_compare'),
    [int]$RepeatCount = 1,
    [switch]$FullBarnase
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$fixtures = if ($FullBarnase) {
    @('barnase_linear_salt0','barnase_linear_salt015','barnase_linear_salt2')
} else {
    @('barnase_linear_salt0')
}

& (Join-Path $PSScriptRoot 'pbe_tabi_vs_grid_compare.ps1') `
    -TitanExe $TitanExe `
    -WorkRoot $WorkRoot `
    -Fixtures $fixtures `
    -RepeatCount $RepeatCount `
    -Mode accuracy `
    -ProfileMemory
if ($LASTEXITCODE -ne 0) {
    throw 'Barnase TABI-vs-grid compare failed.'
}

$csv = Join-Path $WorkRoot 'pbe_tabi_vs_grid_compare.csv'
if (-not (Test-Path -LiteralPath $csv)) {
    throw "Missing barnase TABI-vs-grid compare CSV: $csv"
}
$rows = @(Import-Csv -LiteralPath $csv)
foreach ($fixture in $fixtures) {
    $grid = @($rows | Where-Object { $_.fixture -eq $fixture -and $_.backend -eq 'grid' })
    $tabi = @($rows | Where-Object { $_.fixture -eq $fixture -and $_.backend -eq 'tabi' })
    if ($grid.Count -eq 0 -or $tabi.Count -eq 0) {
        throw "Missing grid/TABI rows for $fixture"
    }
    foreach ($row in @($grid[0], $tabi[0])) {
        foreach ($column in @('energy_kJmol','total_seconds',
                'memory_estimated_bytes_total','memory_peak_working_set_bytes')) {
            if (-not ($row.PSObject.Properties.Name -contains $column) -or
                [string]::IsNullOrWhiteSpace([string]$row.$column) -or
                [string]$row.$column -eq 'NA') {
                throw "Missing $column for $fixture/$($row.backend)"
            }
        }
    }
}

Write-Host "PBE TABI-vs-grid barnase compare PASS: $csv"
