param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_1a63_paper_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'
$fixture = Join-Path $InputRoot '1a63.pqr'
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$outCsv = Join-Path $WorkRoot 'pbe_tabi_1a63_paper_smoke.csv'

if (-not (Test-Path -LiteralPath $fixture)) {
    [PSCustomObject]@{
        fixture = '1A63'
        status = 'missing_fixture'
        detail = "Missing $fixture; paper smoke is report-only until the local PQR/mesh fixture is cached."
    } | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
    Write-Host "PBE TABI 1A63 paper smoke SKIP: missing fixture ($fixture)"
    exit 0
}

& (Join-Path $PSScriptRoot 'pbe_tabi_vs_grid_compare.ps1') -TitanExe $TitanExe `
    -WorkRoot $WorkRoot -Mode accuracy -RepeatCount 1 -Fixtures @('1a63') -ProfileMemory
if ($LASTEXITCODE -ne 0) { throw '1A63 TABI-vs-grid compare failed' }

$rows = @(Import-Csv -LiteralPath (Join-Path $WorkRoot 'pbe_tabi_vs_grid_compare.csv'))
$tabi = @($rows | Where-Object { $_.backend -eq 'tabi' })[0]
if ($null -eq $tabi) { throw 'Missing TABI 1A63 row.' }
if ($tabi.status -ne 'ok' -or $tabi.reaction_energy_source -ne 'tabi_boundary_integral') {
    throw "Unexpected 1A63 TABI status/source: status=$($tabi.status) source=$($tabi.reaction_energy_source)"
}
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE TABI 1A63 paper smoke PASS: $outCsv"
