param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_reference_cache')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI reference-cache root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'
$pqr1a63 = Join-Path $InputRoot '1a63.pqr'
$status = if (Test-Path -LiteralPath $pqr1a63) { 'ready' } else { 'missing_local_pqr' }
$detail = if ($status -eq 'ready') {
    'run pbe_tabi_vs_grid_compare.ps1 -Fixtures 1a63'
} else {
    'place a clean-room generated 1A63 PQR at test/pbe_inputs_local/1a63.pqr before enabling numeric 1A63 gates'
}

New-CleanDir $WorkRoot
$rows = @(
    [PSCustomObject]@{
        pdb_id = '1A63'
        role = 'TABI/APBS computational protein reference'
        pqr_path = $pqr1a63
        status = $status
        source = 'TABI JCP 2013 and TABI-PB 2.0 literature reference; no runtime dependency'
        detail = $detail
    },
    [PSCustomObject]@{
        pdb_id = '1STN'
        role = 'experimental pKa scaffold candidate'
        pqr_path = 'NA'
        status = 'report_only'
        source = 'RCSB 1STN plus PKAD/literature pKa rows'
        detail = 'handled by pbe_tabi_pka_experimental_smoke.ps1 until titration-state workflow exists'
    }
)
$outCsv = Join-Path $WorkRoot 'pbe_tabi_reference_cache.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE TABI reference cache report: $outCsv"
