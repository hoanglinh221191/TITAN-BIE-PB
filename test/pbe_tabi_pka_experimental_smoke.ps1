param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_pka_experimental'),
    [string]$Subset = '1STN'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI pKa root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

New-CleanDir $WorkRoot

$subsetRows = @(
    [PSCustomObject]@{
        pdb_id = '1STN'
        residue = 'ASP/Glu/Lys/His subset'
        experimental_pka = 'PKAD/literature'
        computed_shift_grid = 'NA'
        computed_shift_tabi = 'NA'
        diff_grid = 'NA'
        diff_tabi = 'NA'
        status = 'report_only_no_titration_workflow'
        source = 'PKAD 1350 pKa in 157 wild-type proteins; RCSB 1STN structure'
    }
)

if ($Subset -ne '1STN') {
    $subsetRows += [PSCustomObject]@{
        pdb_id = $Subset
        residue = 'candidate'
        experimental_pka = 'NA'
        computed_shift_grid = 'NA'
        computed_shift_tabi = 'NA'
        diff_grid = 'NA'
        diff_tabi = 'NA'
        status = 'candidate_not_ingested'
        source = 'add PKAD residue rows before enabling numeric pKa comparison'
    }
}

$outCsv = Join-Path $WorkRoot 'pbe_tabi_pka_compare.csv'
$subsetRows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII

$metadata = Join-Path $WorkRoot 'pbe_tabi_pka_sources.txt'
Set-Content -LiteralPath $metadata -Encoding ASCII -Value @(
    'TABI pKa experimental scaffold is report-only in this phase.',
    'Primary experimental database: https://academic.oup.com/database/article/doi/10.1093/database/baz024/5359213',
    'Default structural scaffold: https://www.rcsb.org/structure/1stn',
    'No production pKa claim is made until TITAN has a titration-state workflow.'
)

Write-Host "PBE TABI pKa experimental scaffold PASS: $outCsv"
