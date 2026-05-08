param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_accuracy_matrix')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI accuracy matrix root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Invoke-Step([string]$Name, [scriptblock]$Body) {
    $start = Get-Date
    try {
        & $Body
        $elapsed = ((Get-Date) - $start).TotalSeconds
        Write-Host ("PBE TABI accuracy matrix PASS: {0} seconds={1:N3}" -f $Name, $elapsed)
        return [PSCustomObject]@{ step = $Name; status = 'ok'; seconds = ('{0:N6}' -f $elapsed); detail = 'ok' }
    } catch {
        $elapsed = ((Get-Date) - $start).TotalSeconds
        return [PSCustomObject]@{ step = $Name; status = 'failed'; seconds = ('{0:N6}' -f $elapsed); detail = $_.Exception.Message }
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
New-CleanDir $WorkRoot

$rows = @()
$rows += Invoke-Step 'kirkwood_analytic' {
    & (Join-Path $PSScriptRoot 'pbe_tabi_kirkwood_analytic.ps1') -TitanExe $TitanExe `
        -WorkRoot (Join-Path $WorkRoot 'kirkwood')
    if ($LASTEXITCODE -ne 0) { throw 'kirkwood analytic failed' }
}
$rows += Invoke-Step 'tabi_vs_grid_accuracy_subset' {
    & (Join-Path $PSScriptRoot 'pbe_tabi_vs_grid_compare.ps1') -TitanExe $TitanExe `
        -WorkRoot (Join-Path $WorkRoot 'vs_grid') -Mode accuracy -RepeatCount 1 `
        -Fixtures @('sphere','mini-protein')
    if ($LASTEXITCODE -ne 0) { throw 'tabi vs grid accuracy subset failed' }
}
$rows += Invoke-Step 'direct_vs_treecode' {
    & (Join-Path $PSScriptRoot 'pbe_tabi_treecode_accuracy.ps1') `
        -WorkRoot (Join-Path $WorkRoot 'treecode')
    if ($LASTEXITCODE -ne 0) { throw 'direct vs treecode failed' }
}
$rows += Invoke-Step 'pka_report_only' {
    & (Join-Path $PSScriptRoot 'pbe_tabi_pka_experimental_smoke.ps1') `
        -WorkRoot (Join-Path $WorkRoot 'pka')
    if ($LASTEXITCODE -ne 0) { throw 'pKa scaffold failed' }
}

$rows = @($rows | Where-Object {
    $null -ne $_ -and ($_.PSObject.Properties.Name -contains 'status') -and
    ($_.PSObject.Properties.Name -contains 'step')
})
$failed = @($rows | Where-Object { $_.status -ne 'ok' })
$outCsv = Join-Path $WorkRoot 'pbe_tabi_accuracy_matrix.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
if ($failed.Count -gt 0) {
    $messages = ($failed | ForEach-Object { "$($_.step): $($_.detail)" }) -join '; '
    throw "PBE TABI accuracy matrix failed: $messages"
}

Write-Host "PBE TABI accuracy matrix COMPLETE: $outCsv"
