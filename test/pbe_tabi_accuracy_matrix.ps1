param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_accuracy_matrix'),
    [string[]]$SphereMeshModes = @('latlong')
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
$rows += Invoke-Step 'near_integration_smoke' {
    & (Join-Path $PSScriptRoot 'pbe_tabi_near_integration_smoke.ps1') -TitanExe $TitanExe `
        -WorkRoot (Join-Path $WorkRoot 'near_integration')
    if ($LASTEXITCODE -ne 0) { throw 'near integration smoke failed' }
}
$rows += Invoke-Step 'near_kernel_blocks' {
    & (Join-Path $PSScriptRoot 'pbe_tabi_near_kernel_blocks_smoke.ps1') `
        -WorkRoot (Join-Path $WorkRoot 'near_kernel_blocks')
    if ($LASTEXITCODE -ne 0) { throw 'near kernel block smoke failed' }
}
$rows += Invoke-Step 'near_block_rank' {
    & (Join-Path $PSScriptRoot 'pbe_tabi_near_block_rank.ps1') -TitanExe $TitanExe `
        -WorkRoot (Join-Path $WorkRoot 'near_block_rank')
    if ($LASTEXITCODE -ne 0) { throw 'near block rank failed' }
}
$rows += Invoke-Step 'asi_self_one_ring' {
    & (Join-Path $PSScriptRoot 'pbe_tabi_asi_self_one_ring_smoke.ps1') -TitanExe $TitanExe `
        -WorkRoot (Join-Path $WorkRoot 'asi_self_one_ring')
    if ($LASTEXITCODE -ne 0) { throw 'ASI self/one-ring smoke failed' }
}
$rows += Invoke-Step 'born_ion_matrix' {
    & (Join-Path $PSScriptRoot 'pbe_born_ion_matrix.ps1') -TitanExe $TitanExe `
        -WorkRoot (Join-Path $WorkRoot 'born_ion') -RadiiA @(3.0) `
        -Charges @(1.0,-1.0) -SaltM @(0.0,0.15) `
        -SurfacePointCounts @(96,192) -SphereMeshModes $SphereMeshModes `
        -GridScales @(0.60) -SkipGrid
    if ($LASTEXITCODE -ne 0) { throw 'Born ion matrix failed' }
}
$rows += Invoke-Step 'kirkwood_oracle_selfcheck' {
    & (Join-Path $PSScriptRoot 'pbe_kirkwood_oracle_selfcheck.ps1') `
        -WorkRoot (Join-Path $WorkRoot 'kirkwood_oracle_selfcheck')
    if ($LASTEXITCODE -ne 0) { throw 'Kirkwood oracle selfcheck failed' }
}
$rows += Invoke-Step 'kirkwood_offcenter_matrix' {
    & (Join-Path $PSScriptRoot 'pbe_kirkwood_offcenter_matrix.ps1') -TitanExe $TitanExe `
        -WorkRoot (Join-Path $WorkRoot 'kirkwood_offcenter') -RadiiA @(5.0) `
        -Charges @(1.0,-1.0) -OffsetFractions @(0.0,0.35) `
        -SaltM @(0.0,0.15) -SurfacePointCounts @(96,384) `
        -SphereMeshModes $SphereMeshModes
    if ($LASTEXITCODE -ne 0) { throw 'Kirkwood off-center matrix failed' }
}
$rows += Invoke-Step 'kirkwood_multicharge_matrix' {
    & (Join-Path $PSScriptRoot 'pbe_kirkwood_multicharge_matrix.ps1') -TitanExe $TitanExe `
        -WorkRoot (Join-Path $WorkRoot 'kirkwood_multicharge') `
        -SaltM @(0.0,0.15) -SurfacePointCounts @(96,384) `
        -SphereMeshModes $SphereMeshModes
    if ($LASTEXITCODE -ne 0) { throw 'Kirkwood multi-charge matrix failed' }
}
$rows += Invoke-Step 'grid_born_convergence_diagnostic' {
    & (Join-Path $PSScriptRoot 'pbe_grid_born_convergence.ps1') -TitanExe $TitanExe `
        -WorkRoot (Join-Path $WorkRoot 'grid_born') -RadiiA @(3.0) `
        -Charges @(1.0) -SaltM @(0.0,0.15) -GridScales @(0.60,1.0) `
        -FocusModes @('off') -GridAnchors @('auto')
    if ($LASTEXITCODE -ne 0) { throw 'grid Born convergence diagnostic failed' }
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
$rows += Invoke-Step 'treecode_interaction_list' {
    & (Join-Path $PSScriptRoot 'pbe_tabi_treecode_interaction_list_smoke.ps1') `
        -WorkRoot (Join-Path $WorkRoot 'treecode_interaction_list')
    if ($LASTEXITCODE -ne 0) { throw 'treecode interaction-list smoke failed' }
}
$rows += Invoke-Step 'jcp2013_treecode_vs_direct' {
    & (Join-Path $PSScriptRoot 'pbe_tabi_jcp2013_treecode_full_smoke.ps1') `
        -TitanExe $TitanExe -WorkRoot (Join-Path $WorkRoot 'jcp2013_treecode')
    if ($LASTEXITCODE -ne 0) { throw 'JCP2013 treecode vs direct failed' }
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
