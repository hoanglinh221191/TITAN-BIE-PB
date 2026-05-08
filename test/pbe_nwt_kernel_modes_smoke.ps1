param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_nwt_kernel_modes_smoke'),
    [string]$GoldenCsv = (Join-Path $PSScriptRoot 'pbe_inputs_local\pbe_delphi_golden.csv'),
    [string[]]$Modes = @('packed-list', 'row-span', 'avx2', 'salt-tiled')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$ReferenceSmoke = Join-Path $PSScriptRoot 'pbe_delphi_reference_smoke.ps1'

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

function Set-RowValue([object]$Row, [string]$Name, [string]$Value) {
    if ($Row.PSObject.Properties.Name -contains $Name) {
        $Row.$Name = $Value
    } else {
        $Row | Add-Member -NotePropertyName $Name -NotePropertyValue $Value
    }
}

function Get-RowValue([object]$Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}

function New-CleanDir([string]$Path) {
    $full = [System.IO.Path]::GetFullPath($Path)
    if (Test-Path -LiteralPath $full) { Remove-Item -LiteralPath $full -Recurse -Force }
    New-Item -ItemType Directory -Path $full -Force | Out-Null
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $GoldenCsv)) { throw "Missing golden CSV: $GoldenCsv" }
if (-not (Test-Path -LiteralPath $ReferenceSmoke)) { throw "Missing reference smoke: $ReferenceSmoke" }

$Modes = @($Modes | ForEach-Object {
    ([string]$_ -split ',') | ForEach-Object { $_.Trim() } | Where-Object { $_ }
})
if ($Modes.Count -eq 0) { throw 'No NWT homo kernel modes requested.' }
if ($Modes -notcontains 'packed-list') { $Modes = @('packed-list') + $Modes }

New-CleanDir $WorkRoot
$goldRows = @(Import-Csv -LiteralPath $GoldenCsv | Where-Object { $_.fixture -eq 'barnase_nwt' })
if ($goldRows.Count -eq 0) { throw "Missing barnase_nwt row in $GoldenCsv" }

$summaryRows = New-Object 'System.Collections.Generic.List[object]'
foreach ($mode in $Modes) {
    if ($mode -notin @('packed-list', 'row-span', 'avx2', 'salt-tiled')) {
        throw "Unsupported NWT homo kernel mode: $mode"
    }
    $modeRoot = Join-Path $WorkRoot $mode
    New-Item -ItemType Directory -Path $modeRoot -Force | Out-Null
    $rows = @($goldRows | ForEach-Object { $_.PSObject.Copy() })
    foreach ($row in $rows) {
        Set-RowValue $row 'run_smoke' 'true'
        Set-RowValue $row 'nl_solver' 'newton'
        Set-RowValue $row 'nl_inner' 'sor'
        Set-RowValue $row 'nl_max_iter' '5000'
        Set-RowValue $row 'nl_inner_iter' '50'
        Set-RowValue $row 'nl_block' '50'
        Set-RowValue $row 'nl_dphi' '1e-6'
        Set-RowValue $row 'nl_residual' '1e-10'
        Set-RowValue $row 'nl_continuation' 'off'
        Set-RowValue $row 'nl_line_search' 'off'
        Set-RowValue $row 'parity_dump' 'summary'
        Set-RowValue $row 'pbe_nwt_kernel' 'delphi-half'
        Set-RowValue $row 'pbe_nwt_homo_kernel' $mode
        Set-RowValue $row 'pbe_nwt_progress' 'off'
        Set-RowValue $row 'pbe_timing' 'summary'
        Set-RowValue $row 'test_operator_convention' 'delphi'
        Set-RowValue $row 'test_boundary_convention' 'delphi-grid'
        Set-RowValue $row 'test_field_units' 'native'
        Set-RowValue $row 'temp_K' '297.15'
    }
    $tempGolden = Join-Path $modeRoot 'golden.csv'
    $rows | Export-Csv -LiteralPath $tempGolden -NoTypeInformation -Encoding ASCII
    & $ReferenceSmoke -TitanExe $TitanExe -WorkRoot $modeRoot -GoldenCsv $tempGolden | Out-Host
    if ($LASTEXITCODE -ne 0) { throw "Reference smoke failed for NWT homo kernel mode $mode" }

    $runDir = Join-Path $modeRoot 'barnase_nwt'
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    $timing = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_timing.csv'))[0]
    $result = @(Import-Csv -LiteralPath (Join-Path $modeRoot 'pbe_delphi_reference_results.csv') |
        Where-Object { $_.fixture -eq 'barnase_nwt' })[0]
    $summaryRows.Add([PSCustomObject]@{
        mode = $mode
        status = $energy.status
        nl_status = Get-RowValue $result 'nl_status' 'NA'
        pass = Get-RowValue $result 'pass' 'NA'
        iterations = Get-RowValue $timing 'iterations' $energy.iterations_solvent
        blocks = Get-RowValue $timing 'blocks' $energy.nl_blocks
        total_seconds = Get-RowValue $timing 'total_seconds' 'NA'
        solver_seconds = Get-RowValue $timing 'solver_seconds' 'NA'
        nwt_half_update_seconds = Get-RowValue $timing 'nwt_half_update_seconds' 'NA'
        nwt_homo_update_seconds = Get-RowValue $timing 'nwt_homo_update_seconds' 'NA'
        nwt_hetero_update_seconds = Get-RowValue $timing 'nwt_hetero_update_seconds' 'NA'
        nwt_homo_stencil_seconds = Get-RowValue $timing 'nwt_homo_stencil_seconds' 'NA'
        nwt_homo_salt_math_seconds = Get-RowValue $timing 'nwt_homo_salt_math_seconds' 'NA'
        nwt_homo_write_seconds = Get-RowValue $timing 'nwt_homo_write_seconds' 'NA'
        nwt_thread_wait_seconds = Get-RowValue $timing 'nwt_thread_wait_seconds' 'NA'
        nwt_block_overhead_seconds = Get-RowValue $timing 'nwt_block_overhead_seconds' 'NA'
        nwt_homo_cells = Get-RowValue $timing 'nwt_homo_cells' 'NA'
        nwt_homo_linear_cells = Get-RowValue $timing 'nwt_homo_linear_cells' 'NA'
        nwt_homo_salt_cells = Get-RowValue $timing 'nwt_homo_salt_cells' 'NA'
        nwt_homo_updates = Get-RowValue $timing 'nwt_homo_updates' 'NA'
        nwt_dispatch_calls = Get-RowValue $timing 'nwt_dispatch_calls' 'NA'
        nwt_homo_ns_per_update = Get-RowValue $timing 'nwt_homo_ns_per_update' 'NA'
        final_energy_kJmol = $energy.E_rxn_kJmol
        oracle_tolerance_kJmol = Get-RowValue $result 'tolerance_kJmol' 'NA'
    })
}

$baseline = @($summaryRows | Where-Object { $_.mode -eq 'packed-list' })[0]
$baselineEnergy = Convert-ToDouble $baseline.final_energy_kJmol
$baselineBlocks = [int](Convert-ToDouble $baseline.blocks)
$baselineIterations = [int](Convert-ToDouble $baseline.iterations)
$failures = New-Object 'System.Collections.Generic.List[string]'
foreach ($row in $summaryRows) {
    $energy = Convert-ToDouble $row.final_energy_kJmol
    $blocks = [int](Convert-ToDouble $row.blocks)
    $iterations = [int](Convert-ToDouble $row.iterations)
    $tolerance = Convert-ToDouble $row.oracle_tolerance_kJmol
    $delta = [math]::Abs($energy - $baselineEnergy)
    $row | Add-Member -NotePropertyName delta_vs_packed_kJmol -NotePropertyValue ($delta.ToString('R', $InvariantCulture))
    $row | Add-Member -NotePropertyName compare_status -NotePropertyValue 'ok'
    if ($blocks -ne $baselineBlocks -or $iterations -ne $baselineIterations) {
        $row.compare_status = 'block_or_iteration_mismatch'
        $failures.Add("$($row.mode): block/iteration mismatch")
    } elseif ($delta -gt $tolerance) {
        $row.compare_status = 'energy_mismatch'
        $failures.Add("$($row.mode): energy delta $delta > tolerance $tolerance")
    } elseif ([string]$row.pass -ne 'True') {
        $row.compare_status = 'oracle_gate_failed'
        $failures.Add("$($row.mode): oracle gate failed")
    }
}

$summaryCsv = Join-Path $WorkRoot 'pbe_nwt_kernel_modes_summary.csv'
$summaryRows | Export-Csv -LiteralPath $summaryCsv -NoTypeInformation -Encoding ASCII
if ($failures.Count -gt 0) {
    throw "PBE NWT kernel modes smoke FAIL: $($failures -join '; ')"
}
Write-Host "PBE NWT kernel modes smoke PASS: $summaryCsv"
