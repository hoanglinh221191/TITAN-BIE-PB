param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_timing_compare'),
    [string]$GoldenCsv = (Join-Path $PSScriptRoot 'pbe_inputs_local\pbe_delphi_golden.csv'),
    [string]$ManifestPath = (Join-Path $PSScriptRoot 'pbe_oracle_manifest.csv'),
    [string]$PythonExe = 'C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe',
    [string]$PyDelPhiRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) '.codex_tmp\pyDelPhi'),
    [string[]]$Fixtures = @(),
    [int]$RepeatCount = 2,
    [ValidateSet('auto','packed-list','row-span','avx2','salt-tiled')]
    [string]$NwtHomoKernel = 'auto',
    [double]$MinimumComparableSpeedRatio = 1.05
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
function Convert-ToDoubleOrZero([object]$Value) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return 0.0 }
    return [double]::Parse($text, $InvariantCulture)
}
function Convert-ToDoubleOrNull([object]$Value) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $null }
    return [double]::Parse($text, $InvariantCulture)
}
function Get-RowValue([object]$Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}
function Get-BottleneckBucket([object]$Row) {
    if ($null -eq $Row) { return 'NA' }
    $bestName = 'NA'
    $bestValue = -1.0
    foreach ($name in @(
        'setup_seconds','surface_seconds','epsilon_seconds',
        'boundary_seconds','solver_seconds','energy_seconds',
        'map_io_seconds','rpb_source_seconds','focus_setup_seconds',
        'probe_seconds','focus_parent_solve_seconds',
        'focus_child_setup_seconds','focus_boundary_interp_seconds',
        'focus_child_surface_seconds','focus_child_epsilon_seconds',
        'focus_child_solver_seconds','probe_energy_seconds',
        'audit_output_seconds','nwt_half_update_seconds',
        'nwt_rms_reduce_seconds','nwt_operator_residual_seconds',
        'nwt_source_update_seconds','nwt_progress_io_seconds',
        'nwt_coeff_pack_seconds','nwt_homo_update_seconds',
        'nwt_hetero_update_seconds','nwt_nonlin_term_seconds',
        'nwt_reduce_seconds','nwt_homo_stencil_seconds',
        'nwt_homo_salt_math_seconds','nwt_homo_write_seconds',
        'nwt_thread_wait_seconds','nwt_block_overhead_seconds',
        'linear_half_update_seconds',
        'epsilon_center_density_seconds','epsilon_midpoint_density_seconds',
        'epsilon_ion_exclusion_seconds','epsilon_finalize_seconds'
    )) {
        $value = Convert-ToDoubleOrNull (Get-RowValue $Row $name 'NA')
        if ($null -ne $value -and $value -gt $bestValue) {
            $bestValue = $value
            $bestName = $name -replace '_seconds$',''
        }
    }
    return $bestName
}
function Set-RowValue([object]$Row, [string]$Name, [string]$Value) {
    if ($Row.PSObject.Properties.Name -contains $Name) {
        $Row.$Name = $Value
    } else {
        $Row | Add-Member -NotePropertyName $Name -NotePropertyValue $Value
    }
}
function New-CleanDir([string]$Path) {
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}
function Get-HashText([string]$Path) {
    if ([string]::IsNullOrWhiteSpace($Path) -or -not (Test-Path -LiteralPath $Path)) { return 'NA' }
    try { return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant() }
    catch { return 'hash_error' }
}
function Get-PyDelPhiCommit([string]$Root) {
    if ([string]::IsNullOrWhiteSpace($Root) -or -not (Test-Path -LiteralPath $Root)) { return 'NA' }
    try {
        $commit = & git -C $Root rev-parse HEAD 2>$null
        if ($LASTEXITCODE -eq 0 -and -not [string]::IsNullOrWhiteSpace($commit)) {
            return ([string]$commit).Trim()
        }
    } catch {}
    return 'unknown'
}
function Add-OracleManifestRow([hashtable]$Values) {
    if ([string]::IsNullOrWhiteSpace($ManifestPath)) { return }
    $columns = @(
        'timestamp_utc','producer_script','command_class','fixture',
        'oracle_engine','python_exe','pyroot','pyroot_git_commit',
        'workroot','input_hash','output_hash','output_path',
        'thread_count','status','note'
    )
    $ordered = [ordered]@{}
    foreach ($column in $columns) {
        $ordered[$column] = if ($Values.ContainsKey($column)) { [string]$Values[$column] } else { 'NA' }
    }
    $row = [PSCustomObject]$ordered
    $dir = Split-Path -Parent $ManifestPath
    if (-not [string]::IsNullOrWhiteSpace($dir)) {
        New-Item -ItemType Directory -Path $dir -Force | Out-Null
    }
    if (Test-Path -LiteralPath $ManifestPath) {
        $row | Export-Csv -LiteralPath $ManifestPath -NoTypeInformation -Encoding ASCII -Append
    } else {
        $row | Export-Csv -LiteralPath $ManifestPath -NoTypeInformation -Encoding ASCII
    }
}
function Get-MedianText([object[]]$Values) {
    $numbers = New-Object System.Collections.Generic.List[double]
    foreach ($value in $Values) {
        $number = Convert-ToDoubleOrNull $value
        if ($null -ne $number) { $numbers.Add($number) }
    }
    if ($numbers.Count -eq 0) { return 'NA' }
    $sorted = @($numbers | Sort-Object)
    $middle = [int]($sorted.Count / 2)
    $median = if (($sorted.Count % 2) -eq 1) {
        $sorted[$middle]
    } else {
        0.5 * ($sorted[$middle - 1] + $sorted[$middle])
    }
    return ([math]::Round($median, 9)).ToString($InvariantCulture)
}
function Add-NullableSeconds([double]$Current, [object]$Value) {
    $number = Convert-ToDoubleOrNull $Value
    if ($null -eq $number) { return $Current }
    return $Current + $number
}
function New-MedianRows([object[]]$Rows) {
    $numericColumns = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)
    foreach ($name in @(
        'total_seconds','setup_seconds','surface_seconds','epsilon_seconds',
        'boundary_seconds','solver_seconds','energy_seconds','map_io_seconds',
        'rpb_source_seconds','focus_setup_seconds','probe_seconds',
        'focus_parent_solve_seconds','focus_child_setup_seconds',
        'focus_boundary_interp_seconds','focus_child_surface_seconds',
        'focus_child_epsilon_seconds','focus_child_solver_seconds',
        'probe_energy_seconds','audit_output_seconds',
        'nwt_half_update_seconds','nwt_rms_reduce_seconds',
        'nwt_operator_residual_seconds','nwt_source_update_seconds',
        'nwt_progress_io_seconds',
        'nwt_coeff_pack_seconds','nwt_homo_update_seconds',
        'nwt_hetero_update_seconds','nwt_nonlin_term_seconds',
        'nwt_reduce_seconds','nwt_homo_stencil_seconds',
        'nwt_homo_salt_math_seconds','nwt_homo_write_seconds',
        'nwt_thread_wait_seconds','nwt_block_overhead_seconds',
        'linear_half_update_seconds',
        'nwt_homo_cells','nwt_homo_linear_cells','nwt_homo_salt_cells',
        'nwt_homo_updates','nwt_dispatch_calls','nwt_homo_ns_per_update',
        'epsilon_center_density_seconds',
        'epsilon_midpoint_density_seconds','epsilon_ion_exclusion_seconds',
        'epsilon_finalize_seconds',
        'iterations','blocks','speed_ratio_titan_vs_pydelphi',
        'final_energy_kJmol','titan_threads','pydelphi_threads'
    )) { [void]$numericColumns.Add($name) }

    $medianRows = New-Object System.Collections.Generic.List[object]
    $groups = $Rows | Group-Object -Property engine,fixture,mode
    foreach ($group in $groups) {
        $first = $group.Group[0]
        $ordered = [ordered]@{}
        foreach ($property in $first.PSObject.Properties.Name) {
            if ($property -eq 'repeat') { continue }
            if ($property -eq 'speed_ratio_titan_vs_pydelphi') {
                $ordered[$property] = 'NA'
            } elseif ($numericColumns.Contains($property)) {
                $ordered[$property] = Get-MedianText @($group.Group | ForEach-Object { $_.$property })
            } else {
                $ordered[$property] = $first.$property
            }
        }
        $ordered['timing_stat'] = 'median'
        $ordered['repeat_count'] = $group.Count
        $medianRows.Add([PSCustomObject]$ordered)
    }

    $pydelphiTotals = @{}
    foreach ($row in $medianRows) {
        if ($row.engine -ne 'pydelphi') { continue }
        $value = Convert-ToDoubleOrNull $row.total_seconds
        if ($null -ne $value -and $value -gt 0.0) {
            $pydelphiTotals[$row.fixture] = $value
        }
    }
    foreach ($row in $medianRows) {
        if ($row.engine -ne 'titan' -and $row.engine -ne 'titan_child_scope') { continue }
        $titanSeconds = Convert-ToDoubleOrNull $row.total_seconds
        if ($null -ne $titanSeconds -and $titanSeconds -gt 0.0 -and
            $pydelphiTotals.ContainsKey($row.fixture)) {
            $row.speed_ratio_titan_vs_pydelphi =
                [math]::Round($pydelphiTotals[$row.fixture] / $titanSeconds, 6)
        }
    }
    return $medianRows
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $GoldenCsv)) { throw "Missing golden CSV: $GoldenCsv" }
$Fixtures = @($Fixtures | ForEach-Object {
    ([string]$_ -split ',') | ForEach-Object { $_.Trim() } | Where-Object { $_ }
})
$RepeatCount = [math]::Max(1, $RepeatCount)
if ($RepeatCount -gt 1) {
    New-CleanDir $WorkRoot
    $rawRows = New-Object System.Collections.Generic.List[object]
    for ($repeat = 1; $repeat -le $RepeatCount; ++$repeat) {
        $repeatRoot = Join-Path $WorkRoot ("repeat_{0}" -f $repeat)
        & $PSCommandPath -TitanExe $TitanExe -WorkRoot $repeatRoot `
            -GoldenCsv $GoldenCsv -ManifestPath $ManifestPath -PythonExe $PythonExe `
            -PyDelPhiRoot $PyDelPhiRoot -Fixtures $Fixtures -RepeatCount 1 `
            -NwtHomoKernel $NwtHomoKernel `
            -MinimumComparableSpeedRatio 0 | Out-Null
        $repeatCsv = Join-Path $repeatRoot 'pbe_timing_compare.csv'
        foreach ($row in @(Import-Csv -LiteralPath $repeatCsv)) {
            $row | Add-Member -NotePropertyName repeat -NotePropertyValue $repeat
            $rawRows.Add($row)
        }
    }
    $rawCsv = Join-Path $WorkRoot 'pbe_timing_compare_repeats.csv'
    $rawRows | Export-Csv -LiteralPath $rawCsv -NoTypeInformation -Encoding ASCII
    $medianRows = New-MedianRows $rawRows.ToArray()
    $outCsv = Join-Path $WorkRoot 'pbe_timing_compare.csv'
    $medianRows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
    Add-OracleManifestRow @{
        timestamp_utc = (Get-Date).ToUniversalTime().ToString('o')
        producer_script = 'test/pbe_timing_compare.ps1'
        command_class = 'timing_compare_median'
        fixture = if ($Fixtures.Count -gt 0) { $Fixtures -join ';' } else { 'all_default' }
        oracle_engine = 'pyDelPhi'
        python_exe = $PythonExe
        pyroot = $PyDelPhiRoot
        pyroot_git_commit = Get-PyDelPhiCommit $PyDelPhiRoot
        workroot = $WorkRoot
        input_hash = Get-HashText $GoldenCsv
        output_hash = (Get-HashText $outCsv) + ';repeats=' + (Get-HashText $rawCsv)
        output_path = $outCsv
        thread_count = '2'
        status = 'median_written'
        note = "repeat_count=$RepeatCount"
    }
    if ($MinimumComparableSpeedRatio -gt 0.0) {
        $speedFailures = New-Object System.Collections.Generic.List[string]
        foreach ($row in $medianRows) {
            if ($row.engine -ne 'titan') { continue }
            if ($row.mode -eq 'focus-probe') { continue }
            $ratio = Convert-ToDoubleOrNull $row.speed_ratio_titan_vs_pydelphi
            if ($null -ne $ratio -and $ratio -lt $MinimumComparableSpeedRatio) {
                $speedFailures.Add(("{0}: ratio={1} target={2}" -f $row.fixture, $ratio, $MinimumComparableSpeedRatio))
            }
        }
        if ($speedFailures.Count -gt 0) {
            throw ("Comparable TITAN/PyDelPhi speed gate failed: " + ($speedFailures -join '; '))
        }
    }
    Write-Host "PBE timing compare median COMPLETE: $outCsv"
    Write-Host "PBE timing compare raw repeats: $rawCsv"
    return
}
New-CleanDir $WorkRoot

$workloads = @(
    @{ fixture = 'barnase_linear_salt015'; mode = 'linear'; salt = '0.15' },
    @{ fixture = 'barnase_nwt'; mode = 'nwt'; salt = '0.15' },
    @{ fixture = '5tif_gaussian_cutoff_salt015'; mode = 'gaussian-cutoff'; salt = '0.15' },
    @{ fixture = '5tif_rpbe_salt015'; mode = 'rpbe'; salt = '0.15' },
    @{ fixture = 'barnase_focus_probe'; mode = 'focus-probe'; salt = '0.15' }
)
if ($Fixtures.Count -gt 0) {
    $wantedFixtures = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)
    foreach ($fixtureName in $Fixtures) {
        if (-not [string]::IsNullOrWhiteSpace($fixtureName)) {
            [void]$wantedFixtures.Add($fixtureName.Trim())
        }
    }
    $workloads = @($workloads | Where-Object { $wantedFixtures.Contains([string]$_.fixture) })
} else {
    $wantedFixtures = $null
}

$rowsOut = New-Object System.Collections.Generic.List[object]
$referenceSmoke = Join-Path $PSScriptRoot 'pbe_delphi_reference_smoke.ps1'
$goldenRows = @(Import-Csv -LiteralPath $GoldenCsv)
foreach ($workload in $workloads) {
    $gold = @($goldenRows | Where-Object { $_.fixture -eq $workload.fixture })
    if ($gold.Count -eq 0) { continue }
    foreach ($row in $gold) {
        Set-RowValue $row 'run_smoke' 'true'
        Set-RowValue $row 'pbe_timing' 'summary'
        Set-RowValue $row 'temp_K' '297.15'
        if ($workload.mode -eq 'linear') {
            Set-RowValue $row 'test_energy_convention' 'delphi-iscm'
        }
        if ($workload.mode -eq 'nwt') {
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
            Set-RowValue $row 'pbe_nwt_kernel' 'auto'
            Set-RowValue $row 'pbe_nwt_homo_kernel' $NwtHomoKernel
            Set-RowValue $row 'pbe_nwt_progress' 'off'
            Set-RowValue $row 'test_operator_convention' 'delphi'
            Set-RowValue $row 'test_boundary_convention' 'delphi-grid'
            Set-RowValue $row 'test_field_units' 'native'
            Set-RowValue $row 'temp_K' '297.15'
        }
    }
    $caseRoot = Join-Path $WorkRoot $workload.fixture
    New-Item -ItemType Directory -Path $caseRoot -Force | Out-Null
    $tmpGolden = Join-Path $caseRoot 'golden.csv'
    $gold | Export-Csv -LiteralPath $tmpGolden -NoTypeInformation -Encoding ASCII

    $measure = Measure-Command {
        & $referenceSmoke -TitanExe $TitanExe -WorkRoot (Join-Path $caseRoot 'titan') -GoldenCsv $tmpGolden -PyDelPhiRoot $PyDelPhiRoot | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "TITAN reference smoke failed for $($workload.fixture)" }
    }
    $runDir = Join-Path (Join-Path $caseRoot 'titan') $workload.fixture
    $energyCsv = Join-Path $runDir 'pbe_energy.csv'
    $energy = @(Import-Csv -LiteralPath $energyCsv)[0]
    $timingCsv = Join-Path $runDir 'pbe_timing.csv'
    $timing = $null
    if (Test-Path -LiteralPath $timingCsv) {
        $timing = @(Import-Csv -LiteralPath $timingCsv)[0]
    }
    $titanTotal = Get-RowValue $timing 'total_seconds' ([string][math]::Round($measure.TotalSeconds, 6))
    $resultCsv = Join-Path (Join-Path $caseRoot 'titan') 'pbe_delphi_reference_results.csv'
    $parityStatus = 'report_only'
    if (Test-Path -LiteralPath $resultCsv) {
        $resultRows = @(Import-Csv -LiteralPath $resultCsv | Where-Object { $_.fixture -eq $workload.fixture })
        if ($resultRows.Count -gt 0) {
            $parityStatus = if ($resultRows[0].pass -eq 'True') { 'pass' } else { 'mismatch_report_only' }
        }
    }
    $rowsOut.Add([PSCustomObject]@{
        engine = 'titan'
        fixture = $workload.fixture
        mode = $workload.mode
        parity_status = $parityStatus
        salt = $workload.salt
        grid_shape = "$($energy.grid_nx)x$($energy.grid_ny)x$($energy.grid_nz)"
        atom_count = $energy.atom_count
        titan_threads = Get-RowValue $timing 'num_threads' (Get-RowValue $gold[0] 'nt' '2')
        pydelphi_threads = 'NA'
        scope = Get-RowValue $timing 'scope' $workload.mode
        total_seconds = $titanTotal
        setup_seconds = Get-RowValue $timing 'setup_seconds' 'NA'
        surface_seconds = Get-RowValue $timing 'surface_seconds' 'NA'
        epsilon_seconds = Get-RowValue $timing 'epsilon_seconds' 'NA'
        boundary_seconds = Get-RowValue $timing 'boundary_seconds' 'NA'
        solver_seconds = Get-RowValue $timing 'solver_seconds' 'NA'
        energy_seconds = Get-RowValue $timing 'energy_seconds' 'NA'
        map_io_seconds = Get-RowValue $timing 'map_io_seconds' 'NA'
        rpb_source_seconds = Get-RowValue $timing 'rpb_source_seconds' 'NA'
        focus_setup_seconds = Get-RowValue $timing 'focus_setup_seconds' 'NA'
        probe_seconds = Get-RowValue $timing 'probe_seconds' 'NA'
        focus_parent_solve_seconds = Get-RowValue $timing 'focus_parent_solve_seconds' 'NA'
        focus_child_setup_seconds = Get-RowValue $timing 'focus_child_setup_seconds' 'NA'
        focus_boundary_interp_seconds = Get-RowValue $timing 'focus_boundary_interp_seconds' 'NA'
        focus_child_surface_seconds = Get-RowValue $timing 'focus_child_surface_seconds' 'NA'
        focus_child_epsilon_seconds = Get-RowValue $timing 'focus_child_epsilon_seconds' 'NA'
        focus_child_solver_seconds = Get-RowValue $timing 'focus_child_solver_seconds' 'NA'
        probe_energy_seconds = Get-RowValue $timing 'probe_energy_seconds' 'NA'
        audit_output_seconds = Get-RowValue $timing 'audit_output_seconds' 'NA'
        nwt_half_update_seconds = Get-RowValue $timing 'nwt_half_update_seconds' 'NA'
        nwt_rms_reduce_seconds = Get-RowValue $timing 'nwt_rms_reduce_seconds' 'NA'
        nwt_operator_residual_seconds = Get-RowValue $timing 'nwt_operator_residual_seconds' 'NA'
        nwt_source_update_seconds = Get-RowValue $timing 'nwt_source_update_seconds' 'NA'
        nwt_progress_io_seconds = Get-RowValue $timing 'nwt_progress_io_seconds' 'NA'
        nwt_coeff_pack_seconds = Get-RowValue $timing 'nwt_coeff_pack_seconds' 'NA'
        nwt_homo_update_seconds = Get-RowValue $timing 'nwt_homo_update_seconds' 'NA'
        nwt_hetero_update_seconds = Get-RowValue $timing 'nwt_hetero_update_seconds' 'NA'
        nwt_nonlin_term_seconds = Get-RowValue $timing 'nwt_nonlin_term_seconds' 'NA'
        nwt_reduce_seconds = Get-RowValue $timing 'nwt_reduce_seconds' 'NA'
        nwt_homo_stencil_seconds = Get-RowValue $timing 'nwt_homo_stencil_seconds' 'NA'
        nwt_homo_salt_math_seconds = Get-RowValue $timing 'nwt_homo_salt_math_seconds' 'NA'
        nwt_homo_write_seconds = Get-RowValue $timing 'nwt_homo_write_seconds' 'NA'
        nwt_thread_wait_seconds = Get-RowValue $timing 'nwt_thread_wait_seconds' 'NA'
        nwt_block_overhead_seconds = Get-RowValue $timing 'nwt_block_overhead_seconds' 'NA'
        linear_half_update_seconds = Get-RowValue $timing 'linear_half_update_seconds' 'NA'
        nwt_homo_cells = Get-RowValue $timing 'nwt_homo_cells' 'NA'
        nwt_homo_linear_cells = Get-RowValue $timing 'nwt_homo_linear_cells' 'NA'
        nwt_homo_salt_cells = Get-RowValue $timing 'nwt_homo_salt_cells' 'NA'
        nwt_homo_updates = Get-RowValue $timing 'nwt_homo_updates' 'NA'
        nwt_dispatch_calls = Get-RowValue $timing 'nwt_dispatch_calls' 'NA'
        nwt_homo_ns_per_update = Get-RowValue $timing 'nwt_homo_ns_per_update' 'NA'
        epsilon_center_density_seconds = Get-RowValue $timing 'epsilon_center_density_seconds' 'NA'
        epsilon_midpoint_density_seconds = Get-RowValue $timing 'epsilon_midpoint_density_seconds' 'NA'
        epsilon_ion_exclusion_seconds = Get-RowValue $timing 'epsilon_ion_exclusion_seconds' 'NA'
        epsilon_finalize_seconds = Get-RowValue $timing 'epsilon_finalize_seconds' 'NA'
        iterations = Get-RowValue $timing 'iterations' $energy.iterations_solvent
        blocks = Get-RowValue $timing 'blocks' $energy.nl_blocks
        speed_ratio_titan_vs_pydelphi = 'NA'
        bottleneck_bucket = Get-BottleneckBucket $timing
        final_energy_kJmol = $energy.E_rxn_kJmol
        status = $energy.status
    })
    if ($workload.mode -eq 'focus-probe' -and $null -ne $timing) {
        $childTotal = 0.0
        foreach ($name in @(
            'focus_child_setup_seconds','focus_boundary_interp_seconds',
            'focus_child_surface_seconds','focus_child_epsilon_seconds',
            'focus_child_solver_seconds','probe_energy_seconds'
        )) {
            $childTotal = Add-NullableSeconds $childTotal (Get-RowValue $timing $name 'NA')
        }
        $rowsOut.Add([PSCustomObject]@{
            engine = 'titan_child_scope'
            fixture = $workload.fixture
            mode = $workload.mode
            parity_status = 'report_only_child_scope'
            salt = $workload.salt
            grid_shape = "$($energy.grid_nx)x$($energy.grid_ny)x$($energy.grid_nz)"
            atom_count = Get-RowValue $timing 'focus_child_atom_count' $energy.atom_count
            titan_threads = Get-RowValue $timing 'num_threads' (Get-RowValue $gold[0] 'nt' '2')
            pydelphi_threads = 'NA'
            scope = 'focus_child_report_only_probe'
            total_seconds = [math]::Round($childTotal, 9).ToString($InvariantCulture)
            setup_seconds = 'NA'
            surface_seconds = Get-RowValue $timing 'focus_child_surface_seconds' 'NA'
            epsilon_seconds = Get-RowValue $timing 'focus_child_epsilon_seconds' 'NA'
            boundary_seconds = Get-RowValue $timing 'focus_boundary_interp_seconds' 'NA'
            solver_seconds = Get-RowValue $timing 'focus_child_solver_seconds' 'NA'
            energy_seconds = Get-RowValue $timing 'probe_energy_seconds' 'NA'
            map_io_seconds = 'NA'
            rpb_source_seconds = 'NA'
            focus_setup_seconds = Get-RowValue $timing 'focus_child_setup_seconds' 'NA'
            probe_seconds = Get-RowValue $timing 'probe_seconds' 'NA'
            focus_parent_solve_seconds = 'NA'
            focus_child_setup_seconds = Get-RowValue $timing 'focus_child_setup_seconds' 'NA'
            focus_boundary_interp_seconds = Get-RowValue $timing 'focus_boundary_interp_seconds' 'NA'
            focus_child_surface_seconds = Get-RowValue $timing 'focus_child_surface_seconds' 'NA'
            focus_child_epsilon_seconds = Get-RowValue $timing 'focus_child_epsilon_seconds' 'NA'
            focus_child_solver_seconds = Get-RowValue $timing 'focus_child_solver_seconds' 'NA'
            probe_energy_seconds = Get-RowValue $timing 'probe_energy_seconds' 'NA'
            audit_output_seconds = Get-RowValue $timing 'audit_output_seconds' 'NA'
            nwt_half_update_seconds = 'NA'
            nwt_rms_reduce_seconds = 'NA'
            nwt_operator_residual_seconds = 'NA'
            nwt_source_update_seconds = 'NA'
            nwt_progress_io_seconds = 'NA'
            nwt_coeff_pack_seconds = 'NA'
            nwt_homo_update_seconds = 'NA'
            nwt_hetero_update_seconds = 'NA'
            nwt_nonlin_term_seconds = 'NA'
            nwt_reduce_seconds = 'NA'
            nwt_homo_stencil_seconds = 'NA'
            nwt_homo_salt_math_seconds = 'NA'
            nwt_homo_write_seconds = 'NA'
            nwt_thread_wait_seconds = 'NA'
            nwt_block_overhead_seconds = 'NA'
            linear_half_update_seconds = 'NA'
            nwt_homo_cells = 'NA'
            nwt_homo_linear_cells = 'NA'
            nwt_homo_salt_cells = 'NA'
            nwt_homo_updates = 'NA'
            nwt_dispatch_calls = 'NA'
            nwt_homo_ns_per_update = 'NA'
            epsilon_center_density_seconds = 'NA'
            epsilon_midpoint_density_seconds = 'NA'
            epsilon_ion_exclusion_seconds = 'NA'
            epsilon_finalize_seconds = 'NA'
            iterations = Get-RowValue $timing 'iterations' $energy.iterations_solvent
            blocks = Get-RowValue $timing 'blocks' $energy.nl_blocks
            speed_ratio_titan_vs_pydelphi = 'NA'
            bottleneck_bucket = 'focus_child_solver'
            final_energy_kJmol = $energy.E_rxn_kJmol
            status = 'report_only_child_scope'
        })
    }
}

$exportScript = Join-Path $PSScriptRoot 'pbe_delphi_export_maps.ps1'
if (Test-Path -LiteralPath $PythonExe) {
    $oracleRoot = Join-Path $WorkRoot 'pydelphi'
    $oracleJobs = @(
        @{ fixture = 'barnase_linear_salt015'; mode = 'linear'; flags = @() },
        @{ fixture = 'barnase_nwt'; mode = 'nwt'; flags = @('-FullOperatorOracle', '-FullFieldOracle') },
        @{ fixture = '5tif_gaussian_cutoff_salt015'; mode = 'gaussian-cutoff'; flags = @('-FullGaussianOracle') },
        @{ fixture = '5tif_rpbe_salt015'; mode = 'rpbe'; flags = @('-FullRpbOracle') },
        @{ fixture = 'barnase_focus_probe'; mode = 'focus-probe'; flags = @('-FullFocusProbeOracle') }
    )
    if ($null -ne $wantedFixtures) {
        $oracleJobs = @($oracleJobs | Where-Object { $wantedFixtures.Contains([string]$_.fixture) })
    }
    foreach ($job in $oracleJobs) {
        $jobRoot = Join-Path $oracleRoot $job.fixture
        $exportMeasure = Measure-Command {
            $exportParams = @{
                PythonExe = $PythonExe
                PyDelPhiRoot = $PyDelPhiRoot
                WorkRoot = $jobRoot
                ManifestPath = $ManifestPath
                FullMidpointSalt = '0.15'
                Fixture = @($job.fixture)
            }
            foreach ($flag in $job.flags) {
                $name = ([string]$flag).TrimStart('-')
                $exportParams[$name] = $true
            }
            & $exportScript @exportParams | Out-Null
        }
        $stagePath = Join-Path $jobRoot 'barnase_stage_oracle.tsv'
        if (-not (Test-Path -LiteralPath $stagePath)) {
            $rowsOut.Add([PSCustomObject]@{
                engine = 'pydelphi'
                fixture = $job.fixture
                mode = $job.mode
                parity_status = 'report_only_no_stage'
                salt = 'NA'
                grid_shape = 'NA'
                atom_count = 'NA'
                titan_threads = 'NA'
                pydelphi_threads = 'NA'
                scope = $job.mode
                total_seconds = [math]::Round($exportMeasure.TotalSeconds, 6)
                setup_seconds = 'NA'
                surface_seconds = 'NA'
                epsilon_seconds = 'NA'
                boundary_seconds = 'NA'
                solver_seconds = 'NA'
                energy_seconds = 'NA'
                map_io_seconds = 'NA'
                rpb_source_seconds = 'NA'
                focus_setup_seconds = 'NA'
                probe_seconds = 'NA'
                focus_parent_solve_seconds = 'NA'
                focus_child_setup_seconds = 'NA'
                focus_boundary_interp_seconds = 'NA'
                focus_child_surface_seconds = 'NA'
                focus_child_epsilon_seconds = 'NA'
                focus_child_solver_seconds = 'NA'
                probe_energy_seconds = 'NA'
                audit_output_seconds = 'NA'
                nwt_half_update_seconds = 'NA'
                nwt_rms_reduce_seconds = 'NA'
                nwt_operator_residual_seconds = 'NA'
                nwt_source_update_seconds = 'NA'
                nwt_progress_io_seconds = 'NA'
                nwt_coeff_pack_seconds = 'NA'
                nwt_homo_update_seconds = 'NA'
                nwt_hetero_update_seconds = 'NA'
                nwt_nonlin_term_seconds = 'NA'
                nwt_reduce_seconds = 'NA'
                nwt_homo_stencil_seconds = 'NA'
                nwt_homo_salt_math_seconds = 'NA'
                nwt_homo_write_seconds = 'NA'
                nwt_thread_wait_seconds = 'NA'
                nwt_block_overhead_seconds = 'NA'
                linear_half_update_seconds = 'NA'
                nwt_homo_cells = 'NA'
                nwt_homo_linear_cells = 'NA'
                nwt_homo_salt_cells = 'NA'
                nwt_homo_updates = 'NA'
                nwt_dispatch_calls = 'NA'
                nwt_homo_ns_per_update = 'NA'
                epsilon_center_density_seconds = 'NA'
                epsilon_midpoint_density_seconds = 'NA'
                epsilon_ion_exclusion_seconds = 'NA'
                epsilon_finalize_seconds = 'NA'
                iterations = 'NA'
                blocks = 'NA'
                speed_ratio_titan_vs_pydelphi = 'NA'
                bottleneck_bucket = 'NA'
                final_energy_kJmol = 'NA'
                status = 'report_only_no_stage'
            })
            continue
        }
        $stageRows = @(Import-Csv -LiteralPath $stagePath -Delimiter "`t")
        foreach ($stage in $stageRows) {
            $fixture = [string]$stage.fixture
            $stageSalt = if ($stage.PSObject.Properties.Name -contains 'salt_M') { $stage.salt_M } else { 'NA' }
            $rowsOut.Add([PSCustomObject]@{
                engine = 'pydelphi'
                fixture = $fixture
                mode = $job.mode
                parity_status = 'oracle'
                salt = $stageSalt
                grid_shape = "$($stage.grid_nx)x$($stage.grid_ny)x$($stage.grid_nz)"
                atom_count = $stage.pydelphi_atom_count
                titan_threads = 'NA'
                pydelphi_threads = Get-RowValue $stage 'pydelphi_threads' '2'
                scope = if ($job.mode -eq 'focus-probe') { 'focus_child_oracle_probe' } else { $job.mode }
                total_seconds = $stage.pydelphi_total_seconds
                setup_seconds = $stage.pydelphi_setup_seconds
                surface_seconds = $stage.pydelphi_surface_seconds
                epsilon_seconds = $stage.pydelphi_epsilon_seconds
                boundary_seconds = 'NA'
                solver_seconds = $stage.pydelphi_solver_seconds
                energy_seconds = $stage.pydelphi_energy_seconds
                map_io_seconds = $stage.pydelphi_map_io_seconds
                rpb_source_seconds = 'NA'
                focus_setup_seconds = 'NA'
                probe_seconds = 'NA'
                focus_parent_solve_seconds = 'NA'
                focus_child_setup_seconds = 'NA'
                focus_boundary_interp_seconds = 'NA'
                focus_child_surface_seconds = 'NA'
                focus_child_epsilon_seconds = 'NA'
                focus_child_solver_seconds = 'NA'
                probe_energy_seconds = 'NA'
                audit_output_seconds = 'NA'
                nwt_half_update_seconds = 'NA'
                nwt_rms_reduce_seconds = 'NA'
                nwt_operator_residual_seconds = 'NA'
                nwt_source_update_seconds = 'NA'
                nwt_progress_io_seconds = 'NA'
                nwt_coeff_pack_seconds = 'NA'
                nwt_homo_update_seconds = 'NA'
                nwt_hetero_update_seconds = 'NA'
                nwt_nonlin_term_seconds = 'NA'
                nwt_reduce_seconds = 'NA'
                nwt_homo_stencil_seconds = 'NA'
                nwt_homo_salt_math_seconds = 'NA'
                nwt_homo_write_seconds = 'NA'
                nwt_thread_wait_seconds = 'NA'
                nwt_block_overhead_seconds = 'NA'
                linear_half_update_seconds = 'NA'
                nwt_homo_cells = 'NA'
                nwt_homo_linear_cells = 'NA'
                nwt_homo_salt_cells = 'NA'
                nwt_homo_updates = 'NA'
                nwt_dispatch_calls = 'NA'
                nwt_homo_ns_per_update = 'NA'
                epsilon_center_density_seconds = 'NA'
                epsilon_midpoint_density_seconds = 'NA'
                epsilon_ion_exclusion_seconds = 'NA'
                epsilon_finalize_seconds = 'NA'
                iterations = $stage.pydelphi_iterations
                blocks = 'NA'
                speed_ratio_titan_vs_pydelphi = 'NA'
                bottleneck_bucket = Get-BottleneckBucket ([PSCustomObject]@{
                    setup_seconds = $stage.pydelphi_setup_seconds
                    surface_seconds = $stage.pydelphi_surface_seconds
                    epsilon_seconds = $stage.pydelphi_epsilon_seconds
                    boundary_seconds = 'NA'
                    solver_seconds = $stage.pydelphi_solver_seconds
                    energy_seconds = $stage.pydelphi_energy_seconds
                    map_io_seconds = $stage.pydelphi_map_io_seconds
                    rpb_source_seconds = 'NA'
                    focus_setup_seconds = 'NA'
                    probe_seconds = 'NA'
                })
                final_energy_kJmol = $stage.final_E_rxn_kJmol
                status = 'ok'
            })
        }
    }
} else {
    $rowsOut.Add([PSCustomObject]@{
        engine = 'pydelphi'
        fixture = 'missing_python'
        mode = 'oracle'
        parity_status = 'missing'
        salt = 'NA'
        grid_shape = 'NA'
        atom_count = 'NA'
        titan_threads = 'NA'
        pydelphi_threads = 'NA'
        scope = 'oracle'
        total_seconds = 'NA'
        setup_seconds = 'NA'
        surface_seconds = 'NA'
        epsilon_seconds = 'NA'
        boundary_seconds = 'NA'
        solver_seconds = 'NA'
        energy_seconds = 'NA'
        map_io_seconds = 'NA'
        rpb_source_seconds = 'NA'
        focus_setup_seconds = 'NA'
        probe_seconds = 'NA'
        focus_parent_solve_seconds = 'NA'
        focus_child_setup_seconds = 'NA'
        focus_boundary_interp_seconds = 'NA'
        focus_child_surface_seconds = 'NA'
        focus_child_epsilon_seconds = 'NA'
        focus_child_solver_seconds = 'NA'
        probe_energy_seconds = 'NA'
        audit_output_seconds = 'NA'
        nwt_half_update_seconds = 'NA'
        nwt_rms_reduce_seconds = 'NA'
        nwt_operator_residual_seconds = 'NA'
        nwt_source_update_seconds = 'NA'
        nwt_progress_io_seconds = 'NA'
        nwt_coeff_pack_seconds = 'NA'
        nwt_homo_update_seconds = 'NA'
        nwt_hetero_update_seconds = 'NA'
        nwt_nonlin_term_seconds = 'NA'
        nwt_reduce_seconds = 'NA'
        nwt_homo_stencil_seconds = 'NA'
        nwt_homo_salt_math_seconds = 'NA'
        nwt_homo_write_seconds = 'NA'
        nwt_thread_wait_seconds = 'NA'
        nwt_block_overhead_seconds = 'NA'
        linear_half_update_seconds = 'NA'
        nwt_homo_cells = 'NA'
        nwt_homo_linear_cells = 'NA'
        nwt_homo_salt_cells = 'NA'
        nwt_homo_updates = 'NA'
        nwt_dispatch_calls = 'NA'
        nwt_homo_ns_per_update = 'NA'
        epsilon_center_density_seconds = 'NA'
        epsilon_midpoint_density_seconds = 'NA'
        epsilon_ion_exclusion_seconds = 'NA'
        epsilon_finalize_seconds = 'NA'
        iterations = 'NA'
        blocks = 'NA'
        speed_ratio_titan_vs_pydelphi = 'NA'
        bottleneck_bucket = 'NA'
        final_energy_kJmol = 'NA'
        status = "missing $PythonExe"
    })
}

$pydelphiTotals = @{}
foreach ($row in $rowsOut) {
    if ($row.engine -ne 'pydelphi') { continue }
    $value = Convert-ToDoubleOrNull $row.total_seconds
    if ($null -ne $value -and $value -gt 0.0) {
        $pydelphiTotals[$row.fixture] = $value
    }
}
foreach ($row in $rowsOut) {
    if ($row.engine -ne 'titan' -and $row.engine -ne 'titan_child_scope') { continue }
    $titanSeconds = Convert-ToDoubleOrNull $row.total_seconds
    if ($null -ne $titanSeconds -and $titanSeconds -gt 0.0 -and
        $pydelphiTotals.ContainsKey($row.fixture)) {
        $row.speed_ratio_titan_vs_pydelphi =
            [math]::Round($pydelphiTotals[$row.fixture] / $titanSeconds, 6)
    }
}

$outCsv = Join-Path $WorkRoot 'pbe_timing_compare.csv'
$rowsOut | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Add-OracleManifestRow @{
    timestamp_utc = (Get-Date).ToUniversalTime().ToString('o')
    producer_script = 'test/pbe_timing_compare.ps1'
    command_class = 'timing_compare_single'
    fixture = if ($Fixtures.Count -gt 0) { $Fixtures -join ';' } else { 'all_default' }
    oracle_engine = 'pyDelPhi'
    python_exe = $PythonExe
    pyroot = $PyDelPhiRoot
    pyroot_git_commit = Get-PyDelPhiCommit $PyDelPhiRoot
    workroot = $WorkRoot
    input_hash = Get-HashText $GoldenCsv
    output_hash = Get-HashText $outCsv
    output_path = $outCsv
    thread_count = '2'
    status = 'written'
    note = 'single_run'
}
if ($MinimumComparableSpeedRatio -gt 0.0) {
    $speedFailures = New-Object System.Collections.Generic.List[string]
    foreach ($row in $rowsOut) {
        if ($row.engine -ne 'titan') { continue }
        if ($row.mode -eq 'focus-probe') { continue }
        $ratio = Convert-ToDoubleOrNull $row.speed_ratio_titan_vs_pydelphi
        if ($null -ne $ratio -and $ratio -lt $MinimumComparableSpeedRatio) {
            $speedFailures.Add(("{0}: ratio={1} target={2}" -f $row.fixture, $ratio, $MinimumComparableSpeedRatio))
        }
    }
    if ($speedFailures.Count -gt 0) {
        throw ("Comparable TITAN/PyDelPhi speed gate failed: " + ($speedFailures -join '; '))
    }
}

Write-Host "PBE timing compare COMPLETE: $outCsv"
