param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_snapshot_stability_scan'),
    [string[]]$PqrPaths = @(),
    [ValidateSet('off', 'diagnose', 'anchor', 'hysteresis', 'ensemble')]
    [string]$StabilityMode = 'off',
    [ValidateSet('auto', 'per-frame', 'first-frame', 'fixed-center')]
    [string]$GridAnchor = 'auto',
    [switch]$JitterProbe
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootTrimmed = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootTrimmed + '\'
    if (-not ($full.Equals($rootTrimmed, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside PBE snapshot stability scan root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    $lines = @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       2.050   1.350   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1       1.450   2.350   0.000  1.00  0.00           O',
        'ATOM      5  CB  ALA A   1       1.950  -0.750   1.250  1.00  0.00           C',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Write-SyntheticPqr([string]$Path, [double]$X, [double]$Y = 0.0, [double]$Z = 0.0) {
    $line = ('ATOM      1  CA  GLY     1 {0,10:F7} {1,10:F7} {2,10:F7} {3,8:F4} {4,8:F4}       C' -f $X, $Y, $Z, 0.25, 1.70)
    Set-Content -LiteralPath $Path -Value @($line, 'END') -Encoding ASCII
}

function Write-JitteredPqr([string]$Source, [string]$Dest, [double]$Dx) {
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($line in Get-Content -LiteralPath $Source) {
        if ($line -match '^(ATOM|HETATM)\s+') {
            $fields = @($line -split '\s+' | Where-Object { $_ -ne '' })
            if ($fields.Count -ge 10) {
                $x = [double]::Parse($fields[5], $InvariantCulture) + $Dx
                $fields[5] = $x.ToString('F7', $InvariantCulture)
                $out.Add(($fields -join ' '))
                continue
            }
        }
        $out.Add($line)
    }
    Set-Content -LiteralPath $Dest -Value $out -Encoding ASCII
}

function Get-RowValue([object]$Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}

function Convert-ToDoubleOrZero([object]$Value) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return 0.0 }
    return [double]::Parse($text, $InvariantCulture)
}

function Get-Norm3([double]$X, [double]$Y, [double]$Z) {
    return [math]::Sqrt($X * $X + $Y * $Y + $Z * $Z)
}

function Classify-SnapshotPair([object]$Row, [object]$Previous,
                               [double]$ScanEnergyDelta,
                               [double]$ScanCenterDelta,
                               [double]$ScanOriginDelta) {
    if ($null -eq $Previous) { return 'first_frame' }
    $nearInteger = [int](Convert-ToDoubleOrZero (Get-RowValue $Row 'near_integer_atoms_1e-3' '0'))
    $nearHalf = [int](Convert-ToDoubleOrZero (Get-RowValue $Row 'near_half_atoms_1e-3' '0'))
    $soluteDelta = [int](Convert-ToDoubleOrZero (Get-RowValue $Row 'solute_delta' '0'))
    $iscmDelta = [int](Convert-ToDoubleOrZero (Get-RowValue $Row 'iscm_boundary_delta' '0'))
    $energyDelta = [math]::Abs($ScanEnergyDelta)
    if ($ScanCenterDelta -gt 1.0e-8 -or $ScanOriginDelta -gt 1.0e-8) { return 'grid_center_shift' }
    if ($soluteDelta -ne 0) { return 'surface_midpoint_shift' }
    if ($iscmDelta -ne 0) { return 'iscm_boundary_shift' }
    if (($nearInteger + $nearHalf) -gt 0 -and $energyDelta -gt 1.0e-3) {
        return 'grid_boundary_candidate'
    }
    return 'not_grid_boundary'
}

function Invoke-StabilityCase([string]$Name, [string]$PqrPath, [string]$Anchor = 'auto',
                              [string]$Mode = 'off') {
    $runDir = Join-Path $WorkRoot $Name
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe',
        "-pbe-test-pqr=$(Convert-ToTitanPath $PqrPath)",
        '-pbe-test-grid-size=33',
        '-pbe-test-center=0,0,0',
        '-pbe-scale=1.0',
        '-pbe-margin=2.0',
        '-pbe-surface=vdw',
        '-pbe-probe=1.4',
        '-pbe-max-iter=120',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-3',
        '-pbe-max-grid=200000',
        '-pbe-map=off',
        '-pbe-stability-audit=detail',
        "-pbe-stability-mode=$Mode",
        "-pbe-grid-anchor=$Anchor",
        '-nt=1'
    )
    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII

    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for $Name"
        }
    } finally {
        Pop-Location
    }
    $auditCsv = Join-Path $runDir 'pbe_stability_audit.csv'
    if (-not (Test-Path -LiteralPath $auditCsv)) {
        throw "Missing pbe_stability_audit.csv for $Name"
    }
    $row = @(Import-Csv -LiteralPath $auditCsv)[0]
    $row | Add-Member -NotePropertyName fixture -NotePropertyValue $Name -Force
    $row | Add-Member -NotePropertyName pqr_path -NotePropertyValue $PqrPath -Force
    return $row
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}

New-CleanDir $WorkRoot
$inputRows = New-Object System.Collections.Generic.List[object]
if ($PqrPaths.Count -eq 0) {
    $syntheticRoot = Join-Path $WorkRoot 'synthetic_inputs'
    New-Item -ItemType Directory -Path $syntheticRoot -Force | Out-Null
    $cases = @(
        @{ name = 'atom_on_integer_grid_plane'; x = 0.0; y = 0.0; z = 0.0 },
        @{ name = 'atom_plus_1e-7_A'; x = 1.0e-7; y = 0.0; z = 0.0 },
        @{ name = 'atom_plus_1e-5_A'; x = 1.0e-5; y = 0.0; z = 0.0 },
        @{ name = 'atom_mid_cell_control'; x = 0.5; y = 0.5; z = 0.5 }
    )
    foreach ($case in $cases) {
        $pqrPath = Join-Path $syntheticRoot "$($case.name).pqr"
        Write-SyntheticPqr -Path $pqrPath -X ([double]$case.x) -Y ([double]$case.y) -Z ([double]$case.z)
        $inputRows.Add([PSCustomObject]@{ name = $case.name; path = $pqrPath })
    }
} else {
    $index = 0
    foreach ($pqrPathRaw in $PqrPaths) {
        foreach ($item in ([string]$pqrPathRaw -split ',')) {
            $trimmed = $item.Trim()
            if ([string]::IsNullOrWhiteSpace($trimmed)) { continue }
            if (-not (Test-Path -LiteralPath $trimmed)) { throw "Missing PQR snapshot: $trimmed" }
            ++$index
            $name = ('snapshot_{0:000}' -f $index)
            $inputRows.Add([PSCustomObject]@{ name = $name; path = (Resolve-Path -LiteralPath $trimmed).Path })
        }
    }
}

if ($JitterProbe) {
    $jitterRoot = Join-Path $WorkRoot 'jitter_inputs'
    New-Item -ItemType Directory -Path $jitterRoot -Force | Out-Null
    $extra = New-Object System.Collections.Generic.List[object]
    foreach ($item in $inputRows) {
        foreach ($dx in @(-1.0e-7, 1.0e-7)) {
            $tag = if ($dx -lt 0.0) { 'minus_1e-7_A' } else { 'plus_1e-7_A' }
            $dest = Join-Path $jitterRoot "$($item.name)_$tag.pqr"
            Write-JitteredPqr -Source $item.path -Dest $dest -Dx $dx
            $extra.Add([PSCustomObject]@{ name = "$($item.name)_$tag"; path = $dest })
        }
    }
    foreach ($item in $extra) { $inputRows.Add($item) }
}

$rows = New-Object System.Collections.Generic.List[object]
$previous = $null
foreach ($item in $inputRows) {
    $row = Invoke-StabilityCase -Name $item.name -PqrPath $item.path -Anchor $GridAnchor -Mode $StabilityMode
    $scanEnergyDelta = 0.0
    $scanCenterDelta = 0.0
    $scanOriginDelta = 0.0
    if ($null -ne $previous) {
        $scanEnergyDelta =
            (Convert-ToDoubleOrZero (Get-RowValue $row 'energy_kJmol' '0')) -
            (Convert-ToDoubleOrZero (Get-RowValue $previous 'energy_kJmol' '0'))
        $scanCenterDelta = Get-Norm3 `
            ((Convert-ToDoubleOrZero (Get-RowValue $row 'grid_center_x_A' '0')) -
             (Convert-ToDoubleOrZero (Get-RowValue $previous 'grid_center_x_A' '0'))) `
            ((Convert-ToDoubleOrZero (Get-RowValue $row 'grid_center_y_A' '0')) -
             (Convert-ToDoubleOrZero (Get-RowValue $previous 'grid_center_y_A' '0'))) `
            ((Convert-ToDoubleOrZero (Get-RowValue $row 'grid_center_z_A' '0')) -
             (Convert-ToDoubleOrZero (Get-RowValue $previous 'grid_center_z_A' '0')))
        $scanOriginDelta = Get-Norm3 `
            ((Convert-ToDoubleOrZero (Get-RowValue $row 'grid_origin_x_A' '0')) -
             (Convert-ToDoubleOrZero (Get-RowValue $previous 'grid_origin_x_A' '0'))) `
            ((Convert-ToDoubleOrZero (Get-RowValue $row 'grid_origin_y_A' '0')) -
             (Convert-ToDoubleOrZero (Get-RowValue $previous 'grid_origin_y_A' '0'))) `
            ((Convert-ToDoubleOrZero (Get-RowValue $row 'grid_origin_z_A' '0')) -
             (Convert-ToDoubleOrZero (Get-RowValue $previous 'grid_origin_z_A' '0')))
    }
    $scanLabel = Classify-SnapshotPair -Row $row -Previous $previous `
        -ScanEnergyDelta $scanEnergyDelta -ScanCenterDelta $scanCenterDelta `
        -ScanOriginDelta $scanOriginDelta
    $row | Add-Member -NotePropertyName scan_energy_delta_kJmol -NotePropertyValue ($scanEnergyDelta.ToString('R', $InvariantCulture)) -Force
    $row | Add-Member -NotePropertyName scan_grid_center_delta_A -NotePropertyValue ($scanCenterDelta.ToString('R', $InvariantCulture)) -Force
    $row | Add-Member -NotePropertyName scan_grid_origin_delta_A -NotePropertyValue ($scanOriginDelta.ToString('R', $InvariantCulture)) -Force
    $row | Add-Member -NotePropertyName scan_label -NotePropertyValue $scanLabel -Force
    $rows.Add($row)
    $previous = $row
}

$summaryCsv = Join-Path $WorkRoot 'pbe_snapshot_stability_summary.csv'
$rows | Select-Object fixture,pqr_path,time_ps,energy_kJmol,energy_delta_kJmol,
    stability_mode,grid_anchor_effective,anchor_status,anchor_margin_min_A,
    hysteresis_cells_kept,hysteresis_cells_released,
    native_energy_kJmol,stable_energy_kJmol,energy_delta_vs_native_kJmol,
    ensemble_samples,ensemble_mean_kJmol,ensemble_stddev_kJmol,
    ensemble_min_kJmol,ensemble_max_kJmol,
    scan_energy_delta_kJmol,
    grid_center_delta_A,grid_origin_delta_A,
    scan_grid_center_delta_A,scan_grid_origin_delta_A,
    min_atom_to_integer_grid_units,min_atom_to_half_grid_units,
    min_atom_to_grid_plane_A,near_integer_atoms_1e-6,
    near_integer_atoms_1e-4,near_integer_atoms_1e-3,
    near_half_atoms_1e-6,near_half_atoms_1e-4,near_half_atoms_1e-3,
    charge_cell_flip_candidates,surface_flip_candidates,
    midpoint_flip_candidates,ion_access_flip_candidates,
    solute_delta,solvent_delta,iscm_boundary_delta,status,suspect_stage,
    scan_label |
    Export-Csv -LiteralPath $summaryCsv -NoTypeInformation -Encoding ASCII

if ($PqrPaths.Count -eq 0) {
    $exact = @($rows | Where-Object { $_.fixture -eq 'atom_on_integer_grid_plane' })[0]
    $near = @($rows | Where-Object { $_.fixture -eq 'atom_plus_1e-7_A' })[0]
    $loose = @($rows | Where-Object { $_.fixture -eq 'atom_plus_1e-5_A' })[0]
    $control = @($rows | Where-Object { $_.fixture -eq 'atom_mid_cell_control' })[0]
    if ([int]$exact.'near_integer_atoms_1e-6' -lt 1) {
        throw 'Synthetic exact-boundary atom was not flagged at 1e-6 grid units.'
    }
    if ([int]$near.'near_integer_atoms_1e-6' -lt 1) {
        throw 'Synthetic 1e-7 A near-boundary atom was not flagged at 1e-6 grid units.'
    }
    if ([int]$loose.'near_integer_atoms_1e-4' -lt 1) {
        throw 'Synthetic 1e-5 A near-boundary atom was not flagged at 1e-4 grid units.'
    }
    if ([int]$control.'near_integer_atoms_1e-3' -ne 0 -or
        [int]$control.'near_half_atoms_1e-3' -lt 1) {
        throw 'Synthetic mid-cell control did not isolate half-grid proximity.'
    }
}

Write-Host "PBE snapshot stability scan PASS: $summaryCsv"
