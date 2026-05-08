param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_vs_grid_compare'),
    [string[]]$Fixtures = @('sphere','two-atoms','arg','mini-protein','charged-peptide','compact-mini-protein'),
    [int]$RepeatCount = 1,
    [ValidateSet('accuracy','perf')]
    [string]$Mode = 'accuracy'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

$FixtureMap = @{
    'sphere' = 'sphere.pqr'
    'two-atoms' = 'two-atoms.pqr'
    'arg' = 'arg.pqr'
    'mini-protein' = 'mini-protein.pqr'
    'charged-peptide' = 'charged-peptide.pqr'
    'compact-mini-protein' = 'compact-mini-protein.pqr'
    '1a63' = '1a63.pqr'
}

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDoubleOrNull([object]$Value) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') {
        return $null
    }
    return [double]::Parse($text, $InvariantCulture)
}

function Get-RowValue([object]$Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) {
            return $value
        }
    }
    return $Default
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootFull = [System.IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    if (-not $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside TABI-vs-grid root: $full"
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
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function New-ArgFile([string]$RunDir, [string]$Backend, [string]$PqrPath, [string]$RunMode) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        "-pbe=$Backend",
        "-pbe-test-pqr=$(Convert-ToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-probe=1.4',
        '-pbe-scale=0.35',
        '-pbe-margin=2.0',
        '-pbe-max-grid=300000',
        '-pbe-map=off',
        '-pbe-max-iter=800',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-timing=summary',
        '-nt=1'
    )
    if ($Backend -eq 'tabi') {
        if ($RunMode -eq 'accuracy') {
            $argLines += '-pbe-tabi-sum=direct'
            $argLines += '-pbe-tabi-formulation=jcp2013'
            $argLines += '-pbe-tabi-energy-projection=boundary'
        } else {
            $argLines += '-pbe-tabi-sum=auto'
            $argLines += '-pbe-tabi-formulation=titan-v1'
            $argLines += '-pbe-tabi-energy-projection=boundary'
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Invoke-PbeCase([int]$Repeat, [string]$Fixture, [string]$Backend) {
    $pqrName = $FixtureMap[$Fixture]
    if ([string]::IsNullOrWhiteSpace($pqrName)) {
        throw "Unknown fixture: $Fixture"
    }
    $pqrPath = Join-Path $InputRoot $pqrName
    if (-not (Test-Path -LiteralPath $pqrPath)) {
        throw "Missing PQR fixture: $pqrPath"
    }

    $runDir = Join-Path $WorkRoot ("repeat_{0}\{1}_{2}" -f $Repeat, $Fixture, $Backend)
    New-CleanDir $runDir
    $argPath = New-ArgFile -RunDir $runDir -Backend $Backend -PqrPath $pqrPath -RunMode $Mode
    $measure = Measure-Command {
        Push-Location $runDir
        try {
            & $TitanExe -a $argPath | Out-Null
            if ($LASTEXITCODE -ne 0) {
                throw "TITAN exited with code $LASTEXITCODE for $Fixture/$Backend"
            }
        } finally {
            Pop-Location
        }
    }

    $energyPath = Join-Path $runDir 'pbe_energy.csv'
    $timingPath = Join-Path $runDir 'pbe_timing.csv'
    foreach ($required in @($energyPath, $timingPath)) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Missing PBE output for ${Fixture}/${Backend}: $required"
        }
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $timing = @(Import-Csv -LiteralPath $timingPath)[0]
    $audit = $null
    $auditPath = Join-Path $runDir 'pbe_tabi_formulation_audit.csv'
    if ($Backend -eq 'tabi') {
        if (-not (Test-Path -LiteralPath $auditPath)) {
            throw "Missing TABI formulation audit for ${Fixture}/${Backend}: $auditPath"
        }
        $audit = @(Import-Csv -LiteralPath $auditPath)[0]
    }
    if ($energy.status -ne 'ok') {
        throw "PBE case did not finish ok for ${Fixture}/${Backend}: $($energy.status)"
    }

    $energyKJ = Convert-ToDoubleOrNull $energy.E_rxn_kJmol
    if ($null -eq $energyKJ -or [double]::IsNaN($energyKJ) -or [double]::IsInfinity($energyKJ)) {
        throw "Non-finite PBE energy for ${Fixture}/${Backend}: $($energy.E_rxn_kJmol)"
    }

    return [PSCustomObject]@{
        repeat = $Repeat
        fixture = $Fixture
        backend = $Backend
        mode = Get-RowValue $timing 'mode' (Get-RowValue $energy 'solver_mode' 'NA')
        atom_count = $energy.atom_count
        grid_points = $energy.grid_points
        panels = if ($Backend -eq 'tabi') { Get-RowValue $timing 'tabi_surface_panels' '0' } else { '0' }
        energy_kJmol = $energy.E_rxn_kJmol
        delta_vs_grid_abs = 'NA'
        delta_vs_grid_rel = 'NA'
        energy_ratio_vs_grid = 'NA'
        kirkwood_error_rel = 'NA'
        oracle_engine = 'none'
        oracle_status = 'not_requested'
        formulation = if ($Backend -eq 'tabi') { Get-RowValue $audit 'formulation' 'NA' } else { 'grid' }
        energy_projection = if ($Backend -eq 'tabi') { Get-RowValue $audit 'energy_projection' 'NA' } else { 'grid' }
        fallback_used = if ($Backend -eq 'tabi') { Get-RowValue $audit 'fallback_used' 'NA' } else { 'false' }
        total_seconds = Get-RowValue $timing 'total_seconds' ([math]::Round($measure.TotalSeconds, 6).ToString($InvariantCulture))
        setup_seconds = Get-RowValue $timing 'setup_seconds' 'NA'
        surface_seconds = Get-RowValue $timing 'surface_seconds' 'NA'
        solver_seconds = Get-RowValue $timing 'solver_seconds' 'NA'
        energy_seconds = Get-RowValue $timing 'energy_seconds' 'NA'
        tabi_surface_build_seconds = Get-RowValue $timing 'tabi_surface_build_seconds' 'NA'
        tabi_surface_neighbor_seconds = Get-RowValue $timing 'tabi_surface_neighbor_seconds' 'NA'
        tabi_surface_mesh_seconds = Get-RowValue $timing 'tabi_surface_mesh_seconds' 'NA'
        tabi_solid_angle_seconds = Get-RowValue $timing 'tabi_solid_angle_seconds' 'NA'
        tabi_rhs_seconds = Get-RowValue $timing 'tabi_rhs_seconds' 'NA'
        tabi_gmres_seconds = Get-RowValue $timing 'tabi_gmres_seconds' 'NA'
        tabi_matvec_seconds_total = Get-RowValue $timing 'tabi_matvec_seconds_total' 'NA'
        tabi_matvec_calls = Get-RowValue $timing 'tabi_matvec_calls' 'NA'
        tabi_tree_build_seconds = Get-RowValue $timing 'tabi_tree_build_seconds' 'NA'
        tabi_tree_moment_seconds = Get-RowValue $timing 'tabi_tree_moment_seconds' 'NA'
        tabi_tree_traversal_seconds = Get-RowValue $timing 'tabi_tree_traversal_seconds' 'NA'
        tabi_tree_nodes = Get-RowValue $timing 'tabi_tree_nodes' 'NA'
        tabi_tree_leaf_count = Get-RowValue $timing 'tabi_tree_leaf_count' 'NA'
        tabi_near_interactions = Get-RowValue $timing 'tabi_near_interactions' 'NA'
        tabi_far_interactions = Get-RowValue $timing 'tabi_far_interactions' 'NA'
        tabi_status = Get-RowValue $timing 'tabi_status' 'NA'
        tabi_convergence_reason = Get-RowValue $timing 'tabi_convergence_reason' 'NA'
        reaction_energy_source = $energy.reaction_energy_source
        status = $energy.status
    }
}

function Get-MedianText([object[]]$Values) {
    $numbers = New-Object System.Collections.Generic.List[double]
    foreach ($value in $Values) {
        $parsed = Convert-ToDoubleOrNull $value
        if ($null -ne $parsed) {
            $numbers.Add($parsed)
        }
    }
    if ($numbers.Count -eq 0) {
        return 'NA'
    }
    $sorted = @($numbers | Sort-Object)
    $mid = [int]($sorted.Count / 2)
    $median = if (($sorted.Count % 2) -eq 1) {
        $sorted[$mid]
    } else {
        0.5 * ($sorted[$mid - 1] + $sorted[$mid])
    }
    return ([math]::Round($median, 9)).ToString($InvariantCulture)
}

function New-MedianRows([object[]]$Rows) {
    $numericColumns = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)
    foreach ($name in @(
        'atom_count','grid_points','panels','energy_kJmol',
        'delta_vs_grid_abs','delta_vs_grid_rel','energy_ratio_vs_grid',
        'kirkwood_error_rel','total_seconds',
        'setup_seconds','surface_seconds','solver_seconds','energy_seconds',
        'tabi_surface_build_seconds','tabi_surface_neighbor_seconds',
        'tabi_surface_mesh_seconds','tabi_solid_angle_seconds',
        'tabi_rhs_seconds','tabi_gmres_seconds','tabi_matvec_seconds_total',
        'tabi_matvec_calls','tabi_tree_build_seconds',
        'tabi_tree_moment_seconds','tabi_tree_traversal_seconds',
        'tabi_tree_nodes','tabi_tree_leaf_count','tabi_near_interactions',
        'tabi_far_interactions')) {
        [void]$numericColumns.Add($name)
    }
    $out = New-Object System.Collections.Generic.List[object]
    foreach ($group in ($Rows | Group-Object -Property fixture,backend)) {
        $first = $group.Group[0]
        $ordered = [ordered]@{}
        foreach ($property in $first.PSObject.Properties.Name) {
            if ($property -eq 'repeat') {
                continue
            }
            if ($numericColumns.Contains($property)) {
                $ordered[$property] = Get-MedianText @($group.Group | ForEach-Object { $_.$property })
            } else {
                $ordered[$property] = $first.$property
            }
        }
        $ordered['timing_stat'] = 'median'
        $ordered['repeat_count'] = $group.Count
        $out.Add([PSCustomObject]$ordered)
    }
    return $out
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $InputRoot)) {
    throw "Missing PBE fixture directory: $InputRoot"
}

$Fixtures = @($Fixtures | ForEach-Object {
    ([string]$_ -split ',') | ForEach-Object { $_.Trim() } | Where-Object { $_ }
})
if ($Fixtures.Count -eq 0) {
    throw 'No fixtures selected.'
}
$RepeatCount = [math]::Max(1, $RepeatCount)
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$rawRows = New-Object System.Collections.Generic.List[object]
for ($repeat = 1; $repeat -le $RepeatCount; ++$repeat) {
    foreach ($fixture in $Fixtures) {
        $grid = Invoke-PbeCase -Repeat $repeat -Fixture $fixture -Backend 'grid'
        $tabi = Invoke-PbeCase -Repeat $repeat -Fixture $fixture -Backend 'tabi'
        $gridEnergy = Convert-ToDoubleOrNull $grid.energy_kJmol
        $tabiEnergy = Convert-ToDoubleOrNull $tabi.energy_kJmol
        $grid.delta_vs_grid_abs = '0'
        $grid.delta_vs_grid_rel = '0'
        $grid.energy_ratio_vs_grid = '1'
        $delta = [math]::Abs($tabiEnergy - $gridEnergy)
        $tabi.delta_vs_grid_abs = $delta.ToString($InvariantCulture)
        $tabi.delta_vs_grid_rel = ($delta / [math]::Max(1.0, [math]::Abs($gridEnergy))).ToString($InvariantCulture)
        $ratioDenom = if ([math]::Abs($gridEnergy) -gt 1.0e-12) { $gridEnergy } else { 1.0e-12 }
        $tabi.energy_ratio_vs_grid = ($tabiEnergy / $ratioDenom).ToString($InvariantCulture)
        if ($fixture -eq 'sphere') {
            $tabiAudit = @(Import-Csv -LiteralPath (Join-Path $WorkRoot ("repeat_{0}\{1}_{2}\pbe_tabi_formulation_audit.csv" -f $repeat, $fixture, 'tabi')))[0]
            $refEnergy = Convert-ToDoubleOrNull $tabiAudit.reference_energy_kJmol
            if ($null -ne $refEnergy) {
                $tabi.kirkwood_error_rel = ([math]::Abs($tabiEnergy - $refEnergy) / [math]::Max(1.0e-12, [math]::Abs($refEnergy))).ToString($InvariantCulture)
                $grid.kirkwood_error_rel = ([math]::Abs($gridEnergy - $refEnergy) / [math]::Max(1.0e-12, [math]::Abs($refEnergy))).ToString($InvariantCulture)
            }
        }
        $rawRows.Add($grid)
        $rawRows.Add($tabi)
        Write-Host ("PBE TABI-vs-grid {0} mode={1} repeat={2}: grid={3} kJ/mol tabi={4} kJ/mol abs_delta={5}" -f
            $fixture, $Mode, $repeat, $grid.energy_kJmol, $tabi.energy_kJmol, $tabi.delta_vs_grid_abs)
    }
}

$rawCsv = Join-Path $WorkRoot 'pbe_tabi_vs_grid_compare_repeats.csv'
$rawRows | Export-Csv -LiteralPath $rawCsv -NoTypeInformation -Encoding ASCII
$medianRows = New-MedianRows $rawRows.ToArray()
$outCsv = Join-Path $WorkRoot 'pbe_tabi_vs_grid_compare.csv'
$medianRows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII

Write-Host "PBE TABI-vs-grid compare COMPLETE: $outCsv"
Write-Host "PBE TABI-vs-grid raw repeats: $rawCsv"
