param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_vs_grid_compare'),
    [string[]]$Fixtures = @('sphere','two-atoms','arg','mini-protein','charged-peptide','compact-mini-protein'),
    [ValidateSet('default','barnase','paper2013','all')]
    [string]$FixtureSet = 'default',
    [string]$GoldenCsv = (Join-Path $PSScriptRoot 'pbe_inputs_local\pbe_delphi_golden.csv'),
    [string]$PyDelPhiRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) '.codex_tmp\pyDelPhi'),
    [int]$RepeatCount = 1,
    [ValidateSet('accuracy','perf')]
    [string]$Mode = 'accuracy',
    [switch]$ProfileMemory,
    [int]$BarnaseTabiSurfacePoints = 24
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
    'barnase_linear_salt0' = 'barnase.pqr'
    'barnase_linear_salt015' = 'barnase.pqr'
    'barnase_linear_salt2' = 'barnase.pqr'
}

$BarnaseFixtures = @('barnase_linear_salt0','barnase_linear_salt015','barnase_linear_salt2')
$PaperFixtures = @('sphere')

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

function Get-GoldenRow([string]$Fixture) {
    if (-not (Test-Path -LiteralPath $GoldenCsv)) { return $null }
    $rows = @(Import-Csv -LiteralPath $GoldenCsv | Where-Object { $_.fixture -eq $Fixture } | Select-Object -First 1)
    if ($rows.Count -eq 0) { return $null }
    return $rows[0]
}

function Read-DelPhiValueTable([string]$Path) {
    $values = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        $trim = $line.Trim()
        if ([string]::IsNullOrWhiteSpace($trim) -or $trim.StartsWith('!') -or $trim.StartsWith('#')) { continue }
        $parts = @($trim -split '\s+' | Where-Object { $_ })
        if ($parts.Count -lt 2) { continue }
        $key = $parts[0].ToUpperInvariant()
        $number = 0.0
        for ($i = $parts.Count - 1; $i -ge 1; --$i) {
            if ([double]::TryParse($parts[$i], [System.Globalization.NumberStyles]::Float, $InvariantCulture, [ref]$number)) {
                $values[$key] = $number
                break
            }
        }
    }
    return $values
}

function Convert-DelPhiPdbToPqr([string]$PdbPath, [string]$CrgPath, [string]$SizPath, [string]$OutPath) {
    $charges = Read-DelPhiValueTable $CrgPath
    $radii = Read-DelPhiValueTable $SizPath
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($line in Get-Content -LiteralPath $PdbPath) {
        if (-not ($line.StartsWith('ATOM') -or $line.StartsWith('HETATM'))) { continue }
        $atomName = $line.Substring(12, [Math]::Min(4, $line.Length - 12)).Trim().ToUpperInvariant()
        $charge = if ($charges.ContainsKey($atomName)) { $charges[$atomName] } else { 0.0 }
        $radius = if ($radii.ContainsKey($atomName)) { $radii[$atomName] } else { 1.5 }
        $prefix = if ($line.Length -ge 54) { $line.Substring(0, 54) } else { $line.PadRight(54) }
        $element = if ($line.Length -ge 78) { $line.Substring(76, 2).Trim() } else { $atomName.Substring(0, 1) }
        $out.Add(("{0}{1,9:F4}{2,8:F4}       {3}" -f $prefix, $charge, $radius, $element))
    }
    $out.Add('END')
    Set-Content -LiteralPath $OutPath -Value $out -Encoding ASCII
}

function Resolve-ComparePqr([string]$Fixture, [string]$RunDir) {
    $pqrName = $FixtureMap[$Fixture]
    if ([string]::IsNullOrWhiteSpace($pqrName)) { throw "Unknown fixture: $Fixture" }
    $localPqr = Join-Path $InputRoot $pqrName
    if (Test-Path -LiteralPath $localPqr) { return $localPqr }
    if ($Fixture -notlike 'barnase*') { return $localPqr }
    $fixtureDir = Join-Path $PyDelPhiRoot 'pydelphi\data\test_examples\barnase'
    $directPqr = Join-Path $fixtureDir $pqrName
    if (Test-Path -LiteralPath $directPqr) { return $directPqr }
    $pdbPath = Join-Path $fixtureDir 'barnase.pdb'
    $crgPath = Join-Path $fixtureDir 'amber.crg'
    $sizPath = Join-Path $fixtureDir 'amber.siz'
    if ((Test-Path -LiteralPath $pdbPath) -and
        (Test-Path -LiteralPath $crgPath) -and
        (Test-Path -LiteralPath $sizPath)) {
        $generated = Join-Path $RunDir $pqrName
        Convert-DelPhiPdbToPqr -PdbPath $pdbPath -CrgPath $crgPath -SizPath $sizPath -OutPath $generated
        return $generated
    }
    return $localPqr
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

function New-ArgFile([string]$RunDir, [string]$Backend, [string]$PqrPath, [string]$RunMode, [string]$Fixture) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $gold = Get-GoldenRow $Fixture
    $salt = Get-RowValue $gold 'salt_M' '0.15'
    $scale = if ($Fixture -like 'barnase*') { Get-RowValue $gold 'scale' '2.0' } else { '0.35' }
    $gridSize = if ($Fixture -like 'barnase*') { Get-RowValue $gold 'grid_size' '' } else { '' }
    $surface = if ($Fixture -like 'barnase*') { 'vdw' } else { 'vdw' }
    $margin = if ($Fixture -like 'barnase*') { '0.0' } else { '2.0' }
    $maxGrid = if ($Fixture -like 'barnase*') { '4000000' } else { '300000' }
    $maxIter = if ($Fixture -like 'barnase*') { '1200' } else { '800' }
    $rms = if ($Fixture -like 'barnase*') { '5e-4' } else { '1e-4' }
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        "-pbe=$Backend",
        "-pbe-test-pqr=$(Convert-ToTitanPath $PqrPath)",
        "-pbe-surface=$surface",
        '-pbe-probe=1.4',
        "-pbe-scale=$scale",
        "-pbe-margin=$margin",
        "-pbe-max-grid=$maxGrid",
        '-pbe-map=off',
        "-pbe-max-iter=$maxIter",
        "-pbe-rms=$rms",
        '-pbe-dphi=1e-4',
        '-pbe-timing=summary',
        "-pbe-salt=$salt",
        '-nt=1'
    )
    if (-not [string]::IsNullOrWhiteSpace($gridSize)) {
        $argLines += "-pbe-test-grid-size=$gridSize"
    }
    if ($Backend -eq 'tabi') {
        $argLines += '-pbe-tabi-derivatives=analytic'
        if ($Fixture -like 'barnase*') {
            $argLines += '-pbe-tabi-sum=treecode'
            $argLines += '-pbe-tabi-formulation=jcp2013-treecode'
            $argLines += "-pbe-tabi-surface-points=$BarnaseTabiSurfacePoints"
            $argLines += '-pbe-tabi-energy-projection=boundary'
        } elseif ($RunMode -eq 'accuracy') {
            $argLines += '-pbe-tabi-sum=direct'
            $argLines += '-pbe-tabi-formulation=jcp2013'
            $argLines += '-pbe-tabi-energy-projection=boundary'
        } else {
            $argLines += '-pbe-tabi-sum=treecode'
            $argLines += '-pbe-tabi-formulation=jcp2013-treecode'
            $argLines += '-pbe-tabi-energy-projection=boundary'
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Invoke-PbeCase([int]$Repeat, [string]$Fixture, [string]$Backend) {
    $runDir = Join-Path $WorkRoot ("repeat_{0}\{1}_{2}" -f $Repeat, $Fixture, $Backend)
    New-CleanDir $runDir
    $pqrPath = Resolve-ComparePqr -Fixture $Fixture -RunDir $runDir
    if (-not (Test-Path -LiteralPath $pqrPath)) {
        throw "Missing PQR fixture: $pqrPath"
    }
    $argPath = New-ArgFile -RunDir $runDir -Backend $Backend -PqrPath $pqrPath -RunMode $Mode -Fixture $Fixture
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
    $memory = $null
    $memoryPath = Join-Path $runDir 'pbe_memory_profile.csv'
    if ($ProfileMemory) {
        if (-not (Test-Path -LiteralPath $memoryPath)) {
            throw "Missing PBE memory profile for ${Fixture}/${Backend}: $memoryPath"
        }
        $memory = @(Import-Csv -LiteralPath $memoryPath)[0]
    }
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
        derivative_mode = if ($Backend -eq 'tabi') { Get-RowValue $audit 'derivative_mode' (Get-RowValue $timing 'tabi_derivative_mode' 'NA') } else { 'finite_difference' }
        near_integration = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_integration' (Get-RowValue $timing 'tabi_near_integration' 'NA') } else { 'grid' }
        near_kernel_policy = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_kernel_policy' (Get-RowValue $timing 'tabi_near_kernel_policy' 'NA') } else { 'grid' }
        near_k1_policy = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_k1_policy' (Get-RowValue $timing 'tabi_near_k1_policy' 'NA') } else { 'grid' }
        near_k2_policy = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_k2_policy' (Get-RowValue $timing 'tabi_near_k2_policy' 'NA') } else { 'grid' }
        near_k3_policy = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_k3_policy' (Get-RowValue $timing 'tabi_near_k3_policy' 'NA') } else { 'grid' }
        near_k4_policy = if ($Backend -eq 'tabi') { Get-RowValue $audit 'near_k4_policy' (Get-RowValue $timing 'tabi_near_k4_policy' 'NA') } else { 'grid' }
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
        tabi_singular_near_interactions = Get-RowValue $timing 'tabi_singular_near_interactions' 'NA'
        tabi_leaf_near_interactions = Get-RowValue $timing 'tabi_leaf_near_interactions' 'NA'
        tabi_far_node_interactions = Get-RowValue $timing 'tabi_far_node_interactions' 'NA'
        tabi_leaf_panel_interactions = Get-RowValue $timing 'tabi_leaf_panel_interactions' 'NA'
        tabi_leaf_tile_interactions = Get-RowValue $timing 'tabi_leaf_tile_interactions' 'NA'
        tabi_leaf_tile_panel_interactions = Get-RowValue $timing 'tabi_leaf_tile_panel_interactions' 'NA'
        tabi_leaf_residual_panel_interactions = Get-RowValue $timing 'tabi_leaf_residual_panel_interactions' 'NA'
        tabi_smooth_leaf_node_interactions = Get-RowValue $timing 'tabi_smooth_leaf_node_interactions' 'NA'
        tabi_leaf_pruned_interactions = Get-RowValue $timing 'tabi_leaf_pruned_interactions' 'NA'
        tabi_leaf_near_policy = Get-RowValue $timing 'tabi_leaf_near_policy' 'NA'
        tabi_leaf_panel_kernel = Get-RowValue $timing 'tabi_leaf_panel_kernel' 'NA'
        tabi_leaf_tile_status = Get-RowValue $timing 'tabi_leaf_tile_status' 'NA'
        tabi_leaf_panel_eval_seconds = Get-RowValue $timing 'tabi_leaf_panel_eval_seconds' 'NA'
        tabi_leaf_tile_eval_seconds = Get-RowValue $timing 'tabi_leaf_tile_eval_seconds' 'NA'
        tabi_leaf_residual_eval_seconds = Get-RowValue $timing 'tabi_leaf_residual_eval_seconds' 'NA'
        tabi_smooth_leaf_node_eval_seconds = Get-RowValue $timing 'tabi_smooth_leaf_node_eval_seconds' 'NA'
        tabi_leaf_panel_pairs_per_second = Get-RowValue $timing 'tabi_leaf_panel_pairs_per_second' 'NA'
        tabi_singular_near_eval_seconds = Get-RowValue $timing 'tabi_singular_near_eval_seconds' 'NA'
        tabi_leaf_near_eval_seconds = Get-RowValue $timing 'tabi_leaf_near_eval_seconds' 'NA'
        tabi_near_self_pairs = Get-RowValue $timing 'tabi_near_self_pairs' 'NA'
        tabi_near_one_ring_pairs = Get-RowValue $timing 'tabi_near_one_ring_pairs' 'NA'
        tabi_near_duffy_pairs = Get-RowValue $timing 'tabi_near_duffy_pairs' 'NA'
        tabi_near_centroid_pairs = Get-RowValue $timing 'tabi_near_centroid_pairs' 'NA'
        tabi_near_fallbacks = Get-RowValue $timing 'tabi_near_fallbacks' 'NA'
        tabi_near_k1_accepted = Get-RowValue $timing 'tabi_near_k1_accepted' 'NA'
        tabi_near_k1_rejected = Get-RowValue $timing 'tabi_near_k1_rejected' 'NA'
        tabi_near_k2_accepted = Get-RowValue $timing 'tabi_near_k2_accepted' 'NA'
        tabi_near_k2_rejected = Get-RowValue $timing 'tabi_near_k2_rejected' 'NA'
        tabi_near_k3_accepted = Get-RowValue $timing 'tabi_near_k3_accepted' 'NA'
        tabi_near_k3_rejected = Get-RowValue $timing 'tabi_near_k3_rejected' 'NA'
        tabi_near_k4_accepted = Get-RowValue $timing 'tabi_near_k4_accepted' 'NA'
        tabi_near_k4_rejected = Get-RowValue $timing 'tabi_near_k4_rejected' 'NA'
        tabi_near_one_ring_k1_accepted = Get-RowValue $timing 'tabi_near_one_ring_k1_accepted' 'NA'
        tabi_near_one_ring_k1_rejected = Get-RowValue $timing 'tabi_near_one_ring_k1_rejected' 'NA'
        tabi_near_one_ring_k2_accepted = Get-RowValue $timing 'tabi_near_one_ring_k2_accepted' 'NA'
        tabi_near_one_ring_k2_rejected = Get-RowValue $timing 'tabi_near_one_ring_k2_rejected' 'NA'
        tabi_near_one_ring_k3_accepted = Get-RowValue $timing 'tabi_near_one_ring_k3_accepted' 'NA'
        tabi_near_one_ring_k3_rejected = Get-RowValue $timing 'tabi_near_one_ring_k3_rejected' 'NA'
        tabi_near_one_ring_k4_accepted = Get-RowValue $timing 'tabi_near_one_ring_k4_accepted' 'NA'
        tabi_near_one_ring_k4_rejected = Get-RowValue $timing 'tabi_near_one_ring_k4_rejected' 'NA'
        tabi_near_integration_seconds = Get-RowValue $timing 'tabi_near_integration_seconds' 'NA'
        tabi_far_interactions = Get-RowValue $timing 'tabi_far_interactions' 'NA'
        tabi_status = Get-RowValue $timing 'tabi_status' 'NA'
        tabi_convergence_reason = Get-RowValue $timing 'tabi_convergence_reason' 'NA'
        reaction_energy_source = $energy.reaction_energy_source
        memory_estimated_bytes_total = Get-RowValue $memory 'estimated_bytes_total' 'NA'
        memory_estimated_bytes_surface = Get-RowValue $memory 'estimated_bytes_surface' 'NA'
        memory_estimated_bytes_solver = Get-RowValue $memory 'estimated_bytes_solver' 'NA'
        memory_estimated_bytes_operator = Get-RowValue $memory 'estimated_bytes_operator' 'NA'
        memory_estimated_bytes_tree = Get-RowValue $memory 'estimated_bytes_tree' 'NA'
        memory_rss_before_bytes = Get-RowValue $memory 'rss_before_bytes' 'NA'
        memory_rss_after_bytes = Get-RowValue $memory 'rss_after_bytes' 'NA'
        memory_rss_delta_bytes = Get-RowValue $memory 'rss_delta_bytes' 'NA'
        memory_rss_delta_signed_bytes = Get-RowValue $memory 'rss_delta_signed_bytes' 'NA'
        memory_peak_working_set_bytes = Get-RowValue $memory 'peak_working_set_bytes' 'NA'
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
        'tabi_singular_near_interactions','tabi_leaf_near_interactions',
        'tabi_far_node_interactions','tabi_leaf_panel_interactions',
        'tabi_leaf_tile_interactions','tabi_leaf_tile_panel_interactions',
        'tabi_leaf_residual_panel_interactions',
        'tabi_smooth_leaf_node_interactions','tabi_leaf_pruned_interactions',
        'tabi_leaf_panel_eval_seconds',
        'tabi_leaf_tile_eval_seconds','tabi_leaf_residual_eval_seconds',
        'tabi_smooth_leaf_node_eval_seconds',
        'tabi_leaf_panel_pairs_per_second',
        'tabi_singular_near_eval_seconds',
        'tabi_leaf_near_eval_seconds',
        'tabi_near_self_pairs','tabi_near_one_ring_pairs',
        'tabi_near_duffy_pairs','tabi_near_centroid_pairs',
        'tabi_near_fallbacks','tabi_near_k1_accepted',
        'tabi_near_k1_rejected','tabi_near_k2_accepted',
        'tabi_near_k2_rejected','tabi_near_k3_accepted',
        'tabi_near_k3_rejected','tabi_near_k4_accepted',
        'tabi_near_k4_rejected','tabi_near_one_ring_k1_accepted',
        'tabi_near_one_ring_k1_rejected',
        'tabi_near_one_ring_k2_accepted',
        'tabi_near_one_ring_k2_rejected',
        'tabi_near_one_ring_k3_accepted',
        'tabi_near_one_ring_k3_rejected',
        'tabi_near_one_ring_k4_accepted',
        'tabi_near_one_ring_k4_rejected',
        'tabi_near_integration_seconds',
        'tabi_far_interactions','memory_estimated_bytes_total',
        'memory_estimated_bytes_surface','memory_estimated_bytes_solver',
        'memory_estimated_bytes_operator','memory_estimated_bytes_tree',
        'memory_rss_before_bytes','memory_rss_after_bytes',
        'memory_rss_delta_bytes','memory_rss_delta_signed_bytes',
        'memory_peak_working_set_bytes')) {
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

if ($FixtureSet -eq 'barnase') {
    $Fixtures = $BarnaseFixtures
} elseif ($FixtureSet -eq 'paper2013') {
    $Fixtures = $PaperFixtures
    if (Test-Path -LiteralPath (Join-Path $InputRoot '1a63.pqr')) {
        $Fixtures += '1a63'
    } else {
        Write-Warning "Skipping paper2013 1A63 fixture: missing test\pbe_inputs_local\1a63.pqr"
    }
} elseif ($FixtureSet -eq 'all') {
    $Fixtures = @('sphere','two-atoms','arg','mini-protein','charged-peptide','compact-mini-protein') + $BarnaseFixtures + $PaperFixtures
    if (Test-Path -LiteralPath (Join-Path $InputRoot '1a63.pqr')) {
        $Fixtures += '1a63'
    }
} elseif (($FixtureSet -eq 'default') -and (-not $PSBoundParameters.ContainsKey('Fixtures')) -and
          (Test-Path -LiteralPath (Join-Path $InputRoot '1a63.pqr'))) {
    $Fixtures += '1a63'
}
$Fixtures = @($Fixtures | ForEach-Object {
    ([string]$_ -split ',') | ForEach-Object { $_.Trim() } | Where-Object { $_ }
} | Select-Object -Unique)
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
