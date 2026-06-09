param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_perf_matrix'),
    [string[]]$Fixtures = @('sphere','mini-protein','charged-peptide','compact-mini-protein'),
    [ValidateSet('default','barnase','paper2013','all')]
    [string]$FixtureSet = 'default',
    [string[]]$Modes = @('grid','tabi_direct','tabi_treecode'),
    [int[]]$SurfacePointCounts = @(24,48,96),
    [int]$MaxDirectSurfacePoints = 48,
    [int[]]$TreeLeafSizes = @(32),
    [double[]]$OpeningThetas = @(0.45),
    [string[]]$LeafNearPolicies = @('auto'),
    [string[]]$LeafPanelKernels = @('auto'),
    [string[]]$AdaptiveLevels = @('strict'),
    [string[]]$ClusterPolicies = @('off'),
    [string[]]$InteractionModes = @('cached'),
    [string[]]$Threads = @('1'),
    [ValidateSet('summary','detail')]
    [string]$TimingMode = 'summary',
    [int]$TabiProgressIntervalSeconds = 5,
    [int]$GmresRestart = 32,
    [switch]$TuneTreecode,
    [switch]$ProfileMemory,
    [int]$CaseTimeoutSeconds = 0,
    [string]$NearIntegrationMode = 'centroid',
    [string]$Preconditioner = 'block',
    [switch]$AllowFailedStatus
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'
$PyDelPhiRoot = Join-Path $RepoRoot '.codex_tmp\pyDelPhi'

$FixtureMap = @{
    'sphere' = 'sphere.pqr'
    'two-atoms' = 'two-atoms.pqr'
    'arg' = 'arg.pqr'
    'mini-protein' = 'mini-protein.pqr'
    'charged-peptide' = 'charged-peptide.pqr'
    'compact-mini-protein' = 'compact-mini-protein.pqr'
    'research-medium' = 'rec_027.pqr'
    'research-large' = 'rec_0.pqr'
    '1a63' = '1a63.pqr'
    'barnase' = 'barnase.pqr'
}

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

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

function Format-Number([object]$Value) {
    if ($null -eq $Value) { return 'NA' }
    if ($Value -is [string]) { return $Value }
    return ([double]$Value).ToString('G17', $InvariantCulture)
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI perf matrix root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

function Read-DelPhiValueTable([string]$Path) {
    $values = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        $trim = $line.Trim()
        if ([string]::IsNullOrWhiteSpace($trim) -or $trim.StartsWith('!') -or $trim.StartsWith('#')) { continue }
        $parts = @($trim -split '\s+' | Where-Object { $_ })
        if ($parts.Count -lt 2) { continue }
        $number = 0.0
        for ($i = $parts.Count - 1; $i -ge 1; --$i) {
            if ([double]::TryParse($parts[$i], [System.Globalization.NumberStyles]::Float, $InvariantCulture, [ref]$number)) {
                $values[$parts[0].ToUpperInvariant()] = $number
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

function Convert-ResearchPqr([string]$SourcePath, [string]$RunDir, [string]$OutputName) {
    if (-not (Test-Path -LiteralPath $SourcePath)) { return $SourcePath }
    $outPath = Join-Path $RunDir $OutputName
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($line in Get-Content -LiteralPath $SourcePath) {
        $trim = $line.Trim()
        if ($trim.StartsWith('ATOM') -or $trim.StartsWith('HETATM')) {
            $tokens = @($trim -split '\s+' | Where-Object { $_ })
            if ($tokens.Count -ge 10) {
                $radius = 0.0
                if ([double]::TryParse($tokens[$tokens.Count - 1],
                        [System.Globalization.NumberStyles]::Float,
                        $InvariantCulture, [ref]$radius) -and $radius -le 0.0) {
                    $tokens[$tokens.Count - 1] = '0.600'
                    $out.Add(($tokens -join ' '))
                    continue
                }
            }
        }
        $out.Add($line)
    }
    Set-Content -LiteralPath $outPath -Value $out -Encoding ASCII
    return $outPath
}

function Resolve-PqrFixture([string]$Fixture, [string]$RunDir) {
    if (-not $FixtureMap.ContainsKey($Fixture)) { throw "Unknown fixture: $Fixture" }
    $pqrName = $FixtureMap[$Fixture]
    $localPqr = Join-Path $InputRoot $pqrName
    if (Test-Path -LiteralPath $localPqr) { return $localPqr }
    if ($Fixture -eq 'research-medium') {
        return Convert-ResearchPqr -SourcePath (Join-Path $RepoRoot 'example\1\rec_027.pqr') -RunDir $RunDir -OutputName 'research-medium_normalized.pqr'
    }
    if ($Fixture -eq 'research-large') {
        return Convert-ResearchPqr -SourcePath (Join-Path $RepoRoot 'example\4\rec_0.pqr') -RunDir $RunDir -OutputName 'research-large_normalized.pqr'
    }
    if ($Fixture -eq 'barnase') {
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
    }
    return $null
}

function New-ArgFile([string]$RunDir, [string]$Fixture, [string]$Mode,
                     [string]$PqrPath, [int]$SurfacePoints,
                     [int]$LeafSize, [double]$OpeningTheta,
                     [string]$LeafNearPolicy,
                     [string]$LeafPanelKernel,
                     [string]$AdaptiveLevel,
                     [string]$ClusterPolicy,
                     [string]$InteractionMode,
                     [int]$ThreadCount) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $isBarnase = $Fixture -eq 'barnase'
    $isResearch = $Fixture -eq 'research-medium' -or $Fixture -eq 'research-large'
    $scale = if ($isBarnase -or $isResearch) { '0.35' } else { '0.35' }
    $margin = if ($isBarnase) { '0.0' } elseif ($isResearch) { '2.0' } else { '2.0' }
    $maxGrid = if ($isBarnase -or $isResearch) { '4000000' } else { '300000' }
    $maxIter = if ($isBarnase -or $isResearch) { '1200' } else { '800' }
    $rms = if ($isBarnase -or $isResearch) { '5e-4' } else { '1e-4' }
    $backend = if ($Mode -eq 'grid') { 'grid' } else { 'tabi' }
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        "-pbe=$backend",
        "-pbe-test-pqr=$(Convert-ToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-probe=1.4',
        "-pbe-scale=$scale",
        "-pbe-margin=$margin",
        "-pbe-max-grid=$maxGrid",
        '-pbe-map=off',
        "-pbe-max-iter=$maxIter",
        "-pbe-rms=$rms",
        '-pbe-dphi=1e-4',
        "-pbe-timing=$TimingMode",
        '-pbe-salt=0.15',
        '-pbe-temp=298.15',
        "-nt=$ThreadCount"
    )
    if ($Mode -eq 'tabi_direct') {
        $argLines += '-pbe-tabi-sum=direct'
        $argLines += '-pbe-tabi-formulation=jcp2013'
        $argLines += '-pbe-tabi-precond=diagonal'
    } elseif ($Mode -eq 'tabi_treecode') {
        $argLines += '-pbe-tabi-sum=treecode'
        $argLines += '-pbe-tabi-formulation=jcp2013-treecode'
        $argLines += "-pbe-tabi-precond=$Preconditioner"
        $argLines += "-pbe-tabi-tree-leaf-size=$LeafSize"
        $argLines += "-pbe-tabi-opening-theta=$($OpeningTheta.ToString('G17', $InvariantCulture))"
        $argLines += "-pbe-tabi-leaf-near-policy=$LeafNearPolicy"
        $argLines += "-pbe-tabi-leaf-panel-kernel=$LeafPanelKernel"
        $argLines += "-pbe-tabi-adaptive-level=$AdaptiveLevel"
        $argLines += "-pbe-tabi-cluster-policy=$ClusterPolicy"
        $argLines += "-pbe-tabi-interaction-mode=$InteractionMode"
    }
    if ($Mode -like 'tabi*') {
        $argLines += '-pbe-tabi-derivatives=analytic'
        $argLines += '-pbe-tabi-energy-projection=boundary'
        $argLines += "-pbe-tabi-near-integration=$NearIntegrationMode"
        $argLines += "-pbe-tabi-surface-points=$SurfacePoints"
        $argLines += "-pbe-tabi-gmres-restart=$GmresRestart"
        if ($TimingMode -eq 'detail') {
            $argLines += '-pbe-tabi-progress=auto'
            $argLines += "-pbe-tabi-progress-interval=$TabiProgressIntervalSeconds"
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function New-SkippedRow([string]$Fixture, [string]$Mode, [int]$SurfacePoints,
                        [int]$LeafSize, [double]$OpeningTheta,
                        [string]$LeafNearPolicy, [string]$LeafPanelKernel,
                        [string]$AdaptiveLevel,
                        [string]$ClusterPolicy,
                        [string]$InteractionMode,
                        [int]$ThreadCount,
                        [string]$Reason) {
    return [PSCustomObject]@{
        fixture = $Fixture; mode = $Mode; surface_points = $SurfacePoints
        threads = $ThreadCount
        tree_leaf_size = $LeafSize; opening_theta = Format-Number $OpeningTheta
        leaf_near_policy = if ($Mode -eq 'tabi_treecode') { $LeafNearPolicy } else { 'NA' }
        requested_leaf_panel_kernel = if ($Mode -eq 'tabi_treecode') { $LeafPanelKernel } else { 'NA' }
        adaptive_level = if ($Mode -eq 'tabi_treecode') { $AdaptiveLevel } else { 'NA' }
        cluster_policy = if ($Mode -eq 'tabi_treecode') { $ClusterPolicy } else { 'NA' }
        requested_interaction_mode = if ($Mode -eq 'tabi_treecode') { $InteractionMode } else { 'NA' }
        interaction_mode = 'NA'
        atom_count = 'NA'; grid_points = 'NA'; panels = 'NA'
        surface_point_count_requested = if ($Mode -like 'tabi*') { $SurfacePoints } else { 'NA' }
        panel_density_per_A2 = 'NA'; panel_budget_status = 'NA'
        energy_kJmol = 'NA'
        total_seconds = 'NA'; peak_ram_bytes = 'NA'; estimated_bytes = 'NA'
        gmres_iterations = 'NA'; matvec_calls = 'NA'
        gmres_restart_requested = if ($Mode -like 'tabi*') { $GmresRestart } else { 'NA' }
        gmres_restart_effective = 'NA'
        best_residual = 'NA'; best_residual_iter = 'NA'; final_residual = 'NA'
        near_interactions = 'NA'; far_interactions = 'NA'
        singular_near_pairs = 'NA'; leaf_near_pairs = 'NA'; far_node_pairs = 'NA'
        leaf_panel_interactions = 'NA'; smooth_leaf_node_interactions = 'NA'
        leaf_tile_interactions = 'NA'; leaf_tile_panel_interactions = 'NA'
        leaf_residual_panel_interactions = 'NA'
        leaf_pruned_interactions = 'NA'
        adaptive_smooth_leaf_accepted = 'NA'
        adaptive_smooth_leaf_rejected = 'NA'
        adaptive_span_tile_accepted = 'NA'
        adaptive_span_tile_rejected = 'NA'
        adaptive_exact_leaf_fallbacks = 'NA'
        adaptive_reject_one_ring_overlap = 'NA'
        adaptive_reject_error_bound_failed = 'NA'
        adaptive_reject_p4_check_failed = 'NA'
        adaptive_reject_too_close = 'NA'
        leaf_work_reduction_ratio = 'NA'
        target_cluster_interactions = 'NA'
        target_cluster_targets = 'NA'
        target_cluster_panels_pruned = 'NA'
        target_cluster_rejected = 'NA'
        target_cluster_eval_seconds = 'NA'
        singular_near_eval_seconds = 'NA'; leaf_near_eval_seconds = 'NA'
        leaf_panel_eval_seconds = 'NA'; smooth_leaf_node_eval_seconds = 'NA'
        leaf_tile_eval_seconds = 'NA'; leaf_residual_eval_seconds = 'NA'
        leaf_panel_kernel = 'NA'; leaf_tile_status = 'NA'
        leaf_panel_pairs_per_second = 'NA'
        smooth_node_eval_limited = 'false'
        tile_fast = 'false'; tile_math_limited = 'false'
        tile_memory_limited = 'false'
        tile_avx2_fast = 'false'; tile_avx2_math_limited = 'false'
        tile_avx2_memory_limited = 'false'
        near_per_matvec = 'NA'; far_per_matvec = 'NA'
        singular_near_per_matvec = 'NA'; leaf_near_per_matvec = 'NA'
        leaf_panel_per_matvec = 'NA'; smooth_leaf_node_per_matvec = 'NA'
        leaf_tile_per_matvec = 'NA'; leaf_tile_panel_per_matvec = 'NA'
        leaf_residual_panel_per_matvec = 'NA'
        leaf_pruned_per_matvec = 'NA'
        seconds_per_million_near = 'NA'; time_ratio_vs_grid = 'NA'
        speedup_vs_nt1 = 'NA'; parallel_efficiency = 'NA'
        serial_bottleneck_class = 'NA'
        ram_ratio_vs_grid = 'NA'; status = 'skipped'
        performance_classification = $Reason; run_dir = 'NA'
    }
}

function New-TimeoutRow([string]$Fixture, [string]$Mode, [int]$SurfacePoints,
                        [int]$LeafSize, [double]$OpeningTheta,
                        [string]$LeafNearPolicy, [string]$LeafPanelKernel,
                        [string]$AdaptiveLevel,
                        [string]$ClusterPolicy,
                        [string]$InteractionMode,
                        [int]$ThreadCount,
                        [string]$RunDir, [double]$ElapsedSeconds) {
    $row = New-SkippedRow -Fixture $Fixture -Mode $Mode -SurfacePoints $SurfacePoints `
        -LeafSize $LeafSize -OpeningTheta $OpeningTheta -LeafNearPolicy $LeafNearPolicy `
        -LeafPanelKernel $LeafPanelKernel -AdaptiveLevel $AdaptiveLevel `
        -ClusterPolicy $ClusterPolicy -InteractionMode $InteractionMode `
        -ThreadCount $ThreadCount `
        -Reason 'timeout_no_profile'
    $row.status = 'timeout'
    $row.total_seconds = Format-Number $ElapsedSeconds
    $row.run_dir = $RunDir
    if ($Fixture -eq 'barnase' -and $Mode -eq 'tabi_treecode' -and
        $ElapsedSeconds -lt 300.0) {
        $row.performance_classification = 'barnase_probe_timeout'
    }
    $timingPath = Join-Path $RunDir 'pbe_timing.csv'
    if (Test-Path -LiteralPath $timingPath) {
        $timing = @(Import-Csv -LiteralPath $timingPath | Select-Object -First 1)
        if ($timing.Count -gt 0) {
            $row.performance_classification =
                if ($Fixture -eq 'barnase' -and $Mode -eq 'tabi_treecode' -and
                    $ElapsedSeconds -lt 300.0) {
                    'barnase_probe_timeout'
                } else {
                    'timeout_operator_limited'
                }
            $row.atom_count = Get-RowValue $timing[0] 'atom_count' 'NA'
            $row.panels = Get-RowValue $timing[0] 'tabi_surface_panels' 'NA'
            $row.surface_point_count_requested = Get-RowValue $timing[0] 'tabi_surface_point_count_requested' $SurfacePoints
            $row.panel_density_per_A2 = Get-RowValue $timing[0] 'tabi_panel_density_per_A2' 'NA'
            $row.panel_budget_status = Get-RowValue $timing[0] 'tabi_panel_budget_status' 'NA'
            $row.matvec_calls = Get-RowValue $timing[0] 'tabi_matvec_calls' 'NA'
            $row.gmres_restart_requested = Get-RowValue $timing[0] 'tabi_gmres_restart_requested' $GmresRestart
            $row.gmres_restart_effective = Get-RowValue $timing[0] 'tabi_gmres_restart_effective' 'NA'
            $row.best_residual = Get-RowValue $timing[0] 'tabi_best_residual' 'NA'
            $row.best_residual_iter = Get-RowValue $timing[0] 'tabi_best_residual_iter' 'NA'
            $row.final_residual = Get-RowValue $timing[0] 'tabi_final_residual' 'NA'
            $row.near_interactions = Get-RowValue $timing[0] 'tabi_near_interactions' 'NA'
            $row.far_interactions = Get-RowValue $timing[0] 'tabi_far_interactions' 'NA'
            $row.leaf_panel_interactions = Get-RowValue $timing[0] 'tabi_leaf_panel_interactions' 'NA'
            $row.leaf_tile_interactions = Get-RowValue $timing[0] 'tabi_leaf_tile_interactions' 'NA'
            $row.leaf_tile_panel_interactions = Get-RowValue $timing[0] 'tabi_leaf_tile_panel_interactions' 'NA'
            $row.leaf_residual_panel_interactions = Get-RowValue $timing[0] 'tabi_leaf_residual_panel_interactions' 'NA'
            $row.adaptive_smooth_leaf_accepted = Get-RowValue $timing[0] 'tabi_adaptive_smooth_leaf_accepted' 'NA'
            $row.adaptive_smooth_leaf_rejected = Get-RowValue $timing[0] 'tabi_adaptive_smooth_leaf_rejected' 'NA'
            $row.adaptive_span_tile_accepted = Get-RowValue $timing[0] 'tabi_adaptive_span_tile_accepted' 'NA'
            $row.adaptive_span_tile_rejected = Get-RowValue $timing[0] 'tabi_adaptive_span_tile_rejected' 'NA'
            $row.adaptive_exact_leaf_fallbacks = Get-RowValue $timing[0] 'tabi_adaptive_exact_leaf_fallbacks' 'NA'
            $row.adaptive_reject_one_ring_overlap = Get-RowValue $timing[0] 'tabi_adaptive_reject_one_ring_overlap' 'NA'
            $row.adaptive_reject_error_bound_failed = Get-RowValue $timing[0] 'tabi_adaptive_reject_error_bound_failed' 'NA'
            $row.adaptive_reject_p4_check_failed = Get-RowValue $timing[0] 'tabi_adaptive_reject_p4_check_failed' 'NA'
            $row.adaptive_reject_too_close = Get-RowValue $timing[0] 'tabi_adaptive_reject_too_close' 'NA'
            $row.leaf_work_reduction_ratio = Get-RowValue $timing[0] 'tabi_leaf_work_reduction_ratio' 'NA'
            $row.target_cluster_interactions = Get-RowValue $timing[0] 'tabi_target_cluster_interactions' 'NA'
            $row.target_cluster_targets = Get-RowValue $timing[0] 'tabi_target_cluster_targets' 'NA'
            $row.target_cluster_panels_pruned = Get-RowValue $timing[0] 'tabi_target_cluster_panels_pruned' 'NA'
            $row.target_cluster_rejected = Get-RowValue $timing[0] 'tabi_target_cluster_rejected' 'NA'
            $row.target_cluster_eval_seconds = Get-RowValue $timing[0] 'tabi_target_cluster_eval_seconds' 'NA'
            $row.adaptive_level = Get-RowValue $timing[0] 'tabi_adaptive_level' $AdaptiveLevel
            $row.interaction_mode = Get-RowValue $timing[0] 'tabi_interaction_mode' $InteractionMode
            $row.leaf_panel_eval_seconds = Get-RowValue $timing[0] 'tabi_leaf_panel_eval_seconds' 'NA'
            $row.leaf_tile_eval_seconds = Get-RowValue $timing[0] 'tabi_leaf_tile_eval_seconds' 'NA'
            $row.leaf_residual_eval_seconds = Get-RowValue $timing[0] 'tabi_leaf_residual_eval_seconds' 'NA'
            $row.leaf_panel_kernel = Get-RowValue $timing[0] 'tabi_leaf_panel_kernel' 'NA'
            $row.leaf_tile_status = Get-RowValue $timing[0] 'tabi_leaf_tile_status' 'NA'
            $row.leaf_panel_pairs_per_second = Get-RowValue $timing[0] 'tabi_leaf_panel_pairs_per_second' 'NA'
        }
    }
    return $row
}

function Convert-ToProcessArgument([string]$Value) {
    if ($null -eq $Value) { return '""' }
    $escaped = $Value -replace '"', '\"'
    return '"' + $escaped + '"'
}

function Invoke-TitanWithTimeout([string]$ArgPath, [string]$RunDir, [int]$TimeoutSeconds) {
    $psi = [System.Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = $TitanExe
    $psi.Arguments = '-a ' + (Convert-ToProcessArgument $ArgPath)
    $psi.WorkingDirectory = $RunDir
    $psi.UseShellExecute = $false
    $psi.CreateNoWindow = $true
    $process = [System.Diagnostics.Process]::new()
    $process.StartInfo = $psi
    $watch = [System.Diagnostics.Stopwatch]::StartNew()
    [void]$process.Start()
    $completed = $process.WaitForExit($TimeoutSeconds * 1000)
    if (-not $completed) {
        try { $process.Kill() } catch {}
        $process.WaitForExit()
    } else {
        $process.WaitForExit()
    }
    $watch.Stop()
    return [pscustomobject]@{
        Completed = $completed
        ExitCode = if ($completed) { $process.ExitCode } else { $null }
        ElapsedSeconds = $watch.Elapsed.TotalSeconds
    }
}

function Invoke-PerfCase([string]$Fixture, [string]$Mode, [int]$SurfacePoints,
                         [int]$LeafSize, [double]$OpeningTheta,
                         [string]$LeafNearPolicy,
                         [string]$LeafPanelKernel,
                         [string]$AdaptiveLevel,
                         [string]$ClusterPolicy,
                         [string]$InteractionMode,
                         [int]$ThreadCount) {
    if ($Mode -eq 'tabi_direct' -and $SurfacePoints -gt $MaxDirectSurfacePoints) {
        return New-SkippedRow -Fixture $Fixture -Mode $Mode -SurfacePoints $SurfacePoints `
            -LeafSize $LeafSize -OpeningTheta $OpeningTheta -LeafNearPolicy $LeafNearPolicy `
            -LeafPanelKernel $LeafPanelKernel -AdaptiveLevel $AdaptiveLevel `
            -ClusterPolicy $ClusterPolicy -InteractionMode $InteractionMode `
            -ThreadCount $ThreadCount `
            -Reason 'skipped_slow_direct_guard'
    }
    $safeTheta = $OpeningTheta.ToString('G17', $InvariantCulture) -replace '\.', 'p'
    $safePolicy = ($LeafNearPolicy -replace '[^A-Za-z0-9_-]', '_')
    $safeKernel = ($LeafPanelKernel -replace '[^A-Za-z0-9_-]', '_')
    $safeAdaptive = ($AdaptiveLevel -replace '[^A-Za-z0-9_-]', '_')
    $safeCluster = ($ClusterPolicy -replace '[^A-Za-z0-9_-]', '_')
    $safeInteraction = ($InteractionMode -replace '[^A-Za-z0-9_-]', '_')
    $runName = "{0}_{1}_sp{2}_leaf{3}_theta{4}_lnp{5}_lpk{6}_adl{7}_clp{8}_im{9}_nt{10}" -f $Fixture,$Mode,$SurfacePoints,$LeafSize,$safeTheta,$safePolicy,$safeKernel,$safeAdaptive,$safeCluster,$safeInteraction,$ThreadCount
    $runDir = Join-Path $WorkRoot $runName
    New-CleanDir $runDir
    $pqrPath = Resolve-PqrFixture -Fixture $Fixture -RunDir $runDir
    if ($null -eq $pqrPath -or -not (Test-Path -LiteralPath $pqrPath)) {
        return New-SkippedRow -Fixture $Fixture -Mode $Mode -SurfacePoints $SurfacePoints `
            -LeafSize $LeafSize -OpeningTheta $OpeningTheta -LeafNearPolicy $LeafNearPolicy `
            -LeafPanelKernel $LeafPanelKernel -AdaptiveLevel $AdaptiveLevel `
            -ClusterPolicy $ClusterPolicy -InteractionMode $InteractionMode `
            -ThreadCount $ThreadCount `
            -Reason 'skipped_missing_fixture'
    }
    $argPath = New-ArgFile -RunDir $runDir -Fixture $Fixture -Mode $Mode `
        -PqrPath $pqrPath -SurfacePoints $SurfacePoints -LeafSize $LeafSize `
        -OpeningTheta $OpeningTheta -LeafNearPolicy $LeafNearPolicy `
        -LeafPanelKernel $LeafPanelKernel -AdaptiveLevel $AdaptiveLevel `
        -ClusterPolicy $ClusterPolicy -InteractionMode $InteractionMode `
        -ThreadCount $ThreadCount
    $elapsedSeconds = 0.0
    if ($CaseTimeoutSeconds -gt 0) {
        $runResult = Invoke-TitanWithTimeout -ArgPath $argPath -RunDir $runDir -TimeoutSeconds $CaseTimeoutSeconds
        if (-not $runResult.Completed) {
            return New-TimeoutRow -Fixture $Fixture -Mode $Mode -SurfacePoints $SurfacePoints `
                -LeafSize $LeafSize -OpeningTheta $OpeningTheta -LeafNearPolicy $LeafNearPolicy `
                -LeafPanelKernel $LeafPanelKernel -AdaptiveLevel $AdaptiveLevel `
                -ClusterPolicy $ClusterPolicy -InteractionMode $InteractionMode `
                -ThreadCount $ThreadCount `
                -RunDir $runDir `
                -ElapsedSeconds $runResult.ElapsedSeconds
        }
        $elapsedSeconds = $runResult.ElapsedSeconds
        if ($runResult.ExitCode -ne 0) {
            throw "TITAN exited with code $($runResult.ExitCode) for $runName"
        }
    } else {
        $measure = Measure-Command {
            Push-Location $runDir
            try {
                & $TitanExe -a $argPath | Out-Null
                if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE for $runName" }
            } finally {
                Pop-Location
            }
        }
        $elapsedSeconds = $measure.TotalSeconds
    }
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    $timing = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_timing.csv'))[0]
    $memoryPath = Join-Path $runDir 'pbe_memory_profile.csv'
    if ($ProfileMemory -and -not (Test-Path -LiteralPath $memoryPath)) {
        throw "Missing memory profile for $runName"
    }
    $memory = if (Test-Path -LiteralPath $memoryPath) { @(Import-Csv -LiteralPath $memoryPath)[0] } else { $null }
    $auditPath = Join-Path $runDir 'pbe_tabi_formulation_audit.csv'
    $audit = if (Test-Path -LiteralPath $auditPath) { @(Import-Csv -LiteralPath $auditPath)[0] } else { $null }
    $matvecCalls = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_matvec_calls' '0')
    $near = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_near_interactions' '0')
    $far = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_far_interactions' '0')
    $selfPairs = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_near_self_pairs' '0')
    $oneRingPairs = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_near_one_ring_pairs' '0')
    $emittedSingular = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_singular_near_interactions' 'NA')
    $emittedLeafNear = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_near_interactions' 'NA')
    $farNodePairs = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_far_node_interactions' (Get-RowValue $timing 'tabi_far_interactions' '0'))
    $leafPanelPairs = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_panel_interactions' 'NA')
    $leafTileNodes = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_tile_interactions' 'NA')
    $leafTilePanels = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_tile_panel_interactions' 'NA')
    $leafResidualPanels = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_residual_panel_interactions' 'NA')
    $smoothLeafNodes = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_smooth_leaf_node_interactions' 'NA')
    $leafPrunedPairs = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_pruned_interactions' 'NA')
    $adaptiveSmoothAccepted = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_adaptive_smooth_leaf_accepted' 'NA')
    $adaptiveSmoothRejected = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_adaptive_smooth_leaf_rejected' 'NA')
    $adaptiveSpanAccepted = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_adaptive_span_tile_accepted' 'NA')
    $adaptiveSpanRejected = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_adaptive_span_tile_rejected' 'NA')
    $adaptiveExactFallbacks = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_adaptive_exact_leaf_fallbacks' 'NA')
    $leafWorkReductionRatio = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_work_reduction_ratio' 'NA')
    $singular = if ($null -ne $emittedSingular) { $emittedSingular } else { $selfPairs + $oneRingPairs }
    $leafNear = if ($null -ne $emittedLeafNear) { $emittedLeafNear } else { [math]::Max(0.0, $near - $singular) }
    $leafPanel = if ($null -ne $leafPanelPairs) { $leafPanelPairs } else { $leafNear }
    $leafTileNodeCount = if ($null -ne $leafTileNodes) { $leafTileNodes } else { 0.0 }
    $leafTilePanelCount = if ($null -ne $leafTilePanels) { $leafTilePanels } else { 0.0 }
    $leafResidualPanelCount = if ($null -ne $leafResidualPanels) { $leafResidualPanels } else { [math]::Max(0.0, $leafPanel - $leafTilePanelCount) }
    $smoothLeaf = if ($null -ne $smoothLeafNodes) { $smoothLeafNodes } else { 0.0 }
    $leafPruned = if ($null -ne $leafPrunedPairs) { $leafPrunedPairs } else { 0.0 }
    $matvecSeconds = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_matvec_seconds_total' '0')
    $secondsPerMillionNear = if ($near -gt 0.0) { $matvecSeconds * 1000000.0 / $near } else { $null }
    $leafPanelEvalSeconds = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_panel_eval_seconds' 'NA')
    $leafTileEvalSeconds = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_tile_eval_seconds' 'NA')
    $leafResidualEvalSeconds = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_residual_eval_seconds' 'NA')
    $leafPanelPairsPerSecond = Convert-ToDoubleOrNull (Get-RowValue $timing 'tabi_leaf_panel_pairs_per_second' 'NA')
    if ($energy.status -ne 'ok' -and -not $AllowFailedStatus) {
        throw "PBE perf case failed: $runName status=$($energy.status)"
    }
    return [PSCustomObject]@{
        fixture = $Fixture
        mode = $Mode
        surface_points = if ($Mode -eq 'grid') { 0 } else { $SurfacePoints }
        threads = $ThreadCount
        tree_leaf_size = if ($Mode -eq 'tabi_treecode') { $LeafSize } else { 'NA' }
        opening_theta = if ($Mode -eq 'tabi_treecode') { Format-Number $OpeningTheta } else { 'NA' }
        leaf_near_policy = if ($Mode -eq 'tabi_treecode') { Get-RowValue $timing 'tabi_leaf_near_policy' $LeafNearPolicy } else { 'NA' }
        requested_leaf_panel_kernel = if ($Mode -eq 'tabi_treecode') { $LeafPanelKernel } else { 'NA' }
        adaptive_level = if ($Mode -eq 'tabi_treecode') { Get-RowValue $timing 'tabi_adaptive_level' $AdaptiveLevel } else { 'NA' }
        cluster_policy = if ($Mode -eq 'tabi_treecode') { $ClusterPolicy } else { 'NA' }
        requested_interaction_mode = if ($Mode -eq 'tabi_treecode') { $InteractionMode } else { 'NA' }
        interaction_mode = if ($Mode -eq 'tabi_treecode') { Get-RowValue $timing 'tabi_interaction_mode' $InteractionMode } else { 'NA' }
        atom_count = $energy.atom_count
        grid_points = $energy.grid_points
        panels = if ($Mode -like 'tabi*') { Get-RowValue $timing 'tabi_surface_panels' 'NA' } else { '0' }
        surface_point_count_requested = if ($Mode -like 'tabi*') { Get-RowValue $timing 'tabi_surface_point_count_requested' $SurfacePoints } else { 'NA' }
        panel_density_per_A2 = if ($Mode -like 'tabi*') { Get-RowValue $timing 'tabi_panel_density_per_A2' 'NA' } else { 'NA' }
        panel_budget_status = if ($Mode -like 'tabi*') { Get-RowValue $timing 'tabi_panel_budget_status' 'NA' } else { 'NA' }
        energy_kJmol = $energy.E_rxn_kJmol
        total_seconds = Get-RowValue $timing 'total_seconds' ([math]::Round($elapsedSeconds, 6).ToString($InvariantCulture))
        peak_ram_bytes = Get-RowValue $memory 'peak_working_set_bytes' 'NA'
        estimated_bytes = Get-RowValue $memory 'estimated_bytes_total' 'NA'
        gmres_iterations = if ($Mode -like 'tabi*') { Get-RowValue $audit 'iterations' (Get-RowValue $timing 'iterations' 'NA') } else { Get-RowValue $timing 'iterations' 'NA' }
        matvec_calls = Format-Number $matvecCalls
        gmres_restart_requested = if ($Mode -like 'tabi*') { Get-RowValue $timing 'tabi_gmres_restart_requested' $GmresRestart } else { 'NA' }
        gmres_restart_effective = if ($Mode -like 'tabi*') { Get-RowValue $timing 'tabi_gmres_restart_effective' 'NA' } else { 'NA' }
        best_residual = if ($Mode -like 'tabi*') { Get-RowValue $timing 'tabi_best_residual' 'NA' } else { 'NA' }
        best_residual_iter = if ($Mode -like 'tabi*') { Get-RowValue $timing 'tabi_best_residual_iter' 'NA' } else { 'NA' }
        final_residual = if ($Mode -like 'tabi*') { Get-RowValue $timing 'tabi_final_residual' 'NA' } else { 'NA' }
        near_interactions = Format-Number $near
        far_interactions = Format-Number $far
        singular_near_pairs = Format-Number $singular
        leaf_near_pairs = Format-Number $leafNear
        far_node_pairs = Format-Number $farNodePairs
        leaf_panel_interactions = Format-Number $leafPanel
        leaf_tile_interactions = Format-Number $leafTileNodeCount
        leaf_tile_panel_interactions = Format-Number $leafTilePanelCount
        leaf_residual_panel_interactions = Format-Number $leafResidualPanelCount
        smooth_leaf_node_interactions = Format-Number $smoothLeaf
        leaf_pruned_interactions = Format-Number $leafPruned
        adaptive_smooth_leaf_accepted = Format-Number $adaptiveSmoothAccepted
        adaptive_smooth_leaf_rejected = Format-Number $adaptiveSmoothRejected
        adaptive_span_tile_accepted = Format-Number $adaptiveSpanAccepted
        adaptive_span_tile_rejected = Format-Number $adaptiveSpanRejected
        adaptive_exact_leaf_fallbacks = Format-Number $adaptiveExactFallbacks
        adaptive_reject_one_ring_overlap = Get-RowValue $timing 'tabi_adaptive_reject_one_ring_overlap' 'NA'
        adaptive_reject_error_bound_failed = Get-RowValue $timing 'tabi_adaptive_reject_error_bound_failed' 'NA'
        adaptive_reject_p4_check_failed = Get-RowValue $timing 'tabi_adaptive_reject_p4_check_failed' 'NA'
        adaptive_reject_too_close = Get-RowValue $timing 'tabi_adaptive_reject_too_close' 'NA'
        leaf_work_reduction_ratio = Format-Number $leafWorkReductionRatio
        target_cluster_interactions = Get-RowValue $timing 'tabi_target_cluster_interactions' 'NA'
        target_cluster_targets = Get-RowValue $timing 'tabi_target_cluster_targets' 'NA'
        target_cluster_panels_pruned = Get-RowValue $timing 'tabi_target_cluster_panels_pruned' 'NA'
        target_cluster_rejected = Get-RowValue $timing 'tabi_target_cluster_rejected' 'NA'
        target_cluster_eval_seconds = Get-RowValue $timing 'tabi_target_cluster_eval_seconds' 'NA'
        singular_near_eval_seconds = Get-RowValue $timing 'tabi_singular_near_eval_seconds' 'NA'
        leaf_near_eval_seconds = Get-RowValue $timing 'tabi_leaf_near_eval_seconds' 'NA'
        leaf_panel_eval_seconds = Get-RowValue $timing 'tabi_leaf_panel_eval_seconds' 'NA'
        leaf_tile_eval_seconds = Format-Number $leafTileEvalSeconds
        leaf_residual_eval_seconds = Format-Number $leafResidualEvalSeconds
        smooth_leaf_node_eval_seconds = Get-RowValue $timing 'tabi_smooth_leaf_node_eval_seconds' 'NA'
        leaf_panel_kernel = Get-RowValue $timing 'tabi_leaf_panel_kernel' 'NA'
        leaf_tile_status = Get-RowValue $timing 'tabi_leaf_tile_status' 'NA'
        leaf_panel_pairs_per_second = Format-Number $leafPanelPairsPerSecond
        preconditioner_requested = if ($Mode -like 'tabi*') { Get-RowValue $audit 'preconditioner_requested' 'NA' } else { 'NA' }
        preconditioner_effective = if ($Mode -like 'tabi*') { Get-RowValue $audit 'preconditioner_effective' 'NA' } else { 'NA' }
        smooth_node_eval_limited = 'false'
        tile_fast = 'false'
        tile_math_limited = 'false'
        tile_memory_limited = 'false'
        tile_avx2_fast = 'false'
        tile_avx2_math_limited = 'false'
        tile_avx2_memory_limited = 'false'
        near_per_matvec = if ($matvecCalls -gt 0.0) { Format-Number ($near / $matvecCalls) } else { 'NA' }
        far_per_matvec = if ($matvecCalls -gt 0.0) { Format-Number ($far / $matvecCalls) } else { 'NA' }
        singular_near_per_matvec = if ($matvecCalls -gt 0.0) { Format-Number ($singular / $matvecCalls) } else { 'NA' }
        leaf_near_per_matvec = if ($matvecCalls -gt 0.0) { Format-Number ($leafNear / $matvecCalls) } else { 'NA' }
        leaf_panel_per_matvec = if ($matvecCalls -gt 0.0) { Format-Number ($leafPanel / $matvecCalls) } else { 'NA' }
        leaf_tile_per_matvec = if ($matvecCalls -gt 0.0) { Format-Number ($leafTileNodeCount / $matvecCalls) } else { 'NA' }
        leaf_tile_panel_per_matvec = if ($matvecCalls -gt 0.0) { Format-Number ($leafTilePanelCount / $matvecCalls) } else { 'NA' }
        leaf_residual_panel_per_matvec = if ($matvecCalls -gt 0.0) { Format-Number ($leafResidualPanelCount / $matvecCalls) } else { 'NA' }
        smooth_leaf_node_per_matvec = if ($matvecCalls -gt 0.0) { Format-Number ($smoothLeaf / $matvecCalls) } else { 'NA' }
        leaf_pruned_per_matvec = if ($matvecCalls -gt 0.0) { Format-Number ($leafPruned / $matvecCalls) } else { 'NA' }
        seconds_per_million_near = Format-Number $secondsPerMillionNear
        time_ratio_vs_grid = 'NA'
        speedup_vs_nt1 = 'NA'
        parallel_efficiency = 'NA'
        serial_bottleneck_class = 'NA'
        ram_ratio_vs_grid = 'NA'
        status = $energy.status
        performance_classification = 'pending'
        run_dir = $runDir
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if ($TuneTreecode) {
    $TreeLeafSizes = @(32,64,128)
    $OpeningThetas = @(0.6,0.8,1.0)
    if ($LeafNearPolicies.Count -eq 1 -and $LeafNearPolicies[0] -eq 'auto') {
        $LeafNearPolicies = @('exact','smooth-node')
    }
    if ($LeafPanelKernels.Count -eq 1 -and $LeafPanelKernels[0] -eq 'auto') {
        $LeafPanelKernels = @('scalar','tile-scalar','tile-avx2')
    }
}
if ($FixtureSet -eq 'barnase') {
    $Fixtures = @('barnase')
} elseif ($FixtureSet -eq 'paper2013') {
    $Fixtures = @('sphere')
    if (Test-Path -LiteralPath (Join-Path $InputRoot '1a63.pqr')) { $Fixtures += '1a63' }
} elseif ($FixtureSet -eq 'all') {
    $Fixtures = @('sphere','mini-protein','charged-peptide','compact-mini-protein','barnase')
    if (Test-Path -LiteralPath (Join-Path $InputRoot '1a63.pqr')) { $Fixtures += '1a63' }
}
$Fixtures = @($Fixtures | ForEach-Object {
    ([string]$_ -split ',') | ForEach-Object { $_.Trim().ToLowerInvariant() } | Where-Object { $_ }
} | Select-Object -Unique)
$maxThreads = [Math]::Max(1, [Environment]::ProcessorCount)
$Threads = @($Threads | ForEach-Object {
    ([string]$_ -split ',') | ForEach-Object {
        $text = $_.Trim()
        if ($text) { [Math]::Max(1, [Math]::Min($maxThreads, [int]$text)) }
    }
} | Select-Object -Unique | Sort-Object)
$InteractionModes = @($InteractionModes | ForEach-Object {
    ([string]$_ -split ',') | ForEach-Object {
        $text = $_.Trim().ToLowerInvariant()
        if ($text) { $text }
    }
} | Select-Object -Unique)
foreach ($interactionMode in $InteractionModes) {
    if ($interactionMode -notin @('auto','cached','onthefly')) {
        throw "Invalid -InteractionModes value '$interactionMode'. Expected auto, cached, or onthefly."
    }
}
New-CleanDir $WorkRoot

$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixture in $Fixtures) {
    foreach ($threadCount in $Threads) {
        if ($Modes -contains 'grid') {
            $rows.Add((Invoke-PerfCase -Fixture $fixture -Mode 'grid' -SurfacePoints 0 -LeafSize 0 -OpeningTheta 0.0 -LeafNearPolicy 'NA' -LeafPanelKernel 'NA' -AdaptiveLevel 'NA' -ClusterPolicy 'NA' -InteractionMode 'NA' -ThreadCount $threadCount))
        }
        foreach ($sp in $SurfacePointCounts) {
            if ($Modes -contains 'tabi_direct') {
                $rows.Add((Invoke-PerfCase -Fixture $fixture -Mode 'tabi_direct' -SurfacePoints $sp -LeafSize 0 -OpeningTheta 0.0 -LeafNearPolicy 'NA' -LeafPanelKernel 'NA' -AdaptiveLevel 'strict' -ClusterPolicy 'NA' -InteractionMode 'NA' -ThreadCount $threadCount))
            }
            if ($Modes -contains 'tabi_treecode') {
                foreach ($leaf in $TreeLeafSizes) {
                    foreach ($theta in $OpeningThetas) {
                        foreach ($policy in $LeafNearPolicies) {
                            foreach ($kernel in $LeafPanelKernels) {
                                foreach ($adaptive in $AdaptiveLevels) {
                                    foreach ($clusterPolicy in $ClusterPolicies) {
                                        foreach ($interactionMode in $InteractionModes) {
                                            $rows.Add((Invoke-PerfCase -Fixture $fixture -Mode 'tabi_treecode' -SurfacePoints $sp -LeafSize $leaf -OpeningTheta $theta -LeafNearPolicy $policy -LeafPanelKernel $kernel -AdaptiveLevel $adaptive -ClusterPolicy $clusterPolicy -InteractionMode $interactionMode -ThreadCount $threadCount))
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

foreach ($group in ($rows | Group-Object fixture)) {
    $grid = @($group.Group | Where-Object { $_.mode -eq 'grid' -and $_.status -eq 'ok' } | Select-Object -First 1)
    if ($grid.Count -eq 0) { continue }
    $gridSeconds = Convert-ToDoubleOrNull $grid[0].total_seconds
    $gridRam = Convert-ToDoubleOrNull $grid[0].peak_ram_bytes
    foreach ($row in $group.Group) {
        $seconds = Convert-ToDoubleOrNull $row.total_seconds
        $ram = Convert-ToDoubleOrNull $row.peak_ram_bytes
        if ($null -ne $seconds -and $null -ne $gridSeconds -and $gridSeconds -gt 0.0) {
            $row.time_ratio_vs_grid = Format-Number ($seconds / $gridSeconds)
        }
        if ($null -ne $ram -and $null -ne $gridRam -and $gridRam -gt 0.0) {
            $row.ram_ratio_vs_grid = Format-Number ($ram / $gridRam)
        }
        $rowThreads = Convert-ToDoubleOrNull $row.threads
        if ($row.status -eq 'ok' -and $null -ne $seconds -and
            $null -ne $rowThreads -and $rowThreads -gt 1.0) {
            $nt1 = @($group.Group | Where-Object {
                $_.status -eq 'ok' -and
                $_.mode -eq $row.mode -and
                $_.surface_points -eq $row.surface_points -and
                $_.tree_leaf_size -eq $row.tree_leaf_size -and
                $_.opening_theta -eq $row.opening_theta -and
                $_.leaf_near_policy -eq $row.leaf_near_policy -and
                $_.requested_leaf_panel_kernel -eq $row.requested_leaf_panel_kernel -and
                $_.adaptive_level -eq $row.adaptive_level -and
                $_.cluster_policy -eq $row.cluster_policy -and
                $_.threads -eq 1
            } | Select-Object -First 1)
            if ($nt1.Count -gt 0) {
                $nt1Seconds = Convert-ToDoubleOrNull $nt1[0].total_seconds
                if ($null -ne $nt1Seconds -and $seconds -gt 0.0) {
                    $speedup = $nt1Seconds / $seconds
                    $row.speedup_vs_nt1 = Format-Number $speedup
                    $row.parallel_efficiency = Format-Number ($speedup / $rowThreads)
                    $row.serial_bottleneck_class =
                        if ($speedup / $rowThreads -ge 0.65) {
                            'scales_ok'
                        } elseif ($speedup -gt 1.05) {
                            'partial_scaling'
                        } else {
                            'serial_or_overhead_limited'
                        }
                }
            }
        } elseif ($row.status -eq 'ok' -and $null -ne $rowThreads -and
                  $rowThreads -eq 1.0) {
            $row.speedup_vs_nt1 = '1'
            $row.parallel_efficiency = '1'
            $row.serial_bottleneck_class = 'nt1_baseline'
        }
        if ($row.status -eq 'skipped') { continue }
        if ($row.status -eq 'timeout') { continue }
        if ($row.mode -eq 'grid') {
            $row.performance_classification = 'grid_reference'
        } elseif ($row.mode -eq 'tabi_direct') {
            $row.performance_classification = 'direct_reference'
        } else {
            $near = Convert-ToDoubleOrNull $row.near_per_matvec
            $leafNear = Convert-ToDoubleOrNull $row.leaf_panel_per_matvec
            if ($null -eq $leafNear) { $leafNear = Convert-ToDoubleOrNull $row.leaf_near_per_matvec }
            $leafPruned = Convert-ToDoubleOrNull $row.leaf_pruned_per_matvec
            $leafPanelPairsPerSecond = Convert-ToDoubleOrNull $row.leaf_panel_pairs_per_second
            $far = Convert-ToDoubleOrNull $row.far_per_matvec
            $ratio = Convert-ToDoubleOrNull $row.time_ratio_vs_grid
            $leafTilePanels = Convert-ToDoubleOrNull $row.leaf_tile_panel_interactions
            $leafResidualPanels = Convert-ToDoubleOrNull $row.leaf_residual_panel_interactions
            $leafWorkReductionRatio = Convert-ToDoubleOrNull $row.leaf_work_reduction_ratio
            $targetClusterPanelsPruned = Convert-ToDoubleOrNull $row.target_cluster_panels_pruned
            $targetClusterRejected = Convert-ToDoubleOrNull $row.target_cluster_rejected
            $targetClusterEvalSeconds = Convert-ToDoubleOrNull $row.target_cluster_eval_seconds
            $largeProteinClass = $null
            $adaptiveClass = $null
            $barnaseClass = $null
            if ($row.leaf_near_policy -eq 'large-protein' -or
                $row.leaf_near_policy -eq 'adaptive-large') {
                $matchingExactLarge = @($group.Group | Where-Object {
                    $_.mode -eq 'tabi_treecode' -and
                    $_.surface_points -eq $row.surface_points -and
                    $_.tree_leaf_size -eq $row.tree_leaf_size -and
                    $_.opening_theta -eq $row.opening_theta -and
                    $_.leaf_panel_kernel -eq $row.leaf_panel_kernel -and
                    $_.cluster_policy -eq $row.cluster_policy -and
                    $_.leaf_near_policy -eq 'exact' -and
                    $_.status -eq 'ok'
                } | Select-Object -First 1)
                if ($matchingExactLarge.Count -gt 0) {
                    $exactSeconds = Convert-ToDoubleOrNull $matchingExactLarge[0].total_seconds
                    $exactLeafWork = Convert-ToDoubleOrNull $matchingExactLarge[0].leaf_panel_interactions
                    $exactResidualPanels = Convert-ToDoubleOrNull $matchingExactLarge[0].leaf_residual_panel_interactions
                    $exactGmres = Convert-ToDoubleOrNull $matchingExactLarge[0].gmres_iterations
                    $rowLeafWork = Convert-ToDoubleOrNull $row.leaf_panel_interactions
                    $rowGmres = Convert-ToDoubleOrNull $row.gmres_iterations
                    $timeRegressed = $false
                    $workRegressed = $false
                    $residualPruned = $false
                    $gmresRegressed = $false
                    if ($null -ne $seconds -and $null -ne $exactSeconds) {
                        $timeRegressed = $seconds -gt 1.02 * $exactSeconds
                    }
                    if ($null -ne $rowLeafWork -and $null -ne $exactLeafWork) {
                        $workRegressed = $rowLeafWork -gt 1.10 * $exactLeafWork
                    }
                    if ($null -ne $leafResidualPanels -and $null -ne $exactResidualPanels) {
                        $residualPruned = $leafResidualPanels -le ($exactResidualPanels / 3.0)
                    }
                    if ($null -ne $rowGmres -and $null -ne $exactGmres) {
                        $gmresRegressed = $rowGmres -gt $exactGmres
                    }
                    if ($row.leaf_near_policy -eq 'adaptive-large') {
                        if ($timeRegressed) {
                            $adaptiveClass =
                                if ($row.adaptive_level -eq 'balanced' -or
                                    $row.adaptive_level -eq 'aggressive') {
                                    'adaptive_level_too_aggressive'
                                } else {
                                    'adaptive_error_limited'
                                }
                        } elseif ($workRegressed) {
                            $adaptiveClass = 'span_tile_work_regressed'
                        } elseif ($null -ne $leafWorkReductionRatio -and
                                  $leafWorkReductionRatio -ge 0.25) {
                            $adaptiveClass =
                                if ($null -ne $targetClusterPanelsPruned -and
                                    $targetClusterPanelsPruned -gt 0.0) {
                                    'cluster_leaf_pruned'
                                } else {
                                    'adaptive_leaf_pruned'
                                }
                        } elseif ($null -ne $targetClusterPanelsPruned -and
                                  $targetClusterPanelsPruned -gt 0.0 -and
                                  $null -ne $targetClusterEvalSeconds -and
                                  $targetClusterEvalSeconds -gt 0.0) {
                            $adaptiveClass = 'target_cluster_error_limited'
                        } elseif ($null -ne $targetClusterRejected -and
                                  $targetClusterRejected -gt 0.0) {
                            $adaptiveClass = 'target_cluster_error_limited'
                        } elseif ($residualPruned) {
                            $adaptiveClass = 'span_gate_fallback'
                        } else {
                            $adaptiveClass = 'adaptive_error_limited'
                        }
                    } elseif ($timeRegressed) {
                        $largeProteinClass = 'span_tile_regressed'
                    } elseif ($workRegressed) {
                        $largeProteinClass = 'span_tile_work_regressed'
                    } elseif ($gmresRegressed -and $row.preconditioner_requested -eq 'block') {
                        $largeProteinClass = 'block_precond_regressed'
                    } elseif ($residualPruned) {
                        $largeProteinClass = 'large_protein_candidate'
                    } else {
                        $largeProteinClass = 'span_tile_not_enough'
                    }
                } elseif ($null -ne $leafResidualPanels -and
                          $null -ne $leafTilePanels -and
                          $leafResidualPanels -le 0.1 * [Math]::Max(1.0, $leafTilePanels)) {
                    $largeProteinClass = 'large_protein_residual_pruned'
                }
            }
            if ($row.fixture -eq 'barnase' -and $row.leaf_panel_kernel -eq 'tile-avx2') {
                $leafTileEval = Convert-ToDoubleOrNull $row.leaf_tile_eval_seconds
                $leafResidualEval = Convert-ToDoubleOrNull $row.leaf_residual_eval_seconds
                if ($null -ne $seconds -and $null -ne $ram -and
                    $seconds -lt 124.7 -and $ram -lt 67108864.0) {
                    $barnaseClass = 'barnase_tile_avx2_candidate'
                } elseif ($null -ne $leafResidualPanels -and $null -ne $leafTilePanels -and
                          $leafResidualPanels -gt $leafTilePanels) {
                    $barnaseClass = 'barnase_residual_panel_limited'
                } elseif ($null -ne $leafTileEval -and $null -ne $leafResidualEval -and
                          $leafTileEval -gt $leafResidualEval) {
                    $barnaseClass = 'barnase_tile_eval_limited'
                } else {
                    $barnaseClass = 'barnase_tile_avx2_slow'
                }
            }
            if ($row.leaf_panel_kernel -like 'tile*') {
                $matchingScalar = @($group.Group | Where-Object {
                    $_.mode -eq 'tabi_treecode' -and
                    $_.surface_points -eq $row.surface_points -and
                    $_.tree_leaf_size -eq $row.tree_leaf_size -and
                    $_.opening_theta -eq $row.opening_theta -and
                    $_.leaf_near_policy -eq $row.leaf_near_policy -and
                    $_.cluster_policy -eq $row.cluster_policy -and
                    $_.leaf_panel_kernel -eq 'scalar' -and
                    $_.status -eq 'ok'
                } | Select-Object -First 1)
                if ($matchingScalar.Count -gt 0) {
                    $scalarSeconds = Convert-ToDoubleOrNull $matchingScalar[0].total_seconds
                    if ($null -ne $seconds -and $null -ne $scalarSeconds -and
                        $seconds -lt 0.95 * $scalarSeconds) {
                        $row.tile_fast = 'true'
                    } elseif ($null -ne $leafPanelPairsPerSecond -and
                              $leafPanelPairsPerSecond -gt 100000000.0) {
                        $row.tile_memory_limited = 'true'
                    } else {
                        $row.tile_math_limited = 'true'
                    }
                }
            }
            if ($row.leaf_panel_kernel -eq 'tile-avx2') {
                $matchingTileScalar = @($group.Group | Where-Object {
                    $_.mode -eq 'tabi_treecode' -and
                    $_.surface_points -eq $row.surface_points -and
                    $_.tree_leaf_size -eq $row.tree_leaf_size -and
                    $_.opening_theta -eq $row.opening_theta -and
                    $_.leaf_near_policy -eq $row.leaf_near_policy -and
                    $_.cluster_policy -eq $row.cluster_policy -and
                    $_.leaf_panel_kernel -eq 'tile-scalar' -and
                    $_.status -eq 'ok'
                } | Select-Object -First 1)
                if ($matchingTileScalar.Count -gt 0) {
                    $tileScalarSeconds = Convert-ToDoubleOrNull $matchingTileScalar[0].total_seconds
                    $tileScalarRate = Convert-ToDoubleOrNull $matchingTileScalar[0].leaf_panel_pairs_per_second
                    if ($null -ne $seconds -and $null -ne $tileScalarSeconds -and
                        $seconds -lt 0.95 * $tileScalarSeconds) {
                        $row.tile_avx2_fast = 'true'
                    } elseif ($null -ne $leafPanelPairsPerSecond -and
                              $null -ne $tileScalarRate -and
                              $leafPanelPairsPerSecond -ge 1.5 * $tileScalarRate) {
                        $row.tile_avx2_fast = 'true'
                    } elseif ($null -ne $leafPanelPairsPerSecond -and
                              $leafPanelPairsPerSecond -gt 100000000.0) {
                        $row.tile_avx2_memory_limited = 'true'
                    } else {
                        $row.tile_avx2_math_limited = 'true'
                    }
                }
            }
            if ($row.leaf_near_policy -eq 'smooth-node' -and
                $null -ne $leafPruned -and $leafPruned -gt 0.0) {
                $matchingExact = @($group.Group | Where-Object {
                    $_.mode -eq 'tabi_treecode' -and
                    $_.surface_points -eq $row.surface_points -and
                    $_.tree_leaf_size -eq $row.tree_leaf_size -and
                    $_.opening_theta -eq $row.opening_theta -and
                    $_.leaf_near_policy -eq 'exact' -and
                    $_.cluster_policy -eq $row.cluster_policy -and
                    $_.status -eq 'ok'
                } | Select-Object -First 1)
                if ($matchingExact.Count -gt 0) {
                    $exactSeconds = Convert-ToDoubleOrNull $matchingExact[0].total_seconds
                    if ($null -ne $seconds -and $null -ne $exactSeconds -and
                        $seconds -gt $exactSeconds) {
                        $row.smooth_node_eval_limited = 'true'
                    }
                }
            }
            if (($null -ne $leafNear -and $leafNear -le 420000.0) -or
                ($null -ne $leafPruned -and $leafPruned -gt 0.0 -and
                 $null -ne $leafNear -and $leafNear -gt 0.0 -and
                 $leafPruned -ge 4.0 * $leafNear)) {
                $row.performance_classification = 'interaction_pruned'
            } elseif ($null -ne $largeProteinClass) {
                $row.performance_classification = $largeProteinClass
            } elseif ($null -ne $adaptiveClass) {
                $row.performance_classification = $adaptiveClass
            } elseif ($null -ne $barnaseClass) {
                $row.performance_classification = $barnaseClass
            } elseif ($row.tile_avx2_fast -eq 'true') {
                $row.performance_classification = 'tile_avx2_fast'
            } elseif ($row.tile_avx2_math_limited -eq 'true') {
                $row.performance_classification = 'tile_avx2_math_limited'
            } elseif ($row.tile_avx2_memory_limited -eq 'true') {
                $row.performance_classification = 'tile_avx2_memory_limited'
            } elseif ($row.tile_fast -eq 'true') {
                $row.performance_classification = 'tile_fast'
            } elseif ($row.tile_math_limited -eq 'true') {
                $row.performance_classification = 'tile_math_limited'
            } elseif ($row.tile_memory_limited -eq 'true') {
                $row.performance_classification = 'tile_memory_limited'
            } elseif ($null -ne $near -and $near -gt 1000000.0 -and
                ($null -eq $far -or $far -lt 0.1 * $near)) {
                $row.performance_classification =
                    if ($row.smooth_node_eval_limited -eq 'true') {
                        'smooth_node_eval_limited'
                    } elseif ($null -ne $leafNear -and $leafNear -gt 1000000.0 -and
                        $null -ne $leafPanelPairsPerSecond -and
                        $leafPanelPairsPerSecond -gt 10000000.0) {
                        'leaf_kernel_memory_limited'
                    } else {
                        'operator_limited'
                    }
            } elseif ($null -ne $ratio -and $ratio -le 1.0) {
                $row.performance_classification = 'treecode_candidate'
            } else {
                $row.performance_classification = 'treecode_slow'
            }
        }
    }
}

$outCsv = Join-Path $WorkRoot 'pbe_tabi_perf_matrix.csv'
$rows.ToArray() | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE TABI perf matrix COMPLETE: $outCsv"
