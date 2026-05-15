param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_barnase_truth_matrix'),
    [int]$CaseTimeoutSeconds = 900,
    [string]$NearIntegration = 'centroid',
    [string]$QuadBuffer = 'auto',
    [string]$TrueResidualProbe = 'interval',
    [string]$TrueResidualProbeScope = 'full',
    [int]$TrueResidualProbeInterval = 4,
    [int]$TrueResidualSampleCount = 1024,
    [string[]]$DirectPreconditioners = @('none','diagonal','block'),
    [int[]]$DirectRestarts = @(32,64,100),
    [switch]$IncludeBufferAwareBlock,
    [string[]]$TreeInteractionModes = @('cached','onthefly'),
    [int[]]$TreeLeafSizes = @(32,64),
    [double[]]$TreeOpeningThetas = @(0.5),
    [string[]]$NodeBoundsModes = @('moment'),
    [string[]]$TreeBuilders = @('grouped-morton'),
    [string[]]$MacSafety = @('off'),
    [double]$MacGapFactor = 1.5,
    [string]$MomentAudit = 'off',
    [switch]$ProfileMemory,
    [switch]$Quick,
    [switch]$PrecondOnly,
    [switch]$TreeProbeOnly
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$MacSafety = @($MacSafety | ForEach-Object {
    ([string]$_ -split ',') | Where-Object { -not [string]::IsNullOrWhiteSpace($_) }
})

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'
$PyDelPhiRoot = Join-Path $RepoRoot '.codex_tmp\pyDelPhi'

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI Barnase truth matrix root: $full"
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
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

function Read-DelPhiValueTable([string]$Path) {
    $values = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        $trim = $line.Trim()
        if (-not $trim -or $trim.StartsWith('!') -or $trim.StartsWith('#')) { continue }
        $parts = $trim -split '\s+'
        if ($parts.Count -lt 2) { continue }
        $parsed = 0.0
        if ([double]::TryParse($parts[-1], [System.Globalization.NumberStyles]::Float,
                $InvariantCulture, [ref]$parsed)) {
            $values[$parts[0].ToUpperInvariant()] = $parsed
        }
    }
    return $values
}

function Convert-DelPhiPdbToPqr([string]$PdbPath, [string]$CrgPath,
                                [string]$SizPath, [string]$OutPath) {
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

function Resolve-BarnasePqr([string]$RunDir) {
    $localPqr = Join-Path $InputRoot 'barnase.pqr'
    if (Test-Path -LiteralPath $localPqr) { return $localPqr }
    $fixtureDir = Join-Path $PyDelPhiRoot 'pydelphi\data\test_examples\barnase'
    $directPqr = Join-Path $fixtureDir 'barnase.pqr'
    if (Test-Path -LiteralPath $directPqr) { return $directPqr }
    $pdbPath = Join-Path $fixtureDir 'barnase.pdb'
    $crgPath = Join-Path $fixtureDir 'amber.crg'
    $sizPath = Join-Path $fixtureDir 'amber.siz'
    if ((Test-Path -LiteralPath $pdbPath) -and
        (Test-Path -LiteralPath $crgPath) -and
        (Test-Path -LiteralPath $sizPath)) {
        $generated = Join-Path $RunDir 'barnase.pqr'
        Convert-DelPhiPdbToPqr -PdbPath $pdbPath -CrgPath $crgPath `
            -SizPath $sizPath -OutPath $generated
        return $generated
    }
    return $null
}

function Convert-ToDoubleOrNull($Value) {
    if ($null -eq $Value) { return $null }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $null }
    $parsed = 0.0
    if ([double]::TryParse($text, [System.Globalization.NumberStyles]::Float,
            $InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $null
}

function Get-ValueOrNA($Object, [string]$Name) {
    if ($null -ne $Object -and $Object.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Object.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return 'NA'
}

function Format-Number($Value) {
    if ($null -eq $Value) { return 'NA' }
    if ($Value -is [string]) { return $Value }
    return ([double]$Value).ToString('G17', $InvariantCulture)
}

function Convert-ToProcessArgument([string]$Value) {
    if ($null -eq $Value) { return '""' }
    return '"' + ($Value -replace '"', '\"') + '"'
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
    return [PSCustomObject]@{
        Completed = $completed
        ExitCode = if ($completed) { $process.ExitCode } else { $null }
        ElapsedSeconds = $watch.Elapsed.TotalSeconds
    }
}

function New-TruthArgFile($Case, [string]$RunDir, [string]$PqrPath) {
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
        "-pbe=$($Case.backend)",
        "-pbe-test-pqr=$(Convert-ToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-probe=1.4',
        '-pbe-scale=2.0',
        '-pbe-margin=0.0',
        '-pbe-max-grid=4000000',
        '-pbe-map=off',
        '-pbe-max-iter=1200',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-timing=detail',
        '-pbe-salt=0.15',
        '-pbe-temp=298.15',
        "-nt=$($Case.threads)"
    )
    if ($Case.backend -eq 'tabi') {
        $argLines += "-pbe-tabi-sum=$($Case.sum_mode)"
        $argLines += "-pbe-tabi-formulation=$($Case.formulation)"
        $argLines += "-pbe-tabi-precond=$($Case.preconditioner)"
        if ($Case.PSObject.Properties.Name -contains 'precond_block_near') {
            $argLines += "-pbe-tabi-precond-block-near=$($Case.precond_block_near)"
        }
        $argLines += '-pbe-tabi-derivatives=analytic'
        $argLines += '-pbe-tabi-energy-projection=boundary'
        $argLines += "-pbe-tabi-near-integration=$NearIntegration"
        $argLines += "-pbe-tabi-quad-buffer=$QuadBuffer"
        $argLines += '-pbe-tabi-surface-points=8'
        $argLines += "-pbe-tabi-gmres-restart=$($Case.restart)"
        $argLines += "-pbe-tabi-true-residual-probe=$TrueResidualProbe"
        $argLines += "-pbe-tabi-true-residual-probe-scope=$TrueResidualProbeScope"
        $argLines += "-pbe-tabi-true-residual-probe-interval=$TrueResidualProbeInterval"
        $argLines += "-pbe-tabi-true-residual-sample-count=$TrueResidualSampleCount"
        $argLines += '-pbe-tabi-progress=auto'
        $argLines += '-pbe-tabi-progress-interval=5'
        if ($Case.sum_mode -eq 'treecode') {
            $argLines += "-pbe-tabi-tree-leaf-size=$($Case.leaf)"
            $argLines += "-pbe-tabi-opening-theta=$($Case.theta.ToString('G17', $InvariantCulture))"
            $argLines += '-pbe-tabi-leaf-near-policy=exact'
            $argLines += '-pbe-tabi-leaf-panel-kernel=tile-avx2'
            $argLines += '-pbe-tabi-adaptive-level=strict'
            $argLines += '-pbe-tabi-cluster-policy=off'
            $argLines += "-pbe-tabi-interaction-mode=$($Case.interaction_mode)"
            if ($Case.PSObject.Properties.Name -contains 'node_bounds') {
                $argLines += "-pbe-tabi-node-bounds=$($Case.node_bounds)"
            }
            if ($Case.PSObject.Properties.Name -contains 'tree_builder') {
                $argLines += "-pbe-tabi-tree-builder=$($Case.tree_builder)"
            }
            if ($Case.PSObject.Properties.Name -contains 'mac_safety') {
                $argLines += "-pbe-tabi-mac-safety=$($Case.mac_safety)"
            }
            if ($Case.PSObject.Properties.Name -contains 'mac_gap_factor') {
                $argLines += "-pbe-tabi-mac-gap-factor=$($Case.mac_gap_factor.ToString('G17', $InvariantCulture))"
            }
            if ($Case.PSObject.Properties.Name -contains 'moment_audit') {
                $argLines += "-pbe-tabi-moment-audit=$($Case.moment_audit)"
            }
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function New-Case([string]$Name, [string]$Backend, [string]$SumMode,
                  [string]$Formulation, [string]$Preconditioner, [int]$Restart,
                  [int]$Leaf, [double]$Theta, [int]$Threads,
                  [string]$InteractionMode = 'NA',
                  [string]$PrecondBlockNear = 'singular',
                  [string]$NodeBounds = 'moment',
                  [string]$TreeBuilder = 'grouped-morton',
                  [string]$MacSafety = 'off',
                  [double]$MacGapFactor = 1.5,
                  [string]$MomentAudit = 'off') {
    return [PSCustomObject]@{
        name = $Name
        backend = $Backend
        sum_mode = $SumMode
        formulation = $Formulation
        preconditioner = $Preconditioner
        precond_block_near = $PrecondBlockNear
        restart = $Restart
        leaf = $Leaf
        theta = $Theta
        threads = $Threads
        interaction_mode = $InteractionMode
        node_bounds = $NodeBounds
        tree_builder = $TreeBuilder
        mac_safety = $MacSafety
        mac_gap_factor = $MacGapFactor
        moment_audit = $MomentAudit
    }
}

function Get-LastProgress([string]$RunDir) {
    $path = Join-Path $RunDir 'pbe_tabi_longrun_progress.csv'
    if (Test-Path -LiteralPath $path) {
        $rows = @(Import-Csv -LiteralPath $path)
        if ($rows.Count -gt 0) { return $rows[-1] }
    }
    return $null
}

function Invoke-TruthCase($Case) {
    $safeName = ($Case.name -replace '[^A-Za-z0-9_.-]', '_')
    $runDir = Join-Path $WorkRoot $safeName
    New-CleanDir $runDir
    $pqrPath = Resolve-BarnasePqr -RunDir $runDir
    if ($null -eq $pqrPath -or -not (Test-Path -LiteralPath $pqrPath)) {
        throw "Missing Barnase PQR fixture/cache."
    }
    $argPath = New-TruthArgFile -Case $Case -RunDir $runDir -PqrPath $pqrPath
    $runResult = Invoke-TitanWithTimeout -ArgPath $argPath -RunDir $runDir `
        -TimeoutSeconds $CaseTimeoutSeconds

    $timingPath = Join-Path $runDir 'pbe_timing.csv'
    $energyPath = Join-Path $runDir 'pbe_energy.csv'
    $memoryPath = Join-Path $runDir 'pbe_memory_profile.csv'
    $timing = if (Test-Path -LiteralPath $timingPath) { @(Import-Csv -LiteralPath $timingPath)[0] } else { $null }
    $energy = if (Test-Path -LiteralPath $energyPath) { @(Import-Csv -LiteralPath $energyPath)[0] } else { $null }
    $memory = if (Test-Path -LiteralPath $memoryPath) { @(Import-Csv -LiteralPath $memoryPath)[0] } else { $null }
    if ($ProfileMemory -and $runResult.Completed -and $Case.backend -eq 'tabi' -and $null -eq $memory) {
        throw "Missing memory profile for $($Case.name)"
    }
    $progress = Get-LastProgress -RunDir $runDir
    $status = if (-not $runResult.Completed) {
        'timeout'
    } elseif ($null -ne $energy) {
        $energy.status
    } elseif ($runResult.ExitCode -ne 0) {
        'failed_exit'
    } else {
        'unknown'
    }
    $matvec = Convert-ToDoubleOrNull (Get-ValueOrNA $progress 'matvec_calls')
    if ($null -eq $matvec) { $matvec = Convert-ToDoubleOrNull (Get-ValueOrNA $timing 'tabi_matvec_calls') }
    $leafPanel = Convert-ToDoubleOrNull (Get-ValueOrNA $progress 'leaf_panel_interactions')
    if ($null -eq $leafPanel) { $leafPanel = Convert-ToDoubleOrNull (Get-ValueOrNA $timing 'tabi_leaf_panel_interactions') }
    $leafPerMatvec = if ($null -ne $leafPanel -and $null -ne $matvec -and $matvec -gt 0.0) {
        $leafPanel / $matvec
    } else { $null }
    return [PSCustomObject]@{
        case_id = $Case.name
        backend = $Case.backend
        sum_mode = if ($Case.backend -eq 'tabi') { $Case.sum_mode } else { 'grid' }
        theta = if ($Case.sum_mode -eq 'treecode') { Format-Number $Case.theta } else { 'NA' }
        leaf = if ($Case.sum_mode -eq 'treecode') { $Case.leaf } else { 'NA' }
        preconditioner = if ($Case.backend -eq 'tabi') { $Case.preconditioner } else { 'NA' }
        preconditioner_block_near_requested = if ($Case.backend -eq 'tabi') { $Case.precond_block_near } else { 'NA' }
        restart = if ($Case.backend -eq 'tabi') { $Case.restart } else { 'NA' }
        near_integration = if ($Case.backend -eq 'tabi') { $NearIntegration } else { 'NA' }
        quad_buffer = if ($Case.backend -eq 'tabi') { $QuadBuffer } else { 'NA' }
        true_residual_probe_mode_requested = if ($Case.backend -eq 'tabi') { $TrueResidualProbe } else { 'NA' }
        true_residual_probe_scope_requested = if ($Case.backend -eq 'tabi') { $TrueResidualProbeScope } else { 'NA' }
        true_residual_probe_interval_requested = if ($Case.backend -eq 'tabi') { $TrueResidualProbeInterval } else { 'NA' }
        true_residual_sample_count_requested = if ($Case.backend -eq 'tabi') { $TrueResidualSampleCount } else { 'NA' }
        requested_interaction_mode = if ($Case.sum_mode -eq 'treecode') { $Case.interaction_mode } else { 'NA' }
        requested_node_bounds = if ($Case.sum_mode -eq 'treecode') { $Case.node_bounds } else { 'NA' }
        requested_tree_builder = if ($Case.sum_mode -eq 'treecode') { $Case.tree_builder } else { 'NA' }
        threads = $Case.threads
        status = $status
        exit_code = if ($runResult.Completed) { $runResult.ExitCode } else { 'timeout' }
        total_seconds = Format-Number $runResult.ElapsedSeconds
        energy_kJmol = Get-ValueOrNA $energy 'E_rxn_kJmol'
        panels = Get-ValueOrNA $progress 'panels'
        unknown_count = Get-ValueOrNA $progress 'unknown_count'
        gmres_iter = Get-ValueOrNA $progress 'gmres_iter'
        matvec_calls = Format-Number $matvec
        best_residual = Get-ValueOrNA $progress 'best_residual'
        best_residual_iter = Get-ValueOrNA $progress 'best_residual_iter'
        final_residual = Get-ValueOrNA $progress 'final_residual'
        timing_best_residual = Get-ValueOrNA $timing 'tabi_best_residual'
        timing_final_residual = Get-ValueOrNA $timing 'tabi_final_residual'
        best_preconditioned_residual = Get-ValueOrNA $progress 'best_preconditioned_residual'
        final_preconditioned_residual = Get-ValueOrNA $progress 'final_preconditioned_residual'
        true_residual_abs = Get-ValueOrNA $progress 'true_residual_abs'
        true_residual_rel = Get-ValueOrNA $progress 'true_residual_rel'
        timing_best_preconditioned_residual = Get-ValueOrNA $timing 'tabi_best_preconditioned_residual'
        timing_final_preconditioned_residual = Get-ValueOrNA $timing 'tabi_final_preconditioned_residual'
        timing_true_residual_abs = Get-ValueOrNA $timing 'tabi_true_residual_abs'
        timing_true_residual_rel = Get-ValueOrNA $timing 'tabi_true_residual_rel'
        true_residual_probe_mode = Get-ValueOrNA $progress 'true_residual_probe_mode'
        true_residual_probe_scope = Get-ValueOrNA $progress 'true_residual_probe_scope'
        true_residual_probe_interval = Get-ValueOrNA $progress 'true_residual_probe_interval'
        true_residual_sample_count = Get-ValueOrNA $progress 'true_residual_sample_count'
        true_residual_probe_calls = Get-ValueOrNA $progress 'true_residual_probe_calls'
        true_residual_probe_seconds = Get-ValueOrNA $progress 'true_residual_probe_seconds'
        sampled_true_residual_rel = Get-ValueOrNA $progress 'sampled_true_residual_rel'
        sampled_true_residual_rows = Get-ValueOrNA $progress 'sampled_true_residual_rows'
        sampled_probe_calls = Get-ValueOrNA $progress 'sampled_probe_calls'
        sampled_probe_seconds = Get-ValueOrNA $progress 'sampled_probe_seconds'
        full_true_residual_confirm_calls = Get-ValueOrNA $progress 'full_true_residual_confirm_calls'
        full_probe_seconds = Get-ValueOrNA $progress 'full_probe_seconds'
        timing_true_residual_probe_mode = Get-ValueOrNA $timing 'tabi_true_residual_probe_mode'
        timing_true_residual_probe_scope = Get-ValueOrNA $timing 'tabi_true_residual_probe_scope'
        timing_true_residual_probe_interval = Get-ValueOrNA $timing 'tabi_true_residual_probe_interval'
        timing_true_residual_sample_count = Get-ValueOrNA $timing 'tabi_true_residual_sample_count'
        timing_true_residual_probe_calls = Get-ValueOrNA $timing 'tabi_true_residual_probe_calls'
        timing_true_residual_probe_seconds = Get-ValueOrNA $timing 'tabi_true_residual_probe_seconds'
        timing_sampled_true_residual_rel = Get-ValueOrNA $timing 'tabi_sampled_true_residual_rel'
        timing_sampled_true_residual_rows = Get-ValueOrNA $timing 'tabi_sampled_true_residual_rows'
        timing_sampled_probe_calls = Get-ValueOrNA $timing 'tabi_sampled_probe_calls'
        timing_sampled_probe_seconds = Get-ValueOrNA $timing 'tabi_sampled_probe_seconds'
        timing_full_true_residual_confirm_calls = Get-ValueOrNA $timing 'tabi_full_true_residual_confirm_calls'
        timing_full_probe_seconds = Get-ValueOrNA $timing 'tabi_full_probe_seconds'
        timing_convergence_reason = Get-ValueOrNA $timing 'tabi_convergence_reason'
        interaction_mode = Get-ValueOrNA $timing 'tabi_interaction_mode'
        node_bounds_mode = Get-ValueOrNA $timing 'tabi_node_bounds_mode'
        tree_builder = Get-ValueOrNA $timing 'tabi_tree_builder'
        mac_safety = Get-ValueOrNA $timing 'tabi_mac_safety'
        mac_gap_factor = Get-ValueOrNA $timing 'tabi_mac_gap_factor'
        tree_depth_max = Get-ValueOrNA $timing 'tabi_tree_depth_max'
        mac_reject_radius_limited = Get-ValueOrNA $timing 'tabi_mac_reject_radius_limited'
        mac_accept_rate = Get-ValueOrNA $timing 'tabi_mac_accept_rate'
        node_radius_max = Get-ValueOrNA $timing 'tabi_node_radius_max'
        node_radius_mean = Get-ValueOrNA $timing 'tabi_node_radius_mean'
        leaf_panel_work_per_matvec = Get-ValueOrNA $timing 'tabi_leaf_panel_work_per_matvec'
        far_node_work_per_matvec = Get-ValueOrNA $timing 'tabi_far_node_work_per_matvec'
        far_geometric_unsafe_count = Get-ValueOrNA $timing 'tabi_far_geometric_unsafe_count'
        far_geometric_unsafe_panels = Get-ValueOrNA $timing 'tabi_far_geometric_unsafe_panels'
        far_geometric_min_gap = Get-ValueOrNA $timing 'tabi_far_geometric_min_gap'
        far_geometric_gap_ratio_min = Get-ValueOrNA $timing 'tabi_far_geometric_gap_ratio_min'
        far_geometric_gap_ratio_mean = Get-ValueOrNA $timing 'tabi_far_geometric_gap_ratio_mean'
        far_geometric_audit_samples = Get-ValueOrNA $timing 'tabi_far_geometric_audit_samples'
        mac_physical_gap_rejects = Get-ValueOrNA $timing 'tabi_mac_physical_gap_rejects'
        moment_audit = Get-ValueOrNA $timing 'tabi_moment_audit'
        moment_audit_nodes = Get-ValueOrNA $timing 'tabi_moment_audit_nodes'
        moment_translation_rel_error_max = Get-ValueOrNA $timing 'tabi_moment_translation_rel_error_max'
        moment_translation_rel_error_mean = Get-ValueOrNA $timing 'tabi_moment_translation_rel_error_mean'
        leaf_panel_interactions = Format-Number $leafPanel
        leaf_panel_per_matvec = Format-Number $leafPerMatvec
        rss_bytes = Get-ValueOrNA $progress 'rss_bytes'
        peak_working_set_bytes = Get-ValueOrNA $memory 'peak_working_set_bytes'
        progress_phase = Get-ValueOrNA $progress 'phase'
        preconditioner_effective = if ((Get-ValueOrNA $timing 'tabi_preconditioner_effective') -ne 'NA') {
            Get-ValueOrNA $timing 'tabi_preconditioner_effective'
        } else {
            Get-ValueOrNA $progress 'preconditioner_effective'
        }
        preconditioner_block_near = Get-ValueOrNA $timing 'tabi_preconditioner_block_near'
        decision_classification = 'pending'
        run_dir = $runDir
    }
}

function Convert-RowResidual($Row) {
    $value = Convert-ToDoubleOrNull $Row.timing_final_residual
    if ($null -eq $value) { $value = Convert-ToDoubleOrNull $Row.final_residual }
    return $value
}

function Convert-RowBestResidual($Row) {
    $value = Convert-ToDoubleOrNull $Row.timing_best_residual
    if ($null -eq $value) { $value = Convert-ToDoubleOrNull $Row.best_residual }
    return $value
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
New-CleanDir $WorkRoot

$cases = New-Object System.Collections.Generic.List[object]
$cases.Add((New-Case -Name 'grid_reference' -Backend 'grid' -SumMode 'grid' `
    -Formulation 'grid' -Preconditioner 'NA' -Restart 0 -Leaf 0 -Theta 0.0 -Threads 8))
foreach ($precond in $DirectPreconditioners) {
    foreach ($restart in $DirectRestarts) {
        $cases.Add((New-Case -Name ("direct_{0}_r{1}" -f $precond,$restart) `
            -Backend 'tabi' -SumMode 'direct' -Formulation 'jcp2013' `
            -Preconditioner $precond -Restart $restart -Leaf 0 -Theta 0.0 -Threads 8))
    }
}
if ($IncludeBufferAwareBlock) {
    foreach ($restart in $DirectRestarts) {
        $cases.Add((New-Case -Name ("direct_block_buffer_aware_r{0}" -f $restart) `
            -Backend 'tabi' -SumMode 'direct' -Formulation 'jcp2013' `
            -Preconditioner 'block' -PrecondBlockNear 'buffer-aware' `
            -Restart $restart -Leaf 0 -Theta 0.0 -Threads 8))
    }
}
foreach ($precond in @('diagonal','block')) {
    foreach ($leaf in @(32,64)) {
        foreach ($restart in @(64,100)) {
            foreach ($interactionMode in $TreeInteractionModes) {
                $safeInteraction = ($interactionMode -replace '[^A-Za-z0-9_-]', '_')
                $cases.Add((New-Case -Name ("tree_theta0p5_leaf{0}_{1}_r{2}_im{3}" -f $leaf,$precond,$restart,$safeInteraction) `
                    -Backend 'tabi' -SumMode 'treecode' -Formulation 'jcp2013-treecode' `
                    -Preconditioner $precond -Restart $restart -Leaf $leaf -Theta 0.5 -Threads 8 `
                    -InteractionMode $interactionMode))
            }
        }
    }
}
foreach ($restart in @(64,100)) {
    foreach ($interactionMode in $TreeInteractionModes) {
        $safeInteraction = ($interactionMode -replace '[^A-Za-z0-9_-]', '_')
        $cases.Add((New-Case -Name ("tree_theta0p8_leaf64_block_r{0}_im{1}" -f $restart,$safeInteraction) `
            -Backend 'tabi' -SumMode 'treecode' -Formulation 'jcp2013-treecode' `
            -Preconditioner 'block' -Restart $restart -Leaf 64 -Theta 0.8 -Threads 8 `
            -InteractionMode $interactionMode))
    }
}
if ($TreeProbeOnly) {
    $probeCases = New-Object System.Collections.Generic.List[object]
    $probeCases.Add((New-Case -Name 'grid_reference' -Backend 'grid' -SumMode 'grid' `
        -Formulation 'grid' -Preconditioner 'NA' -Restart 0 -Leaf 0 -Theta 0.0 -Threads 8))
    foreach ($leaf in $TreeLeafSizes) {
        foreach ($theta in $TreeOpeningThetas) {
            $safeTheta = ($theta.ToString('G17', $InvariantCulture) -replace '\.','p')
            foreach ($bounds in $NodeBoundsModes) {
                $safeBounds = ($bounds -replace '[^A-Za-z0-9_-]', '_')
                foreach ($builder in $TreeBuilders) {
                    $safeBuilder = ($builder -replace '[^A-Za-z0-9_-]', '_')
                    foreach ($mac in $MacSafety) {
                        $safeMac = ($mac -replace '[^A-Za-z0-9_-]', '_')
                        $probeCases.Add((New-Case -Name ("tree_mac_leaf{0}_theta{1}_{2}_{3}_{4}" -f $leaf,$safeTheta,$safeBounds,$safeBuilder,$safeMac) `
                            -Backend 'tabi' -SumMode 'treecode' -Formulation 'jcp2013-treecode' `
                            -Preconditioner 'diagonal' -Restart 64 -Leaf $leaf -Theta $theta -Threads 8 `
                            -InteractionMode 'onthefly' -NodeBounds $bounds -TreeBuilder $builder `
                            -MacSafety $mac -MacGapFactor $MacGapFactor -MomentAudit $MomentAudit))
                    }
                }
            }
        }
    }
    $cases = $probeCases
} elseif ($PrecondOnly) {
    $precondCases = New-Object System.Collections.Generic.List[object]
    foreach ($name in @('grid_reference',
                        'direct_none_r64','direct_none_r100',
                        'direct_diagonal_r64','direct_diagonal_r100',
                        'direct_block_r64','direct_block_r100',
                        'direct_block_buffer_aware_r64',
                        'direct_block_buffer_aware_r100')) {
        $match = @($cases | Where-Object { $_.name -eq $name } | Select-Object -First 1)
        if ($match.Count -gt 0) { $precondCases.Add($match[0]) }
    }
    $cases = $precondCases
} elseif ($Quick) {
    $quickCases = New-Object System.Collections.Generic.List[object]
    foreach ($name in @('grid_reference','direct_none_r64','direct_diagonal_r64','direct_block_r64','direct_block_buffer_aware_r64')) {
        $match = @($cases | Where-Object { $_.name -eq $name } | Select-Object -First 1)
        if ($match.Count -gt 0) { $quickCases.Add($match[0]) }
    }
    foreach ($match in @($cases | Where-Object {
            $_.name -like 'tree_theta0p5_leaf32_diagonal_r64_im*'
        })) {
        $quickCases.Add($match)
    }
    $cases = $quickCases
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($case in $cases) {
    $rows.Add((Invoke-TruthCase -Case $case))
}

$directRows = @($rows | Where-Object { $_.backend -eq 'tabi' -and $_.sum_mode -eq 'direct' })
$treeRows = @($rows | Where-Object { $_.backend -eq 'tabi' -and $_.sum_mode -eq 'treecode' })
$directBest = @($directRows | ForEach-Object {
    [PSCustomObject]@{ Row = $_; Final = Convert-RowResidual $_; Best = Convert-RowBestResidual $_ }
} | Where-Object { $null -ne $_.Best } | Sort-Object Best | Select-Object -First 1)
$treeSafe = @($treeRows | Where-Object {
    $_.theta -like '0.5*' -and $_.status -eq 'ok' -and
    (Convert-RowResidual $_) -le 1.0e-4
} | Select-Object -First 1)
$probeConvergedRows = @($rows | Where-Object {
    $_.timing_convergence_reason -eq 'converged_true_residual_probe' -or
    ($_.timing_true_residual_probe_calls -ne 'NA' -and
     (Convert-ToDoubleOrNull $_.timing_true_residual_probe_calls) -gt 0.0 -and
     (Convert-RowBestResidual $_) -ne $null -and
     (Convert-RowBestResidual $_) -le 1.0e-4)
})
$probeTooExpensiveRows = @($rows | Where-Object {
    $probeSeconds = Convert-ToDoubleOrNull $_.timing_true_residual_probe_seconds
    $totalSeconds = Convert-ToDoubleOrNull $_.total_seconds
    $probeSeconds -ne $null -and $totalSeconds -ne $null -and
        $totalSeconds -gt 0.0 -and ($probeSeconds / $totalSeconds) -gt 0.35
})
$sampledOnlyRows = @($directRows | Where-Object {
    $sampleCalls = Convert-ToDoubleOrNull $_.timing_sampled_probe_calls
    $fullCalls = Convert-ToDoubleOrNull $_.timing_full_true_residual_confirm_calls
    $sampleCalls -ne $null -and $sampleCalls -gt 0.0 -and
        ($fullCalls -eq $null -or $fullCalls -le 0.0)
})
$decision = 'inconclusive'
if ($probeConvergedRows.Count -gt 0) {
    $decision = 'true_probe_converged'
} elseif ($sampledOnlyRows.Count -gt 0 -and $directBest.Count -gt 0 -and
          $directBest[0].Best -ge 0.999) {
    $decision = 'sampled_probe_not_predictive'
} elseif ($directBest.Count -gt 0 -and $directBest[0].Best -gt 1.0e-3) {
    $decision = 'true_residual_failed'
} elseif ($directBest.Count -gt 0 -and $directBest[0].Best -le 1.0e-4 -and $treeSafe.Count -eq 0) {
    $decision = 'treecode_noise_limited'
} elseif ($treeSafe.Count -gt 0) {
    $decision = 'safe_treecode_candidate'
}
$restart32 = @($directRows | Where-Object { $_.restart -eq 32 } | ForEach-Object {
    [PSCustomObject]@{ Best = Convert-RowBestResidual $_ }
} | Where-Object { $null -ne $_.Best } | Sort-Object Best | Select-Object -First 1)
$restartLong = @($directRows | Where-Object { $_.restart -ne 32 } | ForEach-Object {
    [PSCustomObject]@{ Best = Convert-RowBestResidual $_ }
} | Where-Object { $null -ne $_.Best } | Sort-Object Best | Select-Object -First 1)
if ($restart32.Count -gt 0 -and $restartLong.Count -gt 0 -and
    $restartLong[0].Best -lt 0.5 * $restart32[0].Best) {
    $decision = 'restart_limited'
}
$phantomRows = @($rows | Where-Object {
    $pre = Convert-ToDoubleOrNull $_.timing_best_preconditioned_residual
    $trueResidualValue = Convert-RowBestResidual $_
    $pre -ne $null -and $trueResidualValue -ne $null -and
        $pre -le 1.0e-4 -and $trueResidualValue -gt 1.0e-4
})
if ($phantomRows.Count -gt 0) {
    $decision = 'phantom_preconditioned'
}
$directNoneBest = @($directRows | Where-Object { $_.preconditioner -eq 'none' } | ForEach-Object {
    [PSCustomObject]@{ Best = Convert-RowBestResidual $_ }
} | Where-Object { $null -ne $_.Best } | Sort-Object Best | Select-Object -First 1)
$directDiagBest = @($directRows | Where-Object { $_.preconditioner -eq 'diagonal' } | ForEach-Object {
    [PSCustomObject]@{ Best = Convert-RowBestResidual $_ }
} | Where-Object { $null -ne $_.Best } | Sort-Object Best | Select-Object -First 1)
$directBlockBest = @($directRows | Where-Object {
    $_.preconditioner -eq 'block' -and $_.preconditioner_block_near_requested -eq 'singular'
} | ForEach-Object {
    [PSCustomObject]@{ Best = Convert-RowBestResidual $_ }
} | Where-Object { $null -ne $_.Best } | Sort-Object Best | Select-Object -First 1)
$directBufferBlockBest = @($directRows | Where-Object {
    $_.preconditioner -eq 'block' -and $_.preconditioner_block_near_requested -eq 'buffer-aware'
} | ForEach-Object {
    [PSCustomObject]@{ Best = Convert-RowBestResidual $_ }
} | Where-Object { $null -ne $_.Best } | Sort-Object Best | Select-Object -First 1)
if ($directBlockBest.Count -gt 0 -and $directDiagBest.Count -gt 0 -and
    $directBlockBest[0].Best -gt [Math]::Max(1.0e-3, 2.0 * $directDiagBest[0].Best)) {
    $decision = 'block_precond_destructive'
} elseif ($directBufferBlockBest.Count -gt 0 -and $directDiagBest.Count -gt 0 -and
          $directBufferBlockBest[0].Best -lt $directDiagBest[0].Best) {
    $decision = 'block_buffer_aware_candidate'
} elseif ($directDiagBest.Count -gt 0 -and $directDiagBest[0].Best -lt 1.0e-3) {
    $decision = 'diagonal_truth_candidate'
} elseif ($directNoneBest.Count -gt 0 -and $directNoneBest[0].Best -lt 1.0e-3) {
    $decision = 'none_truth_candidate'
}
if ($probeTooExpensiveRows.Count -gt 0 -and $decision -eq 'inconclusive') {
    $decision = 'probe_too_expensive'
}
$bestRows = @($rows | Where-Object {
    $_.timing_convergence_reason -eq 'converged_true_residual_best'
})
if ($bestRows.Count -gt 0 -and $decision -eq 'inconclusive') {
    $decision = 'converged_true_residual_best'
}

foreach ($row in $rows) { $row.decision_classification = $decision }

$outCsv = Join-Path $WorkRoot 'pbe_tabi_barnase_truth_matrix.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE TABI Barnase truth matrix COMPLETE: $outCsv decision=$decision"
