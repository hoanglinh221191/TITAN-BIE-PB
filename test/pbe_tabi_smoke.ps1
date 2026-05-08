param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootFull = [System.IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
    if (-not $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean path outside PBE TABI smoke root: $full"
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
        'CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -3.200   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -1.800   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1      -1.160   0.500  -1.320  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -1.720   1.260  -2.100  1.00  0.00           O',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Read-SummaryMap([string]$Path) {
    $map = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ([string]::IsNullOrWhiteSpace($line) -or $line.StartsWith('#')) {
            continue
        }
        $parts = $line -split '\s+', 2
        if ($parts.Count -eq 2) {
            $map[$parts[0]] = $parts[1]
        }
    }
    return $map
}

function New-ArgFile([string]$RunDir, [string]$Backend, [string]$SumMode, [string]$Formulation = 'jcp2013') {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $pqrPath = Join-Path $InputRoot 'mini-protein.pqr'
    if (-not (Test-Path -LiteralPath $pqrPath)) {
        throw "Missing mini-protein PQR fixture: $pqrPath"
    }
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        "-pbe=$Backend",
        "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
        '-pbe-surface=vdw',
        '-pbe-max-iter=800',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-timing=summary',
        "-pbe-tabi-sum=$SumMode",
        "-pbe-tabi-formulation=$Formulation",
        '-nt=1'
    )
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Invoke-TabiCase([string]$Name, [string]$SumMode, [string]$Formulation = 'jcp2013') {
    $runDir = Join-Path $WorkRoot $Name
    New-CleanDir $runDir
    $argPath = New-ArgFile -RunDir $runDir -Backend 'tabi' -SumMode $SumMode -Formulation $Formulation

    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for $Name"
        }
    } finally {
        Pop-Location
    }

    $energyPath = Join-Path $runDir 'pbe_energy.csv'
    $summaryPath = Join-Path $runDir 'pbe_summary.dat'
    $timingPath = Join-Path $runDir 'pbe_timing.csv'
    foreach ($required in @($energyPath, $summaryPath, $timingPath)) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Missing TABI PBE output for ${Name}: $required"
        }
    }

    $rows = @(Import-Csv -LiteralPath $energyPath)
    if ($rows.Count -ne 1) {
        throw "Expected one TABI PBE row for ${Name}, observed $($rows.Count)"
    }
    $row = $rows[0]
    if ($row.pbe_backend -ne 'tabi') {
        throw "Expected pbe_backend=tabi for ${Name}, observed $($row.pbe_backend)"
    }
    if ($row.status -ne 'ok') {
        throw "TABI PBE case did not converge for ${Name}: $($row.status)"
    }
    if (-not $row.solver_mode.StartsWith('tabi_')) {
        throw "TABI PBE row did not report TABI solver_mode for ${Name}: $($row.solver_mode)"
    }
    $energyKJ = Convert-ToDouble $row.E_rxn_kJmol
    if ([double]::IsNaN($energyKJ) -or [double]::IsInfinity($energyKJ)) {
        throw "Non-finite TABI PBE energy for ${Name}: $energyKJ"
    }

    $summary = Read-SummaryMap $summaryPath
    if ($summary['pbe_backend'] -ne 'tabi') {
        throw "TABI summary did not report pbe_backend=tabi for $Name"
    }

    $timingRows = @(Import-Csv -LiteralPath $timingPath)
    if ($timingRows.Count -lt 1 -or $timingRows[0].pbe_backend -ne 'tabi') {
        throw "TABI timing row did not carry pbe_backend=tabi for $Name"
    }
    $timingRow = $timingRows[0]
    foreach ($column in @('tabi_surface_build_seconds','tabi_surface_panels',
            'tabi_surface_exposed_points','tabi_surface_buried_points',
            'tabi_tree_build_seconds','tabi_tree_nodes','tabi_tree_leaf_count',
            'tabi_near_interactions','tabi_far_interactions','tabi_matvec_mode',
            'tabi_surface_geometry','tabi_mesh_vertices','tabi_mesh_triangles',
            'tabi_boundary_nodes','tabi_avg_valence','tabi_nonmanifold_edges',
            'tabi_solid_angle_fallbacks','tabi_diag_min','tabi_diag_max',
            'tabi_taylor_order','tabi_interaction_build_seconds',
            'tabi_tree_eval_seconds','tabi_status',
            'tabi_convergence_reason','tabi_reaction_energy_source',
            'tabi_rhs_seconds','tabi_gmres_seconds',
            'tabi_matvec_seconds_total','tabi_matvec_seconds_avg',
            'tabi_matvec_calls','tabi_orthogonalization_seconds',
            'tabi_preconditioner_seconds','tabi_surface_neighbor_seconds',
            'tabi_surface_mesh_seconds','tabi_solid_angle_seconds',
            'tabi_tree_moment_seconds','tabi_tree_traversal_seconds',
            'tabi_near_eval_seconds','tabi_far_eval_seconds',
            'tabi_formulation','tabi_unknown_count','tabi_bie_system',
            'tabi_lhs_diag_phi','tabi_lhs_diag_phin',
            'tabi_energy_equation','tabi_kernel_convention',
            'tabi_source_norm_s1','tabi_source_norm_s2',
            'tabi_phi_norm','tabi_phin_norm')) {
        if (-not ($timingRow.PSObject.Properties.Name -contains $column)) {
            throw "TABI timing row missing column ${column} for $Name"
        }
    }
    if ($timingRow.tabi_status -ne 'ok') {
        throw "TABI timing status was not ok for ${Name}: $($timingRow.tabi_status)"
    }
    if ($timingRow.tabi_reaction_energy_source -ne 'tabi_boundary_integral') {
        throw "TABI unexpectedly used fallback for ${Name}: $($timingRow.tabi_reaction_energy_source)"
    }
    if ([int64]$timingRow.tabi_matvec_calls -le 0) {
        throw "TABI did not report matvec calls for ${Name}"
    }
    $panelCount = [int64]$timingRow.tabi_surface_panels
    if ($panelCount -le 476) {
        throw "TABI did not build the expected triangle boundary panels for ${Name}: panels=$panelCount"
    }
    if ($timingRow.tabi_surface_geometry -ne 'triangle') {
        throw "TABI did not use triangle surface geometry for ${Name}: $($timingRow.tabi_surface_geometry)"
    }
    if ([int64]$timingRow.tabi_mesh_vertices -le 0 -or [int64]$timingRow.tabi_mesh_triangles -ne $panelCount) {
        throw "TABI triangle diagnostics invalid for ${Name}: vertices=$($timingRow.tabi_mesh_vertices) triangles=$($timingRow.tabi_mesh_triangles) panels=$panelCount"
    }
    if ([int64]$timingRow.tabi_boundary_nodes -ne $panelCount) {
        throw "TABI boundary node count did not match panels for ${Name}: nodes=$($timingRow.tabi_boundary_nodes) panels=$panelCount"
    }
    $diagMin = Convert-ToDouble $timingRow.tabi_diag_min
    $diagMax = Convert-ToDouble $timingRow.tabi_diag_max
    if ($diagMin -le 0.0 -or $diagMax -lt $diagMin) {
        throw "TABI solid-angle diagonal range invalid for ${Name}: min=$diagMin max=$diagMax"
    }
    $expectedTaylor = if ($Formulation -eq 'jcp2013') { 0 } else { 3 }
    if ([int]$timingRow.tabi_taylor_order -ne $expectedTaylor) {
        throw "TABI default Taylor order drifted for ${Name}: $($timingRow.tabi_taylor_order)"
    }
    if ([int64]$row.grid_points -ne $panelCount) {
        throw "TABI energy grid_points did not match surface panel count for ${Name}: grid=$($row.grid_points) panels=$panelCount"
    }
    if ($SumMode -eq 'treecode') {
        if ($timingRow.tabi_matvec_mode -ne 'treecode') {
            throw "TABI treecode case reported matvec_mode=$($timingRow.tabi_matvec_mode)"
        }
        if ([int64]$timingRow.tabi_tree_nodes -le 0 -or [int64]$timingRow.tabi_tree_leaf_count -le 0) {
            throw "TABI treecode case did not report tree nodes/leaves."
        }
    } elseif ($SumMode -eq 'direct' -and $timingRow.tabi_matvec_mode -ne 'direct') {
        throw "TABI direct case reported matvec_mode=$($timingRow.tabi_matvec_mode)"
    }

    Write-Host ("PBE TABI smoke PASS: {0} sum={1} panels={2} E_kJ={3}" -f
        $Name, $SumMode, $row.grid_points, $row.E_rxn_kJmol)
    return [pscustomobject]@{
        name = $Name
        sum_mode = $SumMode
        formulation = $Formulation
        energy = $row
        timing = $timingRow
    }
}

function Invoke-InvalidBackendCase() {
    $runDir = Join-Path $WorkRoot 'invalid_backend'
    New-CleanDir $runDir
    $argPath = New-ArgFile -RunDir $runDir -Backend 'invalid' -SumMode 'direct' -Formulation 'jcp2013'
    $outPath = Join-Path $runDir 'invalid_backend.out'
    Push-Location $runDir
    try {
        try {
            & $TitanExe -a $argPath *> $outPath
        } catch {
            # PowerShell can surface nonzero native exits as NativeCommandError.
        }
        if ($LASTEXITCODE -eq 0) {
            throw 'Invalid PBE backend was unexpectedly accepted.'
        }
    } finally {
        Pop-Location
    }
    Write-Host 'PBE TABI smoke PASS: invalid backend rejected'
}

function Invoke-JcpTreecodeRejectCase() {
    $runDir = Join-Path $WorkRoot 'jcp_treecode_reject'
    New-CleanDir $runDir
    $argPath = New-ArgFile -RunDir $runDir -Backend 'tabi' -SumMode 'treecode' -Formulation 'jcp2013'
    $outPath = Join-Path $runDir 'jcp_treecode_reject.out'
    Push-Location $runDir
    try {
        try {
            & $TitanExe -a $argPath *> $outPath
        } catch {
            # PowerShell can surface nonzero native exits as NativeCommandError.
        }
        if ($LASTEXITCODE -eq 0) {
            throw 'jcp2013 TABI treecode was unexpectedly accepted.'
        }
    } finally {
        Pop-Location
    }
    Write-Host 'PBE TABI smoke PASS: jcp2013 treecode rejected'
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $InputRoot)) {
    throw "Missing PBE fixture directory: $InputRoot"
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$direct = Invoke-TabiCase -Name 'mini_protein_direct' -SumMode 'direct' -Formulation 'jcp2013'
$legacyDirect = Invoke-TabiCase -Name 'mini_protein_legacy_direct' -SumMode 'direct' -Formulation 'titan-v1'
$treecode = Invoke-TabiCase -Name 'mini_protein_legacy_treecode' -SumMode 'treecode' -Formulation 'titan-v1'
$directEnergy = Convert-ToDouble $legacyDirect.energy.E_rxn_kJmol
$treecodeEnergy = Convert-ToDouble $treecode.energy.E_rxn_kJmol
$agreementTol = [math]::Max(1.0e-5, 2.5e-1 * [math]::Max(1.0, [math]::Abs($directEnergy)))
if ([math]::Abs($directEnergy - $treecodeEnergy) -gt $agreementTol) {
    throw "TABI direct/treecode disagreement: direct=$directEnergy treecode=$treecodeEnergy"
}

$benchmarkPath = Join-Path $WorkRoot 'pbe_tabi_benchmark.csv'
@($direct, $legacyDirect, $treecode) | ForEach-Object {
    [pscustomobject]@{
        case = $_.name
        sum_mode = $_.sum_mode
        formulation = $_.formulation
        panel_count = $_.energy.grid_points
        energy_kj_mol = $_.energy.E_rxn_kJmol
        surface_build_seconds = $_.timing.tabi_surface_build_seconds
        tree_build_seconds = $_.timing.tabi_tree_build_seconds
        tree_nodes = $_.timing.tabi_tree_nodes
        tree_leaf_count = $_.timing.tabi_tree_leaf_count
        near_interactions = $_.timing.tabi_near_interactions
        far_interactions = $_.timing.tabi_far_interactions
        matvec_mode = $_.timing.tabi_matvec_mode
        surface_geometry = $_.timing.tabi_surface_geometry
        mesh_vertices = $_.timing.tabi_mesh_vertices
        mesh_triangles = $_.timing.tabi_mesh_triangles
        boundary_nodes = $_.timing.tabi_boundary_nodes
        solid_angle_fallbacks = $_.timing.tabi_solid_angle_fallbacks
        diag_min = $_.timing.tabi_diag_min
        diag_max = $_.timing.tabi_diag_max
        taylor_order = $_.timing.tabi_taylor_order
        tree_eval_seconds = $_.timing.tabi_tree_eval_seconds
    }
} | Export-Csv -LiteralPath $benchmarkPath -NoTypeInformation -Encoding ASCII
Invoke-JcpTreecodeRejectCase
Invoke-InvalidBackendCase
Write-Host 'PBE TABI smoke PASS: all cases'
