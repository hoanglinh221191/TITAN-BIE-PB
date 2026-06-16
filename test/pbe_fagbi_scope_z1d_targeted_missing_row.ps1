param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z1d_targeted'),
    [string]$LibFfRoot = $env:TITAN_LIB_FF,
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 5400,
    [int]$Threads = 1
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z1d_targeted_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z1d_targeted_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z1d_targeted_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z1D_TARGETED_VERIFY_REPORT.txt'
$TemperatureK = 298.15
$InternalDielectric = 1.0
$ExternalDielectric = 80.0
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function F([double]$Value) { return $Value.ToString('G17', $InvariantCulture) }

function Get-Z1dTargetProp([object]$Object, [string]$Name, [object]$Default = '') {
    if ($null -eq $Object) { return $Default }
    if ($Object.PSObject.Properties.Name -contains $Name) {
        $value = $Object.$Name
        if ($null -ne $value -and -not [string]::IsNullOrWhiteSpace([string]$value)) { return $value }
    }
    return $Default
}

function Read-Z1dTargetFirstCsv([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    $rows = @(Import-Csv -LiteralPath $Path)
    if ($rows.Count -eq 0) { return $null }
    return $rows[0]
}

function Write-Z1dTargetArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [int]$RequestedPoints,
    [double]$SaltM,
    [string]$LibPathForArg
) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdbPath
    $argLines = @(
        "-f=$(Convert-FmmCliToTitanPath $pdbPath)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$LibPathForArg",
        "-path-n=$LibPathForArg",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-FmmCliToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-probe=0.0',
        "-pbe-salt=$(F $SaltM)",
        "-pbe-temp=$(F $TemperatureK)",
        "-pbe-int-diel=$(F $InternalDielectric)",
        "-pbe-ext-diel=$(F $ExternalDielectric)",
        '-pbe-timing=detail',
        '-pbe-map=off',
        '-pbe-margin=2.0',
        '-pbe-max-iter=1200',
        '-pbe-rms=1e-6',
        '-pbe-dphi=1e-6',
        '-pbe-tabi-sum=direct',
        '-pbe-tabi-formulation=jcp2013',
        '-pbe-tabi-fmm-m2l=off',
        '-pbe-tabi-energy-projection=boundary',
        '-pbe-tabi-geometry=triangle',
        '-pbe-tabi-near-integration=analytic_singular',
        '-pbe-tabi-quad-buffer=off',
        "-pbe-tabi-surface-points=$RequestedPoints",
        '-pbe-tabi-highres-surface=true',
        '-pbe-tabi-sphere-mesh=icosphere',
        '-pbe-tabi-precond=diagonal',
        "-nt=$Threads"
    )
    $badOptions = @($argLines | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
    })
    if ($badOptions.Count -ne 0) {
        throw "Generated arg.in contains non equals-syntax options: $($badOptions -join '; ')"
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Read-Z1dTargetCaseResult([string]$RunDir) {
    $energy = Read-Z1dTargetFirstCsv (Join-Path $RunDir 'pbe_energy.csv')
    $profile = Read-Z1dTargetFirstCsv (Join-Path $RunDir 'pbe_tabi_profile.csv')
    $quality = Read-Z1dTargetFirstCsv (Join-Path $RunDir 'pbe_tabi_surface_quality.csv')
    $memory = Read-Z1dTargetFirstCsv (Join-Path $RunDir 'pbe_memory_profile.csv')
    if ($null -eq $energy) { throw "Missing pbe_energy.csv in $RunDir" }
    if ($null -eq $profile) { throw "Missing pbe_tabi_profile.csv in $RunDir" }
    if ($null -eq $quality) { throw "Missing pbe_tabi_surface_quality.csv in $RunDir" }
    $nearMode = [string](Get-Z1dTargetProp $profile 'near_integration' '')
    if ($nearMode -ne 'analytic_singular') { throw "Expected analytic_singular, got '$nearMode'" }
    $meshMode = [string](Get-Z1dTargetProp $quality 'sphere_mesh_mode' '')
    if ($meshMode -ne 'icosphere') { throw "Expected icosphere, got '$meshMode'" }
    $status = if (
        [string](Get-Z1dTargetProp $energy 'status' 'unknown') -eq 'ok' -and
        [string](Get-Z1dTargetProp $profile 'status' 'unknown') -eq 'ok' -and
        [string](Get-Z1dTargetProp $quality 'status' 'unknown') -eq 'ok') {
        'ok'
    } else {
        "energy=$(Get-Z1dTargetProp $energy 'status' 'unknown');profile=$(Get-Z1dTargetProp $profile 'status' 'unknown');quality=$(Get-Z1dTargetProp $quality 'status' 'unknown')"
    }
    return [pscustomobject]@{
        status = $status
        energy_kj_mol = Convert-FmmCliDouble (Get-Z1dTargetProp $energy 'E_rxn_kJmol' 'NA')
        actual_panels = Convert-FmmCliInt64 (Get-Z1dTargetProp $quality 'panel_count' 0)
        iterations = Convert-FmmCliInt64 (Get-Z1dTargetProp $profile 'iterations' 0)
        final_residual = Convert-FmmCliDouble (Get-Z1dTargetProp $profile 'final_residual' 'NA')
        true_residual = Convert-FmmCliDouble (Get-Z1dTargetProp $profile 'true_residual_rel' 'NA')
        total_seconds = Convert-FmmCliDouble (Get-Z1dTargetProp $profile 'total_seconds' 'NA')
        near_fallbacks = Convert-FmmCliInt64 (Get-Z1dTargetProp $profile 'near_fallbacks' 0)
        planar_area_rel_error = Convert-FmmCliDouble (Get-Z1dTargetProp $quality 'planar_area_rel_error_vs_4piR2' 'NA')
        enclosed_volume_rel_error = Convert-FmmCliDouble (Get-Z1dTargetProp $quality 'enclosed_volume_rel_error_vs_4piR3_over3' 'NA')
        triangle_area_cv = Convert-FmmCliDouble (Get-Z1dTargetProp $quality 'triangle_area_cv' 'NA')
        triangle_aspect_max = Convert-FmmCliDouble (Get-Z1dTargetProp $quality 'triangle_aspect_max' 'NA')
        peak_memory_mb = Convert-FmmCliDouble (Get-Z1dTargetProp $memory 'peak_working_set_mb' 'NA')
    }
}

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope Z1D targeted work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) { Remove-Item -LiteralPath $workRootFull -Recurse -Force }
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

$resolution = Get-FagbiHighresSurfaceResolution -TargetPanels 20480
$runDir = Join-Path $workRootFull 'kirkwood_offcenter_salt0p15_target20480'
New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
$pqrPath = Join-Path $runDir 'input.pqr'
Write-FagbiKirkwoodPqr -Path $pqrPath -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0

$defaultLibFfRoot = Join-Path $RepoRoot 'lib_ff'
if ([string]::IsNullOrWhiteSpace($LibFfRoot)) { $LibFfRoot = $defaultLibFfRoot }
$libPathForArg = Convert-FmmCliToTitanPath ($LibFfRoot.TrimEnd('\') + '\')
$argPath = Write-Z1dTargetArgFile `
    -RunDir $runDir `
    -PqrPath $pqrPath `
    -RequestedPoints $resolution.requested_points `
    -SaltM 0.15 `
    -LibPathForArg $libPathForArg

$oracle = Get-FagbiKirkwoodEnergyKJ `
    -SphereRadiusA 10.0 `
    -Charge 1.0 `
    -OffsetA 3.0 `
    -SaltM 0.15 `
    -TemperatureK $TemperatureK `
    -InternalDielectric $InternalDielectric `
    -ExternalDielectric $ExternalDielectric

$result = [pscustomobject]@{
    status = 'quick_arg_shape_only'
    energy_kj_mol = [double]::NaN
    actual_panels = 0
    iterations = 0
    final_residual = [double]::NaN
    true_residual = [double]::NaN
    total_seconds = [double]::NaN
    near_fallbacks = 0
    planar_area_rel_error = [double]::NaN
    enclosed_volume_rel_error = [double]::NaN
    triangle_area_cv = [double]::NaN
    triangle_aspect_max = [double]::NaN
    peak_memory_mb = [double]::NaN
}
$run = [pscustomobject]@{ exit_code = 0; timed_out = $false; wall_seconds = 0.0; peak_working_set_mb = [double]::NaN }
$readError = ''

if (-not $Quick.IsPresent) {
    if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
    $TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
    if (-not (Test-Path -LiteralPath $LibFfRoot)) {
        throw "Missing lib_ff directory: $LibFfRoot. Set -LibFfRoot or TITAN_LIB_FF for linked worktrees."
    }
    $run = Invoke-FmmCliTitanCase `
        -TitanExe $TitanExe `
        -RunDir $runDir `
        -ArgPath $argPath `
        -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -eq 0 -and -not $run.timed_out) {
        try {
            $result = Read-Z1dTargetCaseResult -RunDir $runDir
        } catch {
            $readError = $_.Exception.Message
            $result = [pscustomobject]@{
                status = 'read_error'
                energy_kj_mol = [double]::NaN
                actual_panels = 0
                iterations = 0
                final_residual = [double]::NaN
                true_residual = [double]::NaN
                total_seconds = [double]::NaN
                near_fallbacks = 0
                planar_area_rel_error = [double]::NaN
                enclosed_volume_rel_error = [double]::NaN
                triangle_area_cv = [double]::NaN
                triangle_aspect_max = [double]::NaN
                peak_memory_mb = [double]::NaN
            }
        }
    } elseif ($run.timed_out) {
        $result.status = 'timeout'
    } else {
        $result.status = "exit_$($run.exit_code)"
    }
}

$oracleRelativeError = if ($result.status -eq 'ok' -and [math]::Abs($oracle) -gt 1.0e-12) {
    [math]::Abs([double]$result.energy_kj_mol - $oracle) / [math]::Abs($oracle)
} else {
    [double]::NaN
}

$decision = if ($Quick.IsPresent) {
    'z1d_targeted_quick_ready'
} elseif ($result.status -eq 'ok' -and (Test-FmmCliFinite $oracleRelativeError) -and $oracleRelativeError -le 5.0e-3 -and [int64]$result.near_fallbacks -eq 0) {
    'z1d_resolution_converged'
} elseif ($result.status -eq 'ok') {
    'z1d_operator_or_resolution_gap'
} elseif ($result.status -eq 'timeout') {
    'z1d_runtime_limited_confirmed'
} else {
    'z1d_targeted_failed'
}

$matrixRow = [pscustomobject]@{
    fixture = 'kirkwood_offcenter'
    salt_M = 0.15
    mesh_mode = 'icosphere'
    target_panels = 20480
    requested_points = $resolution.requested_points
    expected_panels = $resolution.expected_panels
    actual_panels = $result.actual_panels
    status = $result.status
    decision = $decision
    read_error = $readError
    exit_code = $run.exit_code
    timed_out = $run.timed_out
    energy_kj_mol = $result.energy_kj_mol
    oracle_energy_kj_mol = $oracle
    oracle_relative_error = $oracleRelativeError
    iterations = $result.iterations
    final_residual = $result.final_residual
    true_residual = $result.true_residual
    wall_seconds = $run.wall_seconds
    solver_total_seconds = $result.total_seconds
    peak_memory_mb = $result.peak_memory_mb
    process_peak_memory_mb = $run.peak_working_set_mb
    near_fallbacks = $result.near_fallbacks
    planar_area_rel_error = $result.planar_area_rel_error
    enclosed_volume_rel_error = $result.enclosed_volume_rel_error
    triangle_area_cv = $result.triangle_area_cv
    triangle_aspect_max = $result.triangle_aspect_max
    arg_path = $argPath
    run_dir = $runDir
}

$matrixRow | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
[pscustomobject]@{
    decision = $decision
    quick = $Quick.IsPresent
    status = $result.status
    oracle_relative_error = $oracleRelativeError
    wall_seconds = $run.wall_seconds
    matrix = $MatrixCsv
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    "decision: $decision",
    "quick: $($Quick.IsPresent)",
    "fixture: kirkwood_offcenter",
    "salt_M: 0.15",
    "target_panels: 20480",
    "requested_points: $($resolution.requested_points)",
    "status: $($result.status)",
    "oracle_relative_error: $oracleRelativeError",
    "wall_seconds: $($run.wall_seconds)",
    "matrix: $MatrixCsv"
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    'TITAN FMM/M2L Scope Z1D.1 Targeted Missing-Row Verification Report',
    "Decision: $decision",
    "Generated: $((Get-Date).ToString('s'))",
    '',
    'Target row',
    '- Fixture: off-center Kirkwood sphere.',
    '- Salt: 0.15 M.',
    '- Mesh: icosphere, target 20480 panels, requested surface points 5121.',
    '- Backend: direct TABI, analytic_singular near integration, FMM/M2L off, nt=1.',
    '',
    'Measured row',
    "- Status: $($result.status).",
    "- Oracle relative error: $oracleRelativeError.",
    "- Wall seconds: $($run.wall_seconds).",
    "- Near fallbacks: $($result.near_fallbacks).",
    "- Actual panels: $($result.actual_panels).",
    '',
    'Interpretation'
)
switch ($decision) {
    'z1d_resolution_converged' {
        $report += 'The missing high-resolution row crossed the 0.5% oracle-error gate with no near fallback; Z1D can be upgraded after syncing ledgers.'
    }
    'z1d_operator_or_resolution_gap' {
        $report += 'The missing row completed but did not satisfy the 0.5%/no-fallback gate; keep Z1D as an operator or resolution gap.'
    }
    'z1d_runtime_limited_confirmed' {
        $report += 'The missing row timed out again; high-resolution direct runtime remains the Z1 blocker.'
    }
    'z1d_targeted_quick_ready' {
        $report += 'Quick mode generated artifacts and arg.in shape only; no expensive solve was run.'
    }
    default {
        $report += 'The targeted row failed before producing a scientific classification; inspect stdout/stderr in the run directory.'
    }
}
$report += @(
    '',
    "Matrix CSV: $MatrixCsv",
    "Decision CSV: $DecisionCsv",
    "Summary: $SummaryPath"
)
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "Scope Z1D.1 decision: $decision"
Write-Host "Matrix: $MatrixCsv"
