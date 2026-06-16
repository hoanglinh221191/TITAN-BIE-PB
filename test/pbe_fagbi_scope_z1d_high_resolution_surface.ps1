param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z1d'),
    [string]$LibFfRoot = $env:TITAN_LIB_FF,
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 1800,
    [int]$Threads = 1
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z1d_matrix.csv'
$ConvergenceCsv = Join-Path $WorkRoot 'fagbi_scope_z1d_convergence.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z1d_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z1d_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z1d_highres_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z1D_HIGH_RESOLUTION_REPORT.txt'
$TemperatureK = 298.15
$InternalDielectric = 1.0
$ExternalDielectric = 80.0
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$defaultLibFfRoot = Join-Path $RepoRoot 'lib_ff'
if ([string]::IsNullOrWhiteSpace($LibFfRoot)) { $LibFfRoot = $defaultLibFfRoot }
if (-not (Test-Path -LiteralPath $LibFfRoot)) {
    throw "Missing lib_ff directory: $LibFfRoot. Set -LibFfRoot or TITAN_LIB_FF for linked worktrees."
}
$LibFfRoot = (Resolve-Path -LiteralPath $LibFfRoot).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope Z1D work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) { Remove-Item -LiteralPath $workRootFull -Recurse -Force }
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function F([double]$Value) { return $Value.ToString('G17', $InvariantCulture) }

function Get-Z1dProp([object]$Object, [string]$Name, [object]$Default = '') {
    if ($null -eq $Object) { return $Default }
    if ($Object.PSObject.Properties.Name -contains $Name) {
        $value = $Object.$Name
        if ($null -ne $value -and -not [string]::IsNullOrWhiteSpace([string]$value)) { return $value }
    }
    return $Default
}

function Read-Z1dFirstCsv([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    $rows = @(Import-Csv -LiteralPath $Path)
    if ($rows.Count -eq 0) { return $null }
    return $rows[0]
}

function Write-Z1dArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [int]$RequestedPoints,
    [double]$SaltM
) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdbPath
    $libPath = Convert-FmmCliToTitanPath ($LibFfRoot.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-FmmCliToTitanPath $pdbPath)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
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

function New-Z1dEmptyResult([string]$Status) {
    return [pscustomobject]@{
        status=$Status; energy_kj_mol=[double]::NaN; actual_panels=0
        iterations=0; final_residual=[double]::NaN; true_residual=[double]::NaN
        total_seconds=[double]::NaN; near_fallbacks=0
        planar_area_rel_error=[double]::NaN; enclosed_volume_rel_error=[double]::NaN
        triangle_area_cv=[double]::NaN; triangle_aspect_max=[double]::NaN
        peak_memory_mb=[double]::NaN
    }
}

function Read-Z1dCaseResult([string]$RunDir) {
    $energy = Read-Z1dFirstCsv (Join-Path $RunDir 'pbe_energy.csv')
    $profile = Read-Z1dFirstCsv (Join-Path $RunDir 'pbe_tabi_profile.csv')
    $quality = Read-Z1dFirstCsv (Join-Path $RunDir 'pbe_tabi_surface_quality.csv')
    $memory = Read-Z1dFirstCsv (Join-Path $RunDir 'pbe_memory_profile.csv')
    if ($null -eq $energy) { throw "Missing pbe_energy.csv in $RunDir" }
    if ($null -eq $profile) { throw "Missing pbe_tabi_profile.csv in $RunDir" }
    if ($null -eq $quality) { throw "Missing pbe_tabi_surface_quality.csv in $RunDir" }
    $nearMode = [string](Get-Z1dProp $profile 'near_integration' '')
    if ($nearMode -ne 'analytic_singular') { throw "Expected analytic_singular, got '$nearMode'" }
    $meshMode = [string](Get-Z1dProp $quality 'sphere_mesh_mode' '')
    if ($meshMode -ne 'icosphere') { throw "Expected icosphere, got '$meshMode'" }
    $energyStatus = [string](Get-Z1dProp $energy 'status' 'unknown')
    $profileStatus = [string](Get-Z1dProp $profile 'status' 'unknown')
    $qualityStatus = [string](Get-Z1dProp $quality 'status' 'unknown')
    $status = if ($energyStatus -eq 'ok' -and $profileStatus -eq 'ok' -and $qualityStatus -eq 'ok') { 'ok' } else { "energy=$energyStatus;profile=$profileStatus;quality=$qualityStatus" }
    return [pscustomobject]@{
        status=$status
        energy_kj_mol=Convert-FmmCliDouble (Get-Z1dProp $energy 'E_rxn_kJmol' 'NA')
        actual_panels=Convert-FmmCliInt64 (Get-Z1dProp $quality 'panel_count' 0)
        iterations=Convert-FmmCliInt64 (Get-Z1dProp $profile 'iterations' 0)
        final_residual=Convert-FmmCliDouble (Get-Z1dProp $profile 'final_residual' 'NA')
        true_residual=Convert-FmmCliDouble (Get-Z1dProp $profile 'true_residual_rel' 'NA')
        total_seconds=Convert-FmmCliDouble (Get-Z1dProp $profile 'total_seconds' 'NA')
        near_fallbacks=Convert-FmmCliInt64 (Get-Z1dProp $profile 'near_fallbacks' 0)
        planar_area_rel_error=Convert-FmmCliDouble (Get-Z1dProp $quality 'planar_area_rel_error_vs_4piR2' 'NA')
        enclosed_volume_rel_error=Convert-FmmCliDouble (Get-Z1dProp $quality 'enclosed_volume_rel_error_vs_4piR3_over3' 'NA')
        triangle_area_cv=Convert-FmmCliDouble (Get-Z1dProp $quality 'triangle_area_cv' 'NA')
        triangle_aspect_max=Convert-FmmCliDouble (Get-Z1dProp $quality 'triangle_aspect_max' 'NA')
        peak_memory_mb=Convert-FmmCliDouble (Get-Z1dProp $memory 'peak_working_set_mb' 'NA')
    }
}

function Invoke-Z1dCase([object]$Case) {
    $safeSalt = (F $Case.salt_M) -replace '-', 'm' -replace '\.', 'p'
    $caseName = "{0}_salt{1}_target{2}" -f $Case.fixture,$safeSalt,$Case.target_panels
    $runDir = Join-Path $workRootFull (Get-FmmCliSafeName $caseName)
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $pqrPath = Join-Path $runDir 'input.pqr'
    if ($Case.fixture -eq 'born_centered') {
        Write-FagbiBornPqr -Path $pqrPath -Charge 1.0 -RadiusA 3.0
    } else {
        Write-FagbiKirkwoodPqr -Path $pqrPath -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0
    }
    $argPath = Write-Z1dArgFile -RunDir $runDir -PqrPath $pqrPath -RequestedPoints $Case.requested_points -SaltM $Case.salt_M
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $readError = ''
    $result = $null
    if ($run.exit_code -eq 0 -and -not $run.timed_out) {
        try { $result = Read-Z1dCaseResult -RunDir $runDir } catch { $readError = $_.Exception.Message }
    }
    if ($null -eq $result) {
        $status = if ($run.timed_out) { 'timeout' } elseif (-not [string]::IsNullOrWhiteSpace($readError)) { 'read_error' } else { "exit_$($run.exit_code)" }
        $result = New-Z1dEmptyResult -Status $status
    }
    return [pscustomobject]@{ case=$Case; run=$run; result=$result; read_error=$readError; run_dir=$runDir }
}

function Get-Z1dImprovement([double]$Coarse, [double]$Fine) {
    if (-not (Test-FmmCliFinite $Coarse) -or -not (Test-FmmCliFinite $Fine) -or $Fine -le 0.0) { return [double]::NaN }
    return $Coarse / $Fine
}

function Get-Z1dMedian([double[]]$Values) {
    $finite = @($Values | Where-Object { Test-FmmCliFinite $_ } | Sort-Object)
    if ($finite.Count -eq 0) { return [double]::NaN }
    $middle = [int][math]::Floor($finite.Count / 2)
    if (($finite.Count % 2) -eq 1) { return [double]$finite[$middle] }
    return 0.5 * ([double]$finite[$middle - 1] + [double]$finite[$middle])
}

function Write-Z1dSvgPlot([object[]]$Rows, [string]$Path) {
    $width = 860; $height = 520; $left = 80; $top = 60; $plotWidth = 680; $plotHeight = 340
    $validRows = @($Rows | Where-Object { $_.status -eq 'ok' -and (Test-FmmCliFinite ([double]$_.oracle_relative_error)) })
    if ($validRows.Count -eq 0) {
        Set-Content -LiteralPath $Path -Value "<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height'><text x='20' y='40'>No valid Z1D rows</text></svg>" -Encoding ASCII
        return
    }
    $minX = [math]::Log10(5120.0); $maxX = [math]::Log10(20480.0)
    $minY = -3.0; $maxY = -1.7
    function X([double]$panels) { return $left + (($([math]::Log10($panels)) - $minX) / ($maxX - $minX)) * $plotWidth }
    function Y([double]$err) {
        $log = [math]::Log10([math]::Max(1.0e-12, $err))
        return $top + (1.0 - (($log - $minY) / ($maxY - $minY))) * $plotHeight
    }
    $colors = @{ 'born_centered|0'='#0072b2'; 'born_centered|0.15'='#56b4e9'; 'kirkwood_offcenter|0'='#d55e00'; 'kirkwood_offcenter|0.15'='#cc79a7' }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add("<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
    $svg.Add("<rect width='100%' height='100%' fill='white'/>")
    $svg.Add("<text x='24' y='32' font-size='18' font-family='Arial'>Scope Z1D high-resolution surface convergence</text>")
    $svg.Add("<line x1='$left' y1='$top' x2='$left' y2='$($top+$plotHeight)' stroke='black'/>")
    $svg.Add("<line x1='$left' y1='$($top+$plotHeight)' x2='$($left+$plotWidth)' y2='$($top+$plotHeight)' stroke='black'/>")
    $svg.Add("<text x='360' y='475' font-size='13' font-family='Arial'>actual panels</text>")
    $svg.Add("<text x='12' y='250' font-size='13' font-family='Arial' transform='rotate(-90 12,250)'>relative oracle error</text>")
    foreach ($group in ($validRows | Group-Object fixture,salt_M)) {
        $ordered = @($group.Group | Sort-Object {[double]$_.actual_panels})
        $key = "$($ordered[0].fixture)|$($ordered[0].salt_M)"
        $color = if ($colors.ContainsKey($key)) { $colors[$key] } else { '#333333' }
        $points = @($ordered | ForEach-Object { "$(X ([double]$_.actual_panels)),$(Y ([double]$_.oracle_relative_error))" })
        if ($points.Count -ge 2) { $svg.Add("<polyline fill='none' stroke='$color' stroke-width='2' points='$($points -join ' ')'/>") }
        foreach ($row in $ordered) {
            $svg.Add("<circle cx='$(X ([double]$row.actual_panels))' cy='$(Y ([double]$row.oracle_relative_error))' r='4' fill='$color'/>")
        }
    }
    $legendY = 80
    foreach ($key in $colors.Keys) {
        $color = $colors[$key]
        $svg.Add("<rect x='705' y='$legendY' width='12' height='12' fill='$color'/>")
        $svg.Add("<text x='724' y='$($legendY+11)' font-size='11' font-family='Arial'>$key</text>")
        $legendY += 20
    }
    $svg.Add("</svg>")
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$caseSpecs = @(Get-FagbiZ1dCaseSpecs -Quick:$Quick.IsPresent)
$rows = New-Object System.Collections.Generic.List[object]
foreach ($case in $caseSpecs) {
    $caseResult = Invoke-Z1dCase -Case $case
    $oracle = if ($case.fixture -eq 'born_centered') {
        Get-FagbiBornEnergyKJ -Charge 1.0 -RadiusA 3.0 -SaltM $case.salt_M -TemperatureK $TemperatureK -InternalDielectric $InternalDielectric -ExternalDielectric $ExternalDielectric
    } else {
        Get-FagbiKirkwoodEnergyKJ -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0 -SaltM $case.salt_M -TemperatureK $TemperatureK -InternalDielectric $InternalDielectric -ExternalDielectric $ExternalDielectric
    }
    $result = $caseResult.result
    $energyErr = if ($result.status -eq 'ok' -and [math]::Abs($oracle) -gt 1.0e-12) { [math]::Abs([double]$result.energy_kj_mol - $oracle) / [math]::Abs($oracle) } else { [double]::NaN }
    $rows.Add([pscustomobject]@{
        fixture=$case.fixture; salt_M=$case.salt_M; mesh_mode='icosphere'; target_panels=$case.target_panels
        requested_points=$case.requested_points; expected_panels=$case.expected_panels; actual_panels=$result.actual_panels
        status=$result.status; read_error=$caseResult.read_error; exit_code=$caseResult.run.exit_code; timed_out=$caseResult.run.timed_out
        energy_kj_mol=$result.energy_kj_mol; oracle_energy_kj_mol=$oracle; oracle_relative_error=$energyErr
        iterations=$result.iterations; final_residual=$result.final_residual; true_residual=$result.true_residual
        wall_seconds=$caseResult.run.wall_seconds; solver_total_seconds=$result.total_seconds; peak_memory_mb=$result.peak_memory_mb
        near_fallbacks=$result.near_fallbacks; planar_area_rel_error=$result.planar_area_rel_error
        enclosed_volume_rel_error=$result.enclosed_volume_rel_error; triangle_area_cv=$result.triangle_area_cv
        triangle_aspect_max=$result.triangle_aspect_max; run_dir=$caseResult.run_dir
    }) | Out-Null
}

$convergenceRows = New-Object System.Collections.Generic.List[object]
foreach ($group in ($rows | Group-Object fixture,salt_M)) {
    $ordered = @($group.Group | Sort-Object {[double]$_.actual_panels})
    for ($i=1; $i -lt $ordered.Count; ++$i) {
        $coarse = $ordered[$i-1]; $fine = $ordered[$i]
        $convergenceRows.Add([pscustomobject]@{
            fixture=$fine.fixture; salt_M=$fine.salt_M; coarse_target=$coarse.target_panels; fine_target=$fine.target_panels
            coarse_actual=$coarse.actual_panels; fine_actual=$fine.actual_panels
            energy_improvement=Get-Z1dImprovement ([double]$coarse.oracle_relative_error) ([double]$fine.oracle_relative_error)
            planar_area_improvement=Get-Z1dImprovement ([double]$coarse.planar_area_rel_error) ([double]$fine.planar_area_rel_error)
            volume_improvement=Get-Z1dImprovement ([double]$coarse.enclosed_volume_rel_error) ([double]$fine.enclosed_volume_rel_error)
            fine_oracle_relative_error=$fine.oracle_relative_error; status=$fine.status
        }) | Out-Null
    }
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$convergenceRows | Export-Csv -LiteralPath $ConvergenceCsv -NoTypeInformation
Write-Z1dSvgPlot -Rows $rows -Path $PlotPath

$invalidRows = @($rows | Where-Object { $_.status -ne 'ok' })
$finestRows = @($rows | Where-Object { $_.target_panels -eq 20480 -and $_.status -eq 'ok' })
$finestExpectedRows = @($rows | Where-Object { $_.target_panels -eq 20480 })
$transitionRows = @($convergenceRows | Where-Object { $_.fine_target -eq 20480 -and $_.status -eq 'ok' })
$allFinestBelowHalfPercent = $invalidRows.Count -eq 0 -and $finestRows.Count -eq (@($caseSpecs | Group-Object fixture,salt_M)).Count -and @($finestRows | Where-Object { [double]$_.oracle_relative_error -gt 5.0e-3 }).Count -eq 0
$resolutionLimitedGroups = @($transitionRows | Where-Object { [double]$_.energy_improvement -ge 1.5 -and [double]$_.fine_oracle_relative_error -gt 5.0e-3 })
$plateauGroups = @($transitionRows | Where-Object {
    [double]$_.planar_area_improvement -ge 2.0 -and [double]$_.volume_improvement -ge 2.0 -and
    [double]$_.energy_improvement -lt 1.25 -and [double]$_.fine_oracle_relative_error -gt 1.0e-2
})
$accuracyClass = if ($invalidRows.Count -gt 0) {
    'z1d_runtime_limited'
} elseif ($allFinestBelowHalfPercent) {
    'z1d_resolution_converged'
} elseif ($transitionRows.Count -gt 0 -and $resolutionLimitedGroups.Count -ge [math]::Ceiling(0.75*$transitionRows.Count)) {
    'z1d_resolution_limited'
} elseif ($transitionRows.Count -gt 0 -and $plateauGroups.Count -ge [math]::Ceiling(0.5*$transitionRows.Count)) {
    'z1d_operator_plateau'
} else {
    'z1d_inconclusive'
}
$medianFinestError = Get-Z1dMedian @($finestRows | ForEach-Object { $_.oracle_relative_error })

[pscustomobject]@{
    accuracy_class=$accuracyClass; quick=$Quick.IsPresent; rows=$rows.Count; invalid_rows=$invalidRows.Count
    finest_target_panels=20480; resolution_limited_groups=$resolutionLimitedGroups.Count
    plateau_groups=$plateauGroups.Count; median_finest_oracle_error=$medianFinestError
    matrix=$MatrixCsv; convergence=$ConvergenceCsv; plot=$PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    "accuracy_class: $accuracyClass",
    "quick: $($Quick.IsPresent)",
    "rows: $($rows.Count)",
    "invalid_rows: $($invalidRows.Count)",
    "median_finest_oracle_error: $medianFinestError",
    "matrix: $MatrixCsv",
    "convergence: $ConvergenceCsv",
    "plot: $PlotPath"
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = New-Object System.Collections.Generic.List[string]
$report.Add('TITAN FMM/M2L Scope Z1D High-Resolution Surface Report')
$report.Add("Accuracy class: $accuracyClass")
$report.Add("Generated: $((Get-Date).ToString('s'))")
$report.Add('')
$report.Add('Purpose')
$report.Add('Extend Z1C beyond approximately 5120 panels using direct TABI and analytic_singular near integration.')
$report.Add('')
$report.Add('Matrix')
$report.Add('- Fixtures: centered Born and off-center Kirkwood.')
$report.Add('- Backend: TABI Direct; FMM/M2L disabled; analytic_singular near integration.')
$report.Add('- Mesh: icosphere only; target panels: 5120,20480.')
$report.Add('- Salts: 0,0.15 M.')
$report.Add('')
$report.Add('Measured classification')
$report.Add("- Invalid rows: $($invalidRows.Count).")
$report.Add("- Resolution-limited transitions: $($resolutionLimitedGroups.Count)/$($transitionRows.Count).")
$report.Add("- Operator-plateau transitions: $($plateauGroups.Count)/$($transitionRows.Count).")
$report.Add("- Median finest oracle error: $medianFinestError.")
$report.Add("- Completed 20480-panel rows: $($finestRows.Count)/$($finestExpectedRows.Count).")
foreach ($row in ($finestRows | Sort-Object fixture,salt_M)) {
    $report.Add(("- Completed 20480-panel oracle error: {0} salt {1} = {2:P4}." -f $row.fixture, $row.salt_M, [double]$row.oracle_relative_error))
}
foreach ($row in ($finestExpectedRows | Where-Object { $_.status -ne 'ok' } | Sort-Object fixture,salt_M)) {
    $report.Add(("- Runtime-limited 20480-panel row: {0} salt {1}, status {2}, wall_seconds {3}." -f $row.fixture, $row.salt_M, $row.status, $row.wall_seconds))
}
$report.Add('')
$report.Add('Next decision')
switch ($accuracyClass) {
    'z1d_resolution_converged' { $report.Add('Proceed to Scope Z2A or Z3 with geometry error no longer the dominant observed blocker.') }
    'z1d_resolution_limited' { $report.Add('Consider one more resolution extension or a cheaper high-resolution verifier before changing operator formulation.') }
    'z1d_operator_plateau' { $report.Add('Review collocation, diagonal formulation, and panel basis before FMM near-field work.') }
    'z1d_runtime_limited' { $report.Add('Treat direct high-resolution solve cost as the next blocker; use geometry-only or sampled verifier before more direct solves.') }
    default { $report.Add('Keep Z1D evidence as inconclusive and inspect matrix rows before planning Z2.') }
}
$report.Add('')
$report.Add("Matrix CSV: $MatrixCsv")
$report.Add("Convergence CSV: $ConvergenceCsv")
$report.Add("SVG: $PlotPath")
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "Scope Z1D accuracy class: $accuracyClass"
Write-Host "Matrix: $MatrixCsv"
