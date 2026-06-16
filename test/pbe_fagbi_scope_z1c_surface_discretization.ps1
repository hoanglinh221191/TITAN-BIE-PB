param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z1c'),
    [string]$LibFfRoot = $env:TITAN_LIB_FF,
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 1
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z1c_matrix.csv'
$ConvergenceCsv = Join-Path $WorkRoot 'fagbi_scope_z1c_convergence.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z1c_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z1c_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z1c_geometry_energy_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z1C_SURFACE_DISCRETIZATION_REPORT.txt'
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
    throw "Refusing to clean Scope Z1C work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) { Remove-Item -LiteralPath $workRootFull -Recurse -Force }
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function F([double]$Value) { return $Value.ToString('G17', $InvariantCulture) }

function Get-Z1cProp([object]$Object, [string]$Name, [object]$Default = '') {
    if ($null -eq $Object) { return $Default }
    if ($Object.PSObject.Properties.Name -contains $Name) {
        $value = $Object.$Name
        if ($null -ne $value -and -not [string]::IsNullOrWhiteSpace([string]$value)) { return $value }
    }
    return $Default
}

function Read-Z1cFirstCsv([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    $rows = @(Import-Csv -LiteralPath $Path)
    if ($rows.Count -eq 0) { return $null }
    return $rows[0]
}

function Write-Z1cArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [int]$RequestedPoints,
    [double]$SaltM,
    [ValidateSet('icosphere','latlong')][string]$MeshMode
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
        "-pbe-tabi-sphere-mesh=$MeshMode",
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

function New-Z1cEmptyResult([string]$Status) {
    return [pscustomobject]@{
        status=$Status; energy_kj_mol=[double]::NaN; actual_panels=0
        iterations=0; final_residual=[double]::NaN; true_residual=[double]::NaN
        total_seconds=[double]::NaN; near_fallbacks=0; mesh_mode=''
        spherical_area=[double]::NaN; spherical_area_rel_error=[double]::NaN
        planar_area=[double]::NaN; planar_area_rel_error=[double]::NaN
        enclosed_volume=[double]::NaN; enclosed_volume_rel_error=[double]::NaN
        triangle_area_min=[double]::NaN; triangle_area_mean=[double]::NaN
        triangle_area_max=[double]::NaN; triangle_area_cv=[double]::NaN
        triangle_aspect_max=[double]::NaN; normal_radial_dot_min=[double]::NaN
        normal_radial_dot_mean=[double]::NaN; centroid_radial_error_max=[double]::NaN
    }
}

function Read-Z1cCaseResult([string]$RunDir, [string]$ExpectedMeshMode) {
    $energy = Read-Z1cFirstCsv (Join-Path $RunDir 'pbe_energy.csv')
    $profile = Read-Z1cFirstCsv (Join-Path $RunDir 'pbe_tabi_profile.csv')
    $quality = Read-Z1cFirstCsv (Join-Path $RunDir 'pbe_tabi_surface_quality.csv')
    if ($null -eq $energy) { throw "Missing pbe_energy.csv in $RunDir" }
    if ($null -eq $profile) { throw "Missing pbe_tabi_profile.csv in $RunDir" }
    if ($null -eq $quality) { throw "Missing pbe_tabi_surface_quality.csv in $RunDir" }
    $nearMode = [string](Get-Z1cProp $profile 'near_integration' '')
    if ($nearMode -ne 'analytic_singular') { throw "Expected analytic_singular, got '$nearMode'" }
    $meshMode = [string](Get-Z1cProp $quality 'sphere_mesh_mode' '')
    if ($meshMode -ne $ExpectedMeshMode) { throw "Expected mesh $ExpectedMeshMode, got '$meshMode'" }
    $energyStatus = [string](Get-Z1cProp $energy 'status' 'unknown')
    $profileStatus = [string](Get-Z1cProp $profile 'status' 'unknown')
    $qualityStatus = [string](Get-Z1cProp $quality 'status' 'unknown')
    $status = if ($energyStatus -eq 'ok' -and $profileStatus -eq 'ok' -and $qualityStatus -eq 'ok') { 'ok' } else { "energy=$energyStatus;profile=$profileStatus;quality=$qualityStatus" }
    return [pscustomobject]@{
        status=$status
        energy_kj_mol=Convert-FmmCliDouble (Get-Z1cProp $energy 'E_rxn_kJmol' 'NA')
        actual_panels=Convert-FmmCliInt64 (Get-Z1cProp $quality 'panel_count' 0)
        iterations=Convert-FmmCliInt64 (Get-Z1cProp $profile 'iterations' 0)
        final_residual=Convert-FmmCliDouble (Get-Z1cProp $profile 'final_residual' 'NA')
        true_residual=Convert-FmmCliDouble (Get-Z1cProp $profile 'true_residual_rel' 'NA')
        total_seconds=Convert-FmmCliDouble (Get-Z1cProp $profile 'total_seconds' 'NA')
        near_fallbacks=Convert-FmmCliInt64 (Get-Z1cProp $profile 'near_fallbacks' 0)
        mesh_mode=$meshMode
        spherical_area=Convert-FmmCliDouble (Get-Z1cProp $quality 'total_area' 'NA')
        spherical_area_rel_error=Convert-FmmCliDouble (Get-Z1cProp $quality 'area_rel_error_vs_4piR2' 'NA')
        planar_area=Convert-FmmCliDouble (Get-Z1cProp $quality 'planar_area' 'NA')
        planar_area_rel_error=Convert-FmmCliDouble (Get-Z1cProp $quality 'planar_area_rel_error_vs_4piR2' 'NA')
        enclosed_volume=Convert-FmmCliDouble (Get-Z1cProp $quality 'enclosed_volume' 'NA')
        enclosed_volume_rel_error=Convert-FmmCliDouble (Get-Z1cProp $quality 'enclosed_volume_rel_error_vs_4piR3_over3' 'NA')
        triangle_area_min=Convert-FmmCliDouble (Get-Z1cProp $quality 'triangle_area_min' 'NA')
        triangle_area_mean=Convert-FmmCliDouble (Get-Z1cProp $quality 'triangle_area_mean' 'NA')
        triangle_area_max=Convert-FmmCliDouble (Get-Z1cProp $quality 'triangle_area_max' 'NA')
        triangle_area_cv=Convert-FmmCliDouble (Get-Z1cProp $quality 'triangle_area_cv' 'NA')
        triangle_aspect_max=Convert-FmmCliDouble (Get-Z1cProp $quality 'triangle_aspect_max' 'NA')
        normal_radial_dot_min=Convert-FmmCliDouble (Get-Z1cProp $quality 'normal_radial_dot_min' 'NA')
        normal_radial_dot_mean=Convert-FmmCliDouble (Get-Z1cProp $quality 'normal_radial_dot_mean' 'NA')
        centroid_radial_error_max=Convert-FmmCliDouble (Get-Z1cProp $quality 'centroid_radial_error_max' 'NA')
    }
}

function Invoke-Z1cCase([object]$Case) {
    $safeSalt = (F $Case.salt_M) -replace '-', 'm' -replace '\.', 'p'
    $caseName = "{0}_salt{1}_{2}_target{3}" -f $Case.fixture,$safeSalt,$Case.mesh_mode,$Case.target_panels
    $runDir = Join-Path $workRootFull (Get-FmmCliSafeName $caseName)
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $pqrPath = Join-Path $runDir 'input.pqr'
    if ($Case.fixture -eq 'born_centered') {
        Write-FagbiBornPqr -Path $pqrPath -Charge 1.0 -RadiusA 3.0
    } else {
        Write-FagbiKirkwoodPqr -Path $pqrPath -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0
    }
    $argPath = Write-Z1cArgFile -RunDir $runDir -PqrPath $pqrPath -RequestedPoints $Case.requested_points -SaltM $Case.salt_M -MeshMode $Case.mesh_mode
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $readError = ''
    $result = $null
    if ($run.exit_code -eq 0 -and -not $run.timed_out) {
        try { $result = Read-Z1cCaseResult -RunDir $runDir -ExpectedMeshMode $Case.mesh_mode } catch { $readError = $_.Exception.Message }
    }
    if ($null -eq $result) {
        $status = if ($run.timed_out) { 'timeout' } elseif (-not [string]::IsNullOrWhiteSpace($readError)) { 'read_error' } else { "exit_$($run.exit_code)" }
        $result = New-Z1cEmptyResult -Status $status
    }
    return [pscustomobject]@{ case=$Case; run=$run; result=$result; read_error=$readError; run_dir=$runDir }
}

function Get-Z1cImprovement([double]$Coarse, [double]$Fine) {
    if (-not (Test-FmmCliFinite $Coarse) -or -not (Test-FmmCliFinite $Fine) -or $Fine -le 0.0) { return [double]::NaN }
    return $Coarse / $Fine
}

function Get-Z1cSlope([double]$CoarseError, [double]$FineError, [double]$CoarsePanels, [double]$FinePanels) {
    if ($CoarseError -le 0.0 -or $FineError -le 0.0 -or $FinePanels -le $CoarsePanels) { return [double]::NaN }
    return [math]::Log($FineError / $CoarseError) / [math]::Log($FinePanels / $CoarsePanels)
}

function Get-Z1cMedian([double[]]$Values) {
    $finite = @($Values | Where-Object { Test-FmmCliFinite $_ } | Sort-Object)
    if ($finite.Count -eq 0) { return [double]::NaN }
    $middle = [int][math]::Floor($finite.Count / 2)
    if (($finite.Count % 2) -eq 1) { return [double]$finite[$middle] }
    return 0.5 * ([double]$finite[$middle - 1] + [double]$finite[$middle])
}

function Write-Z1cSvgPlot([object[]]$Rows, [string]$Path) {
    $width = 1120; $height = 760; $left = 80; $top = 70; $panelWidth = 450; $panelHeight = 260
    $colors = @{ icosphere='#0072b2'; latlong='#d55e00' }
    $fixtures = @('born_centered','kirkwood_offcenter'); $salts = @(0.0,0.15)
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add("<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
    $svg.Add("<rect width='100%' height='100%' fill='white'/>")
    $svg.Add("<text x='560' y='30' text-anchor='middle' font-family='sans-serif' font-size='20'>FAGBI Z1C Geometry and Energy Convergence</text>")
    for ($fi=0; $fi -lt $fixtures.Count; ++$fi) {
        for ($si=0; $si -lt $salts.Count; ++$si) {
            $x0=$left+$si*520; $y0=$top+$fi*330
            $svg.Add("<rect x='$x0' y='$y0' width='$panelWidth' height='$panelHeight' fill='none' stroke='#222'/>")
            $svg.Add("<text x='$($x0+$panelWidth/2)' y='$($y0-12)' text-anchor='middle' font-family='sans-serif' font-size='14'>$($fixtures[$fi]), salt=$($salts[$si]) M</text>")
            foreach ($mesh in @('icosphere','latlong')) {
                $series=@($Rows | Where-Object { $_.fixture -eq $fixtures[$fi] -and [math]::Abs([double]$_.salt_M-$salts[$si]) -lt 1.0e-12 -and $_.mesh_mode -eq $mesh -and $_.status -eq 'ok' } | Sort-Object target_panels)
                $energyPoints=New-Object System.Collections.Generic.List[string]
                $geometryPoints=New-Object System.Collections.Generic.List[string]
                foreach($row in $series) {
                    $x=$x0+35+([math]::Log([double]$row.target_panels/80.0)/[math]::Log(64.0))*($panelWidth-70)
                    $energy=[math]::Max(1.0e-6,[double]$row.oracle_relative_error)
                    $geometry=[math]::Max(1.0e-6,[double]$row.planar_area_rel_error)
                    $energyY=$y0+20+((0.0-[math]::Log10($energy))/6.0)*($panelHeight-55)
                    $geometryY=$y0+20+((0.0-[math]::Log10($geometry))/6.0)*($panelHeight-55)
                    $energyPoints.Add("$([math]::Round($x,2)),$([math]::Round($energyY,2))")
                    $geometryPoints.Add("$([math]::Round($x,2)),$([math]::Round($geometryY,2))")
                }
                if($energyPoints.Count -gt 0) {
                    $svg.Add("<polyline points='$($energyPoints -join ' ')' fill='none' stroke='$($colors[$mesh])' stroke-width='2'/>")
                    $svg.Add("<polyline points='$($geometryPoints -join ' ')' fill='none' stroke='$($colors[$mesh])' stroke-width='1.5' stroke-dasharray='5,4'/>")
                }
            }
        }
    }
    $svg.Add("<text x='160' y='735' font-family='sans-serif' font-size='12'>solid: energy error; dashed: planar-area error</text>")
    $svg.Add("<text x='650' y='735' font-family='sans-serif' font-size='12' fill='$($colors.icosphere)'>icosphere</text>")
    $svg.Add("<text x='760' y='735' font-family='sans-serif' font-size='12' fill='$($colors.latlong)'>latlong</text>")
    $svg.Add('</svg>')
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$caseSpecs = @(Get-FagbiZ1cCaseSpecs -Quick:$Quick)
foreach ($case in $caseSpecs) {
    $case | Add-Member -NotePropertyName oracle_energy_kj_mol -NotePropertyValue $(if ($case.fixture -eq 'born_centered') {
        Get-FagbiBornEnergyKJ -Charge 1.0 -RadiusA 3.0 -SaltM $case.salt_M -TemperatureK $TemperatureK -InternalDielectric $InternalDielectric -ExternalDielectric $ExternalDielectric
    } else {
        Get-FagbiKirkwoodEnergyKJ -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0 -SaltM $case.salt_M -TemperatureK $TemperatureK -InternalDielectric $InternalDielectric -ExternalDielectric $ExternalDielectric
    })
}
$records = foreach ($case in $caseSpecs) { Invoke-Z1cCase -Case $case }
$rows = foreach ($record in $records) {
    $case=$record.case; $result=$record.result; $energy=[double]$result.energy_kj_mol; $oracle=[double]$case.oracle_energy_kj_mol
    $oracleError=if((Test-FmmCliFinite $energy)-and(Test-FmmCliFinite $oracle)){[math]::Abs($energy-$oracle)/[math]::Max(1.0e-12,[math]::Abs($oracle))}else{[double]::NaN}
    $status=if($result.status -eq 'ok' -and $result.actual_panels -eq $case.expected_panels -and $result.near_fallbacks -eq 0 -and (Test-FmmCliFinite $oracleError)){'ok'}elseif($result.status -ne 'ok'){$result.status}elseif($result.actual_panels -ne $case.expected_panels){'panel_count_mismatch'}elseif($result.near_fallbacks -ne 0){'near_fallback'}else{'nonfinite'}
    [pscustomobject]@{
        fixture=$case.fixture; salt_M=$case.salt_M; mesh_mode=$case.mesh_mode; target_panels=$case.target_panels
        requested_points=$case.requested_points; expected_panels=$case.expected_panels; actual_panels=$result.actual_panels
        status=$status; read_error=$record.read_error; exit_code=$record.run.exit_code; timed_out=$record.run.timed_out
        energy_kj_mol=$energy; oracle_energy_kj_mol=$oracle; oracle_relative_error=$oracleError
        iterations=$result.iterations; final_residual=$result.final_residual; true_residual=$result.true_residual
        wall_seconds=$record.run.wall_seconds; solver_total_seconds=$result.total_seconds; near_fallbacks=$result.near_fallbacks
        spherical_area=$result.spherical_area; spherical_area_rel_error=$result.spherical_area_rel_error
        planar_area=$result.planar_area; planar_area_rel_error=$result.planar_area_rel_error
        enclosed_volume=$result.enclosed_volume; enclosed_volume_rel_error=$result.enclosed_volume_rel_error
        triangle_area_min=$result.triangle_area_min; triangle_area_mean=$result.triangle_area_mean; triangle_area_max=$result.triangle_area_max
        triangle_area_cv=$result.triangle_area_cv; triangle_aspect_max=$result.triangle_aspect_max
        normal_radial_dot_min=$result.normal_radial_dot_min; normal_radial_dot_mean=$result.normal_radial_dot_mean
        centroid_radial_error_max=$result.centroid_radial_error_max
    }
}

$convergenceRows=New-Object System.Collections.Generic.List[object]
foreach($group in ($rows | Group-Object fixture,salt_M,mesh_mode)) {
    $previous=$null
    foreach($row in @($group.Group | Sort-Object target_panels)) {
        $energyGain=if($null -eq $previous){[double]::NaN}else{Get-Z1cImprovement ([double]$previous.oracle_relative_error) ([double]$row.oracle_relative_error)}
        $planarGain=if($null -eq $previous){[double]::NaN}else{Get-Z1cImprovement ([double]$previous.planar_area_rel_error) ([double]$row.planar_area_rel_error)}
        $volumeGain=if($null -eq $previous){[double]::NaN}else{Get-Z1cImprovement ([double]$previous.enclosed_volume_rel_error) ([double]$row.enclosed_volume_rel_error)}
        $slope=if($null -eq $previous){[double]::NaN}else{Get-Z1cSlope ([double]$previous.oracle_relative_error) ([double]$row.oracle_relative_error) ([double]$previous.actual_panels) ([double]$row.actual_panels)}
        $convergenceRows.Add([pscustomobject]@{
            fixture=$row.fixture; salt_M=$row.salt_M; mesh_mode=$row.mesh_mode; target_panels=$row.target_panels
            actual_panels=$row.actual_panels; oracle_relative_error=$row.oracle_relative_error
            planar_area_rel_error=$row.planar_area_rel_error; enclosed_volume_rel_error=$row.enclosed_volume_rel_error
            energy_improvement=$energyGain; planar_area_improvement=$planarGain; volume_improvement=$volumeGain
            energy_convergence_slope=$slope; status=$row.status
        }) | Out-Null
        $previous=$row
    }
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$convergenceRows | Export-Csv -LiteralPath $ConvergenceCsv -NoTypeInformation
Write-Z1cSvgPlot -Rows $rows -Path $PlotPath

$invalidRows=@($rows | Where-Object { $_.status -ne 'ok' })
$finestTarget=($caseSpecs.target_panels | Measure-Object -Maximum).Maximum
$finestRows=@($rows | Where-Object { $_.target_panels -eq $finestTarget -and $_.status -eq 'ok' })
$bestByPhysics=New-Object System.Collections.Generic.List[object]
foreach($group in ($finestRows | Group-Object fixture,salt_M)) {
    $bestByPhysics.Add((@($group.Group | Sort-Object oracle_relative_error)[0])) | Out-Null
}
$allFinestBelowHalfPercent=$bestByPhysics.Count -eq (@($caseSpecs | Group-Object fixture,salt_M)).Count -and @($bestByPhysics | Where-Object { [double]$_.oracle_relative_error -gt 5.0e-3 }).Count -eq 0

$transitionRows=@($convergenceRows | Where-Object { $_.target_panels -eq $finestTarget -and $_.status -eq 'ok' })
$resolutionLimitedGroups=@($transitionRows | Where-Object { [double]$_.energy_improvement -ge 1.5 })
$plateauGroups=@($transitionRows | Where-Object {
    [double]$_.planar_area_improvement -ge 2.0 -and [double]$_.volume_improvement -ge 2.0 -and
    [double]$_.energy_improvement -lt 1.25 -and [double]$_.oracle_relative_error -gt 1.0e-2
})
$accuracyClass=if($invalidRows.Count -gt 0){'inconclusive'}elseif($allFinestBelowHalfPercent){'resolution_converged'}elseif($transitionRows.Count -gt 0 -and $resolutionLimitedGroups.Count -ge [math]::Ceiling(0.75*$transitionRows.Count)){'resolution_limited'}elseif($transitionRows.Count -gt 0 -and $plateauGroups.Count -ge [math]::Ceiling(0.5*$transitionRows.Count)){'operator_discretization_plateau'}else{'inconclusive'}

$meshStats=@{}
foreach($mesh in @('icosphere','latlong')) {
    $meshRows=@($finestRows | Where-Object { $_.mesh_mode -eq $mesh })
    $meshErrors=@($meshRows | ForEach-Object { $_.oracle_relative_error })
    $meshAspects=@($meshRows | ForEach-Object { $_.triangle_aspect_max })
    $meshCvs=@($meshRows | ForEach-Object { $_.triangle_area_cv })
    $meshStats[$mesh]=[pscustomobject]@{
        error=Get-Z1cMedian $meshErrors
        aspect=Get-Z1cMedian $meshAspects
        cv=Get-Z1cMedian $meshCvs
    }
}
$preferredMesh='none'; $meshComparisonClass='mesh_quality_tie'
if((Test-FmmCliFinite $meshStats.icosphere.error)-and(Test-FmmCliFinite $meshStats.latlong.error)) {
    $winner=if($meshStats.icosphere.error -le $meshStats.latlong.error){'icosphere'}else{'latlong'}
    $loser=if($winner -eq 'icosphere'){'latlong'}else{'icosphere'}
    $ratio=$meshStats[$loser].error/[math]::Max(1.0e-15,$meshStats[$winner].error)
    if($ratio -ge 1.5 -and $meshStats[$winner].aspect -le $meshStats[$loser].aspect -and $meshStats[$winner].cv -le $meshStats[$loser].cv) {
        $preferredMesh=$winner; $meshComparisonClass='mesh_quality_dominant'
    }
}

[pscustomobject]@{
    accuracy_class=$accuracyClass; preferred_mesh=$preferredMesh; mesh_comparison_class=$meshComparisonClass
    quick=$Quick.IsPresent; rows=$rows.Count; invalid_rows=$invalidRows.Count; finest_target_panels=$finestTarget
    resolution_limited_groups=$resolutionLimitedGroups.Count; plateau_groups=$plateauGroups.Count
    matrix=$MatrixCsv; convergence=$ConvergenceCsv; plot=$PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary=@(
    'TITAN FAGBI Scope Z1C Surface Discretization Summary',
    "accuracy_class: $accuracyClass",
    "preferred_mesh: $preferredMesh",
    "mesh_comparison_class: $meshComparisonClass",
    "quick: $($Quick.IsPresent)",
    "rows: $($rows.Count)",
    "invalid_rows: $($invalidRows.Count)",
    'backend: tabi_direct',
    'near_integration: analytic_singular',
    'scientific_note: This gate isolates surface geometry and direct-operator discretization; FMM/M2L is disabled.',
    "matrix: $MatrixCsv",
    "convergence: $ConvergenceCsv",
    "plot: $PlotPath"
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report=New-Object System.Collections.Generic.List[string]
$report.Add('TITAN FMM/M2L Scope Z1C Surface Discretization Report')
$report.Add("Accuracy class: $accuracyClass")
$report.Add("Preferred mesh: $preferredMesh")
$report.Add("Mesh comparison: $meshComparisonClass")
$report.Add("Generated: $(Get-Date -Format s)")
$report.Add('')
$report.Add('Purpose')
$report.Add('Separate surface-resolution and mesh-quality error from direct operator discretization before Scope Z2 near-field integration.')
$report.Add('')
$report.Add('Matrix')
$report.Add('- Fixtures: centered Born and off-center Kirkwood.')
$report.Add('- Backend: TABI Direct; FMM/M2L disabled; analytic_singular near integration.')
$report.Add("- Mesh modes: icosphere and latlong; target panels: $(@($caseSpecs.target_panels | Sort-Object -Unique) -join ',').")
$report.Add("- Salts: $(@($caseSpecs.salt_M | Sort-Object -Unique) -join ',') M.")
$report.Add('')
$report.Add('Measured classification')
$report.Add("- Invalid rows: $($invalidRows.Count).")
$report.Add("- Resolution-limited transitions: $($resolutionLimitedGroups.Count)/$($transitionRows.Count).")
$report.Add("- Operator-plateau transitions: $($plateauGroups.Count)/$($transitionRows.Count).")
$report.Add("- Icosphere median finest energy error: $($meshStats.icosphere.error).")
$report.Add("- Latlong median finest energy error: $($meshStats.latlong.error).")
$report.Add('')
$report.Add('Next decision')
switch($accuracyClass) {
    'resolution_converged' { $report.Add('Proceed to Scope Z2A on matched geometry and integrate the demonstrated near policy into FMM/M2L.') }
    'resolution_limited' { $report.Add('Plan Scope Z1D with a surface-resolution range beyond the current approximately 5120-panel cap.') }
    'operator_discretization_plateau' { $report.Add('Review collocation, panel basis, and diagonal formulation before Scope Z2.') }
    default { $report.Add('Keep Scope Z2 blocked and inspect the geometry/energy convergence matrix before selecting the next scientific change.') }
}
if($meshComparisonClass -eq 'mesh_quality_dominant') { $report.Add("Use $preferredMesh as the controlled mesh for subsequent FAGBI benchmark gates.") }
$report.Add('')
$report.Add("Matrix CSV: $MatrixCsv")
$report.Add("Convergence CSV: $ConvergenceCsv")
$report.Add("SVG: $PlotPath")
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

$expectedRows=if($Quick){8}else{32}
if($rows.Count -ne $expectedRows){throw "Scope Z1C expected $expectedRows rows, got $($rows.Count)"}
if(-not(Test-Path -LiteralPath $PlotPath)){throw "Scope Z1C did not create SVG plot: $PlotPath"}
if($invalidRows.Count -gt 0){Write-Warning "Scope Z1C completed with $($invalidRows.Count) invalid rows."}
Write-Host "Scope Z1C accuracy class: $accuracyClass"
Write-Host "Scope Z1C preferred mesh: $preferredMesh"
Write-Host "Matrix: $MatrixCsv"
