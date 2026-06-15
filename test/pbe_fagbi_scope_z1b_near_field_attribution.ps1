param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z1b'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 1
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z1b_matrix.csv'
$ConvergenceCsv = Join-Path $WorkRoot 'fagbi_scope_z1b_convergence.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z1b_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z1b_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z1b_near_field_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z1B_NEAR_FIELD_ATTRIBUTION_REPORT.txt'
$TemperatureK = 298.15
$InternalDielectric = 1.0
$ExternalDielectric = 80.0
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope Z1B work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) { Remove-Item -LiteralPath $workRootFull -Recurse -Force }
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function F([double]$Value) { return $Value.ToString('G17', $InvariantCulture) }

function Get-Z1bProp([object]$Object, [string]$Name, [object]$Default = '') {
    if ($null -eq $Object) { return $Default }
    if ($Object.PSObject.Properties.Name -contains $Name) {
        $value = $Object.$Name
        if ($null -ne $value -and -not [string]::IsNullOrWhiteSpace([string]$value)) { return $value }
    }
    return $Default
}

function Read-Z1bFirstCsv([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    $rows = @(Import-Csv -LiteralPath $Path)
    if ($rows.Count -eq 0) { return $null }
    return $rows[0]
}

function Write-Z1bArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [int]$SurfacePoints,
    [double]$SaltM,
    [string]$NearMode
) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdbPath
    $libPath = Convert-FmmCliToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
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
        "-pbe-tabi-near-integration=$NearMode",
        "-pbe-tabi-surface-points=$SurfacePoints",
        '-pbe-tabi-sphere-mesh=icosphere',
        '-pbe-tabi-precond=diagonal',
        "-nt=$Threads"
    )
    if ($NearMode -eq 'duffy_adaptive') {
        $argLines += '-pbe-tabi-quad-buffer=adaptive'
    } else {
        $argLines += '-pbe-tabi-quad-buffer=off'
    }
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

function New-Z1bEmptyResult([string]$Status) {
    return [pscustomobject]@{
        status=$Status; energy_kj_mol=[double]::NaN; actual_panels=0
        iterations=0; final_residual=[double]::NaN; true_residual=[double]::NaN
        total_seconds=[double]::NaN; near_integration=''; quad_buffer=''
        near_self_pairs=0; near_one_ring_pairs=0; near_duffy_pairs=0
        near_centroid_pairs=0; near_fallbacks=0; near_integration_seconds=[double]::NaN
        duffy_subdivisions_max=0; duffy_subdivisions_mean=[double]::NaN
        near_kernel_policy=''; near_k1_accepted=0; near_k1_rejected=0
        near_k2_accepted=0; near_k2_rejected=0; near_k3_accepted=0
        near_k3_rejected=0; near_k4_accepted=0; near_k4_rejected=0
        near_one_ring_k1_accepted=0; near_one_ring_k1_rejected=0
        near_one_ring_k2_accepted=0; near_one_ring_k2_rejected=0
        near_one_ring_k3_accepted=0; near_one_ring_k3_rejected=0
        near_one_ring_k4_accepted=0; near_one_ring_k4_rejected=0
    }
}

function Read-Z1bCaseResult([string]$RunDir, [string]$ExpectedNearMode) {
    $energy = Read-Z1bFirstCsv (Join-Path $RunDir 'pbe_energy.csv')
    $profile = Read-Z1bFirstCsv (Join-Path $RunDir 'pbe_tabi_profile.csv')
    if ($null -eq $energy) { throw "Missing pbe_energy.csv in $RunDir" }
    if ($null -eq $profile) { throw "Missing pbe_tabi_profile.csv in $RunDir" }
    $energyStatus = [string](Get-Z1bProp $energy 'status' 'unknown')
    $profileStatus = [string](Get-Z1bProp $profile 'status' 'unknown')
    $nearMode = [string](Get-Z1bProp $profile 'near_integration' '')
    if ($nearMode -ne $ExpectedNearMode) {
        throw "Expected near integration $ExpectedNearMode, got '$nearMode'"
    }
    $status = if ($energyStatus -eq 'ok' -and $profileStatus -eq 'ok') { 'ok' } else { "energy=$energyStatus;profile=$profileStatus" }
    return [pscustomobject]@{
        status=$status
        energy_kj_mol=Convert-FmmCliDouble (Get-Z1bProp $energy 'E_rxn_kJmol' 'NA')
        actual_panels=Convert-FmmCliInt64 (Get-Z1bProp $profile 'panels' 0)
        iterations=Convert-FmmCliInt64 (Get-Z1bProp $profile 'iterations' 0)
        final_residual=Convert-FmmCliDouble (Get-Z1bProp $profile 'final_residual' 'NA')
        true_residual=Convert-FmmCliDouble (Get-Z1bProp $profile 'true_residual_rel' 'NA')
        total_seconds=Convert-FmmCliDouble (Get-Z1bProp $profile 'total_seconds' 'NA')
        near_integration=$nearMode
        quad_buffer=[string](Get-Z1bProp $profile 'quad_buffer' '')
        near_self_pairs=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_self_pairs' 0)
        near_one_ring_pairs=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_one_ring_pairs' 0)
        near_duffy_pairs=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_duffy_pairs' 0)
        near_centroid_pairs=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_centroid_pairs' 0)
        near_fallbacks=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_fallbacks' 0)
        near_integration_seconds=Convert-FmmCliDouble (Get-Z1bProp $profile 'near_integration_seconds' 'NA')
        duffy_subdivisions_max=Convert-FmmCliInt64 (Get-Z1bProp $profile 'duffy_subdivisions_max' 0)
        duffy_subdivisions_mean=Convert-FmmCliDouble (Get-Z1bProp $profile 'duffy_subdivisions_mean' 'NA')
        near_kernel_policy=[string](Get-Z1bProp $profile 'near_kernel_policy' '')
        near_k1_accepted=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_k1_accepted' 0)
        near_k1_rejected=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_k1_rejected' 0)
        near_k2_accepted=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_k2_accepted' 0)
        near_k2_rejected=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_k2_rejected' 0)
        near_k3_accepted=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_k3_accepted' 0)
        near_k3_rejected=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_k3_rejected' 0)
        near_k4_accepted=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_k4_accepted' 0)
        near_k4_rejected=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_k4_rejected' 0)
        near_one_ring_k1_accepted=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_one_ring_k1_accepted' 0)
        near_one_ring_k1_rejected=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_one_ring_k1_rejected' 0)
        near_one_ring_k2_accepted=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_one_ring_k2_accepted' 0)
        near_one_ring_k2_rejected=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_one_ring_k2_rejected' 0)
        near_one_ring_k3_accepted=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_one_ring_k3_accepted' 0)
        near_one_ring_k3_rejected=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_one_ring_k3_rejected' 0)
        near_one_ring_k4_accepted=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_one_ring_k4_accepted' 0)
        near_one_ring_k4_rejected=Convert-FmmCliInt64 (Get-Z1bProp $profile 'near_one_ring_k4_rejected' 0)
    }
}

function Invoke-Z1bCase([object]$Case) {
    $safeSalt = (F $Case.salt_M) -replace '-', 'm' -replace '\.', 'p'
    $caseName = "{0}_salt{1}_p{2}_{3}" -f $Case.fixture,$safeSalt,$Case.requested_panels,$Case.near_mode
    $runDir = Join-Path $workRootFull (Get-FmmCliSafeName $caseName)
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $pqrPath = Join-Path $runDir 'input.pqr'
    if ($Case.fixture -eq 'born_centered') {
        Write-FagbiBornPqr -Path $pqrPath -Charge 1.0 -RadiusA 3.0
    } else {
        Write-FagbiKirkwoodPqr -Path $pqrPath -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0
    }
    $argPath = Write-Z1bArgFile -RunDir $runDir -PqrPath $pqrPath -SurfacePoints $Case.requested_panels -SaltM $Case.salt_M -NearMode $Case.near_mode
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $readError = ''
    $result = $null
    if ($run.exit_code -eq 0 -and -not $run.timed_out) {
        try { $result = Read-Z1bCaseResult -RunDir $runDir -ExpectedNearMode $Case.near_mode } catch { $readError = $_.Exception.Message }
    }
    if ($null -eq $result) {
        $status = if ($run.timed_out) { 'timeout' } elseif (-not [string]::IsNullOrWhiteSpace($readError)) { 'read_error' } else { "exit_$($run.exit_code)" }
        $result = New-Z1bEmptyResult -Status $status
    }
    return [pscustomobject]@{ case=$Case; run=$run; result=$result; read_error=$readError; run_dir=$runDir }
}

function Get-Z1bMeshSlope([object]$Previous, [object]$Current) {
    if ($null -eq $Previous -or $null -eq $Current) { return [double]::NaN }
    $e0 = [double]$Previous.oracle_relative_error
    $e1 = [double]$Current.oracle_relative_error
    $p0 = [double]$Previous.actual_panels
    $p1 = [double]$Current.actual_panels
    if (-not (Test-FmmCliFinite $e0) -or -not (Test-FmmCliFinite $e1) -or $e0 -le 0.0 -or $e1 -le 0.0 -or $p1 -le $p0) {
        return [double]::NaN
    }
    return [math]::Log($e1 / $e0) / [math]::Log($p1 / $p0)
}

function Write-Z1bSvgPlot([object[]]$Rows, [string]$Path) {
    $width = 1100
    $height = 760
    $left = 80
    $top = 70
    $panelWidth = 450
    $panelHeight = 260
    $colors = @{ centroid='#666666'; analytic_singular='#0072b2'; duffy_adaptive='#d55e00' }
    $fixtures = @('born_centered','kirkwood_offcenter')
    $salts = @(0.0,0.15)
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add("<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
    $svg.Add("<rect width='100%' height='100%' fill='white'/>")
    $svg.Add("<text x='550' y='30' text-anchor='middle' font-family='sans-serif' font-size='20'>FAGBI Z1B Near-Field Attribution</text>")
    for ($fi = 0; $fi -lt $fixtures.Count; ++$fi) {
        for ($si = 0; $si -lt $salts.Count; ++$si) {
            $x0 = $left + $si * 520
            $y0 = $top + $fi * 330
            $svg.Add("<rect x='$x0' y='$y0' width='$panelWidth' height='$panelHeight' fill='none' stroke='#222'/>")
            $svg.Add("<text x='$($x0 + $panelWidth / 2)' y='$($y0 - 12)' text-anchor='middle' font-family='sans-serif' font-size='14'>$($fixtures[$fi]), salt=$($salts[$si]) M</text>")
            foreach ($tick in @(20,80,320,1280)) {
                $x = $x0 + 35 + ([math]::Log($tick / 20.0) / [math]::Log(64.0)) * ($panelWidth - 70)
                $svg.Add("<line x1='$x' y1='$($y0 + $panelHeight - 25)' x2='$x' y2='$($y0 + $panelHeight - 20)' stroke='#222'/>")
                $svg.Add("<text x='$x' y='$($y0 + $panelHeight - 5)' text-anchor='middle' font-family='sans-serif' font-size='10'>$tick</text>")
            }
            foreach ($mode in @('centroid','analytic_singular','duffy_adaptive')) {
                $series = @($Rows | Where-Object { $_.fixture -eq $fixtures[$fi] -and [math]::Abs([double]$_.salt_M - $salts[$si]) -lt 1.0e-12 -and $_.near_mode -eq $mode -and $_.status -eq 'ok' } | Sort-Object actual_panels)
                $points = New-Object System.Collections.Generic.List[string]
                foreach ($row in $series) {
                    $err = [math]::Max(1.0e-6, [double]$row.oracle_relative_error)
                    $x = $x0 + 35 + ([math]::Log([double]$row.actual_panels / 20.0) / [math]::Log(64.0)) * ($panelWidth - 70)
                    $y = $y0 + 20 + ((0.0 - [math]::Log10($err)) / 6.0) * ($panelHeight - 55)
                    $points.Add("$([math]::Round($x,2)),$([math]::Round($y,2))")
                }
                if ($points.Count -gt 0) {
                    $svg.Add("<polyline points='$($points -join ' ')' fill='none' stroke='$($colors[$mode])' stroke-width='2'/>")
                }
            }
        }
    }
    $legendX = 750
    $legendY = 735
    $offset = 0
    foreach ($mode in @('centroid','analytic_singular','duffy_adaptive')) {
        $svg.Add("<line x1='$($legendX + $offset)' y1='$legendY' x2='$($legendX + $offset + 25)' y2='$legendY' stroke='$($colors[$mode])' stroke-width='3'/>")
        $svg.Add("<text x='$($legendX + $offset + 30)' y='$($legendY + 4)' font-family='sans-serif' font-size='11'>$mode</text>")
        $offset += 120
    }
    $svg.Add('</svg>')
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$requestedPanels = if ($Quick) { @(20,80) } else { @(20,80,320,1280) }
$salts = if ($Quick) { @(0.15) } else { @(0.0,0.15) }
$nearModes = @('centroid','analytic_singular','duffy_adaptive')
$cases = New-Object System.Collections.Generic.List[object]
foreach ($fixture in @('born_centered','kirkwood_offcenter')) {
    foreach ($salt in $salts) {
        foreach ($panels in $requestedPanels) {
            foreach ($nearMode in $nearModes) {
                $oracle = if ($fixture -eq 'born_centered') {
                    Get-FagbiBornEnergyKJ -Charge 1.0 -RadiusA 3.0 -SaltM $salt -TemperatureK $TemperatureK -InternalDielectric $InternalDielectric -ExternalDielectric $ExternalDielectric
                } else {
                    Get-FagbiKirkwoodEnergyKJ -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0 -SaltM $salt -TemperatureK $TemperatureK -InternalDielectric $InternalDielectric -ExternalDielectric $ExternalDielectric
                }
                $cases.Add([pscustomobject]@{ fixture=$fixture; salt_M=$salt; requested_panels=$panels; near_mode=$nearMode; oracle_energy_kj_mol=$oracle }) | Out-Null
            }
        }
    }
}

$records = foreach ($case in $cases) { Invoke-Z1bCase -Case $case }
$baseRows = foreach ($record in $records) {
    $case = $record.case
    $result = $record.result
    $energy = [double]$result.energy_kj_mol
    $oracle = [double]$case.oracle_energy_kj_mol
    $oracleError = if ((Test-FmmCliFinite $energy) -and (Test-FmmCliFinite $oracle)) { [math]::Abs($energy - $oracle) / [math]::Max(1.0e-12, [math]::Abs($oracle)) } else { [double]::NaN }
    [pscustomobject]@{
        fixture=$case.fixture; salt_M=$case.salt_M; near_mode=$case.near_mode
        requested_panels=$case.requested_panels; actual_panels=$result.actual_panels
        status=$result.status; read_error=$record.read_error; exit_code=$record.run.exit_code
        timed_out=$record.run.timed_out; energy_kj_mol=$energy; oracle_energy_kj_mol=$oracle
        oracle_relative_error=$oracleError; improvement_vs_centroid=[double]::NaN
        mesh_convergence_slope=[double]::NaN; iterations=$result.iterations
        final_residual=$result.final_residual; true_residual=$result.true_residual
        wall_seconds=$record.run.wall_seconds; solver_total_seconds=$result.total_seconds
        near_integration_seconds=$result.near_integration_seconds
        near_self_pairs=$result.near_self_pairs; near_one_ring_pairs=$result.near_one_ring_pairs
        near_duffy_pairs=$result.near_duffy_pairs; near_centroid_pairs=$result.near_centroid_pairs
        near_fallbacks=$result.near_fallbacks; duffy_subdivisions_max=$result.duffy_subdivisions_max
        duffy_subdivisions_mean=$result.duffy_subdivisions_mean; near_kernel_policy=$result.near_kernel_policy
        near_k1_accepted=$result.near_k1_accepted; near_k1_rejected=$result.near_k1_rejected
        near_k2_accepted=$result.near_k2_accepted; near_k2_rejected=$result.near_k2_rejected
        near_k3_accepted=$result.near_k3_accepted; near_k3_rejected=$result.near_k3_rejected
        near_k4_accepted=$result.near_k4_accepted; near_k4_rejected=$result.near_k4_rejected
        near_one_ring_k1_accepted=$result.near_one_ring_k1_accepted; near_one_ring_k1_rejected=$result.near_one_ring_k1_rejected
        near_one_ring_k2_accepted=$result.near_one_ring_k2_accepted; near_one_ring_k2_rejected=$result.near_one_ring_k2_rejected
        near_one_ring_k3_accepted=$result.near_one_ring_k3_accepted; near_one_ring_k3_rejected=$result.near_one_ring_k3_rejected
        near_one_ring_k4_accepted=$result.near_one_ring_k4_accepted; near_one_ring_k4_rejected=$result.near_one_ring_k4_rejected
        attribution_status='pending'
    }
}

$rows = foreach ($row in $baseRows) {
    $centroid = @($baseRows | Where-Object { $_.fixture -eq $row.fixture -and [double]$_.salt_M -eq [double]$row.salt_M -and $_.actual_panels -eq $row.actual_panels -and $_.near_mode -eq 'centroid' }) | Select-Object -First 1
    if ($null -ne $centroid -and (Test-FmmCliFinite ([double]$centroid.oracle_relative_error)) -and (Test-FmmCliFinite ([double]$row.oracle_relative_error)) -and [double]$row.oracle_relative_error -gt 0.0) {
        $row.improvement_vs_centroid = [double]$centroid.oracle_relative_error / [double]$row.oracle_relative_error
    }
    if ($row.status -ne 'ok') { $row.attribution_status = 'run_failed' }
    elseif (-not (Test-FmmCliFinite ([double]$row.oracle_relative_error)) -or $row.near_fallbacks -gt 0) { $row.attribution_status = 'invalid_or_fallback' }
    else { $row.attribution_status = 'ok' }
    $row
}

$convergenceRows = New-Object System.Collections.Generic.List[object]
foreach ($group in ($rows | Group-Object fixture,salt_M,near_mode)) {
    $previous = $null
    foreach ($row in @($group.Group | Sort-Object actual_panels)) {
        $slope = Get-Z1bMeshSlope -Previous $previous -Current $row
        $row.mesh_convergence_slope = $slope
        $convergenceRows.Add([pscustomobject]@{
            fixture=$row.fixture; salt_M=$row.salt_M; near_mode=$row.near_mode
            requested_panels=$row.requested_panels; actual_panels=$row.actual_panels
            oracle_relative_error=$row.oracle_relative_error
            improvement_vs_centroid=$row.improvement_vs_centroid
            convergence_slope=$slope; status=$row.attribution_status
        }) | Out-Null
        $previous = $row
    }
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$convergenceRows | Export-Csv -LiteralPath $ConvergenceCsv -NoTypeInformation
Write-Z1bSvgPlot -Rows $rows -Path $PlotPath

$requiredRows = @($rows | Where-Object { [math]::Abs([double]$_.salt_M - 0.15) -lt 1.0e-12 })
$invalidRows = @($rows | Where-Object { $_.attribution_status -ne 'ok' })
$requiredGroups = @($requiredRows | Group-Object fixture,actual_panels)
$strongGroups = New-Object System.Collections.Generic.List[object]
$lowGainGroups = New-Object System.Collections.Generic.List[object]
$regressionGroups = New-Object System.Collections.Generic.List[object]
foreach ($group in $requiredGroups) {
    $centroid = @($group.Group | Where-Object { $_.near_mode -eq 'centroid' }) | Select-Object -First 1
    $nonCentroid = @($group.Group | Where-Object { $_.near_mode -ne 'centroid' })
    if ($null -eq $centroid -or $nonCentroid.Count -eq 0) { continue }
    $best = @($nonCentroid | Sort-Object oracle_relative_error) | Select-Object -First 1
    $gain = if ([double]$best.oracle_relative_error -gt 0.0) { [double]$centroid.oracle_relative_error / [double]$best.oracle_relative_error } else { [double]::PositiveInfinity }
    if ($gain -ge 2.0) { $strongGroups.Add($group) | Out-Null }
    if ($gain -lt 1.5) { $lowGainGroups.Add($group) | Out-Null }
    if (@($nonCentroid | Where-Object { [double]$_.oracle_relative_error -gt 1.25 * [double]$centroid.oracle_relative_error }).Count -gt 0) {
        $regressionGroups.Add($group) | Out-Null
    }
}

$strongByFixture = @{}
foreach ($fixture in @('born_centered','kirkwood_offcenter')) {
    $strongByFixture[$fixture] = @($strongGroups | Where-Object { $_.Group[0].fixture -eq $fixture }).Count
}
$strongBoth = $strongByFixture['born_centered'] -ge 2 -and $strongByFixture['kirkwood_offcenter'] -ge 2
$finestBestErrors = @()
foreach ($fixture in @('born_centered','kirkwood_offcenter')) {
    $fixtureRows = @($requiredRows | Where-Object { $_.fixture -eq $fixture -and $_.near_mode -ne 'centroid' -and $_.attribution_status -eq 'ok' })
    if ($fixtureRows.Count -gt 0) {
        $maxPanels = ($fixtureRows | Measure-Object actual_panels -Maximum).Maximum
        $best = @($fixtureRows | Where-Object { $_.actual_panels -eq $maxPanels } | Sort-Object oracle_relative_error) | Select-Object -First 1
        $finestBestErrors += [double]$best.oracle_relative_error
    }
}
$finestAboveOnePercent = @($finestBestErrors | Where-Object { $_ -gt 1.0e-2 }).Count -gt 0
$meshDrops = 0
foreach ($fixture in @('born_centered','kirkwood_offcenter')) {
    $centroidRows = @($requiredRows | Where-Object { $_.fixture -eq $fixture -and $_.near_mode -eq 'centroid' -and $_.attribution_status -eq 'ok' } | Sort-Object actual_panels)
    if ($centroidRows.Count -ge 2 -and [double]$centroidRows[-1].oracle_relative_error -le 0.75 * [double]$centroidRows[0].oracle_relative_error) { ++$meshDrops }
}
$majorityLowGain = $requiredGroups.Count -gt 0 -and $lowGainGroups.Count -gt ($requiredGroups.Count / 2.0)

$decision = if ($invalidRows.Count -gt 0 -or $regressionGroups.Count -ge 2) {
    'z1b_near_method_regression'
} elseif ($strongBoth -and $finestAboveOnePercent) {
    'z1b_mixed_quadrature_mesh'
} elseif ($strongBoth) {
    'z1b_quadrature_dominant'
} elseif ($majorityLowGain -and $meshDrops -eq 2) {
    'z1b_mesh_dominant'
} else {
    'z1b_attribution_inconclusive'
}

$winnerRows = @($requiredRows | Where-Object { $_.near_mode -ne 'centroid' -and $_.attribution_status -eq 'ok' } | Sort-Object oracle_relative_error)
$winner = if ($winnerRows.Count -gt 0) { $winnerRows[0].near_mode } else { 'none' }
[pscustomobject]@{
    scope_z1b_decision=$decision; quick=$Quick.IsPresent; rows=$rows.Count
    invalid_rows=$invalidRows.Count; strong_groups=$strongGroups.Count
    low_gain_groups=$lowGainGroups.Count; regression_groups=$regressionGroups.Count
    mesh_drop_fixtures=$meshDrops; best_non_centroid_mode=$winner
    matrix=$MatrixCsv; convergence=$ConvergenceCsv; plot=$PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'TITAN FAGBI Scope Z1B Near-Field Attribution Summary',
    "scope_z1b_decision: $decision",
    "quick: $($Quick.IsPresent)",
    "rows: $($rows.Count)",
    "salts_M: $($salts -join ',')",
    "requested_panel_targets: $($requestedPanels -join ',')",
    'near_modes: centroid,analytic_singular,duffy_adaptive',
    "best_non_centroid_mode: $winner",
    "invalid_rows: $($invalidRows.Count)",
    "strong_groups: $($strongGroups.Count)",
    "regression_groups: $($regressionGroups.Count)",
    'attribution_note: Direct TABI isolates shared near quadrature and surface discretization without FMM approximation.',
    "matrix: $MatrixCsv",
    "convergence: $ConvergenceCsv",
    "plot: $PlotPath"
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = New-Object System.Collections.Generic.List[string]
$report.Add('TITAN FMM/M2L Scope Z1B Near-Field Attribution Report')
$report.Add("Decision: $decision")
$report.Add("Generated: $(Get-Date -Format s)")
$report.Add('')
$report.Add('Purpose')
$report.Add('Separate shared TABI near-field quadrature error from surface-discretization error before wiring Duffy/cubature into FMM/M2L.')
$report.Add('')
$report.Add('Matrix')
$report.Add("- Fixtures: centered Born and off-center Kirkwood.")
$report.Add("- Salts: $($salts -join ',') M.")
$report.Add("- Requested icosphere panel targets: $($requestedPanels -join ',').")
$report.Add('- Near modes: centroid, analytic_singular, duffy_adaptive.')
$report.Add('- Backend: TABI direct only; FMM/M2L is disabled.')
$report.Add('')
$report.Add('Measured attribution')
$report.Add("- Required 0.15 M groups with >=2x best non-centroid improvement: $($strongGroups.Count)/$($requiredGroups.Count).")
$report.Add("- Required groups with a >25% non-centroid regression: $($regressionGroups.Count).")
$report.Add("- Fixtures with clear centroid mesh refinement: $meshDrops/2.")
$report.Add("- Best non-centroid mode by minimum required-row error: $winner.")
$report.Add('')
$report.Add('Interpretation')
switch ($decision) {
    'z1b_quadrature_dominant' { $report.Add('Near quadrature dominates and the winning non-centroid method reaches <=1% at the finest tested mesh. Proceed to FMM near integration Scope Z2A.') }
    'z1b_mixed_quadrature_mesh' { $report.Add('Near quadrature materially improves accuracy, but the finest mesh remains above 1%. Scope Z2A should integrate the winning method while retaining a mesh-resolution study.') }
    'z1b_mesh_dominant' { $report.Add('Mesh refinement explains most of the oracle error. Prioritize surface-resolution validation before modifying the FMM near path.') }
    'z1b_near_method_regression' { $report.Add('At least one non-centroid path regressed or produced invalid/fallback rows. Repair direct near integration before FMM integration.') }
    default { $report.Add('The current matrix does not isolate one dominant source. Keep Z2 blocked and expand diagnostics deliberately.') }
}
$report.Add('')
$report.Add('Scientific boundary')
$report.Add('This report attributes the direct TABI discretization error. It does not claim that FMM/M2L near integration has been fixed or validated.')
$report.Add('')
$report.Add("Matrix CSV: $MatrixCsv")
$report.Add("Convergence CSV: $ConvergenceCsv")
$report.Add("SVG: $PlotPath")
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

if ($Quick -and $rows.Count -ne 12) { throw "Scope Z1B quick expected 12 rows, got $($rows.Count)" }
if (-not $Quick -and $rows.Count -ne 48) { throw "Scope Z1B full expected 48 rows, got $($rows.Count)" }
if (@($rows | Where-Object { $_.actual_panels -notin $requestedPanels -and $_.status -eq 'ok' }).Count -gt 0) {
    throw 'Scope Z1B observed an unexpected actual icosphere panel count.'
}
if (-not (Test-Path -LiteralPath $PlotPath)) { throw "Scope Z1B did not create SVG plot: $PlotPath" }

Write-Host "Scope Z1B decision: $decision"
Write-Host "Matrix: $MatrixCsv"
