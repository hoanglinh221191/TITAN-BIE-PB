param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_born_ion_matrix'),
    [double[]]$RadiiA = @(1.5, 2.0, 3.0, 5.0),
    [double[]]$Charges = @(1.0, -1.0, 2.0, -2.0),
    [double[]]$SaltM = @(0.0, 0.05, 0.15, 0.5),
    [int[]]$SurfacePointCounts = @(96, 192, 384, 768),
    [string[]]$SphereMeshModes = @('latlong'),
    [double[]]$GridScales = @(0.60, 0.45, 0.35),
    [string[]]$NearIntegrationModes = @('analytic_singular','duffy_adaptive'),
    [switch]$SkipGrid
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Convert-ToDoubleOrNull([object]$Value) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $null }
    return [double]::Parse($text, $InvariantCulture)
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside Born ion matrix root: $full"
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
        'CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

function Write-BornPqr([string]$Path, [double]$Charge, [double]$RadiusA) {
    $chargeText = $Charge.ToString('F4', $InvariantCulture)
    $radiusText = $RadiusA.ToString('F4', $InvariantCulture)
    Set-Content -LiteralPath $Path -Value @(
        ("ATOM      1  ION ION     1       0.000   0.000   0.000{0,8}{1,8}       X" -f $chargeText, $radiusText),
        'END'
    ) -Encoding ASCII
}

function Get-DebyeKappaInvA([double]$Salt, [double]$TempK) {
    if ($Salt -le 0.0) { return 0.0 }
    return 0.3291 * [math]::Sqrt($Salt) * [math]::Sqrt(298.15 / [math]::Max(1.0, $TempK))
}

function Get-BornKirkwoodEnergyKJ([double]$Charge, [double]$RadiusA,
                                  [double]$EpsIn, [double]$EpsOut,
                                  [double]$Salt, [double]$TempK) {
    $kappa = Get-DebyeKappaInvA -Salt $Salt -TempK $TempK
    $outside = 1.0 / ($EpsOut * (1.0 + $kappa * $RadiusA))
    $inside = 1.0 / $EpsIn
    $kcal = 0.5 * 332.063713299 * $Charge * $Charge * ($outside - $inside) / $RadiusA
    return $kcal * 4.184
}

function Get-TabiThreshold([int]$SurfacePoints) {
    if ($SurfacePoints -ge 768) { return 1.0e-2 }
    if ($SurfacePoints -ge 384) { return 1.5e-2 }
    if ($SurfacePoints -ge 192) { return 2.5e-2 }
    return 5.0e-2
}

function Get-LinearSlope([double[]]$X, [double[]]$Y) {
    $pairs = @()
    for ($i = 0; $i -lt $X.Count; ++$i) {
        if ($X[$i] -gt 0.0 -and $Y[$i] -gt 0.0) {
            $pairs += [PSCustomObject]@{ x = [math]::Log($X[$i]); y = [math]::Log($Y[$i]) }
        }
    }
    if ($pairs.Count -lt 2) { return 'NA' }
    $meanX = ($pairs | Measure-Object -Property x -Average).Average
    $meanY = ($pairs | Measure-Object -Property y -Average).Average
    $num = 0.0
    $den = 0.0
    foreach ($p in $pairs) {
        $dx = $p.x - $meanX
        $num += $dx * ($p.y - $meanY)
        $den += $dx * $dx
    }
    if ($den -le 0.0) { return 'NA' }
    return ($num / $den).ToString('G17', $InvariantCulture)
}

function New-ArgFile([string]$RunDir, [string]$Backend, [string]$PqrPath,
                     [double]$Salt, [int]$SurfacePoints, [double]$GridScale,
                     [string]$NearMode, [string]$SphereMeshMode) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $saltText = $Salt.ToString('G17', $InvariantCulture)
    $scaleText = $GridScale.ToString('G17', $InvariantCulture)
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
        '-pbe-probe=0.0',
        "-pbe-salt=$saltText",
        '-pbe-temp=298.15',
        '-pbe-timing=summary',
        '-pbe-map=off',
        '-pbe-margin=2.0',
        "-pbe-scale=$scaleText",
        '-pbe-max-grid=550000',
        '-pbe-max-iter=1200',
        '-pbe-rms=1e-5',
        '-nt=1'
    )
    if ($Backend -eq 'tabi') {
        $argLines += '-pbe-tabi-sum=direct'
        $argLines += '-pbe-tabi-formulation=jcp2013'
        $argLines += '-pbe-tabi-energy-projection=boundary'
        $argLines += "-pbe-tabi-near-integration=$NearMode"
        $argLines += "-pbe-tabi-surface-points=$SurfacePoints"
        $argLines += "-pbe-tabi-sphere-mesh=$SphereMeshMode"
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Invoke-BornCase([string]$Backend, [double]$Charge, [double]$RadiusA,
                         [double]$Salt, [int]$SurfacePoints,
                         [double]$GridScale, [string]$NearMode,
                         [string]$SphereMeshMode) {
    $safeCharge = ($Charge.ToString('G17', $InvariantCulture) -replace '-', 'm' -replace '\.', 'p')
    $safeRadius = ($RadiusA.ToString('G17', $InvariantCulture) -replace '\.', 'p')
    $safeSalt = ($Salt.ToString('G17', $InvariantCulture) -replace '\.', 'p')
    $safeScale = ($GridScale.ToString('G17', $InvariantCulture) -replace '\.', 'p')
    $safeNear = if ([string]::IsNullOrWhiteSpace($NearMode)) { 'grid' } else { $NearMode }
    $safeMesh = if ($Backend -eq 'tabi') { $SphereMeshMode } else { 'grid' }
    $suffix = if ($Backend -eq 'tabi') { "sp${SurfacePoints}_${safeNear}_${safeMesh}" } else { "gs$safeScale" }
    $runDir = Join-Path $WorkRoot ("{0}_q{1}_r{2}_s{3}_{4}" -f $Backend, $safeCharge, $safeRadius, $safeSalt, $suffix)
    New-CleanDir $runDir
    $pqrPath = Join-Path $runDir 'born_ion.pqr'
    Write-BornPqr -Path $pqrPath -Charge $Charge -RadiusA $RadiusA
    $argPath = New-ArgFile -RunDir $runDir -Backend $Backend -PqrPath $pqrPath -Salt $Salt -SurfacePoints $SurfacePoints -GridScale $GridScale -NearMode $NearMode -SphereMeshMode $SphereMeshMode
    $measure = Measure-Command {
        Push-Location $runDir
        try {
            & $TitanExe -a $argPath | Out-Null
            if ($LASTEXITCODE -ne 0) {
                throw "TITAN exited with code $LASTEXITCODE for $Backend q=$Charge r=$RadiusA salt=$Salt near=$NearMode"
            }
        } finally {
            Pop-Location
        }
    }
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    $timing = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_timing.csv'))[0]
    $memoryPath = Join-Path $runDir 'pbe_memory_profile.csv'
    $memory = if (Test-Path -LiteralPath $memoryPath) { @(Import-Csv -LiteralPath $memoryPath)[0] } else { $null }
    $auditPath = Join-Path $runDir 'pbe_tabi_formulation_audit.csv'
    $audit = if (Test-Path -LiteralPath $auditPath) { @(Import-Csv -LiteralPath $auditPath)[0] } else { $null }
    $expected = Get-BornKirkwoodEnergyKJ -Charge $Charge -RadiusA $RadiusA -EpsIn 1.0 -EpsOut 80.0 -Salt $Salt -TempK 298.15
    $actual = Convert-ToDoubleOrNull $energy.E_rxn_kJmol
    $rel = [math]::Abs($actual - $expected) / [math]::Max(1.0e-12, [math]::Abs($expected))
    $gridPointsOrPanels = if ($Backend -eq 'tabi') { $timing.tabi_surface_panels } else { $energy.grid_points }
    return [PSCustomObject]@{
        backend = $Backend
        charge_e = $Charge.ToString('G17', $InvariantCulture)
        radius_A = $RadiusA.ToString('G17', $InvariantCulture)
        salt_M = $Salt.ToString('G17', $InvariantCulture)
        surface_points = if ($Backend -eq 'tabi') { $SurfacePoints } else { 0 }
        sphere_mesh_mode = if ($Backend -eq 'tabi') { $SphereMeshMode } else { 'grid' }
        grid_scale_A = if ($Backend -eq 'grid') { $GridScale.ToString('G17', $InvariantCulture) } else { 'NA' }
        near_integration = if ($null -ne $audit) { $audit.near_integration } else { 'grid' }
        near_kernel_policy = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_kernel_policy')) { $audit.near_kernel_policy } else { if ($Backend -eq 'tabi') { 'unknown' } else { 'grid' } }
        near_k1_accepted = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k1_accepted')) { $audit.near_k1_accepted } else { 'NA' }
        near_k1_rejected = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k1_rejected')) { $audit.near_k1_rejected } else { 'NA' }
        near_k2_accepted = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k2_accepted')) { $audit.near_k2_accepted } else { 'NA' }
        near_k2_rejected = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k2_rejected')) { $audit.near_k2_rejected } else { 'NA' }
        near_k3_accepted = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k3_accepted')) { $audit.near_k3_accepted } else { 'NA' }
        near_k3_rejected = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k3_rejected')) { $audit.near_k3_rejected } else { 'NA' }
        near_k4_accepted = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k4_accepted')) { $audit.near_k4_accepted } else { 'NA' }
        near_k4_rejected = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k4_rejected')) { $audit.near_k4_rejected } else { 'NA' }
        near_one_ring_k1_accepted = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_one_ring_k1_accepted')) { $audit.near_one_ring_k1_accepted } else { 'NA' }
        near_one_ring_k1_rejected = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_one_ring_k1_rejected')) { $audit.near_one_ring_k1_rejected } else { 'NA' }
        near_one_ring_k2_accepted = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_one_ring_k2_accepted')) { $audit.near_one_ring_k2_accepted } else { 'NA' }
        near_one_ring_k2_rejected = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_one_ring_k2_rejected')) { $audit.near_one_ring_k2_rejected } else { 'NA' }
        near_one_ring_k3_accepted = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_one_ring_k3_accepted')) { $audit.near_one_ring_k3_accepted } else { 'NA' }
        near_one_ring_k3_rejected = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_one_ring_k3_rejected')) { $audit.near_one_ring_k3_rejected } else { 'NA' }
        near_one_ring_k4_accepted = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_one_ring_k4_accepted')) { $audit.near_one_ring_k4_accepted } else { 'NA' }
        near_one_ring_k4_rejected = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_one_ring_k4_rejected')) { $audit.near_one_ring_k4_rejected } else { 'NA' }
        near_k1_policy = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k1_policy')) { $audit.near_k1_policy } else { 'NA' }
        near_k2_policy = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k2_policy')) { $audit.near_k2_policy } else { 'NA' }
        near_k3_policy = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k3_policy')) { $audit.near_k3_policy } else { 'NA' }
        near_k4_policy = if ($null -ne $audit -and ($audit.PSObject.Properties.Name -contains 'near_k4_policy')) { $audit.near_k4_policy } else { 'NA' }
        analytic_kJmol = $expected.ToString('G17', $InvariantCulture)
        energy_kJmol = $energy.E_rxn_kJmol
        rel_error_vs_analytic = $rel.ToString('G17', $InvariantCulture)
        born_correctness_status = 'pending'
        max_rel_error_by_density = 'NA'
        convergence_slope = 'NA'
        status = $energy.status
        reaction_energy_source = $energy.reaction_energy_source
        fallback_used = if ($null -ne $audit) { $audit.fallback_used } else { 'false' }
        near_fallbacks = if ($null -ne $audit) { $audit.near_fallbacks } else { 'NA' }
        grid_points = $energy.grid_points
        panels = if ($Backend -eq 'tabi') { $timing.tabi_surface_panels } else { '0' }
        grid_points_or_panels = $gridPointsOrPanels
        total_seconds = if ($timing.PSObject.Properties.Name -contains 'total_seconds') { $timing.total_seconds } else { $measure.TotalSeconds.ToString('G17', $InvariantCulture) }
        solver_seconds = if ($timing.PSObject.Properties.Name -contains 'solver_seconds') { $timing.solver_seconds } else { 'NA' }
        peak_working_set_bytes = if ($null -ne $memory) { $memory.peak_working_set_bytes } else { 'NA' }
        run_dir = $runDir
    }
}

function Set-Note([object]$Row, [string]$Name, [object]$Value) {
    $Row | Add-Member -NotePropertyName $Name -NotePropertyValue $Value -Force
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
New-CleanDir $WorkRoot

$rows = New-Object System.Collections.Generic.List[object]
$TabiRowsByKey = @{}

foreach ($radius in $RadiiA) {
    foreach ($charge in $Charges) {
        foreach ($salt in $SaltM) {
            foreach ($points in $SurfacePointCounts) {
                foreach ($meshMode in $SphereMeshModes) {
                    foreach ($nearMode in $NearIntegrationModes) {
                        if ([string]::IsNullOrWhiteSpace($nearMode)) { continue }
                        $row = Invoke-BornCase -Backend 'tabi' -Charge $charge -RadiusA $radius -Salt $salt -SurfacePoints $points -GridScale $GridScales[0] -NearMode $nearMode -SphereMeshMode $meshMode
                        $rows.Add($row)
                        $key = "$charge|$radius|$salt|$points|$nearMode|$meshMode"
                        $TabiRowsByKey[$key] = $row
                        if ($row.status -ne 'ok' -or $row.reaction_energy_source -eq 'tabi_born_fallback' -or $row.fallback_used -eq 'true') {
                            throw "TABI Born ion matrix used fallback or failed: status=$($row.status) source=$($row.reaction_energy_source) fallback=$($row.fallback_used)"
                        }
                        if ($row.near_integration -ne $nearMode) {
                            throw "TABI Born ion matrix expected near=$nearMode, got $($row.near_integration)"
                        }
                        if ($nearMode -eq 'analytic_singular' -and $row.near_kernel_policy -ne 'self_asi') {
                            throw "TABI analytic_singular did not report self_asi policy: $($row.near_kernel_policy)"
                        }
                        if ($nearMode -eq 'duffy_adaptive' -and $row.near_kernel_policy -ne 'one_ring_cpv') {
                            throw "TABI duffy_adaptive did not report one_ring_cpv policy: $($row.near_kernel_policy)"
                        }
                        if ($nearMode -eq 'duffy_adaptive' -and $salt -gt 0.0) {
                            foreach ($block in @('k1','k2','k3','k4')) {
                                $acceptedName = "near_${block}_accepted"
                                $rejectedName = "near_${block}_rejected"
                                if (([int64]$row.$acceptedName + [int64]$row.$rejectedName) -le 0) {
                                    throw "TABI duffy_adaptive did not record any $block near-kernel decisions."
                                }
                                $oneRingAcceptedName = "near_one_ring_${block}_accepted"
                                $oneRingRejectedName = "near_one_ring_${block}_rejected"
                                if (([int64]$row.$oneRingAcceptedName + [int64]$row.$oneRingRejectedName) -le 0) {
                                    throw "TABI duffy_adaptive did not record any one-ring $block near-kernel decisions."
                                }
                                if ([int64]$row.$rejectedName -lt 0) {
                                    throw "TABI duffy_adaptive reported invalid $block rejection count."
                                }
                            }
                        }
                        if ([int64]$row.near_fallbacks -ne 0) {
                            throw "TABI Born ion matrix used near fallback count=$($row.near_fallbacks)"
                        }
                        $threshold = Get-TabiThreshold $points
                        $rel = Convert-ToDoubleOrNull $row.rel_error_vs_analytic
                        if ($nearMode -eq 'analytic_singular' -and $rel -gt $threshold) {
                            throw "TABI Born ion analytic mismatch q=$charge r=$radius salt=$salt points=$points mesh=$meshMode rel=$($row.rel_error_vs_analytic) threshold=$threshold"
                        }
                    }
                    if ($TabiRowsByKey.ContainsKey("$charge|$radius|$salt|$points|analytic_singular|$meshMode") -and
                        $TabiRowsByKey.ContainsKey("$charge|$radius|$salt|$points|duffy_adaptive|$meshMode")) {
                        $analytic = $TabiRowsByKey["$charge|$radius|$salt|$points|analytic_singular|$meshMode"]
                        $duffy = $TabiRowsByKey["$charge|$radius|$salt|$points|duffy_adaptive|$meshMode"]
                        $analyticRel = Convert-ToDoubleOrNull $analytic.rel_error_vs_analytic
                        $duffyRel = Convert-ToDoubleOrNull $duffy.rel_error_vs_analytic
                        if ($duffyRel -gt ($analyticRel + 1.0e-6)) {
                            throw "Duffy one-ring CPV worsened Born ion q=$charge r=$radius salt=$salt points=$points mesh=$meshMode analytic_rel=$analyticRel duffy_rel=$duffyRel"
                        }
                    }
                }
            }
            if (-not $SkipGrid) {
                foreach ($scale in $GridScales) {
                    $grid = Invoke-BornCase -Backend 'grid' -Charge $charge -RadiusA $radius -Salt $salt -SurfacePoints $SurfacePointCounts[0] -GridScale $scale -NearMode '' -SphereMeshMode 'grid'
                    $rows.Add($grid)
                    if ($grid.status -ne 'ok') {
                        throw "Grid Born ion diagnostic failed: status=$($grid.status)"
                    }
                }
            }
        }
    }
}

$rowsArray = @($rows.ToArray())
foreach ($densityGroup in ($rowsArray | Group-Object backend,near_integration,surface_points,sphere_mesh_mode,grid_scale_A)) {
    $maxRel = ($densityGroup.Group | ForEach-Object { Convert-ToDoubleOrNull $_.rel_error_vs_analytic } |
        Where-Object { $null -ne $_ } | Measure-Object -Maximum).Maximum
    foreach ($row in $densityGroup.Group) {
        if ($null -ne $maxRel) {
            Set-Note $row 'max_rel_error_by_density' ($maxRel.ToString('G17', $InvariantCulture))
        }
    }
}

foreach ($caseGroup in ($rowsArray | Group-Object backend,near_integration,sphere_mesh_mode,charge_e,radius_A,salt_M)) {
    $x = New-Object System.Collections.Generic.List[double]
    $y = New-Object System.Collections.Generic.List[double]
    foreach ($row in $caseGroup.Group) {
        $rel = Convert-ToDoubleOrNull $row.rel_error_vs_analytic
        if ($null -eq $rel -or $rel -le 0.0) { continue }
        if ($row.backend -eq 'tabi') {
            $x.Add([double]$row.surface_points)
        } elseif ($row.grid_scale_A -ne 'NA') {
            $x.Add(1.0 / [double]::Parse($row.grid_scale_A, $InvariantCulture))
        }
        $y.Add($rel)
    }
    $slope = Get-LinearSlope -X $x.ToArray() -Y $y.ToArray()
    foreach ($row in $caseGroup.Group) {
        Set-Note $row 'convergence_slope' $slope
    }
}

foreach ($row in $rowsArray) {
    if ($row.backend -eq 'tabi') {
        $threshold = Get-TabiThreshold ([int]$row.surface_points)
        $rel = Convert-ToDoubleOrNull $row.rel_error_vs_analytic
        $status = if ($row.near_integration -eq 'duffy_adaptive') {
            'duffy_not_worse'
        } elseif ($rel -le $threshold) {
            'pass'
        } else {
            'fail'
        }
        Set-Note $row 'born_correctness_status' $status
    } else {
        Set-Note $row 'born_correctness_status' 'grid_diagnostic'
    }
}

$outCsv = Join-Path $WorkRoot 'pbe_born_ion_matrix.csv'
$rowsArray | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE Born ion analytic matrix PASS: $outCsv"
