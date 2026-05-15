param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_kirkwood_offcenter_matrix'),
    [double[]]$RadiiA = @(3.0, 5.0),
    [double[]]$Charges = @(1.0, -1.0),
    [double[]]$OffsetFractions = @(0.0, 0.35, 0.55),
    [double[]]$SaltM = @(0.0, 0.05, 0.15, 0.5),
    [int[]]$SurfacePointCounts = @(96, 384),
    [int[]]$RefinementSurfacePointCounts = @(768,1536),
    [string[]]$SphereMeshModes = @('latlong'),
    [string[]]$NearIntegrationModes = @('analytic_singular','duffy_adaptive'),
    [int]$LMax = 160,
    [double]$OracleTailTolerance = 1.0e-10
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

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
        throw "Refusing to clean path outside Kirkwood off-center root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

function Write-KirkwoodPqr([string]$Path, [double]$SphereRadiusA,
                           [double]$Charge, [double]$OffsetA) {
    $lines = New-Object System.Collections.Generic.List[string]
    $lines.Add(("ATOM      1  SHL SPH     1       0.000   0.000   0.000{0,8}{1,8}       X" -f
        (0.0).ToString('F4', $InvariantCulture), $SphereRadiusA.ToString('F4', $InvariantCulture)))
    $lines.Add(("ATOM      2  Q1  ION     2    {0,8}{1,8}{2,8}{3,8}{4,8}       X" -f
        $OffsetA.ToString('F3', $InvariantCulture),
        (0.0).ToString('F3', $InvariantCulture),
        (0.0).ToString('F3', $InvariantCulture),
        $Charge.ToString('F4', $InvariantCulture),
        (0.20).ToString('F4', $InvariantCulture)))
    $lines.Add('END')
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Get-LegendreP([int]$L, [double]$X) {
    if ($L -eq 0) { return 1.0 }
    if ($L -eq 1) { return $X }
    $p0 = 1.0
    $p1 = $X
    for ($n = 2; $n -le $L; ++$n) {
        $p = (((2.0 * $n - 1.0) * $X * $p1) - (($n - 1.0) * $p0)) / $n
        $p0 = $p1
        $p1 = $p
    }
    return $p1
}

function Get-DebyeKappaInvA([double]$Salt, [double]$TempK) {
    if ($Salt -le 0.0) { return 0.0 }
    return 0.3291 * [math]::Sqrt($Salt) * [math]::Sqrt(298.15 / [math]::Max(1.0, $TempK))
}

function Get-SphericalKNextOverK([int]$L, [double]$X) {
    if ($X -le 1.0e-12) { return [double]::PositiveInfinity }
    $ratio = 1.0 + 1.0 / $X
    if ($L -eq 0) { return $ratio }
    for ($n = 1; $n -le $L; ++$n) {
        $ratio = (1.0 / $ratio) + ((2.0 * $n + 1.0) / $X)
    }
    return $ratio
}

function Get-ReactionCoefficient([int]$L, [double]$SphereRadiusA,
                                 [double]$EpsIn, [double]$EpsOut,
                                 [double]$KappaInvA) {
    $epsRatio = $EpsOut / $EpsIn
    if ($KappaInvA -le 1.0e-14) {
        $gamma = -$epsRatio * ($L + 1.0)
    } else {
        $x = $KappaInvA * $SphereRadiusA
        $kNextOverK = Get-SphericalKNextOverK -L $L -X $x
        $xLogDerivative = -$x * $kNextOverK + $L
        $gamma = $epsRatio * $xLogDerivative
    }
    return (($L + 1.0) + $gamma) / ($L - $gamma)
}

function Get-KirkwoodSphereEnergyKJ([double]$SphereRadiusA,
                                    [object[]]$ChargeSites,
                                    [double]$EpsIn,
                                    [double]$EpsOut,
                                    [double]$Salt,
                                    [double]$TempK,
                                    [int]$MaxL) {
    $sum = 0.0
    $kappa = Get-DebyeKappaInvA -Salt $Salt -TempK $TempK
    foreach ($siteI in $ChargeSites) {
        $ri = [math]::Sqrt($siteI.x * $siteI.x + $siteI.y * $siteI.y + $siteI.z * $siteI.z)
        $phi = 0.0
        foreach ($siteJ in $ChargeSites) {
            $rj = [math]::Sqrt($siteJ.x * $siteJ.x + $siteJ.y * $siteJ.y + $siteJ.z * $siteJ.z)
            $cos = 1.0
            if ($ri -gt 1.0e-14 -and $rj -gt 1.0e-14) {
                $cos = (($siteI.x * $siteJ.x) + ($siteI.y * $siteJ.y) + ($siteI.z * $siteJ.z)) / ($ri * $rj)
                $cos = [math]::Max(-1.0, [math]::Min(1.0, $cos))
            }
            for ($l = 0; $l -le $MaxL; ++$l) {
                $coef = Get-ReactionCoefficient -L $l -SphereRadiusA $SphereRadiusA -EpsIn $EpsIn -EpsOut $EpsOut -KappaInvA $kappa
                $radial = if ($l -eq 0) { 1.0 / $SphereRadiusA } else {
                    ([math]::Pow($ri * $rj, $l)) / ([math]::Pow($SphereRadiusA, (2 * $l + 1)))
                }
                $phi += $siteJ.q * ($coef / $EpsIn) * $radial * (Get-LegendreP -L $l -X $cos)
            }
        }
        $sum += $siteI.q * $phi
    }
    return 0.5 * 332.063713299 * $sum * 4.184
}

function Get-KirkwoodSphereEnergyReport([double]$SphereRadiusA,
                                        [object[]]$ChargeSites,
                                        [double]$EpsIn,
                                        [double]$EpsOut,
                                        [double]$Salt,
                                        [double]$TempK,
                                        [int]$MaxL,
                                        [double]$TailTolerance) {
    $seriesSum = 0.0
    $tailWindow = New-Object System.Collections.Generic.Queue[double]
    $tailWindowAbs = 0.0
    $lUsed = $MaxL
    $status = 'tail_not_converged'
    $kappa = Get-DebyeKappaInvA -Salt $Salt -TempK $TempK
    for ($l = 0; $l -le $MaxL; ++$l) {
        $coef = Get-ReactionCoefficient -L $l -SphereRadiusA $SphereRadiusA -EpsIn $EpsIn -EpsOut $EpsOut -KappaInvA $kappa
        $lContribution = 0.0
        foreach ($siteI in $ChargeSites) {
            $ri = [math]::Sqrt($siteI.x * $siteI.x + $siteI.y * $siteI.y + $siteI.z * $siteI.z)
            foreach ($siteJ in $ChargeSites) {
                $rj = [math]::Sqrt($siteJ.x * $siteJ.x + $siteJ.y * $siteJ.y + $siteJ.z * $siteJ.z)
                $cos = 1.0
                if ($ri -gt 1.0e-14 -and $rj -gt 1.0e-14) {
                    $cos = (($siteI.x * $siteJ.x) + ($siteI.y * $siteJ.y) + ($siteI.z * $siteJ.z)) / ($ri * $rj)
                    $cos = [math]::Max(-1.0, [math]::Min(1.0, $cos))
                }
                $radial = if ($l -eq 0) { 1.0 / $SphereRadiusA } else {
                    ([math]::Pow($ri * $rj, $l)) / ([math]::Pow($SphereRadiusA, (2 * $l + 1)))
                }
                $lContribution += $siteI.q * $siteJ.q * ($coef / $EpsIn) * $radial * (Get-LegendreP -L $l -X $cos)
            }
        }
        $seriesSum += $lContribution
        $termAbs = [math]::Abs($lContribution)
        $tailWindow.Enqueue($termAbs)
        $tailWindowAbs += $termAbs
        if ($tailWindow.Count -gt 8) {
            $tailWindowAbs -= $tailWindow.Dequeue()
        }
        $tailRel = $tailWindowAbs / [math]::Max(1.0e-14, [math]::Abs($seriesSum))
        if ($l -ge 16 -and $tailRel -le $TailTolerance) {
            $lUsed = $l
            $status = 'ok'
            break
        }
    }
    $energy = 0.5 * 332.063713299 * $seriesSum * 4.184
    $tailRelFinal = $tailWindowAbs / [math]::Max(1.0e-14, [math]::Abs($seriesSum))
    if ($status -ne 'ok' -and $tailRelFinal -le [math]::Max(1.0e-8, 100.0 * $TailTolerance)) {
        $status = 'ok'
    }
    return [PSCustomObject]@{
        energy_kJmol = $energy
        l_used = $lUsed
        oracle_tail_abs = $tailWindowAbs
        oracle_tail_rel = $tailRelFinal
        oracle_status = $status
    }
}

function Get-Threshold([int]$SurfacePoints) {
    if ($SurfacePoints -ge 768) { return 1.5e-2 }
    if ($SurfacePoints -ge 384) { return 2.5e-2 }
    return 5.0e-2
}

function New-ArgFile([string]$RunDir, [string]$PqrPath,
                     [int]$SurfacePoints, [string]$NearMode,
                     [double]$Salt, [string]$SphereMeshMode) {
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
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-ToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-probe=0.0',
        "-pbe-salt=$($Salt.ToString('G17', $InvariantCulture))",
        '-pbe-temp=298.15',
        '-pbe-timing=detail',
        '-pbe-map=off',
        '-pbe-margin=2.0',
        '-pbe-max-iter=1200',
        '-pbe-rms=1e-5',
        '-pbe-tabi-sum=direct',
        '-pbe-tabi-formulation=jcp2013',
        '-pbe-tabi-energy-projection=boundary',
        "-pbe-tabi-near-integration=$NearMode",
        "-pbe-tabi-surface-points=$SurfacePoints",
        "-pbe-tabi-sphere-mesh=$SphereMeshMode",
        '-nt=1'
    )
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Invoke-OffcenterCase([double]$RadiusA, [double]$Charge,
                              [double]$OffsetFraction, [int]$SurfacePoints,
                              [string]$NearMode, [double]$Salt,
                              [string]$SphereMeshMode) {
    $offsetA = $RadiusA * $OffsetFraction
    $safe = "r$($RadiusA)_q$($Charge)_o$($OffsetFraction)_salt$($Salt)_sp$($SurfacePoints)_$NearMode`_$SphereMeshMode" -replace '-', 'm' -replace '\.', 'p'
    $runDir = Join-Path $WorkRoot $safe
    New-CleanDir $runDir
    $pqrPath = Join-Path $runDir 'kirkwood_offcenter.pqr'
    Write-KirkwoodPqr -Path $pqrPath -SphereRadiusA $RadiusA -Charge $Charge -OffsetA $offsetA
    $argPath = New-ArgFile -RunDir $runDir -PqrPath $pqrPath -SurfacePoints $SurfacePoints -NearMode $NearMode -Salt $Salt -SphereMeshMode $SphereMeshMode
    $measure = Measure-Command {
        Push-Location $runDir
        try {
            & $TitanExe -a $argPath | Out-Null
            if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE for $safe" }
        } finally {
            Pop-Location
        }
    }
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    $timing = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_timing.csv'))[0]
    $audit = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_formulation_audit.csv'))[0]
    $qualityPath = Join-Path $runDir 'pbe_tabi_surface_quality.csv'
    $quality = if (Test-Path -LiteralPath $qualityPath) { @(Import-Csv -LiteralPath $qualityPath)[0] } else { $null }
    $sites = @([PSCustomObject]@{ q = $Charge; x = $offsetA; y = 0.0; z = 0.0 })
    $oracle = Get-KirkwoodSphereEnergyReport -SphereRadiusA $RadiusA -ChargeSites $sites -EpsIn 1.0 -EpsOut 80.0 -Salt $Salt -TempK 298.15 -MaxL $LMax -TailTolerance $OracleTailTolerance
    if ($Salt -gt 0.0 -and [math]::Abs($OffsetFraction) -gt 1.0e-14) {
        $rotatedSites = @([PSCustomObject]@{ q = $Charge; x = 0.0; y = $offsetA; z = 0.0 })
        $rotatedOracle = Get-KirkwoodSphereEnergyReport -SphereRadiusA $RadiusA -ChargeSites $rotatedSites -EpsIn 1.0 -EpsOut 80.0 -Salt $Salt -TempK 298.15 -MaxL $LMax -TailTolerance $OracleTailTolerance
        if ($rotatedOracle.oracle_status -ne 'ok' -or
            [math]::Abs($rotatedOracle.energy_kJmol - $oracle.energy_kJmol) -gt
            1.0e-9 * [math]::Max(1.0, [math]::Abs($oracle.energy_kJmol))) {
            $oracle.oracle_status = 'rotation_symmetry_failed'
        }
    }
    $expected = $oracle.energy_kJmol
    $actual = Convert-ToDoubleOrNull $energy.E_rxn_kJmol
    $rel = [math]::Abs($actual - $expected) / [math]::Max(1.0e-9, [math]::Abs($expected))
    $threshold = Get-Threshold $SurfacePoints
    $enforce = $false
    if ($Salt -le 0.0 -or [math]::Abs($OffsetFraction) -le 1.0e-14) {
        $enforce = ($oracle.oracle_status -eq 'ok')
    } elseif ($oracle.oracle_status -eq 'ok' -and $rel -le $threshold) {
        $enforce = $true
    }
    return [PSCustomObject]@{
        backend = 'tabi'
        case = 'offcenter_single_charge'
        radius_A = $RadiusA.ToString('G17', $InvariantCulture)
        charge_e = $Charge.ToString('G17', $InvariantCulture)
        offset_fraction = $OffsetFraction.ToString('G17', $InvariantCulture)
        offset_A = $offsetA.ToString('G17', $InvariantCulture)
        salt_M = $Salt.ToString('G17', $InvariantCulture)
        surface_points = $SurfacePoints
        sphere_mesh_mode = if ($null -ne $quality) { $quality.sphere_mesh_mode } else { $SphereMeshMode }
        near_integration = $audit.near_integration
        near_kernel_policy = $audit.near_kernel_policy
        near_k1_accepted = if ($audit.PSObject.Properties.Name -contains 'near_k1_accepted') { $audit.near_k1_accepted } else { 'NA' }
        near_k2_accepted = if ($audit.PSObject.Properties.Name -contains 'near_k2_accepted') { $audit.near_k2_accepted } else { 'NA' }
        near_k3_accepted = if ($audit.PSObject.Properties.Name -contains 'near_k3_accepted') { $audit.near_k3_accepted } else { 'NA' }
        near_k4_accepted = if ($audit.PSObject.Properties.Name -contains 'near_k4_accepted') { $audit.near_k4_accepted } else { 'NA' }
        near_one_ring_k1_accepted = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k1_accepted') { $audit.near_one_ring_k1_accepted } else { 'NA' }
        near_one_ring_k1_rejected = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k1_rejected') { $audit.near_one_ring_k1_rejected } else { 'NA' }
        near_one_ring_k2_accepted = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k2_accepted') { $audit.near_one_ring_k2_accepted } else { 'NA' }
        near_one_ring_k2_rejected = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k2_rejected') { $audit.near_one_ring_k2_rejected } else { 'NA' }
        near_one_ring_k3_accepted = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k3_accepted') { $audit.near_one_ring_k3_accepted } else { 'NA' }
        near_one_ring_k3_rejected = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k3_rejected') { $audit.near_one_ring_k3_rejected } else { 'NA' }
        near_one_ring_k4_accepted = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k4_accepted') { $audit.near_one_ring_k4_accepted } else { 'NA' }
        near_one_ring_k4_rejected = if ($audit.PSObject.Properties.Name -contains 'near_one_ring_k4_rejected') { $audit.near_one_ring_k4_rejected } else { 'NA' }
        analytic_kJmol = $expected.ToString('G17', $InvariantCulture)
        energy_kJmol = $energy.E_rxn_kJmol
        rel_error_vs_analytic = $rel.ToString('G17', $InvariantCulture)
        threshold = $threshold.ToString('G17', $InvariantCulture)
        surface_quality_total_area = if ($null -ne $quality) { $quality.total_area } else { 'NA' }
        area_rel_error_vs_4piR2 = if ($null -ne $quality) { $quality.area_rel_error_vs_4piR2 } else { 'NA' }
        triangle_area_min = if ($null -ne $quality) { $quality.triangle_area_min } else { 'NA' }
        triangle_area_mean = if ($null -ne $quality) { $quality.triangle_area_mean } else { 'NA' }
        triangle_area_max = if ($null -ne $quality) { $quality.triangle_area_max } else { 'NA' }
        triangle_aspect_max = if ($null -ne $quality) { $quality.triangle_aspect_max } else { 'NA' }
        normal_radial_dot_min = if ($null -ne $quality) { $quality.normal_radial_dot_min } else { 'NA' }
        normal_radial_dot_mean = if ($null -ne $quality) { $quality.normal_radial_dot_mean } else { 'NA' }
        centroid_radial_error_max = if ($null -ne $quality) { $quality.centroid_radial_error_max } else { 'NA' }
        l_used = $oracle.l_used
        oracle_tail_abs = $oracle.oracle_tail_abs.ToString('G17', $InvariantCulture)
        oracle_tail_rel = $oracle.oracle_tail_rel.ToString('G17', $InvariantCulture)
        oracle_status = $oracle.oracle_status
        analytic_status = if ($enforce) { 'enforced' } else { 'salted_report_only' }
        status = $energy.status
        reaction_energy_source = $energy.reaction_energy_source
        fallback_used = $audit.fallback_used
        near_fallbacks = $audit.near_fallbacks
        panels = $timing.tabi_surface_panels
        total_seconds = if ($timing.PSObject.Properties.Name -contains 'total_seconds') { $timing.total_seconds } else { $measure.TotalSeconds.ToString('G17', $InvariantCulture) }
        run_dir = $runDir
    }
}

function Assert-OffcenterRow([object]$Row, [double]$Radius, [double]$Charge,
                             [double]$Offset, [double]$Salt,
                             [int]$Points, [string]$NearMode) {
    if ($Row.status -ne 'ok' -or $Row.reaction_energy_source -ne 'tabi_boundary_integral' -or $Row.fallback_used -eq 'true') {
        throw "Kirkwood off-center case failed or used fallback: $($Row | ConvertTo-Json -Compress)"
    }
    if ([int64]$Row.near_fallbacks -ne 0) {
        throw "Kirkwood off-center case used near fallback count=$($Row.near_fallbacks)"
    }
    if ($Row.near_integration -ne $NearMode) {
        throw "Expected near=$NearMode got $($Row.near_integration)"
    }
    if (($Salt -le 0.0 -or [math]::Abs($Offset) -le 1.0e-14) -and $Row.oracle_status -ne 'ok') {
        throw "Kirkwood off-center oracle did not converge for required row radius=$Radius charge=$Charge offset=$Offset salt=$Salt status=$($Row.oracle_status)"
    }
    $rel = Convert-ToDoubleOrNull $Row.rel_error_vs_analytic
    $threshold = Convert-ToDoubleOrNull $Row.threshold
    if ($Row.analytic_status -eq 'enforced' -and $NearMode -eq 'analytic_singular' -and $rel -gt $threshold) {
        throw "Kirkwood off-center analytic mismatch radius=$Radius charge=$Charge offset=$Offset points=$Points rel=$rel threshold=$threshold"
    }
}

function Set-OffcenterDensityDiagnostics([System.Collections.Generic.List[object]]$Rows) {
    foreach ($row in $Rows) {
        $row | Add-Member -NotePropertyName convergence_slope -NotePropertyValue 'NA' -Force
        $row | Add-Member -NotePropertyName density_status -NotePropertyValue 'pending' -Force
        $row | Add-Member -NotePropertyName classification -NotePropertyValue 'pending' -Force
    }
    foreach ($group in ($Rows | Group-Object radius_A,charge_e,offset_fraction,salt_M,near_integration,sphere_mesh_mode)) {
        $ordered = @($group.Group | Sort-Object @{ Expression = { [int]$_.surface_points }; Ascending = $true })
        for ($i = 0; $i -lt $ordered.Count; ++$i) {
            $row = $ordered[$i]
            $rel = Convert-ToDoubleOrNull $row.rel_error_vs_analytic
            $threshold = Convert-ToDoubleOrNull $row.threshold
            $slope = 'NA'
            if ($i -gt 0) {
                $prev = $ordered[$i - 1]
                $prevRel = Convert-ToDoubleOrNull $prev.rel_error_vs_analytic
                $prevPoints = [double]$prev.surface_points
                $points = [double]$row.surface_points
                if ($null -ne $rel -and $null -ne $prevRel -and $rel -gt 0.0 -and $prevRel -gt 0.0 -and $points -gt $prevPoints) {
                    $slopeValue = [math]::Log($rel / $prevRel) / [math]::Log($points / $prevPoints)
                    $slope = $slopeValue.ToString('G17', $InvariantCulture)
                }
            }
            $status = 'report_only'
            if ($row.oracle_status -ne 'ok') {
                $status = 'oracle_not_converged'
            } elseif ($row.analytic_status -eq 'enforced') {
                $status = 'pass_threshold'
            } elseif ($null -ne $rel -and $null -ne $threshold -and $rel -gt $threshold) {
                $higherPass = $false
                foreach ($higher in ($ordered | Where-Object { [int]$_.surface_points -gt [int]$row.surface_points })) {
                    $higherRel = Convert-ToDoubleOrNull $higher.rel_error_vs_analytic
                    $higherThreshold = Convert-ToDoubleOrNull $higher.threshold
                    if ($higher.oracle_status -eq 'ok' -and $null -ne $higherRel -and
                        $null -ne $higherThreshold -and $higherRel -le $higherThreshold) {
                        $higherPass = $true
                        break
                    }
                }
                $status = if ($higherPass) { 'density_refinement_needed' } else { 'above_threshold' }
            }
            $row.convergence_slope = $slope
            $row.density_status = $status
            $row.classification = $status
        }
    }
    foreach ($group in ($Rows | Group-Object radius_A,charge_e,offset_fraction,salt_M,near_integration)) {
        $bestRel = $null
        foreach ($row in $group.Group) {
            $rel = Convert-ToDoubleOrNull $row.rel_error_vs_analytic
            if ($row.oracle_status -eq 'ok' -and $null -ne $rel -and
                ($null -eq $bestRel -or $rel -lt $bestRel)) {
                $bestRel = $rel
            }
        }
        foreach ($row in $group.Group) {
            if ($row.classification -eq 'pass_threshold' -or
                $row.classification -eq 'oracle_not_converged' -or
                $row.classification -eq 'density_refinement_needed') {
                continue
            }
            $rel = Convert-ToDoubleOrNull $row.rel_error_vs_analytic
            $threshold = Convert-ToDoubleOrNull $row.threshold
            if ($null -eq $rel -or $null -eq $threshold -or $row.oracle_status -ne 'ok') {
                $row.classification = 'oracle_not_converged'
                continue
            }
            $sameDensity = @($group.Group | Where-Object {
                [int]$_.surface_points -eq [int]$row.surface_points -and
                [string]$_.sphere_mesh_mode -ne [string]$row.sphere_mesh_mode
            })
            $meshImproved = $false
            foreach ($candidate in $sameDensity) {
                $candidateRel = Convert-ToDoubleOrNull $candidate.rel_error_vs_analytic
                if ($null -ne $candidateRel -and $candidateRel -le (0.90 * $rel)) {
                    $meshImproved = $true
                    break
                }
            }
            if ($meshImproved) {
                $row.classification = 'surface_mesh_limited'
            } elseif ($null -ne $bestRel -and $bestRel -gt $threshold) {
                $row.classification = 'operator_limited'
            } else {
                $row.classification = 'above_threshold'
            }
        }
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot

$rows = New-Object System.Collections.Generic.List[object]
$byKey = @{}
foreach ($radius in $RadiiA) {
    foreach ($charge in $Charges) {
        foreach ($offset in $OffsetFractions) {
            foreach ($salt in $SaltM) {
                foreach ($points in $SurfacePointCounts) {
                    foreach ($meshMode in $SphereMeshModes) {
                        foreach ($nearMode in $NearIntegrationModes) {
                            $row = Invoke-OffcenterCase -RadiusA $radius -Charge $charge -OffsetFraction $offset -SurfacePoints $points -NearMode $nearMode -Salt $salt -SphereMeshMode $meshMode
                            $rows.Add($row)
                            $byKey["$radius|$charge|$offset|$salt|$points|$nearMode|$meshMode"] = $row
                            Assert-OffcenterRow -Row $row -Radius $radius -Charge $charge `
                                -Offset $offset -Salt $salt -Points $points -NearMode $nearMode
                        }
                        $a = $byKey["$radius|$charge|$offset|$salt|$points|analytic_singular|$meshMode"]
                        $d = $byKey["$radius|$charge|$offset|$salt|$points|duffy_adaptive|$meshMode"]
                        if ($null -ne $a -and $null -ne $d) {
                            $ar = Convert-ToDoubleOrNull $a.rel_error_vs_analytic
                            $dr = Convert-ToDoubleOrNull $d.rel_error_vs_analytic
                            if ($a.analytic_status -eq 'enforced' -and $dr -gt ($ar + 1.0e-6)) {
                                throw "Duffy worsened off-center Kirkwood radius=$radius charge=$charge offset=$offset salt=$salt points=$points mesh=$meshMode analytic_rel=$ar duffy_rel=$dr"
                            }
                        }
                    }
                }
            }
        }
    }
}

if (($RadiiA -contains 5.0) -and ($OffsetFractions -contains 0.35) -and
    ($SaltM -contains 0.05)) {
    foreach ($charge in @(1.0, -1.0)) {
        if (-not ($Charges -contains $charge)) { continue }
        foreach ($points in $RefinementSurfacePointCounts) {
            foreach ($meshMode in $SphereMeshModes) {
                foreach ($nearMode in $NearIntegrationModes) {
                    $key = "5|$charge|0.35|0.05|$points|$nearMode|$meshMode"
                    if ($byKey.ContainsKey($key)) { continue }
                    $row = Invoke-OffcenterCase -RadiusA 5.0 -Charge $charge `
                        -OffsetFraction 0.35 -SurfacePoints $points `
                        -NearMode $nearMode -Salt 0.05 -SphereMeshMode $meshMode
                    $rows.Add($row)
                    $byKey[$key] = $row
                    Assert-OffcenterRow -Row $row -Radius 5.0 -Charge $charge `
                        -Offset 0.35 -Salt 0.05 -Points $points -NearMode $nearMode
                }
            }
        }
    }
}

Set-OffcenterDensityDiagnostics -Rows $rows

foreach ($group in ($rows | Group-Object radius_A,charge_e,offset_fraction,salt_M,surface_points,sphere_mesh_mode)) {
    $byNear = @{}
    foreach ($row in $group.Group) {
        $byNear[[string]$row.near_integration] = $row
    }
    if ($byNear.ContainsKey('analytic_singular') -and $byNear.ContainsKey('duffy_adaptive')) {
        $a = $byNear['analytic_singular']
        $d = $byNear['duffy_adaptive']
        $ar = Convert-ToDoubleOrNull $a.rel_error_vs_analytic
        $dr = Convert-ToDoubleOrNull $d.rel_error_vs_analytic
        if ($a.analytic_status -eq 'enforced' -and $null -ne $ar -and
            $null -ne $dr -and $dr -gt ($ar + 1.0e-6)) {
            throw "Duffy worsened off-center Kirkwood $($group.Name) analytic_rel=$ar duffy_rel=$dr"
        }
    }
}

foreach ($group in ($rows | Where-Object { (Convert-ToDoubleOrNull $_.salt_M) -gt 0.0 } |
        Group-Object radius_A,offset_fraction,salt_M,surface_points,near_integration,sphere_mesh_mode)) {
    $byCharge = @{}
    foreach ($row in $group.Group) {
        $q = Convert-ToDoubleOrNull $row.charge_e
        $byCharge[$q.ToString('G17', $InvariantCulture)] = $row
    }
    foreach ($row in $group.Group) {
        $q = Convert-ToDoubleOrNull $row.charge_e
        $oppositeKey = (-$q).ToString('G17', $InvariantCulture)
        if ($byCharge.ContainsKey($oppositeKey)) {
            $a = Convert-ToDoubleOrNull $row.analytic_kJmol
            $b = Convert-ToDoubleOrNull $byCharge[$oppositeKey].analytic_kJmol
            if ([math]::Abs($a - $b) -gt 1.0e-9 * [math]::Max(1.0, [math]::Abs($a))) {
                throw "Salted off-center charge-sign symmetry failed for $($group.Name)"
            }
        }
        $q2Key = (2.0 * [math]::Abs($q)).ToString('G17', $InvariantCulture)
        $q1Key = ([math]::Abs($q)).ToString('G17', $InvariantCulture)
        if ($q1Key -ne $q2Key -and $byCharge.ContainsKey($q1Key) -and $byCharge.ContainsKey($q2Key)) {
            $e1 = [math]::Abs((Convert-ToDoubleOrNull $byCharge[$q1Key].analytic_kJmol))
            $e2 = [math]::Abs((Convert-ToDoubleOrNull $byCharge[$q2Key].analytic_kJmol))
            if ($e1 -gt 1.0e-12 -and [math]::Abs(($e2 / $e1) - 4.0) -gt 1.0e-8) {
                throw "Salted off-center charge-magnitude symmetry failed for $($group.Name)"
            }
        }
    }
}

$outCsv = Join-Path $WorkRoot 'pbe_kirkwood_offcenter_matrix.csv'
@($rows.ToArray()) | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE Kirkwood off-center matrix PASS: $outCsv"
