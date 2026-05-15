param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_kirkwood_multicharge_matrix'),
    [int[]]$SurfacePointCounts = @(96, 384),
    [string[]]$SphereMeshModes = @('latlong'),
    [double[]]$SaltM = @(0.0, 0.05, 0.15, 0.5),
    [string[]]$NearIntegrationModes = @('analytic_singular','duffy_adaptive'),
    [int]$LMax = 180,
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
        throw "Refusing to clean path outside Kirkwood multi-charge root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   70.000   70.000   70.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

function New-ChargeSite([double]$Q, [double]$X, [double]$Y, [double]$Z) {
    return [PSCustomObject]@{ q = $Q; x = $X; y = $Y; z = $Z }
}

function Get-ChargeSets() {
    return @(
        [PSCustomObject]@{
            name = 'same_sign_pair'
            radius_A = 5.0
            sites = @(
                (New-ChargeSite 1.0 1.0 0.0 0.0),
                (New-ChargeSite 1.0 0.0 1.1 0.0)
            )
        },
        [PSCustomObject]@{
            name = 'opposite_pair'
            radius_A = 5.0
            sites = @(
                (New-ChargeSite 1.0 1.0 0.0 0.0),
                (New-ChargeSite -1.0 -1.2 0.0 0.0)
            )
        },
        [PSCustomObject]@{
            name = 'asymmetric_triple'
            radius_A = 5.0
            sites = @(
                (New-ChargeSite 1.0 1.0 0.0 0.0),
                (New-ChargeSite -0.5 0.0 1.2 0.0),
                (New-ChargeSite 0.75 0.0 0.0 -1.4)
            )
        }
    )
}

function Write-KirkwoodPqr([string]$Path, [double]$SphereRadiusA,
                           [object[]]$Sites) {
    $lines = New-Object System.Collections.Generic.List[string]
    $lines.Add(("ATOM      1  SHL SPH     1       0.000   0.000   0.000{0,8}{1,8}       X" -f
        (0.0).ToString('F4', $InvariantCulture), $SphereRadiusA.ToString('F4', $InvariantCulture)))
    $serial = 2
    foreach ($site in $Sites) {
        $lines.Add(("ATOM  {0,5}  Q{1,-2} ION  {0,4}    {2,8}{3,8}{4,8}{5,8}{6,8}       X" -f
            $serial,
            ($serial - 1),
            ([double]$site.x).ToString('F3', $InvariantCulture),
            ([double]$site.y).ToString('F3', $InvariantCulture),
            ([double]$site.z).ToString('F3', $InvariantCulture),
            ([double]$site.q).ToString('F4', $InvariantCulture),
            (0.20).ToString('F4', $InvariantCulture)))
        ++$serial
    }
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

function Convert-Sites([object[]]$Sites, [scriptblock]$Mapper) {
    $out = New-Object System.Collections.Generic.List[object]
    foreach ($site in $Sites) {
        $mapped = & $Mapper $site
        $out.Add($mapped)
    }
    return $out.ToArray()
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
        '-pbe-timing=summary',
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

function Invoke-MultichargeCase([object]$Set, [int]$SurfacePoints,
                                [string]$NearMode, [double]$Salt,
                                [string]$SphereMeshMode) {
    $safe = "$($Set.name)_salt$($Salt)_sp$($SurfacePoints)_$NearMode`_$SphereMeshMode" -replace '-', 'm' -replace '\.', 'p'
    $runDir = Join-Path $WorkRoot $safe
    New-CleanDir $runDir
    $pqrPath = Join-Path $runDir 'kirkwood_multicharge.pqr'
    Write-KirkwoodPqr -Path $pqrPath -SphereRadiusA ([double]$Set.radius_A) -Sites $Set.sites
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
    $oracle = Get-KirkwoodSphereEnergyReport -SphereRadiusA ([double]$Set.radius_A) -ChargeSites $Set.sites -EpsIn 1.0 -EpsOut 80.0 -Salt $Salt -TempK 298.15 -MaxL $LMax -TailTolerance $OracleTailTolerance
    if ($Salt -gt 0.0) {
        $rotatedSites = Convert-Sites $Set.sites { param($site) New-ChargeSite ([double]$site.q) ([double]$site.y) ([double]$site.z) ([double]$site.x) }
        $flippedSites = Convert-Sites $Set.sites { param($site) New-ChargeSite (-[double]$site.q) ([double]$site.x) ([double]$site.y) ([double]$site.z) }
        $scaledSites = Convert-Sites $Set.sites { param($site) New-ChargeSite (2.0 * [double]$site.q) ([double]$site.x) ([double]$site.y) ([double]$site.z) }
        $rotatedOracle = Get-KirkwoodSphereEnergyReport -SphereRadiusA ([double]$Set.radius_A) -ChargeSites $rotatedSites -EpsIn 1.0 -EpsOut 80.0 -Salt $Salt -TempK 298.15 -MaxL $LMax -TailTolerance $OracleTailTolerance
        $flippedOracle = Get-KirkwoodSphereEnergyReport -SphereRadiusA ([double]$Set.radius_A) -ChargeSites $flippedSites -EpsIn 1.0 -EpsOut 80.0 -Salt $Salt -TempK 298.15 -MaxL $LMax -TailTolerance $OracleTailTolerance
        $scaledOracle = Get-KirkwoodSphereEnergyReport -SphereRadiusA ([double]$Set.radius_A) -ChargeSites $scaledSites -EpsIn 1.0 -EpsOut 80.0 -Salt $Salt -TempK 298.15 -MaxL $LMax -TailTolerance $OracleTailTolerance
        $scale = [math]::Max(1.0, [math]::Abs($oracle.energy_kJmol))
        if ($rotatedOracle.oracle_status -ne 'ok' -or
            [math]::Abs($rotatedOracle.energy_kJmol - $oracle.energy_kJmol) -gt 1.0e-9 * $scale) {
            $oracle.oracle_status = 'rotation_symmetry_failed'
        }
        if ($flippedOracle.oracle_status -ne 'ok' -or
            [math]::Abs($flippedOracle.energy_kJmol - $oracle.energy_kJmol) -gt 1.0e-9 * $scale) {
            $oracle.oracle_status = 'charge_sign_symmetry_failed'
        }
        if ($scaledOracle.oracle_status -ne 'ok' -or
            [math]::Abs(($scaledOracle.energy_kJmol / $oracle.energy_kJmol) - 4.0) -gt 1.0e-8) {
            $oracle.oracle_status = 'charge_magnitude_symmetry_failed'
        }
    }
    $expected = $oracle.energy_kJmol
    $actual = Convert-ToDoubleOrNull $energy.E_rxn_kJmol
    $rel = [math]::Abs($actual - $expected) / [math]::Max(1.0, [math]::Abs($expected))
    $threshold = Get-Threshold $SurfacePoints
    $enforce = $false
    if ($Salt -le 0.0) {
        $enforce = ($oracle.oracle_status -eq 'ok')
    } elseif ($oracle.oracle_status -eq 'ok' -and $rel -le $threshold) {
        $enforce = $true
    }
    return [PSCustomObject]@{
        backend = 'tabi'
        case = $Set.name
        radius_A = ([double]$Set.radius_A).ToString('G17', $InvariantCulture)
        charge_count = @($Set.sites).Count
        salt_M = $Salt.ToString('G17', $InvariantCulture)
        surface_points = $SurfacePoints
        sphere_mesh_mode = $SphereMeshMode
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

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot

$rows = New-Object System.Collections.Generic.List[object]
$byKey = @{}
foreach ($set in (Get-ChargeSets)) {
    foreach ($salt in $SaltM) {
        foreach ($points in $SurfacePointCounts) {
            foreach ($meshMode in $SphereMeshModes) {
                foreach ($nearMode in $NearIntegrationModes) {
                    $row = Invoke-MultichargeCase -Set $set -SurfacePoints $points -NearMode $nearMode -Salt $salt -SphereMeshMode $meshMode
                    $rows.Add($row)
                    $byKey["$($set.name)|$salt|$points|$nearMode|$meshMode"] = $row
                    if ($row.status -ne 'ok' -or $row.reaction_energy_source -ne 'tabi_boundary_integral' -or $row.fallback_used -eq 'true') {
                        throw "Kirkwood multi-charge case failed or used fallback: $($row | ConvertTo-Json -Compress)"
                    }
                    if ([int64]$row.near_fallbacks -ne 0) {
                        throw "Kirkwood multi-charge case used near fallback count=$($row.near_fallbacks)"
                    }
                    if ($row.near_integration -ne $nearMode) {
                        throw "Expected near=$nearMode got $($row.near_integration)"
                    }
                    if ($salt -le 0.0 -and $row.oracle_status -ne 'ok') {
                        throw "Kirkwood multi-charge oracle did not converge for required row case=$($set.name) salt=$salt status=$($row.oracle_status)"
                    }
                    $rel = Convert-ToDoubleOrNull $row.rel_error_vs_analytic
                    $threshold = Convert-ToDoubleOrNull $row.threshold
                    if ($row.analytic_status -eq 'enforced' -and $nearMode -eq 'analytic_singular' -and $rel -gt $threshold) {
                        throw "Kirkwood multi-charge analytic mismatch case=$($set.name) salt=$salt points=$points mesh=$meshMode rel=$rel threshold=$threshold"
                    }
                }
                $a = $byKey["$($set.name)|$salt|$points|analytic_singular|$meshMode"]
                $d = $byKey["$($set.name)|$salt|$points|duffy_adaptive|$meshMode"]
                if ($null -ne $a -and $null -ne $d) {
                    $ar = Convert-ToDoubleOrNull $a.rel_error_vs_analytic
                    $dr = Convert-ToDoubleOrNull $d.rel_error_vs_analytic
                    if ($a.analytic_status -eq 'enforced' -and $dr -gt ($ar + 1.0e-6)) {
                        throw "Duffy worsened multi-charge Kirkwood case=$($set.name) salt=$salt points=$points mesh=$meshMode analytic_rel=$ar duffy_rel=$dr"
                    }
                }
            }
        }
    }
}

$outCsv = Join-Path $WorkRoot 'pbe_kirkwood_multicharge_matrix.csv'
@($rows.ToArray()) | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE Kirkwood multi-charge matrix PASS: $outCsv"
