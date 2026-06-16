$script:FagbiInvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Write-FagbiBornPqr(
    [string]$Path,
    [double]$Charge,
    [double]$RadiusA
) {
    $line = [string]::Format(
        $script:FagbiInvariantCulture,
        'ATOM      1  ION ION     1       0.000   0.000   0.000{0,8:F4}{1,8:F4}       X',
        $Charge,
        $RadiusA
    )
    Set-Content -LiteralPath $Path -Value @($line, 'END') -Encoding ASCII
}

function Write-FagbiKirkwoodPqr(
    [string]$Path,
    [double]$SphereRadiusA,
    [double]$Charge,
    [double]$OffsetA
) {
    $lines = @(
        [string]::Format(
            $script:FagbiInvariantCulture,
            'ATOM      1  SHL SPH     1       0.000   0.000   0.000{0,8:F4}{1,8:F4}       X',
            0.0,
            $SphereRadiusA
        ),
        [string]::Format(
            $script:FagbiInvariantCulture,
            'ATOM      2  Q1  ION     2    {0,8:F3}{1,8:F3}{2,8:F3}{3,8:F4}{4,8:F4}       X',
            $OffsetA,
            0.0,
            0.0,
            $Charge,
            0.20
        ),
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Get-FagbiIcosphereTargetPoints([int]$Level) {
    if ($Level -lt 0 -or $Level -gt 6) {
        throw "FAGBI icosphere level must be in 0..6, got $Level"
    }
    return [int](20 * [math]::Pow(4, $Level))
}

function Get-FagbiSurfaceResolution(
    [ValidateSet('icosphere','latlong')][string]$MeshMode,
    [ValidateSet(80,320,1280,5120)][int]$TargetPanels
) {
    if ($MeshMode -eq 'icosphere') {
        $requestedPoints = if ($TargetPanels -eq 5120) { 1281 } else { $TargetPanels }
        return [pscustomobject]@{
            mesh_mode = $MeshMode
            target_panels = $TargetPanels
            requested_points = $requestedPoints
            expected_panels = $TargetPanels
        }
    }
    $mapping = @{
        80 = @{ request=16; actual=80 }
        320 = @{ request=81; actual=360 }
        1280 = @{ request=289; actual=1224 }
        5120 = @{ request=1225; actual=5040 }
    }
    $entry = $mapping[$TargetPanels]
    return [pscustomobject]@{
        mesh_mode = $MeshMode
        target_panels = $TargetPanels
        requested_points = [int]$entry.request
        expected_panels = [int]$entry.actual
    }
}

function Get-FagbiZ1cCaseSpecs([switch]$Quick) {
    $salts = if ($Quick) { @(0.15) } else { @(0.0,0.15) }
    $targets = if ($Quick) { @(320,1280) } else { @(80,320,1280,5120) }
    foreach ($fixture in @('born_centered','kirkwood_offcenter')) {
        foreach ($salt in $salts) {
            foreach ($meshMode in @('icosphere','latlong')) {
                foreach ($target in $targets) {
                    $resolution = Get-FagbiSurfaceResolution `
                        -MeshMode $meshMode `
                        -TargetPanels $target
                    [pscustomobject]@{
                        fixture = $fixture
                        salt_M = [double]$salt
                        mesh_mode = $meshMode
                        target_panels = [int]$target
                        requested_points = [int]$resolution.requested_points
                        expected_panels = [int]$resolution.expected_panels
                    }
                }
            }
        }
    }
}

function Get-FagbiHighresSurfaceResolution(
    [ValidateSet(5120,20480,81920)][int]$TargetPanels
) {
    $mapping = @{
        5120 = @{ request=1281; actual=5120 }
        20480 = @{ request=5121; actual=20480 }
        81920 = @{ request=20481; actual=81920 }
    }
    $entry = $mapping[$TargetPanels]
    return [pscustomobject]@{
        mesh_mode = 'icosphere'
        target_panels = $TargetPanels
        requested_points = [int]$entry.request
        expected_panels = [int]$entry.actual
    }
}

function Get-FagbiZ1dCaseSpecs([switch]$Quick) {
    $fixtures = if ($Quick) { @('born_centered') } else { @('born_centered','kirkwood_offcenter') }
    $salts = if ($Quick) { @(0.15) } else { @(0.0,0.15) }
    foreach ($fixture in $fixtures) {
        foreach ($salt in $salts) {
            foreach ($target in @(5120,20480)) {
                $resolution = Get-FagbiHighresSurfaceResolution -TargetPanels $target
                [pscustomobject]@{
                    fixture = $fixture
                    salt_M = [double]$salt
                    mesh_mode = 'icosphere'
                    target_panels = [int]$target
                    requested_points = [int]$resolution.requested_points
                    expected_panels = [int]$resolution.expected_panels
                }
            }
        }
    }
}

function Get-FagbiDebyeKappaInvA([double]$SaltM, [double]$TemperatureK = 298.15) {
    if ($SaltM -le 0.0) { return 0.0 }
    return 0.3291 * [math]::Sqrt($SaltM) *
        [math]::Sqrt(298.15 / [math]::Max(1.0, $TemperatureK))
}

function Get-FagbiBornEnergyKJ(
    [double]$Charge,
    [double]$RadiusA,
    [double]$SaltM,
    [double]$TemperatureK = 298.15,
    [double]$InternalDielectric = 1.0,
    [double]$ExternalDielectric = 80.0
) {
    $kappa = Get-FagbiDebyeKappaInvA -SaltM $SaltM -TemperatureK $TemperatureK
    $outside = 1.0 / ($ExternalDielectric * (1.0 + $kappa * $RadiusA))
    $inside = 1.0 / $InternalDielectric
    return 0.5 * 332.063713299 * $Charge * $Charge *
        ($outside - $inside) / $RadiusA * 4.184
}

function Get-FagbiLegendreP([int]$L, [double]$X) {
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

function Get-FagbiSphericalKNextOverK([int]$L, [double]$X) {
    if ($X -le 1.0e-12) { return [double]::PositiveInfinity }
    $ratio = 1.0 + 1.0 / $X
    if ($L -eq 0) { return $ratio }
    for ($n = 1; $n -le $L; ++$n) {
        $ratio = (1.0 / $ratio) + ((2.0 * $n + 1.0) / $X)
    }
    return $ratio
}

function Get-FagbiReactionCoefficient(
    [int]$L,
    [double]$SphereRadiusA,
    [double]$SaltM,
    [double]$TemperatureK,
    [double]$InternalDielectric,
    [double]$ExternalDielectric
) {
    $kappa = Get-FagbiDebyeKappaInvA -SaltM $SaltM -TemperatureK $TemperatureK
    $epsRatio = $ExternalDielectric / $InternalDielectric
    if ($kappa -le 1.0e-14) {
        $gamma = -$epsRatio * ($L + 1.0)
    } else {
        $x = $kappa * $SphereRadiusA
        $kNextOverK = Get-FagbiSphericalKNextOverK -L $L -X $x
        $gamma = $epsRatio * (-$x * $kNextOverK + $L)
    }
    return (($L + 1.0) + $gamma) / ($L - $gamma)
}

function Get-FagbiKirkwoodEnergyKJ(
    [double]$SphereRadiusA,
    [double]$Charge,
    [double]$OffsetA,
    [double]$SaltM,
    [double]$TemperatureK = 298.15,
    [double]$InternalDielectric = 1.0,
    [double]$ExternalDielectric = 80.0,
    [int]$LMax = 160
) {
    $ri = [math]::Abs($OffsetA)
    $sum = 0.0
    for ($l = 0; $l -le $LMax; ++$l) {
        $coef = Get-FagbiReactionCoefficient `
            -L $l `
            -SphereRadiusA $SphereRadiusA `
            -SaltM $SaltM `
            -TemperatureK $TemperatureK `
            -InternalDielectric $InternalDielectric `
            -ExternalDielectric $ExternalDielectric
        $radial = if ($l -eq 0) {
            1.0 / $SphereRadiusA
        } else {
            [math]::Pow($ri * $ri, $l) /
                [math]::Pow($SphereRadiusA, (2 * $l + 1))
        }
        $sum += $Charge * $Charge * ($coef / $InternalDielectric) *
            $radial * (Get-FagbiLegendreP -L $l -X 1.0)
    }
    return 0.5 * 332.063713299 * $sum * 4.184
}
