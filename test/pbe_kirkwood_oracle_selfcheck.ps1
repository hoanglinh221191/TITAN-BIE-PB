param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_kirkwood_oracle_selfcheck'),
    [double[]]$RadiiA = @(1.5, 2.0, 3.0, 5.0),
    [double[]]$Charges = @(1.0, -1.0, 2.0, -2.0),
    [double[]]$SaltM = @(0.0, 0.05, 0.15, 0.5),
    [double]$TempK = 298.15,
    [double]$Tolerance = 1.0e-12
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside Kirkwood oracle root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
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

function Get-SeriesCenteredEnergyKJ([double]$SphereRadiusA, [double]$Charge,
                                    [double]$EpsIn, [double]$EpsOut,
                                    [double]$Salt, [double]$TempK) {
    $kappa = Get-DebyeKappaInvA -Salt $Salt -TempK $TempK
    $coef0 = Get-ReactionCoefficient -L 0 -SphereRadiusA $SphereRadiusA `
        -EpsIn $EpsIn -EpsOut $EpsOut -KappaInvA $kappa
    $phi = $Charge * ($coef0 / $EpsIn) / $SphereRadiusA
    return 0.5 * 332.063713299 * $Charge * $phi * 4.184
}

function Get-ClosedCenteredEnergyKJ([double]$SphereRadiusA, [double]$Charge,
                                    [double]$EpsIn, [double]$EpsOut,
                                    [double]$Salt, [double]$TempK) {
    $epsRatio = $EpsOut / $EpsIn
    $kappa = Get-DebyeKappaInvA -Salt $Salt -TempK $TempK
    if ($kappa -le 1.0e-14) {
        $coef0 = (1.0 - $epsRatio) / $epsRatio
    } else {
        $x = $kappa * $SphereRadiusA
        $gamma = -$epsRatio * (1.0 + $x)
        $coef0 = (1.0 + $gamma) / (-$gamma)
    }
    $phi = $Charge * ($coef0 / $EpsIn) / $SphereRadiusA
    return 0.5 * 332.063713299 * $Charge * $phi * 4.184
}

New-CleanDir $WorkRoot
$rows = New-Object System.Collections.Generic.List[object]
foreach ($radius in $RadiiA) {
    foreach ($charge in $Charges) {
        foreach ($salt in $SaltM) {
            $series = Get-SeriesCenteredEnergyKJ -SphereRadiusA $radius -Charge $charge `
                -EpsIn 1.0 -EpsOut 80.0 -Salt $salt -TempK $TempK
            $closed = Get-ClosedCenteredEnergyKJ -SphereRadiusA $radius -Charge $charge `
                -EpsIn 1.0 -EpsOut 80.0 -Salt $salt -TempK $TempK
            $rel = [math]::Abs($series - $closed) / [math]::Max(1.0e-12, [math]::Abs($closed))
            $status = if ($rel -le $Tolerance) { 'ok' } else { 'fail' }
            $rows.Add([PSCustomObject]@{
                radius_A = $radius.ToString('G17', $InvariantCulture)
                charge_e = $charge.ToString('G17', $InvariantCulture)
                salt_M = $salt.ToString('G17', $InvariantCulture)
                kappa_inv_A = (Get-DebyeKappaInvA -Salt $salt -TempK $TempK).ToString('G17', $InvariantCulture)
                series_kJmol = $series.ToString('G17', $InvariantCulture)
                closed_kJmol = $closed.ToString('G17', $InvariantCulture)
                rel_error = $rel.ToString('G17', $InvariantCulture)
                tolerance = $Tolerance.ToString('G17', $InvariantCulture)
                status = $status
            })
            if ($status -ne 'ok') {
                throw "Centered salted Kirkwood oracle selfcheck failed radius=$radius charge=$charge salt=$salt rel=$rel"
            }
        }
    }
}

$outCsv = Join-Path $WorkRoot 'pbe_kirkwood_oracle_selfcheck.csv'
@($rows.ToArray()) | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE Kirkwood oracle selfcheck PASS: $outCsv"
