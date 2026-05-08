param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_kirkwood_analytic'),
    [double[]]$SaltM = @(0.0, 0.15)
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
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI Kirkwood root: $full"
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
        'CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

function Get-DebyeKappaInvA([double]$Salt, [double]$TempK) {
    if ($Salt -le 0.0) {
        return 0.0
    }
    return 0.3291 * [math]::Sqrt($Salt) * [math]::Sqrt(298.15 / [math]::Max(1.0, $TempK))
}

function Get-KirkwoodEnergyKJ([double]$Charge, [double]$RadiusA, [double]$EpsIn,
                              [double]$EpsOut, [double]$Salt, [double]$TempK) {
    $kappa = Get-DebyeKappaInvA -Salt $Salt -TempK $TempK
    $outside = 1.0 / ($EpsOut * (1.0 + $kappa * $RadiusA))
    $inside = 1.0 / $EpsIn
    $kcal = 0.5 * 332.063713299 * $Charge * $Charge * ($outside - $inside) / $RadiusA
    return $kcal * 4.184
}

function New-ArgFile([string]$RunDir, [double]$Salt, [string]$Projection) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $pqrPath = Join-Path $InputRoot 'sphere.pqr'
    if (-not (Test-Path -LiteralPath $pqrPath)) {
        throw "Missing sphere PQR fixture: $pqrPath"
    }
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $saltText = $Salt.ToString('G17', $InvariantCulture)
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
        '-pbe-surface=vdw',
        '-pbe-timing=summary',
        '-pbe-tabi-sum=direct',
        '-pbe-tabi-formulation=jcp2013',
        "-pbe-tabi-energy-projection=$Projection",
        "-pbe-salt=$saltText",
        '-nt=1'
    )
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Invoke-AnalyticCase([double]$Salt, [string]$Projection) {
    $safeSalt = ($Salt.ToString('G17', $InvariantCulture) -replace '\.', 'p')
    $runDir = Join-Path $WorkRoot ("salt_{0}_{1}" -f $safeSalt, $Projection)
    New-CleanDir $runDir
    $argPath = New-ArgFile -RunDir $runDir -Salt $Salt -Projection $Projection
    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for salt=$Salt projection=$Projection"
        }
    } finally {
        Pop-Location
    }

    $energyPath = Join-Path $runDir 'pbe_energy.csv'
    $auditPath = Join-Path $runDir 'pbe_tabi_formulation_audit.csv'
    foreach ($required in @($energyPath, $auditPath)) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Missing TABI Kirkwood output: $required"
        }
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $audit = @(Import-Csv -LiteralPath $auditPath)[0]
    $actual = Convert-ToDouble $energy.E_rxn_kJmol
    $expected = Get-KirkwoodEnergyKJ -Charge 1.0 -RadiusA 3.0 -EpsIn 1.0 -EpsOut 80.0 -Salt $Salt -TempK 298.15
    $rel = [math]::Abs($actual - $expected) / [math]::Max(1.0e-12, [math]::Abs($expected))
    $boundary = Convert-ToDouble $audit.boundary_energy_kJmol
    $boundaryRel = [math]::Abs($boundary - $expected) / [math]::Max(1.0e-12, [math]::Abs($expected))
    return [PSCustomObject]@{
        salt_M = $Salt.ToString('G17', $InvariantCulture)
        projection = $Projection
        energy_kJmol = $actual.ToString('G17', $InvariantCulture)
        kirkwood_kJmol = $expected.ToString('G17', $InvariantCulture)
        rel_error = $rel.ToString('G17', $InvariantCulture)
        boundary_energy_kJmol = $boundary.ToString('G17', $InvariantCulture)
        boundary_rel_error = $boundaryRel.ToString('G17', $InvariantCulture)
        reaction_energy_source = $energy.reaction_energy_source
        fallback_used = $audit.fallback_used
        status = $energy.status
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
New-CleanDir $WorkRoot

$rows = New-Object System.Collections.Generic.List[object]
foreach ($salt in $SaltM) {
    $boundary = Invoke-AnalyticCase -Salt $salt -Projection 'boundary'
    $check = Invoke-AnalyticCase -Salt $salt -Projection 'born-check'
    $rows.Add($boundary)
    $rows.Add($check)
    if ($check.reaction_energy_source -ne 'tabi_born_check') {
        throw "Born-check projection did not report tabi_born_check for salt=$salt"
    }
    if ((Convert-ToDouble $check.rel_error) -gt 1.0e-8) {
        throw "TABI Kirkwood born-check mismatch for salt=$salt rel=$($check.rel_error)"
    }
    if ($boundary.reaction_energy_source -eq 'tabi_born_fallback' -or $boundary.fallback_used -eq 'true') {
        throw "TABI boundary run used fallback for salt=$salt"
    }
    if ((Convert-ToDouble $boundary.boundary_rel_error) -gt 5.0e-2) {
        throw "TABI Kirkwood boundary mismatch for salt=$salt rel=$($boundary.boundary_rel_error)"
    }
    Write-Host ("PBE TABI Kirkwood PASS salt={0} born_check_rel={1} boundary_rel_report={2}" -f
        $check.salt_M, $check.rel_error, $boundary.boundary_rel_error)
}

$outCsv = Join-Path $WorkRoot 'pbe_tabi_kirkwood_analytic.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE TABI Kirkwood analytic report: $outCsv"
