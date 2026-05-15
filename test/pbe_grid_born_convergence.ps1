param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_grid_born_convergence'),
    [double[]]$RadiiA = @(1.5, 2.0, 3.0, 5.0),
    [double[]]$Charges = @(1.0, -1.0),
    [double[]]$SaltM = @(0.0, 0.15, 0.5),
    [double[]]$GridScales = @(0.35, 0.60, 1.0),
    [ValidateSet('off','auto')]
    [string[]]$FocusModes = @('off'),
    [ValidateSet('auto','fixed-center')]
    [string[]]$GridAnchors = @('auto','fixed-center')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside grid Born convergence root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Get-RowValue([object]$Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) {
            return $value
        }
    }
    return $Default
}

function Convert-ToDoubleOrNull([object]$Value) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') {
        return $null
    }
    return [double]::Parse($text, $InvariantCulture)
}

function Get-BornReferenceKJmol([double]$Charge, [double]$RadiusA, [double]$SaltM) {
    $epsIn = 1.0
    $epsOut = 80.0
    $kappa = if ($SaltM -le 0.0) { 0.0 } else { 0.3291 * [Math]::Sqrt($SaltM) }
    $outside = 1.0 / ($epsOut * (1.0 + $kappa * $RadiusA))
    $inside = 1.0 / $epsIn
    $kcal = 0.5 * 332.063713299 * $Charge * $Charge * ($outside - $inside) / $RadiusA
    return $kcal * 4.184
}

function Get-LinearSlope([object[]]$Rows) {
    $points = @($Rows | Where-Object {
        (Convert-ToDoubleOrNull $_.rel_error_vs_analytic) -ne $null -and
        (Convert-ToDoubleOrNull $_.grid_spacing_A) -ne $null -and
        (Convert-ToDoubleOrNull $_.rel_error_vs_analytic) -gt 0.0 -and
        (Convert-ToDoubleOrNull $_.grid_spacing_A) -gt 0.0
    })
    if ($points.Count -lt 2) {
        return 'NA'
    }
    $xs = @($points | ForEach-Object { [Math]::Log((Convert-ToDoubleOrNull $_.grid_spacing_A)) })
    $ys = @($points | ForEach-Object { [Math]::Log((Convert-ToDoubleOrNull $_.rel_error_vs_analytic)) })
    $n = [double]$xs.Count
    $sx = ($xs | Measure-Object -Sum).Sum
    $sy = ($ys | Measure-Object -Sum).Sum
    $sxx = ($xs | ForEach-Object { $_ * $_ } | Measure-Object -Sum).Sum
    $sxy = 0.0
    for ($i = 0; $i -lt $xs.Count; ++$i) {
        $sxy += $xs[$i] * $ys[$i]
    }
    $den = $n * $sxx - $sx * $sx
    if ([Math]::Abs($den) -lt 1.0e-14) {
        return 'NA'
    }
    return (($n * $sxy - $sx * $sy) / $den).ToString($InvariantCulture)
}

function Write-DummyPdb([string]$Path) {
    $lines = @(
        'CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Write-BornPqr([string]$Path, [double]$Charge, [double]$RadiusA) {
    $chargeText = $Charge.ToString('F4', $InvariantCulture)
    $radiusText = $RadiusA.ToString('F4', $InvariantCulture)
    Set-Content -LiteralPath $Path -Value @(
        ("ATOM      1  ION ION     1       0.000   0.000   0.000{0,8}{1,8}       X" -f $chargeText, $radiusText),
        'END'
    ) -Encoding ASCII
}

function New-ArgFile([string]$RunDir, [double]$Charge, [double]$RadiusA,
                     [double]$SaltM, [double]$GridScale,
                     [string]$FocusMode, [string]$GridAnchor) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    $pqrPath = Join-Path $RunDir 'born_ion.pqr'
    Write-DummyPdb $pdbPath
    Write-BornPqr -Path $pqrPath -Charge $Charge -RadiusA $RadiusA
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=grid',
        "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
        '-pbe-surface=vdw',
        '-pbe-probe=0.0',
        "-pbe-salt=$($SaltM.ToString($InvariantCulture))",
        '-pbe-temp=298.15',
        '-pbe-timing=summary',
        '-pbe-map=off',
        '-pbe-margin=2.0',
        "-pbe-scale=$($GridScale.ToString($InvariantCulture))",
        '-pbe-max-grid=550000',
        '-pbe-max-iter=1200',
        '-pbe-rms=1e-5',
        '-nt=1'
    )
    if ($GridAnchor -eq 'fixed-center') {
        $argLines += '-pbe-grid-anchor=fixed-center'
        $argLines += '-pbe-test-center=0,0,0'
    }
    if ($FocusMode -eq 'auto') {
        $argLines += '-pbe-focus=auto'
        $argLines += '-pbe-focus-levels=2'
        $argLines += '-pbe-focus-margin=4.0'
        if ($GridAnchor -eq 'fixed-center') {
            $argLines += '-pbe-focus-center=0,0,0'
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    return $argPath
}

function Invoke-GridCase([double]$Charge, [double]$RadiusA, [double]$SaltM,
                         [double]$GridScale, [string]$FocusMode,
                         [string]$GridAnchor) {
    $safeSalt = ([string]$SaltM).Replace('.', 'p').Replace('-', 'm')
    $safeScale = ([string]$GridScale).Replace('.', 'p').Replace('-', 'm')
    $safeCharge = ([string]$Charge).Replace('.', 'p').Replace('-', 'm')
    $safeRadius = ([string]$RadiusA).Replace('.', 'p').Replace('-', 'm')
    $runDir = Join-Path $WorkRoot ("grid_q{0}_r{1}_s{2}_scale{3}_{4}_{5}" -f $safeCharge, $safeRadius, $safeSalt, $safeScale, $FocusMode, $GridAnchor)
    New-CleanDir $runDir
    $argPath = New-ArgFile -RunDir $runDir -Charge $Charge -RadiusA $RadiusA `
        -SaltM $SaltM -GridScale $GridScale -FocusMode $FocusMode `
        -GridAnchor $GridAnchor
    $measure = Measure-Command {
        Push-Location $runDir
        try {
            & $TitanExe -a $argPath | Out-Null
            if ($LASTEXITCODE -ne 0) {
                throw "TITAN exited with code $LASTEXITCODE"
            }
        } finally {
            Pop-Location
        }
    }
    $energyPath = Join-Path $runDir 'pbe_energy.csv'
    $timingPath = Join-Path $runDir 'pbe_timing.csv'
    if (-not (Test-Path -LiteralPath $energyPath)) {
        throw "Missing pbe_energy.csv: $energyPath"
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $timing = if (Test-Path -LiteralPath $timingPath) { @(Import-Csv -LiteralPath $timingPath)[0] } else { $null }
    $memoryPath = Join-Path $runDir 'pbe_memory_profile.csv'
    $memory = if (Test-Path -LiteralPath $memoryPath) { @(Import-Csv -LiteralPath $memoryPath)[0] } else { $null }
    if ($energy.status -ne 'ok') {
        throw "Grid Born case did not finish ok: $($energy.status)"
    }
    $analytic = Get-BornReferenceKJmol -Charge $Charge -RadiusA $RadiusA -SaltM $SaltM
    $energyKJ = Convert-ToDoubleOrNull $energy.E_rxn_kJmol
    if ($null -eq $energyKJ -or [double]::IsNaN($energyKJ) -or [double]::IsInfinity($energyKJ)) {
        throw "Non-finite grid Born energy: $($energy.E_rxn_kJmol)"
    }
    $rel = [Math]::Abs($energyKJ - $analytic) / [Math]::Max(1.0e-12, [Math]::Abs($analytic))
    $vacuum = Convert-ToDoubleOrNull (Get-RowValue $energy 'vacuum_phase_rf_kJmol' 'NA')
    $vacStatus = if ($null -ne $vacuum) { 'reported' } else { 'missing' }
    return [PSCustomObject]@{
        backend = 'grid'
        charge_e = $Charge.ToString($InvariantCulture)
        radius_A = $RadiusA.ToString($InvariantCulture)
        salt_M = $SaltM.ToString($InvariantCulture)
        grid_scale = $GridScale.ToString($InvariantCulture)
        grid_spacing_A = (1.0 / $GridScale).ToString($InvariantCulture)
        focus_mode = $FocusMode
        grid_anchor = $GridAnchor
        analytic_kJmol = $analytic.ToString($InvariantCulture)
        energy_kJmol = $energy.E_rxn_kJmol
        rel_error_vs_analytic = $rel.ToString($InvariantCulture)
        convergence_slope = 'NA'
        vacuum_reference_status = $vacStatus
        reaction_energy_source = $energy.reaction_energy_source
        seconds = (Get-RowValue $timing 'total_seconds' ([Math]::Round($measure.TotalSeconds, 6).ToString($InvariantCulture)))
        solver_seconds = (Get-RowValue $timing 'solver_seconds' 'NA')
        peak_ram = (Get-RowValue $memory 'peak_working_set_bytes' 'NA')
        estimated_bytes = (Get-RowValue $memory 'estimated_bytes_total' 'NA')
        grid_points_or_panels = $energy.grid_points
        grid_points = $energy.grid_points
        panels = '0'
        status = $energy.status
        run_dir = $runDir
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
New-CleanDir $WorkRoot

$rows = New-Object System.Collections.Generic.List[object]
foreach ($charge in $Charges) {
    foreach ($radius in $RadiiA) {
        foreach ($salt in $SaltM) {
            foreach ($focus in $FocusModes) {
                foreach ($anchor in $GridAnchors) {
                    foreach ($scale in $GridScales) {
                        $rows.Add((Invoke-GridCase -Charge $charge -RadiusA $radius `
                            -SaltM $salt -GridScale $scale -FocusMode $focus `
                            -GridAnchor $anchor))
                    }
                }
            }
        }
    }
}

foreach ($group in ($rows | Group-Object -Property charge_e,radius_A,salt_M,focus_mode,grid_anchor)) {
    $slope = Get-LinearSlope @($group.Group)
    foreach ($row in $group.Group) {
        $row.convergence_slope = $slope
    }
}

$outCsv = Join-Path $WorkRoot 'pbe_grid_born_convergence.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE grid Born convergence diagnostic COMPLETE: $outCsv"
