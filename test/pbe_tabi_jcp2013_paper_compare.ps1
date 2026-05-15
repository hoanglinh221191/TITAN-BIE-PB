param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_jcp2013_paper_compare'),
    [int[]]$SurfacePoints = @(96, 384),
    [double[]]$SaltM = @(0.0, 0.15)
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }
function To-Double([object]$Value) { return [double]::Parse([string]$Value, $InvariantCulture) }
function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    if (-not ($full.Equals($rootBase, [StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootBase + '\', [StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside JCP paper compare root: $full"
    }
}
function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}
function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Encoding ASCII -Value @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    )
}
function Get-Kappa([double]$Salt, [double]$TempK) {
    if ($Salt -le 0.0) { return 0.0 }
    return 0.3291 * [math]::Sqrt($Salt) * [math]::Sqrt(298.15 / [math]::Max(1.0, $TempK))
}
function Get-KirkwoodEnergyKJ([double]$Salt) {
    $radius = 3.0
    $kappa = Get-Kappa -Salt $Salt -TempK 298.15
    $outside = 1.0 / (80.0 * (1.0 + $kappa * $radius))
    $inside = 1.0
    return 0.5 * 332.063713299 * ($outside - $inside) / $radius * 4.184
}
function Invoke-JcpSphere([int]$Points, [double]$Salt) {
    $saltTag = ($Salt.ToString('G17', $InvariantCulture) -replace '\.', 'p')
    $runDir = Join-Path $WorkRoot ("sphere_points_{0}_salt_{1}" -f $Points, $saltTag)
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $pqrPath = Join-Path $InputRoot 'sphere.pqr'
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Encoding ASCII -Value @(
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
        '-pbe-tabi-formulation=jcp2013-centroid',
        "-pbe-tabi-surface-points=$Points",
        ("-pbe-salt={0}" -f $Salt.ToString('G17', $InvariantCulture)),
        '-nt=1'
    )
    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE" }
    } finally {
        Pop-Location
    }
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    $audit = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_formulation_audit.csv'))[0]
    $memory = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_memory_profile.csv'))[0]
    $expected = Get-KirkwoodEnergyKJ -Salt $Salt
    $actual = To-Double $energy.E_rxn_kJmol
    $rel = [math]::Abs($actual - $expected) / [math]::Max(1.0e-12, [math]::Abs($expected))
    return [PSCustomObject]@{
        case = 'kirkwood_sphere'
        surface_points = $Points
        salt_M = $Salt.ToString('G17', $InvariantCulture)
        panels = $audit.panels
        unknown_count = $audit.unknown_count
        bie_system = $audit.bie_system
        formulation = $audit.formulation
        energy_kJmol = $actual.ToString('G17', $InvariantCulture)
        kirkwood_kJmol = $expected.ToString('G17', $InvariantCulture)
        rel_error = $rel.ToString('G17', $InvariantCulture)
        total_seconds = $memory.total_seconds
        estimated_bytes_total = $memory.estimated_bytes_total
        status = $energy.status
        reaction_energy_source = $energy.reaction_energy_source
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot
$rows = New-Object System.Collections.Generic.List[object]
foreach ($points in $SurfacePoints) {
    foreach ($salt in $SaltM) {
        $row = Invoke-JcpSphere -Points $points -Salt $salt
        $rows.Add($row)
        $rel = To-Double $row.rel_error
        $threshold = if ($points -ge 384) { 2.0e-2 } else { 5.0e-2 }
        if ($rel -gt $threshold) {
            throw "JCP Kirkwood sphere error too high points=$points salt=$salt rel=$rel threshold=$threshold"
        }
        if ($row.reaction_energy_source -ne 'tabi_boundary_integral') {
            throw "JCP paper compare used non-boundary source: $($row.reaction_energy_source)"
        }
    }
}
$outCsv = Join-Path $WorkRoot 'pbe_tabi_jcp2013_paper_compare.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE TABI JCP 2013 paper compare PASS: $outCsv"
