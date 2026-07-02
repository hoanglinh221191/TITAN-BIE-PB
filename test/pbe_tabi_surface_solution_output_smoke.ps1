param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_surface_solution_output_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

function Write-DriverPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   GLY A   1      -1.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  GLY A   1       0.300   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   GLY A   1       0.900   1.450   0.000  1.00  0.00           C',
        'END'
    ) -Encoding ASCII
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI surface solution output root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-SpherePqr([string]$Path) {
    Set-Content -LiteralPath $Path -Value @(
        'ATOM      1  SHL SPH     1       0.000   0.000   0.000  0.0000  4.0000       X',
        'ATOM      2  Q1  ION     2       1.250   0.000   0.000  1.0000  0.2000       X',
        'END'
    ) -Encoding ASCII
}

function Write-MinimalForceField([string]$Root) {
    New-Item -ItemType Directory -Path $Root -Force | Out-Null
    Set-Content -LiteralPath (Join-Path $Root 'ffAMBER99SBnonbonded.dat') -Value @(
        '; atom type  symbol  mass  charge  ptype  c6  c12',
        'C C 12.011 0.0 A 0.001 0.000001',
        'N N 14.007 0.0 A 0.001 0.000001'
    ) -Encoding ASCII
    Set-Content -LiteralPath (Join-Path $Root 'AMBER99SB_aminoacids.rtp') -Value @(
        '[ GLY ]',
        ' [ atoms ]',
        ' N  N  0.0  1',
        ' CA C  0.0  1',
        ' C  C  0.0  1',
        ' [ bonds ]'
    ) -Encoding ASCII
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-CleanDir $WorkRoot

$driver = Join-Path $WorkRoot 'driver.pdb'
$pqr = Join-Path $WorkRoot 'sphere.pqr'
$ffRoot = Join-Path $WorkRoot 'lib_ff'
Write-DriverPdb $driver
Write-SpherePqr $pqr
Write-MinimalForceField $ffRoot
$libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath $ffRoot).Path.TrimEnd('\') + '\')
$argLines = @(
    "-f=$(Convert-ToTitanPath $driver)",
    '-1pdb=true',
    '-ff=AMBER99SB',
    '-g=1&1',
    "-path-a=$libPath",
    "-path-n=$libPath",
    '-pbe=tabi',
    "-pbe-test-pqr=$(Convert-ToTitanPath $pqr)",
    '-pbe-surface=vdw',
    '-pbe-probe=0.0',
    '-pbe-salt=0.05',
    '-pbe-temp=298.15',
    '-pbe-timing=detail',
    '-pbe-map=off',
    '-pbe-max-iter=300',
    '-pbe-rms=1e-4',
    '-pbe-tabi-sum=direct',
    '-pbe-tabi-formulation=jcp2013',
    '-pbe-tabi-energy-projection=boundary',
    '-pbe-tabi-near-integration=analytic_singular',
    '-pbe-tabi-surface-points=80',
    '-pbe-tabi-sphere-mesh=icosphere',
    '-pbe-tabi-surface-output=all',
    '-nt=1'
)
$argPath = Join-Path $WorkRoot 'arg.in'
Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII

$badArgLines = $argLines -match '^-{1,2}[^=]+$'
if ($badArgLines.Count -gt 0) {
    throw "arg.in contains bare option lines: $($badArgLines -join ', ')"
}

Push-Location $WorkRoot
try {
    & $TitanExe -a $argPath | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE" }
} finally {
    Pop-Location
}

$qualityPath = Join-Path $WorkRoot 'pbe_tabi_surface_quality.csv'
$csvPath = Join-Path $WorkRoot 'pbe_tabi_surface_solution.csv'
$vtkPath = Join-Path $WorkRoot 'pbe_tabi_surface_solution.vtk'
$summaryPath = Join-Path $WorkRoot 'pbe_summary.dat'
foreach ($path in @($qualityPath, $csvPath, $vtkPath, $summaryPath)) {
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing expected output: $path" }
}

$quality = @(Import-Csv -LiteralPath $qualityPath)[0]
$rows = @(Import-Csv -LiteralPath $csvPath)
$panelCount = [int64]$quality.panel_count
if ($rows.Count -ne $panelCount) {
    throw "Expected $panelCount surface solution rows, got $($rows.Count)"
}
foreach ($column in @('panel_id','x_A','y_A','z_A','nx','ny','nz','area_A2','phi','phin','normal_flux_proxy')) {
    if (-not ($rows[0].PSObject.Properties.Name -contains $column)) {
        throw "Missing CSV column: $column"
    }
}
$finiteCount = 0
foreach ($row in $rows) {
    $phi = Convert-ToDouble $row.phi
    $phin = Convert-ToDouble $row.phin
    $flux = Convert-ToDouble $row.normal_flux_proxy
    $area = Convert-ToDouble $row.area_A2
    if (-not ([double]::IsNaN($phi) -or [double]::IsInfinity($phi) -or
              [double]::IsNaN($phin) -or [double]::IsInfinity($phin) -or
              [double]::IsNaN($flux) -or [double]::IsInfinity($flux))) {
        $finiteCount++
    }
    if ([Math]::Abs($flux - $phin) -gt 1e-12) {
        throw "normal_flux_proxy does not match phin in surface solution CSV"
    }
    if ($area -le 0.0) { throw "Non-positive panel area in surface output" }
}
if ($finiteCount -ne $rows.Count) {
    throw "Non-finite phi/phin values in surface solution CSV"
}

$vtkText = Get-Content -LiteralPath $vtkPath -Raw
if ($vtkText -notmatch 'DATASET POLYDATA' -or
    $vtkText -notmatch "CELL_DATA $panelCount" -or
    $vtkText -notmatch 'SCALARS phi double 1' -or
    $vtkText -notmatch 'SCALARS phin double 1' -or
    $vtkText -notmatch 'SCALARS normal_flux_proxy double 1') {
    throw "VTK surface output does not contain the expected PolyData fields"
}
$summary = Get-Content -LiteralPath $summaryPath -Raw
if ($summary -notmatch 'tabi_surface_solution_csv pbe_tabi_surface_solution.csv' -or
    $summary -notmatch 'tabi_surface_solution_vtk pbe_tabi_surface_solution.vtk') {
    throw "pbe_summary.dat does not advertise TABI surface solution outputs"
}

Write-Host "PBE TABI surface solution output smoke PASS: $csvPath"
