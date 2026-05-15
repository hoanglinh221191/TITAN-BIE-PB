param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_low_density_surface_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI low-density smoke root: $full"
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
    $lines = @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    )
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Run-Case([int]$SurfacePoints, [int]$Threads) {
    $runDir = Join-Path $WorkRoot ("sp{0}_nt{1}" -f $SurfacePoints, $Threads)
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $pqrPath = Join-Path $InputRoot 'sphere.pqr'
    if (-not (Test-Path -LiteralPath $pqrPath)) {
        throw "Missing sphere PQR fixture: $pqrPath"
    }
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argPath = Join-Path $runDir 'arg.in'
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
        '-pbe-salt=0.15',
        '-pbe-max-iter=400',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-timing=detail',
        '-pbe-tabi-progress=auto',
        '-pbe-tabi-progress-interval=1',
        '-pbe-tabi-sum=direct',
        '-pbe-tabi-formulation=jcp2013',
        '-pbe-tabi-near-integration=centroid',
        "-pbe-tabi-surface-points=$SurfacePoints",
        "-nt=$Threads"
    )
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for sp=$SurfacePoints nt=$Threads"
        }
    } finally {
        Pop-Location
    }
    $timingPath = Join-Path $runDir 'pbe_timing.csv'
    $progressPath = Join-Path $runDir 'pbe_tabi_longrun_progress.csv'
    if (-not (Test-Path -LiteralPath $timingPath)) {
        throw "Missing timing CSV for sp=$SurfacePoints nt=$Threads"
    }
    if (-not (Test-Path -LiteralPath $progressPath)) {
        throw "Missing progress CSV for sp=$SurfacePoints nt=$Threads"
    }
    $timing = @(Import-Csv -LiteralPath $timingPath)[0]
    if ($timing.tabi_status -ne 'ok') {
        throw "Low-density TABI did not finish ok for sp=$SurfacePoints nt=${Threads}: $($timing.tabi_status)"
    }
    if ([int64]$timing.tabi_surface_panels -le 0) {
        throw "Non-positive TABI panel count for sp=$SurfacePoints nt=$Threads"
    }
    if ([int]$timing.tabi_surface_point_count_requested -ne $SurfacePoints) {
        throw "Requested surface point count was not reported correctly."
    }
    if ([double]::Parse($timing.tabi_panel_density_per_A2, $InvariantCulture) -le 0.0) {
        throw "Non-positive panel density for sp=$SurfacePoints nt=$Threads"
    }
    $progress = @(Import-Csv -LiteralPath $progressPath)
    if ($progress.Count -lt 1 -or $progress[-1].phase -ne 'solve_done') {
        throw "Progress CSV did not record solve_done for sp=$SurfacePoints nt=$Threads"
    }
    return [PSCustomObject]@{
        surface_points = $SurfacePoints
        threads = $Threads
        panels = [int64]$timing.tabi_surface_panels
        density = [double]::Parse($timing.tabi_panel_density_per_A2, $InvariantCulture)
        status = $timing.tabi_panel_budget_status
        energy = [double]::Parse((@(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]).E_rxn_kJmol, $InvariantCulture)
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}

New-CleanDir $WorkRoot
$rows = @()
foreach ($points in @(8,12,16,24)) {
    $one = Run-Case -SurfacePoints $points -Threads 1
    $four = Run-Case -SurfacePoints $points -Threads 4
    if ($one.panels -ne $four.panels) {
        throw "Panel count is not deterministic for surface_points=$points"
    }
    if ([math]::Abs($one.energy - $four.energy) -gt 1.0e-9) {
        throw "Energy is not deterministic for surface_points=$points"
    }
    $rows += $one
    $rows += $four
}

$outCsv = Join-Path $WorkRoot 'pbe_tabi_low_density_surface_smoke.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE TABI low-density surface smoke PASS: $outCsv"
