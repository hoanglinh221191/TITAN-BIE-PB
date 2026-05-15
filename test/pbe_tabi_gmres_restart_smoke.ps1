param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_gmres_restart_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI GMRES restart smoke root: $full"
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

function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}

function Run-Case([int]$Restart, [int]$Threads) {
    $runDir = Join-Path $WorkRoot ("restart{0}_nt{1}" -f $Restart, $Threads)
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $pqrPath = Join-Path $InputRoot 'sphere.pqr'
    if (-not (Test-Path -LiteralPath $pqrPath)) {
        throw "Missing sphere PQR fixture: $pqrPath"
    }
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-ToTitanPath $pqrPath)",
        '-pbe-surface=vdw',
        '-pbe-probe=1.4',
        '-pbe-salt=0.15',
        '-pbe-max-iter=200',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-timing=detail',
        '-pbe-tabi-progress=auto',
        '-pbe-tabi-progress-interval=1',
        '-pbe-tabi-sum=direct',
        '-pbe-tabi-formulation=jcp2013',
        '-pbe-tabi-near-integration=centroid',
        '-pbe-tabi-surface-points=24',
        "-pbe-tabi-gmres-restart=$Restart",
        "-nt=$Threads"
    ) -Encoding ASCII
    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE for restart=$Restart nt=$Threads"
        }
    } finally {
        Pop-Location
    }
    $timing = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_timing.csv'))[0]
    $progress = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_longrun_progress.csv'))
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    foreach ($col in @('tabi_gmres_restart_requested','tabi_gmres_restart_effective',
                       'tabi_best_residual','tabi_best_residual_iter','tabi_final_residual')) {
        if (-not ($timing.PSObject.Properties.Name -contains $col)) {
            throw "Missing timing column $col"
        }
    }
    foreach ($col in @('gmres_restart_requested','gmres_restart_effective',
                       'best_residual','best_residual_iter','final_residual')) {
        if (-not ($progress[0].PSObject.Properties.Name -contains $col)) {
            throw "Missing progress column $col"
        }
    }
    if ([int]$timing.tabi_gmres_restart_requested -ne $Restart) {
        throw "Restart request was not reported correctly."
    }
    if ([int]$timing.tabi_gmres_restart_effective -ne $Restart) {
        throw "Restart effective value was not reported correctly."
    }
    if ($progress[-1].phase -ne 'solve_done') {
        throw "Progress CSV did not finish with solve_done."
    }
    return [PSCustomObject]@{
        restart = $Restart
        threads = $Threads
        status = $energy.status
        energy_kJmol = Convert-ToDouble $energy.E_rxn_kJmol
        best_residual = Convert-ToDouble $timing.tabi_best_residual
        best_iter = [int]$timing.tabi_best_residual_iter
        final_residual = Convert-ToDouble $timing.tabi_final_residual
        run_dir = $runDir
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
New-CleanDir $WorkRoot

$rows = @()
foreach ($restart in @(32,64,100)) {
    $rows += Run-Case -Restart $restart -Threads 1
}
$rows += Run-Case -Restart 64 -Threads 4

$nt1 = @($rows | Where-Object { $_.restart -eq 64 -and $_.threads -eq 1 })[0]
$nt4 = @($rows | Where-Object { $_.restart -eq 64 -and $_.threads -eq 4 })[0]
if ([math]::Abs($nt1.energy_kJmol - $nt4.energy_kJmol) -gt 1.0e-9) {
    throw "Restart smoke is not deterministic across nt=1/nt=4."
}

$outCsv = Join-Path $WorkRoot 'pbe_tabi_gmres_restart_smoke.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host "PBE TABI GMRES restart smoke PASS: $outCsv"
