param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_presolver_progress_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside pre-solver progress smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Value @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    ) -Encoding ASCII
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}

New-CleanDir $WorkRoot
$pdbPath = Join-Path $WorkRoot 'driver.pdb'
$argPath = Join-Path $WorkRoot 'arg.in'
$pqrPath = Join-Path $InputRoot 'mini-protein.pqr'
$libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
Write-DummyPdb $pdbPath
if (-not (Test-Path -LiteralPath $pqrPath)) {
    throw "Missing mini-protein PQR fixture: $pqrPath"
}

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
    '-pbe-probe=1.4',
    '-pbe-scale=0.35',
    '-pbe-margin=2.0',
    '-pbe-max-grid=300000',
    '-pbe-map=off',
    '-pbe-max-iter=60',
    '-pbe-rms=1e-4',
    '-pbe-dphi=1e-4',
    '-pbe-timing=detail',
    '-pbe-memory-profile=on',
    '-pbe-int-diel=1.0',
    '-pbe-ext-diel=80.0',
    '-pbe-salt=0.15',
    '-pbe-temp=298.15',
    '-nt=2',
    '-pbe-tabi-sum=treecode',
    '-pbe-tabi-formulation=jcp2013-treecode',
    '-pbe-tabi-precond=diagonal',
    '-pbe-tabi-surface-points=16',
    '-pbe-tabi-taylor-order=3',
    '-pbe-tabi-tree-leaf-size=64',
    '-pbe-tabi-opening-theta=0.6',
    '-pbe-tabi-node-bounds=tight-aabb',
    '-pbe-tabi-tree-builder=topdown-aabb',
    '-pbe-tabi-mac-safety=physical-gap',
    '-pbe-tabi-expansion-center=area-centroid',
    '-pbe-tabi-fragmentation-safety=split',
    '-pbe-tabi-near-integration=duffy_adaptive',
    '-pbe-tabi-quad-buffer=adaptive',
    '-pbe-tabi-interaction-mode=onthefly',
    '-pbe-tabi-progress=auto',
    '-pbe-tabi-progress-interval=5',
    '-pbe-tabi-pre-solver-progress=on',
    '-pbe-tabi-convergence-estimator=hybrid',
    '-pbe-tabi-convergence-safety=split',
    '-pbe-tabi-convergence-calibration=consecutive-aspect',
    '-pbe-tabi-convergence-dielectric-amp=on',
    '-pbe-tabi-convergence-estimator-model=channel-tail-eps',
    '-pbe-tabi-convergence-channel-calibration=k3-sqrt-eps',
    '-pbe-tabi-expansion-audit=sampled',
    '-pbe-tabi-expansion-audit-channel=summary',
    '-pbe-tabi-expansion-audit-sample-count=64',
    '-pbe-tabi-expansion-audit-max-panels=512'
)
Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII

$stdoutPath = Join-Path $WorkRoot 'stdout.log'
$stderrPath = Join-Path $WorkRoot 'stderr.log'
$psi = [System.Diagnostics.ProcessStartInfo]::new()
$psi.FileName = $TitanExe
$psi.WorkingDirectory = $WorkRoot
$psi.UseShellExecute = $false
$psi.CreateNoWindow = $true
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$psi.Arguments = ('-a "{0}"' -f ($argPath -replace '"','\"'))
$proc = [System.Diagnostics.Process]::new()
$proc.StartInfo = $psi
[void]$proc.Start()
$stdoutTask = $proc.StandardOutput.ReadToEndAsync()
$stderrTask = $proc.StandardError.ReadToEndAsync()
if (-not $proc.WaitForExit(120000)) {
    $proc.Kill()
    $proc.WaitForExit()
    throw "TABI pre-solver progress smoke timed out."
}
$stdoutTask.Wait()
$stderrTask.Wait()
[System.IO.File]::WriteAllText($stdoutPath, $stdoutTask.Result)
[System.IO.File]::WriteAllText($stderrPath, $stderrTask.Result)
if ($proc.ExitCode -ne 0) {
    throw "TABI pre-solver progress smoke failed with exit code $($proc.ExitCode)."
}

$progressPath = Join-Path $WorkRoot 'pbe_tabi_prepare_progress.csv'
if (-not (Test-Path -LiteralPath $progressPath)) {
    throw "Missing pbe_tabi_prepare_progress.csv."
}
$rows = @(Import-Csv -LiteralPath $progressPath)
$phases = @($rows | ForEach-Object { $_.phase })
foreach ($required in @('tabi_prepare_start','pqr_loaded','tabi_options_built','solver_enter')) {
    if ($phases -notcontains $required) {
        throw "Missing pre-solver progress phase: $required"
    }
}
$solverRow = @($rows | Where-Object { $_.phase -eq 'solver_enter' } | Select-Object -First 1)
if ($solverRow.Count -eq 0 -or $solverRow[0].tabi_convergence_channel_calibration -ne 'k3-sqrt-eps') {
    throw "Pre-solver progress did not preserve Stage 3N channel calibration."
}

Write-Host ("PBE TABI pre-solver progress smoke PASS phases={0} rss={1}" -f (($phases -join ';')), $solverRow[0].rss_bytes)
