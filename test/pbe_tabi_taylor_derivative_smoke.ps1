param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_taylor_derivative_smoke'),
    [string]$Fixture = 'arg',
    [double]$RelativeTolerance = 2.0e-3
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot
$InputRoot = Join-Path $PSScriptRoot 'pbe_inputs_local'

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }
function Convert-ToDouble([object]$Value) {
    return [double]::Parse([string]$Value, $InvariantCulture)
}
function New-CleanDir([string]$Path) {
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}
function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Encoding ASCII -Value @(
        'CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    )
}
function Invoke-TreecodeDerivative([string]$DerivativeMode) {
    $runDir = Join-Path $WorkRoot $DerivativeMode
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $pqrPath = Join-Path $InputRoot "$Fixture.pqr"
    if (-not (Test-Path -LiteralPath $pqrPath)) { throw "Missing fixture PQR: $pqrPath" }
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
        '-pbe-timing=detail',
        '-pbe-max-iter=800',
        '-pbe-rms=1e-4',
        '-pbe-tabi-taylor-order=3',
        "-pbe-tabi-derivatives=$DerivativeMode",
        '-pbe-tabi-formulation=jcp2013-treecode',
        '-pbe-tabi-sum=treecode',
        '-nt=1'
    )
    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE for $DerivativeMode" }
    } finally {
        Pop-Location
    }
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    $timing = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_timing.csv'))[0]
    $audit = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_formulation_audit.csv'))[0]
    if ($energy.status -ne 'ok') { throw "TABI derivative mode $DerivativeMode failed: $($energy.status)" }
    if ($audit.derivative_mode -ne $DerivativeMode) {
        throw "Expected derivative mode $DerivativeMode, got $($audit.derivative_mode)"
    }
    return [PSCustomObject]@{
        derivative_mode = $DerivativeMode
        energy_kJmol = Convert-ToDouble $energy.E_rxn_kJmol
        matvec_mode = $timing.tabi_matvec_mode
        far_interactions = [int64]$timing.tabi_far_interactions
        total_seconds = Convert-ToDouble $timing.total_seconds
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot
$analytic = Invoke-TreecodeDerivative -DerivativeMode 'analytic'
$finite = Invoke-TreecodeDerivative -DerivativeMode 'finite-diff-debug'
if ($analytic.matvec_mode -ne 'jcp2013_treecode') { throw "Expected analytic JCP treecode, got $($analytic.matvec_mode)" }
if ($analytic.far_interactions -le 0) { throw 'Analytic treecode reported no far interactions.' }
$rel = [math]::Abs($analytic.energy_kJmol - $finite.energy_kJmol) / [math]::Max(1.0e-12, [math]::Abs($finite.energy_kJmol))
if ($rel -gt $RelativeTolerance) {
    throw "Analytic-vs-finite-diff Taylor energy mismatch rel=$rel tolerance=$RelativeTolerance"
}
$outCsv = Join-Path $WorkRoot 'pbe_tabi_taylor_derivative_smoke.csv'
@(
    [PSCustomObject]@{
        fixture = $Fixture
        analytic_energy_kJmol = $analytic.energy_kJmol.ToString('G17', $InvariantCulture)
        finite_diff_energy_kJmol = $finite.energy_kJmol.ToString('G17', $InvariantCulture)
        relative_error = $rel.ToString('G17', $InvariantCulture)
        analytic_seconds = $analytic.total_seconds.ToString('G17', $InvariantCulture)
        finite_diff_seconds = $finite.total_seconds.ToString('G17', $InvariantCulture)
        far_interactions = $analytic.far_interactions
    }
) | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host ("PBE TABI Taylor derivative smoke PASS: rel={0} report={1}" -f
    $rel.ToString('G6', $InvariantCulture), $outCsv)
