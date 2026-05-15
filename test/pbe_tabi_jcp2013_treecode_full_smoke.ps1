param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_jcp2013_treecode_full_smoke'),
    [string]$Fixture = 'arg',
    [double]$RelativeTolerance = 5.0e-4
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
function Invoke-TabiJcp([string]$Mode) {
    $runDir = Join-Path $WorkRoot $Mode
    New-CleanDir $runDir
    $pdbPath = Join-Path $runDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $pqrPath = Join-Path $InputRoot "$Fixture.pqr"
    if (-not (Test-Path -LiteralPath $pqrPath)) { throw "Missing fixture PQR: $pqrPath" }
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argPath = Join-Path $runDir 'arg.in'
    $formulation = if ($Mode -eq 'treecode') { 'jcp2013-treecode' } else { 'jcp2013-centroid' }
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
        '-pbe-tabi-derivatives=analytic',
        "-pbe-tabi-formulation=$formulation",
        "-pbe-tabi-sum=$Mode",
        '-nt=1'
    )
    Push-Location $runDir
    try {
        & $TitanExe -a $argPath | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE for $Mode" }
    } finally {
        Pop-Location
    }
    $energy = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_energy.csv'))[0]
    $timing = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_timing.csv'))[0]
    $audit = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_formulation_audit.csv'))[0]
    if ($energy.status -ne 'ok') { throw "TABI JCP $Mode failed: $($energy.status)" }
    return [PSCustomObject]@{
        mode = $Mode
        energy_kJmol = Convert-ToDouble $energy.E_rxn_kJmol
        matvec_mode = $timing.tabi_matvec_mode
        derivative_mode = Get-Member -InputObject $audit -Name derivative_mode -MemberType NoteProperty -ErrorAction SilentlyContinue | ForEach-Object { $audit.derivative_mode }
        tree_nodes = [int64]$timing.tabi_tree_nodes
        near = [int64]$timing.tabi_near_interactions
        far = [int64]$timing.tabi_far_interactions
        total_seconds = Convert-ToDouble $timing.total_seconds
        run_dir = $runDir
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot
$direct = Invoke-TabiJcp -Mode 'direct'
$tree = Invoke-TabiJcp -Mode 'treecode'
$rel = [math]::Abs($tree.energy_kJmol - $direct.energy_kJmol) / [math]::Max(1.0e-12, [math]::Abs($direct.energy_kJmol))
if ($tree.matvec_mode -ne 'jcp2013_treecode') { throw "Expected JCP treecode matvec, got $($tree.matvec_mode)" }
if ($tree.derivative_mode -ne 'analytic') { throw "Expected analytic JCP treecode derivatives, got $($tree.derivative_mode)" }
if ($tree.tree_nodes -le 0 -or $tree.far -le 0) { throw "JCP treecode did not report tree/far interactions." }
if ($rel -gt $RelativeTolerance) {
    throw "JCP treecode/direct energy mismatch rel=$rel tolerance=$RelativeTolerance"
}
$outCsv = Join-Path $WorkRoot 'pbe_tabi_jcp2013_treecode_full_smoke.csv'
@(
    [PSCustomObject]@{
        fixture = $Fixture
        direct_energy_kJmol = $direct.energy_kJmol.ToString('G17', $InvariantCulture)
        treecode_energy_kJmol = $tree.energy_kJmol.ToString('G17', $InvariantCulture)
        relative_error = $rel.ToString('G17', $InvariantCulture)
        tree_nodes = $tree.tree_nodes
        far_interactions = $tree.far
        direct_seconds = $direct.total_seconds.ToString('G17', $InvariantCulture)
        treecode_seconds = $tree.total_seconds.ToString('G17', $InvariantCulture)
    }
) | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII
Write-Host ("PBE TABI JCP2013 treecode full smoke PASS: rel={0} far={1} report={2}" -f
    $rel.ToString('G6', $InvariantCulture), $tree.far, $outCsv)
