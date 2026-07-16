param(
    [string]$TitanExe = (Join-Path $PSScriptRoot '..\..\titan.exe'),
    [string]$ResultsDir = (Join-Path $PSScriptRoot 'results'),
    [ValidateRange(1, 256)] [int]$Threads = 4,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Read-SummaryMap([string]$Path) {
    $map = @{}
    foreach ($line in Get-Content -LiteralPath $Path) {
        if ([string]::IsNullOrWhiteSpace($line) -or $line.StartsWith('#')) { continue }
        $parts = $line -split '\s+', 2
        if ($parts.Count -eq 2) { $map[$parts[0]] = $parts[1] }
    }
    return $map
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing TITAN executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$ResultsDir = [System.IO.Path]::GetFullPath($ResultsDir)

if (Test-Path -LiteralPath $ResultsDir) {
    $existing = @(Get-ChildItem -LiteralPath $ResultsDir -Force)
    if ($existing.Count -gt 0 -and -not $Force) {
        throw "Results directory is not empty: $ResultsDir. Use -Force to replace it."
    }
    if ($Force) { Remove-Item -LiteralPath $ResultsDir -Recurse -Force }
}
New-Item -ItemType Directory -Path $ResultsDir -Force | Out-Null

$argSource = Join-Path $PSScriptRoot 'arg_tabi_mm_pbsa.in'
$argPath = Join-Path $ResultsDir 'arg.in'
$argLines = @(Get-Content -LiteralPath $argSource)
$argLines = @($argLines | ForEach-Object {
    if ($_ -match '^-nt=') { "-nt=$Threads" } else { $_ }
})
if (-not ($argLines -contains '-pbcn')) {
    throw 'The manuscript trajectory example must enable FixBox-style PBC with -pbcn.'
}
$argLines | Set-Content -LiteralPath $argPath -Encoding ASCII

$trajectoryPath = Join-Path $PSScriptRoot 'input\complex_trajectory.pdb'
if (-not (Test-Path -LiteralPath $trajectoryPath)) {
    throw "Missing bundled trajectory: $trajectoryPath"
}
$modelRecords = @(Select-String -LiteralPath $trajectoryPath -Pattern '^MODEL').Count
if ($modelRecords -lt 2) {
    throw "Expected a multi-model trajectory, found $modelRecords MODEL records."
}
# TITAN uses the first model as the coordinate/topology seed. Analysis starts
# from the next timestamp, so N MODEL records produce N-1 sampled frames.
$expectedFrames = $modelRecords - 1

Push-Location $ResultsDir
try {
    & $TitanExe -a $argPath
    if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE" }
} finally {
    Pop-Location
}

$summaryPath = Join-Path $ResultsDir 'mm_pbsa_summary.dat'
if (-not (Test-Path -LiteralPath $summaryPath)) { throw 'Missing mm_pbsa_summary.dat' }
$summary = Read-SummaryMap $summaryPath
if ($summary['pbsa_backend'] -ne 'tabi') { throw 'The example did not use the TABI backend.' }
if ([int]$summary['frames_total'] -ne $expectedFrames) {
    throw "Expected $expectedFrames MM/PBSA frames, found $($summary['frames_total'])."
}
if ($summary['interaction_entropy_status'] -eq 'disabled') {
    throw 'Interaction entropy was unexpectedly disabled for the full trajectory.'
}

$energyPath = Join-Path $ResultsDir 'mm_pbsa_energy.csv'
$energyRows = @(Import-Csv -LiteralPath $energyPath)
if ($energyRows.Count -ne $expectedFrames) {
    throw "Expected $expectedFrames energy rows, found $($energyRows.Count)."
}
if (@($energyRows | Where-Object { $_.status -ne 'ok' }).Count -gt 0) {
    throw 'At least one trajectory frame failed the TABI MM/PBSA cycle.'
}

$logPath = Join-Path $ResultsDir 'log.txt'
if (-not (Test-Path -LiteralPath $logPath) -or
    -not (Select-String -LiteralPath $logPath -Pattern 'Solving PBC' -Quiet)) {
    throw 'The run log does not confirm that the PBC stage was executed.'
}

$decomp = Read-SummaryMap (Join-Path $ResultsDir 'mm_pbsa_decomposition.dat')
if ([math]::Abs([double]$decomp['unit_closure_total_max_abs_kJmol']) -gt 1e-5) {
    throw "MM/PBSA unit closure failed: $($decomp['unit_closure_total_max_abs_kJmol']) kJ/mol"
}

Write-Host ''
Write-Host 'TABI MM/PBSA completed:'
Write-Host "  Trajectory MODEL records:  $modelRecords"
Write-Host "  PBC-corrected frames:      $expectedFrames"
Write-Host "  Delta G_bind (no entropy): $($summary['DeltaG_no_entropy_kcalmol_mean']) kcal/mol"
Write-Host "  -T Delta S (IE):           $($summary['minus_TdS_IE_kcalmol']) kcal/mol"
Write-Host "  Delta G_bind (with IE):    $($summary['DeltaG_binding_IE_kcalmol']) kcal/mol"
Write-Host "  Polar contribution:        $($summary['DeltaDeltaG_solv_PB_kcalmol_mean']) kcal/mol"
Write-Host "  Residue/atom table:         $(Join-Path $ResultsDir 'mm_pbsa_unit_summary.csv')"
Write-Host "  Residue ranking:            $(Join-Path $ResultsDir 'edecomp_mmpbsa_residue_ranking.csv')"
Write-Host "  Atom ranking:               $(Join-Path $ResultsDir 'edecomp_mmpbsa_atom_ranking.csv')"
