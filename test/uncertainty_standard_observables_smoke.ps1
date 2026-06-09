param(
  [string]$TitanExe = ".\titan.exe"
)

$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $RepoRoot
if (-not [System.IO.Path]::IsPathRooted($TitanExe)) {
  $TitanExe = (Resolve-Path -LiteralPath (Join-Path $RepoRoot $TitanExe)).Path
}

$unitExe = Join-Path $RepoRoot "titan_uncertainty_smoke.exe"
& g++ -std=c++23 -O2 -Isrc test\titan_uncertainty_smoke.cpp src\titan_uncertainty.cpp -o $unitExe
if ($LASTEXITCODE -ne 0) { throw "failed to build titan_uncertainty_smoke" }
& $unitExe
if ($LASTEXITCODE -ne 0) { throw "titan_uncertainty_smoke failed" }

function New-WorkDir([string]$Name) {
  $dir = Join-Path $RepoRoot "test_uncertainty_$Name"
  if (Test-Path -LiteralPath $dir) { Remove-Item -LiteralPath $dir -Recurse -Force }
  New-Item -ItemType Directory -Path $dir | Out-Null
  return $dir
}

function Write-ArgFile([string]$Dir, [string[]]$Lines) {
  Set-Content -Path (Join-Path $Dir "arg.in") -Value $Lines -Encoding ASCII
}

function Run-Titan([string]$Dir) {
  Push-Location $Dir
  try {
    & $TitanExe | Out-Host
    if ($LASTEXITCODE -ne 0 -and (Test-Path -LiteralPath "index.ndx")) {
      & $TitanExe | Out-Host
    }
    if ($LASTEXITCODE -ne 0) { throw "TITAN failed in $Dir" }
  } finally {
    Pop-Location
  }
}

function Assert-Exists([string]$Dir, [string]$File) {
  if (-not (Test-Path -LiteralPath (Join-Path $Dir $File))) {
    throw "missing expected file $File in $Dir"
  }
}

function Assert-NotExists([string]$Dir, [string]$File) {
  if (Test-Path -LiteralPath (Join-Path $Dir $File)) {
    throw "unexpected file $File in $Dir"
  }
}

function Assert-Contains([string]$Dir, [string]$File, [string]$Pattern) {
  $path = Join-Path $Dir $File
  Assert-Exists $Dir $File
  if (-not (Select-String -Path $path -Pattern $Pattern -Quiet)) {
    throw "$File did not contain expected pattern $Pattern"
  }
}

$inputPdb = Join-Path $RepoRoot "example\2\md1_short.pdb"
$lib = Join-Path $RepoRoot "lib_ff"

$blockDir = New-WorkDir "block"
Write-ArgFile $blockDir @(
  "-f=$inputPdb",
  "-1pdb=true",
  "-nocal=true",
  "-path-a=$lib",
  "-path-n=$lib",
  "-cmap=true",
  "-cmap-min-seqsep=0",
  "-qnative=true",
  "-qnative-min-seqsep=0",
  "-dihedral=true",
  "-dihedral-kind=both",
  "-uncertainty=block"
)
Run-Titan $blockDir
foreach ($file in @("uncertainty_summary.dat", "cmap_uncertainty.csv", "qnative_uncertainty.csv", "dihedral_uncertainty.csv")) {
  Assert-Exists $blockDir $file
}
if (-not (Select-String -Path (Join-Path $blockDir "uncertainty_summary.dat") -Pattern "uncertainty_mode: block" -Quiet)) {
  throw "summary did not record block mode"
}
if (-not (Select-String -Path (Join-Path $blockDir "cmap_uncertainty.csv") -Pattern "contact_count" -Quiet)) {
  throw "cmap sidecar missing contact_count row"
}
Assert-Contains $blockDir "uncertainty_summary.dat" "implementation_version: 2"
Assert-Contains $blockDir "uncertainty_summary.dat" "sidecar_status_table:"

$msdDir = New-WorkDir "msd"
$singleFrame = [System.IO.File]::ReadAllText($inputPdb)
$multiFrame = Join-Path $msdDir "md1_short_two_frames.pdb"
$frameTwo = $singleFrame.Replace("MODEL        1", "MODEL        2").Replace(
  "t= 100000.00000", "t= 100001.00000")
[System.IO.File]::WriteAllText(
  $multiFrame, $singleFrame + "`r`n" + $frameTwo, [System.Text.Encoding]::ASCII)
Write-ArgFile $msdDir @(
  "-f=$multiFrame",
  "-1pdb=true",
  "-nocal=true",
  "-path-a=$lib",
  "-path-n=$lib",
  "-msd=true",
  "-uncertainty=block"
)
Run-Titan $msdDir
Assert-Contains $msdDir "msd_uncertainty.csv" "D_A2_per_ps"
Assert-Contains $msdDir "msd_uncertainty.csv" "alpha"

$s2Dir = New-WorkDir "s2"
Write-ArgFile $s2Dir @(
  "-f=$(Join-Path $RepoRoot 'test\fixtures\s2_synthetic.pdb')",
  "-n=$(Join-Path $RepoRoot 'test\fixtures\s2_synthetic.ndx')",
  "-1pdb=true",
  "-nocal=true",
  "-path-a=$lib",
  "-path-n=$lib",
  "-s2=true",
  "-uncertainty=block"
)
Run-Titan $s2Dir
Assert-Contains $s2Dir "s2_uncertainty.csv" "plateau_C_raw"
Assert-Contains $s2Dir "s2_uncertainty.csv" "S2_clipped"

$membraneDir = New-WorkDir "membrane"
Write-ArgFile $membraneDir @(
  "-f=$(Join-Path $RepoRoot 'test\fixtures\membrane_synthetic.pdb')",
  "-n=$(Join-Path $RepoRoot 'test\fixtures\membrane_synthetic.ndx')",
  "-1pdb=true",
  "-nocal=true",
  "-path-a=$lib",
  "-path-n=$lib",
  "-membrane=true",
  "-membrane-group=auto",
  "-membrane-normal=z",
  "-membrane-grid=2,2",
  "-uncertainty=block"
)
Run-Titan $membraneDir
Assert-Contains $membraneDir "membrane_uncertainty.csv" "mean_apl_A2"
Assert-Contains $membraneDir "membrane_uncertainty.csv" "mean_thickness_A"
Assert-Contains $membraneDir "membrane_uncertainty.csv" "interdigitation"

$bootDirA = New-WorkDir "bootstrap_a"
Write-ArgFile $bootDirA @(
  "-f=$inputPdb",
  "-1pdb=true",
  "-nocal=true",
  "-path-a=$lib",
  "-path-n=$lib",
  "-qnative=true",
  "-qnative-min-seqsep=0",
  "-uncertainty=bootstrap",
  "-uncertainty-bootstrap=50",
  "-uncertainty-seed=7"
)
Run-Titan $bootDirA
$bootDirB = New-WorkDir "bootstrap_b"
Copy-Item -LiteralPath (Join-Path $bootDirA "arg.in") -Destination (Join-Path $bootDirB "arg.in")
Run-Titan $bootDirB
$a = Get-Content -Raw -Path (Join-Path $bootDirA "qnative_uncertainty.csv")
$b = Get-Content -Raw -Path (Join-Path $bootDirB "qnative_uncertainty.csv")
if ($a -ne $b) { throw "bootstrap uncertainty output is not deterministic" }

$edecompDir = New-WorkDir "edecomp"
Write-ArgFile $edecompDir @(
  "-f=$inputPdb",
  "-1pdb=true",
  "-path-a=$lib",
  "-path-n=$lib",
  "-g=1-2",
  "-edecomp=true",
  "-edecomp-source=mm",
  "-edecomp-level=both",
  "-edecomp-edge-top=5",
  "-uncertainty=block"
)
Run-Titan $edecompDir
Assert-Exists $edecompDir "edecomp_mm_residue_ranking.csv"
Assert-Exists $edecompDir "edecomp_uncertainty.csv"
if (-not (Select-String -Path (Join-Path $edecompDir "edecomp_uncertainty.csv") -Pattern "abs_mean" -Quiet)) {
  throw "EDecomp uncertainty missing abs_mean row"
}

$offDir = New-WorkDir "off"
Write-ArgFile $offDir @(
  "-f=$inputPdb",
  "-1pdb=true",
  "-nocal=true",
  "-path-a=$lib",
  "-path-n=$lib",
  "-qnative=true",
  "-qnative-min-seqsep=0"
)
Run-Titan $offDir
Assert-NotExists $offDir "uncertainty_summary.dat"
Assert-NotExists $offDir "qnative_uncertainty.csv"

Write-Host "PASS uncertainty_standard_observables_smoke"
