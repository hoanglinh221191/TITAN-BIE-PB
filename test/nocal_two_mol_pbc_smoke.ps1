param(
  [string]$TitanExe = ".\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\nocal_two_mol_pbc_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null
$output = Join-Path $repo "nocal_two_mol_pbc_output.pdb"
$arg = Join-Path $work "arg.in"
$index = Join-Path $repo "example\2\index.ndx"
$singleFramePdb = Join-Path $repo "example\2\md1_short.pdb"
$multiFramePdb = Join-Path $work "md1_short_two_frames.pdb"
$frameOne = [System.IO.File]::ReadAllText($singleFramePdb)
$frameTwo = $frameOne.Replace("MODEL        1", "MODEL        2").Replace(
  "t= 100000.00000", "t= 100001.00000")
[System.IO.File]::WriteAllText(
  $multiFramePdb, $frameOne + "`r`n" + $frameTwo, [System.Text.Encoding]::ASCII)

$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

function Read-NdxGroupSizes([string]$Path) {
  $sizes = [System.Collections.Generic.List[int]]::new()
  $active = -1
  foreach ($line in [System.IO.File]::ReadAllLines($Path)) {
    if ($line -match "^\s*\[") {
      $sizes.Add(0)
      $active++
      continue
    }
    if ($active -ge 0) {
      $sizes[$active] += @($line -split "\s+" | Where-Object { $_ }).Count
    }
  }
  return $sizes
}

$argsText = @"
-f=$multiFramePdb
-n=..\..\example\2\index.ndx
-g=3-4
-1pdb=true
-nocal=true
-opdb=nocal_two_mol_pbc_output.pdb
-path-a=..\..\example\2\lib_ff
-path-n=..\..\example\2\lib_ff
"@

[System.IO.File]::WriteAllText($arg, $argsText, [System.Text.Encoding]::ASCII)
if (Test-Path -LiteralPath $output) {
  Remove-Item -LiteralPath $output -Force
}

Push-Location $repo
try {
  & $titanPath -a $arg
  if ($LASTEXITCODE -ne 0) {
    throw "TITAN exited with code $LASTEXITCODE"
  }
  if (-not (Test-Path -LiteralPath $output)) {
    throw "missing PBC-only output PDB"
  }

  $sizes = Read-NdxGroupSizes $index
  $expectedAtoms = $sizes[3] + $sizes[4]
  $actualAtoms = (Select-String -Path $output -Pattern "^ATOM").Count
  if ($actualAtoms -ne $expectedAtoms) {
    throw "output has $actualAtoms atoms; expected mol1+mol2=$expectedAtoms"
  }
  if ((Select-String -Path $output -Pattern "^TER$").Count -ne 2) {
    throw "output PDB does not contain separate mol1 and mol2 records"
  }

  Write-Host "PASS nocal_two_mol_pbc_smoke"
}
finally {
  Pop-Location
  if (Test-Path -LiteralPath $output) {
    Remove-Item -LiteralPath $output -Force
  }
}
