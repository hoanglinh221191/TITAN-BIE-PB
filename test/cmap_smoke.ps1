param(
  [string]$TitanExe = ".\test\bin\native\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\cmap_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null

$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$arg = Join-Path $work "arg_cmap.in"
$argText = @"
-f=..\..\example\2\md1_short.pdb
-1pdb
-nocal
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
-cmap
-cmap-mode=sidechain-heavy
-cmap-cut=4.5
-cmap-domain=intra
-cmap-min-seqsep=3
"@
[System.IO.File]::WriteAllText($arg, $argText, [System.Text.Encoding]::ASCII)

$expectedOutputs = @(
  "contact_probability_matrix.csv",
  "contact_probability_edges.csv",
  "contact_probability_summary.dat"
)

Push-Location $repo
try {
  foreach ($file in $expectedOutputs) {
    if (Test-Path -LiteralPath $file) {
      Remove-Item -LiteralPath $file -Force
    }
  }

  & $titanPath -a $arg
  if ($LASTEXITCODE -ne 0) {
    throw "TITAN exited with code $LASTEXITCODE"
  }

  foreach ($file in $expectedOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "Missing expected output $file"
    }
  }

  if (-not (Select-String -Path "contact_probability_summary.dat" -Pattern "module: titan_cmap" -Quiet)) {
    throw "Summary does not contain module: titan_cmap"
  }
  if (-not (Select-String -Path "contact_probability_summary.dat" -Pattern "mode: sidechain-heavy" -Quiet)) {
    throw "Summary does not contain mode: sidechain-heavy"
  }

  $edgesHeader = Get-Content -LiteralPath "contact_probability_edges.csv" -TotalCount 1
  if ($edgesHeader -notmatch "(^|,)probability(,|$)") {
    throw "Edges header does not contain probability column"
  }

  Write-Host "PASS cmap_smoke"
}
finally {
  Pop-Location
}
