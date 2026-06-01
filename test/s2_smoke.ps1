param(
  [string]$TitanExe = ".\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\s2_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null
$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$trajectoryOutputs = @(
  "s2_order_parameters.csv",
  "s2_correlation_functions.dat",
  "s2_summary.dat",
  "s2_bfactor.pdb"
)
$aggregateOutputs = @(
  "s2_order_parameters_average.csv",
  "s2_summary_average.dat",
  "s2_bfactor_average.pdb"
)
$outputs = $trajectoryOutputs + $aggregateOutputs

function Clear-S2Outputs {
  foreach ($file in $outputs) {
    if (Test-Path -LiteralPath (Join-Path $repo $file)) {
      Remove-Item -LiteralPath (Join-Path $repo $file) -Force
    }
  }
  Get-ChildItem -LiteralPath $repo -Filter "s2_*_[0-9]*.*" -ErrorAction SilentlyContinue |
    Remove-Item -Force
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-S2Outputs
  & $titanPath -a $arg
  if ($LASTEXITCODE -ne 0) {
    throw "${Name}: TITAN exited with code $LASTEXITCODE"
  }
}

$base = @"
-f=..\..\test\fixtures\s2_synthetic.pdb
-n=..\..\test\fixtures\s2_synthetic.ndx
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
"@

Push-Location $repo
try {
  Run-TitanCase "default" ($base + "`r`n" + @"
-s2=true
-s2-bfactor-pdb=true
"@)
  foreach ($file in $trajectoryOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "default: missing expected output $file"
    }
  }

  Run-TitanCase "disabled" ($base + "`r`n-s2=false`r`n")
  foreach ($file in $outputs) {
    if (Test-Path -LiteralPath $file) {
      throw "disabled: S2 output generated with -s2=false: $file"
    }
  }

  Run-TitanCase "explicit" ($base + "`r`n" + @"
-s2=true
-s2-vector=nh
-s2-align=none
-s2-maxlag=4
-s2-fit-window=1,4
-s2-residues=A:13-16
-s2-bfactor-pdb=true
"@)
  foreach ($file in $trajectoryOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "explicit: missing expected output $file"
    }
  }

  $listFile = Join-Path $work "same_pdb_twice.txt"
  $fixture = (Resolve-Path -LiteralPath (Join-Path $repo "test\fixtures\s2_synthetic.pdb")).Path
  [System.IO.File]::WriteAllText(
    $listFile, $fixture + "`r`n" + $fixture + "`r`n",
    [System.Text.Encoding]::ASCII)
  $multi = @"
-f=..\..\test\fixtures\s2_synthetic.pdb
-i=same_pdb_twice.txt
-n=..\..\test\fixtures\s2_synthetic.ndx
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
-s2=true
-s2-bfactor-pdb=true
"@
  Run-TitanCase "multi" $multi
  foreach ($file in @(
      "s2_order_parameters_0.csv",
      "s2_order_parameters_1.csv",
      "s2_correlation_functions_0.dat",
      "s2_correlation_functions_1.dat",
      "s2_summary_0.dat",
      "s2_summary_1.dat",
      "s2_bfactor_0.pdb",
      "s2_bfactor_1.pdb",
      "s2_order_parameters_average.csv",
      "s2_summary_average.dat",
      "s2_bfactor_average.pdb")) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "multi: missing indexed or aggregate output $file"
    }
  }

  Write-Host "PASS s2_smoke"
}
finally {
  Pop-Location
}
