param(
  [string]$TitanExe = ".\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\membrane_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null
$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$trajectoryOutputs = @(
  "membrane_leaflets.csv",
  "membrane_apl_timeseries.csv",
  "membrane_apl_per_lipid.csv",
  "membrane_thickness_timeseries.csv",
  "membrane_thickness_map.csv",
  "membrane_summary.dat",
  "membrane_order_parameters.csv",
  "membrane_interdigitation.dat",
  "membrane_flipflop_events.csv"
)
$aggregateOutputs = @(
  "membrane_summary_average.dat",
  "membrane_thickness_map_average.csv",
  "membrane_order_parameters_average.csv"
)
$outputs = $trajectoryOutputs + $aggregateOutputs

function Clear-MembraneOutputs {
  foreach ($file in $outputs) {
    if (Test-Path -LiteralPath (Join-Path $repo $file)) {
      Remove-Item -LiteralPath (Join-Path $repo $file) -Force
    }
  }
  Get-ChildItem -LiteralPath $repo -Filter "membrane_*_[0-9]*.*" -ErrorAction SilentlyContinue |
    Remove-Item -Force
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-MembraneOutputs
  & $titanPath -a $arg
  if ($LASTEXITCODE -ne 0) {
    throw "${Name}: TITAN exited with code $LASTEXITCODE"
  }
}

$base = @"
-f=..\..\test\fixtures\membrane_synthetic.pdb
-n=..\..\test\fixtures\membrane_synthetic.ndx
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
"@

Push-Location $repo
try {
  Run-TitanCase "default" ($base + "`r`n" + @"
-membrane=true
-membrane-group=auto
-membrane-normal=z
-membrane-grid=2,2
"@)
  foreach ($file in $trajectoryOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "default: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "membrane_summary.dat" -Pattern "reference_note:" -Quiet) -or
      -not (Select-String -Path "membrane_summary.dat" -Pattern "lipids: 8" -Quiet) -or
      -not (Select-String -Path "membrane_apl_timeseries.csv" -Pattern "mean_apl_A2" -Quiet)) {
    throw "default: membrane output does not match geometry contract"
  }

  Run-TitanCase "disabled" $base
  foreach ($file in $outputs) {
    if (Test-Path -LiteralPath $file) {
      throw "disabled: membrane output generated without -membrane=true: $file"
    }
  }

  $listFile = Join-Path $work "same_pdb_twice.txt"
  $fixture = (Resolve-Path -LiteralPath (Join-Path $repo "test\fixtures\membrane_synthetic.pdb")).Path
  [System.IO.File]::WriteAllText(
    $listFile, $fixture + "`r`n" + $fixture + "`r`n",
    [System.Text.Encoding]::ASCII)
  $multi = @"
-f=..\..\test\fixtures\membrane_synthetic.pdb
-i=same_pdb_twice.txt
-n=..\..\test\fixtures\membrane_synthetic.ndx
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
-membrane=true
-membrane-group=auto
-membrane-normal=z
-membrane-grid=2,2
"@
  Run-TitanCase "multi" $multi
  foreach ($file in @(
      "membrane_leaflets_0.csv",
      "membrane_leaflets_1.csv",
      "membrane_summary_0.dat",
      "membrane_summary_1.dat",
      "membrane_summary_average.dat",
      "membrane_thickness_map_average.csv",
      "membrane_order_parameters_average.csv")) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "multi: missing indexed or aggregate output $file"
    }
  }
  if (-not (Select-String -Path "membrane_summary_average.dat" -Pattern "replicates: 2" -Quiet) -or
      -not (Select-String -Path "membrane_thickness_map_average.csv" -Pattern "sem_thickness_A" -Quiet)) {
    throw "multi: aggregate outputs do not record replicate SEM"
  }

  Write-Host "PASS membrane_smoke"
}
finally {
  Pop-Location
}
