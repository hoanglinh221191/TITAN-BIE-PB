param(
  [string]$TitanExe = ".\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\output_manifest_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null

$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$outputs = @(
  "titan_outputs_manifest.csv",
  "titan_run_config_summary.dat",
  "titan_schema.json",
  "contact_probability_matrix.csv",
  "contact_probability_edges.csv",
  "contact_probability_summary.dat",
  "contact_count_timeseries.csv",
  "contact_residue_occupancy.csv"
)

function Clear-ManifestSmokeOutputs {
  foreach ($file in $outputs) {
    $path = Join-Path $repo $file
    if (Test-Path -LiteralPath $path) {
      Remove-Item -LiteralPath $path -Force
    }
  }
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-ManifestSmokeOutputs
  Push-Location $repo
  try {
    & $titanPath -a $arg
    if ($LASTEXITCODE -ne 0) {
      throw "${Name}: TITAN exited with code $LASTEXITCODE"
    }
  }
  finally {
    Pop-Location
  }
}

$base = @"
-f=..\..\example\2\md1_short.pdb
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
"@

Run-TitanCase "normal_manifest" ($base + "`r`n" + @"
-output-manifest=true
-cmap=true
-cmap-min-seqsep=0
"@)

Push-Location $repo
try {
  foreach ($file in @("titan_outputs_manifest.csv", "titan_run_config_summary.dat", "contact_probability_matrix.csv", "contact_probability_edges.csv", "contact_probability_summary.dat")) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "normal_manifest: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "titan_outputs_manifest.csv" -Pattern "^CMap,contact_probability_matrix.csv" -Quiet)) {
    throw "normal_manifest: manifest lacks CMap matrix row"
  }
  if (-not (Select-String -Path "titan_outputs_manifest.csv" -Pattern "^Core,titan_run_config_summary.dat" -Quiet)) {
    throw "normal_manifest: manifest lacks run config row"
  }
  if (-not (Select-String -Path "titan_run_config_summary.dat" -Pattern "active_modules: CMap" -Quiet)) {
    throw "normal_manifest: run config summary lacks active module"
  }
}
finally {
  Pop-Location
}

Write-Host "PASS output_manifest_smoke"
