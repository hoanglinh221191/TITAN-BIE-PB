param(
  [string]$TitanExe = ".\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\dry_run_inputs_local"
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
  "cal.txt",
  "cal.xvg"
)

function Clear-DryRunOutputs {
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
  Clear-DryRunOutputs
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

function Expect-TitanReject([string]$Name, [string]$ArgText, [string]$ExpectedMessage) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  $stdout = Join-Path $work ("stdout_" + $Name + ".txt")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-DryRunOutputs
  Push-Location $repo
  $savedErrorActionPreference = $ErrorActionPreference
  try {
    $ErrorActionPreference = "Continue"
    & $titanPath -a $arg *> $stdout
    $exitCode = $LASTEXITCODE
  }
  finally {
    $ErrorActionPreference = $savedErrorActionPreference
    Pop-Location
  }
  if ($exitCode -eq 0) {
    throw "${Name}: TITAN unexpectedly accepted invalid arguments"
  }
  if (-not (Select-String -Path $stdout -Pattern $ExpectedMessage -Quiet)) {
    throw "${Name}: expected validation message was not emitted"
  }
}

$base = @"
-f=..\..\example\2\md1_short.pdb
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
"@

Run-TitanCase "dry_manifest_schema" ($base + "`r`n" + @"
-dry-run=true
-schema-export=json
-cmap=true
-cmap-min-seqsep=0
"@)

Push-Location $repo
try {
  foreach ($file in @("titan_outputs_manifest.csv", "titan_run_config_summary.dat", "titan_schema.json")) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "dry_manifest_schema: missing expected dry-run output $file"
    }
  }
  foreach ($file in @("contact_probability_matrix.csv", "contact_probability_edges.csv", "contact_probability_summary.dat", "cal.txt", "cal.xvg")) {
    if (Test-Path -LiteralPath $file) {
      throw "dry_manifest_schema: dry-run unexpectedly wrote runtime output $file"
    }
  }
  if (-not (Select-String -Path "titan_run_config_summary.dat" -Pattern "dry_run: true" -Quiet)) {
    throw "dry_manifest_schema: summary does not mark dry_run"
  }
  if (-not (Select-String -Path "titan_outputs_manifest.csv" -Pattern "^CMap,contact_probability_matrix.csv" -Quiet)) {
    throw "dry_manifest_schema: manifest lacks expected CMap output row"
  }
  if (-not (Select-String -Path "titan_schema.json" -Pattern '"name": "-dry-run"' -Quiet)) {
    throw "dry_manifest_schema: schema lacks dry-run option"
  }
}
finally {
  Pop-Location
}

Expect-TitanReject "bad_schema" ($base + "`r`n-schema-export=bad") "schema-export must be one of"

Write-Host "PASS dry_run_validation_smoke"
