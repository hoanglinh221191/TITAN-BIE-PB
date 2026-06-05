param(
  [string]$TitanExe = ".\titan.exe",
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\performance_observer_smoke")
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\performance_observer_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

$unitExe = Join-Path $WorkRoot "titan_performance_observer_smoke.exe"
Push-Location $repo
try {
  & $Compiler -std=c++23 -O2 -Isrc `
    test\titan_performance_observer_smoke.cpp src\titan_performance_observer.cpp `
    -o $unitExe
  if ($LASTEXITCODE -ne 0) {
    throw "performance_observer_smoke: failed to compile C++ observer smoke"
  }
}
finally {
  Pop-Location
}
& $unitExe
if ($LASTEXITCODE -ne 0) {
  throw "performance_observer_smoke: C++ observer smoke failed"
}

$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$perfOutputs = @(
  "titan_performance_summary.dat",
  "titan_performance_modules.csv",
  "titan_performance_frames.csv"
)

function Clear-PerfOutputs {
  foreach ($file in $perfOutputs) {
    $path = Join-Path $repo $file
    if (Test-Path -LiteralPath $path) {
      Remove-Item -LiteralPath $path -Force
    }
  }
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-PerfOutputs
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

Run-TitanCase "summary" ($base + "`r`n" + @"
-cmap=true
-cmap-min-seqsep=0
-qnative=true
-qnative-min-seqsep=0
-dihedral=true
-perf-report=summary
"@)
Push-Location $repo
try {
  foreach ($file in @("titan_performance_summary.dat", "titan_performance_modules.csv")) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "summary: missing expected output $file"
    }
  }
  if (Test-Path -LiteralPath "titan_performance_frames.csv") {
    throw "summary: detail frames file should not be written"
  }
  foreach ($module in @("CMap", "QNative", "Dihedral")) {
    if (-not (Select-String -Path "titan_performance_modules.csv" -Pattern "^$module," -Quiet)) {
      throw "summary: missing module row $module"
    }
  }
  if (-not (Select-String -Path "titan_performance_summary.dat" -Pattern "compatibility_note" -Quiet)) {
    throw "summary: shared performance summary lacks compatibility note"
  }
}
finally {
  Pop-Location
}

Run-TitanCase "detail" ($base + "`r`n" + @"
-cmap=true
-cmap-min-seqsep=0
-qnative=true
-qnative-min-seqsep=0
-dihedral=true
-perf-report=detail
"@)
Push-Location $repo
try {
  foreach ($file in $perfOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "detail: missing expected output $file"
    }
  }
  foreach ($module in @("CMap", "QNative", "Dihedral")) {
    if (-not (Select-String -Path "titan_performance_frames.csv" -Pattern "^$module," -Quiet)) {
      throw "detail: missing frame row $module"
    }
  }
}
finally {
  Pop-Location
}

Run-TitanCase "off_default" ($base + "`r`n" + @"
-cmap=true
-cmap-min-seqsep=0
"@)
Push-Location $repo
try {
  foreach ($file in $perfOutputs) {
    if (Test-Path -LiteralPath $file) {
      throw "off_default: unexpected shared performance output $file"
    }
  }
}
finally {
  Pop-Location
}

Write-Host "PASS performance_observer_smoke"
