param(
  [string]$TitanExe = ".\test\bin\native\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\msd_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null

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

$translationalOutputs = @(
  "msd_timeseries.dat",
  "msd_diffusion_summary.dat",
  "msd_alpha_fit.dat"
)
$rotationalOutputs = @(
  "rotational_msd.dat",
  "rotational_diffusion_summary.dat",
  "rotational_quaternion_timeseries.csv",
  "rotational_covariance_fit.csv"
)
$outputs = $translationalOutputs + $rotationalOutputs

function Clear-MsdOutputs {
  foreach ($file in $outputs) {
    $path = Join-Path $repo $file
    if (Test-Path -LiteralPath $path) {
      Remove-Item -LiteralPath $path -Force
    }
  }
  Get-ChildItem -LiteralPath $repo -Filter "msd_*_[0-9]*.dat" -ErrorAction SilentlyContinue |
    Remove-Item -Force
  Get-ChildItem -LiteralPath $repo -Filter "rotational_*_[0-9]*.dat" -ErrorAction SilentlyContinue |
    Remove-Item -Force
  Get-ChildItem -LiteralPath $repo -Filter "rotational_*_[0-9]*.csv" -ErrorAction SilentlyContinue |
    Remove-Item -Force
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-MsdOutputs
  & $titanPath -a $arg
  if ($LASTEXITCODE -ne 0) {
    throw "${Name}: TITAN exited with code $LASTEXITCODE"
  }
}

function Expect-TitanReject([string]$Name, [string]$ArgText, [string]$ExpectedMessage) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  $stdout = Join-Path $work ("stdout_" + $Name + ".txt")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-MsdOutputs
  $savedErrorActionPreference = $ErrorActionPreference
  try {
    $ErrorActionPreference = "Continue"
    & $titanPath -a $arg *> $stdout
    $exitCode = $LASTEXITCODE
  }
  finally {
    $ErrorActionPreference = $savedErrorActionPreference
  }
  if ($exitCode -eq 0) {
    throw "${Name}: TITAN unexpectedly accepted invalid arguments"
  }
  if (-not (Select-String -Path $stdout -Pattern $ExpectedMessage -Quiet)) {
    throw "${Name}: expected validation message was not emitted"
  }
}

$base = @"
-f=md1_short_two_frames.pdb
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
"@

Push-Location $repo
try {
  Run-TitanCase "default" ($base + "`r`n" + @"
-msd=true
"@)
  foreach ($file in $translationalOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "default: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "msd_timeseries.dat" -Pattern "lag_index lag_time_ps msd_A2 samples" -Quiet) -or
      -not (Select-String -Path "msd_diffusion_summary.dat" -Pattern "selection: mol1" -Quiet) -or
      -not (Select-String -Path "msd_diffusion_summary.dat" -Pattern "mode: com" -Quiet)) {
    throw "default: MSD output does not match translational COM contract"
  }

  Run-TitanCase "atom_mode" ($base + "`r`n" + @"
-msd=true
-msd-selection=atom:CA
-msd-mode=atom
-msd-fit-window=auto
"@)
  if (-not (Select-String -Path "msd_diffusion_summary.dat" -Pattern "mode: atom" -Quiet) -or
      -not (Select-String -Path "msd_diffusion_summary.dat" -Pattern "selection: atom:CA" -Quiet)) {
    throw "atom_mode: MSD summary did not record atom selection"
  }

  Run-TitanCase "rotational" ($base + "`r`n" + @"
-msd=true
-msd-rot=true
-msd-rot-fit-window=1,3
"@)
  foreach ($file in $rotationalOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "rotational: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "rotational_msd.dat" -Pattern "lag_index lag_time_ps rotational_msd_rad2 C1 C2 samples" -Quiet) -or
      -not (Select-String -Path "rotational_diffusion_summary.dat" -Pattern "reference_note:" -Quiet) -or
      -not (Select-String -Path "rotational_diffusion_summary.dat" -Pattern "fit_window_ps: 1.00000000,3.00000000" -Quiet)) {
    throw "rotational: output does not match rotational MSD contract"
  }

  $listFile = Join-Path $work "same_pdb_twice.txt"
  [System.IO.File]::WriteAllText(
    $listFile, $multiFramePdb + "`r`n" + $multiFramePdb + "`r`n",
    [System.Text.Encoding]::ASCII)
  $multi = @"
-f=md1_short_two_frames.pdb
-i=same_pdb_twice.txt
-1pdb=true
-nocal=true
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
-msd=true
-msd-rot=true
"@
  Run-TitanCase "multi" $multi
  foreach ($file in @("msd_timeseries_0.dat", "msd_timeseries_1.dat", "msd_diffusion_summary_0.dat", "msd_diffusion_summary_1.dat", "msd_alpha_fit_0.dat", "msd_alpha_fit_1.dat", "rotational_msd_0.dat", "rotational_msd_1.dat", "rotational_diffusion_summary_0.dat", "rotational_diffusion_summary_1.dat", "rotational_quaternion_timeseries_0.csv", "rotational_quaternion_timeseries_1.csv", "rotational_covariance_fit_0.csv", "rotational_covariance_fit_1.csv")) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "multi: missing indexed output $file"
    }
  }

  Run-TitanCase "disabled" $base
  foreach ($file in $outputs) {
    if (Test-Path -LiteralPath $file) {
      throw "disabled: MSD output generated without -msd=true: $file"
    }
  }

  Expect-TitanReject "invalid_rot_fit_window" ($base + "`r`n" + @"
-msd=true
-msd-rot=true
-msd-rot-fit-window=3,1
"@) "MSD: -msd-rot-fit-window must be auto or tmin,tmax with tmax > tmin."

  Write-Host "PASS msd_smoke"
}
finally {
  Pop-Location
}
