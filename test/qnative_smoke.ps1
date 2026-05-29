param(
  [string]$TitanExe = ".\test\bin\native\titan.exe"
)

$ErrorActionPreference = "Stop"
$env:TITAN_NO_PAUSE = "1"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$work = Join-Path $repo "test\qnative_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null
$singleFramePdb = Join-Path $repo "example\2\md1_short.pdb"
$singleFramePdbArg = "example\2\md1_short.pdb"
$multiFramePdb = Join-Path $work "md1_short_two_frames.pdb"
$shiftedRefPdb = Join-Path $work "md1_short_ref_shifted.pdb"
$frameOne = [System.IO.File]::ReadAllText($singleFramePdb)
$frameTwo = $frameOne.Replace("MODEL        1", "MODEL        2").Replace(
  "t= 100000.00000", "t= 100001.00000")
[System.IO.File]::WriteAllText(
  $multiFramePdb, $frameOne + "`r`n" + $frameTwo, [System.Text.Encoding]::ASCII)

$culture = [System.Globalization.CultureInfo]::InvariantCulture
$shiftedLines = New-Object System.Collections.Generic.List[string]
$shiftedAtomCount = 0
foreach ($line in [System.IO.File]::ReadLines($singleFramePdb)) {
  if (($line.StartsWith("ATOM  ") -or $line.StartsWith("HETATM")) -and
      $line.Length -ge 38 -and $shiftedAtomCount -lt 1000) {
    $x = [double]::Parse($line.Substring(30, 8).Trim(), $culture) + 100.0
    $xText = $x.ToString("F3", $culture).PadLeft(8)
    $shiftedLines.Add($line.Substring(0, 30) + $xText + $line.Substring(38))
    $shiftedAtomCount += 1
  } else {
    $shiftedLines.Add($line)
  }
}
[System.IO.File]::WriteAllLines($shiftedRefPdb, $shiftedLines, [System.Text.Encoding]::ASCII)

$titanPath = $TitanExe
if (-not [System.IO.Path]::IsPathRooted($titanPath)) {
  $titanPath = Join-Path $repo $titanPath
}
$titanPath = (Resolve-Path -LiteralPath $titanPath).Path

$outputs = @(
  "qnative_timeseries.csv",
  "qnative_native_contacts.csv",
  "qnative_contact_occupancy.csv",
  "qnative_summary.dat"
)
$averageOutputs = @(
  "qnative_timeseries_average.csv",
  "qnative_native_contacts_average.csv",
  "qnative_contact_occupancy_average.csv",
  "qnative_summary_average.dat"
)

function Clear-QNativeOutputs {
  foreach ($file in ($outputs + $averageOutputs)) {
    $path = Join-Path $repo $file
    if (Test-Path -LiteralPath $path) {
      Remove-Item -LiteralPath $path -Force
    }
  }
  Get-ChildItem -LiteralPath $repo -Filter "qnative_*_[0-9]*.*" -ErrorAction SilentlyContinue |
    Remove-Item -Force
}

function Run-TitanCase([string]$Name, [string]$ArgText) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-QNativeOutputs
  & $titanPath -a $arg
  if ($LASTEXITCODE -ne 0) {
    throw "${Name}: TITAN exited with code $LASTEXITCODE"
  }
}

function Expect-TitanReject([string]$Name, [string]$ArgText, [string]$ExpectedMessage) {
  $arg = Join-Path $work ("arg_" + $Name + ".in")
  $stdout = Join-Path $work ("stdout_" + $Name + ".txt")
  [System.IO.File]::WriteAllText($arg, $ArgText, [System.Text.Encoding]::ASCII)
  Clear-QNativeOutputs
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

function Get-QNativeTotalNativeContacts([string]$Path) {
  $line = Select-String -Path $Path -Pattern "^total_native_contacts:\s+([0-9]+)" |
    Select-Object -First 1
  if ($null -eq $line -or $line.Line -notmatch "^total_native_contacts:\s+([0-9]+)") {
    throw "could not parse total_native_contacts from $Path"
  }
  return [int]$Matches[1]
}

$base = @"
-f=..\..\example\2\md1_short.pdb
-1pdb
-nocal
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
"@

Push-Location $repo
try {
  Run-TitanCase "default" ($base + "`r`n" + @"
-qnative
-qnative-min-seqsep=0
"@)
  foreach ($file in $outputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "default: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "qnative_summary.dat" -Pattern "mode: sidechain-com" -Quiet) -or
      -not (Select-String -Path "qnative_summary.dat" -Pattern "cutoff_A: 6.500000" -Quiet) -or
      -not (Select-String -Path "qnative_timeseries.csv" -Pattern "Q,formed_native_contacts,total_native_contacts" -Quiet)) {
    throw "default: qnative summary/timeseries does not match sidechain-com contract"
  }
  $defaultNativeContacts = Get-QNativeTotalNativeContacts "qnative_summary.dat"

  Run-TitanCase "explicit_ref" ($base + "`r`n" + @"
-qnative
-qnative-ref=$shiftedRefPdb
-qnative-min-seqsep=0
"@)
  foreach ($file in $outputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "explicit_ref: missing expected output $file"
    }
  }
  if (-not (Select-String -Path "qnative_summary.dat" -Pattern ([regex]::Escape("reference_source: $shiftedRefPdb")) -Quiet)) {
    throw "explicit_ref: qnative summary did not record -qnative-ref source"
  }
  $explicitRefNativeContacts = Get-QNativeTotalNativeContacts "qnative_summary.dat"
  if ($explicitRefNativeContacts -eq $defaultNativeContacts) {
    throw "explicit_ref: shifted -qnative-ref did not change native-contact count"
  }

  $interface = @"
-f=$multiFramePdb
-1pdb
-nocal
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
-n=..\..\example\2\index.ndx
-g=3-4
-qnative
-qnative-scope=interface
"@
  Run-TitanCase "interface" $interface
  if (-not (Select-String -Path "qnative_summary.dat" -Pattern "scope: interface" -Quiet) -or
      -not (Select-String -Path "qnative_native_contacts.csv" -Pattern "mol2\|" -Quiet)) {
    throw "interface: qnative did not use mol1/mol2 native-contact roles"
  }

  Expect-TitanReject "interface_domain_rejected" ($interface + "`r`n-qnative-domain=both") "QNative: -qnative-domain and -qnative-min-seqsep apply only to -qnative-scope=mol1."
  Expect-TitanReject "beta_without_soft_rejected" ($base + "`r`n-qnative-beta=2.0") "QNative: -qnative-beta is valid only with -qnative-soft=true."
  Expect-TitanReject "bad_tpt_qnative_range_rejected" ($base + "`r`n-tpt`r`n-tpt-source=qnative_range:0.8,0.2`r`n-tpt-sink=micro:0") "TPT: -tpt-source qnative_range minimum must be <= maximum."

  $listFile = Join-Path $work "same_pdb_twice.txt"
  [System.IO.File]::WriteAllText(
    $listFile, $multiFramePdb + "`r`n" + $multiFramePdb + "`r`n",
    [System.Text.Encoding]::ASCII)
  $pooled = @"
-f=$multiFramePdb
-i=$listFile
-1pdb
-nocal
-path-a=..\..\lib_ff
-path-n=..\..\lib_ff
-qnative
-qnative-min-seqsep=0
"@
  Run-TitanCase "pooled" $pooled
  foreach ($file in @("qnative_timeseries_0.csv", "qnative_timeseries_1.csv") + $averageOutputs) {
    if (-not (Test-Path -LiteralPath $file)) {
      throw "pooled: missing expected output $file"
    }
  }

  Run-TitanCase "disabled" $base
  foreach ($file in $outputs) {
    if (Test-Path -LiteralPath $file) {
      throw "disabled: qnative output generated without -qnative: $file"
    }
  }

  Write-Host "PASS qnative_smoke"
}
finally {
  Pop-Location
}
