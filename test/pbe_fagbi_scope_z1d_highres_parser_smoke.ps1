param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z1d_parser'),
    [string]$LibFfRoot = $env:TITAN_LIB_FF
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if ([string]::IsNullOrWhiteSpace($LibFfRoot)) { $LibFfRoot = Join-Path $RepoRoot 'lib_ff' }
if (-not (Test-Path -LiteralPath $LibFfRoot)) {
    throw "Missing lib_ff directory: $LibFfRoot. Set -LibFfRoot or TITAN_LIB_FF for linked worktrees."
}
$LibFfRoot = (Resolve-Path -LiteralPath $LibFfRoot).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope Z1D parser work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) { Remove-Item -LiteralPath $workRootFull -Recurse -Force }
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Write-Z1dParserCase(
    [string]$Name,
    [string[]]$ExtraLines
) {
    $runDir = Join-Path $workRootFull (Get-FmmCliSafeName $Name)
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $pdbPath = Join-Path $runDir 'driver.pdb'
    $pqrPath = Join-Path $runDir 'input.pqr'
    Write-FmmCliDummyPdb $pdbPath
    Write-FagbiBornPqr -Path $pqrPath -Charge 1.0 -RadiusA 3.0
    $libPath = Convert-FmmCliToTitanPath ($LibFfRoot.TrimEnd('\') + '\')
    $baseLines = @(
        "-f=$(Convert-FmmCliToTitanPath $pdbPath)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-FmmCliToTitanPath $pqrPath)",
        '-pbe-tabi-sum=direct',
        '-pbe-tabi-fmm-m2l=off',
        '-pbe-tabi-surface-points=5121',
        '-pbe-tabi-sphere-mesh=icosphere',
        '-nt=1'
    )
    $argPath = Join-Path $runDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value ($baseLines + $ExtraLines) -Encoding ASCII
    return [pscustomobject]@{ run_dir=$runDir; arg=$argPath }
}

function Invoke-Z1dParserReject(
    [string]$Name,
    [string[]]$ExtraLines,
    [string]$ExpectedMessage
) {
    $case = Write-Z1dParserCase -Name $Name -ExtraLines $ExtraLines
    $stdoutPath = Join-Path $case.run_dir 'stdout.log'
    $stderrPath = Join-Path $case.run_dir 'stderr.log'
    $process = Start-Process `
        -FilePath $TitanExe `
        -ArgumentList @('-a','arg.in') `
        -WorkingDirectory $case.run_dir `
        -RedirectStandardOutput $stdoutPath `
        -RedirectStandardError $stderrPath `
        -NoNewWindow `
        -Wait `
        -PassThru
    $exitCode = $process.ExitCode
    $output = ''
    if (Test-Path -LiteralPath $stdoutPath) { $output += Get-Content -LiteralPath $stdoutPath -Raw }
    if (Test-Path -LiteralPath $stderrPath) { $output += Get-Content -LiteralPath $stderrPath -Raw }
    if ($exitCode -eq 0) {
        throw "$Name expected rejection but exited 0."
    }
    if ($output -notmatch [regex]::Escape($ExpectedMessage)) {
        throw "$Name missing expected message '$ExpectedMessage'. Output: $output"
    }
}

Invoke-Z1dParserReject `
    -Name 'no_highres_rejects_5121' `
    -ExtraLines @() `
    -ExpectedMessage 'must be between 8 and 4096'

Invoke-Z1dParserReject `
    -Name 'grid_highres_rejects' `
    -ExtraLines @('-pbe=grid','-pbe-tabi-highres-surface=true') `
    -ExpectedMessage 'restricted to direct TABI -pbe-test-pqr benchmark runs'

Invoke-Z1dParserReject `
    -Name 'treecode_highres_rejects' `
    -ExtraLines @('-pbe-tabi-sum=treecode','-pbe-tabi-highres-surface=true') `
    -ExpectedMessage 'restricted to direct TABI -pbe-test-pqr benchmark runs'

Invoke-Z1dParserReject `
    -Name 'highres_over_max_rejects' `
    -ExtraLines @('-pbe-tabi-highres-surface=true','-pbe-tabi-surface-points=65537') `
    -ExpectedMessage 'must be between 8 and 65536'

Write-Host 'pbe_fagbi_scope_z1d_highres_parser_smoke ok'
