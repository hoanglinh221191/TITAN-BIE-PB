param(
    [Parameter(Mandatory = $true)]
    [string]$ExampleDirectory,
    [string]$TitanExe,
    [string[]]$PqrFiles,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
$example = (Resolve-Path -LiteralPath $ExampleDirectory).Path
$pbeRoot = (Resolve-Path -LiteralPath (Join-Path $example '..')).Path
$repoRoot = (Resolve-Path -LiteralPath (Join-Path $pbeRoot '..\..')).Path
if ([string]::IsNullOrWhiteSpace($TitanExe)) {
    $TitanExe = Join-Path $repoRoot 'titan.exe'
}
if (-not [IO.Path]::IsPathRooted($TitanExe)) {
    $TitanExe = Join-Path $repoRoot $TitanExe
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$runDir = Join-Path $example 'run'
if (Test-Path -LiteralPath $runDir) {
    if (-not $Force) {
        throw "Run directory already exists: $runDir. Re-run with -Force to replace it."
    }
    Remove-Item -LiteralPath $runDir -Recurse -Force
}
New-Item -ItemType Directory -Path $runDir | Out-Null

Copy-Item -LiteralPath (Join-Path $example 'arg.in') -Destination $runDir
Copy-Item -LiteralPath (Join-Path $pbeRoot 'common\driver.pdb') -Destination $runDir
Copy-Item -LiteralPath (Join-Path $pbeRoot 'common\driver.top') -Destination $runDir
Copy-Item -LiteralPath (Join-Path $pbeRoot 'common\index.ndx') -Destination $runDir
Copy-Item -LiteralPath (Join-Path $pbeRoot 'common\cal.txt') -Destination $runDir
Copy-Item -LiteralPath (Join-Path $pbeRoot 'common\lib_ff') `
    -Destination (Join-Path $runDir 'lib_ff') -Recurse
foreach ($name in $PqrFiles) {
    Copy-Item -LiteralPath (Join-Path $example $name) -Destination $runDir
}

Push-Location $runDir
try {
    & $TitanExe -a arg.in
    if ($LASTEXITCODE -ne 0) {
        throw "TITAN exited with code $LASTEXITCODE."
    }
    if (-not (Test-Path -LiteralPath 'pbe_summary.dat')) {
        throw 'TITAN completed without pbe_summary.dat.'
    }
    Write-Host "Completed: $runDir"
} finally {
    Pop-Location
}

