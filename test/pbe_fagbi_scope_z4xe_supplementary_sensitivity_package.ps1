param(
    [string]$PythonExe = "",
    [switch]$AllowFailedDecision
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir
$builder = Join-Path $scriptDir "pbe_fagbi_scope_z4xe_build_supplement.py"
$matrix = Join-Path $repoRoot "test\tmp\fagbi_scope_z4w\fagbi_scope_z4w_matrix.csv"

if (-not (Test-Path -LiteralPath $matrix)) {
    throw "Missing Z4W matrix: $matrix. Run Z4W full before packaging Z4X-E."
}

if ([string]::IsNullOrWhiteSpace($PythonExe)) {
    $bundledPython = Join-Path $env:USERPROFILE ".cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
    if (Test-Path -LiteralPath $bundledPython) {
        $PythonExe = $bundledPython
    } else {
        $cmd = Get-Command python -ErrorAction SilentlyContinue
        if ($null -eq $cmd) {
            $cmd = Get-Command py -ErrorAction SilentlyContinue
        }
        if ($null -eq $cmd) {
            throw "Python was not found. Pass -PythonExe <path-to-python>."
        }
        $PythonExe = $cmd.Source
    }
}

& $PythonExe $builder
$exitCode = $LASTEXITCODE
if ($exitCode -ne 0) {
    if ($AllowFailedDecision) {
        Write-Warning "Z4X-E packaging exited with code $exitCode"
        exit 0
    }
    exit $exitCode
}
