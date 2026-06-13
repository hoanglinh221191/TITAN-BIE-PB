param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_cli_helpers_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$repoRoot = Split-Path -Parent $PSScriptRoot
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$fixture = Resolve-FmmCliFixture -Fixture 'research-alt' -RepoRoot $repoRoot
if (-not $fixture.EndsWith('example\3\rec_00.pqr')) {
    throw "Unexpected research-alt fixture: $fixture"
}

$sampleA = Join-Path $WorkRoot 'sample_a.pqr'
$sampleB = Join-Path $WorkRoot 'sample_b.pqr'
$infoA = Write-FmmCliSampledPqr -InputPath $fixture -OutputPath $sampleA -Cap 8
$infoB = Write-FmmCliSampledPqr -InputPath $fixture -OutputPath $sampleB -Cap 8
if ($infoA.sampled_atoms -ne 8 -or $infoB.sampled_atoms -ne 8) {
    throw 'Deterministic sampler did not respect the atom cap.'
}
if ((Get-FileHash -Algorithm SHA256 $sampleA).Hash -ne
    (Get-FileHash -Algorithm SHA256 $sampleB).Hash) {
    throw 'Deterministic sampler produced different outputs.'
}

$runDir = Join-Path $WorkRoot 'args'
New-Item -ItemType Directory -Path $runDir -Force | Out-Null
$argPath = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $repoRoot `
    -PqrPath $sampleA -Mode treecode -FmmMode single-matvec-analytic `
    -SurfacePoints 8 -MaxIter 20 -Restart 8 -Tolerance 1.0e-5 `
    -Threads 1 -SaltM 0.15 -InternalDielectric 1.0 `
    -ExternalDielectric 80.0 -TemperatureK 298.15 -AuditSampleCount 128
$badOptions = @(Get-Content -LiteralPath $argPath | Where-Object {
    -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
})
if ($badOptions.Count -ne 0) {
    throw "Generated arg.in contains non equals-syntax options: $($badOptions -join '; ')"
}
foreach ($expected in @(
    '-pbe-salt=0.15',
    '-pbe-int-diel=1',
    '-pbe-ext-diel=80',
    '-pbe-temp=298.15',
    '-pbe-tabi-opening-theta=0.8',
    '-pbe-tabi-matvec-audit=sampled',
    '-pbe-tabi-matvec-audit-sample-count=128')) {
    if (-not (Select-String -LiteralPath $argPath -SimpleMatch $expected -Quiet)) {
        throw "Generated arg.in missing $expected"
    }
}

$gridRunDir = Join-Path $WorkRoot 'grid_args'
New-Item -ItemType Directory -Path $gridRunDir -Force | Out-Null
$gridArgPath = Write-FmmCliArgFile -RunDir $gridRunDir -RepoRoot $repoRoot `
    -PqrPath $sampleA -Backend grid -Mode direct -FmmMode off `
    -SurfacePoints 8 -MaxIter 20 -Restart 8 -Tolerance 1.0e-5 `
    -Threads 1 -SaltM 0.15 -InternalDielectric 1.0 `
    -ExternalDielectric 80.0 -TemperatureK 298.15
if (-not (Select-String -LiteralPath $gridArgPath -SimpleMatch '-pbe=grid' -Quiet)) {
    throw 'Generated grid arg.in missing -pbe=grid'
}
if (Select-String -LiteralPath $gridArgPath -SimpleMatch '-pbe-tabi-fmm-m2l=' -Quiet) {
    throw 'Generated grid arg.in unexpectedly contains TABI FMM/M2L option'
}

$timeout = Get-FmmCliFailureClass -DirectExit 124 -FmmExit 0 `
    -Direct $null -Fmm $null -Fixture research-medium -RequireSpeed $false
if ($timeout -ne 'timeout_not_verified') {
    throw "Timeout classification drifted: $timeout"
}

Write-Host 'PASS pbe_fmm_m2l_cli_helpers_smoke'
