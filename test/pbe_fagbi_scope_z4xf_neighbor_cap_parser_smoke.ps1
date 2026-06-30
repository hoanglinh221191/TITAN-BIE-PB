param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4xf_parser')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$RuntimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
        $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
    }
    $RuntimeRepoRoot = Split-Path -Parent $commonGitDir
}
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)

$pqr = Join-Path $WorkRoot 'generated.pqr'
@(
    'ATOM      1  C   GLY A   1       0.000   0.000   0.000  1.0000  1.7000',
    'ATOM      2  O   GLY A   1       4.000   0.000   0.000 -1.0000  1.5000'
) | Set-Content -LiteralPath $pqr -Encoding ASCII

function New-Z4xfArg(
    [string]$Name,
    [int]$IluCap,
    [bool]$EnableSweep = $false
) {
    $dir = Join-Path $WorkRoot $Name
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    return Write-FmmCliArgFile -RunDir $dir -RepoRoot $RuntimeRepoRoot -PqrPath $pqr `
        -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
        -SurfacePoints 8 -MaxIter 40 -Restart 16 -Tolerance 1e-6 -Threads 1 `
        -SaltM 0.15 -InternalDielectric 1 -ExternalDielectric 80 `
        -TemperatureK 298.15 -OpeningTheta 0.6 -AuditSampleCount 128 `
        -AuditMode sampled-detail -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap $IluCap -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular -TaylorOrder 4 -FmmCartesianAudit off `
        -FmmLocal gmres-p4 -FmmLocalRefine operator-tail-static `
        -FmmLocalRefineTolerance 0.002 -FmmLocalCache green-jet-compressed `
        -FmmLocalCacheBudgetMb 256 -FmmLocalFallback packed-avx2 `
        -FmmLocalPlanBuilder metric-cache -FmmLocalRefineProbeCount 4 `
        -FmmLocalRefineProbeMode rademacher `
        -FmmLocalGeometryGuard strict-geometry-mac -FmmLocalGeometryTheta 0.5 `
        -FmmLocalIluCapSweep $EnableSweep
}

function Invoke-Z4xfParserCase(
    [string]$Name,
    [string[]]$Lines,
    [bool]$ExpectSuccess,
    [string]$ExpectedMessage = ''
) {
    $dir = Join-Path $WorkRoot $Name
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    $arg = Join-Path $dir 'arg.in'
    $Lines | Set-Content -LiteralPath $arg -Encoding ASCII
    $bare = @(Get-Content -LiteralPath $arg | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) {
        throw "$Name generated bare arg.in options: $($bare -join ', ')"
    }
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $dir -ArgPath $arg `
        -CaseTimeoutSeconds 180
    $log = (Get-Content -Raw (Join-Path $dir 'stdout.log')) +
        (Get-Content -Raw (Join-Path $dir 'stderr.log'))
    if ($ExpectSuccess -and $run.exit_code -ne 0) {
        throw "$Name unexpectedly failed: $log"
    }
    if (-not $ExpectSuccess -and $run.exit_code -eq 0) {
        throw "$Name unexpectedly succeeded"
    }
    if (-not $ExpectSuccess -and $log -notmatch $ExpectedMessage) {
        throw "$Name error was unclear: $log"
    }
}

$cap32NoSweep = @(Get-Content -LiteralPath (New-Z4xfArg 'cap32_no_sweep_seed' 32 $false))
Invoke-Z4xfParserCase 'cap32_no_sweep_rejects' $cap32NoSweep $false `
    'gmres-p4.*neighbor cap 64'

foreach ($cap in @(32, 64, 128)) {
    $arg = @(Get-Content -LiteralPath (New-Z4xfArg "cap${cap}_sweep_seed" $cap $true))
    Invoke-Z4xfParserCase "cap${cap}_sweep_accepts" $arg $true
}

$invalidMode = @(Get-Content -LiteralPath (New-Z4xfArg 'invalid_mode_seed' 32 $true)) -replace `
    '^-pbe-tabi-fmm-local=gmres-p4$',
    '-pbe-tabi-fmm-local=off'
Invoke-Z4xfParserCase 'sweep_flag_rejects_non_gmres' $invalidMode $false `
    'ilu-cap-sweep.*gmres-p4'

Write-Host 'PASS pbe_fagbi_scope_z4xf_neighbor_cap_parser_smoke'
