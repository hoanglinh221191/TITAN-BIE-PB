param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4n_parser')
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
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)

$pqr = Join-Path $WorkRoot 'generated.pqr'
@(
    'ATOM      1  C   GLY A   1       0.000   0.000   0.000  1.0000  1.7000',
    'ATOM      2  O   GLY A   1       4.000   0.000   0.000 -1.0000  1.5000'
) | Set-Content -LiteralPath $pqr -Encoding ASCII

function New-Z4nArgs([string]$Name) {
    $dir = Join-Path $WorkRoot $Name
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    return Write-FmmCliArgFile -RunDir $dir -RepoRoot $RuntimeRepoRoot -PqrPath $pqr `
        -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
        -SurfacePoints 8 -MaxIter 40 -Restart 16 -Tolerance 1e-6 -Threads 1 `
        -SaltM 0.15 -InternalDielectric 1 -ExternalDielectric 80 `
        -TemperatureK 298.15 -OpeningTheta 0.6 -AuditSampleCount 128 `
        -AuditMode sampled-detail -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular -TaylorOrder 4 -FmmCartesianAudit off `
        -FmmLocal gmres-p4 -FmmLocalRefine operator-tail-static `
        -FmmLocalRefineTolerance 0.002 -FmmLocalCache green-jet-compressed `
        -FmmLocalCacheBudgetMb 256 -FmmLocalFallback packed-avx2 `
        -FmmLocalPlanBuilder metric-cache
}

function Invoke-Z4nParserCase([string]$Name,[string[]]$Lines,
        [bool]$ExpectSuccess,[string]$ExpectedMessage='') {
    $dir = Join-Path $WorkRoot $Name
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    $arg = Join-Path $dir 'arg.in'
    $Lines | Set-Content -LiteralPath $arg -Encoding ASCII
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

$acceptedArg = New-Z4nArgs 'accepted_seed'
$accepted = @(Get-Content -LiteralPath $acceptedArg)
Invoke-Z4nParserCase 'accepted' $accepted $true
$profile = @(Import-Csv -LiteralPath (Join-Path $WorkRoot `
    'accepted\pbe_tabi_profile.csv'))[0]
foreach ($column in @('fmm_local_gmres_active','fmm_local_prepare_seconds',
        'fmm_local_plan_build_count','fmm_local_cache_build_count',
        'fmm_local_gmres_matvec_calls','fmm_local_gmres_apply_seconds_total',
        'fmm_local_gmres_apply_seconds_mean','fmm_local_moment_seconds_total',
        'fmm_local_m2l_seconds_total','fmm_local_fallback_seconds_total',
        'fmm_local_l2p_seconds_total',
        'fmm_local_sampled_direct_true_residual','fmm_local_operator_finite',
        'fmm_local_refine_probe_mode','fmm_local_refine_probe_count',
        'fmm_local_solution_probe_used')) {
    if (-not ($profile.PSObject.Properties.Name -contains $column)) {
        throw "accepted profile is missing Z4N column $column"
    }
}
if ($profile.fmm_local_gmres_active -ne 'true' -or
    [int64]$profile.fmm_local_plan_build_count -ne 1 -or
    [int64]$profile.fmm_local_cache_build_count -ne 1 -or
    [int64]$profile.fmm_local_gmres_matvec_calls -le 0 -or
    $profile.fmm_local_operator_finite -ne 'true') {
    throw 'accepted profile did not report a reused finite p4 GMRES operator'
}

$wrongCache = $accepted -replace
    '^-pbe-tabi-fmm-local-cache=green-jet-compressed$',
    '-pbe-tabi-fmm-local-cache=kernel-jet'
Invoke-Z4nParserCase 'wrong_cache' $wrongCache $false `
    'gmres-p4.*green-jet-compressed'

$wrongAudit = $accepted -replace
    '^-pbe-tabi-matvec-audit=sampled-detail$',
    '-pbe-tabi-matvec-audit=sampled'
Invoke-Z4nParserCase 'wrong_audit' $wrongAudit $false `
    'gmres-p4.*sampled-detail'

$wrongPreconditioner = $accepted -replace
    '^-pbe-tabi-precond=fmm-nearfield-ilu0$',
    '-pbe-tabi-precond=diagonal'
Invoke-Z4nParserCase 'wrong_preconditioner' $wrongPreconditioner $false `
    'gmres-p4.*fmm-nearfield-ilu0'

$mixedProbe = $accepted + @(
    '-pbe-tabi-fmm-local-refine-probe-count=8',
    '-pbe-tabi-fmm-local-refine-probe-mode=mixed'
)
Invoke-Z4nParserCase 'mixed_probe' $mixedProbe $true
$mixedProfile = @(Import-Csv -LiteralPath (Join-Path $WorkRoot `
    'mixed_probe\pbe_tabi_profile.csv'))[0]
if ($mixedProfile.fmm_local_refine_probe_mode -ne 'mixed' -or
    [int]$mixedProfile.fmm_local_refine_probe_count -ne 8 -or
    $mixedProfile.fmm_local_solution_probe_used -ne 'true') {
    throw 'mixed probe row did not report requested Z4N.1 probe telemetry'
}

$badProbeCount = $accepted + '-pbe-tabi-fmm-local-refine-probe-count=12'
Invoke-Z4nParserCase 'bad_probe_count' $badProbeCount $false `
    'fmm-local-refine-probe-count.*4.*8.*16'

$badProbeMode = $accepted + '-pbe-tabi-fmm-local-refine-probe-mode=adaptive'
Invoke-Z4nParserCase 'bad_probe_mode' $badProbeMode $false `
    'fmm-local-refine-probe-mode.*rademacher.*solution.*mixed'

$bareOptions = @($accepted | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
if ($bareOptions.Count) { throw "generated arg.in contains bare options" }
Write-Host 'PASS pbe_fagbi_scope_z4n_parser_smoke'
