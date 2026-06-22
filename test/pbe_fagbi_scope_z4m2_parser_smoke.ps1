param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4m2_parser')
)

$ErrorActionPreference = 'Stop'
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$baseline = Join-Path $WorkRoot 'baseline'
& (Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4m1_parser_smoke.ps1') `
    -TitanExe $TitanExe -WorkRoot $baseline

function Invoke-Case {
    param([string]$Name,[string[]]$Lines,[bool]$ExpectSuccess,
          [string]$ExpectedMessage='')
    $dir = Join-Path $WorkRoot $Name
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    $argPath = Join-Path $dir 'arg.in'
    $Lines | Set-Content -LiteralPath $argPath -Encoding ASCII
    $stdout = Join-Path $dir 'stdout.log'
    $stderr = Join-Path $dir 'stderr.log'
    $process = Start-Process -FilePath $TitanExe -ArgumentList @('-a',$argPath) `
        -WorkingDirectory $dir -RedirectStandardOutput $stdout `
        -RedirectStandardError $stderr -Wait -PassThru -NoNewWindow
    $log = (Get-Content -Raw $stdout) + (Get-Content -Raw $stderr)
    if ($ExpectSuccess -and $process.ExitCode -ne 0) {
        throw "$Name unexpectedly failed: $log"
    }
    if (-not $ExpectSuccess -and $process.ExitCode -eq 0) {
        throw "$Name unexpectedly succeeded"
    }
    if (-not $ExpectSuccess -and $log -notmatch $ExpectedMessage) {
        throw "$Name error was unclear: $log"
    }
}

$baseArgs = @(Get-Content (Join-Path $baseline 'accepted\arg.in') |
    Where-Object { $_ -notmatch '^-pbe-tabi-fmm-local-plan-builder=' })

Invoke-Case 'accepted' ($baseArgs +
    '-pbe-tabi-fmm-local-plan-builder=metric-cache') $true
$summary = @(Import-Csv (Join-Path $WorkRoot `
    'accepted\pbe_tabi_operator_attribution_summary.csv'))
$row = @($summary | Where-Object {
    $_.approximation_mode -eq 'local_m2l_p4_static_cached' -and
    $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
    $_.region -eq 'total' -and $_.channel -eq 'aggregate'
}) | Select-Object -First 1
if ($null -eq $row -or $row.fmm_local_plan_builder -ne 'metric-cache') {
    throw 'accepted row did not report metric-cache plan builder'
}
foreach ($column in @('fmm_local_plan_probe_moment_seconds',
        'fmm_local_plan_pair_metric_seconds',
        'fmm_local_plan_refinement_seconds',
        'fmm_local_plan_pair_metric_evaluations',
        'fmm_local_plan_pair_metric_cache_hits',
        'fmm_local_plan_refinement_iterations',
        'fmm_local_plan_metric_cache_peak_bytes',
        'fmm_local_plan_parallel_target_leaves')) {
    if (-not ($row.PSObject.Properties.Name -contains $column)) {
        throw "accepted row is missing telemetry column $column"
    }
}

Invoke-Case 'legacy_default' $baseArgs $true
Invoke-Case 'bad_mode' ($baseArgs +
    '-pbe-tabi-fmm-local-plan-builder=unknown') $false `
    'fmm-local-plan-builder.*legacy-rescan.*metric-cache'
$wrongRefine = (($baseArgs -replace `
    '^-pbe-tabi-fmm-local=single-matvec-p4$',
    '-pbe-tabi-fmm-local=audit-p4') -replace `
    '^-pbe-tabi-fmm-local-fallback=packed-avx2$',
    '-pbe-tabi-fmm-local-fallback=scalar') -replace `
    '^-pbe-tabi-fmm-local-refine=operator-tail-static$',
    '-pbe-tabi-fmm-local-refine=p3-p4-tail' -replace `
    '^-pbe-tabi-fmm-local-cache=kernel-jet$',
    '-pbe-tabi-fmm-local-cache=off'
Invoke-Case 'wrong_refine' ($wrongRefine +
    '-pbe-tabi-fmm-local-plan-builder=metric-cache') $false `
    'metric-cache.*operator-tail-static'
$wrongLocal = ($baseArgs -replace `
    '^-pbe-tabi-fmm-local=single-matvec-p4$',
    '-pbe-tabi-fmm-local=off') -replace `
    '^-pbe-tabi-fmm-local-fallback=packed-avx2$',
    '-pbe-tabi-fmm-local-fallback=scalar'
Invoke-Case 'wrong_local' ($wrongLocal +
    '-pbe-tabi-fmm-local-plan-builder=metric-cache') $false `
    'operator-tail-static.*audit-p4.*single-matvec-p4'

Write-Host 'PASS pbe_fagbi_scope_z4m2_parser_smoke'
