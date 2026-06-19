param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4k_parser')
)

$ErrorActionPreference = 'Stop'
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$baseline = Join-Path $WorkRoot 'baseline'
& (Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4j_parser_smoke.ps1') `
    -TitanExe $TitanExe -WorkRoot $baseline

$candidate = Join-Path $WorkRoot 'candidate'
New-Item -ItemType Directory -Force -Path $candidate | Out-Null
$argPath = Join-Path $candidate 'arg.in'
@(Get-Content -LiteralPath (Join-Path $baseline 'candidate\arg.in')) +
    '-pbe-tabi-fmm-local=audit-p4' |
    Set-Content -LiteralPath $argPath -Encoding ASCII
$stdout = Join-Path $candidate 'stdout.log'
$stderr = Join-Path $candidate 'stderr.log'
$process = Start-Process -FilePath $TitanExe -ArgumentList @('-a',$argPath) `
    -WorkingDirectory $candidate -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr -Wait -PassThru -NoNewWindow
if ($process.ExitCode -ne 0) {
    throw "Z4K audit-p4 candidate failed: $((Get-Content -Raw $stdout) + (Get-Content -Raw $stderr))"
}
$summary = @(Import-Csv (Join-Path $candidate 'pbe_tabi_operator_attribution_summary.csv'))
if ('local_m2l_p4' -notin @($summary.approximation_mode | Sort-Object -Unique)) {
    throw 'Z4K sampled-detail output is missing local_m2l_p4 rows'
}
$localAggregate = @($summary | Where-Object {
        $_.approximation_mode -eq 'local_m2l_p4' -and $_.probe -eq 'solution' -and
        $_.row_component -eq 'both' -and $_.region -eq 'total' -and
        $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
foreach ($column in @('fmm_local_mode','fmm_local_source_order',
        'fmm_local_target_order','fmm_local_kernel_derivative_order',
        'fmm_local_target_leaves','fmm_local_m2l_pairs',
        'fmm_local_moment_seconds','fmm_local_m2l_seconds',
        'fmm_local_l2p_seconds','fmm_local_finite_outputs')) {
    if ($null -eq $localAggregate -or
        -not ($localAggregate.PSObject.Properties.Name -contains $column)) {
        throw "Z4K sampled-detail summary is missing telemetry column $column"
    }
}
if ($localAggregate.fmm_local_mode -ne 'audit-p4' -or
    [int]$localAggregate.fmm_local_source_order -ne 4 -or
    [int]$localAggregate.fmm_local_target_order -ne 4 -or
    [int]$localAggregate.fmm_local_kernel_derivative_order -ne 8) {
    throw 'Z4K local-order telemetry is inconsistent'
}

$unsupported = Join-Path $WorkRoot 'unsupported'
New-Item -ItemType Directory -Force -Path $unsupported | Out-Null
$unsupportedArg = Join-Path $unsupported 'arg.in'
Get-Content $argPath | ForEach-Object {
    ($_ -replace '^-pbe-tabi-matvec-audit=sampled-detail$',
        '-pbe-tabi-matvec-audit=sampled') -replace
        '^-pbe-tabi-fmm-cartesian-audit=on$',
        '-pbe-tabi-fmm-cartesian-audit=off'
} | Set-Content $unsupportedArg -Encoding ASCII
$unsupportedOut = Join-Path $unsupported 'stdout.log'
$unsupportedErr = Join-Path $unsupported 'stderr.log'
$unsupportedProcess = Start-Process -FilePath $TitanExe `
    -ArgumentList @('-a',$unsupportedArg) -WorkingDirectory $unsupported `
    -RedirectStandardOutput $unsupportedOut -RedirectStandardError $unsupportedErr `
    -Wait -PassThru -NoNewWindow
if ($unsupportedProcess.ExitCode -eq 0) {
    throw 'audit-p4 unexpectedly accepted matvec-audit=sampled'
}
$unsupportedLog = (Get-Content -Raw $unsupportedOut) +
    (Get-Content -Raw $unsupportedErr)
if ($unsupportedLog -notmatch 'audit-p4.*requires.*sampled-detail') {
    throw "audit-p4 unsupported-path error was unclear: $unsupportedLog"
}

$gmres = Join-Path $WorkRoot 'gmres_fail_closed'
New-Item -ItemType Directory -Force -Path $gmres | Out-Null
$gmresArg = Join-Path $gmres 'arg.in'
Get-Content $argPath | ForEach-Object {
    $_ -replace '^-pbe-tabi-fmm-local=audit-p4$',
        '-pbe-tabi-fmm-local=gmres-p4'
} | Set-Content $gmresArg -Encoding ASCII
$gmresOut = Join-Path $gmres 'stdout.log'
$gmresErr = Join-Path $gmres 'stderr.log'
$gmresProcess = Start-Process -FilePath $TitanExe -ArgumentList @('-a',$gmresArg) `
    -WorkingDirectory $gmres -RedirectStandardOutput $gmresOut `
    -RedirectStandardError $gmresErr -Wait -PassThru -NoNewWindow
if ($gmresProcess.ExitCode -eq 0) {
    throw 'gmres-p4 unexpectedly bypassed the fail-closed gate'
}
if (((Get-Content -Raw $gmresOut) + (Get-Content -Raw $gmresErr)) -notmatch
    'gmres-p4.*fail-closed') {
    throw 'gmres-p4 fail-closed error was unclear'
}

Write-Host 'PASS pbe_fagbi_scope_z4k_parser_smoke'
