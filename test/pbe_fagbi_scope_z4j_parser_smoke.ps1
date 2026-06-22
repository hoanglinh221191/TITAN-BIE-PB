param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4j_parser')
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$baseline = Join-Path $WorkRoot 'baseline'
& (Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4i_parser_smoke.ps1') `
    -TitanExe $TitanExe -WorkRoot $baseline

$candidate = Join-Path $WorkRoot 'candidate'
New-Item -ItemType Directory -Force -Path $candidate | Out-Null
$candidateArg = Join-Path $candidate 'arg.in'
@(Get-Content -LiteralPath (Join-Path $baseline 'arg.in')) +
    '-pbe-tabi-surface-points=24' +
    '-pbe-tabi-taylor-order=4' +
    '-pbe-tabi-fmm-cartesian-audit=on' |
    Set-Content -LiteralPath $candidateArg -Encoding ASCII
$stdout = Join-Path $candidate 'stdout.log'
$stderr = Join-Path $candidate 'stderr.log'
$process = Start-Process -FilePath $TitanExe -ArgumentList @('-a',$candidateArg) `
    -WorkingDirectory $candidate -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr -Wait -PassThru -NoNewWindow
if ($process.ExitCode -ne 0) {
    $combined = (Get-Content -Raw $stdout) + (Get-Content -Raw $stderr)
    throw "Cartesian attribution candidate failed (exit=$($process.ExitCode)): $combined"
}
$summaryPath = Join-Path $candidate 'pbe_tabi_operator_attribution_summary.csv'
if (-not (Test-Path -LiteralPath $summaryPath)) {
    throw 'Cartesian attribution candidate did not write the summary CSV'
}
$summary = @(Import-Csv -LiteralPath $summaryPath)
$modes = @($summary.approximation_mode | Sort-Object -Unique)
foreach ($mode in @('legacy_degree1','cartesian_p4')) {
    if ($mode -notin $modes) {
        throw "Missing approximation mode in sampled-detail summary: $mode"
    }
}
$legacyAggregate = @($summary | Where-Object {
        $_.approximation_mode -eq 'legacy_degree1' -and $_.probe -eq 'solution' -and
        $_.row_component -eq 'both' -and $_.region -eq 'total' -and
        $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
$cartesianAggregate = @($summary | Where-Object {
        $_.approximation_mode -eq 'cartesian_p4' -and $_.probe -eq 'solution' -and
        $_.row_component -eq 'both' -and $_.region -eq 'total' -and
        $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
if ($null -eq $legacyAggregate -or $null -eq $cartesianAggregate) {
    throw 'Cartesian attribution did not report both aggregate operator rows'
}
if ([double]$legacyAggregate.aggregate_audit_rel -gt 1e-14 -and
    [double]$cartesianAggregate.aggregate_audit_rel -ge
        [double]$legacyAggregate.aggregate_audit_rel) {
    throw 'Cartesian attribution did not report its independently improved operator RMS'
}

$unsupported = Join-Path $WorkRoot 'unsupported_sampled'
New-Item -ItemType Directory -Force -Path $unsupported | Out-Null
$unsupportedArg = Join-Path $unsupported 'arg.in'
Get-Content -LiteralPath $candidateArg | ForEach-Object {
    $_ -replace '^-pbe-tabi-matvec-audit=sampled-detail$',
        '-pbe-tabi-matvec-audit=sampled'
} | Set-Content -LiteralPath $unsupportedArg -Encoding ASCII
$unsupportedOut = Join-Path $unsupported 'stdout.log'
$unsupportedErr = Join-Path $unsupported 'stderr.log'
$unsupportedProcess = Start-Process -FilePath $TitanExe `
    -ArgumentList @('-a',$unsupportedArg) -WorkingDirectory $unsupported `
    -RedirectStandardOutput $unsupportedOut -RedirectStandardError $unsupportedErr `
    -Wait -PassThru -NoNewWindow
if ($unsupportedProcess.ExitCode -eq 0) {
    throw 'Cartesian attribution unexpectedly accepted matvec-audit=sampled'
}
$unsupportedLog = (Get-Content -Raw $unsupportedOut) +
    (Get-Content -Raw $unsupportedErr)
if ($unsupportedLog -notmatch 'requires.*sampled-detail') {
    throw "Cartesian attribution unsupported-path error was unclear: $unsupportedLog"
}

Write-Host 'PASS pbe_fagbi_scope_z4j_parser_smoke'
