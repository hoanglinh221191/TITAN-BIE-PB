param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4m3_parser')
)

$ErrorActionPreference = 'Stop'
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$baseline = Join-Path $WorkRoot 'baseline'
& (Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4m2_parser_smoke.ps1') `
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
    Where-Object { $_ -notmatch '^-pbe-tabi-fmm-local-cache=' })
Invoke-Case 'accepted' ($baseArgs +
    '-pbe-tabi-fmm-local-cache=green-jet-compressed') $true
$summary = @(Import-Csv (Join-Path $WorkRoot `
    'accepted\pbe_tabi_operator_attribution_summary.csv'))
$row = @($summary | Where-Object {
    $_.approximation_mode -eq 'local_m2l_p4_static_green_cached' -and
    $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
    $_.region -eq 'total' -and $_.channel -eq 'aggregate'
}) | Select-Object -First 1
if ($null -eq $row -or
    $row.fmm_local_cache_mode -ne 'green-jet-compressed' -or
    $row.fmm_local_cache_layout -ne 'green-jet-compressed') {
    throw 'accepted row did not report compressed Green-jet cache layout'
}
foreach ($column in @('fmm_local_cache_operator_bytes',
        'fmm_local_cache_uncompressed_estimated_bytes',
        'fmm_local_cache_compression_ratio',
        'fmm_local_cache_far_pair_offsets',
        'fmm_local_cache_build_seconds','fmm_local_cache_hits')) {
    if (-not ($row.PSObject.Properties.Name -contains $column)) {
        throw "accepted row is missing telemetry column $column"
    }
}
Invoke-Case 'bad_mode' ($baseArgs +
    '-pbe-tabi-fmm-local-cache=compressed-unknown') $false `
    'fmm-local-cache.*green-jet-compressed'
$gmres = ($baseArgs -replace '^-pbe-tabi-fmm-local=single-matvec-p4$',
    '-pbe-tabi-fmm-local=gmres-p4')
Invoke-Case 'gmres_rejected' ($gmres +
    '-pbe-tabi-fmm-local-cache=green-jet-compressed') $false `
    'gmres-p4.*fail-closed'

Write-Host 'PASS pbe_fagbi_scope_z4m3_parser_smoke'
