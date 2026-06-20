param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4m1_parser')
)

$ErrorActionPreference = 'Stop'
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$baseline = Join-Path $WorkRoot 'baseline'
& (Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4m_parser_smoke.ps1') `
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
    Where-Object {
        $_ -notmatch '^-pbe-tabi-near-integration=' -and
        $_ -notmatch '^-pbe-tabi-quad-buffer=' -and
        $_ -notmatch '^-pbe-tabi-fmm-local-fallback='
    }) + @('-pbe-tabi-near-integration=analytic_singular',
           '-pbe-tabi-quad-buffer=off')

Invoke-Case 'accepted' ($baseArgs + '-pbe-tabi-fmm-local-fallback=packed-avx2') $true
$summary = @(Import-Csv (Join-Path $WorkRoot `
    'accepted\pbe_tabi_operator_attribution_summary.csv'))
$row = @($summary | Where-Object {
    $_.approximation_mode -eq 'local_m2l_p4_static_cached' -and
    $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
    $_.region -eq 'total' -and $_.channel -eq 'aggregate'
}) | Select-Object -First 1
if ($null -eq $row -or $row.fmm_local_fallback_effective -ne 'packed-avx2') {
    throw 'accepted packed row did not report effective packed-avx2 mode'
}
foreach ($column in @('fmm_local_fallback_worklist_hash',
        'fmm_local_fallback_worklist_build_count','fmm_local_fallback_tasks',
        'fmm_local_fallback_metadata_bytes','fmm_local_fallback_base_near_pairs',
        'fmm_local_fallback_refined_pairs','fmm_local_fallback_vectorized_lanes',
        'fmm_local_fallback_scalar_self_pairs','fmm_local_fallback_seconds')) {
    if (-not ($row.PSObject.Properties.Name -contains $column)) {
        throw "accepted packed row is missing telemetry column $column"
    }
}

Invoke-Case 'scalar' ($baseArgs + '-pbe-tabi-fmm-local-fallback=scalar') $true
Invoke-Case 'bad_mode' ($baseArgs + '-pbe-tabi-fmm-local-fallback=unknown') $false `
    'fmm-local-fallback.*scalar.*packed-avx2'
Invoke-Case 'wrong_local' (($baseArgs -replace `
    '^-pbe-tabi-fmm-local=single-matvec-p4$',
    '-pbe-tabi-fmm-local=audit-p4') +
    '-pbe-tabi-fmm-local-fallback=packed-avx2') $false `
    'packed-avx2.*single-matvec-p4'
Invoke-Case 'wrong_near' (($baseArgs -replace `
    '^-pbe-tabi-near-integration=analytic_singular$',
    '-pbe-tabi-near-integration=centroid') +
    '-pbe-tabi-fmm-local-fallback=packed-avx2') $false `
    'packed-avx2.*analytic_singular'
Invoke-Case 'wrong_quad' (($baseArgs -replace `
    '^-pbe-tabi-quad-buffer=off$',
    '-pbe-tabi-quad-buffer=adaptive') +
    '-pbe-tabi-fmm-local-fallback=packed-avx2') $false `
    'packed-avx2.*quad-buffer=off'

Write-Host 'PASS pbe_fagbi_scope_z4m1_parser_smoke'
