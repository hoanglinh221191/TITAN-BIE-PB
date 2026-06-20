param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4m_parser')
)

$ErrorActionPreference = 'Stop'
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$baseline = Join-Path $WorkRoot 'baseline'
& (Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4k_parser_smoke.ps1') `
    -TitanExe $TitanExe -WorkRoot $baseline

function Invoke-Case {
    param([string]$Name, [string[]]$Lines, [bool]$ExpectSuccess,
          [string]$ExpectedMessage = '')
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

$baseArgs = @(Get-Content -LiteralPath (Join-Path $baseline 'candidate\arg.in'))
$baseArgs = @($baseArgs | Where-Object { $_ -notmatch '^-pbe-tabi-fmm-local=' })
$accepted = $baseArgs + @(
    '-pbe-tabi-fmm-local=single-matvec-p4',
    '-pbe-tabi-fmm-local-refine=operator-tail-static',
    '-pbe-tabi-fmm-local-refine-tol=0.002',
    '-pbe-tabi-fmm-local-cache=kernel-jet',
    '-pbe-tabi-fmm-local-cache-budget-mb=256')
Invoke-Case -Name 'accepted' -Lines $accepted -ExpectSuccess $true

$summary = @(Import-Csv (Join-Path $WorkRoot `
    'accepted\pbe_tabi_operator_attribution_summary.csv'))
$row = @($summary | Where-Object {
        $_.approximation_mode -eq 'local_m2l_p4_static_cached' -and
        $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
        $_.region -eq 'total' -and $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
if ($null -eq $row) { throw 'Z4M accepted row is missing' }
foreach ($column in @('fmm_local_plan_hash','fmm_local_plan_build_count',
        'fmm_local_static_probe_count','fmm_local_cache_mode',
        'fmm_local_cache_operators','fmm_local_cache_bytes',
        'fmm_local_cache_hits','fmm_local_candidate_apply_seconds')) {
    if (-not ($row.PSObject.Properties.Name -contains $column)) {
        throw "Z4M accepted row is missing telemetry column $column"
    }
}

Invoke-Case -Name 'missing_static' -Lines ($baseArgs + @(
        '-pbe-tabi-fmm-local=single-matvec-p4',
        '-pbe-tabi-fmm-local-cache=kernel-jet')) -ExpectSuccess $false `
    -ExpectedMessage 'kernel-jet.*operator-tail-static'
Invoke-Case -Name 'missing_cache' -Lines ($baseArgs + @(
        '-pbe-tabi-fmm-local=single-matvec-p4',
        '-pbe-tabi-fmm-local-refine=operator-tail-static')) -ExpectSuccess $false `
    -ExpectedMessage 'single-matvec-p4.*kernel-jet'
Invoke-Case -Name 'bad_cache' -Lines ($accepted -replace `
        '^-pbe-tabi-fmm-local-cache=kernel-jet$',
        '-pbe-tabi-fmm-local-cache=unknown') -ExpectSuccess $false `
    -ExpectedMessage 'fmm-local-cache.*off.*kernel-jet'
Invoke-Case -Name 'bad_budget' -Lines ($accepted -replace `
        '^-pbe-tabi-fmm-local-cache-budget-mb=256$',
        '-pbe-tabi-fmm-local-cache-budget-mb=32') -ExpectSuccess $false `
    -ExpectedMessage 'cache-budget-mb.*64.*1024'
Invoke-Case -Name 'gmres_closed' -Lines ($accepted -replace `
        '^-pbe-tabi-fmm-local=single-matvec-p4$',
        '-pbe-tabi-fmm-local=gmres-p4') -ExpectSuccess $false `
    -ExpectedMessage 'gmres-p4.*fail-closed'

Write-Host 'PASS pbe_fagbi_scope_z4m_parser_smoke'
