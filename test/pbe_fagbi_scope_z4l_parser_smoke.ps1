param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4l_parser')
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

$baseArgs = @(Get-Content -LiteralPath (Join-Path $baseline 'candidate\arg.in')) +
    '-pbe-tabi-fmm-local=audit-p4'
Invoke-Case -Name 'accepted' -Lines ($baseArgs + @(
        '-pbe-tabi-fmm-local-refine=p3-p4-tail',
        '-pbe-tabi-fmm-local-refine-tol=0.004')) -ExpectSuccess $true
$acceptedSummary = @(Import-Csv (Join-Path $WorkRoot `
    'accepted\pbe_tabi_operator_attribution_summary.csv'))
$acceptedAggregate = @($acceptedSummary | Where-Object {
        $_.approximation_mode -eq 'local_m2l_p4' -and
        $_.probe -eq 'solution' -and $_.row_component -eq 'both' -and
        $_.region -eq 'total' -and $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
if ($null -eq $acceptedAggregate -or
    $acceptedAggregate.fmm_local_refine_mode -ne 'p3-p4-tail' -or
    [Math]::Abs([double]$acceptedAggregate.fmm_local_refine_tolerance - 0.004) -gt 1e-12) {
    throw 'Z4L accepted row did not report refinement mode/tolerance'
}
foreach ($column in @('fmm_local_tail_ratio_before','fmm_local_tail_ratio_after',
        'fmm_local_refined_cells','fmm_local_split_children',
        'fmm_local_direct_fallback_panels','fmm_local_direct_fallback_pairs',
        'fmm_local_far_pairs_before',
        'fmm_local_far_pairs_after','fmm_local_refinement_seconds')) {
    if (-not ($acceptedAggregate.PSObject.Properties.Name -contains $column)) {
        throw "Z4L accepted row is missing telemetry column $column"
    }
}

Invoke-Case -Name 'bad_mode' -Lines ($baseArgs +
    '-pbe-tabi-fmm-local-refine=unknown') -ExpectSuccess $false `
    -ExpectedMessage 'fmm-local-refine.*off.*p3-p4-tail'

Invoke-Case -Name 'bad_tolerance' -Lines ($baseArgs + @(
        '-pbe-tabi-fmm-local-refine=p3-p4-tail',
        '-pbe-tabi-fmm-local-refine-tol=0.1')) -ExpectSuccess $false `
    -ExpectedMessage 'fmm-local-refine-tol.*0.0001.*0.05'

$unsupported = $baseArgs | ForEach-Object {
    $_ -replace '^-pbe-tabi-fmm-local=audit-p4$',
        '-pbe-tabi-fmm-local=off'
}
Invoke-Case -Name 'unsupported' -Lines ($unsupported +
    '-pbe-tabi-fmm-local-refine=p3-p4-tail') -ExpectSuccess $false `
    -ExpectedMessage 'p3-p4-tail.*requires.*audit-p4'

Write-Host 'PASS pbe_fagbi_scope_z4l_parser_smoke'
