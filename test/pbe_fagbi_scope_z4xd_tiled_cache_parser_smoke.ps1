param()

$ErrorActionPreference = 'Stop'

$Root = Split-Path -Parent $PSScriptRoot
$Tmp = Join-Path $Root 'test/tmp/fagbi_scope_z4xd_parser'
New-Item -ItemType Directory -Path $Tmp -Force | Out-Null

function Write-Z4xdParserArg {
    param(
        [Parameter(Mandatory=$true)][string]$Dir,
        [Parameter(Mandatory=$true)][string[]]$Lines
    )
    New-Item -ItemType Directory -Path $Dir -Force | Out-Null
    $path = Join-Path $Dir 'arg.in'
    $Lines | Set-Content -LiteralPath $path -Encoding ASCII
    $bare = @(Get-Content -LiteralPath $path | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) {
        throw "Generated arg.in contains bare options: $($bare -join ', ')"
    }
    return $path
}

[void](Write-Z4xdParserArg -Dir (Join-Path $Tmp 'valid') -Lines @(
    '-pbe=tabi',
    '-pbe-tabi-sum=treecode',
    '-pbe-tabi-fmm-m2l=single-matvec-analytic',
    '-pbe-tabi-fmm-local=single-matvec-p4',
    '-pbe-tabi-fmm-local-refine=operator-tail-static',
    '-pbe-tabi-fmm-local-plan-builder=metric-cache',
    '-pbe-tabi-fmm-local-cache=green-jet-tiled',
    '-pbe-tabi-fmm-local-fallback=packed-avx2',
    '-pbe-tabi-taylor-order=5',
    '-pbe-tabi-matvec-audit=sampled-detail',
    '-pbe-tabi-near-integration=analytic_singular',
    '-pbe-tabi-quad-buffer=off',
    '-pbe-test-pqr=test/input/placeholder.pqr'
))

[void](Write-Z4xdParserArg -Dir (Join-Path $Tmp 'invalid_default') -Lines @(
    '-pbe=tabi',
    '-pbe-tabi-sum=treecode',
    '-pbe-tabi-fmm-local=off',
    '-pbe-tabi-fmm-local-cache=green-jet-tiled'
))

[void](Write-Z4xdParserArg -Dir (Join-Path $Tmp 'invalid_audit') -Lines @(
    '-pbe=tabi',
    '-pbe-tabi-sum=treecode',
    '-pbe-tabi-fmm-m2l=single-matvec-analytic',
    '-pbe-tabi-fmm-local=single-matvec-p4',
    '-pbe-tabi-fmm-local-cache=green-jet-tiled',
    '-pbe-tabi-matvec-audit=sampled'
))

Write-Output 'parse smoke generated'
