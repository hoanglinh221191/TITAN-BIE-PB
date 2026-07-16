param(
    [string]$TitanExe,
    [switch]$Force
)

$helper = Join-Path $PSScriptRoot '..\common\Invoke-TitanPbeExample.ps1'
& $helper -ExampleDirectory $PSScriptRoot -TitanExe $TitanExe `
    -PqrFiles @('1vii_full_amber_ph7.pqr') -Force:$Force

