param(
    [string]$TitanExe,
    [switch]$Force
)

$helper = Join-Path $PSScriptRoot '..\common\Invoke-TitanPbeExample.ps1'
& $helper -ExampleDirectory $PSScriptRoot -TitanExe $TitanExe `
    -PqrFiles @('6lzg_complex_charmm_ph7.pqr', '6lzg_ace2_chain_A_charmm_ph7.pqr', '6lzg_rbd_chain_B_charmm_ph7.pqr') -Force:$Force
