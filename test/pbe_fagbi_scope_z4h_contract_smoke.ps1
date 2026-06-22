param(
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$scriptPath = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4h_full_protein_theta_recovery.ps1'
if (-not (Test-Path -LiteralPath $scriptPath)) {
    throw "Missing Scope Z4H harness: $scriptPath"
}
$text = Get-Content -Raw -LiteralPath $scriptPath
if ($text -match 'Write-FmmCliSampledPqr') {
    throw 'Z4H must run complete prepared PQR structures without atom sampling.'
}
foreach ($theta in @('0.6','0.5','0.45','0.4','0.35')) {
    if ($text -notmatch [regex]::Escape($theta)) {
        throw "Z4H harness lacks required theta sweep value $theta"
    }
}
foreach ($token in @('AuditSampleCount 512','AuditSampleCount 2048','analytic_singular','fmm-nearfield-ilu0','TaylorOrder 4')) {
    if ($text -notmatch [regex]::Escape($token)) {
        throw "Z4H harness lacks required numerical policy: $token"
    }
}
foreach ($token in @('Full-protein confirmation:','Repeatability: $RepeatabilityCsv','z4h_fixture_specific_policy')) {
    if ($text -notmatch [regex]::Escape($token)) {
        throw "Z4H harness lacks required reporting policy: $token"
    }
}

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing public fixture manifest: $manifestPath" }
$manifest = @(Import-Csv -LiteralPath $manifestPath)
foreach ($fixture in @('barnase','1a63','1sh1')) {
    $row = @($manifest | Where-Object { $_.fixture -eq $fixture })
    if ($row.Count -ne 1) { throw "Z4H fixture manifest lacks exactly one $fixture row." }
}

Write-Host 'PASS pbe_fagbi_scope_z4h_contract_smoke'
