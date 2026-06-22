param(
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public')
)

$ErrorActionPreference = 'Stop'
$prepScript = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4g_prepare_public_fixtures.ps1'
$text = Get-Content -Raw -LiteralPath $prepScript
if ($text -match 'ph-calc-method') {
    throw 'PDB2PQR 3.7.1 preparation must not use the obsolete --ph-calc-method option.'
}
if ($text -notmatch 'titration-state-method=propka') {
    throw 'PDB2PQR 3.7.1 preparation must pin --titration-state-method=propka.'
}

$manifest = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
if (Test-Path -LiteralPath $manifest) {
    $rows = @(Import-Csv -LiteralPath $manifest)
    if ($rows.Count -ne 3) { throw 'Fixture manifest must contain three full-protein rows.' }
    foreach ($row in $rows) {
        foreach ($field in @('raw_pdb_sha256','prepared_pdb_sha256','pqr_sha256')) {
            if ($row.PSObject.Properties.Name -notcontains $field -or [string]$row.$field -notmatch '^[0-9a-f]{64}$') {
                throw "Fixture manifest lacks a valid $field provenance hash: $($row.fixture)"
            }
        }
        if ($row.preparation -match 'ph-calc-method' -or
            $row.preparation -notmatch 'titration-state-method=propka') {
            throw "Fixture preparation metadata is incompatible with PDB2PQR 3.7.1: $($row.fixture)"
        }
        if ($row.PSObject.Properties.Name -notcontains 'radius_floor_angstrom' -or
            [double]$row.radius_floor_angstrom -ne 0.6) {
            throw "Fixture manifest lacks the TITAN-compatible 0.6 A radius-floor provenance: $($row.fixture)"
        }
        $pqrPath = Join-Path $FixtureRoot $row.pqr
        foreach ($line in Get-Content -LiteralPath $pqrPath) {
            if (-not ($line.StartsWith('ATOM') -or $line.StartsWith('HETATM'))) { continue }
            $fields = @($line -split '\s+' | Where-Object { $_ -ne '' })
            if ($fields.Count -lt 10 -or [double]$fields[-1] -le 0.0) {
                throw "Fixture contains a non-positive or malformed PQR radius: $($row.fixture) :: $line"
            }
        }
    }
}

Write-Host 'PASS pbe_fagbi_scope_z4g_fixture_contract_smoke'
