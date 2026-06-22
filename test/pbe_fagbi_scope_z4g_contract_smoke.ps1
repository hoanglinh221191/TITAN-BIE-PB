param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4g_contract')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$prepScript = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4g_prepare_public_fixtures.ps1'
$benchmarkScript = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4g_full_protein_benchmark.ps1'
$builderScript = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4g_build_docx.py'
$orchestratorScript = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4g_external_manuscript.ps1'

foreach ($required in @($prepScript, $benchmarkScript, $builderScript, $orchestratorScript)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing Scope Z4G component: $required"
    }
}

$prepText = Get-Content -Raw -LiteralPath $prepScript
$benchmarkText = Get-Content -Raw -LiteralPath $benchmarkScript
if ($prepText -notmatch '1A63' -or $prepText -notmatch '1SH1' -or
    $prepText -notmatch '3\.7\.1' -or $prepText -notmatch 'AMBER') {
    throw 'Fixture preparation contract lacks public IDs or pinned PDB2PQR policy.'
}
if ($benchmarkText -match 'Write-FmmCliSampledPqr') {
    throw 'Z4G full-protein benchmark must not call Write-FmmCliSampledPqr.'
}

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$matrix = Join-Path $WorkRoot 'matrix.csv'
$provenance = Join-Path $WorkRoot 'provenance.csv'
$trace = Join-Path $WorkRoot 'trace.csv'
$runtimeFigure = Join-Path $WorkRoot 'runtime.svg'
$convergenceFigure = Join-Path $WorkRoot 'convergence.svg'
$docx = Join-Path $WorkRoot 'titan_fagbi_external_manuscript_z4g.docx'
$source = Join-Path $WorkRoot 'manuscript.md'
$claimTrace = Join-Path $WorkRoot 'claim_trace.csv'
$qa = Join-Path $WorkRoot 'qa.csv'

@(
    [pscustomobject]@{ protein_name='Barnase'; pdb_id='PyDelPhi benchmark'; atom_count=1100; panels=42000; energy_kj_mol=-1000; true_residual=1e-8; audit_rms=0.002; iterations=90; wall_seconds=40; peak_rss_mb=80; status='ok' },
    [pscustomobject]@{ protein_name='RNA-binding domain of E. coli Rho factor'; pdb_id='1A63'; atom_count=2000; panels=90000; energy_kj_mol=-2000; true_residual=1e-8; audit_rms=0.003; iterations=100; wall_seconds=80; peak_rss_mb=120; status='ok' },
    [pscustomobject]@{ protein_name='Stichodactyla helianthus neurotoxin I'; pdb_id='1SH1'; atom_count=700; panels=28000; energy_kj_mol=-700; true_residual=1e-9; audit_rms=0.0025; iterations=70; wall_seconds=25; peak_rss_mb=60; status='ok' }
) | Export-Csv -LiteralPath $matrix -NoTypeInformation

@(
    [pscustomobject]@{ pdb_id='1A63'; protein_name='RNA-binding domain of E. coli Rho factor'; model='1'; chain='A'; pdb_sha256='a'; pqr_sha256='b'; pdb2pqr_version='3.7.1'; preparation='--ff=AMBER --ffout=AMBER --keep-chain --titration-state-method=propka --with-ph=7.0' },
    [pscustomobject]@{ pdb_id='1SH1'; protein_name='Stichodactyla helianthus neurotoxin I'; model='1'; chain='A'; pdb_sha256='c'; pqr_sha256='d'; pdb2pqr_version='3.7.1'; preparation='--ff=AMBER --ffout=AMBER --keep-chain --titration-state-method=propka --with-ph=7.0' }
) | Export-Csv -LiteralPath $provenance -NoTypeInformation

@(
    [pscustomobject]@{ preconditioner='block_jacobi'; iteration=0; residual=1.0 },
    [pscustomobject]@{ preconditioner='block_jacobi'; iteration=100; residual=1e-3 },
    [pscustomobject]@{ preconditioner='sparse_ilu0'; iteration=0; residual=1.0 },
    [pscustomobject]@{ preconditioner='sparse_ilu0'; iteration=80; residual=1e-9 }
) | Export-Csv -LiteralPath $trace -NoTypeInformation

Set-Content -LiteralPath $runtimeFigure -Encoding ASCII -Value '<svg xmlns="http://www.w3.org/2000/svg" width="800" height="400"><rect width="800" height="400" fill="white"/></svg>'
Set-Content -LiteralPath $convergenceFigure -Encoding ASCII -Value '<svg xmlns="http://www.w3.org/2000/svg" width="800" height="400"><rect width="800" height="400" fill="white"/></svg>'

$python = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $python $builderScript --matrix $matrix --provenance $provenance --trace $trace `
    --runtime-figure $runtimeFigure --convergence-figure $convergenceFigure `
    --output-docx $docx --output-source $source --claim-trace $claimTrace --qa $qa
if ($LASTEXITCODE -ne 0) { throw 'Z4G DOCX builder failed.' }

foreach ($artifact in @($docx, $source, $claimTrace, $qa)) {
    if (-not (Test-Path -LiteralPath $artifact)) { throw "Missing Z4G artifact: $artifact" }
}

$sourceText = Get-Content -Raw -LiteralPath $source
foreach ($forbidden in @('cap-128','cap-256','research-medium','research-alt','money plot','Z1-Z5','Z4 reference policy','audit-clean')) {
    if ($sourceText -match [regex]::Escape($forbidden)) {
        throw "Internal manuscript jargon leaked into Z4G source: $forbidden"
    }
}
foreach ($requiredText in @('1A63','1SH1','N_atoms','N_panels','block-CSR','2 x 2','no fill','forward substitution','backward substitution')) {
    if ($sourceText -notmatch [regex]::Escape($requiredText)) {
        throw "Z4G manuscript source lacks required public/mathematical text: $requiredText"
    }
}

$qaRows = @(Import-Csv -LiteralPath $qa)
if (@($qaRows | Where-Object { $_.status -ne 'pass' }).Count -ne 0) {
    throw 'Z4G DOCX QA contains failing rows.'
}

Write-Host 'PASS pbe_fagbi_scope_z4g_contract_smoke'
