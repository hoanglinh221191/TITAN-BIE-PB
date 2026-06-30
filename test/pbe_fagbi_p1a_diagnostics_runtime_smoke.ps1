param(
  [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
  [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\p1a_diagnostics_runtime')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Join-Path $PSScriptRoot 'tmp')

$pqr = Join-Path $WorkRoot 'two_atom.pqr'
Set-Content -LiteralPath $pqr -Encoding ASCII -Value @(
  'ATOM      1  N   GLY A   1      -1.000   0.000   0.000    0.500  1.500',
  'ATOM      2  C   GLY A   1       1.000   0.000   0.000   -0.500  1.500',
  'END'
)

$argPath = Write-FmmCliArgFile `
  -RunDir $WorkRoot `
  -RepoRoot $RepoRoot `
  -PqrPath $pqr `
  -Backend tabi `
  -Mode direct `
  -FmmMode off `
  -SurfacePoints 8 `
  -MaxIter 120 `
  -Restart 16 `
  -Tolerance 1.0e-4 `
  -Threads 1 `
  -SaltM 0.15 `
  -InternalDielectric 1.0 `
  -ExternalDielectric 80.0 `
  -TemperatureK 298.15 `
  -OpeningTheta 0.6 `
  -AuditSampleCount 0 `
  -AuditMode off `
  -Preconditioner diagonal `
  -NearIntegration analytic_singular `
  -TaylorOrder 4

Add-Content -LiteralPath $argPath -Encoding ASCII -Value '-pbe-tabi-diagnostics=audit'
$badOptions = @(Get-Content -LiteralPath $argPath | Where-Object {
    -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
  })
if ($badOptions.Count -ne 0) {
  throw "Generated arg.in contains non equals-syntax options: $($badOptions -join '; ')"
}

$run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $WorkRoot `
  -ArgPath $argPath -CaseTimeoutSeconds 120
if ($run.exit_code -ne 0) {
  $stderr = Get-Content -Raw -LiteralPath (Join-Path $WorkRoot 'stderr.log')
  $stdout = Get-Content -Raw -LiteralPath (Join-Path $WorkRoot 'stdout.log')
  throw "TITAN diagnostics runtime smoke failed with exit $($run.exit_code).`nSTDOUT:`n$stdout`nSTDERR:`n$stderr"
}

$auditPath = Join-Path $WorkRoot 'pbe_tabi_audit_profile.csv'
if (-not (Test-Path -LiteralPath $auditPath)) {
  throw "Missing pbe_tabi_audit_profile.csv"
}
$rows = @(Import-Csv -LiteralPath $auditPath)
if ($rows.Count -ne 1) {
  throw "Expected exactly one compact diagnostics row, got $($rows.Count)."
}
$required = @(
  'fixture','panels','energy_kJmol','iterations','internal_residual',
  'sampled_direct_residual','audit_rows','wall_seconds','rss_mb',
  'p_order','theta','guard','tail_tol','ilu_cap','far_pairs',
  'fallback_pairs','status'
)
foreach ($column in $required) {
  if (-not ($rows[0].PSObject.Properties.Name -contains $column)) {
    throw "Compact diagnostics row missing column: $column"
  }
}
if ([string]::IsNullOrWhiteSpace([string]$rows[0].fixture)) {
  throw 'Compact diagnostics fixture label is empty.'
}

Write-Host 'PASS pbe_fagbi_p1a_diagnostics_runtime_smoke'
