[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string]$WorkRoot = 'test\tmp\pbe_tabi2_nodepatch_analytic_ses_control',
    [ValidateRange(30,1800)][int]$CaseTimeoutSeconds = 300
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
$TitanExe = if ([IO.Path]::IsPathRooted($TitanExe)) {
    [IO.Path]::GetFullPath($TitanExe)
} else { [IO.Path]::GetFullPath((Join-Path $RepoRoot $TitanExe)) }
$WorkRoot = if ([IO.Path]::IsPathRooted($WorkRoot)) {
    [IO.Path]::GetFullPath($WorkRoot)
} else { [IO.Path]::GetFullPath((Join-Path $RepoRoot $WorkRoot)) }
if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing TITAN executable: $TitanExe" }
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)

$pqr = Join-Path $WorkRoot 'analytic_triple.pqr'
$root3 = [Math]::Sqrt(3.0)
@(
    ('ATOM      1  C1  MOL A   1      {0,8:F3} {1,8:F3}   0.000   0.2500 1.0000' -f -1.5,(-0.5*$root3)),
    ('ATOM      2  C2  MOL A   1      {0,8:F3} {1,8:F3}   0.000  -0.2500 1.0000' -f 1.5,(-0.5*$root3)),
    ('ATOM      3  C3  MOL A   1       0.000 {0,8:F3}   0.000   0.0000 1.0000' -f $root3),
    'END'
) | Set-Content -LiteralPath $pqr -Encoding ASCII

$argPath = Write-FmmCliArgFile -RunDir $WorkRoot -RepoRoot $RepoRoot `
    -PqrPath $pqr -Backend tabi -Mode direct -FmmMode off `
    -Formulation tabi2-nodepatch -SurfacePoints 8 -MaxIter 200 -Restart 64 `
    -Tolerance 1e-4 -Threads 2 -SaltM 0.15 -InternalDielectric 1.0 `
    -ExternalDielectric 80.0 -TemperatureK 298.15 -Preconditioner diagonal `
    -NearIntegration analytic_singular -PrimarySurface analytic-ses `
    -SesGridSpacing 0.5

$run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $WorkRoot `
    -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
$profilePath = Join-Path $WorkRoot 'pbe_tabi_profile.csv'
$auditPath = Join-Path $WorkRoot 'pbe_tabi_formulation_audit.csv'
$profile = if (Test-Path -LiteralPath $profilePath) {
    @(Import-Csv -LiteralPath $profilePath) | Select-Object -Last 1
} else { $null }
$audit = if (Test-Path -LiteralPath $auditPath) {
    @(Import-Csv -LiteralPath $auditPath) | Select-Object -Last 1
} else { $null }
$status = if ($run.timed_out) { 'timeout' } elseif ($run.exit_code -ne 0) {
    'failed'
} elseif ($null -eq $profile -or $null -eq $audit) { 'profile_missing' } else {
    [string]$profile.status
}
$contractPass = $status -eq 'ok' -and
    [string]$profile.surface_generator -eq 'titan_analytic_patch_ray_ses' -and
    [string]$audit.formulation -eq 'tabi2-nodepatch' -and
    [string]$audit.discretization -eq 'node_patch'
$output = [pscustomobject]@{
    status = $status
    contract_pass = $contractPass
    formulation = if ($null -ne $audit) { [string]$audit.formulation } else { '' }
    discretization = if ($null -ne $audit) { [string]$audit.discretization } else { '' }
    surface_generator = if ($null -ne $profile) { [string]$profile.surface_generator } else { '' }
    panels = if ($null -ne $profile) { [int64]$profile.panels } else { 0 }
    iterations = if ($null -ne $profile) { [int]$profile.iterations } else { 0 }
    true_residual = if ($null -ne $profile) { [double]$profile.true_residual_rel } else { [double]::NaN }
    implementation_scope = 'direct_report_path_no_vertex_patch_fmm_preconditioner'
    ilu_comparison_available = $false
    run_dir = $WorkRoot
}
$output | Export-Csv -LiteralPath (Join-Path $WorkRoot 'nodepatch_analytic_ses_control.csv') `
    -NoTypeInformation -Encoding UTF8
if (-not $contractPass) {
    throw "TABI2 node-patch/analytic-SES control failed: status=$status"
}
Write-Host ("nodepatch_analytic_ses_control=pass panels={0} iterations={1}; ILU comparison unavailable by contract" -f `
    $output.panels,$output.iterations)
