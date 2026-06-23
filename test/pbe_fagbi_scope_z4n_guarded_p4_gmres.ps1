param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4n'),
    [string]$RootReport = '',
    [int]$CaseTimeoutSeconds = 3600,
    [switch]$Quick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
if (-not [IO.Path]::IsPathRooted($FixtureRoot)) {
    $FixtureRoot = Join-Path $RepoRoot $FixtureRoot
}
$FixtureRoot = [IO.Path]::GetFullPath($FixtureRoot)
if (-not [IO.Path]::IsPathRooted($WorkRoot)) {
    $WorkRoot = Join-Path $RepoRoot $WorkRoot
}
$WorkRoot = [IO.Path]::GetFullPath($WorkRoot)
if (-not [string]::IsNullOrWhiteSpace($RootReport) -and
    -not [IO.Path]::IsPathRooted($RootReport)) {
    $RootReport = Join-Path $RepoRoot $RootReport
}

$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4n_matrix.csv'
$convergenceCsv = Join-Path $WorkRoot 'fagbi_scope_z4n_convergence.csv'
$auditCsv = Join-Path $WorkRoot 'fagbi_scope_z4n_operator_audit.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4n_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4n_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4n_gmres_plot.svg'
$rootReport = if ([string]::IsNullOrWhiteSpace($RootReport)) {
    Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4N_GUARDED_GMRES_REPORT.txt'
} else { $RootReport }
$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) { throw "Missing fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)

$RuntimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $RuntimeRepoRoot 'lib_ff'))) {
    $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($commonGitDir)) {
        $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir))
    }
    $RuntimeRepoRoot = Split-Path -Parent $commonGitDir
}

$generatedPqr = Join-Path $WorkRoot 'generated_direct_reference.pqr'
@(
    'ATOM      1  C   GLY A   1      -6.000   0.000   0.000  1.0000  1.7000',
    'ATOM      2  O   GLY A   1      -2.000   0.000   0.000 -1.0000  1.5000',
    'ATOM      3  N   GLY A   2       2.000   0.000   0.000  1.0000  1.5500',
    'ATOM      4  O   GLY A   2       6.000   0.000   0.000 -1.0000  1.5000'
) | Set-Content -LiteralPath $generatedPqr -Encoding ASCII

$manifest = @(Import-Csv -LiteralPath $manifestPath)
function Get-Z4nFixture([string]$Name) {
    if ($Name -eq 'generated') {
        return [pscustomobject]@{fixture='generated';protein_name='generated direct reference'
            pdb_id='generated';pqr=$generatedPqr}
    }
    $row = @($manifest | Where-Object fixture -eq $Name)
    if ($row.Count -ne 1) { throw "Expected one public fixture row for $Name" }
    return $row[0]
}

function Get-Z4nDouble([object]$Value,[double]$Default=[double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Invoke-Z4nCase([pscustomobject]$Fixture,[string]$Variant) {
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_$Variant"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $isP4 = $Variant -eq 'p4_gmres'
    $pqr = if ($Fixture.fixture -eq 'generated') { $Fixture.pqr } else {
        Join-Path $FixtureRoot $Fixture.pqr
    }
    $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RuntimeRepoRoot `
        -PqrPath $pqr -Backend tabi -Mode treecode `
        -FmmMode single-matvec-analytic -SurfacePoints 8 -MaxIter 720 `
        -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount 256 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular `
        -TaylorOrder 4 -FmmCartesianAudit off `
        -FmmLocal $(if ($isP4) { 'gmres-p4' } else { 'off' }) `
        -FmmLocalRefine $(if ($isP4) { 'operator-tail-static' } else { 'off' }) `
        -FmmLocalRefineTolerance 0.002 `
        -FmmLocalCache $(if ($isP4) { 'green-jet-compressed' } else { 'off' }) `
        -FmmLocalCacheBudgetMb 256 `
        -FmmLocalFallback $(if ($isP4) { 'packed-avx2' } else { 'scalar' }) `
        -FmmLocalPlanBuilder $(if ($isP4) { 'metric-cache' } else { 'legacy-rescan' })
    $bare = @(Get-Content -LiteralPath $arg | Where-Object { $_ -match '^-[A-Za-z0-9-]+$' })
    if ($bare.Count) { throw "Generated arg.in contains bare options: $($bare -join ', ')" }
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) {
        return [pscustomobject]@{fixture=$Fixture.fixture;protein_name=$Fixture.protein_name
            pdb_id=$Fixture.pdb_id;variant=$Variant;status='timeout'
            wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb
            panels=0;iterations=0;internal_true_residual=[double]::NaN
            final_residual=[double]::NaN;energy_kj_mol=[double]::NaN
            audit_rms=[double]::NaN;matvec_seconds_avg=[double]::NaN
            fmm_local_gmres_active='false';prepare_seconds=0.0
            plan_build_count=0;cache_build_count=0;gmres_matvec_calls=0
            gmres_apply_seconds_total=0.0;gmres_apply_seconds_mean=0.0
            moment_seconds_total=0.0;m2l_seconds_total=0.0
            fallback_seconds_total=0.0;l2p_seconds_total=0.0
            sampled_direct_true_residual=[double]::NaN;operator_finite='false'
            residual_history='';p4_vs_degree1_matvec_ratio=[double]::NaN
            run_dir=$runDir}
    }
    if ($run.exit_code -ne 0) {
        $log = (Get-Content -Raw (Join-Path $runDir 'stdout.log')) +
            (Get-Content -Raw (Join-Path $runDir 'stderr.log'))
        $failure = if ($log -match 'cache_budget_exceeded') {
            'cache_budget_exceeded'
        } elseif ($log -match 'non_finite') { 'non_finite' } else { 'failed' }
        return [pscustomobject]@{fixture=$Fixture.fixture;protein_name=$Fixture.protein_name
            pdb_id=$Fixture.pdb_id;variant=$Variant;status=$failure
            wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb
            panels=0;iterations=0;internal_true_residual=[double]::NaN
            final_residual=[double]::NaN;energy_kj_mol=[double]::NaN
            audit_rms=[double]::NaN;matvec_seconds_avg=[double]::NaN
            fmm_local_gmres_active='false';prepare_seconds=0.0
            plan_build_count=0;cache_build_count=0;gmres_matvec_calls=0
            gmres_apply_seconds_total=0.0;gmres_apply_seconds_mean=0.0
            moment_seconds_total=0.0;m2l_seconds_total=0.0
            fallback_seconds_total=0.0;l2p_seconds_total=0.0
            sampled_direct_true_residual=[double]::NaN;operator_finite='false'
            residual_history='';p4_vs_degree1_matvec_ratio=[double]::NaN
            failure_log=$log;run_dir=$runDir}
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
    return [pscustomobject]@{
        fixture=$Fixture.fixture;protein_name=$Fixture.protein_name;pdb_id=$Fixture.pdb_id
        variant=$Variant;status=$case.status;panels=$case.panels
        iterations=$case.iterations;internal_true_residual=$case.true_residual
        final_residual=$case.final_residual;energy_kj_mol=$case.energy_kj_mol
        audit_rms=$case.tree_vs_direct_matvec_rel;matvec_seconds_avg=$case.matvec_seconds_avg
        wall_seconds=$run.wall_seconds;peak_rss_mb=$run.peak_working_set_mb
        fmm_local_gmres_active=[string]$profile.fmm_local_gmres_active
        prepare_seconds=Get-Z4nDouble $profile.fmm_local_prepare_seconds 0.0
        plan_build_count=[int64]$profile.fmm_local_plan_build_count
        cache_build_count=[int64]$profile.fmm_local_cache_build_count
        gmres_matvec_calls=[int64]$profile.fmm_local_gmres_matvec_calls
        gmres_apply_seconds_total=Get-Z4nDouble `
            $profile.fmm_local_gmres_apply_seconds_total 0.0
        gmres_apply_seconds_mean=Get-Z4nDouble `
            $profile.fmm_local_gmres_apply_seconds_mean 0.0
        moment_seconds_total=Get-Z4nDouble $profile.fmm_local_moment_seconds_total 0.0
        m2l_seconds_total=Get-Z4nDouble $profile.fmm_local_m2l_seconds_total 0.0
        fallback_seconds_total=Get-Z4nDouble `
            $profile.fmm_local_fallback_seconds_total 0.0
        l2p_seconds_total=Get-Z4nDouble $profile.fmm_local_l2p_seconds_total 0.0
        sampled_direct_true_residual=Get-Z4nDouble `
            $profile.fmm_local_sampled_direct_true_residual 0.0
        operator_finite=[string]$profile.fmm_local_operator_finite
        residual_history=$case.residual_history
        p4_vs_degree1_matvec_ratio=[double]::NaN;run_dir=$runDir
    }
}

$fixtureNames = if ($Quick) { @('generated','1sh1') } else {
    @('barnase','1a63','1sh1')
}
$rows = New-Object System.Collections.Generic.List[object]
foreach ($fixtureName in $fixtureNames) {
    $fixture = Get-Z4nFixture $fixtureName
    $rows.Add((Invoke-Z4nCase $fixture 'degree1_baseline')) | Out-Null
    $rows.Add((Invoke-Z4nCase $fixture 'p4_gmres')) | Out-Null
}
$matrix = @($rows.ToArray())

foreach ($fixtureName in $fixtureNames) {
    $baseline = @($matrix | Where-Object {
        $_.fixture -eq $fixtureName -and $_.variant -eq 'degree1_baseline'
    }) | Select-Object -First 1
    $p4 = @($matrix | Where-Object {
        $_.fixture -eq $fixtureName -and $_.variant -eq 'p4_gmres'
    }) | Select-Object -First 1
    if ($null -ne $p4 -and $p4.PSObject.Properties.Name -contains 'matvec_seconds_avg') {
        $ratio = if ($null -ne $baseline -and
            $baseline.PSObject.Properties.Name -contains 'matvec_seconds_avg' -and
            [double]$baseline.matvec_seconds_avg -gt 0) {
            [double]$p4.matvec_seconds_avg / [double]$baseline.matvec_seconds_avg
        } else { [double]::NaN }
        $p4.p4_vs_degree1_matvec_ratio = $ratio
    }
}

$trace = New-Object System.Collections.Generic.List[object]
foreach ($row in $matrix) {
    if (-not ($row.PSObject.Properties.Name -contains 'residual_history')) { continue }
    $values = @(([string]$row.residual_history) -split ';' | Where-Object { $_ })
    for ($iteration = 0; $iteration -lt $values.Count; ++$iteration) {
        $trace.Add([pscustomobject]@{fixture=$row.fixture;variant=$row.variant
            iteration=$iteration;residual=Get-Z4nDouble $values[$iteration]}) | Out-Null
    }
}

$requiredP4 = @($matrix | Where-Object variant -eq 'p4_gmres')
$decision = 'z4n_guarded_gmres_pass'
if (@($requiredP4 | Where-Object status -eq 'timeout').Count) {
    $decision = 'z4n_runtime_limited'
} elseif (@($requiredP4 | Where-Object {
        $_.status -eq 'cache_budget_exceeded' -or
        ($_.PSObject.Properties.Name -contains 'peak_rss_mb' -and
         [double]$_.peak_rss_mb -ge 512)
    }).Count) {
    $decision = 'z4n_memory_limited'
} elseif (@($requiredP4 | Where-Object {
        $_.status -ne 'ok' -or [double]$_.internal_true_residual -gt 1e-6 -or
        -not (Test-FmmCliFinite ([double]$_.energy_kj_mol))
    }).Count) {
    $decision = 'z4n_gmres_convergence_gap'
} elseif (@($requiredP4 | Where-Object {
        [double]$_.audit_rms -gt 8e-3 -or
        [double]$_.sampled_direct_true_residual -gt 1e-2 -or
        $_.operator_finite -ne 'true'
    }).Count) {
    $decision = 'z4n_operator_accuracy_gap'
} elseif (@($requiredP4 | Where-Object {
        [int64]$_.plan_build_count -ne 1 -or
        [int64]$_.cache_build_count -ne 1 -or
        [int64]$_.gmres_matvec_calls -le 0 -or
        $_.fmm_local_gmres_active -ne 'true'
    }).Count) {
    $decision = 'z4n_cache_lifecycle_gap'
} elseif (@($requiredP4 | Where-Object {
        $_.PSObject.Properties.Name -contains 'p4_vs_degree1_matvec_ratio' -and
        [double]$_.p4_vs_degree1_matvec_ratio -gt 2.0
    }).Count) {
    $decision = 'z4n_accuracy_pass_runtime_gap'
}

$matrix | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$trace.ToArray() | Export-Csv -LiteralPath $convergenceCsv -NoTypeInformation
$requiredP4 | Select-Object fixture,audit_rms,sampled_direct_true_residual,
    internal_true_residual,operator_finite,plan_build_count,cache_build_count |
    Export-Csv -LiteralPath $auditCsv -NoTypeInformation
[pscustomobject]@{decision=$decision;quick=[bool]$Quick;rows=$matrix.Count
    required_p4_rows=$requiredP4.Count} |
    Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
@('scope: Z4N guarded p4 GMRES integration',"decision: $decision",
  "quick: $([bool]$Quick)","rows: $($matrix.Count)",
  'gmres_p4: hidden explicit only','tabi_treecode_default: degree1',
  'production_default: grid') | Set-Content -LiteralPath $summaryPath -Encoding ASCII

$svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="1080" height="640" viewBox="0 0 1080 640">',
    '<rect width="1080" height="640" fill="white"/>',
    '<text x="35" y="38" font-family="Arial" font-size="20">Guarded p4 GMRES convergence and operator audit</text>')
$y = 78
foreach ($row in $matrix) {
    $audit = if ($row.PSObject.Properties.Name -contains 'audit_rms') { $row.audit_rms } else { 'n/a' }
    $svg += ('<text x="35" y="{0}" font-family="Arial" font-size="13">{1} {2}: status={3}, residual={4}, audit={5}, matvec={6}s, RSS={7} MB</text>' -f
        $y,$row.fixture,$row.variant,$row.status,$row.internal_true_residual,
        $audit,$row.matvec_seconds_avg,$row.peak_rss_mb)
    $y += 34
}
$svg += '</svg>'
$svg | Set-Content -LiteralPath $plotPath -Encoding ASCII

$report = @('# TITAN FMM/M2L Scope Z4N Guarded p4 GMRES Report','',
    "Decision: $decision",'',
    'The p4 operator is explicit-only. Static refinement, compressed Green-jet cache, packed AVX2 fallback, and FMM-native ILU(0) are prepared once before GMRES.','')
foreach ($row in $matrix) {
    $report += ('- {0} {1}: status={2}, internal_residual={3}, direct_sampled_residual={4}, audit={5}, matvec={6}s, wall={7}s, RSS={8} MB' -f
        $row.fixture,$row.variant,$row.status,$row.internal_true_residual,
        $row.sampled_direct_true_residual,$row.audit_rms,$row.matvec_seconds_avg,
        $row.wall_seconds,$row.peak_rss_mb)
}
$report += @('','Grid remains the production PBE default. TABI treecode remains degree-1 unless gmres-p4 is explicitly selected.')
$report | Set-Content -LiteralPath $rootReport -Encoding ASCII

if ($decision -ne 'z4n_guarded_gmres_pass' -and -not $AllowFailedDecision) {
    throw "Scope Z4N gate did not pass: $decision"
}
Write-Host "PASS pbe_fagbi_scope_z4n_guarded_p4_gmres decision=$decision"
