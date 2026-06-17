param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 300,
    [int]$Threads = 2,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z2b')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z2b_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z2b_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z2b_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z2b_near_policy_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z2B_DUFFY_GATE_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$TitanDir = Split-Path -Parent $TitanExe
$LibFfRoot = Join-Path $RepoRoot 'lib_ff'
if (-not (Test-Path -LiteralPath $LibFfRoot)) {
    $LibFfRoot = Join-Path $TitanDir 'lib_ff'
}
if (-not (Test-Path -LiteralPath $LibFfRoot)) {
    throw "Missing lib_ff beside repo or titan executable: $LibFfRoot"
}
$LibFfRoot = (Resolve-Path -LiteralPath $LibFfRoot).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Z2B work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z2bDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or $text -eq 'nan') {
        return $Default
    }
    return [double]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Convert-Z2bInt64([object]$Value, [int64]$Default = 0) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $Default }
    return [int64]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Write-Z2bArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [ValidateSet('direct','treecode')][string]$Mode,
    [ValidateSet('off','single-matvec-analytic')][string]$FmmMode,
    [string]$Preconditioner,
    [string]$NearIntegration,
    [string]$QuadBuffer,
    [int]$SurfacePoints,
    [int]$MaxIter,
    [int]$Restart,
    [double]$SaltM,
    [int]$AuditSampleCount
) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdbPath
    $libPath = Convert-FmmCliToTitanPath($LibFfRoot.TrimEnd('\') + '\')
    $format = {
        param([double]$Value)
        return $Value.ToString('G15', [System.Globalization.CultureInfo]::InvariantCulture)
    }
    $argLines = @(
        "-f=$(Convert-FmmCliToTitanPath $pdbPath)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-FmmCliToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        "-pbe-max-iter=$MaxIter",
        '-pbe-rms=1e-6',
        '-pbe-dphi=1e-6',
        '-pbe-int-diel=1',
        '-pbe-ext-diel=80',
        "-pbe-salt=$(& $format $SaltM)",
        '-pbe-temp=298.15',
        '-pbe-timing=detail',
        "-nt=$Threads",
        "-pbe-tabi-sum=$Mode",
        '-pbe-tabi-formulation=jcp2013-treecode',
        "-pbe-tabi-fmm-m2l=$FmmMode",
        "-pbe-tabi-surface-points=$SurfacePoints",
        '-pbe-tabi-tree-leaf-size=64',
        '-pbe-tabi-opening-theta=0.6',
        "-pbe-tabi-near-integration=$NearIntegration",
        "-pbe-tabi-quad-buffer=$QuadBuffer",
        "-pbe-tabi-precond=$Preconditioner",
        '-pbe-tabi-precond-ilu-neighbor-cap=64',
        "-pbe-tabi-gmres-restart=$Restart",
        '-pbe-tabi-true-residual-probe=restart',
        '-pbe-tabi-true-residual-probe-scope=sampled',
        '-pbe-tabi-true-residual-sample-count=256'
    )
    if ($AuditSampleCount -gt 0) {
        $argLines += '-pbe-tabi-matvec-audit=sampled'
        $argLines += "-pbe-tabi-matvec-audit-sample-count=$AuditSampleCount"
    } else {
        $argLines += '-pbe-tabi-matvec-audit=off'
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    $badOptions = @($argLines | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
    })
    if ($badOptions.Count -ne 0) {
        throw "Generated arg.in contains non equals-syntax options: $($badOptions -join '; ')"
    }
    return $argPath
}

function Read-Z2bResult([string]$RunDir, [string]$ExpectedFmmMode) {
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $energyPath)) {
        throw "Missing pbe_energy.csv in $RunDir"
    }
    if (-not (Test-Path -LiteralPath $profilePath)) {
        throw "Missing pbe_tabi_profile.csv in $RunDir"
    }
    $energy = @(Import-Csv -LiteralPath $energyPath)[0]
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    $mode = [string]$profile.fmm_m2l_mode
    if ($ExpectedFmmMode -eq 'single-matvec-analytic' -and
        $mode -ne 'single-matvec-analytic') {
        throw "Expected FMM/M2L mode single-matvec-analytic, got $mode"
    }
    return [pscustomobject]@{
        status = if ($energy.status -eq 'ok' -and $profile.status -eq 'ok') { 'ok' } else { "energy=$($energy.status);profile=$($profile.status)" }
        energy_kj_mol = Convert-Z2bDouble $energy.E_rxn_kJmol
        panels = Convert-Z2bInt64 $profile.panels
        iterations = Convert-Z2bInt64 $profile.iterations
        true_residual = Convert-Z2bDouble $profile.true_residual_rel
        near_integration = [string]$profile.near_integration
        quad_buffer = [string]$profile.quad_buffer
        near_kernel_policy = [string]$profile.near_kernel_policy
        near_interactions = Convert-Z2bInt64 $profile.near_interactions
        near_self_pairs = Convert-Z2bInt64 $profile.near_self_pairs
        near_one_ring_pairs = Convert-Z2bInt64 $profile.near_one_ring_pairs
        near_duffy_pairs = Convert-Z2bInt64 $profile.near_duffy_pairs
        near_centroid_pairs = Convert-Z2bInt64 $profile.near_centroid_pairs
        near_fallbacks = Convert-Z2bInt64 $profile.near_fallbacks
        near_integration_seconds = Convert-Z2bDouble $profile.near_integration_seconds 0.0
        fmm_m2l_mode = $mode
        fmm_near_pairs = Convert-Z2bInt64 $profile.fmm_m2l_near_pairs
        fmm_far_pairs = Convert-Z2bInt64 $profile.fmm_m2l_far_pairs
        fmm_finite_outputs = [string]$profile.fmm_m2l_finite_outputs
        tree_vs_direct_matvec_rel = Convert-Z2bDouble $profile.tree_vs_direct_matvec_rel
        matvec_audit_rows = Convert-Z2bInt64 $profile.matvec_audit_rows
    }
}

function Invoke-Z2bCase(
    [string]$Name,
    [string]$PqrPath,
    [ValidateSet('direct','treecode')][string]$Mode,
    [ValidateSet('off','single-matvec-analytic')][string]$FmmMode,
    [string]$Preconditioner,
    [string]$NearIntegration,
    [string]$QuadBuffer,
    [int]$SurfacePoints,
    [int]$MaxIter,
    [int]$Restart,
    [double]$SaltM,
    [int]$AuditSampleCount
) {
    $runDir = Join-Path $workRootFull $Name
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $argPath = Write-Z2bArgFile -RunDir $runDir -PqrPath $PqrPath `
        -Mode $Mode -FmmMode $FmmMode -Preconditioner $Preconditioner `
        -NearIntegration $NearIntegration -QuadBuffer $QuadBuffer `
        -SurfacePoints $SurfacePoints -MaxIter $MaxIter -Restart $Restart `
        -SaltM $SaltM -AuditSampleCount $AuditSampleCount
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        throw "Z2B case $Name failed with exit code $($run.exit_code). See $runDir"
    }
    $result = Read-Z2bResult -RunDir $runDir -ExpectedFmmMode $FmmMode
    return [pscustomobject]@{
        condition = $Name; exit_code = $run.exit_code; timed_out = $run.timed_out
        status = $result.status
        near_integration = $result.near_integration
        quad_buffer = $result.quad_buffer
        near_kernel_policy = $result.near_kernel_policy
        fmm_m2l_mode = $result.fmm_m2l_mode
        fmm_near_pairs = $result.fmm_near_pairs
        fmm_far_pairs = $result.fmm_far_pairs
        fmm_finite_outputs = $result.fmm_finite_outputs
        tree_vs_direct_matvec_rel = $result.tree_vs_direct_matvec_rel
        matvec_audit_rows = $result.matvec_audit_rows
        energy_kj_mol = $result.energy_kj_mol
        panels = $result.panels
        iterations = $result.iterations
        true_residual = $result.true_residual
        near_interactions = $result.near_interactions
        near_self_pairs = $result.near_self_pairs
        near_one_ring_pairs = $result.near_one_ring_pairs
        near_duffy_pairs = $result.near_duffy_pairs
        near_centroid_pairs = $result.near_centroid_pairs
        near_fallbacks = $result.near_fallbacks
        near_integration_seconds = $result.near_integration_seconds
        wall_seconds = $run.wall_seconds
        peak_working_set_mb = $run.peak_working_set_mb
    }
}

function Write-Z2bSvg([object[]]$Rows, [string]$Path) {
    $points = @($Rows | Where-Object { $_.condition -match '_fmm_' })
    $maxY = [Math]::Max(1.0, ($points | Measure-Object -Property wall_seconds -Maximum).Maximum)
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="760" height="300" viewBox="0 0 760 300">') | Out-Null
    $svg.Add('<rect width="760" height="300" fill="white"/>') | Out-Null
    $svg.Add('<text x="24" y="28" font-family="Arial" font-size="16">Scope Z2B FMM Duffy gate</text>') | Out-Null
    $x = 60
    foreach ($row in $points) {
        $h = [Math]::Max(2.0, 210.0 * ([double]$row.wall_seconds / $maxY))
        $y = 250.0 - $h
        $svg.Add(('<rect x="{0}" y="{1:F1}" width="80" height="{2:F1}" fill="#d55e00"/>' -f $x, $y, $h)) | Out-Null
        $svg.Add(('<text x="{0}" y="270" font-family="Arial" font-size="10" transform="rotate(30 {0},270)">{1}</text>' -f $x, $row.condition)) | Out-Null
        $x += 120
    }
    $svg.Add('<text x="24" y="292" font-family="Arial" font-size="11">Bar height = wall seconds for FMM Duffy rows</text>') | Out-Null
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$fixtureDir = Join-Path $workRootFull 'fixtures'
New-Item -ItemType Directory -Path $fixtureDir -Force | Out-Null
$bornPqr = Join-Path $fixtureDir 'born_centered.pqr'
$kirkwoodPqr = Join-Path $fixtureDir 'kirkwood_offcenter.pqr'
Write-FagbiBornPqr -Path $bornPqr -Charge 1.0 -RadiusA 3.0
Write-FagbiKirkwoodPqr -Path $kirkwoodPqr -SphereRadiusA 10.0 -Charge 1.0 -OffsetA 3.0

$rows = @()
foreach ($case in @(
    @{ name='born'; pqr=$bornPqr; max_iter=80 },
    @{ name='kirkwood'; pqr=$kirkwoodPqr; max_iter=120 }
)) {
    $rows += Invoke-Z2bCase -Name "$($case.name)_direct_duffy_adaptive" `
        -PqrPath $case.pqr -Mode direct -FmmMode off -Preconditioner diagonal `
        -NearIntegration duffy_adaptive -QuadBuffer adaptive -SurfacePoints 8 `
        -MaxIter $case.max_iter -Restart 64 -SaltM 0.15 -AuditSampleCount 0
    $rows += Invoke-Z2bCase -Name "$($case.name)_fmm_duffy_adaptive" `
        -PqrPath $case.pqr -Mode treecode -FmmMode single-matvec-analytic `
        -Preconditioner fmm-nearfield-ilu0 -NearIntegration duffy_adaptive `
        -QuadBuffer adaptive -SurfacePoints 8 -MaxIter $case.max_iter `
        -Restart 64 -SaltM 0.15 -AuditSampleCount 128
}

if (-not $Quick) {
    $sampleDir = Join-Path $workRootFull 'samples'
    New-Item -ItemType Directory -Path $sampleDir -Force | Out-Null
    $fixture = Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot `
        -RunDir $sampleDir
    if ([string]::IsNullOrWhiteSpace($fixture)) {
        $fixture = Resolve-FmmCliFixture -Fixture 'research-medium' -RepoRoot $RepoRoot
    }
    $proteinPqr = Join-Path $sampleDir 'z2b_protein_cap128.pqr'
    $sampleInfo = Write-FmmCliSampledPqr -InputPath $fixture -OutputPath $proteinPqr -Cap 128
    if ($sampleInfo.sampled_atoms -lt 64) {
        throw "Z2B protein smoke expected at least 64 sampled atoms, got $($sampleInfo.sampled_atoms)"
    }
    $rows += Invoke-Z2bCase -Name 'protein_fmm_duffy_adaptive' `
        -PqrPath $proteinPqr -Mode treecode -FmmMode single-matvec-analytic `
        -Preconditioner fmm-nearfield-ilu0 -NearIntegration duffy_adaptive `
        -QuadBuffer adaptive -SurfacePoints 8 -MaxIter 160 -Restart 64 `
        -SaltM 0.15 -AuditSampleCount 128
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation

$failures = @()
foreach ($row in $rows) {
    if ($row.status -ne 'ok') {
        $failures += "$($row.condition) status=$($row.status)"
    }
    if ($row.near_integration -ne 'duffy_adaptive') {
        $failures += "$($row.condition) near mode drifted to $($row.near_integration)"
    }
    if ($row.quad_buffer -ne 'adaptive') {
        $failures += "$($row.condition) quad buffer drifted to $($row.quad_buffer)"
    }
    if ($row.near_kernel_policy -ne 'one_ring_cpv') {
        $failures += "$($row.condition) near policy expected one_ring_cpv, got $($row.near_kernel_policy)"
    }
    if ($row.condition -match '_fmm_' -and $row.fmm_m2l_mode -ne 'single-matvec-analytic') {
        $failures += "$($row.condition) missing FMM/M2L telemetry"
    }
    if ($row.condition -match '_fmm_' -and
        ($row.fmm_near_pairs -le 0 -or $row.fmm_finite_outputs -ne 'true')) {
        $failures += "$($row.condition) FMM counters invalid"
    }
    if ($row.near_duffy_pairs -le 0 -or $row.near_fallbacks -ne 0) {
        $failures += "$($row.condition) Duffy near counters invalid: duffy=$($row.near_duffy_pairs) fallback=$($row.near_fallbacks)"
    }
    if ($row.condition -match '_fmm_' -and $row.matvec_audit_rows -gt 0 -and
        $row.tree_vs_direct_matvec_rel -gt 1.0e-2) {
        $failures += "$($row.condition) audit too high: $($row.tree_vs_direct_matvec_rel)"
    }
}

$decision = if ($failures.Count -eq 0) {
    'z2b_duffy_explicit_pass'
} elseif ($rows.Count -gt 0 -and
    (@($rows | Where-Object { $_.status -eq 'ok' }).Count -gt 0)) {
    'z2b_duffy_not_promoted'
} else {
    'z2b_duffy_regression'
}

[pscustomobject]@{
    decision = $decision
    failures = ($failures -join '; ')
    rows = $rows.Count
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'scope: Z2B',
    "decision: $decision",
    'near_policy: explicit duffy_adaptive enabled for FMM/M2L',
    'auto_policy: unchanged',
    'grid_default: unchanged'
)
if ($failures.Count -gt 0) {
    $summary += "failures: $($failures -join '; ')"
}
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII
Write-Z2bSvg -Rows $rows -Path $PlotPath

$report = @(
    '# TITAN FMM/M2L Scope Z2B Duffy Gate Report',
    '',
    'Z2B verifies that explicit duffy_adaptive near integration can run through the FMM/M2L near-source path while keeping auto/default policies unchanged.',
    '',
    "Decision: $decision",
    '',
    '## Rows',
    ''
)
foreach ($row in $rows) {
    $report += ('- {0}: status={1}, near={2}, quad={3}, policy={4}, fmm={5}, duffy_pairs={6}, fallback={7}, audit={8}' -f `
        $row.condition, $row.status, $row.near_integration, $row.quad_buffer, `
        $row.near_kernel_policy, $row.fmm_m2l_mode, $row.near_duffy_pairs, `
        $row.near_fallbacks, $row.tree_vs_direct_matvec_rel)
}
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

if ($failures.Count -gt 0) {
    throw "Z2B Duffy gate failed: $($failures -join '; ')"
}

Write-Host 'PASS pbe_fagbi_scope_z2b_fmm_duffy_gate'
