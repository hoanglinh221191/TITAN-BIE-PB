param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4g'),
    [int]$CaseTimeoutSeconds = 2400,
    [int]$Repeats = 3,
    [switch]$Quick
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4g_matrix.csv'
$RepeatCsv = Join-Path $WorkRoot 'fagbi_scope_z4g_repeatability.csv'
$TraceCsv = Join-Path $WorkRoot 'fagbi_scope_z4g_convergence_trace.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4g_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4g_summary.dat'
$RuntimePlot = Join-Path $WorkRoot 'fagbi_scope_z4g_runtime_memory.svg'
$ConvergencePlot = Join-Path $WorkRoot 'fagbi_scope_z4g_preconditioner_convergence.svg'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4G_EXTERNAL_MANUSCRIPT_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $manifestPath)) {
    throw "Missing public fixture manifest. Run test\pbe_fagbi_scope_z4g_prepare_public_fixtures.ps1 -RefreshFixtures"
}
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$manifest = @(Import-Csv -LiteralPath $manifestPath)
$wanted = if ($Quick) { @('1sh1') } else { @('barnase','1a63','1sh1') }
$repeatCount = if ($Quick) { 1 } else { [Math]::Max(1, $Repeats) }

function Get-Z4gInputAtomCount([string]$Path) {
    return @(Get-Content -LiteralPath $Path | Where-Object {
            $_.StartsWith('ATOM') -or $_.StartsWith('HETATM')
        }).Count
}

function Get-Z4gProfileValue([string]$RunDir, [string]$Name, [string]$Default = '') {
    $path = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $path)) { return $Default }
    $profile = @(Import-Csv -LiteralPath $path)[0]
    if ($profile.PSObject.Properties.Name -notcontains $Name) { return $Default }
    return [string]$profile.$Name
}

function Invoke-Z4gCase(
    [pscustomobject]$Fixture,
    [string]$Backend,
    [int]$Repeat,
    [string]$Preconditioner = 'fmm-nearfield-ilu0'
) {
    $pqr = Join-Path $FixtureRoot $Fixture.pqr
    $inputAtoms = Get-Z4gInputAtomCount $pqr
    if ($inputAtoms -ne [int]$Fixture.atom_count) {
        throw "Full-protein atom count drift for $($Fixture.fixture)"
    }
    $condition = "{0}_{1}_r{2}" -f $Fixture.fixture, $Backend, $Repeat
    $runDir = Join-Path $WorkRoot $condition
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    $isGrid = $Backend -eq 'grid_context'
    $isDirect = $Backend -eq 'tabi_direct'
    $isFmm = -not $isGrid -and -not $isDirect
    $argPath = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot -PqrPath $pqr `
        -Backend $(if ($isGrid) { 'grid' } else { 'tabi' }) `
        -Mode $(if ($isDirect) { 'direct' } else { 'treecode' }) `
        -FmmMode $(if ($isFmm) { 'single-matvec-analytic' } else { 'off' }) `
        -SurfacePoints 8 -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 `
        -SaltM 0.15 -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount $(if ($isFmm) { 512 } else { 0 }) `
        -Preconditioner $(if ($isFmm) { $Preconditioner } else { 'diagonal' }) `
        -PreconditionerBlockPanelCap 32 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular -TaylorOrder 4
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $base = [ordered]@{
        fixture=$Fixture.fixture; protein_name=$Fixture.protein_name; pdb_id=$Fixture.pdb_id;
        backend=$Backend; repeat=$Repeat; atom_count=$inputAtoms; input_atoms=$inputAtoms; analyzed_atoms=$inputAtoms;
        panels=0; status='failed'; failure_class='exit_code'; energy_kj_mol=[double]::NaN;
        true_residual=[double]::NaN; audit_rms=[double]::NaN; audit_rows=0; iterations=0;
        wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb;
        matvec_seconds=[double]::NaN; preconditioner_requested=$(if ($isFmm) { $Preconditioner } else { 'diagonal' });
        preconditioner_effective=''; fmm_mode=$(if ($isFmm) { 'single-matvec-analytic' } else { 'off' });
        near_pairs=0; far_pairs=0; residual_history=''; arg_file=$argPath
    }
    if ($run.timed_out) {
        $base.status='timeout'; $base.failure_class='timeout_not_verified'
        return [pscustomobject]$base
    }
    if ($run.exit_code -ne 0) {
        $base.failure_class="exit_code_$($run.exit_code)"
        return [pscustomobject]$base
    }
    if ($isGrid) {
        $energyPath = Join-Path $runDir 'pbe_energy.csv'
        if (-not (Test-Path -LiteralPath $energyPath)) {
            $base.failure_class='missing_grid_energy'
            return [pscustomobject]$base
        }
        $energy = @(Import-Csv -LiteralPath $energyPath)[0]
        $base.status=[string]$energy.status
        $base.failure_class=if ($base.status -eq 'ok') { 'context_ok' } else { 'context_not_converged' }
        $base.energy_kj_mol=Convert-FmmCliDouble $energy.E_rxn_kJmol
        return [pscustomobject]$base
    }
    try {
        $result = Read-FmmCliCaseResult -RunDir $runDir `
            -ExpectedFmmMode $(if ($isFmm) { 'single-matvec-analytic' } else { 'off' }) `
            -RequireAudit $isFmm
        $base.panels=$result.panels; $base.status=$result.status
        $base.energy_kj_mol=$result.energy_kj_mol; $base.true_residual=$result.true_residual
        $base.audit_rms=$result.tree_vs_direct_matvec_rel; $base.audit_rows=$result.matvec_audit_rows
        $base.iterations=$result.iterations; $base.matvec_seconds=$result.matvec_seconds_avg
        $base.near_pairs=$result.fmm_near_pairs; $base.far_pairs=$result.fmm_far_pairs
        $base.residual_history=$result.residual_history
        $base.preconditioner_effective=Get-Z4gProfileValue $runDir 'preconditioner_effective'
        $base.failure_class=if ($result.status -eq 'ok') { 'pass' } else { 'not_converged' }
    } catch {
        $base.failure_class='profile_parse_error'
        $base.status='failed'
        Set-Content -LiteralPath (Join-Path $runDir 'parse_error.txt') -Value $_.Exception.Message
    }
    return [pscustomobject]$base
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($name in $wanted) {
    $fixture = @($manifest | Where-Object { $_.fixture -eq $name }) | Select-Object -First 1
    if ($null -eq $fixture) { throw "Missing required manifest row: $name" }
    for ($repeat = 1; $repeat -le $repeatCount; ++$repeat) {
        $rows.Add((Invoke-Z4gCase -Fixture $fixture -Backend 'fmm_ilu0' -Repeat $repeat)) | Out-Null
    }
    $rows.Add((Invoke-Z4gCase -Fixture $fixture -Backend 'grid_context' -Repeat 1)) | Out-Null
    $firstFmm = @($rows | Where-Object { $_.fixture -eq $name -and $_.backend -eq 'fmm_ilu0' })[0]
    if (-not $Quick -and $firstFmm.panels -gt 0 -and $firstFmm.panels -le 50000) {
        $rows.Add((Invoke-Z4gCase -Fixture $fixture -Backend 'tabi_direct' -Repeat 1)) | Out-Null
    }
    if (-not $Quick -and $name -eq 'barnase') {
        $rows.Add((Invoke-Z4gCase -Fixture $fixture -Backend 'block_jacobi' -Repeat 1 `
                -Preconditioner 'fmm-block-jacobi')) | Out-Null
    }
}

$rowArray = $rows.ToArray()
$rowArray | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
$repeatRows = New-Object System.Collections.Generic.List[object]
foreach ($name in $wanted) {
    $set = @($rowArray | Where-Object { $_.fixture -eq $name -and $_.backend -eq 'fmm_ilu0' })
    $energies = @($set | ForEach-Object { [double]$_.energy_kj_mol } | Where-Object { Test-FmmCliFinite $_ })
    $audits = @($set | ForEach-Object { [double]$_.audit_rms } | Where-Object { Test-FmmCliFinite $_ })
    $mean = if ($energies.Count) { ($energies | Measure-Object -Average).Average } else { [double]::NaN }
    $spread = if ($energies.Count) { ($energies | Measure-Object -Maximum).Maximum - ($energies | Measure-Object -Minimum).Minimum } else { [double]::NaN }
    $repeatRows.Add([pscustomobject]@{
            fixture=$name; repeats=$set.Count; energy_mean=$mean; energy_spread=$spread;
            audit_max=$(if ($audits.Count) { ($audits | Measure-Object -Maximum).Maximum } else { [double]::NaN });
            all_status_ok=(@($set | Where-Object { $_.status -ne 'ok' }).Count -eq 0)
        }) | Out-Null
}
$repeatRows.ToArray() | Export-Csv -LiteralPath $RepeatCsv -NoTypeInformation

$traceRows = New-Object System.Collections.Generic.List[object]
foreach ($row in @($rowArray | Where-Object { $_.fixture -eq 'barnase' -and $_.backend -in @('fmm_ilu0','block_jacobi') })) {
    foreach ($point in @(Convert-FmmCliResidualHistory $row.residual_history)) {
        $traceRows.Add([pscustomobject]@{
                fixture='barnase'; preconditioner=$(if ($row.backend -eq 'fmm_ilu0') { 'sparse_ilu0' } else { 'block_jacobi' });
                iteration=$point.iteration; residual=$point.residual
            }) | Out-Null
    }
}
$traceRows.ToArray() | Export-Csv -LiteralPath $TraceCsv -NoTypeInformation

function Write-Z4gRuntimeSvg([object[]]$Rows, [string]$Path) {
    $fmm = @($Rows | Where-Object { $_.backend -eq 'fmm_ilu0' -and $_.repeat -eq 1 })
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="420" viewBox="0 0 900 420">') | Out-Null
    $svg.Add('<rect width="900" height="420" fill="white"/>') | Out-Null
    $svg.Add('<text x="32" y="32" font-family="Arial" font-size="18">Full-protein runtime and memory</text>') | Out-Null
    $x=80
    foreach ($row in $fmm) {
        $height=[Math]::Min(280,[Math]::Max(3,[double]$row.wall_seconds))
        $y=340-$height
        $svg.Add(('<rect x="{0}" y="{1}" width="70" height="{2}" fill="#0072b2"/>' -f $x,$y,$height)) | Out-Null
        $svg.Add(('<text x="{0}" y="365" font-family="Arial" font-size="12">{1}</text>' -f ($x-4),$row.pdb_id)) | Out-Null
        $x+=190
    }
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

function Write-Z4gConvergenceSvg([object[]]$Trace, [string]$Path) {
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="520" viewBox="0 0 900 520">') | Out-Null
    $svg.Add('<rect width="900" height="520" fill="white"/>') | Out-Null
    $svg.Add('<text x="32" y="32" font-family="Arial" font-size="18">Barnase preconditioner convergence comparison</text>') | Out-Null
    foreach ($spec in @(@{name='block_jacobi';color='#d55e00'},@{name='sparse_ilu0';color='#0072b2'})) {
        $points = @($Trace | Where-Object { $_.preconditioner -eq $spec.name } | Sort-Object {[int]$_.iteration})
        $pathPoints = @()
        foreach ($p in $points) {
            $x=70+([double]$p.iteration/720.0)*760
            $log=[Math]::Log10([Math]::Max([double]$p.residual,1e-16))
            $y=70+((-1*$log)/16.0)*390
            $pathPoints += ("{0:F1},{1:F1}" -f $x,$y)
        }
        if ($pathPoints.Count -gt 0) {
            $svg.Add(('<polyline fill="none" stroke="{0}" stroke-width="3" points="{1}"/>' -f $spec.color,($pathPoints -join ' '))) | Out-Null
        }
    }
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

Write-Z4gRuntimeSvg -Rows $rowArray -Path $RuntimePlot
Write-Z4gConvergenceSvg -Trace $traceRows.ToArray() -Path $ConvergencePlot

$requiredFmm = @($rowArray | Where-Object { $_.backend -eq 'fmm_ilu0' })
$badStatus = @($requiredFmm | Where-Object { $_.status -ne 'ok' })
$badAccuracy = @($requiredFmm | Where-Object {
        -not (Test-FmmCliFinite ([double]$_.true_residual)) -or [double]$_.true_residual -gt 1e-6 -or
        -not (Test-FmmCliFinite ([double]$_.audit_rms)) -or [double]$_.audit_rms -gt 8e-3
    })
$badRepeat = @($repeatRows | Where-Object {
        -not $_.all_status_ok -or -not (Test-FmmCliFinite ([double]$_.energy_spread)) -or
        [double]$_.energy_spread -gt [Math]::Max(0.5, 1e-4*[Math]::Abs([double]$_.energy_mean))
    })
$decision = if ($badStatus.Count -gt 0) {
    if (@($badStatus | Where-Object { $_.status -eq 'timeout' }).Count) { 'z4g_full_system_runtime_limited' } else { 'z4g_operator_accuracy_gap' }
} elseif ($badAccuracy.Count -gt 0 -or $badRepeat.Count -gt 0) {
    'z4g_operator_accuracy_gap'
} else {
    'z4g_benchmark_ready'
}

[pscustomobject]@{
    decision=$decision; quick=[bool]$Quick; rows=$rowArray.Count; required_fmm_rows=$requiredFmm.Count;
    fixtures=($wanted -join ','); matrix=$MatrixCsv; repeatability=$RepeatCsv; trace=$TraceCsv
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation
@(
    'scope: Z4G public full-protein benchmark',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "fixtures: $($wanted -join ',')",
    "rows: $($rowArray.Count)",
    'sampling_policy: full_pqr_no_atom_sampling',
    'fmm_policy: theta=0.6 order=4 near=analytic_singular ilu_neighbor_cap=64'
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII
$evidenceLines = @($requiredFmm | ForEach-Object {
        '- {0} ({1}): atoms={2}, panels={3}, status={4}, residual={5}, audit_rms={6}, wall_s={7}, peak_rss_mb={8}' -f `
            $_.protein_name,$_.pdb_id,$_.input_atoms,$_.panels,$_.status,$_.true_residual,$_.audit_rms,$_.wall_seconds,$_.peak_rss_mb
    })
$reportLines = @(
    '# TITAN FMM/M2L Scope Z4G External Manuscript Report', '',
    "Decision: $decision", '',
    'This benchmark uses complete prepared PQR structures. Atom-subsampled development rows are excluded.', '',
    'Measured FMM rows:'
)
$reportLines += $evidenceLines
$reportLines += @(
    '',
    'Fixture preparation:',
    '- PDB2PQR 3.7.1, AMBER, pH 7.0, model 1 and chain A.',
    '- Non-positive PDB2PQR hydrogen radii are deterministically floored to 0.6 A and counted in the provenance manifest.',
    '',
    "Matrix: $MatrixCsv", "Fixture provenance: $manifestPath", "Runtime figure: $RuntimePlot",
    "Preconditioner convergence figure: $ConvergencePlot", '',
    'Grid remains the production Poisson-Boltzmann default.'
)
$reportLines | Set-Content -LiteralPath $RootReport -Encoding ASCII

if ($decision -ne 'z4g_benchmark_ready') {
    throw "Scope Z4G benchmark decision: $decision"
}
Write-Host "PASS Scope Z4G full-protein benchmark decision=$decision"
