param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4h'),
    [int]$CaseTimeoutSeconds = 2400,
    [int]$Repeats = 3,
    [switch]$Quick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$ManifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$SweepCsv = Join-Path $WorkRoot 'fagbi_scope_z4h_sweep.csv'
$ConfirmationCsv = Join-Path $WorkRoot 'fagbi_scope_z4h_confirmation.csv'
$RepeatabilityCsv = Join-Path $WorkRoot 'fagbi_scope_z4h_repeatability.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4h_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4h_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4h_accuracy_runtime.svg'
$RootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4H_FULL_PROTEIN_ACCURACY_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $ManifestPath)) {
    throw "Missing public fixture manifest. Run test\pbe_fagbi_scope_z4g_prepare_public_fixtures.ps1 -RefreshFixtures"
}
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$manifest = @(Import-Csv -LiteralPath $ManifestPath)

function Get-Z4hFixture([string]$Name) {
    $rows = @($manifest | Where-Object { $_.fixture -eq $Name })
    if ($rows.Count -ne 1) { throw "Expected exactly one fixture manifest row for $Name" }
    return $rows[0]
}

function Get-Z4hAtomCount([string]$Path) {
    return @(Get-Content -LiteralPath $Path | Where-Object {
            $_.StartsWith('ATOM') -or $_.StartsWith('HETATM')
        }).Count
}

function Get-Z4hProfileValue([string]$RunDir, [string]$Name, [string]$Default = '') {
    $path = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $path)) { return $Default }
    $rows = @(Import-Csv -LiteralPath $path)
    if ($rows.Count -lt 1 -or $rows[0].PSObject.Properties.Name -notcontains $Name) { return $Default }
    return [string]$rows[0].$Name
}

function Get-Z4hMedian([double[]]$Values) {
    $finite = @($Values | Where-Object { Test-FmmCliFinite $_ } | Sort-Object)
    if ($finite.Count -eq 0) { return [double]::NaN }
    $middle = [int][Math]::Floor($finite.Count / 2)
    if (($finite.Count % 2) -eq 1) { return [double]$finite[$middle] }
    return 0.5 * ([double]$finite[$middle - 1] + [double]$finite[$middle])
}

function Convert-Z4hThetaLabel([double]$Theta) {
    return $Theta.ToString('0.00', [Globalization.CultureInfo]::InvariantCulture).Replace('.', 'p')
}

function Invoke-Z4hFmmCase(
    [pscustomobject]$Fixture,
    [string]$Stage,
    [double]$Theta,
    [int]$Repeat,
    [int]$AuditSampleCount
) {
    $pqr = Join-Path $FixtureRoot $Fixture.pqr
    $inputAtoms = Get-Z4hAtomCount $pqr
    if ($inputAtoms -ne [int]$Fixture.atom_count) { throw "Atom-count drift for $($Fixture.fixture)" }
    $condition = '{0}_{1}_theta{2}_r{3}' -f $Stage,$Fixture.fixture,(Convert-Z4hThetaLabel $Theta),$Repeat
    $runDir = Join-Path $WorkRoot $condition
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    $argPath = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot -PqrPath $pqr `
        -Backend tabi -Mode treecode -FmmMode single-matvec-analytic `
        -SurfacePoints 8 -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 `
        -SaltM 0.15 -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta $Theta -AuditSampleCount $AuditSampleCount `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerBlockPanelCap 32 `
        -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular -TaylorOrder 4
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $row = [ordered]@{
        stage=$Stage; fixture=$Fixture.fixture; protein_name=$Fixture.protein_name; pdb_id=$Fixture.pdb_id;
        theta=$Theta; repeat=$Repeat; audit_sample_count=$AuditSampleCount;
        input_atoms=$inputAtoms; analyzed_atoms=$inputAtoms; panels=0; status='failed'; failure_class='exit_code';
        energy_kj_mol=[double]::NaN; true_residual=[double]::NaN; audit_rms=[double]::NaN;
        iterations=0; matvec_seconds=[double]::NaN; wall_seconds=$run.wall_seconds;
        peak_rss_mb=$run.peak_working_set_mb; near_pairs=0; far_pairs=0;
        preconditioner_effective=''; fmm_mode='single-matvec-analytic'; expansion_order=4;
        residual_history=''; arg_file=$argPath
    }
    if ($run.timed_out) {
        $row.status='timeout'; $row.failure_class='timeout_not_verified'
        return [pscustomobject]$row
    }
    if ($run.exit_code -ne 0) {
        $row.failure_class="exit_code_$($run.exit_code)"
        return [pscustomobject]$row
    }
    try {
        $result = Read-FmmCliCaseResult -RunDir $runDir `
            -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
        $row.panels=$result.panels; $row.status=$result.status
        $row.energy_kj_mol=$result.energy_kj_mol; $row.true_residual=$result.true_residual
        $row.audit_rms=$result.tree_vs_direct_matvec_rel; $row.iterations=$result.iterations
        $row.matvec_seconds=$result.matvec_seconds_avg; $row.near_pairs=$result.fmm_near_pairs
        $row.far_pairs=$result.fmm_far_pairs; $row.residual_history=$result.residual_history
        $row.preconditioner_effective=Get-Z4hProfileValue $runDir 'preconditioner_effective'
        $row.expansion_order=[int](Get-Z4hProfileValue $runDir 'fmm_m2l_expansion_order' '0')
        $row.failure_class=if ($result.status -eq 'ok') { 'pass' } else { 'not_converged' }
    } catch {
        $row.status='failed'; $row.failure_class='profile_parse_error'
        Set-Content -LiteralPath (Join-Path $runDir 'parse_error.txt') -Value $_.Exception.Message
    }
    return [pscustomobject]$row
}

function Test-Z4hAccuracyRow([pscustomobject]$Row) {
    return $Row.status -eq 'ok' -and
        (Test-FmmCliFinite ([double]$Row.true_residual)) -and [double]$Row.true_residual -le 1e-6 -and
        (Test-FmmCliFinite ([double]$Row.audit_rms)) -and [double]$Row.audit_rms -le 8e-3 -and
        [int]$Row.input_atoms -eq [int]$Row.analyzed_atoms -and
        $Row.preconditioner_effective -eq 'fmm_nearfield_ilu0' -and [int]$Row.expansion_order -eq 4
}

$sweepThetas = if ($Quick) { @(0.6,0.45,0.4) } else { @(0.6,0.5,0.45,0.4,0.35) }
$sweepRows = New-Object System.Collections.Generic.List[object]
$oneSh1 = Get-Z4hFixture '1sh1'
foreach ($theta in $sweepThetas) {
    $sweepRows.Add((Invoke-Z4hFmmCase -Fixture $oneSh1 -Stage 'sweep' -Theta $theta `
            -Repeat 1 -AuditSampleCount 512)) | Out-Null
}
$sweepArray = $sweepRows.ToArray()
$sweepArray | Export-Csv -LiteralPath $SweepCsv -NoTypeInformation
$candidateRows = @($sweepArray | Where-Object { Test-Z4hAccuracyRow $_ } | Sort-Object {[double]$_.theta} -Descending)
$candidateTheta = if ($candidateRows.Count -gt 0) { [double]$candidateRows[0].theta } else { [double]::NaN }

$confirmationRows = New-Object System.Collections.Generic.List[object]
if (-not $Quick -and (Test-FmmCliFinite $candidateTheta)) {
    foreach ($fixtureName in @('barnase','1a63','1sh1')) {
        $fixture = Get-Z4hFixture $fixtureName
        $confirmationRows.Add((Invoke-Z4hFmmCase -Fixture $fixture -Stage 'baseline' -Theta 0.6 `
                -Repeat 1 -AuditSampleCount 2048)) | Out-Null
        for ($repeat=1; $repeat -le [Math]::Max(1,$Repeats); ++$repeat) {
            $confirmationRows.Add((Invoke-Z4hFmmCase -Fixture $fixture -Stage 'candidate' `
                    -Theta $candidateTheta -Repeat $repeat -AuditSampleCount 2048)) | Out-Null
        }
    }
}
$confirmationArray = $confirmationRows.ToArray()
$confirmationArray | Export-Csv -LiteralPath $ConfirmationCsv -NoTypeInformation

$repeatabilityRows = New-Object System.Collections.Generic.List[object]
if ($confirmationArray.Count -gt 0) {
    foreach ($fixtureName in @('barnase','1a63','1sh1')) {
        $candidateSet = @($confirmationArray | Where-Object { $_.fixture -eq $fixtureName -and $_.stage -eq 'candidate' })
        $baseline = @($confirmationArray | Where-Object { $_.fixture -eq $fixtureName -and $_.stage -eq 'baseline' }) | Select-Object -First 1
        $energies = @($candidateSet | ForEach-Object { [double]$_.energy_kj_mol } | Where-Object { Test-FmmCliFinite $_ })
        $energyMean = if ($energies.Count) { ($energies | Measure-Object -Average).Average } else { [double]::NaN }
        $energySpread = if ($energies.Count) { ($energies | Measure-Object -Maximum).Maximum - ($energies | Measure-Object -Minimum).Minimum } else { [double]::NaN }
        $candidateWall = Get-Z4hMedian ([double[]]@($candidateSet | ForEach-Object { [double]$_.wall_seconds }))
        $baselineWall = if ($null -ne $baseline) { [double]$baseline.wall_seconds } else { [double]::NaN }
        $slowdown = if ((Test-FmmCliFinite $candidateWall) -and (Test-FmmCliFinite $baselineWall) -and $baselineWall -gt 0) { $candidateWall/$baselineWall } else { [double]::NaN }
        $auditValues = @($candidateSet | ForEach-Object { [double]$_.audit_rms } | Where-Object { Test-FmmCliFinite $_ })
        $repeatabilityRows.Add([pscustomobject]@{
                fixture=$fixtureName; repeats=$candidateSet.Count;
                all_accuracy_pass=(@($candidateSet | Where-Object { -not (Test-Z4hAccuracyRow $_) }).Count -eq 0);
                energy_mean=$energyMean; energy_spread=$energySpread;
                energy_spread_limit=$(if (Test-FmmCliFinite $energyMean) { [Math]::Max(0.5,1e-4*[Math]::Abs($energyMean)) } else { [double]::NaN });
                audit_max=$(if ($auditValues.Count) { ($auditValues | Measure-Object -Maximum).Maximum } else { [double]::NaN });
                candidate_wall_median=$candidateWall; baseline_wall=$baselineWall; runtime_slowdown=$slowdown
            }) | Out-Null
    }
}
$repeatArray = $repeatabilityRows.ToArray()
$repeatArray | Export-Csv -LiteralPath $RepeatabilityCsv -NoTypeInformation

$decision = if (-not (Test-FmmCliFinite $candidateTheta)) {
    if (@($sweepArray | Where-Object { $_.status -eq 'timeout' }).Count -gt 0) { 'z4h_runtime_limited' } else { 'z4h_no_safe_theta' }
} elseif ($Quick) {
    'z4h_quick_candidate_found'
} elseif (@($confirmationArray | Where-Object { $_.status -eq 'timeout' }).Count -gt 0) {
    'z4h_runtime_limited'
} else {
    $failedConfirmation = @($confirmationArray | Where-Object { $_.stage -eq 'candidate' -and -not (Test-Z4hAccuracyRow $_) })
    $repeatFailures = @($repeatArray | Where-Object {
            -not $_.all_accuracy_pass -or -not (Test-FmmCliFinite ([double]$_.energy_spread)) -or
            [double]$_.energy_spread -gt [double]$_.energy_spread_limit
        })
    if ($failedConfirmation.Count -gt 0 -or $repeatFailures.Count -gt 0) {
        'z4h_fixture_specific_policy'
    } else {
        $slowdowns = [double[]]@($repeatArray | ForEach-Object { [double]$_.runtime_slowdown })
        $medianSlowdown = Get-Z4hMedian $slowdowns
        if ((Test-FmmCliFinite $medianSlowdown) -and $medianSlowdown -le 3.0) {
            'z4h_theta_recovered'
        } else {
            'z4h_accuracy_recovered_runtime_costly'
        }
    }
}

function Write-Z4hSvg([object[]]$Rows, [string]$Path) {
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="960" height="520" viewBox="0 0 960 520">') | Out-Null
    $svg.Add('<rect width="960" height="520" fill="white"/>') | Out-Null
    $svg.Add('<text x="32" y="34" font-family="Arial" font-size="20">Full 1SH1 opening-theta accuracy and runtime</text>') | Out-Null
    $svg.Add('<line x1="70" y1="420" x2="900" y2="420" stroke="#333"/>') | Out-Null
    $svg.Add('<line x1="70" y1="70" x2="70" y2="420" stroke="#333"/>') | Out-Null
    $ordered = @($Rows | Sort-Object {[double]$_.theta})
    $points = @()
    foreach ($row in $ordered) {
        $x = 100 + (([double]$row.theta - 0.35)/0.25)*740
        $audit = [Math]::Min(0.025,[Math]::Max(0,[double]$row.audit_rms))
        $y = 420 - ($audit/0.025)*330
        $points += ('{0:F1},{1:F1}' -f $x,$y)
        $svg.Add(('<circle cx="{0:F1}" cy="{1:F1}" r="6" fill="#0072b2"/>' -f $x,$y)) | Out-Null
        $svg.Add(('<text x="{0:F1}" y="450" font-family="Arial" font-size="12">{1}</text>' -f ($x-15),$row.theta)) | Out-Null
        $svg.Add(('<text x="{0:F1}" y="{1:F1}" font-family="Arial" font-size="11">{2:F3}s</text>' -f ($x-20),($y-12),[double]$row.matvec_seconds)) | Out-Null
    }
    if ($points.Count -gt 0) { $svg.Add(('<polyline fill="none" stroke="#0072b2" stroke-width="3" points="{0}"/>' -f ($points -join ' '))) | Out-Null }
    $gateY = 420 - (0.008/0.025)*330
    $svg.Add(('<line x1="70" y1="{0:F1}" x2="900" y2="{0:F1}" stroke="#d55e00" stroke-dasharray="8 6"/>' -f $gateY)) | Out-Null
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}
Write-Z4hSvg -Rows $sweepArray -Path $PlotPath

[pscustomobject]@{
    decision=$decision; quick=[bool]$Quick; candidate_theta=$candidateTheta;
    sweep_rows=$sweepArray.Count; confirmation_rows=$confirmationArray.Count;
    sweep=$SweepCsv; confirmation=$ConfirmationCsv; repeatability=$RepeatabilityCsv
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation
@(
    'scope: Z4H full-protein FMM accuracy recovery gate',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "candidate_theta: $candidateTheta",
    "sweep_rows: $($sweepArray.Count)",
    "confirmation_rows: $($confirmationArray.Count)",
    'audit_gate: 8e-3',
    'confirmation_audit_samples: 2048',
    'default_policy_note: no production or TABI treecode default is changed by Z4H'
) | Set-Content -LiteralPath $SummaryPath -Encoding ASCII

$evidence = @($sweepArray | ForEach-Object {
        '- theta={0}: status={1}, residual={2}, audit={3}, matvec_s={4}, wall_s={5}' -f `
            $_.theta,$_.status,$_.true_residual,$_.audit_rms,$_.matvec_seconds,$_.wall_seconds
    })
$confirmationEvidence = @($repeatArray | ForEach-Object {
        '- {0}: candidate_audit_max={1}, energy_spread={2}, median_wall_s={3}, slowdown_vs_theta06={4}, all_accuracy_pass={5}' -f `
            $_.fixture,$_.audit_max,$_.energy_spread,$_.candidate_wall_median,$_.runtime_slowdown,$_.all_accuracy_pass
    })
$report = @('# TITAN FMM/M2L Scope Z4H Full-Protein Accuracy Report','',"Decision: $decision",'',
    'Z4H tightens the FMM opening criterion on complete prepared protein fixtures without changing the production or TABI defaults.','',
    'Full 1SH1 sweep:')
$report += $evidence
$report += @('',"Selected candidate theta: $candidateTheta",'','Full-protein confirmation:')
$report += $confirmationEvidence
$report += @('','Interpretation:',
    '- A theta that passes 1SH1 is not a transferable full-protein policy when Barnase or 1A63 fails the operator-audit gate.',
    '- Converged GMRES residuals do not override failed direct-operator parity; Z4G remains fail-closed.',
    '',"Sweep: $SweepCsv","Confirmation: $ConfirmationCsv","Repeatability: $RepeatabilityCsv","Plot: $PlotPath",'',
    'Grid remains the production Poisson-Boltzmann default. Z4G remains fail-closed until a full Z4H confirmation passes.')
$report | Set-Content -LiteralPath $RootReport -Encoding ASCII

$passDecisions = @('z4h_quick_candidate_found','z4h_theta_recovered','z4h_accuracy_recovered_runtime_costly')
if (-not $AllowFailedDecision -and $decision -notin $passDecisions) {
    throw "Scope Z4H decision: $decision"
}
Write-Host "PASS Scope Z4H decision=$decision candidate_theta=$candidateTheta"
