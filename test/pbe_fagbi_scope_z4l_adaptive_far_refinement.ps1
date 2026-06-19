param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4l'),
    [int]$CaseTimeoutSeconds = 2400,
    [switch]$Quick,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$manifestPath = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4l_matrix.csv'
$refinementCsv = Join-Path $WorkRoot 'fagbi_scope_z4l_refinement.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4l_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4l_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4l_accuracy_work_plot.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4L_ADAPTIVE_FAR_REFINEMENT_REPORT.txt'

if (-not (Test-Path $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path $manifestPath)) { throw "Missing public fixture manifest: $manifestPath" }
$TitanExe = (Resolve-Path $TitanExe).Path
New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$manifest = @(Import-Csv $manifestPath)
$auditTargets = if ($Quick) { 128 } else { 256 }

function Get-Z4lFixture([string]$Name) {
    $rows = @($manifest | Where-Object fixture -eq $Name)
    if ($rows.Count -ne 1) { throw "Expected exactly one fixture row for $Name" }
    return $rows[0]
}

function Get-Z4lDouble([object]$Value, [double]$Default = [double]::NaN) {
    $parsed = 0.0
    if ($null -eq $Value -or -not [double]::TryParse([string]$Value,
            [Globalization.NumberStyles]::Float,
            [Globalization.CultureInfo]::InvariantCulture,[ref]$parsed)) {
        return $Default
    }
    return $parsed
}

function Invoke-Z4lCase([pscustomobject]$Fixture,[string]$RefineMode,
                        [double]$Tolerance) {
    $tolLabel = if ($RefineMode -eq 'off') { 'off' } else {
        $Tolerance.ToString('0.0000',[Globalization.CultureInfo]::InvariantCulture).Replace('.','p')
    }
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_$tolLabel"
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $arg = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
        -PqrPath (Join-Path $FixtureRoot $Fixture.pqr) -Backend tabi `
        -Mode treecode -FmmMode single-matvec-analytic -SurfacePoints 8 `
        -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount $auditTargets `
        -AuditMode sampled-detail -Preconditioner fmm-nearfield-ilu0 `
        -PreconditionerIluNeighborCap 64 -PreconditionerIluDropTol 0.0 `
        -NearIntegration analytic_singular -TaylorOrder 4 `
        -FmmCartesianAudit on -FmmLocal audit-p4 `
        -FmmLocalRefine $RefineMode -FmmLocalRefineTolerance $Tolerance
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $arg -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.timed_out) {
        return [pscustomobject]@{fixture=$Fixture.fixture;protein_name=$Fixture.protein_name
            pdb_id=$Fixture.pdb_id;refine_mode=$RefineMode;refine_tolerance=$Tolerance
            status='timeout';run_dir=$runDir;wall_seconds=$run.wall_seconds}
    }
    if ($run.exit_code -ne 0) {
        $log = (Get-Content -Raw (Join-Path $runDir 'stdout.log')) +
            (Get-Content -Raw (Join-Path $runDir 'stderr.log'))
        throw "Z4L case failed fixture=$($Fixture.fixture) mode=$RefineMode tol=$Tolerance exit=$($run.exit_code): $log"
    }
    $case = Read-FmmCliCaseResult -RunDir $runDir `
        -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    $summary = @(Import-Csv (Join-Path $runDir 'pbe_tabi_operator_attribution_summary.csv'))
    $aggregate = @($summary | Where-Object {
        $_.approximation_mode -eq 'local_m2l_p4' -and $_.probe -eq 'solution' -and
        $_.row_component -eq 'both' -and $_.region -eq 'total' -and
        $_.channel -eq 'aggregate'
    }) | Select-Object -First 1
    if ($null -eq $aggregate) { throw 'Missing Z4L local aggregate row' }
    foreach ($column in @('fmm_local_refine_mode','fmm_local_refine_tolerance',
            'fmm_local_tail_ratio_before','fmm_local_tail_ratio_after',
            'fmm_local_refined_cells','fmm_local_split_children',
            'fmm_local_direct_fallback_panels','fmm_local_direct_fallback_pairs',
            'fmm_local_far_pairs_before',
            'fmm_local_far_pairs_after','fmm_local_refinement_seconds')) {
        if (-not ($aggregate.PSObject.Properties.Name -contains $column)) {
            throw "Z4L summary is missing telemetry column $column"
        }
    }
    $beforeValue = Get-Z4lDouble $aggregate.fmm_local_far_pairs_before 0.0
    $before = [Math]::Max(1.0, $beforeValue)
    return [pscustomobject]@{
        fixture=$Fixture.fixture;protein_name=$Fixture.protein_name;pdb_id=$Fixture.pdb_id
        refine_mode=$RefineMode;refine_tolerance=$Tolerance;status=$case.status
        panels=$case.panels;audit_rms=Get-Z4lDouble $aggregate.aggregate_audit_rel
        recomposition_error=Get-Z4lDouble $aggregate.recomposition_rel_error
        reference_gap=Get-Z4lDouble $aggregate.reference_gap 0.0
        tail_ratio_before=Get-Z4lDouble $aggregate.fmm_local_tail_ratio_before 0.0
        tail_ratio_after=Get-Z4lDouble $aggregate.fmm_local_tail_ratio_after 0.0
        refined_cells=[int64]$aggregate.fmm_local_refined_cells
        split_children=[int64]$aggregate.fmm_local_split_children
        direct_fallback_panels=[int64]$aggregate.fmm_local_direct_fallback_panels
        direct_fallback_pairs=[int64]$aggregate.fmm_local_direct_fallback_pairs
        far_pairs_before=[int64]$aggregate.fmm_local_far_pairs_before
        far_pairs_after=[int64]$aggregate.fmm_local_far_pairs_after
        far_work_growth=(Get-Z4lDouble $aggregate.fmm_local_far_pairs_after 0.0)/$before
        refinement_seconds=Get-Z4lDouble $aggregate.fmm_local_refinement_seconds 0.0
        m2l_seconds=Get-Z4lDouble $aggregate.fmm_local_m2l_seconds 0.0
        l2p_seconds=Get-Z4lDouble $aggregate.fmm_local_l2p_seconds 0.0
        true_residual=$case.true_residual;wall_seconds=$run.wall_seconds
        peak_rss_mb=$run.peak_working_set_mb;run_dir=$runDir
    }
}

$matrix = New-Object System.Collections.Generic.List[object]
$oneA63 = Get-Z4lFixture '1a63'
$sweep = @(
    [pscustomobject]@{mode='off';tolerance=0.004},
    [pscustomobject]@{mode='p3-p4-tail';tolerance=0.008},
    [pscustomobject]@{mode='p3-p4-tail';tolerance=0.004},
    [pscustomobject]@{mode='p3-p4-tail';tolerance=0.002}
)
foreach ($entry in $sweep) {
    $matrix.Add((Invoke-Z4lCase $oneA63 $entry.mode $entry.tolerance)) | Out-Null
}

$rows = $matrix.ToArray()
$baseline = @($rows | Where-Object refine_mode -eq 'off') | Select-Object -First 1
$candidates = @($rows | Where-Object {
    $_.refine_mode -eq 'p3-p4-tail' -and $_.status -ne 'timeout'
})
$passing = @($candidates | Where-Object {
    [double]$_.audit_rms -le 8e-3 -and [double]$_.recomposition_error -le 1e-10 -and
    [double]$_.reference_gap -le 1e-10
} | Sort-Object @{Expression={[double]$_.refine_tolerance};Descending=$true},
    @{Expression={[double]$_.far_work_growth};Descending=$false})
$selected = if ($passing.Count) { $passing[0] } else {
    @($candidates | Sort-Object @{Expression={[double]$_.audit_rms};Descending=$false},
        @{Expression={[double]$_.far_work_growth};Descending=$false}) | Select-Object -First 1
}

if (-not $Quick -and $null -ne $selected) {
    foreach ($fixtureName in @('barnase','1sh1')) {
        $matrix.Add((Invoke-Z4lCase (Get-Z4lFixture $fixtureName) `
            'p3-p4-tail' ([double]$selected.refine_tolerance))) | Out-Null
    }
    $rows = $matrix.ToArray()
}

$baselineWall = [Math]::Max(1e-12, [double]$baseline.wall_seconds)
foreach ($row in $rows) {
    $runtimeGrowth = if ($row.fixture -eq '1a63') {
        [double]$row.wall_seconds / $baselineWall
    } else { [double]::NaN }
    $row | Add-Member -NotePropertyName runtime_growth_vs_baseline `
        -NotePropertyValue $runtimeGrowth -Force
}

$rows | Export-Csv $matrixCsv -NoTypeInformation
$rows | Select-Object fixture,refine_mode,refine_tolerance,audit_rms,
    tail_ratio_before,tail_ratio_after,refined_cells,split_children,
    direct_fallback_panels,direct_fallback_pairs,far_pairs_before,far_pairs_after,far_work_growth,
    refinement_seconds,m2l_seconds,l2p_seconds,runtime_growth_vs_baseline |
    Export-Csv $refinementCsv -NoTypeInformation

$decision = 'z4l_theta06_accuracy_pass'
if (@($rows | Where-Object status -eq 'timeout').Count) {
    $decision = 'z4l_runtime_limited'
} elseif (@($rows | Where-Object {
        $_.refine_mode -eq 'p3-p4-tail' -and
        ([double]$_.reference_gap -gt 1e-10 -or
         [double]$_.recomposition_error -gt 1e-10)
    }).Count) {
    $decision = 'z4l_reference_gap'
} elseif ($null -eq $selected -or [double]$selected.audit_rms -gt 8e-3) {
    $strictest = @($candidates | Sort-Object refine_tolerance | Select-Object -First 1)
    $improvement = if ($strictest.Count -and [double]$baseline.audit_rms -gt 0) {
        1.0 - [double]$strictest[0].audit_rms/[double]$baseline.audit_rms
    } else { 0.0 }
    $decision = if ($improvement -lt 0.2) {
        'z4l_estimator_uncorrelated'
    } else { 'z4l_accuracy_not_recovered' }
} elseif (-not $Quick) {
    $required = @($rows | Where-Object {
        $_.fixture -in @('barnase','1a63','1sh1') -and
        $_.refine_mode -eq 'p3-p4-tail' -and
        [double]$_.refine_tolerance -eq [double]$selected.refine_tolerance
    })
    if ($required.Count -ne 3 -or @($required | Where-Object {
            [double]$_.audit_rms -gt 8e-3
        }).Count) {
        $decision = 'z4l_accuracy_not_recovered'
    } elseif ([double]$selected.runtime_growth_vs_baseline -gt 3.0 -or
        @($required | Where-Object {
            [double]$_.far_work_growth -gt 3.0
        }).Count) {
        $decision = 'z4l_accuracy_pass_costly'
    }
} elseif ([double]$selected.far_work_growth -gt 3.0 -or
    [double]$selected.runtime_growth_vs_baseline -gt 3.0) {
    $decision = 'z4l_accuracy_pass_costly'
}

$selectedTolerance = if ($null -eq $selected) { [double]::NaN } else {
    [double]$selected.refine_tolerance
}
[pscustomobject]@{decision=$decision;quick=[bool]$Quick;rows=$rows.Count
    selected_tolerance=$selectedTolerance;audit_targets=$auditTargets} |
    Export-Csv $decisionCsv -NoTypeInformation
@('scope: Z4L adaptive far-cell refinement',"decision: $decision",
  "quick: $([bool]$Quick)","rows: $($rows.Count)",
  "selected_tolerance: $selectedTolerance",'theta: 0.6',
  'production_default: grid remains default','gmres_p4: disabled') |
    Set-Content $summaryPath -Encoding ASCII

$svg = @('<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="620" viewBox="0 0 1000 620">',
    '<rect width="1000" height="620" fill="white"/>',
    '<text x="35" y="38" font-family="Arial" font-size="20">Adaptive p3/p4 far-cell refinement at theta 0.6</text>')
$y = 82
foreach ($row in $rows) {
    $svg += ('<text x="35" y="{0}" font-family="Arial" font-size="13">{1} {2} tol={3}: audit={4}, tail={5}->{6}, far growth={7}x</text>' -f
        $y,$row.fixture,$row.refine_mode,$row.refine_tolerance,$row.audit_rms,
        $row.tail_ratio_before,$row.tail_ratio_after,$row.far_work_growth)
    $y += 38
}
$svg += '</svg>'
$svg | Set-Content $plotPath -Encoding ASCII

$report = @('# TITAN FMM/M2L Scope Z4L Adaptive Far-Cell Refinement Report','',
    "Decision: $decision",'',
    'Z4L uses a p3/p4 local-tail estimator to split only risky accepted far cells. The production degree-one FMM operator and GMRES remain unchanged.','')
foreach ($row in $rows) {
    $report += ('- {0} {1} tol={2}: audit={3}, tail={4}->{5}, refined={6}, direct_fallback_pairs={7}, far_pairs={8}->{9}, far_growth={10}x, runtime_growth={11}x' -f
        $row.fixture,$row.refine_mode,$row.refine_tolerance,$row.audit_rms,
        $row.tail_ratio_before,$row.tail_ratio_after,$row.refined_cells,
        $row.direct_fallback_pairs,$row.far_pairs_before,$row.far_pairs_after,
        $row.far_work_growth,$row.runtime_growth_vs_baseline)
}
$report += @('','Grid remains the production Poisson-Boltzmann default.',
    'GMRES-p4 remains fail-closed in Scope Z4L.')
$report | Set-Content $rootReport -Encoding ASCII

$bad = @(Get-ChildItem $WorkRoot -Recurse -Filter arg.in | ForEach-Object {
    $path = $_.FullName
    Get-Content $path | Where-Object {
        $_ -match '^\s*-' -and $_ -notmatch '^\s*-[^=\s]+='
    } | ForEach-Object { "$path :: $_" }
})
if ($bad.Count) { throw "Generated arg.in contains bare options: $($bad -join '; ')" }
if (-not $AllowFailedDecision -and $decision -ne 'z4l_theta06_accuracy_pass') {
    throw "Scope Z4L decision: $decision"
}
Write-Host "PASS Scope Z4L decision=$decision rows=$($rows.Count)"
