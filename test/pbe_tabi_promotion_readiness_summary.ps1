param(
    [string]$MediumCsv = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_topdown_physical_gap_matrix\pbe_tabi_topdown_physical_gap_matrix.csv'),
    [string]$BarnaseCsv = (Join-Path $PSScriptRoot 'tmp\b3s_300\pbe_tabi_barnase_truth_matrix.csv'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_promotion_readiness_summary'),
    [string]$MediumCandidateCaseId = 'stage3u_promotion_freeze_o4_theta0p60_surface_margin0p25_poly_abs1e7_cluster_off',
    [switch]$AllowFailedStatus
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Convert-ToDoubleOrNull([object]$Value) {
    if ($null -eq $Value) { return $null }
    if ($Value -is [double] -or $Value -is [single] -or $Value -is [decimal] -or
        $Value -is [int] -or $Value -is [long]) {
        return [double]$Value
    }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $null }
    $parsed = 0.0
    if ([double]::TryParse($text, [System.Globalization.NumberStyles]::Float,
            $InvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $null
}

function Format-Number([object]$Value) {
    if ($null -eq $Value) { return 'NA' }
    if ($Value -is [string]) { return $Value }
    return ([double]$Value).ToString('G17', $InvariantCulture)
}

function Get-RowValue([object]$Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}

function Get-FirstDouble($Row, [string[]]$Names) {
    foreach ($name in $Names) {
        $value = Convert-ToDoubleOrNull (Get-RowValue $Row $name)
        if ($null -ne $value) { return $value }
    }
    return $null
}

function Get-Product($Left, $Right) {
    $leftValue = Convert-ToDoubleOrNull $Left
    $rightValue = Convert-ToDoubleOrNull $Right
    if ($null -eq $leftValue -or $null -eq $rightValue) { return $null }
    return $leftValue * $rightValue
}

function Get-Ratio($Numerator, $Denominator) {
    $numValue = Convert-ToDoubleOrNull $Numerator
    $denValue = Convert-ToDoubleOrNull $Denominator
    if ($null -eq $numValue -or $null -eq $denValue -or [Math]::Abs($denValue) -le 1.0e-300) {
        return $null
    }
    return $numValue / $denValue
}

function Get-ChannelMissTotal($Row) {
    $total = 0.0
    $seen = $false
    foreach ($name in @('convergence_k1_high_error_missed',
                        'convergence_k2_high_error_missed',
                        'convergence_k3_high_error_missed',
                        'convergence_k4_high_error_missed')) {
        $value = Convert-ToDoubleOrNull (Get-RowValue $Row $name)
        if ($null -ne $value) {
            $seen = $true
            $total += $value
        }
    }
    if ($seen) { return $total }
    return Convert-ToDoubleOrNull (Get-RowValue $Row 'convergence_estimator_high_error_missed')
}

function Get-Reason([hashtable]$Gates) {
    $failed = New-Object System.Collections.Generic.List[string]
    foreach ($key in @($Gates.Keys | Sort-Object)) {
        if (-not $Gates[$key]) { $failed.Add($key) }
    }
    if ($failed.Count -eq 0) { return 'pass' }
    return ($failed -join ';')
}

function Get-BarnaseTimingRow($Row) {
    $runDir = Get-RowValue $Row 'run_dir'
    if ($runDir -eq 'NA') { return $null }
    $timingPath = Join-Path $runDir 'pbe_timing.csv'
    if (-not (Test-Path -LiteralPath $timingPath)) { return $null }
    return @(Import-Csv -LiteralPath $timingPath | Select-Object -First 1)[0]
}

function Get-NormalizedBarnaseMetric($Row, $Timing, [string]$SummaryName, [string]$TimingName) {
    $value = Get-RowValue $Row $SummaryName
    if ($value -ne 'NA') { return $value }
    if ($null -ne $Timing) {
        return Get-RowValue $Timing $TimingName
    }
    return 'NA'
}

function New-MediumSummaryRow($Row, [double]$GroupedResidual) {
    $matvecRel = Convert-ToDoubleOrNull (Get-RowValue $Row 'tree_vs_direct_matvec_rel')
    $energyDrift = Convert-ToDoubleOrNull (Get-RowValue $Row 'energy_drift_rel')
    $sampledResidual = Convert-ToDoubleOrNull (Get-RowValue $Row 'sampled_true_residual_rel')
    $farMean = Convert-ToDoubleOrNull (Get-RowValue $Row 'far_expansion_rel_error_mean')
    $missedTotal = Get-ChannelMissTotal $Row
    $statusPass = (Get-RowValue $Row 'status') -eq 'ok'
    $matvecPass = $null -ne $matvecRel -and $matvecRel -lt 2.5e-3
    $energyPass = $null -ne $energyDrift -and $energyDrift -lt 5.0e-3
    $residualPass = $null -ne $sampledResidual -and $GroupedResidual -gt 0.0 -and
        $sampledResidual -le (2.0 * [Math]::Max(1.0e-16, $GroupedResidual))
    $farMeanPass = $null -ne $farMean -and $farMean -lt 5.0e-5
    $missedPass = $null -ne $missedTotal -and $missedTotal -eq 0.0
    $gates = @{
        status_pass = $statusPass
        matvec_pass = $matvecPass
        energy_pass = $energyPass
        residual_pass = $residualPass
        far_mean_pass = $farMeanPass
        accepted_missed_pass = $missedPass
    }
    $overall = -not ($gates.Values -contains $false)

    return [PSCustomObject]@{
        scope = 'medium'
        fixture = Get-RowValue $Row 'fixture'
        candidate_role = 'promotion_candidate'
        case_id = Get-RowValue $Row 'case_id'
        production_default_changed = 'false'
        status = Get-RowValue $Row 'status'
        performance_classification = Get-RowValue $Row 'performance_classification'
        cluster_policy = Get-RowValue $Row 'cluster_policy'
        accepted_miss_absolute_floor = Get-RowValue $Row 'accepted_miss_absolute_floor'
        matvec_rel = Format-Number $matvecRel
        energy_drift_rel = Format-Number $energyDrift
        residual_rel = Format-Number $sampledResidual
        grouped_residual_rel = Format-Number $GroupedResidual
        far_expansion_rel_error_mean = Format-Number $farMean
        accepted_missed_total = Format-Number $missedTotal
        phantom_k3_abs_bound = Get-RowValue $Row 'phantom_k3_abs_bound'
        phantom_k3_abs_over_phi_norm = Get-RowValue $Row 'phantom_k3_abs_over_phi_norm'
        phantom_k3_abs_over_phin_norm = Get-RowValue $Row 'phantom_k3_abs_over_phin_norm'
        status_pass = [string]$statusPass
        matvec_pass = [string]$matvecPass
        energy_pass = [string]$energyPass
        residual_pass = [string]$residualPass
        far_mean_pass = [string]$farMeanPass
        accepted_missed_pass = [string]$missedPass
        phantom_ratio_pass = 'NA'
        phantom_bound_pass = 'NA'
        overall_pass = [string]$overall
        reason = Get-Reason $gates
        source_csv = $MediumCsv
    }
}

function New-MediumReportRow($Row, [string]$Role) {
    return [PSCustomObject]@{
        scope = 'medium_report'
        fixture = Get-RowValue $Row 'fixture'
        candidate_role = $Role
        case_id = Get-RowValue $Row 'case_id'
        production_default_changed = 'false'
        status = Get-RowValue $Row 'status'
        performance_classification = Get-RowValue $Row 'performance_classification'
        cluster_policy = Get-RowValue $Row 'cluster_policy'
        accepted_miss_absolute_floor = Get-RowValue $Row 'accepted_miss_absolute_floor'
        matvec_rel = Get-RowValue $Row 'tree_vs_direct_matvec_rel'
        energy_drift_rel = Get-RowValue $Row 'energy_drift_rel'
        residual_rel = Get-RowValue $Row 'sampled_true_residual_rel'
        grouped_residual_rel = 'NA'
        far_expansion_rel_error_mean = Get-RowValue $Row 'far_expansion_rel_error_mean'
        accepted_missed_total = Format-Number (Get-ChannelMissTotal $Row)
        phantom_k3_abs_bound = Get-RowValue $Row 'phantom_k3_abs_bound'
        phantom_k3_abs_over_phi_norm = Get-RowValue $Row 'phantom_k3_abs_over_phi_norm'
        phantom_k3_abs_over_phin_norm = Get-RowValue $Row 'phantom_k3_abs_over_phin_norm'
        status_pass = 'report'
        matvec_pass = 'report'
        energy_pass = 'report'
        residual_pass = 'report'
        far_mean_pass = 'report'
        accepted_missed_pass = 'report'
        phantom_ratio_pass = 'report'
        phantom_bound_pass = 'report'
        overall_pass = 'report'
        reason = 'report_only'
        source_csv = $MediumCsv
    }
}

function New-BarnaseSummaryRow($Row) {
    $timing = Get-BarnaseTimingRow $Row
    $phiNorm = Get-NormalizedBarnaseMetric $Row $timing 'phi_norm' 'tabi_phi_norm'
    $phinNorm = Get-NormalizedBarnaseMetric $Row $timing 'phin_norm' 'tabi_phin_norm'
    $k3Count = Get-RowValue $Row 'phantom_high_error_k3_count'
    $k3AbsMax = Get-RowValue $Row 'phantom_high_error_k3_abs_max'
    $k3Bound = Convert-ToDoubleOrNull (Get-RowValue $Row 'phantom_k3_abs_bound')
    if ($null -eq $k3Bound) { $k3Bound = Get-Product $k3Count $k3AbsMax }
    $k3OverPhi = Convert-ToDoubleOrNull (Get-RowValue $Row 'phantom_k3_abs_over_phi_norm')
    if ($null -eq $k3OverPhi) { $k3OverPhi = Get-Ratio $k3AbsMax $phiNorm }
    $k3OverPhin = Convert-ToDoubleOrNull (Get-RowValue $Row 'phantom_k3_abs_over_phin_norm')
    if ($null -eq $k3OverPhin) { $k3OverPhin = Get-Ratio $k3AbsMax $phinNorm }

    $residual = Get-FirstDouble $Row @('timing_final_residual',
                                       'final_residual',
                                       'timing_sampled_true_residual_rel',
                                       'sampled_true_residual_rel')
    $farMean = Convert-ToDoubleOrNull (Get-RowValue $Row 'far_expansion_rel_error_mean')
    $missedTotal = Get-ChannelMissTotal $Row
    $statusPass = (Get-RowValue $Row 'status') -eq 'ok'
    $residualPass = $null -ne $residual -and $residual -le 1.0e-4
    $farMeanPass = $null -ne $farMean -and $farMean -lt 5.0e-5
    $missedPass = $null -ne $missedTotal -and $missedTotal -eq 0.0
    $phantomRatioPass = $null -ne $k3OverPhi -and $k3OverPhi -lt 1.0e-4
    $phantomBoundPass = $null -ne $k3Bound -and $k3Bound -lt 1.0e-2
    $gates = @{
        status_pass = $statusPass
        residual_pass = $residualPass
        far_mean_pass = $farMeanPass
        accepted_missed_pass = $missedPass
        phantom_ratio_pass = $phantomRatioPass
        phantom_bound_pass = $phantomBoundPass
    }
    $overall = -not ($gates.Values -contains $false)

    return [PSCustomObject]@{
        scope = 'barnase_300s'
        fixture = 'barnase'
        candidate_role = 'promotion_candidate'
        case_id = Get-RowValue $Row 'case_id'
        production_default_changed = 'false'
        status = Get-RowValue $Row 'status'
        performance_classification = Get-RowValue $Row 'decision_classification'
        cluster_policy = Get-RowValue $Row 'cluster_policy'
        accepted_miss_absolute_floor = Get-RowValue $Row 'accepted_miss_absolute_floor'
        matvec_rel = 'NA'
        energy_drift_rel = 'NA'
        residual_rel = Format-Number $residual
        grouped_residual_rel = 'NA'
        far_expansion_rel_error_mean = Format-Number $farMean
        accepted_missed_total = Format-Number $missedTotal
        phantom_k3_abs_bound = Format-Number $k3Bound
        phantom_k3_abs_over_phi_norm = Format-Number $k3OverPhi
        phantom_k3_abs_over_phin_norm = Format-Number $k3OverPhin
        status_pass = [string]$statusPass
        matvec_pass = 'NA'
        energy_pass = 'NA'
        residual_pass = [string]$residualPass
        far_mean_pass = [string]$farMeanPass
        accepted_missed_pass = [string]$missedPass
        phantom_ratio_pass = [string]$phantomRatioPass
        phantom_bound_pass = [string]$phantomBoundPass
        overall_pass = [string]$overall
        reason = Get-Reason $gates
        source_csv = $BarnaseCsv
    }
}

if (-not (Test-Path -LiteralPath $MediumCsv)) {
    throw "Missing medium matrix CSV: $MediumCsv"
}
if (-not (Test-Path -LiteralPath $BarnaseCsv)) {
    throw "Missing Barnase truth CSV: $BarnaseCsv"
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$mediumRows = @(Import-Csv -LiteralPath $MediumCsv)
$barnaseRows = @(Import-Csv -LiteralPath $BarnaseCsv)
$summary = New-Object System.Collections.Generic.List[object]

$candidateRows = @($mediumRows | Where-Object { $_.case_id -eq $MediumCandidateCaseId })
if ($candidateRows.Count -eq 0) {
    throw "No medium candidate rows found for case_id=$MediumCandidateCaseId in $MediumCsv"
}
foreach ($fixture in @($candidateRows | Select-Object -ExpandProperty fixture -Unique)) {
    $row = @($candidateRows | Where-Object { $_.fixture -eq $fixture } | Select-Object -First 1)[0]
    $grouped = @($mediumRows | Where-Object {
        $_.fixture -eq $fixture -and $_.case_id -eq 'grouped_morton_off'
    } | Select-Object -First 1)
    $groupedResidual = if ($grouped.Count -gt 0) {
        Convert-ToDoubleOrNull (Get-RowValue $grouped[0] 'sampled_true_residual_rel')
    } else {
        $null
    }
    $summary.Add((New-MediumSummaryRow $row ([double]$groupedResidual)))
}

foreach ($row in @($mediumRows | Where-Object {
        $_.case_id -like 'stage3u_timing_*' -or $_.case_id -like 'stage3u_floor_sensitivity_*'
    })) {
    $role = if ($row.case_id -like 'stage3u_timing_*') { 'timing_report' } else { 'floor_sensitivity' }
    $summary.Add((New-MediumReportRow $row $role))
}

$barnaseCandidate = @($barnaseRows | Where-Object {
    $_.backend -eq 'tabi' -and $_.sum_mode -eq 'treecode' -and
        $_.status -eq 'ok' -and
        ((Get-RowValue $_ 'promotion_readiness') -eq 'pass' -or
         ((Get-RowValue $_ 'requested_tree_builder') -eq 'topdown-aabb' -and
          (Get-RowValue $_ 'mac_distance_policy') -eq 'surface-margin' -and
          (Get-RowValue $_ 'convergence_channel_calibration') -eq 'polynomial-eps'))
} | Select-Object -First 1)
if ($barnaseCandidate.Count -eq 0) {
    $barnaseCandidate = @($barnaseRows | Where-Object {
        $_.backend -eq 'tabi' -and $_.sum_mode -eq 'treecode' -and $_.status -eq 'ok'
    } | Select-Object -First 1)
}
if ($barnaseCandidate.Count -eq 0) {
    throw "No completed Barnase treecode candidate row found in $BarnaseCsv"
}
$summary.Add((New-BarnaseSummaryRow $barnaseCandidate[0]))

$mediumPass = @($summary | Where-Object { $_.scope -eq 'medium' -and $_.overall_pass -ne 'True' }).Count -eq 0
$barnasePass = @($summary | Where-Object { $_.scope -eq 'barnase_300s' -and $_.overall_pass -eq 'True' }).Count -gt 0
$overallPass = $mediumPass -and $barnasePass
$summary.Add([PSCustomObject]@{
    scope = 'overall'
    fixture = 'all'
    candidate_role = 'promotion_readiness_package'
    case_id = $MediumCandidateCaseId
    production_default_changed = 'false'
    status = if ($overallPass) { 'ok' } else { 'failed' }
    performance_classification = if ($overallPass) { 'promotion_ready_hidden_candidate' } else { 'promotion_readiness_failed' }
    cluster_policy = 'off'
    accepted_miss_absolute_floor = '1e-7'
    matvec_rel = 'NA'
    energy_drift_rel = 'NA'
    residual_rel = 'NA'
    grouped_residual_rel = 'NA'
    far_expansion_rel_error_mean = 'NA'
    accepted_missed_total = 'NA'
    phantom_k3_abs_bound = 'NA'
    phantom_k3_abs_over_phi_norm = 'NA'
    phantom_k3_abs_over_phin_norm = 'NA'
    status_pass = [string]$overallPass
    matvec_pass = 'NA'
    energy_pass = 'NA'
    residual_pass = 'NA'
    far_mean_pass = 'NA'
    accepted_missed_pass = 'NA'
    phantom_ratio_pass = 'NA'
    phantom_bound_pass = 'NA'
    overall_pass = [string]$overallPass
    reason = if ($overallPass) { 'pass' } else { 'medium_or_barnase_gate_failed' }
    source_csv = "$MediumCsv;$BarnaseCsv"
})

$outCsv = Join-Path $WorkRoot 'pbe_tabi_promotion_readiness_summary.csv'
$summary | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII

if (-not $overallPass -and -not $AllowFailedStatus) {
    throw "PBE TABI promotion readiness summary failed; see $outCsv"
}

Write-Host "PBE TABI promotion readiness summary COMPLETE: $outCsv overall_pass=$overallPass production_default_changed=false"
