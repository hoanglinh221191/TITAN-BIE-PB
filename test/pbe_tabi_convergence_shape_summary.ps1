param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_topdown_physical_gap_matrix'),
    [string]$OutCsv = '',
    [string[]]$Fixtures = @('mini-protein', 'charged-peptide', 'compact-mini-protein'),
    [string]$CaseId = 'stage3ah_shape_audit_detail'
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

if ([string]::IsNullOrWhiteSpace($OutCsv)) {
    $OutCsv = Join-Path $WorkRoot 'pbe_tabi_convergence_shape_summary.csv'
}

function As-Double([object]$Value) {
    [double]::Parse([string]$Value, $InvariantCulture)
}

function As-Long([object]$Value) {
    [long]::Parse([string]$Value, $InvariantCulture)
}

function Get-Percentile([double[]]$Values, [double]$Fraction) {
    $items = @($Values)
    if ($items.Count -eq 0) { return 0.0 }
    $ordered = @($items | Sort-Object)
    $rank = [Math]::Ceiling($Fraction * $ordered.Count) - 1
    $rank = [Math]::Max(0, [Math]::Min($ordered.Count - 1, $rank))
    return $ordered[$rank]
}

function Get-Mean([double[]]$Values) {
    $items = @($Values)
    if ($items.Count -eq 0) { return 0.0 }
    return ($items | Measure-Object -Average).Average
}

function Get-Min([double[]]$Values) {
    $items = @($Values)
    if ($items.Count -eq 0) { return 0.0 }
    return ($items | Measure-Object -Minimum).Minimum
}

function Get-Max([double[]]$Values) {
    $items = @($Values)
    if ($items.Count -eq 0) { return 0.0 }
    return ($items | Measure-Object -Maximum).Maximum
}

function Get-PanelSum([object[]]$Rows) {
    $sum = 0L
    foreach ($row in $Rows) {
        $sum += [long]$row.node_panels
    }
    return $sum
}

function Format-Double([double]$Value) {
    return $Value.ToString('G17', $InvariantCulture)
}

$channels = @('global', 'K1_phi', 'K2_phi_normal', 'K3_phin', 'K4_mixed')
$features = @(
    'rho',
    'rho_power',
    'aspect_ratio',
    'source_norm',
    'tail_norm_p',
    'tail_norm_pm1',
    'tail_ratio',
    'effective_amp',
    'estimated_rel_error',
    'estimated_abs_error',
    'oracle_abs_over_source_norm',
    'oracle_abs_over_tail_norm',
    'oracle_abs_over_tail_rho'
)

$allRows = @()
foreach ($fixture in $Fixtures) {
    $csv = Join-Path $WorkRoot (Join-Path $fixture (Join-Path $CaseId 'pbe_tabi_convergence_shape_calibration.csv'))
    if (-not (Test-Path -LiteralPath $csv)) {
        throw "Missing Stage 3AH shape artifact for ${fixture}: $csv"
    }
    foreach ($row in @(Import-Csv -LiteralPath $csv)) {
        $record = [ordered]@{
            fixture = $fixture
            channel = [string]$row.channel
            classification = [string]$row.classification
            node_panels = As-Long $row.node_panels
        }
        foreach ($feature in $features) {
            $record[$feature] = As-Double $row.$feature
        }
        $allRows += [PSCustomObject]$record
    }
}

if ($allRows.Count -eq 0) {
    throw 'No convergence shape calibration rows were loaded.'
}

$matrixCsv = Join-Path $WorkRoot 'pbe_tabi_topdown_physical_gap_matrix.csv'
$safeLeafPanels = 0.0
if (Test-Path -LiteralPath $matrixCsv) {
    $safeRows = @(Import-Csv -LiteralPath $matrixCsv | Where-Object {
        $_.case_id -eq 'stage3ag_safe_baseline_audit_off' -and
        $Fixtures -contains $_.fixture -and $_.status -eq 'ok'
    })
    foreach ($row in $safeRows) {
        $safeLeafPanels += As-Double $row.leaf_panel_work_per_matvec
    }
}

$summary = @()
$bestProjection = $null
foreach ($channel in $channels) {
    $channelRows = @($allRows | Where-Object { $_.channel -eq $channel })
    $caught = @($channelRows | Where-Object { $_.classification -eq 'caught' })
    $overreject = @($channelRows | Where-Object { $_.classification -eq 'overreject' })
    foreach ($feature in $features) {
        $caughtValues = [double[]]@($caught | ForEach-Object { $_.$feature })
        $overValues = [double[]]@($overreject | ForEach-Object { $_.$feature })
        $minCaught = Get-Min $caughtValues
        $maxOver = Get-Max $overValues
        $threshold = if ($caughtValues.Count -gt 0) { $minCaught } else { [double]::NaN }
        $safeOver = if ($caughtValues.Count -gt 0) {
            @($overreject | Where-Object { $_.$feature -lt $threshold })
        } else {
            @()
        }
        $safePanels = Get-PanelSum $safeOver
        $safeFraction = if ($safeLeafPanels -gt 0.0) { $safePanels / $safeLeafPanels } else { 0.0 }
        $separability = if ($maxOver -gt 0.0) { [double]$minCaught / [double]$maxOver } else { [double]::NaN }
        $row = [PSCustomObject]@{
            channel = $channel
            feature = $feature
            caught_count = $caught.Count
            overreject_count = $overreject.Count
            caught_median = Format-Double (Get-Percentile $caughtValues 0.50)
            caught_p95 = Format-Double (Get-Percentile $caughtValues 0.95)
            caught_p99 = Format-Double (Get-Percentile $caughtValues 0.99)
            overreject_median = Format-Double (Get-Percentile $overValues 0.50)
            overreject_p95 = Format-Double (Get-Percentile $overValues 0.95)
            overreject_p99 = Format-Double (Get-Percentile $overValues 0.99)
            min_caught = Format-Double $minCaught
            max_overreject = Format-Double $maxOver
            separability_min_caught_over_max_overreject = if ([double]::IsNaN($separability)) { '' } else { Format-Double $separability }
            overreject_below_min_caught = @($safeOver).Count
            overreject_below_min_caught_panels = $safePanels
            overreject_below_min_caught_frac_vs_safe_leaf_work = Format-Double $safeFraction
            safe_threshold_exclusive = if ([double]::IsNaN($threshold)) { '' } else { Format-Double $threshold }
            next_branch_recommendation = ''
        }
        if ($safeFraction -ge 0.10 -and $caught.Count -gt 0) {
            if ($null -eq $bestProjection -or
                $safeFraction -gt [double]$bestProjection.overreject_below_min_caught_frac_vs_safe_leaf_work) {
                $bestProjection = $row
            }
        }
        $summary += $row
    }
}

$recommendation =
    if ($null -ne $bestProjection) {
        "implement-estimator-v2:$($bestProjection.channel):$($bestProjection.feature)"
    } else {
        "revisit-channel-source-moment-semantics"
    }

$summary = $summary | ForEach-Object {
    $_.next_branch_recommendation = $recommendation
    $_
}

$outDir = Split-Path -Parent $OutCsv
if (-not [string]::IsNullOrWhiteSpace($outDir)) {
    New-Item -ItemType Directory -Path $outDir -Force | Out-Null
}
$summary | Export-Csv -LiteralPath $OutCsv -NoTypeInformation -Encoding ASCII

Write-Host "PBE TABI convergence shape summary COMPLETE: $OutCsv"
Write-Host "  recommendation=$recommendation"
if ($null -ne $bestProjection) {
    Write-Host "  best channel=$($bestProjection.channel) feature=$($bestProjection.feature) projected_fraction=$($bestProjection.overreject_below_min_caught_frac_vs_safe_leaf_work)"
}
