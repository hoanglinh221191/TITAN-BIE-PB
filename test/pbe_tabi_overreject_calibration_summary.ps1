param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_topdown_physical_gap_matrix'),
    [string]$OutCsv = '',
    [string[]]$Fixtures = @('mini-protein', 'charged-peptide', 'compact-mini-protein'),
    [string]$CaseId = 'stage3ag_overreject_calibration_detail',
    [double]$AbsoluteFloor = 1.0e-7
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

if ([string]::IsNullOrWhiteSpace($OutCsv)) {
    $OutCsv = Join-Path $WorkRoot 'pbe_tabi_overreject_calibration_summary.csv'
}

function As-Double([object]$Value) {
    [double]::Parse([string]$Value, $InvariantCulture)
}

function As-Long([object]$Value) {
    [long]::Parse([string]$Value, $InvariantCulture)
}

function Get-Mean([double[]]$Values) {
    if ($Values.Count -eq 0) { return 0.0 }
    return ($Values | Measure-Object -Average).Average
}

function Get-Max([double[]]$Values) {
    if ($Values.Count -eq 0) { return 0.0 }
    return ($Values | Measure-Object -Maximum).Maximum
}

function Get-Percentile([double[]]$Values, [double]$Fraction) {
    if ($Values.Count -eq 0) { return 0.0 }
    $ordered = @($Values | Sort-Object)
    $rank = [Math]::Ceiling($Fraction * $ordered.Count) - 1
    $rank = [Math]::Max(0, [Math]::Min($ordered.Count - 1, $rank))
    return $ordered[$rank]
}

function Get-PanelSum([object[]]$Rows) {
    $sum = 0L
    foreach ($row in $Rows) {
        $sum += [long]$row.node_panels
    }
    return $sum
}

$scales = @(1.0e-6, 3.0e-6, 1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4,
            1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2, 1.0e-1, 3.0e-1, 1.0)
$channels = @('global', 'K1_phi', 'K2_phi_normal', 'K3_phin', 'K4_mixed')
$allRows = @()

foreach ($fixture in $Fixtures) {
    $csv = Join-Path $WorkRoot (Join-Path $fixture (Join-Path $CaseId 'pbe_tabi_convergence_overreject_calibration.csv'))
    if (-not (Test-Path -LiteralPath $csv)) {
        throw "Missing Stage 3AG detail artifact for ${fixture}: $csv"
    }
    foreach ($row in @(Import-Csv -LiteralPath $csv)) {
        $allRows += [PSCustomObject]@{
            fixture = $fixture
            channel = [string]$row.channel
            classification = [string]$row.classification
            node_panels = As-Long $row.node_panels
            estimated_abs_error = As-Double $row.estimated_abs_error
            oracle_abs_error = As-Double $row.oracle_abs_error
            required_scale = As-Double $row.required_scale
        }
    }
}

if ($allRows.Count -eq 0) {
    throw 'No overreject calibration detail rows were loaded.'
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

$recommended = @{}
foreach ($channel in $channels) {
    $channelRows = @($allRows | Where-Object { $_.channel -eq $channel })
    $caught = @($channelRows | Where-Object { $_.classification -eq 'caught' })
    $safeScale = $null
    foreach ($scale in $scales) {
        $wouldReacceptCaught = @($caught | Where-Object {
            $_.estimated_abs_error * $scale -le $AbsoluteFloor
        })
        if ($wouldReacceptCaught.Count -eq 0) {
            $safeScale = $scale
            break
        }
    }
    if ($null -ne $safeScale) {
        $safeIndex = [Array]::IndexOf($scales, [double]$safeScale)
        $marginIndex = [Math]::Min($safeIndex + 1, $scales.Count - 1)
        $recommended[$channel] = $scales[$marginIndex]
    } else {
        $recommended[$channel] = [double]::NaN
    }
}

$summary = @()
foreach ($channel in $channels) {
    $channelRows = @($allRows | Where-Object { $_.channel -eq $channel })
    $caught = @($channelRows | Where-Object { $_.classification -eq 'caught' })
    $overreject = @($channelRows | Where-Object { $_.classification -eq 'overreject' })
    $oracleValues = [double[]]@($channelRows | ForEach-Object { $_.oracle_abs_error })
    $estimatedValues = [double[]]@($channelRows | ForEach-Object { $_.estimated_abs_error })
    $requiredValues = [double[]]@($channelRows | ForEach-Object { $_.required_scale })
    foreach ($scale in $scales) {
        $caughtWouldReaccept = @($caught | Where-Object {
            $_.estimated_abs_error * $scale -le $AbsoluteFloor
        })
        $overrejectWouldReaccept = @($overreject | Where-Object {
            $_.estimated_abs_error * $scale -le $AbsoluteFloor
        })
        $caughtPanels = Get-PanelSum $caughtWouldReaccept
        $overrejectPanels = Get-PanelSum $overrejectWouldReaccept
        $summary += [PSCustomObject]@{
            channel = $channel
            scale = $scale.ToString('G17', $InvariantCulture)
            absolute_floor = $AbsoluteFloor.ToString('G17', $InvariantCulture)
            sampled_count = $channelRows.Count
            caught_count = $caught.Count
            overreject_count = $overreject.Count
            caught_would_reaccept = $caughtWouldReaccept.Count
            caught_would_reaccept_panels = $caughtPanels
            overreject_would_reaccept = $overrejectWouldReaccept.Count
            overreject_would_reaccept_panels = $overrejectPanels
            caught_would_reaccept_frac_vs_safe_leaf_work = if ($safeLeafPanels -gt 0.0) { ($caughtPanels / $safeLeafPanels).ToString('G17', $InvariantCulture) } else { '' }
            overreject_would_reaccept_frac_vs_safe_leaf_work = if ($safeLeafPanels -gt 0.0) { ($overrejectPanels / $safeLeafPanels).ToString('G17', $InvariantCulture) } else { '' }
            safe_on_sample = ($caughtWouldReaccept.Count -eq 0)
            oracle_abs_mean = (Get-Mean $oracleValues).ToString('G17', $InvariantCulture)
            oracle_abs_max = (Get-Max $oracleValues).ToString('G17', $InvariantCulture)
            oracle_abs_p95 = (Get-Percentile $oracleValues 0.95).ToString('G17', $InvariantCulture)
            oracle_abs_p99 = (Get-Percentile $oracleValues 0.99).ToString('G17', $InvariantCulture)
            estimated_abs_mean = (Get-Mean $estimatedValues).ToString('G17', $InvariantCulture)
            estimated_abs_max = (Get-Max $estimatedValues).ToString('G17', $InvariantCulture)
            required_scale_max = (Get-Max $requiredValues).ToString('G17', $InvariantCulture)
            required_scale_p95 = (Get-Percentile $requiredValues 0.95).ToString('G17', $InvariantCulture)
            required_scale_p99 = (Get-Percentile $requiredValues 0.99).ToString('G17', $InvariantCulture)
            recommended_scale_with_margin = if ([double]::IsNaN($recommended[$channel])) { '' } else { ([double]$recommended[$channel]).ToString('G17', $InvariantCulture) }
        }
    }
}

$outDir = Split-Path -Parent $OutCsv
if (-not [string]::IsNullOrWhiteSpace($outDir)) {
    New-Item -ItemType Directory -Path $outDir -Force | Out-Null
}
$summary | Export-Csv -LiteralPath $OutCsv -NoTypeInformation -Encoding ASCII

Write-Host "PBE TABI overreject calibration summary COMPLETE: $OutCsv"
foreach ($channel in $channels) {
    $value = $recommended[$channel]
    if ([double]::IsNaN($value)) {
        Write-Host "  $channel`: no sampled-safe scale"
    } else {
        Write-Host "  $channel`: recommended scale with margin=$($value.ToString('G17', $InvariantCulture))"
    }
}
