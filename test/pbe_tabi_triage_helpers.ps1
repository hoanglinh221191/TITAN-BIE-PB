$TABI_TRIAGE_HELPERS_LOADED = $true

$Script:TabiTriageInvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Convert-TabiTriageDoubleOrNull([object]$Value) {
    if ($null -eq $Value) { return $null }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return $null }
    $parsed = 0.0
    if ([double]::TryParse($text, [System.Globalization.NumberStyles]::Float,
            $Script:TabiTriageInvariantCulture, [ref]$parsed)) {
        return $parsed
    }
    return $null
}

function Format-TabiTriageNumber([object]$Value) {
    if ($null -eq $Value) { return 'NA' }
    if ($Value -is [string]) { return $Value }
    return ([double]$Value).ToString('G17', $Script:TabiTriageInvariantCulture)
}

function Get-TabiTriageValue($Object, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Object -and $Object.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Object.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}

function Get-TabiProgressRows([string]$RunDir) {
    if ([string]::IsNullOrWhiteSpace($RunDir) -or $RunDir -eq 'NA') {
        return @()
    }
    $path = Join-Path $RunDir 'pbe_tabi_longrun_progress.csv'
    if (-not (Test-Path -LiteralPath $path)) {
        return @()
    }
    return @(Import-Csv -LiteralPath $path)
}

function Get-TabiPrepareProgressRows([string]$RunDir) {
    if ([string]::IsNullOrWhiteSpace($RunDir) -or $RunDir -eq 'NA') {
        return @()
    }
    $path = Join-Path $RunDir 'pbe_tabi_prepare_progress.csv'
    if (-not (Test-Path -LiteralPath $path)) {
        return @()
    }
    return @(Import-Csv -LiteralPath $path)
}

function Get-LastTabiProgressRow([string]$RunDir) {
    $rows = @(Get-TabiProgressRows -RunDir $RunDir)
    if ($rows.Count -eq 0) { return $null }
    return $rows[-1]
}

function Get-LastTabiPrepareProgressRow([string]$RunDir) {
    $rows = @(Get-TabiPrepareProgressRows -RunDir $RunDir)
    if ($rows.Count -eq 0) { return $null }
    return $rows[-1]
}

function Get-TabiTriageRequestedPreconditioner($Row, $ProgressRow) {
    foreach ($name in @('preconditioner_requested','tabi_preconditioner_requested')) {
        $value = Get-TabiTriageValue $Row $name
        if ($value -ne 'NA') { return $value }
    }
    $effective = Get-TabiTriageValue $ProgressRow 'preconditioner_effective'
    if ($effective -like 'diagonal_block*' -or $effective -like 'block*') {
        return 'block'
    }
    if ($effective -like 'diagonal_coarse*' -or $effective -like 'coarse*') {
        return 'coarse'
    }
    return 'NA'
}

function Get-TabiTriageEffectivePreconditioner($Row, $ProgressRow) {
    foreach ($name in @('preconditioner_effective','tabi_preconditioner_effective')) {
        $value = Get-TabiTriageValue $Row $name
        if ($value -ne 'NA') { return $value }
    }
    return Get-TabiTriageValue $ProgressRow 'preconditioner_effective'
}

function Test-TabiPreconditionerDowngraded([string]$Requested, [string]$Effective) {
    if ([string]::IsNullOrWhiteSpace($Requested) -or $Requested -eq 'NA') {
        return $false
    }
    if ([string]::IsNullOrWhiteSpace($Effective) -or $Effective -eq 'NA') {
        return $false
    }
    $requestedLower = $Requested.ToLowerInvariant()
    $effectiveLower = $Effective.ToLowerInvariant()
    $requestedStrong = $requestedLower -match 'block|coarse|budget'
    return ($requestedStrong -and $effectiveLower.StartsWith('diagonal_'))
}

function Get-TabiPreconditionerMetric($Row, $ProgressRow, [string[]]$Names, [string]$Default = 'NA') {
    foreach ($name in $Names) {
        $value = Get-TabiTriageValue $ProgressRow $name
        if ($value -ne 'NA') { return $value }
    }
    foreach ($name in $Names) {
        $value = Get-TabiTriageValue $Row $name
        if ($value -ne 'NA') { return $value }
    }
    return $Default
}

function Get-TabiResidualDeltaPerMatvec($ProgressRows, $FallbackRow, [double]$MatvecCalls) {
    $gmresRows = @($ProgressRows | Where-Object {
        (Get-TabiTriageValue $_ 'phase') -like 'gmres_*' -and
        $null -ne (Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $_ 'residual'))
    })
    if ($gmresRows.Count -ge 2) {
        $first = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $gmresRows[0] 'residual')
        $last = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $gmresRows[-1] 'residual')
        $firstMatvec = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $gmresRows[0] 'matvec_calls')
        $lastMatvec = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $gmresRows[-1] 'matvec_calls')
        $denom = if ($null -ne $firstMatvec -and $null -ne $lastMatvec -and
                     $lastMatvec -gt $firstMatvec) {
            $lastMatvec - $firstMatvec
        } elseif ($MatvecCalls -gt 0.0) {
            $MatvecCalls
        } else {
            1.0
        }
        return ($first - $last) / $denom
    }
    $best = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $FallbackRow 'best_residual')
    $final = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $FallbackRow 'final_residual')
    if ($null -ne $best -and $null -ne $final -and $MatvecCalls -gt 0.0) {
        return ($best - $final) / $MatvecCalls
    }
    return $null
}

function Get-TabiSecondsPerMatvec($ProgressRows, $FallbackRow, [double]$MatvecCalls) {
    if ($MatvecCalls -le 0.0) { return $null }
    $gmresStart = @($ProgressRows | Where-Object { (Get-TabiTriageValue $_ 'phase') -eq 'gmres_start' } | Select-Object -First 1)
    $lastGmres = @($ProgressRows | Where-Object { (Get-TabiTriageValue $_ 'phase') -like 'gmres_*' } | Select-Object -Last 1)
    if ($gmresStart.Count -gt 0 -and $lastGmres.Count -gt 0) {
        $startSeconds = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $gmresStart[0] 'elapsed_seconds')
        $lastSeconds = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $lastGmres[0] 'elapsed_seconds')
        if ($null -ne $startSeconds -and $null -ne $lastSeconds -and $lastSeconds -gt $startSeconds) {
            return ($lastSeconds - $startSeconds) / $MatvecCalls
        }
    }
    $totalSeconds = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $FallbackRow 'total_seconds')
    if ($null -ne $totalSeconds -and $totalSeconds -gt 0.0) {
        return $totalSeconds / $MatvecCalls
    }
    return $null
}

function New-TabiTriageRecord($Row, [string]$RunDir, [string]$ExistingClassification = 'NA') {
    $progressRows = @(Get-TabiProgressRows -RunDir $RunDir)
    $progress = if ($progressRows.Count -gt 0) { $progressRows[-1] } else { $null }
    $prepareRows = @(Get-TabiPrepareProgressRows -RunDir $RunDir)
    $prepare = if ($prepareRows.Count -gt 0) { $prepareRows[-1] } else { $null }
    $mode = Get-TabiTriageValue $Row 'mode'
    if ($mode -eq 'grid') {
        return [PSCustomObject]@{
            triage_classification = if ($ExistingClassification -ne 'NA') { $ExistingClassification } else { 'grid_reference' }
            seconds_per_matvec = 'NA'
            leaf_panel_interactions_per_matvec = 'NA'
            leaf_vs_n2_ratio = 'NA'
            residual_delta_per_matvec = 'NA'
            preconditioner_requested = 'NA'
            preconditioner_effective = 'grid'
            preconditioner_downgrade = 'false'
            preconditioner_block_panel_cap = 'NA'
            preconditioner_memory_budget_bytes = 'NA'
            preconditioner_estimated_bytes = 'NA'
            preconditioner_valid_blocks = 'NA'
            preconditioner_valid_panels = 'NA'
            preconditioner_ill_conditioned_blocks = 'NA'
            preconditioner_quality_rejected_blocks = 'NA'
            preconditioner_quality_probe_max_gain = 'NA'
            preconditioner_quality_status = 'NA'
            preconditioner_setup_seconds = 'NA'
            preconditioner_apply_seconds = 'NA'
            preconditioner_apply_calls = 'NA'
            leaf_mac_candidates = 'NA'
            leaf_mac_accepted = 'NA'
            leaf_mac_rejected = 'NA'
            leaf_mac_far_panels = 'NA'
            leaf_mac_direct_panels_saved_estimate = 'NA'
            prepare_progress_phase = 'grid'
            prepare_progress_status = 'grid'
            prepare_progress_rss_bytes = 'NA'
            triage_source = 'grid'
        }
    }

    $matvec = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $progress 'matvec_calls')
    if ($null -eq $matvec) {
        $matvec = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $Row 'matvec_calls')
    }
    if ($null -eq $matvec) { $matvec = 0.0 }

    $leaf = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $progress 'leaf_panel_interactions')
    if ($null -eq $leaf) {
        $leaf = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $Row 'leaf_panel_interactions')
    }
    $panels = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $progress 'panels')
    if ($null -eq $panels) {
        $panels = Convert-TabiTriageDoubleOrNull (Get-TabiTriageValue $Row 'panels')
    }
    $secondsPerMatvec = Get-TabiSecondsPerMatvec -ProgressRows $progressRows -FallbackRow $Row -MatvecCalls $matvec
    $leafPerMatvec = if ($null -ne $leaf -and $matvec -gt 0.0) { $leaf / $matvec } else { $null }
    $leafVsN2 = if ($null -ne $leaf -and $null -ne $panels -and $panels -gt 0.0) {
        $leaf / [Math]::Max(1.0, $panels * $panels)
    } else {
        $null
    }
    $residualDelta = Get-TabiResidualDeltaPerMatvec -ProgressRows $progressRows -FallbackRow $Row -MatvecCalls $matvec
    $requested = Get-TabiTriageRequestedPreconditioner -Row $Row -ProgressRow $progress
    $effective = Get-TabiTriageEffectivePreconditioner -Row $Row -ProgressRow $progress
    $downgraded = Test-TabiPreconditionerDowngraded -Requested $requested -Effective $effective
    $precondBlockPanelCap = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_block_panel_cap','tabi_preconditioner_block_panel_cap')
    $precondMemoryBudget = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_memory_budget_bytes','tabi_preconditioner_memory_budget_bytes')
    $precondEstimatedBytes = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_estimated_bytes','tabi_preconditioner_estimated_bytes')
    $precondValidBlocks = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_valid_blocks','tabi_preconditioner_valid_blocks')
    $precondValidPanels = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_valid_panels','tabi_preconditioner_valid_panels')
    $precondIllConditioned = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_ill_conditioned_blocks','tabi_preconditioner_ill_conditioned_blocks')
    $precondQualityRejected = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_quality_rejected_blocks','tabi_preconditioner_quality_rejected_blocks')
    $precondQualityGain = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_quality_probe_max_gain','tabi_preconditioner_quality_probe_max_gain')
    $precondQualityStatus = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_quality_status','tabi_preconditioner_quality_status')
    $precondSetupSeconds = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_setup_seconds','tabi_preconditioner_setup_seconds')
    $precondApplySeconds = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_apply_seconds','tabi_preconditioner_apply_seconds')
    $precondApplyCalls = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('preconditioner_apply_calls','tabi_preconditioner_apply_calls')
    $leafMacCandidates = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('leaf_mac_candidates','tabi_leaf_mac_candidates') -Default '0'
    $leafMacAccepted = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('leaf_mac_accepted','tabi_leaf_mac_accepted') -Default '0'
    $leafMacRejected = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('leaf_mac_rejected','tabi_leaf_mac_rejected') -Default '0'
    $leafMacFarPanels = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('leaf_mac_far_panels','tabi_leaf_mac_far_panels') -Default '0'
    $leafMacSavedPanels = Get-TabiPreconditionerMetric -Row $Row -ProgressRow $progress -Names @('leaf_mac_direct_panels_saved_estimate','tabi_leaf_mac_direct_panels_saved_estimate') -Default '0'
    $status = Get-TabiTriageValue $Row 'status'
    $classification = if ($ExistingClassification -ne 'NA') {
        $ExistingClassification
    } else {
        Get-TabiTriageValue $Row 'performance_classification'
    }

    $stalled = $false
    if ($null -ne $residualDelta) {
        $stalled = [Math]::Abs($residualDelta) -le 1.0e-12
    }
    if ($downgraded) {
        $classification = 'preconditioner_downgraded'
    } elseif ($null -ne $leafVsN2 -and $leafVsN2 -gt 1.0) {
        $classification = 'leaf_work_super_n2'
    } elseif ($stalled -and $matvec -gt 0.0) {
        $classification = 'convergence_stalled'
    } elseif ($status -eq 'timeout' -and $matvec -gt 0.0) {
        $classification = 'matvec_limited'
    } elseif ($status -eq 'timeout' -and $progressRows.Count -eq 0 -and
              $prepareRows.Count -gt 0) {
        $preparePhase = Get-TabiTriageValue $prepare 'phase'
        if ($preparePhase -eq 'solver_enter') {
            $classification = 'solver_entry_no_longrun_profile'
        } elseif ($preparePhase -eq 'tabi_options_built') {
            $classification = 'solver_options_built_timeout'
        } elseif ($preparePhase -eq 'pqr_loaded') {
            $classification = 'prepare_after_pqr_timeout'
        } else {
            $classification = 'prepare_limited'
        }
    } elseif ($status -eq 'timeout' -and $progressRows.Count -eq 0 -and
              $classification -eq 'timeout_no_profile' -and
              $prepareRows.Count -gt 0) {
        $classification = 'solver_entry_no_longrun_profile'
    } elseif ($classification -eq 'NA') {
        $classification = 'unclassified'
    }

    [PSCustomObject]@{
        triage_classification = $classification
        seconds_per_matvec = Format-TabiTriageNumber $secondsPerMatvec
        matvec_calls = Format-TabiTriageNumber $matvec
        leaf_panel_interactions_per_matvec = Format-TabiTriageNumber $leafPerMatvec
        leaf_vs_n2_ratio = Format-TabiTriageNumber $leafVsN2
        residual_delta_per_matvec = Format-TabiTriageNumber $residualDelta
        preconditioner_requested = $requested
        preconditioner_effective = $effective
        preconditioner_downgrade = if ($downgraded) { 'true' } else { 'false' }
        preconditioner_block_panel_cap = $precondBlockPanelCap
        preconditioner_memory_budget_bytes = $precondMemoryBudget
        preconditioner_estimated_bytes = $precondEstimatedBytes
        preconditioner_valid_blocks = $precondValidBlocks
        preconditioner_valid_panels = $precondValidPanels
        preconditioner_ill_conditioned_blocks = $precondIllConditioned
        preconditioner_quality_rejected_blocks = $precondQualityRejected
        preconditioner_quality_probe_max_gain = $precondQualityGain
        preconditioner_quality_status = $precondQualityStatus
        preconditioner_setup_seconds = $precondSetupSeconds
        preconditioner_apply_seconds = $precondApplySeconds
        preconditioner_apply_calls = $precondApplyCalls
        leaf_mac_candidates = $leafMacCandidates
        leaf_mac_accepted = $leafMacAccepted
        leaf_mac_rejected = $leafMacRejected
        leaf_mac_far_panels = $leafMacFarPanels
        leaf_mac_direct_panels_saved_estimate = $leafMacSavedPanels
        prepare_progress_phase = Get-TabiTriageValue $prepare 'phase'
        prepare_progress_status = Get-TabiTriageValue $prepare 'status'
        prepare_progress_rss_bytes = Get-TabiTriageValue $prepare 'rss_bytes'
        triage_source = if ($progressRows.Count -gt 0) {
            'progress'
        } elseif ($prepareRows.Count -gt 0) {
            'prepare_progress'
        } else {
            'matrix'
        }
    }
}
