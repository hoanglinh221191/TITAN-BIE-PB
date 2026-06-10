param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_leaf_kernel_transition_probe'),
    [Alias('Fixture')]
    [string[]]$Fixtures = @('barnase','research-medium','research-large'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 90,
    [int[]]$SurfacePointCounts = @(24),
    [int[]]$TreeLeafSizes = @(16,32,64,128),
    [double[]]$OpeningThetas = @(0.45,0.6,0.8,1.0),
    [ValidateSet('baseline','deferred-worklist','target-leaf-block','soa-avx2','target-avx2-v2','cost-tile-4x4','leaf-mac')]
    [string[]]$SchedulerVariants = @('baseline','deferred-worklist','target-leaf-block','soa-avx2','target-avx2-v2','cost-tile-4x4','leaf-mac'),
    [int]$Threads = 8,
    [string]$Preconditioner = 'block-budgeted',
    [string]$PreconditionerBlockNear = 'operator-near',
    [int]$PreconditionerBlockPanelCap = 32,
    [int]$PreconditionerMemoryBudgetMb = 256,
    [switch]$ProfileMemory
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
$RepoRoot = Split-Path -Parent $PSScriptRoot
$PerfMatrixScript = Join-Path $PSScriptRoot 'pbe_tabi_perf_matrix.ps1'
$TriageHelperScript = Join-Path $PSScriptRoot 'pbe_tabi_triage_helpers.ps1'
$OutCsv = Join-Path $WorkRoot 'pbe_tabi_leaf_kernel_transition_probe.csv'
$SummaryPath = Join-Path $WorkRoot 'pbe_tabi_leaf_kernel_transition_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_TABI_LEAF_KERNEL_OPTIMIZATION_REPORT.txt'
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI leaf-kernel probe root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Convert-ToDoubleOrNull([object]$Value) {
    if ($null -eq $Value) { return $null }
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

function Get-SortDoubleOrLarge([object]$Value) {
    $parsed = Convert-ToDoubleOrNull $Value
    if ($null -eq $parsed -or [double]::IsNaN($parsed) -or [double]::IsInfinity($parsed)) {
        return [double]::PositiveInfinity
    }
    return $parsed
}

function Get-RowValue($Object, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Object -and $Object.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Object.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}

function Get-FirstCsvRow([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    $rows = @(Import-Csv -LiteralPath $Path)
    if ($rows.Count -eq 0) { return $null }
    return $rows[0]
}

function Get-LatestCsvValue([string]$Path, [string[]]$Names) {
    if (-not (Test-Path -LiteralPath $Path)) { return '' }
    $rows = @(Import-Csv -LiteralPath $Path)
    for ($i = $rows.Count - 1; $i -ge 0; --$i) {
        foreach ($name in $Names) {
            $value = Get-RowValue $rows[$i] $name ''
            if (-not [string]::IsNullOrWhiteSpace($value) -and $value -ne 'NA') {
                return $value
            }
        }
    }
    return ''
}

function Get-RunMetric($Row, [string[]]$Names, [string]$Default = 'NA') {
    foreach ($name in $Names) {
        $value = Get-RowValue $Row $name ''
        if (-not [string]::IsNullOrWhiteSpace($value) -and $value -ne 'NA') {
            return $value
        }
    }
    $runDir = Get-RowValue $Row 'run_dir'
    if ($runDir -ne 'NA') {
        foreach ($fileName in @('pbe_timing.csv','pbe_tabi_longrun_progress.csv')) {
            $path = Join-Path $runDir $fileName
            $latestValue = Get-LatestCsvValue -Path $path -Names $Names
            if (-not [string]::IsNullOrWhiteSpace($latestValue)) {
                return $latestValue
            }
        }
    }
    return $Default
}

function Get-SchedulerExtraArgs([string]$Scheduler) {
    switch ($Scheduler) {
        'baseline' { return @() }
        'deferred-worklist' {
            return @(
                '-pbe-tabi-leaf-pair-scheduler=deferred-worklist',
                '-pbe-tabi-leaf-pair-worklist-max-active=2048',
                '-pbe-tabi-leaf-pair-worklist-verify=off',
                '-pbe-tabi-leaf-kernel-profile=summary'
            )
        }
        'target-leaf-block' {
            return @(
                '-pbe-tabi-leaf-pair-scheduler=target-leaf-block',
                '-pbe-tabi-leaf-block-verify=off',
                '-pbe-tabi-leaf-kernel-profile=summary'
            )
        }
        'soa-avx2' {
            return @(
                '-pbe-tabi-leaf-kernel-profile=summary',
                '-pbe-tabi-jcp-direct-kernel=soa-avx2'
            )
        }
        'target-avx2-v2' {
            return @(
                '-pbe-tabi-interaction-mode=onthefly',
                '-pbe-tabi-cluster-policy=leaf-pair-balanced',
                '-pbe-tabi-leaf-pair-kernel=target-avx2-v2',
                '-pbe-tabi-target-batch-size=4',
                '-pbe-tabi-leaf-pair-min-source-panels=64',
                '-pbe-tabi-leaf-pair-min-targets=4',
                '-pbe-tabi-leaf-kernel-profile=summary'
            )
        }
        'cost-tile-4x4' {
            return @(
                '-pbe-tabi-interaction-mode=onthefly',
                '-pbe-tabi-cluster-policy=leaf-pair-balanced',
                '-pbe-tabi-leaf-pair-scheduler=cost-tile-worklist',
                '-pbe-tabi-leaf-pair-kernel=target-avx2-4x4',
                '-pbe-tabi-target-batch-size=4',
                '-pbe-tabi-leaf-pair-min-source-panels=32',
                '-pbe-tabi-leaf-pair-min-targets=4',
                '-pbe-tabi-leaf-block-verify=off',
                '-pbe-tabi-leaf-kernel-profile=summary'
            )
        }
        'leaf-mac' {
            return @(
                '-pbe-tabi-leaf-mac=on',
                '-pbe-tabi-accepted-miss-diagnostics=summary',
                '-pbe-tabi-expansion-audit=sampled',
                '-pbe-tabi-expansion-audit-channel=summary',
                '-pbe-tabi-leaf-work-profile=summary',
                '-pbe-tabi-leaf-kernel-profile=summary'
            )
        }
    }
    throw "Unknown scheduler variant: $Scheduler"
}

function Get-SafeName([string]$Text) {
    return ($Text -replace '[^A-Za-z0-9_-]', '_')
}

function Get-MatrixRows([string]$Scheduler, [int[]]$LeafSizes, [double[]]$Thetas) {
    $caseRoot = Join-Path $WorkRoot (Get-SafeName $Scheduler)
    $matrixArgs = @{
        TitanExe = $TitanExe
        WorkRoot = $caseRoot
        Fixtures = $runFixtures
        Modes = @('tabi_treecode')
        SurfacePointCounts = $runSurfacePoints
        TreeLeafSizes = $LeafSizes
        OpeningThetas = $Thetas
        LeafNearPolicies = @('exact')
        LeafPanelKernels = @('tile-avx2')
        AdaptiveLevels = @('strict')
        ClusterPolicies = @('off')
        InteractionModes = @('cached')
        Threads = @([string]$Threads)
        TimingMode = 'detail'
        TabiProgressIntervalSeconds = 5
        CaseTimeoutSeconds = $timeout
        Preconditioner = $Preconditioner
        PreconditionerBlockNear = $PreconditionerBlockNear
        PreconditionerBlockPanelCap = $PreconditionerBlockPanelCap
        PreconditionerMemoryBudgetMb = $PreconditionerMemoryBudgetMb
        TabiExtraArgs = @(Get-SchedulerExtraArgs $Scheduler)
        AllowFailedStatus = $true
    }
    if ($ProfileMemory) { $matrixArgs.ProfileMemory = $true }
    & $PerfMatrixScript @matrixArgs
    $matrixCsv = Join-Path $caseRoot 'pbe_tabi_perf_matrix.csv'
    if (-not (Test-Path -LiteralPath $matrixCsv)) {
        throw "Missing TABI leaf-kernel matrix output: $matrixCsv"
    }
    foreach ($row in @(Import-Csv -LiteralPath $matrixCsv)) {
        $triage = New-TabiTriageRecord -Row $row -RunDir (Get-RowValue $row 'run_dir') -ExistingClassification (Get-RowValue $row 'performance_classification')
        $acceptedMiss = Get-RunMetric $row @(
            'tabi_accepted_missed_total',
            'tabi_accepted_miss_total',
            'accepted_miss',
            'accepted_missed_total'
        )
        $verifyMismatch = Get-RunMetric $row @(
            'tabi_leaf_pair_worklist_verify_mismatches',
            'tabi_leaf_block_verify_mismatches',
            'tabi_leaf_pair_directness_verify_mismatches',
            'leaf_pair_worklist_verify_mismatches',
            'leaf_block_verify_mismatches'
        ) '0'
        $verifyStatus = 'clean'
        $acceptedNumber = Convert-ToDoubleOrNull $acceptedMiss
        $mismatchNumber = Convert-ToDoubleOrNull $verifyMismatch
        if (($null -ne $acceptedNumber -and $acceptedNumber -gt 0.0) -or
            ($null -ne $mismatchNumber -and $mismatchNumber -gt 0.0)) {
            $verifyStatus = 'mismatch'
        } elseif ($acceptedMiss -eq 'NA' -and $verifyMismatch -eq 'NA') {
            $verifyStatus = 'unknown'
        }
        $resolvedKernel = Get-RowValue $row 'leaf_panel_kernel'
        if ($resolvedKernel -eq 'NA') { $resolvedKernel = Get-RowValue $row 'requested_leaf_panel_kernel' }
        $leafPairBatchKernel = Get-RunMetric $row @(
            'tabi_leaf_pair_batch_kernel',
            'leaf_pair_batch_kernel'
        )
        $leafPairVectorSources = Get-RunMetric $row @(
            'tabi_leaf_pair_batch_vectorized_sources',
            'leaf_pair_batch_vectorized_sources'
        )
        $leafPairVectorLanes = Get-RunMetric $row @(
            'tabi_leaf_pair_batch_vector_lanes',
            'leaf_pair_batch_vector_lanes'
        )
        $leafBlockSourceGroups = Get-RunMetric $row @(
            'tabi_leaf_block_source_groups',
            'leaf_block_source_groups'
        )
        $leafPairBatchSourcePanels = Get-RunMetric $row @(
            'tabi_leaf_pair_batch_source_panels',
            'leaf_pair_batch_source_panels'
        )
        $leafPairBatchTargets = Get-RunMetric $row @(
            'tabi_leaf_pair_batch_targets',
            'leaf_pair_batch_targets'
        )
        $leafBlockFlushSeconds = Get-RunMetric $row @(
            'tabi_leaf_block_flush_seconds',
            'leaf_block_flush_seconds'
        )
        $leafMacCandidates = Get-RunMetric $row @(
            'tabi_leaf_mac_candidates',
            'leaf_mac_candidates'
        ) '0'
        $leafMacAccepted = Get-RunMetric $row @(
            'tabi_leaf_mac_accepted',
            'leaf_mac_accepted'
        ) '0'
        $leafMacRejected = Get-RunMetric $row @(
            'tabi_leaf_mac_rejected',
            'leaf_mac_rejected'
        ) '0'
        $leafMacFarPanels = Get-RunMetric $row @(
            'tabi_leaf_mac_far_panels',
            'leaf_mac_far_panels'
        ) '0'
        $leafMacSavedPanels = Get-RunMetric $row @(
            'tabi_leaf_mac_direct_panels_saved_estimate',
            'leaf_mac_direct_panels_saved_estimate'
        ) '0'
        $kernel4x4Sources = '0'
        if ($leafPairBatchKernel -eq 'target-avx2-4x4') {
            $kernel4x4Sources = $leafPairVectorSources
        }
        [PSCustomObject]@{
            fixture = Get-RowValue $row 'fixture'
            scheduler = $Scheduler
            status = Get-RowValue $row 'status'
            surface_points = Get-RowValue $row 'surface_points'
            leaf_size = Get-RowValue $row 'tree_leaf_size'
            theta = Get-RowValue $row 'opening_theta'
            leaf_panel_kernel = $resolvedKernel
            leaf_pair_batch_kernel = $leafPairBatchKernel
            leaf_pair_batch_vectorized_sources = $leafPairVectorSources
            leaf_pair_batch_vector_lanes = $leafPairVectorLanes
            leaf_pair_cost_tile_tasks = if ($Scheduler -eq 'cost-tile-4x4') { $leafBlockSourceGroups } else { '0' }
            leaf_pair_cost_tile_source_panels = if ($Scheduler -eq 'cost-tile-4x4') { $leafPairBatchSourcePanels } else { '0' }
            leaf_pair_cost_tile_target_panels = if ($Scheduler -eq 'cost-tile-4x4') { $leafPairBatchTargets } else { '0' }
            leaf_pair_cost_tile_seconds = if ($Scheduler -eq 'cost-tile-4x4') { $leafBlockFlushSeconds } else { '0' }
            leaf_pair_kernel_4x4_sources = $kernel4x4Sources
            leaf_mac_candidates = $leafMacCandidates
            leaf_mac_accepted = $leafMacAccepted
            leaf_mac_rejected = $leafMacRejected
            leaf_mac_far_panels = $leafMacFarPanels
            leaf_mac_direct_panels_saved_estimate = $leafMacSavedPanels
            leaf_near_policy = Get-RowValue $row 'leaf_near_policy'
            adaptive_level = Get-RowValue $row 'adaptive_level'
            cluster_policy = Get-RowValue $row 'cluster_policy'
            seconds_per_matvec = $triage.seconds_per_matvec
            matvec_calls = $triage.matvec_calls
            total_seconds = Get-RowValue $row 'total_seconds'
            leaf_panel_interactions = Get-RowValue $row 'leaf_panel_interactions'
            leaf_panel_interactions_per_matvec = if ($triage.leaf_panel_interactions_per_matvec -ne 'NA') { $triage.leaf_panel_interactions_per_matvec } else { Get-RowValue $row 'leaf_panel_per_matvec' }
            leaf_vs_n2_ratio = $triage.leaf_vs_n2_ratio
            leaf_panel_pairs_per_second = Get-RowValue $row 'leaf_panel_pairs_per_second'
            peak_ram_bytes = Get-RowValue $row 'peak_ram_bytes'
            accepted_miss = $acceptedMiss
            verify_mismatches = $verifyMismatch
            verify_status = $verifyStatus
            preconditioner_effective = $triage.preconditioner_effective
            preconditioner_quality_status = $triage.preconditioner_quality_status
            preconditioner_setup_seconds = $triage.preconditioner_setup_seconds
            triage_classification = $triage.triage_classification
            performance_classification = Get-RowValue $row 'performance_classification'
            decision_class = 'pending'
            promote_candidate = 'false'
            run_dir = Get-RowValue $row 'run_dir'
            source_matrix_csv = $matrixCsv
            scheduler_extra_args = ((Get-SchedulerExtraArgs $Scheduler) -join ';')
        }
    }
}

function Set-DecisionClasses($Rows) {
    $baselineByFixture = @{}
    foreach ($fixture in @($Rows | Select-Object -ExpandProperty fixture -Unique)) {
        $baseline = @($Rows | Where-Object {
            $_.fixture -eq $fixture -and
            $_.scheduler -eq 'baseline' -and
            $_.leaf_size -eq '128' -and
            $_.theta -eq '0.80000000000000004'
        } | Select-Object -First 1)
        if ($baseline.Count -eq 0) {
            $baseline = @($Rows | Where-Object {
                $_.fixture -eq $fixture -and
                $_.scheduler -eq 'baseline' -and
                $_.leaf_size -eq '128'
            } | Select-Object -First 1)
        }
        if ($baseline.Count -gt 0) { $baselineByFixture[$fixture] = $baseline[0] }
    }
    foreach ($row in $Rows) {
        $baseline = if ($baselineByFixture.ContainsKey($row.fixture)) { $baselineByFixture[$row.fixture] } else { $null }
        $rowSeconds = Convert-ToDoubleOrNull $row.seconds_per_matvec
        $rowLeaf = Convert-ToDoubleOrNull $row.leaf_panel_interactions_per_matvec
        $rowRam = Convert-ToDoubleOrNull $row.peak_ram_bytes
        $baseSeconds = if ($null -ne $baseline) { Convert-ToDoubleOrNull $baseline.seconds_per_matvec } else { $null }
        $baseLeaf = if ($null -ne $baseline) { Convert-ToDoubleOrNull $baseline.leaf_panel_interactions_per_matvec } else { $null }
        $baseRam = if ($null -ne $baseline) { Convert-ToDoubleOrNull $baseline.peak_ram_bytes } else { $null }
        $accepted = Convert-ToDoubleOrNull $row.accepted_miss
        $mismatch = Convert-ToDoubleOrNull $row.verify_mismatches
        $clean = (($null -eq $accepted -or $accepted -eq 0.0) -and
                  ($null -eq $mismatch -or $mismatch -eq 0.0) -and
                  $row.verify_status -ne 'mismatch')
        if (-not $clean) {
            $row.decision_class = 'mac_accuracy_fail'
            continue
        }
        if ($null -ne $rowRam -and $null -ne $baseRam -and $baseRam -gt 0.0 -and
            $rowRam -gt (2.0 * $baseRam)) {
            $row.decision_class = 'memory_limited'
            continue
        }
        $treeWin = $false
        if ($null -ne $rowLeaf -and $null -ne $baseLeaf -and $baseLeaf -gt 0.0 -and
            $rowLeaf -le ($baseLeaf / 3.0)) {
            $treeWin = $true
        }
        $schedulerWin = $false
        if ($null -ne $rowSeconds -and $null -ne $baseSeconds -and $baseSeconds -gt 0.0 -and
            $rowSeconds -le ($baseSeconds / 1.5) -and
            ($null -eq $rowLeaf -or $null -eq $baseLeaf -or $rowLeaf -le (1.2 * $baseLeaf))) {
            $schedulerWin = $true
        }
        if ($treeWin) {
            $row.decision_class = 'tree_depth_win'
        } elseif ($schedulerWin) {
            $row.decision_class = 'scheduler_win'
        } else {
            $row.decision_class = 'no_win'
        }
        if (($row.fixture -eq 'barnase') -and
            (($null -ne $rowSeconds -and $rowSeconds -lt 3.0) -or $treeWin)) {
            $row.promote_candidate = 'true'
        }
    }
}

function Write-LeafKernelReport($Rows) {
    $candidateRows = @($Rows | Where-Object { $_.promote_candidate -eq 'true' } |
        Sort-Object @{Expression = { Get-SortDoubleOrLarge $_.seconds_per_matvec }; Ascending = $true },
                    @{Expression = { Get-SortDoubleOrLarge $_.leaf_panel_interactions_per_matvec }; Ascending = $true })
    $bestRow = $null
    if ($candidateRows.Count -gt 0) {
        $bestRow = $candidateRows[0]
    } else {
        $fallbackRows = @($Rows | Where-Object { $_.fixture -eq 'barnase' } |
            Sort-Object @{Expression = { Get-SortDoubleOrLarge $_.seconds_per_matvec }; Ascending = $true },
                        @{Expression = { Get-SortDoubleOrLarge $_.leaf_panel_interactions_per_matvec }; Ascending = $true } |
            Select-Object -First 1)
        if ($fallbackRows.Count -gt 0) { $bestRow = $fallbackRows[0] }
    }
    $lines = New-Object System.Collections.Generic.List[string]
    $lines.Add('TITAN TABI Leaf-Kernel Optimization Report')
    $lines.Add('==========================================')
    $lines.Add('')
    $lines.Add("quick_mode: $([string][bool]$Quick)")
    $lines.Add("fixtures: $($runFixtures -join ',')")
    $lines.Add("case_timeout_seconds: $timeout")
    $lines.Add("preconditioner: $Preconditioner")
    $lines.Add("preconditioner_block_near: $PreconditionerBlockNear")
    $lines.Add('production_default_note: grid remains the production PBE default; TABI remains experimental.')
    $lines.Add('theta_policy_note: theta/MAC direction is treated empirically; only clean measured rows can be promoted.')
    $lines.Add('')
    $lines.Add('Baseline')
    $lines.Add('--------')
    $lines.Add('baseline_row: leaf_size=128 theta=0.8 scheduler=baseline leaf_near_policy=exact leaf_panel_kernel=tile-avx2 adaptive=strict cluster=off')
    $lines.Add('')
    $lines.Add('Decision Counts')
    $lines.Add('---------------')
    foreach ($group in @($Rows | Group-Object decision_class | Sort-Object Name)) {
        $lines.Add(("decision_{0}: {1}" -f $group.Name, $group.Count))
    }
    $lines.Add('')
    $lines.Add('Top Barnase Rows')
    $lines.Add('----------------')
    $lines.Add('| Scheduler | leaf_size | theta | status | sec/matvec | leaf_work/matvec | verify | decision | promote |')
    $lines.Add('| --- | ---: | ---: | --- | ---: | ---: | --- | --- | --- |')
    foreach ($row in @($Rows | Where-Object { $_.fixture -eq 'barnase' } |
            Sort-Object @{Expression = { Get-SortDoubleOrLarge $_.seconds_per_matvec }; Ascending = $true },
                        @{Expression = { Get-SortDoubleOrLarge $_.leaf_panel_interactions_per_matvec }; Ascending = $true } |
            Select-Object -First 12)) {
        $lines.Add(('| {0} | {1} | {2} | {3} | {4} | {5} | {6} | {7} | {8} |' -f
            $row.scheduler, $row.leaf_size, $row.theta, $row.status,
            $row.seconds_per_matvec, $row.leaf_panel_interactions_per_matvec,
            $row.verify_status, $row.decision_class, $row.promote_candidate))
    }
    $lines.Add('')
    $lines.Add('Leaf-Pair Kernel Rows')
    $lines.Add('---------------------')
    $lines.Add('| Scheduler | leaf_size | theta | batch_kernel | vector_sources | vector_lanes | cost_tile_tasks | cost_tile_seconds |')
    $lines.Add('| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |')
    foreach ($row in @($Rows | Where-Object { $_.fixture -eq 'barnase' } |
            Select-Object -First 12)) {
        $lines.Add(('| {0} | {1} | {2} | {3} | {4} | {5} | {6} | {7} |' -f
            $row.scheduler, $row.leaf_size, $row.theta,
            $row.leaf_pair_batch_kernel,
            $row.leaf_pair_batch_vectorized_sources,
            $row.leaf_pair_batch_vector_lanes,
            $row.leaf_pair_cost_tile_tasks,
            $row.leaf_pair_cost_tile_seconds))
    }
    $lines.Add('')
    $lines.Add('Leaf-MAC Rows')
    $lines.Add('-------------')
    $lines.Add('| Scheduler | leaf_size | theta | candidates | accepted | far_panels | saved_panels | decision |')
    $lines.Add('| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |')
    foreach ($row in @($Rows | Where-Object { $_.fixture -eq 'barnase' -and $_.scheduler -eq 'leaf-mac' } |
            Select-Object -First 12)) {
        $lines.Add(('| {0} | {1} | {2} | {3} | {4} | {5} | {6} | {7} |' -f
            $row.scheduler, $row.leaf_size, $row.theta,
            $row.leaf_mac_candidates, $row.leaf_mac_accepted,
            $row.leaf_mac_far_panels,
            $row.leaf_mac_direct_panels_saved_estimate,
            $row.decision_class))
    }
    $leafMacNoMatvec = @($Rows | Where-Object {
        $_.fixture -eq 'barnase' -and $_.scheduler -eq 'leaf-mac' -and
        ($_.seconds_per_matvec -eq 'NA' -or $_.matvec_calls -eq '0' -or $_.matvec_calls -eq 'NA')
    }).Count
    if ($leafMacNoMatvec -gt 0) {
        $lines.Add('')
        $lines.Add('leaf_mac_note: one or more Barnase Leaf-MAC rows timed out before producing a nonzero matvec/progress work counter; treat these rows as setup/first-matvec stalls, not as evidence that no leaf MAC candidates exist.')
    }
    $lines.Add('')
    $lines.Add('Recommendation')
    $lines.Add('--------------')
    if ($candidateRows.Count -gt 0) {
        $lines.Add(('recommendation: promote the best clean candidate into Phase 7 experimental gate: fixture={0}, scheduler={1}, leaf_size={2}, theta={3}, sec_per_matvec={4}, leaf_work_per_matvec={5}.' -f
            $candidateRows[0].fixture, $candidateRows[0].scheduler,
            $candidateRows[0].leaf_size, $candidateRows[0].theta,
            $candidateRows[0].seconds_per_matvec,
            $candidateRows[0].leaf_panel_interactions_per_matvec))
    } elseif ($null -ne $bestRow) {
        $hasLeafMac = @($Rows | Where-Object { $_.scheduler -eq 'leaf-mac' }).Count -gt 0
        $hasCostTile4x4 = @($Rows | Where-Object { $_.scheduler -eq 'cost-tile-4x4' }).Count -gt 0
        if ($hasLeafMac) {
            $lines.Add(('recommendation: no clean promote candidate yet; best Barnase row is scheduler={0}, leaf_size={1}, theta={2}, sec_per_matvec={3}, decision={4}. Leaf-MAC was exercised as a traversal/work-count fix; inspect leaf_mac_saved_panels and accuracy diagnostics before returning to micro-kernel work.' -f
                $bestRow.scheduler, $bestRow.leaf_size, $bestRow.theta,
                $bestRow.seconds_per_matvec, $bestRow.decision_class))
        } elseif ($hasCostTile4x4) {
            $lines.Add(('recommendation: no clean promote candidate yet; best Barnase row is scheduler={0}, leaf_size={1}, theta={2}, sec_per_matvec={3}, decision={4}. V2B cost-tile-4x4 was exercised and remains no_win, so the next decision should be algorithmic ROI or targeted profiling of 4x4 overhead, not another wrapper-only tuning pass.' -f
                $bestRow.scheduler, $bestRow.leaf_size, $bestRow.theta,
                $bestRow.seconds_per_matvec, $bestRow.decision_class))
        } else {
            $lines.Add(('recommendation: no clean promote candidate yet; best Barnase row is scheduler={0}, leaf_size={1}, theta={2}, sec_per_matvec={3}, decision={4}. Next optimization is V2B 4x4/4x8 register blocking plus static cost-tile worklist, not more wrapper tuning.' -f
                $bestRow.scheduler, $bestRow.leaf_size, $bestRow.theta,
                $bestRow.seconds_per_matvec, $bestRow.decision_class))
        }
    } else {
        $lines.Add('recommendation: no Barnase rows were produced; inspect script/runtime failures before changing solver code.')
    }
    $lines.Add('')
    $lines.Add("probe_csv: $OutCsv")
    $lines.Add("summary_dat: $SummaryPath")
    Set-Content -LiteralPath $ReportPath -Value $lines -Encoding ASCII
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
if (-not (Test-Path -LiteralPath $PerfMatrixScript)) { throw "Missing TABI perf matrix script: $PerfMatrixScript" }
if (-not (Test-Path -LiteralPath $TriageHelperScript)) { throw "Missing TABI triage helper script: $TriageHelperScript" }
. $TriageHelperScript

New-CleanDir $WorkRoot

$runFixtures = if ($Quick) { @($Fixtures | Select-Object -First 1) } else { $Fixtures }
$timeout = if ($Quick -and -not $PSBoundParameters.ContainsKey('CaseTimeoutSeconds')) { 90 } else { $CaseTimeoutSeconds }
$runSurfacePoints = if ($Quick) { @(24) } else { $SurfacePointCounts }
$runLeafSizes = if ($Quick) { @(64,128) } else { $TreeLeafSizes }
$runThetas = if ($Quick) { @(0.8) } else { $OpeningThetas }
$runSchedulers = if ($Quick) { @('baseline','target-avx2-v2','cost-tile-4x4','leaf-mac') } else { $SchedulerVariants }

$rows = New-Object System.Collections.Generic.List[object]
foreach ($scheduler in $runSchedulers) {
    foreach ($row in Get-MatrixRows -Scheduler $scheduler -LeafSizes $runLeafSizes -Thetas $runThetas) {
        $rows.Add($row)
    }
}

Set-DecisionClasses $rows
$rows.ToArray() | Export-Csv -LiteralPath $OutCsv -NoTypeInformation -Encoding ASCII

$summary = @(
    'title: TABI Leaf-Kernel Transition Probe',
    "quick_mode: $([string][bool]$Quick)",
    "fixtures: $($runFixtures -join ',')",
    "case_timeout_seconds: $timeout",
    "threads: $Threads",
    "surface_points: $($runSurfacePoints -join ',')",
    "leaf_sizes: $($runLeafSizes -join ',')",
    "opening_thetas: $($runThetas -join ',')",
    "schedulers: $($runSchedulers -join ',')",
    "rows: $($rows.Count)",
    'note: diagnostics are benchmark-only and do not change production solver defaults'
)
foreach ($group in @($rows | Group-Object decision_class | Sort-Object Name)) {
    $summary += "decision_$($group.Name): $($group.Count)"
}
$summary += "probe_csv: $OutCsv"
$summary += "report_txt: $ReportPath"
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII
Write-LeafKernelReport $rows

Write-Host "PBE TABI leaf-kernel transition probe COMPLETE: $OutCsv"
