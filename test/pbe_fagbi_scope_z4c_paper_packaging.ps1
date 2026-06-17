param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [switch]$RefreshZ4B,
    [int]$CaseTimeoutSeconds = 900,
    [string]$Z4BWorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4b'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4c')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MainTableCsv = Join-Path $WorkRoot 'fagbi_scope_z4c_main_table.csv'
$SupplementTableCsv = Join-Path $WorkRoot 'fagbi_scope_z4c_supplement_table.csv'
$ReferencePolicyCsv = Join-Path $WorkRoot 'fagbi_scope_z4c_reference_policy_table.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4c_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4c_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z4c_runtime_memory_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4C_PAPER_PACKAGE_REPORT.txt'

if ($RefreshZ4B) {
    $z4bScript = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4b_reference_policy_gate.ps1'
    if (-not (Test-Path -LiteralPath $z4bScript)) {
        throw "Missing Z4B refresh script: $z4bScript"
    }
    $z4bRootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4B_REFERENCE_POLICY_REPORT.txt'
    $z4bReportSnapshot = $null
    $z4bReportExisted = Test-Path -LiteralPath $z4bRootReport
    if ($z4bReportExisted) {
        $z4bReportSnapshot = Get-Content -LiteralPath $z4bRootReport
    }
    $z4bArgs = @{
        TitanExe = $TitanExe
        CaseTimeoutSeconds = $CaseTimeoutSeconds
        WorkRoot = $Z4BWorkRoot
    }
    if ($Quick) { $z4bArgs['Quick'] = $true }
    try {
        & $z4bScript @z4bArgs
    } finally {
        if ($z4bReportExisted) {
            Set-Content -LiteralPath $z4bRootReport -Value $z4bReportSnapshot -Encoding ASCII
        } elseif (Test-Path -LiteralPath $z4bRootReport) {
            Remove-Item -LiteralPath $z4bRootReport -Force
        }
    }
}

$z4bMatrixCsv = Join-Path $Z4BWorkRoot 'fagbi_scope_z4b_matrix.csv'
$z4bReferenceCsv = Join-Path $Z4BWorkRoot 'fagbi_scope_z4b_reference_policy.csv'
$z4bDecisionCsv = Join-Path $Z4BWorkRoot 'fagbi_scope_z4b_decision.csv'
foreach ($required in @($z4bMatrixCsv, $z4bReferenceCsv, $z4bDecisionCsv)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing Z4B artifact: $required. Run test\pbe_fagbi_scope_z4b_reference_policy_gate.ps1 without -Quick before packaging, or pass -RefreshZ4B."
    }
}

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Z4C work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z4cDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or $text -eq 'nan') {
        return $Default
    }
    return [double]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Test-Z4cBool([object]$Value) {
    return ([string]$Value).Equals('true', [System.StringComparison]::OrdinalIgnoreCase)
}

function Assert-Z4cColumns([object[]]$Rows, [string[]]$Columns, [string]$Name) {
    if ($Rows.Count -eq 0) {
        throw "$Name has no rows."
    }
    $available = @($Rows[0].PSObject.Properties.Name)
    foreach ($column in $Columns) {
        if ($available -notcontains $column) {
            throw "$Name missing required column: $column"
        }
    }
}

function Get-Z4cDecision([object[]]$MainRows, [bool]$QuickMode) {
    if ($MainRows.Count -eq 0) { return 'z4c_reference_policy_gap' }
    if (-not $QuickMode -and $MainRows.Count -ne 6) {
        return 'z4c_reference_policy_gap'
    }
    if (@($MainRows | Where-Object { -not (Test-Z4cBool $_.paper_usable) }).Count -gt 0) {
        return 'z4c_reference_policy_gap'
    }
    foreach ($row in $MainRows) {
        $cap = [int]$row.atom_cap
        if ($cap -lt 256 -and $row.reference_kind -ne 'full_direct_tabi_energy') {
            return 'z4c_reference_policy_gap'
        }
        if ($cap -ge 256 -and $row.reference_kind -ne 'sampled_direct_operator_audit') {
            return 'z4c_reference_policy_gap'
        }
        if ($cap -lt 256) {
            $rel = Convert-Z4cDouble $row.direct_relative_energy_delta
            $delta = [Math]::Abs((Convert-Z4cDouble $row.direct_energy_delta_kj_mol))
            if (-not (Test-FmmCliFinite $rel) -or ($rel -gt 0.05 -and $delta -gt 2.0)) {
                return 'z4c_reference_policy_gap'
            }
        }
        if ($cap -ge 256) {
            $audit = Convert-Z4cDouble $row.sampled_audit_rms
            if (-not (Test-FmmCliFinite $audit) -or $audit -gt 8.0e-3) {
                return 'z4c_reference_policy_gap'
            }
        }
    }
    if (-not (Test-Path -LiteralPath $PlotPath)) {
        return 'z4c_plot_gap'
    }
    return 'z4c_paper_package_ready'
}

function Write-Z4cSvg([object[]]$MainRows, [object[]]$SupplementRows, [string]$Path) {
    $allRows = @($MainRows + $SupplementRows | Where-Object {
            (Convert-Z4cDouble $_.wall_seconds) -gt 0.0
        })
    $maxWall = ($allRows | ForEach-Object { Convert-Z4cDouble $_.wall_seconds } |
        Where-Object { Test-FmmCliFinite $_ } | Measure-Object -Maximum).Maximum
    $maxRss = ($allRows | ForEach-Object { Convert-Z4cDouble $_.peak_working_set_mb } |
        Where-Object { Test-FmmCliFinite $_ } | Measure-Object -Maximum).Maximum
    if (-not (Test-FmmCliFinite $maxWall) -or $maxWall -le 0.0) { $maxWall = 1.0 }
    if (-not (Test-FmmCliFinite $maxRss) -or $maxRss -le 0.0) { $maxRss = 1.0 }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="1040" height="420" viewBox="0 0 1040 420">') | Out-Null
    $svg.Add('<rect width="1040" height="420" fill="white"/>') | Out-Null
    $svg.Add('<text x="24" y="28" font-family="Arial" font-size="16">Scope Z4C paper package: runtime and memory context</text>') | Out-Null
    $svg.Add('<text x="24" y="48" font-family="Arial" font-size="11">Wall-time bars (solid) and peak RSS ticks (black) for paper FMM rows plus context rows</text>') | Out-Null
    $plotBase = 330.0
    $x = 36
    foreach ($row in $allRows) {
        $wall = Convert-Z4cDouble $row.wall_seconds 0.0
        $rss = Convert-Z4cDouble $row.peak_working_set_mb 0.0
        $h = [Math]::Max(2.0, 240.0 * $wall / $maxWall)
        $tickY = $plotBase - [Math]::Max(2.0, 240.0 * $rss / $maxRss)
        $fill = switch ($row.backend) {
            'fmm_ilu0' { '#0072b2' }
            'tabi_direct_full' { '#d55e00' }
            'grid_context' { '#009e73' }
            default { '#999999' }
        }
        $svg.Add(('<rect x="{0}" y="{1:F1}" width="20" height="{2:F1}" fill="{3}"/>' -f $x, ($plotBase - $h), $h, $fill)) | Out-Null
        $svg.Add(('<line x1="{0}" y1="{1:F1}" x2="{2}" y2="{1:F1}" stroke="black" stroke-width="2"/>' -f ($x - 2), $tickY, ($x + 22))) | Out-Null
        $label = ('{0} c{1} {2}' -f $row.fixture, $row.atom_cap, $row.backend)
        $svg.Add(('<text x="{0}" y="354" font-family="Arial" font-size="8" transform="rotate(45 {0},354)">{1}</text>' -f $x, $label)) | Out-Null
        $x += 32
        if ($x -gt 990) { break }
    }
    $svg.Add('<rect x="790" y="24" width="12" height="12" fill="#0072b2"/><text x="808" y="34" font-family="Arial" font-size="10">FMM ILU0 paper row</text>') | Out-Null
    $svg.Add('<rect x="790" y="44" width="12" height="12" fill="#d55e00"/><text x="808" y="54" font-family="Arial" font-size="10">direct TABI reference</text>') | Out-Null
    $svg.Add('<rect x="790" y="64" width="12" height="12" fill="#009e73"/><text x="808" y="74" font-family="Arial" font-size="10">grid context</text>') | Out-Null
    $svg.Add('</svg>') | Out-Null
    Set-Content -LiteralPath $Path -Value $svg -Encoding ASCII
}

$matrixRows = @(Import-Csv -LiteralPath $z4bMatrixCsv)
$referenceRows = @(Import-Csv -LiteralPath $z4bReferenceCsv)
$decisionRows = @(Import-Csv -LiteralPath $z4bDecisionCsv)
Assert-Z4cColumns -Rows $matrixRows -Name 'Z4B matrix' -Columns @(
    'fixture','atom_cap','backend','status','failure_class','energy_kj_mol',
    'reference_kind','reference_backend','energy_delta_kj_mol',
    'relative_energy_delta','reference_usable_for_paper',
    'tree_vs_direct_matvec_rel','matvec_audit_rows','wall_seconds',
    'peak_working_set_mb')
Assert-Z4cColumns -Rows $referenceRows -Name 'Z4B reference policy' -Columns @(
    'fixture','atom_cap','reference_kind','reference_backend',
    'reference_status','reference_usable_for_paper',
    'tree_vs_direct_matvec_rel','matvec_audit_rows')

$mainRows = @($matrixRows | Where-Object { $_.backend -eq 'fmm_ilu0' } |
    ForEach-Object {
        [pscustomobject]@{
            fixture = $_.fixture
            atom_cap = [int]$_.atom_cap
            backend = $_.backend
            status = $_.status
            failure_class = $_.failure_class
            energy_kj_mol = $_.energy_kj_mol
            reference_kind = $_.reference_kind
            reference_backend = $_.reference_backend
            direct_energy_delta_kj_mol = if ([int]$_.atom_cap -lt 256) { $_.energy_delta_kj_mol } else { 'nan' }
            direct_relative_energy_delta = if ([int]$_.atom_cap -lt 256) { $_.relative_energy_delta } else { 'nan' }
            sampled_audit_rms = $_.tree_vs_direct_matvec_rel
            matvec_audit_rows = $_.matvec_audit_rows
            wall_seconds = $_.wall_seconds
            peak_working_set_mb = $_.peak_working_set_mb
            paper_usable = $_.reference_usable_for_paper
        }
    })

$supplementRows = @($matrixRows | Where-Object { $_.backend -ne 'fmm_ilu0' } |
    ForEach-Object {
        $role = if ($_.backend -eq 'grid_context') {
            if ($_.status -eq 'ok') { 'grid_context_ok' } else { 'grid_context_only_not_converged' }
        } elseif ($_.backend -eq 'tabi_direct_full') {
            'full_direct_tabi_reference'
        } else {
            'context'
        }
        [pscustomobject]@{
            fixture = $_.fixture
            atom_cap = [int]$_.atom_cap
            backend = $_.backend
            table_role = $role
            status = $_.status
            failure_class = $_.failure_class
            energy_kj_mol = $_.energy_kj_mol
            wall_seconds = $_.wall_seconds
            peak_working_set_mb = $_.peak_working_set_mb
            paper_oracle = ($_.backend -eq 'tabi_direct_full' -and $_.status -eq 'ok')
        }
    })

$referencePolicyRows = @($referenceRows | ForEach-Object {
        [pscustomobject]@{
            fixture = $_.fixture
            atom_cap = [int]$_.atom_cap
            reference_kind = $_.reference_kind
            reference_backend = $_.reference_backend
            reference_status = $_.reference_status
            paper_usable = $_.reference_usable_for_paper
            energy_delta_kj_mol = $_.energy_delta_kj_mol
            relative_energy_delta = $_.relative_energy_delta
            sampled_audit_rms = $_.tree_vs_direct_matvec_rel
            matvec_audit_rows = $_.matvec_audit_rows
        }
    })

Write-Z4cSvg -MainRows $mainRows -SupplementRows $supplementRows -Path $PlotPath
$decision = Get-Z4cDecision -MainRows $mainRows -QuickMode ([bool]$Quick)
$maxAudit = ($mainRows | ForEach-Object { Convert-Z4cDouble $_.sampled_audit_rms } |
    Where-Object { Test-FmmCliFinite $_ } | Measure-Object -Maximum).Maximum
if (-not (Test-FmmCliFinite $maxAudit)) { $maxAudit = [double]::NaN }
$gridContextWarnings = @($supplementRows | Where-Object {
        $_.backend -eq 'grid_context' -and $_.status -ne 'ok'
    }).Count

$mainRows | Export-Csv -LiteralPath $MainTableCsv -NoTypeInformation
$supplementRows | Export-Csv -LiteralPath $SupplementTableCsv -NoTypeInformation
$referencePolicyRows | Export-Csv -LiteralPath $ReferencePolicyCsv -NoTypeInformation
[pscustomobject]@{
    decision = $decision
    quick = [bool]$Quick
    refresh_z4b = [bool]$RefreshZ4B
    z4b_decision = if ($decisionRows.Count -gt 0) { $decisionRows[0].decision } else { '' }
    main_rows = $mainRows.Count
    supplement_rows = $supplementRows.Count
    max_sampled_audit = $maxAudit
    grid_context_warnings = $gridContextWarnings
    artifact_main_table = $MainTableCsv
    artifact_supplement_table = $SupplementTableCsv
    artifact_reference_policy_table = $ReferencePolicyCsv
    artifact_plot = $PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'scope: Z4C paper-ready protein panel packaging',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "refresh_z4b: $([bool]$RefreshZ4B)",
    "z4b_decision: $(if ($decisionRows.Count -gt 0) { $decisionRows[0].decision } else { '' })",
    "main_rows: $($mainRows.Count)",
    "supplement_rows: $($supplementRows.Count)",
    "max_sampled_audit: $maxAudit",
    "grid_context_warnings: $gridContextWarnings",
    'reference_policy: cap-128 FMM rows cite full direct TABI energy; cap-256 FMM rows cite sampled direct-operator audit',
    'default_policy_note: grid remains production default; Z4C does not change public behavior'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z4C Paper Package Report',
    '',
    "Decision: $decision",
    '',
    'Z4C packages the Z4B reference-policy gate into paper-ready tables and a compact runtime/memory figure. It does not run solver rows unless -RefreshZ4B is explicitly passed.',
    '',
    'Artifacts:',
    "- Main table: $MainTableCsv",
    "- Supplement table: $SupplementTableCsv",
    "- Reference policy table: $ReferencePolicyCsv",
    "- Decision: $DecisionCsv",
    "- Summary: $SummaryPath",
    "- Plot: $PlotPath",
    '',
    'Main FMM rows:'
)
foreach ($row in $mainRows) {
    $report += ('- fixture={0}, cap={1}, energy={2}, ref_kind={3}, direct_delta={4}, audit={5}, wall_s={6}, rss_mb={7}, paper_usable={8}' -f `
        $row.fixture, $row.atom_cap, $row.energy_kj_mol, $row.reference_kind, `
        $row.direct_energy_delta_kj_mol, $row.sampled_audit_rms, `
        $row.wall_seconds, $row.peak_working_set_mb, $row.paper_usable)
}
$report += ''
$report += 'Supplement/context rows:'
foreach ($row in $supplementRows) {
    $report += ('- fixture={0}, cap={1}, backend={2}, role={3}, status={4}, energy={5}, wall_s={6}, rss_mb={7}' -f `
        $row.fixture, $row.atom_cap, $row.backend, $row.table_role, `
        $row.status, $row.energy_kj_mol, $row.wall_seconds, `
        $row.peak_working_set_mb)
}
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "PASS Scope Z4C paper package: $decision"
