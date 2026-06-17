param(
    [switch]$RefreshZ4C,
    [string]$Z4CWorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4c'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4d')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MainCsv = Join-Path $WorkRoot 'fagbi_scope_z4d_main_table_paper.csv'
$MainMd = Join-Path $WorkRoot 'fagbi_scope_z4d_main_table_paper.md'
$SupplementCsv = Join-Path $WorkRoot 'fagbi_scope_z4d_supplement_table_paper.csv'
$MethodsMd = Join-Path $WorkRoot 'fagbi_scope_z4d_methods_box.md'
$CaptionMd = Join-Path $WorkRoot 'fagbi_scope_z4d_figure_caption.md'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4d_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4d_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4D_MANUSCRIPT_EXPORT_REPORT.txt'

if ($RefreshZ4C) {
    $z4cScript = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4c_paper_packaging.ps1'
    if (-not (Test-Path -LiteralPath $z4cScript)) {
        throw "Missing Z4C packaging script: $z4cScript"
    }
    & $z4cScript -WorkRoot $Z4CWorkRoot
}

$z4cMainCsv = Join-Path $Z4CWorkRoot 'fagbi_scope_z4c_main_table.csv'
$z4cSupplementCsv = Join-Path $Z4CWorkRoot 'fagbi_scope_z4c_supplement_table.csv'
$z4cReferenceCsv = Join-Path $Z4CWorkRoot 'fagbi_scope_z4c_reference_policy_table.csv'
foreach ($required in @($z4cMainCsv, $z4cSupplementCsv, $z4cReferenceCsv)) {
    if (-not (Test-Path -LiteralPath $required)) {
        throw "Missing Z4C artifact: $required. Run test\pbe_fagbi_scope_z4c_paper_packaging.ps1 before Z4D."
    }
}

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Z4D work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z4dDouble([object]$Value, [double]$Default = [double]::NaN) {
    if ($null -eq $Value) { return $Default }
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA' -or $text -eq 'nan') {
        return $Default
    }
    return [double]::Parse($text, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Test-Z4dFinite([double]$Value) {
    return (-not [double]::IsNaN($Value)) -and (-not [double]::IsInfinity($Value))
}

function Test-Z4dBool([object]$Value) {
    return ([string]$Value).Equals('true', [System.StringComparison]::OrdinalIgnoreCase)
}

function Format-Z4dNumber([object]$Value, [string]$Format = 'F3') {
    $number = Convert-Z4dDouble $Value
    if (-not (Test-Z4dFinite $number)) { return 'NA' }
    return $number.ToString($Format, [System.Globalization.CultureInfo]::InvariantCulture)
}

function Assert-Z4dColumns([object[]]$Rows, [string[]]$Columns, [string]$Name) {
    if ($Rows.Count -eq 0) { throw "$Name has no rows." }
    $available = @($Rows[0].PSObject.Properties.Name)
    foreach ($column in $Columns) {
        if ($available -notcontains $column) {
            throw "$Name missing required column: $column"
        }
    }
}

function Convert-Z4dReferenceLabel([string]$ReferenceKind) {
    switch ($ReferenceKind) {
        'full_direct_tabi_energy' { return 'Direct TABI energy' }
        'sampled_direct_operator_audit' { return 'Sampled direct-operator audit' }
        default { return $ReferenceKind }
    }
}

function Write-Z4dMarkdownTable([object[]]$Rows, [string]$Path) {
    $lines = New-Object System.Collections.Generic.List[string]
    $lines.Add('| Fixture | Cap | E_rxn (kJ/mol) | Reference | Delta E (kJ/mol) | Audit RMS | Wall (s) | Peak RSS (MB) |') | Out-Null
    $lines.Add('| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |') | Out-Null
    foreach ($row in $Rows) {
        $lines.Add(('| {0} | {1} | {2} | {3} | {4} | {5} | {6} | {7} |' -f `
            $row.fixture, $row.atom_cap, $row.energy_kj_mol_display, `
            $row.reference_label, $row.direct_energy_delta_kj_mol_display, `
            $row.sampled_audit_rms_display, $row.wall_seconds_display, `
            $row.peak_working_set_mb_display)) | Out-Null
    }
    Set-Content -LiteralPath $Path -Value $lines -Encoding ASCII
}

function Get-Z4dDecision([object[]]$Rows, [string]$MethodsText, [string]$CaptionText) {
    if ($Rows.Count -ne 6) { return 'z4d_traceability_gap' }
    if (@($Rows | Where-Object { -not (Test-Z4dBool $_.paper_usable) }).Count -gt 0) {
        return 'z4d_traceability_gap'
    }
    if (@($Rows | Where-Object {
                [int]$_.atom_cap -lt 256 -and $_.reference_kind -ne 'full_direct_tabi_energy'
            }).Count -gt 0) {
        return 'z4d_traceability_gap'
    }
    if (@($Rows | Where-Object {
                [int]$_.atom_cap -ge 256 -and $_.reference_kind -ne 'sampled_direct_operator_audit'
            }).Count -gt 0) {
        return 'z4d_traceability_gap'
    }
    foreach ($needle in @('grid remains the production default',
            'cap-256 direct TABI energy was not used as the energy oracle',
            'Sampled direct-operator audit')) {
        if ($MethodsText -notlike "*$needle*" -and $CaptionText -notlike "*$needle*") {
            return 'z4d_text_gap'
        }
    }
    return 'z4d_manuscript_assets_ready'
}

$mainRowsRaw = @(Import-Csv -LiteralPath $z4cMainCsv)
$supplementRowsRaw = @(Import-Csv -LiteralPath $z4cSupplementCsv)
Assert-Z4dColumns -Rows $mainRowsRaw -Name 'Z4C main table' -Columns @(
    'fixture','atom_cap','energy_kj_mol','reference_kind','direct_energy_delta_kj_mol',
    'direct_relative_energy_delta','sampled_audit_rms','wall_seconds','peak_working_set_mb',
    'paper_usable')
Assert-Z4dColumns -Rows $supplementRowsRaw -Name 'Z4C supplement table' -Columns @(
    'fixture','atom_cap','backend','table_role','status','energy_kj_mol',
    'wall_seconds','peak_working_set_mb')

$paperRows = @($mainRowsRaw | ForEach-Object {
        $rel = Convert-Z4dDouble $_.direct_relative_energy_delta
        [pscustomobject]@{
            fixture = $_.fixture
            atom_cap = [int]$_.atom_cap
            energy_kj_mol = $_.energy_kj_mol
            energy_kj_mol_display = Format-Z4dNumber $_.energy_kj_mol 'F3'
            reference_kind = $_.reference_kind
            reference_label = Convert-Z4dReferenceLabel $_.reference_kind
            direct_energy_delta_kj_mol = $_.direct_energy_delta_kj_mol
            direct_energy_delta_kj_mol_display = Format-Z4dNumber $_.direct_energy_delta_kj_mol 'F3'
            direct_relative_energy_delta_percent = if (Test-Z4dFinite $rel) {
                (100.0 * $rel).ToString('F4', [System.Globalization.CultureInfo]::InvariantCulture)
            } else {
                'NA'
            }
            sampled_audit_rms = $_.sampled_audit_rms
            sampled_audit_rms_display = Format-Z4dNumber $_.sampled_audit_rms 'G5'
            wall_seconds = $_.wall_seconds
            wall_seconds_display = Format-Z4dNumber $_.wall_seconds 'F1'
            peak_working_set_mb = $_.peak_working_set_mb
            peak_working_set_mb_display = Format-Z4dNumber $_.peak_working_set_mb 'F1'
            paper_usable = $_.paper_usable
        }
    })

$supplementRows = @($supplementRowsRaw | ForEach-Object {
        [pscustomobject]@{
            fixture = $_.fixture
            atom_cap = [int]$_.atom_cap
            backend = $_.backend
            table_role = $_.table_role
            status = $_.status
            energy_kj_mol = $_.energy_kj_mol
            energy_kj_mol_display = Format-Z4dNumber $_.energy_kj_mol 'F3'
            wall_seconds = $_.wall_seconds
            wall_seconds_display = Format-Z4dNumber $_.wall_seconds 'F1'
            peak_working_set_mb = $_.peak_working_set_mb
            peak_working_set_mb_display = Format-Z4dNumber $_.peak_working_set_mb 'F1'
            paper_oracle = $_.paper_oracle
        }
    })

$methodsText = @(
    '# Methods Box: TITAN TABI Treecode Protein Panel',
    '',
    'Protein-panel rows used the TABI treecode FMM/M2L operator with sparse near-field ILU(0), opening theta 0.6, Taylor order 4, analytic-singular near integration, and no Duffy quad buffer.',
    'For cap-128 rows, FMM/M2L reaction energies are compared against full direct TABI energy. For cap-256 rows, cap-256 direct TABI energy was not used as the energy oracle; instead, the paper reference is the sampled direct-operator audit from pbe_tabi_profile.csv.',
    'Grid rows are reported as runtime and memory context only. The grid remains the production default, and this manuscript package does not change TITAN defaults.'
)
$captionText = @(
    '# Figure Caption: Protein Panel Runtime And Verification',
    '',
    'Runtime and memory context for the Z4 protein panel. Blue FMM/M2L rows are paper-usable under the Z4B reference policy. Cap-128 rows use Direct TABI energy as the reference, while cap-256 rows use Sampled direct-operator audit. Grid rows are context only when they do not converge.'
)
$methodsFlat = ($methodsText -join "`n")
$captionFlat = ($captionText -join "`n")
$decision = Get-Z4dDecision -Rows $paperRows -MethodsText $methodsFlat -CaptionText $captionFlat

$paperRows | Export-Csv -LiteralPath $MainCsv -NoTypeInformation
Write-Z4dMarkdownTable -Rows $paperRows -Path $MainMd
$supplementRows | Export-Csv -LiteralPath $SupplementCsv -NoTypeInformation
Set-Content -LiteralPath $MethodsMd -Value $methodsText -Encoding ASCII
Set-Content -LiteralPath $CaptionMd -Value $captionText -Encoding ASCII

$maxAudit = ($paperRows | ForEach-Object { Convert-Z4dDouble $_.sampled_audit_rms } |
    Where-Object { Test-Z4dFinite $_ } | Measure-Object -Maximum).Maximum
if (-not (Test-Z4dFinite $maxAudit)) { $maxAudit = [double]::NaN }
[pscustomobject]@{
    decision = $decision
    main_rows = $paperRows.Count
    supplement_rows = $supplementRows.Count
    max_sampled_audit = $maxAudit
    methods_box = $MethodsMd
    figure_caption = $CaptionMd
    main_table_md = $MainMd
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'scope: Z4D manuscript methods and table export',
    "decision: $decision",
    "main_rows: $($paperRows.Count)",
    "supplement_rows: $($supplementRows.Count)",
    "max_sampled_audit: $maxAudit",
    'reference_policy: cap-128 FMM rows cite Direct TABI energy; cap-256 FMM rows cite Sampled direct-operator audit',
    'default_policy_note: grid remains production default; Z4D does not change public behavior'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z4D Manuscript Export Report',
    '',
    "Decision: $decision",
    '',
    'Z4D converts the Z4C paper package into manuscript-ready Markdown/CSV tables plus methods and caption text.',
    '',
    'Artifacts:',
    "- Main CSV: $MainCsv",
    "- Main Markdown: $MainMd",
    "- Supplement CSV: $SupplementCsv",
    "- Methods box: $MethodsMd",
    "- Figure caption: $CaptionMd",
    "- Decision: $DecisionCsv",
    "- Summary: $SummaryPath",
    '',
    'Main table preview:'
)
$report += Get-Content -LiteralPath $MainMd
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "PASS Scope Z4D manuscript export: $decision"
