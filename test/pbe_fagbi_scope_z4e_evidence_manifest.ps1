param(
    [switch]$RefreshZ4D,
    [string]$Z4DWorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4d'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4e')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$ArtifactManifestCsv = Join-Path $WorkRoot 'fagbi_scope_z4e_artifact_manifest.csv'
$ProvenanceCsv = Join-Path $WorkRoot 'fagbi_scope_z4e_command_provenance.csv'
$ClaimMapCsv = Join-Path $WorkRoot 'fagbi_scope_z4e_claim_map.csv'
$LimitationsMapCsv = Join-Path $WorkRoot 'fagbi_scope_z4e_limitations_map.csv'
$ChecksumsCsv = Join-Path $WorkRoot 'fagbi_scope_z4e_checksums.csv'
$BundleMd = Join-Path $WorkRoot 'fagbi_scope_z4e_evidence_bundle.md'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4e_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4e_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4E_EVIDENCE_MANIFEST_REPORT.txt'

if ($RefreshZ4D) {
    $z4dScript = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4d_manuscript_export.ps1'
    if (-not (Test-Path -LiteralPath $z4dScript)) {
        throw "Missing Z4D manuscript export script: $z4dScript"
    }
    & $z4dScript -WorkRoot $Z4DWorkRoot -RefreshZ4C
}

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Z4E work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z4eBool([object]$Value) {
    return ([string]$Value).Equals('true', [System.StringComparison]::OrdinalIgnoreCase)
}

function Get-Z4eReportPath([string]$ReportPath, [string]$Label) {
    if (-not (Test-Path -LiteralPath $ReportPath)) { return $null }
    $escaped = [regex]::Escape($Label)
    foreach ($line in Get-Content -LiteralPath $ReportPath) {
        if ($line -match ("^\s*$escaped\s*:\s*(.+)\s*$")) {
            return $Matches[1].Trim()
        }
    }
    return $null
}

function New-Z4eArtifact(
    [string]$Id,
    [string]$Scope,
    [string]$Kind,
    [string]$RelativePath,
    [string]$RerunCommand,
    [string]$FallbackPath = ''
) {
    $candidate = Join-Path $RepoRoot $RelativePath
    $resolved = $null
    $source = 'worktree'
    if (Test-Path -LiteralPath $candidate) {
        $resolved = [System.IO.Path]::GetFullPath($candidate)
    } elseif (-not [string]::IsNullOrWhiteSpace($FallbackPath) -and
        (Test-Path -LiteralPath $FallbackPath)) {
        $resolved = [System.IO.Path]::GetFullPath($FallbackPath)
        $source = 'report_pointer'
    }
    if ($null -eq $resolved) {
        return [pscustomobject]@{
            artifact_id = $Id
            scope = $Scope
            artifact_type = $Kind
            relative_path = $RelativePath
            resolved_path = ''
            source_location = 'missing'
            required = 'true'
            sha256 = ''
            size_bytes = ''
            last_write_utc = ''
            rerun_command = $RerunCommand
        }
    }

    $item = Get-Item -LiteralPath $resolved
    $hash = (Get-FileHash -LiteralPath $resolved -Algorithm SHA256).Hash.ToLowerInvariant()
    return [pscustomobject]@{
        artifact_id = $Id
        scope = $Scope
        artifact_type = $Kind
        relative_path = $RelativePath
        resolved_path = $resolved
        source_location = $source
        required = 'true'
        sha256 = $hash
        size_bytes = $item.Length
        last_write_utc = $item.LastWriteTimeUtc.ToString('o', [System.Globalization.CultureInfo]::InvariantCulture)
        rerun_command = $RerunCommand
    }
}

function Assert-Z4eDecision([string]$Path, [string]$ExpectedDecision, [string]$Name) {
    if (-not (Test-Path -LiteralPath $Path)) {
        throw "Missing decision artifact for ${Name}: $Path"
    }
    $rows = @(Import-Csv -LiteralPath $Path)
    if ($rows.Count -eq 0 -or @($rows[0].PSObject.Properties.Name) -notcontains 'decision') {
        throw "${Name} decision artifact has no decision column: $Path"
    }
    if ($rows[0].decision -ne $ExpectedDecision) {
        throw "${Name} decision is '$($rows[0].decision)', expected '$ExpectedDecision'."
    }
}

function Write-Z4eMarkdownTable([System.Collections.IEnumerable]$Rows, [string[]]$Columns) {
    $lines = New-Object System.Collections.Generic.List[string]
    $lines.Add('| ' + ($Columns -join ' | ') + ' |') | Out-Null
    $lines.Add('| ' + (($Columns | ForEach-Object { '---' }) -join ' | ') + ' |') | Out-Null
    foreach ($row in $Rows) {
        $values = foreach ($column in $Columns) {
            $value = [string]$row.$column
            $value = $value.Replace('|', '/')
            $value
        }
        $lines.Add('| ' + ($values -join ' | ') + ' |') | Out-Null
    }
    return $lines
}

function Add-Z4eLines([System.Collections.Generic.List[string]]$Target, [System.Collections.IEnumerable]$Lines) {
    foreach ($line in $Lines) {
        $Target.Add([string]$line) | Out-Null
    }
}

$z5Report = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z5_PRECONDITIONER_MONEY_PLOT_REPORT.txt'
$z5cReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z5C_MULTI_FIXTURE_MONEY_PLOT_REPORT.txt'
$z5Matrix = Get-Z4eReportPath -ReportPath $z5Report -Label 'Matrix'
$z5Trace = Get-Z4eReportPath -ReportPath $z5Report -Label 'Trace'
$z5Plot = Get-Z4eReportPath -ReportPath $z5Report -Label 'Plot'
$z5Dir = if ($z5Matrix) { Split-Path -Parent $z5Matrix } else { '' }
$z5Decision = if ($z5Dir) { Join-Path $z5Dir 'fmm_m2l_scope_z5_decision.csv' } else { '' }
$z5Summary = if ($z5Dir) { Join-Path $z5Dir 'fmm_m2l_scope_z5_summary.dat' } else { '' }
$z5cMatrix = Get-Z4eReportPath -ReportPath $z5cReport -Label 'Matrix'
$z5cTrace = Get-Z4eReportPath -ReportPath $z5cReport -Label 'Trace'
$z5cPlot = Get-Z4eReportPath -ReportPath $z5cReport -Label 'Plot'
$z5cDir = if ($z5cMatrix) { Split-Path -Parent $z5cMatrix } else { '' }
$z5cDecision = if ($z5cDir) { Join-Path $z5cDir 'fmm_m2l_scope_z5c_decision.csv' } else { '' }
$z5cSummary = if ($z5cDir) { Join-Path $z5cDir 'fmm_m2l_scope_z5c_summary.dat' } else { '' }

$artifacts = @(
    (New-Z4eArtifact 'z1d_highres_report' 'Z1D' 'root_report' 'TITAN_FMM_M2L_SCOPE_Z1D_HIGH_RESOLUTION_REPORT.txt' 'test\pbe_fagbi_scope_z1d_high_resolution_surface.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z1d_targeted_report' 'Z1D' 'root_report' 'TITAN_FMM_M2L_SCOPE_Z1D_TARGETED_VERIFY_REPORT.txt' 'test\pbe_fagbi_scope_z1d_targeted_missing_row.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z2a_matrix' 'Z2A' 'matrix_csv' 'test\tmp\fagbi_scope_z2a\fagbi_scope_z2a_matrix.csv' 'test\pbe_fagbi_scope_z2a_fmm_near_parity.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z2a_decision' 'Z2A' 'decision_csv' 'test\tmp\fagbi_scope_z2a\fagbi_scope_z2a_decision.csv' 'test\pbe_fagbi_scope_z2a_fmm_near_parity.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z2a_report' 'Z2A' 'root_report' 'TITAN_FMM_M2L_SCOPE_Z2A_FMM_NEAR_PARITY_REPORT.txt' 'test\pbe_fagbi_scope_z2a_fmm_near_parity.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z2b_matrix' 'Z2B' 'matrix_csv' 'test\tmp\fagbi_scope_z2b\fagbi_scope_z2b_matrix.csv' 'test\pbe_fagbi_scope_z2b_fmm_duffy_gate.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z2b_plot' 'Z2B' 'svg_plot' 'test\tmp\fagbi_scope_z2b\fagbi_scope_z2b_near_policy_plot.svg' 'test\pbe_fagbi_scope_z2b_fmm_duffy_gate.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z2b_report' 'Z2B' 'root_report' 'TITAN_FMM_M2L_SCOPE_Z2B_DUFFY_GATE_REPORT.txt' 'test\pbe_fagbi_scope_z2b_fmm_duffy_gate.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z3d_matrix' 'Z3D' 'matrix_csv' 'test\tmp\fagbi_scope_z3d\fagbi_scope_z3d_matrix.csv' 'test\pbe_fagbi_scope_z3d_theta_cap256_repeatability.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z3d_repeatability' 'Z3D' 'repeatability_csv' 'test\tmp\fagbi_scope_z3d\fagbi_scope_z3d_repeatability.csv' 'test\pbe_fagbi_scope_z3d_theta_cap256_repeatability.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z3d_plot' 'Z3D' 'svg_plot' 'test\tmp\fagbi_scope_z3d\fagbi_scope_z3d_theta_cap256_plot.svg' 'test\pbe_fagbi_scope_z3d_theta_cap256_repeatability.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z3d_report' 'Z3D' 'root_report' 'TITAN_FMM_M2L_SCOPE_Z3D_THETA_CAP256_REPEATABILITY_REPORT.txt' 'test\pbe_fagbi_scope_z3d_theta_cap256_repeatability.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z4b_matrix' 'Z4B' 'matrix_csv' 'test\tmp\fagbi_scope_z4b\fagbi_scope_z4b_matrix.csv' 'test\pbe_fagbi_scope_z4b_reference_policy_gate.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z4b_reference_policy' 'Z4B' 'reference_policy_csv' 'test\tmp\fagbi_scope_z4b\fagbi_scope_z4b_reference_policy.csv' 'test\pbe_fagbi_scope_z4b_reference_policy_gate.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z4b_report' 'Z4B' 'root_report' 'TITAN_FMM_M2L_SCOPE_Z4B_REFERENCE_POLICY_REPORT.txt' 'test\pbe_fagbi_scope_z4b_reference_policy_gate.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z4c_main_table' 'Z4C' 'paper_table_csv' 'test\tmp\fagbi_scope_z4c\fagbi_scope_z4c_main_table.csv' 'test\pbe_fagbi_scope_z4c_paper_packaging.ps1'),
    (New-Z4eArtifact 'z4c_supplement_table' 'Z4C' 'paper_table_csv' 'test\tmp\fagbi_scope_z4c\fagbi_scope_z4c_supplement_table.csv' 'test\pbe_fagbi_scope_z4c_paper_packaging.ps1'),
    (New-Z4eArtifact 'z4c_reference_table' 'Z4C' 'reference_policy_csv' 'test\tmp\fagbi_scope_z4c\fagbi_scope_z4c_reference_policy_table.csv' 'test\pbe_fagbi_scope_z4c_paper_packaging.ps1'),
    (New-Z4eArtifact 'z4c_plot' 'Z4C' 'svg_plot' 'test\tmp\fagbi_scope_z4c\fagbi_scope_z4c_runtime_memory_plot.svg' 'test\pbe_fagbi_scope_z4c_paper_packaging.ps1'),
    (New-Z4eArtifact 'z4c_decision' 'Z4C' 'decision_csv' 'test\tmp\fagbi_scope_z4c\fagbi_scope_z4c_decision.csv' 'test\pbe_fagbi_scope_z4c_paper_packaging.ps1'),
    (New-Z4eArtifact 'z4c_report' 'Z4C' 'root_report' 'TITAN_FMM_M2L_SCOPE_Z4C_PAPER_PACKAGE_REPORT.txt' 'test\pbe_fagbi_scope_z4c_paper_packaging.ps1'),
    (New-Z4eArtifact 'z4d_main_csv' 'Z4D' 'manuscript_csv' 'test\tmp\fagbi_scope_z4d\fagbi_scope_z4d_main_table_paper.csv' 'test\pbe_fagbi_scope_z4d_manuscript_export.ps1'),
    (New-Z4eArtifact 'z4d_main_md' 'Z4D' 'manuscript_markdown' 'test\tmp\fagbi_scope_z4d\fagbi_scope_z4d_main_table_paper.md' 'test\pbe_fagbi_scope_z4d_manuscript_export.ps1'),
    (New-Z4eArtifact 'z4d_methods' 'Z4D' 'manuscript_markdown' 'test\tmp\fagbi_scope_z4d\fagbi_scope_z4d_methods_box.md' 'test\pbe_fagbi_scope_z4d_manuscript_export.ps1'),
    (New-Z4eArtifact 'z4d_caption' 'Z4D' 'manuscript_markdown' 'test\tmp\fagbi_scope_z4d\fagbi_scope_z4d_figure_caption.md' 'test\pbe_fagbi_scope_z4d_manuscript_export.ps1'),
    (New-Z4eArtifact 'z4d_decision' 'Z4D' 'decision_csv' 'test\tmp\fagbi_scope_z4d\fagbi_scope_z4d_decision.csv' 'test\pbe_fagbi_scope_z4d_manuscript_export.ps1'),
    (New-Z4eArtifact 'z4d_report' 'Z4D' 'root_report' 'TITAN_FMM_M2L_SCOPE_Z4D_MANUSCRIPT_EXPORT_REPORT.txt' 'test\pbe_fagbi_scope_z4d_manuscript_export.ps1'),
    (New-Z4eArtifact 'z5_matrix' 'Z5B' 'matrix_csv' 'test\tmp\fmm_m2l_scope_z5\fmm_m2l_scope_z5_matrix.csv' 'test\pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1 -TitanExe .\titan.exe' $z5Matrix),
    (New-Z4eArtifact 'z5_trace' 'Z5B' 'residual_trace_csv' 'test\tmp\fmm_m2l_scope_z5\fmm_m2l_scope_z5_residual_trace.csv' 'test\pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1 -TitanExe .\titan.exe' $z5Trace),
    (New-Z4eArtifact 'z5_plot' 'Z5B' 'svg_plot' 'test\tmp\fmm_m2l_scope_z5\fmm_m2l_scope_z5_money_plot.svg' 'test\pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1 -TitanExe .\titan.exe' $z5Plot),
    (New-Z4eArtifact 'z5_decision' 'Z5B' 'decision_csv' 'test\tmp\fmm_m2l_scope_z5\fmm_m2l_scope_z5_decision.csv' 'test\pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1 -TitanExe .\titan.exe' $z5Decision),
    (New-Z4eArtifact 'z5_report' 'Z5B' 'root_report' 'TITAN_FMM_M2L_SCOPE_Z5_PRECONDITIONER_MONEY_PLOT_REPORT.txt' 'test\pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1 -TitanExe .\titan.exe'),
    (New-Z4eArtifact 'z5c_matrix' 'Z5C' 'matrix_csv' 'test\tmp\fmm_m2l_scope_z5c\fmm_m2l_scope_z5c_matrix.csv' 'test\pbe_fmm_m2l_scope_z5c_multi_fixture_money_plot.ps1 -TitanExe .\titan.exe' $z5cMatrix),
    (New-Z4eArtifact 'z5c_trace' 'Z5C' 'residual_trace_csv' 'test\tmp\fmm_m2l_scope_z5c\fmm_m2l_scope_z5c_residual_trace.csv' 'test\pbe_fmm_m2l_scope_z5c_multi_fixture_money_plot.ps1 -TitanExe .\titan.exe' $z5cTrace),
    (New-Z4eArtifact 'z5c_plot' 'Z5C' 'svg_plot' 'test\tmp\fmm_m2l_scope_z5c\fmm_m2l_scope_z5c_money_plot.svg' 'test\pbe_fmm_m2l_scope_z5c_multi_fixture_money_plot.ps1 -TitanExe .\titan.exe' $z5cPlot),
    (New-Z4eArtifact 'z5c_decision' 'Z5C' 'decision_csv' 'test\tmp\fmm_m2l_scope_z5c\fmm_m2l_scope_z5c_decision.csv' 'test\pbe_fmm_m2l_scope_z5c_multi_fixture_money_plot.ps1 -TitanExe .\titan.exe' $z5cDecision),
    (New-Z4eArtifact 'z5c_report' 'Z5C' 'root_report' 'TITAN_FMM_M2L_SCOPE_Z5C_MULTI_FIXTURE_MONEY_PLOT_REPORT.txt' 'test\pbe_fmm_m2l_scope_z5c_multi_fixture_money_plot.ps1 -TitanExe .\titan.exe')
)

$missing = @($artifacts | Where-Object { $_.source_location -eq 'missing' })
if ($missing.Count -gt 0) {
    $lines = New-Object System.Collections.Generic.List[string]
    $lines.Add('Missing required Z4E artifacts. Re-run the listed commands, then run Z4E again:') | Out-Null
    foreach ($row in $missing) {
        $lines.Add(('- {0}: {1} | rerun: {2}' -f $row.artifact_id, $row.relative_path, $row.rerun_command)) | Out-Null
    }
    throw ($lines -join "`n")
}

Assert-Z4eDecision -Path (Join-Path $RepoRoot 'test\tmp\fagbi_scope_z4c\fagbi_scope_z4c_decision.csv') -ExpectedDecision 'z4c_paper_package_ready' -Name 'Z4C'
Assert-Z4eDecision -Path (Join-Path $RepoRoot 'test\tmp\fagbi_scope_z4d\fagbi_scope_z4d_decision.csv') -ExpectedDecision 'z4d_manuscript_assets_ready' -Name 'Z4D'

$provenanceRows = @(
    [pscustomobject]@{ scope='Z1D'; script='test/pbe_fagbi_scope_z1d_high_resolution_surface.ps1'; expected_input_artifact='analytic Born/Kirkwood generated PQR'; output_family='high-resolution direct TABI reports'; run_type='solver-generated'; rerun_command='test\pbe_fagbi_scope_z1d_high_resolution_surface.ps1 -TitanExe .\titan.exe' },
    [pscustomobject]@{ scope='Z1D-targeted'; script='test/pbe_fagbi_scope_z1d_targeted_missing_row.ps1'; expected_input_artifact='missing off-center Kirkwood 20480 row'; output_family='targeted high-resolution verification'; run_type='solver-generated'; rerun_command='test\pbe_fagbi_scope_z1d_targeted_missing_row.ps1 -TitanExe .\titan.exe' },
    [pscustomobject]@{ scope='Z2A'; script='test/pbe_fagbi_scope_z2a_fmm_near_parity.ps1'; expected_input_artifact='Z1 near-field policy'; output_family='FMM analytic-singular near parity'; run_type='solver-generated'; rerun_command='test\pbe_fagbi_scope_z2a_fmm_near_parity.ps1 -TitanExe .\titan.exe' },
    [pscustomobject]@{ scope='Z2B'; script='test/pbe_fagbi_scope_z2b_fmm_duffy_gate.ps1'; expected_input_artifact='Z2A FMM near parity'; output_family='explicit Duffy near-field gate'; run_type='solver-generated'; rerun_command='test\pbe_fagbi_scope_z2b_fmm_duffy_gate.ps1 -TitanExe .\titan.exe' },
    [pscustomobject]@{ scope='Z3D'; script='test/pbe_fagbi_scope_z3d_theta_cap256_repeatability.ps1'; expected_input_artifact='Z3B/Z3C theta matrices'; output_family='theta 0.6 cap-256 repeatability'; run_type='solver-generated'; rerun_command='test\pbe_fagbi_scope_z3d_theta_cap256_repeatability.ps1 -TitanExe .\titan.exe' },
    [pscustomobject]@{ scope='Z4B'; script='test/pbe_fagbi_scope_z4b_reference_policy_gate.ps1'; expected_input_artifact='Z4A protein panel matrix'; output_family='reference policy table'; run_type='solver-generated'; rerun_command='test\pbe_fagbi_scope_z4b_reference_policy_gate.ps1 -TitanExe .\titan.exe' },
    [pscustomobject]@{ scope='Z4C'; script='test/pbe_fagbi_scope_z4c_paper_packaging.ps1'; expected_input_artifact='test/tmp/fagbi_scope_z4b'; output_family='paper tables and runtime/memory plot'; run_type='packaging-only'; rerun_command='test\pbe_fagbi_scope_z4c_paper_packaging.ps1' },
    [pscustomobject]@{ scope='Z4D'; script='test/pbe_fagbi_scope_z4d_manuscript_export.ps1'; expected_input_artifact='test/tmp/fagbi_scope_z4c'; output_family='manuscript tables, methods, caption'; run_type='packaging-only'; rerun_command='test\pbe_fagbi_scope_z4d_manuscript_export.ps1' },
    [pscustomobject]@{ scope='Z5B'; script='test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1'; expected_input_artifact='Barnase cap-256 sampled fixture'; output_family='preconditioner money plot'; run_type='solver-generated'; rerun_command='test\pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1 -TitanExe .\titan.exe' },
    [pscustomobject]@{ scope='Z5C'; script='test/pbe_fmm_m2l_scope_z5c_multi_fixture_money_plot.ps1'; expected_input_artifact='Z5B money plot harness'; output_family='multi-fixture repeatability/generalization'; run_type='solver-generated'; rerun_command='test\pbe_fmm_m2l_scope_z5c_multi_fixture_money_plot.ps1 -TitanExe .\titan.exe' }
)

$claimRows = @(
    [pscustomobject]@{ claim_id='high_resolution_analytic_convergence'; claim='Direct TABI analytic Born/Kirkwood high-resolution error falls below the Z1D threshold on the verified rows.'; evidence_artifacts='z1d_highres_report;z1d_targeted_report'; paper_use='analytic validation'; caveat='Direct high-resolution rows are expensive and remain benchmark-only.' },
    [pscustomobject]@{ claim_id='fmm_analytic_singular_near_parity'; claim='FMM/M2L can use analytic-singular near integration with finite parity/audit behavior.'; evidence_artifacts='z2a_matrix;z2a_report'; paper_use='near-field implementation validation'; caveat='Duffy is measured separately and remains explicit-only.' },
    [pscustomobject]@{ claim_id='duffy_explicit_only'; claim='Duffy-adaptive FMM near-field rows are safe enough to document as explicit, but not promoted to near=auto.'; evidence_artifacts='z2b_matrix;z2b_plot;z2b_report'; paper_use='near-field policy limitation'; caveat='No automatic Duffy promotion follows from Z2B.' },
    [pscustomobject]@{ claim_id='theta06_order4_policy'; claim='The paper protein panel uses fixed Taylor order 4 and opening theta 0.6.'; evidence_artifacts='z3d_matrix;z3d_repeatability;z3d_report'; paper_use='treecode policy'; caveat='Adaptive order remains future work.' },
    [pscustomobject]@{ claim_id='protein_panel_reference_policy'; claim='Cap-128 protein rows use full direct TABI energy; cap-256 rows use sampled direct-operator audit.'; evidence_artifacts='z4b_reference_policy;z4c_main_table;z4d_main_md;z4d_methods'; paper_use='main protein table'; caveat='Grid rows are runtime/memory context only when not converged.' },
    [pscustomobject]@{ claim_id='preconditioner_money_plot'; claim='Sparse near-field ILU(0) breaks through the Barnase dense/high-coupling cap-256 case where block-diagonal stalls.'; evidence_artifacts='z5_matrix;z5_trace;z5_plot;z5_report;z5c_report'; paper_use='TITAN novelty figure'; caveat='Z5C narrows the generalization claim to Barnase-style dense/high-coupling meshes.' }
)

$limitationRows = @(
    [pscustomobject]@{ limitation_id='grid_default'; statement='Grid remains the production PBE default.'; evidence_artifacts='z4d_methods;z4d_report'; manuscript_note='Do not state that FMM/M2L replaces grid as default.' },
    [pscustomobject]@{ limitation_id='cap256_direct_energy_absent'; statement='Cap-256 direct TABI energy is not used as the paper oracle.'; evidence_artifacts='z4b_reference_policy;z4d_methods'; manuscript_note='Use sampled direct-operator audit for cap-256 verification.' },
    [pscustomobject]@{ limitation_id='z5_barnase_framing'; statement='The Z5 money plot is Barnase dense/high-coupling evidence, not broad protein-generalization evidence.'; evidence_artifacts='z5_report;z5c_report'; manuscript_note='Frame as targeted breakthrough first.' },
    [pscustomobject]@{ limitation_id='duffy_explicit_only'; statement='Duffy-adaptive is explicit-only and not promoted to near=auto.'; evidence_artifacts='z2b_report'; manuscript_note='Do not imply Duffy is the default FMM near policy.' },
    [pscustomobject]@{ limitation_id='z4e_no_new_numerics'; statement='Z4E is packaging and audit only, not new numerical evidence.'; evidence_artifacts='z4e_artifact_manifest'; manuscript_note='Cite upstream Z1-Z5 artifacts for numerical claims.' }
)

$decision = 'z4e_evidence_bundle_ready'
if ($claimRows.Count -lt 6 -or $limitationRows.Count -lt 5) {
    $decision = 'z4e_traceability_gap'
}

$artifacts | Export-Csv -LiteralPath $ArtifactManifestCsv -NoTypeInformation
$provenanceRows | Export-Csv -LiteralPath $ProvenanceCsv -NoTypeInformation
$claimRows | Export-Csv -LiteralPath $ClaimMapCsv -NoTypeInformation
$limitationRows | Export-Csv -LiteralPath $LimitationsMapCsv -NoTypeInformation
$artifacts | Select-Object artifact_id, scope, artifact_type, resolved_path, sha256, size_bytes, last_write_utc |
    Export-Csv -LiteralPath $ChecksumsCsv -NoTypeInformation

$bundle = New-Object System.Collections.Generic.List[string]
$bundle.Add('# TITAN FAGBI Paper Evidence Bundle') | Out-Null
$bundle.Add('') | Out-Null
$bundle.Add("Decision: $decision") | Out-Null
$bundle.Add('') | Out-Null
$bundle.Add('This bundle maps paper-facing claims to the measured artifacts generated across Scope Z1-Z5. It is packaging and audit only; it does not rerun solver rows or change TITAN defaults.') | Out-Null
$bundle.Add('') | Out-Null
$bundle.Add('## Claim Map') | Out-Null
Add-Z4eLines -Target $bundle -Lines (Write-Z4eMarkdownTable -Rows $claimRows -Columns @('claim_id','paper_use','evidence_artifacts','caveat'))
$bundle.Add('') | Out-Null
$bundle.Add('## Limitations') | Out-Null
Add-Z4eLines -Target $bundle -Lines (Write-Z4eMarkdownTable -Rows $limitationRows -Columns @('limitation_id','statement','manuscript_note'))
$bundle.Add('') | Out-Null
$bundle.Add('## Command Provenance') | Out-Null
Add-Z4eLines -Target $bundle -Lines (Write-Z4eMarkdownTable -Rows $provenanceRows -Columns @('scope','script','run_type','output_family'))
$bundle.Add('') | Out-Null
$bundle.Add('## Required Artifact Checksums') | Out-Null
Add-Z4eLines -Target $bundle -Lines (Write-Z4eMarkdownTable -Rows $artifacts -Columns @('artifact_id','scope','artifact_type','source_location','sha256'))
Set-Content -LiteralPath $BundleMd -Value $bundle -Encoding ASCII

[pscustomobject]@{
    decision = $decision
    artifact_count = $artifacts.Count
    claim_count = $claimRows.Count
    limitation_count = $limitationRows.Count
    report = $ReportPath
    bundle = $BundleMd
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'scope: Z4E paper evidence manifest and claim traceability',
    "decision: $decision",
    "artifact_count: $($artifacts.Count)",
    "claim_count: $($claimRows.Count)",
    "limitation_count: $($limitationRows.Count)",
    'default_policy_note: grid remains production default; Z4E does not change public behavior',
    'numerical_evidence_note: Z4E is packaging/audit only and introduces no new solver results'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = New-Object System.Collections.Generic.List[string]
$report.Add('# TITAN FMM/M2L Scope Z4E Evidence Manifest Report') | Out-Null
$report.Add('') | Out-Null
$report.Add("Decision: $decision") | Out-Null
$report.Add('') | Out-Null
$report.Add('Artifacts:') | Out-Null
$report.Add("- Artifact manifest: $ArtifactManifestCsv") | Out-Null
$report.Add("- Command provenance: $ProvenanceCsv") | Out-Null
$report.Add("- Claim map: $ClaimMapCsv") | Out-Null
$report.Add("- Limitations map: $LimitationsMapCsv") | Out-Null
$report.Add("- Checksums: $ChecksumsCsv") | Out-Null
$report.Add("- Evidence bundle: $BundleMd") | Out-Null
$report.Add("- Decision: $DecisionCsv") | Out-Null
$report.Add("- Summary: $SummaryPath") | Out-Null
$report.Add('') | Out-Null
$report.Add('Claims:') | Out-Null
foreach ($claim in $claimRows) {
    $report.Add(('- {0}: {1}' -f $claim.claim_id, $claim.claim)) | Out-Null
}
$report.Add('') | Out-Null
$report.Add('Limitations:') | Out-Null
foreach ($limitation in $limitationRows) {
    $report.Add(('- {0}: {1}' -f $limitation.limitation_id, $limitation.statement)) | Out-Null
}
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "PASS Scope Z4E evidence manifest: $decision"
