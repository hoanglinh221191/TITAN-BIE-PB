param(
    [switch]$RefreshZ4E,
    [switch]$SkipVisualQa,
    [switch]$VisualQaPassed,
    [string]$PythonExe = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [string]$EdgeExe = 'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe',
    [string]$Z4DWorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4d'),
    [string]$Z4EWorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4e'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4f')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Builder = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4f_build_docx.py'
$DocxPath = Join-Path $WorkRoot 'titan_fagbi_manuscript_skeleton.docx'
$SourcePath = Join-Path $WorkRoot 'fagbi_scope_z4f_manuscript_source.md'
$ClaimTracePath = Join-Path $WorkRoot 'fagbi_scope_z4f_claim_trace.csv'
$CitationGapsPath = Join-Path $WorkRoot 'fagbi_scope_z4f_citation_gaps.csv'
$QaPath = Join-Path $WorkRoot 'fagbi_scope_z4f_docx_qa.csv'
$DecisionPath = Join-Path $WorkRoot 'fagbi_scope_z4f_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z4f_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4F_MANUSCRIPT_SKELETON_REPORT.txt'
$RenderRoot = Join-Path $WorkRoot 'render_qa'
$PdfPath = Join-Path $RenderRoot 'titan_fagbi_manuscript_skeleton.pdf'

if ($SkipVisualQa -and $VisualQaPassed) {
    throw '-SkipVisualQa and -VisualQaPassed cannot be used together.'
}
if (-not (Test-Path -LiteralPath $PythonExe)) { throw "Missing bundled Python: $PythonExe" }
if (-not (Test-Path -LiteralPath $Builder)) { throw "Missing Z4F DOCX builder: $Builder" }

if ($RefreshZ4E) {
    $z4eScript = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4e_evidence_manifest.ps1'
    if (-not (Test-Path -LiteralPath $z4eScript)) { throw "Missing Z4E script: $z4eScript" }
    & $z4eScript -RefreshZ4D -Z4DWorkRoot $Z4DWorkRoot -WorkRoot $Z4EWorkRoot
}

$required = @{
    z4d_main = Join-Path $Z4DWorkRoot 'fagbi_scope_z4d_main_table_paper.csv'
    z4d_supplement = Join-Path $Z4DWorkRoot 'fagbi_scope_z4d_supplement_table_paper.csv'
    z4d_methods = Join-Path $Z4DWorkRoot 'fagbi_scope_z4d_methods_box.md'
    z4d_caption = Join-Path $Z4DWorkRoot 'fagbi_scope_z4d_figure_caption.md'
    z4d_decision = Join-Path $Z4DWorkRoot 'fagbi_scope_z4d_decision.csv'
    z4e_manifest = Join-Path $Z4EWorkRoot 'fagbi_scope_z4e_artifact_manifest.csv'
    z4e_checksums = Join-Path $Z4EWorkRoot 'fagbi_scope_z4e_checksums.csv'
    z4e_claims = Join-Path $Z4EWorkRoot 'fagbi_scope_z4e_claim_map.csv'
    z4e_limitations = Join-Path $Z4EWorkRoot 'fagbi_scope_z4e_limitations_map.csv'
    z4e_decision = Join-Path $Z4EWorkRoot 'fagbi_scope_z4e_decision.csv'
}
$missing = @($required.GetEnumerator() | Where-Object { -not (Test-Path -LiteralPath $_.Value) })
if ($missing.Count -gt 0) {
    $details = $missing | ForEach-Object { "$($_.Key)=$($_.Value)" }
    throw "z4f_missing_artifacts: $($details -join '; '). Run test\pbe_fagbi_scope_z4e_evidence_manifest.ps1 -RefreshZ4D."
}

$z4dDecision = @(Import-Csv -LiteralPath $required.z4d_decision)
$z4eDecision = @(Import-Csv -LiteralPath $required.z4e_decision)
if ($z4dDecision.Count -ne 1 -or $z4dDecision[0].decision -ne 'z4d_manuscript_assets_ready') {
    throw "z4f_evidence_gap: Z4D decision is not z4d_manuscript_assets_ready."
}
if ($z4eDecision.Count -ne 1 -or $z4eDecision[0].decision -ne 'z4e_evidence_bundle_ready') {
    throw "z4f_evidence_gap: Z4E decision is not z4e_evidence_bundle_ready."
}

$checksums = @(Import-Csv -LiteralPath $required.z4e_checksums)
if ($checksums.Count -eq 0) { throw 'z4f_evidence_gap: empty Z4E checksum manifest.' }
$drift = New-Object System.Collections.Generic.List[string]
foreach ($row in $checksums) {
    if (-not (Test-Path -LiteralPath $row.resolved_path)) {
        $drift.Add("missing $($row.artifact_id): $($row.resolved_path)") | Out-Null
        continue
    }
    $actual = (Get-FileHash -LiteralPath $row.resolved_path -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($actual -ne ([string]$row.sha256).ToLowerInvariant()) {
        $drift.Add("checksum drift $($row.artifact_id): expected $($row.sha256), got $actual") | Out-Null
    }
}
if ($drift.Count -gt 0) {
    throw ("z4f_checksum_drift: " + ($drift -join '; '))
}

$manifest = @(Import-Csv -LiteralPath $required.z4e_manifest)
function Get-Z4fArtifactPath([string]$ArtifactId) {
    $rows = @($manifest | Where-Object { $_.artifact_id -eq $ArtifactId })
    if ($rows.Count -ne 1 -or -not (Test-Path -LiteralPath $rows[0].resolved_path)) {
        throw "z4f_missing_artifacts: unresolved artifact '$ArtifactId'."
    }
    return [string]$rows[0].resolved_path
}
$runtimeSvg = Get-Z4fArtifactPath 'z4c_plot'
$moneySvg = Get-Z4fArtifactPath 'z5_plot'

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Z4F work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Convert-Z4fSvgToPng([string]$SvgPath, [string]$PngPath, [string]$ProfileName) {
    if (-not (Test-Path -LiteralPath $EdgeExe)) { throw "Missing Microsoft Edge: $EdgeExe" }
    $profile = Join-Path $WorkRoot "edge_$ProfileName"
    New-Item -ItemType Directory -Path $profile -Force | Out-Null
    $uri = ([System.Uri]([System.IO.Path]::GetFullPath($SvgPath))).AbsoluteUri
    $svgHead = Get-Content -LiteralPath $SvgPath -TotalCount 1
    if ($svgHead -notmatch 'width="(\d+)"\s+height="(\d+)"') {
        throw "Cannot resolve SVG dimensions: $SvgPath"
    }
    $windowSize = "$($Matches[1]),$($Matches[2])"
    $edgeArgs = @(
        '--headless=new', '--disable-gpu', '--disable-sync', '--hide-scrollbars', '--no-first-run',
        "--user-data-dir=$profile", "--window-size=$windowSize",
        "--screenshot=$PngPath", $uri
    )
    & $EdgeExe @edgeArgs | Out-Null
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $PngPath)) {
        throw "Failed to rasterize SVG: $SvgPath"
    }
}

$runtimePng = Join-Path $WorkRoot 'figure1_runtime_memory.png'
$moneyPng = Join-Path $WorkRoot 'figure2_convergence.png'
Convert-Z4fSvgToPng -SvgPath $runtimeSvg -PngPath $runtimePng -ProfileName 'runtime'
Convert-Z4fSvgToPng -SvgPath $moneySvg -PngPath $moneyPng -ProfileName 'money'

$builderArgs = @(
    $Builder,
    '--main-csv', $required.z4d_main,
    '--supplement-csv', $required.z4d_supplement,
    '--claim-map', $required.z4e_claims,
    '--limitations-map', $required.z4e_limitations,
    '--runtime-figure', $runtimePng,
    '--money-figure', $moneyPng,
    '--output-docx', $DocxPath,
    '--output-source', $SourcePath,
    '--output-claim-trace', $ClaimTracePath,
    '--output-citation-gaps', $CitationGapsPath,
    '--output-qa', $QaPath
)
& $PythonExe @builderArgs
if ($LASTEXITCODE -ne 0) { throw "Z4F DOCX builder failed with exit code $LASTEXITCODE" }

function Export-Z4fDocxToPdf([string]$InputDocx, [string]$OutputPdf) {
    New-Item -ItemType Directory -Path (Split-Path -Parent $OutputPdf) -Force | Out-Null
    $word = $null
    $doc = $null
    try {
        $word = New-Object -ComObject Word.Application
        $word.Visible = $false
        $word.DisplayAlerts = 0
        $doc = $word.Documents.Open([System.IO.Path]::GetFullPath($InputDocx), $false, $true)
        $doc.ExportAsFixedFormat([System.IO.Path]::GetFullPath($OutputPdf), 17)
    } finally {
        if ($null -ne $doc) { $doc.Close($false); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($doc) }
        if ($null -ne $word) { $word.Quit(); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($word) }
        [GC]::Collect()
        [GC]::WaitForPendingFinalizers()
    }
    if (-not (Test-Path -LiteralPath $OutputPdf)) { throw "Word did not export PDF: $OutputPdf" }
}

function Get-Z4fPdfPageCount([string]$Path) {
    $code = 'from pypdf import PdfReader; import sys; print(len(PdfReader(sys.argv[1]).pages))'
    $value = (& $PythonExe -c $code $Path | Select-Object -Last 1).Trim()
    return [int]$value
}

function Capture-Z4fPdfPages([string]$Path, [int]$PageCount) {
    $renderScript = Join-Path $RenderRoot 'render_pdf_pages.py'
    $renderCode = @'
from pathlib import Path
import sys
import pypdfium2 as pdfium

pdf_path = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
output_dir.mkdir(parents=True, exist_ok=True)
pdf = pdfium.PdfDocument(str(pdf_path))
for index in range(len(pdf)):
    page = pdf[index]
    bitmap = page.render(scale=2.0)
    image = bitmap.to_pil()
    image.save(output_dir / f"page-{index + 1:03d}.png")
    page.close()
'@
    Set-Content -LiteralPath $renderScript -Value $renderCode -Encoding ASCII
    & $PythonExe $renderScript $Path $RenderRoot
    if ($LASTEXITCODE -ne 0) { throw 'pypdfium2 failed to render the manuscript PDF.' }
    $captures = @(Get-ChildItem -LiteralPath $RenderRoot -Filter 'page-*.png' | Sort-Object Name)
    if ($captures.Count -ne $PageCount) {
        throw "Expected $PageCount rendered pages, found $($captures.Count)."
    }
    return @($captures.FullName)
}

$pageCount = 0
$renderStatus = 'skipped'
if (-not $SkipVisualQa) {
    Export-Z4fDocxToPdf -InputDocx $DocxPath -OutputPdf $PdfPath
    $pageCount = Get-Z4fPdfPageCount $PdfPath
    if ($pageCount -lt 1) { throw 'z4f_render_gap: exported PDF has no pages.' }
    $captures = @(Capture-Z4fPdfPages -Path $PdfPath -PageCount $pageCount)
    if ($captures.Count -ne $pageCount) { throw 'z4f_render_gap: incomplete page captures.' }
    $renderStatus = if ($VisualQaPassed) { 'pass' } else { 'pending_visual_review' }
}

$qa = @(Import-Csv -LiteralPath $QaPath)
$qa += [pscustomobject]@{ check='checksum_verification'; status='pass'; value=$checksums.Count; detail='all Z4E SHA256 values match current artifacts' }
$qa += [pscustomobject]@{ check='render_page_count'; status=if($SkipVisualQa){'skipped'}else{'pass'}; value=$pageCount; detail='Word COM PDF export and pypdfium2 page render' }
$qa += [pscustomobject]@{ check='visual_review'; status=$renderStatus; value=$renderStatus; detail='all pages require 100% visual inspection' }
$qa | Export-Csv -LiteralPath $QaPath -NoTypeInformation

$claimCount = @(Import-Csv -LiteralPath $required.z4e_claims).Count
$limitationCount = @(Import-Csv -LiteralPath $required.z4e_limitations).Count
$proteinRowCount = @(Import-Csv -LiteralPath $required.z4d_main).Count
$structuralFailures = @($qa | Where-Object { $_.status -eq 'fail' })
$decision = if ($structuralFailures.Count -gt 0 -or $claimCount -ne 6 -or $limitationCount -ne 5 -or $proteinRowCount -ne 6) {
    'z4f_evidence_gap'
} elseif ($renderStatus -ne 'pass') {
    'z4f_render_gap'
} else {
    'z4f_manuscript_skeleton_ready'
}

[pscustomobject]@{
    decision = $decision
    claim_count = $claimCount
    limitation_count = $limitationCount
    protein_row_count = $proteinRowCount
    figure_count = 2
    render_page_count = $pageCount
    visual_review = $renderStatus
    manuscript_docx = $DocxPath
} | Export-Csv -LiteralPath $DecisionPath -NoTypeInformation

$summary = @(
    'scope: Z4F evidence-backed Word manuscript skeleton',
    "decision: $decision",
    "claim_count: $claimCount",
    "limitation_count: $limitationCount",
    "protein_row_count: $proteinRowCount",
    'figure_count: 2',
    "render_page_count: $pageCount",
    "visual_review: $renderStatus",
    'reference_policy: cap-128 direct TABI energy; cap-256 sampled direct-operator audit',
    'default_policy_note: grid remains production default; Z4F adds no numerical evidence'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Z4F Manuscript Skeleton Report',
    '',
    "Decision: $decision",
    '',
    'Z4F converts the validated Z4D manuscript assets and Z4E evidence manifest into a journal-neutral Word manuscript skeleton.',
    '',
    'Artifacts:',
    "- DOCX: $DocxPath",
    "- Markdown source: $SourcePath",
    "- Claim trace: $ClaimTracePath",
    "- Citation gaps: $CitationGapsPath",
    "- DOCX QA: $QaPath",
    "- PDF render: $PdfPath",
    "- Decision: $DecisionPath",
    "- Summary: $SummaryPath",
    '',
    'Policy:',
    '- Grid remains the production PBE default.',
    '- Cap-128 rows cite full direct TABI energy.',
    '- Cap-256 rows cite sampled direct-operator audit.',
    '- Duffy remains explicit-only; the paper policy uses theta 0.6 and Taylor order 4.',
    '- The Barnase money plot is framed as dense/high-coupling evidence, not broad generalization.',
    '',
    'Citation status:',
    '- Bibliographic gaps remain explicit in the separate citation-gaps CSV; no references were invented.'
)
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "PASS Scope Z4F decision=$decision docx=$DocxPath"
