param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$Z4oWorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4o'),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4p'),
    [string]$PythonExe = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [switch]$SkipRender
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$builder = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4p_build_docx.py'
$matrixCsv = Join-Path $Z4oWorkRoot 'fagbi_scope_z4o_matrix.csv'
$repeatCsv = Join-Path $Z4oWorkRoot 'fagbi_scope_z4o_repeatability.csv'
$z4oDecisionCsv = Join-Path $Z4oWorkRoot 'fagbi_scope_z4o_decision.csv'
$snapshotMatrixCsv = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4p_z4o_snapshot_matrix.csv'
$snapshotRepeatCsv = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4p_z4o_snapshot_repeatability.csv'
$snapshotDecisionCsv = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4p_z4o_snapshot_decision.csv'
$z4oPlot = Join-Path $Z4oWorkRoot 'fagbi_scope_z4o_repeatability_plot.png'
$z4oSvg = Join-Path $Z4oWorkRoot 'fagbi_scope_z4o_repeatability_plot.svg'
$manifestCsv = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$docxPath = Join-Path $WorkRoot 'titan_fagbi_guarded_p4_manuscript.docx'
$sourcePath = Join-Path $WorkRoot 'fagbi_scope_z4p_manuscript_source.md'
$claimTraceCsv = Join-Path $WorkRoot 'fagbi_scope_z4p_claim_trace.csv'
$qaCsv = Join-Path $WorkRoot 'fagbi_scope_z4p_docx_qa.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4p_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4p_summary.dat'
$renderRoot = Join-Path $WorkRoot 'render_qa'
$pdfPath = Join-Path $renderRoot 'titan_fagbi_guarded_p4_manuscript.pdf'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4P_MANUSCRIPT_GUARDED_P4_REPORT.txt'

foreach ($path in @($PythonExe,$builder,$manifestCsv,$snapshotMatrixCsv,$snapshotRepeatCsv,$snapshotDecisionCsv)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "z4p_evidence_gap: missing required artifact $path"
    }
}

function Test-Z4pFullZ4oArtifact([string]$MatrixPath,[string]$RepeatPath,[string]$DecisionPath) {
    if (-not ((Test-Path -LiteralPath $MatrixPath) -and
            (Test-Path -LiteralPath $RepeatPath) -and
            (Test-Path -LiteralPath $DecisionPath))) {
        return $false
    }
    $matrixRows = @(Import-Csv -LiteralPath $MatrixPath)
    $repeatRows = @(Import-Csv -LiteralPath $RepeatPath)
    $decisionRows = @(Import-Csv -LiteralPath $DecisionPath)
    $fixtures = @($matrixRows | Where-Object { $_.variant -eq 'p4_guarded' } |
        Select-Object -ExpandProperty fixture -Unique)
    return ($fixtures.Count -ge 3 -and $repeatRows.Count -ge 3 -and
        $decisionRows.Count -eq 1 -and $decisionRows[0].decision -eq 'z4o_repeatability_pass')
}

$evidenceSource = 'z4o_current_artifacts'
if (-not (Test-Z4pFullZ4oArtifact $matrixCsv $repeatCsv $z4oDecisionCsv)) {
    $matrixCsv = $snapshotMatrixCsv
    $repeatCsv = $snapshotRepeatCsv
    $z4oDecisionCsv = $snapshotDecisionCsv
    $z4oPlot = ''
    $z4oSvg = ''
    $evidenceSource = 'z4o_full_snapshot'
}

$workFull = [IO.Path]::GetFullPath($WorkRoot)
$tmpRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workFull.StartsWith($tmpRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Z4P work root outside test/tmp: $workFull"
}
if (Test-Path -LiteralPath $workFull) { Remove-Item -LiteralPath $workFull -Recurse -Force }
New-Item -ItemType Directory -Path $workFull -Force | Out-Null
New-Item -ItemType Directory -Path $renderRoot -Force | Out-Null

if ((-not [string]::IsNullOrWhiteSpace($z4oSvg)) -and
    (-not [string]::IsNullOrWhiteSpace($z4oPlot)) -and
    (Test-Path -LiteralPath $z4oSvg) -and
    -not (Test-Path -LiteralPath $z4oPlot)) {
    $edge = 'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe'
    if (Test-Path -LiteralPath $edge) {
        $profile = Join-Path $WorkRoot 'edge_z4o_plot'
        New-Item -ItemType Directory -Path $profile -Force | Out-Null
        $head = Get-Content -LiteralPath $z4oSvg -TotalCount 1
        $window = if ($head -match 'width="(\d+)"\s+height="(\d+)"') { "$($Matches[1]),$($Matches[2])" } else { '900,420' }
        $uri = ([Uri]([IO.Path]::GetFullPath($z4oSvg))).AbsoluteUri
        & $edge '--headless=new' '--disable-gpu' '--disable-sync' '--hide-scrollbars' '--no-first-run' `
            "--user-data-dir=$profile" "--window-size=$window" "--screenshot=$z4oPlot" $uri | Out-Null
    }
}
$plotArg = if ((-not [string]::IsNullOrWhiteSpace($z4oPlot)) -and
    (Test-Path -LiteralPath $z4oPlot)) { $z4oPlot } else { '' }

$builderArgs = @(
    $builder,
    '--matrix', $matrixCsv,
    '--repeatability', $repeatCsv,
    '--decision', $z4oDecisionCsv,
    '--provenance', $manifestCsv,
    '--output-docx', $docxPath,
    '--output-source', $sourcePath,
    '--claim-trace', $claimTraceCsv,
    '--qa', $qaCsv
)
if (-not [string]::IsNullOrWhiteSpace($plotArg)) {
    $builderArgs += @('--repeatability-plot', $plotArg)
}
& $PythonExe @builderArgs
if ($LASTEXITCODE -ne 0) { throw "Z4P DOCX builder failed with exit code $LASTEXITCODE" }

function Export-Z4pDocxToPdf([string]$InputDocx, [string]$OutputPdf) {
    New-Item -ItemType Directory -Path (Split-Path -Parent $OutputPdf) -Force | Out-Null
    $word = $null
    $doc = $null
    try {
        $word = New-Object -ComObject Word.Application
        $word.Visible = $false
        $word.DisplayAlerts = 0
        $doc = $word.Documents.Open([IO.Path]::GetFullPath($InputDocx), $false, $true)
        $doc.ExportAsFixedFormat([IO.Path]::GetFullPath($OutputPdf), 17)
    } finally {
        if ($null -ne $doc) { $doc.Close($false); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($doc) }
        if ($null -ne $word) { $word.Quit(); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($word) }
        [GC]::Collect(); [GC]::WaitForPendingFinalizers()
    }
    if (-not (Test-Path -LiteralPath $OutputPdf)) { throw "Word did not export PDF: $OutputPdf" }
}

function Render-Z4pPdfPages([string]$InputPdf) {
    $script = Join-Path $renderRoot 'render_pdf_pages.py'
    $code = @'
from pathlib import Path
import sys
import pypdfium2 as pdfium
pdf = pdfium.PdfDocument(str(Path(sys.argv[1])))
out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
for index in range(len(pdf)):
    page = pdf[index]
    image = page.render(scale=2.0).to_pil()
    image.save(out / f"page-{index + 1:03d}.png")
    page.close()
print(len(pdf))
'@
    Set-Content -LiteralPath $script -Value $code -Encoding ASCII
    $count = (& $PythonExe $script $InputPdf $renderRoot | Select-Object -Last 1).Trim()
    if ($LASTEXITCODE -ne 0) { throw 'z4p_render_gap: PDF page rendering failed.' }
    $captures = @(Get-ChildItem -LiteralPath $renderRoot -Filter 'page-*.png' | Sort-Object Name)
    if ($captures.Count -ne [int]$count -or $captures.Count -lt 1) { throw 'z4p_render_gap: incomplete page captures.' }
    return $captures.Count
}

$renderStatus = 'skipped'
$pageCount = 0
if (-not $SkipRender) {
    try {
        Export-Z4pDocxToPdf $docxPath $pdfPath
        $pageCount = Render-Z4pPdfPages $pdfPath
        $renderStatus = 'pass'
    } catch {
        $renderStatus = 'z4p_render_gap'
        Set-Content -LiteralPath (Join-Path $WorkRoot 'render_error.txt') -Value $_.Exception.Message -Encoding UTF8
    }
}

$qaRows = @(Import-Csv -LiteralPath $qaCsv)
$sourceText = Get-Content -Raw -LiteralPath $sourcePath
$forbidden = @('cap-128','cap-256','research-medium','research-alt','money plot','Z1-Z5','Z4 reference policy','audit-clean')
$leaks = @($forbidden | Where-Object { $sourceText -match [Regex]::Escape($_) })
$qaRows += [pscustomobject]@{ check='render_pages'; status=if($renderStatus -eq 'pass' -or $SkipRender){'pass'}else{'fail'}; value=$pageCount }
$qaRows += [pscustomobject]@{ check='forbidden_terms'; status=if($leaks.Count -eq 0){'pass'}else{'fail'}; value=($leaks -join ';') }
$qaRows | Export-Csv -LiteralPath $qaCsv -NoTypeInformation

$decisionValue = if (@($qaRows | Where-Object status -eq 'fail').Count -gt 0) {
    if ($renderStatus -eq 'z4p_render_gap') { 'z4p_render_gap' } else { 'z4p_evidence_gap' }
} else {
    'z4p_manuscript_ready'
}

[pscustomobject]@{
    decision=$decisionValue; docx=$docxPath; pages=$pageCount; render_status=$renderStatus;
    evidence_source=$evidenceSource; source=$sourcePath; claim_trace=$claimTraceCsv; qa=$qaCsv
} | Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
@(
    'scope: Z4P guarded p4 manuscript refresh',
    "decision: $decisionValue",
    "evidence_source: $evidenceSource",
    "docx: $docxPath",
    "pages: $pageCount",
    "render_status: $renderStatus",
    'policy: guarded p4 is advanced explicit option only',
    'default_policy: grid remains production default; TABI treecode default remains degree-one'
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII
@(
    '# TITAN FMM/M2L Scope Z4P Manuscript Guarded p4 Report',
    '',
    "Decision: $decisionValue",
    '',
    'Z4P refreshes the external manuscript with Scope Z4O repeatability evidence.',
    "Evidence source: $evidenceSource.",
    '',
    "- DOCX: $docxPath",
    "- Source: $sourcePath",
    "- Claim trace: $claimTraceCsv",
    "- QA: $qaCsv",
    '',
    'Policy:',
    '- Guarded p4 is described as a repeatable advanced explicit TABI treecode option.',
    '- Grid remains the production Poisson-Boltzmann default.',
    '- TABI treecode default remains degree-one.',
    '- Large-system verification is sampled direct-operator based.'
) | Set-Content -LiteralPath $rootReport -Encoding ASCII

if ($decisionValue -ne 'z4p_manuscript_ready') {
    throw "Scope Z4P did not pass: $decisionValue"
}
Write-Host "PASS Scope Z4P manuscript refresh: $docxPath"
