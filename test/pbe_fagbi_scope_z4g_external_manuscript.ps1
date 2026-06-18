param(
    [string]$BenchmarkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4g'),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4g_manuscript'),
    [string]$PythonExe = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [string]$EdgeExe = 'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe',
    [switch]$SkipVisualQa,
    [switch]$VisualQaPassed
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Builder = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4g_build_docx.py'
$Matrix = Join-Path $BenchmarkRoot 'fagbi_scope_z4g_matrix.csv'
$Trace = Join-Path $BenchmarkRoot 'fagbi_scope_z4g_convergence_trace.csv'
$BenchmarkDecision = Join-Path $BenchmarkRoot 'fagbi_scope_z4g_decision.csv'
$RuntimeSvg = Join-Path $BenchmarkRoot 'fagbi_scope_z4g_runtime_memory.svg'
$ConvergenceSvg = Join-Path $BenchmarkRoot 'fagbi_scope_z4g_preconditioner_convergence.svg'
$Provenance = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$Docx = Join-Path $WorkRoot 'titan_fagbi_external_manuscript_z4g.docx'
$Source = Join-Path $WorkRoot 'fagbi_scope_z4g_manuscript_source.md'
$ClaimTrace = Join-Path $WorkRoot 'fagbi_scope_z4g_claim_trace.csv'
$Qa = Join-Path $WorkRoot 'fagbi_scope_z4g_docx_qa.csv'
$Checksums = Join-Path $WorkRoot 'fagbi_scope_z4g_checksums.csv'
$Decision = Join-Path $WorkRoot 'fagbi_scope_z4g_decision.csv'
$Summary = Join-Path $WorkRoot 'fagbi_scope_z4g_summary.dat'
$RenderRoot = Join-Path $WorkRoot 'render_qa'
$Pdf = Join-Path $RenderRoot 'titan_fagbi_external_manuscript_z4g.pdf'
$Report = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4G_EXTERNAL_MANUSCRIPT_REPORT.txt'

if ($SkipVisualQa -and $VisualQaPassed) {
    throw '-SkipVisualQa and -VisualQaPassed cannot be used together.'
}
foreach ($path in @($PythonExe, $Builder, $Matrix, $Trace, $BenchmarkDecision, $RuntimeSvg, $ConvergenceSvg, $Provenance)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "z4g_missing_artifacts: $path. Run test\pbe_fagbi_scope_z4g_full_protein_benchmark.ps1 first."
    }
}
$benchmarkRows = @(Import-Csv -LiteralPath $BenchmarkDecision)
if ($benchmarkRows.Count -ne 1 -or $benchmarkRows[0].decision -ne 'z4g_benchmark_ready') {
    throw "z4g_operator_accuracy_gap: benchmark decision is '$($benchmarkRows[0].decision)'."
}

$matrixRows = @(Import-Csv -LiteralPath $Matrix)
$fmmRows = @($matrixRows | Where-Object { $_.backend -eq 'fmm_ilu0' })
$representative = @($fmmRows | Group-Object fixture | ForEach-Object { $_.Group | Sort-Object { [int]$_.repeat } | Select-Object -First 1 })
if ($representative.Count -ne 3) {
    throw "z4g_fixture_preparation_gap: expected three complete-protein FMM rows, found $($representative.Count)."
}
foreach ($row in $fmmRows) {
    if ([int]$row.input_atoms -ne [int]$row.analyzed_atoms -or [int]$row.input_atoms -lt 1) {
        throw "z4g_fixture_preparation_gap: atom sampling detected for $($row.fixture) repeat $($row.repeat)."
    }
}

$workFull = [System.IO.Path]::GetFullPath($WorkRoot)
$tmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workFull.StartsWith($tmpRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Z4G manuscript work root outside test/tmp: $workFull"
}
if (Test-Path -LiteralPath $workFull) { Remove-Item -LiteralPath $workFull -Recurse -Force }
New-Item -ItemType Directory -Path $workFull -Force | Out-Null

function Convert-Z4gSvgToPng([string]$SvgPath, [string]$PngPath, [string]$ProfileName) {
    if (-not (Test-Path -LiteralPath $EdgeExe)) { throw "Missing Microsoft Edge: $EdgeExe" }
    $profile = Join-Path $WorkRoot "edge_$ProfileName"
    New-Item -ItemType Directory -Path $profile -Force | Out-Null
    $head = Get-Content -LiteralPath $SvgPath -TotalCount 1
    if ($head -notmatch 'width="(\d+)"\s+height="(\d+)"') { throw "Cannot resolve SVG dimensions: $SvgPath" }
    $uri = ([System.Uri]([System.IO.Path]::GetFullPath($SvgPath))).AbsoluteUri
    $args = @('--headless=new','--disable-gpu','--disable-sync','--hide-scrollbars','--no-first-run',
        "--user-data-dir=$profile", "--window-size=$($Matches[1]),$($Matches[2])", "--screenshot=$PngPath", $uri)
    & $EdgeExe @args | Out-Null
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $PngPath)) { throw "Failed to rasterize SVG: $SvgPath" }
}

$RuntimePng = Join-Path $WorkRoot 'figure1_runtime_memory.png'
$ConvergencePng = Join-Path $WorkRoot 'figure2_preconditioner_convergence.png'
Convert-Z4gSvgToPng $RuntimeSvg $RuntimePng 'runtime'
Convert-Z4gSvgToPng $ConvergenceSvg $ConvergencePng 'convergence'

& $PythonExe $Builder --matrix $Matrix --provenance $Provenance --trace $Trace `
    --runtime-figure $RuntimePng --convergence-figure $ConvergencePng `
    --output-docx $Docx --output-source $Source --claim-trace $ClaimTrace --qa $Qa
if ($LASTEXITCODE -ne 0) { throw "Z4G DOCX builder failed with exit code $LASTEXITCODE" }

function Export-Z4gDocxToPdf([string]$InputDocx, [string]$OutputPdf) {
    New-Item -ItemType Directory -Path (Split-Path -Parent $OutputPdf) -Force | Out-Null
    $word = $null; $document = $null
    try {
        $word = New-Object -ComObject Word.Application
        $word.Visible = $false; $word.DisplayAlerts = 0
        $document = $word.Documents.Open([System.IO.Path]::GetFullPath($InputDocx), $false, $true)
        $document.ExportAsFixedFormat([System.IO.Path]::GetFullPath($OutputPdf), 17)
    } finally {
        if ($null -ne $document) { $document.Close($false); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($document) }
        if ($null -ne $word) { $word.Quit(); [void][Runtime.InteropServices.Marshal]::ReleaseComObject($word) }
        [GC]::Collect(); [GC]::WaitForPendingFinalizers()
    }
    if (-not (Test-Path -LiteralPath $OutputPdf)) { throw "Word did not export PDF: $OutputPdf" }
}

function Render-Z4gPdfPages([string]$InputPdf) {
    $script = Join-Path $RenderRoot 'render_pdf_pages.py'
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
    $countText = (& $PythonExe $script $InputPdf $RenderRoot | Select-Object -Last 1).Trim()
    if ($LASTEXITCODE -ne 0) { throw 'z4g_manuscript_qa_gap: PDF page rendering failed.' }
    $captures = @(Get-ChildItem -LiteralPath $RenderRoot -Filter 'page-*.png' | Sort-Object Name)
    if ($captures.Count -ne [int]$countText -or $captures.Count -lt 1) { throw 'z4g_manuscript_qa_gap: incomplete page captures.' }
    return @($captures.FullName)
}

$renderStatus = 'skipped'; $pages = @()
if (-not $SkipVisualQa) {
    Export-Z4gDocxToPdf $Docx $Pdf
    $pages = @(Render-Z4gPdfPages $Pdf)
    $renderStatus = if ($VisualQaPassed) { 'pass' } else { 'pending_visual_review' }
}

$qaRows = @(Import-Csv -LiteralPath $Qa)
$qaRows += [pscustomobject]@{ check='full_protein_atom_identity'; status='pass'; value=$fmmRows.Count }
$qaRows += [pscustomobject]@{ check='render_page_count'; status=if($SkipVisualQa){'skipped'}else{'pass'}; value=$pages.Count }
$qaRows += [pscustomobject]@{ check='visual_review'; status=$renderStatus; value=$renderStatus }
$qaRows | Export-Csv -LiteralPath $Qa -NoTypeInformation

$checksumRows = foreach ($item in @($Matrix,$Trace,$Provenance,$RuntimeSvg,$ConvergenceSvg,$Docx,$Source,$ClaimTrace,$Qa)) {
    $file = Get-Item -LiteralPath $item
    [pscustomobject]@{ artifact=$file.Name; path=$file.FullName; bytes=$file.Length; sha256=(Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant() }
}
$checksumRows | Export-Csv -LiteralPath $Checksums -NoTypeInformation

$failures = @($qaRows | Where-Object { $_.status -eq 'fail' })
$decisionValue = if ($failures.Count -gt 0) {
    'z4g_manuscript_qa_gap'
} elseif ($renderStatus -ne 'pass') {
    'z4g_manuscript_qa_gap'
} else {
    'z4g_external_manuscript_ready'
}
[pscustomobject]@{
    decision=$decisionValue; complete_proteins=$representative.Count; fmm_rows=$fmmRows.Count;
    page_count=$pages.Count; visual_review=$renderStatus; manuscript_docx=$Docx
} | Export-Csv -LiteralPath $Decision -NoTypeInformation
@(
    'scope: Z4G public full-protein external manuscript',
    "decision: $decisionValue",
    "complete_proteins: $($representative.Count)",
    "fmm_rows: $($fmmRows.Count)",
    "page_count: $($pages.Count)",
    "visual_review: $renderStatus",
    'sampling_policy: complete prepared PQR structures; no atom sampling',
    'default_policy: grid remains production default'
) | Set-Content -LiteralPath $Summary -Encoding ASCII
@(
    '# TITAN FMM/M2L Scope Z4G External Manuscript Report','',
    "Decision: $decisionValue",'',
    'Z4G replaces internal sampled development labels with reproducible complete-protein evidence for Barnase, 1A63, and 1SH1.',
    '', 'Artifacts:',
    "- DOCX: $Docx", "- PDF render: $Pdf", "- Claim trace: $ClaimTrace", "- Checksums: $Checksums", "- Decision: $Decision",
    '', 'Policy:',
    '- Grid remains the production Poisson-Boltzmann default.',
    '- The manuscript reports complete prepared PQR structures; atom-sampled development energies are excluded.',
    '- Sparse near-field ILU(0) is described with its 2 x 2 block graph and no-fill triangular solve.',
    '- The Barnase preconditioner statement is limited to the measured complete-protein convergence profile.'
) | Set-Content -LiteralPath $Report -Encoding ASCII

Write-Host "PASS Scope Z4G external manuscript decision=$decisionValue docx=$Docx"
