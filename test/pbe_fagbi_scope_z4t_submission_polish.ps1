param(
  [string]$PackageRoot = "docs/fagbi_guarded_p4_package",
  [string]$OutputRoot = "docs/fagbi_guarded_p4_package/submission",
  [string]$PythonExe = "C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe",
  [switch]$SkipRender
)

$ErrorActionPreference = "Stop"
$repoRoot = (Resolve-Path ".").Path
$packagePath = Resolve-Path -LiteralPath $PackageRoot
$outPath = Join-Path $repoRoot $OutputRoot
$builder = Join-Path $repoRoot "test/pbe_fagbi_scope_z4t_build_docx.py"
$renderRoot = Join-Path $outPath "render_qa"
$rootReport = Join-Path $repoRoot "TITAN_FMM_M2L_SCOPE_Z4T_SUBMISSION_POLISH_REPORT.txt"

function Require-File([string]$Path) {
  if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
    throw "z4t_missing_artifact: $Path"
  }
  return (Resolve-Path -LiteralPath $Path).Path
}

function Export-DocxToPdf([string]$InputDocx, [string]$OutputPdf) {
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $OutputPdf) | Out-Null
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
  if (-not (Test-Path -LiteralPath $OutputPdf -PathType Leaf)) { throw "z4t_render_gap: PDF export failed" }
}

function Render-PdfPages([string]$PdfPath, [string]$RenderDir) {
  New-Item -ItemType Directory -Force -Path $RenderDir | Out-Null
  $script = Join-Path $RenderDir "render_pdf_pages.py"
  $code = @'
from pathlib import Path
import sys
import pypdfium2 as pdfium
pdf = pdfium.PdfDocument(str(Path(sys.argv[1])))
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
for index in range(len(pdf)):
    page = pdf[index]
    image = page.render(scale=2.0).to_pil()
    image.save(out / f"page-{index + 1:03d}.png")
    page.close()
print(len(pdf))
'@
  Set-Content -LiteralPath $script -Value $code -Encoding ASCII
  $countText = (& $PythonExe $script $PdfPath $RenderDir | Select-Object -Last 1)
  if ($LASTEXITCODE -ne 0) { throw "z4t_render_gap: PDF page render failed" }
  $count = [int]$countText
  $pages = @(Get-ChildItem -LiteralPath $RenderDir -Filter "page-*.png")
  if ($pages.Count -ne $count -or $count -lt 1) { throw "z4t_render_gap: expected $count pages, got $($pages.Count)" }
  return $count
}

function Convert-SvgToPng([string]$SvgPath, [string]$PngPath) {
  $edge = "C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
  if (-not (Test-Path -LiteralPath $edge -PathType Leaf)) {
    throw "z4t_render_gap: Microsoft Edge headless renderer is unavailable"
  }
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $PngPath) | Out-Null
  $head = Get-Content -LiteralPath $SvgPath -TotalCount 1
  $window = if ($head -match 'width="(\d+)"\s+height="(\d+)"') { "$($Matches[1]),$($Matches[2])" } else { "900,420" }
  $profile = Join-Path ([IO.Path]::GetTempPath()) ("z4t_edge_" + [Guid]::NewGuid().ToString("N"))
  New-Item -ItemType Directory -Force -Path $profile | Out-Null
  $uri = ([Uri]([IO.Path]::GetFullPath($SvgPath))).AbsoluteUri
  try {
    & $edge "--headless=new" "--disable-gpu" "--disable-sync" "--hide-scrollbars" "--no-first-run" `
      "--user-data-dir=$profile" "--window-size=$window" "--screenshot=$PngPath" $uri | Out-Null
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $PngPath -PathType Leaf)) {
      throw "z4t_render_gap: SVG conversion failed for $SvgPath"
    }
  } finally {
    if (Test-Path -LiteralPath $profile) {
      Remove-Item -LiteralPath $profile -Recurse -Force
    }
  }
}

$inputs = @{
  z4oMatrix = Require-File (Join-Path $packagePath "tables/fagbi_scope_z4p_z4o_snapshot_matrix.csv")
  z4oRepeatability = Require-File (Join-Path $packagePath "tables/fagbi_scope_z4p_z4o_snapshot_repeatability.csv")
  z4oDecision = Require-File (Join-Path $packagePath "tables/fagbi_scope_z4p_z4o_snapshot_decision.csv")
  z4qMatrix = Require-File (Join-Path $packagePath "tables/fagbi_scope_z4q_matrix.csv")
  z4qRepeatability = Require-File (Join-Path $packagePath "tables/fagbi_scope_z4q_repeatability.csv")
  z4qManifest = Require-File (Join-Path $packagePath "provenance/public_fixture_expansion_manifest.csv")
  z4qPlot = Require-File (Join-Path $packagePath "figures/fagbi_scope_z4q_repeatability_plot.svg")
  z5Plot = Require-File "test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_money_plot.svg"
}

if (Test-Path -LiteralPath $outPath) {
  Remove-Item -LiteralPath $outPath -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $outPath | Out-Null
$figureRoot = Join-Path $outPath "figures"
$z4qPng = Join-Path $figureRoot "fagbi_scope_z4q_repeatability_plot.png"
$z5Png = Join-Path $figureRoot "fagbi_scope_z5_preconditioner_profile.png"
Convert-SvgToPng $inputs.z4qPlot $z4qPng
Convert-SvgToPng $inputs.z5Plot $z5Png

& $PythonExe $builder `
  --z4o-matrix $inputs.z4oMatrix `
  --z4o-repeatability $inputs.z4oRepeatability `
  --z4o-decision $inputs.z4oDecision `
  --z4q-matrix $inputs.z4qMatrix `
  --z4q-repeatability $inputs.z4qRepeatability `
  --z4q-manifest $inputs.z4qManifest `
  --z4q-plot $z4qPng `
  --z5-plot $z5Png `
  --output-dir $outPath
if ($LASTEXITCODE -ne 0) { throw "z4t_docx_builder_failed" }

$docx = Join-Path $outPath "titan_fagbi_submission_polished.docx"
$pdf = Join-Path $outPath "titan_fagbi_submission_polished.pdf"
$renderStatus = "skipped"
$pageCount = 0
if (-not $SkipRender) {
  Export-DocxToPdf $docx $pdf
  $pageCount = Render-PdfPages $pdf $renderRoot
  $renderStatus = "pass"
}

$qa = Import-Csv (Join-Path $outPath "fagbi_scope_z4t_docx_qa.csv")
$qa += [pscustomobject]@{ check = "render_pages"; status = if ($SkipRender -or $renderStatus -eq "pass") { "pass" } else { "fail" }; value = $pageCount }
$qa | Export-Csv -LiteralPath (Join-Path $outPath "fagbi_scope_z4t_docx_qa.csv") -NoTypeInformation

$refs = Import-Csv (Join-Path $outPath "fagbi_scope_z4t_references.csv")
$gaps = Import-Csv (Join-Path $outPath "fagbi_scope_z4t_citation_gaps.csv")
$decision = if ((@($qa | Where-Object { $_.status -eq "fail" }).Count -eq 0) -and
    $refs.Count -ge 8 -and $gaps[0].status -eq "closed") {
  "z4t_submission_polish_ready"
} else {
  "z4t_submission_polish_gap"
}

@([pscustomobject]@{
  decision = $decision
  docx = $docx
  pdf = $pdf
  pages = $pageCount
  references = $refs.Count
  citation_gaps = $gaps[0].status
  render_status = $renderStatus
  grid_default = "unchanged"
  guarded_p4_policy = "advanced_explicit"
}) | Export-Csv -LiteralPath (Join-Path $outPath "fagbi_scope_z4t_decision.csv") -NoTypeInformation

@(
  "scope: Z4T submission polish",
  "decision: $decision",
  "docx: $docx",
  "pdf: $pdf",
  "pages: $pageCount",
  "references: $($refs.Count)",
  "citation_gaps: $($gaps[0].status)",
  "policy: grid remains production default; guarded p4 remains explicit advanced path"
) | Set-Content -LiteralPath (Join-Path $outPath "fagbi_scope_z4t_summary.dat") -Encoding ASCII

@(
  "TITAN FMM/M2L Scope Z4T Submission Polish Report",
  "Decision: $decision",
  "DOCX: $docx",
  "PDF: $pdf",
  "Pages: $pageCount",
  "References locked: $($refs.Count)",
  "Citation gaps: $($gaps[0].status)",
  "Render status: $renderStatus",
  "Policy: grid remains production default; guarded p4 remains explicit advanced path"
) | Set-Content -LiteralPath $rootReport -Encoding ASCII

if ($decision -ne "z4t_submission_polish_ready") {
  throw "Scope Z4T did not pass: $decision"
}

Write-Host "PASS Scope Z4T submission polish: $docx"
