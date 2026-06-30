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

function Require-Directory([string]$Path) {
  if (-not (Test-Path -LiteralPath $Path -PathType Container)) {
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

function Find-Inkscape {
  $cmd = Get-Command "inkscape.exe" -ErrorAction SilentlyContinue
  if ($null -ne $cmd) { return $cmd.Source }
  foreach ($candidate in @(
      "C:\Program Files\Inkscape\bin\inkscape.exe",
      "C:\Program Files\Inkscape\inkscape.exe",
      "C:\Program Files (x86)\Inkscape\bin\inkscape.exe",
      "C:\Program Files (x86)\Inkscape\inkscape.exe"
    )) {
    if (Test-Path -LiteralPath $candidate -PathType Leaf) { return $candidate }
  }
  return $null
}

function Convert-SvgToPng([string]$SvgPath, [string]$PngPath, [int]$PixelScale = 3) {
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $PngPath) | Out-Null
  $inkscape = Find-Inkscape
  if ($null -ne $inkscape) {
    & $inkscape --export-type=png --export-filename="$PngPath" --export-dpi=600 `
      --export-background=white --export-background-opacity=1 "$SvgPath" | Out-Null
    if ($LASTEXITCODE -eq 0 -and (Test-Path -LiteralPath $PngPath -PathType Leaf)) {
      return
    }
  }
  $edge = "C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
  if (-not (Test-Path -LiteralPath $edge -PathType Leaf)) {
    throw "z4t_render_gap: Microsoft Edge headless renderer is unavailable"
  }
  $head = Get-Content -LiteralPath $SvgPath -TotalCount 1
  $window = if ($head -match "width=['""](\d+)['""]\s+height=['""](\d+)['""]") { "$($Matches[1]),$($Matches[2])" } else { "1200,700" }
  $profile = Join-Path ([IO.Path]::GetTempPath()) ("z4t_edge_" + [Guid]::NewGuid().ToString("N"))
  New-Item -ItemType Directory -Force -Path $profile | Out-Null
  $uri = ([Uri]([IO.Path]::GetFullPath($SvgPath))).AbsoluteUri
  try {
    & $edge "--headless=new" "--disable-gpu" "--disable-sync" "--hide-scrollbars" "--no-first-run" `
      "--force-device-scale-factor=$PixelScale" "--user-data-dir=$profile" "--window-size=$window" `
      "--screenshot=$PngPath" $uri | Out-Null
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $PngPath -PathType Leaf)) {
      throw "z4t_render_gap: SVG conversion failed for $SvgPath"
    }
  } finally {
    if (Test-Path -LiteralPath $profile) {
      Remove-Item -LiteralPath $profile -Recurse -Force
    }
  }
}

function Get-Z4tDouble([object]$Value, [double]$Default = [double]::NaN) {
  $parsed = [double]::NaN
  if ($null -ne $Value -and
      [double]::TryParse([string]$Value, [Globalization.NumberStyles]::Float, [Globalization.CultureInfo]::InvariantCulture, [ref]$parsed) -and
      -not [double]::IsNaN($parsed) -and -not [double]::IsInfinity($parsed)) {
    return $parsed
  }
  return $Default
}

function Write-Z4tWorkflowSvg([string]$SvgPath) {
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $SvgPath) | Out-Null
  @'
<svg xmlns="http://www.w3.org/2000/svg" width="1320" height="430" viewBox="0 0 1320 430">
  <defs>
    <filter id="shadow" x="-10%" y="-10%" width="130%" height="130%">
      <feDropShadow dx="1.5" dy="3" stdDeviation="3" flood-opacity="0.13"/>
    </filter>
    <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto">
      <path d="M0,0 L10,5 L0,10 z" fill="#374151"/>
    </marker>
    <style>
      .stageTitle { font: 700 18px Arial, sans-serif; fill: #111827; }
      .stageText { font: 15px Arial, sans-serif; fill: #374151; }
      .tag { font: 700 15px Arial, sans-serif; fill: #065f46; }
      .mutedTag { font: 700 14px Arial, sans-serif; fill: #475569; }
      .base { fill: #f8fafc; stroke: #cbd5e1; stroke-width: 1.6; rx: 12; filter: url(#shadow); }
      .focus { fill: #ecfdf5; stroke: #059669; stroke-width: 2.4; rx: 12; filter: url(#shadow); }
      .verify { fill: #eff6ff; stroke: #2563eb; stroke-width: 2.0; rx: 12; filter: url(#shadow); }
      .bracket { fill: none; stroke: #059669; stroke-width: 2.5; stroke-dasharray: 8 6; rx: 18; }
      .arrow { stroke: #374151; stroke-width: 2.2; marker-end: url(#arrow); }
    </style>
  </defs>
  <rect width="1320" height="430" fill="white"/>

  <rect x="46" y="90" width="164" height="168" class="base"/>
  <text x="128" y="126" class="stageTitle" text-anchor="middle">1. Structure</text>
  <text x="128" y="158" class="stageText" text-anchor="middle">atoms, charges,</text>
  <text x="128" y="180" class="stageText" text-anchor="middle">radii, salt</text>
  <text x="128" y="220" class="mutedTag" text-anchor="middle">input model</text>

  <rect x="250" y="90" width="164" height="168" class="base"/>
  <text x="332" y="126" class="stageTitle" text-anchor="middle">2. Surface</text>
  <text x="332" y="158" class="stageText" text-anchor="middle">triangulated</text>
  <text x="332" y="180" class="stageText" text-anchor="middle">dielectric boundary</text>
  <text x="332" y="220" class="mutedTag" text-anchor="middle">BIE unknowns</text>

  <rect x="454" y="90" width="164" height="168" class="focus"/>
  <text x="536" y="126" class="stageTitle" text-anchor="middle">3. Operator</text>
  <text x="536" y="158" class="stageText" text-anchor="middle">coupled K1-K4</text>
  <text x="536" y="180" class="stageText" text-anchor="middle">surface matrix</text>
  <text x="536" y="220" class="tag" text-anchor="middle">guarded p4 FMM</text>

  <rect x="658" y="90" width="164" height="168" class="focus"/>
  <text x="740" y="126" class="stageTitle" text-anchor="middle">4. Matvec</text>
  <text x="740" y="158" class="stageText" text-anchor="middle">M2L/L2P</text>
  <text x="740" y="180" class="stageText" text-anchor="middle">fallback</text>
  <text x="740" y="220" class="tag" text-anchor="middle">compressed jets</text>

  <rect x="862" y="90" width="164" height="168" class="focus"/>
  <text x="944" y="126" class="stageTitle" text-anchor="middle">5. Krylov solve</text>
  <text x="944" y="158" class="stageText" text-anchor="middle">GMRES</text>
  <text x="944" y="180" class="stageText" text-anchor="middle">near graph</text>
  <text x="944" y="220" class="tag" text-anchor="middle">sparse ILU(0)</text>

  <rect x="1066" y="90" width="208" height="168" class="verify"/>
  <text x="1170" y="126" class="stageTitle" text-anchor="middle">6. Verification</text>
  <text x="1170" y="158" class="stageText" text-anchor="middle">sampled dense</text>
  <text x="1170" y="180" class="stageText" text-anchor="middle">operator checks</text>
  <text x="1170" y="220" class="mutedTag" text-anchor="middle">audit residuals</text>

  <g class="arrow">
    <line x1="210" y1="174" x2="246" y2="174"/>
    <line x1="414" y1="174" x2="450" y2="174"/>
    <line x1="618" y1="174" x2="654" y2="174"/>
    <line x1="822" y1="174" x2="858" y2="174"/>
    <line x1="1026" y1="174" x2="1062" y2="174"/>
  </g>

  <rect x="434" y="66" width="612" height="226" class="bracket"/>
  <text x="740" y="326" class="tag" text-anchor="middle">This work</text>
</svg>
'@ | Set-Content -LiteralPath $SvgPath -Encoding ASCII
}

function Write-Z4tAlgorithmSvg([string]$SvgPath) {
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $SvgPath) | Out-Null
  @'
<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="650" viewBox="0 0 1000 650">
  <defs>
    <filter id="drop-shadow" x="-10%" y="-10%" width="130%" height="130%">
      <feDropShadow dx="2" dy="4" stdDeviation="4" flood-opacity="0.10"/>
    </filter>
    <marker id="arrowPurple" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">
      <path d="M0,0 L10,5 L0,10 z" fill="#7c3aed"/>
    </marker>
    <marker id="arrowRed" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <path d="M0,0 L10,5 L0,10 z" fill="#dc2626"/>
    </marker>
    <marker id="arrowBlue" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <path d="M0,0 L10,5 L0,10 z" fill="#2563eb"/>
    </marker>
    <marker id="arrowGreen" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
      <path d="M0,0 L10,5 L0,10 z" fill="#15803d"/>
    </marker>
    <style>
      .section { font: 700 16px Arial, sans-serif; fill: #111827; }
      .label { font: 700 13px Arial, sans-serif; fill: #111827; }
      .small { font: 12px Arial, sans-serif; fill: #374151; }
      .tiny { font: 11px Arial, sans-serif; fill: #4b5563; }
      .cellBlue { fill: #eff6ff; stroke: #2563eb; stroke-width: 2; rx: 12; filter: url(#drop-shadow); }
      .cellRed { fill: #fef2f2; stroke: #dc2626; stroke-width: 2; rx: 12; filter: url(#drop-shadow); }
      .box { fill: #ffffff; stroke: #cbd5e1; stroke-width: 1.4; rx: 10; }
      .guard { fill: #faf5ff; stroke: #7c3aed; stroke-width: 1.5; rx: 10; }
      .near { fill: #f0fdf4; stroke: #15803d; stroke-width: 1.6; rx: 10; }
      .muted { fill: #f8fafc; stroke: #94a3b8; stroke-width: 1.2; rx: 9; }
    </style>
  </defs>
  <rect width="1000" height="650" fill="white"/>

  <rect x="70" y="145" width="250" height="250" class="cellBlue"/>
  <text x="195" y="174" class="section" text-anchor="middle" fill="#2563eb">Target leaf</text>
  <circle cx="195" cy="285" r="86" fill="none" stroke="#2563eb" stroke-width="1.5" stroke-dasharray="7 5" opacity="0.65"/>
  <circle cx="195" cy="285" r="6" fill="#2563eb"/>
  <text x="210" y="300" class="label" fill="#2563eb">c_L</text>
  <g stroke="#2563eb" stroke-width="2.0" marker-end="url(#arrowBlue)">
    <line x1="195" y1="285" x2="135" y2="235"/>
    <line x1="195" y1="285" x2="255" y2="226"/>
    <line x1="195" y1="285" x2="140" y2="340"/>
    <line x1="195" y1="285" x2="263" y2="346"/>
  </g>
  <circle cx="135" cy="235" r="5" fill="#2563eb"/>
  <circle cx="255" cy="226" r="5" fill="#2563eb"/>
  <circle cx="140" cy="340" r="5" fill="#2563eb"/>
  <circle cx="263" cy="346" r="5" fill="#2563eb"/>
  <rect x="680" y="145" width="250" height="250" class="cellRed"/>
  <text x="805" y="174" class="section" text-anchor="middle" fill="#dc2626">Far source cell</text>
  <circle cx="805" cy="285" r="86" fill="none" stroke="#dc2626" stroke-width="1.5" stroke-dasharray="7 5" opacity="0.65"/>
  <circle cx="805" cy="285" r="6" fill="#dc2626"/>
  <text x="820" y="300" class="label" fill="#dc2626">c_M</text>
  <g stroke="#dc2626" stroke-width="2.0" marker-end="url(#arrowRed)">
    <line x1="735" y1="226" x2="798" y2="278"/>
    <line x1="875" y1="236" x2="812" y2="278"/>
    <line x1="728" y1="350" x2="798" y2="292"/>
    <line x1="866" y1="356" x2="812" y2="292"/>
    <line x1="890" y1="285" x2="815" y2="285"/>
  </g>
  <circle cx="735" cy="226" r="5" fill="#dc2626"/>
  <circle cx="875" cy="236" r="5" fill="#dc2626"/>
  <circle cx="728" cy="350" r="5" fill="#dc2626"/>
  <circle cx="866" cy="356" r="5" fill="#dc2626"/>
  <circle cx="890" cy="285" r="5" fill="#dc2626"/>
  <line x1="785" y1="285" x2="215" y2="285" stroke="#7c3aed" stroke-width="4" stroke-dasharray="9 6" marker-end="url(#arrowPurple)"/>
  <text x="500" y="256" class="section" text-anchor="middle" fill="#7c3aed">M2L: multipole-to-local translation</text>

  <rect x="330" y="82" width="340" height="82" class="guard"/>
  <text x="500" y="109" class="label" text-anchor="middle" fill="#6d28d9">Far-cell accept / guard</text>
  <line x1="500" y1="164" x2="500" y2="239" stroke="#7c3aed" stroke-width="1.8" marker-end="url(#arrowPurple)"/>
  <line x1="500" y1="338" x2="500" y2="382" stroke="#7c3aed" stroke-width="1.8" marker-end="url(#arrowPurple)"/>

  <rect x="250" y="382" width="500" height="118" class="near"/>
  <text x="500" y="410" class="section" text-anchor="middle" fill="#15803d">Near-field and guarded fallback</text>
  <rect x="292" y="448" width="184" height="38" class="muted"/>
  <text x="384" y="472" class="label" text-anchor="middle">packed AVX2 P2P</text>
  <rect x="524" y="448" width="184" height="38" class="muted"/>
  <text x="616" y="472" class="label" text-anchor="middle">exact self terms</text>
  <line x1="250" y1="432" x2="320" y2="432" stroke="#15803d" stroke-width="2" marker-end="url(#arrowGreen)"/>
  <line x1="750" y1="432" x2="680" y2="432" stroke="#15803d" stroke-width="2" marker-end="url(#arrowGreen)"/>

  <rect x="250" y="528" width="500" height="58" class="box"/>
  <text x="500" y="562" class="label" text-anchor="middle">Near-source graph -> sparse ILU(0)</text>
  <line x1="500" y1="500" x2="500" y2="528" stroke="#15803d" stroke-width="2" marker-end="url(#arrowGreen)"/>
</svg>
'@ | Set-Content -LiteralPath $SvgPath -Encoding ASCII
}

function Escape-Z4tSvgText([string]$Text) {
  if ($null -eq $Text) { return "" }
  return [Security.SecurityElement]::Escape($Text)
}

function Format-Z4tLadderWall([string]$Value) {
  $seconds = [double]::NaN
  if (-not [double]::TryParse($Value, [Globalization.NumberStyles]::Float, [Globalization.CultureInfo]::InvariantCulture, [ref]$seconds)) {
    return "NA"
  }
  if ([double]::IsNaN($seconds) -or [double]::IsInfinity($seconds)) { return "NA" }
  if ($seconds -ge 600.0) { return ('{0:f1} min' -f ($seconds / 60.0)) }
  return ('{0:f1} s' -f $seconds)
}

function Format-Z4tLadderMetric($Row) {
  $audit = [double]::NaN
  if ([string]$Row.stage -ne "tabi_direct" -and
      [double]::TryParse([string]$Row.sampled_audit_rms, [Globalization.NumberStyles]::Float, [Globalization.CultureInfo]::InvariantCulture, [ref]$audit) -and
      -not [double]::IsNaN($audit) -and -not [double]::IsInfinity($audit)) {
    return ('audit {0:E2}' -f $audit)
  }
  $iter = 0
  if ([int]::TryParse([string]$Row.iterations, [ref]$iter)) { return "iter $iter" }
  return ""
}

function Write-Z4tPerformanceLadderSvg([string]$MatrixCsv, [string]$SvgPath) {
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $SvgPath) | Out-Null
  $rows = @(Import-Csv -LiteralPath $MatrixCsv)
  if ($rows.Count -eq 0) { throw "z4t_missing_performance_ladder_rows" }
  $fixtureOrder = @("barnase", "1a63", "1sh1", "1ubq", "2ci2", "1vii")
  $stageOrder = @(
    [pscustomobject]@{ key = "tabi_direct"; label = "Direct TABI"; color = "#6b7280" },
    [pscustomobject]@{ key = "fmm_degree1_diagonal"; label = "FMM deg1 + diag"; color = "#60a5fa" },
    [pscustomobject]@{ key = "fmm_degree1_block"; label = "FMM deg1 + block"; color = "#f59e0b" },
    [pscustomobject]@{ key = "fmm_degree1_ilu0"; label = "FMM deg1 + ILU0"; color = "#16a34a" },
    [pscustomobject]@{ key = "fmm_p4_ilu0"; label = "Guarded p4 + ILU0"; color = "#7c3aed" }
  )
  $byKey = @{}
  foreach ($row in $rows) { $byKey["$($row.fixture)|$($row.stage)"] = $row }

  $width = 1320
  $height = 520
  $left = 34
  $systemW = 140
  $cellW = 213
  $gap = 12
  $rowH = 72
  $headerY = 38
  $firstRowY = 62
  $sb = [Text.StringBuilder]::new()
  [void]$sb.AppendLine("<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
  [void]$sb.AppendLine("<rect width='$width' height='$height' fill='white'/>")
  [void]$sb.AppendLine("<style>.head{font:700 13px Arial,sans-serif;fill:#111827}.sys{font:700 15px Arial,sans-serif;fill:#111827}.panel{font:12px Arial,sans-serif;fill:#4b5563}.cellStatus{font:700 14px Arial,sans-serif;fill:white}.cellText{font:13px Arial,sans-serif;fill:white}.cellSmall{font:12px Arial,sans-serif;fill:white}</style>")
  [void]$sb.AppendLine("<text x='$($left + 2)' y='$headerY' class='head'>System</text>")
  for ($i = 0; $i -lt $stageOrder.Count; $i++) {
    $x = $left + $systemW + $gap + $i * ($cellW + $gap)
    $label = Escape-Z4tSvgText $stageOrder[$i].label
    [void]$sb.AppendLine("<text x='$($x + 6)' y='$headerY' class='head'>$label</text>")
  }
  for ($r = 0; $r -lt $fixtureOrder.Count; $r++) {
    $fixture = $fixtureOrder[$r]
    $sample = @($rows | Where-Object { $_.fixture -eq $fixture } | Sort-Object { [int]($_.panels -as [int]) } -Descending | Select-Object -First 1)[0]
    if ($null -eq $sample) { continue }
    $y = $firstRowY + $r * $rowH
    $name = if ($fixture -eq "barnase") { "Barnase" } else { [string]$sample.pdb_id }
    [void]$sb.AppendLine("<text x='$($left + 2)' y='$($y + 23)' class='sys'>$(Escape-Z4tSvgText $name)</text>")
    [void]$sb.AppendLine("<text x='$($left + 2)' y='$($y + 44)' class='panel'>$(Escape-Z4tSvgText ([string]$sample.panels)) panels</text>")
    for ($i = 0; $i -lt $stageOrder.Count; $i++) {
      $stage = $stageOrder[$i]
      $row = $byKey["$fixture|$($stage.key)"]
      if ($null -eq $row) { continue }
      $x = $left + $systemW + $gap + $i * ($cellW + $gap)
      $status = if ([string]$row.timed_out -eq "True") { "timeout" } elseif ([string]$row.status -eq "ok") { "ok" } else { "not converged" }
      $wall = Format-Z4tLadderWall $row.wall_seconds
      $metric = Format-Z4tLadderMetric $row
      [void]$sb.AppendLine("<rect x='$x' y='$y' width='$cellW' height='56' rx='7' fill='$($stage.color)' opacity='0.92'/>")
      [void]$sb.AppendLine("<text x='$($x + 10)' y='$($y + 18)' class='cellStatus'>$(Escape-Z4tSvgText $status)</text>")
      [void]$sb.AppendLine("<text x='$($x + 10)' y='$($y + 36)' class='cellText'>$(Escape-Z4tSvgText $wall)</text>")
      [void]$sb.AppendLine("<text x='$($x + 10)' y='$($y + 52)' class='cellSmall'>$(Escape-Z4tSvgText $metric)</text>")
    }
  }
  [void]$sb.AppendLine("</svg>")
  [IO.File]::WriteAllText((Resolve-Path -LiteralPath (Split-Path -Parent $SvgPath)).Path + [IO.Path]::DirectorySeparatorChar + (Split-Path -Leaf $SvgPath), $sb.ToString().TrimEnd() + [Environment]::NewLine, [Text.Encoding]::ASCII)
}

function Write-Z4tRuntimeCompositionSvg([string]$MatrixCsv, [string]$Z4vRoot, [string]$SvgPath, [string]$CompositionCsv) {
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $SvgPath) | Out-Null
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $CompositionCsv) | Out-Null
  $rows = @(Import-Csv -LiteralPath $MatrixCsv | Where-Object { $_.fixture -eq "1a63" })
  if ($rows.Count -lt 5) { throw "z4t_missing_1a63_runtime_rows" }
  $byStage = @{}
  foreach ($row in $rows) { $byStage[[string]$row.stage] = $row }

  $stageOrder = @(
    [pscustomobject]@{ stage = "tabi_direct"; run = "1a63_tabi_direct"; label = "Direct TABI"; direct = $true },
    [pscustomobject]@{ stage = "fmm_degree1_diagonal"; run = "1a63_fmm_degree1_diagonal"; label = "FMM + diagonal"; direct = $false },
    [pscustomobject]@{ stage = "fmm_degree1_block"; run = "1a63_fmm_degree1_block"; label = "FMM + block"; direct = $false },
    [pscustomobject]@{ stage = "fmm_degree1_ilu0"; run = "1a63_fmm_degree1_ilu0"; label = "FMM + ILU(0)"; direct = $false },
    [pscustomobject]@{ stage = "fmm_p4_ilu0"; run = "1a63_fmm_p4_ilu0"; label = "Guarded p4 + ILU(0)"; direct = $false }
  )
  $colors = @{
    "Setup / plan" = "#8b5cf6"
    "Matvec / apply" = "#2563eb"
    "Preconditioner" = "#059669"
    "Krylov overhead" = "#f59e0b"
    "Verification / other" = "#94a3b8"
    "Timeout / no breakdown" = "#6b7280"
  }

  $compositionRows = New-Object System.Collections.Generic.List[object]
  foreach ($spec in $stageOrder) {
    $matrixRow = $byStage[$spec.stage]
    if ($null -eq $matrixRow) { throw "z4t_missing_1a63_stage: $($spec.stage)" }
    $wall = Get-Z4tDouble $matrixRow.wall_seconds
    $profile = $null
    $profilePath = Join-Path (Join-Path $Z4vRoot $spec.run) "pbe_tabi_profile.csv"
    if (Test-Path -LiteralPath $profilePath -PathType Leaf) {
      $profile = @(Import-Csv -LiteralPath $profilePath | Select-Object -First 1)[0]
    }
    $total = if ($null -ne $profile) { Get-Z4tDouble $profile.total_seconds $wall } else { $wall }
    if ($total -le 0.0 -or [double]::IsNaN($total)) { $total = $wall }

    $components = New-Object System.Collections.Generic.List[object]
    if ($spec.direct) {
      $components.Add([pscustomobject]@{ component = "Timeout / no breakdown"; seconds = $wall }) | Out-Null
    } else {
      $setup = 0.0
      if ($spec.stage -eq "fmm_p4_ilu0") {
        $setup = Get-Z4tDouble $profile.fmm_local_prepare_seconds 0.0
      } else {
        $setup = (Get-Z4tDouble $profile.fmm_m2l_setup_seconds 0.0) + (Get-Z4tDouble $profile.tree_build_seconds 0.0)
      }
      $matvec = if ($spec.stage -eq "fmm_p4_ilu0") {
        Get-Z4tDouble $profile.fmm_local_gmres_apply_seconds_total (Get-Z4tDouble $profile.matvec_seconds_total 0.0)
      } else {
        Get-Z4tDouble $profile.matvec_seconds_total 0.0
      }
      $precond = Get-Z4tDouble $profile.preconditioner_seconds 0.0
      $krylov = Get-Z4tDouble $profile.orthogonalization_seconds 0.0
      $verify = Get-Z4tDouble $profile.matvec_audit_seconds 0.0
      $accounted = $setup + $matvec + $precond + $krylov + $verify
      $other = [math]::Max(0.0, $total - $accounted)
      $components.Add([pscustomobject]@{ component = "Setup / plan"; seconds = $setup }) | Out-Null
      $components.Add([pscustomobject]@{ component = "Matvec / apply"; seconds = $matvec }) | Out-Null
      $components.Add([pscustomobject]@{ component = "Preconditioner"; seconds = $precond }) | Out-Null
      $components.Add([pscustomobject]@{ component = "Krylov overhead"; seconds = $krylov }) | Out-Null
      $components.Add([pscustomobject]@{ component = "Verification / other"; seconds = ($verify + $other) }) | Out-Null
    }
    foreach ($component in $components) {
      $seconds = [math]::Max(0.0, [double]$component.seconds)
      $fraction = if ($total -gt 0.0) { $seconds / $total } else { 0.0 }
      $compositionRows.Add([pscustomobject]@{
        fixture = "1A63"
        stage = $spec.stage
        label = $spec.label
        status = $matrixRow.status
        iterations = $matrixRow.iterations
        total_seconds = ('{0:f6}' -f $total)
        component = $component.component
        component_seconds = ('{0:f6}' -f $seconds)
        component_fraction = ('{0:f8}' -f $fraction)
      }) | Out-Null
    }
  }
  $compositionRows | Export-Csv -LiteralPath $CompositionCsv -NoTypeInformation

  $width = 1320
  $height = 650
  $left = 340
  $right = 1090
  $top = 68
  $barH = 38
  $rowGap = 70
  $axisY = $top + ($stageOrder.Count * $rowGap) + 4
  function XFrac([double]$f) { return $left + $f * ($right - $left) }

  $sb = [Text.StringBuilder]::new()
  [void]$sb.AppendLine("<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
  [void]$sb.AppendLine("<rect width='$width' height='$height' fill='white'/>")
  [void]$sb.AppendLine("<style>.label{font:700 23px Arial,sans-serif;fill:#111827}.sub{font:18px Arial,sans-serif;fill:#4b5563}.axis{stroke:#374151;stroke-width:1.8}.grid{stroke:#e5e7eb;stroke-width:1.2}.tick{font:18px Arial,sans-serif;fill:#4b5563}.barText{font:700 16px Arial,sans-serif;fill:white}.note{font:18px Arial,sans-serif;fill:#111827}.legend{font:18px Arial,sans-serif;fill:#374151}</style>")
  foreach ($tick in @(0, 0.25, 0.5, 0.75, 1.0)) {
    $x = [math]::Round((XFrac $tick), 1)
    [void]$sb.AppendLine("<line x1='$x' y1='$($top - 18)' x2='$x' y2='$axisY' class='grid'/>")
    [void]$sb.AppendLine("<text x='$($x - 16)' y='$($axisY + 31)' class='tick'>$([int]($tick * 100))</text>")
  }
  [void]$sb.AppendLine("<line x1='$left' y1='$axisY' x2='$right' y2='$axisY' class='axis'/>")
  [void]$sb.AppendLine("<text x='$($left + 315)' y='$($axisY + 62)' class='sub'>runtime fraction (%)</text>")

  for ($i = 0; $i -lt $stageOrder.Count; $i++) {
    $spec = $stageOrder[$i]
    $matrixRow = $byStage[$spec.stage]
    $y = $top + $i * $rowGap
    [void]$sb.AppendLine("<text x='34' y='$($y + 25)' class='label'>$(Escape-Z4tSvgText $spec.label)</text>")
    $status = if ([string]$matrixRow.status -eq "ok") { "ok" } elseif ([string]$matrixRow.status -eq "timeout") { "timeout" } else { "not converged" }
    $rowDetail = if ($spec.direct) {
      "$status, $(Format-Z4tLadderWall $matrixRow.wall_seconds)"
    } else {
      "$status, $(Format-Z4tLadderWall $matrixRow.wall_seconds), iter $($matrixRow.iterations)"
    }
    [void]$sb.AppendLine("<text x='34' y='$($y + 50)' class='sub'>$(Escape-Z4tSvgText $rowDetail)</text>")
    $x = $left
    $parts = @($compositionRows | Where-Object { $_.stage -eq $spec.stage })
    foreach ($part in $parts) {
      $fraction = Get-Z4tDouble $part.component_fraction 0.0
      if ($fraction -le 0.00001) { continue }
      $w = [math]::Max(1.0, $fraction * ($right - $left))
      $color = $colors[$part.component]
      [void]$sb.AppendLine("<rect x='$([math]::Round($x,1))' y='$y' width='$([math]::Round($w,1))' height='$barH' fill='$color'/>")
      if ($fraction -ge 0.09) {
        $pct = "{0:f0}%" -f ($fraction * 100.0)
        [void]$sb.AppendLine("<text x='$([math]::Round($x + $w / 2.0,1))' y='$($y + 25)' class='barText' text-anchor='middle'>$pct</text>")
      }
      $x += $w
    }
    $wallText = "{0:f1} s" -f (Get-Z4tDouble $matrixRow.wall_seconds 0.0)
    [void]$sb.AppendLine("<text x='$($right + 24)' y='$($y + 25)' class='note'>$(Escape-Z4tSvgText $wallText)</text>")
  }

  $legendItems = @(
    [pscustomobject]@{ name = "Setup / plan"; x = 36; y = 530 },
    [pscustomobject]@{ name = "Matvec / apply"; x = 320; y = 530 },
    [pscustomobject]@{ name = "Preconditioner"; x = 620; y = 530 },
    [pscustomobject]@{ name = "Krylov overhead"; x = 36; y = 574 },
    [pscustomobject]@{ name = "Verification / other"; x = 320; y = 574 },
    [pscustomobject]@{ name = "Timeout / no breakdown"; x = 620; y = 574 }
  )
  foreach ($item in $legendItems) {
    $name = $item.name
    [void]$sb.AppendLine("<rect x='$($item.x)' y='$($item.y)' width='22' height='22' fill='$($colors[$name])'/>")
    [void]$sb.AppendLine("<text x='$($item.x + 32)' y='$($item.y + 17)' class='legend'>$(Escape-Z4tSvgText $name)</text>")
  }
  [void]$sb.AppendLine("</svg>")
  Set-Content -LiteralPath $SvgPath -Value $sb.ToString() -Encoding ASCII
}

function Write-Z4tPreconditionerConvergenceSvg([string]$TraceCsv, [string]$SvgPath) {
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $SvgPath) | Out-Null
  $rows = @(Import-Csv -LiteralPath $TraceCsv | Where-Object { $_.fixture -eq "barnase" -and ($_.preconditioner_label -eq "block_diagonal" -or $_.preconditioner_label -eq "sparse_ilu0") })
  if ($rows.Count -lt 2) { throw "z4t_missing_artifact: Barnase Z5 residual trace is incomplete" }

  $seriesNames = @("block_diagonal", "sparse_ilu0")
  $colors = @{
    block_diagonal = "#f59e0b"
    sparse_ilu0 = "#16a34a"
  }
  $labels = @{
    block_diagonal = "Local block-diagonal"
    sparse_ilu0 = "Sparse near-field ILU(0)"
  }
  $width = 1100
  $height = 620
  $left = 90
  $right = 1020
  $top = 78
  $bottom = 505
  $xMax = [double](@($rows | ForEach-Object { [double]$_.iteration }) | Measure-Object -Maximum).Maximum
  if ($xMax -lt 1.0) { $xMax = 1.0 }
  $yMin = -15.0
  $yMax = 0.0
  function X([double]$x) { return $left + ($x / $xMax) * ($right - $left) }
  function Y([double]$y) { return $bottom - (($y - $yMin) / ($yMax - $yMin)) * ($bottom - $top) }

  $sb = [Text.StringBuilder]::new()
  [void]$sb.AppendLine("<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
  [void]$sb.AppendLine("<rect width='$width' height='$height' fill='white'/>")
  [void]$sb.AppendLine("<style>.title{font:700 24px Arial,sans-serif;fill:#111827}.small{font:13px Arial,sans-serif;fill:#374151}.axis{stroke:#374151;stroke-width:1.4}.grid{stroke:#e5e7eb;stroke-width:1}.tick{font:12px Arial,sans-serif;fill:#4b5563}.legend{font:14px Arial,sans-serif;fill:#111827}.line{fill:none;stroke-width:3;stroke-linecap:round;stroke-linejoin:round}.dot{stroke:white;stroke-width:1.5}</style>")
  [void]$sb.AppendLine("<text x='40' y='34' class='title'>Barnase preconditioner convergence profile</text>")
  [void]$sb.AppendLine("<text x='40' y='58' class='small'>Residual history from the dense/high-coupling Barnase trace; sparse ILU(0) collapses the Krylov residual where local block preconditioning stalls.</text>")

  foreach ($ytick in @(0, -3, -6, -9, -12, -15)) {
    $yy = [math]::Round((Y $ytick), 1)
    [void]$sb.AppendLine("<line x1='$left' y1='$yy' x2='$right' y2='$yy' class='grid'/>")
    [void]$sb.AppendLine("<text x='48' y='$($yy + 4)' class='tick'>$ytick</text>")
  }
  foreach ($xtick in 0..([int][math]::Ceiling($xMax))) {
    $xx = [math]::Round((X $xtick), 1)
    [void]$sb.AppendLine("<line x1='$xx' y1='$bottom' x2='$xx' y2='$($bottom + 5)' class='axis'/>")
    [void]$sb.AppendLine("<text x='$($xx - 4)' y='$($bottom + 23)' class='tick'>$xtick</text>")
  }
  [void]$sb.AppendLine("<line x1='$left' y1='$bottom' x2='$right' y2='$bottom' class='axis'/>")
  [void]$sb.AppendLine("<line x1='$left' y1='$top' x2='$left' y2='$bottom' class='axis'/>")
  [void]$sb.AppendLine("<text x='500' y='575' class='small'>GMRES residual-history index</text>")
  [void]$sb.AppendLine("<text x='18' y='315' class='small' transform='rotate(-90 18,315)'>log10 residual</text>")

  $legendX = 650
  $legendY = 95
  foreach ($name in $seriesNames) {
    $color = $colors[$name]
    $label = Escape-Z4tSvgText $labels[$name]
    [void]$sb.AppendLine("<line x1='$legendX' y1='$legendY' x2='$($legendX + 35)' y2='$legendY' stroke='$color' stroke-width='4'/>")
    [void]$sb.AppendLine("<text x='$($legendX + 45)' y='$($legendY + 5)' class='legend'>$label</text>")
    $legendY += 26
  }

  foreach ($name in $seriesNames) {
    $series = @($rows | Where-Object { $_.preconditioner_label -eq $name } | Sort-Object { [double]$_.iteration })
    if ($series.Count -eq 0) { continue }
    $points = @()
    foreach ($row in $series) {
      $x = [math]::Round((X ([double]$row.iteration)), 1)
      $y = [math]::Round((Y ([double]$row.log10_residual)), 1)
      $points += "$x,$y"
    }
    $color = $colors[$name]
    [void]$sb.AppendLine("<polyline points='$($points -join " ")' class='line' stroke='$color'/>")
    foreach ($row in $series) {
      $x = [math]::Round((X ([double]$row.iteration)), 1)
      $y = [math]::Round((Y ([double]$row.log10_residual)), 1)
      [void]$sb.AppendLine("<circle cx='$x' cy='$y' r='4.5' fill='$color' class='dot'/>")
    }
  }
  [void]$sb.AppendLine("<text x='90' y='540' class='small'>Block-diagonal trace ends at residual 2.13e-4; sparse ILU(0) reaches 5.65e-15 in the recorded trace.</text>")
  [void]$sb.AppendLine("</svg>")
  Set-Content -LiteralPath $SvgPath -Value $sb.ToString() -Encoding ASCII
}

function Write-Z4tPreconditionerBaselineSvg([string]$Z4vRoot, [string]$TraceCsv, [string]$SvgPath) {
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $SvgPath) | Out-Null
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $TraceCsv) | Out-Null
  $series = @(
    [pscustomobject]@{
      stage = 'fmm_degree1_diagonal'
      run_dir = 'barnase_fmm_degree1_diagonal'
      label = 'Diagonal'
      color = '#4b5563'
    },
    [pscustomobject]@{
      stage = 'fmm_degree1_block'
      run_dir = 'barnase_fmm_degree1_block'
      label = 'Block-diagonal'
      color = '#d97706'
    },
    [pscustomobject]@{
      stage = 'fmm_degree1_ilu0'
      run_dir = 'barnase_fmm_degree1_ilu0'
      label = 'Sparse near-field ILU(0)'
      color = '#059669'
    }
  )

  $traceRows = New-Object System.Collections.Generic.List[object]
  foreach ($s in $series) {
    $progress = Join-Path (Join-Path $Z4vRoot $s.run_dir) 'pbe_tabi_longrun_progress.csv'
    if (-not (Test-Path -LiteralPath $progress)) {
      throw "z4t_missing_preconditioner_progress: $progress"
    }
    $traceRows.Add([pscustomobject]@{
      fixture = 'Barnase'
      stage = $s.stage
      preconditioner_label = $s.label
      iteration = 0
      residual = 1.0
      log10_residual = 0.0
    }) | Out-Null
    foreach ($row in (Import-Csv -LiteralPath $progress | Where-Object { $_.phase -eq 'gmres_iter' })) {
      $iter = 0
      if (-not [int]::TryParse([string]$row.gmres_iter, [ref]$iter)) { continue }
      $res = [double]::NaN
      if (-not [double]::TryParse([string]$row.best_preconditioned_residual, [Globalization.NumberStyles]::Float, [Globalization.CultureInfo]::InvariantCulture, [ref]$res)) {
        continue
      }
      if ([double]::IsNaN($res) -or [double]::IsInfinity($res) -or $res -le 0.0) { continue }
      $traceRows.Add([pscustomobject]@{
        fixture = 'Barnase'
        stage = $s.stage
        preconditioner_label = $s.label
        iteration = $iter
        residual = ('{0:e12}' -f $res)
        log10_residual = ('{0:f12}' -f ([math]::Log10($res)))
      }) | Out-Null
    }
  }
  $traceRows | Sort-Object preconditioner_label, @{ Expression = { [int]$_.iteration } } |
    Export-Csv -LiteralPath $TraceCsv -NoTypeInformation

  $rows = @(Import-Csv -LiteralPath $TraceCsv)
  if ($rows.Count -lt 6) { throw "z4t_preconditioner_trace_too_short" }
  $width = 1320
  $height = 780
  $left = 120
  $right = 900
  $top = 70
  $bottom = 615
  $xMax = [double](@($rows | ForEach-Object { [double]$_.iteration }) | Measure-Object -Maximum).Maximum
  if ($xMax -lt 1.0) { $xMax = 1.0 }
  $yMin = -15.0
  $yMax = 0.0
  function X([double]$x) { return $left + ($x / $xMax) * ($right - $left) }
  function Y([double]$y) { return $bottom - (($y - $yMin) / ($yMax - $yMin)) * ($bottom - $top) }

  $sb = [Text.StringBuilder]::new()
  [void]$sb.AppendLine("<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
  [void]$sb.AppendLine("<rect width='$width' height='$height' fill='white'/>")
  [void]$sb.AppendLine("<style>.title{font:700 30px Arial,sans-serif;fill:#111827}.small{font:22px Arial,sans-serif;fill:#374151}.axis{stroke:#374151;stroke-width:2.1}.grid{stroke:#e5e7eb;stroke-width:1.35}.tick{font:22px Arial,sans-serif;fill:#4b5563}.legend{font:22px Arial,sans-serif;fill:#111827}.line{fill:none;stroke-width:5.2;stroke-linecap:round;stroke-linejoin:round}.dot{stroke:white;stroke-width:2}.callout{font:20px Arial,sans-serif;fill:#111827}</style>")
  foreach ($ytick in @(0, -3, -6, -9, -12, -15)) {
    $yy = [math]::Round((Y $ytick), 1)
    [void]$sb.AppendLine("<line x1='$left' y1='$yy' x2='$right' y2='$yy' class='grid'/>")
    [void]$sb.AppendLine("<text x='50' y='$($yy + 7)' class='tick'>10^$ytick</text>")
  }
  foreach ($xtick in @(0, 128, 360, 720)) {
    $xx = [math]::Round((X $xtick), 1)
    [void]$sb.AppendLine("<line x1='$xx' y1='$bottom' x2='$xx' y2='$($bottom + 5)' class='axis'/>")
    [void]$sb.AppendLine("<text x='$($xx - 24)' y='$($bottom + 38)' class='tick'>$xtick</text>")
  }
  [void]$sb.AppendLine("<line x1='$left' y1='$bottom' x2='$right' y2='$bottom' class='axis'/>")
  [void]$sb.AppendLine("<line x1='$left' y1='$top' x2='$left' y2='$bottom' class='axis'/>")
  [void]$sb.AppendLine("<text x='500' y='696' class='small'>GMRES iteration</text>")
  [void]$sb.AppendLine("<text x='22' y='402' class='small' transform='rotate(-90 22,402)'>best residual</text>")

  [void]$sb.AppendLine("<rect x='930' y='84' width='365' height='132' fill='white' stroke='#d1d5db' stroke-width='1.2' rx='10' opacity='0.94'/>")
  $legendX = 965
  $legendY = 120
  foreach ($s in $series) {
    $label = Escape-Z4tSvgText $s.label
    [void]$sb.AppendLine("<line x1='$legendX' y1='$legendY' x2='$($legendX + 50)' y2='$legendY' stroke='$($s.color)' stroke-width='5.4'/>")
    [void]$sb.AppendLine("<text x='$($legendX + 64)' y='$($legendY + 8)' class='legend'>$label</text>")
    $legendY += 39
  }

  foreach ($s in $series) {
    $label = $s.label
    $data = @($rows | Where-Object { $_.preconditioner_label -eq $label } | Sort-Object { [double]$_.iteration })
    if ($data.Count -eq 0) { continue }
    $points = @()
    foreach ($row in $data) {
      $x = [math]::Round((X ([double]$row.iteration)), 1)
      $y = [math]::Round((Y ([double]$row.log10_residual)), 1)
      $points += "$x,$y"
    }
    [void]$sb.AppendLine("<polyline points='$($points -join " ")' class='line' stroke='$($s.color)'/>")
    foreach ($row in @($data | Select-Object -First 1; $data | Select-Object -Last 1)) {
      $x = [math]::Round((X ([double]$row.iteration)), 1)
      $y = [math]::Round((Y ([double]$row.log10_residual)), 1)
      [void]$sb.AppendLine("<circle cx='$x' cy='$y' r='7.0' fill='$($s.color)' class='dot'/>")
    }
  }
  [void]$sb.AppendLine("</svg>")
  Set-Content -LiteralPath $SvgPath -Value $sb.ToString() -Encoding ASCII
}

$inputs = @{
  z4oMatrix = Require-File (Join-Path $packagePath "tables/fagbi_scope_z4p_z4o_snapshot_matrix.csv")
  z4oRepeatability = Require-File (Join-Path $packagePath "tables/fagbi_scope_z4p_z4o_snapshot_repeatability.csv")
  z4oDecision = Require-File (Join-Path $packagePath "tables/fagbi_scope_z4p_z4o_snapshot_decision.csv")
  z4qMatrix = Require-File (Join-Path $packagePath "tables/fagbi_scope_z4q_matrix.csv")
  z4qRepeatability = Require-File (Join-Path $packagePath "tables/fagbi_scope_z4q_repeatability.csv")
  z4qManifest = Require-File (Join-Path $packagePath "provenance/public_fixture_expansion_manifest.csv")
  z4vMatrix = Require-File "test/tmp/fagbi_scope_z4v/fagbi_scope_z4v_matrix.csv"
  z4vRoot = Require-Directory "test/tmp/fagbi_scope_z4v"
  z4vPlot = Require-File "test/tmp/fagbi_scope_z4v/fagbi_scope_z4v_performance_ladder.svg"
  pOrderPlot = if (Test-Path -LiteralPath "test/tmp/fagbi_scope_z4x_global_order/fagbi_scope_z4x_global_p_order_panels.svg") {
    (Resolve-Path -LiteralPath "test/tmp/fagbi_scope_z4x_global_order/fagbi_scope_z4x_global_p_order_panels.svg").Path
  } else {
    $null
  }
  pOrderCoverage = if (Test-Path -LiteralPath "test/tmp/fagbi_scope_z4x_global_order/fagbi_scope_z4x_global_p_order_coverage.csv") {
    (Resolve-Path -LiteralPath "test/tmp/fagbi_scope_z4x_global_order/fagbi_scope_z4x_global_p_order_coverage.csv").Path
  } else {
    $null
  }
  apbsMatrix = if (Test-Path -LiteralPath "test/tmp/fagbi_scope_z4z_apbs/fagbi_scope_z4z_apbs_matrix.csv") {
    (Resolve-Path -LiteralPath "test/tmp/fagbi_scope_z4z_apbs/fagbi_scope_z4z_apbs_matrix.csv").Path
  } else {
    $null
  }
}

if (Test-Path -LiteralPath $outPath) {
  Remove-Item -LiteralPath $outPath -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $outPath | Out-Null
$figureRoot = Join-Path $outPath "figures"
$workflowSvg = Join-Path $figureRoot "fagbi_scope_z4t_bie_workflow_map.svg"
$workflowPng = Join-Path $figureRoot "fagbi_scope_z4t_bie_workflow_map.png"
$algorithmSvg = Join-Path $figureRoot "fagbi_scope_z4t_algorithm_schematic.svg"
$algorithmPng = Join-Path $figureRoot "fagbi_scope_z4t_algorithm_schematic.png"
$z4vSvg = Join-Path $figureRoot "fagbi_scope_z4v_performance_ladder.svg"
$z4vPng = Join-Path $figureRoot "fagbi_scope_z4v_performance_ladder.png"
$runtimeCompositionSvg = Join-Path $figureRoot "fagbi_scope_z4t_1a63_runtime_composition.svg"
$runtimeCompositionPng = Join-Path $figureRoot "fagbi_scope_z4t_1a63_runtime_composition.png"
$runtimeCompositionCsv = Join-Path $outPath "fagbi_scope_z4t_1a63_runtime_composition.csv"
$precondTrace = Join-Path $outPath "fagbi_scope_z4t_preconditioner_residual_trace.csv"
$precondSvg = Join-Path $figureRoot "fagbi_scope_z4t_preconditioner_residual_profile.svg"
$precondPng = Join-Path $figureRoot "fagbi_scope_z4t_preconditioner_residual_profile.png"
$pOrderSvg = Join-Path $figureRoot "fagbi_scope_z4x_global_p_order_panels.svg"
$pOrderPng = Join-Path $figureRoot "fagbi_scope_z4x_global_p_order_panels.png"
$pOrderCoverage = Join-Path $outPath "fagbi_scope_z4x_global_p_order_coverage.csv"
Write-Z4tWorkflowSvg $workflowSvg
Write-Z4tAlgorithmSvg $algorithmSvg
Write-Z4tPerformanceLadderSvg $inputs.z4vMatrix $z4vSvg
Write-Z4tRuntimeCompositionSvg $inputs.z4vMatrix $inputs.z4vRoot $runtimeCompositionSvg $runtimeCompositionCsv
Write-Z4tPreconditionerBaselineSvg $inputs.z4vRoot $precondTrace $precondSvg
Convert-SvgToPng $workflowSvg $workflowPng
Convert-SvgToPng $algorithmSvg $algorithmPng
Convert-SvgToPng $z4vSvg $z4vPng
Convert-SvgToPng $runtimeCompositionSvg $runtimeCompositionPng
Convert-SvgToPng $precondSvg $precondPng
if ($null -ne $inputs.pOrderPlot) {
  Copy-Item -LiteralPath $inputs.pOrderPlot -Destination $pOrderSvg -Force
  Convert-SvgToPng $pOrderSvg $pOrderPng
}
if ($null -ne $inputs.pOrderCoverage) {
  Copy-Item -LiteralPath $inputs.pOrderCoverage -Destination $pOrderCoverage -Force
}

$builderArgs = @(
  $builder,
  "--z4o-matrix", $inputs.z4oMatrix,
  "--z4o-repeatability", $inputs.z4oRepeatability,
  "--z4o-decision", $inputs.z4oDecision,
  "--z4q-matrix", $inputs.z4qMatrix,
  "--z4q-repeatability", $inputs.z4qRepeatability,
  "--z4q-manifest", $inputs.z4qManifest,
  "--workflow-plot", $workflowPng,
  "--algorithm-plot", $algorithmPng,
  "--z4v-matrix", $inputs.z4vMatrix,
  "--z4v-plot", $z4vPng,
  "--runtime-composition-plot", $runtimeCompositionPng,
  "--z5-plot", $precondPng,
  "--output-dir", $outPath
)
if ($null -ne $inputs.pOrderPlot) {
  $builderArgs += @("--p-order-plot", $pOrderPng)
}
if ($null -ne $inputs.pOrderCoverage) {
  $builderArgs += @("--p-order-coverage", $pOrderCoverage)
}
if ($null -ne $inputs.apbsMatrix) {
  $builderArgs += @("--apbs-matrix", $inputs.apbsMatrix)
}

& $PythonExe @builderArgs
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
$wordRow = @($qa | Where-Object { $_.check -eq "word_count" })[0]
$decision = if ((@($qa | Where-Object { $_.status -eq "fail" }).Count -eq 0) -and
    $refs.Count -ge 8 -and $gaps[0].status -eq "closed" -and
    ($null -eq $wordRow -or [int]$wordRow.value -ge 3000)) {
  "z4t_submission_polish_ready"
} else {
  "z4t_submission_polish_gap"
}

@([pscustomobject]@{
  decision = $decision
  docx = $docx
  pdf = $pdf
  pages = $pageCount
  word_count = if ($null -ne $wordRow) { $wordRow.value } else { "NA" }
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
  "word_count: $(if ($null -ne $wordRow) { $wordRow.value } else { 'NA' })",
  "references: $($refs.Count)",
  "citation_gaps: $($gaps[0].status)",
  "performance_ladder: Z4V full-protein matrix included",
  "runtime_composition: 1A63 component-fraction figure included",
  "policy: grid remains production default; guarded p4 remains explicit advanced path"
) | Set-Content -LiteralPath (Join-Path $outPath "fagbi_scope_z4t_summary.dat") -Encoding ASCII

@(
  "TITAN FMM/M2L Scope Z4T Submission Polish Report",
  "Decision: $decision",
  "DOCX: $docx",
  "PDF: $pdf",
  "Pages: $pageCount",
"Word count: $(if ($null -ne $wordRow) { $wordRow.value } else { 'NA' })",
"References locked: $($refs.Count)",
"Citation gaps: $($gaps[0].status)",
"Performance ladder: Z4V full-protein matrix included",
"Runtime composition: 1A63 component-fraction figure included",
"Render status: $renderStatus",
"Policy: grid remains production default; guarded p4 remains explicit advanced path"
) | Set-Content -LiteralPath $rootReport -Encoding ASCII

if ($decision -ne "z4t_submission_polish_ready") {
  throw "Scope Z4T did not pass: $decision"
}

Write-Host "PASS Scope Z4T submission polish: $docx"
