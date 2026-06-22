param(
    [string]$PythonExe = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4f_render_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$Harness = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4f_manuscript_skeleton.ps1'
& $Harness -PythonExe $PythonExe -WorkRoot $WorkRoot

$decision = @(Import-Csv -LiteralPath (Join-Path $WorkRoot 'fagbi_scope_z4f_decision.csv'))
if ($decision.Count -ne 1 -or $decision[0].decision -ne 'z4f_render_gap') {
    throw "Render smoke expected pending visual review, got '$($decision[0].decision)'"
}
$pages = @(Get-ChildItem -LiteralPath (Join-Path $WorkRoot 'render_qa') -Filter 'page-*.png')
if ($pages.Count -lt 1) { throw 'Render smoke produced no page PNGs.' }

$pixelCheck = @'
from PIL import Image, ImageStat
import sys

for path in sys.argv[1:]:
    image = Image.open(path).convert("L")
    lo, hi = image.getextrema()
    stat = ImageStat.Stat(image)
    if hi - lo < 40 or stat.stddev[0] < 8.0:
        raise SystemExit(f"blank or low-information render: {path} range={hi-lo} stddev={stat.stddev[0]:.3f}")
print(f"rendered pages carry content: {len(sys.argv)-1}")
'@
$pixelCheckPath = Join-Path $WorkRoot 'check_render_pixels.py'
Set-Content -LiteralPath $pixelCheckPath -Value $pixelCheck -Encoding ASCII
& $PythonExe $pixelCheckPath @($pages.FullName)
if ($LASTEXITCODE -ne 0) { throw 'Rendered page content check failed.' }

Write-Host 'PASS pbe_fagbi_scope_z4f_render_smoke'
