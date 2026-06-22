param(
    [string]$PythonExe = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4f_contract')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Harness = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4f_manuscript_skeleton.ps1'
$Builder = Join-Path $PSScriptRoot 'pbe_fagbi_scope_z4f_build_docx.py'
if (-not (Test-Path -LiteralPath $Harness)) { throw "Missing Z4F harness: $Harness" }
if (-not (Test-Path -LiteralPath $Builder)) { throw "Missing Z4F DOCX builder: $Builder" }
if (-not (Test-Path -LiteralPath $PythonExe)) { throw "Missing bundled Python: $PythonExe" }

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$packageRoot = Join-Path $WorkRoot 'package'
& $Harness -PythonExe $PythonExe -WorkRoot $packageRoot -SkipVisualQa

$expected = @(
    'titan_fagbi_manuscript_skeleton.docx',
    'fagbi_scope_z4f_manuscript_source.md',
    'fagbi_scope_z4f_claim_trace.csv',
    'fagbi_scope_z4f_citation_gaps.csv',
    'fagbi_scope_z4f_docx_qa.csv',
    'fagbi_scope_z4f_decision.csv',
    'fagbi_scope_z4f_summary.dat'
)
foreach ($name in $expected) {
    $path = Join-Path $packageRoot $name
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing Z4F artifact: $path" }
}

$decision = @(Import-Csv -LiteralPath (Join-Path $packageRoot 'fagbi_scope_z4f_decision.csv'))
if ($decision.Count -ne 1 -or $decision[0].decision -ne 'z4f_render_gap') {
    throw "SkipVisualQa must classify z4f_render_gap, got '$($decision[0].decision)'"
}
if ([int]$decision[0].claim_count -ne 6 -or [int]$decision[0].limitation_count -ne 5) {
    throw 'Z4F claim/limitation coverage drifted.'
}
if ([int]$decision[0].protein_row_count -ne 6 -or [int]$decision[0].figure_count -ne 2) {
    throw 'Z4F manuscript table/figure contract drifted.'
}

Add-Type -AssemblyName System.IO.Compression.FileSystem
$docx = Join-Path $packageRoot 'titan_fagbi_manuscript_skeleton.docx'
$zip = [System.IO.Compression.ZipFile]::OpenRead($docx)
try {
    $documentEntry = $zip.GetEntry('word/document.xml')
    if ($null -eq $documentEntry) { throw 'DOCX has no word/document.xml' }
    $reader = New-Object System.IO.StreamReader($documentEntry.Open())
    try { $documentXml = $reader.ReadToEnd() } finally { $reader.Dispose() }
    foreach ($needle in @(
            'TITAN: Analytical FMM/M2L Boundary-Integral Poisson-Boltzmann Solves',
            'Abstract','Introduction','Methods','Results','Discussion','Limitations',
            'Conclusion','Code and Data Availability')) {
        if ($documentXml -notlike "*$needle*") { throw "DOCX missing manuscript text: $needle" }
    }
    foreach ($forbidden in @('TBD','high_resolution_analytic_convergence','preconditioner_money_plot')) {
        if ($documentXml -like "*$forbidden*") { throw "DOCX leaked forbidden token: $forbidden" }
    }
    $pageBreakCount = ([regex]::Matches($documentXml, '<w:br w:type="page"')).Count
    if ($pageBreakCount -lt 4) { throw "DOCX expected at least four deliberate page breaks, found $pageBreakCount" }
    $cantSplitCount = ([regex]::Matches($documentXml, '<w:cantSplit')).Count
    if ($cantSplitCount -lt 17) { throw "DOCX table rows are not protected from splitting: $cantSplitCount" }
    $mediaCount = @($zip.Entries | Where-Object { $_.FullName -like 'word/media/*' }).Count
    if ($mediaCount -ne 2) { throw "DOCX expected two embedded figures, found $mediaCount" }
} finally {
    $zip.Dispose()
}

Add-Type -AssemblyName System.Drawing
$figureContract = @(
    @{ name='figure1_runtime_memory.png'; width=1040; height=420 },
    @{ name='figure2_convergence.png'; width=900; height=520 }
)
foreach ($expectedFigure in $figureContract) {
    $figurePath = Join-Path $packageRoot $expectedFigure.name
    $image = [System.Drawing.Image]::FromFile($figurePath)
    try {
        if ($image.Width -ne $expectedFigure.width -or $image.Height -ne $expectedFigure.height) {
            throw "$($expectedFigure.name) expected $($expectedFigure.width)x$($expectedFigure.height), got $($image.Width)x$($image.Height)"
        }
    } finally {
        $image.Dispose()
    }
}

$z4eCopy = Join-Path $WorkRoot 'z4e_tampered'
Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4e') `
    -Destination $z4eCopy -Recurse
$checksumPath = Join-Path $z4eCopy 'fagbi_scope_z4e_checksums.csv'
$checksums = @(Import-Csv -LiteralPath $checksumPath)
$checksums[0].sha256 = ('0' * 64)
$checksums | Export-Csv -LiteralPath $checksumPath -NoTypeInformation
$driftRoot = Join-Path $WorkRoot 'drift'
$driftStdout = Join-Path $WorkRoot 'drift_stdout.txt'
$driftStderr = Join-Path $WorkRoot 'drift_stderr.txt'
$driftProcess = Start-Process -FilePath 'powershell.exe' -WindowStyle Hidden -Wait -PassThru `
    -RedirectStandardOutput $driftStdout -RedirectStandardError $driftStderr `
    -ArgumentList @('-NoProfile','-ExecutionPolicy','Bypass','-File',$Harness,
        '-PythonExe',$PythonExe,'-WorkRoot',$driftRoot,'-Z4EWorkRoot',$z4eCopy,
        '-SkipVisualQa')
$driftOutput = @()
if (Test-Path -LiteralPath $driftStdout) { $driftOutput += Get-Content -LiteralPath $driftStdout }
if (Test-Path -LiteralPath $driftStderr) { $driftOutput += Get-Content -LiteralPath $driftStderr }
if ($driftProcess.ExitCode -eq 0) { throw 'Tampered checksum manifest unexpectedly passed.' }
if (($driftOutput -join "`n") -notmatch 'checksum drift') {
    throw "Checksum drift failure was not explicit: $($driftOutput -join ' ')"
}

Write-Host 'PASS pbe_fagbi_scope_z4f_contract_smoke'
