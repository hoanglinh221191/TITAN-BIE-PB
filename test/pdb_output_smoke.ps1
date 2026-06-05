param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pdb_output_smoke'),
    [string]$LibRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) 'lib_ff')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$env:TITAN_NO_PAUSE = '1'

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

function Write-ArgFile([string]$Path, [string]$Pdb, [string]$Index, [string]$OutputFlag) {
    $lib = Convert-ToTitanPath ((Resolve-Path -LiteralPath $LibRoot).Path.TrimEnd('\') + '\')
    @(
        "-f=$(Convert-ToTitanPath (Split-Path -Leaf $Pdb))",
        "-n=$(Convert-ToTitanPath $Index)",
        '-g=3-4',
        '-1pdb=true',
        '-nocal=true',
        "-path-a=$lib",
        "-path-n=$lib",
        '-nt=1',
        $OutputFlag
    ) | Set-Content -LiteralPath $Path -Encoding ASCII
}

function Invoke-Titan([string]$CaseDir, [string]$ArgPath) {
    Push-Location $CaseDir
    try {
        & $TitanExe -a (Split-Path -Leaf $ArgPath) | Out-Null
        if ($LASTEXITCODE -ne 0) {
            throw "TITAN exited with code $LASTEXITCODE"
        }
    }
    finally {
        Pop-Location
    }
}

function Assert-PdbContent([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) {
        throw "Missing PDB output: $Path"
    }
    $modelCount = @(Select-String -Path $Path -Pattern '^MODEL').Count
    $atomCount = @(Select-String -Path $Path -Pattern '^ATOM').Count
    $terCount = @(Select-String -Path $Path -Pattern '^TER$').Count
    $endmdlCount = @(Select-String -Path $Path -Pattern '^ENDMDL$').Count
    if ($modelCount -lt 1) {
        throw "$Path has no MODEL records"
    }
    if ($atomCount -le 0) {
        throw "$Path has no ATOM records"
    }
    if ($terCount -lt 2) {
        throw "$Path does not contain separate mol1/mol2 TER records"
    }
    if ($endmdlCount -lt 1) {
        throw "$Path has no ENDMDL records"
    }
    $firstAtom = (Select-String -Path $Path -Pattern '^ATOM' | Select-Object -First 1).Line
    $x = [double]($firstAtom.Substring(30, 8).Trim())
    $occ = [double]($firstAtom.Substring(54, 6).Trim())
    $bf = [double]($firstAtom.Substring(60, 6).Trim())
    if ([double]::IsNaN($x) -or [double]::IsInfinity($x) -or
        [math]::Abs($occ - 1.00) -gt 0.001 -or
        [math]::Abs($bf - 0.00) -gt 0.001) {
        throw "$Path has unexpected coordinate/occupancy/B-factor formatting"
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
if (-not (Test-Path -LiteralPath $LibRoot)) {
    throw "Missing TITAN force-field library directory: $LibRoot"
}
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$sourcePdb = Join-Path $repo 'example\2\md1_short.pdb'
$index = Join-Path $repo 'example\2\index.ndx'
if (-not (Test-Path -LiteralPath $sourcePdb)) {
    throw "Missing fixture PDB: $sourcePdb"
}
if (-not (Test-Path -LiteralPath $index)) {
    throw "Missing fixture index: $index"
}

$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
$WorkRoot = [System.IO.Path]::GetFullPath($WorkRoot)
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$appendDir = Join-Path $WorkRoot 'append'
New-Item -ItemType Directory -Path $appendDir -Force | Out-Null
$appendPdb = Join-Path $appendDir 'fixture.pdb'
$frameOne = [System.IO.File]::ReadAllText($sourcePdb)
$frameTwo = $frameOne.Replace('MODEL        1', 'MODEL        2').Replace(
    't= 100000.00000', 't= 100001.00000')
[System.IO.File]::WriteAllText(
    $appendPdb, $frameOne + "`r`n" + $frameTwo, [System.Text.Encoding]::ASCII)
$appendArg = Join-Path $appendDir 'arg.in'
Write-ArgFile $appendArg $appendPdb $index '-opdb=out.pdb'
Invoke-Titan $appendDir $appendArg
Assert-PdbContent (Join-Path $appendDir 'out.pdb')

$sepDir = Join-Path $WorkRoot 'separate'
New-Item -ItemType Directory -Path $sepDir -Force | Out-Null
$sepPdb = Join-Path $sepDir 'fixture.pdb'
[System.IO.File]::WriteAllText(
    $sepPdb, $frameOne + "`r`n" + $frameTwo, [System.Text.Encoding]::ASCII)
$sepArg = Join-Path $sepDir 'arg.in'
Write-ArgFile $sepArg $sepPdb $index '-opdbsep=out.pdb'
Invoke-Titan $sepDir $sepArg
$sepOutputs = @(Get-ChildItem -LiteralPath $sepDir -Recurse -Filter 'out_*.pdb')
if ($sepOutputs.Count -lt 1) {
    throw 'Expected at least one separated PDB output'
}
foreach ($file in $sepOutputs) {
    Assert-PdbContent $file.FullName
}

Write-Host 'PDB output smoke PASS'
