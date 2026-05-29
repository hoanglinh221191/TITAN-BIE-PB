param(
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$OutRoot = (Join-Path $PSScriptRoot 'tmp\pbe_clean_room_audit')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function New-CleanDir([string]$Path) {
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Get-RelativePath([string]$Base, [string]$Path) {
    $baseFull = [System.IO.Path]::GetFullPath($Base)
    if (-not $baseFull.EndsWith([System.IO.Path]::DirectorySeparatorChar)) {
        $baseFull += [System.IO.Path]::DirectorySeparatorChar
    }
    $pathFull = [System.IO.Path]::GetFullPath($Path)
    $baseUri = [System.Uri]::new($baseFull)
    $pathUri = [System.Uri]::new($pathFull)
    return [System.Uri]::UnescapeDataString(
        $baseUri.MakeRelativeUri($pathUri).ToString()
    ).Replace('/', [System.IO.Path]::DirectorySeparatorChar)
}

function Read-Lines([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return @() }
    return @(Get-Content -LiteralPath $Path)
}

New-CleanDir $OutRoot

$productionFiles = New-Object System.Collections.Generic.List[string]
$explicitProduction = @(
    'src\titan_pbe.h',
    'src\arguments.cpp',
    'src\arguments.h',
    'src\parse_arguments.cpp',
    'src\impl_engine.cpp',
    'Makefile'
)
foreach ($name in $explicitProduction) {
    $path = Join-Path $RepoRoot $name
    if (Test-Path -LiteralPath $path) { $productionFiles.Add($path) }
}
$pbeDir = Join-Path $RepoRoot 'src\tabi'
if (Test-Path -LiteralPath $pbeDir) {
    Get-ChildItem -LiteralPath $pbeDir -Recurse -File -Include *.cpp,*.h |
        ForEach-Object { $productionFiles.Add($_.FullName) }
}

$forbidden = @(
    @{ name = 'codex_tmp_pyDelPhi'; pattern = '(?i)\.codex_tmp[\\/]+pyDelPhi' },
    @{ name = 'pydelphi_env'; pattern = '(?i)pydelphi_env' },
    @{ name = 'pyDelPhi_runtime'; pattern = '(?i)(PyDelPhiRoot|pyDelPhiRoot|pyDelPhi[\\/]|from\s+.*pyDelPhi|import\s+.*pyDelPhi|#include\s+[<"].*pyDelPhi)' },
    @{ name = 'python_oracle_runtime'; pattern = '(?i)\b(PythonExe|python\.exe|python3?)(\b|[\\/:])' },
    @{ name = 'delphi_executable'; pattern = '(?i)\b(delphi\.exe|DelPhi\.exe)\b|[\\/](delphi|DelPhi)(\.exe)?\b' },
    @{ name = 'tabi_external_runtime'; pattern = '(?i)(TABI-PB|NanoShaper|nanoshaper)' },
    @{ name = 'subprocess_runtime'; pattern = '(?i)\b(system|popen|CreateProcess|ShellExecute|std::system)\s*\(' }
)

$allowedDelphiLabel = '(?i)\bdelphi[-_a-z0-9]*\b'
$oracleCsv = '(?i)pbe_oracle_[A-Za-z0-9_]*\.(csv|tsv)'

$failures = New-Object System.Collections.Generic.List[object]
$labelInventory = New-Object System.Collections.Generic.List[object]
$oracleInventory = New-Object System.Collections.Generic.List[object]

foreach ($file in $productionFiles) {
    $rel = Get-RelativePath $RepoRoot $file
    $lineNo = 0
    foreach ($line in Read-Lines $file) {
        ++$lineNo
        foreach ($rule in $forbidden) {
            if ($line -match $rule.pattern) {
                $failures.Add([PSCustomObject]@{
                    file = $rel
                    line = $lineNo
                    rule = $rule.name
                    text = $line.Trim()
                })
            }
        }
        if ($line -match $allowedDelphiLabel) {
            $labelInventory.Add([PSCustomObject]@{
                file = $rel
                line = $lineNo
                text = $line.Trim()
            })
        }
        if ($line -match $oracleCsv) {
            $oracleInventory.Add([PSCustomObject]@{
                file = $rel
                line = $lineNo
                text = $line.Trim()
                classification = if ($line -match '(?i)test_|parity|oracle|hidden') { 'hidden_test_oracle_reference' } else { 'review_required' }
            })
        }
    }
}

foreach ($row in @($oracleInventory | Where-Object { $_.classification -eq 'review_required' })) {
    $failures.Add([PSCustomObject]@{
        file = $row.file
        line = $row.line
        rule = 'oracle_csv_not_hidden_test_path'
        text = $row.text
    })
}

$testInventory = New-Object System.Collections.Generic.List[object]
$testDir = Join-Path $RepoRoot 'test'
if (Test-Path -LiteralPath $testDir) {
    Get-ChildItem -LiteralPath $testDir -File -Include *.ps1,*.py,*.md,*.csv -Recurse |
        Where-Object {
            $relPath = Get-RelativePath $RepoRoot $_.FullName
            $relPath -notmatch '^test[\\/](tmp|bin|\.codex_tmp|benchmark_run_|graph_smoke_|llps_|codex_)'
        } |
        ForEach-Object {
            $file = $_.FullName
            $rel = Get-RelativePath $RepoRoot $file
            $lineNo = 0
            foreach ($line in Read-Lines $file) {
                ++$lineNo
                if ($line -match '(?i)pyDelPhi|pydelphi|\.codex_tmp|pbe_oracle_|DelPhi') {
                    $testInventory.Add([PSCustomObject]@{
                        file = $rel
                        line = $lineNo
                        text = $line.Trim()
                    })
                }
            }
        }
}

$failureCsv = Join-Path $OutRoot 'production_forbidden_references.csv'
$labelCsv = Join-Path $OutRoot 'production_delphi_label_inventory.csv'
$oracleCsvPath = Join-Path $OutRoot 'production_oracle_reference_inventory.csv'
$testCsv = Join-Path $OutRoot 'test_oracle_usage_inventory.csv'
$summaryCsv = Join-Path $OutRoot 'pbe_clean_room_audit_summary.csv'

$failures | Export-Csv -LiteralPath $failureCsv -NoTypeInformation -Encoding ASCII
$labelInventory | Export-Csv -LiteralPath $labelCsv -NoTypeInformation -Encoding ASCII
$oracleInventory | Export-Csv -LiteralPath $oracleCsvPath -NoTypeInformation -Encoding ASCII
$testInventory | Export-Csv -LiteralPath $testCsv -NoTypeInformation -Encoding ASCII

$summary = [PSCustomObject]@{
    status = if ($failures.Count -eq 0) { 'pass' } else { 'fail' }
    production_files_scanned = $productionFiles.Count
    forbidden_references = $failures.Count
    delphi_label_inventory_rows = $labelInventory.Count
    production_oracle_reference_rows = $oracleInventory.Count
    test_oracle_usage_rows = $testInventory.Count
    output_root = $OutRoot
}
$summary | Export-Csv -LiteralPath $summaryCsv -NoTypeInformation -Encoding ASCII

if ($failures.Count -gt 0) {
    Write-Error "PBE clean-room audit FAILED: $($failures.Count) forbidden production reference(s). See $failureCsv"
}

Write-Host "PBE clean-room audit PASS: $summaryCsv"
Write-Host "PBE clean-room inventories: $labelCsv ; $oracleCsvPath ; $testCsv"
