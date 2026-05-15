param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_near_block_rank'),
    [string]$InputCsv = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI near block rank root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Convert-ToInt64([object]$Value) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return [int64]0 }
    return [int64]::Parse($text, $InvariantCulture)
}

function Convert-ToDouble([object]$Value) {
    $text = [string]$Value
    if ([string]::IsNullOrWhiteSpace($text) -or $text -eq 'NA') { return 0.0 }
    return [double]::Parse($text, $InvariantCulture)
}

function Get-Value([object]$Row, [string]$Name, [string]$Default = '0') {
    if ($Row.PSObject.Properties.Name -contains $Name) { return $Row.$Name }
    return $Default
}

function Get-Priority([object]$Row) {
    $block = [string]$Row.block
    $pair = [string]$Row.pair_type
    $orientation = [string]$Row.orientation
    if ($block -eq 'K2' -and ($pair -eq 'vertex' -or $pair -eq 'edge') -and
        $orientation -eq 'same') { return 'P0_k2_source_normal_same_edge_vertex' }
    if ($block -eq 'K1' -and
        ($pair -eq 'vertex' -or $pair -eq 'edge') -and
        $orientation -eq 'same') { return 'P0_k1_k2_same_edge_vertex' }
    if ($block -eq 'K4' -and ($pair -eq 'vertex' -or $pair -eq 'edge')) {
        return 'P1_k4_edge_vertex'
    }
    if ($block -eq 'K3' -and ($pair -eq 'vertex' -or $pair -eq 'edge')) {
        return 'P2_k3_strict_guard'
    }
    if ($pair -eq 'self') { return 'self_reference' }
    return 'background'
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}

New-CleanDir $WorkRoot
if ([string]::IsNullOrWhiteSpace($InputCsv)) {
    $profileRoot = Join-Path $WorkRoot 'profile_source'
    & (Join-Path $PSScriptRoot 'pbe_tabi_profile_smoke.ps1') -TitanExe $TitanExe -WorkRoot $profileRoot
    if ($LASTEXITCODE -ne 0) { throw 'TABI profile smoke failed while generating near block audit input.' }
    $found = @(Get-ChildItem -LiteralPath $profileRoot -Recurse -Filter 'pbe_tabi_near_block_audit.csv' |
        Sort-Object LastWriteTime -Descending)
    if ($found.Count -eq 0) {
        throw "Could not find generated pbe_tabi_near_block_audit.csv under $profileRoot"
    }
    $InputCsv = $found[0].FullName
}

if (-not (Test-Path -LiteralPath $InputCsv)) {
    throw "Missing near block audit CSV: $InputCsv"
}

$inputRows = @(Import-Csv -LiteralPath $InputCsv)
if ($inputRows.Count -eq 0) { throw "Empty near block audit CSV: $InputCsv" }

$ranked = New-Object System.Collections.Generic.List[object]
foreach ($row in $inputRows) {
    $accepted = Convert-ToInt64 (Get-Value $row 'accepted')
    $rejected = Convert-ToInt64 (Get-Value $row 'rejected')
    $worse = Convert-ToInt64 (Get-Value $row 'candidate_worse_than_centroid')
    $refNonconv = Convert-ToInt64 (Get-Value $row 'reference_nonconverged')
    $orient = Convert-ToInt64 (Get-Value $row 'orientation_sign_mismatch')
    $strict = Convert-ToInt64 (Get-Value $row 'strict_ref_reject')
    $maxCand = Convert-ToDouble (Get-Value $row 'max_candidate_ref_rel')
    $maxCent = Convert-ToDouble (Get-Value $row 'max_centroid_ref_rel')
    $meanCand = Convert-ToDouble (Get-Value $row 'mean_candidate_ref_rel')
    $k2SignOk = Convert-ToInt64 (Get-Value $row 'k2_laplace_sign_ok')
    $k2SmoothRel = Convert-ToDouble (Get-Value $row 'k2_yukawa_smooth_ref_rel')
    $k2OrientPairRel = Convert-ToDouble (Get-Value $row 'k2_orientation_pair_ref_rel')
    $priority = Get-Priority $row
    $k2DiagnosticScore = 0.0
    if ([string]$row.block -eq 'K2') {
        $k2DiagnosticScore = (1000.0 * $k2OrientPairRel) + (100.0 * $k2SmoothRel)
    }
    $rankScore = (1000000.0 * $worse) + (100000.0 * $orient) +
        (10000.0 * $refNonconv) + $k2DiagnosticScore +
        (100.0 * $strict) + $maxCand
    $ranked.Add([PSCustomObject]@{
        rank = 0
        priority = $priority
        block = $row.block
        pair_type = $row.pair_type
        orientation = $row.orientation
        accepted = $accepted
        rejected = $rejected
        candidate_worse_than_centroid = $worse
        reference_nonconverged = $refNonconv
        orientation_sign_mismatch = $orient
        strict_ref_reject = $strict
        max_candidate_ref_rel = $maxCand.ToString('G17', $InvariantCulture)
        mean_candidate_ref_rel = $meanCand.ToString('G17', $InvariantCulture)
        max_centroid_ref_rel = $maxCent.ToString('G17', $InvariantCulture)
        k2_laplace_sign_ok = $k2SignOk
        k2_yukawa_smooth_ref_rel = $k2SmoothRel.ToString('G17', $InvariantCulture)
        k2_orientation_pair_ref_rel = $k2OrientPairRel.ToString('G17', $InvariantCulture)
        rank_score = $rankScore.ToString('G17', $InvariantCulture)
    })
}

$ordered = @($ranked.ToArray() | Sort-Object `
    @{ Expression = { [double]$_.rank_score }; Descending = $true },
    @{ Expression = { [double]$_.k2_orientation_pair_ref_rel }; Descending = $true },
    @{ Expression = { [int64]$_.candidate_worse_than_centroid }; Descending = $true },
    @{ Expression = { [double]$_.max_candidate_ref_rel }; Descending = $true },
    block, pair_type, orientation)

for ($i = 0; $i -lt $ordered.Count; ++$i) {
    $ordered[$i].rank = $i + 1
}

$outCsv = Join-Path $WorkRoot 'pbe_tabi_near_block_rank.csv'
$ordered | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII

$top = $ordered[0]
if ((Convert-ToInt64 $top.rejected) -eq 0 -and
    (Convert-ToInt64 $top.candidate_worse_than_centroid) -eq 0 -and
    (Convert-ToInt64 $top.strict_ref_reject) -eq 0) {
    throw "Near block rank did not find any guarded bucket to prioritize."
}

Write-Host ("PBE TABI near block rank PASS: top={0}/{1}/{2} priority={3} out={4}" -f `
    $top.block, $top.pair_type, $top.orientation, $top.priority, $outCsv)
