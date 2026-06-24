param(
  [string]$ArtifactRoot = "test/tmp",
  [string]$PackageRoot = "docs/fagbi_guarded_p4_package",
  [switch]$AllowMissingOptional
)

$ErrorActionPreference = "Stop"

function Resolve-ExistingFile {
  param(
    [string]$Path,
    [bool]$Required = $true
  )
  if (Test-Path -LiteralPath $Path -PathType Leaf) {
    return (Resolve-Path -LiteralPath $Path).Path
  }
  if ($Required) {
    throw "Missing required artifact: $Path"
  }
  return $null
}

function Copy-Artifact {
  param(
    [string]$Source,
    [string]$Destination
  )
  $destDir = Split-Path -Parent $Destination
  New-Item -ItemType Directory -Force -Path $destDir | Out-Null
  Copy-Item -LiteralPath $Source -Destination $Destination -Force
  return (Resolve-Path -LiteralPath $Destination).Path
}

function Write-SimpleSvg {
  param(
    [array]$RepeatRows,
    [string]$Destination
  )
  $width = 900
  $height = 360
  $left = 80
  $bottom = 300
  $barWidth = 58
  $gap = 42
  $maxAudit = [Math]::Max(0.008, ($RepeatRows | ForEach-Object { [double]$_.audit_max } | Measure-Object -Maximum).Maximum)
  $maxRatio = [Math]::Max(1.0, ($RepeatRows | ForEach-Object {
        if ($_.median_p4_vs_degree1_matvec_ratio) { [double]$_.median_p4_vs_degree1_matvec_ratio } else { 0.0 }
      } | Measure-Object -Maximum).Maximum)

  $parts = New-Object System.Collections.Generic.List[string]
  $parts.Add("<svg xmlns=""http://www.w3.org/2000/svg"" width=""$width"" height=""$height"" viewBox=""0 0 $width $height"">")
  $parts.Add("<rect x=""0"" y=""0"" width=""$width"" height=""$height"" fill=""white""/>")
  $parts.Add("<text x=""$left"" y=""34"" font-family=""Arial"" font-size=""18"" font-weight=""700"">Guarded p4 repeatability supplement</text>")
  $parts.Add("<text x=""$left"" y=""57"" font-family=""Arial"" font-size=""12"" fill=""#555"">Audit RMS bars use the left scale; matvec ratio markers use the right scale.</text>")
  $parts.Add("<line x1=""$left"" y1=""78"" x2=""$left"" y2=""$bottom"" stroke=""#444""/>")
  $parts.Add("<line x1=""$left"" y1=""$bottom"" x2=""830"" y2=""$bottom"" stroke=""#444""/>")
  $parts.Add("<text x=""18"" y=""84"" font-family=""Arial"" font-size=""11"" fill=""#333"">audit</text>")
  $parts.Add("<text x=""812"" y=""84"" font-family=""Arial"" font-size=""11"" fill=""#333"">ratio</text>")
  $parts.Add("<line x1=""$left"" y1=""86"" x2=""830"" y2=""86"" stroke=""#e6e6e6""/>")
  $parts.Add("<text x=""24"" y=""91"" font-family=""Arial"" font-size=""10"" fill=""#777"">$('{0:0.###E+0}' -f $maxAudit)</text>")
  $x = $left + 55
  foreach ($row in $RepeatRows) {
    $audit = [double]$row.audit_max
    $ratio = if ($row.median_p4_vs_degree1_matvec_ratio) { [double]$row.median_p4_vs_degree1_matvec_ratio } else { 0.0 }
    $barHeight = [Math]::Max(2, [Math]::Round(($audit / $maxAudit) * 190))
    $y = $bottom - $barHeight
    $ratioY = $bottom - [Math]::Round(($ratio / $maxRatio) * 190)
    $fixture = [System.Security.SecurityElement]::Escape($row.fixture.ToUpperInvariant())
    $parts.Add("<rect x=""$x"" y=""$y"" width=""$barWidth"" height=""$barHeight"" fill=""#2f6f9f""/>")
    $parts.Add("<circle cx=""$($x + [int]($barWidth / 2))"" cy=""$ratioY"" r=""6"" fill=""#b4432f""/>")
    $parts.Add("<text x=""$($x + [int]($barWidth / 2))"" y=""$($bottom + 20)"" text-anchor=""middle"" font-family=""Arial"" font-size=""12"">$fixture</text>")
    $parts.Add("<text x=""$($x + [int]($barWidth / 2))"" y=""$($y - 7)"" text-anchor=""middle"" font-family=""Arial"" font-size=""10"" fill=""#2f6f9f"">$('{0:0.##E+0}' -f $audit)</text>")
    $parts.Add("<text x=""$($x + [int]($barWidth / 2))"" y=""$($ratioY - 10)"" text-anchor=""middle"" font-family=""Arial"" font-size=""10"" fill=""#b4432f"">$('{0:0.###}' -f $ratio)</text>")
    $x += $barWidth + $gap
  }
  $parts.Add("<rect x=""600"" y=""37"" width=""14"" height=""10"" fill=""#2f6f9f""/><text x=""620"" y=""46"" font-family=""Arial"" font-size=""11"">max audit RMS</text>")
  $parts.Add("<circle cx=""607"" cy=""62"" r=""5"" fill=""#b4432f""/><text x=""620"" y=""66"" font-family=""Arial"" font-size=""11"">median p4/degree-one matvec ratio</text>")
  $parts.Add("</svg>")
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $Destination) | Out-Null
  Set-Content -LiteralPath $Destination -Value ($parts -join [Environment]::NewLine) -Encoding UTF8
}

$repoRoot = (Resolve-Path ".").Path
$artifactRootPath = (Resolve-Path -LiteralPath $ArtifactRoot).Path
$packageRootPath = Join-Path $repoRoot $PackageRoot
$z4p = Join-Path $artifactRootPath "fagbi_scope_z4p"
$z4q = Join-Path $artifactRootPath "fagbi_scope_z4q"

$required = @(
  @{ Source = Join-Path $z4p "titan_fagbi_guarded_p4_manuscript.docx"; Dest = "manuscript/titan_fagbi_guarded_p4_manuscript.docx" },
  @{ Source = Join-Path $z4p "render_qa/titan_fagbi_guarded_p4_manuscript.pdf"; Dest = "manuscript/titan_fagbi_guarded_p4_manuscript.pdf" },
  @{ Source = Join-Path $z4p "fagbi_scope_z4p_manuscript_source.md"; Dest = "manuscript/fagbi_scope_z4p_manuscript_source.md" },
  @{ Source = Join-Path $z4p "fagbi_scope_z4p_claim_trace.csv"; Dest = "provenance/fagbi_scope_z4p_claim_trace.csv" },
  @{ Source = Join-Path $z4p "fagbi_scope_z4p_docx_qa.csv"; Dest = "provenance/fagbi_scope_z4p_docx_qa.csv" },
  @{ Source = Join-Path $z4p "fagbi_scope_z4p_decision.csv"; Dest = "tables/fagbi_scope_z4p_decision.csv" },
  @{ Source = Join-Path $z4p "fagbi_scope_z4p_summary.dat"; Dest = "tables/fagbi_scope_z4p_summary.dat" },
  @{ Source = Join-Path $z4q "fagbi_scope_z4q_matrix.csv"; Dest = "tables/fagbi_scope_z4q_matrix.csv" },
  @{ Source = Join-Path $z4q "fagbi_scope_z4q_repeatability.csv"; Dest = "tables/fagbi_scope_z4q_repeatability.csv" },
  @{ Source = Join-Path $z4q "fagbi_scope_z4q_decision.csv"; Dest = "tables/fagbi_scope_z4q_decision.csv" },
  @{ Source = Join-Path $z4q "fagbi_scope_z4q_summary.dat"; Dest = "tables/fagbi_scope_z4q_summary.dat" },
  @{ Source = "test/pbe_inputs_public_z4q/public_fixture_expansion_manifest.csv"; Dest = "provenance/public_fixture_expansion_manifest.csv" },
  @{ Source = "test/pbe_fagbi_scope_z4p_z4o_snapshot_matrix.csv"; Dest = "tables/fagbi_scope_z4p_z4o_snapshot_matrix.csv" },
  @{ Source = "test/pbe_fagbi_scope_z4p_z4o_snapshot_repeatability.csv"; Dest = "tables/fagbi_scope_z4p_z4o_snapshot_repeatability.csv" },
  @{ Source = "test/pbe_fagbi_scope_z4p_z4o_snapshot_decision.csv"; Dest = "tables/fagbi_scope_z4p_z4o_snapshot_decision.csv" },
  @{ Source = "TITAN_FMM_M2L_SCOPE_Z4O_P4_REPEATABILITY_REPORT.txt"; Dest = "reports/TITAN_FMM_M2L_SCOPE_Z4O_P4_REPEATABILITY_REPORT.txt" },
  @{ Source = "TITAN_FMM_M2L_SCOPE_Z4P_MANUSCRIPT_GUARDED_P4_REPORT.txt"; Dest = "reports/TITAN_FMM_M2L_SCOPE_Z4P_MANUSCRIPT_GUARDED_P4_REPORT.txt" },
  @{ Source = "TITAN_FMM_M2L_SCOPE_Z4Q_PUBLIC_FIXTURE_EXPANSION_REPORT.txt"; Dest = "reports/TITAN_FMM_M2L_SCOPE_Z4Q_PUBLIC_FIXTURE_EXPANSION_REPORT.txt" }
)

$copied = New-Object System.Collections.Generic.List[object]
foreach ($item in $required) {
  $source = Resolve-ExistingFile -Path $item.Source -Required $true
  $dest = Join-Path $packageRootPath $item.Dest
  $copiedPath = Copy-Artifact -Source $source -Destination $dest
  $copied.Add([pscustomobject]@{ source = $source; package_path = (Resolve-Path -LiteralPath $copiedPath).Path.Substring($repoRoot.Length + 1) })
}

$z4qPlot = Join-Path $z4q "fagbi_scope_z4q_repeatability_plot.svg"
$plotDest = Join-Path $packageRootPath "figures/fagbi_scope_z4q_repeatability_plot.svg"
if ((Test-Path -LiteralPath $z4qPlot -PathType Leaf) -and ((Get-Item -LiteralPath $z4qPlot).Length -gt 0)) {
  $copiedPath = Copy-Artifact -Source $z4qPlot -Destination $plotDest
} else {
  $repeatRows = Import-Csv (Join-Path $z4q "fagbi_scope_z4q_repeatability.csv")
  $matrixRows = Import-Csv (Join-Path $z4q "fagbi_scope_z4q_matrix.csv")
  foreach ($row in $repeatRows) {
    $p4 = $matrixRows | Where-Object { $_.fixture -eq $row.fixture -and $_.variant -eq "p4_guarded" }
    $ratios = @($p4 | ForEach-Object { [double]$_.p4_vs_degree1_matvec_ratio } | Sort-Object)
    if ($ratios.Count -gt 0) {
      $row | Add-Member -NotePropertyName median_p4_vs_degree1_matvec_ratio -NotePropertyValue $ratios[[int][Math]::Floor(($ratios.Count - 1) / 2)] -Force
    }
  }
  Write-SimpleSvg -RepeatRows $repeatRows -Destination $plotDest
  $copiedPath = (Resolve-Path -LiteralPath $plotDest).Path
}
$copied.Add([pscustomobject]@{ source = $z4qPlot; package_path = (Resolve-Path -LiteralPath $copiedPath).Path.Substring($repoRoot.Length + 1) })

$z4qDecision = Import-Csv (Join-Path $z4q "fagbi_scope_z4q_decision.csv") | Select-Object -First 1
$z4qRepeat = Import-Csv (Join-Path $z4q "fagbi_scope_z4q_repeatability.csv")
$decision = if ($z4qDecision.decision -eq "z4q_expansion_pass" -and (($z4qRepeat | Where-Object { $_.repeatability_pass -ne "True" }).Count -eq 0)) {
  "z4s_paper_package_ready"
} else {
  "z4s_evidence_gap"
}

$decisionPath = Join-Path $packageRootPath "tables/fagbi_scope_z4s_decision.csv"
@([pscustomobject]@{
  decision = $decision
  package_root = $PackageRoot
  z4q_decision = $z4qDecision.decision
  packaged_files = ""
  grid_default = "unchanged"
  tabi_treecode_default = "degree_one"
  guarded_p4_policy = "advanced_explicit"
}) | Export-Csv -LiteralPath $decisionPath -NoTypeInformation

$readme = @"
# TITAN Guarded p4 FAGBI Paper Package

This package is a tracked, paper-facing bundle for the guarded p4 TABI treecode evidence.

## Scope

- Production PBE default remains `-pbe=grid`.
- TABI treecode default remains the degree-one engine.
- Guarded p4 remains an advanced explicit TABI treecode path.
- Z4Q supplement rows cover full public proteins 1UBQ, 2CI2, and 1VII.

## Contents

- `manuscript/`: refreshed DOCX, PDF render, and source Markdown.
- `tables/`: Z4O snapshot tables, Z4P decisions, and full Z4Q supplement matrix/repeatability.
- `figures/`: supplement repeatability figure.
- `provenance/`: claim trace, DOCX QA, public fixture manifest, and package checksums.
- `reports/`: root Z4O/Z4P/Z4Q evidence reports copied at packaging time.

## Decision

Package decision: `$decision`.

The package was generated from artifact root `$artifactRootPath` and keeps generated scratch files out of `test/tmp`.
"@
$readmePath = Join-Path $packageRootPath "README_PAPER_PACKAGE.md"
Set-Content -LiteralPath $readmePath -Value $readme -Encoding UTF8

$hashPath = Join-Path $packageRootPath "provenance/package_checksums.csv"
$hashRows = Get-ChildItem -LiteralPath $packageRootPath -Recurse -File |
  Where-Object { $_.FullName -ne $hashPath } |
  Sort-Object FullName |
  ForEach-Object {
    $hash = Get-FileHash -Algorithm SHA256 -LiteralPath $_.FullName
    [pscustomobject]@{
      relative_path = $_.FullName.Substring($repoRoot.Length + 1)
      bytes = $_.Length
      sha256 = $hash.Hash.ToLowerInvariant()
    }
  }
$hashRows | Export-Csv -LiteralPath $hashPath -NoTypeInformation

$decisionRows = Import-Csv -LiteralPath $decisionPath
$decisionRows[0].packaged_files = ($hashRows.Count + 1).ToString()
$decisionRows | Export-Csv -LiteralPath $decisionPath -NoTypeInformation

$rootReport = "TITAN_FMM_M2L_SCOPE_Z4S_PAPER_PACKAGE_REPORT.txt"
$report = @"
Scope Z4S: Paper package closeout
decision: $decision
package_root: $PackageRoot
artifact_root: $artifactRootPath
z4q_decision: $($z4qDecision.decision)
supplement_fixtures: $((($z4qRepeat | ForEach-Object { $_.fixture }) -join ","))
default_policy: grid remains production default; TABI treecode remains degree-one default; guarded p4 is explicit advanced path
checksums: $($hashPath.Substring($repoRoot.Length + 1))
"@
Set-Content -LiteralPath $rootReport -Value $report -Encoding UTF8

Write-Host "Z4S decision: $decision"
Write-Host "Package root: $PackageRoot"
Write-Host "Packaged files: $($hashRows.Count + 1)"
