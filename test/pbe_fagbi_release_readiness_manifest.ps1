param(
    [string]$OutputRoot = "docs\fagbi_guarded_p4_package",
    [string]$PublicRepositoryUrl = "https://github.com/hoanglinh221191/TITAN-fagbi-submission",
    [string]$ReleaseTag = "fagbi-submission-v0.1",
    [string]$ZenodoDoi = "pending"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Get-RepoRoot {
    return (Resolve-Path -LiteralPath ".").Path
}

function Get-RelativePath {
    param([string]$Path)
    $root = Get-RepoRoot
    $full = (Resolve-Path -LiteralPath $Path).Path
    if ($full.StartsWith($root, [System.StringComparison]::OrdinalIgnoreCase)) {
        return Convert-ToRepoSlash -Path $full.Substring($root.Length).TrimStart("\", "/")
    }
    return Convert-ToRepoSlash -Path $full
}

function Convert-ToRepoSlash {
    param([string]$Path)
    return ($Path -replace "\\", "/")
}

function New-CleanDirectory {
    param([string]$Path)
    $root = Get-RepoRoot
    if (Test-Path -LiteralPath $Path) {
        $resolved = (Resolve-Path -LiteralPath $Path).Path
        if (-not $resolved.StartsWith($root, [System.StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing to remove a path outside the repository: $resolved"
        }
        Remove-Item -LiteralPath $resolved -Recurse -Force
    }
    New-Item -ItemType Directory -Force -Path $Path | Out-Null
}

function Copy-IfPresent {
    param(
        [string]$Source,
        [string]$Destination
    )
    if (-not (Test-Path -LiteralPath $Source -PathType Leaf)) {
        return $false
    }
    $parent = Split-Path -Parent $Destination
    New-Item -ItemType Directory -Force -Path $parent | Out-Null
    Copy-Item -LiteralPath $Source -Destination $Destination -Force
    return $true
}

function Get-FileSha256Shared {
    param([string]$Path)
    $resolved = (Resolve-Path -LiteralPath $Path).Path
    try {
        $stream = [System.IO.File]::Open(
            $resolved,
            [System.IO.FileMode]::Open,
            [System.IO.FileAccess]::Read,
            [System.IO.FileShare]::ReadWrite
        )
        try {
            $sha = [System.Security.Cryptography.SHA256]::Create()
            try {
                $bytes = $sha.ComputeHash($stream)
                return ([System.BitConverter]::ToString($bytes) -replace "-", "").ToLowerInvariant()
            } finally {
                $sha.Dispose()
            }
        } finally {
            $stream.Dispose()
        }
    } catch [System.IO.IOException] {
        return Get-GitBlobSha256 -Path $Path
    }
}

function Get-GitBlobSha256 {
    param([string]$Path)
    $relative = Get-RelativePath -Path $Path
    $ls = (& git ls-files -s -- $relative)
    if (-not $ls) {
        throw "Cannot read locked file and no Git blob is available: $relative"
    }
    $parts = $ls -split "\s+"
    if ($parts.Count -lt 2) {
        throw "Unexpected git ls-files output for ${relative}: $ls"
    }
    $blob = $parts[1]
    $psi = [System.Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = "git"
    $psi.ArgumentList.Add("cat-file")
    $psi.ArgumentList.Add("blob")
    $psi.ArgumentList.Add($blob)
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $proc = [System.Diagnostics.Process]::Start($psi)
    try {
        $sha = [System.Security.Cryptography.SHA256]::Create()
        try {
            $bytes = $sha.ComputeHash($proc.StandardOutput.BaseStream)
            $proc.WaitForExit()
            if ($proc.ExitCode -ne 0) {
                $err = $proc.StandardError.ReadToEnd()
                throw "git cat-file failed for ${relative}: $err"
            }
            return ([System.BitConverter]::ToString($bytes) -replace "-", "").ToLowerInvariant()
        } finally {
            $sha.Dispose()
        }
    } finally {
        $proc.Dispose()
    }
}

function Add-ManifestEntry {
    param(
        [System.Collections.Generic.List[object]]$Entries,
        [string]$Path,
        [string]$Category,
        [string]$ReleaseRole,
        [string]$ManuscriptTarget = "",
        [string]$SourceCommand = ""
    )
    $relative = Convert-ToRepoSlash -Path $Path
    if (Test-Path -LiteralPath $Path -PathType Leaf) {
        $item = Get-Item -LiteralPath $Path
        $hash = Get-FileSha256Shared -Path $Path
        $Entries.Add([pscustomobject]@{
            category = $Category
            manuscript_target = $ManuscriptTarget
            release_role = $ReleaseRole
            relative_path = Get-RelativePath -Path $Path
            status = "present"
            bytes = [int64]$item.Length
            sha256 = $hash
            last_write_utc = $item.LastWriteTimeUtc.ToString("o")
            source_command = $SourceCommand
        })
    } else {
        $Entries.Add([pscustomobject]@{
            category = $Category
            manuscript_target = $ManuscriptTarget
            release_role = $ReleaseRole
            relative_path = $relative
            status = "missing"
            bytes = ""
            sha256 = ""
            last_write_utc = ""
            source_command = $SourceCommand
        })
    }
}

function Add-GlobEntries {
    param(
        [System.Collections.Generic.List[object]]$Entries,
        [string]$BasePath,
        [string[]]$Include,
        [string]$Category,
        [string]$ReleaseRole,
        [string]$ManuscriptTarget = ""
    )
    if (-not (Test-Path -LiteralPath $BasePath -PathType Container)) {
        return
    }
    Get-ChildItem -LiteralPath $BasePath -Recurse -File |
        Where-Object {
            foreach ($pattern in $Include) {
                if ($_.Name -like $pattern) { return $true }
            }
            return $false
        } |
        Sort-Object FullName |
        ForEach-Object {
            Add-ManifestEntry -Entries $Entries -Path $_.FullName -Category $Category -ReleaseRole $ReleaseRole -ManuscriptTarget $ManuscriptTarget
        }
}

function Export-Objects {
    param(
        [object[]]$Rows,
        [string]$Path
    )
    $parent = Split-Path -Parent $Path
    New-Item -ItemType Directory -Force -Path $parent | Out-Null
    $Rows | Export-Csv -LiteralPath $Path -NoTypeInformation
}

$outRoot = New-Item -ItemType Directory -Force -Path $OutputRoot
$payloadRoot = Join-Path $outRoot.FullName "release_payload"
$referenceRoot = Join-Path $payloadRoot "reference_outputs"
$apbsRoot = Join-Path $payloadRoot "apbs_context"
New-CleanDirectory -Path $payloadRoot
New-Item -ItemType Directory -Force -Path $referenceRoot, $apbsRoot | Out-Null

$entries = [System.Collections.Generic.List[object]]::new()
$traceRows = [System.Collections.Generic.List[object]]::new()
$copiedOutputs = [System.Collections.Generic.List[object]]::new()

$sourceMap = @(
    @{
        target = "Table 3"
        title = "Full-protein algorithmic ladder"
        source = "test\tmp\fagbi_scope_z4v\fagbi_scope_z4v_matrix.csv"
        dest = "tables\table_3_full_protein_algorithmic_ladder.csv"
        script = "test\pbe_fagbi_scope_z4v_performance_ladder_matrix.ps1"
    },
    @{
        target = "Table 4"
        title = "Guarded local FMM/M2L p-order and resource rows"
        source = "test\tmp\fagbi_scope_z4xd\fagbi_scope_z4xd_matrix.csv"
        dest = "tables\table_4_guarded_p_order_matrix.csv"
        script = "test\pbe_fagbi_scope_z4xd_p_order_tiled_cache_gate.ps1"
    },
    @{
        target = "Table 5"
        title = "Guarded p-order summary and coverage"
        source = "test\tmp\fagbi_scope_z4xd\fagbi_scope_z4xd_order_summary.csv"
        dest = "tables\table_5_guarded_p_order_summary.csv"
        script = "test\pbe_fagbi_scope_z4xd_p_order_tiled_cache_gate.ps1"
    },
    @{
        target = "Table 6"
        title = "ILU neighbor-cap sensitivity"
        source = "test\tmp\fagbi_scope_z4xf\fagbi_scope_z4xf_matrix.csv"
        dest = "tables\table_6_ilu_neighbor_cap_sweep.csv"
        script = "test\pbe_fagbi_scope_z4xf_ilu_neighbor_cap_sweep.ps1"
    },
    @{
        target = "Supplementary Table S1"
        title = "Supplementary public-protein rows"
        source = "docs\fagbi_guarded_p4_package\tables\fagbi_scope_z4q_matrix.csv"
        dest = "supplement\supplementary_table_s1_public_protein_rows.csv"
        script = "test\pbe_fagbi_scope_z4q_public_fixture_expansion_gate.ps1"
    },
    @{
        target = "Supplementary Table S2"
        title = "Aggregated theta, guard, and tail sensitivity"
        source = "test\tmp\fagbi_scope_z4w\fagbi_scope_z4w_matrix.csv"
        dest = "supplement\supplementary_table_s2_ablation_sensitivity.csv"
        script = "test\pbe_fagbi_scope_z4w_ablation_sensitivity_matrix.ps1"
    },
    @{
        target = "Supplementary Table S3"
        title = "Per-fixture sensitivity rows"
        source = "test\tmp\fagbi_scope_z4xd\fagbi_scope_z4xd_budget.csv"
        dest = "supplement\supplementary_table_s3_p_order_budget.csv"
        script = "test\pbe_fagbi_scope_z4xd_p_order_tiled_cache_gate.ps1"
    },
    @{
        target = "Supplementary Table S4"
        title = "ILU neighbor-cap sensitivity rows"
        source = "test\tmp\fagbi_scope_z4xf\fagbi_scope_z4xf_neighbor_cap_summary.csv"
        dest = "supplement\supplementary_table_s4_ilu_neighbor_cap_summary.csv"
        script = "test\pbe_fagbi_scope_z4xf_ilu_neighbor_cap_sweep.ps1"
    },
    @{
        target = "Supplementary Table S5"
        title = "Small dense direct calibration"
        source = "test\tmp\fagbi_scope_z4y_mini_dense\fagbi_scope_z4y_dense_oracle.csv"
        dest = "supplement\supplementary_table_s5_dense_small_oracle.csv"
        script = "test\pbe_fagbi_scope_z4y_mini_dense_oracle.ps1"
    },
    @{
        target = "Supplementary Table S6"
        title = "Mesh/surface-density sensitivity"
        source = "test\tmp\fagbi_scope_z4aa\fagbi_scope_z4aa_mesh_sensitivity.csv"
        dest = "supplement\supplementary_table_s6_mesh_density_sensitivity.csv"
        script = "test\pbe_fagbi_scope_z4aa_mesh_density_sensitivity.ps1"
    },
    @{
        target = "Supplementary Figure S1"
        title = "APBS grid-size context rows"
        source = "test\tmp\fagbi_scope_z4z_apbs\fagbi_scope_z4z_apbs_matrix.csv"
        dest = "supplement\supplementary_figure_s1_apbs_grid_matrix.csv"
        script = "test\pbe_fagbi_scope_z4z_apbs_grid_sweep.ps1"
    },
    @{
        target = "Supplementary Figure S1"
        title = "APBS grid-size stability settings"
        source = "test\tmp\fagbi_scope_z4z_apbs\fagbi_scope_z4z_apbs_stability.csv"
        dest = "supplement\supplementary_figure_s1_apbs_grid_stability.csv"
        script = "test\pbe_fagbi_scope_z4z_apbs_grid_sweep.ps1"
    }
)

foreach ($row in $sourceMap) {
    $destPath = Join-Path $referenceRoot $row.dest
    $copied = Copy-IfPresent -Source $row.source -Destination $destPath
    $archivedRel = Convert-ToRepoSlash -Path $destPath.Substring((Get-RepoRoot).Length).TrimStart("\", "/")
    $copiedOutputs.Add([pscustomobject]@{
        manuscript_target = $row.target
        title = $row.title
        source_path = Convert-ToRepoSlash -Path $row.source
        archived_path = $archivedRel
        status = $(if ($copied) { "copied" } else { "missing_source" })
        source_script = Convert-ToRepoSlash -Path $row.script
    })
    $traceRows.Add([pscustomobject]@{
        manuscript_target = $row.target
        title = $row.title
        archived_output = $archivedRel
        source_output = Convert-ToRepoSlash -Path $row.source
        reproduction_script = Convert-ToRepoSlash -Path $row.script
        notes = "Archived output copied from the current verified local run tree."
    })
}

$apbsMatrix = "test\tmp\fagbi_scope_z4z_apbs\fagbi_scope_z4z_apbs_matrix.csv"
if (Test-Path -LiteralPath $apbsMatrix) {
    $apbsSettings = Import-Csv -LiteralPath $apbsMatrix |
        Where-Object { $_.backend -eq "apbs" } |
        Select-Object case_id, fixture, protein_name, pdb_id, dime, fine_scale, cglen_x, cglen_y, cglen_z, fglen_x, fglen_y, fglen_z, apbs_exe, status, energy_kj_mol, wall_seconds, peak_rss_mb, run_dir
    Export-Objects -Rows $apbsSettings -Path (Join-Path $apbsRoot "apbs_grid_size_settings_supplementary_figure_s1.csv")
}

if (Test-Path -LiteralPath "test\tmp\fagbi_scope_z4z_apbs" -PathType Container) {
    Get-ChildItem -LiteralPath "test\tmp\fagbi_scope_z4z_apbs" -Directory |
        Where-Object { $_.Name -like "apbs_*" } |
        Sort-Object Name |
        ForEach-Object {
            $caseDest = Join-Path $apbsRoot $_.Name
            New-Item -ItemType Directory -Force -Path $caseDest | Out-Null
            foreach ($name in @("apbs.in", "stdout.log", "stderr.log", "io.mc")) {
                Copy-IfPresent -Source (Join-Path $_.FullName $name) -Destination (Join-Path $caseDest $name) | Out-Null
            }
            Get-ChildItem -LiteralPath $_.FullName -Filter "*.pqr" -File -ErrorAction SilentlyContinue |
                ForEach-Object {
                    Copy-IfPresent -Source $_.FullName -Destination (Join-Path $caseDest $_.Name) | Out-Null
                }
        }
}

$requiredFiles = @(
    @{ path = "LICENSE"; category = "license"; role = "MIT license for the public release"; target = "" },
    @{ path = "CITATION.cff"; category = "citation"; role = "citation metadata template for the archived release"; target = "" },
    @{ path = ".zenodo.json"; category = "archive_metadata"; role = "Zenodo metadata template for DOI minting"; target = "" },
    @{ path = "README.md"; category = "documentation"; role = "repository overview"; target = "" },
    @{ path = "docs\fagbi_guarded_p4_package\README_PAPER_PACKAGE.md"; category = "documentation"; role = "paper package overview and release-payload map"; target = "Code and Data Availability" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\titan_fagbi_jcp_apbs_bounded_references_clean.docx"; category = "manuscript"; role = "main manuscript submitted to JCP"; target = "main manuscript" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\titan_fagbi_jcp_si_apbs_bounded_cap96_mesh.docx"; category = "supplement"; role = "Supporting Information submitted to JCP"; target = "Supporting Information" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\jcp_cover_letter_titan_fagbi.md"; category = "submission_document"; role = "JCP cover letter draft"; target = "cover letter" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\fagbi_scope_z4t_references.bib"; category = "references"; role = "BibTeX references"; target = "references" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\fagbi_scope_z4t_references.csv"; category = "references"; role = "reference table"; target = "references" },
    @{ path = "test\pbe_inputs_public\public_fixture_manifest.csv"; category = "fixture_manifest"; role = "primary public-protein fixture provenance"; target = "input metadata" },
    @{ path = "test\pbe_inputs_public_z4q\public_fixture_expansion_manifest.csv"; category = "fixture_manifest"; role = "supplementary public-protein fixture provenance"; target = "input metadata" },
    @{ path = "test\pbe_fagbi_release_readiness_manifest.ps1"; category = "release_script"; role = "release manifest builder"; target = "release archive" }
)

foreach ($file in $requiredFiles) {
    Add-ManifestEntry -Entries $entries -Path $file.path -Category $file.category -ReleaseRole $file.role -ManuscriptTarget $file.target
}

foreach ($row in $sourceMap) {
    Add-ManifestEntry -Entries $entries -Path $row.script -Category "benchmark_script" -ReleaseRole "reported table or figure reproduction script" -ManuscriptTarget $row.target -SourceCommand ("powershell -NoProfile -ExecutionPolicy Bypass -File " + (Convert-ToRepoSlash -Path $row.script))
}

Add-GlobEntries -Entries $entries -BasePath "docs\fagbi_guarded_p4_package\submission_preview_synced\submission_figures" -Include @("*.svg", "*.png", "*.csv") -Category "figure" -ReleaseRole "submission figure files" -ManuscriptTarget "Figures 1-6 and Supplementary Figures S1-S2"
Add-GlobEntries -Entries $entries -BasePath $referenceRoot -Include @("*.csv", "*.dat", "*.svg", "*.png") -Category "reference_output" -ReleaseRole "archived reference output copied from benchmark run" -ManuscriptTarget "Tables 3-6 and Supplementary Tables S1-S6"
Add-GlobEntries -Entries $entries -BasePath $apbsRoot -Include @("*.csv", "*.in", "*.log", "*.mc", "*.pqr") -Category "apbs_context" -ReleaseRole "APBS input, log, and grid-size settings" -ManuscriptTarget "Supplementary Figure S1"
Add-GlobEntries -Entries $entries -BasePath "src" -Include @("*.cpp", "*.h", "*.hpp", "*.cxx") -Category "source_code" -ReleaseRole "TITAN source code"
Add-GlobEntries -Entries $entries -BasePath "test\pbe_inputs_public" -Include @("*.pdb", "*.pqr", "*.csv") -Category "input_metadata" -ReleaseRole "public protein PDB/PQR inputs and provenance"
Add-GlobEntries -Entries $entries -BasePath "test\pbe_inputs_public_z4q" -Include @("*.pdb", "*.pqr", "*.csv") -Category "input_metadata" -ReleaseRole "expanded public protein PDB/PQR inputs and provenance"

$csvPath = Join-Path $outRoot.FullName "release_artifact_manifest.csv"
$jsonPath = Join-Path $outRoot.FullName "release_artifact_manifest.json"
$summaryPath = Join-Path $outRoot.FullName "release_readiness_summary.dat"
$tracePath = Join-Path $payloadRoot "table_figure_traceability.csv"
$copiedPath = Join-Path $payloadRoot "reference_output_copy_manifest.csv"
$checksumPath = Join-Path $payloadRoot "checksums_sha256.csv"
$readmePath = Join-Path $payloadRoot "README_RELEASE_ARCHIVE.md"
$buildInstructionsPath = Join-Path $payloadRoot "BUILD_INSTRUCTIONS.md"
$rootReport = "TITAN_FAGBI_RELEASE_READINESS_REPORT.txt"

$gitCommit = "unavailable"
try {
    $gitCommit = (& git rev-parse HEAD).Trim()
} catch {
    $gitCommit = "unavailable"
}

$codeAvailability = "The TITAN source code corresponding to this manuscript has been prepared for release under the MIT License as the tagged release $ReleaseTag at $PublicRepositoryUrl. The persistent Zenodo DOI will be inserted after deposition of the tagged release and before journal submission. The archived release will include source code, build instructions, benchmark manifests, input metadata, scripts, figures, APBS input files, APBS logs, APBS grid-size settings used for Supplementary Figure S1, and reference outputs for Tables 3-6 and Supplementary Tables S1-S6. Public protein structures are identified by their RCSB PDB accession codes; derived surface metadata, run manifests, and tabulated outputs corresponding to the manuscript tables are included in the release manifest."

@(
    "# TITAN FAGBI manuscript release payload",
    "",
    "This directory prepares the manuscript-associated release payload for the tagged repository archive $ReleaseTag.",
    "",
    "## Code and data availability statement",
    "",
    $codeAvailability,
    "",
    "## Payload contents",
    "",
    "- reference_outputs/: CSV outputs copied from the verified local benchmark runs for Tables 3-6 and Supplementary Tables S1-S6.",
    "- apbs_context/: APBS input files, APBS stdout/stderr logs, APBS grid-size settings, and the PQR inputs used for Supplementary Figure S1.",
    "- table_figure_traceability.csv: table/figure to source-output/script mapping.",
    "- reference_output_copy_manifest.csv: source-to-archive copy manifest for tabulated reference outputs.",
    "- checksums_sha256.csv: SHA256 checksums for release-facing source, scripts, figures, inputs, and outputs.",
    "",
    "The release manifest at ../release_artifact_manifest.csv is the top-level audit table."
) | Set-Content -LiteralPath $readmePath -Encoding ascii

@(
    "# TITAN manuscript-release build instructions",
    "",
    "These commands document the build and smoke-test environment used by the manuscript release.",
    "",
    "## Windows native build",
    "",
    "Command:",
    "mingw32-make PROFILE=native TARGET=titan.exe all -j2",
    "",
    "The manuscript benchmarks were executed with native double-precision CPU code and AVX2-capable kernels when available.",
    "",
    "## Minimal guarded p4 smoke",
    "",
    "Commands:",
    "powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_y_treecode_default_probe.ps1 -Quick -TitanExe .\titan.exe",
    "powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fagbi_scope_z4xd_p_order_tiled_cache_gate.ps1 -TitanExe .\titan.exe -Quick -AllowFailedDecision",
    "",
    "## Release manifest regeneration",
    "",
    "Command:",
    "powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fagbi_release_readiness_manifest.ps1",
    "",
    "## APBS context",
    "",
    "Supplementary Figure S1 uses APBS input and log files archived under docs/fagbi_guarded_p4_package/release_payload/apbs_context/. APBS itself is external to TITAN; the archived .in files, stdout/stderr logs, and grid-size settings CSV document the exact context rows used in the manuscript supplement.",
    "",
    "## Notes",
    "",
    "- Production PBE default remains -pbe=grid.",
    "- Guarded local p4 TABI/FMM is an explicit advanced path and is not enabled by default.",
    "- Full-protein dense direct TABI is not required for the manuscript tables; sampled dense direct-operator audits and small dense-oracle rows are archived instead."
) | Set-Content -LiteralPath $buildInstructionsPath -Encoding ascii

Add-ManifestEntry -Entries $entries -Path $readmePath -Category "release_document" -ReleaseRole "release payload README and Code/Data statement" -ManuscriptTarget "Code and Data Availability"
Add-ManifestEntry -Entries $entries -Path $buildInstructionsPath -Category "build_instructions" -ReleaseRole "build and smoke-test instructions for the archived release" -ManuscriptTarget "Code and Data Availability"

$entriesSorted = @($entries | Sort-Object category, manuscript_target, relative_path)
$entriesSorted | Export-Csv -LiteralPath $csvPath -NoTypeInformation
$traceRows | Export-Csv -LiteralPath $tracePath -NoTypeInformation
$copiedOutputs | Export-Csv -LiteralPath $copiedPath -NoTypeInformation

$checksums = $entriesSorted |
    Where-Object { $_.status -eq "present" } |
    Select-Object relative_path, bytes, sha256, category, manuscript_target
$checksums | Export-Csv -LiteralPath $checksumPath -NoTypeInformation

$missing = @($entriesSorted | Where-Object { $_.status -ne "present" })
$metadata = [pscustomobject]@{
    release_status = $(if ($missing.Count -eq 0) { "ready_for_tag_and_deposition" } else { "missing_artifacts" })
    generated_utc = (Get-Date).ToUniversalTime().ToString("o")
    git_commit = $gitCommit
    release_tag = $ReleaseTag
    license = "MIT"
    public_repository_url = $PublicRepositoryUrl
    archive_doi = $ZenodoDoi
    entry_count = $entriesSorted.Count
    missing_count = $missing.Count
    table_figure_traceability = Get-RelativePath -Path $tracePath
    checksum_manifest = Get-RelativePath -Path $checksumPath
    entries = $entriesSorted
}
$metadata | ConvertTo-Json -Depth 7 | Set-Content -LiteralPath $jsonPath -Encoding utf8

@(
    "release_status=$(if ($missing.Count -eq 0) { "ready_for_tag_and_deposition" } else { "missing_artifacts" })",
    "license=MIT",
    "release_tag=$ReleaseTag",
    "public_repository_url=$PublicRepositoryUrl",
    "archive_doi=$ZenodoDoi",
    "git_commit=$gitCommit",
    "entry_count=$($entriesSorted.Count)",
    "missing_count=$($missing.Count)",
    "manifest_csv=$(Get-RelativePath -Path $csvPath)",
    "manifest_json=$(Get-RelativePath -Path $jsonPath)",
    "traceability_csv=$(Get-RelativePath -Path $tracePath)",
    "checksums_csv=$(Get-RelativePath -Path $checksumPath)"
) | Set-Content -LiteralPath $summaryPath -Encoding ascii

$reportLines = @()
$reportLines += "TITAN FAGBI release-readiness manifest"
$reportLines += "======================================"
$reportLines += ""
$reportLines += "Status: $(if ($missing.Count -eq 0) { "ready_for_tag_and_deposition" } else { "missing_artifacts" })"
$reportLines += "License: MIT"
$reportLines += "Release tag: $ReleaseTag"
$reportLines += "Public repository URL: $PublicRepositoryUrl"
$reportLines += "Archive DOI: $ZenodoDoi"
$reportLines += "Git commit: $gitCommit"
$reportLines += "Manifest CSV: $(Get-RelativePath -Path $csvPath)"
$reportLines += "Manifest JSON: $(Get-RelativePath -Path $jsonPath)"
$reportLines += "Traceability CSV: $(Get-RelativePath -Path $tracePath)"
$reportLines += "Checksums CSV: $(Get-RelativePath -Path $checksumPath)"
$reportLines += "Entries: $($entriesSorted.Count)"
$reportLines += "Missing required artifacts: $($missing.Count)"
if ($missing.Count -gt 0) {
    $reportLines += ""
    $reportLines += "Missing artifacts:"
    foreach ($m in $missing) {
        $reportLines += "- $($m.relative_path) [$($m.category); $($m.manuscript_target)]"
    }
}
$reportLines += ""
$reportLines += "Prepared Code and Data Availability text:"
$reportLines += $codeAvailability
$reportLines | Set-Content -LiteralPath $rootReport -Encoding ascii

Write-Host "Release readiness manifest: $csvPath"
Write-Host "Release payload: $payloadRoot"
Write-Host "Release readiness report: $rootReport"
if ($missing.Count -gt 0) {
    Write-Host "Missing required artifacts: $($missing.Count)"
    exit 2
}
