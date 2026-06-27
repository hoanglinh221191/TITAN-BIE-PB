param(
    [string]$OutputRoot = "docs\fagbi_guarded_p4_package"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Get-RelativePath {
    param([string]$Path)
    $full = (Resolve-Path -LiteralPath $Path).Path
    $root = (Resolve-Path -LiteralPath ".").Path
    if ($full.StartsWith($root, [System.StringComparison]::OrdinalIgnoreCase)) {
        return $full.Substring($root.Length).TrimStart("\", "/") -replace "\\", "/"
    }
    return $full
}

function Add-FileEntry {
    param(
        [System.Collections.Generic.List[object]]$Entries,
        [string]$Path,
        [string]$Category,
        [string]$Role
    )
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        $Entries.Add([pscustomobject]@{
            category = $Category
            release_role = $Role
            relative_path = ($Path -replace "\\", "/")
            status = "missing"
            bytes = ""
            sha256 = ""
            last_write_utc = ""
        })
        return
    }
    $item = Get-Item -LiteralPath $Path
    $hash = Get-FileHash -LiteralPath $Path -Algorithm SHA256
    $Entries.Add([pscustomobject]@{
        category = $Category
        release_role = $Role
        relative_path = Get-RelativePath -Path $Path
        status = "present"
        bytes = [int64]$item.Length
        sha256 = $hash.Hash.ToLowerInvariant()
        last_write_utc = $item.LastWriteTimeUtc.ToString("o")
    })
}

function Add-GlobEntries {
    param(
        [System.Collections.Generic.List[object]]$Entries,
        [string]$BasePath,
        [string[]]$Include,
        [string]$Category,
        [string]$Role
    )
    if (-not (Test-Path -LiteralPath $BasePath -PathType Container)) {
        return
    }
    Get-ChildItem -LiteralPath $BasePath -Recurse -File |
        Where-Object {
            $name = $_.Name
            foreach ($pattern in $Include) {
                if ($name -like $pattern) { return $true }
            }
            return $false
        } |
        Sort-Object FullName |
        ForEach-Object {
            Add-FileEntry -Entries $Entries -Path $_.FullName -Category $Category -Role $Role
        }
}

$outRoot = New-Item -ItemType Directory -Force -Path $OutputRoot
$entries = [System.Collections.Generic.List[object]]::new()

$requiredFiles = @(
    @{ path = "LICENSE"; category = "license"; role = "MIT license for the public release" },
    @{ path = "CITATION.cff"; category = "citation"; role = "citation metadata template for the archived release" },
    @{ path = ".zenodo.json"; category = "archive_metadata"; role = "Zenodo metadata template for DOI minting" },
    @{ path = "README.md"; category = "documentation"; role = "repository overview" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\titan_fagbi_reader_friendly_p2_cap_revised.docx"; category = "manuscript"; role = "latest user-edited main manuscript" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\titan_fagbi_reader_friendly_supplement_revised.docx"; category = "supplement"; role = "latest user-edited supplement" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\fagbi_scope_z4t_manuscript_source.md"; category = "manuscript_source"; role = "markdown source for manuscript text" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\fagbi_scope_z4t_references.bib"; category = "references"; role = "BibTeX references" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\fagbi_scope_z4t_references.csv"; category = "references"; role = "reference table" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\fagbi_scope_z4t_performance_ladder_table.csv"; category = "table_input"; role = "algorithmic ladder table input" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\fagbi_scope_z4t_1a63_runtime_composition.csv"; category = "table_input"; role = "1A63 runtime composition input" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\fagbi_scope_z4t_preconditioner_residual_trace.csv"; category = "figure_input"; role = "GMRES residual profile trace" },
    @{ path = "docs\fagbi_guarded_p4_package\submission_preview_synced\fagbi_scope_z4x_global_p_order_coverage.csv"; category = "table_input"; role = "p-order coverage table input" },
    @{ path = "test\pbe_inputs_public\public_fixture_manifest.csv"; category = "fixture_manifest"; role = "primary public-protein fixture provenance" },
    @{ path = "test\pbe_inputs_public_z4q\public_fixture_expansion_manifest.csv"; category = "fixture_manifest"; role = "supplementary public-protein fixture provenance" }
)

foreach ($file in $requiredFiles) {
    Add-FileEntry -Entries $entries -Path $file.path -Category $file.category -Role $file.role
}

Add-GlobEntries -Entries $entries -BasePath "docs\fagbi_guarded_p4_package\submission_preview_synced\figures" -Include @("*.svg", "*.png") -Category "figure" -Role "manuscript and supplement figures"
Add-GlobEntries -Entries $entries -BasePath "src" -Include @("*.cpp", "*.h", "*.hpp", "*.cxx") -Category "source_code" -Role "TITAN source code"

$scriptList = @(
    "test\pbe_fagbi_scope_z4v_performance_ladder_matrix.ps1",
    "test\pbe_fagbi_scope_z4w_ablation_sensitivity_matrix.ps1",
    "test\pbe_fagbi_scope_z4x_dense_small_oracle.ps1",
    "test\pbe_fagbi_scope_z4x_p_order_protocol_sweep.ps1",
    "test\pbe_fagbi_scope_z4x_p5_feasibility_gate.ps1",
    "test\pbe_fagbi_scope_z4xd_p_order_tiled_cache_gate.ps1",
    "test\pbe_fagbi_scope_z4x_global_p_order_figure.ps1",
    "test\pbe_fagbi_scope_z4xe_supplementary_sensitivity_package.ps1",
    "test\pbe_fagbi_scope_z4xf_ilu_neighbor_cap_sweep.ps1",
    "test\pbe_fagbi_release_readiness_manifest.ps1"
)
foreach ($script in $scriptList) {
    Add-FileEntry -Entries $entries -Path $script -Category "benchmark_script" -Role "reported-table or release manifest reproduction script"
}

$reportList = @(
    "TITAN_FMM_M2L_SCOPE_Z4T_SUBMISSION_POLISH_REPORT.txt",
    "TITAN_FMM_M2L_SCOPE_Z4W_ABLATION_SENSITIVITY_REPORT.txt",
    "TITAN_FMM_M2L_SCOPE_Z4X_DENSE_SMALL_ORACLE_REPORT.txt",
    "TITAN_FMM_M2L_SCOPE_Z4XB_P_ORDER_SWEEP_REPORT.txt",
    "TITAN_FMM_M2L_SCOPE_Z4XD_P5_TILED_CACHE_REPORT.txt",
    "TITAN_FMM_M2L_SCOPE_Z4XE_SUPPLEMENTARY_SENSITIVITY_REPORT.txt",
    "TITAN_FMM_M2L_SCOPE_Z4XF_ILU_NEIGHBOR_CAP_SWEEP_REPORT.txt"
)
foreach ($report in $reportList) {
    Add-FileEntry -Entries $entries -Path $report -Category "report" -Role "benchmark gate report"
}

$csvPath = Join-Path $outRoot.FullName "release_artifact_manifest.csv"
$jsonPath = Join-Path $outRoot.FullName "release_artifact_manifest.json"
$summaryPath = Join-Path $outRoot.FullName "release_readiness_summary.dat"
$rootReport = "TITAN_FAGBI_RELEASE_READINESS_REPORT.txt"

$entries | Sort-Object category, relative_path | Export-Csv -LiteralPath $csvPath -NoTypeInformation

$gitCommit = ""
try {
    $gitCommit = (& git rev-parse HEAD).Trim()
} catch {
    $gitCommit = "unavailable"
}

$missing = @($entries | Where-Object { $_.status -ne "present" })
$metadata = [pscustomobject]@{
    release_status = "pre_submission_doi_pending"
    generated_utc = (Get-Date).ToUniversalTime().ToString("o")
    git_commit = $gitCommit
    license = "MIT"
    public_repository_url = "pending"
    archive_doi = "pending"
    entry_count = $entries.Count
    missing_count = $missing.Count
    entries = @($entries | Sort-Object category, relative_path)
}
$metadata | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $jsonPath -Encoding utf8

@(
    "release_status=pre_submission_doi_pending",
    "license=MIT",
    "public_repository_url=pending",
    "archive_doi=pending",
    "git_commit=$gitCommit",
    "entry_count=$($entries.Count)",
    "missing_count=$($missing.Count)",
    "manifest_csv=$(Get-RelativePath -Path $csvPath)",
    "manifest_json=$(Get-RelativePath -Path $jsonPath)"
) | Set-Content -LiteralPath $summaryPath -Encoding ascii

$reportLines = @()
$reportLines += "TITAN FAGBI release-readiness manifest"
$reportLines += "======================================"
$reportLines += ""
$reportLines += "Status: pre-submission; public repository URL and DOI pending."
$reportLines += "License: MIT."
$reportLines += "Git commit: $gitCommit"
$reportLines += "Manifest CSV: $(Get-RelativePath -Path $csvPath)"
$reportLines += "Manifest JSON: $(Get-RelativePath -Path $jsonPath)"
$reportLines += "Entries: $($entries.Count)"
$reportLines += "Missing required artifacts: $($missing.Count)"
if ($missing.Count -gt 0) {
    $reportLines += ""
    $reportLines += "Missing artifacts:"
    foreach ($m in $missing) {
        $reportLines += "- $($m.relative_path) [$($m.category)]"
    }
}
$reportLines += ""
$reportLines += "Manuscript statement prepared:"
$reportLines += "The TITAN source code, benchmark manifests, scripts, and reference outputs corresponding to the reported tables will be released under the MIT License at a versioned public repository and archived with a persistent DOI before submission."
$reportLines | Set-Content -LiteralPath $rootReport -Encoding ascii

Write-Host "Release readiness manifest: $csvPath"
Write-Host "Release readiness report: $rootReport"
if ($missing.Count -gt 0) {
    Write-Host "Missing required artifacts: $($missing.Count)"
    exit 2
}
