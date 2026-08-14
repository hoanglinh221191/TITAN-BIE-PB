# =============================================================================
# pbe_precond_cohort_n5_gate.ps1
#
# Extends the PBE preconditioner-robustness cohort from three test geometries
# (barnase, 1a63, 1sh1) to five by adding 1a7m and 1ubq, then re-runs the ablation
# on all three surface representations.
#
# Fixture preparation reuses the exact z4q protocol so that all five proteins
# share one preparation policy:
#     pdb2pqr 3.7.1 --ff=AMBER --ffout=AMBER --keep-chain
#                   --titration-state-method=propka --with-ph=7.0
#     postprocess-radius-floor = 0.6000 A
#
# The ablation driver caches completed rows only when titan.exe, solver contract,
# PQR input, and numerical protocol hashes all match. A new build therefore
# deliberately invalidates rows produced by another executable.
#
# Usage (from repo root, PowerShell):
#     .\test\pbe_precond_cohort_n5_gate.ps1 -RefreshFixtures
#     .\test\pbe_precond_cohort_n5_gate.ps1                  # fixtures already prepared
#
# Outputs (consumed by pbe_precond_cohort_n5_tables.py):
#     test\tmp\fagbi_graph_preconditioner_ablation\fagbi_graph_preconditioner_matrix.csv
#     test\tmp\fagbi_nodepatch_analytic_ses_preconditioner\analytic_ses_ablation\
#         fagbi_graph_preconditioner_matrix.csv
# =============================================================================
[CmdletBinding()]
param(
    [string]$TitanExe = '.\titan.exe',
    [string]$RepoRoot = '',
    [string]$FixtureRoot = '',
    [string[]]$Fixtures = @('barnase', '1a63', '1sh1', '1a7m', '1ubq'),
    [string[]]$NewFixtures = @('1a7m', '1ubq'),
    [ValidateRange(0.25, 4.0)][double]$SesGridSpacing = 1.25,
    [ValidateRange(1, 86400)][int]$CaseTimeoutSeconds = 3600,
    [ValidateRange(256,8192)][int]$RasMemoryBudgetMb = 1024,
    [string]$PythonExe = 'python',
    [string]$Pdb2PqrTarget = '',
    [switch]$RefreshFixtures,
    [switch]$SkipSpherePatch,
    [switch]$SkipGridSes,
    [switch]$SkipAnalyticSes,
    [switch]$RecoverVerifiedProfiles,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$Inv = [Globalization.CultureInfo]::InvariantCulture

if ([string]::IsNullOrWhiteSpace($RepoRoot)) { $RepoRoot = Split-Path -Parent $PSScriptRoot }
if ([string]::IsNullOrWhiteSpace($FixtureRoot)) { $FixtureRoot = Join-Path $PSScriptRoot 'pbe_inputs_public' }
if ([string]::IsNullOrWhiteSpace($Pdb2PqrTarget)) { $Pdb2PqrTarget = Join-Path $RepoRoot '.codex_tmp\pdb2pqr-3.7.1' }
$RepoRoot = [IO.Path]::GetFullPath($RepoRoot)
$FixtureRoot = [IO.Path]::GetFullPath($FixtureRoot)
if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$pdb2pqrVersion = '3.7.1'
$pdb2pqrArguments = '--ff=AMBER --ffout=AMBER --keep-chain --titration-state-method=propka --with-ph=7.0'
$preparationPolicy = "$pdb2pqrArguments; postprocess-radius-floor=0.6000A"
$radiusFloorAngstrom = 0.6
$manifestCsv = Join-Path $FixtureRoot 'public_fixture_manifest.csv'
$legacyZ4qRoot = Join-Path $PSScriptRoot 'pbe_inputs_public_z4q'

# --- fixture specifications for the two additions ---------------------------
# 1a7m widens the cohort upward in size; 1ubq is a standard mid-size control
# whose AMBER PQR already exists under pbe_inputs_public_z4q and is reused
# verbatim when present.
$specs = @(
    [pscustomobject]@{
        fixture = '1a7m'; pdb_id = '1A7M'; protein_name = 'Human interleukin-6'
        model = 1; chain = 'A'
        source_url = 'https://www.rcsb.org/structure/1A7M'
        download_url = 'https://files.rcsb.org/download/1A7M.pdb'
        raw_name = '1a7m_source.pdb'; pdb_name = '1a7m_model1_chainA.pdb'
        pqr_name = '1a7m_full_amber_ph7.pqr'
        reuse_from = ''
    },
    [pscustomobject]@{
        fixture = '1ubq'; pdb_id = '1UBQ'; protein_name = 'Ubiquitin'
        model = 1; chain = 'A'
        source_url = 'https://www.rcsb.org/structure/1UBQ'
        download_url = 'https://files.rcsb.org/download/1UBQ.pdb'
        raw_name = '1ubq_source.pdb'; pdb_name = '1ubq_model1_chainA.pdb'
        pqr_name = '1ubq_full_amber_ph7.pqr'
        reuse_from = $legacyZ4qRoot
    }
)

# --- helpers (byte-for-byte equivalent to the z4q preparation gate) ----------
function Get-Sha256([string]$Path) {
    return (Get-FileHash -Algorithm SHA256 -LiteralPath $Path).Hash.ToLowerInvariant()
}

function Get-AtomCount([string]$Path) {
    return @(Get-Content -LiteralPath $Path |
        Where-Object { $_.StartsWith('ATOM') -or $_.StartsWith('HETATM') }).Count
}

function Write-ModelChainPdb([string]$InputPath, [string]$OutputPath, [string]$Chain, [int]$Model) {
    $hasModels = @(Select-String -LiteralPath $InputPath -Pattern '^MODEL\s').Count -gt 0
    $activeModel = if ($hasModels) { 0 } else { 1 }
    $out = New-Object System.Collections.Generic.List[string]
    foreach ($line in Get-Content -LiteralPath $InputPath) {
        if ($line.StartsWith('MODEL')) {
            $n = 0
            if ([int]::TryParse($line.Substring(5).Trim(), [ref]$n)) { $activeModel = $n }
            continue
        }
        if ($line.StartsWith('ENDMDL')) {
            if ($activeModel -eq $Model) { break }
            $activeModel = 0
            continue
        }
        if ($activeModel -ne $Model -or -not $line.StartsWith('ATOM')) { continue }
        $lineChain = if ($line.Length -gt 21) { [string]$line[21] } else { '' }
        if ($lineChain -eq $Chain) { $out.Add($line) | Out-Null }
    }
    if ($out.Count -eq 0) { throw "No ATOM rows for model $Model chain $Chain in $InputPath" }
    $out.Add('END') | Out-Null
    Set-Content -LiteralPath $OutputPath -Value $out -Encoding ASCII
}

function Install-Pdb2Pqr {
    if (-not (Test-Path -LiteralPath $PythonExe)) { throw "Missing bundled Python: $PythonExe" }
    New-Item -ItemType Directory -Path $Pdb2PqrTarget -Force | Out-Null
    $marker = Join-Path $Pdb2PqrTarget 'pdb2pqr\__init__.py'
    if (-not (Test-Path -LiteralPath $marker)) {
        & $PythonExe -m pip install --disable-pip-version-check --target $Pdb2PqrTarget ("pdb2pqr=={0}" -f $pdb2pqrVersion)
        if ($LASTEXITCODE -ne 0) { throw 'Unable to install pinned PDB2PQR 3.7.1.' }
    }
}

function Invoke-Pdb2Pqr([string]$InputPath, [string]$OutputPath) {
    $oldPythonPath = $env:PYTHONPATH
    try {
        $env:PYTHONPATH = if ([string]::IsNullOrWhiteSpace($oldPythonPath)) { $Pdb2PqrTarget } else { "$Pdb2PqrTarget;$oldPythonPath" }
        & $PythonExe -m pdb2pqr --ff=AMBER --ffout=AMBER --keep-chain `
            --titration-state-method=propka --with-ph=7.0 $InputPath $OutputPath
        if ($LASTEXITCODE -ne 0) { throw "PDB2PQR failed for $InputPath" }
    } finally {
        $env:PYTHONPATH = $oldPythonPath
    }
}

function Normalize-PqrRadii([string]$Path, [double]$MinimumRadius) {
    $culture = [Globalization.CultureInfo]::InvariantCulture
    $output = New-Object System.Collections.Generic.List[string]
    $clamped = 0
    foreach ($line in Get-Content -LiteralPath $Path) {
        if (-not ($line.StartsWith('ATOM') -or $line.StartsWith('HETATM'))) {
            $output.Add($line) | Out-Null
            continue
        }
        $fields = @($line -split '\s+' | Where-Object { $_ -ne '' })
        if ($fields.Count -lt 10) { throw "Malformed PDB2PQR atom row: $line" }
        $radius = 0.0
        if (-not [double]::TryParse($fields[-1], [Globalization.NumberStyles]::Float, $culture, [ref]$radius)) {
            throw "Malformed PDB2PQR radius: $line"
        }
        if ($radius -le 0.0) {
            $tokenStart = $line.LastIndexOf($fields[-1], [StringComparison]::Ordinal)
            if ($tokenStart -lt 0) { throw "Unable to locate radius token: $line" }
            $line = $line.Substring(0, $tokenStart) + $MinimumRadius.ToString('F4', $culture)
            ++$clamped
        }
        $output.Add($line) | Out-Null
    }
    Set-Content -LiteralPath $Path -Value $output -Encoding ASCII
    return $clamped
}

# --- step 1: prepare the two new fixtures -----------------------------------
if ($RefreshFixtures) {
    Write-Host '=== step 1: preparing new fixtures ==='
    New-Item -ItemType Directory -Path $FixtureRoot -Force | Out-Null
    $needsPdb2Pqr = @($specs | Where-Object {
        if ($_.fixture -notin $NewFixtures) { return $false }
        if ([string]::IsNullOrWhiteSpace($_.reuse_from)) { return $true }
        return -not (Test-Path -LiteralPath (Join-Path $_.reuse_from $_.pqr_name))
    }).Count -gt 0
    if ($needsPdb2Pqr) { Install-Pdb2Pqr }

    $newRows = New-Object System.Collections.Generic.List[object]
    foreach ($spec in $specs) {
        if ($spec.fixture -notin $NewFixtures) { continue }
        $rawPath = Join-Path $FixtureRoot $spec.raw_name
        $pdbPath = Join-Path $FixtureRoot $spec.pdb_name
        $pqrPath = Join-Path $FixtureRoot $spec.pqr_name

        $reusePqr = if ([string]::IsNullOrWhiteSpace($spec.reuse_from)) { '' } else { Join-Path $spec.reuse_from $spec.pqr_name }
        if ($reusePqr -and (Test-Path -LiteralPath $reusePqr)) {
            # Identical preparation policy already on disk: copy, do not re-run.
            Write-Host ("  {0}: reusing prepared AMBER PQR from {1}" -f $spec.fixture, $spec.reuse_from)
            Copy-Item -LiteralPath $reusePqr -Destination $pqrPath -Force
            foreach ($n in @($spec.raw_name, $spec.pdb_name)) {
                $src = Join-Path $spec.reuse_from $n
                if (Test-Path -LiteralPath $src) { Copy-Item -LiteralPath $src -Destination (Join-Path $FixtureRoot $n) -Force }
            }
            $reuseManifest = Join-Path $spec.reuse_from 'public_fixture_manifest.csv'
            $reuseRow = @()
            if (Test-Path -LiteralPath $reuseManifest) {
                $reuseRow = @(Import-Csv -LiteralPath $reuseManifest |
                    Where-Object fixture -eq $spec.fixture)
            }
            $radiusClampedAtoms = if ($reuseRow.Count -eq 1) {
                [int]$reuseRow[0].radius_clamped_atoms
            } else { 0 }
        }
        else {
            Write-Host ("  {0}: downloading {1}" -f $spec.fixture, $spec.download_url)
            Invoke-WebRequest -UseBasicParsing -Uri $spec.download_url -OutFile $rawPath
            Write-ModelChainPdb -InputPath $rawPath -OutputPath $pdbPath -Chain $spec.chain -Model $spec.model
            Write-Host ("  {0}: running PDB2PQR {1}" -f $spec.fixture, $pdb2pqrVersion)
            Invoke-Pdb2Pqr -InputPath $pdbPath -OutputPath $pqrPath
            $radiusClampedAtoms = Normalize-PqrRadii -Path $pqrPath -MinimumRadius $radiusFloorAngstrom
        }

        $atomCount = Get-AtomCount $pqrPath
        if ($atomCount -lt 20) { throw "Unexpectedly small PQR fixture: $pqrPath" }
        Write-Host ("  {0}: {1} atoms, radius-clamped {2}" -f $spec.fixture, $atomCount, $radiusClampedAtoms)

        $newRows.Add([pscustomobject]@{
            fixture               = $spec.fixture
            pdb_id                = $spec.pdb_id
            protein_name          = $spec.protein_name
            model                 = $spec.model
            chain                 = $spec.chain
            source_url            = $spec.source_url
            raw_pdb               = if (Test-Path -LiteralPath $rawPath) { $spec.raw_name } else { '' }
            prepared_pdb          = if (Test-Path -LiteralPath $pdbPath) { $spec.pdb_name } else { '' }
            pqr                   = $spec.pqr_name
            raw_pdb_sha256        = if (Test-Path -LiteralPath $rawPath) { Get-Sha256 $rawPath } else { '' }
            prepared_pdb_sha256   = if (Test-Path -LiteralPath $pdbPath) { Get-Sha256 $pdbPath } else { '' }
            pdb_sha256            = if (Test-Path -LiteralPath $pdbPath) { Get-Sha256 $pdbPath } else { '' }
            pqr_sha256            = Get-Sha256 $pqrPath
            atom_count            = $atomCount
            pdb2pqr_version       = $pdb2pqrVersion
            preparation           = $preparationPolicy
            radius_floor_angstrom = $radiusFloorAngstrom
            radius_clamped_atoms  = $radiusClampedAtoms
        }) | Out-Null
    }

    # merge into the cohort manifest, replacing any pre-existing row of the same name
    if (-not (Test-Path -LiteralPath $manifestCsv)) { throw "Missing cohort manifest: $manifestCsv" }
    Copy-Item -LiteralPath $manifestCsv -Destination "$manifestCsv.bak" -Force
    $columns = @((Import-Csv -LiteralPath $manifestCsv)[0].PSObject.Properties.Name)
    $merged = New-Object System.Collections.Generic.List[object]
    foreach ($row in @(Import-Csv -LiteralPath $manifestCsv)) {
        if ([string]$row.fixture -notin $NewFixtures) { $merged.Add($row) | Out-Null }
    }
    foreach ($row in $newRows) {
        $ordered = [ordered]@{}
        foreach ($c in $columns) {
            $prop = $row.PSObject.Properties[$c]
            $ordered[$c] = if ($null -eq $prop) { '' } else { $prop.Value }
        }
        $merged.Add([pscustomobject]$ordered) | Out-Null
    }
    $merged.ToArray() | Export-Csv -LiteralPath $manifestCsv -NoTypeInformation
    Write-Host ("  manifest updated: {0} rows (backup at {1}.bak)" -f $merged.Count, $manifestCsv)
}

# --- step 2: verify every requested fixture resolves ------------------------
Write-Host '=== step 2: verifying cohort manifest ==='
$manifest = @(Import-Csv -LiteralPath $manifestCsv)
foreach ($name in $Fixtures) {
    $rows = @($manifest | Where-Object fixture -eq $name)
    if ($rows.Count -ne 1) { throw "Expected exactly one manifest row for '$name' (found $($rows.Count)). Re-run with -RefreshFixtures." }
    $pqr = Join-Path $FixtureRoot $rows[0].pqr
    if (-not (Test-Path -LiteralPath $pqr)) { throw "Missing fixture PQR: $pqr" }
    Write-Host ("  ok  {0,-8} {1,6} atoms  {2}" -f $name, $rows[0].atom_count, $rows[0].pqr)
}

# --- step 3: run the ablation on all three surface representations ----------
$ablationScript = Join-Path $PSScriptRoot 'pbe_fagbi_graph_preconditioner_ablation.ps1'
$analyticGate = Join-Path $PSScriptRoot 'pbe_fagbi_nodepatch_analytic_ses_preconditioner_gate.ps1'
foreach ($p in @($ablationScript, $analyticGate)) {
    if (-not (Test-Path -LiteralPath $p)) { throw "Missing driver script: $p" }
}

$common = @{
    TitanExe            = $TitanExe
    Fixtures            = $Fixtures
    SesGridSpacing      = $SesGridSpacing
    CaseTimeoutSeconds  = $CaseTimeoutSeconds
    RasMemoryBudgetMb   = $RasMemoryBudgetMb
    AllowFailedDecision = $true
}
if ($Force) { $common.Force = $true }
if ($RecoverVerifiedProfiles) { $common.RecoverVerifiedProfiles = $true }

if (-not $SkipSpherePatch) {
    Write-Host '=== step 3a: sphere-patch stress representation ==='
    & $ablationScript @common -PrimarySurface 'sphere-patch'
}
if (-not $SkipGridSes) {
    Write-Host '=== step 3b: watertight Grid-SES ==='
    & $ablationScript @common -PrimarySurface 'grid-ses'
}
if (-not $SkipAnalyticSes) {
    Write-Host '=== step 3c: analytic patch-ray SES (node-patch) ==='
    $analyticArgs = @{
        TitanExe            = $TitanExe
        Fixtures            = $Fixtures
        SesGridSpacing      = $SesGridSpacing
        CaseTimeoutSeconds  = $CaseTimeoutSeconds
        RasMemoryBudgetMb   = $RasMemoryBudgetMb
        AllowFailedDecision = $true
    }
    if ($Force) { $analyticArgs.Force = $true }
    if ($RecoverVerifiedProfiles) { $analyticArgs.RecoverVerifiedProfiles = $true }
    & $analyticGate @analyticArgs
}

# Fail closed on incomplete, mixed-provenance, fallback, convergence, or
# physical-surface operator-audit gaps before declaring the cohort complete.
$tableScript = Join-Path $PSScriptRoot 'pbe_precond_cohort_n5_tables.py'
if (-not (Test-Path -LiteralPath $tableScript)) { throw "Missing table validator: $tableScript" }
Write-Host '=== step 4: validating complete n=5 PBE cohort ==='
& $PythonExe $tableScript --validate-only
if ($LASTEXITCODE -ne 0) { throw "n=5 cohort validation failed (exit $LASTEXITCODE)." }

Write-Host ''
Write-Host '=== done ==='
Write-Host 'Matrix CSVs:'
Write-Host ("  {0}" -f (Join-Path $RepoRoot 'test\tmp\fagbi_graph_preconditioner_ablation\fagbi_graph_preconditioner_matrix.csv'))
Write-Host ("  {0}" -f (Join-Path $RepoRoot 'test\tmp\fagbi_nodepatch_analytic_ses_preconditioner\analytic_ses_ablation\fagbi_graph_preconditioner_matrix.csv'))
Write-Host ''
Write-Host 'Next: python .\test\pbe_precond_cohort_n5_tables.py --out tables_n5.md'
