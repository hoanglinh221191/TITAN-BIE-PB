param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$RepoRoot = (Split-Path -Parent $PSScriptRoot),
    [string]$FixtureRoot = (Join-Path $PSScriptRoot 'pbe_inputs_public_z4q'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4q'),
    [string]$PythonExe = 'C:\Users\ADMIN\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [string]$Pdb2PqrTarget = '',
    [int]$CaseTimeoutSeconds = 1800,
    [switch]$Quick,
    [switch]$RefreshFixtures,
    [switch]$AllowFailedDecision
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

if ([string]::IsNullOrWhiteSpace($Pdb2PqrTarget)) {
    $Pdb2PqrTarget = Join-Path $RepoRoot '.codex_tmp\pdb2pqr-3.7.1'
}
$pdb2pqrVersion = '3.7.1'
$pdb2pqrArguments = '--ff=AMBER --ffout=AMBER --keep-chain --titration-state-method=propka --with-ph=7.0'
$preparationPolicy = "$pdb2pqrArguments; postprocess-radius-floor=0.6000A"
$radiusFloorAngstrom = 0.6
$manifestCsv = Join-Path $FixtureRoot 'public_fixture_expansion_manifest.csv'
$matrixCsv = Join-Path $WorkRoot 'fagbi_scope_z4q_matrix.csv'
$repeatCsv = Join-Path $WorkRoot 'fagbi_scope_z4q_repeatability.csv'
$decisionCsv = Join-Path $WorkRoot 'fagbi_scope_z4q_decision.csv'
$summaryPath = Join-Path $WorkRoot 'fagbi_scope_z4q_summary.dat'
$plotPath = Join-Path $WorkRoot 'fagbi_scope_z4q_repeatability_plot.svg'
$rootReport = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z4Q_PUBLIC_FIXTURE_EXPANSION_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$specs = @(
    [pscustomobject]@{ fixture='1ubq'; pdb_id='1UBQ'; protein_name='Ubiquitin'; model=1; chain='A'; source_url='https://www.rcsb.org/structure/1UBQ'; download_url='https://files.rcsb.org/download/1UBQ.pdb'; raw_name='1ubq_source.pdb'; pdb_name='1ubq_model1_chainA.pdb'; pqr_name='1ubq_full_amber_ph7.pqr' },
    [pscustomobject]@{ fixture='2ci2'; pdb_id='2CI2'; protein_name='Chymotrypsin inhibitor 2'; model=1; chain='I'; source_url='https://www.rcsb.org/structure/2CI2'; download_url='https://files.rcsb.org/download/2CI2.pdb'; raw_name='2ci2_source.pdb'; pdb_name='2ci2_model1_chainI.pdb'; pqr_name='2ci2_full_amber_ph7.pqr' },
    [pscustomobject]@{ fixture='1vii'; pdb_id='1VII'; protein_name='Villin headpiece subdomain'; model=1; chain='A'; source_url='https://www.rcsb.org/structure/1VII'; download_url='https://files.rcsb.org/download/1VII.pdb'; raw_name='1vii_source.pdb'; pdb_name='1vii_model1_chainA.pdb'; pqr_name='1vii_full_amber_ph7.pqr' }
)

function Get-Z4qSha256([string]$Path) {
    return (Get-FileHash -Algorithm SHA256 -LiteralPath $Path).Hash.ToLowerInvariant()
}

function Get-Z4qAtomCount([string]$Path) {
    return @(Get-Content -LiteralPath $Path | Where-Object { $_.StartsWith('ATOM') -or $_.StartsWith('HETATM') }).Count
}

function Write-Z4qModelChainPdb([string]$InputPath,[string]$OutputPath,[string]$Chain,[int]$Model) {
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

function Install-Z4qPdb2Pqr {
    if (-not (Test-Path -LiteralPath $PythonExe)) { throw "Missing bundled Python: $PythonExe" }
    New-Item -ItemType Directory -Path $Pdb2PqrTarget -Force | Out-Null
    $marker = Join-Path $Pdb2PqrTarget 'pdb2pqr\__init__.py'
    if (-not (Test-Path -LiteralPath $marker)) {
        & $PythonExe -m pip install --disable-pip-version-check --target $Pdb2PqrTarget ("pdb2pqr=={0}" -f $pdb2pqrVersion)
        if ($LASTEXITCODE -ne 0) { throw 'Unable to install pinned PDB2PQR for Z4Q.' }
    }
}

function Invoke-Z4qPdb2Pqr([string]$InputPath,[string]$OutputPath) {
    $oldPythonPath = $env:PYTHONPATH
    try {
        $env:PYTHONPATH = if ([string]::IsNullOrWhiteSpace($oldPythonPath)) { $Pdb2PqrTarget } else { "$Pdb2PqrTarget;$oldPythonPath" }
        & $PythonExe -m pdb2pqr --ff=AMBER --ffout=AMBER --keep-chain --titration-state-method=propka --with-ph=7.0 $InputPath $OutputPath
        if ($LASTEXITCODE -ne 0) { throw "PDB2PQR failed for $InputPath" }
    } finally {
        $env:PYTHONPATH = $oldPythonPath
    }
}

function Normalize-Z4qPqrRadii([string]$Path,[double]$MinimumRadius) {
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

if ($RefreshFixtures) {
    New-Item -ItemType Directory -Path $FixtureRoot -Force | Out-Null
    Install-Z4qPdb2Pqr
    $rows = New-Object System.Collections.Generic.List[object]
    foreach ($spec in $specs) {
        $rawPath = Join-Path $FixtureRoot $spec.raw_name
        $pdbPath = Join-Path $FixtureRoot $spec.pdb_name
        $pqrPath = Join-Path $FixtureRoot $spec.pqr_name
        Invoke-WebRequest -UseBasicParsing -Uri $spec.download_url -OutFile $rawPath
        Write-Z4qModelChainPdb -InputPath $rawPath -OutputPath $pdbPath -Chain $spec.chain -Model $spec.model
        Invoke-Z4qPdb2Pqr -InputPath $pdbPath -OutputPath $pqrPath
        $radiusClampedAtoms = Normalize-Z4qPqrRadii -Path $pqrPath -MinimumRadius $radiusFloorAngstrom
        $atomCount = Get-Z4qAtomCount $pqrPath
        if ($atomCount -lt 20) { throw "Unexpectedly small PQR fixture: $pqrPath" }
        $rows.Add([pscustomobject]@{
            fixture=$spec.fixture; pdb_id=$spec.pdb_id; protein_name=$spec.protein_name;
            model=$spec.model; chain=$spec.chain; source_url=$spec.source_url;
            raw_pdb=$spec.raw_name; prepared_pdb=$spec.pdb_name; pqr=$spec.pqr_name;
            raw_pdb_sha256=(Get-Z4qSha256 $rawPath); prepared_pdb_sha256=(Get-Z4qSha256 $pdbPath);
            pqr_sha256=(Get-Z4qSha256 $pqrPath); atom_count=$atomCount;
            pdb2pqr_version=$pdb2pqrVersion; preparation=$preparationPolicy;
            radius_floor_angstrom=$radiusFloorAngstrom; radius_clamped_atoms=$radiusClampedAtoms
        }) | Out-Null
    }
    $rows.ToArray() | Export-Csv -LiteralPath $manifestCsv -NoTypeInformation
}

if (-not (Test-Path -LiteralPath $manifestCsv)) {
    $msg = "z4q_fixture_preparation_gap: missing $manifestCsv. Run this script with -RefreshFixtures after confirming network/PDB2PQR availability."
    New-Item -ItemType Directory -Path (Split-Path -Parent $rootReport) -Force | Out-Null
    Set-Content -LiteralPath $rootReport -Value @('# TITAN FMM/M2L Scope Z4Q Public Fixture Expansion Report','', "Decision: z4q_fixture_preparation_gap", '', $msg) -Encoding ASCII
    throw $msg
}

$manifest = @(Import-Csv -LiteralPath $manifestCsv)
$wanted = @(if ($Quick) { '1ubq' } else { '1ubq'; '2ci2'; '1vii' })
foreach ($name in $wanted) {
    $row = @($manifest | Where-Object fixture -eq $name) | Select-Object -First 1
    if ($null -eq $row) { throw "z4q_fixture_preparation_gap: missing manifest row $name" }
    foreach ($file in @($row.raw_pdb,$row.prepared_pdb,$row.pqr)) {
        if (-not (Test-Path -LiteralPath (Join-Path $FixtureRoot $file))) { throw "z4q_fixture_preparation_gap: missing cached file $file" }
    }
    if ((Get-Z4qSha256 (Join-Path $FixtureRoot $row.pqr)) -ne $row.pqr_sha256) { throw "z4q_fixture_preparation_gap: checksum drift for $name" }
    if ((Get-Z4qAtomCount (Join-Path $FixtureRoot $row.pqr)) -ne [int]$row.atom_count) { throw "z4q_fixture_preparation_gap: atom-count drift for $name" }
}

New-FmmCliCleanDir -Path $WorkRoot -WorkRoot (Split-Path -Parent $WorkRoot)
$runtimeRepoRoot = $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $runtimeRepoRoot 'lib_ff'))) {
    $commonGitDir = (& git -C $RepoRoot rev-parse --git-common-dir).Trim()
    if (-not [IO.Path]::IsPathRooted($commonGitDir)) { $commonGitDir = [IO.Path]::GetFullPath((Join-Path $RepoRoot $commonGitDir)) }
    $runtimeRepoRoot = Split-Path -Parent $commonGitDir
}

function Get-Z4qProfileValue([object]$Profile,[string]$Name,[object]$Default='') {
    if ($null -eq $Profile) { return $Default }
    $prop = $Profile.PSObject.Properties[$Name]
    if ($null -eq $prop) { return $Default }
    return $prop.Value
}

function Invoke-Z4qCase([pscustomobject]$Fixture,[string]$Variant,[int]$Repeat) {
    $isP4 = $Variant -eq 'p4_guarded'
    $runDir = Join-Path $WorkRoot "$($Fixture.fixture)_${Variant}_r$Repeat"
    New-Item -ItemType Directory -Path $runDir -Force | Out-Null
    $pqr = Join-Path $FixtureRoot $Fixture.pqr
    $argPath = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $runtimeRepoRoot -PqrPath $pqr `
        -Backend tabi -Mode treecode -FmmMode single-matvec-analytic -SurfacePoints 8 `
        -MaxIter 720 -Restart 128 -Tolerance 1e-6 -Threads 4 -SaltM 0.15 `
        -InternalDielectric 1 -ExternalDielectric 80 -TemperatureK 298.15 `
        -OpeningTheta 0.6 -AuditSampleCount 256 -AuditMode sampled-detail `
        -Preconditioner fmm-nearfield-ilu0 -PreconditionerIluNeighborCap 64 `
        -PreconditionerIluDropTol 0.0 -NearIntegration analytic_singular -TaylorOrder 4 `
        -FmmLocal $(if($isP4){'gmres-p4'}else{'off'}) `
        -FmmLocalRefine $(if($isP4){'operator-tail-static'}else{'off'}) `
        -FmmLocalRefineTolerance 0.002 -FmmLocalCache $(if($isP4){'green-jet-compressed'}else{'off'}) `
        -FmmLocalCacheBudgetMb 256 -FmmLocalFallback $(if($isP4){'packed-avx2'}else{'scalar'}) `
        -FmmLocalPlanBuilder $(if($isP4){'metric-cache'}else{'legacy-rescan'}) `
        -FmmLocalGeometryGuard $(if($isP4){'strict-geometry-mac'}else{'off'}) `
        -FmmLocalGeometryTheta $(if($isP4){0.5}else{0.0})
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $base = [ordered]@{
        fixture=$Fixture.fixture; protein_name=$Fixture.protein_name; pdb_id=$Fixture.pdb_id;
        variant=$Variant; repeat=$Repeat; status='failed'; atom_count=$Fixture.atom_count; panels=0;
        iterations=0; internal_true_residual=[double]::NaN; final_residual=[double]::NaN;
        energy_kj_mol=[double]::NaN; audit_rms=[double]::NaN; sampled_direct_true_residual=[double]::NaN;
        matvec_seconds_avg=[double]::NaN; wall_seconds=$run.wall_seconds; peak_rss_mb=$run.peak_working_set_mb;
        fmm_local_gmres_active='false'; plan_build_count=0; cache_build_count=0; gmres_matvec_calls=0;
        operator_finite='false'; p4_vs_degree1_matvec_ratio=[double]::NaN; failure_log=''; run_dir=$runDir
    }
    if ($run.timed_out) { $base.status='timeout'; $base.failure_log='timeout_not_verified'; return [pscustomobject]$base }
    if ($run.exit_code -ne 0) { $base.failure_log="exit_code_$($run.exit_code)"; return [pscustomobject]$base }
    try {
        $result = Read-FmmCliCaseResult -RunDir $runDir -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
        $profile = @(Import-Csv -LiteralPath (Join-Path $runDir 'pbe_tabi_profile.csv'))[0]
        $base.status=$result.status; $base.panels=$result.panels; $base.iterations=$result.iterations
        $base.internal_true_residual=$result.true_residual; $base.final_residual=$result.final_residual
        $base.energy_kj_mol=$result.energy_kj_mol; $base.audit_rms=$result.tree_vs_direct_matvec_rel
        $base.sampled_direct_true_residual=Convert-FmmCliDouble (Get-Z4qProfileValue $profile 'fmm_local_sampled_direct_true_residual' '0')
        $base.matvec_seconds_avg=$result.matvec_seconds_avg
        $base.fmm_local_gmres_active=Get-Z4qProfileValue $profile 'fmm_local_gmres_active' 'false'
        $base.plan_build_count=Convert-FmmCliInt64 (Get-Z4qProfileValue $profile 'fmm_local_plan_build_count' '0')
        $base.cache_build_count=Convert-FmmCliInt64 (Get-Z4qProfileValue $profile 'fmm_local_cache_build_count' '0')
        $base.gmres_matvec_calls=Convert-FmmCliInt64 (Get-Z4qProfileValue $profile 'fmm_local_gmres_matvec_calls' '0')
        $base.operator_finite=Get-Z4qProfileValue $profile 'fmm_local_operator_finite' 'true'
        $base.failure_log=if($result.status -eq 'ok'){'pass'}else{'not_converged'}
    } catch {
        $base.failure_log='profile_parse_error'
        Set-Content -LiteralPath (Join-Path $runDir 'parse_error.txt') -Value $_.Exception.Message -Encoding UTF8
    }
    return [pscustomobject]$base
}

$rows = New-Object System.Collections.Generic.List[object]
foreach ($name in $wanted) {
    $fixture = @($manifest | Where-Object fixture -eq $name)[0]
    $baseline = Invoke-Z4qCase -Fixture $fixture -Variant 'degree1_baseline' -Repeat 0
    $rows.Add($baseline) | Out-Null
    $repeatCount = if ($Quick) { 1 } else { 2 }
    for ($repeat = 1; $repeat -le $repeatCount; ++$repeat) {
        $row = Invoke-Z4qCase -Fixture $fixture -Variant 'p4_guarded' -Repeat $repeat
        if ((Test-FmmCliFinite ([double]$baseline.matvec_seconds_avg)) -and [double]$baseline.matvec_seconds_avg -gt 0 -and
            (Test-FmmCliFinite ([double]$row.matvec_seconds_avg))) {
            $row.p4_vs_degree1_matvec_ratio = [double]$row.matvec_seconds_avg / [double]$baseline.matvec_seconds_avg
        }
        $rows.Add($row) | Out-Null
    }
}

$rowArray = $rows.ToArray()
$rowArray | Export-Csv -LiteralPath $matrixCsv -NoTypeInformation
$repeatRows = New-Object System.Collections.Generic.List[object]
foreach ($name in $wanted) {
    $set = @($rowArray | Where-Object { $_.fixture -eq $name -and $_.variant -eq 'p4_guarded' })
    $energies = @($set | ForEach-Object { [double]$_.energy_kj_mol } | Where-Object { Test-FmmCliFinite $_ })
    $audits = @($set | ForEach-Object { [double]$_.audit_rms } | Where-Object { Test-FmmCliFinite $_ })
    $directs = @($set | ForEach-Object { [double]$_.sampled_direct_true_residual } | Where-Object { Test-FmmCliFinite $_ })
    $mean = if ($energies.Count) { ($energies | Measure-Object -Average).Average } else { [double]::NaN }
    $spread = if ($energies.Count) { ($energies | Measure-Object -Maximum).Maximum - ($energies | Measure-Object -Minimum).Minimum } else { [double]::NaN }
    $rowPasses = @($set | Where-Object {
        $_.status -eq 'ok' -and [double]$_.internal_true_residual -le 1e-6 -and [double]$_.audit_rms -le 8e-3 -and
        [double]$_.sampled_direct_true_residual -le 8e-3 -and [double]$_.peak_rss_mb -lt 512 -and
        [int64]$_.plan_build_count -eq 1 -and [int64]$_.cache_build_count -eq 1
    }).Count
    $repeatRows.Add([pscustomobject]@{
        fixture=$name; repeats=$set.Count; row_passes=$rowPasses; mean_energy=$mean; energy_spread=$spread;
        audit_max=if($audits.Count){($audits | Measure-Object -Maximum).Maximum}else{[double]::NaN};
        sampled_direct_true_residual_max=if($directs.Count){($directs | Measure-Object -Maximum).Maximum}else{[double]::NaN};
        repeatability_pass=($rowPasses -eq $set.Count -and $set.Count -gt 0 -and (Test-FmmCliFinite $spread) -and $spread -le [Math]::Max(0.5,1e-4*[Math]::Abs($mean)))
    }) | Out-Null
}
$repeatArray = $repeatRows.ToArray()
$repeatArray | Export-Csv -LiteralPath $repeatCsv -NoTypeInformation

$passCount = @($repeatArray | Where-Object { $_.repeatability_pass }).Count
$decision = if ($passCount -eq $wanted.Count) {
    'z4q_expansion_pass'
} elseif ($passCount -gt 0) {
    'z4q_partial_expansion'
} elseif (@($rowArray | Where-Object { $_.status -eq 'timeout' }).Count -gt 0) {
    'z4q_runtime_limited'
} else {
    'z4q_operator_accuracy_gap'
}

$svg = @(
    '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="420" viewBox="0 0 900 420">',
    '<rect width="900" height="420" fill="white"/>',
    '<text x="32" y="32" font-family="Arial" font-size="18">Z4Q public fixture expansion repeatability</text>'
)
$x = 80
foreach ($row in $repeatArray) {
    $audit = [double]$row.audit_max
    $height = if (Test-FmmCliFinite $audit) { [Math]::Min(280, [Math]::Max(3, -[Math]::Log10([Math]::Max($audit,1e-12))*32)) } else { 3 }
    $y = 340 - $height
    $fill = if ($row.repeatability_pass) { '#0072b2' } else { '#d55e00' }
    $svg += ('<rect x="{0}" y="{1}" width="70" height="{2}" fill="{3}"/>' -f $x,$y,$height,$fill)
    $svg += ('<text x="{0}" y="365" font-family="Arial" font-size="12">{1}</text>' -f ($x-4),$row.fixture)
    $x += 190
}
$svg += '</svg>'
Set-Content -LiteralPath $plotPath -Value $svg -Encoding ASCII

[pscustomobject]@{
    decision=$decision; quick=[bool]$Quick; fixtures=($wanted -join ','); rows=$rowArray.Count;
    repeatability_passes=$passCount; manifest=$manifestCsv; matrix=$matrixCsv; repeatability=$repeatCsv
} | Export-Csv -LiteralPath $decisionCsv -NoTypeInformation
@(
    'scope: Z4Q public fixture expansion gate',
    "decision: $decision",
    "quick: $([bool]$Quick)",
    "fixtures: $($wanted -join ',')",
    "repeatability_passes: $passCount",
    'policy: supplement-first evidence; no default change'
) | Set-Content -LiteralPath $summaryPath -Encoding ASCII
@(
    '# TITAN FMM/M2L Scope Z4Q Public Fixture Expansion Report',
    '',
    "Decision: $decision",
    '',
    'Z4Q prepares additional public full-protein fixtures and tests guarded p4 repeatability as supplement-first evidence.',
    '',
    "- Manifest: $manifestCsv",
    "- Matrix: $matrixCsv",
    "- Repeatability: $repeatCsv",
    "- Plot: $plotPath",
    '',
    'Policy:',
    '- Passing Z4Q rows may be added to supplement tables.',
    '- Main manuscript guarded p4 claim remains based on Z4O unless all required Z4Q rows pass.',
    '- Grid remains production default and guarded p4 remains explicit.'
) | Set-Content -LiteralPath $rootReport -Encoding ASCII

if ($decision -notin @('z4q_expansion_pass','z4q_partial_expansion') -and -not $AllowFailedDecision) {
    throw "Scope Z4Q did not pass: $decision"
}
Write-Host "Scope Z4Q public fixture expansion decision: $decision"
