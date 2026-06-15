param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z1a'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 4
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')
. (Join-Path $PSScriptRoot 'pbe_fagbi_validation_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fagbi_scope_z1a_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fagbi_scope_z1a_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fagbi_scope_z1a_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fagbi_scope_z1a_accuracy_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FAGBI_SCOPE_Z1A_PARITY_VALIDATION_REPORT.txt'
$TemperatureK = 298.15
$InternalDielectric = 1.0
$ExternalDielectric = 80.0
$SaltM = 0.15
$OpeningTheta = 0.6
$IluNeighborCap = 64

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope Z1A work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) { Remove-Item -LiteralPath $workRootFull -Recurse -Force }
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

function Get-Z1aProp([object]$Object, [string]$Name, [object]$Default = '') {
    if ($null -eq $Object) { return $Default }
    if ($Object.PSObject.Properties.Name -contains $Name) {
        $value = $Object.$Name
        if ($null -ne $value -and -not [string]::IsNullOrWhiteSpace([string]$value)) { return $value }
    }
    return $Default
}

function F([double]$Value) { return $Value.ToString('G17', $script:FmmCliInvariantCulture) }

function Write-Z1aBornPqr([string]$Path, [double]$Charge, [double]$RadiusA) {
    Write-FagbiBornPqr -Path $Path -Charge $Charge -RadiusA $RadiusA
}

function Write-Z1aKirkwoodPqr([string]$Path, [double]$SphereRadiusA, [double]$Charge, [double]$OffsetA) {
    Write-FagbiKirkwoodPqr -Path $Path -SphereRadiusA $SphereRadiusA -Charge $Charge -OffsetA $OffsetA
}

function Get-Z1aDebyeKappaInvA([double]$Salt, [double]$TempK) {
    return Get-FagbiDebyeKappaInvA -SaltM $Salt -TemperatureK $TempK
}

function Get-Z1aBornEnergyKJ([double]$Charge, [double]$RadiusA) {
    return Get-FagbiBornEnergyKJ -Charge $Charge -RadiusA $RadiusA -SaltM $SaltM -TemperatureK $TemperatureK -InternalDielectric $InternalDielectric -ExternalDielectric $ExternalDielectric
}

function Get-Z1aLegendreP([int]$L, [double]$X) {
    return Get-FagbiLegendreP -L $L -X $X
}

function Get-Z1aSphericalKNextOverK([int]$L, [double]$X) {
    return Get-FagbiSphericalKNextOverK -L $L -X $X
}

function Get-Z1aReactionCoefficient([int]$L, [double]$SphereRadiusA) {
    return Get-FagbiReactionCoefficient -L $L -SphereRadiusA $SphereRadiusA -SaltM $SaltM -TemperatureK $TemperatureK -InternalDielectric $InternalDielectric -ExternalDielectric $ExternalDielectric
}

function Get-Z1aKirkwoodEnergyKJ([double]$SphereRadiusA, [double]$Charge, [double]$OffsetA) {
    return Get-FagbiKirkwoodEnergyKJ -SphereRadiusA $SphereRadiusA -Charge $Charge -OffsetA $OffsetA -SaltM $SaltM -TemperatureK $TemperatureK -InternalDielectric $InternalDielectric -ExternalDielectric $ExternalDielectric -LMax 160
}

function Get-Z1aAnalyticThreshold([int]$SurfacePoints) {
    if ($SurfacePoints -ge 24) { return 1.0e-2 }
    if ($SurfacePoints -ge 12) { return 2.0e-2 }
    return 5.0e-2
}

function Write-Z1aArgFile([string]$RunDir, [string]$PqrPath, [object]$Backend, [int]$SurfacePoints, [int]$MaxIter, [int]$Restart, [double]$Tolerance, [string]$SphereMeshMode) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdbPath
    $libPath = Convert-FmmCliToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-FmmCliToTitanPath $pdbPath)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        "-pbe=$($Backend.backend)",
        "-pbe-test-pqr=$(Convert-FmmCliToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-probe=0.0',
        "-pbe-salt=$(F $SaltM)",
        "-pbe-temp=$(F $TemperatureK)",
        "-pbe-int-diel=$(F $InternalDielectric)",
        "-pbe-ext-diel=$(F $ExternalDielectric)",
        '-pbe-timing=detail',
        '-pbe-map=off',
        '-pbe-margin=2.0',
        '-pbe-scale=0.50',
        '-pbe-max-grid=800000',
        "-pbe-max-iter=$MaxIter",
        "-pbe-rms=$(F $Tolerance)",
        "-pbe-dphi=$(F $Tolerance)",
        "-nt=$Threads"
    )
    if ($Backend.backend -eq 'tabi') {
        $argLines += @(
            "-pbe-tabi-sum=$($Backend.mode)",
            '-pbe-tabi-formulation=jcp2013-treecode',
            '-pbe-tabi-energy-projection=boundary',
            "-pbe-tabi-fmm-m2l=$($Backend.fmm_mode)",
            "-pbe-tabi-surface-points=$SurfacePoints",
            '-pbe-tabi-tree-leaf-size=64',
            "-pbe-tabi-opening-theta=$(F $OpeningTheta)",
            '-pbe-tabi-near-integration=centroid',
            '-pbe-tabi-quad-buffer=off',
            "-pbe-tabi-precond=$($Backend.preconditioner)",
            '-pbe-tabi-precond-block-panel-cap=32',
            '-pbe-tabi-precond-overlap-panel-cap=64',
            "-pbe-tabi-precond-ilu-neighbor-cap=$IluNeighborCap",
            '-pbe-tabi-precond-ilu-drop-tol=0.0',
            "-pbe-tabi-gmres-restart=$Restart",
            '-pbe-tabi-true-residual-probe=restart',
            '-pbe-tabi-true-residual-probe-scope=sampled',
            '-pbe-tabi-true-residual-sample-count=512'
        )
        if (-not [string]::IsNullOrWhiteSpace($SphereMeshMode)) { $argLines += "-pbe-tabi-sphere-mesh=$SphereMeshMode" }
        if ($Backend.require_audit) {
            $argLines += '-pbe-tabi-matvec-audit=sampled'
            $argLines += '-pbe-tabi-matvec-audit-sample-count=512'
        } else {
            $argLines += '-pbe-tabi-matvec-audit=off'
        }
    }
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    $badOptions = @($argLines | Where-Object { -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*=' })
    if ($badOptions.Count -ne 0) { throw "Generated arg.in contains non equals-syntax options: $($badOptions -join '; ')" }
    return $argPath
}

function Read-Z1aFirstCsv([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) { return $null }
    $rows = @(Import-Csv -LiteralPath $Path)
    if ($rows.Count -eq 0) { return $null }
    return $rows[0]
}

function New-Z1aEmptyResult([string]$Status) {
    return [pscustomobject]@{ status=$Status; energy_kj_mol=[double]::NaN; reaction_energy_source=''; panels=0; grid_points=0; iterations=0; final_residual=[double]::NaN; true_residual=[double]::NaN; residual_history=''; matvec_seconds_avg=[double]::NaN; total_seconds=[double]::NaN; audit_rel=[double]::NaN; matvec_audit_rows=0; fmm_m2l_mode=''; fmm_near_pairs=0; fmm_far_pairs=0; fmm_finite_outputs=''; preconditioner_effective=''; preconditioner_ilu_neighbor_cap=0; memory_bytes=[double]::NaN }
}

function Read-Z1aCaseResult([string]$RunDir, [object]$Backend) {
    $energy = Read-Z1aFirstCsv (Join-Path $RunDir 'pbe_energy.csv')
    $timing = Read-Z1aFirstCsv (Join-Path $RunDir 'pbe_timing.csv')
    $profile = Read-Z1aFirstCsv (Join-Path $RunDir 'pbe_tabi_profile.csv')
    if ($null -eq $energy) { throw "Missing pbe_energy.csv in $RunDir" }
    if ($Backend.backend -eq 'tabi' -and $null -eq $profile) { throw "Missing pbe_tabi_profile.csv in $RunDir" }
    $fmmMode = [string](Get-Z1aProp $profile 'fmm_m2l_mode' '')
    if ($Backend.expected_fmm -eq 'single-matvec-analytic' -and $fmmMode -ne 'single-matvec-analytic') { throw "Expected FMM/M2L mode single-matvec-analytic, got '$fmmMode'" }
    $auditRows = Convert-FmmCliInt64 (Get-Z1aProp $profile 'matvec_audit_rows' 0)
    $auditRel = Convert-FmmCliDouble (Get-Z1aProp $profile 'tree_vs_direct_matvec_rel' 'NA')
    if ($Backend.require_audit -and ($auditRows -le 0 -or -not (Test-FmmCliFinite $auditRel))) { throw "Missing sampled FMM audit in $RunDir" }
    $energyStatus = [string](Get-Z1aProp $energy 'status' 'unknown')
    $profileStatus = if ($null -eq $profile) { 'grid' } else { [string](Get-Z1aProp $profile 'status' 'unknown') }
    $status = if ($energyStatus -eq 'ok' -and ($Backend.backend -eq 'grid' -or $profileStatus -eq 'ok')) { 'ok' } else { "energy=$energyStatus;profile=$profileStatus" }
    return [pscustomobject]@{ status=$status; energy_kj_mol=Convert-FmmCliDouble (Get-Z1aProp $energy 'E_rxn_kJmol' 'NA'); reaction_energy_source=[string](Get-Z1aProp $energy 'reaction_energy_source' ''); panels=Convert-FmmCliInt64 (Get-Z1aProp $profile 'panels' (Get-Z1aProp $timing 'tabi_surface_panels' 0)); grid_points=Convert-FmmCliInt64 (Get-Z1aProp $energy 'grid_points' 0); iterations=Convert-FmmCliInt64 (Get-Z1aProp $profile 'iterations' 0); final_residual=Convert-FmmCliDouble (Get-Z1aProp $profile 'final_residual' 'NA'); true_residual=Convert-FmmCliDouble (Get-Z1aProp $profile 'true_residual_rel' 'NA'); residual_history=[string](Get-Z1aProp $profile 'residual_history' ''); matvec_seconds_avg=Convert-FmmCliDouble (Get-Z1aProp $profile 'matvec_seconds_avg' 'NA'); total_seconds=Convert-FmmCliDouble (Get-Z1aProp $profile 'total_seconds' 'NA'); audit_rel=$auditRel; matvec_audit_rows=$auditRows; fmm_m2l_mode=$fmmMode; fmm_near_pairs=Convert-FmmCliInt64 (Get-Z1aProp $profile 'fmm_m2l_near_pairs' 0); fmm_far_pairs=Convert-FmmCliInt64 (Get-Z1aProp $profile 'fmm_m2l_far_pairs' 0); fmm_finite_outputs=[string](Get-Z1aProp $profile 'fmm_m2l_finite_outputs' ''); preconditioner_effective=[string](Get-Z1aProp $profile 'preconditioner_effective' ''); preconditioner_ilu_neighbor_cap=Convert-FmmCliInt64 (Get-Z1aProp $profile 'preconditioner_ilu_neighbor_cap' 0); memory_bytes=Convert-FmmCliDouble (Get-Z1aProp $profile 'memory_bytes' 'NA') }
}

function Invoke-Z1aCase([object]$Case, [object]$Backend) {
    $caseSafe = Get-FmmCliSafeName "$($Case.fixture)_$($Case.fixture_class)_sp$($Case.surface_points)_$($Backend.label)_r$($Case.repeat)"
    $runDir = Join-Path $workRootFull "runs\$caseSafe"
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $pqrPath = Join-Path $runDir 'input.pqr'
    if ($Case.fixture_class -like 'analytic*') {
        if ($Case.analytic_type -eq 'born_centered') { Write-Z1aBornPqr -Path $pqrPath -Charge $Case.charge -RadiusA $Case.radius_A } else { Write-Z1aKirkwoodPqr -Path $pqrPath -SphereRadiusA $Case.radius_A -Charge $Case.charge -OffsetA $Case.offset_A }
    } else {
        Copy-Item -LiteralPath $Case.pqr_path -Destination $pqrPath -Force
    }
    $argPath = Write-Z1aArgFile -RunDir $runDir -PqrPath $pqrPath -Backend $Backend -SurfacePoints $Case.surface_points -MaxIter $Case.max_iter -Restart $Case.restart -Tolerance $Case.tolerance -SphereMeshMode $Case.sphere_mesh
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $readError = ''
    $result = $null
    if ($run.exit_code -eq 0 -and -not $run.timed_out) {
        try { $result = Read-Z1aCaseResult -RunDir $runDir -Backend $Backend } catch { $readError = $_.Exception.Message }
    }
    if ($null -eq $result) {
        $status = if ($run.timed_out) { 'timeout' } elseif (-not [string]::IsNullOrWhiteSpace($readError)) { 'read_error' } else { "exit_$($run.exit_code)" }
        $result = New-Z1aEmptyResult -Status $status
    }
    return [pscustomobject]@{ case=$Case; backend=$Backend; run=$run; result=$result; read_error=$readError; run_dir=$runDir }
}

function Get-Z1aFailureClass([object]$Row) {
    if ($Row.status -eq 'skipped_missing_fixture') { return 'skipped_missing_fixture' }
    if ($Row.status -eq 'timeout') { return 'timeout_not_verified' }
    if ($Row.status -eq 'read_error') { return 'blocked_profile_schema' }
    if ($Row.status -ne 'ok') { return 'status_fail' }
    if (-not (Test-FmmCliFinite ([double]$Row.energy_kj_mol))) { return 'nonfinite_energy' }
    if ($Row.backend -eq 'fmm_ilu0') {
        if ($Row.fixture_class -like 'analytic*') {
            if (-not (Test-FmmCliFinite ([double]$Row.relative_error))) { return 'oracle_error_missing' }
            if ([double]$Row.relative_error -gt [double]$Row.acceptance_threshold) { return 'oracle_error_fail' }
        } else {
            if (-not (Test-FmmCliFinite ([double]$Row.reference_energy_kj_mol))) { return 'direct_reference_missing' }
            if ([math]::Abs([double]$Row.energy_delta_kj_mol) -gt [double]$Row.energy_gate_kj_mol) { return 'energy_delta_fail' }
        }
        if ($Row.matvec_audit_rows -le 0 -or -not (Test-FmmCliFinite ([double]$Row.audit_rel)) -or [double]$Row.audit_rel -gt 1.0e-2) { return 'audit_fail' }
    }
    return 'ok'
}

function Write-Z1aSvgPlot([object[]]$Rows, [string]$Path) {
    $points = @($Rows | Where-Object { $_.backend -eq 'fmm_ilu0' -and (Test-FmmCliFinite ([double]$_.relative_error)) -and [double]$_.relative_error -gt 0.0 })
    if ($points.Count -eq 0) {
        Set-Content -LiteralPath $Path -Encoding ASCII -Value '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="520"><text x="20" y="40">No finite FMM relative-error rows</text></svg>'
        return
    }
    $logs = @($points | ForEach-Object { [math]::Log10([double]$_.relative_error) })
    $minLog = [double](@($logs | Measure-Object -Minimum)[0].Minimum)
    $maxLog = [double](@($logs | Measure-Object -Maximum)[0].Maximum)
    if ([math]::Abs($maxLog - $minLog) -lt 1.0e-12) { $maxLog += 1.0; $minLog -= 1.0 }
    $width = 900; $height = 520; $left = 90; $right = 40; $top = 40; $bottom = 110
    $barW = [math]::Max(10.0, ($width - $left - $right) / [math]::Max(1, $points.Count) - 8.0)
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="520" viewBox="0 0 900 520">')
    $svg.Add('<rect width="100%" height="100%" fill="white"/>')
    $svg.Add('<line x1="90" y1="410" x2="860" y2="410" stroke="black"/>')
    $svg.Add('<line x1="90" y1="40" x2="90" y2="410" stroke="black"/>')
    for ($i = 0; $i -lt $points.Count; ++$i) {
        $row = $points[$i]
        $log = [math]::Log10([double]$row.relative_error)
        $x = $left + 8 + $i * ($barW + 8)
        $y = $top + (($height - $top - $bottom) * (($maxLog - $log) / ($maxLog - $minLog)))
        $barH = ($height - $bottom) - $y
        $color = if ([string]$row.failure_class -eq 'ok') { '#009e73' } else { '#d55e00' }
        $svg.Add(('<rect x="{0:F1}" y="{1:F1}" width="{2:F1}" height="{3:F1}" fill="{4}"/>' -f $x,$y,$barW,$barH,$color))
        $label = "$($row.fixture):$($row.surface_points)"
        $svg.Add(('<text x="{0:F1}" y="426" font-family="Arial" font-size="10" transform="rotate(45 {0:F1},426)">{1}</text>' -f $x,$label))
    }
    $svg.Add('<text x="305" y="500" font-family="Arial" font-size="16">FMM/M2L + ILU0 rows</text>')
    $svg.Add('<text x="25" y="275" font-family="Arial" font-size="16" transform="rotate(-90 25,275)">log10 relative error</text>')
    $svg.Add('<text x="90" y="25" font-family="Arial" font-size="18" font-weight="bold">Scope Z1A FAGBI-style parity accuracy</text>')
    $svg.Add('</svg>')
    Set-Content -LiteralPath $Path -Encoding ASCII -Value $svg
}

$surfacePoints = if ($Quick) { @(8) } else { @(8, 12, 24) }
$proteinCap = if ($Quick) { 64 } else { 128 }
$maxIter = if ($Quick) { 180 } else { 240 }
$restart = if ($Quick) { 48 } else { 64 }
$tolerance = 1.0e-6

$backendSpecs = @(
    [pscustomobject]@{ label='grid'; backend='grid'; mode='direct'; fmm_mode='off'; expected_fmm='off'; preconditioner='diagonal'; require_audit=$false },
    [pscustomobject]@{ label='tabi_direct'; backend='tabi'; mode='direct'; fmm_mode='off'; expected_fmm='off'; preconditioner='diagonal'; require_audit=$false },
    [pscustomobject]@{ label='fmm_ilu0'; backend='tabi'; mode='treecode'; fmm_mode='auto'; expected_fmm='single-matvec-analytic'; preconditioner='auto'; require_audit=$true },
    [pscustomobject]@{ label='fmm_block_jacobi'; backend='tabi'; mode='treecode'; fmm_mode='single-matvec-analytic'; expected_fmm='single-matvec-analytic'; preconditioner='fmm-block-jacobi'; require_audit=$true }
)

$caseSpecs = New-Object System.Collections.Generic.List[object]
foreach ($sp in $surfacePoints) {
    $caseSpecs.Add([pscustomobject]@{ fixture='born_centered'; fixture_class='analytic_born'; analytic_type='born_centered'; charge=1.0; radius_A=3.0; offset_A=0.0; surface_points=$sp; cap=0; repeat=1; pqr_path=''; max_iter=$maxIter; restart=$restart; tolerance=$tolerance; sphere_mesh='icosphere'; optional=$false; sampled_atoms=1; oracle_energy_kj_mol=(Get-Z1aBornEnergyKJ -Charge 1.0 -RadiusA 3.0); acceptance_threshold=(Get-Z1aAnalyticThreshold -SurfacePoints $sp); missing_fixture=$false }) | Out-Null
    $offset = 3.0
    $caseSpecs.Add([pscustomobject]@{ fixture='kirkwood_offcenter'; fixture_class='analytic_kirkwood'; analytic_type='kirkwood_offcenter'; charge=1.0; radius_A=10.0; offset_A=$offset; surface_points=$sp; cap=0; repeat=1; pqr_path=''; max_iter=$maxIter; restart=$restart; tolerance=$tolerance; sphere_mesh='icosphere'; optional=$false; sampled_atoms=2; oracle_energy_kj_mol=(Get-Z1aKirkwoodEnergyKJ -SphereRadiusA 10.0 -Charge 1.0 -OffsetA $offset); acceptance_threshold=(Get-Z1aAnalyticThreshold -SurfacePoints $sp); missing_fixture=$false }) | Out-Null
}

$sampleRoot = Join-Path $workRootFull 'samples'
New-Item -ItemType Directory -Path $sampleRoot -Force | Out-Null
$fixtureNames = if ($Quick) { @('barnase') } else { @('barnase', 'research-medium', 'research-alt') }
foreach ($fixtureName in $fixtureNames) {
    $sourceDir = Join-Path $sampleRoot "source_$fixtureName"
    New-Item -ItemType Directory -Path $sourceDir -Force | Out-Null
    $source = Resolve-FmmCliFixture -Fixture $fixtureName -RepoRoot $RepoRoot -RunDir $sourceDir
    if ([string]::IsNullOrWhiteSpace($source) -or -not (Test-Path -LiteralPath $source)) {
        $caseSpecs.Add([pscustomobject]@{ fixture=$fixtureName; fixture_class='protein'; analytic_type='none'; charge=0.0; radius_A=0.0; offset_A=0.0; surface_points=8; cap=$proteinCap; repeat=1; pqr_path=''; max_iter=$maxIter; restart=$restart; tolerance=$tolerance; sphere_mesh=''; optional=$false; sampled_atoms=0; oracle_energy_kj_mol=[double]::NaN; acceptance_threshold=[double]::NaN; missing_fixture=$true }) | Out-Null
    } else {
        $samplePath = Join-Path $sampleRoot ($fixtureName + "_cap$proteinCap.pqr")
        $sampleInfo = Write-FmmCliSampledPqr -InputPath $source -OutputPath $samplePath -Cap $proteinCap
        $caseSpecs.Add([pscustomobject]@{ fixture=$fixtureName; fixture_class='protein'; analytic_type='none'; charge=0.0; radius_A=0.0; offset_A=0.0; surface_points=8; cap=$proteinCap; repeat=1; pqr_path=$samplePath; max_iter=$maxIter; restart=$restart; tolerance=$tolerance; sphere_mesh=''; optional=$false; sampled_atoms=$sampleInfo.sampled_atoms; oracle_energy_kj_mol=[double]::NaN; acceptance_threshold=[double]::NaN; missing_fixture=$false }) | Out-Null
    }
}

$oneA63 = Join-Path $RepoRoot 'test\pbe_inputs_local\1a63.pqr'
if (Test-Path -LiteralPath $oneA63) {
    $samplePath = Join-Path $sampleRoot "1A63_cap$proteinCap.pqr"
    $sampleInfo = Write-FmmCliSampledPqr -InputPath $oneA63 -OutputPath $samplePath -Cap $proteinCap
    $caseSpecs.Add([pscustomobject]@{ fixture='1A63'; fixture_class='protein_optional'; analytic_type='none'; charge=0.0; radius_A=0.0; offset_A=0.0; surface_points=8; cap=$proteinCap; repeat=1; pqr_path=$samplePath; max_iter=$maxIter; restart=$restart; tolerance=$tolerance; sphere_mesh=''; optional=$true; sampled_atoms=$sampleInfo.sampled_atoms; oracle_energy_kj_mol=[double]::NaN; acceptance_threshold=[double]::NaN; missing_fixture=$false }) | Out-Null
} else {
    $caseSpecs.Add([pscustomobject]@{ fixture='1A63'; fixture_class='protein_optional'; analytic_type='none'; charge=0.0; radius_A=0.0; offset_A=0.0; surface_points=8; cap=$proteinCap; repeat=1; pqr_path=''; max_iter=$maxIter; restart=$restart; tolerance=$tolerance; sphere_mesh=''; optional=$true; sampled_atoms=0; oracle_energy_kj_mol=[double]::NaN; acceptance_threshold=[double]::NaN; missing_fixture=$true }) | Out-Null
}

$records = New-Object System.Collections.Generic.List[object]
foreach ($case in $caseSpecs) {
    if ($case.missing_fixture) {
        foreach ($backend in $backendSpecs) {
            $records.Add([pscustomobject]@{ case=$case; backend=$backend; run=$null; result=(New-Z1aEmptyResult -Status 'skipped_missing_fixture'); read_error=''; run_dir='' }) | Out-Null
        }
    } else {
        foreach ($backend in $backendSpecs) { $records.Add((Invoke-Z1aCase -Case $case -Backend $backend)) | Out-Null }
    }
}

$baseRows = foreach ($record in $records) {
    $runStatus = if ($null -eq $record.run) { 'skipped_missing_fixture' } elseif ($record.run.timed_out) { 'timeout' } elseif (-not [string]::IsNullOrWhiteSpace($record.read_error)) { 'read_error' } else { $record.result.status }
    $case = $record.case
    $backend = $record.backend
    $referenceEnergy = if ($case.fixture_class -like 'analytic*') { [double]$case.oracle_energy_kj_mol } else { [double]::NaN }
    $energy = [double]$record.result.energy_kj_mol
    $relError = if ((Test-FmmCliFinite $referenceEnergy) -and (Test-FmmCliFinite $energy)) { [math]::Abs($energy - $referenceEnergy) / [math]::Max(1.0e-12, [math]::Abs($referenceEnergy)) } else { [double]::NaN }
    [pscustomobject]@{ fixture=$case.fixture; fixture_class=$case.fixture_class; backend=$backend.label; backend_engine=$backend.backend; preconditioner=$backend.preconditioner; cap=$case.cap; surface_points=$case.surface_points; repeat=$case.repeat; sampled_atoms=$case.sampled_atoms; panels=$record.result.panels; grid_points=$record.result.grid_points; status=$runStatus; read_error=$record.read_error; exit_code=if ($null -eq $record.run) { 0 } else { $record.run.exit_code }; timed_out=if ($null -eq $record.run) { $false } else { $record.run.timed_out }; energy_kj_mol=$energy; oracle_energy_kj_mol=$case.oracle_energy_kj_mol; reference_energy_kj_mol=$referenceEnergy; relative_error=$relError; acceptance_threshold=$case.acceptance_threshold; energy_delta_kj_mol=[double]::NaN; energy_gate_kj_mol=[double]::NaN; true_residual=$record.result.true_residual; final_residual=$record.result.final_residual; audit_rel=$record.result.audit_rel; matvec_audit_rows=$record.result.matvec_audit_rows; iterations=$record.result.iterations; matvec_seconds_avg=$record.result.matvec_seconds_avg; wall_seconds=if ($null -eq $record.run) { [double]::NaN } else { $record.run.wall_seconds }; peak_rss_mb=if ($null -eq $record.run) { [double]::NaN } else { $record.run.peak_working_set_mb }; fmm_m2l_mode=$record.result.fmm_m2l_mode; fmm_near_pairs=$record.result.fmm_near_pairs; fmm_far_pairs=$record.result.fmm_far_pairs; fmm_finite_outputs=$record.result.fmm_finite_outputs; preconditioner_effective=$record.result.preconditioner_effective; preconditioner_ilu_neighbor_cap=$record.result.preconditioner_ilu_neighbor_cap; reaction_energy_source=$record.result.reaction_energy_source; residual_history=$record.result.residual_history; failure_class='pending' }
}

$rows = foreach ($row in $baseRows) {
    $direct = @($baseRows | Where-Object { $_.fixture -eq $row.fixture -and $_.fixture_class -eq $row.fixture_class -and $_.surface_points -eq $row.surface_points -and $_.cap -eq $row.cap -and $_.repeat -eq $row.repeat -and $_.backend -eq 'tabi_direct' }) | Select-Object -First 1
    if (-not ($row.fixture_class -like 'analytic*') -and $null -ne $direct -and (Test-FmmCliFinite ([double]$direct.energy_kj_mol)) -and (Test-FmmCliFinite ([double]$row.energy_kj_mol)) -and $row.status -ne 'skipped_missing_fixture') {
        $row.reference_energy_kj_mol = [double]$direct.energy_kj_mol
        $row.energy_delta_kj_mol = [double]$row.energy_kj_mol - [double]$direct.energy_kj_mol
        $row.relative_error = [math]::Abs([double]$row.energy_delta_kj_mol) / [math]::Max(1.0e-12, [math]::Abs([double]$direct.energy_kj_mol))
        $row.energy_gate_kj_mol = [math]::Max(2.0, 0.05 * [math]::Max(1.0, [math]::Abs([double]$direct.energy_kj_mol)))
        $row.acceptance_threshold = 5.0e-2
    }
    $row.failure_class = Get-Z1aFailureClass -Row $row
    $row
}

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
Write-Z1aSvgPlot -Rows $rows -Path $PlotPath

$requiredRows = @($rows | Where-Object { $_.fixture -ne '1A63' })
$fmmRows = @($requiredRows | Where-Object { $_.backend -eq 'fmm_ilu0' })
$blockedMissing = @($requiredRows | Where-Object { $_.failure_class -eq 'skipped_missing_fixture' }).Count -gt 0
$profileBlocked = @($requiredRows | Where-Object { $_.failure_class -eq 'blocked_profile_schema' }).Count -gt 0
$accuracyFails = @($fmmRows | Where-Object { $_.failure_class -in @('oracle_error_fail','energy_delta_fail','audit_fail','oracle_error_missing','direct_reference_missing','nonfinite_energy') })
$runtimeFails = New-Object System.Collections.Generic.List[object]
foreach ($fmm in @($fmmRows | Where-Object { $_.fixture_class -like 'protein*' -and $_.failure_class -eq 'ok' })) {
    $direct = @($rows | Where-Object { $_.fixture -eq $fmm.fixture -and $_.backend -eq 'tabi_direct' -and $_.cap -eq $fmm.cap -and $_.surface_points -eq $fmm.surface_points }) | Select-Object -First 1
    if ($null -ne $direct) {
        $directTime = if (Test-FmmCliFinite ([double]$direct.matvec_seconds_avg)) { [double]$direct.matvec_seconds_avg } else { [double]$direct.wall_seconds }
        $fmmTime = if (Test-FmmCliFinite ([double]$fmm.matvec_seconds_avg)) { [double]$fmm.matvec_seconds_avg } else { [double]$fmm.wall_seconds }
        if ((Test-FmmCliFinite $directTime) -and (Test-FmmCliFinite $fmmTime) -and $fmmTime -gt $directTime) { $runtimeFails.Add($fmm) | Out-Null }
    }
}

$z1aDecision = if ($blockedMissing) { 'blocked_missing_fixture' } elseif ($profileBlocked) { 'blocked_profile_schema' } elseif ($accuracyFails.Count -gt 0) { 'z1a_accuracy_gap' } elseif ($runtimeFails.Count -gt 0) { 'z1a_runtime_gap' } elseif (@($fmmRows | Where-Object { $_.failure_class -ne 'ok' }).Count -gt 0) { 'z1a_harness_ready_with_gaps' } else { 'z1a_parity_pass' }

[pscustomobject]@{ scope_z1a_decision=$z1aDecision; quick=$Quick.IsPresent; rows=$rows.Count; fmm_rows=$fmmRows.Count; accuracy_fail_rows=$accuracyFails.Count; runtime_fail_rows=$runtimeFails.Count; optional_1a63_status=(@($rows | Where-Object { $_.fixture -eq '1A63' }) | Select-Object -First 1).status; matrix=$MatrixCsv; plot=$PlotPath } | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation

$summary = @(
    'TITAN FAGBI Scope Z1A Parity Validation Summary',
    "scope_z1a_decision: $z1aDecision",
    "quick: $($Quick.IsPresent)",
    "rows: $($rows.Count)",
    "analytic_surface_points: $($surfacePoints -join ',')",
    "protein_cap: $proteinCap",
    'backends: grid,tabi_direct,fmm_ilu0,fmm_block_jacobi',
    'fmm_default_note: fmm_ilu0 rows use -pbe-tabi-fmm-m2l=auto and -pbe-tabi-precond=auto to test Scope Y TABI treecode defaults.',
    "optional_1a63_status: $((@($rows | Where-Object { $_.fixture -eq '1A63' }) | Select-Object -First 1).status)",
    "matrix: $MatrixCsv",
    "plot: $PlotPath"
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = New-Object System.Collections.Generic.List[string]
$report.Add('TITAN FAGBI Scope Z1A Parity Validation Report')
$report.Add("Decision: $z1aDecision")
$report.Add("Generated: $(Get-Date -Format s)")
$report.Add('')
$report.Add('Purpose')
$report.Add('Measure TITAN TABI treecode against FAGBI-style validation axes before implementing Duffy/cubature or adaptive FMM order work.')
$report.Add('')
$report.Add('Key gates')
$report.Add('- Analytic FMM rows compare against Born/Kirkwood oracle energy.')
$report.Add('- Protein FMM rows compare against TABI direct energy and sampled matvec audit.')
$report.Add('- Grid remains the production default; this report does not promote TABI over grid.')
$report.Add('')
$report.Add('Decision counts')
$report.Add("accuracy_fail_rows: $($accuracyFails.Count)")
$report.Add("runtime_fail_rows: $($runtimeFails.Count)")
$report.Add("blocked_missing_fixture: $blockedMissing")
$report.Add("blocked_profile_schema: $profileBlocked")
$report.Add('')
$report.Add('FMM/ILU0 rows')
foreach ($row in $fmmRows) {
    $report.Add(('{0} class={1} sp={2} cap={3} status={4} rel={5:G6} audit={6:G6} iters={7} failure={8}' -f $row.fixture, $row.fixture_class, $row.surface_points, $row.cap, $row.status, [double]$row.relative_error, [double]$row.audit_rel, $row.iterations, $row.failure_class))
}
$report.Add('')
$report.Add("Matrix: $MatrixCsv")
$report.Add("Decision: $DecisionCsv")
$report.Add("Plot: $PlotPath")
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host "Scope Z1A decision: $z1aDecision"
Write-Host "Matrix: $MatrixCsv"
