param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 180,
    [int]$Threads = 4,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_y')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$MatrixCsv = Join-Path $WorkRoot 'fmm_m2l_scope_y_matrix.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_y_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Y_TREECODE_DEFAULT_REPORT.txt'

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$workRootFull = [System.IO.Path]::GetFullPath($WorkRoot)
$allowedTmpRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot 'tmp')).TrimEnd('\') + '\'
if (-not $workRootFull.StartsWith($allowedTmpRoot,
        [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to clean Scope Y work root outside test/tmp: $workRootFull"
}
if (Test-Path -LiteralPath $workRootFull) {
    Remove-Item -LiteralPath $workRootFull -Recurse -Force
}
New-Item -ItemType Directory -Path $workRootFull -Force | Out-Null

$sampleDir = Join-Path $workRootFull 'samples'
New-Item -ItemType Directory -Path $sampleDir -Force | Out-Null
$fixture = Resolve-FmmCliFixture -Fixture 'barnase' -RepoRoot $RepoRoot `
    -RunDir $sampleDir
if ([string]::IsNullOrWhiteSpace($fixture)) {
    $fixture = Resolve-FmmCliFixture -Fixture 'research-medium' -RepoRoot $RepoRoot
}
$pqrPath = Join-Path $sampleDir 'barnase_cap128.pqr'
$sampleInfo = Write-FmmCliSampledPqr -InputPath $fixture -OutputPath $pqrPath -Cap 128
if ($sampleInfo.sampled_atoms -ne 128) {
    throw "Scope Y expected Barnase cap-128 sample, got $($sampleInfo.sampled_atoms)"
}

function Write-ScopeYArgFile(
    [string]$RunDir,
    [string]$PqrPath,
    [ValidateSet('direct','treecode')][string]$Mode,
    [string]$Formulation,
    [string]$FmmMode,
    [string]$Preconditioner,
    [object]$IluNeighborCap
) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-FmmCliDummyPdb $pdbPath
    $libPath = Convert-FmmCliToTitanPath(
        (Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argLines = @(
        "-f=$(Convert-FmmCliToTitanPath $pdbPath)",
        '-1pdb=true',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-FmmCliToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-max-iter=1',
        '-pbe-rms=1e-6',
        '-pbe-dphi=1e-6',
        '-pbe-int-diel=1',
        '-pbe-ext-diel=80',
        '-pbe-salt=0.15',
        '-pbe-temp=298.15',
        '-pbe-timing=detail',
        "-nt=$Threads",
        "-pbe-tabi-sum=$Mode",
        '-pbe-tabi-surface-points=8',
        '-pbe-tabi-tree-leaf-size=64',
        '-pbe-tabi-opening-theta=0.6',
        '-pbe-tabi-near-integration=centroid',
        '-pbe-tabi-quad-buffer=off',
        '-pbe-tabi-gmres-restart=64',
        '-pbe-tabi-matvec-audit=off',
        '-pbe-tabi-true-residual-probe=off'
    )
    if (-not [string]::IsNullOrWhiteSpace($Formulation)) {
        $argLines += "-pbe-tabi-formulation=$Formulation"
    }
    if (-not [string]::IsNullOrWhiteSpace($FmmMode)) {
        $argLines += "-pbe-tabi-fmm-m2l=$FmmMode"
    }
    if (-not [string]::IsNullOrWhiteSpace($Preconditioner)) {
        $argLines += "-pbe-tabi-precond=$Preconditioner"
    }
    if ($null -ne $IluNeighborCap) {
        $argLines += "-pbe-tabi-precond-ilu-neighbor-cap=$([int]$IluNeighborCap)"
    }

    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Value $argLines -Encoding ASCII
    $badOptions = @(Get-Content -LiteralPath $argPath | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_) -and $_ -notmatch '^-.*='
    })
    if ($badOptions.Count -ne 0) {
        throw "Generated arg.in contains non equals-syntax options: $($badOptions -join '; ')"
    }
    return $argPath
}

function Read-ScopeYProfile([string]$RunDir) {
    $profilePath = Join-Path $RunDir 'pbe_tabi_profile.csv'
    if (-not (Test-Path -LiteralPath $profilePath)) {
        throw "Missing pbe_tabi_profile.csv in $RunDir"
    }
    $profile = @(Import-Csv -LiteralPath $profilePath)[0]
    return [pscustomobject]@{
        fmm_m2l_mode = [string]$profile.fmm_m2l_mode
        fmm_far_pairs = Convert-FmmCliInt64 $profile.fmm_m2l_far_pairs
        preconditioner_requested = [string]$profile.preconditioner_requested
        preconditioner_effective = [string]$profile.preconditioner_effective
        ilu_neighbor_cap = Convert-FmmCliInt64 $profile.preconditioner_ilu_neighbor_cap
        sparse_nnz_blocks = Convert-FmmCliInt64 $profile.preconditioner_sparse_nnz_blocks
    }
}

function Invoke-ScopeYCase(
    [string]$Name,
    [ValidateSet('direct','treecode')][string]$Mode,
    [string]$Formulation,
    [string]$FmmMode,
    [string]$Preconditioner,
    [object]$IluNeighborCap
) {
    $runDir = Join-Path $workRootFull $Name
    New-FmmCliCleanDir -Path $runDir -WorkRoot $workRootFull
    $argPath = Write-ScopeYArgFile -RunDir $runDir -PqrPath $pqrPath `
        -Mode $Mode -Formulation $Formulation -FmmMode $FmmMode `
        -Preconditioner $Preconditioner -IluNeighborCap $IluNeighborCap
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    if ($run.exit_code -ne 0) {
        throw "Scope Y case $Name failed with exit code $($run.exit_code)"
    }
    $profile = Read-ScopeYProfile -RunDir $runDir
    return [pscustomobject]@{
        condition = $Name
        mode = $Mode
        arg_has_explicit_formulation = -not [string]::IsNullOrWhiteSpace($Formulation)
        arg_has_explicit_fmm = -not [string]::IsNullOrWhiteSpace($FmmMode)
        arg_has_explicit_precond = -not [string]::IsNullOrWhiteSpace($Preconditioner)
        wall_seconds = $run.wall_seconds
        peak_working_set_mb = $run.peak_working_set_mb
        fmm_m2l_mode = $profile.fmm_m2l_mode
        fmm_far_pairs = $profile.fmm_far_pairs
        preconditioner_requested = $profile.preconditioner_requested
        preconditioner_effective = $profile.preconditioner_effective
        ilu_neighbor_cap = $profile.ilu_neighbor_cap
        sparse_nnz_blocks = $profile.sparse_nnz_blocks
    }
}

$rows = @()
$rows += Invoke-ScopeYCase -Name 'explicit_fmm_ilu_reference' -Mode treecode `
    -Formulation 'jcp2013-treecode' -FmmMode 'single-matvec-analytic' `
    -Preconditioner 'fmm-nearfield-ilu0' -IluNeighborCap 64
$rows += Invoke-ScopeYCase -Name 'default_treecode_auto_fmm_ilu' -Mode treecode `
    -Formulation '' -FmmMode '' -Preconditioner '' -IluNeighborCap $null
$rows += Invoke-ScopeYCase -Name 'legacy_treecode_escape_off' -Mode treecode `
    -Formulation '' -FmmMode 'off' -Preconditioner '' -IluNeighborCap $null
$rows += Invoke-ScopeYCase -Name 'explicit_diagonal_keeps_override' -Mode treecode `
    -Formulation '' -FmmMode '' -Preconditioner 'diagonal' -IluNeighborCap $null
$rows += Invoke-ScopeYCase -Name 'direct_tabi_no_auto_fmm' -Mode direct `
    -Formulation '' -FmmMode '' -Preconditioner '' -IluNeighborCap $null

$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation

$byName = @{}
foreach ($row in $rows) { $byName[$row.condition] = $row }
$default = $byName['default_treecode_auto_fmm_ilu']
$explicit = $byName['explicit_fmm_ilu_reference']
$legacy = $byName['legacy_treecode_escape_off']
$diagonal = $byName['explicit_diagonal_keeps_override']
$direct = $byName['direct_tabi_no_auto_fmm']

if ($default.fmm_m2l_mode -ne 'single-matvec-analytic') {
    throw "Default TABI treecode did not auto-enable FMM/M2L: $($default.fmm_m2l_mode)"
}
if ($default.arg_has_explicit_formulation) {
    throw 'Default TABI treecode case should test omitted formulation resolution.'
}
if ($default.preconditioner_effective -ne 'fmm_nearfield_ilu0') {
    throw "Default TABI treecode did not auto-enable ILU0: $($default.preconditioner_effective)"
}
if ($default.ilu_neighbor_cap -ne 64) {
    throw "Default TABI treecode ILU neighbor cap drifted: $($default.ilu_neighbor_cap)"
}
if ($default.fmm_m2l_mode -ne $explicit.fmm_m2l_mode -or
    $default.preconditioner_effective -ne $explicit.preconditioner_effective -or
    $default.ilu_neighbor_cap -ne $explicit.ilu_neighbor_cap) {
    throw 'Default TABI treecode does not match explicit FMM+ILU0 reference behavior.'
}
if ($legacy.fmm_m2l_mode -ne 'off') {
    throw "Explicit -pbe-tabi-fmm-m2l=off did not keep legacy treecode: $($legacy.fmm_m2l_mode)"
}
if ($diagonal.fmm_m2l_mode -ne 'single-matvec-analytic' -or
    $diagonal.preconditioner_effective -ne 'diagonal') {
    throw "Explicit diagonal override failed: mode=$($diagonal.fmm_m2l_mode) precond=$($diagonal.preconditioner_effective)"
}
if ($direct.fmm_m2l_mode -ne 'off') {
    throw "Direct TABI unexpectedly auto-enabled FMM/M2L: $($direct.fmm_m2l_mode)"
}

$summary = @(
    'scope: Y',
    'decision: controlled_tabi_treecode_default',
    'grid_default: unchanged',
    'treecode_default_engine: fmm_m2l_single_matvec_analytic',
    'treecode_default_preconditioner: fmm_nearfield_ilu0',
    'treecode_default_ilu_neighbor_cap: 64',
    'legacy_escape: -pbe-tabi-fmm-m2l=off'
)
Set-Content -LiteralPath $SummaryPath -Value $summary -Encoding ASCII

$report = @(
    '# TITAN FMM/M2L Scope Y Treecode Default Report',
    '',
    'Scope Y promotes analytical FMM/M2L plus sparse near-field ILU0 to the default engine only when users explicitly choose TABI treecode.',
    '',
    'Grid remains the production PBE default. The legacy TABI treecode remains available with `-pbe-tabi-fmm-m2l=off`.',
    '',
    '## Smoke Matrix',
    ''
)
foreach ($row in $rows) {
    $report += ('- {0}: mode={1}, preconditioner={2}, ilu_cap={3}, far_pairs={4}' -f `
        $row.condition, $row.fmm_m2l_mode, $row.preconditioner_effective, `
        $row.ilu_neighbor_cap, $row.fmm_far_pairs)
}
Set-Content -LiteralPath $ReportPath -Value $report -Encoding ASCII

Write-Host 'PASS pbe_fmm_m2l_scope_y_treecode_default_probe'
