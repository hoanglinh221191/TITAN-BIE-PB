param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_preconditioner_controlled_probe'),
    [int]$CaseTimeoutSeconds = 45,
    [int]$SurfacePoints = 24,
    [int]$Threads = 2
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture
$RepoRoot = Split-Path -Parent $PSScriptRoot

. (Join-Path $PSScriptRoot 'pbe_tabi_triage_helpers.ps1')

function Convert-ToTitanPath([string]$Path) { return ($Path -replace '\\', '/') }

function Get-RowValue($Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI controlled probe root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Write-DummyPdb([string]$Path) {
    Set-Content -LiteralPath $Path -Encoding ASCII -Value @(
        'CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1           1',
        'ATOM      1  N   ALA A   1      -2.000   0.000   0.000  1.00  0.00           N',
        'ATOM      2  CA  ALA A   1      -0.600   0.100   0.000  1.00  0.00           C',
        'ATOM      3  C   ALA A   1       0.000   1.450   0.000  1.00  0.00           C',
        'ATOM      4  O   ALA A   1      -0.600   2.430   0.000  1.00  0.00           O',
        'END'
    )
}

function New-PqrAtomLine([int]$Serial, [string]$Atom, [string]$Res, [string]$Chain,
                         [int]$Resnr, [double]$X, [double]$Y, [double]$Z,
                         [double]$Charge, [double]$Radius) {
    $chargeText = $Charge.ToString('F4', $InvariantCulture)
    $radiusText = $Radius.ToString('F4', $InvariantCulture)
    return ("ATOM  {0,5} {1,-4} {2,3} {3}{4,4}    {5,8:F3}{6,8:F3}{7,8:F3}{8,8}{9,8}" -f
        $Serial, $Atom, $Res, $Chain, $Resnr, $X, $Y, $Z, $chargeText, $radiusText)
}

function Write-SyntheticProteinPqr([string]$Path) {
    $atoms = @(
        @(1,  'N',  'ALA', 'A', 1, -4.000,  0.000,  0.000,  0.300, 1.550),
        @(2,  'CA', 'ALA', 'A', 1, -2.650,  0.120,  0.050, -0.100, 1.700),
        @(3,  'CB', 'ALA', 'A', 1, -2.200, -1.180,  0.600,  0.050, 1.700),
        @(4,  'C',  'ALA', 'A', 1, -1.750,  1.150, -0.200,  0.550, 1.700),
        @(5,  'O',  'ALA', 'A', 1, -2.100,  2.250, -0.450, -0.550, 1.520),
        @(6,  'N',  'LYS', 'A', 2, -0.500,  0.850, -0.100, -0.300, 1.550),
        @(7,  'CA', 'LYS', 'A', 2,  0.450,  1.850, -0.350,  0.100, 1.700),
        @(8,  'CB', 'LYS', 'A', 2,  1.800,  1.150, -0.100,  0.100, 1.700),
        @(9,  'CG', 'LYS', 'A', 2,  2.900,  2.100,  0.250,  0.100, 1.700),
        @(10, 'NZ', 'LYS', 'A', 2,  4.100,  1.450,  0.500,  1.000, 1.550),
        @(11, 'C',  'LYS', 'A', 2,  0.700,  2.950,  0.650,  0.550, 1.700),
        @(12, 'O',  'LYS', 'A', 2, -0.150,  3.750,  0.950, -0.550, 1.520),
        @(13, 'N',  'ASP', 'A', 3,  1.900,  3.000,  1.050, -0.300, 1.550),
        @(14, 'CA', 'ASP', 'A', 3,  2.350,  4.100,  1.900,  0.100, 1.700),
        @(15, 'CB', 'ASP', 'A', 3,  3.850,  4.000,  2.050, -0.200, 1.700),
        @(16, 'CG', 'ASP', 'A', 3,  4.450,  5.150,  2.850,  0.600, 1.700),
        @(17, 'OD1','ASP', 'A', 3,  5.650,  5.180,  3.050, -0.800, 1.520),
        @(18, 'OD2','ASP', 'A', 3,  3.750,  6.050,  3.250, -0.800, 1.520)
    )
    $lines = foreach ($a in $atoms) {
        New-PqrAtomLine -Serial $a[0] -Atom $a[1] -Res $a[2] -Chain $a[3] -Resnr $a[4] `
            -X $a[5] -Y $a[6] -Z $a[7] -Charge $a[8] -Radius $a[9]
    }
    $lines += 'END'
    Set-Content -LiteralPath $Path -Encoding ASCII -Value $lines
}

function New-ArgFile([string]$RunDir, [string]$CaseName, [string]$SumMode,
                     [string]$Formulation, [string]$Preconditioner,
                     [string]$BlockNear,
                     [string]$PqrPath) {
    $pdbPath = Join-Path $RunDir 'driver.pdb'
    Write-DummyPdb $pdbPath
    $libPath = Convert-ToTitanPath ((Resolve-Path -LiteralPath (Join-Path $RepoRoot 'lib_ff')).Path.TrimEnd('\') + '\')
    $argPath = Join-Path $RunDir 'arg.in'
    Set-Content -LiteralPath $argPath -Encoding ASCII -Value @(
        "-f=$(Convert-ToTitanPath $pdbPath)",
        '-1pdb',
        '-ff=AMBER99SB',
        '-g=1&1',
        "-path-a=$libPath",
        "-path-n=$libPath",
        '-pbe=tabi',
        "-pbe-test-pqr=$(Convert-ToTitanPath $PqrPath)",
        '-pbe-surface=vdw',
        '-pbe-probe=0.0',
        '-pbe-salt=0.15',
        '-pbe-max-iter=40',
        '-pbe-rms=1e-4',
        '-pbe-dphi=1e-4',
        '-pbe-timing=detail',
        "-pbe-tabi-sum=$SumMode",
        "-pbe-tabi-formulation=$Formulation",
        "-pbe-tabi-precond=$Preconditioner",
        "-pbe-tabi-precond-block-near=$BlockNear",
        '-pbe-tabi-precond-block-panel-cap=32',
        '-pbe-tabi-precond-memory-budget-mb=256',
        '-pbe-tabi-near-integration=duffy_adaptive',
        '-pbe-tabi-quad-buffer=adaptive',
        "-pbe-tabi-surface-points=$SurfacePoints",
        '-pbe-tabi-gmres-restart=8',
        '-pbe-tabi-progress=detail',
        '-pbe-tabi-progress-interval=1',
        '-pbe-tabi-pre-solver-progress=on',
        "-nt=$Threads"
    )
    return $argPath
}

function Invoke-TitanCase([string]$RunDir, [string]$ArgPath) {
    $stdout = Join-Path $RunDir 'stdout.log'
    $stderr = Join-Path $RunDir 'stderr.log'
    Push-Location $RunDir
    try {
        $process = Start-Process -FilePath $TitanExe -ArgumentList @('-a', $ArgPath) `
            -NoNewWindow -PassThru -RedirectStandardOutput $stdout -RedirectStandardError $stderr
        if (-not $process.WaitForExit($CaseTimeoutSeconds * 1000)) {
            Stop-Process -Id $process.Id -Force
            $process.WaitForExit()
            return 'timeout'
        }
        $process.Refresh()
        if ($null -ne $process.ExitCode -and $process.ExitCode -ne 0) {
            return "exit_$($process.ExitCode)"
        }
        return 'ok'
    } finally {
        Pop-Location
    }
}

function Read-CaseMetrics([string]$RunDir, [string]$CaseName, [string]$SumMode,
                          [string]$Formulation, [string]$Preconditioner,
                          [string]$BlockNear,
                          [string]$Status) {
    $timingPath = Join-Path $RunDir 'pbe_timing.csv'
    $timing = $null
    if (Test-Path -LiteralPath $timingPath) {
        $timingRows = @(Import-Csv -LiteralPath $timingPath)
        if ($timingRows.Count -gt 0) { $timing = $timingRows[0] }
    }
    $energyPath = Join-Path $RunDir 'pbe_energy.csv'
    $energy = $null
    if (Test-Path -LiteralPath $energyPath) {
        $energyRows = @(Import-Csv -LiteralPath $energyPath)
        if ($energyRows.Count -gt 0) { $energy = $energyRows[0] }
    }
    $progressRows = @(Get-TabiProgressRows -RunDir $RunDir)
    $lastProgress = if ($progressRows.Count -gt 0) { $progressRows[-1] } else { $null }
    $fallback = [PSCustomObject]@{
        mode = $CaseName
        status = $Status
        total_seconds = Get-RowValue $timing 'total_seconds'
        matvec_calls = Get-RowValue $timing 'tabi_matvec_calls' (Get-RowValue $lastProgress 'matvec_calls')
        panels = Get-RowValue $timing 'tabi_panels' (Get-RowValue $lastProgress 'panels')
        preconditioner_effective = Get-RowValue $timing 'tabi_preconditioner_effective' (Get-RowValue $lastProgress 'preconditioner_effective')
        preconditioner_requested = $Preconditioner
        best_residual = Get-RowValue $timing 'tabi_best_residual'
        final_residual = Get-RowValue $timing 'tabi_final_residual' (Get-RowValue $lastProgress 'residual')
        leaf_panel_interactions = Get-RowValue $timing 'tabi_leaf_panel_interactions' (Get-RowValue $lastProgress 'leaf_panel_interactions')
    }
    $triage = New-TabiTriageRecord -Row $fallback -RunDir $RunDir -ExistingClassification 'NA'
    [PSCustomObject]@{
        case = $CaseName
        sum_mode = $SumMode
        formulation = $Formulation
        preconditioner_requested = $Preconditioner
        preconditioner_block_near_requested = $BlockNear
        status = $Status
        total_seconds = Get-RowValue $timing 'total_seconds'
        panels = Get-RowValue $timing 'tabi_panels' (Get-RowValue $lastProgress 'panels')
        matvec_calls = Get-RowValue $timing 'tabi_matvec_calls' (Get-RowValue $lastProgress 'matvec_calls')
        gmres_iterations = Get-RowValue $timing 'tabi_gmres_iterations' (Get-RowValue $timing 'iterations')
        final_residual = Get-RowValue $timing 'tabi_final_residual' (Get-RowValue $lastProgress 'residual')
        best_residual = Get-RowValue $timing 'tabi_best_residual'
        final_preconditioned_residual = Get-RowValue $timing 'tabi_final_preconditioned_residual' (Get-RowValue $lastProgress 'preconditioned_residual')
        best_preconditioned_residual = Get-RowValue $timing 'tabi_best_preconditioned_residual'
        preconditioner_effective = $triage.preconditioner_effective
        preconditioner_setup_seconds = $triage.preconditioner_setup_seconds
        preconditioner_apply_seconds = $triage.preconditioner_apply_seconds
        preconditioner_apply_calls = $triage.preconditioner_apply_calls
        preconditioner_valid_blocks = $triage.preconditioner_valid_blocks
        preconditioner_valid_panels = $triage.preconditioner_valid_panels
        preconditioner_ill_conditioned_blocks = $triage.preconditioner_ill_conditioned_blocks
        preconditioner_quality_rejected_blocks = Get-RowValue $timing 'tabi_preconditioner_quality_rejected_blocks' (Get-RowValue $lastProgress 'preconditioner_quality_rejected_blocks')
        preconditioner_quality_probe_max_gain = Get-RowValue $timing 'tabi_preconditioner_quality_probe_max_gain' (Get-RowValue $lastProgress 'preconditioner_quality_probe_max_gain')
        preconditioner_quality_status = Get-RowValue $timing 'tabi_preconditioner_quality_status' (Get-RowValue $lastProgress 'preconditioner_quality_status')
        seconds_per_matvec = $triage.seconds_per_matvec
        leaf_vs_n2_ratio = $triage.leaf_vs_n2_ratio
        triage_classification = $triage.triage_classification
        energy_status = Get-RowValue $energy 'status'
        energy_kJmol = Get-RowValue $energy 'E_rxn_kJmol'
        has_timing = if ($null -ne $timing) { 'true' } else { 'false' }
        run_dir = $RunDir
    }
}

if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
New-CleanDir $WorkRoot

$pqrPath = Join-Path $WorkRoot 'controlled_synthetic.pqr'
Write-SyntheticProteinPqr -Path $pqrPath

$cases = @(
    [PSCustomObject]@{ name = 'direct_diagonal'; sum = 'direct'; formulation = 'jcp2013'; preconditioner = 'diagonal'; blockNear = 'buffer-aware' },
    [PSCustomObject]@{ name = 'direct_block_buffer_aware'; sum = 'direct'; formulation = 'jcp2013'; preconditioner = 'block-budgeted'; blockNear = 'buffer-aware' },
    [PSCustomObject]@{ name = 'direct_block_operator_near'; sum = 'direct'; formulation = 'jcp2013'; preconditioner = 'block-budgeted'; blockNear = 'operator-near' },
    [PSCustomObject]@{ name = 'treecode_diagonal'; sum = 'treecode'; formulation = 'jcp2013-treecode'; preconditioner = 'diagonal'; blockNear = 'buffer-aware' },
    [PSCustomObject]@{ name = 'treecode_block_buffer_aware'; sum = 'treecode'; formulation = 'jcp2013-treecode'; preconditioner = 'block-budgeted'; blockNear = 'buffer-aware' },
    [PSCustomObject]@{ name = 'treecode_block_operator_near'; sum = 'treecode'; formulation = 'jcp2013-treecode'; preconditioner = 'block-budgeted'; blockNear = 'operator-near' }
)

$rows = @()
foreach ($case in $cases) {
    $runDir = Join-Path $WorkRoot $case.name
    New-CleanDir $runDir
    $argPath = New-ArgFile -RunDir $runDir -CaseName $case.name -SumMode $case.sum `
        -Formulation $case.formulation -Preconditioner $case.preconditioner `
        -BlockNear $case.blockNear -PqrPath $pqrPath
    $status = Invoke-TitanCase -RunDir $runDir -ArgPath $argPath
    $rows += Read-CaseMetrics -RunDir $runDir -CaseName $case.name -SumMode $case.sum `
        -Formulation $case.formulation -Preconditioner $case.preconditioner `
        -BlockNear $case.blockNear -Status $status
}

$outCsv = Join-Path $WorkRoot 'pbe_tabi_preconditioner_controlled_probe.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII

$blockRows = @($rows | Where-Object { $_.preconditioner_requested -eq 'block-budgeted' })
$activeBlockRows = @($blockRows | Where-Object {
    ([string]$_.preconditioner_effective).StartsWith('block_budgeted')
})
$summaryPath = Join-Path $WorkRoot 'pbe_tabi_preconditioner_controlled_probe_summary.dat'
Set-Content -LiteralPath $summaryPath -Encoding ASCII -Value @(
    'TABI controlled preconditioner probe',
    "surface_points: $SurfacePoints",
    "threads: $Threads",
    "case_timeout_seconds: $CaseTimeoutSeconds",
    "cases: $($cases.Count)",
    "block_budgeted_cases: $($blockRows.Count)",
    "block_budgeted_active_cases: $($activeBlockRows.Count)",
    "report_csv: $outCsv",
    'pass_policy: direct operator-near block-budgeted must activate and write pbe_timing.csv; treecode activation is diagnostic.'
)

$directBlock = @($rows | Where-Object { $_.case -eq 'direct_block_operator_near' })[0]
$directDiagonal = @($rows | Where-Object { $_.case -eq 'direct_diagonal' })[0]
if ($directBlock.has_timing -ne 'true') {
    throw "direct_block_operator_near did not write pbe_timing.csv; see $($directBlock.run_dir)"
}
if (-not ([string]$directBlock.preconditioner_effective).StartsWith('block_budgeted')) {
    throw "direct_block_operator_near did not activate block-budgeted preconditioner: $($directBlock.preconditioner_effective)"
}
$directBlockMatvec = Convert-TabiTriageDoubleOrNull $directBlock.matvec_calls
$directDiagonalMatvec = Convert-TabiTriageDoubleOrNull $directDiagonal.matvec_calls
if ($null -eq $directBlockMatvec -or $null -eq $directDiagonalMatvec -or
    $directBlockMatvec -gt $directDiagonalMatvec) {
    throw "direct_block_operator_near did not improve matvecs vs diagonal: block=$($directBlock.matvec_calls) diagonal=$($directDiagonal.matvec_calls)"
}

Write-Host "PBE TABI controlled preconditioner probe PASS: $outCsv"
