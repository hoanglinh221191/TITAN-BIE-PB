param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fagbi_scope_z4i_parser')
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
if (-not (Test-Path -LiteralPath $TitanExe)) { throw "Missing titan executable: $TitanExe" }
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
if (Test-Path -LiteralPath $WorkRoot) { Remove-Item -LiteralPath $WorkRoot -Recurse -Force }
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$pqr = Join-Path $WorkRoot 'born.pqr'
$driver = Join-Path $WorkRoot 'driver.pdb'
Set-Content -LiteralPath $pqr -Encoding ASCII -Value @(
    'ATOM      1  ION ION     1       0.000   0.000   0.000  1.0000  3.0000       X',
    'END'
)
Set-Content -LiteralPath $driver -Encoding ASCII -Value @(
    'ATOM      1  ION ION A   1       0.000   0.000   0.000  1.00  0.00           X',
    'END'
)
$argPath = Join-Path $WorkRoot 'arg.in'
$lines = @(
    "-f=$($driver.Replace('\','/'))",
    '-1pdb=true',
    '-ff=AMBER99SB',
    '-g=1&1',
    "-path-a=$((Join-Path $RepoRoot 'lib_ff').Replace('\','/'))/",
    "-path-n=$((Join-Path $RepoRoot 'lib_ff').Replace('\','/'))/",
    '-pbe=tabi',
    "-pbe-test-pqr=$($pqr.Replace('\','/'))",
    '-pbe-surface=vdw',
    '-pbe-max-iter=8',
    '-pbe-rms=1e-4',
    '-pbe-dphi=1e-4',
    '-pbe-int-diel=1',
    '-pbe-ext-diel=80',
    '-pbe-salt=0.15',
    '-pbe-temp=298.15',
    '-pbe-timing=detail',
    '-nt=1',
    '-pbe-tabi-sum=treecode',
    '-pbe-tabi-formulation=jcp2013-treecode',
    '-pbe-tabi-fmm-m2l=single-matvec-analytic',
    '-pbe-tabi-surface-points=8',
    '-pbe-tabi-tree-leaf-size=8',
    '-pbe-tabi-opening-theta=0.45',
    '-pbe-tabi-near-integration=analytic_singular',
    '-pbe-tabi-quad-buffer=off',
    '-pbe-tabi-precond=diagonal',
    '-pbe-tabi-gmres-restart=8',
    '-pbe-tabi-matvec-audit=sampled-detail',
    '-pbe-tabi-matvec-audit-sample-count=128'
)
$lines | Set-Content -LiteralPath $argPath -Encoding ASCII

Push-Location $WorkRoot
try {
    $stdout = Join-Path $WorkRoot 'stdout.log'
    $stderr = Join-Path $WorkRoot 'stderr.log'
    $process = Start-Process -FilePath $TitanExe -ArgumentList @('-a',$argPath) `
        -WorkingDirectory $WorkRoot -RedirectStandardOutput $stdout `
        -RedirectStandardError $stderr -Wait -PassThru -NoNewWindow
    if ($process.ExitCode -ne 0) {
        $combined = (Get-Content -Raw $stdout) + (Get-Content -Raw $stderr)
        throw "sampled-detail parser/integration row failed (exit=$($process.ExitCode)): $combined"
    }
    foreach ($name in @('pbe_tabi_operator_attribution_samples.csv',
            'pbe_tabi_operator_attribution_summary.csv',
            'pbe_tabi_operator_attribution_levels.csv')) {
        if (-not (Test-Path -LiteralPath $name)) { throw "Missing sampled-detail output: $name" }
    }
    $levelRows = @(Import-Csv -LiteralPath 'pbe_tabi_operator_attribution_levels.csv')
    if ($levelRows.Count -gt 0) {
        $expectedColumns = @('configured_expansion_order','implemented_far_moment_degree',
            'probe','row_component','channel','source_depth','mac_bin','samples',
            'error_rms','reference_rms','relative_rms','max_abs_error')
        $actualColumns = @($levelRows[0].PSObject.Properties.Name)
        if (($actualColumns -join ',') -ne ($expectedColumns -join ',')) {
            throw "Malformed attribution level CSV columns: $($actualColumns -join ',')"
        }
        if (@($levelRows | Where-Object { $_.mac_bin -match ',' }).Count -gt 0) {
            throw 'Attribution MAC bin labels must not contain unquoted CSV delimiters'
        }
    }
} finally {
    Pop-Location
}

$unsupported = Join-Path $WorkRoot 'unsupported_direct'
New-Item -ItemType Directory -Force -Path $unsupported | Out-Null
$unsupportedArg = Join-Path $unsupported 'arg.in'
$lines | ForEach-Object {
    $_ -replace '^-pbe-tabi-sum=treecode$', '-pbe-tabi-sum=direct' `
       -replace '^-pbe-tabi-fmm-m2l=single-matvec-analytic$', '-pbe-tabi-fmm-m2l=off'
} | Set-Content -LiteralPath $unsupportedArg -Encoding ASCII
$unsupportedOut = Join-Path $unsupported 'stdout.log'
$unsupportedErr = Join-Path $unsupported 'stderr.log'
$unsupportedProcess = Start-Process -FilePath $TitanExe `
    -ArgumentList @('-a',$unsupportedArg) -WorkingDirectory $unsupported `
    -RedirectStandardOutput $unsupportedOut -RedirectStandardError $unsupportedErr `
    -Wait -PassThru -NoNewWindow
if ($unsupportedProcess.ExitCode -eq 0) {
    throw 'sampled-detail unexpectedly accepted an unsupported direct TABI path'
}
$unsupportedLog = (Get-Content -Raw $unsupportedOut) + (Get-Content -Raw $unsupportedErr)
if ($unsupportedLog -notmatch 'requires analytical FMM/M2L TABI treecode') {
    throw "sampled-detail unsupported-path error was unclear: $unsupportedLog"
}

Write-Host 'PASS pbe_fagbi_scope_z4i_parser_smoke'
