param(
    [switch]$Quick,
    [string[]]$Fixtures = @('research-medium','research-large'),
    [int]$AtomCap = 128,
    [int]$SurfacePoints = 8,
    [int]$MaxIter = 120,
    [int]$Restart = 24,
    [double]$Tolerance = 1.0e-6,
    [int]$Threads = 4,
    [int]$Repeats = 3,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_l'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutCsv = Join-Path $WorkRoot 'fmm_m2l_scope_l_matrix.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_l_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_L_REPORT.txt'

if ($Quick) {
    $Fixtures = @('research-medium')
    $AtomCap = [Math]::Min($AtomCap, 64)
    $MaxIter = [Math]::Min($MaxIter, 80)
    $Repeats = [Math]::Min($Repeats, 1)
}

New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if ([string]::IsNullOrWhiteSpace($Cxx)) {
    $Clang = 'C:\msys64\ucrt64\bin\clang++.exe'
    $Cxx = if (Test-Path -LiteralPath $Clang) { $Clang } else { 'g++' }
}

$ExePath = Join-Path $WorkRoot 'titan_fmm_m2l_scope_l_probe.exe'
$Sources = @(
    'test\titan_fmm_m2l_scope_l_probe.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_preconditioner.cpp',
    'src\tabi\tabi_solver.cpp',
    'src\tabi\tabi_surface.cpp',
    'src\titan_surface.cpp',
    'src\titan_executor.cpp',
    'src\tabi\tabi_near.cpp',
    'src\tabi\tabi_treecode.cpp',
    'src\tabi\tabi_treecode_internal.cpp',
    'src\tabi\tabi_treecode_geometry.cpp',
    'src\tabi\tabi_treecode_series.cpp',
    'src\tabi\tabi_treecode_policy.cpp',
    'src\tabi\tabi_treecode_leaf.cpp',
    'src\tabi\tabi_treecode_operator.cpp',
    'src\tabi\tabi_preconditioner.cpp'
)
$BuildArgs = @(
    '-std=c++23', '-O2', '-mavx2', '-mfma',
    '-DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1', '-Isrc'
) + $Sources + @('-o', $ExePath)

function Resolve-Fixture([string]$Fixture) {
    switch ($Fixture) {
        'research-medium' { return Join-Path $RepoRoot 'example\1\rec_027.pqr' }
        'research-large' { return Join-Path $RepoRoot 'example\4\rec_0.pqr' }
        default {
            if (Test-Path -LiteralPath $Fixture) {
                return [System.IO.Path]::GetFullPath($Fixture)
            }
            throw "Unknown Scope L fixture: $Fixture"
        }
    }
}

Push-Location $RepoRoot
try {
    & $Cxx @BuildArgs
    if ($LASTEXITCODE -ne 0) {
        throw 'failed to build FMM/M2L Scope L probe executable'
    }

    & $ExePath
    if ($LASTEXITCODE -ne 0) {
        throw 'generated Scope L solve smoke failed'
    }

    $Rows = New-Object System.Collections.Generic.List[string]
    $Header = 'fixture,input_atoms,sampled_atoms,direct_status,fmm_status,direct_panels,fmm_panels,direct_iterations,fmm_iterations,direct_final_residual,fmm_final_residual,direct_true_residual,fmm_true_residual,direct_energy_kj_mol,fmm_energy_kj_mol,energy_delta_kj_mol,energy_delta_abs_kj_mol,energy_delta_rel,direct_matvec_ms_avg,fmm_matvec_ms_avg,matvec_speedup,fmm_m2l_mode,fmm_near_pairs,fmm_far_pairs,fmm_setup_ms,fmm_apply_ms,direct_seconds,fmm_seconds,scope_l_gate'
    $Rows.Add($Header)
    $MediumPassed = $false

    foreach ($Fixture in $Fixtures) {
        if ($Fixture -eq 'research-large' -and -not $MediumPassed) {
            continue
        }
        $PqrPath = Resolve-Fixture $Fixture
        if (-not (Test-Path -LiteralPath $PqrPath)) {
            throw "Missing Scope L PQR fixture: $PqrPath"
        }
        $ProbeArgs = @(
            '--fixture-probe',
            "--fixture=$Fixture",
            "--pqr=$((Resolve-Path -LiteralPath $PqrPath).Path)",
            "--atom-cap=$AtomCap",
            "--surface-points=$SurfacePoints",
            "--max-iter=$MaxIter",
            "--restart=$Restart",
            "--tolerance=$Tolerance",
            "--threads=$Threads"
        )
        $Output = @(& $ExePath @ProbeArgs)
        $ExitCode = $LASTEXITCODE
        $DataLine = @($Output | Where-Object { $_ -and $_ -notmatch '^fixture,' }) | Select-Object -Last 1
        if ([string]::IsNullOrWhiteSpace($DataLine)) {
            throw "Scope L probe emitted no CSV row for $Fixture"
        }
        $Rows.Add($DataLine)
        $Parsed = $DataLine | ConvertFrom-Csv -Header ($Header -split ',')
        if ($Fixture -eq 'research-medium') {
            $MediumPassed = $Parsed.scope_l_gate -eq 'pass'
        }
        if ($ExitCode -ne 0 -and $Parsed.scope_l_gate -eq 'pass') {
            throw "Scope L probe returned exit code $ExitCode for passing fixture $Fixture"
        }
    }

    Set-Content -LiteralPath $OutCsv -Value $Rows -Encoding ASCII
    $Data = @(Import-Csv -LiteralPath $OutCsv)
    $PassRows = @($Data | Where-Object { $_.scope_l_gate -eq 'pass' })

    $Summary = New-Object System.Collections.Generic.List[string]
    $Summary.Add('TITAN FMM/M2L Scope L Probe Summary')
    $Summary.Add('==================================')
    $Summary.Add('')
    $Summary.Add("quick_mode: $($Quick.IsPresent)")
    $Summary.Add('scope: hidden_end_to_end_pb_solve_energy_gate')
    $Summary.Add('hidden_option: -pbe-tabi-fmm-m2l=single-matvec-analytic')
    $Summary.Add('production_default: unchanged_grid_default')
    $Summary.Add("atom_cap: $AtomCap")
    $Summary.Add("surface_points: $SurfacePoints")
    $Summary.Add("max_iter: $MaxIter")
    $Summary.Add("restart: $Restart")
    $Summary.Add("tolerance: $Tolerance")
    $Summary.Add("rows: $($Data.Count)")
    $Summary.Add("scope_l_gate_pass_rows: $($PassRows.Count)")
    foreach ($Row in $Data) {
        $Summary.Add("row: fixture=$($Row.fixture) atoms=$($Row.sampled_atoms) direct_status=$($Row.direct_status) fmm_status=$($Row.fmm_status) direct_energy=$($Row.direct_energy_kj_mol) fmm_energy=$($Row.fmm_energy_kj_mol) energy_delta_abs=$($Row.energy_delta_abs_kj_mol) direct_residual=$($Row.direct_true_residual) fmm_residual=$($Row.fmm_true_residual) speedup=$($Row.matvec_speedup) fmm_mode=$($Row.fmm_m2l_mode) fmm_far_pairs=$($Row.fmm_far_pairs) gate=$($Row.scope_l_gate)")
    }
    $Summary.Add('')
    $Summary.Add('acceptance: direct/FMM status ok, finite energy and residual, FMM true residual bounded against direct, energy delta within fixture gate, and matvec speedup >= 1.2 on medium.')
    $Summary.Add('next_gate: If Scope L passes, plan Scope M controlled titan.exe hidden CLI smoke and broader energy/regression matrix.')
    $Summary.Add("matrix_csv: $OutCsv")
    Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
    Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

    Write-Host "FMM/M2L Scope L probe complete: $OutCsv"
} finally {
    Pop-Location
}
