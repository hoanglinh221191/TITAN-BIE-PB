param(
    [switch]$Quick,
    [string[]]$Fixtures = @('research-medium','research-large'),
    [int]$SurfacePoints = 4,
    [int]$MaxPanels = 1024,
    [int]$LeafSize = 64,
    [double]$Theta = 0.8,
    [int]$Threads = 4,
    [int]$Restart = 24,
    [int]$MaxIterations = 72,
    [double]$Tolerance = 1.0e-6,
    [int]$Repeats = 3,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_k'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutCsv = Join-Path $WorkRoot 'fmm_m2l_scope_k_matrix.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_k_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_K_REPORT.txt'

if ($Quick) {
    $Fixtures = @('research-medium')
    $MaxPanels = [Math]::Max(1024, $MaxPanels)
    $MaxIterations = [Math]::Min($MaxIterations, 48)
}

New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if ([string]::IsNullOrWhiteSpace($Cxx)) {
    $Clang = 'C:\msys64\ucrt64\bin\clang++.exe'
    $Cxx = if (Test-Path -LiteralPath $Clang) { $Clang } else { 'g++' }
}

$ExePath = Join-Path $WorkRoot 'titan_fmm_m2l_scope_k_probe.exe'
$Sources = @(
    'test\titan_fmm_m2l_scope_k_probe.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_fallback.cpp',
    'src\tabi\tabi_fmm_m2l_cartesian.cpp',
    'src\tabi\tabi_fmm_m2l_local.cpp',
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
    'src\tabi\tabi_treecode_operator.cpp'
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
            throw "Unknown Scope K fixture: $Fixture"
        }
    }
}

Push-Location $RepoRoot
try {
    & $Cxx @BuildArgs
    if ($LASTEXITCODE -ne 0) {
        throw 'failed to build FMM/M2L Scope K probe executable'
    }

    & $ExePath
    if ($LASTEXITCODE -ne 0) {
        throw 'generated Scope K GMRES smoke failed'
    }

    $Rows = New-Object System.Collections.Generic.List[string]
    $Header = 'fixture,atoms,surface_panels,sampled_panels,surface_seconds,leaf_size,theta,restart,max_iterations,tolerance,direct_iterations,adapter_iterations,direct_final_residual,adapter_final_residual,direct_true_residual,solution_rel_error,direct_matvec_ms,adapter_matvec_ms,speedup_vs_direct_matvec,near_pairs,far_pairs,adapter_setup_ms,adapter_apply_ms,finite_outputs,scope_k_gate'
    $Rows.Add($Header)
    $MediumPassed = $false

    foreach ($Fixture in $Fixtures) {
        if ($Fixture -eq 'research-large' -and -not $MediumPassed) {
            continue
        }
        $PqrPath = Resolve-Fixture $Fixture
        if (-not (Test-Path -LiteralPath $PqrPath)) {
            throw "Missing Scope K PQR fixture: $PqrPath"
        }
        $ProbeArgs = @(
            '--fixture-probe',
            "--fixture=$Fixture",
            "--pqr=$((Resolve-Path -LiteralPath $PqrPath).Path)",
            "--surface-points=$SurfacePoints",
            "--max-panels=$MaxPanels",
            "--leaf-size=$LeafSize",
            "--theta=$Theta",
            "--threads=$Threads",
            "--restart=$Restart",
            "--max-iterations=$MaxIterations",
            "--tolerance=$Tolerance"
        )
        $Output = @(& $ExePath @ProbeArgs)
        $ExitCode = $LASTEXITCODE
        $DataLine = @($Output | Where-Object { $_ -and $_ -notmatch '^fixture,' }) | Select-Object -Last 1
        if ([string]::IsNullOrWhiteSpace($DataLine)) {
            throw "Scope K probe emitted no CSV row for $Fixture"
        }
        $Rows.Add($DataLine)
        $Parsed = $DataLine | ConvertFrom-Csv -Header ($Header -split ',')
        if ($Fixture -eq 'research-medium') {
            $MediumPassed = $Parsed.scope_k_gate -eq 'pass'
        }
        if ($ExitCode -ne 0 -and $Parsed.scope_k_gate -eq 'pass') {
            throw "Scope K probe returned exit code $ExitCode for passing fixture $Fixture"
        }
    }

    Set-Content -LiteralPath $OutCsv -Value $Rows -Encoding ASCII
    $Data = @(Import-Csv -LiteralPath $OutCsv)
    $PassRows = @($Data | Where-Object { $_.scope_k_gate -eq 'pass' })

    $Summary = New-Object System.Collections.Generic.List[string]
    $Summary.Add('TITAN FMM/M2L Scope K Probe Summary')
    $Summary.Add('==================================')
    $Summary.Add('')
    $Summary.Add("quick_mode: $($Quick.IsPresent)")
    $Summary.Add('scope: guarded_gmres_experiment')
    $Summary.Add('hidden_option: -pbe-tabi-fmm-m2l=single-matvec-analytic')
    $Summary.Add('production_default: unchanged_grid_default')
    $Summary.Add('rhs_policy: b = A_direct * deterministic_x_true')
    $Summary.Add("surface_points: $SurfacePoints")
    $Summary.Add("max_sampled_panels: $MaxPanels")
    $Summary.Add("leaf_size: $LeafSize")
    $Summary.Add("theta: $Theta")
    $Summary.Add("restart: $Restart")
    $Summary.Add("max_iterations: $MaxIterations")
    $Summary.Add("tolerance: $Tolerance")
    $Summary.Add("rows: $($Data.Count)")
    $Summary.Add("scope_k_gate_pass_rows: $($PassRows.Count)")
    foreach ($Row in $Data) {
        $Summary.Add("row: fixture=$($Row.fixture) sampled_panels=$($Row.sampled_panels) adapter_residual=$($Row.adapter_final_residual) direct_true_residual=$($Row.direct_true_residual) solution_rel_error=$($Row.solution_rel_error) direct_matvec_ms=$($Row.direct_matvec_ms) adapter_matvec_ms=$($Row.adapter_matvec_ms) speedup=$($Row.speedup_vs_direct_matvec) near_pairs=$($Row.near_pairs) far_pairs=$($Row.far_pairs) gate=$($Row.scope_k_gate)")
    }
    $Summary.Add('')
    $Summary.Add('acceptance: finite outputs, adapter residual <= 1e-6, direct true residual <= 5e-3 on medium or <= 2e-2 on large, solution relative error <= 5e-2, speedup >= 1.2 on medium.')
    $Summary.Add('next_gate: If Scope K passes medium, plan Scope L hidden end-to-end PB solve and energy comparison. If it fails, do not wire production GMRES.')
    $Summary.Add("matrix_csv: $OutCsv")
    Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
    Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

    Write-Host "FMM/M2L Scope K probe complete: $OutCsv"
} finally {
    Pop-Location
}
