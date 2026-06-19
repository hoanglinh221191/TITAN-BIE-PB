param(
    [switch]$Quick,
    [string[]]$Fixtures = @('research-medium','research-large'),
    [int]$SurfacePoints = 4,
    [int]$MaxPanels = 2048,
    [int]$LeafSize = 64,
    [double]$Theta = 0.8,
    [int]$Threads = 4,
    [int]$Repeats = 5,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_j'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutCsv = Join-Path $WorkRoot 'fmm_m2l_scope_j_matrix.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_j_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_J_REPORT.txt'

if ($Quick) {
    $Fixtures = @('research-medium')
    $MaxPanels = [Math]::Min($MaxPanels, 1024)
    $Repeats = [Math]::Min($Repeats, 3)
}

New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if ([string]::IsNullOrWhiteSpace($Cxx)) {
    $Clang = 'C:\msys64\ucrt64\bin\clang++.exe'
    $Cxx = if (Test-Path -LiteralPath $Clang) { $Clang } else { 'g++' }
}

$ExePath = Join-Path $WorkRoot 'titan_fmm_m2l_scope_j_probe.exe'
$Sources = @(
    'test\titan_fmm_m2l_scope_j_probe.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
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
            throw "Unknown Scope J fixture: $Fixture"
        }
    }
}

Push-Location $RepoRoot
try {
    & $Cxx @BuildArgs
    if ($LASTEXITCODE -ne 0) {
        throw 'failed to build FMM/M2L Scope J probe executable'
    }

    & $ExePath
    if ($LASTEXITCODE -ne 0) {
        throw 'generated production-panel Scope J smoke failed'
    }

    $Rows = New-Object System.Collections.Generic.List[string]
    $Header = 'fixture,atoms,surface_panels,sampled_panels,surface_seconds,production_direct_seconds,leaf_size,theta,repeats,adapter_phi_rel_rms,adapter_phin_rel_rms,adapter_aggregate_rel_rms,adapter_max_abs,near_pairs,far_pairs,setup_ms,apply_ms,adapter_total_seconds,speedup_vs_production_direct,scope_j_gate'
    $Rows.Add($Header)
    $MediumPassed = $false

    foreach ($Fixture in $Fixtures) {
        if ($Fixture -eq 'research-large' -and -not $MediumPassed) {
            continue
        }
        $PqrPath = Resolve-Fixture $Fixture
        if (-not (Test-Path -LiteralPath $PqrPath)) {
            throw "Missing Scope J PQR fixture: $PqrPath"
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
            "--repeats=$Repeats"
        )
        $Output = @(& $ExePath @ProbeArgs)
        $ExitCode = $LASTEXITCODE
        $DataLine = @($Output | Where-Object { $_ -and $_ -notmatch '^fixture,' }) | Select-Object -Last 1
        if ([string]::IsNullOrWhiteSpace($DataLine)) {
            throw "Scope J probe emitted no CSV row for $Fixture"
        }
        $Rows.Add($DataLine)
        $Parsed = $DataLine | ConvertFrom-Csv -Header ($Header -split ',')
        if ($Fixture -eq 'research-medium') {
            $MediumPassed = $Parsed.scope_j_gate -eq 'pass'
        }
        if ($ExitCode -ne 0 -and $Parsed.scope_j_gate -eq 'pass') {
            throw "Scope J probe returned exit code $ExitCode for passing fixture $Fixture"
        }
    }

    Set-Content -LiteralPath $OutCsv -Value $Rows -Encoding ASCII
    $Data = @(Import-Csv -LiteralPath $OutCsv)
    $PassRows = @($Data | Where-Object { $_.scope_j_gate -eq 'pass' })

    $Summary = New-Object System.Collections.Generic.List[string]
    $Summary.Add('TITAN FMM/M2L Scope J Probe Summary')
    $Summary.Add('==================================')
    $Summary.Add('')
    $Summary.Add("quick_mode: $($Quick.IsPresent)")
    $Summary.Add('scope: hidden_production_panel_single_matvec_adapter')
    $Summary.Add('hidden_option: -pbe-tabi-fmm-m2l=single-matvec-analytic')
    $Summary.Add('production_default: unchanged_grid_default')
    $Summary.Add('sign_convention: output equals production apply_jcp_treecode_operator convention lhs*unknown - integral_sum')
    $Summary.Add("surface_points: $SurfacePoints")
    $Summary.Add("max_sampled_panels: $MaxPanels")
    $Summary.Add("leaf_size: $LeafSize")
    $Summary.Add("theta: $Theta")
    $Summary.Add("timing_repeats: $Repeats")
    $Summary.Add("rows: $($Data.Count)")
    $Summary.Add("scope_j_gate_pass_rows: $($PassRows.Count)")
    foreach ($Row in $Data) {
        $Summary.Add("row: fixture=$($Row.fixture) atoms=$($Row.atoms) surface_panels=$($Row.surface_panels) sampled_panels=$($Row.sampled_panels) aggregate=$($Row.adapter_aggregate_rel_rms) phi=$($Row.adapter_phi_rel_rms) phin=$($Row.adapter_phin_rel_rms) near_pairs=$($Row.near_pairs) far_pairs=$($Row.far_pairs) setup_ms=$($Row.setup_ms) apply_ms=$($Row.apply_ms) speedup=$($Row.speedup_vs_production_direct) gate=$($Row.scope_j_gate)")
    }
    $Summary.Add('')
    $Summary.Add('acceptance: finite outputs, aggregate relative RMS <= 1e-2 on real fixtures, speedup >= 1.2 on research-medium; research-large is gated behind medium pass.')
    $Summary.Add('next_gate: If Scope J passes, plan a guarded GMRES experiment. If it fails, optimize adapter internals or pause FMM/M2L integration.')
    $Summary.Add("matrix_csv: $OutCsv")
    Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
    Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

    Write-Host "FMM/M2L Scope J probe complete: $OutCsv"
} finally {
    Pop-Location
}
