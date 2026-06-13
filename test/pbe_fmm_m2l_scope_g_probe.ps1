param(
    [switch]$Quick,
    [string[]]$Fixtures = @('research-medium','research-large'),
    [int]$SurfacePoints = 4,
    [int]$MaxPanels = 2048,
    [int]$Order = 2,
    [double]$Theta = 0.8,
    [int]$Threads = 4,
    [int]$Repeats = 5,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_g'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutCsv = Join-Path $WorkRoot 'fmm_m2l_scope_g_matrix.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_g_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_G_REPORT.txt'

if ($Quick) {
    $Fixtures = @('research-medium')
    $MaxPanels = [Math]::Min($MaxPanels, 1024)
}

New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if ([string]::IsNullOrWhiteSpace($Cxx)) {
    $Clang = 'C:\msys64\ucrt64\bin\clang++.exe'
    $Cxx = if (Test-Path -LiteralPath $Clang) { $Clang } else { 'g++' }
}

$ExePath = Join-Path $WorkRoot 'titan_fmm_m2l_scope_g_probe.exe'
$Sources = @(
    'test\titan_fmm_m2l_scope_g_probe.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
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
            throw "Unknown Scope G fixture: $Fixture"
        }
    }
}

Push-Location $RepoRoot
try {
    & $Cxx @BuildArgs
    if ($LASTEXITCODE -ne 0) {
        throw 'failed to build FMM/M2L Scope G probe executable'
    }

    & $ExePath
    if ($LASTEXITCODE -ne 0) {
        throw 'generated production-surface Scope G smoke failed'
    }

    $Rows = New-Object System.Collections.Generic.List[string]
    $Header = 'fixture,atoms,surface_panels,sampled_panels,surface_seconds,production_direct_seconds,order,theta,repeats,parity_phi_rel_rms,parity_phin_rel_rms,parity_aggregate_rel_rms,k1_rel_rms,derivative_aggregate_rel_rms,derivative_max_channel_rel_rms,derivative_significant_max_rel_rms,k2_rel_rms,k3_rel_rms,k4_rel_rms,k2_reference_rms,k3_reference_rms,k4_reference_rms,near_pairs,far_pairs,setup_seconds,apply_seconds,prototype_direct_seconds,speedup_vs_direct,scope_g_gate'
    $Rows.Add($Header)
    $MediumPassed = $false

    foreach ($Fixture in $Fixtures) {
        if ($Fixture -eq 'research-large' -and -not $MediumPassed) {
            continue
        }
        $PqrPath = Resolve-Fixture $Fixture
        if (-not (Test-Path -LiteralPath $PqrPath)) {
            throw "Missing Scope G PQR fixture: $PqrPath"
        }
        $ProbeArgs = @(
            '--fixture-probe',
            "--fixture=$Fixture",
            "--pqr=$((Resolve-Path -LiteralPath $PqrPath).Path)",
            "--surface-points=$SurfacePoints",
            "--max-panels=$MaxPanels",
            "--order=$Order",
            "--theta=$Theta",
            "--threads=$Threads",
            "--repeats=$Repeats"
        )
        $Output = @(& $ExePath @ProbeArgs)
        $ExitCode = $LASTEXITCODE
        $DataLine = @($Output | Where-Object { $_ -and $_ -notmatch '^fixture,' }) | Select-Object -Last 1
        if ([string]::IsNullOrWhiteSpace($DataLine)) {
            throw "Scope G probe emitted no CSV row for $Fixture"
        }
        $Rows.Add($DataLine)
        $Parsed = $DataLine | ConvertFrom-Csv -Header ($Header -split ',')
        if ($Fixture -eq 'research-medium') {
            $MediumPassed = $Parsed.scope_g_gate -eq 'pass'
        }
        if ($ExitCode -ne 0 -and $Parsed.scope_g_gate -eq 'pass') {
            throw "Scope G probe returned exit code $ExitCode for passing fixture $Fixture"
        }
    }

    Set-Content -LiteralPath $OutCsv -Value $Rows -Encoding ASCII
    $Data = @(Import-Csv -LiteralPath $OutCsv)
    $PassRows = @($Data | Where-Object { $_.scope_g_gate -eq 'pass' })

    $Summary = New-Object System.Collections.Generic.List[string]
    $Summary.Add('TITAN FMM/M2L Scope G Probe Summary')
    $Summary.Add('==================================')
    $Summary.Add('')
    $Summary.Add("quick_mode: $($Quick.IsPresent)")
    $Summary.Add('surface_source: TitanPBE::Tabi::build_surface')
    $Summary.Add("surface_points: $SurfacePoints")
    $Summary.Add("max_sampled_panels: $MaxPanels")
    $Summary.Add("order: $Order")
    $Summary.Add("theta: $Theta")
    $Summary.Add("timing_repeats: $Repeats")
    $Summary.Add("rows: $($Data.Count)")
    $Summary.Add("scope_g_gate_pass_rows: $($PassRows.Count)")
    foreach ($Row in $Data) {
        $Summary.Add("row: fixture=$($Row.fixture) atoms=$($Row.atoms) surface_panels=$($Row.surface_panels) sampled_panels=$($Row.sampled_panels) parity=$($Row.parity_aggregate_rel_rms) k1=$($Row.k1_rel_rms) derivative_aggregate=$($Row.derivative_aggregate_rel_rms) derivative_significant_max=$($Row.derivative_significant_max_rel_rms) near_pairs=$($Row.near_pairs) far_pairs=$($Row.far_pairs) speedup=$($Row.speedup_vs_direct) gate=$($Row.scope_g_gate)")
    }
    $Summary.Add('')
    $Summary.Add('channel_policy: max-relative gate applies only to channels whose reference RMS is at least 0.1% of the strongest derivative channel; aggregate relative RMS always remains mandatory.')
    $Summary.Add('decision_note: Scope G is an offline production-surface comparison. It does not run GMRES or change the production PBE backend.')
    $Summary.Add("matrix_csv: $OutCsv")
    Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
    Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

    Write-Host "FMM/M2L Scope G probe complete: $OutCsv"
} finally {
    Pop-Location
}
