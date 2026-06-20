param(
    [switch]$Quick,
    [int[]]$PanelCounts = @(1024, 2048),
    [int[]]$Orders = @(2, 4),
    [double[]]$Theta = @(0.6, 0.7, 0.8),
    [int]$Threads = 8,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_f'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutCsv = Join-Path $WorkRoot 'fmm_m2l_scope_f_matrix.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_f_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_F_REPORT.txt'

if ($Quick) {
    $PanelCounts = @(512, 1024)
    $Orders = @(2)
    $Theta = @(0.6, 0.7)
}

New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if ([string]::IsNullOrWhiteSpace($Cxx)) {
    $Clang = 'C:\msys64\ucrt64\bin\clang++.exe'
    if (Test-Path $Clang) {
        $Cxx = $Clang
    } else {
        $Cxx = 'g++'
    }
}

$ExePath = Join-Path $WorkRoot 'titan_fmm_m2l_prototype_smoke.exe'
$BuildArgs = @(
    '-std=c++23',
    '-O2',
    '-Isrc',
    'test\titan_fmm_m2l_prototype_smoke.cpp',
    'src\tabi\tabi_fmm_m2l_prototype.cpp',
    'src\tabi\tabi_fmm_m2l_fallback.cpp',
    'src\titan_executor.cpp',
    '-o',
    $ExePath
)

Push-Location $RepoRoot
try {
    & $Cxx @BuildArgs
    if ($LASTEXITCODE -ne 0) {
        throw "failed to build FMM/M2L Scope F probe executable"
    }

    $Rows = New-Object System.Collections.Generic.List[string]
    $Header = 'panel_count,order,theta,status,setup_seconds,apply_seconds,direct_seconds,channels,aggregate_rel_rms_error,max_channel_rel_rms_error,k2_rel_rms_error,k3_rel_rms_error,k4_rel_rms_error,max_abs_error,near_pairs,far_pairs,speedup_vs_direct,scope_f_gate'
    $Rows.Add($Header)

    foreach ($PanelCount in $PanelCounts) {
        foreach ($Order in $Orders) {
            foreach ($ThetaValue in $Theta) {
                $ProbeArgs = @(
                    '--jcp-derivative-bridge-probe',
                    "--panels=$PanelCount",
                    '--kappa=0.2',
                    '--leaf=64',
                    "--order=$Order",
                    "--theta=$ThetaValue",
                    "--threads=$Threads",
                    '--seed=73'
                )
                $Output = & $ExePath @ProbeArgs
                if ($LASTEXITCODE -ne 0) {
                    throw "FMM/M2L Scope F probe failed for panels=$PanelCount order=$Order theta=$ThetaValue"
                }
                foreach ($Line in $Output) {
                    if ($Line -and $Line -notmatch '^panel_count,') {
                        $Rows.Add($Line)
                    }
                }
            }
        }
    }

    Set-Content -LiteralPath $OutCsv -Value $Rows -Encoding ASCII
    $Data = @(Import-Csv -LiteralPath $OutCsv)
    $PassRows = @($Data | Where-Object { $_.scope_f_gate -eq 'pass' })
    $Best = $Data |
        Sort-Object @{Expression = { if ($_.scope_f_gate -eq 'pass') { 0 } else { 1 } }; Ascending = $true},
                    @{Expression = {[double]$_.aggregate_rel_rms_error}; Ascending = $true},
                    @{Expression = {[double]$_.apply_seconds}; Ascending = $true} |
        Select-Object -First 1

    $Summary = New-Object System.Collections.Generic.List[string]
    $Summary.Add('TITAN FMM/M2L Scope F Probe Summary')
    $Summary.Add('==================================')
    $Summary.Add("")
    $Summary.Add("quick_mode: $($Quick.IsPresent)")
    $Summary.Add("fixture: tabi-panel-like-jcp-derivative-bridge")
    $Summary.Add("threads: $Threads")
    $Summary.Add("rows: $($Data.Count)")
    $Summary.Add("scope_f_gate_pass_rows: $($PassRows.Count)")
    if ($Best) {
        $Summary.Add("best_row: panels=$($Best.panel_count) order=$($Best.order) theta=$($Best.theta) aggregate_rel_rms_error=$($Best.aggregate_rel_rms_error) max_channel_rel_rms_error=$($Best.max_channel_rel_rms_error) k2=$($Best.k2_rel_rms_error) k3=$($Best.k3_rel_rms_error) k4=$($Best.k4_rel_rms_error) near_pairs=$($Best.near_pairs) far_pairs=$($Best.far_pairs) speedup=$($Best.speedup_vs_direct) gate=$($Best.scope_f_gate)")
    } else {
        $Summary.Add('best_row: none')
    }
    $Summary.Add("")
    $Summary.Add('decision_note: Scope F is a finite-difference derivative-channel bridge for K2/K3/K4. It is not production derivative M2L and does not run GMRES.')
    $Summary.Add("matrix_csv: $OutCsv")
    Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
    Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

    Write-Host "FMM/M2L Scope F probe complete: $OutCsv"
} finally {
    Pop-Location
}
