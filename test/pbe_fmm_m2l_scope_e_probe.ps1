param(
    [switch]$Quick,
    [int[]]$PanelCounts = @(1024, 2048, 4096),
    [int[]]$Orders = @(2, 4),
    [double[]]$Theta = @(0.6, 0.7),
    [int]$Threads = 8,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_e'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutCsv = Join-Path $WorkRoot 'fmm_m2l_scope_e_matrix.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_e_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_E_REPORT.txt'

if ($Quick) {
    $PanelCounts = @(768, 1024)
    $Orders = @(2)
    $Theta = @(0.7)
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
        throw "failed to build FMM/M2L Scope E probe executable"
    }

    $Rows = New-Object System.Collections.Generic.List[string]
    $Header = 'panel_count,order,theta,status,setup_seconds,apply_seconds,direct_seconds,rel_rms_error,rms_abs_error,max_abs_error,near_pairs,far_pairs,speedup_vs_direct,scope_e_gate'
    $Rows.Add($Header)

    foreach ($PanelCount in $PanelCounts) {
        foreach ($Order in $Orders) {
            foreach ($ThetaValue in $Theta) {
                $ProbeArgs = @(
                    '--jcp-k1-bridge-probe',
                    "--panels=$PanelCount",
                    '--kappa=0.2',
                    '--leaf=64',
                    "--order=$Order",
                    "--theta=$ThetaValue",
                    "--threads=$Threads",
                    '--seed=71'
                )
                $Output = & $ExePath @ProbeArgs
                if ($LASTEXITCODE -ne 0) {
                    throw "FMM/M2L Scope E probe failed for panels=$PanelCount order=$Order theta=$ThetaValue"
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
    $PassRows = @($Data | Where-Object { $_.scope_e_gate -eq 'pass' })
    $Best = $Data |
        Sort-Object @{Expression = {[double]$_.rel_rms_error}; Ascending = $true},
                    @{Expression = {[double]$_.apply_seconds}; Ascending = $true} |
        Select-Object -First 1

    $Summary = New-Object System.Collections.Generic.List[string]
    $Summary.Add('TITAN FMM/M2L Scope E Probe Summary')
    $Summary.Add('==================================')
    $Summary.Add("")
    $Summary.Add("quick_mode: $($Quick.IsPresent)")
    $Summary.Add("fixture: tabi-panel-like-jcp-k1-bridge")
    $Summary.Add("threads: $Threads")
    $Summary.Add("rows: $($Data.Count)")
    $Summary.Add("scope_e_gate_pass_rows: $($PassRows.Count)")
    if ($Best) {
        $Summary.Add("best_row: panels=$($Best.panel_count) order=$($Best.order) theta=$($Best.theta) rel_rms_error=$($Best.rel_rms_error) near_pairs=$($Best.near_pairs) far_pairs=$($Best.far_pairs) speedup=$($Best.speedup_vs_direct) gate=$($Best.scope_e_gate)")
    } else {
        $Summary.Add('best_row: none')
    }
    $Summary.Add("")
    $Summary.Add('decision_note: Scope E validates only the real JCP K1 scalar bridge. K2/K3/K4 derivative M2L remains future Scope F.')
    $Summary.Add("matrix_csv: $OutCsv")
    Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
    Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

    Write-Host "FMM/M2L Scope E probe complete: $OutCsv"
} finally {
    Pop-Location
}
