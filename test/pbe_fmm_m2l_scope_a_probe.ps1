param(
    [string]$TitanExe = ".\titan.exe",
    [switch]$Quick,
    [int[]]$PanelCounts = @(2048, 8192, 32768),
    [int[]]$Orders = @(2, 4, 6),
    [double[]]$Theta = @(0.45, 0.6, 0.8),
    [int]$Threads = 8,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_a'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutCsv = Join-Path $WorkRoot 'fmm_m2l_scope_a_matrix.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_a_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_A_REPORT.txt'

if ($Quick) {
    $PanelCounts = @(512, 1024)
    $Orders = @(2, 4)
    $Theta = @(0.35, 0.6)
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
    'src\titan_executor.cpp',
    '-o',
    $ExePath
)

Push-Location $RepoRoot
try {
    & $Cxx @BuildArgs
    if ($LASTEXITCODE -ne 0) {
        throw "failed to build FMM/M2L Scope A probe executable"
    }

    $Rows = New-Object System.Collections.Generic.List[string]
    $Header = 'panel_count,kernel,order,theta,status,setup_seconds,apply_seconds,direct_seconds,near_pairs,far_pairs,m2l_operators,rms_abs_error,max_abs_error,rel_rms_error,speedup_vs_direct,decision_class'
    $Rows.Add($Header)

    foreach ($Kernel in @('laplace', 'yukawa')) {
        foreach ($PanelCount in $PanelCounts) {
            foreach ($Order in $Orders) {
                foreach ($ThetaValue in $Theta) {
                    $ProbeArgs = @(
                        '--probe',
                        "--panels=$PanelCount",
                        "--kernel=$Kernel",
                        '--kappa=0.2',
                        '--leaf=64',
                        "--order=$Order",
                        "-theta=$ThetaValue".Replace('-theta', '--theta'),
                        "--threads=$Threads",
                        '--seed=29'
                    )
                    $Output = & $ExePath @ProbeArgs
                    if ($LASTEXITCODE -ne 0) {
                        throw "FMM/M2L probe failed for kernel=$Kernel panels=$PanelCount order=$Order theta=$ThetaValue"
                    }
                    foreach ($Line in $Output) {
                        if ($Line -and $Line -notmatch '^panel_count,') {
                            $Rows.Add($Line)
                        }
                    }
                }
            }
        }
    }

    Set-Content -LiteralPath $OutCsv -Value $Rows -Encoding ASCII
    $Data = Import-Csv -LiteralPath $OutCsv
    $Best = $Data |
        Where-Object { $_.status -eq 'ok' } |
        Sort-Object @{Expression = {[double]$_.rel_rms_error}; Ascending = $true},
                    @{Expression = {[double]$_.apply_seconds}; Ascending = $true} |
        Select-Object -First 1
    $Promoted = @($Data | Where-Object { $_.decision_class -eq 'promote_to_scope_b' })
    $Promising = @($Data | Where-Object { $_.decision_class -eq 'promising_scalar' })

    $Summary = New-Object System.Collections.Generic.List[string]
    $Summary.Add('TITAN FMM/M2L Scope A Probe Summary')
    $Summary.Add('==================================')
    $Summary.Add("")
    $Summary.Add("quick_mode: $($Quick.IsPresent)")
    $Summary.Add("threads: $Threads")
    $Summary.Add("rows: $($Data.Count)")
    $Summary.Add("promote_to_scope_b_rows: $($Promoted.Count)")
    $Summary.Add("promising_scalar_rows: $($Promising.Count)")
    if ($Best) {
        $Summary.Add("best_row: panels=$($Best.panel_count) kernel=$($Best.kernel) order=$($Best.order) theta=$($Best.theta) rel_rms_error=$($Best.rel_rms_error) speedup=$($Best.speedup_vs_direct) decision=$($Best.decision_class)")
    } else {
        $Summary.Add('best_row: none')
    }
    $Summary.Add("")
    $Summary.Add('decision_note: Scope A is research-only. Do not wire M2L into TABI JCP unless a row reaches promote_to_scope_b or the user explicitly overrides the gate.')
    $Summary.Add("matrix_csv: $OutCsv")
    Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
    Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

    Write-Host "FMM/M2L Scope A probe complete: $OutCsv"
} finally {
    Pop-Location
}
