param(
    [switch]$Quick,
    [int[]]$PanelCounts = @(1024, 2048, 4096),
    [int[]]$Orders = @(2, 4),
    [double[]]$Theta = @(0.45, 0.6),
    [int]$Threads = 8,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_b'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutCsv = Join-Path $WorkRoot 'fmm_m2l_scope_b_matrix.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_b_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_B_REPORT.txt'

if ($Quick) {
    $PanelCounts = @(512, 1024)
    $Orders = @(2)
    $Theta = @(0.55)
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
        throw "failed to build FMM/M2L Scope B probe executable"
    }

    $Rows = New-Object System.Collections.Generic.List[string]
    $Header = 'fixture,panel_count,kernel,order,theta,status,setup_seconds,apply_seconds,direct_seconds,near_pairs,far_pairs,m2l_operators,rms_abs_error,max_abs_error,rel_rms_error,speedup_vs_direct,decision_class,scope_b_gate'
    $Rows.Add($Header)

    foreach ($Kernel in @('laplace', 'yukawa')) {
        foreach ($PanelCount in $PanelCounts) {
            foreach ($Order in $Orders) {
                foreach ($ThetaValue in $Theta) {
                    $ProbeArgs = @(
                        '--probe',
                        '--fixture=mini-protein-like',
                        "--panels=$PanelCount",
                        "--kernel=$Kernel",
                        '--kappa=0.2',
                        '--leaf=64',
                        "--order=$Order",
                        "--theta=$ThetaValue",
                        "--threads=$Threads",
                        '--seed=37'
                    )
                    $Output = & $ExePath @ProbeArgs
                    if ($LASTEXITCODE -ne 0) {
                        throw "FMM/M2L Scope B probe failed for kernel=$Kernel panels=$PanelCount order=$Order theta=$ThetaValue"
                    }
                    foreach ($Line in $Output) {
                        if ($Line -and $Line -notmatch '^panel_count,') {
                            $Fields = $Line -split ','
                            $Rel = [double]$Fields[13]
                            $Apply = [double]$Fields[6]
                            $Direct = [double]$Fields[7]
                            $Setup = [double]$Fields[5]
                            $Gate = if ($Rel -le 5.0e-2 -and $Apply -lt $Direct -and $Setup -le $Direct -and $Fields[4] -eq 'ok') {
                                'pass'
                            } else {
                                'fail'
                            }
                            $Rows.Add("mini-protein-like,$Line,$Gate")
                        }
                    }
                }
            }
        }
    }

    Set-Content -LiteralPath $OutCsv -Value $Rows -Encoding ASCII
    $Data = Import-Csv -LiteralPath $OutCsv
    $PassRows = @($Data | Where-Object { $_.scope_b_gate -eq 'pass' })
    $Best = $Data |
        Sort-Object @{Expression = {[double]$_.rel_rms_error}; Ascending = $true},
                    @{Expression = {[double]$_.apply_seconds}; Ascending = $true} |
        Select-Object -First 1

    $Summary = New-Object System.Collections.Generic.List[string]
    $Summary.Add('TITAN FMM/M2L Scope B Probe Summary')
    $Summary.Add('==================================')
    $Summary.Add("")
    $Summary.Add("quick_mode: $($Quick.IsPresent)")
    $Summary.Add("fixture: mini-protein-like")
    $Summary.Add("threads: $Threads")
    $Summary.Add("rows: $($Data.Count)")
    $Summary.Add("scope_b_gate_pass_rows: $($PassRows.Count)")
    if ($Best) {
        $Summary.Add("best_row: panels=$($Best.panel_count) kernel=$($Best.kernel) order=$($Best.order) theta=$($Best.theta) rel_rms_error=$($Best.rel_rms_error) speedup=$($Best.speedup_vs_direct) gate=$($Best.scope_b_gate)")
    } else {
        $Summary.Add('best_row: none')
    }
    $Summary.Add("")
    $Summary.Add('decision_note: Scope B uses an explicit point-cloud API and a deterministic mini-protein-like fixture. It is not wired into production TABI/JCP.')
    $Summary.Add("matrix_csv: $OutCsv")
    Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
    Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

    Write-Host "FMM/M2L Scope B probe complete: $OutCsv"
} finally {
    Pop-Location
}
