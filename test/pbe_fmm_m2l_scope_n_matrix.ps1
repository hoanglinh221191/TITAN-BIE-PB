param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 360,
    [int]$MaxIter = 120,
    [int]$Restart = 24,
    [double]$Tolerance = 1.0e-6,
    [int]$Threads = 4,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_n')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutCsv = Join-Path $WorkRoot 'fmm_m2l_scope_n_matrix.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_n_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_n_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_N_REPORT.txt'
$TemperatureK = 298.15
$OpeningTheta = 0.6
$AuditSampleCount = if ($Quick) { 128 } else { 512 }
if ($Quick) { $MaxIter = [Math]::Min($MaxIter, 80) }

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path

$Conditions = @(
    [pscustomobject]@{ id='medium_cap64_baseline'; tier='scale'; fixture='research-medium'; atom_cap=64; surface_points=8; salt_M=0.15; int_diel=1.0; ext_diel=80.0; require_speed=$false },
    [pscustomobject]@{ id='medium_cap128_baseline'; tier='scale'; fixture='research-medium'; atom_cap=128; surface_points=8; salt_M=0.15; int_diel=1.0; ext_diel=80.0; require_speed=$true },
    [pscustomobject]@{ id='large_cap128_baseline'; tier='scale'; fixture='research-large'; atom_cap=128; surface_points=8; salt_M=0.15; int_diel=1.0; ext_diel=80.0; require_speed=$true },
    [pscustomobject]@{ id='alt_cap128_baseline'; tier='fixture'; fixture='research-alt'; atom_cap=128; surface_points=8; salt_M=0.15; int_diel=1.0; ext_diel=80.0; require_speed=$false },
    [pscustomobject]@{ id='medium_cap64_salt0'; tier='physics'; fixture='research-medium'; atom_cap=64; surface_points=8; salt_M=0.0; int_diel=1.0; ext_diel=80.0; require_speed=$false },
    [pscustomobject]@{ id='medium_cap64_salt05'; tier='physics'; fixture='research-medium'; atom_cap=64; surface_points=8; salt_M=0.5; int_diel=1.0; ext_diel=80.0; require_speed=$false },
    [pscustomobject]@{ id='medium_cap64_eps4'; tier='physics'; fixture='research-medium'; atom_cap=64; surface_points=8; salt_M=0.15; int_diel=4.0; ext_diel=80.0; require_speed=$false },
    [pscustomobject]@{ id='medium_cap64_surface16'; tier='resolution'; fixture='research-medium'; atom_cap=64; surface_points=16; salt_M=0.15; int_diel=1.0; ext_diel=80.0; require_speed=$false }
)
if ($Quick) {
    $Conditions = @($Conditions | Where-Object {
        $_.id -in @('medium_cap64_baseline','medium_cap64_salt0')
    })
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
$samplesRoot = Join-Path $WorkRoot 'samples'
New-FmmCliCleanDir -Path $samplesRoot -WorkRoot $WorkRoot
$sampleCache = @{}
$Rows = New-Object System.Collections.Generic.List[object]

function Get-EmptyCase([string]$Status) {
    return [pscustomobject]@{
        status=$Status; panels=0; iterations=0; final_residual=[double]::NaN
        true_residual=[double]::NaN; energy_kj_mol=[double]::NaN
        matvec_seconds_avg=[double]::NaN; total_seconds=[double]::NaN
        fmm_m2l_mode=''; fmm_near_pairs=0; fmm_far_pairs=0
        fmm_setup_seconds=0.0; fmm_apply_seconds=0.0
        fmm_finite_outputs='false'; tree_vs_direct_matvec_rel=[double]::NaN
        matvec_audit_rows=0; matvec_audit_seconds=0.0
        reaction_energy_source=''; timing_backend=''
    }
}

foreach ($Condition in $Conditions) {
    Write-Host "Scope N condition: $($Condition.id)"
    $sampleKey = "$($Condition.fixture)_cap$($Condition.atom_cap)"
    if (-not $sampleCache.ContainsKey($sampleKey)) {
        $sourcePqr = Resolve-FmmCliFixture -Fixture $Condition.fixture `
            -RepoRoot $RepoRoot
        $sampledPqr = Join-Path $samplesRoot "$sampleKey.pqr"
        $sampleInfo = Write-FmmCliSampledPqr -InputPath $sourcePqr `
            -OutputPath $sampledPqr -Cap $Condition.atom_cap
        $sampleCache[$sampleKey] = [pscustomobject]@{
            path = $sampledPqr
            input_atoms = $sampleInfo.input_atoms
            sampled_atoms = $sampleInfo.sampled_atoms
        }
    }
    $sample = $sampleCache[$sampleKey]
    $conditionRoot = Join-Path $WorkRoot $Condition.id
    New-FmmCliCleanDir -Path $conditionRoot -WorkRoot $WorkRoot
    $directDir = Join-Path $conditionRoot 'direct'
    $fmmDir = Join-Path $conditionRoot 'fmm'
    New-FmmCliCleanDir -Path $directDir -WorkRoot $WorkRoot
    New-FmmCliCleanDir -Path $fmmDir -WorkRoot $WorkRoot
    $common = @{
        RepoRoot = $RepoRoot; PqrPath = $sample.path
        SurfacePoints = $Condition.surface_points; MaxIter = $MaxIter
        Restart = $Restart; Tolerance = $Tolerance; Threads = $Threads
        SaltM = $Condition.salt_M; InternalDielectric = $Condition.int_diel
        ExternalDielectric = $Condition.ext_diel; TemperatureK = $TemperatureK
        OpeningTheta = $OpeningTheta
    }
    $directArg = Write-FmmCliArgFile -RunDir $directDir -Mode direct `
        -FmmMode off -AuditSampleCount 0 @common
    $fmmArg = Write-FmmCliArgFile -RunDir $fmmDir -Mode treecode `
        -FmmMode single-matvec-analytic -AuditSampleCount $AuditSampleCount @common

    $directRun = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $directDir `
        -ArgPath $directArg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $fmmRun = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $fmmDir `
        -ArgPath $fmmArg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $direct = if ($directRun.exit_code -eq 0) {
        Read-FmmCliCaseResult -RunDir $directDir -ExpectedFmmMode off `
            -RequireAudit $false
    } else {
        Get-EmptyCase "exit_$($directRun.exit_code)"
    }
    $fmm = if ($fmmRun.exit_code -eq 0) {
        Read-FmmCliCaseResult -RunDir $fmmDir `
            -ExpectedFmmMode single-matvec-analytic -RequireAudit $true
    } else {
        Get-EmptyCase "exit_$($fmmRun.exit_code)"
    }
    $gate = Get-FmmCliFailureClass -DirectExit $directRun.exit_code `
        -FmmExit $fmmRun.exit_code -Direct $direct -Fmm $fmm `
        -Fixture $Condition.fixture -RequireSpeed $Condition.require_speed `
        -RequireAudit $true
    $delta = $fmm.energy_kj_mol - $direct.energy_kj_mol
    $deltaAbs = [Math]::Abs($delta)
    $deltaRel = $deltaAbs / [Math]::Max(1.0, [Math]::Abs($direct.energy_kj_mol))
    $speedup = $direct.matvec_seconds_avg / $fmm.matvec_seconds_avg
    $performanceWarning = if (-not $Condition.require_speed -and
        (Test-FmmCliFinite $speedup) -and $speedup -lt 0.8) {
        'speedup_below_0.8'
    } else { '' }
    $Rows.Add([pscustomobject]@{
        condition = $Condition.id; tier = $Condition.tier
        fixture = $Condition.fixture; required = $true
        atom_cap = $Condition.atom_cap; input_atoms = $sample.input_atoms
        sampled_atoms = $sample.sampled_atoms
        surface_points = $Condition.surface_points
        salt_M = $Condition.salt_M; int_diel = $Condition.int_diel
        ext_diel = $Condition.ext_diel; temperature_K = $TemperatureK
        opening_theta = $OpeningTheta
        audit_sample_count = $AuditSampleCount
        direct_exit_code = $directRun.exit_code; fmm_exit_code = $fmmRun.exit_code
        direct_timed_out = $directRun.timed_out; fmm_timed_out = $fmmRun.timed_out
        direct_status = $direct.status; fmm_status = $fmm.status
        direct_panels = $direct.panels; fmm_panels = $fmm.panels
        direct_iterations = $direct.iterations; fmm_iterations = $fmm.iterations
        direct_final_residual = $direct.final_residual
        fmm_final_residual = $fmm.final_residual
        direct_true_residual = $direct.true_residual
        fmm_true_residual = $fmm.true_residual
        direct_energy_kj_mol = $direct.energy_kj_mol
        fmm_energy_kj_mol = $fmm.energy_kj_mol
        energy_delta_kj_mol = $delta; energy_delta_abs_kj_mol = $deltaAbs
        energy_delta_rel = $deltaRel
        direct_matvec_seconds_avg = $direct.matvec_seconds_avg
        fmm_matvec_seconds_avg = $fmm.matvec_seconds_avg
        matvec_speedup = $speedup
        tree_vs_direct_matvec_rel = $fmm.tree_vs_direct_matvec_rel
        matvec_audit_rows = $fmm.matvec_audit_rows
        matvec_audit_seconds = $fmm.matvec_audit_seconds
        fmm_m2l_mode = $fmm.fmm_m2l_mode
        fmm_near_pairs = $fmm.fmm_near_pairs; fmm_far_pairs = $fmm.fmm_far_pairs
        fmm_setup_seconds = $fmm.fmm_setup_seconds
        fmm_apply_seconds = $fmm.fmm_apply_seconds
        direct_solver_seconds = $direct.total_seconds
        fmm_solver_seconds = $fmm.total_seconds
        direct_wall_seconds = $directRun.wall_seconds
        fmm_wall_seconds = $fmmRun.wall_seconds
        performance_warning = $performanceWarning
        scope_n_gate = $gate
    })
}

$Rows | Export-Csv -LiteralPath $OutCsv -NoTypeInformation -Encoding ASCII
$passRows = @($Rows | Where-Object { $_.scope_n_gate -eq 'pass' })
$warningRows = @($Rows | Where-Object {
    -not [string]::IsNullOrWhiteSpace($_.performance_warning)
})
$timeoutRows = @($Rows | Where-Object {
    $_.scope_n_gate -eq 'timeout_not_verified'
})
$failureRows = @($Rows | Where-Object { $_.scope_n_gate -ne 'pass' })
$overallGate = if ($failureRows.Count -eq 0 -and $Rows.Count -eq $Conditions.Count) {
    'pass'
} elseif ($timeoutRows.Count -gt 0 -and
    $timeoutRows.Count -eq $failureRows.Count) {
    'timeout_not_verified'
} else {
    'fail'
}
$failureClasses = @($failureRows | ForEach-Object { $_.scope_n_gate } |
    Sort-Object -Unique) -join ';'
$decision = [pscustomobject]@{
    scope_n_gate = $overallGate
    quick_mode = $Quick.IsPresent
    required_rows = $Conditions.Count
    rows_written = $Rows.Count
    pass_rows = $passRows.Count
    warning_rows = $warningRows.Count
    timeout_rows = $timeoutRows.Count
    failure_classes = $failureClasses
    production_default = 'grid'
    fmm_status = 'hidden_experimental'
    next_scope = if ($overallGate -eq 'pass' -and -not $Quick) {
        'scope_o_promotion_readiness_repeatability_cap256'
    } else {
        'hold_and_diagnose_scope_n_failures'
    }
}
$decision | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation -Encoding ASCII

$Summary = New-Object System.Collections.Generic.List[string]
$Summary.Add('TITAN FMM/M2L Scope N Physics, Resolution, And Scale Matrix')
$Summary.Add('============================================================')
$Summary.Add('')
$Summary.Add("quick_mode: $($Quick.IsPresent)")
$Summary.Add('scope: hidden_cli_physics_resolution_scale_matrix')
$Summary.Add('hidden_option: -pbe-tabi-fmm-m2l=single-matvec-analytic')
$Summary.Add('production_default: grid')
$Summary.Add('fmm_status: hidden_experimental')
$Summary.Add("audit_mode: sampled")
$Summary.Add("audit_sample_count: $AuditSampleCount")
$Summary.Add("opening_theta: $OpeningTheta")
$Summary.Add("case_timeout_seconds: $CaseTimeoutSeconds")
$Summary.Add("required_rows: $($Conditions.Count)")
$Summary.Add("rows_written: $($Rows.Count)")
$Summary.Add("pass_rows: $($passRows.Count)")
$Summary.Add("warning_rows: $($warningRows.Count)")
$Summary.Add("timeout_rows: $($timeoutRows.Count)")
$Summary.Add("scope_n_gate: $overallGate")
$Summary.Add("failure_classes: $failureClasses")
$Summary.Add('')
foreach ($row in $Rows) {
    $Summary.Add("row: condition=$($row.condition) fixture=$($row.fixture) atoms=$($row.sampled_atoms) surface=$($row.surface_points) salt=$($row.salt_M) eps=$($row.int_diel)/$($row.ext_diel) panels=$($row.fmm_panels) direct_energy=$($row.direct_energy_kj_mol) fmm_energy=$($row.fmm_energy_kj_mol) energy_delta_abs=$($row.energy_delta_abs_kj_mol) fmm_residual=$($row.fmm_true_residual) audit_rel=$($row.tree_vs_direct_matvec_rel) audit_rows=$($row.matvec_audit_rows) speedup=$($row.matvec_speedup) warning=$($row.performance_warning) gate=$($row.scope_n_gate)")
}
$Summary.Add('')
$Summary.Add('acceptance: all required rows must pass status, finite telemetry, residual, sampled direct-operator audit, and fixture energy gates; cap-128 medium and large additionally require matvec speedup >= 1.2.')
$Summary.Add('timeout_policy: timeout_not_verified is evidence missing, not a solver failure classification.')
$Summary.Add('next_gate: Scope O repeatability and cap-256 stress only after the full Scope N matrix passes.')
$Summary.Add("matrix_csv: $OutCsv")
$Summary.Add("decision_csv: $DecisionCsv")
Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

if ($overallGate -ne 'pass') {
    throw "Scope N gate did not pass: $overallGate ($failureClasses)"
}
Write-Host "FMM/M2L Scope N matrix complete: $OutCsv"
