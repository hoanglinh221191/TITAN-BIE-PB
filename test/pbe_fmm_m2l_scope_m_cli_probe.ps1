param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [switch]$Quick,
    [string[]]$Fixtures = @('research-medium','research-large'),
    [int]$AtomCap = 128,
    [int]$SurfacePoints = 8,
    [int]$MaxIter = 120,
    [int]$Restart = 24,
    [double]$Tolerance = 1.0e-6,
    [int]$Threads = 4,
    [int]$CaseTimeoutSeconds = 0,
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_m')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutCsv = Join-Path $WorkRoot 'fmm_m2l_scope_m_matrix.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_m_summary.dat'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_M_REPORT.txt'

if ($Quick) {
    $Fixtures = @('research-medium')
    $AtomCap = [Math]::Min($AtomCap, 64)
    $MaxIter = [Math]::Min($MaxIter, 80)
}
if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}
$TitanExe = (Resolve-Path -LiteralPath $TitanExe).Path
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$Rows = New-Object System.Collections.Generic.List[object]
$MediumPassed = $false
foreach ($Fixture in $Fixtures) {
    if ($Fixture -eq 'research-large' -and -not $MediumPassed) {
        continue
    }
    $fixtureRoot = Join-Path $WorkRoot (Get-FmmCliSafeName $Fixture)
    New-FmmCliCleanDir -Path $fixtureRoot -WorkRoot $WorkRoot
    $sourcePqr = Resolve-FmmCliFixture -Fixture $Fixture -RepoRoot $RepoRoot
    $sampledPqr = Join-Path $fixtureRoot 'sampled.pqr'
    $sample = Write-FmmCliSampledPqr -InputPath $sourcePqr `
        -OutputPath $sampledPqr -Cap $AtomCap
    $directDir = Join-Path $fixtureRoot 'direct'
    $fmmDir = Join-Path $fixtureRoot 'fmm'
    New-FmmCliCleanDir -Path $directDir -WorkRoot $WorkRoot
    New-FmmCliCleanDir -Path $fmmDir -WorkRoot $WorkRoot
    $common = @{
        RepoRoot = $RepoRoot; PqrPath = $sampledPqr
        SurfacePoints = $SurfacePoints; MaxIter = $MaxIter
        Restart = $Restart; Tolerance = $Tolerance; Threads = $Threads
        SaltM = 0.15; InternalDielectric = 1.0
        ExternalDielectric = 80.0; TemperatureK = 298.15
    }
    $directArg = Write-FmmCliArgFile -RunDir $directDir -Mode direct `
        -FmmMode off -AuditSampleCount 0 @common
    $fmmArg = Write-FmmCliArgFile -RunDir $fmmDir -Mode treecode `
        -FmmMode single-matvec-analytic -AuditSampleCount 0 @common
    $directRun = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $directDir `
        -ArgPath $directArg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $fmmRun = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $fmmDir `
        -ArgPath $fmmArg -CaseTimeoutSeconds $CaseTimeoutSeconds
    $direct = if ($directRun.exit_code -eq 0) {
        Read-FmmCliCaseResult -RunDir $directDir -ExpectedFmmMode off `
            -RequireAudit $false
    } else { $null }
    $fmm = if ($fmmRun.exit_code -eq 0) {
        Read-FmmCliCaseResult -RunDir $fmmDir `
            -ExpectedFmmMode single-matvec-analytic -RequireAudit $false
    } else { $null }
    $gate = Get-FmmCliFailureClass -DirectExit $directRun.exit_code `
        -FmmExit $fmmRun.exit_code -Direct $direct -Fmm $fmm `
        -Fixture $Fixture -RequireSpeed $true -RequireAudit $false
    if ($Fixture -eq 'research-medium') { $MediumPassed = $gate -eq 'pass' }
    if ($null -eq $direct -or $null -eq $fmm) {
        throw "Scope M run did not produce readable outputs for $Fixture ($gate)."
    }
    $delta = $fmm.energy_kj_mol - $direct.energy_kj_mol
    $deltaAbs = [Math]::Abs($delta)
    $deltaRel = $deltaAbs / [Math]::Max(1.0, [Math]::Abs($direct.energy_kj_mol))
    $speedup = $direct.matvec_seconds_avg / $fmm.matvec_seconds_avg
    $Rows.Add([pscustomobject]@{
        fixture = $Fixture; atom_cap = $AtomCap
        input_atoms = $sample.input_atoms; sampled_atoms = $sample.sampled_atoms
        direct_exit_code = $directRun.exit_code; fmm_exit_code = $fmmRun.exit_code
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
        matvec_speedup = $speedup; fmm_m2l_mode = $fmm.fmm_m2l_mode
        fmm_near_pairs = $fmm.fmm_near_pairs; fmm_far_pairs = $fmm.fmm_far_pairs
        fmm_setup_seconds = $fmm.fmm_setup_seconds
        fmm_apply_seconds = $fmm.fmm_apply_seconds
        direct_seconds = $direct.total_seconds; fmm_seconds = $fmm.total_seconds
        scope_m_gate = $gate
    })
}

$Rows | Export-Csv -LiteralPath $OutCsv -NoTypeInformation -Encoding ASCII
$passRows = @($Rows | Where-Object { $_.scope_m_gate -eq 'pass' })
$Summary = New-Object System.Collections.Generic.List[string]
$Summary.Add('TITAN FMM/M2L Scope M CLI Probe Summary')
$Summary.Add('=======================================')
$Summary.Add('')
$Summary.Add("quick_mode: $($Quick.IsPresent)")
$Summary.Add('scope: hidden_cli_pb_solve_energy_gate')
$Summary.Add('hidden_option: -pbe-tabi-fmm-m2l=single-matvec-analytic')
$Summary.Add('production_default: unchanged_grid_default')
$Summary.Add("atom_cap: $AtomCap")
$Summary.Add("surface_points: $SurfacePoints")
$Summary.Add("rows: $($Rows.Count)")
$Summary.Add("scope_m_gate_pass_rows: $($passRows.Count)")
foreach ($row in $Rows) {
    $Summary.Add("row: fixture=$($row.fixture) atoms=$($row.sampled_atoms) direct_status=$($row.direct_status) fmm_status=$($row.fmm_status) direct_energy=$($row.direct_energy_kj_mol) fmm_energy=$($row.fmm_energy_kj_mol) energy_delta_abs=$($row.energy_delta_abs_kj_mol) direct_residual=$($row.direct_true_residual) fmm_residual=$($row.fmm_true_residual) speedup=$($row.matvec_speedup) fmm_mode=$($row.fmm_m2l_mode) fmm_far_pairs=$($row.fmm_far_pairs) gate=$($row.scope_m_gate)")
}
$Summary.Add('')
$Summary.Add('acceptance: direct/FMM titan.exe runs exit 0, status ok, finite energy/residual, FMM telemetry present, energy delta within fixture gate, and matvec speedup >= 1.2 on medium.')
$Summary.Add('next_gate: Scope N hidden broader energy/regression matrix before any production exposure.')
$Summary.Add("matrix_csv: $OutCsv")
Set-Content -LiteralPath $SummaryPath -Value $Summary -Encoding ASCII
Set-Content -LiteralPath $ReportPath -Value $Summary -Encoding ASCII

if (@($Rows | Where-Object {
        $_.fixture -eq 'research-medium' -and $_.scope_m_gate -eq 'pass'
    }).Count -eq 0) {
    throw 'Scope M medium gate did not pass.'
}
Write-Host "FMM/M2L Scope M CLI probe complete: $OutCsv"
