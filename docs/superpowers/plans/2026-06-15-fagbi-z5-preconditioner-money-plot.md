# FAGBI Z5 Preconditioner Money Plot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the Scope Z5 benchmark/report pipeline that demonstrates paper-style block-diagonal preconditioning stalling on dense TABI/FMM systems while TITAN sparse near-field ILU(0) reaches strict convergence.

**Architecture:** Reuse the existing FMM/M2L CLI helper and Scope X harness patterns. The new Z5 script runs paired preconditioner cases on the same sampled Barnase cap-256 fixture, parses `pbe_tabi_profile.csv` including `residual_history`, exports matrix/trace CSVs, and writes a paper-facing report plus a publication-ready SVG plot. No production solver defaults change.

**Tech Stack:** PowerShell harnesses, existing `titan.exe -a arg.in` CLI, `pbe_tabi_profile.csv`, `pbe_energy.csv`, `pbe_timing.csv`, ASCII/CSV report artifacts, optional inline SVG generation from PowerShell.

---

## File Structure

- Create `test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1`
  - Scope Z5 harness, fixture setup, case execution, trace parsing, matrix export, decision report, and SVG plot generation.
- Modify `test/pbe_fmm_m2l_cli_helpers.ps1`
  - Add a helper to parse `residual_history` into iteration/value pairs if one does not already exist.
- Modify `TITAN_2110_13778_FAGBI_ANALYSIS.md`
  - Mark Scope Z5 implementation/report path once the harness exists.
- Modify `TITAN_FMM_M2L_ROADMAP.md`
  - Add Scope Z5 result ledger after quick/full verification.
- Modify `TITAN_PBE_TABI_ROADMAP.md`
  - Mirror the Scope Z5 scientific conclusion for TABI/PBE.
- Modify `.gitignore`
  - Whitelist the new Scope Z5 script, report, and this plan if needed.
- Generated outputs:
  - `test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_matrix.csv`
  - `test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_residual_trace.csv`
  - `test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_decision.csv`
  - `test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_summary.dat`
  - `test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_money_plot.svg`
  - `TITAN_FMM_M2L_SCOPE_Z5_PRECONDITIONER_MONEY_PLOT_REPORT.txt`

---

### Task 1: Add Residual-History Parser Helper

**Files:**
- Modify: `test/pbe_fmm_m2l_cli_helpers.ps1`
- Test: existing profile outputs via Scope Z5 quick script

- [ ] **Step 1: Add a parser function near the other `Convert-FmmCli*` helpers**

Add:

```powershell
function Convert-FmmCliResidualHistory([string]$History) {
    $rows = New-Object System.Collections.Generic.List[object]
    if ([string]::IsNullOrWhiteSpace($History)) {
        return @()
    }
    $tokens = @($History -split '[;| ]+' | Where-Object {
        -not [string]::IsNullOrWhiteSpace($_)
    })
    $iteration = 0
    foreach ($token in $tokens) {
        $valueText = $token
        if ($token -match '^[0-9]+[:=](.+)$') {
            $iteration = [int]($token -replace '[:=].+$', '')
            $valueText = $Matches[1]
        }
        $value = Convert-FmmCliDouble $valueText
        if (Test-FmmCliFinite $value) {
            $rows.Add([pscustomobject]@{
                iteration = $iteration
                residual = $value
            }) | Out-Null
            $iteration += 1
        }
    }
    return @($rows)
}
```

- [ ] **Step 2: Add residual history to `Read-FmmCliCaseResult` return object**

Inside the returned `[pscustomobject]`, add:

```powershell
residual_history = [string]$profile.residual_history
```

- [ ] **Step 3: Verify helper syntax**

Run:

```powershell
powershell -NoProfile -Command "& { $null = [scriptblock]::Create((Get-Content -Raw test\pbe_fmm_m2l_cli_helpers.ps1)); 'parse ok' }"
```

Expected: `parse ok`.

---

### Task 2: Create Scope Z5 Harness

**Files:**
- Create: `test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1`
- Uses: `test/pbe_fmm_m2l_cli_helpers.ps1`

- [ ] **Step 1: Create the script header and fixed defaults**

Create the file with:

```powershell
param(
    [string]$TitanExe = (Join-Path $PSScriptRoot '..\titan.exe'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\fmm_m2l_scope_z5'),
    [switch]$Quick,
    [int]$CaseTimeoutSeconds = 900,
    [int]$Threads = 4
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'pbe_fmm_m2l_cli_helpers.ps1')

$RepoRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$WorkRoot = (New-FmmCliCleanDir $WorkRoot $RepoRoot)
$MatrixCsv = Join-Path $WorkRoot 'fmm_m2l_scope_z5_matrix.csv'
$TraceCsv = Join-Path $WorkRoot 'fmm_m2l_scope_z5_residual_trace.csv'
$DecisionCsv = Join-Path $WorkRoot 'fmm_m2l_scope_z5_decision.csv'
$SummaryPath = Join-Path $WorkRoot 'fmm_m2l_scope_z5_summary.dat'
$PlotPath = Join-Path $WorkRoot 'fmm_m2l_scope_z5_money_plot.svg'
$ReportPath = Join-Path $RepoRoot 'TITAN_FMM_M2L_SCOPE_Z5_PRECONDITIONER_MONEY_PLOT_REPORT.txt'
```

- [ ] **Step 2: Resolve and sample Barnase**

Add:

```powershell
$fixture = Resolve-FmmCliFixture -Name 'barnase' -RepoRoot $RepoRoot
if ($fixture.status -ne 'ok') {
    [pscustomobject]@{
        scope_z5_decision = 'skipped_missing_fixture'
        reason = $fixture.status
    } | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation
    Set-Content -LiteralPath $SummaryPath -Encoding ASCII -Value @(
        'scope_z5_decision=skipped_missing_fixture',
        "fixture_status=$($fixture.status)"
    )
    throw "Barnase fixture could not be resolved: $($fixture.status)"
}

$atomCap = if ($Quick) { 128 } else { 256 }
$sampleDir = Join-Path $WorkRoot 'fixture'
New-Item -ItemType Directory -Force -Path $sampleDir | Out-Null
$sampleInfo = Write-FmmCliSampledPqr -InputPath $fixture.path `
    -OutputPath (Join-Path $sampleDir 'barnase_sampled.pqr') `
    -AtomCap $atomCap
```

- [ ] **Step 3: Define preconditioner cases**

Add:

```powershell
$cases = @(
    [pscustomobject]@{
        label='diagonal'; precond='diagonal'; ilu_cap=64; block_cap=32;
        overlap_cap=64; max_iter=720; restart=128; tolerance=1.0e-6
    },
    [pscustomobject]@{
        label='block_diagonal'; precond='fmm-block-jacobi'; ilu_cap=64;
        block_cap=32; overlap_cap=64; max_iter=720; restart=128;
        tolerance=1.0e-6
    },
    [pscustomobject]@{
        label='overlap_schwarz'; precond='fmm-overlap-schwarz'; ilu_cap=64;
        block_cap=32; overlap_cap=64; max_iter=720; restart=128;
        tolerance=1.0e-6
    },
    [pscustomobject]@{
        label='sparse_ilu0'; precond='fmm-nearfield-ilu0'; ilu_cap=64;
        block_cap=32; overlap_cap=64; max_iter=720; restart=128;
        tolerance=1.0e-6
    }
)
if ($Quick) {
    $cases = @($cases | Where-Object { $_.label -in @('block_diagonal','sparse_ilu0') })
}
```

- [ ] **Step 4: Add case runner**

Add:

```powershell
function Invoke-Z5Case([object]$Case) {
    $runDir = Join-Path $WorkRoot $Case.label
    New-FmmCliCleanDir $runDir $WorkRoot | Out-Null
    $argPath = Write-FmmCliArgFile -RunDir $runDir -RepoRoot $RepoRoot `
        -PqrPath $sampleInfo.path -Backend 'tabi' -Mode 'treecode' `
        -FmmMode 'single-matvec-analytic' -SurfacePoints 8 `
        -MaxIter $Case.max_iter -Restart $Case.restart `
        -Tolerance $Case.tolerance -Threads $Threads -SaltM 0.15 `
        -InternalDielectric 1.0 -ExternalDielectric 80.0 `
        -TemperatureK 298.15 -OpeningTheta 0.6 -AuditSampleCount 512 `
        -Preconditioner $Case.precond `
        -PreconditionerBlockPanelCap $Case.block_cap `
        -PreconditionerOverlapPanelCap $Case.overlap_cap `
        -PreconditionerIluNeighborCap $Case.ilu_cap `
        -PreconditionerIluDropTol 0.0
    $run = Invoke-FmmCliTitanCase -TitanExe $TitanExe -RunDir $runDir `
        -ArgPath $argPath -CaseTimeoutSeconds $CaseTimeoutSeconds
    $result = $null
    if ($run.exit_code -eq 0 -and -not $run.timed_out) {
        $result = Read-FmmCliCaseResult -RunDir $runDir `
            -ExpectedFmmMode 'single-matvec-analytic' -RequireAudit $true
    }
    return [pscustomobject]@{
        case = $Case
        run = $run
        result = $result
        run_dir = $runDir
    }
}
```

- [ ] **Step 5: Run cases and export matrix**

Add:

```powershell
$records = New-Object System.Collections.Generic.List[object]
foreach ($case in $cases) {
    $records.Add((Invoke-Z5Case $case)) | Out-Null
}

$rows = foreach ($record in $records) {
    $r = $record.result
    [pscustomobject]@{
        fixture='barnase'
        atom_cap=$atomCap
        sampled_atoms=$sampleInfo.sampled_atoms
        preconditioner_label=$record.case.label
        preconditioner_requested=$record.case.precond
        status=if ($null -ne $r) { $r.status } elseif ($record.run.timed_out) { 'timeout' } else { 'failed' }
        exit_code=$record.run.exit_code
        timed_out=$record.run.timed_out
        wall_seconds=$record.run.wall_seconds
        peak_working_set_mb=$record.run.peak_working_set_mb
        panels=if ($null -ne $r) { $r.panels } else { 0 }
        iterations=if ($null -ne $r) { $r.iterations } else { 0 }
        true_residual=if ($null -ne $r) { $r.true_residual } else { [double]::NaN }
        energy_kj_mol=if ($null -ne $r) { $r.energy_kj_mol } else { [double]::NaN }
        matvec_seconds_avg=if ($null -ne $r) { $r.matvec_seconds_avg } else { [double]::NaN }
        audit_rel=if ($null -ne $r) { $r.tree_vs_direct_matvec_rel } else { [double]::NaN }
        fmm_far_pairs=if ($null -ne $r) { $r.fmm_far_pairs } else { 0 }
        residual_history=if ($null -ne $r) { $r.residual_history } else { '' }
    }
}
$rows | Export-Csv -LiteralPath $MatrixCsv -NoTypeInformation
```

---

### Task 3: Export Residual Trace And SVG Money Plot

**Files:**
- Modify: `test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1`

- [ ] **Step 1: Convert residual histories to trace rows**

Append:

```powershell
$traceRows = New-Object System.Collections.Generic.List[object]
foreach ($row in $rows) {
    foreach ($point in (Convert-FmmCliResidualHistory ([string]$row.residual_history))) {
        $traceRows.Add([pscustomobject]@{
            fixture=$row.fixture
            atom_cap=$row.atom_cap
            preconditioner_label=$row.preconditioner_label
            iteration=$point.iteration
            residual=$point.residual
            log10_residual=if ($point.residual -gt 0.0) {
                [Math]::Log10($point.residual)
            } else {
                [double]::NaN
            }
        }) | Out-Null
    }
}
$traceRows | Export-Csv -LiteralPath $TraceCsv -NoTypeInformation
```

- [ ] **Step 2: Generate an inline SVG plot**

Append:

```powershell
function Write-Z5SvgPlot([object[]]$TraceRows, [string]$Path) {
    $finite = @($TraceRows | Where-Object {
        Test-FmmCliFinite ([double]$_.log10_residual)
    })
    if ($finite.Count -eq 0) {
        Set-Content -LiteralPath $Path -Encoding ASCII -Value '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="520"><text x="20" y="40">No finite residual trace</text></svg>'
        return
    }
    $width = 900; $height = 520
    $left = 80; $right = 30; $top = 30; $bottom = 70
    $maxIter = [double](@($finite | Measure-Object -Property iteration -Maximum)[0].Maximum)
    $minLog = [double](@($finite | Measure-Object -Property log10_residual -Minimum)[0].Minimum)
    $maxLog = [double](@($finite | Measure-Object -Property log10_residual -Maximum)[0].Maximum)
    if ($maxIter -le 0) { $maxIter = 1 }
    if ([Math]::Abs($maxLog - $minLog) -lt 1.0e-12) { $maxLog += 1; $minLog -= 1 }
    $colors = @{
        diagonal = '#666666'
        block_diagonal = '#d55e00'
        overlap_schwarz = '#0072b2'
        sparse_ilu0 = '#009e73'
    }
    $svg = New-Object System.Collections.Generic.List[string]
    $svg.Add("<svg xmlns=`"http://www.w3.org/2000/svg`" width=`"$width`" height=`"$height`" viewBox=`"0 0 $width $height`">")
    $svg.Add('<rect width="100%" height="100%" fill="white"/>')
    $svg.Add("<line x1=`"$left`" y1=`"$($height-$bottom)`" x2=`"$($width-$right)`" y2=`"$($height-$bottom)`" stroke=`"black`"/>")
    $svg.Add("<line x1=`"$left`" y1=`"$top`" x2=`"$left`" y2=`"$($height-$bottom)`" stroke=`"black`"/>")
    foreach ($label in @('diagonal','block_diagonal','overlap_schwarz','sparse_ilu0')) {
        $points = @($finite | Where-Object { $_.preconditioner_label -eq $label } |
            Sort-Object {[int]$_.iteration})
        if ($points.Count -eq 0) { continue }
        $coords = foreach ($p in $points) {
            $x = $left + (($width - $left - $right) * ([double]$p.iteration / $maxIter))
            $y = $top + (($height - $top - $bottom) * (($maxLog - [double]$p.log10_residual) / ($maxLog - $minLog)))
            ('{0:F2},{1:F2}' -f $x,$y)
        }
        $color = $colors[$label]
        $svg.Add("<polyline fill=`"none`" stroke=`"$color`" stroke-width=`"2.4`" points=`"$($coords -join ' ')`"/>")
    }
    $legendY = 60
    foreach ($label in @('diagonal','block_diagonal','overlap_schwarz','sparse_ilu0')) {
        $color = $colors[$label]
        $svg.Add("<line x1=`"610`" y1=`"$legendY`" x2=`"650`" y2=`"$legendY`" stroke=`"$color`" stroke-width=`"3`"/>")
        $svg.Add("<text x=`"660`" y=`"$($legendY+5)`" font-family=`"Arial`" font-size=`"14`">$label</text>")
        $legendY += 24
    }
    $svg.Add('<text x="390" y="500" font-family="Arial" font-size="16">GMRES iteration</text>')
    $svg.Add('<text x="20" y="280" font-family="Arial" font-size="16" transform="rotate(-90 20,280)">log10 residual</text>')
    $svg.Add('<text x="80" y="24" font-family="Arial" font-size="18" font-weight="bold">Scope Z5 preconditioner convergence</text>')
    $svg.Add('</svg>')
    Set-Content -LiteralPath $Path -Encoding ASCII -Value $svg
}

Write-Z5SvgPlot -TraceRows @($traceRows) -Path $PlotPath
```

---

### Task 4: Decision Classification And Report

**Files:**
- Modify: `test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1`
- Create output: `TITAN_FMM_M2L_SCOPE_Z5_PRECONDITIONER_MONEY_PLOT_REPORT.txt`

- [ ] **Step 1: Add classification logic**

Append:

```powershell
$ilu = @($rows | Where-Object { $_.preconditioner_label -eq 'sparse_ilu0' }) | Select-Object -First 1
$block = @($rows | Where-Object { $_.preconditioner_label -eq 'block_diagonal' }) | Select-Object -First 1
$z5Decision = 'no_decision'
if ($null -eq $ilu -or $null -eq $block) {
    $z5Decision = 'missing_required_rows'
} elseif ($ilu.status -ne 'ok') {
    $z5Decision = 'ilu_not_ok'
} elseif (([double]$ilu.true_residual) -gt 1.0e-6) {
    $z5Decision = 'ilu_not_strict'
} elseif ($block.status -eq 'ok' -and ([double]$block.true_residual) -le 1.0e-6) {
    $z5Decision = 'block_also_strict'
} else {
    $z5Decision = 'money_plot_pass'
}

[pscustomobject]@{
    scope_z5_decision=$z5Decision
    atom_cap=$atomCap
    sampled_atoms=$sampleInfo.sampled_atoms
    matrix=$MatrixCsv
    trace=$TraceCsv
    plot=$PlotPath
} | Export-Csv -LiteralPath $DecisionCsv -NoTypeInformation
```

- [ ] **Step 2: Write the report**

Append:

```powershell
$lines = New-Object System.Collections.Generic.List[string]
$lines.Add('TITAN FMM/M2L Scope Z5 Preconditioner Money Plot Report')
$lines.Add("Decision: $z5Decision")
$lines.Add("Fixture: Barnase sampled cap $atomCap ($($sampleInfo.sampled_atoms) atoms)")
$lines.Add('')
$lines.Add('Rows:')
foreach ($row in $rows) {
    $lines.Add(('{0}: status={1} iter={2} residual={3:G6} audit={4:G6} wall={5:G5}s peak={6:G5}MB' -f
        $row.preconditioner_label, $row.status, [int]$row.iterations,
        [double]$row.true_residual, [double]$row.audit_rel,
        [double]$row.wall_seconds, [double]$row.peak_working_set_mb))
}
$lines.Add('')
$lines.Add('Publication framing:')
$lines.Add('Purely local block-diagonal preconditioning is expected to stall on dense boundary-integral PB systems, while TITAN sparse near-field ILU(0) captures cross-leaf near-field coupling and restores Krylov convergence.')
$lines.Add('')
$lines.Add("Matrix: $MatrixCsv")
$lines.Add("Trace: $TraceCsv")
$lines.Add("Plot: $PlotPath")
Set-Content -LiteralPath $ReportPath -Encoding ASCII -Value $lines
Set-Content -LiteralPath $SummaryPath -Encoding ASCII -Value @(
    "scope_z5_decision=$z5Decision",
    "fixture=barnase",
    "atom_cap=$atomCap",
    "sampled_atoms=$($sampleInfo.sampled_atoms)",
    "matrix=$MatrixCsv",
    "trace=$TraceCsv",
    "plot=$PlotPath",
    "report=$ReportPath"
)
Write-Host "Scope Z5 report: $ReportPath"
```

---

### Task 5: Whitelist Artifacts And Sync Roadmaps

**Files:**
- Modify: `.gitignore`
- Modify: `TITAN_2110_13778_FAGBI_ANALYSIS.md`
- Modify: `TITAN_FMM_M2L_ROADMAP.md`
- Modify: `TITAN_PBE_TABI_ROADMAP.md`

- [ ] **Step 1: Whitelist the script and report**

Add to `.gitignore` near the FMM scope report whitelist:

```gitignore
!/TITAN_FMM_M2L_SCOPE_Z5_PRECONDITIONER_MONEY_PLOT_REPORT.txt
!/test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1
!/docs/superpowers/plans/2026-06-15-fagbi-z5-preconditioner-money-plot.md
```

- [ ] **Step 2: Add implementation status to `TITAN_2110_13778_FAGBI_ANALYSIS.md`**

Under Scope Z5, after the acceptance list, add:

```markdown
Implementation status:

- Planned in `docs/superpowers/plans/2026-06-15-fagbi-z5-preconditioner-money-plot.md`.
- Harness target: `test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1`.
- Report target: `TITAN_FMM_M2L_SCOPE_Z5_PRECONDITIONER_MONEY_PLOT_REPORT.txt`.
```

- [ ] **Step 3: Add a Scope Z5 ledger entry to `TITAN_FMM_M2L_ROADMAP.md` and `TITAN_PBE_TABI_ROADMAP.md`**

Add:

```markdown
## Scope Z5 - Preconditioner Money Plot

Status: planned.

Scope Z5 turns the sparse near-field ILU(0) result into the paper-facing
preconditioner benchmark. It compares diagonal, paper-analogue Block-Jacobi,
Restricted Additive-Schwarz, and TITAN sparse near-field ILU(0) on the same
Barnase cap-256 FMM/M2L operator and exports matrix, residual trace, and SVG
plot artifacts. Grid remains the production PBE default.
```

---

### Task 6: Verification

**Files:**
- Test: `test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1`

- [ ] **Step 1: Static parse checks**

Run:

```powershell
powershell -NoProfile -Command "& { $null = [scriptblock]::Create((Get-Content -Raw test\pbe_fmm_m2l_cli_helpers.ps1)); $null = [scriptblock]::Create((Get-Content -Raw test\pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1)); 'parse ok' }"
```

Expected: `parse ok`.

- [ ] **Step 2: Quick Z5 smoke**

Run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1 -TitanExe .\titan.exe -Quick -CaseTimeoutSeconds 900
```

Expected:

- `test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_matrix.csv` exists.
- `test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_residual_trace.csv` exists.
- `test/tmp/fmm_m2l_scope_z5/fmm_m2l_scope_z5_money_plot.svg` exists.
- `TITAN_FMM_M2L_SCOPE_Z5_PRECONDITIONER_MONEY_PLOT_REPORT.txt` exists.

- [ ] **Step 3: Full Z5 gate**

Run:

```powershell
mingw32-make PROFILE=native TARGET=titan.exe all -j1
powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1 -TitanExe .\titan.exe -CaseTimeoutSeconds 2400
```

Expected:

- Decision is `money_plot_pass`, `ilu_not_strict`, `block_also_strict`, or `missing_required_rows`; no silent pass.
- If `money_plot_pass`, matrix shows `sparse_ilu0` strict and `block_diagonal` non-strict.
- Plot is regenerated from CSV without manual steps.

- [ ] **Step 4: Whitespace check**

Run:

```powershell
git diff --check -- test/pbe_fmm_m2l_cli_helpers.ps1 test/pbe_fmm_m2l_scope_z5_preconditioner_money_plot.ps1 TITAN_2110_13778_FAGBI_ANALYSIS.md TITAN_FMM_M2L_ROADMAP.md TITAN_PBE_TABI_ROADMAP.md .gitignore docs/superpowers/plans/2026-06-15-fagbi-z5-preconditioner-money-plot.md
```

Expected: no whitespace errors. CRLF warnings are acceptable on this Windows checkout.
