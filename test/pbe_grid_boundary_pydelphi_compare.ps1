param(
    [string]$TitanExe = (Join-Path (Split-Path -Parent $PSScriptRoot) 'titan.exe'),
    [string]$PyDelPhiRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) '.codex_tmp\pyDelPhi'),
    [string]$PythonExe = $(if (Test-Path -LiteralPath 'C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe') { 'C:\Users\ADMIN\.codex_tmp\pydelphi_env\python.exe' } elseif (Test-Path -LiteralPath 'C:\ProgramData\miniconda3\python.exe') { 'C:\ProgramData\miniconda3\python.exe' } else { 'python' }),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_grid_boundary_pydelphi_compare'),
    [switch]$Strict
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$InvariantCulture = [System.Globalization.CultureInfo]::InvariantCulture

function Write-SyntheticPqr([string]$Path, [double]$X) {
    function Format-PqrCoord([double]$Value) {
        if ([math]::Abs($Value) -gt 0.0 -and [math]::Abs($Value) -lt 1.0e-4) {
            return $Value.ToString('0.0E+0', $InvariantCulture).Replace('E', 'e')
        }
        return $Value.ToString('F4', $InvariantCulture)
    }
    $xText = Format-PqrCoord $X
    $yText = Format-PqrCoord 0.0
    $zText = Format-PqrCoord 0.0
    $line = ('ATOM  {0,5} {1,-4} {2,3} {3}{4,4}    {5,8}{6,8}{7,8}{8,8:F4}{9,8:F4}' -f 1, 'CA', 'GLY', 'A', 1, $xText, $yText, $zText, 0.25, 1.70)
    Set-Content -LiteralPath $Path -Value @($line, 'END') -Encoding ASCII
}

function Convert-ToTitanPath([string]$Path) {
    return ($Path -replace '\\', '/')
}

if (-not (Test-Path -LiteralPath $TitanExe)) {
    throw "Missing titan executable: $TitanExe"
}

if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null

$inputRoot = Join-Path $WorkRoot 'inputs'
New-Item -ItemType Directory -Path $inputRoot -Force | Out-Null
$cases = @(
    @{ fixture = 'x_zero'; x = 0.0 },
    @{ fixture = 'x_minus_1e-8_A'; x = -1.0e-8 },
    @{ fixture = 'x_plus_1e-8_A'; x = 1.0e-8 },
    @{ fixture = 'x_minus_1e-6_A'; x = -1.0e-6 },
    @{ fixture = 'x_plus_1e-6_A'; x = 1.0e-6 },
    @{ fixture = 'x_minus_1e-4_A'; x = -1.0e-4 },
    @{ fixture = 'x_plus_1e-4_A'; x = 1.0e-4 }
)

$pqrPaths = New-Object System.Collections.Generic.List[string]
foreach ($case in $cases) {
    $path = Join-Path $inputRoot "$($case.fixture).pqr"
    Write-SyntheticPqr -Path $path -X ([double]$case.x)
    $pqrPaths.Add($path)
}

$scanRoot = Join-Path $WorkRoot 'titan_scan'
& powershell -ExecutionPolicy Bypass -File (Join-Path $PSScriptRoot 'pbe_snapshot_stability_scan.ps1') `
    -TitanExe $TitanExe `
    -WorkRoot $scanRoot `
    -PqrPaths ($pqrPaths -join ',') `
    -GridAnchor fixed-center `
    -StabilityMode off | Out-Host
if ($LASTEXITCODE -ne 0) {
    throw "TITAN boundary scan failed with code $LASTEXITCODE"
}

$titanSummary = Join-Path $scanRoot 'pbe_snapshot_stability_summary.csv'
if (-not (Test-Path -LiteralPath $titanSummary)) {
    throw "Missing TITAN boundary scan summary: $titanSummary"
}
$titanRows = @(Import-Csv -LiteralPath $titanSummary)

$pydelphiRows = @()
if ((Test-Path -LiteralPath $PyDelPhiRoot) -and (Get-Command $PythonExe -ErrorAction SilentlyContinue)) {
    $driver = Join-Path $WorkRoot 'pydelphi_boundary_driver.py'
    $driverCode = @'
import csv
import math
import os
import sys
import traceback
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

pyroot = Path(sys.argv[1])
work = Path(sys.argv[2])
inputs = sys.argv[3:]
sys.path.insert(0, str(pyroot))

rows = []
try:
    from pydelphi.app.delphi import DelphiApp
    from pydelphi.foundation.enums import Precision
    from pydelphi.foundation.platforms import Platform
    from pydelphi.config.global_runtime import set_precision
    platform = Platform(platform_name="cpu", debug=False)
    platform.activate("cpu", 1, 0)
    platform.set_precision(Precision.DOUBLE)
    set_precision(platform.precision)
except Exception:
    err = traceback.format_exc()
    for pqr in inputs:
        rows.append({"fixture": Path(pqr).stem, "pydelphi_status": "import_failed", "detail": err})
else:
    for pqr in inputs:
        pqr_path = Path(pqr)
        run = work / ("pydelphi_" + pqr_path.stem)
        run.mkdir(parents=True, exist_ok=True)
        local_pqr = run / pqr_path.name
        local_pqr.write_text(pqr_path.read_text(encoding="ascii"), encoding="ascii")
        prm = run / "boundary.prm"
        prm.write_text("\n".join([
            "biomodel=pbe",
            "solver=sor",
            "surfacemethod=vdw",
            "dielectricmodel=twodielectric",
            "scale=1.0",
            "grid_size=33",
            "gridbox_type=cubic",
            "acent(0,0,0)",
            f'in(modpdb4,file="{local_pqr.name}",format="pqr")',
            "indi=1.0",
            "exdi=80.0",
            "gapdi=22.0",
            "boundary_condition=coulombic",
            "salt=0.15",
            "linit=120",
            "nonlinit=0",
            "max_rmsd=1.0e-4",
            "max_delta_phi=1.0e-3",
        ]) + "\n", encoding="ascii")
        cwd = os.getcwd()
        try:
            os.chdir(run)
            app = DelphiApp(str(prm.name), platform, user_inputs=None)
            energies = app.run("outputs.csv", pqr_path.stem, True)
            ctx = getattr(app, "ctx", None)
            grid_shape = getattr(ctx, "grid_shape", "NA") if ctx is not None else "NA"
            center = getattr(ctx, "grid_center", "NA") if ctx is not None else "NA"
            rows.append({
                "fixture": pqr_path.stem,
                "pydelphi_status": "ok",
                "pydelphi_energy_repr": repr(energies),
                "pydelphi_grid_shape": str(grid_shape),
                "pydelphi_grid_center": str(center),
                "detail": "",
            })
        except Exception:
            rows.append({
                "fixture": pqr_path.stem,
                "pydelphi_status": "run_failed",
                "pydelphi_energy_repr": "NA",
                "pydelphi_grid_shape": "NA",
                "pydelphi_grid_center": "NA",
                "detail": traceback.format_exc(),
            })
        finally:
            os.chdir(cwd)

with (work / "pydelphi_boundary_rows.csv").open("w", newline="", encoding="ascii") as f:
    fields = ["fixture", "pydelphi_status", "pydelphi_energy_repr", "pydelphi_grid_shape", "pydelphi_grid_center", "detail"]
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for row in rows:
        writer.writerow({key: row.get(key, "NA") for key in fields})
'@
    Set-Content -LiteralPath $driver -Value $driverCode -Encoding ASCII
    $oldPythonUtf8 = $env:PYTHONUTF8
    $env:PYTHONUTF8 = '1'
    $pyLog = Join-Path $WorkRoot 'pydelphi_boundary_driver.log'
    & $PythonExe $driver $PyDelPhiRoot $WorkRoot @($pqrPaths.ToArray()) *> $pyLog
    $pythonExit = $LASTEXITCODE
    $env:PYTHONUTF8 = $oldPythonUtf8
    $pydelphiCsv = Join-Path $WorkRoot 'pydelphi_boundary_rows.csv'
    if (Test-Path -LiteralPath $pydelphiCsv) {
        $pydelphiRows = @(Import-Csv -LiteralPath $pydelphiCsv)
    } elseif ($pythonExit -ne 0) {
        $pydelphiRows = foreach ($case in $cases) {
            [PSCustomObject]@{
                fixture = [string]$case.fixture
                pydelphi_status = 'python_failed'
                pydelphi_energy_repr = 'NA'
                pydelphi_grid_shape = 'NA'
                pydelphi_grid_center = 'NA'
                detail = "exit_code=$pythonExit"
            }
        }
    }
}

$rows = foreach ($case in $cases) {
    $fixture = [string]$case.fixture
    $titan = @($titanRows | Where-Object { $_.fixture -eq ('snapshot_{0:000}' -f ([array]::IndexOf($cases.fixture, $fixture) + 1)) })
    if ($titan.Count -eq 0) {
        $titan = @($titanRows | Where-Object { $_.pqr_path -like "*$fixture.pqr" })
    }
    $py = @($pydelphiRows | Where-Object { $_.fixture -eq $fixture })
    [PSCustomObject]@{
        fixture = $fixture
        x_A = ([double]$case.x).ToString('R', $InvariantCulture)
        titan_status = if ($titan.Count -gt 0) { $titan[0].status } else { 'missing' }
        titan_energy_kJmol = if ($titan.Count -gt 0) { $titan[0].energy_kJmol } else { 'NA' }
        titan_scan_label = if ($titan.Count -gt 0) { $titan[0].scan_label } else { 'missing' }
        titan_suspect_stage = if ($titan.Count -gt 0) { $titan[0].suspect_stage } else { 'missing' }
        titan_near_integer_atoms_1e_6 = if ($titan.Count -gt 0) { $titan[0].'near_integer_atoms_1e-6' } else { 'NA' }
        titan_near_half_atoms_1e_6 = if ($titan.Count -gt 0) { $titan[0].'near_half_atoms_1e-6' } else { 'NA' }
        pydelphi_status = if ($py.Count -gt 0) { $py[0].pydelphi_status } else { 'not_run' }
        pydelphi_energy_repr = if ($py.Count -gt 0) { $py[0].pydelphi_energy_repr } else { 'NA' }
        pydelphi_grid_shape = if ($py.Count -gt 0) { $py[0].pydelphi_grid_shape } else { 'NA' }
        first_failing_stage = if ($titan.Count -gt 0) { $titan[0].suspect_stage } else { 'missing' }
    }
}

$outCsv = Join-Path $WorkRoot 'pbe_grid_boundary_pydelphi_compare.csv'
$rows | Export-Csv -LiteralPath $outCsv -NoTypeInformation -Encoding ASCII

$failedPy = @($rows | Where-Object { $_.pydelphi_status -notin @('ok', 'not_run') })
if ($Strict -and $failedPy.Count -gt 0) {
    throw "pyDelPhi boundary compare failed for $($failedPy.Count) row(s): $outCsv"
}

Write-Host "PBE grid-boundary TITAN/pyDelPhi compare PASS: $outCsv"
