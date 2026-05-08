param(
    [string]$PyDelPhiRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) '.codex_tmp\pyDelPhi'),
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_delphi_export_maps'),
    [string]$ManifestPath = (Join-Path $PSScriptRoot 'pbe_oracle_manifest.csv'),
    [string]$GoldenCsv = (Join-Path $PSScriptRoot 'pbe_inputs_local\pbe_delphi_golden.csv'),
    [string]$TitanSampleRoot = (Join-Path $PSScriptRoot 'tmp\pbe_barnase_map_oracle_smoke\runs'),
    [string]$PythonExe = $(if (Test-Path -LiteralPath 'C:\ProgramData\miniconda3\python.exe') { 'C:\ProgramData\miniconda3\python.exe' } else { 'python' }),
    [string]$FullMidpointSalt = '0.15',
    [string[]]$Fixture = @(),
    [switch]$FullMidpointEpsilon,
    [switch]$FullVdwmsState,
    [switch]$FullBoundaryOracle,
    [switch]$FullOperatorOracle,
    [switch]$FullChargeGrid,
    [switch]$FullFieldOracle,
    [switch]$FullIscmEnergyOracle,
    [switch]$FullGaussianOracle,
    [switch]$FullRpbOracle,
    [switch]$FullFocusProbeOracle,
    [switch]$FullFieldMap,
    [switch]$NwtOnly,
    [switch]$Strict
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Get-HashText([string]$Path) {
    if ([string]::IsNullOrWhiteSpace($Path) -or -not (Test-Path -LiteralPath $Path)) { return 'NA' }
    try { return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant() }
    catch { return 'hash_error' }
}

function Get-PyDelPhiCommit([string]$Root) {
    if ([string]::IsNullOrWhiteSpace($Root) -or -not (Test-Path -LiteralPath $Root)) { return 'NA' }
    try {
        $commit = & git -C $Root rev-parse HEAD 2>$null
        if ($LASTEXITCODE -eq 0 -and -not [string]::IsNullOrWhiteSpace($commit)) {
            return ([string]$commit).Trim()
        }
    } catch {}
    return 'unknown'
}

function Get-RowValue([object]$Row, [string]$Name, [string]$Default = 'NA') {
    if ($null -ne $Row -and $Row.PSObject.Properties.Name -contains $Name) {
        $value = [string]$Row.$Name
        if (-not [string]::IsNullOrWhiteSpace($value)) { return $value }
    }
    return $Default
}

function Add-OracleManifestRow([hashtable]$Values) {
    if ([string]::IsNullOrWhiteSpace($ManifestPath)) { return }
    $columns = @(
        'timestamp_utc','producer_script','command_class','fixture',
        'oracle_engine','python_exe','pyroot','pyroot_git_commit',
        'workroot','input_hash','output_hash','output_path',
        'thread_count','status','note'
    )
    $ordered = [ordered]@{}
    foreach ($column in $columns) {
        $ordered[$column] = if ($Values.ContainsKey($column)) { [string]$Values[$column] } else { 'NA' }
    }
    $row = [PSCustomObject]$ordered
    $dir = Split-Path -Parent $ManifestPath
    if (-not [string]::IsNullOrWhiteSpace($dir)) {
        New-Item -ItemType Directory -Path $dir -Force | Out-Null
    }
    if (Test-Path -LiteralPath $ManifestPath) {
        $row | Export-Csv -LiteralPath $ManifestPath -NoTypeInformation -Encoding ASCII -Append
    } else {
        $row | Export-Csv -LiteralPath $ManifestPath -NoTypeInformation -Encoding ASCII
    }
}

function Write-ReportOnly([string]$Status, [string]$Detail) {
    New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
    $report = Join-Path $WorkRoot 'pbe_delphi_export_report.csv'
    $rows = foreach ($fixture in @('barnase_linear_salt0', 'barnase_linear_salt015', 'barnase_linear_salt2', 'barnase_nwt', 'barnase_nwt_salt2')) {
        [PSCustomObject]@{
            fixture = $fixture
            salt_M = if ($fixture -eq 'barnase_linear_salt0') { '0' } elseif ($fixture -in @('barnase_linear_salt015', 'barnase_nwt')) { '0.15' } else { '2.0' }
            status = $Status
            detail = $Detail
            charge_grid = 'missing'
            surface_map = 'missing'
            midpoint_epsilon = 'missing'
            boundary_map = 'missing'
            operator_stage = 'missing'
            boundary_samples = 'missing'
            field_samples = 'missing'
            field_source = 'missing'
            field_precision_status = 'missing'
            nwt_trace = 'missing'
            nwt_solver_summary = 'missing'
            iscm_energy_oracle = 'missing'
            gaussian_oracle = 'missing'
            rpbe_oracle = 'missing'
            focus_probe_oracle = 'missing'
            stage_oracle = 'missing'
            surface_semantics = 'NA'
        }
    }
    $rows | Export-Csv -LiteralPath $report -NoTypeInformation -Encoding ASCII
    Add-OracleManifestRow @{
        timestamp_utc = (Get-Date).ToUniversalTime().ToString('o')
        producer_script = 'test/pbe_delphi_export_maps.ps1'
        command_class = 'export_report_only'
        fixture = 'NA'
        oracle_engine = 'pyDelPhi'
        python_exe = $PythonExe
        pyroot = $PyDelPhiRoot
        pyroot_git_commit = Get-PyDelPhiCommit $PyDelPhiRoot
        workroot = $WorkRoot
        input_hash = Get-HashText $GoldenCsv
        output_hash = Get-HashText $report
        output_path = $report
        thread_count = 'NA'
        status = $Status
        note = $Detail
    }
    if ($Strict) { throw "DelPhi map export strict failed: $Status - $Detail" }
    Write-Warning "DelPhi map export is report-only: $Status - $Detail"
    Write-Host "PBE DelPhi map export report: $report"
}

if (-not (Test-Path -LiteralPath $GoldenCsv)) {
    throw "Missing DelPhi golden CSV: $GoldenCsv"
}
if (-not (Test-Path -LiteralPath $PyDelPhiRoot)) {
    Write-ReportOnly -Status 'pyDelPhi_missing' -Detail $PyDelPhiRoot
    return
}

New-Item -ItemType Directory -Path $WorkRoot -Force | Out-Null
function Add-OutputsToClean([System.Collections.Generic.HashSet[string]]$Set, [string[]]$Names) {
    foreach ($name in $Names) { [void]$Set.Add($name) }
}

$outputsToClean = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)
Add-OutputsToClean $outputsToClean @('pbe_delphi_export_report.csv')
if ($FullFieldOracle) {
    Add-OutputsToClean $outputsToClean @(
        'pbe_oracle_field_samples.csv',
        'pbe_oracle_field_samples_fullprecision.csv'
    )
    if ($FullFieldMap) {
        Add-OutputsToClean $outputsToClean @('pbe_oracle_field_map_fullprecision.csv')
    }
}
if ($FullBoundaryOracle) {
    Add-OutputsToClean $outputsToClean @(
        'pbe_oracle_boundary_samples.csv',
        'pbe_oracle_boundary_map.csv',
        'pbe_oracle_boundary_inputs.csv',
        'pbe_oracle_boundary_formula_samples.csv'
    )
}
if ($FullOperatorOracle -or $FullChargeGrid) {
    Add-OutputsToClean $outputsToClean @(
        'pbe_oracle_operator_stage.csv',
        'pbe_oracle_sor_trace.csv',
        'pbe_oracle_solver_summary.csv',
        'pbe_oracle_charge_grid.csv'
    )
}
if ($FullChargeGrid) {
    Add-OutputsToClean $outputsToClean @('pbe_oracle_charge_grid_full.csv')
}
if ($NwtOnly -or $FullOperatorOracle -or $FullFieldOracle) {
    Add-OutputsToClean $outputsToClean @(
        'pbe_oracle_nwt_trace.csv',
        'pbe_oracle_nwt_solver_summary.csv'
    )
}
if ($FullIscmEnergyOracle) {
    Add-OutputsToClean $outputsToClean @(
        'pbe_oracle_iscm_transition.csv',
        'pbe_oracle_iscm_surface_charge_samples.csv',
        'pbe_oracle_iscm_full.csv',
        'pbe_oracle_iscm_energy_terms.csv',
        'pbe_oracle_energy_components.csv'
    )
}
if ($FullGaussianOracle) {
    Add-OutputsToClean $outputsToClean @(
        'pbe_oracle_gaussian_density_summary.csv',
        'pbe_oracle_gaussian_samples.csv',
        'pbe_oracle_gaussian_charge_grid.csv'
    )
}
if ($FullRpbOracle) {
    Add-OutputsToClean $outputsToClean @(
        'pbe_oracle_rpbe_source.csv',
        'pbe_oracle_rpbe_samples.csv',
        'pbe_oracle_rpbe_top_source.csv'
    )
}
if ($FullFocusProbeOracle) {
    Add-OutputsToClean $outputsToClean @(
        'pbe_oracle_focus_boundary.csv',
        'pbe_oracle_focus_probe.csv'
    )
}
if ($FullMidpointEpsilon) {
    Add-OutputsToClean $outputsToClean @(
        'pbe_oracle_midpoint_epsilon.csv',
        'pbe_oracle_midpoint_epsilon_samples.csv'
    )
}
if ($FullVdwmsState) {
    Add-OutputsToClean $outputsToClean @(
        'pbe_oracle_midpoint_state.csv',
        'pbe_oracle_vdwms_audit.csv',
        'pbe_oracle_exposed_points.csv',
        'pbe_oracle_atom_accessibility.csv',
        'pbe_oracle_surface_map.csv'
    )
}
if ($FullBoundaryOracle -or $FullOperatorOracle -or $FullChargeGrid -or
    $FullFieldOracle -or $FullIscmEnergyOracle -or $FullMidpointEpsilon -or
    $FullVdwmsState -or $NwtOnly) {
    Add-OutputsToClean $outputsToClean @('pbe_oracle_atom_radius_audit.csv')
}

foreach ($name in $outputsToClean) {
    $path = Join-Path $WorkRoot $name
    if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path -Force }
}
$barnaseRows = @(Import-Csv -LiteralPath $GoldenCsv | Where-Object {
    $_.fixture -in @('barnase_linear_salt0', 'barnase_linear_salt015', 'barnase_linear_salt2')
})
if ($barnaseRows.Count -eq 0) {
    throw "No barnase linear rows found in $GoldenCsv"
}

$driver = Join-Path $WorkRoot 'pbe_delphi_export_maps_driver.py'
$driverCode = @'
import argparse
import csv
import json
import math
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")


LINEAR_FIXTURES = [
    ("barnase_linear_salt0", "0"),
    ("barnase_linear_salt015", "0.15"),
    ("barnase_linear_salt2", "2.0"),
]
NWT_FIXTURES = [
    ("barnase_nwt", "0.15"),
    ("barnase_nwt_salt2", "2.0"),
]
GAUSSIAN_FIXTURES = [
    ("5tif_gaussian_cutoff_salt0", "0"),
    ("5tif_gaussian_cutoff_salt015", "0.15"),
    ("5tif_gaussian_cutoff_salt2", "2.0"),
    ("5tif_gaussian_cutoff_surface_cutoff_salt015", "0.15"),
    ("5tif_gaussian_cutoff_density050_salt015", "0.15"),
    ("5tif_gaussian_salt0", "0"),
    ("5tif_gaussian_salt015", "0.15"),
    ("5tif_gaussian_salt2", "2.0"),
    ("5tif_gaussian_charge_salt0", "0"),
    ("5tif_gaussian_charge_salt015", "0.15"),
    ("5tif_gaussian_charge_salt2", "2.0"),
]
RPBE_FIXTURES = [
    ("5tif_rpbe_salt0", "0"),
    ("5tif_rpbe_salt015", "0.15"),
    ("5tif_rpbe_salt2", "2.0"),
    ("5tif_rpbe_charge_salt015", "0.15"),
    ("5tif_rpbe_surface_cutoff_salt015", "0.15"),
]
FOCUS_FIXTURES = [
    ("barnase_focus_probe", "0.15"),
    ("barnase_focus_probe_residue_only", "0.15"),
    ("barnase_focus_probe_site_p", "0.15"),
    ("barnase_focus_probe_site_q", "0.15"),
    ("barnase_focus_probe_site_all", "0.15"),
    ("barnase_focus_probe_nomap", "0.15"),
]
FIXTURES = LINEAR_FIXTURES + NWT_FIXTURES + GAUSSIAN_FIXTURES + RPBE_FIXTURES + FOCUS_FIXTURES
SURFACE_SEMANTICS = "pydelphi_vdw_surfmap_1_solvent_0_solute"
MIDPOINT_LAYOUT = "ijk_coords_with_pydelphi_flat_index_i_yz_j_z_k_axis_xyz"
VDWMS_CAPTURES = []
VDWMS_FIND_CAPTURES = []
VDWMS_PATCHED = False
BOUNDARY_CAPTURES = []
BOUNDARY_PATCHED = False
FIELD_CAPTURES = []
FIELD_PATCHED = False
ISCM_CAPTURES = []
ISCM_PATCHED = False
GAUSSIAN_CAPTURES = []
GAUSSIAN_PATCHED = False
RPB_CAPTURES = []
RPB_PATCHED = False
ENERGY_COMPONENT_CAPTURES = []
OPERATOR_CAPTURES = []
SOR_TRACE_CAPTURES = []
NWT_TRACE_CAPTURES = []
NWT_SOLVER_SUMMARIES = []
OPERATOR_PATCHED = False
SOR_PATCHED = False
NWT_PATCHED = False
CURRENT_OPERATOR_SAMPLE_POINTS = set()


def read_csv(path):
    if not path.exists():
        return []
    with path.open(newline="", encoding="ascii", errors="replace") as f:
        return list(csv.DictReader(f))


def read_tsv(path):
    if not path.exists():
        return []
    with path.open(newline="", encoding="ascii", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "NA") for k in fields})


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "NA") for k in fields})


def install_vdwms_capture(enabled):
    global VDWMS_PATCHED
    if not enabled or VDWMS_PATCHED:
        return
    import pydelphi.space.vdwms as vdwms_mod

    original_create = vdwms_mod.SurfaceMolecularVdW.create_vdw_molecular_surfaces
    original_find = vdwms_mod._find_boundary_grid_points

    def wrapped_find(*args, **kwargs):
        result = original_find(*args, **kwargs)
        try:
            VDWMS_FIND_CAPTURES.append({
                "boundary_initial": int(result[0]),
                "external_boundary_initial": int(result[1]),
            })
        except Exception:
            VDWMS_FIND_CAPTURES.append({
                "boundary_initial": "NA",
                "external_boundary_initial": "NA",
            })
        return result

    def wrapped_create(self, *args, **kwargs):
        capture = {}
        capture["epsilon_dimension"] = getattr(self, "epsilon_dimension", "NA")
        capture["num_atoms"] = getattr(self, "num_atoms", "NA")
        capture["num_objects"] = getattr(self, "num_objects", "NA")
        raw = getattr(self, "discrete_epsilon_index_map_1d", None)
        if raw is not None:
            capture["raw_index"] = np.array(raw, copy=True)
        find_start = len(VDWMS_FIND_CAPTURES)
        result = original_create(self, *args, **kwargs)
        final = getattr(self, "discrete_epsilon_index_map_1d", None)
        if final is not None:
            capture["final_index"] = np.array(final, copy=True)
        find_rows = VDWMS_FIND_CAPTURES[find_start:]
        if find_rows:
            capture.update(find_rows[0])
        capture["boundary_final"] = getattr(self, "num_boundary_grid_points", "NA")
        capture["external_boundary_final"] = getattr(self, "num_external_boundary_points", "NA")
        capture["exposed_points"] = getattr(self, "num_exposed_grids", "NA")
        exposed = getattr(self, "exposed_grids_coords", None)
        if exposed is not None:
            capture["exposed_grids_coords"] = np.array(exposed, copy=True)
        atom_accessibility = getattr(self, "atom_accessibility", None)
        if atom_accessibility is not None:
            capture["atom_accessibility"] = np.array(atom_accessibility, copy=True)
        boundary = getattr(self, "boundary_grid_points", None)
        if boundary is not None:
            capture["boundary_grid_points"] = np.array(boundary, copy=True)
        VDWMS_CAPTURES.append(capture)
        return result

    vdwms_mod._find_boundary_grid_points = wrapped_find
    vdwms_mod.SurfaceMolecularVdW.create_vdw_molecular_surfaces = wrapped_create
    VDWMS_PATCHED = True


def install_boundary_capture(enabled):
    global BOUNDARY_PATCHED
    if not enabled or BOUNDARY_PATCHED:
        return
    import pydelphi.solver.pb.common_pb as common_pb
    import pydelphi.solver.pb.sor.nonlinear_pb as sor_pb
    try:
        import pydelphi.solver.pb.nwt.nonlinear_pb as nwt_pb
    except Exception:
        nwt_pb = None

    original_coulomb = sor_pb._cpu_setup_coulombic_boundary_condition
    original_dipolar = sor_pb._cpu_setup_dipolar_boundary_condition

    def grid_boundary_rows(grid_shape, phimap_1d):
        shape = tuple(int(x) for x in np.asarray(grid_shape).reshape(-1)[:3])
        nx, ny, nz = shape
        yz = ny * nz
        arr = np.asarray(phimap_1d).reshape(-1)
        rows = []
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if 0 < i < nx - 1 and 0 < j < ny - 1 and 0 < k < nz - 1:
                        continue
                    flat = i * yz + j * nz + k
                    rows.append((i, j, k, flat, float(arr[flat])))
        return shape, rows

    def atom_rows(atoms_data, limit=12):
        try:
            from pydelphi.constants import (
                ATOMFIELD_GRID_X, ATOMFIELD_GRID_Y, ATOMFIELD_GRID_Z,
                ATOMFIELD_CHARGE,
            )
            arr = np.asarray(atoms_data)
            rows = []
            for idx in range(min(arr.shape[0], limit)):
                rows.append({
                    "atom_index": idx + 1,
                    "atom_grid_x": float(arr[idx, ATOMFIELD_GRID_X]),
                    "atom_grid_y": float(arr[idx, ATOMFIELD_GRID_Y]),
                    "atom_grid_z": float(arr[idx, ATOMFIELD_GRID_Z]),
                    "atom_charge": float(arr[idx, ATOMFIELD_CHARGE]),
                })
            return rows
        except Exception:
            return []

    def capture_common(kind, vacuum, grid_spacing, exdi_scaled, indi_scaled,
                       debye_length, non_zero_salt, grid_shape, atoms_data,
                       phimap_1d):
        shape, rows = grid_boundary_rows(grid_shape, phimap_1d)
        BOUNDARY_CAPTURES.append({
            "kind": kind,
            "vacuum": bool(vacuum),
            "grid_spacing": float(grid_spacing),
            "exdi_scaled": float(exdi_scaled),
            "indi_scaled": float(indi_scaled),
            "debye_length": float(debye_length),
            "non_zero_salt": bool(non_zero_salt),
            "shape": shape,
            "rows": rows,
            "atom_rows": atom_rows(atoms_data),
        })

    def wrapped_coulomb(vacuum, grid_spacing, exdi_scaled, indi_scaled,
                        debye_length, non_zero_salt, grid_shape, atoms_data,
                        phimap_1d):
        result = original_coulomb(vacuum, grid_spacing, exdi_scaled,
                                  indi_scaled, debye_length, non_zero_salt,
                                  grid_shape, atoms_data, phimap_1d)
        capture_common("coulombic", vacuum, grid_spacing, exdi_scaled,
                       indi_scaled, debye_length, non_zero_salt, grid_shape,
                       atoms_data, phimap_1d)
        return result

    def wrapped_dipolar(vacuum, has_pve_charges, has_nve_charges, grid_spacing,
                        exdi_scaled, indi_scaled, debye_length, non_zero_salt,
                        total_pve_charge, total_nve_charge, grid_shape,
                        grid_centroid_pve_charge, grid_centroid_nve_charge,
                        phimap_1d):
        result = original_dipolar(vacuum, has_pve_charges, has_nve_charges,
                                  grid_spacing, exdi_scaled, indi_scaled,
                                  debye_length, non_zero_salt,
                                  total_pve_charge, total_nve_charge,
                                  grid_shape, grid_centroid_pve_charge,
                                  grid_centroid_nve_charge, phimap_1d)
        capture_common("dipolar", vacuum, grid_spacing, exdi_scaled,
                       indi_scaled, debye_length, non_zero_salt, grid_shape,
                       np.zeros((0, 8)), phimap_1d)
        return result

    common_pb._cpu_setup_coulombic_boundary_condition = wrapped_coulomb
    common_pb._cpu_setup_dipolar_boundary_condition = wrapped_dipolar
    sor_pb._cpu_setup_coulombic_boundary_condition = wrapped_coulomb
    sor_pb._cpu_setup_dipolar_boundary_condition = wrapped_dipolar
    if nwt_pb is not None:
        nwt_pb._cpu_setup_coulombic_boundary_condition = wrapped_coulomb
        nwt_pb._cpu_setup_dipolar_boundary_condition = wrapped_dipolar
    BOUNDARY_PATCHED = True


def install_field_capture(enabled):
    global FIELD_PATCHED
    if not enabled or FIELD_PATCHED:
        return
    import pydelphi.app.delphi as delphi_mod

    original_process = delphi_mod.DelphiApp._process_phase_pbe

    def wrapped_process(self, vacuum, n_threads, space_obj):
        phase = "vacuum" if bool(vacuum) else "water"
        result = original_process(self, vacuum, n_threads, space_obj)
        try:
            phimap = self.ctx.phimap_in_vacuum if bool(vacuum) else self.ctx.phimap_in_water
            if phimap is not None:
                FIELD_CAPTURES.append({
                    "phase": phase,
                    "vacuum": bool(vacuum),
                    "grid_shape": tuple(int(x) for x in np.asarray(self.ctx.grid_shape).reshape(-1)[:3]),
                    "phimap": np.array(phimap, dtype=np.float64, copy=True),
                })
        except Exception as exc:
            FIELD_CAPTURES.append({
                "phase": phase,
                "vacuum": bool(vacuum),
                "grid_shape": "NA",
                "phimap": None,
                "error": f"{type(exc).__name__}: {exc}",
            })
        return result

    delphi_mod.DelphiApp._process_phase_pbe = wrapped_process
    FIELD_PATCHED = True


def install_iscm_capture(enabled):
    global ISCM_PATCHED
    if not enabled or ISCM_PATCHED:
        return
    import pydelphi.energy.reactionfield_iscm as iscm_mod
    try:
        import pydelphi.energy.calculator as calc_mod
    except Exception:
        calc_mod = None

    original_calc = iscm_mod.calc_induced_charge_rf_energy

    def wrapped_calc(platform, unique_charged_gridpoints,
                     dielectric_boundary_grids, dielectric_boundary_map_1d,
                     atoms_data, induced_surf_charge_positions,
                     phimap_in_solvent, grid_shape, scale, indi_scaled,
                     debye_factor, epkt, dump_arrays=False):
        energy, induced = original_calc(
            platform, unique_charged_gridpoints, dielectric_boundary_grids,
            dielectric_boundary_map_1d, atoms_data, induced_surf_charge_positions,
            phimap_in_solvent, grid_shape, scale, indi_scaled, debye_factor,
            epkt, dump_arrays
        )
        try:
            induced_arr = np.asarray(induced, dtype=np.float64).reshape(-1, 4)
            boundary_arr = np.asarray(dielectric_boundary_grids, dtype=np.float64)
            charged_arr = np.asarray(unique_charged_gridpoints, dtype=np.float64)
            positions_arr = np.asarray(induced_surf_charge_positions, dtype=np.float64)
            atoms_arr = np.asarray(atoms_data, dtype=np.float64)
            thread_sums = None
            try:
                if positions_arr.ndim == 2 and positions_arr.shape[0] == induced_arr.shape[0] and atoms_arr.ndim == 2:
                    thread_count = int(getattr(platform, "names", {}).get("cpu", {}).get("num_threads", 0))
                    thread_sums = np.zeros(max(1, thread_count), dtype=np.float64)
                    try:
                        import pydelphi.energy.reactionfield_iscm as _iscm_for_terms
                        _ = _iscm_for_terms._calculate_rf_energy_unified(
                            induced.reshape(-1), positions_arr, atoms_arr,
                            float(epkt), max(1, thread_count), thread_sums)
                    except Exception:
                        thread_sums = None
            except Exception:
                thread_sums = None
            ISCM_CAPTURES.append({
                "energy": float(energy),
                "boundary_count": int(boundary_arr.shape[0]) if boundary_arr.ndim >= 2 else int(boundary_arr.size),
                "charged_grid_count": int(charged_arr.shape[0]) if charged_arr.ndim >= 2 else int(charged_arr.size),
                "induced": np.array(induced_arr, dtype=np.float64, copy=True),
                "positions": np.array(positions_arr, dtype=np.float64, copy=True),
                "atoms": np.array(atoms_arr, dtype=np.float64, copy=True),
                "thread_sums": None if thread_sums is None else np.array(thread_sums, dtype=np.float64, copy=True),
                "grid_shape": tuple(int(x) for x in np.asarray(grid_shape).reshape(-1)[:3]),
                "scale": float(scale),
                "indi_scaled": float(indi_scaled),
                "debye_factor": float(debye_factor),
                "epkt": float(epkt),
            })
        except Exception as exc:
            ISCM_CAPTURES.append({
                "energy": float(energy),
                "boundary_count": "NA",
                "charged_grid_count": "NA",
                "induced": None,
                "error": f"{type(exc).__name__}: {exc}",
            })
        return energy, induced

    iscm_mod.calc_induced_charge_rf_energy = wrapped_calc
    if calc_mod is not None:
        calc_mod.calc_induced_charge_rf_energy = wrapped_calc
    ISCM_PATCHED = True


def install_gaussian_capture(enabled):
    global GAUSSIAN_PATCHED
    if not enabled or GAUSSIAN_PATCHED:
        return
    import pydelphi.space.space as space_mod

    original_update = space_mod.Space.update_runtime_context

    def wrapped_update(self, ctx):
        result = original_update(self, ctx)
        try:
            density = getattr(self, "gauss_density_map_1d", None)
            density_mid = getattr(self, "gauss_density_map_midpoints_1d", None)
            if density is not None or density_mid is not None:
                GAUSSIAN_CAPTURES.append({
                    "grid_shape": tuple(int(x) for x in np.asarray(getattr(self, "grid_shape", [])).reshape(-1)[:3]),
                    "grid_origin": np.array(getattr(self, "grid_origin", []), dtype=np.float64, copy=True),
                    "surface_method": str(getattr(getattr(self, "surface_method", None), "name", getattr(self, "surface_method", "NA"))),
                    "dielectric_model": str(getattr(getattr(self, "dielectric_model", None), "name", getattr(self, "dielectric_model", "NA"))),
                    "density": None if density is None else np.array(density, dtype=np.float64, copy=True),
                    "density_midpoints": None if density_mid is None else np.array(density_mid, dtype=np.float64, copy=True),
                    "ion_exclusion": None if getattr(self, "ion_exclusion_map_1d", None) is None else np.array(getattr(self, "ion_exclusion_map_1d"), copy=True),
                })
        except Exception as exc:
            GAUSSIAN_CAPTURES.append({
                "grid_shape": "NA",
                "density": None,
                "density_midpoints": None,
                "error": f"{type(exc).__name__}: {exc}",
            })
        return result

    space_mod.Space.update_runtime_context = wrapped_update
    GAUSSIAN_PATCHED = True


def install_rpb_capture(enabled):
    global RPB_PATCHED
    if not enabled or RPB_PATCHED:
        return
    import pydelphi.solver.rpb.sor.linear_rpb as rpb_mod

    original_set_gridpoint_charges = rpb_mod.RPBESolver._set_gridpoint_charges

    def wrapped_set_gridpoint_charges(self, *args, **kwargs):
        vacuum = kwargs.get("vacuum", False)
        if "vacuum" not in kwargs and args:
            try:
                vacuum = bool(args[0])
            except Exception:
                vacuum = False
        result = original_set_gridpoint_charges(self, *args, **kwargs)
        try:
            source = getattr(self, "eps_dot_coul_map_1d", None)
            RPB_CAPTURES.append({
                "vacuum": bool(vacuum),
                "phase": "vacuum" if bool(vacuum) else "water",
                "grid_shape": tuple(int(x) for x in np.asarray(getattr(self, "grid_shape", [])).reshape(-1)[:3]),
                "coulomb": None if getattr(self, "coulomb_map_1d", None) is None else np.array(getattr(self, "coulomb_map_1d"), dtype=np.float64, copy=True),
                "grad_coulomb": None if getattr(self, "grad_coulomb_map_1d", None) is None else np.array(getattr(self, "grad_coulomb_map_1d"), dtype=np.float64, copy=True),
                "grad_epsilon": None if getattr(self, "grad_epsmap_1d", None) is None else np.array(getattr(self, "grad_epsmap_1d"), dtype=np.float64, copy=True),
                "grad_epsilon_internal": None if getattr(self, "grad_epsin_map_1d", None) is None else np.array(getattr(self, "grad_epsin_map_1d"), dtype=np.float64, copy=True),
                "source": None if source is None else np.array(source, dtype=np.float64, copy=True),
            })
        except Exception as exc:
            RPB_CAPTURES.append({
                "vacuum": bool(vacuum),
                "phase": "vacuum" if bool(vacuum) else "water",
                "source": None,
                "error": f"{type(exc).__name__}: {exc}",
            })
        return result

    rpb_mod.RPBESolver._set_gridpoint_charges = wrapped_set_gridpoint_charges
    RPB_PATCHED = True


def capture_energy_components(ctx, energies, fixture, salt, iscm_captures=None):
    try:
        results = getattr(ctx, "energy_results", None)
        if results is None:
            results = energies
        def get_energy(phase, name):
            try:
                if hasattr(results, "get_energy"):
                    value = results.get_energy(phase, name)
                    return "NA" if value is None else float(value)
            except Exception:
                pass
            try:
                if isinstance(results, dict):
                    value = results.get(phase, {}).get(name)
                    return "NA" if value is None else float(value)
            except Exception:
                pass
            try:
                data = getattr(results, "energies", {})
                if phase in data and name in data[phase]:
                    return float(data[phase][name])
            except Exception:
                pass
            return "NA"
        water_rf = get_energy("water", "reactionfield_energy")
        if water_rf == "NA":
            captures = list(iscm_captures or [])
            if not captures:
                captures = list(ISCM_CAPTURES)
            for cap in reversed(captures):
                try:
                    water_rf = float(cap.get("energy"))
                    break
                except Exception:
                    pass
        vacuum_rf = get_energy("vacuum", "reactionfield_energy")
        total_rf = get_energy("total", "corrected_reaction_field_energy")
        if total_rf == "NA" and water_rf != "NA" and vacuum_rf == "NA":
            total_rf = water_rf
        nonpolar = get_energy("phase_independent", "nonpolar_energy")
        overall = get_energy("total", "overall_solvation_energy")
        if water_rf != "NA" and vacuum_rf != "NA":
            mode = "water_minus_vacuum"
        elif water_rf != "NA":
            mode = "water_only"
        else:
            mode = "missing"
        epkt = safe_attr(ctx, "epkt")
        temp = safe_attr(ctx, "temperature")
        try:
            temp_value = float(temp) if temp is not None else 297.15
        except Exception:
            temp_value = 297.15
        kt_kj = 0.00831446261815324 * temp_value
        def to_kj(value):
            return "NA" if value == "NA" else float(value) * kt_kj
        return {
            "fixture": fixture,
            "salt_M": salt,
            "water_rf_kT": water_rf,
            "water_rf_kJmol": to_kj(water_rf),
            "vacuum_rf_kT": vacuum_rf,
            "vacuum_rf_kJmol": to_kj(vacuum_rf),
            "total_reactionfield_kT": total_rf,
            "total_reactionfield_kJmol": to_kj(total_rf),
            "nonpolar_kT": nonpolar,
            "overall_solvation_kT": overall,
            "final_selected_kT": total_rf,
            "final_selected_kJmol": to_kj(total_rf),
            "finalization_mode": mode,
            "epkt": "NA" if epkt is None else float(epkt),
            "scale": "NA" if safe_attr(ctx, "scale") is None else float(safe_attr(ctx, "scale")),
            "indi_scaled": "NA" if safe_attr(ctx, "internal_dielectric_scaled") is None else float(safe_attr(ctx, "internal_dielectric_scaled")),
            "exdi_scaled": "NA" if safe_attr(ctx, "external_dielectric_scaled") is None else float(safe_attr(ctx, "external_dielectric_scaled")),
            "debye_length": "NA" if safe_attr(ctx, "debye_length") is None else float(safe_attr(ctx, "debye_length")),
            "temperature_K": temp_value,
        }
    except Exception as exc:
        return {
            "fixture": fixture,
            "salt_M": salt,
            "water_rf_kT": "NA",
            "water_rf_kJmol": "NA",
            "vacuum_rf_kT": "NA",
            "vacuum_rf_kJmol": "NA",
            "total_reactionfield_kT": "NA",
            "total_reactionfield_kJmol": "NA",
            "nonpolar_kT": "NA",
            "overall_solvation_kT": "NA",
            "final_selected_kT": "NA",
            "final_selected_kJmol": "NA",
            "finalization_mode": f"error:{type(exc).__name__}:{str(exc).replace(',', ';')}",
            "epkt": "NA",
            "scale": "NA",
            "indi_scaled": "NA",
            "exdi_scaled": "NA",
            "debye_length": "NA",
            "temperature_K": "NA",
        }


def install_operator_capture(enabled):
    global OPERATOR_PATCHED
    if not enabled or OPERATOR_PATCHED:
        return
    import pydelphi.solver.pb.sor.nonlinear_pb as sor_pb

    original_prepare = sor_pb.NLPBESolver._prepare_to_iterate

    def wrapped_prepare(self, vacuum, exdi, is_gaussian_diel_model,
                        grid_spacing, debye_length, non_zero_salt, four_pi,
                        epkt, grid_shape, ion_exclusion_map_1d,
                        epsilon_map_1d, epsmap_midpoints_1d, charge_map_1d,
                        eps_midpoint_neighs_sum_plus_salt_screening_1d,
                        boundary_gridpoints_1d):
        charge_before = np.asarray(charge_map_1d, dtype=float).copy()
        result = original_prepare(
            self, vacuum, exdi, is_gaussian_diel_model, grid_spacing,
            debye_length, non_zero_salt, four_pi, epkt, grid_shape,
            ion_exclusion_map_1d, epsilon_map_1d, epsmap_midpoints_1d,
            charge_map_1d, eps_midpoint_neighs_sum_plus_salt_screening_1d,
            boundary_gridpoints_1d)
        OPERATOR_CAPTURES.append({
            "vacuum": bool(vacuum),
            "exdi": float(exdi),
            "grid_spacing": float(grid_spacing),
            "debye_length": float(debye_length),
            "non_zero_salt": bool(non_zero_salt),
            "four_pi": float(four_pi),
            "epkt": float(epkt),
            "grid_shape": tuple(int(x) for x in np.asarray(grid_shape).reshape(-1)[:3]),
            "charge_before": charge_before,
            "charge_scaled": np.asarray(charge_map_1d, dtype=float).copy(),
            "denom": np.asarray(eps_midpoint_neighs_sum_plus_salt_screening_1d, dtype=float).copy(),
            "boundary_flags": np.asarray(boundary_gridpoints_1d).copy(),
            "ion_exclusion": np.asarray(ion_exclusion_map_1d, dtype=float).copy(),
            "eps_midpoints": np.asarray(epsmap_midpoints_1d, dtype=float).copy(),
        })
        return result

    sor_pb.NLPBESolver._prepare_to_iterate = wrapped_prepare
    OPERATOR_PATCHED = True


def install_sor_trace_capture(enabled):
    global SOR_PATCHED
    if not enabled or SOR_PATCHED:
        return
    import pydelphi.solver.shared.sor.base as base
    import pydelphi.solver.pb.sor.nonlinear_pb as sor_pb

    original_plain = base._cpu_iterate_SOR
    original_rms = base._cpu_iterate_SOR_odd_with_dphi_rmsd
    call_counter = {"value": 0}

    def sample_update_rows(even_odd, omega_sor, grid_shape,
                           read_half, write_half, eps_midpoints, denom_map,
                           boundary_flags, charge_map, phase):
        if len(SOR_TRACE_CAPTURES) > 256:
            return []
        nx, ny, nz = tuple(int(x) for x in np.asarray(grid_shape).reshape(-1)[:3])
        y_stride = nz
        x_stride = ny * y_stride
        y_stride_half = y_stride // 2
        x_stride_half = x_stride // 2
        if int(even_odd) == 0:
            z_lead, z_lag = 1, 0
            y_lead, y_lag = y_stride_half, y_stride_half + 1
            x_lead, x_lag = x_stride_half, x_stride_half + 1
        else:
            z_lead, z_lag = 0, 1
            y_lead, y_lag = y_stride_half + 1, y_stride_half
            x_lead, x_lag = x_stride_half + 1, x_stride_half
        y_stride_x_3 = 3 * y_stride
        x_stride_x_3 = 3 * x_stride
        y_stride_x_3_minus_1 = y_stride_x_3 - 1
        rows = []
        eps_arr = np.asarray(eps_midpoints).reshape(-1)
        denom_arr = np.asarray(denom_map).reshape(-1)
        charge_arr = np.asarray(charge_map).reshape(-1)
        flags_arr = np.asarray(boundary_flags).reshape(-1)
        read_arr = np.asarray(read_half).reshape(-1)
        write_arr = np.asarray(write_half).reshape(-1)
        for i, j, k in sorted(CURRENT_OPERATOR_SAMPLE_POINTS):
            if i <= 0 or j <= 0 or k <= 0 or i >= nx - 1 or j >= ny - 1 or k >= nz - 1:
                continue
            flat = (i * ny + j) * nz + k
            if flat % 2 != int(even_odd):
                continue
            if int(flags_arr[flat]) & 1:
                continue
            half = flat // 2
            base_idx = flat * 3
            eps_k_minus = eps_arr[base_idx - 1]
            eps_k_plus = eps_arr[base_idx + 2]
            eps_j_minus = eps_arr[base_idx - y_stride_x_3_minus_1]
            eps_j_plus = eps_arr[base_idx + 1]
            eps_i_minus = eps_arr[base_idx - x_stride_x_3]
            eps_i_plus = eps_arr[base_idx]
            numerator = (
                read_arr[half - z_lag] * eps_k_minus +
                read_arr[half + z_lead] * eps_k_plus +
                read_arr[half - y_lag] * eps_j_minus +
                read_arr[half + y_lead] * eps_j_plus +
                read_arr[half - x_lag] * eps_i_minus +
                read_arr[half + x_lead] * eps_i_plus
            )
            old_phi = float(write_arr[half])
            charge = float(charge_arr[flat])
            denom = float(denom_arr[flat])
            new_phi = (1.0 - float(omega_sor)) * old_phi + float(omega_sor) * (numerator + charge) / denom
            rows.append({
                "call": call_counter["value"],
                "phase": phase,
                "even_odd": int(even_odd),
                "i": i, "j": j, "k": k,
                "omega": float(omega_sor),
                "old_phi": old_phi,
                "neighbor_sum": float(numerator),
                "charge_density": charge,
                "denom": denom,
                "new_phi_expected": float(new_phi),
            })
        return rows

    def wrapped_plain(even_odd, omega_sor, approx_zero, grid_shape,
                      phi_map_current_half_1d, phi_map_next_half_1d,
                      epsilon_map_midpoints_1d,
                      epsilon_sum_neighbors_plus_salt_screening_1d,
                      is_boundary_gridpoint_1d, charge_map_1d):
        call_counter["value"] += 1
        rows = sample_update_rows(even_odd, omega_sor, grid_shape,
                                  phi_map_current_half_1d,
                                  phi_map_next_half_1d,
                                  epsilon_map_midpoints_1d,
                                  epsilon_sum_neighbors_plus_salt_screening_1d,
                                  is_boundary_gridpoint_1d, charge_map_1d,
                                  "plain")
        result = original_plain(even_odd, omega_sor, approx_zero, grid_shape,
                                phi_map_current_half_1d,
                                phi_map_next_half_1d,
                                epsilon_map_midpoints_1d,
                                epsilon_sum_neighbors_plus_salt_screening_1d,
                                is_boundary_gridpoint_1d, charge_map_1d)
        for row in rows:
            flat = (row["i"] * int(grid_shape[1]) + row["j"]) * int(grid_shape[2]) + row["k"]
            row["new_phi_actual"] = float(np.asarray(phi_map_next_half_1d).reshape(-1)[flat // 2])
            row["dphi"] = row["new_phi_actual"] - row["old_phi"]
            SOR_TRACE_CAPTURES.append(row)
        return result

    def wrapped_rms(num_cpu_threads, even_odd, omega_sor, approx_zero,
                    grid_shape, phi_map_current_half_1d,
                    phi_map_next_half_1d, epsilon_map_midpoints_1d,
                    epsilon_sum_neighbors_plus_salt_screening_1d,
                    is_boundary_gridpoint_1d, charge_map_1d):
        call_counter["value"] += 1
        rows = sample_update_rows(even_odd, omega_sor, grid_shape,
                                  phi_map_current_half_1d,
                                  phi_map_next_half_1d,
                                  epsilon_map_midpoints_1d,
                                  epsilon_sum_neighbors_plus_salt_screening_1d,
                                  is_boundary_gridpoint_1d, charge_map_1d,
                                  "rms")
        result = original_rms(num_cpu_threads, even_odd, omega_sor,
                              approx_zero, grid_shape,
                              phi_map_current_half_1d,
                              phi_map_next_half_1d,
                              epsilon_map_midpoints_1d,
                              epsilon_sum_neighbors_plus_salt_screening_1d,
                              is_boundary_gridpoint_1d, charge_map_1d)
        for row in rows:
            flat = (row["i"] * int(grid_shape[1]) + row["j"]) * int(grid_shape[2]) + row["k"]
            row["new_phi_actual"] = float(np.asarray(phi_map_next_half_1d).reshape(-1)[flat // 2])
            row["dphi"] = row["new_phi_actual"] - row["old_phi"]
            row["rms_sum_sq"] = float(result[0])
            row["max_dphi_block"] = float(result[1])
            SOR_TRACE_CAPTURES.append(row)
        return result

    base._cpu_iterate_SOR = wrapped_plain
    base._cpu_iterate_SOR_odd_with_dphi_rmsd = wrapped_rms
    sor_pb._cpu_iterate_SOR = wrapped_plain
    sor_pb._cpu_iterate_SOR_odd_with_dphi_rmsd = wrapped_rms
    SOR_PATCHED = True


def install_nwt_trace_capture(enabled):
    global NWT_PATCHED
    if not enabled or NWT_PATCHED:
        return
    import pydelphi.solver.pb.nwt.base as nwt_base
    import pydelphi.solver.pb.nwt.nonlinear_pb as nwt_pb
    from pydelphi.constants import BOX_BOUNDARY, BOX_HOMO_EPSILON, BOX_ION_ACCESSIBLE

    original_plain = nwt_base._cpu_iterate_nwt
    original_rms = nwt_base._cpu_iterate_nwt_with_dphi_rmsd
    original_cpu_solve = nwt_pb.NLNewtonPBESolver._cpu_solve_nonlinear_pb_nwt
    call_counter = {"value": 0}

    def safe_sinh_cosh(phi):
        clip = 50.0
        sinh_clip = math.sinh(clip)
        if phi > clip:
            sinh_phi = sinh_clip + sinh_clip * (phi - clip)
            return sinh_phi, math.sqrt(1.0 + sinh_phi * sinh_phi)
        if phi < -clip:
            sinh_phi = -sinh_clip + sinh_clip * (phi + clip)
            return sinh_phi, math.sqrt(1.0 + sinh_phi * sinh_phi)
        return math.sinh(phi), math.cosh(phi)

    def sample_update_rows(vacuum, even_odd, non_zero_salt, approx_zero,
                           omega_adaptive, grid_shape, read_half, write_half,
                           salt_penalty, eps_midpoints, eps_sum_only,
                           boundary_flags, charge_map, ion_exclusion,
                           phase):
        if len(NWT_TRACE_CAPTURES) > 2048:
            return []
        nx, ny, nz = tuple(int(x) for x in np.asarray(grid_shape).reshape(-1)[:3])
        y_stride = nz
        x_stride = ny * y_stride
        y_stride_half = y_stride // 2
        x_stride_half = x_stride // 2
        if int(even_odd) == 0:
            z_lead, z_lag = 0, 1
            y_lead, y_lag = y_stride_half, y_stride_half + 1
            x_lead, x_lag = x_stride_half, x_stride_half + 1
        else:
            z_lead, z_lag = 1, 0
            y_lead, y_lag = y_stride_half + 1, y_stride_half
            x_lead, x_lag = x_stride_half + 1, x_stride_half
        y_stride_x_3 = 3 * y_stride
        x_stride_x_3 = 3 * x_stride
        y_stride_x_3_minus_1 = y_stride_x_3 - 1
        eps_arr = np.asarray(eps_midpoints).reshape(-1)
        eps_sum_arr = np.asarray(eps_sum_only).reshape(-1)
        salt_arr = np.asarray(salt_penalty).reshape(-1)
        charge_arr = np.asarray(charge_map).reshape(-1)
        flags_arr = np.asarray(boundary_flags).reshape(-1)
        read_arr = np.asarray(read_half).reshape(-1)
        write_arr = np.asarray(write_half).reshape(-1)
        sample_points = set(CURRENT_OPERATOR_SAMPLE_POINTS)
        if charge_arr.size:
            nonzero_charge = np.flatnonzero(np.abs(charge_arr) > float(approx_zero))
            for flat_charge in nonzero_charge[:96]:
                ii = int(flat_charge // (ny * nz))
                rem = int(flat_charge - ii * ny * nz)
                jj = int(rem // nz)
                kk = int(rem - jj * nz)
                sample_points.add((ii, jj, kk))
        rows = []
        for i, j, k in sorted(sample_points):
            if i <= 0 or j <= 0 or k <= 0 or i >= nx - 1 or j >= ny - 1 or k >= nz - 1:
                continue
            flat = (i * ny + j) * nz + k
            if flat % 2 != int(even_odd):
                continue
            if int(flags_arr[flat]) & int(BOX_BOUNDARY):
                continue
            half = flat // 2
            base_idx = flat * 3
            eps_sum = float(eps_sum_arr[flat])
            if (int(flags_arr[flat]) & int(BOX_HOMO_EPSILON)) == int(BOX_HOMO_EPSILON):
                eps = float(eps_arr[base_idx])
                neighbor_sum = eps * (
                    read_arr[half - z_lag] + read_arr[half + z_lead] +
                    read_arr[half - y_lag] + read_arr[half + y_lead] +
                    read_arr[half - x_lag] + read_arr[half + x_lead]
                )
            else:
                eps_k_minus = eps_arr[base_idx - 1]
                eps_k_plus = eps_arr[base_idx + 2]
                eps_j_minus = eps_arr[base_idx - y_stride_x_3_minus_1]
                eps_j_plus = eps_arr[base_idx + 1]
                eps_i_minus = eps_arr[base_idx - x_stride_x_3]
                eps_i_plus = eps_arr[base_idx]
                neighbor_sum = (
                    read_arr[half - z_lag] * eps_k_minus +
                    read_arr[half + z_lead] * eps_k_plus +
                    read_arr[half - y_lag] * eps_j_minus +
                    read_arr[half + y_lead] * eps_j_plus +
                    read_arr[half - x_lag] * eps_i_minus +
                    read_arr[half + x_lead] * eps_i_plus
                )
            old_phi = float(write_arr[half])
            charge = float(charge_arr[flat])
            screen = 0.0
            sinh_phi = 0.0
            cosh_phi = 1.0
            linearized_rhs = 0.0
            denom = eps_sum
            if (not bool(vacuum)) and bool(non_zero_salt) and (
                (int(flags_arr[flat]) & int(BOX_ION_ACCESSIBLE)) == int(BOX_ION_ACCESSIBLE)
            ):
                sinh_phi, cosh_phi = safe_sinh_cosh(old_phi)
                screen = float(salt_arr[flat])
                linearized_rhs = screen * (old_phi * cosh_phi - sinh_phi)
                neighbor_sum = float(neighbor_sum) + linearized_rhs
                denom = eps_sum + screen * cosh_phi
            if abs(charge) > float(approx_zero):
                target = (float(neighbor_sum) + charge) / denom
            else:
                target = float(neighbor_sum) / denom
            new_phi = (1.0 - float(omega_adaptive)) * old_phi + float(omega_adaptive) * target
            rows.append({
                "call": call_counter["value"],
                "phase": phase,
                "vacuum": int(bool(vacuum)),
                "iteration": call_counter["value"],
                "color": int(even_odd),
                "i": i, "j": j, "k": k,
                "omega": float(omega_adaptive),
                "old_phi": old_phi,
                "neighbor_sum": float(neighbor_sum) - linearized_rhs,
                "charge_density": charge,
                "eps_sum": eps_sum,
                "screen": screen,
                "sinh_phi": float(sinh_phi),
                "cosh_phi": float(cosh_phi),
                "linearized_rhs": float(linearized_rhs),
                "denom": float(denom),
                "target_phi": float(target),
                "new_phi_expected": float(new_phi),
                "field_units": "pydelphi_raw_epkt",
                "operator_convention": "pydelphi_nwt",
            })
        return rows

    def wrapped_plain(vacuum, even_odd, non_zero_salt, approx_zero,
                      omega_adaptive, grid_shape, phi_map_current_half_1d,
                      phi_map_next_half_1d,
                      salt_ions_solvation_penalty_map_1d,
                      epsilon_map_midpoints_1d,
                      epsilon_sum_neighbors_only_1d,
                      boundary_flags_1d, charge_map_1d,
                      ion_exclusion_map_1d):
        call_counter["value"] += 1
        rows = sample_update_rows(
            vacuum, even_odd, non_zero_salt, approx_zero, omega_adaptive,
            grid_shape, phi_map_current_half_1d, phi_map_next_half_1d,
            salt_ions_solvation_penalty_map_1d, epsilon_map_midpoints_1d,
            epsilon_sum_neighbors_only_1d, boundary_flags_1d, charge_map_1d,
            ion_exclusion_map_1d, "plain")
        result = original_plain(
            vacuum, even_odd, non_zero_salt, approx_zero, omega_adaptive,
            grid_shape, phi_map_current_half_1d, phi_map_next_half_1d,
            salt_ions_solvation_penalty_map_1d, epsilon_map_midpoints_1d,
            epsilon_sum_neighbors_only_1d, boundary_flags_1d, charge_map_1d,
            ion_exclusion_map_1d)
        ny, nz = int(grid_shape[1]), int(grid_shape[2])
        write_arr = np.asarray(phi_map_next_half_1d).reshape(-1)
        for row in rows:
            flat = (int(row["i"]) * ny + int(row["j"])) * nz + int(row["k"])
            row["new_phi"] = float(write_arr[flat // 2])
            row["dphi"] = row["new_phi"] - row["old_phi"]
            row["rms_sum_sq"] = "NA"
            row["max_dphi_block"] = "NA"
            NWT_TRACE_CAPTURES.append(row)
        return result

    def wrapped_rms(vacuum, even_odd, non_zero_salt, approx_zero,
                    omega_adaptive, grid_shape, phi_map_current_half_1d,
                    phi_map_next_half_1d,
                    salt_ions_solvation_penalty_map_1d,
                    epsilon_map_midpoints_1d,
                    epsilon_sum_neighbors_only_1d,
                    boundary_flags_1d, charge_map_1d,
                    ion_exclusion_map_1d, num_cpu_threads):
        call_counter["value"] += 1
        rows = sample_update_rows(
            vacuum, even_odd, non_zero_salt, approx_zero, omega_adaptive,
            grid_shape, phi_map_current_half_1d, phi_map_next_half_1d,
            salt_ions_solvation_penalty_map_1d, epsilon_map_midpoints_1d,
            epsilon_sum_neighbors_only_1d, boundary_flags_1d, charge_map_1d,
            ion_exclusion_map_1d, "rms")
        result = original_rms(
            vacuum, even_odd, non_zero_salt, approx_zero, omega_adaptive,
            grid_shape, phi_map_current_half_1d, phi_map_next_half_1d,
            salt_ions_solvation_penalty_map_1d, epsilon_map_midpoints_1d,
            epsilon_sum_neighbors_only_1d, boundary_flags_1d, charge_map_1d,
            ion_exclusion_map_1d, num_cpu_threads)
        ny, nz = int(grid_shape[1]), int(grid_shape[2])
        write_arr = np.asarray(phi_map_next_half_1d).reshape(-1)
        for row in rows:
            flat = (int(row["i"]) * ny + int(row["j"])) * nz + int(row["k"])
            row["new_phi"] = float(write_arr[flat // 2])
            row["dphi"] = row["new_phi"] - row["old_phi"]
            row["rms_sum_sq"] = float(result[0])
            row["max_dphi_block"] = float(result[1])
            NWT_TRACE_CAPTURES.append(row)
        NWT_TRACE_CAPTURES.append({
            "call": call_counter["value"],
            "phase": "summary",
            "vacuum": int(bool(vacuum)),
            "iteration": call_counter["value"],
            "color": "summary",
            "i": 0, "j": 0, "k": 0,
            "omega": float(omega_adaptive),
            "old_phi": 0.0,
            "neighbor_sum": 0.0,
            "charge_density": 0.0,
            "eps_sum": 0.0,
            "screen": 0.0,
            "sinh_phi": 0.0,
            "cosh_phi": 1.0,
            "linearized_rhs": 0.0,
            "denom": 0.0,
            "target_phi": 0.0,
            "new_phi_expected": 0.0,
            "new_phi": 0.0,
            "dphi": 0.0,
            "rms_sum_sq": float(result[0]),
            "max_dphi_block": float(result[1]),
            "field_units": "pydelphi_raw_epkt",
            "operator_convention": "pydelphi_nwt",
        })
        return result

    def wrapped_cpu_solve(self, vacuum, phimap_current_1d, non_zero_salt,
                          approx_zero, omega_adaptive, grid_shape,
                          salt_ions_solvation_penalty_map_1d,
                          epsilon_map_midpoints_1d,
                          epsilon_sum_neighbors_sum_only_1d,
                          boundary_flags_1d, charge_map_1d,
                          ion_exclusion_map_1d, itr_block_size,
                          max_nonlinear_iters, rms_tol, dphi_tol,
                          check_dphi, num_cpu_threads, verbose=True):
        start = len(NWT_TRACE_CAPTURES)
        result = original_cpu_solve(
            self, vacuum, phimap_current_1d, non_zero_salt, approx_zero,
            omega_adaptive, grid_shape, salt_ions_solvation_penalty_map_1d,
            epsilon_map_midpoints_1d, epsilon_sum_neighbors_sum_only_1d,
            boundary_flags_1d, charge_map_1d, ion_exclusion_map_1d,
            itr_block_size, max_nonlinear_iters, rms_tol, dphi_tol,
            check_dphi, num_cpu_threads, verbose)
        NWT_SOLVER_SUMMARIES.append({
            "vacuum": int(bool(vacuum)),
            "omega": float(omega_adaptive),
            "block_size": int(itr_block_size),
            "max_nonlinear_iters": int(max_nonlinear_iters),
            "rms_tol": float(rms_tol),
            "dphi_tol": float(dphi_tol),
            "check_dphi": int(bool(check_dphi)),
            "num_cpu_threads": int(num_cpu_threads),
            "rmsd": float(result[0]) if len(result) > 0 else "NA",
            "max_dphi": float(result[1]) if len(result) > 1 else "NA",
            "iterations": int(result[2]) if len(result) > 2 else "NA",
            "stop_reason": (
                "dphi_tol"
                if len(result) > 1 and float(result[1]) <= float(dphi_tol)
                else (
                    "rms_tol"
                    if len(result) > 0 and float(result[0]) <= float(rms_tol)
                    else (
                        "work_budget"
                        if len(result) > 2 and int(result[2]) >= int(max_nonlinear_iters)
                        else "stagnation_or_iteration_control"
                    )
                )
            ),
            "trace_rows": len(NWT_TRACE_CAPTURES) - start,
            "field_units": "pydelphi_raw_epkt",
        })
        return result

    nwt_base._cpu_iterate_nwt = wrapped_plain
    nwt_base._cpu_iterate_nwt_with_dphi_rmsd = wrapped_rms
    nwt_pb._cpu_iterate_nwt = wrapped_plain
    nwt_pb._cpu_iterate_nwt_with_dphi_rmsd = wrapped_rms
    nwt_pb.NLNewtonPBESolver._cpu_solve_nonlinear_pb_nwt = wrapped_cpu_solve
    NWT_PATCHED = True


def flat_index(i, j, k, shape):
    nx, ny, nz = shape
    return (int(i) * ny + int(j)) * nz + int(k)


def sample_array(array, i, j, k, shape):
    if array is None:
        return "NA"
    arr = np.asarray(array)
    idx = flat_index(i, j, k, shape)
    if arr.ndim == 1:
        if idx >= arr.size:
            return "NA"
        return float(arr[idx])
    if arr.ndim == 3:
        return float(arr[int(i), int(j), int(k)])
    return "NA"


def surface_value_to_inside(value):
    if value == "NA" or value is None:
        return "NA"
    try:
        return 1 if float(value) < 0.5 else 0
    except Exception:
        return "NA"


def to_float_or_none(value):
    try:
        if value is None:
            return None
        text = str(value).strip()
        if not text or text.upper() == "NA":
            return None
        return float(text)
    except Exception:
        return None


def salt_matches(actual, wanted):
    if wanted is None or str(wanted).strip() == "":
        return True
    try:
        return abs(float(actual) - float(wanted)) <= 1.0e-12
    except Exception:
        return str(actual).strip() == str(wanted).strip()


def has_numeric_column(rows, name):
    for row in rows:
        if to_float_or_none(row.get(name)) is not None:
            return True
    return False


def timing_value(app, *needles):
    timings = getattr(app, "timings", None)
    if not timings:
        return "NA"
    if hasattr(timings, "items"):
        items = list(timings.items())
    else:
        items = []
        try:
            items = [(str(i), getattr(timings, str(i))) for i in dir(timings)]
        except Exception:
            return "NA"
    needles = [n.lower() for n in needles]
    total = 0.0
    found = False
    for key, value in items:
        k = str(key).lower()
        if not any(n in k for n in needles):
            continue
        parsed = to_float_or_none(value)
        if parsed is None:
            continue
        total += parsed
        found = True
    return total if found else "NA"


def sample_midpoints(array, i, j, k, shape):
    if array is None:
        return ("NA", "NA", "NA", "NA")
    if isinstance(array, (list, tuple)) and len(array) >= 3:
        vals = []
        for axis in range(3):
            v = sample_array(array[axis], i, j, k, shape)
            if v == "NA":
                return ("NA", "NA", "NA", "NA")
            vals.append(float(v))
        return (vals[0], vals[1], vals[2], sum(vals))
    arr = np.asarray(array)
    idx = flat_index(i, j, k, shape)
    if arr.ndim == 1 and arr.size >= (idx + 1) * 3:
        vals = [float(arr[idx * 3 + axis]) for axis in range(3)]
    elif arr.ndim == 2 and arr.shape[1] >= 3 and idx < arr.shape[0]:
        vals = [float(arr[idx, axis]) for axis in range(3)]
    elif arr.ndim == 4 and arr.shape[-1] >= 3:
        vals = [float(arr[int(i), int(j), int(k), axis]) for axis in range(3)]
    else:
        return ("NA", "NA", "NA", "NA")
    return (vals[0], vals[1], vals[2], sum(vals))


def sample_index_triplet(index_array, i, j, k, shape):
    if index_array is None:
        return ("NA", "NA", "NA")
    arr = np.asarray(index_array).reshape(-1)
    base = flat_index(i, j, k, shape) * 3
    if base + 2 >= arr.size:
        return ("NA", "NA", "NA")
    return (int(arr[base]), int(arr[base + 1]), int(arr[base + 2]))


def media_from_index(value, epsdim):
    if value == "NA":
        return "NA"
    try:
        return int(abs(int(value)) // int(epsdim))
    except Exception:
        return "NA"


def index_axis_counts(index_array, epsdim):
    if index_array is None:
        return {"x": "NA", "y": "NA", "z": "NA"}
    arr = np.asarray(index_array).reshape(-1)
    if arr.size < 3:
        return {"x": "NA", "y": "NA", "z": "NA"}
    return {
        "x": int(np.count_nonzero((np.abs(arr[0::3]) // epsdim) > 0)),
        "y": int(np.count_nonzero((np.abs(arr[1::3]) // epsdim) > 0)),
        "z": int(np.count_nonzero((np.abs(arr[2::3]) // epsdim) > 0)),
    }


def read_cube_data(path):
    path = Path(path)
    if not path.exists():
        return None
    with path.open("r", encoding="ascii", errors="replace") as f:
        first = f.readline().split()
        if len(first) < 2:
            return None
        nx = int(float(first[1]))
        f.readline()
        natom = int(float(f.readline().split()[0]))
        ny = int(float(f.readline().split()[0]))
        nz = int(float(f.readline().split()[0]))
        f.readline()
        for _ in range(abs(natom)):
            f.readline()
        data = np.fromstring(f.read(), sep=" ", dtype=float)
    expected = nx * ny * nz
    if data.size < expected:
        return None
    return data[:expected].reshape((nx, ny, nz))


def read_energy_summary(path):
    path = Path(path)
    if not path.exists():
        return {}
    lines = [line.strip() for line in path.read_text(encoding="ascii", errors="replace").splitlines() if line.strip()]
    data_lines = [line for line in lines if not line.startswith("#")]
    if len(data_lines) < 2:
        return {}
    header = data_lines[0].split("\t")
    values = data_lines[1].split("\t")
    return dict(zip(header, values))


def read_frc_total_energy(path):
    path = Path(path)
    if not path.exists():
        return "NA"
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if "Total energy" not in line:
            continue
        try:
            return str(float(line.split("=")[1].split()[0]))
        except Exception:
            return "NA"
    return "NA"


def safe_attr(obj, name):
    return getattr(obj, name, None)


def as_shape(ctx):
    shape = safe_attr(ctx, "grid_shape")
    if shape is None:
        return None
    shape = tuple(int(x) for x in list(shape)[:3])
    if len(shape) != 3:
        return None
    return shape


def as_origin_center_spacing(ctx):
    origin = safe_attr(ctx, "grid_origin")
    center = safe_attr(ctx, "grid_center")
    scale = safe_attr(ctx, "scale")
    try:
        h = 1.0 / float(scale)
    except Exception:
        h = "NA"
    def vec(v):
        if v is None:
            return ("NA", "NA", "NA")
        vv = list(np.asarray(v).reshape(-1))
        return tuple(float(x) for x in vv[:3])
    return vec(origin), vec(center), h


def stats(values):
    if values is None:
        return ("NA", "NA", "NA")
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 0:
        return ("NA", "NA", "NA")
    return (float(np.min(arr)), float(np.max(arr)), float(np.mean(arr)))


def epsilon_six_face_stats(eps_mid, shape):
    if not (isinstance(eps_mid, (list, tuple)) and len(eps_mid) >= 3 and all(x is not None for x in eps_mid[:3])):
        return ("NA", "NA", "NA")
    arrays = []
    expected = int(shape[0]) * int(shape[1]) * int(shape[2])
    for axis in range(3):
        arr = np.asarray(eps_mid[axis], dtype=float)
        if arr.size != expected:
            return ("NA", "NA", "NA")
        arrays.append(arr.reshape(shape))
    ex, ey, ez = arrays
    summed = ex + ey + ez
    summed[1:, :, :] += ex[:-1, :, :]
    summed[:, 1:, :] += ey[:, :-1, :]
    summed[:, :, 1:] += ez[:, :, :-1]
    return stats(summed)


def epsilon_axis_stats(eps_mid):
    if not (isinstance(eps_mid, (list, tuple)) and len(eps_mid) >= 3 and all(x is not None for x in eps_mid[:3])):
        return {}
    result = {}
    labels = ("x", "y", "z")
    for axis, label in enumerate(labels):
        arr = np.asarray(eps_mid[axis], dtype=float).reshape(-1)
        mn, mx, mean = stats(arr)
        result[f"eps_mid_{label}_min"] = mn
        result[f"eps_mid_{label}_max"] = mx
        result[f"eps_mid_{label}_mean"] = mean
        result[f"eps_mid_{label}_count"] = int(arr.size)
        result[f"eps_mid_{label}_internal_count"] = int(np.count_nonzero(arr <= 40.5))
    return result


def finite_sum(values):
    if values is None:
        return "NA"
    arr = np.asarray(values, dtype=float).reshape(-1)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return "NA"
    return float(np.sum(arr))


def finite_nonzero(values, threshold=0.0):
    if values is None:
        return "NA"
    arr = np.asarray(values, dtype=float).reshape(-1)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return "NA"
    if threshold > 0.0:
        return int(np.count_nonzero(np.abs(arr) > threshold))
    return int(np.count_nonzero(arr))


def finite_l2_norm(values):
    if values is None:
        return "NA"
    arr = np.asarray(values, dtype=float).reshape(-1)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return "NA"
    return float(np.sqrt(np.sum(arr * arr)))


def finite_rms(values):
    if values is None:
        return "NA"
    arr = np.asarray(values, dtype=float).reshape(-1)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return "NA"
    return float(np.sqrt(np.mean(arr * arr)))


def finite_max_abs(values):
    if values is None:
        return "NA"
    arr = np.asarray(values, dtype=float).reshape(-1)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return "NA"
    return float(np.max(np.abs(arr)))


def density_cutoff_count(values, cutoff):
    if values is None:
        return "NA"
    arr = np.asarray(values, dtype=float).reshape(-1)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return "NA"
    return int(np.count_nonzero(arr >= float(cutoff)))


def atom_radius_audit(ctx, fixture):
    atoms = safe_attr(ctx, "atoms_data")
    if atoms is None:
        return [], {}
    try:
        from pydelphi.constants import (
            ATOMFIELD_X,
            ATOMFIELD_Y,
            ATOMFIELD_Z,
            ATOMFIELD_CHARGE,
            ATOMFIELD_RADIUS,
            ATOMFIELD_MEDIA_ID,
        )
    except Exception:
        return [], {}
    arr = np.asarray(atoms, dtype=float)
    if arr.ndim != 2 or arr.shape[0] == 0:
        return [], {}
    radii = arr[:, ATOMFIELD_RADIUS]
    rows = []
    for idx in range(arr.shape[0]):
        rows.append({
            "fixture": fixture,
            "atom_index": idx + 1,
            "x_A": float(arr[idx, ATOMFIELD_X]),
            "y_A": float(arr[idx, ATOMFIELD_Y]),
            "z_A": float(arr[idx, ATOMFIELD_Z]),
            "q": float(arr[idx, ATOMFIELD_CHARGE]),
            "radius_A": float(arr[idx, ATOMFIELD_RADIUS]),
            "media_id": int(round(float(arr[idx, ATOMFIELD_MEDIA_ID]))),
        })
    stats_row = {
        "pydelphi_atom_count": int(arr.shape[0]),
        "pydelphi_radius_sum_A": float(np.sum(radii)),
        "pydelphi_radius_min_A": float(np.min(radii)),
        "pydelphi_radius_max_A": float(np.max(radii)),
        "pydelphi_radius_mean_A": float(np.mean(radii)),
    }
    return rows, stats_row


def build_prm(path, pyroot, outdir, salt, fixture, write_maps=True,
              phi_name="oracle_phimap.cube", eps_name="oracle_epsmap.cube",
              surf_name="oracle_surfmap.cube"):
    fixture_text = str(fixture).lower()
    phi = phi_name
    eps = eps_name
    surf = surf_name

    if fixture_text.startswith("5tif"):
        bdir = pyroot / "pydelphi" / "data" / "test_examples" / "5tif"
        shutil.copy2(bdir / "5tif.pqr", outdir / "5tif.pqr")
        biomodel = "rpbe" if "rpbe" in fixture_text else "pbe"
        surf_method = "gaussiancutoff" if "cutoff" in fixture_text else "gaussian"
        density_cutoff = "0.50" if "density050" in fixture_text else "0.759"
        lines = [
            f"biomodel={biomodel}",
            "solver=sor",
            f"surfacemethod={surf_method}",
            "dielectricmodel=gaussian",
            "gaussian_exponent=2",
            "gaussian_sigma=1.0",
            f"density_cutoff={density_cutoff}",
            "scale=2.0",
            "grid_size=139",
            "gridbox_type=cubic",
            'in(modpdb4,file="5tif.pqr",format="pqr")',
            "indi=1.0",
            "exdi=80.0",
            "gapdi=22.0",
            "boundary_condition=coulombic",
            f"salt={salt}",
            "linit=5000",
            "max_rmsd=1.0e-10",
            "max_delta_phi=1.0e-6",
        ]
        if write_maps:
            lines[11:11] = [
                f'out(phi,file="{phi}",media="both",format="cube")',
                f'out(eps,file="{eps}",media="both",point="both",format="cube")',
                f'out(surf,file="{surf}",format="cube")',
            ]
        path.write_text("\n".join(lines) + "\n", encoding="ascii")
        return

    bdir = pyroot / "pydelphi" / "data" / "test_examples" / "barnase"
    for src_name, dst_name in [
        ("barnase.pdb", "barnase.pdb"),
        ("amber.crg", "amber.crg"),
        ("amber.siz", "amber.siz"),
    ]:
        shutil.copy2(bdir / src_name, outdir / dst_name)
    if "focus" in fixture_text:
        shutil.copy2(bdir / "barstar-ASP39.pdb", outdir / "barstar-ASP39.pdb")
    pdb = "barnase.pdb"
    crg = "amber.crg"
    siz = "amber.siz"
    solver = "nwt" if "nwt" in fixture_text else "sor"
    boundary = "focusing" if "focus" in fixture_text else "coulombic"
    grid_size = "65" if "focus" in fixture_text else "137"
    scale = "4.0" if "focus" in fixture_text else "2.0"
    lines = [
        "biomodel=pbe",
        f"solver={solver}",
        "surfacemethod=vdw",
        "dielectricmodel=twodielectric",
        f"scale={scale}",
        f"grid_size={grid_size}",
        "gridbox_type=cubic",
        f'in(pdb,file="{pdb}")',
        f'in(crg,file="{crg}")',
        f'in(siz,file="{siz}")',
        "indi=1.0",
        "exdi=80.0",
        "gapdi=22.0",
        f"boundary_condition={boundary}",
        f"salt={salt}",
        "linit=5000",
        "nonlinit=0" if ("focus" in fixture_text or "linear" in fixture_text) else "nonlinit=5000",
        "max_rmsd=1.0e-10",
        "max_delta_phi=1.0e-6",
        "nwt_adaptive_omega=1.0",
    ]
    if "focus" in fixture_text:
        site_value = "a,q,p,f"
        if "site_p" in fixture_text:
            site_value = "p"
        elif "site_q" in fixture_text:
            site_value = "q"
        elif "residue_only" in fixture_text:
            site_value = ""
        lines.extend([
            'in(phi,file="parent_phimap_water.cube",format="cube")',
            'in(frc,file="barstar-ASP39.pdb")',
            'out(frc,file="oracle_frc.out")',
            "acent(29.3265,39.6727,8.45517)",
        ])
        if site_value:
            lines.append(f"site({site_value})")
    if write_maps:
        lines[10:10] = [
            f'out(phi,file="{phi}",media="both",format="cube")',
            f'out(eps,file="{eps}",media="both",point="both",format="cube")',
            f'out(surf,file="{surf}",format="cube")',
        ]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def run_pydelphi(pyroot, run_dir, fixture, salt, capture_vdwms=False,
                 capture_boundary=False, capture_operator=False,
                 capture_field=False, capture_iscm_energy=False,
                 capture_gaussian=False, capture_rpb=False, write_maps=True):
    sys.path.insert(0, str(pyroot))
    from pydelphi.app.delphi import DelphiApp
    from pydelphi.foundation.platforms import Platform
    from pydelphi.foundation.enums import Precision
    from pydelphi.config.global_runtime import set_precision
    install_vdwms_capture(capture_vdwms)
    install_boundary_capture(capture_boundary)

    platform = Platform(platform_name="cpu", debug=False)
    platform.activate("cpu", 2, None)
    platform.set_precision(Precision.DOUBLE)
    set_precision(platform.precision)
    install_field_capture(capture_field)
    install_iscm_capture(capture_field or capture_iscm_energy)
    install_operator_capture(capture_operator)
    install_sor_trace_capture(capture_operator)
    install_nwt_trace_capture(capture_operator)
    install_gaussian_capture(capture_gaussian)
    install_rpb_capture(capture_rpb)
    if capture_operator:
        try:
            import numba
            numba.set_num_threads(2)
        except Exception:
            pass
    fixture_text = str(fixture).lower()
    prm = run_dir / "barnase_oracle.prm"
    if "focus" in fixture_text:
        parent_prm = run_dir / "barnase_parent.prm"
        build_prm(parent_prm, pyroot, run_dir, salt, "barnase_linear_salt015",
                  write_maps=True,
                  phi_name="parent_phimap.cube",
                  eps_name="parent_epsmap.cube",
                  surf_name="parent_surfmap.cube")
        old_cwd_parent = Path.cwd()
        try:
            os.chdir(run_dir)
            parent_app = DelphiApp("barnase_parent.prm", platform, user_inputs=None)
            parent_app.run("pydelphi_parent_energy.csv", f"{fixture}_parent", True)
        finally:
            os.chdir(old_cwd_parent)
    build_prm(prm, pyroot, run_dir, salt, fixture, write_maps=write_maps)
    old_cwd = Path.cwd()
    try:
        os.chdir(run_dir)
        wall0 = time.perf_counter()
        app = DelphiApp("barnase_oracle.prm", platform, user_inputs=None)
        wall1 = time.perf_counter()
        energies = app.run("pydelphi_energy.csv", fixture, True)
        wall2 = time.perf_counter()
        app._titan_wall_setup_seconds = wall1 - wall0
        app._titan_wall_total_seconds = wall2 - wall0
        app._titan_wall_run_seconds = wall2 - wall1
    finally:
        os.chdir(old_cwd)
    return app, energies


def export_fixture(args, fixture, salt):
    global CURRENT_OPERATOR_SAMPLE_POINTS
    fixture_text = str(fixture).lower()
    pyroot = Path(args.pyroot)
    work = Path(args.workroot)
    run_dir = work / fixture
    run_dir.mkdir(parents=True, exist_ok=True)
    sample_root = Path(args.samples) / fixture if args.samples else None
    sample_field = read_csv(sample_root / "pbe_field_samples.csv") if sample_root else []
    sample_boundary = read_csv(sample_root / "pbe_boundary_samples.csv") if sample_root else []
    sample_charge = read_csv(sample_root / "pbe_charge_grid_audit.csv") if sample_root else []
    CURRENT_OPERATOR_SAMPLE_POINTS = set()
    for row in sample_field:
        try:
            CURRENT_OPERATOR_SAMPLE_POINTS.add((int(float(row["i"])), int(float(row["j"])), int(float(row["k"]))))
        except Exception:
            pass
    for row in sample_charge:
        try:
            CURRENT_OPERATOR_SAMPLE_POINTS.add((int(float(row["i"])), int(float(row["j"])), int(float(row["k"]))))
        except Exception:
            pass
    if getattr(args, "full_operator_oracle", False) and not CURRENT_OPERATOR_SAMPLE_POINTS:
        # NWT oracle fixtures may not have prior TITAN sample files. Keep a
        # deterministic interior stencil so trace export cannot pass empty.
        if "focus" in str(fixture).lower():
            default_points = [(32, 32, 32), (31, 32, 32), (32, 31, 32),
                              (32, 32, 31), (32, 32, 1), (32, 32, 63),
                              (1, 32, 32), (63, 32, 32)]
        else:
            default_points = [(68, 68, 68), (67, 68, 68), (68, 67, 68),
                              (68, 68, 67), (68, 68, 1), (68, 68, 135),
                              (1, 68, 68), (135, 68, 68)]
        for point in default_points:
            CURRENT_OPERATOR_SAMPLE_POINTS.add(point)

    capture_start = len(VDWMS_CAPTURES)
    boundary_start = len(BOUNDARY_CAPTURES)
    field_start = len(FIELD_CAPTURES)
    iscm_start = len(ISCM_CAPTURES)
    gaussian_start = len(GAUSSIAN_CAPTURES)
    rpb_start = len(RPB_CAPTURES)
    operator_start = len(OPERATOR_CAPTURES)
    if getattr(args, "full_operator_oracle", False):
        SOR_TRACE_CAPTURES.clear()
        NWT_TRACE_CAPTURES.clear()
        NWT_SOLVER_SUMMARIES.clear()
    sor_start = len(SOR_TRACE_CAPTURES)
    nwt_start = len(NWT_TRACE_CAPTURES)
    nwt_summary_start = len(NWT_SOLVER_SUMMARIES)
    app, energies = run_pydelphi(
        pyroot, run_dir, fixture, salt,
        capture_vdwms=getattr(args, "full_vdwms_state", False),
        capture_boundary=getattr(args, "full_boundary_oracle", False),
        capture_operator=getattr(args, "full_operator_oracle", False),
        capture_field=getattr(args, "full_field_oracle", False),
        capture_iscm_energy=getattr(args, "full_iscm_energy_oracle", False),
        capture_gaussian=getattr(args, "full_gaussian_oracle", False),
        capture_rpb=getattr(args, "full_rpb_oracle", False),
        write_maps=not getattr(args, "nwt_only", False),
    )
    boundary_captures = BOUNDARY_CAPTURES[boundary_start:]
    field_captures = FIELD_CAPTURES[field_start:]
    iscm_captures = ISCM_CAPTURES[iscm_start:]
    gaussian_captures = GAUSSIAN_CAPTURES[gaussian_start:]
    rpb_captures = RPB_CAPTURES[rpb_start:]
    operator_captures = OPERATOR_CAPTURES[operator_start:]
    sor_captures = SOR_TRACE_CAPTURES[sor_start:]
    nwt_captures = NWT_TRACE_CAPTURES[nwt_start:]
    nwt_summary_captures = NWT_SOLVER_SUMMARIES[nwt_summary_start:]
    vdwms_capture = VDWMS_CAPTURES[capture_start] if len(VDWMS_CAPTURES) > capture_start else {}
    ctx = app.ctx
    shape = as_shape(ctx)
    if shape is None:
        raise RuntimeError("pyDelPhi run completed but ctx.grid_shape is unavailable")
    origin, center, h = as_origin_center_spacing(ctx)

    rho = safe_attr(ctx, "charged_gridpoints_1d")
    surface = safe_attr(ctx, "dielectric_boundary_map_1d")
    eps_grid = safe_attr(ctx, "epsilon_map_1d")
    eps_mid = safe_attr(ctx, "epsilon_map_midpoints_water_1d")
    idx_mid = safe_attr(ctx, "index_discrete_epsilon_map_1d")
    if idx_mid is None:
        idx_mid = safe_attr(ctx, "discrete_epsilon_index_map_1d")
    if idx_mid is None and vdwms_capture:
        idx_mid = vdwms_capture.get("final_index")
    captured_water_phi = None
    captured_vacuum_phi = None
    for cap in field_captures:
        if cap.get("phimap") is None:
            continue
        if cap.get("phase") == "water":
            captured_water_phi = cap.get("phimap")
        elif cap.get("phase") == "vacuum":
            captured_vacuum_phi = cap.get("phimap")
    water_phi = captured_water_phi if captured_water_phi is not None else safe_attr(ctx, "phimap_in_water")
    vacuum_phi = captured_vacuum_phi if captured_vacuum_phi is not None else safe_attr(ctx, "phimap_in_vacuum")
    water_field_source = "solver_capture" if captured_water_phi is not None else ("ctx" if water_phi is not None else "cube_fallback")
    vacuum_field_source = "solver_capture" if captured_vacuum_phi is not None else ("ctx" if vacuum_phi is not None else "cube_fallback")
    field_precision_status = "full_precision" if water_field_source in ("solver_capture", "ctx") else "cube_rounded"
    if surface is None:
        surface = read_cube_data(run_dir / "oracle_surfmap.cube")
    if eps_grid is None:
        eps_grid = read_cube_data(run_dir / "oracle_epsmap_water.cube")
    if eps_mid is None:
        mids = [
            read_cube_data(run_dir / "oracle_epsmap_water_mids_4th_dim0.cube"),
            read_cube_data(run_dir / "oracle_epsmap_water_mids_4th_dim1.cube"),
            read_cube_data(run_dir / "oracle_epsmap_water_mids_4th_dim2.cube"),
        ]
        if all(x is not None for x in mids):
            eps_mid = mids
    if water_phi is None:
        water_phi = read_cube_data(run_dir / "oracle_phimap_water.cube")
        water_field_source = "cube_fallback"
        field_precision_status = "cube_rounded"
    if vacuum_phi is None:
        vacuum_phi = read_cube_data(run_dir / "oracle_phimap_vacuum.cube")
        vacuum_field_source = "cube_fallback"

    field_rows = []
    for row in sample_field:
        i, j, k = row["i"], row["j"], row["k"]
        water = sample_array(water_phi, i, j, k, shape)
        vac = sample_array(vacuum_phi, i, j, k, shape)
        reaction = "NA"
        if water != "NA" and vac != "NA":
            reaction = float(water) - float(vac)
        ex, ey, ez, es = sample_midpoints(eps_mid, i, j, k, shape)
        surf = sample_array(surface, i, j, k, shape)
        inside = surface_value_to_inside(surf)
        eps_sum6 = "NA" if es == "NA" else float(es) * 2.0
        field_rows.append({
            "fixture": fixture,
            "point_kind": row.get("point_kind", "sample"),
            "rank": row.get("rank", "NA"),
            "i": i, "j": j, "k": k,
            "rho": sample_array(rho, i, j, k, shape),
            "inside": inside,
            "surface": surf,
            "surface_semantics": SURFACE_SEMANTICS,
            "epsilon": sample_array(eps_grid, i, j, k, shape),
            "eps_x": ex, "eps_y": ey, "eps_z": ez, "eps_sum": es,
            "eps_sum6": eps_sum6,
            "solvent_phi": water,
            "water_phi": water,
            "vacuum_phi": vac,
            "reaction_phi": reaction,
            "field_source": water_field_source,
            "vacuum_field_source": vacuum_field_source,
            "field_precision_status": field_precision_status,
        })

    field_map_rows = []
    if getattr(args, "full_field_map", False) and salt_matches(salt, getattr(args, "full_midpoint_salt", "0.15")):
        if water_phi is not None or vacuum_phi is not None:
            nx, ny, nz = shape
            rank = 0
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        water = sample_array(water_phi, i, j, k, shape)
                        vac = sample_array(vacuum_phi, i, j, k, shape)
                        reaction = "NA"
                        if water != "NA" and vac != "NA":
                            reaction = float(water) - float(vac)
                        rank += 1
                        field_map_rows.append({
                            "fixture": fixture,
                            "rank": rank,
                            "i": i, "j": j, "k": k,
                            "solvent_phi": water,
                            "water_phi": water,
                            "vacuum_phi": vac,
                            "reaction_phi": reaction,
                            "field_source": water_field_source,
                            "vacuum_field_source": vacuum_field_source,
                            "field_precision_status": field_precision_status,
                        })

    boundary_rows = []
    for row in sample_boundary:
        i, j, k = row["i"], row["j"], row["k"]
        water = sample_array(water_phi, i, j, k, shape)
        vac = sample_array(vacuum_phi, i, j, k, shape)
        ex, ey, ez, es = sample_midpoints(eps_mid, i, j, k, shape)
        surf = sample_array(surface, i, j, k, shape)
        inside = surface_value_to_inside(surf)
        eps_sum6 = "NA" if es == "NA" else float(es) * 2.0
        boundary_rows.append({
            "fixture": fixture,
            "boundary_kind": row.get("boundary_kind", "sample_boundary"),
            "rank": row.get("rank", "NA"),
            "i": i, "j": j, "k": k,
            "inside": inside,
            "surface": surf,
            "surface_semantics": SURFACE_SEMANTICS,
            "epsilon": sample_array(eps_grid, i, j, k, shape),
            "eps_x": ex, "eps_y": ey, "eps_z": ez, "eps_sum": es,
            "eps_sum6": eps_sum6,
            "solvent_phi": water,
            "water_phi": water,
            "vacuum_phi": vac,
        })

    charge_rows = []
    for row in sample_charge:
        i, j, k = row["i"], row["j"], row["k"]
        ex, ey, ez, es = sample_midpoints(eps_mid, i, j, k, shape)
        surf = sample_array(surface, i, j, k, shape)
        inside = surface_value_to_inside(surf)
        eps_sum6 = "NA" if es == "NA" else float(es) * 2.0
        charge_rows.append({
            "fixture": fixture,
            "charge_kind": row.get("charge_kind", "sample_rho"),
            "rank": row.get("rank", "NA"),
            "i": i, "j": j, "k": k,
            "rho": sample_array(rho, i, j, k, shape),
            "inside": inside,
            "surface": surf,
            "surface_semantics": SURFACE_SEMANTICS,
            "epsilon": sample_array(eps_grid, i, j, k, shape),
            "eps_x": ex, "eps_y": ey, "eps_z": ez, "eps_sum": es,
            "eps_sum6": eps_sum6,
        })

    iscm_transition_rows = []
    iscm_sample_rows = []
    iscm_full_rows = []
    iscm_energy_term_rows = []
    for capture_index, cap in enumerate(iscm_captures):
        induced = cap.get("induced")
        positions = cap.get("positions")
        atoms_arr = cap.get("atoms")
        induced_count = int(induced.shape[0]) if induced is not None else 0
        induced_sum = float(np.sum(induced[:, 3])) if induced_count > 0 else 0.0
        induced_l1 = float(np.sum(np.abs(induced[:, 3]))) if induced_count > 0 else 0.0
        induced_l2 = float(np.sqrt(np.sum(induced[:, 3] * induced[:, 3]))) if induced_count > 0 else 0.0
        induced_max = float(np.max(np.abs(induced[:, 3]))) if induced_count > 0 else 0.0
        thread_sums = cap.get("thread_sums")
        thread_sum_raw = "NA" if thread_sums is None else float(np.sum(thread_sums))
        thread_sum_kT = "NA"
        if thread_sum_raw != "NA":
            thread_sum_kT = thread_sum_raw * float(cap.get("epkt", 0.0)) * 0.5
        iscm_transition_rows.append({
            "fixture": fixture,
            "salt_M": salt,
            "capture_index": capture_index,
            "boundary_count": cap.get("boundary_count", "NA"),
            "charged_grid_count": cap.get("charged_grid_count", "NA"),
            "induced_count": induced_count,
            "induced_charge_sum": induced_sum,
            "induced_charge_l1": induced_l1,
            "induced_charge_l2": induced_l2,
            "induced_charge_max_abs": induced_max,
            "water_rf_energy": cap.get("energy", "NA"),
            "thread_sum_raw": thread_sum_raw,
            "thread_sum_kT": thread_sum_kT,
            "scale": cap.get("scale", "NA"),
            "indi_scaled": cap.get("indi_scaled", "NA"),
            "debye_factor": cap.get("debye_factor", "NA"),
            "epkt": cap.get("epkt", "NA"),
        })
        if induced_count > 0:
            write_full = getattr(args, "full_iscm_energy_oracle", False) and salt_matches(salt, getattr(args, "full_midpoint_salt", "0.15"))
            sample_indices = list(range(min(64, induced_count)))
            top = np.argsort(-np.abs(induced[:, 3]))[:512]
            sample_indices.extend(int(x) for x in top)
            seen = set()
            rank = 0
            atom_field_x = atom_field_y = atom_field_z = atom_field_charge = None
            try:
                from pydelphi.constants import ATOMFIELD_X, ATOMFIELD_Y, ATOMFIELD_Z, ATOMFIELD_CHARGE
                atom_field_x, atom_field_y, atom_field_z, atom_field_charge = ATOMFIELD_X, ATOMFIELD_Y, ATOMFIELD_Z, ATOMFIELD_CHARGE
            except Exception:
                pass
            def potential_sum_for(idx):
                if positions is None or atoms_arr is None or atom_field_x is None:
                    return "NA"
                if positions.ndim != 2 or atoms_arr.ndim != 2 or idx >= positions.shape[0]:
                    return "NA"
                sx, sy, sz = float(positions[idx, 0]), float(positions[idx, 1]), float(positions[idx, 2])
                total = 0.0
                for atom_idx in range(atoms_arr.shape[0]):
                    dx = sx - float(atoms_arr[atom_idx, atom_field_x])
                    dy = sy - float(atoms_arr[atom_idx, atom_field_y])
                    dz = sz - float(atoms_arr[atom_idx, atom_field_z])
                    r2 = dx * dx + dy * dy + dz * dz
                    if r2 > 1.0e-24:
                        total += float(atoms_arr[atom_idx, atom_field_charge]) / math.sqrt(r2)
                return total
            top_energy_pairs = []
            if positions is not None and atoms_arr is not None and atom_field_x is not None:
                for idx in range(induced_count):
                    psum = potential_sum_for(idx)
                    if psum == "NA":
                        continue
                    energy_kt = 0.5 * float(cap.get("epkt", 0.0)) * float(induced[idx, 3]) * psum
                    top_energy_pairs.append((abs(energy_kt), idx, psum, energy_kt))
                top_energy_pairs.sort(reverse=True, key=lambda x: x[0])
                sample_indices.extend(int(x[1]) for x in top_energy_pairs[:512])
            for idx in sample_indices:
                if idx in seen or idx < 0 or idx >= induced_count:
                    continue
                seen.add(idx)
                rank += 1
                sx = sy = sz = "NA"
                if positions is not None and positions.ndim == 2 and idx < positions.shape[0]:
                    sx, sy, sz = float(positions[idx, 0]), float(positions[idx, 1]), float(positions[idx, 2])
                psum = potential_sum_for(idx)
                energy_kt = "NA" if psum == "NA" else 0.5 * float(cap.get("epkt", 0.0)) * float(induced[idx, 3]) * float(psum)
                iscm_sample_rows.append({
                    "fixture": fixture,
                    "salt_M": salt,
                    "capture_index": capture_index,
                    "rank": rank,
                    "i": int(induced[idx, 0]),
                    "j": int(induced[idx, 1]),
                    "k": int(induced[idx, 2]),
                    "surface_x": sx,
                    "surface_y": sy,
                    "surface_z": sz,
                    "induced_charge": float(induced[idx, 3]),
                    "atom_potential_sum": psum,
                    "energy_contribution_kT": energy_kt,
                })
            if write_full:
                for idx in range(induced_count):
                    sx = sy = sz = "NA"
                    if positions is not None and positions.ndim == 2 and idx < positions.shape[0]:
                        sx, sy, sz = float(positions[idx, 0]), float(positions[idx, 1]), float(positions[idx, 2])
                    psum = potential_sum_for(idx)
                    energy_kt = "NA" if psum == "NA" else 0.5 * float(cap.get("epkt", 0.0)) * float(induced[idx, 3]) * float(psum)
                    iscm_full_rows.append({
                        "fixture": fixture,
                        "salt_M": salt,
                        "capture_index": capture_index,
                        "rank": idx + 1,
                        "i": int(induced[idx, 0]),
                        "j": int(induced[idx, 1]),
                        "k": int(induced[idx, 2]),
                        "surface_x": sx,
                        "surface_y": sy,
                        "surface_z": sz,
                        "induced_charge": float(induced[idx, 3]),
                        "atom_potential_sum": psum,
                        "energy_contribution_kT": energy_kt,
                    })
                for term_rank, (_, idx, psum, energy_kt) in enumerate(top_energy_pairs[:512], start=1):
                    iscm_energy_term_rows.append({
                        "fixture": fixture,
                        "salt_M": salt,
                        "capture_index": capture_index,
                        "rank": term_rank,
                        "i": int(induced[idx, 0]),
                        "j": int(induced[idx, 1]),
                        "k": int(induced[idx, 2]),
                        "induced_charge": float(induced[idx, 3]),
                        "atom_potential_sum": psum,
                        "energy_contribution_kT": energy_kt,
                    })

    energy_component_rows = [capture_energy_components(ctx, energies, fixture, salt, iscm_captures)]

    boundary_map_rows = []
    boundary_inputs_rows = []
    boundary_formula_rows = []
    epkt = float(getattr(ctx, "epkt", float("nan")))
    if getattr(args, "full_boundary_oracle", False):
        for capture_index, cap in enumerate(boundary_captures):
            vacuum = bool(cap.get("vacuum", False))
            boundary_inputs_rows.append({
                "fixture": fixture,
                "salt_M": salt,
                "capture_index": capture_index,
                "kind": cap.get("kind", "coulombic"),
                "vacuum": int(vacuum),
                "grid_spacing": cap.get("grid_spacing", "NA"),
                "exdi_scaled": cap.get("exdi_scaled", "NA"),
                "indi_scaled": cap.get("indi_scaled", "NA"),
                "debye_length": cap.get("debye_length", "NA"),
                "non_zero_salt": int(bool(cap.get("non_zero_salt", False))),
                "epkt": epkt,
                "boundary_points": len(cap.get("rows", [])),
            })
            atom_lookup = cap.get("atom_rows", [])
            atom0 = atom_lookup[0] if atom_lookup else {}
            for i, j, k, flat, phi in cap.get("rows", []):
                boundary_map_rows.append({
                    "fixture": fixture,
                    "salt_M": salt,
                    "capture_index": capture_index,
                    "kind": cap.get("kind", "coulombic"),
                    "vacuum": int(vacuum),
                    "i": i, "j": j, "k": k,
                    "flat_id": flat,
                    "phi": phi,
                    "grid_spacing": cap.get("grid_spacing", "NA"),
                    "exdi_scaled": cap.get("exdi_scaled", "NA"),
                    "indi_scaled": cap.get("indi_scaled", "NA"),
                    "debye_length": cap.get("debye_length", "NA"),
                    "non_zero_salt": int(bool(cap.get("non_zero_salt", False))),
                    "epkt": epkt,
                    "atom_grid_x": atom0.get("atom_grid_x", "NA"),
                    "atom_grid_y": atom0.get("atom_grid_y", "NA"),
                    "atom_grid_z": atom0.get("atom_grid_z", "NA"),
                })
            sample_boundary_keys = {(int(r["i"]), int(r["j"]), int(r["k"])) for r in sample_boundary}
            for i, j, k, flat, phi in cap.get("rows", []):
                if (int(i), int(j), int(k)) not in sample_boundary_keys:
                    continue
                boundary_formula_rows.append({
                    "fixture": fixture,
                    "salt_M": salt,
                    "capture_index": capture_index,
                    "kind": cap.get("kind", "coulombic"),
                    "vacuum": int(vacuum),
                    "i": i, "j": j, "k": k,
                    "flat_id": flat,
                    "phi": phi,
                    "grid_spacing": cap.get("grid_spacing", "NA"),
                    "exdi_scaled": cap.get("exdi_scaled", "NA"),
                    "indi_scaled": cap.get("indi_scaled", "NA"),
                    "debye_length": cap.get("debye_length", "NA"),
                    "non_zero_salt": int(bool(cap.get("non_zero_salt", False))),
                    "epkt": epkt,
                    "atom_grid_x": atom0.get("atom_grid_x", "NA"),
                    "atom_grid_y": atom0.get("atom_grid_y", "NA"),
                    "atom_grid_z": atom0.get("atom_grid_z", "NA"),
                })

    operator_rows = []
    sor_trace_rows = []
    nwt_trace_rows = []
    nwt_summary_rows = []
    solver_summary_rows = []
    charge_full_rows = []
    full_operator_for_fixture = (
        getattr(args, "full_operator_oracle", False)
        and salt_matches(salt, getattr(args, "full_midpoint_salt", "0.15"))
    )
    full_charge_for_fixture = (
        (getattr(args, "full_charge_grid", False) or getattr(args, "full_gaussian_oracle", False))
        and salt_matches(salt, getattr(args, "full_midpoint_salt", "0.15"))
    )
    if getattr(args, "full_operator_oracle", False):
        sample_keys = []
        if full_operator_for_fixture:
            nx, ny, nz = shape
            sample_keys = [(i, j, k, "full_operator", (i * ny + j) * nz + k + 1)
                           for i in range(nx) for j in range(ny) for k in range(nz)]
        else:
            for row in sample_field:
                try:
                    sample_keys.append((int(float(row["i"])), int(float(row["j"])), int(float(row["k"])),
                                        row.get("point_kind", "sample_operator"), row.get("rank", "NA")))
                except Exception:
                    pass
        focus_fixture = "focus" in str(fixture).lower()
        for capture_index, cap in enumerate(operator_captures):
            cap_shape = tuple(cap.get("grid_shape", shape))
            if focus_fixture and tuple(int(x) for x in cap_shape) != tuple(int(x) for x in shape):
                continue
            nx, ny, nz = cap_shape
            y_stride = nz
            x_stride = ny * y_stride
            y_stride_x_3 = 3 * y_stride
            x_stride_x_3 = 3 * x_stride
            y_stride_x_3_minus_1 = y_stride_x_3 - 1
            eps_mid_arr = np.asarray(cap.get("eps_midpoints", []), dtype=float).reshape(-1)
            denom_arr = np.asarray(cap.get("denom", []), dtype=float).reshape(-1)
            charge_before_arr = np.asarray(cap.get("charge_before", []), dtype=float).reshape(-1)
            charge_scaled_arr = np.asarray(cap.get("charge_scaled", []), dtype=float).reshape(-1)
            flags_arr = np.asarray(cap.get("boundary_flags", []))
            ion_arr = np.asarray(cap.get("ion_exclusion", []), dtype=float).reshape(-1)
            phi_arr = np.asarray(water_phi if not cap.get("vacuum", False) else vacuum_phi) if (water_phi is not None or vacuum_phi is not None) else None
            if phi_arr is not None:
                phi_arr = phi_arr.reshape(-1)
            for i, j, k, point_kind, rank in sample_keys:
                if i < 0 or j < 0 or k < 0 or i >= nx or j >= ny or k >= nz:
                    continue
                flat = (i * ny + j) * nz + k
                base_idx = 3 * flat
                def eps_at(offset):
                    idx = base_idx + offset
                    if idx < 0 or idx >= eps_mid_arr.size:
                        return "NA"
                    return float(eps_mid_arr[idx])
                if 0 < i < nx - 1 and 0 < j < ny - 1 and 0 < k < nz - 1:
                    eps_i_plus = eps_at(0)
                    eps_j_plus = eps_at(1)
                    eps_k_plus = eps_at(2)
                    eps_i_minus = eps_at(-x_stride_x_3)
                    eps_j_minus = eps_at(-y_stride_x_3_minus_1)
                    eps_k_minus = eps_at(-1)
                    eps_values = [eps_i_plus, eps_i_minus, eps_j_plus, eps_j_minus, eps_k_plus, eps_k_minus]
                    eps_sum_no_screen = sum(float(v) for v in eps_values if v != "NA")
                else:
                    eps_i_plus = eps_i_minus = eps_j_plus = eps_j_minus = eps_k_plus = eps_k_minus = "NA"
                    eps_sum_no_screen = "NA"
                denom = float(denom_arr[flat]) if flat < denom_arr.size else "NA"
                if denom != "NA" and eps_sum_no_screen != "NA":
                    screen = float(denom) - float(eps_sum_no_screen)
                else:
                    screen = "NA"
                operator_rows.append({
                    "fixture": fixture,
                    "salt_M": salt,
                    "capture_index": capture_index,
                    "vacuum": int(bool(cap.get("vacuum", False))),
                    "point_kind": point_kind,
                    "rank": rank,
                    "i": i, "j": j, "k": k,
                    "rho": float(charge_before_arr[flat]) if flat < charge_before_arr.size else "NA",
                    "charge_map_scaled": float(charge_scaled_arr[flat]) if flat < charge_scaled_arr.size else "NA",
                    "rhs": float(charge_scaled_arr[flat]) if flat < charge_scaled_arr.size else "NA",
                    "screen": screen,
                    "eps_i_plus": eps_i_plus,
                    "eps_i_minus": eps_i_minus,
                    "eps_j_plus": eps_j_plus,
                    "eps_j_minus": eps_j_minus,
                    "eps_k_plus": eps_k_plus,
                    "eps_k_minus": eps_k_minus,
                    "eps_sum_no_screen": eps_sum_no_screen,
                    "denom": denom,
                    "boundary_flag": int(flags_arr.reshape(-1)[flat]) if flat < flags_arr.size else "NA",
                    "ion_exclusion": float(ion_arr[flat]) if flat < ion_arr.size else "NA",
                    "phi_initial": 0.0 if (0 < i < nx - 1 and 0 < j < ny - 1 and 0 < k < nz - 1) else sample_array(water_phi, i, j, k, shape),
                    "phi_water": float(phi_arr[flat]) if phi_arr is not None and flat < phi_arr.size else "NA",
                    "grid_spacing": cap.get("grid_spacing", "NA"),
                    "debye_length": cap.get("debye_length", "NA"),
                    "epkt": cap.get("epkt", "NA"),
                    "operator_units": "pydelphi_raw_epkt",
                })
        for row in sor_captures:
            sor_trace_rows.append({
                "fixture": fixture,
                "salt_M": salt,
                "call": row.get("call", "NA"),
                "phase": row.get("phase", "NA"),
                "even_odd": row.get("even_odd", "NA"),
                "i": row.get("i", "NA"),
                "j": row.get("j", "NA"),
                "k": row.get("k", "NA"),
                "omega": row.get("omega", "NA"),
                "old_phi": row.get("old_phi", "NA"),
                "neighbor_sum": row.get("neighbor_sum", "NA"),
                "charge_density": row.get("charge_density", "NA"),
                "denom": row.get("denom", "NA"),
                "new_phi": row.get("new_phi_actual", row.get("new_phi_expected", "NA")),
                "new_phi_expected": row.get("new_phi_expected", "NA"),
                "dphi": row.get("dphi", "NA"),
                "rms_sum_sq": row.get("rms_sum_sq", "NA"),
                "max_dphi_block": row.get("max_dphi_block", "NA"),
            })
        for row in nwt_captures:
            nwt_trace_rows.append({
                "fixture": fixture,
                "salt_M": salt,
                "call": row.get("call", "NA"),
                "phase": row.get("phase", "NA"),
                "vacuum": row.get("vacuum", "NA"),
                "iteration": row.get("iteration", "NA"),
                "color": row.get("color", "NA"),
                "i": row.get("i", "NA"),
                "j": row.get("j", "NA"),
                "k": row.get("k", "NA"),
                "omega": row.get("omega", "NA"),
                "old_phi": row.get("old_phi", "NA"),
                "neighbor_sum": row.get("neighbor_sum", "NA"),
                "charge_density": row.get("charge_density", "NA"),
                "eps_sum": row.get("eps_sum", "NA"),
                "screen": row.get("screen", "NA"),
                "sinh_phi": row.get("sinh_phi", "NA"),
                "cosh_phi": row.get("cosh_phi", "NA"),
                "linearized_rhs": row.get("linearized_rhs", "NA"),
                "denom": row.get("denom", "NA"),
                "target_phi": row.get("target_phi", "NA"),
                "new_phi": row.get("new_phi", row.get("new_phi_expected", "NA")),
                "new_phi_expected": row.get("new_phi_expected", "NA"),
                "dphi": row.get("dphi", "NA"),
                "rms_sum_sq": row.get("rms_sum_sq", "NA"),
                "max_dphi_block": row.get("max_dphi_block", "NA"),
                "field_units": row.get("field_units", "pydelphi_raw_epkt"),
                "operator_convention": row.get("operator_convention", "pydelphi_nwt"),
            })
        for row in nwt_summary_captures:
            nwt_summary_rows.append({
                "fixture": fixture,
                "salt_M": salt,
                "vacuum": row.get("vacuum", "NA"),
                "temp_K": "297.15" if "nwt" in str(fixture).lower() else "298.15",
                "omega": row.get("omega", "NA"),
                "block_size": row.get("block_size", "NA"),
                "max_nonlinear_iters": row.get("max_nonlinear_iters", "NA"),
                "rms_tol": row.get("rms_tol", "NA"),
                "dphi_tol": row.get("dphi_tol", "NA"),
                "check_dphi": row.get("check_dphi", "NA"),
                "num_cpu_threads": row.get("num_cpu_threads", "NA"),
                "rmsd": row.get("rmsd", "NA"),
                "max_dphi": row.get("max_dphi", "NA"),
                "iterations": row.get("iterations", "NA"),
                "stop_reason": row.get("stop_reason", "NA"),
                "trace_rows": row.get("trace_rows", "NA"),
                "field_units": row.get("field_units", "pydelphi_raw_epkt"),
            })
        for capture_index, cap in enumerate(operator_captures):
            solver_summary_rows.append({
                "fixture": fixture,
                "salt_M": salt,
                "capture_index": capture_index,
                "vacuum": int(bool(cap.get("vacuum", False))),
                "grid_spacing": cap.get("grid_spacing", "NA"),
                "debye_length": cap.get("debye_length", "NA"),
                "non_zero_salt": int(bool(cap.get("non_zero_salt", False))),
                "epkt": cap.get("epkt", "NA"),
                "operator_units": "pydelphi_raw_epkt",
                "operator_rows": len(operator_rows),
                "sor_trace_rows": len(sor_trace_rows),
            })

    if full_charge_for_fixture:
        focus_fixture = "focus" in str(fixture).lower()
        for capture_index, cap in enumerate(operator_captures):
            if bool(cap.get("vacuum", False)):
                continue
            cap_shape = tuple(cap.get("grid_shape", shape))
            if focus_fixture and tuple(int(x) for x in cap_shape) != tuple(int(x) for x in shape):
                continue
            nx, ny, nz = cap_shape
            charge_before_arr = np.asarray(cap.get("charge_before", []), dtype=float).reshape(-1)
            charge_scaled_arr = np.asarray(cap.get("charge_scaled", []), dtype=float).reshape(-1)
            total_size = max(charge_before_arr.size, charge_scaled_arr.size)
            if total_size == 0:
                continue
            rank = 0
            for flat in range(total_size):
                rho_value = float(charge_before_arr[flat]) if flat < charge_before_arr.size else 0.0
                rhs_value = float(charge_scaled_arr[flat]) if flat < charge_scaled_arr.size else 0.0
                if rho_value == 0.0 and rhs_value == 0.0:
                    continue
                ii = int(flat // (ny * nz))
                rem = int(flat - ii * ny * nz)
                jj = int(rem // nz)
                kk = int(rem - jj * nz)
                ex, ey, ez, es = sample_midpoints(eps_mid, ii, jj, kk, shape)
                surf = sample_array(surface, ii, jj, kk, shape)
                inside = surface_value_to_inside(surf)
                eps_sum6 = "NA" if es == "NA" else float(es) * 2.0
                rank += 1
                charge_full_rows.append({
                    "fixture": fixture,
                    "salt_M": salt,
                    "capture_index": capture_index,
                    "charge_kind": "full_charge_grid",
                    "rank": rank,
                    "i": ii, "j": jj, "k": kk,
                    "rho": rho_value,
                    "charge_map_scaled": rhs_value,
                    "rhs": rhs_value,
                    "inside": inside,
                    "surface": surf,
                    "surface_semantics": SURFACE_SEMANTICS,
                    "epsilon": sample_array(eps_grid, ii, jj, kk, shape),
                    "eps_x": ex, "eps_y": ey, "eps_z": ez,
                    "eps_sum": es,
                    "eps_sum6": eps_sum6,
                })
            break

    density_center = safe_attr(ctx, "gauss_density_map_1d")
    density_mid = safe_attr(ctx, "gauss_density_map_midpoints_1d")
    gaussian_capture_source = "ctx"
    if density_center is None or density_mid is None:
        for cap in reversed(gaussian_captures):
            cap_density = cap.get("density")
            cap_density_mid = cap.get("density_midpoints")
            if density_center is None and cap_density is not None:
                density_center = cap_density
                gaussian_capture_source = "space_update_runtime_context"
            if density_mid is None and cap_density_mid is not None:
                density_mid = cap_density_mid
                gaussian_capture_source = "space_update_runtime_context"
            if density_center is not None and density_mid is not None:
                break
    gaussian_summary_rows = []
    gaussian_sample_rows = []
    gaussian_charge_rows = []
    if getattr(args, "full_gaussian_oracle", False) and (
        str(fixture).startswith("5tif_gaussian") or str(fixture).startswith("5tif_rpbe")
    ):
        dens_min, dens_max, dens_mean = stats(density_center)
        dens_mid_min, dens_mid_max, dens_mid_mean = stats(density_mid)
        eps_axis_gaussian = epsilon_axis_stats(eps_mid)
        charge_arr = np.asarray(rho, dtype=float).reshape(-1) if rho is not None else np.array([])
        density_cutoff_value = 0.50 if "density050" in str(fixture).lower() else 0.759
        gaussian_summary_rows.append({
            "fixture": fixture,
            "salt_M": salt,
            "density_source": f"{gaussian_capture_source}.gauss_density_map_1d" if density_center is not None else "missing",
            "midpoint_density_source": f"{gaussian_capture_source}.gauss_density_map_midpoints_1d" if density_mid is not None else "missing",
            "grid_nx": shape[0],
            "grid_ny": shape[1],
            "grid_nz": shape[2],
            "density_sum": finite_sum(density_center),
            "gaussian_density_sum": finite_sum(density_center),
            "density_min": dens_min,
            "density_max": dens_max,
            "density_mean": dens_mean,
            "midpoint_density_sum": finite_sum(density_mid),
            "midpoint_density_min": dens_mid_min,
            "midpoint_density_max": dens_mid_max,
            "midpoint_density_mean": dens_mid_mean,
            "density_cutoff": density_cutoff_value,
            "density_cutoff_points": density_cutoff_count(density_center, density_cutoff_value),
            "gaussian_cutoff_points": density_cutoff_count(density_center, density_cutoff_value),
            "midpoint_cutoff_points": density_cutoff_count(density_mid, density_cutoff_value),
            "eps_mid_x_mean": eps_axis_gaussian.get("eps_mid_x_mean", "NA"),
            "eps_mid_y_mean": eps_axis_gaussian.get("eps_mid_y_mean", "NA"),
            "eps_mid_z_mean": eps_axis_gaussian.get("eps_mid_z_mean", "NA"),
            "eps_mid_min": stats(np.concatenate([np.asarray(x, dtype=float).reshape(-1) for x in eps_mid[:3]]))[0]
                if isinstance(eps_mid, (list, tuple)) and len(eps_mid) >= 3 and all(x is not None for x in eps_mid[:3]) else "NA",
            "gaussian_eps_min": stats(np.concatenate([np.asarray(x, dtype=float).reshape(-1) for x in eps_mid[:3]]))[0]
                if isinstance(eps_mid, (list, tuple)) and len(eps_mid) >= 3 and all(x is not None for x in eps_mid[:3]) else "NA",
            "eps_mid_max": stats(np.concatenate([np.asarray(x, dtype=float).reshape(-1) for x in eps_mid[:3]]))[1]
                if isinstance(eps_mid, (list, tuple)) and len(eps_mid) >= 3 and all(x is not None for x in eps_mid[:3]) else "NA",
            "gaussian_eps_max": stats(np.concatenate([np.asarray(x, dtype=float).reshape(-1) for x in eps_mid[:3]]))[1]
                if isinstance(eps_mid, (list, tuple)) and len(eps_mid) >= 3 and all(x is not None for x in eps_mid[:3]) else "NA",
            "eps_mid_mean": stats(np.concatenate([np.asarray(x, dtype=float).reshape(-1) for x in eps_mid[:3]]))[2]
                if isinstance(eps_mid, (list, tuple)) and len(eps_mid) >= 3 and all(x is not None for x in eps_mid[:3]) else "NA",
            "gaussian_eps_mean": stats(np.concatenate([np.asarray(x, dtype=float).reshape(-1) for x in eps_mid[:3]]))[2]
                if isinstance(eps_mid, (list, tuple)) and len(eps_mid) >= 3 and all(x is not None for x in eps_mid[:3]) else "NA",
            "charge_nonzero_cells": finite_nonzero(charge_arr),
            "charge_sum": finite_sum(charge_arr),
            "charge_grid_source": "ctx.charged_gridpoints_1d" if rho is not None else "missing",
            "status": "ok" if density_center is not None and density_mid is not None else "missing_density",
        })
        sample_keys = []
        if sample_field:
            for row in sample_field:
                try:
                    sample_keys.append((int(float(row["i"])), int(float(row["j"])), int(float(row["k"])), row.get("point_kind", "sample"), row.get("rank", "NA")))
                except Exception:
                    pass
        if not sample_keys:
            nx, ny, nz = shape
            sample_keys = [
                (nx // 2, ny // 2, nz // 2, "center", 1),
                (nx // 2 + 1, ny // 2, nz // 2, "center_x_plus", 2),
                (nx // 2, ny // 2 + 1, nz // 2, "center_y_plus", 3),
                (nx // 2, ny // 2, nz // 2 + 1, "center_z_plus", 4),
            ]
        density_center_flat = np.asarray(density_center, dtype=float).reshape(-1) if density_center is not None else np.array([])
        density_mid_flat = np.asarray(density_mid, dtype=float).reshape(-1) if density_mid is not None else np.array([])
        for i, j, k, kind, rank in sample_keys[:128]:
            flat = (int(i) * int(shape[1]) + int(j)) * int(shape[2]) + int(k)
            mid_base = 3 * flat
            gaussian_sample_rows.append({
                "fixture": fixture,
                "salt_M": salt,
                "point_kind": kind,
                "rank": rank,
                "i": i, "j": j, "k": k,
                "density_center": float(density_center_flat[flat]) if flat < density_center_flat.size else "NA",
                "density_x": float(density_mid_flat[mid_base]) if mid_base < density_mid_flat.size else "NA",
                "density_y": float(density_mid_flat[mid_base + 1]) if mid_base + 1 < density_mid_flat.size else "NA",
                "density_z": float(density_mid_flat[mid_base + 2]) if mid_base + 2 < density_mid_flat.size else "NA",
                "eps_x": sample_midpoints(eps_mid, i, j, k, shape)[0],
                "eps_y": sample_midpoints(eps_mid, i, j, k, shape)[1],
                "eps_z": sample_midpoints(eps_mid, i, j, k, shape)[2],
                "rho": sample_array(rho, i, j, k, shape),
            })
        charge_source_rows = charge_full_rows if charge_full_rows else charge_rows
        for row in charge_source_rows[:100000]:
            gaussian_charge_rows.append({
                "fixture": fixture,
                "salt_M": salt,
                "rank": row.get("rank", "NA"),
                "i": row.get("i", "NA"),
                "j": row.get("j", "NA"),
                "k": row.get("k", "NA"),
                "rho": row.get("rho", "NA"),
                "charge_map_scaled": row.get("charge_map_scaled", "NA"),
                "rhs": row.get("rhs", "NA"),
            })

    rpbe_source_rows = []
    rpbe_sample_rows = []
    rpbe_top_rows = []
    if getattr(args, "full_rpb_oracle", False) and str(fixture).startswith("5tif_rpbe"):
        eps_mn, eps_mx, eps_mean = stats(eps_grid)
        charge_mn, charge_mx, charge_mean = stats(rho)
        source_capture = None
        for cap in reversed(rpb_captures):
            if cap.get("source") is not None and not bool(cap.get("vacuum", False)):
                source_capture = cap
                break
        if source_capture is not None:
            source_arr = np.asarray(source_capture.get("source"), dtype=float).reshape(-1)
            coulomb_arr = source_capture.get("coulomb")
            grad_coul_arr = source_capture.get("grad_coulomb")
            grad_eps_arr = source_capture.get("grad_epsilon")
            status = "ok"
            coulomb_norm = finite_l2_norm(coulomb_arr)
            grad_coulomb_norm = finite_l2_norm(grad_coul_arr)
            grad_epsilon_norm = finite_l2_norm(grad_eps_arr)
            source_rms = finite_rms(source_arr)
            source_max = finite_max_abs(source_arr)
            finite_source_ids = np.flatnonzero(np.isfinite(source_arr))
            if finite_source_ids.size:
                ranked_ids = finite_source_ids[np.argsort(np.abs(source_arr[finite_source_ids]))[::-1]][:32]
                nx, ny, nz = shape
                yz = ny * nz
                coulomb_flat = None if coulomb_arr is None else np.asarray(coulomb_arr, dtype=float).reshape(-1)
                grad_coul_flat = None if grad_coul_arr is None else np.asarray(grad_coul_arr, dtype=float).reshape(-1)
                grad_eps_flat = None if grad_eps_arr is None else np.asarray(grad_eps_arr, dtype=float).reshape(-1)
                density_flat = None if density_center is None else np.asarray(density_center, dtype=float).reshape(-1)
                surface_flat = None if surface is None else np.asarray(surface, dtype=float).reshape(-1)
                eps_flat = None if eps_grid is None else np.asarray(eps_grid, dtype=float).reshape(-1)

                def sample_flat(arr, flat, component=None):
                    if arr is None:
                        return "NA"
                    idx = int(flat) if component is None else int(flat) * 3 + int(component)
                    if idx < 0 or idx >= arr.size:
                        return "NA"
                    value = arr[idx]
                    return float(value) if np.isfinite(value) else "NA"

                for rank, flat in enumerate(ranked_ids, start=1):
                    ii = int(flat // yz)
                    rem = int(flat - ii * yz)
                    jj = int(rem // nz)
                    kk = int(rem - jj * nz)
                    rpbe_top_rows.append({
                        "fixture": fixture,
                        "salt_M": salt,
                        "rank": rank,
                        "i": ii,
                        "j": jj,
                        "k": kk,
                        "source": float(source_arr[flat]),
                        "coulomb": sample_flat(coulomb_flat, flat),
                        "grad_coulomb_x": sample_flat(grad_coul_flat, flat, 0),
                        "grad_coulomb_y": sample_flat(grad_coul_flat, flat, 1),
                        "grad_coulomb_z": sample_flat(grad_coul_flat, flat, 2),
                        "grad_epsilon_x": sample_flat(grad_eps_flat, flat, 0),
                        "grad_epsilon_y": sample_flat(grad_eps_flat, flat, 1),
                        "grad_epsilon_z": sample_flat(grad_eps_flat, flat, 2),
                        "density_center": sample_flat(density_flat, flat),
                        "surface": sample_flat(surface_flat, flat),
                        "epsilon_center": sample_flat(eps_flat, flat),
                        "source_status": status,
                    })
        else:
            source_arr = np.array([])
            status = "ingredient_oracle"
            coulomb_norm = "NA"
            grad_coulomb_norm = "NA"
            grad_epsilon_norm = "NA"
            source_rms = "NA"
            source_max = "NA"
        rpbe_source_rows.append({
            "fixture": fixture,
            "salt_M": salt,
            "source_status": status,
            "coulomb_norm": coulomb_norm,
            "grad_coulomb_norm": grad_coulomb_norm,
            "grad_epsilon_norm": grad_epsilon_norm,
            "source_rms": source_rms,
            "source_max": source_max,
            "epsilon_min": eps_mn,
            "epsilon_max": eps_mx,
            "epsilon_mean": eps_mean,
            "charge_min": charge_mn,
            "charge_max": charge_mx,
            "charge_mean": charge_mean,
            "charge_nonzero_cells": finite_nonzero(rho),
        })
        for row in gaussian_sample_rows[:128]:
            out_row = dict(row)
            try:
                flat = (int(float(row["i"])) * int(shape[1]) + int(float(row["j"]))) * int(shape[2]) + int(float(row["k"]))
            except Exception:
                flat = -1
            out_row["source"] = float(source_arr[flat]) if 0 <= flat < source_arr.size else "NA"
            out_row["source_status"] = status
            rpbe_sample_rows.append(out_row)

    focus_boundary_rows = []
    focus_probe_rows = []
    if getattr(args, "full_focus_probe_oracle", False) and "focus" in str(fixture).lower():
        for row in field_rows[:256]:
            if row.get("i") in ("0", 0) or row.get("j") in ("0", 0) or row.get("k") in ("0", 0):
                focus_boundary_rows.append({
                    "fixture": fixture,
                    "salt_M": salt,
                    "i": row.get("i", "NA"),
                    "j": row.get("j", "NA"),
                    "k": row.get("k", "NA"),
                    "boundary_phi": row.get("water_phi", "NA"),
                    "field_source": row.get("field_source", "NA"),
                    "parent_grid": "NA",
                    "child_grid": f"{shape[0]}x{shape[1]}x{shape[2]}",
                    "interpolation_status": "sampled",
                })
        if not focus_boundary_rows and water_phi is not None:
            nx, ny, nz = shape
            rank = 0
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        if not (i == 0 or j == 0 or k == 0 or i == nx - 1 or j == ny - 1 or k == nz - 1):
                            continue
                        focus_boundary_rows.append({
                            "fixture": fixture,
                            "salt_M": salt,
                            "i": i,
                            "j": j,
                            "k": k,
                            "boundary_phi": sample_array(water_phi, i, j, k, shape),
                            "field_source": water_field_source,
                            "parent_grid": "137x137x137",
                            "child_grid": f"{shape[0]}x{shape[1]}x{shape[2]}",
                            "interpolation_status": f"oracle_from_{water_field_source}",
                        })
                        rank += 1
                        if rank >= 512:
                            break
                    if rank >= 512:
                        break
                if rank >= 512:
                    break
        energy = read_energy_summary(run_dir / "pydelphi_energy.csv")
        probe_energy_kt = energy.get("E_total_probe", "NA")
        if probe_energy_kt == "NA":
            probe_energy_kt = read_frc_total_energy(run_dir / "oracle_frc.out")
        probe_site = "a,p,q,f"
        if "site_p" in fixture_text:
            probe_site = "p"
        elif "site_q" in fixture_text:
            probe_site = "q"
        elif "residue_only" in fixture_text:
            probe_site = ""
        focus_probe_rows.append({
            "fixture": fixture,
            "salt_M": salt,
            "probe_residue": "ASP39",
            "probe_site": probe_site,
            "matched_atoms": "NA",
            "E_probe_kT": probe_energy_kt,
            "E_probe_kJmol": "NA",
            "E_probe_kcalmol": "NA",
            "status": "ok" if probe_energy_kt != "NA" else "missing_probe_energy",
            "parent_grid": "NA",
            "child_grid": f"{shape[0]}x{shape[1]}x{shape[2]}",
        })

    eps_min, eps_max, eps_mean = stats(eps_grid)
    surf_arr = np.asarray(surface).reshape(-1) if surface is not None else np.array([])
    rho_arr = np.asarray(rho).reshape(-1) if rho is not None else np.array([])
    surface_rows = []
    export_full_surface = (
        not str(fixture).startswith("5tif")
        or getattr(args, "full_vdwms_state", False)
        or getattr(args, "full_midpoints", False)
        or getattr(args, "full_boundary_oracle", False)
    )
    if surf_arr.size and export_full_surface:
        solute_ids = np.flatnonzero(surf_arr < 0.5)
        nx, ny, nz = shape
        yz = ny * nz
        for idx in solute_ids:
            ii = int(idx // yz)
            rem = int(idx - ii * yz)
            jj = int(rem // nz)
            kk = int(rem - jj * nz)
            surface_rows.append({
                "fixture": fixture,
                "i": ii,
                "j": jj,
                "k": kk,
                "inside": 1,
                "surface": float(surf_arr[idx]),
                "surface_semantics": SURFACE_SEMANTICS,
            })
    eps_mid_stats = ("NA", "NA", "NA")
    eps_six_stats = ("NA", "NA", "NA")
    eps_axis = {}
    if isinstance(eps_mid, (list, tuple)) and len(eps_mid) >= 3 and all(x is not None for x in eps_mid[:3]):
        mids = np.concatenate([np.asarray(x, dtype=float).reshape(-1) for x in eps_mid[:3]])
        eps_mid_stats = stats(mids)
        eps_six_stats = epsilon_six_face_stats(eps_mid, shape)
        eps_axis = epsilon_axis_stats(eps_mid)
    full_midpoint_for_fixture = (
        getattr(args, "full_midpoints", False)
        and salt_matches(salt, getattr(args, "full_midpoint_salt", "0.15"))
        and eps_mid is not None
    )
    atom_radius_rows, atom_radius_stats = atom_radius_audit(ctx, fixture)
    atom_accessibility = safe_attr(ctx, "atom_accessibility")
    if atom_accessibility is None and vdwms_capture:
        atom_accessibility = vdwms_capture.get("atom_accessibility")
    try:
        atom_accessible_count = int(np.count_nonzero(np.asarray(atom_accessibility).reshape(-1) == 0)) if atom_accessibility is not None else "NA"
    except Exception:
        atom_accessible_count = "NA"
    energy = read_energy_summary(run_dir / "pydelphi_energy.csv")
    e_rxn_kT = energy.get("E_rxn_corr_tot", "NA")
    try:
        e_rxn_kJ = float(e_rxn_kT) * 0.00831446261815324 * 297.15
    except Exception:
        e_rxn_kJ = "NA"
    stage = {
        "fixture": fixture,
        "salt_M": salt,
        "midpoint_layout": MIDPOINT_LAYOUT,
        "midpoint_positive_face_order": "eps_x=i+1/2,eps_y=j+1/2,eps_z=k+1/2",
        "midpoint_full_export": bool(full_midpoint_for_fixture),
        "grid_nx": shape[0],
        "grid_ny": shape[1],
        "grid_nz": shape[2],
        "grid_origin_x_A": origin[0],
        "grid_origin_y_A": origin[1],
        "grid_origin_z_A": origin[2],
        "grid_center_x_A": center[0],
        "grid_center_y_A": center[1],
        "grid_center_z_A": center[2],
        "grid_spacing_A": h,
        "total_charge": getattr(ctx, "total_charge", "NA"),
        "charge_assigned": float(getattr(ctx, "total_charge", "nan")),
        "pydelphi_atom_count": atom_radius_stats.get("pydelphi_atom_count", "NA"),
        "pydelphi_threads": 2,
        "pydelphi_radius_sum_A": atom_radius_stats.get("pydelphi_radius_sum_A", "NA"),
        "pydelphi_radius_min_A": atom_radius_stats.get("pydelphi_radius_min_A", "NA"),
        "pydelphi_radius_max_A": atom_radius_stats.get("pydelphi_radius_max_A", "NA"),
        "pydelphi_radius_mean_A": atom_radius_stats.get("pydelphi_radius_mean_A", "NA"),
        "surface_semantics": SURFACE_SEMANTICS,
        "surface_points": int(np.count_nonzero(surf_arr)) if surf_arr.size else "NA",
        "surface_zero_points": int(np.count_nonzero(surf_arr < 0.5)) if surf_arr.size else "NA",
        "surface_one_points": int(np.count_nonzero(surf_arr >= 0.5)) if surf_arr.size else "NA",
        "oracle_solute_points": int(np.count_nonzero(surf_arr < 0.5)) if surf_arr.size else "NA",
        "oracle_solvent_points": int(np.count_nonzero(surf_arr >= 0.5)) if surf_arr.size else "NA",
        "eps_sum_min": eps_six_stats[0],
        "eps_sum_max": eps_six_stats[1],
        "eps_sum_mean": eps_six_stats[2],
        "eps_mid_min": eps_mid_stats[0],
        "eps_mid_max": eps_mid_stats[1],
        "eps_mid_mean": eps_mid_stats[2],
        "eps_mid_x_min": eps_axis.get("eps_mid_x_min", "NA"),
        "eps_mid_x_max": eps_axis.get("eps_mid_x_max", "NA"),
        "eps_mid_x_mean": eps_axis.get("eps_mid_x_mean", "NA"),
        "eps_mid_x_count": eps_axis.get("eps_mid_x_count", "NA"),
        "eps_mid_x_internal_count": eps_axis.get("eps_mid_x_internal_count", "NA"),
        "eps_mid_y_min": eps_axis.get("eps_mid_y_min", "NA"),
        "eps_mid_y_max": eps_axis.get("eps_mid_y_max", "NA"),
        "eps_mid_y_mean": eps_axis.get("eps_mid_y_mean", "NA"),
        "eps_mid_y_count": eps_axis.get("eps_mid_y_count", "NA"),
        "eps_mid_y_internal_count": eps_axis.get("eps_mid_y_internal_count", "NA"),
        "eps_mid_z_min": eps_axis.get("eps_mid_z_min", "NA"),
        "eps_mid_z_max": eps_axis.get("eps_mid_z_max", "NA"),
        "eps_mid_z_mean": eps_axis.get("eps_mid_z_mean", "NA"),
        "eps_mid_z_count": eps_axis.get("eps_mid_z_count", "NA"),
        "eps_mid_z_internal_count": eps_axis.get("eps_mid_z_internal_count", "NA"),
        "vdwms_boundary_initial": vdwms_capture.get("boundary_initial", "NA") if vdwms_capture else "NA",
        "vdwms_boundary_final": vdwms_capture.get("boundary_final", getattr(ctx, "num_boundary_grid_points", "NA")) if vdwms_capture else getattr(ctx, "num_boundary_grid_points", "NA"),
        "vdwms_added": "NA",
        "vdwms_removed": "NA",
        "vdwms_exposed_points": vdwms_capture.get("exposed_points", getattr(ctx, "num_exposed_grids", "NA")) if vdwms_capture else getattr(ctx, "num_exposed_grids", "NA"),
        "vdwms_accessible_atoms": atom_accessible_count,
        "vdwms_external_boundary_points": vdwms_capture.get("external_boundary_final", getattr(ctx, "num_external_boundary_points", "NA")) if vdwms_capture else getattr(ctx, "num_external_boundary_points", "NA"),
        "pydelphi_total_seconds": getattr(app, "_titan_wall_total_seconds", "NA"),
        "pydelphi_setup_seconds": getattr(app, "_titan_wall_setup_seconds", "NA"),
        "pydelphi_surface_seconds": timing_value(app, "space running", "surf"),
        "pydelphi_epsilon_seconds": timing_value(app, "epsilon"),
        "pydelphi_solver_seconds": timing_value(app, "solving pbe", "solver"),
        "pydelphi_energy_seconds": timing_value(app, "energy"),
        "pydelphi_map_io_seconds": timing_value(app, "writing"),
        "pydelphi_iterations": getattr(ctx, "iteration", getattr(ctx, "iterations", "NA")),
        "final_E_rxn_kT": e_rxn_kT,
        "final_E_rxn_kJmol": e_rxn_kJ,
        "iscm_water_rf_kT": energy_component_rows[0].get("water_rf_kT", "NA") if energy_component_rows else "NA",
        "iscm_water_rf_kJmol": energy_component_rows[0].get("water_rf_kJmol", "NA") if energy_component_rows else "NA",
        "iscm_vacuum_rf_kT": energy_component_rows[0].get("vacuum_rf_kT", "NA") if energy_component_rows else "NA",
        "iscm_vacuum_rf_kJmol": energy_component_rows[0].get("vacuum_rf_kJmol", "NA") if energy_component_rows else "NA",
        "iscm_corrected_rxn_kT": energy_component_rows[0].get("total_reactionfield_kT", "NA") if energy_component_rows else "NA",
        "iscm_corrected_rxn_kJmol": energy_component_rows[0].get("total_reactionfield_kJmol", "NA") if energy_component_rows else "NA",
        "iscm_finalization_mode": energy_component_rows[0].get("finalization_mode", "NA") if energy_component_rows else "NA",
        "rtol": "1e-6",
        "atol": "1e-8",
    }
    midpoint_rows = []
    midpoint_seen = set()
    def sample_midpoint_index_values(i, j, k):
        if idx_mid is None:
            return ("NA", "NA", "NA")
        arr = np.asarray(idx_mid).reshape(-1)
        base = flat_index(i, j, k, shape) * 3
        if base + 2 >= arr.size:
            return ("NA", "NA", "NA")
        return (int(arr[base]), int(arr[base + 1]), int(arr[base + 2]))

    def add_midpoint(kind, rank, i, j, k):
        key = (int(i), int(j), int(k))
        if key in midpoint_seen:
            return
        midpoint_seen.add(key)
        ex, ey, ez, es = sample_midpoints(eps_mid, i, j, k, shape)
        if ex == "NA":
            return
        ix, iy, iz = sample_midpoint_index_values(i, j, k)
        midpoint_rows.append({
            "fixture": fixture,
            "point_kind": kind,
            "rank": rank,
            "i": int(i), "j": int(j), "k": int(k),
            "eps_x": ex, "eps_y": ey, "eps_z": ez, "eps_sum": es,
            "idx_x": ix, "idx_y": iy, "idx_z": iz,
        })

    if full_midpoint_for_fixture:
        nx, ny, nz = shape
        rank = 0
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz):
                    rank += 1
                    add_midpoint("full_grid", rank, ii, jj, kk)
    else:
        rank = 0
        for row in list(sample_field) + list(sample_boundary) + list(sample_charge):
            rank += 1
            add_midpoint(row.get("point_kind", row.get("boundary_kind", row.get("charge_kind", "sample"))), rank, row["i"], row["j"], row["k"])

    vdwms_audit_rows = []
    midpoint_state_rows = []
    exposed_rows = []
    accessibility_rows = []
    if getattr(args, "full_vdwms_state", False) and vdwms_capture:
        epsdim_value = vdwms_capture.get("epsilon_dimension", getattr(ctx, "epsilon_dimension", "NA"))
        try:
            epsdim = int(epsdim_value)
        except Exception:
            epsdim = 0
        if epsdim <= 1:
            atom_count = atom_radius_stats.get("pydelphi_atom_count", "NA")
            try:
                epsdim = int(atom_count) + 2
            except Exception:
                epsdim = int(shape[0] * shape[1] * shape[2] + 2)
        raw_index = vdwms_capture.get("raw_index")
        final_index = vdwms_capture.get("final_index")
        raw_counts = index_axis_counts(raw_index, epsdim)
        final_counts = index_axis_counts(final_index, epsdim)
        atom_accessibility = vdwms_capture.get("atom_accessibility")
        try:
            accessible_atoms = int(np.count_nonzero(np.asarray(atom_accessibility).reshape(-1) == 0)) if atom_accessibility is not None else "NA"
        except Exception:
            accessible_atoms = "NA"
        vdwms_audit_rows.append({
            "fixture": fixture,
            "epsilon_dimension": epsdim,
            "raw_internal_face_count_x": raw_counts["x"],
            "raw_internal_face_count_y": raw_counts["y"],
            "raw_internal_face_count_z": raw_counts["z"],
            "final_internal_face_count_x": final_counts["x"],
            "final_internal_face_count_y": final_counts["y"],
            "final_internal_face_count_z": final_counts["z"],
            "boundary_initial": vdwms_capture.get("boundary_initial", "NA"),
            "boundary_final": vdwms_capture.get("boundary_final", "NA"),
            "external_boundary_initial": vdwms_capture.get("external_boundary_initial", "NA"),
            "external_boundary_final": vdwms_capture.get("external_boundary_final", "NA"),
            "exposed_points": vdwms_capture.get("exposed_points", "NA"),
            "accessible_atoms": accessible_atoms,
        })
        if full_midpoint_for_fixture and final_index is not None:
            nx, ny, nz = shape
            rank = 0
            for ii in range(nx):
                for jj in range(ny):
                    for kk in range(nz):
                        rank += 1
                        rx, ry, rz = sample_index_triplet(raw_index, ii, jj, kk, shape)
                        fx, fy, fz = sample_index_triplet(final_index, ii, jj, kk, shape)
                        ex, ey, ez, es = sample_midpoints(eps_mid, ii, jj, kk, shape)
                        midpoint_state_rows.append({
                            "fixture": fixture,
                            "rank": rank,
                            "i": ii, "j": jj, "k": kk,
                            "raw_idx_x": rx, "raw_idx_y": ry, "raw_idx_z": rz,
                            "idx_x": fx, "idx_y": fy, "idx_z": fz,
                            "raw_media_x": media_from_index(rx, epsdim),
                            "raw_media_y": media_from_index(ry, epsdim),
                            "raw_media_z": media_from_index(rz, epsdim),
                            "media_x": media_from_index(fx, epsdim),
                            "media_y": media_from_index(fy, epsdim),
                            "media_z": media_from_index(fz, epsdim),
                            "eps_x": ex, "eps_y": ey, "eps_z": ez,
                        })
        exposed = vdwms_capture.get("exposed_grids_coords")
        if exposed is not None:
            arr = np.asarray(exposed)
            limit = min(arr.shape[0], 200)
            for idx in range(limit):
                exposed_rows.append({
                    "fixture": fixture,
                    "rank": idx,
                    "x": float(arr[idx, 0]),
                    "y": float(arr[idx, 1]),
                    "z": float(arr[idx, 2]),
                })
        if atom_accessibility is not None:
            arr = np.asarray(atom_accessibility).reshape(-1)
            for idx, value in enumerate(arr[: min(arr.size, 5000)]):
                accessibility_rows.append({
                    "fixture": fixture,
                    "atom_index": idx,
                    "atom_accessibility": int(value),
                })

    return field_rows, field_map_rows, boundary_rows, charge_rows, charge_full_rows, surface_rows, midpoint_rows, atom_radius_rows, full_midpoint_for_fixture, stage, vdwms_audit_rows, midpoint_state_rows, exposed_rows, accessibility_rows, boundary_map_rows, boundary_inputs_rows, boundary_formula_rows, operator_rows, sor_trace_rows, nwt_trace_rows, solver_summary_rows, nwt_summary_rows, iscm_transition_rows, iscm_sample_rows, iscm_full_rows, iscm_energy_term_rows, energy_component_rows, gaussian_summary_rows, gaussian_sample_rows, gaussian_charge_rows, rpbe_source_rows, rpbe_sample_rows, rpbe_top_rows, focus_boundary_rows, focus_probe_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pyroot", required=True)
    ap.add_argument("--workroot", required=True)
    ap.add_argument("--samples", default="")
    ap.add_argument("--full-midpoints", action="store_true")
    ap.add_argument("--full-midpoint-salt", default="0.15")
    ap.add_argument("--full-vdwms-state", action="store_true")
    ap.add_argument("--full-boundary-oracle", action="store_true")
    ap.add_argument("--full-operator-oracle", action="store_true")
    ap.add_argument("--full-charge-grid", action="store_true")
    ap.add_argument("--full-field-oracle", action="store_true")
    ap.add_argument("--full-iscm-energy-oracle", action="store_true")
    ap.add_argument("--full-gaussian-oracle", action="store_true")
    ap.add_argument("--full-rpb-oracle", action="store_true")
    ap.add_argument("--full-focus-probe-oracle", action="store_true")
    ap.add_argument("--full-field-map", action="store_true")
    ap.add_argument("--fixture", action="append", default=[])
    ap.add_argument("--nwt-only", action="store_true")
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    work = Path(args.workroot)
    report_rows = []
    all_field, all_field_map, all_boundary, all_charge, all_charge_full, all_surface, all_atom_radius = [], [], [], [], [], [], []
    all_midpoint_full, all_midpoint_samples, all_stage = [], [], []
    all_vdwms_audit, all_midpoint_state = [], []
    all_exposed, all_accessibility = [], []
    all_boundary_map, all_boundary_inputs, all_boundary_formula = [], [], []
    all_operator_stage, all_sor_trace, all_nwt_trace, all_solver_summary, all_nwt_solver_summary = [], [], [], [], []
    all_iscm_transition, all_iscm_samples = [], []
    all_iscm_full, all_iscm_energy_terms, all_energy_components = [], [], []
    all_gaussian_summary, all_gaussian_samples, all_gaussian_charge = [], [], []
    all_rpbe_source, all_rpbe_samples, all_rpbe_top = [], [], []
    all_focus_boundary, all_focus_probe = [], []
    status = 0
    try:
        fixtures_to_run = []
        linear_requested = (
            args.full_midpoints or args.full_vdwms_state or
            args.full_boundary_oracle or args.full_charge_grid or
            args.full_iscm_energy_oracle
        )
        if args.nwt_only:
            fixtures_to_run.extend(NWT_FIXTURES)
        elif linear_requested:
            fixtures_to_run.extend(LINEAR_FIXTURES)
        elif not (args.full_gaussian_oracle or args.full_rpb_oracle or args.full_focus_probe_oracle):
            fixtures_to_run.extend(LINEAR_FIXTURES)
        if (args.full_operator_oracle or args.full_field_oracle) and not args.nwt_only:
            fixtures_to_run.extend(NWT_FIXTURES)
        if args.full_gaussian_oracle and not args.nwt_only:
            fixtures_to_run.extend(GAUSSIAN_FIXTURES)
        if args.full_rpb_oracle and not args.nwt_only:
            fixtures_to_run.extend(RPBE_FIXTURES)
        if args.full_focus_probe_oracle and not args.nwt_only:
            fixtures_to_run.extend(FOCUS_FIXTURES)
        seen_fixtures = set()
        fixtures_to_run = [row for row in fixtures_to_run if not (row in seen_fixtures or seen_fixtures.add(row))]
        if args.fixture:
            selected = set(str(x) for x in args.fixture)
            fixtures_to_run = [row for row in fixtures_to_run if row[0] in selected]
            missing = sorted(selected.difference(row[0] for row in fixtures_to_run))
            if missing:
                raise RuntimeError("Requested fixture(s) not available for selected oracle flags: " + ",".join(missing))
        for fixture, salt in fixtures_to_run:
            try:
                field, field_map, boundary, charge, charge_full, surface, midpoint, atom_radius, midpoint_is_full, stage, vdwms_audit, midpoint_state, exposed, accessibility, boundary_map, boundary_inputs, boundary_formula, operator_stage, sor_trace, nwt_trace, solver_summary, nwt_summary, iscm_transition, iscm_samples, iscm_full, iscm_energy_terms, energy_components, gaussian_summary, gaussian_samples, gaussian_charge, rpbe_source, rpbe_samples, rpbe_top, focus_boundary, focus_probe = export_fixture(args, fixture, salt)
                all_field.extend(field)
                all_field_map.extend(field_map)
                all_boundary.extend(boundary)
                all_charge.extend(charge)
                all_charge_full.extend(charge_full)
                all_surface.extend(surface)
                all_atom_radius.extend(atom_radius)
                all_vdwms_audit.extend(vdwms_audit)
                all_midpoint_state.extend(midpoint_state)
                all_exposed.extend(exposed)
                all_accessibility.extend(accessibility)
                all_boundary_map.extend(boundary_map)
                all_boundary_inputs.extend(boundary_inputs)
                all_boundary_formula.extend(boundary_formula)
                all_operator_stage.extend(operator_stage)
                all_sor_trace.extend(sor_trace)
                all_nwt_trace.extend(nwt_trace)
                all_solver_summary.extend(solver_summary)
                all_nwt_solver_summary.extend(nwt_summary)
                all_iscm_transition.extend(iscm_transition)
                all_iscm_samples.extend(iscm_samples)
                all_iscm_full.extend(iscm_full)
                all_iscm_energy_terms.extend(iscm_energy_terms)
                all_energy_components.extend(energy_components)
                all_gaussian_summary.extend(gaussian_summary)
                all_gaussian_samples.extend(gaussian_samples)
                all_gaussian_charge.extend(gaussian_charge)
                all_rpbe_source.extend(rpbe_source)
                all_rpbe_samples.extend(rpbe_samples)
                all_rpbe_top.extend(rpbe_top)
                all_focus_boundary.extend(focus_boundary)
                all_focus_probe.extend(focus_probe)
                if midpoint_is_full:
                    all_midpoint_full.extend(midpoint)
                else:
                    all_midpoint_samples.extend(midpoint)
                all_stage.append(stage)
                charge_numeric = has_numeric_column(charge, "rho") or has_numeric_column(charge_full, "rho")
                field_source = field[0].get("field_source", "missing") if field else "missing"
                field_precision = field[0].get("field_precision_status", "missing") if field else "missing"
                report_rows.append({
                    "fixture": fixture,
                    "salt_M": salt,
                    "status": "ok",
                    "detail": "exported",
                    "charge_grid": "present" if charge_numeric else "missing",
                    "surface_map": "present" if surface else "missing",
                    "midpoint_epsilon": "present" if midpoint else "missing",
                    "boundary_map": "present" if boundary_map else "missing",
                    "operator_stage": "present" if operator_stage else "missing",
                    "boundary_samples": "present" if boundary else "missing",
                    "field_samples": "present" if field else "missing",
                    "field_source": field_source,
                    "field_precision_status": field_precision,
                    "nwt_trace": "present" if nwt_trace else "missing",
                    "nwt_solver_summary": "present" if nwt_summary else "missing",
                    "iscm_energy_oracle": "present" if energy_components else "missing",
                    "gaussian_oracle": "present" if gaussian_summary else "missing",
                    "rpbe_oracle": "present" if rpbe_source else "missing",
                    "focus_probe_oracle": "present" if focus_probe or focus_boundary else "missing",
                    "stage_oracle": "present",
                    "surface_semantics": SURFACE_SEMANTICS,
                })
            except Exception as exc:
                status = 2
                report_rows.append({
                    "fixture": fixture,
                    "salt_M": salt,
                    "status": "failed",
                    "detail": f"{type(exc).__name__}: {exc}",
                    "charge_grid": "missing",
                    "surface_map": "missing",
                    "midpoint_epsilon": "missing",
                    "boundary_map": "missing",
                    "operator_stage": "missing",
                    "boundary_samples": "missing",
                    "field_samples": "missing",
                    "field_source": "missing",
                    "field_precision_status": "missing",
                    "nwt_trace": "missing",
                    "nwt_solver_summary": "missing",
                    "iscm_energy_oracle": "missing",
                    "gaussian_oracle": "missing",
                    "rpbe_oracle": "missing",
                    "focus_probe_oracle": "missing",
                    "stage_oracle": "missing",
                    "surface_semantics": "NA",
                })
                (work / fixture / "pydelphi_exception.txt").parent.mkdir(parents=True, exist_ok=True)
                (work / fixture / "pydelphi_exception.txt").write_text(traceback.format_exc(), encoding="utf-8")
                if args.strict:
                    raise
    except Exception:
        if args.strict:
            raise

    if all_field:
        field_columns = [
            "fixture","point_kind","rank","i","j","k","rho","surface","epsilon",
            "inside","surface_semantics","eps_x","eps_y","eps_z","eps_sum","eps_sum6",
            "solvent_phi","water_phi","vacuum_phi","reaction_phi",
            "field_source","vacuum_field_source","field_precision_status"
        ]
        write_csv(work / "pbe_oracle_field_samples.csv", all_field, field_columns)
        if getattr(args, "full_field_oracle", False):
            write_csv(work / "pbe_oracle_field_samples_fullprecision.csv", all_field, field_columns)
    if all_field_map:
        write_csv(work / "pbe_oracle_field_map_fullprecision.csv", all_field_map, [
            "fixture","rank","i","j","k","solvent_phi","water_phi",
            "vacuum_phi","reaction_phi","field_source",
            "vacuum_field_source","field_precision_status"
        ])
    if all_boundary:
        write_csv(work / "pbe_oracle_boundary_samples.csv", all_boundary, [
            "fixture","boundary_kind","rank","i","j","k","surface","epsilon",
            "inside","surface_semantics","eps_x","eps_y","eps_z","eps_sum","eps_sum6",
            "solvent_phi","water_phi","vacuum_phi"
        ])
    if all_boundary_map:
        boundary_fields = [
            "fixture","salt_M","capture_index","kind","vacuum","i","j","k",
            "flat_id","phi","grid_spacing","exdi_scaled","indi_scaled",
            "debye_length","non_zero_salt","epkt",
            "atom_grid_x","atom_grid_y","atom_grid_z"
        ]
        write_csv(work / "pbe_oracle_boundary_map.csv", all_boundary_map, boundary_fields)
    if all_boundary_inputs:
        write_csv(work / "pbe_oracle_boundary_inputs.csv", all_boundary_inputs, [
            "fixture","salt_M","capture_index","kind","vacuum","grid_spacing",
            "exdi_scaled","indi_scaled","debye_length","non_zero_salt",
            "epkt","boundary_points"
        ])
    if all_boundary_formula:
        write_csv(work / "pbe_oracle_boundary_formula_samples.csv", all_boundary_formula, [
            "fixture","salt_M","capture_index","kind","vacuum","i","j","k",
            "flat_id","phi","grid_spacing","exdi_scaled","indi_scaled",
            "debye_length","non_zero_salt","epkt",
            "atom_grid_x","atom_grid_y","atom_grid_z"
        ])
    if all_surface:
        write_csv(work / "pbe_oracle_surface_map.csv", all_surface, [
            "fixture","i","j","k","inside","surface","surface_semantics"
        ])
    if all_midpoint_samples:
        write_csv(work / "pbe_oracle_midpoint_epsilon_samples.csv", all_midpoint_samples, [
            "fixture","point_kind","rank","i","j","k","eps_x","eps_y","eps_z","eps_sum","idx_x","idx_y","idx_z"
        ])
    if all_midpoint_full:
        write_csv(work / "pbe_oracle_midpoint_epsilon.csv", all_midpoint_full, [
            "fixture","point_kind","rank","i","j","k","eps_x","eps_y","eps_z","eps_sum","idx_x","idx_y","idx_z"
        ])
    if all_atom_radius:
        write_csv(work / "pbe_oracle_atom_radius_audit.csv", all_atom_radius, [
            "fixture","atom_index","x_A","y_A","z_A","q","radius_A","media_id"
        ])
    charge_for_replay = all_charge_full if all_charge_full else all_charge
    if charge_for_replay and has_numeric_column(charge_for_replay, "rho"):
        write_csv(work / "pbe_oracle_charge_grid.csv", charge_for_replay, [
            "fixture","salt_M","capture_index","charge_kind","rank","i","j","k",
            "rho","charge_map_scaled","rhs","surface","epsilon",
            "inside","surface_semantics","eps_x","eps_y","eps_z","eps_sum","eps_sum6"
        ])
    if all_charge_full and has_numeric_column(all_charge_full, "rho"):
        write_csv(work / "pbe_oracle_charge_grid_full.csv", all_charge_full, [
            "fixture","salt_M","capture_index","charge_kind","rank","i","j","k",
            "rho","charge_map_scaled","rhs","surface","epsilon",
            "inside","surface_semantics","eps_x","eps_y","eps_z","eps_sum","eps_sum6"
        ])
    if all_vdwms_audit:
        write_csv(work / "pbe_oracle_vdwms_audit.csv", all_vdwms_audit, [
            "fixture","epsilon_dimension",
            "raw_internal_face_count_x","raw_internal_face_count_y","raw_internal_face_count_z",
            "final_internal_face_count_x","final_internal_face_count_y","final_internal_face_count_z",
            "boundary_initial","boundary_final",
            "external_boundary_initial","external_boundary_final",
            "exposed_points","accessible_atoms"
        ])
    if all_midpoint_state:
        write_csv(work / "pbe_oracle_midpoint_state.csv", all_midpoint_state, [
            "fixture","rank","i","j","k",
            "raw_idx_x","raw_idx_y","raw_idx_z","idx_x","idx_y","idx_z",
            "raw_media_x","raw_media_y","raw_media_z",
            "media_x","media_y","media_z","eps_x","eps_y","eps_z"
        ])
    if all_exposed:
        write_csv(work / "pbe_oracle_exposed_points.csv", all_exposed, [
            "fixture","rank","x","y","z"
        ])
    if all_accessibility:
        write_csv(work / "pbe_oracle_atom_accessibility.csv", all_accessibility, [
            "fixture","atom_index","atom_accessibility"
        ])
    if all_operator_stage:
        write_csv(work / "pbe_oracle_operator_stage.csv", all_operator_stage, [
            "fixture","salt_M","capture_index","vacuum","point_kind","rank",
            "i","j","k","rho","charge_map_scaled","rhs","screen",
            "eps_i_plus","eps_i_minus","eps_j_plus","eps_j_minus",
            "eps_k_plus","eps_k_minus","eps_sum_no_screen","denom",
            "boundary_flag","ion_exclusion","phi_initial","phi_water",
            "grid_spacing","debye_length","epkt","operator_units"
        ])
    if all_sor_trace:
        write_csv(work / "pbe_oracle_sor_trace.csv", all_sor_trace, [
            "fixture","salt_M","call","phase","even_odd","i","j","k",
            "omega","old_phi","neighbor_sum","charge_density","denom",
            "new_phi","new_phi_expected","dphi","rms_sum_sq","max_dphi_block"
        ])
    if all_nwt_trace:
        write_csv(work / "pbe_oracle_nwt_trace.csv", all_nwt_trace, [
            "fixture","salt_M","call","phase","vacuum","iteration","color",
            "i","j","k","omega","old_phi","neighbor_sum","charge_density",
            "eps_sum","screen","sinh_phi","cosh_phi","linearized_rhs",
            "denom","target_phi","new_phi","new_phi_expected","dphi",
            "rms_sum_sq","max_dphi_block","field_units","operator_convention"
        ])
    if all_solver_summary:
        write_csv(work / "pbe_oracle_solver_summary.csv", all_solver_summary, [
            "fixture","salt_M","capture_index","vacuum","grid_spacing",
            "debye_length","non_zero_salt","epkt","operator_units",
            "operator_rows","sor_trace_rows"
        ])
    if all_nwt_solver_summary:
        write_csv(work / "pbe_oracle_nwt_solver_summary.csv", all_nwt_solver_summary, [
            "fixture","salt_M","vacuum","temp_K","omega","block_size",
            "max_nonlinear_iters","rms_tol","dphi_tol","check_dphi",
            "num_cpu_threads","rmsd","max_dphi","iterations",
            "stop_reason","trace_rows","field_units"
        ])
    if all_iscm_transition:
        write_csv(work / "pbe_oracle_iscm_transition.csv", all_iscm_transition, [
            "fixture","salt_M","capture_index","boundary_count",
            "charged_grid_count","induced_count","induced_charge_sum",
            "induced_charge_l1","induced_charge_l2",
            "induced_charge_max_abs","water_rf_energy","thread_sum_raw",
            "thread_sum_kT","scale",
            "indi_scaled","debye_factor","epkt"
        ])
    if all_iscm_samples:
        write_csv(work / "pbe_oracle_iscm_surface_charge_samples.csv", all_iscm_samples, [
            "fixture","salt_M","capture_index","rank","i","j","k",
            "surface_x","surface_y","surface_z","induced_charge",
            "atom_potential_sum","energy_contribution_kT"
        ])
    if all_iscm_full:
        write_csv(work / "pbe_oracle_iscm_full.csv", all_iscm_full, [
            "fixture","salt_M","capture_index","rank","i","j","k",
            "surface_x","surface_y","surface_z","induced_charge",
            "atom_potential_sum","energy_contribution_kT"
        ])
    if all_iscm_energy_terms:
        write_csv(work / "pbe_oracle_iscm_energy_terms.csv", all_iscm_energy_terms, [
            "fixture","salt_M","capture_index","rank","i","j","k",
            "induced_charge","atom_potential_sum","energy_contribution_kT"
        ])
    if all_energy_components:
        write_csv(work / "pbe_oracle_energy_components.csv", all_energy_components, [
            "fixture","salt_M","water_rf_kT","water_rf_kJmol",
            "vacuum_rf_kT","vacuum_rf_kJmol","total_reactionfield_kT",
            "total_reactionfield_kJmol","nonpolar_kT","overall_solvation_kT",
            "final_selected_kT","final_selected_kJmol","finalization_mode",
            "epkt","scale","indi_scaled","exdi_scaled","debye_length",
            "temperature_K"
        ])
    if all_gaussian_summary:
        write_csv(work / "pbe_oracle_gaussian_density_summary.csv", all_gaussian_summary, [
            "fixture","salt_M","density_source","midpoint_density_source",
            "grid_nx","grid_ny","grid_nz","density_sum","density_min",
            "gaussian_density_sum",
            "density_max","density_mean","midpoint_density_sum",
            "midpoint_density_min","midpoint_density_max","midpoint_density_mean",
            "density_cutoff","density_cutoff_points","gaussian_cutoff_points",
            "midpoint_cutoff_points",
            "eps_mid_x_mean","eps_mid_y_mean","eps_mid_z_mean",
            "eps_mid_min","eps_mid_max","eps_mid_mean",
            "gaussian_eps_min","gaussian_eps_max","gaussian_eps_mean",
            "charge_nonzero_cells","charge_sum","charge_grid_source","status"
        ])
    if all_gaussian_samples:
        write_csv(work / "pbe_oracle_gaussian_samples.csv", all_gaussian_samples, [
            "fixture","salt_M","point_kind","rank","i","j","k",
            "density_center","density_x","density_y","density_z",
            "eps_x","eps_y","eps_z","rho"
        ])
    if all_gaussian_charge:
        write_csv(work / "pbe_oracle_gaussian_charge_grid.csv", all_gaussian_charge, [
            "fixture","salt_M","rank","i","j","k","rho","charge_map_scaled","rhs"
        ])
    if all_rpbe_source:
        write_csv(work / "pbe_oracle_rpbe_source.csv", all_rpbe_source, [
            "fixture","salt_M","source_status","coulomb_norm",
            "grad_coulomb_norm","grad_epsilon_norm","source_rms","source_max",
            "epsilon_min","epsilon_max","epsilon_mean","charge_min",
            "charge_max","charge_mean","charge_nonzero_cells"
        ])
    if all_rpbe_samples:
        write_csv(work / "pbe_oracle_rpbe_samples.csv", all_rpbe_samples, [
            "fixture","salt_M","point_kind","rank","i","j","k",
            "density_center","density_x","density_y","density_z",
            "eps_x","eps_y","eps_z","rho","source","source_status"
        ])
    if all_rpbe_top:
        write_csv(work / "pbe_oracle_rpbe_top_source.csv", all_rpbe_top, [
            "fixture","salt_M","rank","i","j","k","source","coulomb",
            "grad_coulomb_x","grad_coulomb_y","grad_coulomb_z",
            "grad_epsilon_x","grad_epsilon_y","grad_epsilon_z",
            "density_center","surface","epsilon_center","source_status"
        ])
    if all_focus_boundary:
        write_csv(work / "pbe_oracle_focus_boundary.csv", all_focus_boundary, [
            "fixture","salt_M","i","j","k","boundary_phi","field_source",
            "parent_grid","child_grid","interpolation_status"
        ])
    if all_focus_probe:
        write_csv(work / "pbe_oracle_focus_probe.csv", all_focus_probe, [
            "fixture","salt_M","probe_residue","probe_site","matched_atoms",
            "E_probe_kT","E_probe_kJmol","E_probe_kcalmol","status",
            "parent_grid","child_grid"
        ])
    if all_stage:
        stage_path = work / "barnase_stage_oracle.tsv"
        stage_rows_by_fixture = {}
        for row in read_tsv(stage_path):
            key = str(row.get("fixture", "")).strip().lower()
            if key:
                stage_rows_by_fixture[key] = row
        for row in all_stage:
            key = str(row.get("fixture", "")).strip().lower()
            if key:
                stage_rows_by_fixture[key] = row
        write_tsv(stage_path, list(stage_rows_by_fixture.values()), [
            "fixture","salt_M","midpoint_layout","midpoint_positive_face_order",
            "midpoint_full_export","grid_nx","grid_ny","grid_nz",
            "grid_origin_x_A","grid_origin_y_A","grid_origin_z_A",
            "grid_center_x_A","grid_center_y_A","grid_center_z_A",
            "grid_spacing_A","total_charge","charge_assigned","surface_points",
            "pydelphi_atom_count","pydelphi_threads",
            "pydelphi_radius_sum_A","pydelphi_radius_min_A",
            "pydelphi_radius_max_A","pydelphi_radius_mean_A",
            "surface_semantics","surface_zero_points","surface_one_points",
            "oracle_solute_points","oracle_solvent_points",
            "eps_sum_min","eps_sum_max","eps_sum_mean",
            "eps_mid_min","eps_mid_max","eps_mid_mean",
            "eps_mid_x_min","eps_mid_x_max","eps_mid_x_mean",
            "eps_mid_x_count","eps_mid_x_internal_count",
            "eps_mid_y_min","eps_mid_y_max","eps_mid_y_mean",
            "eps_mid_y_count","eps_mid_y_internal_count",
            "eps_mid_z_min","eps_mid_z_max","eps_mid_z_mean",
            "eps_mid_z_count","eps_mid_z_internal_count",
            "vdwms_boundary_initial","vdwms_boundary_final","vdwms_added",
            "vdwms_removed","vdwms_exposed_points","vdwms_accessible_atoms",
            "vdwms_external_boundary_points",
            "pydelphi_total_seconds","pydelphi_setup_seconds",
            "pydelphi_surface_seconds","pydelphi_epsilon_seconds",
            "pydelphi_solver_seconds","pydelphi_energy_seconds",
            "pydelphi_map_io_seconds","pydelphi_iterations",
            "final_E_rxn_kT","final_E_rxn_kJmol",
            "iscm_water_rf_kT","iscm_water_rf_kJmol",
            "iscm_vacuum_rf_kT","iscm_vacuum_rf_kJmol",
            "iscm_corrected_rxn_kT","iscm_corrected_rxn_kJmol",
            "iscm_finalization_mode","rtol","atol"
        ])
    write_csv(work / "pbe_delphi_export_report.csv", report_rows, [
        "fixture","salt_M","status","detail","charge_grid","surface_map",
        "midpoint_epsilon","boundary_map","operator_stage",
        "boundary_samples","field_samples","field_source",
        "field_precision_status","nwt_trace","nwt_solver_summary",
        "iscm_energy_oracle","gaussian_oracle","rpbe_oracle",
        "focus_probe_oracle","stage_oracle","surface_semantics"
    ])
    return status if args.strict else 0


if __name__ == "__main__":
    raise SystemExit(main())
'@
Set-Content -LiteralPath $driver -Value $driverCode -Encoding ASCII

try {
    $driverArgs = @($driver, '--pyroot', $PyDelPhiRoot, '--workroot', $WorkRoot, '--samples', $TitanSampleRoot, '--full-midpoint-salt', $FullMidpointSalt)
    if ($Strict) { $driverArgs += '--strict' }
    if ($FullMidpointEpsilon) { $driverArgs += '--full-midpoints' }
    if ($FullVdwmsState) { $driverArgs += '--full-vdwms-state' }
    if ($FullBoundaryOracle) { $driverArgs += '--full-boundary-oracle' }
    if ($FullOperatorOracle -or $FullChargeGrid) { $driverArgs += '--full-operator-oracle' }
    if ($FullChargeGrid) { $driverArgs += '--full-charge-grid' }
    if ($FullFieldOracle) { $driverArgs += '--full-field-oracle' }
    if ($FullIscmEnergyOracle) { $driverArgs += '--full-iscm-energy-oracle' }
    if ($FullGaussianOracle) { $driverArgs += '--full-gaussian-oracle' }
    if ($FullRpbOracle) { $driverArgs += '--full-rpb-oracle' }
    if ($FullFocusProbeOracle) { $driverArgs += '--full-focus-probe-oracle' }
    if ($FullFieldMap) { $driverArgs += '--full-field-map' }
    if ($NwtOnly) { $driverArgs += '--nwt-only' }
    foreach ($name in $Fixture) {
        if (-not [string]::IsNullOrWhiteSpace($name)) {
            $driverArgs += @('--fixture', $name)
        }
    }
    & $PythonExe @driverArgs | Out-Host
    $exit = $LASTEXITCODE
} catch {
    Write-ReportOnly -Status 'python_launch_failed' -Detail $_.Exception.Message
    return
}

$report = Join-Path $WorkRoot 'pbe_delphi_export_report.csv'
if (-not (Test-Path -LiteralPath $report)) {
    Write-ReportOnly -Status 'export_driver_failed' -Detail "driver did not write $report"
    return
}

$rows = @(Import-Csv -LiteralPath $report)
$failed = @($rows | Where-Object { $_.status -ne 'ok' })
if ($Strict -and ($exit -ne 0 -or $failed.Count -gt 0)) {
    throw "Strict DelPhi map export failed for $($failed.Count) row(s); see $report"
}
if ($failed.Count -gt 0) {
    Write-Warning "DelPhi map export completed report-only with $($failed.Count) failed row(s); see $report"
} else {
    Write-Host "PBE DelPhi map export PASS: $report"
}

$stage = Join-Path $WorkRoot 'barnase_stage_oracle.tsv'
$stageRowsForManifest = @()
if (Test-Path -LiteralPath $stage) {
    $stageRowsForManifest = @(Import-Csv -LiteralPath $stage -Delimiter "`t")
}
if ($stageRowsForManifest.Count -eq 0) {
    $stageRowsForManifest = @($rows)
}
$commandClassParts = New-Object System.Collections.Generic.List[string]
if ($FullMidpointEpsilon) { $commandClassParts.Add('full_midpoint') }
if ($FullVdwmsState) { $commandClassParts.Add('full_vdwms') }
if ($FullBoundaryOracle) { $commandClassParts.Add('full_boundary') }
if ($FullOperatorOracle) { $commandClassParts.Add('full_operator') }
if ($FullChargeGrid) { $commandClassParts.Add('full_charge_grid') }
if ($FullFieldOracle) { $commandClassParts.Add('full_field') }
if ($FullIscmEnergyOracle) { $commandClassParts.Add('full_iscm_energy') }
if ($FullGaussianOracle) { $commandClassParts.Add('full_gaussian') }
if ($FullRpbOracle) { $commandClassParts.Add('full_rpb') }
if ($FullFocusProbeOracle) { $commandClassParts.Add('full_focus_probe') }
if ($FullFieldMap) { $commandClassParts.Add('full_field_map') }
if ($NwtOnly) { $commandClassParts.Add('nwt_only') }
if ($commandClassParts.Count -eq 0) { $commandClassParts.Add('stage_summary') }
$commandClass = $commandClassParts -join '+'
$manifestOutputHash = (Get-HashText $report) + ';stage=' + (Get-HashText $stage)
foreach ($manifestRow in $stageRowsForManifest) {
    $fixtureName = Get-RowValue $manifestRow 'fixture' 'NA'
    $threadCount = Get-RowValue $manifestRow 'pydelphi_threads' 'NA'
    $statusValue = Get-RowValue $manifestRow 'status' 'ok'
    Add-OracleManifestRow @{
        timestamp_utc = (Get-Date).ToUniversalTime().ToString('o')
        producer_script = 'test/pbe_delphi_export_maps.ps1'
        command_class = $commandClass
        fixture = $fixtureName
        oracle_engine = 'pyDelPhi'
        python_exe = $PythonExe
        pyroot = $PyDelPhiRoot
        pyroot_git_commit = Get-PyDelPhiCommit $PyDelPhiRoot
        workroot = $WorkRoot
        input_hash = Get-HashText $GoldenCsv
        output_hash = $manifestOutputHash
        output_path = if (Test-Path -LiteralPath $stage) { $stage } else { $report }
        thread_count = $threadCount
        status = $statusValue
        note = 'black_box_oracle_output'
    }
}

if (Test-Path -LiteralPath $stage) {
    Write-Host "PBE DelPhi stage oracle: $stage"
} else {
    Write-Warning "No barnase_stage_oracle.tsv was exported; strict barnase stays disabled."
}
