#!/usr/bin/env python3
"""Shared runner for the TITAN BIE-PB revision campaign (RSOS, v0.6).

One case = one PQR solved by titan.exe in its own directory under the paper
protocol plus a named variant.  run_case() writes arg.in, runs the solver and
returns one flat dict of the quantities the manuscript tables report, read
from the CSVs TITAN writes.  The Kirkwood harness and the protein matrices
both build on it, so every number in the revision comes through one path.

As a script it runs a protein matrix:

    python test/tabi_rsos_runner.py --pqr a.pqr,b.pqr --h 1.25,1.0 \
        --variants paper_np_curved --work-root D:/data_MD/tabi/rsos_v06/x
"""

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIB_FF = ROOT / '.codex_tmp/g4c2_package_preflight/inputs/lib_ff'

DRIVER_PDB = """CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   VAL A   1      -2.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  VAL A   1      -0.600   0.100   0.000  1.00  0.00           C
ATOM      3  C   VAL A   1       0.000   1.450   0.000  1.00  0.00           C
ATOM      4  O   VAL A   1      -0.600   2.430   0.000  1.00  0.00           O
END
"""

DRIVER_TOP = """[ defaults ]
1 2 yes 0.5 0.8333333333
[ atomtypes ]
T 12.0 0.0 A 0.3 0.1
[ moleculetype ]
TEST 3
[ atoms ]
1 T 1 VAL N 1 0.0 12.0
2 T 1 VAL CA 2 0.0 12.0
3 T 1 VAL C 3 0.0 12.0
4 T 1 VAL O 4 0.0 12.0
[ system ]
Neutral four-atom CLI bootstrap; the PB source is the PQR
[ molecules ]
TEST 1
"""

# The paper protocol (manuscript section 3.5) with the SES engine chosen for
# the revision.  Variants override or extend it; a later option wins.
BASE = [
    '-pbe=bie',
    '-pbe-surface=vdw',
    '-pbe-bie-primary-surface=analytic-ses',
    '-pbe-bie-analytic-ses-engine=patch-ray-resolution',
    '-pbe-bie-sum=treecode',
    '-pbe-bie-formulation=jcp2013-treecode',
    '-pbe-bie-fast-operator=p1',
    '-pbe-bie-tree-leaf-size=64',
    '-pbe-bie-opening-theta=0.45',
    '-pbe-bie-near-integration=analytic_singular',
    '-pbe-bie-quad-buffer=off',
    '-pbe-bie-precond=fmm-nearfield-ilu0',
    '-pbe-bie-precond-ilu-neighbor-cap=64',
    '-pbe-bie-precond-ilu-drop-tol=0',
    '-pbe-bie-precond-ilu-regularization=adaptive',
    '-pbe-bie-precond-ilu-regularization-max=0.01',
    '-pbe-bie-precond-memory-budget-mb=4096',
    '-pbe-bie-gmres-restart=128',
    '-pbe-max-iter=720',
    '-pbe-rms=1e-4',
    '-pbe-dphi=1e-4',
    '-pbe-int-diel=1',
    '-pbe-ext-diel=80',
    '-pbe-salt=0.15',
    '-pbe-temp=298.15',
    '-pbe-probe=1.4',
    '-pbe-timing=detail',
    # The operator audit of manuscript section 2.6: 512 target panels at
    # index midpoints, both equations, on the solution.
    '-pbe-bie-matvec-audit=targets',
    '-pbe-bie-matvec-audit-sample-count=512',
]

# name -> (extra arguments, extra environment)
VARIANTS = {
    # What the paper protocol runs on main before the fix: the curved cell is
    # refused on the FMM route, so node-patch is the flat legacy cell.
    'paper_np_flat': (['-pbe-bie-discretization=node-patch'],
                      {'TITAN_CURVED_NODE_PATCH': 'off'}),
    'paper_centroid': (['-pbe-bie-discretization=centroid'], {}),
    # The curved node-patch as it exists on the classic treecode.
    'curved_treecode': (['-pbe-bie-discretization=node-patch',
                         '-pbe-bie-fmm-m2l=off',
                         '-pbe-bie-tree-builder=topdown-aabb',
                         '-pbe-bie-tree-leaf-size=8',
                         '-pbe-bie-leaf-mac=on',
                         '-pbe-bie-opening-theta=0.5',
                         '-pbe-bie-taylor-order=4',
                         '-pbe-bie-precond=classic-nearfield-ilu0',
                         '-pbe-bie-curved-sample-order=4',
                         '-pbe-bie-energy-far-eta=0'], {}),
    # The paper protocol with the curved cell on the FMM route.
    'paper_np_curved': (['-pbe-bie-discretization=node-patch'], {}),
    # The revised protocol: the curved cell plus charge-adaptive remeshing.
    'paper_protocol': (['-pbe-bie-discretization=node-patch',
                        '-pbe-bie-ses-charge-refine=0.35',
                        '-pbe-bie-ses-charge-refine-min-edge=0.1'], {}),
    # The revised protocol under the guarded p4 far field.
    'paper_protocol_p4': (['-pbe-bie-discretization=node-patch',
                           '-pbe-bie-ses-charge-refine=0.35',
                           '-pbe-bie-ses-charge-refine-min-edge=0.1',
                           '-pbe-bie-fmm-local=gmres-p4',
                           '-pbe-bie-fmm-local-refine=operator-tail-static',
                           '-pbe-bie-fmm-local-refine-tol=0.002',
                           '-pbe-bie-fmm-local-cache=green-jet-compressed',
                           '-pbe-bie-fmm-local-cache-budget-mb=256',
                           '-pbe-bie-fmm-local-fallback=packed-avx2',
                           '-pbe-bie-fmm-local-plan-builder=metric-cache',
                           '-pbe-bie-fmm-local-refine-probe-count=4',
                           '-pbe-bie-fmm-local-refine-probe-mode=rademacher',
                           '-pbe-bie-fmm-local-geometry-guard=strict-geometry-mac',
                           '-pbe-bie-fmm-local-geometry-theta=0.5',
                           '-pbe-bie-taylor-order=4',
                           # the solution audit of every other row (section 2.6)
                           '-pbe-bie-matvec-audit=targets'], {}),
    # The guarded degree-four accuracy control (manuscript Table S11): the
    # same operator with the local p4 far field.
    'paper_np_curved_p4': (['-pbe-bie-discretization=node-patch',
                            '-pbe-bie-fmm-local=gmres-p4',
                            '-pbe-bie-fmm-local-refine=operator-tail-static',
                            '-pbe-bie-fmm-local-refine-tol=0.002',
                            '-pbe-bie-fmm-local-cache=green-jet-compressed',
                            '-pbe-bie-fmm-local-cache-budget-mb=256',
                            '-pbe-bie-fmm-local-fallback=packed-avx2',
                            '-pbe-bie-fmm-local-plan-builder=metric-cache',
                            '-pbe-bie-fmm-local-refine-probe-count=4',
                            '-pbe-bie-fmm-local-refine-probe-mode=rademacher',
                            '-pbe-bie-fmm-local-geometry-guard=strict-geometry-mac',
                            '-pbe-bie-fmm-local-geometry-theta=0.5',
                            '-pbe-bie-taylor-order=4',
                            '-pbe-bie-matvec-audit=targets'], {}),
}

FIELDS = [
    'case', 'variant', 'h_A', 'status', 'exit', 'energy_kcalmol', 'panels',
    'iterations', 'true_residual_rel', 'op_audit_rel', 'direct_residual_rel',
    'precond_effective',
    'ilu_shift', 'precond_setup_s', 'precond_apply_s', 'gmres_s', 'total_s',
    'wall_s', 'peak_rss_mb', 'node_patch_representation',
    'node_patch_normal_source', 'near_kernel_policy', 'precond_requested',
    'precond_quality', 'audit_rows', 'binding_dgpol_kcalmol',
    'gpol_complex_kcalmol', 'gpol_receptor_kcalmol', 'gpol_ligand_kcalmol',
    'ses_vertices', 'ses_refined_vertices', 'ses_edge_min_A',
    'ses_min_angle_below_5_10_20', 'certificate_source',
    'certificate_components', 'certificate_complete',
    'pn_star_normal_vertices', 'pn_turned_normal_vertices', 'pn_repair_rounds',
    'pn_flattened_triangles', 'pn_separated_triangles', 'pn_straightened_edges',
    'pn_corner_vertices', 'pn_max_collocation_offset_A',
    'pn_flat_contacts', 'pn_max_shared_edge_gap_A', 'self_max_order',
    'self_fallback_blocks', 'self_fallback_step', 'q_rms', 'q_max_median',
    'q_checkerboard_pairs', 'near_adaptive_pairs',
    'near_unresolved_pairs', 'near_min_gap_A', 'near_unconverged_pairs',
    'near_unconverged_worst', 'progress_iter', 'progress_residual',
    'progress_elapsed_s', 'residual_history', 'run_dir',
]


def _first_row(path):
    if not path.exists():
        return {}
    with open(path, newline='', encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    return rows[0] if rows else {}


def _float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def base_for(extra):
    """BASE, without the analytic-SES engine when extra picks another surface
    (the CLI refuses the engine option there)."""
    other = any(o.startswith('-pbe-bie-primary-surface=')
                and o != '-pbe-bie-primary-surface=analytic-ses' for o in extra)
    return [o for o in BASE if not (other and (
        o.startswith('-pbe-bie-analytic-ses-engine=')
        or o.startswith('-pbe-bie-primary-surface=')))]


def run_case(titan, run_dir, pqr, h, extra=(), env_extra=None, nt=4,
             timeout=7200, keep=False):
    """Solve one PQR in run_dir and return the reported quantities."""
    titan = Path(titan).resolve()  # the solver runs with cwd=run_dir
    pqr = Path(pqr).resolve()
    run_dir = Path(run_dir)
    if run_dir.exists() and not keep:
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / 'driver.pdb').write_text(DRIVER_PDB, encoding='ascii')
    (run_dir / 'driver.top').write_text(DRIVER_TOP, encoding='ascii')
    lines = [
        f'-f={(run_dir / "driver.pdb").as_posix()}',
        '-1pdb=true',
        f'-common_topo={(run_dir / "driver.top").as_posix()}',
        '-ff=AMBER99SB',
        '-g=1&1',
        f'-path-a={LIB_FF.as_posix()}/',
        f'-path-n={LIB_FF.as_posix()}/',
        f'-pbe-test-pqr={Path(pqr).as_posix()}',
        f'-pbe-bie-ses-grid-spacing={h}',
        f'-nt={nt}',
    ] + base_for(extra) + list(extra)
    (run_dir / 'arg.in').write_text('\n'.join(lines) + '\n', encoding='ascii')
    env = dict(os.environ)
    env.update(env_extra or {})
    start = time.perf_counter()
    with open(run_dir / 'run.log', 'w') as log:
        try:
            code = subprocess.run([str(titan), '-a', 'arg.in'], cwd=run_dir,
                                  env=env, stdout=log, stderr=subprocess.STDOUT,
                                  timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            code = 'timeout'
    wall = time.perf_counter() - start
    (run_dir / 'case_meta.json').write_text(
        json.dumps({'exit': code, 'wall_s': wall}), encoding='ascii')
    return collect_case(run_dir)


CERTIFICATE_KEYS = {
    # The lift since the normal repair, and before it (flattened, separated).
    'pn_star_normals': 'pn_star_normal_vertices', 'pn_turned_normals': 'pn_turned_normal_vertices',
    'pn_repair_rounds': 'pn_repair_rounds', 'pn_residual_flat': 'pn_flattened_triangles',
    'pn_corner_vertices': 'pn_corner_vertices', 'pn_max_offset_A': 'pn_max_collocation_offset_A',
    'pn_flattened': 'pn_flattened_triangles', 'pn_separated': 'pn_separated_triangles',
    'pn_straightened_edges': 'pn_straightened_edges', 'pn_flat_contacts': 'pn_flat_contacts',
    'pn_edge_gap_A': 'pn_max_shared_edge_gap_A', 'near_adaptive': 'near_adaptive_pairs',
    'near_unresolved': 'near_unresolved_pairs', 'near_min_gap_A': 'near_min_gap_A',
    'self_max_order': 'self_max_order', 'self_fallback': 'self_fallback_blocks',
    'self_fallback_step': 'self_fallback_step',
    # The shape of the solved normal derivative (a spurious oscillating mode).
    'q_rms': 'q_rms', 'q_max_median': 'q_max_median',
    'q_checkerboard_pairs': 'q_checkerboard_pairs',
}
LARGEST = {'pn_edge_gap_A', 'pn_max_offset_A', 'self_max_order',
           'self_fallback_step', 'q_rms', 'q_max_median'}


def node_patch_certificate(run_dir, log):
    """The PN-lift and near-pair certificate of a run, every value None when
    the run did not record it (an older binary, or not a curved surface).

    Newer binaries write one profile row per frame, with separate certificate
    columns for complex, receptor and ligand. Counts are summed over every
    recorded component; gaps are maximised (edge) or minimised (near).
    Older ones only logged: the flattened-triangle count, and near pairs that
    kept an unconverged order-256 value with the worst last change."""
    out = {name: None for name in CERTIFICATE_KEYS.values()}
    out.update(certificate_source='none', certificate_components='',
               certificate_complete=False, near_unconverged_pairs=None,
               near_unconverged_worst=None)
    rows = []
    profile = Path(run_dir) / 'pbe_tabi_profile.csv'
    if profile.exists():
        with open(profile, newline='', encoding='utf-8', errors='replace') as f:
            rows = list(csv.DictReader(f))
    columns = {'complex': 'node_patch_certificate',
               'receptor': 'node_patch_certificate_receptor',
               'ligand': 'node_patch_certificate_ligand'}
    binding = (Path(run_dir) / 'pbe_binding_electrostatics.csv').exists()
    expected = set(columns) if binding else {'complex'}
    records, components = [], set()
    complete = bool(rows)
    for row in rows:
        present = set()
        for component, column in columns.items():
            value = row.get(column, '')
            if value and value != 'n/a':
                records.append(value)
                present.add(component)
        components.update(present)
        complete = complete and expected.issubset(present)
    if records:
        out['certificate_source'] = 'profile'
        out['certificate_components'] = ';'.join(c for c in columns if c in components)
        out['certificate_complete'] = complete
        for row in records:
            for item in row.split(';'):
                key, _, value = item.partition('=')
                name = CERTIFICATE_KEYS.get(key)
                if name is None or value in ('', 'n/a'):
                    continue
                number = float(value)
                old = out[name]
                if old is None:
                    out[name] = number
                elif key in LARGEST:
                    out[name] = max(old, number)
                elif key == 'near_min_gap_A':
                    out[name] = min(old, number)
                else:
                    out[name] = old + number
        return out
    flattened = re.findall(r'TABI curved PN: (\d+) of', log)
    if flattened:
        out['certificate_source'] = 'log'
        out['pn_flattened_triangles'] = sum(int(v) for v in flattened)
        kept = re.findall(r'TABI curved near: (\d+) near-pair evaluation\(s\) took', log)
        out['near_unconverged_pairs'] = sum(int(v) for v in kept)
        worst = [float(v) for v in re.findall(r'worst last relative change (\S+)', log)]
        out['near_unconverged_worst'] = max(worst) if worst else None
    return out


def collect_case(run_dir):
    """The quantities of a finished case, read back from its run directory."""
    run_dir = Path(run_dir)
    meta = json.loads((run_dir / 'case_meta.json').read_text(encoding='ascii'))
    code, wall = meta['exit'], meta['wall_s']
    energy = _first_row(run_dir / 'pbe_energy.csv')
    timing = _first_row(run_dir / 'pbe_timing.csv')
    profile = _first_row(run_dir / 'pbe_tabi_profile.csv')
    memory = _first_row(run_dir / 'pbe_memory_profile.csv')
    peak = _float(memory.get('peak_working_set_bytes'))
    binding = _first_row(run_dir / 'pbe_binding_electrostatics.csv')
    kcal = lambda key: _float(binding.get(key)) / 4.184
    log = (run_dir / 'run.log').read_text(encoding='utf-8', errors='replace') \
        if (run_dir / 'run.log').exists() else ''

    def logged(pattern, default=''):
        found = re.findall(pattern, log)
        return found[-1] if found else default

    certificate = node_patch_certificate(run_dir, log)

    # Where GMRES stood when the run ended, which for a run stopped at its
    # time budget is the only record of how far it got.
    last_iter = {}
    progress = run_dir / 'pbe_tabi_longrun_progress.csv'
    if progress.exists():
        with open(progress, newline='', encoding='utf-8', errors='replace') as f:
            for row in csv.DictReader(f):
                if row.get('phase') == 'gmres_iter':
                    last_iter = row

    return {
        'progress_iter': last_iter.get('gmres_iter', ''),
        'progress_residual': last_iter.get('residual', ''),
        'progress_elapsed_s': last_iter.get('elapsed_seconds', ''),
        'ses_vertices': logged(r'patch-ray-resolution: h=\S+ A vertices=(\d+)'),
        'ses_refined_vertices': logged(r'refined_vertices=(\d+)'),
        'ses_edge_min_A': logged(r'edge_min_A=(\S+)'),
        'ses_min_angle_below_5_10_20': logged(r'min_angle_below_5/10/20_deg=(\S+)'),
        **certificate,
        'precond_requested': profile.get('preconditioner_requested', ''),
        'precond_quality': profile.get('preconditioner_quality_status', ''),
        'audit_rows': profile.get('matvec_audit_rows', ''),
        # r_dir,S of the manuscript (with -pbe-bie-matvec-audit=targets).
        'direct_residual_rel': _float(
            profile.get('fmm_local_sampled_direct_true_residual')),
        'binding_dgpol_kcalmol': kcal('Delta_Gpol_bind_kJmol'),
        'gpol_complex_kcalmol': kcal('Gpol_complex_kJmol'),
        'gpol_receptor_kcalmol': kcal('Gpol_receptor_kJmol'),
        'gpol_ligand_kcalmol': kcal('Gpol_ligand_kJmol'),
        'residual_history': profile.get('residual_history', ''),
        'status': energy.get('status', 'no_output'),
        'exit': code,
        'energy_kcalmol': _float(energy.get('E_rxn_kcalmol')),
        'panels': int(_float(energy.get('solute_points')) or 0)
                  if energy.get('solute_points') else 0,
        'iterations': int(_float(energy.get('iterations_solvent')) or 0)
                      if energy.get('iterations_solvent') else 0,
        'true_residual_rel': _float(timing.get('tabi_true_residual_rel')),
        'op_audit_rel': _float(timing.get('tabi_tree_vs_direct_matvec_rel')),
        'precond_effective': profile.get('preconditioner_effective', ''),
        'ilu_shift': profile.get('preconditioner_ilu_regularization_selected', ''),
        'precond_setup_s': _float(profile.get('preconditioner_setup_seconds')),
        'precond_apply_s': _float(profile.get('preconditioner_apply_seconds')),
        'gmres_s': _float(profile.get('gmres_seconds')),
        'total_s': _float(profile.get('total_seconds')),
        'wall_s': wall,
        'peak_rss_mb': peak / 1048576.0 if math.isfinite(peak) else math.nan,
        'node_patch_representation': profile.get('node_patch_representation', ''),
        'node_patch_normal_source': profile.get('node_patch_normal_source', ''),
        'near_kernel_policy': profile.get('near_kernel_policy', ''),
        'run_dir': run_dir.as_posix(),
    }


def write_rows(path, rows, fields=FIELDS):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f'{v:.10g}' if isinstance(v, float) else v)
                             for k, v in row.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--titan', type=Path, default=ROOT / 'titan.exe')
    parser.add_argument('--pqr', required=True,
                        help='comma-separated PQR files; the case name is the stem')
    parser.add_argument('--h', default='1.25')
    parser.add_argument('--variants', default='paper_np_curved')
    parser.add_argument('--extra', default='',
                        help='semicolon-separated options appended to every case')
    parser.add_argument('--work-root', type=Path, required=True)
    parser.add_argument('--nt', type=int, default=4)
    parser.add_argument('--timeout', type=float, default=7200)
    args = parser.parse_args()
    extra_all = [o for o in args.extra.split(';') if o]
    args.work_root.mkdir(parents=True, exist_ok=True)
    out = args.work_root / 'matrix.csv'
    rows = []
    for pqr in args.pqr.split(','):
        name = Path(pqr).stem
        for variant in args.variants.split(','):
            extra, env = VARIANTS[variant]
            for h in (float(v) for v in args.h.split(',')):
                tag = f'{name}_{variant}_h{h:g}'.replace('.', 'p')
                row = {'case': name, 'variant': variant, 'h_A': h}
                row.update(run_case(args.titan, args.work_root / tag, pqr, h,
                                    list(extra) + extra_all, env, args.nt,
                                    args.timeout))
                rows.append(row)
                print(f"{name:24s} {variant:16s} h={h:<5g} n={row['panels']:<7d} "
                      f"it={row['iterations']:<4d} E={row['energy_kcalmol']:12.4f} "
                      f"eop={row['op_audit_rel']:.3e} rss={row['peak_rss_mb']:.0f}MB "
                      f"{row['wall_s']:.1f}s [{row['status']}|{row['exit']}] "
                      f"{row['node_patch_representation']} {row['node_patch_normal_source']}",
                      flush=True)
                write_rows(out, rows)
    print(f'wrote {out}')


if __name__ == '__main__':
    sys.exit(main())
