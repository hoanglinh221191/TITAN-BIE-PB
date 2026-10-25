#!/usr/bin/env python3
"""APBS grid sweep beside the TITAN BIE-PB convergence study (RSOS revision).

For each PQR and fine-grid spacing, writes the APBS 3.4.1 input of the
manuscript protocol (mg-auto, lpbe, bcfl sdh, chgm spl2, sdens 10, srad 1.4,
pdie/sdie 1/80, 0.15 M, 298.15 K; fine box = radius-expanded PQR bounds plus
12 A per side, coarse box 16 A per side larger; solvated minus a matched
sdie = 1, ion-free reference), runs it with four OpenMP threads, and records
the energy, wall time and peak working set.

Physics modes:
  paper    srfm smol, swin 0.3, ions of radius 2.0 A - the manuscript setting.
  matched  srfm mol (sharp SES dielectric boundary), ions of radius 0.  The
           ion-accessible region then starts at the atoms, inside the SES
           crevices too, while the boundary-integral model screens exactly
           outside the SES; APBS has no setting that places ions there.
  lowsalt  srfm mol, ions of radius 0 at 1e-5 M (Debye length about 960 A).
           Where the ions start no longer matters, so APBS and the
           boundary-integral model at the same ionic strength solve one
           boundary-value problem (sharp SES, eps 1/80) and converge to one
           value.  1e-5 M rather than none because the manuscript's FMM-near
           ILU(0) works on the coupled, salted system.
  vdw_lowsalt  as lowsalt with a zero probe radius, so srfm mol is the van
           der Waals surface (the union of atomic spheres), the surface the
           BIE sphere-patch discretization represents.

    python test/tabi_rsos_apbs_sweep.py --pqr a.pqr --spacing 0.6,0.5,0.4 \
        --mode matched --work-root D:/data_MD/tabi/rsos_v06/apbs
"""

import argparse
import csv
import ctypes
import math
import os
import re
import subprocess
import sys
import time
from ctypes import wintypes
from pathlib import Path

DEFAULT_APBS = Path(os.environ.get('USERPROFILE', '')) / \
    '.conda/envs/titan_apbs/Library/bin/apbs.exe'
KJ_PER_KCAL = 4.184


def read_atoms(pqr):
    atoms = []
    for line in Path(pqr).read_text().splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            fields = line[30:].split()
            x, y, z, q, r = (float(v) for v in fields[:5])
            atoms.append((x, y, z, q, r))
    return atoms


def box(atoms):
    lo = [min(a[k] - a[4] for a in atoms) for k in range(3)]
    hi = [max(a[k] + a[4] for a in atoms) for k in range(3)]
    fglen = [hi[k] - lo[k] + 24.0 for k in range(3)]
    cglen = [f + 32.0 for f in fglen]
    return fglen, cglen


def dime_for(length, spacing):
    points = max(65, math.ceil(length / spacing) + 1)
    return 32 * math.ceil((points - 1) / 32) + 1


def write_input(path, pqr_name, fglen, cglen, spacing, mode):
    dime = [dime_for(l, spacing) for l in fglen]
    if mode == 'paper':
        surface, ion_radius = 'srfm smol\n    swin 0.3', 2.0
    else:
        surface, ion_radius = 'srfm mol', 0.0
    conc = 1e-5 if mode in ('lowsalt', 'vdw_lowsalt') else 0.15
    # vdw_lowsalt: a zero probe makes srfm mol the van der Waals surface, the
    # union of atomic spheres that the BIE sphere-patch surface represents.
    srad = 0.0 if mode == 'vdw_lowsalt' else 1.4

    def elec(name, sdie, ions):
        ion_lines = (f'    ion charge 1 conc {conc} radius {ion_radius}\n'
                     f'    ion charge -1 conc {conc} radius {ion_radius}\n') if ions else ''
        return (f'elec name {name}\n    mg-auto\n'
                f'    dime {dime[0]} {dime[1]} {dime[2]}\n'
                f'    cglen {cglen[0]:.3f} {cglen[1]:.3f} {cglen[2]:.3f}\n'
                f'    fglen {fglen[0]:.3f} {fglen[1]:.3f} {fglen[2]:.3f}\n'
                '    cgcent mol 1\n    fgcent mol 1\n    mol 1\n    lpbe\n    bcfl sdh\n'
                f'{ion_lines}'
                f'    pdie 1.0\n    sdie {sdie}\n    {surface}\n    chgm spl2\n'
                f'    sdens 10.0\n    srad {srad}\n    temp 298.15\n'
                '    calcenergy total\n    calcforce no\nend\n')

    text = (f'read\n    mol pqr {pqr_name}\nend\n' + elec('solvated', 80.0, True)
            + elec('reference', 1.0, False)
            + 'print elecEnergy solvated - reference end\nquit\n')
    path.write_text(text)
    return dime


class ProcessMemoryCounters(ctypes.Structure):
    _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD),
                ('PeakWorkingSetSize', ctypes.c_size_t),
                ('WorkingSetSize', ctypes.c_size_t),
                ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                ('PagefileUsage', ctypes.c_size_t),
                ('PeakPagefileUsage', ctypes.c_size_t)]


def peak_working_set_mb(proc):
    """Peak working set of a finished child, read from its process handle."""
    if os.name != 'nt':
        return math.nan
    counters = ProcessMemoryCounters()
    counters.cb = ctypes.sizeof(counters)
    ok = ctypes.windll.psapi.GetProcessMemoryInfo(
        wintypes.HANDLE(int(proc._handle)), ctypes.byref(counters), counters.cb)
    return counters.PeakWorkingSetSize / 1048576.0 if ok else math.nan


def run_case(apbs, run_dir, pqr, spacing, mode, threads, timeout):
    run_dir.mkdir(parents=True, exist_ok=True)
    local_pqr = run_dir / Path(pqr).name
    local_pqr.write_bytes(Path(pqr).read_bytes())
    fglen, cglen = box(read_atoms(pqr))
    dime = write_input(run_dir / 'apbs.in', local_pqr.name, fglen, cglen,
                       spacing, mode)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    start = time.perf_counter()
    with open(run_dir / 'stdout.log', 'w') as out, \
            open(run_dir / 'stderr.log', 'w') as err:
        proc = subprocess.Popen([str(apbs), 'apbs.in'], cwd=run_dir, env=env,
                                stdout=out, stderr=err)
        try:
            code = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            code = 'timeout'
        wall = time.perf_counter() - start
        peak = peak_working_set_mb(proc)
    text = (run_dir / 'stdout.log').read_text(errors='replace')
    found = re.findall(r'Global net ELEC energy = ([-+0-9.Ee]+) kJ/mol', text)
    energy = float(found[-1]) / KJ_PER_KCAL if found else math.nan
    return {'dime': 'x'.join(map(str, dime)),
            'fine_spacing_A': max(f / (d - 1) for f, d in zip(fglen, dime)),
            'energy_kcalmol': energy, 'exit': code, 'wall_s': wall,
            'peak_rss_mb': peak, 'run_dir': run_dir.as_posix()}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--apbs', type=Path, default=DEFAULT_APBS)
    parser.add_argument('--pqr', required=True, help='comma-separated PQR files')
    parser.add_argument('--spacing', default='0.6,0.5,0.4,0.33,0.25')
    parser.add_argument('--mode', choices=['paper', 'matched', 'lowsalt', 'vdw_lowsalt'],
                        default='matched')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--timeout', type=float, default=7200)
    parser.add_argument('--work-root', type=Path, required=True)
    args = parser.parse_args()
    args.work_root.mkdir(parents=True, exist_ok=True)
    out = args.work_root / f'apbs_{args.mode}.csv'
    fields = ['case', 'mode', 'spacing_A', 'fine_spacing_A', 'dime',
              'energy_kcalmol', 'exit', 'wall_s', 'peak_rss_mb', 'run_dir']
    rows = []
    for pqr in args.pqr.split(','):
        name = Path(pqr).stem
        for spacing in (float(v) for v in args.spacing.split(',')):
            tag = f'{name}_{args.mode}_g{spacing:g}'.replace('.', 'p')
            row = {'case': name, 'mode': args.mode, 'spacing_A': spacing}
            row.update(run_case(args.apbs, args.work_root / tag, pqr, spacing,
                                args.mode, args.threads, args.timeout))
            rows.append(row)
            print(f"{name:24s} {args.mode:8s} g={spacing:<5g} dime={row['dime']:<13s} "
                  f"E={row['energy_kcalmol']:12.4f} {row['wall_s']:7.1f}s "
                  f"{row['peak_rss_mb']:8.0f}MB exit={row['exit']}", flush=True)
            with open(out, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
                writer.writeheader()
                writer.writerows(rows)
    print(f'wrote {out}')


if __name__ == '__main__':
    sys.exit(main())
