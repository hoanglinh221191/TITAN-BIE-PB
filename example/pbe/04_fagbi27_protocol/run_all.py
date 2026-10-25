#!/usr/bin/env python3
"""Run the paper's protocol on its 27 proteins and compare with the paper's results.

  python run_all.py --titan ../../../titan.exe                 # all 27, four threads
  python run_all.py --titan ../../../titan.exe --only 1SH1,1A63

Each protein runs in run/<PDB id>/ from protocol_arg.in.template, with the neutral four-atom driver pair and the
minimal force-field library of ../common (the PB source is the PQR, -pbe-test-pqr). The results go to
run/results.csv and are compared with reference_protocol.csv (Figure 5 and Supplementary Table S7 of the
paper): the reaction-field energy, the GMRES iterations and the node-patch count.  With the release source the
energies agree to the printed digits; other compilers or CPUs may differ in the last digits."""
import argparse
import csv
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMMON = HERE.parent / 'common'
GATE = 8e-3      # sampled operator discrepancy gate of the paper (Section 2.6)


def last_row(path):
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f))
    return rows[-1] if rows else {}


def run_one(titan, pdb_id, pqr, threads, root):
    d = root / pdb_id
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    for name in ('driver.pdb', 'driver.top', 'index.ndx', 'cal.txt'):
        shutil.copy2(COMMON / name, d / name)
    shutil.copytree(COMMON / 'lib_ff', d / 'lib_ff')
    shutil.copy2(HERE / 'inputs' / pqr, d / pqr)
    args = (HERE / 'protocol_arg.in.template').read_text().replace('{PQR}', pqr)
    args = '\n'.join(f'-nt={threads}' if l.startswith('-nt=') else l for l in args.splitlines()) + '\n'
    (d / 'arg.in').write_text(args)
    t0 = time.perf_counter()
    with open(d / 'run.log', 'w') as log:
        code = subprocess.call([str(titan), '-a', 'arg.in'], cwd=d, stdout=log, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - t0
    energy = last_row(d / 'pbe_energy.csv') if (d / 'pbe_energy.csv').exists() else {}
    timing = last_row(d / 'pbe_timing.csv') if (d / 'pbe_timing.csv').exists() else {}
    return dict(pdb_id=pdb_id, exit=code, status=energy.get('status', 'no_output'),
                node_patches=energy.get('solute_points', ''), iterations=energy.get('iterations_solvent', ''),
                dGsolv_kcalmol=energy.get('E_rxn_kcalmol', ''),
                true_residual_rel=timing.get('tabi_true_residual_rel', ''),
                e_op_S=timing.get('tabi_tree_vs_direct_matvec_rel', ''), wall_s=f'{wall:.1f}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--titan', required=True, help='path to the TITAN executable')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--only', default='', help='comma-separated PDB ids')
    ap.add_argument('--out', default=str(HERE / 'run'))
    a = ap.parse_args()
    titan = Path(a.titan).resolve()
    root = Path(a.out).resolve()
    root.mkdir(parents=True, exist_ok=True)
    manifest = list(csv.DictReader(open(HERE / 'inputs' / 'fagbi27_fixture_manifest.csv')))
    only = {x.strip().upper() for x in a.only.split(',') if x.strip()}
    ref = {r['pdb_id']: r for r in csv.DictReader(open(HERE / 'reference_protocol.csv'))}
    results = []
    print(f"{'PDB':5s} {'patches':>8s} {'it':>3s} {'dGsolv (kcal/mol)':>18s} {'paper':>14s} {'diff (%)':>9s} "
          f"{'e_op,S':>9s} {'gate':>5s} {'wall (s)':>8s}")
    for r in manifest:
        if only and r['pdb_id'] not in only:
            continue
        res = run_one(titan, r['pdb_id'], r['pqr_file'], a.threads, root)
        results.append(res)
        p = ref[r['pdb_id']]
        try:
            e, pe = float(res['dGsolv_kcalmol']), float(p['dGsolv_kcalmol'])
            diff = f'{(e - pe) / abs(pe) * 100:+.1e}'
            gate = 'pass' if float(res['e_op_S']) <= GATE else 'fail'
        except ValueError:
            e, pe, diff, gate = float('nan'), float(p['dGsolv_kcalmol']), 'n/a', 'n/a'
        print(f"{r['pdb_id']:5s} {res['node_patches']:>8s} {res['iterations']:>3s} {e:18.6f} {pe:14.6f} "
              f"{diff:>9s} {res['e_op_S'][:9]:>9s} {gate:>5s} {res['wall_s']:>8s}", flush=True)
    with open(root / 'results.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(results[0]))
        w.writeheader()
        w.writerows(results)
    print(f'results: {root / "results.csv"}')
    print('451C is reported with its operator gate failure at h = 1.25 A, as in the paper (Section 3.4).')
    return 0 if all(x['exit'] == 0 for x in results) else 1


if __name__ == '__main__':
    sys.exit(main())
