#!/usr/bin/env python3
"""Run one PBE example in <example>/run and print its main results.

  python run_example.py 01_1vii_solvation --titan ../../titan.exe
  python run_example.py 03_6lzg_binding_cycle --titan ../../titan

The example directory holds arg.in and its PQR files; common/ holds the neutral driver pair and the minimal
force-field library that every run needs (MANUAL_BIE_PB.md, section 4)."""
import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def last_row(path):
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f))
    return rows[-1] if rows else {}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('example', help='example directory name, e.g. 01_1vii_solvation')
    ap.add_argument('--titan', required=True, help='path to the TITAN executable')
    a = ap.parse_args()
    example = (HERE / a.example).resolve()
    titan = Path(a.titan).resolve()
    run = example / 'run'
    if run.exists():
        shutil.rmtree(run)
    run.mkdir()
    for name in ('driver.pdb', 'driver.top', 'index.ndx', 'cal.txt'):
        shutil.copy2(HERE / 'common' / name, run / name)
    shutil.copytree(HERE / 'common' / 'lib_ff', run / 'lib_ff')
    shutil.copy2(example / 'arg.in', run / 'arg.in')
    for pqr in example.glob('*.pqr'):
        shutil.copy2(pqr, run / pqr.name)
    with open(run / 'run.log', 'w') as log:
        code = subprocess.call([str(titan), '-a', 'arg.in'], cwd=run, stdout=log, stderr=subprocess.STDOUT)
    if code != 0 or not (run / 'pbe_energy.csv').exists():
        print(f'TITAN exited with code {code}; see {run / "run.log"}')
        return 1
    e, t = last_row(run / 'pbe_energy.csv'), last_row(run / 'pbe_timing.csv')
    print(f"node patches {e.get('solute_points')}, GMRES iterations {e.get('iterations_solvent')}, "
          f"true residual {t.get('tabi_true_residual_rel')}, e_op,S {t.get('tabi_tree_vs_direct_matvec_rel')}")
    print(f"dG_solv = {float(e['E_rxn_kcalmol']):.6f} kcal/mol")
    b = run / 'pbe_binding_electrostatics.csv'
    if b.exists():
        r = last_row(b)
        print('dG_pol,bind = {:.3f} kcal/mol (complex {:.3f}, receptor {:.3f}, ligand {:.3f})'.format(
            *(float(r[k]) / 4.184 for k in ('Delta_Gpol_bind_kJmol', 'Gpol_complex_kJmol', 'Gpol_receptor_kJmol',
                                             'Gpol_ligand_kJmol'))))
    print(f'outputs: {run}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
