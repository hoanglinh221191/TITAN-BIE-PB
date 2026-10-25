#!/usr/bin/env python3
"""Prepare a protonated PQR from a PDB file the way the paper prepared its 27 proteins.

Step 1 keeps the first model and its ATOM records (protein chains; HETATM records are dropped).
  --policy as-deposited  keeps the ATOM records as deposited (26 of the 27 proteins);
  --policy clean         also removes input hydrogens and keeps one conformer per residue, ranked by
                         occupancy, then altloc A (used for 1CBN, whose alternate conformers PDB2PQR rejects).
Step 2 runs PDB2PQR 3.7.1, which adds hydrogens, assigns protonation states with PROPKA at pH 7.0, and writes
CHARMM charges and radii:
  pdb2pqr --ff=CHARMM --ffout=CHARMM --keep-chain --titration-state-method=propka --with-ph=7.0 in.pdb out.pqr
The PQR is the input of TITAN (-pbe-test-pqr); TITAN reads charges and radii from it.

  python prepare_pqr.py 1sh1_source.pdb 1sh1.pqr                     # policy as-deposited
  python prepare_pqr.py 1cbn_source.pdb 1cbn.pqr --policy clean
  python prepare_pqr.py --verify                                     # rebuild the 27 shipped PQRs and compare

PDB2PQR must be importable by the Python that runs this script (pip install pdb2pqr==3.7.1)."""
import argparse
import csv
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PDB2PQR = ['--ff=CHARMM', '--ffout=CHARMM', '--keep-chain', '--titration-state-method=propka', '--with-ph=7.0']


def first_model_atoms(lines):
    has_models = any(l.startswith('MODEL') for l in lines)
    in_first, seen = not has_models, False
    for line in lines:
        if line.startswith('MODEL'):
            in_first, seen = not seen, True
            continue
        if line.startswith('ENDMDL'):
            if in_first:
                return
            continue
        if in_first and line.startswith('ATOM'):
            yield line.rstrip('\n')


def clean(rows):
    """Drop hydrogens and keep one conformer per residue (highest occupancy, then altloc A)."""
    atoms = []
    for i, line in enumerate(rows):
        fixed = line.ljust(80)
        name, element = fixed[12:16].strip(), fixed[76:78].strip().upper()
        if element in ('H', 'D') or name.lstrip('0123456789')[:1] in ('H', 'D'):
            continue
        try:
            occ = float(fixed[54:60])
        except ValueError:
            occ = 0.0
        atoms.append(dict(line=fixed[:16] + ' ' + fixed[17:], res=(fixed[21], fixed[22:26], fixed[26]),
                          resname=fixed[17:20], name=name, alt=fixed[16].strip(), occ=occ, i=i))
    choice = {}
    for res in {a['res'] for a in atoms}:
        alts = [a for a in atoms if a['res'] == res and a['alt']]
        if not alts:
            continue
        keys = {}
        for a in alts:
            k = (a['alt'], a['resname'])
            keys[k] = max(keys.get(k, -1.0), a['occ'])
        choice[res] = sorted(keys, key=lambda k: (-keys[k], 0 if k[0] == 'A' else 1, k))[0]
    best = {}
    for a in atoms:
        if a['alt'] and (a['alt'], a['resname']) != choice[a['res']]:
            continue
        rank = 2 if not a['alt'] else (1 if a['alt'] == 'A' else 0)
        key = (a['res'], a['name'])
        old = best.get(key)
        if old is None or a['occ'] > old['occ'] or (a['occ'] == old['occ'] and rank > old['rank']):
            best[key] = dict(a, rank=rank)
    return [a['line'].rstrip() for a in sorted(best.values(), key=lambda a: a['i'])]


def prepare(pdb, pqr, policy):
    rows = list(first_model_atoms(Path(pdb).read_text().splitlines()))
    if policy == 'clean':
        rows = clean(rows)
    if not rows:
        raise SystemExit(f'no model-1 ATOM records in {pdb}')
    with tempfile.TemporaryDirectory() as tmp:
        prepared = Path(tmp) / 'prepared.pdb'
        prepared.write_text('\n'.join(rows + ['TER', 'END']) + '\n')
        log = Path(str(pqr) + '.pdb2pqr.log')
        with open(log, 'w') as f:
            code = subprocess.call([sys.executable, '-m', 'pdb2pqr', *PDB2PQR, str(prepared), str(pqr)],
                                   stdout=f, stderr=subprocess.STDOUT)
    if code != 0:
        raise SystemExit(f'PDB2PQR failed for {pdb}; see {log}')
    normalize(pqr)


def normalize(pqr, floor=0.6):
    """As in the paper: a non-positive radius becomes 0.6 A (none of the 27 proteins has one), and the file is
    written with CRLF line ends and a final line end, so the bytes match on every platform."""
    out = []
    for line in Path(pqr).read_text().splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            fields = line.split()
            if float(fields[-1]) <= 0.0:
                line = line[:line.rfind(fields[-1])] + f'{floor:.4f}'
        out.append(line)
    Path(pqr).write_bytes(('\r\n'.join(out) + '\r\n').encode('ascii'))


def verify():
    inputs = HERE / 'inputs'
    out = HERE / 'run' / 'prepared_pqr'
    out.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(inputs / 'fagbi27_fixture_manifest.csv')))
    same = 0
    for r in rows:
        policy = 'clean' if r['conformer_policy'] != 'as deposited' else 'as-deposited'
        pqr = out / r['pqr_file']
        prepare(inputs / r['raw_pdb_file'], pqr, policy)
        ok = hashlib.sha256(pqr.read_bytes()).hexdigest() == r['pqr_sha256']
        same += ok
        print(f"{r['pdb_id']:5s} {policy:13s} {'identical' if ok else 'DIFFERENT'}")
    print(f'{same} of {len(rows)} PQR files identical to the shipped ones')
    return 0 if same == len(rows) else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('pdb', nargs='?')
    ap.add_argument('pqr', nargs='?')
    ap.add_argument('--policy', choices=('as-deposited', 'clean'), default='as-deposited')
    ap.add_argument('--verify', action='store_true', help='rebuild the 27 shipped PQR files and compare them')
    a = ap.parse_args()
    if a.verify:
        return verify()
    if not (a.pdb and a.pqr):
        ap.error('give an input PDB and an output PQR, or --verify')
    prepare(a.pdb, a.pqr, a.policy)
    return 0


if __name__ == '__main__':
    sys.exit(main())
