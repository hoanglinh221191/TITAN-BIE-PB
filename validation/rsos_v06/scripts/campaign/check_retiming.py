"""After a re-timing night: every stage's non-timing fields against the parked rows of the previous binary
(energy to 1e-9 relative, iterations, status), and the status of every row against what the paper expects
(controls that are meant to fail).  Prints one line per stage and every unexpected row.

    py -3 check_retiming.py
"""
import csv
from pathlib import Path

OUT = Path('D:/data_MD/tabi/rsos_v06')
MAIN = OUT / 'campaign_main_idle'
OLD = MAIN / '_campaign_binary_timing'
STAGES = ['proteins', 'binding', 'hsweep_lowsalt', 'hsweep', 'presets', 'presets_lowsalt', 'presets27', 'threads',
          'xgen', 'cap', 'spacing', 'coarse', 'matrix', 'control451c']


def expected_failure(stage, case):
    """Rows the paper reports as failing: diagonal scaling on Grid-SES, and every control other than ILU(0)
    on the sphere-patch stress surfaces (Table 4); nothing else."""
    if stage != 'matrix':
        return False
    return case.endswith('_grid_ses_diagonal') or ('_sphere_patch_' in case and 'ilu0' not in case)


def rows(path):
    return {r['case']: r for r in csv.DictReader(open(path))} if path.exists() else {}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float('nan')


def main():
    bad = []
    for st in STAGES:
        new, old = rows(MAIN / st / f'{st}.csv'), rows(OLD / st / f'{st}.csv')
        diffs, its, unexpected = [], 0, []
        for case, n in new.items():
            ok = n.get('status') == 'ok'
            if not ok and not expected_failure(st, case):
                unexpected.append(f"{case} [{n.get('status')}]")
            o = old.get(case)
            if o and o.get('status') == 'ok' and ok:
                eo, en = num(o['energy_kcalmol']), num(n['energy_kcalmol'])
                diffs.append(abs(en - eo) / abs(eo))
                its += o.get('iterations') != n.get('iterations')
        mx = max(diffs) if diffs else float('nan')
        print(f'{st:16s} rows {len(new):3d}  compared {len(diffs):3d}  max dE/E {mx:.1e}  iteration changes {its}  '
              f'unexpected status {len(unexpected)}')
        for u in unexpected:
            print('    ', u)
        if (diffs and mx > 1e-9) or its or unexpected:
            bad.append(st)
    print('ALL CONSISTENT' if not bad else f'CHECK: {bad}')


if __name__ == '__main__':
    main()
