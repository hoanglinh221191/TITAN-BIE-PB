"""Internal review F4: the p4 spacing sweeps used as the protein reference (hsweep_p4, hsweep_lowsalt_p4)
were audited with 128 rows by an older binary; campaign_verify_p4 reruns them on the release binary with the
solution audit (512 targets, both equations).  Check that every case agrees in energy to 1e-6 relative,
converged, and carries the new audit; with --apply, park the old sweeps and put the new ones in their place,
so that compare_idle, presets_curve_idle and the manuscript builders read the re-audited reference.

    py -3 swap_p4_reference.py            # check only
    py -3 swap_p4_reference.py --apply    # check, then swap
"""
import csv
import shutil
import sys
from pathlib import Path

OUT = Path('D:/data_MD/tabi/rsos_v06')
MAIN = OUT / 'campaign_main_idle'
NEW = OUT / 'campaign_verify_p4'
PARK = MAIN / '_campaign_binary_timing'
STAGES = ('hsweep_p4', 'hsweep_lowsalt_p4')


def rows(path):
    return {r['case']: r for r in csv.DictReader(open(path))}


def check(stage):
    old, new = rows(MAIN / stage / f'{stage}.csv'), rows(NEW / stage / f'{stage}.csv')
    bad = []
    assert set(old) == set(new), (stage, sorted(set(old) ^ set(new)))
    for case, o in sorted(old.items()):
        n = new[case]
        eo, en = float(o['energy_kcalmol']), float(n['energy_kcalmol'])
        rel = abs(en - eo) / abs(eo)
        audit = int(float(n.get('audit_rows') or 0))
        ok = n['status'] == 'ok' and rel <= 1e-6 and audit >= min(1024, 2 * int(float(n['panels'])))
        print(f'{stage:18s} {case:24s} dE/E={rel:.1e} audit_rows {o.get("audit_rows")}->{audit} '
              f'e_op={n.get("op_audit_rel")} [{n["status"]}] {"ok" if ok else "CHECK"}')
        if not ok:
            bad.append(case)
    return bad


def main():
    bad = {s: check(s) for s in STAGES}
    if any(bad.values()):
        sys.exit(f'not swapped: {bad}')
    if '--apply' not in sys.argv:
        print('all cases agree; rerun with --apply to swap')
        return
    for s in STAGES:
        dest = PARK / f'{s}_audit128'
        assert not dest.exists(), dest
        shutil.move(str(MAIN / s), str(dest))
        shutil.copytree(NEW / s, MAIN / s)
        print('swapped', s, '(old in', dest, ')')
    print('next: rerun compare_idle and presets_curve_idle (run_night_v06_timing.sh lines after "afmpb bench"), '
          'then the manuscript builders')


if __name__ == '__main__':
    main()
