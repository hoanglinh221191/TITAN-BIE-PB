"""Does a new TITAN binary give the same PBE results as the release binary on the paper's paths?
Cases: protocol p1 on 1SH1 and 1A63, p4 on 1SH1, the 1BRS binding cycle, a Kirkwood sphere.  Energies must
agree to 1e-9 relative (iterations exactly); timings are not compared.

    py -3 binary_equivalence.py <old.exe> <new.exe> <work_dir>
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, 'C:/Users/ADMIN/OneDrive/Documents/antigravity/elec_vdw/test')
import tabi_rsos_campaign as c  # noqa: E402
import tabi_rsos_runner as runner  # noqa: E402


def cases():
    yield '1sh1_protocol', c.charmm('1sh1'), 1.25, c.ANALYTIC
    yield '1a63_protocol', c.charmm('1a63'), 1.25, c.ANALYTIC
    yield '1sh1_p4', c.charmm('1sh1'), 1.25, c.ANALYTIC + c.P4
    bb = (c.BINDING_JPCB / '1brs_complex_amber_ph7.pqr', c.BINDING_JPCB / '1brs_barnase_chain_A_amber_ph7.pqr',
          c.BINDING_JPCB / '1brs_barstar_chain_D_amber_ph7.pqr')
    yield '1brs_binding', bb[0], 1.25, c.ANALYTIC + c.binding_args(*bb)


def main():
    old, new, work = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]).resolve()
    rows, bad = [], []
    for name, pqr, h, extra in cases():
        r = {}
        for tag, exe in (('old', old), ('new', new)):
            r[tag] = runner.run_case(exe, work / f'{name}_{tag}', pqr, h, extra)
        eo, en = r['old']['energy_kcalmol'], r['new']['energy_kcalmol']
        rel = abs(en - eo) / abs(eo)
        same = rel <= 1e-9 and r['old']['iterations'] == r['new']['iterations'] and r['new']['status'] == 'ok'
        print(f"{name:16s} old {eo:.9f} ({r['old']['iterations']} it)  new {en:.9f} ({r['new']['iterations']} it)  "
              f"dE/E {rel:.1e}  {'same' if same else 'DIFFERENT'}", flush=True)
        rows.append(dict(case=name, old_energy=eo, new_energy=en, rel=rel, old_it=r['old']['iterations'],
                         new_it=r['new']['iterations'], status=r['new']['status'], same=same))
        if not same:
            bad.append(name)
    with open(work / 'binary_equivalence.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    sys.exit(f'different: {bad}' if bad else 0)


if __name__ == '__main__':
    main()
