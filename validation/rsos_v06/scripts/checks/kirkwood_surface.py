"""Boundary potential on the screened Kirkwood sphere against the exact series.

One sphere of radius R (dielectric 1 inside, 80 outside, 0.15 M) carries a
unit charge at t R from its centre (the PQR of test/tabi_node_patch_kirkwood.py).
For one sphere the SES is the sphere, so the protocol's boundary potential
phi at each node can be compared with the linearised-PB solution on r = R:

  phi(R, theta) = (C q / R) sum_l t^l (2l+1) / (eps_in l - eps_out gamma_l) P_l(cos theta)

with C = 332.063713299 kcal A/(mol e^2), gamma_l = l - z k_{l+1}(z)/k_l(z),
z = kappa R (TITAN's kappa), and theta measured from the charge axis.  The
continuity of phi and of eps dphi/dr at r = R gives the coefficients; at
kappa = 0, gamma_l = -(l+1) and the series is Kirkwood's.

    py -3 kirkwood_surface.py <titan.exe> [out_dir] [variant: paper_protocol | paper_protocol_p4]
"""
import csv
import math
import sys
from pathlib import Path

REPO = Path('C:/Users/ADMIN/OneDrive/Documents/antigravity/elec_vdw')
sys.path.insert(0, str(REPO / 'test'))
import tabi_node_patch_kirkwood as kw  # noqa: E402
import tabi_rsos_runner as runner  # noqa: E402

COULOMB = 332.063713299
CASES = [(10.0, t, h) for t in (0.0, 0.5, 0.9, 0.95) for h in (1.25, 1.0, 0.75, 0.5)] + \
        [(2.0, t, h) for t in (0.0, 0.5, 0.8) for h in (1.25, 1.0, 0.75, 0.5, 0.35, 0.25)]
VARIANT = sys.argv[3] if len(sys.argv) > 3 else 'paper_protocol'
if VARIANT != 'paper_protocol':      # the p4 control: the R = 10 A sphere, where p1's far field shows
    CASES = [c for c in CASES if c[0] == 10.0]


def coefficients(t, radius, epsin=1.0, epsout=80.0, salt=0.15, lmax=4000):
    """c_l with phi(R, theta) = sum_l c_l P_l(cos theta), in kcal/(mol e)."""
    kappa = 0.3291 * math.sqrt(salt) * math.sqrt(78.5 / epsout) if salt > 0 else 0.0
    z = kappa * radius
    out = []
    r = 1.0 + 1.0 / z if z > 0 else None          # k_1/k_0
    power = 1.0
    for l in range(lmax):
        gamma = (l - z * r) if z > 0 else -(l + 1.0)
        out.append(COULOMB / radius * power * (2 * l + 1) / (epsin * l - epsout * gamma))
        power *= t
        if z > 0:
            r = (2 * l + 3) / z + 1.0 / r
        if power < 1e-18:
            break
    return out


def exact(cosines, coeff):
    """Sum c_l P_l(x) by the three-term recurrence, for each x."""
    vals = []
    for x in cosines:
        p0, p1 = 1.0, x
        s = coeff[0] * p0 + (coeff[1] * p1 if len(coeff) > 1 else 0.0)
        for l in range(1, len(coeff) - 1):
            p0, p1 = p1, ((2 * l + 1) * x * p1 - l * p0) / (l + 1)
            s += coeff[l + 1] * p1
        vals.append(s)
    return vals


def main():
    titan = Path(sys.argv[1])
    root = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / 'runs'
    root.mkdir(parents=True, exist_ok=True)
    rows = []
    for radius, t, h in CASES:
        pqr = root / f'kirkwood_R{radius:g}_t{t:g}.pqr'
        kw.write_pqr(pqr, radius, t)
        tag = f'R{radius:g}_t{t:g}_h{h:g}'.replace('.', 'p')
        extra = kw.TIGHT + list(runner.VARIANTS[VARIANT][0]) + \
            ['-pbe-salt=0.15', '-pbe-bie-surface-output=all']
        res = runner.run_case(titan, root / tag, pqr, h, extra)
        pts = list(csv.DictReader(open(root / tag / 'pbe_tabi_surface_solution.csv')))
        coeff = coefficients(t, radius)
        cos = []
        for p in pts:
            x, y, zc = float(p['x_A']), float(p['y_A']), float(p['z_A'])
            cos.append(max(-1.0, min(1.0, x / math.sqrt(x * x + y * y + zc * zc))))
        ref = exact(cos, coeff)
        got = [float(p['total_potential_kcal_mol_e']) for p in pts]
        area = [float(p['area_A2']) for p in pts]
        err = [g - e for g, e in zip(got, ref)]
        rms = math.sqrt(sum(a * d * d for a, d in zip(area, err)) / sum(a * e * e for a, e in zip(area, ref)))
        mx = max(abs(d) for d in err) / max(abs(e) for e in ref)
        radial = max(abs(math.sqrt(float(p['x_A']) ** 2 + float(p['y_A']) ** 2 + float(p['z_A']) ** 2) - radius)
                     for p in pts)
        row = dict(radius_A=radius, offset_fraction=t, h_A=h, nodes=len(pts), iterations=res['iterations'],
                   energy_kcalmol=res['energy_kcalmol'], exact_energy_kcalmol=kw.screened_kirkwood(t, radius),
                   phi_exact_min=min(ref), phi_exact_max=max(ref), phi_rel_rms_error=rms, phi_rel_max_error=mx,
                   max_node_radius_offset_A=radial, status=res['status'])
        rows.append(row)
        print(f"R={radius:<4g} t={t:<4g} h={h:<5g} n={len(pts):<6d} it={res['iterations']:<3d} "
              f"phi range {min(ref):8.3f}..{max(ref):8.3f}  rms {rms:.2e}  max {mx:.2e}  "
              f"|r-R|max {radial:.1e}  E rel {(res['energy_kcalmol'] - row['exact_energy_kcalmol']) / abs(row['exact_energy_kcalmol']):+.2e}",
              flush=True)
        with open(root.parent / f'kirkwood_surface_{VARIANT}.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == '__main__':
    main()
