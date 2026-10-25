"""Internal review F1: does the low-salt APBS calculation solve the TITAN problem?
(A) cavity convention: TITAN fills solvent cavities not connected to the exterior ("external-only");
    APBS srfm mol gives any probe-accessible region the solvent dielectric.  Fill the enclosed
    probe-accessible regions with uncharged probe-sized pseudo-atoms and compare APBS energies.
(B) outer boundary: bcfl sdh on a coarse box 16 A beyond the fine box; rerun with 48 A.
Energies only (no timing), lowsalt mode, one spacing.

    py -3 apbs_cavity_box_check.py <out_dir> [spacing]
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, 'C:/Users/ADMIN/OneDrive/Documents/antigravity/elec_vdw/test')
import tabi_rsos_apbs_sweep as ap  # noqa: E402
import tabi_rsos_campaign as c  # noqa: E402

PROBE = 1.4


def cavity_points(atoms, grid=0.15):
    """Probe-centre positions in regions not connected to the exterior.  26-connectivity on a 0.15 A grid:
    on a 0.3 A grid with 6-connectivity, grooves joined to the exterior by channels narrower than the grid
    came out as cavities (1SH1: 3 regions against 1 here; 1A63: 19 against 3), and filling them changed the
    surface itself."""
    xyz = np.array([a[:3] for a in atoms])
    rad = np.array([a[4] for a in atoms]) + PROBE
    lo = xyz.min(0) - rad.max() - 2 * grid
    hi = xyz.max(0) + rad.max() + 2 * grid
    n = np.ceil((hi - lo) / grid).astype(int) + 1
    free = np.ones(n, dtype=bool)                       # probe-centre accessible
    axes = [lo[k] + grid * np.arange(n[k]) for k in range(3)]
    for p, r in zip(xyz, rad):
        i0 = np.maximum(np.floor((p - r - lo) / grid).astype(int), 0)
        i1 = np.minimum(np.ceil((p + r - lo) / grid).astype(int) + 1, n)
        gx, gy, gz = np.meshgrid(*(axes[k][i0[k]:i1[k]] for k in range(3)), indexing='ij')
        inside = (gx - p[0]) ** 2 + (gy - p[1]) ** 2 + (gz - p[2]) ** 2 < r * r
        free[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] &= ~inside
    lab, nlab = ndimage.label(free, structure=ndimage.generate_binary_structure(3, 3))
    edge = set(np.unique(np.concatenate([lab[0].ravel(), lab[-1].ravel(), lab[:, 0].ravel(), lab[:, -1].ravel(),
                                         lab[:, :, 0].ravel(), lab[:, :, -1].ravel()])))
    cav = [k for k in range(1, nlab + 1) if k not in edge]
    pts = [np.argwhere(lab == k) * grid + lo for k in cav]
    return pts


def write_filled(pqr, out, pts):
    lines = [ln for ln in Path(pqr).read_text().splitlines() if not ln.startswith('END')]
    k = 9000
    for region in pts:
        for x, y, z in region:
            k += 1
            # APBS splits PQR records on whitespace: the serial must not touch the record name
            lines.append(f'ATOM  {k:5d}  CAV CAV X 999    {x:9.3f} {y:9.3f} {z:9.3f} {0.0:8.4f} {PROBE:7.4f}')
    lines.append('END')
    out.write_text('\n'.join(lines) + '\n')


def main():
    out = Path(sys.argv[1]).resolve()
    spacing = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
    out.mkdir(parents=True, exist_ok=True)
    apbs = ap.DEFAULT_APBS
    rows = []
    box0 = ap.box
    for name in ('1sh1', '1a63'):
        pqr = c.charmm(name)
        atoms = ap.read_atoms(pqr)
        pts = cavity_points(atoms)
        filled = out / f'{name}_cavfilled.pqr'
        write_filled(pqr, filled, pts)
        print(name, 'cavities', len(pts), 'probe-centre points', [len(p) for p in pts], flush=True)
        for variant, src, extra in (('base', pqr, 0.0), ('cavity_filled', filled, 0.0), ('coarse_box_48A', pqr, 64.0)):
            ap.box = lambda a, extra=extra: (lambda f, cg: (f, [x + extra for x in cg]))(*box0(a))
            r = ap.run_case(apbs, out / f'{name}_{variant}', src, spacing, 'lowsalt', 4, 7200)
            r.update(protein=name.upper(), variant=variant, cavities=len(pts),
                     cavity_points=sum(len(p) for p in pts))
            rows.append(r)
            print(name, variant, f"{r['energy_kcalmol']:.4f}", r['dime'], f"{r['wall_s']:.0f}s", flush=True)
            with open(out / 'apbs_cavity_box_check.csv', 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
    ap.box = box0


if __name__ == '__main__':
    main()
