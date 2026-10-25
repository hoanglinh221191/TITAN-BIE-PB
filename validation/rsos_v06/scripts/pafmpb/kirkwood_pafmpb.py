"""PAFMPB on the screened Kirkwood spheres, on the TITAN protocol meshes written by
kirkwood_surface/kirkwood_surface.py, against the exact energy, beside TITAN on the same nodes.
The exact solution makes this the reference-independent part of the PAFMPB comparison.

The mesh goes in as an OFF file (meshfmt 2): PAFMPB computes its own normals, and no coupling is
dropped.  In MSMS layout (meshfmt 1, analytic normals supplied) PAFMPB zeroes every coupling between
nodes closer than 0.4 A, which on these charge-graded meshes (edges down to 0.1 A) keeps the energy
from converging; `msms` reruns that layout for the record.

    py -3 kirkwood_pafmpb.py <afmpb.exe> [off|msms]   -> kirkwood_pafmpb_<layout>.csv
"""
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path('C:/Users/ADMIN/OneDrive/Documents/antigravity/elec_vdw')
sys.path.insert(0, str(REPO / 'test'))
sys.path.insert(0, str(Path(__file__).parent))
import tabi_node_patch_kirkwood as kw  # noqa: E402
from make_meshes import INPUT, INPUT_OFF, read_vtk  # noqa: E402

KS = Path('D:/data_MD/tabi/rsos_v06/kirkwood_surface')
RUNS = KS / 'runs'
ROOT = Path(__file__).parent


def titan_rows():
    """TITAN protocol energies on the same meshes, keyed by (radius, offset fraction, h)."""
    return {(float(r['radius_A']), float(r['offset_fraction']), float(r['h_A'])): r
            for r in csv.DictReader(open(KS / 'kirkwood_surface_paper_protocol.csv'))}


def main():
    exe = Path(sys.argv[1])
    layout = sys.argv[2] if len(sys.argv) > 2 else 'off'
    out_dir = ROOT / f'kirkwood_{layout}'
    out_dir.mkdir(exist_ok=True)
    T = titan_rows()
    rows = []
    for d in sorted(RUNS.glob('R*_t*_h*')):
        m = re.match(r'R(\d+)_t([0-9p]+)_h([0-9p]+)', d.name)
        radius, t, h = float(m.group(1)), float(m.group(2).replace('p', '.')), float(m.group(3).replace('p', '.'))
        vtk = d / 'pbe_tabi_surface_solution.vtk'
        if not vtk.exists() or (radius, t, h) not in T:
            continue
        pts, tris, normals = read_vtk(vtk)
        tr = T[(radius, t, h)]
        assert int(tr['nodes']) == len(pts), (d.name, tr['nodes'], len(pts))   # the same nodes
        run = out_dir / d.name
        if run.exists():
            shutil.rmtree(run)
        run.mkdir()
        tag = 'k'
        shutil.copyfile(RUNS / f'kirkwood_R{radius:g}_t{t:g}.pqr', run / f'{tag}.pqr')
        if layout == 'off':
            (run / 'input.txt').write_text(INPUT_OFF.format(tag=tag), newline='\n')
            with open(run / f'{tag}.off', 'w', newline='\n') as f:
                f.write(f'OFF\n{len(pts)} {len(tris)} 0\n')
                for x, y, z in pts:
                    f.write(f'{x:.6f} {y:.6f} {z:.6f}\n')
                for a, b, c in tris:
                    f.write(f'3 {a} {b} {c}\n')
        else:
            (run / 'input.txt').write_text(INPUT.format(tag=tag), newline='\n')
            with open(run / f'{tag}.msh', 'w', newline='\n') as f:
                f.write(f'{len(pts)} {len(tris)}\n')
                for (x, y, z), (nx, ny, nz) in zip(pts, normals):
                    f.write(f'{x:.6f} {y:.6f} {z:.6f} {nx:.6f} {ny:.6f} {nz:.6f}\n')
                for a, b, c in tris:
                    f.write(f'{a + 1} {b + 1} {c + 1}\n')
        subprocess.run([str(exe), '-i', 'input.txt', '-o', 'output.txt', '-s', 's.dat', '-t', '4'], cwd=run,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        (run / 's.dat').unlink(missing_ok=True)
        out = (run / 'output.txt').read_text(errors='replace')
        e = float(re.findall(r'\.\.\.polar part:\s+([-+0-9.eE]+)', out)[-1])
        its = re.findall(r'Iteration\s+(\d+) residual norm', out)
        exact = kw.screened_kirkwood(t, radius)
        et = float(tr['energy_kcalmol'])
        assert abs(float(tr['exact_energy_kcalmol']) - exact) < 1e-9 * abs(exact)
        rows.append(dict(radius_A=radius, offset_fraction=t, h_A=h, nodes=len(pts), exact_kcalmol=exact,
                         titan_kcalmol=et, pafmpb_kcalmol=e, titan_rel_pct=100 * (et - exact) / abs(exact),
                         pafmpb_rel_pct=100 * (e - exact) / abs(exact), pafmpb_iterations=int(its[-1]) if its else -1,
                         layout=layout))
        r = rows[-1]
        print(f"R={radius:<4g} t={t:<4g} h={h:<5g} nodes={len(pts):<6d} TITAN {r['titan_rel_pct']:+.3f}%  "
              f"PAFMPB {r['pafmpb_rel_pct']:+.3f}%", flush=True)
    with open(ROOT / f'kirkwood_pafmpb_{layout}.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == '__main__':
    main()
