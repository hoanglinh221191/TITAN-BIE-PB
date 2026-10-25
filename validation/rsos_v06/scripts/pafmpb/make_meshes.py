"""Shared-mesh benchmark against PAFMPB (Reviewer 3): the TITAN protocol surface of two of the 27
proteins, written as a PAFMPB mesh in MSMS layout (meshfmt 1: vertices with their analytic SES
normals, then 1-based triangles), so both programs solve on the same nodes.  Each mesh directory
also gets the PQR and a PAFMPB input file (eps 1/80, 150 mM, 298.15 K, no nonpolar terms), and the
same triangles as an OFF file with its input (meshfmt 2: PAFMPB computes the normals), which the
benchmark uses.  The 6LZG complex and its two chains (Table 6) are named 6lzg, 6lzg_ace2, 6lzg_rbd.

    py -3 make_meshes.py <titan.exe> [spacings, default 1.25] [uniform|protocol] [names, default 1sh1,1a63]
    (needs the ucrt64 runtime on PATH)
"""
import shutil
import sys
from pathlib import Path

REPO = Path('C:/Users/ADMIN/OneDrive/Documents/antigravity/elec_vdw')
sys.path.insert(0, str(REPO / 'test'))
import tabi_rsos_campaign as c  # noqa: E402
import tabi_rsos_runner as runner  # noqa: E402

ROOT = Path(__file__).parent
NAMES = ['1sh1', '1a63']
INPUT = ('#di, de, ion concentration(mM), temperature, surface tension, pressure\n'
         '1.0 80.0 150 298.15 0.0 0.0\n'
         '# meshfmt: 0 auto-mesh, 1 msms, 2 off\n'
         '1 0 0\n'
         '#pqr and mesh files\n'
         '{tag}.pqr {tag}.msh\n')
INPUT_OFF = ('#di de conc T gamma beta\n'
             '1.0 80.0 150 298.15 0.0 0.0\n'
             '# meshfmt: 2 off\n'
             '2 0 0\n'
             '#pqr and mesh files\n'
             '{tag}.pqr {tag}.off\n')
SIXLZG = {'6lzg': '6lzg_complex_charmm_ph7.pqr', '6lzg_ace2': '6lzg_ace2_chain_A_charmm_ph7.pqr',
          '6lzg_rbd': '6lzg_rbd_chain_B_charmm_ph7.pqr'}


def pqr_of(name):
    return c.SIXLZG / SIXLZG[name] if name in SIXLZG else c.charmm(name)


def read_vtk(path):
    lines = path.read_text().splitlines()
    pts, tris, normals = [], [], []
    i = 0
    while i < len(lines):
        f = lines[i].split()
        if f[:1] == ['POINTS']:
            n = int(f[1])
            pts = [tuple(map(float, r.split()[:3])) for r in lines[i + 1:i + 1 + n]]
            i += n + 1
            continue
        if f[:1] == ['POLYGONS']:
            m = int(f[1])
            assert all(len(r.split()) == 4 and r.split()[0] == '3' for r in lines[i + 1:i + 1 + m])
            tris = [tuple(map(int, r.split()[1:4])) for r in lines[i + 1:i + 1 + m]]
            i += m + 1
            continue
        if f[:1] == ['VECTORS'] or f[:1] == ['NORMALS']:
            normals = [tuple(map(float, r.split()[:3])) for r in lines[i + 1:i + 1 + len(pts)]]
            i += len(pts) + 1
            continue
        i += 1
    assert pts and tris and len(normals) == len(pts)
    return pts, tris, normals


def main():
    titan = Path(sys.argv[1])
    spacings = [float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1.25]
    # 'uniform': the analytic SES without the charge-graded refinement, a quasi-uniform mesh like MSMS's
    uniform = len(sys.argv) > 3 and sys.argv[3] == 'uniform'
    names = sys.argv[4].split(',') if len(sys.argv) > 4 else NAMES
    for name in names:
        for h in spacings:
            pqr = pqr_of(name)
            tag = name if h == 1.25 and not uniform else f"{name}{'_uni' if uniform else ''}_h{h:g}".replace('.', 'p')
            run = ROOT / 'titan_mesh' / tag
            surface = c.NODE_PATCH if uniform else c.ANALYTIC
            row = runner.run_case(titan, run, pqr, h, surface + ['-pbe-bie-surface-output=all'])
            pts, tris, normals = read_vtk(run / 'pbe_tabi_surface_solution.vtk')
            out = ROOT / 'pafmpb' / tag
            out.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(pqr, out / f'{tag}.pqr')
            (out / 'input.txt').write_text(INPUT.format(tag=tag), newline='\n')
            with open(out / f'{tag}.msh', 'w', newline='\n') as f:
                f.write(f'{len(pts)} {len(tris)}\n')
                for (x, y, z), (nx, ny, nz) in zip(pts, normals):
                    line = f'{x:.6f} {y:.6f} {z:.6f} {nx:.6f} {ny:.6f} {nz:.6f}\n'
                    assert len(line) < 99
                    f.write(line)
                for a, b, d in tris:
                    f.write(f'{a + 1} {b + 1} {d + 1}\n')
            (out / 'input_off.txt').write_text(INPUT_OFF.format(tag=tag), newline='\n')
            with open(out / f'{tag}.off', 'w', newline='\n') as f:
                f.write(f'OFF\n{len(pts)} {len(tris)} 0\n')
                for x, y, z in pts:
                    f.write(f'{x:.6f} {y:.6f} {z:.6f}\n')
                for a, b, d in tris:
                    f.write(f'3 {a} {b} {d}\n')
            print(tag, 'h', h, 'nodes', len(pts), 'triangles', len(tris), 'TITAN E', row['energy_kcalmol'],
                  'iterations', row['iterations'], flush=True)


if __name__ == '__main__':
    main()
