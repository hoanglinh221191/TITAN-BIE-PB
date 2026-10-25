"""Show the TITAN surface potential in PyMOL, coloured red (negative) - white - blue (positive) as in Figure 6 of
the paper.

  pymol view_surface_pymol.py -- run/pbe_tabi_surface_solution.vtk 1vii_full_charmm_ph7.pqr
  pymol -cq view_surface_pymol.py -- run/pbe_tabi_surface_solution.vtk 1vii_full_charmm_ph7.pqr 5 map.png

Arguments after --: the VTK file written with -pbe-bie-surface-output=all, an optional PQR or PDB file of the
structure (shown as cartoon), the colour limit in kT/e (default 5, the scale of Figure 6), and an optional PNG to
write. The surface is the TITAN mesh itself; the colour is total_potential_kT_e at the nodes, interpolated
linearly across each triangle."""
import sys
from pathlib import Path

from pymol import cmd
from pymol.cgo import BEGIN, COLOR, END, NORMAL, TRIANGLES, VERTEX


def read_vtk(path):
    lines = Path(path).read_text().splitlines()
    points, triangles, normals, fields = [], [], [], {}
    n = i = 0
    while i < len(lines):
        f = lines[i].split()
        if f[:1] == ['POINTS']:
            n = int(f[1])
            points = [tuple(map(float, r.split()[:3])) for r in lines[i + 1:i + 1 + n]]
            i += n + 1
            continue
        if f[:1] == ['POLYGONS']:
            m = int(f[1])
            triangles = [tuple(map(int, r.split()[1:4])) for r in lines[i + 1:i + 1 + m]]
            i += m + 1
            continue
        if f[:1] == ['SCALARS']:        # followed by a LOOKUP_TABLE line
            fields[f[1]] = [float(v) for v in lines[i + 2:i + 2 + n]]
            i += n + 2
            continue
        if f[:1] == ['VECTORS']:
            normals = [tuple(map(float, r.split()[:3])) for r in lines[i + 1:i + 1 + n]]
            i += n + 1
            continue
        i += 1
    return points, triangles, normals, fields


def bwr(value, limit):
    x = max(-1.0, min(1.0, value / limit))
    return (1.0, 1.0 + x, 1.0 + x) if x < 0.0 else (1.0 - x, 1.0 - x, 1.0)


def surface_cgo(points, triangles, normals, values, limit):
    obj = [BEGIN, TRIANGLES]
    for tri in triangles:
        for k in tri:
            obj += [COLOR, *bwr(values[k], limit), NORMAL, *normals[k], VERTEX, *points[k]]
    obj.append(END)
    return obj


args = sys.argv[1:]
if not args:
    raise SystemExit(__doc__)
vtk = args[0]
structure = args[1] if len(args) > 1 and args[1] else None
limit = float(args[2]) if len(args) > 2 else 5.0
png = args[3] if len(args) > 3 else None

points, triangles, normals, fields = read_vtk(vtk)
values = fields['total_potential_kT_e']
cmd.load_cgo(surface_cgo(points, triangles, normals, values, limit), 'titan_potential')
if structure:
    cmd.load(structure, 'structure')
    cmd.hide('everything', 'structure')
    cmd.show('cartoon', 'structure')
    cmd.color('gray70', 'structure')
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.zoom('titan_potential', complete=1)
print(f'{len(points)} nodes, {len(triangles)} triangles; potential {min(values):.2f} to {max(values):.2f} kT/e, '
      f'colour limit +/-{limit:g} kT/e (red negative, blue positive)')
if png:
    cmd.png(png, width=1600, height=1200, dpi=300, ray=1)
    print(f'wrote {png}')
