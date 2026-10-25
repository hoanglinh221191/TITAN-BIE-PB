from pathlib import Path
from pymol import cmd
from pymol.cgo import BEGIN, COLOR, END, NORMAL, TRIANGLES, VERTEX

ROOT = Path.cwd()
if not (ROOT / '6lzg_titan_vs_apbs_surface.vtk').exists():
    ROOT = Path('D:\\data_MD\\tabi\\rsos_v06\\figure6\\common')
VTK = ROOT / '6lzg_titan_vs_apbs_surface.vtk'
PQR = ROOT / '6lzg_complex_charmm_ph7.pqr'
TITAN_LIMIT = 5.0
APBS_LIMIT = 5.0

def read_vtk(path):
    lines = path.read_text().splitlines()
    points, triangles, normals, fields = [], [], [], {}
    n = 0
    i = 0
    while i < len(lines):
        f = lines[i].split()
        if f[:1] == ['POINTS']:
            n = int(f[1]); points = [tuple(map(float, r.split()[:3])) for r in lines[i+1:i+1+n]]; i += n + 1; continue
        if f[:1] == ['POLYGONS']:
            m = int(f[1]); triangles = [tuple(map(int, r.split()[1:4])) for r in lines[i+1:i+1+m]]; i += m + 1; continue
        if f[:1] == ['SCALARS']:
            name = f[1]; fields[name] = [float(v) for v in lines[i+2:i+2+n]]; i += n + 2; continue
        if f[:1] == ['VECTORS']:
            normals = [tuple(map(float, r.split()[:3])) for r in lines[i+1:i+1+n]]; i += n + 1; continue
        i += 1
    return points, triangles, normals, fields

def bwr(value, limit):
    x = max(-1.0, min(1.0, value / limit))
    if x < 0.0:
        return (1.0, 1.0 + x, 1.0 + x)
    return (1.0 - x, 1.0 - x, 1.0)

def mesh_cgo(points, triangles, normals, values, shift, limit):
    obj = [BEGIN, TRIANGLES]
    for tri in triangles:
        for idx in tri:
            x, y, z = points[idx]
            nx, ny, nz = normals[idx]
            r, g, b = bwr(values[idx], limit)
            obj += [COLOR, r, g, b, NORMAL, nx, ny, nz, VERTEX, x + shift, y, z]
    obj.append(END)
    return obj

points, triangles, normals, fields = read_vtk(VTK)
span = max(p[0] for p in points) - min(p[0] for p in points)
offset = 0.60 * span + 8.0
cmd.load_cgo(mesh_cgo(points, triangles, normals, fields['titan_total_potential_kT_e'], -offset, TITAN_LIMIT), 'TITAN_analytic_SES')
cmd.load_cgo(mesh_cgo(points, triangles, normals, fields['apbs_total_potential_kT_e'], offset, APBS_LIMIT), 'APBS_on_TITAN_SES')
cmd.load(str(PQR), 'structure_left')
cmd.create('structure_right', 'structure_left')
cmd.translate([-offset, 0.0, 0.0], 'structure_left')
cmd.translate([offset, 0.0, 0.0], 'structure_right')
cmd.hide('everything', 'structure_left or structure_right')
cmd.show('cartoon', 'structure_left or structure_right')
cmd.color('gray70', 'structure_left or structure_right')
cmd.set('cgo_transparency', 0.0)
cmd.set('orthoscopic', 1)
cmd.set('ray_opaque_background', 0)
cmd.bg_color('white')
cmd.zoom('all', buffer=4.0)
