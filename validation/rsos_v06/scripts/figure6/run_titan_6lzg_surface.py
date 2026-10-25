"""Figure 6: the 6LZG complex under the protocol, with the boundary solution
written on the surface (pbe_tabi_surface_solution.vtk) for the potential map.
    py -3 run_titan_6lzg_surface.py <titan.exe>"""
import sys
from pathlib import Path

REPO = Path('C:/Users/ADMIN/OneDrive/Documents/antigravity/elec_vdw')
sys.path.insert(0, str(REPO / 'test'))
import tabi_rsos_campaign as c  # noqa: E402
import tabi_rsos_runner as runner  # noqa: E402

titan = Path(sys.argv[1])
pqr = c.SIXLZG / '6lzg_complex_charmm_ph7.pqr'
out = Path(__file__).parent / 'titan'
row = runner.run_case(titan, out, pqr, 1.25, c.ANALYTIC + ['-pbe-bie-surface-output=all'])
print({k: row[k] for k in ('status', 'energy_kcalmol', 'iterations', 'op_audit_rel', 'panels', 'wall_s')})
