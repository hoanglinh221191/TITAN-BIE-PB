"""TITAN against PAFMPB (Reviewer 3, major comment 1), two of the 27 proteins.

TITAN runs the protocol (its own analytic SES with charge-graded refinement and curved node
patches).  PAFMPB runs (a) on the same nodes, read from the protocol mesh that make_meshes.py wrote
(as an OFF file: PAFMPB computes the normals itself), and (b) on quasi-uniform analytic SES meshes without
the refinement, like the MSMS meshes it was built for, at several spacings.  All runs use eps 1/80,
0.15 M, 298.15 K and four threads, and each energy is compared with TITAN p4 at h = 0.5 A.  The
configurations of a protein run in a rotated order in each repeat, so that every program sees the
same machine state.  PAFMPB's wall time and peak working set are measured from outside, as for
APBS; TITAN's come from its own records.  TITAN's time also contains building the surface, which
PAFMPB receives ready-made, so the time without the surface build is reported too.  TITAN's protocol
system is also solved without a preconditioner and with diagonal scaling in place of its sparse ILU(0),
in the same rotation, since PAFMPB's GMRES runs without one.  Both programs also report the time of the
GMRES loop alone, after the operator is set up.

    py -3 afmpb_bench.py <titan.exe> <afmpb.exe> [repeats] [out.csv]
"""
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path('C:/Users/ADMIN/OneDrive/Documents/antigravity/elec_vdw')
sys.path.insert(0, str(REPO / 'test'))
import tabi_rsos_apbs_sweep as apbs_driver  # noqa: E402  (peak_working_set_mb)
import tabi_rsos_campaign as c  # noqa: E402
import tabi_rsos_runner as runner  # noqa: E402
import make_meshes as meshes  # noqa: E402  (pqr_of: the 27-protein and 6LZG PQRs)

ROOT = Path(__file__).parent
NAMES = ['1sh1', '1a63']
UNIFORM = ['1p25', '0p75', '0p5', '0p35']
REF = Path('D:/data_MD/tabi/rsos_v06/campaign_main_idle/hsweep_p4/hsweep_p4.csv')
# TITAN's protocol system also solved without a preconditioner and with diagonal scaling, for comparison
# with PAFMPB, whose GMRES runs unpreconditioned on equations scaled to equal diagonal blocks
TITAN_PRECOND = ['none', 'diagonal']


def run_pafmpb(exe, src, run_dir, threads=4, timeout=3600):
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    # OFF input (meshfmt 2): PAFMPB computes the normals and node-patch geometry from the triangles.
    # With the analytic normals supplied in MSMS layout (meshfmt 1) its Kirkwood-sphere energies did
    # not converge with refinement (R = 2 A stayed about 5% low), while the OFF input converged.
    shutil.copyfile(src / 'input_off.txt', run_dir / 'input.txt')
    for f in src.iterdir():
        if f.suffix in ('.pqr', '.off'):
            shutil.copyfile(f, run_dir / f.name)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    start = time.perf_counter()
    with open(run_dir / 'run.log', 'w') as log:
        proc = subprocess.Popen([str(exe), '-i', 'input.txt', '-o', 'output.txt', '-s', 'surfp.dat',
                                 '-t', str(threads)], cwd=run_dir, env=env, stdout=log, stderr=subprocess.STDOUT)
        try:
            code = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            code = 'timeout'
        wall = time.perf_counter() - start
        peak = apbs_driver.peak_working_set_mb(proc)
    out = (run_dir / 'output.txt').read_text(errors='replace') if (run_dir / 'output.txt').exists() else ''

    def grab(label):
        m = re.findall(re.escape(label) + r'\s+([-+0-9.eE]+)', out)
        return float(m[-1]) if m else math.nan
    its = re.findall(r'Iteration\s+(\d+) residual norm', out)
    status = re.findall(r'GMRES solver exits with status:\s+(-?\d+)', out)
    (run_dir / 'surfp.dat').unlink(missing_ok=True)
    # '...time for GMRES' is the Windows build's wall time of the whole GMRES loop (omp_get_wtime); the
    # program's own matrix-vector total sums coarse per-call gettimeofday() differences
    return dict(exit=code, wall_s=wall, peak_rss_mb=peak, energy_kcalmol=grab('...polar part:'),
                iterations=int(its[-1]) if its else -1, solve_s=grab('...time for SolvePoissonBoltzmann'),
                gmres_s=grab('...time for GMRES'), internal_total_s=grab('... Total time'),
                status='ok' if status and status[-1] == '0' else 'failed')


def run_titan(titan, name, run_dir, precond=None):
    """The protocol; precond 'none' or 'diagonal' replaces its sparse ILU(0) (the later flag wins)."""
    extra = c.ANALYTIC + ([f'-pbe-bie-precond={precond}'] if precond else [])
    row = runner.run_case(titan, run_dir, meshes.pqr_of(name), 1.25, extra)
    t = list(csv.DictReader(open(run_dir / 'pbe_timing.csv')))[-1]
    num = lambda k: float(t.get(k) or 'nan')
    return dict(exit=row['exit'], wall_s=row['wall_s'], peak_rss_mb=row['peak_rss_mb'],
                energy_kcalmol=row['energy_kcalmol'], iterations=int(float(t['iterations'])),
                surface_build_s=num('tabi_surface_build_seconds'), gmres_s=num('tabi_gmres_seconds'),
                status=row['status'], op_audit_rel=row['op_audit_rel'],
                true_residual_rel=num('tabi_true_residual_rel'), precond_setup_s=num('tabi_preconditioner_setup_seconds'),
                precond_effective=t.get('tabi_preconditioner_effective', ''))


def main():
    titan, exe = Path(sys.argv[1]), Path(sys.argv[2])
    repeats = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    out_csv = Path(sys.argv[4]) if len(sys.argv) > 4 else ROOT / 'afmpb_bench.csv'
    ref = {r['case']: float(r['energy_kcalmol']) for r in csv.DictReader(open(REF))}
    rows = []

    def record(name, prog, setting, mesh, k, e_ref):
        tag = f"{name}_{prog.lower()}_{(mesh or 'protocol')}_r{k}"
        d = ROOT / 'bench' / tag
        if prog == 'TITAN':
            r = run_titan(titan, name, d, mesh)
        else:
            src = ROOT / 'pafmpb' / mesh
            r = run_pafmpb(exe, src, d)
            r['nodes'] = int(open(src / f'{mesh}.msh').readline().split()[0])
        r.update(protein=name.upper(), program=prog, setting=setting, repeat=k, reference_kcalmol=e_ref,
                 rel_error_pct=100 * (r['energy_kcalmol'] - e_ref) / abs(e_ref), run_dir=d.as_posix())
        rows.append(r)
        print(f"{name} r{k} {prog:6s} {setting:16s} E={r['energy_kcalmol']:.3f} ({r['rel_error_pct']:+.3f}%) "
              f"it={r['iterations']} gmres={r['gmres_s']:.2f}s wall={r['wall_s']:.2f}s peak={r['peak_rss_mb']:.0f}MB "
              f"[{r['status']}]", flush=True)
        fields = sorted({f for row in rows for f in row})
        with open(out_csv, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    for name in NAMES:
        e_ref = ref[f'{name}_charmm_p4_h0p5']
        configs = [('TITAN', 'protocol', None), ('PAFMPB', 'protocol mesh', name)] + \
                  [('PAFMPB', f'uniform h={u.replace("p", ".")}', f'{name}_uni_h{u}') for u in UNIFORM] + \
                  [('TITAN', f'protocol, {p}', p) for p in TITAN_PRECOND]
        for k in range(repeats):
            rot = k % len(configs)
            for prog, setting, mesh in configs[rot:] + configs[:rot]:
                record(name, prog, setting, mesh, k, e_ref)


if __name__ == '__main__':
    main()
