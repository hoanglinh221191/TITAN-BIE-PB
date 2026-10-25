#!/usr/bin/env python3
"""The TITAN BIE-PB revision campaign (RSOS): every manuscript table rerun.

Protocol (manuscript section 3.5, revised): analytic SES from the
patch-ray-resolution engine at h = 1.25 A with charge-adaptive remeshing
(alpha 0.35, smallest edge 0.1 A), curved node-patch cells, p1 Cartesian
multipole operator (theta 0.45, leaf 64), FMM-near ILU(0) cap 64, GMRES
restart 128, true residual 1e-4, four threads, 0.15 M.  Each stage keeps the
settings the submitted manuscript used for that table; only the surface
engine, the node-patch cell and the refinement differ.  Grid-SES rows use the
same curved node patches (no refinement: it lives in the analytic engine);
sphere-patch rows stay centroid, the only valid choice on an open surface.

Stages (manuscript item):
  proteins  Figure 5, Table S7: the 27 CHARMM proteins
  hsweep    spacing convergence (reviewer 2): five CHARMM proteins, h 1.5-0.5
  hsweep_lowsalt  the same at 1e-5 M, for the accuracy-matched APBS cost
  hsweep_p4, hsweep_lowsalt_p4  the same under the guarded p4 (converged reference)
  matrix    Tables 3-5, S8, S9, Figures 2-3: five AMBER geometries x three
            surfaces x six preconditioners
  binding   Table 6: 6LZG ACE2-RBD (standard, and the published cap-8 coarse
            / 1 GiB configuration) and Barnase-Barstar polar cycles
  p1p4      Table S11: p1 against the guarded p4 at theta 0.50
  spacing   Table S12: three AMBER geometries across h
  cap       Table S3: ILU neighbour cap 32/64/128, under p4 as published and p1
  coarse    Table S13: cap-8 coarse aggregation on 1A63
  s1        Table S1: 1UBQ, 1VII, 2CI2 (AMBER) under p1 and p4
  dense     Table S4: small whole proteins, dense (theta 0.05) against p1 and p4
  xgen      Figure 4A: 1A63, Grid-SES against analytic SES under one discretization
  threads   Figure 4B: 1A63, protocol and Grid-SES node-patch, one and four threads
  s2        Table S2: guarded p4 settings on the sphere-patch stress geometry
  s5        Table S5: 1VII spacing under p4 (protocol surface)
  s6        Table S6: direct, p1 + diagonal/block/ILU(0), p4 + ILU(0) on sphere-patch

    python -u test/tabi_rsos_campaign.py --resume \
        --titan D:/data_MD/tabi/rsos_v06/bin/titan_453f1e9f.exe \
        --work-root D:/data_MD/tabi/rsos_v06/campaign
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tabi_rsos_runner as runner  # noqa: E402

ROOT = runner.ROOT
TMP = ROOT / 'test/tmp'
PUBLIC = ROOT / 'test/pbe_inputs_public'
FAGBI27 = TMP / 'pbe_tabi_fagbi27_node_patch'
SIXLZG = TMP / 'pbe_tabi_6lzg_ace2_rbd_binding_apbs/6lzg_titan_binding'
BINDING_JPCB = ROOT / 'test/pbe_inputs_binding_jpcb'
MINI = TMP / 'fagbi_scope_z4y_mini_dense'
S1_INPUTS = TMP / 'fagbi_scope_s1_analytic_nodepatch_p4_apbs'

NODE_PATCH = ['-pbe-bie-discretization=node-patch']
REFINE = ['-pbe-bie-ses-charge-refine=0.35',
          '-pbe-bie-ses-charge-refine-min-edge=0.1']
ANALYTIC = NODE_PATCH + REFINE
# Grid-SES rows of the preconditioner matrix stay centroid, as submitted;
# node-patch on Grid-SES is the flat cell (the curved cell needs the
# remeshed analytic SES).
GRID = ['-pbe-bie-primary-surface=grid-ses', '-pbe-bie-discretization=centroid']
GRID_NP = ['-pbe-bie-primary-surface=grid-ses'] + NODE_PATCH
FLAT_CELL = {'TITAN_CURVED_NODE_PATCH': 'off'}
# The sphere-patch geometry of the solver tests: the clipped union of atomic
# spheres (the van der Waals surface) at the submitted 8 points per sphere.
# The coarse mesh is several percent from the converged energy but poses the
# same linear system: on 1VII the ILU(0) needs 34 iterations here and 38 at
# 96 points.  The submitted whole-triangle surface (-pbe-bie-sphere-mesh=
# latlong -pbe-bie-sphere-patch-clip=off) overlaps itself and has no physical
# energy.
SPHERE = ['-pbe-bie-primary-surface=sphere-patch',
          '-pbe-bie-surface-points=8',
          '-pbe-bie-discretization=centroid']
P4 = runner.VARIANTS['paper_np_curved_p4'][0][1:]  # without the discretization
# The dense reference: at theta 0.05 the tree accepts no far-field
# interaction, so every pair is summed directly (sampled audit ~1e-16 against
# the direct operator, p1_curved_operator_theta0.05 in
# pbe_bie_resolution_ses_smoke.py) while the ILU(0) still preconditions it.
# Near zone 16 adds the effect of the one-point rule beyond the zone.
DENSE = ['-pbe-bie-opening-theta=0.05']
DENSE16 = DENSE + ['-pbe-bie-curved-near-zone=16']
TIGHT = ['-pbe-rms=1e-6', '-pbe-dphi=1e-6', '-pbe-max-iter=1440']
COARSE8 = ['-pbe-bie-precond=fmm-nearfield-ilu0-coarse',
           '-pbe-bie-precond-coarse-aggregate-cap=8',
           '-pbe-bie-precond-coarse-verify=off']

# The submitted matrix's settings per preconditioner (Figure_2_source_matrix).
PRECONDITIONERS = {
    'diagonal': ['-pbe-bie-precond=diagonal'],
    'leaf_block_jacobi': ['-pbe-bie-precond=fmm-block-jacobi',
                          '-pbe-bie-precond-block-panel-cap=64'],
    'fmm_near_graph_ras': ['-pbe-bie-precond=fmm-overlap-schwarz',
                           '-pbe-bie-precond-block-panel-cap=32',
                           '-pbe-bie-precond-overlap-panel-cap=64'],
    'distance_knn_graph_ras': ['-pbe-bie-precond=distance-knn-overlap-schwarz',
                               '-pbe-bie-precond-block-panel-cap=32',
                               '-pbe-bie-precond-overlap-panel-cap=64'],
    'distance_knn_ilu0': ['-pbe-bie-precond=distance-knn-ilu0'],
    'fmm_near_graph_ilu0': ['-pbe-bie-precond=fmm-nearfield-ilu0'],
}

COHORT = ['barnase', '1a63', '1sh1', '1a7m', '1ubq']
FAGBI27_IDS = ['1ajj', '2erl', '1cbn', '1vii', '1fca', '1bbl', '2pde', '1sh1',
               '1vjw', '1uxc', '1ptq', '1bor', '1fxd', '1r69', '1mbg', '1bpi',
               '1hpt', '451c', '1svr', '1frd', '1a2s', '1neq', '1a63', '1a7m',
               '2go0', '1uv0', '4mth']


def amber(name):
    return PUBLIC / f'{name}_full_amber_ph7.pqr'


def charmm(name):
    return FAGBI27 / f'{name}_node_patch' / f'{name}_full_charmm_ph7.pqr'


def tag(h):
    return f'h{h:g}'.replace('.', 'p')


def stage_proteins():
    for name in FAGBI27_IDS:
        yield f'{name}_charmm', charmm(name), 1.25, ANALYTIC


# 451c at the protocol h = 1.25 A fails the operator gate of section 2.6: a
# spurious oscillating q on 40 cells refined in a deep pocket (the certificate's
# q_max_median, q_checkerboard_pairs).  The same protocol at h = 0.75 and 0.5,
# where the mode does not form, is its control.
def stage_control451c():
    for h in (0.75, 0.5):
        yield f'451c_charmm_{tag(h)}', charmm('451c'), h, ANALYTIC


HSWEEP_IDS = ['1vii', '1sh1', '1bpi', '1a63', '1a7m']


def stage_hsweep():
    for name in HSWEEP_IDS:
        for h in (1.5, 1.25, 1.0, 0.75, 0.6, 0.5):
            yield f'{name}_charmm_{tag(h)}', charmm(name), h, ANALYTIC


def stage_hsweep_lowsalt():
    # 1e-5 M: APBS (srfm mol) and the BIE then solve one boundary-value
    # problem, which the accuracy-matched cost comparison needs.
    for name in HSWEEP_IDS:
        for h in (1.25, 1.0, 0.75, 0.6, 0.5):
            yield (f'{name}_charmm_{tag(h)}', charmm(name), h,
                   ANALYTIC + ['-pbe-salt=1e-5'])


# The guarded p4 keeps its Green-jet cache in a 256 MB budget by default,
# which 1A63 and 1A7M overrun below h = 0.75 A (677 MB at h = 0.6).
P4_REF = P4 + ['-pbe-bie-fmm-local-cache-budget-mb=4096']


def stage_hsweep_p4():
    # The p1 far field grows with refinement (fixed 64-panel leaves cover a
    # smaller patch), so the converged reference is the guarded p4.
    for name in HSWEEP_IDS:
        for h in (1.25, 1.0, 0.75, 0.6, 0.5):
            yield f'{name}_charmm_p4_{tag(h)}', charmm(name), h, ANALYTIC + P4_REF


def stage_hsweep_lowsalt_p4():
    for name in HSWEEP_IDS:
        for h in (1.25, 1.0, 0.75, 0.6, 0.5):
            yield (f'{name}_charmm_p4_{tag(h)}', charmm(name), h,
                   ANALYTIC + P4_REF + ['-pbe-salt=1e-5'])


# Cheaper settings than the protocol, each adding one relaxation to the
# previous, for the accuracy/cost curve of the BIE (the balanced setting is
# chosen from it); a few single changes isolate their own effect.
NO_CHECKS = ['-pbe-bie-matvec-audit=off', '-pbe-bie-true-residual-probe=off']
LOOSE = ['-pbe-rms=1e-3', '-pbe-dphi=1e-3']
THETA060 = ['-pbe-bie-opening-theta=0.6']


def REFINE_AT(alpha):
    return NODE_PATCH + [f'-pbe-bie-ses-charge-refine={alpha}',
                         '-pbe-bie-ses-charge-refine-min-edge=0.1']


PRESETS = {
    'p0_protocol': ANALYTIC,
    'p1_nochecks': ANALYTIC + NO_CHECKS,
    'p2_rms1e-3': ANALYTIC + NO_CHECKS + LOOSE,
    'p3_theta060': ANALYTIC + NO_CHECKS + LOOSE + THETA060,
    'p4_remesh3': ANALYTIC + NO_CHECKS + LOOSE + THETA060 + ['-pbe-bie-ses-remesh-iterations=3'],
    'p5_remesh2': ANALYTIC + NO_CHECKS + LOOSE + THETA060 + ['-pbe-bie-ses-remesh-iterations=2'],
    'p6_remesh1': ANALYTIC + NO_CHECKS + LOOSE + THETA060 + ['-pbe-bie-ses-remesh-iterations=1'],
    'p7_alpha05': [a for a in ANALYTIC if 'charge-refine=' not in a]
                  + ['-pbe-bie-ses-charge-refine=0.5'] + NO_CHECKS + LOOSE + THETA060
                  + ['-pbe-bie-ses-remesh-iterations=2'],
    's_remesh2_only': ANALYTIC + ['-pbe-bie-ses-remesh-iterations=2'],
    's_theta060_only': ANALYTIC + THETA060,
    'q1_alpha05': REFINE_AT(0.5) + NO_CHECKS + LOOSE,
    'q2_alpha05_theta060': REFINE_AT(0.5) + NO_CHECKS + LOOSE + THETA060,
    'q3_alpha05_remesh3': REFINE_AT(0.5) + NO_CHECKS + LOOSE + ['-pbe-bie-ses-remesh-iterations=3'],
    'q4_alpha07': REFINE_AT(0.7) + NO_CHECKS + LOOSE,
    'q5_uniform': NODE_PATCH + NO_CHECKS + LOOSE,
    'q6_alpha06_remesh3': REFINE_AT(0.6) + NO_CHECKS + LOOSE + ['-pbe-bie-ses-remesh-iterations=3'],
    'q7_alpha05_remesh3_theta060': REFINE_AT(0.5) + NO_CHECKS + LOOSE + THETA060
                                   + ['-pbe-bie-ses-remesh-iterations=3'],
}


def stage_presets():
    for name in HSWEEP_IDS:
        for preset, extra in PRESETS.items():
            yield f'{name}_charmm_{preset}', charmm(name), 1.25, extra


def stage_presets27():
    # The candidate balanced settings on the 27-protein set (0.15 M).
    for name in FAGBI27_IDS:
        for preset in ('p7_alpha05', 'q3_alpha05_remesh3', 'q4_alpha07',
                       'q6_alpha06_remesh3', 'q7_alpha05_remesh3_theta060'):
            yield f'{name}_charmm_{preset}', charmm(name), 1.25, PRESETS[preset]


def stage_presets_lowsalt():
    for name in HSWEEP_IDS:
        for preset, extra in PRESETS.items():
            yield f'{name}_charmm_{preset}', charmm(name), 1.25, extra + ['-pbe-salt=1e-5']


# The sphere-patch surface as a physical model: the clipped union of atomic
# spheres is the van der Waals surface.  At 1e-5 M APBS with srfm mol and a
# zero probe (tabi_rsos_apbs_sweep.py --mode vdw_lowsalt) solves the same
# boundary-value problem.
VDW_IDS = ['1vii', '1sh1', '1bpi']
SPHERE_VDW = ['-pbe-bie-primary-surface=sphere-patch',
              '-pbe-bie-discretization=centroid', '-pbe-salt=1e-5']


def stage_vdw():
    for name in VDW_IDS:
        for points in (8, 16):
            for route, extra in (('p1', []), ('p4', P4_REF)):
                yield (f'{name}_charmm_vdw_n{points}_{route}', charmm(name), 1.25,
                       SPHERE_VDW + [f'-pbe-bie-surface-points={points}'] + extra)


def stage_matrix():
    surfaces = {'analytic_ses': ANALYTIC, 'grid_ses': GRID, 'sphere_patch': SPHERE}
    for name in COHORT:
        for surface, surface_args in surfaces.items():
            for label, precond in PRECONDITIONERS.items():
                # The 512-target audit of section 2.6, as every other row
                # (the submitted matrix had sampled 256).
                yield (f'{name}_{surface}_{label}', amber(name), 1.25,
                       surface_args + precond)


def binding_args(complex_pqr, receptor, ligand):
    return ['-pbe-binding-cycle=pqr-components',
            f'-pbe-binding-receptor-pqr={Path(receptor).as_posix()}',
            f'-pbe-binding-ligand-pqr={Path(ligand).as_posix()}']


def stage_binding():
    lzg = (SIXLZG / '6lzg_complex_charmm_ph7.pqr',
           SIXLZG / '6lzg_ace2_chain_A_charmm_ph7.pqr',
           SIXLZG / '6lzg_rbd_chain_B_charmm_ph7.pqr')
    yield '6lzg_charmm', lzg[0], 1.25, ANALYTIC + binding_args(*lzg)
    # Table 6 as submitted: cap-8 coarse aggregation in a 1 GiB budget.
    yield ('6lzg_charmm_coarse8_1g', lzg[0], 1.25,
           ANALYTIC + binding_args(*lzg) + COARSE8
           + ['-pbe-bie-precond-memory-budget-mb=1024'])
    bb = (BINDING_JPCB / '1brs_complex_amber_ph7.pqr',
          BINDING_JPCB / '1brs_barnase_chain_A_amber_ph7.pqr',
          BINDING_JPCB / '1brs_barstar_chain_D_amber_ph7.pqr')
    yield '1brs_amber', bb[0], 1.25, ANALYTIC + binding_args(*bb)


def stage_p1p4():
    theta = ['-pbe-bie-opening-theta=0.50']
    for name in ['1sh1', 'barnase', '1a63']:
        yield f'{name}_p1_theta050', amber(name), 1.25, ANALYTIC + theta
        yield f'{name}_p4_theta050', amber(name), 1.25, ANALYTIC + P4 + theta


def stage_spacing():
    for name in ['1a63', '1sh1', 'barnase']:
        for h in (1.25, 1.0, 0.8, 0.75, 0.6, 0.5):
            yield f'{name}_{tag(h)}', amber(name), h, ANALYTIC


def stage_cap():
    for name in ['barnase', '1a63', '1sh1']:
        # The guarded p4 accepts caps other than 64 only as a declared sweep.
        for route, route_args in (
                ('p4', P4 + ['-pbe-bie-fmm-local-ilu-cap-sweep=true']), ('p1', [])):
            for cap in (32, 64, 128):
                yield (f'{name}_{route}_cap{cap}', amber(name), 1.25,
                       ANALYTIC + route_args
                       + [f'-pbe-bie-precond-ilu-neighbor-cap={cap}'])


def stage_coarse():
    yield '1a63_ilu0', amber('1a63'), 1.25, ANALYTIC
    yield '1a63_ilu0_coarse8', amber('1a63'), 1.25, ANALYTIC + COARSE8


def stage_s1():
    for name in ['1ubq', '1vii', '2ci2']:
        pqr = (S1_INPUTS / f'titan_guarded_p4_analytic_nodepatch_p4_s1_v1_{name}'
               / f'{name}_full_amber_ph7.pqr')
        yield f'{name}_p1', pqr, 1.25, ANALYTIC
        yield f'{name}_p4', pqr, 1.25, ANALYTIC + P4


def stage_dense():
    # Whole small proteins: the submitted table's mini fragments are cut out of
    # larger chains and their SES carries thin plates and slivers that no
    # whole protein has, which makes them a poor calibration of the protocol.
    for name in ['1vii', '1cbn', '2erl', '1sh1']:
        pqr = charmm(name)
        yield f'{name}_dense', pqr, 1.25, ANALYTIC + DENSE + TIGHT
        yield f'{name}_dense_zone16', pqr, 1.25, ANALYTIC + DENSE16 + TIGHT
        yield f'{name}_p1', pqr, 1.25, ANALYTIC + TIGHT
        yield f'{name}_p4', pqr, 1.25, ANALYTIC + P4 + TIGHT


def stage_xgen():
    # Figure 4A: the two surface generators under one discretization each.
    # The analytic rows are uniform (no charge refinement, which only the
    # analytic engine has); the protocol row is added for reference.
    for h in (1.25, 1.0, 0.75, 0.6, 0.5):
        for label, extra, env in (
                ('grid_centroid', GRID, None),
                ('analytic_centroid', ['-pbe-bie-discretization=centroid'], None),
                ('grid_np_flat', GRID_NP, None),
                ('analytic_np_flat', NODE_PATCH, FLAT_CELL),
                ('analytic_np_curved', NODE_PATCH, None),
                ('analytic_protocol', ANALYTIC, None)):
            yield f'1a63_{label}_{tag(h)}', amber('1a63'), h, extra, 4, env


def stage_threads():
    # Figure 4B: one and four threads.
    for surface, surface_args in (('analytic_protocol', ANALYTIC),
                                  ('grid_np_flat', GRID_NP)):
        for h in (1.25, 1.0, 0.75, 0.6, 0.5):
            for nt in (1, 4):
                yield (f'1a63_{surface}_{tag(h)}_nt{nt}', amber('1a63'), h,
                       surface_args, nt)


# The sphere-patch stress studies (Tables S2, S6) and the p4 spacing table
# (S5), with the submitted settings: theta 0.6, rms 1e-6 within 720
# iterations, 1800 s per FMM run and 3600 s for the direct sum, four threads.
STRESS_IDS = {'barnase': PUBLIC, '1a63': PUBLIC, '1sh1': PUBLIC, '1ubq': PUBLIC,
              '2ci2': ROOT / 'test/pbe_inputs_public_z4q',
              '1vii': ROOT / 'test/pbe_inputs_public_z4q'}
STRESS_TOL = ['-pbe-rms=1e-6', '-pbe-dphi=1e-6', '-pbe-max-iter=720']
THETA06 = ['-pbe-bie-opening-theta=0.6']
# The Green-jet cache of the clipped surface overruns the 256 MB budget that
# held the submitted whole-triangle one (1SH1: 910 MB), so the stress runs use
# the reference budget.
P4_STRESS = SPHERE + P4_REF + THETA06 + STRESS_TOL


def stress_pqr(name):
    return STRESS_IDS[name] / f'{name}_full_amber_ph7.pqr'


def stage_s2():
    variants = {
        'baseline': [],
        'geometry_off': ['-pbe-bie-fmm-local-geometry-guard=off',
                         '-pbe-bie-fmm-local-geometry-theta=0'],
        'theta0p4': ['-pbe-bie-opening-theta=0.4'],
        'theta0p5': ['-pbe-bie-opening-theta=0.5'],
        'theta0p7': ['-pbe-bie-opening-theta=0.7'],
        # The guarded p4 refuses to run without its tail estimator.
        'tail_off': ['-pbe-bie-fmm-local-refine=off'],
    }
    for variant, extra in variants.items():
        for name in STRESS_IDS:
            yield (f'{name}_{variant}', stress_pqr(name), 1.25,
                   P4_STRESS + extra, 4, None, 1800)


def stage_s5():
    pqr = stress_pqr('1vii')
    for h in (1.5, 1.25, 1.0, 0.75, 0.5):
        yield (f'1vii_p4_{tag(h)}', pqr, h,
               ANALYTIC + P4 + THETA06 + STRESS_TOL, 4, None, 1800)


def stage_s6():
    p1 = SPHERE + THETA06 + STRESS_TOL
    stages = {
        'tabi_direct': (p1 + ['-pbe-bie-sum=direct', '-pbe-bie-fmm-m2l=off',
                              '-pbe-bie-precond=diagonal',
                              '-pbe-bie-matvec-audit=off'], 3600),
        'p1_diagonal': (p1 + ['-pbe-bie-precond=diagonal'], 1800),
        'p1_block': (p1 + ['-pbe-bie-precond=fmm-block-jacobi',
                           '-pbe-bie-precond-block-panel-cap=32'], 1800),
        'p1_ilu0': (p1, 1800),
        'p4_ilu0': (P4_STRESS, 1800),
    }
    for stage, (extra, timeout) in stages.items():
        for name in STRESS_IDS:
            yield f'{name}_{stage}', stress_pqr(name), 1.25, extra, 4, None, timeout


STAGES = {'proteins': stage_proteins, 'control451c': stage_control451c,
          'hsweep': stage_hsweep,
          'hsweep_lowsalt': stage_hsweep_lowsalt, 'hsweep_p4': stage_hsweep_p4,
          'hsweep_lowsalt_p4': stage_hsweep_lowsalt_p4, 'presets': stage_presets,
          'presets_lowsalt': stage_presets_lowsalt, 'presets27': stage_presets27, 'vdw': stage_vdw,
          'matrix': stage_matrix,
          'binding': stage_binding, 'p1p4': stage_p1p4,
          'spacing': stage_spacing, 'cap': stage_cap, 'coarse': stage_coarse,
          's1': stage_s1, 'dense': stage_dense, 'xgen': stage_xgen,
          'threads': stage_threads, 's2': stage_s2, 's5': stage_s5,
          's6': stage_s6}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--titan', type=Path, required=True)
    parser.add_argument('--stages', default=','.join(STAGES))
    parser.add_argument('--only', default='',
                        help='comma-separated case names to run within the stages')
    parser.add_argument('--exclude', default='',
                        help='comma-separated substrings; matching cases are skipped')
    parser.add_argument('--work-root', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=10800)
    parser.add_argument('--resume', action='store_true',
                        help='read back cases that already finished instead of rerunning')
    args = parser.parse_args()
    only = set(filter(None, args.only.split(',')))
    exclude = [e for e in args.exclude.split(',') if e]
    for stage in args.stages.split(','):
        out_dir = args.work_root / stage
        out_dir.mkdir(parents=True, exist_ok=True)
        rows = []
        for spec in STAGES[stage]():
            name, pqr, h, extra = spec[:4]
            if only and name not in only:
                continue
            if any(e in name for e in exclude):
                continue
            nt = spec[4] if len(spec) > 4 else 4
            env = spec[5] if len(spec) > 5 else None
            timeout = spec[6] if len(spec) > 6 else args.timeout
            run_dir = out_dir / name
            row = {'case': name, 'variant': stage, 'h_A': h}
            if args.resume and (run_dir / 'case_meta.json').exists():
                row.update(runner.collect_case(run_dir))
                note = 'resumed'
            else:
                row.update(runner.run_case(args.titan, run_dir, pqr, h, extra,
                                           env, nt, timeout))
                note = ''
            rows.append(row)
            print(f"{stage:8s} {name:40s} n={row['panels']:<7d} "
                  f"it={row['iterations']:<4d} E={row['energy_kcalmol']:12.4f} "
                  f"eop={row['op_audit_rel']:.2e} {row['wall_s']:7.1f}s "
                  f"{row['peak_rss_mb']:6.0f}MB [{row['status']}|{row['exit']}] "
                  f"{row['precond_effective']} {note}", flush=True)
            runner.write_rows(out_dir / f'{stage}.csv', rows)
    return 0


if __name__ == '__main__':
    sys.exit(main())
