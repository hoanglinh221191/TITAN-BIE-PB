#!/usr/bin/env python3
"""TITAN BIE-PB against APBS, two ways (RSOS revision, reviewer 2, point 9c).

1. Equal accuracy (work-precision).  At vanishing salt (1e-5 M) the BIE and
   APBS with srfm mol solve one boundary-value problem, so both are measured
   against one reference: TITAN under the guarded p4 far field at its finest
   spacing.  For each target error (1, 0.5, 0.2%) the cheapest level of each
   program that reaches it is reported with its wall time and peak memory.

2. Balanced settings.  Each program at the setting a user would run: TITAN at
   its protocol (h = 1.25 A, p1) and APBS at the 0.5 A fine grid of the
   manuscript.  At 1e-5 M both errors are against the common reference; at
   0.15 M, where the two use different ion models, each program is measured
   against its own converged value (TITAN p4 finest; APBS finest grid), and
   the cost ratios are reported.  With --balanced-preset the TITAN side is
   that preset of tabi_rsos_campaign.py (stages presets_lowsalt / presets)
   instead of the protocol.

    python test/tabi_rsos_compare.py \
        --titan-lowsalt campaign_timing/hsweep_lowsalt/hsweep_lowsalt.csv \
        --titan-lowsalt-ref campaign/hsweep_lowsalt_p4/hsweep_lowsalt_p4.csv \
        --apbs-lowsalt apbs_sweep/apbs_lowsalt.csv \
        --titan-salt campaign_timing/hsweep/hsweep.csv \
        --titan-salt-ref campaign/hsweep_p4/hsweep_p4.csv \
        --apbs-salt apbs_sweep/apbs_paper.csv --out compare
"""

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

TARGETS = (1e-2, 5e-3, 2e-3)


def load(path, level_key):
    """case -> [{level, energy, wall, rss}] sorted coarse to fine."""
    data = defaultdict(list)
    if path is None or not Path(path).exists():
        return data
    for row in csv.DictReader(open(path, newline='')):
        energy = float(row['energy_kcalmol'])
        if not math.isfinite(energy) or str(row.get('exit', '0')) != '0':
            continue
        data[row['case'].split('_')[0]].append({
            'level': float(row[level_key]), 'energy': energy,
            'wall': float(row['wall_s']), 'rss': float(row['peak_rss_mb'])})
    for rows in data.values():
        rows.sort(key=lambda r: -r['level'])
    return data


def reference(ref_sweep, fallback_sweep, case):
    """The finest guarded-p4 value, else the finest p1 value (flagged)."""
    if ref_sweep.get(case):
        return ref_sweep[case][-1]['energy'], 'p4', ref_sweep[case][-1]['level']
    return fallback_sweep[case][-1]['energy'], 'p1-fallback', fallback_sweep[case][-1]['level']


def load_preset(path, preset):
    """case -> {energy, wall, rss} for one preset of a presets stage."""
    data = {}
    if path is None or not Path(path).exists():
        return data
    for row in csv.DictReader(open(path, newline='')):
        energy = float(row['energy_kcalmol'])
        if (not row['case'].endswith('_' + preset) or not math.isfinite(energy)
                or str(row.get('exit', '0')) != '0'):
            continue
        data[row['case'].split('_')[0]] = {
            'level': float(row['h_A']), 'energy': energy,
            'wall': float(row['wall_s']), 'rss': float(row['peak_rss_mb'])}
    return data


def rel(e, ref):
    return abs(e - ref) / abs(ref)


def at_level(rows, level):
    found = [r for r in rows if abs(r['level'] - level) < 1e-9]
    return found[0] if found else None


def write(path, rows):
    if not rows:
        return
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--titan-lowsalt', type=Path, required=True)
    ap.add_argument('--titan-lowsalt-ref', type=Path)
    ap.add_argument('--apbs-lowsalt', type=Path, required=True)
    ap.add_argument('--titan-salt', type=Path, required=True)
    ap.add_argument('--titan-salt-ref', type=Path)
    ap.add_argument('--apbs-salt', type=Path, required=True)
    ap.add_argument('--titan-h', type=float, default=1.25)
    ap.add_argument('--apbs-spacing', type=float, default=0.5)
    ap.add_argument('--balanced-preset', default='')
    ap.add_argument('--titan-presets-lowsalt', type=Path)
    ap.add_argument('--titan-presets-salt', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    t_low, t_low_ref = load(a.titan_lowsalt, 'h_A'), load(a.titan_lowsalt_ref, 'h_A')
    p_low = load(a.apbs_lowsalt, 'spacing_A')
    t_salt, t_salt_ref = load(a.titan_salt, 'h_A'), load(a.titan_salt_ref, 'h_A')
    p_salt = load(a.apbs_salt, 'spacing_A')
    cases = sorted(set(t_low) & set(p_low))

    # 1. Equal accuracy at 1e-5 M.
    levels, targets = [], []
    for case in cases:
        ref, ref_kind, ref_h = reference(t_low_ref, t_low, case)
        for program, rows in (('TITAN', [r for r in t_low[case] if r['level'] <= a.titan_h + 1e-9]),
                              ('APBS', p_low[case])):
            for r in rows:
                levels.append({'case': case, 'program': program, 'level_A': r['level'],
                               'energy_kcalmol': r['energy'], 'reference_kcalmol': ref,
                               'reference': f'{ref_kind} h={ref_h:g}',
                               'rel_error': rel(r['energy'], ref),
                               'wall_s': r['wall'], 'peak_rss_mb': r['rss']})
            for target in TARGETS:
                ok = [r for r in rows if rel(r['energy'], ref) <= target]
                best = min(ok, key=lambda r: r['wall']) if ok else None
                targets.append({'case': case, 'program': program, 'target': target,
                                'reached': bool(best),
                                'level_A': best['level'] if best else '',
                                'wall_s': round(best['wall'], 2) if best else '',
                                'peak_rss_mb': round(best['rss']) if best else ''})
    write(a.out / 'equal_accuracy_levels.csv', levels)
    write(a.out / 'equal_accuracy_targets.csv', targets)

    # 2. Balanced settings.
    if a.balanced_preset:
        b_low = load_preset(a.titan_presets_lowsalt, a.balanced_preset)
        b_salt = load_preset(a.titan_presets_salt, a.balanced_preset)
        balanced_low = lambda case: b_low.get(case)
        balanced_salt = lambda case: b_salt.get(case)
        setting = a.balanced_preset
    else:
        balanced_low = lambda case: at_level(t_low[case], a.titan_h)
        balanced_salt = lambda case: at_level(t_salt[case], a.titan_h)
        setting = f'h={a.titan_h:g} A p1'
    balanced = []
    for case in cases:
        ref, ref_kind, _ = reference(t_low_ref, t_low, case)
        t, p = balanced_low(case), at_level(p_low[case], a.apbs_spacing)
        if t and p:
            balanced.append({'case': case, 'salt': '1e-5 M (common reference)',
                             'titan_rel_error': rel(t['energy'], ref),
                             'apbs_rel_error': rel(p['energy'], ref),
                             'titan_wall_s': t['wall'], 'apbs_wall_s': p['wall'],
                             'titan_rss_mb': t['rss'], 'apbs_rss_mb': p['rss'],
                             'wall_ratio_t_over_a': t['wall'] / p['wall'],
                             'rss_ratio_t_over_a': t['rss'] / p['rss'],
                             'reference': ref_kind})
        if case in t_salt and case in p_salt:
            ref_t, ref_kind_t, _ = reference(t_salt_ref, t_salt, case)
            ref_p = p_salt[case][-1]['energy']
            t, p = balanced_salt(case), at_level(p_salt[case], a.apbs_spacing)
            if t and p:
                balanced.append({'case': case, 'salt': '0.15 M (own references)',
                                 'titan_rel_error': rel(t['energy'], ref_t),
                                 'apbs_rel_error': rel(p['energy'], ref_p),
                                 'titan_wall_s': t['wall'], 'apbs_wall_s': p['wall'],
                                 'titan_rss_mb': t['rss'], 'apbs_rss_mb': p['rss'],
                                 'wall_ratio_t_over_a': t['wall'] / p['wall'],
                                 'rss_ratio_t_over_a': t['rss'] / p['rss'],
                                 'reference': f'TITAN {ref_kind_t}; APBS finest'})
    write(a.out / 'balanced.csv', balanced)

    print('Equal accuracy (1e-5 M): cheapest level within target')
    for target in TARGETS:
        for program in ('TITAN', 'APBS'):
            rows = [t for t in targets if t['program'] == program and t['target'] == target]
            reached = [t for t in rows if t['reached']]
            walls = [t['wall_s'] for t in reached]
            rss = [t['peak_rss_mb'] for t in reached]
            print(f'  {target:.1%} {program:5s} reached {len(reached)}/{len(rows)}'
                  + (f'  wall {min(walls):.1f}-{max(walls):.1f} s, RSS {min(rss)}-{max(rss)} MB'
                     if reached else ''))
    print(f'Balanced (TITAN {setting}, APBS {a.apbs_spacing:g} A):')
    for salt in sorted({b['salt'] for b in balanced}):
        rows = [b for b in balanced if b['salt'] == salt]
        med = lambda k: statistics.median(b[k] for b in rows)
        print(f'  {salt}: error TITAN {med("titan_rel_error"):.2%} APBS {med("apbs_rel_error"):.2%};'
              f' wall T/A {med("wall_ratio_t_over_a"):.2f}; RSS T/A {med("rss_ratio_t_over_a"):.3f}')
    if not t_low_ref or not t_salt_ref:
        print('note: guarded-p4 reference missing for some sweeps; p1 finest used', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
