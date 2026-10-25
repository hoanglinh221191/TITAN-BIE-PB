#!/usr/bin/env python3
"""Accuracy against cost of the TITAN settings (presets) and of APBS grids.

Each TITAN preset (tabi_rsos_campaign.py stages presets / presets_lowsalt)
and each APBS grid spacing (tabi_rsos_apbs_sweep.py) is measured against the
converged reference, the guarded-p4 TITAN value at its finest spacing
(stages hsweep_p4 / hsweep_lowsalt_p4).  At 1e-5 M both programs solve one
boundary-value problem; at 0.15 M APBS (smol, ion radius 2 A) is measured
against its own finest grid instead.  Writes a per-preset summary (median
and largest error, median wall time and peak memory, and their ratio to
APBS at 0.5 A) and the error-cost figure.

    python test/tabi_rsos_presets_curve.py --final D:/data_MD/tabi/rsos_v06/campaign_final \
        --campaign D:/data_MD/tabi/rsos_v06/campaign --out presets_curve
"""

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


def rows(path):
    return list(csv.DictReader(open(path, newline=''))) if Path(path).exists() else []


def finest_reference(path):
    best = {}
    for r in rows(path):
        case = r['case'].split('_')[0]
        h = float(r['h_A'])
        e = float(r['energy_kcalmol'])
        if math.isfinite(e) and (case not in best or h < best[case][0]):
            best[case] = (h, e)
    return {k: v[1] for k, v in best.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--final', type=Path, required=True)
    ap.add_argument('--campaign', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    salts = {
        '1e-5 M': ('presets_lowsalt', 'hsweep_lowsalt_p4', 'apbs_lowsalt.csv', True),
        '0.15 M': ('presets', 'hsweep_p4', 'apbs_paper.csv', False),
    }
    summary, points = [], defaultdict(list)
    for salt, (stage, ref_stage, apbs_file, common) in salts.items():
        ref = finest_reference(a.campaign / ref_stage / f'{ref_stage}.csv')
        apbs = defaultdict(list)
        for r in rows(a.final / 'apbs_sweep' / apbs_file):
            apbs[r['case'].split('_')[0]].append(r)
        apbs05 = {c: [r for r in rs if abs(float(r['spacing_A']) - 0.5) < 1e-9][0]
                  for c, rs in apbs.items()}
        by_preset = defaultdict(list)
        for r in rows(a.final / stage / f'{stage}.csv'):
            case, _, preset = r['case'].split('_', 2)
            e = float(r['energy_kcalmol'])
            if case not in ref or not math.isfinite(e):
                continue
            err = abs(e - ref[case]) / abs(ref[case])
            wall, rss = float(r['wall_s']), float(r['peak_rss_mb'])
            by_preset[preset].append({'case': case, 'err': err, 'wall': wall, 'rss': rss,
                                      'wall_ratio': wall / float(apbs05[case]['wall_s']),
                                      'rss_ratio': rss / float(apbs05[case]['peak_rss_mb'])})
            points[(salt, 'TITAN ' + preset)].append((wall, err))
        for preset, rs in sorted(by_preset.items()):
            med = lambda k: statistics.median(x[k] for x in rs)
            summary.append({'salt': salt, 'program': 'TITAN', 'setting': preset,
                            'cases': len(rs), 'median_rel_error': med('err'),
                            'max_rel_error': max(x['err'] for x in rs),
                            'median_wall_s': med('wall'), 'median_rss_mb': med('rss'),
                            'median_wall_over_apbs05': med('wall_ratio'),
                            'median_rss_over_apbs05': med('rss_ratio')})
        for spacing in sorted({float(r['spacing_A']) for rs in apbs.values() for r in rs},
                              reverse=True):
            errs, walls, rsss = [], [], []
            for case, rs in apbs.items():
                r = [x for x in rs if abs(float(x['spacing_A']) - spacing) < 1e-9]
                if not r or case not in ref:
                    continue
                r = r[0]
                own = ref[case] if common else float(
                    min(rs, key=lambda x: float(x['spacing_A']))['energy_kcalmol'])
                err = abs(float(r['energy_kcalmol']) - own) / abs(own)
                errs.append(err)
                walls.append(float(r['wall_s']))
                rsss.append(float(r['peak_rss_mb']))
                points[(salt, 'APBS')].append((float(r['wall_s']), err))
            if errs:
                summary.append({'salt': salt, 'program': 'APBS', 'setting': f'{spacing:g} A',
                                'cases': len(errs), 'median_rel_error': statistics.median(errs),
                                'max_rel_error': max(errs),
                                'median_wall_s': statistics.median(walls),
                                'median_rss_mb': statistics.median(rsss),
                                'median_wall_over_apbs05': '', 'median_rss_over_apbs05': ''})
    with open(a.out / 'presets_summary.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)
    for s in summary:
        ratio = (f"  wall/APBS0.5 {s['median_wall_over_apbs05']:.2f}"
                 f" RSS/APBS0.5 {s['median_rss_over_apbs05']:.2f}"
                 if s['program'] == 'TITAN' else '')
        print(f"{s['salt']:7s} {s['program']:5s} {s['setting']:16s} err median "
              f"{100 * s['median_rel_error']:.3f}% max {100 * s['max_rel_error']:.3f}%  "
              f"wall {s['median_wall_s']:.1f}s RSS {s['median_rss_mb']:.0f}MB{ratio}")

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for ax, salt in zip(axes, salts):
        for s in [x for x in summary if x['salt'] == salt]:
            marker = 'o' if s['program'] == 'TITAN' else 's'
            color = '#1f77b4' if s['program'] == 'TITAN' else '#d62728'
            ax.loglog(s['median_wall_s'], 100 * s['median_rel_error'], marker, color=color)
            ax.annotate(s['setting'], (s['median_wall_s'], 100 * s['median_rel_error']),
                        fontsize=7, xytext=(3, 3), textcoords='offset points')
        ax.set_xlabel('median wall time (s, 4 threads)')
        ax.set_ylabel('median relative error (%)')
        ax.set_title(f'{salt}: TITAN settings (o) and APBS grids (s)')
    fig.savefig(a.out / 'presets_curve.png', dpi=200)
    fig.savefig(a.out / 'presets_curve.svg')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
