#!/bin/bash
# After run_night_v06_timing.sh ("v06 timing done"), for the internal review of 21 Oct:
#  F4  the re-audited p4 sweeps (campaign_verify_p4) replace the 128-row-audit ones, if every case agrees;
#  F1  the low-salt APBS rows of the proteins with enclosed cavities come from the cavity-filled sweep, so
#      APBS solves TITAN's problem (external-only cavities);
# then compare_idle and presets_curve_idle are derived again from both.  Stops at the first failed check.
#   bash morning_after_review.sh
set -eu
REPO=/c/Users/ADMIN/OneDrive/Documents/antigravity/elec_vdw
OUT=D:/data_MD/tabi/rsos_v06
MAIN=$OUT/campaign_main_idle
FINAL=$OUT/campaign_final
OLD=$MAIN/_campaign_binary_timing
cd "$OUT/campaign"
py -3 swap_p4_reference.py --apply
py -3 - <<'PY'
import csv
from pathlib import Path
F = Path('D:/data_MD/tabi/rsos_v06/campaign_final')
base = list(csv.DictReader(open(F / 'apbs_sweep/apbs_lowsalt.csv')))
cav = {(r['case'], r['spacing_A']): r for r in csv.DictReader(open(F / 'apbs_sweep_cavfilled/apbs_lowsalt.csv'))}
assert all(r['exit'] == '0' for r in cav.values()), 'a cavity-filled APBS run failed'
out = [dict(cav.get((r['case'], r['spacing_A']), r)) for r in base]
assert sum((r['case'], r['spacing_A']) in cav for r in base) == len(cav)
with open(F / 'apbs_sweep/apbs_lowsalt_matched.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(base[0]))
    w.writeheader()
    w.writerows(out)
print('apbs_lowsalt_matched.csv:', len(cav), 'of', len(base), 'rows from the cavity-filled sweep')
PY
cd "$REPO"
mv "$OUT/compare_idle" "$OLD/compare_idle_before_review"
py -3 -I test/tabi_rsos_compare.py --titan-lowsalt $MAIN/hsweep_lowsalt/hsweep_lowsalt.csv \
  --titan-lowsalt-ref $MAIN/hsweep_lowsalt_p4/hsweep_lowsalt_p4.csv \
  --apbs-lowsalt $FINAL/apbs_sweep/apbs_lowsalt_matched.csv --titan-salt $MAIN/hsweep/hsweep.csv \
  --titan-salt-ref $MAIN/hsweep_p4/hsweep_p4.csv --apbs-salt $FINAL/apbs_sweep/apbs_paper.csv \
  --out $OUT/compare_idle > "$MAIN/compare_idle.log" 2>&1
echo "compare (matched cavities, re-audited p4) exit $?"
mv "$OUT/presets_curve_idle" "$OLD/presets_curve_idle_before_review"
T=$OUT/presets_curve_idle/inputs
mkdir -p $T/presets $T/presets_lowsalt $T/apbs_sweep
cp $MAIN/presets/presets.csv $T/presets/
cp $MAIN/presets_lowsalt/presets_lowsalt.csv $T/presets_lowsalt/
cp $FINAL/apbs_sweep/apbs_lowsalt_matched.csv $T/apbs_sweep/apbs_lowsalt.csv
cp $FINAL/apbs_sweep/apbs_paper.csv $T/apbs_sweep/
py -3 -I test/tabi_rsos_presets_curve.py --final $T --campaign $MAIN --out $OUT/presets_curve_idle \
  > "$MAIN/presets_curve_idle.log" 2>&1
echo "presets curve exit $?"
echo "morning after review done"
