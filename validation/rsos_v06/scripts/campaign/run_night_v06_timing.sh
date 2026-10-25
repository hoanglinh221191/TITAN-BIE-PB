#!/bin/bash
# Re-times on the release binary every manuscript wall time still measured on
# the campaign binary, with APBS rerun in the same block so both programs are
# timed under one condition.  Each APBS block follows the TITAN stage it is
# compared with.  Results other than times are bit-identical between the two
# binaries (checked on 85 rows), so only times and RSS change.
#   proteins + APBS 27 at 0.5 A          Figure 5, Table S7
#   binding  + APBS binding components   Table 6
#   hsweep_lowsalt, hsweep + APBS sweeps (lowsalt, paper; 0.6-0.2 A)
#                                        Abstract, 3.4, S14, S15 (compare_idle)
#   presets, presets_lowsalt, presets27  3.4 cost-accuracy curve, Table R4
#   threads, xgen                        Figure 4
#   cap, spacing, coarse                 Tables S3, S12, S13
#   matrix Grid-SES and analytic-SES rows (Tables 4, 5, S9; Figure 3)
#   control451c                          Table S7 note
#   PAFMPB against TITAN (1SH1, 1A63)     response to Reviewer 3, comment 1
# Then the comparisons are derived again, the Kirkwood stages are rerun, and,
# last, the two p4 spacing sweeps are rerun into a separate work root as an
# energy check (their energies are references only; no time is reported).
# The previous rows move to campaign_main_idle/_campaign_binary_timing/<stage>/,
# the previous APBS runs to campaign_final/_apbs_16oct/; nothing is deleted.
#   bash run_night_v06_timing.sh D:/data_MD/tabi/rsos_v06/bin/titan_f54a0507.exe "2026-10-21 23:00"
set -u
BIN=$1
START=${2:-}
REPO=/c/Users/ADMIN/OneDrive/Documents/antigravity/elec_vdw
OUT=D:/data_MD/tabi/rsos_v06
MAIN=$OUT/campaign_main_idle
OLD=$MAIN/_campaign_binary_timing
FINAL=$OUT/campaign_final
AOLD=$FINAL/_apbs_16oct
if [ -n "$START" ]; then
  target=$(date -d "$START" +%s)
  while [ "$(date +%s)" -lt "$target" ]; do sleep 60; done
fi
park() {  # park <stage> <glob> [<glob> ...]: move matching run folders aside, keep a CSV copy
  local stage=$1; shift
  mkdir -p "$OLD/$stage"
  [ -f "$MAIN/$stage/$stage.csv" ] && cp -n "$MAIN/$stage/$stage.csv" "$OLD/$stage/$stage.csv"
  for g in "$@"; do
    for d in $MAIN/$stage/$g; do
      [ -d "$d" ] && mv "$d" "$OLD/$stage/" || true
    done
  done
}
for stage in proteins control451c binding threads xgen hsweep_lowsalt hsweep presets presets_lowsalt \
             presets27 cap spacing coarse; do
  park $stage '*'
done
park matrix '*_grid_ses_*' '*_analytic_ses_*'
for k in kirkwood_R10 kirkwood_R2 kirkwood_R10_p4 kirkwood_R2_p4; do
  mkdir -p "$OLD/kirkwood"
  [ -d "$MAIN/$k" ] && mv "$MAIN/$k" "$OLD/kirkwood/"
  [ -f "$MAIN/$k.log" ] && mv "$MAIN/$k.log" "$OLD/kirkwood/"
done
mkdir -p "$AOLD/apbs_sweep"
for d in apbs_27 apbs_binding; do
  [ -d "$FINAL/$d" ] && mv "$FINAL/$d" "$AOLD/"
  [ -f "$FINAL/$d.log" ] && mv "$FINAL/$d.log" "$AOLD/"
done
for mode in lowsalt paper; do
  for d in $FINAL/apbs_sweep/*_${mode}_g*; do
    case "$d" in *_vdw_lowsalt_g*) continue ;; esac
    [ -d "$d" ] && mv "$d" "$AOLD/apbs_sweep/"
  done
  [ -f "$FINAL/apbs_sweep/apbs_$mode.csv" ] && mv "$FINAL/apbs_sweep/apbs_$mode.csv" "$AOLD/apbs_sweep/"
  [ -f "$FINAL/apbs_sweep_$mode.log" ] && mv "$FINAL/apbs_sweep_$mode.log" "$AOLD/"
done

cd "$REPO" || exit 1
export PATH=/c/msys64/ucrt64/bin:$PATH
echo "v06 timing start $BIN $(sha256sum "$BIN" | cut -c1-16) $(git rev-parse --short=8 HEAD) $(date)"
SKIP=""
for name in 1a63 1a7m; do
  for precond in diagonal leaf_block_jacobi fmm_near_graph_ras distance_knn_graph_ras; do
    SKIP="$SKIP${SKIP:+,}${name}_sphere_patch_${precond}"
  done
done
campaign() {
  python -u test/tabi_rsos_campaign.py --resume --titan "$BIN" --work-root "$MAIN" \
    --stages "$1" ${2:+--exclude "$2"} >> "$MAIN/campaign_main.log" 2>&1
  echo "$1 exit $? $(date)"
}
LISTS=$(python - <<'PY'
import sys
sys.path.insert(0, 'test')
import tabi_rsos_campaign as c
p27 = ','.join(c.charmm(n).as_posix() for n in c.FAGBI27_IDS)
lzg = ','.join(p.as_posix() for p in (c.SIXLZG / '6lzg_complex_charmm_ph7.pqr',
                                      c.SIXLZG / '6lzg_ace2_chain_A_charmm_ph7.pqr',
                                      c.SIXLZG / '6lzg_rbd_chain_B_charmm_ph7.pqr'))
bb = ','.join(p.as_posix() for p in (c.BINDING_JPCB / '1brs_complex_amber_ph7.pqr',
                                     c.BINDING_JPCB / '1brs_barnase_chain_A_amber_ph7.pqr',
                                     c.BINDING_JPCB / '1brs_barstar_chain_D_amber_ph7.pqr'))
sweep = ','.join(c.charmm(n).as_posix() for n in c.HSWEEP_IDS)
print(p27 + ' ' + lzg + ',' + bb + ' ' + sweep)
PY
)
read -r P27 BINDING SWEEP <<< "$LISTS"
apbs() {  # apbs <name> <pqr list> <spacings> <mode> <work root>
  py -3 -u test/tabi_rsos_apbs_sweep.py --pqr "$2" --spacing "$3" --mode "$4" \
    --work-root "$5" > "$FINAL/$1.log" 2>&1
  echo "apbs $1 exit $? $(date)"
}
campaign proteins
apbs apbs_27 "$P27" 0.5 paper "$FINAL/apbs_27"
campaign binding
apbs apbs_binding "$BINDING" 0.5 paper "$FINAL/apbs_binding"
campaign hsweep_lowsalt
campaign hsweep
apbs apbs_sweep_lowsalt "$SWEEP" 0.6,0.5,0.4,0.3,0.25,0.2 lowsalt "$FINAL/apbs_sweep"
apbs apbs_sweep_paper "$SWEEP" 0.6,0.5,0.4,0.3,0.25,0.2 paper "$FINAL/apbs_sweep"
# Internal review F1: srfm mol gives enclosed cavities the solvent dielectric, while TITAN fills them
# (external-only); the same low-salt sweep with those cavities filled by uncharged probe-sized atoms
# (review_checks/apbs_cavity_box_check.py, 0.15 A grid, 26-connectivity) solves TITAN's problem.
# Same file names as the originals, so the APBS case names match; 1VII has no enclosed cavity.
CAV=$OUT/review_checks/cavfilled_pqr
apbs apbs_sweep_lowsalt_cavfilled "$CAV/1sh1_full_charmm_ph7.pqr,$CAV/1bpi_full_charmm_ph7.pqr,$CAV/1a63_full_charmm_ph7.pqr,$CAV/1a7m_full_charmm_ph7.pqr" \
  0.6,0.5,0.4,0.3,0.25,0.2 lowsalt "$FINAL/apbs_sweep_cavfilled"
campaign presets
campaign presets_lowsalt
campaign presets27
campaign threads
campaign xgen
campaign cap
campaign spacing
campaign coarse
campaign matrix "$SKIP"
campaign control451c
echo "timing block done $(date)"

# Reviewer 3: PAFMPB against TITAN on 1SH1 and 1A63 (shared protocol mesh and uniform SES meshes from
# afmpb_cmp/make_meshes.py), three rotated repeats, the same idle machine
cd "$OUT/afmpb_cmp" || exit 1
py -3 -u afmpb_bench.py "$BIN" D:/softs/pafmpb/AFMPB_win/afmpb.exe 3 "$OUT/afmpb_cmp/afmpb_bench.csv" \
  > "$OUT/afmpb_cmp/afmpb_bench.log" 2>&1
echo "afmpb bench exit $? $(date)"
cd "$REPO" || exit 1

R=$OUT
mv "$R/compare_idle" "$OLD/compare_idle" 2>/dev/null
py -3 -I test/tabi_rsos_compare.py --titan-lowsalt $MAIN/hsweep_lowsalt/hsweep_lowsalt.csv \
  --titan-lowsalt-ref $MAIN/hsweep_lowsalt_p4/hsweep_lowsalt_p4.csv \
  --apbs-lowsalt $FINAL/apbs_sweep/apbs_lowsalt.csv --titan-salt $MAIN/hsweep/hsweep.csv \
  --titan-salt-ref $MAIN/hsweep_p4/hsweep_p4.csv --apbs-salt $FINAL/apbs_sweep/apbs_paper.csv \
  --out $R/compare_idle > "$MAIN/compare_idle.log" 2>&1
echo "compare exit $? $(date)"
mv "$R/presets_curve_idle" "$OLD/presets_curve_idle" 2>/dev/null
T=$R/presets_curve_idle/inputs
mkdir -p $T/presets $T/presets_lowsalt $T/apbs_sweep
cp $MAIN/presets/presets.csv $T/presets/
cp $MAIN/presets_lowsalt/presets_lowsalt.csv $T/presets_lowsalt/
cp $FINAL/apbs_sweep/apbs_lowsalt.csv $FINAL/apbs_sweep/apbs_paper.csv $T/apbs_sweep/
py -3 -I test/tabi_rsos_presets_curve.py --final $T --campaign $MAIN --out $R/presets_curve_idle \
  > "$MAIN/presets_curve_idle.log" 2>&1
echo "presets curve exit $? $(date)"

kirkwood() {
  python -u test/tabi_node_patch_kirkwood.py --titan "$BIN" --work-root "$MAIN/$1" \
    --radius "$2" --offsets "$3" --h "$4" --variants "$5" > "$MAIN/$1.log" 2>&1
  echo "$1 exit $? $(date)"
}
kirkwood kirkwood_R10 10.0 0,0.5,0.9,0.95 1.25,1.0,0.75,0.5 paper_protocol
kirkwood kirkwood_R2 2.0 0,0.5,0.8 1.25,1.0,0.75,0.5,0.35,0.25 paper_protocol
kirkwood kirkwood_R10_p4 10.0 0,0.5,0.9,0.95 1.25,1.0,0.75,0.5 paper_protocol_p4
kirkwood kirkwood_R2_p4 2.0 0,0.5,0.8 1.25,1.0,0.75,0.5,0.35,0.25 paper_protocol_p4
echo "kirkwood done $(date)"

python -u test/tabi_rsos_campaign.py --resume --titan "$BIN" --work-root "$OUT/campaign_verify_p4" \
  --stages hsweep_lowsalt_p4,hsweep_p4 > "$OUT/campaign_verify_p4.log" 2>&1
echo "p4 sweeps exit $? $(date)"
echo "v06 timing done $(date)"
