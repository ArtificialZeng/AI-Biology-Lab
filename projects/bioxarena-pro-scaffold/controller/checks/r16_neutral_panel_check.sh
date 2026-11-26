#!/bin/bash
# Controller: independent recompute of the PREREG_I §4.4 neutral-batch robustness panel (STATE_SE; main line only,
# not in A-T3's pass 2). For each line and mode in {bypass_zero, mean_batch_emb, multi_batch_ens8}, the main line's
# exported arrays are linked into the layout r16_check.py reads (mode arrays under the fixed_index_zero names), scored
# by the unchanged controller code, and compared with summary_r16.json lines/<line>/secondary/neutral_batch_robustness_panel.
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cd "$P" || exit 9
C=controller/checks; OUT=$C/r16_results/neutral_panel; mkdir -p $OUT
for ln in nadig_hepg2 nadig_jurkat a549_ko dld1 cd4t_rest; do
  seen=$(python3 -I -c "import json; print(','.join(json.load(open('round16/results/summary_r16.json'))['lines']['$ln']['strata_sources']['seen_by_checkpoint']))")
  unseen=$(python3 -I -c "import json; print(','.join(json.load(open('round16/results/summary_r16.json'))['lines']['$ln']['strata_sources']['unseen_by_checkpoint']))")
  for mode in bypass_zero mean_batch_emb multi_batch_ens8; do
    L=$OUT/layout/$ln/$mode; rm -rf $L; mkdir -p $L
    ln -s "$P/round16/results/arrays/$ln/observed.npz" $L/observed.npz
    for f in round16/results/arrays/$ln/preds_STATE_SE_${mode}_*.npz; do
      b=$(basename $f); ln -s "$P/$f" "$L/${b/preds_STATE_SE_${mode}_/preds_STATE_SE_fixed_index_zero_}"
    done
    python3 -I $C/r16_check.py $L "$seen" "$unseen" $OUT/${ln}__${mode}.json 4000 20261601 >/dev/null || echo "r16_check failed $ln $mode"
  done
done
python3 -I - <<'EOF'
import json
S = json.load(open('round16/results/summary_r16.json'))['lines']
worst, bad, n = 0.0, [], 0
def num(tag, a, b):
    global worst, n
    n += 1; d = abs(float(a) - float(b)); worst = max(worst, d)
    if not d <= 1e-5: bad.append(f'{tag}: main={a} controller={b}')
for ln, L in S.items():
    nb = L['secondary']['neutral_batch_robustness_panel']
    for mode, v in nb.items():
        c = json.load(open(f'controller/checks/r16_results/neutral_panel/{ln}__{mode}.json'))
        m = v['matched_orth']
        num(f'{ln}.{mode}.matched.est', m['estimate'], c['arms']['matched']['point'])
        num(f'{ln}.{mode}.matched.lo', m['ci975_lo'], c['arms']['matched']['ci975'][0]); num(f'{ln}.{mode}.matched.hi', m['ci975_hi'], c['arms']['matched']['ci975'][1])
        f = v['vs_fct']
        num(f'{ln}.{mode}.vsfct.est', f['estimate'], c['H_I2']['point'])
        num(f'{ln}.{mode}.vsfct.lo', f['ci975_lo'], c['H_I2']['ci975'][0]); num(f'{ln}.{mode}.vsfct.hi', f['ci975_hi'], c['H_I2']['ci975'][1])
        for s, ps in v['paired_vs_mismatched']['per_source'].items():
            num(f'{ln}.{mode}.paired.{s}.est', ps['estimate'], c['pairs'][s]['point'])
            num(f'{ln}.{mode}.paired.{s}.lo', ps['ci975_lo'], c['pairs'][s]['ci975'][0])
        print(f"{ln:13s} {mode:16s} matched {c['arms']['matched']['point']:.4f}  STATE-FCT {c['H_I2']['point']:+.4f} [{c['H_I2']['ci975'][0]:+.4f},{c['H_I2']['ci975'][1]:+.4f}] {c['H_I2']['call']}")
print(json.dumps({'n_checks': n, 'n_fail': len(bad), 'worst': worst}))
for b in bad[:10]: print('FAIL', b)
print('R16_NEUTRAL_PANEL_MATCH' if not bad else 'R16_NEUTRAL_PANEL_DIFF')
EOF
