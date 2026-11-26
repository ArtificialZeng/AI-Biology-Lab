#!/bin/bash
# Controller: verification of A-T3 pass 2 (round16/audit/results/secondary/), run after R16_AT3_SECONDARY_DONE.
#  a) pass-2 arrays (both models) vs main-line arrays: same check as step 2 (r16_compare_arrays.py)
#  b) A-T3's own STATE_HVG metrics vs main summary (r16_compare_at3_metrics.py)
#  c) G4 decomposition, symmetry check, single-index spread (STATE_SE) and HVG Perm100 null (r16_compare_at3_pass2.py)
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cd "$P" || exit 9
C=controller/checks; A=round16/audit/results/secondary; fail=0
for ln in nadig_hepg2 nadig_jurkat a549_ko dld1 cd4t_rest; do
  for m in STATE_SE STATE_HVG; do
    f=$(ls $A/${m}_*__${ln}.npz 2>/dev/null | head -1)
    [ -z "$f" ] && { echo "MISSING pass-2 ${m} file for ${ln}"; fail=1; continue; }
    out=$(python3 -I $C/r16_compare_arrays.py "$f" round16/results/arrays/$ln)
    echo "$out" | grep -E "FAIL|^R16_ARRAYS"
    echo "$out" | grep -q "^R16_ARRAYS_MATCH" || fail=1
  done
done
out=$(python3 -I $C/r16_compare_at3_metrics.py round16/results/summary_r16.json $A STATE_HVG)
echo "$out" | tail -4; echo "$out" | grep -q "^R16_AT3_METRICS_MATCH" || fail=1
out=$(python3 -I $C/r16_compare_at3_pass2.py round16/results/summary_r16.json $A)
echo "$out" | head -1; echo "$out" | grep -E "^(NOTE|FAIL)" | head -20; echo "$out" | tail -1
echo "$out" | grep -q "^R16_AT3_PASS2_MATCH" || fail=1
[ $fail = 0 ] && echo "R16_PASS2_CHECK_PASS" || echo "R16_PASS2_CHECK_FAIL"
