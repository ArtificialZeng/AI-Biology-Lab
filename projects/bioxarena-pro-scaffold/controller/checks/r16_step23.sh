#!/bin/bash
# Controller: staged verification steps 2 and 3 on A-T3 pass 1 (Groups 1-2), run after R16_AT3_REAL_DONE.
#  step 2: A-T3 arrays vs main-line arrays (r16_compare_arrays.py; predicted Δ <= 1e-4, observed side <= 1e-5)
#  step 3: controller recompute on A-T3's arrays (r16_check.py, same code as for the main line) compared with the
#          main line's summary_r16.json (r16_compare_main.py; every estimate and CI bound <= 1e-5, calls identical)
# A-T3's own summary is compared separately once its schema is known.
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cd "$P" || exit 9
C=controller/checks; R=$C/r16_results; mkdir -p $R/at3 $R/at3_layout
fail=0
for ln in nadig_hepg2 nadig_jurkat a549_ko dld1 cd4t_rest; do
  for m in STATE_SE STATE_HVG; do
    f=$(ls round16/audit/results/${m}_*__${ln}.npz 2>/dev/null | head -1)
    [ -z "$f" ] && { echo "MISSING A-T3 ${m} file for ${ln}"; fail=1; continue; }
    python3 -I $C/r16_compare_arrays.py "$f" round16/results/arrays/$ln | grep -E "FAIL|R16_ARRAYS" || fail=1
    python3 -I $C/r16_compare_arrays.py "$f" round16/results/arrays/$ln | grep -q "^R16_ARRAYS_MATCH" || fail=1
  done
  f=$(ls round16/audit/results/STATE_SE_*__${ln}.npz 2>/dev/null | head -1); [ -z "$f" ] && continue
  meta=$(python3 -I $C/r16_at3_to_check_layout.py "$f" $R/at3_layout/$ln)
  seen=$(echo "$meta" | python3 -I -c "import json,sys; print(','.join(json.load(sys.stdin)['seen']))")
  unseen=$(echo "$meta" | python3 -I -c "import json,sys; print(','.join(json.load(sys.stdin)['unseen']))")
  echo "$ln seen=$seen unseen=$unseen"
  python3 -I $C/r16_check.py $R/at3_layout/$ln "$seen" "$unseen" $R/at3/$ln.json 4000 20261601 >/dev/null || { echo "r16_check failed on $ln"; fail=1; }
done
python3 -I $C/r16_compare_main.py round16/results/summary_r16.json $R/at3 | tail -12
python3 -I $C/r16_compare_main.py round16/results/summary_r16.json $R/at3 | grep -q "^R16_MAIN_VS_CONTROLLER_MATCH" || fail=1
[ $fail = 0 ] && echo "R16_STEP23_PASS" || echo "R16_STEP23_FAIL"
