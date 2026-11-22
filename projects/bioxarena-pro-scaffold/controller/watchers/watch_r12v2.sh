#!/bin/bash
# bash 3.2. Exits when PREREG_E v2 and plan v2 are both ready; alerts on big data files in round12 before PREREG_E_COMMIT. 4 h timeout.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
seen=""
for i in $(seq 1 480); do
  if [ ! -s round12/status/PREREG_E_COMMIT.txt ]; then
    big=$(find round12 -type f \( -name '*.mtx*' -o -name '*.h5ad' -o -name '*.h5' \) -size +50M 2>/dev/null | head -1)
    if [ -n "$big" ]; then echo "$(date +%H:%M) DATA_BEFORE_LOCK_ALERT $big"; exit 0; fi
  fi
  for pair in A-T3.txt:PAIRED_RORTH_DONE A-T3.txt:PREREG_E_V2_READY A-T1.txt:R12_PLAN_V2_READY A-T2.txt:R12_CODE_DRAFT_READY A-T5.txt:PRIOR_ART_F_DONE; do
    f=${pair%%:*}; m=${pair##*:}
    case " $seen " in *" $m "*) ;; *) if grep -q "$m" round12/status/$f 2>/dev/null; then seen="$seen $m"; echo "$(date +%H:%M) $m"; fi;; esac
  done
  case "$seen" in *PREREG_E_V2_READY*) case "$seen" in *R12_PLAN_V2_READY*) echo "$(date +%H:%M) V2_READY:$seen"; exit 0;; esac;; esac
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT:$seen"
