#!/bin/bash
# bash 3.2. Prints each R12 marker once; exits when the PREREG_E draft is ready AND the download plan is ready; alerts if any big data file appears in round12 before a PREREG_E_COMMIT exists. 4 h timeout.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
seen=""
for i in $(seq 1 480); do
  if [ ! -s round12/status/PREREG_E_COMMIT.txt ]; then
    big=$(find round12 -type f \( -name '*.mtx*' -o -name '*.h5ad' -o -name '*.h5' \) -size +50M 2>/dev/null | head -1)
    if [ -n "$big" ]; then echo "$(date +%H:%M) DATA_BEFORE_LOCK_ALERT $big"; exit 0; fi
  fi
  for pair in A-T3.txt:PREREG_E_DRAFT_READY A-T1.txt:R12_PLAN_READY A-T2.txt:R12_CODE_DRAFT_READY A-T5.txt:PRIOR_ART_F_DONE; do
    f=${pair%%:*}; m=${pair##*:}
    case " $seen " in *" $m "*) ;; *) if grep -q "$m" round12/status/$f 2>/dev/null; then seen="$seen $m"; echo "$(date +%H:%M) $m"; fi;; esac
  done
  case "$seen" in *PREREG_E_DRAFT_READY*) case "$seen" in *R12_PLAN_READY*) echo "$(date +%H:%M) PREREG_INPUTS_READY:$seen"; exit 0;; esac;; esac
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT:$seen"
