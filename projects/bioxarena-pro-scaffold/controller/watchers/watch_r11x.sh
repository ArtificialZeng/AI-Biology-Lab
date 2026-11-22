#!/bin/bash
# bash 3.2 compatible. Prints each done marker once; exits when A-T5 and A-T3 are both done, or on frozen-file change. 4 h timeout.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
seen=""
for i in $(seq 1 480); do
  if ! shasum -a 256 -c round11/status/FREEZE_MANIFEST_R11.txt >/dev/null 2>&1; then echo "$(date +%H:%M) FROZEN_R11_FILE_CHANGED_ALERT"; exit 0; fi
  for pair in main.txt:R11_VERDICT_PUSHED A-T5.txt:PRIOR_ART_E_DONE A-T3.txt:EXPLORATORY_CHECK_DONE A-T4.txt:OUTLINE_R11_DONE; do
    f=${pair%%:*}; m=${pair##*:}
    case " $seen " in *" $m "*) ;; *) if grep -q "$m" round11/status/$f 2>/dev/null; then seen="$seen $m"; echo "$(date +%H:%M) $m"; fi;; esac
  done
  case "$seen" in *PRIOR_ART_E_DONE*) case "$seen" in *EXPLORATORY_CHECK_DONE*) echo "$(date +%H:%M) KEY_TASKS_DONE:$seen"; exit 0;; esac;; esac
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT:$seen"
