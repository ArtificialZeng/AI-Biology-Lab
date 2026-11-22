#!/bin/bash
# bash 3.2. Prints each marker once; exits when the R12 data screen is done (or frozen R11 files change). 4 h timeout.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
seen=""
for i in $(seq 1 480); do
  if ! shasum -a 256 -c round11/status/FREEZE_MANIFEST_R11.txt >/dev/null 2>&1; then echo "$(date +%H:%M) FROZEN_R11_FILE_CHANGED_ALERT"; exit 0; fi
  for pair in main.txt:MEDIA_BRIEF_REMOVED main.txt:R11_REVIEW_D_PUSHED A-T4.txt:OUTLINE_EXPLORATORY_DONE A-T1.txt:R12_SCREEN_DONE; do
    f=${pair%%:*}; m=${pair##*:}
    case " $seen " in *" $m "*) ;; *) if grep -q "$m" round11/status/$f 2>/dev/null; then seen="$seen $m"; echo "$(date +%H:%M) $m"; fi;; esac
  done
  case "$seen" in *R12_SCREEN_DONE*) echo "$(date +%H:%M) SCREEN_DONE:$seen"; exit 0;; esac
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT:$seen"
