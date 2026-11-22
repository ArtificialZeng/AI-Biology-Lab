#!/bin/bash
# bash 3.2. R12 data phase: prints markers once; alerts if the locked prereg/plan changes, if download exceeds 15 GB, or if any round12/results file appears before CODE_FREEZE_R12. Exits when QC done AND both code-ready markers exist. 6 h timeout.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
seen=""
for i in $(seq 1 720); do
  [ "$(shasum -a 256 round12/prereg/PREREG_E.md | cut -c1-64)" = d10f92ced16a4519c43ea2b92b8b7ee881ab3d1764e39ffd10580473e4108924 ] || { echo "$(date +%H:%M) PREREG_E_CHANGED_ALERT"; exit 0; }
  [ "$(shasum -a 256 round12/shared/R12_DOWNLOAD_PLAN.md | cut -c1-64)" = e3826e09dd4326b1724d25732bb05d3d4dc72681ffc283b4a5efba57bd1555ed ] || { echo "$(date +%H:%M) PLAN_CHANGED_ALERT"; exit 0; }
  kb=$(du -sk round12/shared/raw 2>/dev/null | cut -f1); [ -n "$kb" ] && [ "$kb" -gt 14700000 ] && { echo "$(date +%H:%M) DOWNLOAD_OVER_CAP_ALERT ${kb}KB"; exit 0; }
  if [ ! -s round12/status/CODE_FREEZE_R12.txt ] && [ -n "$(find round12/results -type f -not -path "*/dev_*" 2>/dev/null | head -1)" ]; then echo "$(date +%H:%M) RESULTS_BEFORE_FREEZE_ALERT"; exit 0; fi
  [ -s round12/prereg/DEVIATIONS_R12.md ] && case " $seen " in *" DEV "*) ;; *) seen="$seen DEV"; echo "$(date +%H:%M) DEVIATIONS_R12_WRITTEN";; esac
  for pair in A-T1.txt:R12_QC_DONE A-T2.txt:R12_CODE_READY A-T3.txt:R12_AUDIT_CODE_READY; do
    f=${pair%%:*}; m=${pair##*:}
    case " $seen " in *" $m "*) ;; *) if grep -q "$m" round12/status/$f 2>/dev/null; then seen="$seen $m"; echo "$(date +%H:%M) $m"; fi;; esac
  done
  case "$seen" in *R12_QC_DONE*) case "$seen" in *R12_CODE_READY*) case "$seen" in *R12_AUDIT_CODE_READY*) echo "$(date +%H:%M) READY_FOR_FREEZE:$seen"; exit 0;; esac;; esac;; esac
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT:$seen"
