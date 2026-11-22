#!/bin/bash
# bash 3.2. R12 pre-freeze fixes: markers must stand alone on a line (avoids the 01:17 false trigger).
# Alerts: prereg/plan hash change, raw > 14.7 GB, non-dev results before CODE_FREEZE_R12, A549 line files changed (rebuild forbidden).
# Exits READY_FOR_FREEZE when R12_CODE_PATCHED, R12_AUDIT_PATCHED, R12_QC_DONE all present. 3 h timeout.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
sig() { find round12/shared/lines/a549_ko round12/shared/lines/a549_ko_cas9pos -type f -exec stat -f "%N %z %m" {} \; 2>/dev/null | sort | shasum -a 256 | cut -c1-64; }
A549SIG=$(sig); echo "$(date +%H:%M) start a549sig=${A549SIG:0:12}"
seen=""
for i in $(seq 1 360); do
  [ "$(shasum -a 256 round12/prereg/PREREG_E.md | cut -c1-64)" = d10f92ced16a4519c43ea2b92b8b7ee881ab3d1764e39ffd10580473e4108924 ] || { echo "$(date +%H:%M) PREREG_E_CHANGED_ALERT"; exit 0; }
  [ "$(shasum -a 256 round12/shared/R12_DOWNLOAD_PLAN.md | cut -c1-64)" = e3826e09dd4326b1724d25732bb05d3d4dc72681ffc283b4a5efba57bd1555ed ] || { echo "$(date +%H:%M) PLAN_CHANGED_ALERT"; exit 0; }
  kb=$(du -sk round12/shared/raw 2>/dev/null | cut -f1); [ -n "$kb" ] && [ "$kb" -gt 14700000 ] && { echo "$(date +%H:%M) DOWNLOAD_OVER_CAP_ALERT ${kb}KB"; exit 0; }
  if [ ! -s round12/status/CODE_FREEZE_R12.txt ] && [ -n "$(find round12/results -type f -not -path "*/dev_*" 2>/dev/null | head -1)" ]; then echo "$(date +%H:%M) RESULTS_BEFORE_FREEZE_ALERT"; exit 0; fi
  [ "$(sig)" = "$A549SIG" ] || { echo "$(date +%H:%M) A549_LINES_CHANGED_ALERT"; exit 0; }
  [ -s round12/prereg/DEVIATIONS_R12.md ] && case " $seen " in *" DEV "*) ;; *) seen="$seen DEV"; echo "$(date +%H:%M) DEVIATIONS_R12_WRITTEN";; esac
  [ -s round12/status/CONFIRMATORY_LINES_R12.txt ] && case " $seen " in *" CONF "*) ;; *) seen="$seen CONF"; echo "$(date +%H:%M) CONFIRMATORY_LINES_R12_WRITTEN";; esac
  for pair in A-T2.txt:R12_CODE_PATCHED A-T3.txt:R12_AUDIT_PATCHED A-T1.txt:R12_QC_DONE; do
    f=${pair%%:*}; m=${pair##*:}
    case " $seen " in *" $m "*) ;; *) if grep -qE "^[[:space:]\`*]*${m}[[:space:]\`*]*$" round12/status/$f 2>/dev/null; then seen="$seen $m"; echo "$(date +%H:%M) $m"; fi;; esac
  done
  case "$seen" in *R12_CODE_PATCHED*) case "$seen" in *R12_AUDIT_PATCHED*) case "$seen" in *R12_QC_DONE*) echo "$(date +%H:%M) READY_FOR_FREEZE:$seen"; exit 0;; esac;; esac;; esac
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT:$seen"
