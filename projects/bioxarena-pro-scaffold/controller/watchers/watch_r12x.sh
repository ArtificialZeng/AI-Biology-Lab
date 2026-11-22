#!/bin/bash
# bash 3.2. Waits for A-T3's reader cross-check (whole-line marker). Alerts on prereg change, shared line data change, results before freeze. 2.5 h timeout.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
sig() { find round12/shared/lines -type f -not -path "*/.cache*" -exec stat -f "%N %z %m" {} \; 2>/dev/null | sort | shasum -a 256 | cut -c1-64; }
S0=$(sig); echo "$(date +%H:%M) start linesig=${S0:0:12}"
for i in $(seq 1 300); do
  [ "$(shasum -a 256 round12/prereg/PREREG_E.md | cut -c1-64)" = d10f92ced16a4519c43ea2b92b8b7ee881ab3d1764e39ffd10580473e4108924 ] || { echo "$(date +%H:%M) PREREG_E_CHANGED_ALERT"; exit 0; }
  if [ ! -s round12/status/CODE_FREEZE_R12.txt ] && [ -n "$(find round12/results -type f -not -path "*/dev_*" 2>/dev/null | head -1)" ]; then echo "$(date +%H:%M) RESULTS_BEFORE_FREEZE_ALERT"; exit 0; fi
  [ "$(sig)" = "$S0" ] || { echo "$(date +%H:%M) SHARED_LINES_CHANGED_ALERT"; exit 0; }
  for m in R12_READER_XCHECK_OK R12_READER_XCHECK_MISMATCH; do
    grep -qE "^[[:space:]\`*]*${m}[[:space:]\`*]*$" round12/status/A-T3.txt 2>/dev/null && { echo "$(date +%H:%M) $m"; exit 0; }
  done
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT"
