#!/bin/bash
# bash 3.2. Fix for watch_r12u race: baselines are FIXED (known state at 07:08), not sampled at start.
# Exits on: R12_QC_DONE2 (A-T1), R12_READER_XCHECK_OK / _MISMATCH (A-T3), any alert, or A-T1 stall (30 min no activity).
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
sig() { find round12/shared/lines/a549_ko round12/shared/lines/a549_ko_cas9pos -type f -exec stat -f "%N %z %m" {} \; 2>/dev/null | sort | shasum -a 256 | cut -c1-64; }
cnt() { grep -cE "^[[:space:]\`*]*${2}[[:space:]\`*]*$" round12/status/$1 2>/dev/null; }
A0=$(sig); T0=$(date +%s)
B_QC2=${B_QC2:-0}; B_OK=${B_OK:-0}; B_MM=${B_MM:-1}
echo "$(date +%H:%M) start baselines QC2=$B_QC2 OK=$B_OK MM=$B_MM"
for i in $(seq 1 360); do
  now=$(date +%s)
  [ "$(shasum -a 256 round12/prereg/PREREG_E.md | cut -c1-64)" = d10f92ced16a4519c43ea2b92b8b7ee881ab3d1764e39ffd10580473e4108924 ] || { echo "$(date +%H:%M) PREREG_E_CHANGED_ALERT"; exit 0; }
  [ "$(sig)" = "$A0" ] || { echo "$(date +%H:%M) A549_LINES_CHANGED_ALERT"; exit 0; }
  if [ ! -s round12/status/CODE_FREEZE_R12.txt ] && [ -n "$(find round12/results -type f -not -path "*/dev_*" 2>/dev/null | head -1)" ]; then echo "$(date +%H:%M) RESULTS_BEFORE_FREEZE_ALERT"; exit 0; fi
  [ "$(cnt A-T3.txt R12_READER_XCHECK_OK)" -gt "$B_OK" ] && { echo "$(date +%H:%M) R12_READER_XCHECK_OK"; exit 0; }
  [ "$(cnt A-T3.txt R12_READER_XCHECK_MISMATCH)" -gt "$B_MM" ] && { echo "$(date +%H:%M) R12_READER_XCHECK_MISMATCH"; exit 0; }
  [ "$(cnt A-T1.txt R12_QC_DONE2)" -gt "$B_QC2" ] && { echo "$(date +%H:%M) R12_QC_DONE2"; exit 0; }
  # activity = newest mtime among the watched agent's status file and the HEK293 line dir
  WATCH=${WATCH:-A-T1}
  if [ "$WATCH" = A-T1 ]; then last=$( { stat -f %m round12/status/A-T1.txt; find round12/shared/lines/hek293_xucao round12/shared -maxdepth 1 -newer round12/status/A-T2.txt -exec stat -f %m {} \; ; } 2>/dev/null | sort -n | tail -1)
  else last=$( { stat -f %m round12/status/A-T3.txt; find round12/audit -newer round12/status/A-T1.txt -exec stat -f %m {} \; ; } 2>/dev/null | sort -n | tail -1); fi
  [ "$last" -lt "$T0" ] && last=$T0
  [ $(( now - last )) -gt 1800 ] && { echo "$(date +%H:%M) ${WATCH}_STALL_ALERT"; exit 0; }
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT"
