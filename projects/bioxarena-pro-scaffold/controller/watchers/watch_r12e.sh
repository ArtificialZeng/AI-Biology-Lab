#!/bin/bash
# bash 3.2. Evaluation phase. Fixed baselines (all 0 at 07:18). Exits on any listed marker not yet in $SEEN, frozen-file change, or 40-min stall.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cnt() { grep -cE "^[[:space:]\`*]*${2}[[:space:]\`*]*$" round12/status/$1 2>/dev/null; }
frozen_ok() { awk '/^# Section 3/{exit} !/^#/' round12/status/FREEZE_MANIFEST_R12.txt | shasum -a 256 -c --quiet >/dev/null 2>&1; }
T0=$(date +%s); SEEN=" ${SEEN:-} "
echo "$(date +%H:%M) start seen=[$SEEN]"
for i in $(seq 1 480); do
  now=$(date +%s)
  frozen_ok || { echo "$(date +%H:%M) FROZEN_FILE_CHANGED_ALERT"; exit 0; }
  for spec in A-T2.txt:R12_EVAL_DONE A-T3.txt:R12_AUDIT_DONE A-T3.txt:R12_TRIPLE_OK_A-T3 A-T3.txt:R12_TRIPLE_MISMATCH; do
    f=${spec%%:*}; m=${spec#*:}
    case "$SEEN" in *" $m "*) ;; *) [ "$(cnt $f $m)" -gt 0 ] && { echo "$(date +%H:%M) $m"; exit 0; };; esac
  done
  last=$( { stat -f %m round12/status/A-T2.txt round12/status/A-T3.txt; find round12/results round12/audit -newer round12/status/CODE_FREEZE_R12.txt -exec stat -f %m {} \; ; } 2>/dev/null | sort -n | tail -1)
  [ "$last" -lt "$T0" ] && last=$T0
  [ $(( now - last )) -gt 2400 ] && { echo "$(date +%H:%M) EVAL_STALL_ALERT"; exit 0; }
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT"
