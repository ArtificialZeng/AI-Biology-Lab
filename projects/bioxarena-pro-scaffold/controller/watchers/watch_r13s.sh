#!/bin/bash
# bash 3.2. Round-13 screen + prior art. Fixed baselines (0 at 07:53). Exits on a marker not in $SEEN, R12 frozen-file change, big download, or 60-min stall.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cnt() { local n; n=$(grep -cE "^[[:space:]\`*]*${2}[[:space:]\`*]*$" round13/status/$1 2>/dev/null); echo ${n:-0}; }
frozen_ok() { awk '/^# Section 3/{exit} !/^#/' round12/status/FREEZE_MANIFEST_R12.txt | shasum -a 256 -c --quiet >/dev/null 2>&1; }
T0=$(date +%s); SEEN=" ${SEEN:-} "; echo "$(date +%H:%M) start seen=[$SEEN]"
for i in $(seq 1 480); do
  now=$(date +%s)
  frozen_ok || { echo "$(date +%H:%M) R12_FROZEN_FILE_CHANGED_ALERT"; exit 0; }
  big=$(find round13 -type f -size +100M 2>/dev/null | head -1); [ -n "$big" ] && { echo "$(date +%H:%M) BIG_DOWNLOAD_ALERT $big"; exit 0; }
  for spec in A-T1.txt:R13_SCREEN_DONE A-T5.txt:R13_PRIOR_ART_DONE A-T5.txt:R13_PRIOR_ART_V2_DONE A-T3.txt:PREREG_F_DRAFT_READY A-T2.txt:R13_CODE_READY A-T4.txt:R13_OUTLINE_DONE; do
    f=${spec%%:*}; m=${spec#*:}
    case "$SEEN" in *" $m "*) ;; *) [ "$(cnt $f $m)" -gt 0 ] && { echo "$(date +%H:%M) $m"; exit 0; };; esac
  done
  last=$( { find round13 -type f -exec stat -f %m {} \; ; } 2>/dev/null | sort -n | tail -1); [ -z "$last" ] && last=$T0
  [ "$last" -lt "$T0" ] && last=$T0
  [ $(( now - last )) -gt 3600 ] && { echo "$(date +%H:%M) R13_STALL_ALERT"; exit 0; }
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT"
