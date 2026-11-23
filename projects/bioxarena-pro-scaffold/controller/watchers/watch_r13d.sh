#!/bin/bash
# bash 3.2. R13 download + reader-fix phase. Exits on a marker not in $SEEN, an integrity alert, or a 60-min stall.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cnt() { local n; n=$(grep -cE "^[[:space:]\`*]*${2}[[:space:]\`*]*$" round13/status/$1 2>/dev/null); echo ${n:-0}; }
frozen_ok() { awk '/^# Section 3/{exit} !/^#/' round12/status/FREEZE_MANIFEST_R12.txt | shasum -a 256 -c --quiet >/dev/null 2>&1; }
T0=$(date +%s); SEEN=" ${SEEN:-} "; echo "$(date +%H:%M) start seen=[$SEEN]"
for i in $(seq 1 480); do
  now=$(date +%s)
  frozen_ok || { echo "$(date +%H:%M) R12_FROZEN_FILE_CHANGED_ALERT"; exit 0; }
  [ "$(shasum -a 256 round13/prereg/PREREG_F.md | cut -c1-64)" = 0408ada0c0b88b9b2e1ddf97ba5f3e861dff0448f1bad36c6dd2347f5ac708ec ] || { echo "$(date +%H:%M) PREREG_F_CHANGED_ALERT"; exit 0; }
  [ ! -s round13/status/R13_LINES_AUTHORIZED.txt ] && [ -n "$(find round13/shared/lines -type f 2>/dev/null | head -1)" ] && { echo "$(date +%H:%M) LINES_BUILT_BEFORE_AUTH_ALERT"; exit 0; }
  big=$(find round13 -type f -size +100M -not -path "round13/shared/raw/gse327644/*" -not -path "round13/shared/lines/*" 2>/dev/null | head -1); [ -n "$big" ] && { echo "$(date +%H:%M) UNEXPECTED_BIG_FILE_ALERT $big"; exit 0; }
  for spec in A-T1.txt:R13_DOWNLOAD_DONE A-T2.txt:R13_CODE_PATCHED A-T3.txt:R13_AUDIT_CODE_READY A-T1.txt:R13_QC_DONE A-T3.txt:R13_READER_XCHECK_OK A-T3.txt:R13_READER_XCHECK_MISMATCH; do
    f=${spec%%:*}; m=${spec#*:}
    case "$SEEN" in *" $m "*) ;; *) [ "$(cnt $f $m)" -gt 0 ] && { echo "$(date +%H:%M) $m"; exit 0; };; esac
  done
  last=$( { find round13 -type f -exec stat -f %m {} \; ; } 2>/dev/null | sort -n | tail -1); [ -z "$last" ] && last=$T0
  [ "$last" -lt "$T0" ] && last=$T0
  [ $(( now - last )) -gt 3600 ] && { echo "$(date +%H:%M) R13_STALL_ALERT"; exit 0; }
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT"
