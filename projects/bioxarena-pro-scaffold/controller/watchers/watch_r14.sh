#!/bin/bash
# bash 3.2. Round 14. MARKERS="file:MARK file:MARK ..." (files under round14/status). SEEN = markers to ignore.
# Alerts: R12/R13 frozen code/docs changed; big file under round14 outside authorized dirs; 60-min stall.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cnt() { local n; n=$(grep -cE "^[[:space:]\`*]*${2}[[:space:]\`*]*$" round14/status/$1 2>/dev/null); echo ${n:-0}; }
frozen_ok() { for m in round12/status/FREEZE_MANIFEST_R12.txt round13/status/FREEZE_MANIFEST_R13.txt; do awk '/^# Section 3/{exit} !/^#/' $m | shasum -a 256 -c --quiet >/dev/null 2>&1 || return 1; done; return 0; }
T0=$(date +%s); SEEN=" ${SEEN:-} "; echo "$(date +%H:%M) start markers=[$MARKERS] seen=[$SEEN]"
for i in $(seq 1 480); do
  now=$(date +%s)
  frozen_ok || { echo "$(date +%H:%M) FROZEN_FILE_CHANGED_ALERT"; exit 0; }
  if [ -n "$PREG_SHA" ]; then [ "$(shasum -a 256 round14/prereg/PREREG_G.md | cut -c1-64)" = "$PREG_SHA" ] || { echo "$(date +%H:%M) PREREG_G_CHANGED_ALERT"; exit 0; }; fi
  [ -s round14/status/R14_GO.txt ] || { big=$(find round14 -type f -size +100M 2>/dev/null | head -1); [ -n "$big" ] && { echo "$(date +%H:%M) BIG_FILE_BEFORE_GO_ALERT $big"; exit 0; }; }
  [ ! -s round14/status/R14_LINES_AUTHORIZED.txt ] && [ -n "$(find round14/shared/lines -type f 2>/dev/null | head -1)" ] && { echo "$(date +%H:%M) LINES_BUILT_BEFORE_AUTH_ALERT"; exit 0; }
  for spec in $MARKERS; do
    f=${spec%%:*}; m=${spec#*:}
    case "$SEEN" in *" $m "*) ;; *) [ "$(cnt $f $m)" -gt 0 ] && { echo "$(date +%H:%M) $m"; exit 0; };; esac
  done
  last=$( { find round14 -type f -exec stat -f %m {} \; ; } 2>/dev/null | sort -n | tail -1); [ -z "$last" ] && last=$T0; [ "$last" -lt "$T0" ] && last=$T0
  [ $(( now - last )) -gt 3600 ] && { echo "$(date +%H:%M) R14_STALL_ALERT"; exit 0; }
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT"
