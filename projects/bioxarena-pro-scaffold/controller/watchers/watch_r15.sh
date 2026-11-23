#!/bin/bash
# bash 3.2. Round 15. MARKERS="file:MARK ..." (files under round15/status). SEEN = markers to ignore.
# Alerts: R12-R14 frozen files changed; >100MB file under round15; paper/ in git; VCC 2026 bucket non-empty; 60-min stall.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cnt() { local n; n=$(grep -cE "^[[:space:]\`*]*${2}[[:space:]\`*]*$" round15/status/$1 2>/dev/null); echo ${n:-0}; }
frozen_ok() { for m in round12/status/FREEZE_MANIFEST_R12.txt round13/status/FREEZE_MANIFEST_R13.txt; do awk '/^# Section 3/{exit} !/^#/' $m | shasum -a 256 -c --quiet >/dev/null 2>&1 || return 1; done; return 0; }
gcs() { curl -s --max-time 20 "https://storage.googleapis.com/storage/v1/b/arc-institute-virtual-cell-atlas/o?prefix=virtual-cell-challenge/2026/&maxResults=5&fields=items(name)" | grep -c '"name"'; }
T0=$(date +%s); SEEN=" ${SEEN:-} "; echo "$(date +%H:%M) start markers=[$MARKERS] seen=[$SEEN]"
for i in $(seq 1 720); do
  now=$(date +%s)
  frozen_ok || { echo "$(date +%H:%M) FROZEN_FILE_CHANGED_ALERT"; exit 0; }
  big=$(find round15 -type f -size +100M 2>/dev/null | head -1); [ -n "$big" ] && { echo "$(date +%H:%M) BIG_FILE_ALERT $big"; exit 0; }
  [ -n "$(git -c http.proxy= -c https.proxy= ls-files paper 2>/dev/null | head -1)$(git -c http.proxy= -c https.proxy= diff --cached --name-only -- paper 2>/dev/null | head -1)" ] && { echo "$(date +%H:%M) PAPER_IN_GIT_ALERT"; exit 0; }
  [ $((i % 20)) -eq 1 ] && [ "$(gcs)" -gt 0 ] && { echo "$(date +%H:%M) VCC2026_BUCKET_NONEMPTY_ALERT"; exit 0; }
  for spec in $MARKERS; do
    f=${spec%%:*}; m=${spec#*:}
    case "$SEEN" in *" $m "*) ;; *) [ "$(cnt $f $m)" -gt 0 ] && { echo "$(date +%H:%M) $m"; exit 0; };; esac
  done
  last=$( { find round15 -type f -exec stat -f %m {} \; ; } 2>/dev/null | sort -n | tail -1); [ -z "$last" ] && last=$T0; [ "$last" -lt "$T0" ] && last=$T0
  [ $(( now - last )) -gt 3600 ] && { echo "$(date +%H:%M) R15_STALL_ALERT"; exit 0; }
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT"
