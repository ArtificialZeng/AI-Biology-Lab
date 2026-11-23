#!/bin/bash
# bash 3.2. Waits for the main agent's freeze commit. Exits on marker, big-data-in-git alert, results-before-freeze, or 20-min stall.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
G="git -c http.proxy= -c https.proxy="
T0=$(date +%s); H0=$($G rev-parse HEAD)
echo "$(date +%H:%M) start HEAD=$H0"
for i in $(seq 1 240); do
  now=$(date +%s)
  [ -n "$($G ls-files round12/shared/raw round12/shared/lines round12/audit/reader_xcheck round13/shared/raw | head -1)" ] && { echo "$(date +%H:%M) BIG_DATA_STAGED_OR_COMMITTED_ALERT"; exit 0; }
  [ -n "$($G diff --cached --name-only -- round12/shared/raw round12/shared/lines round12/audit/reader_xcheck round13/shared/raw | head -1)" ] && { echo "$(date +%H:%M) BIG_DATA_STAGED_OR_COMMITTED_ALERT"; exit 0; }
  if [ ! -s round12/status/CODE_FREEZE_R12.txt ] && [ -n "$(find round12/results -type f -not -path "*/dev_*" 2>/dev/null | head -1)" ]; then echo "$(date +%H:%M) RESULTS_BEFORE_FREEZE_ALERT"; exit 0; fi
  grep -qE "^[[:space:]\`*]*${MARK:-R12_FREEZE_COMMITTED}[[:space:]\`*]*$" ${MAINF:-round12/status/main.txt} 2>/dev/null && { echo "$(date +%H:%M) ${MARK:-R12_FREEZE_COMMITTED} HEAD=$($G rev-parse HEAD)"; exit 0; }
  last=$(stat -f %m ${MAINF:-round12/status/main.txt}); [ "$($G rev-parse HEAD)" != "$H0" ] && last=$now
  [ "$last" -lt "$T0" ] && last=$T0
  [ $(( now - last )) -gt 1200 ] && { echo "$(date +%H:%M) MAIN_STALL_ALERT HEAD=$($G rev-parse HEAD)"; exit 0; }
  sleep 20
done
echo "$(date +%H:%M) TIMEOUT"
