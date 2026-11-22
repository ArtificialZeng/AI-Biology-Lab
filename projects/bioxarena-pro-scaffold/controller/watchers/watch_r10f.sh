#!/bin/zsh
setopt nullglob
# Waits for the main agent's freeze-point push (round10/status/CODE_FREEZE_PUSHED.txt). Alerts on early scoring. 1 h limit.
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10
start=$(date +%s)
while true; do
  grep -q QC_LOCK $P/status/CODE_FREEZE_PUSHED.txt 2>/dev/null && { echo "$(date +%H:%M) QC_LOCK_PUSHED: $(grep QC_LOCK $P/status/CODE_FREEZE_PUSHED.txt)"; exit 0; }
  early=$(find $P/results $P/audit -type f -newer $P/status/CONFIRMATORY_LINES.txt 2>/dev/null | grep -v -E "__pycache__|\.py$|r9dev|r9_dev|_dev_" | head -3)
  [ -n "$early" ] && { echo "$(date +%H:%M) RESULTS_BEFORE_FREEZE_ALERT: $early"; exit 0; }
  [ $(( $(date +%s) - start )) -gt 3600 ] && { echo "$(date +%H:%M) TIMEOUT 1h"; exit 0; }
  sleep 60
done
