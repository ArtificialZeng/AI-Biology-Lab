#!/bin/zsh
setopt nullglob
# Waits for the main agent's PREREG_D push (round11/status/PREREG_D_PUSHED.txt). Alerts on any round-11 download before R11_GO. 2 h limit.
B=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
start=$(date +%s)
while true; do
  [ -s $B/round11/status/PREREG_D_PUSHED.txt ] && { echo "$(date +%H:%M) PREREG_D_PUSHED: $(head -1 $B/round11/status/PREREG_D_PUSHED.txt)"; exit 0; }
  if [ ! -f $B/round11/status/R11_GO.txt ]; then
    n=$(find $B/round11/shared/raw $B/round11/shared/lines -type f 2>/dev/null | head -1)
    [ -n "$n" ] && { echo "$(date +%H:%M) R11_DOWNLOAD_BEFORE_GO_ALERT: $n"; exit 0; }
  fi
  [ $(( $(date +%s) - start )) -gt 7200 ] && { echo "$(date +%H:%M) TIMEOUT 2h"; exit 0; }
  sleep 60
done
