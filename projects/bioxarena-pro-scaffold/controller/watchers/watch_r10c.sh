#!/bin/zsh
# Round-10 phase-c watcher: wait for the PREREG_C push; also catches early downloads, the paper revision, and data tracked.
C=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller
P9=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10
R=/Users/mac/AI-Biology-Lab
start=$(date +%s); head0=$(git -C $R rev-parse HEAD)
while true; do
  [ -f $P/status/PREREG_C_PUSHED.txt ] && { echo "$(date +%H:%M) PREREG_C_PUSHED"; exit 0; }
  git -C $R log --name-only --format= $head0..HEAD 2>/dev/null | grep -q "round10/prereg/PREREG_C.md" && { echo "$(date +%H:%M) PREREG_C_COMMITTED"; exit 0; }
  if [ ! -f $P/status/R10_GO.txt ]; then
    big=$(find $P ~/.cache -type f -size +50M -newer $P/status/PREREG_C_AUDIT_R1.md 2>/dev/null | grep -v -E "r7_gears_venv|pip|huggingface/hub" | head -3)
    [ -n "$big" ] && { echo "$(date +%H:%M) DOWNLOAD_BEFORE_GO_ALERT: $big"; exit 0; }
  fi
  git -C $R ls-files | grep -q -E "round(9/p2|10)/shared/(lines|raw)/" && { echo "$(date +%H:%M) DATA_TRACKED_ALERT"; exit 0; }
  [ $(( $(date +%s) - start )) -gt 10800 ] && { echo "$(date +%H:%M) TIMEOUT 3h"; exit 0; }
  sleep 60
done
