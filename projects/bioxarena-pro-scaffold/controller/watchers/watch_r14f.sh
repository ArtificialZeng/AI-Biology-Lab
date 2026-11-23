#!/bin/bash
# bash 3.2. R14 feasibility. Exits on R14_FEASIBILITY_DONE, a big file under round14 (no matrix download allowed yet), or 60-min stall.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
T0=$(date +%s); echo "$(date +%H:%M) start"
for i in $(seq 1 240); do
  now=$(date +%s)
  big=$(find round14 -type f -size +100M 2>/dev/null | head -1); [ -n "$big" ] && { echo "$(date +%H:%M) BIG_FILE_BEFORE_LOCK_ALERT $big"; exit 0; }
  grep -qE "^[[:space:]\`*]*R14_FEASIBILITY_DONE[[:space:]\`*]*$" round14/status/A-T1.txt 2>/dev/null && { echo "$(date +%H:%M) R14_FEASIBILITY_DONE"; exit 0; }
  last=$( { find round14 -type f -exec stat -f %m {} \; ; } 2>/dev/null | sort -n | tail -1); [ -z "$last" ] && last=$T0; [ "$last" -lt "$T0" ] && last=$T0
  [ $(( now - last )) -gt 3600 ] && { echo "$(date +%H:%M) R14_STALL_ALERT"; exit 0; }
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT"
