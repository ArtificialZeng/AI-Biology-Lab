#!/bin/zsh
setopt nullglob
# Round-10 phase-e watcher. Code from A-T2/A-T3 is already reviewed. Wakes controller when CONFIRMATORY_LINES.txt exists
# and the line-preparation process has finished (time to freeze); on scores before CODE_FREEZE_R10; raw downloads > 60 GB;
# data tracked in git; preparation stopping without QC output; or 6 h. Nudges A-T1 only when idle AND no prep process runs.
C=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10
R=/Users/mac/AI-Biology-Lab
start=$(date +%s); nud=0; idle=0; prep_seen=0
while true; do
  source $C/agy_env.sh
  prep=0; pgrep -f prepare_r10_lines.py >/dev/null && prep=1 && prep_seen=1
  qc=0; [ -f $P/status/CONFIRMATORY_LINES.txt ] && qc=1
  [ $qc -eq 1 ] && [ $prep -eq 0 ] && { echo "$(date +%H:%M) READY_FOR_CODE_FREEZE"; exit 0; }
  if [ ! -f $P/status/CODE_FREEZE_R10.txt ]; then
    early=$(find $P/results $P/audit -type f -newer $P/status/R10_GO.txt 2>/dev/null | grep -v -E "test_|\.py$|__pycache__|README|/dev/|_dev_|r9dev|r9_dev|\.md$" | head -3)
    [ -n "$early" ] && { echo "$(date +%H:%M) RESULTS_BEFORE_FREEZE_ALERT: $early"; exit 0; }
  fi
  sz=$(du -sk $P/shared/raw 2>/dev/null | awk '{print int($1/1048576)}')
  [ "${sz:-0}" -gt 60 ] && { echo "$(date +%H:%M) DOWNLOAD_OVER_BUDGET_ALERT: raw ${sz} GB"; exit 0; }
  git -C $R ls-files | grep -q -E "round(9/p2|10)/shared/(lines|raw)/|\.h5$|\.h5ad$|\.npz$|\.bin$" && { echo "$(date +%H:%M) DATA_TRACKED_ALERT"; exit 0; }
  if [ $prep -eq 0 ] && [ $qc -eq 0 ]; then
    st=$(python3 -c "
import sqlite3,os
c=sqlite3.connect('file:'+os.path.expanduser('~/.gemini/antigravity/conversation_summaries.db')+'?mode=ro',uri=True)
print(c.execute(\"select status from conversation_summaries where conversation_id='54060608-790c-477b-baeb-364849e4144e'\").fetchone()[0])" 2>/dev/null)
    if [ "$st" = "CASCADE_RUN_STATUS_IDLE" ]; then idle=$((idle+1)); else idle=0; fi
    if [ $idle -ge 2 ] && [ $nud -lt 4 ]; then
      nud=$((nud+1))
      $AGY send-message 54060608-790c-477b-baeb-364849e4144e "【Claude 主控 · 提醒 #$nud】cd4t_rest 的预处理进程已结束（或未在运行）。请按主控 11:50 的指令：检查 cd4t_rest 产物完整，删除 _tmp_shard_* 临时文件，写 round10/shared/R10_QC_COUNTS.md 和 round10/status/CONFIRMATORY_LINES.txt（只含计数，不算任何指标），在 round9/p2/status/A-T1.txt 末行写 R10_QC_DONE。如果进程是异常退出，先写明原因再重跑。" >/dev/null 2>&1 && echo "$(date +%H:%M) nudged A-T1 #$nud"
      idle=-3
    fi
  fi
  [ $(( $(date +%s) - start )) -gt 21600 ] && { echo "$(date +%H:%M) TIMEOUT 6h"; exit 0; }
  sleep 120
done
