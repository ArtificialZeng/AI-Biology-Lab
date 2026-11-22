#!/bin/zsh
setopt nullglob
# Round-10 evaluation watcher. Wakes controller when A-T2 writes R10_EVAL_DONE and A-T3 writes R10_AUDIT_DONE; on a deviation
# report; if any frozen .py file changes (manifest check); or after 6 h. Nudges idle A-T2/A-T3 (max 4 each).
C=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller
B=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
P=$B/round10
start=$(date +%s); typeset -A nud idle
last(){ tail -1 "$1" 2>/dev/null | tr -d '[:space:]'; }
while true; do
  source $C/agy_env.sh
  e2=0; [ "$(last $P/status/A-T2.txt)" = "R10_EVAL_DONE" ] && e2=1
  e3=0; [ "$(last $P/status/A-T3.txt)" = "R10_AUDIT_DONE" ] && e3=1
  [ $e2 -eq 1 ] && [ $e3 -eq 1 ] && { echo "$(date +%H:%M) EVAL_AND_AUDIT_DONE"; exit 0; }
  [ -f $P/DEVIATIONS_R10.md ] && { echo "$(date +%H:%M) DEVIATION_REPORTED"; exit 0; }
  (cd $B && shasum -a 256 -c $P/status/FREEZE_MANIFEST_R10.txt --quiet >/dev/null 2>&1) || { echo "$(date +%H:%M) FROZEN_CODE_CHANGED_ALERT"; exit 0; }
  for s cid in A-T2 4a6da262-f01a-45d3-8e7e-9ae4511c43a4 A-T3 6b165fe5-445b-46f7-9202-9efcc6b18385; do
    case $s in
      A-T2) [ $e2 -eq 1 ] && continue; msg="请继续第十轮主评测（见主控 12:22 指令），完成后在 round10/status/A-T2.txt 末行写 R10_EVAL_DONE。若进程仍在后台运行，只需确认并等待，不要重复启动。";;
      A-T3) [ $e3 -eq 1 ] && continue; msg="请继续第十轮独立复算和 C3 交叉核对（见主控 12:22 指令），完成后在 round10/status/A-T3.txt 末行写 R10_AUDIT_DONE。若进程仍在后台运行，只需确认并等待，不要重复启动。";;
    esac
    pgrep -f "run_eval_r10.py --phase holdout_r10" >/dev/null && [ $s = A-T2 ] && { idle[$s]=0; continue; }
    pgrep -f "run_r10_audit.py --mode holdout_r10" >/dev/null && [ $s = A-T3 ] && { idle[$s]=0; continue; }
    st=$(python3 -c "
import sqlite3,os
c=sqlite3.connect('file:'+os.path.expanduser('~/.gemini/antigravity/conversation_summaries.db')+'?mode=ro',uri=True)
print(c.execute(\"select status from conversation_summaries where conversation_id='$cid'\").fetchone()[0])" 2>/dev/null)
    if [ "$st" = "CASCADE_RUN_STATUS_IDLE" ]; then idle[$s]=$(( ${idle[$s]:-0} + 1 )); else idle[$s]=0; fi
    if [ "${idle[$s]}" -ge 3 ] && [ "${nud[$s]:-0}" -lt 4 ]; then
      nud[$s]=$(( ${nud[$s]:-0} + 1 ))
      $AGY send-message "$cid" "【Claude 主控 · 续跑提醒 #${nud[$s]}】你已空闲且没有评测进程在跑。$msg" >/dev/null 2>&1 && echo "$(date +%H:%M) nudged $s #${nud[$s]}"
      idle[$s]=-3
    fi
  done
  [ $(( $(date +%s) - start )) -gt 21600 ] && { echo "$(date +%H:%M) TIMEOUT 6h"; exit 0; }
  sleep 120
done
