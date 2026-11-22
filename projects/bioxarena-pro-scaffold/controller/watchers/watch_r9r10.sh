#!/bin/zsh
# Watcher (round-9 wrap-up + round-10 screening). Auto-detects the live Antigravity daemon port.
# Wakes the controller on: preprint draft ready, R10 screening done (A-T5 and A-T1), data tracked in git, or 8 h.
C=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2
R=/Users/mac/AI-Biology-Lab
start=$(date +%s); typeset -A nud idle
last(){ tail -1 "$1" 2>/dev/null | tr -d '[:space:]'; }
while true; do
  source $C/agy_env.sh
  [ -f $P/status/PAPER_DRAFT_READY.txt ] && { echo "$(date +%H:%M) PAPER_DRAFT_READY"; exit 0; }
  [ "$(last $P/status/A-T5.txt)" = "R10_SCREEN_DONE" ] && [ "$(last $P/status/A-T1.txt)" = "R10_SCREEN_DONE" ] && { echo "$(date +%H:%M) R10_SCREEN_READY"; exit 0; }
  git -C $R ls-files | grep -q "round9/p2/shared/lines/" && { echo "$(date +%H:%M) DATA_TRACKED_ALERT"; exit 0; }
  for s cid in A-T2 4a6da262-f01a-45d3-8e7e-9ae4511c43a4 A-T3 6b165fe5-445b-46f7-9202-9efcc6b18385 A-T5 c4d0971d-8452-497f-b5d2-466ba007691c A-T1 54060608-790c-477b-baeb-364849e4144e; do
    case $s in
      A-T2) [ "$(last $P/status/A-T2.txt)" = "TASK_DONE" ] && continue; msg="请继续重跑 S1、S3、S4a、S4b（C2 固定随机数），结果写进 results/*_sens_*，完成后在 status/A-T2.txt 末行写 TASK_DONE。";;
      A-T3) msg="请继续：重跑 S1、S3、S4 写进 audit/sensitivity/，写最终 AUDIT_REPORT.md，再起草 $P/paper/PREPRINT_DRAFT.md，写完创建 $P/status/PAPER_DRAFT_READY.txt。";;
      A-T5) [ "$(last $P/status/A-T5.txt)" = "R10_SCREEN_DONE" ] && continue; msg="请按 $C/directives/R10_第十轮行动指令.md 完成 PRIOR_ART_C.md，末行写 R10_SCREEN_DONE。";;
      A-T1) [ "$(last $P/status/A-T1.txt)" = "R10_SCREEN_DONE" ] && continue; msg="请按 $C/directives/R10_第十轮行动指令.md 完成 shared/R10_DATA_SCREEN.md（不下载表达矩阵），末行写 R10_SCREEN_DONE。";;
    esac
    st=$(python3 -c "
import sqlite3,os
c=sqlite3.connect('file:'+os.path.expanduser('~/.gemini/antigravity/conversation_summaries.db')+'?mode=ro',uri=True)
print(c.execute(\"select status from conversation_summaries where conversation_id='$cid'\").fetchone()[0])" 2>/dev/null)
    if [ "$st" = "CASCADE_RUN_STATUS_IDLE" ]; then idle[$s]=$(( ${idle[$s]:-0} + 1 )); else idle[$s]=0; fi
    if [ "${idle[$s]}" -ge 4 ] && [ "${nud[$s]:-0}" -lt 4 ]; then
      nud[$s]=$(( ${nud[$s]:-0} + 1 ))
      $AGY send-message "$cid" "【Claude 主控 · 续跑提醒 #${nud[$s]}】你已空闲。$msg" >/dev/null 2>&1 && echo "$(date +%H:%M) nudged $s #${nud[$s]}"
      idle[$s]=-4
    fi
  done
  [ $(( $(date +%s) - start )) -gt 28800 ] && { echo "$(date +%H:%M) TIMEOUT 8h"; exit 0; }
  sleep 120
done
