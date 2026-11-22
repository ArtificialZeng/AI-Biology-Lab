#!/bin/zsh
# Round-10 phase-b watcher. Wakes controller on: R9 preprint draft ready, PREREG_C commit request, any large download before R10_GO, data tracked, or 8 h.
C=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller
P9=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10
R=/Users/mac/AI-Biology-Lab
start=$(date +%s); typeset -A nud idle; draft_seen=0
last(){ tail -1 "$1" 2>/dev/null | tr -d '[:space:]'; }
while true; do
  source $C/agy_env.sh
  grep -q R1_REVISED $P9/status/PAPER_DRAFT_READY.txt 2>/dev/null && { echo "$(date +%H:%M) PAPER_R1_REVISED"; exit 0; }
  [ -f $P/status/COMMIT_REQUEST_A-T3.txt ] && [ "$(last $P9/status/A-T1.txt)" = "R10_PLAN_DONE" ] && { echo "$(date +%H:%M) PREREG_C_READY"; exit 0; }
  if [ ! -f $P/status/R10_GO.txt ]; then
    big=$(find $P /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/shared -type f -size +50M -newer $C/directives/R10_第十轮行动指令.md 2>/dev/null | head -3)
    [ -n "$big" ] && { echo "$(date +%H:%M) DOWNLOAD_BEFORE_GO_ALERT: $big"; exit 0; }
  fi
  git -C $R ls-files | grep -q -E "round(9/p2|10)/shared/(lines|raw)/" && { echo "$(date +%H:%M) DATA_TRACKED_ALERT"; exit 0; }
  for s cid in A-T2 4a6da262-f01a-45d3-8e7e-9ae4511c43a4 A-T3 6b165fe5-445b-46f7-9202-9efcc6b18385 A-T1 54060608-790c-477b-baeb-364849e4144e A-T5 c4d0971d-8452-497f-b5d2-466ba007691c; do
    case $s in
      A-T2) [ "$(last $P/status/A-T2.txt)" = "R10_CODE_READY" ] && continue; msg="请继续：先完成 DEVIATIONS.md 和 S2b（秩 1 投影）补跑，再按 $C/directives/R10_第十轮行动指令.md 在 round10/src/ 实现 C 组管线，完成后在 round10/status/A-T2.txt 末行写 R10_CODE_READY。";;
      A-T3) msg="请继续：按 round9/p2/status/PAPER_REVIEW_R1.md 修改初稿（完成后 PAPER_DRAFT_READY.txt 写 R1_REVISED），按 round10/status/PREREG_C_AUDIT_R1.md 修改 PREREG_C 并写 COMMIT_REQUEST_A-T3.txt。";;
      A-T1) [ "$(last $P9/status/A-T1.txt)" = "R10_PLAN_DONE" ] && continue; msg="请按 $C/directives/R10_第十轮行动指令.md 写 round10/shared/R10_DOWNLOAD_PLAN.md（R10_GO.txt 出现前不下载），末行写 R10_PLAN_DONE。";;
      A-T5) [ "$(last $P9/status/A-T5.txt)" = "R10_LIT2_DONE" ] && continue; msg="请把 arXiv:2608.00152 补进 PRIOR_ART_C.md，末行写 R10_LIT2_DONE。";;
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
