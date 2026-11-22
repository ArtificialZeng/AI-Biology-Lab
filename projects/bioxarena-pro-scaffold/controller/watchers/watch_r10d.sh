#!/bin/zsh
setopt nullglob
# Round-10 phase-d watcher (download/QC + code readiness). Wakes controller when QC counts and both code bases are ready (time to freeze),
# on early evaluation (results before CODE_FREEZE_R10), download over 60 GB, data tracked, paper revised, or 8 h. Nudges idle agents.
C=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller
P9=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10
R=/Users/mac/AI-Biology-Lab
start=$(date +%s); typeset -A nud idle; paper_said=0
last(){ tail -1 "$1" 2>/dev/null | tr -d '[:space:]'; }
while true; do
  source $C/agy_env.sh
  qc=0; [ -f $P/status/CONFIRMATORY_LINES.txt ] && qc=1
  c2=0; [ "$(last $P/status/A-T2.txt)" = "R10_CODE_READY" ] && c2=1
  c3=0; grep -q R10_AUDIT_CODE_READY $P/status/A-T3.txt 2>/dev/null && c3=1
  [ $qc -eq 1 ] && [ $c2 -eq 1 ] && [ $c3 -eq 1 ] && { echo "$(date +%H:%M) READY_FOR_CODE_FREEZE"; exit 0; }
  if [ ! -f $P/status/CODE_FREEZE_R10.txt ]; then
    early=$(find $P/results $P/audit -type f -newer $P/status/R10_GO.txt 2>/dev/null | grep -v -E "test_|\.py$|__pycache__|README|/dev/|_dev_|r9dev|\.md$" | head -3)
    [ -n "$early" ] && { echo "$(date +%H:%M) RESULTS_BEFORE_FREEZE_ALERT: $early"; exit 0; }
  fi
  sz=$(du -sk $P/shared ~/.cache/r10* 2>/dev/null | awk '{s+=$1} END{print int(s/1048576)}')
  [ "${sz:-0}" -gt 62 ] && { echo "$(date +%H:%M) DOWNLOAD_OVER_BUDGET_ALERT: ${sz} GB"; exit 0; }
  git -C $R ls-files | grep -q -E "round(9/p2|10)/shared/(lines|raw)/|\.h5$|\.h5ad$" && { echo "$(date +%H:%M) DATA_TRACKED_ALERT"; exit 0; }

  for s cid in A-T1 54060608-790c-477b-baeb-364849e4144e A-T2 4a6da262-f01a-45d3-8e7e-9ae4511c43a4 A-T3 6b165fe5-445b-46f7-9202-9efcc6b18385; do
    case $s in
      A-T1) [ $qc -eq 1 ] && continue; msg="请按 $P/status/R10_GO.txt 继续下载、单 sgRNA 过滤和 QC（只报告计数），完成后写 R10_QC_COUNTS.md 和 CONFIRMATORY_LINES.txt，并在 round9/p2/status/A-T1.txt 末行写 R10_QC_DONE。";;
      A-T2) [ $c2 -eq 1 ] && continue; msg="请继续：第九轮 DEVIATIONS.md 和 S2b（秩 1 投影）补跑，以及 round10/src/ C 组管线（只用第九轮 6 条线测试），完成后在 round10/status/A-T2.txt 末行写 R10_CODE_READY。";;
      A-T3) [ $c3 -eq 1 ] && continue; msg="请继续：在 round10/audit/ 独立实现 C 组复算，完成后在 round10/status/A-T3.txt 末行写 R10_AUDIT_CODE_READY。";;
    esac
    st=$(python3 -c "
import sqlite3,os
c=sqlite3.connect('file:'+os.path.expanduser('~/.gemini/antigravity/conversation_summaries.db')+'?mode=ro',uri=True)
print(c.execute(\"select status from conversation_summaries where conversation_id='$cid'\").fetchone()[0])" 2>/dev/null)
    if [ "$st" = "CASCADE_RUN_STATUS_IDLE" ]; then idle[$s]=$(( ${idle[$s]:-0} + 1 )); else idle[$s]=0; fi
    if [ "${idle[$s]}" -ge 4 ] && [ "${nud[$s]:-0}" -lt 5 ]; then
      nud[$s]=$(( ${nud[$s]:-0} + 1 ))
      $AGY send-message "$cid" "【Claude 主控 · 续跑提醒 #${nud[$s]}】你已空闲。$msg" >/dev/null 2>&1 && echo "$(date +%H:%M) nudged $s #${nud[$s]}"
      idle[$s]=-4
    fi
  done
  [ $(( $(date +%s) - start )) -gt 28800 ] && { echo "$(date +%H:%M) TIMEOUT 8h"; exit 0; }
  sleep 120
done
