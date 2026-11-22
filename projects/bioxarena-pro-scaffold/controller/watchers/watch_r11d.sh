#!/bin/zsh
setopt nullglob
# R11 phase-d watcher (download/QC + code readiness). Wakes controller when CONFIRMATORY_LINES_R11.txt exists, no round-11
# preparation process runs, and both code bases report ready; on scores before CODE_FREEZE_R11; raw downloads > 80 GB; data in git;
# frozen R10 code changed; or 8 h. Nudges idle agents with nothing running (max 5 each).
C=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller
B=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
P=$B/round11
start=$(date +%s); typeset -A nud idle
has(){ grep -q "^$2\$" "$1" 2>/dev/null; }
while true; do
  source $C/agy_env.sh
  prep=0; pgrep -f "round11/shared/scripts|prepare_r11|download_r11" >/dev/null && prep=1
  qc=0; [ -s $P/status/CONFIRMATORY_LINES_R11.txt ] && qc=1
  c2=0; has $P/status/A-T2.txt R11_CODE_READY && c2=1
  c3=0; has $P/status/A-T3.txt R11_AUDIT_CODE_READY && c3=1
  [ $qc -eq 1 ] && [ $prep -eq 0 ] && [ $c2 -eq 1 ] && [ $c3 -eq 1 ] && { echo "$(date +%H:%M) READY_FOR_CODE_FREEZE_R11"; exit 0; }
  if [ ! -f $P/status/CODE_FREEZE_R11.txt ]; then
    early=$(find $P/results $P/audit -type f -newer $P/status/R11_GO.txt 2>/dev/null | grep -v -E "__pycache__|\.py$|README|dev|smoke|test_|\.md$" | head -3)
    [ -n "$early" ] && { echo "$(date +%H:%M) RESULTS_BEFORE_FREEZE_ALERT: $early"; exit 0; }
  fi
  sz=$(du -sk $P/shared/raw 2>/dev/null | awk '{print int($1/1000000)}')
  [ "${sz:-0}" -gt 80 ] && { echo "$(date +%H:%M) DOWNLOAD_OVER_BUDGET_ALERT: raw ${sz} GB"; exit 0; }
  git -C /Users/mac/AI-Biology-Lab ls-files | grep -q -E "shared/(lines|raw)/|\.npz$|\.h5$|\.h5ad$|\.hdf5$|\.bin$" && { echo "$(date +%H:%M) DATA_TRACKED_ALERT"; exit 0; }
  (cd $B && shasum -a 256 -c round10/status/FREEZE_MANIFEST_R10.txt --quiet >/dev/null 2>&1) || { echo "$(date +%H:%M) FROZEN_R10_CODE_CHANGED_ALERT"; exit 0; }
  for s cid flag msg in \
    A-T1 54060608-790c-477b-baeb-364849e4144e $qc "请按 round11/status/R11_GO.txt 继续下载、可行性检查、单 sgRNA 过滤和只写计数的 QC，完成后写 R11_QC_COUNTS.md 和 CONFIRMATORY_LINES_R11.txt，并在 round9/p2/status/A-T1.txt 末行写 R11_QC_DONE。" \
    A-T2 4a6da262-f01a-45d3-8e7e-9ae4511c43a4 $c2 "请继续把 round11/src/ 对齐 PREREG_D v2 和补充 A（只用第九、十轮的线和合成数据测试），完成后在 round11/status/A-T2.txt 末行写 R11_CODE_READY。" \
    A-T3 6b165fe5-445b-46f7-9202-9efcc6b18385 $c3 "请继续在 round11/audit/ 独立实现第十一轮复算（只用第九、十轮的线和合成数据测试），完成后在 round11/status/A-T3.txt 末行写 R11_AUDIT_CODE_READY。"; do
    [ $flag -eq 1 ] && continue
    [ $s = A-T1 ] && [ $prep -eq 1 ] && { idle[$s]=0; continue; }
    st=$(python3 -c "
import sqlite3,os
c=sqlite3.connect('file:'+os.path.expanduser('~/.gemini/antigravity/conversation_summaries.db')+'?mode=ro',uri=True)
print(c.execute(\"select status from conversation_summaries where conversation_id='$cid'\").fetchone()[0])" 2>/dev/null)
    if [ "$st" = "CASCADE_RUN_STATUS_IDLE" ]; then idle[$s]=$(( ${idle[$s]:-0} + 1 )); else idle[$s]=0; fi
    if [ "${idle[$s]}" -ge 3 ] && [ "${nud[$s]:-0}" -lt 5 ]; then
      nud[$s]=$(( ${nud[$s]:-0} + 1 ))
      $AGY send-message "$cid" "【Claude 主控 · 续跑提醒 #${nud[$s]}】你已空闲且没有相关进程在跑。$msg 若有进程在后台运行，只需确认，不要重复启动。" >/dev/null 2>&1 && echo "$(date +%H:%M) nudged $s #${nud[$s]}"
      idle[$s]=-3
    fi
  done
  [ $(( $(date +%s) - start )) -gt 28800 ] && { echo "$(date +%H:%M) TIMEOUT 8h"; exit 0; }
  sleep 120
done
