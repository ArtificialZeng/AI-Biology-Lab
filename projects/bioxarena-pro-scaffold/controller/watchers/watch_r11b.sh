#!/bin/zsh
setopt nullglob
# R11 phase-b watcher. Wakes controller when A-T3 has revised PREREG_D (COMMIT_REQUEST sha no longer 87e428b4...) and A-T1 has
# written R11_PLAN_DONE; on any round-11 download before R11_GO.txt; on data tracked in git; or after 4 h. Nudges idle A-T3/A-T1.
C=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller
B=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
start=$(date +%s); typeset -A nud idle
while true; do
  source $C/agy_env.sh
  r=0; [ -s $B/round11/status/COMMIT_REQUEST_A-T3.txt ] && ! grep -q 87e428b48ee7d87d0ecd90887d0a0e77a664cda764bb5ac4d1fbdc418b0c9f1f $B/round11/status/COMMIT_REQUEST_A-T3.txt && r=1
  pl=0; grep -q "^R11_PLAN_DONE$" $B/round9/p2/status/A-T1.txt 2>/dev/null && pl=1
  [ $r -eq 1 ] && [ $pl -eq 1 ] && { echo "$(date +%H:%M) PREREG_D_REVISED_AND_PLAN_READY"; exit 0; }
  if [ ! -f $B/round11/status/R11_GO.txt ]; then
    n=$(find $B/round11/shared/raw $B/round11/shared/lines -type f 2>/dev/null | head -1)
    [ -n "$n" ] && { echo "$(date +%H:%M) R11_DOWNLOAD_BEFORE_GO_ALERT: $n"; exit 0; }
  fi
  git -C /Users/mac/AI-Biology-Lab ls-files | grep -q -E "shared/(lines|raw)/|\.npz$|\.h5$|\.h5ad$|\.hdf5$|\.bin$" && { echo "$(date +%H:%M) DATA_TRACKED_ALERT"; exit 0; }
  for s cid flag msg in \
    A-T3 6b165fe5-445b-46f7-9202-9efcc6b18385 $r "请按 round11/status/PREREG_D_AUDIT_R1.md 的 Q1–Q10 修改 PREREG_D.md，并更新 COMMIT_REQUEST_A-T3.txt 的 sha256。" \
    A-T1 54060608-790c-477b-baeb-364849e4144e $pl "请写 round11/shared/R11_DOWNLOAD_PLAN.md（只看元数据），完成后在 round9/p2/status/A-T1.txt 末行写 R11_PLAN_DONE。"; do
    [ $flag -eq 1 ] && continue
    st=$(python3 -c "
import sqlite3,os
c=sqlite3.connect('file:'+os.path.expanduser('~/.gemini/antigravity/conversation_summaries.db')+'?mode=ro',uri=True)
print(c.execute(\"select status from conversation_summaries where conversation_id='$cid'\").fetchone()[0])" 2>/dev/null)
    if [ "$st" = "CASCADE_RUN_STATUS_IDLE" ]; then idle[$s]=$(( ${idle[$s]:-0} + 1 )); else idle[$s]=0; fi
    if [ "${idle[$s]}" -ge 3 ] && [ "${nud[$s]:-0}" -lt 4 ]; then
      nud[$s]=$(( ${nud[$s]:-0} + 1 ))
      $AGY send-message "$cid" "【Claude 主控 · 续跑提醒 #${nud[$s]}】你已空闲。$msg" >/dev/null 2>&1 && echo "$(date +%H:%M) nudged $s #${nud[$s]}"
      idle[$s]=-3
    fi
  done
  [ $(( $(date +%s) - start )) -gt 14400 ] && { echo "$(date +%H:%M) TIMEOUT 4h"; exit 0; }
  sleep 120
done
