#!/bin/zsh
setopt nullglob
# R11 phase-a watcher. Wakes controller when the PREREG_D draft (round11/status/COMMIT_REQUEST_A-T3.txt), the R11 data screen
# (R11_SCREEN_DONE) and the Q4/Q5 prior-art update (R11_LIT_DONE) are all in; on any round-11 download before R11_GO.txt;
# on data tracked in git or frozen R10 code changed; or after 6 h. Nudges idle agents (max 4 each).
C=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller
B=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
D=$C/directives/R11_第十一轮行动指令.md
start=$(date +%s); typeset -A nud idle
has(){ grep -q "^$2\$" "$1" 2>/dev/null; }
while true; do
  source $C/agy_env.sh
  p=0; [ -s $B/round11/status/COMMIT_REQUEST_A-T3.txt ] && p=1
  s1=0; has $B/round10/shared/R11_DATA_SCREEN.md R11_SCREEN_DONE && s1=1
  l5=0; has $B/round9/p2/vcc_linespecific/audit/LIT/PRIOR_ART_D.md R11_LIT_DONE && l5=1; has $B/round9/p2/status/A-T5.txt R11_LIT_DONE && l5=1
  o4=0; has $B/round10/status/A-T4.txt R11_OUTLINE_DONE && o4=1
  c2=0; has $B/round11/status/A-T2.txt R11_CODE_DRAFT_READY && c2=1
  [ $p -eq 1 ] && [ $s1 -eq 1 ] && [ $l5 -eq 1 ] && { echo "$(date +%H:%M) R11_PREREG_INPUTS_READY (outline $o4, code $c2)"; exit 0; }
  if [ ! -f $B/round11/status/R11_GO.txt ]; then
    n=$(find $B/round11/shared/raw $B/round11/shared/lines -type f 2>/dev/null | head -1)
    [ -n "$n" ] && { echo "$(date +%H:%M) R11_DOWNLOAD_BEFORE_GO_ALERT: $n"; exit 0; }
  fi
  (cd $B && shasum -a 256 -c round10/status/FREEZE_MANIFEST_R10.txt --quiet >/dev/null 2>&1) || { echo "$(date +%H:%M) FROZEN_R10_CODE_CHANGED_ALERT"; exit 0; }
  git -C /Users/mac/AI-Biology-Lab ls-files | grep -q -E "shared/(lines|raw)/|\.npz$|\.h5$|\.h5ad$|\.bin$" && { echo "$(date +%H:%M) DATA_TRACKED_ALERT"; exit 0; }
  for s cid flag in A-T3 6b165fe5-445b-46f7-9202-9efcc6b18385 $p A-T1 54060608-790c-477b-baeb-364849e4144e $s1 A-T5 c4d0971d-8452-497f-b5d2-466ba007691c $l5 A-T4 9046ae1f-74af-4e33-90d8-21cc389a8108 $o4 A-T2 4a6da262-f01a-45d3-8e7e-9ae4511c43a4 $c2; do
    [ $flag -eq 1 ] && continue
    st=$(python3 -c "
import sqlite3,os
c=sqlite3.connect('file:'+os.path.expanduser('~/.gemini/antigravity/conversation_summaries.db')+'?mode=ro',uri=True)
print(c.execute(\"select status from conversation_summaries where conversation_id='$cid'\").fetchone()[0])" 2>/dev/null)
    if [ "$st" = "CASCADE_RUN_STATUS_IDLE" ]; then idle[$s]=$(( ${idle[$s]:-0} + 1 )); else idle[$s]=0; fi
    if [ "${idle[$s]}" -ge 3 ] && [ "${nud[$s]:-0}" -lt 4 ]; then
      nud[$s]=$(( ${nud[$s]:-0} + 1 ))
      $AGY send-message "$cid" "【Claude 主控 · 续跑提醒 #${nud[$s]}】你已空闲。请继续 $D 里分给 $s 的任务，完成后按要求写完成标记。若有进程在后台运行，只需确认，不要重复启动。" >/dev/null 2>&1 && echo "$(date +%H:%M) nudged $s #${nud[$s]}"
      idle[$s]=-3
    fi
  done
  [ $(( $(date +%s) - start )) -gt 21600 ] && { echo "$(date +%H:%M) TIMEOUT 6h"; exit 0; }
  sleep 120
done
