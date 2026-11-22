#!/bin/zsh
setopt nullglob
# R10-E watcher. Wakes controller when A-T2 (R10E_ABLATION_DONE), A-T3 (R10E_AUDIT_DONE) and A-T5 (R10E_LIT_DONE) are all done;
# on frozen-code change or data tracked in git; or after 6 h. Nudges idle agents with nothing running (max 4 each).
C=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller
B=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
start=$(date +%s); typeset -A nud idle
D=$C/directives/R10E_探索性消融与合稿指令.md
done_(){ grep -q "$2" "$1" 2>/dev/null && [ "$(tail -1 $1 | tr -d '[:space:]')" = "$2" ] || grep -q "^$2\$" "$1" 2>/dev/null; }
while true; do
  source $C/agy_env.sh
  a2=0; done_ $B/round10/status/A-T2.txt R10E_ABLATION_DONE && a2=1
  a3=0; done_ $B/round10/status/A-T3.txt R10E_AUDIT_DONE && a3=1
  a5=0; done_ $B/round9/p2/status/A-T5.txt R10E_LIT_DONE && a5=1
  a1=0; done_ $B/round10/shared/R11_DATA_SCREEN.md R11_SCREEN_DONE && a1=1
  a4=0; done_ $B/round10/status/A-T4.txt R10E_OUTLINE_DONE && a4=1
  [ $a2 -eq 1 ] && [ $a3 -eq 1 ] && [ $a5 -eq 1 ] && { echo "$(date +%H:%M) R10E_CORE_DONE (A-T1 $a1, A-T4 $a4)"; exit 0; }
  (cd $B && shasum -a 256 -c round10/status/FREEZE_MANIFEST_R10.txt --quiet >/dev/null 2>&1) || { echo "$(date +%H:%M) FROZEN_CODE_CHANGED_ALERT"; exit 0; }
  git -C /Users/mac/AI-Biology-Lab ls-files | grep -q -E "shared/(lines|raw)/|\.npz$|\.h5$|\.h5ad$|\.bin$" && { echo "$(date +%H:%M) DATA_TRACKED_ALERT"; exit 0; }
  for s cid flag in A-T2 4a6da262-f01a-45d3-8e7e-9ae4511c43a4 $a2 A-T3 6b165fe5-445b-46f7-9202-9efcc6b18385 $a3 A-T5 c4d0971d-8452-497f-b5d2-466ba007691c $a5 A-T1 54060608-790c-477b-baeb-364849e4144e $a1 A-T4 9046ae1f-74af-4e33-90d8-21cc389a8108 $a4; do
    [ $flag -eq 1 ] && continue
    pgrep -f "round10/(exploratory|audit/exploratory)" >/dev/null && [[ $s = A-T2 || $s = A-T3 ]] && { idle[$s]=0; continue; }
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
