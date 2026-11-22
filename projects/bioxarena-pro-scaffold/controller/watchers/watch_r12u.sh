#!/bin/bash
# bash 3.2. Ruling-2 round: only NEW lines in status files count; markers must stand alone on a line.
# Alerts: prereg change, A549 line files changed, results before freeze, a downstream agent idle > 20 min after its trigger.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
sig() { find round12/shared/lines/a549_ko round12/shared/lines/a549_ko_cas9pos -type f -exec stat -f "%N %z %m" {} \; 2>/dev/null | sort | shasum -a 256 | cut -c1-64; }
# counts standalone marker lines; a marker is new when its count rises above the start count (robust to file rewrites)
cnt() { grep -cE "^[[:space:]\`*]*${2}[[:space:]\`*]*$" round12/status/$1 2>/dev/null; }
A0=$(sig); n1=$(cnt A-T1.txt R12_QC_DONE2); n2=$(cnt A-T2.txt R12_CODE_PATCHED2); n3=$(cnt A-T3.txt R12_AUDIT_PATCHED2)
nok=$(cnt A-T3.txt R12_READER_XCHECK_OK); nmm=$(cnt A-T3.txt R12_READER_XCHECK_MISMATCH)
has() { [ "$(cnt $1 $3)" -gt "$2" ]; }
echo "$(date +%H:%M) start"; seen=""; t_p2=""; t_q2=""
for i in $(seq 1 360); do
  now=$(date +%s)
  [ "$(shasum -a 256 round12/prereg/PREREG_E.md | cut -c1-64)" = d10f92ced16a4519c43ea2b92b8b7ee881ab3d1764e39ffd10580473e4108924 ] || { echo "$(date +%H:%M) PREREG_E_CHANGED_ALERT"; exit 0; }
  [ "$(sig)" = "$A0" ] || { echo "$(date +%H:%M) A549_LINES_CHANGED_ALERT"; exit 0; }
  if [ ! -s round12/status/CODE_FREEZE_R12.txt ] && [ -n "$(find round12/results -type f -not -path "*/dev_*" 2>/dev/null | head -1)" ]; then echo "$(date +%H:%M) RESULTS_BEFORE_FREEZE_ALERT"; exit 0; fi
  for spec in "A-T2.txt:$n2:R12_CODE_PATCHED2" "A-T3.txt:$n3:R12_AUDIT_PATCHED2" "A-T1.txt:$n1:R12_QC_DONE2" "A-T3.txt:$nok:R12_READER_XCHECK_OK" "A-T3.txt:$nmm:R12_READER_XCHECK_MISMATCH"; do
    f=${spec%%:*}; rest=${spec#*:}; n=${rest%%:*}; m=${rest#*:}
    case " $seen " in *" $m "*) ;; *) if has "$f" "$n" "$m"; then seen="$seen $m"; echo "$(date +%H:%M) $m";
      [ "$m" = R12_CODE_PATCHED2 ] && t_p2=$now; [ "$m" = R12_QC_DONE2 ] && t_q2=$now; fi;; esac
  done
  case " $seen " in *" R12_READER_XCHECK_OK "*|*" R12_READER_XCHECK_MISMATCH "*) echo "$(date +%H:%M) DONE:$seen"; exit 0;; esac
  if [ -n "$t_p2" ] && [ -z "$t_q2" ] && [ $(( now - t_p2 )) -gt 1200 ] && [ $(( now - $(stat -f %m round12/status/A-T1.txt) )) -gt 1200 ]; then echo "$(date +%H:%M) A-T1_STALL_ALERT:$seen"; exit 0; fi
  if [ -n "$t_q2" ] && [ $(( now - t_q2 )) -gt 1200 ] && [ $(( now - $(stat -f %m round12/status/A-T3.txt) )) -gt 1200 ]; then echo "$(date +%H:%M) A-T3_STALL_ALERT:$seen"; exit 0; fi
  sleep 30
done
echo "$(date +%H:%M) TIMEOUT:$seen"
