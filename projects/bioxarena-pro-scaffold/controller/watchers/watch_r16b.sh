#!/bin/bash
# R16 watcher, phase b (after approved weight download). Fixed baselines; exits on the first event.
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cd "$P" || exit 9
BASE_HEAD=80f7c17
M=round16/shared/models
REF=$(mktemp); touch -t 202611241105 "$REF"
NFILES=$(find $M -type f | wc -l | tr -d ' ')
MSIG=$(find $M -type f -exec stat -f '%N %z %m' {} \; | sort | md5)
last_change=$(date +%s); prev_sig=""
while true; do
  for m in "round16/status/A-T5.txt:PRIOR_ART_I_DONE" "round16/status/A-T3.txt:PREREG_I_AT3_AUDIT_DONE" "round16/status/MAIN.txt:R16_GATE2_SCHEMA_READY"; do
    f=${m%%:*}; k=${m##*:}
    [ -f "$f" ] && grep -qx "$k" "$f" && { echo "EVENT marker $k in $f at $(date '+%F %T')"; exit 0; }
  done
  w=$(find /Users/mac/AI-Biology-Lab ~/.cache/huggingface/hub -newer "$REF" \( -name '*.ckpt' -o -name '*.safetensors' -o -name 'adata_pred*' -o -name 'adata_real*' -o -name 'eval_final*' \) -not -path "$P/$M/*" -not -path "$P/round16/.venv/*" 2>/dev/null | head -3)
  [ -n "$w" ] && { echo "GUARD model/eval file outside approved dir: $w"; exit 2; }
  [ -n "$(ls -d ~/.cache/huggingface/hub/models--arcinstitute* 2>/dev/null)" ] && { echo "GUARD arcinstitute HF cache created"; exit 2; }
  [ "$(find $M -type f | wc -l | tr -d ' ')" != "$NFILES" ] && { echo "GUARD file count changed in $M"; exit 2; }
  [ "$(find $M -type f -exec stat -f '%N %z %m' {} \; | sort | md5)" != "$MSIG" ] && { echo "GUARD approved model files modified"; exit 2; }
  [ -n "$(ls -d /opt/anaconda3/lib/python3*/site-packages/arc_state* 2>/dev/null)" ] && { echo "GUARD arc-state installed in base env"; exit 2; }
  c=$(find /Users/mac/AI-Biology-Lab -maxdepth 4 -type d -iname 'CellularAI*' -newer "$REF" 2>/dev/null | head -1)
  [ -n "$c" ] && { echo "GUARD CellularAI clone: $c"; exit 2; }
  h=$(git -C /Users/mac/AI-Biology-Lab rev-parse --short HEAD 2>/dev/null)
  [ "$h" != "$BASE_HEAD" ] && { echo "GUARD new git commit HEAD=$h"; exit 2; }
  git -C /Users/mac/AI-Biology-Lab ls-files | grep -q -E '^projects/bioxarena-pro-scaffold/(paper|round16/shared/models)/' && { echo "GUARD paper/ or models tracked in git"; exit 2; }
  big=$(find "$P" -newer "$REF" -type f -size +100M -not -path '*/.git/*' -not -path "$P/$M/*" -not -path "$P/round16/.venv/*" 2>/dev/null | head -3)
  [ -n "$big" ] && { echo "GUARD file >100MB: $big"; exit 2; }
  sig=$(ls -l round16/status round16/prereg round16/src round16/audit 2>/dev/null | md5)
  if [ "$sig" != "$prev_sig" ]; then prev_sig=$sig; last_change=$(date +%s); fi
  [ $(( $(date +%s) - last_change )) -gt 5400 ] && { echo "STALL no round16 changes for 90 min at $(date '+%F %T')"; exit 3; }
  sleep 60
done
