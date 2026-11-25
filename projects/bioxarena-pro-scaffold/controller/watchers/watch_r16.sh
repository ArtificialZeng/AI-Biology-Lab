#!/bin/bash
# R16 watcher (controller). Fixed baselines from the known state at 2026-11-24 11:05 +08:00.
# Exits on the first event: a whole-line marker, a guard violation, or a 90-minute stall.
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cd "$P" || exit 9
BASE_HEAD=80f7c17
REF=$(mktemp); touch -t 202611241105 "$REF"
last_change=$(date +%s)
prev_sig=""
while true; do
  for m in "round16/status/A-T5.txt:PRIOR_ART_I_DONE" "round16/status/A-T3.txt:PREREG_I_AT3_AUDIT_DONE"; do
    f=${m%%:*}; k=${m##*:}
    [ -f "$f" ] && grep -qx "$k" "$f" && { echo "EVENT marker $k in $f at $(date '+%F %T')"; exit 0; }
  done
  w=$(find /Users/mac/AI-Biology-Lab ~/.cache/huggingface/hub -newer "$REF" \( -name '*.ckpt' -o -name '*.safetensors' -o -name 'pert_onehot_map.pt' \) 2>/dev/null | head -3)
  [ -n "$w" ] && { echo "GUARD model weights appeared: $w"; exit 2; }
  d=$(ls -d ~/.cache/huggingface/hub/models--arcinstitute* 2>/dev/null | head -1)
  [ -n "$d" ] && { echo "GUARD arcinstitute HF cache created: $d"; exit 2; }
  s=$(ls -d /opt/anaconda3/lib/python3*/site-packages/arc_state* /opt/homebrew/lib/python3*/site-packages/arc_state* ~/.local/lib/python3*/site-packages/arc_state* 2>/dev/null | head -1)
  [ -n "$s" ] && { echo "GUARD arc-state installed: $s"; exit 2; }
  c=$(find /Users/mac/AI-Biology-Lab -maxdepth 4 -type d -iname 'CellularAI*' -newer "$REF" 2>/dev/null | head -1)
  [ -n "$c" ] && { echo "GUARD CellularAI clone: $c"; exit 2; }
  h=$(git -C /Users/mac/AI-Biology-Lab rev-parse --short HEAD 2>/dev/null)
  [ "$h" != "$BASE_HEAD" ] && { echo "GUARD new git commit HEAD=$h"; exit 2; }
  git -C /Users/mac/AI-Biology-Lab ls-files | grep -q '^projects/bioxarena-pro-scaffold/paper/' && { echo "GUARD paper/ tracked in git"; exit 2; }
  big=$(find "$P" -newer "$REF" -type f -size +100M -not -path '*/.git/*' 2>/dev/null | head -3)
  [ -n "$big" ] && { echo "GUARD file >100MB: $big"; exit 2; }
  sig=$(ls -l round16/status round16/prereg round16/src round16/audit 2>/dev/null | md5)
  if [ "$sig" != "$prev_sig" ]; then prev_sig=$sig; last_change=$(date +%s); fi
  [ $(( $(date +%s) - last_change )) -gt 5400 ] && { echo "STALL no round16 changes for 90 min at $(date '+%F %T')"; exit 3; }
  sleep 60
done
