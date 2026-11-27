#!/bin/bash
# R16 run-phase watcher (controller). Exits on the first event or guard. Baselines fixed at 2026-11-26 14:55 (after R16_GO 7106a3b).
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cd "$P" || exit 9
BASE_HEAD=${BASE_HEAD:-7106a3b}
M=round16/shared/models
S=round16/status
REF=$(mktemp); touch -t ${REF_TIME:-202611261455} "$REF"
MAN=${MAN:-FREEZE_MANIFEST_R16.txt}
HB=${HB:-/private/tmp/claude-503/r16e_heartbeat.txt}
mk() { local c; c=$(grep -c "^$2" "$S/$1" 2>/dev/null); echo "${c:-0}"; }
MARKS="${EXTRA_MARKS:-} MAIN.txt:R16_MAIN_REAL_DONE MAIN.txt:R16_MAIN_DEVIATION MAIN.txt:R16_ADAPTER_CHECK_FAIL MAIN.txt:R16_MAIN_SENSITIVITY_DONE A-T3.txt:R16_AT3_REAL_DONE A-T3.txt:R16_AT3_DEVIATION A-T3.txt:R16_AT3_SECONDARY_DONE"
BASEF=$(mktemp); for m in $MARKS; do echo "$m $(mk "${m%%:*}" "${m##*:}")" >> "$BASEF"; done
MSIG=$(find $M -type f -exec stat -f '%N %z %m' {} \; | sort | md5)
DSIG() { grep -v '^#' $S/FREEZE_MANIFEST_R16_DATA.txt | awk '{print $2}' | while read -r f; do [ -f "$f" ] || f="$P/$f"; stat -f '%N %z %m' "$f" 2>/dev/null || echo "MISSING $f"; done | md5; }
D0=$(DSIG)
seen_main=0; seen_at3=0; gone_main=0; gone_at3=0; i=0; T0=$(date +%s)
while true; do
  for m in $MARKS; do
    c=$(mk "${m%%:*}" "${m##*:}")
    [ "$c" -gt "$(awk -v k="$m" '$1==k{print $2}' "$BASEF")" ] && { echo "EVENT marker ${m##*:} in ${m%%:*} at $(date '+%F %T')"; exit 0; }
  done
  pm=$(pgrep -f "${MAIN_PAT:-round16/src/run_eval_r16.py}" | head -1); pa=$(pgrep -f "${AT3_PAT:-round16/audit/run_audit_r16.py}" | head -1)
  [ -n "$pm" ] && { seen_main=1; gone_main=0; }; [ -n "$pa" ] && { seen_at3=1; gone_at3=0; }
  [ "$seen_main" = 1 ] && [ -z "$pm" ] && gone_main=$((gone_main + 1))
  [ "$seen_at3" = 1 ] && [ -z "$pa" ] && gone_at3=$((gone_at3 + 1))
  [ "$gone_main" -ge 10 ] && { echo "EVENT main-line runner exited 10 min ago with no new marker"; exit 0; }
  [ "$gone_at3" -ge 10 ] && { echo "EVENT A-T3 runner exited 10 min ago with no new marker"; exit 0; }
  el=$(( $(date +%s) - T0 ))
  [ "$el" -ge 2700 ] && [ "$seen_at3" = 0 ] && { echo "EVENT A-T3 runner not started 45 min after dispatch"; exit 0; }
  [ -n "${MAIN_SENT:-}" ] && [ "$el" -ge 2700 ] && [ "$seen_main" = 0 ] && { echo "EVENT main-line runner not started 45 min after dispatch"; exit 0; }
  if [ $((i % 5)) -eq 0 ]; then
    bad=$(grep -E '^[0-9a-f]{64}  ' $S/$MAN | grep -v -E "${MAN_EXCLUDE:-^$}" | shasum -a 256 -c 2>/dev/null | grep -vc ': OK$')
    [ "$bad" != "0" ] && { echo "GUARD frozen file changed ($bad entries fail $MAN)"; exit 2; }
    [ "$(DSIG)" != "$D0" ] && { echo "GUARD frozen data file size/mtime changed (FREEZE_MANIFEST_R16_DATA)"; exit 2; }
  fi
  [ "$(find $M -type f -exec stat -f '%N %z %m' {} \; | sort | md5)" != "$MSIG" ] && { echo "GUARD approved model files modified"; exit 2; }
  w=$(find /Users/mac/AI-Biology-Lab ~/.cache/huggingface/hub -newer "$REF" \( -name '*.ckpt' -o -name '*.safetensors' -o -name 'adata_pred*.h5ad' -o -name 'adata_real*.h5ad' -o -name '*_de.csv' -o -name 'GSE295214*' \) -not -path "$P/$M/*" -not -path "$P/round16/.venv/*" 2>/dev/null | head -3)
  [ -n "$w" ] && { echo "GUARD model/archive/new-data file appeared: $w"; exit 2; }
  stray=$(find round16 -newer "$REF" -type f -not -path 'round16/.venv/*' -not -path 'round16/results/*' -not -path 'round16/audit/results/*' -not -path '*/__pycache__/*' -not -path 'round16/status/MAIN.txt' -not -path 'round16/status/A-T*.txt' -not -path 'round16/status/R16_CONTROLLER_NOTES.md' -not -path 'round16/status/R16_ARC_ARCHIVE_MANIFEST.txt' -not -path "${STRAY_EXCLUDE:-/nonexistent}" 2>/dev/null | head -3)
  [ -n "$stray" ] && { echo "GUARD write outside the allowed output areas: $stray"; exit 2; }
  h=$(git -C /Users/mac/AI-Biology-Lab rev-parse --short HEAD 2>/dev/null)
  [ "$h" != "$BASE_HEAD" ] && { echo "GUARD new git commit HEAD=$h ($(git -C /Users/mac/AI-Biology-Lab log -1 --format='%an: %s' | head -c 120))"; exit 2; }
  git -C /Users/mac/AI-Biology-Lab ls-files | grep -q -E '^projects/bioxarena-pro-scaffold/(paper|round16/shared/models)/' && { echo "GUARD paper/ or models tracked in git"; exit 2; }
  for f in $(find "$P" -newer "$REF" -type f -size +100M -not -path "${BIG_EXCLUDE:-/nonexistent}" -not -path '*/.git/*' -not -path "$P/$M/*" -not -path "$P/round16/.venv/*" 2>/dev/null | head -20); do
    git -C "$P" check-ignore -q "$f" || { echo "GUARD new non-ignored file >100MB: $f"; exit 2; }
  done
  echo "$(date '+%F %T') main_pid=${pm:-none} at3_pid=${pa:-none} res=$(find round16/results -type f 2>/dev/null | wc -l | tr -d ' ') audit_res=$(find round16/audit/results -type f 2>/dev/null | wc -l | tr -d ' ')" > "$HB"
  i=$((i + 1)); sleep 60
done
