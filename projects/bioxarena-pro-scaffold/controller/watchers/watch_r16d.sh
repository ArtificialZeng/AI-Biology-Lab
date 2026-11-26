#!/bin/bash
# R16 post-freeze watcher (controller). Exits on the first event or guard. Fixed baselines from 2026-11-26 02:15.
P=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cd "$P" || exit 9
BASE_HEAD=3c357c7
M=round16/shared/models
NADIG="/Volumes/雷雳传输/A顶会论文同步盘/A_BIG_MAC顶会论文写稿/AI自主研究/gemini/external_data/nadig2025/GSE264667_hepg2_raw_singlecell_01.h5ad"
REF=$(mktemp); touch -t 202611260215 "$REF"
NFILES=$(find $M -type f | wc -l | tr -d ' ')
MSIG=$(find $M -type f -exec stat -f '%N %z %m' {} \; | sort | md5)
i=0
while true; do
  [ -f "$NADIG" ] && { echo "EVENT external drive mounted, Nadig raw present at $(date '+%F %T')"; exit 0; }
  [ -e round16/status/R16_GO.txt ] && { echo "GUARD R16_GO.txt appeared (only the controller may create it)"; exit 2; }
  [ -d round16/results ] && [ -n "$(find round16/results -type f 2>/dev/null | head -1)" ] && { echo "GUARD round16/results has files before R16_GO"; exit 2; }
  if [ $((i % 5)) -eq 0 ]; then
    bad=$(grep -v '^#' round16/status/FREEZE_MANIFEST_R16.txt | shasum -a 256 -c 2>/dev/null | grep -vc ': OK$')
    [ "$bad" != "0" ] && { echo "GUARD frozen file changed ($bad entries fail FREEZE_MANIFEST_R16)"; exit 2; }
  fi
  [ "$(find $M -type f | wc -l | tr -d ' ')" != "$NFILES" ] && { echo "GUARD file count changed in $M"; exit 2; }
  [ "$(find $M -type f -exec stat -f '%N %z %m' {} \; | sort | md5)" != "$MSIG" ] && { echo "GUARD approved model files modified"; exit 2; }
  w=$(find /Users/mac/AI-Biology-Lab ~/.cache/huggingface/hub -newer "$REF" \( -name '*.ckpt' -o -name '*.safetensors' -o -name 'adata_pred*' -o -name 'adata_real*' -o -name 'eval_final*' -o -name 'GSE295214*' \) -not -path "$P/$M/*" -not -path "$P/round16/.venv/*" 2>/dev/null | head -3)
  [ -n "$w" ] && { echo "GUARD model/eval/new-data file appeared: $w"; exit 2; }
  [ -n "$(ls -d /opt/anaconda3/lib/python3*/site-packages/arc_state* 2>/dev/null)" ] && { echo "GUARD arc-state installed in base env"; exit 2; }
  h=$(git -C /Users/mac/AI-Biology-Lab rev-parse --short HEAD 2>/dev/null)
  [ "$h" != "$BASE_HEAD" ] && { echo "GUARD new git commit HEAD=$h"; exit 2; }
  git -C /Users/mac/AI-Biology-Lab ls-files | grep -q -E '^projects/bioxarena-pro-scaffold/(paper|round16/shared/models)/' && { echo "GUARD paper/ or models tracked in git"; exit 2; }
  big=$(find "$P" -newer "$REF" -type f -size +100M -not -path '*/.git/*' -not -path "$P/$M/*" -not -path "$P/round16/.venv/*" 2>/dev/null | head -3)
  [ -n "$big" ] && { echo "GUARD new file >100MB: $big"; exit 2; }
  i=$((i + 1)); sleep 60
done
