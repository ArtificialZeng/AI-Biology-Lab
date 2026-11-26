#!/bin/bash
# Controller-only fetch of the ONE approved Arc archive file (Addendum A section 3: nadig_hepg2 adata_pred.h5ad),
# after R16_GO.txt. Pinned revision; parallel byte-range segments; exact size + Hugging Face LFS sha256 verified.
set -u
URL="https://huggingface.co/arcinstitute/ST-SE-Replogle/resolve/0fa8874b6e3bdf909b09c53c97ad5dc58b18a08f/zeroshot/hepg2/eval_final.ckpt/adata_pred.h5ad"
SIZE=3151145504
SHA=505fcfd3127caf9cfff2697b7322b64372433e0455cc510b7f20dbf736dfee22
DST=/Users/mac/.cache/r16_arc_archive/ST-SE-Replogle/zeroshot/hepg2/eval_final.ckpt/adata_pred.h5ad
SEGDIR=$(dirname "$DST")/.segments
SEG=$((64 * 1024 * 1024)); N=$(( (SIZE + SEG - 1) / SEG ))
mkdir -p "$SEGDIR"
fetch() {  # $1 = segment index
  local i=$1 s=$(( $1 * SEG )) e=$(( ($1 + 1) * SEG - 1 )) f="$SEGDIR/$(printf %04d $1)"
  [ $e -ge $SIZE ] && e=$((SIZE - 1))
  local want=$(( e - s + 1 ))
  for a in $(seq 1 20); do
    [ -f "$f" ] && [ "$(stat -f%z "$f")" -eq $want ] && return 0
    curl -sSLf --noproxy '*' --retry 5 --retry-delay 3 --retry-all-errors -r "$s-$e" -o "$f.part" "$URL" && mv "$f.part" "$f"
  done
  echo "SEGMENT_FAILED $i"; return 1
}
export -f fetch; export URL SIZE SEG SEGDIR
seq 0 $((N - 1)) | xargs -P 12 -I{} bash -c 'fetch {}' || { echo FETCH_FAILED; exit 1; }
cat $(for i in $(seq 0 $((N - 1))); do printf "%s/%04d " "$SEGDIR" $i; done) > "$DST.tmp" || exit 1
[ "$(stat -f%z "$DST.tmp")" -eq $SIZE ] || { echo "SIZE_BAD $(stat -f%z "$DST.tmp")"; exit 1; }
H=$(shasum -a 256 "$DST.tmp" | cut -d' ' -f1)
[ "$H" = "$SHA" ] || { echo "SHA_BAD $H"; exit 1; }
mv "$DST.tmp" "$DST" && chmod a-w "$DST" && rm -rf "$SEGDIR"
echo "$H  $DST  $SIZE"; echo "ARC_ARCHIVE_OK $(date '+%F %T')"
