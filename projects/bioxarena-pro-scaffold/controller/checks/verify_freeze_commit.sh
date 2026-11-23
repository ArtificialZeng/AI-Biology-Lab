#!/bin/bash
# Usage: bash controller/checks/verify_freeze_commit.sh <manifest_path>
# Checks HEAD == origin/main, and that every Section 1–2 manifest entry matches the committed blob at HEAD.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
G="git -c http.proxy= -c https.proxy="
$G fetch -q origin
echo "HEAD=$($G rev-parse HEAD) origin=$($G rev-parse origin/main)"
P=projects/bioxarena-pro-scaffold; ok=0; bad=""
while read h f; do c=$($G show "HEAD:$P/$f" 2>/dev/null | shasum -a 256 | cut -c1-64); if [ "$c" = "$h" ]; then ok=$((ok+1)); else bad="$bad $f"; fi; done < <(awk '/^# Section 3/{exit} !/^#/{print $1, $2}' "$1")
echo "committed-match OK=$ok BAD=[$bad]"
echo "data/paper in git: [$($G ls-files round14/shared/raw round14/shared/lines round14/audit/reader_xcheck paper | head -3)]"
echo "manifest sha256: $(shasum -a 256 "$1" | cut -c1-64); entries $(grep -vc '^#' "$1")"
