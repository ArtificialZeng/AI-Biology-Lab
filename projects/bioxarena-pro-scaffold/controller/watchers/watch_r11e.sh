#!/bin/bash
# Wait for A-T2 R11_EVAL_DONE and A-T3 R11_AUDIT_DONE; alert on frozen-code changes. 4 h timeout.
B=/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
cd $B
for i in $(seq 1 960); do
  if ! grep -v '^#' round11/status/FREEZE_MANIFEST_R11.txt | grep -E '\.py$' | shasum -a 256 -c --quiet >/dev/null 2>&1; then echo "$(date +%H:%M) FROZEN_R11_CODE_CHANGED_ALERT"; grep -v '^#' round11/status/FREEZE_MANIFEST_R11.txt | grep -E '\.py$' | shasum -a 256 -c 2>&1 | grep -v ': OK$'; exit 0; fi
  a=$(grep -c '^R11_EVAL_DONE' round11/status/A-T2.txt 2>/dev/null); b=$(grep -c '^R11_AUDIT_DONE' round11/status/A-T3.txt 2>/dev/null)
  if [ "${a:-0}" -ge 1 ] && [ "${b:-0}" -ge 1 ]; then echo "$(date +%H:%M) R11_EVAL_AND_AUDIT_DONE"; exit 0; fi
  sleep 15
done
echo "$(date +%H:%M) TIMEOUT a=$a b=$b"
