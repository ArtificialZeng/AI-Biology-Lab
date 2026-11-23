#!/bin/bash
# Writes round13/status/FREEZE_MANIFEST_R13.txt (sha256 of code, locked documents, gate files and prepared line data).
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
OUT=round13/status/FREEZE_MANIFEST_R13.txt
{
  echo "# FREEZE_MANIFEST_R13 — written by Claude (master controller) $(date '+%Y-%m-%d %H:%M %z')"
  echo "# Section 1: code (sha256)"
  find round13/src round13/audit -maxdepth 1 -type f -name "*.py" | sort | xargs shasum -a 256
  shasum -a 256 controller/checks/r13_check.py
  echo "# Section 2: locked documents, authorizations, QC and gate files"
  shasum -a 256 round13/prereg/PREREG_F.md round13/status/PREREG_F_COMMIT.txt round13/status/R13_GO.txt \
    round13/status/R13_LINES_AUTHORIZED.txt round13/shared/R13_QC_COUNTS.md round13/status/CONFIRMATORY_LINES_R13.txt \
    round13/shared/R13_CHECKSUMS.sha256 round13/shared/download_manifest_r13.json
  [ -s round13/audit/READER_XCHECK_R13.md ] && shasum -a 256 round13/audit/READER_XCHECK_R13.md
  [ -s round13/status/RULING_R13_READER.md ] && shasum -a 256 round13/status/RULING_R13_READER.md
  echo "# Section 3: prepared line data (inputs to evaluation)"
  find round13/shared/lines -type f -not -path "*/.cache*" -not -name "_tmp_*" | sort | xargs shasum -a 256
} > "$OUT"
grep -vc "^#" "$OUT"
