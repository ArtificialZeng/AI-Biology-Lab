#!/bin/bash
# Writes round12/status/FREEZE_MANIFEST_R12.txt (sha256 of code, locked documents, gate files and prepared line data).
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
OUT=round12/status/FREEZE_MANIFEST_R12.txt
{
  echo "# FREEZE_MANIFEST_R12 — written by Claude (master controller) $(date '+%Y-%m-%d %H:%M %z')"
  echo "# Section 1: code (sha256)"
  find round12/src round12/audit round12/shared/scripts -maxdepth 1 -type f -name "*.py" | sort | xargs shasum -a 256
  shasum -a 256 controller/checks/r12_check.py controller/checks/r12_paired_check.py
  echo "# Section 2: locked documents, rulings, QC and gate files"
  shasum -a 256 round12/prereg/PREREG_E.md round12/shared/R12_DOWNLOAD_PLAN.md round12/status/PREREG_E_COMMIT.txt \
    round12/status/RULING_R12_PREFREEZE.md round12/shared/R12_QC_COUNTS.md round12/status/CONFIRMATORY_LINES_R12.txt \
    round12/shared/R12_CHECKSUMS.sha256 round12/shared/hek293_feasibility.json \
    round12/status/RULING_R12_HEK293_UNITS.md round12/prereg/DEVIATIONS_R12.md \
    round12/status/CONFIRMATORY_LINES_R12_SUPERSEDED_v1.txt
  [ -s round12/audit/READER_XCHECK_R12.md ] && shasum -a 256 round12/audit/READER_XCHECK_R12.md
  echo "# Section 3: prepared line data (inputs to evaluation)"
  find round12/shared/lines -type f -not -path "*/.cache*" -not -name "_tmp_*" | sort | xargs shasum -a 256
} > "$OUT"
grep -vc "^#" "$OUT"
