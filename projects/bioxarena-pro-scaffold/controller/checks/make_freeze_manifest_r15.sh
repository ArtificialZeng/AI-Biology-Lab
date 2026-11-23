#!/bin/bash
# Writes round15/status/FREEZE_MANIFEST_R15.txt (sha256 of code and locked documents). Frozen BEFORE any VCC 2026 data exists:
# there is no Section 3 (prepared line data); line data will be hashed in a post-download addendum before any metric.
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
OUT=round15/status/FREEZE_MANIFEST_R15.txt
{
  echo "# FREEZE_MANIFEST_R15 — written by Claude (master controller) $(date '+%Y-%m-%d %H:%M %z'); pre-release freeze (no VCC 2026 data exists)"
  echo "# Section 1: code (sha256)"
  find round15/src round15/audit -maxdepth 1 -type f -name "*.py" | sort | xargs shasum -a 256
  shasum -a 256 controller/checks/r15_check.py controller/checks/r15_reader_xcheck/make_fixtures.py controller/checks/r15_reader_xcheck/run_xcheck.py
  shasum -a 256 round10/shared/scripts/prepare_r10_lines.py src/prepare_nadig2025_cache.py
  echo "# Section 2: locked documents, audits and pre-freeze verification records"
  shasum -a 256 round15/prereg/PREREG_H.md round15/prereg/PREREG_H_ADDENDUM_A.md round15/status/PREREG_H_COMMIT.txt \
    round15/status/PREREG_H_LOCK_CHECK.txt round15/status/PREREG_H_ADDENDUM_A_CHECK.txt round15/status/PREREG_H_AUDIT_R1.md \
    round15/shared/R15_FEASIBILITY_VCC2026.md round15/shared/PRIOR_ART_H.md round15/status/PRIOR_ART_H_CONTROLLER_FAN2026.md \
    round15/status/R15_READER_XCHECK_CONTROLLER_R1.md controller/checks/r15_reader_xcheck/XCHECK_REPORT_R2.json \
    round15/results/regression_r12_r14/REGRESSION_DIFF_SUMMARY.json round15/audit/regression_r12_r14/REGRESSION_REPORT_R12_R14.json \
    controller/checks/r15_dev/dev12_a549_ko.json controller/checks/r15_dev/dev14_a549_myllia.json
  echo "# Reference outside the repo (not verified against git): $(shasum -a 256 /Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt | cut -c1-64)  hgnc_complete_set.txt"
  echo "# Section 3: prepared line data — none at freeze (data not yet released); to be hashed in FREEZE_MANIFEST_R15_DATA before any metric"
} > "$OUT"
echo "entries: $(grep -vc '^#' "$OUT")"
echo "verify OK: $(grep -v '^#' "$OUT" | shasum -a 256 -c 2>&1 | grep -c ': OK$')"
echo "manifest sha256: $(shasum -a 256 "$OUT" | cut -c1-64)"
