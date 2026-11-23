#!/bin/bash
# Writes round14/status/FREEZE_MANIFEST_R14.txt (sha256 of code, locked documents, gate files and prepared line data).
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
OUT=round14/status/FREEZE_MANIFEST_R14.txt
{
  echo "# FREEZE_MANIFEST_R14 — written by Claude (master controller) $(date '+%Y-%m-%d %H:%M %z')"
  echo "# Section 1: code (sha256)"
  find round14/src round14/audit -maxdepth 1 -type f -name "*.py" | sort | xargs shasum -a 256
  shasum -a 256 controller/checks/r14_check.py
  echo "# Section 2: locked documents, authorizations, QC and gate files"
  shasum -a 256 round14/prereg/PREREG_G.md round14/status/PREREG_G_COMMIT.txt round14/status/R14_GO.txt \
    round14/status/R14_LINES_AUTHORIZED.txt round14/shared/R14_QC_COUNTS.md round14/status/CONFIRMATORY_LINES_R14.txt \
    round14/shared/R14_CHECKSUMS.sha256 round14/shared/download_manifest_r14.json round14/shared/R14_FEASIBILITY_GSE337804.md \
    round14/audit/READER_XCHECK_R14.md
  echo "# Section 3: prepared line data (inputs to evaluation)"
  find round14/shared/lines -type f -not -path "*/.cache*" -not -name "_tmp_*" | sort | xargs shasum -a 256
} > "$OUT"
echo "entries: $(grep -vc '^#' "$OUT")"
echo "verify OK: $(grep -v '^#' "$OUT" | shasum -a 256 -c 2>&1 | grep -c ': OK$')"
echo "manifest sha256: $(shasum -a 256 "$OUT" | cut -c1-64)"
echo "data in git: [$(git -c http.proxy= -c https.proxy= ls-files round14/shared/raw round14/shared/lines round14/audit/reader_xcheck paper | head -3)]"
