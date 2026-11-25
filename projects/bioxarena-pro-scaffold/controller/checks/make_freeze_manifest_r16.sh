#!/bin/bash
# Writes round16/status/FREEZE_MANIFEST_R16.txt (sha256 of code, locked documents, environment and weight checksums).
# Line data are hashed separately in FREEZE_MANIFEST_R16_DATA.txt before R16_GO (the Nadig raw files sit on an external drive).
cd /Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold
OUT=round16/status/FREEZE_MANIFEST_R16.txt
{
  echo "# FREEZE_MANIFEST_R16 — written by Claude (master controller) $(date '+%Y-%m-%d %H:%M %z'); before any real-cell inference"
  echo "# Section 1: Round 16 code (main line round16/src, independent audit round16/audit, controller checks)"
  find round16/src round16/audit -maxdepth 1 -type f \( -name "*.py" -o -name "*.json" \) | sort | xargs shasum -a 256
  shasum -a 256 controller/checks/r16_check.py controller/checks/r16_check_adjudicate.py controller/checks/r16_download_models.py controller/checks/r15_check.py
  echo "# Section 2: frozen upstream code and gene catalog imported read-only (Rounds 9-12)"
  shasum -a 256 round9/p2/shared/COMMON_GENES.txt round10/src/data_loader_r10.py round10/src/fct_pipeline.py round10/src/metrics_r10.py \
    round10/shared/scripts/prepare_r10_lines.py round11/src/metrics_r11.py round12/src/metrics_r12.py
  echo "# Section 3: locked documents, approvals, schema, audits and environment"
  shasum -a 256 round16/prereg/PREREG_I.md round16/prereg/PREREG_I_ADDENDUM_A.md round16/status/PREREG_I_COMMIT.txt \
    round16/status/PREREG_I_LOCK_CHECK.txt round16/status/PREREG_I_ADDENDUM_A_CHECK.txt round16/status/R16_HUMAN_APPROVAL.txt \
    round16/status/R16_GATE2_SCHEMA.md round16/status/PREREG_I_AUDIT_R1.md round16/status/PREREG_I_AUDIT_AT3.md \
    round16/status/PREREG_I_V5_AUDIT_AT3.md round16/status/PRIOR_ART_I.md round16/status/PRIOR_ART_I_CONTROLLER_CHECK.md \
    round16/status/PRIOR_ART_I_ADD1.md round16/status/R16_VENV_FREEZE.txt round16/shared/R16_MODEL_CHECKSUMS.sha256 \
    round16/shared/R16_MODEL_DOWNLOAD_PLAN.json round16/shared/candidates/R16_CANDIDATES.md
  echo "# Reference outside the repo (not in git): $(shasum -a 256 /Users/mac/.cache/r9_p2_raw/hgnc_complete_set.txt | cut -c1-64)  hgnc_complete_set.txt"
  echo "# Model weights: 77 files under round16/shared/models (not in git), verified by R16_MODEL_CHECKSUMS.sha256: $(shasum -a 256 -c round16/shared/R16_MODEL_CHECKSUMS.sha256 2>/dev/null | grep -c ': OK$')/77 OK"
  echo "# Section 4: line data — hashed in FREEZE_MANIFEST_R16_DATA.txt before R16_GO"
} > "$OUT"
echo "entries: $(grep -vc '^#' "$OUT")"
echo "verify OK: $(grep -v '^#' "$OUT" | shasum -a 256 -c 2>&1 | grep -c ': OK$')"
echo "manifest sha256: $(shasum -a 256 "$OUT" | cut -c1-64)"
