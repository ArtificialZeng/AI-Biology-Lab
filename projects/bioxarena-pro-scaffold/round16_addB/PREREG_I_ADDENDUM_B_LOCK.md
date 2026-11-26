# PREREG_I Addendum B (Occam): controller lock record

Claude (master controller), 2026-11-27 01:15 +0800.

## What is locked
- **Locked text:** `round16_addB/PREREG_I_ADDENDUM_B_DRAFT_v2.md`, sha256 `537467db51b00b8e178ad5d7638e25b5d809f350c8feae229d86cdb96acb281e`. From this commit on it is Addendum B. The earlier draft and the review stay as history.
- **Companion prior art:** `round16_addB/PRIOR_ART_I_ADD2_v2.md`, sha256 `2e432d4779a0cb1ee463951cb435fc93732ddd05dd3a6365cadc968dd9faeaa1`.
- **State at lock:**
  - No Occam weight, asset, code archive or interpreter has been downloaded.
  - No Occam inference has run, and no GSE295214 file exists locally.
  - No R16 metric appears in either file (controller scan).
  - The R16 primary analysis is unchanged. Addendum B adds a secondary model only, and does not touch the PREREG_I Bonferroni family.

## Controller review of v2
- v2 addresses review R1 items P1–P6 (`ADDB_CONTROLLER_REVIEW_R1.md`).
- **Controller spot-check:** `xinyizhanglab/perturbation-decomposition/models/state/predict.py` (GitHub raw) evaluates saved STATE predictions produced by the original STATE pipeline. It has no control-cell swap. This matches v2 §2.1.
- **Still not read:** the Molina & Zhang preprint text (Cloudflare, not bypassed). v2 treats it as a possibly substantial overlap. That stays true until the text is read.

## Binding clarifications (part of the lock)
- **L1. Nadig readers.** In Gate 5B and in every Nadig read, the perturbation labels (`obs["gene"]`) and gene names are legacy AnnData categoricals. Decode them with the same legacy-categorical branch as DEVIATION_R16_D3; the integer codes are never labels. Gate 5B must find 611 cohort perturbations and 4,976 non-targeting controls, or it halts.
- **L2. Approval scope.** The human approval (`R16_ADDB_HUMAN_APPROVAL.txt`) must name:
  - the weights bundle (minimal 2.55 GB, or the full 2.93 GB) under licence "other / see-repository";
  - the source archive of commit `9ddb890e…`;
  - the Python interpreter install for `round16_addB/.venv` (Python 3.8.19, or the v2 §6.1 fallback chain).

  Nothing is downloaded or installed before that file exists.

## Next gates (unchanged from v2 §6)
1. Human approval.
2. Downloads, hashed against v2 §1.3.
3. Both parties build the venv and write their adapters.
4. Gate 5A, then Gate 5B.
5. Code freeze.
6. `R16_ADDB_GO.txt` (controller only).
7. Run, but only after A-T3's R16 run has finished: no concurrent MPS.
