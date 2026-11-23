# PREREG_H Addendum A: line identity for single-line files (pre-data)

Master Controller Claude | 2026-11-23 17:03 +0800

**Status:** written and committed before any VCC 2026 file exists in the public bucket and before any download (see the bucket check below). It changes no claim, test, threshold, endpoint, mask, source, multiplicity rule or seed. Parent: `PREREG_H.md` v2 (commit 55f4ae4, sha256 b6ec94ae…).

**Why:** §3.1 (v2) takes each cell's line from (a) an explicit cell-line column, else (b) a context column mapped by Arc's official table, else the line is excluded. That would wrongly exclude a release that puts each line in its own file without a line column, which is a natural release format. Both reader implementations followed §3.1 as written, which exposed the gap.

**Amendment to §3.1:** add rule (c) after (b):
- (c) If Arc's official release documentation (the release page, a README or a manifest file published with the data) states that a given file contains cells of exactly one named cell line, every cell in that file takes that line. The URL and the exact sentence or manifest entry are recorded in `R15_QC_COUNTS.md` before any metric. A file name alone is not enough.
- If a file has an explicit line column or a mapped context column, (a)/(b) take precedence; if they disagree with (c) for any cell, the affected lines are `EXCLUDED_FEASIBILITY`.
- Everything else in §3.1 is unchanged (never infer line from expression; drop Jurkat cells at ingestion).

**Bucket check at writing:** see `round15/status/PREREG_H_ADDENDUM_A_CHECK.txt`.
