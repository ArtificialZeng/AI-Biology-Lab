# Controller review R1: PREREG_I_ADDENDUM_B_DRAFT.md and PRIOR_ART_I_ADD2.md

Claude (master controller), 2026-11-27 00:50 +0800.
Files reviewed:
- `PREREG_I_ADDENDUM_B_DRAFT.md` (sha256 6117709a…)
- `PRIOR_ART_I_ADD2.md` (sha256 9362e83c…)

## Verdict
The design is sound and close to lockable:
- Targets and strata are right: Group 2 only; the Nadig lines are seen sources.
- The one control-swap arm per source fits a model with no cell-line or batch ID.
- Reusing the R16 barcodes and the full (N, G) parity with STATE_SE make the paired contrast clean.
- The adapter gate against the authors' own HepG2 numbers is the right kind of gate.

The prior-art file missed close work. That must be fixed before lock, together with the items below. Send a v2. Do not lock it yourself.

## Facts the controller checked
- Hardware is Apple M3 Pro with 36 GB (`sysctl`), as the draft says.
- Local HepG2 file: `/Users/mac/.cache/r9_p2_raw/nadig2025/GSE264667_hepg2_raw_singlecell_01.h5ad`, sha256 `e1ad7c3c…dadb52e`. This matches the draft.

## Required changes

### P1. Prior art: three close items are missing and one is unsupported
1. **Molina & Zhang 2026.** "Perturbation response decomposition enables biologically aligned generalization to unseen perturbations and cellular contexts", bioRxiv 10.64898/2026.07.24.740459, v1 2026-07-27, CC BY (AITHYRA, Vienna). It is the closest item found so far.
   - From the abstract (bioRxiv API): "the perturbation and cell-line specific components are high-dimensional and cannot be recovered from control expression".
   - A search-engine summary says the paper also evaluates STATE zero-shot and finds it recovers the global template but almost no perturbation-specific or target-interaction component. **The controller could not verify this:** bioRxiv full text gives HTTP 403, and Europe PMC (PPR1286577) has no full text.
   - Try only legitimate routes: the author's or lab's PDF, a code repository, a JATS copy outside Cloudflare. Do not bypass any challenge page. If you cannot read it, record 未能读取 and treat it as a possibly substantial overlap.
   - Either way, add a section that sets it against R16. Their claim is observational (a decomposition of measured responses). R16 is an intervention on frozen model inputs with FCT as the positive control. R12 K1 (target controls beat mismatched controls for FCT on A549) sits in tension with "cannot be recovered from control expression". Say how the two can both hold (metric, orthogonalization, what "recover" means), or say that they conflict.
   - It is cited by COMPASS (Liang & Singh 2026, Europe PMC PPR1293149). Classify COMPASS too.
2. **Li et al. 2026, CELLAUDIT.** "Discover, Falsify, Revise: Auditing Input-Use Claims from Source Code to Predictive Contribution in Agent-Discovered Cell Models", arXiv 2609.27234, v1 2026-09-23.
   - It does frozen-checkpoint input replacement: compound and dose replacement, control-profile replacement within a plate family, and refits on disjoint control-well banks. It covers agent-generated predictors on BBBC047 and sci-Plex.
   - Abstract: "predictive generalization need not imply generalization of input-use claims".
   - No cross-cell-line control swap and no STATE or Occam, as far as the controller's reader could read. Read the whole paper to confirm.
   - It is a direct methodological precedent for "frozen model, replace one input, see whether the prediction depends on it". Cite it in §7.2 and the overlap verdict, and drop any wording that implies the replacement design itself is new.
3. **"Explicit Mechanistic Operators in Predictive Virtual Cells: Evidence, Failure Modes, and a Minimum Falsification Framework for Single-Cell Perturbation Prediction"** (MDPI Cells, doi 10.3390/cells15181711, 2026). Read and classify.
4. **"The Construct Is Part of the Target" (§4 item 3).** It is described as testing "reagent/construct-level sensitivity of released STATE predictions" from secondary search summaries. That breaks the file's own standard, and the controller's search found no such title. Give a verifiable identifier (OpenReview id or URL, authors) or remove the characterization. The same applies to "FalsifyVC" and "PerturbCausal" (no source found under that name).

### P2. Environment and code provenance
- Do not install anything into `round16/.venv`; it is frozen for R16.
- Pre-specify a separate `round16_addB/.venv` from Occam's own pins (pyproject: Python 3.8, torch 2.4.1, and the rest of its locked dependencies). If that cannot be built on macOS arm64, give the fallback now: the exact Python and torch versions, plus the cross-version tolerance gate you cite (`CROSS_VERSION_ATOL 1e-4`, `RTOL 1e-3`).
- Occam code comes in only after approval, as the source archive of the pinned commit `9ddb890e…` (code only, no LFS), with its sha256 recorded.
- Both parties use Occam's own model classes from that commit, unchanged. Each party writes its own adapter independently: input preparation, gene and perturbation mapping, delta, export. Say this in §6.2.
- Add both the weights (licence "other / see-repository") and the code archive to the human approval item.

### P3. Gate 5B fidelity
- Reconcile the HepG2 held-out count: §2.1 says 621 test perturbations, §5.2 and cohort.json say 611.
- If `splithalf/splits_manifest.json` or `cohort.json` holds the authors' cell-level split assignments, Gate 5B must use them exactly. You have read these files, so state now which case applies. Generate new splits only if the authors' assignments are not available; the ±0.015 tolerance is meant to cover only that case.
- Read HepG2 from the local copy above, not from the sync share.

### P4. Gene-coverage confound in the swap
- cd4t_rest is missing 413 of the 6,983 RHJ axis genes and a549_ko is missing 33; both are zero-filled.
- When either line is the target or a mismatched source, swapping controls changes gene coverage as well as control content. That can make a mismatched arm look worse, or better, for a reason unrelated to context use.
- Add a pre-registered sensitivity arm, "pairwise common axis". For each (target, source) pair, zero in both arms every input gene missing in either pool, then run matched and mismatched on that common axis. The primary analysis stays as drafted. Report this arm next to every pair that involves cd4t_rest or a549_ko.

### P5. Exact forward definition
- §3.1 Eq. 13 writes `E(x̄_ctrl)` (an encoding of a mean), but §3.4 averages per-cell outputs.
- State which one the official API does for the matched prediction and use exactly that. If the encoding is per cell, drop the bar.
- State that inference is deterministic: no sampling, the softplus mean output, dropout off. If any stochastic step exists, give its seed.

### P6. GPU rule
- Keep the no-concurrent-MPS rule from §6.1.
- Add: if an MPS run shows no progress for 10 minutes (no new log line and no CPU time), both parties switch to CPU float32 for the whole Addendum B run. Fix this rule now, not after the fact.

## What does not change
- The R16 primary analysis, PREREG_I, and Addendum A do not change.
- Nothing is downloaded, and Addendum B is not locked, until the human project lead approves and the controller locks.
