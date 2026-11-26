# Round 16 verdict — Groups 1–2 (pre-registered retrospective audit of Arc STATE)

Written by Claude (master controller), 2026-11-27 07:30 +0800. This covers PREREG_I v5 + Addendum A, Groups 1–2 and the feng_ipsc sensitivity line.

Not yet covered:
- Group 3 (gse295214_h1, prospective). It is running under `R16_G3_GO.txt`, and A-T3's blind run follows.
- Addendum B (Occam, secondary model). Its CPU gates have passed; it still needs the MPS parity check, A-T3's adapter, the code freeze and its GO.

## 1. What was tested
Arc STATE (ST-SE-Replogle with SE-600M embeddings, primary; ST-HVG-Replogle, secondary), used zero-shot.

**H-I1 (control specificity).** Do predictions on a target line depend on the line's own control cells? Predictions with matched controls are compared against predictions with control cells swapped in from other lines (TOST δ = 0.04, 97.5% CIs).

**H-I2 (baseline comparison).** Is the model more accurate than the FCT baseline? FCT is the closed-form transfer of mean responses from other lines.

**Metric.** r_orth_cal: correlation with the observed perturbation response, orthogonalized against the transfer mean and calibrated by the half-split noise ceiling. Joint perturbation bootstrap, 4,000 draws, seed 20261601.

**Batch token.** Primary neutral setting fixed_index_zero (batch_idx = 0 for every cell). This is the arc-state 0.11.1 fallback for unseen batches.

**Cohort.**
- Group 1: nadig_hepg2 and nadig_jurkat. CRISPRi; these contexts are partly seen by the zero-shot checkpoints.
- Group 2: a549_ko (CRISPRko), dld1 and cd4t_rest (CRISPRi). All unseen.
- feng_ipsc: sensitivity line only.

These cohorts had been analysed with closed-form baselines in Rounds 9–14. This is therefore a **pre-registered retrospective model audit**, not a prospective test.

## 2. Verification (all passed)

| Layer | Parties | Result |
|---|---|---|
| Control-cell barcodes (primary + resample pools) | main line, A-T3 (blind), controller dry run | exact |
| SE-600M embeddings, 11 pools | main line vs A-T3 | bitwise identical |
| Predicted Δ, all arms, both models, 5 lines | main line vs A-T3 (pass 1, STATE_HVG rerun after D8, pass 2) | bitwise identical (max \|diff\| 0) |
| Primary metrics, per-source categories, H-I1/H-I2 calls | controller recompute (`r16_check.py`) on A-T3 arrays vs main summary | worst \|diff\| 5.6e-16 |
| STATE_HVG metrics and calls (A-T3's own scoring code) | A-T3 vs main | 490/490, worst 1.4e-17 |
| G4 2×2, (17, 0) symmetry check, single-index spread, HVG Perm100 | A-T3 pass 2 vs main | 1,213/1,213, worst 5.6e-17 |
| Neutral-batch panel (main line only) | controller recompute vs main | 240/240, worst 4.7e-16 |
| feng_ipsc sensitivity (main line only) | controller recompute vs main | match |

Deviations D1–D8 (records in `round16/status/DEVIATION_R16_*.md`) were all implementation fixes found before any metric was compared. None changed the preregistration, the data or the models. D8 (A-T3's Nadig HVG normalization totals) is the only one found by array comparison, and it changed no call.

These are computational checks between independent implementations. They are not external peer review.

## 3. Primary results (STATE_SE, fixed_index_zero)

| Line | Group | N_eval | Matched r_orth_cal [97.5% CI] | FCT | STATE_SE − FCT [97.5% CI] | H-I1 | H-I2 |
|---|---|---:|---|---:|---|---|---|
| nadig_hepg2 | 1 | 740 | 0.1053 [0.0980, 0.1128] | 0.2015 | −0.0961 [−0.1099, −0.0832] | LINE_MATCHED_SPECIFIC | FCT_SUPERIOR_TO_MODEL |
| nadig_jurkat | 1 | 873 | 0.0383 [0.0324, 0.0440] | 0.1475 | −0.1092 [−0.1250, −0.0931] | MIXED_OR_INCONCLUSIVE | FCT_SUPERIOR_TO_MODEL |
| a549_ko | 2 (KO) | 403 | 0.1143 [0.0980, 0.1309] | 0.2134 | −0.0991 [−0.1245, −0.0728] | MIXED_OR_INCONCLUSIVE | FCT_SUPERIOR_TO_MODEL |
| dld1 | 2 | 75 | 0.1024 [0.0811, 0.1238] | 0.1682 | −0.0658 [−0.1049, −0.0264] | MIXED_OR_INCONCLUSIVE | FCT_SUPERIOR_TO_MODEL |
| cd4t_rest | 2 | 245 | 0.0366 [0.0105, 0.0630] | 0.0583 | −0.0217 [−0.0491, 0.0060] | MIXED_OR_INCONCLUSIVE | INDISTINGUISHABLE |
| feng_ipsc (sensitivity) | — | 88 | 0.0158 [−0.0123, 0.0458] | 0.1302 | −0.1144 [−0.1529, −0.0742] | (sensitivity) | (sensitivity) |

**Group adjudication** (by modality, k/K, no pooling):
- H-I1:
  - Group 1 CRISPRi: LINE_MATCHED_SPECIFIC 1/2.
  - Group 2 CRISPRi: MIXED_OR_INCONCLUSIVE 2/2.
  - Group 2 CRISPRko: MIXED_OR_INCONCLUSIVE 1/1.
- H-I2: under the primary neutral batch setting, FCT is superior to STATE_SE on 4/5 lines and indistinguishable on cd4t_rest. STATE_SE is superior on no line.

**Reading of the primary results.**
- STATE_SE's predictions depend on the target line's own controls in a line-specific way only on nadig_hepg2, a context partly seen in training.
- On unseen lines, swapping in another line's controls gives mixed effects that differ by source. On a549_ko, for example, controls from nadig_hepg2 predict *better* than a549_ko's own controls (paired −0.087), while controls from nadig_jurkat predict worse (+0.146).
- Under the official unseen-batch fallback, a closed-form transfer baseline beats the model on 4 of 5 lines.

## 4. Secondary results (pre-registered PREREG_I §4; exploratory in status)

**STATE_HVG.**
- Its matched accuracy exceeds STATE_SE's on 4 of 5 lines; on a549_ko it is 0.2373 against 0.1143.
- On a549_ko it is LINE_MATCHED_SPECIFIC, and STATE_HVG − FCT is INDISTINGUISHABLE (+0.0239 [−0.0128, 0.0635]).
- FCT is superior on the other 4 lines.
- The HVG Perm100 gene-label null is exceeded by the matched 97.5% lower bound on 4/5 lines, but not on cd4t_rest (q975 0.0024, lower bound −0.0210).

**Neutral-batch panel** (STATE_SE; identical batch input for matched and mismatched arms). This changes how the primary H-I2 result should be read:

| Line | fixed_index_zero (primary) | bypass_zero | mean_batch_emb | multi_batch_ens8 | FCT |
|---|---:|---:|---:|---:|---:|
| nadig_hepg2 | 0.1053 (FCT sup.) | 0.2341 (**model sup.**, +0.0326 [0.0139, 0.0508]) | 0.2022 (indist.) | 0.2008 (indist.) | 0.2015 |
| nadig_jurkat | 0.0383 (FCT sup.) | 0.0419 (FCT sup.) | 0.0763 (FCT sup.) | 0.0576 (FCT sup.) | 0.1475 |
| a549_ko | 0.1143 (FCT sup.) | 0.2085 (indist.) | 0.1846 (indist.) | 0.2874 (**model sup.**, +0.0740 [0.0412, 0.1090]) | 0.2134 |
| dld1 | 0.1024 (FCT sup.) | 0.0922 (FCT sup.) | 0.1069 (FCT sup.) | 0.1131 (FCT sup.) | 0.1682 |
| cd4t_rest | 0.0366 (indist.) | 0.0084 (FCT sup.) | 0.0094 (FCT sup.) | 0.0563 (indist.) | 0.0583 |

**Single trained batch index**, matched controls, Group 2:
- a549_ko: r_orth_cal 0.0863 (index 17) to 0.2426 (index 35), with non-overlapping 97.5% CIs; index 0 gives 0.1143.
- dld1: 0.0425 to 0.1024.
- cd4t_rest: −0.0008 to 0.0366.

**G4 2×2 (control cells × batch token).**
- Group 1, native gem_group batches: almost all of the matched-control advantage comes from the control cells. Δ_ctrl is +0.2904 (hepg2) and +0.2333 (jurkat); the batch share is about 2–3%.
- With native batches, STATE_SE on nadig_hepg2 reaches 0.3522, against FCT 0.2015.
- Group 2 (B_M = 0 vs B_S = 17): the batch-token swap alone moves accuracy by +0.0280 (a549_ko), +0.0598 (dld1) and +0.0374 (cd4t_rest). That is often as large as or larger than the control swap; the batch share is 0.16–0.98 across sources.

**Secondary reading.**
- On unseen contexts, STATE's zero-shot accuracy depends strongly on an arbitrary input with no correct value for an unseen line, the batch index. On two lines that dependence is large enough to reverse the comparison with a simple baseline.
- On partly seen contexts, the control cells carry the line-specific signal.
- These are secondary, retrospective, post-Gate-2 analyses on 3 unseen lines and 4 indices. **They are hypothesis-generating.** Round 17 (PREREG_J, locked by hash commitment in commit 9789b0b before any Round 17 inference) tests them confirmatorily:
  - all 56 batch indices;
  - three fresh confirmatory lines: cd4t_stim48, gse295214_h1, pfizer_telohaec_basal;
  - cross-fitted estimators and pre-specified label-free index rules.

## 5. Limitations
- Retrospective cohorts (see §1). feng_ipsc is a sensitivity line only.
- One model family. Addendum B (Occam) is pending, and Occam has no batch encoder.
- The neutral-batch panel and the G4 decomposition are STATE_SE only.
- Group 1 contexts are partly seen by 3 of the 4 zero-shot checkpoints (PREREG_I M5.2 disclosure).
- dld1 has N_eval = 75, so its CIs are wide.
- Results concern zero-shot use of the released checkpoints with the official unseen-batch fallback. They are not a statement about fine-tuned or few-shot use.
- Close prior work (PRIOR_ART_I, PRIOR_ART_I_ADD2_v2) includes Molina & Zhang 2026, CELLAUDIT, and Zhang et al. arXiv:2611.05849 §C.5 (control + batch probe). The Molina & Zhang full text has not been read (Cloudflare challenge, not bypassed), so overlap with it may be substantial. No priority claim is made.
- The source model and data are credited in PRIOR_ART_I and the preregistration: Adduri et al. 2025 for STATE; Nadig et al. 2025; and the data sources listed for each line.
