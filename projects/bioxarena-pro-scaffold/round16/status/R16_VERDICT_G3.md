# Round 16 verdict — Group 3 (prospective holdout gse295214_h1)

Written by Claude (master controller), 2026-11-27 09:30 +0800. Covers PREREG_I §2.4, Group 3.

**The line.** GSE295214, H1 hESC, CRISPRi 10x 5′; Pan, Saunders, Replogle, Weissman & Zhuang, bioRxiv 10.64898/2026.05.16.725005.
- The line was pre-registered before download.
- Download and counts-only QC: PASS_GROUP3_PRIMARY, N_eval 406, |G_eval^R16| 1,837.
- Data freeze: `FREEZE_MANIFEST_R16_DATA_G3.txt`.
- The GO came after the PREREG_J lock, so no STATE inference on this line preceded either preregistration.

It is reported strictly separately from Groups 1–2. Disclosures:
- it shares the H1 lineage with R9 vcc2025_h1, which is excluded from the primary FCT pool;
- it shares laboratory provenance with the Replogle training lines.

## Verification (all passed)

| Layer | Parties | Result |
|---|---|---|
| Predicted Δ, all arms, both models | main line vs A-T3 (blind) | bitwise identical (max \|diff\| 0) |
| Primary metrics and calls | controller recompute on A-T3 arrays vs main summary | worst \|diff\| 5.6e-16 |
| STATE_HVG metrics and calls (A-T3's own code) | A-T3 vs main | 109/109, worst 1.4e-17 |
| G4, symmetry, single-index spread, Perm100 | A-T3 vs main | 433/433, worst 1.4e-17 |
| Neutral-batch panel (main line only) | controller recompute vs main | 33/33, worst 1.1e-16 |

## Results

| Model / setting | Matched r_orth_cal [97.5% CI] | Model − FCT [97.5% CI] | Call |
|---|---|---|---|
| STATE_SE, fixed_index_zero (primary) | 0.0220 [0.0130, 0.0313] | −0.1297 [−0.1459, −0.1139] | H-I2 FCT_SUPERIOR_TO_MODEL; H-I1 MIXED_OR_INCONCLUSIVE |
| STATE_SE, bypass_zero | 0.0424 | −0.1094 [−0.1271, −0.0911] | FCT superior |
| STATE_SE, mean_batch_emb | 0.0369 | −0.1149 [−0.1318, −0.0974] | FCT superior |
| STATE_SE, multi_batch_ens8 | 0.0351 | −0.1167 [−0.1340, −0.0990] | FCT superior |
| STATE_HVG (secondary) | 0.0349 [0.0229, 0.0475] | −0.1169 [−0.1340, −0.1003] | FCT superior; H-I1 MIXED |

- **FCT** (primary pool, excluding vcc2025_h1) is 0.1518 [0.1397, 0.1646]. With vcc2025_h1 included (sensitivity), STATE_SE − FCT is −0.1356.
- **H-I1 per source** (STATE_SE): within δ = 0.04 on all six sources.
  - matched better: nadig_hepg2 (+0.0174), a549_ko (+0.0066), cd4t_rest (+0.0053);
  - equivalent: nadig_jurkat and dld1;
  - mismatched better: feng_ipsc (−0.0041).
- **Single trained batch index** (STATE_SE): index 0 gives 0.0220, index 1 0.0249, index 17 0.0083, index 35 0.0364.

## Reading
- On a prospective, never-analysed cell context (H1 hESC), zero-shot STATE has little signal: r_orth_cal is about 0.02–0.04 under every batch setting, against the closed-form transfer baseline's 0.15. No neutral batch setting changes the H-I2 call.
- Batch-index dependence is present in relative terms (index 17 at 0.008 vs index 35 at 0.036), but it sits on a low-signal floor. This is the situation the Round 17 J2 detectable-signal guard was written for. Round 17's 56-index analysis of gse295214_h1 uses a fresh control sample (seed 20261701) and the locked PREREG_J rules. These Round 16 numbers do not enter Round 17.
- With Group 3, the Round 16 primary H-I2 result is: FCT superior on 5 of 6 primary lines (Groups 1–3), indistinguishable on cd4t_rest, and STATE superior on none, under the official batch fallback.

## Limitations
Those of `R16_VERDICT_G12.md` §5 apply. In addition, H1 is a pluripotent context far from the K562/RPE1/HepG2/Jurkat training lines, so the weak signal may reflect domain distance rather than anything batch-specific.
