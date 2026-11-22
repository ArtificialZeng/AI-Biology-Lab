# Controller independent check, round 10 (after CODE_FREEZE_R10 f8d2ebf)
Code: controller/checks/r10_fct_check.py (written from PREREG_C text). 1,000 joint bootstraps, controller seed 777 (official runs use 4,000 and seed 20261122), 100 permutations default_rng(k).
Not shared with A-T2/A-T3 before they report, to keep their runs independent.

| Line | N | FCT r_orth_cal [97.5% CI] | max mismatched CI hi | perm 97.5th pct | H-C1 | H-C2 |
|---|---|---|---|---|---|---|
| dld1 | 86 | +0.1401 [+0.0881, +0.1906] | +0.1073 (vcc2025_h1) | +0.1179 | pass | fail |
| cd4t_rest | 302 (284 half-valid) | +0.0560 [+0.0178, +0.0929] | +0.0728 (nadig_jurkat) | +0.0554 | pass | fail |

Controller reading (provisional until official runs and C3 cross-check): H-C1 2/2, H-C2 0/2 -> target-line-specific claim NOT CONFIRMED.
FCT keeps a positive amplitude-removed signal on both fresh lines, but its lower bound does not clear the mismatched-control and permutation thresholds.
