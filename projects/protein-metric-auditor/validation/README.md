# Archived formula check

`archived_formula_review.json` preserves a separate computational audit: 2,331 array/cutoff cases, 240,051 enumerated row orders, and comparison to scikit-learn's tie-aware NDCG with the same normalized linear gains. It records the exact earlier source hash. Its observed finite-span-overflow robustness issue was subsequently guarded in the census evaluator and is also explicitly rejected by this project's sample script.

This archived check is evidence about those formulas and test cases, not a claim that the packaged scripts are unchanged from its source or have undergone external human peer review. The packaged sample independently compares against the pinned official evaluator and checks finite-data guards and order invariance.
