# Round 8 Checkpoint A — R8-0 Preprint Mandatory Revisions & Governance Pause

**Date**: 2026-11-21  
**R8-0 Commit (Local Only)**: `9c4aa6c` (`docs(R8-0): revise PREPRINT_DRAFT for all 6 Round 7 final audit items`)  
**Remote Push Policy**: **PAUSED** — zero `git push` to `origin` (`ArtificialZeng/AI-Biology-Lab`). Local commits only; local clone synchronization performed directly from SSD clone to Thunderbolt Sync Volume clone without touching GitHub.

---

## Verification of All 6 Mandatory R8-0 Fixes (`docs/PREPRINT_DRAFT.md`)

1. **Human-Only Authorship + AI-Assistance Statement**:
   - Author line replaced with `[Human PI]` only (per ICMJE and *Nature*/*Science*/*Cell* editorial policies).
   - Antigravity (Gemini) and Claude roles moved to **Section 2.4 (AI-Assistance & Verification Statement)** and **Section 6 (Acknowledgments)**.
2. **Title & Abstract Scope Calibration**:
   - Title revised to: *Pre-Registered Three-Atlas Evidence That Metric Choice Reverses Method Rankings in Single-Gene Perturbation Prediction* (removed "Outperform Linear and Graph-Neural Baselines").
   - Abstract phrasing changed from `"reflects a metric inversion"` to `"is partly explained by a metric inversion"`.
   - PC1 amplitude mechanism explicitly marked as a **hypothesis** pending direct decomposition testing (`R8-2`) across Abstract, Section 1, and Section 3.1.
3. **Human-AI Collaboration Category Qualification in Abstract**:
   - Directly next to the Codabench score (`0.0998` vs `0.0825`), Abstract now explicitly states: *"multi-round human-AI collaboration, which is not in the same comparison category as 2-hour single-shot autonomous agents"*.
4. **Reactome Enrichment Primary Background Switched to GMT-Annotated Universe (`M_bg = 1,452`)**:
   - Main results in Abstract and Section 3.5 now report `M_bg = 1,452` (`reports/t5_local_reactome_enrichment_summary.json`):
     - `S7 − S0` Top-150 (`128` in GMT): **mRNA Splicing - Major Pathway**, `38/145`, `2.97×`, **`BH FDR = 2.90e-08`** (`0/50` random controls).
     - `S5 − S2` Top-150 (`131` in GMT): **Mitochondrial translation**, `20/82`, `2.70×`, **`BH FDR = 6.98e-03`** (`0/50` random controls).
   - Full evaluated background (`M_bg = 1,784`: `3.12×`, `FDR = 1.18e-08` and `2.90×`, `FDR = 2.76e-03`) moved to Supplementary sensitivity note (Supplementary Table S1).
5. **GEARS Pseudobulk Adaptation Labeling & Methods Disclosure**:
   - Renamed throughout text and Table 1 to **`GEARS (official package, pseudobulk-adapted)`**.
   - Added explicit protocol details in **Section 2.3 (item 4)** and **Section 3.3 / Section 5**: `1` pseudobulk profile per perturbation, `10` jittered copies ($\mathcal{N}(0, 10^{-4})$) of the training control centroid as `"ctrl"`, and `skip_calc_de=True` (non-native single-cell training); removed unqualified claims that official GEARS generally underperforms pending native single-cell evaluation (`R8-1`).
   - Also aligned Table 1 `scGPT_FrozenEmb_Ridge` `de50_r`/`de50_sign` values directly to `r72_scgpt_probe_{replogle,nadig}_metrics.csv` (`0.2646`/`0.6536` Replogle; `0.3064`/`0.6792` Nadig).
6. **Codabench SHA-256 File Attribution**:
   - Explicitly specified in Abstract, Section 3.6, and Section 6 that `e2613a9b6616368ba5ca29edf132fe023761f0f72072b25176724207daa6a207` is the SHA-256 of **`submission.csv`** (`CODABENCH_FREEZE_MANIFEST.json`).

---

## Frozen Solver & Sync Invariants
- `git diff 3609811 -- projects/bioxarena-pro-scaffold/src/xatlas_panhuman_solver.py`: `0` lines.
- Zero `git push` executed since `cb3935f`.

CHECKPOINT_A_DONE
