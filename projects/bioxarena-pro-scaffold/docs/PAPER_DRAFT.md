# Benchmarking Out-of-Distribution Generalization in Biomedical ML Agents: From Baseline Reproduction to Cellular Perturbation Modeling

**Working Draft**  
**Current Milestone**: Level 1 (Rigorous Baseline Reproduction & Infrastructure Alignment)  
**Primary Target Task**: `perturbation/xatlas-orion-perturb`  

---

## 1. Baseline Engineering & Verification Audit

We implemented and verified multi-domain biomedical pipelines across chemical biology, disease diagnosis, and cellular perturbations.
Auditing against the official BioXArena-Pro graders confirms:
- Baseline scores match published reference baselines within statistical noise (gutcrc: 0.4453 vs 0.444; solubility: 0.7150 vs 0.711; cpjump: 0.6112 vs 0.613).
- Custom metrics can introduce substantial evaluation variance (e.g., GDSC scored 0.1199 under official `grade.py` vs 0.148 reference baseline). All evaluations are strictly mandated to execute official `tasks/*/grade.py` and `strict_group_mean_from_row_scores`.

---

## 2. Multi-Split Benchmarking on xatlas-orion-perturb

In the 64-task BioXArena-Pro benchmark, `xatlas-orion-perturb` evaluates genome-wide CRISPRi transcriptional response to held-out target genes in HCT116 and HEK293T cells.
Across 10 independent random splits (Seeds: 42, 0, 1, 2, 3, 7, 11, 99, 2024, 31337):
- **Constant Mean**: $0.0090 \pm 0.0040$ (Official reference: $0.009$).
- **30-kNN Co-expression**: $0.0413 \pm 0.0054$ (Official reference: $0.047$).
- **100-SVD + Ridge ($\alpha=10$)**: $0.0613 \pm 0.0040$ (Official reference: $0.051$).
- **Cosine Similarity Ridge ($\alpha=5.0$)**: $0.0655 \pm 0.0044$.
  - Mean delta over SVD+Ridge: $+0.0041 \pm 0.0023$.
  - With an official standard error of $\approx 0.01$ and noise threshold of $0.02$, the $+0.0041$ gain remains within statistical noise.
  - While seed 42 yielded 0.0804, multi-seed evaluation reveals that seed 42 represents an extreme upper-tail split.

We conclude that current performance represents Level 1 (Engineering Baseline Reproduction). Genuine Level 2 milestone requires demonstrating a 10-split mean gain exceeding the $0.02$ noise threshold over our internal baseline ($0.0613$), reaching $\ge 0.0813$.
