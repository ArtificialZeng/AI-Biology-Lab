# Benchmarking Out-of-Distribution Generalization in Biomedical ML Agents: From Baseline Engineering to Transcriptional Perturbation SOTA

**Working Document**  
**Current Milestone**: Level 1 (Official Baseline Reproduction & Infrastructure Alignment)  
**Target Focus**: `perturbation/xatlas-orion-perturb`  

---

## 1. Baseline Engineering & Verification Audit

We implemented and verified multi-domain biomedical pipelines across chemical biology, disease diagnosis, and cellular perturbations.
Auditing against the official BioXArena-Pro graders confirms:
- Baseline scores match published reference baselines (gutcrc: 0.4453 vs 0.444; solubility: 0.7150 vs 0.711; cpjump: 0.6112 vs 0.613).
- Custom metrics can introduce substantial evaluation variance (e.g., GDSC scored 0.1199 under official `grade.py` vs 0.148 reference baseline). All evaluations are strictly mandated to execute official `tasks/*/grade.py`.

## 2. Core Opportunity: xatlas-orion-perturb

In the 64-task BioXArena-Pro benchmark, `xatlas-orion-perturb` evaluates genome-wide CRISPRi transcriptional response to held-out target genes.
- Current 22-agent ceiling: 8.25 / 100 ($r \approx 0.083$, median 6.71).
- Reference baselines: 0.047 (kNN co-expression) and 0.051 (100-SVD + Ridge).
- We focus research effort on uncovering unsupervised target column structure and pathway co-perturbation manifold representations to surpass the 0.083 frontier ceiling towards > 0.12.
