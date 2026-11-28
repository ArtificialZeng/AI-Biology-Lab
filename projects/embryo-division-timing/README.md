# Embryonic cell division timing

This observational study tests how division labels and imaging cadence change inferred timing associations in published wild-type *Caenorhabditis elegans* lineage recordings. Antigravity performed the archival analyses and drafted the manuscript. Codex supervised the analysis and checked the final reproduction package.

## Manuscript and figures

| File | Content |
| --- | --- |
| [Manuscript v7](v7/manuscript_main.pdf) | Main results and methods |
| [Supplementary information](v7/supplementary_info.pdf) | Estimators, cohort definitions and robustness analyses |
| [Figure 1](v7/figures/Fig1_optical_boundary_control.pdf) | Matched optical boundary comparisons in two laboratories |
| [Figure 2](v7/figures/Fig2_cadence_estimand_and_correction.pdf) | Sampling cadence, estimator sensitivity and correction limits |
| [Supplementary Figure S1](v7/figures/FigS1_richards_transport.pdf) | Reanalysis of published timing comparisons |

## Main findings

- In the RIKEN archive, 52 embryos contributed 634 matched lineage quartets. Holding the outer endpoints fixed and swapping the central mother-end label raised the mean slope by 0.673 (95% CI 0.277 to 1.070). With nuclear outer endpoints, the corresponding contrast was 0.082 (95% CI −0.165 to 0.329). The CShaper comparison came from an independently acquired 17-embryo archive.
- In one archive, two admissible assignments of event times yield the same first and last coarse frames and opposite slope signs. This establishes a limit of the fitted association on the evaluated branch set.
- Recalculating the independence benchmark in three archives weakens the evidence for timing compensation based on the older sum-of-standard-deviations comparison. The published 18-embryo cohort's individual timestamps remain unavailable. [Richards et al. (2013)](https://doi.org/10.1016/j.ydbio.2012.11.034)

The manuscript is a research draft prepared for specialist review. Both optical analyses use membrane segmentation seeded by nuclear tracking. The inference concerns the recorded labels and observation cadence, and biological feedback requires separate evidence.

## Data and reproduction

Start with the Methods and Data and Code Availability sections of the manuscript, then use the branch inclusion and analysis specifications in the supplement. The optical comparison fixes the outer endpoints for matched sister branches and uses embryos as the paired units. The cadence checks use the archived timestamps and reproduce the first and last coarse frames on fixed branch support.

- [CShaper study and data availability](https://doi.org/10.1038/s41467-020-19863-x) · [CShaper data record](https://doi.org/10.6084/m9.figshare.12839315.v7)
- [Azuma et al. study](https://doi.org/10.3389/fbinf.2023.1082531) · [RIKEN SSBD archive](https://ssbd.riken.jp/repository/ssbd-repos-000236/)
- [AceTree curated timing archive](https://doi.org/10.5281/zenodo.16878148)
