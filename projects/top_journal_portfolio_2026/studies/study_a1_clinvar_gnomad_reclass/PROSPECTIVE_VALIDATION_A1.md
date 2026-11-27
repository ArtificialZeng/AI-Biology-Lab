# Prospective Unseen Validation Protocol — Study A-1 (Locked Prior to Data Extraction)

- **Lock Timestamp**: `2026-11-25 09:51 (+08:00)` (strictly locked BEFORE parsing ClinVar or fetching gnomAD v4.1 for any of the 20 unseen validation genes)
- **Purpose**: Provide a genuine, never-inspected out-of-sample replication of Study A-1's post-audit findings outside the 84-gene ACMG SF v3.3 discovery panel.
- **Unseen 20-Gene ClinGen Definitive Autosomal Dominant Validation Panel (0/20 overlap with ACMG SF v3.3)**:
  - *Neurogenetics (6)*: `SCN1A`, `KCNQ2`, `STXBP1`, `SLC2A1`, `NOTCH3`, `PSEN1`
  - *Non-ACMG Hereditary Cancer/Tumor (8)*: `CDH1`, `CDKN2A`, `CHEK2`, `ATM`, `BAP1`, `DICER1`, `NF1`, `PTCH1`
  - *Renal / Connective / Metabolic / Hematologic (6)*: `PKD1`, `PKD2`, `COL1A1`, `ANK1`, `GCK`, `ELN`

## Pre-Specified Falsifiable Predictions (Locked Before Data Access)
1. **PV-1 (Multiallelic `rsID` Collision Mirage)**: Across all 20 unseen genes, naive `rsID`-join false-positive collision rate for `P/LP` (`FAF95 >= 1e-4`) will exceed `50%`.
2. **PV-2 (Denominator-Clean VUS Disparity & Submitter-Velocity Mediation)**: Restricting strictly to high-callability (`AN >= 80%`) coding/splice variants at `FAF95 >= 1e-4` (zero singletons) and Exact Hypergeometric Downsampling:
   - Non-bottleneck Non-NFE alleles will show `VUS CMH OR > 1.5` ($p < 0.01$) and `Conflicting CMH OR < 0.7` ($p < 0.01$);
   - Non-NFE alleles will have significantly fewer ClinVar submitters (Mann-Whitney $p < 10^{-10}$);
   - Controlling for `log(1 + n_submitters)` will attenuate $\ge 50\%$ of both the VUS and Conflicting log-ORs;
   - Replicating our post-audit falsification, aggregate `P/LP` odds ratio at `FAF95 >= 1e-4` will NOT exhibit large spurious inflation (unlike raw `AF >= 1e-4` with `MID` singletons).
3. **PV-3 (NFE-Invisibility of Actionable VUS-to-Benign Liberation Debt)**: Among high-callability (`AN >= 80%`) VUS exceeding `FAF95_grpmax5 >= 5e-4`, $\ge 75\%$ will be invisible in European `NFE` reference data (`FAF95_NFE < 5e-4`).
