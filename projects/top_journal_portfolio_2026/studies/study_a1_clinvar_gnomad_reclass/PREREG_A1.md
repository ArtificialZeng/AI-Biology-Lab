# Post-Pilot Analysis Plan & Audit Amendment Log — Study A-1 (formerly `PREREG_A1`)

> **PROTOCOL PROVENANCE DISCLOSURE (Added 2026-11-25 10:05 +08:00)**:
> This document was originally drafted at **`2026-11-25 08:45 (+08:00)`** *after* an 8-gene exploratory pilot (`BRCA1/2`, `MYBPC3`, `MYH7`, `LDLR`, `TTR`, `SCN5A`, `KCNQ1`) had already examined ClinVar × gnomAD v4.1 rsID inflation and ancestry-stratified frequencies, and *before* full 84-gene ACMG SF v3.3 extraction and hypothesis testing. Accordingly, it is **NOT a pre-pilot preregistration**; it is accurately classified as a **Post-Pilot Analysis Plan** with original timestamp preserved verbatim below.
> By contrast, [`PROSPECTIVE_VALIDATION_A1.md`](file:///Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/PROSPECTIVE_VALIDATION_A1.md) was **genuinely prospectively locked at `2026-11-25 09:51 (+08:00)`** *prior* to querying any ClinVar or gnomAD v4.1 data for the 20 unseen non-ACMG autosomal dominant validation genes (`99,249` ClinVar × `157,021` gnomAD variants).
>
> **Substantive Post-Audit Revisions & Retractions Logged Transparently (`2026-11-25 09:42–10:05 +08:00`)**:
> 1. **Retraction of H1c Aggregate `P/LP` Disparity (Singleton Denominator Bias)**: Audit (`audit/a1_forensic_audit.json`) revealed that raw `AF >= 1e-4` induced `64.7%` `AC == 1` singletons in small-AN `MID` (`median AN = 6,006`). Under denominator-clean `FAF95 >= 1e-4`, `AC >= 10 & AF >= 5e-4`, and exact hypergeometric downsampling (`AN_NFE == AN_p`), the aggregate `P/LP` disparity collapses to null/reverses (`CMH OR = 1.49, p = 0.434`; replicated null `OR = 0.36, p = 0.169` in 20 unseen genes). H1c is retracted and reported as a methodological warning; **H1a (VUS persistence) and Conflicting reversal survive every denominator-clean test**.
> 2. **Complete 55/55 Live gnomAD QC & Gene-Specific ClinGen VCEP BS1/BA1 Re-Audit (C06)**: Live gnomAD v4.1 QC filters caught stealth `AS_VQSR` failure in `APC 5-112775626-T-A` despite `93.9%` AN (`8/55` QC/callability artifacts total). Applying published ClinGen VCEP/Whiffin BS1 rules (`7/55` below BS1) and ClinVar submission-level mode-of-inheritance parsing (`4/55` AR myopathy alleles in dual-mode `RYR1`/`TTN`) leaves **`0/55` clear frequency-only P/LP misclassifications**! Similarly, uniform VUS `FAF95` "immediately benign" claims were retracted and replaced with **gene-specific ClinGen VCEP (29 genes: `122` BS1 / `16` BA1) and VCEP+Whiffin (84 genes: `313` BS1 / `44` BA1, excluding Ghosh 2018 BA1 exceptions) frequency-conflict VUS candidates prioritized for formal expert review**, plus a 20-variant live QC spot-check (`90%` pass).
> 3. **Non-Causal Submitter Attenuation & Collapsibility Fix (C07)**: Conditioning on `n_submitters` is strictly downgraded from "causal mediation / curation velocity" to **descriptive cross-sectional covariate attenuation** (citing Imai et al. 2010 and Mood 2010), evaluated via collapsible LPM risk differences, Logistic AMEs, non-parametric submitter-tier stratification (`1, 2–3, 4–9, 10–19, 20+`), and `LastEvaluated` recency.
> 4. **Citation & Prior-Art Corrections (C05)**: `Chen et al. Nature 625:92–100 (2024)` cited accurately for gnomAD v3.1.2 alongside Broad gnomAD v4.0/v4.1 release docs; explicit side-by-side differentiation added vs. Gudmundsson et al. (*Nat Commun* 2025), Dario (*Ann Hum Genet* 2026), Lai et al. (*Front Genet* 2026), Hotakainen et al. (*Genet Med* 2025), and Ernest & Ansah (*Research Square* 2026).

---

*Original Post-Pilot Protocol Text Preserved Verbatim Below (`2026-11-25 08:45 +08:00`):*

# Pre-Registered Analysis Protocol — Study A-1 (PREREG_A1)

- **Study ID**: `A-1` (Track A Flagship — Target Journals: *Nature Genetics*, *Nature Medicine*, *American Journal of Human Genetics*)
- **Title**: *Allele-Exact Cross-Ancestry Reclassification Debt in ClinVar Across All 84 ACMG SF v3.3 Medically Actionable Genes Under gnomAD v4.1 (807,162 Individuals)*
- **Protocol Lock Timestamp**: `2026-11-25 08:45 (+08:00)` (locked prior to allele-exact 84-gene ClinVar × gnomAD v4.1 join and hypothesis testing)
- **Data Sources & Versions**:
  - **ClinVar GRCh38 Release**: NCBI FTP `clinvar.vcf.gz` (`Last-Modified: 2026-11-24`, `SHA-256: 5aaf944cad89c8b9d409ad74e4eae8b414c14c7f10fd1247fd677a6f445d68d5`) and `variant_summary.txt.gz` (`Last-Modified: 2026-11-23`, `SHA-256: 1d910abfb5b576d92abcf042e118dd7ae1c6cf9e91382931093ae5555816206a`). Terms: `https://www.ncbi.nlm.nih.gov/clinvar/`.
  - **gnomAD v4.1 (`gnomad_r4`)**: Broad Institute gnomAD v4.1 (`730,947` exomes + `76,215` genomes = `807,162` individuals; 8 genetic ancestry groups: `nfe`, `afr`, `amr`, `asj`, `eas`, `fin`, `mid`, `sas`), queried exclusively via the official `gnomad-database` skill wrapper (`search_variants.py` / `get_variant_frequency.py`). Terms: `https://gnomad.broadinstitute.org/policies`.
  - **Gene Panel**: Complete **84-gene ACMG SF v3.3** medically actionable secondary findings list (Lee et al., *Genetics in Medicine* 2025; 27(8):101454; 72 Autosomal Dominant, 9 Autosomal Recessive, 3 X-linked).

---

## 1. Scientific Problem & Pre-Specified Hypotheses

Clinical sequencing laboratories rely on ClinVar classifications across the 84 ACMG SF v3.3 genes to trigger life-altering interventions (prophylactic mastectomy/oophorectomy, implantable cardioverter-defibrillators, cascade family testing). Earlier cross-ancestry variant misclassification audits relied on ExAC (`60k`) or gnomAD v2/v3 (`141k`–`76k`), often suffered from **multiallelic rsID join artifacts**, and predated both the **5.7× expansion of non-European individuals in gnomAD v4.1 (`807k`)** and the expansion of ACMG actionable genes to **84 genes (SF v3.3)**.

**Pre-Specified Primary Hypotheses**:
1. **H1 — Frequency-Matched Cross-Ancestry Classification Asymmetry (Allele-Exact)**:
   Among allele-exact `(CHROM, POS, REF, ALT)` matched coding and canonical-splice variants across all 84 ACMG SF v3.3 genes with max-ancestry filtering allele frequency (`FAF95_popmax >= 1e-4` or `AF_popmax >= 1e-4`), variants enriched in non-European genetic ancestry groups (`AF_nonNFE >= 1e-4` and `AF_nonNFE / max(AF_NFE, 1e-6) >= 5x`) will exhibit:
   - **(H1a)** Significantly **higher odds of being classified as VUS (`Uncertain_significance`) or Conflicting** vs. frequency-matched European-enriched (`NFE`) variants (conditional on exact matched allele frequency bin, gene, and molecular consequence — Mantel-Haenszel & conditional logistic regression!);
   - **(H1b)** Significantly **lower odds of being resolved as Benign/Likely Benign (`B/LB`)**;
   - **(H1c)** Significantly **higher odds of remaining classified as Pathogenic/Likely Pathogenic (`P/LP`)** despite exceeding modern population-frequency ceilings.
2. **H2 — High-Confidence (`>= 2-Star`) ClinVar "Reclassification Debt" Under Exact Whiffin/ClinGen `FAF95` Ceilings**:
   Quantify the exact number, clinical domain concentration (Cardiovascular vs. Cancer vs. Metabolic/Misc), inheritance mode (`AD` vs. `AR` vs. `XL`), and review staleness (`LastEvaluated`) of **allele-exact `2–4 star` ClinVar `P/LP` variants** in the 72 Autosomal Dominant ACMG SF v3.3 genes whose exact gnomAD v4.1 joint (`exome + genome`) **Grpmax 95% Filtering Allele Frequency (`FAF95_max`)** exceeds:
   - Tier 1 (Strict AD BS1 threshold): `FAF95_max >= 1.0e-4` (1 in 5,000 carriers — exceeds maximum credible allele frequency for most high-penetrance AD cardiomyopathies, channelopathies, and cancer syndromes under Whiffin et al. *Genet Med* 2017);
   - Tier 2 (High-Alert AD BS1 threshold): `FAF95_max >= 5.0e-4` (1 in 1,000 carriers);
   - Tier 3 (Stand-alone BA1/Extreme threshold): `FAF95_max >= 1.0e-3` or `AF_max >= 5.0e-3`, OR presence of `>= 2` gnomAD v4.1 homozygotes in a severe AD gene.
   - Measure **Ancestry Attribution**: What fraction of these `FAF95`-flagged `>=2-star` AD `P/LP` variants would have been **invisible in European (`NFE`) reference panels** (`FAF95_NFE < 1e-4` while `FAF95_nonNFE >= 1e-4`)?
3. **H3 — Review-Staleness & Expert-Panel Governance Lag**:
   Test whether `FAF95`-flagged `P/LP` variants (`>=2-star` and `1-star`) have systematically **older `LastEvaluated` dates** (pre-gnomAD v4.0 release, Nov 2023) vs. non-flagged `P/LP` variants — AND test whether a substantial fraction were **actively re-evaluated in 2024–2026 AFTER gnomAD v4.0** yet still retain unconflicted `2–3 star` `P/LP` status!
4. **H4 — Methodological Audit of Multiallelic `rsID` Join Inflation**:
   Compare naive `rsID`-based ClinVar×gnomAD joins against exact GRCh38 `(CHROM, POS, REF, ALT)` joins across all 84 ACMG SF v3.3 genes to quantify the exact false-positive inflation rate caused by multiallelic `rsID` collisions (which our 8-gene pilot flagged at `61.1%`!).

---

## 2. Exact Statistical & Bioinformatics Methods (Pre-Specified)

- **Exact Variant Key**: Normalized GRCh38 `chrom-pos-ref-alt` (matching gnomAD v4.1 VCF minimal representation directly to ClinVar `clinvar.vcf.gz` `CHROM-POS-REF-ALT` with exact left-aligned/trimmed check). Zero `rsID` joins in primary analysis.
- **Exact Joint Exome+Genome Counts & `FAF95`**:
  For each ancestry group $p \in \{\text{nfe}, \text{afr}, \text{amr}, \text{asj}, \text{eas}, \text{fin}, \text{mid}, \text{sas}\}$:
  - $AC_p = AC_{p,\text{exome}} + AC_{p,\text{genome}}$, $AN_p = AN_{p,\text{exome}} + AN_{p,\text{genome}}$, $AF_p = AC_p / AN_p$ (requiring $AN_p \ge 200$ for reliable popmax).
  - Exact **Filtering Allele Frequency (`FAF95`)** per Whiffin et al. (*Genetics in Medicine* 2017): for $AC_p \ge 2$, $\text{FAF95}_p = \frac{1}{2} \chi^2_{2(AC_p - 1),\, 0.05} / AN_p$ (the exact 95% lower confidence bound on the Poisson rate after subtracting 1 allele for single-case allowance; $0$ if $AC_p \le 1$), cross-validated against gnomAD's native `joint.faf95` API output on a random sample of top hits via `get_variant_frequency.py`!
- **Frequency-Matched Stratification (Eliminating Frequency Confounding)**:
  - Bin variants by exact $\log_{10}(\text{AF}_{\max})$ deciles / half-log bins (`[1e-4, 3e-4)`, `[3e-4, 1e-3)`, `[1e-3, 3e-3)`, `[3e-3, 1e-2)`, `[1e-2, 5e-2)`, `[5e-2, 1.0]`).
  - Estimate **Cochran-Mantel-Haenszel (CMH) stratified odds ratios** and **multivariable logistic regression with Gene FE + Consequence FE + $\log_{10}(\text{AF}_{\max})$ spline** (plus separate sensitivity excluding `asj` & `fin` founder populations to separate non-European under-representation from European founder bottlenecks!).
- **Individual-Level Clinical Exposure Estimation**:
  For each ancestry group $p$, compute the **cumulative population carrier burden** of `FAF95`-flagged `>=2-star` AD `P/LP` variants: $C_p = 1 - \prod_{v} (1 - 2 \cdot AF_{v,p} (1 - AF_{v,p}) - AF_{v,p}^2)$ (exact Hardy-Weinberg probability that a sequenced individual of ancestry $p$ carries $\ge 1$ flagged high-confidence ClinVar `P/LP` allele across the 72 AD ACMG SF v3.3 genes!).

---

## 3. Verification & Reporting Commitments
- All 84 genes, exact ClinVar Variation IDs, exact GRCh38 `<chrom>:<pos>:<ref>><alt>` coordinates, exact `FAF95` values, and independent gnomAD `get_variant_frequency.py` spot-checks of top hits will be published transparently.
- Any founder-effect pathogenic variants (e.g., `TTR p.Val142Ile` in `AFR`, `HFE p.Cys282Tyr` in `NFE`, `BRCA1/2` Ashkenazi founder alleles, `MYBPC3` South Asian 25-bp deletion, `KCNQ1` Finnish founder) will be explicitly separated and annotated as **true high-frequency founders / reduced-penetrance alleles vs. likely benign misclassification debt**, ensuring gold-standard clinical-genetics rigor!
