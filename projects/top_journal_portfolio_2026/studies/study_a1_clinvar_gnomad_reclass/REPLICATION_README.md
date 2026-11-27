# Study A-1 Replication & Verification Package

**Manuscript Title:** *Allele-Exact, Ploidy- and Denominator-Controlled Cross-Ancestry Audit of ClinVar Across All 84 ACMG SF v3.3 Genes Under gnomAD v4.1: Pipeline Artifacts, Recessive Return Hazards, and Verified ClinGen CSpec / DAFT VUS Frequency Conflicts*

**Directory:** `/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/`

---

## 1. Environment & Primary Data Provenance

- **Python:** `/opt/anaconda3/bin/python3` (`numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `pyarrow`, `beautifulsoup4`).
- **LaTeX:** TeX Live (`/Library/TeX/texbin/pdflatex`, `bibtex`, `/opt/homebrew/bin/gs`, `pdfinfo`).
- **Primary Public Data Snapshots (Archived Locally):**
  1. **ClinVar GRCh38 Bulk FTP:** `clinvar.vcf.gz` (`2026-11-24`, `4,621,342` rows; `343,328` in ACMG SF v3.3) and `variant_summary.txt.gz` (`2026-11-23`). Public domain (`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/`).
  2. **gnomAD v4.1 (`gnomad_r4`, `807,162` individuals: `730,947` exomes [`367,323` XX / `363,624` XY] + `76,215` genomes [`38,942` XX / `37,273` XY]; joint `406,265` XX / `400,897` XY):** Broad Institute GraphQL API (`https://gnomad.broadinstitute.org/api`), cached per gene (`545,493` variants across 84 ACMG SF v3.3 genes; `157,021` variants across 20 unseen autosomal AD genes). Official v4.1 XX/XY counts verified from gnomAD team post (`2024-06-28`, `https://discuss.gnomad.broadinstitute.org/t/xx-and-xy-sample-numbers-for-v4-exomes-and-genomes/328`), yielding official non-PAR chrX `MAX_AN = N_XY + 2*N_XX = 1,213,427` (`80% = 970,741.6`) vs autosomal `2N = 1,614,324` (`80% = 1,291,459.2`).
  3. **ClinGen Criteria Specification Registry (CSpec) Live API:** All `250` SVI documents retrieved `2026-11-25` (`https://cspec.genome.network/cspec/api/SequenceVariantInterpretation/id/`), archived at `audit/cspec_live_acmg84_all_verified.json` (`30` officially `Released` ACMG SF v3.3 VCEP specs, incl. *APC* `GN089` v2.1 `BS1=1e-5`, *PMS2* `GN139` v2.0 `BS1=2.8e-4`, *ABCD1* `GN105` v1.0, and *OTC* `GN156` v1.0).
  4. **Johnston et al. (2026) Published DAFT Baseline:** *Genetics in Medicine* 28(5):102552, DOI `10.1016/j.gim.2026.102552`, `PMC13170090`, archived at `audit/johnston2026_PMC13170090_parsed.json` (`6` preliminary VCEP genes + `47` published DAFT genes).

---

## 2. Protocol & Audit Trail

- **Post-Pilot Analysis Plan (Relabeled transparently):** `PREREG_A1.md` (original exploratory plan timestamp `2026-11-25 08:45 CST` preserved verbatim; documents every forensic retraction and revision).
- **Genuinely Prospective Out-of-Sample Validation Lock:** `PROSPECTIVE_VALIDATION_A1.md` (locked `2026-11-25 09:51 CST` *prior* to extracting any gnomAD or ClinVar data for the 20 unseen autosomal AD genes).

---

## 3. End-to-End Runnable Replication Commands

Run from `/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/`:

```bash
PY=/opt/anaconda3/bin/python3

# Step 1: Core allele-exact ClinVar x gnomAD v4.1 join & baseline H1-H4 statistics (uses cached data/)
$PY src/analyze_a1_reclassification.py

# Step 2: Forensic denominator-bias audit (equal-AN hypergeometric downsampling, FAF95/AC>=10, live API FAF95 r=0.9997)
$PY audit/forensic_audit_a1.py

# Step 3: Prospective 20-gene unseen non-ACMG AD validation + complete 55/55 live gnomAD v4.1 QC audit
$PY audit/run_a1_prospective_unseen_20genes.py

# Step 4: Collapsible LPM risk differences, logistic AMEs, non-parametric submitter-tier stratification, LastEvaluated recency
$PY audit/audit_a1_vcep_bs1_and_submitters.py

# Step 5a: T05 Definitive Live ClinGen CSpec (30 Released VCEPs) + Johnston et al. (2026) DAFT (6 Prelim + 47 DAFT)
$PY audit/audit_a1_t05_cspec_and_johnston2026.py

# Step 5b: T08 Official chrX Non-PAR Ploidy Denominator Fix (1,213,427; 406/441 X-linked VUS; OTC >=10 hemizygote BA1; 18/18 live chrX QC)
$PY audit/audit_a1_t08_chrx_ploidy_fix.py

# Step 6: Regenerate all 4 publication-grade 300-DPI figures
$PY src/make_a1_figures_v2.py

# Step 7: Compile zero-overflow LaTeX manuscript PDF (requires TeX Live)
cd manuscript
pdflatex -interaction=nonstopmode manuscript_a1.tex
bibtex manuscript_a1
pdflatex -interaction=nonstopmode manuscript_a1.tex
pdflatex -interaction=nonstopmode manuscript_a1.tex
```

---

## 4. Key Output Files & Supplementary Tables

| File | Description |
| :--- | :--- |
| `manuscript/manuscript_a1.pdf` | Compiled 18-page manuscript (`0` undefined refs/cites, `0` overfull hboxes/vboxes; `pypdf` effective body font: Tables 1, 2, 4 native `7.97 pt`, Table 3 `6.997 pt`; SHA-256 `73e625710e011f09fb33c8a124ef9e7ed3d8e1ee3697615b436ec03610e7bceb`) |
| `manuscript/manuscript_a1.tex` & `references_a1.bib` | Complete LaTeX source & verified 27-entry BibTeX |
| `MANUSCRIPT_A1_DRAFT.md` | Synchronized Markdown manuscript |
| `audit/a1_acmg84_verified_cspec_johnston2026_rules.csv` | **Supplementary Table S1:** 100% primary-source-verified 84-gene BS1/BA1 catalog (CSpec ID, version, date, metric, exceptions, verbatim text) |
| `audit/a1_ad55_verified_t05_triage.csv` | **Supplementary Table S2:** Complete 55/55 allele-exact 2–4★ AD P/LP live-QC + CSpec/DAFT + MOI audit |
| `audit/a1_vus_verified_t05_t08_bs1_candidates.csv` | **Supplementary Table S3:** All `1,388` ploidy-aware high-callability VUS exceeding verified BS1/BA1 across 83 genes (`327` Tier 1 Released CSpec incl. `11` *ABCD1* & `5` *OTC* hemizygote BA1; `99` Tier 2 Prelim; `962` Tier 3 Johnston DAFT incl. `1` *GLA*) |
| `audit/a1_t08_chrx_ploidy_closure.json` | Authoritative T08 JSON (official `1,213,427` denominator, `406/441` chrX VUS recovery, ploidy-aware `n=10,058` H1, `18/18` live chrX QC) |
| `audit/a1_t05_verified_closure.json` | Machine-readable JSON of T05 verified triage, `AC>=5` sensitivities, and `20/20` autosomal live QC spot-check |
| `audit/a1_prospective_unseen20_and_qc55.json` | Prospective 20-gene unseen AD validation + `55/55` live gnomAD QC JSON |
| `audit/a1_vcep_bs1_and_collapsible_submitters.json` | Collapsible LPM, AME, non-parametric submitter-tier stratification, and `LastEvaluated` audit JSON |
| `figures/fig{1..4}_a1_*.png` | All 4 300-DPI figures |
