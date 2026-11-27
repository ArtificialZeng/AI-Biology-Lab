# Study B-2 Replication & Verification Package

**Manuscript Title:** *Selective Enforcement, Organizational Governance, and Clinical-Trial Transparency: Quasi-Experimental Evidence and Identification Boundaries from the March 30, 2026 FDA Delinquency Sweep*

**Directory:** `/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep/`

---

## 1. Environments & Primary Data Provenance

- **Python:** `/opt/anaconda3/bin/python3` (`numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `pyarrow`).
- **Official CRAN R `HonestDiD` v0.2.8 (Isolated Conda Env):** `/opt/anaconda3/envs/r_honestdid/bin/Rscript` (`R 4.5.3`, `HonestDiD 0.2.8`, `doParallel`, `foreach`, `jsonlite`). Verified on Rambachan & Roth's official `BCdata_EventStudy` vignette (`audit/verify_honestdid_install.R`).
- **LaTeX:** TeX Live (`/Library/TeX/texbin/pdflatex`, `bibtex`, `/opt/homebrew/bin/gs`, `pdfinfo`).
- **Primary Public Sources (Archived Locally):**
  1. **ClinicalTrials.gov REST API v2 Snapshot (`2026-11-25`):** `data/ctgov_b2_master_cohort.parquet` (`53,736` post-Final-Rule trials; `16,333` strict-`ACTUAL` `COMPLETED|TERMINATED` U.S. FDA-regulated non-Phase-1 ACTs) and `audit/b2_t00_unposted_and_cert_overlay.parquet` (`22,532` unique overlay NCTs covering `8,936` unposted `estimatedResultsFirstSubmitDate`, `3,663` delayed certifications `dispFirstSubmitDate`, and `11,201` posted U.S. FDA-regulated `mcp` records).
  2. **FDA Primary Enforcement Record:** FDA GovDelivery Bulletin (sent `2026-04-13`, documenting compliance reminder messages issued `2026-03-30` to `>2,200` sponsors / `>3,000` trials): `https://content.govdelivery.com/accounts/USFDA/bulletins/412bbc8`.
  3. **NLM QC-Posting Policy:** NLM Technical Bulletin `2020-01-29`: `https://www.nlm.nih.gov/pubs/techbull/jf20/jf20_clinicaltrials_qc.html`.

---

## 2. Protocol & Forensic Amendment Trail

- **Post-Pilot Analysis Plan (Relabeled transparently):** `PREREG_B2.md` (original exploratory plan timestamp `2026-11-25 08:38 CST` preserved verbatim; documents T00–T10 forensic revisions and retractions, including retraction of raw-April Phase-1 spillover and broad efficacy "lemons falsification", plus explicit post-hoc labeling of `±5 pp` / `±10 pp` TOST margins).

---

## 3. End-to-End Runnable Replication Commands

Run from `/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep/`:

```bash
PY=/opt/anaconda3/bin/python3
RS=/opt/anaconda3/envs/r_honestdid/bin/Rscript

# Step 1: T00 API schema & unposted/certification overlay verification
$PY audit/audit_b2_t00_nlm_qc_api.py

# Step 2: T00-T04 Definitive estimation (Combined Tracking Proxy, T01 date-verified zero-censoring H=90/120/150d,
#         Clean-Exempt DDD, Placebos 2024/2025, H4 within-sponsor FE + MDE80 + post-hoc exploratory TOST + QC selection)
$PY audit/audit_b2_definitive_t00_to_t04.py

# Step 3: Quarterly event studies, heterogeneity decompositions (portfolio, cohort, backlog), Phase-1 spillover, 1,000-perm RI
$PY audit/audit_b2_quarterly_honestdid_and_mechanisms.py

# Step 4: T00 proxy-date forensic audit (0 fallback, 0 cert mixing), delayed-cert exclusion sensitivity (n=14,114),
#         T07 fixed pre-window (2025-04..2026-02, df=10) invariance proof, and 15-spec HonestDiD payload generation
$PY audit/audit_b2_t00_t06_t07_closure.py

# Step 5: T06/T10 Official CRAN R HonestDiD v0.2.8 (BCdata benchmark + all 15 payload specs + expanded +-40*SE grid fix on all 5 edge-hitting rows)
$RS audit/verify_honestdid_install.R
$RS audit/run_official_honestdid.R
$RS audit/append_missing_2_honestdid.R
$RS audit/fix_t10_honestdid_grid_boundaries.R

# Step 6: Regenerate all 4 publication-grade 300-DPI figures
$PY src/make_b2_figures_v2.py

# Step 7: Compile zero-overflow, >=9pt-table-font LaTeX manuscript PDF (requires TeX Live)
cd manuscript
pdflatex -interaction=nonstopmode manuscript_b2.tex
bibtex manuscript_b2
pdflatex -interaction=nonstopmode manuscript_b2.tex
pdflatex -interaction=nonstopmode manuscript_b2.tex
```

---

## 4. Key Output Files & Supplementary Tables

| File | Description |
| :--- | :--- |
| `manuscript/manuscript_b2.pdf` | Compiled 16-page manuscript (`0` undefined refs/cites, `0` overfull hboxes/vboxes, `0` `\resizebox`; `pypdf` min effective body font `8.9664 pt` across all pages; SHA-256 `8ab4219b0171d2e4522c8d09904156cd3dc1d310d74c4482f73b59255ee70b47`) |
| `manuscript/manuscript_b2.tex` & `references_b2.bib` | Complete LaTeX source & verified 22-entry BibTeX |
| `MANUSCRIPT_B2_DRAFT.md` | Synchronized Markdown manuscript |
| `audit/b2_supp_monthly_combined_proxy_hazards.csv` | **Supplementary Table S1:** Monthly combined-proxy submission hazards across all 4 treatment cells (`19` months, `2025-04..2026-10`) |
| `audit/b2_supp_t00_t07_master_estimates.csv` | **Supplementary Table S2:** Master table of `61` econometric estimates across T00 proxy, cert-exclusion, T01 zero-censoring, Clean-Exempt DDD, Phase-1 spillover, placebos, and H4 FE/MDE/exploratory TOST |
| `audit/b2_supp_t06_official_r_honestdid_all15.csv` | **Supplementary Table S3:** Complete `180`-row official CRAN R `HonestDiD` v0.2.8 output across all `15/15` specs with `grid_lb`, `grid_ub`, `gridPoints`, `grid_step`, `hits_grid_boundary` (`0` hits after T10 fix), and `expanded_grid_recomputed` flags |
| `audit/b2_t10_grid_boundary_and_convergence_audit.json` | T10 audit JSON documenting all 5 expanded `±40*SE` (`1,001`-pt) recomputations, dual-boundary rejection, and `1,001`-vs.-`2,001` convergence check |
| `audit/b2_t00_proxy_date_forensic_audit.json` | T00 forensic proof (`100.0%` `est == releaseDate`, `0` fallbacks, `0` cert mixing across `53,736` trials) |
| `audit/b2_t07_invariance_and_es_clean.json` | T07 calendar-locked pre-window invariance proof (`max diff <= 4.9e-12`) & delayed-cert exclusion sensitivity |
| `audit/b2_t06_official_r_honestdid_results.json` | Authoritative JSON of `BCdata` vignette benchmark + all `15/15` official CRAN R `HonestDiD` v0.2.8 specs (post-T10 fix) |
| `audit/b2_definitive_t00_to_t04_results.json` | Authoritative JSON of T00–T04 DiD/DDD, T01 zero-censoring horizons, and H4 FE/TOST/selection checks |
| `audit/b2_t00_quarterly_honestdid_and_mechanisms.json` | Quarterly ES, heterogeneity decompositions, Phase-1 spillover, and 1,000-perm RI |
| `figures/fig{1..4}_b2_*.png` | All 4 300-DPI figures |
