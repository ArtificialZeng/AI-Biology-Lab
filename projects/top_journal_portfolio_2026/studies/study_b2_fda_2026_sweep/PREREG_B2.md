# Post-Pilot Analysis Plan & Audit Amendment Log — Study B-2 (formerly `PREREG_B2`)

> **PROTOCOL PROVENANCE DISCLOSURE (Added 2026-11-25 10:05 +08:00)**:
> This document was originally drafted at **`2026-11-25 08:38 (+08:00)`** *after* an exploratory feasibility pilot had already examined aggregate 2026 ClinicalTrials.gov monthly submission counts, and *before* full trial-level universe extraction (`53,927` unique NCTs) and econometric estimation. Accordingly, it is **NOT a pre-pilot preregistration**; it is accurately classified as a **Post-Pilot Analysis Plan** with original timestamp preserved verbatim below.
>
> **Substantive Post-Audit Revisions Logged Transparently (`2026-11-25 09:42–10:15 +08:00`)**:
> 1. **T00 (NLM 2020-01-29 QC-Posting Policy Correction)**: The original Section 2.3 text below erroneously stated that ClinicalTrials.gov hides submission dates until QC completes. Live API v2 audit (`audit/b2_t00_nlm_qc_api_audit.json`) confirmed that while `protocolSection.statusModule.resultsFirstSubmitDate` is indeed populated only after QC clearance (`hasResults=true`, 0/500k+ exceptions), unposted/in-QC and 30-day QC-comment (`MCP`) submissions ARE exposed in `derivedSection.miscInfoModule.submissionTracking.estimatedResultsFirstSubmitDate`. All primary models were therefore re-estimated using **uncensored true first submissions (`resultsFirstSubmitDate` $\cup$ `estimatedResultsFirstSubmitDate`)** alongside separate QC-cleared posting and delayed-certification (`dispFirstSubmitDate`) analyses (`audit/b2_definitive_t00_to_t04_results.json`).
> 2. **T01 (Strict Calendar-Date Observability Assertions)**: Uniform-horizon QC-clearance cutoffs were corrected so that `end_date + H <= 2026-11-25` holds with zero-day violation, and explicitly reframed as composite submission+QC-clearance outcomes.
> 3. **T02 / C02 (Pre-Trend Failure & Rambachan–Roth HonestDiD)**: Monthly joint pre-trend Wald rejections are reported transparently alongside exact Rambachan & Roth (2023) $\Delta^{RM}(\bar{M})$ and $\Delta^{SD}(M)$ partial-identification bounds; any informal claim that a negative pre-trend implies a causal lower bound was retracted.
> 4. **T03 / T04 / C03 (H4 Lemons Power & Honest Downgrade)**: Wide within-sponsor CIs (e.g., `[-17.54, +14.86] pp` for primary $p \ge 0.05$) are reported with exact 80%-power MDEs, TOST bounds, and unposted-vs-posted termination selection checks; "falsification of lemons" claims were retracted.
> 5. **Symmetric Sample & Unseen Validation**: Primary samples restricted strictly to symmetric `COMPLETED|TERMINATED` + `ACTUAL` PCD; new risk-set Overdue-vs-On-Time DDD and unseen `2012–2017` pre-Final-Rule legacy cohort (`9,394` trials) added transparently.

---

*Original Post-Pilot Protocol Text Preserved Verbatim Below (`2026-11-25 08:38 +08:00`):*

# Pre-Registered Analysis Protocol — Study B-2 (PREREG_B2)

- **Study ID**: `B-2` (Track B Flagship — Target Journals: *Management Science*, *Journal of Financial Economics*, *Research Policy*)
- **Title**: *Selective Enforcement, Organizational Governance, and Clinical-Trial Transparency: Quasi-Experimental Evidence from the March 30, 2026 FDA Mass Delinquency Sweep*
- **Protocol Lock Timestamp**: `2026-11-25 08:38 (+08:00)` (locked prior to full-universe trial-level data extraction and estimation)
- **Data Source & Terms**: ClinicalTrials.gov REST API v2 (`https://clinicaltrials.gov/`, accessed exclusively via `clinical-trials-database` skill wrapper `clinical_trials_api.py`).

---

## 1. Institutional Shock & Research Question

On **March 30, 2026**, U.S. Food and Drug Administration (FDA) Commissioner Dr. Marty Makary launched an unprecedented agency-wide enforcement outreach campaign sending formal delinquency reminders to **more than 2,200** sponsors covering **more than 3,000** overdue clinical trials subject to mandatory results reporting under **FDAAA Section 801** (42 U.S.C. 282(j)) and its Final Rule (**42 CFR Part 11**, effective January 18, 2017). The outreach explicitly cited an internal FDA audit showing **29.6%** non-compliance among likely applicable trials and warned of escalating Pre-Notices of Noncompliance, formal Notices of Noncompliance, and inflation-adjusted Civil Money Penalties (CMPs) up to **$15,107 per day** for uncorrected violations beyond 30 days.

**Core Questions**:
1. **Primary Causal Effect (H1 — Enforcement Efficacy)**: Did the March 30, 2026 mass reminder sweep increase the weekly/monthly hazard of overdue results submission among FDA-regulated Applicable Clinical Trials (ACTs) relative to statutorily and jurisdictionally exempt clinical trials?
2. **Primary Governance Asymmetry (H2 — Corporate vs. Academic Responsiveness)**: Did corporate/industry sponsors (`leadSponsor.class == INDUSTRY`) respond significantly more strongly to the FDA sweep than academic/nonprofit/government sponsors (`OTHER`, `NIH`, `FED`, `NETWORK`, `OTHER_GOV`) in a Triple-Difference (DDD) design?
3. **Mechanism Tests (H3 — Organizational Capacity & Repeat-Player Liability)**: Within industry and academic sponsors, is enforcement responsiveness concentrated among (a) multi-trial repeat-player sponsors with larger delinquent backlogs and ongoing FDA regulatory pipelines, (b) more recently overdue vs. deeply stale cohorts, and (c) drug/biologic vs. device trials?
4. **Adverse Selection / "Market for Lemons" Disclosure Test (H4)**: Conditional on sponsor, phase, and therapeutic area, do delinquent trials forced into disclosure *after* March 30, 2026 exhibit higher early-termination rates, higher primary-outcome statistical non-significance ($p \ge 0.05$) or p-value omission, and/or higher Serious Adverse Event (SAE) and mortality rates than pre-sweep or on-time disclosed trials?

---

## 2. Exact Universe & Risk-Set Definition (Pre-Specified)

### 2.1 Primary Cohort (`2017–2025` Statutory Post-Final-Rule Universe)
All `INTERVENTIONAL` trials on ClinicalTrials.gov with:
- `OverallStatus` in `{COMPLETED, TERMINATED}`
- `PrimaryCompletionDate` in `[2017-01-18, 2025-02-28]` (post-42 CFR Part 11 Final Rule effective date, and $\ge 13$ months before March 30, 2026 so every unsubmitted trial is strictly past the statutory 1-year submission deadline $d_i = c_i + 365\text{ days}$ on the event date).
- For partial `YYYY-MM` dates, day defaults conservatively to the last day of the month (biasing against falsely classifying on-time trials as delinquent).

### 2.2 Treatment vs. Exogenous Exempt Control Arms
1. **Treated Arm (`ACT_FDA = 1`)**: Trials with `has_us_site == True` AND (`isFdaRegulatedDrug == True` OR `isFdaRegulatedDevice == True`) AND NOT strictly Phase-1-only (`phases != ['PHASE1']` and `phases != ['EARLY_PHASE1']`).
   - *Alternative statutory proxy (robustness)*: US Phase 2/3 drug/biologic/device interventional trials (`us_p23_act`).
2. **Control Arm 1 (`CTRL_US_NONFDA`, Primary Domestic Control)**: `has_us_site == True` AND `isFdaRegulatedDrug == False` AND `isFdaRegulatedDevice == False` (e.g., behavioral, surgical, diagnostic, dietary US trials sharing US institutional calendars and NLM QC pipeline, but statutorily exempt from FDA FDAAA 801 CMPs).
3. **Control Arm 2 (`CTRL_US_P1`, Statutory Phase-1 Control)**: `has_us_site == True` AND strictly `PHASE1` / `EARLY_PHASE1` (statutorily exempt from mandatory results submission under 42 CFR 11.22(b)).
4. **Control Arm 3 (`CTRL_NONUS_P23`, Jurisdictional Foreign Control)**: `has_us_site == False` AND `Phase2` or `Phase3` AND `isFdaRegulatedDrug == False` AND `isFdaRegulatedDevice == False` (outside U.S. FDA jurisdiction).

### 2.3 Right-Censoring Cutoff for NLM Quality Control (QC) Lag
ClinicalTrials.gov only exposes `resultsFirstSubmitDate` publicly *after* a submission clears NLM QC review (~30–60 days). Because the data pull occurs on `2026-11-25`, submissions made in late October–November 2026 are severely right-censored.
- **Primary Analysis Window**: **18 calendar months** (`2025-04` through `2026-09`) — 12 pre-sweep months (`2025-04..2026-03`, with March 2026 split or treated conservatively as pre-sweep since the sweep occurred on **March 30, 2026**; weekly models split precisely at the Monday `2026-03-30` boundary!) and **6 post-sweep months** (`2026-04..2026-09`, all $\ge 55$ days before data pull).
- **Sensitivity Cutoffs**: `2026-08-31` (86-day QC buffer) and `2026-10-15` (41-day buffer).

---

## 3. Pre-Specified Econometric Models

### 3.1 Exact Delinquent Trial-Month & Trial-Week Risk-Set Panel
For every trial $i$ and time period $t$ (month or week):
- Trial $i$ enters the delinquent risk set ($R_{it} = 1$) iff $t > c_i + 365\text{ days}$ AND trial $i$ has not submitted results strictly before period $t$.
- Outcome $Y_{it} = 1$ iff trial $i$ submits first results (`resultsFirstSubmitDate`) during period $t$; $0$ otherwise. Trial $i$ exits the risk set after period $t$ once $Y_{it} = 1$.

### 3.2 Primary Specifications
1. **Difference-in-Differences (DiD — H1)** on delinquent risk-set panel (`ACT_FDA` vs. `CTRL_US_NONFDA`, and separately vs. `CTRL_US_P1` and `CTRL_NONUS_P23`):
   $$Y_{it} = \beta \, (\text{ACT}_i \times \text{Post}_t) + \gamma_i \text{ or } \alpha_{g(i), \text{cohort}(i)} + \lambda_t + \varepsilon_{it}$$
   Estimated via: (a) Linear Probability Model (OLS) with two-way fixed effects (`Arm × Completion-Year × Phase` + `Period` FE, and `Sponsor FE + Period FE`), (b) Poisson / Log-linear hazard model on cell-level submission counts with $\log(\text{RiskSet}_{gt})$ offset, and (c) Cox Proportional Hazards model with time-varying $\text{ACT}_i \times \text{Post}_t$. Standard errors clustered at the **Lead Sponsor** level.
2. **Triple-Difference (DDD — H2, Main Governance Specification)**:
   $$Y_{it} = \beta_{\text{DDD}} \, (\text{ACT}_i \times \text{Industry}_i \times \text{Post}_t) + \beta_{\text{ACT}\times\text{Post}} \, (\text{ACT}_i \times \text{Post}_t) + \beta_{\text{Ind}\times\text{Post}} \, (\text{Industry}_i \times \text{Post}_t) + \mu_{\text{Arm}\times\text{Ind}\times\text{Cohort}} + \lambda_t + \varepsilon_{it}$$
   Testing $H_0: \beta_{\text{DDD}} = 0$ (two-sided, $\alpha = 0.05$) with sponsor-clustered SEs AND exact event-study monthly/biweekly leads & lags ($k = -12, \dots, -2, 0, +1, \dots, +5$, reference $k = -1$, i.e., `2026-02` / pre-sweep week `-1`) with joint Wald pre-trend F-test on $k \in \{-12, \dots, -2\}$.
3. **Placebo Sweep Falsification Tests**:
   Exact identical risk-set and DDD estimation run on two historical placebo dates when NO FDA mass sweep occurred:
   - **Placebo 2025**: Pseudo-sweep on `2025-03-30` (`2024-04` to `2025-09`).
   - **Placebo 2024**: Pseudo-sweep on `2024-03-30` (`2023-04` to `2024-09`).
   - Falsification pass criterion: Placebo $\hat{\beta}_{\text{DDD}}$ near zero / insignificant and strictly smaller than true 2026 $\hat{\beta}_{\text{DDD}}$.
4. **Lemons / Adverse-Selection Regressions (H4)**:
   Among all FDA-regulated ACT trials submitting results in `2024-01-01..2026-09-30` (both overdue and on-time comparator):
   - Compare (i) `Post-Sweep Overdue` (`2026-03-30..2026-09-30`, delay $>365$d) vs. (ii) `Pre-Sweep Overdue` (`2024-01-01..2026-03-29`, delay $>365$d) and (iii) `On-Time` (delay $\le 365$d), controlling for Sponsor Class (and Sponsor FE in within-sponsor subsample), Phase, Enrollment ($\log N$), and Therapeutic Area.
   - Outcomes: `Terminated` status, `NonSignificant_Primary` ($p_{\min} \ge 0.05$), `Omitted_PValue`, `SAE_Rate`, `AllCause_Mortality_Rate`, and `Submission_Delay_Days`.

---

## 4. Reporting & Transparency Commitments
- Every single test (including null or opposite-sign secondary/lemons results) will be reported verbatim without cherry-picking.
- All raw JSON caches, cleaned Parquet tables, econometric scripts, and verification checksums will be preserved in `studies/study_b2_fda_2026_sweep/`.
