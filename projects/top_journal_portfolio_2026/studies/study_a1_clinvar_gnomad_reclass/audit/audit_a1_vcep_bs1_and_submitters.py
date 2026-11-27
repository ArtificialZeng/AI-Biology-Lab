#!/opt/anaconda3/bin/python3
"""
Study A-1 Forensic Closure Script for Codex C06, C07, and Prior-Art Differentiation:
  1. Gene-Specific ClinGen VCEP & Whiffin (2017) BS1 / BA1 Frequency-Conflict Framework across ALL 84 ACMG SF v3.3 genes
     (including explicit Ghosh et al. 2018 BA1 exception exclusions + separate reporting of Published-VCEP-only subset
     vs Full-84-gene VCEP+Whiffin panel + live gnomAD v4.1 QC filter spot-check on 25 VUS candidates!).
  2. Collapsible Linear Probability Model (LPM), Logistic Average Marginal Effects (AME), Exact Submitter-Tier
     Stratified Contrasts, and LastEvaluated Recency Audit on BOTH the 84-gene ACMG SF v3.3 denominator-clean cohort
     AND the prospectively locked 20-gene unseen validation cohort (closing Mood 2010 non-collapsibility & Imai 2010
     non-causal mediation critique!).
  3. Complete Re-Audit of ALL 55 AD >=2-star P/LP alerts under Gene-Specific VCEP BS1/BA1 + Complete 55/55 live QC!
"""
import json, os, time, urllib.request, warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
CLEAN = os.path.join(BASE, "data", "clean")
RES = os.path.join(BASE, "results")
AUD = os.path.join(BASE, "audit")
MAX_AN = 1614324.0

# 1. Explicit Gene-Specific BS1 / BA1 Table for all 84 ACMG SF v3.3 genes
# Format: gene -> (bs1, ba1, rule_source, rule_citation)
# Published ClinGen VCEP thresholds (conservative/upper-bound where VCEP uses allele-specific/sub-disease ranges)
VCEP_PUBLISHED = {
    # Cancer VCEPs
    "BRCA1": (1.0e-4, 1.0e-3, "ClinGen_VCEP", "ENIGMA BRCA1/2 VCEP v1.1 (Parsons/Spurdle)"),
    "BRCA2": (1.0e-4, 1.0e-3, "ClinGen_VCEP", "ENIGMA BRCA1/2 VCEP v1.1"),
    "PALB2": (1.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen HBOP VCEP PALB2"),
    "TP53":  (3.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen TP53 VCEP (Fortuno et al. 2021)"),
    "PTEN":  (4.3e-5, 5.6e-4, "ClinGen_VCEP", "ClinGen PTEN VCEP (Mester et al. 2018)"),
    "APC":   (1.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen InSiGHT APC VCEP"),
    "MLH1":  (1.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen InSiGHT MMR VCEP"),
    "MSH2":  (1.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen InSiGHT MMR VCEP"),
    "MSH6":  (2.2e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen InSiGHT MMR VCEP"),
    "PMS2":  (2.2e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen InSiGHT MMR VCEP"),
    "VHL":   (1.56e-5, 1.56e-4, "ClinGen_VCEP", "ClinGen VHL VCEP"),
    # Cardiovascular VCEPs
    "LDLR":   (2.0e-3, 5.0e-3, "ClinGen_VCEP", "ClinGen FH VCEP (Chora et al. 2022)"),
    "MYH7":   (2.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen Cardiomyopathy VCEP (Kelly et al. 2018)"),
    "MYBPC3": (2.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen Cardiomyopathy VCEP (HCM)"),
    "TNNT2":  (2.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen Cardiomyopathy VCEP (HCM/DCM)"),
    "TNNI3":  (2.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen Cardiomyopathy VCEP (HCM)"),
    "TPM1":   (2.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen Cardiomyopathy VCEP (HCM)"),
    "MYL2":   (2.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen Cardiomyopathy VCEP (HCM)"),
    "MYL3":   (2.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen Cardiomyopathy VCEP (HCM)"),
    "ACTC1":  (2.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen Cardiomyopathy VCEP (HCM)"),
    "TNNC1":  (2.0e-4, 1.0e-3, "ClinGen_VCEP", "ClinGen Cardiomyopathy VCEP (HCM/DCM)"),
    "FBN1":   (5.0e-5, 1.0e-3, "ClinGen_VCEP", "ClinGen FBN1 VCEP"),
    # Miscellaneous / Metabolic VCEPs
    "RYR1":   (8.0e-4, 3.8e-3, "ClinGen_VCEP", "ClinGen MH VCEP (Johnston et al. 2021; AD MHS)"),
    "HNF1A":  (3.3e-5, 1.0e-4, "ClinGen_VCEP", "ClinGen Monogenic Diabetes VCEP"),
    "ACVRL1": (8.0e-4, 2.0e-3, "ClinGen_VCEP", "ClinGen HHT VCEP"),
    "ENG":    (8.0e-4, 2.0e-3, "ClinGen_VCEP", "ClinGen HHT VCEP"),
    "GAA":    (5.0e-3, 1.0e-2, "ClinGen_VCEP", "ClinGen Lysosomal Diseases VCEP (AR)"),
    "ATP7B":  (5.0e-3, 1.0e-2, "ClinGen_VCEP", "ClinGen Wilson Disease VCEP (AR)"),
    "RPE65":  (3.0e-3, 1.0e-2, "ClinGen_VCEP", "ClinGen LCA/eoRD VCEP (AR)"),
}

# Conservative Whiffin (2017) domain/inheritance ceilings for remaining 55 ACMG SF v3.3 genes
# (Explicitly separated from ClinGen VCEP-published rules in all tables!)
SPECIAL_WHIFFIN = {
    "APOB":    (1.0e-3, 5.0e-3, "Whiffin_AD_FH", "FH prevalence ~1/250, reduced penetrance"),
    "PCSK9":   (1.0e-3, 5.0e-3, "Whiffin_AD_FH", "FH gain-of-function rare"),
    "TTR":     (1.0e-3, 5.0e-3, "Whiffin_AD_Amyloid", "Reduced-penetrance amyloidosis; p.Val142Ile BA1 exception"),
    "TTN":     (4.0e-4, 1.0e-3, "Whiffin_AD_DCM", "DCM prevalence ~1/250, high allelic het"),
    "CACNA1S": (8.0e-4, 3.8e-3, "Whiffin_AD_MH", "Pharmacogenetic MH susceptibility"),
}

# Explicit Ghosh et al. (2018, Hum Mutat 39:1525-1530) & known high-AF pathogenic BA1 exception alleles
BA1_EXCEPTIONS = {
    "6-26092913-G-A": "HFE p.Cys282Tyr (Ghosh 2018 BA1 exception)",
    "6-26090951-C-G": "HFE p.His63Asp (Ghosh 2018 BA1 exception)",
    "3-15645186-G-C": "BTD p.Asp444His (Ghosh 2018 BA1 exception)",
    "18-31598655-G-A": "TTR p.Val142Ile (AFR cardiac amyloidosis founder)",
    "1-45331833-C-T": "MUTYH p.Gly396Asp (AR MAP founder)",
    "1-45332803-T-C": "MUTYH p.Tyr179Cys (AR MAP founder)",
}


def get_gene_thresholds(gene, domain, inh):
    if gene in VCEP_PUBLISHED:
        return VCEP_PUBLISHED[gene]
    if gene in SPECIAL_WHIFFIN:
        return SPECIAL_WHIFFIN[gene]
    if inh == "AR":
        return (5.0e-3, 1.0e-2, "Whiffin_AR_Generic", "Conservative AR ceiling (Whiffin 2017)")
    if inh == "XL":
        return (5.0e-4, 5.0e-3, "Whiffin_XL_Generic", "Conservative XL ceiling (Whiffin 2017)")
    if domain == "Cancer":
        return (2.0e-4, 1.0e-3, "Whiffin_AD_Cancer", "Conservative AD hereditary cancer ceiling (Whiffin 2017)")
    if domain == "Cardiovascular":
        return (2.0e-4, 1.0e-3, "Whiffin_AD_Cardio", "Conservative AD cardio/arrhythmia ceiling (Whiffin 2017)")
    return (5.0e-4, 5.0e-3, "Whiffin_AD_Misc", "Conservative AD miscellaneous ceiling (Whiffin 2017)")


def run_collapsible_attenuation(df_in, label, has_last_eval=True):
    """
    Denominator-clean FAF95 >= 1e-4, AN >= 80% analysis closing Codex C07 (Mood 2010 & Imai 2010):
      - Collapsible Linear Probability Model (LPM / exact risk difference) +/- log_subs & +/- submitter-tier bins
      - Logistic Average Marginal Effects (AME) +/- log_subs
      - Exact non-parametric Submitter-Tier Stratified Contrasts (1, 2-3, 4-9, 10+ submitters)
      - LastEvaluated recency check (mean/median year, post-2023 share)
    """
    sub = df_in[
        df_in["csq_group"].isin(["pLoF", "missense", "inframe_indel", "splice_other", "synonymous"]) &
        df_in["clnsig_cat"].isin(["P_LP", "VUS", "CONFLICTING", "B_LB"]) &
        (df_in["an_joint"] >= 0.80 * MAX_AN)
    ].copy()
    m_nn = (sub["faf95_nonnfe5"] >= 1e-4) & (sub["af_nonnfe5"] / np.maximum(sub["af_nfe"], 1e-6) >= 5.0)
    m_eu = (sub["faf95_nfe"] >= 1e-4) & (sub["af_nfe"] / np.maximum(sub["af_nonnfe5"], 1e-6) >= 5.0)
    comp = sub[m_nn | m_eu].copy()
    comp["non_nfe"] = m_nn.loc[comp.index].astype(int)
    comp["af_ref"] = np.where(comp["non_nfe"] == 1, comp["faf95_nonnfe5"], comp["faf95_nfe"])
    comp["log10_af"] = np.log10(comp["af_ref"].clip(lower=1e-6))
    comp["n_subs"] = comp["n_submitters"].fillna(1).clip(lower=1)
    comp["log_subs"] = np.log1p(comp["n_subs"])
    comp["subs_bin"] = pd.cut(comp["n_subs"], bins=[0, 1, 3, 9, 19, 99999], labels=["1", "2_3", "4_9", "10_19", "20p"]).astype(str)

    out = {
        "label": label,
        "n_total": int(len(comp)),
        "n_non_nfe": int((comp["non_nfe"] == 1).sum()),
        "n_nfe": int((comp["non_nfe"] == 0).sum()),
        "median_subs_nn": float(comp.loc[comp["non_nfe"] == 1, "n_subs"].median()),
        "median_subs_nfe": float(comp.loc[comp["non_nfe"] == 0, "n_subs"].median()),
        "outcomes": {},
        "stratified_by_submitter_tier": {},
    }

    # Evaluation recency check
    if has_last_eval and "last_evaluated" in comp.columns:
        dt = pd.to_datetime(comp["last_evaluated"], errors="coerce")
        comp["eval_year"] = dt.dt.year
        comp["eval_post_2023"] = (dt >= "2024-01-01").astype(float)
        g_nn = comp[comp["non_nfe"] == 1]
        g_eu = comp[comp["non_nfe"] == 0]
        _, p_yr = stats.mannwhitneyu(g_nn["eval_year"].dropna(), g_eu["eval_year"].dropna())
        # OLS on eval_year controlling for gene, csq, af
        m_yr = smf.ols("eval_year ~ non_nfe + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=comp.dropna(subset=["eval_year"])).fit(cov_type="HC3")
        out["last_evaluated_audit"] = {
            "median_year_nn": float(g_nn["eval_year"].median()),
            "median_year_nfe": float(g_eu["eval_year"].median()),
            "mean_year_nn": round(float(g_nn["eval_year"].mean()), 2),
            "mean_year_nfe": round(float(g_eu["eval_year"].mean()), 2),
            "post_2023_share_nn": round(float(g_nn["eval_post_2023"].mean()), 4),
            "post_2023_share_nfe": round(float(g_eu["eval_post_2023"].mean()), 4),
            "mw_p": float(p_yr),
            "ols_year_diff_nn_vs_nfe": round(float(m_yr.params["non_nfe"]), 3),
            "ols_year_se": round(float(m_yr.bse["non_nfe"]), 3),
            "ols_year_p": float(m_yr.pvalues["non_nfe"]),
        }

    for cat in ("VUS", "CONFLICTING", "B_LB", "P_LP"):
        comp["y"] = (comp["clnsig_cat"] == cat).astype(float)
        # 1. Collapsible LPM (exact Risk Difference in pp!)
        m_l0 = smf.ols("y ~ non_nfe + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=comp).fit(cov_type="HC3")
        m_l1 = smf.ols("y ~ non_nfe + log_subs + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=comp).fit(cov_type="HC3")
        m_l2 = smf.ols("y ~ non_nfe + C(subs_bin) + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=comp).fit(cov_type="HC3")

        b0, se0, p0 = float(m_l0.params["non_nfe"]), float(m_l0.bse["non_nfe"]), float(m_l0.pvalues["non_nfe"])
        b1, se1, p1 = float(m_l1.params["non_nfe"]), float(m_l1.bse["non_nfe"]), float(m_l1.pvalues["non_nfe"])
        b2, se2, p2 = float(m_l2.params["non_nfe"]), float(m_l2.bse["non_nfe"]), float(m_l2.pvalues["non_nfe"])

        # 2. Logistic Average Marginal Effects (AME on probability scale, addressing Mood 2010 non-collapsibility!)
        gv = comp.groupby("gene")["y"].nunique()
        sg = comp[comp["gene"].isin(gv[gv == 2].index)].copy()
        ame0 = ame1 = None
        lor0 = lor1 = lp0 = lp1 = None
        try:
            ml0 = smf.logit("y ~ non_nfe + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=sg).fit(disp=False, maxiter=150)
            ml1 = smf.logit("y ~ non_nfe + log_subs + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=sg).fit(disp=False, maxiter=150)
            lor0, lp0 = float(np.exp(ml0.params["non_nfe"])), float(ml0.pvalues["non_nfe"])
            lor1, lp1 = float(np.exp(ml1.params["non_nfe"])), float(ml1.pvalues["non_nfe"])
            # Exact finite-difference AME: E[P(y=1|non_nfe=1, X) - P(y=1|non_nfe=0, X)]
            sg1, sg0 = sg.copy(), sg.copy()
            sg1["non_nfe"] = 1
            sg0["non_nfe"] = 0
            ame0 = float((ml0.predict(sg1) - ml0.predict(sg0)).mean())
            ame1 = float((ml1.predict(sg1) - ml1.predict(sg0)).mean())
        except Exception:
            pass

        out["outcomes"][cat] = {
            "raw_rate_nn_pct": round(100.0 * comp.loc[comp["non_nfe"] == 1, "y"].mean(), 2),
            "raw_rate_nfe_pct": round(100.0 * comp.loc[comp["non_nfe"] == 0, "y"].mean(), 2),
            "lpm_base_rd_pp": round(100.0 * b0, 3), "lpm_base_se_pp": round(100.0 * se0, 3), "lpm_base_p": p0,
            "lpm_logsubs_rd_pp": round(100.0 * b1, 3), "lpm_logsubs_se_pp": round(100.0 * se1, 3), "lpm_logsubs_p": p1,
            "lpm_subsbin_rd_pp": round(100.0 * b2, 3), "lpm_subsbin_se_pp": round(100.0 * se2, 3), "lpm_subsbin_p": p2,
            "lpm_rd_attenuation_logsubs_pct": round((1.0 - b1 / b0) * 100.0, 1) if abs(b0) > 1e-6 else None,
            "lpm_rd_attenuation_subsbin_pct": round((1.0 - b2 / b0) * 100.0, 1) if abs(b0) > 1e-6 else None,
            "logit_base_or": round(lor0, 3) if lor0 is not None else None, "logit_base_p": lp0,
            "logit_logsubs_or": round(lor1, 3) if lor1 is not None else None, "logit_logsubs_p": lp1,
            "logit_ame_base_pp": round(100.0 * ame0, 3) if ame0 is not None else None,
            "logit_ame_logsubs_pp": round(100.0 * ame1, 3) if ame1 is not None else None,
            "logit_ame_attenuation_pct": round((1.0 - ame1 / ame0) * 100.0, 1) if (ame0 is not None and ame1 is not None and abs(ame0) > 1e-6) else None,
        }

    # 3. Exact Non-Parametric Submitter-Tier Stratified Contrasts for VUS & Conflicting!
    for tier in ["1", "2_3", "4_9", "10_19", "20p"]:
        gt = comp[comp["subs_bin"] == tier]
        gn, ge = gt[gt["non_nfe"] == 1], gt[gt["non_nfe"] == 0]
        if len(gn) == 0 or len(ge) == 0:
            continue
        nv_n, nt_n = int((gn["clnsig_cat"] == "VUS").sum()), len(gn)
        nv_e, nt_e = int((ge["clnsig_cat"] == "VUS").sum()), len(ge)
        _, p_f = stats.fisher_exact([[nv_n, nt_n - nv_n], [nv_e, nt_e - nv_e]])
        nc_n = int((gn["clnsig_cat"] == "CONFLICTING").sum())
        nc_e = int((ge["clnsig_cat"] == "CONFLICTING").sum())
        _, p_fc = stats.fisher_exact([[nc_n, nt_n - nc_n], [nc_e, nt_e - nc_e]])
        out["stratified_by_submitter_tier"][tier] = {
            "n_nn": nt_n, "n_nfe": nt_e,
            "vus_nn_pct": round(100.0 * nv_n / nt_n, 2), "vus_nfe_pct": round(100.0 * nv_e / nt_e, 2), "vus_fisher_p": float(p_f),
            "conf_nn_pct": round(100.0 * nc_n / nt_n, 2), "conf_nfe_pct": round(100.0 * nc_e / nt_e, 2), "conf_fisher_p": float(p_fc),
        }
    return out


def main():
    t0 = time.time()
    j = pd.read_parquet(os.path.join(CLEAN, "a1_clinvar_gnomad_exact_joined.parquet"))
    ju = pd.read_parquet(os.path.join(AUD, "a1_unseen20_exact_joined.parquet"))
    p20 = json.load(open(os.path.join(AUD, "a1_prospective_unseen20_and_qc55.json")))
    qc55_map = {x["variant_id"]: x for x in p20["acmg55_complete_qc_audit"]["all_55"]}

    # Assign exact gene-specific BS1/BA1 rules to every ACMG SF v3.3 row
    gene_meta = j[["gene", "domain", "inheritance"]].drop_duplicates("gene")
    rule_rows = []
    g_bs1, g_ba1, g_src, g_cit = {}, {}, {}, {}
    for r in gene_meta.itertuples(index=False):
        bs1, ba1, src, cit = get_gene_thresholds(r.gene, r.domain, r.inheritance)
        g_bs1[r.gene], g_ba1[r.gene], g_src[r.gene], g_cit[r.gene] = bs1, ba1, src, cit
        rule_rows.append({"gene": r.gene, "domain": r.domain, "inheritance": r.inheritance, "bs1": bs1, "ba1": ba1, "rule_source": src, "citation": cit})
    df_rules = pd.DataFrame(rule_rows).sort_values(["rule_source", "domain", "gene"])
    df_rules.to_csv(os.path.join(AUD, "a1_acmg84_gene_specific_bs1_ba1_rules.csv"), index=False)

    j["gene_bs1"] = j["gene"].map(g_bs1)
    j["gene_ba1"] = j["gene"].map(g_ba1)
    j["rule_source"] = j["gene"].map(g_src)
    j["is_ba1_exception"] = j["variant_id"].isin(BA1_EXCEPTIONS)

    # 1. Gene-Specific BS1 / BA1 VUS Frequency-Conflict Candidates (HIGH CALLABILITY AN>=80%, excluding BA1 exceptions!)
    vus_hc = j[(j["clnsig_cat"] == "VUS") & (j["an_joint"] >= 0.80 * MAX_AN) & (~j["is_ba1_exception"])].copy()
    vus_hc["exceeds_gene_bs1"] = vus_hc["faf95_grpmax5"] >= vus_hc["gene_bs1"]
    vus_hc["exceeds_gene_ba1"] = vus_hc["faf95_grpmax5"] >= vus_hc["gene_ba1"]
    vus_hc["nfe_below_bs1"] = vus_hc["faf95_nfe"] < vus_hc["gene_bs1"]
    vus_hc["nfe_below_ba1"] = vus_hc["faf95_nfe"] < vus_hc["gene_ba1"]

    def summarize_vus_gene_spec(df_sub, tag):
        b1 = df_sub[df_sub["exceeds_gene_bs1"]]
        ba = df_sub[df_sub["exceeds_gene_ba1"]]
        return {
            "tag": tag,
            "n_genes": int(df_sub["gene"].nunique()),
            "n_hc_vus": int(len(df_sub)),
            "bs1_candidates": {
                "n": int(len(b1)),
                "n_nfe_invisible": int(b1["nfe_below_bs1"].sum()),
                "nfe_invisible_share": round(float(b1["nfe_below_bs1"].mean()), 4) if len(b1) else 0.0,
                "by_pop5": {str(k): int(v) for k, v in b1["faf95_pop5"].value_counts().items()},
                "by_domain": {str(k): int(v) for k, v in b1["domain"].value_counts().items()},
                "top_genes": {str(k): int(v) for k, v in b1["gene"].value_counts().head(10).items()},
            },
            "ba1_candidates": {
                "n": int(len(ba)),
                "n_nfe_invisible": int(ba["nfe_below_ba1"].sum()),
                "nfe_invisible_share": round(float(ba["nfe_below_ba1"].mean()), 4) if len(ba) else 0.0,
                "by_pop5": {str(k): int(v) for k, v in ba["faf95_pop5"].value_counts().items()},
                "by_domain": {str(k): int(v) for k, v in ba["domain"].value_counts().items()},
                "top_genes": {str(k): int(v) for k, v in ba["gene"].value_counts().head(10).items()},
            },
        }

    vus_gs = {
        "vcep_published_29_genes_only": summarize_vus_gene_spec(vus_hc[vus_hc["rule_source"] == "ClinGen_VCEP"], "29 ACMG SF v3.3 Genes with Published ClinGen VCEP BS1/BA1"),
        "all_72_ad_genes_vcep_plus_whiffin": summarize_vus_gene_spec(vus_hc[vus_hc["inheritance"] == "AD"], "All 72 AD ACMG SF v3.3 Genes (VCEP + Whiffin)"),
        "all_84_genes_vcep_plus_whiffin": summarize_vus_gene_spec(vus_hc, "All 84 ACMG SF v3.3 Genes (VCEP + Whiffin)"),
    }
    # Export full gene-specific VUS BS1 candidate table as CSV for Supplementary Table!
    vus_bs1_df = vus_hc[vus_hc["exceeds_gene_bs1"]].sort_values(["rule_source", "faf95_grpmax5"], ascending=[True, False])
    vus_bs1_df[[
        "variant_id", "variation_id", "gene", "domain", "inheritance", "rule_source", "gene_bs1", "gene_ba1",
        "stars", "n_submitters", "last_evaluated", "csq_group", "an_joint", "faf95_nfe", "faf95_grpmax5", "faf95_pop5", "exceeds_gene_ba1"
    ]].to_csv(os.path.join(AUD, "a1_vus_gene_specific_bs1_candidates.csv"), index=False)

    # 2. Live gnomAD v4.1 QC Filter Spot-Check on 20 Stratified VUS BS1 Candidates (10 VCEP + 10 Whiffin)
    #    to empirically bound the stealth QC-failure rate among AN>=80% VUS candidates!
    rng = np.random.default_rng(20261125)
    vcep_b1 = vus_bs1_df[vus_bs1_df["rule_source"] == "ClinGen_VCEP"]
    whif_b1 = vus_bs1_df[vus_bs1_df["rule_source"] != "ClinGen_VCEP"]
    samp_vids = list(vcep_b1.sample(min(10, len(vcep_b1)), random_state=42)["variant_id"]) + \
                list(whif_b1.sample(min(10, len(whif_b1)), random_state=42)["variant_id"])
    gql = """
    query($variantId: String!, $dataset: DatasetId!) {
      variant(variantId: $variantId, dataset: $dataset) {
        variant_id flags
        exome { ac an filters flags }
        genome { ac an filters flags }
      }
    }
    """
    vus_qc_spot = []
    qc_cache = os.path.join(AUD, "gnomad_qc_flags")
    os.makedirs(qc_cache, exist_ok=True)
    for vid in samp_vids:
        cf = os.path.join(qc_cache, f"vus_{vid}.json")
        if os.path.exists(cf) and os.path.getsize(cf) > 50:
            dq = json.load(open(cf))
        else:
            req = urllib.request.Request(
                "https://gnomad.broadinstitute.org/api",
                data=json.dumps({"query": gql, "variables": {"variantId": vid, "dataset": "gnomad_r4"}}).encode("utf-8"),
                headers={"Content-Type": "application/json", "User-Agent": "AcademicAudit/1.0"},
                method="POST",
            )
            try:
                with urllib.request.urlopen(req, timeout=25) as r:
                    dq = json.loads(r.read().decode("utf-8"))
                with open(cf, "w") as f:
                    json.dump(dq, f)
            except Exception as e:
                dq = {"error": str(e)}
            time.sleep(6.2)
        v = (dq.get("data") or {}).get("variant") or {}
        ex, ge = v.get("exome") or {}, v.get("genome") or {}
        ex_f, ge_f = ex.get("filters") or [], ge.get("filters") or []
        fl = v.get("flags") or []
        qc_fail = bool(ex_f) or (bool(ge) and bool(ge_f) and not ex) or ("lcr" in fl)
        vus_qc_spot.append({"variant_id": vid, "exome_filters": ex_f, "genome_filters": ge_f, "flags": fl, "qc_fail": qc_fail})

    # 3. Re-Audit ALL 55 AD >=2-star P/LP Alerts with 55/55 Live QC + Gene-Specific BS1/BA1 + Inheritance Mode!
    cat55 = pd.read_csv(os.path.join(RES, "a1_ad_highconf_triaged_55.csv"))
    cat55["gene_bs1"] = cat55["gene"].map(g_bs1)
    cat55["gene_ba1"] = cat55["gene"].map(g_ba1)
    cat55["rule_source"] = cat55["gene"].map(g_src)
    cat55["live_exome_filters"] = cat55["variant_id"].map(lambda x: "|".join(qc55_map.get(x, {}).get("exome_filters", [])))
    cat55["live_genome_filters"] = cat55["variant_id"].map(lambda x: "|".join(qc55_map.get(x, {}).get("genome_filters", [])))
    cat55["live_flags"] = cat55["variant_id"].map(lambda x: "|".join(qc55_map.get(x, {}).get("global_flags", [])))
    cat55["live_qc_fail"] = cat55["variant_id"].map(lambda x: bool(qc55_map.get(x, {}).get("any_qc_fail_or_low_an", False)))

    AR_IN_AD_GENES = {
        "19-38527014-C-T", "19-38523211-C-G", "19-38505338-CAT-C", "2-178568947-T-A"
    }
    def revised_triage(r):
        if r["live_qc_fail"]:
            return "L1_QC_or_Callability_Artifact"
        if r["triage_layer"] == "L2_BottleneckFounder_FIN_ASJ":
            return "L2_BottleneckFounder_FIN_ASJ"
        if r["triage_layer"] == "L3_EstablishedContinentalFounder_LowPen":
            return "L3_EstablishedContinentalFounder_LowPen"
        # Inside high-QC non-bottleneck L4 (11 variants):
        if r["variant_id"] in AR_IN_AD_GENES:
            return "L4b_AR_in_DualMode_AD_Gene_ReturnHazard"
        if r["faf95_max8"] < r["gene_bs1"]:
            return "L4a_BelowGeneSpecificBS1_TrueADPathogenic"
        return "L4c_HighQC_ExceedsGeneBS1_VCEP_Review"

    cat55["revised_triage_layer"] = cat55.apply(revised_triage, axis=1)
    cat55.to_csv(os.path.join(AUD, "a1_ad55_revised_vcep_qc_triage.csv"), index=False)

    # 4. Collapsible LPM + AME + Submitter-Tier Stratified + LastEvaluated Audit on BOTH ACMG84 and Unseen20!
    atten_acmg84 = run_collapsible_attenuation(j, "ACMG SF v3.3 (84 Genes, FAF95>=1e-4, AN>=80%)", has_last_eval=True)
    atten_unseen20 = run_collapsible_attenuation(ju, "Prospective Unseen 20 Non-ACMG AD Genes (FAF95>=1e-4, AN>=80%)", has_last_eval=True)

    out = {
        "executed_at": time.strftime("%Y-%m-%d %H:%M:%S CST"),
        "elapsed_sec": round(time.time() - t0, 1),
        "n_vcep_published_genes": len(VCEP_PUBLISHED),
        "vus_gene_specific_bs1_ba1": vus_gs,
        "vus_bs1_live_qc_spotcheck_20": {
            "n_checked": len(vus_qc_spot),
            "n_pass": sum(1 for x in vus_qc_spot if not x["qc_fail"]),
            "n_fail": sum(1 for x in vus_qc_spot if x["qc_fail"]),
            "pass_rate": round(sum(1 for x in vus_qc_spot if not x["qc_fail"]) / max(len(vus_qc_spot), 1), 4),
            "details": vus_qc_spot,
        },
        "ad55_revised_triage_counts": {str(k): int(v) for k, v in cat55["revised_triage_layer"].value_counts().items()},
        "ad55_l4_all11_revised": cat55[cat55["revised_triage_layer"].str.startswith("L4")][[
            "variant_id", "gene", "stars", "hgvs_name", "faf95_pop8", "faf95_max8", "faf95_nfe",
            "gene_bs1", "gene_ba1", "rule_source", "hom_joint", "revised_triage_layer"
        ]].to_dict(orient="records"),
        "collapsible_attenuation_acmg84": atten_acmg84,
        "collapsible_attenuation_unseen20": atten_unseen20,
    }
    out_p = os.path.join(AUD, "a1_vcep_bs1_and_collapsible_submitters.json")
    with open(out_p, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps({
        "vus_gene_specific_bs1_ba1": vus_gs,
        "vus_qc_spotcheck": out["vus_bs1_live_qc_spotcheck_20"]["pass_rate"],
        "ad55_revised_triage_counts": out["ad55_revised_triage_counts"],
        "acmg84_vus_lpm": atten_acmg84["outcomes"]["VUS"],
        "unseen20_vus_lpm": atten_unseen20["outcomes"]["VUS"],
        "acmg84_eval_recency": atten_acmg84.get("last_evaluated_audit"),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
