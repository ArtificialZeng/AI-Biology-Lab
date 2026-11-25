#!/opt/anaconda3/bin/python3
"""
Full Publication-Grade Allele-Exact Reclassification Debt & Cross-Ancestry Equity Engine for Study A-1:
ClinVar GRCh38 (Nov 2026) x gnomAD v4.1 (807,162 individuals) across ALL 84 ACMG SF v3.3 Genes.

Implements pre-registered protocol PREREG_A1.md:
  - 100% allele-exact (CHROM-POS-REF-ALT) GRCh38 join (ZERO rsID collisions)
  - Exact Joint Exome+Genome counts & exact Whiffin/ClinGen 95% Filtering Allele Frequency (FAF95)
    across all 8 ancestry groups (nfe, afr, amr, eas, sas, mid, asj, fin) and official 5 non-bottleneck Grpmax
  - H1: Frequency-matched CMH stratified ORs + Gene/Consequence/AF-spline Logistic Regressions (VUS, B/LB, Conflicting, P/LP)
        + Ancestry-by-ancestry diagnostic equity profiles
  - H2: High-Confidence (>=2-star) AD Reclassification Debt across 72 AD genes (Tier 1/2/3 FAF95 & homozygote ceilings),
        NFE-invisibility rate, founder-excluded sensitivity, and individual-level carrier burden by ancestry
  - H3: Review-staleness & post-gnomAD-v4 (2024-2026) governance audit (including 3-star Expert Panel variants)
  - H4: Quantitative audit of multiallelic rsID join false-positive inflation across all 84 genes
  - Spot cross-validation of computed FAF95 against gnomAD native joint.faf95 via official get_variant_frequency.py
"""
import glob
import json
import os
import subprocess
import time
import warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
RAW = os.path.join(BASE, "data", "raw")
GNOMAD_DIR = os.path.join(RAW, "gnomad_v4_acmg84")
CLEAN = os.path.join(BASE, "data", "clean")
RES = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

GNOMAD_FREQ_SCRIPT = "/Users/mac/.gemini/config/plugins/genomics/skills/gnomad_database/scripts/get_variant_frequency.py"

POPS_ALL8 = ["nfe", "afr", "amr", "eas", "sas", "mid", "asj", "fin"]
POPS_GRPMAX5 = ["nfe", "afr", "amr", "eas", "sas"]
POPS_NONNFE5 = ["afr", "amr", "eas", "sas", "mid"]
POPS_NONNFE7 = ["afr", "amr", "eas", "sas", "mid", "asj", "fin"]

# Well-established true high-frequency pathogenic founder or low-penetrance alleles in ACMG SF v3.3
# (flagged transparently so we report BOTH Gross and Founder-Excluded Reclassification Debt!)
KNOWN_FOUNDER_KEYWORDS = [
    ("TTR", "Val142Ile"), ("TTR", "V142I"), ("TTR", "18-31598655-G-A"),  # AFR cardiac amyloidosis founder (~1.6% AFR)
    ("APOB", "Arg3527Gln"), ("APOB", "2-21006288-C-T"),                 # European FH founder
    ("APC", "Ile1307Lys"), ("APC", "5-112839514-T-A"),                  # ASJ reduced-penetrance risk allele
    ("BRCA1", "17-43124027-ACT-A"), ("BRCA1", "17-43057051-T-TG"),      # ASJ 185delAG, 5382insC
    ("BRCA2", "13-32340300-GT-G"),                                       # ASJ 6174delT
    ("HFE", "6-26092913-G-A"),                                           # NFE p.Cys282Tyr
    ("KCNQ1", "11-2571161-G-A"),                                         # Finnish founder G589D
    ("MYBPC3", "11-47333228"),                                           # SAS 25bp deletion
    ("FH", ""),
]


def calc_faf95_vec(ac, an):
    """Exact Whiffin et al. (Genet Med 2017) 95% lower-bound Filtering Allele Frequency."""
    faf = np.zeros_like(ac, dtype=float)
    mask = (ac >= 2) & (an >= 200)
    if mask.any():
        # Exact Poisson lower 95% CI on (AC - 1) divided by AN: 0.5 * chi2.ppf(0.05, 2*(ac-1)) / an
        faf[mask] = 0.5 * stats.chi2.ppf(0.05, 2.0 * (ac[mask] - 1.0)) / an[mask]
    return faf


def classify_clnsig(sig_raw, conf_raw=""):
    s = (sig_raw or "").lower()
    if not s:
        return "OTHER"
    if "conflicting" in s:
        return "CONFLICTING"
    # Separate explicit low_penetrance / established_risk_allele
    is_low_pen = ("low_penetrance" in s) or ("established_risk" in s) or ("risk_factor" in s)
    if ("pathogenic" in s) and ("benign" not in s) and ("uncertain" not in s):
        return "P_LP_LOWPEN" if is_low_pen else "P_LP"
    if "uncertain_significance" in s:
        return "VUS"
    if ("benign" in s) and ("pathogenic" not in s) and ("uncertain" not in s):
        return "B_LB"
    return "OTHER"


def classify_consequence(mc_raw, gnomad_csq=""):
    txt = f"{mc_raw or ''} {gnomad_csq or ''}".lower()
    if any(k in txt for k in ("nonsense", "stop_gained", "frameshift", "splice_donor", "splice_acceptor", "start_lost", "initiator_codon")):
        return "pLoF"
    if "missense" in txt:
        return "missense"
    if "inframe" in txt:
        return "inframe_indel"
    if "splice_region" in txt or "splice" in txt:
        return "splice_other"
    if "synonymous" in txt:
        return "synonymous"
    if "utr" in txt or "intron" in txt or "non_coding" in txt:
        return "noncoding"
    return "other"


def is_known_founder(row):
    g = str(row.get("gene", ""))
    vid = str(row.get("variant_id", ""))
    hgvs = str(row.get("hgvs_name", ""))
    for kg, pat in KNOWN_FOUNDER_KEYWORDS:
        if kg == g and pat and (pat in vid or pat in hgvs):
            return True
    # Also flag if ClinVar explicitly mentions founder/low penetrance in CLNSIG
    if row.get("clnsig_cat") == "P_LP_LOWPEN":
        return True
    return False


def load_gnomad_acmg84():
    """Flatten all 84 gnomAD v4.1 JSON files into a compact allele-exact DataFrame."""
    cache_pq = os.path.join(CLEAN, "gnomad_v4_acmg84_flat.parquet")
    if os.path.exists(cache_pq) and os.path.getsize(cache_pq) > 100000:
        print("Loading cached gnomad_v4_acmg84_flat.parquet...", flush=True)
        return pd.read_parquet(cache_pq)

    t0 = time.time()
    files = sorted(glob.glob(os.path.join(GNOMAD_DIR, "*_gnomad_r4.json")))
    print(f"Flattening {len(files)} gnomAD v4.1 gene files...", flush=True)
    rows = []
    seen = set()
    for fp in files:
        gene = os.path.basename(fp).replace("_gnomad_r4.json", "")
        with open(fp) as f:
            d = json.load(f)
        vs = (d.get("data", {}).get("gene") or {}).get("variants")
        if vs is None:
            vs = d.get("variants_chunked") or []
        for v in vs:
            vid = v.get("variant_id")
            if not vid or vid in seen:
                continue
            seen.add(vid)
            ex = v.get("exome") or {}
            ge = v.get("genome") or {}
            rec = {
                "variant_id": vid,
                "gnomad_gene": gene,
                "gnomad_csq": v.get("consequence") or "",
                "gnomad_rsids": "|".join(v.get("rsids") or []),
                "ac_joint": int(ex.get("ac") or 0) + int(ge.get("ac") or 0),
                "an_joint": int(ex.get("an") or 0) + int(ge.get("an") or 0),
                "hom_joint": int(ex.get("homozygote_count") or 0) + int(ge.get("homozygote_count") or 0),
                "hemi_joint": int(ex.get("hemizygote_count") or 0) + int(ge.get("hemizygote_count") or 0),
            }
            # Population dicts
            p_ac = {p: 0 for p in POPS_ALL8}
            p_an = {p: 0 for p in POPS_ALL8}
            p_hom = {p: 0 for p in POPS_ALL8}
            for src in (ex, ge):
                for pop in (src.get("populations") or []):
                    pid = pop.get("id")
                    if pid in p_ac:
                        p_ac[pid] += int(pop.get("ac") or 0)
                        p_an[pid] += int(pop.get("an") or 0)
                        p_hom[pid] += int(pop.get("homozygote_count") or 0)
            for p in POPS_ALL8:
                rec[f"ac_{p}"] = p_ac[p]
                rec[f"an_{p}"] = p_an[p]
                rec[f"hom_{p}"] = p_hom[p]
            rows.append(rec)

    gdf = pd.DataFrame(rows)
    gdf["af_joint"] = np.where(gdf["an_joint"] > 0, gdf["ac_joint"] / gdf["an_joint"], 0.0)
    for p in POPS_ALL8:
        ac = gdf[f"ac_{p}"].values.astype(float)
        an = gdf[f"an_{p}"].values.astype(float)
        gdf[f"af_{p}"] = np.where(an >= 200, ac / an, 0.0)
        gdf[f"faf95_{p}"] = calc_faf95_vec(ac, an)

    # Vectorized popmax calculations
    af_all8 = gdf[[f"af_{p}" for p in POPS_ALL8]].values
    faf_all8 = gdf[[f"faf95_{p}" for p in POPS_ALL8]].values
    af_g5 = gdf[[f"af_{p}" for p in POPS_GRPMAX5]].values
    faf_g5 = gdf[[f"faf95_{p}" for p in POPS_GRPMAX5]].values
    af_nn5 = gdf[[f"af_{p}" for p in POPS_NONNFE5]].values
    faf_nn5 = gdf[[f"faf95_{p}" for p in POPS_NONNFE5]].values
    af_nn7 = gdf[[f"af_{p}" for p in POPS_NONNFE7]].values
    faf_nn7 = gdf[[f"faf95_{p}" for p in POPS_NONNFE7]].values

    gdf["af_max8"] = af_all8.max(axis=1)
    gdf["pop_max8"] = [POPS_ALL8[i] for i in af_all8.argmax(axis=1)]
    gdf["faf95_max8"] = faf_all8.max(axis=1)
    gdf["faf95_pop8"] = [POPS_ALL8[i] for i in faf_all8.argmax(axis=1)]

    gdf["af_grpmax5"] = af_g5.max(axis=1)
    gdf["pop_grpmax5"] = [POPS_GRPMAX5[i] for i in af_g5.argmax(axis=1)]
    gdf["faf95_grpmax5"] = faf_g5.max(axis=1)
    gdf["faf95_pop5"] = [POPS_GRPMAX5[i] for i in faf_g5.argmax(axis=1)]

    gdf["af_nonnfe5"] = af_nn5.max(axis=1)
    gdf["pop_nonnfe5"] = [POPS_NONNFE5[i] for i in af_nn5.argmax(axis=1)]
    gdf["faf95_nonnfe5"] = faf_nn5.max(axis=1)

    gdf["af_nonnfe7"] = af_nn7.max(axis=1)
    gdf["pop_nonnfe7"] = [POPS_NONNFE7[i] for i in af_nn7.argmax(axis=1)]
    gdf["faf95_nonnfe7"] = faf_nn7.max(axis=1)

    gdf.to_parquet(cache_pq, index=False)
    print(f"Flattened {len(gdf):,} unique gnomAD v4.1 variants across {len(files)} genes in {time.time()-t0:.1f}s", flush=True)
    return gdf


def cmh_test(df_sub, expo_col, out_col, strata_col):
    """Exact Cochran-Mantel-Haenszel stratified OR, 95% CI (Robins-Breslow-Greenland), and p-value."""
    num = 0.0
    den = 0.0
    var_num = 0.0
    exp_a = 0.0
    obs_a = 0.0
    # RBG variance components
    P_sum, Q_sum, R_sum, S_sum = 0.0, 0.0, 0.0, 0.0
    PR, PS_QR, QS = 0.0, 0.0, 0.0

    for _, g in df_sub.groupby(strata_col):
        a = float(((g[expo_col] == 1) & (g[out_col] == 1)).sum())
        b = float(((g[expo_col] == 1) & (g[out_col] == 0)).sum())
        c = float(((g[expo_col] == 0) & (g[out_col] == 1)).sum())
        d = float(((g[expo_col] == 0) & (g[out_col] == 0)).sum())
        n = a + b + c + d
        if n <= 1:
            continue
        r1, r0, c1, c0 = a + b, c + d, a + c, b + d
        if r1 == 0 or r0 == 0 or c1 == 0 or c0 == 0:
            continue
        obs_a += a
        exp_a += r1 * c1 / n
        var_num += (r1 * r0 * c1 * c0) / (n * n * (n - 1.0))
        R = a * d / n
        S = b * c / n
        P = (a + d) / n
        Q = (b + c) / n
        R_sum += R
        S_sum += S
        PR += P * R
        PS_QR += P * S + Q * R
        QS += Q * S

    if S_sum <= 0 or R_sum <= 0 or var_num <= 0:
        return {"or": np.nan, "ci95": [np.nan, np.nan], "p": np.nan, "z": np.nan}
    or_mh = R_sum / S_sum
    var_log_or = (
        PR / (2.0 * R_sum * R_sum)
        + PS_QR / (2.0 * R_sum * S_sum)
        + QS / (2.0 * S_sum * S_sum)
    )
    se_log = np.sqrt(max(var_log_or, 1e-12))
    ci95 = [float(np.exp(np.log(or_mh) - 1.96 * se_log)), float(np.exp(np.log(or_mh) + 1.96 * se_log))]
    z = (obs_a - exp_a) / np.sqrt(var_num)
    p = float(2.0 * stats.norm.sf(abs(z)))
    return {"or": float(or_mh), "ci95": ci95, "p": p, "z": float(z)}


def run_h1_equity(joined):
    """H1: Frequency-matched cross-ancestry classification asymmetry on allele-exact coding/splice variants."""
    cs_mask = joined["csq_group"].isin(["pLoF", "missense", "inframe_indel", "splice_other", "synonymous"])
    sub = joined[cs_mask & joined["clnsig_cat"].isin(["P_LP", "VUS", "CONFLICTING", "B_LB"])].copy()

    # Define two enrichment comparisons:
    # (A) Primary Non-Bottleneck 5 Non-NFE (afr, amr, eas, sas, mid) vs NFE (AF >= 1e-4, >= 5x enrichment)
    # (B) All-7 Non-NFE (incl asj, fin) vs NFE (AF >= 1e-4, >= 5x enrichment)
    out = {"n_coding_splice_joined": int(len(sub))}

    for mode_name, nn_af_col, nn_pop_col in [
        ("nonbottleneck5_vs_nfe", "af_nonnfe5", "pop_nonnfe5"),
        ("all7_nonnfe_vs_nfe", "af_nonnfe7", "pop_nonnfe7"),
    ]:
        is_nn_enr = (sub[nn_af_col] >= 1e-4) & (sub[nn_af_col] / np.maximum(sub["af_nfe"], 1e-6) >= 5.0)
        is_eur_enr = (sub["af_nfe"] >= 1e-4) & (sub["af_nfe"] / np.maximum(sub[nn_af_col], 1e-6) >= 5.0)
        comp = sub[is_nn_enr | is_eur_enr].copy()
        comp["non_nfe_enr"] = is_nn_enr.loc[comp.index].astype(int)
        comp["af_ref"] = np.where(comp["non_nfe_enr"] == 1, comp[nn_af_col], comp["af_nfe"])
        comp["log10_af"] = np.log10(comp["af_ref"].clip(lower=1e-6))
        comp["af_bin"] = pd.cut(
            comp["af_ref"],
            bins=[1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 5e-2, 1.01],
            labels=["1e-4_3e-4", "3e-4_1e-3", "1e-3_3e-3", "3e-3_1e-2", "1e-2_5e-2", "ge_5e-2"],
            right=False
        )
        comp["strata_af_csq"] = comp["af_bin"].astype(str) + "_" + comp["csq_group"].astype(str)

        res_m = {
            "n_non_nfe": int((comp["non_nfe_enr"] == 1).sum()),
            "n_nfe": int((comp["non_nfe_enr"] == 0).sum()),
            "median_af_non_nfe": float(comp.loc[comp["non_nfe_enr"] == 1, "af_ref"].median()),
            "median_af_nfe": float(comp.loc[comp["non_nfe_enr"] == 0, "af_ref"].median()),
            "outcomes": {},
        }
        for label, target_cats in [
            ("VUS", ["VUS"]),
            ("VUS_or_CONFLICTING", ["VUS", "CONFLICTING"]),
            ("B_LB", ["B_LB"]),
            ("P_LP", ["P_LP"]),
            ("CONFLICTING", ["CONFLICTING"]),
        ]:
            comp["y"] = comp["clnsig_cat"].isin(target_cats).astype(int)
            r_nn = float(comp.loc[comp["non_nfe_enr"] == 1, "y"].mean())
            r_eu = float(comp.loc[comp["non_nfe_enr"] == 0, "y"].mean())
            # Unadjusted Fisher exact
            a = int(((comp["non_nfe_enr"] == 1) & (comp["y"] == 1)).sum())
            b = int(((comp["non_nfe_enr"] == 1) & (comp["y"] == 0)).sum())
            c = int(((comp["non_nfe_enr"] == 0) & (comp["y"] == 1)).sum())
            d = int(((comp["non_nfe_enr"] == 0) & (comp["y"] == 0)).sum())
            f_or, f_p = stats.fisher_exact([[a, b], [c, d]])
            # Frequency-matched CMH & Frequency+Consequence matched CMH
            cmh_af = cmh_test(comp, "non_nfe_enr", "y", "af_bin")
            cmh_af_csq = cmh_test(comp, "non_nfe_enr", "y", "strata_af_csq")
            # Multivariable Logistic Regression (Gene FE + Consequence FE + quadratic log10 AF)
            try:
                mlog = smf.logit("y ~ non_nfe_enr + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=comp).fit(disp=False, maxiter=100)
                b_l = float(mlog.params["non_nfe_enr"])
                se_l = float(mlog.bse["non_nfe_enr"])
                logit_res = {
                    "or": float(np.exp(b_l)),
                    "ci95": [float(np.exp(b_l - 1.96 * se_l)), float(np.exp(b_l + 1.96 * se_l))],
                    "p": float(mlog.pvalues["non_nfe_enr"]),
                }
            except Exception as e:
                # Fallback without gene FE if rare separation (e.g., P_LP)
                mlog = smf.logit("y ~ non_nfe_enr + log10_af + I(log10_af**2) + C(csq_group) + C(domain)", data=comp).fit(disp=False, maxiter=100)
                b_l = float(mlog.params["non_nfe_enr"])
                se_l = float(mlog.bse["non_nfe_enr"])
                logit_res = {
                    "or": float(np.exp(b_l)),
                    "ci95": [float(np.exp(b_l - 1.96 * se_l)), float(np.exp(b_l + 1.96 * se_l))],
                    "p": float(mlog.pvalues["non_nfe_enr"]),
                    "note": "domain_fe_fallback",
                }
            res_m["outcomes"][label] = {
                "rate_non_nfe": r_nn, "rate_nfe": r_eu,
                "counts_2x2": {"nn_pos": a, "nn_neg": b, "nfe_pos": c, "nfe_neg": d},
                "fisher_or": float(f_or), "fisher_p": float(f_p),
                "cmh_af_matched": cmh_af,
                "cmh_af_csq_matched": cmh_af_csq,
                "logit_gene_csq_af_adj": logit_res,
            }

        # Also per-AF-bin breakdown for Figure!
        bin_rows = []
        for bname, bg in comp.groupby("af_bin", observed=True):
            for grp_val, gname in [(0, "NFE_enriched"), (1, "NonNFE_enriched")]:
                gg = bg[bg["non_nfe_enr"] == grp_val]
                bin_rows.append({
                    "af_bin": str(bname), "group": gname, "n": int(len(gg)),
                    "vus_rate": float((gg["clnsig_cat"] == "VUS").mean()) if len(gg) else np.nan,
                    "blb_rate": float((gg["clnsig_cat"] == "B_LB").mean()) if len(gg) else np.nan,
                    "conf_rate": float((gg["clnsig_cat"] == "CONFLICTING").mean()) if len(gg) else np.nan,
                    "plp_rate": float((gg["clnsig_cat"] == "P_LP").mean()) if len(gg) else np.nan,
                })
        res_m["by_af_bin"] = bin_rows
        out[mode_name] = res_m

    # Ancestry-by-Ancestry Profile (each ancestry p enriched >=5x vs NFE & AF_p >= 1e-4, compared to NFE-enriched!)
    eur_base = sub[(sub["af_nfe"] >= 1e-4) & (sub["af_nfe"] / np.maximum(sub["af_nonnfe5"], 1e-6) >= 5.0)].copy()
    eur_base["is_target"] = 0
    eur_base["af_ref"] = eur_base["af_nfe"]
    anc_profile = {}
    for p in POPS_ALL8:
        if p == "nfe":
            anc_profile["nfe"] = {
                "n": int(len(eur_base)),
                "vus_rate": float((eur_base["clnsig_cat"] == "VUS").mean()),
                "blb_rate": float((eur_base["clnsig_cat"] == "B_LB").mean()),
                "conf_rate": float((eur_base["clnsig_cat"] == "CONFLICTING").mean()),
                "plp_rate": float((eur_base["clnsig_cat"] == "P_LP").mean()),
                "cmh_vus_or_vs_nfe": 1.0, "cmh_vus_ci95": [1.0, 1.0], "cmh_vus_p": 1.0,
                "cmh_blb_or_vs_nfe": 1.0, "cmh_blb_ci95": [1.0, 1.0], "cmh_blb_p": 1.0,
            }
            continue
        psub = sub[(sub[f"af_{p}"] >= 1e-4) & (sub[f"af_{p}"] / np.maximum(sub["af_nfe"], 1e-6) >= 5.0) & (sub["pop_max8"] == p)].copy()
        psub["is_target"] = 1
        psub["af_ref"] = psub[f"af_{p}"]
        comb = pd.concat([eur_base, psub], ignore_index=True)
        comb["af_bin"] = pd.cut(comb["af_ref"], bins=[1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 5e-2, 1.01], right=False)
        comb["y_vus"] = (comb["clnsig_cat"] == "VUS").astype(int)
        comb["y_blb"] = (comb["clnsig_cat"] == "B_LB").astype(int)
        comb["y_plp"] = (comb["clnsig_cat"] == "P_LP").astype(int)
        c_vus = cmh_test(comb, "is_target", "y_vus", "af_bin")
        c_blb = cmh_test(comb, "is_target", "y_blb", "af_bin")
        c_plp = cmh_test(comb, "is_target", "y_plp", "af_bin")
        anc_profile[p] = {
            "n": int(len(psub)),
            "vus_rate": float((psub["clnsig_cat"] == "VUS").mean()) if len(psub) else np.nan,
            "blb_rate": float((psub["clnsig_cat"] == "B_LB").mean()) if len(psub) else np.nan,
            "conf_rate": float((psub["clnsig_cat"] == "CONFLICTING").mean()) if len(psub) else np.nan,
            "plp_rate": float((psub["clnsig_cat"] == "P_LP").mean()) if len(psub) else np.nan,
            "cmh_vus_or_vs_nfe": c_vus["or"], "cmh_vus_ci95": c_vus["ci95"], "cmh_vus_p": c_vus["p"],
            "cmh_blb_or_vs_nfe": c_blb["or"], "cmh_blb_ci95": c_blb["ci95"], "cmh_blb_p": c_blb["p"],
            "cmh_plp_or_vs_nfe": c_plp["or"], "cmh_plp_ci95": c_plp["ci95"], "cmh_plp_p": c_plp["p"],
        }
    out["by_ancestry_profile"] = anc_profile
    return out


def run_h2_h3_reclass_debt(joined):
    """H2 & H3: Exact High-Confidence (>=2-star) & 1-star Reclassification Debt, Carrier Burden & Staleness Audit."""
    # Focus on non-conflicting, non-low-penetrance P/LP variants
    plp = joined[joined["clnsig_cat"] == "P_LP"].copy()
    plp["is_founder"] = plp.apply(is_known_founder, axis=1)
    plp["eval_year"] = pd.to_numeric(plp["last_evaluated"].astype(str).str.extract(r"(\d{4})")[0], errors="coerce")
    plp["post_gnomad_v4_eval"] = (plp["eval_year"] >= 2024).astype(int)

    ad_plp = plp[plp["inheritance"] == "AD"].copy()
    ad_hc = ad_plp[ad_plp["stars"] >= 2].copy()  # 2-star, 3-star (Expert Panel), 4-star!
    ad_1s = ad_plp[ad_plp["stars"] == 1].copy()

    def summarize_debt(df_d, label):
        t1_g5 = df_d[df_d["faf95_grpmax5"] >= 1e-4]
        t1_a8 = df_d[df_d["faf95_max8"] >= 1e-4]
        t2_g5 = df_d[df_d["faf95_grpmax5"] >= 5e-4]
        t3_g5 = df_d[(df_d["faf95_grpmax5"] >= 1e-3) | (df_d["af_grpmax5"] >= 5e-3) | (df_d["hom_joint"] >= 2)]
        # NFE invisibility among Tier 1 (all8 and grpmax5)
        nfe_invis_g5 = t1_g5[t1_g5["faf95_nfe"] < 1e-4]
        nfe_invis_a8 = t1_a8[t1_a8["faf95_nfe"] < 1e-4]
        # Founder-excluded Tier 1
        t1_g5_nf = t1_g5[~t1_g5["is_founder"]]
        t1_a8_nf = t1_a8[~t1_a8["is_founder"]]
        nfe_invis_g5_nf = t1_g5_nf[t1_g5_nf["faf95_nfe"] < 1e-4]
        nfe_invis_a8_nf = t1_a8_nf[t1_a8_nf["faf95_nfe"] < 1e-4]
        return {
            "label": label,
            "n_total_plp_in_gnomad": int(len(df_d)),
            "tier1_faf95_grpmax5_ge_1e4": int(len(t1_g5)),
            "tier1_faf95_grpmax5_founder_excl": int(len(t1_g5_nf)),
            "tier1_faf95_all8_ge_1e4": int(len(t1_a8)),
            "tier1_faf95_all8_founder_excl": int(len(t1_a8_nf)),
            "tier2_faf95_grpmax5_ge_5e4": int(len(t2_g5)),
            "tier3_extreme_or_hom2": int(len(t3_g5)),
            "nfe_invisible_grpmax5": int(len(nfe_invis_g5)),
            "nfe_invisible_grpmax5_share": float(len(nfe_invis_g5) / len(t1_g5)) if len(t1_g5) else 0.0,
            "nfe_invisible_grpmax5_nf_share": float(len(nfe_invis_g5_nf) / len(t1_g5_nf)) if len(t1_g5_nf) else 0.0,
            "nfe_invisible_all8": int(len(nfe_invis_a8)),
            "nfe_invisible_all8_share": float(len(nfe_invis_a8) / len(t1_a8)) if len(t1_a8) else 0.0,
            "nfe_invisible_all8_nf_share": float(len(nfe_invis_a8_nf) / len(t1_a8_nf)) if len(t1_a8_nf) else 0.0,
            "stars_breakdown_tier1_all8": {str(k): int(v) for k, v in t1_a8["stars"].value_counts().sort_index().items()},
            "domain_breakdown_tier1_all8": {str(k): int(v) for k, v in t1_a8["domain"].value_counts().items()},
            "domain_breakdown_tier1_grpmax5": {str(k): int(v) for k, v in t1_g5["domain"].value_counts().items()},
            "pop_driver_tier1_all8": {str(k): int(v) for k, v in t1_a8["faf95_pop8"].value_counts().items()},
            "pop_driver_tier1_grpmax5": {str(k): int(v) for k, v in t1_g5["faf95_pop5"].value_counts().items()},
            "post_2024_eval_count_tier1_all8": int(t1_a8["post_gnomad_v4_eval"].sum()),
            "post_2024_eval_share_tier1_all8": float(t1_a8["post_gnomad_v4_eval"].mean()) if len(t1_a8) else 0.0,
            "median_eval_year_flagged": float(t1_a8["eval_year"].median()) if len(t1_a8) else np.nan,
            "median_eval_year_unflagged": float(df_d.loc[df_d["faf95_max8"] < 1e-4, "eval_year"].median()) if len(df_d) else np.nan,
        }

    out = {
        "ad_72genes_highconf_2to4star": summarize_debt(ad_hc, "AD 72 genes, 2-4 star P/LP"),
        "ad_72genes_1star": summarize_debt(ad_1s, "AD 72 genes, 1-star P/LP"),
        "ad_72genes_all_stars": summarize_debt(ad_plp, "AD 72 genes, all stars P/LP"),
        "ar_9genes_highconf_2to4star_ge_1e3": int(len(plp[(plp["inheritance"] == "AR") & (plp["stars"] >= 2) & (plp["faf95_max8"] >= 1e-3)])),
        "xl_3genes_highconf_2to4star_ge_1e4": int(len(plp[(plp["inheritance"] == "XL") & (plp["stars"] >= 2) & (plp["faf95_max8"] >= 1e-4)])),
    }

    # H3 Staleness statistical test (Mann-Whitney & Fisher exact on pre-2024 vs >=2024 eval year)
    flagged_hc = ad_hc[ad_hc["faf95_max8"] >= 1e-4]
    unflagged_hc = ad_hc[ad_hc["faf95_max8"] < 1e-4]
    ev_f = flagged_hc["eval_year"].dropna()
    ev_u = unflagged_hc["eval_year"].dropna()
    mw_stat, mw_p = stats.mannwhitneyu(ev_f, ev_u, alternative="two-sided") if len(ev_f) > 2 else (np.nan, np.nan)
    out["h3_staleness_test_ad_hc"] = {
        "flagged_n_with_date": int(len(ev_f)),
        "unflagged_n_with_date": int(len(ev_u)),
        "flagged_mean_year": float(ev_f.mean()) if len(ev_f) else np.nan,
        "unflagged_mean_year": float(ev_u.mean()) if len(ev_u) else np.nan,
        "flagged_pre2024_share": float((ev_f < 2024).mean()) if len(ev_f) else np.nan,
        "unflagged_pre2024_share": float((ev_u < 2024).mean()) if len(ev_u) else np.nan,
        "flagged_3star_expert_panel_total": int((flagged_hc["stars"] == 3).sum()),
        "flagged_3star_evaluated_2024_2026": int(((flagged_hc["stars"] == 3) & (flagged_hc["eval_year"] >= 2024)).sum()),
        "flagged_2star_evaluated_2024_2026": int(((flagged_hc["stars"] == 2) & (flagged_hc["eval_year"] >= 2024)).sum()),
        "mann_whitney_p": float(mw_p),
    }

    # Individual-level Cumulative Carrier Burden by Ancestry (exact Hardy-Weinberg probability)
    # Computed on AD >=2-star Tier 1 (FAF95_max8 >= 1e-4): (a) ALL flagged, (b) Founder-Excluded, (c) NFE-Invisible only!
    carrier_burden = {}
    for subset_name, sdf in [
        ("all_flagged_ad_hc", flagged_hc),
        ("founder_excluded_ad_hc", flagged_hc[~flagged_hc["is_founder"]]),
        ("nfe_invisible_founder_excl_ad_hc", flagged_hc[(~flagged_hc["is_founder"]) & (flagged_hc["faf95_nfe"] < 1e-4)]),
        ("all_flagged_ad_allstars", ad_plp[(ad_plp["faf95_max8"] >= 1e-4) & (~ad_plp["is_founder"])]),
    ]:
        cb = {}
        for p in POPS_ALL8:
            afs = sdf[f"af_{p}"].values.clip(0.0, 1.0)
            p_car = 1.0 - np.prod(1.0 - (2.0 * afs * (1.0 - afs) + afs * afs))
            cb[p] = {
                "carrier_prob": float(p_car),
                "per_1000_individuals": round(float(p_car * 1000.0), 3),
                "expected_variants_per_1000": round(float((2.0 * afs).sum() * 1000.0), 3),
            }
        # Ratio vs NFE
        nfe_c = cb["nfe"]["carrier_prob"]
        for p in POPS_ALL8:
            cb[p]["ratio_vs_nfe"] = round(cb[p]["carrier_prob"] / nfe_c, 2) if nfe_c > 0 else np.nan
        carrier_burden[subset_name] = cb
    out["carrier_burden_by_ancestry"] = carrier_burden

    # Top gene-level breakdown of AD >=2-star Tier 1 debt
    gene_debt = (
        flagged_hc.groupby(["gene", "domain"], as_index=False)
        .agg(
            n_flagged_hc=("variant_id", "count"),
            n_founder_excl=("is_founder", lambda s: int((~s).sum())),
            n_3star=("stars", lambda s: int((s == 3).sum())),
            n_nfe_invis=("faf95_nfe", lambda s: int((s < 1e-4).sum())),
            max_faf95=("faf95_max8", "max"),
        )
        .sort_values(["n_founder_excl", "n_flagged_hc"], ascending=False)
    )
    out["top_genes_ad_hc_debt"] = gene_debt.to_dict(orient="records")

    # Export complete allele-exact catalog of ALL AD >=2-star Tier 1 (FAF95_max8 >= 1e-4) variants!
    cols_keep = [
        "variant_id", "chrom", "pos", "ref", "alt", "variation_id", "gene", "domain", "inheritance",
        "hgvs_name", "clnsig", "clnrevstat", "stars", "last_evaluated", "n_submitters", "csq_group",
        "is_founder", "ac_joint", "an_joint", "hom_joint", "af_joint",
        "af_nfe", "faf95_nfe", "af_grpmax5", "pop_grpmax5", "faf95_grpmax5", "faf95_pop5",
        "af_max8", "pop_max8", "faf95_max8", "faf95_pop8",
        "af_afr", "af_amr", "af_eas", "af_sas", "af_mid", "af_asj", "af_fin",
    ]
    cat_df = flagged_hc[cols_keep].sort_values(["is_founder", "stars", "faf95_max8"], ascending=[True, False, False])
    cat_df.to_csv(os.path.join(RES, "a1_ad_highconf_reclassification_catalog.csv"), index=False)
    out["top30_founder_excluded_ad_hc"] = cat_df[~cat_df["is_founder"]].head(30).to_dict(orient="records")
    out["known_founders_detected"] = cat_df[cat_df["is_founder"]].to_dict(orient="records")
    return out, cat_df


def run_h4_rsid_collision_audit(cv_df, g_df):
    """
    H4 Methodological Audit: Quantify exact false-positive inflation when joining ClinVar P/LP
    variants to gnomAD v4.1 by rsID vs exact GRCh38 (CHROM-POS-REF-ALT) across ALL 84 ACMG SF v3.3 genes!
    """
    plp = cv_df[cv_df["clnsig_cat"] == "P_LP"].copy()
    plp_rs = plp[plp["rs"].str.startswith("rs", na=False)].copy()

    # Build gnomAD rsID -> list of (variant_id, faf95_max8, af_max8, pop_max8)
    rs_map = {}
    for row in g_df[["variant_id", "gnomad_rsids", "faf95_max8", "af_max8", "faf95_grpmax5", "af_nfe"]].itertuples(index=False):
        if not row.gnomad_rsids:
            continue
        for rs in row.gnomad_rsids.split("|"):
            if rs:
                rs_map.setdefault(rs, []).append((row.variant_id, row.faf95_max8, row.af_max8, row.faf95_grpmax5, row.af_nfe))

    exact_lookup = set(g_df["variant_id"])
    exact_faf = dict(zip(g_df["variant_id"], g_df["faf95_max8"]))
    exact_af = dict(zip(g_df["variant_id"], g_df["af_max8"]))

    naive_flagged_af1e4 = 0
    exact_confirmed_af1e4 = 0
    false_pos_collision_af1e4 = 0
    naive_flagged_faf1e4 = 0
    exact_confirmed_faf1e4 = 0
    false_pos_collision_faf1e4 = 0
    naive_hc_faf1e4 = 0
    exact_hc_faf1e4 = 0
    false_hc_faf1e4 = 0

    for row in plp_rs[["variant_id", "rs", "stars", "inheritance"]].itertuples(index=False):
        hits = rs_map.get(row.rs)
        if not hits:
            continue
        max_naive_af = max(h[2] for h in hits)
        max_naive_faf = max(h[1] for h in hits)
        ex_af = exact_af.get(row.variant_id, 0.0)
        ex_faf = exact_faf.get(row.variant_id, 0.0)
        if max_naive_af >= 1e-4:
            naive_flagged_af1e4 += 1
            if ex_af >= 1e-4:
                exact_confirmed_af1e4 += 1
            else:
                false_pos_collision_af1e4 += 1
        if max_naive_faf >= 1e-4:
            naive_flagged_faf1e4 += 1
            if ex_faf >= 1e-4:
                exact_confirmed_faf1e4 += 1
            else:
                false_pos_collision_faf1e4 += 1
            if row.stars >= 2 and row.inheritance == "AD":
                naive_hc_faf1e4 += 1
                if ex_faf >= 1e-4:
                    exact_hc_faf1e4 += 1
                else:
                    false_hc_faf1e4 += 1

    return {
        "n_clinvar_plp_with_rsid": int(len(plp_rs)),
        "af_ge_1e4": {
            "naive_rsid_flagged": naive_flagged_af1e4,
            "allele_exact_confirmed": exact_confirmed_af1e4,
            "multiallelic_false_positives": false_pos_collision_af1e4,
            "false_positive_rate": float(false_pos_collision_af1e4 / naive_flagged_af1e4) if naive_flagged_af1e4 else 0.0,
        },
        "faf95_ge_1e4": {
            "naive_rsid_flagged": naive_flagged_faf1e4,
            "allele_exact_confirmed": exact_confirmed_faf1e4,
            "multiallelic_false_positives": false_pos_collision_faf1e4,
            "false_positive_rate": float(false_pos_collision_faf1e4 / naive_flagged_faf1e4) if naive_flagged_faf1e4 else 0.0,
        },
        "ad_highconf_2to4star_faf95_ge_1e4": {
            "naive_rsid_flagged": naive_hc_faf1e4,
            "allele_exact_confirmed": exact_hc_faf1e4,
            "multiallelic_false_positives": false_hc_faf1e4,
            "false_positive_rate": float(false_hc_faf1e4 / naive_hc_faf1e4) if naive_hc_faf1e4 else 0.0,
        },
    }


def spot_validate_faf95_api(cat_df, n_check=10):
    """Cross-validate our computed FAF95 against gnomAD's native joint.faf95 via get_variant_frequency.py."""
    nf = cat_df[~cat_df["is_founder"]].head(n_check)
    val_dir = os.path.join(RES, "faf95_api_validation")
    os.makedirs(val_dir, exist_ok=True)
    checks = []
    for row in nf.itertuples(index=False):
        vid = row.variant_id
        out_f = os.path.join(val_dir, f"{vid}.json")
        if not os.path.exists(out_f):
            p = subprocess.run(
                ["uv", "run", GNOMAD_FREQ_SCRIPT, "--variant_id", vid, "--dataset", "gnomad_r4", "--output", out_f],
                capture_output=True, text=True, timeout=60
            )
            time.sleep(6.2)
            if p.returncode != 0 or not os.path.exists(out_f):
                continue
        d = json.load(open(out_f))
        v = (d.get("data") or {}).get("variant") or {}
        jf = (v.get("joint") or {}).get("faf95") or {}
        api_faf = jf.get("popmax")
        api_pop = jf.get("popmax_population")
        checks.append({
            "variant_id": vid, "gene": row.gene, "stars": int(row.stars),
            "computed_faf95_grpmax5": round(float(row.faf95_grpmax5), 7),
            "computed_pop5": row.faf95_pop5,
            "gnomad_api_joint_faf95": round(float(api_faf), 7) if api_faf is not None else None,
            "gnomad_api_joint_pop": api_pop,
        })
    valid_pairs = [(c["computed_faf95_grpmax5"], c["gnomad_api_joint_faf95"]) for c in checks if c["gnomad_api_joint_faf95"] is not None]
    r_val = float(stats.pearsonr([x[0] for x in valid_pairs], [x[1] for x in valid_pairs])[0]) if len(valid_pairs) >= 3 else np.nan
    return {"n_validated": len(valid_pairs), "pearson_r": r_val, "details": checks}


def main():
    t0 = time.time()
    cv_df = pd.read_parquet(os.path.join(CLEAN, "clinvar_acmg84_grch38.parquet"))
    cv_df["clnsig_cat"] = [classify_clnsig(a, b) for a, b in zip(cv_df["clnsig"], cv_df["clnsigconf"])]
    g_df = load_gnomad_acmg84()

    # Exact GRCh38 (CHROM-POS-REF-ALT) inner join
    joined = cv_df.merge(g_df, on="variant_id", how="inner")
    joined["csq_group"] = [classify_consequence(a, b) for a, b in zip(joined["mc"], joined["gnomad_csq"])]
    joined.to_parquet(os.path.join(CLEAN, "a1_clinvar_gnomad_exact_joined.parquet"), index=False)

    print(f"ClinVar ACMG84: {len(cv_df):,} | gnomAD v4.1 ACMG84: {len(g_df):,} | Allele-Exact Joined: {len(joined):,} ({joined['gene'].nunique()} genes)", flush=True)

    print("Running H1 Frequency-Matched Cross-Ancestry Equity analysis...", flush=True)
    h1 = run_h1_equity(joined)

    print("Running H2 & H3 High-Confidence Reclassification Debt & Staleness Audit...", flush=True)
    h2_h3, cat_df = run_h2_h3_reclass_debt(joined)

    print("Running H4 Multiallelic rsID Collision Audit...", flush=True)
    h4 = run_h4_rsid_collision_audit(cv_df, g_df)

    print("Spot-validating top 10 AD >=2-star hits against gnomAD native joint.faf95 API...", flush=True)
    api_val = spot_validate_faf95_api(cat_df, n_check=10)

    master = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_sec": round(time.time() - t0, 1),
        "overview": {
            "n_acmg_genes_panel": 84,
            "n_clinvar_acmg84_variants": int(len(cv_df)),
            "n_clinvar_genes_present": int(cv_df["gene"].nunique()),
            "n_gnomad_v4_acmg84_variants": int(len(g_df)),
            "n_gnomad_genes_present": int(g_df["gnomad_gene"].nunique()),
            "n_allele_exact_joined": int(len(joined)),
            "joined_clnsig_counts": {str(k): int(v) for k, v in joined["clnsig_cat"].value_counts().items()},
            "joined_domain_counts": {str(k): int(v) for k, v in joined["domain"].value_counts().items()},
        },
        "h1_equity": h1,
        "h2_h3_reclassification_debt": h2_h3,
        "h4_rsid_collision_audit": h4,
        "gnomad_api_faf95_validation": api_val,
    }
    out_json = os.path.join(RES, "a1_reclassification_results.json")
    with open(out_json, "w") as f:
        json.dump(master, f, indent=2)
    print(f"DONE in {master['elapsed_sec']}s! Saved {out_json}", flush=True)

    # Crisp headline summary
    print("\n=== STUDY A-1 HEADLINE RESULTS ===")
    nb = h1["nonbottleneck5_vs_nfe"]
    print(f"H1 (Non-Bottleneck 5 Non-NFE n={nb['n_non_nfe']} vs NFE n={nb['n_nfe']}):")
    for o in ("VUS", "B_LB", "P_LP", "CONFLICTING"):
        r = nb["outcomes"][o]
        print(f"  {o:12s}: NonNFE={r['rate_non_nfe']*100:5.2f}% vs NFE={r['rate_nfe']*100:5.2f}% | "
              f"CMH_AF_OR={r['cmh_af_matched']['or']:.2f} (p={r['cmh_af_matched']['p']:.2e}) | "
              f"Logit_GeneFE_OR={r['logit_gene_csq_af_adj']['or']:.2f} (p={r['logit_gene_csq_af_adj']['p']:.2e})")
    hc = h2_h3["ad_72genes_highconf_2to4star"]
    print(f"\nH2 (72 AD genes, >=2-star P/LP total={hc['n_total_plp_in_gnomad']}):")
    print(f"  Tier1 FAF95_grpmax5>=1e-4: {hc['tier1_faf95_grpmax5_ge_1e4']} (founder-excl: {hc['tier1_faf95_grpmax5_founder_excl']}), NFE-invisible: {hc['nfe_invisible_grpmax5']} ({hc['nfe_invisible_grpmax5_share']*100:.1f}%)")
    print(f"  Tier1 FAF95_all8>=1e-4:    {hc['tier1_faf95_all8_ge_1e4']} (founder-excl: {hc['tier1_faf95_all8_founder_excl']}), NFE-invisible: {hc['nfe_invisible_all8']} ({hc['nfe_invisible_all8_share']*100:.1f}%)")
    print(f"  Tier2 FAF95_grpmax5>=5e-4: {hc['tier2_faf95_grpmax5_ge_5e4']} | Tier3 Extreme/Hom>=2: {hc['tier3_extreme_or_hom2']}")
    print(f"  Post-2024 reviewed among flagged >=2-star: {hc['post_2024_eval_count_tier1_all8']} ({hc['post_2024_eval_share_tier1_all8']*100:.1f}%)")
    print(f"\nH4 (rsID Collision Audit, AD >=2-star FAF95>=1e-4): naive={h4['ad_highconf_2to4star_faf95_ge_1e4']['naive_rsid_flagged']}, exact={h4['ad_highconf_2to4star_faf95_ge_1e4']['allele_exact_confirmed']}, FP rate={h4['ad_highconf_2to4star_faf95_ge_1e4']['false_positive_rate']*100:.1f}%")
    print(f"API FAF95 spot validation: n={api_val['n_validated']}, Pearson r={api_val['pearson_r']:.5f}")


if __name__ == "__main__":
    main()
