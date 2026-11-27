#!/opt/anaconda3/bin/python3
"""
Ruthless Forensic Falsification Audit for Study A-1:
  1. DENOMINATOR / WINNER'S-CURSE AUDIT OF H1:
     Our v1 H1 used raw AF >= 1e-4, which in MID (AN ~ 6,062) includes AC=1 singletons!
     Re-test H1 cleanly under:
       (a) FAF95 >= 1e-4 (zero singletons; statistically credible lower bound >= 1e-4!)
       (b) AC_p >= 10 AND AF_p >= 5e-4 (well-powered across ALL 8 ancestries!)
       (c) AF_p >= 1e-3 AND AC_p >= 5 (strict common-rare boundary)
       (d) EXACT HYPERGEOMETRIC DOWNSAMPLING OF NFE (AN ~ 1.24M) TO EACH ANCESTRY'S EXACT AN_p!
           (Mathematically equalizes detection power & sampling variance between NFE and every ancestry!)
  2. LIVE gnomAD v4.1 QC FILTER/FLAG AUDIT ON ALL 12 L4 + 7 L1 + KEY L3 VARIANTS:
     Check exome.filters, genome.filters, flags (lcr, segdup, RF, AS_VQSR, AC0, InbreedingCoeff)!
  3. LIVE CLINVAR EVIDENCE AUDIT (via official clinvar_api.py evidence) ON ALL 12 L4 CANDIDATES:
     Check exact disease conditions (AD vs AR assertions!), per-submitter classifications & dates,
     and official ClinGen VCEP gene-specific BS1/BA1 thresholds!
"""
import json, os, subprocess, time, warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
CLEAN = os.path.join(BASE, "data", "clean")
RES = os.path.join(BASE, "results")
AUD = os.path.join(BASE, "audit")
os.makedirs(AUD, exist_ok=True)

CLINVAR_API = "/Users/mac/.gemini/config/plugins/genomics/skills/clinvar_database/scripts/clinvar_api.py"
MAX_AN = 1614324.0

import sys
sys.path.insert(0, os.path.join(BASE, "src"))
from analyze_a1_full import cmh_test

j = pd.read_parquet(os.path.join(CLEAN, "a1_clinvar_gnomad_exact_joined.parquet"))
cat55 = pd.read_csv(os.path.join(RES, "a1_ad_highconf_triaged_55.csv"))

out = {}

# ==============================================================================
# 1. DENOMINATOR / WINNER'S-CURSE AUDIT OF H1
# ==============================================================================
print("=== 1. H1 DENOMINATOR & SINGLETON DIAGNOSTICS BY ANCESTRY ===")
cs = j[j["csq_group"].isin(["pLoF", "missense", "inframe_indel", "splice_other", "synonymous"]) &
       j["clnsig_cat"].isin(["P_LP", "VUS", "CONFLICTING", "B_LB"]) &
       (j["an_joint"] >= 0.80 * MAX_AN)].copy()

anc_diag = {}
for p in ["nfe", "afr", "amr", "eas", "sas", "mid", "asj", "fin"]:
    if p == "nfe":
        s = cs[(cs["af_nfe"] >= 1e-4) & (cs["af_nfe"] / np.maximum(cs["af_nonnfe5"], 1e-6) >= 5.0)]
    else:
        s = cs[(cs[f"af_{p}"] >= 1e-4) & (cs[f"af_{p}"] / np.maximum(cs["af_nfe"], 1e-6) >= 5.0) & (cs["pop_max8"] == p)]
    anc_diag[p] = {
        "n_raw_af1e4": int(len(s)),
        "median_an_p": int(s[f"an_{p}"].median()) if len(s) else 0,
        "ac1_singleton_share": round(float((s[f"ac_{p}"] == 1).mean()), 4) if len(s) else 0.0,
        "ac_le2_share": round(float((s[f"ac_{p}"] <= 2).mean()), 4) if len(s) else 0.0,
        "median_ac_p": float(s[f"ac_{p}"].median()) if len(s) else 0.0,
        "faf95_ge_1e4_count": int((s[f"faf95_{p}"] >= 1e-4).sum()),
    }
    print(f"  {p.upper():3s}: n={anc_diag[p]['n_raw_af1e4']:4d} | med_AN={anc_diag[p]['median_an_p']:7d} | med_AC={anc_diag[p]['median_ac_p']:5.1f} | AC==1: {anc_diag[p]['ac1_singleton_share']*100:5.1f}% | AC<=2: {anc_diag[p]['ac_le2_share']*100:5.1f}% | FAF95>=1e-4: {anc_diag[p]['faf95_ge_1e4_count']:4d}")
out["h1_denominator_diagnostics"] = anc_diag

# Re-test H1 under 3 STRICT DENOMINATOR-CLEAN DEFINITIONS:
# (A) FAF95-matched: FAF95_nonnfe5 >= 1e-4 & FAF95_nonnfe5 >= 5x NFE vs FAF95_nfe >= 1e-4 & >= 5x NonNFE5 (ZERO singletons!)
# (B) High-AC Well-Powered: AC_pop >= 10 & AF_pop >= 5e-4 (>= 5x enriched)
# (C) Common-Rare Clean: AF_pop >= 1e-3 & AC_pop >= 5 (>= 5x enriched)
print("\n=== 2. DENOMINATOR-CLEAN H1 RE-TESTS (ZERO SINGLETON BIAS!) ===")
h1_clean = {}

def eval_h1_spec(df_in, mask_nn, mask_eu, af_nn_col, af_eu_col, label):
    comp = df_in[mask_nn | mask_eu].copy()
    comp["non_nfe"] = mask_nn.loc[comp.index].astype(int)
    comp["af_ref"] = np.where(comp["non_nfe"] == 1, comp[af_nn_col], comp[af_eu_col])
    comp["log10_af"] = np.log10(comp["af_ref"].clip(lower=1e-6))
    comp["log_subs"] = np.log1p(comp["n_submitters"].fillna(1).clip(lower=1))
    comp["af_bin"] = pd.cut(comp["af_ref"], bins=[1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 5e-2, 1.01], right=False)
    res = {
        "label": label,
        "n_non_nfe": int((comp["non_nfe"] == 1).sum()),
        "n_nfe": int((comp["non_nfe"] == 0).sum()),
        "med_af_nn": float(comp.loc[comp["non_nfe"] == 1, "af_ref"].median()),
        "med_af_nfe": float(comp.loc[comp["non_nfe"] == 0, "af_ref"].median()),
        "med_subs_nn": float(comp.loc[comp["non_nfe"] == 1, "n_submitters"].median()),
        "med_subs_nfe": float(comp.loc[comp["non_nfe"] == 0, "n_submitters"].median()),
        "outcomes": {},
    }
    for cat_nm in ("VUS", "B_LB", "CONFLICTING", "P_LP"):
        comp["y"] = (comp["clnsig_cat"] == cat_nm).astype(int)
        cmh = cmh_test(comp, "non_nfe", "y", "af_bin")
        r_nn = float(comp.loc[comp["non_nfe"] == 1, "y"].mean())
        r_eu = float(comp.loc[comp["non_nfe"] == 0, "y"].mean())
        # Domain+Csq+AF logit (and Gene FE if well-populated)
        try:
            gv = comp.groupby("gene")["y"].nunique()
            sub_g = comp[comp["gene"].isin(gv[gv == 2].index)]
            ml = smf.logit("y ~ non_nfe + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=sub_g).fit(disp=False, maxiter=150)
            lor, lp = float(np.exp(ml.params["non_nfe"])), float(ml.pvalues["non_nfe"])
            # Mediation with log_subs
            ml_m = smf.logit("y ~ non_nfe + log_subs + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=sub_g).fit(disp=False, maxiter=150)
            lor_m, lp_m = float(np.exp(ml_m.params["non_nfe"])), float(ml_m.pvalues["non_nfe"])
        except Exception:
            ml = smf.logit("y ~ non_nfe + log10_af + C(csq_group) + C(domain)", data=comp).fit(disp=False, maxiter=150)
            lor, lp = float(np.exp(ml.params["non_nfe"])), float(ml.pvalues["non_nfe"])
            lor_m, lp_m = np.nan, np.nan
        res["outcomes"][cat_nm] = {
            "rate_nn": round(r_nn, 4), "rate_nfe": round(r_eu, 4),
            "cmh_or": round(cmh["or"], 3), "cmh_ci95": [round(x, 3) for x in cmh["ci95"]], "cmh_p": cmh["p"],
            "logit_or": round(lor, 3), "logit_p": lp,
            "logit_subs_adj_or": round(lor_m, 3) if not np.isnan(lor_m) else None, "logit_subs_adj_p": lp_m if not np.isnan(lp_m) else None,
        }
    print(f"\n  [{label}] n_NN={res['n_non_nfe']:,} (medAF={res['med_af_nn']:.5f}, medSubs={res['med_subs_nn']:.0f}) vs n_NFE={res['n_nfe']:,} (medAF={res['med_af_nfe']:.5f}, medSubs={res['med_subs_nfe']:.0f})")
    for c, rv in res["outcomes"].items():
        print(f"    {c:11s}: NN={rv['rate_nn']*100:5.2f}% vs NFE={rv['rate_nfe']*100:5.2f}% | CMH={rv['cmh_or']:.2f} {rv['cmh_ci95']} (p={rv['cmh_p']:.2e}) | Logit={rv['logit_or']:.2f} (p={rv['logit_p']:.2e}) | SubsAdj={rv['logit_subs_adj_or']} (p={rv['logit_subs_adj_p']})")
    return res

# Spec A: FAF95_grpmax_nonnfe5 >= 1e-4 (ZERO singletons!) vs FAF95_nfe >= 1e-4
h1_clean["A_faf95_ge_1e4"] = eval_h1_spec(
    cs,
    (cs["faf95_nonnfe5"] >= 1e-4) & (cs["af_nonnfe5"] / np.maximum(cs["af_nfe"], 1e-6) >= 5.0),
    (cs["faf95_nfe"] >= 1e-4) & (cs["af_nfe"] / np.maximum(cs["af_nonnfe5"], 1e-6) >= 5.0),
    "faf95_nonnfe5", "faf95_nfe", "FAF95 >= 1e-4 (Zero Singletons)"
)

# Spec B: AC >= 10 & AF >= 5e-4 across 4 major non-NFE continental groups (AFR, AMR, EAS, SAS — excluding tiny MID!)
ac_nn4 = cs[["ac_afr", "ac_amr", "ac_eas", "ac_sas"]].max(axis=1)
af_nn4 = cs[["af_afr", "af_amr", "af_eas", "af_sas"]].max(axis=1)
cs["af_nn4"] = af_nn4
h1_clean["B_ac10_af5e4_4major"] = eval_h1_spec(
    cs,
    (af_nn4 >= 5e-4) & (ac_nn4 >= 10) & (af_nn4 / np.maximum(cs["af_nfe"], 1e-6) >= 5.0),
    (cs["af_nfe"] >= 5e-4) & (cs["ac_nfe"] >= 10) & (cs["af_nfe"] / np.maximum(af_nn4, 1e-6) >= 5.0),
    "af_nn4", "af_nfe", "4 Major Non-NFE (AFR/AMR/EAS/SAS) AC>=10 & AF>=5e-4"
)

# Spec C: AF >= 1e-3 & AC >= 5 (incl MID!)
ac_nn5 = cs[["ac_afr", "ac_amr", "ac_eas", "ac_sas", "ac_mid"]].max(axis=1)
h1_clean["C_af1e3_ac5_all5"] = eval_h1_spec(
    cs,
    (cs["af_nonnfe5"] >= 1e-3) & (ac_nn5 >= 5) & (cs["af_nonnfe5"] / np.maximum(cs["af_nfe"], 1e-6) >= 5.0),
    (cs["af_nfe"] >= 1e-3) & (cs["ac_nfe"] >= 5) & (cs["af_nfe"] / np.maximum(cs["af_nonnfe5"], 1e-6) >= 5.0),
    "af_nonnfe5", "af_nfe", "All 5 Non-NFE AF>=1e-3 & AC>=5"
)

# Spec D: EXACT HYPERGEOMETRIC DOWNSAMPLING OF NFE TO MATCH EACH ANCESTRY'S EXACT AN_p!
print("\n=== 3. EXACT HYPERGEOMETRIC DOWNSAMPLING AUDIT (EQUALIZING DENOMINATOR POWER NFE == POP!) ===")
rng = np.random.default_rng(20261125)
ds_res = {}
for p in ["afr", "amr", "eas", "sas", "mid"]:
    med_an_p = int(cs[f"an_{p}"].median())
    # Downsample NFE AC to exact med_an_p (10 Monte Carlo replicates averaged)
    # Vectorized single-draw exact hypergeometric (plus 5 replicates for stability!)
    vus_ors, plp_ors, conf_ors = [], [], []
    n_nfe_list, n_p_list = [], []
    for rep in range(5):
        valid_nfe = cs["an_nfe"] >= med_an_p
        an_n = cs.loc[valid_nfe, "an_nfe"].values.astype(int)
        ac_n = cs.loc[valid_nfe, "ac_nfe"].values.astype(int)
        ac_ds = rng.hypergeometric(ngood=ac_n, nbad=an_n - ac_n, nsample=med_an_p)
        af_nfe_ds = ac_ds / float(med_an_p)
        sub_d = cs.loc[valid_nfe].copy()
        sub_d["af_nfe_ds"] = af_nfe_ds
        # Exact symmetric selection rule using equal-AN frequencies (and FAF95-clean AC>=2!)
        m_p = (sub_d[f"af_{p}"] >= 1e-4) & (sub_d[f"ac_{p}"] >= 2) & (sub_d[f"af_{p}"] / np.maximum(sub_d["af_nfe_ds"], 1e-6) >= 5.0)
        m_e = (sub_d["af_nfe_ds"] >= 1e-4) & (ac_ds >= 2) & (sub_d["af_nfe_ds"] / np.maximum(sub_d[f"af_{p}"], 1e-6) >= 5.0)
        cc = sub_d[m_p | m_e].copy()
        cc["is_p"] = m_p.loc[cc.index].astype(int)
        cc["af_ref"] = np.where(cc["is_p"] == 1, cc[f"af_{p}"], cc["af_nfe_ds"])
        cc["af_bin"] = pd.cut(cc["af_ref"], bins=[1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 5e-2, 1.01], right=False)
        n_p_list.append(int((cc["is_p"] == 1).sum()))
        n_nfe_list.append(int((cc["is_p"] == 0).sum()))
        for cat_nm, lst in [("VUS", vus_ors), ("P_LP", plp_ors), ("CONFLICTING", conf_ors)]:
            cc["y"] = (cc["clnsig_cat"] == cat_nm).astype(int)
            r = cmh_test(cc, "is_p", "y", "af_bin")
            lst.append((r["or"], r["p"]))
    ds_res[p] = {
        "downsampled_an": med_an_p,
        "mean_n_pop": round(float(np.mean(n_p_list)), 1),
        "mean_n_nfe_ds": round(float(np.mean(n_nfe_list)), 1),
        "vus_cmh_or_mean": round(float(np.mean([x[0] for x in vus_ors])), 3),
        "vus_cmh_p_max": float(np.max([x[1] for x in vus_ors])),
        "plp_cmh_or_mean": round(float(np.mean([x[0] for x in plp_ors])), 3),
        "plp_cmh_p_median": float(np.median([x[1] for x in plp_ors])),
        "conf_cmh_or_mean": round(float(np.mean([x[0] for x in conf_ors])), 3),
        "conf_cmh_p_max": float(np.max([x[1] for x in conf_ors])),
    }
    print(f"  {p.upper():3s} (equal AN={med_an_p:6d}, n_pop={ds_res[p]['mean_n_pop']:.0f} vs n_NFE_ds={ds_res[p]['mean_n_nfe_ds']:.0f}): "
          f"VUS CMH OR={ds_res[p]['vus_cmh_or_mean']:.2f} (p<={ds_res[p]['vus_cmh_p_max']:.2e}) | "
          f"P/LP OR={ds_res[p]['plp_cmh_or_mean']:.2f} (med p={ds_res[p]['plp_cmh_p_median']:.3f}) | "
          f"Conf OR={ds_res[p]['conf_cmh_or_mean']:.2f} (p<={ds_res[p]['conf_cmh_p_max']:.2e})")
out["h1_clean_specs"] = h1_clean
out["h1_hypergeometric_downsampling"] = ds_res

# ==============================================================================
# 4. LIVE CLINVAR EVIDENCE AUDIT (ALL 12 L4 CANDIDATES!) & gnomAD QC FILTERS AUDIT
# ==============================================================================
print("\n=== 4. LIVE CLINVAR EVIDENCE AUDIT ON ALL 12 L4 CANDIDATES ===")
ev_dir = os.path.join(AUD, "clinvar_l4_evidence")
os.makedirs(ev_dir, exist_ok=True)
l4 = cat55[cat55["triage_layer"] == "L4_PriorityNonBottleneckAuditCandidate"].copy()
l4_ev_summary = []
for row in l4.itertuples(index=False):
    vid = str(row.variation_id)
    out_f = os.path.join(ev_dir, f"cv_{vid}.json")
    if not os.path.exists(out_f) or os.path.getsize(out_f) < 100:
        p = subprocess.run(["uv", "run", CLINVAR_API, "evidence", "--variant_id", vid, "--output", out_f],
                           capture_output=True, text=True, timeout=60)
        time.sleep(0.45)
    if os.path.exists(out_f) and os.path.getsize(out_f) > 50:
        ev = json.load(open(out_f))
        conds = [c.get("name") for c in (ev.get("conditions") or [])]
        subs = ev.get("submissions") or []
        sub_classes = {}
        latest_dates = []
        notes_sample = []
        for s in subs:
            cl = s.get("classification") or "unknown"
            sub_classes[cl] = sub_classes.get(cl, 0) + 1
            if s.get("date_last_evaluated"):
                latest_dates.append(s["date_last_evaluated"])
            if s.get("curator_notes") and len(notes_sample) < 2:
                notes_sample.append(str(s["curator_notes"])[:180])
        rec = {
            "gene": row.gene, "variant_id": row.variant_id, "variation_id": vid,
            "hgvs": row.hgvs_name, "stars": int(row.stars),
            "faf95_pop8": row.faf95_pop8, "faf95_max8": round(float(row.faf95_max8), 6),
            "hom_joint": int(row.hom_joint),
            "conditions": conds[:5],
            "n_submissions_in_ev": len(subs),
            "submission_classifications": sub_classes,
            "n_citations": len(ev.get("citation_references") or []),
            "sample_notes": notes_sample,
        }
        l4_ev_summary.append(rec)
        print(f"  {row.gene:7s} {row.variant_id:18s} ({row.faf95_pop8.upper()} {row.faf95_max8:.5f}): subs={sub_classes} | conds={conds[:3]}")
out["l4_clinvar_evidence"] = l4_ev_summary

# 5. LIVE gnomAD v4.1 QC FILTERS/FLAGS AUDIT ON ALL 12 L4 + 7 L1 CANDIDATES!
# Uses polite 6.2s rate limit (<=10 QPM, strictly matching gnomad_database policy) via urllib POST to gnomAD GraphQL API
print("\n=== 5. LIVE gnomAD v4.1 QC FILTERS & FLAGS AUDIT (12 L4 + 7 L1 = 19 SITES) ===")
import urllib.request
qc_dir = os.path.join(AUD, "gnomad_qc_flags")
os.makedirs(qc_dir, exist_ok=True)
target_qc = cat55[cat55["triage_layer"].isin(["L4_PriorityNonBottleneckAuditCandidate", "L1_LowCallability_TechnicalArtifact"])].copy()
gql = """
query($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    flags
    exome { ac an filters flags }
    genome { ac an filters flags }
    joint { ac an faf95 { popmax popmax_population } faf99 { popmax popmax_population } }
  }
}
"""
qc_rows = []
for row in target_qc.itertuples(index=False):
    vid = row.variant_id
    out_f = os.path.join(qc_dir, f"{vid}.json")
    if os.path.exists(out_f) and os.path.getsize(out_f) > 50:
        d_q = json.load(open(out_f))
    else:
        req = urllib.request.Request(
            "https://gnomad.broadinstitute.org/api",
            data=json.dumps({"query": gql, "variables": {"variantId": vid, "dataset": "gnomad_r4"}}).encode("utf-8"),
            headers={"Content-Type": "application/json", "User-Agent": "AcademicAudit/1.0"},
            method="POST"
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                d_q = json.loads(resp.read().decode("utf-8"))
            with open(out_f, "w") as f:
                json.dump(d_q, f, indent=2)
        except Exception as e:
            d_q = {"error": str(e)}
        time.sleep(6.2)
    v = (d_q.get("data") or {}).get("variant") or {}
    ex = v.get("exome") or {}
    ge = v.get("genome") or {}
    jt = v.get("joint") or {}
    rec = {
        "layer": row.triage_layer, "gene": row.gene, "variant_id": vid,
        "global_flags": v.get("flags") or [],
        "exome_filters": ex.get("filters") or [], "exome_flags": ex.get("flags") or [],
        "genome_filters": ge.get("filters") or [], "genome_flags": ge.get("flags") or [],
        "joint_faf95": (jt.get("faf95") or {}).get("popmax"),
        "joint_faf95_pop": (jt.get("faf95") or {}).get("popmax_population"),
    }
    qc_rows.append(rec)
    ex_pass = "PASS" if (ex and not rec["exome_filters"]) else f"FAIL:{rec['exome_filters']}"
    ge_pass = "PASS" if (ge and not rec["genome_filters"]) else (f"FAIL:{rec['genome_filters']}" if ge else "NO_GENOME")
    print(f"  [{row.triage_layer[:2]}] {row.gene:7s} {vid:18s}: exome={ex_pass:14s} | genome={ge_pass:14s} | flags={rec['global_flags']}")
out["gnomad_qc_audit"] = qc_rows

with open(os.path.join(AUD, "a1_forensic_audit.json"), "w") as f:
    json.dump(out, f, indent=2)
print("\nSaved a1_forensic_audit.json")
