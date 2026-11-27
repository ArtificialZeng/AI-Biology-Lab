#!/opt/anaconda3/bin/python3
import gzip, json, math, os, re, time, urllib.request
import numpy as np
import pandas as pd

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
OUT_DIR = os.path.join(BASE, "falsification_ar_in_ad_pipelines")
CLONE_DIR = os.path.join(OUT_DIR, "cloned_pipelines")

df_j = pd.read_parquet(os.path.join(BASE, "data/clean/a1_clinvar_gnomad_exact_joined.parquet"))
print("Columns in joined parquet:", list(df_j.columns))

og_tsv = pd.read_csv(os.path.join(CLONE_DIR, "OwnGenome/assets/acmg_sf.tsv"), sep="\t", comment="#")
ttn_limit = og_tsv.loc[og_tsv["gene"] == "TTN", "limit"].iloc[0]
ttn_hi_psi = [(int(m.group(1)), int(m.group(2))) for m in re.finditer(r"chr2:(\d+)-(\d+)", ttn_limit)]
def in_ttn_hi_psi(pos):
    p = int(pos)
    return any(s <= p <= e for s, e in ttn_hi_psi)
lof_no_ad_genes = set(og_tsv.loc[(og_tsv["inheritance"] == "AD") & (og_tsv["lof"] == "no"), "gene"])

def is_plp(s):
    if pd.isna(s): return False
    s = str(s)
    return ("Conflicting" not in s) and (("Pathogenic" in s) or ("Likely_pathogenic" in s))

ad_j = df_j[(df_j["inheritance"] == "AD") & (df_j["clnsig"].apply(is_plp))].copy()

with open(os.path.join(OUT_DIR, "genome_wide_burden_and_patch_audit.json")) as fh:
    gw = json.load(fh)
ryr1_cm_ar_vids = set(gw["ryr1_3star_ep_audit"]["cm_vcep_ar_only_vids"])

with open(os.path.join(OUT_DIR, "ryr1_all_3star_live_vcep_audit.json")) as fh:
    ryr1_ep_audit = json.load(fh)
ryr1_mh_ad_vids = set()
for vid, info in ryr1_ep_audit.items():
    if any("malignant hyperthermia" in p["panel"].lower() for p in info.get("ep_panels", [])):
        ryr1_mh_ad_vids.add(str(vid))

vids_all = set(ad_j["variation_id"].astype(str))
vcf_map = {}
with gzip.open(os.path.join(BASE, "data/raw/clinvar.vcf.gz"), "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if p[2] in vids_all:
            info = dict(kv.split("=", 1) for kv in p[7].split(";") if "=" in kv)
            vcf_map[p[2]] = (info.get("CLNDN", ""), info.get("CLNDISDB", ""), info.get("MC", ""))

vs = ad_j["variation_id"].astype(str)
ad_j["vcf_clndn"] = vs.map(lambda x: vcf_map.get(x, ("", "", ""))[0])
ad_j["vcf_clndisdb"] = vs.map(lambda x: vcf_map.get(x, ("", "", ""))[1])
ad_j["vcf_mc"] = vs.map(lambda x: vcf_map.get(x, ("", "", ""))[2])
ad_j["is_plof_mc"] = ad_j["vcf_mc"].str.lower().apply(
    lambda s: any(k in s for k in ["nonsense", "frameshift", "splice_donor", "splice_acceptor"])
)
ad_j["in_ttn_hi_psi"] = np.where(ad_j["gene"] == "TTN", ad_j["pos"].apply(in_ttn_hi_psi), False)

def strict_ground_truth(r):
    g = r["gene"]
    dn = str(r["vcf_clndn"]).lower()
    dis = str(r["vcf_clndisdb"]).lower()
    plof = bool(r["is_plof_mc"])
    vid = str(int(r["variation_id"]))
    if g == "RYR1":
        if vid in ryr1_mh_ad_vids: return "AD_ELIGIBLE"
        if vid in ryr1_cm_ar_vids: return "RYR1_ClinGen_CM_VCEP_AR"
        if plof: return "RYR1_pLoF_AR_Myopathy"
        if ("malignant_hyperthermia" not in dn) and ("omim:145600" not in dis): return "RYR1_NonMH_Myopathy_NonLoF"
    elif g == "TTN":
        if not plof: return "TTN_NonTruncating"
        if not bool(r["in_ttn_hi_psi"]): return "TTN_LowPSI_Truncating"
    elif g == "APOB":
        if plof: return "APOB_pLoF_FHBL_LowLDL"
    elif g == "PCSK9":
        if plof: return "PCSK9_pLoF_LowLDL"
    elif g == "RET":
        has_men2 = any(k in dn for k in ["multiple_endocrine_neoplasia", "medullary_thyroid"]) or any(k in dis for k in ["omim:171400", "omim:155240", "omim:162300"])
        if plof or (("hirschsprung" in dn or "omim:142623" in dis) and not has_men2): return "RET_pLoF_or_Hirschsprung"
    elif g == "VHL":
        has_vhl = ("von_hippel" in dn) or ("omim:193300" in dis)
        if (("erythrocytosis" in dn or "polycythemia" in dn or "chuvash" in dn or "omim:263400" in dis) and not has_vhl): return "VHL_AR_Erythrocytosis"
    elif g == "LMNA":
        has_card = ("cardiomyopathy" in dn) or ("cardiac" in dn) or ("omim:115200" in dis)
        if any(k in dn for k in ["charcot-marie-tooth", "mandibuloacral", "lipodystrophy", "progeria"]) and not has_card: return "LMNA_NonCardiac_or_AR"
    elif g in ("DES", "CACNA1S", "MYH7", "MYL2", "MYL3", "TNNI3", "TNNT2", "ACTA2", "MYH11", "TTR"):
        if plof:
            # Exclude verified 3* AD TNNT2 p.Trp297Ter (VID=177636) from discordant!
            if vid == "177636":
                return "AD_ELIGIBLE"
            return "Other_NonLoF_AD_Gene_pLoF"
    return "AD_ELIGIBLE"

ad_j["gt_cat"] = ad_j.apply(strict_ground_truth, axis=1)
ad_j["is_discordant"] = ad_j["gt_cat"] != "AD_ELIGIBLE"

MAX_AN_AUTO = 1614324.0
ad_j["an_frac"] = ad_j["an_joint"] / MAX_AN_AUTO
ad_j["het_joint"] = (ad_j["ac_joint"] - 2 * ad_j["hom_joint"]).clip(lower=0)
ad_j["het_rate_joint"] = (ad_j["het_joint"] / (ad_j["an_joint"] / 2.0).replace(0, np.nan)).clip(upper=1.0)
for pop in ["nfe", "afr", "amr", "eas", "sas"]:
    ad_j[f"het_{pop}"] = (ad_j[f"ac_{pop}"] - 2 * ad_j[f"hom_{pop}"]).clip(lower=0)
    ad_j[f"het_rate_{pop}"] = (ad_j[f"het_{pop}"] / (ad_j[f"an_{pop}"] / 2.0).replace(0, np.nan)).clip(upper=1.0)

def emits_tier1(r): return True
def emits_tier2(r): return not (r["gene"] == "TTN" and not r["is_plof_mc"])
def emits_tier3(r): return not (r["gene"] == "TTN" and ((not r["is_plof_mc"]) or (not r["in_ttn_hi_psi"])))

ad_j["tier1"] = ad_j.apply(emits_tier1, axis=1)
ad_j["tier2"] = ad_j.apply(emits_tier2, axis=1)
ad_j["tier3"] = ad_j.apply(emits_tier3, axis=1)

# 1. Coverage/callset regime breakdown of ALL discordant alleles
def an_regime(an):
    if an < 160000:
        return "genome_only_or_low_exome (<160k AN; <=10% joint)"
    elif an < 0.8 * MAX_AN_AUTO:
        return "partial_exome_capture (160k-1.29M AN; 10-80%)"
    elif an < 0.95 * MAX_AN_AUTO:
        return "high_callability_80_95pct (1.29M-1.53M AN)"
    else:
        return "near_complete_joint_>=95pct (>=1.53M AN)"

ad_j["an_regime"] = ad_j["an_joint"].apply(an_regime)
regime_summary = []
for reg, grp in ad_j[ad_j["is_discordant"]].groupby("an_regime"):
    regime_summary.append({
        "an_regime": reg,
        "n_alleles_all_stars": int(len(grp)),
        "n_alleles_2to4star": int((grp["stars"] >= 2).sum()),
        "sum_het_instances_all_stars": int(grp["het_joint"].sum()),
        "sum_het_instances_2to4star": int(grp.loc[grp["stars"] >= 2, "het_joint"].sum()),
        "boole_upper_per1k_2to4star": round(float(grp.loc[grp["stars"] >= 2, "het_rate_joint"].sum() * 1000.0), 4),
        "genes_top": grp["gene"].value_counts().head(5).to_dict(),
    })

# 2. Compute Boole upper bound (sum p_i), Independence-model estimate (1 - prod(1-p_i)), and Frechet lower bound (max p_i)
# across FOUR callability/detection-scope regimes:
#   (a) All stars, all AN
#   (b) 2-4*, all AN
#   (c) 2-4*, >=80% AN (exome+genome high callability)
#   (d) 2-4*, >=95% AN (near-complete joint exome+genome capture)
def rate_bounds(series):
    p = series.dropna().values
    p = p[(p > 0) & (p <= 1.0)]
    if len(p) == 0:
        return {"boole_upper_per1k": 0.0, "indep_model_per1k": 0.0, "max_single_lower_per1k": 0.0, "max_overlap_gap_pct": 0.0}
    boole = float(np.sum(p))
    # log1p(-p) accurate for small p
    indep = float(-np.expm1(np.sum(np.log1p(-p))))
    lb = float(np.max(p))
    gap_pct = 100.0 * (boole - indep) / boole if boole > 0 else 0.0
    return {
        "boole_upper_per1k": round(boole * 1000.0, 4),
        "indep_model_per1k": round(indep * 1000.0, 4),
        "max_single_lower_per1k": round(lb * 1000.0, 4),
        "boole_vs_indep_diff_pct": round(gap_pct, 4),
    }

bounds_table = []
for subset_lbl, mask in [
    ("ALL_STARS_ALL_AN", pd.Series(True, index=ad_j.index)),
    ("2to4star_ALL_AN", ad_j["stars"] >= 2),
    ("2to4star_GE80pctAN", (ad_j["stars"] >= 2) & (ad_j["an_frac"] >= 0.80)),
    ("2to4star_GE95pctAN", (ad_j["stars"] >= 2) & (ad_j["an_frac"] >= 0.95)),
]:
    sub = ad_j[mask]
    for tier in ["tier1", "tier2", "tier3"]:
        fp = sub[sub["is_discordant"] & sub[tier]]
        em = sub[sub[tier]]
        row = {"subset": subset_lbl, "tier": tier, "n_fp_alleles": int(len(fp)), "sum_fp_het_instances": int(fp["het_joint"].sum())}
        for pop in ["joint", "nfe", "afr", "amr", "eas", "sas"]:
            b_fp = rate_bounds(fp[f"het_rate_{pop}"])
            b_em = rate_bounds(em[f"het_rate_{pop}"])
            row[f"{pop}_fp_boole_per1k"] = b_fp["boole_upper_per1k"]
            row[f"{pop}_fp_indep_per1k"] = b_fp["indep_model_per1k"]
            row[f"{pop}_fp_max_lb_per1k"] = b_fp["max_single_lower_per1k"]
            row[f"{pop}_fp_boole_vs_indep_gap_pct"] = b_fp["boole_vs_indep_diff_pct"]
            row[f"{pop}_emitted_boole_per1k"] = b_em["boole_upper_per1k"]
            row[f"{pop}_emitted_indep_per1k"] = b_em["indep_model_per1k"]
            row[f"{pop}_fp_share_boole_pct"] = round(100.0 * b_fp["boole_upper_per1k"] / b_em["boole_upper_per1k"], 2) if b_em["boole_upper_per1k"] > 0 else 0.0
            row[f"{pop}_fp_share_indep_pct"] = round(100.0 * b_fp["indep_model_per1k"] / b_em["indep_model_per1k"], 2) if b_em["indep_model_per1k"] > 0 else 0.0
        row["sas_nfe_ratio_boole"] = round(row["sas_fp_boole_per1k"] / row["nfe_fp_boole_per1k"], 2) if row["nfe_fp_boole_per1k"] > 0 else None
        row["sas_nfe_ratio_indep"] = round(row["sas_fp_indep_per1k"] / row["nfe_fp_indep_per1k"], 2) if row["nfe_fp_indep_per1k"] > 0 else None
        bounds_table.append(row)

# 3. Within-gene bounds (2-4* all AN and >=80% AN, plus 3* RYR1!)
gene_bounds = []
for g in ["RYR1", "APOB", "TTN", "DES", "CACNA1S"]:
    for s_lbl, sub in [("2to4star_allAN", ad_j[(ad_j["gene"] == g) & (ad_j["stars"] >= 2)]),
                       ("2to4star_ge80AN", ad_j[(ad_j["gene"] == g) & (ad_j["stars"] >= 2) & (ad_j["an_frac"] >= 0.80)])]:
        fp = sub[sub["is_discordant"]]
        b_fp = rate_bounds(fp["het_rate_joint"])
        b_all = rate_bounds(sub["het_rate_joint"])
        # Also note raw allele-instance share vs rate-normalized share!
        raw_fp_inst = int(fp["het_joint"].sum())
        raw_all_inst = int(sub["het_joint"].sum())
        gene_bounds.append({
            "gene": g,
            "subset": s_lbl,
            "n_fp_alleles": int(len(fp)),
            "n_all_alleles": int(len(sub)),
            "sum_fp_het_instances": raw_fp_inst,
            "sum_all_het_instances": raw_all_inst,
            "raw_instance_share_pct": round(100.0 * raw_fp_inst / raw_all_inst, 2) if raw_all_inst > 0 else 0.0,
            "fp_boole_per1k": b_fp["boole_upper_per1k"],
            "fp_indep_per1k": b_fp["indep_model_per1k"],
            "fp_max_lb_per1k": b_fp["max_single_lower_per1k"],
            "all_boole_per1k": b_all["boole_upper_per1k"],
            "all_indep_per1k": b_all["indep_model_per1k"],
            "rate_share_boole_pct": round(100.0 * b_fp["boole_upper_per1k"] / b_all["boole_upper_per1k"], 2) if b_all["boole_upper_per1k"] > 0 else 0.0,
            "rate_share_indep_pct": round(100.0 * b_fp["indep_model_per1k"] / b_all["indep_model_per1k"], 2) if b_all["indep_model_per1k"] > 0 else 0.0,
        })

# 4. Query live gnomAD v4.1 GraphQL API for EXACT exome vs genome vs joint breakdown on:
#    - ALL 10 3* ClinGen Congenital Myopathies VCEP AR-only RYR1 alleles
#    - Top 5 APOB pLoF FHBL alleles
#    - Top 3 DES/CACNA1S pLoF alleles
# This directly inspects exome vs genome detection scope AND confirms whether joint AC/AN == exome + genome (deduplicated)!
opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
ryr1_10 = pd.read_csv(os.path.join(OUT_DIR, "ryr1_10_3star_cm_vcep_ar_alleles_gnomad_v41.csv"))
apob_top5 = ad_j[(ad_j["gene"] == "APOB") & (ad_j["is_discordant"]) & (ad_j["stars"] >= 2)].sort_values("het_joint", ascending=False).head(5)
targets = [(r["variant_id"], "RYR1_3star_CM_VCEP_AR", int(r["het_joint"]), int(r["an_joint"])) for _, r in ryr1_10.iterrows()]
targets += [(r["variant_id"], "APOB_2to4star_FHBL_pLoF", int(r["het_joint"]), int(r["an_joint"])) for _, r in apob_top5.iterrows()]

gql_query = """
query VariantCallsets($vid: String!) {
  variant(variantId: $vid, dataset: gnomad_r4) {
    variant_id
    exome { ac an ac_hom filters }
    genome { ac an ac_hom filters }
    joint { ac an homozygote_count filters }
  }
}
"""
gql_results = []
for vid, grp_lbl, local_het, local_an in targets:
    try:
        payload = json.dumps({"query": gql_query, "variables": {"vid": vid}}).encode("utf-8")
        req = urllib.request.Request(
            "https://gnomad.broadinstitute.org/api",
            data=payload,
            headers={"Content-Type": "application/json", "User-Agent": "A1X-CallsetAudit/1.0"},
        )
        with opener.open(req, timeout=20) as r:
            d = json.loads(r.read().decode("utf-8"))
        v = (d.get("data") or {}).get("variant")
        if v:
            ex = v.get("exome") or {}
            gn = v.get("genome") or {}
            jt = v.get("joint") or {}
            ex_ac, ex_an = ex.get("ac", 0) or 0, ex.get("an", 0) or 0
            gn_ac, gn_an = gn.get("ac", 0) or 0, gn.get("an", 0) or 0
            jt_ac, jt_an = jt.get("ac", 0) or 0, jt.get("an", 0) or 0
            gql_results.append({
                "variant_id": vid,
                "group": grp_lbl,
                "exome_ac": ex_ac, "exome_an": ex_an, "exome_het": max(0, ex_ac - 2*(ex.get("ac_hom",0) or 0)), "exome_filters": ex.get("filters"),
                "genome_ac": gn_ac, "genome_an": gn_an, "genome_het": max(0, gn_ac - 2*(gn.get("ac_hom",0) or 0)), "genome_filters": gn.get("filters"),
                "joint_ac": jt_ac, "joint_an": jt_an, "joint_het": max(0, jt_ac - 2*(jt.get("homozygote_count",0) or 0)),
                "sum_ex_gn_ac": ex_ac + gn_ac,
                "sum_ex_gn_an": ex_an + gn_an,
                "joint_an_minus_sum": jt_an - (ex_an + gn_an),
                "joint_ac_minus_sum": jt_ac - (ex_ac + gn_ac),
                "exome_only": (ex_an > 0) and (gn_an == 0),
                "genome_only": (gn_an > 0) and (ex_an == 0),
                "both_exome_and_genome": (ex_an > 0) and (gn_an > 0),
                "local_parquet_an_matches_joint": (local_an == jt_an),
            })
        else:
            gql_results.append({"variant_id": vid, "group": grp_lbl, "error": "variant_null", "raw": d.get("errors")})
    except Exception as e:
        gql_results.append({"variant_id": vid, "group": grp_lbl, "error": str(e)})
    time.sleep(0.25)

out = {
    "methodological_boundary": (
        "gnomAD v4.1 provides variant-level aggregate AC/AN/HOM/HEMI counts without individual-level joint genotypes. "
        "Therefore cross-variant totals represent cumulative variant-level heterozygous allele-carrier instances (Boole/Frechet upper bound on unique carriers) "
        "rather than exact deduplicated unique individuals. We additionally compute the independence-assumption Poisson-binomial prevalence estimate "
        "(1 - prod(1 - p_i)) and the single-variant Frechet lower bound (max p_i), alongside AN detection-scope strata (>=80%, >=95%, genome-only <160k) "
        "and live gnomAD v4.1 GraphQL exome-vs-genome-vs-joint decomposition."
    ),
    "an_coverage_regimes_discordant": regime_summary,
    "tier_burden_bounds_table": bounds_table,
    "within_gene_bounds_table": gene_bounds,
    "live_gnomad_exome_genome_joint_audit": gql_results,
}

out_json = os.path.join(OUT_DIR, "gnomad_burden_bounds_and_callsets_audit.json")
with open(out_json, "w") as fh:
    json.dump(out, fh, indent=2)
pd.DataFrame(bounds_table).to_csv(os.path.join(OUT_DIR, "gnomad_burden_bounds_by_callability.csv"), index=False)
print("WROTE:", out_json)
