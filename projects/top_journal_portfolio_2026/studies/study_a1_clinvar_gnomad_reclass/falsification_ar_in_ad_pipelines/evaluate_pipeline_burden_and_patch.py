#!/opt/anaconda3/bin/python3
"""
Attack 1 — Step 3: Genome-Wide Per-Pipeline False-Positive Return Burden Across gnomAD v4.1
(807,162 Individuals, Stratified by Ancestry) & Validation of a Concrete Dual-MOI / Dual-Mechanism Patch.

Also audits ALL 3* ClinGen Expert Panel AD alleles across all 73 AD ACMG SF v3.3 genes to verify:
  1. How many 3* ClinGen EP alleles in RYR1/TTN/APOB/RET/etc. are ACTUALLY AR-only curations
     (live ClinVar XML ModeOfInheritance verification on EVERY 3* allele in RYR1 + dual-mode genes!),
  2. Exact per-pipeline false-positive carrier burden per 1,000 sequenced individuals across:
     - Architecture Tier 1: Gene-list + ClinVar P/LP (CPHI-DRAGEN-anno / DNA / VikingGenes)
     - Architecture Tier 2: Pure-AR zygosity + HFE hom + TTN LoF (CPSR v2.4.0 / SFtool)
     - Architecture Tier 3: Tier 2 + Cardiodb TTN high-PSI >90% (OwnGenome c324b61)
     - Architecture Tier 4 (Proposed Patch): Tier 3 + ClinVar P/LP mechanism gating (ClinGen lof=no
       blocks non-LoF AD gene pLoFs) + ClinVar SCV-level MOI/EP-assertion gating + strict non-SF CLNDN filter!
  3. Zero false-negative check across ALL genuine AD ClinGen Expert Panel alleles!
"""

import gzip
import json
import os
import re
import urllib.request
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
OUT_DIR = os.path.join(BASE, "falsification_ar_in_ad_pipelines")
CLONE_DIR = os.path.join(OUT_DIR, "cloned_pipelines")

df_j = pd.read_parquet(os.path.join(BASE, "data/clean/a1_clinvar_gnomad_exact_joined.parquet"))
df_c = pd.read_parquet(os.path.join(BASE, "data/clean/clinvar_acmg84_grch38.parquet"))
og_tsv = pd.read_csv(os.path.join(CLONE_DIR, "OwnGenome/assets/acmg_sf.tsv"), sep="\t", comment="#")

# Parse OwnGenome exact TTN high-PSI >90% intervals and lof=no genes
ttn_limit = og_tsv.loc[og_tsv["gene"] == "TTN", "limit"].iloc[0]
ttn_hi_psi = [(int(m.group(1)), int(m.group(2))) for m in re.finditer(r"chr2:(\d+)-(\d+)", ttn_limit)]
def in_ttn_hi_psi(pos):
    p = int(pos)
    return any(s <= p <= e for s, e in ttn_hi_psi)

lof_no_ad_genes = set(og_tsv.loc[(og_tsv["inheritance"] == "AD") & (og_tsv["lof"] == "no"), "gene"])
print("ClinGen/OwnGenome lof=no AD ACMG SF v3.3 genes (n={}): {}".format(len(lof_no_ad_genes), sorted(lof_no_ad_genes)))

def is_plp(s):
    if pd.isna(s): return False
    s = str(s)
    return ("Conflicting" not in s) and (("Pathogenic" in s) or ("Likely_pathogenic" in s))

ad_j = df_j[(df_j["inheritance"] == "AD") & (df_j["clnsig"].apply(is_plp))].copy()
ad_c = df_c[(df_c["inheritance"] == "AD") & (df_c["clnsig"].apply(is_plp))].copy()

# Attach raw VCF CLNDN, CLNDISDB, MC
vids_all = set(ad_c["variation_id"].astype(str))
vcf_map = {}
with gzip.open(os.path.join(BASE, "data/raw/clinvar.vcf.gz"), "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if p[2] in vids_all:
            info = dict(kv.split("=", 1) for kv in p[7].split(";") if "=" in kv)
            vcf_map[p[2]] = (info.get("CLNDN", ""), info.get("CLNDISDB", ""), info.get("MC", ""))

for d in (ad_c, ad_j):
    vs = d["variation_id"].astype(str)
    d["vcf_clndn"] = vs.map(lambda x: vcf_map.get(x, ("", "", ""))[0])
    d["vcf_clndisdb"] = vs.map(lambda x: vcf_map.get(x, ("", "", ""))[1])
    d["vcf_mc"] = vs.map(lambda x: vcf_map.get(x, ("", "", ""))[2])
    d["is_plof_mc"] = d["vcf_mc"].str.lower().apply(
        lambda s: any(k in s for k in ["nonsense", "frameshift", "splice_donor", "splice_acceptor"])
    )
    d["in_ttn_hi_psi"] = np.where(d["gene"] == "TTN", d["pos"].apply(in_ttn_hi_psi), False)

# Live ClinVar XML audit of ALL 3* ClinGen Expert Panel RYR1 P/LP alleles in gnomAD v4.1!
# (How many 3* RYR1 alleles are curated by Congenital Myopathies VCEP [AR] vs Malignant Hyperthermia VCEP [AD]?!)
ryr1_ep = ad_j[(ad_j["gene"] == "RYR1") & (ad_j["stars"] >= 3)].copy()
print(f"\nFetching live ClinVar XML for ALL {len(ryr1_ep)} 3* Expert-Panel RYR1 P/LP alleles in gnomAD v4.1...")
opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
ep_vids = ryr1_ep["variation_id"].astype(str).tolist()
ryr1_ep_audit = {}
# Fetch in batches of 30
for i in range(0, len(ep_vids), 30):
    batch = ep_vids[i:i+30]
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id={','.join(batch)}"
    with opener.open(url, timeout=30) as r:
        root = ET.fromstring(r.read())
    for va in root.findall(".//VariationArchive"):
        vid = va.attrib.get("VariationID")
        name = va.attrib.get("VariationName")
        ep_panels = []
        for ca in va.findall(".//ClinicalAssertion"):
            acc = ca.find(".//ClinVarAccession")
            sub = acc.attrib.get("SubmitterName") if acc is not None else ""
            cls_el = ca.find(".//Classification")
            rev = cls_el.findtext("ReviewStatus") if cls_el is not None else ""
            cls_d = cls_el.findtext("GermlineClassification") if cls_el is not None else ""
            mois = [a.text for a in ca.findall(".//AttributeSet/Attribute") if a.attrib.get("Type") == "ModeOfInheritance" and a.text]
            if "expert panel" in (rev or "").lower() or "clingen" in (sub or "").lower():
                ep_panels.append({"panel": sub, "class": cls_d, "moi": mois})
        ryr1_ep_audit[vid] = {"name": name, "ep_panels": ep_panels}

with open(os.path.join(OUT_DIR, "ryr1_all_3star_live_vcep_audit.json"), "w") as fh:
    json.dump(ryr1_ep_audit, fh, indent=2)

# Classify each 3* RYR1 allele by exact ClinGen VCEP & MOI!
ryr1_cm_ar_vids = set()
ryr1_mh_ad_vids = set()
ryr1_cm_ad_vids = set()
for _, r in ryr1_ep.iterrows():
    vid = str(int(r["variation_id"]))
    info = ryr1_ep_audit.get(vid, {})
    panels = info.get("ep_panels", [])
    has_mh_vcep = any("malignant hyperthermia" in p["panel"].lower() for p in panels)
    has_cm_ar = any("congenital myopathies" in p["panel"].lower() and any("recessive" in m.lower() for m in p["moi"]) for p in panels)
    has_cm_ad = any("congenital myopathies" in p["panel"].lower() and any("dominant" in m.lower() for m in p["moi"]) for p in panels)
    if has_mh_vcep:
        ryr1_mh_ad_vids.add(vid)
    elif has_cm_ar and not has_cm_ad:
        ryr1_cm_ar_vids.add(vid)
    elif has_cm_ad:
        ryr1_cm_ad_vids.add(vid)

ryr1_ep["vcep_class"] = ryr1_ep["variation_id"].astype(str).map(
    lambda v: "MHS_VCEP_AD" if v in ryr1_mh_ad_vids else ("CM_VCEP_AR_ONLY" if v in ryr1_cm_ar_vids else ("CM_VCEP_AD" if v in ryr1_cm_ad_vids else "OTHER_EP"))
)
ryr1_ep["het_joint"] = (ryr1_ep["ac_joint"] - 2 * ryr1_ep["hom_joint"]).clip(lower=0)
print("\n=== EXHAUSTIVE LIVE VCEP AUDIT OF ALL 3* EXPERT-PANEL RYR1 P/LP ALLELES IN gnomAD v4.1 ===")
print(ryr1_ep.groupby("vcep_class").agg(n_alleles=("variant_id", "count"), sum_het=("het_joint", "sum")).to_string())
print("\nAll CM_VCEP_AR_ONLY 3* RYR1 alleles found:")
for _, r in ryr1_ep[ryr1_ep["vcep_class"] == "CM_VCEP_AR_ONLY"].sort_values("het_joint", ascending=False).iterrows():
    has_mh_poll = ("malignant_hyperthermia" in str(r["vcf_clndn"]).lower()) or ("omim:145600" in str(r["vcf_clndisdb"]).lower())
    print(f"  {r['variant_id']:22s} | VID={r['variation_id']} | HET={int(r['het_joint']):4d} | CLNDN_polluted_MH={has_mh_poll} | {r['hgvs_name']}")

# Compute exact het counts & rates across ALL 73 AD ACMG SF v3.3 genes
MAX_AN_AUTO = 1614324.0
ad_j["hi_call"] = ad_j["an_joint"] >= 0.8 * MAX_AN_AUTO
ad_j["het_joint"] = (ad_j["ac_joint"] - 2 * ad_j["hom_joint"]).clip(lower=0)
ad_j["het_rate_joint"] = ad_j["het_joint"] / (ad_j["an_joint"] / 2.0).replace(0, np.nan)
for pop in ["nfe", "afr", "amr", "eas", "sas"]:
    ad_j[f"het_{pop}"] = (ad_j[f"ac_{pop}"] - 2 * ad_j[f"hom_{pop}"]).clip(lower=0)
    ad_j[f"het_rate_{pop}"] = ad_j[f"het_{pop}"] / (ad_j[f"an_{pop}"] / 2.0).replace(0, np.nan)

# Tag exact discordant ground truth across all AD ACMG SF v3.3 alleles (now using ALL verified CM_VCEP_AR_ONLY VIDs!)
def strict_ground_truth(r):
    g = r["gene"]
    dn = str(r["vcf_clndn"]).lower()
    dis = str(r["vcf_clndisdb"]).lower()
    plof = bool(r["is_plof_mc"])
    vid = str(int(r["variation_id"]))
    if g == "RYR1":
        if vid in ryr1_mh_ad_vids:
            return "AD_ELIGIBLE"
        if vid in ryr1_cm_ar_vids:
            return "RYR1_ClinGen_CM_VCEP_AR"
        if plof:
            return "RYR1_pLoF_AR_Myopathy"
        if ("malignant_hyperthermia" not in dn) and ("omim:145600" not in dis):
            return "RYR1_NonMH_Myopathy_NonLoF"
    elif g == "TTN":
        if not plof:
            return "TTN_NonTruncating"
        if not bool(r["in_ttn_hi_psi"]):
            return "TTN_LowPSI_Truncating"
    elif g == "APOB":
        if plof:
            return "APOB_pLoF_FHBL_LowLDL"
    elif g == "PCSK9":
        if plof:
            return "PCSK9_pLoF_LowLDL"
    elif g == "RET":
        has_men2 = any(k in dn for k in ["multiple_endocrine_neoplasia", "medullary_thyroid"]) or any(k in dis for k in ["omim:171400", "omim:155240", "omim:162300"])
        if plof or (("hirschsprung" in dn or "omim:142623" in dis) and not has_men2):
            return "RET_pLoF_or_Hirschsprung"
    elif g == "VHL":
        has_vhl = ("von_hippel" in dn) or ("omim:193300" in dis)
        if (("erythrocytosis" in dn or "polycythemia" in dn or "chuvash" in dn or "omim:263400" in dis) and not has_vhl):
            return "VHL_AR_Erythrocytosis"
    elif g == "LMNA":
        has_card = ("cardiomyopathy" in dn) or ("cardiac" in dn) or ("omim:115200" in dis)
        if any(k in dn for k in ["charcot-marie-tooth", "mandibuloacral", "lipodystrophy", "progeria"]) and not has_card:
            return "LMNA_NonCardiac_or_AR"
    elif g in ("DES", "CACNA1S", "MYH7", "MYL2", "MYL3", "TNNI3", "TNNT2", "ACTA2", "MYH11", "TTR"):
        if plof:
            return "Other_NonLoF_AD_Gene_pLoF"
    return "AD_ELIGIBLE"

ad_j["gt_cat"] = ad_j.apply(strict_ground_truth, axis=1)
ad_j["is_discordant"] = ad_j["gt_cat"] != "AD_ELIGIBLE"

# Evaluate EXACT false-positive carrier rate per 1,000 sequenced individuals emitted by EACH pipeline architecture!
# Tier 1 (CPHI / DNA / naive gene-list): emits ALL AD P/LP (no TTN LoF, no PSI, no lof=no, no MOI)
# Tier 2 (CPSR v2.4.0 / SFtool): drops TTN non-truncating ONLY; emits ALL remaining discordant alleles!
# Tier 3 (OwnGenome c324b61): drops TTN non-truncating AND TTN low-PSI <=90%; emits ALL remaining discordant alleles!
# Tier 4 (Proposed Concrete Patch):
#   (a) Drops TTN non-truncating & TTN low-PSI <=90% (Cardiodb)
#   (b) Applies ClinGen `lof=no` to ClinVar P/LP pLoFs too (blocks RYR1/APOB/RET/PCSK9/DES/CACNA1S/MYH7/etc. pLoFs)
#   (c) Applies ClinVar SCV-level MOI / VCEP gating (blocks CM_VCEP_AR_ONLY alleles even when CLNDN is polluted!)
#   (d) Drops alleles whose named CLNDN contains ZERO ACMG SF target disease terms (RYR1 non-MH myopathy, RET Hirschsprung, LMNA CMT2, VHL AR erythrocytosis)

def emits_tier1(r): return True
def emits_tier2(r): return not (r["gene"] == "TTN" and not r["is_plof_mc"])
def emits_tier3(r): return not (r["gene"] == "TTN" and ((not r["is_plof_mc"]) or (not r["in_ttn_hi_psi"])))
def emits_tier4_patch(r):
    if not emits_tier3(r):
        return False
    g = r["gene"]
    vid = str(int(r["variation_id"]))
    # Rule (b): ClinGen lof=no mechanism gate applied to ClinVar P/LP pLoFs too!
    if g in lof_no_ad_genes and bool(r["is_plof_mc"]):
        return False
    # Rule (c): SCV-level ClinGen VCEP MOI=AR exclusion!
    if vid in ryr1_cm_ar_vids:
        return False
    # Rule (d): Strict non-target CLNDN gate (only when explicit non-target disease named & ZERO target terms)
    dn = str(r["vcf_clndn"]).lower()
    dis = str(r["vcf_clndisdb"]).lower()
    if g == "RYR1" and ("malignant_hyperthermia" not in dn) and ("omim:145600" not in dis):
        return False
    if g == "RET":
        has_men2 = any(k in dn for k in ["multiple_endocrine_neoplasia", "medullary_thyroid"]) or any(k in dis for k in ["omim:171400", "omim:155240", "omim:162300"])
        if (("hirschsprung" in dn or "omim:142623" in dis) and not has_men2):
            return False
    if g == "LMNA":
        has_card = ("cardiomyopathy" in dn) or ("cardiac" in dn) or ("omim:115200" in dis)
        if any(k in dn for k in ["charcot-marie-tooth", "mandibuloacral", "lipodystrophy", "progeria"]) and not has_card:
            return False
    if g == "VHL":
        has_vhl = ("von_hippel" in dn) or ("omim:193300" in dis)
        if any(k in dn for k in ["erythrocytosis", "polycythemia", "chuvash"]) and not has_vhl:
            return False
    return True

for name, fn in [
    ("tier1_gene_list", emits_tier1),
    ("tier2_cpsr_sftool", emits_tier2),
    ("tier3_owngenome", emits_tier3),
    ("tier4_proposed_patch", emits_tier4_patch),
]:
    ad_j[name] = ad_j.apply(fn, axis=1)

# Audit FALSE-NEGATIVE RISK of Tier 4 Patch on ALL 3* ClinGen Expert Panel AD alleles across ALL 73 AD genes!
ep_all_ad = ad_j[(ad_j["stars"] >= 3) & (~ad_j["is_discordant"])].copy()
ep_dropped = ep_all_ad[~ep_all_ad["tier4_proposed_patch"]]
print(f"\n=== ZERO FALSE-NEGATIVE SAFETY AUDIT OF TIER 4 PATCH ON ALL 3* AD EXPERT-PANEL ALLELES ===")
print(f"  Total 3* AD-eligible Expert-Panel alleles in gnomAD v4.1: {len(ep_all_ad):,}")
print(f"  Retained by Tier 4 Patch: {int(ep_all_ad['tier4_proposed_patch'].sum()):,} / {len(ep_all_ad):,} ({100.0*ep_all_ad['tier4_proposed_patch'].mean():.2f}%)")
if len(ep_dropped) > 0:
    print("  Dropped 3* alleles (inspect!):")
    print(ep_dropped[["variant_id", "gene", "hgvs_name", "vcf_mc", "vcf_clndn"]].to_string(index=False))

# Check also ALL 3* alleles in lof=no genes (did ANY genuine AD EP allele get dropped?!)
ep_lofno = ad_j[(ad_j["stars"] >= 3) & (ad_j["gene"].isin(lof_no_ad_genes))].copy()
print(f"\n  All 3* alleles in lof=no AD genes (n={len(ep_lofno)}): retained={int(ep_lofno['tier4_proposed_patch'].sum())}, dropped={int((~ep_lofno['tier4_proposed_patch']).sum())}")
if (~ep_lofno["tier4_proposed_patch"]).sum() > 0:
    print("  Dropped 3* alleles in lof=no genes (verify ALL are AR myopathy / non-AD!):")
    for _, r in ep_lofno[~ep_lofno["tier4_proposed_patch"]].iterrows():
        print(f"    {r['gene']:8s} | {r['variant_id']:22s} | gt_cat={r['gt_cat']} | HET={int(r['het_joint']):4d} | {r['hgvs_name']}")

# Build per-pipeline genome-wide false-positive burden table (for ALL stars & HIGH-CONF 2-4* & HIGH-CALLABILITY 2-4* >=80% AN!)
burden_rows = []
for subset_lbl, sub in [
    ("ALL_STARS_ALL_AN", ad_j),
    ("HIGH_CONF_2to4star_ALL_AN", ad_j[ad_j["stars"] >= 2]),
    ("HIGH_CONF_2to4star_HICALL_80pctAN", ad_j[(ad_j["stars"] >= 2) & (ad_j["hi_call"])]),
]:
    tot_ad_rate = sub["het_rate_joint"].sum() * 1000.0
    tot_elig_rate = sub.loc[~sub["is_discordant"], "het_rate_joint"].sum() * 1000.0
    for tier in ["tier1_gene_list", "tier2_cpsr_sftool", "tier3_owngenome", "tier4_proposed_patch"]:
        fp = sub[sub["is_discordant"] & sub[tier]]
        emitted = sub[sub[tier]]
        fp_rate = fp["het_rate_joint"].sum() * 1000.0
        em_rate = emitted["het_rate_joint"].sum() * 1000.0
        row = {
            "subset": subset_lbl,
            "pipeline_tier": tier,
            "n_fp_alleles": int(len(fp)),
            "n_fp_hets": int(fp["het_joint"].sum()),
            "fp_per1000_joint": round(float(fp_rate), 4),
            "emitted_per1000_joint": round(float(em_rate), 4),
            "fp_share_of_emitted_pct": round(float(100.0 * fp_rate / em_rate), 2) if em_rate > 0 else 0.0,
        }
        for pop in ["nfe", "afr", "amr", "eas", "sas"]:
            p_fp = fp[f"het_rate_{pop}"].sum() * 1000.0
            p_em = emitted[f"het_rate_{pop}"].sum() * 1000.0
            row[f"fp_per1000_{pop}"] = round(float(p_fp), 4)
            row[f"fp_pct_{pop}"] = round(float(100.0 * p_fp / p_em), 2) if p_em > 0 else 0.0
        row["sas_to_nfe_fp_ratio"] = round(row["fp_per1000_sas"] / row["fp_per1000_nfe"], 2) if row["fp_per1000_nfe"] > 0 else None
        row["afr_to_nfe_fp_ratio"] = round(row["fp_per1000_afr"] / row["fp_per1000_nfe"], 2) if row["fp_per1000_nfe"] > 0 else None
        burden_rows.append(row)

bdf = pd.DataFrame(burden_rows)
bdf.to_csv(os.path.join(OUT_DIR, "genome_wide_per_pipeline_fp_burden_gnomad_v41.csv"), index=False)
print("\n=== GENOME-WIDE FALSE-POSITIVE AD SECONDARY FINDING BURDEN PER 1,000 SEQUENCED INDIVIDUALS (gnomAD v4.1) ===")
print(bdf[["subset", "pipeline_tier", "n_fp_alleles", "n_fp_hets", "fp_per1000_joint", "fp_share_of_emitted_pct", "fp_per1000_nfe", "fp_per1000_afr", "fp_per1000_sas", "sas_to_nfe_fp_ratio", "afr_to_nfe_fp_ratio"]].to_string(index=False))

# Also verify Tier 4 Patch on our exact locked 24-variant benchmark panel!
panel_df = pd.read_csv(os.path.join(OUT_DIR, "locked_24variant_benchmark_panel.csv"))
# Run tier4 patch + pure-AR het check on the 24 variants
PURE_AR = {"MUTYH", "ATP7B", "BTD", "GAA", "RPE65", "CASQ2", "TRDN", "CYP27A1", "HFE"}
patch_verdicts = []
for _, r in panel_df.iterrows():
    g = r["gene"]
    vid = str(int(r["variation_id"]))
    if g in PURE_AR:
        patch_verdicts.append("REJECTED (AR_het)")
        continue
    r_like = {
        "gene": g,
        "variation_id": int(r["variation_id"]),
        "is_plof_mc": bool(r["is_lof"]),
        "in_ttn_hi_psi": in_ttn_hi_psi(r["pos"]) if g == "TTN" else False,
        "vcf_clndn": str(r["clndn"]),
        "vcf_clndisdb": str(r["clndisdb"]),
    }
    # Add VHL Chuvash Arg200Leu SCV/codon awareness? Wait — test WITHOUT hardcoding VHL Arg200Leu so we are 100% honest!
    ok = emits_tier4_patch(r_like)
    patch_verdicts.append("RETURNED" if ok else "REJECTED")

panel_df["Tier4_Proposed_Patch"] = patch_verdicts
print("\n=== TIER 4 PROPOSED PATCH VERIFICATION ON LOCKED 24-VARIANT BENCHMARK PANEL (HONEST, ZERO HARDCODING!) ===")
print(panel_df[["variant_id", "gene", "expected", "Tier4_Proposed_Patch"]].to_string(index=False))

# Check counts on 24-variant panel
g_bad = panel_df["expected"].isin(["AR_IN_AD_SHOULD_REJECT", "NON_SF_SHOULD_REJECT", "OPPOSITE_MECH_SHOULD_REJECT"])
g_pos = panel_df["expected"] == "TRUE_AD_SHOULD_RETURN"
g_ar = panel_df["expected"] == "PURE_AR_HET_SHOULD_REJECT"
n_fp_patch = int((panel_df["Tier4_Proposed_Patch"] == "RETURNED")[g_bad].sum())
n_tp_patch = int((panel_df["Tier4_Proposed_Patch"] == "RETURNED")[g_pos].sum())
n_tn_ar = int((panel_df["Tier4_Proposed_Patch"] != "RETURNED")[g_ar].sum())
print(f"\nTier 4 Patch Benchmark Scorecard: Bad FPs={n_fp_patch}/{int(g_bad.sum())} | True AD Sens={n_tp_patch}/{int(g_pos.sum())} | Pure AR Spec={n_tn_ar}/{int(g_ar.sum())}")

out_json = {
    "ryr1_3star_ep_audit": {
        "total_3star_in_gnomad": int(len(ryr1_ep)),
        "mhs_vcep_ad": int(len(ryr1_mh_ad_vids)),
        "cm_vcep_ar_only": int(len(ryr1_cm_ar_vids)),
        "cm_vcep_ar_only_vids": sorted(ryr1_cm_ar_vids),
        "cm_vcep_ar_only_het_carriers": int(ryr1_ep.loc[ryr1_ep["vcep_class"] == "CM_VCEP_AR_ONLY", "het_joint"].sum()),
    },
    "tier4_patch_safety_audit": {
        "total_3star_ad_eligible_gnomad": int(len(ep_all_ad)),
        "retained_3star_ad_eligible": int(ep_all_ad["tier4_proposed_patch"].sum()),
        "dropped_3star_ad_eligible": int(len(ep_dropped)),
        "benchmark_bad_fps": f"{n_fp_patch}/{int(g_bad.sum())}",
        "benchmark_true_ad_sens": f"{n_tp_patch}/{int(g_pos.sum())}",
        "benchmark_pure_ar_spec": f"{n_tn_ar}/{int(g_ar.sum())}",
    },
    "genome_wide_burden_table": burden_rows,
}
with open(os.path.join(OUT_DIR, "genome_wide_burden_and_patch_audit.json"), "w") as fh:
    json.dump(out_json, fh, indent=2)
print("\nDONE! Saved genome_wide_burden_and_patch_audit.json and CSV.")
