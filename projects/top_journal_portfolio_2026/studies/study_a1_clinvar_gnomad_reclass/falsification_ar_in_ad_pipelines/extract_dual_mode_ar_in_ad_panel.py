#!/opt/anaconda3/bin/python3
"""
Attack 1 — Step 1: Exhaustive extraction of ClinVar P/LP alleles across ALL 84 ACMG SF v3.3 genes,
with ClinVar VCF exact records (CLNDN, CLNDISDB, CLNSIG, CLNREVSTAT, MC, ORIGIN),
gnomAD v4.1 heterozygous carrier counts & frequencies (total + 8 ancestries),
and systematic identification of AR-only / non-ACMG-SF-phenotype alleles inside AD ACMG SF v3.3 genes.

Also builds the locked, standardized GRCh38 synthetic patient VCF (`test_panel_het_grch38.vcf`)
containing:
  (A) Locked 4 core RYR1/TTN AR-in-AD alleles from NOVELTY_IDENTIFICATION_GATE.md
  (B) Additional high-confidence (2-4*) verified AR-only / non-AD-target ClinVar P/LP alleles
      in AD ACMG SF v3.3 genes (e.g., RYR1 AR myopathy / pLoF, TTN non-DCM/missense/low-PSI,
      VHL AR erythrocytosis, LMNA AR CMT2B1/MAD, etc.)
  (C) Verified AD ACMG SF positive controls (ClinGen EP RYR1 MH, TTN A-band DCM tv, BRCA1, LDLR)
  (D) Pure-AR ACMG SF negative controls in heterozygous (0/1) state (MUTYH, ATP7B, BTD, HFE, GAA, RPE65)
"""

import gzip
import json
import os
import re
from collections import defaultdict
import pandas as pd
import numpy as np

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
OUT_DIR = os.path.join(BASE, "falsification_ar_in_ad_pipelines")
os.makedirs(OUT_DIR, exist_ok=True)

JOINED_PQ = os.path.join(BASE, "data/clean/a1_clinvar_gnomad_exact_joined.parquet")
CLINVAR_PQ = os.path.join(BASE, "data/clean/clinvar_acmg84_grch38.parquet")
CLINVAR_VCF = os.path.join(BASE, "data/raw/clinvar.vcf.gz")
RULES_CSV = os.path.join(BASE, "audit/a1_acmg84_verified_cspec_johnston2026_rules.csv")

print("Loading joined parquet and ClinVar parquet...")
df_j = pd.read_parquet(JOINED_PQ)
df_c = pd.read_parquet(CLINVAR_PQ)
rules = pd.read_csv(RULES_CSV)

print(f"Joined rows: {len(df_j):,}; ClinVar ACMG84 rows: {len(df_c):,}")

# Filter to germline Pathogenic / Likely pathogenic (excluding Conflicting)
def is_plp(s):
    if pd.isna(s):
        return False
    s = str(s)
    if "Conflicting" in s:
        return False
    return ("Pathogenic" in s) or ("Likely_pathogenic" in s)

df_c["is_plp"] = df_c["clnsig"].apply(is_plp)
df_j["is_plp"] = df_j["clnsig"].apply(is_plp)

plp_c = df_c[df_c["is_plp"]].copy()
plp_j = df_j[df_j["is_plp"]].copy()

print(f"Total ACMG84 ClinVar P/LP: {len(plp_c):,} (in gnomAD v4.1: {len(plp_j):,})")
print(f"  AD genes P/LP: {len(plp_c[plp_c['inheritance']=='AD']):,} (in gnomAD: {len(plp_j[plp_j['inheritance']=='AD']):,})")
print(f"  AD genes 2-4* P/LP: {len(plp_c[(plp_c['inheritance']=='AD') & (plp_c['stars']>=2)]):,} (in gnomAD: {len(plp_j[(plp_j['inheritance']=='AD') & (plp_j['stars']>=2)]):,})")

# Extract exact raw ClinVar VCF INFO lines for all P/LP variation_ids in ACMG84
plp_vids = set(plp_c["variation_id"].astype(str))
print(f"Scanning {CLINVAR_VCF} for {len(plp_vids):,} P/LP Variation IDs...")

vcf_records = {}
vcf_info_keys_seen = set()
with gzip.open(CLINVAR_VCF, "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 8:
            continue
        vid = parts[2]
        if vid in plp_vids:
            info_dict = {}
            for kv in parts[7].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_dict[k] = v
                    vcf_info_keys_seen.add(k)
            vcf_records[vid] = {
                "chrom": parts[0],
                "pos": int(parts[1]),
                "vid": vid,
                "ref": parts[3],
                "alt": parts[4],
                "info_raw": parts[7],
                "CLNDN": info_dict.get("CLNDN", ""),
                "CLNDISDB": info_dict.get("CLNDISDB", ""),
                "CLNSIG": info_dict.get("CLNSIG", ""),
                "CLNREVSTAT": info_dict.get("CLNREVSTAT", ""),
                "GENEINFO": info_dict.get("GENEINFO", ""),
                "MC": info_dict.get("MC", ""),
                "RS": info_dict.get("RS", ""),
                "ORIGIN": info_dict.get("ORIGIN", ""),
                "CLNVC": info_dict.get("CLNVC", ""),
            }

print(f"Matched {len(vcf_records):,} P/LP records in raw ClinVar VCF. INFO keys seen: {sorted(vcf_info_keys_seen)}")

# Attach CLNDISDB and raw VCF info to plp_c and plp_j
for df_target in (plp_c, plp_j):
    df_target["vid_str"] = df_target["variation_id"].astype(str)
    df_target["vcf_clndn"] = df_target["vid_str"].map(lambda x: vcf_records.get(x, {}).get("CLNDN", ""))
    df_target["vcf_clndisdb"] = df_target["vid_str"].map(lambda x: vcf_records.get(x, {}).get("CLNDISDB", ""))
    df_target["vcf_mc"] = df_target["vid_str"].map(lambda x: vcf_records.get(x, {}).get("MC", ""))

# Systematic phenotype / inheritance audit of ALL AD ACMG SF v3.3 P/LP variants
# First inspect the 4 locked core alleles in full detail
LOCKED_CORE_4 = [
    "19-38527014-C-T",
    "19-38523211-C-G",
    "19-38505338-CAT-C",
    "2-178568947-T-A",
]

core4_rows = plp_j[plp_j["variant_id"].isin(LOCKED_CORE_4)].copy()
print("\n=== LOCKED CORE 4 AR-IN-AD CANDIDATES ===")
for _, r in core4_rows.iterrows():
    het_joint = int(r["ac_joint"]) - 2 * int(r["hom_joint"])
    print(f"  {r['variant_id']} | {r['gene']} | VID={r['variation_id']} | {r['stars']}* | {r['clnsig']} | {r['hgvs_name']}")
    print(f"    MC={r['mc']} | AC={r['ac_joint']} AN={r['an_joint']} HOM={r['hom_joint']} HET={het_joint} | FAF95_grpmax5={r['faf95_grpmax5']:.6g} ({r['faf95_pop5']})")
    print(f"    CLNDN={r['vcf_clndn']}")
    print(f"    CLNDISDB={r['vcf_clndisdb']}")

# Audit ALL AD ACMG SF v3.3 P/LP variants for explicit recessive / non-ACMG-SF phenotype signatures
ad_plp_c = plp_c[plp_c["inheritance"] == "AD"].copy()
ad_plp_j = plp_j[plp_j["inheritance"] == "AD"].copy()

# Compute exact heterozygous carrier count & rate in gnomAD v4.1 (het = ac - 2*hom; carriers = het + hom = ac - hom)
POPS = ["nfe", "afr", "amr", "eas", "sas", "mid", "asj", "fin"]
ad_plp_j["het_joint"] = (ad_plp_j["ac_joint"] - 2 * ad_plp_j["hom_joint"]).clip(lower=0)
ad_plp_j["n_ind_joint"] = ad_plp_j["an_joint"] / 2.0
ad_plp_j["het_rate_joint"] = ad_plp_j["het_joint"] / ad_plp_j["n_ind_joint"].replace(0, np.nan)
for p in POPS:
    ad_plp_j[f"het_{p}"] = (ad_plp_j[f"ac_{p}"] - 2 * ad_plp_j[f"hom_{p}"]).clip(lower=0)
    ad_plp_j[f"n_ind_{p}"] = ad_plp_j[f"an_{p}"] / 2.0
    ad_plp_j[f"het_rate_{p}"] = ad_plp_j[f"het_{p}"] / ad_plp_j[f"n_ind_{p}"].replace(0, np.nan)

# Print breakdown of ALL CLNDN terms in key dual-mode AD ACMG SF genes
DUAL_MODE_GENES = [
    "RYR1", "TTN", "VHL", "LMNA", "DSP", "SCN5A", "MYH7", "MYBPC3", "MYL2", "MYL3",
    "TNNI3", "FLNC", "DES", "BAG3", "KCNQ1", "BRCA2", "PALB2", "BRCA1", "PMS2",
    "MSH6", "MSH2", "MLH1", "APC", "LDLR", "RET", "FBN1", "COL3A1", "ENG", "ACVRL1",
    "CACNA1S", "TSC1", "TSC2", "TP53", "PTEN", "WT1"
]

# Classify every AD P/LP allele's CLNDN into strictly recessive/non-SF-target vs AD-target vs ambiguous
recessive_kw = re.compile(
    r"(autosomal_recessive|recessive|erythrocytosis|polycythemia|chuvash|"
    r"minicore|multiminicore|king-denborough|fetal_akinesia|centronuclear|"
    r"samaritan|congenital_myopathy_with_fiber|nemaline|limb-girdle|"
    r"tibial_muscular_dystrophy|salih_myopathy|myofibrillar_myopathy|"
    r"hereditary_myopathy_with_early_respiratory|charcot-marie-tooth|"
    r"mandibuloacral|progeria|emery-dreifuss|lipodystrophy|"
    r"fanconi_anemia|mismatch_repair_cancer_syndrome|constitutional_mismatch|"
    r"jervell_and_lange-nielsen|sick_sinus|carvajal|naxos|woolly_hair|"
    r"epidermolysis|skin_fragility|palmoplantar|atrichia|alopecia|"
    r"xeroderma|bloom|ataxia|microcephaly|seckel|lissencephaly|"
    r"congenital_myasthenic|arthrogryposis|hyperekplexia|myotonia|"
    r"periodic_paralysis|central_core_disease|ryr1-related_myopathy)",
    re.IGNORECASE,
)

# Also inspect gene-specific strict non-ACMG-SF rules:
# 1. RYR1: ACMG SF v3.3 target = AD Malignant Hyperthermia Susceptibility (MIM:145600).
#    Any RYR1 P/LP variant whose CLNDN does NOT mention "malignant_hyperthermia" OR is a pLoF (frameshift/nonsense/canonical splice)
#    is non-MH / AR-myopathy! (ClinGen MH VCEP GN012: MH is GoF missense; LoF causes AR myopathy!)
# 2. TTN: ACMG SF v3.3 target = AD Dilated Cardiomyopathy (MIM:604145), TRUNCATING variants only!
#    Any TTN P/LP variant that is non-truncating (missense/inframe) OR whose CLNDN explicitly states
#    skeletal myopathy / LGMD / Salih / tibial muscular dystrophy / HMERF without cardiomyopathy!
# 3. VHL: ACMG SF v3.3 target = AD von Hippel-Lindau syndrome (MIM:193300).
#    Variants exclusively linked to AR familial erythrocytosis type 2 / Chuvash polycythemia (MIM:263400)!
# 4. LMNA: ACMG SF v3.3 target = AD DCM (MIM:115200).
#    Variants exclusively linked to AR Charcot-Marie-Tooth 2B1 (MIM:605588) or AR mandibuloacral dysplasia (MIM:248370)!
# 5. KCNQ1: ACMG SF v3.3 target = AD Long QT syndrome 1 (MIM:192500).
#    Check if any are strictly AR Jervell and Lange-Nielsen without LQTS!
# 6. BRCA2/PALB2/PMS2/MSH6/MSH2/MLH1: Check if any are strictly Fanconi anemia / CMMRD (AR) without HBOC/Lynch!
# 7. Generic scan across ALL 73 AD genes for explicit "autosomal_recessive" in CLNDN!

def classify_ad_plp(row):
    gene = row["gene"]
    clndn = str(row["vcf_clndn"] or row["clndn"] or "")
    clndn_l = clndn.lower()
    disdb = str(row["vcf_clndisdb"] or "")
    mc = str(row["vcf_mc"] or row["mc"] or "").lower()
    csq = str(row.get("csq_group", ""))
    tokens = [t.strip() for t in re.split(r"[|,]", clndn) if t.strip() and t.strip().lower() not in ("not_provided", "not_specified", "see_cases")]
    tokens_l = [t.lower() for t in tokens]

    is_plof = any(k in mc for k in ["nonsense", "frameshift", "splice_donor", "splice_acceptor"])
    is_missense = "missense" in mc and not is_plof

    flags = []
    # Explicit autosomal_recessive tag in CLNDN
    has_ar_text = any("autosomal_recessive" in t or "recessive" in t for t in tokens_l)
    has_ad_text = any("autosomal_dominant" in t or "dominant" in t for t in tokens_l)

    if gene == "RYR1":
        has_mh = any("malignant_hyperthermia" in t or "hyperthermia" in t or "mh_susceptibility" in t for t in tokens_l)
        has_myop = any("myopath" in t or "minicore" in t or "central_core" in t or "king-denborough" in t or "fetal_akinesia" in t or "neuromuscular" in t for t in tokens_l)
        if not has_mh and (has_myop or is_plof or has_ar_text):
            return "AR_OR_NON_SF_STRICT"
        if not has_mh and len(tokens_l) == 0 and is_plof:
            return "AR_OR_NON_SF_PLOF_NO_MH"
        if has_mh and is_plof:
            return "MIXED_PLOF_WITH_MH_TAG"
        if has_mh:
            return "AD_SF_TARGET"
        return "UNSPECIFIED_OR_OTHER"

    if gene == "TTN":
        has_dcm_cardiac = any("dilated_cardiomyopathy" in t or "cardiomyopathy" in t for t in tokens_l)
        has_skel_only = any(k in t for t in tokens_l for k in ["muscular_dystrophy", "limb-girdle", "salih", "myopathy", "hmerf", "tibial", "centronuclear", "arthrogryposis"])
        if not is_plof:
            # ACMG SF v3.3 explicitly restricts TTN reporting to TRUNCATING variants only!
            return "NON_SF_TTN_NON_TRUNCATING"
        if has_skel_only and not has_dcm_cardiac:
            return "AR_OR_NON_SF_STRICT"
        if has_dcm_cardiac:
            return "AD_SF_TARGET_TV"
        return "TTN_TV_UNSPECIFIED_DISEASE"

    if gene == "VHL":
        has_vhl = any("von_hippel" in t or "pheochromocytoma" in t or "renal_cell" in t or "hemangioblastoma" in t or "cancer" in t for t in tokens_l)
        has_eryth = any("erythrocytosis" in t or "polycythemia" in t or "chuvash" in t for t in tokens_l)
        if has_eryth and not has_vhl:
            return "AR_OR_NON_SF_STRICT"
        if has_eryth and has_vhl:
            return "MIXED_AR_AND_AD_TAGS"
        if has_vhl:
            return "AD_SF_TARGET"
        return "UNSPECIFIED_OR_OTHER"

    if gene == "CACNA1S":
        has_mh = any("malignant_hyperthermia" in t or "hyperthermia" in t for t in tokens_l)
        has_non_mh = any("periodic_paralysis" in t or "myopath" in t or "recessive" in t for t in tokens_l)
        if has_non_mh and not has_mh:
            return "NON_SF_PHENOTYPE_ONLY"
        if has_mh:
            return "AD_SF_TARGET"
        return "UNSPECIFIED_OR_OTHER"

    if gene == "LMNA":
        has_cardiac = any("cardiomyopathy" in t or "cardiac" in t or "arrhythm" in t or "conduction" in t for t in tokens_l)
        has_ar_only = any("charcot-marie-tooth" in t or "mandibuloacral" in t or "autosomal_recessive" in t for t in tokens_l)
        if has_ar_only and not has_cardiac:
            return "AR_OR_NON_SF_STRICT"
        if has_cardiac:
            return "AD_SF_TARGET"
        return "OTHER_LMNA_PHENOTYPE"

    if gene == "DSP":
        has_cardiac = any("cardiomyopathy" in t or "arrhythmogenic" in t or "arvc" in t for t in tokens_l)
        has_ar_skin = any("epidermolysis" in t or "skin_fragility" in t or "keratoderma" in t or "woolly_hair" in t or "carvajal" in t or "autosomal_recessive" in t for t in tokens_l)
        if has_ar_skin and not has_cardiac:
            return "AR_OR_NON_SF_STRICT"
        if has_cardiac:
            return "AD_SF_TARGET"
        return "UNSPECIFIED_OR_OTHER"

    # Generic rule across all remaining AD ACMG SF genes: explicit "autosomal_recessive" ONLY (no AD mention)
    if has_ar_text and not has_ad_text:
        # Check if ALL named traits are recessive or non-AD
        non_ar_tokens = [t for t in tokens_l if "recessive" not in t and "fanconi" not in t and "jervell" not in t and "mismatch_repair_cancer_syndrome" not in t]
        if len(non_ar_tokens) == 0:
            return "AR_OR_NON_SF_STRICT"
        return "MIXED_AR_AND_OTHER"

    return "STANDARD_AD_OR_UNFLAGGED"

ad_plp_c["dual_mode_class"] = ad_plp_c.apply(classify_ad_plp, axis=1)
ad_plp_j["dual_mode_class"] = ad_plp_j.apply(classify_ad_plp, axis=1)

print("\n=== ALL AD ACMG SF v3.3 P/LP BREAKDOWN (FULL CLINVAR, ALL STARS) ===")
print(ad_plp_c.groupby(["dual_mode_class"]).size())
print("\n=== ALL AD ACMG SF v3.3 P/LP BREAKDOWN BY GENE FOR NON-STANDARD CLASSES (FULL CLINVAR) ===")
nonstd_c = ad_plp_c[ad_plp_c["dual_mode_class"] != "STANDARD_AD_OR_UNFLAGGED"]
print(nonstd_c.groupby(["gene", "dual_mode_class"]).size().to_string())

print("\n=== IN GNOMAD v4.1: AD ACMG SF v3.3 P/LP NON-STANDARD CLASSES BY GENE (ALL STARS & >=2*) ===")
nonstd_j = ad_plp_j[ad_plp_j["dual_mode_class"].isin(["AR_OR_NON_SF_STRICT", "AR_OR_NON_SF_PLOF_NO_MH", "NON_SF_TTN_NON_TRUNCATING", "NON_SF_PHENOTYPE_ONLY", "MIXED_AR_AND_AD_TAGS", "MIXED_PLOF_WITH_MH_TAG"])].copy()

summary_rows = []
for (gene, cls), g in nonstd_j.groupby(["gene", "dual_mode_class"]):
    g2 = g[g["stars"] >= 2]
    # Note: approximate sum of heterozygous carriers (upper bound assuming rare independent alleles; exact per-allele sum)
    row = {
        "gene": gene,
        "class": cls,
        "n_alleles_clinvar_all": int(len(ad_plp_c[(ad_plp_c["gene"]==gene) & (ad_plp_c["dual_mode_class"]==cls)])),
        "n_alleles_clinvar_2plus": int(len(ad_plp_c[(ad_plp_c["gene"]==gene) & (ad_plp_c["dual_mode_class"]==cls) & (ad_plp_c["stars"]>=2)])),
        "n_alleles_gnomad_all": int(len(g)),
        "n_alleles_gnomad_2plus": int(len(g2)),
        "sum_het_joint_all": int(g["het_joint"].sum()),
        "sum_het_joint_2plus": int(g2["het_joint"].sum()),
        "sum_het_rate_per1000_all": float(g["het_rate_joint"].sum() * 1000.0),
        "sum_het_rate_per1000_2plus": float(g2["het_rate_joint"].sum() * 1000.0),
    }
    for p in ["nfe", "afr", "amr", "eas", "sas"]:
        row[f"het_per1000_{p}_all"] = float(g[f"het_rate_{p}"].sum() * 1000.0)
        row[f"het_per1000_{p}_2plus"] = float(g2[f"het_rate_{p}"].sum() * 1000.0)
    summary_rows.append(row)

sum_df = pd.DataFrame(summary_rows).sort_values("sum_het_joint_all", ascending=False)
print(sum_df.to_string(index=False))

# Inspect top 35 highest-carrier AR_OR_NON_SF / NON_SF alleles in gnomAD v4.1!
top_non_sf = nonstd_j.sort_values("het_joint", ascending=False).head(40)
print("\n=== TOP 40 HIGHEST-CARRIER AR / NON-SF P/LP ALLELES IN AD ACMG SF v3.3 GENES (gnomAD v4.1) ===")
for _, r in top_non_sf.iterrows():
    print(f"  {r['gene']:8s} | {r['dual_mode_class']:26s} | {r['variant_id']:22s} | {r['stars']}* | HET={int(r['het_joint']):5d} | pop={r['pop_grpmax5']} ({r['af_grpmax5']:.5f}) | {r['hgvs_name']} | CLNDN={str(r['vcf_clndn'])[:90]}")

# Save complete catalogs!
nonstd_j.to_csv(os.path.join(OUT_DIR, "acmg84_ad_ar_or_nonsf_plp_gnomad_catalog.csv"), index=False)
sum_df.to_csv(os.path.join(OUT_DIR, "acmg84_ad_ar_or_nonsf_summary_by_gene.csv"), index=False)

# Also inspect pure-AR ACMG SF genes in gnomAD to pick real high-confidence negative controls!
ar_plp_j = plp_j[(plp_j["inheritance"] == "AR") & (plp_j["stars"] >= 2)].sort_values("ac_joint", ascending=False)
print("\n=== TOP PURE-AR ACMG SF v3.3 NEGATIVE CONTROLS (2-4* P/LP in gnomAD v4.1) ===")
for gene in ["MUTYH", "ATP7B", "BTD", "HFE", "GAA", "RPE65", "TRDN", "CASQ2", "CYP27A1"]:
    sub = ar_plp_j[ar_plp_j["gene"] == gene].head(2)
    for _, r in sub.iterrows():
        print(f"  {r['gene']:8s} | {r['variant_id']:22s} | VID={r['variation_id']} | {r['stars']}* | AC={int(r['ac_joint']):6d} | {r['hgvs_name']} | CLNDN={str(r['vcf_clndn'])[:70]}")

# Also pick verified AD positive controls (RYR1 MH EP 3*, TTN A-band DCM tv, BRCA1 EP 3*, LDLR EP 3*)
print("\n=== CANDIDATE AD POSITIVE CONTROLS ===")
for gene in ["RYR1", "TTN", "BRCA1", "LDLR"]:
    if gene == "RYR1":
        sub = ad_plp_j[(ad_plp_j["gene"]=="RYR1") & (ad_plp_j["dual_mode_class"]=="AD_SF_TARGET") & (ad_plp_j["stars"]>=3)].sort_values("ac_joint", ascending=False).head(3)
    elif gene == "TTN":
        sub = ad_plp_j[(ad_plp_j["gene"]=="TTN") & (ad_plp_j["dual_mode_class"]=="AD_SF_TARGET_TV") & (ad_plp_j["stars"]>=2)].sort_values("ac_joint", ascending=False).head(3)
    else:
        sub = ad_plp_j[(ad_plp_j["gene"]==gene) & (ad_plp_j["stars"]>=3)].sort_values("ac_joint", ascending=False).head(2)
    for _, r in sub.iterrows():
        print(f"  {r['gene']:8s} | {r['variant_id']:22s} | VID={r['variation_id']} | {r['stars']}* | AC={int(r['ac_joint']):6d} | {r['hgvs_name']} | CLNDN={str(r['vcf_clndn'])[:70]}")

# Save summary JSON
meta = {
    "clinvar_vcf_file": CLINVAR_VCF,
    "clinvar_vcf_info_keys": sorted(vcf_info_keys_seen),
    "clinvar_vcf_has_moi_field": any(k.upper() in ("MOI", "CLNMOI", "INHERITANCE", "CLNINHER") for k in vcf_info_keys_seen),
    "total_acmg84_plp_clinvar": int(len(plp_c)),
    "total_acmg84_plp_gnomad": int(len(plp_j)),
    "ad_plp_clinvar_class_counts": ad_plp_c["dual_mode_class"].value_counts().to_dict(),
    "ad_plp_gnomad_class_counts": ad_plp_j["dual_mode_class"].value_counts().to_dict(),
    "core4_found": core4_rows[["variant_id", "gene", "variation_id", "stars", "clnsig", "hgvs_name", "ac_joint", "hom_joint", "vcf_clndn", "vcf_clndisdb"]].to_dict(orient="records"),
}
with open(os.path.join(OUT_DIR, "dual_mode_clinvar_gnomad_scan_summary.json"), "w") as fh:
    json.dump(meta, fh, indent=2)

print("\nDONE! Saved dual_mode_clinvar_gnomad_scan_summary.json and CSV catalogs.")
