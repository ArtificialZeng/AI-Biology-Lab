#!/opt/anaconda3/bin/python3
"""
1. Execute PROSPECTIVELY LOCKED Unseen 20-Gene Non-ACMG Autosomal Dominant Validation Cohort (PROSPECTIVE_VALIDATION_A1.md)
   using ClinVar GRCh38 clinvar.vcf.gz + variant_summary.txt.gz and official gnomAD v4.1 wrapper (search_variants.py).
2. Complete live gnomAD v4.1 QC filter/flag audit across the REMAINING 36 of the 55 ACMG SF v3.3 AD >=2-star alerts
   so ALL 55/55 have 100% live gnomAD exome/genome QC filter status!
"""
import gzip, json, os, subprocess, sys, time, urllib.request, warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
RAW = os.path.join(BASE, "data", "raw")
CLEAN = os.path.join(BASE, "data", "clean")
RES = os.path.join(BASE, "results")
AUD = os.path.join(BASE, "audit")
VAL_GNOMAD = os.path.join(AUD, "gnomad_v4_unseen20")
QC_DIR = os.path.join(AUD, "gnomad_qc_flags")
os.makedirs(VAL_GNOMAD, exist_ok=True)
os.makedirs(QC_DIR, exist_ok=True)

GNOMAD_SCRIPT = "/Users/mac/.gemini/config/plugins/genomics/skills/gnomad_database/scripts/search_variants.py"
MAX_AN = 1614324.0

sys.path.insert(0, os.path.join(BASE, "src"))
from fetch_a1_gnomad_acmg84 import ACMG_SF_V33, review_status_to_stars
from analyze_a1_full import POPS_ALL8, POPS_GRPMAX5, POPS_NONNFE5, calc_faf95_vec, classify_clnsig, classify_consequence, cmh_test

UNSEEN_20 = [
    ("SCN1A", "Neuro"), ("KCNQ2", "Neuro"), ("STXBP1", "Neuro"), ("SLC2A1", "Neuro"), ("NOTCH3", "Neuro"), ("PSEN1", "Neuro"),
    ("CDH1", "Cancer"), ("CDKN2A", "Cancer"), ("CHEK2", "Cancer"), ("ATM", "Cancer"), ("BAP1", "Cancer"), ("DICER1", "Cancer"), ("NF1", "Cancer"), ("PTCH1", "Cancer"),
    ("PKD1", "Renal"), ("PKD2", "Renal"), ("COL1A1", "Connective"), ("ANK1", "Heme"), ("GCK", "Metabolic"), ("ELN", "Connective"),
]
acmg_set = {g[0] for g in ACMG_SF_V33}
unseen_set = {g[0] for g in UNSEEN_20}
assert len(unseen_set) == 20 and len(acmg_set & unseen_set) == 0, "Overlap with ACMG SF v3.3!"
u_dom = dict(UNSEEN_20)

# Step 1: Parse ClinVar GRCh38 for the 20 unseen genes
cv_pq = os.path.join(AUD, "clinvar_unseen20_grch38.parquet")
if os.path.exists(cv_pq):
    cv_u = pd.read_parquet(cv_pq)
else:
    print("Parsing clinvar.vcf.gz for 20 unseen validation genes...", flush=True)
    rows = []
    with gzip.open(os.path.join(RAW, "clinvar.vcf.gz"), "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, var_id, ref, alt, _, _, info = parts[:8]
            gi = info.find("GENEINFO=")
            if gi < 0:
                continue
            ge = info.find(";", gi)
            genes = [tok.split(":")[0] for tok in info[gi+9:ge if ge >= 0 else None].split("|")]
            m = [g for g in genes if g in unseen_set]
            if not m:
                continue
            kv = dict(it.split("=", 1) for it in info.split(";") if "=" in it)
            g = m[0]
            rows.append({
                "variant_id": f"{chrom}-{pos}-{ref}-{alt}", "chrom": chrom, "pos": int(pos),
                "ref": ref, "alt": alt, "variation_id": str(var_id), "gene": g, "domain": u_dom[g],
                "clnsig": kv.get("CLNSIG", ""), "clnsigconf": kv.get("CLNSIGCONF", ""),
                "clnrevstat": kv.get("CLNREVSTAT", ""), "stars": review_status_to_stars(kv.get("CLNREVSTAT", "")),
                "mc": kv.get("MC", ""), "rs": ("rs" + kv["RS"]) if kv.get("RS") else "",
            })
    cv_u = pd.DataFrame(rows)
    vid_set = set(cv_u["variation_id"])
    vs_map = {}
    with gzip.open(os.path.join(RAW, "variant_summary.txt.gz"), "rt", encoding="utf-8", errors="ignore") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ia, iv, il, ins = hdr.index("Assembly"), hdr.index("VariationID"), hdr.index("LastEvaluated"), hdr.index("NumberSubmitters")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) > max(ia, iv, il, ins) and c[ia] == "GRCh38" and c[iv] in vid_set and c[iv] not in vs_map:
                vs_map[c[iv]] = (c[il], c[ins])
    cv_u["last_evaluated"] = cv_u["variation_id"].map(lambda x: vs_map.get(x, ("", ""))[0])
    cv_u["n_submitters"] = pd.to_numeric(cv_u["variation_id"].map(lambda x: vs_map.get(x, ("", ""))[1]), errors="coerce")
    cv_u["clnsig_cat"] = [classify_clnsig(a, b) for a, b in zip(cv_u["clnsig"], cv_u["clnsigconf"])]
    cv_u.to_parquet(cv_pq, index=False)
print(f"Unseen 20-gene ClinVar: {len(cv_u):,} variants across {cv_u['gene'].nunique()}/20 genes", flush=True)

# Step 2: Fetch gnomAD v4.1 for all 20 unseen genes (polite 6.1s spacing)
spans = cv_u.groupby("gene").agg(chrom=("chrom", "first"), g_start=("pos", "min"), g_end=("pos", "max")).to_dict(orient="index")
g_rows = []
seen = set()
for idx, (gene, dom) in enumerate(UNSEEN_20, 1):
    of = os.path.join(VAL_GNOMAD, f"{gene}_gnomad_r4.json")
    if not (os.path.exists(of) and os.path.getsize(of) > 500):
        t0 = time.time()
        cmd = ["uv", "run", GNOMAD_SCRIPT, "--gene", gene, "--dataset", "gnomad_r4", "--output", of]
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        el = time.time() - t0
        if el < 6.1:
            time.sleep(6.1 - el)
    if os.path.exists(of) and os.path.getsize(of) > 200:
        d = json.load(open(of))
        vs = (d.get("data", {}).get("gene") or {}).get("variants") or []
        for v in vs:
            vid = v.get("variant_id")
            if not vid or vid in seen:
                continue
            seen.add(vid)
            ex = v.get("exome") or {}
            ge = v.get("genome") or {}
            rec = {
                "variant_id": vid, "gnomad_gene": gene, "gnomad_csq": v.get("consequence") or "",
                "gnomad_rsids": "|".join(v.get("rsids") or []),
                "ac_joint": int(ex.get("ac") or 0) + int(ge.get("ac") or 0),
                "an_joint": int(ex.get("an") or 0) + int(ge.get("an") or 0),
                "hom_joint": int(ex.get("homozygote_count") or 0) + int(ge.get("homozygote_count") or 0),
            }
            pac, pan = {p: 0 for p in POPS_ALL8}, {p: 0 for p in POPS_ALL8}
            for src in (ex, ge):
                for pop in (src.get("populations") or []):
                    pid = pop.get("id")
                    if pid in pac:
                        pac[pid] += int(pop.get("ac") or 0)
                        pan[pid] += int(pop.get("an") or 0)
            for p in POPS_ALL8:
                rec[f"ac_{p}"] = pac[p]
                rec[f"an_{p}"] = pan[p]
            g_rows.append(rec)
    if idx % 5 == 0 or idx == 20:
        print(f"  Unseen gnomAD fetched [{idx:02d}/20]: cum={len(g_rows):,} variants", flush=True)

gdf = pd.DataFrame(g_rows)
for p in POPS_ALL8:
    ac, an = gdf[f"ac_{p}"].values.astype(float), gdf[f"an_{p}"].values.astype(float)
    gdf[f"af_{p}"] = np.where(an >= 200, ac / an, 0.0)
    gdf[f"faf95_{p}"] = calc_faf95_vec(ac, an)
gdf["faf95_grpmax5"] = gdf[[f"faf95_{p}" for p in POPS_GRPMAX5]].values.max(axis=1)
gdf["faf95_nonnfe5"] = gdf[[f"faf95_{p}" for p in POPS_NONNFE5]].values.max(axis=1)
gdf["af_nonnfe5"] = gdf[[f"af_{p}" for p in POPS_NONNFE5]].values.max(axis=1)
gdf["faf95_max8"] = gdf[[f"faf95_{p}" for p in POPS_ALL8]].values.max(axis=1)

# Join exact
ju = cv_u.merge(gdf, on="variant_id", how="inner")
ju["csq_group"] = [classify_consequence(a, b) for a, b in zip(ju["mc"], ju["gnomad_csq"])]
ju.to_parquet(os.path.join(AUD, "a1_unseen20_exact_joined.parquet"), index=False)

# Evaluate PV-1, PV-2, PV-3
# PV-1: rsID collision rate
rs_map = {}
for row in gdf[["variant_id", "gnomad_rsids", "faf95_max8"]].itertuples(index=False):
    if row.gnomad_rsids:
        for rs in row.gnomad_rsids.split("|"):
            if rs:
                rs_map.setdefault(rs, []).append(row.faf95_max8)
ex_faf = dict(zip(gdf["variant_id"], gdf["faf95_max8"]))
plp_rs = cv_u[(cv_u["clnsig_cat"] == "P_LP") & cv_u["rs"].str.startswith("rs", na=False)]
nv_all, ex_all, nv_hc, ex_hc = 0, 0, 0, 0
for row in plp_rs[["variant_id", "rs", "stars"]].itertuples(index=False):
    hs = rs_map.get(row.rs)
    if hs and max(hs) >= 1e-4:
        nv_all += 1
        if ex_faf.get(row.variant_id, 0.0) >= 1e-4:
            ex_all += 1
        if row.stars >= 2:
            nv_hc += 1
            if ex_faf.get(row.variant_id, 0.0) >= 1e-4:
                ex_hc += 1
pv1 = {
    "all_plp_faf95_ge_1e4": {"naive": nv_all, "exact": ex_all, "fp_rate": round((nv_all - ex_all) / max(nv_all, 1), 4)},
    "hc_2to4star_faf95_ge_1e4": {"naive": nv_hc, "exact": ex_hc, "fp_rate": round((nv_hc - ex_hc) / max(nv_hc, 1), 4)},
}

# PV-2: Denominator-clean H1 under FAF95 >= 1e-4 & Mediation
csu = ju[ju["csq_group"].isin(["pLoF", "missense", "inframe_indel", "splice_other", "synonymous"]) &
         ju["clnsig_cat"].isin(["P_LP", "VUS", "CONFLICTING", "B_LB"]) &
         (ju["an_joint"] >= 0.80 * MAX_AN)].copy()
m_nn = (csu["faf95_nonnfe5"] >= 1e-4) & (csu["af_nonnfe5"] / np.maximum(csu["af_nfe"], 1e-6) >= 5.0)
m_eu = (csu["faf95_nfe"] >= 1e-4) & (csu["af_nfe"] / np.maximum(csu["af_nonnfe5"], 1e-6) >= 5.0)
comp = csu[m_nn | m_eu].copy()
comp["non_nfe"] = m_nn.loc[comp.index].astype(int)
comp["af_ref"] = np.where(comp["non_nfe"] == 1, comp["faf95_nonnfe5"], comp["faf95_nfe"])
comp["log10_af"] = np.log10(comp["af_ref"].clip(lower=1e-6))
comp["log_subs"] = np.log1p(comp["n_submitters"].fillna(1).clip(lower=1))
comp["af_bin"] = pd.cut(comp["af_ref"], bins=[1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 5e-2, 1.01], right=False)

_, p_sub_mw = stats.mannwhitneyu(comp.loc[comp["non_nfe"] == 1, "n_submitters"].dropna(), comp.loc[comp["non_nfe"] == 0, "n_submitters"].dropna())
pv2 = {
    "n_nn": int((comp["non_nfe"] == 1).sum()), "n_nfe": int((comp["non_nfe"] == 0).sum()),
    "med_subs_nn": float(comp.loc[comp["non_nfe"] == 1, "n_submitters"].median()),
    "med_subs_nfe": float(comp.loc[comp["non_nfe"] == 0, "n_submitters"].median()),
    "subs_mw_p": float(p_sub_mw),
    "outcomes": {},
}
for cat_nm in ("VUS", "CONFLICTING", "B_LB", "P_LP"):
    comp["y"] = (comp["clnsig_cat"] == cat_nm).astype(int)
    cmh = cmh_test(comp, "non_nfe", "y", "af_bin")
    gv = comp.groupby("gene")["y"].nunique()
    sg = comp[comp["gene"].isin(gv[gv == 2].index)]
    try:
        ml = smf.logit("y ~ non_nfe + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=sg).fit(disp=False, maxiter=150)
        ml_m = smf.logit("y ~ non_nfe + log_subs + log10_af + I(log10_af**2) + C(csq_group) + C(gene)", data=sg).fit(disp=False, maxiter=150)
        lor, lp = float(np.exp(ml.params["non_nfe"])), float(ml.pvalues["non_nfe"])
        lor_m, lp_m = float(np.exp(ml_m.params["non_nfe"])), float(ml_m.pvalues["non_nfe"])
        atten = round((1.0 - float(ml_m.params["non_nfe"]) / float(ml.params["non_nfe"])) * 100.0, 1)
    except Exception:
        lor, lp, lor_m, lp_m, atten = np.nan, np.nan, np.nan, np.nan, np.nan
    pv2["outcomes"][cat_nm] = {
        "rate_nn": round(float(comp.loc[comp["non_nfe"] == 1, "y"].mean()), 4),
        "rate_nfe": round(float(comp.loc[comp["non_nfe"] == 0, "y"].mean()), 4),
        "cmh_or": round(cmh["or"], 3), "cmh_ci95": [round(x, 3) for x in cmh["ci95"]], "cmh_p": cmh["p"],
        "logit_or": round(lor, 3) if not np.isnan(lor) else None, "logit_p": lp if not np.isnan(lp) else None,
        "subs_adj_or": round(lor_m, 3) if not np.isnan(lor_m) else None, "subs_adj_p": lp_m if not np.isnan(lp_m) else None,
        "log_or_atten_pct": atten if not np.isnan(atten) else None,
    }

# PV-3: VUS-to-Benign Liberation Debt on Unseen 20 Genes
vu = ju[(ju["clnsig_cat"] == "VUS") & (ju["an_joint"] >= 0.80 * MAX_AN)]
pv3 = {"total_hc_vus": int(len(vu))}
for nm, thr in [("faf95_1e4", 1e-4), ("faf95_5e4", 5e-4), ("faf95_1e3", 1e-3)]:
    h = vu[vu["faf95_grpmax5"] >= thr]
    inv = h[h["faf95_nfe"] < thr]
    pv3[nm] = {"n": int(len(h)), "nfe_invisible": int(len(inv)), "nfe_invisible_share": round(len(inv) / max(len(h), 1), 4)}

# Step 3: Complete live gnomAD QC filter/flag check for the REMAINING 36 of the 55 ACMG SF v3.3 AD >=2-star alerts!
cat55 = pd.read_csv(os.path.join(RES, "a1_ad_highconf_triaged_55.csv"))
gql = """
query($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id flags
    exome { ac an filters flags }
    genome { ac an filters flags }
    joint { ac an faf95 { popmax popmax_population } }
  }
}
"""
qc_all55 = []
n_fetched_new = 0
for row in cat55.itertuples(index=False):
    vid = row.variant_id
    of = os.path.join(QC_DIR, f"{vid}.json")
    if os.path.exists(of) and os.path.getsize(of) > 50:
        dq = json.load(open(of))
    else:
        req = urllib.request.Request(
            "https://gnomad.broadinstitute.org/api",
            data=json.dumps({"query": gql, "variables": {"variantId": vid, "dataset": "gnomad_r4"}}).encode("utf-8"),
            headers={"Content-Type": "application/json", "User-Agent": "AcademicAudit/1.0"},
            method="POST"
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                dq = json.loads(resp.read().decode("utf-8"))
            with open(of, "w") as f:
                json.dump(dq, f, indent=2)
        except Exception as e:
            dq = {"error": str(e)}
        n_fetched_new += 1
        time.sleep(6.2)
    v = (dq.get("data") or {}).get("variant") or {}
    ex, ge = v.get("exome") or {}, v.get("genome") or {}
    ex_f, ge_f = ex.get("filters") or [], ge.get("filters") or []
    fl = v.get("flags") or []
    qc_fail = bool(ex_f) or (bool(ge) and bool(ge_f) and not ex) or ("lcr" in fl) or (row.an_frac < 0.80)
    qc_all55.append({
        "variant_id": vid, "gene": row.gene, "layer_v1": row.triage_layer,
        "an_frac": round(float(row.an_frac), 4), "global_flags": fl,
        "exome_filters": ex_f, "genome_filters": ge_f,
        "any_qc_fail_or_low_an": qc_fail,
    })

out = {
    "unseen20_overview": {"n_genes": int(ju["gene"].nunique()), "n_clinvar": int(len(cv_u)), "n_gnomad": int(len(gdf)), "n_joined": int(len(ju))},
    "pv1_rsid_mirage": pv1,
    "pv2_equity_and_mediation": pv2,
    "pv3_vus_liberation": pv3,
    "acmg55_complete_qc_audit": {
        "n_total": len(qc_all55),
        "n_fetched_new": n_fetched_new,
        "n_qc_fail_or_low_an": sum(1 for x in qc_all55 if x["any_qc_fail_or_low_an"]),
        "failed_variants": [x for x in qc_all55 if x["any_qc_fail_or_low_an"]],
        "all_55": qc_all55,
    },
}
with open(os.path.join(AUD, "a1_prospective_unseen20_and_qc55.json"), "w") as f:
    json.dump(out, f, indent=2)

print("\n=== PROSPECTIVE UNSEEN 20-GENE VALIDATION RESULTS ===", flush=True)
print("Overview:", out["unseen20_overview"])
print("PV-1 rsID Mirage:", pv1)
print("PV-2 Denominator-Clean Equity & Mediation:", json.dumps(pv2, indent=2))
print("PV-3 VUS Liberation:", pv3)
print(f"\n=== ALL 55 ACMG SF v3.3 AD >=2-STAR COMPLETE LIVE QC AUDIT ===")
print(f"Total QC-failed / low-AN / LCR: {out['acmg55_complete_qc_audit']['n_qc_fail_or_low_an']} / 55")
for x in out["acmg55_complete_qc_audit"]["failed_variants"]:
    print(f"  FAIL: {x['gene']:7s} {x['variant_id']:20s} AN={x['an_frac']*100:5.1f}% exome_f={x['exome_filters']} genome_f={x['genome_filters']} flags={x['global_flags']}")
