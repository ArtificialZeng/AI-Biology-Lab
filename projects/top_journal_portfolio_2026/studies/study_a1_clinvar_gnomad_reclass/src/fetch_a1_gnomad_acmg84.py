#!/opt/anaconda3/bin/python3
"""
Step 1 for Study A-1:
  1. Parse NCBI ClinVar GRCh38 clinvar.vcf.gz + variant_summary.txt.gz for ALL 84 ACMG SF v3.3 genes
     (preserving exact GRCh38 chrom-pos-ref-alt, VariationID, CLNSIG, CLNREVSTAT, stars, MC, CLNDN, RS, LastEvaluated, NumberSubmitters).
  2. Fetch ALL gnomAD v4.1 (gnomad_r4) variants across ALL 84 ACMG SF v3.3 genes using the official
     gnomad-database skill wrapper (search_variants.py) with automatic resume & region-chunking fallback
     for ultra-large genes (e.g., TTN).
"""
import gzip
import json
import os
import subprocess
import sys
import time
import hashlib
import pandas as pd

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
RAW = os.path.join(BASE, "data", "raw")
GNOMAD_DIR = os.path.join(RAW, "gnomad_v4_acmg84")
CLEAN = os.path.join(BASE, "data", "clean")
os.makedirs(GNOMAD_DIR, exist_ok=True)
os.makedirs(CLEAN, exist_ok=True)

GNOMAD_SCRIPT = "/Users/mac/.gemini/config/plugins/genomics/skills/gnomad_database/scripts/search_variants.py"

# Complete authoritative 84-gene ACMG SF v3.3 list (Lee et al., Genetics in Medicine 2025; 27(8):101454)
# Tuple: (gene, domain, inheritance, version_added)
ACMG_SF_V33 = [
    # Cancer (28 genes: 25 v2.0 + 3 v3.0)
    ("BRCA1", "Cancer", "AD", "v2.0"), ("BRCA2", "Cancer", "AD", "v2.0"),
    ("TP53", "Cancer", "AD", "v2.0"), ("STK11", "Cancer", "AD", "v2.0"),
    ("MLH1", "Cancer", "AD", "v2.0"), ("MSH2", "Cancer", "AD", "v2.0"),
    ("MSH6", "Cancer", "AD", "v2.0"), ("PMS2", "Cancer", "AD", "v2.0"),
    ("APC", "Cancer", "AD", "v2.0"), ("MUTYH", "Cancer", "AR", "v2.0"),
    ("VHL", "Cancer", "AD", "v2.0"), ("MEN1", "Cancer", "AD", "v2.0"),
    ("RET", "Cancer", "AD", "v2.0"), ("PTEN", "Cancer", "AD", "v2.0"),
    ("RB1", "Cancer", "AD", "v2.0"), ("SDHD", "Cancer", "AD", "v2.0"),
    ("SDHAF2", "Cancer", "AD", "v2.0"), ("SDHC", "Cancer", "AD", "v2.0"),
    ("SDHB", "Cancer", "AD", "v2.0"), ("TSC1", "Cancer", "AD", "v2.0"),
    ("TSC2", "Cancer", "AD", "v2.0"), ("WT1", "Cancer", "AD", "v2.0"),
    ("NF2", "Cancer", "AD", "v2.0"), ("BMPR1A", "Cancer", "AD", "v2.0"),
    ("SMAD4", "Cancer", "AD", "v2.0"),
    ("PALB2", "Cancer", "AD", "v3.0"), ("MAX", "Cancer", "AD", "v3.0"),
    ("TMEM127", "Cancer", "AD", "v3.0"),
    # Cardiovascular (43 genes: 30 v2.0 + 4 v3.0 + 5 v3.1 + 3 v3.2 + 1 v3.3)
    ("COL3A1", "Cardiovascular", "AD", "v2.0"), ("FBN1", "Cardiovascular", "AD", "v2.0"),
    ("TGFBR1", "Cardiovascular", "AD", "v2.0"), ("TGFBR2", "Cardiovascular", "AD", "v2.0"),
    ("SMAD3", "Cardiovascular", "AD", "v2.0"), ("ACTA2", "Cardiovascular", "AD", "v2.0"),
    ("MYH11", "Cardiovascular", "AD", "v2.0"), ("MYBPC3", "Cardiovascular", "AD", "v2.0"),
    ("MYH7", "Cardiovascular", "AD", "v2.0"), ("TNNT2", "Cardiovascular", "AD", "v2.0"),
    ("TNNI3", "Cardiovascular", "AD", "v2.0"), ("TPM1", "Cardiovascular", "AD", "v2.0"),
    ("MYL3", "Cardiovascular", "AD", "v2.0"), ("ACTC1", "Cardiovascular", "AD", "v2.0"),
    ("PRKAG2", "Cardiovascular", "AD", "v2.0"), ("GLA", "Cardiovascular", "XL", "v2.0"),
    ("MYL2", "Cardiovascular", "AD", "v2.0"), ("LMNA", "Cardiovascular", "AD", "v2.0"),
    ("RYR2", "Cardiovascular", "AD", "v2.0"), ("PKP2", "Cardiovascular", "AD", "v2.0"),
    ("DSP", "Cardiovascular", "AD", "v2.0"), ("DSC2", "Cardiovascular", "AD", "v2.0"),
    ("TMEM43", "Cardiovascular", "AD", "v2.0"), ("DSG2", "Cardiovascular", "AD", "v2.0"),
    ("KCNQ1", "Cardiovascular", "AD", "v2.0"), ("KCNH2", "Cardiovascular", "AD", "v2.0"),
    ("SCN5A", "Cardiovascular", "AD", "v2.0"), ("LDLR", "Cardiovascular", "AD", "v2.0"),
    ("APOB", "Cardiovascular", "AD", "v2.0"), ("PCSK9", "Cardiovascular", "AD", "v2.0"),
    ("TTN", "Cardiovascular", "AD", "v3.0"), ("FLNC", "Cardiovascular", "AD", "v3.0"),
    ("TRDN", "Cardiovascular", "AR", "v3.0"), ("CASQ2", "Cardiovascular", "AR", "v3.0"),
    ("BAG3", "Cardiovascular", "AD", "v3.1"), ("DES", "Cardiovascular", "AD", "v3.1"),
    ("RBM20", "Cardiovascular", "AD", "v3.1"), ("TNNC1", "Cardiovascular", "AD", "v3.1"),
    ("TTR", "Cardiovascular", "AD", "v3.1"),
    ("CALM1", "Cardiovascular", "AD", "v3.2"), ("CALM2", "Cardiovascular", "AD", "v3.2"),
    ("CALM3", "Cardiovascular", "AD", "v3.2"),
    ("PLN", "Cardiovascular", "AD", "v3.3"),
    # Inborn Errors of Metabolism (5 genes: 1 v2.0 + 2 v3.0 + 2 v3.3)
    ("OTC", "Metabolic", "XL", "v2.0"),
    ("GAA", "Metabolic", "AR", "v3.0"), ("BTD", "Metabolic", "AR", "v3.0"),
    ("ABCD1", "Metabolic", "XL", "v3.3"), ("CYP27A1", "Metabolic", "AR", "v3.3"),
    # Miscellaneous (8 genes: 3 v2.0 + 5 v3.0)
    ("RYR1", "Miscellaneous", "AD", "v2.0"), ("CACNA1S", "Miscellaneous", "AD", "v2.0"),
    ("ATP7B", "Miscellaneous", "AR", "v2.0"),
    ("HFE", "Miscellaneous", "AR", "v3.0"), ("HNF1A", "Miscellaneous", "AD", "v3.0"),
    ("RPE65", "Miscellaneous", "AR", "v3.0"), ("ACVRL1", "Miscellaneous", "AD", "v3.0"),
    ("ENG", "Miscellaneous", "AD", "v3.0"),
]

assert len(ACMG_SF_V33) == 84, f"Expected 84 genes, got {len(ACMG_SF_V33)}"
assert len({g[0] for g in ACMG_SF_V33}) == 84, "Duplicate gene in ACMG_SF_V33!"
GENE_SET = {g[0] for g in ACMG_SF_V33}
GENE_META = {g[0]: {"domain": g[1], "inheritance": g[2], "version": g[3]} for g in ACMG_SF_V33}


def review_status_to_stars(rs):
    s = (rs or "").lower().replace("_", " ")
    if "practice guideline" in s:
        return 4
    if "reviewed by expert panel" in s:
        return 3
    if "criteria provided, multiple submitters, no conflicts" in s:
        return 2
    if "criteria provided, single submitter" in s or "criteria provided, conflicting" in s:
        return 1
    return 0


def parse_clinvar_acmg84():
    """Extract all GRCh38 ClinVar records for the 84 ACMG SF v3.3 genes from clinvar.vcf.gz + variant_summary.txt.gz."""
    out_pq = os.path.join(CLEAN, "clinvar_acmg84_grch38.parquet")
    if os.path.exists(out_pq) and os.path.getsize(out_pq) > 100000:
        print("ClinVar ACMG84 parquet already exists, skipping parse.", flush=True)
        return pd.read_parquet(out_pq)

    t0 = time.time()
    vcf_path = os.path.join(RAW, "clinvar.vcf.gz")
    vs_path = os.path.join(RAW, "variant_summary.txt.gz")

    print("Parsing clinvar.vcf.gz for all 84 ACMG SF v3.3 genes...", flush=True)
    rows = []
    total_vcf = 0
    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            total_vcf += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, var_id, ref, alt, _, _, info = parts[:8]
            # Fast gene pre-filter
            gi_idx = info.find("GENEINFO=")
            if gi_idx < 0:
                continue
            gi_end = info.find(";", gi_idx)
            gi_val = info[gi_idx + 9: gi_end if gi_end >= 0 else None]
            genes_in_rec = [tok.split(":")[0] for tok in gi_val.split("|")]
            matched = [g for g in genes_in_rec if g in GENE_SET]
            if not matched:
                continue
            # Parse INFO dict
            kv = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    kv[k] = v
            gene = matched[0]
            vid = f"{chrom}-{pos}-{ref}-{alt}"
            rev = kv.get("CLNREVSTAT", "")
            rows.append({
                "variant_id": vid,
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "variation_id": str(var_id),
                "gene": gene,
                "domain": GENE_META[gene]["domain"],
                "inheritance": GENE_META[gene]["inheritance"],
                "acmg_version": GENE_META[gene]["version"],
                "clnsig": kv.get("CLNSIG", ""),
                "clnsigconf": kv.get("CLNSIGCONF", ""),
                "clnrevstat": rev,
                "stars": review_status_to_stars(rev),
                "mc": kv.get("MC", ""),
                "clndn": kv.get("CLNDN", ""),
                "rs": ("rs" + kv["RS"]) if kv.get("RS") else "",
                "clnvc": kv.get("CLNVC", ""),
            })

    df_vcf = pd.DataFrame(rows)
    print(f"Parsed {total_vcf:,} total ClinVar VCF rows -> {len(df_vcf):,} rows in 84 ACMG genes ({df_vcf['gene'].nunique()} genes) in {time.time()-t0:.1f}s", flush=True)

    # Enrich with LastEvaluated, NumberSubmitters, OriginSimple from variant_summary.txt.gz (GRCh38 rows)
    print("Enriching with LastEvaluated & NumberSubmitters from variant_summary.txt.gz...", flush=True)
    t1 = time.time()
    vid_set = set(df_vcf["variation_id"])
    vs_map = {}
    with gzip.open(vs_path, "rt", encoding="utf-8", errors="ignore") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_asm = header.index("Assembly")
        idx_vid = header.index("VariationID")
        idx_le = header.index("LastEvaluated")
        idx_ns = header.index("NumberSubmitters")
        idx_orig = header.index("OriginSimple")
        idx_name = header.index("Name")
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(idx_asm, idx_vid, idx_le, idx_ns, idx_orig, idx_name):
                continue
            if cols[idx_asm] != "GRCh38":
                continue
            v = cols[idx_vid]
            if v in vid_set and v not in vs_map:
                vs_map[v] = (cols[idx_le], cols[idx_ns], cols[idx_orig], cols[idx_name])

    df_vcf["last_evaluated"] = df_vcf["variation_id"].map(lambda x: vs_map.get(x, ("", "", "", ""))[0])
    df_vcf["n_submitters"] = pd.to_numeric(df_vcf["variation_id"].map(lambda x: vs_map.get(x, ("", "", "", ""))[1]), errors="coerce")
    df_vcf["origin_simple"] = df_vcf["variation_id"].map(lambda x: vs_map.get(x, ("", "", "", ""))[2])
    df_vcf["hgvs_name"] = df_vcf["variation_id"].map(lambda x: vs_map.get(x, ("", "", "", ""))[3])
    print(f"Enriched {len(vs_map):,}/{len(vid_set):,} VariationIDs in {time.time()-t1:.1f}s", flush=True)

    df_vcf.to_parquet(out_pq, index=False)
    return df_vcf


def fetch_gnomad_gene(gene, chrom, g_start, g_end):
    """Fetch gnomAD v4.1 variants for a single gene using official search_variants.py wrapper, with chunked-region fallback."""
    out_file = os.path.join(GNOMAD_DIR, f"{gene}_gnomad_r4.json")
    if os.path.exists(out_file) and os.path.getsize(out_file) > 500:
        try:
            d = json.load(open(out_file))
            vs = (d.get("data", {}).get("gene") or {}).get("variants")
            if vs is not None and len(vs) > 0:
                return len(vs), "cached"
            vs2 = d.get("variants_chunked")
            if vs2 is not None and len(vs2) > 0:
                return len(vs2), "cached_chunked"
        except Exception:
            pass

    # Try standard --gene call first (official wrapper handles 10 QPM rate limit)
    for attempt in range(2):
        cmd = ["uv", "run", GNOMAD_SCRIPT, "--gene", gene, "--dataset", "gnomad_r4", "--output", out_file]
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if p.returncode == 0 and os.path.exists(out_file) and os.path.getsize(out_file) > 200:
            try:
                d = json.load(open(out_file))
                vs = (d.get("data", {}).get("gene") or {}).get("variants")
                if vs is not None and len(vs) > 0:
                    return len(vs), "gene_ok"
            except Exception:
                pass
        time.sleep(6.5)

    # Fallback: region chunking (for ultra-large genes like TTN, RYR1, RYR2, FBN1, DMD)
    print(f"  [{gene}] falling back to chunked region query ({chrom}:{g_start}-{g_end})...", flush=True)
    chunk_size = 50000
    all_vars = {}
    cur = max(1, int(g_start) - 1000)
    end = int(g_end) + 1000
    c_idx = 0
    while cur <= end:
        c_idx += 1
        c_end = min(end, cur + chunk_size - 1)
        tmp_f = os.path.join(GNOMAD_DIR, f"{gene}_chunk_{c_idx:03d}.json")
        cmd = ["uv", "run", GNOMAD_SCRIPT, "--chrom", str(chrom), "--start", str(cur), "--end", str(c_end), "--dataset", "gnomad_r4", "--output", tmp_f]
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if p.returncode == 0 and os.path.exists(tmp_f):
            try:
                dc = json.load(open(tmp_f))
                for v in (dc.get("data", {}).get("region") or {}).get("variants") or []:
                    all_vars[v["variant_id"]] = v
            except Exception:
                pass
        time.sleep(6.2)
        cur = c_end + 1

    merged = {"gene": gene, "variants_chunked": list(all_vars.values())}
    with open(out_file, "w") as f:
        json.dump(merged, f)
    return len(all_vars), "chunked_ok"


def main():
    t0 = time.time()
    df_cv = parse_clinvar_acmg84()
    # Compute gene genomic spans from ClinVar coordinates for chunking fallback
    spans = df_cv.groupby("gene").agg(chrom=("chrom", "first"), g_start=("pos", "min"), g_end=("pos", "max")).to_dict(orient="index")

    # Fetch gnomAD v4.1 for all 84 ACMG SF v3.3 genes (polite 6.1s spacing -> <=10 QPM)
    summary = []
    for i, (gene, dom, inh, ver) in enumerate(ACMG_SF_V33, 1):
        sp = spans.get(gene, {"chrom": "", "g_start": 0, "g_end": 0})
        t_g = time.time()
        n_v, mode = fetch_gnomad_gene(gene, sp["chrom"], sp["g_start"], sp["g_end"])
        summary.append({"gene": gene, "domain": dom, "inheritance": inh, "version": ver, "gnomad_n": n_v, "mode": mode})
        if mode != "cached" and mode != "cached_chunked":
            elapsed_g = time.time() - t_g
            if elapsed_g < 6.1:
                time.sleep(6.1 - elapsed_g)
        if i % 10 == 0 or i == 84 or n_v == 0:
            print(f"[{i:02d}/84] {gene:8s} ({dom[:6]}/{inh}): {n_v:6d} gnomAD v4.1 variants ({mode})", flush=True)

    man = {
        "completed_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_sec": round(time.time() - t0, 1),
        "n_genes_total": len(ACMG_SF_V33),
        "n_genes_nonzero": sum(1 for x in summary if x["gnomad_n"] > 0),
        "total_gnomad_variants": sum(x["gnomad_n"] for x in summary),
        "genes": summary,
    }
    with open(os.path.join(CLEAN, "gnomad_acmg84_manifest.json"), "w") as f:
        json.dump(man, f, indent=2)
    print(f"ALL 84 GENES FETCHED: {man['n_genes_nonzero']}/84 nonzero, {man['total_gnomad_variants']:,} variants in {man['elapsed_sec']}s", flush=True)


if __name__ == "__main__":
    main()
