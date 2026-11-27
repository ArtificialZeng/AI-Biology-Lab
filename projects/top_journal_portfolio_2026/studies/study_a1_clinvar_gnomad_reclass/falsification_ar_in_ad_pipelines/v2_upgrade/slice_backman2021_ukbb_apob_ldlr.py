#!/opt/anaconda3/bin/python3
"""
1. Audit the 7 non-RYR1/DES/CACNA1S lof=no pLoF alleles in unbiased_7427_ad_gold_labels.csv
   (ACTA2, MYH7, PLN, TNNI3, TNNT2) to confirm none is an AD-valid pLoF.
2. Slice the official EBI GWAS Catalog tabix-indexed harmonized summary statistics of
   Backman et al. (Nature 2021, PMID 34662886; UK Biobank 454,787 exomes) across the
   APOB (chr2:21001000-21045000) and LDLR (chr19:11089000-11134000) loci for:
     - GCST90079033: Direct LDL-C (UKB field 30780, N=411,351)
     - GCST90079024: Total Cholesterol (UKB field 30690, N=412,075)
     - GCST90079018: Apolipoprotein B (UKB field 30640, N=410,113)
     - GCST90081493: Lipid-lowering medication use (N=320,758)
     - GCST90077883: High cholesterol (N=331,754)
   Yielding SAME-STUDY, SAME-SCALE UK Biobank measured lipid & statin effects across ALL
   ClinVar 2-4star APOB pLoFs vs APOB/LDLR true FH variants!
"""
import gzip, json, os, re, struct, urllib.request, zlib
import pandas as pd

v2 = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/v2_upgrade"
df = pd.read_csv(f"{v2}/unbiased_7427_ad_gold_labels.csv")

print("=== Audit OTHER_LOFNO_PLOF rows (2-4star >=80% AN) ===")
other = df[(df["stars"] >= 2) & (df["hi_call"]) & (df["gold_cat"] == "NOT_MEETING_OTHER_LOFNO_PLOF")]
print(other[["variant_id", "variation_id", "gene", "stars", "clnsig", "hgvs_name", "het_joint"]].to_string(index=False))

opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
UA = {"User-Agent": "A1Xv2-BackmanSlicer/1.0"}

def fetch_bytes(url, start=None, end=None, timeout=45):
    hdrs = dict(UA)
    if start is not None and end is not None:
        hdrs["Range"] = f"bytes={start}-{end}"
    req = urllib.request.Request(url, headers=hdrs)
    with opener.open(req, timeout=timeout) as r:
        return r.read()

def reg2bins(beg, end):
    end -= 1
    bins = [0]
    for k in range(1 + (beg >> 26), 2 + (end >> 26)): bins.append(k)
    for k in range(9 + (beg >> 23), 10 + (end >> 23)): bins.append(k)
    for k in range(73 + (beg >> 20), 74 + (end >> 20)): bins.append(k)
    for k in range(585 + (beg >> 17), 586 + (end >> 17)): bins.append(k)
    for k in range(4681 + (beg >> 14), 4682 + (end >> 14)): bins.append(k)
    return set(bins)

def parse_tbi(tbi_raw_gz):
    data = gzip.decompress(tbi_raw_gz)
    assert data[0:4] == b"TBI\x01", f"Bad TBI magic: {data[0:4]}"
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from("<8i", data, 4)
    off = 36
    names = [n.decode() for n in data[off:off+l_nm].split(b"\x00") if n]
    off += l_nm
    refs = {}
    for i in range(n_ref):
        n_bin, = struct.unpack_from("<i", data, off); off += 4
        bin_map = {}
        for _ in range(n_bin):
            b, n_chunk = struct.unpack_from("<Ii", data, off); off += 8
            chunks = []
            for _ in range(n_chunk):
                c_beg, c_end = struct.unpack_from("<QQ", data, off); off += 16
                chunks.append((c_beg, c_end))
            bin_map[b] = chunks
        n_intv, = struct.unpack_from("<i", data, off); off += 4
        intvs = list(struct.unpack_from(f"<{n_intv}Q", data, off)) if n_intv > 0 else []
        off += 8 * n_intv
        refs[names[i]] = (bin_map, intvs)
    return refs, col_seq, col_beg

def decompress_bgzf(raw_bytes):
    out = []
    pos = 0
    n = len(raw_bytes)
    while pos + 18 <= n:
        if raw_bytes[pos] != 31 or raw_bytes[pos+1] != 139:
            break
        xlen, = struct.unpack_from("<H", raw_bytes, pos+10)
        bsize = None
        epos = pos + 12
        eend = epos + xlen
        while epos + 4 <= eend:
            si1, si2, slen = struct.unpack_from("<BBH", raw_bytes, epos)
            if si1 == 66 and si2 == 67 and slen == 2:
                bsize, = struct.unpack_from("<H", raw_bytes, epos+4)
                break
            epos += 4 + slen
        if bsize is None:
            break
        block_len = bsize + 1
        if pos + block_len > n:
            break
        block = raw_bytes[pos:pos+block_len]
        out.append(zlib.decompress(block[12 + xlen : block_len - 8], -15))
        pos += block_len
    return b"".join(out)

def discover_harmonized_url(gcst):
    num = int(gcst.replace("GCST", ""))
    lo = ((num - 1) // 1000) * 1000 + 1
    hi = lo + 999
    bucket = f"GCST{lo:06d}-GCST{hi:06d}"
    hdir = f"http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{bucket}/{gcst}/harmonised/"
    html = fetch_bytes(hdir).decode("utf-8", errors="replace")
    tbi_files = re.findall(r'href="([^"]+\.h\.tsv\.gz\.tbi)"', html)
    if tbi_files:
        tsv_name = tbi_files[0][:-4]
        return hdir + tsv_name, hdir + tbi_files[0]
    rdir = f"http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{bucket}/{gcst}/"
    html2 = fetch_bytes(rdir).decode("utf-8", errors="replace")
    tbi2 = re.findall(r'href="([^"]+\.tsv\.gz\.tbi)"', html2)
    if tbi2:
        return rdir + tbi2[0][:-4], rdir + tbi2[0]
    return None, None

def slice_gwas_region(tsv_url, tbi_url, chrom_num, beg1, end1):
    # Fetch header from first 64KB
    hdr_raw = fetch_bytes(tsv_url, 0, 65535)
    hdr_lines = decompress_bgzf(hdr_raw).decode("utf-8", errors="replace").splitlines()
    header = hdr_lines[0].split("\t")
    tbi_refs, col_seq, col_beg = parse_tbi(fetch_bytes(tbi_url))
    cname = str(chrom_num) if str(chrom_num) in tbi_refs else f"chr{chrom_num}"
    if cname not in tbi_refs:
        return header, []
    bin_map, intvs = tbi_refs[cname]
    beg0, end0 = max(0, beg1 - 1), end1
    min_off = intvs[beg0 >> 14] if (beg0 >> 14) < len(intvs) else 0
    ranges = []
    # Use leaf bins (4681..) + immediate 128kb parent bins (585..) to avoid bin-0 blowup
    for b in range(4681 + (beg0 >> 14), 4682 + (end0 >> 14)):
        for c_beg, c_end in bin_map.get(b, []):
            if c_end > min_off:
                ranges.append((c_beg >> 16, (c_end >> 16) + 65536))
    for b in range(585 + (beg0 >> 17), 586 + (end0 >> 17)):
        for c_beg, c_end in bin_map.get(b, []):
            if c_end > min_off and ((c_end >> 16) - (c_beg >> 16)) < 2_000_000:
                ranges.append((c_beg >> 16, (c_end >> 16) + 65536))
    if not ranges:
        return header, []
    ranges.sort()
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1] + 65536:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    rows = []
    seen = set()
    for s, e in merged:
        dec = decompress_bgzf(fetch_bytes(tsv_url, s, e)).decode("utf-8", errors="replace")
        for line in dec.splitlines():
            if not line or line.startswith("hm_") or line.startswith("chromosome"):
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                continue
            ch_val = parts[col_seq - 1].replace("chr", "")
            if ch_val != str(chrom_num):
                continue
            try:
                pos_val = int(parts[col_beg - 1])
            except ValueError:
                continue
            if beg1 <= pos_val <= end1:
                k = (pos_val, line[:80])
                if k not in seen:
                    seen.add(k)
                    rows.append(dict(zip(header, parts)))
    return header, rows

# Build lookup of ALL ClinVar P/LP alleles in APOB and LDLR (from unbiased_7427_ad_gold_labels.csv)
target_df = df[df["gene"].isin(["APOB", "LDLR"])].copy()
target_lkp = {}
for _, r in target_df.iterrows():
    ch, pos, ref, alt = str(r["variant_id"]).split("-")
    target_lkp[(ch, int(pos), ref.upper(), alt.upper())] = {
        "variant_id": r["variant_id"],
        "gene": r["gene"],
        "stars": int(r["stars"]),
        "hgvs": r["hgvs_name"],
        "gold_cat": r["gold_cat"],
        "is_plof": bool(r["is_plof_mc"]),
        "het_joint": int(r["het_joint"]),
    }

studies = [
    ("GCST90079033", "Direct_LDL_UKB30780", 411351, "SD"),
    ("GCST90079024", "Total_Cholesterol_UKB30690", 412075, "SD"),
    ("GCST90079018", "ApoB_UKB30640", 410113, "SD"),
    ("GCST90081493", "Lipid_Lowering_Meds", 320758, "logOR"),
    ("GCST90077883", "High_Cholesterol_Binary", 331754, "logOR"),
]
loci = [
    ("APOB", 2, 21001000, 21045000),
    ("LDLR", 19, 11089000, 11134000),
]

matched_rows = []
for gcst, trait_short, n_samp, scale in studies:
    tsv_u, tbi_u = discover_harmonized_url(gcst)
    print(f"\n[{gcst}] {trait_short} -> tsv={tsv_u}")
    if not tsv_u:
        continue
    for gene, ch, b1, e1 in loci:
        hdr, rows = slice_gwas_region(tsv_u, tbi_u, ch, b1, e1)
        m_cnt = 0
        for row in rows:
            # Harmonized EBI columns: hm_chrom/chromosome, hm_pos/base_pair_location, effect_allele, other_allele, beta, standard_error, p_value
            c_str = (row.get("hm_chrom") or row.get("chromosome") or "").replace("chr", "")
            p_str = row.get("hm_pos") or row.get("base_pair_location") or ""
            ea = (row.get("hm_effect_allele") or row.get("effect_allele") or "").upper()
            oa = (row.get("hm_other_allele") or row.get("other_allele") or "").upper()
            beta_s = row.get("hm_beta") or row.get("beta")
            se_s = row.get("standard_error")
            pv_s = row.get("p_value")
            ef_s = row.get("hm_effect_allele_frequency") or row.get("effect_allele_frequency")
            if not p_str or not beta_s or beta_s == "NA":
                continue
            pos = int(float(p_str))
            # Check orientation: ClinVar REF=oa, ALT=ea -> sign=+1; ClinVar REF=ea, ALT=oa -> sign=-1
            meta = target_lkp.get((c_str, pos, oa, ea))
            sign = 1.0
            if meta is None:
                meta = target_lkp.get((c_str, pos, ea, oa))
                sign = -1.0
            if meta is None:
                continue
            try:
                beta_val = float(beta_s) * sign
                se_val = float(se_s) if se_s and se_s != "NA" else None
                pv_val = float(pv_s) if pv_s and pv_s != "NA" else None
                eaf_val = float(ef_s) if ef_s and ef_s != "NA" else None
            except ValueError:
                continue
            m_cnt += 1
            matched_rows.append({
                "study_id": gcst, "trait": trait_short, "n_samples": n_samp, "scale": scale,
                "variant_id": meta["variant_id"], "gene": meta["gene"], "stars": meta["stars"],
                "gold_cat": meta["gold_cat"], "is_plof": meta["is_plof"], "gnomad_het": meta["het_joint"],
                "hgvs": meta["hgvs"], "alt_beta": beta_val, "se": se_val, "p_value": pv_val,
                "eaf_ukb": (eaf_val if sign > 0 else (1.0 - eaf_val if eaf_val is not None else None)),
            })
        print(f"  {gene} (chr{ch}:{b1}-{e1}): scanned {len(rows)} GWAS rows -> matched {m_cnt} ClinVar P/LP alleles!")

mdf = pd.DataFrame(matched_rows)
mdf.to_csv(f"{v2}/backman2021_ukbb_exact_apob_ldlr_clinvar_matches.csv", index=False)
print(f"\nTOTAL matched Backman 2021 UKBB rows: {len(mdf)} across {mdf['variant_id'].nunique() if len(mdf) else 0} unique ClinVar P/LP alleles!")
if len(mdf):
    for tr in ["Direct_LDL_UKB30780", "Total_Cholesterol_UKB30690", "ApoB_UKB30640", "Lipid_Lowering_Meds"]:
        sub = mdf[(mdf["trait"] == tr) & (mdf["stars"] >= 2)].sort_values("alt_beta")
        print(f"\n--- {tr} (2-4star ClinVar P/LP, n={len(sub)}) ---")
        print(sub[["gene", "gold_cat", "variant_id", "hgvs", "gnomad_het", "alt_beta", "se", "p_value"]].to_string(index=False))
