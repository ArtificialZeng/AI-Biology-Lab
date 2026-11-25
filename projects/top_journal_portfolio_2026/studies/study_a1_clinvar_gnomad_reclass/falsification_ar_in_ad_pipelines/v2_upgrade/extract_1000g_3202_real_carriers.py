#!/opt/anaconda3/bin/python3
"""
Pure-Python Tabix (.tbi) + BGZF HTTP Range Slicer for the 1000 Genomes Project
3,202 High-Coverage (30x) NYGC GRCh38 Phased WGS Cohort (2,504 unrelated + 698 related).
Downloads ONLY the tiny .tbi indexes and the exact BGZF byte ranges overlapping all 73 AD
ACMG SF v3.3 gene intervals (zero external C/bcftools dependencies!), then extracts EVERY
real individual-level heterozygous/homozygous carrier of ANY ClinVar P/LP allele across all
2,504 unrelated and 3,202 total real human genomes across AFR, AMR, EAS, EUR, SAS!
"""
import gzip, io, json, os, struct, urllib.request, zlib
from collections import defaultdict
import pandas as pd

BASE = "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass"
OUT_DIR = os.path.join(BASE, "falsification_ar_in_ad_pipelines/v2_upgrade")
os.makedirs(OUT_DIR, exist_ok=True)

opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
UA = {"User-Agent": "A1Xv2-1000G-Slicer/1.0"}

def fetch_bytes(url, start=None, end=None, timeout=60):
    hdrs = dict(UA)
    if start is not None and end is not None:
        hdrs["Range"] = f"bytes={start}-{end}"
    req = urllib.request.Request(url, headers=hdrs)
    with opener.open(req, timeout=timeout) as r:
        return r.read()

def reg2bins(beg, end):
    """Standard CSI/Tabix UCSC binning scheme (0-based half-open [beg, end))."""
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
    off = 0
    magic = data[0:4]
    assert magic == b"TBI\x01", f"Bad TBI magic: {magic}"
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from("<8i", data, 4)
    off = 36
    names = data[off:off+l_nm].split(b"\x00")
    names = [n.decode() for n in names if n]
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
    return refs

def decompress_bgzf_stream(raw_bytes):
    """Decompress concatenated BGZF blocks starting from a block boundary."""
    out = []
    pos = 0
    n = len(raw_bytes)
    while pos + 18 <= n:
        if raw_bytes[pos] != 31 or raw_bytes[pos+1] != 139:
            break
        flg = raw_bytes[pos+3]
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
        # Inflate raw deflate payload (skip header 12+xlen, strip 8-byte trailer)
        cdata = block[12 + xlen : block_len - 8]
        d = zlib.decompress(cdata, -15)
        out.append(d)
        pos += block_len
    return b"".join(out)

def fetch_vcf_header_samples(vcf_url):
    # Fetch first 128 KB (enough for VCF header + #CHROM line)
    raw = fetch_bytes(vcf_url, 0, 262143)
    txt = decompress_bgzf_stream(raw).decode("utf-8", errors="replace")
    for line in txt.splitlines():
        if line.startswith("#CHROM"):
            cols = line.split("\t")
            return cols[9:]
    raise RuntimeError(f"No #CHROM header in first 256KB of {vcf_url}")

def query_tabix_regions(vcf_url, tbi_refs, chrom_name, intervals):
    """Fetch exact BGZF byte ranges covering sorted list of (beg1, end1) 1-based intervals."""
    if chrom_name not in tbi_refs:
        return []
    bin_map, intvs = tbi_refs[chrom_name]
    byte_ranges = []
    for beg1, end1 in intervals:
        beg0 = max(0, beg1 - 1)
        end0 = end1
        min_off = intvs[beg0 >> 14] if (beg0 >> 14) < len(intvs) else 0
        bins = reg2bins(beg0, end0)
        for b in bins:
            for c_beg, c_end in bin_map.get(b, []):
                if c_end > min_off:
                    b_start = c_beg >> 16
                    b_end = (c_end >> 16) + 65536  # include full terminal BGZF block (max 64KB)
                    byte_ranges.append((b_start, b_end))
    if not byte_ranges:
        return []
    # Merge adjacent/nearby byte ranges (within 64KB) to minimize HTTP round-trips!
    byte_ranges.sort()
    merged = [list(byte_ranges[0])]
    for s, e in byte_ranges[1:]:
        if s <= merged[-1][1] + 65536:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    lines = []
    seen = set()
    for s, e in merged:
        raw = fetch_bytes(vcf_url, s, e)
        dec = decompress_bgzf_stream(raw).decode("utf-8", errors="replace")
        for line in dec.splitlines():
            if not line or line.startswith("#"):
                continue
            p = line.split("\t", 6)
            if len(p) < 6 or p[0] != chrom_name:
                continue
            try:
                pos = int(p[1])
            except ValueError:
                continue
            # Check if pos falls inside any requested interval
            if any(b1 <= pos <= e1 for b1, e1 in intervals):
                key = (p[0], pos, p[3], p[4])
                if key not in seen:
                    seen.add(key)
                    lines.append(line)
    return lines

# 1. Fetch official 1000G 3,202 pedigree / superpopulation table (identifies the exact 2,504 unrelated individuals + AFR/AMR/EAS/EUR/SAS!)
ped_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
ped_txt = fetch_bytes(ped_url).decode("utf-8")
ped_df = pd.read_csv(io.StringIO(ped_txt), sep=r"\s+")
ped_df.to_csv(os.path.join(OUT_DIR, "g1k_3202_ped_population.csv"), index=False)
# Unrelated = FatherID == '0' & MotherID == '0' (the standard 2,504 unrelated 1000G panel!)
unrel_set = set(ped_df.loc[(ped_df["FatherID"].astype(str) == "0") & (ped_df["MotherID"].astype(str) == "0"), "SampleID"])
pop_map = dict(zip(ped_df["SampleID"], ped_df["Superpopulation"]))
print(f"Loaded 1000G ped: total={len(ped_df)}, unrelated={len(unrel_set)}, superpops={ped_df['Superpopulation'].value_counts().to_dict()}")

# 2. Load ALL 73 AD ACMG SF v3.3 ClinVar P/LP alleles observed in gnomAD v4.1 (and ALL ClinVar P/LP across 73 AD genes!)
df_c = pd.read_parquet(os.path.join(BASE, "data/clean/clinvar_acmg84_grch38.parquet"))
def is_plp(s):
    if pd.isna(s): return False
    s = str(s)
    return ("Conflicting" not in s) and (("Pathogenic" in s) or ("Likely_pathogenic" in s))

ad_c = df_c[(df_c["inheritance"] == "AD") & (df_c["clnsig"].apply(is_plp))].copy()
# Group exact variant lookup by chr: (pos, ref, alt) -> row metadata!
# Use gnomAD-observed + all ClinVar P/LP positions on each autosome!
by_chr_lookup = defaultdict(dict)
by_chr_intervals = defaultdict(list)
for _, r in ad_c.iterrows():
    ch = str(r["chrom"]).replace("chr", "")
    if ch not in [str(i) for i in range(1, 23)]:
        continue
    pos = int(r["pos"])
    ref = str(r["ref"]).upper()
    alt = str(r["alt"]).upper()
    by_chr_lookup[f"chr{ch}"][(pos, ref, alt)] = {
        "variant_id": r["variant_id"],
        "variation_id": str(int(r["variation_id"])),
        "gene": r["gene"],
        "stars": int(r["stars"]),
        "clnsig": r["clnsig"],
        "hgvs": r["hgvs_name"],
        "mc": r["mc"],
    }
    by_chr_intervals[f"chr{ch}"].append((pos, pos))

# Merge gene-level intervals per chromosome (gap <= 100kb so we only fetch exact gene spans!)
merged_gene_iv = {}
for ch, ivs in by_chr_intervals.items():
    ivs.sort()
    m = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= m[-1][1] + 50000:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    merged_gene_iv[ch] = [(max(1, s - 10), e + 10) for s, e in m]

print(f"Chromosomes with AD ACMG SF v3.3 ClinVar P/LP alleles: {len(merged_gene_iv)} (total gene spans: {sum(len(v) for v in merged_gene_iv.values())})")

# 3. Slice ALL 22 autosomes of the 1000 Genomes 3,202 30x GRCh38 phased VCF!
base_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
samples_ref = None
carrier_records = []
total_matched_alleles = 0

for ch in sorted(merged_gene_iv.keys(), key=lambda x: int(x.replace("chr", ""))):
    vcf_url = f"{base_url}/1kGP_high_coverage_Illumina.{ch}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    tbi_url = vcf_url + ".tbi"
    tbi_raw = fetch_bytes(tbi_url)
    tbi_refs = parse_tbi(tbi_raw)
    if samples_ref is None:
        samples_ref = fetch_vcf_header_samples(vcf_url)
        print(f"1000G VCF sample count: {len(samples_ref)}")
    lines = query_tabix_regions(vcf_url, tbi_refs, ch, merged_gene_iv[ch])
    lkp = by_chr_lookup[ch]
    ch_matches = 0
    for line in lines:
        parts = line.split("\t")
        pos = int(parts[1])
        ref = parts[3].upper()
        alts = parts[4].upper().split(",")
        for alt_idx, alt in enumerate(alts, start=1):
            meta = lkp.get((pos, ref, alt))
            if meta is None:
                continue
            ch_matches += 1
            total_matched_alleles += 1
            alt_str = str(alt_idx)
            # Scan all 3,202 sample genotypes!
            for s_i, gt_field in enumerate(parts[9:]):
                gt = gt_field.split(":", 1)[0]
                alleles_gt = gt.replace("|", "/").split("/")
                n_alt = sum(1 for a in alleles_gt if a == alt_str)
                if n_alt > 0:
                    sid = samples_ref[s_i]
                    carrier_records.append({
                        "sample_id": sid,
                        "is_unrelated_2504": sid in unrel_set,
                        "superpop": pop_map.get(sid, "UNK"),
                        "zygosity": "hom" if n_alt >= 2 else "het",
                        "variant_id": meta["variant_id"],
                        "variation_id": meta["variation_id"],
                        "gene": meta["gene"],
                        "stars": meta["stars"],
                        "clnsig": meta["clnsig"],
                        "hgvs": meta["hgvs"],
                        "mc": meta["mc"],
                        "chrom": ch,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                    })
    print(f"  {ch:5s}: scanned {len(lines):6d} VCF rows -> matched {ch_matches:3d} ClinVar P/LP alleles (cum carrier records: {len(carrier_records)})")

cdf = pd.DataFrame(carrier_records)
cdf.to_csv(os.path.join(OUT_DIR, "g1k_3202_real_individual_ad_acmg_sf_carriers.csv"), index=False)
print(f"\nDONE! Extracted {len(cdf)} individual-level carrier calls across {total_matched_alleles} matched AD ACMG SF ClinVar P/LP alleles in 1000G 3,202 real WGS genomes!")
