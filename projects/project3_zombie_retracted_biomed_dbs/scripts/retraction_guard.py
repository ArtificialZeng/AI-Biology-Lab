#!/opt/anaconda3/bin/python3
"""
RetractionGuard: Lightweight, Zero-Cost Automated Retraction Synchronization & Filtering Gate
for Biomedical Databases, AI Knowledge Graphs (PubTator3 / Open Targets), and JATS XML LLM Corpora.

Usage (CLI or Python import):
    from retraction_guard import RetractionGuard
    rg = RetractionGuard.from_master_registry("data/master_retracted_pmids_100pct.json",
                                              "data/retraction_watch_clean_retractions.parquet")
    rg.is_retracted_pmid(32076465)       # -> True (plus mechanism metadata)
    rg.check_jats_xml(xml_bytes)         # -> Detects BOTH JATS <related-article> AND PMID/DOI registry match!
    rg.filter_tsv_stream(in_gz, out_gz)  # -> Filters PubTator3 / BioC / ChEMBL / ClinVar tables in streaming O(1) RAM!
"""
import argparse, gzip, json, os, re, time, xml.etree.ElementTree as ET
from typing import Dict, Optional, Set, Tuple


class RetractionGuard:
    """Fast O(1) hash-set retraction filter combining Crossref Retraction Watch + PubMed Retracted[pt],
    strictly excluding Reinstatements."""

    def __init__(self, retracted_pmids: Set[int], retracted_dois: Optional[Set[str]] = None,
                 pmid_reason: Optional[Dict[int, str]] = None):
        self.retracted_pmids = set(retracted_pmids)
        self.retracted_dois = {d.strip().lower() for d in (retracted_dois or set()) if d}
        self.pmid_reason = pmid_reason or {}

    @classmethod
    def from_master_registry(cls, master_json_path: str, rw_parquet_path: Optional[str] = None) -> "RetractionGuard":
        with open(master_json_path) as f:
            reg = json.load(f)
        pmids = set(int(x) for x in reg["union_pmids"])
        dois: Set[str] = set()
        reasons: Dict[int, str] = {}
        if rw_parquet_path and os.path.exists(rw_parquet_path):
            import pandas as pd
            rw = pd.read_parquet(rw_parquet_path, columns=["pmid_clean", "doi_clean", "reason_class"])
            for p, d, rc in zip(rw["pmid_clean"], rw["doi_clean"], rw["reason_class"]):
                if pd.notna(p) and int(p) in pmids:
                    reasons[int(p)] = str(rc)
                if isinstance(d, str) and d.startswith("10."):
                    dois.add(d.strip().lower())
        return cls(pmids, dois, reasons)

    def is_retracted_pmid(self, pmid: int) -> bool:
        return int(pmid) in self.retracted_pmids

    def is_retracted_doi(self, doi: str) -> bool:
        if not doi:
            return False
        d = re.sub(r"^https?://(dx\.)?doi\.org/", "", str(doi).strip().lower())
        return d in self.retracted_dois

    def check_jats_xml(self, xml_bytes: bytes) -> Tuple[bool, str]:
        """Inspect a PMC Open Access JATS XML document BEFORE extracting <abstract>/<body> text.
        Catches BOTH (1) JATS `<related-article related-article-type='retracted-article'>` metadata
        and (2) external Retraction Watch + PubMed PMID/DOI registry matches (including the 2.1%
        missing JATS tags or 1,865 RW-only articles unflagged in PubMed!)."""
        try:
            root = ET.fromstring(xml_bytes)
        except ET.ParseError:
            return False, "parse_error"

        # Gate 1: Explicit JATS XML article-type or related-article retraction metadata
        atype = (root.attrib.get("article-type") or "").lower()
        if "retract" in atype:
            return True, f"jats_article_type:{atype}"
        for ra in root.findall(".//related-article"):
            rat = (ra.attrib.get("related-article-type") or "").lower()
            if "retract" in rat or "retract" in "".join(ra.itertext()).lower():
                return True, f"jats_related_article:{rat}"

        # Gate 2: External Master Registry lookup by PMID / DOI
        for aid in root.findall(".//article-meta/article-id"):
            idtype = (aid.attrib.get("pub-id-type") or "").lower()
            val = (aid.text or "").strip()
            if idtype == "pmid" and val.isdigit() and self.is_retracted_pmid(int(val)):
                return True, f"registry_pmid:{self.pmid_reason.get(int(val), 'PUBMED_RETRACTED')}"
            if idtype == "doi" and self.is_retracted_doi(val):
                return True, "registry_doi"

        return False, "clean"

    def filter_tsv_stream(self, in_path: str, pmid_col: int = 0, out_path: Optional[str] = None,
                          has_header: bool = False) -> Dict[str, int]:
        """Streaming O(1)-memory retraction filter for PubTator3 / ClinVar / GOA / IntAct TSV/GZ dumps."""
        t0 = time.time()
        opener_in = gzip.open if in_path.endswith(".gz") else open
        kept, blocked, total = 0, 0, 0
        out_f = None
        if out_path:
            out_f = gzip.open(out_path, "wt", encoding="utf-8") if out_path.endswith(".gz") else open(out_path, "w", encoding="utf-8")
        try:
            with opener_in(in_path, "rt", encoding="utf-8", errors="replace") as f:
                if has_header:
                    hdr = f.readline()
                    if out_f:
                        out_f.write(hdr)
                for line in f:
                    total += 1
                    parts = line.split("\t")
                    is_ret = False
                    if len(parts) > pmid_col:
                        raw_id = parts[pmid_col].strip().replace("PMID:", "").replace("pubmed:", "")
                        if raw_id.isdigit() and int(raw_id) in self.retracted_pmids:
                            is_ret = True
                    if is_ret:
                        blocked += 1
                    else:
                        kept += 1
                        if out_f:
                            out_f.write(line)
        finally:
            if out_f:
                out_f.close()
        dt = time.time() - t0
        return {
            "total_rows": total,
            "blocked_retracted_rows": blocked,
            "kept_clean_rows": kept,
            "elapsed_sec": round(dt, 3),
            "throughput_rows_per_sec": int(total / max(1e-6, dt)),
        }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="RetractionGuard CLI benchmark & filter")
    ap.add_argument("--registry", required=True, help="Path to master_retracted_pmids_100pct.json")
    ap.add_argument("--rw-parquet", default=None, help="Optional path to retraction_watch_clean_retractions.parquet")
    ap.add_argument("--filter-tsv", default=None, help="Optional TSV/GZ file to benchmark streaming filter")
    ap.add_argument("--pmid-col", type=int, default=0, help="0-based PMID column index")
    ap.add_argument("--has-header", action="store_true")
    args = ap.parse_args()

    t_load = time.time()
    rg = RetractionGuard.from_master_registry(args.registry, args.rw_parquet)
    print(json.dumps({
        "loaded_pmids": len(rg.retracted_pmids),
        "loaded_dois": len(rg.retracted_dois),
        "load_sec": round(time.time() - t_load, 3),
    }, indent=2))
    if args.filter_tsv:
        res = rg.filter_tsv_stream(args.filter_tsv, pmid_col=args.pmid_col, has_header=args.has_header)
        print(json.dumps(res, indent=2))
