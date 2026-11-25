#!/opt/anaconda3/bin/python3
"""
Analysis 2 (C2-ZOMBIE):
  Large-Scale Live AWS S3 `pmc-oa-opendata` Full-Text (.txt & .xml) Parser Audit!
  Tests a stratified sample of N = 240 retracted PMC Open Access articles across 2010-2024:
    - Downloads actual live `.txt` and `.xml` files directly from AWS S3 `pmc-oa-opendata` (the official open corpus powering biomedical LLMs!).
    - Inspects EXACTLY WHERE any retraction notice lives:
      (a) `.txt` raw file header (`==== Front`) vs `.txt` body (`==== Body`!)
      (b) `.xml` JATS standard fields extracted by LLM pipelines:
          - `<article-title>` (does the title say Retracted?)
          - `<abstract>` (does the abstract mention retraction?)
          - `<body>` (does the main scientific body text — methods/results/discussion — mention retraction, or is it 100% CLEAN unflagged scientific text?!)
          - `<article article-type="...">` attribute & `<related-article>` metadata tag!
    - Quantifies exact clean body characters & estimated LLM tokens (`chars / 4`) silently ingested by standard JATS XML `<abstract>+<body>` parsers (like S2ORC, The Pile, Dolma, BioMedLM, PMC-LLaMA)!
"""
import json, os, re, time, xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd
import requests

ROOT = "/Users/mac/AI-Biology-Lab/projects/project3_zombie_retracted_biomed_dbs"
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results")

S = requests.Session()
S.trust_env = False
S.headers.update({"User-Agent": "C2ZombieS3Audit/1.0 (mailto:research@example.org)"})

# Stratified sampling: 20 retracted PMC OA PMIDs per year bucket across 2012-2024 (12 buckets * 20 = 240 articles!)
print("=== [1] Drawing Stratified Sample of N=240 Retracted PMC OA Articles (2012-2024) ===", flush=True)
sampled_pmids = []
year_buckets = [
    ("2012_2014", "2012:2014[dp]", 30),
    ("2015_2016", "2015:2016[dp]", 30),
    ("2017_2018", "2017:2018[dp]", 30),
    ("2019", "2019[dp]", 25),
    ("2020", "2020[dp]", 25),
    ("2021", "2021[dp]", 25),
    ("2022", "2022[dp]", 25),
    ("2023", "2023[dp]", 25),
    ("2024", "2024[dp]", 25),
]
rng = np.random.default_rng(20261125)
for label, dp, k in year_buckets:
    term = f'Retracted Publication[pt] AND "pubmed pmc open access"[filter] AND {dp}'
    r = S.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
              params={"db": "pubmed", "term": term, "retmax": 1000, "retmode": "json"}, timeout=25).json()["esearchresult"]
    ids = r["idlist"]
    pick = rng.choice(ids, size=min(k, len(ids)), replace=False).tolist()
    for pid in pick:
        sampled_pmids.append((int(pid), label))
    time.sleep(0.35)

print(f"  Sampled {len(sampled_pmids)} stratified PMIDs; resolving PMCIDs via PubMed esummary...", flush=True)
pmc_items = []
for i in range(0, len(sampled_pmids), 100):
    chunk = sampled_pmids[i:i+100]
    ids_str = ",".join(str(p) for p, _ in chunk)
    r = S.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
              params={"db": "pubmed", "id": ids_str, "retmode": "json"}, timeout=25).json()["result"]
    for pid, bucket in chunk:
        it = r.get(str(pid), {})
        pmcid = next((x["value"] for x in it.get("articleids", []) if x.get("idtype") == "pmc"), None)
        if pmcid and pmcid.startswith("PMC"):
            pmc_items.append({
                "pmid": pid, "pmcid": pmcid, "year_bucket": bucket,
                "pubdate": it.get("pubdate", ""), "journal": it.get("source", ""),
                "pubmed_title": it.get("title", ""),
            })
    time.sleep(0.35)

print(f"  Resolved {len(pmc_items)} PMCIDs. Auditing live AWS S3 `pmc-oa-opendata` (.txt & .xml) in parallel (12 workers)...", flush=True)

def audit_one(item):
    sess = requests.Session()
    sess.trust_env = False
    pmcid = item["pmcid"]
    rec = dict(item)
    try:
        rl = sess.get(f"https://pmc-oa-opendata.s3.amazonaws.com/?list-type=2&prefix={pmcid}.", timeout=15)
        keys = re.findall(r"<Key>([^<]+)</Key>", rl.text)
        txt_key = next((k for k in keys if k.endswith(".txt")), None)
        xml_key = next((k for k in keys if k.endswith(".xml")), None)
        rec["s3_found"] = bool(keys)
        rec["has_txt"] = txt_key is not None
        rec["has_xml"] = xml_key is not None
        rec["txt_key"] = txt_key or ""
        rec["xml_key"] = xml_key or ""

        # Inspect TXT: split into ==== Front, ==== Body, ==== Refs
        if txt_key:
            rt = sess.get(f"https://pmc-oa-opendata.s3.amazonaws.com/{txt_key}", timeout=20)
            if rt.status_code == 200:
                txt = rt.text
                rec["txt_total_chars"] = len(txt)
                parts = re.split(r"====\s*(Front|Body|Refs)\s*", txt)
                sec_map = {}
                for j in range(1, len(parts)-1, 2):
                    sec_map[parts[j]] = parts[j+1]
                front_t = sec_map.get("Front", txt[:2000])
                body_t = sec_map.get("Body", "")
                rec["txt_front_chars"] = len(front_t)
                rec["txt_body_chars"] = len(body_t)
                rec["txt_front_has_retract"] = bool(re.search(r"retract", front_t, re.I))
                rec["txt_body_has_retract"] = bool(re.search(r"retract", body_t, re.I))
                # Capture exact retraction banner snippet if present in Front
                m = re.search(r".{0,60}retract.{0,120}", front_t, re.I | re.S)
                rec["txt_front_retract_snippet"] = re.sub(r"\s+", " ", m.group(0))[:180] if m else ""

        # Inspect XML using standard JATS LLM parser (<article-title>, <abstract>, <body>, <related-article>)
        if xml_key:
            rx = sess.get(f"https://pmc-oa-opendata.s3.amazonaws.com/{xml_key}", timeout=20)
            if rx.status_code == 200:
                xml_str = rx.text
                rec["xml_total_chars"] = len(xml_str)
                rec["xml_anywhere_has_retract"] = bool(re.search(r"retract", xml_str, re.I))
                try:
                    root = ET.fromstring(rx.content)
                    rec["xml_article_type"] = root.attrib.get("article-type", "")
                    # Related article tags (standard JATS metadata for retractions)
                    rel_arts = root.findall(".//related-article")
                    rec["xml_has_related_retraction_tag"] = any(
                        "retract" in (str(ra.attrib.get("related-article-type", "")) + " " + "".join(ra.itertext())).lower()
                        for ra in rel_arts
                    )
                    # Standard LLM text extraction fields: <article-title>, <abstract>, <body>
                    title_txt = " ".join("".join(el.itertext()) for el in root.findall(".//title-group/article-title"))
                    abs_txt = " ".join("".join(el.itertext()) for el in root.findall(".//abstract"))
                    body_txt = " ".join("".join(el.itertext()) for el in root.findall(".//body"))
                    rec["xml_title_chars"] = len(title_txt)
                    rec["xml_abstract_chars"] = len(abs_txt)
                    rec["xml_body_chars"] = len(body_txt)
                    rec["xml_title_has_retract"] = bool(re.search(r"retract", title_txt, re.I))
                    rec["xml_abstract_has_retract"] = bool(re.search(r"retract", abs_txt, re.I))
                    rec["xml_body_has_retract"] = bool(re.search(r"retract", body_txt, re.I))
                    # Standard LLM ingestion (title + abstract + body): is it COMPLETELY SILENT on retraction?!
                    rec["llm_extracted_text_silent"] = (
                        (len(body_txt) > 500 or len(abs_txt) > 100) and
                        (not rec["xml_title_has_retract"]) and
                        (not rec["xml_abstract_has_retract"]) and
                        (not rec["xml_body_has_retract"])
                    )
                except Exception as e:
                    rec["xml_parse_error"] = str(e)[:120]
    except Exception as e:
        rec["error"] = str(e)[:120]
    return rec

rows = []
with ThreadPoolExecutor(max_workers=12) as ex:
    futs = [ex.submit(audit_one, it) for it in pmc_items]
    for fut in as_completed(futs):
        rows.append(fut.result())

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "pmc_oa_s3_stratified_240_audit.csv"), index=False)

# Compute publication-grade summary
has_xml = df[df["has_xml"] == True].copy()
has_txt = df[df["has_txt"] == True].copy()

summary = {
    "audit_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "n_sampled_pmcids": len(df),
    "n_live_on_s3": int(df["s3_found"].sum()),
    "n_has_txt": len(has_txt),
    "n_has_xml": len(has_xml),
    "txt_front_mentions_retract_n": int(has_txt["txt_front_has_retract"].sum()),
    "txt_front_mentions_retract_pct": round(100 * has_txt["txt_front_has_retract"].mean(), 2) if len(has_txt) else 0,
    "txt_body_silent_no_retract_n": int((has_txt["txt_body_has_retract"] == False).sum()),
    "txt_body_silent_no_retract_pct": round(100 * (has_txt["txt_body_has_retract"] == False).mean(), 2) if len(has_txt) else 0,
    "xml_has_related_retraction_tag_n": int(has_xml["xml_has_related_retraction_tag"].sum()) if "xml_has_related_retraction_tag" in has_xml else 0,
    "xml_has_related_retraction_tag_pct": round(100 * has_xml["xml_has_related_retraction_tag"].mean(), 2) if len(has_xml) else 0,
    "xml_article_type_counts": has_xml["xml_article_type"].value_counts().head(10).to_dict() if "xml_article_type" in has_xml else {},
    "xml_title_silent_pct": round(100 * (has_xml["xml_title_has_retract"] == False).mean(), 2) if len(has_xml) else 0,
    "xml_abstract_silent_pct": round(100 * (has_xml["xml_abstract_has_retract"] == False).mean(), 2) if len(has_xml) else 0,
    "xml_body_silent_pct": round(100 * (has_xml["xml_body_has_retract"] == False).mean(), 2) if len(has_xml) else 0,
    "llm_standard_jats_extractor_completely_silent_n": int(has_xml["llm_extracted_text_silent"].sum()) if "llm_extracted_text_silent" in has_xml else 0,
    "llm_standard_jats_extractor_completely_silent_pct": round(100 * has_xml["llm_extracted_text_silent"].mean(), 2) if len(has_xml) else 0,
    "median_xml_body_chars": int(has_xml["xml_body_chars"].median()) if len(has_xml) else 0,
    "mean_xml_body_chars": int(has_xml["xml_body_chars"].mean()) if len(has_xml) else 0,
    "median_txt_body_chars": int(has_txt["txt_body_chars"].median()) if len(has_txt) else 0,
    "example_txt_front_snippets": has_txt["txt_front_retract_snippet"].dropna().head(5).tolist(),
}
with open(os.path.join(OUT, "analysis2_s3_llm_corpus_audit_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)

print("\n=== ANALYSIS 2 (LIVE AWS S3 PMC-OA-OPENDATA PARSER AUDIT) SUMMARY ===", flush=True)
print(json.dumps(summary, indent=2), flush=True)
