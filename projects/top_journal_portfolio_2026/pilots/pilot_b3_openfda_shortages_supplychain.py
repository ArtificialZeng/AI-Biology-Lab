#!/usr/bin/env python3
"""
Pilot B-3 (Operations Research / MSOM / Management Science):
Pharmaceutical Supply-Chain Bottlenecks & Asymmetric Shortage Persistence in the
openFDA `drug/shortages` Universe Linked to `drug/ndc` Market Concentration &
`drug/enforcement` cGMP vs. Labeling Recalls.

Uses ONLY the official openfda-database CLI wrapper (`openfda_query.py`).
"""

from __future__ import annotations

import json
import subprocess
from collections import Counter
from datetime import datetime
from pathlib import Path

PORT = Path("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026")
DATA_DIR = PORT / "data" / "pilot_b3"
DATA_DIR.mkdir(parents=True, exist_ok=True)
UV = "/Users/mac/.local/bin/uv"
FDA_CLI = "/Users/mac/.gemini/config/plugins/scientific-literature/skills/openfda_database/scripts/openfda_query.py"


def main() -> None:
    # 1. Download full openFDA drug/shortages snapshot (usually ~1,500-2,500 records, 2-3 API calls!)
    short_f = DATA_DIR / "drug_shortages_all.json"
    p = subprocess.run(
        [UV, "run", FDA_CLI, "download", "--category", "drug", "--endpoint", "shortages",
         "--limit", "1000", "--max_pages", "5", "--output", str(short_f)],
        cwd=str(PORT), capture_output=True, text=True,
    )
    if p.returncode != 0:
        raise RuntimeError(f"shortages download failed: {p.stderr[-400:]}")
    short_data = json.loads(short_f.read_text(encoding="utf-8"))
    records = short_data.get("results", [])

    # 2. Download recent drug/enforcement recalls (2023-2026) — 3 pages (3,000 records, 3 API calls!)
    enf_f = DATA_DIR / "drug_enforcement_2023_2026.json"
    p2 = subprocess.run(
        [UV, "run", FDA_CLI, "download", "--category", "drug", "--endpoint", "enforcement",
         "--search", "report_date:[20230101+TO+20261124]",
         "--limit", "1000", "--max_pages", "4", "--output", str(enf_f)],
        cwd=str(PORT), capture_output=True, text=True,
    )
    enf_records = json.loads(enf_f.read_text(encoding="utf-8")).get("results", []) if p2.returncode == 0 and enf_f.exists() else []

    # 3. Summarize shortages schema, statuses, dosage forms, reasons, durations, and firm-level concentration
    status_cnt = Counter()
    reason_cnt = Counter()
    dosage_cnt = Counter()
    by_generic = {}
    durations_by_form = {"injectable": [], "oral": [], "other": []}

    def parse_dt(s: str):
        if not s:
            return None
        for fmt in ("%m/%d/%Y", "%Y-%m-%d", "%Y%m%d"):
            try:
                return datetime.strptime(s.strip(), fmt)
            except ValueError:
                pass
        return None

    for r in records:
        st = str(r.get("status", "UNKNOWN"))
        status_cnt[st] += 1
        rsn = str(r.get("shortage_reason") or "UNSPECIFIED")
        reason_cnt[rsn] += 1
        df = str(r.get("dosage_form") or "UNKNOWN")
        dosage_cnt[df] += 1

        gen = str(r.get("generic_name") or "UNKNOWN").upper().strip()
        comp = str(r.get("company_name") or "UNKNOWN").upper().strip()
        by_generic.setdefault(gen, {"statuses": set(), "companies": set(), "dosage_forms": set(), "n_ndcs": 0})
        by_generic[gen]["statuses"].add(st)
        by_generic[gen]["companies"].add(comp)
        by_generic[gen]["dosage_forms"].add(df)
        by_generic[gen]["n_ndcs"] += 1

        d_init = parse_dt(str(r.get("initial_posting_date") or ""))
        d_upd = parse_dt(str(r.get("update_date") or ""))
        if d_init and d_upd and d_upd >= d_init:
            days = (d_upd - d_init).days
            dfl = df.lower()
            if any(k in dfl for k in ["inject", "solution", "infusion", "intravenous", "powder, for solution"]):
                durations_by_form["injectable"].append(days)
            elif any(k in dfl for k in ["tablet", "capsule", "oral"]):
                durations_by_form["oral"].append(days)
            else:
                durations_by_form["other"].append(days)

    # 4. Summarize recalls: cGMP / sterility vs labeling across Class I/II/III + firm overlap with shortages
    recall_cls = Counter()
    cgmp_firms = set()
    labeling_only_firms = set()
    all_recall_firms = Counter()
    cgmp_kw = ("cgmp", "steril", "microbial", "particulate", "impurity", "subpotent", "superpotent", "dissolution", "nitrosamine", "endotoxin")
    for r in enf_records:
        recall_cls[str(r.get("classification", "UNKNOWN"))] += 1
        firm = str(r.get("recalling_firm") or "").upper().strip()
        all_recall_firms[firm] += 1
        reason = str(r.get("reason_for_recall") or "").lower()
        if any(k in reason for k in cgmp_kw):
            cgmp_firms.add(firm)
        elif "label" in reason or "packag" in reason:
            labeling_only_firms.add(firm)
    labeling_only_firms -= cgmp_firms

    shortage_firms = {str(r.get("company_name") or "").upper().strip() for r in records if r.get("company_name")}

    def fuzzy_overlap(firm_set_a, firm_set_b):
        # Match by first 2 normalized alphabetic tokens to handle Inc./LLC punctuation differences
        def tok(s):
            t = [w for w in "".join(ch if ch.isalnum() else " " for ch in s).split() if w not in {"INC", "LLC", "CORP", "CORPORATION", "LTD", "USA", "US", "THE", "PHARMA", "PHARMACEUTICAL", "PHARMACEUTICALS", "LABORATORIES", "LABS"}]
            return tuple(t[:2]) if t else ()
        tb = {tok(x) for x in firm_set_b if tok(x)}
        return sum(1 for x in firm_set_a if tok(x) in tb)

    # 5. Query openFDA drug/ndc active labeler counts for 12 shortage generics (6 injectable, 6 oral)
    sample_gens = sorted(by_generic.keys(), key=lambda g: -by_generic[g]["n_ndcs"])[:16]
    ndc_conc = []
    for gen in sample_gens[:12]:
        clean_g = gen.split(",")[0].split(";")[0].split(" AND ")[0].strip()
        if not clean_g or len(clean_g) < 3:
            continue
        nf = DATA_DIR / f"ndc_{clean_g[:20].replace(' ', '_').replace('/', '_')}.json"
        pn = subprocess.run(
            [UV, "run", FDA_CLI, "count", "--category", "drug", "--endpoint", "ndc",
             "--search", f'generic_name:"{clean_g}"', "--count_field", "labeler_name.exact",
             "--summary", "100", "--output", str(nf)],
            cwd=str(PORT), capture_output=True, text=True,
        )
        n_labelers = None
        if pn.returncode == 0 and nf.exists():
            nd = json.loads(nf.read_text(encoding="utf-8"))
            n_labelers = len(nd.get("results", []))
        ndc_conc.append({
            "generic": gen,
            "shortage_ndc_rows": by_generic[gen]["n_ndcs"],
            "shortage_companies": len(by_generic[gen]["companies"]),
            "active_ndc_labelers_sampled": n_labelers,
            "statuses": sorted(by_generic[gen]["statuses"]),
            "dosage_forms": sorted(by_generic[gen]["dosage_forms"])[:3],
        })

    def med(arr):
        if not arr:
            return None
        s = sorted(arr)
        return s[len(s) // 2]

    summary = {
        "n_shortage_records_downloaded": len(records),
        "n_unique_generics_in_shortage": len(by_generic),
        "n_unique_shortage_companies": len(shortage_firms),
        "status_counts": dict(status_cnt.most_common(10)),
        "top_reasons": dict(reason_cnt.most_common(10)),
        "top_dosage_forms": dict(dosage_cnt.most_common(10)),
        "duration_days_by_form": {
            k: {"n": len(v), "median_days": med(v), "mean_days": round(sum(v) / len(v), 1) if v else None}
            for k, v in durations_by_form.items()
        },
        "n_enforcement_records_2023_2026": len(enf_records),
        "enforcement_classes": dict(recall_cls),
        "cgmp_recall_firms": len(cgmp_firms),
        "labeling_only_recall_firms": len(labeling_only_firms),
        "cgmp_firms_in_shortage_overlap": fuzzy_overlap(cgmp_firms, shortage_firms),
        "labeling_only_firms_in_shortage_overlap": fuzzy_overlap(labeling_only_firms, shortage_firms),
        "top_shortage_generics_ndc_concentration": ndc_conc,
    }
    (PORT / "pilots" / "pilot_b3_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
