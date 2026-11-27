#!/usr/bin/env python3
"""
Pilot B-2 (Marketing Science / ISR / Management Science):
Viral Consumer Demand, Legal Salience, and Algorithmic Pharmacovigilance Distortion
in openFDA FAERS (Semaglutide & Tirzepatide vs. Matched Clinical Controls, 2019-2026).

Uses ONLY the official openfda-database CLI wrapper (`openfda_query.py`), keeping total
API calls strictly < 60 (well below the keyless 1,000/day ceiling).
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

PORT = Path("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026")
DATA_DIR = PORT / "data" / "pilot_b2"
DATA_DIR.mkdir(parents=True, exist_ok=True)
UV = "/Users/mac/.local/bin/uv"
FDA_CLI = "/Users/mac/.gemini/config/plugins/scientific-literature/skills/openfda_database/scripts/openfda_query.py"


def fda_count(search: str, count_field: str, tag: str) -> dict:
    out_f = DATA_DIR / f"cnt_{tag}.json"
    cmd = [
        UV, "run", FDA_CLI, "count",
        "--category", "drug", "--endpoint", "event",
        "--search", search,
        "--count_field", count_field,
        "--summary", "25",
        "--output", str(out_f),
    ]
    p = subprocess.run(cmd, cwd=str(PORT), capture_output=True, text=True)
    if p.returncode != 0:
        return {"error": p.stderr[-300:], "results": []}
    return json.loads(out_f.read_text(encoding="utf-8"))


def fda_total(search: str, tag: str) -> int:
    out_f = DATA_DIR / f"tot_{tag}.json"
    cmd = [
        UV, "run", FDA_CLI, "search",
        "--category", "drug", "--endpoint", "event",
        "--search", search,
        "--limit", "1",
        "--output", str(out_f),
    ]
    p = subprocess.run(cmd, cwd=str(PORT), capture_output=True, text=True)
    if p.returncode != 0:
        return 0
    d = json.loads(out_f.read_text(encoding="utf-8"))
    return int((d.get("meta") or {}).get("results", {}).get("total", 0))


def main() -> None:
    # Viral shock drugs vs active-comparator T2D/metabolic controls
    drugs = {
        "SEMAGLUTIDE": "viral",
        "TIRZEPATIDE": "viral",
        "EMPAGLIFLOZIN": "control",
        "SITAGLIPTIN": "control",
        "DULAGLUTIDE": "control",
    }
    periods = {
        "pre_2019_2021": "[20190101+TO+20211231]",
        "post_2023_2026": "[20230101+TO+20261124]",
    }

    # 1. Qualification breakdown (1=Physician, 2=Pharmacist, 3=Other HP, 4=Lawyer, 5=Consumer)
    qual_map = {"1": "physician", "2": "pharmacist", "3": "other_hp", "4": "lawyer", "5": "consumer"}
    panel = {}
    for drug, role in drugs.items():
        panel[drug] = {"role": role, "periods": {}}
        for pname, drange in periods.items():
            q = f'patient.drug.openfda.generic_name.exact:"{drug}"+AND+receivedate:{drange}'
            res_q = fda_count(q, "primarysource.qualification", f"{drug}_{pname}_qual")
            counts = {qual_map.get(str(r["term"]), str(r["term"])): int(r["count"]) for r in res_q.get("results", [])}
            tot = sum(counts.values())
            # Serious outcome (death) by consumer vs physician
            death_all = fda_total(f"{q}+AND+seriousnessdeath:1", f"{drug}_{pname}_death_all")
            death_md = fda_total(f"{q}+AND+primarysource.qualification:1+AND+seriousnessdeath:1", f"{drug}_{pname}_death_md")
            death_cons = fda_total(f"{q}+AND+primarysource.qualification:5+AND+seriousnessdeath:1", f"{drug}_{pname}_death_cons")
            # Media/litigation-salient reaction (Gastroparesis) by lawyer+consumer vs physician
            gp_md = fda_total(f'{q}+AND+primarysource.qualification:1+AND+patient.reaction.reactionmeddrapt.exact:"Gastroparesis"', f"{drug}_{pname}_gp_md")
            gp_law = fda_total(f'{q}+AND+primarysource.qualification:4+AND+patient.reaction.reactionmeddrapt.exact:"Gastroparesis"', f"{drug}_{pname}_gp_law")
            gp_cons = fda_total(f'{q}+AND+primarysource.qualification:5+AND+patient.reaction.reactionmeddrapt.exact:"Gastroparesis"', f"{drug}_{pname}_gp_cons")
            # Severe clinical event (Pancreatitis) by physician vs pooled
            panc_md = fda_total(f'{q}+AND+primarysource.qualification:1+AND+patient.reaction.reactionmeddrapt.exact:"Pancreatitis"', f"{drug}_{pname}_panc_md")
            panc_all = fda_total(f'{q}+AND+patient.reaction.reactionmeddrapt.exact:"Pancreatitis"', f"{drug}_{pname}_panc_all")

            n_md = counts.get("physician", 0)
            n_cons = counts.get("consumer", 0)
            n_law = counts.get("lawyer", 0)
            panel[drug]["periods"][pname] = {
                "qual_counts": counts,
                "total_with_qual": tot,
                "consumer_plus_lawyer_share": round((n_cons + n_law) / tot, 4) if tot else None,
                "lawyer_share": round(n_law / tot, 4) if tot else None,
                "physician_share": round(n_md / tot, 4) if tot else None,
                "death_all": death_all,
                "death_md": death_md,
                "death_cons": death_cons,
                "death_rate_all": round(death_all / tot, 5) if tot else None,
                "death_rate_md": round(death_md / n_md, 5) if n_md else None,
                "death_rate_cons": round(death_cons / n_cons, 5) if n_cons else None,
                "gastroparesis_md": gp_md,
                "gastroparesis_lawyer": gp_law,
                "gastroparesis_consumer": gp_cons,
                "gastroparesis_rate_md": round(gp_md / n_md, 5) if n_md else None,
                "gastroparesis_rate_lawyer": round(gp_law / n_law, 5) if n_law else None,
                "pancreatitis_all": panc_all,
                "pancreatitis_md": panc_md,
                "pancreatitis_rate_all": round(panc_all / tot, 5) if tot else None,
                "pancreatitis_rate_md": round(panc_md / n_md, 5) if n_md else None,
            }

    (PORT / "pilots" / "pilot_b2_summary.json").write_text(
        json.dumps(panel, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(panel, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
