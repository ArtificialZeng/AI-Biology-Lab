#!/usr/bin/env python3
"""
Pilot A-5 (Nature Medicine / Science / Management Science Healthcare Operations):
The Trial-to-Real-World Demographic Exclusion Gap — Linking Pivotal Phase III
Enrolled Older-Adult Age (ClinicalTrials.gov `baselineCharacteristicsModule`)
to Real-World FAERS Older-Adult Utilization Share & Within-Drug Excluded-Stratum
(`>=65` vs. `18-64`) Mortality Disparity (openFDA FAERS).

Uses ONLY official CLI wrappers (`clinical_trials_api.py` and `openfda_query.py`).
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import numpy as np
from scipy import stats

PORT = Path("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026")
DATA_DIR = PORT / "data" / "pilot_a5"
DATA_DIR.mkdir(parents=True, exist_ok=True)
UV = "/Users/mac/.local/bin/uv"
CT_CLI = "/Users/mac/.gemini/config/plugins/scientific-literature/skills/clinical_trials_database/scripts/clinical_trials_api.py"
FDA_CLI = "/Users/mac/.gemini/config/plugins/scientific-literature/skills/openfda_database/scripts/openfda_query.py"

DRUGS = [
    ("PEMBROLIZUMAB", "Oncology"),
    ("OSIMERTINIB", "Oncology"),
    ("PALBOCICLIB", "Oncology"),
    ("IBRUTINIB", "Oncology"),
    ("VENETOCLAX", "Oncology"),
    ("APIXABAN", "Cardiology"),
    ("RIVAROXABAN", "Cardiology"),
    ("EMPAGLIFLOZIN", "Cardiometabolic"),
    ("DUPILUMAB", "Immunology"),
    ("UPADACITINIB", "Immunology"),
    ("ADALIMUMAB", "Immunology"),
    ("OCRELIZUMAB", "Neurology"),
]


def fda_total(search: str, tag: str) -> int:
    out_f = DATA_DIR / f"fda_{tag}.json"
    if not out_f.exists():
        p = subprocess.run(
            [UV, "run", FDA_CLI, "search", "--category", "drug", "--endpoint", "event",
             "--search", search, "--limit", "1", "--output", str(out_f)],
            cwd=str(PORT), capture_output=True, text=True,
        )
        if p.returncode != 0:
            return 0
    return int((json.loads(out_f.read_text()).get("meta") or {}).get("results", {}).get("total", 0))


def extract_trial_age_stats(studies: list) -> dict:
    mean_ages = []
    pct_65_vals = []
    n_max_age_caps = 0
    n_trials_parsed = 0
    for s in studies:
        ps = s.get("protocolSection", {})
        em = ps.get("eligibilityModule", {})
        if em.get("maximumAge"):
            n_max_age_caps += 1

        rs = s.get("resultsSection", {})
        bm = rs.get("baselineCharacteristicsModule", {})
        measures = bm.get("measures", [])
        got_age = False
        for m in measures:
            title = str(m.get("title") or "").lower()
            if "age" not in title:
                continue
            unit = str(m.get("unitOfMeasure") or "").lower()
            param = str(m.get("paramType") or "").upper()
            # 1. Continuous mean or median age in years
            if ("year" in unit or param in ("MEAN", "MEDIAN")) and not got_age:
                vals = []
                for cls in m.get("classes", []):
                    for cat in cls.get("categories", []):
                        for meas in cat.get("measurements", []):
                            v = meas.get("value")
                            if v is not None:
                                try:
                                    fv = float(v)
                                    if 15.0 <= fv <= 95.0:
                                        vals.append(fv)
                                except ValueError:
                                    pass
                if vals:
                    mean_ages.append(float(np.mean(vals)))
                    got_age = True
            # 2. Categorical >=65 share
            if "categor" in title or param in ("COUNT_OF_PARTICIPANTS", "NUMBER"):
                tot_n = 0
                ge65_n = 0
                for cls in m.get("classes", []):
                    for cat in cls.get("categories", []):
                        ct = str(cat.get("title") or "").lower()
                        # Sum across groups (avoiding Total group double-count by taking last if labeled Total or all non-total)
                        cvals = []
                        for meas in cat.get("measurements", []):
                            v = meas.get("value")
                            if v is not None:
                                try:
                                    cvals.append(float(v))
                                except ValueError:
                                    pass
                        if not cvals:
                            continue
                        cat_sum = max(cvals) if len(cvals) > 1 and max(cvals) >= sum(cvals) * 0.45 else sum(cvals)
                        tot_n += cat_sum
                        if re.search(r"(>=|≥|>|65\s*to|65\s*-|65\+|over\s*65|65\s*years\s*and\s*older)", ct):
                            ge65_n += cat_sum
                if tot_n >= 20 and 0 < ge65_n < tot_n:
                    pct_65_vals.append(ge65_n / tot_n)
        if got_age:
            n_trials_parsed += 1
    return {
        "n_p3_fetched": len(studies),
        "n_with_age": n_trials_parsed,
        "n_with_max_age_cap": n_max_age_caps,
        "mean_trial_age": round(float(np.mean(mean_ages)), 2) if mean_ages else None,
        "n_with_cat65": len(pct_65_vals),
        "mean_trial_65plus_share": round(float(np.mean(pct_65_vals)), 4) if pct_65_vals else None,
    }


def main() -> None:
    rows = []
    for drug, area in DRUGS:
        # 1. ClinicalTrials.gov completed Phase 3 interventional trials with results
        ct_f = DATA_DIR / f"ct_p3_{drug}.json"
        if not ct_f.exists():
            p = subprocess.run(
                [UV, "run", CT_CLI, "search",
                 "--intervention", drug,
                 "--study-type", "INTERVENTIONAL",
                 "--phase", "PHASE3",
                 "--status", "COMPLETED",
                 "--has-results",
                 "--fields", "NCTId,BriefTitle,EligibilityModule,resultsSection.baselineCharacteristicsModule",
                 "--limit", "25",
                 "--count-total",
                 "--output", str(ct_f)],
                cwd=str(PORT), capture_output=True, text=True,
            )
            if p.returncode != 0:
                continue
        ct_d = json.loads(ct_f.read_text(encoding="utf-8"))
        t_stats = extract_trial_age_stats(ct_d.get("studies", []))

        # 2. openFDA FAERS 2019-2026: 18-64 vs >=65 reports and within-drug death rates!
        base = f'patient.drug.openfda.substance_name.exact:"{drug}"+AND+receivedate:[20190101+TO+20261124]+AND+patient.patientonsetageunit:801'
        n_adult = fda_total(f"{base}+AND+patient.patientonsetage:[18+TO+64]", f"{drug}_18_64")
        n_older = fda_total(f"{base}+AND+patient.patientonsetage:[65+TO+110]", f"{drug}_65p")
        d_adult = fda_total(f"{base}+AND+patient.patientonsetage:[18+TO+64]+AND+seriousnessdeath:1", f"{drug}_18_64_death")
        d_older = fda_total(f"{base}+AND+patient.patientonsetage:[65+TO+110]+AND+seriousnessdeath:1", f"{drug}_65p_death")

        tot_aged = n_adult + n_older
        faers_65p_share = (n_older / tot_aged) if tot_aged > 0 else None
        dr_adult = (d_adult / n_adult) if n_adult > 0 else None
        dr_older = (d_older / n_older) if n_older > 0 else None
        mort_rr = (dr_older / dr_adult) if (dr_older is not None and dr_adult and dr_adult > 0) else None
        gap_65p = (
            round(faers_65p_share - t_stats["mean_trial_65plus_share"], 4)
            if (faers_65p_share is not None and t_stats["mean_trial_65plus_share"] is not None)
            else None
        )

        rows.append({
            "drug": drug,
            "area": area,
            "ct_p3_total": ct_d.get("totalCount"),
            **t_stats,
            "faers_n_18_64": n_adult,
            "faers_n_65plus": n_older,
            "faers_65plus_share": round(faers_65p_share, 4) if faers_65p_share is not None else None,
            "trial_to_faers_65p_gap": gap_65p,
            "faers_death_rate_18_64": round(dr_adult, 4) if dr_adult is not None else None,
            "faers_death_rate_65plus": round(dr_older, 4) if dr_older is not None else None,
            "within_drug_older_vs_adult_mortality_rr": round(mort_rr, 3) if mort_rr is not None else None,
        })

    # Correlation between Trial mean age / 65+ exclusion gap and within-drug mortality RR
    valid_age = [(r["mean_trial_age"], r["within_drug_older_vs_adult_mortality_rr"]) for r in rows
                 if r["mean_trial_age"] is not None and r["within_drug_older_vs_adult_mortality_rr"] is not None]
    valid_gap = [(r["trial_to_faers_65p_gap"], r["within_drug_older_vs_adult_mortality_rr"]) for r in rows
                 if r["trial_to_faers_65p_gap"] is not None and r["within_drug_older_vs_adult_mortality_rr"] is not None]

    corrs = {}
    if len(valid_age) >= 5:
        r_sp, p_sp = stats.spearmanr([x[0] for x in valid_age], [x[1] for x in valid_age])
        corrs["trial_mean_age_vs_mortality_rr_spearman"] = {"n": len(valid_age), "rho": round(float(r_sp), 3), "p": round(float(p_sp), 4)}
    if len(valid_gap) >= 5:
        r_sp, p_sp = stats.spearmanr([x[0] for x in valid_gap], [x[1] for x in valid_gap])
        corrs["gap_65p_vs_mortality_rr_spearman"] = {"n": len(valid_gap), "rho": round(float(r_sp), 3), "p": round(float(p_sp), 4)}

    out = {"rows": rows, "correlations": corrs}
    (PORT / "pilots" / "pilot_a5_summary.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
