#!/usr/bin/env python3
"""
Pilot B-1 (UTD24 / FT50): Causal Impact of the March 30, 2026 FDA Mass Clinical-Trial
Delinquency Enforcement Sweep (2,200+ Sponsors, 3,000+ Trials) on Results Disclosure &
Selective Reporting (ClinicalTrials.gov API v2).

Uses ONLY the official clinical-trials-database CLI wrapper (`clinical_trials_api.py`).
"""

from __future__ import annotations

import calendar
import json
import subprocess
from pathlib import Path

PORT = Path("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026")
DATA_DIR = PORT / "data" / "pilot_b1"
DATA_DIR.mkdir(parents=True, exist_ok=True)
UV = "/Users/mac/.local/bin/uv"
CT_CLI = "/Users/mac/.gemini/config/plugins/scientific-literature/skills/clinical_trials_database/scripts/clinical_trials_api.py"


def ct_count(advanced: str, tag: str) -> int:
    out_f = DATA_DIR / f"cnt_{tag}.json"
    if not out_f.exists():
        cmd = [
            UV, "run", CT_CLI, "search",
            "--study-type", "INTERVENTIONAL",
            "--has-results",
            "--advanced", advanced,
            "--fields", "NCTId",
            "--limit", "1",
            "--count-total",
            "--output", str(out_f),
        ]
        p = subprocess.run(cmd, cwd=str(PORT), capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"CT count failed for {tag}: {p.stderr[-400:]}")
    d = json.loads(out_f.read_text(encoding="utf-8"))
    return int(d.get("totalCount", 0))


def ct_fetch(advanced: str, fields: str, limit: int, tag: str) -> dict:
    out_f = DATA_DIR / f"fetch_{tag}.json"
    if not out_f.exists():
        cmd = [
            UV, "run", CT_CLI, "search",
            "--study-type", "INTERVENTIONAL",
            "--has-results",
            "--advanced", advanced,
            "--fields", fields,
            "--limit", str(limit),
            "--count-total",
            "--sort", "LastUpdatePostDate:desc",
            "--output", str(out_f),
        ]
        p = subprocess.run(cmd, cwd=str(PORT), capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"CT fetch failed for {tag}: {p.stderr[-400:]}")
    return json.loads(out_f.read_text(encoding="utf-8"))


def month_range(year: int, month: int) -> str:
    last_day = calendar.monthrange(year, month)[1]
    return f"{month:02d}/01/{year}, {month:02d}/{last_day:02d}/{year}"


def main() -> None:
    # Overdue cohort definition: PrimaryCompletionDate strictly <= 12/31/2023 (already >24 months overdue
    # by early 2026 and >12 months overdue throughout 2025-2026 — eliminates mechanical aging confounding!)
    overdue_filter = "AREA[PrimaryCompletionDate]RANGE[01/01/2017, 12/31/2023]"

    arms = {
        "treated_p23_us": f"{overdue_filter} AND (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3) AND AREA[LocationCountry]United States",
        "ctrl_p1_us": f"{overdue_filter} AND AREA[Phase]PHASE1 AND NOT (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3) AND AREA[LocationCountry]United States",
        "ctrl_p23_nonus": f"{overdue_filter} AND (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3) AND NOT AREA[LocationCountry]United States",
        "treated_p23_us_industry": f"{overdue_filter} AND (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3) AND AREA[LocationCountry]United States AND AREA[LeadSponsorClass]INDUSTRY",
        "treated_p23_us_other": f"{overdue_filter} AND (AREA[Phase]PHASE2 OR AREA[Phase]PHASE3) AND AREA[LocationCountry]United States AND AREA[LeadSponsorClass]OTHER",
    }

    # Monthly series: Jan 2025 through Nov 2026 (23 months, centered on March 30, 2026 FDA sweep!)
    months = [(2025, m) for m in range(1, 13)] + [(2026, m) for m in range(1, 12)]
    series = {arm: [] for arm in arms}
    for y, m in months:
        mr = month_range(y, m)
        ym = f"{y}-{m:02d}"
        for arm_name, base_q in arms.items():
            q = f"{base_q} AND AREA[ResultsFirstSubmitDate]RANGE[{mr}]"
            cnt = ct_count(q, f"{arm_name}_{ym}")
            series[arm_name].append({"month": ym, "count": cnt})

    # Fetch 250 overdue Phase 2/3 US trials disclosed POST-sweep (Apr-Nov 2026) vs PRE-sweep (Apr-Nov 2025)
    # to test Lemons / Selective-Hiding Hypothesis (p-values & serious adverse event rates!)
    fields = (
        "NCTId,BriefTitle,OverallStatus,Phase,LeadSponsorName,LeadSponsorClass,"
        "PrimaryCompletionDate,ResultsFirstSubmitDate,EnrollmentCount,"
        "resultsSection.adverseEventsModule.eventGroups,"
        "resultsSection.outcomeMeasuresModule.outcomeMeasures"
    )
    post_q = (
        f"{arms['treated_p23_us']} AND AREA[ResultsFirstSubmitDate]RANGE[04/01/2026, 11/24/2026]"
    )
    pre_q = (
        f"{arms['treated_p23_us']} AND AREA[ResultsFirstSubmitDate]RANGE[04/01/2025, 11/24/2025]"
    )
    post_data = ct_fetch(post_q, fields, 250, "post_2026_sample")
    pre_data = ct_fetch(pre_q, fields, 250, "pre_2025_sample")

    def summarize_sample(studies: list) -> dict:
        n = len(studies)
        sponsor_cls = {}
        terminated = 0
        sae_rates = []
        pvals = []
        sig_p = 0
        non_sig_p = 0
        for s in studies:
            ps = s.get("protocolSection", {})
            st = ps.get("statusModule", {}).get("overallStatus", "")
            if st == "TERMINATED":
                terminated += 1
            sc = ps.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}).get("class", "UNKNOWN")
            sponsor_cls[sc] = sponsor_cls.get(sc, 0) + 1

            rs = s.get("resultsSection", {})
            # Serious AE proportion across event groups
            ae_groups = rs.get("adverseEventsModule", {}).get("eventGroups", [])
            tot_risk = sum(int(g.get("seriousNumAtRisk", 0) or 0) for g in ae_groups)
            tot_aff = sum(int(g.get("seriousNumAffected", 0) or 0) for g in ae_groups)
            if tot_risk > 0:
                sae_rates.append(tot_aff / tot_risk)

            # Primary outcome p-values
            oms = rs.get("outcomeMeasuresModule", {}).get("outcomeMeasures", [])
            trial_ps = []
            for om in oms:
                if om.get("type") != "PRIMARY":
                    continue
                for an in om.get("analyses", []):
                    pv_raw = str(an.get("pValue", "")).strip()
                    if not pv_raw:
                        continue
                    cleaned = pv_raw.replace("<", "").replace(">", "").replace("=", "").strip()
                    try:
                        val = float(cleaned)
                        if 0.0 <= val <= 1.0:
                            trial_ps.append(val)
                    except ValueError:
                        pass
            if trial_ps:
                min_p = min(trial_ps)
                pvals.append(min_p)
                if min_p < 0.05:
                    sig_p += 1
                else:
                    non_sig_p += 1

        return {
            "n_fetched": n,
            "sponsor_class_counts": sponsor_cls,
            "terminated_share": round(terminated / n, 4) if n else None,
            "n_with_sae": len(sae_rates),
            "mean_sae_rate": round(sum(sae_rates) / len(sae_rates), 4) if sae_rates else None,
            "n_with_primary_pval": len(pvals),
            "sig_p_lt_005_share": round(sig_p / len(pvals), 4) if pvals else None,
            "non_sig_p_ge_005_share": round(non_sig_p / len(pvals), 4) if pvals else None,
        }

    # Compute DiD summary: Pre (2025-04..2025-11 vs 2025-01..2026-02) vs Post (2026-04..2026-10 complete months)
    def avg_months(arm_list, ym_set):
        vals = [x["count"] for x in arm_list if x["month"] in ym_set]
        return sum(vals) / len(vals) if vals else 0.0

    pre_yoy = {f"2025-{m:02d}" for m in range(4, 11)}   # Apr-Oct 2025 (7 full months)
    post_yoy = {f"2026-{m:02d}" for m in range(4, 11)}  # Apr-Oct 2026 (7 full months post-March 30 sweep!)
    pre_imm = {"2025-11", "2025-12", "2026-01", "2026-02"}  # 4 months immediately before sweep

    did_summary = {}
    for arm_name, arm_list in series.items():
        m_pre_yoy = avg_months(arm_list, pre_yoy)
        m_pre_imm = avg_months(arm_list, pre_imm)
        m_post = avg_months(arm_list, post_yoy)
        peak_post = max((x for x in arm_list if x["month"] >= "2026-04"), key=lambda x: x["count"])
        did_summary[arm_name] = {
            "mean_apr_oct_2025": round(m_pre_yoy, 2),
            "mean_nov25_feb26": round(m_pre_imm, 2),
            "mean_apr_oct_2026": round(m_post, 2),
            "ratio_post_vs_imm_pre": round(m_post / m_pre_imm, 3) if m_pre_imm > 0 else None,
            "ratio_post_vs_yoy_2025": round(m_post / m_pre_yoy, 3) if m_pre_yoy > 0 else None,
            "peak_post_month": peak_post,
        }

    out = {
        "monthly_series": series,
        "did_summary": did_summary,
        "post_2026_total_matching": post_data.get("totalCount"),
        "pre_2025_total_matching": pre_data.get("totalCount"),
        "post_2026_sample_stats": summarize_sample(post_data.get("studies", [])),
        "pre_2025_sample_stats": summarize_sample(pre_data.get("studies", [])),
    }
    (PORT / "pilots" / "pilot_b1_summary.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(out["did_summary"], indent=2, ensure_ascii=False))
    print("POST_SAMPLE:", json.dumps(out["post_2026_sample_stats"], indent=2))
    print("PRE_SAMPLE:", json.dumps(out["pre_2025_sample_stats"], indent=2))


if __name__ == "__main__":
    main()
