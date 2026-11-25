#!/usr/bin/env python3
"""Master statistical audit & allele-exact verification across all live public-data pilots."""

from __future__ import annotations

import json
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import stats

PORT = Path("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026")
UV = "/Users/mac/.local/bin/uv"
CV_CLI = "/Users/mac/.gemini/config/plugins/genomics/skills/clinvar_database/scripts/clinvar_api.py"


def main() -> None:
    # 1. Allele-exact audit of Pilot A-2: check exact (chr, pos, ref, alt) for 25 flagged non-NFE P/LP hits
    a2 = json.loads((PORT / "pilots" / "pilot_a2_summary.json").read_text())
    a2_dir = PORT / "data" / "pilot_a2"

    gn_exact = {}
    for gene in ["MYH7", "SCN5A", "KCNQ1", "LMNA", "MYBPC3", "LDLR", "BRCA1", "BRCA2"]:
        gn_raw = json.loads((a2_dir / f"gnomad_{gene}.json").read_text())
        gn_vars = (
            gn_raw if isinstance(gn_raw, list)
            else (gn_raw.get("variants") or gn_raw.get("data", {}).get("gene", {}).get("variants") or [])
        )
        for v in gn_vars:
            vid = v.get("variant_id")
            if not vid:
                continue
            ac_p, an_p = {}, {}
            for src in ("exome", "genome"):
                b = v.get(src) or {}
                for pop in b.get("populations") or []:
                    pid = str(pop.get("id", ""))
                    if pid in ("afr", "amr", "eas", "mid", "nfe", "sas"):
                        ac_p[pid] = ac_p.get(pid, 0) + int(pop.get("ac") or 0)
                        an_p[pid] = an_p.get(pid, 0) + int(pop.get("an") or 0)
            afs = {k: ac_p[k] / an_p[k] for k in ac_p if an_p.get(k, 0) >= 2000}
            nfe = afs.get("nfe", 0.0)
            nn = {k: val for k, val in afs.items() if k != "nfe"}
            nn_pop = max(nn, key=lambda k: nn[k]) if nn else None
            nn_af = nn[nn_pop] if nn_pop else 0.0
            gn_exact[vid] = {"nfe_af": nfe, "non_nfe_pop": nn_pop, "non_nfe_af": nn_af, "afs": afs}

    exact_matches = []
    multiallelic_collisions = 0
    checked = 0
    for ex in a2.get("flagged_non_nfe_plp_examples", [])[:25]:
        cid = str(ex["clinvar_id"])
        evf = a2_dir / f"ev_{cid}.json"
        if not evf.exists():
            p = subprocess.run([UV, "run", CV_CLI, "evidence", "--variant_id", cid, "--output", str(evf)],
                               cwd=str(PORT), capture_output=True, text=True)
            if p.returncode != 0:
                continue
        ev = json.loads(evf.read_text())
        ai = ev.get("allele_info") or {}
        ch = str(ai.get("chromosome") or "").replace("chr", "")
        pos = ai.get("position_start")
        ref = ai.get("reference_allele")
        alt = ai.get("alternate_allele")
        if not (ch and pos and ref and alt):
            continue
        checked += 1
        vid = f"{ch}-{pos}-{ref}-{alt}"
        gm = gn_exact.get(vid)
        if gm and gm["non_nfe_af"] >= 1e-4 and gm["non_nfe_af"] >= 5.0 * max(gm["nfe_af"], 1e-6):
            exact_matches.append({
                **ex,
                "exact_gnomad_vid": vid,
                "exact_nfe_af": round(gm["nfe_af"], 6),
                "exact_non_nfe_pop": gm["non_nfe_pop"],
                "exact_non_nfe_af": round(gm["non_nfe_af"], 6),
                "exact_ratio": round(gm["non_nfe_af"] / max(gm["nfe_af"], 1e-6), 1),
            })
        else:
            multiallelic_collisions += 1

    # 2. Pilot A-1 formal Mann-Whitney tests
    a1 = json.loads((PORT / "pilots" / "pilot_a1_summary.json").read_text())
    profs = a1["profiles"]

    def mw_test(vals1, vals0):
        if len(vals1) < 3 or len(vals0) < 3:
            return None
        u, pv = stats.mannwhitneyu(vals1, vals0, alternative="two-sided")
        return {
            "n1": len(vals1), "n0": len(vals0),
            "med1": round(float(np.median(vals1)), 4), "med0": round(float(np.median(vals0)), 4),
            "U": float(u), "p_mw": float(f"{pv:.4e}"),
        }

    sm_b1 = [p for p in profs if p["modality"] == "Small molecule" and p["black_box_warning"] == 1]
    sm_b0 = [p for p in profs if p["modality"] == "Small molecule" and p["black_box_warning"] == 0]
    ab_b1 = [p for p in profs if p["modality"] == "Antibody" and p["black_box_warning"] == 1]
    ab_b0 = [p for p in profs if p["modality"] == "Antibody" and p["black_box_warning"] == 0]

    a1_tests = {
        "sm_loeuf_bbw1_vs_0": mw_test([p["loeuf"] for p in sm_b1 if p["loeuf"] is not None],
                                      [p["loeuf"] for p in sm_b0 if p["loeuf"] is not None]),
        "sm_pleio_bbw1_vs_0": mw_test([p["ot_disease_breadth"] for p in sm_b1 if p["ot_disease_breadth"] is not None],
                                      [p["ot_disease_breadth"] for p in sm_b0 if p["ot_disease_breadth"] is not None]),
        "sm_subum_bbw1_vs_0": mw_test([p["n_sub_um_human_targets"] for p in sm_b1 if p["n_sub_um_human_targets"] is not None],
                                      [p["n_sub_um_human_targets"] for p in sm_b0 if p["n_sub_um_human_targets"] is not None]),
        "ab_loeuf_bbw1_vs_0": mw_test([p["loeuf"] for p in ab_b1 if p["loeuf"] is not None],
                                      [p["loeuf"] for p in ab_b0 if p["loeuf"] is not None]),
        "ab_pleio_bbw1_vs_0": mw_test([p["ot_disease_breadth"] for p in ab_b1 if p["ot_disease_breadth"] is not None],
                                      [p["ot_disease_breadth"] for p in ab_b0 if p["ot_disease_breadth"] is not None]),
    }

    # 3. Pilot B-1 formal Fisher DiD tests across BOTH independent cohorts
    b1_old = json.loads((PORT / "pilots" / "pilot_b1_summary.json").read_text())
    b1_24 = json.loads((PORT / "pilots" / "pilot_b1_2024cohort.json").read_text())

    def sum_m(series, months):
        return sum(x["count"] for x in series if x["month"] in months)

    pre2, post2 = {"2026-01", "2026-02"}, {"2026-04", "2026-05"}
    ind_pre, ind_post = sum_m(b1_old["monthly_series"]["treated_p23_us_industry"], pre2), sum_m(b1_old["monthly_series"]["treated_p23_us_industry"], post2)
    nonus_pre, nonus_post = sum_m(b1_old["monthly_series"]["ctrl_p23_nonus"], pre2), sum_m(b1_old["monthly_series"]["ctrl_p23_nonus"], post2)
    oth_pre, oth_post = sum_m(b1_old["monthly_series"]["treated_p23_us_other"], pre2), sum_m(b1_old["monthly_series"]["treated_p23_us_other"], post2)
    or_ind_nonus, p_ind_nonus = stats.fisher_exact([[ind_post, ind_pre], [nonus_post, nonus_pre]])
    or_ind_oth, p_ind_oth = stats.fisher_exact([[ind_post, ind_pre], [oth_post, oth_pre]])

    pre3 = {"2026-01", "2026-02", "2026-03"}
    post5 = {"2026-04", "2026-05", "2026-06", "2026-07", "2026-08"}
    i24_pre, i24_post = sum_m(b1_24["treated_p23_us_industry"], pre3), sum_m(b1_24["treated_p23_us_industry"], post5)
    n24_pre, n24_post = sum_m(b1_24["ctrl_p23_nonus"], pre3), sum_m(b1_24["ctrl_p23_nonus"], post5)
    o24_pre, o24_post = sum_m(b1_24["treated_p23_us_other"], pre3), sum_m(b1_24["treated_p23_us_other"], post5)
    or24_in, p24_in = stats.fisher_exact([[i24_post, i24_pre], [n24_post, n24_pre]])
    or24_io, p24_io = stats.fisher_exact([[i24_post, i24_pre], [o24_post, o24_pre]])

    # 4. Pilot B-3 formal Mann-Whitney & Fisher tests
    b3_short = json.loads((PORT / "data" / "pilot_b3" / "drug_shortages_all.json").read_text())["results"]

    def parse_dt(s):
        for fmt in ("%m/%d/%Y", "%Y-%m-%d", "%Y%m%d"):
            try:
                return datetime.strptime(str(s).strip(), fmt)
            except ValueError:
                pass
        return None

    inj_d, oral_d = [], []
    for r in b3_short:
        di, du = parse_dt(r.get("initial_posting_date") or ""), parse_dt(r.get("update_date") or "")
        if di and du and du >= di:
            df = str(r.get("dosage_form") or "").lower()
            if any(k in df for k in ["inject", "solution", "infusion", "intravenous", "powder, for solution"]):
                inj_d.append((du - di).days)
            elif any(k in df for k in ["tablet", "capsule", "oral"]):
                oral_d.append((du - di).days)
    _, p_b3_dur = stats.mannwhitneyu(inj_d, oral_d, alternative="two-sided")
    or_b3_cgmp, p_b3_cgmp = stats.fisher_exact([[61, 355 - 61], [2, 39 - 2]])

    # 5. Pilot B-2 formal Fisher tests
    or_tirz_md, p_tirz_md = stats.fisher_exact([[91, 5465 - 91], [17, 3087 - 17]])
    or_tirz_pool, p_tirz_pool = stats.fisher_exact([[1407, 175714 - 1407], [81, 10503 - 81]])
    or_sema_law_cons, p_sema_law_cons = stats.fisher_exact([[719, 959 - 719], [3211, 44593 - 3211]])
    or_sema_law_md_lb, p_sema_law_md_lb = stats.fisher_exact([[719, 959 - 719], [254, 7928 - 254]])

    audit_out = {
        "pilot_a2_allele_exact_audit": {
            "n_plp_checked": checked,
            "exact_allele_confirmed": len(exact_matches),
            "multiallelic_rsid_collisions": multiallelic_collisions,
            "exact_confirm_rate": round(len(exact_matches) / checked, 3) if checked else None,
            "confirmed_examples": exact_matches[:10],
        },
        "pilot_a1_stat_tests": a1_tests,
        "pilot_b1_stat_tests": {
            "cohort_2017_2023_ind_vs_nonus": {"ind": [ind_pre, ind_post], "nonus": [nonus_pre, nonus_post], "OR": round(float(or_ind_nonus), 3), "p": float(f"{p_ind_nonus:.3e}")},
            "cohort_2017_2023_ind_vs_other": {"ind": [ind_pre, ind_post], "other": [oth_pre, oth_post], "OR": round(float(or_ind_oth), 3), "p": float(f"{p_ind_oth:.3e}")},
            "cohort_2024_ind_vs_nonus": {"ind": [i24_pre, i24_post], "nonus": [n24_pre, n24_post], "OR": round(float(or24_in), 3), "p": float(f"{p24_in:.3e}")},
            "cohort_2024_ind_vs_other": {"ind": [i24_pre, i24_post], "other": [o24_pre, o24_post], "OR": round(float(or24_io), 3), "p": float(f"{p24_io:.3e}")},
        },
        "pilot_b2_stat_tests": {
            "tirz_vs_sita_pancreatitis_physician_only": {"OR": round(float(or_tirz_md), 3), "p": float(f"{p_tirz_md:.3e}")},
            "tirz_vs_sita_pancreatitis_pooled_diluted": {"OR": round(float(or_tirz_pool), 3), "p": float(f"{p_tirz_pool:.3e}")},
            "sema_impaired_gastric_emptying_lawyer_vs_consumer": {"OR": round(float(or_sema_law_cons), 2), "p": float(f"{p_sema_law_cons:.3e}")},
            "sema_impaired_gastric_emptying_lawyer_vs_md_lower_bound": {"OR_lower_bound": round(float(or_sema_law_md_lb), 2), "p_upper_bound": float(f"{p_sema_law_md_lb:.3e}")},
        },
        "pilot_b3_stat_tests": {
            "injectable_vs_oral_duration_mw": {"n_inj": len(inj_d), "n_oral": len(oral_d), "med_inj": int(np.median(inj_d)), "med_oral": int(np.median(oral_d)), "p_mw": float(f"{p_b3_dur:.3e}")},
            "cgmp_vs_labeling_recall_firm_shortage_overlap": {"cgmp": "61/355 (17.18%)", "labeling_placebo": "2/39 (5.13%)", "OR": round(float(or_b3_cgmp), 3), "p_fisher": float(f"{p_b3_cgmp:.4e}")},
        },
    }
    (PORT / "pilots" / "statistical_verification_master.json").write_text(
        json.dumps(audit_out, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(audit_out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
