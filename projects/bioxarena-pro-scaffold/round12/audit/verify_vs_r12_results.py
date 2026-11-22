#!/usr/bin/env python3
"""
Clean-Room Numerical Consistency Verifier (`round12/audit/` vs `round12/results/` Outputs).

Compares strictly OLD-DATA (`dev_r9`, `dev_r10`, `dev_r11`) CSV/JSON evaluation artifacts
produced independently by `round12/audit/run_r12_audit.py` against `round12/results/`
outputs across shared primary endpoints (`FCT`, `E1`, `Mis_s`, `H-E1b` paired diffs,
`H-E2a E1(split)`, `H-E2b E1(split)-E1(all)`, `Perm100`, `N_5`, `SA-Src8`), asserting:
  - `max_point_abs_diff <= 1e-10`
  - `max_ci_abs_diff <= 1e-10`
  - `100%` verdict parity across every sub-hypothesis and source contrast.

Strictly NEVER imports or inspects `round12/src/` source code.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

METHOD_NAME_AUDIT_TO_RESULTS = {
    "fct": "FCT",
    "e1_const_med": "E1_ConstMed",
    "source_grand_mean": "Source_Grand_Mean",
    "transfer_mean_var_scaled": "Transfer_Mean_VarScaled",
    "hybrid": "Hybrid_ProjShared",
    "e2_const_p10": "E2_ConstP10",
    "e2_const_p90": "E2_ConstP90",
    "e3_src_mean_ctrl": "E3_SrcMeanCtrl",
    "e5_rank_matched_mu": "E5_RankMatchedMu",
    "e6_log_fc_bar": "E6_LogFCBar",
    "e7_fc_minus_1": "E7_FCBarMinus1",
    "transfer_mean_global_scale": "Transfer_Mean_GlobalScale",
}


def verify_pair(audit_json: Path, results_dir: Path, tol: float = 1e-10) -> Dict[str, Any]:
    a = json.loads(audit_json.read_text(encoding="utf-8"))
    summary_p = next(results_dir.glob("summary_*.json"))
    adj_p = next(results_dir.glob("adjudication_*.json"))
    lm_p = next(results_dir.glob("line_metrics_*.csv"))
    ci_p = next(results_dir.glob("bootstrap_ci_*.csv"))

    r_sum = json.loads(summary_p.read_text(encoding="utf-8"))
    r_adj = json.loads(adj_p.read_text(encoding="utf-8"))

    lm_map: Dict[Tuple[str, str], Dict[str, str]] = {}
    with open(lm_p, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            lm_map[(row["target_line"], row["method"])] = row

    ci_map: Dict[Tuple[str, str], Dict[str, str]] = {}
    with open(ci_p, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            ci_map[(row["target_line"], row["comparison"])] = row

    point_diffs: List[Tuple[str, float]] = []
    ci_diffs: List[Tuple[str, float]] = []
    verdict_mismatches: List[str] = []

    def _chk_pt(label: str, va: Any, vb: Any) -> None:
        if va is None or vb in (None, ""):
            return
        fa, fb = float(va), float(vb)
        point_diffs.append((label, abs(fa - fb)))

    def _chk_ci(label: str, alo: Any, ahi: Any, blo: Any, bhi: Any) -> None:
        if None in (alo, ahi) or blo in (None, "") or bhi in (None, ""):
            return
        ci_diffs.append((f"{label}:lo", abs(float(alo) - float(blo))))
        ci_diffs.append((f"{label}:hi", abs(float(ahi) - float(bhi))))

    for ln, ad in a["per_line"].items():
        # 1. Core `FCT`, `E1`, `Mis_s`, `TMVS`, `Hybrid`, `SGM` under `all_nontarget`
        for mk_a, mk_r in METHOD_NAME_AUDIT_TO_RESULTS.items():
            if mk_a in ad["method_eval"] and (ln, mk_r) in lm_map:
                ma = ad["method_eval"][mk_a]["all_nontarget"]
                mr = lm_map[(ln, mk_r)]
                _chk_pt(f"{ln}.{mk_r}.r_orth_cal", ma.get("r_orth_cal_est"), mr.get("r_orth_cal"))
                _chk_ci(
                    f"{ln}.{mk_r}.r_orth_cal_ci975",
                    ma.get("r_orth_cal_ci975_lo"),
                    ma.get("r_orth_cal_ci975_hi"),
                    mr.get("r_orth_cal_ci975_lo"),
                    mr.get("r_orth_cal_ci975_hi"),
                )

        for src in ad["source_names"]:
            rm = f"NegCtrl_FCT_MismatchedMean({src})"
            if src in ad["mismatched_all"] and (ln, rm) in lm_map:
                ma = ad["mismatched_all"][src]
                mr = lm_map[(ln, rm)]
                _chk_pt(f"{ln}.Mis({src}).r_orth_cal", ma.get("r_orth_cal_est"), mr.get("r_orth_cal"))
                _chk_ci(
                    f"{ln}.Mis({src}).ci975",
                    ma.get("r_orth_cal_ci975_lo"),
                    ma.get("r_orth_cal_ci975_hi"),
                    mr.get("r_orth_cal_ci975_lo"),
                    mr.get("r_orth_cal_ci975_hi"),
                )

            # H-E1b paired diffs
            ck = f"H-E1b_FCT_minus_Mis_{src}_all_nontarget"
            if src in ad["paired_he1b"] and (ln, ck) in ci_map:
                pa = ad["paired_he1b"][src]
                pr = ci_map[(ln, ck)]
                _chk_pt(f"{ln}.HE1b({src}).est", pa.get("diff_est"), pr.get("estimate"))
                _chk_ci(
                    f"{ln}.HE1b({src}).ci975",
                    pa.get("diff_ci975_lo"),
                    pa.get("diff_ci975_hi"),
                    pr.get("ci975_lo"),
                    pr.get("ci975_hi"),
                )

        # 2. SA-8 (`H-E2a` & `H-E2b`) + in-sample `E1(resp_top50)`
        if (ln, "E1_ConstMed") in lm_map:
            r_e1 = lm_map[(ln, "E1_ConstMed")]
            _chk_pt(f"{ln}.HE2a.E1_split", ad["sa8_e1"].get("r_orth_cal_split_est"), r_e1.get("resp_top50_split_r_orth_cal"))
            _chk_ci(
                f"{ln}.HE2a.E1_split_ci975",
                ad["sa8_e1"].get("r_orth_cal_split_ci975_lo"),
                ad["sa8_e1"].get("r_orth_cal_split_ci975_hi"),
                r_e1.get("resp_top50_split_r_orth_cal_ci975_lo"),
                r_e1.get("resp_top50_split_r_orth_cal_ci975_hi"),
            )
        if (ln, "FCT") in lm_map:
            _chk_pt(f"{ln}.SA8.FCT_split", ad["sa8_fct"].get("r_orth_cal_split_est"), lm_map[(ln, "FCT")].get("resp_top50_split_r_orth_cal"))

        if (ln, "H-E2b_E1_split_minus_all") in ci_map:
            cr = ci_map[(ln, "H-E2b_E1_split_minus_all")]
            _chk_pt(f"{ln}.HE2b.est", ad["sa8_e1_split_minus_all"].get("diff_est"), cr.get("estimate"))
            _chk_ci(
                f"{ln}.HE2b.ci975",
                ad["sa8_e1_split_minus_all"].get("diff_ci975_lo"),
                ad["sa8_e1_split_minus_all"].get("diff_ci975_hi"),
                cr.get("ci975_lo"),
                cr.get("ci975_hi"),
            )

        if (ln, "H-D3_FCT_minus_E1_all_nontarget") in ci_map:
            cr = ci_map[(ln, "H-D3_FCT_minus_E1_all_nontarget")]
            _chk_pt(f"{ln}.HD3.est", ad["pair_fct_minus_e1_all"].get("diff_est"), cr.get("estimate"))
            _chk_ci(
                f"{ln}.HD3.ci975",
                ad["pair_fct_minus_e1_all"].get("diff_ci975_lo"),
                ad["pair_fct_minus_e1_all"].get("diff_ci975_hi"),
                cr.get("ci975_lo"),
                cr.get("ci975_hi"),
            )

        # 3. Perm100 (`all_nontarget`) & N_5 (`FCT`, `GlobalScale`, `SGM`, `E1`, `Mis_s`, paired diffs)
        if ln in r_sum.get("fct_perm100", {}):
            rp = r_sum["fct_perm100"][ln]["all_nontarget"]["r_orth_cal"]
            ap = ad["perm100_all_nontarget"]["r_orth_cal"]
            _chk_pt(f"{ln}.perm100.p975", ap.get("p975"), rp.get("p97_5"))
            _chk_pt(f"{ln}.perm100.median", ap.get("median"), rp.get("median"))

        n5_r = r_sum.get("n5_cone_metrics", {}).get(ln, {})
        n5_a = ad.get("n5_all", {})
        if n5_r and n5_a:
            for mk_a, mk_r in METHOD_NAME_AUDIT_TO_RESULTS.items():
                if mk_a in n5_a.get("methods", {}) and mk_r in n5_r.get("methods", {}):
                    ma = n5_a["methods"][mk_a]
                    mr = n5_r["methods"][mk_r]
                    _chk_pt(f"{ln}.N5.{mk_r}.r_cone_sp", ma.get("r_cone_sp_est"), mr["r_cone_sp"].get("estimate"))
                    _chk_pt(f"{ln}.N5.{mk_r}.eta", ma.get("eta_est"), mr.get("eta"))
                    _chk_ci(
                        f"{ln}.N5.{mk_r}.ci975",
                        ma.get("r_cone_sp_ci975_lo"),
                        ma.get("r_cone_sp_ci975_hi"),
                        mr["r_cone_sp"].get("ci975_lo"),
                        mr["r_cone_sp"].get("ci975_hi"),
                    )
            for src, pa in n5_a.get("paired_fct_minus_mis", {}).items():
                if src in n5_r.get("paired_fct_minus_mis", {}):
                    pr = n5_r["paired_fct_minus_mis"][src]
                    _chk_pt(f"{ln}.N5.paired({src}).est", pa.get("diff_est"), pr.get("estimate"))
                    _chk_ci(
                        f"{ln}.N5.paired({src}).ci975",
                        pa.get("diff_ci975_lo"),
                        pa.get("diff_ci975_hi"),
                        pr.get("ci975_lo"),
                        pr.get("ci975_hi"),
                    )

        # 4. Verdict parity across adjudication JSON
        desc = r_adj.get("all_lines_descriptive", {})
        if ln in desc.get("H-E1", {}).get("per_line", {}):
            re1 = desc["H-E1"]["per_line"][ln]
            for a_k, r_k in (("he1a_pass", "H-E1a_pass"), ("he1b_pass", "H-E1b_pass")):
                if bool(ad[a_k]) != bool(re1.get(r_k)):
                    verdict_mismatches.append(f"{ln}.{a_k}: audit={ad[a_k]} results={re1.get(r_k)}")
        if ln in desc.get("H-E2", {}).get("per_line", {}):
            re2 = desc["H-E2"]["per_line"][ln]
            for a_k, r_k in (("he2a_pass", "H-E2a_pass"), ("he2b_pass", "H-E2b_pass")):
                if bool(ad[a_k]) != bool(re2.get(r_k)):
                    verdict_mismatches.append(f"{ln}.{a_k}: audit={ad[a_k]} results={re2.get(r_k)}")

    max_pt_name, max_pt = max(point_diffs, key=lambda x: x[1]) if point_diffs else ("NONE", float("nan"))
    max_ci_name, max_ci = max(ci_diffs, key=lambda x: x[1]) if ci_diffs else ("NONE", float("nan"))
    passed = bool(
        len(point_diffs) > 0
        and len(ci_diffs) > 0
        and max_pt <= tol
        and max_ci <= tol
        and len(verdict_mismatches) == 0
    )
    return {
        "results_dir": str(results_dir),
        "n_point_checks": len(point_diffs),
        "n_ci_checks": len(ci_diffs),
        "max_point_abs_diff": max_pt,
        "max_point_metric": max_pt_name,
        "max_ci_abs_diff": max_ci,
        "max_ci_metric": max_ci_name,
        "verdict_mismatches": verdict_mismatches,
        "pass_tol": passed,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--audit-json", type=Path, required=True)
    p.add_argument("--results-dir", type=Path, required=True)
    p.add_argument("--tol", type=float, default=1e-10)
    args = p.parse_args()
    out = verify_pair(args.audit_json, args.results_dir, tol=args.tol)
    print(json.dumps(out, indent=2))
    if not out["pass_tol"]:
        sys.exit(1)
