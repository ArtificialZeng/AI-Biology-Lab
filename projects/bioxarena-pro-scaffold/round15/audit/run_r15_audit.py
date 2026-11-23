#!/usr/bin/env python3
"""
Independent Round 15 Confirmatory Audit Runner, Release-Wave IUT Adjudicator, & SA-H1/SA-Panel (A-T3).

Strictly follows locked `PREREG_H.md` (`v2`, addendum to `PREREG_E.md` commit `4c606f4`):
  - Reuses frozen `round12/audit/independent_r12_metrics.py::evaluate_target_line_r12`
    with `DEFAULT_R15_SEED = 20261231`, `n_boot = 4,000`, `n_perm = 100`, two-sided
    Bonferroni `97.5%` percentile CIs (`[1.25%, 98.75%]`, one-sided alpha = 0.0125).
  - Enforces strict blinding gates, open-license + `2027-06-30` sunset check, and hard
    `jurkat` exclusion (`PREREG_H §1–§3`); never reads `round15/src/` or `round15/results/`.
  - Release-wave & overall cross-wave IUT adjudication (`PREREG_H v2 §8.2, §8.4`):
      * Evaluates each open release wave immediately and reports per-wave verdicts
        (labelled `INTERIM_MORE_LINES_PENDING` until all 5 candidate lines are resolved
        or sunset `2027-06-30` is reached).
      * Aggregates qualifying lines across ALL waves (`K_conf_total in 0..5`): confirmation
        requires unanimous `K_conf_total / K_conf_total` (`100%`) IUT passage across every
        qualifying line across all waves (later waves can keep or revoke, never rescue).
      * Enforces mandatory `K_conf = 1` limited generalizability notice AND `K_conf >= 2`
        single-laboratory / single-platform (`Arc VCC 2026, 10x Flex CRISPRi`) notice.
  - Pre-registered §9 descriptive secondary suite (`只描述、不判定`):
      * `SA-H1`: Leave-`vcc2025_h1`-out `5`-source recalculation (`FCT` + 5 paired `H-E1b`).
      * `SA-Panel`: `129` shared validation–test targets vs remaining panel targets.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

import numpy as np

AUDIT_R15_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R15_DIR.parents[1]
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R12_AUDIT_DIR = REPO_ROOT / "round12" / "audit"
R14_AUDIT_DIR = REPO_ROOT / "round14" / "audit"
for _p in (str(AUDIT_R15_DIR), str(R14_AUDIT_DIR), str(R12_AUDIT_DIR), str(R11_AUDIT_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from independent_r11_loader import (  # noqa: E402
    R9_LINE_NAMES,
    R9_LINE_SET,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R10_SOURCE_EXT_NAMES,
    AuditLineData,
    CanonicalHGNC,
    load_line_dataset,
    load_r9_common_genes,
)
from independent_r11_predictors import (  # noqa: E402
    deterministic_subsample,
    sample_covariance,
)
from independent_r12_metrics import (  # noqa: E402
    align_source_and_target_r12,
    evaluate_target_line_r12,
)
from independent_r15_loader import (  # noqa: E402
    DEFAULT_R15_SEED,
    DEV_ALLOWED_LINES_BY_PHASE,
    PREREG_H_SUNSET_DATE,
    R15_CANDIDATE_LINES,
    R15_CANDIDATE_SET,
    R15_EXCLUDED_LINES,
    R15_STATUS_DIR,
    _normalize_line_token,
    enforce_r15_gate,
    parse_confirmatory_lines_r15,
    resolve_r15_line_cache_dir,
)
from run_r14_audit import subset_audit_line_targets  # noqa: E402


def check_r15_line_thresholds(
    target_data: AuditLineData,
    source_list: Sequence[AuditLineData],
    n_nt_half_a: int,
    n_nt_half_b: int,
    only_reproducible: bool = True,
    min_targets: int = 50,
    min_genes: int = 3000,
    min_ntc_half: int = 15,
) -> Dict[str, Any]:
    """Verify `PREREG_H §7` per-line gate (`n_eval_targets >= 50`, `|G_eval| >= 3000`, `per-half NTC >= 15`)."""
    norm = _normalize_line_token(target_data.name)
    if norm in R15_EXCLUDED_LINES:
        raise RuntimeError("PREREG_H §3 hard exclusion: jurkat cannot enter threshold check or evaluation.")

    try:
        cohort = align_source_and_target_r12(source_list, target_data, only_reproducible=only_reproducible)
        n_genes = len(cohort.eval_genes)
        n_tgts = len(cohort.eval_targets)
        ev_tgts = list(cohort.eval_targets)
    except Exception as exc:
        n_genes = 0
        n_tgts = 0
        ev_tgts = []
        align_err = str(exc)
    else:
        align_err = ""

    ok_tgts = n_tgts >= int(min_targets)
    ok_genes = n_genes >= int(min_genes)
    ok_nt_a = int(n_nt_half_a) >= int(min_ntc_half)
    ok_nt_b = int(n_nt_half_b) >= int(min_ntc_half)
    passes = bool((not align_err) and ok_tgts and ok_genes and ok_nt_a and ok_nt_b)

    reasons: List[str] = []
    if align_err:
        reasons.append(align_err)
    if not ok_tgts:
        reasons.append(f"n_eval_targets={n_tgts}<{min_targets}")
    if not ok_genes:
        reasons.append(f"n_genes_eval={n_genes}<{min_genes}")
    if not ok_nt_a:
        reasons.append(f"n_nt_half_a={n_nt_half_a}<{min_ntc_half}")
    if not ok_nt_b:
        reasons.append(f"n_nt_half_b={n_nt_half_b}<{min_ntc_half}")

    return {
        "line": target_data.name,
        "n_genes_eval": int(n_genes),
        "n_eval_targets": int(n_tgts),
        "n_nt_half_a": int(n_nt_half_a),
        "n_nt_half_b": int(n_nt_half_b),
        "passes_thresholds": passes,
        "eval_targets": ev_tgts,
        "reason": ";".join(reasons) if reasons else "OK",
    }


def adjudicate_release_wave(
    wave_id: str,
    per_line_summary: Dict[str, Dict[str, Any]],
    gate_info: Dict[str, Any],
    thresholds_by_line: Optional[Dict[str, Dict[str, Any]]] = None,
    all_resolved_lines_so_far: Optional[Sequence[str]] = None,
    check_date: Optional[_dt.date] = None,
    min_targets: int = 50,
    min_genes: int = 3000,
) -> Dict[str, Any]:
    """
    Adjudicate ONE open release wave immediately under `PREREG_H v2 §7, §8, §8.4`:
      - Validates every confirmatory candidate line in this wave against `>= 50` targets,
        `>= 3,000` genes, and `>= 15` NTC/half (demoting failed lines if thresholds fail).
      - Applies IUT (`100%` passage) within the wave across all qualifying wave lines.
      - Labels output `INTERIM_MORE_LINES_PENDING` whenever `< 5` candidate lines have been
        resolved and `check_date <= 2027-06-30`.
    """
    today = check_date or _dt.date.today()
    thresholds_by_line = thresholds_by_line or {}
    role_by_line = dict(gate_info.get("role_by_line", {}))
    wave_conf_lines: List[str] = []
    demoted_lines: Dict[str, str] = {}

    for ln in gate_info.get("confirmatory_lines", []):
        norm = _normalize_line_token(ln)
        if norm in R15_EXCLUDED_LINES:
            raise RuntimeError("PREREG_H §3 violation: jurkat cannot enter confirmatory adjudication.")
        if norm not in per_line_summary:
            raise RuntimeError(f"Wave {wave_id}: missing summary for confirmatory line {norm!r}")
        res = per_line_summary[norm]
        thr = thresholds_by_line.get(norm, {})
        n_t = int(res.get("n_eval_targets", 0))
        n_g = int(res.get("n_genes_eval", 0))
        thr_ok = thr.get("passes_thresholds", True)
        if n_t < int(min_targets) or n_g < int(min_genes) or not thr_ok:
            demoted_lines[norm] = thr.get("reason", f"n_targets={n_t},n_genes={n_g}")
            role_by_line[norm] = "SECONDARY_ONLY_FAILED_THRESHOLD"
        else:
            wave_conf_lines.append(norm)

    k_wave = len(wave_conf_lines)
    c1_per_line = {ln: bool(per_line_summary[ln]["he1_pass"]) for ln in wave_conf_lines}
    c2_per_line = {ln: bool(per_line_summary[ln]["he2_pass"]) for ln in wave_conf_lines}

    c1_wave_ok = bool(k_wave > 0 and all(c1_per_line.values()))
    c2_wave_ok = bool(k_wave > 0 and all(c2_per_line.values()))

    resolved_set = set(all_resolved_lines_so_far or []) | set(per_line_summary.keys()) | {
        ln for ln, r in role_by_line.items() if r.startswith("EXCLUDED")
    }
    remaining = [c for c in R15_CANDIDATE_LINES if c not in resolved_set]
    is_final_study = (len(remaining) == 0) or (today > PREREG_H_SUNSET_DATE)
    interim_label = "FINAL" if is_final_study else "interim, more lines pending"

    if k_wave == 0:
        wave_status = "CONFIRMATORY_WAVE_ABORTED_K0"
    elif k_wave == 1:
        wave_status = "CONFIRMATORY_WAVE_SINGLE_LINE_K1"
    else:
        wave_status = f"CONFIRMATORY_WAVE_MULTI_LINE_K{k_wave}"

    return {
        "wave_id": wave_id,
        "date": today.isoformat(),
        "status": wave_status,
        "interim_label": interim_label,
        "is_final_study_endpoint": is_final_study,
        "k_conf_wave": k_wave,
        "confirmatory_lines": wave_conf_lines,
        "secondary_eval_lines": [ln for ln in gate_info.get("secondary_eval_lines", []) if ln in per_line_summary] + list(demoted_lines.keys()),
        "role_by_line": role_by_line,
        "demoted_lines": demoted_lines,
        "claim1_he1_per_line": c1_per_line,
        "claim2_he2_per_line": c2_per_line,
        "claim1_confirmed_wave": c1_wave_ok,
        "claim2_confirmed_wave": c2_wave_ok,
        "claim1_he1_verdict": "CONFIRMED" if c1_wave_ok else ("ABORTED_K0" if k_wave == 0 else "NOT_CONFIRMED"),
        "claim2_he2_verdict": "CONFIRMED" if c2_wave_ok else ("ABORTED_K0" if k_wave == 0 else "NOT_CONFIRMED"),
        "remaining_candidate_lines": remaining,
    }


def adjudicate_overall_round15(
    wave_results: Sequence[Dict[str, Any]],
    check_date: Optional[_dt.date] = None,
) -> Dict[str, Any]:
    """
    Aggregate ALL release waves into the pre-registered overall Round 15 IUT verdict (`PREREG_H v2 §8.2, §8.4`):
      - Collects every qualifying confirmatory line across all waves (`K_conf_total in 0..5`).
      - Strictly enforces NO duplicate confirmatory/resolved lines across waves and
        unanimous `K_conf_total / K_conf_total` IUT per claim across ALL waves (later waves
        can only maintain or revoke an earlier confirmation, never rescue an earlier failure).
      - Applies exact `§8.2` mandatory scope wording (`单线，推广性有限` for `K=1`; single-lab/platform
        `Arc Institute VCC 2026, 10x Flex CRISPRi` + shared-panel note for `K>=2`).
    """
    today = check_date or _dt.date.today()
    conf_lines_all: List[str] = []
    seen_conf: Set[str] = set()
    c1_all: Dict[str, bool] = {}
    c2_all: Dict[str, bool] = {}
    resolved_all: Dict[str, str] = {}
    wave_ids: List[str] = []

    for w in wave_results:
        wid = str(w.get("wave_id", f"wave_{len(wave_ids) + 1}"))
        wave_ids.append(wid)
        for ln, role in w.get("role_by_line", {}).items():
            norm = _normalize_line_token(ln)
            if norm not in R15_CANDIDATE_SET:
                continue
            if role in ("NOT_RELEASED", "PENDING"):
                continue
            if norm in resolved_all and role != resolved_all[norm]:
                raise RuntimeError(f"Duplicate/conflicting line {norm!r} across waves: {resolved_all[norm]} vs {role}")
            resolved_all[norm] = str(role)

        for ln in w.get("confirmatory_lines", []):
            norm = _normalize_line_token(ln)
            if norm in R15_EXCLUDED_LINES:
                raise RuntimeError("PREREG_H §3: jurkat forbidden in overall adjudication.")
            if norm in seen_conf:
                raise RuntimeError(f"Duplicate confirmatory line {norm!r} across release waves ({wave_ids}).")
            seen_conf.add(norm)
            conf_lines_all.append(norm)
            c1_all[norm] = bool(w["claim1_he1_per_line"][ln])
            c2_all[norm] = bool(w["claim2_he2_per_line"][ln])

    k_total = len(conf_lines_all)
    remaining = [c for c in R15_CANDIDATE_LINES if c not in resolved_all]
    past_sunset = today > PREREG_H_SUNSET_DATE
    is_final = (len(remaining) == 0) or past_sunset
    interim_label = "FINAL" if is_final else "interim, more lines pending"

    if k_total == 0:
        c1_ok = False
        c2_ok = False
        status = "CONFIRMATORY_EXPIRED_UNRELEASED" if (past_sunset and not wave_results) else "CONFIRMATORY_ABORTED_K0"
        c1_verdict = "EXPIRED_NO_OPEN_RELEASE" if status == "CONFIRMATORY_EXPIRED_UNRELEASED" else "ABORTED_K0"
        c2_verdict = c1_verdict
        scope_msg = (
            f"PREREG_H overall adjudication ({interim_label}, K_conf = 0): zero qualifying lines cleared §7; "
            f"remaining_lines={remaining}. Zero confirmatory claims asserted."
        )
    else:
        c1_ok = all(c1_all.values())
        c2_ok = all(c2_all.values())
        c1_verdict = "CONFIRMED" if c1_ok else "NOT_CONFIRMED"
        c2_verdict = "CONFIRMED" if c2_ok else "NOT_CONFIRMED"
        if k_total == 1:
            status = "CONFIRMATORY_SINGLE_LINE_K1"
            scope_msg = (
                f"MANDATORY SCOPE NOTICE (K_conf = 1, 单线，推广性有限; {interim_label}): "
                f"Confirmatory adjudication covers only a single line ({conf_lines_all[0]}). "
                f"Cross-line generalizability is limited; NEVER claim multi-line validation."
            )
        else:
            status = f"CONFIRMATORY_MULTI_LINE_K{k_total}"
            scope_msg = (
                f"MANDATORY SCOPE NOTICE (K_conf = {k_total}; {interim_label}): Any confirmed claim MAY ONLY "
                f"be described as confirmed across {k_total} unseen human cell lines from a single laboratory "
                f"and platform (Arc Institute VCC 2026, 10x Flex CRISPRi) across lines {conf_lines_all}. "
                f"Validation and test panels share 129 targets (not independent target samples)."
            )

    return {
        "status": status,
        "date": today.isoformat(),
        "sunset_date": PREREG_H_SUNSET_DATE.isoformat(),
        "is_final_study_endpoint": is_final,
        "interim_label": interim_label,
        "wave_ids": wave_ids,
        "k_conf_total": k_total,
        "confirmatory_lines": conf_lines_all,
        "resolved_roles": resolved_all,
        "remaining_candidate_lines": remaining,
        "excluded_lines_permanent": sorted(R15_EXCLUDED_LINES),
        "claim1_he1_per_line": c1_all,
        "claim2_he2_per_line": c2_all,
        "claim1_confirmed": c1_ok,
        "claim2_confirmed": c2_ok,
        "claim1_he1_verdict": c1_verdict,
        "claim2_he2_verdict": c2_verdict,
        "mandatory_scope_disclosure": scope_msg,
        "prior_art_disclosure": (
            "PREREG_H §10.6 (Fan et al. NeurIPS 2026 OpenReview lnpzb1HvwX; Kawsar et al. 2026 "
            "CanPert/CASTER doi:10.64898/2026.11.15.766815): Claim 1 builds on Fan et al. and "
            "Kawsar et al. and MUST NEVER be claimed as first."
        ),
    }


def evaluate_sa_h1_leave_vcc2025_out(
    target_data: AuditLineData,
    source_list: Sequence[AuditLineData],
    only_reproducible: bool = True,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    n_boot: int = 4000,
    n_perm: int = 100,
    seed: int = DEFAULT_R15_SEED,
    precomputed_covs: Optional[Dict[str, np.ndarray]] = None,
    precomputed_gene_orders: Optional[Dict[str, Sequence[str]]] = None,
) -> Dict[str, Any]:
    """
    Pre-registered `PREREG_H §9.1` `SA-H1` (Descriptive ONLY / `只描述、不判定`):
    Removes `vcc2025_h1` (same Arc Institute Flex/UG100/CRISPRi platform) from `source_list`
    and recomputes `FCT`, `5` paired `H-E1b` contrasts, `Perm100`, and `H-E2` (`seed = 20261231`).
    """
    filtered_sources = [s for s in source_list if s.name != "vcc2025_h1" and s.name != target_data.name]
    if len(filtered_sources) < 2:
        return {
            "status": "SKIPPED_INSUFFICIENT_NON_VCC_SOURCES",
            "tag": "SA_H1_LEAVE_VCC2025_OUT_DESCRIPTIVE_ONLY",
            "n_sources": len(filtered_sources),
        }
    res = evaluate_target_line_r12(
        target_data=target_data,
        source_list=filtered_sources,
        only_reproducible=only_reproducible,
        shrink=shrink,
        cov_max_cells=cov_max_cells,
        n_boot=n_boot,
        n_perm=n_perm,
        seed=seed,
        include_sa3_anchors=False,
        include_exploratory_e=False,
        include_n5=False,
        include_n5_split=False,
        include_r19=False,
        include_s7=False,
        precomputed_covs=precomputed_covs,
        precomputed_gene_orders=precomputed_gene_orders,
    )
    res["tag"] = "SA_H1_LEAVE_VCC2025_OUT_DESCRIPTIVE_ONLY"
    res["excluded_source"] = "vcc2025_h1"
    return res


def evaluate_sa_panel_overlap(
    target_data: AuditLineData,
    source_list: Sequence[AuditLineData],
    shared_129_targets: Optional[Set[str]] = None,
    only_reproducible: bool = True,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    n_boot: int = 4000,
    n_perm: int = 100,
    seed: int = DEFAULT_R15_SEED,
    precomputed_covs: Optional[Dict[str, np.ndarray]] = None,
    precomputed_gene_orders: Optional[Dict[str, Sequence[str]]] = None,
) -> Dict[str, Any]:
    """
    Pre-registered `PREREG_H §9.2` `SA-Panel` (Descriptive ONLY / `只描述、不判定`):
    Separately evaluates `eval_targets ∩ shared_129_targets` vs `eval_targets \\ shared_129_targets`.
    """
    if not shared_129_targets:
        return {
            "tag": "SA_PANEL_DESCRIPTIVE_ONLY",
            "status": "PENDING_OPEN_SHARED_129_TARGET_LIST",
            "subgroups": {},
        }

    cohort_all = align_source_and_target_r12(source_list, target_data, only_reproducible=only_reproducible)
    eval_set = set(cohort_all.eval_targets)
    s_shared = eval_set & set(shared_129_targets)
    s_unique = eval_set - set(shared_129_targets)

    cov_map = dict(precomputed_covs) if precomputed_covs is not None else {}
    gord_map = dict(precomputed_gene_orders) if precomputed_gene_orders is not None else {}
    out: Dict[str, Any] = {
        "tag": "SA_PANEL_DESCRIPTIVE_ONLY",
        "status": "OK",
        "parent_line": target_data.name,
        "subgroups": {},
    }

    for sg_name, kset in (("shared_val_test_129", s_shared), ("panel_unique_remaining", s_unique)):
        if len(kset) < 2:
            out["subgroups"][sg_name] = {
                "status": "INSUFFICIENT_SUBGROUP_TARGETS",
                "n_eval_targets": len(kset),
                "eval_targets": sorted(kset),
            }
            continue
        sub_tag = f"{target_data.name}::{sg_name}"
        sub_data = subset_audit_line_targets(target_data, kset, subgroup_name=sub_tag)
        if target_data.name in cov_map:
            cov_map[sub_tag] = cov_map[target_data.name]
        if target_data.name in gord_map:
            gord_map[sub_tag] = gord_map[target_data.name]

        res = evaluate_target_line_r12(
            target_data=sub_data,
            source_list=source_list,
            only_reproducible=only_reproducible,
            shrink=shrink,
            cov_max_cells=cov_max_cells,
            n_boot=n_boot,
            n_perm=n_perm,
            seed=seed,
            include_sa3_anchors=False,
            include_exploratory_e=False,
            include_n5=False,
            include_n5_split=False,
            include_r19=False,
            include_s7=False,
            precomputed_covs=cov_map,
            precomputed_gene_orders=gord_map,
        )
        res["status"] = "OK"
        res["tag"] = "SA_PANEL_DESCRIPTIVE_ONLY"
        res["eval_targets"] = sorted(kset)
        out["subgroups"][sg_name] = res
    return out


def _write_r15_outputs(summary: Dict[str, Any], out_dir: Path, tag: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"audit_r15_{tag}.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    line_rows: List[Dict[str, Any]] = []
    pair_rows: List[Dict[str, Any]] = []
    for ln, d in summary.get("per_line", {}).items():
        fct = d["method_eval"]["fct"]["all_nontarget"]
        sa8_e1 = d["sa8_e1"]
        d_sa8 = d["sa8_e1_split_minus_all"]
        perm_all = d["perm100_all_nontarget"]["r_orth_cal"]
        line_rows.append({
            "line": ln,
            "n_genes_eval": d["n_genes_eval"],
            "n_eval_targets": d["n_eval_targets"],
            "fct_r_orth_cal_est": fct.get("r_orth_cal_est"),
            "fct_r_orth_cal_ci975_lo": fct.get("r_orth_cal_ci975_lo"),
            "fct_r_orth_cal_ci975_hi": fct.get("r_orth_cal_ci975_hi"),
            "fct_r_orth_cal_valid_frac": fct.get("r_orth_cal_valid_frac"),
            "fct_rho_ab_orth": fct.get("rho_ab_orth"),
            "perm100_all_p975": perm_all.get("p975"),
            "he1_pass": d["he1_pass"],
            "sa8_e1_split_est": sa8_e1.get("r_orth_cal_split_est"),
            "sa8_e1_split_minus_all_est": d_sa8.get("diff_est"),
            "he2_pass": d["he2_pass"],
        })
        for src, pr in d["paired_he1b"].items():
            mr = d["mismatched_all"][src]
            pair_rows.append({
                "line": ln,
                "source": src,
                "mis_est": mr.get("r_orth_cal_est"),
                "diff_est": pr.get("diff_est"),
                "diff_ci975_lo": pr.get("diff_ci975_lo"),
                "diff_ci975_hi": pr.get("diff_ci975_hi"),
                "passes_he1b_source": pr.get("passes_he1b_source"),
            })

    for fname, rows in (
        (f"audit_r15_lines_{tag}.csv", line_rows),
        (f"audit_r15_paired_he1b_{tag}.csv", pair_rows),
    ):
        if not rows:
            continue
        keys: List[str] = []
        for r in rows:
            for k in r.keys():
                if k not in keys:
                    keys.append(k)
        with (out_dir / fname).open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Independent Round 15 Release-Wave Audit Runner (`PREREG_H v2`, A-T3)")
    parser.add_argument(
        "--phase",
        choices=("synthetic_test", "dev_r9", "dev_r10", "dev_r11", "dev_r12", "dev_r13", "dev_r14", "holdout_r15"),
        default="dev_r10",
    )
    parser.add_argument("--wave-id", default="wave_1")
    parser.add_argument("--source-mode", choices=("r9_6", "r9_r10_8"), default="r9_6")
    parser.add_argument("--lines", nargs="+", default=None)
    parser.add_argument("--prior-wave-jsons", nargs="*", type=Path, default=())
    parser.add_argument("--shared-129-targets-file", type=Path, default=None)
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--seed", type=int, default=DEFAULT_R15_SEED)
    parser.add_argument("--cov-max-cells", type=int, default=None)
    parser.add_argument("--min-targets", type=int, default=50)
    parser.add_argument("--min-genes", type=int, default=3000)
    parser.add_argument("--no-hgnc", action="store_true")
    parser.add_argument("--skip-sa-h1", action="store_true")
    parser.add_argument("--out-dir", type=Path, default=AUDIT_R15_DIR / "dev_smoke")
    args = parser.parse_args()

    is_holdout = args.phase == "holdout_r15"
    if not is_holdout:
        if args.phase == "synthetic_test":
            raise RuntimeError("Use `test_r15_audit.py` or `xcheck_r15_synthetic.py` for synthetic_test.")
        default_targets = list(DEV_ALLOWED_LINES_BY_PHASE[args.phase])[:1]
        targets = list(args.lines or default_targets)
        gate_meta = enforce_r15_gate(
            args.phase,
            target_lines=targets,
            extra_dirs=[resolve_r15_line_cache_dir(t, args.phase) for t in targets],
        )
        gate_info: Dict[str, Any] = {
            "confirmatory_lines": [],
            "k_conf": 0,
            "secondary_eval_lines": targets,
            "role_by_line": {},
        }
    else:
        gate_info = parse_confirmatory_lines_r15(R15_STATUS_DIR / "CONFIRMATORY_LINES_R15.txt")
        targets = list(args.lines) if args.lines else list(gate_info["confirmatory_lines"] + gate_info["secondary_eval_lines"])
        gate_meta = enforce_r15_gate(
            "holdout_r15",
            target_lines=targets,
            status_dir=R15_STATUS_DIR,
            require_confirmatory_file=True,
        )

    shared_129: Set[str] = set()
    if args.shared_129_targets_file and args.shared_129_targets_file.exists():
        shared_129 = {ln.strip() for ln in args.shared_129_targets_file.read_text(encoding="utf-8").splitlines() if ln.strip()}

    hgnc = None if args.no_hgnc else CanonicalHGNC.from_tsv()
    common_genes = load_r9_common_genes(hgnc=hgnc)
    source_pool_names = list(R9_LINE_NAMES) + (list(R10_SOURCE_EXT_NAMES) if args.source_mode == "r9_r10_8" else [])

    loaded_cache: Dict[str, AuditLineData] = {}
    cov_cache: Dict[str, np.ndarray] = {}
    gene_order_cache: Dict[str, Sequence[str]] = {}
    for sn in source_pool_names:
        sdir = R9_LINES_DIR / sn if sn in R9_LINE_SET else R10_LINES_DIR / sn
        sdata = load_line_dataset(sdir, representation="log1p", aligned_genes=common_genes, hgnc=hgnc, min_gene_overlap=3000)
        loaded_cache[sn] = sdata
        cov_cache[sn] = sample_covariance(deterministic_subsample(sdata.ctrl_cells, args.cov_max_cells, seed=42), shrink=0.0)
        gene_order_cache[sn] = sdata.genes

    per_line: Dict[str, Dict[str, Any]] = {}
    thresholds_by_line: Dict[str, Dict[str, Any]] = {}
    sa_h1_by_line: Dict[str, Dict[str, Any]] = {}
    sa_panel_by_line: Dict[str, Dict[str, Any]] = {}

    for tgt in targets:
        norm = _normalize_line_token(tgt)
        if norm in R15_EXCLUDED_LINES:
            raise RuntimeError("PREREG_H §3 violation: jurkat is strictly excluded.")
        if tgt in loaded_cache:
            tdata = loaded_cache[tgt]
        else:
            min_ov = int(args.min_genes) if (is_holdout and norm in gate_info["confirmatory_lines"]) else 50
            tdata = load_line_dataset(
                resolve_r15_line_cache_dir(tgt, args.phase),
                representation="log1p",
                aligned_genes=common_genes,
                hgnc=hgnc,
                min_gene_overlap=min_ov,
            )
            cov_cache[tgt] = sample_covariance(deterministic_subsample(tdata.ctrl_cells, args.cov_max_cells, seed=42), shrink=0.0)
            gene_order_cache[tgt] = tdata.genes

        nt_a = int(tdata.metadata.get("n_nt_half_a", tdata.ctrl_cells.shape[0] // 2))
        nt_b = int(tdata.metadata.get("n_nt_half_b", tdata.ctrl_cells.shape[0] // 2))
        src_list = [loaded_cache[sn] for sn in source_pool_names if sn != tgt]
        thresholds_by_line[norm] = check_r15_line_thresholds(
            tdata,
            src_list,
            n_nt_half_a=nt_a,
            n_nt_half_b=nt_b,
            min_targets=args.min_targets,
            min_genes=args.min_genes,
        )
        per_line[norm] = evaluate_target_line_r12(
            target_data=tdata,
            source_list=src_list,
            only_reproducible=True,
            cov_max_cells=args.cov_max_cells,
            n_boot=args.n_boot,
            n_perm=args.n_perm,
            seed=args.seed,
            include_sa3_anchors=True,
            include_exploratory_e=True,
            include_n5=True,
            include_n5_split=True,
            include_r19=True,
            include_s7=True,
            precomputed_covs=cov_cache,
            precomputed_gene_orders=gene_order_cache,
        )
        if not args.skip_sa_h1:
            sa_h1_by_line[norm] = evaluate_sa_h1_leave_vcc2025_out(
                target_data=tdata,
                source_list=src_list,
                cov_max_cells=args.cov_max_cells,
                n_boot=args.n_boot,
                n_perm=args.n_perm,
                seed=args.seed,
                precomputed_covs=cov_cache,
                precomputed_gene_orders=gene_order_cache,
            )
        sa_panel_by_line[norm] = evaluate_sa_panel_overlap(
            target_data=tdata,
            source_list=src_list,
            shared_129_targets=shared_129,
            cov_max_cells=args.cov_max_cells,
            n_boot=args.n_boot,
            n_perm=args.n_perm,
            seed=args.seed,
            precomputed_covs=cov_cache,
            precomputed_gene_orders=gene_order_cache,
        )

    if not is_holdout:
        wave_verdict: Dict[str, Any] = {
            "wave_id": args.wave_id,
            "status": f"DEV_DRY_RUN_{args.phase.upper()}",
            "k_conf_wave": 0,
            "confirmatory_lines": [],
            "claim1_he1_verdict": "DEV_ONLY_NO_ADJUDICATION",
            "claim2_he2_verdict": "DEV_ONLY_NO_ADJUDICATION",
        }
        overall_verdict = wave_verdict
    else:
        prior_waves = [json.loads(p.read_text(encoding="utf-8"))["wave_verdict"] for p in args.prior_wave_jsons]
        prior_resolved = {ln for pw in prior_waves for ln, r in pw.get("role_by_line", {}).items() if r != "NOT_RELEASED"}
        wave_verdict = adjudicate_release_wave(
            wave_id=args.wave_id,
            per_line_summary=per_line,
            gate_info=gate_info,
            thresholds_by_line=thresholds_by_line,
            all_resolved_lines_so_far=prior_resolved,
            min_targets=args.min_targets,
            min_genes=args.min_genes,
        )
        overall_verdict = adjudicate_overall_round15(prior_waves + [wave_verdict])

    summary = {
        "phase": args.phase,
        "wave_id": args.wave_id,
        "source_mode": args.source_mode,
        "seed": args.seed,
        "n_boot": args.n_boot,
        "n_perm": args.n_perm,
        "gate_meta": gate_meta,
        "gate_info": gate_info,
        "thresholds": thresholds_by_line,
        "targets": targets,
        "wave_verdict": wave_verdict,
        "overall_verdict": overall_verdict,
        "per_line": per_line,
        "sa_h1_leave_vcc2025_out": sa_h1_by_line,
        "sa_panel": sa_panel_by_line,
    }
    tag = f"{args.phase}_{args.wave_id}_{args.source_mode}"
    _write_r15_outputs(summary, args.out_dir, tag)
    print(json.dumps({"wave": wave_verdict["status"], "overall": overall_verdict.get("status"), "targets": targets}, indent=2))


if __name__ == "__main__":
    main()
