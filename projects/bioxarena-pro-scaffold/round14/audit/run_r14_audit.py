#!/usr/bin/env python3
"""
Independent Round 14 Confirmatory Audit Runner & §7 Subgroup Breakdown (A-T3).

Strictly follows locked `PREREG_G.md` (`commit d34560d`):
  - Reuses frozen `round12/audit/independent_r12_metrics.py::evaluate_target_line_r12`
    with `DEFAULT_R14_SEED = 20261224`, `n_boot = 4,000`, `n_perm = 100`, two-sided
    Bonferroni `97.5%` percentile CIs (`[1.25%, 98.75%]`, one-sided alpha = 0.0125).
  - Enforces strict blinding gates (`synthetic_test`, `dev_r9..dev_r13`,
    `reader_xcheck_r14`, `holdout_r14`) and NEVER imports/reads `round14/src/` or
    `round14/results/`.
  - Confirmatory adjudication: strictly `max K_conf = 1` on `a549_myllia` (`GSE337804`,
    same-lineage independent-lab replication in `A549`, mandatory `n_eval_targets >= 50`
    & `n_genes_eval >= 3,000` gate; strictly NO cross-lineage generalization and NO
    cross-round pooling with Round 12 `a549_ko`).
  - §7 Descriptive Subgroup Analysis (`只报告、不判定`):
      * `overlap_r12_a549` (15 targets shared with R12 `a549_ko`)
      * `novel_not_in_r12_a549` (46 novel `u6_rep ∩ R9_COMMON` targets absent from R12 `a549_ko`)
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

AUDIT_R14_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R14_DIR.parents[1]
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R12_AUDIT_DIR = REPO_ROOT / "round12" / "audit"
for _p in (str(AUDIT_R14_DIR), str(R12_AUDIT_DIR), str(R11_AUDIT_DIR)):
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
from independent_r14_loader import (  # noqa: E402
    DEV_ALLOWED_LINES_BY_PHASE,
    R14_STATUS_DIR,
    R14_TARGET_KEY,
    enforce_r14_gate,
    parse_confirmatory_lines_r14,
    resolve_r14_line_cache_dir,
)

DEFAULT_R14_SEED = 20261224

# Exact pre-registered §7 subgroups from locked PREREG_G.md (15 + 46 = 61 u6_rep ∩ R9_COMMON targets)
OVERLAP_R12_A549: Tuple[str, ...] = (
    "BMS1",
    "CMTR1",
    "CNOT10",
    "CREBBP",
    "DHDDS",
    "EIF4G1",
    "GMPPB",
    "HDAC3",
    "KDM6A",
    "PTEN",
    "RIOK2",
    "SEC61A1",
    "SLC35B1",
    "SMC3",
    "SPCS3",
)

NOVEL_NOT_IN_R12_A549: Tuple[str, ...] = (
    "BRPF1",
    "CALR",
    "CBX3",
    "CTCF",
    "DDX1",
    "DHX29",
    "DHX36",
    "DOT1L",
    "EHMT1",
    "EP300",
    "EWSR1",
    "FASTKD5",
    "GSK3B",
    "HDAC1",
    "HSD17B10",
    "INTS14",
    "IRF2BP2",
    "KAT6A",
    "KDM1A",
    "KDM2B",
    "KDM3B",
    "KEAP1",
    "MSL1",
    "MTPAP",
    "NEDD4",
    "PIK3R4",
    "PRORP",
    "RAF1",
    "RPP14",
    "RPUSD3",
    "RPUSD4",
    "RRP9",
    "SHOC2",
    "SIN3B",
    "SMARCA5",
    "SRF",
    "SRFBP1",
    "SUPT20H",
    "SUPV3L1",
    "TFAM",
    "TFPT",
    "TRMT10C",
    "UBA3",
    "UBA5",
    "UTP6",
    "ZNF598",
)


def _fmt(val: Any, digits: int = 4) -> str:
    try:
        fv = float(val)
        return f"{fv:+.{digits}f}" if np.isfinite(fv) else "nan"
    except (TypeError, ValueError):
        return str(val)


def subset_audit_line_targets(
    target_data: AuditLineData,
    keep_targets: Set[str],
    subgroup_name: Optional[str] = None,
) -> AuditLineData:
    """Clone `AuditLineData` keeping ONLY targets in `keep_targets` (preserving `genes` & NTC cells)."""
    row_idx = np.asarray(
        [i for i, t in enumerate(target_data.target_names) if t in keep_targets],
        dtype=np.int64,
    )
    kept_names = [target_data.target_names[int(i)] for i in row_idx]
    return AuditLineData(
        name=subgroup_name or target_data.name,
        genes=list(target_data.genes),
        target_names=kept_names,
        target_indices=np.asarray(target_data.target_indices[row_idx], dtype=np.int64),
        deltas=np.asarray(target_data.deltas[row_idx], dtype=np.float64),
        ctrl_cells=target_data.ctrl_cells,
        ctrl_mean_libnorm=target_data.ctrl_mean_libnorm,
        ctrl_mean_libnorm_full=target_data.ctrl_mean_libnorm_full,
        repr_mask=np.asarray(target_data.repr_mask[row_idx], dtype=bool) if target_data.repr_mask is not None else None,
        half_a_deltas=(
            np.asarray(target_data.half_a_deltas[row_idx], dtype=np.float64)
            if target_data.half_a_deltas is not None
            else None
        ),
        half_b_deltas=(
            np.asarray(target_data.half_b_deltas[row_idx], dtype=np.float64)
            if target_data.half_b_deltas is not None
            else None
        ),
        half_valid=np.asarray(target_data.half_valid[row_idx], dtype=bool) if target_data.half_valid is not None else None,
        n_cells=np.asarray(target_data.n_cells[row_idx], dtype=np.float64) if target_data.n_cells is not None else None,
        kd_depth_log1p=np.asarray(target_data.kd_depth_log1p[row_idx], dtype=np.float64),
        metadata={**dict(target_data.metadata), "subgroup_subset": subgroup_name or "full"},
    )


def check_r14_target_thresholds(
    target_data: AuditLineData,
    source_list: Sequence[AuditLineData],
    only_reproducible: bool = True,
    min_targets: int = 50,
    min_genes: int = 3000,
) -> Dict[str, Any]:
    """Verify pre-registered `PREREG_G §5` eligibility gates (`n_eval_targets >= 50`, `n_genes_eval >= 3000`)."""
    try:
        cohort = align_source_and_target_r12(source_list, target_data, only_reproducible=only_reproducible)
        n_genes = len(cohort.eval_genes)
        n_tgts = len(cohort.eval_targets)
        ev_tgts = list(cohort.eval_targets)
        passes = (n_genes >= int(min_genes)) and (n_tgts >= int(min_targets))
        reason = "OK" if passes else f"below_threshold(n_targets={n_tgts}<{min_targets}_or_n_genes={n_genes}<{min_genes})"
    except Exception as exc:
        n_genes = 0
        n_tgts = 0
        ev_tgts = []
        passes = False
        reason = str(exc)

    return {
        "line": target_data.name,
        "n_genes_eval": int(n_genes),
        "n_eval_targets": int(n_tgts),
        "min_genes": int(min_genes),
        "min_targets": int(min_targets),
        "passes_thresholds": bool(passes),
        "eval_targets": ev_tgts,
        "reason": reason,
    }


def adjudicate_round14(
    per_line_summary: Dict[str, Dict[str, Any]],
    gate_info: Dict[str, Any],
    min_targets: int = 50,
    min_genes: int = 3000,
) -> Dict[str, Any]:
    """
    Apply `PREREG_G.md §0, §5, §6` single-line confirmatory adjudication (`max K_conf = 1` on `a549_myllia`):
      - Aborts confirmatory adjudication (`K_conf = 0`) if `CONFIRMATORY_LINES_R14.txt`
        fails gate, or `n_eval_targets < 50`, or `n_genes_eval < 3,000`.
      - When `K_conf = 1`, adjudicates Claim 1 (`H-E1`: `H-E1a` & all 6 `H-E1b` paired `FCT - Mis_s`
        & `H-E1c`) and Claim 2 (`H-E2`: `H-E2a` & `H-E2b` on `SA-8`).
      - Strictly enforces same-lineage scope (`A549` only, NO cross-lineage claim, NO R12 pooling).
    """
    conf_lines = list(gate_info.get("confirmatory_lines", []))
    role = str(gate_info.get("a549_myllia_role", "FAILED_GATE"))

    if not conf_lines or R14_TARGET_KEY not in conf_lines:
        return {
            "status": "CONFIRMATORY_ABORTED_K0",
            "k_conf": 0,
            "confirmatory_lines": [],
            "a549_myllia_role": role,
            "claim1_he1_verdict": "ABORTED_MANDATORY_A549_MYLLIA_FAILED_GATE",
            "claim2_he2_verdict": "ABORTED_MANDATORY_A549_MYLLIA_FAILED_GATE",
            "claim1_confirmed": False,
            "claim2_confirmed": False,
            "same_line_replication_only": True,
            "limited_generalizability_notice": True,
            "mandatory_disclosure": (
                f"Round 14 confirmatory adjudication aborted (K_conf = 0; role={role}). "
                f"Descriptive numbers (if ingested) carry zero confirmatory claim; never pool with R12 a549_ko."
            ),
        }

    if R14_TARGET_KEY not in per_line_summary:
        raise RuntimeError(f"Missing confirmatory line summary for {R14_TARGET_KEY}")

    res = per_line_summary[R14_TARGET_KEY]
    n_t = int(res.get("n_eval_targets", 0))
    n_g = int(res.get("n_genes_eval", 0))
    if n_t < int(min_targets) or n_g < int(min_genes):
        return {
            "status": "CONFIRMATORY_ABORTED_K0",
            "k_conf": 0,
            "confirmatory_lines": [],
            "a549_myllia_role": "FAILED_COUNT_THRESHOLD",
            "n_eval_targets": n_t,
            "n_genes_eval": n_g,
            "claim1_he1_verdict": "ABORTED_BELOW_TARGET_OR_GENE_MINIMUM",
            "claim2_he2_verdict": "ABORTED_BELOW_TARGET_OR_GENE_MINIMUM",
            "claim1_confirmed": False,
            "claim2_confirmed": False,
            "same_line_replication_only": True,
            "limited_generalizability_notice": True,
            "mandatory_disclosure": (
                f"Round 14 confirmatory adjudication aborted (K_conf = 0): a549_myllia yielded "
                f"n_eval_targets={n_t} (req >= {min_targets}), n_genes_eval={n_g} (req >= {min_genes}). "
                f"Strictly no cross-round rescue pooling with R12 a549_ko."
            ),
        }

    c1_ok = bool(res["he1_pass"])
    c2_ok = bool(res["he2_pass"])
    if c1_ok and c2_ok:
        interp = (
            "REPLICATED_INDEPENDENT_LAB_WITHIN_A549: Replicates Claim 1 (H-E1) & Claim 2 (H-E2) "
            "in an independent-lab A549 CRISPRko cohort (SpCas9 + 10x 5' HT + DepleteX, 75% novel "
            "targets/guides). Scope strictly restricted to A549; does NOT prove cross-lineage generalization."
        )
    else:
        interp = (
            "FAILED_INDEPENDENT_A549_REPLICATION: Failed Claim 1 and/or Claim 2 in independent-lab "
            f"a549_myllia (H-E1={c1_ok}, H-E2={c2_ok}); platform/lab factors or target shift in A549 not replicated."
        )

    return {
        "status": "CONFIRMATORY_SINGLE_LINE_K1",
        "k_conf": 1,
        "confirmatory_lines": [R14_TARGET_KEY],
        "a549_myllia_role": "CONFIRMATORY",
        "n_eval_targets": n_t,
        "n_genes_eval": n_g,
        "claim1_he1_per_line": {R14_TARGET_KEY: c1_ok},
        "claim2_he2_per_line": {R14_TARGET_KEY: c2_ok},
        "claim1_confirmed": c1_ok,
        "claim2_confirmed": c2_ok,
        "claim1_he1_verdict": "CONFIRMED" if c1_ok else "NOT_CONFIRMED",
        "claim2_he2_verdict": "CONFIRMED" if c2_ok else "NOT_CONFIRMED",
        "same_line_replication_only": True,
        "limited_generalizability_notice": True,
        "scientific_interpretation": interp,
        "mandatory_disclosure": (
            "MANDATORY SCOPE NOTICE (K_conf = 1, PREREG_G §0.1/§0.2/§6): Evaluates ONLY independent-lab, "
            "cross-platform (SpCas9 + 10x 5' HT + DepleteX) replication within A549 (CRISPRko). "
            "Cannot extrapolate to multi-lineage generalization; §7 overlap vs novel split is descriptive only; "
            "strictly forbidden from cross-round pooling with Round 12 a549_ko."
        ),
    }


def evaluate_r14_subgroups(
    target_data: AuditLineData,
    source_list: Sequence[AuditLineData],
    only_reproducible: bool = True,
    shrink: float = 0.0,
    cov_max_cells: Optional[int] = None,
    n_boot: int = 4000,
    n_perm: int = 100,
    seed: int = DEFAULT_R14_SEED,
    precomputed_covs: Optional[Dict[str, np.ndarray]] = None,
    precomputed_gene_orders: Optional[Dict[str, Sequence[str]]] = None,
) -> Dict[str, Any]:
    """
    Pre-registered `PREREG_G.md §7` Descriptive Subgroup Breakdown (`只报告、不判定`, NO multiplicity):
      - `overlap_r12_a549`: targets present in R12 `a549_ko` (15 expected in `u6_rep ∩ R9_COMMON`)
      - `novel_not_in_r12_a549`: remaining novel `u6_rep ∩ R9_COMMON` targets (46 expected)
    Computes `FCT r_orth_cal [97.5% CI]`, all 6 `H-E1b` paired `FCT - Mis_s [97.5% CI]`,
    `Perm100`, `E1_ConstMed`, and `SA-8` strictly descriptively under `seed = 20261224`.
    """
    overlap_set = set(OVERLAP_R12_A549)
    novel_set = set(NOVEL_NOT_IN_R12_A549)

    cov_map = dict(precomputed_covs) if precomputed_covs is not None else {}
    gord_map = dict(precomputed_gene_orders) if precomputed_gene_orders is not None else {}

    subgroups_spec: List[Tuple[str, Set[str]]] = [
        ("overlap_r12_a549", overlap_set),
        ("novel_not_in_r12_a549", novel_set),
    ]

    out: Dict[str, Any] = {
        "tag": "DESCRIPTIVE_SUBGROUP_ONLY_NO_ADJUDICATION",
        "parent_line": target_data.name,
        "seed": int(seed),
        "n_boot": int(n_boot),
        "n_perm": int(n_perm),
        "subgroups": {},
    }

    for sub_name, kset in subgroups_spec:
        sub_tag = f"{target_data.name}::{sub_name}"
        sub_data = subset_audit_line_targets(target_data, kset, subgroup_name=sub_tag)
        if target_data.name in cov_map:
            cov_map[sub_tag] = cov_map[target_data.name]
        if target_data.name in gord_map:
            gord_map[sub_tag] = gord_map[target_data.name]

        try:
            cohort = align_source_and_target_r12(source_list, sub_data, only_reproducible=only_reproducible)
            n_eligible = len(cohort.eval_targets)
            elig_names = list(cohort.eval_targets)
        except Exception as exc:
            n_eligible = 0
            elig_names = []
            err_msg = str(exc)
        else:
            err_msg = ""

        if n_eligible < 2:
            out["subgroups"][sub_name] = {
                "status": "INSUFFICIENT_SUBGROUP_TARGETS",
                "tag": "DESCRIPTIVE_SUBGROUP_ONLY_NO_ADJUDICATION",
                "n_eval_targets": int(n_eligible),
                "eval_targets": elig_names,
                "reason": err_msg or f"n_eval_targets={n_eligible} < 2",
            }
            continue

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
        res["tag"] = "DESCRIPTIVE_SUBGROUP_ONLY_NO_ADJUDICATION"
        res["eval_targets"] = elig_names
        out["subgroups"][sub_name] = res

    return out


def _write_r14_outputs(summary: Dict[str, Any], out_dir: Path, tag: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    json_p = out_dir / f"audit_r14_{tag}.json"
    json_p.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    line_rows: List[Dict[str, Any]] = []
    pair_rows: List[Dict[str, Any]] = []
    sub_rows: List[Dict[str, Any]] = []

    for ln, d in summary["per_line"].items():
        fct = d["method_eval"]["fct"]["all_nontarget"]
        e1_all = d["method_eval"]["e1_const_med"]["all_nontarget"]
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
            "he1a_pass": d["he1a_pass"],
            "he1b_pass": d["he1b_pass"],
            "he1b_n_pass": d["he1b_n_pass"],
            "he1c_pass": d["he1c_pass"],
            "he1_pass": d["he1_pass"],
            "sa8_e1_split_est": sa8_e1.get("r_orth_cal_split_est"),
            "sa8_e1_split_ci975_lo": sa8_e1.get("r_orth_cal_split_ci975_lo"),
            "sa8_e1_split_ci975_hi": sa8_e1.get("r_orth_cal_split_ci975_hi"),
            "sa8_e1_split_minus_all_est": d_sa8.get("diff_est"),
            "sa8_e1_split_minus_all_ci975_lo": d_sa8.get("diff_ci975_lo"),
            "sa8_e1_split_minus_all_ci975_hi": d_sa8.get("diff_ci975_hi"),
            "he2_pass": d["he2_pass"],
            "e1_all_est": e1_all.get("r_orth_cal_est"),
            "pair_fct_minus_e1_all_est": d["pair_fct_minus_e1_all"].get("diff_est"),
        })
        for src, pr in d["paired_he1b"].items():
            mr = d["mismatched_all"][src]
            pair_rows.append({
                "line": ln,
                "source": src,
                "mis_est": mr.get("r_orth_cal_est"),
                "mis_ci975_lo": mr.get("r_orth_cal_ci975_lo"),
                "mis_ci975_hi": mr.get("r_orth_cal_ci975_hi"),
                "diff_est": pr.get("diff_est"),
                "diff_ci975_lo": pr.get("diff_ci975_lo"),
                "diff_ci975_hi": pr.get("diff_ci975_hi"),
                "passes_he1b_source": pr.get("passes_he1b_source"),
            })

    for ln, sblk in summary.get("subgroups_§7", summary.get("subgroups_s7", {})).items():
        for sg_name, sg in sblk.get("subgroups", {}).items():
            if sg.get("status") != "OK":
                sub_rows.append({"line": ln, "subgroup": sg_name, "status": sg.get("status"), "n_eval_targets": sg.get("n_eval_targets")})
                continue
            sf = sg["method_eval"]["fct"]["all_nontarget"]
            row = {
                "line": ln,
                "subgroup": sg_name,
                "status": "OK",
                "n_eval_targets": sg["n_eval_targets"],
                "fct_r_orth_cal_est": sf.get("r_orth_cal_est"),
                "fct_r_orth_cal_ci975_lo": sf.get("r_orth_cal_ci975_lo"),
                "fct_r_orth_cal_ci975_hi": sf.get("r_orth_cal_ci975_hi"),
                "fct_r_orth_cal_valid_frac": sf.get("r_orth_cal_valid_frac"),
                "fct_rho_ab_orth": sf.get("rho_ab_orth"),
                "perm100_p975": sg["perm100_all_nontarget"]["r_orth_cal"].get("p975"),
                "he1b_n_pass": sg["he1b_n_pass"],
                "he1b_n_total": sg["he1b_n_total"],
            }
            for src, pr in sg["paired_he1b"].items():
                row[f"diff_{src}_est"] = pr.get("diff_est")
                row[f"diff_{src}_lo"] = pr.get("diff_ci975_lo")
                row[f"diff_{src}_hi"] = pr.get("diff_ci975_hi")
            sub_rows.append(row)

    for fname, rows in (
        (f"audit_r14_lines_{tag}.csv", line_rows),
        (f"audit_r14_paired_he1b_{tag}.csv", pair_rows),
        (f"audit_r14_subgroups_{tag}.csv", sub_rows),
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
    parser = argparse.ArgumentParser(description="Independent Round 14 Audit Runner (`PREREG_G`, A-T3)")
    parser.add_argument(
        "--phase",
        choices=("synthetic_test", "dev_r9", "dev_r10", "dev_r11", "dev_r12", "dev_r13", "holdout_r14"),
        default="dev_r10",
    )
    parser.add_argument("--source-mode", choices=("r9_6", "r9_r10_8"), default="r9_6")
    parser.add_argument("--lines", nargs="+", default=None)
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--seed", type=int, default=DEFAULT_R14_SEED)
    parser.add_argument("--cov-max-cells", type=int, default=None)
    parser.add_argument("--min-targets", type=int, default=50)
    parser.add_argument("--min-genes", type=int, default=3000)
    parser.add_argument("--skip-subgroups", action="store_true")
    parser.add_argument("--no-hgnc", action="store_true")
    parser.add_argument("--out-dir", type=Path, default=AUDIT_R14_DIR / "dev_smoke")
    args = parser.parse_args()

    is_holdout = args.phase == "holdout_r14"
    if not is_holdout:
        if args.phase == "synthetic_test":
            raise RuntimeError("Use `test_r14_audit.py` or pass explicit synthetic AuditLineData for synthetic_test.")
        default_targets = list(DEV_ALLOWED_LINES_BY_PHASE[args.phase])[:1]
        targets = list(args.lines or default_targets)
        gate_meta = enforce_r14_gate(
            args.phase,
            target_lines=targets,
            extra_dirs=[resolve_r14_line_cache_dir(t, args.phase) for t in targets],
        )
        gate_info: Dict[str, Any] = {
            "confirmatory_lines": [],
            "k_conf": 0,
            "a549_myllia_role": "DEV_DRY_RUN",
            "secondary_eval_lines": targets,
        }
    else:
        gate_meta = enforce_r14_gate(
            "holdout_r14",
            target_lines=args.lines or (R14_TARGET_KEY,),
            status_dir=R14_STATUS_DIR,
            require_confirmatory_file=True,
        )
        gate_info = parse_confirmatory_lines_r14(R14_STATUS_DIR / "CONFIRMATORY_LINES_R14.txt")
        if gate_info.get("excluded_feasibility"):
            raise RuntimeError("PREREG_G §4.1: refused evaluating a549_myllia marked EXCLUDED_FEASIBILITY.")
        targets = list(args.lines) if args.lines else list(gate_info["confirmatory_lines"] or gate_info["secondary_eval_lines"])

    hgnc = None if args.no_hgnc else CanonicalHGNC.from_tsv()
    common_genes = load_r9_common_genes(hgnc=hgnc)

    source_pool_names = list(R9_LINE_NAMES)
    if args.source_mode == "r9_r10_8":
        source_pool_names.extend(R10_SOURCE_EXT_NAMES)

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
    subgroups_s7: Dict[str, Any] = {}

    for tgt in targets:
        if tgt in loaded_cache:
            tdata = loaded_cache[tgt]
        else:
            min_ov = int(args.min_genes) if is_holdout else 50
            tdata = load_line_dataset(
                resolve_r14_line_cache_dir(tgt, args.phase),
                representation="log1p",
                aligned_genes=common_genes,
                hgnc=hgnc,
                min_gene_overlap=min_ov,
            )
            cov_cache[tgt] = sample_covariance(deterministic_subsample(tdata.ctrl_cells, args.cov_max_cells, seed=42), shrink=0.0)
            gene_order_cache[tgt] = tdata.genes

        src_list = [loaded_cache[sn] for sn in source_pool_names if sn != tgt]
        thresholds_by_line[tgt] = check_r14_target_thresholds(
            tdata,
            src_list,
            only_reproducible=True,
            min_targets=args.min_targets,
            min_genes=args.min_genes,
        )
        per_line[tgt] = evaluate_target_line_r12(
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
        if not args.skip_subgroups:
            subgroups_s7[tgt] = evaluate_r14_subgroups(
                target_data=tdata,
                source_list=src_list,
                only_reproducible=True,
                cov_max_cells=args.cov_max_cells,
                n_boot=args.n_boot,
                n_perm=args.n_perm,
                seed=args.seed,
                precomputed_covs=cov_cache,
                precomputed_gene_orders=gene_order_cache,
            )

    if not is_holdout:
        adjudication = {
            "status": f"DEV_DRY_RUN_{args.phase.upper()}",
            "k_conf": 0,
            "confirmatory_lines": [],
            "claim1_he1_verdict": "DEV_ONLY_NO_ADJUDICATION",
            "claim2_he2_verdict": "DEV_ONLY_NO_ADJUDICATION",
            "same_line_replication_only": True,
            "limited_generalizability_notice": True,
            "mandatory_disclosure": "Pre-freeze dry-run on unblinded historical R9-R13 lines only.",
        }
    else:
        adjudication = adjudicate_round14(per_line, gate_info, min_targets=args.min_targets, min_genes=args.min_genes)

    summary = {
        "phase": args.phase,
        "source_mode": args.source_mode,
        "seed": args.seed,
        "n_boot": args.n_boot,
        "n_perm": args.n_perm,
        "gate_meta": gate_meta,
        "gate_info": gate_info,
        "thresholds": thresholds_by_line,
        "targets": targets,
        "adjudication": adjudication,
        "per_line": per_line,
        "subgroups_s7": subgroups_s7,
    }
    _write_r14_outputs(summary, args.out_dir, f"{args.phase}_{args.source_mode}")
    print(json.dumps({"status": adjudication["status"], "targets": targets, "out_dir": str(args.out_dir)}, indent=2))


if __name__ == "__main__":
    main()
