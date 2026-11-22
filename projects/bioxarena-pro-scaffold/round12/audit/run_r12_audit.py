#!/usr/bin/env python3
"""
Independent Round 12 Pre-Registration Audit Runner & Report Generator (A-T3).

Strictly follows `PREREG_E.md` (`v2`, commit `4c606f4`):
  - Strictly NEVER reads or imports `round12/src/`.
  - Pre-freeze (`--phase dev_r9 | dev_r10 | dev_r11`): executes ONLY on historical Round 9–11
    standardized caches (`seed=20261201`, `n_boot=4000`, `n_perm=100` defaults) and blocks
    any attempt to access `round12/shared/` or Round 12 datasets.
  - Holdout (`--phase holdout_r12`): blocked at runtime until `PREREG_E_COMMIT.txt`,
    `CODE_FREEZE_R12.txt`, and `CONFIRMATORY_LINES_R12.txt` exist and are non-empty.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np

AUDIT_R12_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R12_DIR.parents[1]
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
for _p in (str(AUDIT_R12_DIR), str(R11_AUDIT_DIR)):
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
from independent_r12_loader import (  # noqa: E402
    DEV_ALLOWED_LINES_BY_PHASE,
    R12_LINES_DIR,
    R12_STATUS_DIR,
    enforce_r12_gate,
    parse_confirmatory_lines_r12,
    resolve_line_cache_dir,
)
from independent_r12_metrics import (  # noqa: E402
    DEFAULT_R12_SEED,
    adjudicate_round12,
    evaluate_target_line_r12,
)


def _fmt(val: Any, digits: int = 4) -> str:
    try:
        fv = float(val)
        return f"{fv:+.{digits}f}" if np.isfinite(fv) else "nan"
    except (TypeError, ValueError):
        return str(val)


def _write_csvs(per_line: Dict[str, Dict[str, Any]], out_dir: Path, tag: str) -> None:
    line_rows: List[Dict[str, Any]] = []
    pair_rows: List[Dict[str, Any]] = []
    n5_rows: List[Dict[str, Any]] = []

    for ln, d in per_line.items():
        fct = d["method_eval"]["fct"]["all_nontarget"]
        e1_all = d["method_eval"]["e1_const_med"]["all_nontarget"]
        e1_top = d["method_eval"]["e1_const_med"]["resp_top50"]
        sa8_e1 = d["sa8_e1"]
        sa8_f = d["sa8_fct"]
        d_sa8 = d["sa8_e1_split_minus_all"]
        perm_all = d["perm100_all_nontarget"]["r_orth_cal"]
        line_rows.append({
            "line": ln,
            "n_genes_eval": d["n_genes_eval"],
            "n_eval_targets": d["n_eval_targets"],
            "n_sources": d["n_sources"],
            "fct_r_orth_cal_est": fct.get("r_orth_cal_est"),
            "fct_r_orth_cal_ci975_lo": fct.get("r_orth_cal_ci975_lo"),
            "fct_r_orth_cal_ci975_hi": fct.get("r_orth_cal_ci975_hi"),
            "fct_r_orth_cal_valid_frac": fct.get("r_orth_cal_valid_frac"),
            "fct_rho_ab_orth": fct.get("rho_ab_orth"),
            "perm100_all_median": perm_all.get("median"),
            "perm100_all_p975": perm_all.get("p975"),
            "he1a_pass": d["he1a_pass"],
            "he1b_pass": d["he1b_pass"],
            "he1b_n_pass": d["he1b_n_pass"],
            "he1b_n_total": d["he1b_n_total"],
            "he1c_pass": d["he1c_pass"],
            "he1_pass": d["he1_pass"],
            "sa8_e1_n_eval": sa8_e1.get("n_eval_perts"),
            "sa8_e1_split_est": sa8_e1.get("r_orth_cal_split_est"),
            "sa8_e1_split_ci975_lo": sa8_e1.get("r_orth_cal_split_ci975_lo"),
            "sa8_e1_split_ci975_hi": sa8_e1.get("r_orth_cal_split_ci975_hi"),
            "sa8_e1_split_minus_all_n": d_sa8.get("n_eval_perts"),
            "sa8_e1_split_minus_all_est": d_sa8.get("diff_est"),
            "sa8_e1_split_minus_all_ci975_lo": d_sa8.get("diff_ci975_lo"),
            "sa8_e1_split_minus_all_ci975_hi": d_sa8.get("diff_ci975_hi"),
            "he2a_pass": d["he2a_pass"],
            "he2b_pass": d["he2b_pass"],
            "he2_pass": d["he2_pass"],
            "sa8_fct_split_est": sa8_f.get("r_orth_cal_split_est"),
            "sa8_perm100_p975": d["sa8_perm100"]["r_orth_cal_split"]["p975"],
            "legacy_hd2_pass": d["legacy_hd2_pass"],
            "pair_fct_minus_e1_all_est": d["pair_fct_minus_e1_all"].get("diff_est"),
            "pair_fct_minus_e1_all_ci975_lo": d["pair_fct_minus_e1_all"].get("diff_ci975_lo"),
            "pair_fct_minus_e1_all_ci975_hi": d["pair_fct_minus_e1_all"].get("diff_ci975_hi"),
            "e1_all_est": e1_all.get("r_orth_cal_est"),
            "e1_top50_est": e1_top.get("r_orth_cal_est"),
            "hd_util_pass": d["hd_util"]["hd_util_pass_line"],
        })

        for src, pr in d["paired_he1b"].items():
            mr = d["mismatched_all"].get(src, {})
            pair_rows.append({
                "line": ln,
                "source_line": src,
                "n_eval": pr.get("n_eval_perts"),
                "fct_est": fct.get("r_orth_cal_est"),
                "mis_est": mr.get("r_orth_cal_est"),
                "mis_ci975_lo": mr.get("r_orth_cal_ci975_lo"),
                "mis_ci975_hi": mr.get("r_orth_cal_ci975_hi"),
                "paired_diff_est": pr.get("diff_est"),
                "paired_diff_ci975_lo": pr.get("diff_ci975_lo"),
                "paired_diff_ci975_hi": pr.get("diff_ci975_hi"),
                "paired_diff_valid_frac": pr.get("diff_valid_frac"),
                "passes_he1b_source": pr.get("passes_he1b_source"),
            })

        for mask_key in ("n5_all", "n5_split"):
            blk = d.get(mask_key) or {}
            for mn, mres in blk.get("methods", {}).items():
                n5_rows.append({
                    "line": ln,
                    "mask": blk.get("mode", mask_key),
                    "n_eval": mres.get("n_eval_perts"),
                    "method": mn,
                    "audit_classification": mres.get("audit_classification"),
                    "r_cone_sp_est": mres.get("r_cone_sp_est"),
                    "r_cone_sp_ci975_lo": mres.get("r_cone_sp_ci975_lo"),
                    "r_cone_sp_ci975_hi": mres.get("r_cone_sp_ci975_hi"),
                    "eta_est": mres.get("eta_est"),
                    "r2_in_n5": mres.get("r2_in_n5"),
                    "rho_ab_n5": mres.get("rho_ab_n5"),
                    "valid_frac": mres.get("r_cone_sp_valid_frac"),
                })

    for fname, rows in (
        (f"audit_r12_per_line_{tag}.csv", line_rows),
        (f"audit_r12_paired_diffs_{tag}.csv", pair_rows),
        (f"audit_r12_n5_{tag}.csv", n5_rows),
    ):
        if rows:
            with open(out_dir / fname, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)


def render_audit_markdown(summary: Dict[str, Any]) -> str:
    adj = summary["adjudication"]
    lines = [
        f"# Round 12 Independent Pre-Registration Audit Report (`{summary['phase']}`)",
        "",
        f"- **Source Mode:** `{summary['source_mode']}` | **Seed:** `{summary['seed']}` | "
        f"**Bootstraps:** `{summary['n_boot']}` (`97.5%` two-sided Bonferroni CIs) | **Perms:** `{summary['n_perm']}`",
        f"- **Adjudication Status:** `{adj['status']}` (`K_conf = {adj.get('k_conf', 0)}`, "
        f"confirmatory lines = `{adj.get('confirmatory_lines', [])}`)",
        f"- **Claim 1 (`H-E1`):** `{adj.get('claim1_he1_verdict', 'DEV')}` | "
        f"**Claim 2 (`H-E2`):** `{adj.get('claim2_he2_verdict', 'DEV')}`",
        "",
        "## Mandatory Pre-Registered Disclosures (`PREREG_E §0.1, §0.2, §7, §8`)",
        "- **Post-Hoc Provenance (`E6`):** Both `H-E1` (paired `FCT - Mis_s`) and `H-E2` (`SA-8`) hypotheses — and the interpretation that prior T-cell unpaired `H-D2` failures reflected low unpaired power — were formulated after inspecting Rounds 9–11 data.",
        "- **Claim 2 Empirical & Metric Scope (`E6`):** `SA-8` is strictly an empirical finding on residual metric `r_orth_cal` (`resp_top50_split`), NOT mechanistic proof (Round 11 Check 4 noiseless simulation predicted the wrong sign on `8/12` lines under `all_nontarget`), and does NOT extrapolate to Nicol et al. (2026) DEG metrics.",
        "- **Claim 1 Novelty Scope (`§8`):** Must cite Kawsar et al. 2026 (`CanPert/CASTER`, `doi:10.64898/2026.11.15.766815`) for prior per-target paired matched-vs-mismatched control comparison; Claim 1 tests prospectively whether baseline-only `FCT` outperforms all `6` mismatched-line baselines on residual `r_orth_cal` in an unseen lineage and modality.",
        f"- **Scope & Modality:** {adj.get('mandatory_disclosure', '')}",
        f"- **CRISPRi -> CRISPRko Modality Boundary:** {adj.get('modality_interpretation', '')}",
        "",
        "## Confirmatory Endpoints (`H-E1` & `H-E2`, `r_orth_cal`)",
        "| Line (`N`) | `FCT` `[97.5% CI]` (`vf`, `rho`) | `H-E1a` | `H-E1b` Paired Pass | `H-E1c` (`> Perm97.5`) | `H-E1` | `E1(split)` `[97.5% CI]` | `E1(split)-E1(all)` `[97.5% CI]` (`N_c`) | `H-E2` |",
        "|---|---|:---:|:---:|:---:|:---:|---|---|:---:|",
    ]
    for ln, d in summary["per_line"].items():
        f = d["method_eval"]["fct"]["all_nontarget"]
        p975 = d["perm100_all_nontarget"]["r_orth_cal"]["p975"]
        se1 = d["sa8_e1"]
        ds = d["sa8_e1_split_minus_all"]
        lines.append(
            f"| `{ln}` (`{d['n_eval_targets']}`) | "
            f"`{_fmt(f.get('r_orth_cal_est'))}` `[{_fmt(f.get('r_orth_cal_ci975_lo'))}, {_fmt(f.get('r_orth_cal_ci975_hi'))}]` "
            f"(`{f.get('r_orth_cal_valid_frac', 0):.3f}`, `{_fmt(f.get('rho_ab_orth'), 3)}`) | "
            f"`{d['he1a_pass']}` | `{d['he1b_n_pass']}/{d['he1b_n_total']}` (`{d['he1b_pass']}`) | "
            f"`{d['he1c_pass']}` (`{_fmt(p975)}`) | **`{d['he1_pass']}`** | "
            f"`{_fmt(se1.get('r_orth_cal_split_est'))}` `[{_fmt(se1.get('r_orth_cal_split_ci975_lo'))}, {_fmt(se1.get('r_orth_cal_split_ci975_hi'))}]` | "
            f"`{_fmt(ds.get('diff_est'))}` `[{_fmt(ds.get('diff_ci975_lo'))}, {_fmt(ds.get('diff_ci975_hi'))}]` (`{ds.get('n_eval_perts')}`) | "
            f"**`{d['he2_pass']}`** |"
        )
    lines.extend(["", "## `H-E1b` Source-by-Source Paired `FCT - Mis_s` (`all_nontarget`, `r_orth_cal`)"])
    lines.append("| Target | Source `s'` | `Mis_s` `[97.5% CI]` | Paired `FCT - Mis_s` `[97.5% CI]` | Pass |")
    lines.append("|---|---|---|---|---|")
    for ln, d in summary["per_line"].items():
        for src, pr in d["paired_he1b"].items():
            mr = d["mismatched_all"][src]
            lines.append(
                f"| `{ln}` | `{src}` | `{_fmt(mr.get('r_orth_cal_est'))}` `[{_fmt(mr.get('r_orth_cal_ci975_lo'))}, {_fmt(mr.get('r_orth_cal_ci975_hi'))}]` | "
                f"`{_fmt(pr.get('diff_est'))}` `[{_fmt(pr.get('diff_ci975_lo'))}, {_fmt(pr.get('diff_ci975_hi'))}]` | `{pr.get('passes_he1b_source')}` |"
            )

    lines.extend([
        "",
        "## §5 Non-Confirmatory Descriptive Summary (`只报告、不判定`)",
        "*Explicit notice (`PREREG_E §5.3`): Within `N_5`, `Transfer_Mean`, `Transfer_Mean_GlobalScale` (`eta=0`), and `Source_Grand_Mean` (`eta=0`) equal zero by exact algebraic construction, and `E1_ConstMed`, `E2_ConstP10/P90`, `E6_LogFCBar`, `E7_FCBarMinus1` are suppressed by cone construction — none constitute empirical validation evidence.*",
    ])
    for ln, d in summary["per_line"].items():
        n5a = (d.get("n5_all") or {}).get("methods", {})
        n5_fct = n5a.get("fct", {})
        n5_gs = n5a.get("transfer_mean_global_scale", {})
        n5_e1 = n5a.get("e1_const_med", {})
        lines.append(
            f"- **`{ln}`**: Legacy unpaired `H-D2`=`{d['legacy_hd2_pass']}`; "
            f"`FCT - E1 (all)`=`{_fmt(d['pair_fct_minus_e1_all'].get('diff_est'))}` "
            f"`[{_fmt(d['pair_fct_minus_e1_all'].get('diff_ci975_lo'))}, {_fmt(d['pair_fct_minus_e1_all'].get('diff_ci975_hi'))}]`; "
            f"`N_5(all) FCT`=`{_fmt(n5_fct.get('r_cone_sp_est'))}` (`eta={n5_fct.get('eta_est', float('nan')):.3f}`), "
            f"`GS`=`{_fmt(n5_gs.get('r_cone_sp_est'))}` (`eta={n5_gs.get('eta_est', float('nan')):.1e}`), "
            f"`E1`=`{_fmt(n5_e1.get('r_cone_sp_est'))}` (`eta={n5_e1.get('eta_est', float('nan')):.4f}`); "
            f"`H-D_util` pass=`{d['hd_util']['hd_util_pass_line']}` (`R^2(FCT orth/all)` approx via unprojected `r_cal^2` reported in JSON)."
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Independent Round 12 Pre-Registration Audit Runner (A-T3)")
    parser.add_argument(
        "--phase",
        choices=("dev_r9", "dev_r10", "dev_r11", "holdout_r12"),
        default="dev_r10",
    )
    parser.add_argument("--source-mode", choices=("r9_6", "r9_r10_8"), default="r9_6")
    parser.add_argument("--lines", nargs="+", default=None, help="Target lines to evaluate")
    parser.add_argument("--n-boot", type=int, default=4000)
    parser.add_argument("--n-perm", type=int, default=100)
    parser.add_argument("--seed", type=int, default=DEFAULT_R12_SEED)
    parser.add_argument("--cov-max-cells", type=int, default=None)
    parser.add_argument("--no-hgnc", action="store_true")
    parser.add_argument("--skip-n5", action="store_true")
    parser.add_argument("--skip-n5-split", action="store_true")
    parser.add_argument("--skip-r19", action="store_true")
    parser.add_argument("--skip-s7", action="store_true")
    parser.add_argument("--out-dir", type=Path, default=AUDIT_R12_DIR / "dev_smoke")
    args = parser.parse_args()

    is_dev = args.phase.startswith("dev_")
    if is_dev:
        default_targets = list(DEV_ALLOWED_LINES_BY_PHASE[args.phase])[:1]
        targets = list(args.lines or default_targets)
        gate_meta = enforce_r12_gate(
            args.phase,
            target_lines=targets,
            extra_dirs=[resolve_line_cache_dir(t, args.phase) for t in targets],
        )
        gate_info: Dict[str, Any] = {
            "confirmatory_lines": [],
            "k_conf": 0,
            "a549_role": "DEV_DRY_RUN",
            "hek293_role": "DEV_DRY_RUN",
            "secondary_eval_lines": targets,
        }
    else:
        gate_meta = enforce_r12_gate(
            "holdout_r12",
            target_lines=args.lines or (),
            status_dir=R12_STATUS_DIR,
            require_confirmatory_file=True,
        )
        gate_info = parse_confirmatory_lines_r12(R12_STATUS_DIR / "CONFIRMATORY_LINES_R12.txt")
        if args.lines is not None:
            # Never permit evaluating HEK293 if marked EXCLUDED_FEASIBILITY_Q3
            if gate_info.get("hek293_excluded_feasibility_q3") and "hek293_xucao" in args.lines:
                raise RuntimeError("PREREG_E §3.2 / E3: refused evaluating hek293_xucao (EXCLUDED_FEASIBILITY_Q3).")
            targets = list(args.lines)
        else:
            targets = list(gate_info["confirmatory_lines"]) + [
                t for t in gate_info["secondary_eval_lines"] if t not in gate_info["confirmatory_lines"]
            ]
            if (R12_LINES_DIR / "a549_ko_cas9pos" / "pseudobulk_delta.npz").exists() and "a549_ko_cas9pos" not in targets:
                targets.append("a549_ko_cas9pos")

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
    for tgt in targets:
        if tgt in loaded_cache:
            tdata = loaded_cache[tgt]
        else:
            # Secondary/dev lines or HEK293 SECONDARY_ONLY may fall below 3,000 overlap
            min_ov = 3000 if (not is_dev and tgt in gate_info["confirmatory_lines"]) else 50
            tdata = load_line_dataset(
                resolve_line_cache_dir(tgt, args.phase),
                representation="log1p",
                aligned_genes=common_genes,
                hgnc=hgnc,
                min_gene_overlap=min_ov,
            )
            cov_cache[tgt] = sample_covariance(deterministic_subsample(tdata.ctrl_cells, args.cov_max_cells, seed=42), shrink=0.0)
            gene_order_cache[tgt] = tdata.genes

        # Enforce zero target-self leakage in LOO dev checks
        src_list = [loaded_cache[sn] for sn in source_pool_names if sn != tgt]
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
            include_n5=not args.skip_n5,
            include_n5_split=not args.skip_n5_split,
            include_r19=not args.skip_r19,
            include_s7=not args.skip_s7,
            precomputed_covs=cov_cache,
            precomputed_gene_orders=gene_order_cache,
        )

    if is_dev:
        adjudication = {
            "status": f"DEV_DRY_RUN_{args.phase.upper()}",
            "k_conf": 0,
            "confirmatory_lines": [],
            "claim1_he1_verdict": "DEV_ONLY_NO_ADJUDICATION",
            "claim2_he2_verdict": "DEV_ONLY_NO_ADJUDICATION",
            "limited_generalizability_notice": False,
            "modality_interpretation": "DEV_R9_TO_R11_VALIDATION",
            "mandatory_disclosure": "Pre-freeze verification strictly on unblinded historical R9-R11 lines.",
        }
    else:
        adjudication = adjudicate_round12(per_line, gate_info)

    summary = {
        "phase": args.phase,
        "source_mode": args.source_mode,
        "seed": args.seed,
        "n_boot": args.n_boot,
        "n_perm": args.n_perm,
        "gate_meta": gate_meta,
        "gate_info": gate_info,
        "targets": targets,
        "adjudication": adjudication,
        "per_line": per_line,
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{args.phase}_{args.source_mode}"
    _write_csvs(per_line, args.out_dir, tag)
    (args.out_dir / f"audit_r12_summary_{tag}.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (args.out_dir / f"AUDIT_R12_REPORT_{tag}.md").write_text(
        render_audit_markdown(summary),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": f"R12_AUDIT_{args.phase.upper()}_COMPLETE",
                "source_mode": args.source_mode,
                "targets": targets,
                "adjudication_status": adjudication["status"],
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
