#!/usr/bin/env python3
"""
Round 15 Prospective Group E Runner (`PREREG_H` v2, Arc VCC 2026 Multi-Line Holdout).

Reuses frozen Round 10 / 11 / 12 / 14 evaluation primitives strictly unmodified:
  - Uses `DEFAULT_R15_SEED = 20261231`, `n_boot = 4,000`, `n_perm = 100`, `r_orth_cal`,
    `all_nontarget`, `resp_top50_split`, `97.5%` two-sided Bonferroni CIs, IUT over
    `K_conf in 0..5` lines (`ALLOWED_R15_LINES`; hard `jurkat` exclusion).
  - Enforces hard `PREREG_H_COMMIT.txt` + `CODE_FREEZE_R15.txt` blinding gate on
    `--phase holdout_r15` BEFORE loading any target expression dataset.
  - Emits strictly descriptive `PREREG_H §9` secondary analyses:
      * `SA-H1`: Leave-`vcc2025_h1`-out sensitivity across remaining 5 `r9_6` sources
      * `SA-Panel`: `129` shared val-test target overlap vs remainder (when panel file available)
  - Supports `PREREG_H §8.4` release-wave evaluation (`--wave-id`, `--prior-wave-jsons`)
    emitting both per-wave adjudication and cumulative cross-wave IUT summary.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

ROUND15_SRC_DIR = Path(__file__).resolve().parent
ROUND15_ROOT = ROUND15_SRC_DIR.parent
REPO_ROOT = ROUND15_ROOT.parent
R10_SRC_DIR = REPO_ROOT / "round10" / "src"
R11_SRC_DIR = REPO_ROOT / "round11" / "src"
R12_SRC_DIR = REPO_ROOT / "round12" / "src"
R14_SRC_DIR = REPO_ROOT / "round14" / "src"
for _p in (str(ROUND15_SRC_DIR), str(R10_SRC_DIR), str(R11_SRC_DIR), str(R12_SRC_DIR), str(R14_SRC_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from adjudicate_r15 import (  # noqa: E402
    adjudicate_r15_results,
    aggregate_r15_wave_adjudications,
    render_r15_markdown_report,
    render_r15_overall_markdown_report,
)
from data_loader_r15 import (  # noqa: E402
    DEFAULT_HGNC_FILE,
    DEFAULT_R15_SEED,
    HGNCBridge,
    LineDatasetR10,
    R9_COMMON_GENES_FILE,
    R9_LINES_DIR,
    R10_LINES_DIR,
    R11_LINES_DIR,
    R15_CONFIRMATORY_LINES_FILE,
    R15_LINES_DIR,
    R15_RESULTS_DIR,
    SA_H1_EXCLUDED_SOURCE,
    check_prereg_h_gate,
    list_standardized_lines,
    load_r9_common_genes,
    load_source_cohort_r15,
    load_standardized_line_r10,
    load_vcc2026_panel_overlap_targets,
    parse_confirmatory_lines_r15,
    validate_r15_candidate_line,
)
from run_eval_r12 import _dev_min_overlap, evaluate_target_against_sources_r12  # noqa: E402
from run_eval_r14 import compute_overlap_subgroup_descriptive_r14  # noqa: E402


def compute_sa_h1_leave_vcc2025_h1_out_r15(
    target_ds: LineDatasetR10,
    source_datasets: Dict[str, LineDatasetR10],
    only_reproducible: bool = True,
    s5_top75_genes: bool = False,
    include_r9_anchors: bool = True,
    include_exploratory_e: bool = True,
    min_gene_overlap: int = 3000,
    n_boot: int = 4000,
    seed: int = DEFAULT_R15_SEED,
    n_perm: int = 100,
) -> Dict[str, Any]:
    """
    `PREREG_H §9.1 SA-H1` (strictly descriptive): recompute `FCT r_orth_cal [97.5% CI]`,
    `5` paired `H-E1b` contrasts `FCT - Mis_s'`, `Perm100`, and `H-E2` (`SA-8`) after
    excluding `vcc2025_h1` from `source_datasets`.
    """
    no_h1_sources = {k: v for k, v in source_datasets.items() if k != SA_H1_EXCLUDED_SOURCE}
    if len(no_h1_sources) == 0:
        raise ValueError("SA-H1 requires >= 1 source remaining after excluding vcc2025_h1.")

    l_rows, b_rows, ex = evaluate_target_against_sources_r12(
        target_ds=target_ds,
        source_datasets=no_h1_sources,
        only_reproducible=only_reproducible,
        s5_top75_genes=s5_top75_genes,
        include_r9_anchors=include_r9_anchors,
        include_exploratory_e=include_exploratory_e,
        min_gene_overlap=min_gene_overlap,
        n_boot=n_boot,
        seed=seed,
        n_perm=n_perm,
    )
    by_m = {str(r["method"]): r for r in l_rows}
    fct_r = by_m["FCT"]
    e1_r = by_m["E1_ConstMed"]

    paired: Dict[str, Dict[str, Any]] = {}
    for br in b_rows:
        comp = str(br.get("comparison", ""))
        if comp.startswith("H-E1b_FCT_minus_Mis_") and comp.endswith("_all_nontarget"):
            src = comp[len("H-E1b_FCT_minus_Mis_") : -len("_all_nontarget")]
            paired[src] = {
                "estimate": float(br["estimate"]),
                "ci975_lo": float(br["ci975_lo"]),
                "ci975_hi": float(br["ci975_hi"]),
                "valid_frac": float(br["valid_frac"]),
                "n_items": int(br["n_items"]),
            }

    h2b = next((r for r in b_rows if r.get("comparison") == "H-E2b_E1_split_minus_all"), None)
    return {
        "target_line": target_ds.name,
        "excluded_source": SA_H1_EXCLUDED_SOURCE,
        "remaining_sources": sorted(no_h1_sources.keys()),
        "n_eval_targets": int(fct_r.get("n_eval_targets", 0)),
        "fct": {
            "r_orth_cal": float(fct_r["r_orth_cal"]),
            "ci975_lo": float(fct_r["r_orth_cal_ci975_lo"]),
            "ci975_hi": float(fct_r["r_orth_cal_ci975_hi"]),
            "valid_frac": float(fct_r.get("r_orth_cal_valid_frac", np.nan)),
            "rho_ab_orth": float(fct_r["rho_half_ab_orth"]),
        },
        "e1_split": {
            "r_orth_cal_split": float(e1_r.get("resp_top50_split_r_orth_cal", np.nan)),
            "ci975_lo": float(e1_r.get("resp_top50_split_r_orth_cal_ci975_lo", np.nan)),
            "ci975_hi": float(e1_r.get("resp_top50_split_r_orth_cal_ci975_hi", np.nan)),
            "h_e2b_diff_est": float(h2b["estimate"]) if h2b is not None else float("nan"),
            "h_e2b_diff_ci975_lo": float(h2b["ci975_lo"]) if h2b is not None else float("nan"),
            "h_e2b_diff_ci975_hi": float(h2b["ci975_hi"]) if h2b is not None else float("nan"),
        },
        "paired_fct_minus_mismatched": paired,
        "perm100": ex.get("perm100", {}),
        "line_rows": l_rows,
        "boot_rows": b_rows,
    }


def compute_sa_panel_overlap_descriptive_r15(
    target_ds: LineDatasetR10,
    source_datasets: Dict[str, LineDatasetR10],
    panel_overlap_targets: Set[str],
    only_reproducible: bool = True,
    s5_top75_genes: bool = False,
    include_r9_anchors: bool = True,
    include_exploratory_e: bool = True,
    n_boot: int = 4000,
    seed: int = DEFAULT_R15_SEED,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    """
    `PREREG_H §9.2 SA-Panel` (strictly descriptive): partition `eval_targets` into
    `in_val_test_overlap_panel` (inside the 129 shared val-test set) vs
    `panel_unique_remainder`, reporting subgroup counts, `FCT r_orth_cal [97.5% CI]`,
    and `H-E1b` paired contrasts cleanly using the frozen R14 subgroup engine.
    """
    raw_d, raw_rows = compute_overlap_subgroup_descriptive_r14(
        target_ds=target_ds,
        source_datasets=source_datasets,
        r12_overlap_targets=panel_overlap_targets,
        only_reproducible=only_reproducible,
        s5_top75_genes=s5_top75_genes,
        include_r9_anchors=include_r9_anchors,
        include_exploratory_e=include_exploratory_e,
        n_boot=n_boot,
        seed=seed,
    )
    key_map = {
        "in_r12_a549_overlap": "in_val_test_overlap_panel",
        "novel_not_in_r12_a549": "panel_unique_remainder",
    }
    out_d = {key_map.get(k, k): v for k, v in raw_d.items()}
    out_rows = [{**r, "subgroup": key_map.get(str(r.get("subgroup", "")), str(r.get("subgroup", "")))} for r in raw_rows]
    return out_d, out_rows


def evaluate_r15_cohort(
    target_ds_map: Dict[str, LineDatasetR10],
    source_datasets: Dict[str, LineDatasetR10],
    confirmatory_lines: Sequence[str],
    excluded_lines: Sequence[str] = (),
    panel_overlap_targets: Optional[Set[str]] = None,
    panel_overlap_meta: Optional[Dict[str, Any]] = None,
    phase: str = "holdout_r15",
    source_mode: str = "r9_6",
    representation: str = "log1p",
    prereg_commit_hash: str = "",
    code_freeze_hash: str = "",
    wave_id: str = "wave1",
    sunset_reached: bool = False,
    only_reproducible: bool = True,
    s5_top75_genes: bool = False,
    include_r9_anchors: bool = True,
    include_exploratory_e: bool = True,
    min_gene_overlap: Optional[int] = None,
    n_boot: int = 4000,
    n_perm: int = 100,
    seed: int = DEFAULT_R15_SEED,
    run_sa_h1: bool = True,
) -> Dict[str, Any]:
    """
    Evaluate a set of standardized target lines cleanly under `PREREG_H` v2.
    """
    all_l_rows: List[Dict[str, Any]] = []
    all_b_rows: List[Dict[str, Any]] = []
    perm_map: Dict[str, Any] = {}
    n5_map: Dict[str, Any] = {}
    n5_sp_map: Dict[str, Any] = {}
    n5_perm_map: Dict[str, Any] = {}
    sa_h1_map: Dict[str, Any] = {}
    sa_h1_line_rows: List[Dict[str, Any]] = []
    sa_h1_boot_rows: List[Dict[str, Any]] = []
    sa_panel_map: Dict[str, Any] = {}
    sa_panel_rows: List[Dict[str, Any]] = []

    for ln, t_ds in target_ds_map.items():
        if phase == "holdout_r15":
            validate_r15_candidate_line(ln)
        m_ov = (
            int(min_gene_overlap)
            if min_gene_overlap is not None
            else (3000 if phase in ("holdout_r15", "regression_legacy") else _dev_min_overlap(phase, ln))
        )
        # Drop target line if accidentally present in source_datasets (e.g. dev_r9 leave-self-out parity)
        src_for_ln = {k: v for k, v in source_datasets.items() if k != ln}

        l_rows, b_rows, extra = evaluate_target_against_sources_r12(
            target_ds=t_ds,
            source_datasets=src_for_ln,
            only_reproducible=only_reproducible,
            s5_top75_genes=s5_top75_genes,
            include_r9_anchors=include_r9_anchors,
            include_exploratory_e=include_exploratory_e,
            min_gene_overlap=m_ov,
            n_boot=n_boot,
            seed=seed,
            n_perm=n_perm,
        )
        all_l_rows.extend(l_rows)
        all_b_rows.extend(b_rows)
        perm_map[ln] = extra.get("perm100", {})
        n5_map[ln] = extra.get("n5_cone_metrics", {})
        n5_sp_map[ln] = extra.get("n5_cone_split_metrics", {})
        n5_perm_map[ln] = extra.get("n5_perm100", {})

        if run_sa_h1 and SA_H1_EXCLUDED_SOURCE in src_for_ln and len(src_for_ln) > 1:
            sa1 = compute_sa_h1_leave_vcc2025_h1_out_r15(
                target_ds=t_ds,
                source_datasets=src_for_ln,
                only_reproducible=only_reproducible,
                s5_top75_genes=s5_top75_genes,
                include_r9_anchors=include_r9_anchors,
                include_exploratory_e=include_exploratory_e,
                min_gene_overlap=m_ov,
                n_boot=n_boot,
                seed=seed,
                n_perm=n_perm,
            )
            sa_h1_line_rows.extend([{**r, "sa_analysis": "SA-H1_no_vcc2025_h1"} for r in sa1.pop("line_rows", [])])
            sa_h1_boot_rows.extend([{**r, "sa_analysis": "SA-H1_no_vcc2025_h1"} for r in sa1.pop("boot_rows", [])])
            sa_h1_map[ln] = sa1

        if panel_overlap_targets is not None:
            pan_d, pan_rows = compute_sa_panel_overlap_descriptive_r15(
                target_ds=t_ds,
                source_datasets=src_for_ln,
                panel_overlap_targets=panel_overlap_targets,
                only_reproducible=only_reproducible,
                s5_top75_genes=s5_top75_genes,
                include_r9_anchors=include_r9_anchors,
                include_exploratory_e=include_exploratory_e,
                n_boot=n_boot,
                seed=seed,
            )
            sa_panel_map[ln] = pan_d
            sa_panel_rows.extend(pan_rows)

    df_lines = pd.DataFrame(all_l_rows)
    df_boot = pd.DataFrame(all_b_rows)
    summary: Dict[str, Any] = {
        "protocol": "PREREG_H_v2",
        "phase": phase,
        "wave_id": wave_id,
        "source_mode": source_mode,
        "representation": representation,
        "source_lines": sorted(source_datasets.keys()),
        "target_lines": sorted(target_ds_map.keys()),
        "prereg_commit_hash": prereg_commit_hash,
        "code_freeze_hash": code_freeze_hash,
        "bootstrap_seed": int(seed),
        "n_boot": int(n_boot),
        "n_perm": int(n_perm),
        "fct_perm100": perm_map,
        "n5_cone_metrics": n5_map,
        "n5_cone_split_metrics": n5_sp_map,
        "n5_perm100": n5_perm_map,
        "sa_h1": sa_h1_map,
        "sa_panel_meta": panel_overlap_meta or {"available": panel_overlap_targets is not None, "status": "OK" if panel_overlap_targets is not None else "PANEL_OVERLAP_LIST_NOT_AVAILABLE"},
        "sa_panel": sa_panel_map,
    }

    adj = adjudicate_r15_results(
        df_lines=df_lines,
        df_boot=df_boot,
        summary=summary,
        confirmatory_lines=confirmatory_lines if phase in ("holdout_r15", "regression_legacy") else [],
        excluded_lines=excluded_lines,
        wave_id=wave_id,
        sunset_reached=sunset_reached,
        allow_legacy_lines=(phase == "regression_legacy"),
    )
    return {
        "df_lines": df_lines,
        "df_boot": df_boot,
        "df_sa_h1_lines": pd.DataFrame(sa_h1_line_rows),
        "df_sa_h1_boot": pd.DataFrame(sa_h1_boot_rows),
        "df_sa_panel": pd.DataFrame(sa_panel_rows),
        "summary": summary,
        "adj": adj,
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Round 15 (`PREREG_H` v2) prospective multi-line evaluator.")
    p.add_argument("--phase", choices=["dev_r9", "dev_r10", "dev_r11", "holdout_r15", "regression_legacy"], required=True)
    p.add_argument("--source-mode", choices=["r9_6", "r9_r10_8"], default="r9_6")
    p.add_argument("--repr", choices=["log1p", "zscore", "libnorm"], default="log1p")
    p.add_argument("--all-targets", action="store_true", default=False)
    p.add_argument("--s5-top75", action="store_true", default=False)
    p.add_argument("--no-r9-anchors", action="store_true", default=False)
    p.add_argument("--no-exploratory-e", action="store_true", default=False)
    p.add_argument("--no-canon-hgnc", action="store_true", default=False)
    p.add_argument("--lines", type=str, default="", help="Comma-separated target lines to evaluate in this wave.")
    p.add_argument("--confirmatory-lines", type=str, default="", help="Optional confirmatory lines override or NONE.")
    p.add_argument("--conf-file", type=Path, default=R15_CONFIRMATORY_LINES_FILE)
    p.add_argument("--excluded-lines", type=str, default="", help="Optional comma-separated lines excluded pre-metric in this wave.")
    p.add_argument("--lines-root", type=Path, default=R15_LINES_DIR)
    p.add_argument("--min-gene-overlap", type=int, default=None, help="Optional min gene overlap override (default 3000 on holdout/regression).")
    p.add_argument("--gate-dir", type=Path, default=None)
    p.add_argument("--panel-overlap-file", type=Path, default=None, help="Optional 129 shared val-test targets list for SA-Panel.")
    p.add_argument("--skip-sa-h1", action="store_true", default=False)
    p.add_argument("--wave-id", type=str, default="wave1")
    p.add_argument("--prior-wave-jsons", type=Path, nargs="*", default=[])
    p.add_argument("--sunset-reached", action="store_true", default=False)
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--n-perm", type=int, default=100)
    p.add_argument("--seed", type=int, default=DEFAULT_R15_SEED)
    p.add_argument("--tag-suffix", type=str, default="")
    p.add_argument("--out-dir", type=Path, default=R15_RESULTS_DIR)
    return p.parse_args()


def main() -> None:
    t0 = time.time()
    args = parse_args()
    only_rep = not args.all_targets
    include_anchors = not args.no_r9_anchors
    hgnc = None if args.no_canon_hgnc else HGNCBridge.from_tsv(DEFAULT_HGNC_FILE)
    common_r9 = load_r9_common_genes(hgnc_bridge=hgnc, common_genes_file=R9_COMMON_GENES_FILE)

    prereg_hash = freeze_hash = ""
    if args.phase in ("dev_r9", "dev_r10", "dev_r11"):
        src_ds_base = load_source_cohort_r15("r9_6" if args.phase == "dev_r9" else args.source_mode, args.repr, common_r9, hgnc)
        lines_root = {"dev_r9": R9_LINES_DIR, "dev_r10": R10_LINES_DIR, "dev_r11": R11_LINES_DIR}[args.phase]
        avail = list_standardized_lines(lines_root)
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail
        target_ds_map = {
            ln: load_standardized_line_r10(
                lines_root / ln,
                representation=args.repr,
                common_genes=common_r9,
                hgnc_bridge=hgnc,
                min_gene_overlap=_dev_min_overlap(args.phase, ln),
            )
            for ln in eval_lines
        }
        conf_lines: List[str] = []
    elif args.phase == "regression_legacy":
        # Non-confirmatory legacy regression mode pointing explicitly at prior-round directories
        src_ds_base = load_source_cohort_r15(args.source_mode, args.repr, common_r9, hgnc)
        lines_root = Path(args.lines_root)
        avail = list_standardized_lines(lines_root)
        eval_lines = [x.strip() for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail
        raw_conf = [x.strip() for x in args.confirmatory_lines.split(",") if x.strip()]
        conf_lines = [] if (not raw_conf or (len(raw_conf) == 1 and raw_conf[0].upper() in {"NONE", "EMPTY", "K0"})) else raw_conf
        m_ov = int(args.min_gene_overlap) if args.min_gene_overlap is not None else 3000
        target_ds_map = {
            ln: load_standardized_line_r10(
                lines_root / ln,
                representation=args.repr,
                common_genes=common_r9,
                hgnc_bridge=hgnc,
                min_gene_overlap=m_ov,
            )
            for ln in eval_lines
        }
    else:
        # Hard PREREG_H + CODE_FREEZE_R15 blinding lock BEFORE loading ANY held-out data
        prereg_hash, freeze_hash = check_prereg_h_gate(args.gate_dir)
        src_ds_base = load_source_cohort_r15(args.source_mode, args.repr, common_r9, hgnc)
        lines_root = Path(args.lines_root)
        avail = [validate_r15_candidate_line(x) for x in list_standardized_lines(lines_root)]
        eval_lines = [validate_r15_candidate_line(x) for x in args.lines.split(",") if x.strip()] if args.lines.strip() else avail
        if args.confirmatory_lines.strip():
            conf_lines = parse_confirmatory_lines_r15(raw_lines=[args.confirmatory_lines])
        else:
            conf_lines = parse_confirmatory_lines_r15(conf_file=args.conf_file)
        m_ov = int(args.min_gene_overlap) if args.min_gene_overlap is not None else 3000
        target_ds_map = {
            ln: load_standardized_line_r10(
                lines_root / ln,
                representation=args.repr,
                common_genes=common_r9,
                hgnc_bridge=hgnc,
                min_gene_overlap=m_ov,
            )
            for ln in eval_lines
        }

    panel_targets, panel_meta = load_vcc2026_panel_overlap_targets(args.panel_overlap_file, hgnc_bridge=hgnc)
    excluded_lines = [x.strip() for x in args.excluded_lines.split(",") if x.strip()]

    res = evaluate_r15_cohort(
        target_ds_map=target_ds_map,
        source_datasets=src_ds_base,
        confirmatory_lines=conf_lines,
        excluded_lines=excluded_lines,
        panel_overlap_targets=panel_targets,
        panel_overlap_meta=panel_meta,
        phase=args.phase,
        source_mode=args.source_mode,
        representation=args.repr,
        prereg_commit_hash=prereg_hash,
        code_freeze_hash=freeze_hash,
        wave_id=args.wave_id,
        sunset_reached=bool(args.sunset_reached),
        only_reproducible=only_rep,
        s5_top75_genes=args.s5_top75,
        include_r9_anchors=include_anchors,
        include_exploratory_e=not args.no_exploratory_e,
        n_boot=args.n_boot,
        n_perm=args.n_perm,
        seed=args.seed,
        run_sa_h1=not args.skip_sa_h1,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_tag = f"{args.phase}_{args.wave_id}_{args.source_mode}_{args.repr}{('_' + args.tag_suffix) if args.tag_suffix else ''}"

    res["df_lines"].to_csv(out_dir / f"lines_{run_tag}.csv", index=False)
    res["df_boot"].to_csv(out_dir / f"boot_{run_tag}.csv", index=False)
    if not res["df_sa_h1_lines"].empty:
        res["df_sa_h1_lines"].to_csv(out_dir / f"sa_h1_lines_{run_tag}.csv", index=False)
        res["df_sa_h1_boot"].to_csv(out_dir / f"sa_h1_boot_{run_tag}.csv", index=False)
    if not res["df_sa_panel"].empty:
        res["df_sa_panel"].to_csv(out_dir / f"sa_panel_{run_tag}.csv", index=False)

    res["summary"]["elapsed_sec"] = round(time.time() - t0, 2)
    (out_dir / f"summary_{run_tag}.json").write_text(json.dumps(res["summary"], indent=2) + "\n", encoding="utf-8")
    (out_dir / f"adjudicate_{run_tag}.json").write_text(json.dumps(res["adj"], indent=2) + "\n", encoding="utf-8")
    md_text = render_r15_markdown_report(res["df_lines"], res["df_boot"], res["summary"], res["adj"], run_tag)
    (out_dir / f"report_{run_tag}.md").write_text(md_text, encoding="utf-8")

    # Optional cross-wave overall IUT aggregation (`PREREG_H §8.4`)
    if args.prior_wave_jsons:
        wave_objs = [json.loads(Path(p).read_text(encoding="utf-8")) for p in args.prior_wave_jsons]
        wave_objs.append(res["adj"])
        overall = aggregate_r15_wave_adjudications(wave_objs, sunset_reached=bool(args.sunset_reached))
        (out_dir / f"adjudicate_overall_{run_tag}.json").write_text(json.dumps(overall, indent=2) + "\n", encoding="utf-8")
        (out_dir / f"report_overall_{run_tag}.md").write_text(
            render_r15_overall_markdown_report(overall, report_tag=run_tag), encoding="utf-8"
        )

    print(
        json.dumps(
            {
                "run_tag": run_tag,
                "k_conf": res["adj"].get("k_conf"),
                "Claim1_H_E1": res["adj"].get("Group_E_Claim1_Verdict"),
                "Claim2_H_E2": res["adj"].get("Group_E_Claim2_Verdict"),
                "statement_label": res["adj"].get("statement_label"),
                "pending_lines": res["adj"].get("pending_lines"),
                "elapsed_sec": res["summary"]["elapsed_sec"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
