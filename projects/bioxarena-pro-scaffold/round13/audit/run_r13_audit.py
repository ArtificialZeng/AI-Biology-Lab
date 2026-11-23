#!/usr/bin/env python3
"""
Independent Round 13 Pre-Registration Audit Runner & Report Generator (A-T3).

Strictly follows locked `PREREG_F.md` (`commit 5b51e5b`) and `PREREG_E.md`:
  - Reuses frozen `round12/audit/independent_r12_metrics.py::evaluate_target_line_r12`
    and `round11/audit/` primitives with Round 13 seed `20261223`, `n_boot = 4,000`,
    `n_perm = 100`, two-sided `97.5%` percentile CIs (`[1.25%, 98.75%]`).
  - Strictly NEVER reads or imports `round13/src/` or `round13/results/`.
  - Runtime blinding lock blocks `--phase holdout_r13` until `PREREG_F_COMMIT.txt`,
    `CODE_FREEZE_R13.txt`, and `CONFIRMATORY_LINES_R13.txt` exist and are valid.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np

AUDIT_R13_DIR = Path(__file__).resolve().parent
REPO_ROOT = AUDIT_R13_DIR.parents[1]
R11_AUDIT_DIR = REPO_ROOT / "round11" / "audit"
R12_AUDIT_DIR = REPO_ROOT / "round12" / "audit"
for _p in (str(AUDIT_R13_DIR), str(R11_AUDIT_DIR), str(R12_AUDIT_DIR)):
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
from independent_r11_predictors import deterministic_subsample, sample_covariance  # noqa: E402
from independent_r12_metrics import evaluate_target_line_r12  # noqa: E402
from independent_r13_loader import (  # noqa: E402
    DEV_ALLOWED_LINES_BY_PHASE,
    R13_LINES_DIR,
    R13_STATUS_DIR,
    enforce_r13_gate,
    parse_confirmatory_lines_r13,
    resolve_line_cache_dir_r13,
)
from run_r12_audit import _fmt, _write_csvs  # noqa: E402

DEFAULT_R13_SEED = 20261223


def adjudicate_round13(
    per_line_summary: Dict[str, Dict[str, Any]],
    gate_info: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Apply `PREREG_F.md §2, §5, §6, §9.3` adjudication:
      - Maximum `K_conf = 1` (`mac_ctrl` ONLY; `mac_lps` is ALWAYS `SECONDARY_ONLY`).
      - If `mac_ctrl` failed feasibility (`EXCLUDED_FEASIBILITY`) or count gates
        (`SECONDARY_ONLY`, `k_conf == 0`): ZERO confirmatory adjudication; report
        all feasible lines strictly descriptively under `§5` (`R12 A549` unchanged).
      - If `mac_ctrl` entered `CONFIRMATORY` (`k_conf == 1`): adjudicate `H-E1`
        (`H-E1a AND H-E1b AND H-E1c`) and `H-E2` (`H-E2a AND H-E2b`) independently
        on `mac_ctrl` at `97.5%` CI; strictly prohibit pooling R12 and R13 data.
    """
    conf_lines = list(gate_info.get("confirmatory_lines", []))
    roles = dict(gate_info.get("roles", {}))
    mac_ctrl_role = roles.get("mac_ctrl", "CONFIRMATORY" if "mac_ctrl" in conf_lines else "SECONDARY_ONLY")
    mac_lps_role = roles.get("mac_lps", "SECONDARY_ONLY")
    if mac_lps_role == "CONFIRMATORY" or "mac_lps" in conf_lines:
        raise ValueError("PREREG_F §2 violation: mac_lps can never enter confirmatory adjudication.")
    if len(conf_lines) > 1 or (len(conf_lines) == 1 and conf_lines[0] != "mac_ctrl"):
        raise ValueError(f"PREREG_F §6 violation: invalid confirmatory_lines {conf_lines!r}")

    if not conf_lines:
        return {
            "status": "DESCRIPTIVE_ONLY_GATE_NOT_PASSED_K0",
            "k_conf": 0,
            "confirmatory_lines": [],
            "mac_ctrl_role": mac_ctrl_role,
            "mac_lps_role": mac_lps_role,
            "claim1_he1_verdict": "NO_CONFIRMATORY_ADJUDICATION_GATE_FAILED",
            "claim2_he2_verdict": "NO_CONFIRMATORY_ADJUDICATION_GATE_FAILED",
            "claim1_confirmed": False,
            "claim2_confirmed": False,
            "no_cross_round_pooling_notice": (
                "Round 13 performs zero confirmatory adjudication (K_conf = 0) because mac_ctrl "
                f"status = {mac_ctrl_role}. Descriptive metrics only; Round 12 (A549) remains unchanged."
            ),
        }

    if "mac_ctrl" not in per_line_summary:
        raise RuntimeError("Missing mac_ctrl summary in adjudicate_round13")

    m_res = per_line_summary["mac_ctrl"]
    c1_ok = bool(m_res["he1_pass"])
    c2_ok = bool(m_res["he2_pass"])

    if c1_ok and c2_ok:
        interp = (
            "BOTH CLAIMS CONFIRMED ON mac_ctrl (K_conf = 1 in Round 13). State ONLY as "
            "'confirmed independently on two separate cell lines in two separately preregistered "
            "prospective tests' (Round 12 A549 + Round 13 mac_ctrl). NEVER pool R12 and R13 data."
        )
    else:
        failed_parts = ([] if c1_ok else ["Claim 1 (H-E1)"]) + ([] if c2_ok else ["Claim 2 (H-E2)"])
        interp = (
            f"FAILED TO REPLICATE ON SECOND CELL LINE (mac_ctrl): {', '.join(failed_parts)} not "
            "confirmed. Report failure honestly; note confounding with macrophage biology, "
            "Flex probe panel, donor RG5581, and power; Round 12 A549 remains unchanged."
        )

    return {
        "status": "CONFIRMATORY_R13_K1_MAC_CTRL",
        "k_conf": 1,
        "confirmatory_lines": ["mac_ctrl"],
        "mac_ctrl_role": mac_ctrl_role,
        "mac_lps_role": mac_lps_role,
        "claim1_he1_per_line": {"mac_ctrl": c1_ok},
        "claim2_he2_per_line": {"mac_ctrl": c2_ok},
        "claim1_confirmed": c1_ok,
        "claim2_confirmed": c2_ok,
        "claim1_he1_verdict": "CONFIRMED" if c1_ok else "FAILED_TO_REPLICATE",
        "claim2_he2_verdict": "CONFIRMED" if c2_ok else "FAILED_TO_REPLICATE",
        "no_cross_round_pooling_notice": interp,
    }


def render_r13_audit_markdown(summary: Dict[str, Any]) -> str:
    """Emit concise clean-room Markdown audit report adhering to `PREREG_F §6–§7`."""
    adj = summary["adjudication"]
    per_line = summary["per_line"]
    lines: List[str] = [
        f"# Round 13 Independent Clean-Room Audit Report (`{summary['phase']}`, `{summary['source_mode']}`)",
        "",
        f"- **PREREG_F Commit:** `{summary['gate_meta'].get('prereg_commit', '')}`",
        f"- **Seed / Boot / Perm:** `seed = {summary['seed']}`, `n_boot = {summary['n_boot']}`, `n_perm = {summary['n_perm']}`, `97.5%` CIs",
        f"- **Adjudication Status (`K_conf = {adj.get('k_conf', 0)}`):** `{adj['status']}`",
        f"- **Claim 1 (`H-E1`) Verdict:** `{adj['claim1_he1_verdict']}`",
        f"- **Claim 2 (`H-E2` / `SA-8`) Verdict:** `{adj['claim2_he2_verdict']}`",
        f"- **Interpretive Lock:** {adj.get('no_cross_round_pooling_notice', adj.get('mandatory_disclosure', ''))}",
        "",
        "| Line | `|G_eval|` | `n_eval_targets` | `FCT r_orth_cal [97.5% CI]` | `H-E1b` (`pass/total`) | `Perm Q97.5` | `H-E1` | `SA-8 E1 split` | `E1(split)-E1(all)` | `H-E2` |",
        "|---|---:|---:|---|---:|---:|:---:|---|---:|:---:|",
    ]
    for ln, d in per_line.items():
        fct = d["method_eval"]["fct"]["all_nontarget"]
        sa8 = d["sa8_e1"]
        d_sa8 = d["sa8_e1_split_minus_all"]
        perm_q975 = d["perm100_all_nontarget"]["r_orth_cal"]["p975"]
        lines.append(
            f"| `{ln}` | `{d['n_genes_eval']}` | `{d['n_eval_targets']}` | "
            f"`{_fmt(fct.get('r_orth_cal_est'))}` `[{_fmt(fct.get('r_orth_cal_ci975_lo'))}, {_fmt(fct.get('r_orth_cal_ci975_hi'))}]` | "
            f"`{d['he1b_n_pass']}/{d['he1b_n_total']}` | `{_fmt(perm_q975)}` | `{'PASS' if d['he1_pass'] else 'FAIL'}` | "
            f"`{_fmt(sa8.get('r_orth_cal_split_est'))}` `[{_fmt(sa8.get('r_orth_cal_split_ci975_lo'))}, {_fmt(sa8.get('r_orth_cal_split_ci975_hi'))}]` | "
            f"`{_fmt(d_sa8.get('diff_est'))}` `[{_fmt(d_sa8.get('diff_ci975_lo'))}, {_fmt(d_sa8.get('diff_ci975_hi'))}]` | `{'PASS' if d['he2_pass'] else 'FAIL'}` |"
        )
    lines.append("")
    lines.append("> Note: `E1_ConstMed` is a target-baseline-free transfer predictor; Claim 2 (`SA-8`) is an empirical observation on `r_orth_cal` with unknown mechanism.")
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Independent Round 13 clean-room audit runner (A-T3).")
    p.add_argument(
        "--phase",
        required=True,
        choices=["dev_r9", "dev_r10", "dev_r11", "dev_r12", "holdout_r13"],
    )
    p.add_argument("--source-mode", default="r9_6", choices=["r9_6", "r9_r10_8"])
    p.add_argument("--targets", nargs="+", default=None)
    p.add_argument("--seed", type=int, default=DEFAULT_R13_SEED)
    p.add_argument("--n-boot", type=int, default=4000)
    p.add_argument("--n-perm", type=int, default=100)
    p.add_argument("--cov-max-cells", type=int, default=20000)
    p.add_argument("--status-dir", type=Path, default=R13_STATUS_DIR)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--skip-n5", action="store_true", default=False)
    p.add_argument("--skip-n5-split", action="store_true", default=False)
    p.add_argument("--skip-r19", action="store_true", default=False)
    p.add_argument("--skip-s7", action="store_true", default=False)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    is_dev = args.phase in DEV_ALLOWED_LINES_BY_PHASE

    if args.targets is not None:
        targets = list(args.targets)
    elif is_dev:
        targets = list(DEV_ALLOWED_LINES_BY_PHASE[args.phase])
    else:
        conf_info = parse_confirmatory_lines_r13(args.status_dir / "CONFIRMATORY_LINES_R13.txt")
        targets = [ln for ln in ("mac_ctrl", "mac_lps") if conf_info["roles"].get(ln) != "EXCLUDED_FEASIBILITY"]

    extra_dirs = [R13_LINES_DIR] if (not is_dev or any(t in {"mac_ctrl", "mac_lps"} for t in targets)) else []
    gate_meta = enforce_r13_gate(
        phase=args.phase,
        target_lines=targets,
        status_dir=args.status_dir,
        extra_dirs=extra_dirs,
        require_confirmatory_file=not is_dev,
    )
    gate_info = (
        {"roles": {}, "confirmatory_lines": [], "k_conf": 0}
        if is_dev
        else parse_confirmatory_lines_r13(args.status_dir / "CONFIRMATORY_LINES_R13.txt")
    )

    hgnc = CanonicalHGNC.from_tsv()
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
            min_ov = 3000 if (not is_dev and tgt in gate_info["confirmatory_lines"]) else 50
            tdata = load_line_dataset(
                resolve_line_cache_dir_r13(tgt, args.phase),
                representation="log1p",
                aligned_genes=common_genes,
                hgnc=hgnc,
                min_gene_overlap=min_ov,
            )
            cov_cache[tgt] = sample_covariance(deterministic_subsample(tdata.ctrl_cells, args.cov_max_cells, seed=42), shrink=0.0)
            gene_order_cache[tgt] = tdata.genes

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
            "mandatory_disclosure": f"Pre-freeze verification strictly on unblinded historical {args.phase} lines.",
        }
    else:
        adjudication = adjudicate_round13(per_line, gate_info)

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
    (args.out_dir / f"audit_r13_summary_{tag}.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (args.out_dir / f"AUDIT_R13_REPORT_{tag}.md").write_text(
        render_r13_audit_markdown(summary),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": f"R13_AUDIT_{args.phase.upper()}_COMPLETE",
                "source_mode": args.source_mode,
                "targets": targets,
                "adjudication_status": adjudication["status"],
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
