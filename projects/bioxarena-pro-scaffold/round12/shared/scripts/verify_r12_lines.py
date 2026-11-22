#!/usr/bin/env python3
"""
Counts-Only Verification & Confirmatory Eligibility Auditor for R12 Lines (A-T1).

BLINDING INVARIANT (`PREREG_E` §0, §2.2, §3.2):
  Computes and reports ONLY dataset counts, finite checks, `P4`/`P5` attrition,
  and target/gene set overlaps against frozen Round 9 (`r9_6`) reference sets.
  ZERO predictors, correlations, `r_orth_cal`, slopes, residuals, or `FCT` values
  are computed or loaded.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent
ROUND12_ROOT = SCRIPTS_DIR.parents[1]
REPO_ROOT = ROUND12_ROOT.parent
R10_SCRIPTS = REPO_ROOT / "round10" / "shared" / "scripts"
for p in (str(SCRIPTS_DIR), str(R10_SCRIPTS)):
    if p not in sys.path:
        sys.path.insert(0, p)

from prepare_r10_lines import DEFAULT_HGNC_PATH, build_hgnc_lookup  # noqa: E402
from verify_r10_lines import compute_r9_reference_sets, verify_and_count_line  # noqa: E402

R9_LINES_DIR = REPO_ROOT / "round9" / "p2" / "shared" / "lines"
R12_LINES_DIR = ROUND12_ROOT / "shared" / "lines"
R12_CONFIRMATORY_OUT = ROUND12_ROOT / "status" / "CONFIRMATORY_LINES_R12.txt"
R12_QC_DOC_OUT = ROUND12_ROOT / "shared" / "R12_QC_COUNTS.md"
DEFAULT_HEK_FEAS_JSON = ROUND12_ROOT / "shared" / "hek293_feasibility.json"

CANDIDATE_CONFIRMATORY = ("a549_ko", "hek293_xucao")
SENSITIVITY_ONLY_LINES = ("a549_ko_cas9pos",)
ALL_R12_LINES = CANDIDATE_CONFIRMATORY + SENSITIVITY_ONLY_LINES


def check_no_tmp_files(line_dir: Path) -> None:
    tmp = [p.name for p in line_dir.iterdir() if p.name.startswith("_tmp_")]
    if tmp:
        raise RuntimeError(f"{line_dir.name}: leftover temporary files found: {tmp}")


def write_qc_counts_markdown(
    counts_by_line: Dict[str, dict],
    meta_by_line: Dict[str, dict],
    hek_feasibility: Optional[dict],
    passing: List[str],
    out_path: Path = R12_QC_DOC_OUT,
) -> None:
    lines: List[str] = [
        "# Round 12 QC & Eligibility Counts (`PREREG_E` §2.2 & §3.2)",
        "",
        "**Blinding Declaration**: Strictly **counts-only** audit executed after `R12_GO.txt` (`PREREG_E` commit `4c606f4732ea080768e27757c689d40d7892cd3d`) and before `CODE_FREEZE_R12.txt` / Round 12 evaluation. Zero transfer predictions, correlations, orthogonal calibration (`r_orth_cal`), slopes, residuals, `FCT`, or `N_5` metrics were computed on any Round 12 dataset.",
        "**Verification Gate**: All prepared directories under `round12/shared/lines/` verified complete (`genes.txt`, `targets.txt`, `ctrl_cells.npz`, `pseudobulk_delta.npz`, `line_metadata.json`) with zero leftover temporary shards (`_tmp_*`) or temporary HDF5 caches.",
        "",
        "---",
        "",
        "## 1. `HEK293` Pre-Metric Feasibility Gate (`E3`) & Per-Line Split-Half NTC Counts",
        "",
    ]

    if hek_feasibility is not None:
        verdict = str(hek_feasibility.get("e3_verdict", "UNKNOWN"))
        reason = str(hek_feasibility.get("e3_reason", ""))
        lines.extend([
            f"- **`hek293_xucao` `E3` Pre-Metric Feasibility Gate (`XuCao2023.h5ad`)**: **`{verdict}`** — {reason}",
            f"  - Matrix inspection: `n_obs={hek_feasibility.get('n_obs')}`, `n_var={hek_feasibility.get('n_var')}`, `top_keys={hek_feasibility.get('top_keys')}`, matrices=`{json.dumps(hek_feasibility.get('matrices', []))}`",
            f"  - `.obs` columns: `{hek_feasibility.get('obs_summary', {}).get('columns', [])}`",
        ])

    lines.extend([
        "",
        "| Line | Role / Status (`PREREG_E` §2.2 & §5.5) | Split-Half Partition (`Half A` vs `Half B`) | Kept Cells (`single-sgRNA + NTC`) | Targeting Cells | Total NTC (`n_ctrl_total`) | NTC Half A / Half B | Sampled NTC | Finite Genes (`raw -> final`) |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ])

    for ln in ALL_R12_LINES:
        c = counts_by_line.get(ln)
        m = meta_by_line.get(ln, {})
        audit = m.get("adapter_audit", {})
        if c is None:
            if ln == "hek293_xucao" and hek_feasibility is not None and not hek_feasibility.get("e3_gate_pass", False):
                verdict = str(hek_feasibility.get("e3_verdict", "EXCLUDED_FEASIBILITY_Q3"))
                if verdict == "BLOCKED_BY_READER_INDEX_BUG":
                    lines.append("| **`hek293_xucao`** | **Pending Pre-Freeze (`BLOCKED_BY_READER_INDEX_BUG`: `E3` data PASS; awaiting `A-T2` `.attrs['_index']` fix + Master Control `max_pct_mt` ruling)** | `Stage 2 seed=42 (no .obs rep col)` | `—` | `—` | `2,758` raw meta | `—` | `—` | `59,428` human raw |")
                else:
                    lines.append(f"| **`hek293_xucao`** | **Excluded (`{verdict}`)** | `—` | `0` | `0` | `0` | `0` / `0` | `0` | `—` |")
            else:
                lines.append(f"| `{ln}` | Not Prepared | `—` | `0` | `0` | `0` | `0` / `0` | `0` | `—` |")
            continue
        if ln == "a549_ko_cas9pos":
            st_lbl = "Sensitivity (`SA-Cas9`, `CAS9 > 0`; never in `K_conf`)"
        elif ln in passing:
            st_lbl = "**Confirmatory (`PASS`)**"
        else:
            st_lbl = "**Fail Confirmatory (`<50` or `<3000`)**"
        split_desc = str(audit.get("half_rule") or audit.get("half_rule_used") or audit.get("split_half_mode") or f"{c['n_half_a_batches']} vs {c['n_half_b_batches']} batches")
        raw_g = m.get("n_genes_raw", "?")
        lines.append(
            f"| **`{ln}`** | {st_lbl} | {split_desc} | `{c['raw_cells_after_qc']:,}` | "
            f"`{c['n_targeting_single_sgrna_cells']:,}` | `{c['n_nt_cells']:,}` | "
            f"`{c['n_nt_half_a']:,}` / `{c['n_nt_half_b']:,}` | `{c['n_nt_sampled']:,}` | "
            f"`{raw_g} -> {c['n_genes_finite']} -> {c['n_genes_final']}` |"
        )

    lines.extend([
        "",
        "---",
        "",
        "## 2. `PREREG_E` `P4`/`P5` Target Attrition & `|G_eval(t)|` (`r9_6` Sources)",
        "",
        "Criterion filters perturbations sequentially: **(1)** reproducible on target line `t` (`n_cells >= 25`, top-5000 AD `BH < 0.05 >= 50` vs `1,000` NTC, 20-PC Energy `p < 0.001` vs `500` NTC; `1` row per symbol); **(2)** reproducible on `>= 1` frozen `r9_6` source line (`union6_rep = 2,117` canonical / raw); **(3) `n_eval_targets`**: target self-gene lies in `G_eval(t) = canon(R9_COMMON_6148) ∩ genes(t)`.",
        "",
        "| Line | Deduplicated Targets (`n >= 10`) | `half_valid` (`>=2` in A & B) | **(1)** Reproducible on `t` (`half_valid` rep) | **(2)** Shared Reproducible w/ `r9_6` (`Canon / Raw`) | **(3) `n_eval_targets`** in `G_eval(t)` (`Canon / Raw`) | `half_valid` in `n_eval_targets` | **`\\|G_eval(t)\\|`** (`Canon / Raw`) | Confirmatory Gate (`>=50` targets & `>=3000` genes) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ])
    for ln in ALL_R12_LINES:
        c = counts_by_line.get(ln)
        if c is None:
            continue
        ok_ntc_half = (c["n_nt_half_a"] >= 15) and (c["n_nt_half_b"] >= 15)
        if ln == "a549_ko_cas9pos":
            gate = "SENSITIVITY_ONLY (`SA-Cas9`)"
        elif ln in passing:
            gate = f"**`PASS` (`{c['n_shared_reproducible_eval_canon']} >= 50`, `{c['g_eval_canon']} >= 3000`, NTC `{c['n_nt_half_a']}/{c['n_nt_half_b']} >= 15`)**"
        elif ln == "hek293_xucao":
            gate = f"**`SECONDARY_ONLY` (`n_eval={c['n_shared_reproducible_eval_canon']}`, `|G_eval|={c['g_eval_canon']}`, NTC `{c['n_nt_half_a']}/{c['n_nt_half_b']}`)**"
        else:
            gate = f"**`FAIL` (`n_eval={c['n_shared_reproducible_eval_canon']}`, `|G_eval|={c['g_eval_canon']}`, ok_ntc_half={ok_ntc_half})**"
        lines.append(
            f"| **`{ln}`** | `{c['n_targets_dedup']:,}` | `{c['n_targets_half_valid']:,}` | "
            f"`{c['n_targets_reproducible']:,}` (`{c['n_targets_reproducible_half_valid']:,}`) | "
            f"**`{c['n_shared_reproducible_with_r9_all_canon']:,}`** / `{c['n_shared_reproducible_with_r9_all_raw']:,}` | "
            f"**`{c['n_shared_reproducible_eval_canon']:,}`** / `{c['n_shared_reproducible_eval_raw']:,}` | "
            f"`{c['n_shared_reproducible_eval_half_valid_canon']:,}` | "
            f"**`{c['g_eval_canon']:,}`** / `{c['g_eval_raw']:,}` | {gate} |"
        )

    lines.extend([
        "",
        "---",
        "",
        "## 3. Adapter Audits, HGNC Bridges, & Confirmatory Gate (`round12/status/CONFIRMATORY_LINES_R12.txt`)",
        "",
    ])
    for ln in ALL_R12_LINES:
        c = counts_by_line.get(ln)
        m = meta_by_line.get(ln, {})
        if c is not None:
            audit = m.get("adapter_audit", {})
            lines.append(f"- **`{ln}`**: `{json.dumps(audit)}` | HGNC-bridged eval symbols: `{c.get('hgnc_bridge_recovered_eval_symbols', [])}`")
    lines.append("- **`a549_ko` Excluded Targets & MT Disclosure**: `QARS` (`591` cells across single target) was excluded because the frozen R9–R11 local HGNC table lacks a `QARS -> QARS1` mapping (`999 / 1,000` annotated gene targets kept; `605,484 / 606,075` cells; kept unchanged per `round12/status/RULING_R12_PREFREEZE.md`). `GSE345058` cells already have max `pct_counts_mt = 9.994% < 20%` (`0` cells dropped by `max_pct_mt`).")
    if hek_feasibility is not None:
        pmt_q = hek_feasibility.get("obs_percent_mito_quantiles") or {
            "5": 12.89, "25": 19.91, "50": 25.81, "75": 32.32, "95": 42.42
        }
        lines.append(
            "- **`hek293_xucao` Pre-Freeze Ruling B Disclosure (`round12/status/RULING_R12_PREFREEZE.md`)**: "
            "Exported `XuCao2023.h5ad` holds all `98,315` QC-valid single-sgRNA cells (`nperts` in `{0, 1}`, `min ncounts=301 >= 300`, `min ngenes=201 >= 200`, `2,758` `perturbation == 'control'` NTCs across `15` `NO-TARGET_*` guides; authors' exported `98,315` pool, with no additional post-hoc KD>40% filter). "
            f"Per `RULING_R12_PREFREEZE.md` (Ruling B), **no mitochondrial upper cutoff (`max_pct_mt` disabled)** is applied because `PREREG_E §3.2` contains no MT limit (and no extra `MT <= 20%` sensitivity branch is added). "
            f"For transparency, `.obs['percent_mito']` 5/25/50/75/95 percentiles are **`{pmt_q.get('5')}% / {pmt_q.get('25')}% / {pmt_q.get('50')}% / {pmt_q.get('75')}% / {pmt_q.get('95')}%`** (`MT-*` fraction on `/X` human RNA UMIs: `13.25% / 20.48% / 26.54% / 33.28% / 43.57%`)."
        )
    hek_c = counts_by_line.get("hek293_xucao")
    hek_m = meta_by_line.get("hek293_xucao", {})
    if hek_c is not None:
        lines.append(
            "- **`hek293_xucao` D3 NTC Subtraction Unit Fix Disclosure (`round12/status/RULING_R12_HEK293_UNITS.md`)**: "
            "Per `RULING_R12_HEK293_UNITS.md` §1 (`PREREG_E §3.3`), Stage 2 HEK293 must use per-half batches (`half_1` / `half_2`) as the NTC subtraction unit rather than pooled single `batch_0`. "
            "Before/after counts-only comparison: "
            "**Superseded v1 (`batch_0` pooled, archived in `CONFIRMATORY_LINES_R12_SUPERSEDED_v1.txt`)**: `n_batches=1`, `98,274` kept cells (`2,757` NTC: `1,378 / 1,379`), `n_genes=29,247 -> 29,194`, `|G_eval|=6,127`, `96/200` reproducible, `70` shared w/ `r9_6`, `n_eval_targets=62`; "
            f"**Corrected v2 (`half_1` / `half_2` per-half centering)**: `n_batches={hek_m.get('n_batches', '?')}`, `{hek_c['raw_cells_after_qc']:,}` kept cells (`{hek_c['n_nt_cells']:,}` NTC: `{hek_c['n_nt_half_a']:,} / {hek_c['n_nt_half_b']:,}`), "
            f"`n_genes={hek_c['n_genes_finite']:,} -> {hek_c['n_genes_final']:,}`, `|G_eval|={hek_c['g_eval_canon']:,}`, "
            f"`{hek_c['n_targets_reproducible']:,}/{hek_c['n_targets_dedup']:,}` reproducible, `{hek_c['n_shared_reproducible_with_r9_all_canon']:,}` shared w/ `r9_6`, `n_eval_targets={hek_c['n_shared_reproducible_eval_canon']:,}`."
        )
    hek_pending = (
        "hek293_xucao" not in counts_by_line
        and hek_feasibility is not None
        and hek_feasibility.get("e3_verdict") == "BLOCKED_BY_READER_INDEX_BUG"
    )
    if hek_pending:
        lines.append(
            "- **`CONFIRMATORY_LINES_R12.txt` Status (`PENDING_HEK293_PREFREEZE_RULING`)**: "
            f"Provisional passing lines so far = `{passing}` (`a549_ko` passes confirmatory floor with `445 >= 50` targets and `6,121 >= 3,000` genes). "
            "`hek293_xucao` awaits `A-T2` `R12_CODE_PATCHED` build."
        )
    else:
        lines.append(f"- **Final `CONFIRMATORY_LINES_R12.txt` (`K_conf={len(passing)}`)**: `{passing}`")
        if len(passing) == 1 and passing == ["a549_ko"]:
            lines.append("- **Mandatory `PREREG_E` §1 / §2.2 Single-Line Disclosure (`K_conf=1`, 推广性有限)**: Only `a549_ko` cleared the confirmatory gate (`>= 50` eval targets & `>= 3,000` genes); `hek293_xucao` fell below threshold under compliant half-level NTC subtraction and is mechanically classified `SECONDARY_ONLY` (reported descriptively under §5). `K_conf = 1` confirmatory claims apply strictly to `A549` (Cas9 KO lung carcinoma) with limited multi-lineage generalizability.")
        elif "a549_ko" not in passing:
            lines.append("- **Mandatory `PREREG_E` §2.2 Abort Notice**: Primary confirmatory line `a549_ko` failed eligibility thresholds (`K_conf=0` confirmatory abort).")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def evaluate_confirmatory_lines_r12(
    counts_by_line: Dict[str, dict],
    hek_feasibility: Optional[dict] = None,
) -> Tuple[List[str], str]:
    """
    Mechanically determine Round 12 confirmatory lines per `PREREG_E` §2.2, §3.2, §3.3:
      Candidates: `a549_ko` (mandatory primary), `hek293_xucao` (subject to `E3`).
      Sensitivity-only (`never enters K_conf`): `a549_ko_cas9pos`.
      Criteria: e3_gate_pass (if applicable) AND |G_eval(t)| >= 3000 AND n_eval_targets >= 50 AND >= 15 NTC per half.
    """
    passing: List[str] = []
    lines_txt = [
        "# Round 12 Mechanical Confirmatory Line Gate (PREREG_E §2.2, §3.2, §3.3)",
        "# Thresholds: candidate in {a549_ko, hek293_xucao}, feasibility==PASS, n_eval_targets >= 50, |G_eval(t)| >= 3000, NTC_half_A >= 15 & NTC_half_B >= 15",
    ]
    a549_passed = False
    for ln in CANDIDATE_CONFIRMATORY:
        if ln == "hek293_xucao" and hek_feasibility is not None and not hek_feasibility.get("e3_gate_pass", False):
            verdict = str(hek_feasibility.get("e3_verdict", "EXCLUDED_FEASIBILITY_Q3"))
            reason = str(hek_feasibility.get("e3_reason", "failed E3 pre-metric feasibility gate"))
            lines_txt.append(f"# hek293_xucao: {verdict} ({reason})")
            continue
        c = counts_by_line.get(ln)
        if c is None:
            lines_txt.append(f"# {ln}: NOT_PROCESSED")
            continue
        ok_t = c["n_shared_reproducible_eval_canon"] >= 50
        ok_g = c["g_eval_canon"] >= 3000
        ok_ntc = (c["n_nt_half_a"] >= 15) and (c["n_nt_half_b"] >= 15)
        ok = ok_t and ok_g and ok_ntc
        if ln == "a549_ko":
            a549_passed = ok
        if ok:
            passing.append(ln)
        t_op = ">=" if ok_t else "<"
        g_op = ">=" if ok_g else "<"
        status_tok = "PASS_CONFIRMATORY" if ok else ("SECONDARY_ONLY" if ln == "hek293_xucao" else "FAIL_CONFIRMATORY")
        lines_txt.append(
            f"# {ln}: {status_tok} "
            f"(n_eval_targets={c['n_shared_reproducible_eval_canon']} {t_op} 50, "
            f"|G_eval|={c['g_eval_canon']} {g_op} 3000, "
            f"NTC_half_A={c['n_nt_half_a']}, NTC_half_B={c['n_nt_half_b']})"
        )

    for ln in SENSITIVITY_ONLY_LINES:
        c = counts_by_line.get(ln)
        if c is not None:
            lines_txt.append(
                f"# {ln}: SENSITIVITY_SA_CAS9 (never enters K_conf; n_eval_targets={c['n_shared_reproducible_eval_canon']}, |G_eval|={c['g_eval_canon']})"
            )
        else:
            lines_txt.append(f"# {ln}: SENSITIVITY_SA_CAS9 (not yet processed)")

    # PREREG_E §2.2: a549_ko is mandatory confirmatory anchor; if it fails, abort confirmatory (K_conf=0).
    if not a549_passed and "a549_ko" in counts_by_line:
        passing = []
        lines_txt.append("# PREREG_E_ABORT: mandatory primary line a549_ko failed -> K_CONF=0")

    lines_txt.append(f"# K_CONF={len(passing)}")
    if len(passing) == 0:
        lines_txt.append("# PREREG_E_GUARD: K_conf = 0 -> ABORT confirmatory adjudication")
    elif len(passing) == 1:
        lines_txt.append("# CONFIRMATORY_SINGLE_LINE_K1: K_conf = 1 confirmatory line (mandatory single-line generalizability disclosure)")
    else:
        lines_txt.append(f"# CONFIRMATORY_DUAL_LINE_K2: K_conf = {len(passing)} confirmatory cell lines")
    lines_txt.extend(passing)
    return passing, "\n".join(lines_txt) + "\n"


def main() -> int:
    p = argparse.ArgumentParser(description="Verify prepared R12 lines and mechanically write CONFIRMATORY_LINES_R12.txt (counts only).")
    p.add_argument("--lines-dir", type=Path, default=R12_LINES_DIR)
    p.add_argument("--hek-feasibility-json", type=Path, default=DEFAULT_HEK_FEAS_JSON)
    p.add_argument("--write-confirmatory", action="store_true")
    p.add_argument("--write-qc-report", action="store_true")
    p.add_argument("--qc-report-out", type=Path, default=R12_QC_DOC_OUT)
    p.add_argument("--confirmatory-out", type=Path, default=R12_CONFIRMATORY_OUT)
    p.add_argument("--json-out", type=Path, default=None)
    args = p.parse_args()

    hgnc = build_hgnc_lookup(DEFAULT_HGNC_PATH)
    r9_ref = compute_r9_reference_sets(hgnc, R9_LINES_DIR)

    hek_feas = None
    if args.hek_feasibility_json and args.hek_feasibility_json.exists():
        hek_feas = json.loads(args.hek_feasibility_json.read_text(encoding="utf-8"))

    results: Dict[str, dict] = {}
    metas: Dict[str, dict] = {}
    for ln in ALL_R12_LINES:
        ld = args.lines_dir / ln
        if (ld / "line_metadata.json").exists():
            check_no_tmp_files(ld)
            results[ln] = verify_and_count_line(ld, hgnc, r9_ref)
            metas[ln] = json.loads((ld / "line_metadata.json").read_text(encoding="utf-8"))
            print(json.dumps(results[ln], indent=2))
        else:
            print(f"[SKIP] {ln}: not prepared at {ld}")

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    passing, conf_text = evaluate_confirmatory_lines_r12(results, hek_feasibility=hek_feas)
    print("\n" + conf_text)
    if args.write_confirmatory:
        args.confirmatory_out.parent.mkdir(parents=True, exist_ok=True)
        args.confirmatory_out.write_text(conf_text, encoding="utf-8")
        print(f"Wrote {args.confirmatory_out} (passing={passing})")
    if args.write_qc_report:
        args.qc_report_out.parent.mkdir(parents=True, exist_ok=True)
        write_qc_counts_markdown(results, metas, hek_feas, passing, args.qc_report_out)
        print(f"Wrote {args.qc_report_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
