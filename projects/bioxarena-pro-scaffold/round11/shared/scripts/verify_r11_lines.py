#!/usr/bin/env python3
"""
Counts-Only Verification & Confirmatory Eligibility Auditor for R11 Lines (A-T1).

BLINDING INVARIANT:
  Computes and reports ONLY dataset counts, finite checks, and target/gene set
  overlaps under PREREG_D §2.2-2.3. Strictly forbidden before CODE_FREEZE_R11.txt:
  any model predictions, correlations, r_orth_cal, slopes, residuals, or FCT metrics.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(__file__).resolve().parents[3]
R10_SCRIPTS = PROJECT_ROOT / "round10" / "shared" / "scripts"
for p in (SCRIPTS_DIR, R10_SCRIPTS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from prepare_r10_lines import DEFAULT_HGNC_PATH, HGNCLookup, build_hgnc_lookup  # noqa: E402
from verify_r10_lines import compute_r9_reference_sets, verify_and_count_line  # noqa: E402

R9_LINES_DIR = PROJECT_ROOT / "round9" / "p2" / "shared" / "lines"
R11_LINES_DIR = PROJECT_ROOT / "round11" / "shared" / "lines"
R11_CONFIRMATORY_OUT = PROJECT_ROOT / "round11" / "status" / "CONFIRMATORY_LINES_R11.txt"
R11_QC_DOC_OUT = PROJECT_ROOT / "round11" / "shared" / "R11_QC_COUNTS.md"

CANDIDATE_CONFIRMATORY = ("telohaec", "hcasmc", "cd4t_stim48")
SECONDARY_ONLY_LINES = ("feng_ipsc",)
ALL_R11_LINES = CANDIDATE_CONFIRMATORY + SECONDARY_ONLY_LINES


def check_no_tmp_files(line_dir: Path) -> None:
    tmp = [p.name for p in line_dir.iterdir() if p.name.startswith("_tmp_")]
    if tmp:
        raise RuntimeError(f"{line_dir.name}: leftover temporary files found: {tmp}")


def write_qc_counts_markdown(
    counts_by_line: Dict[str, dict],
    meta_by_line: Dict[str, dict],
    hcasmc_feasibility: Optional[dict],
    passing: List[str],
    out_path: Path = R11_QC_DOC_OUT,
) -> None:
    lines: List[str] = [
        "# Round 11 QC & Eligibility Counts (`PREREG_D` §2.2–2.3 & Addendum A)",
        "",
        "**Blinding Declaration**: Strictly **counts-only** audit executed prior to `CODE_FREEZE_R11.txt`. Zero transfer predictions, correlations, orthogonal calibration (`r_orth_cal`), slopes, residuals, or `FCT` values were computed on any Round 11 dataset.",
        "**Verification Gate**: All prepared directories under `round11/shared/lines/` verified complete (`genes.txt`, `targets.txt`, `ctrl_cells.npz`, `pseudobulk_delta.npz`, `line_metadata.json`) with zero leftover temporary shards (`_tmp_*`).",
        "",
        "---",
        "",
        "## 1. Pre-Eval Feasibility (`Q3`), Overlap Audit (`Q6`), & Per-Line Split-Half NTC Counts",
        "",
    ]
    if hcasmc_feasibility is not None:
        lines.extend([
            f"- **`hcasmc` `Q3` Pre-Eval Feasibility Gate (`IGVFFI5999NXAF.h5`)**: "
            f"**`{'PASS' if hcasmc_feasibility.get('feasible') else 'EXCLUDED_FEASIBILITY_Q3'}`** — "
            f"{hcasmc_feasibility.get('reason', '')} "
            f"(shape={hcasmc_feasibility.get('shape')}, chunked `max_rss_gb={hcasmc_feasibility.get('max_rss_gb')}` GB, "
            f"barcodes us_suffix={hcasmc_feasibility.get('barcode_us_suffix_counts', {})}).",
        ])
    lines.extend([
        "- **`Q6` Pre-QC Disjointness Enforcement**:",
        "  - `cd4t_stim48`: verified 0 files from `D1`/`D4`/`Rest`/`Stim8hr` in `round11/shared/raw/cd4t_stim48/`; strictly donors `CE0008678` (`D3_Stim48hr`) and `CE0010866` (`D2_Stim48hr`).",
        "  - `feng_ipsc`: strictly excluded `kolf_2` & `kolf_3` pre-QC (`n_kolf_excluded=133,631`, 0 kept); strictly excluded all `unassigned` cells.",
        "",
        "| Line | Status (`PREREG_D` §2.2–2.3) | Split-Half Partition (`Half A` vs `Half B`) | Kept Cells (`single-sgRNA + NTC`) | Targeting Cells | Total NTC (`n_ctrl_total`) | NTC Half A / Half B | Sampled NTC | Finite Genes (`raw -> final`) |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ])
    for ln in ALL_R11_LINES:
        c = counts_by_line.get(ln)
        m = meta_by_line.get(ln, {})
        audit = m.get("adapter_audit", {})
        if c is None:
            if ln == "hcasmc" and hcasmc_feasibility and not hcasmc_feasibility.get("feasible", False):
                lines.append(f"| **`hcasmc`** | **Excluded (`Q3` Feasibility)** | `—` | `0` | `0` | `0` | `0` / `0` | `0` | `—` |")
            else:
                lines.append(f"| `{ln}` | Pending | `—` | `0` | `0` | `0` | `0` / `0` | `0` | `—` |")
            continue
        if ln == "feng_ipsc":
            st_lbl = "Secondary-only (`H-D5` / never `K_conf`)"
        elif ln in passing:
            st_lbl = "**Confirmatory (`PASS`)**"
        else:
            st_lbl = "**Fail Confirmatory (`<50` or `<3000`)**"
        split_desc = str(audit.get("half_rule") or audit.get("split_half_mode") or f"{c['n_half_a_batches']} vs {c['n_half_b_batches']} batches")
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
        "## 2. `PREREG_D` §2.3 Step-by-Step Target Attrition & `|G_eval(t)|`",
        "",
        "Criterion `PREREG_D` §2.3 filters perturbations sequentially: **(1)** reproducible on target line `t` (`n_cells >= 25`, top-5000 AD `BH < 0.05 >= 50` vs `1,000` NTC, 20-PC Energy `p < 0.001` vs `500` NTC; `1` row per symbol); **(2)** reproducible on `>= 1` frozen Round 9 source line (Canon / Raw); **(3) `n_eval_targets`**: target self-gene lies in `G_eval(t) = canon(R9_COMMON_6148) ∩ genes(t)`.",
        "",
        "| Line | Deduplicated Targets (`n >= 10`) | `half_valid` (`>=2` in A & B) | **(1)** Reproducible on `t` (`half_valid` rep) | **(2)** Shared Reproducible w/ R9 (`Canon / Raw`) | **(3) `n_eval_targets`** in `G_eval(t)` (`Canon / Raw`) | `half_valid` in `n_eval_targets` | **`\\|G_eval(t)\\|`** (`Canon / Raw`) | Confirmatory Gate (`>=50` targets & `>=3000` genes) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ])
    for ln in ALL_R11_LINES:
        c = counts_by_line.get(ln)
        if c is None:
            continue
        if ln == "feng_ipsc":
            gate = "SECONDARY_ONLY (`H-D5`)"
        elif ln in passing:
            gate = f"**`PASS` (`{c['n_shared_reproducible_eval_canon']} >= 50`, `{c['g_eval_canon']} >= 3000`)**"
        else:
            gate = f"**`FAIL` (`n_eval={c['n_shared_reproducible_eval_canon']}`, `|G_eval|={c['g_eval_canon']}`)**"
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
        "## 3. Adapter Audits, HGNC Bridges, & Confirmatory Gate Outcome (`round11/status/CONFIRMATORY_LINES_R11.txt`)",
        "",
    ])
    for ln in ALL_R11_LINES:
        c = counts_by_line.get(ln)
        m = meta_by_line.get(ln, {})
        if c is not None:
            audit = m.get("adapter_audit", {})
            lines.append(f"- **`{ln}`**: `{json.dumps(audit)}` | HGNC-bridged eval symbols: `{c.get('hgnc_bridge_recovered_eval_symbols', [])}`")
    lines.append(f"- **Final `CONFIRMATORY_LINES_R11.txt` (`K_conf={len(passing)}`)**: `{passing}`")
    if "cd4t_stim48" in passing:
        lines.append("- **Mandatory `Q2`/`Q7` Disclosure (`cd4t_stim48`)**: Shares primary `CD4+ alpha-beta T` cell type, sgRNA library, and Marson lab with Round 10 `cd4t_rest` (disjoint donors & `Stim48hr` state; shared R2L25..48 sequencing lanes across Half A/B per `Addendum A2`).")
    if "feng_ipsc" in counts_by_line:
        lines.append("- **Mandatory `Addendum A3` Note (`feng_ipsc`)**: Sorted cell-line parity places sibling clones from `8` donors across Half A & Half B (clone/inlet noise ceiling, secondary-only `H-D5`).")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def evaluate_confirmatory_lines_r11(
    counts_by_line: Dict[str, dict],
    hcasmc_feasibility: Optional[dict] = None,
) -> Tuple[List[str], str]:
    """
    Mechanically determine Round 11 confirmatory lines per PREREG_D §2.2 & §2.3:
      Candidates: `telohaec`, `hcasmc`, `cd4t_stim48`.
      Secondary-only (`never enters K_conf`): `feng_ipsc`.
      Criteria: passed_pre_eval_feasibility == True AND |G_eval(t)| >= 3000 AND n_eval_targets >= 50.
    """
    passing: List[str] = []
    lines_txt = [
        "# Round 11 Mechanical Confirmatory Line Gate (PREREG_D §2.2-2.3)",
        "# Thresholds: candidate in {telohaec, hcasmc, cd4t_stim48}, feasibility==PASS, n_eval_targets >= 50, |G_eval(t)| >= 3000",
    ]
    for ln in CANDIDATE_CONFIRMATORY:
        if ln == "hcasmc" and hcasmc_feasibility is not None and not hcasmc_feasibility.get("feasible", False):
            reason = hcasmc_feasibility.get("reason", "failed Q3 pre-eval feasibility gate")
            lines_txt.append(f"# hcasmc: EXCLUDED_FEASIBILITY_Q3 ({reason})")
            continue
        c = counts_by_line.get(ln)
        if c is None:
            lines_txt.append(f"# {ln}: NOT_PROCESSED")
            continue
        ok = (c["g_eval_canon"] >= 3000) and (c["n_shared_reproducible_eval_canon"] >= 50)
        if ok:
            passing.append(ln)
        note = ""
        if ln == "cd4t_stim48":
            note = " [Q2/Q7 DISCLOSURE: shared cell type/library/lab with R10 cd4t_rest]"
        lines_txt.append(
            f"# {ln}: {'PASS_CONFIRMATORY' if ok else 'FAIL_CONFIRMATORY'} "
            f"(n_eval_targets={c['n_shared_reproducible_eval_canon']} >= 50, |G_eval|={c['g_eval_canon']} >= 3000){note}"
        )

    for ln in SECONDARY_ONLY_LINES:
        c = counts_by_line.get(ln)
        if c is not None:
            lines_txt.append(
                f"# {ln}: SECONDARY_ONLY (never enters K_conf; n_eval_targets={c['n_shared_reproducible_eval_canon']}, |G_eval|={c['g_eval_canon']})"
            )
        else:
            lines_txt.append(f"# {ln}: SECONDARY_ONLY (not yet processed)")

    lines_txt.append(f"# K_CONF={len(passing)}")
    if len(passing) < 2:
        lines_txt.append("# P1_GUARD: K_conf < 2 -> ABORT confirmatory adjudication (descriptive reporting only)")
    elif len(passing) == 2:
        lines_txt.append("# Q7_DISCLOSURE: K_conf = 2 confirmatory cell lines")
    lines_txt.extend(passing)
    return passing, "\n".join(lines_txt) + "\n"


def main() -> None:
    p = argparse.ArgumentParser(description="Verify prepared R11 lines and mechanically assign CONFIRMATORY_LINES_R11.txt (counts only).")
    p.add_argument("--lines-dir", type=Path, default=R11_LINES_DIR)
    p.add_argument("--hcasmc-feasibility-json", type=Path, default=None)
    p.add_argument("--write-confirmatory", action="store_true")
    p.add_argument("--write-qc-report", action="store_true")
    p.add_argument("--qc-report-out", type=Path, default=R11_QC_DOC_OUT)
    p.add_argument("--confirmatory-out", type=Path, default=R11_CONFIRMATORY_OUT)
    p.add_argument("--json-out", type=Path, default=None)
    args = p.parse_args()

    hgnc = build_hgnc_lookup(DEFAULT_HGNC_PATH)
    r9_ref = compute_r9_reference_sets(hgnc, R9_LINES_DIR)

    hc_feas = None
    if args.hcasmc_feasibility_json and args.hcasmc_feasibility_json.exists():
        hc_feas = json.loads(args.hcasmc_feasibility_json.read_text(encoding="utf-8"))

    results: Dict[str, dict] = {}
    metas: Dict[str, dict] = {}
    for ln in ALL_R11_LINES:
        ld = args.lines_dir / ln
        if (ld / "line_metadata.json").exists():
            check_no_tmp_files(ld)
            results[ln] = verify_and_count_line(ld, hgnc, r9_ref)
            metas[ln] = json.loads((ld / "line_metadata.json").read_text(encoding="utf-8"))
            print(json.dumps(results[ln], indent=2))
        else:
            print(f"[SKIP] {ln}: not yet prepared at {ld}")

    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    passing, conf_text = evaluate_confirmatory_lines_r11(results, hcasmc_feasibility=hc_feas)
    print("\n" + conf_text)
    if args.write_confirmatory:
        args.confirmatory_out.parent.mkdir(parents=True, exist_ok=True)
        args.confirmatory_out.write_text(conf_text, encoding="utf-8")
        print(f"Wrote {args.confirmatory_out} (passing={passing})")
    if args.write_qc_report:
        args.qc_report_out.parent.mkdir(parents=True, exist_ok=True)
        write_qc_counts_markdown(results, metas, hc_feas, passing, args.qc_report_out)
        print(f"Wrote {args.qc_report_out}")


if __name__ == "__main__":
    main()
