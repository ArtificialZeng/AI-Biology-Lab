#!/usr/bin/env python3
"""
Replogle-only (K562 <-> RPE1) development ablation sweep.

Strictly enforces Rule 3: refuses to run on any non-Replogle cell line.
Used exclusively to benchmark preprocessing representations (`log1p`, `zscore`, `libnorm`),
shrinkage (`0.0` vs `1e-3`), and `mask_self` choices prior to locking `PREREG_A.md`.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

SRC_DIR = Path(__file__).resolve().parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from adjudicate import adjudicate_results
from data_loader import is_replogle_dev_line, list_available_lines, load_line_dataset
from run_eval import RESULTS_DIR, evaluate_datasets


def main() -> None:
    t0 = time.time()
    avail = [ln for ln in list_available_lines() if is_replogle_dev_line(ln)]
    if len(avail) < 2:
        raise RuntimeError(f"Need >= 2 Replogle dev lines in p2/shared/lines/, found: {avail}")

    # Align strictly onto common genes of the two Replogle lines (avoid touching held-out COMMON_GENES.txt)
    g_sets = []
    for ln in avail:
        g_file = SRC_DIR.parents[1] / "shared" / "lines" / ln / "genes.txt"
        g_sets.append([x.strip() for x in g_file.read_text(encoding="utf-8").splitlines() if x.strip()])
    common_replogle_genes = sorted(set(g_sets[0]).intersection(*g_sets[1:]))
    print(f"Replogle dev lines={avail}, common genes={len(common_replogle_genes)}", flush=True)

    configs = [
        {"repr": "log1p", "shrink": 0.0, "mask_self": True, "only_rep": False},
        {"repr": "zscore", "shrink": 0.0, "mask_self": True, "only_rep": False},
        {"repr": "libnorm", "shrink": 0.0, "mask_self": True, "only_rep": False},
        {"repr": "log1p", "shrink": 1e-3, "mask_self": True, "only_rep": False},
        {"repr": "zscore", "shrink": 1e-3, "mask_self": True, "only_rep": False},
        {"repr": "log1p", "shrink": 0.0, "mask_self": True, "only_rep": True},
        {"repr": "zscore", "shrink": 0.0, "mask_self": True, "only_rep": True},
    ]

    sweep_rows: List[Dict[str, Any]] = []
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    for cfg in configs:
        tag = f"dev_{cfg['repr']}_sh{cfg['shrink']}_rep{int(cfg['only_rep'])}"
        datasets = {
            ln: load_line_dataset(ln, representation=cfg["repr"], common_genes=common_replogle_genes)
            for ln in avail
        }
        df_lines, df_boot, summary = evaluate_datasets(
            datasets=datasets,
            shrink=cfg["shrink"],
            mask_self=cfg["mask_self"],
            shared_targets_only=True,
            only_reproducible=cfg["only_rep"],
            n_boot=4000,
            seed=12345,
            require_official_spec=True,
        )
        df_lines.to_csv(RESULTS_DIR / f"line_metrics_{tag}.csv", index=False)
        df_boot.to_csv(RESULTS_DIR / f"bootstrap_ci_{tag}.csv", index=False)
        (RESULTS_DIR / f"summary_{tag}.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8"
        )
        adj = adjudicate_results(df_lines, df_boot, summary, primary_method="Hybrid_ProjShared", min_lines_pass=2)
        (RESULTS_DIR / f"adjudication_{tag}.json").write_text(
            json.dumps(adj, indent=2), encoding="utf-8"
        )
        for m_name, vals in summary["macro_means_by_method"].items():
            sweep_rows.append(
                {
                    "tag": tag,
                    **cfg,
                    "method": m_name,
                    "n_eval_targets": vals["n_eval_targets"],
                    "ls_slope_ols_top50": vals["ls_slope_ols_top50"],
                    "ls_resid_corr_cal_top50": vals["ls_resid_corr_cal_top50"],
                    "cent_disc_raw": vals["cent_disc_raw"],
                    "cent_de50_r": vals["cent_de50_r"],
                    "ctrl_disc_raw": vals["ctrl_disc_raw"],
                    "ctrl_de50_r": vals["ctrl_de50_r"],
                    "mse": vals["mse"],
                }
            )
        print(f"Completed {tag} in {time.time() - t0:.1f}s", flush=True)

    pd.DataFrame(sweep_rows).to_csv(RESULTS_DIR / "dev_replogle_sweep.csv", index=False)
    print(f"Saved dev_replogle_sweep.csv ({len(sweep_rows)} rows)", flush=True)


if __name__ == "__main__":
    main()
