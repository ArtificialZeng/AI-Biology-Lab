#!/usr/bin/env python3
"""
T2 Supplemental Row-Level Stratified Paired Bootstrap 95% CIs (4,000 resamples)
for Perturbation Discrimination (raw Delta space `disc_raw` AND v5-projected space `disc`)
and Standard Delta-Expression Pearson (`delta_r`) comparing S5 vs AE and S7 vs AE
on BOTH xatlas (10 dev splits) and Replogle 2022 (10 strict union splits).

Matches `replogle_rowlevel.py` definitions line-for-line.
"""

import gc
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compute_row_bootstrap_ci import stratified_row_bootstrap_ci
from replogle_external_validation import load_replogle_bundle, make_strict_union_split, run_ahlmann_eltze_bilinear
from round6_t2_multi_metric_bench import load_xatlas_bundle, make_xatlas_split
from xatlas_panhuman_solver import REPORTS_DIR, SEEDS_10, PanPerturbBridgeModel, _project_out


def compute_row_metrics_for_split(preds_by_method, val_df, all_y, bases, axes, score_indices, gene2idx):
    """Return {method: {metric: np.ndarray(len(val_df))}} for delta_r, disc_raw, disc_v5."""
    methods = list(preds_by_method.keys())
    n = len(val_df)
    out = {m: {"delta_r": np.empty(n), "disc_raw": np.empty(n), "disc_v5": np.empty(n)} for m in methods}

    for src in ["HCT116", "HEK293T"]:
        mask = (val_df["source"].values == src)
        pos_idx = np.where(mask)[0]
        sub = val_df.iloc[pos_idx].reset_index(drop=True)
        m_len = len(sub)
        base, ax, s_idx = bases[src], axes[src], score_indices[f"{src}_indices"]

        Q_full, Rm = np.linalg.qr(np.stack([np.ones(len(s_idx)), base[s_idx], ax[s_idx]], axis=1))
        Q_full = Q_full[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())]

        # Common-projection true matrices (exact match to replogle_rowlevel.py lines 66-69)
        Tp = np.stack([_project_out(all_y[r].astype(np.float64)[s_idx] - base[s_idx], Q_full) for r in sub["id"]])
        Tp_c = Tp - Tp.mean(axis=1, keepdims=True)
        Tp_c /= np.linalg.norm(Tp_c, axis=1, keepdims=True) + 1e-12

        Tr = np.stack([all_y[r].astype(np.float64)[s_idx] - base[s_idx] for r in sub["id"]])
        Tr_c = Tr - Tr.mean(axis=1, keepdims=True)
        Tr_c /= np.linalg.norm(Tr_c, axis=1, keepdims=True) + 1e-12

        tids = [gene2idx[src].get(p) for p in sub["perturbation"]]
        for m in methods:
            preds = preds_by_method[m]
            Pp = np.stack([_project_out(preds[r].astype(np.float64)[s_idx] - base[s_idx], Q_full) for r in sub["id"]])
            Pp_c = Pp - Pp.mean(axis=1, keepdims=True)
            Pp_c /= np.linalg.norm(Pp_c, axis=1, keepdims=True) + 1e-12
            S_v5 = Pp_c @ Tp_c.T
            d_v5 = np.diag(S_v5)
            disc_v5 = 1.0 - (S_v5 > d_v5[:, None]).sum(axis=1) / max(m_len - 1, 1)

            Pr = np.stack([preds[r].astype(np.float64)[s_idx] - base[s_idx] for r in sub["id"]])
            Pr_c = Pr - Pr.mean(axis=1, keepdims=True)
            Pr_c /= np.linalg.norm(Pr_c, axis=1, keepdims=True) + 1e-12
            S_raw = Pr_c @ Tr_c.T
            d_raw = np.diag(S_raw)
            disc_raw = 1.0 - (S_raw > d_raw[:, None]).sum(axis=1) / max(m_len - 1, 1)

            dr_arr = np.empty(m_len, dtype=np.float64)
            for i, r in enumerate(sub["id"]):
                t = tids[i]
                idx = s_idx[s_idx != t] if t is not None else s_idx
                a = all_y[r].astype(np.float64)[idx] - base[idx]
                b = preds[r].astype(np.float64)[idx] - base[idx]
                a = a - a.mean(); b = b - b.mean()
                denom = np.linalg.norm(a) * np.linalg.norm(b)
                dr_arr[i] = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0)) if denom > 1e-12 else 0.0

            out[m]["delta_r"][pos_idx] = dr_arr
            out[m]["disc_raw"][pos_idx] = disc_raw
            out[m]["disc_v5"][pos_idx] = disc_v5
    return out


def run_dataset(dataset: str):
    t0 = time.time()
    print(f"=== T2 Discrim/Delta Row-Level 95% CI on {dataset} (10 splits) ===", flush=True)
    if dataset == "xatlas":
        train_df, all_y, bases_pub, axes_pub, sidx, sh0, sh1, g2i_map = load_xatlas_bundle()
    else:
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_replogle_bundle()
    model = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)

    records = []
    for seed in SEEDS_10:
        ts = time.time()
        if dataset == "xatlas":
            val_df, vr_set, vt_line = make_xatlas_split(train_df, g2i_map, seed)
            bases, axes = bases_pub, axes_pub
        else:
            val_df, df_tr, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)

        preds = {"S0": {}, "S5": {}, "S7": {}, "AE": {}}
        for src in ["HCT116", "HEK293T"]:
            oth = "HEK293T" if src == "HCT116" else "HCT116"
            g2i = g2i_map[src]
            if dataset == "xatlas":
                df_own = train_df[(train_df["source"] == src) & (~train_df["id"].isin(vr_set))]
                df_oth = train_df[(train_df["source"] == oth) & (~train_df["perturbation"].isin(vt_line[src]))]
            else:
                df_own = df_tr[df_tr["source"] == src]
                df_oth = df_tr[df_tr["source"] == oth]
            q_df = val_df[val_df["source"] == src]
            res = model.predict_for_source_all_stages(src, df_own, df_oth, q_df, all_y, bases, axes, ("S0", "S5", "S7"))
            for st in ("S0", "S5", "S7"):
                preds[st].update(res[st])

            own_rep = df_own[df_own["reproducible"] == 1]
            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            q_idx = [g2i[p] for p in q_df["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_own["id"]], dtype=np.float32)
            p_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_df["id"].tolist()):
                preds["AE"][rid] = p_ae[i]

        rm = compute_row_metrics_for_split(preds, val_df, all_y, bases, axes, sidx, g2i_map)
        src_arr = val_df["source"].values
        for metric in ["delta_r", "disc_raw", "disc_v5"]:
            for comp, ka, kb in [("S5_vs_AE", "S5", "AE"), ("S7_vs_AE", "S7", "AE"), ("S5_vs_S0", "S5", "S0"), ("S7_vs_S0", "S7", "S0")]:
                st_ci = stratified_row_bootstrap_ci(rm[ka][metric], rm[kb][metric], src_arr)
                records.append({"dataset": dataset, "seed": seed, "metric": metric, "comparison": comp, **st_ci})
        r5 = [r for r in records if r["seed"] == seed and r["metric"] == "disc_raw" and r["comparison"] == "S5_vs_AE"][0]
        r7 = [r for r in records if r["seed"] == seed and r["metric"] == "disc_raw" and r["comparison"] == "S7_vs_AE"][0]
        print(
            f"  [{dataset}] Seed {seed:5d} | disc_raw S5-AE={r5['delta']:+.4f} [{r5['ci95_lo']:+.4f}, {r5['ci95_hi']:+.4f}] "
            f"| S7-AE={r7['delta']:+.4f} [{r7['ci95_lo']:+.4f}, {r7['ci95_hi']:+.4f}] [{time.time()-ts:.1f}s]",
            flush=True,
        )
        gc.collect()
    return pd.DataFrame(records)


def main():
    df_all = pd.concat([run_dataset("xatlas"), run_dataset("replogle2022")], ignore_index=True)
    df_all.to_csv(REPORTS_DIR / "t2_discrim_delta_row_ci.csv", index=False)
    summary = {}
    for (ds, metric, comp), grp in df_all.groupby(["dataset", "metric", "comparison"]):
        summary.setdefault(ds, {}).setdefault(metric, {})[comp] = {
            "mean_a": round(float(grp["score_a"].mean()), 4),
            "mean_b": round(float(grp["score_b"].mean()), 4),
            "mean_delta": round(float(grp["delta"].mean()), 4),
            "median_ci95_lo": round(float(grp["ci95_lo"].median()), 4),
            "median_ci95_hi": round(float(grp["ci95_hi"].median()), 4),
            "mean_boot_se": round(float(grp["boot_se"].mean()), 4),
            "splits_ci_gt_zero": f"{int(grp['ci_gt_zero'].sum())}/10",
            "splits_ci_lt_zero": f"{int(grp['ci_lt_zero'].sum())}/10",
        }
    with open(REPORTS_DIR / "t2_discrim_delta_row_ci_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n=== T2 DISCRIMINATION & DELTA-R ROW-LEVEL 95% CI SUMMARY ===")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
