#!/usr/bin/env python3
"""
Round 8 (R8-3 Confirmatory): Pre-Registered Unseen 4th-Dataset Confirmation of
Parameter-Free Orthogonal Fusion (`OrthoFusion_scGPT_S5`) on
Feng et al. (Cell Genomics 2025, PMID 41330380, Figshare 27989294 TargetedScreen, MIT License;
individual human iPSC donor cell lines `jejf_2` and `zapk_3`).

Strictly enforces Addendum R8-A Post-QC Eligibility Gate BEFORE running any model:
  - Both lines >= 100 reproducible target genes, >= 80 shared reproducible target genes,
    and >= 2,000 shared readout genes. If gate fails, exits cleanly without relaxing thresholds.

Evaluates pre-registered hypotheses (reports/ROUND8_PREREGISTRATION.md Section 3.2):
  - Formula (ZERO tunable weights):
      y_hat_fusion = b_c + P_S(y_hat_scGPT - b_c) + P_perp(y_hat_S5 - b_c),  S = span{1, b_c, a_c}.
  - H4: delta_r (OrthoFusion - scGPT_Probe) <= 2/10 splits CI < 0 AND v5 (OrthoFusion - S5) <= 2/10 splits CI < 0.
  - H5: disc_raw (OrthoFusion - scGPT_Probe) >= 7/10 splits CI > 0 AND (OrthoFusion - S5) >= 7/10 splits CI > 0.
  - Also records secondary consistency checks on H1/H2/H3 (S5 vs AE) and R8-2 PC1 decomposition on Dataset 4.
"""

import gc
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import Ridge

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compute_row_bootstrap_ci import row_scores_v5, stratified_row_bootstrap_ci
from prepare_feng2025_cache import load_feng_bundle, resolve_feng_dir
from replogle_external_validation import make_strict_union_split, run_ahlmann_eltze_bilinear
from round6_t2_discrim_bootstrap_ci import compute_row_metrics_for_split
from round6_t2_multi_metric_bench import evaluate_5metrics
from round8_r82_and_r83_exploratory import (
    compute_r82_decomposition,
    ortho_fusion_predict,
    pearson_1d,
    resolve_scgpt_dir,
)
from xatlas_panhuman_solver import REPORTS_DIR, SEEDS_10, PanPerturbBridgeModel


def main():
    t0 = time.time()
    feng_dir = resolve_feng_dir()
    cache_sum_path = feng_dir / "cache" / "feng2025_cache_summary.json"
    if not cache_sum_path.exists():
        raise FileNotFoundError(f"Missing Feng 2025 cache summary: {cache_sum_path}. Run prepare_feng2025_cache.py first.")

    cache_sum = json.load(open(cache_sum_path))
    gate = cache_sum["eligibility_gate"]
    if not gate["passed"]:
        print("[HARD STOP] Feng 2025 failed pre-registered Eligibility Gate. Setting confirmation to PENDING.", flush=True)
        out = {
            "dataset": cache_sum["dataset"],
            "status": "PENDING_FAILED_ELIGIBILITY_GATE",
            "eligibility_gate": gate,
            "cache_summary": cache_sum,
        }
        with open(REPORTS_DIR / "r83_confirmatory_feng2025_summary.json", "w") as f:
            json.dump(out, f, indent=2)
        return

    print("=== R8-3 PRE-REGISTERED CONFIRMATORY TEST ON UNSEEN FENG 2025 (jejf_2 + zapk_3) ===", flush=True)
    train_df, all_y, shared_genes, sh0, sh1, g2i_map, sidx = load_feng_bundle()

    wdir = resolve_scgpt_dir()
    vocab = json.load(open(wdir / "scgpt_vocab.json"))
    sd = torch.load(wdir / "scgpt_gh_repo_original_model.bin", map_location="cpu", weights_only=False)
    emb_mat = sd["encoder.embedding.weight"].float().numpy()
    pad_id = vocab["<pad>"]
    del sd

    model = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)
    eval_models = ["Constant_Mean", "AE", "S0", "S5", "S7", "scGPT_FrozenEmb_Ridge", "OrthoFusion_scGPT_S5"]

    rows = []
    ci_rows = []
    split_stats = []
    vocab_covs = []

    for seed in SEEDS_10:
        ts = time.time()
        val_df, tr_df, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
        df_hct = tr_df[tr_df["source"] == "HCT116"]
        df_hek = tr_df[tr_df["source"] == "HEK293T"]
        q_hct = val_df[val_df["source"] == "HCT116"]
        q_hek = val_df[val_df["source"] == "HEK293T"]

        split_stats.append({
            "seed": seed,
            "n_val_jejf2": int(len(q_hct)),
            "n_val_zapk3": int(len(q_hek)),
            "n_union_val_targets": int(val_df["perturbation"].nunique()),
            "n_tr_jejf2_rep": int((df_hct["reproducible"] == 1).sum()),
            "n_tr_zapk3_rep": int((df_hek["reproducible"] == 1).sum()),
        })

        p_const, p_ae, p_scgpt = {}, {}, {}
        cov_s = {}
        for src, df_s, q_s in [("HCT116", df_hct, q_hct), ("HEK293T", df_hek, q_hek)]:
            g2i = g2i_map[src]
            own_rep = df_s[df_s["reproducible"] == 1]
            mu = np.mean([all_y[rid] for rid in own_rep["id"].tolist()], axis=0).astype(np.float64)
            for rid in q_s["id"].tolist():
                p_const[rid] = mu

            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            q_idx = [g2i[p] for p in q_s["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_s["id"]], dtype=np.float32)
            arr_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_s["id"].tolist()):
                p_ae[rid] = arr_ae[i]

            tr_syms = own_rep["gene_symbol"].astype(str).tolist()
            q_syms = q_s["gene_symbol"].astype(str).tolist()
            X_tr = np.array([emb_mat[vocab.get(s, pad_id)] for s in tr_syms], dtype=np.float32)
            X_q = np.array([emb_mat[vocab.get(s, pad_id)] for s in q_syms], dtype=np.float32)
            Y_tr_all = np.array([all_y[rid] for rid in own_rep["id"]], dtype=np.float32)
            p_arr = Ridge(alpha=10.0).fit(X_tr, Y_tr_all).predict(X_q)
            for i, rid in enumerate(q_s["id"].tolist()):
                p_scgpt[rid] = p_arr[i].astype(np.float64)
            cov_s[src] = round(sum(1 for s in q_syms if s in vocab) / max(len(q_syms), 1), 4)
        vocab_covs.append({"seed": seed, "jejf_2": cov_s["HCT116"], "zapk_3": cov_s["HEK293T"]})

        r0 = model.predict_for_source_all_stages("HCT116", df_hct, df_hek, q_hct, all_y, bases, axes, ("S0", "S5", "S7"))
        r1 = model.predict_for_source_all_stages("HEK293T", df_hek, df_hct, q_hek, all_y, bases, axes, ("S0", "S5", "S7"))
        p_s5 = {**r0["S5"], **r1["S5"]}
        p_fusion = ortho_fusion_predict(p_scgpt, p_s5, val_df, bases, axes)

        preds_map = {
            "Constant_Mean": p_const,
            "AE": p_ae,
            "S0": {**r0["S0"], **r1["S0"]},
            "S5": p_s5,
            "S7": {**r0["S7"], **r1["S7"]},
            "scGPT_FrozenEmb_Ridge": p_scgpt,
            "OrthoFusion_scGPT_S5": p_fusion,
        }

        rm = compute_row_metrics_for_split(preds_map, val_df, all_y, bases, axes, sidx, g2i_map)
        v5_map = {m: row_scores_v5(preds_map[m], val_df, all_y, bases, axes, sidx, g2i_map) for m in eval_models}
        row_delta, row_off, z_pack = compute_r82_decomposition(preds_map, tr_df, val_df, all_y, sidx, g2i_map)

        max_v5_diff = float(np.max(np.abs(v5_map["OrthoFusion_scGPT_S5"] - v5_map["S5"])))
        assert max_v5_diff < 1e-10, f"Orthogonal fusion v5 parity violated: {max_v5_diff}"

        src_arr = val_df["source"].values
        i0, i1 = np.where(src_arr == "HCT116")[0], np.where(src_arr == "HEK293T")[0]

        for m in eval_models:
            ev = evaluate_5metrics(preds_map[m], val_df, all_y, bases, axes, sidx, g2i_map)
            axr = float(0.5 * (
                pearson_1d(z_pack["z_true"]["HCT116"], z_pack["z_pred"][m]["HCT116"])
                + pearson_1d(z_pack["z_true"]["HEK293T"], z_pack["z_pred"][m]["HEK293T"])
            ))
            offr = float(0.5 * (row_off[m][i0].mean() + row_off[m][i1].mean()))
            rows.append({
                "dataset": "feng2025", "status": "PREREGISTERED_CONFIRMATORY", "seed": seed, "method": m,
                **{k: round(v, 6) for k, v in ev.items()},
                "disc_raw": round(float(0.5 * (rm[m]["disc_raw"][i0].mean() + rm[m]["disc_raw"][i1].mean())), 6),
                "disc_v5": round(float(0.5 * (rm[m]["disc_v5"][i0].mean() + rm[m]["disc_v5"][i1].mean())), 6),
                "axis_r": round(axr, 6),
                "off_axis_r": round(offr, 6),
            })

        comps = [
            ("OrthoFusion_vs_scGPT_Probe", "OrthoFusion_scGPT_S5", "scGPT_FrozenEmb_Ridge"),
            ("OrthoFusion_vs_S5", "OrthoFusion_scGPT_S5", "S5"),
            ("OrthoFusion_vs_AE", "OrthoFusion_scGPT_S5", "AE"),
            ("S5_vs_AE", "S5", "AE"),
            ("S5_vs_scGPT_Probe", "S5", "scGPT_FrozenEmb_Ridge"),
            ("scGPT_Probe_vs_AE", "scGPT_FrozenEmb_Ridge", "AE"),
        ]
        for cname, ma, mb in comps:
            for met, a_v, b_v in [
                ("v5", v5_map[ma], v5_map[mb]),
                ("delta_r", rm[ma]["delta_r"], rm[mb]["delta_r"]),
                ("disc_raw", rm[ma]["disc_raw"], rm[mb]["disc_raw"]),
                ("disc_v5", rm[ma]["disc_v5"], rm[mb]["disc_v5"]),
                ("off_axis_r", row_off[ma], row_off[mb]),
            ]:
                ci = stratified_row_bootstrap_ci(a_v, b_v, src_arr, n_boot=4000, seed=12345)
                ci_rows.append({
                    "dataset": "feng2025", "status": "PREREGISTERED_CONFIRMATORY",
                    "seed": seed, "metric": met, "comparison": cname,
                    **{k: round(v, 6) if isinstance(v, float) else v for k, v in ci.items()},
                })

        print(f"  [feng2025] seed={seed} done in {time.time()-ts:.1f}s (val={len(q_hct)}+{len(q_hek)})", flush=True)

    df_m = pd.DataFrame(rows)
    df_c = pd.DataFrame(ci_rows)
    df_m.to_csv(REPORTS_DIR / "r83_confirmatory_feng2025_metrics.csv", index=False)
    df_c.to_csv(REPORTS_DIR / "r83_confirmatory_feng2025_row_ci.csv", index=False)

    summary = {
        "dataset": cache_sum["dataset"],
        "status": "PREREGISTERED_CONFIRMATORY",
        "formula": "y_hat = b + P_S(y_hat_scGPT - b) + P_perp(y_hat_S5 - b), S = span{1, b, a}, ZERO tunable weights",
        "shared_readout_genes": len(shared_genes),
        "n_splits": len(SEEDS_10),
        "split_stats": split_stats,
        "scgpt_vocab_coverage": vocab_covs,
        "methods": {},
        "comparisons": {},
    }
    for m, g in df_m.groupby("method"):
        summary["methods"][m] = {
            col: {"mean": round(float(g[col].mean()), 6), "std": round(float(g[col].std(ddof=1)), 6)}
            for col in ["v5", "delta_r", "de50_r", "de50_sign", "disc_raw", "disc_v5", "axis_r", "off_axis_r"]
        }
    for (met, comp), g in df_c.groupby(["metric", "comparison"]):
        summary["comparisons"].setdefault(met, {})[comp] = {
            "mean_delta": round(float(g["delta"].mean()), 6),
            "std_delta": round(float(g["delta"].std(ddof=1)), 6),
            "mean_ci95_lo": round(float(g["ci95_lo"].mean()), 6),
            "mean_ci95_hi": round(float(g["ci95_hi"].mean()), 6),
            "splits_ci_gt_zero": f"{int(g['ci_gt_zero'].sum())}/{len(g)}",
            "splits_ci_lt_zero": f"{int(g['ci_lt_zero'].sum())}/{len(g)}",
        }

    # Pre-registered H4 & H5 adjudication
    c_dr = summary["comparisons"]["delta_r"]["OrthoFusion_vs_scGPT_Probe"]
    c_v5 = summary["comparisons"]["v5"]["OrthoFusion_vs_S5"]
    c_d1 = summary["comparisons"]["disc_raw"]["OrthoFusion_vs_scGPT_Probe"]
    c_d2 = summary["comparisons"]["disc_raw"]["OrthoFusion_vs_S5"]

    h4_dr_neg = int(c_dr["splits_ci_lt_zero"].split("/")[0])
    h4_v5_neg = int(c_v5["splits_ci_lt_zero"].split("/")[0])
    h5_sc_pos = int(c_d1["splits_ci_gt_zero"].split("/")[0])
    h5_s5_pos = int(c_d2["splits_ci_gt_zero"].split("/")[0])

    h4_pass = (h4_dr_neg <= 2) and (h4_v5_neg <= 2)
    h5_pass = (h5_sc_pos >= 7) and (c_d1["mean_delta"] > 0) and (h5_s5_pos >= 7) and (c_d2["mean_delta"] > 0)

    summary["preregistered_adjudication"] = {
        "H4_non_inferiority": {
            "criterion": "delta_r(Fusion - scGPT) <= 2/10 CI < 0 AND v5(Fusion - S5) <= 2/10 CI < 0",
            "delta_r_vs_scGPT": c_dr,
            "v5_vs_S5": c_v5,
            "passed": h4_pass,
            "verdict": "确认 (CONFIRMED)" if h4_pass else "未确认 (NOT CONFIRMED)",
        },
        "H5_disc_raw_synergy": {
            "criterion": "disc_raw(Fusion - scGPT) >= 7/10 CI > 0 AND disc_raw(Fusion - S5) >= 7/10 CI > 0",
            "disc_raw_vs_scGPT": c_d1,
            "disc_raw_vs_S5": c_d2,
            "passed": h5_pass,
            "verdict": "确认 (CONFIRMED)" if h5_pass else "未确认 (NOT CONFIRMED)",
        },
    }
    with open(REPORTS_DIR / "r83_confirmatory_feng2025_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== R8-3 CONFIRMATORY VERDICT ON UNSEEN FENG 2025 ===", flush=True)
    for m in ["AE", "S5", "scGPT_FrozenEmb_Ridge", "OrthoFusion_scGPT_S5"]:
        mm = summary["methods"][m]
        print(f"  {m:24s}: delta_r={mm['delta_r']['mean']:.4f}, v5={mm['v5']['mean']:.4f}, disc_raw={mm['disc_raw']['mean']:.4f}, axis_r={mm['axis_r']['mean']:.4f}, off_axis_r={mm['off_axis_r']['mean']:.4f}", flush=True)
    print(f"  H4: {summary['preregistered_adjudication']['H4_non_inferiority']['verdict']}", flush=True)
    print(f"  H5: {summary['preregistered_adjudication']['H5_disc_raw_synergy']['verdict']}", flush=True)
    print(f"Done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
