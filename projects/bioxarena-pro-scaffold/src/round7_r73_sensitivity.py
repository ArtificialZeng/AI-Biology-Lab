#!/usr/bin/env python3
"""
Round 7 (R7-3): One-Factor-At-A-Time (+-50%) Hyperparameter Sensitivity Analysis
for Primary Frozen Method S5 on External Datasets (Replogle 2022 & Nadig 2025).

Purpose:
  Demonstrate robustness of primary method S5 (Cross-Line Bridge) across +-50%
  perturbations of every constant used in S5. Strictly diagnostic / evaluation-only:
  NEVER modifies any frozen constant in `src/xatlas_panhuman_solver.py`.

Reports for each perturbed constant (-50%, Frozen 0%, +50%):
  - v5_orth_r (mean +- SD across SEEDS_10, delta vs AE, # splits with row-level 95% CI > 0)
  - disc_raw  (mean +- SD across SEEDS_10, delta vs AE, # splits with row-level 95% CI > 0)
  - delta_r   (mean +- SD across SEEDS_10, delta vs AE)
"""

import argparse
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
from replogle_external_validation import (
    load_replogle_bundle,
    make_strict_union_split,
    run_ahlmann_eltze_bilinear,
)
from round6_t2_discrim_bootstrap_ci import compute_row_metrics_for_split
from xatlas_panhuman_solver import (
    REPORTS_DIR,
    SEEDS_10,
    PanPerturbBridgeModel,
    build_pyramid_kernel,
    build_spectral_pyramid_Gs,
    compute_hybrid_tom_from_G,
    fit_geneslope_coskrr_translator,
    krr_weighted_solve,
    project_rows_5000,
    score_prediction,
)

# Frozen S5 defaults (exact match to PanPerturbBridgeModel @ commit 3609811)
FROZEN_S5_PARAMS = {
    "soft_p": 1.5,              # HybridTOM soft-power exponent
    "krr_alpha": 0.5,           # KRR regularization alpha
    "pan_other_scale": 1.0,     # Multiplier on PAN_OTHER_WEIGHTS {0:0.6,1:0.8,2:0.8,5:0.6,10:0.5}
    "deflation_scale": 1.0,     # Multiplier on deflated PC weights {1:0.8,2:0.5,5:0.3,10:0.2} (d0=1.0)
    "nr_weight_scale": 0.25,    # Non-reproducible anchor weight base multiplier
    "extra_trans_weight": 0.6,  # Cross-line translated extra target sample weight
    "trans_skip_w": 0.8,        # GeneSlope-CosKRR translator skip weight
    "trans_alpha": 2.0,         # GeneSlope-CosKRR translator ridge alpha
    "trans_pk": 1.5,            # GeneSlope-CosKRR translator signed-power exponent
    "bridge_w": 0.9,            # Symmetric cross-line bridge fusion weight
}

BASE_PC_WEIGHTS = {0: 1.0, 1: 0.8, 2: 0.5, 5: 0.3, 10: 0.2}
BASE_PAN_OTHER_WEIGHTS = {0: 0.6, 1: 0.8, 2: 0.8, 5: 0.6, 10: 0.5}
PC_LEVELS = (0, 1, 2, 5, 10)


def build_sensitivity_configs():
    """Return list of (config_name, param_name, factor_label, params_dict)."""
    configs = [("S5_Frozen", "NONE", "0%", dict(FROZEN_S5_PARAMS))]
    for k, v in FROZEN_S5_PARAMS.items():
        for lbl, mult in [("-50%", 0.5), ("+50%", 1.5)]:
            p = dict(FROZEN_S5_PARAMS)
            p[k] = round(v * mult, 6)
            configs.append((f"{k}_{lbl}", k, lbl, p))
    return configs


def predict_s5_parametric(
    source: str,
    df_own: pd.DataFrame,
    df_other: pd.DataFrame,
    query_df: pd.DataFrame,
    all_train_y: dict,
    baselines: dict,
    axes: dict,
    gene2idx: dict,
    hct_sh_idx: np.ndarray,
    hek_sh_idx: np.ndarray,
    score_indices: dict,
    Gs_own: dict,
    Gs_oth: dict,
    params: dict,
    tom_cache: dict,
) -> dict:
    """Exact reproduction of PanPerturbBridgeModel S5 with parameterized constants + exact TOM kernel cache."""
    other_src = "HEK293T" if source == "HCT116" else "HCT116"
    g2i = gene2idx[source]
    g2i_o = gene2idx[other_src]
    s_idx = score_indices[f"{source}_indices"]
    s_idx_o = score_indices[f"{other_src}_indices"]
    own_sh = hct_sh_idx if source == "HCT116" else hek_sh_idx
    oth_sh = hek_sh_idx if source == "HCT116" else hct_sh_idx

    base_own, ax_own = baselines[source], axes[source]
    base_oth, ax_oth = baselines[other_src], axes[other_src]

    query_perts = query_df["perturbation"].tolist()
    query_idx = [g2i[p] for p in query_perts]
    query_rids = query_df["id"].tolist()

    own_rep = df_own[df_own["reproducible"] == 1]
    oth_rep = df_other[df_other["reproducible"] == 1]
    tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
    tr_perts = tr_sub["perturbation"].tolist()
    tr_idx = [g2i[p] for p in tr_perts]
    tr_Y_raw = np.array([all_train_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
    tr_Y_p = project_rows_5000(tr_Y_raw, base_own, ax_own)

    # Parameterized PC weights and pan-human weights
    pc_weights = {
        k: (1.0 if k == 0 else BASE_PC_WEIGHTS[k] * params["deflation_scale"])
        for k in PC_LEVELS
    }
    pan_scale = params["pan_other_scale"]
    soft_p = params["soft_p"]

    def get_cached_pyramid_kernel(branch_tag, Gs_base, Gs_partner, sh_base, sh_partner, tr_i, q_i):
        K_tr = np.zeros((len(tr_i), len(tr_i)), dtype=np.float64)
        K_val = np.zeros((len(q_i), len(tr_i)), dtype=np.float64)
        for k, w in pc_weights.items():
            if w <= 0:
                continue
            ckey = (source, branch_tag, k, pan_scale, soft_p)
            if ckey not in tom_cache:
                wo = BASE_PAN_OTHER_WEIGHTS[k] * pan_scale
                M = Gs_base[k].copy()
                M[np.ix_(sh_base, sh_base)] += wo * Gs_partner[k][np.ix_(sh_partner, sh_partner)]
                tom_cache[ckey] = compute_hybrid_tom_from_G(M, tr_i, q_i, soft_p=soft_p)
            kt, kv = tom_cache[ckey]
            K_tr += w * kt
            K_val += w * kv
        return K_tr, K_val

    # Non-reproducible anchors
    nr_sub_all = df_own[(df_own["reproducible"] == 0) & (df_own["perturbation"].isin(g2i.keys()))]
    if len(nr_sub_all) > 0:
        nr_Y_all_p = project_rows_5000(
            np.array([all_train_y[rid] for rid in nr_sub_all["id"]], dtype=np.float32),
            base_own,
            ax_own,
        )
        z_eff_all = np.linalg.norm(nr_Y_all_p[:, s_idx], axis=1) * np.sqrt(
            nr_sub_all["n_cells"].values.astype(np.float64)
        )
        med_z = max(float(np.median(z_eff_all)), 1e-8)
        w_nr_all = params["nr_weight_scale"] * np.clip(
            (z_eff_all - 0.9 * med_z) / (0.3 * med_z), 0.1, 1.5
        )
    else:
        nr_Y_all_p = np.zeros((0, 5000), dtype=np.float64)
        w_nr_all = np.zeros((0,), dtype=np.float64)

    # Cross-line translator
    own_p2id = dict(zip(own_rep["perturbation"], own_rep["id"]))
    oth_p2id = dict(zip(oth_rep["perturbation"], oth_rep["id"]))
    paired = sorted(list(set(own_p2id.keys()) & set(oth_p2id.keys())))
    Y_p_own = project_rows_5000(
        np.array([all_train_y[own_p2id[p]] for p in paired], dtype=np.float32), base_own, ax_own
    )
    Y_p_oth = project_rows_5000(
        np.array([all_train_y[oth_p2id[p]] for p in paired], dtype=np.float32), base_oth, ax_oth
    )
    trans_o2s = fit_geneslope_coskrr_translator(
        Y_p_oth,
        Y_p_own,
        oth_sh,
        own_sh,
        skip_w=params["trans_skip_w"],
        alpha=params["trans_alpha"],
        p_k=params["trans_pk"],
    )
    trans_s2o = fit_geneslope_coskrr_translator(
        Y_p_own,
        Y_p_oth,
        own_sh,
        oth_sh,
        skip_w=params["trans_skip_w"],
        alpha=params["trans_alpha"],
        p_k=params["trans_pk"],
    )

    # Model 1: own cell line
    tr_set = set(tr_perts)
    extra_df = oth_rep[(oth_rep["perturbation"].isin(g2i.keys())) & (~oth_rep["perturbation"].isin(tr_set))]
    extra_perts = extra_df["perturbation"].tolist()
    extra_set = set(extra_perts)

    all_idx = list(tr_idx)
    y_blocks = [tr_Y_p]
    sw_blocks = [np.ones(len(tr_idx), dtype=np.float64)]

    if len(extra_perts) > 0:
        Y_ext_trans = trans_o2s(
            project_rows_5000(
                np.array([all_train_y[rid] for rid in extra_df["id"]], dtype=np.float32),
                base_oth,
                ax_oth,
            )
        )
        all_idx += [g2i[p] for p in extra_perts]
        y_blocks.append(Y_ext_trans)
        sw_blocks.append(np.full(len(extra_perts), params["extra_trans_weight"], dtype=np.float64))

    if len(nr_sub_all) > 0:
        keep_nr_mask = np.array([p not in extra_set for p in nr_sub_all["perturbation"].tolist()], dtype=bool)
        if np.any(keep_nr_mask):
            nr_perts_kept = nr_sub_all["perturbation"].values[keep_nr_mask].tolist()
            all_idx += [g2i[p] for p in nr_perts_kept]
            y_blocks.append(nr_Y_all_p[keep_nr_mask])
            sw_blocks.append(w_nr_all[keep_nr_mask])

    all_Y = np.vstack(y_blocks)
    sw = np.concatenate(sw_blocks)
    K_tr, K_val = get_cached_pyramid_kernel("own", Gs_own, Gs_oth, own_sh, oth_sh, all_idx, query_idx)
    pred_own = krr_weighted_solve(K_tr, K_val, all_Y, alpha=params["krr_alpha"], sample_weight=sw)

    # Model 2: other cell line bridge
    oth_meas = oth_rep[oth_rep["perturbation"].isin(g2i_o.keys())]
    oth_tr_perts = oth_meas["perturbation"].tolist()
    oth_tr_idx = [g2i_o[p] for p in oth_tr_perts]
    oth_tr_Y = project_rows_5000(
        np.array([all_train_y[rid] for rid in oth_meas["id"]], dtype=np.float32), base_oth, ax_oth
    )
    oth_sw = np.ones(len(oth_tr_idx), dtype=np.float64)
    oth_set = set(oth_tr_perts)

    ext_o_df = own_rep[(own_rep["perturbation"].isin(g2i_o.keys())) & (~own_rep["perturbation"].isin(oth_set))]
    ext_o_perts = ext_o_df["perturbation"].tolist()
    if len(ext_o_perts) > 0:
        Y_ext_o = trans_s2o(
            project_rows_5000(
                np.array([all_train_y[rid] for rid in ext_o_df["id"]], dtype=np.float32),
                base_own,
                ax_own,
            )
        )
        oth_tr_idx += [g2i_o[p] for p in ext_o_perts]
        oth_tr_Y = np.vstack([oth_tr_Y, Y_ext_o])
        oth_sw = np.concatenate([oth_sw, np.full(len(ext_o_perts), params["extra_trans_weight"])])
        oth_set |= set(ext_o_perts)

    nr_o_sub = df_other[
        (df_other["reproducible"] == 0)
        & (df_other["perturbation"].isin(g2i_o.keys()))
        & (~df_other["perturbation"].isin(oth_set))
    ]
    if len(nr_o_sub) > 0:
        nr_o_Y = project_rows_5000(
            np.array([all_train_y[rid] for rid in nr_o_sub["id"]], dtype=np.float32), base_oth, ax_oth
        )
        z_o = np.linalg.norm(nr_o_Y[:, s_idx_o], axis=1) * np.sqrt(nr_o_sub["n_cells"].values.astype(np.float64))
        med_zo = max(float(np.median(z_o)), 1e-8)
        w_nro = params["nr_weight_scale"] * np.clip((z_o - 0.9 * med_zo) / (0.3 * med_zo), 0.1, 1.5)
        oth_tr_idx += [g2i_o[p] for p in nr_o_sub["perturbation"].tolist()]
        oth_tr_Y = np.vstack([oth_tr_Y, nr_o_Y])
        oth_sw = np.concatenate([oth_sw, w_nro])

    q_in_o = [i for i, p in enumerate(query_perts) if p in g2i_o]
    q_oi = [g2i_o[query_perts[i]] for i in q_in_o]
    out_s5 = pred_own.copy()
    if len(q_oi) > 0:
        K_tr_o, K_val_o = get_cached_pyramid_kernel("oth", Gs_oth, Gs_own, oth_sh, own_sh, oth_tr_idx, q_oi)
        pred_in_other = krr_weighted_solve(K_tr_o, K_val_o, oth_tr_Y, alpha=params["krr_alpha"], sample_weight=oth_sw)
        p_br_s5 = trans_o2s(pred_in_other)
        for j, qi in enumerate(q_in_o):
            out_s5[qi] = out_s5[qi] + params["bridge_w"] * p_br_s5[j]

    return {rid: base_own + out_s5[i] for i, rid in enumerate(query_rids)}


def load_dataset_bundle(dataset: str):
    if dataset == "replogle":
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_replogle_bundle()
        return train_df, all_y, sh0, sh1, g2i_map, sidx
    elif dataset == "nadig":
        from prepare_nadig2025_cache import load_nadig_bundle
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_nadig_bundle()
        return train_df, all_y, sh0, sh1, g2i_map, sidx
    else:
        raise ValueError(f"Unsupported dataset for R7-3 sensitivity: {dataset}")


def run_sensitivity(dataset: str):
    t0 = time.time()
    print(f"=== R7-3 S5 +-50% Hyperparameter Sensitivity on {dataset} (10 splits) ===", flush=True)
    train_df, all_y, sh0, sh1, g2i_map, sidx = load_dataset_bundle(dataset)
    configs = build_sensitivity_configs()

    # Verify S5_Frozen matches exact PanPerturbBridgeModel S5 on first split
    ref_model = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)
    records = []

    for split_i, seed in enumerate(SEEDS_10):
        ts = time.time()
        val_df, tr_df, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
        df_hct = tr_df[tr_df["source"] == "HCT116"]
        df_hek = tr_df[tr_df["source"] == "HEK293T"]
        q_hct = val_df[val_df["source"] == "HCT116"]
        q_hek = val_df[val_df["source"] == "HEK293T"]
        sources = val_df["source"].values

        # Precompute AE baseline and spectral pyramids once per split
        p_ae = {}
        for src, df_s, q_s in [("HCT116", df_hct, q_hct), ("HEK293T", df_hek, q_hek)]:
            g2i = g2i_map[src]
            own_rep = df_s[df_s["reproducible"] == 1]
            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            q_idx = [g2i[p] for p in q_s["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_s["id"]], dtype=np.float32)
            arr_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_s["id"].tolist()):
                p_ae[rid] = arr_ae[i]
        Gs_hct = build_spectral_pyramid_Gs(df_hct, all_y, g2i_map["HCT116"], bases["HCT116"], PC_LEVELS)
        Gs_hek = build_spectral_pyramid_Gs(df_hek, all_y, g2i_map["HEK293T"], bases["HEK293T"], PC_LEVELS)

        # Row-level v5 for AE
        v5_ae = np.array([
            score_prediction(
                all_y[r["id"]].astype(np.float64),
                p_ae[r["id"]].astype(np.float64),
                bases[r["source"]],
                axes[r["source"]],
                sidx[f"{r['source']}_indices"],
                g2i_map[r["source"]].get(r["perturbation"]),
            )
            for _, r in val_df.iterrows()
        ])
        ae_row_metrics = compute_row_metrics_for_split({"AE": p_ae}, val_df, all_y, bases, axes, sidx, g2i_map)["AE"]

        tom_cache = {}
        for cfg_name, pname, flbl, params in configs:
            p0 = predict_s5_parametric(
                "HCT116", df_hct, df_hek, q_hct, all_y, bases, axes, g2i_map, sh0, sh1, sidx, Gs_hct, Gs_hek, params, tom_cache
            )
            p1 = predict_s5_parametric(
                "HEK293T", df_hek, df_hct, q_hek, all_y, bases, axes, g2i_map, sh0, sh1, sidx, Gs_hek, Gs_hct, params, tom_cache
            )
            preds = {**p0, **p1}

            # Exact bit-parity check against frozen PanPerturbBridgeModel S5 on seed 42
            if split_i == 0 and cfg_name == "S5_Frozen":
                r0 = ref_model.predict_for_source_all_stages(
                    "HCT116", df_hct, df_hek, q_hct, all_y, bases, axes, stages_to_run=("S5",)
                )["S5"]
                r1 = ref_model.predict_for_source_all_stages(
                    "HEK293T", df_hek, df_hct, q_hek, all_y, bases, axes, stages_to_run=("S5",)
                )["S5"]
                ref_p = {**r0, **r1}
                max_diff = max(float(np.max(np.abs(preds[rid] - ref_p[rid]))) for rid in preds)
                print(f"  [PARITY CHECK seed={seed}] S5_Frozen vs PanPerturbBridgeModel S5 max_abs_diff = {max_diff:.3e}", flush=True)
                assert max_diff < 1e-10, f"S5_Frozen parity mismatch: {max_diff}"

            v5_cfg = np.array([
                score_prediction(
                    all_y[r["id"]].astype(np.float64),
                    preds[r["id"]].astype(np.float64),
                    bases[r["source"]],
                    axes[r["source"]],
                    sidx[f"{r['source']}_indices"],
                    g2i_map[r["source"]].get(r["perturbation"]),
                )
                for _, r in val_df.iterrows()
            ])
            cfg_row_metrics = compute_row_metrics_for_split({cfg_name: preds}, val_df, all_y, bases, axes, sidx, g2i_map)[cfg_name]

            ci_v5 = stratified_row_bootstrap_ci(v5_cfg, v5_ae, sources)
            ci_disc = stratified_row_bootstrap_ci(cfg_row_metrics["disc_raw"], ae_row_metrics["disc_raw"], sources)
            ci_dr = stratified_row_bootstrap_ci(cfg_row_metrics["delta_r"], ae_row_metrics["delta_r"], sources)

            records.append({
                "dataset": dataset,
                "seed": seed,
                "config": cfg_name,
                "param": pname,
                "perturbation": flbl,
                "v5": ci_v5["score_a"],
                "v5_ae": ci_v5["score_b"],
                "v5_delta_vs_ae": ci_v5["delta"],
                "v5_ci95_lo": ci_v5["ci95_lo"],
                "v5_ci95_hi": ci_v5["ci95_hi"],
                "v5_ci_gt_zero": ci_v5["ci_gt_zero"],
                "disc_raw": ci_disc["score_a"],
                "disc_raw_ae": ci_disc["score_b"],
                "disc_raw_delta_vs_ae": ci_disc["delta"],
                "disc_raw_ci95_lo": ci_disc["ci95_lo"],
                "disc_raw_ci95_hi": ci_disc["ci95_hi"],
                "disc_raw_ci_gt_zero": ci_disc["ci_gt_zero"],
                "delta_r": ci_dr["score_a"],
                "delta_r_ae": ci_dr["score_b"],
                "delta_r_delta_vs_ae": ci_dr["delta"],
            })

        del Gs_hct, Gs_hek, tom_cache
        gc.collect()
        print(f"  [seed={seed}] 21 configs done in {time.time()-ts:.1f}s", flush=True)

    df = pd.DataFrame(records)
    csv_path = REPORTS_DIR / f"r73_sensitivity_{dataset}.csv"
    df.to_csv(csv_path, index=False)

    # Aggregate summary
    frozen_v5 = float(df[df["config"] == "S5_Frozen"]["v5"].mean())
    frozen_disc = float(df[df["config"] == "S5_Frozen"]["disc_raw"].mean())
    summary = {"dataset": dataset, "n_splits": len(SEEDS_10), "configs": {}}
    for cfg_name, pname, flbl, _ in configs:
        sub = df[df["config"] == cfg_name]
        summary["configs"][cfg_name] = {
            "param": pname,
            "perturbation": flbl,
            "v5_mean": round(float(sub["v5"].mean()), 6),
            "v5_std": round(float(sub["v5"].std(ddof=1)), 6),
            "v5_shift_vs_frozen": round(float(sub["v5"].mean() - frozen_v5), 6),
            "v5_delta_vs_ae_mean": round(float(sub["v5_delta_vs_ae"].mean()), 6),
            "v5_ci_gt_zero_splits": f"{int(sub['v5_ci_gt_zero'].sum())}/{len(sub)}",
            "disc_raw_mean": round(float(sub["disc_raw"].mean()), 6),
            "disc_raw_std": round(float(sub["disc_raw"].std(ddof=1)), 6),
            "disc_raw_shift_vs_frozen": round(float(sub["disc_raw"].mean() - frozen_disc), 6),
            "disc_raw_delta_vs_ae_mean": round(float(sub["disc_raw_delta_vs_ae"].mean()), 6),
            "disc_raw_ci_gt_zero_splits": f"{int(sub['disc_raw_ci_gt_zero'].sum())}/{len(sub)}",
            "delta_r_delta_vs_ae_mean": round(float(sub["delta_r_delta_vs_ae"].mean()), 6),
        }

    json_path = REPORTS_DIR / f"r73_sensitivity_{dataset}_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n=== R7-3 Summary ({dataset}, {time.time()-t0:.1f}s) ===", flush=True)
    for k, v in summary["configs"].items():
        print(
            f"  {k:24s} | v5={v['v5_mean']:.4f} (shift={v['v5_shift_vs_frozen']:+.4f}, vsAE CI>0 {v['v5_ci_gt_zero_splits']}) "
            f"| disc_raw={v['disc_raw_mean']:.4f} (shift={v['disc_raw_shift_vs_frozen']:+.4f}, vsAE CI>0 {v['disc_raw_ci_gt_zero_splits']}) "
            f"| dr-AE={v['delta_r_delta_vs_ae_mean']:+.4f}",
            flush=True,
        )
    print(f"Saved: {csv_path} & {json_path}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["replogle", "nadig"], required=True)
    args = ap.parse_args()
    run_sensitivity(args.dataset)
