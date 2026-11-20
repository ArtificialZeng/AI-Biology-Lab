#!/usr/bin/env python3
"""
T4: Leave-One-Module-Out (LOMO) Ablation on BOTH xatlas (10 dev splits)
and Replogle et al. (Cell 2022, 10 strict union splits).

Uses the EXACT frozen hyperparameters of PanPerturbBridgeModel (`src/xatlas_panhuman_solver.py`
is NOT modified). Starting from Full S7, removes ONE component at a time while keeping
all other components fixed:
  1. LOO_S1_SelfMask       : disable self-knockdown spike zeroing in gene graph
  2. LOO_S1_BaseSub        : disable baseline subtraction in gene graph
  3. LOO_S1_SqrtNCells     : replace sqrt(n_cells) row weighting with uniform weights
  4. LOO_S1_HybridTOM      : replace soft-power p=1.5 HybridTOM with raw cosine p=1.0
  5. LOO_S2_Deflation      : replace d0..d10 spectral deflation pyramid with d0 only
  6. LOO_S3_CrossLineGraph : disable cross-line gene graph fusion (PAN_OTHER_WEIGHTS = 0)
  7. LOO_S4_NonRepAnchors  : drop non-reproducible measured anchors from weighted KRR
  8. LOO_S5_CrossLineBridge: disable cross-line translator extra targets & bidirectional bridge
  9. LOO_S6_Reciprocal     : disable Reciprocal Response Borrowing (rec_w = 0)
 10. LOO_S7_WienerSnap     : disable Stage 7 output SVD Wiener denoise & manifold snapping

Reports exact module contribution `Full_S7 - LOO_variant` with row-level stratified
paired bootstrap 95% CI (4,000 resamples) per split.
"""

import gc
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compute_row_bootstrap_ci import row_scores_v5, stratified_row_bootstrap_ci
from replogle_external_validation import load_replogle_bundle, make_strict_union_split
from round6_t2_multi_metric_bench import load_xatlas_bundle, make_xatlas_split
from xatlas_panhuman_solver import (
    REPORTS_DIR,
    SEEDS_10,
    PanPerturbBridgeModel,
    build_pyramid_kernel,
    compute_relu_reciprocal,
    fit_geneslope_coskrr_translator,
    krr_weighted_solve,
    manifold_snap,
    project_rows_5000,
    svd_wiener_denoise,
)


def build_pyramid_configurable(df_s, all_y, g2i, base, pc_levels=(0, 1, 2, 5, 10),
                               self_mask=True, base_sub=True, sqrt_ncells=True):
    Y = np.array([all_y[rid] for rid in df_s["id"]], dtype=np.float32)
    if self_mask:
        for i, p in enumerate(df_s["perturbation"].tolist()):
            if p in g2i:
                Y[i, g2i[p]] = 0.0
    if base_sub:
        Y = Y - base.astype(np.float32)[None, :]
    if sqrt_ncells:
        w = np.sqrt(df_s["n_cells"].values.astype(np.float32))
        Y = Y * (w / np.mean(w))[:, None]
    max_pc = max(pc_levels)
    if max_pc > 0:
        svd = TruncatedSVD(n_components=max_pc, random_state=42)
        U = svd.fit_transform(Y)
        V = svd.components_
    Gs = {}
    for k in pc_levels:
        Yk = Y if k == 0 else (Y - U[:, :k] @ V[:k, :])
        Ck = Yk.T
        Ck_u = Ck / (np.linalg.norm(Ck, axis=1, keepdims=True) + 1e-12)
        Gs[k] = Ck_u @ Ck_u.T
    return Gs


def run_variant_for_source(source, df_own, df_other, query_df, all_y, bases, axes,
                           gene2idx, hct_sh_idx, hek_sh_idx, score_indices, cfg):
    """Exact PanPerturbBridgeModel S7 logic with one boolean flag toggled off for LOMO ablation."""
    other_src = "HEK293T" if source == "HCT116" else "HCT116"
    g2i, g2i_o = gene2idx[source], gene2idx[other_src]
    s_idx, s_idx_o = score_indices[f"{source}_indices"], score_indices[f"{other_src}_indices"]
    own_sh = hct_sh_idx if source == "HCT116" else hek_sh_idx
    oth_sh = hek_sh_idx if source == "HCT116" else hct_sh_idx
    base_own, ax_own = bases[source], axes[source]
    base_oth, ax_oth = bases[other_src], axes[other_src]

    query_perts = query_df["perturbation"].tolist()
    query_idx = [g2i[p] for p in query_perts]
    query_rids = query_df["id"].tolist()

    own_rep = df_own[df_own["reproducible"] == 1]
    oth_rep = df_other[df_other["reproducible"] == 1]
    tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
    tr_perts = tr_sub["perturbation"].tolist()
    tr_idx = [g2i[p] for p in tr_perts]
    tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
    tr_Y_p = project_rows_5000(tr_Y_raw, base_own, ax_own)

    pc_levels = (0, 1, 2, 5, 10) if cfg["use_deflation"] else (0,)
    pc_weights = PanPerturbBridgeModel.PC_WEIGHTS if cfg["use_deflation"] else {0: 1.0}
    soft_p = 1.5 if cfg["use_hybrid_tom"] else 1.0

    Gs_own = build_pyramid_configurable(
        df_own, all_y, g2i, base_own, pc_levels,
        self_mask=cfg["use_self_mask"], base_sub=cfg["use_base_sub"], sqrt_ncells=cfg["use_sqrt_ncells"],
    )
    Gs_oth = build_pyramid_configurable(
        df_other, all_y, g2i_o, base_oth, pc_levels,
        self_mask=cfg["use_self_mask"], base_sub=cfg["use_base_sub"], sqrt_ncells=cfg["use_sqrt_ncells"],
    )

    Gs_pan_own, Gs_pan_oth = {}, {}
    for k in pc_levels:
        wo = PanPerturbBridgeModel.PAN_OTHER_WEIGHTS[k] if cfg["use_cross_graph"] else 0.0
        M1 = Gs_own[k].copy()
        if wo > 0:
            M1[np.ix_(own_sh, own_sh)] += wo * Gs_oth[k][np.ix_(oth_sh, oth_sh)]
        Gs_pan_own[k] = M1
        M2 = Gs_oth[k].copy()
        if wo > 0:
            M2[np.ix_(oth_sh, oth_sh)] += wo * Gs_own[k][np.ix_(own_sh, own_sh)]
        Gs_pan_oth[k] = M2
    del Gs_own, Gs_oth

    nr_sub_all = df_own[(df_own["reproducible"] == 0) & (df_own["perturbation"].isin(g2i.keys()))]
    nr_Y_all_p = project_rows_5000(np.array([all_y[rid] for rid in nr_sub_all["id"]], dtype=np.float32), base_own, ax_own)
    z_eff_all = np.linalg.norm(nr_Y_all_p[:, s_idx], axis=1) * np.sqrt(nr_sub_all["n_cells"].values.astype(np.float64))
    w_nr_all = 0.25 * np.clip((z_eff_all - 0.9 * np.median(z_eff_all)) / (0.3 * np.median(z_eff_all)), 0.1, 1.5)

    # Cross-line translator
    own_p2id = dict(zip(own_rep["perturbation"], own_rep["id"]))
    oth_p2id = dict(zip(oth_rep["perturbation"], oth_rep["id"]))
    paired = sorted(list(set(own_p2id.keys()) & set(oth_p2id.keys())))
    Y_p_own = project_rows_5000(np.array([all_y[own_p2id[p]] for p in paired], dtype=np.float32), base_own, ax_own)
    Y_p_oth = project_rows_5000(np.array([all_y[oth_p2id[p]] for p in paired], dtype=np.float32), base_oth, ax_oth)
    trans_o2s = fit_geneslope_coskrr_translator(Y_p_oth, Y_p_own, oth_sh, own_sh, skip_w=0.8, alpha=2.0, p_k=1.5)
    trans_s2o = fit_geneslope_coskrr_translator(Y_p_own, Y_p_oth, own_sh, oth_sh, skip_w=0.8, alpha=2.0, p_k=1.5)

    # Own-line KRR
    tr_set = set(tr_perts)
    all_idx = list(tr_idx)
    all_Y_parts = [tr_Y_p]
    sw_parts = [np.ones(len(tr_idx))]
    extra_set = set()
    if cfg["use_bridge"]:
        extra_df = oth_rep[(oth_rep["perturbation"].isin(g2i.keys())) & (~oth_rep["perturbation"].isin(tr_set))]
        extra_perts = extra_df["perturbation"].tolist()
        extra_set = set(extra_perts)
        if extra_perts:
            Y_ext = trans_o2s(project_rows_5000(np.array([all_y[rid] for rid in extra_df["id"]], dtype=np.float32), base_oth, ax_oth))
            all_idx += [g2i[p] for p in extra_perts]
            all_Y_parts.append(Y_ext)
            sw_parts.append(np.full(len(extra_perts), 0.6))
    if cfg["use_nonrep"]:
        keep_nr = np.array([p not in extra_set for p in nr_sub_all["perturbation"].tolist()], dtype=bool)
        if np.any(keep_nr):
            all_idx += [g2i[p] for p in nr_sub_all["perturbation"].values[keep_nr].tolist()]
            all_Y_parts.append(nr_Y_all_p[keep_nr])
            sw_parts.append(w_nr_all[keep_nr])

    K_tr, K_val = build_pyramid_kernel(Gs_pan_own, all_idx, query_idx, pc_weights, soft_p=soft_p)
    pred_own = krr_weighted_solve(K_tr, K_val, np.vstack(all_Y_parts), alpha=0.5, sample_weight=np.concatenate(sw_parts))

    # Other-line KRR for bidirectional bridge
    q_in_o = [i for i, p in enumerate(query_perts) if p in g2i_o]
    q_oi = [g2i_o[query_perts[i]] for i in q_in_o]
    pred_in_other = None
    if cfg["use_bridge"] and len(q_oi) > 0:
        oth_meas = oth_rep[oth_rep["perturbation"].isin(g2i_o.keys())]
        oth_tr_perts = oth_meas["perturbation"].tolist()
        oth_tr_idx = [g2i_o[p] for p in oth_tr_perts]
        oth_Y_parts = [project_rows_5000(np.array([all_y[rid] for rid in oth_meas["id"]], dtype=np.float32), base_oth, ax_oth)]
        oth_sw_parts = [np.ones(len(oth_tr_idx), dtype=np.float64)]
        oth_set = set(oth_tr_perts)

        ext_o_df = own_rep[(own_rep["perturbation"].isin(g2i_o.keys())) & (~own_rep["perturbation"].isin(oth_set))]
        ext_o_perts = ext_o_df["perturbation"].tolist()
        if ext_o_perts:
            Y_ext_o = trans_s2o(project_rows_5000(np.array([all_y[rid] for rid in ext_o_df["id"]], dtype=np.float32), base_own, ax_own))
            oth_tr_idx += [g2i_o[p] for p in ext_o_perts]
            oth_Y_parts.append(Y_ext_o)
            oth_sw_parts.append(np.full(len(ext_o_perts), 0.6))
            oth_set |= set(ext_o_perts)
        if cfg["use_nonrep"]:
            nr_o_sub = df_other[(df_other["reproducible"] == 0) & (df_other["perturbation"].isin(g2i_o.keys())) & (~df_other["perturbation"].isin(oth_set))]
            if len(nr_o_sub) > 0:
                nr_o_Y = project_rows_5000(np.array([all_y[rid] for rid in nr_o_sub["id"]], dtype=np.float32), base_oth, ax_oth)
                z_o = np.linalg.norm(nr_o_Y[:, s_idx_o], axis=1) * np.sqrt(nr_o_sub["n_cells"].values.astype(np.float64))
                w_nro = 0.25 * np.clip((z_o - 0.9 * np.median(z_o)) / (0.3 * np.median(z_o)), 0.1, 1.5)
                oth_tr_idx += [g2i_o[p] for p in nr_o_sub["perturbation"].tolist()]
                oth_Y_parts.append(nr_o_Y)
                oth_sw_parts.append(w_nro)

        K_tr_o, K_val_o = build_pyramid_kernel(Gs_pan_oth, oth_tr_idx, q_oi, pc_weights, soft_p=soft_p)
        pred_in_other = krr_weighted_solve(K_tr_o, K_val_o, np.vstack(oth_Y_parts), alpha=0.5, sample_weight=np.concatenate(oth_sw_parts))
    del Gs_pan_own, Gs_pan_oth

    def prep_rep_matrices(rep_df, g2is, base_s, ax_s):
        Y_raw = np.array([all_y[rid] for rid in rep_df["id"]], dtype=np.float32)
        Y_p = project_rows_5000(Y_raw, base_s, ax_s)
        Z = Y_raw.copy()
        for i, p in enumerate(rep_df["perturbation"].tolist()):
            if p in g2is:
                Z[i, g2is[p]] = 0.0
        Q5, Rm = np.linalg.qr(np.stack([np.ones(5000), base_s, ax_s], axis=1))
        Q5 = Q5[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())].astype(np.float32)
        Z = Z - base_s.astype(np.float32)[None, :]
        Z = Z - (Z @ Q5) @ Q5.T
        Z_z = Z / (np.std(Z, axis=0, keepdims=True) + 1e-6)
        return Y_p, Z_z

    Y_rep_own_p, Z_rep_own_z = prep_rep_matrices(own_rep, g2i, base_own, ax_own)

    bridge_w = 0.9
    rec_w = 0.55
    out = pred_own.copy()
    if cfg["use_bridge"] and pred_in_other is not None:
        p_other = pred_in_other.copy()
        if cfg["use_reciprocal"]:
            Y_rep_oth_p, Z_rep_oth_z = prep_rep_matrices(oth_rep, g2i_o, base_oth, ax_oth)
            rec_o, conf_o = compute_relu_reciprocal(Y_rep_oth_p, Z_rep_oth_z, q_oi, thresh=4.0, pow_p=2.0)
            for j in range(len(q_oi)):
                if conf_o[j] > 0:
                    no = np.linalg.norm(p_other[j, s_idx_o]) + 1e-8
                    nr = np.linalg.norm(rec_o[j, s_idx_o]) + 1e-8
                    gate = np.clip(conf_o[j] / 3.0, 0.2, 2.0)
                    p_other[j] += (rec_w * gate) * (no / nr) * rec_o[j]
        p_br = trans_o2s(p_other)
        for j, qi in enumerate(q_in_o):
            out[qi] = out[qi] + bridge_w * p_br[j]

    if cfg["use_reciprocal"]:
        rec_own, conf_own = compute_relu_reciprocal(Y_rep_own_p, Z_rep_own_z, query_idx, thresh=4.0, pow_p=2.0)
        for i in range(len(query_idx)):
            if conf_own[i] > 0:
                no = np.linalg.norm(out[i, s_idx]) + 1e-8
                nr = np.linalg.norm(rec_own[i, s_idx]) + 1e-8
                gate = np.clip(conf_own[i] / 3.0, 0.2, 2.0)
                out[i, s_idx] += (rec_w * gate) * (no / nr) * rec_own[i, s_idx]

    if cfg["use_wiener_snap"]:
        svd_own = TruncatedSVD(n_components=min(250, len(Y_rep_own_p) - 1), random_state=42).fit(Y_rep_own_p[:, s_idx])
        V_own, S_own = svd_own.components_, svd_own.singular_values_
        Y_rep_own_sn = Y_rep_own_p[:, s_idx] / (np.linalg.norm(Y_rep_own_p[:, s_idx], axis=1, keepdims=True) + 1e-8)
        out = svd_wiener_denoise(out, V_own, S_own, s_idx, tau=0.020)
        out = manifold_snap(out, Y_rep_own_p, Y_rep_own_sn, Z_rep_own_z, query_idx, s_idx, snap_w=0.45, temp=18.0, z_bonus=0.08, top_k=40)
        out = svd_wiener_denoise(out, V_own, S_own, s_idx, tau=0.015)

    return {rid: base_own + out[i] for i, rid in enumerate(query_rids)}


FULL_CFG = {
    "use_self_mask": True,
    "use_base_sub": True,
    "use_sqrt_ncells": True,
    "use_hybrid_tom": True,
    "use_deflation": True,
    "use_cross_graph": True,
    "use_nonrep": True,
    "use_bridge": True,
    "use_reciprocal": True,
    "use_wiener_snap": True,
}

LOO_VARIANTS = [
    ("Full_S7", {}),
    ("LOO_S1_SelfMask", {"use_self_mask": False}),
    ("LOO_S1_BaseSub", {"use_base_sub": False}),
    ("LOO_S1_SqrtNCells", {"use_sqrt_ncells": False}),
    ("LOO_S1_HybridTOM", {"use_hybrid_tom": False}),
    ("LOO_S2_Deflation", {"use_deflation": False}),
    ("LOO_S3_CrossLineGraph", {"use_cross_graph": False}),
    ("LOO_S4_NonRepAnchors", {"use_nonrep": False}),
    ("LOO_S5_CrossLineBridge", {"use_bridge": False}),
    ("LOO_S6_Reciprocal", {"use_reciprocal": False}),
    ("LOO_S7_WienerSnap", {"use_wiener_snap": False}),
]


def run_dataset_lomo(dataset: str):
    t0 = time.time()
    print(f"\n=== T4: Leave-One-Module-Out Ablation on {dataset} (10 splits) ===", flush=True)
    if dataset == "xatlas":
        train_df, all_y, bases_pub, axes_pub, sidx, sh0, sh1, g2i_map = load_xatlas_bundle()
    else:
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_replogle_bundle()

    rows = []
    for seed in SEEDS_10:
        ts = time.time()
        if dataset == "xatlas":
            val_df, vr_set, vt_line = make_xatlas_split(train_df, g2i_map, seed)
            bases, axes = bases_pub, axes_pub
        else:
            val_df, df_tr, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)

        src_arr = val_df["source"].values
        row_scores = {}
        for vname, overrides in LOO_VARIANTS:
            cfg = {**FULL_CFG, **overrides}
            preds = {}
            for src in ["HCT116", "HEK293T"]:
                oth = "HEK293T" if src == "HCT116" else "HCT116"
                if dataset == "xatlas":
                    df_own = train_df[(train_df["source"] == src) & (~train_df["id"].isin(vr_set))]
                    df_oth = train_df[(train_df["source"] == oth) & (~train_df["perturbation"].isin(vt_line[src]))]
                else:
                    df_own = df_tr[df_tr["source"] == src]
                    df_oth = df_tr[df_tr["source"] == oth]
                q_df = val_df[val_df["source"] == src]
                preds.update(run_variant_for_source(src, df_own, df_oth, q_df, all_y, bases, axes, g2i_map, sh0, sh1, sidx, cfg))
            row_scores[vname] = row_scores_v5(preds, val_df, all_y, bases, axes, sidx, g2i_map)

        full_sc = row_scores["Full_S7"]
        for vname, _ in LOO_VARIANTS:
            if vname == "Full_S7":
                ov = 0.5 * (full_sc[src_arr == "HCT116"].mean() + full_sc[src_arr == "HEK293T"].mean())
                rows.append({"dataset": dataset, "seed": seed, "variant": vname, "v5_score": round(float(ov), 6),
                             "module_contrib_full_minus_loo": 0.0, "ci95_lo": 0.0, "ci95_hi": 0.0, "boot_se": 0.0, "ci_gt_zero": False, "ci_lt_zero": False})
            else:
                st_ci = stratified_row_bootstrap_ci(full_sc, row_scores[vname], src_arr)
                rows.append({
                    "dataset": dataset, "seed": seed, "variant": vname,
                    "v5_score": st_ci["score_b"],
                    "module_contrib_full_minus_loo": st_ci["delta"],
                    "ci95_lo": st_ci["ci95_lo"], "ci95_hi": st_ci["ci95_hi"],
                    "boot_se": st_ci["boot_se"], "ci_gt_zero": st_ci["ci_gt_zero"], "ci_lt_zero": st_ci["ci_lt_zero"],
                })
        print(f"  [{dataset}] Seed {seed:5d} complete ({len(LOO_VARIANTS)} variants) [{time.time()-ts:.1f}s]", flush=True)
        gc.collect()

    df_out = pd.DataFrame(rows)
    df_out.to_csv(REPORTS_DIR / f"t4_lomo_{dataset}.csv", index=False)
    print(f"Saved t4_lomo_{dataset}.csv in {time.time()-t0:.1f}s!", flush=True)
    return df_out


def main():
    df_x = run_dataset_lomo("xatlas")
    df_r = run_dataset_lomo("replogle2022")

    summary = {}
    for ds, df in [("xatlas", df_x), ("replogle2022", df_r)]:
        summary[ds] = {}
        print(f"\n=== T4 LOMO SUMMARY ({ds}) — Module Contribution = Full_S7 - LOO ===")
        for vname, _ in LOO_VARIANTS:
            sub = df[df["variant"] == vname]
            summary[ds][vname] = {
                "mean_v5": round(float(sub["v5_score"].mean()), 4),
                "mean_contrib": round(float(sub["module_contrib_full_minus_loo"].mean()), 4),
                "sd_contrib": round(float(sub["module_contrib_full_minus_loo"].std()), 4),
                "median_ci95_lo": round(float(sub["ci95_lo"].median()), 4),
                "median_ci95_hi": round(float(sub["ci95_hi"].median()), 4),
                "splits_pos": f"{int((sub['module_contrib_full_minus_loo'] > 0).sum())}/10",
                "splits_ci_gt_zero": f"{int(sub['ci_gt_zero'].sum())}/10",
                "splits_ci_lt_zero": f"{int(sub['ci_lt_zero'].sum())}/10",
            }
            s = summary[ds][vname]
            print(f"  {vname:24s} | v5={s['mean_v5']:.4f} | contrib={s['mean_contrib']:+.4f} (med CI [{s['median_ci95_lo']:+.4f}, {s['median_ci95_hi']:+.4f}]) | pos={s['splits_pos']} | CI>0:{s['splits_ci_gt_zero']} | CI<0:{s['splits_ci_lt_zero']}")

    with open(REPORTS_DIR / "t4_lomo_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\nSaved reports/t4_lomo_summary.json!")


if __name__ == "__main__":
    main()
