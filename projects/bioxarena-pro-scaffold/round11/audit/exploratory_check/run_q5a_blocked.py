#!/usr/bin/env python3
"""
Exploratory Check 3: Direct Q5a Cross-Line Conservation (`logFC` vs `log1p`) & Direct
Transfer (`FCT` vs `log1p`) across all 33 Human Cell-Line Pairs (`>= 25` shared reproducible
targets) and ALL 3 gene universes (`all`, `expr01`, `expr05`).

Replaces invalid naive edge-level Wilcoxon (`N=33` dependent edges on a 10-line graph)
with dependency-aware inference:
  1. Target-line-blocked exact sign-flip permutation (`2^10 = 1,024` over 10 target-line blocks)
     and Source-line-blocked exact sign-flip permutation (`2^10 = 1,024` over 10 source blocks)
     on all `66` directed pairs (`FCT` vs `log1p`).
  2. Line-summarized exact sign-flip permutation (`2^10 = 1,024` on the 10 active lines'
     per-line edge means d_bar_v) for both `logFC` vs `log1p` and `FCT` vs `log1p`.
  3. Two-sided & one-sided Exact Line-Level Binomial Sign Test (`N = 10` lines) +
     LOO-Line Hub-Drop sensitivity analysis (`min_{-v}` across dropping each line).
  4. Dyadic Two-Way Cluster-Robust standard error & t-stat (Cameron-Gelbach-Miller 2011).
  5. Explicit `Bonferroni x 3` correction across the 3 gene universes (`all`, `expr01`, `expr05`).
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
from scipy.stats import binomtest, t as student_t, wilcoxon

EXPLORATORY_CHECK_DIR = Path(__file__).resolve().parent
AUDIT_DIR = EXPLORATORY_CHECK_DIR.parent
for p in (str(EXPLORATORY_CHECK_DIR), str(AUDIT_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

from exploratory_core import LINES_12, load_12_lines  # noqa: E402
from independent_r11_loader import AuditLineData  # noqa: E402
from independent_r11_metrics import EPS, ceiling_from_halves, disattenuate, pearson_vec  # noqa: E402

GENE_UNIVERSES: Tuple[Tuple[str, float], ...] = (
    ("all", -1.0),
    ("expr01", 0.1),
    ("expr05", 0.5),
)


def _log1p_to_logfc(deltas: np.ndarray, mu_ctrl: np.ndarray, floor: float = 1e-8) -> np.ndarray:
    """Invert Delta = log1p(pert) - log1p(mu) -> log(pert / mu) safely."""
    mu = np.maximum(np.asarray(mu_ctrl, dtype=np.float64), floor)
    lmu = np.log1p(mu)[None, :]
    pert = np.maximum(np.expm1(lmu + np.asarray(deltas, dtype=np.float64)), floor)
    return np.log(pert / mu[None, :])


def _fct_map(logfc_src: np.ndarray, mu_tgt: np.ndarray) -> np.ndarray:
    """Delta_FCT(src -> tgt) = log1p(max(mu_tgt, 0) * exp(logFC_src)) - log1p(max(mu_tgt, 0))."""
    mu = np.maximum(np.asarray(mu_tgt, dtype=np.float64), 0.0)
    # Clip exponent safely before exp only against IEEE overflow while preserving dynamic range
    fc = np.exp(np.clip(logfc_src, -40.0, 40.0))
    return np.log1p(mu[None, :] * fc) - np.log1p(mu)[None, :]


def _pool_corr(vecs_a: Sequence[np.ndarray], vecs_b: Sequence[np.ndarray]) -> float:
    if not vecs_a:
        return float("nan")
    a = np.concatenate(vecs_a)
    b = np.concatenate(vecs_b)
    return pearson_vec(a, b)


def _exact_line_sign_flip_p(block_values: np.ndarray) -> Tuple[float, float]:
    """
    Exact 2^V enumeration of block sign-flip permutation test on mean(block_values).
    Returns (p_two_sided, p_one_sided_greater).
    """
    v = int(block_values.size)
    if v == 0 or not np.all(np.isfinite(block_values)):
        return float("nan"), float("nan")
    obs = float(np.mean(block_values))
    # Exhaustive 2^V via binary matrix (V <= 15; here V = 10 -> 1,024 configs)
    signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=v)), dtype=np.float64)
    null_means = (signs @ block_values) / float(v)
    tol = 1e-15
    p_one = float(np.mean(null_means >= obs - tol))
    p_two = float(np.mean(np.abs(null_means) >= abs(obs) - tol))
    return p_two, p_one


def _exact_binomial_sign_test(values: np.ndarray) -> Tuple[int, int, float, float]:
    fin = values[np.isfinite(values)]
    nonzero = fin[np.abs(fin) > 1e-15]
    n = int(nonzero.size)
    k = int(np.sum(nonzero > 0.0))
    if n == 0:
        return 0, 0, float("nan"), float("nan")
    p_two = float(binomtest(k, n, p=0.5, alternative="two-sided").pvalue)
    p_one = float(binomtest(k, n, p=0.5, alternative="greater").pvalue)
    return k, n, p_two, p_one


def _dyadic_cluster_test(
    pairs: Sequence[Tuple[str, str]],
    diffs: np.ndarray,
    active_lines: Sequence[str],
) -> Dict[str, Any]:
    """
    Cameron-Gelbach-Miller (2011) two-way / dyadic cluster variance over shared lines,
    plus Conservative Single-Way / LOO Jackknife checks and LOO-line minimum remaining mean.
    """
    m = len(diffs)
    v_cnt = len(active_lines)
    mu_hat = float(np.mean(diffs))
    resid = diffs - mu_hat

    # Edge covariance: V_hat = sum_{e cap e' != empty} r_e r_{e'}
    # = sum_v (sum_{e ni v} r_e)^2 - sum_e r_e^2
    sum_by_v = np.zeros(v_cnt, dtype=np.float64)
    line_to_idx = {ln: i for i, ln in enumerate(active_lines)}
    deg = np.zeros(v_cnt, dtype=np.int64)
    per_line_vals: List[List[float]] = [[] for _ in range(v_cnt)]
    for (a, b), r_val, d_val in zip(pairs, resid, diffs):
        ia, ib = line_to_idx[a], line_to_idx[b]
        sum_by_v[ia] += r_val
        sum_by_v[ib] += r_val
        deg[ia] += 1
        deg[ib] += 1
        per_line_vals[ia].append(float(d_val))
        per_line_vals[ib].append(float(d_val))

    line_means = np.asarray([float(np.mean(vs)) for vs in per_line_vals], dtype=np.float64)
    s_cgm = float(np.sum(sum_by_v ** 2) - np.sum(resid ** 2))
    # Fallback floor: max single-line-sum variance or naive variance
    s_diag = float(np.sum(resid ** 2))
    # One-way upper-bound variance via each node sum
    var_mean_cgm = s_cgm / float(m * m) if s_cgm > 0.0 else float("nan")
    se_cgm = math.sqrt(var_mean_cgm) if np.isfinite(var_mean_cgm) and var_mean_cgm > 0 else float("nan")
    t_cgm = mu_hat / se_cgm if np.isfinite(se_cgm) and se_cgm > EPS else float("nan")
    # Conservative df = V - 1 (10 - 1 = 9 lines)
    p_cgm_two = float(2.0 * student_t.sf(abs(t_cgm), df=max(v_cnt - 1, 1))) if np.isfinite(t_cgm) else float("nan")
    p_cgm_one = float(student_t.sf(t_cgm, df=max(v_cnt - 1, 1))) if np.isfinite(t_cgm) else float("nan")

    # Line-level inference on the V=10 per-line averages d_bar_v
    p_line_perm_two, p_line_perm_one = _exact_line_sign_flip_p(line_means)
    k_pos_lines, n_lines_nz, p_binom_two, p_binom_one = _exact_binomial_sign_test(line_means)

    # LOO-line drop check (dropping any single line v and all edges incident to v)
    loo_means: Dict[str, float] = {}
    for ln in active_lines:
        keep_d = [float(d) for (a, b), d in zip(pairs, diffs) if a != ln and b != ln]
        loo_means[ln] = float(np.mean(keep_d)) if keep_d else float("nan")
    min_loo_line = min(loo_means, key=lambda k: loo_means[k])
    max_loo_line = max(loo_means, key=lambda k: loo_means[k])

    # Naive edge Wilcoxon for explicit side-by-side audit comparison
    n_pos_edges = int(np.sum(diffs > 0.0))
    try:
        w_stat_two = float(wilcoxon(diffs, alternative="two-sided").pvalue)
        w_stat_one = float(wilcoxon(diffs, alternative="greater").pvalue)
    except Exception:
        w_stat_two, w_stat_one = float("nan"), float("nan")

    return {
        "n_pairs": m,
        "n_active_lines": v_cnt,
        "mean_diff": mu_hat,
        "median_diff": float(np.median(diffs)),
        "n_pos_pairs": n_pos_edges,
        "naive_wilcoxon_p_two": w_stat_two,
        "naive_wilcoxon_p_one": w_stat_one,
        "line_means": {ln: float(lm) for ln, lm in zip(active_lines, line_means)},
        "n_pos_lines": k_pos_lines,
        "n_lines_tested": n_lines_nz,
        "line_binomial_sign_p_two": p_binom_two,
        "line_binomial_sign_p_one": p_binom_one,
        "line_exact_signflip_p_two": p_line_perm_two,
        "line_exact_signflip_p_one": p_line_perm_one,
        "dyadic_cgm_se": se_cgm,
        "dyadic_cgm_t_df9": t_cgm,
        "dyadic_cgm_p_two": p_cgm_two,
        "dyadic_cgm_p_one": p_cgm_one,
        "cgm_var_raw": s_cgm,
        "naive_var_raw": s_diag,
        "loo_line_min_drop": {"line": min_loo_line, "remaining_mean": loo_means[min_loo_line]},
        "loo_line_max_drop": {"line": max_loo_line, "remaining_mean": loo_means[max_loo_line]},
        "loo_line_all_positive": bool(all(val > 0.0 for val in loo_means.values())),
    }


def _directed_block_tests(
    directed_pairs: Sequence[Tuple[str, str]],
    diffs: np.ndarray,
    active_lines: Sequence[str],
) -> Dict[str, Any]:
    """
    Exact Line-Blocked Sign-Flip Permutation over all 2^10 = 1,024 configurations
    grouping all 66 directed `(source -> target)` pairs either by target-line block
    (identical target matrix / ceiling / mu_t) or by source-line block.
    """
    # Group by target line & compute equal-line-weighted block mean AND pair-weighted block sum
    tgt_blocks = {ln: [] for ln in active_lines}
    src_blocks = {ln: [] for ln in active_lines}
    for (s, t), d in zip(directed_pairs, diffs):
        tgt_blocks[t].append(float(d))
        src_blocks[s].append(float(d))

    tgt_means = np.asarray([float(np.mean(tgt_blocks[ln])) for ln in active_lines], dtype=np.float64)
    tgt_sums = np.asarray([float(np.sum(tgt_blocks[ln])) for ln in active_lines], dtype=np.float64)
    src_means = np.asarray([float(np.mean(src_blocks[ln])) for ln in active_lines], dtype=np.float64)
    src_sums = np.asarray([float(np.sum(src_blocks[ln])) for ln in active_lines], dtype=np.float64)

    p_tgt_unw_two, p_tgt_unw_one = _exact_line_sign_flip_p(tgt_means)
    p_tgt_pw_two, p_tgt_pw_one = _exact_line_sign_flip_p(tgt_sums)
    p_src_unw_two, p_src_unw_one = _exact_line_sign_flip_p(src_means)
    p_src_pw_two, p_src_pw_one = _exact_line_sign_flip_p(src_sums)

    k_tgt, n_tgt, p_tgt_bin_two, p_tgt_bin_one = _exact_binomial_sign_test(tgt_means)
    k_src, n_src, p_src_bin_two, p_src_bin_one = _exact_binomial_sign_test(src_means)

    return {
        "n_directed_pairs": int(diffs.size),
        "n_blocks": len(active_lines),
        "target_block_means": {ln: float(m) for ln, m in zip(active_lines, tgt_means)},
        "target_block_n_pos": k_tgt,
        "target_block_binomial_p_two": p_tgt_bin_two,
        "target_block_binomial_p_one": p_tgt_bin_one,
        "target_block_signflip_pairweighted_p_two": p_tgt_pw_two,
        "target_block_signflip_pairweighted_p_one": p_tgt_pw_one,
        "target_block_signflip_lineweighted_p_two": p_tgt_unw_two,
        "target_block_signflip_lineweighted_p_one": p_tgt_unw_one,
        "source_block_means": {ln: float(m) for ln, m in zip(active_lines, src_means)},
        "source_block_n_pos": k_src,
        "source_block_binomial_p_two": p_src_bin_two,
        "source_block_binomial_p_one": p_src_bin_one,
        "source_block_signflip_pairweighted_p_two": p_src_pw_two,
        "source_block_signflip_pairweighted_p_one": p_src_pw_one,
        "source_block_signflip_lineweighted_p_two": p_src_unw_two,
        "source_block_signflip_lineweighted_p_one": p_src_unw_one,
    }


def max_independent_matching_size(pairs: Sequence[Tuple[str, str]], nodes: Sequence[str]) -> int:
    """Exact maximum vertex-disjoint matching size on the V<=12 graph via backtracking."""
    adj: Dict[int, List[int]] = {i: [] for i in range(len(nodes))}
    idx = {ln: i for i, ln in enumerate(nodes)}
    for a, b in pairs:
        u, v = idx[a], idx[b]
        adj[u].append(v)
        adj[v].append(u)

    best = [0]

    def _dfs(u_start: int, used: int, cnt: int) -> None:
        if cnt > best[0]:
            best[0] = cnt
        for u in range(u_start, len(nodes)):
            if used & (1 << u):
                continue
            for v in adj[u]:
                if v > u and not (used & (1 << v)):
                    _dfs(u + 1, used | (1 << u) | (1 << v), cnt + 1)

    _dfs(0, 0, 0)
    return best[0]


def run_check3(out_dir: Path = EXPLORATORY_CHECK_DIR / "outputs") -> Dict[str, Any]:
    t0 = time.time()
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = load_12_lines()
    line_names = [ln for ln, _ in LINES_12]

    # Map target gene -> row index (among reproducible targets only) per line
    rep_map: Dict[str, Dict[str, int]] = {}
    for ln in line_names:
        ds = lines[ln]
        rmap: Dict[str, int] = {}
        for i, tgt in enumerate(ds.target_names):
            if ds.repr_mask is None or bool(ds.repr_mask[i]):
                rmap[tgt] = i
        rep_map[ln] = rmap

    # Identify all undirected pairs with >= 25 shared reproducible targets
    pairs_33: List[Tuple[str, str, List[str]]] = []
    for i, a in enumerate(line_names):
        for b in line_names[i + 1 :]:
            shared = sorted(set(rep_map[a]) & set(rep_map[b]))
            if len(shared) >= 25:
                pairs_33.append((a, b, shared))

    active_lines = sorted({ln for a, b, _ in pairs_33 for ln in (a, b)})
    max_indep = max_independent_matching_size([(a, b) for a, b, _ in pairs_33], active_lines)
    print(f"[Check3 Q5a] Pairs={len(pairs_33)}, Active lines={len(active_lines)}, Max vertex-disjoint pairs={max_indep}")

    # Precompute full logFC and split-half logFC arrays per line cleanly
    logfc_cache: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for ln in active_lines:
        ds = lines[ln]
        lfc = _log1p_to_logfc(ds.deltas, ds.ctrl_mean_libnorm)
        lfc_a = _log1p_to_logfc(ds.half_a_deltas, ds.ctrl_mean_libnorm)
        lfc_b = _log1p_to_logfc(ds.half_b_deltas, ds.ctrl_mean_libnorm)
        logfc_cache[ln] = (lfc, lfc_a, lfc_b)

    all_pair_rows: List[Dict[str, Any]] = []
    universe_results: Dict[str, Any] = {}

    for u_name, mu_thr in GENE_UNIVERSES:
        undirected_pairs_u: List[Tuple[str, str]] = []
        directed_pairs_u: List[Tuple[str, str]] = []

        raw_log1p_v, cal_log1p_v = [], []
        raw_logfc_v, cal_logfc_v = [], []
        raw_fct_avg_v, cal_fct_avg_v = [], []
        raw_fct_dir_v, cal_fct_dir_v = [], []
        raw_log1p_dir_v, cal_log1p_dir_v = [], []

        for a, b, shared in pairs_33:
            da, db = lines[a], lines[b]
            pos_b = {g: j for j, g in enumerate(db.genes)}
            shared_pairs = [(i, pos_b[g], g) for i, g in enumerate(da.genes) if g in pos_b]
            ga = np.asarray([p[0] for p in shared_pairs], dtype=np.int64)
            gb = np.asarray([p[1] for p in shared_pairs], dtype=np.int64)
            shared_g_pos = {p[2]: k for k, p in enumerate(shared_pairs)}
            n_g = len(shared_pairs)

            lfca_full, lfc_ha_full, lfc_hb_full = logfc_cache[a]
            lfcb_full, lfc_hba_full, lfc_hbb_full = logfc_cache[b]
            lfca, lfc_ha, lfc_hb = lfca_full[:, ga], lfc_ha_full[:, ga], lfc_hb_full[:, ga]
            lfcb, lfc_hba, lfc_hbb = lfcb_full[:, gb], lfc_hba_full[:, gb], lfc_hbb_full[:, gb]
            da_del, da_ha, da_hb = da.deltas[:, ga], da.half_a_deltas[:, ga], da.half_b_deltas[:, ga]
            db_del, db_ha, db_hb = db.deltas[:, gb], db.half_a_deltas[:, gb], db.half_b_deltas[:, gb]

            ra_idx = np.asarray([rep_map[a][t] for t in shared], dtype=np.int64)
            rb_idx = np.asarray([rep_map[b][t] for t in shared], dtype=np.int64)

            mu_a = np.asarray(da.ctrl_mean_libnorm[ga], dtype=np.float64)
            mu_b = np.asarray(db.ctrl_mean_libnorm[gb], dtype=np.float64)
            gene_ok = np.isfinite(mu_a) & np.isfinite(mu_b)
            if mu_thr >= 0.0:
                gene_ok &= (mu_a >= mu_thr) & (mu_b >= mu_thr)

            # Build element lists for pair (a, b) excluding self-target column
            va_l1p, vb_l1p = [], []
            va_lfc, vb_lfc = [], []
            vfct_ab, vfct_ba = [], []
            ha_a_l1p, hb_a_l1p, ha_b_l1p, hb_b_l1p = [], [], [], []
            ha_a_lfc, hb_a_lfc, ha_b_lfc, hb_b_lfc = [], [], [], []

            fct_ab_mat = _fct_map(lfca[ra_idx], mu_b)
            fct_ba_mat = _fct_map(lfcb[rb_idx], mu_a)

            for k, (ia, ib, tgt_sym) in enumerate(zip(ra_idx, rb_idx, shared)):
                mask = (
                    gene_ok
                    & np.isfinite(da_del[ia])
                    & np.isfinite(db_del[ib])
                    & np.isfinite(da_ha[ia])
                    & np.isfinite(da_hb[ia])
                    & np.isfinite(db_ha[ib])
                    & np.isfinite(db_hb[ib])
                )
                t_shared_col = shared_g_pos.get(tgt_sym, -1)
                if 0 <= t_shared_col < n_g:
                    mask[t_shared_col] = False
                idx = np.flatnonzero(mask)
                if idx.size == 0:
                    continue
                va_l1p.append(da_del[ia, idx].astype(np.float64))
                vb_l1p.append(db_del[ib, idx].astype(np.float64))
                va_lfc.append(lfca[ia, idx])
                vb_lfc.append(lfcb[ib, idx])
                vfct_ab.append(fct_ab_mat[k, idx])
                vfct_ba.append(fct_ba_mat[k, idx])

                ha_a_l1p.append(da_ha[ia, idx].astype(np.float64))
                hb_a_l1p.append(da_hb[ia, idx].astype(np.float64))
                ha_b_l1p.append(db_ha[ib, idx].astype(np.float64))
                hb_b_l1p.append(db_hb[ib, idx].astype(np.float64))

                ha_a_lfc.append(lfc_ha[ia, idx])
                hb_a_lfc.append(lfc_hb[ia, idx])
                ha_b_lfc.append(lfc_hba[ib, idx])
                hb_b_lfc.append(lfc_hbb[ib, idx])

            r_l1p_raw = _pool_corr(va_l1p, vb_l1p)
            r_lfc_raw = _pool_corr(va_lfc, vb_lfc)
            r_fct_ab_raw = _pool_corr(vfct_ab, vb_l1p)
            r_fct_ba_raw = _pool_corr(vfct_ba, va_l1p)
            r_fct_avg_raw = 0.5 * (r_fct_ab_raw + r_fct_ba_raw)

            c_a_l1p = ceiling_from_halves(_pool_corr(ha_a_l1p, hb_a_l1p))
            c_b_l1p = ceiling_from_halves(_pool_corr(ha_b_l1p, hb_b_l1p))
            c_pair_l1p = c_a_l1p * c_b_l1p

            c_a_lfc = ceiling_from_halves(_pool_corr(ha_a_lfc, hb_a_lfc))
            c_b_lfc = ceiling_from_halves(_pool_corr(ha_b_lfc, hb_b_lfc))
            c_pair_lfc = c_a_lfc * c_b_lfc

            r_l1p_cal = disattenuate(r_l1p_raw, c_pair_l1p)
            r_lfc_cal = disattenuate(r_lfc_raw, c_pair_lfc)
            r_fct_ab_cal = disattenuate(r_fct_ab_raw, c_pair_l1p)
            r_fct_ba_cal = disattenuate(r_fct_ba_raw, c_pair_l1p)
            r_fct_avg_cal = 0.5 * (r_fct_ab_cal + r_fct_ba_cal)

            undirected_pairs_u.append((a, b))
            raw_log1p_v.append(r_l1p_raw)
            cal_log1p_v.append(r_l1p_cal)
            raw_logfc_v.append(r_lfc_raw)
            cal_logfc_v.append(r_lfc_cal)
            raw_fct_avg_v.append(r_fct_avg_raw)
            cal_fct_avg_v.append(r_fct_avg_cal)

            directed_pairs_u.extend([(a, b), (b, a)])
            raw_fct_dir_v.extend([r_fct_ab_raw, r_fct_ba_raw])
            cal_fct_dir_v.extend([r_fct_ab_cal, r_fct_ba_cal])
            raw_log1p_dir_v.extend([r_l1p_raw, r_l1p_raw])
            cal_log1p_dir_v.extend([r_l1p_cal, r_l1p_cal])

            all_pair_rows.append({
                "universe": u_name,
                "line_a": a,
                "line_b": b,
                "n_shared_repro": len(shared),
                "raw_log1p": r_l1p_raw,
                "raw_logfc": r_lfc_raw,
                "raw_fct_avg": r_fct_avg_raw,
                "raw_fct_a_to_b": r_fct_ab_raw,
                "raw_fct_b_to_a": r_fct_ba_raw,
                "cal_log1p": r_l1p_cal,
                "cal_logfc": r_lfc_cal,
                "cal_fct_avg": r_fct_avg_cal,
                "cal_fct_a_to_b": r_fct_ab_cal,
                "cal_fct_b_to_a": r_fct_ba_cal,
                "c_sb_a_log1p": c_a_l1p,
                "c_sb_b_log1p": c_b_l1p,
                "c_sb_a_logfc": c_a_lfc,
                "c_sb_b_logfc": c_b_lfc,
            })

        arr_raw_l1p = np.asarray(raw_log1p_v, dtype=np.float64)
        arr_cal_l1p = np.asarray(cal_log1p_v, dtype=np.float64)
        arr_raw_lfc = np.asarray(raw_logfc_v, dtype=np.float64)
        arr_cal_lfc = np.asarray(cal_logfc_v, dtype=np.float64)
        arr_raw_fct = np.asarray(raw_fct_avg_v, dtype=np.float64)
        arr_cal_fct = np.asarray(cal_fct_avg_v, dtype=np.float64)

        contrasts = {
            "raw_logfc_minus_log1p": _dyadic_cluster_test(undirected_pairs_u, arr_raw_lfc - arr_raw_l1p, active_lines),
            "cal_logfc_minus_log1p": _dyadic_cluster_test(undirected_pairs_u, arr_cal_lfc - arr_cal_l1p, active_lines),
            "raw_fct_minus_log1p_undir": _dyadic_cluster_test(undirected_pairs_u, arr_raw_fct - arr_raw_l1p, active_lines),
            "cal_fct_minus_log1p_undir": _dyadic_cluster_test(undirected_pairs_u, arr_cal_fct - arr_cal_l1p, active_lines),
            "raw_fct_minus_log1p_dir_blocks": _directed_block_tests(
                directed_pairs_u, np.asarray(raw_fct_dir_v) - np.asarray(raw_log1p_dir_v), active_lines
            ),
            "cal_fct_minus_log1p_dir_blocks": _directed_block_tests(
                directed_pairs_u, np.asarray(cal_fct_dir_v) - np.asarray(cal_log1p_dir_v), active_lines
            ),
        }
        # Append explicit Bonferroni x3 across the 3 gene universes
        for cname, cdict in contrasts.items():
            bonf = {}
            for pk, pv in list(cdict.items()):
                if pk.endswith(("_p_two", "_p_one")) and isinstance(pv, (int, float)) and np.isfinite(pv):
                    bonf[f"{pk}_bonf3"] = min(1.0, float(pv) * 3.0)
            cdict.update(bonf)

        summary_means = {
            "raw_log1p_mean": float(np.mean(arr_raw_l1p)),
            "raw_logfc_mean": float(np.mean(arr_raw_lfc)),
            "raw_fct_mean": float(np.mean(arr_raw_fct)),
            "cal_log1p_mean": float(np.mean(arr_cal_l1p)),
            "cal_logfc_mean": float(np.mean(arr_cal_lfc)),
            "cal_fct_mean": float(np.mean(arr_cal_fct)),
        }
        universe_results[u_name] = {
            "means": summary_means,
            "contrasts": contrasts,
        }

        c_lfc = contrasts["cal_logfc_minus_log1p"]
        c_fct = contrasts["cal_fct_minus_log1p_undir"]
        c_dir = contrasts["cal_fct_minus_log1p_dir_blocks"]
        print(
            f"[Check3 Q5a] {u_name:6s} | Cal log1p={summary_means['cal_log1p_mean']:+.4f} "
            f"logFC={summary_means['cal_logfc_mean']:+.4f} FCT={summary_means['cal_fct_mean']:+.4f}\n"
            f"   Cal logFC-log1p: edges={c_lfc['n_pos_pairs']}/33 (naive_w_p2={c_lfc['naive_wilcoxon_p_two']:.2e}) | "
            f"lines={c_lfc['n_pos_lines']}/10 line_perm_p2={c_lfc['line_exact_signflip_p_two']:.4f} (bonf3={c_lfc['line_exact_signflip_p_two_bonf3']:.4f}), "
            f"cgm_p2={c_lfc['dyadic_cgm_p_two']:.4f} (bonf3={c_lfc['dyadic_cgm_p_two_bonf3']:.4f})\n"
            f"   Cal FCT-log1p:   edges={c_fct['n_pos_pairs']}/33 (naive_w_p2={c_fct['naive_wilcoxon_p_two']:.2e}) | "
            f"lines={c_fct['n_pos_lines']}/10 line_perm_p2={c_fct['line_exact_signflip_p_two']:.4f} (bonf3={c_fct['line_exact_signflip_p_two_bonf3']:.4f}), "
            f"tgt_blk_p2={c_dir['target_block_signflip_lineweighted_p_two']:.4f} (bonf3={c_dir['target_block_signflip_lineweighted_p_two_bonf3']:.4f}), "
            f"src_blk_p2={c_dir['source_block_signflip_lineweighted_p_two']:.4f} (bonf3={c_dir['source_block_signflip_lineweighted_p_two_bonf3']:.4f}), "
            f"loo_all_pos={c_fct['loo_line_all_positive']}"
        )

    payload = {
        "check": "Q5a_blocked_and_loo_line_audit",
        "n_pairs": len(pairs_33),
        "active_lines": active_lines,
        "max_vertex_disjoint_pairs": max_indep,
        "bonferroni_factor": len(GENE_UNIVERSES),
        "elapsed_sec": round(time.time() - t0, 2),
        "universes": universe_results,
    }
    with open(out_dir / "q5a_blocked.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    with open(out_dir / "q5a_33pairs_all_universes.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(all_pair_rows[0].keys()))
        w.writeheader()
        w.writerows(all_pair_rows)
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=EXPLORATORY_CHECK_DIR / "outputs")
    args = parser.parse_args()
    run_check3(out_dir=args.out_dir)
