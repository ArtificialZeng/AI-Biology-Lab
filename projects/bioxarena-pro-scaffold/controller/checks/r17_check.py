"""Controller (third-party) recomputation for Round 17 from a party's exported arrays.

Written from PREREG_J v2 (sections 1.2, 3.1-3.5, 5.2), PREREG_J_CLARIFICATIONS_C1/C2 and
round17/src/R17_ARRAY_SCHEMA.md only. It never imports or reads round17/src/*.py or round17/audit/*.py.

Metric: r_orth_cal exactly as the controller's Round 16 code (r16_check.Line / pooled_r / cal, imported unchanged):
uncentered orthogonalization against Transfer_Mean, pooled moment correlation, Spearman-Brown ceiling from the
orthogonalized half A / half B residuals of half-valid targets (half-specific when the weights cover one half),
and joint perturbation bootstrap weights bincount(rng.integers(0, N, N)) per draw from
np.random.default_rng(seed). Bootstraps are evaluated as W @ stats (same formula as r16_check.cal per draw).

Usage:
  python3 -I r17_check.py <results_dir> [--party auto|main|at3] [--tol 1e-5] [--out recompute.json]
         [--n-boot N] [--n-boot-bonf56 N] [--n-perm N] [--n-splits-diag N] [--seed S] [--skip-j5b]
  <results_dir> holds arrays/<line>/{observed.npz, preds_<M>_sweep56.npz, ...} and summary_r17.json.
  n_boot etc. default to the values recorded in the party's summary_r17.json (4000 / 20000 / 10000 / 20 if absent).

Output: one JSON line {"n_checks", "n_fail", "n_not_exported", "worst_abs_diff", ...}, the list of failing fields
(with a diagnosis hint where the controller knows one), then R17_CHECK_MATCH or R17_CHECK_DIFF.
`--self-test` runs the adjudication branch tests only.

Controller readings where the locked text leaves a choice (each is reported, none changes a locked rule):
  * Holm (C2.1): reject when p <= alpha_step (standard Holm; PREREG_J writes p_var < alpha_adj for the base
    rule, which is used with strict < on non-Holm lines). Ties in raw p are ordered by unit (PREREG_J §2.1 order),
    then test (H-J1, H-J3a, H-J3b), then model (pooled family); if a party uses another valid p-sorted tie order,
    Holm is re-run in the party's order and a NOTE is printed. A WARNING is printed if any C2.1 gate depends on the
    tie order (gate passes at one end of its tie group's alpha range but not the other).
  * C1.7: subreasons are compared only for INCONCLUSIVE calls (a party label on other verdicts is ignored).
  * C2.3 lines (smaller half < 2): cross-fitted fields are not compared; the H-J1 call must be INCONCLUSIVE /
    VALIDITY_FAIL; no J5a cross-fitted pair may be reported; the J2 predicate uses L1_ens56 only (oracle undefined).
  * J5b: Spearman by exact integer arithmetic when there are no rank ties, so |rho_k| >= |rho_obs| counts exact ties
    (PREREG_J §3.5 / C2.2, no tolerance); permutation table = rng.permutation(56) per row from default_rng(seed).
  * J5a: r_orth_cal only (r_dorth_cal is not recomputed).
"""
import argparse
import json
import math
import sys
from pathlib import Path
from statistics import NormalDist

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r16_check as R16  # noqa: E402  (controller metric code, unchanged)
import r16_check_adjudicate as ADJ  # noqa: E402

MODELS = ("STATE_SE", "STATE_HVG")
NB_IDX = 56
ALPHA_FW = 0.05
ALPHA_LINE = 0.05 / 3          # PREREG_J §3.1 single-line base level (98.333% CI)
DELTA = 0.04                   # §3.2 practical threshold
TAU_THR = 0.02                 # §3.2: 2 tau > 0.04
MIN_VALID = 0.95
ALPHA_BONF56 = 0.025 / 56      # C1.8 per-index two-sided level
Z_BONF56 = NormalDist().inv_cdf(1 - ALPHA_BONF56 / 2)
CONF_ORDER = ("cd4t_stim48", "gse295214_h1", "pfizer_telohaec_basal")   # PREREG_J §2.1 unit order
HOLM_TESTS = ("H-J1", "H-J3a", "H-J3b")
ENS8 = (0, 7, 14, 21, 28, 35, 42, 49)
NOT_FRESH = ("a549_ko_cas9pos",)  # C1.10: not a fresh line
PFIZER_PAIR = ("pfizer_telohaec_basal", "pfizer_telohaec_cytokine")  # C1.3 (cytokine == treated)

np.seterr(all="ignore")


# ----------------------------------------------------------------------------------------------- metric
def pooled_many(S):
    """Vectorized r16_check.pooled_r over rows of S (B, 6)."""
    n, sx, sy, sxy, sxx, syy = S.T
    cov = sxy - sx * sy / n
    vx = sxx - sx * sx / n
    vy = syy - sy * sy / n
    ok = (vx > 0) & (vy > 0)
    out = np.full(S.shape[0], np.nan)
    out[ok] = cov[ok] / np.sqrt(vx[ok] * vy[ok])
    return out


class Weights:
    """Bootstrap (or indicator) weight matrix W (B, N) with the Spearman-Brown term precomputed."""

    def __init__(self, L, W):
        self.W = np.atleast_2d(np.asarray(W, float))
        self.rab = pooled_many((self.W * L.HV) @ L.cs)
        ok = self.rab > 0
        self.den = np.full(self.W.shape[0], np.nan)
        self.den[ok] = np.sqrt(2 * self.rab[ok] / (1 + self.rab[ok]))

    def r(self, st):
        out = np.clip(pooled_many(self.W @ st) / self.den, -1, 1)
        return out


def boot_weights(N, n_boot, seed, members=None):
    """Joint (members None) or within-half bootstrap weights: per draw rng.integers(0, n, n) over the members."""
    rng = np.random.default_rng(seed)
    if members is None:
        members = np.arange(N)
    members = np.asarray(members, int)
    n = members.size
    W = np.zeros((n_boot, N))
    for u in range(n_boot):
        W[u] = np.bincount(members[rng.integers(0, n, n)], minlength=N)
    return W


def indicator(N, members):
    w = np.zeros(N)
    w[np.asarray(members, int)] = 1.0
    return w


def pct_ci(draws, alpha):
    d = np.asarray(draws, float)
    d = d[np.isfinite(d)]
    if d.size == 0:
        return (math.nan, math.nan)
    lo, hi = np.percentile(d, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return (float(lo), float(hi))


def p_two_sided(draws):
    """C1.5: min(1, 2 (1 + min(#{d<=0}, #{d>=0})) / (B_valid + 1))."""
    d = np.asarray(draws, float)
    d = d[np.isfinite(d)]
    return float(min(1.0, 2 * (1 + min(np.sum(d <= 0), np.sum(d >= 0))) / (d.size + 1)))


def p_pos(draws):
    """C1.5: p_Delta = min(1, 2 (1 + #{Delta <= 0}) / (B_valid + 1))."""
    d = np.asarray(draws, float)
    d = d[np.isfinite(d)]
    return float(min(1.0, 2 * (1 + np.sum(d <= 0)) / (d.size + 1)))


def summ(est, draws):
    d = np.asarray(draws, float)
    fin = np.isfinite(d)
    out = {"est": float(est), "valid_frac": float(fin.mean()) if d.size else 0.0,
           "boot_se": float(np.std(d[fin], ddof=1)) if fin.sum() > 1 else math.nan,
           "p": p_two_sided(d), "draws": d}
    for tag, a in (("95", 0.05), ("975", 0.025), ("9833", ALPHA_LINE)):
        out["ci" + tag] = pct_ci(d, a)
    return out


def h_i2_call(s):
    return ADJ.classify_h_i2({"valid_frac": s["valid_frac"], "ci975": s["ci975"]})


def call_from_ci(valid_frac, lo, hi):
    if valid_frac >= MIN_VALID and math.isfinite(lo) and math.isfinite(hi):
        if lo > 0.0:
            return "MODEL_SUPERIOR_TO_FCT"
        if hi < 0.0:
            return "FCT_SUPERIOR_TO_MODEL"
    return "INDISTINGUISHABLE"


# ----------------------------------------------------------------------------------------------- splits
def split_c14(target_names, half_valid, s, seed):
    """C1.4: sort by target symbol; permute half_valid True then False indices; odd positions -> P1, even -> P2."""
    order = np.argsort(np.asarray(target_names), kind="stable")
    hv = np.asarray(half_valid, bool)[order]
    rng = np.random.default_rng(seed + s)
    perm_t = rng.permutation(order[hv])
    perm_f = rng.permutation(order[~hv])
    perm = np.concatenate([perm_t, perm_f]).astype(int)
    return np.sort(perm[1::2]), np.sort(perm[0::2])


def safe_argmax(v):
    v = np.where(np.isfinite(v), v, -np.inf)
    return int(np.argmax(v))


def safe_argmin(v):
    v = np.where(np.isfinite(v), v, np.inf)
    return int(np.argmin(v))


# ----------------------------------------------------------------------------------------------- L2
def criterion_summary(c):
    c = np.asarray(c, float)
    order = np.argsort(-c, kind="stable")  # descending; ties -> lowest index first
    t1, t2 = int(order[0]), int(order[1])
    return {"top1_idx": t1, "top2_idx": t2, "crit_top1": float(c[t1]), "crit_top2": float(c[t2]),
            "margin": float(c[t1] - c[t2]), "b0_crit": float(c[0]),
            "b0_rank": int(1 + np.sum(c > c[0])), "values": c}


def l2_select(ntc56, obs_ntc):
    pear = np.array([np.corrcoef(ntc56[b], obs_ntc)[0, 1] for b in range(ntc56.shape[0])])
    negl2 = -np.linalg.norm(ntc56 - obs_ntc[None, :], axis=1)
    return criterion_summary(pear), criterion_summary(negl2)


# ----------------------------------------------------------------------------------------------- J5b
def rankdata(x):
    x = np.asarray(x, float)
    order = np.argsort(x, kind="stable")
    r = np.empty(x.size)
    xs = x[order]
    i = 0
    while i < x.size:
        j = i
        while j + 1 < x.size and xs[j + 1] == xs[i]:
            j += 1
        r[order[i:j + 1]] = 0.5 * (i + j) + 1
        i = j + 1
    return r


def spearman_perm(a, b, perm_table):
    """Spearman rho and the §3.5 / C2.2 permutation p. Exact integer arithmetic when there are no ties."""
    ra, rb = rankdata(a), rankdata(b)
    n = ra.size
    ties = (np.unique(ra).size < n) or (np.unique(rb).size < n)
    za = ra - ra.mean()
    zb = rb - rb.mean()
    rho = float(za @ zb / math.sqrt((za @ za) * (zb @ zb)))
    if not ties:
        ia, ib = ra.astype(np.int64), rb.astype(np.int64)
        c = n * (n * n - 1)
        d_obs = abs(c - 6 * int(np.sum((ia - ib) ** 2)))
        dk = np.abs(c - 6 * np.sum((ia[None, :] - ib[perm_table]) ** 2, axis=1))
        cnt = int(np.sum(dk >= d_obs))
        n_tie = int(np.sum(dk == d_obs))   # permutations whose |rho| equals |rho_obs| exactly
        rho_exact = 1 - 6 * int(np.sum((ia - ib) ** 2)) / c
        rho = float(rho_exact)
    else:
        rk = (zb[perm_table] @ za) / math.sqrt((za @ za) * (zb @ zb))
        cnt = int(np.sum(np.abs(rk) >= abs(rho)))
        n_tie = int(np.sum(np.abs(np.abs(rk) - abs(rho)) <= 1e-12))
    return rho, (1 + cnt) / (perm_table.shape[0] + 1), ties, n_tie


def bh(pvals):
    p = np.asarray(pvals, float)
    m = p.size
    order = np.argsort(p, kind="stable")
    q = np.empty(m)
    prev = 1.0
    for k in range(m - 1, -1, -1):
        i = order[k]
        prev = min(prev, p[i] * m / (k + 1))
        q[i] = prev
    return q


def top_bottom(v, k=5):
    v = np.asarray(v, float)
    o = np.argsort(-v, kind="stable")
    return set(o[:k].tolist()), set(o[::-1][:k].tolist())


# ----------------------------------------------------------------------------------------------- one line x model
def analyze(L, obs, sw, cfg, g4=None):
    N, seed, nb = L.N, cfg["seed"], cfg["n_boot"]
    D56 = sw["delta_pred_56"]
    out = {"N": N}
    chk = {}
    # array-level consistency checks (schema §3.2)
    l1_pred = D56.mean(axis=0)
    chk["arr.L1_equals_mean56"] = float(np.max(np.abs(l1_pred - sw["delta_pred_L1_ens56"])))
    chk["arr.ens8_equals_mean8"] = float(np.max(np.abs(D56[list(ENS8)].mean(axis=0) - sw["delta_pred_multi_batch_ens8"])))
    med = int(sw["medoid_batch_idx"])
    chk["arr.medoid_equals_row"] = float(np.max(np.abs(D56[med] - sw["delta_pred_medoid_idx"])))

    st56 = np.stack([L.stats(D56[b]) for b in range(NB_IDX)])
    st = {"FCT": L.stats(obs["fct_delta_pred"]), "L1": L.stats(l1_pred),
          "bypass": L.stats(sw["delta_pred_bypass_zero"]), "meanemb": L.stats(sw["delta_pred_mean_batch_emb"]),
          "ens8": L.stats(D56[list(ENS8)].mean(axis=0)), "medoid": st56[med]}

    # L2 selection (§3.4, C2.5)
    pear, negl2 = l2_select(sw["ntc_pb_ensemble_56"], obs["obs_ntc_mean_official_median_log1p"])
    b_ntc = pear["top1_idx"]
    out["l2"], out["l2neg"] = pear, negl2
    out["l2_fold_idx"] = [l2_select(sw["ntc_pb_per_fold_56"][f], obs["obs_ntc_mean_official_median_log1p"])[0]["top1_idx"]
                          for f in range(sw["ntc_pb_per_fold_56"].shape[0])]
    chk["arr.L2_export_idx_equals_recomputed"] = int(int(sw["l2_selected_batch_idx"]) != b_ntc)
    chk["arr.L2_equals_row"] = float(np.max(np.abs(D56[int(sw["l2_selected_batch_idx"])] - sw["delta_pred_L2_ntc_recon_match"])))
    st["L2"] = st56[b_ntc]

    # point estimates and joint bootstrap (§3.1)
    one = Weights(L, np.ones(N))
    r_pt = np.array([one.r(st56[b])[0] for b in range(NB_IDX)])
    pt = {k: one.r(v)[0] for k, v in st.items()}
    out["rho_ab_orth"] = float(one.rab[0])
    WJ = Weights(L, boot_weights(N, nb, seed))
    R56 = np.stack([WJ.r(st56[b]) for b in range(NB_IDX)], axis=1)  # (B, 56)
    RB = {k: WJ.r(v) for k, v in st.items()}
    out["r56"] = [summ(r_pt[b], R56[:, b]) for b in range(NB_IDX)]
    out["r56_pt"] = r_pt
    for k in st:
        out[k] = summ(pt[k], RB[k])

    # landscape
    out["land"] = {"min": float(np.min(r_pt)), "max": float(np.max(r_pt)), "median": float(np.median(r_pt)),
                   "mean": float(np.mean(r_pt)), "q10": float(np.percentile(r_pt, 10)), "q25": float(np.percentile(r_pt, 25)),
                   "q75": float(np.percentile(r_pt, 75)), "q90": float(np.percentile(r_pt, 90)),
                   "iqr": float(np.percentile(r_pt, 75) - np.percentile(r_pt, 25)), "s_obs": float(np.std(r_pt, ddof=1)),
                   "range": float(np.max(r_pt) - np.min(r_pt)), "argmax": safe_argmax(r_pt), "argmin": safe_argmin(r_pt),
                   "b0_rank_desc": int(1 + np.sum(r_pt > r_pt[0]))}

    # variance component (§3.2.3, C1.6)
    fin_idx = np.isfinite(R56).mean(axis=0)
    joint = np.all(np.isfinite(R56), axis=1)
    Rv = R56[joint]
    s2_obs = float(np.var(r_pt, ddof=1))
    if Rv.shape[0] > 1:
        v_noise = float(np.sum(np.var(Rv - Rv.mean(axis=1, keepdims=True), axis=0, ddof=1)) / (NB_IDX - 1))
        Z = Rv - Rv.mean(axis=0, keepdims=True)
        s2_null = np.var(Z, axis=1, ddof=1)
        p_var = float((1 + np.sum(s2_null >= s2_obs)) / (Rv.shape[0] + 1))
    else:
        v_noise, p_var = math.nan, math.nan
    tau = math.sqrt(max(0.0, s2_obs - v_noise)) if math.isfinite(v_noise) else math.nan
    out["vc"] = {"s2_obs": s2_obs, "v_noise": v_noise, "tau": tau, "p_var": p_var,
                 "joint_valid_frac": float(joint.mean()), "min_index_valid_frac": float(fin_idx.min())}

    # cross-fitting on split s = 0 (§3.2.2, C1.4) with within-half bootstrap (P2 seed, P1 seed + 10000)
    P1, P2 = split_c14(obs["target_names"], obs["half_valid"], 0, seed)
    out["split"] = {"p1": P1.tolist(), "p2": P2.tolist(),
                    "names_sorted": bool(list(obs["target_names"]) == sorted(obs["target_names"].tolist()))}
    chk["arr.split_s0_p1_matches_C1.4"] = int(not np.array_equal(P1, np.asarray(obs["split_s0_p1_indices"])))
    chk["arr.split_s0_p2_matches_C1.4"] = int(not np.array_equal(P2, np.asarray(obs["split_s0_p2_indices"])))
    out["half_lt2"] = bool(min(P1.size, P2.size) < 2)   # C2.3
    cf = cross_fit(L, st56, st, P1, P2, nb, seed, b_ntc)
    out["cf"] = cf
    # multi-split diagnostic (secondary)
    ms = []
    for s in range(cfg["n_splits_diag"]):
        q1, q2 = split_c14(obs["target_names"], obs["half_valid"], s, seed)
        w1, w2 = Weights(L, indicator(N, q1)), Weights(L, indicator(N, q2))
        r1 = np.array([w1.r(st56[b])[0] for b in range(NB_IDX)])
        r2 = np.array([w2.r(st56[b])[0] for b in range(NB_IDX)])
        bs1, bd1, bs2, bd2 = safe_argmax(r1), safe_argmin(r1), safe_argmax(r2), safe_argmin(r2)
        ms.append(0.5 * ((r2[bs1] - r2[bd1]) + (r1[bs2] - r1[bd2])))
    ms = np.array(ms, float)
    msf = ms[np.isfinite(ms)]
    out["ms"] = ({"mean": float(msf.mean()), "median": float(np.median(msf)), "min": float(msf.min()),
                  "max": float(msf.max()), "iqr": float(np.percentile(msf, 75) - np.percentile(msf, 25))}
                 if msf.size else {})

    # J2 detectable-signal predicate
    l1s = out["L1"]
    sig_l1 = bool(l1s["valid_frac"] >= MIN_VALID and l1s["ci9833"][0] > 0.0)
    orc = cf["oracle"]
    sig_or = bool(orc["valid_frac"] >= MIN_VALID and orc["ci9833"][0] > 0.0)
    # C2.3: on a line whose smaller half has < 2 perturbations the cross-fitted oracle is undefined, so only the
    # L1_ens56 arm of the J2 predicate can qualify (controller reading; it feeds the descriptive H-J2 class there)
    out["sig"] = {"has_signal": bool(sig_l1 or (sig_or and not out["half_lt2"])), "l1": sig_l1, "oracle": sig_or,
                  "has_signal_with_oracle": bool(sig_l1 or sig_or)}

    # H-J1 validity (C1.6, §3.2.5.1, C2.3)
    valid = bool(fin_idx.min() >= MIN_VALID and joint.mean() >= MIN_VALID and math.isfinite(out["rho_ab_orth"])
                 and out["rho_ab_orth"] > 0 and cf["vf_half2"] >= MIN_VALID and cf["vf_half1"] >= MIN_VALID
                 and cf["dcf"]["valid_frac"] >= MIN_VALID and not out["half_lt2"])
    out["hj1_valid"] = valid
    out["p_hj1"] = float(max(p_var, cf["p_delta"])) if math.isfinite(p_var) else math.nan
    out["hj1_base"] = hj1_verdict(out, ALPHA_LINE, holm_rejected=None)

    # H-J2 (§3.3, C1.8)
    calls, calls_n, calls_e, vs = [], [], [], []
    WE = Weights(L, boot_weights(N, cfg["n_boot_bonf56"], seed))
    RE56 = np.stack([WE.r(st56[b]) for b in range(NB_IDX)], axis=1)
    REf = WE.r(st["FCT"])
    q = ALPHA_BONF56 / 2
    for b in range(NB_IDX):
        s = summ(r_pt[b] - pt["FCT"], R56[:, b] - RB["FCT"])
        c975 = h_i2_call(s)
        lo_n, hi_n = s["est"] - Z_BONF56 * s["boot_se"], s["est"] + Z_BONF56 * s["boot_se"]
        de = RE56[:, b] - REf
        fe = np.isfinite(de)
        lo_e, hi_e = (np.quantile(de[fe], [q, 1 - q]) if fe.any() else (math.nan, math.nan))
        calls.append(c975)
        calls_n.append(call_from_ci(s["valid_frac"], lo_n, hi_n))
        calls_e.append(call_from_ci(float(fe.mean()), float(lo_e), float(hi_e)))
        vs.append({"est": s["est"], "ci975": s["ci975"], "boot_se": s["boot_se"], "valid_frac": s["valid_frac"],
                   "norm": (lo_n, hi_n), "emp": (float(lo_e), float(hi_e))})
    out["hj2"] = {"vs": vs, "primary": hj2_class(calls, out["sig"]["has_signal"]),
                  "bonf_norm": hj2_class(calls_n, out["sig"]["has_signal"]),
                  "bonf_emp": hj2_class(calls_e, out["sig"]["has_signal"])}

    # H-J3 (§3.4, C1.5)
    for tag, key in (("hj3a", "L1"), ("hj3b", "L2")):
        d = summ(pt[key] - r_pt[0], RB[key] - R56[:, 0])
        tie = (tag == "hj3b" and b_ntc == 0)
        d["p"] = 1.0 if tie else d["p"]
        d["tie"] = tie
        d["base"] = hj3_verdict(d, ALPHA_LINE, None)
        d["vsfct"] = summ(pt[key] - pt["FCT"], RB[key] - RB["FCT"])
        d["vsfct_call"] = h_i2_call(d["vsfct"])
        d["gap"] = cf["gap_" + key]
        d["practical"] = bool(d["est"] > DELTA)
        out[tag] = d
    if g4 is not None:
        out["j5a"], chk_g4 = j5a(L, st56, g4, cf, one, WJ, out["half_lt2"])
        chk.update(chk_g4)
    return out, chk


def j5a(L, st56, g4, cf, one, WJ, half_lt2):
    """§3.5 J5a under r_orth_cal: G4 2x2 per mismatched source, forward (B_M, B_S) and reversed (B_S, B_M).
    Cross-fitted (b*1, b_dag1) on P2 (bootstrap inside P2, seed); auxiliary (b*2, b_dag2) on P1 (seed + 10000, C1.11);
    continuity (0, 17) and (0, 35) on all N (joint bootstrap, seed). C2.3: no cross-fitted pair when a half has < 2."""
    chk = {}
    pos = {int(b): i for i, b in enumerate(g4["g4_batch_indices"])}
    srcs = [str(x) for x in g4["mismatched_sources"]]
    MIS = g4["delta_pred_mismatched"]
    w1, w2, WB1, WB2 = cf["_w"]
    for k in ("b_star_1", "b_dag_1", "b_star_2", "b_dag_2"):
        chk[f"arr.g4_{k}_equals_recomputed"] = int(int(g4[k]) != cf[k]) if not half_lt2 else 0
    pairs = {"c17": (0, 17, one, WJ), "c35": (0, 35, one, WJ)}
    if not half_lt2:
        pairs["cf12"] = (cf["b_star_1"], cf["b_dag_1"], w2, WB2)
        pairs["cf21"] = (cf["b_star_2"], cf["b_dag_2"], w1, WB1)
    res = {}
    cache = {}

    def stS(si, b):
        if (si, b) not in cache:
            cache[(si, b)] = L.stats(MIS[si, pos[b]])
        return cache[(si, b)]

    for pn, (bm0, bs0, wpt, wbs) in pairs.items():
        miss = [b for b in (bm0, bs0) if b not in pos]
        chk[f"arr.g4_{pn}_indices_exported"] = int(bool(miss))
        if miss:
            continue
        for dname, (bm, bs) in (("fwd", (bm0, bs0)), ("rev", (bs0, bm0))):
            for si, src in enumerate(srcs):
                cells = {"r_M_BM": st56[bm], "r_M_BS": st56[bs], "r_S_BM": stS(si, bm), "r_S_BS": stS(si, bs)}
                ptv = {k: wpt.r(v)[0] for k, v in cells.items()}
                bsv = {k: wbs.r(v) for k, v in cells.items()}
                e = {k: summ(ptv[k], bsv[k]) for k in cells}
                for dk, (a, b) in (("d_batch", ("r_M_BM", "r_M_BS")), ("d_ctrl", ("r_M_BM", "r_S_BM")),
                                   ("d_both", ("r_M_BM", "r_S_BS"))):
                    e[dk] = summ(ptv[a] - ptv[b], bsv[a] - bsv[b])
                db, dc, dbo = e["d_batch"]["est"], e["d_ctrl"]["est"], e["d_both"]["est"]
                e["d_int"] = dbo - dc - db
                den = abs(dc) + abs(db)
                e["share"] = abs(db) / den if den > 0 else math.nan
                res[(pn, dname, src)] = e
    return res, chk


def cross_fit(L, st56, st, P1, P2, nb, seed, b_ntc):
    N = L.N
    w1, w2 = Weights(L, indicator(N, P1)), Weights(L, indicator(N, P2))
    r1 = np.array([w1.r(st56[b])[0] for b in range(NB_IDX)])
    r2 = np.array([w2.r(st56[b])[0] for b in range(NB_IDX)])
    bs1, bd1, bs2, bd2 = safe_argmax(r1), safe_argmin(r1), safe_argmax(r2), safe_argmin(r2)
    WB2 = Weights(L, boot_weights(N, nb, seed, P2))           # inside P2, seed
    WB1 = Weights(L, boot_weights(N, nb, seed + 10000, P1))   # inside P1, seed + 10000
    B2 = {b: WB2.r(st56[b]) for b in {bs1, bd1}}
    B1 = {b: WB1.r(st56[b]) for b in {bs2, bd2}}
    d12_pt, d21_pt = r2[bs1] - r2[bd1], r1[bs2] - r1[bd2]
    d12_bs, d21_bs = B2[bs1] - B2[bd1], B1[bs2] - B1[bd2]
    dcf_bs = 0.5 * (d12_bs + d21_bs)
    orc_pt = 0.5 * (r2[bs1] + r1[bs2])
    orc_bs = 0.5 * (B2[bs1] + B1[bs2])
    cf = {"b_star_1": bs1, "b_dag_1": bd1, "b_star_2": bs2, "b_dag_2": bd2,
          "rho_p1": float(w1.rab[0]), "rho_p2": float(w2.rab[0]),
          "d12": summ(d12_pt, d12_bs), "d21": summ(d21_pt, d21_bs),
          "dcf": summ(0.5 * (d12_pt + d21_pt), dcf_bs), "oracle": summ(orc_pt, orc_bs),
          "vf_half2": float(np.isfinite(d12_bs).mean()), "vf_half1": float(np.isfinite(d21_bs).mean()),
          "p_delta": p_pos(dcf_bs), "r1": r1, "r2": r2, "_w": (w1, w2, WB1, WB2)}
    for key, s in (("L1", st["L1"]), ("L2", st56[b_ntc])):
        g_pt = 0.5 * (w2.r(s)[0] + w1.r(s)[0]) - orc_pt
        g_bs = 0.5 * (WB2.r(s) + WB1.r(s)) - orc_bs
        cf["gap_" + key] = summ(g_pt, g_bs)
    return cf


def hj1_verdict(o, alpha, holm_rejected):
    """§3.2.5 with C1.7 subreasons, C2.1 Holm gate and C2.3 half-size rule. holm_rejected None = non-Holm line."""
    cf = o["cf"]
    if o["half_lt2"]:
        return {"verdict": "INCONCLUSIVE", "subreason": "VALIDITY_FAIL", "detail": "HALF_SIZE_LT_2"}
    lo, hi = pct_ci(cf["dcf"]["draws"], alpha)
    hi95 = cf["dcf"]["ci95"][1]
    if not o["hj1_valid"] or not (math.isfinite(lo) and math.isfinite(hi) and math.isfinite(hi95)):
        return {"verdict": "INCONCLUSIVE", "subreason": "VALIDITY_FAIL", "ci_lo": lo}
    p_var, tau, d = o["vc"]["p_var"], o["vc"]["tau"], cf["dcf"]["est"]
    if holm_rejected is None:
        detected = bool(p_var < alpha and lo > 0.0)
    else:
        detected = bool(holm_rejected and lo > 0.0)
    if detected:
        if d > DELTA and tau > TAU_THR:
            v = "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL"
        else:
            v = "BATCH_DEPENDENCE_DETECTED_SUB_DELTA"
        return {"verdict": v, "subreason": "", "ci_lo": lo, "detected": True}
    if not o["sig"]["has_signal"]:
        return {"verdict": "INCONCLUSIVE_LOW_SIGNAL", "subreason": "", "ci_lo": lo, "detected": False}
    if hi95 < DELTA and tau <= TAU_THR:
        return {"verdict": "BATCH_INVARIANT_WITHIN_DELTA", "subreason": "", "ci_lo": lo, "detected": False}
    return {"verdict": "INCONCLUSIVE", "subreason": "INSUFFICIENT_PRECISION", "ci_lo": lo, "detected": False}


def hj1_gate(o, alpha):
    if o["half_lt2"] or not o["hj1_valid"]:
        return False
    return bool(pct_ci(o["cf"]["dcf"]["draws"], alpha)[0] > 0.0)


def hj3_verdict(d, alpha, holm_rejected):
    lo, hi = pct_ci(d["draws"], alpha)
    if d["tie"]:
        return {"verdict": "EXACT_TIE_SELECTED_ZERO", "ci": (lo, hi)}
    ok = d["valid_frac"] >= MIN_VALID and math.isfinite(lo) and math.isfinite(hi)
    if ok and (holm_rejected is None or holm_rejected):
        if lo > 0.0:
            return {"verdict": "RULE_SUPERIOR_TO_ZERO", "ci": (lo, hi)}
        if hi < 0.0:
            return {"verdict": "ZERO_SUPERIOR_TO_RULE", "ci": (lo, hi)}
    return {"verdict": "INDISTINGUISHABLE", "ci": (lo, hi)}


def hj3_gate(d, alpha):
    if d["tie"]:
        return False
    lo, hi = pct_ci(d["draws"], alpha)
    return bool(d["valid_frac"] >= MIN_VALID and math.isfinite(lo) and math.isfinite(hi) and (lo > 0.0 or hi < 0.0))


def hj2_class(calls, has_signal):
    c0 = calls[0]
    flips = [b for b, c in enumerate(calls) if c != c0]
    f = len(flips) / NB_IDX
    nm = sum(c == "MODEL_SUPERIOR_TO_FCT" for c in calls)
    nf = sum(c == "FCT_SUPERIOR_TO_MODEL" for c in calls)
    rev = nm >= 1 and nf >= 1
    if rev:
        k = "CALL_UNSTABLE_SIGN_REVERSAL"
    elif len(flips) >= 5:
        k = "CALL_UNSTABLE_SUBSTANTIAL"
    elif len(flips) >= 1:
        k = "CALL_UNSTABLE_MINOR"
    else:
        k = "CALL_STABLE" if has_signal else "INCONCLUSIVE_LOW_SIGNAL"
    return {"calls": list(calls), "b0_call": c0, "n_flip": len(flips), "f_flip": f, "has_full_reversal": rev,
            "n_model": nm, "n_fct": nf, "n_ind": NB_IDX - nm - nf, "classification": k}


# ----------------------------------------------------------------------------------------------- Holm (C1.5, C2.1)
def holm(tests, order_keys=None):
    """tests: list of (key, p, gate(alpha) -> bool), in tie-break order. Step-down on raw p alone (C2.1).
    order_keys: optional explicit ordering (a party's tie-break); used only if it is a p-sorted permutation of the keys."""
    m = len(tests)
    order = sorted(range(m), key=lambda i: tests[i][1])  # stable: ties keep the listed order
    if order_keys is not None:
        pos = {t[0]: i for i, t in enumerate(tests)}
        if sorted(order_keys) == sorted(pos) and all(tests[pos[a]][1] <= tests[pos[b]][1]
                                                     for a, b in zip(order_keys, order_keys[1:])):
            order = [pos[k] for k in order_keys]
    res, stop, run_max = {}, None, 0.0
    for j, i in enumerate(order):
        key, p, gate = tests[i]
        a = ALPHA_FW / (m - j)
        run_max = max(run_max, min(1.0, (m - j) * p))
        holm_rej = stop is None and p <= a
        if not holm_rej and stop is None:
            stop = a
        eff = a if holm_rej else stop
        res[key] = {"rank": j + 1, "raw_p": p, "step_alpha": a, "holm_adj_p": run_max, "holm_rejected": bool(holm_rej),
                    "effective_alpha": eff}
    for key, p, gate in tests:
        r = res[key]
        r["gate"] = bool(gate(r["effective_alpha"]))
        r["rejected"] = bool(r["holm_rejected"] and r["gate"])
    return res


# ----------------------------------------------------------------------------------------------- driver
def load_line(d):
    obs = dict(np.load(d / "observed.npz", allow_pickle=False))
    sws = {}
    for m in MODELS:
        f = d / f"preds_{m}_sweep56.npz"
        if f.exists():
            sws[m] = dict(np.load(f, allow_pickle=False))
            fg = d / f"preds_{m}_g4_mismatched.npz"
            sws[m]["_g4"] = dict(np.load(fg, allow_pickle=False)) if fg.exists() else None
    return obs, sws


def recompute(results_dir, cfg, log=print):
    adir = Path(results_dir) / "arrays"
    lines = sorted(p.name for p in adir.iterdir() if (p / "observed.npz").exists())
    R = {"cfg": cfg, "lines": {}, "chk": {}}
    for ln in lines:
        obs, sws = load_line(adir / ln)
        L = R16.Line(obs)
        R["lines"][ln] = {"role": str(obs["role"]), "N": L.N, "models": {}}
        for m, sw in sws.items():
            chk_names = (list(map(str, sw["target_names"])) == list(map(str, obs["target_names"]))
                         and list(map(str, sw["gene_names"])) == list(map(str, obs["gene_names"])))
            o, chk = analyze(L, obs, sw, cfg, None if cfg.get("skip_j5a") else sw.get("_g4"))
            chk["arr.names_match_observed"] = int(not chk_names)
            R["lines"][ln]["models"][m] = o
            R["chk"][f"{ln}|{m}"] = chk
            log(f"  recomputed {ln} {m}: N={L.N} b_ntc={o['l2']['top1_idx']} H-J1(base)={o['hj1_base']['verdict']}")
    # confirmatory family: units present with role confirmatory, in PREREG_J order
    conf = [ln for ln in CONF_ORDER if ln in R["lines"] and R["lines"][ln]["role"] == "confirmatory"]
    conf += sorted(ln for ln in R["lines"] if R["lines"][ln]["role"] == "confirmatory" and ln not in conf)
    R["conf"] = conf
    apply_holm(R)
    if not cfg.get("skip_j5b"):
        R["j5b"] = j5b(R, cfg)
    return R


def holm_tests(R, m):
    t = []
    for ln in R["conf"]:
        o = R["lines"][ln]["models"][m]
        t.append((f"{m}|{ln}|H-J1", o["p_hj1"], lambda a, o=o: hj1_gate(o, a)))
        t.append((f"{m}|{ln}|H-J3a", o["hj3a"]["p"], lambda a, o=o: hj3_gate(o["hj3a"], a)))
        t.append((f"{m}|{ln}|H-J3b", o["hj3b"]["p"], lambda a, o=o: hj3_gate(o["hj3b"], a)))
    return t


def tie_material(tests, res):
    """Holm-rejected tie groups (equal raw p): does the C2.1 gate depend on where a test sits inside its group?"""
    out = []
    by_p = {}
    for key, p, gate in tests:
        by_p.setdefault(p, []).append((key, gate))
    for p, grp in by_p.items():
        if len(grp) < 2 or not res[grp[0][0]]["holm_rejected"]:
            continue
        alphas = sorted(res[k]["step_alpha"] for k, _ in grp)
        for k, gate in grp:
            if gate(alphas[0]) != gate(alphas[-1]):
                out.append({"key": k, "raw_p": p, "gate_at_min_alpha": gate(alphas[0]), "gate_at_max_alpha": gate(alphas[-1]),
                            "alpha_range": [alphas[0], alphas[-1]]})
    return out


def apply_holm(R, orders=None):
    """Separate-model Holm (STATE_SE, STATE_HVG) and pooled sensitivity. orders: {family: [keys]} tie-break override."""
    orders = orders or {}
    conf = R["conf"]
    R["holm"], R["cohort"], R["holm_tie_material"] = {}, {}, {}
    for m in MODELS:
        if not conf or not all(m in R["lines"][ln]["models"] for ln in conf):
            continue
        tests = holm_tests(R, m)
        h = holm(tests, orders.get(m))
        R["holm"][m] = h
        R["holm_tie_material"][m] = tie_material(tests, h)
        k = 0
        for ln in conf:
            o = R["lines"][ln]["models"][m]
            hj = h[f"{m}|{ln}|H-J1"]
            o["hj1_holm"] = hj1_verdict(o, hj["effective_alpha"], hj["holm_rejected"])
            o["hj1_holm"]["alpha_eff"] = hj["effective_alpha"]
            for tag, tn in (("hj3a", "H-J3a"), ("hj3b", "H-J3b")):
                ht = h[f"{m}|{ln}|{tn}"]
                o[tag]["holm"] = hj3_verdict(o[tag], ht["effective_alpha"], ht["holm_rejected"])
            k += o["hj1_holm"]["verdict"] == "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL"
        R["cohort"][m] = {"K_conf": len(conf), "k": k, "claim": bool(k >= 2 and len(conf) >= 2)}
    if all(m in R["holm"] for m in MODELS):
        tests = holm_tests(R, MODELS[0]) + holm_tests(R, MODELS[1])
        hp = holm(tests, orders.get("pooled"))
        R["holm_pooled"] = hp
        R["holm_tie_material"]["pooled"] = tie_material(tests, hp)
        R["cohort_pooled"] = {}
        for m in MODELS:
            k = 0
            for ln in conf:
                o = R["lines"][ln]["models"][m]
                hj = hp[f"{m}|{ln}|H-J1"]
                k += hj1_verdict(o, hj["effective_alpha"], hj["holm_rejected"])["verdict"] == "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL"
            R["cohort_pooled"][m] = {"k": k, "claim": bool(k >= 2 and len(conf) >= 2)}


def j5b(R, cfg):
    rng = np.random.default_rng(cfg["seed"])
    perm = np.stack([rng.permutation(NB_IDX) for _ in range(cfg["n_perm"])])  # C1.10 shared table
    out = {"by_model": {}, "cross": {}}
    lines = R["lines"]
    conf = R["conf"]
    fresh16 = sorted(ln for ln in lines if ln not in NOT_FRESH and lines[ln]["N"] >= 16)
    sets = {"conf": [(a, b) for i, a in enumerate(conf) for b in conf[i + 1:]],
            "pfizer": [PFIZER_PAIR] if all(x in lines for x in PFIZER_PAIR) else [],
            "fresh16": [(a, b) for i, a in enumerate(fresh16) for b in fresh16[i + 1:]]}
    for m in MODELS:
        out["by_model"][m] = {}
        for sname, pairs in sets.items():
            res = {}
            for a, b in pairs:
                if m not in lines[a]["models"] or m not in lines[b]["models"]:
                    continue
                va, vb = lines[a]["models"][m]["r56_pt"], lines[b]["models"][m]["r56_pt"]
                rho, p, ties, n_tie = spearman_perm(va, vb, perm)
                ta, ba = top_bottom(va)
                tb, bb = top_bottom(vb)
                res[tuple(sorted((a, b)))] = {"rho": rho, "p_perm": p, "ties": ties, "n_tie": n_tie,
                                              "top5": len(ta & tb), "bot5": len(ba & bb)}
            if res:
                q = bh([v["p_perm"] for v in res.values()])
                for v, qq in zip(res.values(), q):
                    v["bh"] = float(qq)
            out["by_model"][m][sname] = res
    res = {}
    for ln in sorted(lines):
        if all(m in lines[ln]["models"] for m in MODELS):
            va, vb = lines[ln]["models"]["STATE_SE"]["r56_pt"], lines[ln]["models"]["STATE_HVG"]["r56_pt"]
            rho, p, ties, n_tie = spearman_perm(va, vb, perm)
            res[ln] = {"rho": rho, "p_perm": p, "ties": ties, "n_tie": n_tie}
    if res:
        q = bh([v["p_perm"] for v in res.values()])
        for v, qq in zip(res.values(), q):
            v["bh"] = float(qq)
    out["cross"] = res
    return out


# ----------------------------------------------------------------------------------------------- canonical flat view
def flat_controller(R):
    F = {}
    for ln, LD in R["lines"].items():
        for m, o in LD["models"].items():
            p = f"{ln}|{m}|"
            F[p + "rho_ab_orth"] = o["rho_ab_orth"]
            for b in range(NB_IDX):
                s = o["r56"][b]
                F[p + f"r56.est[{b}]"] = s["est"]
                F[p + f"r56.ci975_lo[{b}]"], F[p + f"r56.ci975_hi[{b}]"] = s["ci975"]
                F[p + f"r56.boot_se[{b}]"] = s["boot_se"]
            for k, lab in (("FCT", "fct"), ("L1", "l1"), ("bypass", "sec.bypass_zero"), ("meanemb", "sec.mean_batch_emb"),
                           ("ens8", "sec.multi_batch_ens8"), ("medoid", "sec.medoid_idx")):
                s = o[k]
                F[p + lab + ".est"] = s["est"]
                F[p + lab + ".ci975_lo"], F[p + lab + ".ci975_hi"] = s["ci975"]
                F[p + lab + ".ci9833_lo"], F[p + lab + ".ci9833_hi"] = s["ci9833"]
            for k in ("min", "max", "median", "mean", "q10", "q25", "q75", "q90", "iqr", "s_obs", "range", "argmax", "argmin",
                      "b0_rank_desc"):
                F[p + "land." + k] = o["land"][k]
            for tag in ("l2", "l2neg"):
                for k in ("top1_idx", "top2_idx", "crit_top1", "crit_top2", "margin", "b0_crit", "b0_rank"):
                    F[p + f"{tag}.{k}"] = o[tag][k]
            F[p + "l2.fold_idx"] = list(o["l2_fold_idx"])
            for k in ("s2_obs", "v_noise", "tau", "p_var"):
                F[p + "vc." + k] = o["vc"][k]
            F[p + "split.p1"], F[p + "split.p2"] = o["split"]["p1"], o["split"]["p2"]
            F[p + "has_signal"] = o["sig"]["has_signal"]
            F[p + "hj1.base.verdict"] = o["hj1_base"]["verdict"]
            if o["hj1_base"]["verdict"] == "INCONCLUSIVE":   # C1.7: subreasons belong to INCONCLUSIVE only
                F[p + "hj1.base.subreason"] = o["hj1_base"]["subreason"]
            F[p + "vc.p_var"] = o["vc"]["p_var"]
            cf = o["cf"]
            if not o["half_lt2"]:   # C2.3: cross-fitted quantities undefined on such lines; skip them
                for k in ("b_star_1", "b_dag_1", "b_star_2", "b_dag_2"):
                    F[p + "cf." + k] = cf[k]
                for k in ("d12", "d21"):
                    F[p + f"cf.{k}.est"] = cf[k]["est"]
                for k in ("dcf", "oracle"):
                    F[p + f"cf.{k}.est"] = cf[k]["est"]
                    F[p + f"cf.{k}.ci9833_lo"], F[p + f"cf.{k}.ci9833_hi"] = cf[k]["ci9833"]
                    F[p + f"cf.{k}.ci95_lo"], F[p + f"cf.{k}.ci95_hi"] = cf[k]["ci95"]
                    F[p + f"cf.{k}.boot_se"] = cf[k]["boot_se"]
                F[p + "cf.p_delta"] = cf["p_delta"]
                F[p + "hj1.p_iut"] = o["p_hj1"]
                F[p + "hj1.base.ci_adj_lo"] = o["hj1_base"].get("ci_lo", math.nan)
                for k in ("L1", "L2"):
                    g = cf["gap_" + k]
                    F[p + f"gap.{k}.est"] = g["est"]
                    F[p + f"gap.{k}.ci9833_lo"], F[p + f"gap.{k}.ci9833_hi"] = g["ci9833"]
                for k, v in o["ms"].items():
                    F[p + "ms." + k] = v
            h2 = o["hj2"]
            for tag in ("primary", "bonf_norm", "bonf_emp"):
                for k in ("calls", "b0_call", "n_flip", "f_flip", "has_full_reversal", "n_model", "n_fct", "n_ind",
                          "classification"):
                    F[p + f"hj2.{tag}.{k}"] = h2[tag][k]
            for b in range(NB_IDX):
                v = h2["vs"][b]
                F[p + f"hj2.vsfct.est[{b}]"] = v["est"]
                F[p + f"hj2.vsfct.ci975_lo[{b}]"], F[p + f"hj2.vsfct.ci975_hi[{b}]"] = v["ci975"]
                F[p + f"hj2.vsfct.boot_se[{b}]"] = v["boot_se"]
                F[p + f"hj2.vsfct.bonf_norm_lo[{b}]"], F[p + f"hj2.vsfct.bonf_norm_hi[{b}]"] = v["norm"]
                F[p + f"hj2.vsfct.bonf_emp_lo[{b}]"], F[p + f"hj2.vsfct.bonf_emp_hi[{b}]"] = v["emp"]
            for tag in ("hj3a", "hj3b"):
                d = o[tag]
                F[p + f"{tag}.est"] = d["est"]
                F[p + f"{tag}.ci9833_lo"], F[p + f"{tag}.ci9833_hi"] = d["ci9833"]
                F[p + f"{tag}.boot_se"] = d["boot_se"]
                F[p + f"{tag}.p"] = d["p"]
                F[p + f"{tag}.base.verdict"] = d["base"]["verdict"]
                F[p + f"{tag}.practical_gt_delta"] = d["practical"]
                F[p + f"{tag}.vsfct.est"] = d["vsfct"]["est"]
                F[p + f"{tag}.vsfct.ci975_lo"], F[p + f"{tag}.vsfct.ci975_hi"] = d["vsfct"]["ci975"]
                F[p + f"{tag}.vsfct.call"] = d["vsfct_call"]
                if "holm" in d:
                    F[p + f"{tag}.holm.verdict"] = d["holm"]["verdict"]
                    F[p + f"{tag}.holm.ci_lo"], F[p + f"{tag}.holm.ci_hi"] = d["holm"]["ci"]
            for (pn, dname, src), e in (o.get("j5a") or {}).items():
                q = p + f"j5a.{pn}.{dname}.{src}."
                for k in ("r_M_BM", "r_M_BS", "r_S_BM", "r_S_BS", "d_batch", "d_ctrl", "d_both"):
                    F[q + k + ".est"] = e[k]["est"]
                    F[q + k + ".ci975_lo"], F[q + k + ".ci975_hi"] = e[k]["ci975"]
                F[q + "d_int"], F[q + "share"] = e["d_int"], e["share"]
            if o.get("half_lt2"):
                F[p + "j5a.cf_pair_reported"] = False   # C2.3: no J5a cross-fitted pair is reported
            if "hj1_holm" in o:
                F[p + "hj1.holm.verdict"] = o["hj1_holm"]["verdict"]
                if o["hj1_holm"]["verdict"] == "INCONCLUSIVE":
                    F[p + "hj1.holm.subreason"] = o["hj1_holm"]["subreason"]
                F[p + "hj1.holm.alpha_eff"] = o["hj1_holm"]["alpha_eff"]
                F[p + "hj1.holm.ci_lo"] = o["hj1_holm"].get("ci_lo", math.nan)
    for fam, H in [("holm", R["holm"][m]) for m in R["holm"]] + ([("pholm", R["holm_pooled"])] if "holm_pooled" in R else []):
        for key, r in H.items():
            for k in ("rank", "raw_p", "step_alpha", "effective_alpha", "holm_adj_p", "holm_rejected", "rejected"):
                F[f"{fam}|{key}|{k}"] = r[k]
    for m, c in R["cohort"].items():
        F[f"cohort|{m}|K_conf"], F[f"cohort|{m}|k"], F[f"cohort|{m}|claim"] = c["K_conf"], c["k"], c["claim"]
    for m, c in R.get("cohort_pooled", {}).items():
        F[f"pcohort|{m}|k"], F[f"pcohort|{m}|claim"] = c["k"], c["claim"]
    if "j5b" in R:
        for m, sets in R["j5b"]["by_model"].items():
            for sname, res in sets.items():
                for (a, b), v in res.items():
                    for k in ("rho", "p_perm", "bh", "top5", "bot5"):
                        F[f"j5b|{m}|{sname}|{a}__vs__{b}|{k}"] = v[k]
        for ln, v in R["j5b"]["cross"].items():
            for k in ("rho", "p_perm", "bh"):
                F[f"j5b|cross|{ln}|{k}"] = v[k]
    return F


# ----------------------------------------------------------------------------------------------- party adapters
def _g(d, *path):
    for k in path:
        if d is None:
            return None
        if isinstance(d, dict):
            d = d.get(k)
        elif isinstance(d, list) and isinstance(k, int) and k < len(d):
            d = d[k]
        else:
            return None
    return d


def _put(P, key, v):
    if v is not None:
        P[key] = v


def flat_main(S):
    P = {}
    for ln, LD in S.get("lines", {}).items():
        for m, o in LD.get("models", {}).items():
            p = f"{ln}|{m}|"
            land = o.get("landscape_r_orth_cal", {})
            _put(P, p + "rho_ab_orth", land.get("rho_ab_orth"))
            for b, s in enumerate(land.get("per_index", [])):
                _put(P, p + f"r56.est[{b}]", s.get("estimate"))
                _put(P, p + f"r56.ci975_lo[{b}]", s.get("ci975_lo"))
                _put(P, p + f"r56.ci975_hi[{b}]", s.get("ci975_hi"))
                _put(P, p + f"r56.boot_se[{b}]", s.get("boot_se"))
            for k, kk in (("min", "min"), ("max", "max"), ("median", "median"), ("mean", "mean"), ("q10", "q10"), ("q25", "q25"),
                          ("q75", "q75"), ("q90", "q90"), ("iqr", "iqr"), ("s_obs", "observed_sd"), ("range", "raw_range"),
                          ("argmax", "argmax_idx"), ("argmin", "argmin_idx"), ("b0_rank_desc", "b0_rank_1based_desc")):
                _put(P, p + "land." + k, land.get(kk))
            for lab, src in (("fct", o.get("fct_r_orth_cal")), ("l1", o.get("l1_ens56_r_orth_cal"))):
                if src:
                    for k in ("est", "ci975_lo", "ci975_hi", "ci9833_lo", "ci9833_hi"):
                        _put(P, p + f"{lab}.{k}", src.get("estimate" if k == "est" else k))
            for lab in ("bypass_zero", "mean_batch_emb", "multi_batch_ens8", "medoid_idx"):
                src = _g(o, "secondary_neutral_comparators", lab, "r_orth_cal")
                if src:
                    for k in ("est", "ci975_lo", "ci975_hi", "ci9833_lo", "ci9833_hi"):
                        _put(P, p + f"sec.{lab}.{k}", src.get("estimate" if k == "est" else k))
            l2 = o.get("l2_ntc_selection", {})
            for tag, src in (("l2", l2.get("primary_pearson")), ("l2neg", l2.get("secondary_neg_l2"))):
                if src:
                    for k, kk in (("top1_idx", "top1_idx"), ("top2_idx", "top2_idx"), ("crit_top1", "criterion_top1"),
                                  ("crit_top2", "criterion_top2"), ("margin", "criterion_margin_top1_minus_top2"),
                                  ("b0_crit", "b0_criterion"), ("b0_rank", "b0_rank_1based")):
                        _put(P, p + f"{tag}.{k}", src.get(kk))
            pf = _g(l2, "secondary_per_fold_agreement", "per_fold_selected_idx")
            if pf:
                _put(P, p + "l2.fold_idx", list(pf.values()))
            vc = o.get("variance_component", {})
            for k, kk in (("s2_obs", "s_obs_sq"), ("v_noise", "v_noise"), ("tau", "tau_batch"), ("p_var", "p_var")):
                _put(P, p + "vc." + k, vc.get(kk))
            sp = _g(o, "cross_fitting", "split_s0") or {}
            _put(P, p + "split.p1", sp.get("p1_indices"))
            _put(P, p + "split.p2", sp.get("p2_indices"))
            for k in ("b_star_1", "b_dag_1", "b_star_2", "b_dag_2"):
                _put(P, p + "cf." + k, sp.get(k))
            for k, kk in (("d12", "delta_1to2"), ("d21", "delta_2to1")):
                _put(P, p + f"cf.{k}.est", _g(sp, kk, "estimate"))
            for k, kk in (("dcf", "delta_cf_0"), ("oracle", "r_cf_oracle_0")):
                s = sp.get(kk) or {}
                for f, ff in (("est", "estimate"), ("ci9833_lo", "ci9833_lo"), ("ci9833_hi", "ci9833_hi"), ("ci95_lo", "ci95_lo"),
                              ("ci95_hi", "ci95_hi"), ("boot_se", "boot_se")):
                    _put(P, p + f"cf.{k}.{f}", s.get(ff))
            _put(P, p + "cf.p_delta", _g(sp, "delta_cf_0", "p_pos_two_sided"))
            for k, kk in (("L1", "l1_vs_cf_oracle"), ("L2", "l2_vs_cf_oracle")):
                s = sp.get(kk) or {}
                for f, ff in (("est", "estimate"), ("ci9833_lo", "ci9833_lo"), ("ci9833_hi", "ci9833_hi")):
                    _put(P, p + f"gap.{k}.{f}", s.get(ff))
            ms = _g(o, "cross_fitting", "multi_split_diag") or {}
            for k, kk in (("mean", "delta_cf_mean"), ("median", "delta_cf_median"), ("min", "delta_cf_min"),
                          ("max", "delta_cf_max"), ("iqr", "delta_cf_iqr")):
                _put(P, p + "ms." + k, ms.get(kk))
            _put(P, p + "has_signal", _g(o, "detectable_signal_j2", "has_detectable_signal"))
            hb = o.get("H_J1_line_base", {})
            _put(P, p + "hj1.base.verdict", hb.get("verdict"))
            _put(P, p + "hj1.base.subreason", hb.get("subreason"))
            _put(P, p + "hj1.p_iut", hb.get("p_h_j1_iut"))
            _put(P, p + "hj1.base.ci_adj_lo", hb.get("delta_cf_0_ci_adj_lo"))
            h2 = o.get("H_J2", {})
            for tag, src in (("primary", h2), ("bonf_norm", h2.get("companion_bonf56_normal")),
                             ("bonf_emp", h2.get("companion_bonf56_emp20k"))):
                if not src:
                    continue
                cnt = src.get("counts", {})
                for k, v in (("calls", src.get("calls_56")), ("b0_call", src.get("baseline_b0_call")), ("n_flip", src.get("n_flip")),
                             ("f_flip", src.get("f_flip")), ("has_full_reversal", src.get("has_full_reversal")),
                             ("n_model", cnt.get("MODEL_SUPERIOR_TO_FCT")), ("n_fct", cnt.get("FCT_SUPERIOR_TO_MODEL")),
                             ("n_ind", cnt.get("INDISTINGUISHABLE")), ("classification", src.get("classification"))):
                    _put(P, p + f"hj2.{tag}.{k}", v)
            for v in h2.get("per_index_vs_fct", []):
                b = v["batch_idx"]
                for k, kk in (("est", "estimate"), ("ci975_lo", "ci975_lo"), ("ci975_hi", "ci975_hi"), ("boot_se", "boot_se"),
                              ("bonf_norm_lo", "bonf56_norm_lo"), ("bonf_norm_hi", "bonf56_norm_hi"),
                              ("bonf_emp_lo", "bonf56_emp20k_lo"), ("bonf_emp_hi", "bonf56_emp20k_hi")):
                    _put(P, p + f"hj2.vsfct.{k}[{b}]", v.get(kk))
            g4 = o.get("secondary_j5a_g4") or {}
            for pn, gk, pre in (("cf12", "cross_fitted_primary_p1_to_p2", "cf_p1_to_p2"),
                                ("cf21", "cross_fitted_auxiliary_p2_to_p1", "cf_p2_to_p1"),
                                ("c17", "continuity_b0_vs_b17", "continuity_b0_vs_b17"),
                                ("c35", "continuity_b0_vs_b35", "continuity_b0_vs_b35")):
                ps = _g(g4, gk, "per_source") or {}
                if pn == "cf12":
                    _put(P, p + "j5a.cf_pair_reported", bool(ps))
                for name, v in ps.items():
                    for dname in ("fwd", "rev"):
                        tag = f"{pre}_{dname}_mis_"
                        if name.startswith(tag):
                            src = name[len(tag):]
                            ro = v.get("r_orth_cal") or {}
                            q = p + f"j5a.{pn}.{dname}.{src}."
                            for k, kk in (("r_M_BM", "r_M_BM"), ("r_M_BS", "r_M_BS"), ("r_S_BM", "r_S_BM"), ("r_S_BS", "r_S_BS"),
                                          ("d_batch", "delta_batch"), ("d_ctrl", "delta_ctrl"), ("d_both", "delta_both")):
                                x = ro.get(kk) or {}
                                for f, ff in (("est", "estimate"), ("ci975_lo", "ci975_lo"), ("ci975_hi", "ci975_hi")):
                                    _put(P, q + f"{k}.{f}", x.get(ff))
                            _put(P, q + "d_int", ro.get("delta_interaction_point"))
                            _put(P, q + "share", ro.get("batch_share_point"))
            for tag, kk in (("hj3a", "H_J3a_L1_ens56"), ("hj3b", "H_J3b_L2_ntc_recon_match")):
                s = _g(o, "H_J3_line_base", kk) or {}
                vb = s.get("vs_b0") or {}
                for f, ff in (("est", "estimate"), ("ci9833_lo", "ci9833_lo"), ("ci9833_hi", "ci9833_hi"), ("boot_se", "boot_se"),
                              ("p", "p_two_sided")):
                    _put(P, p + f"{tag}.{f}", vb.get(ff))
                _put(P, p + f"{tag}.base.verdict", s.get("verdict"))
                ma = s.get("magnitude_annotation")
                if ma is not None and s.get("verdict") in ("RULE_SUPERIOR_TO_ZERO",):
                    _put(P, p + f"{tag}.practical_gt_delta", ma == "PRACTICAL_GT_DELTA")
                vf = s.get("vs_fct") or {}
                for f, ff in (("est", "estimate"), ("ci975_lo", "ci975_lo"), ("ci975_hi", "ci975_hi")):
                    _put(P, p + f"{tag}.vsfct.{f}", vf.get(ff))
                _put(P, p + f"{tag}.vsfct.call", _g(s, "vs_fct_h_i2_call", "classification"))
    sep = S.get("separate_model_holm_fwer05", {})
    for m, tab in sep.get("holm_tables", {}).items():
        for key, r in tab.items():
            ck = key.replace("::", "|")
            for k in ("rank_1based", "raw_p", "step_alpha", "effective_alpha", "holm_adj_p", "rejected"):
                _put(P, f"holm|{ck}|{'rank' if k == 'rank_1based' else k}", r.get(k))
    for m, c in sep.get("cohort_adjudications", {}).items():
        _put(P, f"cohort|{m}|K_conf", c.get("K_conf"))
        _put(P, f"cohort|{m}|k", c.get("k_batch_dependence_practically_meaningful"))
        _put(P, f"cohort|{m}|claim", c.get("general_cross_context_claim_permitted"))
        for ln, e in (c.get("per_confirmatory_line") or {}).items():
            p = f"{ln}|{m}|"
            h = e.get("H_J1_holm") or {}
            _put(P, p + "hj1.holm.verdict", h.get("verdict"))
            _put(P, p + "hj1.holm.subreason", h.get("subreason"))
            _put(P, p + "hj1.holm.alpha_eff", h.get("alpha_adj"))
            _put(P, p + "hj1.holm.ci_lo", h.get("delta_cf_0_ci_adj_lo"))
            for tag, kk in (("hj3a", "H_J3a_L1_ens56_holm"), ("hj3b", "H_J3b_L2_ntc_recon_match_holm")):
                s = e.get(kk) or {}
                _put(P, p + f"{tag}.holm.verdict", s.get("verdict"))
                _put(P, p + f"{tag}.holm.ci_lo", _g(s, "vs_b0", "ci_adj_lo"))
                _put(P, p + f"{tag}.holm.ci_hi", _g(s, "vs_b0", "ci_adj_hi"))
    pool = S.get("pooled_model_holm_sensitivity_fwer05", {})
    for key, r in (pool.get("holm_table") or {}).items():
        ck = key.replace("::", "|")
        for k in ("rank_1based", "raw_p", "step_alpha", "effective_alpha", "holm_adj_p", "rejected"):
            _put(P, f"pholm|{ck}|{'rank' if k == 'rank_1based' else k}", r.get(k))
    for m, c in (pool.get("cohort_adjudications") or {}).items():
        if isinstance(c, dict):
            _put(P, f"pcohort|{m}|k", c.get("k_batch_dependence_practically_meaningful"))
            _put(P, f"pcohort|{m}|claim", c.get("general_cross_context_claim_permitted"))
    j = S.get("j5b_cross_line_rank_correlations", {})
    smap = {"confirmatory_pairs": "conf", "pfizer_within_family_pair": "pfizer", "all_fresh_n_ge_16_pairs": "fresh16"}
    for m, sets in (j.get("by_model") or {}).items():
        for sname, sv in sets.items():
            if sname not in smap or not isinstance(sv, dict):
                continue
            for pk, v in (sv.get("pairs") or {}).items():
                a, b = sorted(pk.split("__vs__"))
                base = f"j5b|{m}|{smap[sname]}|{a}__vs__{b}|"
                _put(P, base + "rho", v.get("spearman_rho"))
                _put(P, base + "p_perm", v.get("p_perm_two_sided"))
                _put(P, base + "bh", v.get("bh_fdr_q"))
                _put(P, base + "top5", v.get("top5_overlap_count"))
                _put(P, base + "bot5", v.get("bottom5_overlap_count"))
    cm = j.get("cross_model_within_line") or {}
    cmp_ = cm.get("pairs", cm) if isinstance(cm, dict) else {}
    for pk, v in (cmp_ or {}).items():
        if isinstance(v, dict) and "spearman_rho" in v:
            ln2 = v.get("unit_a") or pk.split("__")[0]
            _put(P, f"j5b|cross|{ln2}|rho", v.get("spearman_rho"))
            _put(P, f"j5b|cross|{ln2}|p_perm", v.get("p_perm_two_sided"))
            _put(P, f"j5b|cross|{ln2}|bh", v.get("bh_fdr_q"))
    return P


def flat_at3(S):
    P = {}
    for m, lines in (S.get("models") or {}).items():
        for ln, o in lines.items():
            p = f"{ln}|{m}|"
            h1 = o.get("h_j1", {})
            _put(P, p + "rho_ab_orth", h1.get("rho_ab_orth"))
            for b, v in enumerate(h1.get("r_obs_56") or []):
                _put(P, p + f"r56.est[{b}]", v)
            for b, v in enumerate(h1.get("r_boot_se_56") or []):
                _put(P, p + f"r56.boot_se[{b}]", v)
            land = h1.get("landscape", {})
            for k, kk in (("min", "min"), ("max", "max"), ("median", "median"), ("mean", "mean"), ("q10", "q10"), ("q25", "q25"),
                          ("q75", "q75"), ("q90", "q90"), ("iqr", "iqr"), ("s_obs", "s_obs"), ("range", "range_obs"),
                          ("b0_rank_desc", "b0_rank_desc")):
                _put(P, p + "land." + k, land.get(kk))
            l1 = land.get("l1_ens56") or {}
            for k, kk in (("est", "estimate"), ("ci9833_lo", "ci9833_lo"), ("ci9833_hi", "ci9833_hi")):
                _put(P, p + "l1." + k, l1.get(kk))
            for k, kk in (("s2_obs", "s2_obs"), ("v_noise", "v_noise"), ("tau", "tau_batch"), ("p_var", "p_var")):
                _put(P, p + "vc." + k, h1.get(kk))
            _put(P, p + "cf.p_delta", h1.get("p_delta"))
            _put(P, p + "hj1.p_iut", h1.get("p_hj1_iu"))
            cfs = h1.get("cross_fit_s0") or {}
            for k in ("b_star_1", "b_dag_1", "b_star_2", "b_dag_2"):
                _put(P, p + "cf." + k, cfs.get(k))
            _put(P, p + "cf.d12.est", cfs.get("delta_1_to_2"))
            _put(P, p + "cf.d21.est", cfs.get("delta_2_to_1"))
            _put(P, p + "cf.dcf.est", cfs.get("delta_cf_0"))
            _put(P, p + "cf.dcf.ci95_lo", h1.get("ci95_lo"))
            _put(P, p + "cf.dcf.ci95_hi", h1.get("ci95_hi"))
            _put(P, p + "cf.dcf.boot_se", h1.get("delta_cf_boot_se"))
            _put(P, p + "hj1.base.ci_adj_lo", h1.get("ci_adj_lo"))
            orc = land.get("cf_oracle_s0") or {}
            for k, kk in (("est", "estimate"), ("ci9833_lo", "ci9833_lo"), ("ci9833_hi", "ci9833_hi"), ("boot_se", "boot_se")):
                _put(P, p + "cf.oracle." + k, orc.get(kk))
            ms = h1.get("multi_split_diag") or {}
            for k in ("mean", "median", "min", "max", "iqr"):
                _put(P, p + "ms." + k, ms.get(k))
            _put(P, p + "has_signal", h1.get("has_detectable_signal"))
            _put(P, p + "hj1.base.verdict", h1.get("classification"))
            _put(P, p + "hj1.base.subreason", h1.get("subreason") or "")
            h2 = o.get("h_j2", {})
            _put(P, p + "fct.est", h2.get("fct_r_orth_cal"))
            for tag, kk in (("primary", "primary_975"), ("bonf_norm", "bonf56_normal"), ("bonf_emp", "bonf56_empirical_20k")):
                src = h2.get(kk)
                if not src:
                    continue
                for k, k2 in (("calls", "calls_56"), ("b0_call", "b0_call"), ("n_flip", "n_flip"), ("f_flip", "f_flip"),
                              ("has_full_reversal", "has_full_reversal"), ("n_model", "n_MODEL_SUPERIOR"),
                              ("n_fct", "n_FCT_SUPERIOR"), ("n_ind", "n_INDISTINGUISHABLE")):
                    _put(P, p + f"hj2.{tag}.{k}", src.get(k2))
            _put(P, p + "hj2.primary.classification", h2.get("classification"))
            for b, v in enumerate(h2.get("delta_vs_fct_56") or []):
                _put(P, p + f"hj2.vsfct.est[{b}]", v)
            for b, v in enumerate(h2.get("ci975_56") or []):
                _put(P, p + f"hj2.vsfct.ci975_lo[{b}]", v[0])
                _put(P, p + f"hj2.vsfct.ci975_hi[{b}]", v[1])
            h3 = o.get("h_j3", {})
            for tag, kk in (("hj3a", "h_j3a_L1_ens56"), ("hj3b", "h_j3b_L2_ntc_recon_match")):
                s = h3.get(kk) or {}
                for f, ff in (("est", "estimate"), ("ci9833_lo", "ci_adj_lo"), ("ci9833_hi", "ci_adj_hi"), ("boot_se", "boot_se"),
                              ("p", "p_value")):
                    _put(P, p + f"{tag}.{f}", s.get(ff))
                _put(P, p + f"{tag}.base.verdict", s.get("classification"))
                if s.get("classification") == "RULE_SUPERIOR_TO_ZERO":
                    _put(P, p + f"{tag}.practical_gt_delta", s.get("practical_gt_0_04"))
                vf = s.get("vs_fct") or {}
                for f, ff in (("est", "estimate"), ("ci975_lo", "ci975_lo"), ("ci975_hi", "ci975_hi"), ("call", "call")):
                    _put(P, p + f"{tag}.vsfct.{f}", vf.get(ff))
                g = s.get("vs_cf_oracle_s0") or {}
                for f, ff in (("est", "estimate"), ("ci9833_lo", "ci9833_lo"), ("ci9833_hi", "ci9833_hi")):
                    _put(P, p + f"gap.{'L1' if tag == 'hj3a' else 'L2'}.{f}", g.get(ff))
            l2 = h3.get("l2_selection_diagnostics") or {}
            for k, kk in (("top1_idx", "top1_idx"), ("top2_idx", "top2_idx"), ("crit_top1", "criterion_top1"),
                          ("crit_top2", "criterion_top2"), ("margin", "criterion_margin_top1_minus_top2"),
                          ("b0_crit", "b0_criterion"), ("b0_rank", "b0_rank_desc")):
                _put(P, p + f"l2.{k}", l2.get(kk))
            neg = l2.get("neg_l2_diag") or {}
            for k, kk in (("top1_idx", "top1_idx"), ("top2_idx", "top2_idx"), ("crit_top1", "criterion_top1"),
                          ("crit_top2", "criterion_top2"), ("margin", "criterion_margin_top1_minus_top2")):
                _put(P, p + f"l2neg.{k}", neg.get(kk))
            _put(P, p + "l2.fold_idx", l2.get("per_fold_top1_indices"))
            for lab in ("bypass_zero", "mean_batch_emb", "multi_batch_ens8", "medoid_idx"):
                s = _g(o, "secondary_comparators", lab) or {}
                for f, ff in (("est", "estimate"), ("ci975_lo", "ci975_lo"), ("ci975_hi", "ci975_hi")):
                    _put(P, p + f"sec.{lab}.{f}", s.get(ff))
            g4 = o.get("j5a_g4") or {}
            for pn, gk in (("cf12", "cross_fit_primary_on_p2"), ("cf21", "cross_fit_auxiliary_on_p1"),
                           ("c17", "continuity_0_17"), ("c35", "continuity_0_35")):
                gg = g4.get(gk) or {}
                if pn == "cf12":
                    _put(P, p + "j5a.cf_pair_reported", bool(gg.get("per_source")))
                for dname, sk in (("fwd", "per_source"), ("rev", "reversed_symmetry")):
                    for src, v in (gg.get(sk) or {}).items():
                        if not isinstance(v, dict):
                            continue
                        q = p + f"j5a.{pn}.{dname}.{src}."
                        for k, kk in (("d_batch", "delta_batch_orth"), ("d_ctrl", "delta_ctrl_orth"), ("d_both", "delta_both_orth")):
                            x = v.get(kk) or {}
                            for f, ff in (("est", "estimate"), ("ci975_lo", "ci975_lo"), ("ci975_hi", "ci975_hi")):
                                _put(P, q + f"{k}.{f}", x.get(ff))
                        _put(P, q + "d_int", v.get("delta_interaction_orth"))
                        _put(P, q + "share", v.get("batch_share_orth"))
            hh = o.get("h_j1_holm")
            if hh:
                _put(P, p + "hj1.holm.verdict", hh.get("classification"))
                _put(P, p + "hj1.holm.subreason", hh.get("subreason") or "")
                _put(P, p + "hj1.holm.alpha_eff", hh.get("alpha_eff"))
                _put(P, p + "hj1.holm.ci_lo", hh.get("ci_adj_lo"))
    tmap = {"H-J1": "H-J1", "H-J3a": "H-J3a", "H-J3b": "H-J3b"}
    for fam, tab in [("holm", t) for m, t in (S.get("holm_by_model") or {}).items()] + \
                    [("pholm", S.get("holm_pooled_sensitivity") or {})]:
        for st in tab.get("steps") or []:
            m, ln, t = st["key"].split(":")
            ck = f"{fam}|{m}|{ln}|{tmap.get(t, t)}"
            for k, kk in (("rank", "rank"), ("raw_p", "p_raw"), ("step_alpha", "alpha_step"), ("effective_alpha", "alpha_eff"),
                          ("holm_adj_p", "p_holm"), ("holm_rejected", "holm_p_rejected"), ("rejected", "rejected")):
                _put(P, f"{ck}|{k}", st.get(kk))
    for m, c in (S.get("claim_rollups") or {}).items():
        _put(P, f"cohort|{m}|K_conf", c.get("K_conf"))
        _put(P, f"cohort|{m}|k", c.get("k_batch_dependence_practically_meaningful"))
        _put(P, f"cohort|{m}|claim", c.get("general_cross_context_claim_permitted"))
    j = S.get("j5b_rank_correlations") or {}
    smap = {"confirmatory_pairs": "conf", "pfizer_within_family": "pfizer", "fresh_n_ge_16_pairs": "fresh16"}
    for m, sets in (j.get("by_model") or {}).items():
        for sname, lst in sets.items():
            if sname not in smap or not isinstance(lst, list):
                continue
            for v in lst:
                a, b = sorted((v["line_1"], v["line_2"]))
                base = f"j5b|{m}|{smap[sname]}|{a}__vs__{b}|"
                for k, kk in (("rho", "spearman_rho"), ("p_perm", "p_perm"), ("bh", "bh_fdr"), ("top5", "top5_shared"),
                              ("bot5", "bot5_shared")):
                    _put(P, base + k, v.get(kk))
    for v in j.get("cross_model_within_line") or []:
        for k, kk in (("rho", "spearman_rho"), ("p_perm", "p_perm"), ("bh", "bh_fdr")):
            _put(P, f"j5b|cross|{v['line']}|{k}", v.get(kk))
    return P


# ----------------------------------------------------------------------------------------------- compare
def _num(v):
    return isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool)


def compare(C, P, tol):
    rows, worst, worst_key, n_checks, missing = [], 0.0, None, 0, []
    compare.worst_matched = (0.0, None)
    for k, cv in C.items():
        if k not in P:
            missing.append(k)
            continue
        pv = P[k]
        n_checks += 1
        if _num(cv) and _num(pv):
            cv, pv = float(cv), float(pv)
            if math.isnan(cv) and math.isnan(pv):
                continue
            diff = abs(cv - pv) if (math.isfinite(cv) and math.isfinite(pv)) else (0.0 if cv == pv else math.inf)
            if diff > worst:
                worst, worst_key = diff, k
            if not diff <= tol:
                rows.append((k, cv, pv, diff))
            elif diff > compare.worst_matched[0]:
                compare.worst_matched = (diff, k)
        else:
            if isinstance(cv, (list, tuple)):
                cv = [float(x) if _num(x) else x for x in cv]
                pv = [float(x) if _num(x) else x for x in pv] if isinstance(pv, (list, tuple)) else pv
            if cv != pv:
                rows.append((k, cv, pv, None))
    # set-membership families: a party pair or Holm test that the controller does not have is itself a difference
    ckeys = {k.rsplit("|", 1)[0] for k in C if k.startswith(("j5b|", "holm|", "pholm|"))}
    for k in sorted({k.rsplit("|", 1)[0] for k in P if k.startswith(("j5b|", "holm|", "pholm|"))} - ckeys):
        n_checks += 1
        rows.append((k + "|membership", "absent", "present", None))
    return rows, worst, worst_key, n_checks, missing


def diagnose(k, cv, pv, R):
    """Controller notes for known disagreement classes (protocol text cited)."""
    parts = k.split("|")
    C15 = ("PARTY: p_Delta computed two-sided; C1.5 defines the one-sided-doubled 'p_Delta = min(1, 2*(1 + "
           "#{Delta_CF^(u) <= 0}) / (B_valid + 1))' (the two-sided form is for H-J3a/b only)")
    if len(parts) >= 3 and parts[0] in R["lines"]:
        o = R["lines"][parts[0]]["models"].get(parts[1], {})
        if o.get("half_lt2") and parts[2].startswith("hj1."):
            return ("PARTY: C2.3 says 'Whatever any code computes, the reported H-J1 call for such a line is INCONCLUSIVE "
                    "(subreason VALIDITY_FAIL, detail HALF_SIZE_LT_2)'; this line's smaller half has < 2 perturbations")
        if o.get("half_lt2") and (parts[2] == "has_signal" or parts[2].endswith(".classification")) \
                and o["sig"]["has_signal"] != o["sig"]["has_signal_with_oracle"]:
            return ("PROTOCOL READING: C2.3 makes the cross-fitted oracle undefined on this line; the controller's J2 "
                    "predicate uses L1_ens56 only, the party's appears to use the oracle as computed")
        if o.get("half_lt2") and parts[2] == "j5a.cf_pair_reported":
            return "PARTY: C2.3 says '... and no J5a cross-fitted pair is reported' on a line whose smaller half has < 2"
        if parts[2] in ("cf.p_delta", "hj1.p_iut") and _num(pv) and "cf" in o:
            alt = o["cf"]["dcf"]["p"]
            alt = alt if parts[2] == "cf.p_delta" else max(o["vc"]["p_var"], alt)
            if abs(float(pv) - alt) <= 1e-12:
                return C15
    if k.startswith(("holm|", "pholm|")) and k.endswith("|raw_p") and parts[-2] == "H-J1" and _num(pv):
        o = R["lines"].get(parts[2], {}).get("models", {}).get(parts[1], {})
        if "cf" in o and abs(float(pv) - max(o["vc"]["p_var"], o["cf"]["dcf"]["p"])) <= 1e-12:
            return C15
    if k.endswith("|membership"):
        return ("party reports a J5b pair / Holm test outside the controller's set (PREREG_J §3.1 family; §3.5 sets "
                "[i]-[iv] with C1.10 / C2.2: a549_ko_cas9pos excluded from [iii], cytokine included)")
    if parts[0] == "j5b" and (k.endswith("|p_perm") or k.endswith("|bh")):
        J = R.get("j5b", {})
        v = (J.get("cross", {}).get(parts[2]) if parts[1] == "cross"
             else J.get("by_model", {}).get(parts[1], {}).get(parts[2], {}).get(tuple(parts[3].split("__vs__"))))
        n_perm = R["cfg"]["n_perm"]
        if v and v.get("n_tie", 0) and _num(cv) and _num(pv) and k.endswith("|p_perm") \
                and abs(float(cv) - float(pv)) <= v["n_tie"] / (n_perm + 1) + 1e-12:
            return (f"PARTY (floating point): {v['n_tie']} permutation(s) give |rho_k| exactly equal to |rho_obs| "
                    f"(exact integer Spearman); PREREG_J §3.5 / C2.2 count them ('>='), the party's float |rho_k| fell "
                    f"below |rho_obs|; difference = {round(abs(float(cv) - float(pv)) * (n_perm + 1))} count(s)")
        if v and v.get("n_tie", 0) and k.endswith("|bh"):
            return "follows from the p_perm exact-tie difference in the same pair set"
        return "check the permutation-table convention (PREREG_J §3.5, C1.10, C2.2)"
    if k.startswith(("holm|", "pholm|")) and (k.endswith("|rank") or k.endswith("|step_alpha")):
        fam = parts[0]
        raw_bad = [r[0] for r in diagnose.rows if r[0].startswith(fam + "|") and r[0].endswith("|raw_p")
                   and (fam == "pholm" or r[0].split("|")[1] == parts[1])]
        if raw_bad:
            return "follows from the raw-p difference(s) in this Holm family: " + ", ".join(x.split("|", 1)[1][:-6] for x in raw_bad)
        return "Holm tie order (PREREG_J §3.1 / C2.1 do not specify how tied p-values are ordered)"
    return ""


def party_holm_orders(P):
    fam = {}
    for k, v in P.items():
        parts = k.split("|")
        if parts[0] in ("holm", "pholm") and parts[-1] == "rank":
            key = "|".join(parts[1:-1])
            f = parts[1] if parts[0] == "holm" else "pooled"
            fam.setdefault(f, []).append((v, key))
    return {f: [k for _, k in sorted(lst)] for f, lst in fam.items()}


def _self_test():
    """Branch tests of the controller's adjudication code on hand-built inputs (no party data)."""
    rng = np.random.default_rng(1)

    def mk(dcf_mu, dcf_sd, p_var, tau, sig, valid=True, half_lt2=False):
        dr = rng.normal(dcf_mu, dcf_sd, 4000)
        return {"half_lt2": half_lt2, "hj1_valid": valid, "vc": {"p_var": p_var, "tau": tau}, "sig": {"has_signal": sig},
                "cf": {"dcf": {"draws": dr, "est": dcf_mu, "ci95": pct_ci(dr, 0.05)}}}

    def d3(mu, sd, tie=False, vf=1.0):
        dr = np.full(4000, 0.0) if tie else rng.normal(mu, sd, 4000)
        return {"draws": dr, "tie": tie, "valid_frac": vf}

    V = lambda o, a=ALPHA_LINE, h=None: hj1_verdict(o, a, h)["verdict"]
    T = [
        (V(mk(0.10, 0.01, 0.001, 0.05, True)), "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL"),
        (V(mk(0.03, 0.005, 0.001, 0.05, True)), "BATCH_DEPENDENCE_DETECTED_SUB_DELTA"),
        (V(mk(0.10, 0.01, 0.001, 0.015, True)), "BATCH_DEPENDENCE_DETECTED_SUB_DELTA"),
        (V(mk(0.00, 0.01, 0.5, 0.0, False)), "INCONCLUSIVE_LOW_SIGNAL"),
        (V(mk(0.00, 0.01, 0.5, 0.0, True)), "BATCH_INVARIANT_WITHIN_DELTA"),
        (V(mk(0.00, 0.03, 0.5, 0.0, True)), "INCONCLUSIVE"),
        (hj1_verdict(mk(0.00, 0.03, 0.5, 0.0, True), ALPHA_LINE, None)["subreason"], "INSUFFICIENT_PRECISION"),
        (V(mk(0.00, 0.01, 0.5, 0.03, True)), "INCONCLUSIVE"),
        (hj1_verdict(mk(0.10, 0.01, 0.001, 0.05, True, valid=False), ALPHA_LINE, None)["subreason"], "VALIDITY_FAIL"),
        (hj1_verdict(mk(0.10, 0.01, 0.001, 0.05, True, half_lt2=True), ALPHA_LINE, None)["detail"], "HALF_SIZE_LT_2"),
        (V(mk(0.10, 0.01, 0.001, 0.05, True), 0.005, False), "INCONCLUSIVE"),           # Holm not rejected
        (V(mk(0.10, 0.01, 0.001, 0.05, True), 0.005, True), "BATCH_DEPENDENCE_PRACTICALLY_MEANINGFUL"),
        (V(mk(0.01, 0.005, 0.001, 0.05, True), 0.005, True), "INCONCLUSIVE"),          # Holm rejects, gate fails
        (hj3_verdict(d3(0, 0, tie=True), ALPHA_LINE, None)["verdict"], "EXACT_TIE_SELECTED_ZERO"),
        (hj3_verdict(d3(0.05, 0.01), ALPHA_LINE, None)["verdict"], "RULE_SUPERIOR_TO_ZERO"),
        (hj3_verdict(d3(-0.05, 0.01), ALPHA_LINE, None)["verdict"], "ZERO_SUPERIOR_TO_RULE"),
        (hj3_verdict(d3(0.0, 0.01), ALPHA_LINE, None)["verdict"], "INDISTINGUISHABLE"),
        (hj3_verdict(d3(0.05, 0.01), ALPHA_LINE, False)["verdict"], "INDISTINGUISHABLE"),
        (hj3_verdict(d3(0.05, 0.01, vf=0.9), ALPHA_LINE, None)["verdict"], "INDISTINGUISHABLE"),
        (p_two_sided([1, 2, 3, -1]), min(1, 2 * 2 / 5)), (p_pos([1, 2, 3, -1]), min(1, 2 * 2 / 5)),
        (p_pos([-1, -2, -3, 1]), 1.0), (p_two_sided([-1, -2, -3, 1]), min(1, 2 * 2 / 5)),
    ]
    I, M, F = "INDISTINGUISHABLE", "MODEL_SUPERIOR_TO_FCT", "FCT_SUPERIOR_TO_MODEL"
    T += [(hj2_class([F] + [M] + [F] * 54, True)["classification"], "CALL_UNSTABLE_SIGN_REVERSAL"),
          (hj2_class([F] + [I] * 5 + [F] * 50, True)["classification"], "CALL_UNSTABLE_SUBSTANTIAL"),
          (hj2_class([F] + [I] * 4 + [F] * 51, True)["classification"], "CALL_UNSTABLE_MINOR"),
          (hj2_class([F] * 56, True)["classification"], "CALL_STABLE"),
          (hj2_class([I] * 56, False)["classification"], "INCONCLUSIVE_LOW_SIGNAL"),
          (hj2_class([I] * 55 + [M], False)["classification"], "CALL_UNSTABLE_MINOR")]
    g_ok, g_bad = (lambda a: True), (lambda a: False)
    h = holm([("a", 0.001, g_ok), ("b", 0.01, g_bad), ("c", 0.02, g_ok), ("d", 0.04, g_ok)])
    T += [([h[k]["rejected"] for k in "abcd"], [True, False, True, True]),
          ([h[k]["effective_alpha"] for k in "abcd"], [0.0125, 0.05 / 3, 0.025, 0.05])]
    h = holm([("a", 0.001, g_ok), ("b", 0.03, g_ok), ("c", 0.04, g_ok)])
    T += [([h[k]["holm_rejected"] for k in "abc"], [True, False, False]),
          ([h[k]["effective_alpha"] for k in "abc"], [0.05 / 3, 0.025, 0.025])]
    h = holm([("a", 0.004, g_ok), ("b", 0.004, g_ok)], order_keys=["b", "a"])
    T += [(h["b"]["rank"], 1), (h["a"]["step_alpha"], 0.05)]
    tm = tie_material([("a", 0.004, lambda a: a > 0.03), ("b", 0.004, g_ok)], holm([("a", 0.004, g_ok), ("b", 0.004, g_ok)]))
    T += [([t["key"] for t in tm], ["a"])]
    x = np.arange(56.0)
    perm = np.stack([np.random.default_rng(3).permutation(56) for _ in range(50)])
    rho, pp, ties, nt = spearman_perm(x, x, perm)
    T += [(rho, 1.0), (pp, 1 / 51), (spearman_perm(x, -x, perm)[0], -1.0)]
    P1, P2 = split_c14(np.array(["b", "a", "c", "d", "e"]), np.array([True, True, False, True, True]), 0, 7)
    T += [((P1.size, P2.size), (2, 3)), (sorted(P1.tolist() + P2.tolist()), [0, 1, 2, 3, 4])]
    T += [(round(Z_BONF56, 4), 3.5110)]
    bad = [(i, got, exp) for i, (got, exp) in enumerate(T) if got != exp]
    print(f"r17_check self-test: {len(T) - len(bad)}/{len(T)} branch tests pass", "" if not bad else bad)
    return not bad


def detect_party(S):
    if "lines" in S and "separate_model_holm_fwer05" in S:
        return "main"
    if "models" in S and "holm_by_model" in S:
        return "at3"
    raise SystemExit("cannot detect the party schema of summary_r17.json")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("results_dir", nargs="?")
    ap.add_argument("--self-test", action="store_true", help="run the adjudication branch tests and exit")
    ap.add_argument("--party", default="auto", choices=["auto", "main", "at3"])
    ap.add_argument("--tol", type=float, default=1e-5)
    ap.add_argument("--n-boot", type=int)
    ap.add_argument("--n-boot-bonf56", type=int)
    ap.add_argument("--n-perm", type=int)
    ap.add_argument("--n-splits-diag", type=int)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--skip-j5b", action="store_true")
    ap.add_argument("--skip-j5a", action="store_true")
    ap.add_argument("--out", help="write the controller recompute (flat) + diff list as JSON")
    ap.add_argument("--max-print", type=int, default=200)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    if a.self_test:
        return 0 if _self_test() else 1
    if not a.results_dir:
        ap.error("results_dir is required (or pass --self-test)")
    rd = Path(a.results_dir)
    S = json.loads((rd / "summary_r17.json").read_text())
    party = detect_party(S) if a.party == "auto" else a.party
    cfg = {"n_boot": a.n_boot or int(S.get("n_boot", 4000)),
           "n_boot_bonf56": a.n_boot_bonf56 or int(S.get("n_boot_bonf56", 20000)),
           "n_perm": a.n_perm or int(S.get("n_perm", 10000)),
           "n_splits_diag": a.n_splits_diag or int(S.get("n_splits_diag", 20)),
           "seed": a.seed or int(S.get("seed", 20261701)), "skip_j5b": a.skip_j5b, "skip_j5a": a.skip_j5a}
    log = (lambda *x: None) if a.quiet else (lambda *x: print(*x, file=sys.stderr))
    log(f"r17_check: party={party} cfg={cfg}")
    R = recompute(rd, cfg, log)
    C = flat_controller(R)
    P = flat_main(S) if party == "main" else flat_at3(S)
    notes = []
    orders = party_holm_orders(P)
    mine = {m: [k for k, _ in sorted(R["holm"][m].items(), key=lambda kv: kv[1]["rank"])] for m in R["holm"]}
    if "holm_pooled" in R:
        mine["pooled"] = [k for k, _ in sorted(R["holm_pooled"].items(), key=lambda kv: kv[1]["rank"])]
    diff_fams = [f for f in orders if f in mine and orders[f] != mine[f]]
    if diff_fams:
        apply_holm(R, orders)
        new = {m: [k for k, _ in sorted(R["holm"][m].items(), key=lambda kv: kv[1]["rank"])] for m in R["holm"]}
        if "holm_pooled" in R:
            new["pooled"] = [k for k, _ in sorted(R["holm_pooled"].items(), key=lambda kv: kv[1]["rank"])]
        for f in diff_fams:
            if new.get(f) == orders[f]:
                notes.append(f"Holm family {f}: party breaks p-value ties in a different order; it is a valid p-sorted "
                             f"order, so the controller re-ran Holm in the party's tie order (PREREG_J / C2.1 do not "
                             f"specify tie order)")
            else:
                notes.append(f"Holm family {f}: party rank order is not p-sorted under the controller's p-values; kept "
                             f"the controller order")
        C = flat_controller(R)
    for f, lst in R["holm_tie_material"].items():
        for t in lst:
            notes.append(f"WARNING Holm family {f}: {t['key']} (raw p {t['raw_p']:.6g}) passes its C2.1 gate at one end of "
                         f"its tie group's alpha range {t['alpha_range']} but not the other: the call depends on the "
                         f"unspecified tie order")
    rows, worst, worst_key, n_checks, missing = compare(C, P, a.tol)
    # array-level checks are controller-only invariants (must be exactly 0 for flags, <= 1e-12 for arrays)
    arr_fail = [(f"arrays|{lm}|{k}", v) for lm, ch in R["chk"].items() for k, v in ch.items() if v > 1e-12]
    rep = {"party": party, "n_checks": n_checks + sum(len(c) for c in R["chk"].values()),
           "n_fail": len(rows) + len(arr_fail), "n_not_exported": len(missing), "worst_abs_diff": worst,
           "worst_key": worst_key, "worst_abs_diff_within_tol": compare.worst_matched[0],
           "worst_key_within_tol": compare.worst_matched[1], "tol": a.tol, "cfg": {k: v for k, v in cfg.items()}, "confirmatory": R["conf"]}
    locked = {"n_boot": 4000, "n_boot_bonf56": 20000, "n_perm": 10000, "n_splits_diag": 20, "seed": 20261701}
    off = {k: cfg[k] for k in locked if cfg[k] != locked[k]}
    if off:
        notes.append(f"configuration differs from the C1.13 locked run configuration {locked}: {off} "
                     f"(acceptable for synthetic smoke runs only)")
    rep["notes"] = notes
    print(json.dumps(rep))
    for n in notes:
        print("NOTE " + n)
    diagnose.rows = rows
    for k, cv, pv, diff in rows[:a.max_print]:
        dg = diagnose(k, cv, pv, R)
        print(f"DIFF {k}: controller={_short(cv)} party={_short(pv)}" + (f" |diff|={diff:.3g}" if diff is not None else "")
              + (f"  [{dg}]" if dg else ""))
    if len(rows) > a.max_print:
        print(f"... {len(rows) - a.max_print} more DIFF rows")
    for k, v in arr_fail:
        print(f"DIFF {k}: {v}")
    if a.out:
        Path(a.out).write_text(json.dumps({"report": rep, "controller": {k: _js(v) for k, v in C.items()},
                                           "diffs": [[k, _js(cv), _js(pv), diff] for k, cv, pv, diff in rows],
                                           "not_exported": missing, "array_checks": R["chk"]}, indent=1))
    print("R17_CHECK_MATCH" if rep["n_fail"] == 0 else "R17_CHECK_DIFF")
    return 0 if rep["n_fail"] == 0 else 1


def _js(v):
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, (list, tuple)):
        return [_js(x) for x in v]
    if isinstance(v, float) and not math.isfinite(v):
        return str(v)
    return v


def _short(v):
    s = json.dumps(_js(v)) if not isinstance(v, str) else v
    return s if len(s) <= 160 else s[:157] + "..."


if __name__ == "__main__":
    raise SystemExit(main())
