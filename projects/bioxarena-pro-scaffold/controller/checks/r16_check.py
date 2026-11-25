"""Controller (third-party) recomputation for Round 16 from exported arrays.

Metric code is carried over unchanged from controller/checks/r15_check.py (frozen in R15, predates R16):
uncentered orthogonalization of residuals against the per-perturbation Transfer_Mean, pooled moment
correlation, Spearman-Brown ceiling from orthogonalized half A / half B residuals on half-valid targets,
and joint perturbation bootstrap weights from np.random.default_rng(seed).integers(0, N, N) per draw.
Adjudication comes from controller/checks/r16_check_adjudicate.py (PREREG_I §3 text only).
Never imports round16/src or round16/audit.

Usage: python3 -I r16_check.py <line_arrays_dir> <seen_csv> <unseen_csv> <out.json> [n_boot] [seed]
  <line_arrays_dir> holds observed.npz and preds_STATE_SE_fixed_index_zero_{matched,mis_<src>}.npz
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r16_check_adjudicate as ADJ  # noqa: E402

LO975, HI975 = 1.25, 98.75
LO95, HI95 = 2.5, 97.5


def orth_row(x, u):
    return x - np.dot(x, u) * u


class Line:
    def __init__(self, obs):
        self.D = obs["delta_obs"]; self.A = obs["delta_half_a"]; self.B = obs["delta_half_b"]
        self.TM = obs["delta_transfer"]; self.HV = obs["half_valid"].astype(bool)
        self.selfg = obs["target_indices"].astype(int)
        self.N, self.G = self.D.shape
        self.M, self.U, self.O = [], [], []
        self.cs = np.zeros((self.N, 6))
        for r in range(self.N):
            m = np.ones(self.G, bool)
            if 0 <= self.selfg[r] < self.G:
                m[self.selfg[r]] = False
            m &= np.isfinite(self.D[r]) & np.isfinite(self.TM[r])
            m = np.flatnonzero(m); self.M.append(m)
            tm = self.TM[r, m]; u = tm / max(np.linalg.norm(tm), 1e-12); self.U.append(u)
            self.O.append(orth_row(self.D[r, m] - tm, u))
            if self.HV[r]:
                a = self.A[r, m]; b = self.B[r, m]; ok = np.isfinite(a) & np.isfinite(b)
                if ok.all():
                    a = orth_row(a - tm, u); b = orth_row(b - tm, u)
                else:  # pairwise-complete fallback, as in R15
                    a = orth_row(np.where(ok, a, tm) - tm, u)[ok]; b = orth_row(np.where(ok, b, tm) - tm, u)[ok]
                self.cs[r] = [a.size, a.sum(), b.sum(), (a * b).sum(), (a * a).sum(), (b * b).sum()]

    def stats(self, P):
        st = np.zeros((self.N, 6))
        for r, m in enumerate(self.M):
            y = orth_row(P[r, m] - self.TM[r, m], self.U[r]); x = self.O[r]
            st[r] = [x.size, x.sum(), y.sum(), (x * y).sum(), (x * x).sum(), (y * y).sum()]
        return st


def pooled_r(st, w):
    n, sx, sy, sxy, sxx, syy = (w @ st)
    cov = sxy - sx * sy / n; vx = sxx - sx * sx / n; vy = syy - sy * sy / n
    return cov / np.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan


def cal(L, st, w):
    rab = pooled_r(L.cs, w * L.HV)
    if not rab > 0:
        return np.nan, rab
    return float(np.clip(pooled_r(st, w) / np.sqrt(2 * rab / (1 + rab)), -1, 1)), rab


def summarize(pt, bs):
    return {"point": float(pt), "valid_frac": float(np.mean(np.isfinite(bs))),
            "ci975": tuple(float(v) for v in np.nanpercentile(bs, [LO975, HI975])),
            "ci95": tuple(float(v) for v in np.nanpercentile(bs, [LO95, HI95]))}


def run(line_dir, seen, unseen, n_boot=4000, seed=20261601):
    d = Path(line_dir)
    obs = dict(np.load(d / "observed.npz", allow_pickle=False))
    L = Line(obs)
    arms = {"matched": np.load(d / "preds_STATE_SE_fixed_index_zero_matched.npz")["delta_pred"]}
    for s in list(seen) + list(unseen):
        arms["mis_" + s] = np.load(d / f"preds_STATE_SE_fixed_index_zero_mis_{s}.npz")["delta_pred"]
    arms["FCT"] = obs["fct_delta_pred"]
    rng = np.random.default_rng(seed)
    W = [np.bincount(rng.integers(0, L.N, L.N), minlength=L.N).astype(float) for _ in range(n_boot)]
    ones = np.ones(L.N)
    ST = {k: L.stats(v) for k, v in arms.items()}
    PT = {k: cal(L, st, ones)[0] for k, st in ST.items()}
    rab = cal(L, ST["matched"], ones)[1]
    BS = {k: np.array([cal(L, st, w)[0] for w in W]) for k, st in ST.items()}
    out = {"N": L.N, "G": L.G, "n_boot": n_boot, "seed": seed, "rho_ab_orth": float(rab),
           "arms": {k: summarize(PT[k], BS[k]) for k in arms}}
    pairs = {s: summarize(PT["matched"] - PT["mis_" + s], BS["matched"] - BS["mis_" + s]) for s in list(seen) + list(unseen)}
    out["pairs"] = {s: dict(v, call=ADJ.classify_pair(v), within_delta=ADJ.within_delta(v)) for s, v in pairs.items()}
    m = out["arms"]["matched"]; matched = {"valid_frac": m["valid_frac"], "rho_ab_orth": out["rho_ab_orth"], "ci975": m["ci975"]}
    out["H_I1_line"] = ADJ.adjudicate_line_h_i1(matched, {s: pairs[s] for s in seen}, {s: pairs[s] for s in unseen})
    h2 = summarize(PT["matched"] - PT["FCT"], BS["matched"] - BS["FCT"])
    out["H_I2"] = dict(h2, call=ADJ.classify_h_i2(h2))
    return out


if __name__ == "__main__":
    line_dir, seen, unseen, outp = sys.argv[1], sys.argv[2].split(","), sys.argv[3].split(","), sys.argv[4]
    nb = int(sys.argv[5]) if len(sys.argv) > 5 else 4000
    sd = int(sys.argv[6]) if len(sys.argv) > 6 else 20261601
    res = run(line_dir, [s for s in seen if s], [s for s in unseen if s], nb, sd)
    Path(outp).write_text(json.dumps(res, indent=1))
    print(json.dumps({"N": res["N"], "G": res["G"], "rho_ab_orth": res["rho_ab_orth"], "H_I1_line": res["H_I1_line"],
                      "H_I2": res["H_I2"]["call"], "matched": res["arms"]["matched"]["point"]}, indent=1))
