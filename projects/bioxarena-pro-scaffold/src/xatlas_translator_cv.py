#!/usr/bin/env python3
"""
Standalone 5-Fold Cross-Validation Verification of the Cross-Line Perturbation Translator
(HEK293T <-> HCT116) on paired reproducible training knockdowns.

Note on D7 (Comparability & Inference-Time Availability):
- This script measures cross-line translation accuracy on paired training targets (N=781)
  solely to select the architecture of the auxiliary translator T_{A->B}.
- Because 0/455 test targets appear in train in either cell line, T_{A->B} is NEVER used to
  copy test targets from the other cell line. Instead, T_{A->B} is used strictly at train/bridge
  time to (1) translate other-line training labels for unseen targets and (2) bridge predicted
  responses across the 4,040 shared genes (385/455 = 84.6% test coverage).
"""

import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

candidate_dirs = [
    *( [Path(os.environ["BIOXARENA_DIR"])] if os.environ.get("BIOXARENA_DIR") else [] ),
    Path.cwd() / "BioXArena-Pro",
    Path(__file__).resolve().parents[4] / "BioXArena-Pro",
]
BENCHMARK_DIR = next((c for c in candidate_dirs if (c / "data/perturbation/xatlas-orion-perturb/public").exists()), candidate_dirs[0])
DATA_DIR = BENCHMARK_DIR / "data/perturbation/xatlas-orion-perturb/public"


def proj5000(Y, base, ax):
    Q, Rm = np.linalg.qr(np.stack([np.ones(5000), base, ax], axis=1))
    Q = Q[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())]
    Ys = Y.astype(np.float64) - base[None, :]
    return Ys - (Ys @ Q) @ Q.T


def eval_pair_cos(Y_true, Y_pred, base, ax, s_idx):
    Q, Rm = np.linalg.qr(np.stack([np.ones(len(s_idx)), base[s_idx], ax[s_idx]], axis=1))
    Q = Q[:, np.abs(np.diag(Rm)) > 1e-10 * max(1.0, np.abs(Rm).max())]
    scores = []
    for i in range(len(Y_true)):
        a = Y_true[i, s_idx]; a = a - Q @ (Q.T @ a)
        b = Y_pred[i, s_idx]; b = b - Q @ (Q.T @ b)
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        scores.append(float(np.dot(a, b) / (na * nb)) if na > 1e-8 and nb > 1e-8 else 0.0)
    return float(np.mean(scores)), float(np.std(scores))


def signed_power(M, p):
    return M if p == 1.0 else np.sign(M) * (np.abs(M) ** p)


def main():
    train_df = pd.read_csv(DATA_DIR / "train.csv")
    genes_npz = np.load(DATA_DIR / "genes.npz")
    baselines_npz = np.load(DATA_DIR / "context_baseline.npz")
    axes_npz = np.load(DATA_DIR / "residual_axis.npz")
    score_idx_npz = np.load(DATA_DIR / "score_gene_indices.npz")
    train_y_npz = np.load(DATA_DIR / "train_y.npz")

    hct_g = list(genes_npz["HCT116_genes"])
    hek_g = list(genes_npz["HEK293T_genes"])
    sh = sorted(list(set(hct_g) & set(hek_g)))
    h_sh = np.array([hct_g.index(g) for g in sh])
    k_sh = np.array([hek_g.index(g) for g in sh])

    baselines = {k: np.asarray(baselines_npz[k], dtype=np.float64) for k in baselines_npz.files}
    axes = {k: np.asarray(axes_npz[k], dtype=np.float64) for k in axes_npz.files}
    score_indices = {k: np.asarray(score_idx_npz[k], dtype=np.int64) for k in score_idx_npz.files}

    h_rep = train_df[(train_df["source"] == "HCT116") & (train_df["reproducible"] == 1)]
    k_rep = train_df[(train_df["source"] == "HEK293T") & (train_df["reproducible"] == 1)]
    hp = dict(zip(h_rep["perturbation"], h_rep["id"]))
    kp = dict(zip(k_rep["perturbation"], k_rep["id"]))
    paired = sorted(list(set(hp.keys()) & set(kp.keys())))

    Yh = proj5000(np.array([train_y_npz[hp[p]] for p in paired], dtype=np.float32), baselines["HCT116"], axes["HCT116"])
    Yk = proj5000(np.array([train_y_npz[kp[p]] for p in paired], dtype=np.float32), baselines["HEK293T"], axes["HEK293T"])

    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    summary = {"n_paired_reproducible": len(paired), "n_shared_genes": len(sh), "directions": {}}

    for src, dst, Xs, Yd, s_sh, d_sh in [
        ("HEK293T", "HCT116", Yk, Yh, k_sh, h_sh),
        ("HCT116", "HEK293T", Yh, Yk, h_sh, k_sh),
    ]:
        s_idx = score_indices[f"{dst}_indices"]
        P0 = np.zeros_like(Yd); P0[:, d_sh] = Xs[:, s_sh]
        m0, sd0 = eval_pair_cos(Yd, P0, baselines[dst], axes[dst], s_idx)

        oof_gs = np.zeros_like(Yd)
        for tr, va in kf.split(Xs):
            xs_tr = Xs[tr][:, s_sh]; yd_tr = Yd[tr][:, d_sh]
            slope = np.clip(np.sum(xs_tr * yd_tr, axis=0) / (np.sum(xs_tr**2, axis=0) + 0.5), 0.0, 3.0)
            Ysk_tr = np.zeros_like(Yd[tr]); Ysk_tr[:, d_sh] = 0.8 * slope[None, :] * Xs[tr][:, s_sh]
            Ysk_va = np.zeros_like(Yd[va]); Ysk_va[:, d_sh] = 0.8 * slope[None, :] * Xs[va][:, s_sh]
            norm_tr = np.linalg.norm(Xs[tr], axis=1, keepdims=True) + 1e-8
            norm_va = np.linalg.norm(Xs[va], axis=1, keepdims=True) + 1e-8
            Xn_tr = Xs[tr] / norm_tr; Xn_va = Xs[va] / norm_va
            K_tr = signed_power(Xn_tr @ Xn_tr.T, 1.5); K_tr[np.diag_indices_from(K_tr)] += 2.0
            K_va = signed_power(Xn_va @ Xn_tr.T, 1.5)
            dual = np.linalg.solve(K_tr, (Yd[tr] - Ysk_tr) / norm_tr)
            oof_gs[va] = Ysk_va + (K_va @ dual) * norm_va
        m_gs, sd_gs = eval_pair_cos(Yd, oof_gs, baselines[dst], axes[dst], s_idx)

        summary["directions"][f"{src}->{dst}"] = {
            "raw_shared_copy_mean": round(m0, 4),
            "raw_shared_copy_sd": round(sd0, 4),
            "geneslope_coskrr_p15_mean": round(m_gs, 4),
            "geneslope_coskrr_p15_sd": round(sd_gs, 4),
        }
        print(f"[{src} -> {dst}] Raw shared copy: {m0:.4f} +/- {sd0:.4f} | GeneSlope+CosKRR(p=1.5): {m_gs:.4f} +/- {sd_gs:.4f}")

    out_path = Path(__file__).resolve().parents[1] / "reports" / "xatlas_translator_5fold_cv.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved translator CV audit to {out_path}")


if __name__ == "__main__":
    main()
