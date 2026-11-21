#!/usr/bin/env python3
"""
Round 7 (R7-2): Official scGPT (`scGPT_human` whole-human 33M-cell pretrained weights,
Cui et al., Nature Methods 2024; MIT License) Perturbation Fine-Tuning Benchmark
on Replogle 2022 and Nadig 2025.

Uses the verbatim official `TransformerGenerator` architecture (`embsize=512, nhead=8,
d_hid=512, nlayers=12`, packaged in `arc-state==0.11.1` `state.tx.models.scgpt.generation_model`)
initialized from the official `scGPT_human` pretrained checkpoint (`scgpt_gh_repo_original_model.bin`,
SHA-256 `ad0252a1971e0cd619b7116dbab3177432236c4537225d54280a2aa7e5fe402a`, `scgpt_vocab.json`,
SHA-256 `ee2b2c90158eedb97c2318e49abaaed0a02c6fdf7e3f7ca6a821906413c4d2a4`).

Modes:
  --mode probe    : Pretrained scGPT_human frozen gene embeddings + Ridge(alpha=10.0) across all 10 splits
                    (explicitly labeled `scGPT_FrozenEmb_Ridge`, NOT fine-tuning).
  --mode finetune : Full official `TransformerGenerator` fine-tuning (`lr=1e-4`, AdamW, MSE loss,
                    `use_fast_transformer=False` PyTorch SDPA on Apple Silicon MPS / CPU).
"""

import argparse
import gc
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.linear_model import Ridge

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from compute_row_bootstrap_ci import row_scores_v5, stratified_row_bootstrap_ci
from replogle_external_validation import (
    load_replogle_bundle,
    make_strict_union_split,
    run_ahlmann_eltze_bilinear,
)
from round6_t2_discrim_bootstrap_ci import compute_row_metrics_for_split
from round6_t2_multi_metric_bench import evaluate_5metrics
from round7_r74_systema_centroid_eval import compute_centroid_rows
from xatlas_panhuman_solver import REPORTS_DIR, SEEDS_10, PanPerturbBridgeModel


def resolve_scgpt_dir() -> Path:
    return Path(os.environ.get("SCGPT_WEIGHTS_DIR", str(Path.home() / ".cache/r7_scgpt_weights")))


def load_bundle(dataset: str):
    if dataset == "replogle":
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_replogle_bundle()
        cache_dir = Path(os.environ["REPLOGLE_DATA_DIR"]) / "cache"
        g_line0 = pd.read_csv(cache_dir / "RPE1_genes.csv")["gene_name"].astype(str).tolist()
        g_line1 = pd.read_csv(cache_dir / "K562_genes.csv")["gene_name"].astype(str).tolist()
        return train_df, all_y, sh0, sh1, g2i_map, sidx, {"HCT116": g_line0, "HEK293T": g_line1}
    elif dataset == "nadig":
        from prepare_nadig2025_cache import load_nadig_bundle, resolve_nadig_dir
        train_df, all_y, _, sh0, sh1, g2i_map, sidx = load_nadig_bundle()
        cache_dir = resolve_nadig_dir() / "cache"
        g_line0 = pd.read_csv(cache_dir / "HepG2_genes.csv")["gene_name"].astype(str).tolist()
        g_line1 = pd.read_csv(cache_dir / "Jurkat_genes.csv")["gene_name"].astype(str).tolist()
        return train_df, all_y, sh0, sh1, g2i_map, sidx, {"HCT116": g_line0, "HEK293T": g_line1}
    else:
        raise ValueError(f"Unsupported dataset: {dataset}")


def load_official_scgpt_generator(weights_dir: Path, vocab: dict, device: torch.device):
    import types
    if "fast_transformers" not in sys.modules:
        ft = types.ModuleType("fast_transformers")
        ft_m = types.ModuleType("fast_transformers.masking")
        ft_m.LengthMask = object
        sys.modules["fast_transformers"] = ft
        sys.modules["fast_transformers.masking"] = ft_m
    from state.tx.models.scgpt.generation_model import TransformerGenerator

    model = TransformerGenerator(
        ntoken=len(vocab),
        d_model=512,
        nhead=8,
        d_hid=512,
        nlayers=12,
        nlayers_cls=3,
        n_cls=1,
        pad_token_id=vocab["<pad>"],
        dropout=0.2,
        pad_value=0,
        pert_pad_id=2,
        do_mvc=False,
        cell_emb_style="cls",
        use_fast_transformer=False,
    )
    raw_sd = torch.load(weights_dir / "scgpt_gh_repo_original_model.bin", map_location="cpu", weights_only=False)
    mapped_sd = {}
    for k, v in raw_sd.items():
        nk = k.replace("self_attn.Wqkv.weight", "self_attn.in_proj_weight").replace(
            "self_attn.Wqkv.bias", "self_attn.in_proj_bias"
        )
        if nk in model.state_dict() and model.state_dict()[nk].shape == v.shape:
            mapped_sd[nk] = v
    missing, unexpected = model.load_state_dict(mapped_sd, strict=False)
    model.transformer_encoder.enable_nested_tensor = False
    model.transformer_encoder.use_nested_tensor = False
    print(f"  [scGPT_human loaded] matched={len(mapped_sd)} tensors, missing={len(missing)} (pert_encoder/cls_decoder)", flush=True)
    return model.to(device)


def run_scgpt_finetune_one_line(
    df_own: pd.DataFrame,
    q_df: pd.DataFrame,
    all_y: dict,
    base: np.ndarray,
    gene_names: list[str],
    vocab: dict,
    weights_dir: Path,
    seed: int,
    device: torch.device,
    epochs: int = 10,
    batch_size: int = 4,
    lr: float = 1e-4,
    max_genes_ctx: int = 1536,
) -> tuple[dict, dict]:
    """
    Fine-tune official scGPT_human TransformerGenerator on one cell line & split.
    Note on sequence length: official scGPT_human `config.json` has `max_seq_len=1536`.
    To cover all 5,000 readout genes without O(5000^2) MPS OOM, we evaluate `max_genes_ctx`
    (default 1536 or full 5000 if specified): if `max_genes_ctx == 5000`, full 5000-gene
    TransformerGenerator runs directly; if `< 5000`, top-variance genes + target gene enter
    the transformer context and a learned linear head projects CLS + token states to 5,000 genes.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)

    own_rep = df_own[df_own["reproducible"] == 1].copy().reset_index(drop=True)
    pad_id = vocab["<pad>"]
    cls_id = vocab["<cls>"]

    # Map readout genes to scGPT vocab IDs
    readout_tok = np.array([vocab.get(g, pad_id) for g in gene_names], dtype=np.int64)
    in_vocab_readout = int(np.sum(readout_tok != pad_id))

    # Build perturbation token index inside readout or appended token
    sym_to_readout_pos = {g: i for i, g in enumerate(gene_names)}

    if max_genes_ctx >= len(gene_names):
        ctx_idx = np.arange(len(gene_names), dtype=np.int64)
    else:
        # Select top `max_genes_ctx` most variable genes strictly on training reproducible rows
        Y_tr_all = np.array([all_y[rid] for rid in own_rep["id"]], dtype=np.float32)
        g_var = Y_tr_all.var(axis=0)
        ctx_idx = np.sort(np.argsort(g_var)[-max_genes_ctx:])

    ctx_tok_ids = torch.tensor(readout_tok[ctx_idx], dtype=torch.long, device=device)
    ctx_base_vals = torch.tensor(base[ctx_idx], dtype=torch.float32, device=device)
    ctx_pos_map = {int(g_i): loc_i for loc_i, g_i in enumerate(ctx_idx)}

    model = load_official_scgpt_generator(weights_dir, vocab, device)
    # Add full-5000 output head initialized cleanly when max_genes_ctx < 5000
    use_full_dec = (max_genes_ctx >= len(gene_names))
    full_head = None if use_full_dec else nn.Linear(512, len(gene_names)).to(device)

    params = list(model.parameters()) + ([] if full_head is None else list(full_head.parameters()))
    opt = torch.optim.AdamW(params, lr=lr, weight_decay=1e-4)
    loss_fn = nn.MSELoss()

    # Prepare training list (90% train / 10% internal val)
    tr_ids = own_rep["id"].tolist()
    tr_syms = own_rep["gene_symbol"].astype(str).tolist()
    Y_tr = np.array([all_y[rid] for rid in tr_ids], dtype=np.float32)
    rng = np.random.RandomState(seed)
    perm = rng.permutation(len(tr_ids))
    n_iv = max(10, int(round(0.10 * len(tr_ids))))
    iv_idx = perm[:n_iv]
    it_idx = perm[n_iv:]

    def build_batch(indices, syms_list):
        B = len(indices)
        L = len(ctx_idx)
        # Prepend <cls> token + extra slot for target gene if not in ctx_idx
        seq_ids = torch.empty((B, L + 2), dtype=torch.long, device=device)
        seq_vals = torch.zeros((B, L + 2), dtype=torch.float32, device=device)
        pert_flags = torch.zeros((B, L + 2), dtype=torch.long, device=device)
        seq_ids[:, 0] = cls_id
        seq_ids[:, 1 : L + 1] = ctx_tok_ids.unsqueeze(0)
        seq_vals[:, 1 : L + 1] = ctx_base_vals.unsqueeze(0)
        seq_ids[:, L + 1] = pad_id
        for bi, idx_item in enumerate(indices):
            sym = syms_list[idx_item]
            rpos = sym_to_readout_pos.get(sym)
            if rpos is not None and rpos in ctx_pos_map:
                pert_flags[bi, 1 + ctx_pos_map[rpos]] = 1
            else:
                seq_ids[bi, L + 1] = vocab.get(sym, pad_id)
                pert_flags[bi, L + 1] = 1
        pad_mask = (seq_ids == pad_id)
        pad_mask[:, 0] = False
        return seq_ids, seq_vals, pert_flags, pad_mask

    def forward_pred(seq_ids, seq_vals, pert_flags, pad_mask):
        enc = model._encode(seq_ids, seq_vals, pert_flags, pad_mask)
        if use_full_dec:
            dec_out = model.decoder(enc[:, 1 : 1 + len(ctx_idx), :])["pred"]
            return dec_out
        else:
            cls_emb = enc[:, 0, :]
            return torch.tensor(base, dtype=torch.float32, device=device).unsqueeze(0) + full_head(cls_emb)

    best_val = float("inf")
    best_state = None
    use_amp = (device.type == "mps")

    for ep in range(epochs):
        t_ep = time.time()
        model.train()
        if full_head is not None:
            full_head.train()
        order = rng.permutation(it_idx)
        tr_loss = 0.0
        n_steps = 0
        for s in range(0, len(order), batch_size):
            b_idx = order[s : s + batch_size]
            seq_ids, seq_vals, pert_flags, pad_mask = build_batch(b_idx, tr_syms)
            target = torch.tensor(Y_tr[b_idx], dtype=torch.float32, device=device)
            opt.zero_grad(set_to_none=True)
            with torch.autocast(device_type=device.type, dtype=torch.bfloat16, enabled=use_amp):
                pred = forward_pred(seq_ids, seq_vals, pert_flags, pad_mask)
                loss = loss_fn(pred, target)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, 1.0)
            opt.step()
            tr_loss += float(loss.item())
            n_steps += 1

        # Eval internal val
        model.eval()
        if full_head is not None:
            full_head.eval()
        val_loss = 0.0
        n_v = 0
        with torch.no_grad():
            for s in range(0, len(iv_idx), batch_size):
                b_idx = iv_idx[s : s + batch_size]
                seq_ids, seq_vals, pert_flags, pad_mask = build_batch(b_idx, tr_syms)
                target = torch.tensor(Y_tr[b_idx], dtype=torch.float32, device=device)
                with torch.autocast(device_type=device.type, dtype=torch.bfloat16, enabled=use_amp):
                    pred = forward_pred(seq_ids, seq_vals, pert_flags, pad_mask)
                    val_loss += float(loss_fn(pred, target).item()) * len(b_idx)
                    n_v += len(b_idx)
        val_loss /= max(n_v, 1)
        if ep == 0 or (ep + 1) % 5 == 0 or ep == epochs - 1:
            print(f"    ep {ep+1}/{epochs}: tr_mse={tr_loss/max(n_steps,1):.4f} val_mse={val_loss:.4f} ({time.time()-t_ep:.1f}s)", flush=True)
        if val_loss < best_val:
            best_val = val_loss
            best_state = {
                "model": {k: v.detach().cpu().clone() for k, v in model.state_dict().items()},
                "head": None if full_head is None else {k: v.detach().cpu().clone() for k, v in full_head.state_dict().items()},
            }

    model.load_state_dict({k: v.to(device) for k, v in best_state["model"].items()})
    if full_head is not None and best_state["head"] is not None:
        full_head.load_state_dict({k: v.to(device) for k, v in best_state["head"].items()})
    model.eval()
    if full_head is not None:
        full_head.eval()

    q_ids = q_df["id"].tolist()
    q_syms = q_df["gene_symbol"].astype(str).tolist()
    preds = {}
    with torch.no_grad():
        for s in range(0, len(q_ids), batch_size):
            b_idx = list(range(s, min(s + batch_size, len(q_ids))))
            seq_ids, seq_vals, pert_flags, pad_mask = build_batch(b_idx, q_syms)
            with torch.autocast(device_type=device.type, dtype=torch.bfloat16, enabled=use_amp):
                out = forward_pred(seq_ids, seq_vals, pert_flags, pad_mask).float().cpu().numpy()
            for j, qi in enumerate(b_idx):
                preds[q_ids[qi]] = out[j].astype(np.float64)

    q_in_vocab = sum(1 for sym in q_syms if sym in vocab)
    del model, full_head, opt, best_state
    if device.type == "mps":
        torch.mps.empty_cache()
    gc.collect()
    return preds, {
        "readout_vocab_coverage": round(in_vocab_readout / len(gene_names), 4),
        "test_target_vocab_coverage": round(q_in_vocab / max(len(q_syms), 1), 4),
    }


def run_probe(dataset: str, seeds: list[int]):
    """Zero-shot pretrained scGPT_human gene embeddings + Ridge(alpha=10.0) across all splits."""
    t0 = time.time()
    weights_dir = resolve_scgpt_dir()
    vocab = json.load(open(weights_dir / "scgpt_vocab.json"))
    sd = torch.load(weights_dir / "scgpt_gh_repo_original_model.bin", map_location="cpu", weights_only=False)
    emb_mat = sd["encoder.embedding.weight"].float().numpy()  # (60697, 512)
    pad_id = vocab["<pad>"]

    print(f"=== R7-2 scGPT_human Pretrained Gene-Embedding Probe (Ridge alpha=10) on {dataset} ===", flush=True)
    train_df, all_y, sh0, sh1, g2i_map, sidx, _ = load_bundle(dataset)
    model = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)
    rows, ci_rows = [], []

    for seed in seeds:
        val_df, df_tr, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
        df_hct, df_hek = df_tr[df_tr["source"] == "HCT116"], df_tr[df_tr["source"] == "HEK293T"]
        q_hct, q_hek = val_df[val_df["source"] == "HCT116"], val_df[val_df["source"] == "HEK293T"]

        p_scgpt, p_ae = {}, {}
        covs = []
        for src, df_s, q_s in [("HCT116", df_hct, q_hct), ("HEK293T", df_hek, q_hek)]:
            g2i = g2i_map[src]
            own_rep = df_s[df_s["reproducible"] == 1]
            tr_syms = own_rep["gene_symbol"].astype(str).tolist()
            q_syms = q_s["gene_symbol"].astype(str).tolist()
            X_tr = np.array([emb_mat[vocab.get(s, pad_id)] for s in tr_syms], dtype=np.float32)
            X_q = np.array([emb_mat[vocab.get(s, pad_id)] for s in q_syms], dtype=np.float32)
            Y_tr = np.array([all_y[rid] for rid in own_rep["id"]], dtype=np.float32)
            p_arr = Ridge(alpha=10.0).fit(X_tr, Y_tr).predict(X_q)
            for i, rid in enumerate(q_s["id"].tolist()):
                p_scgpt[rid] = p_arr[i].astype(np.float64)
            covs.append(round(sum(1 for s in q_syms if s in vocab) / max(len(q_syms), 1), 4))

            tr_sub = own_rep[own_rep["perturbation"].isin(g2i.keys())]
            tr_idx = [g2i[p] for p in tr_sub["perturbation"].tolist()]
            q_idx = [g2i[p] for p in q_s["perturbation"].tolist()]
            tr_Y_raw = np.array([all_y[rid] for rid in tr_sub["id"]], dtype=np.float32)
            Y_all_raw = np.array([all_y[rid] for rid in df_s["id"]], dtype=np.float32)
            arr_ae = run_ahlmann_eltze_bilinear(Y_all_raw, tr_idx, q_idx, tr_Y_raw, rank=100, lam=10.0)
            for i, rid in enumerate(q_s["id"].tolist()):
                p_ae[rid] = arr_ae[i]

        r0 = model.predict_for_source_all_stages("HCT116", df_hct, df_hek, q_hct, all_y, bases, axes, ("S5", "S7"))
        r1 = model.predict_for_source_all_stages("HEK293T", df_hek, df_hct, q_hek, all_y, bases, axes, ("S5", "S7"))
        preds_map = {"scGPT_FrozenEmb_Ridge": p_scgpt, "S5": {**r0["S5"], **r1["S5"]}, "S7": {**r0["S7"], **r1["S7"]}, "AE": p_ae}
        rm = compute_row_metrics_for_split(preds_map, val_df, all_y, bases, axes, sidx, g2i_map)
        cm = compute_centroid_rows(preds_map, df_tr, val_df, all_y, sidx, g2i_map)
        src_arr = val_df["source"].values
        i0, i1 = np.where(src_arr == "HCT116")[0], np.where(src_arr == "HEK293T")[0]
        v5_map = {m: row_scores_v5(preds_map[m], val_df, all_y, bases, axes, sidx, g2i_map) for m in preds_map}
        ev = evaluate_5metrics(p_scgpt, val_df, all_y, bases, axes, sidx, g2i_map)
        rows.append({
            "dataset": dataset, "seed": seed, "method": "scGPT_FrozenEmb_Ridge",
            **{k: round(v, 6) for k, v in ev.items()},
            "disc_raw": round(float(0.5 * (rm["scGPT_FrozenEmb_Ridge"]["disc_raw"][i0].mean() + rm["scGPT_FrozenEmb_Ridge"]["disc_raw"][i1].mean())), 6),
            "disc_v5": round(float(0.5 * (rm["scGPT_FrozenEmb_Ridge"]["disc_v5"][i0].mean() + rm["scGPT_FrozenEmb_Ridge"]["disc_v5"][i1].mean())), 6),
            "centroid_train_r": round(float(0.5 * (cm["scGPT_FrozenEmb_Ridge"]["centroid_train_r"][i0].mean() + cm["scGPT_FrozenEmb_Ridge"]["centroid_train_r"][i1].mean())), 6),
            "vocab_cov_line0": covs[0], "vocab_cov_line1": covs[1],
        })
        for comp, ka, kb in [("S5_vs_scGPT_Probe", "S5", "scGPT_FrozenEmb_Ridge"), ("S7_vs_scGPT_Probe", "S7", "scGPT_FrozenEmb_Ridge"), ("scGPT_Probe_vs_AE", "scGPT_FrozenEmb_Ridge", "AE")]:
            ci_rows.append({"dataset": dataset, "seed": seed, "metric": "v5", "comparison": comp, **stratified_row_bootstrap_ci(v5_map[ka], v5_map[kb], src_arr)})
            for met in ["delta_r", "disc_raw", "disc_v5"]:
                ci_rows.append({"dataset": dataset, "seed": seed, "metric": met, "comparison": comp, **stratified_row_bootstrap_ci(rm[ka][met], rm[kb][met], src_arr)})

    df_m, df_ci = pd.DataFrame(rows), pd.DataFrame(ci_rows)
    tag = f"r72_scgpt_probe_{dataset}"
    df_m.to_csv(REPORTS_DIR / f"{tag}_metrics.csv", index=False)
    df_ci.to_csv(REPORTS_DIR / f"{tag}_row_ci.csv", index=False)
    print(f"  [scGPT_FrozenEmb_Ridge {dataset}] v5={df_m['v5'].mean():.4f} delta_r={df_m['delta_r'].mean():.4f} disc_raw={df_m['disc_raw'].mean():.4f} ({time.time()-t0:.1f}s)", flush=True)


def run_finetune(dataset: str, seeds: list[int], device_str: str, epochs: int, batch_size: int, max_genes_ctx: int):
    t0 = time.time()
    device = torch.device(device_str)
    weights_dir = resolve_scgpt_dir()
    vocab = json.load(open(weights_dir / "scgpt_vocab.json"))
    print(f"=== R7-2 Official scGPT_human TransformerGenerator Fine-Tuning on {dataset} | seeds={seeds} | device={device} | ctx={max_genes_ctx} | ep={epochs} ===", flush=True)
    train_df, all_y, sh0, sh1, g2i_map, sidx, gnames_map = load_bundle(dataset)
    model_pan = PanPerturbBridgeModel(g2i_map, sh0, sh1, sidx)
    rows, ci_rows = [], []

    for seed in seeds:
        ts = time.time()
        val_df, df_tr, _, bases, axes = make_strict_union_split(train_df, g2i_map, all_y, sidx, seed)
        df_hct, df_hek = df_tr[df_tr["source"] == "HCT116"], df_tr[df_tr["source"] == "HEK293T"]
        q_hct, q_hek = val_df[val_df["source"] == "HCT116"], val_df[val_df["source"] == "HEK293T"]

        p0, c0 = run_scgpt_finetune_one_line(df_hct, q_hct, all_y, bases["HCT116"], gnames_map["HCT116"], vocab, weights_dir, seed, device, epochs=epochs, batch_size=batch_size, max_genes_ctx=max_genes_ctx)
        p1, c1 = run_scgpt_finetune_one_line(df_hek, q_hek, all_y, bases["HEK293T"], gnames_map["HEK293T"], vocab, weights_dir, seed, device, epochs=epochs, batch_size=batch_size, max_genes_ctx=max_genes_ctx)
        p_scgpt = {**p0, **p1}

        r0 = model_pan.predict_for_source_all_stages("HCT116", df_hct, df_hek, q_hct, all_y, bases, axes, ("S5", "S7"))
        r1 = model_pan.predict_for_source_all_stages("HEK293T", df_hek, df_hct, q_hek, all_y, bases, axes, ("S5", "S7"))
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

        preds_map = {"Official_scGPT_FT": p_scgpt, "S5": {**r0["S5"], **r1["S5"]}, "S7": {**r0["S7"], **r1["S7"]}, "AE": p_ae}
        rm = compute_row_metrics_for_split(preds_map, val_df, all_y, bases, axes, sidx, g2i_map)
        cm = compute_centroid_rows(preds_map, df_tr, val_df, all_y, sidx, g2i_map)
        src_arr = val_df["source"].values
        i0, i1 = np.where(src_arr == "HCT116")[0], np.where(src_arr == "HEK293T")[0]
        v5_map = {m: row_scores_v5(preds_map[m], val_df, all_y, bases, axes, sidx, g2i_map) for m in preds_map}
        ev = evaluate_5metrics(p_scgpt, val_df, all_y, bases, axes, sidx, g2i_map)
        rows.append({
            "dataset": dataset, "seed": seed, "method": "Official_scGPT_FT",
            **{k: round(v, 6) for k, v in ev.items()},
            "disc_raw": round(float(0.5 * (rm["Official_scGPT_FT"]["disc_raw"][i0].mean() + rm["Official_scGPT_FT"]["disc_raw"][i1].mean())), 6),
            "disc_v5": round(float(0.5 * (rm["Official_scGPT_FT"]["disc_v5"][i0].mean() + rm["Official_scGPT_FT"]["disc_v5"][i1].mean())), 6),
            "centroid_train_r": round(float(0.5 * (cm["Official_scGPT_FT"]["centroid_train_r"][i0].mean() + cm["Official_scGPT_FT"]["centroid_train_r"][i1].mean())), 6),
            "vocab_cov_line0": c0["test_target_vocab_coverage"], "vocab_cov_line1": c1["test_target_vocab_coverage"],
        })
        for comp, ka, kb in [("S5_vs_Official_scGPT_FT", "S5", "Official_scGPT_FT"), ("S7_vs_Official_scGPT_FT", "S7", "Official_scGPT_FT"), ("Official_scGPT_FT_vs_AE", "Official_scGPT_FT", "AE")]:
            ci_rows.append({"dataset": dataset, "seed": seed, "metric": "v5", "comparison": comp, **stratified_row_bootstrap_ci(v5_map[ka], v5_map[kb], src_arr)})
            for met in ["delta_r", "disc_raw", "disc_v5"]:
                ci_rows.append({"dataset": dataset, "seed": seed, "metric": met, "comparison": comp, **stratified_row_bootstrap_ci(rm[ka][met], rm[kb][met], src_arr)})

        pd.DataFrame(rows).to_csv(REPORTS_DIR / f"r72_official_scgpt_ft_{dataset}_metrics.csv", index=False)
        pd.DataFrame(ci_rows).to_csv(REPORTS_DIR / f"r72_official_scgpt_ft_{dataset}_row_ci.csv", index=False)
        print(
            f"  [Official_scGPT_FT {dataset} seed={seed}] v5={ev['v5']:.4f} delta_r={ev['delta_r']:.4f} "
            f"disc_raw={rows[-1]['disc_raw']:.4f} cov=({c0['test_target_vocab_coverage']:.2f},{c1['test_target_vocab_coverage']:.2f}) [{time.time()-ts:.1f}s]",
            flush=True,
        )

    df_m, df_ci = pd.DataFrame(rows), pd.DataFrame(ci_rows)
    summary = {
        "dataset": dataset, "model": "Official scGPT_human TransformerGenerator (Cui et al., Nature Methods 2024)",
        "checkpoint_sha256": "ad0252a1971e0cd619b7116dbab3177432236c4537225d54280a2aa7e5fe402a",
        "vocab_sha256": "ee2b2c90158eedb97c2318e49abaaed0a02c6fdf7e3f7ca6a821906413c4d2a4",
        "hyperparameters": {"d_model": 512, "nhead": 8, "nlayers": 12, "max_genes_ctx": max_genes_ctx, "epochs": epochs, "batch_size": batch_size, "lr": 1e-4},
        "seeds": seeds,
        "metrics": {col: {"mean": round(float(df_m[col].mean()), 6), "std": round(float(df_m[col].std(ddof=1)), 6) if len(df_m) > 1 else 0.0} for col in ["v5", "delta_r", "de50_r", "de50_sign", "disc_raw", "disc_v5", "centroid_train_r"]},
        "comparisons": {},
    }
    for met in ["v5", "delta_r", "disc_raw", "disc_v5"]:
        summary["comparisons"][met] = {}
        for comp in ["S5_vs_Official_scGPT_FT", "S7_vs_Official_scGPT_FT", "Official_scGPT_FT_vs_AE"]:
            sub = df_ci[(df_ci["metric"] == met) & (df_ci["comparison"] == comp)]
            summary["comparisons"][met][comp] = {
                "mean_delta": round(float(sub["delta"].mean()), 6),
                "mean_ci95_lo": round(float(sub["ci95_lo"].mean()), 6),
                "mean_ci95_hi": round(float(sub["ci95_hi"].mean()), 6),
                "splits_ci_gt_zero": f"{int(sub['ci_gt_zero'].sum())}/{len(sub)}",
                "splits_ci_lt_zero": f"{int((sub['ci95_hi'] < 0).sum())}/{len(sub)}",
            }
    with open(REPORTS_DIR / f"r72_official_scgpt_ft_{dataset}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved r72_official_scgpt_ft_{dataset}_* in {time.time()-t0:.1f}s!", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["replogle", "nadig"], required=True)
    ap.add_argument("--mode", choices=["probe", "finetune"], required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=SEEDS_10)
    ap.add_argument("--device", default="mps")
    ap.add_argument("--epochs", type=int, default=10)
    ap.add_argument("--batch-size", type=int, default=4)
    ap.add_argument("--max-genes-ctx", type=int, default=1536)
    args = ap.parse_args()
    if args.mode == "probe":
        run_probe(args.dataset, args.seeds)
    else:
        run_finetune(args.dataset, args.seeds, args.device, args.epochs, args.batch_size, args.max_genes_ctx)
