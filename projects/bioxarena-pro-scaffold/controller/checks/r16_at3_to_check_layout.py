"""Controller: repackage one A-T3 audit npz ({model}_{ckpt}__{line}.npz) into the directory layout r16_check.py reads,
so the same controller code scores both parties. Pure key renaming; arrays are copied bit-for-bit.

Usage: python3 -I r16_at3_to_check_layout.py <at3_line_npz> <out_dir>
Prints the seen and unseen source lists stored by A-T3 (pass them to r16_check.py).
"""
import json
import sys
from pathlib import Path

import numpy as np

src, out = Path(sys.argv[1]), Path(sys.argv[2])
z = np.load(src, allow_pickle=False)
out.mkdir(parents=True, exist_ok=True)
obs = {k: z[k] for k in ["delta_obs", "delta_half_a", "delta_half_b", "half_valid", "delta_transfer", "source_grand_mean"]}
obs["target_indices"] = z["target_gene_idx"]
obs["fct_delta_pred"] = z["delta_fct"]
obs["gene_names"] = z["eval_genes"]
obs["target_names"] = z["eval_perts"]
np.savez(out / "observed.npz", **obs)
np.savez(out / "preds_STATE_SE_fixed_index_zero_matched.npz", delta_pred=z["delta_pred_matched"])
seen, unseen = [str(s) for s in z["seen_sources"]], [str(s) for s in z["unseen_sources"]]
for s in seen + unseen:
    np.savez(out / f"preds_STATE_SE_fixed_index_zero_mis_{s}.npz", delta_pred=z[f"delta_pred_mismatched__{s}"])
if "delta_pred_resample" in z.files:
    np.savez(out / "preds_resample.npz", delta_pred=z["delta_pred_resample"])
print(json.dumps({"model": str(z["model_key"]), "line": str(z["target_line"]), "seen": seen, "unseen": unseen,
                  "N": int(z["delta_obs"].shape[0]), "G": int(z["delta_obs"].shape[1])}))
