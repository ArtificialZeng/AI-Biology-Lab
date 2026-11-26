# Round 16 Group 3 (`GSE295214`, `H1 hESC`) Counts-Only QC Gate

- **Status:** `PASS_GROUP3_PRIMARY`
- **`N_eval^R16` (`r9_6 \ {vcc2025_h1} ∩ P_model`):** `406` (`>= 50` required)
- **`|G_eval^R9|`:** `6131` (`>= 3,000` required); **`|G_eval^R16|`:** `1837`
- **Per-Half NTC (`non-targeting_*`):** Half A = `28929`, Half B = `29224` (`>= 15` each required)

```json
{
  "line": "gse295214_h1",
  "qualifies_group3_primary": true,
  "n_eval_p5_r9_excl_vcc2025_h1": 427,
  "n_eval_r16": 406,
  "n_genes_g_eval_r9": 6131,
  "n_genes_g_eval_r16": 1837,
  "n_kept_half_a_ntc": 28929,
  "n_kept_half_b_ntc": 29224,
  "per_half_ntc_ge15": true,
  "adapter_audit": {
    "n_channels": 85,
    "channels_sorted_sample": [
      "GSM8943713",
      "GSM8943714",
      "GSM8943715",
      "GSM8943716",
      "GSM8943717",
      "...",
      "GSM8943795",
      "GSM8943796",
      "GSM8943797"
    ],
    "half_a_n_channels": 43,
    "half_b_n_channels": 42,
    "n_raw_cells": 1521921,
    "n_hq_rna_cells": 1499727,
    "qc_spec": {
      "min_umi": 300,
      "min_ngenes": 200,
      "max_pct_mt": null
    },
    "n_guide_features": 6638,
    "n_ntc_guides": 370,
    "n_targeting_guides": 6257,
    "n_excluded_guides": 11,
    "median_guide_threshold": 5.0,
    "n_kept_ntc": 58153,
    "n_kept_targeting": 834203,
    "n_excluded_guide_singlet_cells": 1100,
    "n_kept_half_a_ntc": 28929,
    "n_kept_half_b_ntc": 29224,
    "per_half_ntc_ge15": true,
    "n_low_ntc_channel_fallback_cells": 21,
    "half_rule": "PREREG_H_Tier1_lexicographical_GSM_1based_odd_HalfA_even_HalfB"
  }
}
```
