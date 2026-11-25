# Round 16 Pre-Data Addendum A (`PREREG_I_ADDENDUM_A.md`, Items 1–11)

> **STATUS: PRE-DATA ADDENDUM A (DRAFT FOR MASTER CONTROLLER LOCK PRIOR TO `CODE_FREEZE_R16.txt` AND `R16_GO.txt`).**
> **Locked Parent Protocol (Unmodified):** [`round16/prereg/PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md) (`v5-FINAL-PRE-LOCK`, commit `c1fc7c6d24dd128bc5e2d969d8c78224a2cbc118`, SHA-256 `43d242e460fd64db2e4af1e628e886cf1d94fb4ccd1e95925928a99cb96ee4a8`, pushed to `origin/main` `2026-11-26 00:24:20 +08:00`; verified in [`PREREG_I_COMMIT.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PREREG_I_COMMIT.txt) and [`PREREG_I_LOCK_CHECK.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PREREG_I_LOCK_CHECK.txt); GCS `virtual-cell-challenge/2026/` verified empty at `2026-11-26 00:23:53 +08:00`).
> **Pre-Data & Blinding Attestation:** Drafted strictly **before** code freeze (`CODE_FREEZE_R16.txt`), **before** execution release (`R16_GO.txt`), **before** any real-cell inference (`0` real cells inferred), **before** downloading Arc archived `adata_pred.h5ad`, and **before** downloading `GSE295214`. `PREREG_I.md` itself is untouched. Does not weaken or modify primary `H-I1`/`H-I2` hypotheses, Bonferroni $\alpha_{\mathrm{claim}} = 0.025$, or TOST $\delta = 0.04$ (beyond closing the exact-identical-array boundary case in Item 10). Zero priority claims are made.

- **Authorizing Rulings:**
  - Master Controller (`Claude`) rulings `2026-11-26 00:25 +08:00` (Items 1–7, resolving all 7 post-lock items from `A-T3`'s pre-lock final audit [`PREREG_I_V5_AUDIT_AT3.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PREREG_I_V5_AUDIT_AT3.md));
  - Master Controller (`Claude`) supplemental rulings `2026-11-26 00:28 +08:00` (Items 8–10, addressing synthetic verification edge cases observed in `controller/checks/r16_check.py`);
  - Master Controller (`Claude`) pre-freeze reproducibility directive `2026-11-26 01:01 +08:00` (cross-process `PYTHONHASHSEED` determinism gate recorded in `§7`);
  - Master Controller (`Claude`) Item 11 rulings `2026-11-26 01:31 & 01:32 +08:00` (unified official `arc-state 0.11.1` `SE-600M` collator, per-pool seeding, `X_state` embedding export + SHA-256 cross-party comparison `<= 1e-5`, and cross-process embedding determinism recorded in `§7` and `§11`).
- **Scope:** Must be locked on `origin/main` strictly prior to `CODE_FREEZE_R16.txt` and any real-cell inference. Every specification below matches `round16/src/` verbatim.

---

## §1. Group 3 (`GSE295214`) `full r9_6` Secondary Sensitivity Locked to Primary `N_eval^R16` (`PREREG_I §2.4.4`)

1. **Fixed Target and Gene Universe:** The Group 3 secondary lineage-overlap sensitivity (`full r9_6`, including `vcc2025_h1` alongside the 5 primary sources `{nadig_hepg2, nadig_jurkat, nourreddine_kolf21j, replogle_k562_essential, replogle_rpe1}`) MUST be evaluated on the **exact same target set $\mathcal{T}_{\mathrm{eval}}^{\mathrm{R16}}(\text{gse295214\_h1})$** (`N_eval^R16` determined by the primary 5-line pool `r9_6 \ {vcc2025_h1}`) and the **exact same gene universe $G_{\mathrm{eval}}^{\mathrm{R16}}(\text{gse295214\_h1})$**. No extra targets introduced solely by `vcc2025_h1` are added.
2. **Role of `vcc2025_h1`:** Adding `vcc2025_h1` affects ONLY the computation of `FCT` (`FctSourceState.fit` $\to$ $\overline{\mathrm{FC}}_p$, $\hat{\Delta}^{\mathrm{FCT}}_p$), `Transfer_Mean` ($\bar{\Delta}^{\mathrm{src}}_p$), and `Source_Grand_Mean` ($\bar{b}_S$) on those locked `(targets, genes)`.
3. **Evaluation & Audit Export:** Report `FCT`, `STATE_SE`, and paired `STATE_SE - FCT` under both `r_orth_cal` and `r_dorth_cal` using the `full r9_6` `Transfer_Mean` (and `Source_Grand_Mean`) strictly as a secondary sensitivity alongside primary `r9_6 \ {vcc2025_h1}`. Export `fct_delta_pred`, `delta_transfer`, and `source_grand_mean` for `FCT_full_r9_6_incl_vcc2025_h1` in `round16/results/arrays/gse295214_h1/` for independent recomputation.
4. **Code Correspondence:** `round16/src/run_eval_r16.py::_align_line_and_sources_r16(..., fixed_target_names=q_names, fixed_gene_names=g_eval)` and `evaluate_line_from_arrays_r16(..., extra_fct_sensitivity=...)`.

---

## §2. Group 3 (`GSE295214`) Exact Guide Parsing, Duplicate-Feature Deduplication, and `P5` Thresholds (`PREREG_I §2.4.2`)

All rules below reproduce `round16/src/run_eval_r16.py::build_gse295214_parts`, `parse_gse295214_guide_r16`, `_resolve_single_guide_token`, and frozen `round10/shared/scripts/prepare_r10_lines.py::select_canonical_genes` & `process_r10_line` verbatim.

### §2.1 Guide Name Parsing (`parse_gse295214_guide_r16` & `_resolve_single_guide_token`)
For every 10x `.h5` feature with `matrix/features/feature_type != "Gene Expression"`, let `gid = str(features/id).strip()` and `gnm = str(features/name).strip()`:
1. **NTC Classification (Highest Precedence):** If `gid.lower().startswith("non-targeting_")` OR `gnm.lower().startswith("non-targeting_")`, classify strictly as `"ntc"`.
2. **Primary String Selection:** Let `raw = gnm or gid`; if empty, classify as `"exclude"`.
3. **Pipe-Delimited Dual/Multi-Guide Cassettes (`|`):** If `"|"` occurs in `raw`, split on `"|"` and strip whitespace into non-empty substrings `subs`:
   - If `subs` is empty $\to$ `"exclude"`;
   - If every substring in `subs` starts with `"non-targeting_"` (case-insensitive) $\to$ `"ntc"`;
   - Otherwise resolve every substring via `_resolve_single_guide_token(s, hgnc)`. If any substring returns `None` OR the set of resolved canonical HGNC gene symbols across `subs` has size $\neq 1$ (multi-gene cassette), classify as `"exclude"`. Only when all substrings resolve to the **single identical approved HGNC gene** is the feature classified as `"target"` (retaining the leading resolved tuple `resolved[0]`).
4. **Single-Token Resolution (`_resolve_single_guide_token(token, hgnc)`):**
   - Strip whitespace; reject (`None`) if empty, matching multi-target/illegal delimiter regex `[,;\s]` (comma, semicolon, or whitespace), OR matching illegal `+` regex `(?<!_)\+|\+(?![_\d])` (rejecting multi-gene constructs such as `TP53+MYC` or `TP53_1+MYC_1` while preserving Weissman `+` strand tokens such as `TP53_+_7687490.23-P1P2` and signed coordinate offsets such as `TP53_+12345_P1`; any dual-gene token of the form `TP53_+_MYC` is still rejected by the stage-2 multi-gene guard below).
   - Lowercase `low = tok.lower()`:
     - Reject (`None`) if `low.startswith("non-targeting_")` or `low == "non-targeting"`.
     - **Safe-Targeting, Intergenic, Multiplet, & Multi-Gene Blacklist:** Reject (`None`) if `low` contains any substring in `("safe", "intergenic", "unassigned", "multiplet", "doublet", "gene_desert")` OR contains `"-and-"`.
   - Extract optional embedded Ensembl gene ID via case-insensitive regex `(ENSG\d+)` (`ens_raw`).
   - Helper `_try_resolve(sym_cand, ens_cand)` calls frozen `HGNCLookup.resolve(sym_cand, ens_cand)` (precedence: `approved -> ENSG -> unambiguous prev_symbol -> unambiguous alias_symbol`, or equivalent `HGNCBridge.canonicalize_with_source`) and **accepts ONLY symbols present in HGNC's `approved` set**.
   - Three-stage deterministic lookup (returns upon earliest hit):
     1. **Exact full-token lookup:** `_try_resolve(tok, ens_raw)` (protects hyphenated HGNC symbols such as `HLA-A`, `MT-CO1`).
     2. **Underscore `_` delimited tokens (`len(parts) >= 2`, dropping empty segments):**
        - *Multi-gene guard:* Ignore tokens matching the case-insensitive guide-suffix regex
          `^(?:P\d+(?:[-]?P\d+)*|ENST\d+|ENSG\d+|\d+(?:\.\d+)?(?:-P\d+.*)?|sg\d+|TSS\d*|[+-]|[+-]\d+.*|AB|A|B|[ACGT]{15,})$`
          (covering promoter indices `P1`/`P1P2`/`P1-P2` without matching `P2RX7`/`P4HB`-like HGNC gene symbols, `ENST*`/`ENSG*`, integer/Weissman genomic coordinate tokens such as `7687490.23-P1P2`, `sg\d+`, `TSS\d*`, Weissman strand flags `+`/`-`, signed coordinate offsets, `AB`/`A`/`B`, and $\ge 15$ nt protospacer sequences). Resolve every remaining segment via `_try_resolve(p, "")`; if `> 1` distinct canonical HGNC gene symbols match (e.g., `TP53_MYC_1`, `TP53_+_MYC`, `TP53_P2RX7_1`), reject (`None`).
        - *Candidate resolution:* If `parts[0]` starts with `"ENSG"` (case-insensitive), try `_try_resolve(parts[1], parts[0].split(".")[0])`; next try `_try_resolve(parts[0], ens_raw)` (handling `GENE_1`, `GENE_P1P2_1`, `GENE_TSS1_1`, `GENE_+_7687490.23-P1P2`, `HLA-A_P1P2_1`, `ENSG..._GENE_1`).
     3. **Trailing hyphen `-` suffix strip (only if stages 1–2 did not match):** Match case-insensitive regex `^(.+?)-(?:\d+|sg\d+|P\d+(?:[-]?P\d+)*|TSS\d*)$` (stripping trailing `-<int>`, `-sg<int>`, `-P<int>...`, or `-TSS[<int>]` such as `GENE-1`, `GENE-sg1`, `GENE-P1P2`, `GENE-TSS1`) and call `_try_resolve(prefix, ens_raw)`.
5. **ID Fallback:** If `raw` (`gnm`) returns `None` and `gid` is non-empty with `gid != raw`, run `_resolve_single_guide_token(gid, hgnc)` once more; otherwise `"exclude"`.
6. **Per-Cell Guide Key:** Each QC-passing targeting singlet (`n_pos == 1`, `top1 >= 5`, `top2 < 3 or top1 / max(top2, 0.5) >= 3.0`) is tagged with its exact guide key `gkey = crp_names_ref[g_i] or crp_ids_ref[g_i]` and passed into frozen `process_r10_line` for per-guide Addendum-B testing followed by canonical gene deduplication.

### §2.2 Duplicate RNA Feature & Duplicate Target-Gene Deduplication
1. **RNA Feature Deduplication (`prepare_r10_lines.py::select_canonical_genes`, `max_genes = None`):**
   - Retain raw `Gene Expression` features passing finite-expression QC across all batches with `>= 15` NTCs (`mean` finite, `sd > 1e-12`) AND global NTC `libnorm(1e4)` `ctrl_mean_lib > 0` and `sd_global > 1e-12`.
   - Canonicalize each feature via `HGNCLookup.resolve(raw_sym, raw_ens)`. When multiple raw features map to the same canonical gene symbol, **retain the single feature with strictly highest global NTC `libnorm(1e4)` mean `ctrl_mean_lib`** (`score > best_score`, breaking exact ties by smallest original raw feature index). Retained features are sorted by original raw feature index ascending.
   - *Control-tensor extraction for `ST-HVG` (`model_adapters_r16.py::preprocess_hvg_counts_r16`):* Normalization factors (`official_median_log1p` or `libnorm1e4_log1p`) are computed over ALL raw `Gene Expression` features in the sampled control pool; when mapping raw columns via `HGNCBridge` onto `G_model`, duplicate canonical symbols retain the column with strictly highest mean raw UMI count in that sampled control pool (exact ties broken by smallest raw column index).
2. **Duplicate Target-Gene Deduplication (`prepare_r10_lines.py::process_r10_line`):**
   - Per-guide groups with `< 10` total QC-passing cells are discarded upfront (`len(cells) >= 10` retained).
   - When multiple guides map to the same canonical HGNC target gene symbol, retain the single guide maximizing the lexicographic tuple **`(int(reproducible), float(ad_counts), float(n_cells))`** (strict `>` comparison; exact ties retain the lexicographically smallest `pert_key` because `all_pkeys` is sorted ascending).

### §2.3 Frozen `P5` Reproducibility Thresholds (`prepare_r10_lines.py::process_r10_line`)
Inherited unchanged from frozen Rounds 10–12 (`PREREG_E §3.3`):
1. **Minimum Perturbation Cell Count:** `n_cells >= 25` QC-passing single-guide cells (`< 25` automatically receives `reproducible = False`, `ad_counts = 0.0`, `energy_test_p_value = 1.0`).
2. **Feature Space for Addendum-B QC:** Top `5,000` finite genes ranked by global NTC `ctrl_mean_lib`, batch-standardized to z-scores using per-batch NTC `(mu_lib, sd_lib)` (falling back to global NTC SD if batch NTC `< 15` or `sd <= 1e-12`).
3. **Anderson–Darling DE Gene Count:** Two-sample Anderson–Darling test vs. `1,000` batch-stratified NTC cells (`stratified_nt_sample(..., n_target=1000, seed=42)`) across the top `5,000` genes with Benjamini–Hochberg FDR correction: require **`ad_counts = sum(q_AD < 0.05) >= 50`**.
4. **20-PC Energy-Distance Permutation Test:** Project the top `5,000` z-scores onto the leading `20` eigenvectors of global $Z^\top Z$, compare against `500` batch-stratified NTC cells (`stratified_nt_sample(..., n_target=500, seed=42)`) using `n_perm = 2,000` permutations (`RandomState((42 + pi * 1000003) & 0x7FFFFFFF)`): require **`energy_test_p_value < 0.001`**.
5. **Reproducibility & Split-Half Validity:** `reproducible = (n_cells >= 25) AND (ad_counts >= 50) AND (energy_test_p_value < 0.001)`; `half_valid = (half_a_n_cells >= 2) AND (half_b_n_cells >= 2)`.
6. **Qualifying Floor (`R16_GSE295214_QC_COUNTS.md`, logged pre-metric):** `N_eval^R16 >= 50` (`reproducible` on `gse295214_h1`, `reproducible` in $\ge 1$ source of `r9_6 \ {vcc2025_h1}`, self-gene in `G_eval^R9`, canonical symbol in `canon(P_model)`), `|G_eval^R9| >= 3,000`, AND `per-half NTC >= 15` in both Half A and Half B.

---

## §3. Arc Archived `adata_pred.h5ad` Adapter Fidelity Check (`PREREG_I §1.3` & `§5.5`)

1. **Fixed Single Check Line:** Because human download approval covers only one Group 1 archived prediction file, the adapter check is **strictly fixed to `nadig_hepg2` ONLY** (`zeroshot/hepg2/eval_final.ckpt/adata_pred.h5ad`). `nadig_jurkat`'s archived `adata_pred.h5ad` will not be downloaded.
2. **Download Ordering & Control:** `adata_pred.h5ad` is downloaded strictly **after `R16_GO.txt`**, under a download manifest controlled by Master Controller.
3. **Comparison Protocol (`run_eval_r16.py::verify_arc_archived_adata_pred_r16` & `compare_predicted_deltas_vs_arc_archive_r16`):**
   - Run our `ArcStateAdapter("STATE_SE", "zeroshot/hepg2", batch_mode="native_gem_group")` on `nadig_hepg2`'s sampled `512` NTC cells using each cell's native `gem_group` mapped through `batch_onehot_map.torch` (identical to `§4.4` cell `(M, B_M)`).
   - From Arc's archived `adata_pred.h5ad`, compute per-perturbation pseudobulk shift $\hat{\Delta}^{\mathrm{arc}}_p = \overline{Y}_{\mathrm{arc}}(p) - \overline{Y}_{\mathrm{arc}}(\text{"non-targeting"})$.
   - Align strictly on all common non-control perturbations ($P_{\mathrm{common}} = \mathcal{T}_{\mathrm{eval}}^{\mathrm{R16}}(\text{nadig\_hepg2}) \cap P_{\mathrm{arc}}$) and common evaluated genes ($G_{\mathrm{eval}}^{\mathrm{R16}}(\text{nadig\_hepg2})$), and compute:
     - Per-perturbation Pearson $r_p = \operatorname{corr}(\hat{\Delta}^{\mathrm{adapter}}_p, \hat{\Delta}^{\mathrm{arc}}_p)$ across genes, taking the median across $p \in P_{\mathrm{common}}$ (`median_per_pert_pearson`);
     - Pooled Pearson correlation across all flattened `(p, g)` pairs (`pooled_pearson`).
4. **Locked Pass Criteria (Both Required):**
   $$\mathrm{median\_per\_pert\_pearson} \ge 0.90 \quad \text{AND} \quad \mathrm{pooled\_pearson} \ge 0.95.$$
5. **Failure Protocol:** If either criterion fails, **halt immediately prior to final adjudication** and escalate to Master Controller to diagnose the adapter implementation. **No adjudication rule, threshold, or hypothesis may be modified.**

---

## §4. `STATE_HVG` `Perm100` Null Definition and `G_model` Column Mapping (`PREREG_I §2.1.4` & `§4.3.1`)

1. **`G_model` Column Mapping (`HGNCBridge`):** For `STATE_HVG` (`preprocess_hvg_counts_r16`), raw input gene names and the `2,000` ordered `G_model` column names are harmonized through frozen R10 `HGNCBridge` (`hgnc_complete_set.txt`; recovering legacy/updated pairs such as `AARS`/`AARS1`, `DARS`/`DARS1`, `H2AFZ`/`H2AZ1`, `TAZ`/`TAFAZZIN`), while `SE-600M` retains exact gene-symbol string lookup into `protein_embeddings.pt` per `PREREG_I §2.1.4`.
2. **`STATE_HVG` `Perm100` Null Definition (`metrics_r16.py::compute_gene_label_perm_null_r16`):**
   - Follows frozen `PREREG_E §4.1` (`k = 0..99`) adapted to predicted delta matrices within $G_{\mathrm{eval}}^{\mathrm{R16}}(t)$ (let $G = |G_{\mathrm{eval}}^{\mathrm{R16}}(t)|$).
   - For each $k \in \{0, 1, \dots, 99\}$:
     - Generate deterministic permutation $\pi_k = \texttt{np.random.default\_rng(k).permutation(G)}$ (seed strictly equal to `k`, matching `PREREG_E §4.1`);
     - Permute the gene columns of the `STATE_HVG` matched-arm predicted delta matrix (`official_median_log1p` primary normalization): $\hat{\Delta}^{(k)}_{p, g} = \hat{\Delta}^{\mathrm{HVG}}_{\mathrm{matched}, p, \pi_k(g)}$, holding `delta_obs`, `delta_transfer`, `target_indices`, `delta_half_a`, `delta_half_b`, `half_valid`, and `source_grand_mean` fixed;
     - Compute point-estimate `r_orth_cal` (and descriptively `r_dorth_cal`) on `all_nontarget`.
   - Report all `100` permutation scores, `median`, `mean`, `q025`, `q975`, `max`, and one-sided empirical tail proportion `(1 + #{perm >= obs}) / 101` strictly as a secondary descriptive analysis for every evaluated line.

---

## §5. `CellularAI_OpenCode_PublicRefit` Explicitly Dropped (`PREREG_I §1.2` & `§4.3.3`)

- **`CellularAI_OpenCode_PublicRefit` WILL NOT BE EXECUTED (`不做`).**
- This decision is locked pre-data (prior to code freeze and prior to any real-cell execution) to eliminate all operational degrees of freedom associated with the former `"optional and non-blocking"` label in `PREREG_I.md §1.2` and `§4.3.3`.
- In `round16/src/model_adapters_r16.py`, `CellularAIOpenCodePublicRefitAdapter` is marked `role = "dropped_per_addendum_a"` and raises `NotImplementedError` if invoked; `round16/src/run_eval_r16.py` never runs it.

---

## §6. `G4` 2×2 Factorial Pairing and Symmetry Check (`PREREG_I §4.4.3`)

1. **Group 1 (`nadig_hepg2`, `nadig_jurkat`):**
   - Target line $t$ and counterpart line $s$ each supply `512` sampled NTC cells ordered `0..511` in their deterministic `16 × 32` post-shuffle order (`seed = 20261601`), together with their respective native `gem_group` labels $B_M[0..511]$ and $B_S[0..511]$ mapped via `batch_onehot_map.torch`.
   - Cross-cells `(M, B_S)` and `(S, B_M)` align **position-by-position (`i = 0..511`)**: cell `i` of `X_ctrl` receives batch label `B[i]` (`1`-to-`1` positional pairing, zero reshuffling).
2. **Group 2 (`a549_ko`, `dld1`, `cd4t_rest`) & Group 3 (`gse295214_h1`, if approved/qualified):**
   - Primary `G4` pairing fixes **`B_M = 0` (`fixed_index_zero`), `B_S = 17` (`fixed_index:17`)** (`b0_vs_b17_mis_<s'>`), evaluated for every mismatched source $s' \in \mathcal{S}_{\mathrm{mis}}^{\mathrm{seen}}(t) \cup \mathcal{S}_{\mathrm{mis}}^{\mathrm{unseen}}(t)$.
   - **Mandatory `(17, 0)` Symmetry Check:** Additionally report the reversed pairing **`B_M = 17`, `B_S = 0`** (`b17_vs_b0_sym_mis_<s'>`) for every $s'$, using the exact same four already-computed cell matrices (`(M, 17)`, `(M, 0)`, `(S, 17)`, `(S, 0)`), verifying whether `delta_batch`, `delta_ctrl`, `delta_int`, and `batch_share` are invariant to swapping index order.
   - Report per-source `G4` contrasts individually alongside stratum-level summaries (`seen_by_checkpoint` and `unseen_by_checkpoint` separately), plus matched-arm single-index spread across `b ∈ {0, 1, 17, 35}`.

---

## §7. Three-Party Verification Tolerances and Mandatory Control-Barcode Export (`PREREG_I §5.4` & `§5.5`)

Because independent `MPS float32` runs on Apple Silicon are not guaranteed to be bit-for-bit identical across independent implementations (`round16/src/` vs. `round16/audit/`), three-party audit tolerances are locked as follows:

| Audit Layer | Parties Compared | Locked Acceptance Criterion |
|---|---|---|
| **1. Control-cell sampling** | Main Line vs. `A-T3` (for every control source, BOTH primary `seed=20261601` AND Item 8 resample `seed=20261602`) | **Exact `100%` ordered string equality** of sampled cell barcode lists |
| **2. `SE-600M` `X_state` embeddings** | Main Line vs. `A-T3` (independent `MPS float32` encoding per `§11`, checked **prior to** predicted deltas for every primary `seed=20261601` AND resample `seed=20261602` control pool) | $\max_{i, d} \left|X_{\mathrm{state}}^{\mathrm{Main}}[i, d] - X_{\mathrm{state}}^{\mathrm{AT3}}[i, d]\right| \le 10^{-5}$ (with exported `x_state_*.npz` and `embeddings_manifest_sha256.json`) |
| **3. Independent predicted delta arrays** | Main Line vs. `A-T3` (independent `MPS float32` inference) | $\max_{p, g} \left|\hat{\Delta}^{\mathrm{Main}}_{p, g} - \hat{\Delta}^{\mathrm{AT3}}_{p, g}\right| \le 10^{-4}$ |
| **4. Same-array metric recomputation** | Master Controller recomputing from Main Line's or `A-T3`'s exported `.npz` arrays | $\left|\mathrm{metric}^{\mathrm{Controller}} - \mathrm{metric}^{\mathrm{Party}}\right| \le 10^{-10}$ (all point estimates & bootstrap CIs under identical `seed = 20261601`, `n_boot = 4,000`) |
| **5. Cross-party metrics & verdicts** | Main Line vs. `A-T3` (from independent `MPS float32` arrays) | $\left|\mathrm{metric}^{\mathrm{Main}} - \mathrm{metric}^{\mathrm{AT3}}\right| \le 10^{-5}$ **AND `100%` verdict parity** (`H-I1` pair/stratum/line/group; `H-I2` line) |
| **6. Discrepancy reporting** | All three parties | **Every** verdict discrepancy (if any occurs) MUST be reported item-by-item with full root-cause disclosure; selective omission is strictly prohibited. |

### Locked Canonical Cell Barcode Format & Export Locations
To ensure exact `100%` ordered barcode parity and run-to-run determinism across processes and parties:
- **Pre-Freeze Cross-Process Determinism Requirement (`2026-11-26 01:01` & `01:31`):** No production code path may iterate an unordered container (`set`/`frozenset`) to determine sampling, array order, or model outputs; every stochastic operation derives strictly from pre-registered seeds (`20261601`, `20261602`, `Perm100` `k = 0..99`); unit tests write strictly inside `tempfile.TemporaryDirectory()` (never `round16/status/` or `round16/results/`); and independent fresh processes with distinct `PYTHONHASHSEED` (`1` vs. `99991`) MUST produce bit-for-bit identical (`max_abs == 0.0`) `SE-600M` `X_state` embeddings (`§11`, for both `seed=20261601` and `seed=20261602`) and predicted-delta arrays on CPU, and `max_abs <= 1e-6` on MPS (`test_unit_r16.py::test_18_cross_process_pythonhashseed_determinism`).
- Hamilton largest-remainder batch quotas in `sample_stratified_control_sets_r16` MUST use **exact integer arithmetic** (`quota_k = (target_n * n_k) // n_total`, remainder numerator `-(target_n * n_k - quota_k * n_total)` with lexicographic batch-label tie-break `k`), eliminating any `float64` rounding divergence.
- Every control source exports its ordered sampled cell barcode list (length `512`, or `(N_NTC // 32) * 32` when `< 512`) using these exact identifiers:
- **`nadig_hepg2`, `nadig_jurkat`:** exact AnnData `obs` index string (`obs.attrs["_index"]` $\to$ `"_index"` $\to$ `"index"`).
- **`a549_ko`:** exact `cell_id` column string from `GSE345058_cell_annotations.csv.gz`.
- **`feng_ipsc`:** exact `Cell_ID` column header string from `TargetedScreen_RNA-UMI-Counts.csv.gz`.
- **`dld1`:** `"{gsm}:{raw_barcode}"` (`gsm = filename.split("_")[0]`, `raw_barcode` from `matrix/barcodes`).
- **`cd4t_rest`:** `"{donor_id}_{lane_tag}:{raw_barcode}"` (`raw_barcode` from `matrix/barcodes`).
- **`gse295214_h1`:** `"{gsm}:{raw_barcode}"` (`gsm` uppercase `GSM\d+`, `raw_barcode` from `matrix/barcodes`).

Exported identically in: (a) `round16/results/control_barcodes_r16.json`, (b) `round16/results/summary_r16.json` (`control_source_barcodes`, `embeddings_manifest_sha256`, and each line's `ctrl_audit`), (c) `observed.npz` (`matched_ctrl_barcodes`, `resample_ctrl_barcodes`) and each `preds_<arm>.npz` (`control_barcodes`), and (d) `round16/results/embeddings/x_state_<safe_key>.npz` (`control_barcodes`) + `embeddings_manifest_sha256.json` (`§11`).

---

## §8. Within-Line Disjoint Resample Reference Arm (Secondary Descriptive Analysis)

*Motivation:* Master Controller's synthetic verification (`controller/checks/r16_check.py`) demonstrated that when a model exhibits positive target-line signal yet swapping control cells merely injects small isotropic noise (e.g., $\sigma = 0.05$) unrelated to cell-line identity, large $N$ causes paired perturbation bootstrapping to detect a tiny systematic correlation attenuation ($\approx +0.0015$, `95% CI [0.001, 0.002]`) and classify the line as `LINE_MATCHED_SPECIFIC_SUB_DELTA`. Because any stochastic input perturbation slightly attenuates correlation, interpreting small paired losses requires an explicit **within-line resampling null reference**.

1. **Exact Disjoint Sampling Rule (`model_adapters_r16.py::sample_disjoint_resample_control_sets_r16`):**
   - Primary control sampling (`PREREG_I §2.2`, `seed = 20261601`) is executed **completely unchanged** (zero impact on primary `H-I1` or `H-I2`).
   - From the remaining QC-passing NTC cells of target line $t$ **strictly disjoint from the primary sample** (preserving ascending original cell order; let $N_{\mathrm{rem}} = N_{\mathrm{NTC}} - n_{\mathrm{primary}}$):
     - **When $N_{\mathrm{NTC}} \ge 1,024$** ($N_{\mathrm{rem}} \ge 512$, satisfied by all local primary lines): apply the exact `PREREG_I §2.2` batch-stratified sampler with **`seed = 20261602`** (`target_n = 512`, `cell_set_len = 32`, `16 × 32`) on the remaining pool, guaranteeing `0` cell overlap with the primary `512`-cell arm.
     - **When $640 \le N_{\mathrm{NTC}} < 1,024$** ($128 \le N_{\mathrm{rem}} < 512$): use all $N_{\mathrm{rem}}$ disjoint remaining NTC cells without replacement, shuffled once with **`seed = 20261602`** and truncated to `(N_rem // 32) * 32` ($\ge 128$, `0` cell overlap with the locked `512`-cell primary arm).
     - **When $N_{\mathrm{rem}} < 128$** (fewer than `128` cells remain disjoint from the locked `PREREG_I §2.2` primary sample, i.e., $N_{\mathrm{NTC}} < 640$): omit the within-line resample arm (`eligible = False`) and log transparently.
   - Export the second resample arm's ordered barcode list (`control_barcodes_r16.json`, `observed.npz`, `preds_STATE_SE_fixed_index_zero_resample.npz`) and verify `overlap_with_primary == 0`.
2. **Identical Pipeline & Joint Paired Bootstrap (`metrics_r16.py::compute_paired_with_within_line_resample_r16`):**
   - Run the exact same `SE-600M -> ST-SE` (`fixed_index_zero`, with `SE-600M` per-pool seed `= 20261602` per `§11.4`) pipeline (and identically for secondary `STATE_HVG`) on the resample arm to obtain $\hat{\Delta}^{\mathrm{resample}}$.
   - Using the **exact same `n_boot = 4,000` joint perturbation bootstrap draws (`seed = 20261601`)** and noise ceiling as primary `H-I1`:
     - Compute within-line resampling loss $\Delta^{\mathrm{pair}}(\mathrm{resample}) = r(\mathrm{matched}) - r(\mathrm{resample})$ with two-sided `97.5%` and `95.0%` CIs;
     - For every mismatched source $s' \in \mathcal{S}_{\mathrm{mis}}^{\mathrm{seen}}(t) \cup \mathcal{S}_{\mathrm{mis}}^{\mathrm{unseen}}(t)$, compute the paired excess loss over within-line resampling:
       $$\Delta^{\mathrm{excess}}(s') = \Delta^{\mathrm{pair}}(s') - \Delta^{\mathrm{pair}}(\mathrm{resample}) = \big[r(\mathrm{matched}) - r(\mathrm{mis}(s'))\big] - \big[r(\mathrm{matched}) - r(\mathrm{resample})\big]$$
       together with its exact paired `97.5%` CI (`[1.25%, 98.75%]`) and `95.0%` CI (`[2.5%, 97.5%]`).
3. **Locked Descriptive Interpretation Rule (`adjudicate_r16.py::summarize_within_line_resample_r16`, Does NOT Alter Primary `H-I1`/`H-I2`):**
   - For each mismatched source $s'$, bootstrap validity requires `valid_frac >= 0.95` and finite `97.5%` CI; mismatched loss is declared to exceed within-line resampling loss (`exceeds_resample_loss = True`) iff valid AND `excess_ci975_lo > 0.0`. Any pair failing bootstrap validity is reported as invalid/inconclusive (`INCONCLUSIVE_VS_WITHIN_LINE_RESAMPLE`) rather than claiming equivalence to resampling.
   - **Mandatory narrative rule:**
     - When all pairs are valid and mismatched loss does **not** exceed within-line resampling loss: MUST state **「换成其他细胞系的对照，影响不超过重抽本细胞系的对照」**;
     - **ONLY** when all pairs are valid and mismatched loss exceeds within-line resampling loss (`excess_ci975_lo > 0.0` across all sources, alongside positive matched signal `matched_positive_ci975 == True`): may state **「模型用到了目标线特有的对照信息」**;
     - Mixed outcomes across sources/strata (or all pairs exceeding resample loss when `matched_positive_ci975 == False`) MUST be reported source-by-source without making a line-wide target-specificity claim.
   - This analysis is strictly secondary and descriptive; it does NOT change primary `H-I1` or `H-I2` verdicts.

---

## §9. Locked Interpretation of `LINE_MATCHED_SPECIFIC_SUB_DELTA` (`PREREG_I §3.1`)

Whenever any line receives `LINE_MATCHED_SPECIFIC_SUB_DELTA` (or any stratum receives `STRATUM_MATCHED_SPECIFIC_SUB_DELTA`):
1. **Mandatory Wording:** MUST be stated strictly as **「能检测到，但小于事先规定的最小有意义效应 δ=0.04」**.
2. **Mandatory Joint Reporting:** MUST be reported **together with §8's within-line resample contrast** ($\Delta^{\mathrm{pair}}(\mathrm{resample})$ and per-source $\Delta^{\mathrm{pair}}(s') - \Delta^{\mathrm{pair}}(\mathrm{resample})$).
3. **Prohibited Standalone Claim:** `LINE_MATCHED_SPECIFIC_SUB_DELTA` MUST **NEVER** be reported standalone as 「模型用到了对照」.

---

## §10. Exact-Identical Predicted-Delta Boundary Case (`EXACT_INVARIANCE`) (`PREREG_I §3.1`)

1. **Boundary Problem Closed:** Under `PREREG_I v5 §3.1`'s conservative zero-tie rule (`ci975_lo == 0.0 or ci975_hi == 0.0 -> INCONCLUSIVE`), if a model is algebraically or numerically invariant to swapping `X_ctrl` so that matched and mismatched predicted delta arrays are identical ($\max_{p, g} |\hat{\Delta}^{\mathrm{matched}}_{p, g} - \hat{\Delta}^{\mathrm{mis}(s')}_{p, g}| \le 10^{-12}$), every bootstrap paired difference equals `0.0`, producing `97.5% CI = [0.0, 0.0]` and spuriously classifying exact invariance as `INCONCLUSIVE`.
2. **Locked Pair-Level Rule (`adjudicate_r16.py::classify_pair_tost_r16`, `metrics_r16.py::compute_paired_model_vs_all_mismatched_r16`):**
   - Let $d_{\max} = \max_{p, g} |\hat{\Delta}^{\mathrm{matched}}_{p, g} - \hat{\Delta}^{\mathrm{mis}(s')}_{p, g}|$ over the evaluated `(P, G)` predicted delta matrices (requiring all entries finite).
   - Revised strict per-pair priority order:
     1. **`INCONCLUSIVE`** if `valid_frac < 0.95` or `97.5%` CI is non-finite (`within_delta = False`);
     2. **`EXACT_INVARIANCE`** if both predicted delta matrices are finite AND **$d_{\max} \le 10^{-12}$**;
     3. **`INCONCLUSIVE`** (conservative zero-tie rule preserved for non-identical arrays) if `ci975_lo == 0.0` OR `ci975_hi == 0.0`;
     4. **`MATCHED_BETTER`** if `ci975_lo > 0.0`;
     5. **`MISMATCHED_BETTER`** if `ci975_hi < 0.0`;
     6. **`EQUIVALENT`** if `within_delta == True`;
     7. **`INCONCLUSIVE`** for all remaining cases (including exact `±0.04` open-boundary ties).
3. **Locked Stratum- and Line-Level Counting (`adjudicate_r16.py::_summarize_stratum_r16` & `adjudicate_claim_i1_line`):**
   - `EXACT_INVARIANCE` is reported in its own count slot (`counts["EXACT_INVARIANCE"]`) for transparency, and **counts identically as `EQUIVALENT`** in all stratum-level and line-level `H-I1` rules (`n_equiv_effective = counts["EQUIVALENT"] + counts["EXACT_INVARIANCE"]`).
   - Matched-arm validity (`matched_arm_valid`: `valid_frac >= 0.95` AND finite `rho_ab_orth > 0.0` AND finite matched `97.5%` CI) still takes precedence over equivalence, and `LINE_EQUIVALENT_NONZERO_SIGNAL` vs. `LINE_EQUIVALENT_NULL_SIGNAL` remains strictly split by `matched_positive_ci975`.

---

## §11. SE-600M Official Collator, Per-Pool Seeding, & `X_state` Embedding Export (`PREREG_I §2.1.4` & `§5.5`)

*Motivation:* Inside `arc-state 0.11.1`'s official `VCIDatasetSentenceCollator` (`state/emb/data/loader.py:360, 529–576`), sentence construction uses PyTorch's global RNG (`torch.randperm` to break ties among equal-/zero-count genes whenever a cell expresses fewer than `pad_length - 1 = 2,047` mapped genes, plus `torch.randint`/`torch.randperm` for batch `shared_genes` and task tokens). To guarantee exact reproducibility across processes and cross-party parity between Main Line (`round16/src/`) and `A-T3` (`round16/audit/`), BOTH independent implementations MUST follow all six rules below verbatim:

1. **Official `arc-state 0.11.1` Collation (`FilteredGenesCounts` + `VCIDatasetSentenceCollator`):**
   - All `SE-600M` control-pool encoding MUST call `arc-state 0.11.1`'s official `state.emb.data.loader.FilteredGenesCounts` and `state.emb.data.loader.VCIDatasetSentenceCollator(..., is_train=False, precision=torch.float32)` with `collator.training = False` in pure `float32`.
   - Because `FilteredGenesCounts` and `VCIDatasetSentenceCollator` are part of the upstream frozen model library rather than Round 16 analysis code, both parties invoking this identical library step preserves complete two-party implementation independence.
2. **Per-Pool RNG Seeding Immediately Before Encoding (Once Per Pool):**
   - Immediately before the `DataLoader` encoding loop (`for batch in dl:`) begins for each control pool — strictly **after** `StateEmbeddingModel` instantiation and `load_state_dict` so model initialization never perturbs global PyTorch RNG — call the following functions in exact order **once per pool**:
     1. `random.seed(seed)`
     2. `np.random.seed(seed)`
     3. `torch.manual_seed(seed)`
     4. If executing on Apple MPS (`device == "mps"`): `torch.mps.manual_seed(seed)`
   - **Primary control pools (matched target and every mismatched control source):** **`seed = 20261601`** (`DEFAULT_R16_SEED`).
3. **Locked Input Ordering & DataLoader Parameters:**
   - **Cell order:** Exact sampled order produced by `PREREG_I §2.2` / `§8` (ascending sampled position index `0..n_ctrl-1`, e.g., `0..511` in the `16 × 32` layout).
   - **Gene column order & mapping:** Exact raw file column order (`adata.var["gene_name"] = list(raw_gene_names)`, `adata.var_names = list(raw_gene_names)`, `adata.var_names_make_unique()` retaining first occurrences on duplicate raw symbols without duplicate global IDs), mapped strictly by exact gene-symbol lookup into `SE-600M/protein_embeddings.pt` (`19,790` genes; zero HGNC synonym remapping per `PREREG_I §2.1.4`).
   - **DataLoader settings:** `batch_size = 16`, `shuffle = False`, `num_workers = 0`.
4. **Within-Line Disjoint Resample Arm Seeding (`§8`):**
   - Every eligible within-line disjoint resample control pool (`§8`) is encoded as its own pool using **`seed = 20261602`** (`DEFAULT_RESAMPLE_SEED`) via the exact same per-pool seeding sequence and official collator.
5. **Mandatory `X_state` Embedding Export, SHA-256 Manifest, & Staged Three-Party Comparison:**
   - Export every encoded control pool's `2,058`-dim `float32` `X_state` array (`[emb_2048 || ds_emb_10]`) alongside `control_barcodes`, `seed`, `dataset_name`, and `se_audit_json` to `round16/results/embeddings/x_state_<safe_key>.npz` (primary key `<line>`, resample key `<line>_resample`), plus `round16/results/embeddings/embeddings_manifest_sha256.json` (`{filename: sha256}`, mirrored in `summary_r16.json["embeddings_manifest_sha256"]`).
   - **Staged Verification Order:** Master Controller compares Main Line vs. `A-T3` `X_state` embeddings **first** under locked tolerance **$\max_{i, d} |X_{\mathrm{state}}^{\mathrm{Main}}[i, d] - X_{\mathrm{state}}^{\mathrm{AT3}}[i, d]| \le 10^{-5}$** across all primary and resample control pools before comparing downstream predicted delta arrays (`<= 1e-4`, `§7`).
6. **Cross-Process Embedding Determinism Coverage:**
   - Both Main Line (`test_unit_r16.py::test_18_cross_process_pythonhashseed_determinism`) and `A-T3` reproducibility tests MUST cover `SE-600M` control encoding across two independent fresh processes (`PYTHONHASHSEED=1` vs. `99991`) on synthetic integer UMI counts containing zero-count genes (exercising `loader.py:529–533` `torch.randperm` tie-breaking), asserting bit-for-bit identical (`max_abs == 0.0`) `X_state` embeddings on CPU for both `seed = 20261601` and `seed = 20261602` (and `max_abs <= 1e-6` on MPS when available).
   - **Code Correspondence:** `round16/src/model_adapters_r16.py::encode_controls_se600m_r16(..., batch_size=16, seed=...)`, `round16/src/run_eval_r16.py::save_control_embeddings_r16`, and `round16/src/test_unit_r16.py::test_11`, `test_18`.
