# Replication Guide: Dual-MOI & Opposite-Mechanism ACMG SF v3.3 Audit (Study A1-X)

## Prerequisites
- Python 3 (`/opt/anaconda3/bin/python3`) with `pandas`, `pyarrow`, `numpy`, `matplotlib`
- R (`/opt/anaconda3/envs/r_honestdid/bin/Rscript`) with `dplyr`, `stringr`
- Rust toolchain (`cargo`) — used to compile `GeiserX/OwnGenome` (`c324b617b5e95a9abea82d352420dfa62238f094`) binary `og`
- TeX Live (`pdflatex`, `bibtex`, `pdftoppm`) and `pypdf` Python environment for PDF compilation and font-matrix verification

## Pinned External Repositories (`cloned_pipelines_manifest.json`)
All 7 public repositories were shallow-cloned and pinned to exact commits under `cloned_pipelines/`:
1. `GeiserX/OwnGenome` @ `c324b617b5e95a9abea82d352420dfa62238f094` (compiled binary: `cloned_pipelines/OwnGenome/rust/target/release/og`)
2. `sigven/cpsr` @ `v2.4.0` (`7f68f9dd8ca5aeb83c8349a479e53a7bd749c52b`)
3. `babelomics/SFtool` @ `0d7224935e6525f65899cf55ea8469cd6885a5f5`
4. `ccmbioinfo/CPHI-DRAGEN-anno` @ `974860d2b66e502a975b47f8dcd5f9db228806c0`
5. `viking-genes/clinvar_pipeline` @ `3ea7abe9d2ab13315ec5153b47c9084c9a1c586a`
6. `unbalancedparentheses/DNA` @ `f459b22ba9cb2c10f19d3dc2e5a3a242432c342d`
7. `populationgenomics/talos` @ `v12.3.1` (`4aada15adc1f903a8020421c3dce7ad6c2082aea`, diagnostic MOI contrast)

## Execution Order
1. **Build locked 24-variant GRCh38 het VCF & execute all 6 real public pipelines (R, Python, compiled Rust `og`):**
   ```bash
   /opt/anaconda3/bin/python3 build_panel_and_run_pipelines.py
   ```
   Produces `locked_24variant_benchmark_panel.csv`, `test_panel_het_grch38.vcf`, `ar_in_ad_pipeline_verdict_matrix.csv`, `ar_in_ad_real_pipeline_audit.json`, and `pipeline_execution_runs/owngenome_synth_run/summary_group_b.json`.

2. **Live ClinVar XML audit of all 69 3★ *RYR1* Expert-Panel alleles, genome-wide gnomAD v4.1 (`807,162`) per-tier burden, and Tier-4 patch validation:**
   ```bash
   /opt/anaconda3/bin/python3 evaluate_pipeline_burden_and_patch.py
   ```
   Produces `ryr1_all_3star_live_vcep_audit.json`, `ryr1_10_3star_cm_vcep_ar_alleles_gnomad_v41.csv`, `genome_wide_per_pipeline_fp_burden_gnomad_v41.csv`, and `genome_wide_burden_and_patch_audit.json`.

3. **Compile standalone Brief Communication manuscript (`manuscript_a1x/manuscript_a1x.pdf`):**
   ```bash
   cd manuscript_a1x && pdflatex manuscript_a1x.tex && bibtex manuscript_a1x && pdflatex manuscript_a1x.tex && pdflatex manuscript_a1x.tex
   ```
   Verified outputs & Crossref citation audit archived in `manuscript_a1x/manuscript_a1x_audit.json` and `manuscript_a1x/cross_audit_verification.json`.
