# Portfolio 1 Freeze Decision & Deliverable Manifest

**Date:** November 25, 2026  
**Status:** FROZEN AT PUBLIC-DATA CEILING (Ready for archive / branch commit / PR)

## Frozen Manuscripts & Verified Outputs

1. **Study A1-X v2 (Primary Upgraded Deliverable — Target Tier: *Nature Communications* / *Genome Medicine* / *AJHG* / *Genetics in Medicine*):**
   - Manuscript PDF: `studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/v2_upgrade/manuscript_a1x_v2.pdf`
   - Key empirical advances over v1 (directly closing all 6 external reviewer critiques using 100% open public data):
     - **Unbiased complete-universe benchmark:** All `7,427` gnomAD v4.1-observed (`4,799` high-callability `2–4★`) AD ACMG SF v3.3 ClinVar P/LP alleles across 73 AD genes, conservatively curated via live ClinVar VCV XML (explicitly preserving AD single-exon `PLN p.Leu39Ter` and 3★ `TNNT2 p.Trp297Ter`).
     - **Tier A clinical/biobank standards added:** `Exomiser` (`98f4e0b`, Monarch / Genomics England) and `AutoGVP` (`054e2b7`, CHOP / NIH Kids First) added alongside `CPSR v2.4.0`, `SFtool`, `CPHI-DRAGEN-anno`, `viking-genes`, `Talos v12.3.1`, and Tier B tools (`OwnGenome`, `DNA`); confirmed that NCBI ClinVar omits `ModeOfInheritance` across **all three** flat releases (`clinvar.vcf.gz`, `variant_summary.txt.gz`, `submission_summary.txt.gz`).
     - **Real individual-level human genomes (`1000 Genomes Project`, `N = 2,504` unrelated / `3,202` total `30×` GRCh38 WGS):** Zero multi-allele overlaps across all `3,202` real genomes (`dedup ratio = 1.000`); `3/6` (`50.0%`, `6.13/1,000`) high-confidence South Asian AD findings are strictly recessive 3★ `RYR1 p.Arg3550Trp`; direct execution of compiled Rust `OwnGenome`, `CPSR v2.4.0`, and `SFtool` exposed a fourth VCF failure mode — overlapping-gene `GENEINFO` bleed-through misassigning AR `NTHL1 p.Gln82Ter` to AD `TSC2` in `2/503` Europeans (`HG00182`, `HG01613`).
     - **Direct human biomarker phenotype proof (`UK Biobank`, Backman et al. *Nature* 2021, `N = 411,351` exomes):** Across all `13` matched `2–4★` `APOB` pLoFs misreturned as Familial Hypercholesterolemia by 8/8 pipelines, Direct LDL-C is lowered (`IVW β = -1.727 ± 0.058 SD, p = 3.98e-194`), directly opposite to true FH `APOB` missense (`+0.976 SD`, `Δβ = -2.703 SD, Z = -34.67`) and `LDLR` (`+0.311 SD`, `n = 84`, `p = 6.34e-37`), replicating across Total Cholesterol (`N = 412,075`) and ApoB (`N = 410,113`) alongside Nightingale NMR (`N = 434,646`) and clinical odds ratios.
2. **Study A1-X v1 (Archived Baseline Pre-Upgrade Version):**
   - `studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/manuscript_a1x/manuscript_a1x.pdf` (8 pp., preserved bit-identical).
3. **Study A1 (Phase 1 Accepted Manuscript — Preserved Bit-Identical):**
   - `studies/study_a1_clinvar_gnomad_reclass/manuscript_a1/manuscript_a1.pdf` (18 pp., SHA-256 `73e625710e011f09fb33c8a124ef9e7ed3d8e1ee3697615b436ec03610e7bceb`).
4. **Study B-2 (Phase 1 Accepted Manuscript + Phase 2 Adversarial Confounding Audit — Preserved Bit-Identical):**
   - `studies/study_b2_fda_2026_sweep/manuscript_b2/manuscript_b2.pdf` (16 pp., SHA-256 `8ab4219b0171d2e4522c8d09904156cd3dc1d310d74c4482f73b59255ee70b47`).
   - Phase 2 Attack 2D/2E verdict (`ATTACK2_B2_CONFOUNDING_VERDICT.md`): transparently documents that while 2026 Q1 broad-sponsor diffusion differs sharply from 2024 single-sponsor bulk uploads (`ImmunityBio`), March-2026 ACT-vs.-Non-ACT separation remains constrained by pre-trends and unmeasured status-differential shocks; honest scope is descriptive quasi-experimental audit rather than causal identification.

## Honest Ceiling Rationale for Freezing Portfolio 1
Further elevation of A1-X beyond *Nature Communications* / *Genome Medicine* / *AJHG* into *Nature* / *Nature Medicine* main journal would require controlled-access clinical diagnostic lab final sign-off charts (to measure real-world manual curation interception rates and post-disclosure patient harm) or prospective multi-site clinical implementation trials, which are outside open public data. With all open-data proofs exhausted cleanly and rigorously, Portfolio 1 is frozen in `FREEZE_MANIFEST.json` and autonomous discovery advances immediately to **Project #3 (`C2-ZOMBIE`)**.
