
suppressPackageStartupMessages(library(dplyr))
d_exome <- read.delim("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/pipeline_execution_runs/vg_exome_24.tsv", stringsAsFactors=FALSE)
d_clinvar <- read.delim("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/pipeline_execution_runs/vg_clinvar_24.tsv", stringsAsFactors=FALSE)
d_exome$Chromosome <- as.character(d_exome$Chromosome)
d_clinvar$Chromosome <- as.character(d_clinvar$Chromosome)

# Exact lines 88-89 of a03_acmg_gene_clinvar.r: keeping only coding variants
d_exome <- d_exome %>% filter(funct_conseq != "intronic" & funct_conseq != "downstream" & funct_conseq != "upstream" & funct_conseq != "3_prime_UTR" & funct_conseq != "5_prime_UTR")
# Exact line 119: merge by Chromosome, Position, REF, ALT
d_merge <- left_join(d_exome, d_clinvar, by=c("Chromosome","Position","REF","ALT"))
# Exact lines 142-155: clinical_significance_categories
d_merge <- d_merge %>% mutate(clinical_significance_categories = "") %>%
  mutate(clinical_significance_categories = ifelse(grepl("likely pathogenic", ClinicalSignificance, ignore.case=T), "likely pathogenic", clinical_significance_categories)) %>%
  mutate(clinical_significance_categories = ifelse(grepl("likely benign", ClinicalSignificance, ignore.case=T), "likely benign", clinical_significance_categories)) %>%
  mutate(clinical_significance_categories = ifelse(grepl("conflicting", ClinicalSignificance, ignore.case=T), "conflicting", clinical_significance_categories)) %>%
  mutate(clinical_significance_categories = ifelse(grepl("uncertain", ClinicalSignificance, ignore.case=T), "uncertain", clinical_significance_categories)) %>%
  mutate(clinical_significance_categories = ifelse(ClinicalSignificance == "Pathogenic", "pathogenic", clinical_significance_categories)) %>%
  mutate(clinical_significance_categories = ifelse(grepl("Pathogenic[|/,]", ClinicalSignificance, ignore.case=T), "pathogenic", clinical_significance_categories))
# Exact lines 168-178: star_status
d_merge <- d_merge %>% mutate(star_status = 0L) %>%
  mutate(star_status = ifelse(grepl("single submitter", ReviewStatus, ignore.case=T), 1L, star_status)) %>%
  mutate(star_status = ifelse(grepl("multiple submitters", ReviewStatus, ignore.case=T), 2L, star_status)) %>%
  mutate(star_status = ifelse(grepl("reviewed by expert panel", ReviewStatus, ignore.case=T), 3L, star_status)) %>%
  mutate(star_status = ifelse(grepl("practice guideline", ReviewStatus, ignore.case=T), 4L, star_status))
# Exact lines 189-191 of a03_acmg_gene_clinvar.r + a08_summarise.r line 134 (HFE homozygotes only!)
d_filt <- d_merge %>% filter(grepl("pathogenic", clinical_significance_categories)) %>% filter(star_status >= 1)
# HFE special handling in a08_summarise.r (lines 134-135: "only interested in homozygotes" for HFE rs1800562)
d_filt <- d_filt %>% filter(!(hgnc_name == "HFE" & GT != "1/1"))
write.table(d_filt, "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/pipeline_execution_runs/vg_output.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("viking-genes clinvar_pipeline returned %d / 24 variants\n", nrow(d_filt)))
