
suppressPackageStartupMessages({ library(dplyr); library(stringr) })
# Minimal shim ONLY for assertable::assert_colnames (not installed in r_honestdid); verifies exact colnames
assertable_env <- new.env()
assertable_env$assert_colnames <- function(data, colnames, only_colnames=FALSE, quiet=TRUE) {
  missing <- setdiff(colnames, names(data))
  if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse=",")))
  invisible(TRUE)
}
# Source the UNMODIFIED retrieve_secondary_calls function directly from cloned sigven/cpsr v2.4.0 R/utils.R!
utils_lines <- readLines("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/cloned_pipelines/cpsr/R/utils.R")
# Extract lines 127..215 verbatim (replacing only assertable::assert_colnames call with our identical check)
fn_lines <- utils_lines[127:215]
fn_lines <- gsub("assertable::assert_colnames", "assertable_env$assert_colnames", fn_lines, fixed=TRUE)
eval(parse(text = paste(fn_lines, collapse="\n")))

df <- read.delim("/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/pipeline_execution_runs/cpsr_input_24variants.tsv", stringsAsFactors=FALSE)
df$PRIMARY_TARGET <- as.logical(df$PRIMARY_TARGET)
df$LOSS_OF_FUNCTION <- as.logical(df$LOSS_OF_FUNCTION)

res <- retrieve_secondary_calls(df)
write.table(res, "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_a1_clinvar_gnomad_reclass/falsification_ar_in_ad_pipelines/pipeline_execution_runs/cpsr_output_sf.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("CPSR v2.4.0 retrieve_secondary_calls() returned %d / %d variants\n", nrow(res), nrow(df)))
