
suppressPackageStartupMessages({
  library(HonestDiD)
  library(jsonlite)
  library(doParallel)
})
n_cores <- min(8L, max(1L, parallel::detectCores() - 2L))
registerDoParallel(cores = n_cores)
in_path  <- "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep/audit/b2_t06_honestdid_input_payload.json"
out_path <- "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep/audit/b2_t06_official_r_honestdid_results.json"
payload <- fromJSON(in_path, simplifyVector = FALSE)
results <- fromJSON(out_path, simplifyVector = FALSE)

df_to_records <- function(df) {
  df$lb_irr <- round(exp(df$lb), 4)
  df$ub_irr <- round(exp(df$ub), 4)
  df$excludes_zero <- (df$lb > 0) | (df$ub < 0)
  df
}

missing_specs <- setdiff(names(payload$specs), names(results$specs))
cat("Missing specs to append:", paste(missing_specs, collapse=", "), "\n")
for (sname in missing_specs) {
  t_sp <- proc.time()[3]
  sp <- payload$specs[[sname]]
  betahat <- unlist(sp$betahat)
  sigma <- do.call(rbind, lapply(sp$sigma, unlist))
  sigma <- 0.5 * (sigma + t(sigma))
  n_pre <- as.integer(sp$numPrePeriods)
  n_post <- as.integer(sp$numPostPeriods)
  l_vec <- matrix(unlist(sp$l_vec), ncol = 1)
  orig <- constructOriginalCS(betahat=betahat, sigma=sigma, numPrePeriods=n_pre, numPostPeriods=n_post, l_vec=l_vec)
  sd_flci <- tryCatch({
    df_to_records(as.data.frame(createSensitivityResults(betahat=betahat, sigma=sigma, numPrePeriods=n_pre, numPostPeriods=n_post, l_vec=l_vec, method="FLCI", Mvec=c(0, 0.02, 0.05, 0.10, 0.20))))
  }, error=function(e) list(error=as.character(e)))
  rm_clf <- tryCatch({
    df_to_records(as.data.frame(createSensitivityResults_relativeMagnitudes(betahat=betahat, sigma=sigma, numPrePeriods=n_pre, numPostPeriods=n_post, l_vec=l_vec, method="C-LF", Mbarvec=c(0, 0.25, 0.5, 1.0), gridPoints=100L, parallel=TRUE)))
  }, error=function(e) list(error=as.character(e)))
  theta_hat <- as.numeric(t(l_vec) %*% betahat[(n_pre + 1):(n_pre + n_post)])
  se_theta  <- sqrt(as.numeric(t(l_vec) %*% sigma[(n_pre + 1):(n_pre + n_post), (n_pre + 1):(n_pre + n_post)] %*% l_vec))
  results$specs[[sname]] <- list(
    description = sp$description, ref_period = sp$ref_period, periods = unlist(sp$periods),
    numPrePeriods = n_pre, numPostPeriods = n_post, l_vec = as.numeric(l_vec), gridPoints_rm = 100L,
    theta_hat_log_irr = round(theta_hat, 5), theta_hat_irr = round(exp(theta_hat), 4), se_theta = round(se_theta, 5),
    original_cs = df_to_records(as.data.frame(orig)), delta_SD_FLCI = sd_flci, delta_RM_CLF = rm_clf
  )
  cat(sprintf("[%.1fs] Appended %s | theta_irr=%.3f | orig=[%.3f, %.3f]\n",
              proc.time()[3] - t_sp, sname, exp(theta_hat), exp(orig$lb), exp(orig$ub)))
}
results$n_specs_completed <- length(results$specs)
results$completed_all_15_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
write_json(results, out_path, pretty = TRUE, auto_unbox = TRUE, digits = 6)
stopImplicitCluster()
cat("TOTAL SPECS IN OFFICIAL HONESTDID JSON:", length(results$specs), "\n")
