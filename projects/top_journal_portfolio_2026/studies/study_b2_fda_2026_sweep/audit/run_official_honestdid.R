#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(HonestDiD)
  library(jsonlite)
  library(doParallel)
})

n_cores <- min(8L, max(1L, parallel::detectCores() - 2L))
registerDoParallel(cores = n_cores)
cat(sprintf("Official HonestDiD v%s initialized with %d parallel cores\n",
            as.character(packageVersion("HonestDiD")), n_cores))
flush.console()

in_path  <- "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep/audit/b2_t06_honestdid_input_payload.json"
out_path <- "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep/audit/b2_t06_official_r_honestdid_results.json"

payload <- fromJSON(in_path, simplifyVector = FALSE)

# 1. Official Vignette Benchmark Check (BCdata_EventStudy)
data("BCdata_EventStudy")
vig_sd <- createSensitivityResults(
  betahat = BCdata_EventStudy$betahat,
  sigma = BCdata_EventStudy$sigma,
  numPrePeriods = 4, numPostPeriods = 4,
  Mvec = c(0, 0.01, 0.02, 0.05)
)
vig_rm <- createSensitivityResults_relativeMagnitudes(
  betahat = BCdata_EventStudy$betahat,
  sigma = BCdata_EventStudy$sigma,
  numPrePeriods = 4, numPostPeriods = 4,
  Mbarvec = c(0.5, 1.0)
)

results <- list(
  executed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  r_version = R.version.string,
  honestdid_version = as.character(packageVersion("HonestDiD")),
  parallel_cores = n_cores,
  vignette_benchmark_BCdata = list(
    delta_SD_FLCI = vig_sd,
    delta_RM_CLF = vig_rm
  ),
  specs = list()
)

df_to_records <- function(df) {
  df$lb_irr <- round(exp(df$lb), 4)
  df$ub_irr <- round(exp(df$ub), 4)
  df$excludes_zero <- (df$lb > 0) | (df$ub < 0)
  df
}

# Run quarterly specs FIRST (exact default gridPoints=1000), then key monthly specs (gridPoints=100, parallel=TRUE)!
spec_order <- c(
  "quarterly_ot_od_ind_only_ref2026Q1_Q2_impact",
  "quarterly_ot_od_ind_only_ref2026Q1_Q2Q3_avg",
  "quarterly_ot_od_ddd_ref2026Q1_Q2_impact",
  "quarterly_ot_od_ddd_ref2026Q1_Q2Q3_avg",
  "quarterly_within_act_ref2026Q1_Q2_impact",
  "quarterly_within_act_ref2026Q1_Q2Q3_avg",
  "monthly_ot_od_ind_only_ref2602_apr_impact",
  "monthly_ot_od_ind_only_ref2602_apr_sep_avg",
  "monthly_ot_od_ddd_ref2602_apr_impact",
  "monthly_ot_od_ddd_ref2602_apr_sep_avg",
  "monthly_within_act_ref2602_apr_sep_avg",
  "monthly_ot_od_ind_only_ref2603_apr_sep_avg",
  "monthly_ot_od_ddd_ref2603_apr_sep_avg"
)

for (sname in spec_order) {
  t_sp <- proc.time()[3]
  sp <- payload$specs[[sname]]
  betahat <- unlist(sp$betahat)
  sigma <- do.call(rbind, lapply(sp$sigma, unlist))
  sigma <- 0.5 * (sigma + t(sigma))
  n_pre <- as.integer(sp$numPrePeriods)
  n_post <- as.integer(sp$numPostPeriods)
  l_vec <- matrix(unlist(sp$l_vec), ncol = 1)
  is_q <- grepl("^quarterly_", sname)
  gp <- if (is_q) 1000L else 100L
  m_vec <- if (is_q) c(0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30) else c(0, 0.02, 0.05, 0.10, 0.20)
  mbar_vec <- if (is_q) c(0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0) else c(0, 0.25, 0.5, 1.0)

  orig <- constructOriginalCS(
    betahat = betahat, sigma = sigma,
    numPrePeriods = n_pre, numPostPeriods = n_post,
    l_vec = l_vec
  )

  # Delta^SD(M): Official FLCI (Fixed-Length Confidence Interval)
  sd_flci <- tryCatch({
    res <- createSensitivityResults(
      betahat = betahat, sigma = sigma,
      numPrePeriods = n_pre, numPostPeriods = n_post,
      l_vec = l_vec, method = "FLCI",
      Mvec = m_vec
    )
    df_to_records(as.data.frame(res))
  }, error = function(e) list(error = as.character(e)))

  # Delta^RM(Mbar): Official Conditional / Hybrid (C-LF)
  rm_clf <- tryCatch({
    res <- createSensitivityResults_relativeMagnitudes(
      betahat = betahat, sigma = sigma,
      numPrePeriods = n_pre, numPostPeriods = n_post,
      l_vec = l_vec, method = "C-LF",
      Mbarvec = mbar_vec,
      gridPoints = gp,
      parallel = TRUE
    )
    df_to_records(as.data.frame(res))
  }, error = function(e) list(error = as.character(e)))

  theta_hat <- as.numeric(t(l_vec) %*% betahat[(n_pre + 1):(n_pre + n_post)])
  se_theta  <- sqrt(as.numeric(t(l_vec) %*% sigma[(n_pre + 1):(n_pre + n_post), (n_pre + 1):(n_pre + n_post)] %*% l_vec))

  results$specs[[sname]] <- list(
    description = sp$description,
    ref_period = sp$ref_period,
    periods = unlist(sp$periods),
    numPrePeriods = n_pre,
    numPostPeriods = n_post,
    l_vec = as.numeric(l_vec),
    gridPoints_rm = gp,
    theta_hat_log_irr = round(theta_hat, 5),
    theta_hat_irr = round(exp(theta_hat), 4),
    se_theta = round(se_theta, 5),
    original_cs = df_to_records(as.data.frame(orig)),
    delta_SD_FLCI = sd_flci,
    delta_RM_CLF = rm_clf
  )
  write_json(results, out_path, pretty = TRUE, auto_unbox = TRUE, digits = 6)
  dt_s <- proc.time()[3] - t_sp
  cat(sprintf("[%.1fs] %s | theta_irr=%.3f | orig=[%.3f, %.3f] | SD_M0.05=[%.3f,%.3f] | RM_0.5=[%.3f,%.3f]\n",
              dt_s, sname, exp(theta_hat), exp(orig$lb), exp(orig$ub),
              if (is.data.frame(sd_flci)) sd_flci$lb_irr[sd_flci$M == 0.05] else NA,
              if (is.data.frame(sd_flci)) sd_flci$ub_irr[sd_flci$M == 0.05] else NA,
              if (is.data.frame(rm_clf)) rm_clf$lb_irr[rm_clf$Mbar == 0.5] else NA,
              if (is.data.frame(rm_clf)) rm_clf$ub_irr[rm_clf$Mbar == 0.5] else NA))
  flush.console()
}

stopImplicitCluster()
cat("ALL OFFICIAL HONESTDID SPECS COMPLETED AND SAVED!\n")
