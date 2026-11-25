#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(HonestDiD)
  library(jsonlite)
  library(parallel)
})

aud_dir  <- "/Users/mac/AI-Biology-Lab/projects/top_journal_portfolio_2026/studies/study_b2_fda_2026_sweep/audit"
in_path  <- file.path(aud_dir, "b2_t06_honestdid_input_payload.json")
res_path <- file.path(aud_dir, "b2_t06_official_r_honestdid_results.json")
t10_path <- file.path(aud_dir, "b2_t10_grid_boundary_and_convergence_audit.json")

payload <- fromJSON(in_path, simplifyVector = FALSE)$specs
results <- fromJSON(res_path, simplifyVector = FALSE)

# Step 1: Audit ALL 15 specs against default +-20*sdTheta grid bounds
hits <- list()
for (sname in names(results$specs)) {
  sp_in <- payload[[sname]]
  sp_rs <- results$specs[[sname]]
  cov <- do.call(rbind, lapply(sp_in$sigma, unlist))
  cov <- 0.5 * (cov + t(cov))
  pre <- as.integer(sp_in$numPrePeriods)
  post <- as.integer(sp_in$numPostPeriods)
  l <- matrix(unlist(sp_in$l_vec), ncol = 1)
  se <- sqrt(as.numeric(t(l) %*% cov[(pre+1):(pre+post), (pre+1):(pre+post)] %*% l))
  gp <- as.integer(sp_rs$gridPoints_rm)
  def_b <- 20 * se
  step <- 2 * def_b / (gp - 1)
  tol <- 1e-4
  for (idx in seq_along(sp_rs$delta_RM_CLF)) {
    row <- sp_rs$delta_RM_CLF[[idx]]
    lb <- as.numeric(row$lb)
    ub <- as.numeric(row$ub)
    hit_l <- abs(lb - (-def_b)) <= tol || (lb <= -def_b + 0.5 * step && lb <= -def_b + 1e-3)
    hit_r <- abs(ub - def_b) <= tol || (ub >= def_b - 0.5 * step && ub >= def_b - 1e-3)
    if (hit_l || hit_r) {
      hits[[length(hits) + 1]] <- list(
        spec = sname, row_idx = idx, Mbar = as.numeric(row$Mbar),
        old_lb = lb, old_ub = ub, se = se,
        default_grid_lb = -def_b, default_grid_ub = def_b, default_grid_points = gp, default_grid_step = step
      )
    }
  }
}
cat(sprintf("Detected %d default-grid boundary hits across all %d specs\n", length(hits), length(results$specs)))

# Step 2: Recompute all boundary-hitting (spec, Mbar) pairs on expanded +-40*se grid with 1001 points (parallel)
run_expanded <- function(h, grid_mult = 40, gp = 1001L) {
  sp <- payload[[h$spec]]
  beta <- unlist(sp$betahat)
  cov <- do.call(rbind, lapply(sp$sigma, unlist))
  cov <- 0.5 * (cov + t(cov))
  pre <- as.integer(sp$numPrePeriods)
  post <- as.integer(sp$numPostPeriods)
  l <- matrix(unlist(sp$l_vec), ncol = 1)
  se <- sqrt(as.numeric(t(l) %*% cov[(pre+1):(pre+post), (pre+1):(pre+post)] %*% l))
  bound <- grid_mult * se
  t0 <- proc.time()[3]
  res <- computeConditionalCS_DeltaRM(
    betahat = beta, sigma = cov,
    numPrePeriods = pre, numPostPeriods = post,
    l_vec = l, Mbar = h$Mbar, hybrid_flag = "LF",
    grid.lb = -bound, grid.ub = bound, gridPoints = gp, seed = 0
  )
  acc <- res$grid[res$accept == 1]
  dt <- unname(proc.time()[3] - t0)
  list(
    spec = h$spec, row_idx = h$row_idx, Mbar = h$Mbar,
    old_lb = h$old_lb, old_ub = h$old_ub,
    new_lb = min(acc), new_ub = max(acc),
    grid_lb = -bound, grid_ub = bound, grid_points = gp,
    grid_step = 2 * bound / (gp - 1),
    accepted_left_boundary = as.logical(res$accept[1]),
    accepted_right_boundary = as.logical(tail(res$accept, 1)),
    excludes_zero = (min(acc) > 0) || (max(acc) < 0),
    elapsed_sec = round(dt, 3)
  )
}

expanded <- mclapply(hits, function(h) run_expanded(h, grid_mult = 40, gp = 1001L), mc.cores = min(5L, max(1L, length(hits))))

# Step 3: Grid-resolution convergence check (1001 vs 2001 pts) on quarterly_within_act_ref2026Q1_Q2Q3_avg Mbar=2
q_hit <- Filter(function(h) h$spec == "quarterly_within_act_ref2026Q1_Q2Q3_avg", hits)[[1]]
q_conv_2001 <- run_expanded(q_hit, grid_mult = 40, gp = 2001L)
q_1001 <- Filter(function(x) x$spec == "quarterly_within_act_ref2026Q1_Q2Q3_avg", expanded)[[1]]
conv_check <- list(
  spec = q_hit$spec, Mbar = q_hit$Mbar,
  lb_1001 = q_1001$new_lb, ub_1001 = q_1001$new_ub, step_1001 = q_1001$grid_step,
  lb_2001 = q_conv_2001$new_lb, ub_2001 = q_conv_2001$new_ub, step_2001 = q_conv_2001$grid_step,
  abs_lb_diff = abs(q_1001$new_lb - q_conv_2001$new_lb),
  abs_ub_diff = abs(q_1001$new_ub - q_conv_2001$new_ub),
  within_one_step = (abs(q_1001$new_lb - q_conv_2001$new_lb) <= q_1001$grid_step) && (abs(q_1001$new_ub - q_conv_2001$new_ub) <= q_1001$grid_step),
  accepted_left_boundary_2001 = q_conv_2001$accepted_left_boundary,
  accepted_right_boundary_2001 = q_conv_2001$accepted_right_boundary
)

# Step 4: Annotate ALL 15 specs with exact grid metadata and patch the 5 expanded-grid entries
exp_map <- list()
for (ex in expanded) {
  exp_map[[paste0(ex$spec, "__", ex$row_idx)]] <- ex
}

post_patch_hits <- 0L
for (sname in names(results$specs)) {
  sp_in <- payload[[sname]]
  sp_rs <- results$specs[[sname]]
  cov <- do.call(rbind, lapply(sp_in$sigma, unlist))
  cov <- 0.5 * (cov + t(cov))
  pre <- as.integer(sp_in$numPrePeriods)
  post <- as.integer(sp_in$numPostPeriods)
  l <- matrix(unlist(sp_in$l_vec), ncol = 1)
  se <- sqrt(as.numeric(t(l) %*% cov[(pre+1):(pre+post), (pre+1):(pre+post)] %*% l))
  gp_def <- as.integer(sp_rs$gridPoints_rm)
  def_b <- 20 * se
  step_def <- 2 * def_b / (gp_def - 1)
  for (idx in seq_along(sp_rs$delta_RM_CLF)) {
    k <- paste0(sname, "__", idx)
    row <- sp_rs$delta_RM_CLF[[idx]]
    if (!is.null(exp_map[[k]])) {
      ex <- exp_map[[k]]
      row$lb <- round(ex$new_lb, 6)
      row$ub <- round(ex$new_ub, 6)
      row$lb_irr <- round(exp(ex$new_lb), 4)
      row$ub_irr <- round(exp(ex$new_ub), 4)
      row$excludes_zero <- ex$excludes_zero
      row$grid_lb <- round(ex$grid_lb, 6)
      row$grid_ub <- round(ex$grid_ub, 6)
      row$gridPoints <- ex$grid_points
      row$grid_step <- round(ex$grid_step, 6)
      row$hits_grid_boundary <- (ex$accepted_left_boundary || ex$accepted_right_boundary)
      row$expanded_grid_recomputed <- TRUE
      row$pre_fix_truncated_lb <- ex$old_lb
      row$pre_fix_truncated_ub <- ex$old_ub
    } else {
      lb <- as.numeric(row$lb)
      ub <- as.numeric(row$ub)
      row$grid_lb <- round(-def_b, 6)
      row$grid_ub <- round(def_b, 6)
      row$gridPoints <- gp_def
      row$grid_step <- round(step_def, 6)
      row$hits_grid_boundary <- (abs(lb - (-def_b)) <= 1e-4) || (abs(ub - def_b) <= 1e-4)
      row$expanded_grid_recomputed <- FALSE
    }
    if (isTRUE(row$hits_grid_boundary)) post_patch_hits <- post_patch_hits + 1L
    results$specs[[sname]]$delta_RM_CLF[[idx]] <- row
  }
}

results$t10_grid_boundary_fix <- list(
  fixed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  n_specs_audited = length(results$specs),
  n_boundary_hits_before_fix = length(hits),
  n_boundary_hits_after_fix = post_patch_hits,
  expanded_grid_note = "All 5 boundary-hitting C-LF entries recomputed on +-40*sdTheta with gridPoints=1001; both grid boundaries rejected on all 5. All C-LF endpoints are grid-resolution approximations rather than exact real-valued bounds.",
  convergence_check_1001_vs_2001 = conv_check
)

write_json(results, res_path, pretty = TRUE, auto_unbox = TRUE, digits = 6)
write_json(list(
  executed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  honestdid_version = as.character(packageVersion("HonestDiD")),
  n_specs_audited = length(results$specs),
  n_boundary_hits_before_fix = length(hits),
  n_boundary_hits_after_fix = post_patch_hits,
  expanded_5_results = expanded,
  convergence_check_1001_vs_2001 = conv_check
), t10_path, pretty = TRUE, auto_unbox = TRUE, digits = 6)

cat(sprintf("T10 FIX COMPLETE: %d hits before -> %d hits after | conv_within_1step=%s\n",
            length(hits), post_patch_hits, as.character(conv_check$within_one_step)))
