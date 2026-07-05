# parameter_sweep.R
# ============================================================
# Micro parameter sweep for the grant-funding model
# ============================================================
# Usage:
#   Rscript parameter_sweep.R
# ============================================================

cat("=== Parameter Sweep: MICRO ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ---- Source model code (lines 1-818 of app.R, stripped) ----

source_model <- function() {
  app_lines <- readLines("app.R", warn = FALSE)
  app_lines <- app_lines[1:min(818, length(app_lines))]
  strip_patterns <- c(
    "^options\\(",
    "^suppressPackageStartupMessages",
    "^\\s*library\\(",
    "^\\}\\)$"
  )
  keep <- rep(TRUE, length(app_lines))
  for (pat in strip_patterns) {
    keep <- keep & !grepl(pat, app_lines)
  }
  code <- paste(app_lines[keep], collapse = "\n")
  eval(parse(text = code), envir = globalenv())
}

source_model()
cat("Model sourced successfully.\n")

# ---- Output directory ----
out_dir <- file.path(getwd(), "sweep_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ============================================================
# Fixed defaults
# ============================================================

DEFAULTS <- list(
  k_min           = 1.0,
  k_shape         = 2.0,
  r_min           = 1.0,
  r_shape         = 2.0,
  rho_kr          = 0.0,
  gamma           = 1.0,
  epsilon         = 0.1,
  delta           = 1.0,
  tau_r           = 1.0,
  tau_k           = 1.0,
  use_resource_signal = TRUE,
  n_pre_rounds    = 0,
  x_seed          = 0.5,
  n               = 50,
  budget_fraction = 0.10
)

# ============================================================
# Build conditions — Micro sweep only
# ============================================================

build_conditions <- function() {
  conds <- list()

  # ---- Default condition (special case) ----
  # M=500, 30 seeds, strategies 1:9
  conds[[length(conds) + 1]] <- list(
    sweep_type = "micro", block = "default",
    n = 50, M = 500, seeds = 30, strategies = 1:9,
    is_default = TRUE
  )

  # ---- Parameter sweeps (M=200, 7 seeds, per-block strategies) ----
  # Only non-default values (default=50 for n, 2.0 for alpha, etc.)

  # n (pool size): 10, 50*, 100  (* = default, skip)
  for (nv in c(10, 100)) {
    conds[[length(conds) + 1]] <- list(
      sweep_type = "micro", block = "n",
      n = nv, M = 200, seeds = 7, strategies = c(1, 2, 4, 5, 6)
    )
  }

  # alpha (Pareto shape): 1.5, 2.0*, 3.0
  for (a in c(1.5, 3.0)) {
    conds[[length(conds) + 1]] <- list(
      sweep_type = "micro", block = "alpha",
      n = 50, M = 200, seeds = 7, strategies = c(1, 2, 4, 5, 6),
      k_shape = a, r_shape = a
    )
  }

  # r_min/k_min: 0.5, 1.0*, 2.0
  for (ratio in c(0.5, 2.0)) {
    conds[[length(conds) + 1]] <- list(
      sweep_type = "micro", block = "ratio",
      n = 50, M = 200, seeds = 7, strategies = c(1, 2, 4, 5, 6),
      r_min = ratio
    )
  }

  # budget_fraction: 0.01, 0.10*, 1.0
  for (bf in c(0.01, 1.0)) {
    conds[[length(conds) + 1]] <- list(
      sweep_type = "micro", block = "budget",
      n = 50, M = 200, seeds = 7, strategies = c(1, 2, 4, 5),
      budget_fraction = bf
    )
  }

  # tau_k: 0.1, 1.0*, 10.0
  for (tk in c(0.1, 10.0)) {
    conds[[length(conds) + 1]] <- list(
      sweep_type = "micro", block = "tau_k",
      n = 50, M = 200, seeds = 7, strategies = c(1, 4, 5),
      tau_k = tk
    )
  }

  # tau_r: 0.1, 1.0*, 10.0
  for (tr in c(0.1, 10.0)) {
    conds[[length(conds) + 1]] <- list(
      sweep_type = "micro", block = "tau_r",
      n = 50, M = 200, seeds = 7, strategies = c(1, 4, 5),
      tau_r = tr
    )
  }

  # x_seed: 0.1, 0.5*, 1.0
  for (xs in c(0.1, 1.0)) {
    conds[[length(conds) + 1]] <- list(
      sweep_type = "micro", block = "x_seed",
      n = 50, M = 200, seeds = 7, strategies = c(1, 2, 4, 6),
      x_seed = xs
    )
  }

  conds
}

# ============================================================
# Run one condition x seed
# ============================================================

run_one <- function(cond, seed_val) {
  # Resolve parameters: condition overrides > defaults
  n_sim   <- cond$n      %||% DEFAULTS$n
  M_sim   <- cond$M      %||% 200
  k_shape <- cond$k_shape %||% DEFAULTS$k_shape
  r_shape <- cond$r_shape %||% DEFAULTS$r_shape
  r_min   <- cond$r_min   %||% DEFAULTS$r_min
  k_min   <- cond$k_min   %||% DEFAULTS$k_min
  tau_k   <- cond$tau_k   %||% DEFAULTS$tau_k
  tau_r   <- cond$tau_r   %||% DEFAULTS$tau_r
  budget_fraction <- cond$budget_fraction %||% DEFAULTS$budget_fraction
  gamma   <- cond$gamma   %||% DEFAULTS$gamma
  epsilon <- cond$epsilon  %||% DEFAULTS$epsilon
  rho_kr  <- cond$rho_kr  %||% DEFAULTS$rho_kr
  n_pre   <- cond$n_pre_rounds %||% DEFAULTS$n_pre_rounds
  x_seed  <- cond$x_seed  %||% DEFAULTS$x_seed
  delta   <- DEFAULTS$delta
  strategies <- cond$strategies

  # Budget conversion
  alpha_budget <- r_shape
  alpha_budget_sub <- alpha_budget
  if (alpha_budget <= 1.0) {
    alpha_budget_sub <- 1.5
  }
  R_median <- r_min * 2^(1 / alpha_budget_sub)
  B_raw    <- budget_fraction * n_sim * R_median

  # Flags

  degenerate <- B_raw < delta
  ess_risk   <- (!is.infinite(tau_k) && tau_k < 0.25) ||
                (!is.infinite(tau_r) && tau_r < 0.25)

  # Map Inf tau to 1e6 for numerical stability
  tau_k_eff <- if (is.infinite(tau_k)) 1e6 else tau_k
  tau_r_eff <- if (is.infinite(tau_r)) 1e6 else tau_r

  # Run simulation
  sim <- run_simulation(
    seed   = seed_val,
    n      = n_sim,
    k_min  = k_min, k_shape = k_shape,
    r_min  = r_min, r_shape = r_shape,
    rho_kr = rho_kr,
    gamma  = gamma, epsilon = epsilon,
    B      = B_raw,
    delta  = delta,
    tau_r  = tau_r_eff, tau_k = tau_k_eff,
    use_resource_signal = DEFAULTS$use_resource_signal,
    n_pre_rounds = n_pre,
    x_seed = x_seed,
    M      = M_sim,
    strategies = strategies,
    verbose = FALSE
  )

  # ---- Pre-funding population summaries ----
  K_init <- sim$K_at_start
  R_init <- sim$R0_at_start
  D_init <- (K_init - R_init) / (K_init + R_init)
  S_init <- D_init^2

  row <- data.frame(
    sweep_type       = cond$sweep_type,
    block            = cond$block,
    seed             = seed_val,
    n                = n_sim,
    alpha            = k_shape,
    k_min            = k_min,
    k_shape          = k_shape,
    r_min            = r_min,
    r_shape          = r_shape,
    r_min_k_min_ratio = r_min / k_min,
    budget_fraction  = budget_fraction,
    B_raw            = B_raw,
    tau_k            = tau_k,
    tau_r            = tau_r,
    gamma            = gamma,
    epsilon          = epsilon,
    rho_kr           = rho_kr,
    n_pre_rounds     = n_pre,
    x_seed           = x_seed,
    M                = M_sim,
    ess_risk         = ess_risk,
    degenerate       = degenerate,
    alpha_budget_sub = alpha_budget_sub,
    mean_K_init      = mean(K_init),
    mean_R_init      = mean(R_init),
    mean_D_init      = mean(D_init),
    mean_S_init      = mean(S_init),
    cor_KR_init      = if (length(K_init) > 1) cor(K_init, R_init) else NA_real_,
    stringsAsFactors = FALSE
  )

  # ---- Per-strategy columns ----
  for (s in 1:9) {
    pexp <- paste0("total_expected_s", s)
    pout <- paste0("total_output_s", s)
    pK2  <- paste0("mean_K2_s", s)
    pR2  <- paste0("mean_R2_s", s)
    pD2  <- paste0("mean_D2_s", s)
    pS2  <- paste0("mean_S2_s", s)
    pcor <- paste0("cor_K2R2_s", s)

    if (s %in% strategies && !is.null(sim$strategies[[s]])) {
      st <- sim$strategies[[s]]
      row[[pexp]] <- st$total_expected
      row[[pout]] <- st$total_output

      K2 <- st$K2
      R2 <- st$R2
      D2 <- (K2 - R2) / (K2 + R2)
      row[[pK2]]  <- mean(K2)
      row[[pR2]]  <- mean(R2)
      row[[pD2]]  <- mean(D2)
      row[[pS2]]  <- mean(D2^2)
      row[[pcor]] <- if (length(K2) > 1) cor(K2, R2) else NA_real_
    } else {
      row[[pexp]] <- NA_real_
      row[[pout]] <- NA_real_
      row[[pK2]]  <- NA_real_
      row[[pR2]]  <- NA_real_
      row[[pD2]]  <- NA_real_
      row[[pS2]]  <- NA_real_
      row[[pcor]] <- NA_real_
    }
  }

  # Return sim object alongside row for default condition researcher data
  list(row = row, sim = sim)
}

# ============================================================
# Run sweep
# ============================================================

run_sweep <- function() {
  conds <- build_conditions()
  n_conds <- length(conds)
  n_total_runs <- sum(vapply(conds, function(c) as.integer(c$seeds), integer(1)))
  cat(sprintf("Total conditions: %d (%d runs)\n", n_conds, n_total_runs))

  n_cores <- max(1, parallel::detectCores() - 1)
  cat(sprintf("Using %d cores\n\n", n_cores))

  # Storage for default condition researcher data
  default_researchers <- list()

  timing_conds <- min(5, n_conds)

  for (ci in seq_len(n_conds)) {
    cond <- conds[[ci]]
    n_seeds <- cond$seeds
    is_default <- isTRUE(cond$is_default)
    cond_id <- sprintf("micro_%05d", ci)
    result_file <- file.path(out_dir, sprintf("%s.rds", cond_id))

    # Skip if already computed (unless default — need researcher data)
    if (file.exists(result_file) && !is_default) {
      if (ci <= timing_conds) {
        cat(sprintf("[%d/%d] %s — skipped (exists)\n", ci, n_conds, cond$block))
      }
      next
    }

    t0 <- proc.time()["elapsed"]

    # Run all seeds for this condition in parallel
    seed_list <- seq_len(n_seeds)
    results_raw <- parallel::mclapply(seed_list, function(s) {
      tryCatch(
        run_one(cond, s),
        error = function(e) {
          warning(sprintf("Cond %d seed %d failed: %s", ci, s, e$message))
          NULL
        }
      )
    }, mc.cores = n_cores)

    # Filter NULLs
    results_raw <- results_raw[!vapply(results_raw, is.null, logical(1))]

    if (length(results_raw) > 0) {
      # Extract rows
      rows <- lapply(results_raw, `[[`, "row")
      df <- do.call(rbind, rows)
      saveRDS(df, result_file)

      # For default condition: extract researcher-level data
      if (is_default) {
        default_researchers <- lapply(results_raw, function(r) {
          sim <- r$sim
          seed_val <- r$row$seed
          n_res <- sim$params$n

          researcher_df <- data.frame(
            seed = seed_val,
            researcher_id = seq_len(n_res),
            K_init = sim$K_at_start,
            R_init = sim$R0_at_start,
            stringsAsFactors = FALSE
          )

          for (s in 1:9) {
            if (!is.null(sim$strategies[[s]])) {
              st <- sim$strategies[[s]]
              researcher_df[[paste0("K2_s", s)]] <- st$K2
              researcher_df[[paste0("R2_s", s)]] <- st$R2
              researcher_df[[paste0("g1_s", s)]] <- st$g1
              researcher_df[[paste0("g2_s", s)]] <- st$g2
            }
          }

          researcher_df
        })

        saveRDS(default_researchers,
                file.path(out_dir, "default_condition_researchers.rds"))
        cat("  -> Saved default_condition_researchers.rds\n")
      }
    }

    elapsed <- proc.time()["elapsed"] - t0

    if (ci <= timing_conds) {
      cat(sprintf("[%d/%d] %s — %.1fs (%d seeds, M=%d, strats=%s)\n",
                  ci, n_conds, cond$block, elapsed, n_seeds,
                  cond$M, paste(cond$strategies, collapse=",")))
    }

    if (ci == timing_conds) {
      cat(sprintf("\n"))
    }

    if (ci > timing_conds && ci %% 5 == 0) {
      cat(sprintf("[%d/%d] completed\n", ci, n_conds))
    }
  }

  # ---- Merge all results ----
  cat("\nMerging results...\n")
  rds_files <- sort(list.files(out_dir, pattern = "^micro_.*\\.rds$", full.names = TRUE))
  all_results <- lapply(rds_files, readRDS)
  df_all <- do.call(rbind, all_results)

  # ---- Derived columns ----
  df_all$gain_grant_myopic <- df_all$total_expected_s5 - df_all$total_expected_s4
  df_all$gain_grant_forward <- df_all$total_expected_s8 - df_all$total_expected_s7
  df_all$gain_seed_myopic <- df_all$total_expected_s6 - df_all$total_expected_s4
  df_all$gain_seed_vs_grant <- df_all$total_expected_s6 - df_all$total_expected_s5

  # Best strategy
  exp_cols <- grep("^total_expected_s", names(df_all), value = TRUE)
  exp_mat <- as.matrix(df_all[, exp_cols])
  df_all$best_strategy_idx <- apply(exp_mat, 1, function(x) {
    valid <- which(!is.na(x))
    if (length(valid) == 0) return(NA_integer_)
    idx <- valid[which.max(x[valid])]
    as.integer(gsub("total_expected_s", "", exp_cols[idx]))
  })
  df_all$best_strategy_expected <- apply(exp_mat, 1, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })

  # ---- Save ----
  saveRDS(df_all, file.path(out_dir, "sweep_micro.rds"))
  write.csv(df_all, file.path(out_dir, "sweep_micro.csv"), row.names = FALSE)

  cat(sprintf("\nDone! %d rows saved to sweep_results/sweep_micro.{rds,csv}\n",
              nrow(df_all)))
  cat("Finished:", format(Sys.time()), "\n")

  invisible(df_all)
}

# ============================================================
# Go
# ============================================================
run_sweep()
