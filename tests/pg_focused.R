# ============================================================
# tests/pg_focused.R  —  focused PG(T) headline confirmation
# ------------------------------------------------------------
# Locks the §6.10 headline before the full manifest: does the
# forward-vs-myopic gap turn positive as the horizon T grows?
# Grid: T ∈ {1..5} × ε ∈ {0.1,0.3,0.6}, 50 seeds, n=50, M=200,
# all 9 strategies. Reports every key contrast with SE.
# ============================================================
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages(library(parallel))

source("model.R")

T_vals   <- 1:5
eps_vals <- c(0.1, 0.3, 0.6)
seeds    <- 1:50
NCORES   <- max(1, detectCores() - 1)

one_run <- function(seed, Tr, eps) {
  r <- run_simulation_T(seed = seed, T_rounds = Tr, n = 50, epsilon = eps,
                        b = 0.5, tau_k = 1.0, x_seed = 0.25, k_shape = 2.0,
                        M = 200, strategies = 1:9)
  o <- vapply(1:9, function(s) r$strategies[[s]]$total_expected, numeric(1))
  c(PG           = o[8] - o[5],   # forward vs myopic, pubs+grant  (headline)
    fwd_vs_myo_P = o[7] - o[4],   # forward vs myopic, pubs only
    fwd_vs_myo_PS= o[9] - o[6],   # forward vs myopic, pubs+seed
    signal_fwd   = o[8] - o[7],   # grant-signal value (forward)
    seed_fwd     = o[9] - o[7],   # seed value (forward)
    b_idx_S8     = r$strategies[[8]]$b_idx)
}

grid <- expand.grid(T = T_vals, eps = eps_vals)
rows <- list()
for (g in seq_len(nrow(grid))) {
  Tr <- grid$T[g]; eps <- grid$eps[g]
  mats <- mclapply(seeds, function(sd) one_run(sd, Tr, eps), mc.cores = NCORES)
  M <- do.call(rbind, mats)
  mu <- colMeans(M, na.rm = TRUE)
  se <- apply(M, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  rows[[g]] <- data.frame(T = Tr, eps = eps, t(mu),
                          setNames(as.list(se), paste0(names(se), "_se")))
}
res <- do.call(rbind, rows)
saveRDS(res, "sweep_results/pg_focused_summary.rds")

cat("\n=== PG(T) headline: forward vs myopic (pubs+grant), n=50, 50 seeds ===\n")
cat(sprintf("%4s %5s | %18s %14s %14s %10s\n",
            "T", "eps", "PG=S8-S5 (z)", "fwd_myo_P", "signal_fwd", "b_idx_S8"))
for (i in seq_len(nrow(res))) {
  z <- res$PG[i] / res$PG_se[i]
  cat(sprintf("%4d %5.2f | %+8.3f ± %5.3f (%4.1f) %+7.3f±%4.2f %+7.3f±%4.2f %8.3f\n",
              res$T[i], res$eps[i], res$PG[i], res$PG_se[i], z,
              res$fwd_vs_myo_P[i], res$fwd_vs_myo_P_se[i],
              res$signal_fwd[i], res$signal_fwd_se[i], res$b_idx_S8[i]))
}
cat("\nSaved: sweep_results/pg_focused_summary.rds\n")
