# Package D1/D2: particle convergence and heavy-tail stability.
# Cells: tau_K in {0.3, 1, 3} x alpha_K in {1.3, 2, 3.5} (fixed mean 2) x M in {200, 1000, 5000}.
# 100 populations per cell (same seeds across M, so comparisons are paired); medians reported.
suppressPackageStartupMessages(library(parallel))
source("model.R")
ess_of <- function(post) 1 / sum(post$w^2)
cell <- function(tk, ks, M, seeds = 1:100) {
  kmin <- 2 * (ks - 1) / ks
  rows <- mclapply(seeds, function(sd) {
    rp <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks, k_min = kmin, M = M,
                           allocator = "smooth", strategies = c(1, 2, 4, 5), detail = TRUE)
    ro <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks, k_min = kmin,
                           allocator = "smooth", strategies = c(5), oracle = TRUE)
    # ESS of the round-1 posteriors with review (S5), rebuilt here with the run's own signals
    set.seed(sd * 1000)
    posts <- build_posteriors_hist(rp$p_cumul, rp$signals$sigma_r, rp$signals$sigma_k,
                                   M, kmin, ks, 1, 2, 1, 1, tk, TRUE, TRUE)
    ess <- sapply(posts, ess_of)
    c(S1 = rp$strategies[[1]]$total_expected, S2 = rp$strategies[[2]]$total_expected,
      S4 = rp$strategies[[4]]$total_expected, S5 = rp$strategies[[5]]$total_expected,
      OR = ro$strategies[[5]]$total_expected, ess_med = median(ess), ess_q10 = unname(quantile(ess, 0.1)),
      meanK = mean(rp$K_at_start))
  }, mc.cores = 2)
  d <- do.call(rbind, rows)
  rs <- (d[, "S5"] - d[, "S4"]) / (d[, "OR"] - d[, "S4"]); tg <- (d[, "OR"] - d[, "S2"]) / (d[, "S2"] - d[, "S1"])
  rg <- (d[, "S5"] - d[, "S4"]) / (d[, "S5"] - d[, "S1"])
  se <- function(x) sd(x) / sqrt(length(x))
  data.frame(tau_k = tk, k_shape = ks, M = M, rev_share = mean(rs), rev_share_se = se(rs), rev_share_med = median(rs),
             rev_gain = mean(rg), rev_gain_se = se(rg), targ = mean(tg), targ_se = se(tg), targ_med = median(tg),
             ess_med = median(d[, "ess_med"]), ess_q10 = median(d[, "ess_q10"]), meanK = mean(d[, "meanK"]), meanK_med = median(d[, "meanK"]))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
MS <- c(200, 1000, 5000); out <- list(); t0 <- Sys.time()
for (M in MS) for (ks in c(1.3, 2, 3.5)) for (tk in c(0.3, 1, 3)) {
  out[[length(out) + 1]] <- cell(tk, ks, M)
  cat(M, ks, tk, round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
}
write.csv(do.call(rbind, out), "audit_results.csv", row.names = FALSE)
cat("AUDIT DONE\n")
