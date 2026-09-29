# Package E1: regime map. Budget scale b (paper units; code b' = b/2) x capability tail alpha_K
# at fixed mean capability 2. T = 2, 50 populations, myopic funders with and without review at
# tau_K = 1, complete-information benchmark, uniform, no funding. Records the Gini of the
# capability draw and of the optimal grants.
suppressPackageStartupMessages(library(parallel))
source("model.R"); formals(run_simulation_T)$M <- 1000
gap_rule <- function(K, R, budget) {
  if (budget <= 0) return(rep(0, length(K)))
  f <- function(cc) sum(pmax(cc * K - R, 0)) - budget
  cc <- uniroot(f, c(1e-9, 1e7), tol = 1e-10)$root
  pmax(cc * K - R, 0)
}
gini <- function(x) { x <- sort(x); n <- length(x); if (sum(x) == 0) return(0); sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x)) }
cell <- function(b, ks, seeds = 1:50, tk = 1) {
  kmin <- 2 * (ks - 1) / ks
  rows <- mclapply(seeds, function(sd) {
    rp <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks, k_min = kmin, b = b / 2,
                           allocator = "smooth", strategies = c(1, 2, 4, 5), detail = TRUE)
    ro <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks, k_min = kmin, b = b / 2,
                           allocator = "smooth", strategies = c(5), oracle = TRUE)
    K <- rp$K_at_start; R <- rp$R0_at_start
    g <- gap_rule(K, R, rp$params$B_total / 2)
    c(S1 = rp$strategies[[1]]$total_expected, S2 = rp$strategies[[2]]$total_expected,
      S4 = rp$strategies[[4]]$total_expected, S5 = rp$strategies[[5]]$total_expected,
      OR = ro$strategies[[5]]$total_expected, gini_K = gini(K), gini_g = gini(g), n_funded = sum(g > 1e-9),
      meanK = mean(K), meanR = mean(R))
  }, mc.cores = 2)
  d <- do.call(rbind, rows)
  m <- function(x) mean(x, na.rm = TRUE); se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  targ_unif <- (d[, "OR"] - d[, "S2"]) / (d[, "S2"] - d[, "S1"])       # value of targeting, share of uniform's gain
  targ_fund <- (d[, "OR"] - d[, "S2"]) / (d[, "OR"] - d[, "S1"])       # share of the optimum's gain lost by uniform funding
  rev_share <- (d[, "S5"] - d[, "S4"]) / (d[, "OR"] - d[, "S4"])
  rev_gain <- (d[, "S5"] - d[, "S4"]) / (d[, "S5"] - d[, "S1"])       # review's value, share of funding's gain
  short4 <- (d[, "OR"] - d[, "S4"]) / (d[, "OR"] - d[, "S1"])          # records-only shortfall
  data.frame(b = b, k_shape = ks, tau_k = tk, targ_unif = m(targ_unif), targ_unif_se = se(targ_unif),
             targ_fund = m(targ_fund), rev_share = m(rev_share), rev_share_se = se(rev_share),
             rev_gain = m(rev_gain), rev_gain_se = se(rev_gain), short4 = m(short4),
             gini_K = m(d[, "gini_K"]), gini_g = m(d[, "gini_g"]), n_funded = m(d[, "n_funded"]),
             meanK = m(d[, "meanK"]), meanR = m(d[, "meanR"]))
}
out <- list(); t0 <- Sys.time()
for (ks in c(1.2, 1.3, 1.5, 2, 2.5, 3.5, 5)) for (b in c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3)) {
  out[[length(out) + 1]] <- cell(b, ks)
  cat(ks, b, round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
}
write.csv(do.call(rbind, out), "regime_results.csv", row.names = FALSE)
cat("REGIME DONE\n")
