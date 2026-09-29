# Sensitivity runs for the revision (work package C). Harmonic model, smooth allocator,
# T = 2, myopic funders S4 (records + resource signal) and S5 (S4 + review), oracle =
# complete information, S1 no funding, S2 uniform. Per cell we record:
#   rev_share  = (S5 - S4)/(OR - S4)   share of records-only shortfall recovered by review
#   rev_gain   = (S5 - S4)/(S5 - S1)   review's value as share of the gain from funding
#   targ_unif  = (OR - S2)/(S2 - S1)   value of targeting relative to uniform's gain
#   gap_cv     = CV of positive round-1 gaps cK - R0 under complete information
#   gap_gini   = Gini of round-1 optimal grants
#   kr_cor     = Spearman correlation of K0 and R0 in the drawn population
#   sk_cor     = Spearman correlation of the review score and K0 (informativeness)
suppressPackageStartupMessages(library(parallel))
source("model.R")
gap_rule <- function(K, R, budget) {
  f <- function(cc) sum(pmax(cc * K - R, 0)) - budget
  cc <- uniroot(f, c(1e-9, 1e6), tol = 1e-10)$root
  pmax(cc * K - R, 0)
}
gini <- function(x) { x <- sort(x); n <- length(x); if (sum(x) == 0) return(0); sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x)) }
cell <- function(tk, ks, rho = 0, rs = 2, kmin = 1, n = 50, seeds = 1:50, label = "") {
  rows <- mclapply(seeds, function(sd) {
    rp <- run_simulation_T(seed = sd, T_rounds = 2, n = n, tau_k = tk, k_shape = ks, r_shape = rs, rho_kr = rho, k_min = kmin,
                           allocator = "smooth", strategies = c(1, 2, 4, 5), detail = TRUE)
    ro <- run_simulation_T(seed = sd, T_rounds = 2, n = n, tau_k = tk, k_shape = ks, r_shape = rs, rho_kr = rho, k_min = kmin,
                           allocator = "smooth", strategies = c(5), oracle = TRUE)
    K <- rp$K_at_start; R <- rp$R0_at_start
    g <- gap_rule(K, R, rp$params$B_total / 2)
    gaps <- g[g > 0]
    set.seed(sd + 10^6); S <- K + rnorm(length(K), 0, tk)
    c(S1 = rp$strategies[[1]]$total_expected, S2 = rp$strategies[[2]]$total_expected,
      S4 = rp$strategies[[4]]$total_expected, S5 = rp$strategies[[5]]$total_expected,
      OR = ro$strategies[[5]]$total_expected,
      gap_cv = if (length(gaps) > 1) sd(gaps) / mean(gaps) else 0, gap_gini = gini(g),
      kr_cor = cor(K, R, method = "spearman"), sk_cor = cor(S, K, method = "spearman"),
      meanK = mean(K), meanR = mean(R))
  }, mc.cores = 2)
  d <- do.call(rbind, rows)
  m <- function(x) mean(x); se <- function(x) sd(x) / sqrt(length(x))
  rev_share <- (d[, "S5"] - d[, "S4"]) / (d[, "OR"] - d[, "S4"])
  rev_gain <- (d[, "S5"] - d[, "S4"]) / (d[, "S5"] - d[, "S1"])
  targ <- (d[, "OR"] - d[, "S2"]) / (d[, "S2"] - d[, "S1"])
  data.frame(label = label, tau_k = tk, k_shape = ks, rho = rho, r_shape = rs, k_min = kmin, n = n,
             rev_share = m(rev_share), rev_share_se = se(rev_share), rev_gain = m(rev_gain), rev_gain_se = se(rev_gain),
             targ_unif = m(targ), targ_unif_se = se(targ), gap_cv = m(d[, "gap_cv"]), gap_gini = m(d[, "gap_gini"]),
             kr_cor = m(d[, "kr_cor"]), sk_cor = m(d[, "sk_cor"]), meanK = m(d[, "meanK"]), meanR = m(d[, "meanR"]))
}
out <- list(); t0 <- Sys.time()
# C1a: rho grid
for (rho in c(-0.5, 0, 0.5, 0.8)) for (ks in c(1.3, 2, 3.5)) for (tk in c(0.3, 1, 3))
  out[[length(out) + 1]] <- cell(tk, ks, rho = rho, label = "rho")
cat("C1a done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
# C1b: alpha_R grid
for (rs in c(1.3, 3.5)) for (ks in c(1.3, 2, 3.5)) for (tk in c(0.3, 1, 3))
  out[[length(out) + 1]] <- cell(tk, ks, rs = rs, label = "rshape")
cat("C1b done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
# C2a: mean-normalized capability, E[K] = 2 for every alpha, fine tau grid, 100 seeds
for (ks in c(1.3, 2, 3.5)) for (tk in c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10, 20))
  out[[length(out) + 1]] <- cell(tk, ks, kmin = 2 * (ks - 1) / ks, seeds = 1:100, label = "meanK2")
cat("C2a done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
# C2b: the original (scale 1) grid with 100 seeds for error bars
for (ks in c(1.3, 2, 3.5)) for (tk in c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10, 20))
  out[[length(out) + 1]] <- cell(tk, ks, seeds = 1:100, label = "scale1")
cat("C2b done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
# C2c: applicant-pool size
for (n in c(25, 100, 200)) for (ks in c(1.3, 2, 3.5)) for (tk in c(0.3, 1, 3))
  out[[length(out) + 1]] <- cell(tk, ks, n = n, seeds = 1:50, label = "n")
cat("C2c done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
res <- do.call(rbind, out)
write.csv(res, "sens_results.csv", row.names = FALSE)
cat("ALL DONE\n")
