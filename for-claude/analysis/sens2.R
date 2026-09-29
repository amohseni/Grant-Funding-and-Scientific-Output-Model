# Second batch: (1) mean-normalized fine grid with gap-rule correlations and shortfalls
# (replaces the Fig. 2C data); (2) review contaminated by resources, S = K + beta R + eta,
# funder knows beta (C4); (3) T = 20 records-only vs review, output shortfall by round (C5).
suppressPackageStartupMessages(library(parallel))
source("model.R")
gap_rule <- function(K, R, budget) {
  f <- function(cc) sum(pmax(cc * K - R, 0)) - budget
  cc <- uniroot(f, c(1e-9, 1e6), tol = 1e-10)$root
  pmax(cc * K - R, 0)
}
BETA <- 0
# contaminated review: patch the signal draw and both likelihood sites
draw_signals <- function(K0, R0, tau_r, tau_k) {
  n <- length(K0)
  list(sigma_r = R0 + rnorm(n, 0, tau_r), sigma_k = K0 + BETA * R0 + rnorm(n, 0, tau_k))
}
loglik_grant_signal2 <- function(sigma_k, K, R, tau_k) dnorm(sigma_k, mean = K + BETA * R, sd = tau_k, log = TRUE)
for (fn in c("posterior_samples_single", "build_posteriors_hist")) {
  src <- deparse(get(fn))
  src <- gsub("loglik_grant_signal(sigma_k, K_s, tau_k)", "loglik_grant_signal2(sigma_k, K_s, R_s, tau_k)", src, fixed = TRUE)
  src <- gsub("loglik_grant_signal(sigma_k[i], K_s, tau_k)", "loglik_grant_signal2(sigma_k[i], K_s, R_s, tau_k)", src, fixed = TRUE)
  assign(fn, eval(parse(text = src)))
}
cell <- function(tk, ks, kmin = 1, seeds = 1:100, beta = 0, label = "") {
  BETA <<- beta
  rows <- mclapply(seeds, function(sd) {
    rp <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks, k_min = kmin,
                           allocator = "smooth", strategies = c(1, 2, 4, 5), detail = TRUE)
    ro <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks, k_min = kmin,
                           allocator = "smooth", strategies = c(5), oracle = TRUE)
    K <- rp$K_at_start; R <- rp$R0_at_start
    gstar <- gap_rule(K, R, rp$params$B_total / 2)
    g5 <- rp$strategies[[5]]$g_rounds[[1]]; g4 <- rp$strategies[[4]]$g_rounds[[1]]
    c(S1 = rp$strategies[[1]]$total_expected, S2 = rp$strategies[[2]]$total_expected,
      S4 = rp$strategies[[4]]$total_expected, S5 = rp$strategies[[5]]$total_expected,
      OR = ro$strategies[[5]]$total_expected,
      corr_S5 = suppressWarnings(cor(g5, gstar)), corr_S4 = suppressWarnings(cor(g4, gstar)))
  }, mc.cores = 2)
  d <- do.call(rbind, rows)
  m <- function(x) mean(x, na.rm = TRUE); se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  rs <- (d[, "S5"] - d[, "S4"]) / (d[, "OR"] - d[, "S4"])
  data.frame(label = label, tau_k = tk, k_shape = ks, k_min = kmin, beta = beta,
             rev_share = m(rs), rev_share_se = se(rs),
             rev_gain = m((d[, "S5"] - d[, "S4"]) / (d[, "S5"] - d[, "S1"])),
             short4 = m((d[, "OR"] - d[, "S4"]) / (d[, "OR"] - d[, "S1"])), short5 = m((d[, "OR"] - d[, "S5"]) / (d[, "OR"] - d[, "S1"])),
             corr_S5 = m(d[, "corr_S5"]), corr_S4 = m(d[, "corr_S4"]), targ_unif = m((d[, "OR"] - d[, "S2"]) / (d[, "S2"] - d[, "S1"])))
}
out <- list(); t0 <- Sys.time()
for (ks in c(1.3, 2, 3.5)) for (tk in c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10, 20))
  out[[length(out) + 1]] <- cell(tk, ks, kmin = 2 * (ks - 1) / ks, label = "norm_fine")
cat("norm_fine done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
for (beta in c(0.25, 0.5, 1)) for (ks in c(1.3, 2, 3.5)) for (tk in c(0.3, 1, 3))
  out[[length(out) + 1]] <- cell(tk, ks, kmin = 2 * (ks - 1) / ks, beta = beta, seeds = 1:50, label = "contam")
cat("contam done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
BETA <<- 0
write.csv(do.call(rbind, out), "sens2_results.csv", row.names = FALSE)
# C5: T = 20 by-round shortfall, default field and heavy, mean-normalized, 50 seeds
byround <- function(ks, kmin, seeds = 1:50) {
  rows <- mclapply(seeds, function(sd) {
    rp <- run_simulation_T(seed = sd, T_rounds = 20, tau_k = 1, k_shape = ks, k_min = kmin, allocator = "smooth", strategies = c(1, 4, 5), detail = TRUE)
    ro <- run_simulation_T(seed = sd, T_rounds = 20, tau_k = 1, k_shape = ks, k_min = kmin, allocator = "smooth", strategies = c(5), oracle = TRUE)
    l1 <- sapply(rp$strategies[[1]]$lam_rounds, sum); l4 <- sapply(rp$strategies[[4]]$lam_rounds, sum)
    l5 <- sapply(rp$strategies[[5]]$lam_rounds, sum); lo <- sapply(ro$strategies[[5]]$lam_rounds, sum)
    cbind(round = seq_along(l1), sf4 = (lo - l4) / (lo - l1), sf5 = (lo - l5) / (lo - l1))
  }, mc.cores = 2)
  a <- do.call(rbind, rows)
  agg <- aggregate(cbind(sf4, sf5) ~ round, data = as.data.frame(a), FUN = mean)
  agg$k_shape <- ks; agg
}
br <- rbind(byround(2, 1), byround(1.3, 2 * 0.3 / 1.3))
write.csv(br, "byround_results.csv", row.names = FALSE)
cat("byround done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\nALL DONE\n")
