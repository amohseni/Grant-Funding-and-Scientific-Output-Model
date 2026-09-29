# Package C: review-signal robustness. Fixed mean capability 2, b = 1, T = 2 unless stated,
# 50 populations. C1 resource-signal noise tau_R; C2 multiplicative review noise
# log S = log K + e; C3 fresh review each round (the model's review is a single score
# observed at the start and reused in every round; here it is redrawn every round).
suppressPackageStartupMessages(library(parallel))
source("model.R"); formals(run_simulation_T)$M <- 1000
draw_signals_add <- draw_signals; loglik_add <- loglik_grant_signal; build_hist_orig <- build_posteriors_hist
metrics <- function(d) {
  m <- function(x) mean(x, na.rm = TRUE); se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  rs <- (d[, "S5"] - d[, "S4"]) / (d[, "OR"] - d[, "S4"]); rg <- (d[, "S5"] - d[, "S4"]) / (d[, "S5"] - d[, "S1"])
  tg <- (d[, "OR"] - d[, "S2"]) / (d[, "S2"] - d[, "S1"])
  data.frame(rev_share = m(rs), rev_share_se = se(rs), rev_gain = m(rg), rev_gain_se = se(rg), targ_unif = m(tg),
             short4 = m((d[, "OR"] - d[, "S4"]) / (d[, "OR"] - d[, "S1"])), short5 = m((d[, "OR"] - d[, "S5"]) / (d[, "OR"] - d[, "S1"])),
             sk_cor = m(d[, "sk_cor"]))
}
run_cell <- function(ks, tk, tr = 1, seeds = 1:50, T_rounds = 2, use_rs = TRUE) {
  kmin <- 2 * (ks - 1) / ks
  rows <- mclapply(seeds, function(sd) {
    rp <- run_simulation_T(seed = sd, T_rounds = T_rounds, tau_k = tk, tau_r = tr, k_shape = ks, k_min = kmin,
                           use_resource_signal = use_rs, allocator = "smooth", strategies = c(1, 2, 4, 5), detail = TRUE)
    ro <- run_simulation_T(seed = sd, T_rounds = T_rounds, tau_k = tk, tau_r = tr, k_shape = ks, k_min = kmin,
                           use_resource_signal = use_rs, allocator = "smooth", strategies = c(5), oracle = TRUE, detail = TRUE)
    sk <- rp$signals$sigma_k; if (is.matrix(sk)) sk <- sk[, 1]
    c(S1 = rp$strategies[[1]]$total_expected, S2 = rp$strategies[[2]]$total_expected,
      S4 = rp$strategies[[4]]$total_expected, S5 = rp$strategies[[5]]$total_expected,
      OR = ro$strategies[[5]]$total_expected, sk_cor = cor(sk, rp$K_at_start, method = "spearman"))
  }, mc.cores = 2)
  metrics(do.call(rbind, rows))
}
out <- list(); t0 <- Sys.time()
# ---- C1: tau_R sweep at tau_K = 1, and no resource signal at all
for (ks in c(1.3, 2, 3.5)) {
  for (tr in c(0.1, 0.3, 1, 3, 10)) out[[length(out) + 1]] <- cbind(label = "tauR", k_shape = ks, tau_k = 1, par = tr, run_cell(ks, 1, tr = tr))
  out[[length(out) + 1]] <- cbind(label = "noRsignal", k_shape = ks, tau_k = 1, par = NA, run_cell(ks, 1, use_rs = FALSE))
}
cat("C1 done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
# ---- C2: multiplicative review noise, log S = log K + e, e ~ N(0, s^2); funder knows the form
draw_signals <- function(K0, R0, tau_r, tau_k) {
  n <- length(K0); list(sigma_r = R0 + rnorm(n, 0, tau_r), sigma_k = K0 * exp(rnorm(n, 0, tau_k)))
}
loglik_grant_signal <- function(sigma_k, K, tau_k) dnorm(log(sigma_k), mean = log(K), sd = tau_k, log = TRUE)
for (ks in c(1.3, 2, 3.5)) for (s in c(0.1, 0.3, 0.6, 1, 1.5))
  out[[length(out) + 1]] <- cbind(label = "mult", k_shape = ks, tau_k = s, par = s, run_cell(ks, s))
cat("C2 done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
draw_signals <- draw_signals_add; loglik_grant_signal <- loglik_add
# ---- C3: fresh review each round. sigma_k becomes an n x T matrix; the posterior at round t
# uses columns 1..t (the reviews observed so far).
FRESH_T <- 2
draw_signals <- function(K0, R0, tau_r, tau_k) {
  n <- length(K0)
  list(sigma_r = R0 + rnorm(n, 0, tau_r), sigma_k = matrix(K0 + rnorm(n * FRESH_T, 0, tau_k), n, FRESH_T))
}
build_posteriors_hist <- function(p_init, sigma_r, sigma_k, M, k_min, k_shape, r_min, r_shape, gamma,
                                  tau_r, tau_k, use_resource_signal, use_grant_signal,
                                  p_hist = NULL, g_hist = NULL, gamma_i = NULL, latentA = NULL, k_lognormal = NULL) {
  n <- length(p_init); t_now <- length(p_hist) + 1
  posts <- vector("list", n)
  for (i in seq_len(n)) {
    K_s <- rpareto(M, k_min, k_shape); R_s <- rpareto(M, r_min, r_shape)
    ll <- loglik_pubs(p_init[i], K_s, R_s, gamma)
    if (!is.null(p_hist) && length(p_hist) > 0) for (s in seq_along(p_hist)) ll <- ll + loglik_pubs(p_hist[[s]][i], K_s, R_s + g_hist[[s]][i], gamma)
    if (use_resource_signal) ll <- ll + loglik_resource_signal(sigma_r[i], R_s, tau_r)
    if (use_grant_signal) for (s in seq_len(min(t_now, ncol(sigma_k)))) ll <- ll + loglik_grant_signal(sigma_k[i, s], K_s, tau_k)
    ll <- ll - max(ll); w <- exp(ll); w <- w / sum(w)
    posts[[i]] <- list(K0 = K_s, R0 = R_s, w = w, gam = NULL)
  }
  posts
}
for (ks in c(1.3, 2, 3.5)) for (tk in c(0.3, 1, 3))
  out[[length(out) + 1]] <- cbind(label = "fresh", k_shape = ks, tau_k = tk, par = tk, run_cell(ks, tk))
cat("C3 T=2 done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
write.csv(do.call(rbind, out), "signals_results.csv", row.names = FALSE)
# C3 by round, T = 20: shortfall of records-only and of fresh-review funder by round
FRESH_T <- 20
byround <- function(ks, seeds = 1:50) {
  kmin <- 2 * (ks - 1) / ks
  rows <- mclapply(seeds, function(sd) {
    rp <- run_simulation_T(seed = sd, T_rounds = 20, tau_k = 1, k_shape = ks, k_min = kmin, allocator = "smooth", strategies = c(1, 4, 5), detail = TRUE)
    ro <- run_simulation_T(seed = sd, T_rounds = 20, tau_k = 1, k_shape = ks, k_min = kmin, allocator = "smooth", strategies = c(5), oracle = TRUE, detail = TRUE)
    l1 <- sapply(rp$strategies[[1]]$lam_rounds, sum); l4 <- sapply(rp$strategies[[4]]$lam_rounds, sum)
    l5 <- sapply(rp$strategies[[5]]$lam_rounds, sum); lo <- sapply(ro$strategies[[5]]$lam_rounds, sum)
    cbind(round = seq_along(l1), sf4 = (lo - l4) / (lo - l1), sf5 = (lo - l5) / (lo - l1))
  }, mc.cores = 2)
  a <- as.data.frame(do.call(rbind, rows))
  agg <- aggregate(cbind(sf4, sf5) ~ round, data = a, FUN = mean); agg$k_shape <- ks; agg
}
br <- rbind(byround(2), byround(1.3))
write.csv(br, "byround_fresh_results.csv", row.names = FALSE)
cat("C3 T=20 done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\nSIGNALS DONE\n")
